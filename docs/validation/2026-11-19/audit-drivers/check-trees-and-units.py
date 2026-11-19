from pathlib import Path
from datetime import datetime,timezone
import sys,os,json,time,subprocess,hashlib,re,traceback,math
root=Path('/workspaces/meta2data-validation')
source=root/'audit-source-550d440'
sys.path.insert(0,str(source/'tests/integration'))
sys.path.insert(0,str(source/'scripts'))
import validate_runtime as runtime
from taxonomy_provenance import check_taxonomy_provenance
print('Waiting for real workflows to finish before isolated verification',flush=True)
deadline=time.monotonic()+1800
while not (root/'tree-audit-01.exit').exists():
    if time.monotonic()>deadline: raise RuntimeError('Tree workflow has not finished within the verification window')
    time.sleep(2)
assert (root/'tree-audit-01.exit').read_text().strip()=='0'
assert (root/'workflow-audit-01.exit').read_text().strip()=='0'
assert json.loads((root/'taxa-audit-01/resume-audit.json').read_text())['status']=='passed'
import pandas as pd
from biom import Table
from qiime2 import Artifact,Visualization
from skbio import TreeNode
model=root/'results/db/silva-144-full-length-sklearn-1.7.1.qza'
report={'status':'running','scope':'SILVA144 tree output compatibility, not the all-classifier release gate','run_revision':'86dea0b3d5ad3b060588bb5466060138df96d6f3','validator_revision':'550d440e0b127e84b0b412ebd0884fc36b937064','cases':{}}
report_path=root/'tree-audit-01/tree-evidence.json'
assert not report_path.exists()
def save(): report_path.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
try:
    for mode,folder in [('singleV','final-dada2-singleV'),('sepp','final-dada2-multiV')]:
        final=root/'tree-audit-01'/mode/'results/taxa'/folder
        taxonomy_path=final/'silva144Taxonomy.qza'
        oriented=final/('tmp/orientedRepSeqs.qza' if mode=='sepp' else 'orientedRepSeqs.qza')
        identity=check_taxonomy_provenance(taxonomy_path,model,reads_path=oriented,confidence=0.7)
        taxonomy=Artifact.load(taxonomy_path).view(pd.DataFrame)
        ranks=runtime.check_silva144_ranks(taxonomy['Taxon'])
        table_path=final/('treeFilteredTable.qza' if mode=='sepp' else 'orientedTable.qza')
        sequence_path=final/('treeFilteredRepSeqs.qza' if mode=='sepp' else 'orientedRepSeqs.qza')
        tree_path=final/('seppTree.qza' if mode=='sepp' else 'denovoRootedTree.qza')
        table=Artifact.load(table_path).view(Table)
        seqs=Artifact.load(sequence_path).view(pd.Series)
        tree=Artifact.load(tree_path)
        assert str(tree.type)=='Phylogeny[Rooted]'
        tips={tip.name for tip in tree.view(TreeNode).tips()}
        runtime.check_identifiers(table.ids(axis='observation'),seqs.index,taxonomy.index,tips)
        assert len(seqs)>0 and not seqs.index.has_duplicates and not taxonomy.index.has_duplicates
        if mode=='singleV': assert tips==set(seqs.index)==set(taxonomy.index)
        assert all(math.isfinite(float(x)) and 0<=float(x)<=1 for x in taxonomy['Confidence'])
        counts=dict(zip(map(str,table.ids(axis='sample')),map(int,table.sum(axis='sample'))))
        assert set(counts)=={'explicitPE_A','explicitPE_B'} and all(0<n<=8000 for n in counts.values())
        files={}
        for path in sorted(final.rglob('*.qz*')):
            obj=Visualization.load(path) if path.suffix=='.qzv' else Artifact.load(path)
            obj.validate(level='min')
            files[str(path.relative_to(final))]={'uuid':str(obj.uuid),'type':str(obj.type),'bytes':path.stat().st_size}
        for stem in ('mergedTableSummary','treeFilteredTableSummary' if mode=='sepp' else 'orientedTableSummary'):
            assert files['tmp/'+stem+'.qzv']['type']=='Visualization'
            for suffix in ('feature-frequencies','sample-frequencies'):
                assert files['tmp/'+stem+'-'+suffix+'.qza']['type']=='ImmutableMetadata'
        report['cases'][mode]={'status':'passed','retained_features':len(seqs),'taxonomy_features':len(taxonomy),'tree_tips':len(tips),'sample_counts':counts,**identity,**ranks,'files':files}
        save();print('TREE_OUTPUTS_PASSED',mode,len(seqs),flush=True)
    report['status']='passed';save()
except Exception as error:
    report['status']='failed';report['error']=repr(error);save();raise
processes=subprocess.check_output(['ps','-eo','pid,ppid,args'],text=True)
leftovers=[line for line in processes.splitlines() if re.search(r'/workspaces/meta2data-validation/tmp/tmp[^/]+/repo/(bin|scripts)/',line)]
assert not leftovers,('Leftover timed-out test processes require inspection',leftovers)
v=runtime.Validation(root/'unit-audit-03','isolated-unit-and-shell-audit',1)
v.report['source_revision']='550d440e0b127e84b0b412ebd0884fc36b937064'
v.report['prior_parallel_attempt']={'directory':str(root/'unit-audit-02'),'tests':317,'failures':0,'errors':3,'skipped':0,'cause_observed':'Three subprocess.TimeoutExpired errors at the existing 30-second bound while real workflows were running'}
try:
    v.unit_tests()
    v.run('help-AmpliconTAXA',[source/'bin/Meta2Data','AmpliconTAXA','--help'])
    v.report['status']='passed'
except Exception as error:
    v.report['status']='failed';v.report['error']=repr(error);traceback.print_exc()
finally:
    v.report['finished']=datetime.now(timezone.utc).isoformat();v.save()
print('ISOLATED_RESULT',json.dumps({k:v.report.get(k) for k in ('status','unit_tests','error')}),flush=True)
sys.exit(0 if v.report['status']=='passed' else 1)
