from pathlib import Path
import sys,json,hashlib,subprocess,zipfile,math
import pandas as pd
from biom import Table
from qiime2 import Artifact,Visualization
sys.path.insert(0,'/workspaces/Meta2Data/tests/integration')
import validate_runtime as runtime
root=Path('/workspaces/meta2data-validation')
work=root/'taxa-audit-01'
final=work/'results/taxa/final-dada2-notree'
output=work/'first-run-evidence.json'
assert not output.exists()
model=root/'results/db/silva-144-full-length-sklearn-1.7.1.qza'
taxonomy_path=final/'silva144Taxonomy.qza'
identity=runtime.check_taxonomy_classifier(taxonomy_path,model)
taxonomy=Artifact.load(taxonomy_path).view(pd.DataFrame)
ranks=runtime.check_silva144_ranks(taxonomy['Taxon'])
table=Artifact.load(final/'orientedTable.qza').view(Table)
seqs=Artifact.load(final/'orientedRepSeqs.qza').view(pd.Series)
runtime.check_identifiers(table.ids(axis='observation'),seqs.index,taxonomy.index)
assert set(taxonomy.index)==set(seqs.index) and not taxonomy.index.has_duplicates
assert len(seqs)>0 and all(math.isfinite(float(x)) and 0<=float(x)<=1 for x in taxonomy['Confidence'])
state=json.loads((final/'taxa-run-state.json').read_text())
assert state['classifier_id']=='silva-144' and state['db_label']=='silva144' and state['confidence']==0.7
assert not (final/'silvaTaxonomy.qza').exists() and not (final/'gg2Taxonomy.qza').exists()
files={}
for path in sorted(final.rglob('*.qz*')):
    obj=Visualization.load(path) if path.suffix=='.qzv' else Artifact.load(path)
    obj.validate(level='min')
    files[str(path.relative_to(final))]={'uuid':str(obj.uuid),'type':str(obj.type),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'mtime_ns':path.stat().st_mtime_ns}
for summary in ('mergedTableSummary','orientedTableSummary'):
    assert files['tmp/'+summary+'.qzv']['type']=='Visualization'
    for suffix in ('feature-frequencies','sample-frequencies'):
        assert files['tmp/'+summary+'-'+suffix+'.qza']['type']=='ImmutableMetadata'
counts=dict(zip(map(str,table.ids(axis='sample')),map(int,table.sum(axis='sample'))))
assert sum(counts.values())>0
report={'status':'passed','scope':'SILVA144 real no-tree run; not the full release gate','run_revision':'82d41e0b174ac2934e3c614add80629bb159d2a9','validator_revision':subprocess.check_output(['git','-C','/workspaces/Meta2Data','rev-parse','HEAD'],text=True).strip(),'requested_classifier':'silva-144','default_confidence':0.7,'features':len(taxonomy),'sample_counts':counts,**identity,**ranks,'files':files,'state':state}
output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
taxonomy.to_csv(work/'first-run-taxonomy.tsv',sep='\t')
import shutil
shutil.copy2(taxonomy_path,work/'first-run-taxonomy-confidence-0.7.qza')
print(json.dumps({k:v for k,v in report.items() if k not in ('files','state')},sort_keys=True))
print('TAXA_FIRST_CHECK_PASSED',len(files))
