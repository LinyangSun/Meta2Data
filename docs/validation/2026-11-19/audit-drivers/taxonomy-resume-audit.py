from pathlib import Path
from datetime import datetime,timezone
import sys,os,json,subprocess,time,shutil,hashlib,zipfile,yaml,traceback
root=Path('/workspaces/meta2data-validation')
source=root/'audit-source-550d440'
sys.path.insert(0,str(source/'scripts'))
from taxonomy_provenance import check_taxonomy_provenance
work=root/'taxa-audit-01'
final=work/'results/taxa/final-dada2-notree'
report_path=work/'resume-audit.json'
assert not report_path.exists()
model=root/'results/db/silva-144-full-length-sklearn-1.7.1.qza'
taxonomy=final/'silva144Taxonomy.qza'
reads=final/'orientedRepSeqs.qza'
backup=work/'first-run-taxonomy-confidence-0.7.qza'
assert backup.exists()
env=dict(os.environ,PATH=str(source/'bin')+os.pathsep+os.environ['PATH'],TMPDIR=str(root/'tmp'),PYTHONNOUSERSITE='1',MPLBACKEND='agg',UNIFRAC_USE_GPU='N')
base=[str(source/'bin/Meta2Data'),'AmpliconTAXA','-i',str(root/'workflow-audit-01/pip-dada2/results/pip/pe_explicit'),'--dada2','--classifier','silva-144','--notree','-t','1']
report={'status':'running','scope':'Real taxonomy reuse and stale-content rejection; not the full release gate','source_revision':'550d440e0b127e84b0b412ebd0884fc36b937064','started':datetime.now(timezone.utc).isoformat(),'steps':[]}
def save(): report_path.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
def snapshot():
    files={}
    for path in sorted(final.rglob('*.qz*')):
        with zipfile.ZipFile(path) as z:
            name=[n for n in z.namelist() if n.count('/')==1 and n.endswith('/metadata.yaml')][0]
            metadata=yaml.load(z.read(name),Loader=yaml.BaseLoader)
        files[str(path.relative_to(final))]={'uuid':metadata['uuid'],'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    return files
def run(name,confidence):
    command=base+['--confidence',str(confidence)]
    log=root/'logs'/('resume-audit-'+name+'.log')
    started=time.monotonic()
    with log.open('w') as stream:
        result=subprocess.run(command,cwd=work,env=env,stdout=stream,stderr=subprocess.STDOUT)
    step={'name':name,'command':command,'returncode':result.returncode,'seconds':round(time.monotonic()-started,2),'log':str(log)}
    report['steps'].append(step);save()
    assert result.returncode==0,(name,log.read_text(errors='replace')[-3500:])
    step['identity']=check_taxonomy_provenance(taxonomy,model,reads_path=reads,confidence=confidence)
    step['files']=snapshot();save()
    return step,log.read_text(errors='replace')
def assert_only_taxonomy_changed(before,after):
    assert set(before)==set(after)
    changed=[n for n in before if before[n]!=after[n]]
    assert changed==['silva144Taxonomy.qza'],changed
try:
    save(); initial=snapshot();report['initial_files']=initial;save()
    reuse,log=run('identical-0.7',0.7)
    assert reuse['files']==initial,'Identical run changed artifact contents'
    assert '[OK] reuse: silva144Taxonomy.qza' in log
    reuse['checks']=['all artifact UUIDs and SHA256 values unchanged','taxonomy reuse observed'];save()
    changed,log=run('changed-0.8',0.8)
    assert_only_taxonomy_changed(initial,changed['files'])
    changed['checks']=['only taxonomy UUID and SHA256 changed','confidence provenance is 0.8'];save()
    state_path=final/'taxa-run-state.json'
    state_hash=hashlib.sha256(state_path.read_bytes()).hexdigest()
    shutil.copy2(backup,taxonomy)
    assert hashlib.sha256(state_path.read_bytes()).hexdigest()==state_hash
    assert snapshot()['silva144Taxonomy.qza']==initial['silva144Taxonomy.qza']
    repaired,log=run('restore-stale-0.7-request-0.8',0.8)
    assert 'Cached taxonomy rejected:' in log and 'confidence mismatch' in log
    assert_only_taxonomy_changed(changed['files'],repaired['files'])
    repaired['checks']=['stale 0.7 artifact copied over 0.8 output with current state retained','provenance mismatch rejected','only taxonomy regenerated','regenerated confidence provenance is 0.8'];save()
    report['status']='passed'
except Exception as error:
    report['status']='failed';report['error']=repr(error);traceback.print_exc()
finally:
    report['finished']=datetime.now(timezone.utc).isoformat();save()
print(json.dumps({'status':report['status'],'steps':[{'name':s['name'],'returncode':s['returncode'],'checks':s.get('checks',[])} for s in report['steps']]}),flush=True)
sys.exit(0 if report['status']=='passed' else 1)
