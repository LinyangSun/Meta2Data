from pathlib import Path
from datetime import datetime, timezone
import json, subprocess, sys
repo=Path('/workspaces/Meta2Data')
root=Path('/workspaces/meta2data-validation')
sys.path[:0]=[str(repo/'tests/integration'),str(repo/'scripts')]
from validate_runtime import Validation
from generate_fixture import fasta_records
import reference_resources as refs
v=Validation(root/'workflow-audit-01','pip-platform-audit',2)
v.cache=root/'results/db'
v.report.update(source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),scope='PIP platform and output-file compatibility only; full release gate remains pending',case_results={},errors=[])
v.save()
try:
    reference=refs.prepare_resource('gg2-sequences',work_dir=root)
    fixture,expected,metadata=v.fixtures(reference)
    v.report['fixture_metadata']={name:{'csv':str(value[0]),'datasets':value[1]} for name,value in metadata.items()}
    v.save()
    outputs={}
    for method,(csv_path,names) in metadata.items():
        work=v.case_directory('pip-'+method)
        command=['Meta2Data','AmpliconPIP','--local-m',str(csv_path),'--local-datasets-colNAME','datasets','--local-path-colNAME','path','--local-platform-colNAME','platform','--local-primer-f-colNAME','primer_f','--local-primer-r-colNAME','primer_r','--'+method,'-t','2','--max-parallel','1']
        try:
            v.run('pip-'+method,command,cwd=work)
        except Exception as error:
            v.report['errors'].append({'step':'pip-'+method,'error':str(error)})
            print('PIPELINE_ERROR',method,str(error),flush=True)
        output=work/'results/pip'
        outputs[method]=str(output)
        for name in names:
            case=method+'/'+name
            try:
                samples={sample:record['expected_pipeline_reads'] for sample,record in expected[name]['samples'].items()}
                dataset=output/name
                table=dataset/(name+'-'+method+'-final-table.qza')
                sequences_path=dataset/(name+'-'+method+'-final-rep-seqs.qza')
                _,sequences=v.table_and_sequences(table,sequences_path,samples)
                if name=='pacbio' and method=='dada2':
                    inserts={sequence for _,sequence in fasta_records(fixture/'pacbio_inserts.fasta')}
                    assert set(map(str,sequences))<=inserts,'PacBio representatives differ from primer-free inserts'
                v.report['case_results'][case]={'status':'passed','table':str(table),'sequences':str(sequences_path),'features':len(sequences)}
                print('CASE_PASSED',case,len(sequences),flush=True)
            except Exception as error:
                v.report['case_results'][case]={'status':'failed','error':str(error)}
                v.report['errors'].append({'step':case,'error':str(error)})
                print('CASE_FAILED',case,str(error),flush=True)
            v.save()
    v.report['pip_outputs']=outputs
    v.report['status']='failed' if v.report['errors'] else 'passed'
except Exception as error:
    v.report['status']='failed'
    v.report['errors'].append({'step':'setup','error':str(error)})
    print('AUDIT_ERROR',str(error),flush=True)
v.report['finished']=datetime.now(timezone.utc).isoformat()
v.save()
print('PIP_AUDIT_FINISHED',v.report['status'],flush=True)
sys.exit(0 if v.report['status']=='passed' else 1)
