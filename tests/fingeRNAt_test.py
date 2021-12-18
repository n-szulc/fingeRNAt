import subprocess, os, sys, traceback, shutil

sys_sep = os.sep
RNA = ['3d2v.pdb']
ligands = [['redocked.sdf']]
fingerprints = ['FULL']
yaml = 'custom-interactions.yaml'

program_path = '..' + sys_sep + 'code' + sys_sep + 'fingeRNAt.py'
test_inputs_path = 'test_inputs' + sys_sep
test_outputs_path = 'outputs' + sys_sep
test_ex_outputs_path = 'expected_outputs' + sys_sep

def run_test():

    shutil.rmtree('outputs', ignore_errors=True)
    os.mkdir('outputs')

    OK = True

    for i in range(len(RNA)):
        for j in range(len(fingerprints)):
            for l in range(len(ligands[i])):
                if fingerprints[j] == 'FULL':
                    #subprocess.call("> %sdebug.txt" %(test_outputs_path), shell = True)
                    command = 'python %s -r %s -l %s -f %s -h2o -o outputs -custom %s' %(program_path, test_inputs_path + RNA[i], test_inputs_path + ligands[i][l], fingerprints[j], test_inputs_path + yaml)

                    #command = 'python %s -r %s -l %s -f %s -wrapper ACUG,PuPy,Counter -h2o -o outputs -detail' %(program_path, test_inputs_path + RNA[i], test_inputs_path + ligands[i][l], fingerprints[j])
                if subprocess.call(command, shell = True):
                    print ('fingeRNAt had problem running %s on %s  and %s' % (fingerprints[j], test_inputs_path + RNA[i], test_inputs_path + ligands[i][l]))
                    OK = False



    if OK:
        for root, dirs, files in os.walk('outputs'):
            files = [f for f in files if not f[0] == '.']
            for name in files:
            	try:
                    #subprocess.call('cat %s' %(test_outputs_path + name), shell = True)

                    out = subprocess.check_output('comm -3 %s %s' %(test_outputs_path + name, test_ex_outputs_path + name), shell = True)
                    if len(out) != 0:
                        OK = False
                        print ('%s and %s differ!' %(test_outputs_path + name, test_ex_outputs_path + name))
                        subprocess.call('cat %s' %(test_outputs_path + name), shell = True)

            	except:
                		mssg = '# Something is wrong, attention needed! #'
		                print('#'*len(mssg))
		                print(mssg)
		                print('#'*len(mssg))
		                traceback.print_exc()
		                subprocess.call('cat %s' %(test_outputs_path + name), shell = True)
		                sys.exit(3)

    return OK

OK = run_test()

if OK:
    mssg = '# Tests ran successfully, everything is OK! #'
    print('#'*len(mssg))
    print(mssg)
    print('#'*len(mssg))
    sys.exit(0)
else:
    mssg = '# Something is wrong, attention needed! #'
    print('#'*len(mssg))
    print(mssg)
    print('#'*len(mssg))
    sys.exit(1)
