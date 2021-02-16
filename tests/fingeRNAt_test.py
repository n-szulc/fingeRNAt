import subprocess, os, sys, traceback


RNA = ['1aju_model1.pdb', '3d2v.mol2']
ligands = [['ligands.sdf'], ['redocked.sdf', 'various.sdf']]
fingerprints = ['SIMPLE', 'PBS', 'FULL']

def run_test():

    os.system("rm -rf outputs/ > /dev/null 2>&1")

    OK = True

    for i in range(len(RNA)):
        for j in range(len(fingerprints)):
            for l in range(len(ligands[i])):
                if subprocess.call('python ../code/fingeRNAt.py -r test_inputs/%s -l test_inputs/%s -f %s -wrapper ACUG,PuPy,Counter -h2o' %(RNA[i], ligands[i][l], fingerprints[j]), shell = True):
                    print ('fingeRNAt had problem running fingerprint %s on test_inputs/%s  and test_inputs/%s' % (fingerprints[j], RNA[i], ligands[i][l]))
                    OK = False

    if OK:
        for root, dirs, files in os.walk('outputs/'):
            files = [f for f in files if not f[0] == '.']
            for name in files:
            	try:
                    out = subprocess.check_output('comm -3 outputs/%s expected_outputs/%s' %(name, name), shell = True)
                    if len(out) != 0:
                        OK = False
                        print ('outputs/%s and expected_outputs/%s differ!' %(name, name))
            	except:
                		mssg = '# Something is wrong, attention needed! #'
		                print('#'*len(mssg))
		                print(mssg)
		                print('#'*len(mssg))
		                traceback.print_exc()
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
