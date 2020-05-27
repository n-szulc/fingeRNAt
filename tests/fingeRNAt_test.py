import subprocess, os


RNA = ['1aju_model1.pdb']
ligands = ['ligands.sdf']
fingerprints = ['SIMPLE', 'PBS', 'FULL', 'XP']

def run_test():

    os.system("rm -rf outputs/ > /dev/null 2>&1")

    OK = True

    for i in range(len(RNA)):
        for j in range(len(fingerprints)):
            if subprocess.call('python ../code/fingeRNAt.py -r test_inputs/%s -l test_inputs/%s -f %s' %(RNA[i], ligands[i], fingerprints[j]), shell = True):
                print ('fingeRNAt had problem running fingerprint %s on test_inputs/%s  and test_inputs/%s' % (fingerprints[j], RNA[i], ligands[i]))
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
            	except subprocess.CalledProcessError:
                		mssg = '# Something is wrong, attention needed! #'
		                print('#'*len(mssg))
		                print(mssg)
		                print('#'*len(mssg))
		                exit(1)

    return OK


OK = run_test()

if OK:
    mssg = '# Tests ran successfully, everything is OK! #'
    print('#'*len(mssg))
    print(mssg)
    print('#'*len(mssg))
else:
    mssg = '# Something is wrong, attention needed! #'
    print('#'*len(mssg))
    print(mssg)
    print('#'*len(mssg))
    exit(1)
