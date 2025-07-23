import os

os.system('make')
for i in range(32,81):
	os.system('./propagate.exe '+str(i))	
