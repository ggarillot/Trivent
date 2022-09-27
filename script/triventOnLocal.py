#!/usr/bin/env python

import os
import sys
import argparse

import Trivent

if __name__ == '__main__' :

	parser = argparse.ArgumentParser()
	parser.add_argument('runNumber' , help='The run number to process' , type=int)
	parser.add_argument('-e' , '--energy' , help='Energy of the run' , type=float , default=0)
	parser.add_argument('-g' , '--geom' , help='Geometry file' , type=str , default='')
	parser.add_argument('-c' , '--cerenkovDelay' , help='Cerenkov clock delay' , type=int , default=0)
	args = parser.parse_args()

	runNumber = str(args.runNumber) 

	dir = '/home/acqilc/Desktop/Analysis/files/STREAMOUT'

	print ('Searching files in ' + dir)

	#list files
	fileList = []

	for fileName in os.listdir(dir) :
		if f'R{runNumber}.slcio' in fileName :
			fileList.append(dir + '/' + fileName)

	fileList.sort()
	print('File List :')
	print(fileList)

	#os.environ["MARLIN"] = '/home/guillaume/ilcsoft/v02-00-01/Marlin/v01-16'
	#os.environ["PATH"] = os.environ["MARLIN"] + '/bin:' + os.environ["PATH"]
	os.environ["MARLIN_DLL"] = '/home/acqilc/Desktop/Analysis/Trivent/lib/libTrivent.so'

	a = Trivent.Params()
	a.energy = args.energy
	a.geometry = args.geom
	a.cerenkovDelay = 0

	a.outputFileName = 'TDHCAL_' + str(sys.argv[1]) + '.slcio'
	a.outputRootFileName = 'noise' + str(sys.argv[1]) + '.root'

	Trivent.launch(a , fileList)

	outputDir = dir.replace('STREAMOUT' , 'TRIVENT')

	os.system('mkdir -p ' + outputDir)
	os.system('mv ' + a.outputFileName + ' ' + outputDir)
	os.system(f'mv noise{runNumber}.root /home/acqilc/Desktop/Analysis/files/noise')
