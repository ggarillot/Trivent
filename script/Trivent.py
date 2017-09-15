#!/usr/bin/env python

import os
import sys
import argparse

class Params :
	def __init__(self) :
		self.energy = 0
		self.layerCut = 5
		self.noiseCut = 7
		self.timeWin = 2
		self.cerenkovDelay = 0
		self.removeSquareEvents = True
		self.removeRamFullEvents = True
		self.geometry = ""


def launch(a , files) :

	fileList = ''
	for name in files :
		fileList += name + ' '

	pid = os.getpid()

	xmlFileName = str(pid) + '.xml'
	tempOutputFile = str(pid) + '.slcio'

	xml = '''<marlin>

  <execute>
    <processor name="MyTriventProc"/>
  </execute>

  <global>
    <parameter name="LCIOInputFiles">''' + fileList + '''</parameter>
    <parameter name="Verbosity" type="string"> MESSAGE </parameter>
    <parameter name="SupressCheck" value="false"/>
  </global>

  <processor name="MyTriventProc" type="TriventProc">
    <parameter name="HitCollectionName" type="StringVec"> DHCALRawHits </parameter>
    <parameter name="setup_geometry" type="string">''' + a.geometry + '''</parameter>
    <parameter name="beamEnergy" type="double">''' + str(a.energy) + '''</parameter>
	<parameter name="electronic_noise_cut" type="int"> 500000 </parameter>
    <parameter name="cerenkovDelay" type="int">''' + str(a.cerenkovDelay) + '''</parameter>
    <parameter name="RemoveSquareEvents" type="bool">''' + str(a.removeSquareEvents).lower() + '''</parameter>
    <parameter name="RemoveRamFullEvents" type="bool">''' + str(a.removeRamFullEvents).lower() + '''</parameter>
    <parameter name="LayerCut" type="int">''' + str(a.layerCut) + '''</parameter>
    <parameter name="noiseCut" type="int">''' + str(a.noiseCut) + '''</parameter>
    <parameter name="timeWin" type="int">''' + str(a.timeWin) + '''</parameter>
	<parameter name="LCIOOutputFile" type="string" >''' + tempOutputFile + '''</parameter>
  </processor>

</marlin>'''

	xmlFile = open(xmlFileName , 'w')
	xmlFile.write(xml)
	xmlFile.close()

	os.system('Marlin ' + xmlFileName)
	os.system('rm ' + xmlFileName)
	os.system('mv ' + tempOutputFile + ' ' + a.outputFileName)



if __name__ == '__main__' :

	parser = argparse.ArgumentParser()
	parser.add_argument('runNumber' , help='The run number to process' , type=int)
	parser.add_argument('-e' , '--energy' , help='Energy of the run' , type=float , default=0)
	parser.add_argument('-g' , '--geom' , help='Geometry file' , type=str , default='')
	args = parser.parse_args()

	runNumber = str(args.runNumber) 

	dir = '/home/garillot/files/DATA/STREAMOUT'

	print ('Searching files in ' + dir)

	#list files
	fileList = []

	for fileName in os.listdir(dir) :
		if runNumber in fileName :
			fileList.append(dir + '/' + fileName)

	print 'File List :'
	print fileList

	os.environ["MARLIN"] = '/home/garillot/ilcsoft/v01-19-03/Marlin/v01-12'
	os.environ["PATH"] = '/home/garillot/ilcsoft/v01-19-03/Marlin/v01-12/bin' + ':' + os.environ["PATH"]
	os.environ["MARLIN_DLL"] = '/home/garillot/Trivent/lib/libTrivent.so'

	a = Params()
	a.energy = args.energy
	a.geometry = args.geom
	a.cerenkovDelay = 6

	a.outputFileName = '/home/garillot/files/DATA/TDHCAL_' + str(sys.argv[1]) + '.slcio'

	launch(a , fileList)

