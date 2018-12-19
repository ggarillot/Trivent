#ifndef TriventProc_h
#define TriventProc_h

#define  HISTOGRAM_PARSER true

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/RawCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <TTree.h>
#include <TFile.h>

#include "IO/LCWriter.h"
#include <map>
#include <algorithm>
#include <set>
#include "Mapping.h"


class TriventProc : public marlin::Processor
{
	public :

		Processor*  newProcessor() { return new TriventProc ; }

		TriventProc() ;

		~TriventProc() {}

		void init() ;

		void processEvent( LCEvent * evtP );


		void getMaxTime() ;
		void computeTimeSpectrum() ;
		bool isLocalPeak(unsigned int bin) ;

		bool eventBuilder(LCCollection* col_event, unsigned int time_peak, unsigned int prev_time_peak) ;
		int	findTheBifSignal(unsigned int timeStamp) ;

		void end() ;

		TriventProc(const TriventProc&) = delete ;
		void operator=(const TriventProc&) = delete ;


	protected :

		std::string rootFileName = "" ;
		TFile* rootFile = nullptr ;
		TTree* tree = nullptr ;

		unsigned int nHit = 0 ;
		unsigned int trigger = 0 ;
		unsigned int clock = 0 ;

		unsigned long long triggerBeginTime = 0 ;
		unsigned long long timeInTrigger = 0 ;

		unsigned long long timeOfFirstTrigger = 0 ;

		bool isFirstEvent = true ;
		//		double acquisitionTime = 0 ;

		std::map<unsigned int , std::vector<EVENT::RawCalorimeterHit*>> triggerHitMap {{}} ;
		std::vector<unsigned int> timeSpectrum {} ;

		std::string _outFileName = "" ;

		std::string _colName = "" ;
		std::string _fileName = "" ;

		std::string geometryFile = "" ;

		mapping::Mapping mapping{} ;

		std::vector<std::string> _hcalCollections {} ;

		Int_t _elec_noise_cut = 500000 ;


		bool removeSquareEvents = true ;
		bool removeRamFullEvents = true ;
		int maxHitPerLayer = 400 ;

		unsigned int cerenkovBif = 3 ;
		int cerenkovDelay = 0 ;

		unsigned int _noiseCut = 10 ;
		unsigned int _timeWin = 2 ;

		int _cerenkovBifForMarlin = 3 ;
		int _noiseCutForMarlin = 10 ;
		int _timeWinForMarlin = 2 ;

		int _LayerCut = 5 ;
		int _time2prev_event_cut = 0 ;
		float _beamEnergy = 0 ;
		unsigned int _maxtime = 0 ;
		int evtnum  = 0 ;
		int _rejectedNum  = 0 ;
		int _index[3] = {0,0,0} ;
		std::set<int> zcut {} ;
		LCWriter* _lcWriter = nullptr ;
		int _bcid1 = 0 ;
		int _bcid2 = 0 ;

		std::string _outputFormat = "" ;

		std::string normal  = {0x1b,'[','0',';','3','9','m',0} ;
		std::string red     = {0x1b,'[','1',';','3','1','m',0} ;
		std::string green   = {0x1b,'[','1',';','3','2','m',0} ;
		std::string yellow  = {0x1b,'[','1',';','3','3','m',0} ;
		std::string blue    = {0x1b,'[','1',';','3','4','m',0} ;
		std::string magenta = {0x1b,'[','1',';','3','5','m',0} ;
		std::string white   = {0x1b,'[','1',';','3','9','m',0} ;

} ;
#endif //TriventProc_h


