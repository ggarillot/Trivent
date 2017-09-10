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
#include <TH1.h>
#include <TH1F.h>
#include "IO/LCWriter.h"
#include <map>
#include <algorithm>
#include <set>
#include "Mapping.h"

using namespace std;

class TriventProc : public marlin::Processor
{
	public :

		Processor*  newProcessor() { return new TriventProc ; }

		TriventProc() ;

		~TriventProc() {}

		void init() ;


		void    processEvent( LCEvent * evtP );
		void    XMLReader(std::string xmlfile);
		void    readDifGeomFile(std::string geomfile);
		void    printDifGeom();

		uint    getCellDif_id(int cell_id);
		uint    getCellAsic_id(int cell_id);
		uint    getCellChan_id(int cell_id);

		void    getMaxTime();
		std::vector<int> getTimeSpectrum();
		uint*   getPadIndex(uint dif_id, uint asic_id, uint chan_id);
		bool    eventBuilder(LCCollection* col_event,int time_peak, int prev_time_peak);
		bool    findTheBifSignal(int timeStamp);

		void    end();

	protected :
		TH1F *noise_dist;
		TH1F *gain_chan;
		TH1F *mean_hit_dif;
		TH1F *time_hit_dif;
		// xml test
		std::map<std::string,std::string> m_parameters;

		std::vector<EVENT::RawCalorimeterHit*> _trigger_raw_hit;

		std::map<int , std::vector<EVENT::RawCalorimeterHit*>> triggerHitMap ;

		bool GAIN_CORRECTION_MODE;
		std::string _outFileName;
		std::string _noiseFileName;
		std::string _treeName;
		std::string _logrootName;
		std::string _colName;
		std::string _fileName;
		std::string _mappingfile;
		std::string _geomXML;
		std::ostream *_output;
		std::vector<std::string> _hcalCollections;
		int _overwrite;
		TTree *_outputTree;
		unsigned int _eventNr;
		Int_t _Nhit;
		Int_t _elec_noise_cut;

		std::map<int, LayerID  > _mapping;
		std::map<int, double  > _chamber_pos;//camber , pos

		bool removeSquareEvents ;
		bool removeRamFullEvents ;

		int cerenkovBif = 0 ;
		int cerenkovDelay = 0 ;

		int _noiseCut;
		int _timeWin;
		int _LayerCut;
		int _time2prev_event_cut;
		int _cerenkovTime;
		float _beamEnergy;
		int trig_count;
		int _maxtime;
		int evtnum;
		int _rejectedNum;
		uint _index[3];
		std::set<int> zcut ;
		//		uintVec zcut;
		LCWriter* _lcWriter;
		int _bcid1;
		int _bcid2;

		std::string _outputFormat;

		std::string normal  ;
		std::string red     ;
		std::string green   ;
		std::string yellow  ;
		std::string blue    ;
		std::string magenta ;
		std::string white   ;

} ;
#endif //TriventProc_h


