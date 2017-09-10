/**
* Yacine HADDAD
* LLR Ecole polytechnique
* avril 2012
* Trivent v0.3
*/

#include <TriventProc.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <limits.h>
#include <cmath>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>
#include <Rtypes.h>
#include <sstream>
#include <UTIL/CellIDEncoder.h>
#include "Mapping.h"
#include "TObject.h"
#include "TRefArray.h"
#include "TRef.h"
#include <fstream>
#include <algorithm>
#include <lcio.h>
#include "marlin/VerbosityLevels.h"
#include "marlin/tinyxml.h"

#include <limits>

TriventProc a_TriventProc_instance ;

//=========================================================
TriventProc::TriventProc()
	: Processor("TriventProc"),
	  _output(0),
	  _outputTree(0)//,
	//_hitArray(0)
{

	streamlog_out( MESSAGE )<< "Trivent ... begin " << endl;
	_rejectedNum = 0;

	// collection
	std::vector<std::string> hcalCollections;
	hcalCollections.push_back(std::string("DHCALRawHits"));
	registerInputCollections( LCIO::RAWCALORIMETERHIT ,
							  "HCALCollections"       ,
							  "HCAL Collection Names" ,
							  _hcalCollections        ,
							  hcalCollections         );

	// Option of output file with clean events
	_outFileName="LCIO_clean_run.slcio";
	registerProcessorParameter("LCIOOutputFile" ,
							   "LCIO file" ,
							   _outFileName ,
							   _outFileName);
	// Energy
	_beamEnergy = 10;
	registerProcessorParameter("beamEnergy" ,
							   "The beam ",
							   _beamEnergy ,
							   _beamEnergy);

	registerProcessorParameter("cerenkovBif" ,
							   "dif for cerenkov hits" ,
							   cerenkovBif ,
							   3) ;

	registerProcessorParameter("cerenkovDelay" ,
							   "cerenkovDelay",
							   cerenkovDelay ,
							   0) ;


	// Option of output file with noise
	_noiseFileName="noise_run.slcio";
	registerProcessorParameter("NOISEutputFile" ,
							   "NOISE file" ,
							   _noiseFileName ,
							   _noiseFileName);


	registerProcessorParameter("RemoveSquareEvents" ,
							   "Remove events with an ASIC with 64hits" ,
							   removeSquareEvents ,
							   true) ;

	registerProcessorParameter("RemoveRamFullEvents" ,
							   "RemoveRamFullEvents" ,
							   removeRamFullEvents ,
							   true) ;

	// layer cut
	_LayerCut = 10;
	registerProcessorParameter("LayerCut" ,
							   "cut in number of layer 10 in default",
							   _LayerCut ,
							   _LayerCut);


	// noise cut
	_noiseCut = 10;
	registerProcessorParameter("noiseCut" ,
							   "noise cut in time spectrum 10 in default",
							   _noiseCut ,
							   _noiseCut);

	// time windows
	_timeWin = 2;
	registerProcessorParameter("timeWin" ,
							   "time window = 2 in default",
							   _timeWin ,
							   _timeWin);
	//maping on XML file
	_geomXML = "setup_geometry.xml";
	registerProcessorParameter("setup_geometry" ,
							   "Dif geometry and position on the detector XML",
							   _geomXML,
							   _geomXML);

	//maping on txt file
	_mappingfile = "mapping_ps.txt";
	registerProcessorParameter("DIFMapping" ,
							   "dif's mapping file ",
							   _mappingfile,
							   _mappingfile);

	// electronic noise cut
	_elec_noise_cut = 5000;
	registerProcessorParameter("electronic_noise_cut" ,
							   "number of hit max on time stamp",
							   _elec_noise_cut,
							   _elec_noise_cut);

	// electronic noise cut
	_time2prev_event_cut = 0;
	registerProcessorParameter("_time2prev_event_cut" ,
							   "cut on time to previous event (x 200 ns)",
							   _time2prev_event_cut,
							   _time2prev_event_cut);


	//log root file
	_treeName = "TEST";
	registerProcessorParameter("TreeName_logroot" ,
							   "Logroot tree name",
							   _treeName,
							   _treeName);
	// histogram control tree
	_logrootName = "logroot.root";
	registerProcessorParameter("logroot_Name" ,
							   "Logroot name",
							   _logrootName,
							   _logrootName);

	GAIN_CORRECTION_MODE = false;
	registerProcessorParameter("GAIN_CORRECTION_MODE",
							   "GAIN_CORRECTION_MODE",
							   GAIN_CORRECTION_MODE,
							   GAIN_CORRECTION_MODE);

	registerProcessorParameter( "DataFormat" ,
								"Data Format string: it could be M:3,S-1:3,I:9,J:9,K-1:6 (ILD_ENDCAP) or I:9,J:9,K-1:6,Dif_id:8,Asic_id:6,Chan_id:7",
								_outputFormat,
								std::string("M:3,S-1:3,I:9,J:9,K-1:6"));
}

void TriventProc::XMLReader(std::string xmlfile){
	TiXmlDocument doc(xmlfile.c_str());
	bool load_key = doc.LoadFile();
	if(load_key){
		streamlog_out( MESSAGE ) << green << "File : " << xmlfile.c_str() << normal <<std::endl;
		// tout ici
		TiXmlHandle hDoc(&doc);
		TiXmlElement* pElem;
		TiXmlHandle hRoot(0);
		// name block
		{
			pElem=hDoc.FirstChildElement().Element();
			// should always have a valid root but handle gracefully if it does
			if (!pElem) streamlog_out( WARNING ) << red << "error elem" << normal << std::endl;
			streamlog_out( MESSAGE ) << green << pElem->Value() << normal << std::endl;

			// save this for later
			hRoot=TiXmlHandle(pElem);
		}
		// parameters block
		{
			m_parameters.clear();
			pElem=hRoot.FirstChild("parameter").Element();
			std::string key = pElem->Attribute("name");
			streamlog_out( MESSAGE ) << green << key.c_str() << normal << std::endl;
			streamlog_out( DEBUG1 ) << green
									<<"parameter : "
								   << pElem->Attribute("name")
								   << normal
								   << std::endl;

			std::vector<std::string> lines;
			{
				std::string value = pElem->GetText() ;
				std::vector<std::string> lines;
				istringstream iss(value);
				copy(istream_iterator<string>(iss),
					 istream_iterator<string>(),
					 back_inserter<vector<string> >(lines));
				for(unsigned int iline = 0; iline < lines.size(); iline++){
					std::string line = lines.at(iline);
					streamlog_out( MESSAGE ) << red << line << normal << std::endl;

					stringstream ss( line.c_str() );
					vector<string> result;

					LayerID mapp;
					int Dif_id;
					while( ss.good() )
					{
						string substr;
						getline( ss, substr, ',' );
						result.push_back( substr );
					}
					istringstream ( result.at(0) ) >> Dif_id;
					istringstream ( result.at(1) ) >> mapp.K;
					istringstream ( result.at(2) ) >> mapp.DifX;
					istringstream ( result.at(3) ) >> mapp.DifY;
					istringstream ( result.at(4) ) >> mapp.IncX;
					istringstream ( result.at(5) ) >> mapp.IncY;
					_mapping[Dif_id] = mapp;
				}
			}
			pElem = pElem->NextSiblingElement();
			// ChamberGeom  Node.
			{
				streamlog_out( DEBUG1 ) << green
										<<"parameter : "
									   << pElem->Attribute("name")
									   << normal
									   << std::endl;
				std::vector<std::string> lines;
				{
					std::string value = pElem->GetText() ;
					std::vector<std::string> lines;
					istringstream iss(value);
					copy(istream_iterator<string>(iss),
						 istream_iterator<string>(),
						 back_inserter<vector<string> >(lines));
					for(unsigned int iline = 0; iline < lines.size(); iline++){
						std::string line = lines.at(iline);
						streamlog_out( MESSAGE ) << red << line << normal << std::endl;

						stringstream ss( line.c_str() );
						vector<string> result;

						double position;
						int Dif_id;
						while( ss.good() )
						{
							string substr;
							getline( ss, substr, ',' );
							result.push_back( substr );
						}
						istringstream ( result.at(0) ) >> Dif_id;
						istringstream ( result.at(3) ) >> position;

						_chamber_pos[Dif_id] = position;
					}
				}
			}
		}
	}else{
		streamlog_out( WARNING ) << red << "Failed to load file : " << xmlfile.c_str() << normal <<std::endl;
	}
}

void TriventProc::readDifGeomFile(std::string geomfile){

	cout << "read the mapping file .."<< endl;

	LayerID contenu;
	ifstream file(geomfile.c_str(), ios::in);
	if(file){
		while(!file.eof()){
			int Dif_id;
			char co;
			file >> Dif_id >> co
					>> contenu.K >> co
					>> contenu.DifX >> co
					>> contenu.DifY >> co
					>> contenu.IncX >> co
					>> contenu.IncY ;
			_mapping [Dif_id] = contenu;
		}
		file.close();
	}
	else
		cerr << "ERROR ... maping file not correct !" << endl;
}

void TriventProc::printDifGeom(){

	for(std::map<int,LayerID>::iterator itt = _mapping.begin();itt!=_mapping.end();itt++)     {
		streamlog_out( MESSAGE ) << itt->first << "\t" << itt->second.K
								 <<"\t"<<itt->second.DifX
								<<"\t"<<itt->second.DifY
							   <<"\t"<<itt->second.IncX
							  <<"\t"<<itt->second.IncY
							 << std::endl;
	}
}

// ============ decode the cell ids =============
uint TriventProc::getCellDif_id(int cell_id){
	return cell_id & 0xFF;
}
uint TriventProc::getCellAsic_id(int cell_id){
	return (cell_id & 0xFF00)>>8;
}
uint TriventProc::getCellChan_id(int cell_id){
	return (cell_id & 0x3F0000)>>16;
}

uint* TriventProc::getPadIndex(uint dif_id, uint asic_id, uint chan_id)
{
	_index[0]=_index[1]=_index[2]=0;
	double DifY = -1.,DifZ = -1.;
	DifZ = _mapping.find(dif_id)->second.K;
	DifY = _mapping.find(dif_id)->second.DifY;
	_index[0] = (1+MapILargeHR2[chan_id]+AsicShiftI[asic_id]);
	_index[1] = (32-(MapJLargeHR2[chan_id]+AsicShiftJ[asic_id]))+int(DifY);
	_index[2] = abs(int(DifZ));
	streamlog_out( DEBUG0 ) << " Dif_id == " << dif_id
							<< " Asic_id ==" << asic_id
							<< " Chan_id ==" << chan_id
							<< " I == " << _index[0]
							<< " J == " << _index[1]
							<< " K == " << _index[2]
							<< std::endl;
	return _index;
}
//===============================================
void TriventProc::getMaxTime()
{
	_maxtime = 0 ;
	if ( !triggerHitMap.empty() )
		_maxtime = triggerHitMap.rbegin()->first ;
}


std::vector<int> TriventProc::getTimeSpectrum()
{
	std::vector<int> time_spectrum(_maxtime + 1) ;
	for ( const auto& it : triggerHitMap )
	{
		if ( it.first >= 0 )
			time_spectrum.at(it.first) += it.second.size() ;
	}

	return time_spectrum ;
}

int IJKToKey(const int i,const int j,const int k){return 100*100*k+100*j+i;}
int findAsicKey(int i,int j,int k)
{
	if(i>96||i<0||j>96||j<0) return -1;
	int jnum=(j-1)/8;
	int inum=(i-1)/8;
	int num=jnum*12+inum;
	return k*1000+num;
}

bool TriventProc::eventBuilder(LCCollection* col_event , int time_peak , int prev_time_peak)
{
	//	std::cout << "eventBuilder" << std::endl ;
	zcut.clear();
	col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_LONG));
	col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_TIME));
	CellIDEncoder<CalorimeterHitImpl> cd( _outputFormat.c_str() ,col_event) ;
	std::map<int,int> asicMap;

	std::map<int , unsigned int> ramFullDetectorMap ;

	try
	{
		std::vector<int> hitKeys ;
		for ( int i = time_peak - _timeWin ; i <= time_peak + _timeWin ; ++i )
		{
			if ( triggerHitMap.find(i) == triggerHitMap.end() )
				continue ;
			if ( !(i > prev_time_peak + _timeWin) )
				continue ;

			for ( const auto& rawhit : triggerHitMap.at(i) )
			{
				int Dif_id  =  getCellDif_id ( rawhit->getCellID0() ) ;
				int Asic_id =  getCellAsic_id( rawhit->getCellID0() ) ;
				int Chan_id =  getCellChan_id( rawhit->getCellID0() ) ;

				if ( _mapping.find(Dif_id) == _mapping.end() )
				{
					//std::cout << "wtf dif : " << Dif_id << std::endl ;
					continue ;
				}

				if ( Chan_id == 29 || Chan_id == 31 )
					ramFullDetectorMap[Dif_id] ++ ;

				if ( removeRamFullEvents && ramFullDetectorMap[Dif_id] > 70 )
				{
					//					std::cout << "Reject ramFull event of dif : " << Dif_id << std::endl ;
					return false ;
				}

				int I = getPadIndex(Dif_id, Asic_id, Chan_id)[0];
				int J = getPadIndex(Dif_id, Asic_id, Chan_id)[1];
				int K = getPadIndex(Dif_id, Asic_id, Chan_id)[2];

				//find ans remove square events
				int asickey = findAsicKey(I,J,K);
				if(asicMap[asickey])
					asicMap[asickey]++;
				else
					asicMap[asickey]=1;
				if( removeSquareEvents && asicMap[asickey] == 64 )
				{
					zcut.clear();
					hitKeys.clear();
					asicMap.clear();
					//					std::cout << "fullAsic of asicKey : " << asickey << std::endl ;
					return false ;
				}

				int aHitKey = IJKToKey(I,J,K) ;

				float pos[3] ;
				pos[0] = I*10.408f ;
				pos[1] = J*10.408f ;
				pos[2] = K*26.131f ;

				CalorimeterHitImpl* caloHit = new CalorimeterHitImpl();
				caloHit->setTime( float(rawhit->getTimeStamp()) ) ; // done !!

				if(float(rawhit->getAmplitude()&3)>2.5)
					caloHit->setEnergy(float(rawhit->getAmplitude()&3));
				else if(float(rawhit->getAmplitude()&3)>1.5)
					caloHit->setEnergy(float(rawhit->getAmplitude()&3)-1);
				else
					caloHit->setEnergy(float(rawhit->getAmplitude()&3)+1);

				//avoid two hit in the same cell
				if(std::find(hitKeys.begin(),hitKeys.end(),aHitKey)!=hitKeys.end())
				{
					IMPL::CalorimeterHitImpl* hit = dynamic_cast<IMPL::CalorimeterHitImpl*>(col_event->getElementAt(std::distance(hitKeys.begin(),std::find(hitKeys.begin(),hitKeys.end(),aHitKey))));
					float hitTime=hit->getTime();
					if( std::abs(time_peak-hitTime) > std::abs(time_peak-i) )
						hit->setEnergy(caloHit->getEnergy());

					continue ;
				}
				// set the cell id
				if( _outputFormat==std::string("M:3,S-1:3,I:9,J:9,K-1:6") ){
					cd["I"] = I ;
					cd["J"] = J ;
					cd["K-1"] = K-1 ;
					cd["M"] = 0 ;
					cd["S-1"] = 3 ;
				}
				else if( _outputFormat==std::string("I:9,J:9,K-1:6,Dif_id:8,Asic_id:6,Chan_id:7") ){
					cd["I"] = I ;
					cd["J"] = J ;
					cd["K-1"] = K-1 ;
					cd["Dif_id"] = Dif_id ;
					cd["Asic_id"] = Asic_id ;
					cd["Chan_id"] = Chan_id ;
				}
				else{
					streamlog_out( ERROR ) << "WRONG OUTPUT DATA FORMAT -> THROW" << std::endl;
					throw;
				}
				streamlog_out( DEBUG0 ) << " I == " << I
										<< " J == " << J
										<< " K == " << K
										<< std::endl;
				cd.setCellID( caloHit ) ;
				zcut.insert(K) ;

				caloHit->setPosition(pos);
				col_event->addElement(caloHit);
				hitKeys.push_back(aHitKey);

			}//loop over the hit
		}
		hitKeys.clear() ;
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << " collection not available" << std::endl;
	}
	catch(const std::out_of_range& e)
	{
		std::cout << e.what() << std::endl ;
	}

	return true ;
}

bool TriventProc::findTheBifSignal(int timeStamp)
{
//	int time_diff = std::numeric_limits<int>::max() ;

	std::map<int , std::vector<EVENT::RawCalorimeterHit*>>::iterator it ;
	for ( int i = timeStamp - cerenkovDelay - 1 ; i <= timeStamp - cerenkovDelay + 1 ; ++i )
	{
		it = triggerHitMap.find(i) ;

		if ( it != triggerHitMap.end() )
			break ;
	}

	if ( it == triggerHitMap.end() )
		return false ;


	for ( auto& jt : it->second )
	{
		if ( getCellDif_id( jt->getCellID0() ) == cerenkovBif )
		{
			//			std::cout << "Dif " << getCellDif_id( jt->getCellID0() )
			//					  << " , Asic " << getCellAsic_id( jt->getCellID0() )
			//					  << " , Pad " << getCellChan_id( jt->getCellID0())
			//					  << " , Thr " << jt->getAmplitude() << std::endl ;
			_cerenkovTime = it->first - timeStamp ;
			return true ;
		}
	}

	return false ;

	//	for(std::vector<EVENT::RawCalorimeterHit*>::iterator rawhit = _trigger_raw_hit.begin() ; rawhit != _trigger_raw_hit.end() ; rawhit++)
	//	{
	//		if( getCellDif_id((*rawhit)->getCellID0())==3 && fabs(timeStamp-(*rawhit)->getTimeStamp())<time_diff /*&& time_peak-(*rawhit)->getTimeStamp()>0*/ )
	//			time_diff=fabs(timeStamp-(*rawhit)->getTimeStamp());
	//	}
	//	if( time_diff<cerenkovDelay ) {
	//		_cerenkovTime=time_diff;
	//		return true;
	//	}
	//	else return false;
}
//===============================================
void TriventProc::init() {
	trig_count = 0;
	//========================
	//readDifGeomFile(_mappingfile.c_str());

	// ========================

	printParameters();
	// new process

	char cnormal[8] =  {0x1b,'[','0',';','3','9','m',0};
	char cred[8]     = {0x1b,'[','1',';','3','1','m',0};
	char cgreen[8]   = {0x1b,'[','1',';','3','2','m',0};
	char cyellow[8]  = {0x1b,'[','1',';','3','3','m',0};
	char cblue[8]    = {0x1b,'[','1',';','3','4','m',0};
	char cmagenta[8] = {0x1b,'[','1',';','3','5','m',0};
	char cwhite[8]   = {0x1b,'[','1',';','3','9','m',0};

	normal   = cnormal;
	red      = cred;
	green    = cgreen;
	yellow   = cyellow;
	blue     = cblue;
	magenta  = cmagenta;
	white    = cwhite;

	_lcWriter = LCFactory::getInstance()->createLCWriter() ;
	_lcWriter->setCompressionLevel( 0 ) ;
	_lcWriter->open(_outFileName.c_str(),LCIO::WRITE_NEW) ;


	XMLReader(_geomXML.c_str());
	printDifGeom();
	evtnum=0;// event number

}

//==================================================================================
void TriventProc::processEvent( LCEvent* evtP )
{
	if (evtP != nullptr)
	{
		try
		{
			_eventNr = evtP->getEventNumber() ;
			for(unsigned int i=0; i< _hcalCollections.size(); i++)
			{//!loop over collection

				unsigned int currentTriggerNEvents = 0 ;
				try
				{
					LCCollection* col = nullptr ;
					col = evtP ->getCollection(_hcalCollections[i].c_str()) ;
					int numElements = col->getNumberOfElements() ;// hit number


					streamlog_out( MESSAGE ) << yellow << "Trigger number == " << trig_count++ << normal << std::endl;

					if( col == nullptr )
					{
						streamlog_out( WARNING )<< red << "TRIGGER SKIPED ..."<< normal <<std::endl;
						break;
					}

					if(numElements > _elec_noise_cut)
					{
						streamlog_out( MESSAGE ) << red << "TRIGGER SKIPED ..."<< normal <<std::endl;
						break;
					}

					// set raw hits
					_trigger_raw_hit.clear() ;
					triggerHitMap.clear() ;
					std::vector<int> vTrigger;
					for (int ihit(0); ihit < col->getNumberOfElements(); ++ihit)
					{// loop over the hits
						RawCalorimeterHit* raw_hit = dynamic_cast<RawCalorimeterHit*>( col->getElementAt(ihit) ) ;
						if (NULL != raw_hit)
						{
							_trigger_raw_hit.push_back(raw_hit) ;
							triggerHitMap[raw_hit->getTimeStamp()].push_back(raw_hit) ;
							//extract abolute bcid information:
							if(ihit==0)
							{
								unsigned int difid=0;
								difid = raw_hit->getCellID0()&0xFF;
								if (difid==0) return;
								std::stringstream pname("");
								pname <<"DIF"<<difid<<"_Triggers";
								col->getParameters().getIntVals(pname.str(),vTrigger);
								if (vTrigger.size()!=0)
								{
									_bcid1=vTrigger[4] ;
									_bcid2=vTrigger[3] ;
									unsigned long long Shift=16777216ULL;//to shift the value from the 24 first bits
									unsigned long long theBCID_=_bcid1*Shift+_bcid2;
									streamlog_out( DEBUG1 ) << "trigger time : " << theBCID_ << std::endl;
								}
							}
						}
					}
					//					auto sortHitsByTime = [](const RawCalorimeterHit* hit1 , const RawCalorimeterHit* hit2) -> bool { return hit1->getTimeStamp() < hit2->getTimeStamp() ; } ;
					//					std::sort( _trigger_raw_hit.begin() , _trigger_raw_hit.end() , sortHitsByTime) ;
					getMaxTime() ;
					std::vector<int> time_spectrum = getTimeSpectrum();

					//---------------------------------------------------------------
					//! Find the condidate event
					int ibin=0;
					int bin_c_prev = -2 * _timeWin ; //  the previous bin center

					int time_prev = 0;
					while(ibin < (_maxtime+1))
					{
						if(time_spectrum[ibin] >= _noiseCut &&
						   time_spectrum[ibin] >= time_spectrum[ibin+1] &&
						   time_spectrum[ibin] >= time_spectrum[ibin-1] &&
						   time_spectrum[ibin] >= time_spectrum[ibin-2] &&
						   time_spectrum[ibin] >= time_spectrum[ibin+2] )
						{
							LCEventImpl* evt = new LCEventImpl() ;     // create the event

							//---------- set event paramters ------
							const std::string parname_trigger = "trigger";
							const std::string parname_energy  = "beamEnergy";
							const std::string parname_bcid1 = "bcid1";
							const std::string parname_bcid2 = "bcid2";
							evt->parameters().setValue(parname_trigger,evtP->getEventNumber());
							evt->parameters().setValue(parname_energy , _beamEnergy);
							evt->parameters().setValue(parname_bcid1 , _bcid1);
							evt->parameters().setValue(parname_bcid2 , _bcid2);
							evt->parameters().setValue("eventTimeInTrigger" , time_spectrum[ibin]);
							evt->setRunNumber( evtP->getRunNumber()) ;
							//-------------------------------------

							LCCollectionVec* outcol = new LCCollectionVec(LCIO::CALORIMETERHIT) ;
							bool add = TriventProc::eventBuilder(outcol,ibin,bin_c_prev) ;

							streamlog_out( DEBUG1 ) << "zcut.size() = " << zcut.size() << "\t _LayerCut = " << _LayerCut << std::endl ;

							int maxConsecutiveLayers = 0 ;
							int currentCtr = 1 ;
							int previousOk = std::numeric_limits<int>::min() ;
							bool bonusAllowed = false ;
							for ( std::set<int>::iterator it = zcut.begin() ; it != zcut.end() ; ++it )
							{
								if ( *it == previousOk + 1 )
								{
									currentCtr++ ;
									if ( currentCtr > maxConsecutiveLayers )
										maxConsecutiveLayers = currentCtr ;
								}
								else if ( *it == previousOk + 2 && !bonusAllowed )
								{
									bonusAllowed = true ;
									currentCtr++ ;
									if ( currentCtr > maxConsecutiveLayers )
										maxConsecutiveLayers = currentCtr ;
								}
								else
								{
									currentCtr = 1 ;
									bonusAllowed = false ;
								}

								previousOk = *it ;
								if ( maxConsecutiveLayers > _LayerCut )
									break ;
							}

							if( maxConsecutiveLayers >_LayerCut && abs(int(ibin)-time_prev) > _time2prev_event_cut && add)
							{
								streamlog_out( DEBUG5 ) <<green<<" Trivent find event at :==> "<< red << ibin
													   <<green<<"\t :Nhit: ==> "<< magenta
													  <<outcol->getNumberOfElements() << normal <<std::endl;
								evt->setEventNumber( evtnum++ ) ;

								if( findTheBifSignal(ibin) )
								{
									evt->parameters().setValue("cerenkovTag",true) ;
									evt->parameters().setValue("cerenkovTime",_cerenkovTime) ;
									streamlog_out( DEBUG ) << evtnum << " is tagged by Cerenkov " << ibin << std::endl;
								}
								else
									evt->parameters().setValue("cerenkovTag",false) ;
								evt->addCollection(outcol, "SDHCAL_HIT");
								_lcWriter->writeEvent( evt ) ;
								currentTriggerNEvents++ ;
							}
							else
							{
								_rejectedNum++ ;
								delete outcol ;
							}
							time_prev = ibin ;
							delete evt ; evt = nullptr ;

							bin_c_prev = ibin;
							ibin = ibin + _timeWin ;
						}
						else
							ibin++ ;
					}

				}
				catch (lcio::DataNotAvailableException zero) {}

				std::cout << "nEvents : " << currentTriggerNEvents << std::endl ;
			}
		}
		catch (lcio::DataNotAvailableException err) {}
	}


}
//==============================================================
void TriventProc::end()
{
	streamlog_out( MESSAGE )<< "Trivent Selected "<< evtnum <<" events"<<std::endl;
	streamlog_out( MESSAGE )<< "Trivent Rejected "<< _rejectedNum <<" events"<<std::endl;
	streamlog_out( MESSAGE )<< "Trivent end"<<std::endl;
	//cc.StoreHistos("test.root");
	_lcWriter->close();

	if (_outputTree) {
		TFile *_logroot = _outputTree->GetCurrentFile();
		_logroot->Write();
		delete _logroot;
	}
}
//==============================================================
