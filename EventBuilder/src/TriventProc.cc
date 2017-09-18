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
#include <algorithm>
#include <exception>

TriventProc a_TriventProc_instance ;

//=========================================================
TriventProc::TriventProc()
	: Processor("TriventProc")
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
							  hcalCollections         ) ;

	// Option of output file with clean events
	_outFileName="LCIO_clean_run.slcio";
	registerProcessorParameter("LCIOOutputFile" ,
							   "LCIO file" ,
							   _outFileName ,
							   _outFileName) ;

	registerProcessorParameter("beamEnergy" ,
							   "The beam ",
							   _beamEnergy ,
							   0.0f) ;

	registerProcessorParameter("cerenkovBif" ,
							   "dif for cerenkov hits" ,
							   _cerenkovBifForMarlin ,
							   3) ;

	cerenkovBif = static_cast<unsigned int>(_cerenkovBifForMarlin) ;

	registerProcessorParameter("cerenkovDelay" ,
							   "cerenkovDelay",
							   cerenkovDelay ,
							   0) ;



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
	registerProcessorParameter("noiseCut" ,
							   "noise cut in time spectrum 10 in default",
							   _noiseCutForMarlin ,
							   10) ;

	_noiseCut = static_cast<unsigned int>(_noiseCutForMarlin) ;

	// time windows
	registerProcessorParameter("timeWin" ,
							   "time window = 2 in default",
							   _timeWinForMarlin ,
							   2) ;

	_timeWin = static_cast<unsigned int>(_timeWinForMarlin) ;


	//maping on XML file
	registerProcessorParameter("setup_geometry" ,
							   "Dif geometry and position on the detector XML" ,
							   _geomXML ,
							   std::string("") ) ;


	// electronic noise cut
	registerProcessorParameter("electronic_noise_cut" ,
							   "number of hit max on time stamp",
							   _elec_noise_cut,
							   500000) ;

	registerProcessorParameter("_time2prev_event_cut" ,
							   "cut on time to previous event (x 200 ns)",
							   _time2prev_event_cut,
							   0) ;


	registerProcessorParameter( "DataFormat" ,
								"Data Format string: it could be M:3,S-1:3,I:9,J:9,K-1:6 (ILD_ENDCAP) or I:9,J:9,K-1:6,Dif_id:8,Asic_id:6,Chan_id:7",
								_outputFormat,
								std::string("M:3,S-1:3,I:9,J:9,K-1:6"));
}

//===============================================
void TriventProc::init()
{
	printParameters() ;

	if ( _geomXML == std::string("") )
		throw std::invalid_argument("wtf put a geometry plz") ;

	trig_count = 0 ;

	_lcWriter = LCFactory::getInstance()->createLCWriter() ;
	_lcWriter->setCompressionLevel( 0 ) ;
	_lcWriter->open(_outFileName.c_str(),LCIO::WRITE_NEW) ;

	XMLReader(_geomXML.c_str()) ;
	printDifGeom() ;
	evtnum=0;// event number

}

//===============================================
void TriventProc::XMLReader(std::string xmlfile)
{
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

			{
				std::string value = pElem->GetText() ;
				std::vector<std::string> lines;
				istringstream iss(value);
				copy(istream_iterator<string>(iss),
					 istream_iterator<string>(),
					 back_inserter<vector<string> >(lines));
				for(unsigned int iline = 0; iline < lines.size(); iline++)
				{
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
					}
				}
			}
		}
	}
	else
	{
		streamlog_out( WARNING ) << red << "Failed to load file : " << xmlfile.c_str() << normal <<std::endl;
	}
}

void TriventProc::printDifGeom()
{

	for(std::map<int,LayerID>::iterator itt = _mapping.begin();itt!=_mapping.end();itt++)
	{
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


void TriventProc::computeTimeSpectrum()
{
	timeSpectrum.clear() ;
	timeSpectrum.assign(_maxtime + 1 , 0) ;
	for ( const auto& it : triggerHitMap )
		timeSpectrum.at(it.first) += it.second.size() ;
}

bool TriventProc::isLocalPeak(unsigned int bin)
{
	int min = std::max( static_cast<int>(bin - _timeWin) , 0 ) ;
	int max = std::max( static_cast<int>(bin + _timeWin) , 0 ) ;

	for ( unsigned int i = static_cast<unsigned int>(min) ; i <= static_cast<unsigned int>(max) ; ++i )
	{
		if ( i == bin || i > _maxtime )
			continue ;

		if (timeSpectrum.at(bin) <= timeSpectrum.at(i) )
			return false ;
	}
	return true ;
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

bool TriventProc::eventBuilder(LCCollection* col_event , unsigned int time_peak , unsigned int prev_time_peak)
{
	zcut.clear() ;
	col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_LONG));
	col_event->setFlag(col_event->getFlag()|( 1 << LCIO::RCHBIT_TIME));
	CellIDEncoder<CalorimeterHitImpl> cd( _outputFormat.c_str() ,col_event) ;
	std::map<int,int> asicMap;

	std::map<int , unsigned int> ramFullDetectorMap ;

	try
	{
		std::vector<int> hitKeys ;

		int min = std::max( static_cast<int>(time_peak - _timeWin) , 0 ) ;
		int max = std::max( static_cast<int>(time_peak + _timeWin) , 0 ) ;

		for ( unsigned int i = static_cast<unsigned int>(min) ; i <= static_cast<unsigned int>(max) ; ++i )
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

int TriventProc::findTheBifSignal(unsigned int timeStamp)
{
	std::map<unsigned int , std::vector<EVENT::RawCalorimeterHit*>>::iterator it = triggerHitMap.end() ;

	int min = std::max( static_cast<int>(timeStamp) - cerenkovDelay - 1 , 0 ) ;
	int max = std::max( static_cast<int>(timeStamp) - cerenkovDelay + 1 , 0 ) ;

	for ( unsigned int i = static_cast<unsigned int>(min) ; i <= static_cast<unsigned int>(max) ; ++i )
	{
		it = triggerHitMap.find(i) ;

		if ( it != triggerHitMap.end() )
			break ;
	}

	if ( it == triggerHitMap.end() )
		return 0 ;

	for ( auto& jt : it->second )
	{
		if ( getCellDif_id( jt->getCellID0() ) == cerenkovBif )
		{
			//						std::cout << "Dif " << getCellDif_id( jt->getCellID0() )
			//								  << " , Asic " << getCellAsic_id( jt->getCellID0() )
			//								  << " , Pad " << getCellChan_id( jt->getCellID0())
			//								  << " , Thr " << jt->getAmplitude() << std::endl ;
			_cerenkovTime = it->first - timeStamp ;
			return jt->getAmplitude() ;
		}
	}
	return 0 ;
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
					triggerHitMap.clear() ;
					std::vector<int> vTrigger;
					for (int ihit(0); ihit < col->getNumberOfElements(); ++ihit)
					{// loop over the hits
						RawCalorimeterHit* raw_hit = dynamic_cast<RawCalorimeterHit*>( col->getElementAt(ihit) ) ;
						if (NULL != raw_hit)
						{
							if (raw_hit->getTimeStamp() < 0 )
								continue ;
							triggerHitMap[static_cast<unsigned int>(raw_hit->getTimeStamp())].push_back(raw_hit) ;
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

					getMaxTime() ;
					computeTimeSpectrum() ;

					unsigned int ibin = 0 ;
					int bin_c_prev = -2 * _timeWin ; //  the previous bin center

					int time_prev = 0;
					while( ibin < (_maxtime - 2) )
					{
						if( timeSpectrum.at(ibin) <= _noiseCut || !isLocalPeak(ibin) )
						{
							ibin++ ;
							continue ;
						}

						LCEventImpl* evt = new LCEventImpl() ;     // create the event

						//---------- set event paramters ------
						evt->parameters().setValue("trigger" , evtP->getEventNumber()) ;
						if ( _beamEnergy > 0 )
							evt->parameters().setValue("ParticleEnergy" , _beamEnergy) ;
						evt->parameters().setValue("bcid1" , _bcid1) ;
						evt->parameters().setValue("bcid2" , _bcid2) ;
						evt->parameters().setValue("eventTimeInTrigger" , static_cast<int>(timeSpectrum.at(ibin)) ) ;
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

							evt->parameters().setValue("cerenkovTag" , findTheBifSignal(ibin) ) ;
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
						delete evt ;
						evt = nullptr ;

						bin_c_prev = ibin;
						ibin = ibin + _timeWin ;
					}

				}
				catch (lcio::DataNotAvailableException zero) {}
				catch ( std::out_of_range& e) { std::cout << e.what() << std::endl ; }

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
	_lcWriter->close() ;
}
//==============================================================
