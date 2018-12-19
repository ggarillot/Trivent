#include "FromStreamoutToRoot.h"

#include <EVENT/LCCollection.h>

FromStreamoutToRootProcessor aFromStreamoutToRootProcessor ;

FromStreamoutToRootProcessor::FromStreamoutToRootProcessor()
	: Processor("FromStreamoutToRootProcessor")
{
	// collection
	std::vector<std::string> hcalCollections ;
	hcalCollections.push_back(std::string("DHCALRawHits")) ;
	registerInputCollections( LCIO::RAWCALORIMETERHIT ,
							  "HCALCollections"       ,
							  "HCAL Collection Names" ,
							  _hcalCollections        ,
							  hcalCollections         ) ;


	registerProcessorParameter("RootFileName" ,
							   "ROOT file name" ,
							   rootFileName ,
							   std::string("root.root") ) ;


	//maping on JSON file
	registerProcessorParameter("geometry" ,
							   "JSON geometry" ,
							   geometryFile ,
							   std::string("") ) ;


	registerProcessorParameter( "DataFormat" ,
								"Data Format string: it could be M:3,S-1:3,I:9,J:9,K-1:6 (ILD_ENDCAP) or I:9,J:9,K-1:6,Dif_id:8,Asic_id:6,Chan_id:7",
								_outputFormat,
								std::string("M:3,S-1:3,I:9,J:9,K-1:6"));
}

void FromStreamoutToRootProcessor::init()
{
	printParameters() ;

	assert ( geometryFile != std::string("") ) ;

	mapping = mapping::Mapping(geometryFile) ;
	mapping.print() ;

	evtnum = 0 ;

	rootFile = new TFile(rootFileName.c_str() , "RECREATE") ;
	tree = new TTree("tree" , "tree") ;

	tree->Branch("trigger" , &trigger) ;
	tree->Branch("clockInTrigger" , &clockInTrigger) ;
	tree->Branch("I" , &I) ;
	tree->Branch("J" , &J) ;
	tree->Branch("K" , &K) ;
	tree->Branch("thr" , &thr) ;
}

void FromStreamoutToRootProcessor::processEvent( LCEvent* evtP )
{
	triggerHitMap.clear() ;

	if ( evtP == nullptr )
		return ;

	try
	{
		//!loop over collection
		for ( unsigned int i = 0 ; i < _hcalCollections.size() ; i++ )
		{
			try
			{
				LCCollection* col = nullptr ;
				col = evtP->getCollection( _hcalCollections[i].c_str() ) ;

				streamlog_out( MESSAGE ) << "Trigger number == " << evtP->getEventNumber() << std::endl ;

				if( col == nullptr )
					break ;

				trigger = static_cast<unsigned int> ( evtP->getEventNumber() ) ;

				for ( int ihit = 0 ; ihit < col->getNumberOfElements() ; ++ihit )
				{
					auto rawHit = dynamic_cast<RawCalorimeterHit*>( col->getElementAt(ihit) ) ;
					auto clock = static_cast<unsigned int>( rawHit->getTimeStamp() ) ;

					triggerHitMap[clock].push_back(rawHit) ;
				}
			}
			catch ( lcio::DataNotAvailableException zero )
			{}


			std::cout << triggerHitMap.size() << std::endl ;


			for ( const auto clock : triggerHitMap )
			{
				clockInTrigger = clock.first ;

				for ( const auto hit : clock.second )
				{
					auto difID = mapping::getDifID( hit->getCellID0() ) ;
					auto asicID = mapping::getAsicID( hit->getCellID0() ) ;
					auto padID = mapping::getPadID( hit->getCellID0() ) ;

					if ( !mapping.isListedDif(difID) )
						continue ;

					auto padIndex = mapping.getPadIndex(difID, asicID, padID) ;

					I = padIndex.at(0) ;
					J =  padIndex.at(1) ;
					K = padIndex.at(2) ;

					//fix bug threshold 2 <-> 1
					if ( (hit->getAmplitude()&3) > 2.5f )
						thr = 3 ;
					else if ( (hit->getAmplitude()&3) > 1.5f )
						thr = 1 ;
					else
						thr = 2 ;

					tree->Fill() ;
				}
			}
		}
	}
	catch ( lcio::DataNotAvailableException err )
	{}
}

void FromStreamoutToRootProcessor::end()
{
	rootFile->cd() ;
	tree->Write() ;
	rootFile->Purge() ;
	rootFile->Close() ;
}
