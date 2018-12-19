#ifndef FromStreamoutToRoot_h
#define FromStreamoutToRoot_h

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <EVENT/RawCalorimeterHit.h>

#include <TTree.h>
#include <TFile.h>

#include "Mapping.h"


class FromStreamoutToRootProcessor : public marlin::Processor
{
	public :

		Processor* newProcessor() { return new FromStreamoutToRootProcessor ; }

		FromStreamoutToRootProcessor() ;

		~FromStreamoutToRootProcessor() = default ;

		void init() ;

		void processEvent( LCEvent * evtP );

		void end() ;

		FromStreamoutToRootProcessor(const FromStreamoutToRootProcessor&) = delete ;
		void operator=(const FromStreamoutToRootProcessor&) = delete ;


	protected :

		std::string rootFileName = "" ;
		std::vector<std::string> _hcalCollections = {} ;
		std::string geometryFile = "" ;
		std::string _outputFormat = "" ;

		mapping::Mapping mapping{} ;

		std::map<unsigned int , std::vector<EVENT::RawCalorimeterHit*>> triggerHitMap {{}} ;

		TFile* rootFile = nullptr ;
		TTree* tree = nullptr ;

		unsigned int evtnum = 0 ;

		unsigned int trigger = 0 ;
		unsigned long long clockInTrigger = 0 ;

		int I = 0 ;
		int J = 0 ;
		int K = 0 ;
		int thr = 1 ;
} ;


#endif //FromStreamoutToRoot_h
