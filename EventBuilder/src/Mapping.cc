#include "Mapping.h"


unsigned int mapping::getDifID(int cell_id)
{
	return cell_id & 0xFF;
}
unsigned int mapping::getAsicID(int cell_id)
{
	return (cell_id & 0xFF00)>>8;
}
unsigned int mapping::getPadID(int cell_id)
{
	return (cell_id & 0x3F0000)>>16;
}

mapping::Mapping::Mapping(std::string jsonFile)
{
	difList.clear() ;

	std::ifstream file(jsonFile) ;
	auto json = nlohmann::json::parse(file) ;

	auto chambersList = json.at("chambers") ;

	for ( const auto& i : chambersList )
	{
		int slot = i.at("slot") ;
		unsigned int difID = i.at("left") ;

		mapping::Dif temp{ slot , 0 } ;

		//insert the new dif and check if it already exists elsewhere
		auto it = difList.insert( { difID , temp } ) ;
		if ( !it.second )
		{
			std::cout << "ERROR in geometry : dif " << difID << " of layer " << slot << " already exists" << std::endl ;
			std::terminate() ;
		}

		difID = i.at("center") ;
		temp.shiftY = 32 ;

		it = difList.insert( { difID , temp } ) ;
		if ( !it.second )
		{
			std::cout << "ERROR in geometry : dif " << difID << " of layer " << slot << " already exists" << std::endl ;
			std::terminate() ;
		}

		difID = i.at("right") ;
		temp.shiftY = 64 ;

		it = difList.insert( { difID , temp } ) ;
		if ( !it.second )
		{
			std::cout << "ERROR in geometry : dif " << difID << " of layer " << slot << " already exists" << std::endl ;
			std::terminate() ;
		}
	}

	file.close() ;
}

std::array<int,3> mapping::Mapping::getPadIndex(unsigned int dif_id, unsigned int asic_id, unsigned int chan_id)
{
	std::array<int,3> index = {{0,0,0}} ;

	const auto it = difList.find(dif_id) ;
	int DifZ = it->second.K ;
	int DifY = it->second.shiftY ;

	index[0] = ( 1 + MapILargeHR2[chan_id] + AsicShiftI[asic_id] ) ;
	index[1] = ( 32 - (MapJLargeHR2[chan_id] + AsicShiftJ[asic_id]) ) + DifY ;
	index[2] = DifZ ;

	return index ;
}

void mapping::Mapping::print()
{
	for ( const auto& it : difList )
	{
		std::cout << "Dif " << it.first << "\t: "
				  << "Layer " << it.second.K << "\t"
				  << "Shift " << it.second.shiftY << std::endl ;
	}
}

bool mapping::Mapping::isListedDif( unsigned int difID )
{
	if ( difList.find(difID) != difList.end() )
		return true ;
	return false ;
}
