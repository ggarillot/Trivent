#ifndef Mapping_H_
#define Mapping_H_

#include <iostream>
#include <vector>
#include <map>
#include <array>

#include <string>
#include <fstream>

#include "json.hpp"


namespace mapping
{

unsigned int getDifID(int cell_id) ;
unsigned int getAsicID(int cell_id) ;
unsigned int getPadID(int cell_id) ;

struct Dif
{
		int K ;
		int shiftY ;
} ;

//=====================================
const unsigned short MapJLargeHR2[64] =
{1,1,2,2,3,3,4,4,5,5,6,6,7,7,4,3,
 2,0,0,1,0,5,6,7,7,6,5,4,3,2,1,0,
 0,1,2,3,4,5,6,7,7,6,0,0,1,1,0,2,
 3,4,5,7,7,6,6,5,5,4,4,3,3,2,2,1};
const unsigned short MapILargeHR2[64] =
{1,0,1,0,1,0,1,0,0,1,0,1,0,1,2,2,
 2,0,1,2,2,2,2,2,3,3,3,3,3,3,3,3,
 4,4,4,4,4,4,4,4,5,5,5,6,5,6,7,5,
 5,5,5,6,7,6,7,6,7,6,7,7,6,7,6,7};

//      (J Axis)
// 7	|46|63|61|59|58|56|54|52|
// 6	|43|45|62|60|57|55|53|51|
// 5	|42|44|47|48|49|50|41|40|
// 4	|32|33|34|35|36|37|38|39|
// 3	|31|30|29|28|27|26|25|24|  	     TOP VIEW (ASICs SIDE)
// 2	|20|19|16|15|14|21|22|23|
// 1	|18|00|02|04|06|09|11|13|
// 0	|17|01|03|05|07|08|10|12|
//       0  1  2  3  4  5  6  7    (I Axis)  ----->
//				|	 |
//				|DIFF|
//				|____|

const unsigned short AsicShiftI[49]=
{	0,	0,	0,	0,	0,	8,	8,
	8,	8,	16,	16,	16,	16,	24,
	24,	24,	24,	32,	32,	32,	32,
	40,	40,	40,	40,	48,	48,	48,
	48,	56,	56,	56,	56,	64,	64,
	64,	64,	72,	72,	72,	72,	80,
	80,	80,	80,	88,	88,	88,	88};
const unsigned short AsicShiftJ[49]=
{	0,	0,	8,	16,	24,	24,	16,
	8,	0,	0,	8,	16,	24,	24,
	16,	8,	0,	0,	8,	16,	24,
	24,	16,	8,	0,	0,	8,	16,
	24,	24,	16,	8,	0,	0,	8,
	16,	24,	24,	16,	8,	0,	0,
	8,	16,	24,	24,	16,	8,	0};

struct Mapping
{
		Mapping() = default ;
		Mapping(std::string jsonFile) ;

		~Mapping() = default ;

		std::array<int,3> getPadIndex(unsigned int dif_id, unsigned int asic_id, unsigned int chan_id) ;

		void print() ;

		bool isListedDif( unsigned int difID ) ;

	protected :
		std::map<unsigned int, Dif> difList = {} ;
} ;

} //namespace mapping
#endif
