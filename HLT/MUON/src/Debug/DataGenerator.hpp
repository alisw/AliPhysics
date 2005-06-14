////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DEBUG_DATA_GENERATOR_HPP
#define dHLT_DEBUG_DATA_GENERATOR_HPP

#include <math.h>
#include <stdlib.h>
#include "Debug/TriggerSource.hpp"
#include "Debug/ClusterSource.hpp"

namespace dHLT
{
namespace Debug
{


class DataGenerator
{
public:

	DataGenerator();

	void GenerateData(
			const EventID eventid, const UInt number_of_tracks,
			TriggerSource& triggers, ClusterSource& clusters
		);

protected:

	Float Random();
	Float GausRandom();

	void ComputeTrack(Float theta, Float phi, Float momentum, Float charge, ClusterPoint track[14]);


	Float BL;
	Float chamber_z[14];
	Float zf;
};


} // Debug
} // dHLT

#endif // dHLT_DEBUG_DATA_GENERATOR_HPP
