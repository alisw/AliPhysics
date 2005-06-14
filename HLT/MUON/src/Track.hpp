////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_TRACK_HPP
#define dHLT_TRACK_HPP

#include "Point.hpp"
#include "RegionOfInterest.hpp"
#include "TriggerRecord.hpp"

namespace dHLT
{


typedef UInt TrackID;


struct Track
{

	TriggerRecordID triggerid;
	ParticleSign sign;
	Float p;   // momentum.
	Float pt;  // transverse momentum.
	Point point[10];  // Computed track coordinates on the 10 tracking chambers.
	ROI region[10];   // Regions of interest from which clusters were used to compute this track.

};


} // dHLT

#endif // dHLT_TRACK_HPP
