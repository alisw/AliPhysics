////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Clustering/CenterOfGravityFinder.hpp"
#include "Error.hpp"

namespace dHLT
{
namespace Clustering
{


CenterOfGravityFinder::CenterOfGravityFinder() : ClusterFinder()
{
	Diff_Y = 0.5; // 5000 micron slat size in Y
	Diff_X = 1.0; // 10000 micron slat size in X
	X = 56;
	Y = 128;
	DigitMax = 35; // maximum number of padhits in columns or rows.
	DDLMax = 200; // Maximum number of padhits in one ddl;
	DDLTot = 500; // totoal number of padhits in one ddl;
}


CenterOfGravityFinder::~CenterOfGravityFinder()
{
	// TODO
}


void CenterOfGravityFinder::FindClusters(const ADCStream* stream)
{
	// TODO
}


UInt CenterOfGravityFinder::FillClusterData(
		ClusterPoint* clusters, UInt arraysize
	)
{
	// TODO
	return 0;
}


void CenterOfGravityFinder::Reset() {
	    // TODO
}


} // Clustering
} // dHLT

