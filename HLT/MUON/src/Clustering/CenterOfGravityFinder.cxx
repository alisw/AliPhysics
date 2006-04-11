////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Clustering/CenterOfGravityFinder.hpp"
#include "Error.hpp"


AliHLTMUONCoreCenterOfGravityFinder::AliHLTMUONCoreCenterOfGravityFinder()
	: AliHLTMUONCoreClusterFinder()
{
	fDiff_Y = 0.5; // 5000 micron slat size in Y
	fDiff_X = 1.0; // 10000 micron slat size in X
	fX = 56;
	fY = 128;
	fDigitMax = 35; // maximum number of padhits in columns or rows.
	fDDLMax = 200; // Maximum number of padhits in one ddl;
	fDDLTot = 500; // totoal number of padhits in one ddl;
}


AliHLTMUONCoreCenterOfGravityFinder::~AliHLTMUONCoreCenterOfGravityFinder()
{
	// TODO
}


void AliHLTMUONCoreCenterOfGravityFinder::FindClusters(
		const AliHLTMUONCoreADCStream* /*stream*/
	)
{
	// TODO
}


UInt AliHLTMUONCoreCenterOfGravityFinder::FillClusterData(
		AliHLTMUONCoreClusterPoint* /*clusters*/, UInt /*arraysize*/
	)
{
	// TODO
	return 0;
}


void AliHLTMUONCoreCenterOfGravityFinder::Reset()
{
	// TODO
}

