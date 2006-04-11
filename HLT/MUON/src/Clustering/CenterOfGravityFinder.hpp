////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONCORECENTEROFGRAVITYFINDER_H
#define ALIHLTMUONCORECENTEROFGRAVITYFINDER_H

#include "Clustering/ClusterFinder.hpp"


class AliHLTMUONCoreCenterOfGravityFinder : public AliHLTMUONCoreClusterFinder
{
public:
	
	AliHLTMUONCoreCenterOfGravityFinder();
	
	virtual ~AliHLTMUONCoreCenterOfGravityFinder();
	
	// Inherited from ClusterFinder
	virtual void FindClusters(const AliHLTMUONCoreADCStream* stream);
	virtual UInt FillClusterData(AliHLTMUONCoreClusterPoint* clusters, UInt arraysize);
	virtual void Reset();

private:

	FILE* fLut13;
	FILE* fLut14;
	FILE* fLut15;
	FILE* fLut16;
	FILE* fLut17;
	FILE* fLut18;
	FILE* fLut19;
	FILE* fLut20;
	
	Float_t fDiff_Y; // 5000 micron slat size in Y
	Float_t fDiff_X; // 10000 micron slat size in X
	UInt_t fX; // number of columns in bending plane
	UInt_t fY; // number of rows in non-bending plane
	UInt_t fDigitMax; // maximum number of padhits in columns or rows.
	UInt_t fDDLMax; // Maximum number of padhits in one ddl;
	UInt_t fDDLTot; // totoal number of padhits in one ddl;

};


#endif // ALIHLTMUONCORECENTEROFGRAVITYFINDER_H

