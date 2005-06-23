////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_CLUSTERING_CENTER_OF_GRAVITY_FINDER_HPP
#define dHLT_CLUSTERING_CENTER_OF_GRAVITY_FINDER_HPP

#include "Clustering/ClusterFinder.hpp"

namespace dHLT
{
namespace Clustering
{
      
      
class CenterOfGravityFinder : public ClusterFinder
{
public:
	
	CenterOfGravityFinder();
	
	virtual ~ CenterOfGravityFinder();
	
	// Inherited from ClusterFinder
	virtual void FindClusters(const ADCStream * stream);
	virtual UInt FillClusterData(ClusterPoint * clusters,
				     UInt arraysize);
	virtual void Reset();
private:
	FILE* lut13;
	FILE* lut14;
	FILE* lut15;
	FILE* lut16;
	FILE* lut17;
	FILE* lut18;
	FILE* lut19;
	FILE* lut20;
	
	Float_t Diff_Y; // 5000 micron slat size in Y
	Float_t Diff_X; // 10000 micron slat size in X
	UInt_t X; // number of columns in bending plane
	UInt_t Y; // number of rows in non-bending plane
	UInt_t DigitMax; // maximum number of padhits in columns or rows.
	UInt_t DDLMax; // Maximum number of padhits in one ddl;
	UInt_t DDLTot; // totoal number of padhits in one ddl;

};
      

} // Clustering
} // dHLT

#endif // dHLT_CLUSTERING_CENTER_OF_GRAVITY_FINDER_HPP

