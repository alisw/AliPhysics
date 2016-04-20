#ifndef __AliITSUVertexerFast__
#define __AliITSUVertexerFast__

/// \class AliITSUVertexerFast
/// \brief (hopefully) Fast implementation of a primary vertex finder
///
/// This implementation of the primary vertex finder is the first code used to study the performance
/// of a parallel primary vertex finder algorithm (at least for the author).
///
/// The first step is loading the clusters. This is carried out by the LoadClusters() method in 
/// which the standard AliROOT format of the clusters is converted. The information about the space 
/// point is here stored in 4 std::vector for each of the two innermost layer of ITSU.
///
/// The second step in the reconstruction is the so-called trackleting. Using the first to layer of 
/// the ITSU a list of line is compiled. Each of this lines is checked using a veto mask on the 
/// third layer. 
///
/// \author Maximiliano Puccio <maximiliano.puccio@cern.ch>, Università e INFN Torino
/// \date May 7, 2015

#define _TUNE
#define _THCLK

#include <TObject.h>

#include "AliITSUVertexCandidate.h"

#include <omp.h>
#include <math.h>
#include <ctime>
#include <vector>
using std::vector;

struct cluster
{
	float fX;              ///< Position x of the clusters on the layers.
	float fY;              ///< Position y of the clusters on the layers.
    float fZ;              ///< Position z of the clusters on the layers.
	float fP;              ///< Position \f$\phi\f$ of the clusters on the layers.
	float fR;              ///< Position r of the clusters on the layers.
	int   fID;             ///< Id in the LUT

  cluster() : fX(), fY(), fZ(), fR(), fP(), fID(-1) {};
	cluster(float x, float y, float z)
    : fX(x),
    fY(y),
    fZ(z),
    fP(atan2(-fY,-fX)+3.14159265359f),
    fR(sqrt(fX*fX+fY*fY)),
    fID(-1)
	{
	}

	bool operator< (const cluster& clu ) const 
	{
		return (fP < clu.fP);
	}
};

class TClonesArray;
class AliITSUVertexerFast : public TObject {
	public: 
		AliITSUVertexerFast();
		AliITSUVertexerFast(int granularityPhi, int granularityZ, float layerSize);

		void FindVertex(float *xyz);
		void LoadClusters(TClonesArray *clusters[3]);
		void SetRadii(float r0, float r1, float r2) { fRLayer[0] = r0; fRLayer[1] = r1; fRLayer[2] = r2; }
		void SetNumberOfThreads(int n) { omp_set_num_threads(n); }
		int GetMaxNumberOfThreads() { return omp_get_max_threads(); }
		int GetNumberOfThreads() { return omp_get_num_threads(); }

		#ifdef _THCLK
		void SetTimesVecSize() { fTimes.reserve( omp_get_max_threads() ); }
		void StartTime( int idx ) { fTimes[idx] = clock(); }
		void EndTime( int idx ) { fTimes[idx] = clock() - fTimes[idx]; }
		int ElapsedTicks( int idx ) { return fTimes[idx]; }
		float ElapsedTime( int idx ) { return ((float)fTimes[idx])/CLOCKS_PER_SEC; };
 		#endif

	private:

		float IntersectCylinder(Line &l, float r);

		int           fGranularityPhi;    ///< Number of cells in \f$\phi\f$ for the validation mask and for the lookup tables
		int           fGranularityZ;      ///< Number of cells in z for the validation mask
		float         fSizePhi;           ///< Size of the lookup tables cells on \f$\phi\f$ 1/deltaphi
		float         fSizeZ;             ///< Size of the mask cell in z
		float         fHalfLayerLength;   ///< Length of half layer
		float         fRLayer[3];         ///< Radii of the three innermost layers

		#ifdef _THCLK
		vector <clock_t>       fTimes;    ///< Profiling CPU time vector
		#endif

		vector<bool>  fMask;              ///< Validation mask on the third layer  
		                                                          
		vector<cluster> fClusters[3];     ///< Array of clusters vectors. fCluster[j] <==> layer.
		vector<int>   fLUT[3];            ///< Lookup table to browse clusters on layer i

		// vector<Line>  fLines;             ///< Array of tracklets

		#ifdef _TUNE
		vector<int>   fL[3];
		#endif
		ClassDef(AliITSUVertexerFast,1)
};

#endif
