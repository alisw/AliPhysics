#ifndef ALIITSUTRACKERSA_H
#define ALIITSUTRACKERSA_H

//-------------------------------------------------------------------------
//                   The stand-alone ITSU tracker
//     It reads AliITSUClusterPix clusters and writes the tracks to the ESD
//-------------------------------------------------------------------------

#include "AliTracker.h"
#include "AliITSUGeomTGeo.h"
#include <TClonesArray.h>
#include <vector>
#include "AliITSUTrackerSAaux.h"   // Structs and other stuff 

class TTree;
class AliCluster;
class AliESDEvent;

using std::vector;

//-------------------------------------------------------------------------
class AliITSUTrackerSA : public AliTracker {
public:
  AliITSUTrackerSA();
  virtual ~AliITSUTrackerSA() {} ;

  // These functions must be implemented 
  Int_t Clusters2Tracks(AliESDEvent *event);
  Int_t PropagateBack(AliESDEvent *event);
  Int_t RefitInward(AliESDEvent *event);
  Int_t LoadClusters(TTree *ct);
  void UnloadClusters();
  AliCluster *GetCluster(Int_t index) const;

  // Possibly, other public functions


protected:
  AliITSUTrackerSA(const AliITSUTrackerSA&);

  void MakeDoublets();
  void MakeTriplets();
  void CASelection();
  void GlobalFit();
  void ChiSquareSelection();
  // Other protected functions
  // (Sorting, labeling, calculations of "roads", etc)

private:
  AliITSUTrackerSA &operator=(const AliITSUTrackerSA &tr);

  // Data members
  // Internal tracker arrays, layers, modules, etc
  vector<itsCluster> fClusters[7];
  TClonesArray fClustersTC[7];
  vector<nPlets> fDoublets[6];
  Int_t *fIndex[7];
  Int_t fNClusters[7];
  Int_t fNDoublets[6];
  Float_t fPhiCut;
  Float_t fRPhiCut;
  Float_t fZCut;

  ClassDef(AliITSUTrackerSA,1)   //ITSU stand-alone tracker
};

#endif
