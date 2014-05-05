#ifndef ALIITSUTRACKERCooked_H
#define ALIITSUTRACKERCooked_H

//-------------------------------------------------------------------------
//                   The stand-alone ITSU tracker
//    The pattern recongintion based on the "cooked covariance" approach
//-------------------------------------------------------------------------

#include "AliTracker.h"

class TTree;
class TClonesArray;
class TObjArray;

class AliESDEvent;
class AliCluster;
class AliITSUClusterPix;

//-------------------------------------------------------------------------
class AliITSUTrackerCooked : public AliTracker {
public:
  enum{kNLayers=7, kMaxClusterPerLayer=70000};
  AliITSUTrackerCooked();
  virtual ~AliITSUTrackerCooked();

  // These functions must be implemented 
  Int_t Clusters2Tracks(AliESDEvent *event);
  Int_t PropagateBack(AliESDEvent *event);
  Int_t RefitInward(AliESDEvent *event);
  Int_t LoadClusters(TTree *ct);
  void UnloadClusters();
  AliCluster *GetCluster(Int_t index) const;

  // internal helper classes
  class AliITSUlayer;

protected:
  AliITSUTrackerCooked(const AliITSUTrackerCooked&);
  // Other protected functions
  Int_t MakeSeeds(Int_t layer1, Int_t layer2);
  Bool_t AddCookedSeed(const Float_t r1[3], Int_t l1, Int_t i1,
                       const Float_t r2[3], Int_t l2, Int_t i2,
                       const AliCluster *c3,Int_t l3, Int_t i3);

private:
  AliITSUTrackerCooked &operator=(const AliITSUTrackerCooked &tr);

  // Data members
  // Internal tracker arrays, layers, modules, etc
  static AliITSUlayer fgLayers[kNLayers];// Layers
    
  TObjArray *fSeeds; // Track seeds
  
  ClassDef(AliITSUTrackerCooked,1)   //ITSU stand-alone tracker
};



// The helper classes
class AliITSUTrackerCooked::AliITSUlayer {
  public:
    AliITSUlayer();
   ~AliITSUlayer(){DeleteClusters();}

    void InsertClusters(TClonesArray *clusters);
    void SetR(Double_t r) {fR=r;}
    void DeleteClusters();

    Int_t FindClusterIndex(Double_t z) const;
    Double_t GetR() const {return fR;}
    AliCluster *GetCluster(Int_t i) const { return fClusters[i]; } 
    Int_t GetNumberOfClusters() const {return fN;}

  protected:
    AliITSUlayer(const AliITSUlayer&);
    AliITSUlayer &operator=(const AliITSUlayer &tr);  
    Int_t InsertCluster(AliCluster *c);

    Double_t fR;                // mean radius of this layer

    AliCluster *fClusters[kMaxClusterPerLayer]; // array of clusters
    Int_t fN; //number of clusters 

    Int_t fIndex[kMaxClusterPerLayer]; // Indexes of selected clusters 
    Int_t fNsel;                       // number of selected clusters
  };


#endif
