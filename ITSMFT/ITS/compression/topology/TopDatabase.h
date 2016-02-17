#ifndef TOPDATABASE_H
#define TOPDATABASE_H
#include "TObject.h"
#include "TBits.h"
#include "./Topology.h"
#include <Riostream.h>
#include "TArrayI.h"
#include "TArrayF.h"
#include "TObjArray.h"
#include "AliITSMFTClusterPix.h"


class TopDatabase : public TObject {

 public:

  enum SortMode_t{kFrequency=0, kHashes=1};//fMode
  enum FitIndex_t{kDeltaZmean=0, kDeltaZmeanErr=1, kDeltaXmean=2, kDeltaXmeanErr=3, kDeltaZsigma=4, kDeltaZsigmaErr=5,
		  kDeltaXsigma=6, kDeltaXsigmaErr=7, kChi2z=8, kChi2x=9, kNDFx=10, kNDFz=11, kFitLength=12}; //position in Topology fArrFit

  TopDatabase();
  TopDatabase(TopDatabase &ogg);
  ~TopDatabase();
  void AccountTopology(const AliITSMFTClusterPix &cluster, Float_t dX, Float_t dZ, Float_t alpha, Float_t beta);

  Int_t GetN() const {return fN;}
  Int_t GetTotClusters() const {return fTotClusters;}
  Float_t GetThreshold() const {return fThreshold;}
  Int_t GetOverThr() const {return fOverThr;}
  Int_t GetNGroups() const {return fNGroups;}
  Int_t GetNmax() const {return fNmax;}

  void EndAndSort(Int_t mode = kHashes);//to end the database and sort key wrt hashes, in ascending order
  void PrintDB(const char* output = "Database.txt") const; //print the database on a txt file
  void SetThresholdCumulative(Float_t cumulative);
  //Threshold is the frequency for which you have a fraction = cumulative of topology not in groups
  void SetThreshold(Float_t thr);
  void Grouping(Int_t NumberofShiftXbins, Int_t NumberofShiftZbins);//return patterns over threshold
  void SetNmax(Int_t a) { fNmax = a;}
  Int_t FromCluster2GroupID(const AliITSMFTClusterPix &cl) const;
  Bool_t TestChain2Ways(const AliITSMFTClusterPix &cl) const;


 private:
  Int_t fN; //length of arrays
  TObjArray fArrTopologies;//array of topologies (class Topology)
  Int_t fTotClusters;
  Float_t fThreshold;//frequency threshold
  Int_t fOverThr;//number of patterns topologies over threshold
  Int_t fNGroups;
  Int_t fNmax;//patterns above this number (included) belong to a "junk" bin
  TObjArray fArrHisto;
  TObjArray* GetArrTopologies() {return &fArrTopologies;}
  void ExpandDB(const TBits* patt);


ClassDef(TopDatabase,1)

};

#endif
