#ifndef AliTPCCALIBTRACKSCUTS_H
#define AliTPCCALIBTRACKSCUTS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////
//                                                  //
//     Class to specify cuts for track analysis     //
//     with AliTPCcalibTracks                       //
//                                                  //
//////////////////////////////////////////////////////



#include <TNamed.h>
#include <TObjString.h>

class TChain;
class AliTPCseed;
class AliVTrack;

using namespace std;

class AliTPCcalibTracksCuts: public TNamed {

public:
   AliTPCcalibTracksCuts(Int_t minClusters, Float_t minRatio, Float_t max1pt,
      Float_t edgeXZCutNoise, Float_t edgeThetaCutNoise);
   AliTPCcalibTracksCuts(AliTPCcalibTracksCuts *cuts);
   AliTPCcalibTracksCuts();
   virtual ~AliTPCcalibTracksCuts();
   static  AliTPCcalibTracksCuts  *CreateCuts(char* ctype);

   Int_t AcceptTrack(const AliTPCseed * track) const;
   Int_t AcceptTrack(const AliVTrack * track) const;

   void SetMinClusters(Int_t minClusters){fMinClusters = minClusters;}
   void SetMinRatio(Float_t minRatio){fMinRatio = minRatio;}
   void SetMax1pt(Float_t max1pt){fMax1pt = max1pt;}
   void SetEdgeXYCutNoise(Float_t edgeCutNoise){fEdgeYXCutNoise = edgeCutNoise;}
   void SetEdgeThetaCutNoise(Float_t edgeCutNoise){fEdgeThetaCutNoise = edgeCutNoise;}
   Int_t   GetMinClusters() const {return fMinClusters;}
   Float_t GetMinRatio() const {return fMinRatio;}
   Float_t GetMax1pt() const {return fMax1pt;}
   Float_t GetEdgeYXCutNoise() const {return fEdgeYXCutNoise;}
   Float_t GetEdgeThetaCutNoise() const {return fEdgeThetaCutNoise;}
   virtual void Print(Option_t* option = "") const;
   
private:
   Int_t   fMinClusters;         // number of clusters
   Float_t fMinRatio;            // kMinRratio = 0.4
   Float_t fMax1pt;              // kMax1pt = 0.5
   Float_t fEdgeYXCutNoise;      // kEdgeYXCutNoise = 0.13
   Float_t fEdgeThetaCutNoise;   // kEdgeThetaCutNoise = 0.018

protected:         
   ClassDef(AliTPCcalibTracksCuts,1)
};


#endif
