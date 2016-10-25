#ifndef ALICALOTRACKMATCHER_H
#define ALICALOTRACKMATCHER_H

#include "AliAnalysisTaskSE.h"
#include "AliEMCALGeometry.h"
#include "AliPHOSGeometry.h"
#include <vector>
#include <map>
#include <utility>

using namespace std;

class AliCaloTrackMatcher : public AliAnalysisTaskSE {

 public:

   AliCaloTrackMatcher(const char *name, Int_t clusterType);
   //Uncopyable & operator=(const Uncopyable&);

   virtual ~AliCaloTrackMatcher();                            //virtual destructor
   void UserCreateOutputObjects();

   virtual void UserExec(Option_t *option);
   virtual void Terminate(Option_t *);

   void Initialize(Int_t runNumber);
   void ProcessEvent(AliVEvent *event);

   void SetV0ReaderName(TString name)    {fV0ReaderName = name; return;}
   void SetMatchingResidual(Float_t res) {fMatchingResidual = res; return;}
   void SetMatchingWindow(Float_t win) {fMatchingWindow = win; return;}

   Bool_t GetTrackClusterMatchingResidual(Int_t trackID, Int_t clusterID, Float_t &dEta, Float_t &dPhi);
   Int_t GetNMatchedTrackIDsForCluster(Int_t clusterID, Float_t dEtaPos, Float_t dEtaNeg, Float_t dPhiPos, Float_t dPhiNeg);
   Int_t GetNMatchedClusterIDsForTrack(Int_t trackID, Float_t dEtaPos, Float_t dEtaNeg, Float_t dPhiPos, Float_t dPhiNeg);
   vector<Int_t> GetMatchedTrackIDsForCluster(Int_t clusterID, Float_t dEtaPos, Float_t dEtaNeg, Float_t dPhiPos, Float_t dPhiNeg);
   vector<Int_t> GetMatchedClusterIDsForTrack(Int_t trackID, Float_t dEtaPos, Float_t dEtaNeg, Float_t dPhiPos, Float_t dPhiNeg);

 private:
   typedef pair<Int_t, Int_t> pairInt;
   typedef pair<Float_t, Float_t> pairFloat;
   typedef map<pairInt, Int_t> mapT;

   AliCaloTrackMatcher (const AliCaloTrackMatcher&); // not implemented
   AliCaloTrackMatcher & operator=(const AliCaloTrackMatcher&); // not implemented

   Int_t                 fClusterType;            // EMCal(1), PHOS(2) or not running (0)
   TString               fV0ReaderName;           // Name of V0Reader
   Double_t              fMatchingWindow;         // matching window to prevent unnecessary propagations
   Float_t               fMatchingResidual;       // matching residual below which track <-> cluster associations should be stored
   Int_t                 fRunNumber;              // current run number

   AliEMCALGeometry*     fGeomEMCAL;              // pointer to EMCAL geometry
   AliPHOSGeometry*      fGeomPHOS;               // pointer to PHOS geometry

   multimap<Int_t,Int_t> fMapTrackToCluster;      // connects a given track ID with all associated cluster IDs
   multimap<Int_t,Int_t> fMapClusterToTrack;      // connects a given cluster ID with all associated track IDs

   Int_t                 fNEntries;               // number of current TrackID/ClusterID -> Eta/Phi connections
   vector<pairFloat>     fVectorDeltaEtaDeltaPhi; // vector of all matching residuals for a specific
   mapT                  fMap_TrID_ClID_ToIndex;  // map tuple of (trackID,clusterID) to index in vector vecDeltaEtaDeltaPhi

   ClassDef(AliCaloTrackMatcher,1)
};

#endif
