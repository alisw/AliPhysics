#ifndef ALIANALYSISTASKEVENTCOUNT_H
#define ALIANALYSISTASKEVENTCOUNT_H

// analysis task creating basic QA plots for resonance particles
// Author: Ayben Karasu Uysal

class TH1I;
class TH1F;
class TH2F;

class AliESDpid;
class AliESDEvent;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEventCount : public AliAnalysisTaskSE {
public:

   enum EType {
      kAll = 0,
      kPassPhysicsSel,
      kGoodPrimaryVertex,
      kHas1TrackAny,
      kHas1TrackQuality
      
      kTypes
   };

   AliAnalysisTaskEventCount(const char *name = "RsnQA");
   virtual ~AliAnalysisTaskEventCount() {}

   virtual void   UserCreateOutputObjects();
   virtual void   UserExec(Option_t *option);
   virtual void   Terminate(Option_t *);
   
   void               SetVzRange(Double_t vz)                      {fVz = vz;}
   void               SetAcceptTPCPrimaryVertex(Bool_t yn = kTRUE) {fAcceptTPC = yn;}
   static const char* Title(EType type);

private:

   Bool_t            fAcceptTPC;     // accept TPC primary vertex?
   Double_t          fVz;            // symmetric range around 0
   TList            *fOutputList;    // output container
   Bool_t            fCheck[kTypes]; // check if accepted according to criterion in enum
   TH1I             *fHEvent;        // histogram with event counts
   AliESDtrackCuts  *fTrackCuts;     // cut checker
   
   AliAnalysisTaskEventCount(const AliAnalysisTaskEventCount&);
   AliAnalysisTaskEventCount& operator=(const AliAnalysisTaskEventCount&);

   ClassDef(AliAnalysisTaskEventCount, 1);
};

inline const char* Title(AliAnalysisTaskEventCount::EType type)
{
//
// Title string depending on type
//
   
   switch (type) {
      case kAll              : return "All";
      case kPassPhysicsSel   : return "Phys Sel OK";
      case kGoodPrimaryVertex: return "Good primary vertex";
      case kHas1TrackAny     : return "1 track (any)";
      case kHas1TrackQuality : return "1 track (quality OK)";
      default                : return "Unknown";
   }
}

#endif


