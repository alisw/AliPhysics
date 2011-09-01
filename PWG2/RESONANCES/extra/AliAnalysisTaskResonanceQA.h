#ifndef ALIANALYSISTASKRESONANCEQA_H
#define ALIANALYSISTASKRESONANCEQA_H

// analysis task creating basic QA plots for resonance particles
// Author: Ayben Karasu Uysal
// Computes some QA histograms useful for checking productions and data
// both for counting resonances inside MC, and for checking PID performances

class TH1I;
class TH1F;
class TH2F;

class AliESDpid;
class AliESDEvent;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskResonanceQA : public AliAnalysisTaskSE {
public:

   enum ERsn {
      kPhi = 0,
      kKStar0,
      kRho,
      kLambdaStar,
      kXiStar0,
      kSigmaStarP,
      kSigmaStarM,
      kDeltaPP,
      
      // this must be last and counter
      kResonances
   };
   
   enum EEvtType {
      kBadNoVertex,
      kBadVertexOutZ,
      kGoodEvent,
      
      // this must be last and counter
      kEvtTypes
   };

   AliAnalysisTaskResonanceQA(const char *name = "RsnQA");
   virtual ~AliAnalysisTaskResonanceQA() {}

   virtual void   UserCreateOutputObjects();
   virtual void   UserExec(Option_t *option);
   virtual void   Terminate(Option_t *);
   
   void           SetT0(AliESDpid::EStartTimeType_t ftype) {fT0 = ftype;}
   void           SetPrimaryThr(Double_t d) {fPrimaryThr = d;}
   void           SetVz(Double_t vz)        {fVz = vz;}
   
   const char*    EvtName  (Int_t type) const {return EvtName  ((EEvtType)type);}
   const char*    RsnName  (Int_t type) const {return RsnName  ((ERsn)type);}
   const char*    RsnSymbol(Int_t type) const {return RsnSymbol((ERsn)type);}
   Int_t          RsnPDG   (Int_t type) const {return RsnPDG   ((ERsn)type);}
   const char*    EvtName  (EEvtType type) const;
   const char*    RsnName  (ERsn type)  const;
   const char*    RsnSymbol(ERsn type)  const;
   Int_t          RsnPDG   (ERsn type)  const;

private:

   AliESDpid::EStartTimeType_t fT0;       // T0 type for TOF computation
   Double_t fPrimaryThr;                  // maximum DCA for selecting primary particles w.r. to gen primary vertex
   Double_t fVz;                          // maximum VZ for primary vertex (reconstructed)
                                         
   TList *fOutputList;                    // list with output histograms
   
   TH1I  *fSelectedEvts;                  // selected events
   TH1F  *fEventVz;                       // Vz distribution of events
   TH2F  *fdEdxTPC;                       // TPC PID QA
   TH2F  *fdEdxITS;                       // ITS PID QA
   TH2F  *fTOFpid;                        // TOF PID QA
   TH2F  *fDCAXYvsPtBeforeCuts;           // DCA QA r
   TH2F  *fDCAZvsPtBeforeCuts;            // DCA QA z
   TH2F  *fNClusterPtBeforeCuts;          // N cluster TPC
   TH2F  *fNFindableClusterPtBeforeCuts;  // N findable TPC clusters
   TH2F  *fNCrossedRowsTPCPtBeforeCuts;   // crossed rows
   TH3F  *fRsnYPt[2][kResonances];        // rapidity vs pt vs event type distribution of resonances
   TH1I  *fProducedParticles;             // synoptic of all resonances
      
   AliESDEvent      *fESD;                //! temporary object (event)
   AliESDpid        *fESDpid;             //  temporary object (PID)
   AliESDtrackCuts  *fTrackCuts;          //  temporary object (quality track cuts)
   
   AliAnalysisTaskResonanceQA(const AliAnalysisTaskResonanceQA&);                // disabled
   AliAnalysisTaskResonanceQA& operator=(const AliAnalysisTaskResonanceQA&);     // disabled

   ClassDef(AliAnalysisTaskResonanceQA, 1);   // Resonance QA class
};

#endif


