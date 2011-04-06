#ifndef ALIANALYSISTASTRESONANCEQA_H
#define ALIANALYSISTASTRESONANCEQA_H

// analysis task creating basic QA plots for resonance particles
// Author: Ayben Karasu Uysal

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
      kXiStarM,
      kSigmaStarP,
      kSigmaStar0,
      kSigmaStarM,
      kDeltaPP,
      
      // this must be last and counter
      kResonances
   };

   AliAnalysisTaskResonanceQA(const char *name = "RsnQA");
   virtual ~AliAnalysisTaskResonanceQA() {}

   virtual void   UserCreateOutputObjects();
   virtual void   UserExec(Option_t *option);
   virtual void   Terminate(Option_t *);
   
   void           SetT0(AliESDpid::EStartTimeType_t ftype) {fT0 = ftype;}
   void           SetPrimaryThr(Double_t d) {fPrimaryThr = d;}
   
   const char*    RsnName  (Int_t type) {return RsnName  ((ERsn)type);}
   const char*    RsnSymbol(Int_t type) {return RsnSymbol((ERsn)type);}
   Int_t          RsnPDG   (Int_t type) {return RsnPDG   ((ERsn)type);}
   const char*    RsnName  (ERsn type);
   const char*    RsnSymbol(ERsn type);
   Int_t          RsnPDG   (ERsn type);

private:

   AliESDpid::EStartTimeType_t fT0;
   Double_t fPrimaryThr;
   
   TList *fOutputList;
   
   TH1I  *fSelectedEvts;
   TH2F  *fdEdxTPC;
   TH2F  *fdEdxITS;
   TH2F  *fTOFpid;
   TH2F  *fDCAXYvsPtBeforeCuts;
   TH2F  *fDCAZvsPtBeforeCuts;
   TH2F  *fNClusterPtBeforeCuts;
   TH2F  *fNFindableClusterPtBeforeCuts;
   TH2F  *fNCrossedRowsTPCPtBeforeCuts;
   TH2F  *fRsnYPt[2][kResonances];
   TH1I  *fProducedParticles;
   
   AliESDEvent      *fESD;
   AliESDpid        *fESDpid;
   AliESDtrackCuts  *fTrackCuts;
   
   AliAnalysisTaskResonanceQA(const AliAnalysisTaskResonanceQA&);
   AliAnalysisTaskResonanceQA& operator=(const AliAnalysisTaskResonanceQA&);

   ClassDef(AliAnalysisTaskResonanceQA, 1);
};

#endif


