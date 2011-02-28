//
// Header file for implementation of data analysis aft 900 GeV
//
// Author: A. Pulvirenti
//

#ifndef ALIRSNANALYSISPHI900GEV_H
#define ALIRSNANALYSISPHI900GEV_H

#include "AliAnalysisTaskSE.h"
#include "AliRsnTOFT0maker.h"

class TH1I;
class TTree;

class AliStack;
class AliESDEvent;
class AliESDVertex;

class AliRsnAnalysisPhi900GeV : public AliAnalysisTaskSE {
public:

   AliRsnAnalysisPhi900GeV(const char *name = "Phi900GeV");
   AliRsnAnalysisPhi900GeV(const AliRsnAnalysisPhi900GeV& copy);
   AliRsnAnalysisPhi900GeV& operator=(const AliRsnAnalysisPhi900GeV& copy);
   virtual ~AliRsnAnalysisPhi900GeV();

   void            SetUseMC(Bool_t yn = kTRUE) {fUseMC = yn;}
   void            SetMaxDCAr(Double_t v) {fDCAr = v;}
   void            SetMaxDCAz(Double_t v) {fDCAz = v;}
   void            SetMaxChi2(Double_t v) {fChi2 = v;}
   void            SetMinNTPC(Int_t    n) {fNTPC = n;}
   void            SetTPCparams(Bool_t isMC);
   void            SetTPCrange(Double_t min, Double_t max) {fMinTPC = min; fMaxTPC = max;}
   void            SetTPCpar(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4)
   {fTPCpar[0] = p0; fTPCpar[1] = p1; fTPCpar[2] = p2; fTPCpar[3] = p3; fTPCpar[4] = p4;}

   void           SetTOFESD(Bool_t yn = kTRUE) {fTOFESD = yn;}
   void           SetTOFSigma(Double_t v) {fTOFSigma = v;}
   void           SetTOFSettings(AliRsnTOFT0maker::ESettings set) {fTOFSettings = set;}

   virtual void    UserCreateOutputObjects();
   virtual void    UserExec(Option_t *option = "");
   virtual void    Terminate(Option_t *option = "");

private:

   void     ProcessESD(AliESDEvent *esd, const AliESDVertex *v, Double_t time0, AliStack *stack);
   void     ProcessMC(AliStack *stack);
   Double_t AlephBB(Double_t p, Double_t mass = 0.493677) const;
   Double_t RemakeTOFtimeMC(AliESDEvent *& esd);

   Bool_t   fUseMC;  // use MC or data?

   Short_t  fPDG;    // PDG code
   Float_t  fIM;     // inv mass
   Float_t  fPt;     // transv momentum
   Float_t  fY;      // rapidity
   Float_t  fEta;    // pseud-o-rapidity

   Double_t fDCAr;   // transverse DCA
   Double_t fDCAz;   // longitudinal DCA
   Double_t fChi2;   // track chi2
   Int_t    fNTPC;   // number of TPC clusters

   Double_t fTPCpar[5];  // parameters for TPC bethe-Bloch
   Double_t fMinTPC;     // range for TPC de/dx band - min
   Double_t fMaxTPC;     // range for TPC de/dx band - max

   TTree   *fOutTree[2]; // output tree
   TList   *fOutList;    // list for output event counts
   TH1I    *fHEvents;    // histogram of event types

   Bool_t                       fTOFESD;              //  TOF flag to check if ESD data should be used
   Double_t                     fTOFSigma;            //  TOF default resolution
   AliRsnTOFT0maker            *fTOFmaker;            //! TOF time0 computator
   AliRsnTOFT0maker::ESettings  fTOFSettings;         //  TOF settings

   // ROOT dictionary
   ClassDef(AliRsnAnalysisPhi900GeV, 1)
};

#endif
