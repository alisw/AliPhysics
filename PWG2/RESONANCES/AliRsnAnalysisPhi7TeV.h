//
// Header file for implementation of data analysis aft 900 GeV
//
// Author: A. Pulvirenti
//

#ifndef ALIRSNANALYSISPHI7TEV_H
#define ALIRSNANALYSISPHI7TEV_H

#include "AliAnalysisTaskSE.h"
#include "AliRsnTOFT0maker.h"

class TH1I;
class TTree;

class AliStack;
class AliESDEvent;
class AliESDVertex;
class AliESDpid;
class AliTOFT0maker;
class AliTOFcalib;

class AliRsnAnalysisPhi7TeV : public AliAnalysisTaskSE
{
  public:

    AliRsnAnalysisPhi7TeV(const char *name = "Phi7TeV");
    AliRsnAnalysisPhi7TeV(const AliRsnAnalysisPhi7TeV& copy);
    AliRsnAnalysisPhi7TeV& operator=(const AliRsnAnalysisPhi7TeV& copy);
    virtual ~AliRsnAnalysisPhi7TeV();

    void           SetUseMC(Bool_t yn = kTRUE) {fUseMC = yn;}
    
    void           SetMaxDCAr(Double_t v) {fMaxDCAr = v;}
    void           SetMaxDCAz(Double_t v) {fMaxDCAz = v;}
    void           SetMaxChi2(Double_t v) {fMaxChi2 = v;}
    void           SetMinNTPC(Int_t    n) {fMinNTPC = n;}
    
    void           SetTPCpLimit(Double_t v) {fTPCpLimit = v;}
    void           SetTPCrange(Double_t min, Double_t max) {fMinTPCband = min; fMaxTPCband = max;}
    void           SetTPCpar(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4)
                     {fTPCpar[0]=p0;fTPCpar[1]=p1;fTPCpar[2]=p2;fTPCpar[3]=p3;fTPCpar[4]=p4;}

    void           SetTOFcalibrateESD(Bool_t yn = kTRUE)  {fTOFcalibrateESD = yn;}
    void           SetTOFcorrectTExp (Bool_t yn = kTRUE)  {fTOFcorrectTExp = yn;}
    void           SetTOFuseT0       (Bool_t yn = kTRUE)  {fTOFuseT0 = yn;}
    void           SetTOFtuneMC      (Bool_t yn = kTRUE)  {fTOFtuneMC = yn;}
    void           SetTOFresolution  (Double_t v = 100.0) {fTOFresolution = v;}

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option = "");
    virtual void   Terminate(Option_t *option = "");

  private:

    void     ProcessESD(AliESDEvent *esd, const AliESDVertex *v, AliStack *stack);
    void     ProcessMC(AliStack *stack);

    Bool_t   fUseMC;      // use MC or data?
    
    Short_t  fPDG;        // PDG code
    Float_t  fIM;         // inv mass
    Float_t  fPt;         // transv momentum
    Float_t  fY;          // rapidity
    Float_t  fEta;        // pseudo-rapidity
    
    Double_t fMaxDCAr;    // transverse DCA
    Double_t fMaxDCAz;    // longitudinal DCA
    Double_t fMaxChi2;    // track normalized chi2
    Int_t    fMinNTPC;    // number of TPC clusters

    Double_t fTPCpLimit;  // limit to choose what band to apply
    Double_t fTPCpar[5];  // parameters for TPC bethe-Bloch
    Double_t fMinTPCband; // range for TPC de/dx band - min
    Double_t fMaxTPCband; // range for TPC de/dx band - max

    TTree     *fRsnTreeComp;    // output tree of computed pairs
    TTree     *fRsnTreeTrue;    // output tree of true pairs
    TList     *fOutList;        // list for monitoring histograms
    TH1I      *fHEvents;        // histogram of event types
    TH1I      *fHCuts;          // histogram telling how many tracks don't pass each cut
    
    AliESDpid       *fESDpid;           //! PID manager
    AliTOFT0maker   *fTOFmaker;         //! TOF time0 computator
    AliTOFcalib     *fTOFcalib;         //! TOF calibration
    Bool_t           fTOFcalibrateESD;  //  TOF settings
    Bool_t           fTOFcorrectTExp;   //  TOF settings
    Bool_t           fTOFuseT0;         //  TOF settings
    Bool_t           fTOFtuneMC;        //  TOF settings
    Double_t         fTOFresolution;    //  TOF settings

    // ROOT dictionary
    ClassDef(AliRsnAnalysisPhi7TeV,1)
};

#endif
