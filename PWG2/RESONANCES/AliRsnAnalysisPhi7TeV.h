//
// Header file for implementation of data analysis aft 900 GeV
//
// Author: A. Pulvirenti
//

#ifndef ALIRSNANALYSISPHI7TEV_H
#define ALIRSNANALYSISPHI7TEV_H

#include "AliAnalysisTaskSE.h"
#include "AliRsnTOFT0maker.h"
#include "AliESDtrackCuts.h"

class TH1I;
class TH1F;
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
  
    enum
    {
      kGoodTracksPrimaryVertex = 0,
      kGoodSPDPrimaryVertex    = 1,
      kFarTracksPrimaryVertex  = 2,
      kFarSPDPrimaryVertex     = 3,
      kNoGoodPrimaryVertex     = 4,
      kEvaluations             = 5
    };

    AliRsnAnalysisPhi7TeV(const char *name = "Phi7TeV");
    AliRsnAnalysisPhi7TeV(const AliRsnAnalysisPhi7TeV& copy);
    AliRsnAnalysisPhi7TeV& operator=(const AliRsnAnalysisPhi7TeV& copy);
    virtual ~AliRsnAnalysisPhi7TeV();

    void             SetUseMC(Bool_t yn = kTRUE) {fUseMC = yn;}
    void             SetCheckITS(Bool_t yn = kTRUE) {fCheckITS = yn;}
    void             SetCheckTPC(Bool_t yn = kTRUE) {fCheckTPC = yn;}
    void             SetCheckTOF(Bool_t yn = kTRUE) {fCheckTOF = yn;}
    
    void             SetMaxVz(Double_t v)   {fMaxVz = v;}
    
    void             SetITSband(Double_t v) {fMaxITSband = v;}
    
    void             SetTPCpLimit(Double_t v) {fTPCpLimit = v;}
    void             SetTPCrange(Double_t min, Double_t max) {fMinTPCband = min; fMaxTPCband = max;}
    void             SetTPCpar(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4)
                       {fTPCpar[0]=p0;fTPCpar[1]=p1;fTPCpar[2]=p2;fTPCpar[3]=p3;fTPCpar[4]=p4;}

    void             SetTOFcalibrateESD(Bool_t yn = kTRUE)  {fTOFcalibrateESD = yn;}
    void             SetTOFcorrectTExp (Bool_t yn = kTRUE)  {fTOFcorrectTExp = yn;}
    void             SetTOFuseT0       (Bool_t yn = kTRUE)  {fTOFuseT0 = yn;}
    void             SetTOFtuneMC      (Bool_t yn = kTRUE)  {fTOFtuneMC = yn;}
    void             SetTOFresolution  (Double_t v = 100.0) {fTOFresolution = v;}

    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option = "");
    virtual void     Terminate(Option_t *option = "");
    
    Int_t            EventEval(AliESDEvent *esd);
    AliESDtrackCuts* GetCutsTPC() {return &fESDtrackCutsTPC;}
    AliESDtrackCuts* GetCutsITS() {return &fESDtrackCutsITS;}

  private:

    void     ProcessESD(AliESDEvent *esd, AliStack *stack);
    void     ProcessMC(AliStack *stack);
    void     AddEntryFromESD(AliESDEvent *esd, Int_t i1, Int_t i2, Int_t its1, Int_t its2, Short_t charge, AliStack *stack = 0x0);

    Bool_t   fUseMC;      // use MC or data?
    Bool_t   fCheckITS;   // chec ITS PID?
    Bool_t   fCheckTPC;   // chec TPC PID?
    Bool_t   fCheckTOF;   // chec TOF PID?
    
    Short_t  fPDG;        // PDG code
    Short_t  fCh;         // control flag for like/unlike sign
    Short_t  fITS[2];     // check flag to know if one or both candidates are ITS standalone
    Float_t  fIM;         // inv mass
    Float_t  fPt;         // transv momentum
    Float_t  fY;          // rapidity
    Float_t  fEta;        // pseudo-rapidity
    
    Double_t fMaxVz;      // range in Z of primary vertex w.r. to origin
    
    Double_t fMaxITSband; // range for ITS de/dx band

    Double_t fTPCpLimit;  // limit to choose what band to apply
    Double_t fTPCpar[5];  // parameters for TPC bethe-Bloch
    Double_t fMinTPCband; // range for TPC de/dx band - min
    Double_t fMaxTPCband; // range for TPC de/dx band - max

    TTree     *fRsnTreeComp;    // output tree of computed pairs
    TTree     *fRsnTreeTrue;    // output tree of true pairs
    TList     *fOutList;        // list for monitoring histograms
    TH1I      *fHEvents;        // histogram of event types
    TH1F      *fVertexX[2];     // histogram of X coordinate of primary vertex ([0] = tracks, [1] = SPD)
    TH1F      *fVertexY[2];     // histogram of Y coordinate of primary vertex ([0] = tracks, [1] = SPD)
    TH1F      *fVertexZ[2];     // histogram of Z coordinate of primary vertex ([0] = tracks, [1] = SPD)
    
    AliESDtrackCuts  fESDtrackCutsTPC;  //  ESD standard defined track cuts for TPC tracks
    AliESDtrackCuts  fESDtrackCutsITS;  //  ESD standard defined track cuts for ITS-SA tracks
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
