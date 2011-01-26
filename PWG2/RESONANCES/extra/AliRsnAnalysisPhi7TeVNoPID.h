//
// Header file for implementation of data analysis aft 900 GeV
//
// Author: A. Pulvirenti
//

#ifndef ALIRSNANALYSISPHI7TEVNOPID_H
#define ALIRSNANALYSISPHI7TEVNOPID_H

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

class AliRsnAnalysisPhi7TeVNoPID : public AliAnalysisTaskSE
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

    AliRsnAnalysisPhi7TeVNoPID(const char *name = "Phi7TeVNoPID");
    AliRsnAnalysisPhi7TeVNoPID(const AliRsnAnalysisPhi7TeVNoPID& copy);
    AliRsnAnalysisPhi7TeVNoPID& operator=(const AliRsnAnalysisPhi7TeVNoPID& copy);
    virtual ~AliRsnAnalysisPhi7TeVNoPID();

    void             SetUseMC(Bool_t yn = kTRUE) {fUseMC = yn;}
    
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

    Bool_t   fUseMC;        // use MC or data?
    
    Short_t  fPDG;          // PDG code
    Short_t  fCh;           // control flag for like/unlike sign
    Short_t  fITS[2];       // check flag to know if one or both candidates are ITS standalone
    Float_t  fIM;           // inv mass
    Float_t  fPt;           // transv momentum
    Float_t  fY;            // rapidity
    Float_t  fEta;          // pseudo-rapidity
    Float_t  fTPCnsigma[2]; // number of sigma in TPC
    Float_t  fITSnsigma[2]; // number of sigma in ITS
    Float_t  fTOFdiff[2];   // relative PID signal in TOF
    Float_t  fP[2];         // total momentum at vertex
    Float_t  fPTPC[2];      // total momentum at inner TPC wall
    
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
    ClassDef(AliRsnAnalysisPhi7TeVNoPID,1)
};

#endif
