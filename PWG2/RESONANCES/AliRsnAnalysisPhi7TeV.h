//
// Header file for implementation of data analysis aft 900 GeV
//
// Author: A. Pulvirenti
//

#ifndef ALIRSNANALYSISPHI7TEV_H
#define ALIRSNANALYSISPHI7TEV_H

#include "AliRsnDaughter.h"
#include "AliRsnCutESD2010.h"
#include "AliAnalysisTaskSE.h"
#include "AliRsnTOFT0maker.h"
#include "AliESDtrackCuts.h"
#include "AliPID.h"

class TH1I;
class TH1F;
class TH3F;
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
      kEmptyEvent              = 5
    };

    AliRsnAnalysisPhi7TeV(const char *name = "Phi7TeV", Bool_t isMC = kFALSE);
    AliRsnAnalysisPhi7TeV(const AliRsnAnalysisPhi7TeV& copy);
    AliRsnAnalysisPhi7TeV& operator=(const AliRsnAnalysisPhi7TeV& copy);
    virtual ~AliRsnAnalysisPhi7TeV();

    void             SetUseMC   (Bool_t yn = kTRUE);
    void             SetCheckITS(Bool_t yn = kTRUE) {fCheckITS = yn;}
    void             SetCheckTPC(Bool_t yn = kTRUE) {fCheckTPC = yn;}
    void             SetCheckTOF(Bool_t yn = kTRUE) {fCheckTOF = yn;}
    void             SetAddITSSA(Bool_t yn = kTRUE) {fAddITSSA = yn;}
    
    void             SetMaxVz(Double_t v)   {fMaxVz = v;}
    
    void             SetITSband(Double_t v) {fMaxITSband = v;}
    void             SetITSmom (Double_t v) {fMaxITSmom  = v;}
    
    void             SetTPCpLimit(Double_t v) {fTPCpLimit = v;}
    void             SetTPCrange(Double_t min, Double_t max) {fMinTPCband = min; fMaxTPCband = max;}
    void             SetTPCpar(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4)
                       {fTPCpar[0]=p0;fTPCpar[1]=p1;fTPCpar[2]=p2;fTPCpar[3]=p3;fTPCpar[4]=p4;}

    void             SetTOFcalibrateESD(Bool_t yn = kTRUE)  {fTOFcalibrateESD = yn;}
    void             SetTOFcorrectTExp (Bool_t yn = kTRUE)  {fTOFcorrectTExp = yn;}
    void             SetTOFuseT0       (Bool_t yn = kTRUE)  {fTOFuseT0 = yn;}
    void             SetTOFtuneMC      (Bool_t yn = kTRUE)  {fTOFtuneMC = yn;}
    void             SetTOFresolution  (Double_t v = 100.0) {fTOFresolution = v;}
    void             SetMinTOF         (Double_t v)         {fMinTOF = v;}
    void             SetMaxTOF         (Double_t v)         {fMaxTOF = v;}

    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option = "");
    virtual void     Terminate(Option_t *option = "");
    
    Int_t            EventEval(AliESDEvent *esd);
    AliESDtrackCuts* GetCutsTPC() {return &fESDtrackCutsTPC;}
    AliESDtrackCuts* GetCutsITS() {return &fESDtrackCutsITS;}
    
    Bool_t           IsITSTPC (AliESDtrack *track);
    Bool_t           IsITSSA  (AliESDtrack *track);
    Bool_t           MatchTOF (AliESDtrack *track);
    Bool_t           OkQuality(AliESDtrack *track);
    Bool_t           OkITSPID (AliESDtrack *track, AliPID::EParticleType pid);
    Bool_t           OkTPCPID (AliESDtrack *track, AliPID::EParticleType pid);
    Bool_t           OkTOFPID (AliESDtrack *track, AliPID::EParticleType pid);
    Bool_t           OkTrack  (AliESDtrack *track, AliPID::EParticleType pid);

  private:

    void     ProcessESD(AliESDEvent *esd, AliStack *stack);

    Bool_t   fUseMC;      // use MC or data?
    Bool_t   fCheckITS;   // chec ITS PID?
    Bool_t   fCheckTPC;   // chec TPC PID?
    Bool_t   fCheckTOF;   // chec TOF PID?
    Bool_t   fAddITSSA;   // add ITS standalone?
    
    Double_t fMaxVz;      // range in Z of primary vertex w.r. to origin
    
    Double_t fMaxITSband; // range for ITS de/dx band
    Double_t fMaxITSmom;  // maximum momentum for ITS identification

    Double_t fTPCpLimit;  // limit to choose what band to apply
    Double_t fTPCpar[5];  // parameters for TPC bethe-Bloch
    Double_t fMinTPCband; // range for TPC de/dx band - min
    Double_t fMaxTPCband; // range for TPC de/dx band - max
    Double_t fMinTOF;     // TOF range (min)
    Double_t fMaxTOF;     // TOF range (min)

    TList     *fOutList;        // list for monitoring histograms
    TH3F      *fUnlike;         // unlike-sign pairs
    TH3F      *fLikePP;         // unlike-sign pairs
    TH3F      *fLikeMM;         // unlike-sign pairs
    TH3F      *fTrues;          // true pairs
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
    
    AliRsnDaughter   fDaughter;
    AliRsnCutESD2010 fRsnCuts;

    // ROOT dictionary
    ClassDef(AliRsnAnalysisPhi7TeV,1)
};

inline Bool_t AliRsnAnalysisPhi7TeV::IsITSTPC(AliESDtrack *vtrack)
{
//
// Checks if the track has the status flags required for a global track
//

  if (!vtrack)
  {
    AliWarning("NULL argument: impossible to check status");
    return kFALSE;
  }
  
  return vtrack->IsOn(AliESDtrack::kTPCin);
}

inline Bool_t AliRsnAnalysisPhi7TeV::IsITSSA(AliESDtrack *vtrack)
{
//
// Checks if the track has the status flags required for an ITS standalone track
//

  if (!vtrack)
  {
    AliWarning("NULL argument: impossible to check status");
    return kFALSE;
  }
  
  Bool_t isTPCin     = vtrack->IsOn(AliESDtrack::kTPCin);
  Bool_t isITSrefit  = vtrack->IsOn(AliESDtrack::kITSrefit);
  Bool_t isITSpureSA = vtrack->IsOn(AliESDtrack::kITSpureSA);
  Bool_t isITSpid    = vtrack->IsOn(AliESDtrack::kITSpid);
  
  return ( (!isTPCin) && isITSrefit && (!isITSpureSA) && isITSpid );
}


inline Bool_t AliRsnAnalysisPhi7TeV::MatchTOF(AliESDtrack *vtrack)
{
//
// Checks if the track has matched the TOF detector
//

  if (!vtrack)
  {
    AliWarning("NULL argument: impossible to check status");
    return kFALSE;
  }
  
  // require a minimum length to have meaningful match
  if (vtrack->GetIntegratedLength() < 350.) return kFALSE;
  
  Bool_t isTOFout = vtrack->IsOn(AliESDtrack::kTOFout);
  Bool_t isTIME   = vtrack->IsOn(AliESDtrack::kTIME);
  
  return ( isTOFout && isTIME );
}

#endif
