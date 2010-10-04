//
// Header file for implementation of data analysis aft 900 GeV
//
// Author: A. Pulvirenti
//

#ifndef ALIRSNANALYSISMONITORPAIRTASK_H
#define ALIRSNANALYSISMONITORPAIRTASK_H

#include "AliAnalysisTaskSE.h"
#include "AliRsnCutSet.h"

class TH1I;
class TH1F;
class TTree;

class AliStack;
class AliESDEvent;
class AliESDVertex;
class AliESDpid;
class AliTOFT0maker;
class AliTOFcalib;

class AliRsnAnalysisMonitorPairTask : public AliAnalysisTaskSE
{
  public:
  
    AliRsnAnalysisMonitorPairTask(const char *name = "Phi7TeV");
    AliRsnAnalysisMonitorPairTask(const AliRsnAnalysisMonitorPairTask& copy);
    AliRsnAnalysisMonitorPairTask& operator=(const AliRsnAnalysisMonitorPairTask& copy);
    virtual ~AliRsnAnalysisMonitorPairTask();
    
    void             SetMasses(Double_t m1, Double_t m2) {fMass[0] = m1; fMass[1] = m2;}
    void             SetInvMassRange(Double_t m1, Double_t m2) {fRangeMin = m1, fRangeMax = m2;}
    
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
    Bool_t           IsTPCtrack(AliESDtrack *track);
    Bool_t           IsITSSAtrack(AliESDtrack *track);
    void             ProcessESD(AliESDEvent *esd, const AliESDVertex *v);
    Bool_t           ProcessTrack(Int_t myIndex, Int_t esdIndex, AliESDEvent *esd, const AliESDVertex *v);
    
    AliRsnCutSet*    GetEventCuts() {return &fEventCuts;}
    AliRsnCutSet*    GetTrackCuts() {return &fTrackCuts;}

  private:
    
    TTree              *fOut;              //  output TTree
    AliRsnMonitorTrack *fTrack[2];         //  branch objects for output TTree
    Double_t            fMass[2];          //  masses assigned to daughters
    Float_t             fInvMass;          //  pair inv mass (computed with above masses)
    Double_t            fRangeMin;         //  minimum accepted invmass
    Double_t            fRangeMax;         //  maximum accepted invmass
    
    Double_t            fTPCpar[5];        //  parameters for TPC bethe-Bloch
   
    AliESDpid          *fESDpid;           //! PID manager
    
    AliTOFT0maker      *fTOFmaker;         //! TOF time0 computator
    AliTOFcalib        *fTOFcalib;         //! TOF calibration
    Bool_t              fTOFcalibrateESD;  //  TOF settings
    Bool_t              fTOFcorrectTExp;   //  TOF settings
    Bool_t              fTOFuseT0;         //  TOF settings
    Bool_t              fTOFtuneMC;        //  TOF settings
    Double_t            fTOFresolution;    //  TOF settings
    
    AliRsnCutSet        fEventCuts;        //  event cuts
    AliRsnCutSet        fTrackCuts;        //  track cuts

    // ROOT dictionary
    ClassDef(AliRsnAnalysisMonitorPairTask,1)
};

#endif
