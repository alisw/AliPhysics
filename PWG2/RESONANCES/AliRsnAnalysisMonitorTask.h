//
// Header file for implementation of data analysis aft 900 GeV
//
// Author: A. Pulvirenti
//

#ifndef ALIRSNANALYSISMONITORTASK_H
#define ALIRSNANALYSISMONITORTASK_H

#include "AliAnalysisTaskSE.h"
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

class AliRsnAnalysisMonitorTask : public AliAnalysisTaskSE
{
  public:
  
    AliRsnAnalysisMonitorTask(const char *name = "Phi7TeV");
    AliRsnAnalysisMonitorTask(const AliRsnAnalysisMonitorTask& copy);
    AliRsnAnalysisMonitorTask& operator=(const AliRsnAnalysisMonitorTask& copy);
    virtual ~AliRsnAnalysisMonitorTask();
    
    void             SetITSband(Double_t v) {fMaxITSband = v;}
    
    void             SetTPClargeBandLimit(Double_t v)        {fTPCpLimit = v;}
    void             SetTPCbands(Double_t min, Double_t max) {fLargeTPCband = min; fSmallTPCband = max;}
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
    
    void             EventEval(AliESDEvent *esd);
    Bool_t           IsTPCtrack(AliESDtrack *track);
    Bool_t           IsITSSAtrack(AliESDtrack *track);
    AliESDtrackCuts& GetCutsTPC() {return fESDtrackCutsTPC;}
    AliESDtrackCuts& GetCutsITS() {return fESDtrackCutsITS;}
    void             ProcessESD(AliESDEvent *esd, const AliESDVertex *v, AliStack *stack);

  private:
    
    Int_t            fEventType;           // event classification (0 = vertex with tracks, 1 = vertex with SPD, 2 = bad vertex)
    Double_t         fVertex[3];           // primary vertex position
    Int_t            fNTracks;             // counter for tracks
    
    TTree           *fOut;                 // output TTree
    TClonesArray    *fTracks;              // array of data from tracks
    
    Double_t         fMaxITSband;          // range for ITS de/dx band
    Double_t         fTPCpLimit;           // limit to choose what band to apply
    Double_t         fTPCpar[5];           // parameters for TPC bethe-Bloch
    Double_t         fLargeTPCband;        // range for TPC de/dx band - min
    Double_t         fSmallTPCband;        // range for TPC de/dx band - max
   
    AliESDtrackCuts  fESDtrackCutsTPC;     //  ESD standard defined track cuts for TPC tracks
    AliESDtrackCuts  fESDtrackCutsITS;     //  ESD standard defined track cuts for ITS-SA tracks
    AliESDpid       *fESDpid;              //! PID manager
    
    AliTOFT0maker   *fTOFmaker;            //! TOF time0 computator
    AliTOFcalib     *fTOFcalib;            //! TOF calibration
    Bool_t           fTOFcalibrateESD;     //  TOF settings
    Bool_t           fTOFcorrectTExp;      //  TOF settings
    Bool_t           fTOFuseT0;            //  TOF settings
    Bool_t           fTOFtuneMC;           //  TOF settings
    Double_t         fTOFresolution;       //  TOF settings

    // ROOT dictionary
    ClassDef(AliRsnAnalysisMonitorTask,1)
};

#endif
