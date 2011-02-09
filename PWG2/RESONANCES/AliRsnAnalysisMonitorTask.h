//
// Header file for implementation of data analysis aft 900 GeV
//
// Author: A. Pulvirenti
//

#ifndef ALIRSNANALYSISMONITORTASK_H
#define ALIRSNANALYSISMONITORTASK_H

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

class AliRsnAnalysisMonitorTask : public AliAnalysisTaskSE {
public:

   AliRsnAnalysisMonitorTask(const char *name = "Phi7TeV");
   AliRsnAnalysisMonitorTask(const AliRsnAnalysisMonitorTask& copy);
   AliRsnAnalysisMonitorTask& operator=(const AliRsnAnalysisMonitorTask& copy);
   virtual ~AliRsnAnalysisMonitorTask();

   void             SetTPCpar(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4)
   {fTPCpar[0] = p0; fTPCpar[1] = p1; fTPCpar[2] = p2; fTPCpar[3] = p3; fTPCpar[4] = p4;}

   void             SetTOFcalibrateESD(Bool_t yn = kTRUE)  {fTOFcalibrateESD = yn;}
   void             SetTOFcorrectTExp(Bool_t yn = kTRUE)  {fTOFcorrectTExp = yn;}
   void             SetTOFuseT0(Bool_t yn = kTRUE)  {fTOFuseT0 = yn;}
   void             SetTOFtuneMC(Bool_t yn = kTRUE)  {fTOFtuneMC = yn;}
   void             SetTOFresolution(Double_t v = 100.0) {fTOFresolution = v;}

   virtual void     UserCreateOutputObjects();
   virtual void     UserExec(Option_t *option = "");
   virtual void     Terminate(Option_t *option = "");

   Int_t            EventEval(AliESDEvent *esd);
   Bool_t           IsTPCtrack(AliESDtrack *track);
   Bool_t           IsITSSAtrack(AliESDtrack *track);
   void             ProcessESD(AliESDEvent *esd, const AliESDVertex *v, AliStack *stack);

   AliRsnCutSet*    GetEventCuts() {return &fEventCuts;}
   AliRsnCutSet*    GetTrackCuts() {return &fTrackCuts;}

private:

   TTree              *fOut;              //  output TTree
   AliRsnMonitorTrack *fTrack;            //  branch object for output TTree

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
   ClassDef(AliRsnAnalysisMonitorTask, 1)
};

#endif
