#ifndef ALIANALYSISTASKSIGMA1385_H
#define ALIANALYSISTASKSIGMA1385_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskSigma1385 class
//-----------------------------------------------------------------

class TList;
class TH1F;
class TNtuple;
class AliESDcascade;

//-------------------------------- For PID


#include  "AliESDpid.h"
#include  "AliTOFT0maker.h"
#include  "AliTOFcalib.h"
#include  "AliCDBManager.h"
#include  "AliESDtrackCuts.h"

//-------------------------------

#include "TString.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSigma1385 : public AliAnalysisTaskSE {
public:
   AliAnalysisTaskSigma1385();
   AliAnalysisTaskSigma1385(const char *name);
   virtual ~AliAnalysisTaskSigma1385() {}

   virtual void   UserCreateOutputObjects();
   virtual void   UserExec(Option_t *option);
   virtual void   Terminate(Option_t *);

   void SetCollidingSystems(Short_t collidingSystems = 0)     {fCollidingSystems = collidingSystems;}
   void SetAnalysisType(const char* analysisType = "ESD") {fAnalysisType = analysisType;}
   void SetDataType(const char* dataType = "REAL") {fDataType = dataType;}

   virtual const char *GetDataType() {return fDataType;}

   Bool_t *IsSelected(AliESDtrack* track);




//-------------------------------- For Alberto's PID


   void             SetCheckITS(Bool_t yn = kTRUE) {fCheckITS = yn;}
   void             SetCheckTPC(Bool_t yn = kTRUE) {fCheckTPC = yn;}
   void             SetCheckTOF(Bool_t yn = kTRUE) {fCheckTOF = yn;}
   void             SetUseGlobal(Bool_t yn = kTRUE) {fUseGlobal = yn;}
   void             SetUseITSSA(Bool_t yn = kTRUE) {fUseITSSA = yn;}

   void             SetITSband(Double_t v) {fMaxITSband = v;}

   void             SetTPCpLimit(Double_t v) {fTPCpLimit = v;}
   void             SetTPCrange(Double_t min, Double_t max) {fMinTPCband = min; fMaxTPCband = max;}
   void             SetTPCpar(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4)
   {fTPCpar[0] = p0; fTPCpar[1] = p1; fTPCpar[2] = p2; fTPCpar[3] = p3; fTPCpar[4] = p4;}

   void             SetTOFcalibrateESD(Bool_t yn = kTRUE)  {fTOFcalibrateESD = yn;}
   void             SetTOFcorrectTExp(Bool_t yn = kTRUE)  {fTOFcorrectTExp = yn;}
   void             SetTOFuseT0(Bool_t yn = kTRUE)  {fTOFuseT0 = yn;}
   void             SetTOFtuneMC(Bool_t yn = kTRUE)  {fTOFtuneMC = yn;}
   void             SetTOFresolution(Double_t v = 100.0) {fTOFresolution = v;}
   void             SetTOFrange(Double_t v1, Double_t v2) {fMinTOF = v1; fMaxTOF = v2;}


//-------------------------------

private:

   TString fAnalysisType;            // "ESD" or "AOD" analysis type
   Short_t fCollidingSystems;        // 0 = pp collisions or 1 = AA collisions
   TString fDataType;          // "REAL" or "SIM" data type

   TList  *fListHistCascade;      //! List of Cascade histograms
   TH1F    *fHistEventMultiplicity;
   TH1F    *fHistEventMultiplicityRAVS;  //event rejected after vertex selection
   TNtuple *fNtuple1;
   TNtuple *fNtuple2;
   TNtuple *fNtuple3;
   TNtuple *fNtuple4;



//-------------------------------- For PID

protected:

   Bool_t           isMC;
   Bool_t           fIsMC;             //  switch for MC analysis
   Bool_t           fCheckITS;         //  switch for ITS dE/dx check
   Bool_t           fCheckTPC;         //  switch for TPC dE/dx check
   Bool_t           fCheckTOF;         //  switch for TOF time check
   Bool_t           fUseGlobal;        //  switch to use TPC global tracks
   Bool_t           fUseITSSA;         //  switch to use ITS standalone tracks

   Double_t         fMaxITSband;       //  range for ITS de/dx band

   Double_t         fTPCpLimit;        //  limit to choose what band to apply
   Double_t         fTPCpar[5];        //  parameters for TPC bethe-Bloch
   Double_t         fMinTPCband;       //  range for TPC de/dx band - min
   Double_t         fMaxTPCband;       //  range for TPC de/dx band - max


   AliESDpid       *fESDpid;           //! PID manager
   AliTOFT0maker   *fTOFmaker;         //! TOF time0 computator
   AliTOFcalib     *fTOFcalib;         //! TOF calibration
   Bool_t           fTOFcalibrateESD;  //  TOF settings
   Bool_t           fTOFcorrectTExp;   //  TOF settings
   Bool_t           fTOFuseT0;         //  TOF settings
   Bool_t           fTOFtuneMC;        //  TOF settings
   Double_t         fTOFresolution;    //  TOF settings
   Double_t         fMinTOF;           //  range for TOF PID (min)
   Double_t         fMaxTOF;           //  range for TOF PID (max)
   Int_t            fLastRun;          //  last run number

//-------------------------------


   AliAnalysisTaskSigma1385(const AliAnalysisTaskSigma1385&);            // not implemented
   AliAnalysisTaskSigma1385& operator=(const AliAnalysisTaskSigma1385&); // not implemented

   ClassDef(AliAnalysisTaskSigma1385, 3);
};

#endif
