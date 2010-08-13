/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
// Description: Maker to analyze Flow from the generated MC reaction plane.
//              This class is used to get the real value of the flow 
//              to compare the other methods to when analysing simulated events.

/* $Id$ */

#ifndef ALIFLOWANALYSISWITHMCEVENTPLANE_H
#define ALIFLOWANALYSISWITHMCEVENTPLANE_H

class TVector2;
class TString;
class TDirectoryFile;

class AliFlowTrackSimple;
class AliFlowEventSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;

class TH1F;
class TH1D;
class TProfile;
class TProfile2D;
class TObjArray;
class TFile;
class TList;
class TComplex;
class Riostream;


 
class AliFlowAnalysisWithMCEventPlane {

 public:
 
   AliFlowAnalysisWithMCEventPlane();            //default constructor
   virtual  ~AliFlowAnalysisWithMCEventPlane();  //destructor
 
   void      WriteHistograms(TString* outputFileName);
   void      WriteHistograms(TString outputFileName);
   void      WriteHistograms(TDirectoryFile *outputFileName);
   void      Init();                                       //defines variables and histograms
   void      Make(AliFlowEventSimple* anEvent);            //calculates variables and fills histograms
   void      GetOutputHistograms(TList *outputListHistos); //get pointers to all output histograms (called before Finish()) 
   void      Finish();                                     //saves histograms
   
   void      SetDebug(Bool_t kt)          { this->fDebug = kt ; }
   Bool_t    GetDebug() const             { return this->fDebug ; }

   void      SetEventNumber(Int_t n)      { this->fEventNumber = n; }
   Int_t     GetEventNumber() const       { return this->fEventNumber; }

   // Output 
   TList*    GetHistList() const          { return this->fHistList ; }  
   AliFlowCommonHist* GetCommonHists() const  { return this->fCommonHists; }
   void      SetCommonHists(AliFlowCommonHist* const aCommonHist)  
     { this->fCommonHists = aCommonHist; }
   AliFlowCommonHistResults*  GetCommonHistsRes() const { return this->fCommonHistsRes; }
   void      SetCommonHistsRes( AliFlowCommonHistResults* const aCommonHistResult ) 
     { this->fCommonHistsRes = aCommonHistResult; }
   
   //histograms
   TH1F*     GetHistRP() const            {return this->fHistRP; } 
   void      SetHistRP(TH1F* const  aHistRP)     {this->fHistRP = aHistRP; }
   
   TProfile* GetHistProIntFlow() const    {return this->fHistProIntFlow; } 
   void      SetHistProIntFlow(TProfile* const aHistProIntFlow) 
     {this->fHistProIntFlow = aHistProIntFlow; }
     
   TProfile* GetHistProIntFlowVsM() const    {return this->fHistProIntFlowVsM; } 
   void      SetHistProIntFlowVsM(TProfile* const aHistProIntFlowVsM) 
     {this->fHistProIntFlowVsM = aHistProIntFlowVsM; }

   TProfile2D* GetHistProDiffFlowPtEtaRP() const    {return this->fHistProDiffFlowPtEtaRP; } 
   void      SetHistProDiffFlowPtEtaRP(TProfile2D* const aHistProDiffFlowPtEtaRP) 
     {this->fHistProDiffFlowPtEtaRP = aHistProDiffFlowPtEtaRP; }   
   
   TProfile* GetHistProDiffFlowPtRP() const    {return this->fHistProDiffFlowPtRP; } 
   void      SetHistProDiffFlowPtRP(TProfile* const aHistProDiffFlowPtRP) 
     {this->fHistProDiffFlowPtRP = aHistProDiffFlowPtRP; } 
   
   TProfile* GetHistProDiffFlowEtaRP() const   {return this->fHistProDiffFlowEtaRP; } 
   void      SetHistProDiffFlowEtaRP(TProfile* const aHistProDiffFlowEtaRP) 
     {this->fHistProDiffFlowEtaRP = aHistProDiffFlowEtaRP; } 
     
   TProfile2D* GetHistProDiffFlowPtEtaPOI()const     {return this->fHistProDiffFlowPtEtaPOI; } 
   void      SetHistProDiffFlowPtEtaPOI(TProfile2D* const aHistProDiffFlowPtEtaPOI) 
     {this->fHistProDiffFlowPtEtaPOI = aHistProDiffFlowPtEtaPOI; }   
   
   TProfile* GetHistProDiffFlowPtPOI()const    {return this->fHistProDiffFlowPtPOI; } 
   void      SetHistProDiffFlowPtPOI(TProfile* const aHistProDiffFlowPtPOI) 
     {this->fHistProDiffFlowPtPOI = aHistProDiffFlowPtPOI; } 
   
   TProfile* GetHistProDiffFlowEtaPOI()const   {return this->fHistProDiffFlowEtaPOI; } 
   void      SetHistProDiffFlowEtaPOI(TProfile* const aHistProDiffFlowEtaPOI) 
     {this->fHistProDiffFlowEtaPOI = aHistProDiffFlowEtaPOI; } 
     
   TH1D* GetHistSpreadOfFlow()const   {return this->fHistSpreadOfFlow; } 
   void      SetHistSpreadOfFlow(TH1D* const aHistSpreadOfFlow) 
     {this->fHistSpreadOfFlow = aHistSpreadOfFlow; }    
   
   // harmonic:
   void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;};
   Int_t GetHarmonic() const {return this->fHarmonic;};
   
   // mixed harmonics:
   // a) methods:
   virtual void InitalizeArraysForMixedHarmonics();
   virtual void BookObjectsForMixedHarmonics();
   virtual void EvaluateMixedHarmonics(AliFlowEventSimple* anEvent);
   virtual void GetOutputHistoramsForMixedHarmonics(TList *mixedHarmonicsList);
   // b) setters and getters:
   void SetMixedHarmonicsList(TList* const mhl) {this->fMixedHarmonicsList = mhl;}
   TList* GetMixedHarmonicsList() const {return this->fMixedHarmonicsList;}    
   void SetEvaluateMixedHarmonics(Bool_t const emh) {this->fEvaluateMixedHarmonics = emh;};
   Bool_t GetEvalauteMixedHarmonics() const {return this->fEvaluateMixedHarmonics;};   
   void SetMixedHarmonicsSettings(TProfile* const mhs) {this->fMixedHarmonicsSettings = mhs;};
   TProfile* GetMixedHarmonicsSettings() const {return this->fMixedHarmonicsSettings;};  
   void SetPairCorrelator(TProfile* const spc, Int_t const cs) {this->fPairCorrelator[cs] = spc;};
   TProfile* GetPairCorrelator(Int_t const cs) const {return this->fPairCorrelator[cs];};
   void SetPairCorrelatorVsM(TProfile* const spcVsM, Int_t const cs) {this->fPairCorrelatorVsM[cs] = spcVsM;};
   TProfile* GetPairCorrelatorVsM(Int_t const cs) const {return this->fPairCorrelatorVsM[cs];};   
   void SetnBinsMult(Int_t const nbm) {this->fnBinsMult = nbm;};
   Int_t GetnBinsMult() const {return this->fnBinsMult;};  
   void SetMinMult(Double_t const minm) {this->fMinMult = minm;};
   Double_t GetMinMult() const {return this->fMinMult;};
   void SetMaxMult(Double_t const maxm) {this->fMaxMult = maxm;};
   Double_t GetMaxMult() const {return this->fMaxMult;};   
   void SetPairCorrelatorVsPtSumDiff(TProfile* const spcVspsd, Int_t const cs, Int_t const sd) {this->fPairCorrelatorVsPtSumDiff[cs][sd] = spcVspsd;};
   TProfile* GetPairCorrelatorVsPtSumDiff(Int_t const cs, Int_t const sd) const {return this->fPairCorrelatorVsPtSumDiff[cs][sd];};
   void SetNinCorrelator(Int_t const n) {this->fNinCorrelator = n;};
   Int_t GetNinCorrelator() const {return this->fNinCorrelator;};     
   void SetMinCorrelator(Int_t const m) {this->fMinCorrelator = m;};
   Int_t GetMinCorrelator() const {return this->fMinCorrelator;};     
   void SetXinPairAngle(Double_t const xipa) {this->fXinPairAngle = xipa;};
   Double_t GetXinPairAngle() const {return this->fXinPairAngle;};   
  
 private:
 
   AliFlowAnalysisWithMCEventPlane(const AliFlowAnalysisWithMCEventPlane& aAnalysis);             //copy constructor
   AliFlowAnalysisWithMCEventPlane& operator=(const AliFlowAnalysisWithMCEventPlane& aAnalysis);  //assignment operator 

      
#ifndef __CINT__
   TVector2*    fQsum;              // flow vector sum
   Double_t     fQ2sum;             // flow vector sum squared
#endif /*__CINT__*/

   Int_t        fEventNumber;       // event counter
   Bool_t       fDebug ;            //! flag for lyz analysis: more print statements

   TList*       fHistList;          //list to hold all output histograms  
    
   AliFlowCommonHist* fCommonHists;              // hist
   AliFlowCommonHistResults* fCommonHistsRes;    // hist
   
   TH1F*        fHistRP;                  // reaction plane
   TProfile*    fHistProIntFlow;          // profile used to calculate the integrated flow of RP particles
   TProfile*    fHistProIntFlowVsM;       // profile used to calculate the integrated flow of RP particles vs multiplicity
   TProfile2D*  fHistProDiffFlowPtEtaRP;  // profile used to calculate the differential flow (Pt,Eta) of RP particles
   TProfile*    fHistProDiffFlowPtRP;     // profile used to calculate the differential flow (Pt) of RP particles 
   TProfile*    fHistProDiffFlowEtaRP;    // profile used to calculate the differential flow (Eta) of RP particles 
   TProfile2D*  fHistProDiffFlowPtEtaPOI; // profile used to calculate the differential flow (Pt,Eta) of POI particles
   TProfile*    fHistProDiffFlowPtPOI;    // profile used to calculate the differential flow (Pt) of POI particles 
   TProfile*    fHistProDiffFlowEtaPOI;   // profile used to calculate the differential flow (Eta) of POI particles
   TH1D*        fHistSpreadOfFlow;        // histogram filled with reference flow calculated e-b-e    
   Int_t        fHarmonic;                // harmonic 
   
   // mixed harmonics:
   TList *fMixedHarmonicsList; // list to hold all objects relevant for mixed harmonics 
   Bool_t fEvaluateMixedHarmonics; // evaluate and store objects relevant for mixed harmonics
   TProfile *fMixedHarmonicsSettings; // profile used to hold all flags relevant for the mixed harmonics
   TProfile *fPairCorrelator[2]; // profiles used to calculate <cos[m*phi_{pair}-n*RP]> and <sin[m*phi_{pair}-n*RP]> (0 = cos, 1 = sin), where phi_{pair} = x*phi1+(1-x)*phi2
   TProfile *fPairCorrelatorVsM[2]; // <cos[m*phi_{pair}-n*RP]> and <sin[m*phi_{pair}-n*RP]> versus multiplicity (0 = cos, 1 = sin), where phi_{pair} = x*phi1+(1-x)*phi2
   Int_t fnBinsMult; // number of multiplicity bins for mixed harmonics analysis versus multiplicity  
   Double_t fMinMult; // minimal multiplicity for mixed harmonics analysis versus multiplicity  
   Double_t fMaxMult; // maximal multiplicity for mixed harmonics analysis versus multiplicity    
   TProfile *fPairCorrelatorVsPtSumDiff[2][2]; // <cos/sin[m*phi_{pair}-n*RP]> vs (1/2)(pt1+pt2) (0) and |pt1-pt2| (1), where phi_{pair} = x*phi1+(1-x)*phi2
   Int_t fNinCorrelator; // n in <cos[m*phi_{pair}-n*RP]> and <sin[m*phi_{pair}-n*RP]>, where phi_{pair} = x*phi1+(1-x)*phi2
   Int_t fMinCorrelator; // m in <cos[m*phi_{pair}-n*RP]> and <sin[m*phi_{pair}-n*RP]>, where phi_{pair} = x*phi1+(1-x)*phi2   
   Double_t fXinPairAngle; // x in definition phi_{pair} = x*phi1+(1-x)*phi2
                                       
   ClassDef(AliFlowAnalysisWithMCEventPlane,0)  // Analyse particle distribution versus MC reaction plane
  
  };

     
#endif


