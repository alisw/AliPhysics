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
   
   // flow of resonances:
   // a) methods:
   virtual void BookObjectsForFlowOfResonances();
   virtual void FlowOfResonances(AliFlowEventSimple* anEvent);
   // b) setters and getters:
   void SetResonanceList(TList* const rl) {this->fResonanceList = rl;}
   TList* GetResonanceList() const {return this->fResonanceList;}    
   void SetFlowOfResonances(Bool_t const ffor) {this->fFlowOfResonances = ffor;};
   Bool_t GetFlowOfResonances() const {return this->fFlowOfResonances;};   
   void SetResonanceSettings(TProfile* const rs) {this->fResonanceSettings = rs;};
   TProfile* GetResonanceSettings() const {return this->fResonanceSettings;};
   void SetPairCorrelator(TProfile* const spc, Int_t const cs) {this->fPairCorrelator[cs] = spc;};
   TProfile* GetPairCorrelator(Int_t cs) const {return this->fPairCorrelator[cs];};
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
   TProfile2D*  fHistProDiffFlowPtEtaRP;  // profile used to calculate the differential flow (Pt,Eta) of RP particles
   TProfile*    fHistProDiffFlowPtRP;     // profile used to calculate the differential flow (Pt) of RP particles 
   TProfile*    fHistProDiffFlowEtaRP;    // profile used to calculate the differential flow (Eta) of RP particles 
   TProfile2D*  fHistProDiffFlowPtEtaPOI; // profile used to calculate the differential flow (Pt,Eta) of POI particles
   TProfile*    fHistProDiffFlowPtPOI;    // profile used to calculate the differential flow (Pt) of POI particles 
   TProfile*    fHistProDiffFlowEtaPOI;   // profile used to calculate the differential flow (Eta) of POI particles
   TH1D*        fHistSpreadOfFlow;        // histogram filled with reference flow calculated e-b-e    
   Int_t        fHarmonic;                // harmonic 
   
   // objects needed for a study of flow of resonances (and in addition for strong parity violation):
   TList *fResonanceList; // list to hold all objects relevant for a study of flow of resonances 
   Bool_t fFlowOfResonances; // evaluate and store objects relevant for study of flow of resonances
   TProfile *fResonanceSettings; // profile used to hold all flags relevant for the flow of resonances
   TProfile *fPairCorrelator[2]; // profiles used to calculate <cos[n(phi_{pair}-RP)]> and <sin[n(phi_{pair}-RP)]> (0 = cos, 1 = sin), where phi_{pair} = x*phi1+(1-x)*phi2
   Double_t fXinPairAngle; // x in definition phi_{pair} = x*phi1+(1-x)*phi2
                                       
   ClassDef(AliFlowAnalysisWithMCEventPlane,1)  // Analyse particle distribution versus MC reaction plane
     };

     
#endif


