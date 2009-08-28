// This class is creatig TH2D histogramms for Nch - Nch , Nch - Pt , Pt - Pt 
// dirtributions for given ETA windows and some supplementary data for Long Range Correlation (LRC) analysis .  
// Class is designid to work with AliAnalysisTaskLRC

// Author : Andrey Ivanov , St.Peterburg State University
// Email: Andrey.Ivanov@cern.ch

// Version line : 3.5
// Version 3.5.5


#ifndef ALILRCPROCESS_H
#define ALILRCPROCESS_H

// Include
#include "TH1F.h"
#include  "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TList.h"
#include "TString.h"

class AliLRCProcess: public TObject 
{
public:
// Constructors 

AliLRCProcess(Double_t _StartForwardETA=-1.0,Double_t _EndForwardETA=1.0,Double_t _StartBakwardETA=-1.0,Double_t _EndBakwardETA=1.0); // Constructor with window setup makes ready-to-run processor

Bool_t InitDataMembers(); //Is to be called in CreateOutputObjects method

// Destructor 
virtual ~AliLRCProcess(); 

 // Setters 
  void SetForwardWindow(Double_t StartETA,Double_t EndETA);
  void SetBackwardWindow(Double_t StartETA,Double_t EndETA);
  void SetETAWindows(Double_t _StartForwardETA,Double_t _EndForwardETA,Double_t _StartBakwardETA,Double_t _EndBakwardETA);
  void SetPtRange(Double_t LoPt,Double_t HiPt,Int_t PtBins);
  void SetMultRange(Int_t LoMult,Int_t HiMult,Int_t MultBins=0);
  

// Getters
TList* CreateOutput(); // Creates output object
TString GetShortDef(); // Returns fShortDef value

// Track by track event import
void StartEvent();  // Open new Event for track by track event import
void AddTrackPtEta(Double_t Pt, Double_t Eta ); // Imports track to the event
void FinishEvent(); // Close opened event

// 



private:

 void SetShortDef();  // Sets fShortDef according to window paramiters

//Data Init and checks 
//Bool_t StateCheck();   // Check if data is OK
Bool_t fIsEventOpend;  // Indicates opened event 
Bool_t fIsOnline;  // Are data objects created
Bool_t fDisplayInitOnDemandWarning; // Switch warning when InitDataInDemand is called;
//Bool_t InitDataOnDemand(); // Create data objects 

// Statistics
Int_t fEventCount; //Number of processed events

 // Windows paramiters -----------------------------------
  
  Double_t fStartForwardETA;  // Forward windos lover rapidity
  Double_t fEndForwardETA;    // Forward window higer rapidity	
  Double_t fStartBakwardETA;  // Bakward window lover rapidity
  Double_t fEndBakwardETA;    // Bakward window higer rapidity

  Double_t fHiPt;		// Max Pt for histos
  Double_t fLoPt; 		// Min Pt for histos
  Int_t fHiMult;		// Min multiplicity for histos
  Int_t fLoMult;		// Max multiplicity for histos
  
  Int_t fMultBins;		// N bins for multiplicity
  Int_t fPtBins;		// N bins for Pt
  

// Track by track import values
Double_t fSumPtFw; 		// Summary Pt forward
Double_t fSumPtBw; 		// Summary Pt backward
Int_t fNchFw; 			// Number of tracks Forward
Int_t fNchBw; 			// Number of tracks Backward




//Output data
TList* fOutList;    // List of output data

TString fShortDef; // Short desctiption 

  // Total spectras (debugging for TAG selection)
  TH1F        *fHistPt; // Overal Pt spectrum
  TH1F        *fHistEta; // Overal Eta spectrum

 
 // Output histogramms -----------------------------------

 
  TH2D* fHistNN;        // N-N 2D Profile
  TH2D* fHistPtN;	// Pt-N 2D Profile
  TH2D* fHistPtPt;	// Pt-Pt 2D Profile
  TProfile* fProfNberr;	// Nbackward error Profile
  TProfile* fProfdPtB;  // Used to store (in first bin) summary of PtB and its std diviation
  TProfile* fProfTestLRC; // Diognostic LRC Pt - N correlation

  // Supp. info for windows
  //Forward
  TH1D* fHistPtForward;   // Pt spectrum in Forward windows
  TH1D* fHistEtaForward;  // Eta spectrum in Forward windows
  TH1D* fHistNchForward;  // Nch spectrum in Forward windows
  
   //Bakward
  TH1D* fHistPtBakward;   // Pt spectrum in Bakward windows
  TH1D* fHistEtaBakward;  // Eta spectrum in Bakward windows
  TH1D* fHistNchBakward;  // Nch spectrum in Bakward windows
 
 

ClassDef(AliLRCProcess, 1);
};

#endif

