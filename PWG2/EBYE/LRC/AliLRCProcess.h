// This class is creatig TH2D histogramms for Nch - Nch , Nch - Pt , Pt - Pt 
// dirtributions for given ETA windows and some supplementary data for Long Range Correlation (LRC) analysis .  
// Class is designid to work with AliAnalysisTaskLRC

// Author : Andrey Ivanov , St.Peterburg State University
// Email: Andrey.Ivanov@cern.ch

/*  See cxx source for full Copyright notice */


#ifndef ALILRCPROCESS_H
#define ALILRCPROCESS_H



// Include
#include "TString.h"
#include "AliLRCBase.h"

class TH1F;
class TH1D;
class TH2D;
class TProfile;
class TList;


class AliLRCProcess: public AliLRCBase
{
public:
// Constructors 
AliLRCProcess();

AliLRCProcess(Double_t _StartForwardETA,Double_t _EndForwardETA,Double_t _StartBakwardETA,Double_t _EndBakwardETA); // Constructor with window setup makes ready-to-run processor

Bool_t InitDataMembers(); //Is to be called in CreateOutputObjects method

// Destructor 
virtual ~AliLRCProcess(); 

 // Setters 
virtual  void SetForwardWindow(Double_t StartETA,Double_t EndETA);  //Sets Forward ETA window
virtual void SetBackwardWindow(Double_t StartETA,Double_t EndETA);  // Sets Backward window
virtual  void SetETAWindows(Double_t _StartForwardETA,Double_t _EndForwardETA,Double_t _StartBakwardETA,Double_t _EndBakwardETA); // Sets both forward and backward windows
virtual  void SetForwardWindowPhi(Double_t StartForwardPhi,Double_t EndForwardPhi){fStartForwardPhi=StartForwardPhi;fEndForwardPhi=EndForwardPhi;}
virtual void SetBackwardWindowPhi(Double_t StartBackwardPhi,Double_t EndBackwardPhi){fStartBackwardPhi=StartBackwardPhi;fEndBackwardPhi=EndBackwardPhi;}

virtual  void SetHistPtRange(Double_t LoPt,Double_t HiPt,Int_t PtBins);  // Sets range for Pt histos axis
virtual  void SetHistMultRange(Int_t LoMult,Int_t HiMult,Int_t MultBins=0); // Sets range for Nch histos axis
virtual  void SetOutputSlotNumber(Int_t SlotNumber); // Sets number of output slot for LRCProcessor

// Getters
virtual TList* CreateOutput() const ; // Creates output object
virtual TString GetShortDef() const ; // Returns fShortDef value
virtual Int_t GetOutputSlotNumber() const; // Returns number of output slot for LRCProcessor

// Track by track event import
virtual void StartEvent();  // Open new Event for track by track event import
virtual void AddTrackPtEta(Double_t Pt, Double_t Eta , Double_t Phi=0.1); // Imports track to the event
virtual void AddTrackForward(Double_t Pt, Double_t Eta , Double_t Phi=0.1);  // Imports track to the event directly to Forward window
virtual void AddTrackBackward(Double_t Pt, Double_t Eta , Double_t Phi=0.1);  // Imports track to the event directly to Backward window
virtual void FinishEvent(); // Close opened event

// 



private:

virtual  void SetShortDef();  // Sets fShortDef according to window paramiters

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
  Double_t fStartForwardPhi;  // Forward phi angle interval
  Double_t fEndForwardPhi;   // Forward phi angle interval
  
  
  Double_t fStartBakwardETA;  // Backward window lover rapidity
  Double_t fEndBakwardETA;    // Backward window higer rapidity
  
  Double_t fStartBackwardPhi;  // Backward window phi angle interval
  Double_t fEndBackwardPhi;    // Backward window phi angle interval
  
  

  Double_t fHiPt;		// Max Pt for histos
  Double_t fLoPt; 		// Min Pt for histos
  Int_t fHiMult;		// Min multiplicity for histos
  Int_t fLoMult;		// Max multiplicity for histos
  
  Int_t fMultBins;		// N bins for multiplicity
  Int_t fPtBins;		// N bins for Pt
  

// Track by track import values
Double_t fSumPtFw; 		// Summary Pt forward
Double_t fSumPtBw; 		// Summary Pt backward
Double_t fSumPtBw2;		// Summary Pt^2 backward
Int_t fNchFw; 			// Number of tracks Forward
Int_t fNchBw; 			// Number of tracks Backward




//Output data
TList* fOutList;    //! List of output data

TString fShortDef; // Short desctiption 
Int_t fOutputSlot; // Output slot number for this Processor

  // Total spectras (debugging for TAG selection)
  TH1F        *fHistPt; //! Overal Pt spectrum
  TH1F        *fHistEta; //! Overal Eta spectrum

 
 // Output histogramms -----------------------------------

 
  TH2D* fHistNN;        //! N-N 2D Profile
  TH2D* fHistPtN;	//! Pt-N 2D Profile
  TH2D* fHistPtPt;	//! Pt-Pt 2D Profile
  TProfile* fProfNberr;	//! Nbackward error Profile PtN
  TProfile* fProfNberrPtPt;	//! Nbackward error Profile PtPt
  TProfile* fProfdPtB;  //! Used to store (in first bin) summary of PtB and its std diviation
  TProfile* fProfTestLRC; //! Diognostic LRC Pt - N correlation

  // Supp. info for windows
  //Forward
  TH1D* fHistPtForward;   //! Pt spectrum in Forward windows
  TH1D* fHistEtaForward;  //! Eta spectrum in Forward windows
  TH1D* fHistNchForward;  //! Nch spectrum in Forward windows
  TH1D* fHistPhiForward;  //! Phi spectrum in Forward windows
  
   //Bakward
  TH1D* fHistPtBakward;   //! Pt spectrum in Bakward windows
  TH1D* fHistEtaBakward;  //! Eta spectrum in Bakward windows
  TH1D* fHistNchBakward;  //! Nch spectrum in Bakward windows
  TH1D* fHistPhiBakward;  //! Phi spectrum in Bakward windows
  
  
         AliLRCProcess(const AliLRCProcess&); // not implemented
        AliLRCProcess& operator=(const AliLRCProcess&); // not implemented


ClassDef(AliLRCProcess, 1);
};

#endif

