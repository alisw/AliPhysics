//-------------------------------------------------------------------------
//    Description: 
//    This class is abstact (interface) base class for different classes 
//    used for track-by-track data analysis in AliAnalysisTaskLRC
//    Origin: Andrey Ivanov (SPbSU-CERN) anivanov@cern.ch, 
//    Igor Altsebeev (SPbSU-CERN) 
//-------------------------------------------------------------------------
/*  See cxx source for full Copyright notice */

#ifndef ALILRCBASE_H
#define ALILRCBASE_H

#include <TObject.h>

// Include
class TList;
class AliLRCBase: public TObject
{
public:
// Constructors 
	AliLRCBase(): TObject() {}; // Constructor

// Destructor 
	virtual ~AliLRCBase(){}; 

	virtual Bool_t InitDataMembers() = 0; //Is to be called in CreateOutputObjects method

// Setters 
	virtual  void SetOutputSlotNumber(Int_t SlotNumber) = 0; // Sets number of output slot for LRCProcessor
//	virtual  void SetPrintInfo( Bool_t _flagPrintInfo )  = 0; // Print info flag 

// Getters
	virtual TList* CreateOutput() const 	= 0 ; // Creates output object
	virtual TString GetShortDef() const 	= 0  ; // Returns fShortDef value
	virtual Int_t GetOutputSlotNumber() const = 0 ; // Returns number of output slot for LRCProcessor

// Track by track event import
	virtual void StartEvent() 	= 0;  // Open new Event for track by track event import
	virtual void FinishEvent() 	= 0; // Close opened event

	virtual void SetETAWindows( Double_t _StartForwardETA, Double_t _EndForwardETA, Double_t _StartBackwardETA, Double_t _EndBackwardETA ) = 0; // Sets both forward and backward

	virtual void AddTrackPtEta( Double_t Pt, Double_t Eta, Double_t Phi = 0.1 )	= 0; // Imports track from the event
	virtual void AddTrackForward( Double_t Pt, Double_t Eta, Double_t Phi ) 	= 0; // Imports track to the event directly to Forward window
	virtual void AddTrackBackward( Double_t Pt, Double_t Eta, Double_t Phi ) 	= 0; // Imports track to the event directly to Backward window


private:

	virtual  void SetShortDef()	=0;  // Sets fShortDef according to window paramiters
	
	AliLRCBase(const AliLRCBase&); // not implemented
	AliLRCBase& operator=(const AliLRCBase&); // not implemented


ClassDef(AliLRCBase, 0);
};

#endif

