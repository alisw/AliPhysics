/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Class AliRsnAnalysis
//             Reconstruction and analysis of K* Rsn
// ........................................
// ........................................
// ........................................
// ........................................
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef AliRsnAnalysis_H
#define AliRsnAnalysis_H

#include "AliRsnPID.h"

class TTree;
class TArrayI;
class TObjArray;
class AliRsnPair;
class AliRsnPID;

class AliRsnAnalysis : public TObject
{
public:

	AliRsnAnalysis(const char *branchname = "events");
	AliRsnAnalysis(const AliRsnAnalysis &copy) :
      TObject(copy),fSkipUnbalanced(kFALSE),fStep(1000),fBranchName("events"),
      fMixMultCut(10),fMixVzCut(0.5),fNEvents(0),fPID(0x0),
      fPairs(0x0),fMixPairs(0x0),fTree(0x0),fMatches(0x0),fEvents(0x0) { }
	AliRsnAnalysis& operator=(const AliRsnAnalysis & /*copy*/) { return (*this); }
	virtual ~AliRsnAnalysis() {Clear();}
	virtual void Clear(Option_t *option = "C");
	
	/* setters */
	void    SetPID(AliRsnPID *pid) {fPID = pid;}
	void    SetEventsTree(TTree *tree);
	void    SetStep(Int_t step) {fStep = step;}
    void    SetBranchName(const char *name) {fBranchName = name;}
    void    SetMixMultiplicityCut(Int_t cut) {fMixMultCut = cut;}
    void    SetMixVzCut(Double_t cut) {fMixVzCut = cut;}
    void    SetSkipUnbalanced(Bool_t doit = kTRUE) {fSkipUnbalanced = doit;}
	
	/* working routines */
	void    AddPair(AliRsnPair *pair);
	Stat_t  Process();
    Stat_t  EventMixing(Int_t nEventsToMix); 
	void    SaveOutput(const char *fileName, const char *fileOpt) const;
	
private:

	/* service methods (private) */
	AliRsnEvent *Evt(Int_t i);
    void         FindMatches(Int_t nEventsToMatch);
    Bool_t       AdjustPID(AliRsnEvent* &event);
    
    Bool_t        fSkipUnbalanced;  // reject events with no pos or no neg particles
	Int_t         fStep;            // step for progress message
    TString       fBranchName;      // name of event branch
    Int_t         fMixMultCut;      // multiplicity cut for event mixing
    Double_t      fMixVzCut;        // difference in VZ cut for event mixing
    Int_t         fNEvents;         // number of events (it cannot be accessed by user)
    
    AliRsnPID    *fPID;         //! if necessary, a new alternative PID object can be added
	TObjArray    *fPairs;       //! collection of particle pairs to be read for 1-event analysis
	TObjArray    *fMixPairs;    //! collection of particle pairs to be read for event-mixing
	TTree        *fTree;        //! TTree of events (can not be created here, must be passed)
    TArrayI      *fMatches;     //! collection of indexes of all events that match each one for mixing
    TClonesArray *fEvents;      //! collection of events --> branch of fTree
	
	// Rsn analysis implementation
	ClassDef(AliRsnAnalysis,1)
};

#endif
