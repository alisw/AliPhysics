/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Class AliRsnManager
//             Reconstruction and analysis of K* Rsn
// ........................................
// ........................................
// ........................................
// ........................................
// 
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef AliRsnManager_H
#define AliRsnManager_H

#include "AliRsnPID.h"

class TList;
class TObjArray;
class AliRsnPair;

class AliRsnManager : public TObject
{
public:

	AliRsnManager();
	virtual ~AliRsnManager() {Clear();}
	virtual void Clear(Option_t *option = "");
	
	/* setters */
	void    SetEventsTree(TTree *tree);
	void    SetStep(Int_t step) {fStep = step;}
    void    SetMixEvents(Int_t num) {fMixEvents = num;}
    void    SetMixMultiplicityCut(Int_t cut) {fMixMultCut = cut;}
    void    SetMixVzCut(Double_t cut) {fMixVzCut = cut;}
    void    SetQueueSize(Int_t size);
    void    SetUsePID(Bool_t doit) {fUsePID = doit;}
	
	/* working routines */
	void    AddPair(AliRsnPair *pair);
	Stat_t  Process(AliRsnEvent *event);
    Stat_t  Mix(AliRsnEvent *event); 
	Bool_t  CanBeMixed(AliRsnEvent *ev1, AliRsnEvent *ev2);
	TList*  GetOutputList() {return fOutputList;}
	
private:
    
    AliRsnManager(const AliRsnManager &copy) :
      TObject(copy),fUsePID(kTRUE),fStep(100),
      fMixEvents(10),fMixMultCut(10),fMixVzCut(0.5),fQueuePos(-1),
      fPairs(0x0),fMixPairs(0x0),fBuffer(0x0),fOutputList(0x0) { }
	AliRsnManager& operator=(const AliRsnManager & /*copy*/) { return (*this); }
    
    Bool_t        fUsePID;      // flag to switch between PID/noPID analysis
	Int_t         fStep;        // step for progress message
    Int_t         fMixEvents;   // number of events to be mixed (maximum)
    Int_t         fMixMultCut;  // multiplicity cut for event mixing
    Double_t      fMixVzCut;    // difference in VZ cut for event mixing
    Int_t         fQueuePos;    // position in queue for event adding
    
	TObjArray         *fPairs;       //! collection of particle pairs to be read for 1-event analysis
	TObjArray         *fMixPairs;    //! collection of particle pairs to be read for event-mixing
    TObjArray         *fBuffer;      //! buffer for events to mix
    TList             *fOutputList;  //! list of output histograms
	
	// Rsn analysis implementation
	ClassDef(AliRsnManager,1)
};

#endif
