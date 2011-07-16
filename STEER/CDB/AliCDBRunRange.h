#ifndef ALI_CDB_RUN_RANGE_H
#define ALI_CDB_RUN_RANGE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBRunRange					   //
//  defines the run validity range of the object:		   //
//  [fFirstRun, fLastRun] 					   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliCDBRunRange: public TObject {

public:
	AliCDBRunRange();
	AliCDBRunRange(Int_t firstRun, Int_t lastRun);

	virtual ~AliCDBRunRange();

	
	Int_t 	GetFirstRun() const {return fFirstRun;};
	Int_t 	GetLastRun() const {return fLastRun;};

	void 	SetFirstRun(Int_t firstRun) {fFirstRun = firstRun;};
	void 	SetLastRun(Int_t lastRun) {fLastRun = lastRun;};

	void 	SetRunRange(Int_t firstRun, Int_t lastRun) 
	                {fFirstRun=firstRun; fLastRun = lastRun;};


	Bool_t IsValid() const;
	
	Bool_t IsAnyRange() const {return fFirstRun < 0 && fLastRun < 0;};

	Bool_t Overlaps(const AliCDBRunRange& other) const;

	Bool_t Comprises(const AliCDBRunRange& other) const;

	virtual Bool_t IsEqual(const TObject *obj) const;

	static Int_t Infinity() {return fgkInfinity;}
	
private:

	Int_t fFirstRun;	// first valid run
	Int_t fLastRun;		// last valid run	

	static const Int_t fgkInfinity = 999999999; 	 //! Flag for "infinity"

	ClassDef(AliCDBRunRange, 1);	
};


#endif
