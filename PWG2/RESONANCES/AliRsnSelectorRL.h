/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//  Class AliRsnSelectorRL
// ------------------------
// Reader for conversion of ESD output into the internal format
// used for resonance study.
// ---
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
// ---
// adapted for TSelector compliance
// by    : R. Vernet                 (email: renaud.vernet@cern.ch)
//-------------------------------------------------------------------------

#ifndef ALIRSNSELECTORRL_H
#define ALIRSNSELECTORRL_H

#include "AliPID.h"
#include "AliSelectorRL.h"
#include "AliRsnReader.h"

class TH3D;
class TH1D;
class TOrdCollection;
class TTree;
class TBranch;
class AliRunLoader;

class AliRsnSelectorRL : public AliSelectorRL, public AliRsnReader
{
public:

    AliRsnSelectorRL(TTree *tree = 0);
    AliRsnSelectorRL(const AliRsnSelectorRL& obj);
	virtual ~AliRsnSelectorRL();
    AliRsnSelectorRL& operator=(const AliRsnSelectorRL& obj);
	
	// TSelector-inherited member functions
	virtual Int_t   Version() const {return 1;}
	virtual void    Begin(TTree *tree) const;
	virtual void    SlaveBegin(TTree *tree);
	virtual void    Init(TTree *tree);
	virtual Bool_t  Process(Long64_t entry);
	virtual void    SetOption(const char *option) {fOption = option;}
	virtual void    SetObject(TObject *obj) {fObject = obj;}
	virtual void    SetInputList(TList *input) {fInput = input;}
	virtual TList  *GetOutputList() const {return fOutput;}
	virtual void    SlaveTerminate() ;
	virtual void    Terminate();
	
	// Other
	void            Clear(const Option_t *option = "");
	
protected:
 
 	// Parameters/flags
	TString*      fOutputPath;       //! path where output tree will be stored
	
	// Workaround for AliSelectorRL:
	Bool_t        fIsRunLoaderOpen;  //  flag to check if run loader is open
	
	// IO tree
	TTree*        fRsnEventTree;     //  output tree which should contain AliRsnEvents
	AliRsnEvent*  fRsnEvent;         //  pointer on AliRsnEvent to store
	TBranch*      fRsnEventBranch;   //  tree branch to store AliRsnEvents in
	
	// Rsn event reader implementation
	ClassDef(AliRsnSelectorRL,1)
};

#endif

