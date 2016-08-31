#ifndef ALIMCGENHANDLER_H
#define ALIMCGENHANDLER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id: AliMCGenHandler.h 61697 2013-03-26 12:17:26Z agheata $ */

//-------------------------------------------------------------------------
//                          Class AliMCGenHandler
// This class can be used with the analysis framework to generate event on
// the fly and analyse them.
//      
// Origin: Andrei Gheata, Jan Fiete Grosse-Oetringhaus
//-------------------------------------------------------------------------

#include "AliInputEventHandler.h"

class TFile;
class TTree;
class AliMCEvent;
class AliGenerator;
class AliRunLoader;
class AliStack;
class AliHeader;
class TMacro;

class AliMCGenHandler : public AliInputEventHandler
{
public:

    AliMCGenHandler();
    AliMCGenHandler(const char* name, const char* title);
    virtual ~AliMCGenHandler();

    virtual Bool_t       Init(Option_t* /*opt*/);
    virtual Bool_t       Init(TTree* tree, Option_t* opt) { return AliInputEventHandler::Init(tree, opt); }
    virtual Bool_t       BeginEvent(Long64_t entry);
    virtual Bool_t       FinishEvent();

    void                 CreateGenerator();

    AliMCEvent* MCEvent() const {return fMCEvent;}
    const AliGenerator* GetGenerator() const { return fGenerator; }

    void 		 SetGenerator(AliGenerator* generator) { fGenerator = generator; }
    
    void		 SetSeedMode(Int_t mode) { fSeedMode = mode; }
    void		 SetSeed(UInt_t seed) { fSeed = seed; }
    UInt_t		 GetSeed() { return fSeed; }
    
    void		 SetGeneratorMacroPath(const char* macroPath) { fGeneratorMacroPath = macroPath; }
    void		 SetGeneratorMacroParameters(const char* params) { fGeneratorMacroParameters = params; }
    void		 SetGeneratorCustomization(TMacro* macro) { fGeneratorCustomization = macro; }

private:
    AliMCGenHandler(const AliMCGenHandler& handler);             
    AliMCGenHandler& operator=(const AliMCGenHandler& handler);  

    AliMCEvent            *fMCEvent;            //! MC Event
    Long64_t		   fEventNumber;	//! current event number
    AliStack* 		   fStack;		//! current AliStack pointer
    AliHeader*		   fHeader;		//! current AliHeader pointer

    AliGenerator	  *fGenerator;		// generator
    Int_t		   fSeedMode;		// which seed is to be used: 0 (default): nothing/set externally; 1: use fSeed; 2: current time; 3: AliEn job id
    UInt_t		   fSeed;		// can be used to set seed manually (fSeedMode == 1); contains last used seed (fSeedMode == 2 || 3)
    
    TString		   fGeneratorMacroPath; // path to macro creating the generator object
    TString		   fGeneratorMacroParameters; // parameters passed to the creating macro
    TMacro		  *fGeneratorCustomization; // customization macro for generator object

    ClassDef(AliMCGenHandler, 2)  // MC Gen Handler
};
#endif 
