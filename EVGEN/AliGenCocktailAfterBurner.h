#ifndef AliGenCocktailAfterBurner_H
#define AliGenCocktailAfterBurner_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Container class for AliGenerator through recursion.
// (Container is itself an AliGenerator)
// Author: andreas.morsch@cern.ch 
//
#include "AliGenCocktailAfterBurner.h"
#include "AliGenCocktail.h"


class AliGenCocktailEntry;
class AliStack;

class AliGenCocktailAfterBurner : public  AliGenCocktail
{
 public:
    AliGenCocktailAfterBurner();
//    AliGenCocktailAfterBurner(const AliGenCocktailAfterBurner &cocktail){}
     
    virtual ~AliGenCocktailAfterBurner();
    virtual void Init();
    virtual void Generate();
    virtual void SetTracks(Int_t stackno);
    //
    // Add a new generator to the list
    virtual void AddAfterBurner
	(AliGenerator *Generator, char* Name, Float_t RateExp );
    AliGenCocktailAfterBurner & operator=(const AliGenCocktailAfterBurner & rhs);
    
    AliStack* GetStack(Int_t n);
    AliStack* GetActiveStack() {return fActiveStack;}
    
    AliGenerator* GetCurrentGenerator();
    virtual void  SetActiveEventNumber(Int_t actev);
    Int_t GetActiveEventNumber() {return fActiveEvent;}
    virtual void SetNumberOfEvents(Int_t n)   {fNumberOfEvents=n;}
    virtual Int_t GetNumberOfEvents() {return fNumberOfEvents;}

    static AliMCProcess IntToMCProcess(Int_t no);
 protected:
    Int_t fNAfterBurners;       // Number of afterburners  
    TList  *fAfterBurnerEntries;// List of afterburners
    Bool_t fGenerationDone;
    TObjArray *fInternalStacks; //! List of internal stacks
    Int_t fCurrentEvent;        //  Number of current event/stack
    
    Int_t fNumberOfEvents;      // Number of events to process

    AliStack* fActiveStack;   //! pointer to the current stack
    Int_t fActiveEvent;       //HBT Processor needs more then one event to do correlations
                              //Due to complications in fortran, it first calls C routine
                              //that sets the active event to be read. All alihbtp_gettrack
                              //are addressed to this event
    
    AliGenerator *fCurrentGenerator;      // Current event generator 
    ClassDef(AliGenCocktailAfterBurner,1) // Particle cocktail generator a la SHAKER
};

inline  AliGenerator*  
    AliGenCocktailAfterBurner::GetCurrentGenerator()
{
  return fCurrentGenerator;
}


#endif





