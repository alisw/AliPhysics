#ifndef AliGenCocktailAfterBurner_H
#define AliGenCocktailAfterBurner_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Container class for AliGenerator through recursion.
// (Container is itself an AliGenerator)
// Author: piotr.skowronski@cern.ch 
//

#include <TMCProcess.h>

#include "AliGenCocktail.h"
#include "AliRun.h"

class AliGenCocktailEntry;
class AliStack;

// ANDREAS MORSCH ------------------------------------- (
class AliCollisionGeometry;
// ANDREAS MORSCH ------------------------------------- )

class AliGenCocktailAfterBurner : public  AliGenCocktail
{
//container class for other generators
//extends AliGenCocktail functionality
//with possiblity of adding after-burners

 public:
    AliGenCocktailAfterBurner();
    virtual ~AliGenCocktailAfterBurner();
    
    virtual void  Init();
    virtual void  Generate();
    virtual void  SetTracks(Int_t stackno);
    //
    // Add a new generator to the list
    virtual void  AddAfterBurner
	(AliGenerator *Generator, char* Name, Float_t RateExp );
    
    AliStack*     GetStack(Int_t n) const;
    AliStack*     GetActiveStack() const{return fActiveStack;}

// ANDREAS MORSCH ------------------------------------- (
    AliCollisionGeometry* GetCollisionGeometry(Int_t n) const;
// ANDREAS MORSCH ------------------------------------- )
 
    AliGenerator* GetCurrentGenerator() const;
    virtual void  SetActiveEventNumber(Int_t actev);
    Int_t         GetActiveEventNumber() const {return fActiveEvent;}
    virtual Int_t GetNumberOfEvents() const {return gAlice->GetEventsPerRun() + fNBgEvents;}
    void          SetNBgEvents(Int_t nbg=0){fNBgEvents = nbg;}

    static TMCProcess IntToMCProcess(Int_t no);

 protected:
    Int_t fNAfterBurners;       // Number of afterburners  
    TList  *fAfterBurnerEntries;// List of afterburners
    Bool_t fGenerationDone;     // flag if generation is already done 
                                //   during first call of Generate method
                                //   if true just return event to gAlice
                                //   
    TObjArray *fInternalStacks; //! List of internal stacks

// ANDREAS MORSCH ------------------------------------- (
    AliCollisionGeometry** fCollisionGeometries; //! List of Collision Geometries
// ANDREAS MORSCH ------------------------------------- )
    
    Int_t fCurrentEvent;        //  Number of current event/stack
    

    AliStack* fActiveStack;   //! pointer to the current stack
    Int_t fActiveEvent;       //HBT Processor needs more then one event to do correlations
                              //Due to complications in fortran, it first calls C routine
                              //that sets the active event to be read. All alihbtp_gettrack
                              //are addressed to this event
    
    AliGenerator *fCurrentGenerator;      // Current event generator 
    Int_t fNBgEvents;                     //Nuber of backgrouns events 
                                          //(events that are generated only temporarly)
                                          //needed by some afterburners that works better with higher statistics 
                                          //this generates such a artificial one
 private:
    AliGenCocktailAfterBurner(const AliGenCocktailAfterBurner& in);
    AliGenCocktailAfterBurner & operator=(const AliGenCocktailAfterBurner & rhs);

    ClassDef(AliGenCocktailAfterBurner,2) // Particle cocktail generator a la SHAKER
                                          //background events added
};

inline  AliGenerator*  
    AliGenCocktailAfterBurner::GetCurrentGenerator() const
{
  return fCurrentGenerator;
}


#endif





