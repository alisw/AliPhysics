#ifndef ALIGENPILEUP_H
#define ALIGENPILEUP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          Class AliGenPileup
//   This is a generator of beam-beam pileup.
//   It generates interactions within 3 orbits (+-1) around
//   the trigger event. The trigger event itself is chosen
//   randomly among the bunch crossings within the central orbit.
//   The user can decide whenever to include in the simulation the
//   "trigger" interaction or not. This is handled by the
//   GenerateTrigInteraction(Bool_t flag) method.
//   In the case the trigger interaction is included, it is
//   generated using the same settings (vertex smear for example) as
//   the pileup events.
//   In case the trigger simulation is not included, the user can make
//   a cocktail of generator used to produce the trigger interaction and
//   AliGenPileup. In this case in order to avoid a fake increase of the rate around the
//   trigger, the number of background events within the bunch
//   crossing of the trigger is readuced by one.
//   The beam profile (the list of the active bunch crossings) can be
//   controlled via the SetBCMask(const char *mask) method. The syntax
//   follows the one in AliTriggerBCMask class. For example:
//   "3564H" would mean that all the bunch corssings within the orbit
//   are aloowed (which is of course unphysical). In case one wants to simulate
//   one-bunch-crossing-per-orbit scenario, the way to do it is to put something like:
//   "1H3563L" or similar.
//   The SetGenerator(AliGenerator *generator, Float_t rate) method is
//   used in order to define the generator to be used. The second argument is the pileup
//   rate in terms of #_of_interactions/bunch-crossing = sigma_tot * luminosity.
//   The pileup generation time window can be set via
//   AliGenerator::SetPileUpTimeWindow(Float_t pileUpTimeW) method. By the default the
//   window is set to 88micros (= TPC readout window).
//      
// cvetan.cheshkov@cern.ch  9/12/2008
//-------------------------------------------------------------------------

#include "AliGenCocktail.h"
#include "AliTriggerBCMask.h"
class TFormula;
class AliGenPileup : public AliGenCocktail
{
 public:
    AliGenPileup();
    virtual ~AliGenPileup();

    virtual void Generate();
    virtual void SetRandomise(Bool_t flag);
    virtual void UsePerEventRates();
	    
    void         SetGenerator(AliGenerator *generator, Float_t rate, Bool_t flag = kFALSE);
    //void         SetGenerator(AliGenerator *generator, Float_t rate);
    Bool_t       SetBCMask(const char *mask);
    void         GenerateTrigInteraction(Bool_t flag) {fGenTrig = flag;}

 protected:
    virtual void AddGenerator
      (AliGenerator *Generator, const char* Name, Float_t RateExp, TFormula* formula = 0 );

    AliTriggerBCMask fBCMask;    // Mask used to tag the active bunch-crossings within an orbit
    Bool_t           fGenTrig;   // Generate or not the trigger interaction
    Bool_t           fFlag;      // fixed interaction rate (integer)

 private:
    AliGenPileup(const AliGenPileup &gen);
    AliGenPileup & operator=(const AliGenPileup & gen);

    ClassDef(AliGenPileup,1) // Beam-beam pileup generator based on cocktail generator
};

#endif

