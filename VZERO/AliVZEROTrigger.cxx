#include "AliRun.h"
#include "AliRunLoader.h"

#include "AliVZEROTrigger.h"

//______________________________________________________________________
ClassImp(AliVZEROTrigger)
////////////////////////////////////////////////////////////////////////
//
// Version 1
//
// AliVZEROTrigger: 
//
////////////////////////////////////////////////////////////////////////

//______________________________________________________________________
AliVZEROTrigger::AliVZEROTrigger()
  :AliTriggerDetector(),
   fAdcThresHold(0.0),
   fTimeWindowWidth(50.0)
   
{
   SetName("VZERO");
   CreateInputs();

   SetAdcThreshold();
   SetTimeWindowWidth();
}

//______________________________________________________________________
void AliVZEROTrigger::CreateInputs()
{
   // inputs

   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;

   fInputs.AddLast( new AliTriggerInput( "VZERO_LEFT", "At least one hit in the left VZERO", 0x01 ) );
   fInputs.AddLast( new AliTriggerInput( "VZERO_RIGHT","At least one hit in the right VZERO", 0x02 ) );
   fInputs.AddLast( new AliTriggerInput( "VZERO_AND",  "At least one hit in the each (left and right) VZERO", 0x04 ) );
   fInputs.AddLast( new AliTriggerInput( "VZERO_OR",   "At least one hit in the one (left one right) VZERO", 0x08 ) );

   fInputs.AddLast( new AliTriggerInput( "VZERO_BEAMGAS", "Beam gas VZERO trigger ", 0x010 ) );
}

//______________________________________________________________________
void AliVZEROTrigger::Trigger()
{
  
  //  ********** Get run loader for the current event **********
  AliRunLoader* runLoader = gAlice->GetRunLoader();

  AliVZEROLoader* loader = (AliVZEROLoader* )runLoader->GetLoader( "VZEROLoader" );

  loader->LoadDigits("READ");
  TTree* vzeroDigitsTree = loader->TreeD();
  if (!vzeroDigitsTree) return;

  TClonesArray* vzeroDigits = new TClonesArray("AliVZEROdigit",1000);
  TBranch* digitBranch = vzeroDigitsTree->GetBranch("VZERODigit");
  digitBranch->SetAddress(&vzeroDigits);

  // number of hits in left/right
  Int_t nLeftDig  = 0;
  Int_t nRightDig = 0;
  
  // first time 
  Float_t firstTimeLeft  = 9999.0;
  Float_t firstTimeRight = 9999.0;
  Float_t TimeHalfWidth  = fTimeWindowWidth/2.0;
 
  // loop over vzero entries
  Int_t nEntries = (Int_t)vzeroDigitsTree->GetEntries();
  for (Int_t e=0; e<nEntries; e++) {
    vzeroDigitsTree->GetEvent(e);

    Int_t nDigits = vzeroDigits->GetEntriesFast();
    
    for (Int_t d=0; d<nDigits; d++) {
      //      vzeroDigitsTree->GetEvent(d);
      AliVZEROdigit* digit = (AliVZEROdigit*)vzeroDigits->At(d);
      
      Int_t   PMNumber   = digit->PMNumber();
      Float_t adc        = digit->ADC();
      Float_t tdc        = digit->Time(); // in 100 of picoseconds
      
      if (PMNumber<=31 && adc>fAdcThresHold) {
	if (tdc>(30.0-TimeHalfWidth) && tdc<(30.0+TimeHalfWidth)) nRightDig++;
	if (tdc<firstTimeRight) firstTimeRight = tdc;
      }      
      if (PMNumber>=32 && adc>fAdcThresHold) {
 	if (tdc>(114.0-TimeHalfWidth) && tdc<(114.0+TimeHalfWidth)) nLeftDig++;
	if (tdc<firstTimeLeft) firstTimeLeft = tdc;
      }	
      
    } // end of loop over digits
  } // end of loop over events in digits tree
  
  // Beam gas trigger set from the time difference. The time it takes
  // to travel between the two counters is ~14.3 ns = 143 * 100 ps.
  //  NB: this should be defined
  // from time windows relative to the time of the bunch crossing!
  // beam gas comming from the left ...

  if (TMath::Abs(TMath::Abs(firstTimeLeft - firstTimeRight)-143) < 20) // time window of 2 ns
    SetInput( "VZERO_BEAMGAS" );

  if (nLeftDig > 0)
      SetInput( "VZERO_LEFT" );

  if (nRightDig > 0)
      SetInput( "VZERO_RIGHT" );
  
  if (nLeftDig>0 || nRightDig>0) {
      SetInput( "VZERO_OR" );

    if (nLeftDig>0 && nRightDig>0) {
        SetInput( "VZERO_AND" );   
    }
  }
  
  AliDebug(1,Form("VZERO PMs fired: %d (left) %d (right)", nLeftDig, nRightDig));
 
  return;
}


