/*
***********************************************************
  Implementation of reduced ESD information classes for
  quick analysis.
  Contact: i.c.arsene@gsi.de, i.c.arsene@cern.ch
  2012/06/21
  *********************************************************
*/

#ifndef ALIEVENTPLANECUTS_H
#include "AliEventPlaneCuts.h"
#endif

#include <TMath.h>
#include <TClonesArray.h>
#include <TClass.h>

ClassImp(AliEventPlaneCuts)



//_______________________________________________________________________________
AliEventPlaneCuts::AliEventPlaneCuts() :
  TObject(),
  fCuts(),
  fNcuts(0),
  fName("")
{   
  //
  // Constructor
  //

  for(Int_t i=0; i<100; ++i)
    for(Int_t j=0; j<4; ++j){
      fCuts[i][j]=0.0;
    }
}




//_______________________________________________________________________________
AliEventPlaneCuts::~AliEventPlaneCuts()
{
  //
  // De-Constructor
  //
}


//_______________________________________________________________________________
void AliEventPlaneCuts::CopyCuts(AliEventPlaneCuts* cuts){
  for(Int_t i=0; i<100; ++i){
    fCuts[i][0] = cuts->Type(i);
    fCuts[i][1] = cuts->Min(i);
    fCuts[i][2] = cuts->Max(i);
    fCuts[i][3] = cuts->ExcludeRange(i);
  }
  fNcuts = cuts->Ncuts();
}


