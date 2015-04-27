/*
***********************************************************
    Q-vector class for event plane correction framework and analysis
    Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
***********************************************************
*/


//#ifndef ALIEVENTPLANEQVECTOR_H
#include "AliEventPlaneQvector.h"
//#endif

#include <TMath.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TArrayS.h>
#include <iostream>

ClassImp(AliEventPlaneQvector)



//_______________________________________________________________________________
AliEventPlaneQvector::AliEventPlaneQvector() :
  fQvector(),
  fBin(-1),
  fMultiplicity(0.0),
  fEventPlaneStatus()
{
  //
  // Constructor
  //
  for(Int_t ih=0; ih<fgkEPMaxHarmonics; ++ih){
      fEventPlaneStatus[ih] = 0;
      for(Int_t ic=0; ic<2; ++ic)
        fQvector[ih][ic] = 0.0;
  }
        

}

//_______________________________________________________________________________
AliEventPlaneQvector::~AliEventPlaneQvector()
{
  //
  // De-Constructor
  //
}


