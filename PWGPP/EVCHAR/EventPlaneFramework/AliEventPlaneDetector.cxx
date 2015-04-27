
/*
***********************************************************
    Detector class that contains detector information (angles and weights)
    Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
***********************************************************
*/

#ifndef ALIEVENTPLANEDETECTOR_H
#include "AliEventPlaneDetector.h"
#endif


ClassImp(AliEventPlaneDetector)


//_______________________________________________________________________________
AliEventPlaneDetector::AliEventPlaneDetector() :
  TObject(),
  fPhi(0.0),
  fX(0.0),
  fY(0.0),
  fWeight(0.0),
  fEqualizedWeight(),
  fId(0),
  fEventPlaneDetectorMask(0),
  fBin(0)
{   
  //
  // Constructor
  //
  for(Int_t idet=0; idet<fgkEPMaxDetectors; ++idet)
     for(Int_t imeth=0; imeth<2; ++imeth) {
    fEqualizedWeight[idet][imeth]=0.0;
  }

}



//_______________________________________________________________________________
AliEventPlaneDetector::~AliEventPlaneDetector()
{
  //
  // De-Constructor
  //
}


