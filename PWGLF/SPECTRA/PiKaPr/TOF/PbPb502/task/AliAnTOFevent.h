#ifndef AliAnTOFevent_H
#define AliAnTOFevent_H

#include "AliAnTOFtrack.h"
#include "Rtypes.h"
#include <vector>

using namespace AliUtilTOFParams;

class AliESDEvent;
class AliESDtrack;

///////////////////////////////////////////////////////////////////////////////
///                                                                          //
///               Class container for TOF analysis of Pi/K/p.                //
///                                                                          //
///                                                                          //
/// Authors:                                                                 //
/// N. Jacazio,  nicolo.jacazio[AROBASe]bo.infn.it                           //
///////////////////////////////////////////////////////////////////////////////

class AliAnTOFevent {
  public:
  //Constructors and destructor
  AliAnTOFevent();
  virtual ~AliAnTOFevent();

  //Masks
  Short_t fEvtMultBin; //Binned multiplicity
  std::vector<AliAnTOFtrack> fAliAnTOFtracks;

  //******************************
  ////////Utility methods/////////
  //******************************

  ///
  ///Method to get the number of tracks
  Int_t GetNtracks();

  ///
  ///Mehod to reset the event
  void Reset();

  ///
  ///Method to get a track from the event, if the index is not given (or negative) it returns the last element of the vector
  AliAnTOFtrack* GetTrack(Int_t i = -1);

  ClassDef(AliAnTOFevent, 1); //AliAnTOFevent : TOF analysis event container class
};

#endif
