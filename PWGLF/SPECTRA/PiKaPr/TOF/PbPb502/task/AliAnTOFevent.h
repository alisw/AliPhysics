#ifndef AliAnTOFevent_H
#define AliAnTOFevent_H

#include "AliAnTOFtrack.h"
#include "AliESDVertex.h"
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

  ///Masks
  Short_t fEvtMultBin;                        /// Binned multiplicity
  Double32_t fVtxX;                           //[-0.128,0.128,8]  Vtx X
  Double32_t fVtxY;                           //[0.154,0.410,8]  Vtx Y
  Double32_t fVtxZ;                           //[-12.8,12.8,8]  Vtx Z
  std::vector<AliAnTOFtrack> fAliAnTOFtracks; /// Array of AliAnTOFtrack

  //******************************
  ////////Utility methods/////////
  //******************************

  ///
  ///Method to adopt an event vertex
  void AdoptVertex(const AliESDVertex* vtx);

  ///
  ///Method to get the number of tracks
  Int_t GetNtracks() const { return fAliAnTOFtracks.size(); }

  ///
  ///Mehod to reset the event
  void Reset();

  ///
  ///Method to get a track from the event, if the index is not given (or negative) it returns the last element of the vector
  AliAnTOFtrack* GetTrack(Int_t i = -1);


  ///
  /// Status Printer
  void Print() const
  {
    Printf(" fEvtMultBin = %i", fEvtMultBin);
    Printf(" fVtxX = %f", fVtxX);
    Printf(" fVtxY = %f", fVtxY);
    Printf(" fVtxZ = %f", fVtxZ);
  }
};

#endif
