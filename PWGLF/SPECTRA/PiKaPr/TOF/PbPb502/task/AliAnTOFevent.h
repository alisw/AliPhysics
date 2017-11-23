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
  UChar_t fVtxX;                              /// Binned Vtx X
  UChar_t fVtxY;                              /// Binned Vtx Y
  UChar_t fVtxZ;                              /// Binned Vtx Z
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
  /// Setter for the Vtx Z position
  void SetVtxZ(const Double_t z)
  {
    fVtxZ = static_cast<UChar_t>(BinData(z, -12.8, 12.8, 256));
  }

  ///
  /// Setter for the Vtx X position
  void SetVtxX(const Double_t x)
  {
    fVtxX = static_cast<UChar_t>(BinData(x, -0.128, 0.128, 256));
  }

  ///
  /// Setter for the Vtx Y position
  void SetVtxY(const Double_t y)
  {
    fVtxY = static_cast<UChar_t>(BinData(y, -0.128, 0.128, 256));
  }

  ///
  /// Getter for the Vtx Z position
  Double_t GetVtxZ() const { return GetBinnedData(fVtxZ, -12.8, 12.8, 256); }

  ///
  /// Getter for the Vtx X position
  Double_t GetVtxX() const { return GetBinnedData(fVtxX, -0.128, 0.128, 256); }

  ///
  /// Getter for the Vtx Y position
  Double_t GetVtxY() const { return GetBinnedData(fVtxY, -0.128, 0.128, 256); }
};

#endif
