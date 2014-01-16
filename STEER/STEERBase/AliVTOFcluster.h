// -*- mode: C++ -*- 
#ifndef ALIVTOFCLUSTER_H
#define ALIVTOFCLUSTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

//----------------------------------------------------------------------//
//                                                                      //
// AliVTOFcluster Class                                                //
//                                                                      //
//----------------------------------------------------------------------//

#include <TObject.h>
#include <TArrayI.h>
#include <TArrayF.h>

class AliVTOFcluster : public TObject {

 public:

  AliVTOFcluster() { }
  virtual ~AliVTOFcluster() { }
  AliVTOFcluster(const AliVTOFcluster & source);
  AliVTOFcluster & operator=(const AliVTOFcluster & source);

  virtual Int_t GetClusterIndex() const {return 0;} // cluster index
  virtual Int_t GetTOFchannel() const {return 0;}; // TOF channel
  virtual Float_t GetTime() const {return 0;}; // TOF time
  virtual Float_t GetTimeRaw() const {return 0;}; // TOF raw time
  virtual Float_t GetTOT() const {return 0;}; // TOF ToT
  virtual Int_t GetLabel(Int_t ) const {return 0;};
  virtual Int_t GetDeltaBC() const {return 0;};
  virtual Int_t GetL0L1Latency() const {return 0;};
  virtual Bool_t GetStatus() const {return 0;};
  virtual Float_t GetZ() const {return 0;};
  virtual Float_t GetPhi() const {return 0;};
  virtual Float_t GetR() const {return 0;};
  virtual Int_t GetNMatchableTracks() const {return 0;};
  virtual Int_t GetTrackIndex(Int_t ) const {return 0;};
  virtual Float_t GetDistanceInStripPlane(Int_t )   const {return 0;};
  virtual Float_t GetDx(Int_t )  const {return 0;};
  virtual Float_t GetDy(Int_t )  const {return 0;};
  virtual Float_t GetDz(Int_t )  const {return 0;};
  virtual Float_t GetLength(Int_t ) const {return 0;};
  virtual Double_t GetIntegratedTime(Int_t ,Int_t ) const {return 0;};

  ClassDef(AliVTOFcluster, 1) // TOF matchable cluster

}; 

#endif
