// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTracklet_H
#define AliEveTracklet_H

#include "AliEveTrack.h"

class AliVVertex;

//______________________________________________________________________________
// AliEveTracklet is a representation of SPD tracklet.
// It inherits from AliEveTrack to allow for common functionality
// regarding track counting.

class AliEveTracklet : public AliEveTrack
{
public:
  AliEveTracklet(Int_t index, const AliVVertex* pv, Float_t theta, Float_t phi, TEveTrackPropagator* prop=0);
  virtual ~AliEveTracklet() {}

  virtual void MakeTrack(Bool_t recurse=kTRUE);

  virtual void SecSelected(TEveTrack*);              // *SIGNAL*
  virtual void SecSelectedTracklet(AliEveTracklet*); // *SIGNAL*

  // ----------------------------------------------------------------

  static Float_t GetDefaultRadius();
  static void    SetDefaultRadius(Float_t r);

protected:
  static Float_t fgDefaultRadius;

private:
  AliEveTracklet(const AliEveTracklet&);            // Not implemented
  AliEveTracklet& operator=(const AliEveTracklet&); // Not implemented

  ClassDef(AliEveTracklet, 0); // Short description.
};

#endif
