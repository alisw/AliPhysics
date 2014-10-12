/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Comment describing what this class does needed!

#include "AliJJet.h"

AliJJet::AliJJet(): 
    AliJBaseTrack(), 
    fArea(0),
    fConstituents()
{;}


AliJJet::AliJJet(float px,float py, float pz, float e, Int_t id, Short_t ptype, Char_t charge):
    AliJBaseTrack( px, py, pz, e, id, ptype, charge ),
    fArea(0),
    fConstituents()
{;}

AliJJet::AliJJet(const AliJJet& a):
    AliJBaseTrack( a ),
    fArea( a.fArea ),
    fConstituents( a.fConstituents )
{;}

AliJJet::AliJJet(const TLorentzVector & a):
    AliJBaseTrack( a ),
    fArea(0),
    fConstituents()
{;}


AliJJet& AliJJet::operator=(const AliJJet& trk){
  //operator =  
  if(this != &trk){
    AliJBaseTrack::operator=(trk);
    fArea = trk.fArea;
    fConstituents = trk.fConstituents;
  }
  return *this;
}

ClassImp(AliJJet)
