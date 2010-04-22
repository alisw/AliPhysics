/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notifce   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id: AliJMCTrack.cxx,v 1.2 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
//
//  \file AliJMCTrack.cxx
//  \brief
//  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
//  \email: djkim@jyu.fi
//  \version $Revision: 1.2 $
//  \date $Date: 2008/05/08 13:44:45 $
//
//  class which encapsulates track  monte carlo information
////////////////////////////////////////////////////

#include "AliJMCTrack.h"

ClassImp(AliJMCTrack)

//______________________________________________________________________________
AliJMCTrack::AliJMCTrack() : 
  AliPhJBaseTrack(),
  fPdgCode(0),
  fStatus(0),
  fFlag(0),
  fLabel(0),
  fVx(0),
  fVy(0),
  fVz(0),
  fEta(0),
  fPrimary(0),
  fPHOS(0),
  fEMCAL(0),
  fTPC(0),
  fPtHard(0)
{
  // default constructor
  fMother[0] =   fMother[1] = 0;
  fDaughter[0] =   fDaughter[1] = 0;


}

//______________________________________________________________________________
AliJMCTrack::AliJMCTrack(const AliJMCTrack& a):
  AliPhJBaseTrack(a),
  fPdgCode(a.fPdgCode),
  fStatus(a.fStatus),
  fFlag(a.fFlag),
  fLabel(a.fLabel),
  fVx(a.fVx),
  fVy(a.fVy),
  fVz(a.fVz),
  fEta(a.fEta),
  fPrimary(a.fPrimary),
  fPHOS(a.fPHOS),
  fEMCAL(a.fEMCAL),
  fTPC(a.fTPC),
  fPtHard(a.fPtHard)
{
//copy constructor
  for(Int_t i=0;i<2;i++){
    fMother[i]   = a.fMother[i];
    fDaughter[i] = a.fDaughter[i];
  }
}


//______________________________________________________________________________

AliJMCTrack& AliJMCTrack::operator=(const AliJMCTrack& trk){
  //operator= 
  if(this != &trk){
    AliPhJBaseTrack::operator=(trk);
    fPdgCode   =  trk.fPdgCode;
    fStatus    =  trk.fStatus;
    fFlag      =  trk.fFlag;
    fLabel     =  trk.fLabel;
    for(Int_t i=0;i<2;i++){
      fMother[i]   = trk.fMother[i];
      fDaughter[i] = trk.fDaughter[i];
    }
    fVx  = trk.fVx;
    fVy  = trk.fVy;
    fVz  = trk.fVz;
    fEta = trk.fEta;
    fPrimary = trk.fPrimary;
    fPHOS    = trk.fPHOS;
    fEMCAL   = trk.fEMCAL;
    fTPC     = trk.fTPC;
    fPtHard  = trk.fPtHard;
  }
  return *this;
}





