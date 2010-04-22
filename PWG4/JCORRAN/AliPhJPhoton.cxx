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

// $Id: AliPhJPhoton.cxx,v 1.4 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
//
//  \file AliPhJPhoton.cxx
//  \brief
//  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
//  \email: djkim@jyu.fi
//  \version $Revision: 1.4 $
//  \date $Date: 2008/05/08 13:44:45 $
//
//  Class encapsulation photon information from aliroot
////////////////////////////////////////////////////

#include "AliPhJBaseTrack.h"
#include "AliPhJPhoton.h"

ClassImp(AliPhJPhoton)

//______________________________________________________________________________
AliPhJPhoton::AliPhJPhoton() : 
  AliPhJBaseTrack(),
  fChi2(-999),              
  fTof(-999),                   
  fX(-999),               
  fY(-999),                  
  fZ(-999),
  fProbPhot(-999)
{
  // default constructor
}

//_____________________________________________________________________________
AliPhJPhoton::AliPhJPhoton(const AliPhJPhoton& a) : 
  AliPhJBaseTrack(a),
  fChi2(a.fChi2),
  fTof(a.fTof),
  fX(a.fX),
  fY(a.fY),
  fZ(a.fZ),
  fProbPhot(a.fProbPhot)
{
  //copy constructor
}


//_____________________________________________________________________________
AliPhJPhoton& AliPhJPhoton::operator=(const AliPhJPhoton& photon){
  //operator=
  if(this != &photon){
    AliPhJBaseTrack::operator=(photon); 
    fChi2 = photon.fChi2;
    fTof  = photon.fTof;
    fX    = photon.fX;
    fY    = photon.fY;
    fZ    = photon.fZ;
    fProbPhot = photon.fProbPhot;
  }

  return *this;
}
