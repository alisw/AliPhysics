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

// $Id: AliJTrack.cxx,v 1.2 2008/01/21 11:56:39 djkim Exp $

////////////////////////////////////////////////////
//
//  \file AliJTrack.cxx
//  \brief
//  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
//  \email: djkim@jyu.fi
//  \version $Revision: 1.1 $
//  \date $Date: 2008/05/02 11:56:39 $
//  
// class encapsulating aliroot track information
////////////////////////////////////////////////////

#include "AliPhJBaseTrack.h"
#include "AliJTrack.h"

//ClassImp(AliJTrack)

//______________________________________________________________________________
AliJTrack::AliJTrack() : 
  AliPhJBaseTrack(),
  fChi2perNDF(-999),
  fChi2Trig(-999),
  fRecFlags(-999),
  fTPCdEdx(-999),         
  fTPCnClust(-999),
  fImapactXY(-999),
  fImapactZ(-999),
  fTPCDCAXY(-999),
  fTPCDCAZ(-999),
  fTPCClustPerFindClust(-999),
  fTPCChi2PerClust(-999),
  fKinkIndex(-999),
  fstatus(-999)

{
  // default constructor
  SetPID((Double_t*)NULL);
  SetExternalDiaCovariance((Double_t*)NULL);
}

//______________________________________________________________________________
AliJTrack::AliJTrack(const AliJTrack& a):
  AliPhJBaseTrack(a),
  fChi2perNDF(a.fChi2perNDF),
  fChi2Trig(a.fChi2Trig),
  fRecFlags(a.fRecFlags),
  fTPCdEdx(a.fTPCdEdx),
  fTPCnClust(a.fTPCnClust),
  fImapactXY(a.fImapactXY),
  fImapactZ(a.fImapactZ),
  fTPCDCAXY(a.fTPCDCAXY),
  fTPCDCAZ(a.fTPCDCAZ),
  fTPCClustPerFindClust(a.fTPCClustPerFindClust),
  fTPCChi2PerClust(a.fTPCChi2PerClust),
  fKinkIndex(a.fKinkIndex),
  fstatus(a.fstatus)
{ 
  //copy constructor
  for(Int_t i=0;i<10;i++) ftrkPID[i]=a.ftrkPID[i];
  for(Int_t i=0;i<5;i++)  fextDiaCov[i]=a.fextDiaCov[i];
}
//______________________________________________________________________________
void AliJTrack::ConvertAliPID(){

  // Converts AliPID array.
  // The numbering scheme is the same for electrons, muons, pions, kaons, and protons.
  // Everything else has to be set to zero.

  ftrkPID[kDeuteronAli] = 0.;
  ftrkPID[kTritonAli]   = 0.;
  ftrkPID[kHelium3Ali]  = 0.;
  ftrkPID[kAlphaAli]    = 0.;
  ftrkPID[kUnknownAli]  = 0.;
  
  return;

}

//______________________________________________________________________________
AliJTrack&  AliJTrack::operator=(const AliJTrack& trk){
  //operator = 
  if(this != &trk){
    AliPhJBaseTrack::operator=(trk);
    for(Int_t i=0;i<10;i++){
      ftrkPID[i] = trk.ftrkPID[i];
    }
    fChi2perNDF = trk.fChi2perNDF;
    fChi2Trig   = trk.fChi2Trig;
    fRecFlags   = trk.fRecFlags;
    fTPCdEdx    = trk.fTPCdEdx;
    fTPCnClust  = trk.fTPCnClust;
    fImapactXY = trk.fImapactXY;
    fImapactZ  = trk.fImapactZ;
    fTPCDCAXY  = trk.fTPCDCAXY;
    fTPCDCAZ   = trk.fTPCDCAZ;
    fTPCClustPerFindClust = trk.fTPCClustPerFindClust;
    fTPCChi2PerClust      = trk.fTPCChi2PerClust;
    fKinkIndex  = trk.fKinkIndex;
    fstatus     = trk.fstatus;
    for(Int_t i=0; i<5; i++){
      fextDiaCov[i] = trk.fextDiaCov[i];
    }
  }

  return *this;
}


//___________________________________________________________________
void AliJTrack::SetPID(const Double_t *pid) {
  //set pid
  if(pid) for(Int_t i=0; i<10; ++i) ftrkPID[i]=pid[i];
  else {
    for(Int_t i=0; i<10; ftrkPID[i++]=0.){;}
    ftrkPID[kUnknownAli]=1.;
  }
}

//___________________________________________________________________
void AliJTrack::SetExternalDiaCovariance(const Double_t *ecov) {
  //set diagonal elements of the covariance  matrix 
  if(ecov) {
    for(Int_t i=0; i<5; i++) fextDiaCov[i]=ecov[i];
  } else {
    for(Int_t i=0; i<5;i++) fextDiaCov[i]=-999;
  }
}
