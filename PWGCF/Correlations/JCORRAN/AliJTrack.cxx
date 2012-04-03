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

#include "AliJBaseTrack.h"
#include "AliJTrack.h"

//ClassImp(AliJTrack)

//______________________________________________________________________________
AliJTrack::AliJTrack() : 
    AliJBaseTrack(),
    fFilterMap(0),
    fTPCnClust(-1),
    fTPCdEdx(-1), 
    fTOFbeta(0)
{
  // default constructor
  for( int i=0;i<kNAliJTrkPID;i++ ) SetPID( AliJTrkPID(i), 0, kTOF);
  for( int i=0;i<kNAliJTrkPID;i++ ) SetPID( AliJTrkPID(i), 0, kTPC);
  for( int i=0;i<kNAliJTrkPID;i++ ) SetPID( AliJTrkPID(i), 0, kTPCTOF);
  for( int i=0;i<kNAliJTrkPID;i++) fExpTOFbeta[i]= 0;
  for( int i=0;i<kNAliJTrkPID;i++) fExpTPCdEdx[i]= 0;
  for( int i=0;i<kNAliJTrkPID;i++) fTPCsigma[i]= 0;
  for( int i=0;i<kNAliJTrkPID;i++) fTOFsigma[i]= 0;
  for( int i=0;i<3;i++) fTPCTrack[i] = 0;

}

//______________________________________________________________________________
AliJTrack::AliJTrack(const AliJTrack& a):
    AliJBaseTrack(a),
    fFilterMap( a.fFilterMap ),
    fTPCnClust(a.fTPCnClust),
    fTPCdEdx(a.fTPCdEdx), 
    fTOFbeta( a.fTOFbeta )
{ 
  //copy constructor
  for(Int_t i=0;i<kNAliJTrkPID;i++) fTrkPID[i][kTOF] = a.fTrkPID[i][kTOF];
  for(Int_t i=0;i<kNAliJTrkPID;i++) fTrkPID[i][kTPC] = a.fTrkPID[i][kTPC];
  for(Int_t i=0;i<kNAliJTrkPID;i++) fTrkPID[i][kTPCTOF] = a.fTrkPID[i][kTPCTOF];
  for(Int_t i=0;i<kNAliJTrkPID;i++) fExpTOFbeta[i]= a.fExpTOFbeta[i];
  for( int i=0;i<kNAliJTrkPID;i++) fExpTPCdEdx[i]= a.fExpTPCdEdx[i];
  for( int i=0;i<kNAliJTrkPID;i++) fTPCsigma[i]= a.fTPCsigma[i];
  for( int i=0;i<kNAliJTrkPID;i++) fTOFsigma[i]= a.fTOFsigma[i];
  for( int i=0;i<3;i++) fTPCTrack[i] = a.fTPCTrack[i];
}


//______________________________________________________________________________
AliJTrack&  AliJTrack::operator=(const AliJTrack& trk){
  //operator = 
  if(this != &trk){
    AliJBaseTrack::operator=(trk);
    for(Int_t i=0;i<kNAliJTrkPID;i++){
      fTrkPID[i][kTOF] = trk.fTrkPID[i][kTOF];
      fTrkPID[i][kTPC] = trk.fTrkPID[i][kTPC];
      fTrkPID[i][kTPCTOF] = trk.fTrkPID[i][kTPCTOF];
      fExpTOFbeta[i]= trk.fExpTOFbeta[i];
      fExpTPCdEdx[i] = trk.fExpTPCdEdx[i];
      fTPCsigma[i]  = trk.fTPCsigma[i];
      fTOFsigma[i]  = trk.fTOFsigma[i];
    }
    for( int i=0;i<3;i++) fTPCTrack[i] = trk.fTPCTrack[i];
    fFilterMap  = trk.fFilterMap;
    fTPCnClust  = trk.fTPCnClust;
    fTPCdEdx = trk.fTPCdEdx;
    fTOFbeta = trk.fTOFbeta;
  }
  return *this;
}

