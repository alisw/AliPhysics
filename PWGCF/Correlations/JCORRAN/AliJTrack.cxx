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
    fTOFsignal(9999),
    fTPCmom(0)
{
  // default constructor
  for( int i=0;i<kNAliJTrkPID;i++) fExpTOFsignal[i]= 9999;
  for( int i=0;i<3;i++) fTPCTrack[i] = 0;
  for( int i=0;i<3;i++) fGCGTrack[i] = 0;
  for( int i=0;i<3;i++) fTrackPos[i] = 0;

}

//______________________________________________________________________________
AliJTrack::AliJTrack(const AliJTrack& a):
    AliJBaseTrack(a),
    fFilterMap( a.fFilterMap ),
    fTPCnClust(a.fTPCnClust),
    fTPCdEdx(a.fTPCdEdx), 
    fTOFsignal(a.fTOFsignal),
    fTPCmom(a.fTPCmom)
{ 
  //copy constructor
  for(Int_t i=0;i<kNAliJTrkPID;i++) fExpTOFsignal[i]= a.fExpTOFsignal[i];
  for( int i=0;i<3;i++) fTPCTrack[i] = a.fTPCTrack[i];
  for( int i=0;i<3;i++) fGCGTrack[i] = a.fGCGTrack[i];
  for( int i=0;i<3;i++) fTrackPos[i] = a.fTrackPos[i];
}


//______________________________________________________________________________
AliJTrack&  AliJTrack::operator=(const AliJTrack& trk){
  //operator = 
  if(this != &trk){
    AliJBaseTrack::operator=(trk);
    for(Int_t i=0;i<kNAliJTrkPID;i++){
	 fExpTOFsignal[i]= trk.fExpTOFsignal[i];
    }
    for( int i=0;i<3;i++) fTPCTrack[i] = trk.fTPCTrack[i];
    for( int i=0;i<3;i++) fGCGTrack[i] = trk.fGCGTrack[i];
    for( int i=0;i<3;i++) fTrackPos[i] = trk.fTrackPos[i];
    fFilterMap  = trk.fFilterMap;
    fTPCnClust  = trk.fTPCnClust;
    fTPCdEdx = trk.fTPCdEdx;
    fTOFsignal = trk.fTOFsignal;
  }
  return *this;
}

