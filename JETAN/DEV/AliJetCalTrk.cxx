/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//--------------------------------------------------
// CalTrk is used to store Tracks and CaloCells information 
//
// Author: alexandre.shabetai@cern.ch & magali.estienne@subatech.in2p3.fr
//-------------------------------------------------

#include "AliJetCalTrk.h"

#include "AliVCaloCells.h"

ClassImp(AliJetCalTrkTrack)

//////////////////////////////////////////////////////////////////

AliJetCalTrkTrack::AliJetCalTrkTrack():
  TObject(),
  fCalTrkTrackRef(),
  fCalTrkCutFlag(kFALSE),
  fCalTrkSignalFlag(kFALSE), 
  fCalTrkPtCorr(1.)
{
  // Default constructor
}  

//------------------------------------------------------------------------
AliJetCalTrkTrack::AliJetCalTrkTrack(AliVTrack* track, Bool_t cutFlag, Bool_t signalFlag, Float_t ptCorr):
  TObject(),
  fCalTrkTrackRef(track),
  fCalTrkCutFlag(cutFlag),
  fCalTrkSignalFlag(signalFlag),
  fCalTrkPtCorr(ptCorr)
{
  // Constructor 2
}

//------------------------------------------------------------------------
AliJetCalTrkTrack::AliJetCalTrkTrack(AliVParticle* track, Bool_t cutFlag, Bool_t signalFlag, Float_t ptCorr):
  TObject(),
  fCalTrkTrackRef(track),
  fCalTrkCutFlag(cutFlag),
  fCalTrkSignalFlag(signalFlag),
  fCalTrkPtCorr(ptCorr)
{
  // Constructor 3
}

//------------------------------------------------------------------------
AliJetCalTrkTrack::AliJetCalTrkTrack(const AliJetCalTrkTrack& rCalTrk):
  TObject(rCalTrk),
  fCalTrkTrackRef(rCalTrk.fCalTrkTrackRef),
  fCalTrkCutFlag(rCalTrk.fCalTrkCutFlag),
  fCalTrkSignalFlag(rCalTrk.fCalTrkSignalFlag),
  fCalTrkPtCorr(rCalTrk.fCalTrkPtCorr)
{
  // Copy constructor
}

//------------------------------------------------------------------------
AliJetCalTrkTrack& AliJetCalTrkTrack::operator=(const AliJetCalTrkTrack& rhs)
{
  // Assignment
  if (this != &rhs) {
   TObject::operator=(rhs);
   fCalTrkTrackRef   = rhs.fCalTrkTrackRef;
   fCalTrkCutFlag    = rhs.fCalTrkCutFlag;
   fCalTrkSignalFlag = rhs.fCalTrkSignalFlag;
   fCalTrkPtCorr     = rhs.fCalTrkPtCorr;
 }
  return *this;

}

//------------------------------------------------------------------------
void AliJetCalTrkTrack::Clear(Option_t* /*option*/)
{
  // Clear objects
  fCalTrkTrackRef   = 0;
  fCalTrkCutFlag    = 0;
  fCalTrkSignalFlag = 0;
  fCalTrkPtCorr     = 1.;
 
}

//-----------------------------------------------------------------------
void AliJetCalTrkTrack::Print(const Option_t* option) const
{
  cout << "Track: " << option << ", Pt: " << GetPt() << ", Eta: " << GetEta() << ", Phi: " << GetPhi() << endl;

}

//...........................................................................
//***************************************************************************
ClassImp(AliJetCalTrkTrackKine)

  AliJetCalTrkTrackKine::AliJetCalTrkTrackKine():
    AliJetCalTrkTrack(),
    fCalTrkPtReso(1.),
    fCalTrkTrackE(-999.),
    fCalTrkTrackPt(-999.),
    fCalTrkTrackP(-999.),
    fCalTrkTrackPx(-999.),
    fCalTrkTrackPy(-999.),
    fCalTrkTrackPz(-999.)
{
  // Default constructor
}

//------------------------------------------------------------------------
AliJetCalTrkTrackKine::AliJetCalTrkTrackKine(AliVParticle* track, Bool_t cutFlag, Bool_t signalFlag, Float_t ptReso) : 
  AliJetCalTrkTrack(track,cutFlag,signalFlag),
  fCalTrkPtReso(ptReso),
  fCalTrkTrackE(-999.),
  fCalTrkTrackPt(-999.),
  fCalTrkTrackP(-999.),
  fCalTrkTrackPx(-999.),
  fCalTrkTrackPy(-999.),
  fCalTrkTrackPz(-999.)
{
  // Constructor 2
  CalcPx(); CalcPy(); CalcPz(); CalcP(); CalcPt(); CalcE();

}

//------------------------------------------------------------------------
void AliJetCalTrkTrackKine::Clear(Option_t* option)
{
  // Clear objects
  fCalTrkPtReso     = 1.;
  fCalTrkTrackE     = -999;
  fCalTrkTrackPt    = -999;
  fCalTrkTrackP     = -999;
  fCalTrkTrackPx    = -999;
  fCalTrkTrackPy    = -999;
  fCalTrkTrackPz    = -999;
  AliJetCalTrkTrack::Clear(option);

}

//------------------------------------------------------------------------
Float_t AliJetCalTrkTrackKine::CalcE()
{
  // Particle energy
  if(fCalTrkTrackE==-999){
    if ( fCalTrkPtReso != 1 ){
      fCalTrkTrackE = TMath::Sqrt(GetPx()*GetPx()+GetPy()*GetPy()+GetPz()*GetPz()+GetM()*GetM());
    }
    else {fCalTrkTrackE = GetParticle()->E(); }
  }
  
  return fCalTrkTrackE;
  
}

//------------------------------------------------------------------------
Float_t AliJetCalTrkTrackKine::CalcPt()
{
  // Particle transverse momentum
  if(fCalTrkTrackPt==-999){
    if ( fCalTrkPtReso != 1 ){
      fCalTrkTrackPt = TMath::Sqrt(GetPx()*GetPx()+GetPy()*GetPy());
    }
    else {fCalTrkTrackPt = GetParticle()->Pt();}
  }

  return fCalTrkTrackPt;

}

//------------------------------------------------------------------------
Float_t AliJetCalTrkTrackKine::CalcP()
{
  // Particle momentum
  if(fCalTrkTrackP==-999){
    if ( fCalTrkPtReso != 1 ){
      fCalTrkTrackP = TMath::Sqrt(GetPx()*GetPx()+GetPy()*GetPy()+GetPz()*GetPz());
    }
    else {fCalTrkTrackP = GetParticle()->P(); }
  }
  
  return fCalTrkTrackP;

}

//...........................................................................
//***************************************************************************

ClassImp(AliJetCalTrkEvent)

  AliJetCalTrkEvent::AliJetCalTrkEvent():
    TObject(),
    fJetCalTrkTrack(0x0),
    fJetCalTrkCell(0x0),
    fNJetCalTrkTrack(0)
{
  // Default constructor
}

//----------------------------------------------------------------
AliJetCalTrkEvent::AliJetCalTrkEvent(Short_t opt,Bool_t kine,Bool_t kIsHighMult):
  TObject(),
  fJetCalTrkTrack(0x0),
  fJetCalTrkCell(0x0),
  fNJetCalTrkTrack(0)
{
  // Constructor 2
  if (kine==0) {
    // Tracks (real or MC)
    if(opt%2==!0 || opt==0){
      fJetCalTrkTrack = new TClonesArray("AliJetCalTrkTrack", kIsHighMult*3800+200);
    }
  }
  else { // Kine cases
    fJetCalTrkTrack = new TClonesArray("AliJetCalTrkTrackKine", kIsHighMult*3800+200);
  }

}

//----------------------------------------------------------------
AliJetCalTrkEvent::~AliJetCalTrkEvent()
{
 // destructor
 if (fJetCalTrkTrack) delete fJetCalTrkTrack;
 if (fJetCalTrkCell)  delete fJetCalTrkCell;

}
//----------------------------------------------------------------
AliJetCalTrkEvent::AliJetCalTrkEvent(const AliJetCalTrkEvent& rCalTrkEvent):
  TObject(),
  fJetCalTrkTrack(rCalTrkEvent.fJetCalTrkTrack),
  fJetCalTrkCell(rCalTrkEvent.fJetCalTrkCell),
  fNJetCalTrkTrack(rCalTrkEvent.fNJetCalTrkTrack)
{
  // Copy constructor
}

//----------------------------------------------------------------
AliJetCalTrkEvent& AliJetCalTrkEvent::operator=(const AliJetCalTrkEvent& rhs)
{
  // Assignment
  if (this != &rhs) {
   TObject::operator=(rhs);
   if (fJetCalTrkTrack) delete fJetCalTrkTrack; 
   if (fJetCalTrkCell) delete  fJetCalTrkCell;
   fJetCalTrkTrack  = rhs.fJetCalTrkTrack;
   fJetCalTrkCell   = rhs.fJetCalTrkCell;
   fNJetCalTrkTrack = rhs.fNJetCalTrkTrack;
  }
  
  return *this;

}

//----------------------------------------------------------------
AliJetCalTrkTrack* AliJetCalTrkEvent::AddCalTrkTrack(AliVTrack* track, Bool_t cutFlag, Bool_t signalFlag, Float_t ptCorr)
{
  // Add a track to the CalTrkEvent  
  TClonesArray &tJetCalTrkTrack = *fJetCalTrkTrack ;
  AliJetCalTrkTrack *n = new(tJetCalTrkTrack[fNJetCalTrkTrack++]) AliJetCalTrkTrack(track, cutFlag, signalFlag, ptCorr) ;
  return n ;

}

//----------------------------------------------------------------
AliJetCalTrkTrack* AliJetCalTrkEvent::AddCalTrkTrack(AliVParticle* track, Bool_t cutFlag, Bool_t signalFlag, Float_t ptCorr)
{
  // Add a track to the CalTrkEvent
  TClonesArray &tJetCalTrkTrack = *fJetCalTrkTrack ;
  AliJetCalTrkTrack *n = new(tJetCalTrkTrack[fNJetCalTrkTrack++]) AliJetCalTrkTrack(track, cutFlag, signalFlag, ptCorr) ;
  return n ;

}

//_________________________________________________________________
AliJetCalTrkTrackKine* AliJetCalTrkEvent::AddCalTrkTrackKine(AliVParticle* track, Bool_t cutFlag, Bool_t signalFlag, Float_t ptReso)
{
  // Add a track to the CalTrkEvent
  TClonesArray &tJetCalTrkTrack = *fJetCalTrkTrack ;
  AliJetCalTrkTrackKine *n = new(tJetCalTrkTrack[fNJetCalTrkTrack++]) AliJetCalTrkTrackKine(track, cutFlag, signalFlag, ptReso) ;
  return n ;

}


//----------------------------------------------------------------
AliJetCalTrkTrack* AliJetCalTrkEvent::GetCalTrkTrack(Int_t i)
{
  // Get track i
  return (AliJetCalTrkTrack*) fJetCalTrkTrack->At(i);

}

//-----------------------------------------------------------------
void AliJetCalTrkEvent::Clear(Option_t* /*option*/)
{
  // Clear object

  if(fJetCalTrkTrack)  fJetCalTrkTrack->Clear("C"); // array of Tracks
  fNJetCalTrkTrack     = 0; // Number of tracks
}

//________________________________________________________________
void  AliJetCalTrkEvent::Print(const Option_t* /*option*/) const
{
  // prints event information
  cout<< "Number of tracks:" << fNJetCalTrkTrack << endl;

}
