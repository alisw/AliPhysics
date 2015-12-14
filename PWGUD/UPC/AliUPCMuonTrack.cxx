
//_____________________________________________________________________________
//    Class for UPC muon track
//    Author: Jaroslav Adam
//
//    contains parameters of the muon track relevant for UPC analysis
//_____________________________________________________________________________

#include "TLorentzVector.h"
#include "TArrayI.h"
#include "TArrayD.h"

#include "AliUPCMuonTrack.h"

ClassImp(AliUPCMuonTrack)

//_____________________________________________________________________________
AliUPCMuonTrack::AliUPCMuonTrack()
 :TObject(),
  fPt(0), fEta(0), fPhi(0), fCharge(0), fMatchTrigger(0), fRabs(0),
  fChi2perNDF(0), fDca(0), fPdca(0),
  fArrayInt(0x0), fArrayD(0x0),
  fkMuonMass(0.105658)
{

  // Default constructor


}

//_____________________________________________________________________________
AliUPCMuonTrack::~AliUPCMuonTrack()
{
  //destructor

  if(fArrayInt) {delete fArrayInt; fArrayInt = 0x0;}
  if(fArrayD) {delete fArrayD; fArrayD = 0x0;}
}

//_____________________________________________________________________________
void AliUPCMuonTrack::Clear(Option_t *)
{
  // clear the track

  if(fArrayInt) fArrayInt->Reset();
  if(fArrayD) fArrayD->Reset();
}

void AliUPCMuonTrack::GetMomentum(TLorentzVector *v) const
{
  // get track 4-momentum
  v->SetPtEtaPhiM(fPt,fEta,fPhi,fkMuonMass);
}

//_____________________________________________________________________________
Int_t AliUPCMuonTrack::MakeArrayInt(Int_t size)
{
  //create the object of TArrayI, skip if it has been already created
  //the array allows the event to be extended for other integer variables

  if( fArrayInt ) return 999; // already initialized

  fArrayInt = new TArrayI(size);

  return 0;
}

//_____________________________________________________________________________
Int_t AliUPCMuonTrack::MakeArrayD(Int_t size)
{
  //create the object of TArrayD, skip if it has been already created
  //the array allows the event to be extended for other double variables

  if( fArrayD ) return 999; // already initialized

  fArrayD = new TArrayD(size);

  return 0;
}










































