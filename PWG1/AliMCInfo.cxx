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


///////////////////////////////////////////////////////////////////////////
/*

Origin: marian.ivanov@cern.ch
Container classes with MC infomation

The AliMCInfo contains the information about the particles properties 
during transportation throuch ALICE Detector

The base Information :
TParticle - fParticle             -  properties of the particle at creation point
AliTrackReference - fXXXRefernces -  TClonesArray of refernces in differnt detectors
fNXXXRef                          -  number of the track refernces in differnt detectors
AliTPCdigitRow    - fTPCRow       -  the map of the hitted rows - (will be repalced by TBits)
fRowsWith*                        - number of rows hitted by particle
fMCtracks         -               - number of turn over of the track inside of the TPC

++++
some additional information usable for tree draw - TO SPEED UP tree queries
IMPORTANT FOR PROOF FAST PROTOTYPING ANALYSIS 
                                      



*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
//ROOT includes
#include "Rtypes.h"
#include "TClonesArray.h"
//ALIROOT includes
#include "AliTrackReference.h"
#include "AliMCInfo.h" 
#endif

//
// 

ClassImp(AliTPCdigitRow)
ClassImp(AliMCInfo)



////////////////////////////////////////////////////////////////////////
AliMCInfo::AliMCInfo():
  fTrackRef(),
  fTrackRefOut(),
  fTRdecay(),
  fPrimPart(0),
  fParticle(),
  fMass(0),
  fCharge(0),
  fLabel(0),
  fEventNr(),
  fMCtracks(),
  fPdg(0),
  fTPCdecay(0),
  fRowsWithDigitsInn(0),
  fRowsWithDigits(0),
  fRowsTrackLength(0),
  fPrim(0),
  fTPCRow(), 
  fNTPCRef(0),                    // tpc references counter
  fNITSRef(0),                    // ITS references counter
  fNTRDRef(0),                    // TRD references counter
  fNTOFRef(0),                    // TOF references counter
  fTPCReferences(0),
  fITSReferences(0),
  fTRDReferences(0),
  fTOFReferences(0)
{
  //
  // Default constructor
  //
  fTPCReferences  = new TClonesArray("AliTrackReference",10);
  fITSReferences  = new TClonesArray("AliTrackReference",10);
  fTRDReferences  = new TClonesArray("AliTrackReference",10);
  fTOFReferences  = new TClonesArray("AliTrackReference",10);
  fTRdecay.SetTrack(-1);
  fCharge = 0;
}

AliMCInfo::AliMCInfo(const AliMCInfo& info):
  TObject(),
  fTrackRef(info.fTrackRef),
  fTrackRefOut(info.fTrackRefOut),
  fTRdecay(info.GetTRdecay()),
  fPrimPart(info.fPrimPart),
  fParticle(info.fParticle),
  fMass(info.GetMass()),
  fCharge(info.fCharge),
  fLabel(info.fLabel),
  fEventNr(info.fEventNr),
  fMCtracks(info.fMCtracks),
  fPdg(info.fPdg),
  fTPCdecay(info.fTPCdecay),
  fRowsWithDigitsInn(info.fRowsWithDigitsInn),
  fRowsWithDigits(info.fRowsWithDigits),
  fRowsTrackLength(info.fRowsTrackLength),
  fPrim(info.fPrim),
  fTPCRow(info.fTPCRow), 
  fNTPCRef(info.fNTPCRef),                    // tpc references counter
  fNITSRef(info.fNITSRef),                    // ITS references counter
  fNTRDRef(info.fNTRDRef),                    // TRD references counter
  fNTOFRef(info.fNTOFRef),                    // TOF references counter
  fTPCReferences(0),
  fITSReferences(0),
  fTRDReferences(0),
  fTOFReferences(0)
{
  //
  // copy constructor
  //
  fTPCReferences = (TClonesArray*)info.fTPCReferences->Clone();
  fITSReferences = (TClonesArray*)info.fITSReferences->Clone();
  fTRDReferences = (TClonesArray*)info.fTRDReferences->Clone();
  fTOFReferences = (TClonesArray*)info.fTOFReferences->Clone();
}


AliMCInfo::~AliMCInfo()
{
  //
  // Destructor of the class
  //
  if (fTPCReferences) {
    delete fTPCReferences;
  }
  if (fITSReferences){
    delete fITSReferences;
  }
  if (fTRDReferences){    
    delete fTRDReferences;  
  }
  if (fTOFReferences){    
    delete fTOFReferences;  
  }

}



void AliMCInfo::Update()
{
  //
  // Update MC info
  // Calculates some derived variables
  //
  //
  fMCtracks =1;
  if (!fTPCReferences) {
    fNTPCRef =0;
    return;
  }
  Float_t direction=1;
  //Float_t rlast=0;
  fNTPCRef = fTPCReferences->GetEntriesFast();
  fNITSRef = fITSReferences->GetEntriesFast();
  fNTRDRef = fTRDReferences->GetEntriesFast();
  fNTOFRef = fTOFReferences->GetEntriesFast();
  
  for (Int_t iref =0;iref<fTPCReferences->GetEntriesFast();iref++){
    AliTrackReference * ref = (AliTrackReference *) fTPCReferences->At(iref);
    //Float_t r = (ref->X()*ref->X()+ref->Y()*ref->Y());
    Float_t newdirection = ref->X()*ref->Px()+ref->Y()*ref->Py(); //inside or outside
    if (iref==0) direction = newdirection;
    if ( newdirection*direction<0){
      //changed direction
      direction = newdirection;
      fMCtracks+=1;
    }
    //rlast=r;			    
  }
  //
  // decay info
  fTPCdecay=kFALSE;
  if (fTRdecay.GetTrack()>0){
    fDecayCoord[0] = fTRdecay.X();
    fDecayCoord[1] = fTRdecay.Y();
    fDecayCoord[2] = fTRdecay.Z();
    if ( (fTRdecay.R()<250)&&(fTRdecay.R()>85) && (TMath::Abs(fTRdecay.Z())<250) ){
      fTPCdecay=kTRUE;     
    }
    else{
      fDecayCoord[0] = 0;
      fDecayCoord[1] = 0;
      fDecayCoord[2] = 0;
    }
  }
  //
  //
  //digits information update
  fRowsWithDigits    = fTPCRow.RowsOn();    
  fRowsWithDigitsInn = fTPCRow.RowsOn(63); // 63 = number of inner rows
  fRowsTrackLength   = fTPCRow.Last() - fTPCRow.First();
  //
  //
  // calculate primary ionization per cm  
  if (fParticle.GetPDG()){
    fMass = fParticle.GetMass();  
    fCharge = fParticle.GetPDG()->Charge();
    if (fTPCReferences->GetEntriesFast()>0){
      fTrackRef = *((AliTrackReference*)fTPCReferences->At(0));
    }
    if (fMass>0){
      Float_t p = TMath::Sqrt(fTrackRef.Px()*fTrackRef.Px()+
			      fTrackRef.Py()*fTrackRef.Py()+
			      fTrackRef.Pz()*fTrackRef.Pz());    
      if (p>0.001){
	Float_t betagama = p /fMass;
	fPrim = TPCBetheBloch(betagama);
      }else fPrim=0;
    }
  }else{
    fMass =0;
    fPrim =0;
  }  
}




////////////////////////////////////////////////////////////////////////
AliTPCdigitRow::AliTPCdigitRow()
{
  Reset();
}
////////////////////////////////////////////////////////////////////////
AliTPCdigitRow & AliTPCdigitRow::operator=(const AliTPCdigitRow &digOld)
{
  for (Int_t i = 0; i<32; i++) fDig[i] = digOld.fDig[i];
  return (*this);
}
////////////////////////////////////////////////////////////////////////
void AliTPCdigitRow::SetRow(Int_t row) 
{
  //
  // set bit mask for given row
  //
  if (row >= 8*32) {
    //    cerr<<"AliTPCdigitRow::SetRow: index "<<row<<" out of bounds."<<endl;
    return;
  }
  Int_t iC = row/8;
  Int_t iB = row%8;
  SETBIT(fDig[iC],iB);
}

////////////////////////////////////////////////////////////////////////
Bool_t AliTPCdigitRow::TestRow(Int_t row) const
{
//
// return kTRUE if row is on
//
  Int_t iC = row/8;
  Int_t iB = row%8;
  return TESTBIT(fDig[iC],iB);
}
////////////////////////////////////////////////////////////////////////
Int_t AliTPCdigitRow::RowsOn(Int_t upto) const
{
//
// returns number of rows with a digit  
// count only rows less equal row number upto
//
  Int_t total = 0;
  for (Int_t i = 0; i<32; i++) {
    for (Int_t j = 0; j < 8; j++) {
      if (i*8+j > upto) return total;
      if (TESTBIT(fDig[i],j))  total++;
    }
  }
  return total;
}
////////////////////////////////////////////////////////////////////////
void AliTPCdigitRow::Reset()
{
//
// resets all rows to zero
//
  for (Int_t i = 0; i<32; i++) {
    fDig[i] <<= 8;
  }
}
////////////////////////////////////////////////////////////////////////
Int_t AliTPCdigitRow::Last() const
{
//
// returns the last row number with a digit
// returns -1 if now digits 
//
  for (Int_t i = 32-1; i>=0; i--) {
    for (Int_t j = 7; j >= 0; j--) {
      if TESTBIT(fDig[i],j) return i*8+j;
    }
  }
  return -1;
}
////////////////////////////////////////////////////////////////////////
Int_t AliTPCdigitRow::First() const
{
//
// returns the first row number with a digit
// returns -1 if now digits 
//
  for (Int_t i = 0; i<32; i++) {
    for (Int_t j = 0; j < 8; j++) {
      if (TESTBIT(fDig[i],j)) return i*8+j;
    }
  }
  return -1;
}


//_____________________________________________________________________________
Float_t AliMCInfo::TPCBetheBloch(Float_t bg)
{
  //
  // Bethe-Bloch energy loss formula
  //
  const Double_t kp1=0.76176e-1;
  const Double_t kp2=10.632;
  const Double_t kp3=0.13279e-4;
  const Double_t kp4=1.8631;
  const Double_t kp5=1.9479;

  Double_t dbg = (Double_t) bg;

  Double_t beta = dbg/TMath::Sqrt(1.+dbg*dbg);

  Double_t aa = TMath::Power(beta,kp4);
  Double_t bb = TMath::Power(1./dbg,kp5);

  bb=TMath::Log(kp3+bb);
  
  return ((Float_t)((kp2-aa-bb)*kp1/aa));
}

