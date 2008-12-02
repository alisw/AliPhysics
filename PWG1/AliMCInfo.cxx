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
#include "AliMathBase.h" 
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
  fTPCtrackLength(-1),
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
  fTPCtrackLength(info.fTPCtrackLength),
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


AliMCInfo& AliMCInfo::operator=(const AliMCInfo& info) { 
  //
  // Assignment operator
  //
  this->~AliMCInfo();
  new (this) AliMCInfo(info);
  return *this;
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
	fPrim = AliMathBase::BetheBlochAleph(betagama);
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

void AliMCInfo::CalcTPCrows(TClonesArray * runArrayTR){
  //
  // Calculates the numebr of the track references for detectors
  // In case of changing direction - curling tracks - the counter is not increasing
  //
  // Rough calculation 
  // of the first and last point in the TPC  
  //
  fNTPCRef = 0;
  fNITSRef = 0;
  fNTRDRef = 0;
  fNTOFRef = 0;  
  Float_t tpcminRadius=250;
  Float_t tpcmaxRadius=80;
  Float_t dir=0;
  Int_t nover=0;
  
  for (Int_t iTrackRef = 0; iTrackRef < runArrayTR->GetEntriesFast(); iTrackRef++) {
    //
    AliTrackReference *ref = (AliTrackReference*)runArrayTR->At(iTrackRef);  
    Float_t newdirection = (ref->X()*ref->Px()+ref->Y()*ref->Py()>0)? 1.:-1.; //inside or outside
    if (dir*newdirection<0.5) {
      nover++;
      dir = newdirection;
    }
    //
    if (ref->DetectorId()== AliTrackReference::kTRD){
      tpcmaxRadius =250;
      fNTRDRef++;
    }
    if (ref->DetectorId()== AliTrackReference::kITS){
      fNITSRef++;
      tpcminRadius =90;
    }
    if (ref->DetectorId()== AliTrackReference::kITS){
      fNTOFRef++;
      tpcmaxRadius =250;
    }
    //
    if (ref->DetectorId()== AliTrackReference::kTPC){
      fNTPCRef++;
      if (ref->R()>tpcmaxRadius) tpcmaxRadius = ref->R();
      if (ref->R()<tpcminRadius) tpcminRadius = ref->R();
    }
    if (ref->DetectorId()== AliTrackReference::kDisappeared){
      if (TMath::Abs(ref->Z())<250 && TMath::Abs(ref->R()<250))
	tpcmaxRadius = ref->R();
    }
  }
  fTPCtrackLength = tpcmaxRadius-tpcminRadius;
  fMCtracks=nover;
}
