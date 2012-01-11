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
#include "TDatabasePDG.h"
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
  fNTPCRefOut(0),                    // tpc references counter
  fNITSRefOut(0),                    // ITS references counter
  fNTRDRefOut(0),                    // TRD references counter
  fNTOFRefOut(0),                    // TOF references counter
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
  for (Int_t i=0;i<4;i++) fVDist[i]=0;
  for (Int_t i=0;i<3;i++) fDecayCoord[i]=0;
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
  fNTPCRefOut(info.fNTPCRefOut),                    // tpc references counter
  fNITSRefOut(info.fNITSRefOut),                    // ITS references counter
  fNTRDRefOut(info.fNTRDRefOut),                    // TRD references counter
  fNTOFRefOut(info.fNTOFRefOut),                    // TOF references counter
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
  for (Int_t i=0;i<4;i++) fVDist[i]=info.fVDist[i];
  for (Int_t i=0;i<3;i++) fDecayCoord[i]=info.fDecayCoord[i];

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
  for (Int_t i = 0; i<32; i++) { fDig[i] = 0; }
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
    //fDig[i] <<= 8;
    fDig[i] = 0;
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


void AliMCInfo::Clear(Option_t* ){
  //
  // clear the content
  //
  if (fTPCReferences) fTPCReferences->Clear();
  if (fITSReferences) fITSReferences->Clear();
  if (fTRDReferences) fTRDReferences->Clear();
  if (fTOFReferences) fTOFReferences->Clear();
  for (Int_t i=0;i<4;i++) fVDist[i]=0;
  for (Int_t i=0;i<3;i++) fDecayCoord[i]=0;
  fTRdecay.SetTrack(-1);
  fCharge = 0;
  fEventNr=-1;
  fPrimPart=-1;               // index of primary particle in TreeH
  fMass=0;                   // mass of the particle
  fCharge=0;                 // charge of the particle
  fLabel=0;                  // track label
  fEventNr=0;                // event number
  fMCtracks=0;               // indication of how many times the track is retuturned back
  fPdg=0;                        //pdg code
  fTPCdecay=0;                  //indicates decay in TPC
  fRowsWithDigitsInn=0;          // number of rows with digits in the inner sectors
  fRowsWithDigits=0;             // number of rows with digits in the outer sectors
  fRowsTrackLength=0;            // last - first row with digit
  fTPCtrackLength=0;           // distance between first and last track reference
  fPrim=0;                     // theoretical dedx in tpc according particle momenta and mass
  fNTPCRef=0;                    // tpc references counter
  fNITSRef=0;                    // ITS references counter
  fNTRDRef=0;                    // TRD references counter
  fNTOFRef=0;                    // TOF references counter
  //
  fNTPCRefOut=0;                    // tpc references counter - out
  fNITSRefOut=0;                    // ITS references counter - out
  fNTRDRefOut=0;                    // TRD references counter - out
  fNTOFRefOut=0;                    // TOF references counter - out
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
  fNTPCRefOut = 0;
  fNITSRefOut = 0;
  fNTRDRefOut = 0;
  fNTOFRefOut = 0;  
  //
  //
  Float_t tpcminRadius=TMath::Max(90., fParticle.R());
  Float_t tpcmaxRadius=TMath::Min(250.,fParticle.R());  
  Int_t nover=-1;

  if (runArrayTR->GetEntriesFast()>0){  
    nover=0;
    AliTrackReference *ref0 = (AliTrackReference*)runArrayTR->At(0);    
    Float_t dir = ((ref0->X()*ref0->Px()+ref0->Y()*ref0->Py())>0)? 1.:-1.; //inside or outside
    //
    //    
    for (Int_t iTrackRef = 0; iTrackRef < runArrayTR->GetEntriesFast(); iTrackRef++) {
      //
      AliTrackReference *ref = (AliTrackReference*)runArrayTR->At(iTrackRef);  
      Float_t newdirection = ((ref->X()*ref->Px()+ref->Y()*ref->Py())>0)? 1.:-1.; //inside or outside
      //
      //
      if (dir*newdirection<0.0) {
	nover++;
	dir = newdirection;
      }
      //
      if (ref->DetectorId()== AliTrackReference::kTRD){
	if (nover==0) {
	  tpcmaxRadius =250;
	  fNTRDRefOut++;
	}
	fNTRDRef++;
      }
      if (ref->DetectorId()== AliTrackReference::kITS){ 
	if (nover==0) {
	  fNITSRefOut++;
	  tpcminRadius =90;
	}
	fNITSRef++;
      }
      if (ref->DetectorId()== AliTrackReference::kTOF){
	fNTOFRef++;
	if (nover==0) {
	  fNTOFRefOut++;
	  tpcmaxRadius =250;
	}
      }
      //
      if (ref->DetectorId()== AliTrackReference::kTPC){
	fNTPCRef++;
	if (nover==0) {
	  fNTPCRefOut++;
	  fTrackRefOut = (*ref);
	  if (ref->R()>tpcmaxRadius) tpcmaxRadius = ref->R();
	  if (ref->R()<tpcminRadius) tpcminRadius = ref->R();
	}
      }
      if (ref->DetectorId()== AliTrackReference::kDisappeared){
	if (nover==0 &&TMath::Abs(ref->Z())<250){
	  tpcmaxRadius = TMath::Min(Float_t(250.),ref->R());
	}
      }
    }
  }  // if we have track refs
  fTPCtrackLength = tpcmaxRadius-tpcminRadius;
  fMCtracks=nover+1;
}

//
//
//



void AliMCInfo::Update(TParticle * part, TClonesArray * runArrayTR, Double_t pvertex[4], Int_t label){
  //
  //
  //
  Clear();
  fLabel=label;
  fParticle=(*part);
  AliMCInfo *info = this;
  fVDist[0]=part->Vx()-pvertex[0];
  fVDist[1]=part->Vy()-pvertex[1];
  fVDist[2]=part->Vz()-pvertex[2];
  fVDist[3]=part->R()-pvertex[3];
  //
  //
  //
  for (Int_t iTrackRef = 0; iTrackRef < runArrayTR->GetEntriesFast(); iTrackRef++) {
    AliTrackReference *trackRef = (AliTrackReference*)runArrayTR->At(iTrackRef);         
    if (trackRef->DetectorId()== AliTrackReference::kTPC){     
      TClonesArray & arr = *(info->fTPCReferences);
      new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
    }
    if (trackRef->DetectorId()== AliTrackReference::kITS){
      TClonesArray & arr = *(info->fITSReferences);
      new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
    }
    if (trackRef->DetectorId()== AliTrackReference::kTRD){
      TClonesArray & arr = *(info->fTRDReferences);
      new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
    }
    if (trackRef->DetectorId()== AliTrackReference::kTOF){
      TClonesArray & arr = *(info->fTOFReferences);
      new (arr[arr.GetEntriesFast()]) AliTrackReference(*trackRef);     
    }
    //
    // decay case
    //
    if (trackRef->DetectorId()== AliTrackReference::kDisappeared){
      info->fTRdecay = *trackRef;      	
    }
  }
  //
  Update();
  CalcTPCrows(runArrayTR);
}

