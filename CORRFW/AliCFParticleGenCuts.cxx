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

////////////////////////////////////////////////////////////////////////////
//          ----   CORRECTION FRAMEWORK   ----
// class AliCFParticleGenCuts implementation
// Using this class a user may define selections relative to 
// MC particle (AliMCParticle) using generation-level information.
////////////////////////////////////////////////////////////////////////////
// author : R. Vernet (renaud.vernet@cern.ch)
////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliCFParticleGenCuts.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "AliMCEvent.h"
#include "TObject.h"
#include "AliStack.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TBits.h"
#include "TList.h"
#include "TArrayF.h"
#include "TDecayChannel.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"

ClassImp(AliCFParticleGenCuts)

//______________________________
AliCFParticleGenCuts::AliCFParticleGenCuts() : 
  AliCFCutBase(),
  fIsAODMC(0),
  fMCInfo(0x0),
  fRequireIsCharged(0),
  fRequireIsNeutral(0),
  fRequireIsPrimary(0),
  fRequireIsSecondary(0),
  fRequirePdgCode(0),
  fRequireAbsolutePdg(0),
  fProdVtxRange2D(0),
  fPdgCode(0),
  fProdVtxXMin (-1.e+09),
  fProdVtxYMin (-1.e+09),
  fProdVtxZMin (-1.e+09),
  fProdVtxXMax ( 1.e+09),
  fProdVtxYMax ( 1.e+09),
  fProdVtxZMax ( 1.e+09),
  fDecayVtxXMin(-1.e+09),
  fDecayVtxYMin(-1.e+09),
  fDecayVtxZMin(-1.e+09),
  fDecayVtxXMax( 1.e+09),
  fDecayVtxYMax( 1.e+09),
  fDecayVtxZMax( 1.e+09),
  fDecayLengthMin(-1.),
  fDecayLengthMax(1.e+09),
  fDecayRxyMin(-1),
  fDecayRxyMax(1.e+09),
  fDecayChannel(0x0),
  fhCutStatistics(0x0),
  fhCutCorrelation(0x0),
  fCutValues(new TArrayF(kNCuts)),
  fBitmap(new TBits(0))
{
  //
  //ctor
  //
  for (int i=0; i<kNCuts; i++) 
    for (int j=0; j<kNStepQA; j++) 
      fhQA[i][j]=0x0;

  for (int j=0; j<kNStepQA; j++)
    fhProdVtxXY[j]=0x0;
}

//______________________________
AliCFParticleGenCuts::AliCFParticleGenCuts(const Char_t* name, const Char_t* title) : 
  AliCFCutBase(name,title),
  fIsAODMC(0),
  fMCInfo(0x0),
  fRequireIsCharged(0),
  fRequireIsNeutral(0),
  fRequireIsPrimary(0),
  fRequireIsSecondary(0),
  fRequirePdgCode(0),
  fRequireAbsolutePdg(0),
  fProdVtxRange2D(0),
  fPdgCode(0),
  fProdVtxXMin (-1.e+09),
  fProdVtxYMin (-1.e+09),
  fProdVtxZMin (-1.e+09),
  fProdVtxXMax ( 1.e+09),
  fProdVtxYMax ( 1.e+09),
  fProdVtxZMax ( 1.e+09),
  fDecayVtxXMin(-1.e+09),
  fDecayVtxYMin(-1.e+09),
  fDecayVtxZMin(-1.e+09),
  fDecayVtxXMax( 1.e+09),
  fDecayVtxYMax( 1.e+09),
  fDecayVtxZMax( 1.e+09),
  fDecayLengthMin(-1.),
  fDecayLengthMax(1.e+09),
  fDecayRxyMin(-1.),
  fDecayRxyMax(1.e+09),
  fDecayChannel(0x0),
  fhCutStatistics(0x0),
  fhCutCorrelation(0x0),
  fCutValues(new TArrayF(kNCuts)),
  fBitmap(new TBits(0))
{
  //
  //ctor
  //
  for (int i=0; i<kNCuts; i++) 
    for (int j=0; j<kNStepQA; j++) 
      fhQA[i][j]=0x0;

  for (int j=0; j<kNStepQA; j++)
    fhProdVtxXY[j]=0x0;
}

//______________________________
AliCFParticleGenCuts::AliCFParticleGenCuts(const AliCFParticleGenCuts& c) : 
  AliCFCutBase(c),
  fIsAODMC(c.fIsAODMC),
  fMCInfo(c.fMCInfo),
  fRequireIsCharged(c.fRequireIsCharged),
  fRequireIsNeutral(c.fRequireIsNeutral),
  fRequireIsPrimary(c.fRequireIsPrimary),
  fRequireIsSecondary(c.fRequireIsSecondary),
  fRequirePdgCode(c.fRequirePdgCode),
  fRequireAbsolutePdg(c.fRequireAbsolutePdg),
  fProdVtxRange2D(c.fProdVtxRange2D),
  fPdgCode(c.fPdgCode),
  fProdVtxXMin (c.fProdVtxXMin),
  fProdVtxYMin (c.fProdVtxYMin),
  fProdVtxZMin (c.fProdVtxZMin),
  fProdVtxXMax (c.fProdVtxXMax),
  fProdVtxYMax (c.fProdVtxYMax),
  fProdVtxZMax (c.fProdVtxZMax),
  fDecayVtxXMin(c.fDecayVtxXMin),
  fDecayVtxYMin(c.fDecayVtxYMin),
  fDecayVtxZMin(c.fDecayVtxZMin),
  fDecayVtxXMax(c.fDecayVtxXMax),
  fDecayVtxYMax(c.fDecayVtxYMax),
  fDecayVtxZMax(c.fDecayVtxZMax),
  fDecayLengthMin(c.fDecayLengthMin),
  fDecayLengthMax(c.fDecayLengthMin),
  fDecayRxyMin(c.fDecayLengthMin),
  fDecayRxyMax(c.fDecayLengthMin),
  fDecayChannel(c.fDecayChannel),
  fhCutStatistics(new TH1F(*c.fhCutStatistics)),
  fhCutCorrelation(new TH2F(*c.fhCutCorrelation)),
  fCutValues(new TArrayF(*c.fCutValues)),
  fBitmap(new TBits(*c.fBitmap))
{
  //
  //copy ctor
  //
  for (int i=0; i<kNCuts; i++)
    for (int j=0; j<kNStepQA; j++)
      fhQA[i][j]=(TH1F*)c.fhQA[i][j]->Clone();

  for (int j=0; j<kNStepQA; j++)
    fhProdVtxXY[j]=(TH2F*)c.fhProdVtxXY[j]->Clone();
}

//______________________________
AliCFParticleGenCuts& AliCFParticleGenCuts::operator=(const AliCFParticleGenCuts& c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    fIsAODMC=c.fIsAODMC;
    fMCInfo=c.fMCInfo;
    fRequireIsCharged=c.fRequireIsCharged;
    fRequireIsNeutral=c.fRequireIsNeutral;
    fRequireIsPrimary=c.fRequireIsPrimary;
    fRequireIsSecondary=c.fRequireIsSecondary;
    fRequirePdgCode=c.fRequirePdgCode;
    fRequireAbsolutePdg=c.fRequireAbsolutePdg;
    fProdVtxRange2D=c.fProdVtxRange2D;
    fPdgCode=c.fPdgCode;
    fProdVtxXMin=c.fProdVtxXMin;
    fProdVtxYMin=c.fProdVtxYMin;
    fProdVtxZMin=c.fProdVtxZMin;
    fProdVtxXMax=c.fProdVtxXMax;
    fProdVtxYMax=c.fProdVtxYMax;
    fProdVtxZMax=c.fProdVtxZMax;
    fDecayVtxXMin=c.fDecayVtxXMin;
    fDecayVtxYMin=c.fDecayVtxYMin;
    fDecayVtxZMin=c.fDecayVtxZMin;
    fDecayVtxXMax=c.fDecayVtxXMax;
    fDecayVtxYMax=c.fDecayVtxYMax;
    fDecayVtxZMax=c.fDecayVtxZMax;      
    fDecayLengthMin=c.fDecayVtxZMax;
    fDecayLengthMax=c.fDecayLengthMax;
    fDecayRxyMin=c.fDecayRxyMin;
    fDecayRxyMax=c.fDecayRxyMax;
    fDecayChannel=c.fDecayChannel;
    fCutValues=new TArrayF(*c.fCutValues);
    fBitmap=new TBits(*c.fBitmap);
    
    if (fhCutStatistics)  fhCutStatistics =new TH1F(*c.fhCutStatistics) ;
    if (fhCutCorrelation) fhCutCorrelation=new TH2F(*c.fhCutCorrelation);
    
    for (int i=0; i<kNCuts; i++)
      for (int j=0; j<kNStepQA; j++)
	fhQA[i][j]=(TH1F*)c.fhQA[i][j]->Clone();

    for (int j=0; j<kNStepQA; j++)
	fhProdVtxXY[j]=(TH2F*)c.fhProdVtxXY[j]->Clone();
  }
  return *this ;
}

//______________________________
Bool_t AliCFParticleGenCuts::IsSelected(TObject* obj) {
  //
  // check if selections on 'obj' are passed
  // 'obj' must be an AliMCParticle
  //
  
  if (!obj) return kFALSE ;

  if (!fIsAODMC) SelectionBitMap((AliMCParticle*)   obj);
  else           SelectionBitMap((AliAODMCParticle*)obj);

  if (fIsQAOn) FillHistograms(obj,0);

  for (UInt_t icut=0; icut<fBitmap->GetNbits();icut++)
    if (!fBitmap->TestBitNumber(icut)) return kFALSE ; 
  
  if (fIsQAOn) FillHistograms(obj,1);
  return kTRUE;
}
	     
//__________________________________________________________________________________
void AliCFParticleGenCuts::SelectionBitMap(AliMCParticle* mcPart)
{
  //
  // test if the track passes the single cuts
  // and store the information in a bitmap
  //

  for (UInt_t i=0; i<kNCuts; i++) {
    fBitmap->SetBitNumber(i,kFALSE);
    fCutValues->SetAt((Double32_t)0,i) ;
  }

  // fill the cut array
  Double32_t partVx=(Double32_t)mcPart->Xv();
  Double32_t partVy=(Double32_t)mcPart->Yv();
  Double32_t partVz=(Double32_t)mcPart->Zv();

  // calculate the production vertex ellipse
  Double32_t prodVtxXYmin = 0.;
  if (fProdVtxXMin>0 && fProdVtxYMin>0)
	prodVtxXYmin = partVx*partVx/(fProdVtxXMin*fProdVtxXMin) + partVy*partVy/(fProdVtxYMin*fProdVtxYMin);
  Double32_t prodVtxXYmax = 0.;
  if(fProdVtxXMax>0 && fProdVtxYMax>0)
	prodVtxXYmax = partVx*partVx/(fProdVtxXMax*fProdVtxXMax) + partVy*partVy/(fProdVtxYMax*fProdVtxYMax);

  Double32_t decayVx=0.;
  Double32_t decayVy=0.;
  Double32_t decayVz=0.;
  Double32_t decayL=0.;
  Double32_t decayRxy=0.;

  TParticle* part = mcPart->Particle();
  AliStack*  stack = ((AliMCEvent*)fMCInfo)->Stack();
  TParticle* daughter=0x0;
  if ( part->GetNDaughters() > 0 ) {
    daughter = stack->Particle(part->GetFirstDaughter()) ;
    decayVx=(Double32_t)daughter->Vx();
    decayVy=(Double32_t)daughter->Vy();
    decayVz=(Double32_t)daughter->Vz();
    decayL = TMath::Sqrt(TMath::Power(partVx-decayVx,2) + 
			 TMath::Power(partVy-decayVy,2) + 
			 TMath::Power(partVz-decayVz,2) ) ;
    decayRxy = TMath::Sqrt(TMath::Power(decayVx,2) + TMath::Power(decayVy,2) ) ;
  }

  fCutValues->SetAt(partVx  ,kCutProdVtxXMin) ;
  fCutValues->SetAt(partVy  ,kCutProdVtxYMin) ;
  fCutValues->SetAt(partVz  ,kCutProdVtxZMin) ;
  fCutValues->SetAt(partVx  ,kCutProdVtxXMax) ;
  fCutValues->SetAt(partVy  ,kCutProdVtxYMax) ;
  fCutValues->SetAt(partVz  ,kCutProdVtxZMax) ;
  fCutValues->SetAt(decayVx ,kCutDecVtxXMin)  ;
  fCutValues->SetAt(decayVy ,kCutDecVtxYMin)  ;
  fCutValues->SetAt(decayVz ,kCutDecVtxZMin)  ;
  fCutValues->SetAt(decayVx ,kCutDecVtxXMax)  ;
  fCutValues->SetAt(decayVy ,kCutDecVtxYMax)  ;
  fCutValues->SetAt(decayVz ,kCutDecVtxZMax)  ;
  fCutValues->SetAt(decayL  ,kCutDecLgthMin)  ;
  fCutValues->SetAt(decayL  ,kCutDecLgthMax)  ;
  fCutValues->SetAt(decayRxy,kCutDecRxyMin)   ;
  fCutValues->SetAt(decayRxy,kCutDecRxyMax)   ;
  
  // cut on charge
  if ( fRequireIsCharged || fRequireIsNeutral ) {
    if (fRequireIsCharged &&  IsCharged(mcPart)) fCutValues->SetAt((Double32_t)kTRUE,kCutCharge) ;
    if (fRequireIsNeutral && !IsCharged(mcPart)) fCutValues->SetAt((Double32_t)kTRUE,kCutCharge) ;
  } 
  else fCutValues->SetAt((Double32_t)kTRUE,kCutCharge) ;
  
  // cut on primary/secondary
  if ( fRequireIsPrimary || fRequireIsSecondary) {
    if (fRequireIsPrimary   &&  IsPrimary(mcPart)) fCutValues->SetAt((Double32_t)kTRUE,kCutPrimSec);
    if (fRequireIsSecondary && !IsPrimary(mcPart)) fCutValues->SetAt((Double32_t)kTRUE,kCutPrimSec);
  } 
  else fCutValues->SetAt((Double32_t)kTRUE,kCutPrimSec);
  
  // cut on PDG code
  if ( fRequirePdgCode ) {
    if (IsA(mcPart,fPdgCode,fRequireAbsolutePdg)) fCutValues->SetAt((Double32_t)kTRUE,kCutPDGCode);
  }
  else fCutValues->SetAt((Double32_t)kTRUE,kCutPDGCode);
  
  // cut on decay channel
  if ( fDecayChannel ) {
    Bool_t goodDecay = kTRUE ;
    Short_t nDaughters = mcPart->Particle()->GetNDaughters() ;
    if (nDaughters != fDecayChannel->NDaughters()) goodDecay = kFALSE ;
    //now number of daughters is OK
    if (goodDecay) {
      // now check if decay channel is respected
      // first try
      for (Int_t iDaughter = 0; iDaughter<nDaughters; iDaughter++) {
	TParticle* daug = stack->Particle(mcPart->Particle()->GetDaughter(iDaughter)) ;
	if (daug->GetPdgCode() != fDecayChannel->DaughterPdgCode(iDaughter)) {goodDecay = kFALSE; break;}
      }
      if (!goodDecay) {
	//second try inverting the order of the daughters
	goodDecay = kTRUE ;
	for (Int_t iDaughter = 0; iDaughter<nDaughters; iDaughter++) {
	  TParticle* daug = stack->Particle(mcPart->Particle()->GetDaughter(nDaughters-(iDaughter+1))) ;
	  if (daug->GetPdgCode() != fDecayChannel->DaughterPdgCode(iDaughter)) {goodDecay = kFALSE; break;}
	}
      }
      if (!goodDecay && fRequireAbsolutePdg) {
	//now tries inverting the sign of the daughters in case the anti-particle is also looked at
	// third try
	goodDecay = kTRUE ;
	for (Int_t iDaughter = 0; iDaughter<nDaughters; iDaughter++) {
	  TParticle* daug = stack->Particle(mcPart->Particle()->GetDaughter(iDaughter)) ;
	  if (daug->GetPdgCode() != -fDecayChannel->DaughterPdgCode(iDaughter)) {goodDecay = kFALSE; break;}
	}
	if (!goodDecay) {
	  //fourth try inverting the order of the daughters
	  goodDecay = kTRUE ;
	  for (Int_t iDaughter = 0; iDaughter<nDaughters; iDaughter++) {
	    TParticle* daug = stack->Particle(mcPart->Particle()->GetDaughter(nDaughters-(iDaughter+1))) ;
	    if (daug->GetPdgCode() != -fDecayChannel->DaughterPdgCode(iDaughter)) {goodDecay = kFALSE; break;}
	  }
	}
      } //end check anti-particle
    } //end # daughters OK
    fCutValues->SetAt((Double32_t)goodDecay,kCutDecayChannel) ;
  } //end require decay channel
  else fCutValues->SetAt((Double32_t)kTRUE,kCutDecayChannel);
  
  
  // now array of cut is build, fill the bitmap consequently
  Int_t iCutBit = -1;
  if ( fCutValues->At(++iCutBit) !=0 )              fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) !=0 )              fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) !=0 )              fBitmap->SetBitNumber(iCutBit,kTRUE);

  ++iCutBit;
  if ( (!fProdVtxRange2D && fCutValues->At(iCutBit) > fProdVtxXMin)
    || ( fProdVtxRange2D && (fProdVtxXMin>0 && fProdVtxYMin>0) && prodVtxXYmin >= 1)
    || ( fProdVtxRange2D && (fProdVtxXMin<=0 || fProdVtxYMin<=0) ) )
   fBitmap->SetBitNumber(iCutBit,kTRUE);

  ++iCutBit;
  if ( (!fProdVtxRange2D && fCutValues->At(iCutBit) < fProdVtxXMax)
    || ( fProdVtxRange2D && (fProdVtxXMax>0 && fProdVtxYMax>0) && prodVtxXYmax <= 1)
    || ( fProdVtxRange2D && (fProdVtxXMax<=0 || fProdVtxYMax<=0) ) )
  fBitmap->SetBitNumber(iCutBit,kTRUE);

  ++iCutBit;
  if ( (!fProdVtxRange2D && fCutValues->At(iCutBit) > fProdVtxYMin)
    || ( fProdVtxRange2D &&  (fProdVtxXMin>0 && fProdVtxYMin>0) && prodVtxXYmin >= 1)
    || ( fProdVtxRange2D &&  (fProdVtxXMin<=0 || fProdVtxYMin<=0) ) )
  fBitmap->SetBitNumber(iCutBit,kTRUE);

  ++iCutBit;
  if ( (!fProdVtxRange2D && fCutValues->At(iCutBit) < fProdVtxYMax)
    || ( fProdVtxRange2D && (fProdVtxXMax>0 && fProdVtxYMax>0) && prodVtxXYmax <= 1)
    || ( fProdVtxRange2D && (fProdVtxXMax<=0 || fProdVtxYMax<=0) ) )
  fBitmap->SetBitNumber(iCutBit,kTRUE);

  if ( fCutValues->At(++iCutBit) > fProdVtxZMin)    fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) < fProdVtxZMax)    fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) > fDecayVtxXMin)   fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) < fDecayVtxXMax)   fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) > fDecayVtxYMin)   fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) < fDecayVtxYMax)   fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) > fDecayVtxZMin)   fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) < fDecayVtxZMax)   fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) > fDecayLengthMin) fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) < fDecayLengthMax) fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) > fDecayRxyMin)    fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) < fDecayRxyMax)    fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) != 0 )             fBitmap->SetBitNumber(iCutBit,kTRUE);
}

//__________________________________________________________________________________
void AliCFParticleGenCuts::SelectionBitMap(AliAODMCParticle* mcPart)
{
  //
  // test if the track passes the single cuts
  // and store the information in a bitmap
  //
  
  for (UInt_t i=0; i<kNCuts; i++) {
    fBitmap->SetBitNumber(i,kFALSE);
    fCutValues->SetAt((Double32_t)0,i) ;
  }

  // fill the cut array
  Double32_t partVx=(Double32_t)mcPart->Xv();
  Double32_t partVy=(Double32_t)mcPart->Yv();
  Double32_t partVz=(Double32_t)mcPart->Zv();

  // calculate the production vertex ellipse
  Double32_t prodVtxXYmin = 0.;
  if (fProdVtxXMin!=0 && fProdVtxYMin!=0)
	prodVtxXYmin = partVx*partVx/(fProdVtxXMin*fProdVtxXMin) + partVy*partVy/(fProdVtxYMin*fProdVtxYMin);
  Double32_t prodVtxXYmax = 0.;
  if(fProdVtxXMax!=0 && fProdVtxYMax!=0)
	prodVtxXYmax = partVx*partVx/(fProdVtxXMax*fProdVtxXMax) + partVy*partVy/(fProdVtxYMax*fProdVtxYMax);

  Double32_t decayVx=0.;
  Double32_t decayVy=0.;
  Double32_t decayVz=0.;
  Double32_t decayL=0.;
  Double32_t decayRxy=0.;

  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(fMCInfo);
  
  if (!aod) {
    AliError("AOD event casting failed");
    return;
  }
  
  TClonesArray* mcArray = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcArray) {
    AliError("array casting failed");
    return;
  }
  AliAODMCParticle* daughter = 0x0 ;

  if (mcPart->GetDaughter(0)>0) {
    daughter = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcPart->GetDaughter(0)));
    if (!daughter) {
      AliError("daughter casting failed");
      return;
    }

    decayVx=(Double32_t)daughter->Xv();
    decayVy=(Double32_t)daughter->Yv();
    decayVz=(Double32_t)daughter->Zv();
    decayL = TMath::Sqrt(TMath::Power(partVx-decayVx,2) + 
			 TMath::Power(partVy-decayVy,2) + 
			 TMath::Power(partVz-decayVz,2) ) ;
    decayRxy = TMath::Sqrt(TMath::Power(decayVx,2) + TMath::Power(decayVy,2) ) ;
    
  }
  
  fCutValues->SetAt(partVx  ,kCutProdVtxXMin) ;
  fCutValues->SetAt(partVy  ,kCutProdVtxYMin) ;
  fCutValues->SetAt(partVz  ,kCutProdVtxZMin) ;
  fCutValues->SetAt(partVx  ,kCutProdVtxXMax) ;
  fCutValues->SetAt(partVy  ,kCutProdVtxYMax) ;
  fCutValues->SetAt(partVz  ,kCutProdVtxZMax) ;
  fCutValues->SetAt(decayVx ,kCutDecVtxXMin)  ;
  fCutValues->SetAt(decayVy ,kCutDecVtxYMin)  ;
  fCutValues->SetAt(decayVz ,kCutDecVtxZMin)  ;
  fCutValues->SetAt(decayVx ,kCutDecVtxXMax)  ;
  fCutValues->SetAt(decayVy ,kCutDecVtxYMax)  ;
  fCutValues->SetAt(decayVz ,kCutDecVtxZMax)  ;
  fCutValues->SetAt(decayL  ,kCutDecLgthMin)  ;
  fCutValues->SetAt(decayL  ,kCutDecLgthMax)  ;
  fCutValues->SetAt(decayRxy,kCutDecRxyMin)   ;
  fCutValues->SetAt(decayRxy,kCutDecRxyMax)   ;
  
  // cut on charge
  if ( fRequireIsCharged || fRequireIsNeutral ) {
    if (fRequireIsCharged &&  IsCharged(mcPart)) fCutValues->SetAt((Double32_t)kTRUE,kCutCharge) ;
    if (fRequireIsNeutral && !IsCharged(mcPart)) fCutValues->SetAt((Double32_t)kTRUE,kCutCharge) ;
  } 
  else fCutValues->SetAt((Double32_t)kTRUE,kCutCharge) ;
  
  // cut on primary/secondary
  if ( fRequireIsPrimary || fRequireIsSecondary) {
    if (fRequireIsPrimary   &&  IsPrimary(mcPart)) fCutValues->SetAt((Double32_t)kTRUE,kCutPrimSec);
    if (fRequireIsSecondary && !IsPrimary(mcPart)) fCutValues->SetAt((Double32_t)kTRUE,kCutPrimSec);
  } 
  else fCutValues->SetAt((Double32_t)kTRUE,kCutPrimSec);
  
  // cut on PDG code
  if ( fRequirePdgCode ) {
    if (IsA(mcPart,fPdgCode,fRequireAbsolutePdg)) fCutValues->SetAt((Double32_t)kTRUE,kCutPDGCode);
  }
  else fCutValues->SetAt((Double32_t)kTRUE,kCutPDGCode);
  
  // cut on decay channel
  if ( fDecayChannel ) {
    Bool_t goodDecay = kTRUE ;
    Short_t nDaughters = 0 ;
    if (mcPart->GetDaughter(0) >=0) nDaughters = 1 + mcPart->GetDaughter(1) - mcPart->GetDaughter(0);
    
    if (nDaughters != fDecayChannel->NDaughters()) goodDecay = kFALSE ;
    if (goodDecay) {
      // now check if decay channel is respected
      // first try
      for (Int_t iDaughter = 0 ; iDaughter<nDaughters; iDaughter++) {
	AliAODMCParticle* daug = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcPart->GetDaughter(0)+iDaughter));
	if (!daug) {
	  AliError("daughter: casting failed");
	  continue;
	}
	if (daug->GetPdgCode() != fDecayChannel->DaughterPdgCode(iDaughter)) {goodDecay = kFALSE; break;}
      }
      if (!goodDecay) {
	//second try inverting the order of the daughters
	goodDecay = kTRUE ;
	for (Int_t iDaughter = 0 ; iDaughter<nDaughters; iDaughter++) {
	  AliAODMCParticle* daug = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcPart->GetDaughter(1)-iDaughter));
	  if (!daug) AliFatal("");
	  if (daug->GetPdgCode() != fDecayChannel->DaughterPdgCode(iDaughter)) {goodDecay = kFALSE; break;}
	}
      }
    }
    fCutValues->SetAt((Double32_t)goodDecay,kCutDecayChannel) ;
  }
  else fCutValues->SetAt((Double32_t)kTRUE,kCutDecayChannel);
  
  
  // now array of cut is build, fill the bitmap consequently
  Int_t iCutBit = -1;
  if ( fCutValues->At(++iCutBit) !=0 )              fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) !=0 )              fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) !=0 )              fBitmap->SetBitNumber(iCutBit,kTRUE);

  ++iCutBit;
  if ( (!fProdVtxRange2D && fCutValues->At(iCutBit) > fProdVtxXMin)
    || ( fProdVtxRange2D && prodVtxXYmin >= 1))     fBitmap->SetBitNumber(iCutBit,kTRUE);

  ++iCutBit;
  if ( (!fProdVtxRange2D && fCutValues->At(iCutBit) < fProdVtxXMax)
    || ( fProdVtxRange2D && prodVtxXYmax <= 1))     fBitmap->SetBitNumber(iCutBit,kTRUE);

  ++iCutBit;
  if ( (!fProdVtxRange2D && fCutValues->At(iCutBit) > fProdVtxYMin)
    || ( fProdVtxRange2D && prodVtxXYmin >= 1))     fBitmap->SetBitNumber(iCutBit,kTRUE);

  ++iCutBit;
  if ( (!fProdVtxRange2D && fCutValues->At(iCutBit) < fProdVtxYMax)
    || ( fProdVtxRange2D && prodVtxXYmax <= 1))     fBitmap->SetBitNumber(iCutBit,kTRUE);

  if ( fCutValues->At(++iCutBit) > fProdVtxZMin)    fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) < fProdVtxZMax)    fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) > fDecayVtxXMin)   fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) < fDecayVtxXMax)   fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) > fDecayVtxYMin)   fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) < fDecayVtxYMax)   fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) > fDecayVtxZMin)   fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) < fDecayVtxZMax)   fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) > fDecayLengthMin) fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) < fDecayLengthMax) fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) > fDecayRxyMin)    fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) < fDecayRxyMax)    fBitmap->SetBitNumber(iCutBit,kTRUE);
  if ( fCutValues->At(++iCutBit) != 0 )             fBitmap->SetBitNumber(iCutBit,kTRUE);
}


//__________________________________________________________________________________
void AliCFParticleGenCuts::FillHistograms(TObject* /*obj*/, Bool_t afterCuts)
{
  //
  // fill the QA histograms
  //

  for (int iCutNumber = 0; iCutNumber < kNCuts; iCutNumber++) 
    fhQA[iCutNumber][afterCuts]->Fill(fCutValues->At(iCutNumber));

  fhProdVtxXY[afterCuts]->Fill(fCutValues->At(4),fCutValues->At(5));

  // fill cut statistics and cut correlation histograms with information from the bitmap
  if (afterCuts) return;

  // Number of single cuts in this class
  UInt_t ncuts = fBitmap->GetNbits();
  for(UInt_t bit=0; bit<ncuts;bit++) {
    if (!fBitmap->TestBitNumber(bit)) {
      fhCutStatistics->Fill(bit+1);
      for (UInt_t bit2=bit; bit2<ncuts;bit2++) {
	if (!fBitmap->TestBitNumber(bit2)) 
	  fhCutCorrelation->Fill(bit+1,bit2+1);
      }
    }
  }
}

//__________________________________________________________________________________
void AliCFParticleGenCuts::AddQAHistograms(TList *qaList) {
  //
  // saves the histograms in a TList
  //

  DefineHistograms();

  qaList->Add(fhCutStatistics);
  qaList->Add(fhCutCorrelation);

  for (Int_t j=0; j<kNStepQA; j++) {
    qaList->Add(fhProdVtxXY[j]);
    for(Int_t i=0; i<kNCuts; i++)
      qaList->Add(fhQA[i][j]);
  }
}

//__________________________________________________________________________________
void AliCFParticleGenCuts::DefineHistograms() {
  //
  // histograms for cut variables, cut statistics and cut correlations
  //
  Int_t color = 2;

  // book cut statistics and cut correlation histograms
  fhCutStatistics = new TH1F(Form("%s_cut_statistics",GetName()),"",kNCuts,0.5,kNCuts+0.5);
  fhCutStatistics->SetLineWidth(2);
  int k = 1;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"charge")     ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"prim/sec")   ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"PDG")        ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"VtxXMin")    ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"VtxXMax")    ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"VtxYMin")    ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"VtxYMax")    ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"VtxZMin")    ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"DecZMax")    ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"DecXMin")    ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"DecXMax")    ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"DecYMin")    ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"DecYMax")    ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"DecZMin")    ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"DecZMax")    ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"DecLgthMin") ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"DecLgthMax") ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"DecRxyMin")  ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"DecRxyMax")  ; k++;
  fhCutStatistics->GetXaxis()->SetBinLabel(k,"DecChannel") ; k++;


  fhCutCorrelation = new TH2F(Form("%s_cut_correlation",GetName()),"",kNCuts,0.5,kNCuts+0.5,kNCuts,0.5,kNCuts+0.5);
  fhCutCorrelation->SetLineWidth(2);
  for (k=1; k<=kNCuts; k++) {
    fhCutCorrelation->GetXaxis()->SetBinLabel(k,fhCutStatistics->GetXaxis()->GetBinLabel(k));
    fhCutCorrelation->GetYaxis()->SetBinLabel(k,fhCutStatistics->GetXaxis()->GetBinLabel(k));
  }

  Char_t str[5];
  for (int i=0; i<kNStepQA; i++) {
    if (i==0) snprintf(str,5," ");
    else snprintf(str,5,"_cut");
    fhQA[kCutCharge]      [i] = new TH1F(Form("%s_charge%s"      ,GetName(),str),"",2,0,2);
    fhQA[kCutPrimSec]     [i] = new TH1F(Form("%s_primSec%s"     ,GetName(),str),"",2,0,2);
    fhQA[kCutPDGCode]     [i] = new TH1F(Form("%s_pdgCode%s"     ,GetName(),str),"",2,0,2);
    fhQA[kCutProdVtxXMin] [i] = new TH1F(Form("%s_prodVtxXMin%s" ,GetName(),str),"",100,-10,10);
    fhQA[kCutProdVtxXMax] [i] = new TH1F(Form("%s_prodVtxXMax%s" ,GetName(),str),"",100,-10,10);
    fhQA[kCutProdVtxYMin] [i] = new TH1F(Form("%s_prodVtxYMin%s" ,GetName(),str),"",100,-10,10);
    fhQA[kCutProdVtxYMax] [i] = new TH1F(Form("%s_prodVtxYMax%s" ,GetName(),str),"",100,-10,10);
    fhQA[kCutProdVtxZMin] [i] = new TH1F(Form("%s_prodVtxZMin%s" ,GetName(),str),"",100,-10,10);
    fhQA[kCutProdVtxZMax] [i] = new TH1F(Form("%s_prodVtxZMax%s" ,GetName(),str),"",100,-10,10);
    fhQA[kCutDecVtxXMin]  [i] = new TH1F(Form("%s_decVtxXMin%s"  ,GetName(),str),"",100,0,10);
    fhQA[kCutDecVtxXMax]  [i] = new TH1F(Form("%s_decVtxXMax%s"  ,GetName(),str),"",100,0,10);
    fhQA[kCutDecVtxYMin]  [i] = new TH1F(Form("%s_decVtxYMin%s"  ,GetName(),str),"",100,0,10);
    fhQA[kCutDecVtxYMax]  [i] = new TH1F(Form("%s_decVtxYMax%s"  ,GetName(),str),"",100,0,10);
    fhQA[kCutDecVtxZMin]  [i] = new TH1F(Form("%s_decVtxZMin%s"  ,GetName(),str),"",100,0,10);
    fhQA[kCutDecVtxZMax]  [i] = new TH1F(Form("%s_decVtxZMax%s"  ,GetName(),str),"",100,0,10);
    fhQA[kCutDecLgthMin]  [i] = new TH1F(Form("%s_decLengthMin%s",GetName(),str),"",100,0,10);
    fhQA[kCutDecLgthMax]  [i] = new TH1F(Form("%s_decLengthMax%s",GetName(),str),"",100,0,10);
    fhQA[kCutDecRxyMin]   [i] = new TH1F(Form("%s_decRxyMin%s"   ,GetName(),str),"",100,0,10);
    fhQA[kCutDecRxyMax]   [i] = new TH1F(Form("%s_decRxyMax%s"   ,GetName(),str),"",100,0,10);
    fhQA[kCutDecayChannel][i] = new TH1F(Form("%s_decayChannel%s",GetName(),str),"",2,0,2);
    fhProdVtxXY		  [i] = new TH2F(Form("%s_prodVtxXY%s"   ,GetName(),str),"",100,0,10,100,0,10);
    fhProdVtxXY		  [i] ->GetXaxis()->SetTitle("x_{production vertex}");
    fhProdVtxXY		  [i] ->GetYaxis()->SetTitle("y_{production vertex}");
    fhQA[kCutProdVtxXMax] [i] ->GetXaxis()->SetTitle("x_{production vertex}");
    fhQA[kCutProdVtxYMax] [i] ->GetXaxis()->SetTitle("y_{production vertex}");
  }
  for(Int_t i=0; i<kNCuts; i++) fhQA[i][1]->SetLineColor(color);
}


//______________________________
Bool_t AliCFParticleGenCuts::IsCharged(AliVParticle *mcPart) {
  //
  //check if particle is charged.
  //
  if (mcPart->Charge()==0) return kFALSE;
  return kTRUE;
}
//______________________________
Bool_t AliCFParticleGenCuts::IsPrimary(AliMCParticle *mcPart) {
  //
  //check if particle is primary (standard definition)
  //
  
  AliStack* stack = ((AliMCEvent*)fMCInfo)->Stack();

  if (!stack->IsPhysicalPrimary(mcPart->GetLabel())) return kFALSE;
  return kTRUE;
}
//______________________________
Bool_t AliCFParticleGenCuts::IsPrimary(AliAODMCParticle *mcPart) {
  //
  //check if particle is primary (standard definition)
  //
  
  if (!mcPart->IsPhysicalPrimary()) return kFALSE;
  return kTRUE;
}
//______________________________
Bool_t AliCFParticleGenCuts::IsPrimaryCharged(AliVParticle *mcPart) {
  //
  //check if a charged particle is primary (standard definition)
  //

  if (!fIsAODMC) {
    if (!IsPrimary((AliMCParticle*)mcPart) || !IsCharged(mcPart)) return kFALSE ;
  }
  else {
    if (!IsPrimary((AliAODMCParticle*)mcPart) || !IsCharged(mcPart)) return kFALSE ;
  }
  return kTRUE;
}
//______________________________
Bool_t AliCFParticleGenCuts::IsA(AliMCParticle *mcPart, Int_t pdg, Bool_t abs) {
  //
  //Check on the pdg code of the MC particle. if abs=kTRUE then check on the 
  //absolute value. 
  //
  TParticle* part = mcPart->Particle();
  Int_t pdgCode = part->GetPdgCode();

  if (abs) {
    pdgCode = TMath::Abs(pdgCode);
    pdg = TMath::Abs(pdg);
  }
  if (pdgCode != pdg ) return kFALSE;
  return kTRUE;
}
//______________________________
Bool_t AliCFParticleGenCuts::IsA(AliAODMCParticle *mcPart, Int_t pdg, Bool_t abs) {
  //
  //Check on the pdg code of the MC particle. if abs=kTRUE then check on the 
  //absolute value. 
  //
  Int_t pdgCode = mcPart->GetPdgCode();
  
  if (abs) {
    pdgCode = TMath::Abs(pdgCode);
    pdg = TMath::Abs(pdg);
  }
  if (pdgCode != pdg ) return kFALSE;
  return kTRUE;
}
//______________________________
void AliCFParticleGenCuts::SetMCEventInfo(const TObject* mcEvent) {
  //
  // Sets pointer to MC event information (AliMCEvent)
  //

  if (!mcEvent) {
    AliError("Pointer to MC Event is null !");
    return;
  }
  
  TString className(mcEvent->ClassName());
  if (className.CompareTo("AliMCEvent") != 0 && className.CompareTo("AliAODEvent") != 0) {
    AliError("argument must point to an AliMCEvent or an AliAODEvent !");
    return ;
  }

  fMCInfo = (AliVEvent*)mcEvent ;
}
