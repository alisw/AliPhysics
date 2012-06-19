/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// Class for cuts on AOD reconstructed Lc->V0+X
//
// Modified by A.De Caro - decaro@sa.infn.it
//
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>

#include "AliRDHFCutsLctoV0.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAODv0.h"
#include "AliESDv0.h"

ClassImp(AliRDHFCutsLctoV0)

//--------------------------------------------------------------------------
  AliRDHFCutsLctoV0::AliRDHFCutsLctoV0(const char* name, Short_t /*v0channel*/) : 
  AliRDHFCuts(name),
  fPidHFV0pos(new AliAODPidHF()),
  fPidHFV0neg(new AliAODPidHF())
{
  //
  // Default Constructor
  //

  Int_t nvars=9;
  SetNVars(nvars);
  TString varNames[9]={"inv. mass if K0S [GeV/c2]",
		       "inv. mass if Lambda [GeV/c2]",
		       "inv. mass V0 if K0S [GeV/c2]",
		       "inv. mass V0 if Lambda [GeV/c2]",
		       "pT min bachelor track [GeV/c]",
		       "pT min V0-positive track [GeV/c]",
		       "pT min V0-negative track [GeV/c]",
		       "dca cascade cut [cm]",
		       "dca V0 cut [cm]"};

  Bool_t isUpperCut[9]={kTRUE,
			kTRUE, 
			kTRUE, 
			kTRUE,
			kFALSE,
			kFALSE,
			kFALSE,
			kTRUE,
			kTRUE};
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[9]={kFALSE,
		    kFALSE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
		    kTRUE};
  SetVarsForOpt(9,forOpt); // It was 5: why only 5? And which ones?

  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);

  /*
  switch (v0channel) {
  case 0:
    fV0channel = 0x0001;
    break;
  case 1:
    fV0channel = 0x0002;
    break;
  case 2:
    fV0channel = 0x0004;
    break;
  }
  */

}
//--------------------------------------------------------------------------
AliRDHFCutsLctoV0::AliRDHFCutsLctoV0(const AliRDHFCutsLctoV0 &source) :
  AliRDHFCuts(source),
  fPidHFV0pos(new AliAODPidHF(*(source.fPidHFV0pos))),
  fPidHFV0neg(new AliAODPidHF(*(source.fPidHFV0neg)))/*,
						       fV0channel(source.fV0channel)*/
{
  //
  // Copy constructor
  //

}
//--------------------------------------------------------------------------
AliRDHFCutsLctoV0 &AliRDHFCutsLctoV0::operator=(const AliRDHFCutsLctoV0 &source)
{
  //
  // assignment operator
  //

  if (this != &source) {

    AliRDHFCuts::operator=(source);
    delete fPidHFV0pos;
    fPidHFV0pos = new AliAODPidHF(*(source.fPidHFV0pos));
    delete fPidHFV0neg;
    fPidHFV0neg = new AliAODPidHF(*(source.fPidHFV0neg));

    //fV0channel = source.fV0channel;

  }

  return *this;
}


//---------------------------------------------------------------------------
AliRDHFCutsLctoV0::~AliRDHFCutsLctoV0() {
 //
 //  // Default Destructor
 //   

 if (fPidHFV0pos) {
  delete fPidHFV0pos;
  fPidHFV0pos=0;
 }
 if (fPidHFV0neg) {
  delete fPidHFV0neg;
  fPidHFV0neg=0;
 }

}

//---------------------------------------------------------------------------
void AliRDHFCutsLctoV0::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) {
  // 
  // Fills in vars the values of the variables 
  //

  if (pdgdaughters[0]==-9999) return; // dummy

  if (nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsLctoV0::GetCutsVarsForOpt: wrong number of variables\n");
    return;
  }

  AliAODRecoCascadeHF *dd = (AliAODRecoCascadeHF*)d;

  // Get the v0 and all daughter tracks
  AliAODTrack *bachelorTrack = (AliAODTrack*)dd->GetBachelor();
  AliAODv0 *v0 = (AliAODv0*)dd->Getv0();
  AliAODTrack *v0positiveTrack = (AliAODTrack*)dd->Getv0PositiveTrack();
  AliAODTrack *v0negativeTrack = (AliAODTrack*)dd->Getv0NegativeTrack(); 

  Int_t iter=-1;
  // cut on cascade mass, if K0S + p
  if (fVarsForOpt[0]) {
    iter++;
    vars[iter]=dd->InvMassLctoK0sP();
  }
  // cut on cascade mass, if Lambda + pi
  if (fVarsForOpt[1]) {
    iter++;
    vars[iter]=dd->InvMassLctoLambdaPi();
  }

  // cut on V0 mass if K0S
  if (fVarsForOpt[2]) {
    iter++;
    vars[iter]=v0->MassK0Short();
  }

  // cut on V0 mass if Lambda
  // ----------------------------- pb with anti-lambda?? --------->>>>>>>>
  if (fVarsForOpt[3]) {

    if (bachelorTrack->Charge()==1) {
      iter++;
      vars[iter]=v0->MassLambda();
    } else if (bachelorTrack->Charge()==-1) {
      iter++;
      vars[iter]=v0->MassAntiLambda();
    }

  }

  // cut bachelor min pt
  if (fVarsForOpt[4]) {
    iter++;
    vars[iter]=bachelorTrack->Pt();
  }

  // cut on V0-positive min pt
  if (fVarsForOpt[5]) {
    iter++;
    vars[iter]=v0positiveTrack->Pt();
  }

  // cut on V0-negative min pt
  if (fVarsForOpt[6]) {
    iter++;
    vars[iter]=v0negativeTrack->Pt();
  }

  // cut on cascade dca
  if (fVarsForOpt[7]) {
    iter++;
    vars[iter]=dd->GetDCA(); // prong-to-prong DCA
  }
 
  // cut on V0 dca
  if (fVarsForOpt[8]) {
    iter++;
    vars[iter]=v0->GetDCA(); // prong-to-prong DCA
  }

  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctoV0::IsSelected(TObject* obj,Int_t selectionLevel) {
  //
  // Apply selection
  //

  if (!fCutsRD) {
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }

  //PrintAll();

  AliAODRecoCascadeHF* d = (AliAODRecoCascadeHF*)obj;
  if (!d) {
    cout<<"AliAODRecoCascadeHF null"<<endl;
    return 0;
  }

  // selection on daughter tracks 
  if (selectionLevel==AliRDHFCuts::kAll || 
      selectionLevel==AliRDHFCuts::kTracks) {
    if (!AreDaughtersSelected(d)) return 0;
  }

  // Get the bachelor track
  AliAODTrack *bachelorTrack = (AliAODTrack*)d->GetBachelor();
  if (!(bachelorTrack->TestFilterMask(BIT(4)))) return 0;

  // Get the V0 and all daughter tracks
  AliAODv0 *v0 = (AliAODv0*)d->Getv0();
  AliAODTrack *v0positiveTrack = (AliAODTrack*)d->Getv0PositiveTrack();
  AliAODTrack *v0negativeTrack = (AliAODTrack*)d->Getv0NegativeTrack(); 

  // If reading ESDv0, return false
  if ( !d->Getv0() || !d->Getv0PositiveTrack() || !d->Getv0NegativeTrack() ) {
    AliInfo(Form("Not adapted for ESDv0s, return false..."));
    return 0;
  }

  Int_t returnvaluePID = 7;

  // selection on candidate
  if (selectionLevel==AliRDHFCuts::kAll || 
      selectionLevel==AliRDHFCuts::kCandidate || 
      selectionLevel==AliRDHFCuts::kPID)
    returnvaluePID = IsSelectedPID(d);

  //if (fUsePID && returnvaluePID==0) return 0;

  Bool_t okLck0sp=kTRUE, okLcLpi=kTRUE, okLcLBarpi=kTRUE;
  Bool_t okK0spipi=kTRUE, okLppi=kTRUE, okLBarpip=kTRUE;

  // selection on candidate
  if (selectionLevel==AliRDHFCuts::kAll || 
      selectionLevel==AliRDHFCuts::kCandidate) {

    Double_t pt = d->Pt();
    
    Int_t ptbin = PtBin(pt);
    
    Double_t mLcPDG = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
    Double_t mk0sPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();
    Double_t mLPDG = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

    // K0S + p
    Double_t mk0s = v0->MassK0Short();
    Double_t mLck0sp = d->InvMassLctoK0sP();
    // cut on Lc mass with K0S+p hypothesis
    if (TMath::Abs(mLck0sp-mLcPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) {
      okLck0sp = kFALSE;
      AliDebug(4,Form(" cascade mass is %2.2e and does not correspond to Lambda_c into K0S+p cut",mLck0sp));
    }

    // cuts on the V0 mass: K0S case
    if (TMath::Abs(mk0s-mk0sPDG)>fCutsRD[GetGlobalIndex(2,ptbin)]) { 
      okK0spipi = 0;
      AliDebug(4,Form(" V0 mass is %2.2e and does not correspond to K0S cut",mk0s));
    }

    // Lambda + pi 
    Double_t mlambda = v0->MassLambda();
    Double_t malambda = v0->MassAntiLambda();
    Double_t mLcLpi = d->InvMassLctoLambdaPi();
    // cut on Lc mass with Lambda+pi hypothesis
    if (TMath::Abs(mLcLpi-mLcPDG)>fCutsRD[GetGlobalIndex(1,ptbin)]) {
      okLcLpi = kFALSE;
      okLcLBarpi = kFALSE;
      AliDebug(4,Form(" cascade mass is %2.2e and does not correspond to Lambda_c into Lambda+pi cut",mLcLpi));
    }

    // cuts on the V0 mass: Lambda/LambdaBar case
    if ( TMath::Abs(mlambda-mLPDG)>fCutsRD[GetGlobalIndex(3,ptbin)] ) {
      okLppi = kFALSE;
      AliDebug(4,Form(" V0 mass is %2.2e and does not correspond to LambdaBar cut",mlambda));
    }

    if ( TMath::Abs(malambda-mLPDG)>fCutsRD[GetGlobalIndex(3,ptbin)] ) {
      okLBarpip = kFALSE;
      AliDebug(4,Form(" V0 mass is %2.2e and does not correspond to LambdaBar cut",malambda));
    }

    okLck0sp   = okLck0sp   && okK0spipi;
    okLcLpi    = okLcLpi    && okLppi;
    okLcLBarpi = okLcLBarpi && okLBarpip;

    if (!okLck0sp && !okLcLpi && !okLcLBarpi) return 0;

    // cuts on the minimum pt of the tracks 
    if (TMath::Abs(bachelorTrack->Pt()) < fCutsRD[GetGlobalIndex(4,ptbin)]) {
      AliDebug(4,Form(" bachelor track Pt=%2.2e > %2.2e",bachelorTrack->Pt(),fCutsRD[GetGlobalIndex(4,ptbin)]));
      return 0;
    }
    if (TMath::Abs(v0positiveTrack->Pt()) < fCutsRD[GetGlobalIndex(5,ptbin)]) {
      AliDebug(4,Form(" V0-positive track Pt=%2.2e > %2.2e",v0positiveTrack->Pt(),fCutsRD[GetGlobalIndex(5,ptbin)]));
      return 0;
    }
    if (TMath::Abs(v0negativeTrack->Pt()) < fCutsRD[GetGlobalIndex(6,ptbin)]) {
      AliDebug(4,Form(" V0-negative track Pt=%2.2e > %2.2e",v0negativeTrack->Pt(),fCutsRD[GetGlobalIndex(6,ptbin)]));
      return 0;
    }

    // cut on cascade dca
    //if (TMath::Abs(d->GetDCA()) > fCutsRD[GetGlobalIndex(7,ptbin)]) {
    if ( TMath::Abs(d->GetDCA()) > fCutsRD[GetGlobalIndex(7,ptbin)] /*|| // prong-to-prong DCA
	 TMath::Abs(d->Getd0Prong(0)) > fCutsRD[GetGlobalIndex(7,ptbin)] || // rphi impact params w.r.t. Primary Vtx
	 TMath::Abs(d->Getd0Prong(1)) > fCutsRD[GetGlobalIndex(7,ptbin)]*/ ) { // rphi impact params w.r.t. Primary Vtx
      AliDebug(4,Form(" cascade tracks DCA don't pass the cut"));
      return 0;
    }

    // cut on V0 dca
    //if (TMath::Abs(v0->DcaV0Daughters()) > fCutsRD[GetGlobalIndex(8,ptbin)]) { // prong-to-prong DCA
    if ( TMath::Abs(v0->GetDCA()) > fCutsRD[GetGlobalIndex(8,ptbin)] /*|| // prong-to-prong DCA
	 TMath::Abs(v0->DcaV0ToPrimVertex()) > fCutsRD[GetGlobalIndex(8,ptbin)] || // V0-to-vertex DCA
	 TMath::Abs(v0->DcaPosToPrimVertex()) > fCutsRD[GetGlobalIndex(8,ptbin)] || // prong-to-vertex DCA
	 TMath::Abs(v0->DcaNegToPrimVertex()) > fCutsRD[GetGlobalIndex(8,ptbin)]*/ ) { // prong-to-vertex DCA
      AliDebug(4,Form(" V0 DCA don't pass the cut"));
      return 0;
    }

  }

  Int_t returnvalue = okLck0sp+2*okLcLBarpi+4*okLcLpi;
  /*
    retvalue case
          1  Lc->K0S + p
          2  Lc->LambdaBar + pi
          3  Lc->K0S + p AND Lc->LambdaBar + pi
          4  Lc->Lambda + pi
          5  Lc->K0S + p AND Lc->Lambda + pi
          6  Lc->LambdaBar + pi AND Lc->Lambda + pi
          7  Lc->K0S + p AND Lc->LambdaBar + pi AND Lc->Lambda + pi
  */

  Int_t returnvalueTot = 0;
  if (fUsePID )
    returnvalueTot = CombinePIDCuts(returnvalue,returnvaluePID);
  else
    returnvalueTot = returnvalue;

  return returnvalueTot;

}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctoV0::IsSelectedPID(AliAODRecoDecayHF* obj) {

  // fPidHF -> PID object for bachelor
  // fPidHFV0pos -> PID object for positive V0 daughter
  // fPidHFV0neg -> PID object for negative V0 daughter

  if (!fUsePID || !obj) return 7; // all hypothesis are valid

  AliAODRecoCascadeHF *objD = (AliAODRecoCascadeHF*)obj;

  Bool_t isPeriodd = fPidHF->GetOnePad();
  Bool_t isMC = fPidHF->GetMC();

  if (isPeriodd) {
    fPidHFV0pos->SetOnePad(kTRUE);
    fPidHFV0neg->SetOnePad(kTRUE);
  }
  if (isMC) {
    fPidHFV0neg->SetMC(kTRUE);
    fPidHFV0pos->SetMC(kTRUE);
  }

  AliAODTrack *bachelor = (AliAODTrack*)objD->GetBachelor();
  AliAODTrack *v0Pos = (AliAODTrack*)objD->Getv0PositiveTrack();
  AliAODTrack *v0Neg = (AliAODTrack*)objD->Getv0NegativeTrack();

  if (!bachelor || !v0Pos || !v0Neg) return 0;

  // identify bachelor
  Bool_t isBachelorID1 = fPidHF->IsProtonRaw(bachelor,"TPC") && fPidHF->IsProtonRaw(bachelor,"TOF"); // K0S case
  //Bool_t isBachelorID2 = fPidHF->IsPionRaw(bachelor,"TPC") && fPidHF->IsPionRaw(bachelor,"TOF"); // LambdaBar case
  //Bool_t isBachelorID4 = isBachelorID2; // Lambda case

  // identify V0pos
  //Bool_t isV0PosID1 = fPidHFV0pos->IsPionRaw(v0Pos,"TPC") && fPidHFV0pos->IsPionRaw(v0Pos,"TOF"); // K0S case
  //Bool_t isV0PosID2 = isV0PosID1; // LambdaBar case
  Bool_t isV0PosID4 = fPidHFV0pos->IsProtonRaw(v0Pos,"TPC") && fPidHFV0pos->IsProtonRaw(v0Pos,"TOF"); // Lambda case

  // identify V0neg
  //Bool_t isV0NegID1 = fPidHFV0neg->IsPionRaw(v0Neg,"TPC") && fPidHFV0neg->IsPionRaw(v0Neg,"TOF"); // K0S case
  Bool_t isV0NegID2 = fPidHFV0neg->IsProtonRaw(v0Neg,"TPC") && fPidHFV0neg->IsProtonRaw(v0Neg,"TOF"); // LambdaBar case
  //Bool_t isV0NegID4 = isV0NegID1; // Lambda case

  Bool_t okLcK0Sp = isBachelorID1; // K0S case
  Bool_t okLcLambdaBarPi = isV0NegID2; // LambdaBar case
  Bool_t okLcLambdaPi = isV0PosID4; // Lambda case

  Int_t returnvalue = okLcK0Sp+2*okLcLambdaBarPi+4*okLcLambdaPi;

  return returnvalue;
}
//-----------------------

Int_t AliRDHFCutsLctoV0::CombinePIDCuts(Int_t returnvalue, Int_t returnvaluePID) const {
  // combine PID with topological cuts

 Int_t returnvalueTot=returnvalue&returnvaluePID;

 return returnvalueTot;
}

//----------------------------------
Int_t AliRDHFCutsLctoV0::IsSelected(TObject* obj, Int_t selectionLevel, Int_t cutIndex) {
  //
  // Apply selection
  //

  if (!fCutsRD) {
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }

  AliAODRecoCascadeHF* d = (AliAODRecoCascadeHF*)obj;
  if (!d) {
    cout<<"AliAODRecoCascadeHF null"<<endl;
    return 0;
  }

  // selection on daughter tracks 
  if (selectionLevel==AliRDHFCuts::kAll || 
      selectionLevel==AliRDHFCuts::kTracks || 
      selectionLevel==AliRDHFCuts::kCandidate) {
    if (!AreDaughtersSelected(d)) return 0;
  }

  // Get the v0 and all daughter tracks
  AliAODTrack *bachelorTrack = (AliAODTrack*)d->GetBachelor();
  if (!(bachelorTrack->TestFilterMask(BIT(4)))) return 0;

  AliAODv0 *v0 = (AliAODv0*)d->Getv0();
  AliAODTrack *v0positiveTrack = (AliAODTrack*)d->Getv0PositiveTrack();
  AliAODTrack *v0negativeTrack = (AliAODTrack*)d->Getv0NegativeTrack(); 
  if ( !d->Getv0() || !d->Getv0PositiveTrack() || !d->Getv0NegativeTrack() ) {
    AliInfo(Form("Not adapted for ESDv0s, return false..."));
    return 0;
  }


  Bool_t okLck0sp=kFALSE, okLcLpi=kFALSE, okLcLBarpi=kFALSE;

  // selection on candidate
  if (selectionLevel==AliRDHFCuts::kAll || 
      selectionLevel==AliRDHFCuts::kCandidate) {

    Double_t mLcPDG =  TDatabasePDG::Instance()->GetParticle(4122)->Mass();
    Double_t mk0sPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();
    Double_t mLPDG =   TDatabasePDG::Instance()->GetParticle(3122)->Mass();

    // k0s + p
    Double_t mk0s    = v0->MassK0Short();
    Double_t mLck0sp = d->InvMassLctoK0sP();

    // lambda + pi 
    Double_t mlambda  = v0->MassLambda();
    Double_t malambda = v0->MassAntiLambda();
    Double_t mLcLpi   = d->InvMassLctoLambdaPi();

    Double_t pt = d->Pt();
    Int_t ptbin = PtBin(pt);

    switch (cutIndex) {
    case 0:
      okLck0sp   = TMath::Abs(mLck0sp-mLcPDG)<=fCutsRD[GetGlobalIndex(0,ptbin)];
      okLcLpi    = kFALSE;
      okLcLBarpi = kFALSE;
      break;
    case 1:
      okLck0sp   = kFALSE;
      okLcLpi    = TMath::Abs(mLcLpi-mLcPDG)<=fCutsRD[GetGlobalIndex(1,ptbin)];
      okLcLBarpi = okLcLpi;
      break;
    case 2:
      okLck0sp   = TMath::Abs(mk0s-mk0sPDG)<=fCutsRD[GetGlobalIndex(2,ptbin)];
      okLcLpi    = kFALSE;
      okLcLBarpi = kFALSE;
      break;
    case 3:
      okLck0sp   = kFALSE;
      okLcLpi    = TMath::Abs(mlambda-mLPDG)<=fCutsRD[GetGlobalIndex(3,ptbin)];
      okLcLBarpi = TMath::Abs(malambda-mLPDG)<=fCutsRD[GetGlobalIndex(3,ptbin)];
      break;
    case 4:
      okLck0sp   = TMath::Abs(bachelorTrack->Pt())>=fCutsRD[GetGlobalIndex(4,ptbin)];
      okLcLpi    = okLck0sp;
      okLcLBarpi = okLck0sp;
      break;
    case 5:
      okLck0sp   = TMath::Abs(v0positiveTrack->Pt())>=fCutsRD[GetGlobalIndex(5,ptbin)];
      okLcLpi    = okLck0sp;
      okLcLBarpi = okLck0sp;
      break;
    case 6:
      okLck0sp   = TMath::Abs(v0negativeTrack->Pt())>=fCutsRD[GetGlobalIndex(6,ptbin)];
      okLcLpi    = okLck0sp;
      okLcLBarpi = okLck0sp;
      break;
    case 7:
      //okLck0sp = TMath::Abs(d->GetDCA(0))<=fCutsRD[GetGlobalIndex(8,ptbin)];
      okLck0sp   = TMath::Abs(d->GetDCA())<=fCutsRD[GetGlobalIndex(7,ptbin)] /*&&
	TMath::Abs(d->Getd0Prong(0))<=fCutsRD[GetGlobalIndex(7,ptbin)] &&
	TMath::Abs(d->Getd0Prong(1))<=fCutsRD[GetGlobalIndex(7,ptbin)]*/;
      okLcLpi    = okLck0sp;
      okLcLBarpi = okLck0sp;
      break;
    case 8:
      //okLck0sp   = TMath::Abs(v0->DcaV0Daughters())<=fCutsRD[GetGlobalIndex(7,ptbin)];
      okLck0sp   = TMath::Abs(v0->GetDCA())<=fCutsRD[GetGlobalIndex(8,ptbin)] /*&&
	TMath::Abs(v0->DcaV0ToPrimVertex())<=fCutsRD[GetGlobalIndex(8,ptbin)] &&
	TMath::Abs(v0->DcaPosToPrimVertex())<=fCutsRD[GetGlobalIndex(8,ptbin)] &&
	TMath::Abs(v0->DcaNegToPrimVertex())<=fCutsRD[GetGlobalIndex(8,ptbin)]*/;
      okLcLpi    = okLck0sp;
      okLcLBarpi = okLck0sp;
      break;
    }

  }

  Int_t returnvalue = okLck0sp+2*okLcLBarpi+4*okLcLpi;
  /*
    retvalue case
          1  Lc->K0S + p
          2  Lc->LambdaBar + pi
          3  Lc->K0S + p AND Lc->LambdaBar + pi
          4  Lc->Lambda + pi
          5  Lc->K0S + p AND Lc->Lambda + pi
          6  Lc->LambdaBar + pi AND Lc->Lambda + pi
          7  Lc->K0S + p AND Lc->LambdaBar + pi AND Lc->Lambda + pi
  */


  Int_t returnvaluePID = 7;

  // selection on PID
  if (selectionLevel==AliRDHFCuts::kAll || 
      selectionLevel==AliRDHFCuts::kPID)
    returnvaluePID = IsSelectedPID(d);

  //if (fUsePID && returnvaluePID==0) return 0;

  Int_t returnvalueTot = /*0;
  if (fUsePID)
  returnvalueTot =*/ CombinePIDCuts(returnvalue,returnvaluePID);
  /*else
    returnvalueTot = returnvalue;*/

  return returnvalueTot;

}
