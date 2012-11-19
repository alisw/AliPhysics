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
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliAODv0.h"
#include "AliESDv0.h"

using std::cout;
using std::endl;

ClassImp(AliRDHFCutsLctoV0)

//--------------------------------------------------------------------------
  AliRDHFCutsLctoV0::AliRDHFCutsLctoV0(const char* name, Short_t /*v0channel*/) :
  AliRDHFCuts(name),
  fPidSelectionFlag(0),
  fPidHFV0pos(0),
  fPidHFV0neg(0),
  fV0daughtersCuts(0),
  fV0Type(0)
{
  //
  // Default Constructor
  //

  const Int_t nvars=10;
  SetNVars(nvars);
  TString varNames[nvars]={"inv. mass if K0S [GeV/c2]",
		       "inv. mass if Lambda [GeV/c2]",
		       "inv. mass V0 if K0S [GeV/c2]",
		       "inv. mass V0 if Lambda [GeV/c2]",
		       "pT min bachelor track [GeV/c]",
		       "pT min V0-positive track [GeV/c]",
		       "pT min V0-negative track [GeV/c]",
		       "dca cascade cut [cm]",
		       "dca V0 cut [cm]",
               "V0 type"
                        };

  Bool_t isUpperCut[nvars]={kTRUE,
			kTRUE,
			kTRUE,
			kTRUE,
			kFALSE,
			kFALSE,
			kFALSE,
			kTRUE,
			kTRUE,
            kTRUE
                        };
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[nvars]={kFALSE,
		    kFALSE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
            kTRUE
                    };
  SetVarsForOpt(nvars,forOpt); // It was 5: why only 5? And which ones?

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
  fPidSelectionFlag(0),
  fPidHFV0pos(0),
  fPidHFV0neg(0),
  fV0daughtersCuts(0),
  fV0Type(0)
    /*
		       fV0channel(source.fV0channel)*/
{
  //
  // Copy constructor
  //

  if (source.fPidHFV0pos) fPidHFV0pos = new AliAODPidHF(*(source.fPidHFV0pos));
  else fPidHFV0pos = new AliAODPidHF();
  if (source.fPidHFV0neg) fPidHFV0neg = new AliAODPidHF(*(source.fPidHFV0neg));
  else fPidHFV0neg = new AliAODPidHF();

  if (source.fV0daughtersCuts) fV0daughtersCuts = new AliESDtrackCuts(*(source.fV0daughtersCuts));
  else fV0daughtersCuts = new AliESDtrackCuts();

}
//--------------------------------------------------------------------------
AliRDHFCutsLctoV0 &AliRDHFCutsLctoV0::operator=(const AliRDHFCutsLctoV0 &source)
{
  //
  // assignment operator
  //

  if (this != &source) {

    AliRDHFCuts::operator=(source);
    fPidSelectionFlag = source.fPidSelectionFlag;
    delete fPidHFV0pos;
    fPidHFV0pos = new AliAODPidHF(*(source.fPidHFV0pos));
    delete fPidHFV0neg;
    fPidHFV0neg = new AliAODPidHF(*(source.fPidHFV0neg));

    delete fV0daughtersCuts;
    fV0daughtersCuts = new AliESDtrackCuts(*(source.fV0daughtersCuts));

    fV0Type  = source.fV0Type;

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

 if (fV0daughtersCuts) {
  delete fV0daughtersCuts;
  fV0daughtersCuts=0;
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
    AliDebug(2,"AliAODRecoCascadeHF null");
    return 0;
  }

  // Get the bachelor track
  AliAODTrack *bachelorTrack = (AliAODTrack*)d->GetBachelor();
  if (!bachelorTrack) {
    AliDebug(2,"No bachelor object");
    return 0;
  }

  // not used
  //if ( fUseTrackSelectionWithFilterBits &&
  //!(bachelorTrack->TestFilterMask(BIT(4))) ) return 0;

  // Get V0
  AliAODv0 *v0 = (AliAODv0*)d->Getv0();
  if (!v0) {
    AliDebug(2,"No v0 object");
    return 0;
  }



  // Get the V0 daughter tracks
  AliAODTrack *v0positiveTrack = (AliAODTrack*)d->Getv0PositiveTrack();
  AliAODTrack *v0negativeTrack = (AliAODTrack*)d->Getv0NegativeTrack();
  if (!v0positiveTrack || !v0negativeTrack ) {
    AliDebug(2,"No V0 daughters' objects");
    return 0;
  }

  // selection on daughter tracks
  if (selectionLevel==AliRDHFCuts::kAll ||
      selectionLevel==AliRDHFCuts::kTracks) {

    if (fIsCandTrackSPDFirst && d->Pt()<fMaxPtCandTrackSPDFirst) {
      if (!bachelorTrack->HasPointOnITSLayer(0)) return 0;
    }
    if (fTrackCuts && fV0daughtersCuts) {
      AliAODVertex *vAOD = (AliAODVertex*)d->GetPrimaryVtx();
      Double_t pos[3]; vAOD->GetXYZ(pos);
      Double_t cov[6]; vAOD->GetCovarianceMatrix(cov);
      const AliESDVertex vESD(pos,cov,100.,100);
      if ( !(IsDaughterSelected(bachelorTrack,&vESD,fTrackCuts)) ||
	   !(IsDaughterSelected(v0negativeTrack,&vESD,fV0daughtersCuts)) ||
	   !(IsDaughterSelected(v0positiveTrack,&vESD,fV0daughtersCuts)) ) return 0;
    }
    //if (!AreDaughtersSelected(d)) return 0;
  }

  Bool_t okLck0sp=kTRUE, okLcLpi=kTRUE, okLcLBarpi=kTRUE;
  Bool_t okK0spipi=kTRUE, okLppi=kTRUE, okLBarpip=kTRUE;

  // selection on candidate
  if (selectionLevel==AliRDHFCuts::kAll ||
      selectionLevel==AliRDHFCuts::kCandidate) {

    Double_t pt = d->Pt();
    Int_t ptbin = PtBin(pt);

    Double_t mLcPDG  = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
    Double_t mk0sPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();
    Double_t mLPDG   = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

    // K0S + p
    Double_t mk0s    = v0->MassK0Short();
    Double_t mLck0sp = d->InvMassLctoK0sP();

    // Lambda + pi
    Double_t mlambda  = v0->MassLambda();
    Double_t malambda = v0->MassAntiLambda();
    Double_t mLcLpi   = d->InvMassLctoLambdaPi();

    // cut on Lc mass with K0S+p hypothesis
    if (TMath::Abs(mLck0sp-mLcPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) {
      okLck0sp = kFALSE;
      AliDebug(4,Form(" cascade mass is %2.2e and does not correspond to Lambda_c into K0S+p cut",mLck0sp));
    }

    // cuts on the V0 mass: K0S case
    if (TMath::Abs(mk0s-mk0sPDG)>fCutsRD[GetGlobalIndex(2,ptbin)]) {
      okK0spipi = kFALSE;
      AliDebug(4,Form(" V0 mass is %2.2e and does not correspond to K0S cut",mk0s));
    }

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

  Int_t returnvaluePID = 7;

  // selection on candidate
  if (selectionLevel==AliRDHFCuts::kAll ||
      //selectionLevel==AliRDHFCuts::kCandidate ||
      selectionLevel==AliRDHFCuts::kPID )
    returnvaluePID = IsSelectedPID(d);

  //if (fUsePID && returnvaluePID==0) return 0;

  Int_t returnvalueTot = 0;
  if ( fUsePID )
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

  Bool_t okLcK0Sp = kTRUE; // K0S case
  Bool_t okLcLambdaBarPi = kTRUE; // LambdaBar case
  Bool_t okLcLambdaPi = kTRUE; // Lambda case

  CheckPID(bachelor,v0Neg,v0Pos,okLcK0Sp,okLcLambdaBarPi,okLcLambdaPi);

  Int_t returnvalue = okLcK0Sp+2*okLcLambdaBarPi+4*okLcLambdaPi;

  return returnvalue;
}
//-----------------------
void AliRDHFCutsLctoV0::CheckPID(AliAODTrack *bachelor, AliAODTrack *v0Neg, AliAODTrack *v0Pos,
				 Bool_t &isBachelorID1, Bool_t &isV0NegID2, Bool_t &isV0PosID4) {
  // identification strategy

  Int_t trackIDtof = -1;
  Int_t trackIDtpc = -1;

  switch (fPidSelectionFlag) {

  case 0:

    // identify bachelor
    trackIDtof = fPidHF->ApplyPidTOFRaw(bachelor,4);
    trackIDtpc = fPidHF->ApplyPidTPCRaw(bachelor,4);
    AliDebug(1,Form(" fPidHF->ApplyPidTOFRaw(bachelor,4)=%d fPidHF->ApplyPidTPCRaw(bachelor,4)=%d",trackIDtof,trackIDtpc));
    isBachelorID1 = (trackIDtof==4) && (trackIDtpc==4); // K0S case
    //Bool_t isBachelorID2 = (fPidHF->ApplyPidTPCRaw(bachelor,2)==2) && (fPidHF->ApplyPidTOFRaw(bachelor,2)==2); // LambdaBar case
    //Bool_t isBachelorID4 = isBachelorID2; // Lambda case

    // identify V0neg
    //Bool_t isV0NegID1 = (fPidHFV0neg->ApplyPidTPCRaw(v0Neg,2)==2) && (fPidHFV0neg->ApplyPidTOFRaw(v0Neg,2)==2); // K0S case
    trackIDtof = fPidHFV0neg->ApplyPidTOFRaw(v0Neg,4);
    trackIDtpc = fPidHFV0neg->ApplyPidTPCRaw(v0Neg,4);
    AliDebug(1,Form(" fPidHFV0neg->ApplyPidTOFRaw(v0Neg,4)=%d fPidHFV0neg->ApplyPidTPCRaw(v0Neg,4)=%d",trackIDtof,trackIDtpc));
    isV0NegID2 = (trackIDtof==4) && (trackIDtpc==4); // LambdaBar case
    //Bool_t isV0NegID4 = isV0NegID1; // Lambda case

    // identify V0pos
    //Bool_t isV0PosID1 = (fPidHFV0pos->ApplyPidTPCRaw(v0Pos,2)==2) && (fPidHFV0pos->ApplyPidTOFRaw(v0Pos,2)==2); // K0S case
    //Bool_t isV0PosID2 = isV0PosID1; // LambdaBar case
    trackIDtof = fPidHFV0pos->ApplyPidTOFRaw(v0Pos,4);
    trackIDtpc = fPidHFV0pos->ApplyPidTPCRaw(v0Pos,4);
    AliDebug(1,Form(" fPidHFV0pos->ApplyPidTOFRaw(v0Pos,4)=%d fPidHFV0pos->ApplyPidTPCRaw(v0POS,4)=%d",trackIDtof,trackIDtpc));
    isV0PosID4 = (trackIDtof==4) && (trackIDtpc==4); // Lambda case

    break;
  case 1:

    // identify bachelor
    trackIDtof = fPidHF->ApplyPidTOFRaw(bachelor,4);
    trackIDtpc = fPidHF->ApplyPidTPCRaw(bachelor,4);
    AliDebug(1,Form(" fPidHF->ApplyPidTOFRaw(bachelor,4)=%d fPidHFV0->ApplyPidTPCRaw(bachelor,4)=%d",trackIDtof,trackIDtpc));
    isBachelorID1 = ( trackIDtof==4 );
    Bool_t dummy1 = ( !(fPidHF->CheckTOFPIDStatus(bachelor)) && (trackIDtpc==4) &&
		      fPidHF->IsExcluded(bachelor,2,2.,"TPC") && fPidHF->IsExcluded(bachelor,3,2.,"TPC") ); // K0S case
    isBachelorID1 = isBachelorID1 || dummy1;


    // identify V0neg
    trackIDtof = fPidHFV0neg->ApplyPidTOFRaw(v0Neg,4);
    trackIDtpc = fPidHFV0neg->ApplyPidTPCRaw(v0Neg,4);
    AliDebug(1,Form(" fPidHFV0neg->ApplyPidTOFRaw(v0Neg,4)=%d fPidHFV0neg->ApplyPidTPCRaw(v0Neg,4)=%d",trackIDtof,trackIDtpc));
    isV0NegID2    = ( trackIDtof==4 );
    Bool_t dummy2 = ( !(fPidHFV0neg->CheckTOFPIDStatus(v0Neg)) && (trackIDtpc==4) &&
		      fPidHFV0neg->IsExcluded(v0Neg,2,2.,"TPC") && fPidHFV0neg->IsExcluded(v0Neg,3,2.,"TPC") ); // LambdaBar case
    isV0NegID2 = isV0NegID2 || dummy2;


    // identify V0pos
    trackIDtof = fPidHFV0pos->ApplyPidTOFRaw(v0Pos,4);
    trackIDtpc = fPidHFV0pos->ApplyPidTPCRaw(v0Pos,4);
    AliDebug(1,Form(" fPidHFV0pos->ApplyPidTOFRaw(v0Pos,4)=%d fPidHFV0pos->ApplyPidTPCRaw(v0Pos,4)=%d",trackIDtof,trackIDtpc));
    isV0PosID4    = ( trackIDtof==4 );
    Bool_t dummy4 = ( !(fPidHFV0pos->CheckTOFPIDStatus(v0Pos)) && (trackIDtpc==4) &&
		      fPidHFV0pos->IsExcluded(v0Pos,2,2.,"TPC") && fPidHFV0pos->IsExcluded(v0Pos,3,2.,"TPC") ); // Lambda case
    isV0PosID4 = isV0PosID4 || dummy4;


    break;

  }

}
//----------------
Int_t AliRDHFCutsLctoV0::CombinePIDCuts(Int_t returnvalue, Int_t returnvaluePID) const {
  // combine PID with topological cuts

 Int_t returnvalueTot=returnvalue&returnvaluePID;

 return returnvalueTot;
}

//----------------------------------
Int_t AliRDHFCutsLctoV0::IsSelectedSingleCut(TObject* obj, Int_t selectionLevel, Int_t cutIndex) {
  //
  // Apply selection on single cut
  //

  if (!fCutsRD) {
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }

  AliAODRecoCascadeHF* d = (AliAODRecoCascadeHF*)obj;
  if (!d) {
    AliDebug(2,"AliAODRecoCascadeHF null");
    return 0;
  }

  // Get the v0 and all daughter tracks
  AliAODTrack *bachelorTrack = (AliAODTrack*)d->GetBachelor();
  if (!bachelorTrack) {
    AliDebug(2,"No bachelor object");
    return 0;
  }

  // not used
  //if ( fUseTrackSelectionWithFilterBits &&
  //!(bachelorTrack->TestFilterMask(BIT(4))) ) return 0;

  AliAODv0 *v0 = (AliAODv0*)d->Getv0();
  if (!v0) {
    AliDebug(2,"No v0 object");
    return 0;
  }

  // Get the V0 daughter tracks
  AliAODTrack *v0positiveTrack = (AliAODTrack*)d->Getv0PositiveTrack();
  AliAODTrack *v0negativeTrack = (AliAODTrack*)d->Getv0NegativeTrack();
  if (!v0positiveTrack || !v0negativeTrack ) {
    AliDebug(2,"No V0 daughters' objects");
    return 0;
  }

  // selection on daughter tracks
  if (selectionLevel==AliRDHFCuts::kAll ||
      selectionLevel==AliRDHFCuts::kTracks) {

    if (fIsCandTrackSPDFirst && d->Pt()<fMaxPtCandTrackSPDFirst) {
      if (!bachelorTrack->HasPointOnITSLayer(0)) return 0;
    }
    if (fTrackCuts && fV0daughtersCuts) {
      AliAODVertex *vAOD = (AliAODVertex*)d->GetPrimaryVtx();
      Double_t pos[3]; vAOD->GetXYZ(pos);
      Double_t cov[6]; vAOD->GetCovarianceMatrix(cov);
      const AliESDVertex vESD(pos,cov,100.,100);
      if ( !(IsDaughterSelected(bachelorTrack,&vESD,fTrackCuts)) ||
	   !(IsDaughterSelected(v0negativeTrack,&vESD,fV0daughtersCuts)) ||
	   !(IsDaughterSelected(v0positiveTrack,&vESD,fV0daughtersCuts)) ) return 0;
    }
    //if (!AreDaughtersSelected(d)) return 0;
  }

  Bool_t okLck0sp=kFALSE, okLcLpi=kFALSE, okLcLBarpi=kFALSE;

  // selection on candidate
  if (selectionLevel==AliRDHFCuts::kAll ||
      selectionLevel==AliRDHFCuts::kCandidate) {

    Double_t pt = d->Pt();
    Int_t ptbin = PtBin(pt);

    Double_t mLcPDG  = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
    Double_t mk0sPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();
    Double_t mLPDG   = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

    // K0S + p
    Double_t mk0s    = v0->MassK0Short();
    Double_t mLck0sp = d->InvMassLctoK0sP();

    // Lambda + pi
    Double_t mlambda  = v0->MassLambda();
    Double_t malambda = v0->MassAntiLambda();
    Double_t mLcLpi   = d->InvMassLctoLambdaPi();

    switch (cutIndex) {
    case 0:
      // cut on Lc mass with K0S+p hypothesis
      okLck0sp   = TMath::Abs(mLck0sp-mLcPDG)<=fCutsRD[GetGlobalIndex(0,ptbin)];
      okLcLpi    = kFALSE;
      okLcLBarpi = kFALSE;
      break;
    case 1:
      // cut on Lc mass with Lambda+pi hypothesis
      okLck0sp   = kFALSE;
      okLcLpi    = TMath::Abs(mLcLpi-mLcPDG)<=fCutsRD[GetGlobalIndex(1,ptbin)];
      okLcLBarpi = okLcLpi;
      break;
    case 2:
      // cuts on the V0 mass: K0S case
      okLck0sp   = TMath::Abs(mk0s-mk0sPDG)<=fCutsRD[GetGlobalIndex(2,ptbin)];
      okLcLpi    = kFALSE;
      okLcLBarpi = kFALSE;
      break;
    case 3:
      // cuts on the V0 mass: Lambda/LambdaBar case
      okLck0sp   = kFALSE;
      okLcLpi    = TMath::Abs(mlambda-mLPDG)<=fCutsRD[GetGlobalIndex(3,ptbin)];
      okLcLBarpi = TMath::Abs(malambda-mLPDG)<=fCutsRD[GetGlobalIndex(3,ptbin)];
      break;
    case 4:
      // cuts on the minimum pt of bachelor
      okLck0sp   = TMath::Abs(bachelorTrack->Pt())>=fCutsRD[GetGlobalIndex(4,ptbin)];
      okLcLpi    = okLck0sp;
      okLcLBarpi = okLck0sp;
      break;
    case 5:
      // cuts on the minimum pt of positive V0daughter
      okLck0sp   = TMath::Abs(v0positiveTrack->Pt())>=fCutsRD[GetGlobalIndex(5,ptbin)];
      okLcLpi    = okLck0sp;
      okLcLBarpi = okLck0sp;
      break;
    case 6:
      // cuts on the minimum pt of negative V0daughter
      okLck0sp   = TMath::Abs(v0negativeTrack->Pt())>=fCutsRD[GetGlobalIndex(6,ptbin)];
      okLcLpi    = okLck0sp;
      okLcLBarpi = okLck0sp;
      break;
    case 7:
      // cut on cascade dca
      //okLck0sp = TMath::Abs(d->GetDCA(0))<=fCutsRD[GetGlobalIndex(8,ptbin)];
      okLck0sp   = TMath::Abs(d->GetDCA())<=fCutsRD[GetGlobalIndex(7,ptbin)] /*&&
	TMath::Abs(d->Getd0Prong(0))<=fCutsRD[GetGlobalIndex(7,ptbin)] &&
	TMath::Abs(d->Getd0Prong(1))<=fCutsRD[GetGlobalIndex(7,ptbin)]*/;
      okLcLpi    = okLck0sp;
      okLcLBarpi = okLck0sp;
      break;
    case 8:
      // cut on V0 dca
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


  /*
  Int_t returnvaluePID = 7;

  // selection on PID
  if (selectionLevel==AliRDHFCuts::kAll ||
      //selectionLevel==AliRDHFCuts::kCandidate ||
      selectionLevel==AliRDHFCuts::kPID )
    returnvaluePID = IsSelectedPID(d);
  */

  Int_t returnvalueTot = 0;
  //if ( fUsePID )
  //returnvalueTot = CombinePIDCuts(returnvalue,returnvaluePID);
  //else
  returnvalueTot = returnvalue;

  return returnvalueTot;

}
//----------------------------------
void AliRDHFCutsLctoV0::SetStandardCutsPP2010() {

 SetName("LctoV0ProductionCuts");
 SetTitle("Production cuts for Lc->V0+bachelor analysis");

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMinNClustersITS(0);//(4); // default is 5
  esdTrackCuts->SetMinNClustersTPC(70);
  //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
  //				           AliESDtrackCuts::kAny);
  // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3,1.e10);
  //esdTrackCuts->SetEtaRange(-0.8,+0.8);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  AddTrackCuts(esdTrackCuts);


  AliESDtrackCuts* esdTrackCutsV0daughters=new AliESDtrackCuts();
  esdTrackCutsV0daughters->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCutsV0daughters->SetRequireTPCRefit(kTRUE);
  esdTrackCutsV0daughters->SetRequireITSRefit(kFALSE);//(kTRUE);
  esdTrackCutsV0daughters->SetMinNClustersITS(0);//(4); // default is 5
  esdTrackCutsV0daughters->SetMinNClustersTPC(70);
  //esdTrackCutsV0daughters->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
  //					              AliESDtrackCuts::kAny);
  // default is kBoth, otherwise kAny
  esdTrackCutsV0daughters->SetMinDCAToVertexXY(0.);
  esdTrackCutsV0daughters->SetPtRange(0.,1.e10);
  esdTrackCutsV0daughters->SetAcceptKinkDaughters(kFALSE);
  AddTrackCutsV0daughters(esdTrackCutsV0daughters);

  const Int_t nptbins=1;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=99999999.;

  SetPtBins(nptbins+1,ptbins);
  SetPtBins(nptbins+1,ptbins);

  const Int_t nvars=9 ;

  Float_t** prodcutsval;
  prodcutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){prodcutsval[ic]=new Float_t[nptbins];}
  for(Int_t ipt2=0;ipt2<nptbins;ipt2++){
   prodcutsval[0][ipt2]=1.;    // inv. mass if K0S [GeV/c2]
   prodcutsval[1][ipt2]=1.;    // inv. mass if Lambda [GeV/c2]
   prodcutsval[2][ipt2]=0.05;  // inv. mass V0 if K0S [GeV/c2]
   prodcutsval[3][ipt2]=0.05;  // inv. mass V0 if Lambda [GeV/c2]
   prodcutsval[4][ipt2]=0.3;   // pT min bachelor track [GeV/c] // AOD by construction
   prodcutsval[5][ipt2]=0.;    // pT min V0-positive track [GeV/c]
   prodcutsval[6][ipt2]=0.;    // pT min V0-negative track [GeV/c]
   prodcutsval[7][ipt2]=1000.; // dca cascade cut [cm]
   prodcutsval[8][ipt2]=1000.; // dca V0 cut [nSigma] // it's 1.5 x offline V0s
  }
  SetCuts(nvars,nptbins,prodcutsval);

  SetGlobalIndex(nvars,nptbins);
  SetPtBins(nptbins+1,ptbins);


  //pid settings
  //1. bachelor: default one
  AliAODPidHF* pidObjBachelor = new AliAODPidHF();
  Double_t sigmasBac[5]={3.,1.,1.,3.,3.}; // 0, 1(A), 2(A) -> TPC; 3 -> TOF; 4 -> ITS
  pidObjBachelor->SetSigma(sigmasBac);
  pidObjBachelor->SetAsym(kFALSE);
  pidObjBachelor->SetMatch(1);
  pidObjBachelor->SetTPC(kTRUE);
  pidObjBachelor->SetTOF(kTRUE);
  pidObjBachelor->SetTOFdecide(kFALSE);
  SetPidHF(pidObjBachelor);

  //2. V0pos
  AliAODPidHF* pidObjV0pos = new AliAODPidHF();
  Double_t sigmasV0pos[5]={3.,1.,1.,3.,3.}; // 0, 1(A), 2(A) -> TPC; 3 -> TOF; 4 -> ITS
  pidObjV0pos->SetSigma(sigmasV0pos);
  pidObjV0pos->SetAsym(kFALSE);
  pidObjV0pos->SetMatch(1);
  pidObjV0pos->SetTPC(kTRUE);
  pidObjV0pos->SetTOF(kTRUE);
  pidObjV0pos->SetTOFdecide(kFALSE);
  SetPidV0pos(pidObjV0pos);

  //2. V0neg
  AliAODPidHF* pidObjV0neg = new AliAODPidHF();
  Double_t sigmasV0neg[5]={3.,1.,1.,3.,3.}; // 0, 1(A), 2(A) -> TPC; 3 -> TOF; 4 -> ITS
  pidObjV0neg->SetSigma(sigmasV0neg);
  pidObjV0neg->SetAsym(kFALSE);
  pidObjV0neg->SetMatch(1);
  pidObjV0neg->SetTPC(kTRUE);
  pidObjV0neg->SetTOF(kTRUE);
  pidObjV0neg->SetTOFdecide(kFALSE);
  SetPidV0neg(pidObjV0neg);

  SetUsePID(kFALSE);//(kTRUE);

  PrintAll();

 for(Int_t iiv=0;iiv<nvars;iiv++){
  delete [] prodcutsval[iiv];
 }
 delete [] prodcutsval;
 prodcutsval=NULL;
 delete [] ptbins;
 ptbins=NULL;


 delete pidObjBachelor;
 pidObjBachelor=NULL;

 delete pidObjV0pos;
 pidObjV0pos=NULL;

 delete pidObjV0neg;
 pidObjV0neg=NULL;

 return;
}
//------------------
void AliRDHFCutsLctoV0::SetStandardCutsPbPb2010() {

 SetName("LctoV0ProductionCuts");
 SetTitle("Production cuts for Lc->V0+bachelor analysis");

 SetStandardCutsPP2010();

 return;
}
//------------------
void AliRDHFCutsLctoV0::SetStandardCutsPbPb2011() {

  // Default 2010 PbPb cut object
  SetStandardCutsPbPb2010();

  //
  // Enable all 2011 PbPb run triggers
  //
  SetTriggerClass("");
  ResetMaskAndEnableMBTrigger();
  EnableCentralTrigger();
  EnableSemiCentralTrigger();

}
//-----------------------
Int_t AliRDHFCutsLctoV0::GetV0Type(){

  const Int_t nvars = this->GetNVars() ;
  //Float_t *vArray =GetCuts();
  //fV0Type = vArray[nvars-1];
  fV0Type = (this->GetCuts())[nvars-1];
  //this->GetCuts(vArray);
  TString *sVarNames=GetVarNames();

  if(sVarNames[nvars-1].Contains("V0 type")) return (Int_t)fV0Type;
  else {AliInfo("AliRDHFCutsLctoV0 Last variable is not the Vo type!!!"); return -999;}
}
