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
// Class for cuts on AOD reconstructed Lc->pKpi
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
AliRDHFCutsLctoV0::AliRDHFCutsLctoV0(const char* name) : 
AliRDHFCuts(name)
{
  //
  // Default Constructor
  //
  Int_t nvars=9;
  SetNVars(nvars);
  TString varNames[9]={"inv. mass if K0s [GeV]",
		       "inv. mass if Lambda [GeV]",
		       "inv. mass V0 if K0s [GeV]",
		       "inv. mass V0 if Lambda [GeV]",
		       "pT min bachelor track [GeV/c]",
		       "pT min V0-positive track [GeV/c]",
		       "pT min V0-negative track [GeV/c]",
		       "dca cascade cut (cm)",
		       "dca V0 cut (cm)"}; 

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
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kFALSE};
  SetVarsForOpt(5,forOpt);
  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);
}
//--------------------------------------------------------------------------
AliRDHFCutsLctoV0::AliRDHFCutsLctoV0(const AliRDHFCutsLctoV0 &source) :
  AliRDHFCuts(source)
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
  if(&source == this) return *this;

  AliRDHFCuts::operator=(source);

  return *this;
}


//---------------------------------------------------------------------------
void AliRDHFCutsLctoV0::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) {
  // 
  // Fills in vars the values of the variables 
  //
  if(pdgdaughters[0]==-9999) return; // dummy


  if(nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsLctoV0::GetCutsVarsForOpt: wrong number of variables\n");
    return;
  }

  AliAODRecoCascadeHF *dd = (AliAODRecoCascadeHF*)d;

  // Get the v0 and all daughter tracks
  AliAODTrack *bachelor_track = dd->GetBachelor();
  AliAODv0 *v0 = dd->Getv0();
  AliAODTrack *v0positive_track = dd->Getv0PositiveTrack();
  AliAODTrack *v0negative_track = dd->Getv0NegativeTrack(); 

  Int_t iter=-1;
  // cut on cascade mass, if k0s + p
  if(fVarsForOpt[0]){
    iter++;
    vars[iter]=dd->InvMassLctoK0sP();
  }
  // cut on cascade mass, if lambda + pi
  if(fVarsForOpt[1]){
    iter++;
    vars[iter]=dd->InvMassLctoLambdaPi();
  }

  // cut on v0 mass if k0s
  if(fVarsForOpt[2]){
    iter++;
    vars[iter]=v0->MassK0Short();
  }

  // cut on v0 mass if lambda
  // ----------------------------- pb with anti-lambda?? --------->>>>>>>>
  if(fVarsForOpt[3]){
    iter++;
    vars[iter]=v0->MassLambda();
  }

  // cut on v0-positive min pt
  if(fVarsForOpt[4]){
    iter++;
    vars[iter]=v0positive_track->Pt();
  }

  // cut on v0-negative min pt
  if(fVarsForOpt[5]){
    iter++;
    vars[iter]=v0negative_track->Pt();
  }

  // cut bachelor min pt
  if(fVarsForOpt[6]){
    iter++;
    vars[iter]=bachelor_track->Pt();
  }

  // cut on v0 dca
  if(fVarsForOpt[7]){
    iter++;
    vars[iter]=v0->GetDCA();
  }

  // cut on cascade dca
  if(fVarsForOpt[8]){
    iter++;
    vars[iter]=dd->GetDCA();
  }

  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctoV0::IsSelected(TObject* obj,Int_t selectionLevel) {
  //
  // Apply selection
  //

  if(!fCutsRD){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  //PrintAll();
  AliAODRecoCascadeHF* d=(AliAODRecoCascadeHF*)obj;

  if(!d){
    cout<<"AliAODRecoCascadeHF null"<<endl;
    return 0;
  }

  if(d->HasBadDaughters()) return 0;

  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(d)) return 0;
  }


  // Get the v0 and all daughter tracks
  AliAODTrack *bachelor_track = d->GetBachelor();
  AliAODv0 *v0 = d->Getv0();
  AliAODTrack *v0positive_track = d->Getv0PositiveTrack();
  AliAODTrack *v0negative_track = d->Getv0NegativeTrack(); 

//   // If reading ESDv0, return false
//   if ( !d->Getv0() || !d->Getv0PositiveTrack() || !d->Getv0NegativeTrack() )
//     { AliInfo(Form("Not adapted for ESDv0s, return false...")); return false; }

  Int_t returnvalue=1;

  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {

    Double_t pt=d->Pt();
    
    Int_t ptbin=PtBin(pt);
    
    Double_t mLck0sp,mLcLpi;
    Int_t okLck0sp=1,okLcLpi=1;

    Double_t mLcPDG = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
    Double_t mk0sPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();
    Double_t mLPDG = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

    // k0s + p
    double mk0s = v0->MassK0Short();
    mLck0sp = d->InvMassLctoK0sP();

    // lambda + pi 
    double mlambda = v0->MassLambda();
    double malambda = v0->MassAntiLambda();
    mLcLpi = d->InvMassLctoLambdaPi();

    // cut on Lc mass
    //   with k0s p hypothesis
    if(TMath::Abs(mLck0sp-mLcPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) {
      okLck0sp = 0;
      AliDebug(4,Form(" cascade mass is %2.2e and does not correspond to lambda_c into k0s+p cut",mLck0sp));
    }
    //   with Lambda pi hypothesis
    if(TMath::Abs(mLcLpi-mLcPDG)>fCutsRD[GetGlobalIndex(1,ptbin)]) {
      okLcLpi = 0;
      AliDebug(4,Form(" cascade mass is %2.2e and does not correspond to lambda_c into lambda+pi cut",mLcLpi));
    }

    // cuts on the v0 mass
    if(TMath::Abs(mk0s-mk0sPDG)>fCutsRD[GetGlobalIndex(2,ptbin)]) { 
      okLck0sp = 0;
      AliDebug(4,Form(" v0 mass is %2.2e and does not correspond to k0s cut",mk0s));
    }
    if( TMath::Abs(mlambda-mLPDG)>fCutsRD[GetGlobalIndex(3,ptbin)] && 
	TMath::Abs(malambda-mLPDG)>fCutsRD[GetGlobalIndex(3,ptbin)] ) {
      okLcLpi = 0;
      AliDebug(4,Form(" v0 mass is %2.2e and does not correspond to lambda cut",mlambda));
    }

    if(!okLck0sp && !okLcLpi) return 0;

    // cuts on the minimum pt of the tracks 
    if(TMath::Abs(bachelor_track->Pt()) < fCutsRD[GetGlobalIndex(4,ptbin)]) {
      AliDebug(4,Form(" bachelor track Pt=%2.2e > %2.2e",bachelor_track->Pt(),fCutsRD[GetGlobalIndex(4,ptbin)]));
      return 0;
    }
    if(TMath::Abs(v0positive_track->Pt()) < fCutsRD[GetGlobalIndex(5,ptbin)]) {
      AliDebug(4,Form(" v0 positive track Pt=%2.2e > %2.2e",v0positive_track->Pt(),fCutsRD[GetGlobalIndex(5,ptbin)]));
      return 0;
    }
    if(TMath::Abs(v0negative_track->Pt()) < fCutsRD[GetGlobalIndex(6,ptbin)]) {
      AliDebug(4,Form(" v0 negative track Pt=%2.2e > %2.2e",v0negative_track->Pt(),fCutsRD[GetGlobalIndex(6,ptbin)]));
      return 0;
    }

    // cut on the v0 dca
    if(TMath::Abs(v0->DcaV0Daughters()) > fCutsRD[GetGlobalIndex(7,ptbin)]) {
      AliDebug(4,Form(" v0 daughters DCA = %2.2e > %2.2e",v0->DcaV0Daughters(),fCutsRD[GetGlobalIndex(7,ptbin)]));
      return 0;
    }

    // cut on the cascade dca
    if( TMath::Abs(d->GetDCA(0))>fCutsRD[GetGlobalIndex(8,ptbin)] ||
	TMath::Abs(v0->DcaPosToPrimVertex())>fCutsRD[GetGlobalIndex(8,ptbin)] ||
	TMath::Abs(v0->DcaNegToPrimVertex())>fCutsRD[GetGlobalIndex(8,ptbin)] ) {
      AliDebug(4,Form(" cascade tracks DCA at primary vertex don't pass the cut"));
      return 0;
    }

    if(okLck0sp) returnvalue=1; //cuts passed as Lc -> k0s + p
    if(okLcLpi) returnvalue=2; //cuts passed as Lc-> lambda + pi
    if(okLck0sp && okLcLpi) returnvalue=3; //cuts passed as both  Lc -> k0s + p; Lc-> lambda + pi

  }

  return returnvalue;
}
//---------------------------------------------------------------------------
