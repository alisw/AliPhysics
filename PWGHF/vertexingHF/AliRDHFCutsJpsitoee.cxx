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
// Class for cuts on AOD reconstructed Jpsi->ee
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>

#include "AliRDHFCutsJpsitoee.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"

ClassImp(AliRDHFCutsJpsitoee)

//--------------------------------------------------------------------------
AliRDHFCutsJpsitoee::AliRDHFCutsJpsitoee(const char* name) : 
AliRDHFCuts(name)
{
  //
  // Default Constructor
  //
  Int_t nvars=9;
  SetNVars(nvars);
  TString varNames[9]={"inv. mass [GeV]",   
		       "dca [cm]",
		       "cosThetaStar", // negative electron 
		       "pTP [GeV/c]",
		       "pTN [GeV/c]",
		       "d0P [cm]",
		       "d0N [cm]",
		       "d0d0 [cm^2]",
		       "cosThetaPoint"};
  Bool_t isUpperCut[9]={kTRUE,
			kTRUE,
			kTRUE,
			kFALSE,
			kFALSE,
			kTRUE,
			kTRUE,
			kTRUE,
			kFALSE};
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[9]={kFALSE,
		    kTRUE,
		    kTRUE,
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kTRUE,
		    kTRUE};
  SetVarsForOpt(4,forOpt);
  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);
}
//--------------------------------------------------------------------------
AliRDHFCutsJpsitoee::AliRDHFCutsJpsitoee(const AliRDHFCutsJpsitoee &source) :
  AliRDHFCuts(source)
{
  //
  // Copy constructor
  //

}
//--------------------------------------------------------------------------
AliRDHFCutsJpsitoee &AliRDHFCutsJpsitoee::operator=(const AliRDHFCutsJpsitoee &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliRDHFCuts::operator=(source);

  return *this;
}


//---------------------------------------------------------------------------
void AliRDHFCutsJpsitoee::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) {
  // 
  // Fills in vars the values of the variables 
  //

  if(nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsJpsitoee::GetCutsVarsForOpt: wrong number of variables\n");
    return;
  }

  AliAODRecoDecayHF2Prong *dd = (AliAODRecoDecayHF2Prong*)d;
 
  Int_t iter=-1;
  if(fVarsForOpt[0]){
    iter++;
    vars[iter]=dd->InvMassJPSIee();
  }
  if(fVarsForOpt[1]){
    iter++;
    vars[iter]=dd->GetDCA();
  }
  if(fVarsForOpt[2]){
    iter++;
    vars[iter] = dd->CosThetaStarJPSI();
  }
  if(fVarsForOpt[3]){
    iter++;
    if(pdgdaughters[0]==11) vars[iter]=dd->PtProng(0);
  }
  if(fVarsForOpt[4]){
    iter++;
    if(pdgdaughters[1]==11) vars[iter]=dd->PtProng(1);
  }
  if(fVarsForOpt[5]){
    iter++;
    vars[iter]=dd->Getd0Prong(0);
  }
  if(fVarsForOpt[6]){
    iter++;
    vars[iter]=dd->Getd0Prong(1);
  }
  if(fVarsForOpt[7]){
    iter++;
    vars[iter]= dd->Prodd0d0();
  }
  if(fVarsForOpt[8]){
    iter++;
    vars[iter]=dd->CosPointingAngle();
  }
  
  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsJpsitoee::IsSelected(TObject* obj,Int_t selectionLevel) {
  //
  // Apply selection
  //

  if(!fCutsRD){
    cout<<"Cut matrix not inizialized. Exit..."<<endl;
    return 0;
  }
  //PrintAll();
  AliAODRecoDecayHF2Prong* d=(AliAODRecoDecayHF2Prong*)obj;

  if(!d){
    cout<<"AliAODRecoDecayHF2Prong null"<<endl;
    return 0;
  }

  if(d->HasBadDaughters()) return 0;


  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(d)) return 0;
  }


  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    
    Double_t pt=d->Pt();
    Double_t mJPsi,ctsJPsi;

    Double_t mJPSIPDG = TDatabasePDG::Instance()->GetParticle(443)->Mass();
    
    Int_t ptbin=PtBin(pt);

    if(d->PtProng(1) < fCutsRD[GetGlobalIndex(3,ptbin)] || d->PtProng(0) < fCutsRD[GetGlobalIndex(4,ptbin)]) return 0;

    if(TMath::Abs(d->Getd0Prong(1)) > fCutsRD[GetGlobalIndex(5,ptbin)] ||
       TMath::Abs(d->Getd0Prong(0)) > fCutsRD[GetGlobalIndex(6,ptbin)]) return 0;

  
    if(d->GetDCA() > fCutsRD[GetGlobalIndex(1,ptbin)]) return 0;

    mJPsi=d->InvMassJPSIee();
    if(TMath::Abs(mJPsi-mJPSIPDG)    > fCutsRD[GetGlobalIndex(0,ptbin)]) return 0;

    ctsJPsi=d->CosThetaStarJPSI();
    if(TMath::Abs(ctsJPsi)    > fCutsRD[GetGlobalIndex(2,ptbin)]) return 0;

    if(d->Prodd0d0() > fCutsRD[GetGlobalIndex(7,ptbin)]) return 0;
    
    if(d->CosPointingAngle()   < fCutsRD[GetGlobalIndex(8,ptbin)]) return 0;
  }

  return 1;
}
//---------------------------------------------------------------------------
