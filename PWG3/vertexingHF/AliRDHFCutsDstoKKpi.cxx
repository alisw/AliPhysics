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

/////////////////////////////////////////////////////////////
//
// Class for cuts on AOD reconstructed Ds->KKpi
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>

#include "AliRDHFCutsDstoKKpi.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"

ClassImp(AliRDHFCutsDstoKKpi)

//--------------------------------------------------------------------------
AliRDHFCutsDstoKKpi::AliRDHFCutsDstoKKpi(const char* name) : 
AliRDHFCuts(name)
{
  //
  // Default Constructor
  //
  Int_t nvars=14;
  SetNVars(nvars);
  TString varNames[14]={"inv. mass [GeV]",   
			"pTK [GeV/c]",
			"pTPi [GeV/c]",
			"d0K [cm]",
			"d0Pi [cm]",
			"dist12 [cm]",
			"sigmavert [cm]",
			"decLen [cm]",
			"ptMax [GeV/c]",
			"cosThetaPoint",
			"Sum d0^2 (cm^2)",
			"dca [cm]",
			"inv. mass (Mphi-MKK) [GeV]",
			"inv. mass (MKo*-MKpi) [GeV]"};
  Bool_t isUpperCut[14]={kTRUE,
			 kFALSE,
			 kFALSE,
			 kFALSE,
			 kFALSE,
			 kFALSE,
			 kTRUE,
			 kFALSE,
			 kFALSE,
			 kFALSE,
			 kFALSE,
			 kTRUE,
			 kTRUE,
			 kTRUE};
  SetVarNames(14,varNames,isUpperCut);
  Bool_t forOpt[14]={kFALSE,
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
		    kFALSE,
		    kTRUE,
		    kTRUE};
  SetVarsForOpt(7,forOpt);
  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);
}
//--------------------------------------------------------------------------
AliRDHFCutsDstoKKpi::AliRDHFCutsDstoKKpi(const AliRDHFCutsDstoKKpi &source) :
  AliRDHFCuts(source)
{
  //
  // Copy constructor
  //

}
//--------------------------------------------------------------------------
AliRDHFCutsDstoKKpi &AliRDHFCutsDstoKKpi::operator=(const AliRDHFCutsDstoKKpi &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliRDHFCuts::operator=(source);

  return *this;
}


//---------------------------------------------------------------------------
void AliRDHFCutsDstoKKpi::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) {
  // 
  // Fills in vars the values of the variables 
  //

  if(nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsDstoKKpi::GetCutsVarsForOpt: wrong number of variables\n");
    return;
  }

  AliAODRecoDecayHF3Prong *dd = (AliAODRecoDecayHF3Prong*)d;
 
  Int_t iter=-1;
  if(fVarsForOpt[0]){
    iter++;
    if(TMath::Abs(pdgdaughters[0]==321)){
      vars[iter]=dd->InvMassDsKKpi();
    }else{
      vars[iter]=dd->InvMassDspiKK();
    }
  }
  if(fVarsForOpt[1]){
    iter++;
    Float_t minPtDau=99999.;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==321 && 
	 dd->PtProng(iprong)<minPtDau) minPtDau=dd->PtProng(iprong);
    }
    vars[iter]=minPtDau;
  }
  if(fVarsForOpt[2]){
    iter++;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==211) {
	vars[iter]=dd->PtProng(iprong);
      }
    }
  }
  if(fVarsForOpt[3]){
    iter++;
    Float_t minImpParDau=99999.;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==321 &&
	 dd->Getd0Prong(iprong)<minImpParDau) minImpParDau=dd->Getd0Prong(iprong);
    }
    vars[iter]=minImpParDau;
  }
  if(fVarsForOpt[4]){
    iter++;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==211) {
	vars[iter]=dd->Getd0Prong(iprong);
      }
    }
  }
  if(fVarsForOpt[5]){
    iter++;
    Float_t minDistPair=TMath::Min(dd->GetDist12toPrim(),dd->GetDist23toPrim());
    vars[iter]=minDistPair;
  }
  if(fVarsForOpt[6]){
    iter++;
    vars[iter]=dd->GetSigmaVert();
  }
  if(fVarsForOpt[7]){
    iter++;
    vars[iter] = dd->DecayLength();
  }
  if(fVarsForOpt[8]){
    iter++;
    Float_t ptmax=0;
    for(Int_t i=0;i<3;i++){
      if(dd->PtProng(i)>ptmax)ptmax=dd->PtProng(i);
    }
    vars[iter]=ptmax;
  }
  if(fVarsForOpt[9]){
    iter++;
    vars[iter]=dd->CosPointingAngle();
  }
  if(fVarsForOpt[10]){
    iter++;
    vars[iter]=dd->Getd0Prong(0)*dd->Getd0Prong(0)+dd->Getd0Prong(1)*dd->Getd0Prong(1)+dd->Getd0Prong(2)*dd->Getd0Prong(2);
  }
  if(fVarsForOpt[11]){
    iter++;
    Float_t maxDCA=0.;
    for(Int_t i=0;i<3;i++){ 
      if(d->GetDCA(i)>maxDCA) maxDCA=d->GetDCA(i);
    }
    vars[iter]=maxDCA;
  }
  if(fVarsForOpt[12]){
    iter++;
    if(TMath::Abs(pdgdaughters[0]==321)){
      vars[iter]=dd->InvMass2Prongs(0,1,321,321);
    }else{
      vars[iter]=dd->InvMass2Prongs(1,2,321,321);      
    }
  }
  if(fVarsForOpt[13]){
    iter++;
    if(TMath::Abs(pdgdaughters[0]==321)){
      vars[iter]=dd->InvMass2Prongs(1,2,321,211);
    }else{
      vars[iter]=dd->InvMass2Prongs(0,1,211,321);      
    }
  }

  
  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsDstoKKpi::IsSelected(TObject* obj,Int_t selectionLevel) {
  //
  // Apply selection
  //

  if(!fCutsRD){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  //PrintAll();
  AliAODRecoDecayHF3Prong* d=(AliAODRecoDecayHF3Prong*)obj;

  if(!d){
    cout<<"AliAODRecoDecayHF3Prong null"<<endl;
    return 0;
  }


  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(d)) return 0;
  }



  Int_t returnvalue=1;
  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    

    Int_t okDsKKpi=1;
    Int_t okDspiKK=1;
    Int_t okMassPhi=0;
    Int_t okMassK0star=0;

    Double_t pt=d->Pt();
    Int_t ptbin=PtBin(pt);

    Double_t mDsPDG = TDatabasePDG::Instance()->GetParticle(431)->Mass();
    Double_t mDsKKpi=d->InvMassDsKKpi();
    Double_t mDspiKK=d->InvMassDspiKK();
    if(TMath::Abs(mDsKKpi-mDsPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) okDsKKpi = 0;
    if(TMath::Abs(mDspiKK-mDsPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) okDspiKK = 0;
    if(!okDsKKpi && !okDspiKK) return 0;

    //single track
    if(TMath::Abs(d->PtProng(1)) < fCutsRD[GetGlobalIndex(1,ptbin)] || 
       TMath::Abs(d->Getd0Prong(1))<fCutsRD[GetGlobalIndex(3,ptbin)]) return 0;
    if(okDsKKpi){
      if(TMath::Abs(d->PtProng(0)) < fCutsRD[GetGlobalIndex(1,ptbin)] || 
	 TMath::Abs(d->Getd0Prong(0))<fCutsRD[GetGlobalIndex(3,ptbin)]) okDsKKpi=0;
      if(TMath::Abs(d->PtProng(2)) < fCutsRD[GetGlobalIndex(2,ptbin)] || 
	 TMath::Abs(d->Getd0Prong(2))<fCutsRD[GetGlobalIndex(4,ptbin)]) okDsKKpi=0;
    }
    if(okDspiKK){
      if(TMath::Abs(d->PtProng(0)) < fCutsRD[GetGlobalIndex(2,ptbin)] || 
	 TMath::Abs(d->Getd0Prong(0))<fCutsRD[GetGlobalIndex(4,ptbin)]) okDspiKK=0;
      if(TMath::Abs(d->PtProng(2)) < fCutsRD[GetGlobalIndex(1,ptbin)] || 
	 TMath::Abs(d->Getd0Prong(2))<fCutsRD[GetGlobalIndex(3,ptbin)]) okDspiKK=0;
    }
    if(!okDsKKpi && !okDspiKK) return 0;

    // cuts on resonant decays (via Phi or K0*)
    Double_t mPhiPDG = TDatabasePDG::Instance()->GetParticle(333)->Mass();
    Double_t mK0starPDG = TDatabasePDG::Instance()->GetParticle(313)->Mass();
    if(okDsKKpi){
      Double_t mass01phi=d->InvMass2Prongs(0,1,321,321);
      Double_t mass12K0s=d->InvMass2Prongs(1,2,321,211);
      if(TMath::Abs(mass01phi-mPhiPDG)<fCutsRD[GetGlobalIndex(12,ptbin)]) okMassPhi=1;
      if(TMath::Abs(mass12K0s-mK0starPDG)<fCutsRD[GetGlobalIndex(13,ptbin)]) okMassK0star = 1;
      if(!okMassPhi && !okMassK0star) okDsKKpi=0;
    }
    if(okDspiKK){
      Double_t mass01K0s=d->InvMass2Prongs(0,1,211,321);
      Double_t mass12phi=d->InvMass2Prongs(1,2,321,321);
      if(TMath::Abs(mass01K0s-mK0starPDG)<fCutsRD[GetGlobalIndex(13,ptbin)]) okMassK0star = 1;
      if(TMath::Abs(mass12phi-mPhiPDG)<fCutsRD[GetGlobalIndex(12,ptbin)]) okMassPhi=1;
      if(!okMassPhi && !okMassK0star) okDspiKK=0;
    }
    if(!okDsKKpi && !okDspiKK) return 0;

    // Cuts on track pairs
    for(Int_t i=0;i<3;i++) if(d->GetDCA(i)>fCutsRD[GetGlobalIndex(11,ptbin)])  return 0;
    if(d->GetDist12toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)] || 
       d->GetDist23toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]) return 0;


    // Cuts on candidate triplet
    if(d->GetSigmaVert()>fCutsRD[GetGlobalIndex(6,ptbin)]) return 0;
    if(d->DecayLength()<fCutsRD[GetGlobalIndex(7,ptbin)]) return 0;
    if(TMath::Abs(d->PtProng(0))<fCutsRD[GetGlobalIndex(8,ptbin)] && 
       TMath::Abs(d->PtProng(1))<fCutsRD[GetGlobalIndex(8,ptbin)] && 
       TMath::Abs(d->PtProng(2))<fCutsRD[GetGlobalIndex(8,ptbin)]) return 0;
    if(d->CosPointingAngle()< fCutsRD[GetGlobalIndex(9,ptbin)])return 0;
    Double_t sum2=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
    if(sum2<fCutsRD[GetGlobalIndex(10,ptbin)])return 0;

    returnvalue=0;
    if(okDsKKpi) returnvalue+=1;
    if(okDspiKK) returnvalue+=2;
    if(okMassPhi) returnvalue+=4;
    if(okMassK0star) returnvalue+=8;

  }
  return returnvalue;

}
//---------------------------------------------------------------------------
