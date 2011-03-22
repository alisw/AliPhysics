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
  if(fPidHF)delete fPidHF;
  fPidHF=new AliAODPidHF();
  Double_t plim[2]={0.6,0.8};
  Double_t nsigma[5]={2.,1.,2.,3.,0.};
  
  fPidHF->SetPLimit(plim);
  fPidHF->SetAsym(kTRUE);
  fPidHF->SetSigma(nsigma);
  fPidHF->SetMatch(1);
  fPidHF->SetTPC(1);
  fPidHF->SetTOF(1);
  fPidHF->SetITS(0);
  fPidHF->SetTRD(0);
  fPidHF->SetCompat(kTRUE);

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
    Double_t mPDGPhi = TDatabasePDG::Instance()->GetParticle(333)->Mass();
    if(TMath::Abs(pdgdaughters[0]==321)){
      
      Double_t phimass01=d->InvMass2Prongs(0,1,321,321);
       vars[iter]=TMath::Abs(phimass01-mPDGPhi);
       // vars[iter]=dd->InvMass2Prongs(0,1,321,321);
    }else{
      Double_t phimass12=d->InvMass2Prongs(1,2,321,321);
       vars[iter]=TMath::Abs(phimass12-mPDGPhi);
       // vars[iter]=dd->InvMass2Prongs(1,2,321,321);      
    }
  }
  if(fVarsForOpt[13]){
    iter++;
    Double_t mPDGK0star = TDatabasePDG::Instance()->GetParticle(313)->Mass();
    if(TMath::Abs(pdgdaughters[0]==321)){
      
      Double_t mass12kpi=d->InvMass2Prongs(1,2,321,211);
      vars[iter]=TMath::Abs(mass12kpi-mPDGK0star);
      //	      vars[iter]=dd->InvMass2Prongs(1,2,321,211);
    }else{
      Double_t mass01pik=d->InvMass2Prongs(0,1,211,321);
      vars[iter]=TMath::Abs(mass01pik-mPDGK0star);
      //	vars[iter]=dd->InvMass2Prongs(0,1,211,321);      
    }
  }

  
  return;
}
//---------------------------------------------------------------------------
Bool_t AliRDHFCutsDstoKKpi::IsInFiducialAcceptance(Double_t pt, Double_t y) const
{
  //
  // Checking if Ds is in fiducial acceptance region 
  //

  if(pt > 5.) {
    // applying cut for pt > 5 GeV
    AliDebug(2,Form("pt of Ds = %f (> 5), cutting at |y| < 0.8",pt)); 
    if (TMath::Abs(y) > 0.8) return kFALSE;
    
  } else {
    // appliying smooth cut for pt < 5 GeV
    Double_t maxFiducialY = -0.2/15*pt*pt+1.9/15*pt+0.5; 
    Double_t minFiducialY = 0.2/15*pt*pt-1.9/15*pt-0.5;		
    AliDebug(2,Form("pt of Ds = %f (< 5), cutting  according to the fiducial zone [%f, %f]\n",pt,minFiducialY,maxFiducialY)); 
    if (y < minFiducialY || y > maxFiducialY) return kFALSE;    
  }

  return kTRUE;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsDstoKKpi::IsSelectedPID(AliAODRecoDecayHF *rd) {
  // PID selection
  // return values: 0->NOT OK, 1->OK as KKpi, 2->OK as piKK, 3->OK as both 
  Int_t retCode=3;
  Bool_t okKKpi=kTRUE;
  Bool_t okpiKK=kTRUE;
  if(!fUsePID || !rd) return retCode;
  if(!fPidHF){
    AliWarning("AliAODPidHF not created!");
    return retCode;
  }
  Int_t nKaons=0;
  Int_t nNotKaons=0;
  Int_t sign= rd->GetCharge(); 
  for(Int_t iDaught=0; iDaught<3; iDaught++){
    AliAODTrack *track=(AliAODTrack*)rd->GetDaughter(iDaught);
    Int_t isPion=fPidHF->MakeRawPid(track,AliPID::kPion);
    Int_t isKaon=fPidHF->MakeRawPid(track,AliPID::kKaon);
    Int_t isProton=fPidHF->MakeRawPid(track,AliPID::kProton);
    
    if(isProton>0 &&  isKaon<0  && isPion<0) return 0;
    if(sign!=track->Charge()){// must be kaon
      if(isKaon<0) return 0;
    }
    if(isKaon>0 && isPion<0) nKaons++;
    if(isKaon<0) nNotKaons++;
    if(iDaught==0){
      if(isKaon<0) okKKpi=kFALSE;
      if(isPion<0) okpiKK=kFALSE;
    }
    else if(iDaught==2){
      if(isKaon<0) okpiKK=kFALSE;
      if(isPion<0) okKKpi=kFALSE;
    }
  }
  
  if(nKaons>2)return 0;
  if(nNotKaons>1) return 0;
  
  if(!okKKpi) retCode-=1;
  if(!okpiKK) retCode-=2;

  return retCode;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsDstoKKpi::IsSelected(TObject* obj,Int_t selectionLevel, AliAODEvent* aod) {
  //
  // Apply selection
  //

  if(!fCutsRD){
    cout<<"Cut matrix not inizialized. Exit..."<<endl;
    return 0;
  }
  //PrintAll();
  AliAODRecoDecayHF3Prong* d=(AliAODRecoDecayHF3Prong*)obj;

  if(!d){
    cout<<"AliAODRecoDecayHF3Prong null"<<endl;
    return 0;
  }
 
  Double_t ptD=d->Pt();
  if(ptD<fMinPtCand) return 0;
  if(ptD>fMaxPtCand) return 0;

  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(d)) return 0;
  }



 
  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    //recalculate vertex w/o daughters
    AliAODVertex *origownvtx=0x0;
    AliAODVertex *recvtx=0x0;
    if(fRemoveDaughtersFromPrimary) {
      if(!RecalcOwnPrimaryVtx(d,aod,origownvtx,recvtx)) return 0;
    }

    Int_t okDsKKpi=1;
    Int_t okDspiKK=1;
    Int_t okMassPhi=0;
    Int_t okMassK0star=0;

    Double_t pt=d->Pt();
    Int_t ptbin=PtBin(pt);
    if (ptbin==-1) {
      CleanOwnPrimaryVtx(d,origownvtx);
      return 0;
    }
 
    Double_t mDsPDG = TDatabasePDG::Instance()->GetParticle(431)->Mass();
    Double_t mDsKKpi=d->InvMassDsKKpi();
    Double_t mDspiKK=d->InvMassDspiKK();
    if(TMath::Abs(mDsKKpi-mDsPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) okDsKKpi = 0;
    if(TMath::Abs(mDspiKK-mDsPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) okDspiKK = 0;
    if(!okDsKKpi && !okDspiKK){
      CleanOwnPrimaryVtx(d,origownvtx);
      return 0;
    }



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
    if(!okDsKKpi && !okDspiKK){
      CleanOwnPrimaryVtx(d,origownvtx);
      return 0;
    }

    // Cuts on track pairs
    for(Int_t i=0;i<3;i++){
      if(d->GetDCA(i)>fCutsRD[GetGlobalIndex(11,ptbin)]){
	CleanOwnPrimaryVtx(d,origownvtx);
	return 0;
      }
    }
    if(d->GetDist12toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)] || 
       d->GetDist23toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]){
      CleanOwnPrimaryVtx(d,origownvtx);
      return 0;
    }



    //single track
    if(TMath::Abs(d->Pt2Prong(1)) < fCutsRD[GetGlobalIndex(1,ptbin)]*fCutsRD[GetGlobalIndex(1,ptbin)] || 
       TMath::Abs(d->Getd0Prong(1))<fCutsRD[GetGlobalIndex(3,ptbin)]){
      CleanOwnPrimaryVtx(d,origownvtx);
      return 0;
    }

    if(okDsKKpi){
      if(TMath::Abs(d->Pt2Prong(0)) < fCutsRD[GetGlobalIndex(1,ptbin)]*fCutsRD[GetGlobalIndex(1,ptbin)] || 
	 TMath::Abs(d->Getd0Prong(0))<fCutsRD[GetGlobalIndex(3,ptbin)]) okDsKKpi=0;
      if(TMath::Abs(d->Pt2Prong(2)) < fCutsRD[GetGlobalIndex(2,ptbin)]*fCutsRD[GetGlobalIndex(2,ptbin)] || 
	 TMath::Abs(d->Getd0Prong(2))<fCutsRD[GetGlobalIndex(4,ptbin)]) okDsKKpi=0;
    }
    if(okDspiKK){
      if(TMath::Abs(d->Pt2Prong(0)) < fCutsRD[GetGlobalIndex(2,ptbin)]*fCutsRD[GetGlobalIndex(2,ptbin)] || 
	 TMath::Abs(d->Getd0Prong(0))<fCutsRD[GetGlobalIndex(4,ptbin)]) okDspiKK=0;
      if(TMath::Abs(d->Pt2Prong(2)) < fCutsRD[GetGlobalIndex(1,ptbin)]*fCutsRD[GetGlobalIndex(1,ptbin)] || 
	 TMath::Abs(d->Getd0Prong(2))<fCutsRD[GetGlobalIndex(3,ptbin)]) okDspiKK=0;
    }
    if(!okDsKKpi && !okDspiKK){
      CleanOwnPrimaryVtx(d,origownvtx);
      return 0;
    }

       
    // Cuts on candidate triplet

    if(d->GetSigmaVert()>fCutsRD[GetGlobalIndex(6,ptbin)]){
      CleanOwnPrimaryVtx(d,origownvtx);
      return 0;
    }

    if(d->CosPointingAngle()< fCutsRD[GetGlobalIndex(9,ptbin)]){
      CleanOwnPrimaryVtx(d,origownvtx); 
      return 0;
    }
     
    if(d->Pt2Prong(0)<fCutsRD[GetGlobalIndex(8,ptbin)]*fCutsRD[GetGlobalIndex(8,ptbin)] && 
       d->Pt2Prong(1)<fCutsRD[GetGlobalIndex(8,ptbin)]*fCutsRD[GetGlobalIndex(8,ptbin)] && 
       d->Pt2Prong(2)<fCutsRD[GetGlobalIndex(8,ptbin)]*fCutsRD[GetGlobalIndex(8,ptbin)]) {
      CleanOwnPrimaryVtx(d,origownvtx); 
      return 0;
    }

    if(d->DecayLength2()<fCutsRD[GetGlobalIndex(7,ptbin)]*fCutsRD[GetGlobalIndex(7,ptbin)]){
      CleanOwnPrimaryVtx(d,origownvtx);
      return 0;
    }


    Double_t sum2=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
    if(sum2<fCutsRD[GetGlobalIndex(10,ptbin)]){
      CleanOwnPrimaryVtx(d,origownvtx);
      return 0;
    }

     // unset recalculated primary vertex when not needed any more
    CleanOwnPrimaryVtx(d,origownvtx);
   

    // PID selection
    Int_t returnvaluePID=3;  
    if(selectionLevel==AliRDHFCuts::kAll || 
       selectionLevel==AliRDHFCuts::kCandidate ||     
       selectionLevel==AliRDHFCuts::kPID) {
      returnvaluePID = IsSelectedPID(d);
      fIsSelectedPID=returnvaluePID;
    }
    if(returnvaluePID==0)return 0;

    Bool_t okPidDsKKpi=returnvaluePID&1;
    Bool_t okPidDspiKK=returnvaluePID&2;
    if(okPidDsKKpi && !okDsKKpi)  return 0;
    if(okPidDspiKK && !okDspiKK) return 0;

    Int_t returnvalue=0;
    if(okDsKKpi) returnvalue+=1;
    if(okDspiKK) returnvalue+=2;
    if(okMassPhi) returnvalue+=4;
    if(okMassK0star) returnvalue+=8;

    return returnvalue;
  }
  return 15;

}
//---------------------------------------------------------------------------
