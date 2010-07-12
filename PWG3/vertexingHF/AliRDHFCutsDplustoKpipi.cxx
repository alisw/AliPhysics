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
// Class for cuts on AOD reconstructed D+->Kpipi
//
// Author: R. Bala, bala@to.infn.it
//         G. Ortona, ortona@to.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>

#include "AliRDHFCutsDplustoKpipi.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"


ClassImp(AliRDHFCutsDplustoKpipi)

//--------------------------------------------------------------------------
AliRDHFCutsDplustoKpipi::AliRDHFCutsDplustoKpipi(const char* name) : 
  AliRDHFCuts(name)
{
  //
  // Default Constructor
  //
  Int_t nvars=12;
  SetNVars(nvars);
  TString varNames[12]={"inv. mass [GeV]",
			"pTK [GeV/c]",
			"pTPi [GeV/c]",
			"d0K [cm]   lower limit!",
			"d0Pi [cm]  lower limit!",
			"dist12 (cm)",
			"sigmavert (cm)",
			"dist prim-sec (cm)",
			"pM=Max{pT1,pT2,pT3} (GeV/c)",
			"cosThetaPoint",
			"Sum d0^2 (cm^2)",
			"dca cut (cm)"};
  Bool_t isUpperCut[12]={kTRUE,
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
			 kTRUE};
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[12]={kFALSE,
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
		     kFALSE};
  SetVarsForOpt(5,forOpt);
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
AliRDHFCutsDplustoKpipi::AliRDHFCutsDplustoKpipi(const AliRDHFCutsDplustoKpipi &source) :
  AliRDHFCuts(source)
{
  //
  // Copy constructor
  //

}
//--------------------------------------------------------------------------
AliRDHFCutsDplustoKpipi &AliRDHFCutsDplustoKpipi::operator=(const AliRDHFCutsDplustoKpipi &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliRDHFCuts::operator=(source);

  return *this;
}
//


//---------------------------------------------------------------------------
void AliRDHFCutsDplustoKpipi::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) {
  // 
  // Fills in vars the values of the variables 
  //


  if(nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsDplustoKpipi::GetCutsVarsForOpt: wrong number of variables\n");
    return;
  }

  AliAODRecoDecayHF3Prong *dd = (AliAODRecoDecayHF3Prong*)d;

  Int_t iter=-1;
  if(fVarsForOpt[0]){
    iter++;
    vars[iter]=dd->InvMassDplus();
  }
  if(fVarsForOpt[1]){
    iter++;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==321) {
	vars[iter]=dd->PtProng(iprong);
      }
    }
  }
  if(fVarsForOpt[2]){
    iter++;
    Float_t minPtDau=1000000.0;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==211) {
	if(dd->PtProng(iprong)<minPtDau){
	  minPtDau=dd->PtProng(iprong);
	}
      }
    }
    vars[iter]=minPtDau;
  }
  if(fVarsForOpt[3]){
    iter++;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==321) {
	vars[iter]=dd->Getd0Prong(iprong);
      }
    }
  }
  if(fVarsForOpt[4]){
    iter++;
    Float_t minImpParDau=1000000.0;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==211) {
	if(dd->Getd0Prong(iprong)<minImpParDau){
	  minImpParDau=dd->Getd0Prong(iprong);
	}
      }
    }
   vars[iter]=minImpParDau;
  }
  if(fVarsForOpt[5]){
    iter++;
    Float_t dist12 = dd->GetDist12toPrim();
    Float_t dist23 = dd->GetDist23toPrim();
    if(dist12<dist23)vars[iter]=dist12;
    else vars[iter]=dist23;
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
    Float_t maxDCA=0;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(dd->GetDCA(iprong)<maxDCA){
	maxDCA=dd->GetDCA(iprong);
      }
    }
    vars[iter]=maxDCA;
  }
  return;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsDplustoKpipi::IsSelectedPID(AliAODRecoDecayHF *rd)const {
  //PID 
  if(!fUsePID || !rd) return 1;
  //if(fUsePID)printf("i am inside the pid \n");
  Int_t nkaons=0;
  Int_t nNotKaons=0;
  Int_t sign= rd->GetCharge(); 
  for(Int_t daught=0;daught<3;daught++){
    AliAODTrack *track=(AliAODTrack*)rd->GetDaughter(daught);
    Int_t isPion=fPidHF->MakeRawPid(track,AliPID::kPion);
    Int_t isKaon=fPidHF->MakeRawPid(track,AliPID::kKaon);
    Int_t isProton=fPidHF->MakeRawPid(track,AliPID::kProton);
    
    if(isProton>0 &&  isKaon<0  && isPion<0) return 0;
    if(isKaon>0 && isPion<0) nkaons++;
    if(isKaon<0) nNotKaons++;  
    if(sign==track->Charge()){//pions
      if(isPion<0)return 0;
    }
      else{//kaons
	if(isKaon<0)return 0;
      }
    
      
  }
  
  if(nkaons>1)return 0;
  if(nNotKaons==3)return 0;
  
  return 1;   
}



//---------------------------------------------------------------------------
Int_t AliRDHFCutsDplustoKpipi::IsSelected(TObject* obj,Int_t selectionLevel) {
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


 
  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(d)) return 0;
  }
  
  // PID selection
  Int_t returnvaluePID=1;  
                                          

  //if(selectionLevel==AliRDHFCuts::kAll || 
  if(selectionLevel==AliRDHFCuts::kCandidate ||     
 selectionLevel==AliRDHFCuts::kPID) {
    returnvaluePID = IsSelectedPID(d);
  }
  if(returnvaluePID==0)return 0;
  
  
  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    
    Double_t pt=d->Pt();
    
    Int_t ptbin=PtBin(pt);
    
    Double_t mDplusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    Double_t mDplus=d->InvMassDplus();
    if(TMath::Abs(mDplus-mDplusPDG)>fCutsRD[GetGlobalIndex(0,ptbin)])return 0;
    //    if(d->PtProng(1) < fCutsRD[GetGlobalIndex(3,ptbin)] || d->PtProng(0) < fCutsRD[GetGlobalIndex(4,ptbin)]) okD0 = 0;
    if(TMath::Abs(d->PtProng(1)) < fCutsRD[GetGlobalIndex(1,ptbin)] || TMath::Abs(d->Getd0Prong(1))<fCutsRD[GetGlobalIndex(3,ptbin)])return 0;//Kaon
    if(TMath::Abs(d->PtProng(0)) < fCutsRD[GetGlobalIndex(2,ptbin)] || TMath::Abs(d->Getd0Prong(0))<fCutsRD[GetGlobalIndex(4,ptbin)])return 0;//Pion1
    if(TMath::Abs(d->PtProng(2)) < fCutsRD[GetGlobalIndex(2,ptbin)] || TMath::Abs(d->Getd0Prong(2))<fCutsRD[GetGlobalIndex(4,ptbin)])return 0;//Pion2

    

  //2track cuts
  if(d->GetDist12toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]|| d->GetDist23toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)])return 0;
  if(d->Getd0Prong(0)*d->Getd0Prong(1)<0. && d->Getd0Prong(2)*d->Getd0Prong(1)<0.)return 0;

  //sec vert
  if(d->GetSigmaVert()>fCutsRD[GetGlobalIndex(6,ptbin)])return 0;

  if(d->DecayLength()<fCutsRD[GetGlobalIndex(7,ptbin)])return 0;

  if(TMath::Abs(d->PtProng(0))<fCutsRD[GetGlobalIndex(8,ptbin)] && TMath::Abs(d->PtProng(1))<fCutsRD[GetGlobalIndex(8,ptbin)] && TMath::Abs(d->PtProng(2))<fCutsRD[GetGlobalIndex(8,ptbin)])return 0;
  if(d->CosPointingAngle()< fCutsRD[GetGlobalIndex(9,ptbin)])return 0;
  Double_t sum2=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
  if(sum2<fCutsRD[GetGlobalIndex(10,ptbin)])return 0;

  //DCA
  for(Int_t i=0;i<3;i++) if(d->GetDCA(i)>fCutsRD[GetGlobalIndex(11,ptbin)]) return 0;
  
  return 1;
  }
  return 1;
}
//---------------------------------------------------------------------------
