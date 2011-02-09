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
AliRDHFCuts(name),
fUseStrongPid(kFALSE)
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
  AliRDHFCuts(source),
  fUseStrongPid(source.fUseStrongPid)

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
Bool_t AliRDHFCutsDplustoKpipi::IsInFiducialAcceptance(Double_t pt, Double_t y) const
{
  //
  // Checking if Dplus is in fiducial acceptance region 
  //

  if(pt > 5.) {
    // applying cut for pt > 5 GeV
    AliDebug(2,Form("pt of D+ = %f (> 5), cutting at |y| < 0.8",pt)); 
    if (TMath::Abs(y) > 0.8) return kFALSE;
    
  } else {
    // appliying smooth cut for pt < 5 GeV
    Double_t maxFiducialY = -0.2/15*pt*pt+1.9/15*pt+0.5; 
    Double_t minFiducialY = 0.2/15*pt*pt-1.9/15*pt-0.5;		
    AliDebug(2,Form("pt of D+ = %f (< 5), cutting  according to the fiducial zone [%f, %f]\n",pt,minFiducialY,maxFiducialY)); 
    if (y < minFiducialY || y > maxFiducialY) return kFALSE;    
  }

  return kTRUE;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsDplustoKpipi::IsSelectedPID(AliAODRecoDecayHF *rd)
{
  //
  // PID selection, returns 3 if accepted, 0 if not accepted
  // 
  if(!fUsePID || !rd) return 3;
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
      if(rd->Pt()<2. && isPion<=0 && fUseStrongPid)return 0;
    }
      else{//kaons
	if(isKaon<0)return 0;
	if(rd->Pt()<2. && isKaon<=0 && fUseStrongPid)return 0;
      }
    
      
  }
  
  if(nkaons>1)return 0;
  if(nNotKaons==3)return 0;
  
  return 3;   
}


//---------------------------------------------------------------------------
Int_t AliRDHFCutsDplustoKpipi::IsSelected(TObject* obj,Int_t selectionLevel, AliAODEvent* aod) {
  //
  // Apply selection, returns 3 if accepted, 0 if not accepted
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

  // PID selection
  Int_t returnvaluePID=3;
  Int_t returnvalueCuts=3;


  
  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    
    //recalculate vertex w/o daughters
    AliAODVertex *origownvtx=0x0;
    AliAODVertex *recvtx=0x0;
    if(fRemoveDaughtersFromPrimary) {
      if(!RecalcOwnPrimaryVtx(d,aod,origownvtx,recvtx)) return 0;
    }

    Double_t pt=d->Pt();
    
    Int_t ptbin=PtBin(pt);
    if (ptbin==-1) {
      CleanOwnPrimaryVtx(d,origownvtx);
      return 0;
    }
    
    Double_t mDplusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    Double_t mDplus=d->InvMassDplus();
    if(TMath::Abs(mDplus-mDplusPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) {CleanOwnPrimaryVtx(d,origownvtx); return 0;}
    //    if(d->PtProng(1) < fCutsRD[GetGlobalIndex(3,ptbin)] || d->PtProng(0) < fCutsRD[GetGlobalIndex(4,ptbin)]) okD0 = 0;
    if(TMath::Abs(d->PtProng(1)) < fCutsRD[GetGlobalIndex(1,ptbin)] || TMath::Abs(d->Getd0Prong(1))<fCutsRD[GetGlobalIndex(3,ptbin)]) {CleanOwnPrimaryVtx(d,origownvtx); return 0;}//Kaon
    if(TMath::Abs(d->PtProng(0)) < fCutsRD[GetGlobalIndex(2,ptbin)] || TMath::Abs(d->Getd0Prong(0))<fCutsRD[GetGlobalIndex(4,ptbin)]) {CleanOwnPrimaryVtx(d,origownvtx); return 0;}//Pion1
    if(TMath::Abs(d->PtProng(2)) < fCutsRD[GetGlobalIndex(2,ptbin)] || TMath::Abs(d->Getd0Prong(2))<fCutsRD[GetGlobalIndex(4,ptbin)]) {CleanOwnPrimaryVtx(d,origownvtx); return 0;}//Pion2

    

    //2track cuts
    if(d->GetDist12toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]|| d->GetDist23toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]) {CleanOwnPrimaryVtx(d,origownvtx); return 0;}
    if(d->Getd0Prong(0)*d->Getd0Prong(1)<0. && d->Getd0Prong(2)*d->Getd0Prong(1)<0.) {CleanOwnPrimaryVtx(d,origownvtx); return 0;}
    
    //sec vert
    if(d->GetSigmaVert()>fCutsRD[GetGlobalIndex(6,ptbin)]) {CleanOwnPrimaryVtx(d,origownvtx); return 0;}
    
    if(d->DecayLength()<fCutsRD[GetGlobalIndex(7,ptbin)]) {CleanOwnPrimaryVtx(d,origownvtx); return 0;}
    
    if(TMath::Abs(d->PtProng(0))<fCutsRD[GetGlobalIndex(8,ptbin)] && TMath::Abs(d->PtProng(1))<fCutsRD[GetGlobalIndex(8,ptbin)] && TMath::Abs(d->PtProng(2))<fCutsRD[GetGlobalIndex(8,ptbin)]) {CleanOwnPrimaryVtx(d,origownvtx); return 0;}
    if(d->CosPointingAngle()< fCutsRD[GetGlobalIndex(9,ptbin)]) {CleanOwnPrimaryVtx(d,origownvtx); return 0;}
    Double_t sum2=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
    if(sum2<fCutsRD[GetGlobalIndex(10,ptbin)]) {CleanOwnPrimaryVtx(d,origownvtx); return 0;}
    
    //DCA
    for(Int_t i=0;i<3;i++) if(d->GetDCA(i)>fCutsRD[GetGlobalIndex(11,ptbin)]) {CleanOwnPrimaryVtx(d,origownvtx); return 0;}

    // unset recalculated primary vertex when not needed any more
    CleanOwnPrimaryVtx(d,origownvtx);
    
    if(!returnvalueCuts) return 0;
  }
  
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate ||     
     selectionLevel==AliRDHFCuts::kPID) {
    returnvaluePID = IsSelectedPID(d);
  }
  if(returnvaluePID==0)return 0;

  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(d)) return 0;
  }
  
 


  return 3;
}




//---------------------------------------------------------------------------


void AliRDHFCutsDplustoKpipi::SetStandardCutsPP2010() {
  //
  //STANDARD CUTS USED FOR 2010 pp analysis 
  //                                            
  
  SetName("DplustoKpipiCutsStandard");
  SetTitle("Standard Cuts for D+ analysis");
  
  // PILE UP REJECTION
  SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  // EVENT CUTS
  SetMinVtxContr(1);

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4); // default is 5
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny); 
  // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3,1.e10);
  
  AddTrackCuts(esdTrackCuts);
  
 
  const Int_t nptbins =13;
  const Int_t nvars=12;
  Float_t ptbins[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=1;	
  ptbins[2]=2.;
  ptbins[3]=3.;
  ptbins[4]=4.;
  ptbins[5]=5.;
  ptbins[6]=6.;
  ptbins[7]=8.;
  ptbins[8]=10.;
  ptbins[9]=12.;
  ptbins[10]=14.;
  ptbins[11]=16.;
  ptbins[12]=24.;
  ptbins[13]=99999.;
      
    
  Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}

  //Double_t cutsDplus[12]={0.2,0.4,0.4,0.,0.,0.01,0.06,0.02,0.,0.85,0.,10000000000.};
  Int_t ic=0;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.2;
  }
  
  ic=3;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=4;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=5;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.01;
  }
  //ic=6;
  //for(Int_t ipt=0;ipt<nptbins;ipt++){
  // anacutsval[ic][ipt]=0.06;
  
  ic=11;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=10000000000.;
  }
  
  
  anacutsval[1][0]=0.3;
  anacutsval[1][1]=0.3;
  anacutsval[1][2]=0.4;
  
  anacutsval[2][0]=0.3;
  anacutsval[2][1]=0.3;
  anacutsval[2][2]=0.4;
  
  
  ic=1;
  for(Int_t ipt=3;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.4;
  }
  ic=2;
  for(Int_t ipt=3;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.4;
  }
  
  anacutsval[6][0]=0.022100;
  anacutsval[6][1]=0.022100;
  anacutsval[6][2]=0.034;
  anacutsval[6][3]=0.020667;
  anacutsval[6][4]=0.020667;
  anacutsval[6][5]=0.023333;
  anacutsval[6][6]=0.023333;
  anacutsval[6][7]=0.023333;
  anacutsval[6][8]=0.023333;
  anacutsval[6][9]=0.023333;
  anacutsval[6][10]=0.06;
  anacutsval[6][11]=0.06;
  anacutsval[6][12]=0.06;
  
  
  anacutsval[7][0]=0.08;
  anacutsval[7][1]=0.08;
  anacutsval[7][2]=0.09;  
  anacutsval[7][3]=0.095;
  anacutsval[7][4]=0.095;
  anacutsval[7][5]=0.115;
  anacutsval[7][6]=0.115;
  anacutsval[7][7]=0.115;
  anacutsval[7][8]=0.115;
  anacutsval[7][9]=0.115;
  anacutsval[7][10]=0.02;
  anacutsval[7][11]=0.02;
  anacutsval[7][12]=0.02;
  
  
  anacutsval[8][0]=0.5;
  anacutsval[8][1]=0.5;
  anacutsval[8][2]=1.0;
  anacutsval[8][3]=0.5;
  anacutsval[8][4]=0.5;
  anacutsval[8][5]=0.5;
  anacutsval[8][6]=0.5;
  anacutsval[8][7]=0.5;
  anacutsval[8][8]=0.5;
  anacutsval[8][9]=0.5;
  anacutsval[8][10]=0.2;
  anacutsval[8][11]=0.2;
  anacutsval[8][12]=0.2;
    
    
  anacutsval[9][0]=0.95;
  anacutsval[9][1]=0.95;
  anacutsval[9][2]=0.95; 
  anacutsval[9][3]=0.95; 
  anacutsval[9][4]= 0.95;
  anacutsval[9][5]=0.92;
  anacutsval[9][6]=0.92;
  anacutsval[9][7]=0.92;
  anacutsval[9][8]=0.92;
  anacutsval[9][9]=0.90;
  anacutsval[9][10]=0.90;
  anacutsval[9][11]=0.90;
  anacutsval[9][12]=0.90;
  
  anacutsval[10][0]=0.0055;
  anacutsval[10][1]=0.0055;
  anacutsval[10][2]= 0.0028;
  anacutsval[10][3]=0.000883;
  anacutsval[10][4]=0.000883;
  anacutsval[10][5]=0.000883;
  anacutsval[10][6]=0.000883;
  anacutsval[10][7]=0.000883;
  anacutsval[10][8]=0.000883;
  anacutsval[10][9]=0.000883;
  anacutsval[10][10]=0.0;
  anacutsval[10][11]=0.0;
  anacutsval[10][12]=0.0;
  
  SetGlobalIndex(nvars,nptbins);
  SetPtBins(nptbins+1,ptbins);
  SetCuts(nvars,nptbins,anacutsval);
  SetUsePID(kTRUE);
  SetRemoveDaughtersFromPrim(kTRUE);
  
  PrintAll();

  for(Int_t iic=0;iic<nvars;iic++){delete [] anacutsval[iic];}
  delete [] anacutsval;
  anacutsval=NULL;

  return;
}


void AliRDHFCutsDplustoKpipi::SetStandardCutsPbPb2010() {
  //
  //STANDARD CUTS USED FOR 2010 Pb Pb analysis.... not optimized yet 
  //                                              
  
  SetName("DplustoKpipiCutsStandard");
  SetTitle("Standard Cuts for D+ analysis in PbPb2010 run");
  
  // PILE UP REJECTION
  //SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  // EVENT CUTS
  SetMinVtxContr(1);

  
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4); // default is 5
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny); 
  // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.8,1.e10);
    
  AddTrackCuts(esdTrackCuts);
    
  const Int_t nptbins=10;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];

  ptbins[0]=0.;
  ptbins[1]=1.;
  ptbins[2]=2.;
  ptbins[3]=3.;
  ptbins[4]=4.;
  ptbins[5]=5.;
  ptbins[6]=6.;
  ptbins[7]=8.;
  ptbins[8]=12.;
  ptbins[9]=16.;
  ptbins[10]=24.;
  const Int_t nvars=12;
    
    
    
  Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
  //Double_t cutsDplus[12]={0.2,0.4,0.4,0.,0.,0.01,0.06,0.02,0.,0.85,0.,10000000000.};
  Int_t ic=0;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.2;
  }
  
  ic=3;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=4;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.;
  }
  ic=5;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.01;
  }
  //ic=6;
  //for(Int_t ipt=0;ipt<nptbins;ipt++){
  // anacutsval[ic][ipt]=0.06;
  
  ic=11;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=10000000000.;
  }
  
  
  
  
  ic=1;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.8;
  }
  ic=2;
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[ic][ipt]=0.8;
  }
  
    
  anacutsval[6][0]=0.022100;
  anacutsval[6][1]=0.022100;
  anacutsval[6][2]=0.034;
  anacutsval[6][3]=0.020667;
  anacutsval[6][4]=0.020667;
  anacutsval[6][5]=0.023333;
  anacutsval[6][6]=0.023333;
  anacutsval[6][7]=0.023333;
  anacutsval[6][8]=0.023333;
  anacutsval[6][9]=0.023333;
  
  anacutsval[7][0]=0.08;
  anacutsval[7][1]=0.08;
  anacutsval[7][2]=0.09;  
  anacutsval[7][3]=0.095;
  anacutsval[7][4]=0.095;
  anacutsval[7][5]=0.115;
  anacutsval[7][6]=0.115;
  anacutsval[7][7]=0.115;
  anacutsval[7][8]=0.115;
  anacutsval[7][9]=0.115;
  
  anacutsval[8][0]=0.5;
  anacutsval[8][1]=0.5;
  anacutsval[8][2]=1.0;
  anacutsval[8][3]=0.5;
  anacutsval[8][4]=0.5;
  anacutsval[8][5]=0.5;
  anacutsval[8][6]=0.5;
  anacutsval[8][7]=0.5;
  anacutsval[8][8]=0.5;
  anacutsval[8][9]=0.5;
  
  anacutsval[9][0]=0.995;
  anacutsval[9][1]=0.995;
  anacutsval[9][2]=0.997; 
  anacutsval[9][3]=0.989; 
  anacutsval[9][4]= 0.989;
  anacutsval[9][5]=0.997;
  anacutsval[9][6]=0.997;
  anacutsval[9][7]=0.997;
  anacutsval[9][8]=0.997;
  anacutsval[9][9]=0.997;
 
  anacutsval[10][0]=0.0055;
  anacutsval[10][1]=0.0055;
  anacutsval[10][2]= 0.0028;
  anacutsval[10][3]=0.000883;
  anacutsval[10][4]=0.000883;
  anacutsval[10][5]=0.000883;
  anacutsval[10][6]=0.000883;
  anacutsval[10][7]=0.000883;
  anacutsval[10][8]=0.000883;
  anacutsval[10][9]=0.000883;
    
    
  SetPtBins(nptbins+1,ptbins);
  SetCuts(nvars,nptbins,anacutsval);
  SetUsePID(kTRUE);
  SetRemoveDaughtersFromPrim(kFALSE);
    
  PrintAll();

  for(Int_t iic=0;iic<nvars;iic++){delete [] anacutsval[iic];}
  delete [] anacutsval;
  anacutsval=NULL;

  delete [] ptbins;
  ptbins=NULL;

  return;
}
