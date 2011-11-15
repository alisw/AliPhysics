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
  fUseStrongPid(0),
  fMaxPtStrongPid(0.),
  fUseImpParProdCorrCut(kFALSE)
{
  //
  // Default Constructor
  //
  Int_t nvars=14;
  SetNVars(nvars);
  TString varNames[14]={"inv. mass [GeV]",
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
			"dca cut (cm)",
			"dec len XY (cm)",
			"cosThetaPointXY"};
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
			 kFALSE,
			 kFALSE};
  SetVarNames(nvars,varNames,isUpperCut);
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
AliRDHFCutsDplustoKpipi::AliRDHFCutsDplustoKpipi(const AliRDHFCutsDplustoKpipi &source) :
  AliRDHFCuts(source),
  fUseStrongPid(source.fUseStrongPid),
  fMaxPtStrongPid(source.fMaxPtStrongPid),
  fUseImpParProdCorrCut(source.fUseImpParProdCorrCut)
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
void AliRDHFCutsDplustoKpipi::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters,AliAODEvent *aod) {
  // 
  // Fills in vars the values of the variables 
  //


  if(nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsDplustoKpipi::GetCutsVarsForOpt: wrong number of variables\n");
    return;
  }

  AliAODRecoDecayHF3Prong *dd = (AliAODRecoDecayHF3Prong*)d;
 
  //recalculate vertex w/o daughters
  Bool_t cleanvtx=kFALSE;
  AliAODVertex *origownvtx=0x0;
  if(fRemoveDaughtersFromPrimary) {
    if(dd->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*dd->GetOwnPrimaryVtx());
    cleanvtx=kTRUE;
    if(!RecalcOwnPrimaryVtx(dd,aod)) {
      CleanOwnPrimaryVtx(dd,aod,origownvtx);
      cleanvtx=kFALSE;
    }
  }

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
    vars[iter]=dd->GetSigmaVert(aod);
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
  if(fVarsForOpt[12]){
    iter++;
    vars[iter]=dd->NormalizedDecayLengthXY()*dd->P()/dd->Pt();
  }
  if(fVarsForOpt[13]){
    iter++;
    vars[iter]=dd->CosPointingAngleXY();
  }

  if(cleanvtx)CleanOwnPrimaryVtx(dd,aod,origownvtx);

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
      if(rd->Pt()<fMaxPtStrongPid && isPion<=0 && fUseStrongPid>1)return 0;
    }
    else{//kaons
      if(isKaon<0)return 0;
	if(rd->Pt()<fMaxPtStrongPid && isKaon<=0 && fUseStrongPid>0)return 0;
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


  fIsSelectedCuts=0;
  fIsSelectedPID=0;

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

  if(fKeepSignalMC) if(IsSignalMC(d,aod,411)) return 3;

  // PID selection
  Int_t returnvaluePID=3;
  Int_t returnvalueCuts=3;

  Double_t pt=d->Pt();
  if(pt<fMinPtCand) return 0;
  if(pt>fMaxPtCand) return 0;

  if(d->HasBadDaughters()) return 0;
  
  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    
    //recalculate vertex w/o daughters
    AliAODVertex *origownvtx=0x0;
    if(fRemoveDaughtersFromPrimary) {
      if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
      if(!RecalcOwnPrimaryVtx(d,aod)) {
	CleanOwnPrimaryVtx(d,aod,origownvtx);
	return 0;
      }
    }

    Int_t ptbin=PtBin(pt);
    if (ptbin==-1) {
      CleanOwnPrimaryVtx(d,aod,origownvtx);
      return 0;
    }
    
    Double_t mDplusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
    Double_t mDplus=d->InvMassDplus();
    if(TMath::Abs(mDplus-mDplusPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}

    //2track cuts
    if(d->GetDist12toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]|| d->GetDist23toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}

    Double_t sum2=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
    if(sum2<fCutsRD[GetGlobalIndex(10,ptbin)]) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}

    if(fUseImpParProdCorrCut){
      if(d->Getd0Prong(0)*d->Getd0Prong(1)<0. && d->Getd0Prong(2)*d->Getd0Prong(1)<0.) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}
    }


    //DCA
    for(Int_t i=0;i<3;i++) if(d->GetDCA(i)>fCutsRD[GetGlobalIndex(11,ptbin)]) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}

    if(d->Pt2Prong(1) < fCutsRD[GetGlobalIndex(1,ptbin)]*fCutsRD[GetGlobalIndex(1,ptbin)] || TMath::Abs(d->Getd0Prong(1))<fCutsRD[GetGlobalIndex(3,ptbin)]) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}//Kaon

    if(d->Pt2Prong(0) < fCutsRD[GetGlobalIndex(2,ptbin)]*fCutsRD[GetGlobalIndex(2,ptbin)] || TMath::Abs(d->Getd0Prong(0))<fCutsRD[GetGlobalIndex(4,ptbin)]) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}//Pion1

    if(d->Pt2Prong(2) < fCutsRD[GetGlobalIndex(2,ptbin)]*fCutsRD[GetGlobalIndex(2,ptbin)] || TMath::Abs(d->Getd0Prong(2))<fCutsRD[GetGlobalIndex(4,ptbin)]) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}//Pion2
    
    if(d->Pt2Prong(0)<fCutsRD[GetGlobalIndex(8,ptbin)]*fCutsRD[GetGlobalIndex(8,ptbin)] && d->Pt2Prong(1)<fCutsRD[GetGlobalIndex(8,ptbin)]*fCutsRD[GetGlobalIndex(8,ptbin)] && d->Pt2Prong(2)<fCutsRD[GetGlobalIndex(8,ptbin)]*fCutsRD[GetGlobalIndex(8,ptbin)]) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}

    if(d->DecayLength2()<fCutsRD[GetGlobalIndex(7,ptbin)]*fCutsRD[GetGlobalIndex(7,ptbin)]) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}

    if(d->CosPointingAngle()< fCutsRD[GetGlobalIndex(9,ptbin)]) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}

    if(d->NormalizedDecayLengthXY()*d->P()/pt<fCutsRD[GetGlobalIndex(12,ptbin)]){CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}

    if(d->CosPointingAngleXY()<fCutsRD[GetGlobalIndex(13,ptbin)]){CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}

    //sec vert
    Double_t sigmavert=d->GetSigmaVert(aod);
    if(sigmavert>fCutsRD[GetGlobalIndex(6,ptbin)]) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}

    // unset recalculated primary vertex when not needed any more
    CleanOwnPrimaryVtx(d,aod,origownvtx);
    
    fIsSelectedCuts=returnvalueCuts;

    //if(!returnvalueCuts) return 0; // returnvalueCuts cannot be 0 here
  }
  
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate ||     
     selectionLevel==AliRDHFCuts::kPID) {
    returnvaluePID = IsSelectedPID(d);
    fIsSelectedPID=returnvaluePID;
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
  const Int_t nvars=14;
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
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[0][ipt]=0.2;
    anacutsval[3][ipt]=0.;
    anacutsval[4][ipt]=0.;
    anacutsval[5][ipt]=0.01;
    anacutsval[11][ipt]=10000000000.;
    anacutsval[12][ipt]=0.;
    anacutsval[13][ipt]=0.;
  }

  anacutsval[1][0]=0.3;
  anacutsval[1][1]=0.3;
  anacutsval[1][2]=0.4; 
  anacutsval[2][0]=0.3;
  anacutsval[2][1]=0.3;
  anacutsval[2][2]=0.4;  
  for(Int_t ipt=3;ipt<nptbins;ipt++){
    anacutsval[1][ipt]=0.4;
    anacutsval[2][ipt]=0.4;
  }
  
  anacutsval[6][0]=0.022100;
  anacutsval[6][1]=0.022100;
  anacutsval[6][2]=0.034;
  anacutsval[6][3]=0.020667;
  anacutsval[6][4]=0.020667;
  anacutsval[6][5]=0.023333;
    
  
  anacutsval[7][0]=0.08;
  anacutsval[7][1]=0.08;
  anacutsval[7][2]=0.09;  
  anacutsval[7][3]=0.095;
  anacutsval[7][4]=0.095;
   
  anacutsval[8][0]=0.5;
  anacutsval[8][1]=0.5;
  anacutsval[8][2]=1.0;
  anacutsval[8][3]=0.5;
  anacutsval[8][4]=0.5;
     
    
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

  
  for(Int_t ipt=5;ipt<nptbins;ipt++){
    anacutsval[6][ipt]=0.02333;
    anacutsval[7][ipt]=0.115;
    anacutsval[8][ipt]=0.5;
    anacutsval[10][ipt]=0.000883;
    }   


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
  const Int_t nvars=14;
    
  Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
  //Double_t cutsDplus[12]={0.2,0.4,0.4,0.,0.,0.01,0.06,0.02,0.,0.85,0.,10000000000.};

  for(Int_t ipt=0;ipt<nptbins;ipt++){
    anacutsval[0][ipt]=0.2;
    anacutsval[1][ipt]=0.8;
    anacutsval[2][ipt]=0.8;
    anacutsval[3][ipt]=0.;
    anacutsval[4][ipt]=0.;
    anacutsval[5][ipt]=0.01;
    anacutsval[11][ipt]=10000000000.;
    anacutsval[12][ipt]=0.;
    anacutsval[13][ipt]=0.;
  }
  anacutsval[1][5]=0.9;

  anacutsval[6][0]=0.022100;
  anacutsval[6][1]=0.022100;
  anacutsval[6][2]=0.034;
  anacutsval[6][3]=0.020667;
  anacutsval[6][4]=0.020667;
  anacutsval[6][5]=0.023333;
  
  anacutsval[7][0]=0.08;
  anacutsval[7][1]=0.08;
  anacutsval[7][2]=0.17;  
  anacutsval[7][3]=0.14;
  anacutsval[7][4]=0.14;
  anacutsval[7][5]=0.19;
  
  anacutsval[8][0]=0.8;
  anacutsval[8][1]=0.8;
  anacutsval[8][2]=1.1;
  anacutsval[8][3]=0.5;
  anacutsval[8][4]=0.5;
  anacutsval[8][5]=0.5;
  
  anacutsval[9][0]=0.995;
  anacutsval[9][1]=0.995;
  anacutsval[9][2]=0.997; 
  anacutsval[9][3]=0.998; 
  anacutsval[9][4]=0.998;
  anacutsval[9][5]=0.995;
  
  anacutsval[10][0]=0.0055;
  anacutsval[10][1]=0.0055;
  anacutsval[10][2]= 0.0028;
  anacutsval[10][3]=0.000883;
  anacutsval[10][4]=0.000883;
  anacutsval[10][5]=0.000883;

  anacutsval[12][5]=12.;
  anacutsval[13][5]=0.998571;
  anacutsval[12][6]=10.;
  anacutsval[13][6]=0.997143;

  for(Int_t ipt=6;ipt<nptbins;ipt++){
    anacutsval[6][ipt]=0.02333;
    anacutsval[7][ipt]=0.19;
    anacutsval[8][ipt]=2.0;
    anacutsval[9][ipt]=0.997;
    anacutsval[10][ipt]=0.000883;
  }   
  anacutsval[7][6]=0.14;
  anacutsval[9][6]=0.995;

  SetPtBins(nptbins+1,ptbins);
  SetCuts(nvars,nptbins,anacutsval);
  SetUsePID(kTRUE);
  SetMinCentrality(1E-10);
  SetMaxCentrality(20.);
  SetUseCentrality(AliRDHFCuts::kCentV0M);
  SetRemoveDaughtersFromPrim(kFALSE);
    
  PrintAll();

  for(Int_t iic=0;iic<nvars;iic++){delete [] anacutsval[iic];}
  delete [] anacutsval;
  anacutsval=NULL;

  delete [] ptbins;
  ptbins=NULL;

  return;
}
