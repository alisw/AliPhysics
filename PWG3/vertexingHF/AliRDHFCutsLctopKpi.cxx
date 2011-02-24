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
// Author: A.Dainese, andrea.dainese@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>

#include "AliRDHFCutsLctopKpi.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliRDHFCuts.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliKFParticle.h"
#include "AliESDVertex.h"

ClassImp(AliRDHFCutsLctopKpi)

//--------------------------------------------------------------------------
AliRDHFCutsLctopKpi::AliRDHFCutsLctopKpi(const char* name) : 
AliRDHFCuts(name),
fPidObjprot(0),
fPidObjpion(0),
fRecoKF(kFALSE)
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
}
//--------------------------------------------------------------------------
AliRDHFCutsLctopKpi::AliRDHFCutsLctopKpi(const AliRDHFCutsLctopKpi &source) :
  AliRDHFCuts(source),
  fPidObjprot(0),
  fPidObjpion(0),
  fRecoKF(kFALSE)
{
  //
  // Copy constructor
  //
  if(source.fPidObjprot) SetPidprot(source.fPidObjprot);
  if(source.fPidObjpion) SetPidpion(source.fPidObjpion);
}
//--------------------------------------------------------------------------
AliRDHFCutsLctopKpi &AliRDHFCutsLctopKpi::operator=(const AliRDHFCutsLctopKpi &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliRDHFCuts::operator=(source);
  SetPidprot(source.GetPidprot());
  SetPidpion(source.GetPidpion());

  return *this;
}
//---------------------------------------------------------------------------
AliRDHFCutsLctopKpi::~AliRDHFCutsLctopKpi() {
 //
 //  // Default Destructor
 //   
 if(fPidObjpion){
  delete fPidObjpion;
  fPidObjpion=0;
 }
 if(fPidObjprot){
  delete fPidObjprot;
  fPidObjprot=0;
 }

}

//---------------------------------------------------------------------------
void AliRDHFCutsLctopKpi::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) {
  // 
  // Fills in vars the values of the variables 
  //

  if(nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsLctopKpi::GetCutsVarsForOpt: wrong number of variables\n");
    return;
  }

  AliAODRecoDecayHF3Prong *dd = (AliAODRecoDecayHF3Prong*)d;

    Int_t iter=-1;
  if(fVarsForOpt[0]){
    iter++;
    vars[iter]=dd->InvMassLcpKpi();
  }
  if(fVarsForOpt[1]){
    iter++;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==2212) {
	vars[iter]=dd->PtProng(iprong);
      }
    }
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
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==2212) {
	vars[iter]=dd->Getd0Prong(iprong);
      }
    }
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
    vars[iter]=dd->GetDist12toPrim();
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
    vars[iter]=dd->GetDCA();
  }

  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctopKpi::IsSelected(TObject* obj,Int_t selectionLevel,AliAODEvent *aod) {
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



  Int_t returnvalue=3;
  Int_t returnvaluePID=3;

  if(d->Pt()<fMinPtCand) return 0;
  if(d->Pt()>fMaxPtCand) return 0;

  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {

    Double_t pt=d->Pt();
    
    Int_t ptbin=PtBin(pt);
    
    Double_t mLcpKpi,mLcpiKp;
    Int_t okLcpKpi=1,okLcpiKp=1;

    Double_t mLcPDG = TDatabasePDG::Instance()->GetParticle(4122)->Mass();

    mLcpKpi=d->InvMassLcpKpi();
    mLcpiKp=d->InvMassLcpiKp();

    if(TMath::Abs(mLcpKpi-mLcPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) okLcpKpi = 0;
    if(TMath::Abs(mLcpiKp-mLcPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) okLcpiKp = 0;
    if(!okLcpKpi && !okLcpiKp) return 0;

    if(TMath::Abs(d->PtProng(1)) < fCutsRD[GetGlobalIndex(1,ptbin)] || TMath::Abs(d->Getd0Prong(1))<fCutsRD[GetGlobalIndex(3,ptbin)]) return 0;//Kaon
    //if(TMath::Abs(d->PtProng(0)) < fCutsRD[GetGlobalIndex(2,ptbin)] || TMath::Abs(d->Getd0Prong(0))<fCutsRD[GetGlobalIndex(4,ptbin)]) return 0;//Proton
    //if(TMath::Abs(d->PtProng(2)) < fCutsRD[GetGlobalIndex(2,ptbin)] || TMath::Abs(d->Getd0Prong(2))<fCutsRD[GetGlobalIndex(4,ptbin)]) return 0;//Pion
    if((TMath::Abs(d->PtProng(0)) < fCutsRD[GetGlobalIndex(2,ptbin)]) || (TMath::Abs(d->PtProng(2)) < 0.4)) okLcpKpi=0;
    if((TMath::Abs(d->PtProng(2)) < fCutsRD[GetGlobalIndex(2,ptbin)]) || (TMath::Abs(d->PtProng(0)) < 0.4))okLcpiKp=0;
    if(!okLcpKpi && !okLcpiKp) return 0;

    
    if(fRecoKF){
     Int_t valueTmp=3;
     if(okLcpKpi) returnvalue=1; //cuts passed as Lc->pKpi
     if(okLcpiKp) returnvalue=2; //cuts passed as Lc->piKp
     if(okLcpKpi && okLcpiKp) returnvalue=3; //cuts passed as both pKpi and piKp

     Int_t valueTotTmp=CombinePIDCuts(valueTmp,returnvaluePID);
     Int_t pdgs[3]={2212,321,211};
     if(valueTotTmp>=2) {
      pdgs[0]=211;
      pdgs[2]=2212;
     }
     if(!d->GetOwnPrimaryVtx()){
      AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex(); 
      d->SetOwnPrimaryVtx(vtx1);
     }
     Double_t field=aod->GetMagneticField(); 
     ReconstructKF(d,pdgs,field);
    }
    //2track cuts
    if(d->GetDist12toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]|| d->GetDist23toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]) return 0;
    if(d->GetDist12toPrim()>1.) return 0;
    if(d->GetDist23toPrim()>1.) return 0;
    if(d->Getd0Prong(0)*d->Getd0Prong(1)<0. && d->Getd0Prong(2)*d->Getd0Prong(1)<0.) return 0;
    
    //sec vert
    if(d->GetSigmaVert()>fCutsRD[GetGlobalIndex(6,ptbin)]) return 0;

    if(d->DecayLength()<fCutsRD[GetGlobalIndex(7,ptbin)]) return 0;
    if(d->DecayLength()>0.5) return 0;
    
    if(TMath::Abs(d->PtProng(0))<fCutsRD[GetGlobalIndex(8,ptbin)] && TMath::Abs(d->PtProng(1))<fCutsRD[GetGlobalIndex(8,ptbin)] && TMath::Abs(d->PtProng(2))<fCutsRD[GetGlobalIndex(8,ptbin)]) return 0;
    if(d->CosPointingAngle()< fCutsRD[GetGlobalIndex(9,ptbin)]) return 0;
    Double_t sum2=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
    if(sum2<fCutsRD[GetGlobalIndex(10,ptbin)]) return 0;
    
    //DCA
    for(Int_t i=0;i<3;i++) if(d->GetDCA(i)>fCutsRD[GetGlobalIndex(11,ptbin)]) return 0;


    if(okLcpKpi) returnvalue=1; //cuts passed as Lc->pKpi
    if(okLcpiKp) returnvalue=2; //cuts passed as Lc->piKp
    if(okLcpKpi && okLcpiKp) returnvalue=3; //cuts passed as both pKpi and piKp
   
  }

  if(fUsePID || selectionLevel==AliRDHFCuts::kPID) returnvaluePID = IsSelectedPID(d);
  if(returnvaluePID==0) return 0;

  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(d)) return 0;
  }


  Int_t returnvalueTot=CombinePIDCuts(returnvalue,returnvaluePID);
  return returnvalueTot;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctopKpi::IsSelectedPID(AliAODRecoDecayHF* obj) {


    if(!obj) {return 3;}
    Int_t okLcpKpi=0,okLcpiKp=0;
    Int_t returnvalue=0;
    Bool_t isPeriodd=fPidHF->GetOnePad();
    Bool_t isPbPb=fPidHF->GetPbPb();
    Bool_t ispion0=kTRUE,ispion2=kTRUE;
    Bool_t isproton0=kFALSE,isproton2=kFALSE;
    Bool_t iskaon1=kFALSE;
    Double_t sigmaTOF=120.;
    if(isPeriodd) sigmaTOF=160.;
    if(isPbPb) sigmaTOF=160.;
    fPidObjprot->SetTofSigma(sigmaTOF);
    fPidHF->SetTofSigma(sigmaTOF);

    for(Int_t i=0;i<3;i++){
     AliAODTrack *track=(AliAODTrack*)obj->GetDaughter(i);
     if(!track) return 0;
     // identify kaon
     if(i==1) {
      Int_t isKaon=fPidHF->MakeRawPid(track,3); 
      if(isKaon>=1) {
       iskaon1=kTRUE;
      if(fPidHF->MakeRawPid(track,2)>=1) iskaon1=kFALSE;
      }
      if(!iskaon1) return 0;
     
     }else{
     //pion or proton
      
     Int_t isProton=fPidObjprot->MakeRawPid(track,4);
     if(isProton>=1){
      if(fPidHF->MakeRawPid(track,2)>=1) isProton=-1;
      if(fPidHF->MakeRawPid(track,3)>=1) isProton=-1;
     }

     Int_t isPion=fPidObjpion->MakeRawPid(track,2);
     if(fPidHF->MakeRawPid(track,3)>=1) isPion=-1;
     if(fPidObjprot->MakeRawPid(track,4)>=1) isPion=-1;


     if(i==0) {
      if(isPion<0) ispion0=kFALSE;
      if(isProton>=1) isproton0=kTRUE;

     }
      if(!ispion0 && !isproton0) return 0;
     if(i==2) {
      if(isPion<0) ispion2=kFALSE;
      if(isProton>=1) isproton2=kTRUE;
     }

    }
   }

    if(ispion2 && isproton0 && iskaon1) okLcpKpi=1;
    if(ispion0 && isproton2 && iskaon1) okLcpiKp=1;
    if(okLcpKpi) returnvalue=1; //cuts passed as Lc->pKpi
    if(okLcpiKp) returnvalue=2; //cuts passed as Lc->piKp
    if(okLcpKpi && okLcpiKp) returnvalue=3; //cuts passed as both pKpi and piKp

 return returnvalue;
}
//-----------------------
Int_t AliRDHFCutsLctopKpi::CombinePIDCuts(Int_t returnvalue, Int_t returnvaluePID) const {

 Int_t returnvalueTot=0;
 Int_t okLcpKpi=0,okLcpiKp=0;
 if(returnvaluePID==1){
   if(returnvalue==1 || returnvalue==3) okLcpKpi=1;
 }
 if(returnvaluePID==2){
   if(returnvalue>=2) okLcpiKp=1;
 }
 if(returnvaluePID==3 && returnvalue>0){
  if(returnvalue==1 || returnvalue==3) okLcpKpi=1;
  if(returnvalue>=2) okLcpiKp=1;
 } 

 if(okLcpKpi) returnvalueTot=1; //cuts passed as Lc->pKpi
 if(okLcpiKp) returnvalueTot=2; //cuts passed as Lc->piKp
 if(okLcpKpi && okLcpiKp) returnvalueTot=3; //cuts passed as both pKpi and piKp
 return returnvalueTot;
}
//----------------------------------
void AliRDHFCutsLctopKpi::SetStandardCutsPP2010() {

 SetName("LctopKpiProdCuts");
 SetTitle("Production cuts for Lc analysis");

 AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
 esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
 esdTrackCuts->SetRequireTPCRefit(kTRUE);
 esdTrackCuts->SetMinNClustersTPC(70);
 esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                          AliESDtrackCuts::kAny);
 esdTrackCuts->SetRequireITSRefit(kTRUE);
 esdTrackCuts->SetMinNClustersITS(4);
 esdTrackCuts->SetMinDCAToVertexXY(0.);
 esdTrackCuts->SetEtaRange(-0.8,0.8);
 esdTrackCuts->SetPtRange(0.3,1.e10);
 AddTrackCuts(esdTrackCuts);

 const Int_t nptbins=4;
 const Int_t nvars=12;
 Float_t* ptbins;
 ptbins=new Float_t[nptbins+1];
 
 ptbins[0]=0.;
 ptbins[1]=2.;
 ptbins[2]=3.;
 ptbins[3]=4.;
 ptbins[4]=9999.;

 SetGlobalIndex(nvars,nptbins);
 SetPtBins(nptbins+1,ptbins);

 Float_t** prodcutsval;
 prodcutsval=new Float_t*[nvars];
 for(Int_t iv=0;iv<nvars;iv++){
  prodcutsval[iv]=new Float_t[nptbins];
 }

 for(Int_t ipt=0;ipt<nptbins;ipt++){
  prodcutsval[0][ipt]=0.18;
  prodcutsval[1][ipt]=0.4;
  prodcutsval[2][ipt]=0.5;
  prodcutsval[3][ipt]=0.;
  prodcutsval[4][ipt]=0.;
  prodcutsval[5][ipt]=0.01;
  prodcutsval[6][ipt]=0.06;
  prodcutsval[7][ipt]=0.005;
  prodcutsval[8][ipt]=0.;
  prodcutsval[9][ipt]=0.;
  prodcutsval[10][ipt]=0.;
  prodcutsval[11][ipt]=0.05;
 }
 SetCuts(nvars,nptbins,prodcutsval);

 AliAODPidHF* pidObjK=new AliAODPidHF();
 Double_t sigmasK[5]={3.,1.,1.,3.,2.};
 pidObjK->SetSigma(sigmasK);
 pidObjK->SetAsym(kTRUE);
 pidObjK->SetMatch(1);
 pidObjK->SetTPC(kTRUE);
 pidObjK->SetTOF(kTRUE);
 pidObjK->SetITS(kTRUE);
 Double_t plimK[2]={0.5,0.8};
 pidObjK->SetPLimit(plimK,2);

 SetPidHF(pidObjK);

 AliAODPidHF* pidObjpi=new AliAODPidHF();
 pidObjpi->SetTPC(kTRUE);
 Double_t sigmaspi[5]={3.,0.,0.,0.,0.};
 pidObjpi->SetSigma(sigmaspi);
 SetPidpion(pidObjpi);

 AliAODPidHF* pidObjp=new AliAODPidHF();
 Double_t sigmasp[5]={3.,1.,1.,3.,2.};
 pidObjp->SetSigma(sigmasp);
 pidObjp->SetAsym(kTRUE);
 pidObjp->SetMatch(1);
 pidObjp->SetTPC(kTRUE);
 pidObjp->SetTOF(kTRUE);
 pidObjp->SetITS(kTRUE);
 Double_t plimp[2]={1.,2.};
 pidObjp->SetPLimit(plimp,2);

 SetPidprot(pidObjp);

 SetUsePID(kTRUE);

 PrintAll();

 for(Int_t iiv=0;iiv<nvars;iiv++){
  delete [] prodcutsval[iiv];
 }
 delete [] prodcutsval;
 prodcutsval=NULL;
 delete [] ptbins;
 ptbins=NULL;

 delete pidObjp;
 pidObjp=NULL;

 return;
}
//------------------
void AliRDHFCutsLctopKpi::SetStandardCutsPbPb2010() {

 SetName("LctopKpiProdCuts");
 SetTitle("Production cuts for Lc analysis");

 AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");

 esdTrackCuts->SetRequireTPCRefit(kTRUE);
 esdTrackCuts->SetMinNClustersTPC(70);
 esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                          AliESDtrackCuts::kAny);
 esdTrackCuts->SetRequireITSRefit(kTRUE);
 esdTrackCuts->SetMinNClustersITS(4);
 esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0100*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
 esdTrackCuts->SetEtaRange(-0.8,0.8);
 esdTrackCuts->SetMaxDCAToVertexXY(1.);
 esdTrackCuts->SetMaxDCAToVertexZ(1.);
 esdTrackCuts->SetPtRange(0.8,1.e10);
 AddTrackCuts(esdTrackCuts);

 const Int_t nptbins=4;
 const Int_t nvars=12;
 Float_t* ptbins;
 ptbins=new Float_t[nptbins+1];
 
 ptbins[0]=0.;
 ptbins[1]=2.;
 ptbins[2]=3.;
 ptbins[3]=4.;
 ptbins[4]=9999.;

 SetGlobalIndex(nvars,nptbins);
 SetPtBins(nptbins+1,ptbins);

 Float_t** prodcutsval;
 prodcutsval=new Float_t*[nvars];
 for(Int_t iv=0;iv<nvars;iv++){
  prodcutsval[iv]=new Float_t[nptbins];
 }

 for(Int_t ipt=0;ipt<nptbins;ipt++){
  prodcutsval[0][ipt]=0.13;
  prodcutsval[1][ipt]=0.5;
  prodcutsval[2][ipt]=0.6;
  prodcutsval[3][ipt]=0.;
  prodcutsval[4][ipt]=0.;
  prodcutsval[5][ipt]=0.01;
  prodcutsval[6][ipt]=0.04;
  prodcutsval[7][ipt]=0.006;
  prodcutsval[8][ipt]=0.8;
  prodcutsval[9][ipt]=0.3;
  prodcutsval[10][ipt]=0.;
  prodcutsval[11][ipt]=0.05;
 }
 SetCuts(nvars,nptbins,prodcutsval);

 AliAODPidHF* pidObjK=new AliAODPidHF();
 Double_t sigmasK[5]={3.,1.,1.,3.,2.};
 pidObjK->SetSigma(sigmasK);
 pidObjK->SetAsym(kTRUE);
 pidObjK->SetMatch(1);
 pidObjK->SetTPC(kTRUE);
 pidObjK->SetTOF(kTRUE);
 pidObjK->SetITS(kTRUE);
 Double_t plimK[2]={0.5,0.8};
 pidObjK->SetPLimit(plimK,2);

 SetPidHF(pidObjK);

 AliAODPidHF* pidObjpi=new AliAODPidHF();
 pidObjpi->SetTPC(kTRUE);
 Double_t sigmaspi[5]={3.,0.,0.,0.,0.};
 pidObjpi->SetSigma(sigmaspi);
 SetPidpion(pidObjpi);

 AliAODPidHF* pidObjp=new AliAODPidHF();
 Double_t sigmasp[5]={3.,1.,1.,3.,2.};
 pidObjp->SetSigma(sigmasp);
 pidObjp->SetAsym(kTRUE);
 pidObjp->SetMatch(1);
 pidObjp->SetTPC(kTRUE);
 pidObjp->SetTOF(kTRUE);
 pidObjp->SetITS(kTRUE);
 Double_t plimp[2]={1.,2.};
 pidObjp->SetPLimit(plimp,2);

 SetPidprot(pidObjp);

 SetUsePID(kTRUE);

 PrintAll();

 for(Int_t iiv=0;iiv<nvars;iiv++){
  delete [] prodcutsval[iiv];
 }
 delete [] prodcutsval;
 prodcutsval=NULL;
 delete [] ptbins;
 ptbins=NULL;

 delete pidObjp;
 pidObjp=NULL;

 return;
}
//------------------
Bool_t AliRDHFCutsLctopKpi::ReconstructKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field) const{

 Int_t nprongs=d->GetNProngs();
 Int_t iprongs[nprongs];
 for(Int_t i=0;i<nprongs;i++) iprongs[i]=i;

 Double_t mass[2]={0.,0.};

 AliKFParticle *decay=d->ApplyVertexingKF(iprongs,nprongs,pdgs,kTRUE,field,mass);
 if(!decay) return kTRUE;
 AliESDVertex *vertexESD = new AliESDVertex(decay->Parameters(),
                                            decay->CovarianceMatrix(),
					    decay->GetChi2(),
					    nprongs);
 Double_t pos[3],cov[6],chi2perNDF;
 vertexESD->GetXYZ(pos);
 vertexESD->GetCovMatrix(cov);
 chi2perNDF = vertexESD->GetChi2toNDF();
 delete vertexESD; vertexESD=NULL;
 AliAODVertex *vertexAOD = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,nprongs);
 d->SetSecondaryVtx(vertexAOD);
 return kTRUE;
}

