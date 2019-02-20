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

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliRDHFCutsDstoKKpi);
/// \endcond


//--------------------------------------------------------------------------
AliRDHFCutsDstoKKpi::AliRDHFCutsDstoKKpi(const char* name) : 
  AliRDHFCuts(name),
  fCutOnResonances(kTRUE),
  fPidOption(0),
  fMaxPtStrongPid(0.),
  fMaxPStrongPidK(0.),
  fMaxPStrongPidpi(0.),
  fDistToMaxProb(0.01),
  fBayesThreshold(0.05),
  fWeightKKpi(1.),
  fWeightpiKK(1.),
  fPhiMassRef(1.019),
  fUseRefPhiMass(kFALSE),
  fUsed0MeasMinusExpCut(kFALSE),
  fMaxd0MeasMinusExp(0x0),
  fUsed0Cut(kFALSE),
  fMaxd0(0x0),
  fCheckK0star(kTRUE)
{
  //
  // Default Constructor
  //
  Int_t nvars=20;
  SetNVars(nvars);
  TString varNames[20]={"inv. mass [GeV]",
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
			"inv. mass (MKo*-MKpi) [GeV]",
			"Abs(CosineKpiPhiRFrame)^3",
			"CosPiDsLabFrame",
			"decLenXY [cm]"
			"NormdecLen",
			"NormdecLenXY [cm]",
			"cosThetaPointXY"};
			
  Bool_t isUpperCut[20]={kTRUE,
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
			 kTRUE,
			 kFALSE,
			 kTRUE,
			 kFALSE,
			 kFALSE,
			 kFALSE,
			 kFALSE};
  SetVarNames(20,varNames,isUpperCut);
  Bool_t forOpt[20]={kFALSE,
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
		     kTRUE,
		     kFALSE,
		     kFALSE,
		     kTRUE,
		     kTRUE,
		     kTRUE,
		     kTRUE};
		    
  SetVarsForOpt(11,forOpt);
  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);
  if(fPidHF)delete fPidHF;
  fPidHF=new AliAODPidHF();
  Double_t plim[2]={0.6,0.8};
  Double_t nsigma[5]={2.,1.,2.,3.,0.};
  
  fPidHF->SetPLimit(plim,2);
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
  AliRDHFCuts(source),
  fCutOnResonances(source.fCutOnResonances),
  fPidOption(source.fPidOption),
  fMaxPtStrongPid(source.fMaxPtStrongPid),
  fMaxPStrongPidK(source.fMaxPStrongPidK),
  fMaxPStrongPidpi(source.fMaxPStrongPidpi),
  fDistToMaxProb(source.fDistToMaxProb),
  fBayesThreshold(source.fBayesThreshold),
  fWeightKKpi(source.fWeightKKpi),
  fWeightpiKK(source.fWeightpiKK),
  fPhiMassRef(source.fPhiMassRef),
  fUseRefPhiMass(source.fUseRefPhiMass),
  fUsed0MeasMinusExpCut(source.fUsed0MeasMinusExpCut),
  fMaxd0MeasMinusExp(0x0),
  fUsed0Cut(source.fUsed0Cut),
  fMaxd0(0x0),
  fCheckK0star(source.fCheckK0star)
{
  //
  // Copy constructor
  //
  if(source.fMaxd0MeasMinusExp) Setd0MeasMinusExpCut(source.fnPtBins,source.fMaxd0MeasMinusExp);
  if(source.fMaxd0) Setd0Cut(source.fnPtBins,source.fMaxd0);
}
//--------------------------------------------------------------------------
AliRDHFCutsDstoKKpi &AliRDHFCutsDstoKKpi::operator=(const AliRDHFCutsDstoKKpi &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliRDHFCuts::operator=(source);
  
  fCutOnResonances=source.fCutOnResonances;
  fPidOption=source.fPidOption;
  fMaxPtStrongPid=source.fMaxPtStrongPid;
  fMaxPStrongPidK=source.fMaxPStrongPidK;
  fMaxPStrongPidpi=source.fMaxPStrongPidpi;
  fDistToMaxProb=source.fDistToMaxProb;
  fBayesThreshold=source.fBayesThreshold;
  fWeightKKpi=source.fWeightKKpi;
  fWeightpiKK=source.fWeightpiKK;
  fPhiMassRef=source.fPhiMassRef;
  fUseRefPhiMass=source.fUseRefPhiMass;
  fUsed0MeasMinusExpCut=source.fUsed0MeasMinusExpCut;
  fUsed0Cut=source.fUsed0Cut;
  if(source.fMaxd0MeasMinusExp) Setd0MeasMinusExpCut(source.fnPtBins,source.fMaxd0MeasMinusExp);
  if(source.fMaxd0) Setd0Cut(source.fnPtBins,source.fMaxd0);
  if(source.fCheckK0star) fCheckK0star=source.fCheckK0star;

  return *this;
}
//---------------------------------------------------------------------------
AliRDHFCutsDstoKKpi::~AliRDHFCutsDstoKKpi(){
  //
  // Destructor
  //
  if(fMaxd0MeasMinusExp) delete [] fMaxd0MeasMinusExp;
  if(fMaxd0) delete [] fMaxd0;
}
//---------------------------------------------------------------------------
void AliRDHFCutsDstoKKpi::Setd0MeasMinusExpCut(Int_t nPtBins, Float_t *cutval) {
  //
  // store the cuts
  //
  if(nPtBins!=fnPtBins) {
    printf("Wrong number of pt bins: it has to be %d\n",fnPtBins);
    AliFatal("exiting");
  }
  if(!fMaxd0MeasMinusExp)  fMaxd0MeasMinusExp = new Float_t[fnPtBins];
  for(Int_t ib=0; ib<fnPtBins; ib++) fMaxd0MeasMinusExp[ib] = cutval[ib];
  fUsed0MeasMinusExpCut=kTRUE;
  return;
}

//---------------------------------------------------------------------------
void AliRDHFCutsDstoKKpi::Setd0Cut(Int_t nPtBins, Float_t *cutval) {
  //
  // store the cuts
  //
  if(nPtBins!=fnPtBins) {
    printf("Wrong number of pt bins: it has to be %d\n",fnPtBins);
    AliFatal("exiting");
  }
  if(!fMaxd0)  fMaxd0 = new Float_t[fnPtBins];
  for(Int_t ib=0; ib<fnPtBins; ib++) fMaxd0[ib] = cutval[ib];
  fUsed0Cut=kTRUE;
  return;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsDstoKKpi::PreSelect(TObjArray aodTracks){
  //
  // apply pre selection
  //
  if(!fUsePreselect) return 3;
  Int_t retVal=3;
  
  //compute pt
  Double_t px=0, py=0;
  AliAODTrack *track[3];
  for(Int_t iDaught=0; iDaught<3; iDaught++) {
    track[iDaught] = (AliAODTrack*)aodTracks.At(iDaught);
    if(!track[iDaught]) return retVal;
    px += track[iDaught]->Px();
    py += track[iDaught]->Py();
  }
  
  Double_t ptD=TMath::Sqrt(px*px+py*py);
  
  Int_t pidoptmem = fPidOption;
  fPidOption=kConservative;
  retVal=IsSelectedPID(ptD,aodTracks);
  fPidOption=pidoptmem;
  
  if(fUsePreselect==1) return retVal;
  
  Int_t ptbin=PtBin(ptD);
  Double_t mPhiPDG = TDatabasePDG::Instance()->GetParticle(333)->Mass();
  Double_t mK0starPDG = TDatabasePDG::Instance()->GetParticle(313)->Mass();

  Bool_t okDsKKpi = (retVal==1 || retVal==3) ? kTRUE : kFALSE;
  Bool_t okDspiKK = (retVal==2 || retVal==3) ? kTRUE : kFALSE;

  if(okDsKKpi){
    //compute min and max with 20% tolerance
    Double_t minmassphi2=(mPhiPDG-fCutsRD[GetGlobalIndex(12,ptbin)]*1.2)*(mPhiPDG-fCutsRD[GetGlobalIndex(12,ptbin)]*1.2);
    Double_t maxmassphi2=(mPhiPDG+fCutsRD[GetGlobalIndex(12,ptbin)]*1.2)*(mPhiPDG+fCutsRD[GetGlobalIndex(12,ptbin)]*1.2);
    
    Double_t mass01phi2=ComputeInvMass2(track[0],track[1],321,321);
    
    if(mass01phi2<minmassphi2 || mass01phi2>maxmassphi2) okDsKKpi=kFALSE;
    
    if(!okDsKKpi && fCheckK0star){
      //compute min and max with 20% tolerance
      Double_t minmassK0s2=(mK0starPDG-fCutsRD[GetGlobalIndex(13,ptbin)]*1.2)*(mK0starPDG-fCutsRD[GetGlobalIndex(13,ptbin)]*1.2);
      Double_t maxmassK0s2=(mK0starPDG+fCutsRD[GetGlobalIndex(13,ptbin)]*1.2)*(mK0starPDG+fCutsRD[GetGlobalIndex(13,ptbin)]*1.2);
      
      Double_t mass12K0s=ComputeInvMass2(track[1],track[2],321,211);
      
      if(mass12K0s<minmassK0s2 || mass12K0s>maxmassK0s2) okDsKKpi=kFALSE;
    }
  }
  if(okDspiKK){
    //compute min and max with 20% tolerance
    Double_t minmassphi2=(mPhiPDG-fCutsRD[GetGlobalIndex(12,ptbin)]*1.2)*(mPhiPDG-fCutsRD[GetGlobalIndex(12,ptbin)]*1.2);
    Double_t maxmassphi2=(mPhiPDG+fCutsRD[GetGlobalIndex(12,ptbin)]*1.2)*(mPhiPDG+fCutsRD[GetGlobalIndex(12,ptbin)]*1.2);
    
    Double_t mass01phi2=ComputeInvMass2(track[1],track[2],321,321);
    
    if(mass01phi2<minmassphi2 || mass01phi2>maxmassphi2) okDsKKpi=kFALSE;
    
    if(!okDsKKpi && fCheckK0star){
      //compute min and max with 20% tolerance
      Double_t minmassK0s2=(mK0starPDG-fCutsRD[GetGlobalIndex(13,ptbin)]*1.2)*(mK0starPDG-fCutsRD[GetGlobalIndex(13,ptbin)]*1.2);
      Double_t maxmassK0s2=(mK0starPDG+fCutsRD[GetGlobalIndex(13,ptbin)]*1.2)*(mK0starPDG+fCutsRD[GetGlobalIndex(13,ptbin)]*1.2);
      
      Double_t mass12K0s=ComputeInvMass2(track[0],track[1],211,321);
      
      if(mass12K0s<minmassK0s2 || mass12K0s>maxmassK0s2) okDsKKpi=kFALSE;
    }
  }
  
  if(okDsKKpi && okDspiKK) retVal=3;
  else if(okDsKKpi && !okDspiKK) retVal=1;
  else if(!okDsKKpi && okDspiKK) retVal=2;

  return retVal;
}

//---------------------------------------------------------------------------
void AliRDHFCutsDstoKKpi::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters,AliAODEvent *aod) {
  // 
  // Fills in vars the values of the variables 
  //

  if(nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsDstoKKpi::GetCutsVarsForOpt: wrong number of variables\n");
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
    if(TMath::Abs(pdgdaughters[0])==321){
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
    Float_t maxDCA=0.;
    for(Int_t i=0;i<3;i++){ 
      if(d->GetDCA(i)>maxDCA) maxDCA=d->GetDCA(i);
    }
    vars[iter]=maxDCA;
  }
  if(fVarsForOpt[12]){
    iter++;
    Double_t mPDGPhi = TDatabasePDG::Instance()->GetParticle(333)->Mass();
    if(TMath::Abs(pdgdaughters[0])==321){
      
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
    if(TMath::Abs(pdgdaughters[0])==321){
      
      Double_t mass12kpi=d->InvMass2Prongs(1,2,321,211);
      vars[iter]=TMath::Abs(mass12kpi-mPDGK0star);
      //	      vars[iter]=dd->InvMass2Prongs(1,2,321,211);
    }else{
      Double_t mass01pik=d->InvMass2Prongs(0,1,211,321);
      vars[iter]=TMath::Abs(mass01pik-mPDGK0star);
      //	vars[iter]=dd->InvMass2Prongs(0,1,211,321);      
    }
  }
  if(fVarsForOpt[14]){
    iter++;
    if(TMath::Abs(pdgdaughters[0])==321){
      vars[iter]=dd->CosPiKPhiRFrameKKpi();
    }else{
      vars[iter]=dd->CosPiKPhiRFramepiKK();
    }
  }
  if(fVarsForOpt[15]){
    iter++;
    if(TMath::Abs(pdgdaughters[0])==321){
      vars[iter]=dd->CosPiDsLabFrameKKpi();
    }else{
      vars[iter]=dd->CosPiDsLabFramepiKK();
    }
  }
  
  if(fVarsForOpt[16]){
    iter++;
    vars[iter]=dd->DecayLengthXY();
  }
  
  if(fVarsForOpt[17]){
    iter++;
    vars[iter]=dd->NormalizedDecayLength();
  }
  
  if(fVarsForOpt[18]){
    iter++;
    vars[iter]=dd->NormalizedDecayLengthXY();
  }
  
  if(fVarsForOpt[19]){
    iter++;
    vars[iter]=dd->CosPointingAngleXY();
  }

  if(cleanvtx)CleanOwnPrimaryVtx(dd,aod,origownvtx); 
  return;
}
//---------------------------------------------------------------------------
Bool_t AliRDHFCutsDstoKKpi::IsInFiducialAcceptance(Double_t pt, Double_t y) const
{
  //
  // Checking if Ds is in fiducial acceptance region 
  //

  if(fMaxRapidityCand>-998.){
    if(TMath::Abs(y) > fMaxRapidityCand) return kFALSE;
    else return kTRUE;
  }

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
Int_t AliRDHFCutsDstoKKpi::IsSelectedPIDBayes(AliAODRecoDecayHF *rd) {
  Int_t retCode=3;
  if(rd) {
    Double_t Pt = rd->Pt();
    TObjArray aodtracks(3);
    for(Int_t iDaught=0; iDaught<3; iDaught++) {
      aodtracks.AddAt(rd->GetDaughter(iDaught),iDaught);
    }
    retCode = IsSelectedPIDBayes(Pt,aodtracks);
  }
  
  return retCode;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsDstoKKpi::IsSelectedPIDBayes(Double_t Pt, TObjArray aodtracks) {
  Int_t retCode=3;
  Bool_t okKKpi=kTRUE;
  Bool_t okpiKK=kTRUE;
  
  AliAODTrack *track[3];
  for(Int_t iDaught=0; iDaught<3; iDaught++){
    track[iDaught]=(AliAODTrack*)aodtracks.At(iDaught);
  }
  
  if(!fUsePID || !track[0] || !track[1] || !track[2]) return retCode;
  
  Int_t sign = track[0]->Charge();
  if(!fPidHF){
    AliWarning("AliAODPidHF not created!");
    return retCode;
  }
  if(fPidOption!=kBayesianMaxProb && fPidOption!=kBayesianThreshold && fPidOption!=kBayesianWeights){
    AliWarning("Wrong call to Bayesian PID");
    return retCode;
  }
  
  AliPIDCombined* copid=fPidHF->GetPidCombined();
  copid->SetDetectorMask(AliPIDResponse::kDetTPC | AliPIDResponse::kDetTOF);
  AliPIDResponse* pidres=fPidHF->GetPidResponse();
  Double_t bayesProb[AliPID::kSPECIES];
  Int_t nKaons=0;
  Int_t nNotKaons=0;
  fWeightKKpi=1.;
  fWeightpiKK=1.;
  for(Int_t iDaught=0; iDaught<3; iDaught++){
    
    Int_t isPion=0;
    Int_t isKaon=0;
    Int_t isProton=0;
    UInt_t usedDet=copid->ComputeProbabilities(track[iDaught],pidres,bayesProb);
    if(usedDet!=0){
      if(fPidOption==kBayesianMaxProb){
        Double_t maxProb=TMath::MaxElement(AliPID::kSPECIES,bayesProb);
        if(TMath::Abs(maxProb-bayesProb[AliPID::kPion])<fDistToMaxProb) isPion=1;
        else isPion=-1;
        if(TMath::Abs(maxProb-bayesProb[AliPID::kKaon])<fDistToMaxProb) isKaon=1;
        else isKaon=-1;
        if(TMath::Abs(maxProb-bayesProb[AliPID::kProton])<fDistToMaxProb) isProton=1;
        else isProton=-1;
      }
      if(fPidOption==kBayesianThreshold){
        if(bayesProb[AliPID::kPion]>fBayesThreshold) isPion=1;
        else isPion=-1;
        if(bayesProb[AliPID::kKaon]>fBayesThreshold) isKaon=1;
        else isKaon=-1;
        if(bayesProb[AliPID::kProton]>fBayesThreshold) isProton=1;
        else isProton=-1;
      }
    }
    if(fPidOption==kBayesianWeights){ // store the probabilities in the case kBayesianWeights
      if(iDaught==0){
        fWeightKKpi*=bayesProb[AliPID::kKaon];
        fWeightpiKK*=bayesProb[AliPID::kPion];
      }else if(iDaught==1){
        fWeightKKpi*=bayesProb[AliPID::kKaon];
        fWeightpiKK*=bayesProb[AliPID::kKaon];
      }else if(iDaught==2){
        fWeightKKpi*=bayesProb[AliPID::kPion];
        fWeightpiKK*=bayesProb[AliPID::kKaon];
      }
    }else{ // selection for the other cases
      if(isProton>0 &&  isKaon<0  && isPion<0) return 0;
      if(sign!=track[iDaught]->Charge()){// must be kaon
        if(isKaon<0) return 0;
      }
      if(isKaon>0 && isPion<0) nKaons++;
      if(isKaon<0) nNotKaons++;
      if(iDaught==0){
        if(isKaon<0) okKKpi=kFALSE;
        if(isPion<0) okpiKK=kFALSE;
      }else if(iDaught==2){
        if(isKaon<0) okpiKK=kFALSE;
        if(isPion<0) okKKpi=kFALSE;
      }
    }
  }
  if(fPidOption==kBayesianWeights) return retCode;
  
  if(nKaons>2)return 0;
  if(nNotKaons>1) return 0;
  
  if(!okKKpi) retCode-=1;
  if(!okpiKK) retCode-=2;
  
  return retCode;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsDstoKKpi::IsSelectedPID(AliAODRecoDecayHF *rd) {
  // PID selection
  // return values: 0->NOT OK, 1->OK as KKpi, 2->OK as piKK, 3->OK as both
  Int_t retCode=3;
  if(!fUsePID) return retCode;
  if(rd) {
    Double_t Pt = rd->Pt();
    TObjArray aodtracks(3);
    for(Int_t iDaught=0; iDaught<3; iDaught++) {
      aodtracks.AddAt(rd->GetDaughter(iDaught),iDaught);
    }
    retCode = IsSelectedPID(Pt,aodtracks);
  }
  
  return retCode;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsDstoKKpi::IsSelectedPID(Double_t Pt, TObjArray aodtracks) {
  // PID selection
  // return values: 0->NOT OK, 1->OK as KKpi, 2->OK as piKK, 3->OK as both
  
  Int_t retCode=3;
  Bool_t okKKpi=kTRUE;
  Bool_t okpiKK=kTRUE;
  
  AliAODTrack *track[3];
  for(Int_t iDaught=0; iDaught<3; iDaught++){
    track[iDaught]=(AliAODTrack*)aodtracks.At(iDaught);
  }
  
  if(!fUsePID || !track[0] || !track[1] || !track[2]) return retCode;

  Int_t sign = track[0]->Charge();
  if(!fPidHF){
    AliWarning("AliAODPidHF not created!");
    return retCode;
  }
  if(fPidOption==kBayesianMaxProb || fPidOption==kBayesianThreshold || fPidOption==kBayesianWeights){
    // call method for Bayesian probability
    return IsSelectedPIDBayes(Pt,aodtracks);
  }
  
  Double_t origCompatTOF=fPidHF->GetPCompatTOF();
  Double_t origThreshTPC=fPidHF->GetPtThresholdTPC();
  if(fPidOption==kStrong){
    fPidHF->SetPCompatTOF(999999.);
    fPidHF->SetPtThresholdTPC(999999.);
  }
  
  Int_t nKaons=0;
  Int_t nNotKaons=0;
  for(Int_t iDaught=0; iDaught<3; iDaught++){
    
    Int_t isPion=fPidHF->MakeRawPid(track[iDaught],AliPID::kPion);
    Int_t isKaon=fPidHF->MakeRawPid(track[iDaught],AliPID::kKaon);
    Int_t isProton=fPidHF->MakeRawPid(track[iDaught],AliPID::kProton);
    
    if(isProton>0 &&  isKaon<0  && isPion<0){
      fPidHF->SetPCompatTOF(origCompatTOF);
      fPidHF->SetPtThresholdTPC(origThreshTPC);
      return 0;
    }
    if(sign!=track[iDaught]->Charge()){// must be kaon
      if(isKaon<0){
        fPidHF->SetPCompatTOF(origCompatTOF);
        fPidHF->SetPtThresholdTPC(origThreshTPC);
        return 0;
      }
      if(fPidOption==kStrong && Pt<fMaxPtStrongPid && isKaon<=0){
        fPidHF->SetPCompatTOF(origCompatTOF);
        fPidHF->SetPtThresholdTPC(origThreshTPC);
        return 0;
      }
      if(fPidOption==kStrongPDep && Pt<fMaxPtStrongPid){
        if(isKaon<=0 && track[iDaught]->P()<fMaxPStrongPidK) return 0;
      }
    }
    
    if(isKaon>0 && isPion<0) nKaons++;
    if(isKaon<0) nNotKaons++;
    if(iDaught==0){
      if(isKaon<0) okKKpi=kFALSE;
      if(isPion<0) okpiKK=kFALSE;
      if(fPidOption==kStrong && Pt<fMaxPtStrongPid){
        if(isKaon<=0) okKKpi=kFALSE;
        if(isPion<=0) okpiKK=kFALSE;
      }
      if(fPidOption==kStrongPDep && Pt<fMaxPtStrongPid){
        if(isKaon<=0 && track[iDaught]->P()<fMaxPStrongPidK) okKKpi=kFALSE;
        if(isPion<=0 && track[iDaught]->P()<fMaxPStrongPidpi) okpiKK=kFALSE;
      }
    }
    else if(iDaught==2){
      if(isKaon<0) okpiKK=kFALSE;
      if(isPion<0) okKKpi=kFALSE;
      if(fPidOption==kStrong && Pt<fMaxPtStrongPid){
        if(isKaon<=0) okpiKK=kFALSE;
        if(isPion<=0) okKKpi=kFALSE;
      }
      if(fPidOption==kStrongPDep && Pt<fMaxPtStrongPid){
        if(isKaon<=0 && track[iDaught]->P()<fMaxPStrongPidK) okpiKK=kFALSE;
        if(isPion<=0 && track[iDaught]->P()<fMaxPStrongPidpi) okKKpi=kFALSE;
      }
    }
  }
  
  fPidHF->SetPCompatTOF(origCompatTOF);
  fPidHF->SetPtThresholdTPC(origThreshTPC);
  
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
  
  if(fKeepSignalMC) if(IsSignalMC(d,aod,431)) return 3;
 
  Double_t ptD=d->Pt();
  if(ptD<fMinPtCand) return 0;
  if(ptD>fMaxPtCand) return 0;

  if(fUseTrackSelectionWithFilterBits && d->HasBadDaughters()) return 0;
  

  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(d,aod)) return 0;
  }



 
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

    Int_t okDsKKpi=1;
    Int_t okDspiKK=1;
    Int_t okMassPhiKKpi=0;
    Int_t okMassPhipiKK=0;
    Int_t okMassK0starKKpi=0;
    Int_t okMassK0starpiKK=0;
    Int_t okDsPhiKKpi=0;
    Int_t okDsPhipiKK=0;
    Int_t okDsK0starKKpi=0;
    Int_t okDsK0starpiKK=0;

    Double_t pt=d->Pt();
    Int_t ptbin=PtBin(pt);
    if (ptbin==-1) {
      CleanOwnPrimaryVtx(d,aod,origownvtx);
      return 0;
    }
 
    Double_t mDsPDG = TDatabasePDG::Instance()->GetParticle(431)->Mass();
    Double_t mDsKKpi=d->InvMassDsKKpi();
    Double_t mDspiKK=d->InvMassDspiKK();
    if(TMath::Abs(mDsKKpi-mDsPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) okDsKKpi = 0;
    if(TMath::Abs(mDspiKK-mDsPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) okDspiKK = 0;
    if(!okDsKKpi && !okDspiKK){
      CleanOwnPrimaryVtx(d,aod,origownvtx);
      return 0;
    }



    // cuts on resonant decays (via Phi or K0*)
    if(fCutOnResonances){
      Double_t mPhiPDG = TDatabasePDG::Instance()->GetParticle(333)->Mass();
      Double_t mK0starPDG = TDatabasePDG::Instance()->GetParticle(313)->Mass();
      Double_t mPhiRef = mPhiPDG;
      if(fUseRefPhiMass) mPhiRef = fPhiMassRef;
      if(okDsKKpi){
	Double_t mass01phi=d->InvMass2Prongs(0,1,321,321);
	if(TMath::Abs(mass01phi-mPhiRef)<fCutsRD[GetGlobalIndex(12,ptbin)]) okMassPhiKKpi=1;
	if(fCheckK0star){
	  Double_t mass12K0s=d->InvMass2Prongs(1,2,321,211);
	  if(TMath::Abs(mass12K0s-mK0starPDG)<fCutsRD[GetGlobalIndex(13,ptbin)]) okMassK0starKKpi = 1;
	}
	if(!okMassPhiKKpi && !okMassK0starKKpi) okDsKKpi=0;
	if(okMassPhiKKpi) okDsPhiKKpi=1;
	if(okMassK0starKKpi) okDsK0starKKpi=1;
      }
      if(okDspiKK){
	if(fCheckK0star){
	  Double_t mass01K0s=d->InvMass2Prongs(0,1,211,321);
	  if(TMath::Abs(mass01K0s-mK0starPDG)<fCutsRD[GetGlobalIndex(13,ptbin)]) okMassK0starpiKK = 1;
	}
	Double_t mass12phi=d->InvMass2Prongs(1,2,321,321);
	if(TMath::Abs(mass12phi-mPhiRef)<fCutsRD[GetGlobalIndex(12,ptbin)]) okMassPhipiKK=1;
	if(!okMassPhipiKK && !okMassK0starpiKK) okDspiKK=0;
	if(okMassPhipiKK) okDsPhipiKK=1;
	if(okMassK0starpiKK) okDsK0starpiKK=1;
      }
      if(!okDsKKpi && !okDspiKK){
	CleanOwnPrimaryVtx(d,aod,origownvtx);
	return 0;
      }
    }

    // Cuts on track pairs
    for(Int_t i=0;i<3;i++){
      if(d->GetDCA(i)>fCutsRD[GetGlobalIndex(11,ptbin)]){
	CleanOwnPrimaryVtx(d,aod,origownvtx);
	return 0;
      }
    }
    if(d->GetDist12toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)] || 
       d->GetDist23toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]){
      CleanOwnPrimaryVtx(d,aod,origownvtx);
      return 0;
    }



    //single track
    if(TMath::Abs(d->Pt2Prong(1)) < fCutsRD[GetGlobalIndex(1,ptbin)]*fCutsRD[GetGlobalIndex(1,ptbin)] || 
       TMath::Abs(d->Getd0Prong(1))<fCutsRD[GetGlobalIndex(3,ptbin)]){
      CleanOwnPrimaryVtx(d,aod,origownvtx);
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
      CleanOwnPrimaryVtx(d,aod,origownvtx);
      return 0;
    }

    // Cuts on candidate triplet


    if(d->CosPointingAngle()< fCutsRD[GetGlobalIndex(9,ptbin)]){
      CleanOwnPrimaryVtx(d,aod,origownvtx); 
      return 0;
    }
     
    if(d->Pt2Prong(0)<fCutsRD[GetGlobalIndex(8,ptbin)]*fCutsRD[GetGlobalIndex(8,ptbin)] && 
       d->Pt2Prong(1)<fCutsRD[GetGlobalIndex(8,ptbin)]*fCutsRD[GetGlobalIndex(8,ptbin)] && 
       d->Pt2Prong(2)<fCutsRD[GetGlobalIndex(8,ptbin)]*fCutsRD[GetGlobalIndex(8,ptbin)]) {
      CleanOwnPrimaryVtx(d,aod,origownvtx); 
      return 0;
    }

    if(d->DecayLength2()<fCutsRD[GetGlobalIndex(7,ptbin)]*fCutsRD[GetGlobalIndex(7,ptbin)]){
      CleanOwnPrimaryVtx(d,aod,origownvtx);
      return 0;
    }


    Double_t sum2=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
    if(sum2<fCutsRD[GetGlobalIndex(10,ptbin)]){
      CleanOwnPrimaryVtx(d,aod,origownvtx);
      return 0;
    }

   
    //sec vert
    Double_t sigmavert=d->GetSigmaVert(aod);
    if(sigmavert>fCutsRD[GetGlobalIndex(6,ptbin)]){
      CleanOwnPrimaryVtx(d,aod,origownvtx);
      return 0;
    }
    
    // decay length XY
    if(d->DecayLengthXY()<fCutsRD[GetGlobalIndex(16,ptbin)]){
      CleanOwnPrimaryVtx(d,aod,origownvtx); 
      return 0;
    }
    
    //norm decay length
    if(d->NormalizedDecayLength()<fCutsRD[GetGlobalIndex(17,ptbin)]){
      CleanOwnPrimaryVtx(d,aod,origownvtx); 
      return 0;
    }
    
    //norm decay length XY
    if(d->NormalizedDecayLengthXY()<fCutsRD[GetGlobalIndex(18,ptbin)]){
      CleanOwnPrimaryVtx(d,aod,origownvtx); 
      return 0;
    }
    
    //cos pointing XY
    if(d->CosPointingAngleXY()<fCutsRD[GetGlobalIndex(19,ptbin)]){
      CleanOwnPrimaryVtx(d,aod,origownvtx); 
      return 0;
    }
    

    if(okDsKKpi){
      Double_t cosPiKPhiRFKKpi=d->CosPiKPhiRFrameKKpi();
      Double_t kincutPiKPhiKKpi=TMath::Abs(cosPiKPhiRFKKpi*cosPiKPhiRFKKpi*cosPiKPhiRFKKpi);
      if(kincutPiKPhiKKpi<fCutsRD[GetGlobalIndex(14,ptbin)]) okDsKKpi=0;
    }
    if(okDspiKK){
      Double_t cosPiKPhiRFpiKK=d->CosPiKPhiRFramepiKK();
      Double_t kincutPiKPhipiKK=TMath::Abs(cosPiKPhiRFpiKK*cosPiKPhiRFpiKK*cosPiKPhiRFpiKK);
      if(kincutPiKPhipiKK<fCutsRD[GetGlobalIndex(14,ptbin)]) okDspiKK=0;
    }
    if(!okDsKKpi && !okDspiKK){
      CleanOwnPrimaryVtx(d,aod,origownvtx);
      return 0;
    }
    
    
    
    if(okDsKKpi){
      Double_t cosPiDsLabFrameKKpi=d->CosPiDsLabFrameKKpi();
      if(cosPiDsLabFrameKKpi>fCutsRD[GetGlobalIndex(15,ptbin)]) okDsKKpi=0;
    }
    if(okDspiKK){
      Double_t cosPiDsLabFramepiKK=d->CosPiDsLabFramepiKK();
      if(cosPiDsLabFramepiKK>fCutsRD[GetGlobalIndex(15,ptbin)]) okDspiKK=0;
    }
    if(!okDsKKpi && !okDspiKK){
      CleanOwnPrimaryVtx(d,aod,origownvtx);
      return 0;
    }
    
    // d0meas-exp
    if(fUsed0MeasMinusExpCut){
      Double_t dd0max=0;
      for(Int_t ipr=0; ipr<3; ipr++) {
	Double_t diffIP, errdiffIP;
	d->Getd0MeasMinusExpProng(ipr,aod->GetMagneticField(),diffIP,errdiffIP);
	Double_t normdd0=0.;
	if(errdiffIP>0) normdd0=diffIP/errdiffIP;
	if(ipr==0) dd0max=normdd0;
	else if(TMath::Abs(normdd0)>TMath::Abs(dd0max)) dd0max=normdd0;
      }
      if(TMath::Abs(dd0max)>fMaxd0MeasMinusExp[ptbin]) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}
    }
      
    // d0
    if(fUsed0Cut){
      Double_t d0=d->ImpParXY()*10000.;
      if(TMath::Abs(d0)>fMaxd0[ptbin]) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}
    }
      

    // unset recalculated primary vertex when not needed any more
    CleanOwnPrimaryVtx(d,aod,origownvtx);
      
    

    if(!okDsKKpi){
      okDsPhiKKpi=0;
      okDsK0starKKpi=0;
    }
    if(!okDspiKK){
      okDsPhipiKK=0;
      okDsK0starpiKK=0;
    }

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
    if(!okPidDsKKpi){
      okDsPhiKKpi=0;
      okDsK0starKKpi=0;
    }
    if(!okPidDspiKK){
      okDsPhipiKK=0;
      okDsK0starpiKK=0;
    }

    if((okPidDsKKpi && okDsKKpi)||(okPidDspiKK && okDspiKK)){
      Int_t returnvalue=0;
      if(okDsKKpi) returnvalue+=1;
      if(okDspiKK) returnvalue+=2;
      if(okDsPhiKKpi) returnvalue+=4;
      if(okDsPhipiKK) returnvalue+=8;
      if(okDsK0starKKpi) returnvalue+=16;
      if(okDsK0starpiKK) returnvalue+=32;
      return returnvalue;
    }else{
      return 0;
    }
  }
  return 15;

}

//--------------------------------------------------------------------------

UInt_t AliRDHFCutsDstoKKpi::GetPIDTrackTPCTOFBitMap(AliAODTrack *track) const{

  UInt_t bitmap=0;

  Double_t sigmaTPCPionHyp=-999.;
  Double_t sigmaTPCKaonHyp=-999.;
  Double_t sigmaTPCProtonHyp=-999.;
  Double_t sigmaTOFPionHyp=-999.;
  Double_t sigmaTOFKaonHyp=-999.;
  Double_t sigmaTOFProtonHyp=-999.;
  
  Int_t oksigmaTPCPionHyp=fPidHF->GetnSigmaTPC(track,2,sigmaTPCPionHyp);
  Int_t oksigmaTPCKaonHyp=fPidHF->GetnSigmaTPC(track,3,sigmaTPCKaonHyp);
  Int_t oksigmaTPCProtonHyp=fPidHF->GetnSigmaTPC(track,4,sigmaTPCProtonHyp);
  Int_t oksigmaTOFPionHyp=fPidHF->GetnSigmaTOF(track,2,sigmaTOFPionHyp);
  Int_t oksigmaTOFKaonHyp=fPidHF->GetnSigmaTOF(track,3,sigmaTOFKaonHyp);
  Int_t oksigmaTOFProtonHyp=fPidHF->GetnSigmaTOF(track,4,sigmaTOFProtonHyp);
  
  sigmaTPCPionHyp=TMath::Abs(sigmaTPCPionHyp);
  sigmaTPCKaonHyp=TMath::Abs(sigmaTPCKaonHyp);
  sigmaTPCProtonHyp=TMath::Abs(sigmaTPCProtonHyp);
  sigmaTOFPionHyp=TMath::Abs(sigmaTOFPionHyp);
  sigmaTOFKaonHyp=TMath::Abs(sigmaTOFKaonHyp);
  sigmaTOFProtonHyp=TMath::Abs(sigmaTOFProtonHyp);

  if (oksigmaTPCPionHyp && sigmaTPCPionHyp>0.){
    if (sigmaTPCPionHyp<1.) bitmap+=1<<kTPCPionLess1;
    else{
      if (sigmaTPCPionHyp<2.) bitmap+=1<<kTPCPionMore1Less2;
      else { 
        if (sigmaTPCPionHyp<3.) bitmap+=1<<kTPCPionMore2Less3; 
        else bitmap+=1<<kTPCPionMore3;
      }
    }
  }
  
  if (oksigmaTPCKaonHyp && sigmaTPCKaonHyp>0.){
    if (sigmaTPCKaonHyp<1.) bitmap+=1<<kTPCKaonLess1;
    else{
      if (sigmaTPCKaonHyp<2.) bitmap+=1<<kTPCKaonMore1Less2;
      else { 
        if (sigmaTPCKaonHyp<3.) bitmap+=1<<kTPCKaonMore2Less3; 
        else bitmap+=1<<kTPCKaonMore3;
      }
    }
  }
  
  if (oksigmaTPCProtonHyp && sigmaTPCProtonHyp>0.){
    if (sigmaTPCProtonHyp<1.) bitmap+=1<<kTPCProtonLess1;
    else{
      if (sigmaTPCProtonHyp<2.) bitmap+=1<<kTPCProtonMore1Less2;
      else { 
        if (sigmaTPCProtonHyp<3.) bitmap+=1<<kTPCProtonMore2Less3; 
        else bitmap+=1<<kTPCProtonMore3;
      }
    }
  }
  
  if (oksigmaTOFPionHyp && sigmaTOFPionHyp>0.){
    if (sigmaTOFPionHyp<1.) bitmap+=1<<kTOFPionLess1;
    else{
      if (sigmaTOFPionHyp<2.) bitmap+=1<<kTOFPionMore1Less2;
      else { 
        if (sigmaTOFPionHyp<3.) bitmap+=1<<kTOFPionMore2Less3; 
        else bitmap+=1<<kTOFPionMore3;
      }
    }
  }
  
  if (oksigmaTOFKaonHyp && sigmaTOFKaonHyp>0.){
    if (sigmaTOFKaonHyp<1.) bitmap+=1<<kTOFKaonLess1;
    else{
      if (sigmaTOFKaonHyp<2.) bitmap+=1<<kTOFKaonMore1Less2;
      else { 
        if (sigmaTOFKaonHyp<3.) bitmap+=1<<kTOFKaonMore2Less3; 
        else bitmap+=1<<kTOFKaonMore3;
      }
    }
  }
  
  if (oksigmaTOFProtonHyp && sigmaTOFProtonHyp>0.){
    if (sigmaTOFProtonHyp<1.) bitmap+=1<<kTOFProtonLess1;
    else{
      if (sigmaTOFProtonHyp<2.) bitmap+=1<<kTOFProtonMore1Less2;
      else { 
        if (sigmaTOFProtonHyp<3.) bitmap+=1<<kTOFProtonMore2Less3; 
        else bitmap+=1<<kTOFProtonMore3;
      }
    }
  }
  
  
  
  return bitmap;

}


//---------------------------------------------------------------------------

void AliRDHFCutsDstoKKpi::SetStandardCutsPP2010() {
  //
  //STANDARD CUTS USED FOR 2010 pp analysis 
  //                                            
  
  SetName("DstoKKpiCutsStandard");
  SetTitle("Standard Cuts for D+s analysis");
  
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
  delete esdTrackCuts;
  esdTrackCuts=NULL;
 
   
  const Int_t nptbins=4;
  Float_t ptbins[nptbins+1];
  ptbins[0]=2.;
  ptbins[1]=4.;
  ptbins[2]=6.;
  ptbins[3]=8.;
  ptbins[4]=12.;
  
  const Int_t nvars=20;
      
  Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}  
  for(Int_t ipt=0;ipt<nptbins;ipt++){
  
    anacutsval[0][ipt]=0.35;
    anacutsval[1][ipt]=0.3;
    anacutsval[2][ipt]=0.3;
    anacutsval[3][ipt]=0.;
    anacutsval[4][ipt]=0.;
    anacutsval[5][ipt]=0.005;
    anacutsval[8][ipt]=0.;
    anacutsval[10][ipt]=0.;
    anacutsval[11][ipt]=1000.0;
    anacutsval[13][ipt]=0.1;
    anacutsval[16][ipt]=0.;
    anacutsval[17][ipt]=0.;
    anacutsval[18][ipt]=0.;
    anacutsval[19][ipt]=-1.;
    
    
  }
 
  //sigmavertex

  anacutsval[6][0]=0.020;
  anacutsval[6][1]=0.030;
  anacutsval[6][2]=0.030;
  anacutsval[6][3]=0.060;
   
  //Decay length
    
  anacutsval[7][0]=0.035;
  anacutsval[7][1]=0.035;
  anacutsval[7][2]=0.040;
  anacutsval[7][3]=0.040;
 
  //cosThetaPoint
    
  anacutsval[9][0]=0.94;
  anacutsval[9][1]=0.94;
  anacutsval[9][2]=0.94;
  anacutsval[9][3]=0.94;
   
  //phi DeltaMass
    
  anacutsval[12][0]=0.0080;
  anacutsval[12][1]=0.0050;
  anacutsval[12][2]=0.0045;
  anacutsval[12][3]=0.0090;
      
  //Kin1
    
  anacutsval[14][0]=0.10;
  anacutsval[14][1]=0.05;
  anacutsval[14][2]=0.0;
  anacutsval[14][3]=0.05;
      
  //Kin2
  
  anacutsval[15][0]=0.95;
  anacutsval[15][1]=0.95;
  anacutsval[15][2]=1.;
  anacutsval[15][3]=0.95;     
  
  fPidHF->SetOldPid(kTRUE);
  SetUsePID(kTRUE); 
  SetPidOption(1);
  SetMaxPtStrongPid(9999.);
  SetGlobalIndex(nvars,nptbins);
  SetPtBins(nptbins+1,ptbins);
  SetCuts(nvars,nptbins,anacutsval);
  SetRemoveDaughtersFromPrim(kTRUE);
  
  //  PrintAll();

  for(Int_t iic=0;iic<nvars;iic++){delete [] anacutsval[iic];}
  delete [] anacutsval;
  anacutsval=NULL;

  return;
}

Double_t AliRDHFCutsDstoKKpi::ComputeInvMass2(AliAODTrack* track1, AliAODTrack* track2, Int_t pdg1, Int_t pdg2) {
  
  Double_t mass1 = TDatabasePDG::Instance()->GetParticle(pdg1)->Mass();
  Double_t mass2 = TDatabasePDG::Instance()->GetParticle(pdg2)->Mass();

  Double_t px1 = track1->Px();
  Double_t py1 = track1->Py();
  Double_t pz1 = track1->Pz();
  Double_t px2 = track2->Px();
  Double_t py2 = track2->Py();
  Double_t pz2 = track2->Pz();
  
  Double_t E12 = TMath::Sqrt(mass1*mass1+px1*px1+py1*py1+pz1*pz1)+TMath::Sqrt(mass2*mass2+px2*px2+py2*py2+pz2*pz2);
  
  return E12*E12-((px1+px2)*(px1+px2)+(py1+py2)*(py1+py2)+(pz1+pz2)*(pz1+pz2));
}

void AliRDHFCutsDstoKKpi::PrintAll() const {
  //
  // print all cuts values
  // 
  AliRDHFCuts::PrintAll();
  if(fUsed0MeasMinusExpCut){
    printf("Cuts on d0meas-d0exp:\n");
    for(Int_t ib=0;ib<fnPtBins;ib++){
      printf("%f   ",fMaxd0MeasMinusExp[ib]);
    }
    printf("\n");
  }else{
    printf("No cut on d0meas-d0exp:\n");
  }
  if(fUsed0Cut){
    printf("Cuts on d0:\n");
    for(Int_t ib=0;ib<fnPtBins;ib++){
      printf("%f   ",fMaxd0[ib]);
    }
    printf("\n");
  }else{
    printf("No cut on d0\n");
  }
}
