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

/* $Id: AliRDHFCutsDStartoKpipi.cxx 61203 2013-03-02 22:52:17Z fprino $ */

/////////////////////////////////////////////////////////////
//
// Class for cuts on AOD reconstructed DStar->Kpipi
//
// Author: A.Grelli, alessandro.grelli@uu.nl
//
// PID method implemented by   Y.Wang, yifei@physi.uni-heidelberg.de
//           
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAODPid.h"
#include "AliTPCPIDResponse.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "TObjArray.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliRDHFCutsDStartoKpipi);
/// \endcond


//--------------------------------------------------------------------------
AliRDHFCutsDStartoKpipi::AliRDHFCutsDStartoKpipi(const char* name) : 
  AliRDHFCuts(name),
  fTrackCutsSoftPi(0),
  fMaxPtPid(9999.),
  fTPCflag(999.),
  fCircRadius(0.),
  fUseTPCtrackCutsOnD0Daughters(kTRUE),
  fUseTPCtrackCutsOnSoftPion(kFALSE)
{
  //
  // Default Constructor
  //
  
  Int_t nvars=16;
  SetNVars(nvars);
  TString varNames[16]={
    "inv. mass [GeV]",   
    "dca [cm]",
    "cosThetaStar", 
    "pTK [GeV/c]",
    "pTPi [GeV/c]",
    "d0K [cm]",
    "d0Pi [cm]",
    "d0d0 [cm^2]",
    "cosThetaPoint",
    "inv. mass half width of D* [GeV]",
    "half width of (M_Kpipi-M_D0) [GeV]",
    "PtMin of pi_s [GeV/c]",
    "PtMax of pi_s [GeV/c]",
    "theta, angle between the pi_s and decay plane of the D0 [rad]",
    "|cosThetaPointXY|", 
    "NormDecayLenghtXY"};
  Bool_t isUpperCut[16]={
    kTRUE,
    kTRUE,
    kTRUE,
    kFALSE,
    kFALSE,
    kTRUE,
    kTRUE,
    kTRUE,
    kFALSE,
    kTRUE,
    kTRUE,
    kTRUE,
    kTRUE,
    kFALSE,
    kFALSE, 
    kFALSE};
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[16]={
    kFALSE,
    kTRUE,
    kTRUE,
    kFALSE,
    kFALSE,
    kFALSE,
    kFALSE,
    kTRUE,
    kTRUE,
    kFALSE,
    kTRUE,
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
AliRDHFCutsDStartoKpipi::AliRDHFCutsDStartoKpipi(const AliRDHFCutsDStartoKpipi &source) :
  AliRDHFCuts(source),
  fTrackCutsSoftPi(0),
  fMaxPtPid(9999.),
  fTPCflag(999.),
  fCircRadius(0.),
  fUseTPCtrackCutsOnD0Daughters(source.fUseTPCtrackCutsOnD0Daughters),
  fUseTPCtrackCutsOnSoftPion(source.fUseTPCtrackCutsOnSoftPion)
{
  //
  // Copy constructor
  //
  
  if(source.GetTrackCutsSoftPi()) AddTrackCutsSoftPi(source.GetTrackCutsSoftPi());
  
}
//--------------------------------------------------------------------------
AliRDHFCutsDStartoKpipi &AliRDHFCutsDStartoKpipi::operator=(const AliRDHFCutsDStartoKpipi &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliRDHFCuts::operator=(source);
  if(source.GetTrackCutsSoftPi()) {
    delete fTrackCutsSoftPi;
    fTrackCutsSoftPi = new AliESDtrackCuts(*(source.GetTrackCutsSoftPi()));
  }

  return *this;
}
//---------------------------------------------------------------------------
AliRDHFCutsDStartoKpipi::~AliRDHFCutsDStartoKpipi(){
  //
  // Destructor
  //
  if (fTrackCutsSoftPi) { delete fTrackCutsSoftPi; fTrackCutsSoftPi = nullptr;}
}

//---------------------------------------------------------------------------
void AliRDHFCutsDStartoKpipi::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) {
  // 
  // Fills in vars the values of the variables 
  //
  if(nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsDStartoKpipi::GetCutsVarsForOpt: wrong number of variables\n");
    return;
  }
  
 
  AliAODRecoCascadeHF* dstarD0pi = (AliAODRecoCascadeHF*)d;

  AliAODTrack *softPi = (AliAODTrack*)dstarD0pi->GetBachelor();

  AliAODRecoDecayHF2Prong* dd = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();
  

  Int_t iter=-1;
  if(fVarsForOpt[0]){
    iter++;
    if(TMath::Abs(pdgdaughters[0])==211) {
      vars[iter]=dd->InvMassD0();
    } else {
      vars[iter]=dd->InvMassD0bar();
    }
  }
  if(fVarsForOpt[1]){
    iter++;
    vars[iter]=dd->GetDCA();
  }
  if(fVarsForOpt[2]){
    iter++;
    if(TMath::Abs(pdgdaughters[0])==211) {
      vars[iter] = dd->CosThetaStarD0();
    } else {
      vars[iter] = dd->CosThetaStarD0bar();
    }
  }
  if(fVarsForOpt[3]){
    iter++;
   if(TMath::Abs(pdgdaughters[0])==321) {
     vars[iter]=dd->PtProng(0);
   }
   else{
     vars[iter]=dd->PtProng(1);
   }
  }
  if(fVarsForOpt[4]){
    iter++;
   if(TMath::Abs(pdgdaughters[0])==211) {
     vars[iter]=dd->PtProng(0);
   }
   else{
     vars[iter]=dd->PtProng(1);
   }
  }
  if(fVarsForOpt[5]){
    iter++;
    if(TMath::Abs(pdgdaughters[0])==321) {
     vars[iter]=dd->Getd0Prong(0);
   }
   else{
     vars[iter]=dd->Getd0Prong(1);
   }
  }
  if(fVarsForOpt[6]){
    iter++;
     if(TMath::Abs(pdgdaughters[0])==211) {
     vars[iter]=dd->Getd0Prong(0);
   }
   else{
     vars[iter]=dd->Getd0Prong(1);
   }
  }
  if(fVarsForOpt[7]){
    iter++;
    vars[iter]= dd->Prodd0d0();
  }
  if(fVarsForOpt[8]){
    iter++;
    vars[iter]=dd->CosPointingAngle();
  }
  if(fVarsForOpt[9]){
    iter++;
    vars[iter]=dstarD0pi->InvMassDstarKpipi();
  }
  if(fVarsForOpt[10]){
    iter++;
    vars[iter]=dstarD0pi->DeltaInvMass();
  }
  if(fVarsForOpt[11]){
    iter++;
    vars[iter] = softPi->Pt();
  }
  if(fVarsForOpt[12]){
    iter++;
    vars[iter] = softPi->Pt();
  }
  if(fVarsForOpt[13]){
    iter++;
    vars[iter] =dstarD0pi->AngleD0dkpPisoft();
  }
  if(fVarsForOpt[14]){
    iter++;
    vars[iter]=TMath::Abs(dd->CosPointingAngleXY());
  }
  if(fVarsForOpt[15]){
    iter++;
    vars[iter]=(dd->NormalizedDecayLengthXY()*(dd->P()/dd->Pt()));
  }
 
  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsDStartoKpipi::PreSelect(TObjArray aodTracks){
  //
  // apply pre selection
  //
  if(!fUsePreselect)return 3;
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
  
  if(ptD<fMinPtCand) return 0;
  if(ptD>fMaxPtCand) return 0;

  Int_t ptbin=PtBin(ptD);
  if (ptbin==-1) {
    return 0;
  }
  Float_t xy[2],z[2];
  track[1]->GetImpactParameters(xy[0],z[0]);
  track[2]->GetImpactParameters(xy[1],z[1]);
  if(xy[0]*xy[1] > fCutsRD[GetGlobalIndex(7,ptbin)]) return 0;

  retVal=IsSelectedPID(ptD,aodTracks);

  return retVal;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsDStartoKpipi::IsSelected(TObject* obj,Int_t selectionLevel, AliAODEvent* aod) {
  //
  // Apply selection for D*.
  // Added functionality to remove the D0 daughters from primary vertex (not dafult) 

  fIsSelectedCuts=0;
  fIsSelectedPID=0;

  if(!fCutsRD){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  
  AliAODRecoCascadeHF* d = (AliAODRecoCascadeHF*)obj;
  if(!d){
    cout<<"AliAODRecoCascadeHF null"<<endl;
    return 0;
  }
  
  Double_t ptD=d->Pt();
  if(ptD<fMinPtCand) return 0;
  if(ptD>fMaxPtCand) return 0;
  

  AliAODRecoDecayHF2Prong* dd = (AliAODRecoDecayHF2Prong*)d->Get2Prong();  
  if(!dd){
    cout<<"AliAODRecoDecayHF2Prong null"<<endl;
    return 0;
  }

  if(fUseTrackSelectionWithFilterBits && dd->HasBadDaughters()) return 0;

  AliAODTrack *b = (AliAODTrack*)d->GetBachelor();
  if(fTrackCutsSoftPi && fTrackCutsSoftPi->GetRequireTPCRefit()){
    if(!(b->TestFilterMask(BIT(4)))) return 0;
  }
  
  Int_t returnvalue=1;
  Int_t returnvaluePID=3;


  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    
    Double_t pt=d->Pt();
    Int_t ptbin=PtBin(pt);    
    
    // DStarMass and D0mass
    Double_t mDSPDG = TDatabasePDG::Instance()->GetParticle(413)->Mass();
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    // delta mass PDG
    Double_t deltaPDG = mDSPDG-mD0PDG;
   
    // Half width DStar mass
    if(TMath::Abs(mDSPDG - (d->InvMassDstarKpipi()))>fCutsRD[GetGlobalIndex(9,ptbin)]) return 0;
    // Half width Delta mass
    
    if(TMath::Abs(deltaPDG-(d->DeltaInvMass())) > fCutsRD[GetGlobalIndex(10,ptbin)]) return 0;
    
    // cut on soft pion pt
    if(b->Pt() < fCutsRD[GetGlobalIndex(11,ptbin)] || b->Pt() > fCutsRD[GetGlobalIndex(12,ptbin)]) return 0;
    // cut on the angle between D0 decay plane and soft pion
    if(d->AngleD0dkpPisoft() > fCutsRD[GetGlobalIndex(13,ptbin)]) return 0;
  
    // select D0 that passes D* cuts
    returnvalue = IsD0FromDStarSelected(pt,dd,selectionLevel, aod);
    if((b->Charge()==+1 && returnvalue==2) || (b->Charge()==-1 && returnvalue==1)) return 0; 
    
  }

  fIsSelectedCuts = returnvalue;

  // selection on PID 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate ||
     selectionLevel==AliRDHFCuts::kPID) {
    returnvaluePID = IsSelectedPID(d);
    fIsSelectedPID = returnvaluePID;
  }
  if(returnvaluePID!=3) returnvalue =0;


  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {

    if(fUseTPCtrackCutsOnD0Daughters || !fUseTPCtrackCutsOnSoftPion){
    //if by mistake both flags turned off in cutobject => only apply to D0daughters (default option)
      SetUseTPCtrackCutsOnThisDaughter(kTRUE);
    } else { SetUseTPCtrackCutsOnThisDaughter(kFALSE); }

    if(!AreDaughtersSelected(dd,aod)) return 0;
    if(fTrackCutsSoftPi) {
      AliAODVertex *vAOD = d->GetPrimaryVtx();
      Double_t pos[3],cov[6];
      vAOD->GetXYZ(pos);
      vAOD->GetCovarianceMatrix(cov);
      const AliESDVertex vESD(pos,cov,100.,100);

      if(fUseTPCtrackCutsOnSoftPion){ SetUseTPCtrackCutsOnThisDaughter(kTRUE); }
      else{ SetUseTPCtrackCutsOnThisDaughter(kFALSE); }

      if(!IsDaughterSelected(b,&vESD,fTrackCutsSoftPi,aod)) return 0;
    }
  }
  
  return returnvalue;

}
//_________________________________________________________________________________________________
Int_t AliRDHFCutsDStartoKpipi::IsD0FromDStarSelected(Double_t pt, TObject* obj,Int_t selectionLevel, AliAODEvent* aod) const {
  //
  // Apply selection for D0 from D*. The selection in on D0 prongs
  // added functionality to recalculate the primary vertex without D0 prongs (not default)
  // 
  
  if(!fCutsRD){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  
  AliAODRecoDecayHF2Prong* dd = (AliAODRecoDecayHF2Prong*)obj;
  
  if(!dd){
    cout<<"AliAODRecoDecayHF2Prong null"<<endl;
    return 0;
  }

  // selection on daughter tracks is done in IsSelected()
  
  Int_t returnvalue=1;
  
  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    
    // D0 mass
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    // delta mass PDG
 
    // add vertex recalculation without daughters
    AliAODVertex *origownvtx=0x0;
    if(fRemoveDaughtersFromPrimary  && !fUseMCVertex) {
      if(dd->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*dd->GetOwnPrimaryVtx());
      if(!RecalcOwnPrimaryVtx(dd,aod)) { 
	CleanOwnPrimaryVtx(dd,aod,origownvtx);
	return 0;
      }
    }
    
    
    if(fUseMCVertex) {
      if(dd->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*dd->GetOwnPrimaryVtx());
      if(!SetMCPrimaryVtx(dd,aod)) {
	CleanOwnPrimaryVtx(dd,aod,origownvtx);
	return 0;
      }
    }
    
    Int_t ptbin=PtBin(pt);
    if (ptbin==-1) {
      CleanOwnPrimaryVtx(dd,aod,origownvtx);
      return 0;
    }
    
    Double_t mD0,mD0bar,ctsD0,ctsD0bar;
  
    Int_t okD0     =0;
    Int_t okD0bar  =0;
    okD0=1; okD0bar=1;

    if(dd->PtProng(1) < fCutsRD[GetGlobalIndex(3,ptbin)] || dd->PtProng(0) < fCutsRD[GetGlobalIndex(4,ptbin)]) okD0 = 0;
    if(dd->PtProng(0) < fCutsRD[GetGlobalIndex(3,ptbin)] || dd->PtProng(1) < fCutsRD[GetGlobalIndex(4,ptbin)]) okD0bar = 0;
    
    if(!okD0 && !okD0bar) {CleanOwnPrimaryVtx(dd,aod,origownvtx); return 0;}
    
    if(TMath::Abs(dd->Getd0Prong(1)) > fCutsRD[GetGlobalIndex(5,ptbin)] || 
       TMath::Abs(dd->Getd0Prong(0)) > fCutsRD[GetGlobalIndex(6,ptbin)]) okD0 = 0;
    if(TMath::Abs(dd->Getd0Prong(0)) > fCutsRD[GetGlobalIndex(6,ptbin)] ||
       TMath::Abs(dd->Getd0Prong(1)) > fCutsRD[GetGlobalIndex(5,ptbin)]) okD0bar = 0;
    if(!okD0 && !okD0bar) {CleanOwnPrimaryVtx(dd,aod,origownvtx); return 0;}
    
    if(dd->GetDCA() > fCutsRD[GetGlobalIndex(1,ptbin)]) {CleanOwnPrimaryVtx(dd,aod,origownvtx); return 0;}
    
    dd->InvMassD0(mD0,mD0bar);
    if(TMath::Abs(mD0-mD0PDG) > fCutsRD[GetGlobalIndex(0,ptbin)]) okD0 = 0;
    if(TMath::Abs(mD0bar-mD0PDG) > fCutsRD[GetGlobalIndex(0,ptbin)]) okD0bar = 0;
    if(!okD0 && !okD0bar) {CleanOwnPrimaryVtx(dd,aod,origownvtx); return 0;}
    
    dd->CosThetaStarD0(ctsD0,ctsD0bar);
    if(TMath::Abs(ctsD0) > fCutsRD[GetGlobalIndex(2,ptbin)]) okD0 = 0;
    if(TMath::Abs(ctsD0bar) > fCutsRD[GetGlobalIndex(2,ptbin)]) okD0bar = 0;
    if(!okD0 && !okD0bar) {CleanOwnPrimaryVtx(dd,aod,origownvtx); return 0;}
    
    if(dd->Prodd0d0() > fCutsRD[GetGlobalIndex(7,ptbin)]) {CleanOwnPrimaryVtx(dd,aod,origownvtx); return 0;}
    
    if(dd->CosPointingAngle() < fCutsRD[GetGlobalIndex(8,ptbin)]) {CleanOwnPrimaryVtx(dd,aod,origownvtx); return 0;}

    if(TMath::Abs(dd->CosPointingAngleXY()) < fCutsRD[GetGlobalIndex(14,ptbin)]) {CleanOwnPrimaryVtx(dd,aod,origownvtx); return 0;}
	
    Double_t normalDecayLengXY=(dd->NormalizedDecayLengthXY()*(dd->P()/dd->Pt()));
    if (normalDecayLengXY < fCutsRD[GetGlobalIndex(15, ptbin)]) {CleanOwnPrimaryVtx(dd,aod,origownvtx); return 0;}

    if (okD0) returnvalue=1; //cuts passed as D0
    if (okD0bar) returnvalue=2; //cuts passed as D0bar
    if (okD0 && okD0bar) returnvalue=3; //both

    // unset recalculated primary vertex when not needed any more
    CleanOwnPrimaryVtx(dd,aod,origownvtx);

  }
 
  return returnvalue;
}
//----------------------------------------------------------------------------------
Bool_t AliRDHFCutsDStartoKpipi::IsInFiducialAcceptance(Double_t pt, Double_t y) const
{
  //
  // D* fiducial acceptance region 
  //

  if(fMaxRapidityCand>-998.){
    if(TMath::Abs(y) > fMaxRapidityCand) return kFALSE;
    else return kTRUE;
  }

  if(pt > 5.) {
    // applying cut for pt > 5 GeV
    AliDebug(4,Form("pt of D* = %f (> 5), cutting at |y| < 0.8\n",pt)); 
    if (TMath::Abs(y) > 0.8){
      return kFALSE;
    }
  } else {    
    // appliying smooth cut for pt < 5 GeV
    Double_t maxFiducialY = -0.2/15*pt*pt+1.9/15*pt+0.5; 
    Double_t minFiducialY = 0.2/15*pt*pt-1.9/15*pt-0.5;		
    AliDebug(2,Form("pt of D* = %f (< 5), cutting  according to the fiducial zone [%f, %f]\n",pt,minFiducialY,maxFiducialY)); 
    if (y < minFiducialY || y > maxFiducialY){
      return kFALSE;
    }
  }
    
  return kTRUE;
}
//_______________________________________________________________________________-
Int_t AliRDHFCutsDStartoKpipi::IsSelectedPID(AliAODRecoDecayHF* obj)
{
  Int_t retCode=3;
  if(!fUsePID) return retCode;
  
  AliAODRecoCascadeHF* dstar = (AliAODRecoCascadeHF*)obj;
  if(!dstar){
    cout<<"AliAODRecoCascadeHF null"<<endl;
    return 0;
  }
  AliAODRecoDecayHF2Prong* d0 = (AliAODRecoDecayHF2Prong*)dstar->Get2Prong();
  if(!d0){
    cout<<"AliAODRecoDecayHF2Prong null"<<endl;
    return 0;
  }

  if(dstar) {
    Double_t pt = dstar->Pt();
    TObjArray aodTracks(3);
    aodTracks.AddAt(dstar->GetBachelor(),0); //soft pion
    aodTracks.AddAt(d0->GetDaughter(0),1); //D0 prong0
    aodTracks.AddAt(d0->GetDaughter(1),2); //D0 prong1

    retCode = IsSelectedPID(pt,aodTracks);
  }
    
  return retCode;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsDStartoKpipi::IsSelectedPID(Double_t pt, TObjArray aodTracks)
{
  //
  // PID method, n signa approach default
  //
   
  if(!fUsePID || pt > fMaxPtPid) return 3;

  //  here the PID
  AliAODTrack *softPion = (AliAODTrack*)aodTracks.At(0);
  AliAODTrack *pos = (AliAODTrack*)aodTracks.At(1);
  AliAODTrack *neg = (AliAODTrack*)aodTracks.At(2);

  if (softPion->Charge()>0){
    if(!SelectPID(pos,2)) return 0;//pion+
    if(!SelectPID(neg,3)) return 0;//kaon-
  }else{
    if(!SelectPID(pos,3)) return 0;//kaon+
    if(!SelectPID(neg,2)) return 0;//pion-
  }

  if ((fPidHF->GetMatch() == 10 || fPidHF->GetMatch() == 11) && fPidHF->GetITS()) { //ITS n sigma band
    if (fPidHF->CheckBands(AliPID::kPion, AliPIDResponse::kITS, softPion) == -1) {
      return 0;
    }
  }

  return 3;
}

//_______________________________________________________________________________-
Int_t AliRDHFCutsDStartoKpipi::SelectPID(AliAODTrack *track, Int_t type)
{
  //
  //  here the PID
    
  Bool_t isParticle=kTRUE;
  Int_t match = fPidHF->GetMatch();

  if(match == 1){//n-sigma
    Bool_t TPCon=TMath::Abs(2)>1e-4?kTRUE:kFALSE;
    Bool_t TOFon=TMath::Abs(3)>1e-4?kTRUE:kFALSE;
    
    Bool_t isTPC=kTRUE;
    Bool_t isTOF=kTRUE;

    if (TPCon){//TPC
      if(fPidHF->CheckStatus(track,"TPC")){
	if(type==2) isTPC=fPidHF->IsPionRaw(track,"TPC");
	if(type==3) isTPC=fPidHF->IsKaonRaw(track,"TPC");
      }
    }
    if (TOFon){//TOF
      if(fPidHF->CheckStatus(track,"TOF")){
	if(type==2) isTOF=fPidHF->IsPionRaw(track,"TOF");
	if(type==3) isTOF=fPidHF->IsKaonRaw(track,"TOF");
      }
    }

    //--------------------------------
    // cut on high momentum in the TPC
    //--------------------------------
    Double_t pPIDcut = track->P();
    if(pPIDcut>fTPCflag) isTPC=1;
    
    isParticle = isTPC&&isTOF;
  }
  
  if(match == 2){//bayesian
    //Double_t priors[5]={0.01,0.001,0.3,0.3,0.3};
    Double_t prob[5]={1.,1.,1.,1.,1.};
    
    //fPidHF->SetPriors(priors,5);
    //    fPidHF->BayesianProbability(track,prob);
    
    Double_t max=0.;
    Int_t k=-1;
    for (Int_t i=0; i<5; i++) {
      if (prob[i]>max) {k=i; max=prob[i];}
    }
    isParticle = Bool_t(k==type);
  }

  if (match == 10 || match == 11) { //Assymetric PID using histograms
    Int_t checkTPC = fPidHF->CheckBands((AliPID::EParticleType) type, AliPIDResponse::kTPC, track);
    Int_t checkTOF = fPidHF->CheckBands((AliPID::EParticleType) type, AliPIDResponse::kTOF, track);

    isParticle = checkTPC >= 0 && checkTOF >= 0 ? kTRUE : kFALSE; //Standard: both at least compatible
    if (match == 11) { //Extra requirement: at least one identified
      isParticle = isParticle && checkTPC+checkTOF >= 1 ? kTRUE : kFALSE;
    }
  }
  
  if (match == 12) { //Circular cut
    Double_t nSigmaTOF = 0;
    Double_t nSigmaTPC = 0;

    Double_t radius = fCircRadius;

    isParticle = kTRUE;
    if (radius > 0) {
      Int_t TPCok = fPidHF->GetnSigmaTPC(track, type, nSigmaTPC);
      Int_t TOFok = fPidHF->GetnSigmaTOF(track, type, nSigmaTOF);
      if (TPCok != -1 && TOFok != -1) {
        //Both detectors gave PID information
        isParticle = TMath::Sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF) <= radius ? kTRUE : kFALSE;
      }
      else {
        //Maximum one detector gave PID information
        if (TPCok != -1) {
          isParticle = nSigmaTPC <= radius ? kTRUE : kFALSE;
        }
        if (TOFok != -1) {
          isParticle = nSigmaTOF <= radius ? kTRUE : kFALSE;
        }
      }
    }
  }

  return isParticle;
  
}
//__________________________________________________________________________________-
void  AliRDHFCutsDStartoKpipi::SetStandardCutsPP2010() {
  //
  //STANDARD CUTS USED FOR 2010 pp analysis 
  //                                           
  // Need to be updated for the final cut version
  //

  SetName("DStartoD0piCutsStandard");
  SetTitle("Standard Cuts for D* analysis");
  
  // PILE UP REJECTION
  SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  // EVENT CUTS
  SetMinVtxContr(1);
  
  // CUTS ON SINGLE TRACKS
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.3,1.e10);
  
  // CUTS on SOFT PION
  AliESDtrackCuts* esdSoftPicuts=new AliESDtrackCuts();
  esdSoftPicuts->SetRequireSigmaToVertex(kFALSE);
  esdSoftPicuts->SetRequireTPCRefit(kFALSE);
  esdSoftPicuts->SetRequireITSRefit(kFALSE);
  esdSoftPicuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					  AliESDtrackCuts::kAny); 
  esdSoftPicuts->SetPtRange(0.0,1.e10);

  AddTrackCuts(esdTrackCuts);
  AddTrackCutsSoftPi(esdSoftPicuts);
  delete esdTrackCuts;
  esdTrackCuts=NULL;
  delete esdSoftPicuts;
  esdSoftPicuts=NULL;

  const Int_t nptbins =13;
  const Double_t ptmax = 9999.;
  const Int_t nvars=16;
  Float_t ptbins[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=0.5;	
  ptbins[2]=1.;
  ptbins[3]=2.;
  ptbins[4]=3.;
  ptbins[5]=4.;
  ptbins[6]=5.;
  ptbins[7]=6.;
  ptbins[8]=7.;
  ptbins[9]=8.;
  ptbins[10]=12.;
  ptbins[11]=16.;
  ptbins[12]=24.;
  ptbins[13]=ptmax;

  SetGlobalIndex(nvars,nptbins);
  SetPtBins(nptbins+1,ptbins);
  
  Float_t cutsMatrixD0toKpiStand[nptbins][nvars]={{0.7,220.*1E-4,0.7,0.21,0.21,500.*1E-4,500.*1E-4,-2000.*1E-8,0.85,0.3,0.1,0.05,100,0.5,-1.,0.},/* pt<0.5*/
						  {0.7,220.*1E-4,0.7,0.21,0.21,500.*1E-4,500.*1E-4,-16000.*1E-8,0.85,0.3,0.1,0.05,100,0.5,-1.,0.},/* 0.5<pt<1*/
						  {0.7,400.*1E-4,0.8,0.7,0.7,400.*1E-4,400.*1E-4,-36000.*1E-8,0.82,0.3,0.1,0.05,100,0.5,-1.,0.},/* 1<pt<2 */
						  {0.7,200.*1E-4,0.8,0.7,0.7,800.*1E-4,800.*1E-4,-16000.*1E-8,0.9,0.3,0.1,0.05,100,0.5,-1.,0.},/* 2<pt<3 */
						  {0.7,500.*1E-4,0.8,1.0,1.0,420.*1E-4,560.*1E-4,-6500.*1E-8,0.9,0.3,0.1,0.05,100,0.5,-1.,0.},/* 3<pt<4 */
						  {0.7,800.*1E-4,0.9,1.2,1.2,700.*1E-4,700.*1E-4,1000.*1E-8,0.9,0.3,0.1,0.05,100,0.5,-1.,0.},/* 4<pt<5 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,800.*1E-4,800.*1E-4,50000.*1E-8,0.8,0.3,0.1,0.05,100,0.5,-1.,0.},/* 5<pt<6 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,1000.*1E-4,1000.*1E-4,100000.*1E-8,0.7,0.3,0.1,0.05,100,0.5,-1.,0.},/* 6<pt<7 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,1000.*1E-4,1000.*1E-4,100000.*1E-8,0.7,0.3,0.1,0.05,100,0.5,-1.,0.},/* 7<pt<8 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,1000.*1E-4,1000.*1E-4,600000.*1E-8,0.7,0.3,0.1,0.05,100,0.5,-1.,0.},/* 8<pt<12 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,1500.*1E-4,1500.*1E-4,1000000.*1E-8,0.7,0.3,0.1,0.05,100,0.5,-1.,0.},/* 12<pt<16 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,1500.*1E-4,1500.*1E-4,1000000.*1E-8,0.7,0.3,0.1,0.05,100,0.5,-1.,0.},/* 16<pt<20 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,1500.*1E-4,1500.*1E-4,1000000.*1E-8,0.7,0.3,0.1,0.05,100,0.5,-1.,0.}};/* pt>24 */
  
  
  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
  
  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar = 0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];      
    }
  }
  
  SetCuts(nvars,nptbins,cutsMatrixTransposeStand);

  for(Int_t iv=0;iv<nvars;iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand=NULL;

  // PID SETTINGS FOR D* analysis
  AliAODPidHF* pidObj=new AliAODPidHF();
  //pidObj->SetName("pid4DSatr");
  Int_t mode=1;
  Double_t priors[5]={0.01,0.001,0.3,0.3,0.3};
  pidObj->SetPriors(priors,5);
  pidObj->SetMatch(mode);
  pidObj->SetSigma(0,2); // TPC
  pidObj->SetSigma(3,3); // TOF
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  
  SetPidHF(pidObj);
  SetUsePID(kTRUE);

  //  PrintAll();

  delete pidObj;
  pidObj=NULL;

  return;
}
//_____________________________________________________________________________-
void  AliRDHFCutsDStartoKpipi::SetStandardCutsPbPb2010(){  
  //
  // TEMPORARY, WORK IN PROGRESS ... BUT WORKING! 
  //
  //  Lead Lead
  //

  SetName("DStartoD0piCutsStandard");
  SetTitle("Standard Cuts for D* analysis in PbPb 2010");

  // EVENT CUTS
  SetMinVtxContr(1);
  
  // CUTS ON SINGLE TRACKS
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.3,1.e10);
  
  // CUTS on SOFT PION
  AliESDtrackCuts* esdSoftPicuts=new AliESDtrackCuts();
  esdSoftPicuts->SetRequireSigmaToVertex(kFALSE);
  esdSoftPicuts->SetRequireTPCRefit(kTRUE);
  esdSoftPicuts->SetRequireITSRefit(kTRUE);
  esdSoftPicuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					  AliESDtrackCuts::kAny); //test d0 asimmetry
  esdSoftPicuts->SetPtRange(0.25,5);

  AddTrackCuts(esdTrackCuts);
  AddTrackCutsSoftPi(esdSoftPicuts);
  delete esdTrackCuts;
  esdTrackCuts=NULL;
  delete esdSoftPicuts;
  esdSoftPicuts=NULL;

  const Int_t nptbins =13;
  const Double_t ptmax = 9999.;
  const Int_t nvars=16;
  Float_t ptbins[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=0.5;	
  ptbins[2]=1.;
  ptbins[3]=2.;
  ptbins[4]=3.;
  ptbins[5]=4.;
  ptbins[6]=5.;
  ptbins[7]=6.;
  ptbins[8]=7.;
  ptbins[9]=8.;
  ptbins[10]=12.;
  ptbins[11]=16.;
  ptbins[12]=24.;
  ptbins[13]=ptmax;

  SetGlobalIndex(nvars,nptbins);
  SetPtBins(nptbins+1,ptbins);
  
  Float_t cutsMatrixD0toKpiStand[nptbins][nvars]={{0.7,220.*1E-4,0.7,0.21,0.21,500.*1E-4,500.*1E-4,-2000.*1E-8,0.85,0.3,0.1,0.05,100,0.5,-1.,0.},/* pt<0.5*/
						  {0.7,220.*1E-4,0.7,0.21,0.21,500.*1E-4,500.*1E-4,-16000.*1E-8,0.85,0.3,0.1,0.05,100,0.5,-1.,0.},/* 0.5<pt<1*/
						  {0.7,400.*1E-4,0.8,0.7,0.7,800.*1E-4,800.*1E-4,-36000.*1E-8,0.82,0.3,0.1,0.05,100,0.5,-1.,0.},/* 1<pt<2 */
						  {0.7,200.*1E-4,0.8,0.7,0.7,800.*1E-4,800.*1E-4,-16000.*1E-8,0.9,0.3,0.1,0.05,100,0.5,-1.,0.},/* 2<pt<3 */
						  {0.7,500.*1E-4,0.8,1.0,1.0,420.*1E-4,560.*1E-4,-6500.*1E-8,0.9,0.3,0.1,0.05,100,0.5,-1.,0.},/* 3<pt<4 */
						  {0.7,800.*1E-4,0.9,1.2,1.2,700.*1E-4,700.*1E-4,1000.*1E-8,0.9,0.3,0.1,0.05,100,0.5,-1.,0.},/* 4<pt<5 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,800.*1E-4,800.*1E-4,50000.*1E-8,0.8,0.3,0.1,0.05,100,0.5,-1.,0.},/* 5<pt<6 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,1000.*1E-4,1000.*1E-4,100000.*1E-8,0.7,0.3,0.1,0.05,100,0.5,-1.,0.},/* 6<pt<7 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,1000.*1E-4,1000.*1E-4,100000.*1E-8,0.7,0.3,0.1,0.05,100,0.5,-1.,0.},/* 7<pt<8 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,1000.*1E-4,1000.*1E-4,600000.*1E-8,0.7,0.3,0.1,0.05,100,0.5,-1.,0.},/* 8<pt<12 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,1500.*1E-4,1500.*1E-4,1000000.*1E-8,0.7,0.3,0.1,0.05,100,0.5,-1.,0.},/* 12<pt<16 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,1500.*1E-4,1500.*1E-4,1000000.*1E-8,0.7,0.3,0.1,0.05,100,0.5,-1.,0.},/* 16<pt<24 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,1500.*1E-4,1500.*1E-4,1000000.*1E-8,0.7,0.3,0.1,0.05,100,0.5,-1.,0.}};/* pt>24 */
  
  
  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
  
  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar = 0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];      
    }
  }
  
  SetCuts(nvars,nptbins,cutsMatrixTransposeStand);

  for(Int_t iv=0;iv<nvars;iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand=NULL;
  
  // PID SETTINGS
  AliAODPidHF* pidObj=new AliAODPidHF();
  // pidObj->SetName("pid4DSatr");
  Int_t mode=1;
  Double_t priors[5]={0.01,0.001,0.3,0.3,0.3};
  pidObj->SetPriors(priors,5);
  pidObj->SetMatch(mode);
  pidObj->SetSigma(0,2); // TPC
  pidObj->SetSigma(3,3); // TOF
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  
  SetPidHF(pidObj);
  SetUsePID(kTRUE);

  //  PrintAll();

  delete pidObj;
  pidObj=NULL;

  return;

}

//_____________________________________________________________________________
void  AliRDHFCutsDStartoKpipi::SetStandardCutsPbPb2011(){
  //
  // Not implemented !!
  //
  SetStandardCutsPbPb2011DStar(0);
  return;
}
//_________________________________here the PbPb vs pt _________________________________________
void  AliRDHFCutsDStartoKpipi::SetStandardCutsPbPb2011DStar(TH1F *hfl){

  // Default 2010 PbPb cut object

  SetName("DStartoD0piCutsStandard2011");
  SetTitle("Standard Cuts for D* analysis in PbPb 2011");

  // EVENT CUTS
  SetMinVtxContr(1);
  // MAX Z-VERTEX CUT
  SetMaxVtxZ(10.);

  SetTriggerClass("");
  ResetMaskAndEnableMBTrigger();
  EnableCentralTrigger();
  EnableSemiCentralTrigger();

  // CENTRALITY SELECTION
  SetMinCentrality(0.);
  SetMaxCentrality(10.);
  SetUseCentrality(AliRDHFCuts::kCentV0M);

  // CUTS ON SINGLE TRACKS
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.3,1.e10);
  
  esdTrackCuts->SetMaxDCAToVertexXY(1.);  
  esdTrackCuts->SetMaxDCAToVertexZ(1.);


  // CUTS on SOFT PION
  AliESDtrackCuts* esdSoftPicuts=new AliESDtrackCuts();
  esdSoftPicuts->SetRequireSigmaToVertex(kFALSE);
  esdSoftPicuts->SetRequireTPCRefit(kTRUE);
  esdSoftPicuts->SetRequireITSRefit(kTRUE);
  esdSoftPicuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					  AliESDtrackCuts::kAny); //test d0 asimmetry
  esdSoftPicuts->SetPtRange(0.1,10);

  esdSoftPicuts->SetMaxDCAToVertexXY(1.);  
  esdSoftPicuts->SetMaxDCAToVertexZ(1.);

  SetSelectCandTrackSPDFirst(kTRUE, 4);

  //nothing below 3 GeV/c, speed up the calculation
  SetMinPtCandidate(3.0);

  AddTrackCuts(esdTrackCuts);
  AddTrackCutsSoftPi(esdSoftPicuts);
  delete esdTrackCuts;
  esdTrackCuts=NULL;
  delete esdSoftPicuts;
  esdSoftPicuts=NULL;

  const Int_t nptbins =14;
  const Double_t ptmax = 36.;
  const Int_t nvars=16;
  Float_t ptbins[nptbins+1];
  ptbins[0]=0.5;	
  ptbins[1]=1.;
  ptbins[2]=2.;
  ptbins[3]=3.;
  ptbins[4]=4.;
  ptbins[5]=5.;
  ptbins[6]=6.;
  ptbins[7]=7.;
  ptbins[8]=8.;
  ptbins[9]=10.;
  ptbins[10]=12.;
  ptbins[11]=16.;
  ptbins[12]=20.;
  ptbins[13]=24.;
  ptbins[14]=ptmax;

  SetGlobalIndex(nvars,nptbins);
  SetPtBins(nptbins+1,ptbins);
  
  Float_t cutsMatrixD0toKpiStand[nptbins][nvars]={{0.032,220.*1E-4,0.9,0.5,0.5,500.*1E-4,500.*1E-4,-16000.*1E-8,0.85,0.3,0.1,0.05,100,0.5,-1.,9.},/* 0.5<pt<1*/
						  {0.032,350.*1E-4,0.9,0.5,0.5,800.*1E-4,800.*1E-4,-20000.*1E-8,0.96,0.3,0.15,0.05,100,0.5,0.99,10.},/* 1<pt<2 */
						  {0.032,300.*1E-4,0.9,0.5,0.5,900.*1E-4,900.*1E-4,-42000.*1E-8,0.96,0.3,0.15,0.05,100,0.5,0.99,9.},/* 2<pt<3 */
						  {0.036,300.*1E-4,0.8,0.8,0.8,900.*1E-4,900.*1E-4,-39000.*1E-8,0.99,0.3,0.15,0.05,100,0.5,0.998,8.},/* 3<pt<4 */
						  {0.038,225.*1E-4,0.8,1.0,1.0,1000.*1E-4,1000.*1E-4,-30000.*1E-8,0.99,0.3,0.15,0.05,100,0.5,0.998,7.5},/* 4<pt<5 */
						  {0.045,200.*1E-4,1.0,1.0,1.0,1000.*1E-4,1000.*1E-4,-23000.*1E-8,0.99,0.3,0.15,0.05,100,0.5,0.998,7.},/* 5<pt<6 */
						  {0.045,210.*1E-4,1.0,1.0,1.0,1000.*1E-4,1000.*1E-4,-10000.*1E-8,0.982,0.3,0.15,0.05,100,0.5,0.998,6.4},/* 6<pt<7 */
						  {0.050,230.*1E-4,1.0,1.0,1.0,1000.*1E-4,1000.*1E-4,-12700.*1E-8,0.98,0.3,0.15,0.05,100,0.5,0.998,6.4},/* 7<pt<8 */
						  {0.060,200.*1E-4,1.0,0.9,0.9,1000.*1E-4,1000.*1E-4,-7500.*1E-8,0.98,0.3,0.15,0.05,100,0.5,0.998,4.7},/* 8<pt<10 */
						  {0.060,200.*1E-4,1.0,0.9,0.9,1000.*1E-4,1000.*1E-4,-7500.*1E-8,0.97,0.3,0.15,0.05,100,1.0,0.998,4.7},/* 10<pt<12 */
						  {0.074,200.*1E-4,1.0,0.5,0.5,1500.*1E-4,1500.*1E-4,-7500.*1E-8,0.95,0.3,0.15,0.05,100,1.0,0.998,3},/* 12<pt<16 */
						  {0.074,210.*1E-4,1.0,0.5,0.5,1500.*1E-4,1500.*1E-4,-5000.*1E-8,0.95,0.3,0.15,0.05,100,1.0,0.998,2.},/* 16<pt<20 */
                                                  {0.074,220.*1E-4,1.0,0.5,0.5,1500.*1E-4,1500.*1E-4,-5000.*1E-8,0.93,0.3,0.15,0.05,100,1.0,0.995,2.},/* 20<pt<24 */
						  {0.074,400.*1E-4,1.0,0.5,0.5,2000.*1E-4,2000.*1E-4,40000.*1E-8,0.7,0.3,0.15,0.05,100,1.0,0.9,1.}};/* 24<pt<36 */
  
  
  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
  
  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar = 0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];      
    }
  }
  
  SetCuts(nvars,nptbins,cutsMatrixTransposeStand);

  for(Int_t iv=0;iv<nvars;iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand=NULL;
  
  // PID SETTINGS // --- 
  AliAODPidHF* pidObj=new AliAODPidHF();
  // pidObj->SetName("pid4DSatr");
  Int_t mode=1;

  pidObj->SetMatch(mode);
  pidObj->SetSigma(0,3); // TPC -- 2 sigma for pt < 4
  pidObj->SetSigma(3,3); // TOF
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  pidObj->SetOldPid(kFALSE);
  
  SetPidHF(pidObj);
  SetUsePID(kTRUE);

  // PID off for tracks with pt above 4 GeV/c
  SetOffHighPtPIDinTPC(4.0);
  SetRemoveDaughtersFromPrim(kFALSE);
  // flattening
  SetHistoForCentralityFlattening(hfl,0.,10,0.,0);

  //  PrintAll();

  delete pidObj;
  pidObj=NULL;


  return;
}
//-----------------------------------Here the multiplicity pp--------------------------------

void  AliRDHFCutsDStartoKpipi::SetStandardCutsPP2010DStarMult(Bool_t rec){


 //
  // STANDARD CUTS USED FOR 2010 pp analysis (multiplicity)
  //                                           
  //

  SetName("DStartoD0piCutsStandard");
  SetTitle("Standard Cuts for D* analysis pp mult");
  
  // PILE UP REJECTION
  SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  // EVENT CUTS
  SetMinVtxContr(1);
   // MAX Z-VERTEX CUT
  SetMaxVtxZ(10.);

  // CUTS ON SINGLE TRACKS
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.3,1.e10);

  esdTrackCuts->SetMaxDCAToVertexXY(1.);  
  esdTrackCuts->SetMaxDCAToVertexZ(1.);

  
  // CUTS on SOFT PION
  AliESDtrackCuts* esdSoftPicuts=new AliESDtrackCuts();
  esdSoftPicuts->SetRequireSigmaToVertex(kFALSE);
  esdSoftPicuts->SetRequireTPCRefit(kFALSE);
  esdSoftPicuts->SetRequireITSRefit(kTRUE);
  esdSoftPicuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					  AliESDtrackCuts::kAny); 
  esdSoftPicuts->SetPtRange(0.08,1.e10);
  SetUseCentrality(kFALSE);

  AddTrackCuts(esdTrackCuts);
  AddTrackCutsSoftPi(esdSoftPicuts);
  delete esdTrackCuts;
  esdTrackCuts=NULL;
  delete esdSoftPicuts;
  esdSoftPicuts=NULL;

  const Int_t nptbins =14;
  const Double_t ptmax = 9999.;
  const Int_t nvars=16;
  Float_t ptbins[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=0.5;	
  ptbins[2]=1.;
  ptbins[3]=2.;
  ptbins[4]=3.;
  ptbins[5]=4.;
  ptbins[6]=5.;
  ptbins[7]=6.;
  ptbins[8]=7.;
  ptbins[9]=8.;
  ptbins[10]=10.;
  ptbins[11]=12.;
  ptbins[12]=16.;
  ptbins[13]=20.;
  ptbins[14]=ptmax;

  SetGlobalIndex(nvars,nptbins);
  SetPtBins(nptbins+1,ptbins);
  
  Float_t cutsMatrixD0toKpiStand[nptbins][nvars]={{0.026,220.*1E-4,0.7,0.21,0.21,500.*1E-4,500.*1E-4,-2000.*1E-8,0.85,0.3,0.1,0.05,100,0.5,-1.,0.},/* pt<0.5*/
						  {0.039,300.*1E-4,0.7,0.21,0.21,500.*1E-4,500.*1E-4,-16000.*1E-8,0.85,0.3,0.1,0.05,100,0.5,-1.,0.},/* 0.5<pt<1*/
						  {0.022,400.*1E-4,0.8,0.7,0.7,800.*1E-4,800.*1E-4,-20000.*1E-8,0.8,0.3,0.3,0.05,100,0.5,-1.,0.},/* 1<pt<2 */
						  {0.032,350.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-13000.*1E-8,0.9,0.3,0.3,0.05,100,0.8,-1.,0.},/* 2<pt<3 */
						  {0.032,500.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.9,0.3,0.3,0.05,100,0.8,-1.,0.},/* 3<pt<4 */
						  {0.032,700.*1E-4,0.9,0.9,0.9,1000.*1E-4,1000.*1E-4,10000.*1E-8,0.9,0.3,0.3,0.05,100,1.0,-1.,0.},/* 4<pt<5 */
						  {0.036,1000.*1E-4,1.0,0.8,0.8,1000.*1E-4,1000.*1E-4,50000.*1E-8,0.9,0.3,0.3,0.05,100,1.0,-1.,0.},/* 5<pt<6 */
						  {0.036,1000.*1E-4,1.0,0.8,0.8,1000.*1E-4,1000.*1E-4,100000.*1E-8,0.7,0.3,0.3,0.05,100,1.0,-1.,0.},/* 6<pt<7 */
						  {0.036,1000.*1E-4,1.0,0.8,0.8,1200.*1E-4,1200.*1E-4,100000.*1E-8,0.6,0.3,0.3,0.05,100,1.0,-1.,0.},/* 7<pt<8 */
						  {0.065,2000.*1E-4,1.0,0.3,0.3,2000.*1E-4,2000.*1E-4,1000000.*1E-8,0.5,0.3,0.3,0.05,100,1.0,-1.,0.},/* 8<pt<10 */
						  {0.075,2000.*1E-4,1.0,0.3,0.3,2000.*1E-4,2000.*1E-4,1000000.*1E-8,0.3,0.3,0.3,0.05,100,1.0,-1.,0.},/* 10<pt<12 */
						  {0.084,6000.*1E-4,1.0,0.3,0.3,1500.*1E-4,1500.*1E-4,1000000.*1E-8,0.1,0.3,0.1,0.05,100,1.0,-1.,0.},/* 12<pt<16 */
						  {0.084,1000.*1E-4,1.0,1.0,1.0,1500.*1E-4,1500.*1E-4,1000000.*1E-8,0.7,0.3,0.1,0.05,100,1.0,-1.,0.},/* 16<pt<20 */
						  {0.7,1000.*1E-4,1.0,1.0,1.0,1500.*1E-4,1500.*1E-4,1000000.*1E-8,0.7,0.3,0.1,0.05,100,1.0,-1.,0.}};/* pt>24 */
  
  
  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
  
  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar = 0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];      
    }
  }
  
  SetCuts(nvars,nptbins,cutsMatrixTransposeStand);

  for(Int_t iv=0;iv<nvars;iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand=NULL;

  // remove daughters from primary vertex
  SetRemoveDaughtersFromPrim(rec);

  // PID SETTINGS FOR D* analysis
  AliAODPidHF* pidObj=new AliAODPidHF();
  //pidObj->SetName("pid4DSatr");
  Int_t mode=1;
  pidObj->SetMatch(mode);
  pidObj->SetSigma(0,3); // TPC
  pidObj->SetSigma(3,3); // TOF
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  
  SetPidHF(pidObj);
  SetUsePID(kTRUE);
  pidObj->SetOldPid(kTRUE);

  //  PrintAll();

  delete pidObj;
  pidObj=NULL;

  return;


}
