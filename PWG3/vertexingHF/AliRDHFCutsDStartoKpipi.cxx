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
// Class for cuts on AOD reconstructed DStar->Kpipi
//
// Author: A.Grelli, alessandro.grelli@uu.nl
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>

#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"

ClassImp(AliRDHFCutsDStartoKpipi)

//--------------------------------------------------------------------------
AliRDHFCutsDStartoKpipi::AliRDHFCutsDStartoKpipi(const char* name) : 
AliRDHFCuts(name),
fTrackCutsSoftPi(0)
{
  //
  // Default Constructor
  //
 
  Int_t nvars=14;
  SetNVars(nvars);
  TString varNames[14]={
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
    "theta, angle between the pi_s and decay plane of the D0 [rad]"};
  Bool_t isUpperCut[14]={
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
    kFALSE};
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[14]={
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
    kFALSE};
  SetVarsForOpt(5,forOpt);
  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);
}
//--------------------------------------------------------------------------
AliRDHFCutsDStartoKpipi::AliRDHFCutsDStartoKpipi(const AliRDHFCutsDStartoKpipi &source) :
  AliRDHFCuts(source),
  fTrackCutsSoftPi(0)
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

  if(source.GetTrackCutsSoftPi()) AddTrackCutsSoftPi(source.GetTrackCutsSoftPi());
  return *this;
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
 
  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsDStartoKpipi::IsSelected(TObject* obj,Int_t selectionLevel) {
  //
  // Apply selection for D*.
  //
  if(!fCutsRD){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  
  AliAODRecoCascadeHF* d = (AliAODRecoCascadeHF*)obj;
  if(!d){
    cout<<"AliAODRecoCascadeHF null"<<endl;
    return 0;
  }
  
  AliAODRecoDecayHF2Prong* dd = (AliAODRecoDecayHF2Prong*)d->Get2Prong();  
  if(!dd){
    cout<<"AliAODRecoDecayHF2Prong null"<<endl;
    return 0;
  }

  AliAODTrack *b = (AliAODTrack*)d->GetBachelor();

  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(dd)) return 0;
    if(fTrackCutsSoftPi) {
      AliAODVertex *vAOD = d->GetPrimaryVtx();
      Double_t pos[3],cov[6];
      vAOD->GetXYZ(pos);
      vAOD->GetCovarianceMatrix(cov);
      const AliESDVertex vESD(pos,cov,100.,100);
      if(!IsDaughterSelected(b,&vESD,fTrackCutsSoftPi)) return 0;
    }
  }
  
  Int_t returnvalue=1;
  
  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    
    Double_t pt=d->Pt();
    Int_t ptbin=PtBin(pt);

    // select D0 that passes D* cuts
    returnvalue = IsD0FromDStarSelected(pt,dd,selectionLevel);

    // DStarMass and D0mass
    Double_t mDSPDG = TDatabasePDG::Instance()->GetParticle(413)->Mass();
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    // delta mass PDG
    Double_t deltaPDG = mDSPDG-mD0PDG;
    
    Double_t mD0,mD0bar;
  
    Int_t okDStarP =0;
    Int_t okDStarM =0;
   
    okDStarP=1; okDStarM=1; 

    dd->InvMassD0(mD0,mD0bar);
    // Half width DStar mass
    if(TMath::Abs((d->InvMassDstarKpipi()-mDSPDG))>fCutsRD[GetGlobalIndex(9,ptbin)]) return 0;
    // Half width Delta mass
    
    if(TMath::Abs((d->InvMassDstarKpipi()-mD0)-deltaPDG) > fCutsRD[GetGlobalIndex(10,ptbin)]) okDStarP =0;
    if(TMath::Abs((d->InvMassDstarKpipi()-mD0bar)-deltaPDG) > fCutsRD[GetGlobalIndex(10,ptbin)]) okDStarM =0;
    if(!okDStarP && !okDStarM) return 0;
    
    // cut on soft pion pt
    if(b->Pt() < fCutsRD[GetGlobalIndex(11,ptbin)] || b->Pt() > fCutsRD[GetGlobalIndex(12,ptbin)]) return 0;
    // cut on the angle between D0 decay plane and soft pion
    if(d->AngleD0dkpPisoft() > fCutsRD[GetGlobalIndex(13,ptbin)]) return 0;
  
  }

  return returnvalue;
}
//_________________________________________________________________________________________________
Int_t AliRDHFCutsDStartoKpipi::IsD0FromDStarSelected(Double_t pt, TObject* obj,Int_t selectionLevel) const {
  //
  // Apply selection for D0 from D*. The selection in on D0 prongs
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
 
    Int_t ptbin=PtBin(pt);
    
    Double_t mD0,mD0bar,ctsD0,ctsD0bar;
  
    Int_t okD0     =0;
    Int_t okD0bar  =0;
    okD0=1; okD0bar=1;

    if(dd->PtProng(1) < fCutsRD[GetGlobalIndex(3,ptbin)] || dd->PtProng(0) < fCutsRD[GetGlobalIndex(4,ptbin)]) okD0 = 0;
    if(dd->PtProng(0) < fCutsRD[GetGlobalIndex(3,ptbin)] || dd->PtProng(1) < fCutsRD[GetGlobalIndex(4,ptbin)]) okD0bar = 0;
    
    if(!okD0 && !okD0bar) return 0;
    
    if(TMath::Abs(dd->Getd0Prong(1)) > fCutsRD[GetGlobalIndex(5,ptbin)] || 
       TMath::Abs(dd->Getd0Prong(0)) > fCutsRD[GetGlobalIndex(6,ptbin)]) okD0 = 0;
    if(TMath::Abs(dd->Getd0Prong(0)) > fCutsRD[GetGlobalIndex(6,ptbin)] ||
       TMath::Abs(dd->Getd0Prong(1)) > fCutsRD[GetGlobalIndex(5,ptbin)]) okD0bar = 0;
    if(!okD0 && !okD0bar) return 0;
    
    if(dd->GetDCA() > fCutsRD[GetGlobalIndex(1,ptbin)]) return 0;
    
    dd->InvMassD0(mD0,mD0bar);
    if(TMath::Abs(mD0-mD0PDG) > fCutsRD[GetGlobalIndex(0,ptbin)]) okD0 = 0;
    if(TMath::Abs(mD0bar-mD0PDG) > fCutsRD[GetGlobalIndex(0,ptbin)]) okD0bar = 0;
    if(!okD0 && !okD0bar) return 0;
    
    dd->CosThetaStarD0(ctsD0,ctsD0bar);
    if(TMath::Abs(ctsD0) > fCutsRD[GetGlobalIndex(2,ptbin)]) okD0 = 0;
    if(TMath::Abs(ctsD0bar) > fCutsRD[GetGlobalIndex(2,ptbin)]) okD0bar = 0;
    if(!okD0 && !okD0bar) return 0;
    
    if(dd->Prodd0d0() > fCutsRD[GetGlobalIndex(7,ptbin)]) return 0;
    
    if(dd->CosPointingAngle() < fCutsRD[GetGlobalIndex(8,ptbin)]) return 0;

    if (okD0) returnvalue=1; //cuts passed as D0
    if (okD0bar) returnvalue=2; //cuts passed as D0bar
    if (okD0 && okD0bar) returnvalue=3; //both
  }

  return returnvalue;
}


