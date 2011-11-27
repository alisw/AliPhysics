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
// Class for cuts on AOD reconstructed D0->Kpi
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>

#include "AliRDHFCutsD0toKpi.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAODPid.h"
#include "AliTPCPIDResponse.h"
#include "AliAODVertex.h"
#include "AliKFVertex.h"
#include "AliKFParticle.h"

ClassImp(AliRDHFCutsD0toKpi)

//--------------------------------------------------------------------------
AliRDHFCutsD0toKpi::AliRDHFCutsD0toKpi(const char* name) : 
AliRDHFCuts(name),
fUseSpecialCuts(kFALSE),
fLowPt(kTRUE),
fDefaultPID(kFALSE),
fUseKF(kFALSE),
fPtLowPID(2.)
{
  //
  // Default Constructor
  //
  Int_t nvars=11;
  SetNVars(nvars);
  TString varNames[11]={"inv. mass [GeV]",   
			"dca [cm]",
			"cosThetaStar", 
			"pTK [GeV/c]",
			"pTPi [GeV/c]",
			"d0K [cm]",
			"d0Pi [cm]",
			"d0d0 [cm^2]",
			"cosThetaPoint",
			"|cosThetaPointXY|", 
			"NormDecayLenghtXY"};
  Bool_t isUpperCut[11]={kTRUE,
			 kTRUE,
			 kTRUE,
			 kFALSE,
			 kFALSE,
			 kTRUE,
			 kTRUE,
			 kTRUE,
			 kFALSE,
			 kFALSE, 
			 kFALSE};
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[11]={kFALSE,
		     kTRUE,
		     kTRUE,
		     kFALSE,
		     kFALSE,
		     kFALSE,
		     kFALSE,
		     kTRUE,
		     kTRUE,
		     kFALSE,
		     kFALSE};
  SetVarsForOpt(4,forOpt);
  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);

}
//--------------------------------------------------------------------------
AliRDHFCutsD0toKpi::AliRDHFCutsD0toKpi(const AliRDHFCutsD0toKpi &source) :
  AliRDHFCuts(source),   
  fUseSpecialCuts(source.fUseSpecialCuts),
  fLowPt(source.fLowPt),
  fDefaultPID(source.fDefaultPID),
  fUseKF(source.fUseKF),
  fPtLowPID(source.fPtLowPID)
{
  //
  // Copy constructor
  //

}
//--------------------------------------------------------------------------
AliRDHFCutsD0toKpi &AliRDHFCutsD0toKpi::operator=(const AliRDHFCutsD0toKpi &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliRDHFCuts::operator=(source); 
  fUseSpecialCuts=source.fUseSpecialCuts;
  fLowPt=source.fLowPt;
  fDefaultPID=source.fDefaultPID;
  fPtLowPID=source.fPtLowPID;

  return *this;
}


//---------------------------------------------------------------------------
void AliRDHFCutsD0toKpi::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters,AliAODEvent *aod) {
  // 
  // Fills in vars the values of the variables 
  //

  if(nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsD0toKpi::GetCutsVarsForOpt: wrong number of variables\n");
    return;
  }

  AliAODRecoDecayHF2Prong *dd = (AliAODRecoDecayHF2Prong*)d;
 
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
		vars[iter]=TMath::Abs(dd->CosPointingAngleXY());
	}
  
   if(fVarsForOpt[10]){
		iter++;
	   vars[iter]=dd->NormalizedDecayLengthXY();
	}

   if(cleanvtx)CleanOwnPrimaryVtx(dd,aod,origownvtx);

  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpi::IsSelected(TObject* obj,Int_t selectionLevel,AliAODEvent* aod) {
  //
  // Apply selection
  //
 

  fIsSelectedCuts=0;
  fIsSelectedPID=0;

  if(!fCutsRD){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  //PrintAll();
  AliAODRecoDecayHF2Prong* d=(AliAODRecoDecayHF2Prong*)obj;

  if(!d){
    cout<<"AliAODRecoDecayHF2Prong null"<<endl;
    return 0;
  }

  if(fKeepSignalMC) if(IsSignalMC(d,aod,421)) return 3;

  Double_t ptD=d->Pt();
  if(ptD<fMinPtCand) return 0;
  if(ptD>fMaxPtCand) return 0;

  if(d->HasBadDaughters()) return 0;

  // returnvalue: 0 not sel, 1 only D0, 2 only D0bar, 3 both
  Int_t returnvaluePID=3;
  Int_t returnvalueCuts=3;

  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {

    if(!fUseKF) {

      //recalculate vertex w/o daughters
      AliAODVertex *origownvtx=0x0;
      if(fRemoveDaughtersFromPrimary && !fUseMCVertex) {
	if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
	if(!RecalcOwnPrimaryVtx(d,aod)) { 
	  CleanOwnPrimaryVtx(d,aod,origownvtx);
	  return 0;
	}
      }

      if(fUseMCVertex) {
	if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
	if(!SetMCPrimaryVtx(d,aod)) {
	  CleanOwnPrimaryVtx(d,aod,origownvtx);
	  return 0;
	}
      }
      
      Double_t pt=d->Pt();
      
      Int_t okD0=0,okD0bar=0;
      
      Int_t ptbin=PtBin(pt);
      if (ptbin==-1) {
	CleanOwnPrimaryVtx(d,aod,origownvtx);
	return 0;
      }

      Double_t mD0,mD0bar,ctsD0,ctsD0bar;
      okD0=1; okD0bar=1;
      
      Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

      d->InvMassD0(mD0,mD0bar);
      if(TMath::Abs(mD0-mD0PDG) > fCutsRD[GetGlobalIndex(0,ptbin)]) okD0 = 0;
      if(TMath::Abs(mD0bar-mD0PDG) > fCutsRD[GetGlobalIndex(0,ptbin)])  okD0bar = 0;
      if(!okD0 && !okD0bar)  {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}

      if(d->Prodd0d0() > fCutsRD[GetGlobalIndex(7,ptbin)])  {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}
    
      
      if(d->Pt2Prong(1) < fCutsRD[GetGlobalIndex(3,ptbin)]*fCutsRD[GetGlobalIndex(3,ptbin)] || d->Pt2Prong(0) < fCutsRD[GetGlobalIndex(4,ptbin)]*fCutsRD[GetGlobalIndex(4,ptbin)]) okD0 = 0;
      if(d->Pt2Prong(0) < fCutsRD[GetGlobalIndex(3,ptbin)]*fCutsRD[GetGlobalIndex(3,ptbin)] || d->Pt2Prong(1) < fCutsRD[GetGlobalIndex(4,ptbin)]*fCutsRD[GetGlobalIndex(4,ptbin)]) okD0bar = 0;
      if(!okD0 && !okD0bar) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}
      
      
      if(TMath::Abs(d->Getd0Prong(1)) > fCutsRD[GetGlobalIndex(5,ptbin)] || 
	 TMath::Abs(d->Getd0Prong(0)) > fCutsRD[GetGlobalIndex(6,ptbin)]) okD0 = 0;
      if(TMath::Abs(d->Getd0Prong(0)) > fCutsRD[GetGlobalIndex(6,ptbin)] ||
	 TMath::Abs(d->Getd0Prong(1)) > fCutsRD[GetGlobalIndex(5,ptbin)]) okD0bar = 0;
      if(!okD0 && !okD0bar)  {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}
      
      if(d->GetDCA() > fCutsRD[GetGlobalIndex(1,ptbin)])  {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}
      
    
      d->CosThetaStarD0(ctsD0,ctsD0bar);
      if(TMath::Abs(ctsD0) > fCutsRD[GetGlobalIndex(2,ptbin)]) okD0 = 0; 
      if(TMath::Abs(ctsD0bar) > fCutsRD[GetGlobalIndex(2,ptbin)]) okD0bar = 0;
      if(!okD0 && !okD0bar)   {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}
    
      if(d->CosPointingAngle() < fCutsRD[GetGlobalIndex(8,ptbin)])  {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}

      
      if(TMath::Abs(d->CosPointingAngleXY()) < fCutsRD[GetGlobalIndex(9,ptbin)])  {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}
	
      Double_t normalDecayLengXY=d->NormalizedDecayLengthXY();
      if (normalDecayLengXY < fCutsRD[GetGlobalIndex(10, ptbin)]) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}
      
      if (returnvalueCuts!=0) {
	if (okD0) returnvalueCuts=1; //cuts passed as D0
	if (okD0bar) returnvalueCuts=2; //cuts passed as D0bar
	if (okD0 && okD0bar) returnvalueCuts=3; //cuts passed as D0 and D0bar
      }

      // call special cuts
      Int_t special=1;
      if(fUseSpecialCuts) special=IsSelectedSpecialCuts(d);
      if(!special) {CleanOwnPrimaryVtx(d,aod,origownvtx); return 0;}

      // unset recalculated primary vertex when not needed any more
      CleanOwnPrimaryVtx(d,aod,origownvtx);

    } else {
      // go to selection with Kalman vertexing, if requested
      returnvalueCuts = IsSelectedKF(d,aod);
    }

    fIsSelectedCuts=returnvalueCuts;
    if(!returnvalueCuts) return 0;
  }


  // selection on PID 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate ||
     selectionLevel==AliRDHFCuts::kPID) {
    returnvaluePID = IsSelectedPID(d);
    fIsSelectedPID=returnvaluePID;
    if(!returnvaluePID) return 0;
  }

  Int_t returnvalueComb=CombineSelectionLevels(3,returnvalueCuts,returnvaluePID);

  if(!returnvalueComb) return 0;

  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(d)) return 0;
  }
 
  //  cout<<"Pid = "<<returnvaluePID<<endl;
  return returnvalueComb;
}

//------------------------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpi::IsSelectedKF(AliAODRecoDecayHF2Prong *d,
				       AliAODEvent* aod) const {
  //
  // Apply selection using KF-vertexing
  //
	
  AliAODTrack *track0 = (AliAODTrack*)d->GetDaughter(0);
  AliAODTrack *track1 = (AliAODTrack*)d->GetDaughter(1); 
   
  if(!track0 || !track1) {
    cout<<"one or two D0 daughters missing!"<<endl;
    return 0;
  }

  // returnvalue: 0 not sel, 1 only D0, 2 only D0bar, 3 both
  Int_t returnvalueCuts=3;
	 
  // candidate selection with AliKF
  AliKFParticle::SetField(aod->GetMagneticField()); // set the magnetic field
  
  Int_t okD0=0,okD0bar=0;
  okD0=1; okD0bar=1;
  
  // convert tracks into AliKFParticles
  
  AliKFParticle negPiKF(*track1,-211); // neg pion kandidate
  AliKFParticle negKKF(*track1,-321); // neg kaon kandidate
  AliKFParticle posPiKF(*track0,211); // pos pion kandidate
  AliKFParticle posKKF(*track0,321); // pos kaon kandidate
  
  // build D0 candidates
  
  AliKFParticle d0c(negKKF,posPiKF); // D0 candidate
  AliKFParticle ad0c(posKKF,negPiKF); // D0bar candidate
  
  // create kf primary vertices
  
  AliAODVertex *vtx1 = aod->GetPrimaryVertex();
  AliKFVertex primVtx(*vtx1); 
  AliKFVertex aprimVtx(*vtx1);
  
  if(primVtx.GetNContributors()<=0) okD0 = 0;
  if(aprimVtx.GetNContributors()<=0) okD0bar = 0;
  if(!okD0 && !okD0bar) returnvalueCuts=0;
	
  // calculate mass
	
  Double_t d0mass = d0c.GetMass();
  Double_t ad0mass = ad0c.GetMass();
	
  // calculate P of D0 and D0bar
  Double_t d0P = d0c.GetP();
  Double_t d0Px = d0c.GetPx();
  Double_t d0Py = d0c.GetPy();
  Double_t d0Pz = d0c.GetPz();
  Double_t ad0P = ad0c.GetP(); 
  Double_t ad0Px = ad0c.GetPx();
  Double_t ad0Py = ad0c.GetPy();
  Double_t ad0Pz = ad0c.GetPz();
  
  //calculate Pt of D0 and D0bar
	
  Double_t pt=d0c.GetPt(); 
  Double_t apt=ad0c.GetPt();
	
  // remove D0 daughters from primary vertices (if used in vertex fit) and add D0-candidates
  
  if(track0->GetUsedForPrimVtxFit()) {
    primVtx -= posPiKF; 
    aprimVtx -= posKKF;
  }
  
  if(track1->GetUsedForPrimVtxFit()) { 
    primVtx -= negKKF; 
    aprimVtx -= negPiKF;
  }
  
  primVtx += d0c;
  aprimVtx += ad0c;
  
  if(primVtx.GetNContributors()<=0) okD0 = 0;
  if(aprimVtx.GetNContributors()<=0) okD0bar = 0;
  if(!okD0 && !okD0bar) returnvalueCuts=0;
  
  //calculate cut variables
  
  // calculate impact params of daughters w.r.t recalculated vertices
  
  Double_t impactPi = posPiKF.GetDistanceFromVertexXY(primVtx);
  Double_t aimpactPi = negPiKF.GetDistanceFromVertexXY(aprimVtx);
  Double_t impactKa = negKKF.GetDistanceFromVertexXY(primVtx);
  Double_t aimpactKa = posKKF.GetDistanceFromVertexXY(aprimVtx);
	
  // calculate Product of Impact Params
	
  Double_t prodParam = impactPi*impactKa;
  Double_t aprodParam = aimpactPi*aimpactKa;
	
  // calculate cosine of pointing angles
	
  TVector3 mom(d0c.GetPx(),d0c.GetPy(),d0c.GetPz());
  TVector3 fline(d0c.GetX()-primVtx.GetX(),
		 d0c.GetY()-primVtx.GetY(),
		 d0c.GetZ()-primVtx.GetZ());
  Double_t pta = mom.Angle(fline);
  Double_t cosP = TMath::Cos(pta); // cosine of pta for D0 candidate
  
  TVector3 amom(ad0c.GetPx(),ad0c.GetPy(),ad0c.GetPz());
  TVector3 afline(ad0c.GetX()-aprimVtx.GetX(),
		  ad0c.GetY()-aprimVtx.GetY(),
		  ad0c.GetZ()-aprimVtx.GetZ());
  Double_t apta = amom.Angle(afline);
  Double_t acosP = TMath::Cos(apta); // cosine of pta for D0bar candidate
  
  // calculate P of Pions at Decay Position of D0 and D0bar candidates
  negKKF.TransportToParticle(d0c);
  posPiKF.TransportToParticle(d0c);
  posKKF.TransportToParticle(ad0c);
  negPiKF.TransportToParticle(ad0c);
  
  Double_t pxPi =  posPiKF.GetPx();
  Double_t pyPi =  posPiKF.GetPy();
  Double_t pzPi =  posPiKF.GetPz();
  Double_t ptPi =  posPiKF.GetPt();
  
  Double_t apxPi =  negPiKF.GetPx();
  Double_t apyPi =  negPiKF.GetPy();
  Double_t apzPi =  negPiKF.GetPz();
  Double_t aptPi =  negPiKF.GetPt();
  
  // calculate Pt of Kaons at Decay Position of D0 and D0bar candidates
  
  Double_t ptK =  negKKF.GetPt();
  Double_t aptK =  posKKF.GetPt();
	
  //calculate cos(thetastar)
  Double_t massvtx = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t massp[2];
  massp[0] = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  massp[1] = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  Double_t pStar = TMath::Sqrt(TMath::Power(massvtx*massvtx-massp[0]*massp[0]-massp[1]*massp[1],2.)
			       -4.*massp[0]*massp[0]*massp[1]*massp[1])/(2.*massvtx);
  
  // cos(thetastar) for D0 and Pion
  
  Double_t d0E = TMath::Sqrt(massvtx*massvtx + d0P*d0P);
  Double_t beta = d0P/d0E;
  Double_t gamma = d0E/massvtx;
  TVector3 momPi(pxPi,pyPi,pzPi);
  TVector3 momTot(d0Px,d0Py,d0Pz);
  Double_t q1 = momPi.Dot(momTot)/momTot.Mag();
  Double_t cts = (q1/gamma-beta*TMath::Sqrt(pStar*pStar+massp[1]*massp[1]))/pStar;	
    	
  // cos(thetastar) for D0bar and Pion	
  
  Double_t ad0E = TMath::Sqrt(massvtx*massvtx + ad0P*ad0P);
  Double_t abeta = ad0P/ad0E;
  Double_t agamma = ad0E/massvtx;
  TVector3 amomPi(apxPi,apyPi,apzPi);
  TVector3 amomTot(ad0Px,ad0Py,ad0Pz);
  Double_t aq1 = amomPi.Dot(amomTot)/amomTot.Mag();
  Double_t acts = (aq1/agamma-abeta*TMath::Sqrt(pStar*pStar+massp[1]*massp[1]))/pStar;	
  
  // calculate reduced Chi2 for the full D0 fit
  d0c.SetProductionVertex(primVtx);
  ad0c.SetProductionVertex(aprimVtx);
  negKKF.SetProductionVertex(d0c);
  posPiKF.SetProductionVertex(d0c);
  posKKF.SetProductionVertex(ad0c);
  negPiKF.SetProductionVertex(ad0c);
  d0c.TransportToProductionVertex();
  ad0c.TransportToProductionVertex();
	
  // calculate the decay length
  Double_t decayLengthD0 = d0c.GetDecayLength();
  Double_t adecayLengthD0 = ad0c.GetDecayLength();
  
  Double_t chi2D0 = 50.;
  if(d0c.GetNDF() > 0 && d0c.GetChi2() >= 0) {
    chi2D0 = d0c.GetChi2()/d0c.GetNDF();
  }
  
  Double_t achi2D0 = 50.;
  if(ad0c.GetNDF() > 0 && ad0c.GetChi2() >= 0) {
    achi2D0 = ad0c.GetChi2()/ad0c.GetNDF();
  }
	
  // Get the Pt-bins
  Int_t ptbin=PtBin(pt);
  Int_t aptbin=PtBin(apt);

  if(ptbin < 0) okD0 = 0;
  if(aptbin < 0) okD0bar = 0;
  if(!okD0 && !okD0bar) returnvalueCuts=0;
  
  if(ptK < fCutsRD[GetGlobalIndex(3,ptbin)] || ptPi < fCutsRD[GetGlobalIndex(4,ptbin)]) okD0 = 0;
  if(aptK < fCutsRD[GetGlobalIndex(3,aptbin)] || aptPi < fCutsRD[GetGlobalIndex(4,aptbin)]) okD0bar = 0;
  if(!okD0 && !okD0bar) returnvalueCuts=0;
  
  
  if(TMath::Abs(impactKa) > fCutsRD[GetGlobalIndex(5,ptbin)] || 
     TMath::Abs(impactPi) > fCutsRD[GetGlobalIndex(6,ptbin)]) okD0 = 0;
  
  if(TMath::Abs(aimpactKa) > fCutsRD[GetGlobalIndex(5,aptbin)] ||
     TMath::Abs(aimpactPi) > fCutsRD[GetGlobalIndex(6,aptbin)]) okD0bar = 0;
  
  if(!okD0 && !okD0bar)  returnvalueCuts=0;
  
  // for the moment via the standard method due to bug in AliKF
  if(d->GetDCA() > fCutsRD[GetGlobalIndex(1,ptbin)]) okD0 = 0;
  if(d->GetDCA() > fCutsRD[GetGlobalIndex(1,aptbin)]) okD0bar = 0;
  if(!okD0 && !okD0bar) returnvalueCuts=0;
    
    
  if(TMath::Abs(d0mass-massvtx) > fCutsRD[GetGlobalIndex(0,ptbin)]) okD0 = 0;
  if(TMath::Abs(ad0mass-massvtx) > fCutsRD[GetGlobalIndex(0,aptbin)])  okD0bar = 0;
  if(!okD0 && !okD0bar)  returnvalueCuts=0;
  
  
  if(TMath::Abs(cts) > fCutsRD[GetGlobalIndex(2,ptbin)]) okD0 = 0; 
  if(TMath::Abs(acts) > fCutsRD[GetGlobalIndex(2,aptbin)]) okD0bar = 0;
  if(!okD0 && !okD0bar)   returnvalueCuts=0;
  
  if(prodParam  > fCutsRD[GetGlobalIndex(7,ptbin)]) okD0 = 0;
  if(aprodParam > fCutsRD[GetGlobalIndex(7,aptbin)]) okD0bar = 0;
  if(!okD0 && !okD0bar)  returnvalueCuts=0;
    
  if(cosP  < fCutsRD[GetGlobalIndex(8,ptbin)]) okD0 = 0; 
  if(acosP < fCutsRD[GetGlobalIndex(8,aptbin)]) okD0bar = 0;
  if(!okD0 && !okD0bar)  returnvalueCuts=0;
	
  if(chi2D0  > fCutsRD[GetGlobalIndex(10,ptbin)]) okD0 = 0; 
  if(achi2D0 > fCutsRD[GetGlobalIndex(10,aptbin)]) okD0bar = 0;
  if(!okD0 && !okD0bar)  returnvalueCuts=0;
	
  if(decayLengthD0  < fCutsRD[GetGlobalIndex(9,ptbin)]) okD0 = 0; 
  if(adecayLengthD0 < fCutsRD[GetGlobalIndex(9,aptbin)]) okD0bar = 0;
  if(!okD0 && !okD0bar)  returnvalueCuts=0;
    
  if(returnvalueCuts!=0) {
    if(okD0) returnvalueCuts=1; //cuts passed as D0
    if(okD0bar) returnvalueCuts=2; //cuts passed as D0bar
    if(okD0 && okD0bar) returnvalueCuts=3; //cuts passed as D0 and D0bar
  }

  return returnvalueCuts;  
}

//---------------------------------------------------------------------------

Bool_t AliRDHFCutsD0toKpi::IsInFiducialAcceptance(Double_t pt, Double_t y) const
{
  //
  // Checking if D0 is in fiducial acceptance region 
  //

  if(pt > 5.) {
    // applying cut for pt > 5 GeV
    AliDebug(2,Form("pt of D0 = %f (> 5), cutting at |y| < 0.8\n",pt)); 
    if (TMath::Abs(y) > 0.8){
      return kFALSE;
    }
  } else {
    // appliying smooth cut for pt < 5 GeV
    Double_t maxFiducialY = -0.2/15*pt*pt+1.9/15*pt+0.5; 
    Double_t minFiducialY = 0.2/15*pt*pt-1.9/15*pt-0.5;		
    AliDebug(2,Form("pt of D0 = %f (< 5), cutting  according to the fiducial zone [%f, %f]\n",pt,minFiducialY,maxFiducialY)); 
    if (y < minFiducialY || y > maxFiducialY){
      return kFALSE;
    }
  }

  return kTRUE;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpi::IsSelectedPID(AliAODRecoDecayHF* d) 
{
  // ############################################################
  //
  // Apply PID selection
  //
  //
  // ############################################################

  if(!fUsePID) return 3;
  if(fDefaultPID) return IsSelectedPIDdefault(d);
  fWhyRejection=0;
  Int_t isD0D0barPID[2]={1,2};
  Int_t combinedPID[2][2];// CONVENTION: [daught][isK,IsPi]; [0][0]=(prong 1, isK)=value [0][1]=(prong 1, isPi)=value; 
  //                                                                                                 same for prong 2
  //                                               values convention -1 = discarded 
  //                                                                  0 = not identified (but compatible) || No PID (->hasPID flag)
  //                                                                  1 = identified
  // PID search:   pion (TPC) or not K (TOF), Kaon hypothesis for both 
  // Initial hypothesis: unknwon (but compatible) 
  combinedPID[0][0]=0;  // prima figlia, Kaon
  combinedPID[0][1]=0;  // prima figlia, pione
  combinedPID[1][0]=0;  // seconda figlia, Kaon
  combinedPID[1][1]=0;  // seconda figlia, pion
  
  Bool_t checkPIDInfo[2]={kTRUE,kTRUE};
  Double_t sigma_tmp[3]={fPidHF->GetSigma(0),fPidHF->GetSigma(1),fPidHF->GetSigma(2)};
  for(Int_t daught=0;daught<2;daught++){
    //Loop con prongs
    AliAODTrack *aodtrack=(AliAODTrack*)d->GetDaughter(daught);
    if(fPidHF->IsTOFPiKexcluded(aodtrack,5.)) return 0; 
    
    if(!(fPidHF->CheckStatus(aodtrack,"TPC")) && !(fPidHF->CheckStatus(aodtrack,"TOF"))) {
      checkPIDInfo[daught]=kFALSE; 
      continue;
    }

    // identify kaon
    combinedPID[daught][0]=fPidHF->MakeRawPid(aodtrack,3);

    // identify pion

    if(!(fPidHF->CheckStatus(aodtrack,"TPC"))) {
     combinedPID[daught][1]=0;
    }else{
      fPidHF->SetTOF(kFALSE);
      combinedPID[daught][1]=fPidHF->MakeRawPid(aodtrack,2);
      fPidHF->SetTOF(kTRUE);
      fPidHF->SetCompat(kTRUE);
     }


    if(combinedPID[daught][0]<=-1&&combinedPID[daught][1]<=-1){ // if not a K- and not a pi- both D0 and D0bar excluded
      isD0D0barPID[0]=0;
      isD0D0barPID[1]=0;
    }
    else if(combinedPID[daught][0]==2&&combinedPID[daught][1]>=1){
      if(aodtrack->Charge()==-1)isD0D0barPID[1]=0;//if K- D0bar excluded
      else isD0D0barPID[0]=0;// if K+ D0 excluded
    }
    /*    else if(combinedPID[daught][0]==1&&combinedPID[daught][1]>=1){
	  isD0D0barPID[0]=0;
	  isD0D0barPID[1]=0;
	  }
    */
    else if(combinedPID[daught][0]>=1||combinedPID[daught][1]<=-1){ 
      if(aodtrack->Charge()==-1)isD0D0barPID[1]=0;// not a D0bar if K- or if pi- excluded
      else isD0D0barPID[0]=0;//  not a D0 if K+ or if pi+ excluded
    }
    else if(combinedPID[daught][0]<=-1||combinedPID[daught][1]>=1){
      if(aodtrack->Charge()==-1)isD0D0barPID[0]=0;// not a D0 if pi- or if K- excluded
      else isD0D0barPID[1]=0;// not a D0bar if pi+ or if K+ excluded
    }

    if(fLowPt && d->Pt()<fPtLowPID){
     Double_t sigmaTPC[3]={3.,2.,0.};
     fPidHF->SetSigmaForTPC(sigmaTPC);
    // identify kaon
    combinedPID[daught][0]=fPidHF->MakeRawPid(aodtrack,3);

    Double_t ptProng=aodtrack->P();

    if(ptProng<0.6){
     fPidHF->SetCompat(kFALSE);
     combinedPID[daught][0]=fPidHF->MakeRawPid(aodtrack,3);
     fPidHF->SetCompat(kTRUE);
    }

    if(!(fPidHF->CheckStatus(aodtrack,"TPC"))) {
     combinedPID[daught][1]=0;
    }else{
      fPidHF->SetTOF(kFALSE);
      Double_t sigmaTPCpi[3]={3.,3.,0.};
      fPidHF->SetSigmaForTPC(sigmaTPCpi);
      combinedPID[daught][1]=fPidHF->MakeRawPid(aodtrack,2);
      fPidHF->SetTOF(kTRUE);
       if(ptProng<0.8){
        Bool_t isTPCpion=fPidHF->IsPionRaw(aodtrack,"TPC");
        if(isTPCpion){
         combinedPID[daught][1]=1;
        }else{
         combinedPID[daught][1]=-1;
        }
      }
    }

   }
   fPidHF->SetSigmaForTPC(sigma_tmp);
  }// END OF LOOP ON DAUGHTERS

   if(!checkPIDInfo[0] && !checkPIDInfo[1]) {
    if(fLowPt) fPidHF->SetSigmaForTPC(sigma_tmp);
    return 0;
   }


  // FURTHER PID REQUEST (both daughter info is needed)
  if(combinedPID[0][0]<=-1&&combinedPID[1][0]<=-1){
    fWhyRejection=31;// reject cases in which no kaon-compatible tracks are found
    if(fLowPt) fPidHF->SetSigmaForTPC(sigma_tmp);
    return 0;
  }

  if(fLowPt && d->Pt()<fPtLowPID){    
    if(combinedPID[0][0]<=0&&combinedPID[1][0]<=0){
      fWhyRejection=32;// reject cases where the Kaon is not identified
      fPidHF->SetSigmaForTPC(sigma_tmp);
      return 0;
    }
  }
    if(fLowPt) fPidHF->SetSigmaForTPC(sigma_tmp);

  //  cout<<"Why? "<<fWhyRejection<<endl;  
  return isD0D0barPID[0]+isD0D0barPID[1];
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpi::IsSelectedPIDdefault(AliAODRecoDecayHF* d) 
{
  // ############################################################
  //
  // Apply PID selection
  //
  //
  // temporary selection: PID AS USED FOR D0 by Andrea Rossi (up to 28/06/2010)
  //
  // d must be a AliAODRecoDecayHF2Prong object
  // returns 0 if both D0 and D0bar are rejectecd
  //         1 if D0 is accepted while D0bar is rejected
  //         2 if D0bar is accepted while D0 is rejected
  //         3 if both are accepted
  // fWhyRejection variable (not returned for the moment, print it if needed)
  //               keeps some information on why a candidate has been 
  //               rejected according to the following (unfriendly?) scheme 
  //             if more rejection cases are considered interesting, just add numbers
  //
  //      TO BE CONSIDERED WITH A GRAIN OF SALT (the order in which cut are applied is relevant) 
  //              from 20 to 30: "detector" selection (PID acceptance)                                             
  //                                                  26: TPC refit
  //                                                  27: ITS refit
  //                                                  28: no (TOF||TPC) pid information (no kTOFpid,kTOFout,kTIME,kTPCpid,...)
  //
  //              from 30 to 40: PID selection
  //                                                  31: no Kaon compatible tracks found between daughters
  //                                                  32: no Kaon identified tracks found (strong sel. at low momenta)
  //                                                  33: both mass hypotheses are rejected 
  //                  
  // ############################################################

  if(!fUsePID) return 3;
  fWhyRejection=0;
  Int_t isD0D0barPID[2]={1,2};
  Double_t nsigmaTPCpi=-1., nsigmaTPCK=-1.; //used for TPC pid
  Double_t tofSig,times[5];// used fot TOF pid
  Int_t hasPID[2]={2,2};// flag to count how many detectors give PID info for the daughters
  Int_t isKaonPionTOF[2][2],isKaonPionTPC[2][2];
  Int_t combinedPID[2][2];// CONVENTION: [daught][isK,IsPi]; [0][0]=(prong 1, isK)=value [0][1]=(prong 1, isPi)=value; 
  //                                                                                                 same for prong 2
  //                                               values convention -1 = discarded 
  //                                                                  0 = not identified (but compatible) || No PID (->hasPID flag)
  //                                                                  1 = identified
  // PID search:   pion (TPC) or not K (TOF), Kaon hypothesis for both 
  // Initial hypothesis: unknwon (but compatible) 
  isKaonPionTOF[0][0]=0;
  isKaonPionTOF[0][1]=0;
  isKaonPionTOF[1][0]=0;
  isKaonPionTOF[1][1]=0;
  
  isKaonPionTPC[0][0]=0;
  isKaonPionTPC[0][1]=0;
  isKaonPionTPC[1][0]=0;
  isKaonPionTPC[1][1]=0;
  
  combinedPID[0][0]=0;
  combinedPID[0][1]=0;
  combinedPID[1][0]=0;
  combinedPID[1][1]=0;
  
  
 
  for(Int_t daught=0;daught<2;daught++){
    //Loop con prongs
    
    // ########### Step 0- CHECKING minimal PID "ACCEPTANCE" ####################

    AliAODTrack *aodtrack=(AliAODTrack*)d->GetDaughter(daught); 
   
    if(!(aodtrack->GetStatus()&AliESDtrack::kTPCrefit)){
      fWhyRejection=26;
      return 0;
    } 
    if(!(aodtrack->GetStatus()&AliESDtrack::kITSrefit)){
      fWhyRejection=27;
      return 0;
    } 
    
    AliAODPid *pid=aodtrack->GetDetPid();
    if(!pid) {
      //delete esdtrack;
      hasPID[daught]--;
      continue;
    }
  
    // ########### Step 1- Check of TPC and TOF response ####################

    Double_t ptrack=aodtrack->P();
    //#################### TPC PID #######################
     if (!(aodtrack->GetStatus()&AliESDtrack::kTPCpid )){
       // NO TPC PID INFO FOR THIS TRACK 
       hasPID[daught]--;
     }
     else {
       static AliTPCPIDResponse theTPCpid;
       AliAODPid *pidObj = aodtrack->GetDetPid();
       Double_t ptProng=pidObj->GetTPCmomentum();
       nsigmaTPCpi = theTPCpid.GetNumberOfSigmas(ptProng,(Float_t)pid->GetTPCsignal(),(Int_t)aodtrack->GetTPCClusterMap().CountBits(),AliPID::kPion);
       nsigmaTPCK =  theTPCpid.GetNumberOfSigmas(ptProng,(Float_t)pid->GetTPCsignal(),(Int_t)aodtrack->GetTPCClusterMap().CountBits(),AliPID::kKaon);
       //if(ptrack<0.6){
       if(ptProng<0.6){
	 if(TMath::Abs(nsigmaTPCK)<2.)isKaonPionTPC[daught][0]=1;
	 else if(TMath::Abs(nsigmaTPCK)>3.)isKaonPionTPC[daught][0]=-1;
	 if(TMath::Abs(nsigmaTPCpi)<2.)isKaonPionTPC[daught][1]=1;
	 else if(TMath::Abs(nsigmaTPCpi)>3.)isKaonPionTPC[daught][1]=-1;
       }
       //else if(ptrack<.8){
       else if(ptProng<.8){
	 if(TMath::Abs(nsigmaTPCK)<1.)isKaonPionTPC[daught][0]=1;
	 else if(TMath::Abs(nsigmaTPCK)>3.)isKaonPionTPC[daught][0]=-1;
	 if(TMath::Abs(nsigmaTPCpi)<1.)isKaonPionTPC[daught][1]=1;
	 else if(TMath::Abs(nsigmaTPCpi)>3.)isKaonPionTPC[daught][1]=-1;
       }     
       else {
	 //	if(nsigmaTPCK>-2.&&nsigmaTPCK<1.)isKaonPionTPC[daught][0]=1;
	 if(TMath::Abs(nsigmaTPCK)>3.)isKaonPionTPC[daught][0]=-1;
	 //if(nsigmaTPCpi>-1.&&nsigmaTPCpi<2.)isKaonPionTPC[daught][1]=1;
	 if(TMath::Abs(nsigmaTPCpi)>3.)isKaonPionTPC[daught][1]=-1;
       }
     }
    
    
    // ##### TOF PID: do not ask nothing for pion/protons ############
     if(!((aodtrack->GetStatus()&AliESDtrack::kTOFpid)&&(aodtrack->GetStatus()&AliESDtrack::kTOFout)&&(aodtrack->GetStatus()&AliESDtrack::kTIME))){
       // NO TOF PID INFO FOR THIS TRACK      
       hasPID[daught]--;
     }
     else{
       tofSig=pid->GetTOFsignal(); 
       pid->GetIntegratedTimes(times);
       if((tofSig-times[3])>5.*160.)return 0;// PROTON REJECTION
       if(TMath::Abs(tofSig-times[3])>3.*160.){
	 isKaonPionTOF[daught][0]=-1;
       }
       else {	 
	 if(ptrack<1.5){
	   isKaonPionTOF[daught][0]=1;
	 }
       }
     }
     
     //######### Step 2: COMBINE TOF and TPC PID ###############
     // we apply the following convention: if TPC and TOF disagree (discarded Vs identified) -> unknown
     combinedPID[daught][0]=isKaonPionTOF[daught][0]+isKaonPionTPC[daught][0];
     combinedPID[daught][1]=isKaonPionTOF[daught][1]+isKaonPionTPC[daught][1];
     
     
     //######### Step 3:   USE PID INFO     
     
     if(combinedPID[daught][0]<=-1&&combinedPID[daught][1]<=-1){// if not a K- and not a pi- both D0 and D0bar excluded
       isD0D0barPID[0]=0;
       isD0D0barPID[1]=0;
     }
     else if(combinedPID[daught][0]==2&&combinedPID[daught][1]>=1){// if in conflict (both pi- and K-), if k for both TPC and TOF -> is K
       if(aodtrack->Charge()==-1)isD0D0barPID[1]=0;//if K- D0bar excluded
       else isD0D0barPID[0]=0;// if K+ D0 excluded
     }
     else if(combinedPID[daught][0]==1&&combinedPID[daught][1]>=1){// if in conflict (both pi- and K-) and k- only for TPC or TOF -> reject
       isD0D0barPID[0]=0;
       isD0D0barPID[1]=0;
     }
     else if(combinedPID[daught][0]>=1||combinedPID[daught][1]<=-1){
       if(aodtrack->Charge()==-1)isD0D0barPID[1]=0;// not a D0bar if K- or if pi- excluded
       else isD0D0barPID[0]=0;//  not a D0 if K+ or if pi+ excluded
     }
     else if(combinedPID[daught][0]<=-1||combinedPID[daught][1]>=1){
       if(aodtrack->Charge()==-1)isD0D0barPID[0]=0;// not a D0 if pi- or if K- excluded
      else isD0D0barPID[1]=0;// not a D0bar if pi+ or if K+ excluded
     }
     
     // ##########  ALSO DIFFERENT TPC PID REQUEST FOR LOW pt D0: request of K identification      ###############################
     // ########## more tolerant criteria for single particle ID-> more selective criteria for D0   ##############################
     // ###############                     NOT OPTIMIZED YET                                  ###################################
     if(d->Pt()<2.){
       isKaonPionTPC[daught][0]=0;
       isKaonPionTPC[daught][1]=0;
       AliAODPid *pidObj = aodtrack->GetDetPid();
       Double_t ptProng=pidObj->GetTPCmomentum();
       //if(ptrack<0.6){
       if(ptProng<0.6){
	 if(TMath::Abs(nsigmaTPCK)<3.)isKaonPionTPC[daught][0]=1;
	 else if(TMath::Abs(nsigmaTPCK)>3.)isKaonPionTPC[daught][0]=-1;
	 if(TMath::Abs(nsigmaTPCpi)<3.)isKaonPionTPC[daught][1]=1;
	 else if(TMath::Abs(nsigmaTPCpi)>3.)isKaonPionTPC[daught][1]=-1;
     }
       //else if(ptrack<.8){
       else if(ptProng<.8){
	 if(TMath::Abs(nsigmaTPCK)<2.)isKaonPionTPC[daught][0]=1;
	 else if(TMath::Abs(nsigmaTPCK)>3.)isKaonPionTPC[daught][0]=-1;
	 if(TMath::Abs(nsigmaTPCpi)<3.)isKaonPionTPC[daught][1]=1;
	 else if(TMath::Abs(nsigmaTPCpi)>3.)isKaonPionTPC[daught][1]=-1;
       }     
       else {
	 if(TMath::Abs(nsigmaTPCK)>3.)isKaonPionTPC[daught][0]=-1;
	 if(TMath::Abs(nsigmaTPCpi)>3.)isKaonPionTPC[daught][1]=-1;
       }
     }
     
  }// END OF LOOP ON DAUGHTERS
  
  // FURTHER PID REQUEST (both daughter info is needed)
  if(combinedPID[0][0]<=-1&&combinedPID[1][0]<=-1){
    fWhyRejection=31;// reject cases in which no kaon-compatible tracks are found
    return 0;
  }
  else if(hasPID[0]==0&&hasPID[1]==0){
    fWhyRejection=28;// reject cases in which no PID info is available  
    return 0;
  }
  if(d->Pt()<2.){
    // request of K identification at low D0 pt
    combinedPID[0][0]=0;
    combinedPID[0][1]=0;
    combinedPID[1][0]=0;
    combinedPID[1][1]=0;
    
    combinedPID[0][0]=isKaonPionTOF[0][0]+isKaonPionTPC[0][0];
    combinedPID[0][1]=isKaonPionTOF[0][1]+isKaonPionTPC[0][1];
    combinedPID[1][0]=isKaonPionTOF[1][0]+isKaonPionTPC[1][0];
    combinedPID[1][1]=isKaonPionTOF[1][1]+isKaonPionTPC[1][1];
    
    if(combinedPID[0][0]<=0&&combinedPID[1][0]<=0){
      fWhyRejection=32;// reject cases where the Kaon is not identified
      return 0;
    }
  }

  //  cout<<"Why? "<<fWhyRejection<<endl;  
  return isD0D0barPID[0]+isD0D0barPID[1];
}



//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpi::CombineSelectionLevels(Int_t selectionvalTrack,
						 Int_t selectionvalCand,
						 Int_t selectionvalPID) const
{
  //
  // This method combines the tracks, PID and cuts selection results
  //
  if(selectionvalTrack==0) return 0;

  Int_t returnvalue;

  switch(selectionvalPID) {
  case 0:
    returnvalue=0;
    break;
  case 1:
    returnvalue=((selectionvalCand==1 || selectionvalCand==3) ? 1 : 0);
    break;
  case 2:
    returnvalue=((selectionvalCand==2 || selectionvalCand==3) ? 2 : 0);
    break;
  case 3:
    returnvalue=selectionvalCand;
    break;
  default:
    returnvalue=0;
    break;
  }

  return returnvalue;
}
//----------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpi::IsSelectedSpecialCuts(AliAODRecoDecayHF *d) const
{
  //
  // Note: this method is temporary
  // Additional cuts on decay lenght and lower cut for d0 norm are applied using vertex without candidate's daughters
  //

  //apply cuts

  Float_t normDecLengthCut=1.,decLengthCut=TMath::Min(d->P()*0.0066+0.01,0.06/*cm*/), normd0Cut=0.5;
  // "decay length" expo law with tau' = beta*gamma*ctau= p/m*ctau =p*0.0123/1.864~p*0.0066
  // decay lenght > ctau' implies to retain (1-1/e) (for signal without considering detector resolution), 

  Int_t returnvalue=3; //cut passed
  for(Int_t i=0;i<2/*prongs*/;i++){
    if(TMath::Abs(d->Normalizedd0Prong(i))<normd0Cut) return 0; //normd0Cut not passed
  }
  if(d->DecayLength2()<decLengthCut*decLengthCut)  return 0; //decLengthCut not passed
  if(d->NormalizedDecayLength2()<normDecLengthCut*normDecLengthCut)  return 0; //decLengthCut not passed
	
  return returnvalue;
}

//----------------------------------------------
void AliRDHFCutsD0toKpi::SetUseKF(Bool_t useKF)
{
  //switch on candidate selection via AliKFparticle
  if(!useKF) return;
  if(useKF){
    fUseKF=useKF;
    Int_t nvarsKF=11;
    SetNVars(nvarsKF);
    TString varNamesKF[11]={"inv. mass [GeV]",   
		            "dca [cm]",
		            "cosThetaStar", 
		            "pTK [GeV/c]",
		            "pTPi [GeV/c]",
		            "d0K [cm]",
		            "d0Pi [cm]",
		            "d0d0 [cm^2]",
		            "cosThetaPoint"
		            "DecayLength[cm]",
		            "RedChi2"};
    Bool_t isUpperCutKF[11]={kTRUE,
			     kTRUE,
			     kTRUE,
			     kFALSE,
			     kFALSE,
			     kTRUE,
			     kTRUE,
			     kTRUE,
			     kFALSE,
			     kFALSE,
			     kTRUE};
    SetVarNames(nvarsKF,varNamesKF,isUpperCutKF);
    SetGlobalIndex();
    Bool_t forOpt[11]={kFALSE,
		    kTRUE,
		    kTRUE,
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kTRUE,
		    kTRUE,
		    kFALSE,
		    kFALSE};
    SetVarsForOpt(4,forOpt);
  }
  return;
}


void AliRDHFCutsD0toKpi::SetStandardCutsPP2010() {
  //
  //STANDARD CUTS USED FOR 2010 pp analysis 
  //dca cut will be enlarged soon to 400 micron
  //
  
  SetName("D0toKpiCutsStandard");
  SetTitle("Standard Cuts for D0 analysis");
  
  // PILE UP REJECTION
  SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  // EVENT CUTS
  SetMinVtxContr(1);

  
  // TRACKS ON SINGLE TRACKS
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //  esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.3,1.e10);
  
  AddTrackCuts(esdTrackCuts);
  
  const Int_t nptbins =14;
  const Double_t ptmax = 9999.;
  const Int_t nvars=11;
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
  ptbins[12]=20.;
  ptbins[13]=24.;
  ptbins[14]=ptmax;

  SetGlobalIndex(nvars,nptbins);
  SetPtBins(nptbins+1,ptbins);
  
  Float_t cutsMatrixD0toKpiStand[nptbins][nvars]={{0.400,350.*1E-4,0.8,0.5,0.5,1000.*1E-4,1000.*1E-4,-5000.*1E-8,0.80,0.,0.},/* pt<0.5*/
						  {0.400,350.*1E-4,0.8,0.5,0.5,1000.*1E-4,1000.*1E-4,-5000.*1E-8,0.80,0.,0.},/* 0.5<pt<1*/
						  {0.400,300.*1E-4,0.8,0.4,0.4,1000.*1E-4,1000.*1E-4,-25000.*1E-8,0.80,0.,0.},/* 1<pt<2 */
						  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-15000.*1E-8,0.85,0.,0.},/* 2<pt<3 */
						  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.85,0.,0.},/* 3<pt<4 */
						  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.85,0.,0.},/* 4<pt<5 */
						  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.85,0.,0.},/* 5<pt<6 */
						  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-8000.*1E-8,0.85,0.,0.},/* 6<pt<7 */
						  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-7000.*1E-8,0.85,0.,0.},/* 7<pt<8 */
						  {0.400,300.*1E-4,0.9,0.7,0.7,1000.*1E-4,1000.*1E-4,-5000.*1E-8,0.85,0.,0.},/* 8<pt<12 */
						  {0.400,300.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,10000.*1E-8,0.85,0.,0.},/* 12<pt<16 */
						  {0.400,300.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,999999.*1E-8,0.85,0.,0.},/* 16<pt<20 */
						  {0.400,300.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,999999.*1E-8,0.85,0.,0.},/* 20<pt<24 */
						  {0.400,300.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,999999.*1E-8,0.85,0.,0.}};/* pt>24 */
  
  
  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
  
  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar = 0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];      
    }
  }
  
  SetCuts(nvars,nptbins,cutsMatrixTransposeStand);
  SetUseSpecialCuts(kTRUE);
  SetRemoveDaughtersFromPrim(kTRUE);
  
  for(Int_t iv=0;iv<nvars;iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand=NULL;

  // PID SETTINGS
  AliAODPidHF* pidObj=new AliAODPidHF();
  //pidObj->SetName("pid4D0");
  Int_t mode=1;
  const Int_t nlims=2;
  Double_t plims[nlims]={0.6,0.8}; //TPC limits in momentum [GeV/c]
  Bool_t compat=kTRUE; //effective only for this mode
  Bool_t asym=kTRUE;
  Double_t sigmas[5]={2.,1.,0.,3.,0.}; //to be checked and to be modified with new implementation of setters by Rossella
  pidObj->SetAsym(asym);// if you want to use the asymmetric bands in TPC
  pidObj->SetMatch(mode);
  pidObj->SetPLimit(plims,nlims);
  pidObj->SetSigma(sigmas);
  pidObj->SetCompat(compat);
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  pidObj->SetPCompatTOF(1.5);
  pidObj->SetSigmaForTPCCompat(3.);
  pidObj->SetSigmaForTOFCompat(3.);

  SetPidHF(pidObj);
  SetUsePID(kTRUE);
  SetUseDefaultPID(kFALSE);


  PrintAll();

  delete pidObj;
  pidObj=NULL;

  return;

}


void AliRDHFCutsD0toKpi::SetStandardCutsPbPb2010() {
  //
  //PRELIMINARY CUTS USED FOR 2010 PbPb analysis
  //... EVOLVING SOON 
  // 
  
  SetName("D0toKpiCutsStandard");
  SetTitle("Standard Cuts for D0 analysis in PbPb2010 run");
  
  // PILE UP REJECTION
  //SetOptPileup(AliRDHFCuts::kRejectPileupEvent);
  // CENTRALITY SELECTION
  SetMinCentrality(0.);
  SetMaxCentrality(80.);
  SetUseCentrality(AliRDHFCuts::kCentV0M);


  // EVENT CUTS
  SetMinVtxContr(1);
  // MAX Z-VERTEX CUT
  SetMaxVtxZ(10.);
  
  // TRACKS ON SINGLE TRACKS
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //  esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.7,1.e10);

  esdTrackCuts->SetMaxDCAToVertexXY(1.);  
  esdTrackCuts->SetMaxDCAToVertexZ(1.);
  esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0075*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");  


  AddTrackCuts(esdTrackCuts);

  // CANDIDATE CUTS  
  const Int_t nptbins =13;
  const Double_t ptmax = 9999.;
  const Int_t nvars=11;
  Float_t ptbins[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=0.5;	
  ptbins[2]=1.;
  ptbins[3]=2.;
  ptbins[4]=3.;
  ptbins[5]=4.;
  ptbins[6]=5.;
  ptbins[7]=6.;
  ptbins[8]=8.;
  ptbins[9]=12.;
  ptbins[10]=16.;
  ptbins[11]=20.;
  ptbins[12]=24.;
  ptbins[13]=ptmax;

  SetGlobalIndex(nvars,nptbins);
  SetPtBins(nptbins+1,ptbins);
  SetMinPtCandidate(2.);

  Float_t cutsMatrixD0toKpiStand[nptbins][nvars]={{0.400,400.*1E-4,0.8,0.3,0.3,1000.*1E-4,1000.*1E-4,-10000.*1E-8,0.85,0.,5.},/* pt<0.5*/
						  {0.400,400.*1E-4,0.8,0.3,0.3,1000.*1E-4,1000.*1E-4,-35000.*1E-8,0.9,0.,5.},/* 0.5<pt<1*/
						  {0.400,400.*1E-4,0.8,0.4,0.4,1000.*1E-4,1000.*1E-4,-43000.*1E-8,0.85,0.,5.},/* 1<pt<2 */
						  {0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-40000.*1E-8,0.95,0.998,5.},/* 2<pt<3 */
						  {0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-36000.*1E-8,0.95,0.998,5.},/* 3<pt<4 */
						  {0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-27000.*1E-8,0.95,0.998,5.},/* 4<pt<5 */
						  {0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-21000.*1E-8,0.92,0.998,5.},/* 5<pt<6 */
						  {0.400,270.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-14000.*1E-8,0.88,0.998,5.},/* 6<pt<8 */
						  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-5000.*1E-8,0.85,0.998,5.},/* 8<pt<12 */
						  {0.400,350.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,-1000.*1E-8,0.83,0.998,5.},/* 12<pt<16 */
						  {0.400,400.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,-1000.*1E-8,0.82,0.998,5.},/* 16<pt<20 */
						  {0.400,400.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,-1000.*1E-8,0.81,0.998,5.},/* 20<pt<24 */
						  {0.400,400.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,-1000.*1E-8,0.8,0.998,5.}};/* pt>24 */
  
  
  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
  
  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar = 0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];      
    }
  }
  
  SetCuts(nvars,nptbins,cutsMatrixTransposeStand);
  SetUseSpecialCuts(kTRUE);
  SetRemoveDaughtersFromPrim(kFALSE);// THIS IS VERY IMPORTANT! TOO SLOW IN PbPb
  for(Int_t iv=0;iv<nvars;iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand=NULL;
  
  // PID SETTINGS
  AliAODPidHF* pidObj=new AliAODPidHF();
  //pidObj->SetName("pid4D0");
  Int_t mode=1;
  const Int_t nlims=2;
  Double_t plims[nlims]={0.6,0.8}; //TPC limits in momentum [GeV/c]
  Bool_t compat=kTRUE; //effective only for this mode
  Bool_t asym=kTRUE;
  Double_t sigmas[5]={2.,1.,0.,3.,0.}; //to be checked and to be modified with new implementation of setters by Rossella
  pidObj->SetAsym(asym);// if you want to use the asymmetric bands in TPC
  pidObj->SetMatch(mode);
  pidObj->SetPLimit(plims,nlims);
  pidObj->SetSigma(sigmas);
  pidObj->SetCompat(compat);
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  pidObj->SetPCompatTOF(2.);
  pidObj->SetSigmaForTPCCompat(3.);
  pidObj->SetSigmaForTOFCompat(3.);  


  SetPidHF(pidObj);
  SetUsePID(kTRUE);
  SetUseDefaultPID(kFALSE);


  PrintAll();


  delete pidObj;
  pidObj=NULL;

  return;

}
