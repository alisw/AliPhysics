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
// Class for cuts on AOD reconstructed D0->Kpipipi
//
// Author: r.romita@gsi.de, andrea.dainese@pd.infn.it,
//	   fabio.colamaria@ba.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>

#include "AliRDHFCutsD0toKpipipi.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAODPidHF.h"

ClassImp(AliRDHFCutsD0toKpipipi)

//--------------------------------------------------------------------------
AliRDHFCutsD0toKpipipi::AliRDHFCutsD0toKpipipi(const char* name) : 
AliRDHFCuts(name)
{
  //
  // Default Constructor
  //
  Int_t nvars=9;
  SetNVars(nvars);
  TString varNames[9]={"inv. mass [GeV]",   
		       "dca [cm]",
                       "Dist 2-trk Vtx to PrimVtx [cm]",
		       "Dist 3-trk Vtx to PrimVtx [cm]",
		       "Dist 4-trk Vtx to PrimVtx [cm]",
		       "cosThetaPoint",
		       "pt [GeV/c]",
		       "rho mass [GeV]",
		       "PID cut"};
  Bool_t isUpperCut[9]={kTRUE,
			kTRUE,
			kFALSE,
			kFALSE,
			kFALSE,
			kFALSE,
			kFALSE,
			kTRUE,
			kFALSE};
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[9]={kFALSE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
		    kTRUE,
		    kFALSE,
		    kFALSE,
		    kFALSE};
  SetVarsForOpt(5,forOpt);
  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);
}
//--------------------------------------------------------------------------
AliRDHFCutsD0toKpipipi::AliRDHFCutsD0toKpipipi(const AliRDHFCutsD0toKpipipi &source) :
  AliRDHFCuts(source)
{
  //
  // Copy constructor
  //

}
//--------------------------------------------------------------------------
AliRDHFCutsD0toKpipipi &AliRDHFCutsD0toKpipipi::operator=(const AliRDHFCutsD0toKpipipi &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliRDHFCuts::operator=(source);

  return *this;
}


//---------------------------------------------------------------------------
void AliRDHFCutsD0toKpipipi::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters,AliAODEvent* aod) {
  // 
  // Fills in vars the values of the variables 
  //

  if(nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsD0toKpipipi::GetCutsVarsForOpt: wrong number of variables\n");
    return;
  }

  AliAODRecoDecayHF4Prong *dd = (AliAODRecoDecayHF4Prong*)d;

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

  if(fVarsForOpt[0]) {
    iter++;
    Double_t mD0[2],mD0bar[2];
    if(TMath::Abs(pdgdaughters[1])==321 || TMath::Abs(pdgdaughters[3])==321) {
      dd->InvMassD0(mD0);
      if(TMath::Abs(pdgdaughters[1])==321) {
       vars[iter]=mD0[0];
      }else{
       vars[iter]=mD0[1];
      }
    } else {
      dd->InvMassD0bar(mD0bar);
      if(TMath::Abs(pdgdaughters[0])==321) {
       vars[iter]=mD0bar[0];
      }else{
       vars[iter]=mD0bar[1];
      }
   }
  }

  if(fVarsForOpt[1]){
    iter++;
    vars[iter]=dd->GetDCA();
  }

  if(fVarsForOpt[2]){
    iter++;
    vars[iter]=dd->GetDist12toPrim();
  }
  if(fVarsForOpt[3]){
    iter++;
    vars[iter]=dd->GetDist3toPrim();
  }
  if(fVarsForOpt[4]){
    iter++;
    vars[iter]=dd->GetDist4toPrim();
  }
  if(fVarsForOpt[5]){
    iter++;
    vars[iter]=dd->CosPointingAngle();
  }
  if(fVarsForOpt[6]){
    iter++;
    vars[iter]=dd->Pt();
  }
  if(fVarsForOpt[7]){
    iter++;
    vars[iter]=999999999.;
    printf("ERROR: optmization for rho mass cut not implemented\n");
  }
  if(fVarsForOpt[8]){
    iter++;
    vars[iter]=999999999.;
    printf("ERROR: optmization for PID cut not implemented\n");
  }
  
    if(cleanvtx)CleanOwnPrimaryVtx(dd,aod,origownvtx);
  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpipipi::IsSelected(TObject* obj,Int_t selectionLevel) {
  //
  // Apply selection
  //

  if(!fCutsRD){
    cout<<"Cut matrix not inizialized. Exit..."<<endl;
    return 0;
  }
  //PrintAll();
  AliAODRecoDecayHF4Prong* d=(AliAODRecoDecayHF4Prong*)obj;

  if(!d){
    cout<<"AliAODRecoDecayHF4Prong null"<<endl;
    return 0;
  }

  Double_t ptD=d->Pt();
  if(ptD<fMinPtCand) return 0;
  if(ptD>fMaxPtCand) return 0;

  if(d->HasBadDaughters()) return 0;  

  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(d)) return 0;
  }


  Int_t returnvalue=1;

  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {

    Int_t ptbin=PtBin(d->Pt());
  
    Int_t okD0=1,okD0bar=1;       
    Double_t mD0[2],mD0bar[2];
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    d->InvMassD0(mD0);
    if(TMath::Abs(mD0[0]-mD0PDG) > fCutsRD[GetGlobalIndex(0,ptbin)] &&
       TMath::Abs(mD0[1]-mD0PDG) > fCutsRD[GetGlobalIndex(0,ptbin)]) okD0 = 0;
    d->InvMassD0bar(mD0bar);
    if(TMath::Abs(mD0bar[0]-mD0PDG) > fCutsRD[GetGlobalIndex(0,ptbin)] &&
       TMath::Abs(mD0bar[1]-mD0PDG) > fCutsRD[GetGlobalIndex(0,ptbin)]) okD0bar = 0;
    if(!okD0 && !okD0bar) return 0;
    
    if(d->GetDCA() > fCutsRD[GetGlobalIndex(1,ptbin)]
       || d->GetDCA(3) > fCutsRD[GetGlobalIndex(1,ptbin)]
       || d->GetDCA(2) > fCutsRD[GetGlobalIndex(1,ptbin)]
       || d->GetDCA(5) > fCutsRD[GetGlobalIndex(1,ptbin)]) return 0;
    if(d->GetDist12toPrim() < fCutsRD[GetGlobalIndex(2,ptbin)]) return 0;
    if(d->GetDist3toPrim() < fCutsRD[GetGlobalIndex(3,ptbin)]) return 0;
    if(d->GetDist4toPrim() < fCutsRD[GetGlobalIndex(4,ptbin)]) return 0;
    if(d->CosPointingAngle() < fCutsRD[GetGlobalIndex(5,ptbin)]) return 0;
    if(d->Pt() < fCutsRD[GetGlobalIndex(6,ptbin)]) return 0;
    if(!d->CutRhoMass(mD0,mD0bar,fCutsRD[GetGlobalIndex(0,ptbin)],fCutsRD[GetGlobalIndex(7,ptbin)])) return 0;

    if (okD0) returnvalue=1; //cuts passed as D0
    if (okD0bar) returnvalue=2; //cuts passed as D0bar
    if (okD0 && okD0bar) returnvalue=3; //cuts passed as D0 and D0bar
  }

  // selection on PID (from AliAODPidHF)
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kPID) {

    Int_t selD01 = D01Selected(d,AliRDHFCuts::kCandidate);
    Int_t selD02 = D02Selected(d,AliRDHFCuts::kCandidate);  
    Int_t selD0bar1 = D0bar1Selected(d,AliRDHFCuts::kCandidate);
    Int_t selD0bar2 = D0bar2Selected(d,AliRDHFCuts::kCandidate);

    Int_t d01PID = 0, d02PID = 0, d0bar1PID = 0, d0bar2PID = 0;
    returnvalue = IsSelectedFromPID(d, &d01PID, &d02PID, &d0bar1PID, &d0bar2PID);  //This returnvalue is dummy! Now it's modified as it must be!

returnvalue = 0;

    if((selD01 == 1 && d01PID == 1)||(selD02 == 1 && d02PID == 1)||(selD0bar1 == 1 && d0bar1PID == 1)||(selD0bar2 == 1 && d0bar2PID == 1)) returnvalue = 1;
  }

  return returnvalue;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpipipi::IsSelectedFromPID(AliAODRecoDecayHF4Prong *d, Int_t *hyp1, Int_t *hyp2, Int_t *hyp3, Int_t *hyp4) {
  //
  // Apply selection (using AliAODPidHF methods)
  // Mass hypothesis true if each particle is at least compatible with specie of hypothesis
  // 

  Int_t output=0;

  Int_t matchK[4], matchPi[4];
  Double_t ptlimit[2] = {0.6,0.8};
  AliAODTrack* trk[4];
  trk[0] = (AliAODTrack*)d->GetDaughter(0);
  trk[1] = (AliAODTrack*)d->GetDaughter(1);
  trk[2] = (AliAODTrack*)d->GetDaughter(2);
  trk[3] = (AliAODTrack*)d->GetDaughter(3);
//  if(!trk[0] || !trk[1] || !trk[2] || !trk[3]) {          //REMOVED (needs also AliAODEvent to be passed, here and in IsSelected method)
//    trk[0]=aodIn->GetTrack(trkIDtoEntry[d->GetProngID(0)]);
//    trk[1]=aodIn->GetTrack(trkIDtoEntry[d->GetProngID(1)]);
//    trk[2]=aodIn->GetTrack(trkIDtoEntry[d->GetProngID(2)]);
//    trk[3]=aodIn->GetTrack(trkIDtoEntry[d->GetProngID(3)]);}

  AliAODPidHF* PidObj = new AliAODPidHF();
  PidObj->SetAsym(kTRUE);
  PidObj->SetPLimit(ptlimit);
  PidObj->SetSigma(0,2.);  //TPC sigma, in three pT ranges
  PidObj->SetSigma(1,1.);
  PidObj->SetSigma(2,0.);  
  PidObj->SetSigma(3,2.);  //TOF sigma, whole pT range
  PidObj->SetTPC(kTRUE);
  PidObj->SetTOF(kTRUE);
   
  for(Int_t ii=0; ii<4; ii++) {
    PidObj->SetSigma(0,2.);
    matchK[ii] = PidObj->MatchTPCTOF(trk[ii],1,3,kTRUE); //arguments: track, mode, particle#, compatibility allowed
    PidObj->SetSigma(0,2.);
    matchPi[ii] = PidObj->MatchTPCTOF(trk[ii],1,2,kTRUE); 
  }

  //Rho invariant mass under various hypotheses (to be matched with PID infos in order to accet the candidate)
  Int_t d01rho03 = 0, d01rho23 = 0, d02rho01 = 0, d02rho12 = 0, d0bar1rho12 = 0, d0bar1rho23 = 0, d0bar2rho01 = 0, d0bar2rho03 = 0;
  if(TMath::Abs(0.775 - d->InvMassRho(0,3))<0.1)  {d01rho03 = 1; d0bar2rho03 = 1;}
  if(TMath::Abs(0.775 - d->InvMassRho(2,3))<0.1)  {d01rho23 = 1; d0bar1rho23 = 1;}
  if(TMath::Abs(0.775 - d->InvMassRho(0,1))<0.1)  {d02rho01 = 1; d0bar2rho01 = 1;}
  if(TMath::Abs(0.775 - d->InvMassRho(1,2))<0.1)  {d02rho12 = 1; d0bar1rho12 = 1;}
  Int_t d01rho = 0, d02rho = 0, d0bar1rho = 0, d0bar2rho = 0;
  if(d01rho03==1||d01rho23==1) d01rho = 1;
  if(d02rho01==1||d02rho12==1) d02rho = 1;
  if(d0bar1rho12==1||d0bar1rho23==1) d0bar1rho = 1;
  if(d0bar2rho01==1||d0bar2rho03==1) d0bar2rho = 1;

  //This way because there could be multiple hypotheses accepted
  if(d01rho==1 && (matchK[1]>=0 && matchPi[0]>=0 && matchPi[2]>=0 && matchPi[3]>=0)) {*hyp1 = 1; output = 1;} //d01 hyp
  if(d02rho==1 && (matchK[3]>=0 && matchPi[0]>=0 && matchPi[1]>=0 && matchPi[2]>=0)) {*hyp2 = 1; output = 1;} //d02 hyp
  if(d0bar1rho==1 && (matchK[0]>=0 && matchPi[1]>=0 && matchPi[2]>=0 && matchPi[3]>=0)) {*hyp3 = 1; output = 1;} //d0bar1 hyp
  if(d0bar2rho==1 && (matchK[2]>=0 && matchPi[0]>=0 && matchPi[1]>=0 && matchPi[3]>=0)) {*hyp4 = 1; output = 1;} //d0bar2 hyp

  return output;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpipipi::D01Selected(TObject* obj,Int_t selectionLevel) {
  //
  // Apply selection
  //

  if(!fCutsRD){
    cout<<"Cut matrix not inizialized. Exit..."<<endl;
    return 0;
  }
  //PrintAll();
  AliAODRecoDecayHF4Prong* d=(AliAODRecoDecayHF4Prong*)obj;

  if(!d){
    cout<<"AliAODRecoDecayHF4Prong null"<<endl;
    return 0;
  }

  Int_t returnvalue=0;

  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {

    Int_t ptbin=PtBin(d->Pt());
    
    Double_t mD0[2];
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    d->InvMassD0(mD0);
    if(TMath::Abs(mD0[0]-mD0PDG) < fCutsRD[GetGlobalIndex(0,ptbin)]) returnvalue = 1;
  }

  return returnvalue;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpipipi::D02Selected(TObject* obj,Int_t selectionLevel) {
  //
  // Apply selection
  //

  if(!fCutsRD){
    cout<<"Cut matrix not inizialized. Exit..."<<endl;
    return 0;
  }
  //PrintAll();
  AliAODRecoDecayHF4Prong* d=(AliAODRecoDecayHF4Prong*)obj;

  if(!d){
    cout<<"AliAODRecoDecayHF4Prong null"<<endl;
    return 0;
  }

  Int_t returnvalue=0;

  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {

    Int_t ptbin=PtBin(d->Pt());
      
    Double_t mD0[2];
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    d->InvMassD0(mD0);
    if(TMath::Abs(mD0[1]-mD0PDG) < fCutsRD[GetGlobalIndex(0,ptbin)]) returnvalue = 1;
  }

  return returnvalue;
}//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpipipi::D0bar1Selected(TObject* obj,Int_t selectionLevel) {
  //
  // Apply selection
  //

  if(!fCutsRD){
    cout<<"Cut matrix not inizialized. Exit..."<<endl;
    return 0;
  }
  //PrintAll();
  AliAODRecoDecayHF4Prong* d=(AliAODRecoDecayHF4Prong*)obj;

  if(!d){
    cout<<"AliAODRecoDecayHF4Prong null"<<endl;
    return 0;
  }

  Int_t returnvalue=0;

  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {

    Int_t ptbin=PtBin(d->Pt());
 
    Double_t mD0bar[2];
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    d->InvMassD0bar(mD0bar);
    if(TMath::Abs(mD0bar[0]-mD0PDG) < fCutsRD[GetGlobalIndex(0,ptbin)]) returnvalue = 1;
  }

  return returnvalue;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpipipi::D0bar2Selected(TObject* obj,Int_t selectionLevel) {
  //
  // Apply selection
  //

  if(!fCutsRD){
    cout<<"Cut matrix not inizialized. Exit..."<<endl;
    return 0;
  }
  //PrintAll();
  AliAODRecoDecayHF4Prong* d=(AliAODRecoDecayHF4Prong*)obj;

  if(!d){
    cout<<"AliAODRecoDecayHF4Prong null"<<endl;
    return 0;
  }

  Int_t returnvalue=0;

  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {

    Int_t ptbin=PtBin(d->Pt());
    
    Double_t mD0bar[2];
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    d->InvMassD0bar(mD0bar);
    if(TMath::Abs(mD0bar[1]-mD0PDG) < fCutsRD[GetGlobalIndex(0,ptbin)]) returnvalue = 1;
  }

  return returnvalue;
}
//---------------------------------------------------------------------------
Bool_t AliRDHFCutsD0toKpipipi::IsInFiducialAcceptance(Double_t pt, Double_t y) const
{
  //
  // Checking if D0 is in fiducial acceptance region 
  //

  if(pt > 5.) {
    // applying cut for pt > 5 GeV
    AliDebug(4,Form("pt of D0 = %f (> 5), cutting at |y| < 0.8\n",pt)); 
    if (TMath::Abs(y) > 0.8){
      return kFALSE;
    }
  } else {
    // appliying smooth cut for pt < 5 GeV
    Double_t maxFiducialY = -0.2/15*pt*pt+1.9/15*pt+0.5; 
    Double_t minFiducialY = 0.2/15*pt*pt-1.9/15*pt-0.5;		
    AliDebug(4,Form("pt of D0 = %f (< 5), cutting  according to the fiducial zone [%f, %f]\n",pt,minFiducialY,maxFiducialY)); 
    if (y < minFiducialY || y > maxFiducialY){
      return kFALSE;
    }
  }

  return kTRUE;
}

