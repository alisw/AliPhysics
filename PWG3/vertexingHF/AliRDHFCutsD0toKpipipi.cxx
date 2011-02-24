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
// Class for cuts on AOD reconstructed D0->Kpipipi
//
// Author: r.romita@gsi.de, andrea.dainese@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>

#include "AliRDHFCutsD0toKpipipi.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"

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
void AliRDHFCutsD0toKpipipi::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) {
  // 
  // Fills in vars the values of the variables 
  //

  if(nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsD0toKpipipi::GetCutsVarsForOpt: wrong number of variables\n");
    return;
  }

  AliAODRecoDecayHF4Prong *dd = (AliAODRecoDecayHF4Prong*)d;

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
    
    if(d->GetDCA() > fCutsRD[GetGlobalIndex(1,ptbin)]) return 0;
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

  return returnvalue;
}
//---------------------------------------------------------------------------
