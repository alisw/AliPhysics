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
// Class for cuts on AOD reconstructed D0->Kpi
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>

#include "AliRDHFCutsD0toKpipipi.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"

ClassImp(AliRDHFCutsD0toKpipipi)

//--------------------------------------------------------------------------
AliRDHFCutsD0toKpipipi::AliRDHFCutsD0toKpipipi() : 
AliRDHFCuts()
{
  //
  // Default Constructor
  //
  Int_t nvars=9;
  SetNVars(nvars);
  TString varNames[9]={"inv. mass [GeV]",   
		       "dca [cm]",
		       "cosThetaStar", 
		       "pTK [GeV/c]",
		       "pTPi [GeV/c]",
		       "d0K [cm]",
		       "d0Pi [cm]",
		       "d0d0 [cm^2]",
		       "cosThetaPoint"};
  Bool_t isUpperCut[9]={kTRUE,
			kTRUE,
			kTRUE,
			kFALSE,
			kFALSE,
			kTRUE,
			kTRUE,
			kTRUE,
			kFALSE};
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[9]={kFALSE,
		    kTRUE,
		    kTRUE,
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kFALSE,
		    kTRUE,
		    kTRUE};
  SetVarsForOpt(4,forOpt);
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
 
  
  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsD0toKpipipi::IsSelected(TObject* obj,Int_t selectionLevel) {
  //
  // Apply selection
  //

  if(!fCutsRD){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  //PrintAll();
  AliAODRecoDecayHF4Prong* d=(AliAODRecoDecayHF4Prong*)obj;

  if(!d){
    cout<<"AliAODRecoDecayHF4Prong null"<<endl;
    return 0;
  }


  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(d)) return 0;
  }


  Int_t returnvalue=1;

  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    
  }


  return returnvalue;

}
//---------------------------------------------------------------------------
