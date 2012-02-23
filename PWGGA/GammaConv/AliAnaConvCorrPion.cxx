/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Author: Svein Lindal <slindal@fys.uio.no>                      *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliAnaConvCorrPion.cxx
/// @author Svein Lindal
/// @brief  Class used to run conversion gamma/pion - hadron/jet analysis



#include "TH2F.h"
#include "AliAnaConvCorrPion.h"
#include "AliAODConversionParticle.h"
#include <iostream>


using namespace std;
ClassImp(AliAnaConvCorrPion)

//________________________________________________________________________________
AliAnaConvCorrPion::AliAnaConvCorrPion() :
AliAnaConvCorrBase("pion_hadron_corr", "Pion dPhi"),
  fAxisM()
{
  //consctructor
  InitMassAxis();
}
//________________________________________________________________________________
AliAnaConvCorrPion::AliAnaConvCorrPion(TString name, TString title = "Pion Corr") :
  AliAnaConvCorrBase(name, title),
  fAxisM()
{
  //consctructor
  InitMassAxis();
}


//________________________________________________________________________________
AliAnaConvCorrPion::~AliAnaConvCorrPion() {
  //destructor
}

//________________________________________________________________________________
void AliAnaConvCorrPion::InitMassAxis() {
  //Init mass axis
  
  Double_t mbins[7] = {0.1, 0.11, 0.12, 0.15, 0.16, 0.18, 0.2};
  fAxisM.Set(6, mbins);
  fAxisM.SetNameTitle("InvMass", "invariant mass");
  GetAxisList().AddAt(&fAxisM, 5);
  GetTrackAxisList().AddAt(&fAxisM, 5);
  GetTrigAxisList().AddAt(&fAxisM, 4);
}

///________________________________________________________________________________
void AliAnaConvCorrPion::CreateHistograms() {
  //Create histograms
  CreateBaseHistograms();

  for(Int_t i = 0; i < 3; i++) {
  	TH2F * hpthist = new TH2F(Form("%s_iso_%d_ptVsmass", fName.Data(), i), 
  							  Form("%s_iso_%d_ptVsmass", fName.Data(), i), 
  							  400, 0, .400,
  							  GetAxistPt().GetNbins(), GetAxistPt().GetXbins()->GetArray());
  	GetHistograms()->Add(hpthist);
  }
}


///________________________________________________________________________________
void AliAnaConvCorrPion::FillTriggerCounters(const AliAODConversionParticle * particle, Int_t isolated) {
  //Fill histograms counting triggers
  TH2F *hpt = dynamic_cast<TH2F*>(GetHistograms()->FindObject(Form("%s_iso_%d_ptVsmass", fName.Data(), isolated)));
  if(hpt) hpt->Fill(particle->M(), particle->Pt());
}
