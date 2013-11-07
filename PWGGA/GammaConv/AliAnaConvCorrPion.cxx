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



#include "TH2D.h"
#include "AliAnaConvCorrPion.h"
//#include "AliAODTrack.h"
#include "AliAODConversionParticle.h"




using namespace std;
ClassImp(AliAnaConvCorrPion)

//________________________________________________________________________________
AliAnaConvCorrPion::AliAnaConvCorrPion() :
AliAnaConvCorrBase("pion_hadron_corr", "Pion dPhi"),
//hTriggerPtvsMass(NULL),
  hTriggerPtvsMass(NULL),
  fAxisM()  
{
  //consctructor
  InitMassAxis();
}
//________________________________________________________________________________
AliAnaConvCorrPion::AliAnaConvCorrPion(TString name, TString title = "Pion Corr") :
  AliAnaConvCorrBase(name, title),
  //hTriggerPtvsMass(NULL),
  hTriggerPtvsMass(NULL),
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
  Double_t mbins[7] = {0.1, 0.11, 0.12, 0.15, 0.16, 0.18, 0.2};
  fAxisM.Set(6, mbins);
  fAxisM.SetNameTitle("InvMass", "invariant mass");
  GetAxisList().AddAt(&fAxisM, 6);
  GetTrackAxisList().AddAt(&fAxisM, 5);
  GetTrigAxisList().AddAt(&fAxisM, 4);
}

///________________________________________________________________________________
void AliAnaConvCorrPion::CreateHistograms() {
  //Create histograms
  CreateBaseHistograms();
 
  hTriggerPtvsMass = new TH2D(Form("hTriggerPtvsMass_all_%s", GetName()), "Pt vs Mass all pizero", 
			      400, 0, .400, GetAxistPt().GetNbins(), GetAxistPt().GetXbins()->GetArray());
  GetHistograms()->Add(hTriggerPtvsMass);
}


///________________________________________________________________________________
void AliAnaConvCorrPion::FillTriggerCounters(const AliAODConversionParticle * particle) {
  hTriggerPtvsMass->Fill(particle->M(), particle->Pt());
}

