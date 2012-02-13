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
#include "TClonesArray.h"
#include "AliAODConversionParticle.h"
//#include "AliAODConversionMother.h"
//#include "AliAODConversionPhoton.h"
//#include "THnSparse.h"
//#include "TH2F.h"

#include <iostream>


using namespace std;
ClassImp(AliAnaConvCorrPion)

//________________________________________________________________________________
AliAnaConvCorrPion::AliAnaConvCorrPion() :
AliAnaConvCorrBase("pion_hadron_corr", "Pion dPhi"),
  hTriggerPtvsMass(NULL),
  fAxisM()
{
  //consctructor
  InitMassAxis();
}
//________________________________________________________________________________
AliAnaConvCorrPion::AliAnaConvCorrPion(TString name, TString title = "Pion Corr") :
  AliAnaConvCorrBase(name, title),
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
  hTriggerPtvsMass = new TH2D("hTriggerPtvsMass", "Pt vs Mass", 400, 0, .400, GetAxistPt().GetNbins(), GetAxistPt().GetXbins()->GetArray());
  GetHistograms()->Add(hTriggerPtvsMass);
}


///________________________________________________________________________________
void AliAnaConvCorrPion::FillTriggerCounters(const AliAODConversionParticle * particle, Bool_t isolated) {
  //Fill histograms counting triggers
  fHNTriggers[isolated]->Fill(particle->Pt());
  hTriggerPtvsMass->Fill(particle->M(), particle->Pt());
}

//________________________________________________________________________________
// void AliAnaConvCorrPion::Process(TClonesArray * pions, TClonesArray * photons, TClonesArray * tracks) {

//   for(Int_t ip = 0; ip < pions->GetEntriesFast(); ip++) {

// 	AliAODConversionParticle * pion = static_cast<AliAODConversionParticle*>(pions->UncheckedAt(ip));
	
// 	Int_t tIDs[4] = {-1, -1, -1, -1};
// 	AliAODConversionParticle * photon1 = static_cast<AliAODConversionParticle*>(photons->UncheckedAt(pion->GetLabel(0)));
// 	tIDs[0] =  photon1->GetLabel(0);
// 	tIDs[1] =  photon1->GetLabel(1);
// 	AliAODConversionParticle * photon2 = static_cast<AliAODConversionParticle*>(photons->UncheckedAt(pion->GetLabel(1)));
// 	tIDs[2] =  photon2->GetLabel(0);
// 	tIDs[3] =  photon2->GetLabel(1);
	
// 	CorrelateWithTracks(static_cast<AliAODConversionParticle*>(pion), tracks, tIDs, kFALSE);
//   }
// }

