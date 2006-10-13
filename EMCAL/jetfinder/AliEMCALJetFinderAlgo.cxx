/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//_________________________________________________________________________
//  Base Class for JetFinder Algorithms     
// --                  
//*-- Author: Mark Horner (LBL/UCT)
// --



#include "AliEMCALJetFinderAlgo.h"
#include "AliRun.h"
#include "AliEMCALGeometry.h"
#include "AliMagF.h"
ClassImp(AliEMCALJetFinderAlgo)

AliEMCALJetFinderAlgo::AliEMCALJetFinderAlgo() :
  fInputPointer(0),fOutputPointer(0),fOutputAllocated(kFALSE),
  fDebug(0),fPythiaComparison(0)
{
   fDebug =0;
   fOutputPointer=0;
   fInputPointer=0;
   fOutputAllocated=kFALSE;
}


AliEMCALJetFinderAlgo::AliEMCALJetFinderAlgo(const AliEMCALJetFinderAlgo& jfa) 
  : TTask(jfa.GetName(), jfa.GetTitle()),
    fInputPointer(jfa.fInputPointer),fOutputPointer(jfa.fOutputPointer),
    fOutputAllocated(jfa.fOutputAllocated), fDebug(jfa.fDebug),
    fPythiaComparison(jfa.fPythiaComparison)
{
  //copy ctor
}

AliEMCALJetFinderAlgo::~AliEMCALJetFinderAlgo()
{
  if (fOutputAllocated)
    delete fOutputPointer;
}

void AliEMCALJetFinderAlgo::SetOutput(AliEMCALJetFinderOutput* output) {
  if (fOutputAllocated)
    delete fOutputPointer;
  fOutputPointer=output;
  fOutputAllocated=kFALSE;
}

void AliEMCALJetFinderAlgo::InitInput(AliEMCALJetFinderInput* input)
{
// Take input data	
if (fDebug>1) Info("InitInput","Beginning InitInput");		
	fInputPointer = input;
	if (fOutputPointer==0) {
	  
	  if (fDebug>1) Info("InitInput","Allocating output object");		
	  fOutputPointer=new AliEMCALJetFinderOutput();
	  fOutputAllocated=kTRUE;
	}
	fOutputPointer->Reset(kResetAll);	
	// automatically copy parton and particle info to output object
	 
	for (Int_t counter = 0 ; counter < fInputPointer->GetNPartons();counter++)
	{
		fOutputPointer->AddParton(fInputPointer->GetParton(counter));
	}
	for (Int_t counter = 0 ; counter < fInputPointer->GetNParticles();counter++)
	{
		fOutputPointer->AddParticle(fInputPointer->GetParticle(counter));
	}
}

Float_t AliEMCALJetFinderAlgo::PropagatePhi(Float_t pt, Float_t charge, Bool_t& curls)
{
	// Propagates phi angle to EMCAL radius
	// //
 Float_t b = 0.0, rEMCAL = -1.0;
 if(rEMCAL<0) 
 {	
	 b =  gAlice->Field()->SolenoidField();
	 rEMCAL = AliEMCALGeometry::GetInstance()->GetIPDistance();
 }
 Float_t dPhi = 0.;
 Float_t rB = 3335.6 * pt / b;  // [cm]  (case of |charge|=1)
 if (2.*rB < rEMCAL) 
 {
	 curls = kTRUE;
	 return dPhi;
 }
 Float_t phi = TMath::ACos(1.-rEMCAL*rEMCAL/(2.*rB*rB));
 dPhi = TMath::ATan2(1.-TMath::Cos(phi), TMath::Sin(phi));
 dPhi = -TMath::Sign(dPhi, charge);
 return dPhi;
}

