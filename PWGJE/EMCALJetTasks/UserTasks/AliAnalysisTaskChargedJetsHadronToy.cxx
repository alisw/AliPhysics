/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: R. Haake.                                                      *
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
#include <iostream>
#include <TRandom3.h>
#include <AliLog.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <AliEmcalJet.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include <AliAODTrack.h>
#include <TClonesArray.h>
#include "AliAnalysisTaskChargedJetsHadronToy.h"

ClassImp(AliAnalysisTaskChargedJetsHadronToy)

//_____________________________________________________________________________________________________
AliAnalysisTaskChargedJetsHadronToy::AliAnalysisTaskChargedJetsHadronToy() :
  AliAnalysisTaskSE("AliAnalysisTaskChargedJetsHadronToy"), fCreateUE(1), fCreateJets(1), fUEMultDistribution(0), fUEDistribution(0), fUEMultiplicity(1000), fGeneratedJetParticleDistribution(0), fGeneratedJetPtDistribution(0), fGeneratedJetCount(1), fGeneratedJetPtMin(20.), fGeneratedJetPtMax(30.), fGeneratedJetWidthPhi(0.2), fGeneratedJetWidthEta(0.2), fGeneratedJetMinEta(-0.9), fGeneratedJetMaxEta(0.9), fInputArrTracks(0), fInputArrTracksName(""), fOutputArrTracks(0), fOutputArrTracksName(""), fGeneratedJetsArr(0), fGeneratedJetsArrName(""), fDistEtaGaussian(0), fDistPhiGaussian(0), fRandom(), fInitialized()
{
// constructor
}

//_____________________________________________________________________________________________________
AliAnalysisTaskChargedJetsHadronToy::~AliAnalysisTaskChargedJetsHadronToy()
{
// destructor
  if(fUEMultDistribution)
    delete fUEMultDistribution;
  if(fUEDistribution)
    delete fUEDistribution;
  if(fGeneratedJetPtDistribution)
    delete fGeneratedJetPtDistribution;
  if(fGeneratedJetParticleDistribution)
    delete fGeneratedJetParticleDistribution;
  if(fDistPhiGaussian)
    delete fDistPhiGaussian;
  if(fDistEtaGaussian)
    delete fDistEtaGaussian;
}


//_____________________________________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronToy::UserCreateOutputObjects()
{
  fRandom = new TRandom3(0);
}

//_____________________________________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronToy::ExecOnce()
{
  // Check if input array can be loaded
  if(!fInputArrTracksName.IsNull())
  {
    fInputArrTracks = static_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s", fInputArrTracksName.Data())));
    if(!fInputArrTracks)
      AliFatal(Form("Input array '%s' demanded, but not found!", fInputArrTracksName.Data()));
    if(strcmp(fInputArrTracks->GetClass()->GetName(), "AliAODTrack"))
      AliFatal(Form("Input array has track type %s. Only AliAODTracks are supported.", fInputArrTracks->GetClass()->GetName()));
  }

  // Check if output arrays can be created
  if((InputEvent()->FindListObject(Form("%s", fOutputArrTracksName.Data()))))
    AliFatal(Form("Output array '%s' already exists in the event! Rename it.", fOutputArrTracksName.Data()));
  if((InputEvent()->FindListObject(Form("%s", fGeneratedJetsArrName.Data()))))
    AliFatal(Form("Output array '%s' already exists in the event! Rename it.", fGeneratedJetsArrName.Data()));

  // Define functions used in toy model + put the arrays to the event
  // NOTE: fOutputArrTracks must be added in any case since it contains the tracks of UE and jets
  fOutputArrTracks = new TClonesArray("AliAODTrack");
  fOutputArrTracks->SetName(fOutputArrTracksName.Data());
  InputEvent()->AddObject(fOutputArrTracks);

  if(fCreateUE)
  {
    if(!fUEDistribution)
    {
      std::cout << "\n### Distribution for the UE not given -- using default thermal distribution ###\n\n";
      fUEDistribution = new TF1("fUEDistribution","[0]*exp([1]*x)",0,200.);
      fUEDistribution->SetParameters(1.0,-1.5);
    }
    fUEDistribution->SetNpx(400);
  }

  if(fCreateJets)
  {
    fGeneratedJetsArr = new TClonesArray("AliEmcalJet");
    fGeneratedJetsArr->SetName(fGeneratedJetsArrName.Data());
    InputEvent()->AddObject(fGeneratedJetsArr);
    if(!fGeneratedJetParticleDistribution)
    {
      std::cout << "\n### Distribution for the particles in jets not given -- using default power law distribution ###\n\n";
      fGeneratedJetParticleDistribution = new TF1("fGeneratedJetParticleDistribution","[0]*(x^[1])",4.,120.);
      fGeneratedJetParticleDistribution->SetParameters(6.0,-2.5);
    }
    fGeneratedJetParticleDistribution->SetNpx(400);

    if(!fGeneratedJetPtDistribution)
    {
      std::cout << "\n### Distribution for the jet pt not given -- only min jet pT ###\n\n";
    }
    else
    {
      Int_t minBin = fGeneratedJetPtDistribution->GetXaxis()->FindBin(fGeneratedJetPtMin);
      Int_t maxBin = fGeneratedJetPtDistribution->GetXaxis()->FindBin(fGeneratedJetPtMax);
      for(Int_t i=0; i<=fGeneratedJetPtDistribution->GetNbinsX(); i++)
        if(i < minBin || i > maxBin)
        {
          fGeneratedJetPtDistribution->SetBinContent(i,0.0);
          fGeneratedJetPtDistribution->SetBinError(i,0.0);
        }

    }

    fDistEtaGaussian = new TF1("fDistEtaGaussian","gaus(0)",-1.,1.);// gaus(0) is  [0]*exp(-0.5*((x-[1])/[2])**2)
    fDistEtaGaussian->SetParameters(1.0,0.,0.5*fGeneratedJetWidthEta);
    fDistPhiGaussian = new TF1("fDistPhiGaussian","gaus(0)",-1.,1.);
    fDistPhiGaussian->SetParameters(1.0,0.,0.5*fGeneratedJetWidthPhi);
  }

  fInitialized = kTRUE;
}

//_____________________________________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronToy::UserExec(Option_t *)
{
  // Run once the exec once (must be here to have the input event ready)
  if(!fInitialized)
    ExecOnce();

  AssembleEvent();
//    std::cout << fOutputArrTracks->GetName() << ":" << fOutputArrTracks->GetEntries() << std::endl;

}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronToy::AssembleEvent()
{
  // Create the event from the several inputs and run the jet finder
  // Note: those tracks are guaranteed to be AliAODTracks
  if(fInputArrTracks)
    fOutputArrTracks->AddAll(fInputArrTracks);

  // 1. Create a vertex if there is none (needed by correlation task)
  if(dynamic_cast<AliESDEvent*>(InputEvent()))
  {
    if(!(dynamic_cast<AliESDEvent*>(InputEvent()))->GetPrimaryVertexTracks()->GetNContributors())
      static_cast<AliESDEvent*>(fInputEvent)->SetPrimaryVertexTracks(new AliESDVertex(0.,0., 100));
  }
  else if(dynamic_cast<AliAODEvent*>(InputEvent()))
  {
    if( (!(dynamic_cast<AliAODEvent*>(InputEvent()))->GetPrimaryVertex()) || (!(dynamic_cast<AliAODEvent*>(InputEvent()))->GetPrimaryVertex()->GetNContributors()) )
    {
      Double_t* p = new Double_t[3] {0., 0., 0.};
      AliAODVertex* vertex = new AliAODVertex(p,1.);
      vertex->SetNContributors(100);
      vertex->SetName("PrimaryVertex");
      static_cast<AliAODEvent*>(fInputEvent)->AddVertex(vertex);
    }
  }

  // 2. Create underlying event
  if(fCreateUE)
  {
    Double_t etaMin = -0.9;
    Double_t etaMax = +0.9;

    Double_t thrownPt = 0.;
    Int_t count = fOutputArrTracks->GetEntries();

    Int_t multiplicity = fUEMultiplicity;
    if(fUEMultDistribution)
      multiplicity = (Int_t)fUEMultDistribution->GetRandom();

    for(Int_t i=0;i<multiplicity; i++)
    {
      Double_t trackPt = fUEDistribution->GetRandom();
      Double_t trackEta = etaMin + fRandom->Rndm()*(etaMax-etaMin);
      Double_t trackTheta = 2.*atan(exp(-trackEta));
      Double_t trackPhi = fRandom->Rndm()*TMath::TwoPi();
      Double_t trackCharge = fRandom->Rndm() - 0.5;

      if(trackCharge>0) trackCharge = 1; else trackCharge = -1;

      // Add very basic particle to event
      new ((*fOutputArrTracks)[count]) AliAODTrack();
      static_cast<AliAODTrack*>(fOutputArrTracks->At(count))->SetPt(trackPt);
      static_cast<AliAODTrack*>(fOutputArrTracks->At(count))->SetPhi(trackPhi);
      static_cast<AliAODTrack*>(fOutputArrTracks->At(count))->SetTheta(trackTheta); // AliAODTrack cannot set eta directly
      static_cast<AliAODTrack*>(fOutputArrTracks->At(count))->SetCharge(trackCharge);
      count++;
      thrownPt += trackPt;
    }

    std::cout << Form("Created the underlying event using %i particles. Pt per unit area: %5.2f", multiplicity, thrownPt/(TMath::TwoPi()*(etaMax-etaMin))) << std::endl;
  }

  // 3. Embed gaussian jets into event
  if(fCreateJets)
  {
    // Define jets and throw them into the acceptance

    Double_t thrownPt = 0.;
    for(Int_t i=0;i<fGeneratedJetCount; i++)
    {
      Double_t jetEta = fGeneratedJetMinEta + fRandom->Rndm()*(fGeneratedJetMaxEta-fGeneratedJetMinEta);
      Double_t jetPhi = fRandom->Rndm()*TMath::TwoPi();

      Double_t jetPt = fGeneratedJetPtMin;
      if(fGeneratedJetPtDistribution)
        jetPt = fGeneratedJetPtDistribution->GetRandom();

      Int_t count = fOutputArrTracks->GetEntries();
      Int_t particlesInJet = 0;
      while(thrownPt < jetPt)
      {
        Double_t trackPt = fGeneratedJetParticleDistribution->GetRandom();
        Double_t trackEta = jetEta + fDistEtaGaussian->GetRandom();
        Double_t trackTheta = 2.*atan(exp(-trackEta));
        Double_t trackPhi = jetPhi + fDistPhiGaussian->GetRandom();
        Double_t trackCharge = fRandom->Rndm() - 0.5;
        if(trackCharge>0) trackCharge = 1; else trackCharge = -1;

        // Add particle to event
        new ((*fOutputArrTracks)[count]) AliAODTrack();
        static_cast<AliAODTrack*>(fOutputArrTracks->At(count))->SetPt(trackPt);
        static_cast<AliAODTrack*>(fOutputArrTracks->At(count))->SetPhi(trackPhi);
        static_cast<AliAODTrack*>(fOutputArrTracks->At(count))->SetTheta(trackTheta); // AliAODTrack cannot set eta directly
        static_cast<AliAODTrack*>(fOutputArrTracks->At(count))->SetCharge(trackCharge);

        count++;
        particlesInJet++;
        thrownPt += trackPt;
      }
      std::cout << Form("Created a gaussian %5.2f GeV jet using %i particles",thrownPt, particlesInJet) << std::endl;

      // Save the generated jet for later matching
      new ((*fGeneratedJetsArr)[i]) AliEmcalJet(thrownPt, jetEta, jetPhi, 0);
    }
  }
}

