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

//==========================================================================
// AliEPDependentDiHadronOnTheFlyMCTask - 
//    This analysis task constructs an event plane based on the 
//    VZERO detector range and consecutively fills several dihadron 
//    histograms with the dihadron-pair triggers having a certain 
//    angle to that event plane.
//==========================================================================
// It is mend for use with a HIJING or AMPT event generator.
// It divides the event plane into three sectors (in-plane, mid-plane and 
// out-of-plane), which correspond to the P0, P1 and P2 output containers. 
// The P-1 output is reserved for the plane-inclusive results.
//==========================================================================

#include "AliEPDependentDiHadronOnTheFlyMCTask.h"


ClassImp(AliEPDependentDiHadronOnTheFlyMCTask)


AliEPDependentDiHadronOnTheFlyMCTask::AliEPDependentDiHadronOnTheFlyMCTask():
AliAnalysisTaskSE("EPDependentDiHadronOnTheFlyMCTask"),
fListOutputHistograms(0x0),
fHistNoEvents(0x0),
fHistImpactVsMultiplicity(0x0),
fProfileFractionPrimaryTracks(0x0),
fHistReactionPlane(0x0),
fHistEventReactionPlane(0x0),
fHistTrueEventReactionPlane(0x0),
fHistos(0x0),
fHistosIn(0x0),
fHistosMid(0x0),
fHistosOut(0x0),
fEtaMinV0A(2.8),
fEtaMaxV0A(5.1),
fEtaMinV0C(-3.7),
fEtaMaxV0C(-1.7),
fTrackEtaCut(1e9),
fNoSectorsV0(8),
fCentralityEstimator("V0M")
{
  // Default constructor
  DefineOutput(1, TList::Class());
}


AliEPDependentDiHadronOnTheFlyMCTask::AliEPDependentDiHadronOnTheFlyMCTask(const char* name):
AliAnalysisTaskSE(name),
fListOutputHistograms(0x0),
fHistNoEvents(0x0),
fHistImpactVsMultiplicity(0x0),
fProfileFractionPrimaryTracks(0x0),
fHistReactionPlane(0x0),
fHistEventReactionPlane(0x0),
fHistTrueEventReactionPlane(0x0),
fHistos(0x0),
fHistosIn(0x0),
fHistosMid(0x0),
fHistosOut(0x0),
fEtaMinV0A(2.8),
fEtaMaxV0A(5.1),
fEtaMinV0C(-3.7),
fEtaMaxV0C(-1.7),
fTrackEtaCut(1e9),
fNoSectorsV0(8),
fCentralityEstimator("V0M")
{
  // Constructor with name
  DefineOutput(1, TList::Class());
}


AliEPDependentDiHadronOnTheFlyMCTask::~AliEPDependentDiHadronOnTheFlyMCTask()
{
  // Destructor
  // Nothing to destruct... yet
}


void AliEPDependentDiHadronOnTheFlyMCTask::UserCreateOutputObjects()
{
  // Initialize output list of containers
  if(fListOutputHistograms != NULL) {
    delete fListOutputHistograms;
    fListOutputHistograms = NULL;
  }
  if(!fListOutputHistograms){
    fListOutputHistograms = new TList();
    fListOutputHistograms->SetOwner();
  }

  // Create the output containers
  Double_t onebinlimits[2] = {0.0, 1.0};
  fHistNoEvents = new TH1D("noEventsVZEROOnTheFlyMCTask", "noEventsVZEROOnTheFlyMCTask", 1, onebinlimits);
  fListOutputHistograms->Add(fHistNoEvents);

  fHistImpactVsMultiplicity = new TH2D("impactVsLogMultiplicity", "Impact parameter versus 10-Log of Multiplicity", 40, 0, 20, 20, 0, 5);
  fListOutputHistograms->Add(fHistImpactVsMultiplicity);
  
  fProfileFractionPrimaryTracks = new TProfile("fractionPrimaryTracks", "fractionPrimaryTracks", 1, onebinlimits, "s");
  fListOutputHistograms->Add(fProfileFractionPrimaryTracks);

  fHistReactionPlane = new TH1D("ReactionPlane", "Reactionplane", 200, -2*TMath::Pi(), 2*TMath::Pi());
  fListOutputHistograms->Add(fHistReactionPlane);

  fHistEventReactionPlane = new TH1D("EventReactionPlane", "Two-event plane w.r.t. the two-reactionplane", 100, 0.0, TMath::Pi());
  fListOutputHistograms->Add(fHistEventReactionPlane);

  fHistTrueEventReactionPlane = new TH1D("TrueEventReactionPlane", "All particle two-event plane w.r.t. the two-reactionplane", 100, 0.0, TMath::Pi());
  fListOutputHistograms->Add(fHistTrueEventReactionPlane);

  fHistos = new AliUEHistograms("AliUEHistogramsSame-1", "4R", fCustomBinning);
  fHistos->SetTrackEtaCut(fTrackEtaCut);
  fListOutputHistograms->Add(fHistos);

  fHistosIn = new AliUEHistograms("AliUEHistogramsSame0", "4R", fCustomBinning);
  fHistosIn->SetTrackEtaCut(fTrackEtaCut);
  fListOutputHistograms->Add(fHistosIn);

  fHistosMid = new AliUEHistograms("AliUEHistogramsSame1", "4R", fCustomBinning);
  fHistosMid->SetTrackEtaCut(fTrackEtaCut);
  fListOutputHistograms->Add(fHistosMid);

  fHistosOut = new AliUEHistograms("AliUEHistogramsSame2", "4R", fCustomBinning);
  fHistosOut->SetTrackEtaCut(fTrackEtaCut);
  fListOutputHistograms->Add(fHistosOut);

  PostData(1, fListOutputHistograms);
}


void AliEPDependentDiHadronOnTheFlyMCTask::UserExec(Option_t* option)
{
  fHistNoEvents->Fill(0.5);

  AliMCEvent* mcEvent = dynamic_cast<AliMCEvent*>(MCEvent());
  
  Double_t reactionplane = 0, zvtx = 0, centrality = 0;
  GetEventDetails(mcEvent, reactionplane, zvtx, centrality);

  FillEventLevelQA(mcEvent);

  TObjArray* allTracks = SelectTracks(mcEvent);

  Double_t eventplaneV0 = 0;
  FillEventPlaneHistograms(allTracks, reactionplane, eventplaneV0);


//  TODO: bepaal centrality

  // Divide the leading tracks with respect to the event plane. Done in a single loop to save on computation time
  TObjArray* inplaneTracks = new TObjArray; // cutOnEventPlaneAngle(allTracks, 0.0, TMath::Pi()/6);
  TObjArray* midplaneTracks = new TObjArray; // cutOnEventPlaneAngle(allTracks, TMath::Pi()/6, TMath::Pi()/3);
  TObjArray* outplaneTracks = new TObjArray; // cutOnEventPlaneAngle(allTracks, TMath::Pi()/3, TMath::Pi()/2);
  for(Int_t part = 0; part < allTracks->GetEntriesFast(); part++) {
    AliMCParticle* track = dynamic_cast<AliMCParticle*>(allTracks->At(part));
    Double_t dphi = piToHalfPi(twoPiToPi(angleBetween(track->Phi(), reactionplane)));
    if(dphi < TMath::Pi()/6) inplaneTracks->Add(track);
    if((dphi >= TMath::Pi()/6) && (dphi < TMath::Pi()/3)) midplaneTracks->Add(track);
    if(dphi >= TMath::Pi()/3) outplaneTracks->Add(track);
  }

  fHistos->FillCorrelations(centrality, zvtx, AliUEHist::kCFStepAll, allTracks, allTracks, 1);
  fHistosIn->FillCorrelations(centrality, zvtx, AliUEHist::kCFStepAll, inplaneTracks, allTracks, 1);
  fHistosMid->FillCorrelations(centrality, zvtx, AliUEHist::kCFStepAll, midplaneTracks, allTracks, 1);
  fHistosOut->FillCorrelations(centrality, zvtx, AliUEHist::kCFStepAll, outplaneTracks, allTracks, 1);

  delete allTracks; delete inplaneTracks; delete midplaneTracks; delete outplaneTracks;
}


void AliEPDependentDiHadronOnTheFlyMCTask::FillEventLevelQA(AliMCEvent* mcEvent)
{
  AliGenHijingEventHeader* eventHeader = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());
  Double_t impactParameter = (Double_t) eventHeader->ImpactParameter();
  Double_t multiplicity = mcEvent->GetNumberOfTracks();

  if(multiplicity>0)
    fHistImpactVsMultiplicity->Fill(impactParameter, TMath::Log10(multiplicity));
}


TObjArray* AliEPDependentDiHadronOnTheFlyMCTask::SelectTracks(AliMCEvent* mcEvent)
{
  TObjArray* acceptedTracks = new TObjArray;
//  acceptedTracks->SetOwner(kTRUE);

  if(!mcEvent) {
    AliError("No correct mcEvent provided to SelectTracks");
    return 0x0;
  }

  if(mcEvent) {
    Int_t nPrimaryTracks = mcEvent->GetNumberOfPrimaries();
    Int_t nTracks = mcEvent->GetNumberOfTracks();
    fProfileFractionPrimaryTracks->Fill(0.5, nTracks/nPrimaryTracks);

    for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
      AliMCParticle* track = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(iTrack));
      if (!track) {
        AliError(Form("Could not retreive particle %d during selection", iTrack));
        continue;
      }
    
    // Allways exclude non-primary particles
      if( ! mcEvent->IsPhysicalPrimary(iTrack) ) continue;
    
    // TODO: See AliAnalysisTaskBFPsi::GetAcceptedTracks for other idea's about what to exclude

      if( track->Charge() == 0 ) continue;

      acceptedTracks->Add(track);
    }
  }

  return acceptedTracks;
}


void AliEPDependentDiHadronOnTheFlyMCTask::GetEventDetails(AliMCEvent* mcEvent, Double_t& reactionplane, Double_t& zvtx, Double_t& centrality)
{
  AliGenHijingEventHeader* eventHeader = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());

  reactionplane = (Double_t) eventHeader->ReactionPlaneAngle();

  TArrayF vertex = TArrayF(3);
  eventHeader->PrimaryVertex(vertex);
  zvtx = vertex.At(2);

/* TODO: terug als primertask werkt.
  AliMultSelection* multSelection = (AliMultSelection*) mcEvent->FindListObject("MultSelection");
  if (!multSelection)
    AliFatal("MultSelection not found in input event. Did the AliMultSelectionTask run before this task?");
*/

  centrality = 35.0; // TODO: terug als primertask werkt: multSelection->GetMultiplicityPercentile(fCentralityEstimator, kTRUE);
}


void AliEPDependentDiHadronOnTheFlyMCTask::FillEventPlaneHistograms(TObjArray* allTracks, Double_t reactionplane, Double_t& eventplaneV0)
{
  TH1D* V0A = FillVirtualV0(allTracks, fEtaMinV0A, fEtaMaxV0A);
  TH1D* V0C = FillVirtualV0(allTracks, fEtaMinV0C, fEtaMaxV0C);
  TH1D* V0 = new TH1D("V0dist", "SimulatedV0", fNoSectorsV0, 0, 2*TMath::Pi());
  V0->Add(V0A, V0C);
  Double_t trueeventplane = ComputeTrueEventPlane(allTracks, 2);
  Double_t eventplaneV0A = ComputeEventPlane(V0A, 2);
  Double_t eventplaneV0C = ComputeEventPlane(V0C, 2);
  eventplaneV0 = ComputeEventPlane(V0, 2);
  Double_t angleBetweenPlanes = twoPiToPi(angleBetween(reactionplane, eventplaneV0));
  Double_t angleBetweenPlanesTrue = twoPiToPi(angleBetween(trueeventplane, reactionplane));

  fHistReactionPlane->Fill(reactionplane);
  fHistEventReactionPlane->Fill(angleBetweenPlanes);
  fHistTrueEventReactionPlane->Fill(angleBetweenPlanesTrue);

  delete V0A; delete V0C; delete V0;
}


TH1D* AliEPDependentDiHadronOnTheFlyMCTask::FillVirtualV0(TObjArray* allTracks, Double_t etaMin, Double_t etaMax)
{
  TH1D* sectors = new TH1D("V0dist", "SimulatedV0", fNoSectorsV0, 0, 2*TMath::Pi());
  Int_t nTracks = allTracks->GetEntriesFast();
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    AliMCParticle* track = dynamic_cast<AliMCParticle*>(allTracks->At(iTrack));
    Double_t eta = track->Eta();
    if((eta > etaMin) && (eta < etaMax)) 
      sectors->Fill(track->Phi());
  }
  return sectors;
}


Double_t AliEPDependentDiHadronOnTheFlyMCTask::ComputeTrueEventPlane(TObjArray* allTracks, Int_t harmonic)
{
  Int_t nSelected = 0;
  Int_t nTracks = allTracks->GetEntriesFast();
  Double_t qx = 0, qy = 0;
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    AliMCParticle* track = dynamic_cast<AliMCParticle*>(allTracks->At(iTrack));
    nSelected++;
    qx += TMath::Cos(harmonic*track->Phi());
    qy += TMath::Sin(harmonic*track->Phi());
  }
  if (nSelected==0) return -999;
  return (TMath::ATan2(qy, qx)/harmonic);
}


Double_t AliEPDependentDiHadronOnTheFlyMCTask::ComputeEventPlane(TH1D* sectors, Int_t harmonic)
{
  Double_t qx = 0, qy = 0;
  Int_t nSelected = 0;
  for (Int_t sector = 1; sector <= fNoSectorsV0; sector++) {
    Double_t phi = sectors->GetBinCenter(sector);
    Double_t mult = sectors->GetBinContent(sector);
    qx += mult*TMath::Cos(harmonic*phi);
    qy += mult*TMath::Sin(harmonic*phi);
    nSelected += mult;
  }
  if (nSelected==0) return -999;
  return (TMath::ATan2(qy, qx)/harmonic);
}


Double_t AliEPDependentDiHadronOnTheFlyMCTask::angleBetween(Double_t phi1, Double_t phi2)
{
  Double_t phi = phi1-phi2;
  if(phi<0)
    phi += 2*TMath::Pi();
  if(phi>2*TMath::Pi())
    phi -= 2*TMath::Pi();
  return phi;
}


Double_t AliEPDependentDiHadronOnTheFlyMCTask::twoPiToPi(Double_t phi)
{
  if(phi>=TMath::Pi())
    return phi-TMath::Pi();
  else
    return phi;
}

Double_t AliEPDependentDiHadronOnTheFlyMCTask::piToHalfPi(Double_t phi)
{
  if(phi>TMath::Pi()/2)
    return TMath::Pi() - phi;
  else
    return phi;
}


void AliEPDependentDiHadronOnTheFlyMCTask::Terminate(Option_t* option)
{
  // Terminate analysis  // Nothing to terminate... yet
  // Nothing to terminate... yet
  
}
