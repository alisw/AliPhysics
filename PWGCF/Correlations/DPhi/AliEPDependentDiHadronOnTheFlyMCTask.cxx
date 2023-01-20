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
AliAnalysisTaskSE(),
fListOutputHistograms(0x0),
fProfileFractionPrimaryTracks(0x0),
fHistNoEvents(0x0),
fHistImpactVsMultiplicity(0x0),
fHistImpactVsExperimentalMultiplicity(0x0),
fHistReactionPlane(0x0),
fHistEventReactionPlane(0x0),
fHistTrueEventReactionPlane(0x0),
fHistTrueEvent3Plane(0x0),
fHistTrueEvent4Plane(0x0),
fHistSingleParticles(0x0),
fHistos(0x0),
fHistosMixed(0x0),
fHistosIn(0x0),
fHistosInMixed(0x0),
fHistosMid(0x0),
fHistosMidMixed(0x0),
fHistosOut(0x0),
fHistosOutMixed(0x0),
fEtaMinV0A(2.8),
fEtaMaxV0A(5.1),
fEtaMinV0C(-3.7),
fEtaMaxV0C(-1.7),
fTrackEtaCut(1e9),
fMultiplicityEtaCut(0.5),
fMultiplicityPtCut(0.15),
fNoSectorsV0(8),
fFillMixed(kTRUE),
fMixingTracks(5000),
fPoolMgr(0x0)
//fCentralityEstimator("V0M")
{
  // Default Constructor (used for streaming)
}

AliEPDependentDiHadronOnTheFlyMCTask::AliEPDependentDiHadronOnTheFlyMCTask(const char* name):
AliAnalysisTaskSE(name),
fListOutputHistograms(0x0),
fProfileFractionPrimaryTracks(0x0),
fHistNoEvents(0x0),
fHistImpactVsMultiplicity(0x0),
fHistImpactVsExperimentalMultiplicity(0x0),
fHistReactionPlane(0x0),
fHistEventReactionPlane(0x0),
fHistTrueEventReactionPlane(0x0),
fHistTrueEvent3Plane(0x0),
fHistTrueEvent4Plane(0x0),
fHistSingleParticles(0x0),
fHistos(0x0),
fHistosMixed(0x0),
fHistosIn(0x0),
fHistosInMixed(0x0),
fHistosMid(0x0),
fHistosMidMixed(0x0),
fHistosOut(0x0),
fHistosOutMixed(0x0),
fEtaMinV0A(2.8),
fEtaMaxV0A(5.1),
fEtaMinV0C(-3.7),
fEtaMaxV0C(-1.7),
fTrackEtaCut(1e9),
fMultiplicityEtaCut(0.5),
fMultiplicityPtCut(0.15),
fNoSectorsV0(8),
fFillMixed(kTRUE),
fMixingTracks(5000),
fPoolMgr(0x0)
//fCentralityEstimator("V0M")
{
  // Constructor with name
  AliInfo("Constructing Event Plane based on fictional VZERO including NAME\n");
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
  
  fHistImpactVsExperimentalMultiplicity = new TH2D("impactVsLogExperimentalMultiplicity", "Impact parameter versus 10-Log of Multiplicity in experimental conditions", 40, 0, 20, 20, 0, 5);
  fListOutputHistograms->Add(fHistImpactVsExperimentalMultiplicity);
  
  fProfileFractionPrimaryTracks = new TProfile("fractionPrimaryTracks", "fractionPrimaryTracks", 1, onebinlimits, "s");
  fListOutputHistograms->Add(fProfileFractionPrimaryTracks);

  fHistReactionPlane = new TH1D("ReactionPlane", "Reactionplane", 200, -2*TMath::Pi(), 2*TMath::Pi());
  fListOutputHistograms->Add(fHistReactionPlane);

  fHistEventReactionPlane = new TH1D("EventReactionPlane", "Two-event plane w.r.t. the two-reactionplane", 100, 0.0, TMath::Pi());
  fListOutputHistograms->Add(fHistEventReactionPlane);

  fHistTrueEventReactionPlane = new TH1D("TrueEventReactionPlane", "All particle two-event plane w.r.t. the two-reactionplane", 100, 0.0, TMath::Pi());
  fListOutputHistograms->Add(fHistTrueEventReactionPlane);

  fHistTrueEvent3Plane = new TH1D("TrueEvent3Plane", "All particle three-event plane w.r.t. the all particle two-event plane", 200, 0.0, 2*TMath::Pi());
  fListOutputHistograms->Add(fHistTrueEvent3Plane);

  fHistTrueEvent4Plane = new TH1D("TrueEvent4Plane", "All particle three-event plane w.r.t. the all particle two-event plane", 200, 0.0, 2*TMath::Pi());
  fListOutputHistograms->Add(fHistTrueEvent4Plane);

  fHistSingleParticles = new TH2D("SingleParticleDistribution", "Distribution of all particles that are used in the dihadron distributions", 128, -TMath::Pi(), 3*TMath::Pi(), 80, -4.0, 4.0);
  fListOutputHistograms->Add(fHistSingleParticles);

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

  if(fFillMixed) {
    fHistosMixed = new AliUEHistograms("AliUEHistogramsMixed-1", "4R", fCustomBinning);
    fHistosMixed->SetTrackEtaCut(fTrackEtaCut);
    fListOutputHistograms->Add(fHistosMixed);

    fHistosInMixed = new AliUEHistograms("AliUEHistogramsMixed0", "4R", fCustomBinning);
    fHistosInMixed->SetTrackEtaCut(fTrackEtaCut);
    fListOutputHistograms->Add(fHistosInMixed);

    fHistosMidMixed = new AliUEHistograms("AliUEHistogramsMixed1", "4R", fCustomBinning);
    fHistosMidMixed->SetTrackEtaCut(fTrackEtaCut);
    fListOutputHistograms->Add(fHistosMidMixed);

    fHistosOutMixed = new AliUEHistograms("AliUEHistogramsMixed2", "4R", fCustomBinning);
    fHistosOutMixed->SetTrackEtaCut(fTrackEtaCut);
    fListOutputHistograms->Add(fHistosOutMixed);
  }

  if(!fPoolMgr)
  {
    Int_t nCentralityBins = fHistos->GetUEHist(2)->GetEventHist()->GetNBins(1);
    Double_t* centralityBins = (Double_t*) fHistos->GetUEHist(2)->GetEventHist()->GetAxis(1,0)->GetXbins()->GetArray();
    Int_t nZvtxBins = 1;
    Double_t zvtxbins[2] = {-999., 999.};
    Int_t nPsiBins = 1;
    Double_t psibins[2] = {-999., 999.};
    //Int_t nPtBins = 1;
    Int_t nPtBins = fHistos->GetUEHist(2)->GetTrackHist(AliUEHist::kToward)->GetNBins(1);
    //Double_t defaultPtBins[2] = {-9999., 9999.};
    //Double_t* ptbins = defaultPtBins;
    Double_t* ptbins = (Double_t*) fHistos->GetUEHist(2)->GetTrackHist(AliUEHist::kToward)->GetAxis(1,0)->GetXbins()->GetArray();
    fPoolMgr = new AliEventPoolManager(1000, fMixingTracks, nCentralityBins, centralityBins, nZvtxBins, zvtxbins, nPsiBins, psibins, nPtBins, ptbins);
    fPoolMgr->SetTargetValues(fMixingTracks, 0.1, 5);
  }

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


  // Cut on central barrel acceptancy
  TObjArray* etacutTracks = SubSelectTracksEta(allTracks);

  // Divide the leading tracks with respect to the event plane. Done in a single loop to save on computation time
  TObjArray* inplaneTracks = new TObjArray; // cutOnEventPlaneAngle(allTracks, 0.0, TMath::Pi()/6);
  TObjArray* midplaneTracks = new TObjArray; // cutOnEventPlaneAngle(allTracks, TMath::Pi()/6, TMath::Pi()/3);
  TObjArray* outplaneTracks = new TObjArray; // cutOnEventPlaneAngle(allTracks, TMath::Pi()/3, TMath::Pi()/2);
  for(Int_t part = 0; part < etacutTracks->GetEntries(); part++) {
    AliMCParticle* track = dynamic_cast<AliMCParticle*>(etacutTracks->At(part));
    Double_t dphi = piToHalfPi(twoPiToPi(angleBetween(track->Phi(), reactionplane)));
    if(dphi < TMath::Pi()/6) inplaneTracks->Add(track);
    if((dphi >= TMath::Pi()/6) && (dphi < TMath::Pi()/3)) midplaneTracks->Add(track);
    if(dphi >= TMath::Pi()/3) outplaneTracks->Add(track);
  }

  fHistos->FillCorrelations(centrality, zvtx, AliUEHist::kCFStepAll, etacutTracks, etacutTracks, 1);
  fHistosIn->FillCorrelations(centrality, zvtx, AliUEHist::kCFStepAll, inplaneTracks, etacutTracks, 1);
  fHistosMid->FillCorrelations(centrality, zvtx, AliUEHist::kCFStepAll, midplaneTracks, etacutTracks, 1);
  fHistosOut->FillCorrelations(centrality, zvtx, AliUEHist::kCFStepAll, outplaneTracks, etacutTracks, 1);

  if(fFillMixed) {
    for(Int_t iPool=0; iPool<fPoolMgr->GetNumberOfPtBins(); iPool++) {
      AliEventPool* pool = fPoolMgr->GetEventPool(centrality, zvtx, 0., iPool);
      if(pool->IsReady()) {
        for(Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++) {
          fHistosMixed->FillCorrelations(centrality, zvtx, AliUEHist::kCFStepAll, etacutTracks, pool->GetEvent(jMix), 1.0 / pool->GetCurrentNEvents(), (jMix==0));
          fHistosMixed->FillCorrelations(centrality, zvtx, AliUEHist::kCFStepAll, inplaneTracks, pool->GetEvent(jMix), 1.0 / pool->GetCurrentNEvents(), (jMix==0));
          fHistosMixed->FillCorrelations(centrality, zvtx, AliUEHist::kCFStepAll, midplaneTracks, pool->GetEvent(jMix), 1.0 / pool->GetCurrentNEvents(), (jMix==0));
          fHistosMixed->FillCorrelations(centrality, zvtx, AliUEHist::kCFStepAll, outplaneTracks, pool->GetEvent(jMix), 1.0 / pool->GetCurrentNEvents(), (jMix==0));
        }
      }
      pool->UpdatePool(CloneAndReduceTrackList(etacutTracks, pool->GetPtMin(), pool->GetPtMax()));
    }
  }

  delete allTracks; delete etacutTracks; delete inplaneTracks; delete midplaneTracks; delete outplaneTracks;
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

  if(!mcEvent) {
    AliError("No correct mcEvent provided to SelectTracks");
    return 0x0;
  }

  Int_t experimentalMultiplicity = 0;

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
    
      if( track->Charge() == 0 ) continue;

      acceptedTracks->Add(track);

      if((TMath::Abs(track->Eta()) < fMultiplicityEtaCut) & (track->Pt() > fMultiplicityPtCut)) experimentalMultiplicity++;
    }
    
    AliGenHijingEventHeader* eventHeader = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());
    Double_t impactParameter = (Double_t) eventHeader->ImpactParameter();
    
    fHistImpactVsExperimentalMultiplicity->Fill(impactParameter, TMath::Log10(experimentalMultiplicity));
  }

  return acceptedTracks;
}

TObjArray* AliEPDependentDiHadronOnTheFlyMCTask::SubSelectTracksEta(TObjArray* selectedTracks)
{
  TObjArray* subselectedTracks = new TObjArray;

  if(!selectedTracks) {
    AliError("TObjArray passed to SubSelectTracksEta is invalid.");
    return 0x0;
  }

  if(selectedTracks) {
    Int_t nTracks = selectedTracks->GetEntries();
    for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
      AliMCParticle* track = dynamic_cast<AliMCParticle*>(selectedTracks->At(iTrack));
      if(!track) {
        AliError(Form("Object number %d in the selectedTracks was not an AliMCParticle.", iTrack));
        continue;
      }

      if(TMath::Abs(track->Eta()) > fTrackEtaCut) continue;

      fHistSingleParticles->Fill(track->Phi(), track->Eta());
      subselectedTracks->Add(track);
    }
  }

  return subselectedTracks;
}


TObjArray* AliEPDependentDiHadronOnTheFlyMCTask::CloneAndReduceTrackList(TObjArray* tracks, Double_t minPt, Double_t maxPt)
{
  // Clones the track list with AliBasicParticles that use less memory (for event mixing), copied from AliAnalysisTaskPhiCorrelations::CloneAndReduceTrackList

  // Check if we already have a reduced track list. In that case a simple Clone is enough
  // Only possible for a inclusive pt
  if(maxPt-minPt < 0)
    if(tracks->GetEntriesFast() == 0 || tracks->UncheckedAt(0)->InheritsFrom("AliBasicParticle"))
      return (TObjArray*) tracks->Clone();

  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);

  for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
    AliVParticle* particle = (AliVParticle*) tracks->UncheckedAt(i);
    AliBasicParticle* copy = 0;

    if( (maxPt-minPt > 0) && ((particle->Pt()<minPt) || (particle->Pt()>=maxPt)) )
      continue;

    copy = new AliBasicParticle(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge());
    copy->SetUniqueID(particle->GetUniqueID());
    tracksClone->Add(copy);
  }

  return tracksClone;
}


void AliEPDependentDiHadronOnTheFlyMCTask::GetEventDetails(AliMCEvent* mcEvent, Double_t& reactionplane, Double_t& zvtx, Double_t& centrality)
{
  AliGenHijingEventHeader* eventHeader = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());

  reactionplane = (Double_t) eventHeader->ReactionPlaneAngle();

  TArrayF vertex = TArrayF(3);
  eventHeader->PrimaryVertex(vertex);
  zvtx = vertex.At(2);

/* TODO: code to implement centrality in the future
  AliMultSelection* multSelection = (AliMultSelection*) mcEvent->FindListObject("MultSelection");
  if (!multSelection)
    AliFatal("MultSelection not found in input event. Did the AliMultSelectionTask run before this task?");
*/

  centrality = 35.0; // TODO: centrality not implemented: multSelection->GetMultiplicityPercentile(fCentralityEstimator, kTRUE);
}


void AliEPDependentDiHadronOnTheFlyMCTask::FillEventPlaneHistograms(TObjArray* allTracks, Double_t reactionplane, Double_t& eventplaneV0)
{
  TH1D* V0A = FillVirtualV0(allTracks, fEtaMinV0A, fEtaMaxV0A);
  TH1D* V0C = FillVirtualV0(allTracks, fEtaMinV0C, fEtaMaxV0C);
  TH1D* V0 = new TH1D("V0dist", "SimulatedV0", fNoSectorsV0, 0, 2*TMath::Pi());
  V0->Add(V0A, V0C);
  Double_t trueeventplane = ComputeTrueEventPlane(allTracks, 2);
  Double_t trueevent3plane = ComputeTrueEventPlane(allTracks, 3);
  Double_t trueevent4plane = ComputeTrueEventPlane(allTracks, 4);
  Double_t eventplaneV0A = ComputeEventPlane(V0A, 2);
  Double_t eventplaneV0C = ComputeEventPlane(V0C, 2);
  eventplaneV0 = ComputeEventPlane(V0, 2);
  Double_t angleBetweenPlanes = twoPiToPi(angleBetween(reactionplane, eventplaneV0));
  Double_t angleBetweenPlanesTrue = twoPiToPi(angleBetween(trueeventplane, reactionplane));
  Double_t angleBetweenPlanes3True = angleBetween(trueevent3plane, trueeventplane);
  Double_t angleBetweenPlanes4True = angleBetween(trueevent4plane, trueeventplane);

  fHistReactionPlane->Fill(reactionplane);
  fHistEventReactionPlane->Fill(angleBetweenPlanes);
  fHistTrueEventReactionPlane->Fill(angleBetweenPlanesTrue);
  fHistTrueEvent3Plane->Fill(angleBetweenPlanes3True);
  fHistTrueEvent4Plane->Fill(angleBetweenPlanes4True);

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
    if(TMath::Abs(track->Eta()) > fTrackEtaCut) continue;
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
