#include <cstring>
#include <iostream>

#include <THnSparse.h>
#include <TMath.h>
#include <TString.h>

#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include <AliAODMCParticle.h>
#include <AliLog.h>
#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>

#include "AliMCParticleContainer.h"

#include "AliEmcalTrackingQATask.h"

ClassImp(AliEmcalTrackingQATask)

/**
 * Default constructor
 */
AliEmcalTrackingQATask::AliEmcalTrackingQATask() : 
  AliAnalysisTaskEmcalLight("AliEmcalTrackingQA", kTRUE),
  fDoSigma1OverPt(kFALSE),
  fDoSigmaPtOverPtGen(kFALSE),
  fDoSeparateTRDrefit(kFALSE),
  fUseTRDUpdateFlag(kTRUE),
  fUseQOverPtShift(kFALSE),
  fQOverPtShift(0),
  fIsEsd(kFALSE),
  fGeneratorLevel(nullptr),
  fDetectorLevel(nullptr),
  fPtHistBins(),
  fEtaHistBins(),
  fPhiHistBins(),
  fCentHistBins(),
  fPtRelDiffHistBins(),
  fPtResHistBins(),
  f1OverPtResHistBins(),
  fIntegerHistBins(),
  fChargeHistBins(),
  fTracks(nullptr),
  fParticlesPhysPrim(nullptr),
  fParticlesMatched(nullptr)
{
  SetMakeGeneralHistograms(kTRUE);

  GenerateHistoBins();
}

/**
 * Standard constructor
 */
AliEmcalTrackingQATask::AliEmcalTrackingQATask(const char *name) : 
  AliAnalysisTaskEmcalLight("AliEmcalTrackingQA", kTRUE),
  fDoSigma1OverPt(kFALSE),
  fDoSigmaPtOverPtGen(kFALSE),
  fDoSeparateTRDrefit(kFALSE),
  fUseTRDUpdateFlag(kTRUE),
  fUseQOverPtShift(kFALSE),
  fQOverPtShift(0),
  fIsEsd(kFALSE),
  fGeneratorLevel(nullptr),
  fDetectorLevel(nullptr),
  fPtHistBins(),
  fEtaHistBins(),
  fPhiHistBins(),
  fCentHistBins(),
  fPtRelDiffHistBins(),
  fPtResHistBins(),
  f1OverPtResHistBins(),
  fIntegerHistBins(),
  fChargeHistBins(),
  fTracks(nullptr),
  fParticlesPhysPrim(nullptr),
  fParticlesMatched(nullptr)
{
  SetMakeGeneralHistograms(kTRUE);

  GenerateHistoBins();
}

/**
 * Destructor
 */
AliEmcalTrackingQATask::~AliEmcalTrackingQATask()
{
}

/**
 * Generate histogram bins
 */
void AliEmcalTrackingQATask::GenerateHistoBins()
{
  GenerateFixedBinArray(6,   0.0,   0.3, fPtHistBins, false);
  GenerateFixedBinArray(7,   0.3,   1.0, fPtHistBins, false);
  GenerateFixedBinArray(10,  1.0,   3.0, fPtHistBins, false);
  GenerateFixedBinArray(14,  3.0,  10.0, fPtHistBins, false);
  GenerateFixedBinArray(10, 10.0,  20.0, fPtHistBins, false);
  GenerateFixedBinArray(15, 20.0,  50.0, fPtHistBins, false);
  GenerateFixedBinArray(40, 50.0, 250.0, fPtHistBins, false);
  GenerateFixedBinArray(10, 250.0, 350.0, fPtHistBins);

  GenerateFixedBinArray(100, -1.0, 1.0, fEtaHistBins);

  GenerateFixedBinArray(100, 0.0, TMath::TwoPi(), fPhiHistBins);

  fCentHistBins.push_back(0);
  fCentHistBins.push_back(10);
  fCentHistBins.push_back(30);
  fCentHistBins.push_back(50);
  fCentHistBins.push_back(90);

  GenerateFixedBinArray(50, 0.00, 0.05, fPtResHistBins, false);
  GenerateFixedBinArray(25, 0.05, 0.10, fPtResHistBins, false);
  GenerateFixedBinArray(25, 0.10, 0.20, fPtResHistBins, false);
  GenerateFixedBinArray(30, 0.20, 0.50, fPtResHistBins, false);
  GenerateFixedBinArray(25, 0.50, 1.00, fPtResHistBins, false);
  GenerateFixedBinArray(20, 1.00, 2.00, fPtResHistBins);

  GenerateFixedBinArray(200, -2.0, 2.0, fPtRelDiffHistBins);

  GenerateFixedBinArray(100, 0.00, 0.02, f1OverPtResHistBins, false);
  GenerateFixedBinArray( 60, 0.02, 0.05, f1OverPtResHistBins, false);
  GenerateFixedBinArray( 50, 0.05, 0.10, f1OverPtResHistBins, false);
  GenerateFixedBinArray( 50, 0.10, 0.20, f1OverPtResHistBins, false);
  GenerateFixedBinArray( 75, 0.20, 0.50, f1OverPtResHistBins, false);
  GenerateFixedBinArray( 80, 0.50, 1.50, f1OverPtResHistBins);

  GenerateFixedBinArray(10, -0.5, 9.5, fIntegerHistBins);
  
  GenerateFixedBinArray(2, -1.1, 1.1, fChargeHistBins, true);
}

/**
 * Create histograms
 */
void AliEmcalTrackingQATask::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalLight::UserCreateOutputObjects();

  if (fParticleCollArray.empty()) {
    AliErrorStream() << "This task needs at least one particle container! The task won't run." << std::endl;
    fInhibit = kTRUE;
    return;
  }

  AliDebugStream(3) << "Loading the detector track container" << std::endl;
  if (!fDetectorLevel) {
    auto iter = fParticleCollArray.find("detector");
    if (iter == fParticleCollArray.end()) {
      AliErrorStream() << "This task needs at least one particle container named 'detector'! The task won't run." << std::endl;
      fInhibit = kTRUE;
      return;
    }
    else {
      fDetectorLevel = static_cast<AliTrackContainer*>(iter->second);
    }
  }

  AliDebugStream(3) << "Loading the generator particle container" << std::endl;
  if (!fGeneratorLevel) {
    auto iter = fParticleCollArray.find("generator");
    if (iter == fParticleCollArray.end()) {
      AliInfoStream() << "No particle container named 'generator' was found. Assuming this is not a MC production." << std::endl;
    }
    else {
      fGeneratorLevel = static_cast<AliMCParticleContainer*>(iter->second);
    }
  }

  AliDebugStream(3) << "Allocating histograms" << std::endl;
  AllocateDetectorLevelTHnSparse();

  if (fGeneratorLevel) {
    AllocateGeneratorLevelTHnSparse();
    AllocateMatchedParticlesTHnSparse();
  }
}

void AliEmcalTrackingQATask::ExecOnce()
{
  AliAnalysisTaskEmcalLight::ExecOnce();
  if (!fDetectorLevel->GetArray()) {
    AliErrorStream() << "Could not load track array! The task won't run." << std::endl;
    fInhibit = kTRUE;
    return;
  }
  if (fDetectorLevel->GetArray()->GetClass()->InheritsFrom("AliESDtrack")) fIsEsd = kTRUE;
}

/**
 * Generate a THnSparseF based on a list of axis provided
 * @param name Name of the output histogram
 * @param axis Vector containing the axis
 * @return The newly created THnSparse
 */
THnSparse* AliEmcalTrackingQATask::GenerateTHnSparse(const char* name, const std::vector<std::tuple<std::string, std::vector<Double_t>::iterator, std::vector<Double_t>::iterator>>& axis)
{
  std::vector<int> nbins;
  for (auto a : axis) nbins.push_back(int(std::get<2>(a) - std::get<1>(a) - 1));

  THnSparse* h = new THnSparseF(name, name, nbins.size(), &nbins[0]);
  Int_t i = 0;
  for (auto a : axis) {
    h->GetAxis(i)->SetTitle(std::get<0>(a).c_str());
    h->SetBinEdges(i, &(*(std::get<1>(a))));
    i++;
  }

  return h;
}

/**
 * Allocate THnSparse to contain tracks
 */
void AliEmcalTrackingQATask::AllocateDetectorLevelTHnSparse()
{
  typedef std::vector<Double_t>::iterator my_iterator;

  std::vector<std::tuple<std::string, my_iterator, my_iterator>> axis;

  if (fForceBeamType != AliAnalysisTaskEmcalLight::kpp) {
    axis.push_back(std::make_tuple("Centrality %", fCentHistBins.begin(), fCentHistBins.end()));
  }

  axis.push_back(std::make_tuple("#it{p}_{T} (GeV/#it{c})", fPtHistBins.begin(), fPtHistBins.end()));
  axis.push_back(std::make_tuple("#eta", fEtaHistBins.begin(), fEtaHistBins.end()));
  axis.push_back(std::make_tuple("#phi", fPhiHistBins.begin(), fPhiHistBins.end()));
  axis.push_back(std::make_tuple("MC Generator", fIntegerHistBins.begin(), fIntegerHistBins.begin() + 3));
  axis.push_back(std::make_tuple("track type", fIntegerHistBins.begin(), fIntegerHistBins.begin() + (fDoSeparateTRDrefit ? 9 : 5)));

  if (fDoSigma1OverPt) {
    axis.push_back(std::make_tuple("#sigma(1/#it{p}_{T}) (GeV/#it{c})^{-1}", f1OverPtResHistBins.begin(), f1OverPtResHistBins.end()));
  }
  else {
    axis.push_back(std::make_tuple("#sigma(#it{p}_{T}) / #it{p}_{T}", fPtResHistBins.begin(), fPtResHistBins.end()));
  }
  axis.push_back(std::make_tuple("charge", fChargeHistBins.begin(), fChargeHistBins.end()));

  fTracks = GenerateTHnSparse("fTracks", axis);

  fOutput->Add(fTracks);
}

/**
 * Allocate THnSparse to contain particles (generator level)
 */
void AliEmcalTrackingQATask::AllocateGeneratorLevelTHnSparse()
{
  typedef std::vector<Double_t>::iterator my_iterator;

  std::vector<std::tuple<std::string, my_iterator, my_iterator>> axis;

  if (fForceBeamType != AliAnalysisTaskEmcalLight::kpp) {
    axis.push_back(std::make_tuple("Centrality %", fCentHistBins.begin(), fCentHistBins.end()));
  }

  axis.push_back(std::make_tuple("#it{p}_{T} (GeV/#it{c})", fPtHistBins.begin(), fPtHistBins.end()));
  axis.push_back(std::make_tuple("#eta", fEtaHistBins.begin(), fEtaHistBins.end()));
  axis.push_back(std::make_tuple("#phi", fPhiHistBins.begin(), fPhiHistBins.end()));
  axis.push_back(std::make_tuple("MC Generator", fIntegerHistBins.begin(), fIntegerHistBins.begin() + 3));
  axis.push_back(std::make_tuple("Findable", fIntegerHistBins.begin(), fIntegerHistBins.begin() + 3));
  axis.push_back(std::make_tuple("charge", fChargeHistBins.begin(), fChargeHistBins.end()));

  fParticlesPhysPrim = GenerateTHnSparse("fParticlesPhysPrim", axis);

  fOutput->Add(fParticlesPhysPrim);
}

/**
 * Allocate THnSparse to contain tracks matched to particles
 */
void AliEmcalTrackingQATask::AllocateMatchedParticlesTHnSparse()
{
  typedef std::vector<Double_t>::iterator my_iterator;

  std::vector<std::tuple<std::string, my_iterator, my_iterator>> axis;

  if (fForceBeamType != AliAnalysisTaskEmcalLight::kpp) {
    axis.push_back(std::make_tuple("Centrality %", fCentHistBins.begin(), fCentHistBins.end()));
  }

  axis.push_back(std::make_tuple("#it{p}_{T}^{gen} (GeV/#it{c})", fPtHistBins.begin(), fPtHistBins.end()));
  axis.push_back(std::make_tuple("#eta^{gen}", fEtaHistBins.begin(), fEtaHistBins.end()));
  axis.push_back(std::make_tuple("#phi^{gen}", fPhiHistBins.begin(), fPhiHistBins.end()));
  axis.push_back(std::make_tuple("#it{p}_{T}^{det} (GeV/#it{c})", fPtHistBins.begin(), fPtHistBins.end()));
  axis.push_back(std::make_tuple("#eta^{det}", fEtaHistBins.begin(), fEtaHistBins.end()));
  axis.push_back(std::make_tuple("#phi^{det}", fPhiHistBins.begin(), fPhiHistBins.end()));
  axis.push_back(std::make_tuple("track type", fIntegerHistBins.begin(), fIntegerHistBins.begin() + (fDoSeparateTRDrefit ? 9 : 5)));

  if (fDoSigma1OverPt) {
    axis.push_back(std::make_tuple("(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{gen}", fPtRelDiffHistBins.begin(), fPtRelDiffHistBins.end()));
  }
  else {
    axis.push_back(std::make_tuple("(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{det}", fPtRelDiffHistBins.begin(), fPtRelDiffHistBins.end()));
  }
  axis.push_back(std::make_tuple("charge", fChargeHistBins.begin(), fChargeHistBins.end()));

  fParticlesMatched = GenerateTHnSparse("fParticlesMatched", axis);

  fOutput->Add(fParticlesMatched);
}

/**
 * Fill THnSparse with tracks
 */
void AliEmcalTrackingQATask::FillDetectorLevelTHnSparse(Double_t cent, Double_t trackEta, Double_t trackPhi, Double_t trackPt, 
    Double_t sigma1OverPt, Int_t mcGen, Byte_t trackType, Double_t trackCharge) 
{
  AliDebugStream(10) << "Filling detector level THnSparse" << std::endl;
  std::vector<Double_t> contents(fTracks->GetNdimensions());

  for (Int_t i = 0; i < fTracks->GetNdimensions(); i++) {
    TString title(fTracks->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = cent;
    else if (title=="#it{p}_{T} (GeV/#it{c})")
      contents[i] = trackPt;
    else if (title=="#eta")
      contents[i] = trackEta;
    else if (title=="#phi")
      contents[i] = trackPhi;
    else if (title=="#sigma(1/#it{p}_{T}) (GeV/#it{c})^{-1}")
      contents[i] = sigma1OverPt;
    else if (title=="#sigma(#it{p}_{T}) / #it{p}_{T}")
      contents[i] = sigma1OverPt*trackPt;
    else if (title=="MC Generator")
      contents[i] = mcGen;
    else if (title=="track type")
      contents[i] = trackType;
    else if (title=="charge")
      contents[i] = trackCharge;
    else 
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), fTracks->GetName()));
  }

  fTracks->Fill(&contents[0]);
}

/**
 * Fill THnSparse with particles
 */
void AliEmcalTrackingQATask::FillGeneratorLevelTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt, Int_t mcGen, Byte_t findable, Double_t partCharge)
{
  std::vector<Double_t> contents(fParticlesPhysPrim->GetNdimensions());

  for (Int_t i = 0; i < fParticlesPhysPrim->GetNdimensions(); i++) {
    TString title(fParticlesPhysPrim->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = cent;
    else if (title=="#it{p}_{T} (GeV/#it{c})")
      contents[i] = partPt;
    else if (title=="#eta")
      contents[i] = partEta;
    else if (title=="#phi")
      contents[i] = partPhi;
    else if (title=="MC Generator")
      contents[i] = mcGen;
    else if (title=="Findable")
      contents[i] = findable;
    else if (title=="charge")
      contents[i] = partCharge;
    else 
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), fParticlesPhysPrim->GetName()));
  }

  fParticlesPhysPrim->Fill(&contents[0]);
}

/**
 * Fill THnSparse with tracks matched to particles
 */
void AliEmcalTrackingQATask::FillMatchedParticlesTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt,
    Double_t trackEta, Double_t trackPhi, Double_t trackPt, Byte_t trackType, Double_t trackCharge)
{
  std::vector<Double_t> contents(fParticlesMatched->GetNdimensions());

  for (Int_t i = 0; i < fParticlesMatched->GetNdimensions(); i++) {
    TString title(fParticlesMatched->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = cent;
    else if (title=="#it{p}_{T}^{gen} (GeV/#it{c})")
      contents[i] = partPt;
    else if (title=="#eta^{gen}")
      contents[i] = partEta;
    else if (title=="#phi^{gen}")
      contents[i] = partPhi;
    else if (title=="#it{p}_{T}^{det} (GeV/#it{c})")
      contents[i] = trackPt;
    else if (title=="#eta^{det}")
      contents[i] = trackEta;
    else if (title=="#phi^{det}")
      contents[i] = trackPhi;
    else if (title=="(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{gen}")
      contents[i] = (partPt - trackPt) / partPt;
    else if (title=="(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{det}")
      contents[i] = (partPt - trackPt) / trackPt;
    else if (title=="track type")
      contents[i] = (Double_t)trackType;
    else if (title=="charge")
      contents[i] = trackCharge;
    else 
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), fParticlesMatched->GetName()));
  }

  fParticlesMatched->Fill(&contents[0]);
}

/**
 * Fill all the histograms
 */
Bool_t AliEmcalTrackingQATask::FillHistograms()
{
  AliDebugStream(1) << "Called: Tracks [" << fDetectorLevel->GetNTracks() << "], Accepted [" << fDetectorLevel->GetNAcceptedTracks() << "]\n";
  auto iterable = fDetectorLevel->accepted_momentum();
  for (auto trackIterator = iterable.begin(); trackIterator != iterable.end(); trackIterator++) {
    auto track = trackIterator->second;
    Byte_t type = fDetectorLevel->GetTrackType(track);
    AliDebugStream(2) << "Next track of type " << static_cast<Int_t>(type) << std::endl;
    Byte_t ntracklets = 0;
    if (type <= 3) {
      Double_t sigma = 0;
      
      if (fIsEsd) {
        AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(track);
        if (esdTrack){
          sigma = TMath::Sqrt(esdTrack->GetSigma1Pt2());
          ntracklets = esdTrack->GetTRDntracklets();
        }
      }
      else { // AOD
        AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track);
        if(!aodtrack) AliFatal("Not a standard AOD");
        
        AliExternalTrackParam exParam;
        
        //get covariance matrix
        Double_t cov[21] = {0,};
        aodtrack->GetCovMatrix(cov);
        Double_t pxpypz[3] = {0,};
        aodtrack->PxPyPz(pxpypz);
        Double_t xyz[3] = {0,};
        aodtrack->GetXYZ(xyz);
        Short_t sign = aodtrack->Charge();
        exParam.Set(xyz,pxpypz,cov,sign);

        sigma = TMath::Sqrt(exParam.GetSigma1Pt2());
        ntracklets = track->GetTRDntrackletsPID();
      }

      if(fDoSeparateTRDrefit) {
        // Gold condition:
        // - at least 3 TRD tracklets (with this cut track without TRD in global track fit is at % level)
        if(fUseTRDUpdateFlag) {
          if(!(track->GetStatus() & AliVTrack::kTRDupdate)) type += 4;
        } else  {
          if(ntracklets < 3) type += 4;    // failed TRD gold condition
        }
      }

      Int_t label = TMath::Abs(track->GetLabel());
      AliDebugStream(10) << "Track " << trackIterator.current_index() << " with type " << int(type) << " and label " << label <<
          ", pt = " << track->Pt() << std::endl;
      Int_t mcGen = 1;
      // reject particles generated from other generators in the cocktail but keep fake tracks (label == 0)
      if (label == 0 || track->GetGeneratorIndex() <= 0) mcGen = 0;

      double pt = track->Pt();
      if(fUseQOverPtShift) {
        double chargeval = track->Charge() > 0 ? 1. : -1.;
        pt = 1./(fQOverPtShift*chargeval  + 1./pt);
        AliDebugStream(1) <<  "Applying q/pt shift " << fQOverPtShift << ", before shift: " << track->Pt() << ", after shift: " << pt << std::endl;
      }

      FillDetectorLevelTHnSparse(fCent, track->Eta(), track->Phi(), pt, sigma, mcGen, type, track->Charge());

      if (fGeneratorLevel && label > 0) {
        AliAODMCParticle *part =  fGeneratorLevel->GetAcceptMCParticleWithLabel(label);
        if (part) {
          if (part->GetGeneratorIndex() <= 0) {
            Int_t pdg = TMath::Abs(part->PdgCode());
            // select charged pions, protons, kaons , electrons, muons
            if (pdg == 211 || pdg == 2212 || pdg == 321 || pdg == 11 || pdg == 13) {
              FillMatchedParticlesTHnSparse(fCent, part->Eta(), part->Phi(), part->Pt(), track->Eta(), track->Phi(), track->Pt(), type, track->Charge());
            }
          }
        }
      }
    }
    else {
      AliErrorStream() << "Track " << trackIterator.current_index() << " has type " << type << " not recognized!" << std::endl;
    }
  }

  if (fGeneratorLevel) {
    auto iterable = fGeneratorLevel->accepted_momentum();
    for (auto partIterator = iterable.begin(); partIterator != iterable.end(); partIterator++) {
      auto part = partIterator->second;

      Int_t mcGen = 1;
      Byte_t findable = 0;
      Double_t partcharge; // translate to +- 1
      if(part->Charge() > 0) partcharge = 1.;
      else partcharge = -1.;

      if (part->GetGeneratorIndex() <= 0) mcGen = 0;

      Int_t pdg = TMath::Abs(part->PdgCode());
      // select charged pions, protons, kaons , electrons, muons
      if (pdg == 211 || pdg == 2212 || pdg == 321 || pdg == 11 || pdg == 13) findable = 1;

      FillGeneratorLevelTHnSparse(fCent, part->Eta(), part->Phi(), part->Pt(), mcGen, findable, partcharge);
    }
  }

  return kTRUE;
}

/**
 * Add task macro
 * @param isMC Whether it is an MC analysis
 * @return Pointer to the newly created object
 */
AliEmcalTrackingQATask* AliEmcalTrackingQATask::AddTaskTrackingQA(Bool_t isMC, const char *suffix)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    AliErrorClassStream() << "No analysis manager to connect to." << std::endl;
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis manager
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) {
    AliErrorClassStream() << "This task requires an input event handler" << std::endl;
    return nullptr;
  }
  else {
    AliInfoClassStream() << "Event handler ok!" << std::endl;
  }

  EDataType_t dataType = kUnknownDataType;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
    AliInfoClassStream() << "Data type is ESD." << std::endl;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
    AliInfoClassStream() << "Data type is AOD." << std::endl;
  }

  TString track_name("tracks");
  if (dataType == kESD) track_name = "Tracks";

  // Init the task and do settings
  TString name("AliEmcalTrackingQATask");
  if(strlen(suffix)) name += TString::Format("_%s", suffix);
  AliInfoClassStream() << "Allocating task." << std::endl;
  AliEmcalTrackingQATask *qaTask = new AliEmcalTrackingQATask(name);
  qaTask->SetVzRange(-10,10);
  AliInfoClassStream() << "Task allocated, setting containers." << std::endl;
  qaTask->AddParticleContainer(track_name.Data(), "detector");
  if (isMC) qaTask->AddParticleContainer("mcparticles", "generator");

  AliInfoClassStream() << "Containers ok, adding task to the analysis manager." << std::endl;
  // Final settings, pass to manager and set the containers
  mgr->AddTask(qaTask);

  AliInfoClassStream() << "Task added, setting input/output." << std::endl;
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;

  TString contName(name);
  contName += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(),
                  TList::Class(),AliAnalysisManager::kOutputContainer,
                  Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (qaTask, 0,  cinput1 );
  mgr->ConnectOutput (qaTask, 1, coutput1 );


  AliInfoClassStream() << "Task configuration done." << std::endl;
  return qaTask;
}

