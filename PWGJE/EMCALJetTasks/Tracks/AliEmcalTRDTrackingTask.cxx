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

#include "AliEmcalTRDTrackingTask.h"

ClassImp(AliEmcalTRDTrackingTask)

/**
 * Default constructor
 */
AliEmcalTRDTrackingTask::AliEmcalTRDTrackingTask() : 
  AliAnalysisTaskEmcalLight("AliEmcalTrackingQA", kTRUE),
  fIsEsd(kFALSE),
  fGeneratorLevel(nullptr),
  fDetectorLevel(nullptr),
  fPtHistBins(),
  fPtResHistBins(),
  fIntegerHistBins(),
  fChargeHistBins(),
  fTracks(nullptr)
{
  SetMakeGeneralHistograms(kTRUE);

  GenerateHistoBins();
}

/**
 * Standard constructor
 */
AliEmcalTRDTrackingTask::AliEmcalTRDTrackingTask(const char *name) : 
  AliAnalysisTaskEmcalLight("AliEmcalTrackingQA", kTRUE),
  fIsEsd(kFALSE),
  fGeneratorLevel(nullptr),
  fDetectorLevel(nullptr),
  fPtHistBins(),
  fPtResHistBins(),
  fIntegerHistBins(),
  fChargeHistBins(),
  fTracks(nullptr)
{
  SetMakeGeneralHistograms(kTRUE);

  GenerateHistoBins();
}

/**
 * Destructor
 */
AliEmcalTRDTrackingTask::~AliEmcalTRDTrackingTask()
{
}

/**
 * Generate histogram bins
 */
void AliEmcalTRDTrackingTask::GenerateHistoBins()
{
  GenerateFixedBinArray(6,   0.0,   0.3, fPtHistBins, false);
  GenerateFixedBinArray(7,   0.3,   1.0, fPtHistBins, false);
  GenerateFixedBinArray(10,  1.0,   3.0, fPtHistBins, false);
  GenerateFixedBinArray(14,  3.0,  10.0, fPtHistBins, false);
  GenerateFixedBinArray(10, 10.0,  20.0, fPtHistBins, false);
  GenerateFixedBinArray(15, 20.0,  50.0, fPtHistBins, false);
  GenerateFixedBinArray(40, 50.0, 250.0, fPtHistBins, false);
  GenerateFixedBinArray(10, 250.0, 350.0, fPtHistBins);


  GenerateFixedBinArray(50, 0.00, 0.05, fPtResHistBins, false);
  GenerateFixedBinArray(25, 0.05, 0.10, fPtResHistBins, false);
  GenerateFixedBinArray(25, 0.10, 0.20, fPtResHistBins, false);
  GenerateFixedBinArray(30, 0.20, 0.50, fPtResHistBins, false);
  GenerateFixedBinArray(25, 0.50, 1.00, fPtResHistBins, false);
  GenerateFixedBinArray(20, 1.00, 2.00, fPtResHistBins);

  GenerateFixedBinArray(10, -0.5, 9.5, fIntegerHistBins);
  
  GenerateFixedBinArray(2, -1.1, 1.1, fChargeHistBins, true);
}

/**
 * Create histograms
 */
void AliEmcalTRDTrackingTask::UserCreateOutputObjects()
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
}

void AliEmcalTRDTrackingTask::ExecOnce()
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
THnSparse* AliEmcalTRDTrackingTask::GenerateTHnSparse(const char* name, const std::vector<std::tuple<std::string, std::vector<Double_t>::iterator, std::vector<Double_t>::iterator>>& axis)
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
void AliEmcalTRDTrackingTask::AllocateDetectorLevelTHnSparse()
{
  typedef std::vector<Double_t>::iterator my_iterator;

  std::vector<std::tuple<std::string, my_iterator, my_iterator>> axis;


  axis.push_back(std::make_tuple("#it{p}_{T} (GeV/#it{c})", fPtHistBins.begin(), fPtHistBins.end()));
  axis.push_back(std::make_tuple("track type", fIntegerHistBins.begin(), fIntegerHistBins.begin() + 4));
  axis.push_back(std::make_tuple("number of tracklets", fIntegerHistBins.begin(), fIntegerHistBins.begin()+7));
  axis.push_back(std::make_tuple("#sigma(#it{p}_{T}) / #it{p}_{T}", fPtResHistBins.begin(), fPtResHistBins.end()));
  axis.push_back(std::make_tuple("charge", fChargeHistBins.begin(), fChargeHistBins.end()));

  fTracks = GenerateTHnSparse("fTracks", axis);

  fOutput->Add(fTracks);
}

/**
 * Fill THnSparse with tracks
 */
void AliEmcalTRDTrackingTask::FillDetectorLevelTHnSparse(Double_t trackPt, Double_t sigma1OverPt, Byte_t trackType, Int_t ntracklets, Double_t trackCharge) 
{
  AliDebugStream(10) << "Filling detector level THnSparse" << std::endl;
  std::vector<Double_t> contents(fTracks->GetNdimensions());

  for (Int_t i = 0; i < fTracks->GetNdimensions(); i++) {
    TString title(fTracks->GetAxis(i)->GetTitle());
    if (title=="#it{p}_{T} (GeV/#it{c})")
      contents[i] = trackPt;
    else if (title=="#sigma(#it{p}_{T}) / #it{p}_{T}")
      contents[i] = sigma1OverPt*trackPt;
    else if (title=="track type")
      contents[i] = trackType;
    else if (title == "number of tracklets")
      contents[i] = ntracklets;
    else if (title=="charge")
      contents[i] = trackCharge;
    else 
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), fTracks->GetName()));
  }

  fTracks->Fill(&contents[0]);
}

/**
 * Fill all the histograms
 */
Bool_t AliEmcalTRDTrackingTask::FillHistograms()
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

      // Gold condition:
      // - at least 3 TRD tracklets (with this cut track without TRD in global track fit is at % level)
      if(!(track->GetStatus() & AliVTrack::kTRDupdate)) type += 2;

      Int_t label = TMath::Abs(track->GetLabel());
      AliDebugStream(10) << "Track " << trackIterator.current_index() << " with type " << int(type) << " and label " << label <<
          ", pt = " << track->Pt() << std::endl;

      double pt = track->Pt();

      FillDetectorLevelTHnSparse(pt, sigma, type, ntracklets, track->Charge());
    }
    else {
      AliErrorStream() << "Track " << trackIterator.current_index() << " has type " << type << " not recognized!" << std::endl;
    }
  }

  return kTRUE;
}

/**
 * Add task macro
 * @param isMC Whether it is an MC analysis
 * @return Pointer to the newly created object
 */
AliEmcalTRDTrackingTask* AliEmcalTRDTrackingTask::AddTaskTRDTracking(const char *suffix)
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
  TString name("AliEmcalTRDTrackingTask");
  if(strlen(suffix)) name += TString::Format("_%s", suffix);
  AliInfoClassStream() << "Allocating task." << std::endl;
  AliEmcalTRDTrackingTask *qaTask = new AliEmcalTRDTrackingTask(name);
  qaTask->SetVzRange(-10,10);
  AliInfoClassStream() << "Task allocated, setting containers." << std::endl;
  qaTask->AddParticleContainer(track_name.Data(), "detector");

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

