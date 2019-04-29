#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"

#include "AliMCEvent.h"

#include "AliMultDepSpecAnalysisTask.h"


/// \cond CLASSIMP
ClassImp(AliMultDepSpecAnalysisTask);
/// \endcond

/***************************************************************************//**
 * ROOT I/O Constructor.
 ******************************************************************************/
 AliMultDepSpecAnalysisTask::AliMultDepSpecAnalysisTask() : AliAnalysisTaskSE(),
   //General member variables
   fEventCuts(),
   fOutputList(nullptr),
   fEvent(nullptr),
   fMCEvent(nullptr),
   fESDtrackCuts(nullptr),
   fCutMode(100),
   //Toggles
   fIsESD(kTRUE),
   fIsMC(kFALSE),
   // Cut Parameters
   fTriggerMask(AliVEvent::kMB | AliVEvent::kINT7),
   fMinEta(-10),
   fMaxEta(10),
   fMinPt(0.0),
   fMaxPt(50.0),
   //Arrays for Binning
   fBinsMult(nullptr),
   fBinsCent(nullptr),
   fBinsPt(nullptr),
   fBinsEta(nullptr),
   fBinsZv(nullptr),
   fBinsPtReso(nullptr),
   //Event-Histograms
   fHistEventSelection(nullptr),
   fHistEvents(nullptr),
   fHistTracks(nullptr),
   fHistRelPtReso(nullptr),
   fHistMCRelPtReso(nullptr),
   fHistMCMultCorrelMatrix(nullptr),
   fHistMCPtCorrelMatrix(nullptr),
   fHistMCEtaCorrelMatrix(nullptr),
   fHistMCPrimTrue(nullptr),
   fHistMCPrimMeas(nullptr),
   fHistMCSecMeas(nullptr)
{
  // ROOT IO constructor, don't allocate memory here!
}


/***************************************************************************//**
 * Constructor.
 ******************************************************************************/
AliMultDepSpecAnalysisTask::AliMultDepSpecAnalysisTask(const char* name) : AliAnalysisTaskSE(name),
  //General member variables
  fEventCuts(),
  fOutputList(nullptr),
  fEvent(nullptr),
  fMCEvent(nullptr),
  fESDtrackCuts(nullptr),
  fCutMode(100),
  //Toggles
  fIsESD(kTRUE),
  fIsMC(kFALSE),
  // Cut Parameters
  fTriggerMask(AliVEvent::kMB | AliVEvent::kINT7),
  fMinEta(-10),
  fMaxEta(10),
  fMinPt(0.0),
  fMaxPt(50.0),
  //Arrays for Binning
  fBinsMult(nullptr),
  fBinsCent(nullptr),
  fBinsPt(nullptr),
  fBinsEta(nullptr),
  fBinsZv(nullptr),
  fBinsPtReso(nullptr),
  //Event-Histograms
  fHistEventSelection(nullptr),
  fHistEvents(nullptr),
  fHistTracks(nullptr),
  fHistRelPtReso(nullptr),
  fHistMCRelPtReso(nullptr),
  fHistMCMultCorrelMatrix(nullptr),
  fHistMCPtCorrelMatrix(nullptr),
  fHistMCEtaCorrelMatrix(nullptr),
  fHistMCPrimTrue(nullptr),
  fHistMCPrimMeas(nullptr),
  fHistMCSecMeas(nullptr)
{
  // Set default binning
  Double_t binsMultDefault[2] = {0., 10000.};
  Double_t binsCentDefault[2] = {0., 100.};
  Double_t binsPtDefault[53] = {0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,20.0,30.0,40.0,50.0,60.0};
  Double_t binsEtaDefault[19] = {-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
  Double_t binsZvDefault[13] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};

  // binning for relative pT resolution
  const Int_t nBinsPtReso = 300;
  Double_t binsPtReso[nBinsPtReso+1];
  SetFixedBinEdges(binsPtReso, 0., 0.3, nBinsPtReso);
  SetBinsPtReso(nBinsPtReso, binsPtReso);

  SetBinsMult(1,binsMultDefault);
  SetBinsCent(1,binsCentDefault);
  SetBinsPt(52, binsPtDefault);
  SetBinsEta(18,binsEtaDefault);
  SetBinsZv(12,binsZvDefault);

  DefineOutput(1, TList::Class());
}

/***************************************************************************//**
 * Function executed once before the event loop.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::UserCreateOutputObjects(){
  // Create histograms here (function is called once)
  OpenFile(1, "recreate");
  fOutputList = new TList();
  fOutputList->SetOwner();

  // Control histogram to check the effect of event cuts
  fHistEventSelection = new TH1F("fHistEventSelection","fHistEventSelection [all : selected]",2,0.5,2.5);
  fHistEventSelection->GetYaxis()->SetTitle("#it{N}_{events}");
  fHistEventSelection->GetXaxis()->SetBinLabel(1, "all");
  fHistEventSelection->GetXaxis()->SetBinLabel(2, "selected");
  fOutputList->Add(fHistEventSelection);

  fHistEvents = CreateHistogram("fHistEvents", {"mult_meas", "cent"});
  fOutputList->Add(fHistEvents);

  fHistTracks = CreateHistogram("fHistTracks", {"pt_meas", "eta_meas", "mult_meas", "cent"});
  fOutputList->Add(fHistTracks);

  fHistRelPtReso = CreateHistogram("fHistRelPtReso", {"sigmapt", "pt_meas", "cent"});
  fOutputList->Add(fHistRelPtReso);

  if(fIsMC)
  {
    fHistMCRelPtReso = CreateHistogram("fHistMCRelPtReso", {"deltapt", "pt_meas", "cent"});
    fOutputList->Add(fHistMCRelPtReso);

    fHistMCMultCorrelMatrix = CreateHistogram("fHistMCMultCorrelMatrix", {"mult_meas", "mult_true"});
    fOutputList->Add(fHistMCMultCorrelMatrix);

    fHistMCPtCorrelMatrix = CreateHistogram("fHistMCPtCorrelMatrix", {"pt_meas", "pt_true"});
    fOutputList->Add(fHistMCPtCorrelMatrix);

    fHistMCEtaCorrelMatrix = CreateHistogram("fHistMCEtaCorrelMatrix", {"eta_meas", "eta_true"});
    fOutputList->Add(fHistMCEtaCorrelMatrix);

    fHistMCPrimTrue = CreateHistogram("fHistMCPrimTrue", {"pt_true", "eta_true", "mult_true", "cent"});
    fOutputList->Add(fHistMCPrimTrue);

    fHistMCPrimMeas = CreateHistogram("fHistMCPrimMeas", {"pt_true", "eta_true", "mult_true", "cent"});
    fOutputList->Add(fHistMCPrimMeas);

    fHistMCSecMeas = CreateHistogram("fHistMCSecMeas", {"pt_true", "eta_true", "mult_true", "cent"});
    fOutputList->Add(fHistMCSecMeas);
  }


  // override event automatic event selection settings
  //fEventCuts.SetMaxVertexZposition();
  //fEventCuts.SetCentralityRange();
  //fEventCuts.OverrideAutomaticTriggerSelection();

  if(fIsESD) InitESDTrackCuts();

  PostData(1, fOutputList);
}

/// Destructor
AliMultDepSpecAnalysisTask::~AliMultDepSpecAnalysisTask(){
  if(fESDtrackCuts){delete fESDtrackCuts; fESDtrackCuts = nullptr;}
}

/***************************************************************************//**
 * Function executed for each event.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::UserExec(Option_t *){

  fEvent = InputEvent();
  if (!fEvent) {Printf("ERROR: fEvent not available\n"); return;}

  if(fIsMC){
    fMCEvent = MCEvent();
    if (!fMCEvent) {Printf("ERROR: fMCEvent not available\n"); return;}
  }

  fHistEventSelection->Fill(1.0); // all events
  if(!fEventCuts.AcceptEvent(fEvent)) return;
  fHistEventSelection->Fill(2.0); // selected events

  Double_t mult_meas = 0;
  Double_t mult_true = 0;

  Double_t centrality = 50;
  if((fBinsCent->GetSize()-1) > 1) centrality = GetCentrality(fEvent);

  /// ------------------ Count Multiplicities --------------------------------------

  // True Multiplicity mult_true:
  if(fIsMC){
    for(Int_t iGenPart = 0; iGenPart < fMCEvent->GetNumberOfTracks(); iGenPart++) {
      AliMCParticle* mcGenParticle  = (AliMCParticle*)fMCEvent->GetTrack(iGenPart);
      if(!mcGenParticle) {Printf("ERROR: mcGenParticle  not available\n"); continue;}
      if(!AcceptKinematics(mcGenParticle)) continue;
      if(IsChargedPrimary(iGenPart)) mult_true++;
    }
  }

  // Measured Multiplicity mult_meas:
  AliVTrack* track = nullptr;
  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); iTrack++){
    track = fEvent->GetVTrack(iTrack);
    if (!track){Printf("ERROR: Could not receive track %d\n", iTrack); continue;}
    if(!AcceptKinematics(track)) continue;
    if(!AcceptTrackQuality(track)) continue;
    mult_meas++;
  }

  if(fIsMC){
    // Response Matrix
    FillHisto(fHistMCMultCorrelMatrix, {mult_meas, mult_true});
  }

  /// ------------------ Event Histogram ---------------------------------------

  FillHisto(fHistEvents, {mult_meas, centrality});

  ///--------------- Loop over measured Tracks ---------------------------------

  track = NULL;
  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); iTrack++){
    track = fEvent->GetVTrack(iTrack);
    if(!track) {Printf("ERROR: Could not receive track %d\n", iTrack); continue;}
    if(!AcceptKinematics(track)) continue;
    if(!AcceptTrackQuality(track)) continue;


    FillHisto(fHistTracks, {track->Pt(), track->Eta(), mult_meas, centrality});
    // todo make this work for alivtrack..
    FillHisto(fHistRelPtReso, {1./TMath::Abs(dynamic_cast<AliESDtrack*>(track)->GetSigned1Pt())*TMath::Sqrt(dynamic_cast<AliESDtrack*>(track)->GetSigma1Pt2()), track->Pt(), centrality});


    /// Find original particle in MC-Stack
    if(fIsMC){
      Int_t mcLabel = TMath::Abs(track->GetLabel()); // negative label means bad quality track
      AliMCParticle* mcParticle  = (AliMCParticle*)fMCEvent->GetTrack(mcLabel);
      if(!mcParticle) {Printf("ERROR: mcParticle not available\n"); continue;}

      if(!AcceptKinematics(mcParticle)) continue;

      FillHisto(fHistMCRelPtReso, {TMath::Abs(track->Pt() - mcParticle->Pt())/track->Pt(), track->Pt(), centrality});

      if(IsChargedPrimary(mcLabel))
      {
        FillHisto(fHistMCPtCorrelMatrix, {track->Pt(), mcParticle->Pt()});
        FillHisto(fHistMCEtaCorrelMatrix, {track->Eta(), mcParticle->Eta()});
        FillHisto(fHistMCPrimMeas, {mcParticle->Pt(), mcParticle->Eta(), mult_true, centrality});
      }else{
        FillHisto(fHistMCSecMeas, {mcParticle->Pt(), mcParticle->Eta(), mult_true, centrality});
      }
    }
  }

  ///------------------- Loop over Generated Tracks (True MC)------------------------------
  if (fIsMC){

    for(Int_t iGenPart = 0; iGenPart < fMCEvent->GetNumberOfTracks(); iGenPart++) {
      AliMCParticle* mcGenParticle  = (AliMCParticle*)fMCEvent->GetTrack(iGenPart);
      if(!mcGenParticle) {Printf("ERROR: mcGenParticle  not available\n"); continue;}

      if(!AcceptKinematics(mcGenParticle)) continue;

      if(IsChargedPrimary(iGenPart)){
        FillHisto(fHistMCPrimTrue, {mcGenParticle->Pt(), mcGenParticle->Eta(), mult_true, centrality});
      }
    }
  }
  PostData(1, fOutputList);
}

/***************************************************************************//**
 * Function executed after all events were processed.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::Terminate(Option_t*)
{

}

/***************************************************************************//**
 * Function to select primary charged particles.
 ******************************************************************************/
Bool_t AliMultDepSpecAnalysisTask::IsChargedPrimary(Int_t mcLabel)
{
  if(!fMCEvent->IsPhysicalPrimary(mcLabel)) return kFALSE;
  AliMCParticle* mcParticle  = (AliMCParticle*)fMCEvent->GetTrack(mcLabel);
  if(!mcParticle) {Printf("ERROR: mcGenParticle  not available\n"); return kFALSE;}
  if(!(TMath::Abs(mcParticle->Charge()) > 0.01)) return kFALSE;
  return kTRUE;
}

/***************************************************************************//**
 * Function to select if the track or mc particle is in defined kinematic range.
 ******************************************************************************/
Bool_t AliMultDepSpecAnalysisTask::AcceptKinematics(AliVParticle* particle)
{
  if(!particle) return kFALSE;

  Double_t eta = particle->Eta();
  Double_t pt  = particle->Pt();

  if(eta <= fMinEta + PRECISION)  return kFALSE;
  if(eta >= fMaxEta - PRECISION)  return kFALSE;
  if(pt  <= fMinPt  + PRECISION)  return kFALSE;
  if(pt  >= fMaxPt  - PRECISION)  return kFALSE;
  return kTRUE;
}

/***************************************************************************//**
 * Function to select tracks with required quality.
 ******************************************************************************/
Bool_t AliMultDepSpecAnalysisTask::AcceptTrackQuality(AliVTrack* track){
  if(fIsESD){
    AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (track);
    if(!fESDtrackCuts->AcceptTrack(esdTrack)) return kFALSE;
  }
  return kTRUE;
}

/***************************************************************************//**
 * Function to obtain V0M centrality.
 ******************************************************************************/
Double_t AliMultDepSpecAnalysisTask::GetCentrality(AliVEvent* event)
{
  Double_t centrality = -1;
  AliMultSelection* multSelection = (AliMultSelection*) fEvent->FindListObject("MultSelection");
  if(!multSelection){AliInfo("ERROR: No MultSelection found!"); return 999;}
  centrality = multSelection->GetMultiplicityPercentile("V0M"/*, lEmbedEventSelection = kFALSE*/);
  if(centrality > 100) {AliInfo("ERROR: Centrality determination does not work proprely!"); return 999;}
  return centrality;
}

/***************************************************************************//**
 * Function to initialize the ESD track cuts object.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::InitESDTrackCuts(){

  if(fESDtrackCuts) delete fESDtrackCuts;
  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
  if(!fESDtrackCuts) {printf("ERROR: fESDtrackCuts not available\n"); return;}

  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(4.0);
  fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.4);
  fESDtrackCuts->SetRequireITSRefit(kTRUE);
  fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  fESDtrackCuts->SetMaxChi2PerClusterITS(36.);
  fESDtrackCuts->SetDCAToVertex2D(kFALSE);
  fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
  fESDtrackCuts->SetMaxDCAToVertexZ(2.0);
  fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); // 7*(0.0026+0.0050/pt^1.01)
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36.); // tpcc cut
  fESDtrackCuts->SetCutGeoNcrNcl(3,130,1.5,0.85,0.7); // Geometrical-Length Cut

  fESDtrackCuts->SetEtaRange(-0.8, 0.8);
  fESDtrackCuts->SetPtRange(0.15, 50.0);

  // cut-variations:
  if(fCutMode == 101) {fESDtrackCuts->SetMaxChi2PerClusterITS(25.);}
  if(fCutMode == 102) {fESDtrackCuts->SetMaxChi2PerClusterITS(49.);}

  if(fCutMode == 103) {fESDtrackCuts->SetMaxChi2PerClusterTPC(3.0); }
  if(fCutMode == 104) {fESDtrackCuts->SetMaxChi2PerClusterTPC(5.0); }

  if(fCutMode == 105) {fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);}
  if(fCutMode == 106) {fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);}

  if(fCutMode == 107) {fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.2);}
  if(fCutMode == 108) {fESDtrackCuts->SetMaxFractionSharedTPCClusters(1.0);}

  if(fCutMode == 109) {fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(25.);}
  if(fCutMode == 110) {fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(49.);}

  if(fCutMode == 111) {fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0104+0.0200/pt^1.01");}
  if(fCutMode == 112) {fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0260+0.0500/pt^1.01");}

  if(fCutMode == 113) {fESDtrackCuts->SetMaxDCAToVertexZ(1.0);}
  if(fCutMode == 114) {fESDtrackCuts->SetMaxDCAToVertexZ(5.0);}

  if(fCutMode == 115) {fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);}

  if(fCutMode == 116) {fESDtrackCuts->SetCutGeoNcrNcl(3,120,1.5,0.85,0.7);}
  if(fCutMode == 117) {fESDtrackCuts->SetCutGeoNcrNcl(3,140,1.5,0.85,0.7);}

  if(fCutMode == 118) {fESDtrackCuts->SetCutGeoNcrNcl(4,130,1.5,0.85,0.7);}
  if(fCutMode == 119) {fESDtrackCuts->SetCutGeoNcrNcl(2,130,1.5,0.85,0.7);}
}


/***************************************************************************//**
 * Function to get array of equidistant bin edges between lower and upper edge.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::SetFixedBinEdges(Double_t* array, Double_t lowerEdge, Double_t upperEdge, Int_t nBins){
  for(Int_t i = 0; i <= nBins; i++){
    array[i] = lowerEdge + i*(upperEdge - lowerEdge)/nBins;
  }
}

/***************************************************************************//**
 * Function to create THnSparseF histogram with the specified axes.
 ******************************************************************************/
THnSparseF* AliMultDepSpecAnalysisTask::CreateHistogram(string name, vector<string> axes){
  Int_t nAxes = axes.size();
  if(nAxes > MAX_HISTO_DIM) return nullptr;

  Int_t nBins[MAX_HISTO_DIM] = {0};
  Double_t lowerBounds[MAX_HISTO_DIM] = {0.0};
  Double_t upperBounds[MAX_HISTO_DIM] = {0.0};

  string title = name + " [";
  // first figure out number of bins and dimensions
  for(Int_t i = 0; i < nAxes; i++){
    TArrayD* binEdges = GetBinEdges(axes[i]);
    nBins[i] = binEdges->GetSize()-1;
    lowerBounds[i] = binEdges->GetAt(0);
    upperBounds[i] = binEdges->GetAt(binEdges->GetSize()-1);
    title += axes[i];
    if(i < nAxes-1) title += " : "; else title += "]";
  }
  // create histogram
  THnSparseF* histogram = new THnSparseF(name.c_str(), title.c_str(), nAxes, nBins, lowerBounds, upperBounds);

  // set histogram axes
  for(Int_t i = 0; i < nAxes; i++){
    TArrayD* binEdges = GetBinEdges(axes[i]);
    histogram->SetBinEdges(i, binEdges->GetArray());
    histogram->GetAxis(i)->SetTitle(GetAxisTitle(axes[i]).c_str());
  }
  histogram->Sumw2();
  return histogram;
}

/***************************************************************************//**
 * Function to obtain the correct binning for the respective axis.
 ******************************************************************************/
TArrayD* AliMultDepSpecAnalysisTask::GetBinEdges(string& axisName){
  if(axisName.find("sigmapt") != string::npos) return fBinsPtReso;
  else if(axisName.find("deltapt") != string::npos) return fBinsPtReso;
  else if(axisName.find("pt") != string::npos) return fBinsPt;
  else if(axisName.find("eta") != string::npos) return fBinsEta;
  else if(axisName.find("mult") != string::npos) return fBinsMult;
  else if(axisName.find("cent") != string::npos) return fBinsCent;
  else if(axisName.find("zv") != string::npos) return fBinsZv;
  else return nullptr;
}

/***************************************************************************//**
 * Function to get the correct title for each histogram axis.
 ******************************************************************************/
string AliMultDepSpecAnalysisTask::GetAxisTitle(string& axisName){
  if(axisName == "pt") return "#it{p}_{T} (GeV/#it{c})";
  else if(axisName == "deltapt") return "#Delta(#it{p}_{T}) / #it{p}^{ meas}_{T}";
  else if(axisName == "mult") return "Multiplicity";
  else if(axisName == "cent") return "Centrality (%)";
  else if(axisName == "eta_meas") return "#eta^{ meas}";
  else if(axisName == "eta_true") return "#eta^{ true}";
  else if(axisName == "pt_meas") return "#it{p}^{ meas}_{T} (GeV/#it{c})";
  else if(axisName == "pt_true") return "#it{p}^{ true}_{T} (GeV/#it{c})";
  else if(axisName == "mult_meas") return "#it{N}^{ meas}_{ch}";
  else if(axisName == "mult_true") return "#it{N}^{ true}_{ch}";
  else if(axisName == "sigmapt") return "#sigma(#it{p}^{ meas}_{T}) / #it{p}^{ meas}_{T}";
  else return "dummyTitle";
}

/***************************************************************************//**
 * Function to fill a histogram.
 ******************************************************************************/
void AliMultDepSpecAnalysisTask::FillHisto(THnSparseF* histo, array<Double_t, MAX_HISTO_DIM> values){
  histo->Fill(values.data());
}




/***************************************************************************//**
 * Function to add this task to a train.
 ******************************************************************************/
//void AliMultDepSpecAnalysisTask::AddTaskMultDepSpec(){


//}
