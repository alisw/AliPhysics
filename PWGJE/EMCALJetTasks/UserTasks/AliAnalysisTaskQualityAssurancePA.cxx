#ifndef ALIANALYSISTASKSE_H
#include <Riostream.h>
#include <TROOT.h>
#include <TString.h>
#include <TFile.h>
#include <TCint.h>
#include <TChain.h>
#include <TTree.h>
#include <TKey.h>
#include <TProfile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TMath.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TH1.h>
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#endif

#include <TRandom3.h>
#include "AliGenPythiaEventHeader.h"
#include "AliMCEvent.h"
#include "AliLog.h"
#include <AliEmcalJet.h>
#include <AliRhoParameter.h>
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliAnalysisUtils.h"

#include "AliAnalysisTaskQualityAssurancePA.h"


//TODO: FillHistogram can be done better with virtual TH1(?)
ClassImp(AliAnalysisTaskQualityAssurancePA)

// ######################################################################################## DEFINE HISTOGRAMS
void AliAnalysisTaskQualityAssurancePA::Init()
{
  #ifdef DEBUGMODE
    AliInfo("Creating histograms.");
  #endif
  *fRunNumbers = *fRunNumbers + " ALL";
   
  TObjArray* runNumberArr = fRunNumbers->Tokenize(" ");

  for(Int_t i=0; i<runNumberArr->GetEntries();i++)
  {
    const char* tmpRunNum = (static_cast<TObjString*>(runNumberArr->At(i)))->String().Data();
    
    // NOTE: Track & Cluster & QA histograms
    if (fAnalyzeQA)
    {
      AddHistogram1D<TH1D>(tmpRunNum, "hNumberEvents", "Number of events (0 = before, 1 = after vertex cuts)", "", 2, 0, 2, "#Delta z(cm)","N^{Events}/cut");
      AddHistogram1D<TH1D>(tmpRunNum, "hVertexX", "X distribution of the vertex", "", 10000, -1., 1., "#Delta x(cm)","dN^{Events}/dx");
      AddHistogram1D<TH1D>(tmpRunNum, "hVertexY", "Y distribution of the vertex", "", 10000, -1., 1., "#Delta y(cm)","dN^{Events}/dy");
      AddHistogram2D<TH2D>(tmpRunNum, "hVertexXY", "XY distribution of the vertex", "COLZ", 1000, -1., 1., 1000, -1., 1.,"#Delta x(cm)", "#Delta y(cm)","dN^{Events}/dxdy");
      AddHistogram1D<TH1D>(tmpRunNum, "hVertexZ", "Z distribution of the vertex", "", 400, -40., 40., "#Delta z(cm)","dN^{Events}/dz");
      AddHistogram1D<TH1D>(tmpRunNum, "hVertexR", "R distribution of the vertex", "", 100, 0., 1., "#Delta r(cm)","dN^{Events}/dr");
      AddHistogram1D<TH1D>(tmpRunNum, "hCentralityV0M", "Centrality distribution V0M", "", 100, 0., 100., "Centrality","dN^{Events}");
      AddHistogram1D<TH1D>(tmpRunNum, "hCentralityV0A", "Centrality distribution V0A", "", 100, 0., 100., "Centrality","dN^{Events}");
      AddHistogram1D<TH1D>(tmpRunNum, "hCentralityV0C", "Centrality distribution V0C", "", 100, 0., 100., "Centrality","dN^{Events}");

      AddHistogram2D<TH2D>(tmpRunNum, "hTrackCountAcc", "Number of tracks in acceptance vs. centrality", "LEGO2", 750, 0., 750., 100, 0, 100, "N tracks","Centrality", "dN^{Events}/dN^{Tracks}");
      AddHistogram1D<TH1D>(tmpRunNum, "hTrackPt", "Tracks p_{T} distribution", "", 20000, 0., 400., "p_{T} (GeV/c)","dN^{Tracks}/dp_{T}");
      AddHistogram1D<TH1D>(tmpRunNum, "hTrackPtNegEta", "Tracks p_{T} distribution (negative #eta)", "", 2000, 0., 400., "p_{T} (GeV/c)","dN^{Tracks}/dp_{T}");
      AddHistogram1D<TH1D>(tmpRunNum, "hTrackPtPosEta", "Tracks p_{T} distribution (positive #eta)", "", 2000, 0., 400., "p_{T} (GeV/c)","dN^{Tracks}/dp_{T}");      
      AddHistogram1D<TH1D>(tmpRunNum, "hTrackCharge", "Charge", "", 11, -5, 5, "Charge (e)","dN^{Tracks}/dq");
      AddHistogram1D<TH1D>(tmpRunNum, "hTrackPhi", "Track #phi distribution", "", 360, 0, TMath::TwoPi(), "#phi","dN^{Tracks}/d#phi");
      AddHistogram2D<TH2D>(tmpRunNum, "hTrackPhiEta", "Track angular distribution", "LEGO2", 100, 0., 2*TMath::Pi(),100, -2.5, 2.5, "#phi","#eta","dN^{Tracks}/(d#phi d#eta)");

      AddHistogram2D<TH2D>(tmpRunNum, "hTrackPhiPtCut", "Track #phi distribution for different pT cuts", "LEGO2", 360, 0, TMath::TwoPi(), 20, 0, 20, "#phi", "p_{T} lower cut", "dN^{Tracks}/d#phi dp_{T}");      
      AddHistogram2D<TH2D>(tmpRunNum, "hTrackPhiLabel", "Track #phi distribution in different labels", "LEGO2", 360, 0, TMath::TwoPi(), 3, 0, 3, "#phi", "Label", "dN^{Tracks}/d#phi");
      AddHistogram1D<TH1D>(tmpRunNum, "hTrackEta", "Track #eta distribution", "", 180, -fTrackEtaWindow, +fTrackEtaWindow, "#eta","dN^{Tracks}/d#eta");
    }

    // NOTE: Pythia histograms
    if (fAnalyzePythia)
    {
      AddHistogram1D<TH1D>(tmpRunNum, "hPythiaPtHard", "Pythia p_{T} hard distribution", "", 2000, 0, 400, "p_{T} hard","dN^{Events}/dp_{T,hard}");
      AddHistogram1D<TProfile>(tmpRunNum, "hPythiaXSection", "Pythia cross section distribution", "", fNumPtHardBins, 0, fNumPtHardBins, "p_{T} hard bin","dN^{Events}/dp_{T,hard}");
      AddHistogram1D<TH1D>(tmpRunNum, "hPythiaNTrials", "Pythia trials (no correction for manual cuts)", "", fNumPtHardBins, 0, fNumPtHardBins, "p_{T} hard bin", "Trials");
    }

    // NOTE: Jet histograms
    if (fAnalyzeJets)
    {
      // ######## Jet spectra
      AddHistogram1D<TH1D>(tmpRunNum, "hJetPt", "Jets p_{T} distribution", "", 1000, 0., 200., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
      AddHistogram1D<TH1D>(tmpRunNum, "hJetArea", "Jets area distribution", "", 200, 0., 2., "Area","dN^{Jets}/dA");
      AddHistogram2D<TH2D>(tmpRunNum, "hJetAreaVsPt", "Jets area vs. p_{T} distribution", "COLZ", 200, 0., 2., 400, 0., 200., "Area", "p_{T}", "dN^{Jets}/dA dp_{T}");

      AddHistogram2D<TH2D>(tmpRunNum, "hJetPhiEta", "Jets angular distribution", "LEGO2", 360, 0., 2*TMath::Pi(),100, -0.6, 0.6, "#phi","#eta","dN^{Jets}/(d#phi d#eta)");
      AddHistogram2D<TH2D>(tmpRunNum, "hJetPtVsConstituentCount", "Jets number of constituents vs. jet p_{T}", "COLZ", 800, 0., 400., 100, 0., 100., "p_{T}","N^{Tracks}","dN^{Jets}/(dp_{T} dN^{tracks})");
      AddHistogram1D<TH1D>(tmpRunNum, "hJetCountAll", "Number of Jets", "", 200, 0., 200., "N jets","dN^{Events}/dN^{Jets}");
      AddHistogram1D<TH1D>(tmpRunNum, "hJetCountAccepted", "Number of accepted Jets", "", 200, 0., 200., "N jets","dN^{Events}/dN^{Jets}");

    }
  }
  // register Histograms
  for (Int_t i = 0; i < fHistCount; i++)
  {
    fOutputList->Add(fHistList->At(i));
  }
  
  PostData(1,fOutputList); // important for merging

}


//________________________________________________________________________
AliAnalysisTaskQualityAssurancePA::AliAnalysisTaskQualityAssurancePA() : AliAnalysisTaskSE("AliAnalysisTaskQualityAssurancePA"), fOutputList(0), fAnalyzeQA(1), fAnalyzeJets(1), fAnalyzePythia(0), fHasTracks(0), fHasClusters(0), fHasJets(0), fIsMC(0), fJetArray(0), fTrackArray(0), fClusterArray(0), fJetArrayName(0), fTrackArrayName(0), fClusterArrayName(0), fRunNumbers(0), fNumPtHardBins(11), fSignalJetRadius(0.4), fNumberExcludedJets(2), fSignalJetEtaWindow(0.5), fTrackEtaWindow(0.9), fClusterEtaWindow(0.7), fVertexWindow(10.0), fVertexMaxR(1.0), fMinTrackPt(0.150), fMinClusterPt(0.300), fMinJetPt(1.0), fMinJetArea(0.4), fFirstLeadingJet(0), fSecondLeadingJet(0), fNumberSignalJets(0), fCrossSection(0.0), fTrials(0.0), fRandom(0), fHelperClass(0), fInitialized(0), fTaskInstanceCounter(0), fHistList(0), fHistCount(0)
{
  // default constructor
}


//________________________________________________________________________
AliAnalysisTaskQualityAssurancePA::AliAnalysisTaskQualityAssurancePA(const char *name, const char* trackArrayName, const char* clusterArrayName, const char* jetArrayName) : AliAnalysisTaskSE(name), fOutputList(0), fAnalyzeQA(1), fAnalyzeJets(1), fAnalyzePythia(0), fHasTracks(0), fHasClusters(0), fHasJets(0), fIsMC(0), fJetArray(0), fTrackArray(0), fClusterArray(0), fJetArrayName(0), fTrackArrayName(0), fClusterArrayName(0), fRunNumbers(0), fNumPtHardBins(11), fSignalJetRadius(0.4), fNumberExcludedJets(2), fSignalJetEtaWindow(0.5), fTrackEtaWindow(0.9), fClusterEtaWindow(0.7), fVertexWindow(10.0), fVertexMaxR(1.0), fMinTrackPt(0.150), fMinClusterPt(0.300), fMinJetPt(1.0), fMinJetArea(0.4), fFirstLeadingJet(0), fSecondLeadingJet(0), fNumberSignalJets(0), fCrossSection(0.0), fTrials(0.0), fRandom(0), fHelperClass(0), fInitialized(0), fTaskInstanceCounter(0), fHistList(0), fHistCount(0)
{
  #ifdef DEBUGMODE
    AliInfo("Calling constructor.");
  #endif

  // Constructor
  // Define input and output slots here (never in the dummy constructor)
  // Input slot #0 works with a TChain - it is connected to the default input container
  // Output slot #1 writes into a TH1 container
  // Constructor

  // Every instance of this task gets his own number
  static Int_t instance = 0;
  fTaskInstanceCounter = instance;
  instance++;

  fTrackArrayName = new TString(trackArrayName);
  fClusterArrayName = new TString(clusterArrayName);
  fRunNumbers = new TString("");
  if (strcmp(fTrackArrayName->Data(),"") == 0)
    fAnalyzeQA = kFALSE;
  else
  {
    fAnalyzeQA = kTRUE;
    if (fTrackArrayName->Contains("MCParticles")) //TODO: Hardcoded for now
      fIsMC = kTRUE;
  }

  fJetArrayName = new TString(jetArrayName);
  if (strcmp(fJetArrayName->Data(),"") == 0)
    fAnalyzeJets = kFALSE;
  else
    fAnalyzeJets = kTRUE;

  DefineOutput(1, TList::Class());
 
  fHistList = new TList();

  #ifdef DEBUGMODE
    AliInfo("Constructor done.");
  #endif
  
}

//________________________________________________________________________
inline Double_t AliAnalysisTaskQualityAssurancePA::GetPtHard()
{
  Double_t tmpPtHard = -1.0;

  if (!MCEvent())
    AliError("MCEvent not accessible although demanded!");
  else
  {
    AliGenPythiaEventHeader* pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
    if (!pythiaHeader)
      AliError("Pythia Header not accessible!");
    else
      tmpPtHard = pythiaHeader->GetPtHard();
  }
  return tmpPtHard;
}

//________________________________________________________________________
inline Int_t AliAnalysisTaskQualityAssurancePA::GetPtHardBin()
{
  // ########## PT HARD BIN EDGES
  const Int_t localkPtHardLowerEdges[] = { 0, 5,11,21,36,57, 84,117,152,191,234};
  const Int_t localkPtHardHigherEdges[] = { 5,11,21,36,57,84,117,152,191,234,1000000};

  Int_t tmpPtHardBin = 0;
  Double_t tmpPtHard = GetPtHard();
 
  for (tmpPtHardBin = 0; tmpPtHardBin <= fNumPtHardBins; tmpPtHardBin++)
    if (tmpPtHard >= localkPtHardLowerEdges[tmpPtHardBin] && tmpPtHard < localkPtHardHigherEdges[tmpPtHardBin])
      break;

  return tmpPtHardBin;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskQualityAssurancePA::IsTrackInAcceptance(AliVParticle* track)
{
  if (track != 0)
    if (TMath::Abs(track->Eta()) <= fTrackEtaWindow)
      if (track->Pt() >= fMinTrackPt)
        return kTRUE;

  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskQualityAssurancePA::IsClusterInAcceptance(AliVCluster* cluster)
{
  if (cluster != 0)
 //   if (TMath::Abs(cluster->Eta()) <= fClusterEtaWindow)
//      if (cluster->Phi() <= 187.0/360.0 * TMath::TwoPi());
//        if (cluster->Phi() >= 80.0/360.0 * TMath::TwoPi());
          if (cluster->E() >= fMinClusterPt)
            return kTRUE;

  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskQualityAssurancePA::IsSignalJetInAcceptance(AliEmcalJet *jet)
{   
  if (jet != 0)
    if (TMath::Abs(jet->Eta()) <= fSignalJetEtaWindow)
      if (jet->Pt() >= fMinJetPt)
        if (jet->Area() >= fMinJetArea)
          return kTRUE;

  return kFALSE;
}

//________________________________________________________________________
void AliAnalysisTaskQualityAssurancePA::ExecOnce()
{
  #ifdef DEBUGMODE
    AliInfo("Starting ExecOnce.");
  #endif
  fInitialized = kTRUE;

  // Check for track array
  if (strcmp(fTrackArrayName->Data(), "") != 0)
  {
    fTrackArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTrackArrayName->Data()));
    fHasTracks = kTRUE;
    if (!fTrackArray) 
    {
      AliInfo(Form("%s: Could not retrieve tracks %s! This is OK, if tracks are not demanded.", GetName(), fTrackArrayName->Data())); 
      fHasTracks = kFALSE;
    } 
    else
    {
      TClass *cl = fTrackArray->GetClass();
      if (!cl->GetBaseClass("AliVParticle"))
      {
      	AliError(Form("%s: Collection %s does not contain AliVParticle objects!", GetName(), fTrackArrayName->Data())); 
      	fTrackArray = 0;
        fHasTracks = kFALSE;
      }
    }
  }
  // Check for cluster array
  if (strcmp(fClusterArrayName->Data(), "") != 0)
  {
    fClusterArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fClusterArrayName->Data()));
    fHasClusters = kTRUE;
    if (!fClusterArray) 
    {
      AliInfo(Form("%s: Could not retrieve clusters %s! This is OK, if clusters are not demanded.", GetName(), fClusterArrayName->Data())); 
      fHasClusters = kFALSE;
    } 
    else
    {
      TClass *cl = fClusterArray->GetClass();
      if (!cl->GetBaseClass("AliVCluster"))
      {
      	AliError(Form("%s: Collection %s does not contain AliVCluster objects!", GetName(), fClusterArrayName->Data())); 
      	fClusterArray = 0;
        fHasClusters = kFALSE;
      }
    }
  }

  // Check for jet array
  if (strcmp(fJetArrayName->Data(), "") != 0)
  {
    fJetArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetArrayName->Data()));
    fHasJets = kTRUE;

    if (!fJetArray) 
    {
      AliInfo(Form("%s: Could not retrieve jets %s! This is OK, if jets are not demanded.", GetName(), fJetArrayName->Data())); 
      fHasJets = kFALSE;
    } 
    else
    {
      if (!fJetArray->GetClass()->GetBaseClass("AliEmcalJet")) 
      {
        AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fJetArrayName->Data())); 
        fJetArray = 0;
        fHasJets = kFALSE;
      }
    }
  }

  // Look, if initialization is OK
  if (!fHasTracks && fAnalyzeQA)
  {
    AliError(Form("%s: Tracks NOT successfully casted although demanded! Deactivating QA",GetName()));
    fAnalyzeQA = kFALSE;
  }
  if (!fHasJets && fAnalyzeJets)
  {
    AliError(Form("%s: Jets NOT successfully casted although demanded!  Deactivating jet",GetName()));
    fAnalyzeJets = kFALSE;
  }

  // Initialize helper class (for vertex selection)
  fHelperClass = new AliAnalysisUtils();

  // Histogram init
  Init();

  #ifdef DEBUGMODE
    AliInfo("ExecOnce done.");
  #endif

}

//________________________________________________________________________
void AliAnalysisTaskQualityAssurancePA::GetSignalJets()
{
  // Reset vars
  fFirstLeadingJet = NULL;
  fSecondLeadingJet = NULL;
  fNumberSignalJets = 0;

  Float_t maxJetPts[] = {0, 0};
  Int_t jetIDArray[]   = {-1, -1};
  Int_t jetCount = fJetArray->GetEntries();

  // Go through all jets and save signal jets and the two leading ones
  for (Int_t i = 0; i < jetCount; i++)
  {
    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJetArray->At(i));
    if (!jet)
    {
      AliError(Form("%s: Could not receive jet %d", GetName(), i));
      continue;
    }

    if (!IsSignalJetInAcceptance(jet)) continue;
    
    if (jet->Pt() > maxJetPts[0]) 
    {
      maxJetPts[1] = maxJetPts[0];
      jetIDArray[1] = jetIDArray[0];
      maxJetPts[0] = jet->Pt();
      jetIDArray[0] = i;
    }
    else if (jet->Pt() > maxJetPts[1]) 
    {
      maxJetPts[1] = jet->Pt();
      jetIDArray[1] = i;
    }
    fSignalJets[fNumberSignalJets] = jet;
    fNumberSignalJets++;
  }
  
  if (fNumberSignalJets > 0)
    fFirstLeadingJet  = static_cast<AliEmcalJet*>(fJetArray->At(jetIDArray[0]));
  if (fNumberSignalJets > 1)
    fSecondLeadingJet = static_cast<AliEmcalJet*>(fJetArray->At(jetIDArray[1]));

}

//________________________________________________________________________
Int_t AliAnalysisTaskQualityAssurancePA::GetLeadingJets(TClonesArray* jetArray, Int_t* jetIDArray)
{
// Writes first two leading jets into already registered array jetIDArray

  if (!jetArray)
  {
    AliError("Could not get the jet array to get leading jets from it!");
    return 0;
  }

  Float_t maxJetPts[] = {0, 0};
  jetIDArray[0] = -1;
  jetIDArray[1] = -1;

  Int_t jetCount = jetArray->GetEntries();
  Int_t jetCountAccepted = 0;

  for (Int_t i = 0; i < jetCount; i++)
  {
    AliEmcalJet* jet = static_cast<AliEmcalJet*>(jetArray->At(i));
    if (!jet) 
    {
      AliError(Form("%s: Could not receive jet %d", GetName(), i));
      continue;
    }

    if (!IsSignalJetInAcceptance(jet)) continue;

    if (jet->Pt() > maxJetPts[0]) 
    {
      maxJetPts[1] = maxJetPts[0];
      jetIDArray[1] = jetIDArray[0];
      maxJetPts[0] = jet->Pt();
      jetIDArray[0] = i;
    }
    else if (jet->Pt() > maxJetPts[1]) 
    {
      maxJetPts[1] = jet->Pt();
      jetIDArray[1] = i;
    }
    jetCountAccepted++;
  }
  return jetCountAccepted;
}


//________________________________________________________________________
void AliAnalysisTaskQualityAssurancePA::Calculate(AliVEvent* event)
{
  #ifdef DEBUGMODE
    AliInfo("Starting Calculate().");
  #endif
  ////////////////////// NOTE: initialization & casting

  if (!event) {
    AliError("??? Event pointer == 0 ???");
    return;
  }
 
  if (!fInitialized)
    ExecOnce(); // Get tracks, jets from arrays if not already given + Init Histos
  
  // Get run number
  TString tmpRunNum("");
  tmpRunNum += event->GetRunNumber();
  // Additional cuts
  FillHistogram(tmpRunNum.Data(), "hNumberEvents", 0.5); // number of events before manual cuts
  
  if(!fHelperClass->IsVertexSelected2013pA(event))
    return;

  FillHistogram(tmpRunNum.Data(), "hNumberEvents", 1.5); // number of events after manual cuts

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Init done.");
  #endif

  ////////////////////// NOTE: Get Centrality, (Leading)Signal jets

  // Get centrality (V0A)
  AliCentrality* tmpCentrality = NULL;
  tmpCentrality = event->GetCentrality();
  Double_t centralityPercentileV0A = 0.0;
  Double_t centralityPercentileV0C = 0.0;
  Double_t centralityPercentileV0M = 0.0;
  if (tmpCentrality != NULL)
  {
    centralityPercentileV0A = tmpCentrality->GetCentralityPercentile("V0A");
    centralityPercentileV0C = tmpCentrality->GetCentralityPercentile("V0C");
    centralityPercentileV0M = tmpCentrality->GetCentralityPercentile("V0M");
  }
  // Get jets
  if (fAnalyzeJets)
    GetSignalJets();

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Centrality&SignalJets done.");
  #endif
  ////////////////////// NOTE: Pythia histograms
  if(fAnalyzePythia)
  {
    FillHistogram(tmpRunNum.Data(), "hPythiaPtHard", GetPtHard());
    FillHistogram(tmpRunNum.Data(), "hPythiaNTrials", GetPtHardBin()+0.1, fTrials);
    FillHistogram(tmpRunNum.Data(), "hPythiaXSection", GetPtHardBin()+0.1, fCrossSection);

    #ifdef DEBUGMODE
      AliInfo("Calculate()::Pythia done.");
    #endif
  }

  ////////////////////// NOTE: Track & QA histograms

  if (fAnalyzeQA)
  {
    FillHistogram(tmpRunNum.Data(), "hVertexX",event->GetPrimaryVertex()->GetX());
    FillHistogram(tmpRunNum.Data(), "hVertexY",event->GetPrimaryVertex()->GetY());
    FillHistogram(tmpRunNum.Data(), "hVertexXY",event->GetPrimaryVertex()->GetX(), event->GetPrimaryVertex()->GetY());
    FillHistogram(tmpRunNum.Data(), "hVertexZ",event->GetPrimaryVertex()->GetZ());
    FillHistogram(tmpRunNum.Data(), "hVertexR",TMath::Sqrt(event->GetPrimaryVertex()->GetX()*event->GetPrimaryVertex()->GetX() + event->GetPrimaryVertex()->GetY()*event->GetPrimaryVertex()->GetY()));
    FillHistogram(tmpRunNum.Data(), "hCentralityV0M",centralityPercentileV0M);
    FillHistogram(tmpRunNum.Data(), "hCentralityV0A",centralityPercentileV0A);
    FillHistogram(tmpRunNum.Data(), "hCentralityV0C",centralityPercentileV0C);

    Int_t trackCountAcc = 0;
    Int_t nTracks = fTrackArray->GetEntries();
    for (Int_t i = 0; i < nTracks; i++)
    {
      AliVTrack* track = static_cast<AliVTrack*>(fTrackArray->At(i));
      if (IsTrackInAcceptance(track))
      {
        FillHistogram(tmpRunNum.Data(), "hTrackPhiEta", track->Phi(),track->Eta(), 1);
        FillHistogram(tmpRunNum.Data(), "hTrackPt", track->Pt());
        if(track->Eta() >= 0)
          FillHistogram(tmpRunNum.Data(), "hTrackPtPosEta", track->Pt());
        else
          FillHistogram(tmpRunNum.Data(), "hTrackPtNegEta", track->Pt());        
                
        FillHistogram(tmpRunNum.Data(), "hTrackEta", track->Eta());
        FillHistogram(tmpRunNum.Data(), "hTrackPhi", track->Phi());
        FillHistogram(tmpRunNum.Data(), "hTrackPhiLabel", track->Phi(), track->GetLabel());
        for(Int_t j=0;j<20;j++)
          if(track->Pt() > j)
            FillHistogram(tmpRunNum.Data(), "hTrackPhiPtCut", track->Phi(), track->Pt());

        FillHistogram(tmpRunNum.Data(), "hTrackCharge", track->Charge());
        trackCountAcc++;
      }
    }
    FillHistogram(tmpRunNum.Data(), "hTrackCountAcc", trackCountAcc, centralityPercentileV0M);

  }
  #ifdef DEBUGMODE
    AliInfo("Calculate()::QA done.");
  #endif

  ////////////////////// NOTE: Jet analysis and calculations

  if (fAnalyzeJets)
  {
    FillHistogram(tmpRunNum.Data(), "hJetCountAll", fJetArray->GetEntries());
    FillHistogram(tmpRunNum.Data(), "hJetCountAccepted", fNumberSignalJets);
    // SIGNAL JET ANALYSIS
    for (Int_t i = 0; i<fNumberSignalJets; i++)
    {
      AliEmcalJet* tmpJet = fSignalJets[i];

      FillHistogram(tmpRunNum.Data(), "hJetArea", tmpJet->Area());
      FillHistogram(tmpRunNum.Data(), "hJetAreaVsPt", tmpJet->Area(), tmpJet->Pt());
      FillHistogram(tmpRunNum.Data(), "hJetPt", tmpJet->Pt());
      FillHistogram(tmpRunNum.Data(), "hJetPtVsConstituentCount", tmpJet->Pt(),tmpJet->GetNumberOfTracks());
      FillHistogram(tmpRunNum.Data(), "hJetPhiEta", tmpJet->Phi(),tmpJet->Eta());


    }
  } //endif AnalyzeJets

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Jets done.");
  #endif
} 

//________________________________________________________________________
Bool_t AliAnalysisTaskQualityAssurancePA::Notify()
{
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 
  #ifdef DEBUGMODE
    AliInfo("Notify started.");
  #endif

  if(fAnalyzePythia)
  {
    TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
    TFile *currFile = tree->GetCurrentFile();

    TString file(currFile->GetName());

    if(file.Contains("root_archive.zip#")){
      Ssiz_t pos1 = file.Index("root_archive",12,TString::kExact);
      Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
      file.Replace(pos+1,20,"");
    }
    else {
      // not an archive take the basename....
      file.ReplaceAll(gSystem->BaseName(file.Data()),"");
    }
   
    TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
    if(!fxsec){
      // next trial fetch the histgram file
      fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
      if(!fxsec){
          // not a severe condition but inciate that we have no information
        return kFALSE;
      }
      else{
        // find the tlist we want to be independtent of the name so use the Tkey
        TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0); 
        if(!key){
          fxsec->Close();
          return kFALSE;
        }
        TList *list = dynamic_cast<TList*>(key->ReadObj());
        if(!list){
          fxsec->Close();
          return kFALSE;
        }
        fCrossSection = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
        fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
        fxsec->Close();
      }
    } // no tree pyxsec.root
    else {
      TTree *xtree = (TTree*)fxsec->Get("Xsection");
      if(!xtree){
        fxsec->Close();
        return kFALSE;
      }
      UInt_t   ntrials  = 0;
      Double_t  xsection  = 0;
      xtree->SetBranchAddress("xsection",&xsection);
      xtree->SetBranchAddress("ntrials",&ntrials);
      xtree->GetEntry(0);
      fTrials = ntrials;
      fCrossSection = xsection;
      fxsec->Close();
    }
    #ifdef DEBUGMODE
      AliInfo("Notify ended.");
    #endif
  }
  return kTRUE;
}

//________________________________________________________________________
inline void AliAnalysisTaskQualityAssurancePA::FillHistogram(const char* runNumber, const char * key, Double_t x)
{
  if(strcmp(runNumber,"ALL") != 0)
    FillHistogram("ALL", key, x);

  TH1* tmpHist = static_cast<TH1*>(fOutputList->FindObject(GetHistoName(runNumber, key)));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",GetHistoName(runNumber, key))) ;
    return;
  }

  tmpHist->Fill(x);
}

//________________________________________________________________________
inline void AliAnalysisTaskQualityAssurancePA::FillHistogram(const char* runNumber, const char * key, Double_t x, Double_t y)
{
  if(strcmp(runNumber,"ALL") != 0)
    FillHistogram("ALL", key, x, y);

  TH1* tmpHist = static_cast<TH1*>(fOutputList->FindObject(GetHistoName(runNumber, key)));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",GetHistoName(runNumber, key)));
    return;
  }

  if (tmpHist->IsA()->GetBaseClass("TH1"))
    static_cast<TH1*>(tmpHist)->Fill(x,y); // Fill x with y
  else if (tmpHist->IsA()->GetBaseClass("TH2"))
    static_cast<TH2*>(tmpHist)->Fill(x,y); // Fill x,y with 1
}

//________________________________________________________________________
inline void AliAnalysisTaskQualityAssurancePA::FillHistogram(const char* runNumber, const char * key, Double_t x, Double_t y, Double_t add)
{
  if(strcmp(runNumber,"ALL") != 0)
    FillHistogram("ALL", key, x, y, add);


  TH2* tmpHist = static_cast<TH2*>(fOutputList->FindObject(GetHistoName(runNumber, key)));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",GetHistoName(runNumber, key)));
    return;
  }
  
  tmpHist->Fill(x,y,add);
}
//________________________________________________________________________
template <class T> T* AliAnalysisTaskQualityAssurancePA::AddHistogram1D(const char* runNumber, const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, const char* xTitle, const char* yTitle)
{
  T* tmpHist = new T(GetHistoName(runNumber, name), GetHistoName(runNumber, title), xBins, xMin, xMax);

  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fHistList->Add(tmpHist);
  fHistCount++;
  
  return tmpHist;
}

//________________________________________________________________________
template <class T> T* AliAnalysisTaskQualityAssurancePA::AddHistogram2D(const char* runNumber, const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, const char* xTitle, const char* yTitle, const char* zTitle)
{
  T* tmpHist = new T(GetHistoName(runNumber, name), GetHistoName(runNumber, title), xBins, xMin, xMax, yBins, yMin, yMax);
  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->GetZaxis()->SetTitle(zTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fHistList->Add(tmpHist);
  fHistCount++;

  return tmpHist;
}

//________________________________________________________________________
void AliAnalysisTaskQualityAssurancePA::Terminate(Option_t *)
{
  PostData(1, fOutputList);

  // Mandatory
  fOutputList = dynamic_cast<TList*> (GetOutputData(1)); // '1' refers to the output slot
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}

//________________________________________________________________________
AliAnalysisTaskQualityAssurancePA::~AliAnalysisTaskQualityAssurancePA()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputList;
  }
}

//________________________________________________________________________
void AliAnalysisTaskQualityAssurancePA::UserCreateOutputObjects()
{
  // called once to create user defined output objects like histograms, plots etc. 
  // and to put it on the output list.
  // Note: Saving to file with e.g. OpenFile(0) is must be before creating other objects.

  fRandom = new TRandom3(0);
  
  fOutputList = new TList();
  fOutputList->SetOwner(); // otherwise it produces leaks in merging

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskQualityAssurancePA::UserExec(Option_t *) 
{
  #ifdef DEBUGMODE
    AliInfo("UserExec() started.");
  #endif

  Calculate(InputEvent());
        
  PostData(1, fOutputList);
}
