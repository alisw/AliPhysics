#ifndef ALIANALYSISTASKSE_H
#include <Riostream.h>
#include <TROOT.h>
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
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliAnalysisUtils.h"


#include "AliAnalysisTaskChargedJetsPA.h"

//TODO: Not accessing the particles when using MC
//TODO: FillHistogram can be done better with virtual TH1(?)
ClassImp(AliAnalysisTaskChargedJetsPA)

// ######################################################################################## DEFINE HISTOGRAMS
void AliAnalysisTaskChargedJetsPA::Init()
{
  #ifdef DEBUGMODE
    AliInfo("Creating histograms.");
  #endif

  AddHistogram1D<TH1D>("hNumberEvents", "Number of events (0 = before, 1 = after vertex cuts)", "", 2, 0, 2, "#Delta z(cm)","N^{Events}/cut");

  // NOTE: Jet histograms
  if (fAnalyzeJets)
  {
    // ######## Jet spectra
    AddHistogram2D<TH2D>("hJetPt", "Jets p_{T} distribution", "", 500, -50., 200., 5, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedRC", "Jets p_{T} distribution, RC background subtracted", "", 500, -50., 200.,5, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKT", "Jets p_{T} distribution, KT background subtracted", "", 500, -50., 200., 5, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedTR", "Jets p_{T} distribution, TR background subtracted", "", 500, -50., 200.,5, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedRCPosEta", "Jets p_{T} distribution, RC background subtracted, #eta > 0.2", "", 500, -50., 200.,5, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTPosEta", "Jets p_{T} distribution, KT background subtracted, #eta > 0.2", "", 500, -50., 200., 5, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedTRPosEta", "Jets p_{T} distribution, TR background subtracted, #eta > 0.2", "", 500, -50., 200.,5, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedRCNegEta", "Jets p_{T} distribution, RC background subtracted, #eta < -0.2", "", 500, -50., 200.,5, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTNegEta", "Jets p_{T} distribution, KT background subtracted, #eta < -0.2", "", 500, -50., 200., 5, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedTRNegEta", "Jets p_{T} distribution, TR background subtracted, #eta < -0.2", "", 500, -50., 200.,5, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");

    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedRCEtaWindow", "Jets p_{T} distribution, RC background (in #eta window around jet) subtracted", "", 500, -50., 200., 5, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedKTEtaWindow", "Jets p_{T} distribution, KT background (in #eta window around jet) subtracted", "", 500, -50., 200., 5, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hJetPtBgrdSubtractedTREtaWindow", "Jets p_{T} distribution, TR background (in #eta window around jet) subtracted", "", 500, -50., 200., 5, 0, 100, "p_{T} (GeV/c)","Centrality","dN^{Jets}/dp_{T}");

    // ######## Jet stuff
    AddHistogram1D<TH1D>("hJetCountAll", "Number of Jets", "", 200, 0., 200., "N jets","dN^{Events}/dN^{Jets}");
    AddHistogram1D<TH1D>("hJetCountAccepted", "Number of accepted Jets", "", 200, 0., 200., "N jets","dN^{Events}/dN^{Jets}");
    AddHistogram1D<TH1D>("hLeadingJetPt", "Leading jet p_{T}", "", 500,  0, 100, "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hSecondLeadingJetPt", "Second Leading jet p_{T}", "", 500,  0, 100, "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hJetDeltaPhi", "Jets combinatorial #Delta #phi", "", 250, 0., TMath::Pi(), "#Delta #phi","dN^{Jets}/d(#Delta #phi)");
    AddHistogram1D<TH1D>("hLeadingJetDeltaPhi", "1st and 2nd leading jet #Delta #phi", "", 250, 0., TMath::Pi(), "#Delta #phi","dN^{Jets}/d(#Delta #phi)");
  }

  // NOTE: Jet background histograms
  if (fAnalyzeBackground)
  {
    // ########## Different background estimates
    AddHistogram2D<TH2D>("hRCBackground", "RC background density (2 leading jets excluded)", "LEGO2", 400, 0., 40., 5, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hKTBackground", "KT background density (2 leading jets excluded)", "LEGO2", 400, 0., 40., 5, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hTRBackground", "TR background density (2 leading jets excluded)", "LEGO2", 400, 0., 40., 5, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");

    // ########## Delta Pt
    AddHistogram2D<TH2D>("hDeltaPtKT", "Background fluctuations #delta p_{T} (KT)", "", 600, -40., 80., 5, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtRC", "Background fluctuations #delta p_{T} (RC)", "", 600, -40., 80., 5, 0, 100,  "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtTR", "Background fluctuations #delta p_{T} (TR)", "", 600, -40., 80., 5, 0, 100,  "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtKTNoExcl", "Background fluctuations #delta p_{T} (KT, no leading jet correction)", "", 600, -40., 80., 5, 0, 100, "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtRCNoExcl", "Background fluctuations #delta p_{T} (RC, no leading jet correction)", "", 600, -40., 80., 5, 0, 100,  "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");
    AddHistogram2D<TH2D>("hDeltaPtTRNoExcl", "Background fluctuations #delta p_{T} (TR, no leading jet correction)", "", 600, -40., 80., 5, 0, 100,  "#delta p_{T} (GeV/c)","Centrality","dN^{Jets}/d#delta p_{T}");

    // ########## Min bias background in eta bins
    AddHistogram2D<TH2D>("hRCBackgroundEtaBins", "RC background density (2 leading jets excluded)", "LEGO2", 400, 0., 40., 5, -0.5, +0.5, "#rho (GeV/c)","#eta", "dN^{Events}/d#rho d#eta");
    AddHistogram2D<TH2D>("hKTBackgroundEtaBins", "KT background density (2 leading jets excluded)", "LEGO2", 400, 0., 40., 5, -0.5, +0.5, "#rho (GeV/c)","#eta", "dN^{Events}/d#rho d#eta");
    AddHistogram2D<TH2D>("hTRBackgroundEtaBins", "TR background density (2 leading jets excluded)", "LEGO2", 400, 0., 40., 5, -0.5, +0.5, "#rho (GeV/c)","#eta", "dN^{Events}/d#rho d#eta");

    // ########## Min bias background in sliding eta bins
    AddHistogram2D<TH2D>("hRCBackgroundEtaWindow", "RC background density (2 leading jets excluded)", "LEGO2", 400, 0., 40., 5, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho d#eta");
    AddHistogram2D<TH2D>("hKTBackgroundEtaWindow", "KT background density (2 leading jets excluded)", "LEGO2", 400, 0., 40., 5, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho d#eta");
    AddHistogram2D<TH2D>("hTRBackgroundEtaWindow", "TR background density (2 leading jets excluded)", "LEGO2", 400, 0., 40., 5, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho d#eta");


    // ########## Dijet stuff
    AddHistogram1D<TH1D>("hDijetLeadingJetPt", "Dijet leading jet p_{T} distribution", "", 500, 0., 100., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram1D<TH1D>("hDijetConstituentsPt", "Dijet constituents p_{T} distribution", "", 500, 0., 100., "p_{T} (GeV/c)","dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>("hDijetPtCorrelation", "Dijet constituents p_{T} correlation", "COLZ", 500, 5., 100., 500, 5., 100., "1st leading jet p_{T} (GeV/c)","2nd leading jet p_{T} (GeV/c)","dN^{Dijets}/d^{2}p_{T}");
    AddHistogram2D<TH2D>("hDijetBackground", "Background density (dijets excluded)", "", 200, 0., 20., 5, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
    AddHistogram2D<TH2D>("hDijetBackgroundPerpendicular", "Background density (dijets excluded)", "", 200, 0., 20., 5, 0, 100, "#rho (GeV/c)","Centrality", "dN^{Events}/d#rho");
  }

  // NOTE: Pythia histograms
  if (fAnalyzePythia)
  {
    AddHistogram1D<TH1D>("hPythiaPtHard", "Pythia p_{T} hard distribution", "", 2000, 0, 400, "p_{T} hard","dN^{Events}/dp_{T,hard}");
    AddHistogram1D<TProfile>("hPythiaXSection", "Pythia cross section distribution", "", fNumPtHardBins+2, -1, fNumPtHardBins+1, "p_{T} hard bin","dN^{Events}/dp_{T,hard}");
    AddHistogram1D<TH1D>("hPythiaNTrials", "Pythia trials (no correction for manual cuts)", "", fNumPtHardBins+2, -1, fNumPtHardBins+1, "p_{T} hard bin", "Trials");
  }

  // register Histograms
  for (Int_t i = 0; i < fHistCount; i++)
  {
    fOutputList->Add(fHistList->At(i));
  }
  
  PostData(1,fOutputList); // important for merging

}

//________________________________________________________________________
AliAnalysisTaskChargedJetsPA::AliAnalysisTaskChargedJetsPA(const char *name, const char* trackArrayName, const char* jetArrayName, const char* backgroundJetArrayName) : AliAnalysisTaskSE(name), fOutputList(0), fAnalyzeJets(1), fAnalyzeBackground(1), fAnalyzePythia(0), fHasTracks(0), fHasJets(0), fHasBackgroundJets(0), fIsMC(0), fJetArray(0), fTrackArray(0), fBackgroundJetArray(0), fJetArrayName(0), fTrackArrayName(0), fBackgroundJetArrayName(0), fNumPtHardBins(11), fRandConeRadius(0.4), fSignalJetRadius(0.4), fBackgroundJetRadius(0.4), fTRBackgroundConeRadius(0.4), fNumberRandCones(8), fNumberExcludedJets(2), fDijetMaxAngleDeviation(10.0), fSignalJetEtaWindow(0.5), fBackgroundJetEtaWindow(0.5), fTrackEtaWindow(0.9), fVertexWindow(10.0), fVertexMaxR(1.0), fMinTrackPt(0.150), fMinJetPt(1.0), fMinJetArea(0.4), fMinBackgroundJetPt(0.15), fMinDijetLeadingPt(10.0), fCentralityType("V0A"), fFirstLeadingJet(0), fSecondLeadingJet(0), fNumberSignalJets(0), fCrossSection(0.0), fTrials(0.0),  fRandom(0), fHelperClass(0), fInitialized(0), fTaskInstanceCounter(0), fHistList(0), fHistCount(0)
{
  #ifdef DEBUGMODE
    AliInfo("Calling constructor.");
  #endif

  // Every instance of this task gets his own number
  static Int_t instance = 0;
  fTaskInstanceCounter = instance;
  instance++;

  fTrackArrayName = new TString(trackArrayName);
  if (fTrackArrayName->Contains("MCParticles")) //TODO: Not working for now
    fIsMC = kTRUE;

  fJetArrayName = new TString(jetArrayName);
  if (strcmp(fJetArrayName->Data(),"") == 0)
    fAnalyzeJets = kFALSE;
  else
    fAnalyzeJets = kTRUE;
    
  fBackgroundJetArrayName = new TString(backgroundJetArrayName);
  if (strcmp(fBackgroundJetArrayName->Data(),"") == 0)
    fAnalyzeBackground = kFALSE;
  else
    fAnalyzeBackground = kTRUE;

  DefineOutput(1, TList::Class());
 
  fHistList = new TList();

  #ifdef DEBUGMODE
    AliInfo("Constructor done.");
  #endif
  
}

//________________________________________________________________________
inline Double_t AliAnalysisTaskChargedJetsPA::GetConePt(Double_t eta, Double_t phi, Double_t radius)
{
  Double_t tmpConePt = 0.0;

  for (Int_t i = 0; i < fTrackArray->GetEntries(); i++)
  {
    AliVTrack* tmpTrack = static_cast<AliVTrack*>(fTrackArray->At(i));
    if (IsTrackInAcceptance(tmpTrack))
      if(IsTrackInCone(tmpTrack, eta, phi, radius))
        tmpConePt = tmpConePt + tmpTrack->Pt();
  }
  return tmpConePt;
}


//________________________________________________________________________
inline Double_t AliAnalysisTaskChargedJetsPA::GetPtHard()
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
inline Int_t AliAnalysisTaskChargedJetsPA::GetPtHardBin()
{
  // ########## PT HARD BIN EDGES
  const Int_t kPtHardLowerEdges[] =  { 0, 5,11,21,36,57, 84,117,152,191,234};
  const Int_t kPtHardHigherEdges[] = { 5,11,21,36,57,84,117,152,191,234,1000000};

  Int_t tmpPtHardBin = 0;
  Double_t tmpPtHard = GetPtHard();
 
  for (tmpPtHardBin = 0; tmpPtHardBin <= fNumPtHardBins; tmpPtHardBin++)
    if (tmpPtHard >= kPtHardLowerEdges[tmpPtHardBin] && tmpPtHard < kPtHardHigherEdges[tmpPtHardBin])
      break;

  return tmpPtHardBin;
}


//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsTrackInCone(AliVTrack* track, Double_t eta, Double_t phi, Double_t radius)
{
  // This is to use a full cone in phi even at the edges of phi (2pi -> 0) (0 -> 2pi)
  Double_t trackPhi = 0.0;
  if (track->Phi() > (TMath::TwoPi() - (radius-phi)))
    trackPhi = track->Phi() - TMath::TwoPi();
  else if (track->Phi() < (phi+radius - TMath::TwoPi()))
    trackPhi = track->Phi() + TMath::TwoPi();
  else
    trackPhi = track->Phi();
  
  if ( TMath::Abs(trackPhi-phi)*TMath::Abs(trackPhi-phi) + TMath::Abs(track->Eta()-eta)*TMath::Abs(track->Eta()-eta) <= radius*radius)
    return kTRUE;
  
  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsTrackInAcceptance(AliVParticle* track)
{
  if (track != 0)
    if (TMath::Abs(track->Eta()) <= fTrackEtaWindow)
      if (track->Pt() >= fMinTrackPt)
        return kTRUE;

  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsBackgroundJetInAcceptance(AliEmcalJet *jet)
{   
  if (jet != 0)
    if (TMath::Abs(jet->Eta()) <= fBackgroundJetEtaWindow)
      if (jet->Pt() >= fMinBackgroundJetPt)
        return kTRUE;

  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsSignalJetInAcceptance(AliEmcalJet *jet)
{   
  if (jet != 0)
    if (TMath::Abs(jet->Eta()) <= fSignalJetEtaWindow)
      if (jet->Pt() >= fMinJetPt)
        if (jet->Area() >= fMinJetArea)
          return kTRUE;

  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsPA::IsDijet(AliEmcalJet *jet1, AliEmcalJet *jet2)
{   
  // Output from GetDeltaPhi is < pi in any case
  if ((jet1 != 0) && (jet2 != 0))
    if((TMath::Pi() - GetDeltaPhi(jet1->Phi(),jet2->Phi())) < fDijetMaxAngleDeviation)
      if ((jet1->Pt() > fMinDijetLeadingPt) && (jet2->Pt() > fMinDijetLeadingPt)) //TODO: Introduce recoil jet?
        return kTRUE;

  return kFALSE;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::ExecOnce()
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
      AliWarning(Form("%s: Could not retrieve tracks %s! This is OK, if tracks are not demanded.", GetName(), fTrackArrayName->Data())); 
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

  // Check for jet array
  if (strcmp(fJetArrayName->Data(), "") != 0)
  {
    fJetArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetArrayName->Data()));
    fHasJets = kTRUE;

    if (!fJetArray) 
    {
      AliWarning(Form("%s: Could not retrieve jets %s! This is OK, if jets are not demanded.", GetName(), fJetArrayName->Data())); 
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

  // Check for background object
  if (strcmp(fBackgroundJetArrayName->Data(), "") != 0)
  {
    fBackgroundJetArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fBackgroundJetArrayName->Data()));
    fHasBackgroundJets = kTRUE;
    if (!fBackgroundJetArray)
    {
      AliInfo(Form("%s: Could not retrieve background jets %s! This is OK, if background is not demanded.", GetName(), fBackgroundJetArrayName->Data())); 
      fHasBackgroundJets = kFALSE;
    }
  }

  // Look, if initialization is OK
  if (!fHasTracks && fAnalyzeBackground)
  {
    AliError(Form("%s: Tracks NOT successfully casted although demanded! Deactivating background analysis",GetName()));
    fAnalyzeBackground = kFALSE;
  }
  if ((!fHasJets && fAnalyzeJets) || (!fHasJets && fAnalyzeBackground))
  {
    AliError(Form("%s: Jets NOT successfully casted although demanded!  Deactivating jet- and background analysis",GetName()));
    fAnalyzeJets = kFALSE;
    fAnalyzeBackground = kFALSE;
  }
  if (!fHasBackgroundJets && fAnalyzeBackground)
  {
    AliError(Form("%s: Background NOT successfully casted although demanded!  Deactivating background analysis",GetName()));
    fAnalyzeBackground = kFALSE;
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
void AliAnalysisTaskChargedJetsPA::GetSignalJets()
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
Int_t AliAnalysisTaskChargedJetsPA::GetLeadingJets(TClonesArray* jetArray, Int_t* jetIDArray, Bool_t isSignalJets)
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

    if(isSignalJets)
    {
      if (!IsSignalJetInAcceptance(jet)) continue;
    }
    else
    {
      if (!IsBackgroundJetInAcceptance(jet)) continue;
    }    

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
Double_t AliAnalysisTaskChargedJetsPA::GetCorrectedJetPt(AliEmcalJet* jet, Double_t background)
{
  #ifdef DEBUGMODE
    AliInfo("Getting corrected jet spectra.");
  #endif

  if(!jet)
  {
    AliError("Jet pointer passed to GetCorrectedJet() not valid!");
    return -1.0;
  }

  Double_t correctedPt = -1.0;

  // if the passed background is not valid, do not subtract it
  if(background < 0)
    background = 0;

  // Subtract background
  correctedPt = jet->Pt() - background * jet->Area();

  #ifdef DEBUGMODE
    AliInfo("Got corrected jet spectra.");
  #endif 

  return correctedPt;
}



//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetDeltaPt(Double_t& deltaPt, Double_t rho, Bool_t leadingJetExclusion)
{
  #ifdef DEBUGMODE
    AliInfo("Getting Delta Pt.");
  #endif

  // Define an invalid delta pt
  deltaPt = -10000.0;

  // Define eta range
  Double_t etaMin, etaMax;
  etaMin = -(fTrackEtaWindow-fRandConeRadius);
  etaMax = +(fTrackEtaWindow-fRandConeRadius);

  // Define random cone
  Bool_t coneValid = kTRUE;
  Double_t tmpRandConeEta = etaMin + fRandom->Rndm()*(etaMax-etaMin);
  Double_t tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();

  // if there is a jet, check for overlap if demanded
  if(leadingJetExclusion)
  {
    for (Int_t i = 0; i<fNumberSignalJets; i++)
    {
      AliEmcalJet* tmpJet = fSignalJets[i];

      Double_t excludedJetPhi = tmpJet->Phi();
      Double_t excludedJetEta = tmpJet->Eta();
      Double_t tmpDeltaPhi = GetDeltaPhi(tmpRandConePhi, excludedJetPhi);

      // Check, if cone has overlap with jet
      if ( tmpDeltaPhi*tmpDeltaPhi + TMath::Abs(tmpRandConeEta-excludedJetEta)*TMath::Abs(tmpRandConeEta-excludedJetEta) <= fRandConeRadius*fRandConeRadius)
      {
        // Define probability to exclude the RC
        Double_t probability = 1 - (fNumberSignalJets-1)/fNumberSignalJets;

        // Only exclude cone with a given probability
        if (fRandom->Rndm()<=probability)
        {
          coneValid = kFALSE;
          break;
        }
      }
    }
  }

  // Get the cones' pt and calculate delta pt
  if (coneValid)
    deltaPt = GetConePt(tmpRandConeEta,tmpRandConePhi,fRandConeRadius) - (rho*fRandConeRadius*fRandConeRadius*TMath::Pi());

  #ifdef DEBUGMODE
    AliInfo("Got Delta Pt.");
  #endif
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetKTBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoMedian, Double_t& areaMean, Double_t etaMin, Double_t etaMax, Bool_t excludeSignalJets)
{
  #ifdef DEBUGMODE
    AliInfo("Getting KT background density.");
  #endif

  // static declaration. Advantage: more speed. Disadvantage: Problematic for events with more than 1024 jets :)
  static Double_t tmpRhos[1024];
  static Double_t tmpAreas[1024];
  Int_t maxJetIds[]   = {-1, -1}; // Indices for excludes jets (up to two)

  // Setting invalid values
  rhoMedian = -1.0;
  areaMean= -1.0;

  GetLeadingJets(fBackgroundJetArray, &maxJetIds[0], kFALSE);
  // Exclude UP TO numberExcludeLeadingJets
  if (fNumberSignalJets < 2)
    numberExcludeLeadingJets = fNumberSignalJets;

  // Check if etaMin/etaMax is given correctly
  if(etaMin < -fBackgroundJetEtaWindow)
    etaMin = -fBackgroundJetEtaWindow;
  if(etaMax > fBackgroundJetEtaWindow)
    etaMax = fBackgroundJetEtaWindow;

  if ((etaMin == 0) && (etaMax == 0))
  {
    etaMin = -fBackgroundJetEtaWindow;
    etaMax = +fBackgroundJetEtaWindow;
  }

  Int_t jetCountAccepted = 0;
  for (Int_t i = 0; i < fBackgroundJetArray->GetEntries(); i++)
  {
    AliEmcalJet* backgroundJet = static_cast<AliEmcalJet*>(fBackgroundJetArray->At(i));
    if (!backgroundJet)
    {
      AliError(Form("%s: Could not receive jet %d", GetName(), i));
      continue;
    } 

    // Exclude signal jets
    if (excludeSignalJets)
    {
      Bool_t isOverlapping = kFALSE;
      for(Int_t j=0;j<numberExcludeLeadingJets;j++)
      {
        AliEmcalJet* signalJet = NULL;

        if (j==0)
          signalJet = fFirstLeadingJet;
        else if (j==1)
          signalJet = fSecondLeadingJet;
        else
          AliFatal("Trying to exclude more than 2 signal jets in KT background -- not implemented.");


        Double_t tmpDeltaPhi = GetDeltaPhi(backgroundJet->Phi(), signalJet->Phi());
        
        if ( tmpDeltaPhi*tmpDeltaPhi + TMath::Abs(signalJet->Eta()-backgroundJet->Eta())*TMath::Abs(signalJet->Eta()-backgroundJet->Eta()) <= fSignalJetRadius*fSignalJetRadius)
        {
          isOverlapping = kTRUE;
          break;
        }
      }
      if(isOverlapping)
        continue;
    }
    else  // exclude leading background jets
    {
      if (numberExcludeLeadingJets > 0)
        if (i == maxJetIds[0])
          continue;
      if (numberExcludeLeadingJets > 1)
        if (i == maxJetIds[1])
          continue;
    }      

    if (!IsBackgroundJetInAcceptance(backgroundJet))
      continue;
    if (!((backgroundJet->Eta() >= etaMin) && (backgroundJet->Eta() < etaMax)))
      continue;

    
    tmpRhos[jetCountAccepted] = backgroundJet->Pt() / backgroundJet->Area();
    tmpAreas[jetCountAccepted] = backgroundJet->Area();
    jetCountAccepted++;
  }

  if (jetCountAccepted > 0)
  {
    rhoMedian = TMath::Median(jetCountAccepted, tmpRhos);
    areaMean   = TMath::Mean(jetCountAccepted, tmpAreas);
  }
  #ifdef DEBUGMODE
    AliInfo("Got KT background density.");
  #endif
}


//________________________________________________________________________
Int_t AliAnalysisTaskChargedJetsPA::GetRCBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoMean, Double_t& rhoMedian, Double_t etaMin, Double_t etaMax, Int_t numberRandCones)
{
  #ifdef DEBUGMODE
    AliInfo("Getting RC background density.");
  #endif

  if(numberRandCones == 0)
    numberRandCones = fNumberRandCones;

  std::vector<AliEmcalJet> tmpCones(numberRandCones);

  // Setting invalid values
  rhoMean = -1.0;
  rhoMedian = -1.0;

  // Exclude UP TO numberExcludeLeadingJets
  if (fNumberSignalJets < 2)
    numberExcludeLeadingJets = fNumberSignalJets;

  // Search given amount of RCs
  Int_t numAcceptedRCs = 0;
  for(Int_t i=0;i<numberRandCones;i++)
  {
    Double_t tmpRandConeEta = 0.0;
    Double_t tmpRandConePhi = 0.0;
    Double_t excludedJetEta = 0.0;
    Double_t excludedJetPhi = 0.0;

    // Search random cone in acceptance with no overlap with already excluded jets (leading jets and random cones)
    Bool_t coneValid = kTRUE;

    // Check if etaMin/etaMax is given correctly
    if(etaMin < -fSignalJetEtaWindow)
      etaMin = -fSignalJetEtaWindow;
    if(etaMax > fSignalJetEtaWindow)
      etaMax = fSignalJetEtaWindow;

    // Set the random cone position
    if ((etaMin == 0) && (etaMax == 0))
      tmpRandConeEta = (fTrackEtaWindow-fRandConeRadius)*(2.0*fRandom->Rndm()-1.0); // full RC is in acceptance
    else
      tmpRandConeEta = etaMin + fRandom->Rndm()*(etaMax-etaMin);

    tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();

    // Go through all excluded leading jets and check if there's an overlap
     
    for(Int_t j=0;j<numberExcludeLeadingJets;j++)
    {
      AliEmcalJet* tmpJet = NULL;

      if (j==0)
        tmpJet = fFirstLeadingJet;
      else if (j==1)
        tmpJet = fSecondLeadingJet;
      else
        AliFatal("Trying to exclude more than 2 jets in RC background -- not implemented.");

      excludedJetPhi = tmpJet->Phi();
      excludedJetEta = tmpJet->Eta();
      Double_t tmpDeltaPhi = GetDeltaPhi(tmpRandConePhi, excludedJetPhi);
      
      if ( tmpDeltaPhi*tmpDeltaPhi + TMath::Abs(tmpRandConeEta-excludedJetEta)*TMath::Abs(tmpRandConeEta-excludedJetEta) <= fRandConeRadius*fRandConeRadius)
      {
        coneValid = kFALSE;
        break;
      }
    }

    // RC is accepted, so save it
    if(coneValid)
    {
      AliEmcalJet tmpJet(GetConePt(tmpRandConeEta, tmpRandConePhi, fRandConeRadius), tmpRandConeEta, tmpRandConePhi, 0.0);
      tmpCones[numAcceptedRCs] = tmpJet;
      numAcceptedRCs++;
    }
  }

  // Calculate Rho and the mean from the RCs (no excluded jets are considered!)
  if(numAcceptedRCs > 0)
  {
    std::vector<Double_t> tmpRho(numAcceptedRCs);
    for (Int_t i=0; i<numAcceptedRCs;i++)
      tmpRho[i] = tmpCones[i].Pt()/(fRandConeRadius*fRandConeRadius*TMath::Pi());

    rhoMean = TMath::Mean(tmpRho.begin(), tmpRho.end());
    rhoMedian = 0.0; // NOT IMPLEMENTED because TMath::Median is not working with iterators
  }
    
  #ifdef DEBUGMODE
    AliInfo("Got RC background density.");
  #endif
  return numAcceptedRCs;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetTRBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoMean, Double_t& area, Double_t etaMin, Double_t etaMax)
{
  #ifdef DEBUGMODE
    AliInfo("Getting TR background density.");
  #endif

  Double_t summedTracksPt = 0.0;

  // Check if etaMin/etaMax is given correctly
  if(etaMin < -fTrackEtaWindow)
    etaMin = -fTrackEtaWindow;
  if(etaMax > fTrackEtaWindow)
    etaMax = fTrackEtaWindow;

  if ((etaMin == 0) && (etaMax == 0))
  {
    etaMin = -fTrackEtaWindow;
    etaMax = +fTrackEtaWindow;
  }

  // Setting invalid values
  rhoMean = -1.0;
  area = -1.0;
  // Exclude UP TO numberExcludeLeadingJets
  if (fNumberSignalJets < 2)
    numberExcludeLeadingJets = fNumberSignalJets;


  Int_t   trackCount = fTrackArray->GetEntries();
  Int_t   trackCountAccepted = 0;
  for (Int_t i = 0; i < trackCount; i++)
  {
    Bool_t  trackValid = kTRUE;
    AliVTrack* tmpTrack = static_cast<AliVTrack*>(fTrackArray->At(i));
    if (IsTrackInAcceptance(tmpTrack))
      if ((tmpTrack->Eta() >= etaMin) && (tmpTrack->Eta() < etaMax))
      {
        for (Int_t j = 0; j < numberExcludeLeadingJets; j++)
        {
          AliEmcalJet* tmpJet = NULL;
          if (j==0)
            tmpJet = fFirstLeadingJet;
          else if (j==1)
            tmpJet = fSecondLeadingJet;
          else
            AliFatal("Trying to exclude more than 2 jets in track background -- not implemented.");

          if (IsTrackInCone(tmpTrack, tmpJet->Eta(), tmpJet->Phi(), fTRBackgroundConeRadius))
          {
            trackValid = kFALSE;
            break;
          }
        }
        if (trackValid)
        {
          // Add track pt to array
          summedTracksPt = summedTracksPt + tmpTrack->Pt();
          trackCountAccepted++;
        }
      }
  }

  if (trackCountAccepted > 0)
  {
    Double_t tmpArea = 0.0;

    tmpArea = (2.0*fTrackEtaWindow) * TMath::TwoPi() * (etaMax-etaMin)/(2.0*fTrackEtaWindow); // area of the used eta strip
    
    // Now: exclude the part of the leading jet that is in the strip
    if (numberExcludeLeadingJets == 2)
      tmpArea = tmpArea*(1.0-MCGetOverlapCircleRectancle(fFirstLeadingJet->Eta(), fFirstLeadingJet->Phi(), fTRBackgroundConeRadius, etaMin, etaMax, 0., TMath::TwoPi()) -MCGetOverlapCircleRectancle(fSecondLeadingJet->Eta(), fSecondLeadingJet->Phi(), fTRBackgroundConeRadius, etaMin, etaMax, 0., TMath::TwoPi()));
    else if (numberExcludeLeadingJets == 1)
      tmpArea = tmpArea*(1.0-MCGetOverlapCircleRectancle(fFirstLeadingJet->Eta(), fFirstLeadingJet->Phi(), fTRBackgroundConeRadius, etaMin, etaMax, 0., TMath::TwoPi()));
   
    rhoMean = summedTracksPt/tmpArea;
    area  = tmpArea;
  }

  #ifdef DEBUGMODE
    AliInfo("Got TR background density.");
  #endif
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::GetTRBackgroundDensity(Int_t numberExcludeLeadingJets, Double_t& rhoMean, Double_t& area, AliEmcalJet* excludeJet1, AliEmcalJet* excludeJet2, Bool_t doSearchPerpendicular)
{
  #ifdef DEBUGMODE
    AliInfo("Getting TR background density.");
  #endif

  // Setting invalid values
  Double_t summedTracksPt = 0.0;
  rhoMean = -1.0;
  area = -1.0;

  Double_t tmpRadius = 0.0;
  if (doSearchPerpendicular)
    tmpRadius = 0.5*TMath::Pi(); // exclude 90 degrees around jets
  else
    tmpRadius = fSignalJetRadius;
    
  numberExcludeLeadingJets = 2; // dijet is excluded here in any case



  if (!fTrackArray || !fJetArray)
  {
    AliError("Could not get the track/jet array to calculate track rho!");
    return;
  }

  Int_t   trackCount = fTrackArray->GetEntries();
  Int_t   trackCountAccepted = 0;
  for (Int_t i = 0; i < trackCount; i++)
  {
    AliVTrack* tmpTrack = static_cast<AliVTrack*>(fTrackArray->At(i));
    if (IsTrackInAcceptance(tmpTrack))
    {
      if (IsTrackInCone(tmpTrack, excludeJet1->Eta(), excludeJet1->Phi(), tmpRadius))
        continue;

      if (numberExcludeLeadingJets > 1)
        if (IsTrackInCone(tmpTrack, excludeJet2->Eta(), excludeJet2->Phi(), tmpRadius))
          continue;

        // Add track pt to array
        summedTracksPt = summedTracksPt + tmpTrack->Pt();
        trackCountAccepted++;
    }
  }

  if (trackCountAccepted > 0)
  {
    Double_t tmpArea = 2.0*fTrackEtaWindow*TMath::TwoPi()  - 2*(tmpRadius*tmpRadius * TMath::Pi()); //TPC area - excluding jet area
    rhoMean = summedTracksPt/tmpArea;
    area = tmpArea;
  }

  #ifdef DEBUGMODE
    AliInfo("Got TR background density.");
  #endif
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::Calculate(AliVEvent* event)
{
  #ifdef DEBUGMODE
    AliInfo("Starting Calculate().");
  #endif
  ////////////////////// NOTE: initialization & casting

  // Additional cuts
  FillHistogram("hNumberEvents", 0.5); // number of events before manual cuts

  if(!fHelperClass->IsVertexSelected2013pA(event))
    return;

  FillHistogram("hNumberEvents", 1.5); // number of events after manual cuts

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Init done.");
  #endif

  ////////////////////// NOTE: Get Centrality, (Leading)Signal jets and Background

  // Get centrality
  AliCentrality* tmpCentrality = NULL;
  tmpCentrality = event->GetCentrality();
  Double_t centralityPercentile = 0.0;
  if (tmpCentrality != NULL)
    centralityPercentile = tmpCentrality->GetCentralityPercentile(fCentralityType.Data());

  // Get jets
  if (fAnalyzeBackground || fAnalyzeJets)
    GetSignalJets();

  // Get background estimates
  Double_t              backgroundKTMedian = -1.0;
  Double_t              backgroundRCMean = -1.0;
  Double_t              backgroundRCMedian = -1.0;
  Double_t              backgroundTRMean = -1.0;
  Double_t              backgroundKTAreaMean = -1.0;
  Double_t              backgroundTRAreaMean = -1.0;
  Double_t              dijetBackground = -1.0;
  Double_t              dijetBackgroundPerpendicular = -1.0;

  if (fAnalyzeBackground)
  {
    GetRCBackgroundDensity    (fNumberExcludedJets, backgroundRCMean, backgroundRCMedian);
    GetTRBackgroundDensity    (fNumberExcludedJets, backgroundTRMean, backgroundTRAreaMean);
    GetKTBackgroundDensity    (fNumberExcludedJets, backgroundKTMedian, backgroundKTAreaMean);
  }

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Centrality&SignalJets&Background-Calculation done.");
  #endif

  ////////////////////// NOTE: Jet analysis and calculations

  if (fAnalyzeJets)
  {
    // ### SIGNAL JET ANALYSIS
    for (Int_t i = 0; i<fNumberSignalJets; i++)
    {
      AliEmcalJet* tmpJet = fSignalJets[i];

      // Jet spectra
      FillHistogram("hJetPt", tmpJet->Pt(), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedRC", GetCorrectedJetPt(tmpJet, backgroundRCMean), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedKT", GetCorrectedJetPt(tmpJet, backgroundKTMedian), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedTR", GetCorrectedJetPt(tmpJet, backgroundTRMean), centralityPercentile);

      // Jet spectra  (pos+neg eta)
      if(tmpJet->Eta() > 0.2)
      {
        FillHistogram("hJetPtBgrdSubtractedRCPosEta", GetCorrectedJetPt(tmpJet, backgroundRCMean), centralityPercentile);
        FillHistogram("hJetPtBgrdSubtractedKTPosEta", GetCorrectedJetPt(tmpJet, backgroundKTMedian), centralityPercentile);
        FillHistogram("hJetPtBgrdSubtractedTRPosEta", GetCorrectedJetPt(tmpJet, backgroundTRMean), centralityPercentile);
      }
      else if(tmpJet->Eta() < -0.2)
      {
        FillHistogram("hJetPtBgrdSubtractedRCNegEta", GetCorrectedJetPt(tmpJet, backgroundRCMean), centralityPercentile);
        FillHistogram("hJetPtBgrdSubtractedKTNegEta", GetCorrectedJetPt(tmpJet, backgroundKTMedian), centralityPercentile);
        FillHistogram("hJetPtBgrdSubtractedTRNegEta", GetCorrectedJetPt(tmpJet, backgroundTRMean), centralityPercentile);
      }

      // Correct EVERY jet with suitable background
      Double_t tmpKTRho = 0.0;
      Double_t tmpRCRho = 0.0;
      Double_t tmpTRRho = 0.0;
      Double_t dummy = 0.0;
      // Calculate backgrounds
      GetKTBackgroundDensity (fNumberExcludedJets, tmpKTRho, dummy, tmpJet->Eta()-0.1, tmpJet->Eta()+0.1, kTRUE);
      GetRCBackgroundDensity (fNumberExcludedJets, tmpRCRho, dummy, tmpJet->Eta()-0.1, tmpJet->Eta()+0.1);
      GetTRBackgroundDensity (fNumberExcludedJets, tmpTRRho, dummy, tmpJet->Eta()-0.1, tmpJet->Eta()+0.1);

      FillHistogram("hKTBackgroundEtaWindow", tmpKTRho, centralityPercentile);
      FillHistogram("hRCBackgroundEtaWindow", tmpRCRho, centralityPercentile);
      FillHistogram("hTRBackgroundEtaWindow", tmpTRRho, centralityPercentile);

      FillHistogram("hJetPtBgrdSubtractedKTEtaWindow", GetCorrectedJetPt(tmpJet, tmpKTRho), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedRCEtaWindow", GetCorrectedJetPt(tmpJet, tmpRCRho), centralityPercentile);
      FillHistogram("hJetPtBgrdSubtractedTREtaWindow", GetCorrectedJetPt(tmpJet, tmpTRRho), centralityPercentile);


      // Signal jet vs. signal jet - "Combinatorial"
      for (Int_t j = i+1; j<fNumberSignalJets; j++)
        FillHistogram("hJetDeltaPhi", GetDeltaPhi(tmpJet->Phi(), fSignalJets[j]->Phi()));
    }

    // ### DIJETS
    if(fNumberSignalJets >= 2)
    {
      FillHistogram("hLeadingJetDeltaPhi", GetDeltaPhi(fFirstLeadingJet->Phi(), fSecondLeadingJet->Phi()));

      if (IsDijet(fFirstLeadingJet, fSecondLeadingJet))
      {
        FillHistogram("hDijetConstituentsPt", fFirstLeadingJet->Pt());
        FillHistogram("hDijetConstituentsPt", fSecondLeadingJet->Pt());

        FillHistogram("hDijetLeadingJetPt", fFirstLeadingJet->Pt());
        FillHistogram("hDijetPtCorrelation", fFirstLeadingJet->Pt(), fSecondLeadingJet->Pt());
        Double_t dummyArea = 0;
        GetTRBackgroundDensity (2, dijetBackground, dummyArea, fFirstLeadingJet, fSecondLeadingJet, kFALSE);
        GetTRBackgroundDensity (2, dijetBackgroundPerpendicular, dummyArea, fFirstLeadingJet, fSecondLeadingJet, kTRUE);
      }
    }

    // ### SOME JET PLOTS
    FillHistogram("hJetCountAll", fJetArray->GetEntries());
    FillHistogram("hJetCountAccepted", fNumberSignalJets);
    if (fFirstLeadingJet)
      FillHistogram("hLeadingJetPt", fFirstLeadingJet->Pt());
    if (fSecondLeadingJet)
      FillHistogram("hSecondLeadingJetPt", fSecondLeadingJet->Pt());

  } //endif AnalyzeJets

  #ifdef DEBUGMODE
    AliInfo("Calculate()::Jets done.");
  #endif
  ////////////////////// NOTE: Background analysis

  if (fAnalyzeBackground)
  {
    // Calculate (non-eta corrected) background in centrality classes
    FillHistogram("hRCBackground", backgroundRCMean, centralityPercentile);
    FillHistogram("hTRBackground", backgroundTRMean, centralityPercentile);
    FillHistogram("hKTBackground", backgroundKTMedian, centralityPercentile);

    // Calculate backgrounds in eta bins
    for(Int_t i=0;i<5;i++)
    {
      Double_t dummy = 0.0;
      Double_t tmpKTRho = 0.0;
      Double_t tmpRCRho = 0.0;
      Double_t tmpTRRho = 0.0;

      Double_t etaMin = -(fTrackEtaWindow-fSignalJetRadius) + 2*(fTrackEtaWindow-fSignalJetRadius)/5 *  i;
      Double_t etaMax = -(fTrackEtaWindow-fSignalJetRadius) + 2*(fTrackEtaWindow-fSignalJetRadius)/5 * (i+1);

      // Calculate backgrounds
      GetKTBackgroundDensity    (fNumberExcludedJets, tmpKTRho, dummy, etaMin, etaMax);
      GetTRBackgroundDensity    (fNumberExcludedJets, tmpTRRho, dummy, etaMin, etaMax);
      GetRCBackgroundDensity    (fNumberExcludedJets, tmpRCRho, dummy, etaMin, etaMax);

      FillHistogram("hRCBackgroundEtaBins", tmpRCRho, (etaMin+etaMax)/2.0);
      FillHistogram("hTRBackgroundEtaBins", tmpTRRho, (etaMin+etaMax)/2.0);
      FillHistogram("hKTBackgroundEtaBins", tmpKTRho, (etaMin+etaMax)/2.0);
    }

    // In case of dijets -> look at the background
    if (dijetBackground >= 0)
      FillHistogram("hDijetBackground", dijetBackground, centralityPercentile); 
    if (dijetBackgroundPerpendicular >= 0)
      FillHistogram("hDijetBackgroundPerpendicular", dijetBackgroundPerpendicular, centralityPercentile); 


    // Calculate the delta pt
    Double_t tmpDeltaPtKT = 0.0;
    Double_t tmpDeltaPtRC = 0.0;
    Double_t tmpDeltaPtTR = 0.0;
    Double_t tmpDeltaPtKTNoExcl = 0.0;
    Double_t tmpDeltaPtRCNoExcl = 0.0;
    Double_t tmpDeltaPtTRNoExcl = 0.0;
    GetDeltaPt(tmpDeltaPtKT, backgroundKTMedian);
    GetDeltaPt(tmpDeltaPtRC, backgroundRCMean);
    GetDeltaPt(tmpDeltaPtTR, backgroundTRMean);
    GetDeltaPt(tmpDeltaPtKTNoExcl, backgroundKTMedian, kFALSE);
    GetDeltaPt(tmpDeltaPtRCNoExcl, backgroundRCMean, kFALSE);
    GetDeltaPt(tmpDeltaPtTRNoExcl, backgroundTRMean, kFALSE);

    // If valid, fill the delta pt histograms
    if(tmpDeltaPtKT > -10000.0)
      FillHistogram("hDeltaPtKT", tmpDeltaPtKT, centralityPercentile);
    if(tmpDeltaPtRC > -10000.0)
      FillHistogram("hDeltaPtRC", tmpDeltaPtRC, centralityPercentile);
    if(tmpDeltaPtTR > -10000.0)
      FillHistogram("hDeltaPtTR", tmpDeltaPtTR, centralityPercentile);
    if(tmpDeltaPtKTNoExcl > -10000.0)
      FillHistogram("hDeltaPtKTNoExcl", tmpDeltaPtKTNoExcl, centralityPercentile);
    if(tmpDeltaPtRCNoExcl > -10000.0)
      FillHistogram("hDeltaPtRCNoExcl", tmpDeltaPtRCNoExcl, centralityPercentile);
    if(tmpDeltaPtTRNoExcl > -10000.0)
      FillHistogram("hDeltaPtTRNoExcl", tmpDeltaPtTRNoExcl, centralityPercentile);

  }
  
  #ifdef DEBUGMODE
    AliInfo("Calculate()::Background done.");
  #endif
  
  ////////////////////// NOTE: Pythia histograms
  if(fAnalyzePythia)
  {
    FillHistogram("hPythiaPtHard", GetPtHard());
    FillHistogram("hPythiaNTrials", GetPtHardBin()-0.1, fTrials);
    FillHistogram("hPythiaXSection", GetPtHardBin()-0.1, fCrossSection);

    #ifdef DEBUGMODE
      AliInfo("Calculate()::Pythia done.");
    #endif
  }
  #ifdef DEBUGMODE
    AliInfo("Calculate() done.");
  #endif
}

//________________________________________________________________________
Bool_t AliAnalysisTaskChargedJetsPA::Notify()
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
inline Double_t AliAnalysisTaskChargedJetsPA::EtaToTheta(Double_t arg)
  {return 2.*atan(exp(-arg));} 
//________________________________________________________________________
inline Double_t AliAnalysisTaskChargedJetsPA::ThetaToEta(Double_t arg)
{
  if ((arg > TMath::Pi()) || (arg < 0.0))
  {
    AliError(Form("ThetaToEta got wrong input! (%f)", arg));
    return 0.0;
  }
  return -log(tan(arg/2.));
}
//________________________________________________________________________
inline Double_t AliAnalysisTaskChargedJetsPA::GetDeltaPhi(Double_t phi1, Double_t phi2)
  {return min(TMath::Abs(phi1-phi2),TMath::TwoPi() - TMath::Abs(phi1-phi2));}

//________________________________________________________________________
Double_t AliAnalysisTaskChargedJetsPA::MCGetOverlapCircleRectancle(Double_t cPosX, Double_t cPosY, Double_t cRadius, Double_t rPosXmin, Double_t rPosXmax, Double_t rPosYmin, Double_t rPosYmax)
{
  const Int_t kTests = 1000;
  Int_t hits = 0;
  TRandom3 randomGen(0);
 
  // Loop over kTests-many tests
  for (Int_t i=0; i<kTests; i++)
  {
    //Choose random position in rectangle for the tester
    Double_t tmpTestX = randomGen.Uniform(rPosXmin, rPosXmax);
    Double_t tmpTestY = randomGen.Uniform(rPosYmin, rPosYmax);

    //Check, if tester is in circle. If yes, increment circle counter.
    Double_t tmpDistance = TMath::Sqrt( (tmpTestX - cPosX)*(tmpTestX - cPosX) + (tmpTestY - cPosY)*(tmpTestY - cPosY) );
    if(tmpDistance < cRadius)
      hits++;
  }

  // return ratio
  return (static_cast<Double_t>(hits)/static_cast<Double_t>(kTests));
}

//________________________________________________________________________
inline void AliAnalysisTaskChargedJetsPA::FillHistogram(const char * key, Double_t x)
{
  TH1* tmpHist = static_cast<TH1*>(fOutputList->FindObject(GetHistoName(key)));
  if(!tmpHist)
  {
    AliWarning(Form("Cannot find histogram <%s> ",key)) ;
    return;
  }

  tmpHist->Fill(x);
}

//________________________________________________________________________
inline void AliAnalysisTaskChargedJetsPA::FillHistogram(const char * key, Double_t x, Double_t y)
{
  TH1* tmpHist = static_cast<TH1*>(fOutputList->FindObject(GetHistoName(key)));
  if(!tmpHist)
  {
    AliWarning(Form("Cannot find histogram <%s> ",key));
    return;
  }

  if (tmpHist->IsA()->GetBaseClass("TH1"))
    static_cast<TH1*>(tmpHist)->Fill(x,y); // Fill x with y
  else if (tmpHist->IsA()->GetBaseClass("TH2"))
    static_cast<TH2*>(tmpHist)->Fill(x,y); // Fill x,y with 1
}

//________________________________________________________________________
inline void AliAnalysisTaskChargedJetsPA::FillHistogram(const char * key, Double_t x, Double_t y, Double_t add)
{
  TH2* tmpHist = static_cast<TH2*>(fOutputList->FindObject(GetHistoName(key)));
  if(!tmpHist)
  {
    AliWarning(Form("Cannot find histogram <%s> ",key));
    return;
  }
  
  tmpHist->Fill(x,y,add);
}
//________________________________________________________________________
template <class T> T* AliAnalysisTaskChargedJetsPA::AddHistogram1D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, const char* xTitle, const char* yTitle)
{
  T* tmpHist = new T(GetHistoName(name), GetHistoName(title), xBins, xMin, xMax);

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
template <class T> T* AliAnalysisTaskChargedJetsPA::AddHistogram2D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, const char* xTitle, const char* yTitle, const char* zTitle)
{
  T* tmpHist = new T(GetHistoName(name), GetHistoName(title), xBins, xMin, xMax, yBins, yMin, yMax);
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
void AliAnalysisTaskChargedJetsPA::Terminate(Option_t *)
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
AliAnalysisTaskChargedJetsPA::~AliAnalysisTaskChargedJetsPA()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputList;
  }
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsPA::UserCreateOutputObjects()
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
void AliAnalysisTaskChargedJetsPA::UserExec(Option_t *) 
{
  #ifdef DEBUGMODE
    AliInfo("UserExec() started.");
  #endif

  if (!InputEvent())
  {
    AliError("??? Event pointer == 0 ???");
    return;
  }

  if (!fInitialized)
    ExecOnce(); // Get tracks, jets, background from arrays if not already given + Init Histos
  
  Calculate(InputEvent());
        
  PostData(1, fOutputList);
}
