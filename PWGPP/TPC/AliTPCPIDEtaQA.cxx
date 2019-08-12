#include "TChain.h"
#include "TTree.h"
#include "TF1.h"
#include "TAxis.h"
#include "TH1F.h"
#include "THnSparse.h"

#include "AliMCParticle.h"
//#include "AliStack.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"

#include "AliESDtrackCuts.h"
#include "AliVVertex.h"
#include "AliAnalysisFilter.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"

#include "AliTPCPIDEtaQA.h"

/*
This task determines the eta dependence of the TPC signal.
For this purpose, only tracks fulfilling some standard quality cuts are taken into account.
The obtained data can be used to derive the functional behaviour of the eta dependence.
Such a function can be plugged into this task to correct for the eta dependence and to see
if there is then really no eta dependence left.

Class written by Benjamin Hess.
Contact: bhess@cern.ch
*/

ClassImp(AliTPCPIDEtaQA)

//________________________________________________________________________
AliTPCPIDEtaQA::AliTPCPIDEtaQA()
  : AliTPCPIDBase()
  , fPtThresholdForPhiCut(2.0)
  , fPhiCutSecondBranchLow(0x0)
  , fPhiCutSecondBranchHigh(0x0)
  , fCentralityEstimator("")
  , fhPIDdataAll(0x0)
  //OLD clusterQA, fhNumClustersPhiPrimePtBeforeCut(0x0)
  //OLD clusterQA, fhNumClustersPhiPrimePtAfterCut(0x0)
  , fhPhiPrimeCutEfficiency(0x0)
  , fOutputContainer(0x0)
{
  // default Constructor

  // Question: Is this the right place to initialize these functions?
  // Will it work on proof? i.e. will they be streamed to the workers?
  // Also one should add getters and setters
  //TODO fPhiCutLow  = new TF1("StandardPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 50);
  //TODO fPhiCutHigh = new TF1("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 50);
  
  //TODO NEW fPhiCutLow  = new TF1("StandardPhiCutLow",  "0.072/x+pi/18.0-0.035", 0, 50);
  //TODO NEW fPhiCutHigh = new TF1("StandardPhiCutHigh", "0.07/x/x+0.1/x+pi/18.0+0.035", 0, 50);
  
  fPhiCutSecondBranchLow  = new TF1("StandardPhiCutSecondBranchLow",  "NewStandardPhiCutLow - 2.*pi/18.", 0, 50);
  fPhiCutSecondBranchHigh = new TF1("StandardPhiCutSecondBranchHigh", "0.07/x+pi/18.0+0.125 - 2.*pi/18.", 0, 50);
}

//________________________________________________________________________
AliTPCPIDEtaQA::AliTPCPIDEtaQA(const char *name)
  : AliTPCPIDBase(name)
  , fPtThresholdForPhiCut(2.0)
  , fPhiCutSecondBranchLow(0x0)
  , fPhiCutSecondBranchHigh(0x0)
  , fCentralityEstimator("")
  , fhPIDdataAll(0x0)
  //OLD clusterQA, fhNumClustersPhiPrimePtBeforeCut(0x0)
  //OLD clusterQA, fhNumClustersPhiPrimePtAfterCut(0x0)
  , fhPhiPrimeCutEfficiency(0x0)
  , fOutputContainer(0x0)
{
  // Constructor
  
  // Question: Is this the right place to initialize these functions?
  // Will it work on proof? i.e. will they be streamed to the workers?
  // Also one should add getters and setters
  //TODO fPhiCutLow  = new TF1("StandardPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 50);
  //TODO fPhiCutHigh = new TF1("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 50);
  
  //TODO NEW fPhiCutLow  = new TF1("StandardPhiCutLow",  "0.072/x+pi/18.0-0.035", 0, 50);
  //TODO NEW fPhiCutHigh = new TF1("StandardPhiCutHigh", "0.07/x/x+0.1/x+pi/18.0+0.035", 0, 50);
  
  fPhiCutSecondBranchLow  = new TF1("StandardPhiCutSecondBranchLow",  "StandardPhiCutLow - 2.*pi/18.", 0, 50);
  fPhiCutSecondBranchHigh = new TF1("StandardPhiCutSecondBranchHigh", "0.07/x+pi/18.0+0.125 - 2.*pi/18.", 0, 50);

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1, TObjArray::Class());
}


//________________________________________________________________________
AliTPCPIDEtaQA::~AliTPCPIDEtaQA()
{
  // dtor
  
  delete fOutputContainer;
  fOutputContainer = 0;
  
  delete fPhiCutSecondBranchLow;
  fPhiCutSecondBranchLow = 0;
  
  delete fPhiCutSecondBranchHigh;
  fPhiCutSecondBranchHigh = 0;
}


//________________________________________________________________________
void AliTPCPIDEtaQA::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  AliTPCPIDBase::UserCreateOutputObjects();

  OpenFile(1);
  
  fOutputContainer = new TObjArray(1);
  fOutputContainer->SetName(GetName()) ;
  fOutputContainer->SetOwner(kTRUE);
  
  const Int_t nEta = TMath::Nint(40 * fEtaCut);
  
  const Int_t nPtBins = 68;
  Double_t binsPt[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
           0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
           1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
           2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
           4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
           11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
           26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0 };
    
      
  const Int_t nBins = 6;
  // MC PID, SelectSpecies, P(TPC_inner), multiplicity, deltaPrimeSpecies, eta
  
  
  // Binning for multiplicity
  Int_t nMultBins = 40;
  Double_t multLow = 0.;
  Double_t multUp = 20000;
  
  CentralityEstimatorType centType = GetCentralityEstimatorType();
  if (centType == kITSTPCtracklets) {
    // Special pp centrality estimator
    nMultBins = 100;
    multLow = 0;
    multUp = 100;
  }
  else if (centType == kPPmultV0M) {
    // Another special pp centrality estimator
    nMultBins = 100;
    multLow = 0.;
    multUp = 1.;
  }
  
  
  Int_t bins[nBins] = 
    {  9,   4, nPtBins, nMultBins, 201, nEta };
  Double_t xmin[nBins] = 
    {  0., 0.,      0.,   multLow, 0.5, -fEtaCut};
  Double_t xmax[nBins] = 
    {  9., 4.,    50.0,    multUp, 2.0, fEtaCut };
 
  		
  fhPIDdataAll = new THnSparseI("hPIDdataAll","", nBins, bins, xmin, xmax);
  SetUpHist(fhPIDdataAll, binsPt);
  fOutputContainer->Add(fhPIDdataAll);
  
  /*OLD clusterQA
  const Int_t nBinsQA = 3;
  // pT, phiPrime, Ncl
  
  Int_t binsQA[nBinsQA] = 
    {  nPtBins, 50, 90 }; 
  Double_t xminQA[nBinsQA] = 
    {  0., 0., 70.};
  Double_t xmaxQA[nBinsQA] = 
    {  50.0, TMath::Pi() / 9., 160 };
    
  fhNumClustersPhiPrimePtBeforeCut = new THnSparseI("hNumClustersPhiPrimePtBeforeCut", "", nBinsQA, binsQA, xminQA, xmaxQA);
  fhNumClustersPhiPrimePtBeforeCut->SetBinEdges(0, binsPt);
  fhNumClustersPhiPrimePtBeforeCut->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fhNumClustersPhiPrimePtBeforeCut->GetAxis(1)->SetTitle("#phi'");
  fhNumClustersPhiPrimePtBeforeCut->GetAxis(2)->SetTitle("Ncl");
  fOutputContainer->Add(fhNumClustersPhiPrimePtBeforeCut);
  
  fhNumClustersPhiPrimePtAfterCut = new THnSparseI("hNumClustersPhiPrimePtAfterCut", "", nBinsQA, binsQA, xminQA, xmaxQA);
  fhNumClustersPhiPrimePtAfterCut->SetBinEdges(0, binsPt);
  fhNumClustersPhiPrimePtAfterCut->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fhNumClustersPhiPrimePtAfterCut->GetAxis(1)->SetTitle("#phi'");
  fhNumClustersPhiPrimePtAfterCut->GetAxis(2)->SetTitle("Ncl");
  fOutputContainer->Add(fhNumClustersPhiPrimePtAfterCut);
  
  fhPhiPrimeCutEfficiency = new TH1F("hPhiPrimeCutEfficiency", "Efficiency of #phi' cut;p_{T} (GeV/c);Cut efficiency", nPtBins, binsPt);
  fOutputContainer->Add(fhPhiPrimeCutEfficiency);
  */
  
  PostData(1, fOutputContainer);
}


//________________________________________________________________________
void AliTPCPIDEtaQA::UserExec(Option_t *)
{
  // Main loop
  // Called for each event

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  
  if (GetIsPileUp(fESD, kPileUpRejectionSPD))
    return;
  
  if (!fPIDResponse || !fV0KineCuts)
    return;
  
  if (!GetVertexIsOk(fESD))
    return;
  
  const AliVVertex* primaryVertex = fESD->GetPrimaryVertexTracks(); 
  if (!primaryVertex)
    return;
    
  if(primaryVertex->GetNContributors() <= 0) 
    return;
  
  fMC = dynamic_cast<AliMCEvent*>(MCEvent());
  
  Double_t magField = fESD->GetMagneticField();
  
  // Fill V0 arrays with V0s
  FillV0PIDlist();
  
  Double_t multiplicity = -1;
  
  CentralityEstimatorType centType = GetCentralityEstimatorType();
  if (centType == kITSTPCtracklets) {
    // Special pp centrality estimator
    multiplicity = AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);
  }
  else if (centType == kPPmultV0M) {
    // Another special pp centrality estimator
    multiplicity = fAnaUtils->GetMultiplicityPercentile(fESD, "V0M");
  }
  else {
    // Take default, namely number of ESD tracks
    multiplicity = fESD->GetNumberOfESDTracks();
  }
  
  // Track loop to fill a Train spectrum
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    
    if (TMath::Abs(track->Eta()) >= fEtaCut)  continue;
    
    //if (TMath::Abs(track->Eta()) < 0.6 || TMath::Abs(track->Eta()) > 0.8)  continue;
    
    // Do not distinguish positively and negatively charged V0's
    Char_t v0tagAllCharges = TMath::Abs(GetV0tag(iTracks));
    if (v0tagAllCharges == -99) {
      AliError(Form("Problem with V0 tag list (requested V0 track for track %d from %d (list status %d))!", iTracks, fESD->GetNumberOfTracks(),
                    fV0tags != 0x0));
      continue;
    }
    
    Bool_t isV0el = v0tagAllCharges == 14;
    Bool_t isV0pi = v0tagAllCharges == 15;
    Bool_t isV0pr = v0tagAllCharges == 16;
    Bool_t isV0 = isV0el || isV0pi || isV0pr;
    
    // Apply track cut. Accept track, if from V0 or if accepted by cut.
    // Do not apply the track cuts to V0s, since the track cuts are usually
    // cuts on primaries, what will throw away almost all V0s.
    if(!isV0 && (fTrackFilter && !fTrackFilter->IsSelected(track)))
      continue;
    
    if (GetUseTPCCutMIGeo()) {
      // If cut on geometry is active, don't cut on number of clusters, since such a cut is already included.
      if (!isV0) {
        if (!TPCCutMIGeo(track, InputEvent()))
          continue;
      }
      else {
        // One should not cut on geometry for V0's since they have different topology. (Loosely) cut on num of clusters instead.
        if (track->GetTPCsignalN() < 60)
          continue;
      }
    }
    else {
      // If cut on geometry is not active, always cut on num clusters
      if (track->GetTPCsignalN() < 60)
        continue;
    }
    
    Double_t pT = track->Pt();
    //Double_t phiPrime = GetPhiPrime(track->Phi(), magField, track->Charge());
    
    //OLD clusterQA Double_t entryQA[3] = { pT, phiPrime, track->GetTPCsignalN() };
    //OLD clusterQA fhNumClustersPhiPrimePtBeforeCut->Fill(entryQA);

    // Apply PhiPrimeCut only on high-pt tracks, in this case: Tracks with pt >= fPtThresholdForPhiCut
    if(fUsePhiCut && pT >= fPtThresholdForPhiCut) {
      
      if(fUsePhiCut) {
      if (!PhiPrimeCut(track, magField))
        continue; // reject track
    }
      //Use cut for second branch?
      if(!PhiPrimeCut(track, magField))
      //if(!PhiPrimeCut(track, magField)
      //   || (phiPrime < fPhiCutSecondBranchHigh->Eval(pT)  && phiPrime > fPhiCutSecondBranchLow->Eval(pT)))
	      continue; // reject track
    }
    
    //OLD clusterQA fhNumClustersPhiPrimePtAfterCut->Fill(entryQA);
    
    Int_t label = track->GetLabel();

    AliMCParticle* mcTrack = 0x0;
    
    if (fMC) {
      if (label < 0)
        continue; // If MC is available, reject tracks with negative ESD label
      mcTrack = dynamic_cast<AliMCParticle*>(fMC->GetTrack(TMath::Abs(label)));
      if (!mcTrack) {
        Printf("ERROR: Could not receive mcTrack with label %d for track %d", label, iTracks);
        continue;
      }
      
      /*// Only accept MC primaries
      if (!fMC->Stack()->IsPhysicalPrimary(TMath::Abs(label))) {
        continue;
      }*/
    }
    
    // Momentum
    Double_t pTPC = track->GetTPCmomentum();
    
    // Eta
    Float_t eta = track->Eta();
    
    // TPC signal
    Double_t dEdxTPC = track->GetTPCsignal();
    
    Double_t dEdxExpectedEl = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault, 
                                                                              fPIDResponse->UseTPCEtaCorrection(),
                                                                              fPIDResponse->UseTPCMultiplicityCorrection());
    Double_t dEdxExpectedKa = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kKaon, AliTPCPIDResponse::kdEdxDefault, 
                                                                              fPIDResponse->UseTPCEtaCorrection(),
                                                                              fPIDResponse->UseTPCMultiplicityCorrection());
    Double_t dEdxExpectedPi = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kPion, AliTPCPIDResponse::kdEdxDefault, 
                                                                              fPIDResponse->UseTPCEtaCorrection(),
                                                                              fPIDResponse->UseTPCMultiplicityCorrection());
    Double_t dEdxExpectedPr = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault, 
                                                                              fPIDResponse->UseTPCEtaCorrection(),
                                                                              fPIDResponse->UseTPCMultiplicityCorrection());
    
    Double_t deltaPrimeEl = (dEdxExpectedEl > 0) ? (dEdxTPC / dEdxExpectedEl) : (-1);
    Double_t deltaPrimeKa = (dEdxExpectedKa > 0) ? (dEdxTPC / dEdxExpectedKa) : (-1);
    Double_t deltaPrimePi = (dEdxExpectedPi > 0) ? (dEdxTPC / dEdxExpectedPi) : (-1);
    Double_t deltaPrimePr = (dEdxExpectedPr > 0) ? (dEdxTPC / dEdxExpectedPr) : (-1);
    
    Int_t binMC = -1;
    if (fMC) {
      Int_t pdg = mcTrack->PdgCode();
      
      if (TMath::Abs(pdg) == 211) {//Pion
        // If V0, only accept if correctly identified
        if (isV0)   {
          if (isV0pi)
            binMC = 7;
          else
            continue;
        }
        else
          binMC = 3;
      }
      else if (TMath::Abs(pdg) == 321) {//Kaon
        // No kaons from V0 => Reject
        if (isV0)
          continue;
        else
          binMC = 1;
      }
      else if (TMath::Abs(pdg) == 2212) {//Proton
        // If V0, only accept if correctly identified
        if (isV0)   {
          if (isV0pr)
            binMC = 8;
          else
            continue;
        }
        else
          binMC = 4;
      }
      else if (TMath::Abs(pdg) == 11) {//Electron
        // If V0, only accept if correctly identified
        if (isV0)   {
          if (isV0el)
            binMC = 6;
          else
            continue;
        }
        else
          binMC = 0;
      }
      else if (TMath::Abs(pdg) == 13) {//Muon
        // No muons from V0 => Reject
        if (isV0)
          continue;
        else
          binMC = 2;
      }
      else
        continue; // Reject all other particles
    }
    
    Bool_t isKaon = kFALSE;
    Bool_t isPion = kFALSE;
    Bool_t isProton = kFALSE;
    Bool_t isElectron = kFALSE;
    Bool_t isV0plusTOFel = kFALSE;
    
    
    if (!fMC && !isV0) {    
      UInt_t status = track->GetStatus();
      Bool_t hasTOFout  = status&AliESDtrack::kTOFout; 
      Bool_t hasTOFtime = status&AliESDtrack::kTIME;
      Bool_t hasTOF     = kFALSE;
      if (hasTOFout && hasTOFtime)
        hasTOF = kTRUE;
      Float_t length = track->GetIntegratedLength();
      // Length check only for primaries, not for V0's!
      if (length < 350 && !isV0)
        hasTOF = kFALSE;

    
      // Select Kaons, Pions and Protons in 3 sigma band (TOF and TPC). Below some momentum threshold, only use TPC cut,
      // since TOF efficiency is really bad.
      // In case of ambiguity, add entries for corresponding species
      if ((pTPC <= 0.25 && TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)) < 6.0) ||
          (pTPC > 0.25 && pTPC <= 0.3 && TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)) < 4.0) ||
          (pTPC > 0.3 && pTPC <= 0.35 && TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)) < 3.0) ||
          (pTPC > 0.35 && TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)) < 3.0 && hasTOF && 
           TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon)) < 3.0)) {
          isKaon = kTRUE;
      }
      if ((dEdxTPC >= 50. / TMath::Power(pTPC, 1.3)) && // Pattern recognition instead of TPC cut to be ~independent of old TPC expected dEdx
         ((pTPC < 0.6) ||
         (pTPC >= 0.6 && hasTOF && TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton)) < 3.0))) {
        isProton = kTRUE;
      }
      /*
      if ((TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton)) < (pTPC < 0.3 ? 8.0 : 5.0)) &&
          ((pTPC < 0.6) ||
          (pTPC >= 0.6 && hasTOF && TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton)) < 3.0))) {
          isProton = kTRUE;
      }*/
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)) < 3.0 &&
          (pTPC <= 0.3 || 
           (pTPC > 0.3 && hasTOF && TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion)) < 3.0)))  {
          isPion = kTRUE;
      }
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron)) < 3.0 &&
          (pTPC <= 0.3 ||
           (pTPC > 0.3 && hasTOF && TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron)) < 3.0)))  {
          isElectron = kTRUE;
      }
    }
    else if (!fMC && isV0el) {
      // Special treatment for V0 electrons -> Look for V0+TOF electrons
      
      UInt_t status = track->GetStatus();
      Bool_t hasTOFout  = status&AliESDtrack::kTOFout; 
      Bool_t hasTOFtime = status&AliESDtrack::kTIME;
      Bool_t hasTOF     = kFALSE;
      if (hasTOFout && hasTOFtime)
        hasTOF = kTRUE;
      
      // No length check for V0's
      if (hasTOF && TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron)) < 3.0) {
        isV0plusTOFel = kTRUE;
      }      
    }
    
    // MC PID, SelectSpecies, P(TPC_inner), multiplicity, deltaPrimeSpecies, eta
    Double_t entry[6] = { (Double_t)binMC, 0., pTPC, multiplicity, deltaPrimeEl, eta };
    
    if (fMC)  {
      fhPIDdataAll->Fill(entry);
      
      entry[1] = 1;
      entry[4] = deltaPrimeKa;
      fhPIDdataAll->Fill(entry);
      
      entry[1] = 2;
      entry[4] = deltaPrimePi;
      fhPIDdataAll->Fill(entry);
      
      entry[1] = 3;
      entry[4] = deltaPrimePr;
      fhPIDdataAll->Fill(entry);
    }
    else  {      
      for (Int_t i = 0; i < 8; i++) {
        if (i == 0 && isKaon) 
          binMC = 1;
        else if (i == 1 && isPion)  
          binMC = 3;
        else if (i == 2 && isProton)  
          binMC = 4;
        else if (i == 3 && isElectron)
          binMC = 0;
        else if (i == 4 && isV0plusTOFel)
          binMC = 5;
        else if (i == 5 && isV0el)
          binMC = 6;
        else if (i == 6 && isV0pi)
          binMC = 7;
        else if (i == 7 && isV0pr)
          binMC = 8;
        else
          continue;
          
        entry[0] = binMC;
        
        entry[1] = 0;
        entry[4] = deltaPrimeEl;
        fhPIDdataAll->Fill(entry);
        
        entry[1] = 1;
        entry[4] = deltaPrimeKa;
        fhPIDdataAll->Fill(entry);
        
        entry[1] = 2;
        entry[4] = deltaPrimePi;
        fhPIDdataAll->Fill(entry);
        
        entry[1] = 3;
        entry[4] = deltaPrimePr;
        fhPIDdataAll->Fill(entry);
      }
    }
  } //track loop 

  // Post output data.
  PostData(1, fOutputContainer);
  
  // Clear the V0 PID arrays
  ClearV0PIDlist();
}      

//________________________________________________________________________
void AliTPCPIDEtaQA::Terminate(const Option_t *)
{
  // Called once at the end of the query
  /*OLD clusterQA
  TObjArray* output = (TObjArray*)GetOutputData(1);
  if (!output)
    return;
  
  TH1F* hPhiPrimeCutEfficiency = (TH1F*)(output->FindObject("hPhiPrimeCutEfficiency"));
  if (!hPhiPrimeCutEfficiency) 
    return;
  
  THnSparse* hNumClustersPhiPrimePtBeforeCut = (THnSparse*)(output->FindObject("hNumClustersPhiPrimePtBeforeCut"));
  if (!hNumClustersPhiPrimePtBeforeCut)
    return;
  
  THnSparse* hNumClustersPhiPrimePtAfterCut = (THnSparse*)(output->FindObject("hNumClustersPhiPrimePtAfterCut"));
  if (!hNumClustersPhiPrimePtAfterCut)
    return;
  
  TH1D* hNumEntriesBefore = hNumClustersPhiPrimePtBeforeCut->Projection(0);
  TH1D* hNumEntriesAfter = hNumClustersPhiPrimePtAfterCut->Projection(0);
   
  for (Int_t bin = 1; bin <= hNumEntriesBefore->GetXaxis()->GetNbins(); bin++)  {
    Double_t numEntriesBefore = hNumEntriesBefore->GetBinContent(bin);
    
    if (numEntriesBefore > 0) {
      Double_t numEntriesAfter = hNumEntriesAfter->GetBinContent(bin);
      hPhiPrimeCutEfficiency->SetBinContent(bin,  numEntriesAfter / numEntriesBefore);
      // Errors are 100%-correlated: Accepted tracks should be a constant fraction of all tracks in given bin.
      // Thus: Only take statistical fluctuation from accepted tracks as error
      hPhiPrimeCutEfficiency->SetBinError(bin, TMath::Sqrt(numEntriesAfter) / numEntriesBefore);
      //                          (numEntriesAfter / numEntriesBefore) * TMath::Sqrt(1. / numEntriesAfter + 1. / numEntriesBefore));
    }
  }
  
  hPhiPrimeCutEfficiency->SetDrawOption("e");
  */
}


//________________________________________________________________________
void AliTPCPIDEtaQA::SetUpHist(THnSparse* hist, Double_t* binsPt) const
{
  // Sets bin limits for axes which are not standard binned and the axes titles.
  // MC PID, SelectSpecies, P(TPC_inner), multiplicity, deltaPrimeSpecies, eta
  hist->SetBinEdges(2, binsPt);
                          
  // Set axes titles
  hist->GetAxis(0)->SetTitle("MC PID");
  hist->GetAxis(0)->SetBinLabel(1, "e");
  hist->GetAxis(0)->SetBinLabel(2, "K");
  hist->GetAxis(0)->SetBinLabel(3, "#mu");
  hist->GetAxis(0)->SetBinLabel(4, "#pi");
  hist->GetAxis(0)->SetBinLabel(5, "p");
  hist->GetAxis(0)->SetBinLabel(6, "V0+TOF e");
  hist->GetAxis(0)->SetBinLabel(7, "V0 e");
  hist->GetAxis(0)->SetBinLabel(8, "V0 #pi");
  hist->GetAxis(0)->SetBinLabel(9, "V0 p");
  
  hist->GetAxis(1)->SetTitle("Select Species");
  hist->GetAxis(1)->SetBinLabel(1, "e");
  hist->GetAxis(1)->SetBinLabel(2, "K");
  hist->GetAxis(1)->SetBinLabel(3, "#pi");
  hist->GetAxis(1)->SetBinLabel(4, "p");
  
  hist->GetAxis(2)->SetTitle("p_{TPC_inner} (GeV/c)");
  
  hist->GetAxis(3)->SetTitle("Event multiplicity");
  
  hist->GetAxis(4)->SetTitle("TPC #Delta'{species}");
  
  hist->GetAxis(5)->SetTitle("#eta");
     
  //hist->Sumw2();
}
