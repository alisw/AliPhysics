/* $Id$ */

#include "AliMultiplicityTask.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <TObjString.h>
#include <TF1.h>

#include <AliLog.h>
#include <AliESDVertex.h>
#include <AliESDEvent.h>
#include <AliStack.h>
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include <AliMultiplicity.h>
#include <AliAnalysisManager.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliESDInputHandler.h>

#include "AliESDtrackCuts.h"
#include "AliPWG0Helper.h"
#include "multiplicity/AliMultiplicityCorrection.h"
#include "AliCorrection.h"
#include "AliCorrectionMatrix3D.h"
#include "AliPhysicsSelection.h"
#include "AliTriggerAnalysis.h"

ClassImp(AliMultiplicityTask)

AliMultiplicityTask::AliMultiplicityTask(const char* opt) :
  AliAnalysisTask("AliMultiplicityTask", ""),
  fESD(0),
  fOption(opt),
  fAnalysisMode((AliPWG0Helper::AnalysisMode) (AliPWG0Helper::kSPD | AliPWG0Helper::kFieldOn)),
  fTrigger(AliTriggerAnalysis::kMB1),
  fDeltaPhiCut(-1),
  fDiffTreatment(AliPWG0Helper::kMCFlags),
  fReadMC(kFALSE),
  fUseMCVertex(kFALSE),
  fMultiplicity(0),
  fEsdTrackCuts(0),
  fSystSkipParticles(kFALSE),
  fSelectProcessType(0),
  fParticleSpecies(0),
  fdNdpT(0),
  fPtSpectrum(0),
  fTemp1(0),
  fTemp2(0),
  fVertex(0),
  fEtaPhi(0),
  fOutput(0)
{
  //
  // Constructor. Initialization of pointers
  //

  for (Int_t i = 0; i<8; i++)
    fParticleCorrection[i] = 0;
    
  for (Int_t i=0; i<3; i++)
    fEta[i] = 0;

  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());

  if (fOption.Contains("only-process-type-nd"))
    fSelectProcessType = 1;

  if (fOption.Contains("only-process-type-sd"))
    fSelectProcessType = 2;

  if (fOption.Contains("only-process-type-dd"))
    fSelectProcessType = 3;

  if (fSelectProcessType != 0)
    AliInfo(Form("WARNING: Systematic study enabled. Only considering process type %d", fSelectProcessType));
}

AliMultiplicityTask::~AliMultiplicityTask()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

//________________________________________________________________________
void AliMultiplicityTask::ConnectInputData(Option_t *)
{
  // Connect ESD
  // Called once

  Printf("AliMultiplicityTask::ConnectInputData called");

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read tree from input slot 0");
  } else {
    // Disable all branches and enable only the needed ones
    /*
    tree->SetBranchStatus("*", 0);

    tree->SetBranchStatus("AliESDHeader*", 1);
    tree->SetBranchStatus("*Vertex*", 1);

    if (fAnalysisMode & AliPWG0Helper::kSPD) {
      tree->SetBranchStatus("AliMultiplicity*", 1);
    }

    if (fAnalysisMode & AliPWG0Helper::kTPC || fAnalysisMode & AliPWG0Helper::kTPCITS) {
      //AliESDtrackCuts::EnableNeededBranches(tree);
      tree->SetBranchStatus("Tracks*", 1);
    }
    */

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }

  // disable info messages of AliMCEvent (per event)
  AliLog::SetClassDebugLevel("AliMCEvent", AliLog::kWarning - AliLog::kDebug + 1);
}

void AliMultiplicityTask::CreateOutputObjects()
{
  // create result objects and add to output list

  fOutput = new TList;
  fOutput->SetOwner();

  fMultiplicity = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");
  fOutput->Add(fMultiplicity);
  
  fdNdpT = new TH1F("fdNdpT", "fdNdpT", 1000, 0, 10);
  fdNdpT->Sumw2();
  fOutput->Add(fdNdpT);

  if (fOption.Contains("particle-efficiency"))
    for (Int_t i = 0; i<8; i++)
    {
      fParticleCorrection[i] = new AliCorrection(Form("correction_%d", i), Form("correction_%d", i));
      fOutput->Add(fParticleCorrection[i]);
    }

  if (fOption.Contains("pt-spectrum-hist"))
  {
    TFile* file = TFile::Open("ptspectrum_fit.root");
    if (file)
    {
      TString subStr(fOption(fOption.Index("pt-spectrum")+17, 3));
      TString histName(Form("ptspectrum_%s", subStr.Data()));
      AliInfo(Form("Pt-Spectrum modification. Using %s.", histName.Data()));
      fPtSpectrum = (TH1D*) file->Get(histName);
      if (!fPtSpectrum)   
        AliError("Histogram not found");
    }
    else
      AliError("Could not open ptspectrum_fit.root. Pt Spectrum study could not be enabled.");
  }

  if (fOption.Contains("pt-spectrum-func"))
  {
    if (fPtSpectrum)
    {
      Printf("Using function for systematic p_t study");
    }
    else
    {
      Printf("ERROR: Could not find function for systematic p_t study");
      fPtSpectrum = new TH1D("fPtSpectrum", "fPtSpectrum", 1, 0, 100);
      fPtSpectrum->SetBinContent(1, 1);
    }
  }

  if (fPtSpectrum)
    Printf("WARNING: Systematic study enabled. Pt spectrum will be modified");
  
  if (fOption.Contains("particle-species")) {
    fParticleSpecies = new TNtuple("fParticleSpecies", "fParticleSpecies", "vtx:Pi_True:K_True:p_True:o_True:Pi_Rec:K_Rec:p_Rec:o_Rec:nolabel:doublePrim:doubleCount");
    fOutput->Add(fParticleSpecies);
  }

  fTemp1 = new TH2F("fTemp1", "fTemp1", 100, -0.5, 99.5, 100, -0.5, 99.5);
  fOutput->Add(fTemp1);
  
  for (Int_t i=0; i<3; i++)
  {
    fEta[i] = new TH1F(Form("fEta_%d", i), ";#eta", 100, -2, 2);
    fOutput->Add(fEta[i]);
  }
  
  fVertex = new TH3F("vertex_check", "vertex_check;x;y;z", 100, -1, 1, 100, -1, 1, 100, -30, 30);
  fOutput->Add(fVertex);
  
  fEtaPhi = new TH2F("fEtaPhi", "fEtaPhi;#eta;#phi in rad.;count", 80, -4, 4, 18*20, 0, 2 * TMath::Pi());
  fOutput->Add(fEtaPhi);
  
  // TODO set seed for random generator
}

void AliMultiplicityTask::Exec(Option_t*)
{
  // process the event

  // Check prerequisites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD not available");
    return;
  }

  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }
    
  Bool_t eventTriggered = inputHandler->IsEventSelected();

  static AliTriggerAnalysis* triggerAnalysis = 0;
  if (!triggerAnalysis)
  {
    AliPhysicsSelection* physicsSelection = dynamic_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    triggerAnalysis = physicsSelection->GetTriggerAnalysis();
  }
  if (eventTriggered)
    eventTriggered = triggerAnalysis->IsTriggerFired(fESD, fTrigger);
    
  const AliESDVertex* vtxESD = AliPWG0Helper::GetVertex(fESD, fAnalysisMode);
  if (vtxESD && !AliPWG0Helper::TestVertex(vtxESD, fAnalysisMode))
    vtxESD = 0;
    
  // remove vertices outside +- 15 cm
  if (vtxESD && TMath::Abs(vtxESD->GetZv()) > 15)
    vtxESD = 0;
  
  Bool_t eventVertex = (vtxESD != 0);

  Double_t vtx[3];
  if (vtxESD)
  {
    vtxESD->GetXYZ(vtx);
    fVertex->Fill(vtxESD->GetXv(), vtxESD->GetYv(), vtxESD->GetZv());
  }
  
  // post the data already here
  PostData(0, fOutput);
  
  //const Float_t kPtCut = 0.3;

  // create list of (label, eta) tuples
  Int_t inputCount = 0;
  Int_t* labelArr = 0;
  Float_t* etaArr = 0;
  if (fAnalysisMode & AliPWG0Helper::kSPD)
  {
    if (vtxESD)
    {
      // get tracklets
      const AliMultiplicity* mult = fESD->GetMultiplicity();
      if (!mult)
      {
        AliDebug(AliLog::kError, "AliMultiplicity not available");
        return;
      }
  
      labelArr = new Int_t[mult->GetNumberOfTracklets()];
      etaArr = new Float_t[mult->GetNumberOfTracklets()];
      
      Bool_t foundInEta10 = kFALSE;
      
      // get multiplicity from ITS tracklets
      for (Int_t i=0; i<mult->GetNumberOfTracklets(); ++i)
      {
        //printf("%d %f %f %f\n", i, mult->GetTheta(i), mult->GetPhi(i), mult->GetDeltaPhi(i));
  
        Float_t deltaPhi = mult->GetDeltaPhi(i);
        
        if (fDeltaPhiCut > 0 && TMath::Abs(deltaPhi) > fDeltaPhiCut)
          continue;
  
        if (fSystSkipParticles && gRandom->Uniform() < (0.0153))
        {
          Printf("Skipped tracklet!");
          continue;
        }
          
        etaArr[inputCount] = mult->GetEta(i);
        if (mult->GetLabel(i, 0) == mult->GetLabel(i, 1))
        {
          labelArr[inputCount] = mult->GetLabel(i, 0);
        }
        else
          labelArr[inputCount] = -1;
          
        for (Int_t j=0; j<3; j++)
        {
          if (vtx[2] > fMultiplicity->GetVertexBegin(j) && vtx[2] < fMultiplicity->GetVertexEnd(j))
            fEta[j]->Fill(etaArr[inputCount]);
        }
        
        if (vtxESD && TMath::Abs(vtx[2]) < 10)
          fEtaPhi->Fill(etaArr[inputCount], mult->GetPhi(i));
        
          // we have to repeat the trigger here, because the tracklet might have been kicked out fSystSkipParticles
        if (TMath::Abs(etaArr[inputCount]) < 1)
          foundInEta10 = kTRUE;
          
        ++inputCount;
      }
      
      if (fSystSkipParticles && (fTrigger & AliTriggerAnalysis::kOneParticle) && !foundInEta10)
        eventTriggered = kFALSE;
    }
  }
  else if (fAnalysisMode & AliPWG0Helper::kTPC || fAnalysisMode & AliPWG0Helper::kTPCITS)
  {
    if (!fEsdTrackCuts)
    {
      AliDebug(AliLog::kError, "fESDTrackCuts not available");
      return;
    }

    if (vtxESD)
    {
      // get multiplicity from ESD tracks
      TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD, (fAnalysisMode & AliPWG0Helper::kTPC));
      Int_t nGoodTracks = list->GetEntries();
  
      labelArr = new Int_t[nGoodTracks];
      etaArr = new Float_t[nGoodTracks];
  
      // loop over esd tracks
      for (Int_t i=0; i<nGoodTracks; i++)
      {
        AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (list->At(i));
        if (!esdTrack)
        {
          AliDebug(AliLog::kError, Form("ERROR: Could not retrieve track %d.", i));
          continue;
        }
        
        if (esdTrack->Pt() < 0.15)
          continue;
        
        Float_t d0z0[2],covd0z0[3];
        esdTrack->GetImpactParameters(d0z0,covd0z0);
        Float_t sigma= 0.0050+0.0060/TMath::Power(esdTrack->Pt(),0.9);
        Float_t d0max = 7.*sigma;
        if (TMath::Abs(d0z0[0]) > d0max) 
          continue;
  
        if (vtxESD && TMath::Abs(vtx[2]) < 10)
          fEtaPhi->Fill(esdTrack->Eta(), esdTrack->Phi());
        
        etaArr[inputCount] = esdTrack->Eta();
        labelArr[inputCount] = TMath::Abs(esdTrack->GetLabel());
        ++inputCount;
      }
  
      delete list;
    }
  }

  // eta range for nMCTracksSpecies, nESDTracksSpecies
  Float_t etaRange = 1.0;
//   switch (fAnalysisMode) {
//     case AliPWG0Helper::kInvalid: break;
//     case AliPWG0Helper::kTPC:
//     case AliPWG0Helper::kTPCITS:
//     	etaRange = 1.0; break;
//     case AliPWG0Helper::kSPD: etaRange = 1.0; break;
//     default: break;
//   }

  if (!fReadMC) // Processing of ESD information
  {
    Int_t nESDTracks05 = 0;
    Int_t nESDTracks10 = 0;
    Int_t nESDTracks14 = 0;
    
    for (Int_t i=0; i<inputCount; ++i)
    {
      Float_t eta = etaArr[i];

      if (TMath::Abs(eta) < 0.5)
        nESDTracks05++;

      if (TMath::Abs(eta) < 1.0)
        nESDTracks10++;

      if (TMath::Abs(eta) < 1.3)
        nESDTracks14++;
    }
    
    //if (nESDTracks05 >= 20 || nESDTracks10 >= 30 || nESDTracks14 >= 32)
    //  Printf("File: %s, IEV: %d, TRG: ---, Orbit: 0x%x, Period: %d, BC: %d; Tracks: %d %d %d", ((TTree*) GetInputData(0))->GetCurrentFile()->GetName(), fESD->GetEventNumberInFile(), fESD->GetOrbitNumber(),fESD->GetPeriodNumber(),fESD->GetBunchCrossNumber(), nESDTracks05, nESDTracks10, nESDTracks14);

    if (eventTriggered)
      fMultiplicity->FillTriggeredEvent(nESDTracks05, nESDTracks10, nESDTracks14);
    
    if (eventTriggered && eventVertex)
      fMultiplicity->FillMeasured(vtx[2], nESDTracks05, nESDTracks10, nESDTracks14);
  }
  else if (fReadMC)   // Processing of MC information
  {
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }

    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    AliStack* stack = mcEvent->Stack();
    if (!stack)
    {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }
    
    AliHeader* header = mcEvent->Header();
    if (!header)
    {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }

    if (fUseMCVertex)
    {
      Printf("WARNING: Replacing vertex by MC vertex. This is for systematical checks only.");
      // get the MC vertex
      AliGenEventHeader* genHeader = header->GenEventHeader();
      TArrayF vtxMC(3);
      genHeader->PrimaryVertex(vtxMC);

      vtx[2] = vtxMC[2];
    }
    
    // get process information
    AliPWG0Helper::MCProcessType processType = AliPWG0Helper::GetEventProcessType(fESD, header, stack, fDiffTreatment);

    Bool_t processEvent = kTRUE;
    if (fSelectProcessType > 0)
    {
      processEvent = kFALSE;

      // non diffractive
      if (fSelectProcessType == 1 && processType == AliPWG0Helper::kND)
        processEvent = kTRUE;

      // single diffractive
      if (fSelectProcessType == 2 && processType == AliPWG0Helper::kSD)
        processEvent = kTRUE;

      // double diffractive
      if (fSelectProcessType == 3 && processType == AliPWG0Helper::kDD)
        processEvent = kTRUE;

      if (!processEvent)
        Printf("Skipping this event, because it is not of the requested process type (%d)", (Int_t) processType);
    }

    if (processEvent)
    {
      // get the MC vertex
      AliGenEventHeader* genHeader = header->GenEventHeader();
      TArrayF vtxMC(3);
      genHeader->PrimaryVertex(vtxMC);

      // get number of tracks from MC
      // loop over mc particles
      Int_t nPrim  = stack->GetNprimary();
      Int_t nMCPart = stack->GetNtrack();

      // tracks in different eta ranges
      Int_t nMCTracks05 = 0;
      Int_t nMCTracks10 = 0;
      Int_t nMCTracks14 = 0;
      Int_t nMCTracksAll = 0;

      // tracks per particle species, in |eta| < 2 (systematic study)
      Int_t nMCTracksSpecies[4]; // (pi, K, p, other)
      for (Int_t i = 0; i<4; ++i)
        nMCTracksSpecies[i] = 0;

      for (Int_t iMc = 0; iMc < nPrim; ++iMc)
      {
        AliDebug(AliLog::kDebug+1, Form("MC Loop: Processing particle %d.", iMc));

        TParticle* particle = stack->Particle(iMc);

        if (!particle)
        {
          AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (mc loop).", iMc));
          continue;
        }

        Bool_t debug = kFALSE;
        if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim, debug) == kFALSE)
        {
          //printf("%d) DROPPED ", iMC);
          //particle->Print();
          continue;
        }

        //printf("%d) OK      ", iMC);
        //particle->Print();

        //if (particle->Pt() < kPtCut)
        //  continue;

        Int_t particleWeight = 1;

        Float_t pt = particle->Pt();

        // in case of systematic study, weight according to the change of the pt spectrum
        // it cannot be just multiplied because we cannot count "half of a particle"
        // instead a random generator decides if the particle is counted twice (if value > 1)
        // or not (if value < 0)
        if (fPtSpectrum)
        {
          Int_t bin = fPtSpectrum->FindBin(pt);
          if (bin > 0 && bin <= fPtSpectrum->GetNbinsX())
          {
            Float_t factor = fPtSpectrum->GetBinContent(bin);
            if (factor > 0)
            {
              Float_t random = gRandom->Uniform();
              if (factor > 1 && random < factor - 1)
              {
                particleWeight = 2;
              }
              else if (factor < 1 && random < 1 - factor)
                particleWeight = 0;
            }
          }
        }

        //Printf("MC weight is: %d", particleWeight);

        if (TMath::Abs(particle->Eta()) < 0.5)
          nMCTracks05 += particleWeight;

        if (TMath::Abs(particle->Eta()) < 1.0)
          nMCTracks10 += particleWeight;

        if (TMath::Abs(particle->Eta()) < 1.3)
          nMCTracks14 += particleWeight;

        nMCTracksAll += particleWeight;
        
        if (particle->Pt() > 0 && TMath::Abs(particle->Eta()) < 1.0)
          fdNdpT->Fill(particle->Pt(), 1.0 / particle->Pt());

        if (fParticleCorrection[0] || fParticleSpecies)
        {
          Int_t id = -1;
          switch (TMath::Abs(particle->GetPdgCode()))
          {
            case 211: id = 0; break;
            case 321: id = 1; break;
            case 2212: id = 2; break;
            default: id = 3; break;
          }

          if (TMath::Abs(particle->Eta()) < etaRange)
            nMCTracksSpecies[id]++;
            
          if (fParticleCorrection[id])
            fParticleCorrection[id]->GetTrackCorrection()->FillGene(vtxMC[2], particle->Eta(), particle->Pt());
        }
      } // end of mc particle

      fMultiplicity->FillGenerated(vtxMC[2], eventTriggered, eventVertex, processType, (Int_t) nMCTracks05, (Int_t) nMCTracks10, (Int_t) nMCTracks14, (Int_t) nMCTracksAll);

      // ESD processing
      Int_t nESDTracks05 = 0;
      Int_t nESDTracks10 = 0;
      Int_t nESDTracks14 = 0;

      // tracks per particle species, in |eta| < 2 (systematic study)
      Int_t nESDTracksSpecies[7]; // (pi, K, p, other, nolabel, doublecount_prim, doublecount_all)
      for (Int_t i = 0; i<7; ++i)
        nESDTracksSpecies[i] = 0;

      Bool_t* foundPrimaries = new Bool_t[nPrim];   // to prevent double counting
      for (Int_t i=0; i<nPrim; i++)
        foundPrimaries[i] = kFALSE;

      Bool_t* foundPrimaries2 = new Bool_t[nPrim];   // to prevent double counting
      for (Int_t i=0; i<nPrim; i++)
        foundPrimaries2[i] = kFALSE;

      Bool_t* foundTracks = new Bool_t[nMCPart];    // to prevent double counting
      for (Int_t i=0; i<nMCPart; i++)
        foundTracks[i] = kFALSE;

      for (Int_t i=0; i<inputCount; ++i)
      {
        Float_t eta = etaArr[i];
        Int_t label = labelArr[i];

        Int_t particleWeight = 1;

        // in case of systematic study, weight according to the change of the pt spectrum
        if (fPtSpectrum)
        {
          TParticle* mother = 0;

          // preserve label for later
          Int_t labelCopy = label;
          if (labelCopy >= 0)
            labelCopy = AliPWG0Helper::FindPrimaryMotherLabel(stack, labelCopy);
          if (labelCopy >= 0)
            mother = stack->Particle(labelCopy);

          // in case of pt study we do not count particles w/o label, because they cannot be scaled
          if (!mother)
            continue;

          // it cannot be just multiplied because we cannot count "half of a particle"
          // instead a random generator decides if the particle is counted twice (if value > 1) 
          // or not (if value < 0)
          Int_t bin = fPtSpectrum->FindBin(mother->Pt());
          if (bin > 0 && bin <= fPtSpectrum->GetNbinsX())
          {
            Float_t factor = fPtSpectrum->GetBinContent(bin);
            if (factor > 0)
            {
              Float_t random = gRandom->Uniform();
              if (factor > 1 && random < factor - 1)
              {
                particleWeight = 2;
              }
              else if (factor < 1 && random < 1 - factor)
                particleWeight = 0;
            }
          }
        }

        //Printf("ESD weight is: %d", particleWeight);

        if (TMath::Abs(eta) < 0.5)
          nESDTracks05 += particleWeight;

        if (TMath::Abs(eta) < 1.0)
          nESDTracks10 += particleWeight;

        if (TMath::Abs(eta) < 1.3)
          nESDTracks14 += particleWeight;

        if (fParticleSpecies)
        {
          Int_t motherLabel = -1;
          TParticle* mother = 0;

          // find mother
          if (label >= 0)
            motherLabel = AliPWG0Helper::FindPrimaryMotherLabel(stack, label);
          if (motherLabel >= 0)
            mother = stack->Particle(motherLabel);

          if (!mother)
          {
            // count tracks that did not have a label
            if (TMath::Abs(eta) < etaRange)
              nESDTracksSpecies[4]++;
          }
          else
          {
            // get particle type (pion, proton, kaon, other) of mother
            Int_t idMother = -1;
            switch (TMath::Abs(mother->GetPdgCode()))
            {
              case 211: idMother = 0; break;
              case 321: idMother = 1; break;
              case 2212: idMother = 2; break;
              default: idMother = 3; break;
            }

            // double counting is ok for particle ratio study
            if (TMath::Abs(eta) < etaRange)
              nESDTracksSpecies[idMother]++;

            // double counting is not ok for efficiency study

            // check if we already counted this particle, this way distinguishes double counted particles (bug/artefact in tracking) or double counted primaries due to secondaries (physics)
            if (foundTracks[label])
            {
              if (TMath::Abs(eta) < etaRange)
                nESDTracksSpecies[6]++;
            }
            else
            {
              foundTracks[label] = kTRUE;

              // particle (primary) already counted?
              if (foundPrimaries[motherLabel])
              {
                if (TMath::Abs(eta) < etaRange)
                  nESDTracksSpecies[5]++;
              }
              else
                foundPrimaries[motherLabel] = kTRUE;
            }
          }
        }

        if (fParticleCorrection[0])
        {
          if (label >= 0 && stack->IsPhysicalPrimary(label))
          {
            TParticle* particle = stack->Particle(label);

            // get particle type (pion, proton, kaon, other)
            Int_t id = -1;
            switch (TMath::Abs(particle->GetPdgCode()))
            {
              case 211: id = 0; break;
              case 321: id = 1; break;
              case 2212: id = 2; break;
              default: id = 3; break;
            }

            // todo check if values are not completely off??

            // particle (primary) already counted?
            if (!foundPrimaries2[label])
            {
              foundPrimaries2[label] = kTRUE;
              fParticleCorrection[id]->GetTrackCorrection()->FillMeas(vtxMC[2], particle->Eta(), particle->Pt());
            }
          }
        }
      }
        
      if (fParticleCorrection[0])
      {
        // if the particle decays/stops before this radius we do not see it
        // 8cm larger than SPD layer 2
        // 123cm TPC radius where a track has about 50 clusters (cut limit)          
        const Float_t endRadius = (fAnalysisMode & AliPWG0Helper::kSPD) ? 8. : 123;
                
        // loop over all primaries that have not been found
        for (Int_t i=0; i<nPrim; i++)
        {
          // already found
          if (foundPrimaries2[i])
            continue;
            
          TParticle* particle = 0;
          TClonesArray* trackrefs = 0;
          mcEvent->GetParticleAndTR(i, particle, trackrefs);
          
          // true primary and charged
          if (!AliPWG0Helper::IsPrimaryCharged(particle, nPrim))
            continue;              
          
          //skip particles with larger |eta| than 3, to keep the log clean, is anyway not included in correction map
          if (TMath::Abs(particle->Eta()) > 3)
            continue;
          
          // skipping checking of process type of daughter: Neither kPBrem, kPDeltaRay nor kPCerenkov should appear in the event generation
          
          // get particle type (pion, proton, kaon, other)
          Int_t id = -1;
          switch (TMath::Abs(particle->GetPdgCode()))
          {
            case 211: id = 4; break;
            case 321: id = 5; break;
            case 2212: id = 6; break;
            default: id = 7; break;
          }            
          
          if (!fParticleCorrection[id])
            continue;
            
          // get last track reference
          AliTrackReference* trackref = dynamic_cast<AliTrackReference*> (trackrefs->Last());
          
          if (!trackref)
          {
            Printf("ERROR: Could not get trackref of %d (count %d)", i, trackrefs->GetEntries());
            particle->Print();
            continue;
          }
            
          // particle in tracking volume long enough...
          if (trackref->R() > endRadius)
            continue;  
          
          if (particle->GetLastDaughter() >= 0)
          {
            Int_t uID = stack->Particle(particle->GetLastDaughter())->GetUniqueID();
            //if (uID != kPBrem && uID != kPDeltaRay && uID < kPCerenkov)
            if (uID == kPDecay)
            {
              // decayed
              
              Printf("Particle %d (%s) decayed at %f, daugher uniqueID: %d:", i, particle->GetName(), trackref->R(), uID);
              particle->Print();
              Printf("Daugthers:");
              for (Int_t d = particle->GetFirstDaughter(); d <= particle->GetLastDaughter(); d++)
                stack->Particle(d)->Print();
              Printf(" ");
              
              fParticleCorrection[id]->GetTrackCorrection()->FillGene(vtxMC[2], particle->Eta(), particle->Pt());
              continue;
            }
          }
          
          if (trackref->DetectorId() == -1)
          {
            // stopped
            Printf("Particle %d stopped at %f:", i, trackref->R());
            particle->Print();
            Printf(" ");
            
            fParticleCorrection[id]->GetTrackCorrection()->FillMeas(vtxMC[2], particle->Eta(), particle->Pt());
            continue;
          }
          
          Printf("Particle %d simply not tracked", i);
          particle->Print();
          Printf(" ");
        }
      }
      
      delete[] foundTracks;
      delete[] foundPrimaries;
      delete[] foundPrimaries2;

//         if ((Int_t) nMCTracks14 > 10 && nESDTracks14 <= 3)
//         {
//             TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
//             printf("WARNING: Event %lld %s (vtx-z = %f, recon: %f, contrib: %d, res: %f) has %d generated and %d reconstructed...\n", tree->GetReadEntry(), tree->GetCurrentFile()->GetName(), vtxMC[2], vtx[2], vtxESD->GetNContributors(), vtxESD->GetZRes(), nMCTracks14, nESDTracks14);
//         }

      if (eventTriggered)
      {
        fMultiplicity->FillTriggeredEvent(nESDTracks05, nESDTracks10, nESDTracks14);
        fMultiplicity->FillNoVertexEvent(vtxMC[2], eventVertex, nMCTracks05, nMCTracks10, nMCTracks14, nESDTracks05, nESDTracks10, nESDTracks14);
//         if (!eventVertex)
//         {
//           if (nESDTracks05 == 0)
//             fTemp1->Fill(nMCTracks05, nESDTracks05);
//         }
      }
      
      if (eventTriggered && eventVertex)
      {
        // fill response matrix using vtxMC (best guess)
        fMultiplicity->FillCorrection(vtxMC[2],  nMCTracks05,  nMCTracks10,  nMCTracks14,  nMCTracksAll,  nESDTracks05, nESDTracks10, nESDTracks14);

        fMultiplicity->FillMeasured(vtx[2], nESDTracks05, nESDTracks10, nESDTracks14);

        if (fParticleSpecies)
          fParticleSpecies->Fill(vtxMC[2], nMCTracksSpecies[0], nMCTracksSpecies[1], nMCTracksSpecies[2], nMCTracksSpecies[3], nESDTracksSpecies[0], nESDTracksSpecies[1], nESDTracksSpecies[2], nESDTracksSpecies[3], nESDTracksSpecies[4], nESDTracksSpecies[5], nESDTracksSpecies[6]);
      }
    }
  }

  if (etaArr)
    delete[] etaArr;
  if (labelArr)
    delete[] labelArr;
}

void AliMultiplicityTask::Terminate(Option_t *)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput) {
    Printf("ERROR: fOutput not available");
    return;
  }

  fMultiplicity = dynamic_cast<AliMultiplicityCorrection*> (fOutput->FindObject("Multiplicity"));
  for (Int_t i = 0; i < 8; ++i)
    fParticleCorrection[i] = dynamic_cast<AliCorrection*> (fOutput->FindObject(Form("correction_%d", i)));
  fParticleSpecies = dynamic_cast<TNtuple*> (fOutput->FindObject("fParticleSpecies"));
  
  fdNdpT = dynamic_cast<TH1*> (fOutput->FindObject("fdNdpT"));

  if (!fMultiplicity)
  {
    AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p", (void*) fMultiplicity));
    return;
  }

  TString fileName("multiplicity");
  if (fSelectProcessType == 1)
    fileName += "ND";
  if (fSelectProcessType == 2)
    fileName += "SD";
  if (fSelectProcessType == 3)
    fileName += "DD";
  fileName += ".root";
  TFile* file = TFile::Open(fileName, "RECREATE");

  fMultiplicity->SaveHistograms();
  for (Int_t i = 0; i < 8; ++i)
    if (fParticleCorrection[i])
      fParticleCorrection[i]->SaveHistograms();
  if (fParticleSpecies)
    fParticleSpecies->Write();
  if (fdNdpT)
    fdNdpT->Write();

  fTemp1 = dynamic_cast<TH1*> (fOutput->FindObject("fTemp1"));
  if (fTemp1)
    fTemp1->Write();
    
  fEtaPhi = dynamic_cast<TH2F*> (fOutput->FindObject("fEtaPhi"));
  if (fEtaPhi)
    fEtaPhi->Write();
  
  fVertex = dynamic_cast<TH3F*> (fOutput->FindObject("vertex_check"));
  if (fVertex)
    fVertex->Write();
  
  for (Int_t i=0; i<3; i++)
  {
    fEta[i] = dynamic_cast<TH1*> (fOutput->FindObject(Form("fEta_%d", i)));
    if (fEta[i])
    {
      fEta[i]->Sumw2();
      Float_t events = fMultiplicity->GetMultiplicityESD(i)->Integral(1, fMultiplicity->GetMultiplicityESD(i)->GetNbinsX());
      if (events > 0)
        fEta[i]->Scale(1.0 / events);
      fEta[i]->Scale(1.0 / fEta[i]->GetXaxis()->GetBinWidth(1));
      fEta[i]->Write();
    }
  }
  
  TObjString option(fOption);
  option.Write();

  file->Close();

  Printf("Written result to %s", fileName.Data());
}
