/**************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
#include <array>

#include <THistManager.h>
#include <TLinearBinning.h>
#include <TList.h>
#include <TString.h>
#include "AliAnalysisTaskK0toPi0Pi0.h"
#include "AliAODConversionMother.h"
#include "AliAODConversionPhoton.h"
#include "AliCaloPhotonCuts.h"
#include "AliClusterContainer.h"
#include "AliConversionPhotonCuts.h"
#include "AliConversionMesonCuts.h"
#include "AliConvEventCuts.h"
#include "AliLog.h"
#include "AliVCluster.h"
#include "AliV0ReaderV1.h"
#include "AliGammaConversionAODBGHandler.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskK0toPi0Pi0)
/// \endcond

AliAnalysisTaskK0toPi0Pi0::AliAnalysisTaskK0toPi0Pi0():
  AliAnalysisTaskSE(),
  fLocalInitialized(kFALSE),
  fCurrentRun(-1),
  fNewFile(kFALSE),
  fV0Reader(nullptr),
  fV0ReaderName("V0ReaderV1"),
  fClusterContainer(nullptr),
  fIsMC(false),
  fWeightJetJetMC(1.),
  fEventPlaneAngle(0.),
  fEventCuts(nullptr),
  fConvPhotonCuts(nullptr),
  fCaloPhotonCuts(nullptr),
  fPi0CutsConvConv(nullptr),
  fPi0CutsCaloCalo(nullptr),
  fPi0CutsConvCalo(nullptr),
  fK0Cuts(nullptr),
  fBGHandler(nullptr),
  fHistos(nullptr),
  fOutput(nullptr)
{

}

AliAnalysisTaskK0toPi0Pi0::AliAnalysisTaskK0toPi0Pi0(const char *name):
  AliAnalysisTaskSE(name),
  fLocalInitialized(kFALSE),
  fCurrentRun(-1),
  fNewFile(kFALSE),
  fV0Reader(nullptr),
  fV0ReaderName("V0ReaderV1"),
  fClusterContainer(nullptr),
  fIsMC(false),
  fWeightJetJetMC(1.),
  fEventPlaneAngle(0),
  fEventCuts(nullptr),
  fConvPhotonCuts(nullptr),
  fCaloPhotonCuts(nullptr),
  fPi0CutsConvConv(nullptr),
  fPi0CutsCaloCalo(nullptr),
  fPi0CutsConvCalo(nullptr),
  fK0Cuts(nullptr),
  fBGHandler(nullptr), 
  fHistos(nullptr),
  fOutput(nullptr)
{
  DefineOutput(1, TList::Class());
}



AliAnalysisTaskK0toPi0Pi0::~AliAnalysisTaskK0toPi0Pi0() {
  if(fHistos) delete fHistos;
}

//=================================== CREATE OUTPUT OBJECTS ====================================================================
void AliAnalysisTaskK0toPi0Pi0::UserCreateOutputObjects(){
  fOutput        = new TList();
  fOutput->SetOwner(kTRUE);
  
  
  // Connecting input V0 reader
  fV0Reader=dynamic_cast<AliV0ReaderV1*>(AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()));
  if(!fV0Reader){
    AliFatal("Error: No V0 Reader");
  }// GetV0Reader

 // fEventCuts = fV0Reader->GetEventCuts();
 
  
  
  // Define histograms
  fHistos = new THistManager("K0stoPi0Pi0");
  
  AliDebug(2, "************* Defining Event Counter Histograms ********************");
  
  
  // Event counter histograms
  fHistos->CreateTH1("hEventQualityBefore", "V0 reader Event Quality (0 = good)", 13, -0.5, 12.5);
  fHistos->CreateTH1("hEventQualityAfter", "V0 reader Event Quality (0 = good)", 13, -0.5, 12.5);
  fHistos->CreateTH1("hEventSelectionStatusBefore", "Event selection status (0 = good)", 14, -0.5, 13.5);
  fHistos->CreateTH1("hEventSelectionStatusAfter", "Event selection status (0 = good)", 14, -0.5, 13.5 );
  fHistos->CreateTH1("hVertexZ", "z-component of the primary vertex; z (cm); Number of events", 1000, -40., 40.);
  fHistos->CreateTH1("hCaloPhotonsBefore", "Number of Events", 13, -0.5, 12.5);
  fHistos->CreateTH1("hCaloPhotonsAfter", "Number of Events", 13, -0.5, 12.5);
  fHistos->CreateTH1("hConvPhotonsBefore", "Number of Events", 13, -0.5, 12.5);
  fHistos->CreateTH1("hConvPhotonsAfter", "Number of Events", 13, -0.5, 12.5);
  
  AliDebug(2, "************* Defining Photon QA Histograms ********************");
  
  // Photon QA
  fHistos->CreateTH1("hCaloPhotonPt", "p_{t}-distribution of the conversion photons; p_{t} (GeV); Yield", 300, 0., 30.);
  fHistos->CreateTH1("hConvPhotonPt", "p_{t}-distribution of the conversion photons; p_{t} (GeV); Yield", 300, 0., 30.);
  fHistos->CreateTH2("hConvPhotonEtaR", "#eta vs conversion radius of conversion photons; #eta; R (cm)", 200, -1.5, 1.5, 300, 0., 300);

  AliDebug(2, "************* Defining Pi0 Histograms ********************");
  
  // Pi0 invariant mass, alpha and opening angle distributions
  const std::array<TString, 3> pi0rec = {"ConvConv", "ConvCalo", "CaloCalo"}; // aka PCM, EMCAL, PCM-EMCAL
  for(const auto &reccase : pi0rec){
    // before selection
    // for candidates in a wide mass region
    fHistos->CreateTH1("hNPi0CandidatesPerEventBefore" +  reccase, "Number of pi0 candidates in event before selection", 101, -0.5, 100.5);
    fHistos->CreateTH2("hMassvsPtPi0Before" + reccase + "All", "inv. mass vs. p_{t} for all #pi^{0} (" + reccase + ") candidates; inv. mass (GeV/c^{2}); p_{t} (GeV/c)", 500, 0., 0.5, 300, 0.3, 30.);
    // only for candidates in the pi0 mass region
    fHistos->CreateTH2("hMassvsPtPi0Before" + reccase + "Sel", "inv. mass vs. p_{t} for selected #pi^{0} (" + reccase + ") candidates; inv. mass (GeV/c^{2}); p_{t} (GeV/c)", 500, 0., 0.5, 300, 0.3, 30.);
    fHistos->CreateTH2("hAlphavsPtPi0Before" + reccase, "#alpha vs p_{t} for selected #p^i{0} (" + reccase + ") candidates; #alpha; p_{t}", 200, -1., 1., 300, 0., 30.);
    fHistos->CreateTH2("hOpeningAnglevsPtPi0Before" + reccase, "Opening angle vs. p_{t} for selected #pi^{0} (" + reccase + ") candidates; opening angle; p_{t} (GeV/c)", 100, 0., 1., 300., 0.3, 30.);

    fHistos->CreateTH1("hPi0Selection" + reccase, "Pi0 selection status bit for reconstruction case " + reccase, 2, -0.5, 1.5);
  	  // after selection
    fHistos->CreateTH1("hNPi0CandidatesPerEventAfter" +  reccase, "Number of pi0 candidates in event before selection", 101, -0.5, 100.5);
  	  fHistos->CreateTH2("hMassvsPtPi0After" + reccase + "All", "inv. mass vs. p_{t} for all #pi^{0} (" + reccase + ") candidates; inv. mass (GeV/c^{2}); p_{t} (GeV/c)", 500, 0., 0.5, 300, 0.3, 30.);
    // only for candidates in the pi0 mass region
    fHistos->CreateTH2("hMassvsPtPi0After" + reccase + "Sel", "inv. mass vs. p_{t} for selected #pi^{0} (" + reccase + ") candidates; inv. mass (GeV/c^{2}); p_{t} (GeV/c)", 500, 0., 0.5, 300, 0.3, 30.);
    fHistos->CreateTH2("hAlphavsPtPi0After" + reccase, "#alpha vs p_{t} for selected #p^i{0} (" + reccase + ") candidates; #alpha; p_{t}", 200, -1., 1., 300, 0., 30.);
    fHistos->CreateTH2("hOpeningAnglevsPtPi0After" + reccase, "Opening angle vs. p_{t} for selected #pi^{0} (" + reccase + ") candidates; opening angle; p_{t} (GeV/c)", 100, 0., 1., 300., 0.3, 30.);
  
  }

  AliDebug(2, "************* Defining K0 Histograms ********************");
  
  // K0short invariant mass distribution and opening angle distributions
  const std::array<TString, 6> k0Shortrec = {"AllConv", "AllCalo", "DiffMixed", "SameMixed", "ConvoCalo","CaloConvo" }; // aka 6 cases 
  for(const auto &reccase1 : k0Shortrec){
  	// before selection
    fHistos->CreateTH1("hNK0CandidatesPerEventBefore" +  reccase1, "Number of K0 candidates in event before selection", 101, -0.5, 100.5);
    fHistos->CreateTH2("hMassvsPtK0ShortBefore" + reccase1, "inv. mass vs. p_{t} for #k^{0}s (" + reccase1 + ") candidates; inv. mass (GeV/c^{2}); p_{t} (GeV/c)",  500, 0.3, 0.6, 300, 0.3, 30.); 
    fHistos->CreateTH2("hOpeningAnglevsPtK0ShortBefore"+ reccase1, "Opening angle vs. p_{t} for  k0Short (" + reccase1 + ") candidates; opening angle; p_{t} (GeV/c)",  100, 0., 1., 300., 0.3, 30.); 
    
    fHistos->CreateTH1("hK0Selection" + reccase1, "Pi0 selection status bit for reconstruction case " + reccase1, 2, -0.5, 1.5);
    
    // after selection
    fHistos->CreateTH1("hNK0CandidatesPerEventAfter" +  reccase1, "Number of K0 candidates in event after selection", 101, -0.5, 100.5);
    fHistos->CreateTH2("hMassvsPtK0ShortAfter" + reccase1, "inv. mass vs. p_{t} for #k^{0}s (" + reccase1 + ") candidates; inv. mass (GeV/c^{2}); p_{t} (GeV/c)",  500, 0.3, 0.6, 300, 0.3, 30.); 
    fHistos->CreateTH2("hOpeningAnglevsPtK0ShortAfter"+ reccase1, "Opening angle vs. p_{t} for  k0Short (" + reccase1 + ") candidates; opening angle; p_{t} (GeV/c)",  100, 0., 1., 300., 0.3, 30.); 
    
    /*
    //background histograms
    fHistos->CreateTH2("hMassvsPtK0ShortBKG" + reccase1, "inv. mass vs. p_{t} for #k^{0}s (" + reccase1 + ") candidates; inv. mass (GeV/c^{2}); p_{t} (GeV/c)",  500, 0.3, 0.6, 300, 0.3, 30.); 
    fHistos->CreateTH2("hOpeningAnglevsPtK0ShortBKG"+ reccase1, "Opening angle vs. p_{t} for  k0Short (" + reccase1 + ") candidates; opening angle; p_{t} (GeV/c)",  100, 0., 1., 300., 0.3, 30.); 
    */
  }
  
  for(auto hist : *(fHistos->GetListOfHistograms())) fOutput->Add(hist);
  
  // Adding cut QA
  
  TList *qaV0reader = new TList;
  qaV0reader->SetName("QA_V0reader");
  qaV0reader->SetOwner(kTRUE);
  qaV0reader->Add(fV0Reader->GetEventCutHistograms());
  qaV0reader->Add(fV0Reader->GetCutHistograms());
  fOutput->Add(qaV0reader);
  
  // Adding to QA histos to the output - change name in order to
  // better identify the histos
  fOutput->Add(fEventCuts->GetCutHistograms());
  fOutput->Add(fConvPhotonCuts->GetCutHistograms());
  fOutput->Add(fCaloPhotonCuts->GetCutHistograms());
  TList *histlist = fPi0CutsConvConv->GetCutHistograms();
  TString histname = "Pi0CutsConvConv_" + TString(histlist->GetName());
  histlist->SetName(histname);
  fOutput->Add(histlist);
  histlist = fPi0CutsCaloCalo->GetCutHistograms();
  histname = "Pi0CutsCalo_" + TString(histlist->GetName());
  histlist->SetName(histname);
  fOutput->Add(histlist);
  histlist = fPi0CutsConvCalo->GetCutHistograms();
  histname = "Pi0CutsConvCalo_" + TString(histlist->GetName());
  histlist->SetName(histname);
  fOutput->Add(histlist);
  histlist = fK0Cuts->GetCutHistograms();
  histname = "K0cuts_" + TString(histlist->GetName());
  histlist->SetName(histname);
  fOutput->Add(histlist);
  
  
  PostData(1, fOutput);
  
}

//=================================== EXEC ONCE ====================================================================
void AliAnalysisTaskK0toPi0Pi0::ExecOnce() {
  if(!fClusterContainer) fClusterContainer = new AliClusterContainer(fInputEvent->IsA() == AliESDEvent::Class() ? "CaloClusters" : "caloClusters");
  fClusterContainer->SetArray(fInputEvent);
   
  if(fConvPhotonCuts){
    fConvPhotonCuts->InitPIDResponse();   
  }
   
}

//=================================== USER EXEC ====================================================================
void AliAnalysisTaskK0toPi0Pi0::UserExec(Option_t *){
  if(!fLocalInitialized) {
    ExecOnce();
    fLocalInitialized = kTRUE;
  }
  
  if(fCurrentRun != fInputEvent->GetRunNumber()) {
    RunChanged();
    fCurrentRun = fInputEvent->GetRunNumber();
  }
  
  
  // do event selection
  // Use the same event selection as for the v0 reader
  // Good events defined as events with event quality 0
  Int_t selectionStatus = fEventCuts->IsEventAcceptedByCut(fV0Reader->GetEventCuts(), fInputEvent, fMCEvent, false, false);
  Int_t eventQuality = fV0Reader->GetEventCuts()->GetEventQuality();
  fHistos->FillTH1("hEventQualityBefore", eventQuality);
  fHistos->FillTH1("hEventSelectionStatusBefore", selectionStatus);
  if(selectionStatus || eventQuality) return;
  fHistos->FillTH1("hEventQualityAfter", eventQuality);
  fHistos->FillTH1("hEventSelectionStatusAfter", selectionStatus);
  fHistos->FillTH1("hVertexZ", fInputEvent->GetPrimaryVertex()->GetZ());

  
  // get Photon candidates
  fHistos->FillTH1("hConvPhotonsBefore",fV0Reader->GetNReconstructedGammas());
  std::vector<AliAODConversionPhoton> conversionPhotons = MakeConversionPhotonCandidates(*fV0Reader, *fConvPhotonCuts);
  MakePhotonQAConv(conversionPhotons, *fEventCuts);
  Int_t numPhotons = conversionPhotons.size();
  fHistos->FillTH1("hConvPhotonsAfter", numPhotons); 
  
  fHistos->FillTH1("hCaloPhotonsBefore",fClusterContainer->GetNEntries());
  std::vector<AliAODConversionPhoton> caloPhotons = MakeCaloPhotonCandidates(*fClusterContainer, *fCaloPhotonCuts);
  MakePhotonQACalo(caloPhotons, *fEventCuts);
  Int_t numCaloPhotons = caloPhotons.size();
  fHistos->FillTH1("hCaloPhotonsAfter", numCaloPhotons);
  

  
  // get Pi0 candidates
  std::vector<AliAODConversionMother> samePi0PCM = MakePi0Candidates(&conversionPhotons, nullptr, *fPi0CutsConvConv),
                                      samePi0EMCAL = MakePi0Candidates(&caloPhotons, nullptr, *fPi0CutsCaloCalo),
                                      mixedPi0 = MakePi0Candidates(&conversionPhotons, &caloPhotons, *fPi0CutsConvCalo);
  
  /*                                    
  // create the BG Handlers and get the number of events
  // collision system 0 for pp
  // center min and center max don't really matter here
  // nEvents is done by using ->GetNumberOfBGEvents()
  // track mult is done by using ->UseTrackMultiplicity()
  // mode, bins Z, bins multiplicity 0,8,5
  AliGammaConversionBGHandler *samePCMHandler = new AliGammaConversionBGHandler(0,0,0,fPi0Cuts->GetNumberOfBGEvents(),fPi0Cuts->UseTrackMultiplicity(),0,8,5); 
  AliGammaConversionBGHandler *sameEMCALHandler = new AliGammaConversionBGHandler(0,0,0,fPi0CutsCaloCalo->GetNumberOfBGEvents(),fPi0CutsCaloCalo->UseTrackMultiplicity(),0,8,5);
  AliGammaConversionBGHandler *mixedHandler = new AliGammaConversionBGHandler(0,0,0,fPi0Cuts->GetNumberOfBGEvents(),fPi0Cuts->UseTrackMultiplicity(),0,8,5); 
   
  Int_t nEventsSamePCM   = samePCMHandler->GetNBackgroundEventsInBuffer(binZ, multiplicity);
  Int_t nEventsSameEMCAL = sameEMCALHandler->GetNBackgroundEventsInBuffer(binZ, multiplicity);  
  Int_t nEventsMixed     = mixedHandler->GetNBackgroundEventsInBuffer(binZ, multiplicity); 
  
  
  std::vector<AliAODConversionMother> samePi0PCM_BG; 
  std::vector<AliAODConversionMother> samePi0EMCAL_BG;
  std::vector<AliAODConversionMother> mixedPi0_BG;  
  
  
  // loop combine all pi0s with all of the events currently in the buffer
  // add it to the buffer once it is already 
  // three for loops, one for each buffer event
  
  for(int i=0;i<nEventsSamePCM;i++){
     	samePi0PCM_BG.insert(handler->GetBGGoodMesons(binZ, multiplicity, i)); 	
  }
  
  for(int h=0;h<nEventsSameEMCAL;h++){
     	samePi0EMCAL_BG.insert(handler->GetBGGoodMesons(binZ, multiplicity, h)); 	
  }
  
  for(int h=0;h<nEventsSameEMCAL;h++){
     	mixedPi0_BG.insert(handler->GetBGGoodMesons(binZ, multiplicity, h)); 	
  }
  
  // add the events now as a tlist
  handler->AddMesonEvent((TList*) samePi0PCM_BG, samePi0PCM_BG->At(i)->GetProductionX(), samePi0PCM_BG->At(i)->GetProductionY() , samePi0PCM_BG->At(i)->GetProductionZ(),multiplicity);
  */

  
  // fill the QA histograms before selection
  AliDebug(1, "New event - make Pi0 candidates");
  MakePi0QA(samePi0PCM, "ConvConv", "Before");
  MakePi0QA(samePi0EMCAL, "CaloCalo", "Before");
  MakePi0QA(mixedPi0, "ConvCalo", "Before");
  AliDebug(1, "New event - finished pi0 candidates");
  
  //make the selections
  std::vector<AliAODConversionMother> samePi0PCMSelection   = SelectMeson(samePi0PCM, *fPi0CutsConvConv, kPi0, "ConvConv");
  std::vector<AliAODConversionMother> samePi0EMCALSelection = SelectMeson(samePi0EMCAL, *fPi0CutsCaloCalo, kPi0, "CaloCalo");
  std::vector<AliAODConversionMother> mixedPi0Selection     = SelectMeson(mixedPi0, *fPi0CutsConvCalo, kPi0, "ConvCalo");
  
  // Additionally fill the QA histograms after the selection
  MakePi0QA(samePi0PCMSelection, "ConvConv", "After");
  MakePi0QA(samePi0EMCALSelection, "CaloCalo", "After");
  MakePi0QA(mixedPi0Selection, "ConvCalo", "After");
                                    
    
  
  
                                      
  // get K0Short candidates
  std::vector<AliAODConversionMother> allPCM = MakeK0ShortCandidates(&samePi0PCMSelection, nullptr, *fK0Cuts); 
  std::vector<AliAODConversionMother> allEMC = MakeK0ShortCandidates(&samePi0EMCALSelection, nullptr, *fK0Cuts);
  std::vector<AliAODConversionMother> PCMEMC = MakeK0ShortCandidates(&samePi0PCMSelection, &mixedPi0Selection, *fK0Cuts);
  std::vector<AliAODConversionMother> EMCPCM = MakeK0ShortCandidates(&samePi0EMCALSelection, &mixedPi0Selection, *fK0Cuts);
  std::vector<AliAODConversionMother> mixedSame = MakeK0ShortCandidates(&samePi0EMCALSelection, &samePi0PCMSelection, *fK0Cuts);
  std::vector<AliAODConversionMother> mixedDiff = MakeK0ShortCandidates(&mixedPi0Selection, nullptr, *fK0Cuts);
  
  /*
  // add duplicates for the background as well
  std::vector<AliAODConversionMother> allPCM_BG = MakeK0ShortCandidates(&samePi0PCM_BG, nullptr, *fK0Cuts);
  std::vector<AliAODConversionMother> allEMC_BG = MakeK0ShortCandidates(&samePi0EMCAL_BG, nullptr, *fK0Cuts);
  std::vector<AliAODConversionMother> PCMEMC_BG = MakeK0ShortCandidates(&samePi0PCM_BG, &mixedPi0, *fK0Cuts);
  std::vector<AliAODConversionMother> EMCPCM_BG = MakeK0ShortCandidates(&samePi0EMCAL_BG, &mixedPi0, *fK0Cuts);
  std::vector<AliAODConversionMother> mixedSame_BG = MakeK0ShortCandidates(&samePi0EMCAL_BG, &samePi0PCM, *fK0Cuts);
  std::vector<AliAODConversionMother> mixedDiff_BG = MakeK0ShortCandidates(&mixedPi0_BG, nullptr, *fK0Cuts); 
  */
  
  // fill the QA histograms prior to selection
  MakeK0ShortQA(allPCM, "AllConv", "Before");
  MakeK0ShortQA(allEMC,"AllCalo", "Before");
  MakeK0ShortQA(PCMEMC,"ConvoCalo", "Before");
  MakeK0ShortQA(EMCPCM,"CaloConvo", "Before");
  MakeK0ShortQA(mixedSame, "SameMixed", "Before");
  MakeK0ShortQA(mixedDiff, "DiffMixed", "Before");
  
  // make the selections
  std::vector<AliAODConversionMother> allPCMSelection    = SelectMeson(allPCM, *fK0Cuts, kK0, "AllConv");
  std::vector<AliAODConversionMother> allEMCSelection    = SelectMeson(allEMC, *fK0Cuts, kK0, "AllCalo");
  std::vector<AliAODConversionMother> PCMEMCSelection    = SelectMeson(PCMEMC, *fK0Cuts, kK0, "ConvoCalo");
  std::vector<AliAODConversionMother> EMCPCMSelection    = SelectMeson(EMCPCM,*fK0Cuts, kK0, "CaloConvo");
  std::vector<AliAODConversionMother> mixedSameSelection = SelectMeson(mixedSame, *fK0Cuts, kK0, "SameMixed");
  std::vector<AliAODConversionMother> mixedDiffSelection = SelectMeson(mixedDiff, *fK0Cuts, kK0, "DiffMixed");
  
  // fill the QA histograms after selection
  MakeK0ShortQA(allPCMSelection, "AllConv", "After");
  MakeK0ShortQA(allEMCSelection,"AllCalo", "After");
  MakeK0ShortQA(PCMEMCSelection,"ConvoCalo", "After");
  MakeK0ShortQA(EMCPCMSelection,"CaloConvo", "After");
  MakeK0ShortQA(mixedSameSelection, "SameMixed", "After");
  MakeK0ShortQA(mixedDiffSelection, "DiffMixed", "After");
   
  
  PostData(1, fOutput);
}


//=================================== MAKE CALO PHOTON CANDIDATES ====================================================================
std::vector<AliAODConversionPhoton> AliAnalysisTaskK0toPi0Pi0::MakeCaloPhotonCandidates(const AliClusterContainer &inputcont, AliCaloPhotonCuts &cuts){
  std::vector<AliAODConversionPhoton> candidates;
  cuts.FillHistogramsExtendedQA(fInputEvent, fIsMC);

  // vertex
  Double_t vertex[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  int clusterindex = 0;
  for(auto c : inputcont.all()) {
    clusterindex++;
    if(!cuts.ClusterIsSelected(c, fInputEvent, fMCEvent, fIsMC, 1., c->GetID())) continue;   // weight to be added later

    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    c->GetMomentum(clusterVector,vertex);

    // convert to AODConversionPhoton
    AliAODConversionPhoton photonCandidate(&clusterVector);

    // Flag Photon as CaloPhoton
    photonCandidate.SetIsCaloPhoton();
    photonCandidate.SetCaloClusterRef(clusterindex);

    // get MC label
    if(fIsMC>0){
      photonCandidate.SetNCaloPhotonMCLabels(c->GetNLabels());
      for (UInt_t k = 0; k < c->GetNLabels(); k++){
        if(k < 50) photonCandidate.SetCaloPhotonMCLabel(k,c->GetLabels()[k]);
      }
    }

    candidates.push_back(photonCandidate);
  }

  return candidates;
}

//=================================== MAKE CONVERSION PHOTON CANDIDATES ====================================================================
std::vector<AliAODConversionPhoton> AliAnalysisTaskK0toPi0Pi0::MakeConversionPhotonCandidates(const AliV0ReaderV1 &reader, AliConversionPhotonCuts &cuts) {
  std::vector<AliAODConversionPhoton> candidates;
  Int_t nV0 = 0;
  std::vector<AliAODConversionPhoton *> GammaCandidatesStepOne;
  TList GammaCandidatesStepTwo;
  
  // Loop over Photon Candidates allocated by ReaderV1  
  for(auto photon : reader){
    if(!cuts.PhotonIsSelected(photon ,fInputEvent)) continue;
    if(!cuts.InPlaneOutOfPlaneCut(photon->GetPhotonPhi(), fEventPlaneAngle)) continue;
    if(!cuts.UseElecSharingCut() && !cuts.UseToCloseV0sCut()){
      candidates.push_back(*(static_cast<AliAODConversionPhoton *>(photon))); // if no second loop is required add to events good gammas
    }else if(cuts.UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
      cuts.FillElectonLabelArray(static_cast<AliAODConversionPhoton *>(photon), nV0);
      nV0++;
      GammaCandidatesStepOne.push_back(static_cast<AliAODConversionPhoton *>(photon));
    }else if(!cuts.UseElecSharingCut() && cuts.UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
      GammaCandidatesStepTwo.Add(static_cast<AliAODConversionPhoton *>(photon));
    }
  }

  if(cuts.UseElecSharingCut()){
    int iV0 = 0;
    for(auto photon : GammaCandidatesStepOne){
      if(!cuts.RejectSharedElectronV0s(photon,iV0,GammaCandidatesStepOne.size())) continue;
      if(!cuts.UseToCloseV0sCut()){ // To Close v0s cut diabled, step two not needed
        candidates.push_back(*photon);
      } else GammaCandidatesStepTwo.Add(photon); // Close v0s cut enabled -> add to list two
      iV0++;
    }
  }

  if(cuts.UseToCloseV0sCut()){
    for(int i = 0; i < GammaCandidatesStepTwo.GetEntries(); i++){
      AliAODConversionPhoton *photon = static_cast<AliAODConversionPhoton *>(GammaCandidatesStepTwo.At(i));
      if(!cuts.RejectToCloseV0s(photon, &GammaCandidatesStepTwo,i)) continue;
      candidates.push_back(*photon); // Add gamma to current cut TList
    }
  }
  return candidates;
}

//=================================== SELECT MESON ====================================================================
std::vector<AliAODConversionMother> AliAnalysisTaskK0toPi0Pi0::SelectMeson(std::vector<AliAODConversionMother> &candidates, 
																		   AliConversionMesonCuts &cuts, MesonType_t meson, const char *reccase){
  std::vector<AliAODConversionMother> selectedCandidates; 
  TString mesonName;
  switch(meson) {
  case kPi0: mesonName = "Pi0"; break;
  case kK0: mesonName = "K0"; break;
  };
  TString qaname = "h" + mesonName + "Selection" + reccase;
  for(auto candidate: candidates){
    Int_t selectionStatus = cuts.MesonIsSelected(&candidate, kTRUE, 0);
    fHistos->FillTH1(qaname, selectionStatus);
  	  if(!selectionStatus) continue;
    selectedCandidates.push_back(candidate);
  }
  
  return selectedCandidates; 


}

//=================================== MAKE PI0 CANDIDATES ====================================================================
std::vector<AliAODConversionMother> AliAnalysisTaskK0toPi0Pi0::MakePi0Candidates(const std::vector<AliAODConversionPhoton> *primaryLeg,
                                                                                 const std::vector<AliAODConversionPhoton> *secondaryLeg,
                                                                                 AliConversionMesonCuts &cuts){
  // secondary leg: optional argument, if different methods for photon identification (i.e. PCM-EMCAL) is used
  AliDebug(1, Form("Make Pi0 candidates: Number of photons in first leg: %lu, second leg: %lu", primaryLeg->size(), secondaryLeg ? secondaryLeg->size() : 0));
  std::vector<AliAODConversionMother> candidates;
  if(secondaryLeg){
    // Different methods for photon identification identification
    for(auto primphoton : *primaryLeg){
      for(auto secphoton : *secondaryLeg) {
        AliAODConversionMother candidate(&primphoton, &secphoton);
        AliDebug(2, Form("Made Pi0 candidate same - Mass %f [%f - %f] GeV/c2", candidate.M(), cuts.GetSelectionLow(), cuts.GetSelectionHigh()));
        if(!cuts.CheckWhetherInMassRange(candidate.M()))continue;
        AliDebug(2, "Candidate in mass window");
        // Do Pi0 selection
        //if(!cuts.MesonIsSelected(&candidate, kTRUE, 0)) continue;      // Rapidity shift needed when going to asymmetric systems
        candidates.push_back(candidate);
      }
    }
  } else {
    // Same method for photon identification
    for(auto primiter = primaryLeg->begin(); primiter != primaryLeg->end(); ++primiter){
      for(auto seciter = primiter + 1; seciter != primaryLeg->end(); ++seciter){
        AliAODConversionMother candidate(&(*primiter), &(*seciter));
        AliDebug(2, Form("Made Pi0 candidate mixed - Mass %f [%f - %f] GeV/c2", candidate.M(), cuts.GetSelectionLow(), cuts.GetSelectionHigh()));
        if(!cuts.CheckWhetherInMassRange(candidate.M()))continue;
        AliDebug(2, "Candidate in mass window");
        // Do pi0 selection
        //if(!cuts.MesonIsSelected(&candidate, kTRUE, 0)) continue;      // Rapidity shift needed when going to asymmetric systems
        candidates.push_back(candidate);
      }
    }
  }
  AliDebug(1, Form("Found %lu pi0 candidates in event\n", candidates.size()));
  return candidates;
}

// =================================== MAKE K0 SHORT CANDIDATES ====================================================================
std::vector<AliAODConversionMother> AliAnalysisTaskK0toPi0Pi0::MakeK0ShortCandidates(const std::vector<AliAODConversionMother> *primaryLeg,
                                                                                     const std::vector<AliAODConversionMother> *secondaryLeg,
                                                                                    AliConversionMesonCuts &cuts){
  std::vector<AliAODConversionMother> candidates;
  if(secondaryLeg) {
    // Different methods for Pi0 identification (one same, one mixed)
    for(const auto &primpi0 : *primaryLeg) {
      for(const auto &secpi0 : *secondaryLeg) {
        AliAODConversionMother candidate(&primpi0, &secpi0);
        //if(!cuts.MesonIsSelected(&candidate, kTRUE, 0)) continue;
        candidates.push_back(candidate);
      }
    }
  } else {
    // Same methods for Pi0 identification (both same or both mixed)
    for(auto primpi0 = primaryLeg->begin(); primpi0 != primaryLeg->end(); ++primpi0) {
      for(auto secpi0 = primpi0 + 1; secpi0 != primaryLeg->end(); ++secpi0) {
        AliAODConversionMother candidate(&(*primpi0), &(*secpi0));
        //if(!cuts.MesonIsSelected(&candidate, kTRUE, 0)) continue;
        candidates.push_back(candidate);
      }
    }
  }
  return candidates;
}

//=================================== MAKE PHOTON QA CALO ====================================================================
void AliAnalysisTaskK0toPi0Pi0::MakePhotonQACalo(const std::vector<AliAODConversionPhoton> &photons, AliConvEventCuts &cuts) {
  for(const auto &photon : photons) {
    fHistos->FillTH1("hCaloPhotonPt", photon.Pt());
  }
}

//=================================== MAKE PHOTON QA CONV ====================================================================
void AliAnalysisTaskK0toPi0Pi0::MakePhotonQAConv(const std::vector<AliAODConversionPhoton> &photons, AliConvEventCuts &cuts) {
  for(const auto &photon : photons) {
    fHistos->FillTH1("hConvPhotonPt", photon.Pt());
    fHistos->FillTH2("hConvPhotonEtaR", photon.Eta(), photon.GetConversionRadius());
  }
}

//=================================== MAKE PI0 QA ====================================================================
void AliAnalysisTaskK0toPi0Pi0::MakePi0QA(const std::vector<AliAODConversionMother> &pi0s, const char *reccase, TString selectionStatus){
  TString reccaseString = reccase;
  fHistos->FillTH1("hNPi0CandidatesPerEvent" + selectionStatus + reccaseString, pi0s.size());
  for(const auto &pi0 : pi0s) {
  	  fHistos->FillTH2("hMassvsPtPi0" + selectionStatus + reccaseString + "All", pi0.M(), pi0.Pt());
  	  // if in the pi0 mass region  
  	  if((0.08 <=pi0.M()) && (pi0.M()<= 0.145)){
  	    fHistos->FillTH2("hMassvsPtPi0" + selectionStatus + reccaseString + "Sel", pi0.M(),pi0.Pt());
  	    fHistos->FillTH2("hAlphavsPtPi0" + selectionStatus + reccaseString, pi0.GetAlpha(), pi0.Pt());
  	    fHistos->FillTH2("hOpeningAnglevsPtPi0" + selectionStatus + reccaseString, pi0.GetOpeningAngle(), pi0.Pt());
  	  }
  	}
}

//=================================== MAKE KO SHORT QA ====================================================================
void AliAnalysisTaskK0toPi0Pi0::MakeK0ShortQA(const std::vector<AliAODConversionMother> &k0s,const char *reccase, TString selectionStatus){
  TString reccaseString = reccase;
  fHistos->FillTH1("hNK0CandidatesPerEvent" + selectionStatus + reccaseString, k0s.size());
  for(const auto &k0: k0s) {
  	  fHistos->FillTH2("hMassvsPtK0Short" + selectionStatus + reccaseString, k0.M(), k0.Pt());
  	  fHistos->FillTH2("hOpeningAnglevsPtK0Short" + selectionStatus + reccaseString, k0.GetOpeningAngle(), k0.Pt());
  }
}

//=================================== ADD CLUSTER CONTAINER ====================================================================
AliClusterContainer *AliAnalysisTaskK0toPi0Pi0::AddClusterContainer(const char *name) {
  fClusterContainer = new AliClusterContainer(name);
  return fClusterContainer;
}



