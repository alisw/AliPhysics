/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
 ***************************************************************************/
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>

#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliAnalysisTaskEmcalJetEnergyFlow.h"
//#include "AliAODJet.h"
#include "TList.h"
#include "TArrayI.h"
#include "TArrayS.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetEnergyFlow);
/// \endcond

/**
 * Default constructor. Needed by I/O
 */
AliAnalysisTaskEmcalJetEnergyFlow::AliAnalysisTaskEmcalJetEnergyFlow():
	AliAnalysisTaskEmcalJet(),
	fHistManager()// ,fOutput{0}
{
}
/**
 * Standard constructor, intended for use by the the user.
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalJetEnergyFlow::AliAnalysisTaskEmcalJetEnergyFlow(const char* name):
AliAnalysisTaskEmcalJet(name, kTRUE),
fHistManager(name) // ,fOutput{0}
{
	SetMakeGeneralHistograms(kFALSE);
//	DefineInput(0, TChain::Class());
//	DefineOutput(1,TList::Class());

}
/**
 * Destructor
 */

AliAnalysisTaskEmcalJetEnergyFlow::~AliAnalysisTaskEmcalJetEnergyFlow()
{
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalJetEnergyFlow::UserCreateOutputObjects()
{
	AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

	AllocateJetHistograms();
	AllocateTrackHistograms();
	AllocateClusterHistograms();
	AllocateCellHistograms();
	AllocateEnergyflowHistograms();

  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalJetEnergyFlow::ExecOnce()
{
	AliAnalysisTaskEmcalJet::ExecOnce();
}
/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 *       */

Bool_t AliAnalysisTaskEmcalJetEnergyFlow::Run()
{
	return kTRUE;
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */

Bool_t AliAnalysisTaskEmcalJetEnergyFlow::FillHistograms()
{
	DoJetLoop();
	DoTrackLoop();
	DoClusterLoop();
	DoCellLoop();
	FillEFHistograms();	
        return kTRUE;
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmcalJetEnergyFlow::Terminate(Option_t *)
{
}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */

AliAnalysisTaskEmcalJetEnergyFlow* AliAnalysisTaskEmcalJetEnergyFlow::AddTaskEmcalJetEnergyFlow(
              
	 const char *ntracks,
         const char *nclusters,
         const char *ncells,
         const char *suffix              )
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetEnergyFlow", "No analysis manager to connect to.");
    return 0;
  }
 
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalJetEnergyFlow", "This task requires an input event handler");
    return 0;
  }

  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD
  };

  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  //-------------------------------------------------------
  //   Init the task and do settings
  //-------------------------------------------------------

 TString trackName(ntracks);
  TString clusName(nclusters);
  TString cellName(ncells);

  if (trackName == "usedefault") {
    if (dataType == kESD) {
      trackName = "Tracks";
    }
    else if (dataType == kAOD) {
      trackName = "tracks";
    }
    else {
      trackName = "";
    }
  }

  if (clusName == "usedefault") {
    if (dataType == kESD) {
      clusName = "CaloClusters";
    }
    else if (dataType == kAOD) {
      clusName = "caloClusters";
    }
    else {
      clusName = "";
    }
  }

  if (cellName == "usedefault") {
    if (dataType == kESD) {
      cellName = "EMCALCells";
    }
    else if (dataType == kAOD) {
      cellName = "emcalCells";
    }
    else {
      cellName = "";
    }
  }

  TString name("AliAnalysisTaskEmcalJetEnergyFlow");
  if (!trackName.IsNull()) {
    name += "_";
    name += trackName;
  }
  if (!clusName.IsNull()) {
    name += "_";
    name += clusName;
  }
  if (!cellName.IsNull()) {
    name += "_";
    name += cellName;
  }
  if (strcmp(suffix,"") != 0) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskEmcalJetEnergyFlow* EFTask = new AliAnalysisTaskEmcalJetEnergyFlow(name);
  EFTask->SetCaloCellsName(cellName);
  EFTask->SetVzRange(-10,10);

  if (trackName == "mcparticles") {
    EFTask->AddMCParticleContainer(trackName);
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    EFTask->AddTrackContainer(trackName);
  }
  else if (!trackName.IsNull()) {
    EFTask->AddParticleContainer(trackName);
  }
  EFTask->AddClusterContainer(clusName);

 //-------------------------------------------------------
 // Final settings, pass to manager and set the containers
 //-------------------------------------------------------
   mgr->AddTask(EFTask);
 //Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
      TList::Class(),AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput  (EFTask, 0,  cinput1 );
  mgr->ConnectOutput (EFTask, 1, coutput1 );

  return EFTask; 
}

/*
 * This function allocates the histograms for basic jet QA.
 * A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
 * per each jet container and per each centrality bin.
 */

void AliAnalysisTaskEmcalJetEnergyFlow::AllocateJetHistograms(){
  TString histname;
  TString histtitle;
  TString groupname;
  Int_t fNPtBins = 100;
  Double_t fMinPtBin = 0.0;
  Double_t fMaxPtBin = 100.0;
  Int_t fNEtaBins = 100;
  Double_t fMaxEtaBin = 1.0;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    // Protection against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/histJetPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNPtBins, fMinPtBin, fMaxPtBin);

      histname = TString::Format("%s/histJetArea_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{A}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNPtBins / 2, 0, 3);

      histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNPtBins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histJetEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNEtaBins, -fMaxEtaBin, fMaxEtaBin);

      histname = TString::Format("%s/histNJets_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of jets;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 500);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 100, 0, 100);
      }

      if (!jetCont->GetRhoName().IsNull()) {
        histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histname.Data());
        fHistManager.CreateTH1(histname, histtitle, fNPtBins, -fMaxPtBin / 2, fMaxPtBin / 2);
      }
    }
  }

}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetEnergyFlow::DoJetLoop(){

TString histname;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    UInt_t count = 0;
    for(auto jet : jetCont->accepted()) {
      if (!jet) continue;
      count++;

      histname = TString::Format("%s/histJetPt_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Pt());

      histname = TString::Format("%s/histJetArea_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Area());

      histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Phi());

      histname = TString::Format("%s/histJetEta_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Eta());

      if (jetCont->GetRhoParameter()) {
        histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, jet->Pt() - jetCont->GetRhoVal() * jet->Area());
      }
    }
    histname = TString::Format("%s/histNJets_%d", groupname.Data(), fCentBin);
    fHistManager.FillTH1(histname, count);
  }
}

void AliAnalysisTaskEmcalJetEnergyFlow::AllocateEnergyflowHistograms(){

  TString histname;
  TString histtitle;
  TString groupname;
  Double_t Rjet = 0.1;
  Int_t fNPtBins = 100;
  Double_t fMinPtBin = 0.0;
  Double_t fMaxPtBin = 100.0;
  Int_t fNEtaBins = 100;
  Double_t fMaxEtaBin = 1.0;
  Int_t fNDPtBins = 100;
  Double_t fMaxDPtBin = 80.0;
  Double_t fMinDPtBin = -20.0;
  Int_t fNDRBins = 100;
  Double_t fMaxDRBin = 0.4;
  Double_t fMinDRBin = 0.0;
  // This step is done in order to avoid hardcoded jet radii step value
   AliJetContainer* jetCont1= dynamic_cast<AliJetContainer*>(fJetCollArray[0]); // One container for the jets of the lower R (of each comparison pair)
   AliJetContainer* jetCont2= dynamic_cast<AliJetContainer*>(fJetCollArray[1]); // One container for the jets of the higher R (of each comparison pair)
//   Double_t Rjet = jetCont2->GetJetRadius() - jetCont1->GetJetRadius();


  Int_t Pair_number = fJetCollArray.GetEntries()-1;
  for(Int_t i=0;i<Pair_number;i++){
	histname = TString::Format("hJetPtDeltaPt_R%d",int(Rjet*(i+1)*100));
	histtitle = TString::Format("#DeltaP_{t} between %.2f and %.2f jet radii;P_{t,jet(R=%.2f)}(GeV/c);#Delta P_{t}(GeV/c)",Rjet*(i+1),Rjet*(i+2),Rjet*(i+1));
	fHistManager.CreateTH2(histname,histtitle,fNPtBins,fMinPtBin,fMaxPtBin,fNDPtBins,fMinDPtBin,fMaxDPtBin);

        histname = TString::Format("hJetPtDeltaRDeltaPt_R%d",int(Rjet*(i+1)*100));
        histtitle = TString::Format("#DeltaP_{t} between %.2f and %.2f jet radii vs #DeltaR;P_{t,jet(R=%.2f)};#DeltaR;#DeltaP_{t}",Rjet*(i+1),Rjet*(i+2),Rjet*(i+1));
        fHistManager.CreateTH3(histname,histtitle,fNPtBins,fMinPtBin,fMaxPtBin,fNDRBins,fMinDRBin,fMaxDRBin,fNDPtBins,fMinDPtBin,fMaxDPtBin);
	
	histname = TString::Format("hJetPtDeltaR_R%d",int(Rjet*(i+1)*100));
        histtitle = TString::Format("#DeltaR between %.2f and %.2f jet radii;P_{t,R=%.2f};#DeltaR",Rjet*(i+1),Rjet*(i+2),Rjet*(i+1));
        fHistManager.CreateTH2(histname,histtitle,fNPtBins,fMinPtBin,fMaxPtBin,fNDRBins,fMinDRBin,fMaxDRBin);

        histname = TString::Format("hDeltaPt_overPt_R%d",int(Rjet*(i+1)*100));
        histtitle = TString::Format("#DeltaP_{t} between %.2f and %.2f jet radii over P_{t,R=%.2f};P_{t,jet(R=%f)}(GeV/c);#Delta P_{t}/P_{t,jet(R=%f)}",Rjet*(i+1),Rjet*(i+2),Rjet*(i+1),Rjet*(i+1),Rjet*(i+1));
        fHistManager.CreateTH2(histname,histtitle,fNPtBins,fMinPtBin,fMaxPtBin,100,0,2);

	histname = TString::Format("hEtaJetDeltaR_R%d",int(Rjet*(i+1)*100));
        histtitle = TString::Format("#DeltaR between %.2f and %.2f jet radii;#eta_{jet(R=%f)};#Delta R",Rjet*(i+1),Rjet*(i+2),Rjet*(i+1));
        fHistManager.CreateTH2(histname,histtitle,fNEtaBins,-fMaxEtaBin,fMaxEtaBin,fNDRBins,fMinDRBin,fMaxDRBin);

        histname = TString::Format("hDeltaPtvPtvDeltaEta_R%d",int(Rjet*(i+1)*100));
        histtitle = TString::Format("#Delta#eta vs the #DeltaP_{t} between %.2f and %.2f jet radii;#Delta P_{t};P_{t,jet(R=%f)}(GEV/v);#Delta #eta",Rjet*(i+1),Rjet*(i+2),Rjet*(i+1));
        fHistManager.CreateTH3(histname,histtitle,fNDPtBins,fMinDPtBin,fMaxDPtBin,fNPtBins,fMinPtBin,fMaxPtBin,200,0,2);

	histname = TString::Format("hDeltaPtvPtvEta_low_R%d",int(Rjet*(i+1)*100));
        histtitle = TString::Format("#DeltaP_{t} between %.2f and %.2f jet radii vs #Eta vs P_{t} of %.2f;#Delta P_{t};P_{t,jet(R=%f)}(GEV/v); #eta_{R=%f}",Rjet*(i+1),Rjet*(i+2),Rjet*(i+1),Rjet*(i+1),Rjet*(i+1));
        fHistManager.CreateTH3(histname,histtitle,fNDPtBins,fMinDPtBin,fMaxDPtBin,fNPtBins,fMinPtBin,fMaxPtBin,fNEtaBins,-fMaxEtaBin,fMaxEtaBin);

	histname = TString::Format("hDeltaPtvPtvEta_high_R%d",int(Rjet*(i+1)*100));
        histtitle = TString::Format("#DeltaP_{t} between %.2f and %.2f jet radii vs #Eta vs P_{t} of %.2f;#Delta P_{t};P_{t,jet(R=%f)}(GEV/v); #eta_{R=%f}",Rjet*(i+1),Rjet*(i+2),Rjet*(i+1),Rjet*(i+1),Rjet*(i+2));
        fHistManager.CreateTH3(histname,histtitle,fNDPtBins,fMinDPtBin,fMaxDPtBin,fNPtBins,fMinPtBin,fMaxPtBin,fNEtaBins,-fMaxEtaBin,fMaxEtaBin);

	histname = TString::Format("hDeltaPtvPtvMultiplicity_R%d",int(Rjet*(i+1)*100));
	histtitle = TString::Format("#DeltaP_{t} between %.2f and %.2f jet radii vs P_{t} vs Multiplicity;#Delta P_{t};P_{t,jet(R=%f)}(GEV/v); Multiplicity",Rjet*(i+1),Rjet*(i+2),Rjet*(i+1));
	fHistManager.CreateTH3(histname,histtitle,fNDPtBins,fMinDPtBin,fMaxDPtBin,fNPtBins,fMinPtBin,fMaxPtBin,50,0,50);

	histname = TString::Format("hMatchedJetPt_R%d",int(Rjet*(i+1)*100));
	histtitle = TString::Format("Matched jet P_{t} spectrum of R=%.2f",Rjet*(i+1));
	fHistManager.CreateTH1(histname, histtitle,fNPtBins,fMinPtBin,fMaxPtBin);

        histname = TString::Format("hMatchedJetEta_R%d",int(Rjet*(i+1)*100));
        histtitle = TString::Format("Matched jet #eta spectrum of R=%.2f",Rjet*(i+1));
        fHistManager.CreateTH1(histname, histtitle,fNEtaBins,-fMaxEtaBin,fMaxEtaBin);

	if(i==Pair_number-1){
	histname = TString::Format("hMatchedJetPt_R%d",int(Rjet*(i+2)*100));
        histtitle = TString::Format("Matched jet P_{t} spectrum of R=%.2f",Rjet*(i+2));
        fHistManager.CreateTH1(histname, histtitle,fNPtBins,fMinPtBin,fMaxPtBin);

        histname = TString::Format("hMatchedJetEta_R%d",int(Rjet*(i+2)*100));
        histtitle = TString::Format("Matched jet #eta spectrum of R=%.2f",Rjet*(i+2));
        fHistManager.CreateTH1(histname, histtitle,fNEtaBins,-fMaxEtaBin,fMaxEtaBin);
			}

	}
}

void AliAnalysisTaskEmcalJetEnergyFlow::FillEFHistograms(){


Int_t Njets = 200;	//Just a high number so that the matching matrix created by the helper task will always have the size of the input jet lists
Int_t &kLowRJets = Njets;	
Int_t &kHighRJets = Njets;

TList LowRJetsList;	// List of the accepted generated(lower R) jets
TList HighRJetsList;    // List of the accepted rec(higher R) jets

//This array points to the low R jet that matches to each high R jet
TArrayI iLowRIndex;

// This array points to the high R jet that matches to each low R jet
TArrayI iHighRIndex;	

TString histname1;	//This refers to the the Delta Pt histograms
TString histname2;	//This refers to the DeltaR histograms
TString histname3; 	//This refers to the Dpt/pt histograms
TString histname4;	//This refers to the Delta Pt vs Delta R vs Pt histograms 
TString histname5;	//This refers to the Delta R vs Eta histograms
TString histname6;	//This refers to the Dpt vs Pt vs DEta histograms
TString histname7;	//This refers to the Matched Jet Pt histograms
TString histname8;	//This refers to the Matched Jet Eta histograms
TString histname9;	//This refers to the Matched Jet Pt histograms at the last R iteration
TString histname10;	//This refers to the Matched Jet Eta histograms at the last R iteration
TString histname11;	//This refers to the DeltaPt vs Pt vs Eta low Rjet histograms
TString histname12;     //This refers to the DeltaPt vs Pt vs Eta high Rjet histograms
TString histname13;     //This refers to the DeltaPt vs Pt vs Multiplicity histograms

TString groupname;
Double_t Rjet=0.1;		//This is the size of the radius step with which we expand the jet
Double_t DeltaPt = 0.0;
Double_t DeltaEta = 0.0;
Float_t Max_dist = 0.2;
Float_t max_eta = 0.5;
AliJetContainer* jetCont1=0; // One container for the jets of the lower R (of each comparison pair)
AliJetContainer* jetCont2=0; // One container for the jets of the higher R (of each comparison pair)

//Loop over the number of comparison pairs
for (Int_t i=0;i<fJetCollArray.GetEntries()-1;i++)
{	LowRJetsList.Clear();
	HighRJetsList.Clear();	
	// Casting the lower R half of the pair to a container and getting the accepted jets to a list
	jetCont1 = dynamic_cast<AliJetContainer*>(fJetCollArray[i]);
	for(auto jet:jetCont1->accepted()) LowRJetsList.Add(jet);
	
	// Casting the higher R half of the pair to a container and getting the accepted jets to a list
	jetCont2 = dynamic_cast<AliJetContainer*>(fJetCollArray[i+1]);
	for(auto jet:jetCont2->accepted()) HighRJetsList.Add(jet);
//	Rjet = jetCont2->GetJetRadius()- jetCont1->GetJetRadius();

	// Histogram name for the Dpt histograms
	histname1 = TString::Format("hJetPtDeltaPt_R%d",int(Rjet*(i+1)*100));

	//Histogram name for the DR histograms
	histname2 = TString::Format("hJetPtDeltaR_R%d",int(Rjet*(i+1)*100));

	//Histogram name for the Dpt/Pt histograms
	histname3 = TString::Format("hDeltaPt_overPt_R%d",int(Rjet*(i+1)*100));	

        //Histogram name for the Dpt vs Eta histograms
        histname4 = TString::Format("hJetPtDeltaRDeltaPt_R%d",int(Rjet*(i+1)*100));

        //Histogram name for the DR vs Eta histograms
        histname5 = TString::Format("hEtaJetDeltaR_R%d",int(Rjet*(i+1)*100));
	       		
	//Histogram name for the Dpt vs Pt vs DEta histograms
	histname6 = TString::Format("hDeltaPtvPtvDeltaEta_R%d",int(Rjet*(i+1)*100));

	//Histogram name for the Matched Jet Pt histograms
	histname7 = TString::Format("hMatchedJetPt_R%d",int(Rjet*(i+1)*100));

	//Histogram name for the Matched Jet Eta histograms
	histname8 = TString::Format("hMatchedJetEta_R%d",int(Rjet*(i+1)*100));

	if(i==fJetCollArray.GetEntries()-2){
	//Histogram name for the Matched Jet Pt histograms at the last R iteration
	histname9 = TString::Format("hMatchedJetPt_R%d",int(Rjet*(i+2)*100));
	//Histogram name for the Matched Jet Pt histograms at the last R iteration
        histname10 = TString::Format("hMatchedJetEta_R%d",int(Rjet*(i+2)*100));
						}
	//Histogram name for the DeltaPt vs Pt vs Eta low Rjet histograms
	histname11 = TString::Format("hDeltaPtvPtvEta_low_R%d",int(Rjet*(i+1)*100));

	//Histogram name for the DeltaPt vs Pt vs Eta high Rjet histograms
	histname12 = TString::Format("hDeltaPtvPtvEta_high_R%d",int(Rjet*(i+1)*100));

	//Histogram name for the DeltaPt vs Pt vs Multiplicity histograms
	histname13 = TString::Format("hDeltaPtvPtvMultiplicity_R%d",int(Rjet*(i+1)*100));

	if(LowRJetsList.GetEntries()==0||HighRJetsList.GetEntries()==0) continue;
	iLowRIndex.Set(HighRJetsList.GetEntries());
	iHighRIndex.Set(LowRJetsList.GetEntries());
	JetMatcher(&LowRJetsList,kLowRJets,&HighRJetsList,kHighRJets, iLowRIndex,iHighRIndex,0,Max_dist,max_eta);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for (Int_t j=0; j<iHighRIndex.GetSize()-1;j++)
                {
                        if(iHighRIndex[j]>=0){
                        Int_t match_index = iHighRIndex[j];
			if (iLowRIndex[match_index]==j){

                Double_t pt_low = (dynamic_cast<AliEmcalJet*>(LowRJetsList.At(j)))->Pt();
                Double_t pt_high = (dynamic_cast<AliEmcalJet*>(HighRJetsList.At(match_index)))->Pt();
		Double_t eta_low = (dynamic_cast<AliEmcalJet*>(LowRJetsList.At(j)))->Eta(); 
		Double_t eta_high = (dynamic_cast<AliEmcalJet*>(HighRJetsList.At(match_index)))->Eta();
		DeltaPt = pt_high-pt_low;
		Double_t DeltaR = (dynamic_cast<AliEmcalJet*>(LowRJetsList.At(j)))->DeltaR((dynamic_cast<AliEmcalJet*>(HighRJetsList.At(match_index))));
		if(eta_low*eta_high>=0) DeltaEta =fabs(eta_high - eta_low);
		else{ if(eta_low>0) DeltaEta = eta_low-eta_high;
		      else DeltaEta = eta_high -eta_low;
			}
                        fHistManager.FillTH2(histname1,pt_low,DeltaPt);

                        fHistManager.FillTH2(histname2,pt_low,DeltaR);
		
			fHistManager.FillTH2(histname3,pt_low,DeltaPt/pt_low);
			fHistManager.FillTH3(histname4,pt_low,DeltaR,DeltaPt);
			fHistManager.FillTH2(histname5,eta_low,DeltaR);
			fHistManager.FillTH3(histname6,DeltaPt,pt_low,DeltaEta);
			fHistManager.FillTH1(histname7,pt_low);
			fHistManager.FillTH1(histname8,eta_low);
			if(i==fJetCollArray.GetEntries()-2){
			fHistManager.FillTH1(histname9,pt_high);
                        fHistManager.FillTH1(histname10,eta_high);
								}
			fHistManager.FillTH3(histname11,DeltaPt,pt_low,eta_low);
			fHistManager.FillTH3(histname12,DeltaPt,pt_low,eta_high);
			fHistManager.FillTH3(histname13,DeltaPt,pt_low,(dynamic_cast<AliEmcalJet*>(LowRJetsList.At(j)))->GetNumberOfConstituents());
							}
					}
		}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
}			//End of JetCollArray loop

}

/*
 *This is a local implementation of the AliAnalysisHelperJetTasks helper task for the jet matching of the closest jets. There are two reasons for the locality of the implementation
 * a)Better control over the debugging procedure b)Adjust the type of the input lists from AliAODjet to AliEmcalJet.  
 */

void AliAnalysisTaskEmcalJetEnergyFlow::JetMatcher(const TList *genJetsList,const Int_t &kGenJets,
                                           const TList *recJetsList,const Int_t &kRecJets,
                                           TArrayI &iGenIndex,TArrayI &iRecIndex,
                                           Int_t iDebug, Float_t maxDist,Float_t max_eta){

   // Size indepnedendentt Implemenation of jet matching
   // Thepassed TArrayI should be static in the user function an only increased if needed
   //
   // Relate the two input jet Arrays
   // The association has to be unique
   // So check in two directions
   // find the closest rec to a gen
   // and check if there is no other rec which is closer
   // Caveat: Close low energy/split jets may disturb this correlation
  
   // Idea: search in two directions generated e.g (a--e) and rec (1--3)
   // Fill a matrix with Flags (1 for closest rec jet, 2 for closest rec jet
   // in the end we have something like this
   //    1   2   3
   // ------------
   // a| 3   2   0
   // b| 0   1   0
   // c| 0   0   3
   // d| 0   0   1
   // e| 0   0   1
   // Topology
   //   1     2
   //     a         b        
   //
   //  d      c
   //        3     e
   // Only entries with "3" match from both sides

	iGenIndex.Reset(-1);
	iRecIndex.Reset(-1);


	const int kMode = 3;
	const Int_t nGenJets = TMath::Min(genJetsList->GetEntries(),kGenJets);
 	const Int_t nRecJets = TMath::Min(recJetsList->GetEntries(),kRecJets);
        if(nRecJets==0||nGenJets==0)return;
 
        static TArrayS iFlag(nGenJets*nRecJets);
        if(iFlag.GetSize()<(nGenJets*nRecJets)){
        iFlag.Set(nGenJets*nRecJets);
        }
        iFlag.Reset(0);
  
       // find the closest distance to the generated
      for(int ig = 0;ig<nGenJets;++ig){
      AliEmcalJet *genJet = (AliEmcalJet*)genJetsList->At(ig);
    	if(!genJet)continue;
        if(fabs(genJet->Eta())>max_eta)continue;
      Float_t dist = maxDist;
      if(iDebug>1)Printf("Gen (%d) p_T %3.3f eta %3.3f ph %3.3f ",ig,genJet->Pt(),genJet->Eta(),genJet->Phi());
      for(int ir = 0;ir<nRecJets;++ir){
      AliEmcalJet *recJet = (AliEmcalJet*)recJetsList->At(ir);
      if(!recJet)continue;
      if(fabs(recJet->Eta())>max_eta)continue;
      if(iDebug>1){
      printf("Rec (%d) ",ir);
      Printf("p_T %3.3f eta %3.3f ph %3.3f ",recJet->Pt(),recJet->Eta(),recJet->Phi());
          }
        Double_t dR = genJet->DeltaR(recJet);
 	if(iDebug>1)Printf("Distance (%d)--(%d) %g ",ig,ir,dR);
	if(dR<dist){
    iRecIndex[ig] = ir;
    dist = dR;
        }
      }
      if(iRecIndex[ig]>=0)iFlag[ig*nRecJets+iRecIndex[ig]]+=1;
      // reset...
      iRecIndex[ig] = -1;
    }
    // other way around
       for(int ir = 0;ir<nRecJets;++ir){
	       AliEmcalJet *recJet = (AliEmcalJet*)recJetsList->At(ir);
 	      if(!recJet)continue;
	      if(fabs(recJet->Eta())>max_eta)continue;
 	       Float_t dist = maxDist;
 	       for(int ig = 0;ig<nGenJets;++ig){
  	   AliEmcalJet *genJet = (AliEmcalJet*)genJetsList->At(ig);
  	   if(!genJet)continue;
	   if(fabs(genJet->Eta())>max_eta)continue;
  	   Double_t dR = genJet->DeltaR(recJet);
	 if(dR<dist){
  	   iGenIndex[ir] = ig;
  	   dist = dR;
  	       }
  	     }
  	     if(iGenIndex[ir]>=0)iFlag[iGenIndex[ir]*nRecJets+ir]+=2;
 	     // reset...
  	     iGenIndex[ir] = -1;
  	   }
  
     // check for "true" correlations
  
     if(iDebug>1)Printf(">>>>>> Matrix Size %d",iFlag.GetSize());
  
     for(int ig = 0;ig<nGenJets;++ig){
       for(int ir = 0;ir<nRecJets;++ir){
         // Print
         if(iDebug>1)printf("Flag2[%d][%d] %d ",ig,ir,iFlag[ig*nRecJets+ir]);
  
         if(kMode==3){
    // we have a unique correlation
  	   if(iFlag[ig*nRecJets+ir]==3){
 	     iGenIndex[ir] = ig;
       	     iRecIndex[ig] = ir;
    	 }
         }
         else{
     // we just take the correlation from on side
  	   if((iFlag[ig*nRecJets+ir]&2)==2){
  	     iGenIndex[ir] = ig;
  	   }
  	   if((iFlag[ig*nRecJets+ir]&1)==1){
      	 iRecIndex[ig] = ir;
     }
  	       }
 	     }
  	     if(iDebug>1)printf("\n");
    	 }
 }




//*/
/*
 * This function allocates the histograms for basic tracking QA.
 * A set of histograms (pT, eta, phi, difference between kinematic properties
 * at the vertex and at the EMCal surface, number of tracks) is allocated
 * per each particle container and per each centrality bin.
 */

void AliAnalysisTaskEmcalJetEnergyFlow::AllocateTrackHistograms(){
  
  TString histname;
  TString histtitle;
  TString groupname;
  Int_t fNPtBins = 100;
  Double_t fMinPtBin = 0.0;
  Double_t fMaxPtBin = 100.0;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    groupname = partCont->GetName();
	// Protect against creating the histogram twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The track containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNPtBins / 2, fMinPtBin, fMaxPtBin / 2);

      histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{track};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNPtBins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{track};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNPtBins / 6, -1, 1);

      if (TClass(partCont->GetClassName()).InheritsFrom("AliVTrack")) {
        histname = TString::Format("%s/fHistDeltaEtaPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{#eta}_{track}^{vertex} - #it{#eta}_{track}^{EMCal};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNPtBins / 2, fMinPtBin, fMaxPtBin, 50, -0.5, 0.5);

        histname = TString::Format("%s/fHistDeltaPhiPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{#phi}_{track}^{vertex} - #it{#phi}_{track}^{EMCal};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNPtBins / 2, fMinPtBin, fMaxPtBin, 200, -2, 2);

        histname = TString::Format("%s/fHistDeltaPtvsPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{p}_{T,track}^{vertex} - #it{p}_{T,track}^{EMCal} (GeV/#it{c});counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNPtBins / 2, fMinPtBin, fMaxPtBin, fNPtBins / 2, -fMaxPtBin/2, fMaxPtBin/2);

        histname = TString::Format("%s/fHistEoverPvsP_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{P}_{track} (GeV/#it{c});#it{E}_{cluster} / #it{P}_{track} #it{c};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNPtBins / 2, fMinPtBin, fMaxPtBin, fNPtBins / 2, 0, 4);
      }

      histname = TString::Format("%s/histNTracks_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of tracks;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
      }
    }
  }

  histname = "fHistSumNTracks";
  histtitle = TString::Format("%s;Sum of n tracks;events", histname.Data());
  if (fForceBeamType != kpp) {
    fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000);
  }
  else {
    fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
  }
}

/**
 * This function performs a loop over the reconstructed tracks
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetEnergyFlow::DoTrackLoop(){

 AliClusterContainer* clusCont = GetClusterContainer(0);
  TString histname;
  TString groupname;
  UInt_t sumAcceptedTracks = 0;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    groupname = partCont->GetName();
    UInt_t count = 0;
    for(auto part : partCont->accepted()) {
      if (!part) continue;
      count++;

      histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, part->Pt());

      histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, part->Phi());

      histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, part->Eta());

      if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
        const AliVTrack* track = static_cast<const AliVTrack*>(part);

        histname = TString::Format("%s/fHistDeltaEtaPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, track->Pt(), track->Eta() - track->GetTrackEtaOnEMCal());

        histname = TString::Format("%s/fHistDeltaPhiPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, track->Pt(), track->Phi() - track->GetTrackPhiOnEMCal());

        histname = TString::Format("%s/fHistDeltaPtvsPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, track->Pt(), track->Pt() - track->GetTrackPtOnEMCal());

        if (clusCont) {
          Int_t iCluster = track->GetEMCALcluster();
          if (iCluster >= 0) {
            AliVCluster* cluster = clusCont->GetAcceptCluster(iCluster);
            if (cluster) {
              histname = TString::Format("%s/fHistEoverPvsP_%d", groupname.Data(), fCentBin);
              fHistManager.FillTH2(histname, track->P(), cluster->GetNonLinCorrEnergy() / track->P());
            }
          }
        }
      }
    }
    sumAcceptedTracks += count;

    histname = TString::Format("%s/histNTracks_%d", groupname.Data(), fCentBin);
    fHistManager.FillTH1(histname, count);
  }

  histname = "fHistSumNTracks";
  fHistManager.FillTH1(histname, sumAcceptedTracks);
}

/*
 * This function allocates the histograms for basic EMCal QA.
 * One 2D histogram with the cell energy spectra and the number of cells
 * per event is allocated per each centrality bin.
 */

void AliAnalysisTaskEmcalJetEnergyFlow::AllocateCellHistograms(){
 
 TString histname;
  TString histtitle;
  TString groupname(fCaloCellsName);

  fHistManager.CreateHistoGroup(groupname);
  for (Int_t cent = 0; cent < fNcentBins; cent++) {
    histname = TString::Format("%s/histCellEnergy_%d", groupname.Data(), cent);
    histtitle = TString::Format("%s;#it{E}_{cell} (GeV);counts", histname.Data());
    fHistManager.CreateTH1(histname, histtitle, 300, 0, 150);

    histname = TString::Format("%s/histNCells_%d", groupname.Data(), cent);
    histtitle = TString::Format("%s;number of cells;events", histname.Data());
    if (fForceBeamType != kpp) {
      fHistManager.CreateTH1(histname, histtitle, 500, 0, 6000);
    }
    else {
      fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
    }
  }
}

/**
 * This function performs a loop over the reconstructed EMCal cells
 * in the current event and fills the relevant histograms.
 */

void AliAnalysisTaskEmcalJetEnergyFlow::DoCellLoop(){

  if (!fCaloCells) return;

  TString histname;

  const Short_t ncells = fCaloCells->GetNumberOfCells();

  histname = TString::Format("%s/histNCells_%d", fCaloCellsName.Data(), fCentBin);
  fHistManager.FillTH1(histname, ncells);

  histname = TString::Format("%s/histCellEnergy_%d", fCaloCellsName.Data(), fCentBin);
  for (Short_t pos = 0; pos < ncells; pos++) {
    Double_t amp   = fCaloCells->GetAmplitude(pos);

    fHistManager.FillTH1(histname, amp);
  }

}
/*
 * This function allocates the histograms for basic EMCal cluster QA.
 * A set of histograms (energy, eta, phi, number of cluster) is allocated
 * per each cluster container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetEnergyFlow::AllocateClusterHistograms(){

  TString histname;
  TString histtitle;
  TString groupname;
  AliClusterContainer* clusCont = 0;
  Int_t fNPtBins = 100;
  Double_t fMinPtBin = 0.0;
  Double_t fMaxPtBin = 100.0;
  TIter next(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(next()))) {
    groupname = clusCont->GetName();
  // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The cluster containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/histClusterEnergy_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNPtBins / 2, fMinPtBin, fMaxPtBin / 2);

      histname = TString::Format("%s/histClusterEnergyExotic_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster}^{exotic} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNPtBins / 2, fMinPtBin, fMaxPtBin / 2);

      histname = TString::Format("%s/histClusterNonLinCorrEnergy_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster}^{non-lin.corr.} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNPtBins / 2, fMinPtBin, fMaxPtBin / 2);

      histname = TString::Format("%s/histClusterHadCorrEnergy_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster}^{had.corr.} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNPtBins / 2, fMinPtBin, fMaxPtBin / 2);

      histname = TString::Format("%s/histClusterPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{custer};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNPtBins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histClusterEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{custer};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNPtBins / 6, -1, 1);

      histname = TString::Format("%s/histNClusters_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of clusters;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 3000);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
      }
    }
  }

  histname = "fHistSumNClusters";
  histtitle = TString::Format("%s;Sum of n clusters;events", histname.Data());
  if (fForceBeamType != kpp) {
    fHistManager.CreateTH1(histname, histtitle, 500, 0, 3000);
  }
  else {
    fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
  }
}

/**
 * This function performs a loop over the reconstructed EMCal clusters
 * in the current event and fills the relevant histograms.
 */

void AliAnalysisTaskEmcalJetEnergyFlow::DoClusterLoop(){

  TString histname;
  TString groupname;
  UInt_t sumAcceptedClusters = 0; 
  AliClusterContainer* clusCont = 0;
  TIter next(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(next()))) {
    groupname = clusCont->GetName();

    for(auto cluster : clusCont->all()) {
      if (!cluster) continue;

      if (cluster->GetIsExotic()) {
        histname = TString::Format("%s/histClusterEnergyExotic_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, cluster->E());
      }
    }

    UInt_t count = 0;
    for(auto cluster : clusCont->accepted()) {
      if (!cluster) continue;
      count++;

      AliTLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);

      histname = TString::Format("%s/histClusterEnergy_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, cluster->E());

      histname = TString::Format("%s/histClusterNonLinCorrEnergy_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, cluster->GetNonLinCorrEnergy());

      histname = TString::Format("%s/histClusterHadCorrEnergy_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, cluster->GetHadCorrEnergy());

      histname = TString::Format("%s/histClusterPhi_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, nPart.Phi_0_2pi());

      histname = TString::Format("%s/histClusterEta_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, nPart.Eta());
    }
    sumAcceptedClusters += count;

    histname = TString::Format("%s/histNClusters_%d", groupname.Data(), fCentBin);
    fHistManager.FillTH1(histname, count);
  }

  histname = "fHistSumNClusters";
  fHistManager.FillTH1(histname, sumAcceptedClusters);

}



