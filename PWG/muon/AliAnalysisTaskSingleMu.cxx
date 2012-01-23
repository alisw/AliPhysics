/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//-----------------------------------------------------------------------------
/// \class AliAnalysisTaskSingleMu
/// Analysis task for single muons in the spectrometer.
/// The output is a list of histograms.
/// The macro class can run on AODs or ESDs.
/// In the latter case a flag can be activated to produce a tree as output.
/// If Monte Carlo information is present, some basics checks are performed.
///
/// \author Diego Stocco
//-----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//    Implementation of the single muon analysis class
// An example of usage can be found in the macro RunSingleMuAnalysisFromAOD.C.
//----------------------------------------------------------------------------

#define AliAnalysisTaskSingleMu_cxx

// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAxis.h"
#include "TList.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TMap.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TIterator.h"
#include "TParameter.h"
#include "TMCProcess.h"

// STEER includes
#include "AliInputEventHandler.h"
#include "AliVVertex.h"
#include "AliMultiplicity.h"
#include "AliCentrality.h"

//#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
//#include "AliAODVertex.h"

#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"

//#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

// CORRFW includes
#include "AliCFManager.h"
#include "AliCFContainer.h"
#include "AliCFGridSparse.h"
#include "AliCFTrackKineCuts.h"
#include "AliCFParticleGenCuts.h"
#include "AliCFEventRecCuts.h"

#include "AliAnalysisTaskSingleMu.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSingleMu) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliAnalysisTaskSingleMu::AliAnalysisTaskSingleMu(const char *name, Int_t fillTreeScaleDown, Bool_t keepAll) :
  AliAnalysisTaskSE(name),
  fFillTreeScaleDown(fillTreeScaleDown),
  fKeepAll(keepAll),
  fkNvtxContribCut(1),
  fTriggerClasses(0x0),
  fCFManager(0),
  fHistoList(0),
  fHistoListMC(0),
  fHistoListQA(0),
  fTreeSingleMu(0),
  fVarFloat(0),
  fVarInt(0),
  fVarChar(0),
  fVarUInt(0),
  fVarFloatMC(0),
  fVarIntMC(0),
  fAuxObjects(new TMap()),
  fDebugString("")
{
  //
  /// Constructor.
  //
  if ( fFillTreeScaleDown <= 0 )
    fKeepAll = kFALSE;

  DefineOutput(1, AliCFContainer::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());

  if ( fFillTreeScaleDown > 0 )
    DefineOutput(5, TTree::Class());

  fAuxObjects->SetOwner();

  SetTriggerClasses();
}


//________________________________________________________________________
AliAnalysisTaskSingleMu::~AliAnalysisTaskSingleMu()
{
  //
  /// Destructor
  //

  delete fTriggerClasses;
  // For proof: do not delete output containers
  if ( ! AliAnalysisManager::GetAnalysisManager() || ! AliAnalysisManager::GetAnalysisManager()->IsProofMode() ) {
    // In terminate mode fCFManager does not exist!
    // So, check before deleting
    if ( fCFManager )
      delete fCFManager->GetParticleContainer();  // The container is not deleted by framework
    delete fCFManager;
    delete fHistoList;
    delete fHistoListMC;
    delete fHistoListQA;
    delete fTreeSingleMu;
  }
  delete fVarFloat;
  delete fVarInt;
  delete [] fVarChar;
  delete fVarUInt;
  delete fVarFloatMC;
  delete fVarIntMC;
  delete fAuxObjects;
}


//___________________________________________________________________________
void AliAnalysisTaskSingleMu::NotifyRun()
{
  //
  /// Use the event handler information to correctly fill the analysis flags:
  /// - check if Monte Carlo information is present
  //  
  
  if ( fMCEvent ) {
    AliInfo("Monte Carlo information is present");
  }
  else {
    AliInfo("No Monte Carlo information in run");
  }
}


//___________________________________________________________________________
void AliAnalysisTaskSingleMu::FinishTaskOutput()
{
  //
  /// Perform histogram normalization after the last analyzed event
  /// but before merging
  //

  // Set the correct run limits for the histograms
  // vs run number to cope with the merging
  Int_t histoIndex = -1;
  Int_t indexPerRun[2] = {kHistoNeventsPerRun, kHistoNmuonsPerRun};
  for ( Int_t ihisto=0; ihisto<2; ihisto++) {
    histoIndex = GetHistoIndex(indexPerRun[ihisto]);
    TH2F* histo2D = (TH2F*)fHistoList->At(histoIndex);
    Double_t minX = 1e10, maxX = -1e10;
    for (Int_t ibin=1; ibin<=histo2D->GetXaxis()->GetNbins(); ibin++){
      TString runNum = histo2D->GetXaxis()->GetBinLabel(ibin);
      minX = TMath::Min(runNum.Atof()-0.5, minX);
      maxX = TMath::Max(runNum.Atof()+0.5, maxX);
    }
    histo2D->GetXaxis()->SetLimits(minX, maxX);
    AliInfo(Form("Histogram %s run limits (%f, %f)",histo2D->GetName(), minX, maxX));
  }

  // Add stat. info from physics selection
  // (usefull when running on AODs)
  TString histoName = "";
  if ( fInputHandler ) {
    for ( Int_t istat=0; istat<2; istat++ ) {
      TString statType = ( istat == 0 ) ? "ALL" : "BIN0";
      TH2* hStat = dynamic_cast<TH2*>(fInputHandler->GetStatistics(statType.Data()));
      if ( hStat ) {
        histoName = Form("%s_SingleMuon", hStat->GetName());
        TH2* cloneStat = static_cast<TH2*>(hStat->Clone(histoName.Data()));
        cloneStat->SetDirectory(0);
        fHistoList->Add(cloneStat);
      }
      else {
        AliWarning("Stat histogram not available");
        break;
      }
    } // loop on stat type
  }
}


//___________________________________________________________________________
void AliAnalysisTaskSingleMu::UserCreateOutputObjects() 
{
  //
  /// Create output histograms
  //
  AliInfo(Form("   CreateOutputObjects of task %s\n", GetName()));

  // initialize histogram lists
  fHistoList = new TList();
  fHistoList->SetOwner();
  fHistoListMC = new TList();
  fHistoListMC->SetOwner();
  fHistoListQA = new TList();
  fHistoListQA->SetOwner();

  // Init variables
  fVarFloat = new Float_t [kNvarFloat];
  fVarInt = new Int_t [kNvarInt];
  fVarChar = new Char_t *[kNvarChar];
  fVarUInt = new UInt_t [kNvarUInt];
  fVarFloatMC = new Float_t [kNvarFloatMC];
  fVarIntMC = new Int_t [kNvarIntMC];

  const Int_t charWidth[kNvarChar] = {255};
  for(Int_t ivar=0; ivar<kNvarChar; ivar++){
    fVarChar[ivar] = new Char_t [charWidth[ivar]];
  }

  Int_t nPtBins = 60; //200; //60;
  Double_t ptMin = 0., ptMax = 30.; //100.; //30.; // extend range for z
  TString ptName("Pt"), ptTitle("p_{t}"), ptUnits("GeV/c");

  Int_t nRapidityBins = 25;
  Double_t rapidityMin = -4.5, rapidityMax = -2.;
  //TString rapidityName("Rapidity"), rapidityTitle("y"), rapidityUnits("");
  TString rapidityName("Eta"), rapidityTitle("#eta"), rapidityUnits("");

  Int_t nPhiBins = 36;
  Double_t phiMin = 0.; Double_t phiMax = 2*TMath::Pi();
  TString phiName("Phi"), phiTitle("#phi"), phiUnits("rad");
  
  Int_t nDcaBins = 30;
  Double_t dcaMin = 0., dcaMax = 300.;
  TString dcaName("DCA"), dcaTitle("DCA"), dcaUnits("cm");

  Int_t nVzBins = 40;
  Double_t vzMin = -20., vzMax = 20.;
  TString vzName("Vz"), vzTitle("Vz"), vzUnits("cm");
  
  Int_t nThetaAbsEndBins = 4;
  Double_t thetaAbsEndMin = -0.5, thetaAbsEndMax = 3.5;
  TString thetaAbsEndName("ThetaAbsEnd"), thetaAbsEndTitle("#theta_{abs}"), thetaAbsEndUnits("a.u.");  

  Int_t nChargeBins = 2;
  Double_t chargeMin = -2, chargeMax = 2.;
  TString chargeName("Charge"), chargeTitle("charge"), chargeUnits("e");

  Int_t nMatchTrigBins = 4;
  Double_t matchTrigMin = -0.5, matchTrigMax = 3.5;
  TString matchTrigName("MatchTrig"), matchTrigTitle("Trigger match"), matchTrigUnits("");
  
  Int_t nTrigClassBins = fTriggerClasses->GetEntries();
  Double_t trigClassMin = 0.5, trigClassMax = (Double_t)nTrigClassBins + 0.5;
  TString trigClassName("TrigClass"), trigClassTitle("Fired trigger class"), trigClassUnits("");

  Int_t nGoodVtxBins = 3;
  Double_t goodVtxMin = -0.5, goodVtxMax = 2.5;
  TString goodVtxName("GoodVtx"), goodVtxTitle("Vertex flags"), goodVtxUnits("");

  Int_t nMotherTypeBins = kNtrackSources;
  Double_t motherTypeMin = -0.5, motherTypeMax = (Double_t)kNtrackSources - 0.5;
  TString motherType("MotherType"), motherTypeTitle("motherType"), motherTypeUnits("");

  Int_t nCentralityBins = 12;
  Double_t centralityMin = -10., centralityMax = 110.;
  TString centralityName("Centrality"), centralityTitle("centrality"), centralityUnits("%");

  TString trigName[kNtrigCuts];
  trigName[kNoMatchTrig] = "NoMatch";
  trigName[kAllPtTrig]   = "AllPt";
  trigName[kLowPtTrig]   = "LowPt";
  trigName[kHighPtTrig]  = "HighPt";

  TString srcName[kNtrackSources];
  srcName[kCharmMu]     = "Charm";
  srcName[kBeautyMu]    = "Beauty";
  srcName[kPrimaryMu]   = "Decay";
  srcName[kSecondaryMu] = "Secondary";
  srcName[kRecoHadron]  = "Hadrons";
  srcName[kUnknownPart] = "Unidentified";

  TH1F* histo1D = 0x0;
  TH2F* histo2D = 0x0;
  TString histoName, histoTitle;
  Int_t histoIndex = 0;

  // Multi-dimensional histo
  Int_t nbins[kNvars] = {nPtBins, nRapidityBins, nPhiBins, nDcaBins, nVzBins, nThetaAbsEndBins, nChargeBins, nMatchTrigBins, nTrigClassBins, nGoodVtxBins, nMotherTypeBins, nCentralityBins};
  Double_t xmin[kNvars] = {ptMin, rapidityMin, phiMin, dcaMin, vzMin, thetaAbsEndMin, chargeMin, matchTrigMin, trigClassMin, goodVtxMin, motherTypeMin, centralityMin};
  Double_t xmax[kNvars] = {ptMax, rapidityMax, phiMax, dcaMax, vzMax, thetaAbsEndMax, chargeMax, matchTrigMax, trigClassMax, goodVtxMax, motherTypeMax, centralityMax};
  TString axisTitle[kNvars] = {ptTitle, rapidityTitle, phiTitle, dcaTitle, vzTitle, thetaAbsEndTitle, chargeTitle, matchTrigTitle, trigClassTitle, goodVtxTitle, motherTypeTitle, centralityTitle};
  TString axisUnits[kNvars] = {ptUnits, rapidityUnits, phiUnits, dcaUnits, vzUnits, thetaAbsEndUnits, chargeUnits, matchTrigUnits, trigClassUnits, goodVtxUnits, motherTypeUnits, centralityUnits};

  //TString stepTitle[kNsteps] = {"reconstructed", "in acceptance", "generated", "in acceptance (MC)"};
  TString stepTitle[kNsteps] = {"reconstructed", "generated"};

  // Create CF container

  // The framework has problems if the name of the object
  // and the one of container differ
  // To be complaint, get the name from container and set it
  TString containerName = GetOutputSlot(1)->GetContainer()->GetName();
  
  AliCFContainer* container = new AliCFContainer(containerName.Data(),"container for tracks",kNsteps,kNvars,nbins);

  for ( Int_t idim = 0; idim<kNvars; idim++){
    histoTitle = Form("%s (%s)", axisTitle[idim].Data(), axisUnits[idim].Data());
    histoTitle.ReplaceAll("()","");

    container->SetVarTitle(idim, histoTitle.Data());
    container->SetBinLimits(idim, xmin[idim], xmax[idim]);
  }

  for (Int_t istep=0; istep<kNsteps; istep++){
    container->SetStepTitle(istep, stepTitle[istep].Data());
    AliCFGridSparse* gridSparse = container->GetGrid(istep);

    SetAxisLabel(gridSparse->GetAxis(kHvarTrigClass));
    TAxis* isGoodVtxAxis = gridSparse->GetAxis(kHvarIsGoodVtx);
    isGoodVtxAxis->SetBinLabel(1,Form("No vertex (contrib<%i)",fkNvtxContribCut));
    isGoodVtxAxis->SetBinLabel(2,"Good vertex");
    isGoodVtxAxis->SetBinLabel(3,"Pileup");
    TAxis* motherTypeAxis = gridSparse->GetAxis(kHvarMotherType);
    for (Int_t ibin=0; ibin<kNtrackSources; ibin++){
      motherTypeAxis->SetBinLabel(ibin+1,srcName[ibin]);
    }
  }

  // Create cuts

  //Particle-Level cuts:
  AliCFParticleGenCuts* mcGenCuts = new AliCFParticleGenCuts();
  mcGenCuts->SetNameTitle("mcGenCuts","MC particle generation cuts");
  mcGenCuts->SetRequirePdgCode(13,kTRUE);
  mcGenCuts->SetQAOn(fHistoListQA);

  /*
  // MC kinematic cuts
  AliCFTrackKineCuts *mcAccCuts = new AliCFTrackKineCuts();
  mcAccCuts->SetNameTitle("mcAccCuts","MC-level acceptance cuts");
  mcAccCuts->SetEtaRange(-4.,-2.5);
  mcAccCuts->SetQAOn(fHistoListQA);

  // Rec-Level kinematic cuts
  AliCFTrackKineCuts *recAccCuts = new AliCFTrackKineCuts();
  recAccCuts->SetNameTitle("recAccCuts","Reco-level acceptance cuts");
  recAccCuts->SetEtaRange(-4.,-2.5);
  recAccCuts->SetQAOn(fHistoListQA);
  */

  TObjArray* mcGenList = new TObjArray(0) ;
  mcGenList->AddLast(mcGenCuts);

  /*
  TObjArray* mcAccList = new TObjArray(0);
  mcAccList->AddLast(mcAccCuts);

  TObjArray* recAccList = new TObjArray(0);
  recAccList->AddLast(recAccCuts);
  */

  // Create CF manager
  fCFManager = new AliCFManager() ;
  fCFManager->SetParticleContainer(container);

  // Add cuts  
  // Dummy event container
  Int_t dummyBins[1] = {1};
  AliCFContainer* evtContainer = new AliCFContainer("dummyContainer","dummy contaier for events",1,1,dummyBins);
  fCFManager->SetEventContainer(evtContainer);
  fCFManager->SetEventCutsList(0,0x0);
 
  // Init empty cuts (avoid warnings in framework)
  for (Int_t istep=0; istep<kNsteps; istep++) {
    fCFManager->SetParticleCutsList(istep,0x0);    
  }
  //fCFManager->SetParticleCutsList(kStepAcceptance,recAccList);
  fCFManager->SetParticleCutsList(kStepGeneratedMC,mcGenList);
  //fCFManager->SetParticleCutsList(kStepAcceptanceMC,mcAccList);

  // Summary histos
  histoName = "histoNeventsPerTrig";
  histoTitle = "Number of events per trigger class";
  histo2D = new TH2F(histoName.Data(), histoTitle.Data(), nTrigClassBins, trigClassMin, trigClassMax, 5, 0.5, 5.5);
  histo2D->GetXaxis()->SetTitle("Fired class");
  SetAxisLabel(histo2D->GetXaxis());
  histo2D->GetYaxis()->SetBinLabel(1, "All events");
  histo2D->GetYaxis()->SetBinLabel(2, "Good vtx events");
  histo2D->GetYaxis()->SetBinLabel(3, "Pileup events");
  histo2D->GetYaxis()->SetBinLabel(4, "All vertices");
  histo2D->GetYaxis()->SetBinLabel(5, "Pileup vertices");
  histoIndex = GetHistoIndex(kHistoNeventsPerTrig);
  fHistoList->AddAt(histo2D, histoIndex);

  histoName = "histoMuonMultiplicity";
  histoTitle = "Muon track multiplicity";
  histo2D = new TH2F(histoName.Data(), histoTitle.Data(), 15, -0.5, 15-0.5,
		     nTrigClassBins, trigClassMin, trigClassMax);
  histo2D->GetXaxis()->SetTitle("# of muons");
  SetAxisLabel(histo2D->GetYaxis());
  histoIndex = GetHistoIndex(kHistoMuonMultiplicity);
  fHistoList->AddAt(histo2D, histoIndex);

  histoName = "histoEventVz";
  histoTitle = "All events IP Vz distribution";
  histo2D = new TH2F(histoName.Data(), histoTitle.Data(), nTrigClassBins, trigClassMin, trigClassMax, nVzBins, vzMin, vzMax);
  histo2D->GetXaxis()->SetTitle("Fired class");
  SetAxisLabel(histo2D->GetXaxis());
  histoTitle = Form("%s (%s)", vzTitle.Data(), vzUnits.Data());
  histo2D->GetYaxis()->SetTitle(histoTitle.Data());
  histoIndex = GetHistoIndex(kHistoEventVz);
  fHistoList->AddAt(histo2D, histoIndex);

  Int_t hRunIndex[2] = {kHistoNeventsPerRun, kHistoNmuonsPerRun};
  TString hRunName[2] = {"histoNeventsPerRun", "histoNmuonsPerRun"};
  TString hRunTitle[2] = {"Number of events per run", "Number of muons per run"};
  for ( Int_t ihisto=0; ihisto<2; ihisto++){
    histoName = hRunName[ihisto];
    histoTitle = hRunTitle[ihisto];
    histo2D = new TH2F(histoName.Data(), histoTitle.Data(),
		       1, 0., 0.,
		       nTrigClassBins, trigClassMin, trigClassMax);
    histo2D->GetXaxis()->SetTitle("Run number");
    SetAxisLabel(histo2D->GetYaxis());
    histo2D->Sumw2();
    histoIndex = hRunIndex[ihisto];
    fHistoList->AddAt(histo2D, histoIndex);
  }

  // MC histos summary
  Int_t hCheckVzIndex[3] = {kHistoCheckVzMC, kHistoCheckVzHasVtxMC, kHistoCheckVzNoPileupMC};
  TString hCheckVzName[3] = {"histoCheckVz", "histoCheckVzHasVtx", "histoCheckVzIsPileup"};
  TString hCheckVzTitle[3] = {"", " w/ vertex contributors", "w/ pileup SPD"};

  for ( Int_t ihisto=0; ihisto<3; ihisto++ ) {
    histoName = hCheckVzName[ihisto];
    histoTitle = Form("Check IP Vz distribution %s", hCheckVzTitle[ihisto].Data());

    histo2D = new TH2F(histoName.Data(), histoTitle.Data(),
		       nVzBins, vzMin, vzMax,
		       nVzBins, vzMin, vzMax);
    histoTitle = Form("%s (%s)", vzTitle.Data(), vzUnits.Data());
    histo2D->GetXaxis()->SetTitle(histoTitle.Data());
    histoTitle = Form("%s MC (%s)", vzTitle.Data(), vzUnits.Data());
    histo2D->GetYaxis()->SetTitle(histoTitle.Data());
    histoIndex = GetHistoIndex(hCheckVzIndex[ihisto]);
    fHistoListMC->AddAt(histo2D, histoIndex);
  }

  // MC histos
  for (Int_t itrig=0; itrig < kNtrigCuts; itrig++) {
    for (Int_t isrc = 0; isrc < kNtrackSources; isrc++) {
      histoName = Form("%sResolution%sTrig%s", ptName.Data(), trigName[itrig].Data(), srcName[isrc].Data());
      histoTitle = Form("%s resolution. Trigger: %s (%s)", ptTitle.Data(), trigName[itrig].Data(), srcName[isrc].Data());
      histo1D = new TH1F(histoName.Data(), histoTitle.Data(), 100, -5., 5.);
      histo1D->GetXaxis()->SetTitle(Form("%s^{reco} - %s^{MC} (%s)", ptTitle.Data(), ptTitle.Data(), ptUnits.Data()));
      histoIndex = GetHistoIndex(kHistoPtResolutionMC, itrig, isrc);
      fHistoListMC->AddAt(histo1D, histoIndex);
    }
  }

  // Trees
  if ( fFillTreeScaleDown > 0 ) {
    TString leavesFloat[kNvarFloat] = {"Px", "Py", "Pz", "Pt",
				       "PxAtDCA", "PyAtDCA", "PzAtDCA", "PtAtDCA",
				       "PxUncorrected", "PyUncorrected", "PzUncorrected", "PtUncorrected",
				       "XUncorrected", "YUncorrected", "ZUncorrected",
				       "XatDCA", "YatDCA", "DCA",
				       "Eta", "Rapidity", "Charge", "RAtAbsEnd",
				       "IPVx", "IPVy", "IPVz", "Centrality"};
    TString leavesInt[kNvarInt] = {"MatchTrig", "IsMuon", "IsGhost", "LoCircuit", "PassPhysicsSelection", "NVtxContrib", "NspdTracklets", "IsPileupVertex"};
    TString leavesChar[kNvarChar] = {"FiredTrigClass"};
    TString leavesUInt[kNvarUInt] = {"BunchCrossNum", "OrbitNum", "PeriodNum", "RunNum"};
    TString leavesFloatMC[kNvarFloatMC] = {"PxMC", "PyMC", "PzMC", "PtMC",
					   "EtaMC", "RapidityMC",
					   "VxMC", "VyMC", "VzMC",
					   "MotherPxMC", "MotherPyMC", "MotherPzMC",
					   "MotherEtaMC", "MotherRapidityMC",
					   "MotherVxMC", "MotherVyMC", "MotherVzMC",
					   "IPVxMC", "IPVyMC", "IPVzMC"};
    TString leavesIntMC[kNvarIntMC] = {"Pdg", "MotherPdg", "MotherType"};

    TString treeName = GetOutputSlot(5)->GetContainer()->GetName();
    Bool_t hasMC = AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
    TString treeTitle = ( hasMC ) ? "Single Mu" : "Single Mu MC";

    OpenFile(5);
    if ( ! fTreeSingleMu ) fTreeSingleMu = new TTree(treeName.Data(), treeTitle.Data());

    for(Int_t itree=0; itree<2; itree++){
      TParameter<Int_t>* par1 = new TParameter<Int_t>("fillTreeScaleDown",fFillTreeScaleDown);
      TParameter<Int_t>* par2 = new TParameter<Int_t>("keepAllEvents",fKeepAll);
      fTreeSingleMu->GetUserInfo()->Add(par1);
      fTreeSingleMu->GetUserInfo()->Add(par2);
      for(Int_t ivar=0; ivar<kNvarFloat; ivar++){
	fTreeSingleMu->Branch(leavesFloat[ivar].Data(), &fVarFloat[ivar], Form("%s/F", leavesFloat[ivar].Data()));
      }
      for(Int_t ivar=0; ivar<kNvarInt; ivar++){
	fTreeSingleMu->Branch(leavesInt[ivar].Data(), &fVarInt[ivar], Form("%s/I", leavesInt[ivar].Data()));
      }
      for(Int_t ivar=0; ivar<kNvarChar; ivar++){
	TString addString = leavesChar[ivar] + "/C";
	fTreeSingleMu->Branch(leavesChar[ivar].Data(), fVarChar[ivar], addString.Data());
      }
      for(Int_t ivar=0; ivar<kNvarUInt; ivar++){
	fTreeSingleMu->Branch(leavesUInt[ivar].Data(), &fVarUInt[ivar], Form("%s/i", leavesUInt[ivar].Data()));
      }
      if ( hasMC ){
	for(Int_t ivar=0; ivar<kNvarFloatMC; ivar++){
	  fTreeSingleMu->Branch(leavesFloatMC[ivar].Data(), &fVarFloatMC[ivar], Form("%s/F", leavesFloatMC[ivar].Data()));
	}
	for(Int_t ivar=0; ivar<kNvarIntMC; ivar++){
	  fTreeSingleMu->Branch(leavesIntMC[ivar].Data(), &fVarIntMC[ivar], Form("%s/I", leavesIntMC[ivar].Data()));
	}
      }
    } // loop on trees
    PostData(5, fTreeSingleMu);
  } // fillNTuple
  
  PostData(1, fCFManager->GetParticleContainer());
  PostData(2, fHistoList);
  PostData(3, fHistoListQA);
  PostData(4, fHistoListMC);

}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::UserExec(Option_t * /*option*/) 
{
  //
  /// Main loop
  /// Called for each event
  //

  AliESDEvent* esdEvent = 0x0;
  AliAODEvent* aodEvent = 0x0;

  esdEvent = dynamic_cast<AliESDEvent*> (InputEvent());
  if ( ! esdEvent ){
    aodEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  }

  if ( ! aodEvent && ! esdEvent ) {
    AliError ("AOD or ESD event not found. Nothing done!");
    return;
  }

  if ( ! fMCEvent && InputEvent()->GetEventType() != 7 ) return; // Run only on physics events!

  Bool_t fillCurrentEventTree = ( fFillTreeScaleDown == 0 ) ? kFALSE : ( Entry() % fFillTreeScaleDown == 0 );

  Reset(kFALSE);

  //
  // Global event info
  //
  TString firedTrigClasses = ( esdEvent ) ? esdEvent->GetFiredTriggerClasses() : aodEvent->GetFiredTriggerClasses();
  fVarFloat[kVarCentrality] = InputEvent()->GetCentrality()->GetCentralityPercentile("V0M");

  AliVVertex* primaryVertex = ( esdEvent ) ? (AliVVertex*)esdEvent->GetPrimaryVertexSPD() : (AliVVertex*)aodEvent->GetPrimaryVertexSPD();

  fVarFloat[kVarIPVz] = primaryVertex->GetZ();
  fVarInt[kVarNVtxContrib] = primaryVertex->GetNContributors();
  fVarInt[kVarIsPileup] = ( esdEvent ) ? esdEvent->IsPileupFromSPDInMultBins() : aodEvent->IsPileupFromSPDInMultBins();
  fVarUInt[kVarRunNumber] = InputEvent()->GetRunNumber();

  Int_t isGoodVtxBin = ( fVarInt[kVarNVtxContrib] >= fkNvtxContribCut );
  if ( isGoodVtxBin && fVarInt[kVarIsPileup] > 0 )
    isGoodVtxBin = 2;

  if ( fillCurrentEventTree ){
    strncpy(fVarChar[kVarTrigMask], firedTrigClasses.Data(),255);
    fVarInt[kVarPassPhysicsSelection] = fInputHandler->IsEventSelected();

    // Small workaround: in MC the bunch ID are not properly set and the timestamp is in seconds
    // So fill bunchCrossing with the read timestamp
    //    fill the orbit and period number with a timestamp created while reading the run
    TTimeStamp ts;
    fVarUInt[kVarBunchCrossNumber] = ( fMCEvent ) ? (UInt_t)Entry() : InputEvent()->GetBunchCrossNumber();
    fVarUInt[kVarOrbitNumber] = ( fMCEvent ) ? (UInt_t)ts.GetNanoSec() : InputEvent()->GetOrbitNumber();
    fVarUInt[kVarPeriodNumber] = ( fMCEvent ) ? ts.GetTime() : InputEvent()->GetPeriodNumber();

    fVarFloat[kVarIPVx] = primaryVertex->GetX();
    fVarFloat[kVarIPVy] = primaryVertex->GetY();
    if ( esdEvent )
      fVarInt[kVarNspdTracklets] = esdEvent->GetMultiplicity()->GetNumberOfTracklets();
    else if ( aodEvent->GetTracklets() )
      fVarInt[kVarNspdTracklets] = aodEvent->GetTracklets()->GetNumberOfTracklets();
  }

  firedTrigClasses.Append(" ANY");

  // Object declaration
  AliMCParticle* mcPart = 0x0;
  AliVParticle* track = 0x0;

  Double_t containerInput[kNvars];
  Int_t histoIndex = -1;

  fCFManager->SetRecEventInfo(InputEvent());

  //
  // Pure Monte Carlo part
  //
  if ( fMCEvent ) {
    fCFManager->SetMCEventInfo (fMCEvent);
    Int_t nMCtracks = fMCEvent->GetNumberOfTracks();
    if ( nMCtracks > 0 ) {
      fVarFloatMC[kVarIPVxMC] = fMCEvent->Stack()->Particle(0)->Vx();
      fVarFloatMC[kVarIPVyMC] = fMCEvent->Stack()->Particle(0)->Vy();
      fVarFloatMC[kVarIPVzMC] = fMCEvent->Stack()->Particle(0)->Vz();
      containerInput[kHvarVz] = fVarFloatMC[kVarIPVzMC];
      containerInput[kHvarIsGoodVtx] = 1;
      containerInput[kHvarCentrality] = fVarFloat[kVarCentrality];
      histoIndex = GetHistoIndex(kHistoCheckVzMC);
      ((TH2F*)fHistoListMC->At(histoIndex))->Fill(fVarFloat[kVarIPVz], fVarFloatMC[kVarIPVzMC]);
      if ( isGoodVtxBin >= 1 ) {
	histoIndex = GetHistoIndex(kHistoCheckVzHasVtxMC);
	((TH2F*)fHistoListMC->At(histoIndex))->Fill(fVarFloat[kVarIPVz], fVarFloatMC[kVarIPVzMC]);
      }
      if ( isGoodVtxBin == 1 ) {
	histoIndex = GetHistoIndex(kHistoCheckVzNoPileupMC);
	((TH2F*)fHistoListMC->At(histoIndex))->Fill(fVarFloat[kVarIPVz], fVarFloatMC[kVarIPVzMC]);
      }
    }

    for (Int_t ipart=0; ipart<nMCtracks; ipart++) {
      mcPart = (AliMCParticle*)fMCEvent->GetTrack(ipart);

      //check the MC-level cuts
      if ( ! fCFManager->CheckParticleCuts(kStepGeneratedMC,mcPart) )
	continue;

      containerInput[kHvarPt]  = mcPart->Pt();
      //containerInput[kHvarY]   = mcPart->Y();
      containerInput[kHvarEta]   = mcPart->Eta();
      containerInput[kHvarPhi] = mcPart->Phi();
      containerInput[kHvarDCA] = TMath::Sqrt(mcPart->Xv()*mcPart->Xv() +
					     mcPart->Yv()*mcPart->Yv());
      containerInput[kHvarThetaZones] = GetBinThetaAbsEnd(TMath::Pi()-mcPart->Theta(),kTRUE);
      containerInput[kHvarCharge] = mcPart->Charge()/3.;
      containerInput[kHvarMatchTrig] = 1.;
      containerInput[kHvarMotherType] = (Double_t)RecoTrackMother(mcPart);
      //containerInput[kHvarP] = mcPart->P();

      for ( Int_t itrig=0; itrig<fTriggerClasses->GetEntries(); itrig++ ) {
	TString trigName = ((TObjString*)fTriggerClasses->At(itrig))->GetString();
	if ( ! firedTrigClasses.Contains(trigName.Data()) )
	  continue;
	containerInput[kHvarTrigClass] = (Double_t)(itrig+1);
	fCFManager->GetParticleContainer()->Fill(containerInput,kStepGeneratedMC);
	// if ( fCFManager->CheckParticleCuts(kStepAcceptanceMC,mcPart) ) fCFManager->GetParticleContainer()->Fill(containerInput,kStepAcceptanceMC);
      } // loop on trigger classes
      if ( fDebug >= 2 ) printf("AliAnalysisTaskSingleMu: Pure MC. %s. Set mother %i\n", fDebugString.Data(), (Int_t)containerInput[kHvarMotherType]);
    } // loop on MC particles
  } // is MC


  //
  // Reconstruction part
  //
  Int_t trackLabel = -1;
  Bool_t isGhost = kFALSE;
  Int_t nGhosts = 0, nMuons = 0;

  Int_t nTracks = ( esdEvent ) ? esdEvent->GetNumberOfMuonTracks() : aodEvent->GetNTracks();

  for (Int_t itrack = 0; itrack < nTracks; itrack++) {
    if ( esdEvent ){
      track = esdEvent->GetMuonTrack(itrack);
      isGhost = ( ((AliESDMuonTrack*)track)->ContainTriggerData() && ! ((AliESDMuonTrack*)track)->ContainTrackerData() );
    }
    else {
      track = aodEvent->GetTrack(itrack);
      if ( ! ((AliAODTrack*)track)->IsMuonTrack() )
	continue;
    }    

    if ( isGhost ) {
      ++nGhosts;
      // Nothing to do for ghosts if the tree is not filled
      if ( ! fillCurrentEventTree ) continue;
    }
    else ++nMuons;

    fVarFloat[kVarPt] = track->Pt();
    //fVarFloat[kVarRapidity] = ( isGhost ) ? 0. : track->Y();
    fVarFloat[kVarEta] = ( isGhost ) ? 0. : track->Eta();
    fVarFloat[kVarXatDCA] = ( esdEvent ) ? ((AliESDMuonTrack*)track)->GetNonBendingCoorAtDCA() : ((AliAODTrack*)track)->XAtDCA();
    fVarFloat[kVarYatDCA] = ( esdEvent ) ? ((AliESDMuonTrack*)track)->GetBendingCoorAtDCA() : ((AliAODTrack*)track)->YAtDCA();
    fVarFloat[kVarDCA] = 
      TMath::Sqrt( fVarFloat[kVarXatDCA] * fVarFloat[kVarXatDCA] +
		   fVarFloat[kVarYatDCA] * fVarFloat[kVarYatDCA] );
    fVarFloat[kVarCharge] = ( isGhost ) ? 0. : (Float_t)track->Charge();
    fVarFloat[kVarRAtAbsEnd] = ( esdEvent ) ? ((AliESDMuonTrack*)track)->GetRAtAbsorberEnd() : ((AliAODTrack*)track)->GetRAtAbsorberEnd();
    fVarInt[kVarMatchTrig] = ( esdEvent ) ? ((AliESDMuonTrack*)track)->GetMatchTrigger() : ((AliAODTrack*)track)->GetMatchTrigger();

    fVarIntMC[kVarMotherType] = kUnknownPart;

    // Monte Carlo part
    if ( fMCEvent ) {
      trackLabel = track->GetLabel();
      if ( trackLabel >= 0 ) {
	AliMCParticle* matchedMCTrack = (AliMCParticle*)fMCEvent->GetTrack(trackLabel);
	fVarIntMC[kVarMotherType] = RecoTrackMother(matchedMCTrack);
	fVarFloatMC[kVarPtMC] = matchedMCTrack->Pt();
	if ( ! isGhost ) FillTriggerHistos(kHistoPtResolutionMC, fVarInt[kVarMatchTrig], fVarIntMC[kVarMotherType], fVarFloat[kVarPt] - fVarFloatMC[kVarPtMC]);
	if ( fillCurrentEventTree ){
	  fVarFloatMC[kVarPxMC] = matchedMCTrack->Px();
	  fVarFloatMC[kVarPyMC] = matchedMCTrack->Py();
	  fVarFloatMC[kVarPzMC] = matchedMCTrack->Pz();
	  fVarFloatMC[kVarEtaMC] = matchedMCTrack->Eta();
	  fVarFloatMC[kVarRapidityMC] = matchedMCTrack->Y();
	  fVarFloatMC[kVarVxMC] = matchedMCTrack->Xv();
	  fVarFloatMC[kVarVyMC] = matchedMCTrack->Yv();
	  fVarFloatMC[kVarVzMC] = matchedMCTrack->Zv();
	  fVarIntMC[kVarPdg] = matchedMCTrack->PdgCode();

	  Int_t imother = matchedMCTrack->GetMother();
	  if ( imother >= 0 ) {
	    AliMCParticle* motherTrack = (AliMCParticle*)fMCEvent->GetTrack(imother);
	    fVarFloatMC[kVarMotherPxMC] = motherTrack->Px();
	    fVarFloatMC[kVarMotherPyMC] = motherTrack->Py();
	    fVarFloatMC[kVarMotherPzMC] = motherTrack->Pz();
	    fVarFloatMC[kVarMotherEtaMC] = motherTrack->Eta();
	    fVarFloatMC[kVarMotherRapidityMC] = motherTrack->Y();
	    fVarFloatMC[kVarMotherVxMC] = motherTrack->Xv();
	    fVarFloatMC[kVarMotherVyMC] = motherTrack->Yv();
	    fVarFloatMC[kVarMotherVzMC] = motherTrack->Zv();
	    fVarIntMC[kVarMotherPdg] = motherTrack->PdgCode();
	  }
	}
      }
      if ( fDebug >= 1 ) printf("AliAnalysisTaskSingleMu: Reco track. %s. Set mother %i\n", fDebugString.Data(), fVarIntMC[kVarMotherType]);
    } // if use MC

    if ( fillCurrentEventTree ) {
      fVarInt[kVarIsMuon] = ( isGhost ) ? 0 : nMuons;
      fVarInt[kVarIsGhost] = ( isGhost ) ? nGhosts : 0;
      //fVarFloat[kVarEta] = ( isGhost ) ? 0. : track->Eta();
      fVarFloat[kVarRapidity] = ( isGhost ) ? 0.: track->Y();
      fVarFloat[kVarPx] = track->Px();
      fVarFloat[kVarPy] = track->Py();
      fVarFloat[kVarPz] = track->Pz();
      fVarFloat[kVarPxAtDCA] = ( esdEvent ) ? ((AliESDMuonTrack*)track)->PxAtDCA() : ((AliAODTrack*)track)->PxAtDCA();
      fVarFloat[kVarPyAtDCA] = ( esdEvent ) ? ((AliESDMuonTrack*)track)->PyAtDCA() : ((AliAODTrack*)track)->PyAtDCA();
      fVarFloat[kVarPzAtDCA] = ( esdEvent ) ? ((AliESDMuonTrack*)track)->PzAtDCA() : ((AliAODTrack*)track)->PzAtDCA();
      fVarFloat[kVarPtAtDCA] = TMath::Sqrt(fVarFloat[kVarPxAtDCA]*fVarFloat[kVarPxAtDCA] + fVarFloat[kVarPyAtDCA]*fVarFloat[kVarPyAtDCA]);

      // Information present only on ESD tracks
      if ( esdEvent ) {
        fVarFloat[kVarPxUncorrected] = ( isGhost ) ? -TMath::Tan(((AliESDMuonTrack*)track)->GetThetaXUncorrected()) : ((AliESDMuonTrack*)track)->PxUncorrected();
        fVarFloat[kVarPyUncorrected] = ( isGhost ) ? -TMath::Tan(((AliESDMuonTrack*)track)->GetThetaYUncorrected()) : ((AliESDMuonTrack*)track)->PyUncorrected();
        fVarFloat[kVarPzUncorrected] = ( isGhost ) ? -1 : ((AliESDMuonTrack*)track)->PzUncorrected();

        fVarFloat[kVarPtUncorrected] = 
          TMath::Sqrt(fVarFloat[kVarPxUncorrected] * fVarFloat[kVarPxUncorrected] + 
                      fVarFloat[kVarPyUncorrected] * fVarFloat[kVarPyUncorrected]);

        fVarFloat[kVarXUncorrected] = ((AliESDMuonTrack*)track)->GetNonBendingCoorUncorrected();
        fVarFloat[kVarYUncorrected] = ((AliESDMuonTrack*)track)->GetBendingCoorUncorrected();
        fVarFloat[kVarZUncorrected] = ((AliESDMuonTrack*)track)->GetZUncorrected();

        fVarInt[kVarLocalCircuit] = ((AliESDMuonTrack*)track)->LoCircuit();
      }

      fTreeSingleMu->Fill();
    }

    if ( isGhost ) continue;
    
    //
    // Fill container
    //
    containerInput[kHvarPt]  = fVarFloat[kVarPt];
    //containerInput[kHvarY]   = fVarFloat[kVarRapidity];
    containerInput[kHvarEta]   = fVarFloat[kVarEta];
    containerInput[kHvarPhi] = track->Phi();
    containerInput[kHvarDCA] = fVarFloat[kVarDCA];
    containerInput[kHvarVz]  = fVarFloat[kVarIPVz];
    containerInput[kHvarThetaZones] = GetBinThetaAbsEnd(fVarFloat[kVarRAtAbsEnd]);
    containerInput[kHvarCharge] = fVarFloat[kVarCharge];
    containerInput[kHvarMatchTrig] = (Double_t)fVarInt[kVarMatchTrig];
    containerInput[kHvarIsGoodVtx] = (Double_t)isGoodVtxBin;
    containerInput[kHvarMotherType] = (Double_t)fVarIntMC[kVarMotherType];
    containerInput[kHvarCentrality] = fVarFloat[kVarCentrality];
    //containerInput[kHvarP] = track->P();

    histoIndex = GetHistoIndex(kHistoNmuonsPerRun);
    for ( Int_t itrig=0; itrig<fTriggerClasses->GetEntries(); itrig++ ) {
      TString trigName = ((TObjString*)fTriggerClasses->At(itrig))->GetString();
      if ( ! firedTrigClasses.Contains(trigName.Data()) )
	continue;
      containerInput[kHvarTrigClass] = (Double_t)(itrig+1);
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed);
      // if ( fCFManager->CheckParticleCuts(kStepAcceptance,track) ) fCFManager->GetParticleContainer()->Fill(containerInput,kStepAcceptance);
      ((TH2F*)fHistoList->At(histoIndex))->Fill(Form("%u",fVarUInt[kVarRunNumber]), containerInput[kHvarTrigClass], 1.);
    } // loop on trigger classes
  } // loop on tracks

  //
  // Complete global information 
  //
  if ( fillCurrentEventTree && fKeepAll &&  ( ( nMuons + nGhosts ) == 0 ) ) {
    // Fill event also if there is not muon (when explicitely required)
    fTreeSingleMu->Fill();
  }

  for ( Int_t itrig=0; itrig<fTriggerClasses->GetEntries(); itrig++ ) {
    TString trigName = ((TObjString*)fTriggerClasses->At(itrig))->GetString();
    if ( ! firedTrigClasses.Contains(trigName.Data()) )
      continue;
    Double_t trigClassBin = (Double_t)(itrig+1);
    histoIndex = GetHistoIndex(kHistoNeventsPerTrig);
    ((TH2F*)fHistoList->At(histoIndex))->Fill(trigClassBin, 1., 1.); // All events
    ((TH2F*)fHistoList->At(histoIndex))->Fill(trigClassBin, 4., (Double_t)(fVarInt[kVarIsPileup]+1)); // All vertices
    if ( isGoodVtxBin == 1 )
      ((TH2F*)fHistoList->At(histoIndex))->Fill(trigClassBin, 2., 1.); // Good vtx events
    else if ( isGoodVtxBin == 2 ) {
      ((TH2F*)fHistoList->At(histoIndex))->Fill(trigClassBin, 3., 1.); // Pileup events
      ((TH2F*)fHistoList->At(histoIndex))->Fill(trigClassBin, 5., 1.); // Pileup vertices
    }

    histoIndex = GetHistoIndex(kHistoMuonMultiplicity);
    ((TH1F*)fHistoList->At(histoIndex))->Fill(nMuons, trigClassBin);
    
    if ( isGoodVtxBin == 1 ) {
      histoIndex = GetHistoIndex(kHistoEventVz);
      ((TH2F*)fHistoList->At(histoIndex))->Fill(trigClassBin,fVarFloat[kVarIPVz]);
    }

    histoIndex = GetHistoIndex(kHistoNeventsPerRun);
    ((TH2F*)fHistoList->At(histoIndex))->Fill(Form("%u",fVarUInt[kVarRunNumber]), trigClassBin, 1.);
  } // loop on trigger classes


  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fCFManager->GetParticleContainer());
  PostData(2, fHistoList);
  PostData(3, fHistoListQA);
  PostData(4, fHistoListMC);
  if ( fFillTreeScaleDown > 0 )
    PostData(5, fTreeSingleMu);
}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::Terminate(Option_t *) {
  //
  /// Draw some histogram at the end.
  //

  AliCFContainer* container = dynamic_cast<AliCFContainer*> (GetOutputData(1));
  if ( ! container ) {
    AliError("Cannot find container in file");
    return;
  }

  //Int_t histoIndex = -1;

  if ( ! gROOT->IsBatch() ) {
    TString currName = GetName();
    currName.Prepend("c1_");
    TCanvas *c1_SingleMu = new TCanvas(currName.Data(),"Vz vs Pt",10,10,310,310);
    c1_SingleMu->SetFillColor(10); c1_SingleMu->SetHighLightColor(10);
    c1_SingleMu->SetLeftMargin(0.15); c1_SingleMu->SetBottomMargin(0.15);
    TH2* histo = static_cast<TH2*>(container->Project(kStepReconstructed,kHvarPt,kHvarVz));
    currName = GetName();
    currName.Prepend("hPtVz_");
    histo->SetName(currName.Data());
    histo->Draw("COLZ");
  }

}

//________________________________________________________________________
 void AliAnalysisTaskSingleMu::FillTriggerHistos(Int_t histoIndex, Int_t matchTrig, Int_t motherType,
						Float_t var1, Float_t var2, Float_t var3)
{
  //
  /// Fill all histograms passing the trigger cuts
  //

  Int_t nTrigs = TMath::Max(1, matchTrig);
  TArrayI trigToFill(nTrigs);
  trigToFill[0] = matchTrig;
  for(Int_t itrig = 1; itrig < matchTrig; itrig++){
    trigToFill[itrig] = itrig;
  }

  TList* histoList = (motherType < 0 ) ? fHistoList : fHistoListMC;

  TString className;
  for(Int_t itrig = 0; itrig < nTrigs; itrig++){
    Int_t currIndex = GetHistoIndex(histoIndex, trigToFill[itrig], motherType);
    className = histoList->At(currIndex)->ClassName();
    if ( fDebug >= 3 ) printf("AliAnalysisTaskSingleMu: matchTrig %i  Fill %s\n", matchTrig, histoList->At(currIndex)->GetName());
    if (className.Contains("1"))
      ((TH1F*)histoList->At(currIndex))->Fill(var1);
    else if (className.Contains("2"))
      ((TH2F*)histoList->At(currIndex))->Fill(var1, var2);
    else if (className.Contains("3"))
      ((TH3F*)histoList->At(currIndex))->Fill(var1, var2, var3);
    else
      AliWarning("Histogram not filled");
  }
}

//________________________________________________________________________
Int_t AliAnalysisTaskSingleMu::RecoTrackMother(AliMCParticle* mcParticle)
{
  //
  /// Find track mother from kinematics
  //

  Int_t recoPdg = mcParticle->PdgCode();

  fDebugString = Form("History: %i", recoPdg);

  // Track is not a muon
  if ( TMath::Abs(recoPdg) != 13 ) return kRecoHadron;

  Int_t imother = mcParticle->GetMother();

  //Int_t baseFlv[2] = {4,5};
  //Int_t mType[2] = {kCharmMu, kBeautyMu};
  Int_t den[3] = {100, 1000, 1};

  //Bool_t isFirstMotherHF = kFALSE;
  //Int_t step = 0;
  Int_t motherType = kPrimaryMu;
  while ( imother >= 0 ) {
    TParticle* part = ((AliMCParticle*)fMCEvent->GetTrack(imother))->Particle();

    fDebugString += Form(" <- %s (%s)", part->GetName(), TMCProcessName[part->GetUniqueID()]);

    Int_t absPdg = TMath::Abs(part->GetPdgCode());
    //step++;

    if ( imother < fMCEvent->GetNumberOfPrimaries() ) {
      /*
      // Hadronic processes are not possible for "primary" =>
      // either is decay or HF
      if ( absPdg > 100 && absPdg < 400 ) {
	// it is decay mu
	motherType = kPrimaryMu;
	break; // particle loop
      }
      else if ( step == 1 || isFirstMotherHF ){
	// Check if it is HF muon
	// avoid the case when the HF was not the first mother
	// (the mu is produced by a decain chain => decay mu)
	Bool_t foundHF = kFALSE;
	for ( Int_t idec=0; idec<3; idec++ ) {
	  for ( Int_t ihf=0; ihf<2; ihf++ ) {
	    if ( ( absPdg/den[idec] ) != baseFlv[ihf] ) continue;
	    motherType = mType[ihf];
	    foundHF = kTRUE;
	    break;
	  } // loop on hf
	  if ( foundHF ) {
	    if ( step == 1 ) isFirstMotherHF = kTRUE;
	    break; // break loop on pdg code
	    // but continue the global loop to find higher mass HF
	  }
	} // loop on pdg codes
	if ( ! foundHF ) {
	  motherType = kPrimaryMu;
	  break;
	}
	else if ( absPdg < 10 ) {
	  // found HF quark: break particle loop
	  break;
	}
      } // potential HF code
      */
      for ( Int_t idec=0; idec<3; idec++ ) {
	Int_t flv = absPdg/den[idec];
	if ( flv > 0 && flv < 4 ) return kPrimaryMu;
	else if ( flv == 0 || flv > 5 ) continue;
	else {
	  motherType = ( flv == 4 ) ? kCharmMu : kBeautyMu;
	  break; // break loop on pdg code
	  // but continue the global loop to find higher mass HF
	}
      } // loop on pdg code
      if ( absPdg < 10 ) break; // particle loop
    } // is primary
    else {
      // If hadronic process => secondary
      if ( part->GetUniqueID() == kPHadronic ) {
	//motherType = kSecondaryMu;
	//break; // particle loop
	return kSecondaryMu;
      }
    } // is secondary

    imother = part->GetFirstMother();

  } // loop on mothers

  return motherType;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSingleMu::GetHistoIndex(Int_t histoTypeIndex, Int_t trigIndex, Int_t srcIndex)
{
  //
  /// Get histogram index in the list
  //

  if ( srcIndex < 0 ) {
    return histoTypeIndex;
  }

  return
    kNsummaryHistosMC + 
    kNtrackSources * kNtrigCuts * histoTypeIndex  + 
    kNtrackSources * trigIndex  + 
    srcIndex;
}

//________________________________________________________________________
Float_t AliAnalysisTaskSingleMu::GetBinThetaAbsEnd(Float_t RAtAbsEnd, Bool_t isTheta)
{
  //
  /// Get bin of theta at absorber end region
  //
  Float_t thetaDeg = ( isTheta ) ? RAtAbsEnd : TMath::ATan( RAtAbsEnd / 505. );
  thetaDeg *= TMath::RadToDeg();
  if ( thetaDeg < 2. )
    return 0.;
  else if ( thetaDeg < 3. )
    return 1.;
  else if ( thetaDeg < 9. )
    return 2.;

  return 3.;
}


//________________________________________________________________________
void AliAnalysisTaskSingleMu::Reset(Bool_t keepGlobal)
{
  //
  /// Reset variables
  //
  Int_t lastVarFloat = ( keepGlobal ) ? kVarIPVx : kNvarFloat;
  for(Int_t ivar=0; ivar<lastVarFloat; ivar++){
    fVarFloat[ivar] = 0.;
  }
  Int_t lastVarInt = ( keepGlobal ) ? kVarPassPhysicsSelection : kNvarInt;
  for(Int_t ivar=0; ivar<lastVarInt; ivar++){
    fVarInt[ivar] = 0;
  }
  fVarInt[kVarMatchTrig] = -1;

  if ( ! keepGlobal ){
    for(Int_t ivar=0; ivar<kNvarChar; ivar++){
      strncpy(fVarChar[ivar]," ",255);
    }
    for(Int_t ivar=0; ivar<kNvarUInt; ivar++){
      fVarUInt[ivar] = 0;
    }
  }
  if ( fMCEvent ){
    lastVarFloat = ( keepGlobal ) ? kVarIPVxMC : kNvarFloatMC;
    for(Int_t ivar=0; ivar<lastVarFloat; ivar++){
      fVarFloatMC[ivar] = 0.;
    }
    for(Int_t ivar=0; ivar<kNvarIntMC; ivar++){
      fVarIntMC[ivar] = 0;
    }
    fVarIntMC[kVarMotherType] = -1;
  }
}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::SetTriggerClasses(TString triggerClasses)
{
  /// Set trigger classes
  if ( fTriggerClasses )
    delete fTriggerClasses;

  fTriggerClasses = triggerClasses.Tokenize(" ");
  fTriggerClasses->SetOwner();

}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::SetAxisLabel(TAxis* axis)
{
  //
  /// Set the labels of trigger class axis
  //
  for ( Int_t itrig=0; itrig<fTriggerClasses->GetEntries(); itrig++ ) {
    TString trigName = ((TObjString*)fTriggerClasses->At(itrig))->GetString();
    axis->SetBinLabel(itrig+1,trigName.Data());
  }
    axis->SetTitle("Fired class");
}
