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
          fHistManager(),fAnalysisType(kppData),R_jet_step(0.1),Max_match_dist(0.2)// , fOutput{0}
  {
  }
  /**
   * Standard constructor, intended for use by the the user.
   * @param[in] name Name of the task
   */
  AliAnalysisTaskEmcalJetEnergyFlow::AliAnalysisTaskEmcalJetEnergyFlow(const char* name):
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistManager(name),fAnalysisType(kppData),R_jet_step(0.1),Max_match_dist(0.2)// , fOutput{0}
  {
          SetMakeGeneralHistograms(1);
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
    //      AllocateClusterHistograms();
    //      AllocateCellHistograms();
          AllocateEnergyflowHistograms();
  
    TIter next(fHistManager.GetListOfHistograms());
    TObject* obj = 0;
    if(fAnalysisType==kppMC) fOutput->SetUseScaling(kTRUE);
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
      //    DoClusterLoop();
      //    DoCellLoop();
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
           Double_t   Rstep_EF,
           Double_t   Max_match_dr,
           Double_t   Lead_pt_cut,
           AnalysisType fAnType,
           const char *suffix )
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
    TString trackName2("mcparticles"); // If the task runs over MC, then this is used to access the generator level tracks
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
    
    if (fAnType==kppData)name+= "_ppdata";
    else if (fAnType==kPbPbData)name+= "_PbPbdata";
    else if (fAnType==kppMC)name+= "_ppMC";
    else if (fAnType==kEmbedded)name+= "_Embed";
  
    if (strcmp(suffix,"") != 0) {
      name += "_";
      name += suffix;
    }
  
    AliAnalysisTaskEmcalJetEnergyFlow* EFTask = new AliAnalysisTaskEmcalJetEnergyFlow(name);
    EFTask->R_jet_step=Rstep_EF;
    EFTask->SetCaloCellsName(cellName);
    EFTask->Max_match_dist=Max_match_dr;
    EFTask->LeadPtCut=Lead_pt_cut;
    EFTask->SetVzRange(-10,10);
    EFTask->SetAnalysisType(fAnType);
    if((fAnType==kppData)||(fAnType==kppMC))EFTask->SetForceBeamType(AliAnalysisTaskEmcal::kpp);  
    if((fAnType==kPbPbData)||(fAnType==kEmbedded))EFTask->SetForceBeamType(AliAnalysisTaskEmcal::kAA);
   // The first case would be used to run only on generator level of MC productions
    if (trackName == "mcparticles") {           
      EFTask->AddMCParticleContainer(trackName);
    }
    else if (trackName == "tracks" || trackName == "Tracks") {
      EFTask->AddTrackContainer(trackName);
    }
    else if (!trackName.IsNull()) {
      EFTask->AddParticleContainer(trackName);
    }
 if ((fAnType==kppMC)||(fAnType==kEmbedded))EFTask->AddMCParticleContainer(trackName2); //In the case of MC analysis, we add the additional particle container for the generator level analysis
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


void AliAnalysisTaskEmcalJetEnergyFlow::AllocateJetHistograms(){
    TString histname;
    TString histtitle;
    TString groupname;
    Int_t fNPtBins = 40; 
    Double_t fMinPtBin = 0.0; 
    Double_t fMaxPtBin = 200.0;
    Int_t fNEtaBins = 20; 
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
        fHistManager.CreateTH1(histname, histtitle, fNPtBins, 0, 3);
  
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
        fHistManager.CreateTH1(histname, histtitle, fNPtBins+10, -fMaxPtBin /4, fMaxPtBin);
        }
      }
    }
  
  }

void AliAnalysisTaskEmcalJetEnergyFlow::DoJetLoop(){
  
  TString histname;
    TString groupname;
    AliJetContainer* jetCont = 0;
    Double_t LeadTrack_cut = LeadPtCut;   //Pt cut on the leading trach of the jet
    TIter next(&fJetCollArray);
    while ((jetCont = static_cast<AliJetContainer*>(next()))) {
      groupname = jetCont->GetName();
      UInt_t count = 0;
      jetCont->SetJetEtaLimits(-.5,.5);
      for(auto jet : jetCont->accepted()) {
        if (!jet) continue;
        Double_t ptLeading = jetCont->GetLeadingHadronPt(jet);
        if(ptLeading<LeadTrack_cut) continue;
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
    Double_t Rstep = R_jet_step;
    Int_t fNPtBins = 40;
    Double_t fMinPtBin = 0.0;
    Double_t fMaxPtBin = 200.0;
    Int_t fNEtaBins = 20;
    Double_t fMaxEtaBin = 1.0;
    Int_t fNDPtBins = 100;
    Double_t fMaxDPtBin = 79.5;
    Double_t fMinDPtBin = -20.5;
    Int_t fNDRBins = 120; //Maybe this is too fine
    Double_t fMaxDRBin = 0.3;
    Double_t fMinDRBin = 0.0;
    Int_t Bins[4] = {fNPtBins,fNPtBins,fNDPtBins,fNDPtBins}; //These three arrays are needed for the THnSparse Response matrix
    Double_t MaxBin[4] = {fMaxPtBin,fMaxPtBin,fMaxDPtBin,fMaxDPtBin};
    Double_t MinBin[4] = {fMinPtBin,fMinPtBin,fMinDPtBin,fMinDPtBin};

  Int_t Pair_number = 3;
    if((fAnalysisType==kPbPbData)||(fAnalysisType==kEmbedded))fNcentBins=5;
    else fNcentBins=1;

    if ((fAnalysisType==kppMC)||(fAnalysisType==kEmbedded))Pair_number = (fJetCollArray.GetEntries()/2)-1;
    else Pair_number = fJetCollArray.GetEntries()-1;
for (Int_t cent = 0; cent < fNcentBins; cent++){

          histname = TString::Format("hBkgRho_%d",cent);
          histtitle = TString::Format("Background density median #rho distribution;#rho");
          fHistManager.CreateTH1(histname,histtitle,50,0,300);

        for(Int_t i=0;i<Pair_number;i++){
       
          histname = TString::Format("hJetPtDeltaPt_R%03d_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("#DeltaP_{t} between %.2f and %.2f jet radii;P_{t,jet(R=%.2f)}(GeV/c);#Delta P_{t}(GeV/c)",Rstep*(i+1),Rstep*(i+2),Rstep*(i+1));
          fHistManager.CreateTH2(histname,histtitle,fNPtBins,fMinPtBin,fMaxPtBin,fNDPtBins,fMinDPtBin,fMaxDPtBin);

          histname = TString::Format("hJetPtSmallDeltaPt_R%03d_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("Small #DeltaP_{t} between %.2f and %.2f jet radii;P_{t,jet(R=%.2f)}(GeV/c);#Delta P_{t}(GeV/c)",Rstep*(i+1),Rstep*(i+2),Rstep*(i+1));
          fHistManager.CreateTH2(histname,histtitle,fNPtBins,fMinPtBin,fMaxPtBin,2*fNDPtBins,0,2);
  
          histname = TString::Format("hJetPtDeltaRDeltaPt_R%03d_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("#DeltaP_{t} between %.2f and %.2f jet radii vs #DeltaR;P_{t,jet(R=%.2f)};#DeltaR;#DeltaP_{t}",Rstep*(i+1),Rstep*(i+2),Rstep*(i+1));
          fHistManager.CreateTH3(histname,histtitle,fNPtBins,fMinPtBin,fMaxPtBin,fNDRBins,fMinDRBin,fMaxDRBin,fNDPtBins,fMinDPtBin,fMaxDPtBin);
  
          histname = TString::Format("hJetPtDeltaR_R%03d_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("#DeltaR between %.2f and %.2f jet radii;P_{t,R=%.2f};#DeltaR",Rstep*(i+1),Rstep*(i+2),Rstep*(i+1));
          fHistManager.CreateTH2(histname,histtitle,fNPtBins,fMinPtBin,fMaxPtBin,fNDRBins,fMinDRBin,fMaxDRBin);
  
          histname = TString::Format("hDeltaPtvPtvMultiplicity_R%03d_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("#DeltaP_{t} between %.2f and %.2f jet radii vs P_{t} vs Multiplicity;#Delta P_{t};P_{t,jet(R=%f)}(GEV/v); Multiplicity",Rstep*(i+1),Rstep*(i+2),Rstep*(i+1));
          fHistManager.CreateTH3(histname,histtitle,fNDPtBins,fMinDPtBin,fMaxDPtBin,fNPtBins,fMinPtBin,fMaxPtBin,20,0,100);
         
          histname = TString::Format("hMatchedJetPt_R%03d_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("Matched jet P_{t} spectrum of R=%.2f",Rstep*(i+1));
          fHistManager.CreateTH1(histname, histtitle,fNPtBins,fMinPtBin,fMaxPtBin);
  
          histname = TString::Format("hMatchedJetEta_R%03d_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("Matched jet #eta spectrum of R=%.2f",Rstep*(i+1));
          fHistManager.CreateTH1(histname, histtitle,fNEtaBins,-fMaxEtaBin,fMaxEtaBin);

          if(i==Pair_number-1){                                         //Once we arrive at the last pair iteration, make the spectra for the large R jet of the pair
          histname = TString::Format("hMatchedJetPt_R%03d_%d",int(Rstep*(i+2)*100),cent);
          histtitle = TString::Format("Matched jet P_{t} spectrum of R=%.2f",Rstep*(i+2));
          fHistManager.CreateTH1(histname, histtitle,fNPtBins,fMinPtBin,fMaxPtBin);
  
          histname = TString::Format("hMatchedJetEta_R%03d_%d",int(Rstep*(i+2)*100),cent);
          histtitle = TString::Format("Matched jet #eta spectrum of R=%.2f",Rstep*(i+2));
          fHistManager.CreateTH1(histname, histtitle,fNEtaBins,-fMaxEtaBin,fMaxEtaBin);
                          }
         }
 
        if((fAnalysisType==kppMC)||(fAnalysisType==kEmbedded)){
        for(Int_t i=0;i<Pair_number;i++){

          histname = TString::Format("hJetPtMismatches_R%03d_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("Ammount of mismatches between the two levels for R_{jet}= %.2f;P_{t,jet(R=%.2f)};N_{mismatches}",Rstep*(i+2),Rstep*(i+1));
          fHistManager.CreateTH1(histname,histtitle,fNPtBins,fMinPtBin,fMaxPtBin);

          histname = TString::Format("hJetPtDeltaPt_R%03d_gen_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("#DeltaP_{t} between %.2f and %.2f jet radii (Generator level);P_{t,jet(R=%.2f)}(GeV/c);#Delta P_{t}(GeV/c)",Rstep*(i+1),Rstep*(i+2),Rstep*(i+1));
          fHistManager.CreateTH2(histname,histtitle,fNPtBins,fMinPtBin,fMaxPtBin,fNDPtBins,fMinDPtBin,fMaxDPtBin);

          histname = TString::Format("hJetPtSmallDeltaPt_R%03d_gen_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("Small #DeltaP_{t} between %.2f and %.2f jet radii (Generator level);P_{t,jet(R=%.2f)}(GeV/c);#Delta P_{t}(GeV/c)",Rstep*(i+1),Rstep*(i+2),Rstep*(i+1));
          fHistManager.CreateTH2(histname,histtitle,fNPtBins,fMinPtBin,fMaxPtBin,2*fNDPtBins,0,2); 
 
          histname = TString::Format("hJetPtDeltaRDeltaPt_R%03d_gen_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("#DeltaP_{t} between %.2f and %.2f jet radii vs #DeltaR;P_{t,jet(R=%.2f)} (Generator level);#DeltaR;#DeltaP_{t}",Rstep*(i+1),Rstep*(i+2),Rstep*(i+1));
          fHistManager.CreateTH3(histname,histtitle,fNPtBins,fMinPtBin,fMaxPtBin,fNDRBins,fMinDRBin,fMaxDRBin,fNDPtBins,fMinDPtBin,fMaxDPtBin);
       
          histname = TString::Format("hJetPtDeltaR_R%03d_gen_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("#DeltaR between %.2f and %.2f jet radii (Generator level);P_{t,R=%.2f};#DeltaR",  Rstep*(i+1),Rstep*(i+2),Rstep*(i+1));
          fHistManager.CreateTH2(histname,histtitle,fNPtBins,fMinPtBin,fMaxPtBin,fNDRBins,fMinDRBin,fMaxDRBin);
 
          histname = TString::Format("hDeltaPtvPtvMultiplicity_R%03d_gen_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("#DeltaP_{t} between %.2f and %.2f jet radii vs P_{t} vs Multiplicity (Generator level);#Delta P_{t};P_{t,jet(R=%f)}(GEV/v); Multiplicity",Rstep*(i+1),Rstep*(i+2),Rstep*(i+1));
          fHistManager.CreateTH3(histname,histtitle,fNDPtBins,fMinDPtBin,fMaxDPtBin,fNPtBins,fMinPtBin,fMaxPtBin,20,0,100);
  
          histname = TString::Format("hMatchedJetPt_R%03d_gen_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("Matched jet P_{t} spectrum of R=%.2f (Generator level)",Rstep*(i+1));
          fHistManager.CreateTH1(histname, histtitle,fNPtBins,fMinPtBin,fMaxPtBin);
  
          histname = TString::Format("hMatchedJetEta_R%03d_gen_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("Matched jet #eta spectrum of R=%.2f (Generator level)",Rstep*(i+1));
          fHistManager.CreateTH1(histname, histtitle,fNEtaBins,-fMaxEtaBin,fMaxEtaBin);
  
          histname = TString::Format("hJetEnergyResponse_R%03d_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("Relative Jet Energy Scale (Hybrid-Truth)/Truth  of R=%.2f",Rstep*(i+1));
          fHistManager.CreateTH2(histname, histtitle,fNPtBins,fMinPtBin,fMaxPtBin,100,-1,1);

          histname = TString::Format("hMatchedJetEnergyResponse_R%03d_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("Matched Relative Jet Energy Scale (Hybrid-Truth)/Truth  of R=%.2f",Rstep*(i+1));
          fHistManager.CreateTH2(histname, histtitle,fNPtBins,fMinPtBin,fMaxPtBin,100,-1,1);

          histname = TString::Format("ResponseMatrix_R%03d_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("Response Matrix of R %.2f;P_{t} (Generator level); P_{t} (Detector level);#DeltaP_{t} (Generator level); #DeltaP_{t} (Detector level)",Rstep*(i+1));  
          fHistManager.CreateTHnSparse(histname,histtitle,4,Bins,MinBin,MaxBin);

          histname = TString::Format("MismatchResponseMatrix_R%03d_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("Response Matrix of mismatches for R %.2f;P_{t} (Generator level); P_{t} (Detector level);#DeltaP_{t} (Generator level); #DeltaP_{t} (Detector level)",Rstep*(i+1));  
          fHistManager.CreateTHnSparse(histname,histtitle,4,Bins,MinBin,MaxBin);


          histname = TString::Format("DeltaResponseMatrix_R%03d_%d",int(Rstep*(i+1)*100),cent);
          histtitle = TString::Format("DeltaResponseMatrix_R%.2f;P_{t,low} (GeV/c);P_{t,det}-P_{t,gen} (GeV/c);#DeltaP_{t,det}-#DeltaP_{t,gen} (GeV/c)",Rstep*(i+1));
          fHistManager.CreateTH3(histname, histtitle,fNPtBins,fMinPtBin,fMaxPtBin,100,-100,100,2*fNDPtBins,-fMaxDPtBin,fMaxDPtBin);

          if(i==Pair_number-1){
          histname = TString::Format("hMatchedJetPt_R%03d_gen_%d",int(Rstep*(i+2)*100),cent);
          histtitle = TString::Format("Matched jet P_{t} spectrum of R=%.2f (Generator level)",Rstep*(i+2));
          fHistManager.CreateTH1(histname, histtitle,fNPtBins,fMinPtBin,fMaxPtBin);
  
          histname = TString::Format("hMatchedJetEta_R%03d_gen_%d",int(Rstep*(i+2)*100),cent);
          histtitle = TString::Format("Matched jet #eta spectrum of R=%.2f (Generator level)",Rstep*(i+2));
          fHistManager.CreateTH1(histname, histtitle,fNEtaBins,-fMaxEtaBin,fMaxEtaBin);
                          }
                } // End of pairs' loop
        }       //End of if MC clause
        }       //End of centrality loop
} 



void AliAnalysisTaskEmcalJetEnergyFlow::FillEFHistograms(){

  Int_t Njets = 200;      //Just a large number so that the matching matrix created by the matcher task will always have the size of the input jet lists
  Int_t &kLowRJets = Njets;
  Int_t &kHighRJets = Njets;
  
  TList DetLowRJetsList;     // List of the accepted (lower R) jets (Detector level or data)
  TList DetHighRJetsList;    // List of the accepted (higher R) jets (Detector level or data)

  TList MatchGenDetList;     //List of accepted(small R) jets at the generator level that match to detector level (Only used for MC runs)
  TList GenHighRJetsList;    // List of the accepted (higher R) jets (Generator level)  (Only used for MC)

  //This array points to the low R jet that matches to each high R jet
  TArrayI iLowRIndex_det;       //Detector level (MC or data)
  TArrayI iLowRIndex_gen;       //Generator level (Only MC)
  
  // This array points to the high R jet that matches to each low R jet
  TArrayI iHighRIndex_det;      //Detector level (MC or data)
  TArrayI iHighRIndex_gen;      //Generator level (Only MC)
 
  TString histname;
  TString Contname;
  Double_t Rho = 0.0;             //Median background density of the event
  Double_t Rstep = R_jet_step;            //This is the size of the radius step used in the jet expansion
  Double_t DeltaPt_gen = 0.0;           //DeltaPt @ Gen-level
  Double_t DeltaPt_det = 0.0;           //DeltaPt @ Det-level
  Double_t pt_Ldet = 0.0;               //Low R jet pt @ Det level (needed for the R-matrix)
  Double_t pt_Hdet = 0.0;               //High R jet pt @Det level (needed for the R-matrix)
  Double_t DeltaEta = 0.0;
 // Double_t Mismatch_count =0.0;         //Counter for mismatches between truth and detector level
 // Double_t Response_count =0.0;         //Counter for all truth level jets which made it to the response matrix
  Double_t LeadTrack_cut = LeadPtCut;   //Pt cut on the leading trach of the jet
  Float_t Max_dist = Max_match_dist;    //Maximum distance used as a matching criterion for JetMatcher
  Float_t max_eta = 0.5;                //Maximum eta jet acceptance (cross-check within the JetMatcher)
  Double_t ResponseData[4];             //Array needed for filling the THnSparse R-matrix
  Bool_t   Bkg_sub_method =kFALSE;      //Boolean flag in case a constituent subtraction method is chosen instead of rho*A
  AliJetContainer* DetjetCont1=0; // One container for the jets of the lower R (of each comparison pair)
  AliJetContainer* DetjetCont2=0; // One container for the jets of the higher R (of each comparison pair)
  AliJetContainer* GenjetCont1=0; // One container for the jets of the lower R (of each comparison pair)
  AliJetContainer* GenjetCont2=0; // One container for the jets of the higher R (of each comparison pair)

   Int_t Pair_number = 3;
  AliEmcalJet* Jet_genlowR =0; //Convenient pointer for the EF calculations @ Gen-level
  Int_t NumJet= 4;
if ((fAnalysisType==kppMC)||(fAnalysisType==kEmbedded)) NumJet=fJetCollArray.GetEntries()/2;
else NumJet =fJetCollArray.GetEntries();
Pair_number= NumJet-1;

//Loop over the number of comparison pairs
  for (Int_t i=0;i<Pair_number;i++)
  {       DetLowRJetsList.Clear();
          DetHighRJetsList.Clear();
          DeltaPt_gen = 0.0;
          DeltaPt_det = 0.0;
          pt_Ldet = 0.0;
          pt_Hdet = 0.0;
          DeltaEta = 0.0;
        DetjetCont1=0; DetjetCont2=0; GenjetCont1=0; GenjetCont2=0; //Reseting the containers at the start of each pair loop     
        if((fAnalysisType==kppMC)||(fAnalysisType==kEmbedded)){
                MatchGenDetList.Clear();
                GenHighRJetsList.Clear();

           Contname= dynamic_cast<AliJetContainer*>(fJetCollArray[i])->GetName();
                if(Contname.Contains("TruthJet")||Contname.Contains("GenJet"))
                {GenjetCont1 = dynamic_cast<AliJetContainer*>(fJetCollArray[i]);
                 GenjetCont1->SetJetEtaLimits(-max_eta,max_eta); 
                }
                else if(Contname.Contains("DetJet")&&(fAnalysisType==kEmbedded))
                {GenjetCont1 = dynamic_cast<AliJetContainer*>(fJetCollArray[i]);
                 GenjetCont1->SetJetEtaLimits(-max_eta,max_eta);
                }
                else{ 
                DetjetCont1 = dynamic_cast<AliJetContainer*>(fJetCollArray[i]);
                DetjetCont1->SetJetEtaLimits(-max_eta,max_eta); 
                }
                if (Contname.Contains("ConstSub"))Bkg_sub_method=kTRUE;

         Contname= dynamic_cast<AliJetContainer*>(fJetCollArray[i+1])->GetName();
                  if(Contname.Contains("TruthJet")||Contname.Contains("GenJet"))
                  {GenjetCont2 = dynamic_cast<AliJetContainer*>(fJetCollArray[i+1]);
                   GenjetCont2->SetJetEtaLimits(-max_eta,max_eta);
                  }
                  else if(Contname.Contains("DetJet")&&(fAnalysisType==kEmbedded))
                  {GenjetCont2 = dynamic_cast<AliJetContainer*>(fJetCollArray[i+1]);
                   GenjetCont2->SetJetEtaLimits(-max_eta,max_eta);
                  }
                  else{ 
                  DetjetCont2 = dynamic_cast<AliJetContainer*>(fJetCollArray[i+1]);
                  DetjetCont2->SetJetEtaLimits(-max_eta,max_eta);
                  }
                if (Contname.Contains("ConstSub"))Bkg_sub_method=kTRUE;
      
         Contname= dynamic_cast<AliJetContainer*>(fJetCollArray[i+NumJet])->GetName();
                  if(Contname.Contains("TruthJet")||Contname.Contains("GenJet"))
                  {GenjetCont1 = dynamic_cast<AliJetContainer*>(fJetCollArray[i+NumJet]);
                   GenjetCont1->SetJetEtaLimits(-max_eta,max_eta);
                  }
                  else if(Contname.Contains("DetJet")&&(fAnalysisType==kEmbedded))
                  {GenjetCont1 = dynamic_cast<AliJetContainer*>(fJetCollArray[i+NumJet]);
                   GenjetCont1->SetJetEtaLimits(-max_eta,max_eta);
                  }
                  else{ 
                  DetjetCont1 = dynamic_cast<AliJetContainer*>(fJetCollArray[i+NumJet]);
                  DetjetCont1->SetJetEtaLimits(-max_eta,max_eta);
                  }
                if (Contname.Contains("ConstSub"))Bkg_sub_method=kTRUE;

         Contname= dynamic_cast<AliJetContainer*>(fJetCollArray[i+NumJet+1])->GetName();
                  if(Contname.Contains("TruthJet")||Contname.Contains("GenJet"))
                  {GenjetCont2 = dynamic_cast<AliJetContainer*>(fJetCollArray[i+NumJet+1]);
                   GenjetCont2->SetJetEtaLimits(-max_eta,max_eta);
                  }
                  else if(Contname.Contains("DetJet")&&(fAnalysisType==kEmbedded))
                  {GenjetCont2 = dynamic_cast<AliJetContainer*>(fJetCollArray[i+NumJet+1]);
                   GenjetCont2->SetJetEtaLimits(-max_eta,max_eta);
                  }
                  else{ 
                  DetjetCont2 = dynamic_cast<AliJetContainer*>(fJetCollArray[i+NumJet+1]);
                  DetjetCont2->SetJetEtaLimits(-max_eta,max_eta);
                  }
                if (Contname.Contains("ConstSub"))Bkg_sub_method=kTRUE;

                //------------Debugging tool/ check that containers are loaded correctly i.e.All Gen-level followed by all Det-level or vice versa. TLDR no container mixing
            //    Printf(TString(GenjetCont1->GetName()));   Printf(TString(GenjetCont2->GetName()));
           //     Printf(TString(DetjetCont1->GetName()));   Printf(TString(DetjetCont2->GetName()));
                         
                for(auto jet:DetjetCont1->accepted())if(DetjetCont1->GetLeadingHadronPt(jet)>=LeadTrack_cut)DetLowRJetsList.Add(jet);
                for(auto jet:DetjetCont2->accepted())if(DetjetCont2->GetLeadingHadronPt(jet)>=LeadTrack_cut)DetHighRJetsList.Add(jet);
                for(auto jet:GenjetCont1->accepted())if(GenjetCont1->GetLeadingHadronPt(jet)>=LeadTrack_cut)if(jet->ClosestJet())MatchGenDetList.Add(jet);
                for(auto jet:GenjetCont2->accepted())if(GenjetCont2->GetLeadingHadronPt(jet)>=LeadTrack_cut) GenHighRJetsList.Add(jet);
                // First, match and calculate the energy flow on the detector level jets for this R pair iteration
         if(DetLowRJetsList.GetEntries()==0||DetHighRJetsList.GetEntries()==0) continue;
         //     Printf("Point A \n");
            iLowRIndex_det.Set(DetHighRJetsList.GetEntries()); iHighRIndex_det.Set(DetLowRJetsList.GetEntries());
            JetMatcher(&DetLowRJetsList,kLowRJets,&DetHighRJetsList,kHighRJets, iLowRIndex_det,iHighRIndex_det,0,Max_dist,max_eta);
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                for (Int_t j=0; j<iHighRIndex_det.GetSize();j++)      // loop over the low R jets (Detector level)
                    {
                        //Calculating the jet response between 2 levels e.g. truth-detector
           //             Printf(Form("Jet: %d",j));
                        if((dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->ClosestJet()){    
                        Double_t pt_det_JER=(dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->Pt();
                        Double_t pt_tru_JER=((dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->ClosestJet())->Pt();                      
                        if(!Bkg_sub_method){
                        if(DetjetCont1->GetRhoParameter())pt_det_JER = pt_det_JER - DetjetCont1->GetRhoVal() * (dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->Area();
                        if(GenjetCont1->GetRhoParameter())pt_tru_JER = pt_tru_JER- DetjetCont1->GetRhoVal()*((dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->ClosestJet())->Area();
                                                }
                        histname = TString::Format("hJetEnergyResponse_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                        fHistManager.FillTH2(histname,pt_tru_JER,(pt_det_JER-pt_tru_JER));                      
                                                                                                }
                        //Check the match conditions between Ri and Ri+1 at the detector level
                            if(iHighRIndex_det[j]>=0){                  // if there is a match
                            Int_t match_index = iHighRIndex_det[j];
                            if (iLowRIndex_det[match_index]==j){        //And if the match is bijective
  
                    Double_t pt_low = (dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->Pt();
                    Double_t pt_high = (dynamic_cast<AliEmcalJet*>(DetHighRJetsList.At(match_index)))->Pt();

                        //In case of PbPb or Embedded analysis correct the jet pt by rho*area (This may have to be removed if we want to look at const subtraction) 
                    
                if(!Bkg_sub_method){
                    if(DetjetCont1->GetRhoParameter())pt_low = (dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->Pt() - DetjetCont1->GetRhoVal() * (dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->Area();
                    if(DetjetCont2->GetRhoParameter())pt_high = (dynamic_cast<AliEmcalJet*>(DetHighRJetsList.At(match_index)))->Pt() - DetjetCont2->GetRhoVal() *(dynamic_cast<AliEmcalJet*>(DetHighRJetsList.At(match_index)))->Area();}

                    if((pt_low<=0)||(pt_high<=0)) continue;
                    Double_t eta_low = (dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->Eta();
                    Double_t eta_high = (dynamic_cast<AliEmcalJet*>(DetHighRJetsList.At(match_index)))->Eta();
                    if((pt_low>0)&&(pt_high>0)) DeltaPt_det = pt_high-pt_low;
                    Double_t DeltaR = (dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->DeltaR((dynamic_cast<AliEmcalJet*>(DetHighRJetsList.At(match_index))));
                    Double_t DeltaEta = fabs(eta_high - eta_low);
             //           Printf("Point B \n");
                        if((pt_low>0)&&(pt_high>0)){
                    histname = TString::Format("hJetPtDeltaPt_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                    fHistManager.FillTH2(histname,pt_low,DeltaPt_det);

                        if((DeltaPt_det<=2)&&(DeltaPt_det>=0)){
                   histname = TString::Format("hJetPtSmallDeltaPt_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                   fHistManager.FillTH2(histname,pt_low,DeltaPt_det);}
  
//                    histname = TString::Format("hJetPtConstZ_R%02d_%d",int(Rstep*(i+1)*100),fCentBin);
//                    for(auto cont:dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j))->GetParticleConstituents()){
//                     fHistManager.FillTH2(histname,pt_low,cont.Pt()/pt_low);}

                    histname = TString::Format("hJetPtDeltaR_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                    fHistManager.FillTH2(histname,pt_low,DeltaR);

                     histname = TString::Format("hJetPtDeltaRDeltaPt_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                     fHistManager.FillTH3(histname,pt_low,DeltaR,DeltaPt_det);

                    histname = TString::Format("hMatchedJetPt_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                    fHistManager.FillTH1(histname,pt_low);
  
                    histname = TString::Format("hMatchedJetEta_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                    fHistManager.FillTH1(histname,eta_low);
                    
                    histname = TString::Format("hDeltaPtvPtvMultiplicity_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                    fHistManager.FillTH3(histname,DeltaPt_det,pt_low,(dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->GetNumberOfConstituents());
                       
                    if(i== Pair_number-1){
  
                    histname = TString::Format("hMatchedJetPt_R%03d_%d",int(Rstep*(i+2)*100),fCentBin);
                    fHistManager.FillTH1(histname,pt_high);
  
                    histname = TString::Format("hMatchedJetEta_R%03d_%d",int(Rstep*(i+2)*100),fCentBin);
                    fHistManager.FillTH1(histname,eta_high);
                                                         }}
               //      Printf("Point C \n");                                    
                                                             } //And if the match is bijective
                                                  } // if there is a match
              //          Printf("Point D \n");
                  } // loop over the low R jets (Detector level)

                // Repeat the energy flow calculation on the generator level jets for this R pair iteration (BUT ONLY FOR THE LOW R JETS THAT HAVE A MATCH ON DETECTOR LEVEL)
            //            Printf("Point A1 \n");
           if(MatchGenDetList.GetEntries()==0||GenHighRJetsList.GetEntries()==0) continue;
              iLowRIndex_gen.Set(GenHighRJetsList.GetEntries()); iHighRIndex_gen.Set(MatchGenDetList.GetEntries());
              JetMatcher(&MatchGenDetList,kLowRJets,&GenHighRJetsList,kHighRJets, iLowRIndex_gen,iHighRIndex_gen,0,Max_dist,max_eta);
   //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                  for (Int_t j=0; j<iHighRIndex_gen.GetSize();j++) // loop over the low R jets (Generator level that match to Detector level)
                      {
                              if(iHighRIndex_gen[j]>=0){        // if there is a match
                              Int_t match_index = iHighRIndex_gen[j];
                              if (iLowRIndex_gen[match_index]==j){      //And if the match is bijective
                        Jet_genlowR = dynamic_cast<AliEmcalJet*>(MatchGenDetList.At(j)); //Useful definition for the following calculations
                      Double_t pt_low = Jet_genlowR->Pt();
                      Double_t pt_high = (dynamic_cast<AliEmcalJet*>(GenHighRJetsList.At(match_index)))->Pt();
                        if(!Bkg_sub_method){
                      if(GenjetCont1->GetRhoParameter())pt_low = Jet_genlowR->Pt() - GenjetCont1->GetRhoVal() * Jet_genlowR->Area();
                      if(GenjetCont2->GetRhoParameter())pt_high = (dynamic_cast<AliEmcalJet*>(GenHighRJetsList.At(match_index)))->Pt()- GenjetCont2->GetRhoVal()* (dynamic_cast<AliEmcalJet*>(GenHighRJetsList.At(match_index)))->Area();}
                       if((pt_low>0)&&(pt_high>0)){ 
                      Double_t eta_low = Jet_genlowR->Eta();
                      Double_t eta_high = (dynamic_cast<AliEmcalJet*>(GenHighRJetsList.At(match_index)))->Eta();
                      DeltaPt_gen = pt_high-pt_low;
                      Double_t DeltaR = Jet_genlowR->DeltaR((dynamic_cast<AliEmcalJet*>(GenHighRJetsList.At(match_index))));
                      Double_t DeltaEta = fabs(eta_high - eta_low);
                        DeltaPt_det = 0.0;
              //          Printf("Point A2 \n");
                      //------- In order to construct the Rmatrix we need to evaluate the pt & Dpt @ Det level while analysing the corresponding jets @ Gen level
                        pt_Ldet = Jet_genlowR->ClosestJet()->Pt();
                        if(DetjetCont1->GetRhoParameter())pt_Ldet = pt_Ldet - DetjetCont1->GetRhoVal() * Jet_genlowR->ClosestJet()->Area(); 
                        for (Int_t w =0;w<DetLowRJetsList.GetEntries();w++){
                                Int_t MI = iHighRIndex_det[w];
                        if(MI>=0){
                        if(((DetLowRJetsList.At(w) == Jet_genlowR->ClosestJet()))&&(iLowRIndex_det[MI]==w)) {
                                 pt_Hdet = (dynamic_cast<AliEmcalJet*>(DetHighRJetsList.At(MI)))->Pt();
                                
                                if((pt_low<0)||(pt_high<0))continue;
                        if(DetjetCont2->GetRhoParameter())pt_Hdet = pt_Hdet-DetjetCont2->GetRhoVal() * (dynamic_cast<AliEmcalJet*>(DetHighRJetsList.At(MI)))->Area();
                                DeltaPt_det = pt_Hdet -pt_Ldet;
                                 histname = TString::Format("ResponseMatrix_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                                 ResponseData[0]=pt_low;ResponseData[1]=pt_Ldet;ResponseData[2]=DeltaPt_gen;ResponseData[3]=DeltaPt_det;      
                                 fHistManager.FillTHnSparse(histname,ResponseData);

                                 histname = TString::Format("DeltaResponseMatrix_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                                fHistManager.FillTH3(histname,pt_low,pt_Ldet-pt_low,DeltaPt_det-DeltaPt_gen);

                        //Treatment of mismatches
                                histname = TString::Format("hJetPtMismatches_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                                if((!(dynamic_cast<AliEmcalJet*>(GenHighRJetsList.At(match_index)))->ClosestJet())||(dynamic_cast<AliEmcalJet*>(DetHighRJetsList.At(MI))->DeltaR((dynamic_cast<AliEmcalJet*>(GenHighRJetsList.At(match_index)))->ClosestJet())>1e-3)){
                                fHistManager.FillTH1(histname,pt_low);
                                histname = TString::Format("MismatchResponseMatrix_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                                ResponseData[0]=pt_low;ResponseData[1]=pt_Ldet;ResponseData[2]=DeltaPt_gen;ResponseData[3]=DeltaPt_det;
                                 fHistManager.FillTHnSparse(histname,ResponseData);
                                }
                                break;}}
                                          }
                      histname = TString::Format("hJetPtDeltaPt_R%03d_gen_%d",int(Rstep*(i+1)*100),fCentBin);
                      fHistManager.FillTH2(histname,pt_low,DeltaPt_gen);
                //        Printf("Point B1 \n");
                          if((DeltaPt_gen<=2)&&(DeltaPt_gen>=0)){
                     histname = TString::Format("hJetPtSmallDeltaPt_R%03d_gen_%d",int(Rstep*(i+1)*100),fCentBin);
                     fHistManager.FillTH2(histname,pt_low,DeltaPt_gen);}

                      histname = TString::Format("hMatchedJetEnergyResponse_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                        fHistManager.FillTH2(histname,pt_low,(pt_Ldet-pt_low)/pt_low);
  
                      histname = TString::Format("hJetPtDeltaRDeltaPt_R%03d_gen_%d",int(Rstep*(i+1)*100),fCentBin);
                      fHistManager.FillTH3(histname,pt_low,DeltaR,DeltaPt_gen);


                      histname = TString::Format("hJetPtDeltaR_R%03d_gen_%d",int(Rstep*(i+1)*100),fCentBin);
                      fHistManager.FillTH2(histname,pt_low,DeltaR);
  
                      histname = TString::Format("hMatchedJetPt_R%03d_gen_%d",int(Rstep*(i+1)*100),fCentBin);
                      fHistManager.FillTH1(histname,pt_low);
  
                      histname = TString::Format("hMatchedJetEta_R%03d_gen_%d",int(Rstep*(i+1)*100),fCentBin);
                      fHistManager.FillTH1(histname,eta_low);
  
                      histname = TString::Format("hDeltaPtvPtvMultiplicity_R%03d_gen_%d",int(Rstep*(i+1)*100),fCentBin);
                      fHistManager.FillTH3(histname,DeltaPt_gen,pt_low,Jet_genlowR->GetNumberOfConstituents());
  

                        if(i== Pair_number-1){
  
                      histname = TString::Format("hMatchedJetPt_R%03d_gen_%d",int(Rstep*(i+2)*100),fCentBin);
                      fHistManager.FillTH1(histname,pt_high);
  
                      histname = TString::Format("hMatchedJetEta_R%03d_gen_%d",int(Rstep*(i+2)*100),fCentBin);
                      fHistManager.FillTH1(histname,eta_high);
                                                           }}
                  //      Printf("Point C1 \n");                
                                                } //And if the match is bijective
                                                    } // if there is a match
                    } // loop over the low R jets (Generator level that match to detector level)
                 } //End of MC production case
        else{
          // Casting the lower R half of the pair to a container and getting the accepted jets to a list
          DetjetCont1 = dynamic_cast<AliJetContainer*>(fJetCollArray[i]);
          DetjetCont1->SetJetEtaLimits(-max_eta,max_eta);
          for(auto jet:DetjetCont1->accepted())if(DetjetCont1->GetLeadingHadronPt(jet)>=LeadTrack_cut) DetLowRJetsList.Add(jet);
          Contname= DetjetCont1->GetName();
          if (Contname.Contains("ConstSub"))Bkg_sub_method=kTRUE;
          // Casting the higher R half of the pair to a container and getting the accepted jets to a list
          DetjetCont2 = dynamic_cast<AliJetContainer*>(fJetCollArray[i+1]);
          DetjetCont2->SetJetEtaLimits(-max_eta,max_eta);
          for(auto jet:DetjetCont2->accepted())if(DetjetCont2->GetLeadingHadronPt(jet)>=LeadTrack_cut) DetHighRJetsList.Add(jet);
          Contname= DetjetCont2->GetName();
          if (Contname.Contains("ConstSub"))Bkg_sub_method=kTRUE;
        // For the case of data, match and calculate the energy flow only on the detector level jets for this R pair iteration
          if(DetLowRJetsList.GetEntries()==0||DetHighRJetsList.GetEntries()==0) continue;
          iLowRIndex_det.Set(DetHighRJetsList.GetEntries()); iHighRIndex_det.Set(DetLowRJetsList.GetEntries());
          JetMatcher(&DetLowRJetsList,kLowRJets,&DetHighRJetsList,kHighRJets, iLowRIndex_det,iHighRIndex_det,0,Max_dist,max_eta);
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

                for (Int_t j=0; j<iHighRIndex_det.GetSize();j++) // loop over the low R jets
                  {
                          if(iHighRIndex_det[j]>=0){  // if there is a match
                          Int_t match_index = iHighRIndex_det[j];
                          if (iLowRIndex_det[match_index]==j){  //And if the match is bijective
  
                  Double_t pt_low = (dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->Pt();
                  Double_t pt_high = (dynamic_cast<AliEmcalJet*>(DetHighRJetsList.At(match_index)))->Pt();
                if(!Bkg_sub_method){
                  if(DetjetCont1->GetRhoParameter())pt_low = (dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->Pt() - DetjetCont1->GetRhoVal() * (dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->Area();
                  if(DetjetCont2->GetRhoParameter())pt_high = (dynamic_cast<AliEmcalJet*>(DetHighRJetsList.At(match_index)))->Pt() - DetjetCont2->GetRhoVal() *(dynamic_cast<AliEmcalJet*>(DetHighRJetsList.At(match_index)))->Area();}

                  Double_t eta_low = (dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->Eta();
                  Double_t eta_high = (dynamic_cast<AliEmcalJet*>(DetHighRJetsList.At(match_index)))->Eta();
                  if((pt_low<=0)||(pt_high<=0))continue;
                  DeltaPt_det = pt_high-pt_low;
                  Double_t DeltaR = (dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->DeltaR((dynamic_cast<AliEmcalJet*>(DetHighRJetsList.At(match_index))));
                  Double_t DeltaEta = fabs(eta_high - eta_low);
                
                if((pt_low>0)&&(pt_high>0)){
                  histname = TString::Format("hJetPtDeltaPt_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);             
                  fHistManager.FillTH2(histname,pt_low,DeltaPt_det);
       
                if((DeltaPt_det<=2)&&(DeltaPt_det>=0)){
                   histname = TString::Format("hJetPtSmallDeltaPt_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                   fHistManager.FillTH2(histname,pt_low,DeltaPt_det);}
  
                  histname = TString::Format("hJetPtDeltaR_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                  fHistManager.FillTH2(histname,pt_low,DeltaR);

                  histname = TString::Format("hJetPtDeltaRDeltaPt_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                  fHistManager.FillTH3(histname,pt_low,DeltaR,DeltaPt_det);
                  
                  histname = TString::Format("hMatchedJetPt_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                  fHistManager.FillTH1(histname,pt_low);
        
                  histname = TString::Format("hMatchedJetEta_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                  fHistManager.FillTH1(histname,eta_low);

                  histname = TString::Format("hDeltaPtvPtvMultiplicity_R%03d_%d",int(Rstep*(i+1)*100),fCentBin);
                  fHistManager.FillTH3(histname,DeltaPt_det,pt_low,(dynamic_cast<AliEmcalJet*>(DetLowRJetsList.At(j)))->GetNumberOfConstituents());

                          if(i== Pair_number-1){
                
                  histname = TString::Format("hMatchedJetPt_R%03d_%d",int(Rstep*(i+2)*100),fCentBin);          
                  fHistManager.FillTH1(histname,pt_high);
                
                  histname = TString::Format("hMatchedJetEta_R%03d_%d",int(Rstep*(i+2)*100),fCentBin);
                  fHistManager.FillTH1(histname,eta_high);
                                               }}
                                                            } //And if the match is bijective
                                                } // if there is a match
                } // loop over the low R jets
                }// End of data case
} //End of loop over the R pair

            if (DetjetCont1->GetRhoParameter()){
                        Rho = DetjetCont1->GetRhoVal();
                        histname = TString::Format("hBkgRho_%d",fCentBin);
                        fHistManager.FillTH1(histname,Rho);
                                                }
}

void AliAnalysisTaskEmcalJetEnergyFlow::JetMatcher(
                                        const TList *genJetsList,const Int_t &kGenJets, 
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
         //     if (TClass(partCont->GetClassName()).InheritsFrom("AliVTrack")) {
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
 //       }
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
  
      //  if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
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
       // }
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

