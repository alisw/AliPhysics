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
 **************************************************************************/

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

#include "AliAnalysisTaskEmcalJetDijetMass.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetDijetMass);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalJetDijetMass::AliAnalysisTaskEmcalJetDijetMass() : 
    AliAnalysisTaskEmcalJet(),
    fHistManager(),
    fana(NULL),
    //fhistos(NULL),
    //fEmcalHistos(NULL),
    hTest(NULL),
    fJOutput(NULL)
{
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalJetDijetMass::AliAnalysisTaskEmcalJetDijetMass(const char *name) : 
    AliAnalysisTaskEmcalJet(name, kTRUE),
    fHistManager(name),
    fana(NULL),
    //fhistos(NULL),
    //fEmcalHistos(NULL),
    hTest(NULL),
    fJOutput(NULL)
{
    SetMakeGeneralHistograms(kTRUE);
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalJetDijetMass::~AliAnalysisTaskEmcalJetDijetMass()
{
    //delete fhistos;
    //delete fEmcalHistos;
    delete hTest;
    delete fana;
    delete fJOutput;
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalJetDijetMass::UserCreateOutputObjects()
{
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

    // For testing, later connect to addtask macro:
    fjetCone=0.4;
    fktJetCone=0.4;
    fktScheme=1;
    fantiktScheme=1;
    fusePionMass=0;
    fuseDeltaPhiBGSubtr=1;
    fparticleEtaCut=0.9;
    fparticlePtCut=0.15;
    fleadingJetCut=20;
    fsubleadingJetCut=20;
    fMinJetPt=5;
    fconstituentCut=5;
    fdeltaPhiCut=2;
    fmatchingR=0.3;
    ftrackingIneff=0.0;
    std::vector<double> fcentralityBins = {0,10,20,40,60,80};

    OpenFile(1);
    fJOutput = gDirectory;
    fJOutput->cd();

    //fhistos = new AliAnalysisTaskEmcalJetDijetMassHisto();
    /*
    fhistos = new AliJCDijetHistos();
    fhistos->SetName("jcdijet");
    fhistos->SetCentralityBinsHistos(fcentralityBins);
    fhistos->CreateEventTrackHistos();
    fhistos->fHMG->Print();
    */

    /*
    fEmcalHistos = new AliAnalysisTaskEmcalJetDijetMassHistos();
    fEmcalHistos->SetCentralityBinsHistos(fcentralityBins);
    fEmcalHistos->SetName("jcdijet");
    fEmcalHistos->CreateHistos();
    */

    hTest = new TH1D("test","test",4,0.0,10.0);


    //Work in progress
    const int fNCentBin = fcentralityBins.size()-1;
    /*
    TH1D *fh_events[fNCentBin];
    for (int i=0; i<fNCentBin; i++) {
        fh_events[i] = new TH1D(Form("h_eventsCentBin%02d",i), Form("h_eventsCentBin%02d",i), 40, 0.0, 40.0 );
    }
    */

    AllocateTrackHistograms();

    fana = new AliJEmcalDijetAna();
    fana->SetParticleCollArray(fParticleCollArray);
    fana->SetCentBins(fNCentBin);

    TString sktScheme;
    TString santiktScheme;
    switch (fktScheme) {
        case 0:  sktScheme = "E_scheme";
                 break;
        case 1:  sktScheme = "pt_scheme";
                 break;
        case 2:  sktScheme = "pt2_scheme";
                 break;
        case 3:  sktScheme = "Et_scheme";
                 break;
        case 4:  sktScheme = "Et2_scheme";
                 break;
        case 5:  sktScheme = "BIpt_scheme";
                 break;
        case 6:  sktScheme = "BIpt2_scheme";
                 break;
        default: sktScheme = "Unknown, check macro arguments!";
                 break;
    }
    switch (fantiktScheme) {
        case 0:  santiktScheme = "E_scheme";
                 break;
        case 1:  santiktScheme = "pt_scheme";
                 break;
        case 2:  santiktScheme = "pt2_scheme";
                 break;
        case 3:  santiktScheme = "Et_scheme";
                 break;
        case 4:  santiktScheme = "Et2_scheme";
                 break;
        case 5:  santiktScheme = "BIpt_scheme";
                 break;
        case 6:  santiktScheme = "BIpt2_scheme";
                 break;
        default: santiktScheme = "Unknown, check macro arguments!";
                 break;
    }

    cout << endl;
    cout << "===========SETTINGS===========" << endl;
    cout << "MC:                         " << fIsMC << endl;
    cout << "Centrality bins:            ";
    for(unsigned i=0; i< fcentralityBins.size(); i++) cout << fcentralityBins.at(i) << " ";
    cout << endl;
    cout << "Jet cone size:              " << fjetCone << endl;
    cout << "kt-jet cone size:           " << fktJetCone << endl;
    cout << "Using kt-jet scheme:        " << sktScheme.Data() << endl;
    cout << "Using antikt-jet scheme:    " << santiktScheme.Data() << endl;
    cout << "Using pion mass:            " << fusePionMass << endl;
    cout << "Using DeltaPhi in BG subtr: " << fuseDeltaPhiBGSubtr << endl;
    cout << "Particle eta cut:           " << fparticleEtaCut << endl;
    cout << "Particle pt cut:            " << fparticlePtCut << endl;
    cout << "Dijet leading jet cut:      " << fleadingJetCut << endl;
    cout << "Dijet subleading jet cut:   " << fsubleadingJetCut << endl;
    cout << "Jet min pt cut:             " << fMinJetPt << endl;
    cout << "Jet leading const. cut:     " << fconstituentCut << endl;
    cout << "Dijet DeltaPhi cut:         pi/" << fdeltaPhiCut << endl;
    cout << "Matching R for MC:          " << fmatchingR << endl;
    cout << "Tracking ineff for DetMC:   " << ftrackingIneff << endl;
    cout << endl;

    if(fusePionMass && (fktScheme!=0 || fantiktScheme!=0)) {
        cout << "Warning: Using pion mass for jets but not using E_scheme!" << endl;
        cout << endl;
    }

    TIter next(fHistManager.GetListOfHistograms());
    TObject* obj = 0;
    while ((obj = next())) {
        fOutput->Add(obj);
    }


#if !defined(__CINT__) && !defined(__MAKECINT__)
    fana->SetSettings(fDebug,
                      fparticleEtaCut,
                      fparticlePtCut,
                      fjetCone,
                      fktJetCone,
                      fktScheme,
                      fantiktScheme,
                      fusePionMass,
                      fuseDeltaPhiBGSubtr,
                      fconstituentCut,
                      fleadingJetCut,
                      fsubleadingJetCut,
                      fMinJetPt,
                      fdeltaPhiCut,
                      fmatchingR,
                      0.0); //Tracking ineff only for det level.

    //fana->InitHistos(fIsMC);
#endif




    //PostData(1, fOutput); // Post data for ALL output slots > 0 here.
    PostData(1, fJOutput); // Post data for ALL output slots > 0 here.
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetDijetMass::FillHistograms()
{
    cout << "AliAnalysisTaskEmcalJetDijetMass::FillHistograms" << endl;
    //fhistos->fh_eventSel->Fill("events wo/ cuts",1.0);
    //fEmcalHistos->fh_pt[0]->Fill(1.0);
    //fh_events[0]->Fill(1.0);
    hTest->Fill(1.0);
    DoJetLoop();
    DoTrackLoop();
    //DoClusterLoop();
    //DoCellLoop();

    return kTRUE;
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetDijetMass::DoJetLoop()
{
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
            //fhistos->fh_jetPt[0][0]->Fill(jet->Pt());
            //cout << "Filling jet pt of " << jet->Pt() << endl;
            //fHistManager.FillTH1(histname, jet->Pt());

            histname = TString::Format("%s/histJetArea_%d", groupname.Data(), fCentBin);
            //fHistManager.FillTH1(histname, jet->Area());

            histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), fCentBin);
            //fHistManager.FillTH1(histname, jet->Phi());

            histname = TString::Format("%s/histJetEta_%d", groupname.Data(), fCentBin);
            //fHistManager.FillTH1(histname, jet->Eta());

            if (jetCont->GetRhoParameter()) {
                histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), fCentBin);
                //fHistManager.FillTH1(histname, jet->Pt() - jetCont->GetRhoVal() * jet->Area());
            }
        }
        histname = TString::Format("%s/histNJets_%d", groupname.Data(), fCentBin);
        //fHistManager.FillTH1(histname, count);
    }
}

/**
 * This function performs a loop over the reconstructed tracks
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetDijetMass::DoTrackLoop()
{
    AliClusterContainer* clusCont = GetClusterContainer(0);

    TString histname;
    TString groupname;
    UInt_t sumAcceptedTracks = 0;
    AliParticleContainer* partCont = 0;
    TIter next(&fParticleCollArray);
    while ((partCont = static_cast<AliParticleContainer*>(next()))) {
        groupname = partCont->GetName();
        UInt_t count = 0;
        //fana->CalculateJets(fInputList, fCentBin);
        for(auto part : partCont->accepted()) {
            if (!part) continue;
            count++;

            histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), fCentBin);
            //fhistos->fh_pt[0]->Fill(part->Pt());
            //cout << "Filling track pt of " << part->Pt() << endl;
            fHistManager.FillTH1(histname, part->Pt());
            //fhistos->fh_pt[0]->Print();

            histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), fCentBin);
            //fHistManager.FillTH1(histname, part->Phi());

            histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), fCentBin);
            //fHistManager.FillTH1(histname, part->Eta());

            if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
                const AliVTrack* track = static_cast<const AliVTrack*>(part);

                histname = TString::Format("%s/fHistDeltaEtaPt_%d", groupname.Data(), fCentBin);
                //fHistManager.FillTH1(histname, track->Pt(), track->Eta() - track->GetTrackEtaOnEMCal());

                histname = TString::Format("%s/fHistDeltaPhiPt_%d", groupname.Data(), fCentBin);
                //fHistManager.FillTH1(histname, track->Pt(), track->Phi() - track->GetTrackPhiOnEMCal());

                histname = TString::Format("%s/fHistDeltaPtvsPt_%d", groupname.Data(), fCentBin);
                //fHistManager.FillTH1(histname, track->Pt(), track->Pt() - track->GetTrackPtOnEMCal());

                if (clusCont) {
                    Int_t iCluster = track->GetEMCALcluster();
                    if (iCluster >= 0) {
                        AliVCluster* cluster = clusCont->GetAcceptCluster(iCluster);
                        if (cluster) {
                            histname = TString::Format("%s/fHistEoverPvsP_%d", groupname.Data(), fCentBin);
                            //fHistManager.FillTH2(histname, track->P(), cluster->GetNonLinCorrEnergy() / track->P());
                        }
                    }
                }
            }
        }
        sumAcceptedTracks += count;
        //cout << "sumAcceptedTracks: " << sumAcceptedTracks << endl;

        histname = TString::Format("%s/histNTracks_%d", groupname.Data(), fCentBin);
        //fHistManager.FillTH1(histname, count);
    }

    histname = "fHistSumNTracks";
    //fHistManager.FillTH1(histname, sumAcceptedTracks);
}

/*
 * This function allocates the histograms for basic tracking QA.
 * A set of histograms (pT, eta, phi, difference between kinematic properties
 * at the vertex and at the EMCal surface, number of tracks) is allocated
 * per each particle container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetDijetMass::AllocateTrackHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    groupname = partCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The track containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{track};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{track};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

      if (TClass(partCont->GetClassName()).InheritsFrom("AliVTrack")) {
        histname = TString::Format("%s/fHistDeltaEtaPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{#eta}_{track}^{vertex} - #it{#eta}_{track}^{EMCal};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, 50, -0.5, 0.5);

        histname = TString::Format("%s/fHistDeltaPhiPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{#phi}_{track}^{vertex} - #it{#phi}_{track}^{EMCal};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, 200, -2, 2);

        histname = TString::Format("%s/fHistDeltaPtvsPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{p}_{T,track}^{vertex} - #it{p}_{T,track}^{EMCal} (GeV/#it{c});counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, fNbins / 2, -fMaxBinPt/2, fMaxBinPt/2);

        histname = TString::Format("%s/fHistEoverPvsP_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{P}_{track} (GeV/#it{c});#it{E}_{cluster} / #it{P}_{track} #it{c};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, fNbins / 2, 0, 4);
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
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalJetDijetMass::ExecOnce()
{
    AliAnalysisTaskEmcalJet::ExecOnce();
}

/**
 * Run analysis code here, if needed.
% * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetDijetMass::Run()
{
    return kTRUE;
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmcalJetDijetMass::Terminate(Option_t *) 
{
}



/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */

AliAnalysisTaskEmcalJetDijetMass *AliAnalysisTaskEmcalJetDijetMass::AddTaskEmcalJetDijetMass(
                                    const char *ntracks,
                                    const char *nclusters,
                                    const char* ncells,
                                    const char *suffix,
                                    TString taskName,
                                    Bool_t isMC,
                                    TString sJCatalyst,
                                    TString sJCatalystDetMC,
                                    UInt_t flags,
                                    TString centBins,
                                    double jetCone,
                                    double ktjetCone,
                                    int ktScheme,
                                    int antiktScheme,
                                    Bool_t usePionMass,
                                    Bool_t useDeltaPhiBGSubtr,
                                    double particleEtaCut,
                                    double particlePtCut,
                                    double leadingJetCut,
                                    double subleadingJetCut,
                                    double minJetPt,
                                    double constituentCut,
                                    double deltaPhiCut,
                                    double matchingR,
                                    double trackingIneff){
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AddTaskEmcalJetDijetMass", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    AliVEventHandler* handler = mgr->GetInputEventHandler();
    if (!handler)
    {
        ::Error("AddTaskEmcalJetDijetMass", "This task requires an input event handler");
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
    // Init the task and do settings
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

    TString name("AliAnalysisTaskEmcalJetDijetMass");
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

    //-------------------------------------------------------
    // Final settings, pass to manager and set the containers
    //-------------------------------------------------------
    // Load Custom Configuration and parameters
    // override values with parameters

    // flags can manipulate event selection:
    // 0: no additional events rejected.
    // AliAnalysisTask::DIJET_VERTEX13PA: use IsVertexSelected2013pA
    // AliAnalysisTask::DIJET_PILEUPSPD: use IsPileupFromSPD(3,0.6,3,2,5)
    // AliAnalysisTask::DIJET_UTILSPILEUPSPD: use IsPileUpSPD(InputEvent())
    // Combinations of these can be used by giving argument for example:
    // AliAnalysisTask::DIJET_VERTEX13PA|AliAnalysisTask::DIJET_PILEUPSPD

    // jet recombination schemes can be set with ktScheme argument:
    // E_scheme     = 0
    // pt_scheme    = 1
    // pt2_scheme   = 2
    // Et_scheme    = 3
    // Et2_scheme   = 4
    // BIpt_scheme  = 5
    // BIpt2_scheme = 6

    cout<<"AddTaskJCDijetTask::flags = "<<flags<<endl;

    std::stringstream ss( centBins.Data() );
    double binBorder;
    vector<double> vecCentBins;
    ss >> binBorder;
    while (!ss.fail()) {
        vecCentBins.push_back(binBorder);
        ss >> binBorder;
    }

    if (vecCentBins.size() < 2) {
        ::Error("AddTaskJCDijetTask", "Centrality bins are not properly set. At least two bin borders are needed. Terminating.");
        return NULL;
    }

    for (unsigned ivec=0; ivec < vecCentBins.size()-1; ivec++) {
        if(vecCentBins.at(ivec+1) <= vecCentBins.at(ivec)) {
            ::Error("AddTaskJCDijetTask", "Centrality bins are not properly set. Terminating.");
            return NULL;
        }
    }

    if (jetCone > 0.8 || jetCone < 0.0 || ktjetCone > 0.8 || ktjetCone < 0.0) {
        ::Error("AddTaskJCDijetTask", "Jet cones are set to be too small or too big. Terminating.");
        return NULL;
    }

    if (ktScheme < 0 || ktScheme > 6) {
        ::Error("AddTaskJCDijetTask", "Invalid ktScheme set. Please choose a setting from 0 till 6. Terminating.");
        return NULL;
    }
    if (antiktScheme < 0 || antiktScheme > 6) {
        ::Error("AddTaskJCDijetTask", "Invalid antiktScheme set. Please choose a setting from 0 till 6. Terminating.");
        return NULL;
    }
    cout << "MC: " << isMC << endl;

    AliAnalysisTaskEmcalJetDijetMass* dijetMassTask = new AliAnalysisTaskEmcalJetDijetMass(name);
    dijetMassTask->SetCaloCellsName(cellName);
    dijetMassTask->SetVzRange(-10,10);
    //dijetTask->SetDebugLevel(5);
    //dijetTask->SetJCatalystTaskName(sJCatalyst.Data());
    //dijetTask->SetJCatalystTaskNameDetMC(sJCatalystDetMC.Data());
    //dijetTask->SetCentralityBins(vecCentBins);
    //dijetTask->SetJetConeSize(jetCone, ktjetCone);
    //dijetTask->SetBGSubtrSettings(ktScheme, antiktScheme, usePionMass, useDeltaPhiBGSubtr);
    //dijetTask->SetIsMC(isMC);
    //dijetTask->SetCuts(particleEtaCut, particlePtCut, leadingJetCut, subleadingJetCut, constituentCut, deltaPhiCut, matchingR, trackingIneff, minJetPt);
    //dijetTask->AddFlags(flags);
    cout << dijetMassTask->GetName() << endl;

    if (trackName == "mcparticles") {
        dijetMassTask->AddMCParticleContainer(trackName);
    }
    else if (trackName == "tracks" || trackName == "Tracks") {
        dijetMassTask->AddTrackContainer(trackName);
    }
    else if (!trackName.IsNull()) {
        dijetMassTask->AddParticleContainer(trackName);
    }
    dijetMassTask->AddClusterContainer(clusName);


    //==== Set up the dijet task ====

    mgr->AddTask(dijetMassTask);

    // Create containers for input/output
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer()  ;
    TString contname(name);
    contname += "_histos";
    mgr->ConnectInput  (dijetMassTask, 0, cinput );
    AliAnalysisDataContainer *emcalHist = mgr->CreateContainer(Form("%scontainerList",name.Data()),
            TList::Class(), AliAnalysisManager::kOutputContainer,
            Form("%s:%s",AliAnalysisManager::GetCommonFileName(), name.Data()));
    AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",name.Data()),
            TDirectory::Class(), AliAnalysisManager::kOutputContainer,
            Form("%s:%s",AliAnalysisManager::GetCommonFileName(), name.Data()));
    mgr->ConnectOutput (dijetMassTask, 1, emcalHist );
    mgr->ConnectOutput (dijetMassTask, 1, jHist );
    //cout << "AliAnalysisDataContainer *jHist = mgr->CreateContainer(" << Form("%scontainer",dijetMassTask->GetName()) << "," << endl;
    //cout << "       TDirectory::Class(), AliAnalysisManager::kOutputContainer," << endl;
    //cout << "       " << Form("%s:%s",AliAnalysisManager::GetCommonFileName(), dijetMassTask->GetName()) << ");" << endl;

    return dijetMassTask;
}

