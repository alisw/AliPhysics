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
    //fHistManager(),
    fana(NULL),
    fhistos(NULL)
{
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalJetDijetMass::AliAnalysisTaskEmcalJetDijetMass(const char *name) : 
    AliAnalysisTaskEmcalJet(name, kTRUE),
    //fHistManager(name),
    fana(NULL),
    fhistos(NULL)
{
    SetMakeGeneralHistograms(kTRUE);
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalJetDijetMass::~AliAnalysisTaskEmcalJetDijetMass()
{
    delete fhistos;
    delete fana;
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalJetDijetMass::UserCreateOutputObjects()
{
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

    std::vector<double> fcentralityBins = {0,10,20,40,60,80};

    //fhistos = new AliAnalysisTaskEmcalJetDijetMassHisto();
    fhistos = new AliJCDijetHistos();
    fhistos->SetName("jcdijet");
    fhistos->SetCentralityBinsHistos(fcentralityBins);
    fhistos->CreateEventTrackHistos();
    fhistos->fHMG->Print();

    fana = new AliJCDijetAna();

#if !defined(__CINT__) && !defined(__MAKECINT__)
    //fana->SetSettings(fDebug,
    //                  fparticleEtaCut,
    //                  fparticlePtCut,
    //                  fjetCone,
    //                  fktJetCone,
    //                  fktScheme,
    //                  fantiktScheme,
    //                  fusePionMass,
    //                  fuseDeltaPhiBGSubtr,
    //                  fconstituentCut,
    //                  fleadingJetCut,
    //                  fsubleadingJetCut,
    //                  fMinJetPt,
    //                  fdeltaPhiCut,
    //                  fmatchingR,
    //                  0.0); //Tracking ineff only for det level.
#endif




    //TIter next(fHistManager.GetListOfHistograms());
    //TObject* obj = 0;
    //while ((obj = next())) {
    //    fOutput->Add(obj);
    //}

    PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetDijetMass::FillHistograms()
{
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
        for(auto part : partCont->accepted()) {
            if (!part) continue;
            count++;

            histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), fCentBin);
            //fHistManager.FillTH1(histname, part->Pt());

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

        histname = TString::Format("%s/histNTracks_%d", groupname.Data(), fCentBin);
        //fHistManager.FillTH1(histname, count);
    }

    histname = "fHistSumNTracks";
    //fHistManager.FillTH1(histname, sumAcceptedTracks);
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
    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
    TString contname(name);
    contname += "_histos";
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
            TList::Class(),AliAnalysisManager::kOutputContainer,
            Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectInput  (dijetMassTask, 0,  cinput1 );
    mgr->ConnectOutput (dijetMassTask, 1, coutput1 );

    return dijetMassTask;
}
