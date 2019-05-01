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
#include <THnSparse.h>
#include <TList.h>
#include <array>
#include <iostream>
#include <map>
#include <vector>
#include <TClonesArray.h>
#include <TGrid.h>
#include <THistManager.h>
#include <THashList.h>
#include <TLinearBinning.h>
#include <TObjArray.h>
#include <TParameter.h>

#include <AliVCluster.h>
#include <AliVParticle.h>
#include "AliAODTrack.h"
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"


#include "AliAnalysisUtils.h"
#include "AliESDEvent.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "Tracks/AliEmcalTriggerOfflineSelection.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliOADBContainer.h"


#include "AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA);
/// \endcond

using std::cout;
using std::endl;

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA() :
AliAnalysisTaskEmcalJet(),
fUseRecalcPatches(false),
fHistManager(),
//fRecevent(),
//fMCevent(),
//fGeneratorLevel(0),
//fMCJetContainer(nullptr),
//fRecoUtil(0x0),  fClusterEResolution(0x0), fVaryTrkPtRes(),
fUseSumw2(0),
fHistNumbJets(0x0),
fHistJetPt(0x0),
fHistJetJetPatchE(0x0),
fHistJetGammaPatchE(0x0),
fHistJetJetPatchPt(0x0),
fHistJetGammaPatchPt(0x0),
fHistTriggerPatchE(0x0)
//fhnMBJetSpectra(),
//fhnTrkQA(),
//fhnClusQA()

{
    //Array Initiation
    for(int i=0;i<9;i++){
        //fHistEMCalTowerMult[i] = NULL;
        
    }
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA(const char *name) :
AliAnalysisTaskEmcalJet(name, kTRUE),
fUseRecalcPatches(0x0),
fHistManager(name),
//fRecevent(NULL),
//fMCevent(NULL),
//fGeneratorLevel(0),
//fMCJetContainer(nullptr),
//fRecoUtil(0x0),  fClusterEResolution(0x0), fVaryTrkPtRes(0),
fUseSumw2(0),
fHistNumbJets(0x0),
fHistJetPt(0x0),
fHistJetJetPatchE(0x0),
fHistJetGammaPatchE(0x0),
fHistJetJetPatchPt(0x0),
fHistJetGammaPatchPt(0x0),
fHistTriggerPatchE(0x0)
//fhnMBJetSpectra(0x0),
//fhnTrkQA(0x0),
//fhnClusQA(0x0)


{
    //Array Initiation
    for(int i=0;i<=9;i++){
        //fHistEMCalTowerMult[i] = NULL;
        
    }
    SetMakeGeneralHistograms(kTRUE);
    SetCaloTriggerPatchInfoName("EmcalTriggers");
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::~AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA()
{
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::UserCreateOutputObjects()
{
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
    
    TString histName = "";
    TString histTitle = "";
     for(int i=0;i<=9;i++){
         histName = TString::Format("fHistEMCalTowerMult_%d",i);
         histTitle = TString::Format("%s;N_{tower};Counts",histName.Data());
         //fHistEMCalTowerMult[i] = new TH1F(histName.Data(),histTitle.Data(),400,0,400);
         //fOutput->Add(fHistEMCalTowerMult[i]);
     }
    
    // Get the MC particle branch, in case it exists
    
    AliJetContainer* jetMCCont = GetJetContainer();
    TString jetMCContName = jetMCCont->GetName();
    
    //fGeneratorLevel = GetMCParticleContainer("mcparticles");
   // if(fGeneratorLevel)
    //{
    //    Printf("MCparticleContainer: %s", fGeneratorLevel->GetName());
    //}
    //if(jetMCContName.Contains("mcparticles"))
    //{
    //    fMCJetContainer = jetMCCont;
    //}
    
    
    AllocateClusterHistograms();
    AllocateTrackHistograms();
    AllocateJetHistograms();
    AllocateCellHistograms();
    //AllocateParticleHistograms();
    
    
    TIter next(fHistManager.GetListOfHistograms());
    TObject* obj = 0;
    while ((obj = next())) {
        fOutput->Add(obj);
    }
    fOutput->Add(fHistNumbJets);
    fOutput->Add(fHistJetPt);
    
    
   // UInt_t SparseBit = 0;// Bit Code from GetDimParams()
    
 
    //SparseBit = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 |  1<<10 | 1<<11 | 1<<12 | 1<<13 | 1<<14 | 1<<15 | 1<<16;
    //fhnMBJetSpectra = NewTHnSparseF("fhnMBJetSpectra", SparseBit);
    //fhnMBJetSpectra->Sumw2();
   // fOutput->Add(fhnMBJetSpectra);
    
    if(fUseSumw2==0) {
        // =========== Switch on Sumw2 for all histos ===========
        for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
            TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
            if (h1){
                h1->Sumw2();
                continue;
            }
            THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
            if(hn)hn->Sumw2();
        }
    }
    
    PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/*
 * This function allocates the histograms for basic EMCal cluster QA.
 * A set of histograms (energy, eta, phi, number of cluster) is allocated
 * per each cluster container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::AllocateClusterHistograms()
{
    TString histname;
    TString histtitle;
    TString groupname;
    AliClusterContainer* clusCont = 0;
    TIter next(&fClusterCollArray);
    while ((clusCont = static_cast<AliClusterContainer*>(next()))) {
        groupname = clusCont->GetName();
        fHistManager.CreateHistoGroup(groupname);
        
        
        for (Int_t cent = 0; cent < fNcentBins; cent++) {
            histname = TString::Format("%s/histClusterEnergy_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E}_{cluster} (GeV);counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2, "s");
            
            histname = TString::Format("%s/histClusterEnergyExotic_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E}_{cluster}^{exotic} (GeV);counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2, "s");
            
            histname = TString::Format("%s/histClusterNonLinCorrEnergy_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E}_{cluster}^{non-lin.corr.} (GeV);counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2, "s");
            
            histname = TString::Format("%s/histClusterHadCorrEnergy_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E}_{cluster}^{had.corr.} (GeV);counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2, "s");
            
            histname = TString::Format("%s/histClusterPhi_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#phi}_{custer};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi(), "s");
            
            histname = TString::Format("%s/histClusterEta_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#eta}_{custer};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1, "s");
            
            histname = TString::Format("%s/fHistFcrossvEonline_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E}_{cluster} (GeV);#it{F}_{cross}", histname.Data());
            fHistManager.CreateTH2(histname, histtitle, 100, 0, 100, 100, 0, 1.0);
            
            
            histname = TString::Format("%s/histNClusters_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;number of clusters;events", histname.Data());
            if (fForceBeamType != kpp) {
                fHistManager.CreateTH1(histname, histtitle, 500, 0, 3000, "s");
            }
            else {
                fHistManager.CreateTH1(histname, histtitle, 200, 0, 200, "s");
            }
        }
    }
}

/*
 *  This function allocate histograms to the particle level MC event
 *  before the particles are propagated through GEANT.
 *  This function also has adaptability for the jet finder.
 *
 */

/*
void AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::AllocateParticleHistograms()
{
    // If MC, get the MC event
   // const AliMCEvent* mcevent = nullptr;
    //if (fGeneratorLevel) {
     //   mcevent = MCEvent();
    //}
    //else {
     //   return;
    //}
    TString histname;
    TString histtitle;
    TString groupname;
    Double_t pi = TMath::Pi();
    //AliParticleContainer* MCCont = 0;
    //MCCont = AddMCParticleContainer("mcparticles");
    //groupname = MCCont->GetName();
    fHistManager.CreateHistoGroup(groupname);
    
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/fHistParticleLvlpT_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt, "s");
    }
    
}
*/

/*
 * This function allocates the histograms for basic EMCal QA.
 * One 2D histogram with the cell energy spectra and the number of cells
 * per event is allocated per each centrality bin.
 */
void AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::AllocateCellHistograms()
{
    TString histname;
    TString histtitle;
    TString groupname(fCaloCellsName);
    
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
        histname = TString::Format("%s/histCellEnergy_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{E}_{cell} (GeV);counts", histname.Data());
        fHistManager.CreateTH1(histname, histtitle, 600, 0, 150, "s");
        
        histname = TString::Format("%s/histNCells_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;number of cells;events", histname.Data());
        if (fForceBeamType != kpp) {
            fHistManager.CreateTH1(histname, histtitle, 500, 0, 6000,"s");
        }
        else {
            fHistManager.CreateTH1(histname, histtitle, 200, 0, 200,"s");
        }
    }
}

/*
 * This function allocates the histograms for basic tracking QA.
 * A set of histograms (pT, eta, phi, difference between kinematic properties
 * at the vertex and at the EMCal surface, number of tracks) is allocated
 * per each particle container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::AllocateTrackHistograms()
{
    TString histname;
    TString histtitle;
    TString groupname;
    AliParticleContainer* partCont = 0;
    TIter next(&fParticleCollArray);
    while ((partCont = static_cast<AliParticleContainer*>(next()))) {
        groupname = partCont->GetName();
        fHistManager.CreateHistoGroup(groupname);
        
        
        for (Int_t cent = 0; cent < fNcentBins; cent++) {
            histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2, "s");
            
            histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#phi}_{track};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi(), "s");
            
            histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#eta}_{track};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1, "s");
            
            if (TClass(partCont->GetClassName()).InheritsFrom("AliVTrack")) {
                
                histname = TString::Format("%s/fHistTrackResolution_%d", groupname.Data(), cent);
                histtitle = TString::Format("%s;#sigma(#it{p}_{T}) / #it{p}_{T};counts", histname.Data());
                fHistManager.CreateTH1(histname, histtitle, 100, 0, 1);
                
                histname = TString::Format("%s/fHistCovMatrix_%d", groupname.Data(), cent);
                histtitle = TString::Format("%s; #it{p}_{T}; #sigma(#it{p}_{T}) / #it{p}_{T}; counts", histname.Data());
                fHistManager.CreateTH2(histname, histtitle, 150, 0, 150, 100,0,1);
                
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
                fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000, "s");
            }
            else {
                fHistManager.CreateTH1(histname, histtitle, 200, 0, 200, "s");
            }
        }
    }
}

/*
 * This function allocates the histograms for basic jet QA.
 * A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
 * per each jet container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::AllocateJetHistograms()
{
    TString histname;
    TString histtitle;
    TString groupname;
    Double_t pi = TMath::Pi();
    AliJetContainer* jetCont = 0;
    TIter next(&fJetCollArray);
    UInt_t SparseBit = 0;// Bit Code from GetDimParams()
    SparseBit = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 |  1<<10 | 1<<11 | 1<<12 | 1<<13 | 1<<14 | 1<<15 | 1<<16;
    
    while ((jetCont = static_cast<AliJetContainer*>(next()))) {
        groupname = jetCont->GetName();
        fHistManager.CreateHistoGroup(groupname);
        
        histname = TString::Format("%s/fhnJetSparse", groupname.Data());
        histtitle = histname + ";Multiplicity;#it{p}_{T,jet}^{uncorr} (GeV/#it{c});#it{p}_{T,leading} (GeV/#it{c};#it{E}_{leading} (GeV);#eta_{jet};#phi_{jet};m#it{z}_{leading};#it{A}_{jet};#it{NEF};#it{N}_{constit}";
        Int_t nbins13[10]  = {300, fNbins, fNbins, fNbins,24,36,20,20,20,100};
        Double_t min13[10] = {0,fMinBinPt,fMinBinPt,fMinBinPt,-1.2, -0.5*pi,0,0,0,0};
        Double_t max13[10] = {300,fMaxBinPt, fMaxBinPt, fMaxBinPt,1.2,1.5*pi,1,1,1,1};
        fHistManager.CreateTHnSparse(histname.Data(), histtitle.Data(), 10, nbins13, min13, max13, "s");
        
        
        for (Int_t cent = 0; cent < fNcentBins; cent++) {
            
            histname = TString::Format("%s/histJetPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt, "s");
            
            histname = TString::Format("%s/histJetClusterEnergy_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E}_{JetClus} (GeV);counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2, "s");
            
            histname = TString::Format("%s/histJetTrkPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{p}_{T,JetTrk} (GeV/#it{c});counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt, "s");
            
            histname = TString::Format("%s/histJetNEFvJetPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{p}_{T,Jet};NEF", histname.Data());
            fHistManager.CreateTH2(histname, histtitle, 200, 0, 200, 100, 0, 1, "s");
            
            histname = TString::Format("%s/histNumbJetConstvJetPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{p}_{T,Jet};NumberofJetConstit", histname.Data());
            fHistManager.CreateTH2(histname, histtitle, 200, 0, 200, 20,0,20, "s");
            
            histname = TString::Format("%s/histJetFF_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;Z=#it{p}_{T,leading} / #it{p}_{T,Jet};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, 100, 0, 1, "s");
            
            histname = TString::Format("%s/histJetZvJetPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{p}_{T,Jet};Z=#it{p}_{T,Trk} / #it{p}_{T,Jet}", histname.Data());
            fHistManager.CreateTH2(histname, histtitle, 200, 0, 200, 100,0,1);
            
            //histname = TString::Format("%s/histJetNeutralPtvJetPt_%d", groupname.Data(), cent);
            //histtitle = TString::Format("%s;#it{p}_{T,Jet};#it{p}_{T,NeutralConstituents}", histname.Data());
            //fHistManager.CreateTH2(histname, histtitle, 200, 0, 200, 200, 0, 200);
            
            histname = TString::Format("%s/histFCrossvZleading_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{z}_{leading};#it{F}_{Cross}", histname.Data());
            fHistManager.CreateTH2(histname, histtitle, 100, 0, 1, 100, 0, 1);
            
            //Matched Trigger Histos
            histname = TString::Format("%s/fHistNumbJETrigger_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;JE Trigger Normaliztion;counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, 1, 0, 1, "s");
            
            histname = TString::Format("%s/fHistNumbGATrigger_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;GA Trigger Normaliztion;counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, 1, 0, 1, "s");
            
            histname = TString::Format("%s/fHistJetJetPatchE_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E}_{JetGammaPatch} (GeV);counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt, "s");
            
            histname = TString::Format("%s/fHistJetGammaPatchE_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E}_{JetGammaPatch} (GeV);counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt, "s");
            
            histname = TString::Format("%s/fHistJetJetPatchPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{p}_{T,JetJetPatch} (GeV/#it{c});counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt, "s");
            
            histname = TString::Format("%s/fHistJetGammaPatchPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{p}_{T,JetGammaPatch} (GeV/#it{c});counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt, "s");
            
            histname = TString::Format("%s/fHistTriggerPatchE_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E}_{TriggerPatch} (GeV);counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt, "s");
            
            histname = TString::Format("%s/fHistDeltaEtaDeltaPhiJE_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#phi}_{jet} - #it{#phi}_{JE,patch};#it{#eta}_{Jet} - #it{#eta}_{JE,patch};counts", histname.Data());
            fHistManager.CreateTH2(histname, histtitle, 280, -0.7, 0.7, 280, -0.7, 0.7);
            
            histname = TString::Format("%s/fHistDeltaEtaDeltaPhiGA_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#phi}_{jet} - #it{#phi}_{GA,patch};#it{#eta}_{Jet} - #it{#eta}_{GA,patch};counts", histname.Data());
            fHistManager.CreateTH2(histname, histtitle, 280, -0.7, 0.7, 280, -0.7, 0.7);
            
            histname = TString::Format("%s/fHistJetEPatchEJE_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E_{JE,patch}}#it{E_{Jet,UnCorr}};counts", histname.Data());
            fHistManager.CreateTH2(histname, histtitle, 400, 0, 200, 400, 0, 200);
            
            histname = TString::Format("%s/fHistJetEPatchEGA_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{E_{GA,patch}};#it{E_{Jet,UnCorr}};counts", histname.Data());
            fHistManager.CreateTH2(histname, histtitle, 400, 0, 200, 400, 0, 200);
            //  End Matched Trigger Histos
            
            histname = TString::Format("%s/fHistJetTrackResolution_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#sigma(#it{p}_{T}) / #it{p}_{T};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, 100, 0, 1);
            
            histname = TString::Format("%s/fHistJetCovMatrix_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s; #it{p}_{T}; #sigma(#it{p}_{T}) / #it{p}_{T}; counts", histname.Data());
            fHistManager.CreateTH2(histname, histtitle, 150, 0, 150, 100,0,1);
            
            histname = TString::Format("%s/histJetArea_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{A}_{jet};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, 3);
            
            histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#phi}_{jet};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi(), "s");
            
            histname = TString::Format("%s/histJetEta_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#eta}_{jet};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1, "s");
            
            histname = TString::Format("%s/histJetClusterPhi_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#phi}_{JetClus};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi(), "s");
            
            histname = TString::Format("%s/histJetClusterEta_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{#eta}_{JetClus};counts", histname.Data());
            fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1, "s");
            
            histname = TString::Format("%s/histNJets_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;number of jets;events", histname.Data());
            if (fForceBeamType != kpp) {
                fHistManager.CreateTH1(histname, histtitle, 500, 0, 500, "s");
            }
            else {
                fHistManager.CreateTH1(histname, histtitle, 100, 0, 100, "s");
            }
            
            if (!jetCont->GetRhoName().IsNull()) {
                histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), cent);
                histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histname.Data());
                fHistManager.CreateTH1(histname, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
            }
        }
    }
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::FillHistograms()
{
    DoJetLoop();
    DoTrackLoop();
    DoClusterLoop();
    DoCellLoop();
    //DoParticleLoop();
    
    return kTRUE;
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::DoJetLoop()
{
    TString histname;
    TString groupname;
    AliJetContainer* jetCont = 0;
    TIter next(&fJetCollArray);
    
    //TClonesArray *triggerpatches = dynamic_cast<TClonesArray *>(InputEvent()->FindListObject("EmcalTriggers"));
    //if(!triggerpatches)
    //    AliErrorStream() << "Trigger patch container EmcalTriggers not found in task " << GetName() << std::endl;
    
    while ((jetCont = static_cast<AliJetContainer*>(next()))) {
        groupname = jetCont->GetName();
        UInt_t count = 0;
        for(auto jet : jetCont->accepted()) {
            if (!jet) continue;
            count++;

            TLorentzVector leadPart;
            jetCont->GetLeadingHadronMomentum(leadPart, jet);
            Double_t z = GetParallelFraction(leadPart.Vect(), jet);
            if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999;
            
            histname = TString::Format("%s/histJetPt_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, jet->Pt());
            
            histname = TString::Format("%s/histJetClusterEnergy_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, jet->E());
            
            histname = TString::Format("%s/histJetArea_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, jet->Area());
            
            histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, jet->Phi());
            
            histname = TString::Format("%s/histJetEta_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, jet->Eta());
            
            histname = TString::Format("%s/histJetNEFvJetPt_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH2(histname,jet->Pt(),jet->NEF() );
            
            
            histname = TString::Format("%s/histNumbJetConstvJetPt_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH2(histname,jet->Pt(),jet->N());
            
            histname = TString::Format("%s/histJetFF_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH1(histname, z );
            
            Double_t MatchedEta = 0.04, MatchedPhi = 0.04;
           //Look at  emcal trigger patches in event and associate with jet
            this->SetCaloTriggerPatchInfoName("EmcalTriggers");
            
            if (fTriggerPatchInfo) {
                TString objname(fTriggerPatchInfo->GetClass()->GetName());
                TClass cls(objname);
                if (!cls.InheritsFrom("AliEMCALTriggerPatchInfo")) {
                    AliError(Form("%s: Objects of type %s in %s are not inherited from AliEMCALTriggerPatchInfo!",GetName(), cls.GetName(), "EmcalTriggers"));
                    fTriggerPatchInfo = 0;
                }
            }
            
            if (!fTriggerPatchInfo) {
                AliError(Form("%s: Unable to get trigger patch container with name %s. Aborting", GetName(), "EmcalTriggers"));
                continue;
            }
            
            //AliEMCALTriggerPatchInfo *currentpatch(nullptr);
             for(auto p : *fTriggerPatchInfo){
                AliEMCALTriggerPatchInfo *currentpatch = static_cast<AliEMCALTriggerPatchInfo *>(p);
                
                //Select Single Shower Patch
                if(currentpatch->IsGammaLowRecalc() || currentpatch->IsGammaLowSimple()){
                    histname = TString::Format("%s/fHistDeltaEtaDeltaPhiGA_%d", groupname.Data(), fCentBin);
                    fHistManager.FillTH2(histname, jet->Phi() - currentpatch->GetPhiGeo(), jet->Eta() - currentpatch->GetEtaGeo());
                    histname = TString::Format("%s/fHistJetEPatchEGA_%d", groupname.Data(), fCentBin);
                    fHistManager.FillTH2(histname,jet->E(),currentpatch->GetPatchE());
                    
                    if(jet->Phi() - currentpatch->GetPhiGeo() <= MatchedPhi && jet->Eta() - currentpatch->GetEtaGeo() <= MatchedEta){
                        histname = TString::Format("%s/fHistNumbGATrigger_%d", groupname.Data(), fCentBin);
                        fHistManager.FillTH1(histname, 0.5);
                        histname = TString::Format("%s/fHistJetGammaPatchE_%d", groupname.Data(), fCentBin);
                        fHistManager.FillTH1(histname, currentpatch->GetPatchE());
                        histname = TString::Format("%s/fHistJetGammaPatchPt_%d", groupname.Data(), fCentBin);
                        fHistManager.FillTH1(histname, jet->Pt());
                    }
                }
                
                //Select Jet Patch
                if(currentpatch->IsJetLowRecalc() || currentpatch->IsJetLowSimple()){
                    histname = TString::Format("%s/fHistDeltaEtaDeltaPhiJE_%d", groupname.Data(), fCentBin);
                    fHistManager.FillTH2(histname, jet->Phi() - currentpatch->GetPhiGeo(), jet->Eta() - currentpatch->GetEtaGeo());
                    histname = TString::Format("%s/fHistJetEPatchEJE_%d", groupname.Data(), fCentBin);
                    fHistManager.FillTH2(histname,jet->E(),currentpatch->GetPatchE());
                    
                    if(jet->Phi() - currentpatch->GetPhiGeo() <= MatchedPhi && jet->Eta() - currentpatch->GetEtaGeo() <= MatchedEta){
                        histname = TString::Format("%s/fHistNumbJETrigger_%d", groupname.Data(), fCentBin);
                        fHistManager.FillTH1(histname, 0.5);
                        histname = TString::Format("%s/fHistJetJetPatchE_%d", groupname.Data(), fCentBin);
                        fHistManager.FillTH1(histname, currentpatch->GetPatchE());
                        histname = TString::Format("%s/fHistJetJetPatchPt_%d", groupname.Data(), fCentBin);
                        fHistManager.FillTH1(histname, jet->Pt());
                    }
                }
                
            }
            
             Double_t Z_part = 0.0;
            
            //Look at assoc tracks of a given jet
            AliParticleContainer* tracks = jetCont->GetParticleContainer();
            if (tracks) {
                for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
                    AliVParticle *JetTrk = jet->TrackAt(it, tracks->GetArray());
                    if (JetTrk) {
                        
                        AliVTrack *Vtrack = dynamic_cast<AliVTrack*>(JetTrk);
                        const AliAODTrack* AODJetTrack = static_cast<const AliAODTrack*>(Vtrack);
                        
                        if(!AODJetTrack) AliFatal("Not a standard AOD");
                        
                        AliExternalTrackParam exParam;
                        Double_t sigma = 0;
                        
                        //get covariance matrix
                        Double_t cov[21] = {0,};
                        AODJetTrack->GetCovMatrix(cov);
                        Double_t pxpypz[3] = {0,};
                        AODJetTrack->PxPyPz(pxpypz);
                        Double_t xyz[3] = {0,};
                        AODJetTrack->GetXYZ(xyz);
                        Short_t sign = AODJetTrack->Charge();
                        exParam.Set(xyz,pxpypz,cov,sign);
                        sigma = TMath::Sqrt(exParam.GetSigma1Pt2());
                        Z_part = Vtrack->Pt() / jet->Pt();
                        
                        histname = TString::Format("%s/fHistJetTrackResolution_%d", groupname.Data(), fCentBin);
                        fHistManager.FillTH1(histname, sigma);
                        histname = TString::Format("%s/fHistJetCovMatrix_%d", groupname.Data(), fCentBin);
                        fHistManager.FillTH2(histname, Vtrack->Pt(), sigma);
                        histname = TString::Format("%s/histJetTrkPt_%d", groupname.Data(), fCentBin);
                        fHistManager.FillTH1(histname, JetTrk->Pt());
                        histname = TString::Format("%s/histJetZvJetPt_%d", groupname.Data(), fCentBin);
                        fHistManager.FillTH2(histname,jet->Pt(),Z_part);
                        
                        //if(jetContName.Contains("mcparticles")){
                        //  fMCJetContainer
                        //  Double_t jetPID = 0;
                        //  jetPID->PID();
                       //   cout<<"PID of trk in Jet: "<<jetPID<<endl;
                        //}
                    }
                }
            }
            //Loop over Clusters in Jet
            AliClusterContainer* JetCluster = 0;
            JetCluster = jetCont->GetClusterContainer();
            if(JetCluster) {
                for(auto cluster : JetCluster->accepted()) {
                    AliTLorentzVector nPart;
                    cluster->GetMomentum(nPart, fVertex);
                    //Double_t JetFCross = GetFcross(cluster, fCaloCells);
                    Double_t JetFCross = 1.;
                    histname = TString::Format("%s/histJetClusterEnergy_%d", groupname.Data(), fCentBin);
                    fHistManager.FillTH1(histname, cluster->E());
                    histname = TString::Format("%s/histJetClusterPhi_%d", groupname.Data(), fCentBin);
                    fHistManager.FillTH1(histname, nPart.Phi_0_2pi());
                    histname = TString::Format("%s/histJetClusterEta_%d", groupname.Data(), fCentBin);
                    fHistManager.FillTH1(histname, nPart.Eta());
                    histname = TString::Format("%s/histFCrossvZleading_%d", groupname.Data(), fCentBin);
                    fHistManager.FillTH2(histname,jet->MaxTrackPt() / jet->Pt(), JetFCross);
                }
                
                
            }
            
            Double_t leadingtrackpT = 0.0;
            Double_t leadingclusterE = 0.0;
            Double_t JetFCrossLeading = 0.0;
            //AliVParticle *leadingParticle = 0x0;
            AliVCluster *leadingCluster = 0x0;
            
            
            leadingCluster = jet->GetLeadingCluster();
            if(leadingCluster){
                leadingclusterE = leadingCluster->E();
                //JetFCrossLeading = GetFcross(leadingCluster, fCaloCells);
                JetFCrossLeading = 1.;
            }
            
            AliVParticle* leadingTrk = jet->GetLeadingTrack();
            if(leadingTrk){
                leadingtrackpT = leadingTrk->Pt();
           }
            AliTrackContainer * Globaltracks = static_cast<AliTrackContainer * >(GetParticleContainer("tracks"));
            Double_t TrackMultiplicity = Globaltracks->GetNTracks();
            Double_t Numb = jet->N();
            Double_t x[10]={TrackMultiplicity,jet->Pt(),leadingtrackpT,leadingclusterE,jet->Eta(),jet->Phi(),z,jet->Area(),jet->NEF(),Numb};
            histname = TString::Format("%s/fhnJetSparse", groupname.Data());
            fHistManager.FillTHnSparse(histname, x);
            
            if (jetCont->GetRhoParameter()) {
                histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, jet->Pt() - jetCont->GetRhoVal() * jet->Area());
            }//Background Subtaction in PbPb
        }// Loop over each jet in a event
        histname = TString::Format("%s/histNJets_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, count);
    }
}

/**
 * This function performs a loop over the reconstructed tracks
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::DoTrackLoop()
{
    AliClusterContainer* clusCont = GetClusterContainer(0);
    
    AliTrackContainer* trackCont = GetTrackContainer("tracks");
   
    TString histname;
    TString groupname;
    AliParticleContainer* partCont = 0;
    TIter next(&fParticleCollArray);
    
    //for (auto trackIt : trackCont->accepted_momentum() ) {
    //}
    
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
                const AliAODTrack* AODTrack = static_cast<const AliAODTrack*>(track);
                
                if(!AODTrack) AliFatal("Not a standard AOD");
                
                AliExternalTrackParam exParam;
                Double_t sigma = 0;
                
                //get covariance matrix
                Double_t cov[21] = {0,};
                AODTrack->GetCovMatrix(cov);
                Double_t pxpypz[3] = {0,};
                AODTrack->PxPyPz(pxpypz);
                Double_t xyz[3] = {0,};
                AODTrack->GetXYZ(xyz);
                Short_t sign = AODTrack->Charge();
                exParam.Set(xyz,pxpypz,cov,sign);
                sigma = TMath::Sqrt(exParam.GetSigma1Pt2());
                
                histname = TString::Format("%s/fHistTrackResolution_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH1(histname, sigma);
                
                histname = TString::Format("%s/fHistCovMatrix_%d", groupname.Data(), fCentBin);
                fHistManager.FillTH2(histname, track->Pt(), sigma);
                
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
        
        histname = TString::Format("%s/histNTracks_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, count);
    }
}

/**
 * This function performs a loop over the reconstructed EMCal clusters
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::DoClusterLoop()
{
    Double_t TrkPt = 0.0, EovP = 0.0, trkpSum = 0.0;
    TString histname;
    TString groupname;
    AliClusterContainer* clusCont = 0;
    TIter next(&fClusterCollArray);
    while ((clusCont = static_cast<AliClusterContainer*>(next()))) {
        groupname = clusCont->GetName();
        
        //Cluster Exotics
        for(auto cluster : clusCont->all()) {
            if (!cluster) continue;
            
            //Double_t FCross = GetFcross(cluster, fCaloCells);
            Double_t FCross = 1.;
            histname = TString::Format("%s/fHistFcrossvEonline_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH2(histname, cluster->E(), FCross);
            
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
            
            //Int_t nTracksMatched = 0;
            const AliVTrack* track = 0x0;
            
            for (Int_t itrack=0; itrack < cluster->GetNTracksMatched(); itrack++){
              track = dynamic_cast<AliVTrack*>(cluster->GetTrackMatched(itrack));
              if(track){
                TrkPt = track->Pt();
                trkpSum += track->P();
                EovP = cluster->GetNonLinCorrEnergy() / trkpSum;
        
              }
            }
            

            
        }
        
        histname = TString::Format("%s/histNClusters_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, count);
    }
}

/**
 * This function performs a loop over the reconstructed EMCal cells
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::DoCellLoop()
{
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

/**
 *  This function perform a loop ov generator level MC
 *  particles in a current event and fills the relevent histograms.
 *  Functionality still in development
 */


/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::ExecOnce()
{
    
    AliAnalysisTaskEmcalJet::ExecOnce();
}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::Run()
{
    //Int_t runNum = InputEvent()->GetRunNumber();
    
    //Bool_t PileUp = AliVEvent::IsPileupFromSPD();
    
    //Double_t zVertex=fVertex[2];
    
    this->SetCaloTriggerPatchInfoName("EmcalTriggers");
    
    if (fTriggerPatchInfo) {
        TString objname(fTriggerPatchInfo->GetClass()->GetName());
        TClass cls(objname);
        if (!cls.InheritsFrom("AliEMCALTriggerPatchInfo")) {
            AliError(Form("%s: Objects of type %s in %s are not inherited from AliEMCALTriggerPatchInfo!",GetName(), cls.GetName(), "EmcalTriggers"));
            fTriggerPatchInfo = 0;
        }
    }
    
    if (!fTriggerPatchInfo) {
        AliError(Form("%s: Unable to get trigger patch container with name %s. Aborting", GetName(), "EmcalTriggers"));
        return kFALSE;
    }
    for(auto p : *fTriggerPatchInfo){
        AliEMCALTriggerPatchInfo *recpatch = static_cast<AliEMCALTriggerPatchInfo *>(p);
        if (recpatch) {
        }
    }
    
    //fClusterEResolution = new TF1("fClusterEResolution","sqrt([0]^2+[1]^2*x+([2]*x)^2)*0.01");// 2010 EMCal Test Beam Resolution
    //fClusterEResolution->SetParameters(4.35,9.07,1.63);//Fit to test beam
    
    
    
    
    //IsLEDEvent();
    
    return kTRUE;
}

/*
void AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::ExtractMainPatch() {
    
    this->SetCaloTriggerPatchInfoName("EmcalTriggers");
    
    if (fTriggerPatchInfo) {
        TString objname(fTriggerPatchInfo->GetClass()->GetName());
        TClass cls(objname);
        if (!cls.InheritsFrom("AliEMCALTriggerPatchInfo")) {
            AliError(Form("%s: Objects of type %s in %s are not inherited from AliEMCALTriggerPatchInfo!",GetName(), cls.GetName(), "EmcalTriggers"));
            fTriggerPatchInfo = 0;
        }
    }
    
    if (!fTriggerPatchInfo) {
        AliError(Form("%s: Unable to get trigger patch container with name %s. Aborting", GetName(), "EmcalTriggers"));
        return;
    }
    for(auto p : *fTriggerPatchInfo){
        AliEMCALTriggerPatchInfo *recpatch = static_cast<AliEMCALTriggerPatchInfo *>(p);
        if (recpatch) {
            cout<<"Pstchphi:"<<recpatch->GetPhiCM()<<endl;
        }
    }
    //Find main trigger
 
     (fTriggerPatchInfo) {
     TString objname(fTriggerPatchInfo->GetClass()->GetName());
     TClass cls(objname);
     if (!cls.InheritsFrom("AliEMCALTriggerPatchInfo")) {
     AliError(Form("%s: Objects of type %s in %s are not inherited from AliEMCALTriggerPatchInfo!",GetName(), cls.GetName(), "EmcalTriggers"));
     fTriggerPatchInfo = 0;
     }
     if (!fTriggerPatchInfo) {
     AliError(Form("%s: Unable to get trigger patch container with name %s. Aborting", GetName(), "EmcalTriggers"));
     return;
     }
     }
     //number of patches in event
     //Int_t nPatch = fTriggerPatchInfo->GetEntriesFast();
     //extract main trigger patch
     //Double_t emax = -1.;
     //for (Int_t iPatch = 0; iPatch < nPatch; iPatch++) {
     //   AliEMCALTriggerPatchInfo *patch = (AliEMCALTriggerPatchInfo*)fTriggerPatchInfo->At( iPatch );
     //  if (!patch) continue;
     // cout<<"PstchE:"<<patch->GetPatchE()<<endl;
     //}
 
}
*/

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::Terminate(Option_t *)
{
    cout<<"*****************************"<<endl;
    cout<<"******* Task Finished *******"<<endl;
    cout<<"*****************************"<<endl;
}
/*
Double_t AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz, const Double_t jetPx, const Double_t jetPy, const Double_t jetPz) const
{
    //
    // Get the z of a constituent inside of a jet
    //
    
    return (trkPx*jetPx+trkPy*jetPy+trkPz*jetPz)/(jetPx*jetPx+jetPy*jetPy+jetPz*jetPz);
}
 */
/*
Bool_t AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::IsLEDEvent() const
{
    Bool_t isLED = kFALSE;
    //
    // Check if the event is contaminated by LED signal
    //
    if (!fCaloCells) return kFALSE;
    
    TString histname;
    
    Int_t nCellCount[10] = {0};
    const Short_t ncells = fCaloCells->GetNumberOfCells();
    
    
    
    for (Short_t pos = 0; pos < ncells; pos++) {
        Int_t cellId = fCaloCells->GetCellNumber(pos);
        Double_t amp   = fCaloCells->GetAmplitude(pos);
        Int_t sMod = fGeom->GetSuperModuleNumber(cellId);
        
        if(amp<0.1) continue;
        nCellCount[sMod]++;
        
    }
    
    for(Int_t i=0; i<=9; i++){
        fHistEMCalTowerMult[i]->Fill(nCellCount[i]);
    }
    

    
    return isLED;
}
*/
/*
Double_t AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::GetFcross(const AliVCluster *cluster, AliVCaloCells *cells)
{
    Int_t    AbsIdseed  = -1;
    Double_t Eseed      = 0;
    for (Int_t i = 0; i < cluster->GetNCells(); i++) {
        if (cells->GetCellAmplitude(cluster->GetCellAbsId(i)) > Eseed) {
            Eseed     = cells->GetCellAmplitude(cluster->GetCellAbsId(i));
            AbsIdseed = cluster->GetCellAbsId(i);
        }
    }
    
    if (Eseed < 1e-9) {
        return 100;
    }
    
    Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1;
    fGeom->GetCellIndex(AbsIdseed,imod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi,iIeta,iphi,ieta);
    
    //Get close cells index and energy, not in corners
    
    Int_t absID1 = -1;
    Int_t absID2 = -1;
    
    if (iphi < AliEMCALGeoParams::fgkEMCALRows-1) {
        absID1 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta);
    }
    if (iphi > 0) {
        absID2 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta);
    }
    
    // In case of cell in eta = 0 border, depending on SM shift the cross cell index
    
    Int_t absID3 = -1;
    Int_t absID4 = -1;
    
    if (ieta == AliEMCALGeoParams::fgkEMCALCols-1 && !(imod%2)) {
        absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod+1, iphi, 0);
        absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod,   iphi, ieta-1);
    }
    else if (ieta == 0 && imod%2) {
        absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod,   iphi, ieta+1);
        absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod-1, iphi, AliEMCALGeoParams::fgkEMCALCols-1);
    }
    else {
        if (ieta < AliEMCALGeoParams::fgkEMCALCols-1) {
            absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta+1);
        }
        if (ieta > 0) {
            absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta-1);
        }
    }
    
    Double_t  ecell1 = cells->GetCellAmplitude(absID1);
    Double_t  ecell2 = cells->GetCellAmplitude(absID2);
    Double_t  ecell3 = cells->GetCellAmplitude(absID3);
    Double_t  ecell4 = cells->GetCellAmplitude(absID4);
    
    Double_t Ecross = ecell1 + ecell2 + ecell3 + ecell4;
    
    Double_t Fcross = 1 - Ecross/Eseed;
    
    return Fcross;
}

*/
THnSparse* AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::NewTHnSparseF(const char* name, UInt_t entries)
{
    Int_t count = 0;
    UInt_t tmp = entries;
    while(tmp!=0){
        count++;
        tmp = tmp &~ -tmp;  // clear lowest bit
    }
    
    TString hnTitle(name);
    const Int_t dim = count;
    Int_t nbins[dim];
    Double_t xmin[dim];
    Double_t xmax[dim];
    
    Int_t i=0;
    Int_t c=0;
    while(c<dim && i<32){
        if(entries&(1<<i)){
            
            TString label("");
            GetDimParams(i, label, nbins[c], xmin[c], xmax[c]);
            hnTitle += Form(";%s",label.Data());
            c++;
        }
        
        i++;
    }
    hnTitle += ";";
    
    return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
}
/**
 * Stores labels and binning of axes for creating a new THnSparseF.
 *
 * @param[in] iEntry Which axis to pick out of the cases.
 * @param[out] label Label of the axis.
 * @param[out] nbins Number of bins for the axis.
 * @param[out] xmin Minimum value of the axis.
 * @param[out] xmax Maximum value of the axis.
 */
void AliAnalysisTaskEmcalJetSpectra8TeVTriggerQA::GetDimParams(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
    const Double_t pi = TMath::Pi();
    
    switch(iEntry){
            
        case 0:
            label = "V0 centrality (%)";
            nbins = 30;
            xmin = 0.;
            xmax = 100.;
            // Adjust for pp, since we are retrieving multiplicity instead
            if (fForceBeamType == AliAnalysisTaskEmcal::kpp) {
                label = "Multiplicity";
                xmax = 300.;
            }
            break;
            
        case 1:
            label = "Jet p_{T}";
            nbins = 200;
            xmin = 0.;
            xmax = 200.;
            break;
            
        case 2:
            label = "Track p_{T}";
            nbins = 100;
            xmin = 0.;
            xmax = 100;
            break;
            
        case 3:
            label = "Cluster E";
            nbins = 100;
            xmin = 0.;
            xmax = 100;
            break;
            
        case 4:
            label = "Jet#eta";
            nbins = 24;
            xmin = -1.2;
            xmax = 1.2;
            break;
            
        case 5:
            label = "Jet #phi";
            nbins = 72;
            xmin = -0.5*pi;
            xmax = 1.5*pi;
            break;
            
        case 6:
            label = "#it{p}_{T,track}^{leading}";
            nbins = 100;
            xmin = 0;
            xmax = 100;
            break;
            
        case 7:
            label = "Matched Trigger Amp";
            nbins = 10;
            xmin = 0;
            xmax = 50;
            break;
            
        case 8:
            label = "#Delta#eta";
            nbins = 48;
            xmin = -1.2;
            xmax = 1.2;
            break;
            
        case 9:
            label = "#Delta#phi";
            nbins = 72;
            xmin = 0;
            xmax = 2*pi;
            break;
            
        case 10:
            label = "F_{cross}";
            nbins = 20;
            xmin = 0;
            xmax = 1;
            break;
            
        case 11:
            label = "z_{leading}";
            nbins = 20;
            xmin = 0;
            xmax = 1;
            break;
            
        case 12:
            label = "#it{A}_{jet}";
            nbins = 20;
            xmin = 0;
            xmax = 1;
            break;
            
        case 13:
            label = "Jet NEF";
            nbins = 20;
            xmin = 0;
            xmax = 1;
            break;
            
        case 14:
            label = "Numb Constit";
            nbins = 50;
            xmin = 0;
            xmax = 50;
            break;
            
        case 15:
            label = "Numb Chrg Constit";
            nbins = 50;
            xmin = 0;
            xmax = 50;
            break;
            
        case 16:
            label = "Numb Neu Constit";
            nbins = 50;
            xmin = 0;
            xmax = 50;
            break;
    }
}

