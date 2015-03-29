/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
 
// Analysis task for pi0-hadron correlation whis PHOS detector.
// Author:     Daniil Ponomarenko <Daniil.Ponomarenko@cern.ch>
// 20-Feb-2015

#include <Riostream.h>
#include "THashList.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TH3D.h"
#include "TMath.h"
#include "TVector3.h"
#include "TProfile.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliPHOSCorrelations.h"
#include "AliOADBContainer.h"
#include "AliInputEventHandler.h"

#include "AliPHOSGeometry.h"
#include "AliCentrality.h"
#include "AliEventplane.h"

#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"

#include "AliVTrack.h"
#include "AliCaloPhoton.h"

using std::cout;
using std::endl;

ClassImp(AliPHOSCorrelations)

//_______________________________________________________________________________
AliPHOSCorrelations::AliPHOSCorrelations()
:AliAnalysisTaskSE(),
    fPHOSGeo(0x0),
    fOutputContainer(0x0),
    fEvent(0x0),
    fEventESD(0x0),
    fEventAOD(0x0),
    fEventHandler(0),
    fCaloPhotonsPHOS(0x0),
    fTracksTPC(0x0),
    fCaloPhotonsPHOSLists(0x0),
    fTracksTPCLists(0x0),
    fRunNumber(-999),
    fInternalRunNumber(0),
    fPeriod("0x0"),
    fTriggerSelectionPbPb(kCentAndSCent),
    isREALEvent(false),
    fMBEvent(false),
    fNVtxZBins(5),
    fVtxEdges(10),
    fCentEdges(10),
    fCentNMixed(),
    fNEMRPBins(3),
    fAssocBins(),
    fVertexVector(),
    fVtxBin(0),
    fCentralityEstimator("V0M"),
    fCentrality(0.),
    fCentBin(0),
    fHaveTPCRP(0),
    fRP(0.),
    fEMRPBin(0),
    fEventPlaneMethod("V0"),
    fMaxAbsVertexZ(10.),
    fCentralityLowLimit(0.),
    fCentralityHightLimit(90),
    fMinClusterEnergy(0.3),
    fMinBCDistance(0),
    fMinNCells(3),
    fMinM02(0.2),
    fTOFCutEnabled(1),
    fTOFCut(100.e-9),
    fUseMassWindowParametrisation(true),
    fMassInvMeanMin(0.13),
    fMassInvMeanMax(0.145),
    fNSigmaWidth(0.),
    fUseEfficiency(true),
    fESDtrackCuts(0x0),
    fSelectHybridTracks(0),
    fTrackStatus(0),      
    fTrackFilterMask(786),
    fSelectSPDHitTracks(0),
    fSelectFractionTPCSharedClusters(kTRUE),
    fCutTPCSharedClustersFraction(0.4)
{
     // Constructor

    fMassMean[0]  = 1.00796e-05 ;
    fMassMean[1]  = 0.136096    ;
    fMassSigma[0] = 0.00383029 ;
    fMassSigma[1] = 0.0041709 ;
    fMassSigma[2] = 0.00468736 ;

    const Int_t nPtAssoc = 10 ;
    Double_t ptAssocBins[nPtAssoc] = {0.,0.5,1.0,1.5,2.0,3.,5.,7.,10.,16} ;
    fAssocBins.Set(nPtAssoc,ptAssocBins) ;

    const int nbins = 9;
    Double_t edges[nbins+1] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80.};
    TArrayD centEdges( nbins+1, edges );
    Int_t nMixed[nbins] = {4,4,6,10,20,30,50,100,100};
    TArrayI centNMixed(nbins, nMixed);
    SetCentralityBinning(centEdges, centNMixed);
    SetVertexBinning();

    fVertex[0] = 0; 
    fVertex[1] = 0; 
    fVertex[2] = 0;

    //SetGeometry();

    ZeroingVariables();
}

//_______________________________________________________________________________
AliPHOSCorrelations::AliPHOSCorrelations(const char *name)
:AliAnalysisTaskSE(name),
    fPHOSGeo(0x0),
    fOutputContainer(0x0),
    fEvent(0x0),
    fEventESD(0x0),
    fEventAOD(0x0),
    fEventHandler(0),
    fCaloPhotonsPHOS(0x0),
    fTracksTPC(0x0),
    fCaloPhotonsPHOSLists(0x0),
    fTracksTPCLists(0x0),
    fRunNumber(-999),
    fInternalRunNumber(0),
    fPeriod("0x0"),
    fTriggerSelectionPbPb(kCentAndSCent),
    isREALEvent(false),
    fMBEvent(false),
    fNVtxZBins(5),
    fVtxEdges(10),
    fCentEdges(10),
    fCentNMixed(),
    fNEMRPBins(3),
    fAssocBins(),
    fVertexVector(),
    fVtxBin(0),
    fCentralityEstimator("V0M"),
    fCentrality(0.),
    fCentBin(0),
    fHaveTPCRP(0),
    fRP(0.),
    fEMRPBin(0),
    fEventPlaneMethod("V0"),
    fMaxAbsVertexZ(10.),
    fCentralityLowLimit(0.),
    fCentralityHightLimit(90),
    fMinClusterEnergy(0.3),
    fMinBCDistance(0),
    fMinNCells(3),
    fMinM02(0.2),
    fTOFCutEnabled(1),
    fTOFCut(100.e-9),
    fUseMassWindowParametrisation(true),
    fMassInvMeanMin(0.13),
    fMassInvMeanMax(0.145),
    fNSigmaWidth(0.),
    fUseEfficiency(true),
    fESDtrackCuts(0x0),
    fSelectHybridTracks(0),
    fTrackStatus(0),      
    fTrackFilterMask(786),
    fSelectSPDHitTracks(0),
    fSelectFractionTPCSharedClusters(kTRUE),
    fCutTPCSharedClustersFraction(0.4)
{
    // Constructor
    // Output slots #0 write into a TH1 container
    DefineOutput(1,THashList::Class());

    fMassMean[0]  = 1.00796e-05 ;
    fMassMean[1]  = 0.136096    ;
    fMassSigma[0] = 0.00383029 ;
    fMassSigma[1] = 0.0041709 ;
    fMassSigma[2] = 0.00468736 ;

    const Int_t nPtAssoc = 10 ;
    Double_t ptAssocBins[nPtAssoc] = {0.,0.5,1.0,1.5,2.0,3.,5.,7.,10.,16} ;
    fAssocBins.Set(nPtAssoc,ptAssocBins) ;

    const int nbins = 9;
    Double_t edges[nbins+1] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80.};
    TArrayD centEdges( nbins+1, edges );
    Int_t nMixed[nbins] = {4,4,6,10,20,30,50,100,100};
    TArrayI centNMixed(nbins, nMixed);
    SetCentralityBinning(centEdges, centNMixed);
    SetVertexBinning();

    fVertex[0] = 0; 
    fVertex[1] = 0; 
    fVertex[2] = 0;

    //SetGeometry();

    ZeroingVariables();
}

//_______________________________________________________________________________
AliPHOSCorrelations::~AliPHOSCorrelations()
{
    if(fCaloPhotonsPHOS){ 
      delete fCaloPhotonsPHOS;
      fCaloPhotonsPHOS=0x0;
    }
    
    if(fTracksTPC){
      delete fTracksTPC;
      fTracksTPC=0x0;
    }

    if(fCaloPhotonsPHOSLists){
      fCaloPhotonsPHOSLists->SetOwner() ;
      delete fCaloPhotonsPHOSLists;
      fCaloPhotonsPHOSLists=0x0;
    }
    
    if(fTracksTPCLists){
      fTracksTPCLists->SetOwner() ;
      delete fTracksTPCLists;
      fTracksTPCLists=0x0 ;
    }
     
    if( fESDtrackCuts){      
      delete fESDtrackCuts;
      fESDtrackCuts=0x0 ;
    }
          
    if(fOutputContainer){
      delete fOutputContainer;
      fOutputContainer=0x0;
    }      
}

//_______________________________________________________________________________
void AliPHOSCorrelations::UserCreateOutputObjects()
{
    // Create histograms
    // Called once
    const Int_t     nRuns   = 200 ;
    const Int_t     ptMult  = 30 ;
    const Double_t  ptMin   = 0.  ;
    const Double_t  ptMax   = 30. ;

    // Create histograms
    if(fOutputContainer != NULL) { delete fOutputContainer; }
    fOutputContainer = new THashList();
    fOutputContainer->SetOwner(kTRUE);

    // Event selection
    fOutputContainer->Add(new TH1F( "hTriggerPassedEvents", "Event selection passed Cuts",                  60, 0., 60.  )) ;

    fOutputContainer->Add(new TH1F( "hMBTriggerperRunNumber", "MB trigger vs run number",                   200, 0., 200.)) ;
    fOutputContainer->Add(new TH1F( "hCentralTriggerperRunNumber", "Central trigger vs run number",         200, 0., 200.)) ;
    fOutputContainer->Add(new TH1F( "hSemiCentralTriggerperRunNumber", "SemiCentral trigger vs run number", 200, 0., 200.)) ;

    // Analysis event's progress
    fOutputContainer->Add(new TH1F( "hTotSelEvents", "Event selection", 15, 0., 15 )) ;
    fOutputContainer->Add(new TH2F( "hSelEvents", "Event selection", kTotalSelected+1, 0., double(kTotalSelected+1), nRuns,0., float(nRuns) )) ;

    // Centrality, Reaction plane and Vertex selection
    fOutputContainer->Add(new TH2F( "hCentrality", "Event centrality of all events",                100, 0., 100., nRuns,0., float(nRuns)   )) ;
    fOutputContainer->Add(new TH2F( "hCentralityTriggerEvent", "Event centrality trigger events",   100, 0., 100., nRuns,0., float(nRuns)   )) ;
    fOutputContainer->Add(new TH2F( "hCentralityMBEvent", "Event centrality MB events",             100, 0., 100., nRuns,0., float(nRuns)   )) ;
    fOutputContainer->Add(new TH2F( "phiRPflat", "RP distribution with TPC flat",                   100, 0., 2.*TMath::Pi(), 20, 0., 100.   )) ;

    fOutputContainer->Add(new TH1F( "hCentralityBining", "Event centrality bining of all events",            GetNumberOfCentralityBins(), 0., GetNumberOfCentralityBins() )) ;
    fOutputContainer->Add(new TH1F( "hCentralityBiningTrigger", "Event centrality bining of trigger events", GetNumberOfCentralityBins(), 0., GetNumberOfCentralityBins() )) ;
    fOutputContainer->Add(new TH1F( "hCentralityBiningMB", "Event centrality bining of MB events",           GetNumberOfCentralityBins(), 0., GetNumberOfCentralityBins() )) ;

    fOutputContainer->Add(new TH1F( "phiRPflatBining", "Event RP bining of all events",                     GetNumberOfRPBins(), 0., GetNumberOfRPBins() )) ;
    fOutputContainer->Add(new TH1F( "phiRPflatBiningTrigger", "Event RP bining of trigger events",          GetNumberOfRPBins(), 0., GetNumberOfRPBins() )) ;
    fOutputContainer->Add(new TH1F( "phiRPflatBiningMB", "Event RP bining of MB events",                    GetNumberOfRPBins(), 0., GetNumberOfRPBins() )) ;

    fOutputContainer->Add(new TH1F( "hVertexZBining", "Event vertex bining of all events",                  GetNumberOfVertexBins(), 0., GetNumberOfVertexBins() )) ;
    fOutputContainer->Add(new TH1F( "hVertexZBiningTrigger", "Event vertex bining of trigger events",       GetNumberOfVertexBins(), 0., GetNumberOfVertexBins() )) ;
    fOutputContainer->Add(new TH1F( "hVertexZBiningMB", "Event vertex bining of MB events",                 GetNumberOfVertexBins(), 0., GetNumberOfVertexBins() )) ;

    // Mass selection
    fOutputContainer->Add(new TH1F( "massWindowPass", "Mass selection", 10,     0.,     10. )) ;
    fOutputContainer->Add(new TH2F( "massWindow", "mean & sigma",       100.,   0.085,  0.185,  50,    0., 0.05 )) ;
    // Cluster multiplisity
    fOutputContainer->Add(new TH2F( "hCluEvsClu", "ClusterMult vs E",   200,    0.,     10.,    100,    0., 100. )) ;

    // Mixing progress
    // TODO: Fix filling maximum mixing depth value
    fOutputContainer->Add(new TH2F( "hCentralityMixingProgress", "Centrality Mixing Progress",  GetNumberOfCentralityBins(), 0., GetNumberOfCentralityBins(), 110, 0, 110 )) ;
    fOutputContainer->Add(new TH2F( "hVertexMixingProgress", "Vertex Mixing Progress",          GetNumberOfVertexBins(),     0., GetNumberOfVertexBins(),     110, 0, 110 )) ;
    fOutputContainer->Add(new TH2F( "hRPMixingProgress", "RP Mixing Progress",                  GetNumberOfRPBins(),         0., GetNumberOfRPBins(),         110, 0, 110 )) ;
                                              
    // Time Of Flight
    fOutputContainer->Add(new TH1F( "hTOFcut", "TOF PHOS's clastrers distribution", 10000, -5000, 5000 )) ;

    // Set hists, with track's and cluster's angle distributions.
    SetHistPtNumTrigger (ptMult, ptMin, ptMax);
    SetHistEtaPhi       (ptMult, ptMin, ptMax);
    SetHistPHOSClusterMap();
    SetHistMass         (ptMult, ptMin, ptMax);
    SetHistPtAssoc      (ptMult, ptMin, ptMax);

    // Setup photon lists
    Int_t kapacity = fNVtxZBins * GetNumberOfCentralityBins() * fNEMRPBins;
    fCaloPhotonsPHOSLists = new TObjArray(kapacity);
    fCaloPhotonsPHOSLists->SetOwner();

    fTracksTPCLists = new TObjArray(kapacity);
    fTracksTPCLists->SetOwner();

    PostData(1, fOutputContainer);
}

//_______________________________________________________________________________
void AliPHOSCorrelations::SetHistPtNumTrigger(const Int_t& ptMult, const Double_t& ptMin, const Double_t& ptMax) const
{
    TString spid[4]={"all","cpv","disp","both"} ;
    for(Int_t ipid=0; ipid<4; ipid++)    
    {
        fOutputContainer->Add(new TH1F(    Form("nTrigger_%s", spid[ipid].Data()), 
                                        Form("Num of trigger particle %s", spid[ipid].Data()), 
                                        ptMult, ptMin, ptMax ) );
        TH1F *h = static_cast<TH1F*>(fOutputContainer->Last()) ;
        h->Sumw2();
        h->GetXaxis()->SetTitle("Pt [GEV]");
    }
    for(Int_t ipid=0; ipid<4; ipid++)    
    {
        fOutputContainer->Add(new TH1F( Form("nTrigger_%s_MB", spid[ipid].Data()), 
                                        Form("Num of trigger particle %s", spid[ipid].Data()), 
                                        ptMult, ptMin, ptMax ) );
        TH1F *h = static_cast<TH1F*>(fOutputContainer->Last()) ;
        h->Sumw2();
        h->GetXaxis()->SetTitle("Pt [GEV]");
    }
}

//_______________________________________________________________________________
void AliPHOSCorrelations::SetHistEtaPhi(const Int_t& ptMult, const Double_t& ptMin, const Double_t& ptMax) const
{
    // Set hists, with track's and cluster's angle distributions.

    Float_t pi = TMath::Pi();

    //===
    fOutputContainer->Add(new TH2F( "clu_phieta","Cluster's #phi & #eta distribution", 
                            300, double(-1.8), double(-0.6), 
                            300, double(-0.2), double(0.2) ) );
    TH2F * h = static_cast<TH2F*>(fOutputContainer->Last()) ;
    h->GetXaxis()->SetTitle("#phi [rad]");
    h->GetYaxis()->SetTitle("#eta");

    //===
    fOutputContainer->Add(new TH2F( "clusingle_phieta","Cluster's  #phi & #eta distribution", 
                            300, double(-1.8), double(-0.6), 
                            300, double(-0.2), double(0.2) ) );
    h = static_cast<TH2F*>(fOutputContainer->Last()) ;
    h->GetXaxis()->SetTitle("#phi [rad]");
    h->GetYaxis()->SetTitle("#eta");

    //===
    fOutputContainer->Add(new TH3F( "track_phieta","TPC track's  #phi & #eta distribution", 
                            200, double(-pi-0.3), double(pi+0.3), 
                            200, double(-0.9), double(0.9), 
    ptMult, ptMin, ptMax ) );
    TH3F* h3 = static_cast<TH3F*>(fOutputContainer->FindObject("track_phieta")) ;
    h3->GetXaxis()->SetTitle("#phi [rad]");
    h3->GetYaxis()->SetTitle("#eta");
    h3->GetZaxis()->SetTitle("Pt [GeV/c]");
} 

//_______________________________________________________________________________
void AliPHOSCorrelations::SetHistMass(const Int_t& ptMult, const Double_t& ptMin, const Double_t& ptMax) const
{
    // Set mass histograms.

    Double_t binMult = 400;
    Double_t massMin = 0.0;
    Double_t massMax = 0.4;

    TString spid[4]={"all","cpv","disp","both"} ;

    TH2F * h;

    for(Int_t ipid=0; ipid<4; ipid++)    
    {
        // Real ++++++++++++++++++++++++++++++

        fOutputContainer->Add(new TH2F(Form("%s_mpt", spid[ipid].Data() ), "Real", 
                                binMult, massMin, massMax, 
                                ptMult, ptMin, ptMax ) );
        h = static_cast<TH2F*>(fOutputContainer->Last()) ;
        h->Sumw2();
        h->GetXaxis()->SetTitle("Mass [GeV]");
        h->GetYaxis()->SetTitle("Pt [GEV]");

        // MIX +++++++++++++++++++++++++

        fOutputContainer->Add(new TH2F(Form("mix_%s_mpt", spid[ipid].Data() ), "Mix", 
                                binMult, massMin, massMax, 
                                ptMult, ptMin, ptMax ) );
        h = static_cast<TH2F*>(fOutputContainer->Last()) ;
        h->Sumw2();
        h->GetXaxis()->SetTitle("Mass [GeV]");
        h->GetYaxis()->SetTitle("Pt [GEV]");
    }

    // Calibration PHOS Module Pi0peak {REAL}
    for(Int_t mod=1; mod<4; mod++)
    {
        fOutputContainer->Add(new TH2F(Form(  "both%d_mpt",mod), Form("Both cuts (CPV + Disp) mod[%d]",mod), 
                                binMult, massMin, massMax, 
                                ptMult, ptMin, ptMax ) );
        h = static_cast<TH2F*>(fOutputContainer->Last()) ;
        h->Sumw2();
        h->GetXaxis()->SetTitle("Mass [GeV]");
        h->GetYaxis()->SetTitle("Pt [GEV]");

        // Calibration PHOS Module Pi0peak {MIX}
        fOutputContainer->Add(new TH2F(Form(    "mix_both%d_mpt",mod), Form(" Both cuts (CPV + Disp) mod[%d]",mod), 
                                binMult, massMin, massMax, 
                                ptMult, ptMin, ptMax ) );
        h = static_cast<TH2F*>(fOutputContainer->Last()) ;
        h->Sumw2();
        h->GetXaxis()->SetTitle("Mass [GeV]");
        h->GetYaxis()->SetTitle("Pt [GEV]");

    }
}

//_______________________________________________________________________________
void AliPHOSCorrelations::SetHistPtAssoc(const Int_t& ptMult, const Double_t& ptMin, const Double_t& ptMax) const
{
    Double_t pi = TMath::Pi();

    Int_t PhiMult  =  60      ;
    Float_t PhiMin =  -0.5*pi ;
    Float_t PhiMax =  1.5*pi  ;
    Int_t EtaMult  =  20      ; 
    Float_t EtaMin = -1.      ;
    Float_t EtaMax =  1.      ;

    TString spid[4]={"all","cpv","disp","both"} ;

    for (int i = 0; i<fAssocBins.GetSize()-1; i++)
    {
        for(Int_t ipid=0; ipid<4; ipid++)
        {
            // Main histo for ConsiderPi0s().
            fOutputContainer->Add(new TH3F(Form("%s_ptphieta_ptAssoc_%3.1f", spid[ipid].Data(), fAssocBins.At(i+1)),
                                    Form("%s_ptphieta_ptAssoc_%3.1f", spid[ipid].Data(), fAssocBins.At(i+1)), 
                                    ptMult, ptMin, ptMax,  
                                    PhiMult, PhiMin, PhiMax, 
                                    EtaMult, EtaMin, EtaMax ) );
            TH3F * h = static_cast<TH3F*>(fOutputContainer->Last()) ;
            h->Sumw2();
            h->GetXaxis()->SetTitle("Pt_{triger} [GEV]");
            h->GetYaxis()->SetTitle("#phi [rad]");
            h->GetZaxis()->SetTitle("#eta");

            // For ConsiderPi0s_MBSelection().
            fOutputContainer->Add(new TH3F(Form("%s_ptphieta_ptAssoc_%3.1f_MB", spid[ipid].Data(), fAssocBins.At(i+1)),
                                    Form("%s_ptphieta_ptAssoc_%3.1f", spid[ipid].Data(), fAssocBins.At(i+1)), 
                                    ptMult, ptMin, ptMax,  
                                    PhiMult, PhiMin, PhiMax, 
                                    EtaMult, EtaMin, EtaMax ) );
            h = static_cast<TH3F*>(fOutputContainer->Last()) ;
            h->Sumw2();
            h->GetXaxis()->SetTitle("Pt_{triger} [GEV]");
            h->GetYaxis()->SetTitle("#phi [rad]");
            h->GetZaxis()->SetTitle("#eta");

            // For Mixed events in ConsiderTracksMix()
            fOutputContainer->Add(new TH3F(Form("mix_%s_ptphieta_ptAssoc_%3.1f", spid[ipid].Data(), fAssocBins.At(i+1)),
                                    Form("Mixed %s_ptphieta_ptAssoc_%3.1f", spid[ipid].Data(), fAssocBins.At(i+1)),
                                    ptMult, ptMin, ptMax,  
                                    PhiMult, PhiMin, PhiMax, 
                                    EtaMult, EtaMin, EtaMax ) );
            h = static_cast<TH3F*>(fOutputContainer->Last()) ;
            h->Sumw2();
            h->GetXaxis()->SetTitle("Pt_{triger} [GEV]");
            h->GetYaxis()->SetTitle("#phi [rad]");
            h->GetZaxis()->SetTitle("#eta");
        }
    }
}

//_______________________________________________________________________________
void AliPHOSCorrelations::SetHistPHOSClusterMap()
{
    //  Cluster X/Z/E distribution.
    for(int i =  0; i<5; i++)
    {
        fOutputContainer->Add(new TH3F( Form("QA_cluXZE_mod%i", i), Form("PHOS Clusters XZE distribution of module %i", i), 
                                70, 0, 70, 
                                60, 0, 60, 
                                200, 0, 20 ) );
        TH3F *h = static_cast<TH3F*>(fOutputContainer->Last()) ;
        h->GetXaxis()->SetTitle("X");
        h->GetYaxis()->SetTitle("Z");
        h->GetZaxis()->SetTitle("E");
    }    
}

//_______________________________________________________________________________
void AliPHOSCorrelations::UserExec(Option_t *) 
{
    // Main loop, called for each event analyze ESD/AOD 
    // Step 0: Event Objects
    LogProgress(0);

    fEvent = InputEvent();
    if( ! fEvent ) 
    {
        AliError("Event could not be retrieved");
        PostData(1, fOutputContainer);
        return ;
    }
    LogProgress(1);

    // Step 1(done once):  
    if( fRunNumber != fEvent->GetRunNumber() )
    {
        ShowTaskInfo();
        fRunNumber = fEvent->GetRunNumber();
        fInternalRunNumber = ConvertToInternalRunNumber(fRunNumber);
        SetGeometry();
        SetESDTrackCuts();
    }
    LogProgress(2);

    if(GetPeriod().Contains("0x0"))
    {
        AliWarning("Undefined period!");
    }

    // Step 2: Preparation variables for new event
    ZeroingVariables();

    fEventESD = dynamic_cast<AliESDEvent*>(fEvent);
    fEventAOD = dynamic_cast<AliAODEvent*>(fEvent);

    // Get Event-Handler for the trigger information
    fEventHandler= dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!fEventHandler) 
    {
        AliError("Could not get InputHandler");
        PostData(1, fOutputContainer);
        return; // Reject!
    }
    LogProgress(3);

    // Step 3: Event trigger selection
    // isREALEvent, fMBEvent
    if( RejectTriggerMaskSelection() ) 
    {
        PostData(1, fOutputContainer);
        return; // Reject!
    }
    LogProgress(4);

    // Step 4: Vertex
    // fVertex, fVertexVector, fVtxBin
    SetVertex();
    if( RejectEventVertex() ) 
    {
        PostData(1, fOutputContainer);
        return; // Reject!
    }
    LogProgress(5);

    // Step 5: Centrality
    // fCentrality, fCentBin
    SetCentrality(); 
    if( RejectEventCentrality() ) 
    {
        PostData(1, fOutputContainer);
        return; // Reject!
    }
    LogProgress(6);
    if(isREALEvent)  FillHistogram( "hCentralityTriggerEvent",    fCentrality, fInternalRunNumber-0.5 ) ;
    if(fMBEvent)     FillHistogram( "hCentralityMBEvent",        fCentrality, fInternalRunNumber-0.5 ) ;
    FillHistogram( "hCentrality", fCentrality, fInternalRunNumber-0.5 ) ;


    // Step 6: Reaction Plane
    // fHaveTPCRP, fRP, fRPV0A, fRPV0C, fRPBin
    EvalReactionPlane();  
    fEMRPBin = GetRPBin(); 

    // Step 7: Event Photons (PHOS Clusters) selection
    SelectPhotonClusters();
    if( ! fCaloPhotonsPHOS->GetEntriesFast() )    
        LogSelection(kHasPHOSClusters, fInternalRunNumber);

    // Step 8: Event Associated particles (TPC Tracks) selection
    SelectAccosiatedTracks();
    if( ! fTracksTPC->GetEntriesFast() ) LogSelection(kHasTPCTracks, fInternalRunNumber);
    LogSelection(kTotalSelected, fInternalRunNumber);

    // Step 9: Fill TPC's track mask and control bining hists.
    FillTrackEtaPhi();
    FillEventBiningProperties();
    LogProgress(7);

    // Step 10: Extract one most energetic pi0 candidate in this event.   
    SelectTriggerPi0ME();

    // Step 11: Start correlation analysis.
    if (isREALEvent)
    {
        ConsiderPi0s(); // Consider the most energetic Pi0 in this event with all tracks of this event.
        LogProgress(8);
    }

    if(GetPeriod().Contains("13") && fMBEvent)
    {
        ConsiderPi0s_MBSelection();
        LogProgress(9);
    }

     // Filling mixing histograms:
    ConsiderPi0sMix();      // Make background for extracting pi0 mass.
    ConsiderTracksMix();    // Compare only one most energetic pi0 candidate with all tracks from previous MB events.

    // Update pull
    if (fMBEvent)
    {
        // Update pull using MB events only!
        UpdatePhotonLists();    // Updating pull of photons.
        UpdateTrackLists();     // Updating pull of tracks.
        LogProgress(10);
    }

    LogProgress(11);
    // Post output data.
    PostData(1, fOutputContainer);
}

//_______________________________________________________________________________
void AliPHOSCorrelations::SetESDTrackCuts()
{
    if( fEventESD ) 
    {
        // Create ESD track cut
        fESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts() ;
        fESDtrackCuts->SetRequireTPCRefit(kTRUE) ;
    }
}

//_______________________________________________________________________________
Int_t AliPHOSCorrelations::ConvertToInternalRunNumber(const Int_t& run) const
{
    // Manual setup using data from logbook.
    if(GetPeriod().Contains("11h"))
    {
        switch(run)
        {
            case  170593 : return 179 ;
            case  170572 : return 178 ;
            case  170556 : return 177 ;
            case  170552 : return 176 ;
            case  170546 : return 175 ;
            case  170390 : return 174 ;
            case  170389 : return 173 ;
            case  170388 : return 172 ;
            case  170387 : return 171 ;
            case  170315 : return 170 ;
            case  170313 : return 169 ;
            case  170312 : return 168 ;
            case  170311 : return 167 ;
            case  170309 : return 166 ;
            case  170308 : return 165 ;
            case  170306 : return 164 ;
            case  170270 : return 163 ;
            case  170269 : return 162 ;
            case  170268 : return 161 ;
            case  170267 : return 160 ;
            case  170264 : return 159 ;
            case  170230 : return 158 ;
            case  170228 : return 157 ;
            case  170208 : return 156 ;
            case  170207 : return 155 ;
            case  170205 : return 154 ;
            case  170204 : return 153 ;
            case  170203 : return 152 ;
            case  170195 : return 151 ;
            case  170193 : return 150 ;
            case  170163 : return 149 ;
            case  170162 : return 148 ;
            case  170159 : return 147 ;
            case  170155 : return 146 ;
            case  170152 : return 145 ;
            case  170091 : return 144 ;
            case  170089 : return 143 ;
            case  170088 : return 142 ;
            case  170085 : return 141 ;
            case  170084 : return 140 ;
            case  170083 : return 139 ;
            case  170081 : return 138 ;
            case  170040 : return 137 ;
            case  170038 : return 136 ;
            case  170036 : return 135 ;
            case  170027 : return 134 ;
            case  169981 : return 133 ;
            case  169975 : return 132 ;
            case  169969 : return 131 ;
            case  169965 : return 130 ;
            case  169961 : return 129 ;
            case  169956 : return 128 ;
            case  169926 : return 127 ;
            case  169924 : return 126 ;
            case  169923 : return 125 ;
            case  169922 : return 124 ;
            case  169919 : return 123 ;
            case  169918 : return 122 ;
            case  169914 : return 121 ;
            case  169859 : return 120 ;
            case  169858 : return 119 ;
            case  169855 : return 118 ;
            case  169846 : return 117 ;
            case  169838 : return 116 ;
            case  169837 : return 115 ;
            case  169835 : return 114 ;
            case  169683 : return 113 ;
            case  169628 : return 112 ;
            case  169591 : return 111 ;
            case  169590 : return 110 ;
            case  169588 : return 109 ;
            case  169587 : return 108 ;
            case  169586 : return 107 ;
            case  169584 : return 106 ;
            case  169557 : return 105 ;
            case  169555 : return 104 ;
            case  169554 : return 103 ;
            case  169553 : return 102 ;
            case  169550 : return 101 ;
            case  169515 : return 100 ;
            case  169512 : return 99 ;
            case  169506 : return 98 ;
            case  169504 : return 97 ;
            case  169498 : return 96 ;
            case  169475 : return 95 ;
            case  169420 : return 94 ;
            case  169419 : return 93 ;
            case  169418 : return 92 ;
            case  169417 : return 91 ;
            case  169415 : return 90 ;
            case  169411 : return 89 ;
            case  169238 : return 88 ;
            case  169236 : return 87 ;
            case  169167 : return 86 ;
            case  169160 : return 85 ;
            case  169156 : return 84 ;
            case  169148 : return 83 ;
            case  169145 : return 82 ;
            case  169144 : return 81 ;
            case  169143 : return 80 ;
            case  169138 : return 79 ;
            case  169099 : return 78 ;
            case  169094 : return 77 ;
            case  169091 : return 76 ;
            case  169045 : return 75 ;
            case  169044 : return 74 ;
            case  169040 : return 73 ;
            case  169035 : return 72 ;
            case  168992 : return 71 ;
            case  168988 : return 70 ;
            case  168984 : return 69 ;
            case  168826 : return 68 ;
            case  168777 : return 67 ;
            case  168514 : return 66 ;
            case  168512 : return 65 ;
            case  168511 : return 64 ;
            case  168467 : return 63 ;
            case  168464 : return 62 ;
            case  168461 : return 61 ;
            case  168460 : return 60 ;
            case  168458 : return 59 ;
            case  168362 : return 58 ;
            case  168361 : return 57 ;
            case  168356 : return 56 ;
            case  168342 : return 55 ;
            case  168341 : return 54 ;
            case  168325 : return 53 ;
            case  168322 : return 52 ;
            case  168318 : return 51 ;
            case  168311 : return 50 ;
            case  168310 : return 49 ;
            case  168213 : return 48 ;
            case  168212 : return 47 ;
            case  168208 : return 46 ;
            case  168207 : return 45 ;
            case  168206 : return 44 ;
            case  168205 : return 43 ;
            case  168204 : return 42 ;
            case  168203 : return 41 ;
            case  168181 : return 40 ;
            case  168177 : return 39 ;
            case  168175 : return 38 ;
            case  168173 : return 37 ;
            case  168172 : return 36 ;
            case  168171 : return 35 ;
            case  168115 : return 34 ;
            case  168108 : return 33 ;
            case  168107 : return 32 ;
            case  168105 : return 31 ;
            case  168104 : return 30 ;
            case  168103 : return 29 ;
            case  168076 : return 28 ;
            case  168069 : return 27 ;
            case  168068 : return 26 ;
            case  168066 : return 25 ;
            case  167988 : return 24 ;
            case  167987 : return 23 ;
            case  167986 : return 22 ;
            case  167985 : return 21 ;
            case  167921 : return 20 ;
            case  167920 : return 19 ;
            case  167915 : return 18 ;
            case  167909 : return 17 ;
            case  167903 : return 16 ;
            case  167902 : return 15 ;
            case  167818 : return 14 ;
            case  167814 : return 13 ;
            case  167813 : return 12 ;
            case  167808 : return 11 ;
            case  167807 : return 10 ;
            case  167806 : return 9 ;
            case  167713 : return 8 ;
            case  167712 : return 7 ;
            case  167711 : return 6 ;
            case  167706 : return 5 ;
            case  167693 : return 4 ;
            case  166532 : return 3 ;
            case  166530 : return 2 ;
            case  166529 : return 1 ;

            default : return 199;
        }
    }
    
    if(GetPeriod().Contains("10h"))
    {
        switch(run)
        {
            case  139517 : return 137;
            case  139514 : return 136;
            case  139513 : return 135;
            case  139511 : return 134;
            case  139510 : return 133;
            case  139507 : return 132;
            case  139505 : return 131;
            case  139504 : return 130;
            case  139503 : return 129;
            case  139470 : return 128;
            case  139467 : return 127;
            case  139466 : return 126;
            case  139465 : return 125;
            case  139440 : return 124;
            case  139439 : return 123;
            case  139438 : return 122;
            case  139437 : return 121;
            case  139360 : return 120;
            case  139329 : return 119;
            case  139328 : return 118;
            case  139314 : return 117;
            case  139311 : return 116;
            case  139310 : return 115;
            case  139309 : return 114;
            case  139308 : return 113;
            case  139173 : return 112;
            case  139172 : return 111;
            case  139110 : return 110;
            case  139107 : return 109;
            case  139105 : return 108;
            case  139104 : return 107;
            case  139042 : return 106;
            case  139038 : return 105;
            case  139037 : return 104;
            case  139036 : return 103;
            case  139029 : return 102;
            case  139028 : return 101;
            case  138983 : return 100;
            case  138982 : return 99;
            case  138980 : return 98;
            case  138979 : return 97;
            case  138978 : return 96;
            case  138977 : return 95;
            case  138976 : return 94;
            case  138973 : return 93;
            case  138972 : return 92;
            case  138965 : return 91;
            case  138924 : return 90;
            case  138872 : return 89;
            case  138871 : return 88;
            case  138870 : return 87;
            case  138837 : return 86;
            case  138830 : return 85;
            case  138828 : return 84;
            case  138826 : return 83;
            case  138796 : return 82;
            case  138795 : return 81;
            case  138742 : return 80;
            case  138732 : return 79;
            case  138730 : return 78;
            case  138666 : return 77;
            case  138662 : return 76;
            case  138653 : return 75;
            case  138652 : return 74;
            case  138638 : return 73;
            case  138624 : return 72;
            case  138621 : return 71;
            case  138583 : return 70;
            case  138582 : return 69;
            case  138579 : return 68;
            case  138578 : return 67;
            case  138534 : return 66;
            case  138469 : return 65;
            case  138442 : return 64;
            case  138439 : return 63;
            case  138438 : return 62;
            case  138396 : return 61;
            case  138364 : return 60;
            case  138359 : return 59;
            case  138275 : return 58;
            case  138225 : return 57;
            case  138201 : return 56;
            case  138200 : return 55;
            case  138197 : return 54;
            case  138192 : return 53;
            case  138190 : return 52;
            case  138154 : return 51;
            case  138153 : return 50;
            case  138151 : return 49;
            case  138150 : return 48;
            case  138126 : return 47;
            case  138125 : return 46;
            case  137848 : return 45;
            case  137847 : return 44;
            case  137844 : return 43;
            case  137843 : return 42;
            case  137752 : return 41;
            case  137751 : return 40;
            case  137748 : return 39;
            case  137724 : return 38;
            case  137722 : return 37;
            case  137718 : return 36;
            case  137704 : return 35;
            case  137693 : return 34;
            case  137692 : return 33;
            case  137691 : return 32;
            case  137689 : return 31;
            case  137686 : return 30;
            case  137685 : return 29;
            case  137639 : return 28;
            case  137638 : return 27;
            case  137608 : return 26;
            case  137595 : return 25;
            case  137549 : return 24;
            case  137546 : return 23;
            case  137544 : return 22;
            case  137541 : return 21;
            case  137539 : return 20;
            case  137531 : return 19;
            case  137530 : return 18;
            case  137443 : return 17;
            case  137441 : return 16;
            case  137440 : return 15;
            case  137439 : return 14;
            case  137434 : return 13;
            case  137432 : return 12;
            case  137431 : return 11;
            case  137430 : return 10;
            case  137366 : return 9;
            case  137243 : return 8;
            case  137236 : return 7;
            case  137235 : return 6;
            case  137232 : return 5;
            case  137231 : return 4;
            case  137165 : return 3;
            case  137162 : return 2;
            case  137161 : return 1;
            default : return 199;
        }
    }
    
    if(GetPeriod().Contains("13"))
    {
        switch(run)
        {
            case  195344 : return 1;
            case  195346 : return 2;
            case  195351 : return 3;
            case  195389 : return 4;
            case  195390 : return 5;
            case  195391 : return 6;
            case  195478 : return 7;
            case  195479 : return 8;
            case  195480 : return 9;
            case  195481 : return 10;
            case  195482 : return 11;
            case  195483 : return 12;
            case  195529 : return 13;
            case  195531 : return 14;
            case  195532 : return 15;
            case  195566 : return 16;
            case  195567 : return 17;
            case  195568 : return 18;
            case  195592 : return 19;
            case  195593 : return 20;
            case  195596 : return 21;
            case  195633 : return 22;
            case  195635 : return 23;
            case  195644 : return 24;
            case  195673 : return 25;
            case  195675 : return 26;
            case  195676 : return 27;
            case  195677 : return 28;
            case  195681 : return 29;
            case  195682 : return 30;
            case  195720 : return 31;
            case  195721 : return 32;
            case  195722 : return 33;
            case  195724 : return 34;
            case  195725 : return 34;
            case  195726 : return 35;
            case  195727 : return 36;
            case  195760 : return 37;
            case  195761 : return 38;
            case  195765 : return 39;
            case  195767 : return 40;
            case  195783 : return 41;
            case  195787 : return 42;
            case  195826 : return 43;
            case  195827 : return 44;
            case  195829 : return 45;
            case  195830 : return 46;
            case  195831 : return 47;
            case  195867 : return 48;
            case  195869 : return 49;
            case  195871 : return 50;
            case  195872 : return 51;
            case  195873 : return 52;
            case  195935 : return 53;
            case  195949 : return 54;
            case  195950 : return 55;
            case  195954 : return 56;
            case  195955 : return 57;
            case  195958 : return 58;
            case  195989 : return 59;
            case  195994 : return 60;
            case  195998 : return 61;
            case  196000 : return 62;
            case  196006 : return 63;
            case  196085 : return 64;
            case  196089 : return 65;
            case  196090 : return 66;
            case  196091 : return 67;
            case  196099 : return 68;
            case  196105 : return 69;
            case  196107 : return 70;
            case  196185 : return 71;
            case  196187 : return 72;
            case  196194 : return 73;
            case  196197 : return 74;
            case  196199 : return 75;
            case  196200 : return 76;
            case  196201 : return 77;
            case  196203 : return 78;
            case  196208 : return 79;
            case  196214 : return 80;
            case  196308 : return 81;
            case  196309 : return 82;
            case  196310 : return 83;
            case  196311 : return 84;
            case  196433 : return 85;
            case  196474 : return 86;
            case  196475 : return 87;
            case  196477 : return 88;
            case  196528 : return 89;
            case  196533 : return 90;
            case  196535 : return 91;
            case  196563 : return 92;
            case  196564 : return 93;
            case  196566 : return 94;
            case  196568 : return 95;
            case  196601 : return 96;
            case  196605 : return 97;
            case  196608 : return 98;
            case  196646 : return 99;
            case  196648 : return 100;
            case  196701 : return 101;
            case  196702 : return 102;
            case  196703 : return 103;
            case  196706 : return 104;
            case  196714 : return 105;
            case  196720 : return 106;
            case  196721 : return 107;
            case  196722 : return 108;
            case  196772 : return 109;
            case  196773 : return 110;
            case  196774 : return 111;
            case  196869 : return 112;
            case  196870 : return 113;
            case  196874 : return 114;
            case  196876 : return 115;
            case  196965 : return 116;
            case  196967 : return 117;
            case  196972 : return 118;
            case  196973 : return 119;
            case  196974 : return 120;
            case  197003 : return 121;
            case  197011 : return 122;
            case  197012 : return 123;
            case  197015 : return 124;
            case  197027 : return 125;
            case  197031 : return 126;
            case  197089 : return 127;
            case  197090 : return 128;
            case  197091 : return 129;
            case  197092 : return 130;
            case  197094 : return 131;
            case  197098 : return 132;
            case  197099 : return 133;
            case  197138 : return 134;
            case  197139 : return 135;
            case  197142 : return 136;
            case  197143 : return 137;
            case  197144 : return 138;
            case  197145 : return 139;
            case  197146 : return 140;
            case  197147 : return 140;
            case  197148 : return 141;
            case  197149 : return 142;
            case  197150 : return 143;
            case  197152 : return 144;
            case  197153 : return 145;
            case  197184 : return 146;
            case  197189 : return 147;
            case  197247 : return 148;
            case  197248 : return 149;
            case  197254 : return 150;
            case  197255 : return 151;
            case  197256 : return 152;
            case  197258 : return 153;
            case  197260 : return 154;
            case  197296 : return 155;
            case  197297 : return 156;
            case  197298 : return 157;
            case  197299 : return 158;
            case  197300 : return 159;
            case  197302 : return 160;
            case  197341 : return 161;
            case  197342 : return 162;
            case  197348 : return 163;
            case  197349 : return 164;
            case  197351 : return 165;
            case  197386 : return 166;
            case  197387 : return 167;
            case  197388 : return 168;
            default : return 199;
        }
    }
    
    if(GetPeriod().Contains("0x0") && fDebug >= 1)
    {
        AliWarning("Period not defined");
    }
    return 1;
}

//_______________________________________________________________________________
Bool_t AliPHOSCorrelations::RejectTriggerMaskSelection()
{
    // Analyse trigger event and reject it if it not intresting.
    const Bool_t REJECT = true;
    const Bool_t ACCEPT = false;

    if( fDebug >= 2 )
        AliInfo( Form("Event passed offline phos trigger test: %s ", fEvent->GetFiredTriggerClasses().Data() ) );

    const Int_t physSelMask = fEventHandler->IsEventSelected();

    Bool_t isAny            = physSelMask & AliVEvent::kAny;            // to accept any trigger

    Bool_t isPHI1           = physSelMask & AliVEvent::kPHI1;           // PHOS trigger, CINT1 suite
    Bool_t isPHI7           = physSelMask & AliVEvent::kPHI7;           // PHOS trigger, CINT7 suite
    Bool_t isPHI8           = physSelMask & AliVEvent::kPHI8;           // PHOS trigger, CINT8 suite
    Bool_t isCentral        = physSelMask & AliVEvent::kCentral;        // PbPb central collision trigger
    Bool_t isSemiCentral    = physSelMask & AliVEvent::kSemiCentral;    // PbPb semicentral collision trigger
    Bool_t isPHOSPb         = physSelMask & AliVEvent::kPHOSPb;         // PHOS trigger, CINT8 suite for PbPb

    Bool_t isMB         = physSelMask & AliVEvent::kMB;                 // Minimum bias trigger, i.e. interaction trigger, offline SPD or V0 selection
    Bool_t isINT7       = physSelMask & AliVEvent::kINT7;               // V0AND trigger, offline V0 selection
    Bool_t isAnyINT     = physSelMask & AliVEvent::kAnyINT;             // to accept any interaction (aka minimum bias) trigger

    // Other triggers
    Bool_t isMUON           = physSelMask & AliVEvent::kMUON;           // Muon trigger, offline SPD or V0 selection
    Bool_t isHighMult       = physSelMask & AliVEvent::kHighMult;       // High-multiplicity trigger (threshold defined online), offline SPD or V0 selection
    Bool_t isEMC1           = physSelMask & AliVEvent::kEMC1;           // EMCAL trigger
    Bool_t isCINT5          = physSelMask & AliVEvent::kCINT5;          // Minimum bias trigger without SPD. i.e. interaction trigger, offline V0 selection
    Bool_t isCMUS5          = physSelMask & AliVEvent::kCMUS5;          // Muon trigger, offline V0 selection
    Bool_t isMUSPB          = physSelMask & AliVEvent::kMUSPB;          // idem for PbPb
    Bool_t isMUSH7          = physSelMask & AliVEvent::kMUSH7;          // Muon trigger: high pt, single muon, offline V0 selection, CINT7 suite
    Bool_t isMUSHPB         = physSelMask & AliVEvent::kMUSHPB;         // idem for PbPb
    Bool_t isMUL7           = physSelMask & AliVEvent::kMUL7;           // Muon trigger: like sign dimuon, offline V0 selection, CINT7 suite
    Bool_t isMuonLikePB     = physSelMask & AliVEvent::kMuonLikePB;     // idem for PbPb
    Bool_t isMUU7           = physSelMask & AliVEvent::kMUU7;           // Muon trigger, unlike sign dimuon, offline V0 selection, CINT7 suite
    Bool_t isMuonUnlikePB   = physSelMask & AliVEvent::kMuonUnlikePB;   // idem for PbPb
    Bool_t isEMC7           = physSelMask & AliVEvent::kEMC7;           // EMCAL trigger, CINT7 suite
    Bool_t isEMC8           = physSelMask & AliVEvent::kEMC8;           // EMCAL trigger, CINT8 suite
    Bool_t isMUS7           = physSelMask & AliVEvent::kMUS7;           // Muon trigger: low pt, single muon, offline V0 selection, CINT7 suite

    Bool_t isEMCEJE         = physSelMask & AliVEvent::kEMCEJE;                 // EMCAL jet patch trigger
    Bool_t isEMCEGA         = physSelMask & AliVEvent::kEMCEGA;                 // EMCAL gamma trigger

    Bool_t isDG5                = physSelMask & AliVEvent::kDG5;                // Double gap diffractive
    Bool_t isZED                = physSelMask & AliVEvent::kZED;                // ZDC electromagnetic dissociation
    Bool_t isSPI7               = physSelMask & AliVEvent::kSPI7;               // Power interaction trigger
    Bool_t isSPI                = physSelMask & AliVEvent::kSPI;                // Power interaction trigger
    Bool_t isINT8               = physSelMask & AliVEvent::kINT8;               // CINT8 trigger: 0TVX (T0 vertex) triger
    Bool_t isMuonSingleLowPt8   = physSelMask & AliVEvent::kMuonSingleLowPt8;   // Muon trigger : single muon, low pt, T0 selection, CINT8 suite
    Bool_t isMuonSingleHighPt8  = physSelMask & AliVEvent::kMuonSingleHighPt8;  // Muon trigger : single muon, high pt, T0 selection, CINT8 suite
    Bool_t isMuonLikeLowPt8     = physSelMask & AliVEvent::kMuonLikeLowPt8;     // Muon trigger : like sign muon, low pt, T0 selection, CINT8 suite
    Bool_t isMuonUnlikeLowPt8   = physSelMask & AliVEvent::kMuonUnlikeLowPt8;   // Muon trigger : unlike sign muon, low pt, T0 selection, CINT8 suite
    Bool_t isMuonUnlikeLowPt0   = physSelMask & AliVEvent::kMuonUnlikeLowPt0;   // Muon trigger : unlike sign muon, low pt, no additional L0 requirement
    Bool_t isUserDefined        = physSelMask & AliVEvent::kUserDefined;        // Set when custom trigger classes are set in AliPhysicsSelection, offline SPD or V0 selection
    Bool_t isTRD                = physSelMask & AliVEvent::kTRD;                // TRD trigger
    Bool_t isFastOnly           = physSelMask & AliVEvent::kFastOnly;           // The fast cluster fired. This bit is set in to addition another trigger bit, e.g. kMB

    // All input events
    FillHistogram("hTriggerPassedEvents", 0 );
    if ( isAny )        FillHistogram("hTriggerPassedEvents",  1.) ;        

    // PHOS events.
    if ( isPHI1 )       FillHistogram("hTriggerPassedEvents",  2. );    
    if ( isPHI7 )       FillHistogram("hTriggerPassedEvents",  3. );    
    if ( isPHI8 )       FillHistogram("hTriggerPassedEvents",  4. ); 
    if ( isCentral ) 
    {
        FillHistogram("hTriggerPassedEvents",   5. );
        FillHistogram("hCentralTriggerperRunNumber", fInternalRunNumber );
    }         
    if ( isSemiCentral ) 
    {
        FillHistogram("hTriggerPassedEvents",     6. ); 
        FillHistogram("hSemiCentralTriggerperRunNumber", fInternalRunNumber );
    }
    if ( isPHOSPb )     FillHistogram("hTriggerPassedEvents",     7. );    

    // MB events.
    if ( isMB )          
    {   
        FillHistogram("hTriggerPassedEvents",     8. );
        FillHistogram("hMBTriggerperRunNumber", fInternalRunNumber );
    }
    if ( isINT7 )       FillHistogram("hTriggerPassedEvents",     9. );
    if ( isAnyINT )     FillHistogram("hTriggerPassedEvents", 10. );

    Bool_t isTriggerEvent   = false;
    Bool_t isMIXEvent       = false;

    // TODO: Set false by default.
    isREALEvent  = false ;
    fMBEvent    = false ;

    // Choosing triggers for different periods.
    if(GetPeriod().Contains("13")) 
    {
        // Working: MB + TriggerPHOS events in real events; MB only in mixed events.
        isTriggerEvent  =  isPHI7 || isINT7 ;
        isMIXEvent      =  isINT7 ;
    }

    if(GetPeriod().Contains("11h"))
    {
        if(fTriggerSelectionPbPb == kCentAndSCent)
            isTriggerEvent = isCentral || isSemiCentral ;
        if(fTriggerSelectionPbPb == kCent)
             isTriggerEvent = isCentral ;
        if(fTriggerSelectionPbPb == kSCent)
             isTriggerEvent = isSemiCentral ;

        // Working: The same trigger in real and mixed events.
        isMIXEvent     = isCentral || isSemiCentral ;
    }

    // Select events.
    if(isTriggerEvent || isMIXEvent)
    {
        if ( isTriggerEvent )
        {
            FillHistogram("hTriggerPassedEvents", 17.);
            isREALEvent = true;
        }

        if ( isMIXEvent )
        {
            FillHistogram("hTriggerPassedEvents", 18.);
            fMBEvent = true;
        }

        return ACCEPT;
    }

    // other events
    FillHistogram("hTriggerPassedEvents",  19.); 

    //if ( isAny )      FillHistogram("hTriggerPassedEvents",  1.+20) ; // Not necessary
    // PHOS events.
    if ( isPHI1 )       FillHistogram("hTriggerPassedEvents",   2.+20 );
    if ( isPHI7 )       FillHistogram("hTriggerPassedEvents",   3.+20 );
    if ( isPHI8 )       FillHistogram("hTriggerPassedEvents",   4.+20 ); 
    if ( isCentral )    FillHistogram("hTriggerPassedEvents",   5.+20 );
    if ( isSemiCentral )FillHistogram("hTriggerPassedEvents",   6.+20 ); 
    if ( isPHOSPb )     FillHistogram("hTriggerPassedEvents",   7.+20 );
    // MB events.
    if ( isMB )         FillHistogram("hTriggerPassedEvents",   8.+20 );
    if ( isINT7 )       FillHistogram("hTriggerPassedEvents",   9.+20 );
    if ( isAnyINT )     FillHistogram("hTriggerPassedEvents",   10.+20 );
    // Other rejected events
    if ( isMUON )               FillHistogram("hTriggerPassedEvents",     11.+20 );
    if ( isHighMult )           FillHistogram("hTriggerPassedEvents",     12.+20 );
    if ( isEMC1 )               FillHistogram("hTriggerPassedEvents",     13.+20 );
    if ( isCINT5 )              FillHistogram("hTriggerPassedEvents",     14.+20 );
    if ( isCMUS5 )              FillHistogram("hTriggerPassedEvents",     15.+20 );
    if ( isMUSPB )              FillHistogram("hTriggerPassedEvents",     16.+20 );
    if ( isMUSH7 )              FillHistogram("hTriggerPassedEvents",     17.+20 );
    if ( isMUSHPB )             FillHistogram("hTriggerPassedEvents",     18.+20 );
    if ( isMUL7 )               FillHistogram("hTriggerPassedEvents",     19.+20 );
    if ( isMuonLikePB )         FillHistogram("hTriggerPassedEvents",     20.+20 );
    if ( isMUU7 )               FillHistogram("hTriggerPassedEvents",     21.+20 );
    if ( isMuonUnlikePB )       FillHistogram("hTriggerPassedEvents",     22.+20 );
    if ( isEMC7 )               FillHistogram("hTriggerPassedEvents",     23.+20 );
    if ( isEMC8 )               FillHistogram("hTriggerPassedEvents",     24.+20 );
    if ( isMUS7 )               FillHistogram("hTriggerPassedEvents",     25.+20 );
    if ( isEMCEJE )             FillHistogram("hTriggerPassedEvents",     26.+20 );
    if ( isEMCEGA )             FillHistogram("hTriggerPassedEvents",     27.+20 );
    if ( isDG5 )                FillHistogram("hTriggerPassedEvents",     28.+20 );
    if ( isZED )                FillHistogram("hTriggerPassedEvents",     29.+20 );
    if ( isSPI7 )               FillHistogram("hTriggerPassedEvents",     30.+20 );
    if ( isSPI )                FillHistogram("hTriggerPassedEvents",     31.+20 );
    if ( isINT8 )               FillHistogram("hTriggerPassedEvents",     32.+20 );
    if ( isMuonSingleLowPt8 )   FillHistogram("hTriggerPassedEvents",     33.+20 );
    if ( isMuonSingleHighPt8 )  FillHistogram("hTriggerPassedEvents",     34.+20 );
    if ( isMuonLikeLowPt8 )     FillHistogram("hTriggerPassedEvents",     35.+20 );
    if ( isMuonUnlikeLowPt8 )   FillHistogram("hTriggerPassedEvents",     36.+20 );
    if ( isMuonUnlikeLowPt0 )   FillHistogram("hTriggerPassedEvents",     37.+20 );
    if ( isUserDefined )        FillHistogram("hTriggerPassedEvents",     38.+20 );
    if ( isTRD )                FillHistogram("hTriggerPassedEvents",     39.+20 );
    if ( isFastOnly )           FillHistogram("hTriggerPassedEvents",     40.+20 );

    return REJECT;
}

//_______________________________________________________________________________
void AliPHOSCorrelations::SetVertex()
{
    const AliVVertex *primaryVertex = fEvent->GetPrimaryVertex();
    if( primaryVertex ) 
    {
        fVertex[0] = primaryVertex->GetX();
        fVertex[1] = primaryVertex->GetY();
        fVertex[2] = primaryVertex->GetZ();
    }
    else
    {
        //AliError("Event has 0x0 Primary Vertex, defaulting to origo");
        fVertex[0] = 0;
        fVertex[1] = 0;
        fVertex[2] = 0;
    }
    fVertexVector = TVector3(fVertex);

    fVtxBin=GetVertexBin(fVertexVector) ;
}

//_______________________________________________________________________________
Bool_t AliPHOSCorrelations::RejectEventVertex() const
{
    if( ! fEvent->GetPrimaryVertex() ) return true; // reject
    LogSelection(kHasVertex, fInternalRunNumber);

    if ( TMath::Abs(fVertexVector.z()) > fMaxAbsVertexZ ) return true; // reject
    LogSelection(kHasAbsVertex, fInternalRunNumber);

    return false; // accept event.
}

//_______________________________________________________________________________
void AliPHOSCorrelations::SetVertexBinning()
{
    // Define vertex bins by their edges
    const int nbins = fNVtxZBins+1;
    const double binWidth = 2*fMaxAbsVertexZ/fNVtxZBins;
    Double_t edges[nbins];
    for (int i = 0; i < nbins; ++i)
    {
        edges[i] = -1.*fMaxAbsVertexZ + binWidth*(double)i;
    }

    TArrayD vtxEdges(nbins, edges);

    for(int i=0; i<vtxEdges.GetSize()-1; ++i)
    {
        if(vtxEdges.At(i) > vtxEdges.At(i+1)) 
            AliFatal("edges are not sorted");

        fVtxEdges = vtxEdges;
    }
}

//_______________________________________________________________________________
Int_t AliPHOSCorrelations::GetVertexBin(const TVector3&  vertexVector)
{
    int lastBinUpperIndex = fVtxEdges.GetSize() -1;
    if( vertexVector.z() > fVtxEdges[lastBinUpperIndex] ) 
    {
        if( fDebug >= 1 )
            AliWarning( Form("vertex (%f) larger then upper edge of last vertex bin (%f)!", vertexVector.z(), fVtxEdges[lastBinUpperIndex]) );
        return lastBinUpperIndex-1;
    }
    if( vertexVector.z() < fVtxEdges[0] ) 
    {
        if( fDebug >= 1 )
        AliWarning( Form("vertex (%f) smaller then lower edge of first bin (%f)!", vertexVector.z(), fVtxEdges[0]) );
        return 0;
    }

    fVtxBin = TMath::BinarySearch<Double_t> ( GetNumberOfVertexBins(), fVtxEdges.GetArray(), vertexVector.z() );
    return fVtxBin;
}

//_______________________________________________________________________________
void AliPHOSCorrelations::SetCentrality()
{
    AliCentrality *centrality = fEvent->GetCentrality();
    if( centrality ) 
        fCentrality=centrality->GetCentralityPercentile(fCentralityEstimator);
    else 
    {
        AliError("Event has 0x0 centrality");
        fCentrality = -1.;
    }

    fCentBin = GetCentralityBin(fCentrality);
}

//_______________________________________________________________________________
Bool_t AliPHOSCorrelations::RejectEventCentrality() const
{
    if (fCentrality<fCentralityLowLimit)
        return true; //reject
    if(fCentrality>fCentralityHightLimit)
        return true; //reject

    return false;  // accept event.
}

//_______________________________________________________________________________
void AliPHOSCorrelations::SetCentralityBinning(const TArrayD& edges, const TArrayI& nMixed)
{
    // Define centrality bins by their edges
    for(int i=0; i<edges.GetSize()-1; ++i)
        if(edges.At(i) > edges.At(i+1)) AliFatal("edges are not sorted");
    if( edges.GetSize() != nMixed.GetSize()+1) AliFatal("edges and nMixed don't have appropriate relative sizes");
      
    fCentEdges = edges;
    fCentNMixed = nMixed;
}

//_______________________________________________________________________________
Int_t AliPHOSCorrelations::GetCentralityBin(const Float_t& centralityV0M)
{
    int lastBinUpperIndex = fCentEdges.GetSize() -1;
    if( centralityV0M > fCentEdges[lastBinUpperIndex] ) 
    {
        if( fDebug >= 1 )
            AliWarning( Form("centrality (%f) larger then upper edge of last centrality bin (%f)!", centralityV0M, fCentEdges[lastBinUpperIndex]) );
        return lastBinUpperIndex-1;
    }
    if( centralityV0M < fCentEdges[0] ) 
    {
        if( fDebug >= 1 )
        AliWarning( Form("centrality (%f) smaller then lower edge of first bin (%f)!", centralityV0M, fCentEdges[0]) );
        return 0;
    }

    fCentBin = TMath::BinarySearch<Double_t> ( GetNumberOfCentralityBins(), fCentEdges.GetArray(), centralityV0M );
    return fCentBin;
}

//_______________________________________________________________________________
void AliPHOSCorrelations::SetCentralityBorders(const Double_t& downLimit , const Double_t& upLimit )
{
    if (downLimit < 0. || upLimit > 100 || upLimit<=downLimit)
        AliError( Form("Warning. Bad value of centrality borders. Setting as default: fCentralityLowLimit=%2.f, fCentralityHightLimit=%2.f", fCentralityLowLimit, fCentralityHightLimit) );
    else
    {
        fCentralityLowLimit     = downLimit; 
        fCentralityHightLimit     = upLimit;
        AliInfo( Form("Centrality border was set as fCentralityLowLimit=%2.f, fCentralityHightLimit=%2.f", fCentralityLowLimit, fCentralityHightLimit ) );
    }
}

//_______________________________________________________________________________
void AliPHOSCorrelations::EvalReactionPlane()
{
    // assigns: fHaveTPCRP and fRP
    // also does RP histogram fill

    if(fEvent->GetEventplane())
    {
        Float_t ep =  fEvent->GetEventplane()->GetEventplane(GetEventPlaneMethod(), fEvent);

        if(GetEventPlaneMethod()=="Q" && (ep < 0 || ep > TMath::Pi()))
        {
            AliDebug(1,Form("Bad EP for <Q> method : %f\n",ep));
            ep = -1000;
        }
        else 
        if(GetEventPlaneMethod().Contains("V0")  )
        {
            if((ep > TMath::Pi()/2 || ep < -TMath::Pi()/2))
            {
                AliDebug(1,Form("Bad EP for <%s> method : %f\n",GetEventPlaneMethod().Data(), ep));
                ep = -1000;
            }

            ep+=TMath::Pi()/2; // put same range as for <Q> method
        }

        AliDebug(3,Form("Event plane angle %f",ep));
        
        if(ep>=999 || ep < 0.)
        {
            fRP = 0.;
        }
        else
            fRP = ep;
    }
    
    FillHistogram("phiRPflat",fRP,fCentrality) ;
}

//_______________________________________________________________________________
Int_t AliPHOSCorrelations::GetRPBin()
{
    Double_t averageRP;
    averageRP = fRP ;         // If possible, it is better to have EP bin from TPC
                            // to have similar events for miximng (including jets etc)   (fRPV0A+fRPV0C+fRP) /3.;
    fEMRPBin = Int_t(fNEMRPBins*(averageRP)/TMath::Pi());
    if(fEMRPBin > (Int_t)fNEMRPBins-1) 
        fEMRPBin = fNEMRPBins-1 ;
    else 
        if(fEMRPBin < 0) fEMRPBin=0;

    return fEMRPBin;
}

//_______________________________________________________________________________
void AliPHOSCorrelations::SelectPhotonClusters()
{
    //Selects PHOS clusters

    // clear (or create) array for holding events photons/clusters
    if(fCaloPhotonsPHOS)
        fCaloPhotonsPHOS->Clear();
    else
    {
        fCaloPhotonsPHOS = new TClonesArray("AliCaloPhoton",200);
        fCaloPhotonsPHOS->SetOwner();
    }

    Int_t inPHOS = 0 ;

    for (Int_t i = 0;  i < fEvent->GetNumberOfCaloClusters();  i++) 
    {
        AliVCluster *clu = fEvent->GetCaloCluster(i); 
        if (!clu->IsPHOS() || clu->E()< fMinClusterEnergy) continue; // reject cluster

        Float_t  position[3];
        clu->GetPosition(position);
        TVector3 global(position) ;
        Int_t relId[4] ;
        fPHOSGeo->GlobalPos2RelId(global,relId) ;
        Int_t modPHOS  = relId[0] ;
        Int_t cellXPHOS = relId[2];
        Int_t cellZPHOS = relId[3] ;
        
        Double_t distBC=clu->GetDistanceToBadChannel();
        if(distBC<fMinBCDistance)             continue ; // reject cluster
        if(clu->GetNCells() < fMinNCells)     continue ; // reject cluster
        if(clu->GetM02() < fMinM02)           continue ; // reject cluster

        FillHistogram("hTOFcut", clu->GetTOF()*1.e9 );

        if(fTOFCutEnabled)
        {
            Double_t tof = clu->GetTOF();
            if(TMath::Abs(tof) > fTOFCut ) continue ;
        }
        TLorentzVector lorentzMomentum;
        Double_t ecore = clu->GetCoreEnergy();
        //Double_t ecore = clu->E();

        FillHistogram("hCluEvsClu", clu->E(), clu->GetNCells()) ; 

        Double_t origo[3] = {0,0,0}; // don't rely on event vertex, assume (0,0,0) ?
        clu->GetMomentum(lorentzMomentum, origo);
    
        if(inPHOS>=fCaloPhotonsPHOS->GetSize())
            fCaloPhotonsPHOS->Expand(inPHOS+50) ;
        
        AliCaloPhoton * ph =new((*fCaloPhotonsPHOS)[inPHOS]) AliCaloPhoton(lorentzMomentum.X(), lorentzMomentum.Py(), lorentzMomentum.Z(), lorentzMomentum.E());
        inPHOS++ ;
        ph->SetCluster(clu);

        // Manual PHOS module number calculation
        /*Float_t cellId=clu->GetCellAbsId(0) ;
        Int_t mod = (Int_t)TMath:: Ceil(cellId/(56*64) ) ; */
        ph->SetModule(modPHOS) ;

        lorentzMomentum*=ecore/lorentzMomentum.E() ;

        //ph->SetNCells(clu->GetNCells());
        ph->SetMomV2(&lorentzMomentum) ;
        ph->SetDispBit(clu->GetDispersion() < 2.5) ;
        ph->SetCPVBit(clu->GetEmcCpvDistance() > 2.) ;

        FillHistogram(Form("QA_cluXZE_mod%i", modPHOS), cellXPHOS, cellZPHOS, lorentzMomentum.E() ) ;
    }
}

//_______________________________________________________________________________
void AliPHOSCorrelations::SelectAccosiatedTracks()
{
    // clear (or create) array for holding events tracks
    if(fTracksTPC)
        fTracksTPC->Clear();
    else 
    {
        fTracksTPC = new TClonesArray("TLorentzVector",12000);
    }
    Int_t iTracks = 0 ;
    for (Int_t i = 0; i < fEvent->GetNumberOfTracks(); i++) 
    {
      
        AliVTrack * track = (AliVTrack*)fEvent->GetTrack(i);

        //Select tracks under certain conditions, TPCrefit, ITSrefit ...
        ULong_t status = track->GetStatus();
        if (fTrackStatus && !((status & fTrackStatus) == fTrackStatus))
        {
            if( fDebug > 2 ) printf("\t Reject track, status != fTrackStatus\n" );
            continue ;
        }

        if(fEventESD)
        {
            if(!SelectESDTrack((AliESDtrack*)track)) continue ; // reject track
        }
        else
        {
            if(!SelectAODTrack((AliAODTrack*)track)) continue ;    // reject track  
        }

        Double_t px = track->Px();
        Double_t py = track->Py();
        Double_t pz = track->Pz();
        Double_t e  = track->E() ;
        
        if(iTracks >= fTracksTPC->GetSize())
            fTracksTPC->Expand(iTracks+50) ;
        
        new((*fTracksTPC)[iTracks]) TLorentzVector(px, py, pz, e);
        iTracks++ ;
    }
}

//_______________________________________________________________________________
void AliPHOSCorrelations::SelectTriggerPi0ME()
{
    const Int_t nPHOS = fCaloPhotonsPHOS->GetEntriesFast() ;
    for(Int_t i1 = 0; i1 < nPHOS-1; i1++)
    {
        AliCaloPhoton * ph1 = (AliCaloPhoton*)fCaloPhotonsPHOS->At(i1) ;
        for (Int_t i2 = i1+1; i2 < nPHOS; i2++)
        {
            AliCaloPhoton * ph2=(AliCaloPhoton*)fCaloPhotonsPHOS->At(i2) ;
            TLorentzVector p12 = *ph1 + *ph2;

            Double_t phiTrigger = p12.Phi() ;
            Double_t etaTrigger = p12.Eta() ;

            Double_t m      = p12.M() ;
            Double_t pt  = p12.Pt();
            Double_t eff = 1./GetEfficiency(pt);
            int mod1 = ph1->Module() ;
            int mod2 = ph2->Module() ;

            FillHistogram("clu_phieta",       phiTrigger, etaTrigger );
            FillHistogram("clusingle_phieta", ph1->Phi(), ph1->Eta() );
            FillHistogram("clusingle_phieta", ph2->Phi(), ph2->Eta() );

            FillHistogram("all_mpt", m, pt, eff );
      
            if ( ph1->IsCPVOK() && ph2->IsCPVOK() )
            {
                FillHistogram("cpv_mpt", m, pt, eff );
            }

            if ( ph1->IsDispOK() && ph2->IsDispOK() )
            {
                FillHistogram("disp_mpt", m, pt, eff );
                if ( ph1->IsCPVOK() && ph2->IsCPVOK() )
                {
                    FillHistogram("both_mpt", m, pt, eff );
                    if(mod1 == mod2) // for each module
                    {
                        FillHistogram(Form("both%d_mpt", mod1), m, pt, eff );
                    }
                }
            }

            if(!TestMass(m,pt)) continue; //reject this pair

            Int_t modCase = GetModCase(mod1, mod2);

            //Now we choosing most energetic pi0.
            TestPi0ME(kPidAll, p12, modCase);
            if ( ph1->IsCPVOK() && ph2->IsCPVOK() )
                TestPi0ME(kPidCPV, p12, modCase);
            if ( ph1->IsDispOK() && ph2->IsDispOK() )
            {
                TestPi0ME(kPidDisp, p12, modCase);
                if ( ph1->IsCPVOK() && ph2->IsCPVOK() )
                    TestPi0ME(kPidBoth, p12, modCase);
            }
        }
    }
}

//_______________________________________________________________________________
void AliPHOSCorrelations::ConsiderPi0s()
{
    TString spid[4] = {"all","cpv","disp","both"} ;
    // Counting number of trigger particles.
    for (int ipid = 0; ipid < 4; ipid++)
    {
        if (fMEExists[ipid])
            FillHistogram( Form("nTrigger_%s", spid[ipid].Data()), GetMEPt(ipid), 1./GetEfficiency(GetMEPt(ipid)) );
    }

    // Take track's angles and compare with trigger's angles.
    for(Int_t i3 = 0; i3 < fTracksTPC->GetEntriesFast(); i3++)
    {
        TLorentzVector * track = (TLorentzVector*)fTracksTPC->At(i3);

        Double_t phiAssoc = track->Phi();
        Double_t etaAssoc = track->Eta();
        Double_t ptAssoc  = track->Pt() ;

        Double_t ptAssocBin = GetAssocBin(ptAssoc) ;
        Double_t dPhi(0.), dEta(0.);

        for (int ipid = 0; ipid < 4; ipid++)
        {
            if (GetMEExists(ipid))
            {
                dPhi = GetMEPhi(ipid) - phiAssoc;
                while (dPhi > 1.5*TMath::Pi()) dPhi -= 2*TMath::Pi();
                while (dPhi < -.5*TMath::Pi()) dPhi += 2*TMath::Pi();
                dEta = GetMEEta(ipid) - etaAssoc;
                FillHistogram( Form("%s_ptphieta_ptAssoc_%3.1f", spid[ipid].Data(), ptAssocBin), GetMEPt(ipid), dPhi, dEta, 1./GetEfficiency(GetMEPt(ipid)) );
            }    
        }
    } 
}

//_______________________________________________________________________________
void AliPHOSCorrelations::ConsiderPi0s_MBSelection()
{
    TString spid[4] = {"all","cpv","disp","both"} ;
    // Counting number of trigger particles.
    for (int ipid = 0; ipid < 4; ipid++)
    {
        if (GetMEExists(ipid))
        {
            FillHistogram( Form("nTrigger_%s_MB", spid[ipid].Data()), GetMEPt(ipid), 1./GetEfficiency(GetMEPt(ipid)) );
        }
    }

    // Take track's angles and compare with trigger's angles.
    for(Int_t i3 = 0; i3 < fTracksTPC->GetEntriesFast(); i3++)
    {
        TLorentzVector * track = (TLorentzVector*)fTracksTPC->At(i3);

        Double_t phiAssoc = track->Phi();
        Double_t etaAssoc = track->Eta();
        Double_t ptAssoc  = track->Pt();

        Double_t ptAssocBin = GetAssocBin(ptAssoc) ;
        Double_t dPhi(0.), dEta(0.);

        for (int ipid = 0; ipid < 4; ipid++)
        {
            if (GetMEExists(ipid))
            {
                dPhi = GetMEPhi(ipid) - phiAssoc;
                while (dPhi > 1.5*TMath::Pi()) dPhi -= 2*TMath::Pi();
                while (dPhi < -.5*TMath::Pi()) dPhi += 2*TMath::Pi();
                dEta = GetMEEta(ipid) - etaAssoc;
                FillHistogram(Form("%s_ptphieta_ptAssoc_%3.1f_MB", spid[ipid].Data(), ptAssocBin),  GetMEPt(ipid), dPhi, dEta, 1./GetEfficiency(GetMEPt(ipid)) );
            }    
        }
    } 
}

//_______________________________________________________________________________
void AliPHOSCorrelations::ConsiderPi0sMix()
{
    TList * arrayList = GetCaloPhotonsPHOSList(fVtxBin, fCentBin, fEMRPBin);
    FillHistogram("hCentralityMixingProgress", fCentBin, arrayList->GetEntries() );
    FillHistogram("hVertexMixingProgress", fVtxBin, arrayList->GetEntries() );
    FillHistogram("hRPMixingProgress", fEMRPBin, arrayList->GetEntries() );
    for(Int_t evi = 0; evi < arrayList->GetEntries(); evi++)
    {
        TClonesArray * mixPHOS = static_cast<TClonesArray*>(arrayList->At(evi));
        for (Int_t i1 = 0; i1 < fCaloPhotonsPHOS->GetEntriesFast(); i1++)
        {
            AliCaloPhoton * ph1 = (AliCaloPhoton*)fCaloPhotonsPHOS->At(i1) ;
            for(Int_t i2 = 0; i2 < mixPHOS->GetEntriesFast(); i2++)
            {
                AliCaloPhoton * ph2 = (AliCaloPhoton*)mixPHOS->At(i2) ;
                TLorentzVector p12 = *ph1 + *ph2;
                Double_t m      = p12.M() ;
                Double_t pt  = p12.Pt() ;
                Double_t eff = 1./GetEfficiency(pt);
                
                int mod1 = ph1->Module() ;
                int mod2 = ph2->Module() ;

                FillHistogram("mix_all_mpt", m, pt, eff);
                if ( ph1->IsCPVOK() && ph2->IsCPVOK() ) 
                {
                    FillHistogram("mix_cpv_mpt",m, pt, eff);
                }
                if ( ph1->IsDispOK() && ph2->IsDispOK() )
                {
                    FillHistogram("mix_disp_mpt",m, pt, eff);
                    if ( ph1->IsCPVOK() && ph2->IsCPVOK() )
                    {
                        FillHistogram("mix_both_mpt",m, pt, eff);
                        if (mod1 == mod2) // for each module
                        {
                            FillHistogram(Form("mix_both%d_mpt",mod1),m, pt, eff);
                        }
                    }
                }
            }
        }
    }
}

//_______________________________________________________________________________
void AliPHOSCorrelations::ConsiderTracksMix()
{
    TString spid[4] = {"all","cpv","disp","both"} ;

    TList * arrayList = GetTracksTPCList(fVtxBin, fCentBin, fEMRPBin);

    for(Int_t evi = 0; evi < arrayList->GetEntries();evi++)
    {
        TClonesArray * mixTracks = static_cast<TClonesArray*>(arrayList->At(evi));
        for(Int_t i3 = 0; i3 < mixTracks->GetEntriesFast(); i3++)
        {
            TLorentzVector * track = (TLorentzVector*)mixTracks->At(i3);        

            Double_t phiAssoc = track->Phi();
            Double_t etaAssoc = track->Eta();
            Double_t ptAssoc  = track->Pt();

            Double_t ptAssocBin = GetAssocBin(ptAssoc) ;

            Double_t ptTrigger(0.);

            Double_t dPhi(0.), dEta(0.);

            for (int ipid = 0; ipid < 4; ipid++)
            {
                if (GetMEExists(ipid))
                {
                    dPhi = GetMEPhi(ipid) - phiAssoc;
                    while (dPhi > 1.5*TMath::Pi()) dPhi -= 2*TMath::Pi();
                    while (dPhi < -.5*TMath::Pi()) dPhi += 2*TMath::Pi();
                    dEta = GetMEEta(ipid) - etaAssoc;
                    ptTrigger = GetMEPt(ipid);

                    FillHistogram(Form("mix_%s_ptphieta_ptAssoc_%3.1f", spid[ipid].Data(), ptAssocBin), ptTrigger, dPhi, dEta, 1./GetEfficiency(ptTrigger));
                }    
            }
        }
    } 
}

//_______________________________________________________________________________
TList* AliPHOSCorrelations::GetCaloPhotonsPHOSList(const UInt_t vtxBin, const UInt_t centBin, const UInt_t rpBin)
{
    int offset = vtxBin * GetNumberOfCentralityBins() * fNEMRPBins + centBin * fNEMRPBins + rpBin;
    if( fCaloPhotonsPHOSLists->At(offset) ) 
    {
        // list exists
        TList* list = dynamic_cast<TList*> (fCaloPhotonsPHOSLists->At(offset));
        return list;
    }
    else
    { 
        // no list for this bin has been created, yet
        TList* list = new TList();
        fCaloPhotonsPHOSLists->AddAt(list, offset);
        return list;
    }
}

//_______________________________________________________________________________
TList* AliPHOSCorrelations::GetTracksTPCList(const UInt_t vtxBin, const UInt_t centBin, const UInt_t rpBin)
{
    int offset = vtxBin * GetNumberOfCentralityBins() * fNEMRPBins + centBin * fNEMRPBins + rpBin;
    if( fTracksTPCLists->At(offset) ) 
    { 
        // list exists
        TList* list = dynamic_cast<TList*> (fTracksTPCLists->At(offset));
        return list;
    }
    else 
    { 
        // no list for this bin has been created, yet
        TList* list = new TList();
        fTracksTPCLists->AddAt(list, offset);
        return list;
    }
}

//_______________________________________________________________________________
Double_t AliPHOSCorrelations::GetAssocBin(const Double_t& pt) const
{
    //Calculates bin of associated particle pt.
    for(Int_t i=1; i<fAssocBins.GetSize(); i++)
    {
        if(pt>fAssocBins.At(i-1) && pt<fAssocBins.At(i))
            return fAssocBins.At(i) ;
    }

    return fAssocBins.At(fAssocBins.GetSize()-1) ;
}

//_______________________________________________________________________________
void AliPHOSCorrelations::FillTrackEtaPhi() const
{
    // Distribution TPC's tracks by angles.
    for (Int_t i1=0; i1<fTracksTPC->GetEntriesFast(); i1++)
    {
        TLorentzVector * track = (TLorentzVector*)fTracksTPC->At(i1);
        Double_t phiAssoc = track->Phi();
        Double_t etaAssoc = track->Eta();
        Double_t ptAssoc  = track->Pt() ;

        FillHistogram( "track_phieta", phiAssoc, etaAssoc, ptAssoc );
    }
}

//_______________________________________________________________________________
void AliPHOSCorrelations::UpdatePhotonLists()
{
    //Now we either add current events to stack or remove
    //If no photons in current event - no need to add it to mixed

    TList * arrayList = GetCaloPhotonsPHOSList(fVtxBin, fCentBin, fEMRPBin);
    if( fDebug >= 2 )
        AliInfo( Form("fCentBin=%d, fCentNMixed[]=%d",fCentBin,fCentNMixed[fCentBin]) );
    if(fCaloPhotonsPHOS->GetEntriesFast()>0)
    {
        arrayList->AddFirst(fCaloPhotonsPHOS) ;
        fCaloPhotonsPHOS=0x0;
        if(arrayList->GetEntries() > fCentNMixed[fCentBin])
        { 
            // Remove redundant events
            TClonesArray * tmp = static_cast<TClonesArray*>(arrayList->Last()) ;
            arrayList->RemoveLast() ;
            delete tmp; 
        }
    }
}

//_______________________________________________________________________________
void AliPHOSCorrelations::UpdateTrackLists()
{
    //Now we either add current events to stack or remove
    //If no photons in current event - no need to add it to mixed

    TList * arrayList = GetTracksTPCList(fVtxBin, fCentBin, fEMRPBin);

    if( fDebug >= 2 )
        AliInfo( Form("fCentBin=%d, fCentNMixed[]=%d",fCentBin,fCentNMixed[fCentBin]) );
    if(fTracksTPC->GetEntriesFast()>0)
    {
        arrayList->AddFirst(fTracksTPC) ;
        fTracksTPC=0x0;
        if(arrayList->GetEntries() > fCentNMixed[fCentBin])
        { 
            // Remove redundant events
            TClonesArray * tmp = static_cast<TClonesArray*>(arrayList->Last()) ;
            arrayList->RemoveLast() ;
            delete tmp; 
        }
    }
}

//_______________________________________________________________________________
Bool_t AliPHOSCorrelations::SelectESDTrack(const AliESDtrack * t) const
{
    // Estimate if this track can be used for correlation analisys. If all right - return "TRUE".
    Float_t pt = t->Pt();
    if( pt<0.5 || pt>20. ) return kFALSE ;
    if( TMath::Abs( t->Eta() ) > 0.8 ) return kFALSE;
    if( !fESDtrackCuts->AcceptTrack(t) ) return kFALSE ;
    
    return kTRUE ;
}

//_______________________________________________________________________________
Bool_t AliPHOSCorrelations::SelectAODTrack(const AliAODTrack * t) const
{
    // Estimate if this track can be used for correlation analisys. If all right - return "TRUE".
    Float_t pt = t->Pt() ;
    if( pt<0.5 || pt>20.     ) return kFALSE ;
    if( TMath::Abs(t->Eta()) > 0.8 ) return kFALSE ;

    if(fSelectHybridTracks)
    {
        if(!t->IsHybridGlobalConstrainedGlobal())
        {
            if( fDebug > 2 ) printf("\t Reject track, IsHybridGlobalConstrainedGlobal() is FALSE\n ");
            return kFALSE ;
        }
    }
    
    if(fSelectSPDHitTracks) //Not much sense to use with TPC only or Hybrid tracks
    {
        if(!t->HasPointOnITSLayer(0) && !t->HasPointOnITSLayer(1)) 
        {
            if( fDebug > 2 ) printf("\t Reject track, HasPointOnITSLayer(0) is %i and HasPointOnITSLayer(1) is %i\n", t->HasPointOnITSLayer(0), t->HasPointOnITSLayer(1) );
            return kFALSE ; ;
        }
    }

    if(fSelectFractionTPCSharedClusters)
    {
        Double_t frac = Double_t(t->GetTPCnclsS()) / Double_t(t->GetTPCncls());
        if (frac > fCutTPCSharedClustersFraction)
        {
            if( fDebug > 2 ) printf("\t Reject track, shared cluster fraction %f > %f\n",frac, fCutTPCSharedClustersFraction);
            return kFALSE ;
        }
    }

    return kTRUE ;
}

//_______________________________________________________________________________
void AliPHOSCorrelations::LogProgress(int step)
{
    // Fill "step by step" hist
    FillHistogram("hTotSelEvents", step+0.5);
}

//_______________________________________________________________________________
void AliPHOSCorrelations::LogSelection(const int& step, const int& internalRunNumber) const
{
    // the +0.5 is not realy neccisarry, but oh well... -henrik
    FillHistogram("hSelEvents", step+0.5, internalRunNumber-0.5);
}

//_______________________________________________________________________________
Bool_t AliPHOSCorrelations::TestMass(const Double_t& m, const Double_t& pt) const
{
    // If pi0 candidate outside of mass peak, then return false. 
    if (fUseMassWindowParametrisation && fNSigmaWidth)     
    {
        // Parametrization
        FillHistogram("massWindow", MassMeanFunction(pt), MassSigmaFunction(pt)*fNSigmaWidth);
        if ( MassMeanFunction(pt)-MassSigmaFunction(pt)*fNSigmaWidth<=m && m<=MassMeanFunction(pt)+MassSigmaFunction(pt)*fNSigmaWidth )
        {
            FillHistogram("massWindowPass", 1);
            return true;
        }
        else
        {
            FillHistogram("massWindowPass", 2);
            return false;
        }
    }
    else
    {
        if(!fNSigmaWidth)
        {
            AliInfo( Form("fNSigmaWidth equel 0. Class will use default mass window: from %f to %f GeV.", fMassInvMeanMin, fMassInvMeanMax ) ); 
        }

        if(fMassInvMeanMin<=m && m<=fMassInvMeanMax)
        {
            FillHistogram("massWindowPass", 3);
            FillHistogram("massWindow", m, 0.025); // 0.025 just for filling.
            return true;
        }
        else
        {
            FillHistogram("massWindowPass", 4);
            return false;
        }
    }
} 

//_______________________________________________________________________________
Double_t AliPHOSCorrelations::MassMeanFunction(const Double_t &pt) const
{
    // Parametrization mean of mass window
    return ( fMassMean[0]*pt + fMassMean[1] );
}

//_______________________________________________________________________________
Double_t AliPHOSCorrelations::MassSigmaFunction(const Double_t &pt) const
{
    // Parametrization sigma of mass window
    return ( TMath::Sqrt(fMassSigma[0]*fMassSigma[0]/pt + fMassSigma[1]*fMassSigma[1]/pt/pt + fMassSigma[2]*fMassSigma[2]));
}

//_____________________________________________________________________________
void AliPHOSCorrelations::FillHistogram(const char * key, Double_t x) const
{
    //FillHistogram
    TH1 * hist = dynamic_cast<TH1*>(fOutputContainer->FindObject(key)) ;
    if(hist)
        hist->Fill(x) ;
    else
        AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}

//_____________________________________________________________________________
void AliPHOSCorrelations::FillHistogram(const char * key, Double_t x, Double_t y)const
{
    //FillHistogram
    TH1 * th1 = dynamic_cast<TH1*> (fOutputContainer->FindObject(key));
    if(th1)
        th1->Fill(x, y) ;
    else
        AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}

//_____________________________________________________________________________
void AliPHOSCorrelations::FillHistogram(const char * key, Double_t x, Double_t y, Double_t z) const
{
    //Fills 1D histograms with key
    TObject * obj = fOutputContainer->FindObject(key);

    TH2 * th2 = dynamic_cast<TH2*> (obj);
    if(th2) 
    {
        th2->Fill(x, y, z) ;
        return;
    }

    TH3 * th3 = dynamic_cast<TH3*> (obj);
    if(th3) 
    {
        th3->Fill(x, y, z) ;
        return;
    }

    AliError(Form("can not find histogram (of instance TH2) <%s> ",key)) ;
}

//_____________________________________________________________________________
void AliPHOSCorrelations::FillHistogram(const char * key, Double_t x,Double_t y, Double_t z, Double_t w) const
{
    //Fills 1D histograms with key
    TObject * obj = fOutputContainer->FindObject(key);

    TH3 * th3 = dynamic_cast<TH3*> (obj);
    if(th3) 
    {
        th3->Fill(x, y, z, w) ;
        return;
    }

    AliError(Form("can not find histogram (of instance TH3) <%s> ",key)) ;
}

//_____________________________________________________________________________
void AliPHOSCorrelations::SetGeometry()
{
    // Initialize the PHOS geometry
    //Init geometry
    if(!fPHOSGeo)
    {
        AliOADBContainer geomContainer("phosGeo");
        geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
        TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");
        fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP") ;
        for(Int_t mod=0; mod<5; mod++) 
        {
            if(!matrixes->At(mod)) 
            {
                if( fDebug )
                AliInfo(Form("No PHOS Matrix for mod:%d, geo=%p\n", mod, fPHOSGeo));
                continue;
            }
            else 
            {
                fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
                if( fDebug >1 )
                    AliInfo(Form("Adding PHOS Matrix for mod:%d, geo=%p\n", mod, fPHOSGeo));
            }
        }
    } 
}

//_____________________________________________________________________________
Double_t AliPHOSCorrelations::GetEfficiency(const Double_t& pt) const 
{
    Double_t x = pt;
    //Efficiency for Both2core only!
    if (!fUseEfficiency)
        return 1.;

    Double_t e =1.;
     // From 0 to 5 - 11h for different centrality.
     /*0: 0-5%
    1: 5-10%
    2: 10-20%
    3: 20-40%
    4: 40-60%
    5: 60-80%
    6: 0-20%
    7: 0-10%*/
    Double_t par0[9] = { -798863,       339.714,    6407.1,     -457.778,   1283.65,    -117.075,   -19.3764,   0,          0       };
    Double_t par1[9] = { -799344,       -1852.1,    3326.29,    -384.229,   504.046,    562.608,    130.518,    0,          0       };
    Double_t par2[9] = { -858904,       -1923.28,   5350.74,    -568.946,   945.497,    419.647,    101.911,    0,          0       };
    Double_t par3[9] = { -795652,       -1495.97,   2926.46,    -357.804,   478.961,    551.127,    128.86,     0,          0       };
    Double_t par4[9] = { -891951,       279626,     -123110,    -5464.75,   27470.8,    283264,     15355.1,    192762,     44828.6 };
    Double_t par5[9] = { -1.1094e+06,   -986.915,   2127.71,    -268.908,   375.594,    380.791,    89.4053,    0,          0       };
    // Double_t par6[7] = {4.86106e+09, 4.47013e+08, -1.48079e+09, 1.47233e+08, -2.62356e+08, -1.00639e+08, -2.45629e+07, 0, 0};
    // Double_t par7[7] = {-1.36243e+06, -26011.1, 135838, -12161.3, 24956.8, 4985.4, 1285.57, 0, 0};

     // 8 bin for pPb13 and 0-100%
    Double_t par8[9] = { 6.87095e+06,   8.36553e+06,    -3.29572e+06,   2.18688e+06,    -739490,    521666,     106661,     0,  0   };
         
    Double_t* pFitPoint;

    // TODO: fix functions value below 1 GeV/c.
    if(x < 1.) x = 1.;

    if(GetPeriod().Contains("11h"))
    {
        if (fCentrality <= 5)                       pFitPoint = &par0[0];
        if (fCentrality > 5 && fCentrality <= 10)   pFitPoint = &par1[0];
        if (fCentrality > 10 && fCentrality <= 20)  pFitPoint = &par2[0];
        if (fCentrality > 20 && fCentrality <= 40)  pFitPoint = &par3[0];
        if (fCentrality > 40 && fCentrality <= 60)  pFitPoint = &par4[0];
        if (fCentrality > 60)                       pFitPoint = &par5[0];

        Double_t pFit[9];
        for (int i = 0; i < 9; ++i)
         {
             pFit[i] = *(pFitPoint+i);
         }

        if (fCentrality > 40 && fCentrality <= 60)
            e = TMath::Exp(-(((((1.+(pFit[1]*x))+(pFit[2]*(x*x)))+(pFit[5]*(x*(x*x))))+(pFit[7]*(x*(x*(x*x)))))/((((pFit[3]*x)+(pFit[4]*(x*x)))+(pFit[6]*(x*(x*x))))+(pFit[8]*(x*(x*(x*x))))))) ;
        else
            e = TMath::Exp(-((((1.+(pFit[1]*x))+(pFit[2]*(x*x)))+(pFit[5]*(x*(x*x))))/(((pFit[3]*x)+(pFit[4]*(x*x)))+(pFit[6]*(x*(x*x)))))) ;
    }
    else
    if(GetPeriod().Contains("13")) 
    {
        pFitPoint = &par8[0];
        Double_t pFit[9];
        for( int i = 0; i < 9; i++ )
         {
             pFit[i] = *(pFitPoint+i);
         }

        e = TMath::Exp(-((((pFit[0]+(pFit[1]*x))+(pFit[2]*(x*x)))+(pFit[5]*(x*(x*x))))/(((1.+(pFit[3]*x))+(pFit[4]*(x*x)))+(pFit[6]*(x*(x*x)))))) ;
    }
    else
    {
        // No case
        AliWarning(Form("No efficiensy choise. Return 1"));
        e = 1.;
    }

    return e;
}

//_____________________________________________________________________________
Int_t AliPHOSCorrelations::GetModCase(const Int_t &mod1, const Int_t &mod2) const 
{
    // Return modules pair namber.
    if(mod1 == mod2)
    {
        if(mod1 == 1) return 1;
        if(mod1 == 2) return 2;
        if(mod1 == 3) return 3;
    }
    else
    {
        if(mod1 == 1 || mod2 == 1)
            if(mod1 == 2 || mod2 == 2)
                return 12;

        if(mod1 == 1 || mod2 == 1)
            if(mod1 == 3 || mod2 == 3)
                return 13;
        if(mod1 == 2 || mod2 == 2)
            if(mod1 == 3 || mod2 == 3)
                return 23;
    }

    AliError(Form("No choise for mod1 = %i, mod2 = %i", mod1, mod2));
    return 1;
}

//_____________________________________________________________________________
void AliPHOSCorrelations::TestPi0ME(const Int_t& ipid, const TLorentzVector& p12, const Int_t& modCase)
{
    Double_t phiTrigger = p12.Phi() ;
    Double_t etaTrigger = p12.Eta() ;
    Double_t pt         = p12.Pt() ;

    if ( GetMEExists(ipid) )
    {
        if ( pt >= GetMEPt(ipid) )
        {
            SetMEPt(ipid,pt);
            SetMEPhi(ipid, phiTrigger);
            SetMEEta(ipid, etaTrigger);
            SetMEModCase(ipid, modCase);
        }
    }
    else
    {
        SetMEPt(ipid,pt);
        SetMEPhi(ipid, phiTrigger);
        SetMEEta(ipid, etaTrigger);
        SetMEModCase(ipid, modCase);
        SetMEExists(ipid);
    }
}

//_____________________________________________________________________________
void AliPHOSCorrelations::ZeroingVariables()
{
    // Set Phi, Eta, pT, modNumber andtrigger variable of moust energetic trigger particle to zero.
    for (int i = 0; i < 4; ++i)
    {
        fMEExists[i] = false;
        fMEPhi[i] = fMEEta[i] = fMEPt[i] = -99;
        fMEModCase[i] = 1;
    }
}

//_____________________________________________________________________________
void AliPHOSCorrelations::SetMassMeanParametrs(const Double_t par[2])                                      
{ 
    for (int i = 0; i < 2; ++i)
    {
        fMassMean[i] = par[i] ;  
    }                  
} 

//_____________________________________________________________________________
void AliPHOSCorrelations::SetMassSigmaParametrs(const Double_t par[3])                                   
{ 
    for (int i = 0; i < 3; ++i)
    {
        fMassSigma[i] = par[i] ;    
    }                  
}

//_____________________________________________________________________________
void AliPHOSCorrelations::FillEventBiningProperties() const
{
    // Fill fCentBin, fEMRPBin, fVtxBin.
    if(isREALEvent || fMBEvent)
    {
        FillHistogram( "hCentralityBining", fCentBin) ;
        FillHistogram( "phiRPflatBining",   fEMRPBin) ;
        FillHistogram( "hVertexZBining",    fVtxBin)  ;
        if(isREALEvent) 
        {
            FillHistogram( "hCentralityBiningTrigger", fCentBin) ;
            FillHistogram( "phiRPflatBiningTrigger",   fEMRPBin) ;
            FillHistogram( "hVertexZBiningTrigger",    fVtxBin) ;
        }
        if(fMBEvent)   
        {
            FillHistogram( "hCentralityBiningMB", fCentBin) ;
            FillHistogram( "phiRPflatBiningMB",   fEMRPBin) ;
            FillHistogram( "hVertexZBiningMB",    fVtxBin) ;
        }
    }
}

void AliPHOSCorrelations::SetEventMixingVtxBinning(const Int_t nBins) 
{ 
    fNVtxZBins = nBins ; 
    SetVertexBinning();
}

//_____________________________________________________________________________
void AliPHOSCorrelations::ShowTaskInfo()
{
    // Show all info about task settings.
    AliInfo("//________________________________________________");
    AliInfo(Form("Period: %s", fPeriod.Data()));
    
    AliInfo("Bining:");
    AliInfo(Form("Number Of Centrality Bins = %i", GetNumberOfCentralityBins()));
    AliInfo(Form("fNVtxZBins = %i", fNVtxZBins));
    AliInfo(Form("fNEMRPBins = %i", fNEMRPBins));

    AliInfo(Form("fCentralityEstimator = %s", fCentralityEstimator.Data()));
    AliInfo(Form("fEventPlaneMethod = %s", fEventPlaneMethod.Data()));
    
    AliInfo("Global event cuts:");
    if (GetPeriod().Contains("11h")) AliInfo(Form("fTriggerSelectionPbPb = %i (Info: kCentAndSCentCS=0; kCent=1; kSCent=2)", fTriggerSelectionPbPb ));
    AliInfo(Form("fCentralityLowLimit = %f", fCentralityLowLimit));
    AliInfo(Form("fCentralityHightLimit = %f", fCentralityHightLimit));
    AliInfo(Form("fMaxAbsVertexZ = %f", fMaxAbsVertexZ));
    
    AliInfo("PHOS claster cuts:");
    AliInfo(Form("fMinClusterEnergy = %f", fMinClusterEnergy));
    AliInfo(Form("fMinBCDistance = %f", fMinBCDistance));
    AliInfo(Form("fMinNCells = %i", fMinNCells));
    AliInfo(Form("fMinM02 = %f", fMinM02));
    AliInfo(Form("fTOFCutEnabled = %i", fTOFCutEnabled));
    AliInfo(Form("fTOFCut = %f", fTOFCut));
    AliInfo(Form("fUseMassWindowParametrisation = %i", fUseMassWindowParametrisation));
    AliInfo(Form("fMassInvMeanMin = %f", fMassInvMeanMin));
    AliInfo(Form("fMassInvMeanMax = %f", fMassInvMeanMax));
    AliInfo(Form("fNSigmaWidth = %f", fNSigmaWidth));
    AliInfo(Form("fUseEfficiency = %i", fUseEfficiency));

    AliInfo("Track cuts:");
    AliInfo(Form("fSelectHybridTracks = %i", fSelectHybridTracks));
    AliInfo(Form("fSelectSPDHitTracks = %i", fSelectSPDHitTracks));
    AliInfo(Form("fSelectFractionTPCSharedClusters = %i", fSelectFractionTPCSharedClusters));
    AliInfo(Form("fCutTPCSharedClustersFraction = %f", fCutTPCSharedClustersFraction));

    AliInfo("Mass cuts:");
    if (fUseMassWindowParametrisation) AliInfo(Form("Parametrization: Maen [%f, %f], Sigma[%f, %f, %f]", fMassMean[0], fMassMean[1], fMassSigma[0], fMassSigma[1], fMassSigma[2]));
    else  AliInfo( Form("Class will use default mass window: from %f to %f GeV.", fMassInvMeanMin, fMassInvMeanMax ) ); 
    AliInfo(Form("fNSigmaWidth = %f", fNSigmaWidth));
}
