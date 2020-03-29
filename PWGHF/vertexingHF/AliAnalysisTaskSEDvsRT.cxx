/**************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Class AliAnalysisTaskSEDvsRT
// AliAnalysisTaskSE for charmed hadrons vs. transverse activity analysis
// Authors: Jeremy Wilkinson,
/////////////////////////////////////////////////////////////


#include <TClonesArray.h>
#include <TCanvas.h>
#include <TList.h>
#include <TString.h>
#include <TDatabasePDG.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TProfile.h>
#include "AliAnalysisManager.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDvsRT.h"
#include "AliNormalizationCounter.h"
#include "AliVertexingHFUtils.h"
#include "AliAODVZERO.h"
#include "AliESDUtils.h"


//________________________________________________________________________
AliAnalysisTaskSEDvsRT::AliAnalysisTaskSEDvsRT():
   AliAnalysisTaskSE(),
   fOutput(0),
   fListCuts(0),
   fOutputCounters(0),
   fListQAhists(0),
   fUpmasslimit(1.965),
   fLowmasslimit(1.765),
   fNMassBins(200),
   fRDCutsAnalysis(0),
   fHistNEvents(0),
   fGlobalRT(0),
   fHistPtLead(0),
   fRTvsZvtxvsMult(0),
   fCounter(0),
   fReadMC(kFALSE),
   fMCOption(0),
   fUseBit(kTRUE),
   fAODProtection(1),
   fPdgSpecies(411),
   fLctoV0(kFALSE),
   fisPPbData(kFALSE),
   fEtaCut(1.5),
   fLeadMin(6.0),
   fAveMultiInTrans(4.939),
   fPhiLeading(0),
   fPtvsMassvsRTToward(0),
   fPtvsMassvsRTAway(0),
   fPtvsMassvsRTTrans(0),
   fUseNsparse(kFALSE),
   fOutNsparse(0)

   {
      ///default constructor
      for (Int_t i = 0; i < 18; i++) {fTrackFilter[i] = 0;}
   }

AliAnalysisTaskSEDvsRT::AliAnalysisTaskSEDvsRT(const char *name, Int_t pdgSpecies, AliRDHFCuts *cuts):
   AliAnalysisTaskSE(name),
   fOutput(0),
   fListCuts(0),
   fOutputCounters(0),
   fListQAhists(0),
   fUpmasslimit(1.965),
   fLowmasslimit(1.765),
   fNMassBins(200),
   fRDCutsAnalysis(cuts),
   fHistNEvents(0),
   fGlobalRT(0),
   fHistPtLead(0),
   fRTvsZvtxvsMult(0),
   fCounter(0),
   fReadMC(kFALSE),
   fMCOption(0),
   fUseBit(kTRUE),
   fAODProtection(1),
   fPdgSpecies(pdgSpecies),
   fLctoV0(kFALSE),
   fisPPbData(kFALSE),
   fEtaCut(1.5),
   fLeadMin(6.0),
   fAveMultiInTrans(4.939),
   fPhiLeading(0),
   fPtvsMassvsRTToward(0),
   fPtvsMassvsRTAway(0),
   fPtvsMassvsRTTrans(0),
   fUseNsparse(kFALSE),
   fOutNsparse(0)
   {
      ///Default constructor
      if (fPdgSpecies == 413) {
         fNMassBins = 200;
         SetMassLimits(0.12,0.2);
      }else if(fPdgSpecies == 431) {
        Double_t MassBinSize  = 0.002;
        Int_t nInvMassBins = (Int_t)(0.7/MassBinSize+0.5);
        Double_t massDs  = TDatabasePDG::Instance()->GetParticle(431)->Mass();
        Double_t minMass = massDs-0.5*nInvMassBins*MassBinSize;
        Double_t maxMass = massDs+0.5*nInvMassBins*MassBinSize;
        SetMassLimits(minMass,maxMass);
        SetNMassBins(nInvMassBins); 
      }else if(fPdgSpecies == 4122) {
        Double_t massLc  = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
        Int_t nInvMassBins = 500;
        Double_t minMass = massLc-0.180;
        Double_t maxMass = massLc+0.180;
        SetMassLimits(minMass,maxMass);
        SetNMassBins(nInvMassBins);
      }
      else{
        fNMassBins=200;
        SetMassLimits(fPdgSpecies,0.1);
      }

      DefineOutput(1,TList::Class()); //Output slot 1: fOutput TList container
      DefineOutput(2,TList::Class()); //Output slot 2: fListCuts
      DefineOutput(3,TList::Class()); //Output slot 3: normalisation counters
      DefineOutput(4,TList::Class()); // Output slot 4: list for QA plots
      for (Int_t i = 0; i < 18; i++) {fTrackFilter[i] = 0;}
}   
   
//________________________________________________________________________
AliAnalysisTaskSEDvsRT::~AliAnalysisTaskSEDvsRT()
{
   //
   /// Standard destructor
   //
   if (fOutput) {
      delete fOutput;
      fOutput = 0x0;
   }
   if (fListCuts) {
      delete fListCuts;
      fListCuts = 0x0;
   }
  
   if (fOutputCounters) {
      delete fOutputCounters;
      fOutputCounters = 0x0;
   }
   if (fListQAhists) {
      delete fListQAhists;
      fListQAhists = 0x0;
   }
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsRT::SetMassLimits(Double_t lowlimit, Double_t uplimit){
   /// set limits for inv. mass axis between two specific values
   if (uplimit > lowlimit) {
      fLowmasslimit = lowlimit;
      fUpmasslimit = uplimit;
   } else {
      AliError("Wrong mass limits: upper limit should be larger than lower limit");
   }
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsRT::SetMassLimits(Int_t pdg, Double_t range) {
   /// set limits for inv. mass axis for given particle species
   Double_t mass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdg))->Mass();
   SetMassLimits(mass-range,mass+range);
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsRT::Init(){
   //
   /// Initialisation
   //
   
   printf("AliAnalysisTaskSEDvsRT_0::Init() \n");
   
   fListCuts = new TList();
   fListCuts->SetOwner();
   fListCuts->SetName("CutsList");
   
   if(fPdgSpecies==411){
    AliRDHFCutsDplustoKpipi* copycut=new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDplus");
    fListCuts->Add(copycut);
  }else if(fPdgSpecies==421){
    AliRDHFCutsD0toKpi* copycut=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDzero");
    fListCuts->Add(copycut);
  }else if(fPdgSpecies==413){
    AliRDHFCutsDStartoKpipi* copycut=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDStar");
    fListCuts->Add(copycut);
  }else if(fPdgSpecies==431){
    AliRDHFCutsDstoKKpi* copycut=new AliRDHFCutsDstoKKpi(*(static_cast<AliRDHFCutsDstoKKpi*>(fRDCutsAnalysis)));
    copycut->SetName("AnalysisCutsDs");
    fListCuts->Add(copycut);
  }else if(fPdgSpecies==4122){
    if(fLctoV0){
      AliRDHFCutsLctoV0* copycut=new AliRDHFCutsLctoV0(*(static_cast<AliRDHFCutsLctoV0*>(fRDCutsAnalysis)));
      copycut->SetName("AnalysisCutsLc2pK0S");
      fListCuts->Add(copycut);
      }else{
      AliRDHFCutsLctopKpi *copycut=new AliRDHFCutsLctopKpi(*(static_cast<AliRDHFCutsLctopKpi*>(fRDCutsAnalysis)));
      copycut->SetName("LctopKpiProdCuts");
      fListCuts->Add(copycut);
      }
  }
   PostData(2,fListCuts);
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsRT::UserCreateOutputObjects()
{
   /// Create output containers
   
   if (fDebug > 1) printf("AliAnalysisTaskSEDvsRT::UserCreateOutputObjects() \n");
   
   // TList for output
   fOutput = new TList();
   fOutput->SetOwner();
   fOutput->SetName("OutputHistos");

   fCounter = new AliNormalizationCounter("NormCounterCorrMult");
   fCounter->SetStudyMultiplicity(kTRUE,1.);
   fCounter->Init();
   fOutputCounters = new TList();
   fOutputCounters->SetOwner();
   fOutputCounters->SetName("OutputCounters");
   fOutputCounters->Add(fCounter);

   Double_t ptmin = 0.; Double_t ptmax = 24.;
   Int_t nPTbins = 48;

   Double_t firstRTbin = 0.; Double_t lastRTbin = 10.;
   Int_t nRTbins = 20;

   fPtvsMassvsRTToward = new TH3D("hPtvsMassvsRTToward", "D candidates in toward region: p_{T} vs mass vs R_{T};R_{T};Mass [GeV/c^{2}];p_{T} [GeV/c]",nRTbins,firstRTbin,lastRTbin, fNMassBins,fLowmasslimit,fUpmasslimit,nPTbins, ptmin, ptmax);
   fPtvsMassvsRTAway = new TH3D("hPtvsMassvsRTAway", "D candidates in away region: p_{T} vs mass vs R_{T};R_{T};Mass [GeV/c^{2}];p_{T} [GeV/c]",nRTbins,firstRTbin,lastRTbin, fNMassBins,fLowmasslimit,fUpmasslimit,nPTbins, ptmin, ptmax);

   fPtvsMassvsRTTrans = new TH3D("hPtvsMassvsRTTrans", "D candidates in transverse region: p_{T} vs mass vs R_{T};R_{T};Mass [GeV/c^{2}];p_{T} [GeV/c}",nRTbins,firstRTbin,lastRTbin, fNMassBins, fLowmasslimit, fUpmasslimit, nPTbins, ptmin, ptmax);
  
   fOutput->Add(fPtvsMassvsRTToward);
   fOutput->Add(fPtvsMassvsRTAway);
   fOutput->Add(fPtvsMassvsRTTrans);
   
   
   fHistNEvents = new TH1F("fHistNEvents", "number of events ",12,-0.5,11.5);
   fHistNEvents->GetXaxis()->SetBinLabel(1,"nEvents total");
   fHistNEvents->GetXaxis()->SetBinLabel(2,"nEvents with Z vertex");
   fHistNEvents->GetXaxis()->SetBinLabel(3,"nEvents selected");
   fHistNEvents->GetXaxis()->SetBinLabel(4,"Rejected due to trigger");
   fHistNEvents->GetXaxis()->SetBinLabel(5,"Rejected due to phys sel");
   fHistNEvents->GetXaxis()->SetBinLabel(6,"Rejected due to vertex cuts");
   fHistNEvents->GetXaxis()->SetBinLabel(7,"Rejected due to pileup");
   fHistNEvents->GetXaxis()->SetBinLabel(8,"Total no. of candidate");
   fHistNEvents->GetXaxis()->SetBinLabel(9,"no. of cand wo bitmask");
   fHistNEvents->GetXaxis()->SetBinLabel(10,"D after cuts (No PID)");
   fHistNEvents->GetXaxis()->SetBinLabel(11,"D after cuts + PID)");
   fHistNEvents->GetXaxis()->SetBinLabel(12,"nEvents with RT calc");
   
   fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);  
   fHistNEvents->Sumw2();
   fHistNEvents->SetMinimum(0);
   fOutput->Add(fHistNEvents);
   
   
   //THnSparse notation
   //axes: pt, mass, rtval, candidate phi
   Double_t phimin = 0.; Double_t phimax = TMath::TwoPi();
   Int_t nPhibins  = 200;
   
   static const Int_t sparsendim = 4;
   Double_t sparsemin[sparsendim] = {ptmin, fLowmasslimit, firstRTbin, phimin};
   Double_t sparsemax[sparsendim] = {ptmax, fUpmasslimit,  lastRTbin,  phimax};
   Int_t sparsenbins[sparsendim] = {nPTbins, fNMassBins, nRTbins, nPhibins};
   
   
   if(fUseNsparse) { // define NSparse object if set to use
   fOutNsparse = new THnSparseD("hNsparse", ";pT (GeV/c);Mass (MeV/c^2);RT; #Delta#phi(cand)", sparsendim, sparsenbins,  sparsemin, sparsemax);
   fOutput->Add(fOutNsparse);
   }
   
   ///settings for track filter used in RT determination
  AliESDtrackCuts* esdTrackCutsRun2[18] = {0};
	
  for ( int iTc = 0 ; iTc < 18 ; iTc++ )
  {
	// standar parameters ------------------- //
    double maxdcaz = 2.;  
    double minratiocrossrowstpcover = 0.8;
    double maxfraclusterstpcshared = 0.4; 
    double maxchi2perclustertpc = 4.0; 
    double maxchi2perclusterits = 36.; 
    double geowidth = 3.;
    double geolenght = 130.;
    double maxchi2tpcglobal = 36.;
	// ------------------------------------- //
	
	// variations of the track cuts -------- //
    if ( iTc == 1) maxdcaz = 1.0;
    if ( iTc == 2) maxdcaz = 5.0;
    if ( iTc == 5) minratiocrossrowstpcover = 0.7;
    if ( iTc == 6) minratiocrossrowstpcover = 0.9;
    if ( iTc == 7) maxfraclusterstpcshared = 0.2;
    if ( iTc == 8) maxfraclusterstpcshared = 1.0;
    if ( iTc == 9) maxchi2perclustertpc = 3.0;
    if ( iTc == 10) maxchi2perclustertpc = 5.0;
    if ( iTc == 11) maxchi2perclusterits = 25.0;
    if ( iTc == 12) maxchi2perclusterits = 49.0;
    if ( iTc == 14) geowidth = 2.0;
    if ( iTc == 15) geowidth = 4.0;
    if ( iTc == 16) geolenght = 120.0;
    if ( iTc == 17) geolenght = 140.0;
	
	// variations of the track cuts -------- //
    
    
    fTrackFilter[iTc] = new AliAnalysisFilter(Form("fTrackFilter%d",iTc));
    if (iTc != 0) esdTrackCutsRun2[iTc] = new AliESDtrackCuts(Form("esdTrackCutsRun2%d",iTc));
    
    // TPC
    if (iTc != 0) { //variations using ITS-TPC
    esdTrackCutsRun2[iTc]->SetCutGeoNcrNcl(geowidth,geolenght,1.5,0.85,0.7);
    esdTrackCutsRun2[iTc]->SetRequireTPCRefit(kTRUE);
    esdTrackCutsRun2[iTc]->SetMinRatioCrossedRowsOverFindableClustersTPC(minratiocrossrowstpcover); 
    esdTrackCutsRun2[iTc]->SetMaxChi2PerClusterTPC(maxchi2perclustertpc);
    esdTrackCutsRun2[iTc]->SetMaxFractionSharedTPCClusters(maxfraclusterstpcshared); 
    //esdTrackCutsRun2[iTc]->SetMaxChi2TPCConstrainedGlobal(maxchi2tpcglobal); TODO VZ: check this cut   
    // ITS
    esdTrackCutsRun2[iTc]->SetRequireITSRefit(kTRUE);
    if ( iTc != 13 )
    esdTrackCutsRun2[iTc]->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    
    esdTrackCutsRun2[iTc]->SetMaxChi2PerClusterITS(maxchi2perclusterits);

    // primary selection
    esdTrackCutsRun2[iTc]->SetDCAToVertex2D(kFALSE);
    esdTrackCutsRun2[iTc]->SetRequireSigmaToVertex(kFALSE);
    esdTrackCutsRun2[iTc]->SetMaxDCAToVertexZ(maxdcaz); 
    esdTrackCutsRun2[iTc]->SetAcceptKinkDaughters(kFALSE);
    
    if ( iTc == 3 )
      // esdTrackCutsRun2[iTc]->SetMaxDCAToVertexXYPtDep("4*(0.0026+0.0050/pt^1.01)");
    	esdTrackCutsRun2[iTc]->SetMaxDCAToVertexXYPtDep("6.5*(0.0026+0.0050/pt^1.01)");
    else if ( iTc == 4 )
      // esdTrackCutsRun2[iTc]->SetMaxDCAToVertexXYPtDep("10*(0.0026+0.0050/pt^1.01)");
    	esdTrackCutsRun2[iTc]->SetMaxDCAToVertexXYPtDep("7.5*(0.0026+0.0050/pt^1.01)");
    else
      esdTrackCutsRun2[iTc]->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); // (7*(------))
    }
    else { //Default: TPC-only track filter
      esdTrackCutsRun2[iTc] = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    }
       
  
  
    fTrackFilter[iTc]->AddCuts(esdTrackCutsRun2[iTc]);
  }
   
   
   fListQAhists = new TList();
   fListQAhists->SetOwner();
   fListQAhists->SetName("QAhists");
   
   fGlobalRT = new TH1F("fGlobalRT","RT for all events;R_{T};Entries",100,0,10);
   fHistPtLead = new TH1F("fHistPtLead","pT distribution of leading track;p_{T} (GeV/c);Entries",100,0,100);
   fRTvsZvtxvsMult = new TH3F("fRTvsZvtxvsMult","RT vs Zvtx vs mult;R_{T};Z_{vtx} (cm);N_{trk}",100,0,10,200,-10,10,200,0,200);
   
   
   fListQAhists->Add(fGlobalRT);
   fListQAhists->Add(fHistPtLead);
   fListQAhists->Add(fRTvsZvtxvsMult); 
   
   PostData(1,fOutput);
   PostData(3,fOutputCounters);
   PostData(4,fListQAhists);
   return;   
}

//________________________________________________________________________
void AliAnalysisTaskSEDvsRT::UserExec(Option_t */*option*/)
{
   /// Execute analysis for current event:
   /// heavy flavour candidates as function of RT
   
   AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
   AliESDEvent *esd = dynamic_cast<AliESDEvent*> (InputEvent());
   //  AliAODTracklets* tracklets = aod->GetTracklets();
   //Int_t ntracklets = tracklets->GetNumberOfTracklets();
   if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      return;
    }
  }
   
    
  TClonesArray *arrayCand = 0;
  TString arrayName="";
  UInt_t pdgDau[3];
  Int_t nDau=0;
  Int_t selbit=0;
  if(fPdgSpecies==411){
    arrayName="Charm3Prong";
    pdgDau[0]=211; pdgDau[1]=321; pdgDau[2]=211; 
    nDau=3;
    selbit=AliRDHFCuts::kDplusCuts;
  }else if(fPdgSpecies==421){
    arrayName="D0toKpi";
    pdgDau[0]=211; pdgDau[1]=321; pdgDau[2]=0;
    nDau=2;
    selbit=AliRDHFCuts::kD0toKpiCuts;
  }else if(fPdgSpecies==413){
    arrayName="Dstar";
    pdgDau[0]=321; pdgDau[1]=211; pdgDau[2]=0; // Quoting here D0 daughters (D* ones on another variable later)
    nDau=2;
    selbit=AliRDHFCuts::kDstarCuts;
  }else if(fPdgSpecies==431){
    arrayName="Charm3Prong";
    pdgDau[0]=321; pdgDau[1]=321; pdgDau[2]=211;
    nDau=3;
    selbit=AliRDHFCuts::kDsCuts;
  }else if(fPdgSpecies==4122){
    if(fLctoV0){
    arrayName="CascadesHF";
    pdgDau[0]=211; pdgDau[1]=211; pdgDau[2]=0; // Quoting here K0S daughters (Lc ones on another variable later)
    nDau=2;
    selbit=AliRDHFCuts::kLctoV0Cuts;
    }else{
    arrayName="Charm3Prong";
    pdgDau[0]=2212; pdgDau[1]=321; pdgDau[2]=211;
    nDau=3;
    selbit=AliRDHFCuts::kLcCuts;
    }
  }

//Arrays of daughters for inv mass calculation
    UInt_t pdgdaughtersD0[2] = {211,321}; //pi, K
    UInt_t pdgdaughtersD0bar[2] = {321,211}; //K, pi
    UInt_t pdgDsKKpi[3] = {321,321,211};
    UInt_t pdgDspiKK[3] = {311, 321,321};  
  
  
  
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      arrayCand=(TClonesArray*)aodFromExt->GetList()->FindObject(arrayName.Data());
    }
  } else if(aod) {
    arrayCand=(TClonesArray*)aod->GetList()->FindObject(arrayName.Data());
  }

  if(!aod || !arrayCand) {
    printf("AliAnalysisTaskSEDvsMultiplicity::UserExec: Charm3Prong branch not found!\n");
    return;
  }

  if(fisPPbData && fReadMC){
    Int_t runnumber = aod->GetRunNumber();
    if(aod->GetTriggerMask()==0 && 
       (runnumber>=195344 && runnumber<=195677)){
      AliDebug(3,"Event rejected because of null trigger mask");
      return;
    }
  }
  fHistNEvents->Fill(0); // count event


  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex()||TMath::Abs(aod->GetMagneticField())<0.001) return;
   switch (fRDCutsAnalysis->GetWhyRejection()) {
      case 5: fHistNEvents->Fill(3); break;
      case 7: fHistNEvents->Fill(4); break;
      case 6: fHistNEvents->Fill(5); break;
      case 1: fHistNEvents->Fill(6); break;
      
   }  
    Bool_t isEvSel = fRDCutsAnalysis->IsEventSelected(aod); //?
  if(!isEvSel)return;
  fHistNEvents->Fill(2);
  
  
  
  ///!----- l.797-833 for multiplicity estimation; RT determination goes here
  
  Double_t rtval = CalculateRTVal(aod);
  if (rtval < 0.) return;  // event rejected during RT calculation
  fCounter->StoreEvent(aod,fRDCutsAnalysis,fReadMC,rtval);
  fHistNEvents->Fill(12);
  fGlobalRT->Fill(rtval);
  
  AliAODTracklets* tracklets = aod->GetTracklets();
  Int_t nTr = tracklets->GetNumberOfTracklets();
  Int_t countMult = 0;
  Double_t theta, eta;
  for (Int_t iTr = 0; iTr<nTr; iTr++) { //count Ntr in eta < 1 for QA
   theta = tracklets->GetTheta(iTr);
   eta = -TMath::Log(TMath::Tan(theta/2.));
   if (eta > -1.0 && eta < 1.0) countMult++;
  }
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  fRTvsZvtxvsMult->Fill(rtval,vtx1->GetZ(),countMult);  /// QA plot for RT vs Zvtx vs Ntrk
  Double_t weight = 1.; //dummy weight for filling (needed later?)
  //!----l.839-888: multiplicity correction
  



//l.948-1094: MC multiplicity counting/reweighting  
  
  Int_t nCand = arrayCand->GetEntriesFast();
  Int_t nSelectedNoPID=0,nSelectedPID=0,nSelectedInMassPeak=0;
  Double_t mD0PDG    = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  Double_t mDplusPDG = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  Double_t mDstarPDG = TDatabasePDG::Instance()->GetParticle(413)->Mass();
  Double_t mDsPDG    = TDatabasePDG::Instance()->GetParticle(431)->Mass();
  Double_t mLcPDG    = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  
  // PDG of daughters for D*
  UInt_t pdgDgDStartoD0pi[2] = {421, 211};
  
  // PDG of daughters for Lc2pK0
  UInt_t pdgDgLctopK0S[2] = {2212, 310};
 //Load MC info if MC  
TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;

  Double_t nchWeight=1.0;

  // load MC particles
  if(fReadMC){
     
    //arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      printf("AliAnalysisTaskSEDvsMultiplicity::UserExec: MC particles branch not found!\n");
      return;
    }  
    // load MC header
   // mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEDvsMultiplicity::UserExec: MC header branch not found!\n");
      return;
    }
 } 
   // omitting "aveMult" l.1110
  
  //TODO:  Loop on candidates, perform selection, !!determine phi of candidate wrt leading particle, fill corresponding histo (toward/away)
  
  for (Int_t iCand = 0; iCand < nCand; iCand++) {  //Loop over candidates
     
     AliAODRecoDecayHF *d =(AliAODRecoDecayHF*)arrayCand->UncheckedAt(iCand);
     AliAODRecoCascadeHF *dCascade = NULL;
     if (fPdgSpecies == 413 || (fPdgSpecies==4122 && fLctoV0)) dCascade = (AliAODRecoCascadeHF*)d;
     
     fHistNEvents->Fill(7);
     
     if(fPdgSpecies == 4122) { //V0 selection for LctopK0
        if(fLctoV0) {
           AliAODv0 *v0part = (AliAODv0*)dCascade->Getv0();
           Bool_t onFlyV0 = v0part->GetOnFlyStatus();
           if (onFlyV0)  {fHistNEvents->Fill(8); continue;}
        } else if (!fLctoV0 && fUseBit && !d->HasSelectionBit(selbit)) { // check selection bit for LcpKpi
           fHistNEvents->Fill(8); continue;
        }
     }
    
     if (fPdgSpecies != 4122 && fUseBit && !d->HasSelectionBit(selbit)) { //check selection bit for all others
         fHistNEvents->Fill(8); continue;
     }

     Double_t ptCand = d->Pt();
     Double_t rapid = d->Y(fPdgSpecies);
     Bool_t isInFidAcc = fRDCutsAnalysis->IsInFiducialAcceptance(ptCand,rapid);
     if (!isInFidAcc) continue;
     TClonesArray *arrayMC =0; 
     Int_t labD = -1;
     if(fReadMC) {
      if(fPdgSpecies==413){
         labD = dCascade->MatchToMC(fPdgSpecies,421,(Int_t*)pdgDgDStartoD0pi,(Int_t*)pdgDau,arrayMC);
      } else if(fPdgSpecies==4122 && fLctoV0){
         labD = dCascade->MatchToMC(fPdgSpecies,pdgDgLctopK0S[1],(Int_t*)pdgDgLctopK0S,(Int_t*)pdgDau,arrayMC,kTRUE);
      } else {
         labD = d->MatchToMC(fPdgSpecies,arrayMC,nDau,(Int_t*)pdgDau);
      }
      // FillMCMassHistos(arrayMC,labD, countMult,nchWeight);  ////!TODO be implemented
     }
     
     Int_t passAllCuts = 0, passTopolCuts = 0;
     
     if (fPdgSpecies == 4122) {  //special cuts cases for Lc
        if (fLctoV0) { //LctopK0
           passAllCuts = (((fRDCutsAnalysis->IsSelected(d, AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr));
           passTopolCuts = (((fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr));
        }
        else { //LctopKpi
           passAllCuts = fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kAll,aod);
           passTopolCuts= fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kCandidate,aod);
        }
     }
     else { //everything that isn't Lc
         passAllCuts = fRDCutsAnalysis->IsSelected(d,AliRDHFCuts::kAll,aod);
         passTopolCuts = fRDCutsAnalysis->GetIsSelectedCuts();
     }
     
     if (fPdgSpecies != 431 && passTopolCuts == 0) continue;
     nSelectedNoPID++;
     fHistNEvents->Fill(9);
     if (fPdgSpecies == 431 && passAllCuts == 0) continue;
     if (passAllCuts) {
        nSelectedPID++;
        fHistNEvents->Fill(10);
     }
     //daughter subtraction from multiplicity, l.1169-1192
     
 //    Bool_t isPrimary = kTRUE;
 //    Double_t trueImpParXY = 9999.;
 //    Double_t impparXY = d->ImpParXY()*10000.;
 //     Double_t dlen = 0.1;
     

     
     
     Double_t mass[2];  //Calculate inv mass and check if in peak region
     if (fPdgSpecies == 411) {   //D+
         mass[0] = d->InvMass(nDau,pdgDau);
         mass[1] = -1.;
         if (TMath::Abs(mass[0]-mDplusPDG) < 0.02) nSelectedInMassPeak++;
     } else if (fPdgSpecies == 421) { //D0

         mass[0] = d->InvMass(2,pdgdaughtersD0);
         mass[1] = d->InvMass(2,pdgdaughtersD0bar);
         if (TMath::Abs(mass[0]-mD0PDG)<0.02 || TMath::Abs(mass[1]-mD0PDG) < 0.02) nSelectedInMassPeak++;
     } else if (fPdgSpecies == 413) { //D*
         mass[0] = dCascade->DeltaInvMass();
         mass[1] = -1.;
         if (TMath::Abs(mass[0] - (mDstarPDG-mD0PDG)) < 0.0015 ) nSelectedInMassPeak++;
     } else if (fPdgSpecies == 431) { //Ds
         mass[0] = d->InvMass(nDau,pdgDsKKpi);
         mass[1] = d->InvMass(nDau,pdgDspiKK);
         if (TMath::Abs(mass[0]-mDsPDG) < 0.02 || TMath::Abs(mass[1]-mDsPDG) < 0.02) nSelectedInMassPeak++;
     } else if (fPdgSpecies == 4122) { //Lc
         if (fLctoV0) { //pK0 inv mass
            mass[0] = d->InvMass(2,pdgDgLctopK0S);
            mass[1] = -1.;
            if (TMath::Abs(mass[0]-mLcPDG)< 0.02) nSelectedInMassPeak++;
         } else { //pKpi inv mass
            UInt_t pdgpKpi[3] = {2212,321,211};
            UInt_t pdgpiKp[3] = {211, 321,2212};
            if (passTopolCuts == 3 || passTopolCuts == 1) mass[0] = d->InvMass(3,pdgpKpi);
            if (passTopolCuts >= 2 ) mass[1] = d->InvMass(3,pdgpiKp);
            if (TMath::Abs(mass[0] - mLcPDG) < 0.02 ) nSelectedInMassPeak++;
         }
     }
            
     for (Int_t iHyp = 0; iHyp < 2; iHyp++) {
        if (mass[iHyp] < 0.) continue; //for D+ , D*, Lc2pK0, only one mass hypothesis
        Double_t invMass = mass[iHyp];
     
        if (fReadMC) {//TODO: add MC selections (l.1238-1279 in mult task) 
        }
      
        if (fPdgSpecies ==421) {
           if (iHyp == 0 && !(passAllCuts&1)) continue; // candidate not passing as D0
           if (iHyp == 1 && !(passAllCuts&2)) continue; // candidate not passing as D0bar
        }
        
        if (fPdgSpecies == 431) {
           if (iHyp == 0 && !(passAllCuts&4)) continue; // candidate Ds not passing as kk(phi)pi
           if (iHyp == 1 && !(passAllCuts&8)) continue; // candidate Ds not passing as pikk(phi)
        }
        if (fPdgSpecies == 4122 && !fLctoV0) {
           if (iHyp == 0 && !(passTopolCuts&1)) continue; // candidate not passing as Lc->pKpi
           if (iHyp == 1 && !(passTopolCuts&2)) continue; // candidate not passing as Lc->piKp
        }
     
      if (passAllCuts) {
         // if using nsparse, fill it
         
         //select phi region and fill appropriate histo
         Double_t phiCand = d->Phi();
         Double_t candDeltaPhi = phiCand - fPhiLeading; //delta-phi wrt leading particle
         Double_t arrayForSparse[4] = {ptCand, invMass, rtval, candDeltaPhi};  
         if (fUseNsparse) {fOutNsparse->Fill(arrayForSparse,weight); } 
        
           
         if (candDeltaPhi <= -TMath::PiOver2()) candDeltaPhi += TMath::TwoPi();
         if (candDeltaPhi > 3*TMath::PiOver2()) candDeltaPhi-=TMath::TwoPi();
         Double_t fDeltaPhiMinCut = TMath::DegToRad()*60.;
         Double_t fDeltaPhiMaxCut = TMath::DegToRad()*120.;
         Int_t region = 0;
         if ((candDeltaPhi < -fDeltaPhiMinCut) || (candDeltaPhi > 2*fDeltaPhiMaxCut)) region = -1; //left
         if ((candDeltaPhi >  fDeltaPhiMinCut) && (candDeltaPhi < fDeltaPhiMaxCut))   region = 1; //right
         
         if ((candDeltaPhi > -fDeltaPhiMinCut) && (candDeltaPhi < fDeltaPhiMinCut))   region = 2; // toward
         if ((candDeltaPhi > fDeltaPhiMinCut) && (candDeltaPhi < 2*fDeltaPhiMaxCut))   region = -2; //away
         
         if ((region == 1) || (region == -1)) { // fill transverse histo
         fPtvsMassvsRTTrans->Fill(rtval,invMass,ptCand,weight);
        
         }
         if (region == 2) { // fill toward histo
         fPtvsMassvsRTToward->Fill(rtval,invMass,ptCand,weight);
         }
         if (region == -2) { // fill away histo
         fPtvsMassvsRTAway->Fill(rtval,invMass,ptCand,weight);
         }
         
      }     
      }
    } 
  
  
  
  PostData(1, fOutput);
  PostData(2, fListCuts);
  PostData(3, fOutputCounters);
  PostData(4, fListQAhists);
}

Double_t AliAnalysisTaskSEDvsRT::CalculateRTVal(AliAODEvent* esdEvent)
{
   /// Calculate RT value for input event (method ported from AliAnalysisTaskUeSpectraRT)
   /// also sets phi of leading particle (fPhiLeading)
   Int_t runNumber = esdEvent->GetRunNumber();
   Int_t eventId = 0;
   Double_t trackRTval = -1;
   if (esdEvent->GetHeader()) eventId = GetEventIdAsLong(esdEvent->GetHeader());
   
   const Int_t nESDTracks = esdEvent->GetNumberOfTracks();
   TObjArray *fCTSTracks = new TObjArray();
   Double_t nRecTracks = 0;
   AliVParticle* part = 0x0;
   Double_t eta, pt, LeadingPt = -1; 
   for (Int_t iT = 0; iT < nESDTracks; iT++) 
   {
      part = esdEvent->GetTrack(iT);
      eta = part->Eta();
      pt  = part->Pt();
      if (TMath::Abs(eta) > fEtaCut) continue;
      if (!(TMath::Abs(pt) > 0.15)) continue;
      
      // Default track filter (to be checked)
      for ( Int_t i = 0; i < 1; i++)
      {
         UInt_t selectDebug = 0;
         if (fTrackFilter[i])
         {
            selectDebug = fTrackFilter[i]->IsSelected(part);
            if (!selectDebug)
            {
               continue;
            }
            /// fill tracks array
            fCTSTracks->Add(part);
            if (!part) continue;
         }
      }           
   }      
   
   //find leading object
   TObjArray *LeadingTrackReco = FindLeading(fCTSTracks);
   if (LeadingTrackReco) {
      AliVParticle* LeadingReco = 0;
      LeadingReco = (AliVParticle*)LeadingTrackReco->At(0);
      LeadingPt = LeadingReco->Pt();
      fPhiLeading = LeadingReco->Phi();
      if (LeadingPt > fLeadMin && LeadingPt < 300. ) {// calculate only if leading pt is in acceptable range
         //Sorting
         TObjArray *regionSortedParticlesReco = SortRegions((AliVParticle*)LeadingTrackReco->At(0), fCTSTracks);
         // Transverse regions
         TObjArray *regionsMinMaxReco = GetMinMaxRegion((TList*)regionSortedParticlesReco->At(2),(TList*)regionSortedParticlesReco->At(3));
         TList *listMax = (TList*)regionsMinMaxReco->At(0);
         TList *listMin = (TList*)regionsMinMaxReco->At(1);
         
         trackRTval = (listMax->GetEntries() + listMin->GetEntries()) / fAveMultiInTrans; //sum of transverse regions / average
         fHistPtLead->Fill(LeadingPt);
      }
      
   }
   
 
  return trackRTval; 
}


ULong64_t AliAnalysisTaskSEDvsRT::GetEventIdAsLong(AliVHeader* header)
{
//	unique ID for each event
   return ((ULong64_t)header->GetBunchCrossNumber() +
           (ULong64_t)header->GetOrbitNumber()*3564 +
           (ULong64_t)header->GetPeriodNumber()*16777215*3564);
}

void AliAnalysisTaskSEDvsRT::Terminate(Option_t */*option*/)
{
   /// Terminate analysis
   ///
   
   if (fDebug > 1) printf("AliAnalysisTaskSEDvsRT: Terminate() \n");
   fOutput = dynamic_cast<TList*> (GetOutputData(1));
   if (!fOutput) {
      printf("ERROR: fOutput not found\n");
      return;
   }
   
   fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));
   if (!fHistNEvents){
      printf("ERROR: fHistNEvents not found\n");
      return;
   }
   printf("Number of events analysed = %d\n",(Int_t)fHistNEvents->GetBinContent(3));
   return; 
}

TObjArray *AliAnalysisTaskSEDvsRT::FindLeading(TObjArray *array) 
{
   if (!array) return 0;
   Int_t nTracks = array->GetEntriesFast();
   if (!nTracks) return 0;
   
   TObjArray *tracks = new TObjArray(nTracks);
   for (Int_t ipart = 0; ipart < nTracks; ipart++) {
      AliVParticle *part = (AliVParticle*)(array->At(ipart));
      if(!part) continue;
      tracks->AddLast(part);
   }
   QSortTracks(*tracks, 0, tracks->GetEntriesFast());
   nTracks = tracks->GetEntriesFast();
   if (!nTracks) return 0;
   return tracks;
}

void AliAnalysisTaskSEDvsRT::QSortTracks(TObjArray &a, Int_t first, Int_t last)
{
   //Sort array by pT
  static TObject *tmp;
  static int i;           // "static" to save stack space
  int j;

  while (last - first > 1) {
    i = first;
    j = last;
    for (;;) {
      while (++i < last && ((AliVParticle*)a[i])->Pt() > ((AliVParticle*)a[first])->Pt() )
        ;
      while (--j > first && ((AliVParticle*)a[j])->Pt() < ((AliVParticle*)a[first])->Pt() )
        ;
      if (i >= j)
        break;

      tmp  = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
    if (j == first) {
      ++first;
      continue;
    }
    tmp = a[first];
    a[first] = a[j];
    a[j] = tmp;
    if (j - first < last - (j + 1)) {
      QSortTracks(a, first, j);
      first = j + 1;   
    } else {
      QSortTracks(a, j + 1, last);
      last = j;      
    }
  }
}


      
      
   



TObjArray *AliAnalysisTaskSEDvsRT::SortRegions(const AliVParticle* leading, TObjArray *array)
{
   if (!array) return 0;
   static const Double_t k60rad = 60.*TMath::Pi()/180.;
   static const Double_t k120rad = 120.*TMath::Pi()/180.;
   
   // define output lists of particles
   TList *toward = new TList();
   TList *away = new TList();
   TList *transverse1 = new TList();
   TList *transverse2 = new TList();
   TObjArray *regionParticles = new TObjArray;
   regionParticles->SetOwner();
   
   regionParticles->AddLast(toward);
   regionParticles->AddLast(away);
   regionParticles->AddLast(transverse1);
   regionParticles->AddLast(transverse2);
   if (!leading) return regionParticles;
   
   TVector3 leadVect(leading->Px(),leading->Py(),leading->Pz());
   Int_t nTracks = array->GetEntriesFast();
   if (!nTracks) return 0;

   //loop over tracks
   for (Int_t ipart = 0; ipart < nTracks; ipart++) {
      AliVParticle* part = (AliVParticle*)(array->At(ipart));
      if(!part) continue;
      //vector notation for particles
      TVector3 partVect(part->Px(), part->Py(), part->Pz());
      Int_t region = 0;
      Float_t deltaPhi = leadVect.DeltaPhi(partVect);
      if (deltaPhi <= -TMath::PiOver2()) deltaPhi+= TMath::TwoPi();
      if (deltaPhi > 3*TMath::PiOver2()) deltaPhi-= TMath::TwoPi();
      Double_t fUeDeltaPhiMinCut = TMath::DegToRad()*60.;
      Double_t fUeDeltaPhiMaxCut = TMath::DegToRad()*120.;

      //transverse regions

      if((deltaPhi<-fUeDeltaPhiMinCut) || (deltaPhi >2*fUeDeltaPhiMaxCut))region = -1; //left
      if((deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi <fUeDeltaPhiMaxCut)) region = 1;   //right
   
      if(deltaPhi > -fUeDeltaPhiMinCut && deltaPhi < fUeDeltaPhiMinCut) region = 2;    //forward
      if(deltaPhi > fUeDeltaPhiMaxCut && deltaPhi < 2*fUeDeltaPhiMaxCut) region = -2;  //backward
   
      // skip leading particle   
      if(leading == part) continue;
      if(part->Pt() >= leading->Pt()) continue;
      if(!region)continue;
   
      if(region == 1) transverse1->Add(part);
      if(region == -1) transverse2->Add(part);
      if(region == 2) toward->Add(part);
      if(region == -2) away->Add(part);
   }//end loop on tracks

   return regionParticles;   
}


TObjArray* AliAnalysisTaskSEDvsRT::GetMinMaxRegion(TList *transv1, TList *transv2)
{
// Returns two lists of particles, one for MIN and one for MAX region
  Double_t sumpT1 = 0.;
  Double_t sumpT2 = 0.;

  Int_t particles1 = transv1->GetEntries();
  Int_t particles2 = transv2->GetEntries();

// Loop on transverse region 1
  for(Int_t i=0; i<particles1; i++){
   AliVParticle *part = (AliVParticle*)transv1->At(i);
   sumpT1 +=  part->Pt();
   }

// Loop on transverse region 2
  for(Int_t i=0; i<particles2; i++){
   AliVParticle *part = (AliVParticle*)transv2->At(i);
   sumpT2 +=  part->Pt();
   }

  TObjArray *regionParticles = new TObjArray;
  if(sumpT2 >= sumpT1){
   regionParticles->AddLast(transv1); // MIN
   regionParticles->AddLast(transv2); // MAX 
   }
  else{
   regionParticles->AddLast(transv2); // MIN
   regionParticles->AddLast(transv1); // MAX
   }

  return regionParticles;
}

