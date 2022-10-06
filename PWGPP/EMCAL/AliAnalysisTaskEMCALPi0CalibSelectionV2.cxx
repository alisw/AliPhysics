/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *

 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <iostream>
#include <fstream>


// Root 
#include <TRefArray.h>
#include <TList.h>
#include <TH1F.h>
#include <TGeoManager.h>
#include <TFile.h>
#include "TChain.h"

// AliRoot
#include "AliAnalysisTaskEMCALPi0CalibSelectionV2.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliEMCALGeometry.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliEMCALRecoUtils.h"
#include "AliOADBContainer.h"
#include "AliDataFile.h"
#include "AliEmcalCorrectionTask.h"
#include "AliEmcalCorrectionComponent.h"
#include "AliEmcalCorrectionClusterNonLinearity.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEMCALPi0CalibSelectionV2) ;
/// \endcond

///
/// Default constructor. Arrays initialization is done here.
//______________________________________________________________________________________________
AliAnalysisTaskEMCALPi0CalibSelectionV2::AliAnalysisTaskEMCALPi0CalibSelectionV2() :
AliAnalysisTaskSE(),
fInputEvent(NULL), fMCEvent(NULL),  
fEMCALGeo(0x0),           fLoadMatrices(0),
fEMCALGeoName("EMCAL_COMPLETE12SMV1_DCAL_8SM"),
fTriggerName("EMC"),      
fRecoUtils(NULL),
fPeriodName(""),
fCorrTaskSetting(""),
fEMCALInitialized(kFALSE),
fIsMC(0),                 fSaveHistos(kTRUE),       fSaveCells(kFALSE),     
fSaveClusters(kFALSE),    fSaveFullTree(kFALSE),    fIsHeavyIon(kFALSE),
fNContributorsCutEnabled(kFALSE),
fOADBFilePath(""),        
fRecalPosition(kTRUE),
fEMCALCells(0x0),   
fOutputContainer(0x0),
fCheckCentrality(kFALSE), fCentralityClass("V0M"),  fCentWithEventSel(kFALSE),
fCentMin(-1),             fCentMax(10000000), 
fVertex(),                
fImportGeometryFromFile(1), fImportGeometryFilePath(""),
fCellEmin(0.05),
fEmin(0.5),               fEmax(15.),        
fL0min(0.01),             fL0max(0.5),
fOpAnglemin(0.),          fOpAnglemax(3.0),
fDTimeCut(100.),          fTimeMax(1000000),        fTimeMin(-1000000),
fAsyCut(1.),              fMinNCells(2),            fGroupNCells(0),
fSameSM(kFALSE),          
fNMaskCellColumns(11),    fMaskCellColumns(0x0),
fInvMassCutMin(110.),     fInvMassCutMax(160.),
fOfflineTriggerMask(0), isEMC(kFALSE), isDMC(kFALSE),
// Histograms binning
fNbins(300),              fMinBin(0.),              fMaxBin(300.),
fNTimeBins(1000),         fMinTimeBin(0.),          fMaxTimeBin(1000.),
fNEnergybins(1000),       fMinEnergyBin(0.),        fMaxEnergyBin(100.),
// Temporal
fMomentum1(),             fMomentum2(),             fMomentum12(),
// Histograms
fHmgg(0x0),               fHmggDifferentSM(0x0), 
fHmggMaskFrame(0x0),      fHmggDifferentSMMaskFrame(0x0), 
fHOpeningAngle(0x0),      fHOpeningAngleDifferentSM(0x0),  
fHAsymmetry(0x0),         fHAsymmetryDifferentSM(0x0),  
fhNEvents(0x0),
fhCentrality(0x0),        fhCentralitySelected(0x0),
fhClusterTime(0x0),
//Trees
fCellTree(NULL),
fVBuffer_NCells(0),         fBuffer_EventWeight(0),         fBuffer_ptHard(0),
fBuffer_Event_VertexZ(0),   fBuffer_EventNPrimaryTracks(0),
fBuffer_Event_V0Centrality(0),
fVBuffer_Cell_ID(0),        fVBuffer_Cell_E(0),                      
fVBuffer_Cell_t(0),         fVBuffer_Cell_gain(0),      fVBuffer_Cell_MCParticleID(0),
fBuffer_NClusters(0),
fVBuffer_Cluster_E(0),          fVBuffer_Cluster_Eta(0),    fVBuffer_Cluster_Phi(0),
fVBuffer_Cluster_LeadCellId(0), fVBuffer_Cluster_t(0),      fVBuffer_Cluster_M02(0),
fVBuffer_TrueCluster_MCId(0),   fVBuffer_Cluster_NCells(0),
fVFBuffer_Cluster_E(0),         fVFBuffer_Cluster_Eta(0),   fVFBuffer_Cluster_Phi(0),
fVBuffer_Cluster_X(0),          fVBuffer_Cluster_Y(0),      fVBuffer_Cluster_Z(0)
{
  

  for(Int_t iCell=0; iCell < kMaxActiveCells_calib; iCell++){
    fHmpi0[iCell] = 0;
  }
  
  fVertex[0]=fVertex[1]=fVertex[2]=-1000;
  
  fHTpi0[0]= 0 ;
  fHTpi0[1]= 0 ;
  fHTpi0[2]= 0 ;
  fHTpi0[3]= 0 ;
  
  fMaskCellColumns    = new Int_t[fNMaskCellColumns];
  fMaskCellColumns[0] = 6 ;  fMaskCellColumns[1] = 7 ;  fMaskCellColumns[2] = 8 ; 
  fMaskCellColumns[3] = 35;  fMaskCellColumns[4] = 36;  fMaskCellColumns[5] = 37; 
  fMaskCellColumns[6] = 12+AliEMCALGeoParams::fgkEMCALCols; fMaskCellColumns[7] = 13+AliEMCALGeoParams::fgkEMCALCols;
  fMaskCellColumns[8] = 40+AliEMCALGeoParams::fgkEMCALCols; fMaskCellColumns[9] = 41+AliEMCALGeoParams::fgkEMCALCols; 
  fMaskCellColumns[10]= 42+AliEMCALGeoParams::fgkEMCALCols; 

  
  for(Int_t iSM = 0; iSM < AliEMCALGeoParams::fgkEMCALModules; iSM++) {
    fHmggSM[iSM]                          = 0;
    fHmggSMMaskFrame[iSM]                 = 0;
    fHOpeningAngleSM[iSM]                 = 0;
    fHOpeningAnglePairSM[iSM]             = 0;
    fHAsymmetrySM[iSM]                    = 0;
    fHAsymmetryPairSM[iSM]                = 0;
    fhTowerDecayPhotonHit[iSM]            = 0;
    fhTowerDecayPhotonEnergy[iSM]         = 0;
    fhTowerDecayPhotonAsymmetry[iSM]      = 0;
    fhTowerDecayPhotonHitMaskFrame[iSM]   = 0;
    fMatrix[iSM]                          = 0x0;
    fhClusterTimeSM[iSM]                  = 0;
  }

}


/// Constructor with name as option. Arrays initialization is done here.
///
/// \param name: Name of task.
///
//______________________________________________________________________________________________
AliAnalysisTaskEMCALPi0CalibSelectionV2::AliAnalysisTaskEMCALPi0CalibSelectionV2(const char* name) :
AliAnalysisTaskSE(name),  
fInputEvent(NULL), fMCEvent(NULL),  
fEMCALGeo(0x0),           fLoadMatrices(0),
fEMCALGeoName("EMCAL_COMPLETE12SMV1_DCAL_8SM"),
fTriggerName("EMC"),      
fRecoUtils(NULL),
fPeriodName(""),
fCorrTaskSetting(""),
fEMCALInitialized(kFALSE),
fIsMC(0),                 fSaveHistos(kTRUE),        fSaveCells(kFALSE),     
fSaveClusters(kFALSE),    fSaveFullTree(kFALSE),     fIsHeavyIon(kFALSE),
fNContributorsCutEnabled(kFALSE),
fOADBFilePath(""),        
fRecalPosition(kTRUE),
fEMCALCells(0x0),
//fCuts(0x0),               
fOutputContainer(0x0),
fCheckCentrality(kFALSE), fCentralityClass("V0M"),  fCentWithEventSel(kFALSE),
fCentMin(-1),             fCentMax(10000000), 
fVertex(),
fImportGeometryFromFile(1), fImportGeometryFilePath(""),
fCellEmin(0.05),
fEmin(0.5),               fEmax(15.),     
fL0min(0.01),             fL0max(0.5),
fOpAnglemin(0.),          fOpAnglemax(3.0),
fDTimeCut(100.),          fTimeMax(1000000),        fTimeMin(-1000000),
fAsyCut(1.),              fMinNCells(2),            fGroupNCells(0),
fSameSM(kFALSE),          
fNMaskCellColumns(11),    fMaskCellColumns(0x0),
fInvMassCutMin(110.),     fInvMassCutMax(160.),
fOfflineTriggerMask(0), isEMC(kFALSE), isDMC(kFALSE),
// Histograms binning
fNbins(300),              fMinBin(0.),              fMaxBin(300.),
fNTimeBins(1000),         fMinTimeBin(0.),          fMaxTimeBin(1000.),
fNEnergybins(1000),       fMinEnergyBin(0.),        fMaxEnergyBin(100.),
// Temporal
fMomentum1(),             fMomentum2(),             fMomentum12(),
// Histograms
fHmgg(0x0),               fHmggDifferentSM(0x0), 
fHmggMaskFrame(0x0),      fHmggDifferentSMMaskFrame(0x0), 
fHOpeningAngle(0x0),      fHOpeningAngleDifferentSM(0x0),  
fHAsymmetry(0x0),         fHAsymmetryDifferentSM(0x0),  
fhNEvents(0x0),
fhCentrality(0x0),        fhCentralitySelected(0x0),
fhClusterTime(0x0),
//Trees
fCellTree(NULL),
fVBuffer_NCells(0),         fBuffer_EventWeight(0), fBuffer_ptHard(0),
fBuffer_Event_VertexZ(0),   fBuffer_EventNPrimaryTracks(0),
fBuffer_Event_V0Centrality(0),
fVBuffer_Cell_ID(0),        fVBuffer_Cell_E(0),                      
fVBuffer_Cell_t(0),         fVBuffer_Cell_gain(0),      fVBuffer_Cell_MCParticleID(0),
fBuffer_NClusters(0),
fVBuffer_Cluster_E(0),          fVBuffer_Cluster_Eta(0),    fVBuffer_Cluster_Phi(0),
fVBuffer_Cluster_LeadCellId(0), fVBuffer_Cluster_t(0),      fVBuffer_Cluster_M02(0),
fVBuffer_TrueCluster_MCId(0),   fVBuffer_Cluster_NCells(0),
fVFBuffer_Cluster_E(0),         fVFBuffer_Cluster_Eta(0),   fVFBuffer_Cluster_Phi(0),
fVBuffer_Cluster_X(0),          fVBuffer_Cluster_Y(0),      fVBuffer_Cluster_Z(0)
{
  
  for(Int_t iCell=0; iCell < kMaxActiveCells_calib; iCell++){
    fHmpi0[iCell] = 0;
  }
  
  fVertex[0]=fVertex[1]=fVertex[2]=-1000;
  
  fHTpi0[0]= 0 ;
  fHTpi0[1]= 0 ;
  fHTpi0[2]= 0 ;
  fHTpi0[3]= 0 ;
  
  fMaskCellColumns    = new Int_t[fNMaskCellColumns];
  fMaskCellColumns[0] = 6 ;  fMaskCellColumns[1] = 7 ;  fMaskCellColumns[2] = 8 ; 
  fMaskCellColumns[3] = 35;  fMaskCellColumns[4] = 36;  fMaskCellColumns[5] = 37; 
  fMaskCellColumns[6] = 12+AliEMCALGeoParams::fgkEMCALCols; fMaskCellColumns[7] = 13+AliEMCALGeoParams::fgkEMCALCols;
  fMaskCellColumns[8] = 40+AliEMCALGeoParams::fgkEMCALCols; fMaskCellColumns[9] = 41+AliEMCALGeoParams::fgkEMCALCols; 
  fMaskCellColumns[10]= 42+AliEMCALGeoParams::fgkEMCALCols; 

  
  for(Int_t iSM = 0; iSM < AliEMCALGeoParams::fgkEMCALModules; iSM++) {
    fHmggSM[iSM]                     = 0;
    fHmggSMMaskFrame[iSM]            = 0;
    fHOpeningAngleSM[iSM]            = 0;
    fHOpeningAnglePairSM[iSM]        = 0;
    fHAsymmetrySM[iSM]               = 0;
    fHAsymmetryPairSM[iSM]           = 0;
    fhTowerDecayPhotonHit[iSM]       = 0;
    fhTowerDecayPhotonEnergy[iSM]    = 0;
    fhTowerDecayPhotonAsymmetry[iSM] = 0;
    fhTowerDecayPhotonHitMaskFrame[iSM]= 0;
    fMatrix[iSM]                     = 0x0;
    fhClusterTimeSM[iSM]             = 0;
  }
    
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
//DefineOutput(2, TList::Class());  // will contain cuts or local params
}

///
/// Destructor.
//_____________________________________________________________________________
AliAnalysisTaskEMCALPi0CalibSelectionV2::~AliAnalysisTaskEMCALPi0CalibSelectionV2() {
  
  if(fOutputContainer) {
    fOutputContainer->Delete() ; 
    delete fOutputContainer ;
  }
  
  if(fEMCALGeo)         delete fEMCALGeo  ;
  if(fNMaskCellColumns) delete [] fMaskCellColumns;
}

///
/// Init geometry and set the geometry matrix, for the first event, skip the rest.
/// Also set once the run dependent calibrations.
//________________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::InitGeometryMatrices() {
  Int_t runnumber = InputEvent()->GetRunNumber() ;
  
  //
  // Load default geo matrices if requested
  if(fImportGeometryFromFile && !gGeoManager) {
    if(fImportGeometryFilePath=="") { // If not specified, set location depending on run number
      // "$ALICE_ROOT/EVE/alice-data/default_geo.root"
      if     (runnumber <  140000) fImportGeometryFilePath = AliDataFile::GetFileNameOADB("EMCAL/geometry_2010.root").data();
      else if(runnumber <  171000) fImportGeometryFilePath = AliDataFile::GetFileNameOADB("EMCAL/geometry_2011.root").data();
      else if(runnumber <  198000) fImportGeometryFilePath = AliDataFile::GetFileNameOADB("EMCAL/geometry_2012.root").data(); // 2012-2013
      else                         fImportGeometryFilePath = AliDataFile::GetFileNameOADB("EMCAL/geometry_2015.root").data(); // >= 2015
    }
    AliInfo(Form("Import %s",fImportGeometryFilePath.Data()));
    TGeoManager::Import(fImportGeometryFilePath) ; // default need file "geometry.root" in local dir!!!!
  }
  
  //
  if(fLoadMatrices) { // Load matrices
    AliInfo("Load user defined EMCAL geometry matrices");
    // OADB if available
    AliOADBContainer emcGeoMat("AliEMCALgeo");

    if(fOADBFilePath!="")
      emcGeoMat.InitFromFile(Form("%s/EMCALlocal2master.root",fOADBFilePath.Data()),"AliEMCALgeo");
    else
      emcGeoMat.InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALlocal2master.root").data(),"AliEMCALgeo");

    TObjArray *matEMCAL=(TObjArray*)emcGeoMat.GetObject(runnumber,"EmcalMatrices");
    
    for(Int_t mod = 0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++) { //SM loop
      
      if (!fMatrix[mod]) { // Get it from OADB
        AliDebug(1,Form("EMCAL matrices SM %d, %p",mod,((TGeoHMatrix*) matEMCAL->At(mod))));
        fMatrix[mod] = (TGeoHMatrix*) matEMCAL->At(mod) ;
      }        
      
      if(fMatrix[mod]) {  
        if(DebugLevel() > 1) 
          fMatrix[mod]->Print();
        fEMCALGeo->SetMisalMatrix(fMatrix[mod],mod) ;  
      } else if(gGeoManager) {
        AliWarning(Form("Set matrix for SM %d from gGeoManager",mod));
        fEMCALGeo->SetMisalMatrix(fEMCALGeo->GetMatrixForSuperModuleFromGeoManager(mod),mod) ;
      } else {
        AliError(Form("Alignment matrix for SM %d is not available",mod));
      }
    }
  } else if(!gGeoManager) { // Load matrices from Data
    AliInfo("Get geo matrices from data");
    //Still not implemented in AOD, just a workaround to be able to work at least with ESDs	
    if(!strcmp(InputEvent()->GetName(),"AliAODEvent")) {
      AliWarning("Use ideal geometry, values geometry matrix not kept in AODs");
    } else {// ESD
      AliDebug(1,"AliAnalysisTaskEMCALClusterize Load Misaligned matrices");
      
      for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++) {
        if(InputEvent()->GetEMCALMatrix(mod)) {
          if(DebugLevel() > 1) 
            InputEvent()->GetEMCALMatrix(mod)->Print();
          fEMCALGeo->SetMisalMatrix(InputEvent()->GetEMCALMatrix(mod),mod) ;
        } 
      }
    }
  } else if(gGeoManager) { // Load default gGeoManager matrices
    for(Int_t mod = 0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++) {
      AliWarning(Form("Set matrix for SM %d from gGeoManager",mod));
      fEMCALGeo->SetMisalMatrix(fEMCALGeo->GetMatrixForSuperModuleFromGeoManager(mod),mod) ;
    }
  }

}

///
/// Initialize EMCAL
//__________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::InitializeEMCAL(){
    AliEmcalCorrectionTask* emcalCorrTask=0x0;
    //     

    if( !fCorrTaskSetting.CompareTo("") ){
      emcalCorrTask  = (AliEmcalCorrectionTask*) AliAnalysisManager::GetAnalysisManager()->GetTopTasks()->FindObject("AliEmcalCorrectionTask");
    } else {
      emcalCorrTask  = (AliEmcalCorrectionTask*) AliAnalysisManager::GetAnalysisManager()->GetTopTasks()->FindObject( Form("AliEmcalCorrectionTask_%s",fCorrTaskSetting.Data()) );
    }

    if( emcalCorrTask ){
      AliEmcalCorrectionComponent * emcalCorrComponent = 0x0;
      if( !fCorrTaskSetting.CompareTo("") ){
        emcalCorrComponent = emcalCorrTask->GetCorrectionComponent("AliEmcalCorrectionClusterNonLinearity");
      } else {
        emcalCorrComponent = emcalCorrTask->GetCorrectionComponent( Form("AliEmcalCorrectionClusterNonLinearity_%s", fCorrTaskSetting.Data()) );
      }
      if( emcalCorrComponent ){
        fRecoUtils        = emcalCorrComponent->GetRecoUtils();
      } else {
        if( !fCorrTaskSetting.CompareTo("") ){
          emcalCorrComponent = emcalCorrTask->GetCorrectionComponent("AliEmcalCorrectionCellBadChannel");
        } else {
          emcalCorrComponent = emcalCorrTask->GetCorrectionComponent( Form("AliEmcalCorrectionCellBadChannel_%s", fCorrTaskSetting.Data()) );
        }
        if( emcalCorrComponent ){
          fRecoUtils        = emcalCorrComponent->GetRecoUtils();
        }
      }
    }
    
    if (fRecoUtils) {
      fEMCALInitialized = kTRUE;
      fRecoUtils->SetNumberOfCellsFromEMCALBorder(0);
    }


    return;
}

///
/// Create output container, init geometry.
//___________________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::UserCreateOutputObjects() {

  if(!fEMCALGeo)fEMCALGeo =  AliEMCALGeometry::GetInstance(fEMCALGeoName) ;
  Int_t nSM = (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules();
  
  fOutputContainer = new TList();
  const Int_t buffersize = 255;
  char hname[buffersize], htitl[buffersize];
  
  fhNEvents        = new TH1I("hNEvents", "Number of analyzed events"   , 1 , 0 , 1  ) ;
  fOutputContainer->Add(fhNEvents);
  
  if ( fCheckCentrality ) {
    fhCentrality   = new TH1F ("hCentrality",
                                Form("Number of events in centrality bins, |vz|<10 cm, method <%s> ",fCentralityClass.Data()),
                                1000,-1000.,1000.) ;
    fhCentrality->SetXTitle("Centrality bin");
    fOutputContainer->Add(fhCentrality) ;  
    
    fhCentralitySelected   = new TH1F ("hCentralitySelected",
                                        Form("Number of selected events in centrality bin, |vz|<10 cm, method <%s> ",fCentralityClass.Data()),
                                        100,0.,100.) ;
    fhCentralitySelected->SetXTitle("Centrality bin");
    fOutputContainer->Add(fhCentralitySelected) ;  
  }

  if( fSaveHistos ){
    fHmgg = new TH2F("hmgg","2-cluster invariant mass",fNbins,fMinBin,fMaxBin,100,0,10);
    fHmgg->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
    fHmgg->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
    fOutputContainer->Add(fHmgg);

    fHmggDifferentSM = new TH2F("hmggDifferentSM","2-cluster invariant mass, different SM",fNbins,fMinBin,fMaxBin,100,0,10);
    fHmggDifferentSM->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
    fHmggDifferentSM->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
    fOutputContainer->Add(fHmggDifferentSM);

    fHOpeningAngle = new TH2F("hopang","2-cluster opening angle",100,0.,50.,100,0,10);
    fHOpeningAngle->SetXTitle("#alpha_{#gamma #gamma}");
    fHOpeningAngle->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
    fOutputContainer->Add(fHOpeningAngle);

    fHOpeningAngleDifferentSM = new TH2F("hopangDifferentSM","2-cluster opening angle, different SM",100,0,50.,100,0,10);
    fHOpeningAngleDifferentSM->SetXTitle("#alpha_{#gamma #gamma}");
    fHOpeningAngleDifferentSM->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
    fOutputContainer->Add(fHOpeningAngleDifferentSM);

    fHAsymmetry = new TH2F("hasym","2-cluster opening angle",100,0.,1.,100,0,10);
    fHAsymmetry->SetXTitle("a");
    fHAsymmetry->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
    fOutputContainer->Add(fHAsymmetry);

    fHAsymmetryDifferentSM = new TH2F("hasymDifferentSM","2-cluster opening angle, different SM",100,0,1.,100,0,10);
    fHAsymmetryDifferentSM->SetXTitle("a");
    fHAsymmetryDifferentSM->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
    fOutputContainer->Add(fHAsymmetryDifferentSM);

    fHmggMaskFrame = new TH2F("hmggMaskFrame","2-cluster invariant mass, frame masked",fNbins,fMinBin,fMaxBin,100,0,10);
    fHmggMaskFrame->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
    fHmggMaskFrame->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
    fOutputContainer->Add(fHmggMaskFrame);

    fHmggDifferentSMMaskFrame = new TH2F("hmggDifferentSMMaskFrame","2-cluster invariant mass, different SM, frame masked",
                                         fNbins,fMinBin,fMaxBin,100,0,10);
    fHmggDifferentSMMaskFrame->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
    fHmggDifferentSMMaskFrame->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
    fOutputContainer->Add(fHmggDifferentSMMaskFrame);

    for(Int_t iSM = 0; iSM < nSM; iSM++) {
      snprintf(hname, buffersize, "hmgg_SM%d",iSM);
      snprintf(htitl, buffersize, "Two-gamma inv. mass for super mod %d",iSM);
      fHmggSM[iSM] = new TH2F(hname,htitl,fNbins,fMinBin,fMaxBin,100,0,10);
      fHmggSM[iSM]->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
      fHmggSM[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
      fOutputContainer->Add(fHmggSM[iSM]);

      snprintf(hname, buffersize, "hmgg_SM%d_MaskFrame",iSM);
      snprintf(htitl, buffersize, "Two-gamma inv. mass for super mod %d",iSM);
      fHmggSMMaskFrame[iSM] = new TH2F(hname,htitl,fNbins,fMinBin,fMaxBin,100,0,10);
      fHmggSMMaskFrame[iSM]->SetXTitle("m_{#gamma #gamma} (MeV/c^{2})");
      fHmggSMMaskFrame[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
      fOutputContainer->Add(fHmggSMMaskFrame[iSM]);

      snprintf(hname, buffersize, "hopang_SM%d",iSM);
      snprintf(htitl, buffersize, "Opening angle for super mod %d",iSM);
      fHOpeningAngleSM[iSM] = new TH2F(hname,htitl,100,0.,50.,100,0,10);
      fHOpeningAngleSM[iSM]->SetXTitle("#alpha_{#gamma #gamma} (deg)");
      fHOpeningAngleSM[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
      fOutputContainer->Add(fHOpeningAngleSM[iSM]);

      snprintf(hname,buffersize, "hopang_PairSM%d",iSM);
      snprintf(htitl,buffersize, "Opening angle for SM pair: %d",iSM);
      fHOpeningAnglePairSM[iSM] = new TH2F(hname,htitl,100,0.,50.,100,0,10);
      fHOpeningAnglePairSM[iSM]->SetXTitle("#alpha_{#gamma #gamma} (deg)");
      fHOpeningAnglePairSM[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
      fOutputContainer->Add(fHOpeningAnglePairSM[iSM]);    

      snprintf(hname, buffersize, "hasym_SM%d",iSM);
      snprintf(htitl, buffersize, "Asymmetry for super mod %d",iSM);
      fHAsymmetrySM[iSM] = new TH2F(hname,htitl,100,0.,1.,100,0,10);
      fHAsymmetrySM[iSM]->SetXTitle("a");
      fHAsymmetrySM[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
      fOutputContainer->Add(fHAsymmetrySM[iSM]);

      snprintf(hname,buffersize, "hasym_PairSM%d",iSM);
      snprintf(htitl,buffersize, "Asymmetry for SM pair: %d",iSM);
      fHAsymmetryPairSM[iSM] = new TH2F(hname,htitl,100,0.,1.,100,0,10);
      fHAsymmetryPairSM[iSM]->SetXTitle("a");
      fHAsymmetryPairSM[iSM]->SetYTitle("p_{T #gamma #gamma} (GeV/c)");
      fOutputContainer->Add(fHAsymmetryPairSM[iSM]);    

      Int_t colmax = 48;
      Int_t rowmax = 24;

      fhTowerDecayPhotonHit[iSM]  = new TH2F (Form("hTowerDecPhotonHit_Mod%d",iSM),
                                              Form("Entries in grid of cells in Module %d",iSM), 
                                              colmax+2,-1.5,colmax+0.5, rowmax+2,-1.5,rowmax+0.5); 
      fhTowerDecayPhotonHit[iSM]->SetYTitle("row (phi direction)");
      fhTowerDecayPhotonHit[iSM]->SetXTitle("column (eta direction)");
      fOutputContainer->Add(fhTowerDecayPhotonHit[iSM]);

      fhTowerDecayPhotonEnergy[iSM]  = new TH2F (Form("hTowerDecPhotonEnergy_Mod%d",iSM),
                                                 Form("Accumulated energy in grid of cells in Module %d",iSM), 
                                                 colmax+2,-1.5,colmax+0.5, rowmax+2,-1.5,rowmax+0.5); 
      fhTowerDecayPhotonEnergy[iSM]->SetYTitle("row (phi direction)");
      fhTowerDecayPhotonEnergy[iSM]->SetXTitle("column (eta direction)");
      fOutputContainer->Add(fhTowerDecayPhotonEnergy[iSM]);

      fhTowerDecayPhotonAsymmetry[iSM]  = new TH2F (Form("hTowerDecPhotonAsymmetry_Mod%d",iSM),
                                                    Form("Accumulated asymmetry in grid of cells in Module %d",iSM), 
                                                    colmax+2,-1.5,colmax+0.5, rowmax+2,-1.5,rowmax+0.5); 
      fhTowerDecayPhotonAsymmetry[iSM]->SetYTitle("row (phi direction)");
      fhTowerDecayPhotonAsymmetry[iSM]->SetXTitle("column (eta direction)");
      fOutputContainer->Add(fhTowerDecayPhotonAsymmetry[iSM]);

      fhTowerDecayPhotonHitMaskFrame[iSM]  = new TH2F (Form("hTowerDecPhotonHit_Mod%d_MaskFrame",iSM),Form("Entries in grid of cells in Module %d",iSM), 
                                                       colmax+2,-1.5,colmax+0.5, rowmax+2,-1.5,rowmax+0.5); 
      fhTowerDecayPhotonHitMaskFrame[iSM]->SetYTitle("row (phi direction)");
      fhTowerDecayPhotonHitMaskFrame[iSM]->SetXTitle("column (eta direction)");
      fOutputContainer->Add(fhTowerDecayPhotonHitMaskFrame[iSM]);

      fhClusterTimeSM[iSM] = new TH2F(Form("hClusterTime_SM%d",iSM),"cluster time vs E",100,0,10, 200,-1000,1000);
      fhClusterTimeSM[iSM]->SetXTitle("E (GeV)");
      fhClusterTimeSM[iSM]->SetYTitle("t (ns)");
      fOutputContainer->Add(fhClusterTimeSM[iSM]);
    }  

    Int_t nchannels = nSM*AliEMCALGeoParams::fgkEMCALRows*AliEMCALGeoParams::fgkEMCALCols;
    for(Int_t ibc = 0; ibc < 4; ibc++) {
      fHTpi0[ibc] = new TH2F(Form("hTime_BC%d",ibc),Form("Time of cell clusters under pi0 peak, bunch crossing %d",ibc),
                             nchannels,0,nchannels, fNTimeBins,fMinTimeBin,fMaxTimeBin);
      fOutputContainer->Add(fHTpi0[ibc]);       
      fHTpi0[ibc]->SetYTitle("time (ns)");
      fHTpi0[ibc]->SetXTitle("abs. Id. ");
    }

    fhClusterTime = new TH2F("hClusterTime","cluster time vs E",100,0,10, 100,0,1000);
    fhClusterTime->SetXTitle("E (GeV)");
    fhClusterTime->SetYTitle("t (ns)");
    fOutputContainer->Add(fhClusterTime);

    // for(Int_t iCell=0; iCell < kMaxActiveCells_calib; iCell++){
    //   fHmpi0[iCell]    = new TH1F( Form("%d",iCell), Form("Two-gamma inv. mass for cell %d",iCell), fNbins,fMinBin,fMaxBin);
    //   fHmpi0[iCell]->SetXTitle("mass (MeV/c^{2}");
    //   fOutputContainer->Add(fHmpi0[iCell]);
    // }
  }

  fOutputContainer->SetOwner(kTRUE);
  PostData(1,fOutputContainer);


  if( fSaveCells || fSaveClusters ){
    fCellTree = new TTree(Form("EMCALCells_%s",fTriggerName.Data()),Form("EMCALCells_%s",fTriggerName.Data()));

    if( fIsMC > 1 ) {
      fCellTree->Branch("EventWeight",          &fBuffer_EventWeight, "EventWeight/D");
      fCellTree->Branch("ptHard",               &fBuffer_ptHard,      "ptHard/F");
    }

    fCellTree->Branch("VertexZ",                &fBuffer_Event_VertexZ,       "VertexZ/S");
    fCellTree->Branch("PrimaryTracks",          &fBuffer_EventNPrimaryTracks, "PrimaryTracks/s");

    if( fIsHeavyIon ){
      fCellTree->Branch("V0Centrality",         &fBuffer_Event_V0Centrality,  "V0Centrality/F");
    }
  }

  if( fSaveCells ){
      fCellTree->Branch("NCells",                 &fVBuffer_NCells,             "NCells/I");
      fCellTree->Branch("Cell_ID",                "std::vector<UShort_t>",      &fVBuffer_Cell_ID);
      if( fSaveFullTree )
        fCellTree->Branch("Cell_E",               "std::vector<Float_t>",       &fVBuffer_Cell_E);
      else
        fCellTree->Branch("Cell_E",               "std::vector<UShort_t>",      &fVBuffer_Cell_E);
      fCellTree->Branch("Cell_time",              "std::vector<Short_t>",       &fVBuffer_Cell_t);
      fCellTree->Branch("Cell_highGain",          "std::vector<Bool_t>",        &fVBuffer_Cell_gain);

      if( fIsMC ){
        fCellTree->Branch("Cell_MCParticleID",    "std::vector<Short_t>",       &fVBuffer_Cell_MCParticleID);
      }

  }

  if( fSaveClusters ){
    fCellTree->Branch("NClusters",            &fBuffer_NClusters,               "NClusters/s");
    fCellTree->Branch("Cluster_LeadCellId",   "std::vector<UShort_t>",             &fVBuffer_Cluster_LeadCellId);
    fCellTree->Branch("Cluster_NCells",       "std::vector<Int_t>",             &fVBuffer_Cluster_NCells);
    fCellTree->Branch("Cluster_t",            "std::vector<Short_t>",           &fVBuffer_Cluster_t);
    if( fSaveFullTree ){
      fCellTree->Branch("Cluster_E_long",            "std::vector<Float_t>",          &fVFBuffer_Cluster_E);
      fCellTree->Branch("Cluster_Eta_long",          "std::vector<Float_t>",          &fVFBuffer_Cluster_Eta);
      fCellTree->Branch("Cluster_Phi_long",          "std::vector<Float_t>",          &fVFBuffer_Cluster_Phi);
      fCellTree->Branch("Cluster_M02",          "std::vector<Double_t>",         &fVBuffer_Cluster_M02);
      fCellTree->Branch("Cluster_X",            "std::vector<Float_t>",          &fVBuffer_Cluster_X);
      fCellTree->Branch("Cluster_Y",            "std::vector<Float_t>",          &fVBuffer_Cluster_Y);
      fCellTree->Branch("Cluster_Z",            "std::vector<Float_t>",          &fVBuffer_Cluster_Z);
    } else {
      fCellTree->Branch("Cluster_E",            "std::vector<UShort_t>",         &fVBuffer_Cluster_E);
      fCellTree->Branch("Cluster_Eta",          "std::vector<Short_t>",          &fVBuffer_Cluster_Eta);
      fCellTree->Branch("Cluster_Phi",          "std::vector<UShort_t>",         &fVBuffer_Cluster_Phi);
    }
    
    if( fIsMC ){
      fCellTree->Branch("TrueCluster_MCId",    "std::vector<Short_t>",          &fVBuffer_TrueCluster_MCId);
    }
  }

  if( fSaveCells || fSaveClusters ){
    OpenFile(2);
    PostData(2, fCellTree);
  }
}

///
/// Fill the invariant mass analysis per channel with the
/// corrected clusters, and other general histograms.
//__________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::FillHistograms() {

  Int_t absId1   = -1;
  Int_t iSupMod1 = -1;
  Int_t iphi1    = -1;
  Int_t ieta1    = -1;
  Int_t absId2   = -1;
  Int_t iSupMod2 = -1;
  Int_t iphi2    = -1;
  Int_t ieta2    = -1;
  Bool_t shared  = kFALSE;
  
  Int_t bc  = InputEvent()->GetBunchCrossNumber();
  TClonesArray* clusters;

  if( !fCorrTaskSetting.CompareTo("") ){
    clusters = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject("caloClusters"));
  } else {
    clusters = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject( Form("%sClustersBranch", fCorrTaskSetting.Data() ) ));
  }

  if ( !clusters ) return;
  
  for(Int_t iClu=0; iClu<clusters->GetEntries(); iClu++) {
    
    // AliVCluster *c1 = (AliVCluster *) clusters->At(iClu);
    AliVCluster *c1   = NULL;
    if( fInputEvent->IsA()==AliESDEvent::Class() ){
      c1 = new AliESDCaloCluster( *(AliESDCaloCluster*) clusters->At(iClu) );
    } else if (fInputEvent->IsA()==AliAODEvent::Class() ){
      c1 = new AliAODCaloCluster( *(AliAODCaloCluster*) clusters->At(iClu) );
    }
    
    if (!AcceptCluster(c1)) continue;
    
    fRecoUtils->GetMaxEnergyCell(fEMCALGeo, fEMCALCells,c1,absId1,iSupMod1,ieta1,iphi1,shared); //RECOUTILS

    if( isEMC && !(isDMC) &&  iSupMod1 > 11 ) continue;
    if( isDMC && !(isEMC) &&  iSupMod1 < 12 ) continue;
    
    c1->GetMomentum(fMomentum1,fVertex);
    
    // Check if cluster is in fidutial region, not too close to borders
    if(fEMCALGeo == NULL || fEMCALCells == NULL) continue;
    Bool_t in1 = kTRUE;
    fRecoUtils->CheckCellFiducialRegion(fEMCALGeo, c1, fEMCALCells); //RECOUTILS
    
    // Clusters not facing frame structures
    Bool_t mask1 = MaskFrameCluster(iSupMod1, ieta1);
    
    Double_t time1 = c1->GetTOF()*1.e9;
    
    fhClusterTime            ->Fill(c1->E(),time1);
    fhClusterTimeSM[iSupMod1]->Fill(c1->E(),time1);
    
    // Combine cluster with other clusters and get the invariant mass
    for (Int_t jClu=iClu+1; jClu < clusters->GetEntries(); jClu++) {
      
      // AliAODCaloCluster *c2 = (AliAODCaloCluster *) clusters->At(jClu);
      AliVCluster *c2   = NULL;
      if( fInputEvent->IsA()==AliESDEvent::Class() ){
        c2 = new AliESDCaloCluster( *(AliESDCaloCluster*) clusters->At(jClu) );
      } else if (fInputEvent->IsA()==AliAODEvent::Class() ){
        c2 = new AliAODCaloCluster( *(AliAODCaloCluster*) clusters->At(jClu) );
      }
    
      if (!AcceptCluster(c2)) continue;
      
      fRecoUtils->GetMaxEnergyCell(fEMCALGeo, fEMCALCells,c2,absId2,iSupMod2,ieta2,iphi2,shared); //RECOUTILS

      if( isEMC && !(isDMC) &&  iSupMod2 > 11 ) continue;
      if( isDMC && !(isEMC) &&  iSupMod2 < 12 ) continue;
      
      c2->GetMomentum(fMomentum2,fVertex);
      
      fMomentum12 = fMomentum1+fMomentum2;
      Float_t invmass = fMomentum12.M()*1000; 
      
      //Asimetry cut      
      Float_t asym = TMath::Abs(fMomentum1.E()-fMomentum2.E())/(fMomentum1.E()+fMomentum2.E());
      if(asym > fAsyCut) continue;
      
      //Time cut
      Double_t time2 = c2->GetTOF()*1.e9;
      
      // fhClusterPairDiffTime->Fill(fMomentum12.E(),time1-time2);
      if(TMath::Abs(time1-time2) > fDTimeCut) continue;

      //Opening angle cut to reject too close clusters
      if( ((fMomentum1.Angle(fMomentum2.Vect()) < fOpAnglemin) || (fMomentum1.Angle(fMomentum2.Vect()) > fOpAnglemax)))  continue;
      
      if(invmass < fMaxBin && invmass > fMinBin ) {
        //Check if cluster is in fidutial region, not too close to borders
        Bool_t in2 = kTRUE;
        fRecoUtils->CheckCellFiducialRegion(fEMCALGeo, c2, fEMCALCells); //RECOUTILS
        
        // Clusters not facing frame structures
        Bool_t mask2 = MaskFrameCluster(iSupMod2, ieta2);         
        
        if(in1 && in2) { //acceptance cuts
          fHmgg->Fill(invmass,fMomentum12.Pt()); 
          
          if(iSupMod1==iSupMod2) {
            fHmggSM[iSupMod1]->Fill(invmass,fMomentum12.Pt()); 
            // fhClusterPairDiffTimeSameSM[iSupMod1]->Fill(fMomentum12.E(),time1-time2);
          } else {           
            fHmggDifferentSM ->Fill(invmass,fMomentum12.Pt());
          }        
          
          if(!mask1 && !mask2) {  // Pair not facing frame
            fHmggMaskFrame->Fill(invmass,fMomentum12.Pt()); 
            
            if(iSupMod1==iSupMod2) fHmggSMMaskFrame[iSupMod1]->Fill(invmass,fMomentum12.Pt()); 
            else                   fHmggDifferentSMMaskFrame ->Fill(invmass,fMomentum12.Pt());
          }
          
          if(invmass > fInvMassCutMin && invmass < fInvMassCutMax) { //restrict to clusters really close to pi0 peak (pair in 100 < mass < 160)

            // Check time of cells in both clusters, and fill time histogram
            for(Int_t icell = 0; icell < c1->GetNCells(); icell++) {
              Int_t absID = c1->GetCellAbsId(icell);   
              fHTpi0[bc%4]->Fill(absID, fEMCALCells->GetCellTime(absID)*1.e9);  
            }
            
            for(Int_t icell = 0; icell < c2->GetNCells(); icell++) {
              Int_t absID = c2->GetCellAbsId(icell);   
              fHTpi0[bc%4]->Fill(absID, fEMCALCells->GetCellTime(absID)*1.e9);  
            }
            
            //Opening angle of 2 photons
            Float_t opangle = fMomentum1.Angle(fMomentum2.Vect())*TMath::RadToDeg();
            
            fHOpeningAngle ->Fill(opangle,fMomentum12.Pt()); 
            fHAsymmetry    ->Fill(asym,fMomentum12.Pt()); 
            
            if(iSupMod1==iSupMod2) {
              fHOpeningAngleSM[iSupMod1] ->Fill(opangle,fMomentum12.Pt());
              fHAsymmetrySM[iSupMod1]    ->Fill(asym,fMomentum12.Pt());
            } else {      
              fHOpeningAngleDifferentSM  ->Fill(opangle,fMomentum12.Pt());
              fHAsymmetryDifferentSM     ->Fill(asym,fMomentum12.Pt());
            }
            
          }
          
        }
        
        //In case of filling only channels with second cluster in same SM
        // if(fSameSM && iSupMod1!=iSupMod2) continue;
        
        // if (fGroupNCells == 0) {
        //   fHmpi0[absId1]->Fill(invmass);
        //   fHmpi0[absId2]->Fill(invmass);
          
          // if(invmass > fInvMassCutMin && invmass < fInvMassCutMax) { //restrict to clusters really close to pi0 peak
          //   fhTowerDecayPhotonHit      [iSupMod1]->Fill(ieta1,iphi1);
          //   fhTowerDecayPhotonEnergy   [iSupMod1]->Fill(ieta1,iphi1,fMomentum1.E());
          //   fhTowerDecayPhotonAsymmetry[iSupMod1]->Fill(ieta1,iphi1,asym);
            
          //   fhTowerDecayPhotonHit      [iSupMod2]->Fill(ieta2,iphi2);
          //   fhTowerDecayPhotonEnergy   [iSupMod2]->Fill(ieta2,iphi2,fMomentum2.E());
          //   fhTowerDecayPhotonAsymmetry[iSupMod2]->Fill(ieta2,iphi2,asym);
            
          //   if(!mask1)fhTowerDecayPhotonHitMaskFrame[iSupMod1]->Fill(ieta1,iphi1);
          //   if(!mask2)fhTowerDecayPhotonHitMaskFrame[iSupMod2]->Fill(ieta2,iphi2);
          // }
        // }	else { //group cells
        //   //printf("Regroup N %d, eta1 %d, phi1 %d, eta2 %d, phi2 %d \n",fGroupNCells, ieta1, iphi1, ieta2, iphi2);
        //   for (Int_t i = -fGroupNCells; i < fGroupNCells+1; i++) {
        //     for (Int_t j = -fGroupNCells; j < fGroupNCells+1; j++) {              
        //       Int_t absId11 = fEMCALGeo->GetAbsCellIdFromCellIndexes(iSupMod1, iphi1+j, ieta1+i); 
        //       Int_t absId22 = fEMCALGeo->GetAbsCellIdFromCellIndexes(iSupMod2, iphi2+j, ieta2+i); 
              
        //       Bool_t ok1 = kFALSE;
        //       Bool_t ok2 = kFALSE;
              
        //       for(Int_t icell = 0; icell < c1->GetNCells(); icell++) {
        //         if(c1->GetCellsAbsId()[icell] == absId11) ok1=kTRUE;
        //       }
              
        //       for(Int_t icell = 0; icell < c2->GetNCells(); icell++) {
        //         if(c2->GetCellsAbsId()[icell] == absId22) ok2=kTRUE;
        //       }
              
        //       if(ok1) {
        //         fHmpi0[absId11]->Fill(invmass);
        //       }
                
        //       if(ok2 ) {
        //         fHmpi0[absId22]->Fill(invmass);
        //       }
        //     }
        //   }
        // }
        
      }
    }
  } // end of loop over EMCAL clusters
  
}


///_____________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::ProcessCells() {
  Int_t   nCells          = 0;
  nCells    =   fEMCALCells->GetNumberOfCells();
  if( nCells == 0 ) return;
  UShort_t nCellsAboveTh = 0;

  for(Long_t i=0; i<kMaxActiveCells_calib; i++){
    if(fEMCALCells->GetCellAmplitude(i) < fCellEmin ) continue;          // 50 MeV cut on cell energy

    fVBuffer_Cell_ID.push_back( static_cast<UShort_t>(i) );
    fVBuffer_Cell_t.push_back( static_cast<Short_t>(fEMCALCells->GetCellTime(i)*1e9) );
    fVBuffer_Cell_gain.push_back( fEMCALCells->GetCellHighGain(i) );

    fVBuffer_Cell_E.push_back( static_cast<UShort_t>(fEMCALCells->GetCellAmplitude(i)*1000) );

    nCellsAboveTh++;

    if( fIsMC ){
      fVBuffer_Cell_MCParticleID.push_back( static_cast<Short_t>(fEMCALCells->GetCellMCLabel(i)) );
    }
  }

  fVBuffer_NCells = nCellsAboveTh;
  return;
}


///_____________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::ProcessClusters() {
  Int_t nclus   = 0;
  UShort_t nClusAboveTh = 0;

  TClonesArray* clusters  = NULL;

  if( !fCorrTaskSetting.CompareTo("") ){
    clusters = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject("caloClusters"));
  } else {
    clusters = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject( Form("%sClustersBranch", fCorrTaskSetting.Data() ) ));
    if(!clusters)
      // AliFatal(Form("%sClustersBranch was not found! Check the correction framework settings!",fCorrTaskSetting.Data()));
      return;
  }
  nclus   = clusters->GetEntries();
  if(nclus == 0) return;   

  for(Long_t i=0; i<nclus; i++){
    if(fInputEvent->IsA()==AliAODEvent::Class()){
      AliVCluster* clus = NULL;
      clus  = dynamic_cast<AliAODCaloCluster*>(clusters->At(i));
      if( !clus ) continue;
      if( clus->IsPHOS()) continue;
      Float_t     clusPos[3];
      clus->GetPosition(clusPos);
      TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
      Float_t     etaCluster            = (Float_t) (clusterVector.Eta());
      Float_t     phiCluster            = (Float_t) (clusterVector.Phi());
      if( phiCluster < 0 ) phiCluster  += 2*TMath::Pi();


      fVBuffer_Cluster_LeadCellId.push_back( clus->GetCellAbsId(0) );
      fVBuffer_Cluster_t.push_back( static_cast<Short_t>(clus->GetTOF()*1.e9) );
      fVBuffer_Cluster_NCells.push_back( clus->GetNCells() );

      if( fSaveFullTree ){
        fVFBuffer_Cluster_E.push_back( clus->E() );
        fVFBuffer_Cluster_Eta.push_back( etaCluster );
        fVFBuffer_Cluster_Phi.push_back(phiCluster );
        fVBuffer_Cluster_M02.push_back( clus->GetM02() );
        fVBuffer_Cluster_X.push_back( clusPos[0] );
        fVBuffer_Cluster_Y.push_back( clusPos[1] );
        fVBuffer_Cluster_Z.push_back( clusPos[2] );
      } else {
        fVBuffer_Cluster_E.push_back( static_cast<UShort_t>(clus->E()*1000) );
        fVBuffer_Cluster_Eta.push_back( static_cast<Short_t>(etaCluster*1000) );
        fVBuffer_Cluster_Phi.push_back( static_cast<UShort_t>(phiCluster*1000) );
      }

      nClusAboveTh++;

      if ( fIsMC ){
        fVBuffer_TrueCluster_MCId.push_back( static_cast<Short_t>(clus->GetLabel()) );
      }

    } else  if(fInputEvent->IsA()==AliESDEvent::Class()){
      AliVCluster* clus = NULL;
      clus  = dynamic_cast<AliESDCaloCluster*>(clusters->At(i));
      if( !clus ) continue;
      if( clus->IsPHOS()) continue;

      Float_t     clusPos[3];
      clus->GetPosition(clusPos);
      TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
      Float_t     etaCluster            = (Float_t) (clusterVector.Eta());
      Float_t     phiCluster            = (Float_t) (clusterVector.Phi());
      if( phiCluster < 0 ) phiCluster  += 2*TMath::Pi();


      fVBuffer_Cluster_LeadCellId.push_back( clus->GetCellAbsId(0) );
      fVBuffer_Cluster_t.push_back( static_cast<Short_t>(clus->GetTOF()*1.e9) );
      fVBuffer_Cluster_NCells.push_back( clus->GetNCells() );

      if( fSaveFullTree ){
        fVFBuffer_Cluster_E.push_back( clus->E() );
        fVFBuffer_Cluster_Eta.push_back( etaCluster );
        fVFBuffer_Cluster_Phi.push_back(phiCluster );
        fVBuffer_Cluster_M02.push_back( clus->GetM02() );
        fVBuffer_Cluster_X.push_back( clusPos[0] );
        fVBuffer_Cluster_Y.push_back( clusPos[1] );
        fVBuffer_Cluster_Z.push_back( clusPos[2] );
      } else {
        fVBuffer_Cluster_E.push_back( static_cast<UShort_t>(clus->E()*1000) );
        fVBuffer_Cluster_Eta.push_back( static_cast<Short_t>(etaCluster*1000) );
        fVBuffer_Cluster_Phi.push_back( static_cast<UShort_t>(phiCluster*1000) );
      }
      
      nClusAboveTh++;

      if ( fIsMC ){
        fVBuffer_TrueCluster_MCId.push_back( static_cast<Short_t>(clus->GetLabel()) );
      }
    }
  }

  fBuffer_NClusters = nClusAboveTh;
  return;
}


///_____________________________________________________
UShort_t AliAnalysisTaskEMCALPi0CalibSelectionV2::GetPrimaryTracks(){
  Int_t                       fMinClsTPC = 70;  
  Double_t                    fChi2PerClsTPC = 5;   
  Int_t                       fMinClsITS = 0;  
  Double_t                    fEtaCut = 0.9;  
  Double_t                    fPtCut= 0.1;  

  UShort_t prim = 0;
  AliAODTrack *fCurrentTrack = NULL;

  for(Int_t t=0; t<fInputEvent->GetNumberOfTracks(); t++){
    fCurrentTrack = static_cast<AliAODTrack*>(fInputEvent->GetTrack(t));
    if( !fCurrentTrack ) continue;
    if( !fCurrentTrack->IsHybridGlobalConstrainedGlobal()) continue;
    if( fCurrentTrack->GetTPCNcls()<fMinClsTPC ) continue;
    if( fCurrentTrack->GetTPCchi2perCluster()>fChi2PerClsTPC ) continue;
    if( fCurrentTrack->GetITSNcls()<fMinClsITS ) continue;
    if( TMath::Abs(fCurrentTrack->Eta()) > fEtaCut ) continue;
    if( fCurrentTrack->Pt() < fPtCut ) continue;
    prim++;
  }

  return prim;
}

///
/// Main method, do the analysis per event:
/// * first, select the events;
/// * then, correct the clusters if needed;
/// * finally, fill the histograms per channel after recalibration.
//__________________________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::UserExec(Option_t* /* option */) {
  // Get the input event
  fInputEvent   = InputEvent();

  AliESDEvent*  esdev   = 0;
  AliAODEvent*  aodev   = 0;

  if( InputEvent()->IsA() == AliESDEvent::Class() ){
    esdev   = dynamic_cast<AliESDEvent*>(fInputEvent);
  } else if ( InputEvent()->IsA() == AliAODEvent::Class() ){
    aodev   = dynamic_cast<AliAODEvent*>(fInputEvent);
  }

  if( fIsMC>0 ){
    fMCEvent = MCEvent();
  }

  if(!fInputEvent) {
    AliWarning("Input event not available! Returning....");
    return;
  }

  if( fNContributorsCutEnabled ){             
    if( aodev->GetPrimaryVertex() != NULL ){
      if( aodev->GetPrimaryVertex()->GetNContributors() <= 0) {
        AliDebug(1,"No contributors to the vertex! Returning...");
        return;
      }
    }
  }

  if( !fEMCALInitialized ) {
    InitializeEMCAL();
    InitGeometryMatrices();
    if( !fEMCALInitialized ){
      AliWarning("EMCAL not initialized! Returning...");
      return;
    }
  }

  // Event selection
  isEMC=kFALSE;
  isDMC=kFALSE;

  if( !IsTriggerSelected() ){
    AliDebug(1,"Trigger not selected! Returning...");
    return;
  }
  
  if( !(isEMC) && !(isDMC) ){
    AliDebug(1,"Trigger not selected! Returning...");
    return;
  }


  // Centrality selection
  if ( fCheckCentrality ) {
    AliMultSelection* multSelection = (AliMultSelection * ) fInputEvent->FindListObject("MultSelection") ;
    if ( multSelection ) {
      Float_t cent = multSelection->GetMultiplicityPercentile(fCentralityClass, fCentWithEventSel);
      fhCentrality->Fill(cent);
      
      AliDebug(1,Form("Centrality %2.1f for class <%s>, event sel %d\n",cent,fCentralityClass.Data(),fCentWithEventSel));
      
      if ( cent < fCentMin || cent >= fCentMax ) {
        AliWarning("Outside given centrality, returnin....");
        return;
      }
      fhCentralitySelected->Fill(cent);
    }
  }

  
  AliDebug(1,Form("<<< %s: Event %d >>>",fInputEvent->GetName(), (Int_t)Entry()));
  
  // Get the primary vertex
  fInputEvent->GetPrimaryVertex()->GetXYZ(fVertex) ;
  AliDebug(1,Form("Vertex: (%.3f,%.3f,%.3f)",fVertex[0],fVertex[1],fVertex[2]));


  if(fSaveCells || fSaveClusters ){
    fBuffer_Event_VertexZ = fInputEvent->GetPrimaryVertex()->GetZ() * 100;
  }
  fBuffer_EventNPrimaryTracks = GetPrimaryTracks();
  
  fhNEvents->Fill(0); //Count the events to be analyzed

  //Get the list of clusters and cells
  fEMCALCells       = fInputEvent->GetEMCALCells();
  if( fEMCALCells == 0 ) return;

  if( fSaveHistos ) FillHistograms();  
  PostData(1,fOutputContainer); 

  if( fSaveCells ) ProcessCells();
  if( fSaveClusters ) ProcessClusters();

  if( fSaveCells || fSaveClusters) {
    Double_t fWeightJetJetMC = 1;
    if( fIsMC > 1 ){
      Float_t pthard = -1;
      Bool_t isMCJet = IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC, pthard, fInputEvent, fPeriodName );
      if (!isMCJet ) return;
      fBuffer_EventWeight = fWeightJetJetMC;
      fBuffer_ptHard      = pthard;
    }
    fCellTree->Fill();
    PostData(2,fCellTree);

    ResetBufferVectors();
  }
}

///
/// Print settings.
//_____________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::PrintInfo()
{
  printf("Cluster cuts: %2.2f < E < %2.2f GeV; number of cells > %d; Assymetry < %1.2f, pair time diff < %2.2f, %2.2f < t < %2.2f ns\n", 
         fEmin,fEmax, fMinNCells, fAsyCut, fDTimeCut,fTimeMin,fTimeMax) ;
  
  printf("Group %d cells\n", fGroupNCells) ;
    
  printf("Histograms: bins %d; energy range: %2.2f < E < %2.2f MeV\n",fNbins,fMinBin,fMaxBin) ;
  
  printf("OADB path        : %s\n",fOADBFilePath .Data());  
  printf("EMCAL Geometry name: < %s >, Load Matrices %d\n",fEMCALGeoName.Data(), fLoadMatrices) ;

  if ( fCheckCentrality )
    printf("Select centrality bin [%1.1f,%1.1f] for class <%s>, event sel. %d\n",fCentMin,fCentMax,fCentralityClass.Data(),fCentWithEventSel);
  
  if(fLoadMatrices) { 
    for(Int_t ism = 0; ism < AliEMCALGeoParams::fgkEMCALModules; ism++) {
      if(fMatrix[ism]) fMatrix[ism]->Print(); 
      }
  }
}

///
/// Accept cluster routine
//__________________________________________________________
Bool_t AliAnalysisTaskEMCALPi0CalibSelectionV2::AcceptCluster( AliVCluster* c1){
  // Exclude bad channels
  
  Float_t e1i = c1->E();   // cluster energy before correction   
  
  if (!c1->IsEMCAL()) 
    return kFALSE;
  if (c1->GetNCells() < fMinNCells) 
    return kFALSE;
  if (e1i > fEmax) 
    return kFALSE;
  if (e1i < fEmin) 
    return kFALSE;
  if (c1->GetIsExotic())
    return kFALSE;
  if ( fIsMC==0 && ( c1->GetTOF()*1.e9 > fTimeMax || c1->GetTOF()*1.e9 < fTimeMin))
    return kFALSE;
  if (c1->GetM02() < fL0min || c1->GetM02() > fL0max)
    return kFALSE;
  
  return kTRUE;

}

///
/// Check if trigger selected
///
//____________________________________________________________________
Bool_t AliAnalysisTaskEMCALPi0CalibSelectionV2::IsTriggerSelected(){

  UInt_t isSelected = 0;

  if(fTriggerName!="" && fIsMC==0 ) {
      AliInputEventHandler *fInputHandler=(AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()); 
      if (fInputHandler==NULL) return kFALSE;

      if( fInputHandler->GetEventSelection() &&  fTriggerName.CompareTo("CALO_ONLY") != 0) {
        if( fTriggerName.CompareTo("INT7")==0 ){
          fOfflineTriggerMask = AliVEvent::kINT7;
          isSelected = fOfflineTriggerMask & fInputHandler->IsEventSelected();
          if( !isSelected ) return kFALSE;
          isEMC = kTRUE; 
          isDMC = kTRUE;
        } else if ( fTriggerName.CompareTo("EMC7") == 0 ) {
          fOfflineTriggerMask = AliVEvent::kEMC7;

          TString firedTrigClass = fInputEvent->GetFiredTriggerClasses();
          if( (firedTrigClass.Contains("CEMC7-B-") || firedTrigClass.Contains("CDMC7-B-")) && (!firedTrigClass.Contains("EG1") && !firedTrigClass.Contains("EG2") && !firedTrigClass.Contains("DG1") && !firedTrigClass.Contains("DG2")  && !firedTrigClass.Contains("EJ1") && !firedTrigClass.Contains("EJ2") && !firedTrigClass.Contains("DJ1") && !firedTrigClass.Contains("DJ2") && !firedTrigClass.Contains("EGA") && !firedTrigClass.Contains("EJE") )) {
            if( firedTrigClass.Contains("CEMC7-B-") )
              isEMC = kTRUE;
            if( firedTrigClass.Contains("CDMC7-B-") )
              isDMC = kTRUE;
          } else return kFALSE;
        } else if ( fTriggerName.CompareTo("EG2") == 0 ) {
          fOfflineTriggerMask = AliVEvent::kEMCEGA;
          isSelected = fOfflineTriggerMask & fInputHandler->IsEventSelected();
          if( !isSelected ) return kFALSE;
          
          TString firedTrigClass = fInputEvent->GetFiredTriggerClasses();
          if( firedTrigClass.Contains("CEMC7EG2-B") || firedTrigClass.Contains("CDMC7DG2-B") && (!firedTrigClass.Contains("EG1") && !firedTrigClass.Contains("DG1")  && !firedTrigClass.Contains("EJ1") && !firedTrigClass.Contains("EJ2") && !firedTrigClass.Contains("DJ1") && !firedTrigClass.Contains("DJ2") && !firedTrigClass.Contains("EGA") && !firedTrigClass.Contains("EJE") )) {
            if( firedTrigClass.Contains("EMC7") )
              isEMC = kTRUE;
            if( firedTrigClass.Contains("DMC7") )
              isDMC = kTRUE;
          } else return kFALSE;
        } else if ( fTriggerName.CompareTo("EG1") == 0) {
          fOfflineTriggerMask = AliVEvent::kEMCEGA;
          isSelected = fOfflineTriggerMask & fInputHandler->IsEventSelected();
          if( !isSelected ) return kFALSE;
          
          TString firedTrigClass = fInputEvent->GetFiredTriggerClasses();
          if( firedTrigClass.Contains("CEMC7EG1-B") || firedTrigClass.Contains("CDMC7DG1-B")&& (!firedTrigClass.Contains("EJ1") && !firedTrigClass.Contains("EJ2") && !firedTrigClass.Contains("DJ1") && !firedTrigClass.Contains("DJ2") && !firedTrigClass.Contains("EGA") && !firedTrigClass.Contains("EJE") )) {
            if( firedTrigClass.Contains("EMC7") )
              isEMC = kTRUE;
            if( firedTrigClass.Contains("DMC7") )
              isDMC = kTRUE;
          } else return kFALSE;
        } else if ( fTriggerName.CompareTo("EGA") == 0 ) {
          fOfflineTriggerMask = AliVEvent::kEMCEGA;
          isSelected = fOfflineTriggerMask & fInputHandler->IsEventSelected();
          if( !isSelected ) return kFALSE;
          
          TString firedTrigClass = fInputEvent->GetFiredTriggerClasses();
          if( firedTrigClass.Contains("EGA") && (!firedTrigClass.Contains("EJ1") && !firedTrigClass.Contains("EJ2") && !firedTrigClass.Contains("DJ1") && !firedTrigClass.Contains("DJ2") && !firedTrigClass.Contains("EJE") )) {
              isEMC = kTRUE;
              isDMC = kFALSE;
          } else return kFALSE;
        }else if( fTriggerName.CompareTo("EMCAL") == 0 ){            // EMC7+EG1+EG2
          fOfflineTriggerMask = AliVEvent::kEMC7 & AliVEvent::kEMCEGA;
          TString firedTrigClass = fInputEvent->GetFiredTriggerClasses();

          if( firedTrigClass.Contains("CEMC7") || firedTrigClass.Contains("CDMC7") ) {
            if( firedTrigClass.Contains("CEMC7-B-") )
              isEMC = kTRUE;
            if( firedTrigClass.Contains("CDMC7-B-") )
              isDMC = kTRUE;
          } else return kFALSE;
        }
        return kTRUE;
      } else if( fTriggerName.CompareTo("CALO_ONLY") == 0 ) {
        TString firedTrigClass = fInputEvent->GetFiredTriggerClasses();
        // std::cout << fInputEvent->GetFiredTriggerClasses().Data() << std::endl;
        if( firedTrigClass.Contains("CEMC7-B-NOPF-EMCAL") || firedTrigClass.Contains("CDMC7-B-NOPF-EMCAL") ) {
          if( firedTrigClass.Contains("EMC7") )
            isEMC = kTRUE;
          if( firedTrigClass.Contains("DMC7") )
            isDMC = kTRUE;
        } else return kFALSE;
        // std::cout << "found it" << std::endl;
        return kTRUE;
      }
  } else if (fIsMC){
    AliInputEventHandler *fInputHandler=(AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()); 
      if (fInputHandler==NULL) return kFALSE;
      if( fInputHandler->GetEventSelection() ) {
        fOfflineTriggerMask = AliVEvent::kAny;
        if( !(fOfflineTriggerMask) ) return kFALSE;
        isEMC = kTRUE;
        isDMC = kTRUE;
      }
      return kTRUE;
  }

  std::cout << "No trigger selected" << std::endl;
  return kFALSE;

}

///
/// Return JJ MC weights
///
//____________________________________________________________________
Bool_t AliAnalysisTaskEMCALPi0CalibSelectionV2::IsJetJetMCEventAccepted( AliMCEvent *mcEvent, Double_t& weight, Float_t& pthard, AliVEvent* event, TString fPeriodName ){
  if(fPeriodName.CompareTo("LHC16P1JJ") != 0 && fPeriodName.CompareTo("LHC17P1JJ") != 0 && fPeriodName.CompareTo("LHC18P1JJ") != 0 && !fPeriodName.Contains("LHC19d3") ){
    std::cout << "Weights not implemented for given period" << std::endl;
    weight = 1;
    return kFALSE;
  }

  AliGenCocktailEventHeader *cHeader   = 0x0;
  Bool_t headerFound                   = kFALSE;
  weight                               = -1;
  Double_t fMaxPtJetMC                 = 0;
  Bool_t eventAccepted = kTRUE;

  if(mcEvent){
    cHeader           = dynamic_cast<AliGenCocktailEventHeader*>(mcEvent->GenEventHeader());
    if(cHeader) headerFound   = kTRUE;
  } else {
    //no mcEvent available -> not running on MC
    weight = 1;
    return kTRUE;
  }

  if(headerFound) {
    TList *genHeaders         = 0x0;
    if(cHeader) genHeaders    = cHeader->GetHeaders();
    AliGenEventHeader* gh     = 0;

    for(Int_t i = 0; i<genHeaders->GetEntries();i++){
      gh = (AliGenEventHeader*)genHeaders->At(i);
      TString GeneratorName = gh->GetName();
      Int_t nTriggerJets = dynamic_cast<AliGenPythiaEventHeader*>(gh)->NTriggerJets();
      Float_t ptHard = dynamic_cast<AliGenPythiaEventHeader*>(gh)->GetPtHard();
      Float_t tmpjet[]={0,0,0,0};

      for(Int_t ijet = 0; ijet< nTriggerJets; ijet++){
        dynamic_cast<AliGenPythiaEventHeader*>(gh)->TriggerJet(ijet, tmpjet);
        TParticle jet(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
        //Compare jet pT and pt Hard
        if(jet.Pt() > 2.5 * ptHard){
          eventAccepted= kFALSE;
        }
        //set highest jet pT
        if (jet.Pt() > fMaxPtJetMC) fMaxPtJetMC = jet.Pt();
      }

      if (mcEvent){
        for(Long_t i = 0; i < mcEvent->GetNumberOfPrimaries(); i++) {
          AliMCParticle* particle = (AliMCParticle*) mcEvent->GetTrack(i);
          if (!particle) continue;
          // if (TMath::Abs(particle->GetPdgCode()) == 111 || TMath::Abs(particle->GetPdgCode()) == 221){
              if (particle->Pt() > 1.5*ptHard && TMath::Abs(particle->PdgCode()) > 21){
                eventAccepted= kFALSE;
              }
          // }
        }
      }

      Double_t ptHardBinRanges[21]  = { 5,  7,  9, 12, 16,
                                              21, 28, 36, 45, 57,
                                              70, 85, 99, 115, 132,
                                              150, 169, 190, 212, 235,
                                              1000000};
      Double_t weightsBins[20]      = { 43.8654,  13.6215, 6.79856, 2.67526, 0.978794,
                                        0.390797,  0.127769, 0.0465714, 0.0206173, 0.00750282,
                                        0.00318773,  0.00122533, 0.000644385, 0.000321225,  0.00016846,
                                        9.18305e-05, 5.33507e-05, 3.00677e-05, 1.74608e-05, 2.80823e-05};
      Int_t bin = 0;
      if(ptHard >= ptHardBinRanges[0]){
        while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
        if (bin < 20) weight = weightsBins[bin];
      }
    }
  } else {
    AliGenEventHeader * eventHeader = mcEvent->GenEventHeader();
    TString eventHeaderName     = eventHeader->ClassName();
    Bool_t eventAccepted = kFALSE;
    if (eventHeaderName.CompareTo("AliGenPythiaEventHeader") == 0 || eventHeaderName.Contains("Pythia8Jets")){
      eventAccepted = kTRUE;
    }else { //special case for pythia8jets embedded in EPOSLHC for AODs
      if(event->IsA()==AliAODEvent::Class()){
        AliAODMCHeader *mch = NULL;
        AliAODEvent * aod = dynamic_cast<AliAODEvent*> (event);
        if(aod){
          mch = dynamic_cast<AliAODMCHeader*>(aod->FindListObject("mcHeader"));
          if ( mch ){
            Int_t nGenerators = mch->GetNCocktailHeaders();
            if ( nGenerators > 0  ){
              for(Int_t igen = 0; igen < nGenerators; igen++)
              {
                AliGenEventHeader * eventHeaderGen = mch->GetCocktailHeader(igen) ;
                TString name = eventHeaderGen->GetName();
                if (name.CompareTo("AliGenPythiaEventHeader") == 0 || name.Contains("Pythia8Jets") || name.Contains("Pythia8GammaJet")){
                  eventAccepted = kTRUE;
                  eventHeader = eventHeaderGen;
                }
              }
            }
          }
        }
      }
    }
    if(eventAccepted){
      Int_t nTriggerJets =  dynamic_cast<AliGenPythiaEventHeader*>(eventHeader)->NTriggerJets();
      Float_t ptHard = dynamic_cast<AliGenPythiaEventHeader*>(eventHeader)->GetPtHard();
      pthard = ptHard;
      Float_t tmpjet[]={0,0,0,0};

        for(Int_t ijet = 0; ijet< nTriggerJets; ijet++){
          dynamic_cast<AliGenPythiaEventHeader*>(eventHeader)->TriggerJet(ijet, tmpjet);
          TParticle jet(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
          //Compare jet pT and pt Hard
          if(jet.Pt() > 2.5 * ptHard){
            eventAccepted= kFALSE;
          }
          //set highest jet pT
          if (jet.Pt() > fMaxPtJetMC){
            fMaxPtJetMC = jet.Pt();
          }
        }

      if (mcEvent){
        for(Long_t i = 0; i < mcEvent->GetNumberOfPrimaries(); i++) {
          // TParticle* particle = (TParticle *)mcEvent->Particle(i);
          AliMCParticle* particle = (AliMCParticle*) mcEvent->GetTrack(i);
          if (!particle) continue;
          // if (TMath::Abs(particle->GetPdgCode()) == 111 || TMath::Abs(particle->GetPdgCode()) == 221){
              if (particle->Pt() > 1.5*ptHard && TMath::Abs(particle->PdgCode()) > 21){
                eventAccepted= kFALSE;
              }
          // }
        }
      }

        Double_t ptHardBinRanges[21]  = { 5,  7,  9, 12, 16,
                                          21, 28, 36, 45, 57,
                                          70, 85, 99, 115, 132,
                                          150, 169, 190, 212, 235,
                                          1000000};
        Double_t weightsBins[20]      = { 43.8654,  13.6215, 6.79856, 2.67526, 0.978794,
                                          0.390797,  0.127769, 0.0465714, 0.0206173, 0.00750282,
                                          0.00318773,  0.00122533, 0.000644385, 0.000321225,  0.00016846,
                                          9.18305e-05, 5.33507e-05, 3.00677e-05, 1.74608e-05, 2.80823e-05};

        Int_t bin = 0;
        if(ptHard >= ptHardBinRanges[0]){
          while (!((ptHard< ptHardBinRanges[bin+1] && ptHard > ptHardBinRanges[bin]) || (ptHard == ptHardBinRanges[bin]) ) )bin++;
          if (bin < 20) weight = weightsBins[bin];
        }
      }

  }
  if(weight == -1) return kFALSE;
  else return eventAccepted;
}

///
/// Set the total number of columns to be masked in the analysis
/// \param n: number of columns
//_____________________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::SetNMaskCellColumns(Int_t n)
{
    if(n > fNMaskCellColumns) {
        delete [] fMaskCellColumns ;
        fMaskCellColumns = new Int_t[n] ;
    }
    
    fNMaskCellColumns = n ;
}

///
/// Define which column must be masked and its position in the array of columns
/// \param ipos: position in the list of columns.
/// \param icol: column to be masked in EMCal mapping coordinate.
//___________________________________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::SetMaskCellColumn(Int_t ipos, Int_t icol)
{
    if(ipos < fNMaskCellColumns) fMaskCellColumns[ipos] = icol            ;
    else AliWarning("Mask column not set, position larger than allocated set size first") ;
}

///
/// Check if cell is in one of the regions where we have significant amount of material in front of EMCAL.
/// \return True if this cell is in one problematic region
/// \param iSM: supermodule number of the cell.
/// \param ieta: column index of the cell.
//______________________________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPi0CalibSelectionV2::MaskFrameCluster(Int_t iSM, Int_t ieta) const
{
  Int_t icol = ieta;
  if(iSM%2) icol+=48; // Impair SM, shift index [0-47] to [48-96]
  
  if (fNMaskCellColumns && fMaskCellColumns) {
    for (Int_t imask = 0; imask < fNMaskCellColumns; imask++) {
      if(icol==fMaskCellColumns[imask]) return kTRUE;
    }
  }
  
  return kFALSE;
}

Int_t AliAnalysisTaskEMCALPi0CalibSelectionV2::RecalculateRow( Int_t row, Int_t nSupMod ){
    row += nSupMod/2 * (24+1);
    return row;
}

Int_t AliAnalysisTaskEMCALPi0CalibSelectionV2::RecalculateColumn( Int_t col, Int_t nSupMod) {
    if( nSupMod < 12 ){
        col += nSupMod%2 * 48;
    } else {
        col += nSupMod%2 * (48+1);
    }
    return col;
}

///_____________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::ResetBufferVectors() {

  fBuffer_NClusters = 0;

  fVBuffer_Cell_ID.clear();
  fVBuffer_Cell_E.clear();
  fVBuffer_Cell_t.clear();
  fVBuffer_Cell_gain.clear();
  fVBuffer_Cell_MCParticleID.clear();

  fVBuffer_Cell_ID.resize(0);
  fVBuffer_Cell_E.resize(0);
  fVBuffer_Cell_t.resize(0);
  fVBuffer_Cell_gain.resize(0);
  fVBuffer_Cell_MCParticleID.resize(0);

  fVBuffer_Cluster_E.clear();
  fVBuffer_Cluster_Eta.clear();
  fVBuffer_Cluster_Phi.clear();
  fVFBuffer_Cluster_E.clear();
  fVFBuffer_Cluster_Eta.clear();
  fVFBuffer_Cluster_Phi.clear();
  fVBuffer_Cluster_t.clear();
  fVBuffer_TrueCluster_MCId.clear();

  fVBuffer_Cluster_E.resize(0);
  fVBuffer_Cluster_Eta.resize(0);
  fVBuffer_Cluster_Phi.resize(0);
  fVFBuffer_Cluster_E.resize(0);
  fVFBuffer_Cluster_Eta.resize(0);
  fVFBuffer_Cluster_Phi.resize(0);
  fVBuffer_Cluster_t.resize(0);
  fVBuffer_TrueCluster_MCId.resize(0);

  fVBuffer_Cluster_LeadCellId.clear(); fVBuffer_Cluster_LeadCellId.resize(0);

  fVBuffer_Cluster_M02.clear(); fVBuffer_Cluster_M02.resize(0);
  fVBuffer_Cluster_X.clear();   fVBuffer_Cluster_X.resize(0);
  fVBuffer_Cluster_Y.clear();   fVBuffer_Cluster_Y.resize(0);
  fVBuffer_Cluster_Z.clear();   fVBuffer_Cluster_Z.resize(0);
  fVBuffer_Cluster_NCells.clear(); fVBuffer_Cluster_NCells.resize(0);
}

///
/// Create cuts/param objects and publish to slot. Comment out for the moment.
//______________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelectionV2::Terminate(Option_t*)
{
  AliDebug(1,"Not implemented");

}

