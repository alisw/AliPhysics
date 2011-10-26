/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  EMCAL tender, apply corrections to EMCAl clusters                        //
//  and do track matching.                                                   //                                                                           
//  Author: Deepa Thomas (Utrecht University)                                // 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TInterpreter.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliAnalysisManager.h>
#include <AliTender.h>
#include "AliOADBContainer.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliMagF.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALTenderSupply.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALClusterizer.h"
#include "AliEMCALRecParam.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALAfterBurnerUF.h"
#include "AliEMCALClusterizerNxN.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALClusterizerv2.h"
#include "AliEMCALDigit.h"

ClassImp(AliEMCALTenderSupply)

AliEMCALTenderSupply::AliEMCALTenderSupply() :
AliTenderSupply()
,fEMCALGeo(0x0)
,fEMCALGeoName("EMCAL_COMPLETEV1")
,fEMCALRecoUtils(0)
,fConfigName("")
,fDebugLevel(0)
,fNonLinearFunc(AliEMCALRecoUtils::kNoCorrection) 
,fNonLinearThreshold(30)      	
,fReCalibCluster(kFALSE)	
,fReCalibCell(kFALSE)	
,fUpdateCell(kFALSE)  
,fRecalClusPos(kFALSE)
,fFiducial(kFALSE) 
,fNCellsFromEMCALBorder(1)	
,fRecalDistToBadChannels(kFALSE)	
,fInputTree(0)	
,fInputFile(0)
,fFilepass(0) 
,fMass(0.139)
,fStep(50)
,fCutEtaPhiSum(kTRUE)
,fCutEtaPhiSeparate(kFALSE)
,fRcut(0.05)	
,fEtacut(0.025)	
,fPhicut(0.05)	
,fBasePath("")
,fReClusterize(kFALSE)
,fClusterizer(0)
,fGeomMatrixSet(kFALSE)
,fLoadGeomMatrices(kFALSE)
,fRecParam(0x0)
,fRecParamSet(kFALSE)
,fOCDBpath(" ")
,fUnfolder(0)
,fDigitsArr(0)
,fClusterArr(0)
,fMisalignSurvey(kdefault)	
{
  // Default constructor.
  for(Int_t i = 0; i < 10; i++) fEMCALMatrix[i] = 0 ;
}

//_____________________________________________________
AliEMCALTenderSupply::AliEMCALTenderSupply(const char *name, const AliTender *tender) :
AliTenderSupply(name,tender)
,fEMCALGeo(0x0)
,fEMCALGeoName("EMCAL_COMPLETEV1")
,fEMCALRecoUtils(0)
,fConfigName("") 
,fDebugLevel(0)
,fNonLinearFunc(AliEMCALRecoUtils::kNoCorrection)      	
,fNonLinearThreshold(30)      	
,fReCalibCluster(kFALSE)	
,fReCalibCell(kFALSE)	
,fUpdateCell(kFALSE)  
,fRecalClusPos(kFALSE)
,fFiducial(kFALSE) 
,fNCellsFromEMCALBorder(1)	
,fRecalDistToBadChannels(kFALSE)	
,fInputTree(0)	
,fInputFile(0)
,fFilepass(0)
,fMass(0.139)
,fStep(50)
,fCutEtaPhiSum(kTRUE)
,fCutEtaPhiSeparate(kFALSE)
,fRcut(0.05)
,fEtacut(0.025)	
,fPhicut(0.05)	
,fBasePath("")
,fReClusterize(kFALSE)
,fClusterizer(0)
,fGeomMatrixSet(kFALSE)
,fLoadGeomMatrices(kFALSE)
,fRecParam(0x0)
,fRecParamSet(kFALSE)
,fOCDBpath(" ")
,fUnfolder(0)
,fDigitsArr(0)
,fClusterArr(0)
,fMisalignSurvey(kdefault)	
{
  // Named constructor
  
  for(Int_t i = 0; i < 10; i++) fEMCALMatrix[i] = 0 ;
  
  fEMCALRecoUtils  = new AliEMCALRecoUtils();
  fRecParam        = new AliEMCALRecParam;
  fDigitsArr       = new TClonesArray("AliEMCALDigit",1000);
}

//_____________________________________________________
AliEMCALTenderSupply::~AliEMCALTenderSupply()
{
  //Destructor
  
  delete fEMCALRecoUtils;
  delete fClusterizer;
  delete fUnfolder;
  if (fDigitsArr){
    fDigitsArr->Clear("C");
    delete fDigitsArr; 
  }
}

//_____________________________________________________
void AliEMCALTenderSupply::Init()
{
  // Initialise EMCAL tender.

  if (fDebugLevel>0) 
    AliInfo("Init EMCAL Tender supply");	
  
  if(fConfigName.Length()>0 && gROOT->LoadMacro(fConfigName) >=0) {
    AliDebug(1, Form("Loading settings from macro %s", fConfigName.Data()));
    AliEMCALTenderSupply *tender = (AliEMCALTenderSupply*)gInterpreter->ProcessLine("ConfigEMCALTenderSupply()");
    fDebugLevel             = tender->fDebugLevel;
    fEMCALGeoName           = tender->fEMCALGeoName; 
    delete fEMCALRecoUtils;
    fEMCALRecoUtils         = tender->fEMCALRecoUtils; 
    fConfigName             = tender->fConfigName;
    fNonLinearFunc          = tender->fNonLinearFunc;
    fNonLinearThreshold     = tender->fNonLinearThreshold;
    fReCalibCluster         = tender->fReCalibCluster;
    fReCalibCell            = tender->fReCalibCell;
    fRecalClusPos           = tender->fRecalClusPos;
    fFiducial	            = tender->fFiducial;
    fNCellsFromEMCALBorder  = tender->fNCellsFromEMCALBorder;
    fRecalDistToBadChannels = tender->fRecalDistToBadChannels;    
    fMass                   = tender->fMass;
    fStep                   = tender->fStep;
    fRcut                   = tender->fRcut;
    fEtacut                 = tender->fEtacut;
    fPhicut                 = tender->fPhicut;
    fReClusterize           = tender->fReClusterize;
    fLoadGeomMatrices       = tender->fLoadGeomMatrices;
    fRecParam               = tender->fRecParam;
    fOCDBpath               = tender->fOCDBpath;
    for(Int_t i = 0; i < 10; i++) 
      fEMCALMatrix[i] = tender->fEMCALMatrix[i] ;
  }
  
  // Init goemetry	
  fEMCALGeo = AliEMCALGeometry::GetInstance(fEMCALGeoName) ;
  
  // Initialising Non linearity parameters
  fEMCALRecoUtils->SetNonLinearityThreshold(fNonLinearThreshold);
  fEMCALRecoUtils->SetNonLinearityFunction(fNonLinearFunc);
  
  // Setting mass, step size and residual cut
  fEMCALRecoUtils->SetMass(fMass);
  fEMCALRecoUtils->SetStep(fStep);
  if(fCutEtaPhiSum){ 
    fEMCALRecoUtils->SwitchOnCutEtaPhiSum(); 
    fEMCALRecoUtils->SetCutR(fRcut);
  } else if (fCutEtaPhiSeparate) {
    fEMCALRecoUtils->SwitchOnCutEtaPhiSeparate();
    fEMCALRecoUtils->SetCutEta(fEtacut);
    fEMCALRecoUtils->SetCutPhi(fPhicut);
  }
}

//_____________________________________________________
void AliEMCALTenderSupply::ProcessEvent()
{
  // Event loop.
  
  AliESDEvent *event = fTender->GetEvent();
  if (!event) {
    AliError("ESD event ptr = 0, returning");
    return;
  }
/*
  if(!fRecParamSet)
   {
    if(!fRecParam)
    {
      InitRecParam();
      fRecParamSet = kTRUE;
    }
   }
 */ 
  // Initialising parameters once per run number
  if(fTender->RunChanged()){ 
    GetPass();
    // Init bad channels
    Int_t fInitBC=InitBadChannels();
    
    if (fInitBC==0)
      {
	AliError("InitBadChannels returned false, returning");
	return;
      }
    if(fInitBC>1)
      {
	AliInfo(Form("No external hot channel set: %d - %s", event->GetRunNumber(), fFilepass.Data()));
      }

    if (fReCalibCluster || fReCalibCell || fUpdateCell) { 
      Int_t fInitRecalib=InitRecalib();
      if (fInitRecalib==0)
	{
	  AliError("InitRecalib returned false, returning");
	  return;
	}
      if(fInitRecalib >1)
	{
	  AliInfo(Form("No recalibration available: %d - %s", event->GetRunNumber(), fFilepass.Data()));
	  fReCalibCell=kFALSE;
	  fReCalibCluster=kFALSE;
	}
    }

    if (fRecalClusPos || fReClusterize || fUpdateCell) { 
      if (!InitMisalignMatrix()) { 
        AliError("InitMisalignmentMatrix returned false, returning");
        return;
      }
    }

    if (fReClusterize || fUpdateCell) {
      if (!InitClusterization()) {
        AliError("InitClusterization returned false, returning");
        return;
      }
    }
    
    if(fDebugLevel>1) 
      fEMCALRecoUtils->Print("");
      
  }

  // Test if cells present
  AliESDCaloCells *cells= event->GetEMCALCells();
  if (cells->GetNumberOfCells()<=0) {
    AliWarning(Form("Number of EMCAL cells = %d, returning", cells->GetNumberOfCells()));
    return;
  }
  
  // Recalibrate cells
  if (fReCalibCell || fUpdateCell)
    { 
      RecalibrateCells();
      fReCalibCluster=kFALSE;
    }
  if(fDebugLevel>2)
     AliInfo(Form("Re-calibrate cluster %d\n",fReCalibCluster));

  if(fUpdateCell)
    {
      printf("Update cells\n");
      UpdateCells();
      fReClusterize=kTRUE;
    }

    // Reclusterize
  if(fReClusterize)
    {

      FillDigitsArray();
      Clusterize();
      UpdateClusters();
    }

  // Store good clusters
  TClonesArray *clusArr = dynamic_cast<TClonesArray*>(event->FindListObject("caloClusters"));
  if (!clusArr) 
    clusArr = dynamic_cast<TClonesArray*>(event->FindListObject("CaloClusters"));
  if (!clusArr) {
    AliWarning(Form("No cluster array, number of cells in event = %d, returning", cells->GetNumberOfCells()));
    return;
  }
  Int_t nclusters = clusArr->GetEntriesFast();
  
  for (Int_t icluster=0; icluster < nclusters; ++icluster) { 
    AliVCluster *clust = static_cast<AliVCluster*>(clusArr->At(icluster));
    if (!clust) 
      continue;
    if  (!clust->IsEMCAL()) 
      continue;
    
    if (fEMCALRecoUtils->ClusterContainsBadChannel(fEMCALGeo, clust->GetCellsAbsId(), clust->GetNCells())) {
      delete clusArr->RemoveAt(icluster);
      continue; //todo is it really needed to remove it? Or should we flag it?
    }
    
    if (fFiducial){
      if (!fEMCALRecoUtils->CheckCellFiducialRegion(fEMCALGeo, clust, cells)){
        delete clusArr->RemoveAt(icluster);
        continue; // todo it would be nice to store the distance
      }
    }
    
    fEMCALRecoUtils->CorrectClusterEnergyLinearity(clust);
    if(fRecalDistToBadChannels) 
      fEMCALRecoUtils->RecalculateClusterDistanceToBadChannel(fEMCALGeo, cells, clust);  
    if(fReCalibCluster) 
      fEMCALRecoUtils->RecalibrateClusterEnergy(fEMCALGeo, clust, cells);
    if(fRecalClusPos) 
      fEMCALRecoUtils->RecalculateClusterPosition(fEMCALGeo, cells, clust);
  }
  clusArr->Compress();
  
  // Track matching
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    const AliESDRun *erun = event->GetESDRun();
    AliMagF *field = AliMagF::CreateFieldMap(erun->GetCurrentL3(),
                                             erun->GetCurrentDip(),
                                             AliMagF::kConvLHC,
                                             kFALSE,
                                             erun->GetBeamEnergy(),
                                             erun->GetBeamType());
    TGeoGlobalMagField::Instance()->SetField(field);
  }
  
  fEMCALRecoUtils->FindMatches(event,0x0,fEMCALGeo);
 
  SetClusterMatchedToTrack(event);
  SetTracksMatchedToCluster(event);
  
}

//_____________________________________________________
void AliEMCALTenderSupply::SetClusterMatchedToTrack(AliESDEvent *event)
{
  // Checks if tracks are matched to EMC clusters and set the matched EMCAL cluster index to ESD track. 
  
  Int_t nTracks = event->GetNumberOfTracks();
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    AliESDtrack* track = event->GetTrack(iTrack);
    if (!track) {
      AliWarning(Form("Could not receive track %d", iTrack));
      continue;
    }
    Int_t matchClusIndex = fEMCALRecoUtils->GetMatchedClusterIndex(iTrack);		   
    track->SetEMCALcluster(matchClusIndex); //sets -1 if track not matched within residual
    if(matchClusIndex != -1) 
      track->SetStatus(AliESDtrack::kEMCALmatch);
    else
      track->ResetStatus(AliESDtrack::kEMCALmatch);
  }
  if (fDebugLevel>2) 
    AliInfo("Track matched to closest cluster");	
}

//_____________________________________________________
void AliEMCALTenderSupply::SetTracksMatchedToCluster(AliESDEvent *event)
{
  // Checks if EMC clusters are matched to ESD track.
  // Adds track indexes of all the tracks matched to a cluster withing residuals in ESDCalocluster.
  
  for (Int_t iClus=0; iClus < event->GetNumberOfCaloClusters(); ++iClus) {
    AliESDCaloCluster *cluster = event->GetCaloCluster(iClus);
    if (!cluster->IsEMCAL()) 
      continue;
    
    Int_t nTracks = event->GetNumberOfTracks();
    TArrayI arrayTrackMatched(nTracks);
    
    // Get the closest track matched to the cluster
    Int_t nMatched = 0;
    Int_t matchTrackIndex = fEMCALRecoUtils->GetMatchedTrackIndex(iClus);
    if (matchTrackIndex != -1) {
      arrayTrackMatched[nMatched] = matchTrackIndex;
      nMatched++;
    }
    
    // Get all other tracks matched to the cluster
    for(Int_t iTrk=0; iTrk<nTracks; ++iTrk) {
      AliESDtrack* track = event->GetTrack(iTrk);
      if(iTrk == matchTrackIndex) continue;
      if(track->GetEMCALcluster() == iClus){
        arrayTrackMatched[nMatched] = iTrk;
        ++nMatched;
      }
    }
    
    arrayTrackMatched.Set(nMatched);
    cluster->AddTracksMatched(arrayTrackMatched);
    
    Float_t eta= -999, phi = -999;
    if (matchTrackIndex != -1) 
      fEMCALRecoUtils->GetMatchedResiduals(iClus, eta, phi);
    cluster->SetTrackDistance(phi, eta);
  }
  
  if (fDebugLevel>2) 
    AliInfo("Cluster matched to tracks");	
}

//_____________________________________________________
Bool_t AliEMCALTenderSupply::InitMisalignMatrix()
{
  // Initialising misalignment matrices
  
  AliESDEvent *event = fTender->GetEvent();
  if (!event) 
    return kFALSE;
  
  if (fGeomMatrixSet) {
    AliInfo("Misalignment matrix already set");	
    return kTRUE;
  }
  
  if (fDebugLevel>0) 
    AliInfo("Initialising misalignment matrix");	
  
  fEMCALRecoUtils->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);
  
  if (fLoadGeomMatrices) {
    for(Int_t mod=0; mod < fEMCALGeo->GetNumberOfSuperModules(); ++mod) {
      if (fEMCALMatrix[mod]){
        if(fDebugLevel > 2) 
          fEMCALMatrix[mod]->Print();
        fEMCALGeo->SetMisalMatrix(fEMCALMatrix[mod],mod);  
      }
    }
    fGeomMatrixSet = kTRUE;
    return kTRUE;
  }
  
  Int_t runGM = event->GetRunNumber();
  TObjArray *mobj = 0;

 if(fMisalignSurvey == kdefault){ //take default alignment corresponding to run no
    AliOADBContainer emcalgeoCont(Form("emcal"));
    emcalgeoCont.InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALlocal2master.root",Form("AliEMCALgeo"));
    mobj=(TObjArray*)emcalgeoCont.GetObject(runGM,"EmcalMatrices");
 }

 if(fMisalignSurvey == kSurveybyS){ //take alignment at sector level
  if (runGM <= 140000) { //2010 data
    AliOADBContainer emcalgeoCont(Form("emcal2010"));
    emcalgeoCont.InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALlocal2master.root",Form("AliEMCALgeo"));
    mobj=(TObjArray*)emcalgeoCont.GetObject(100,"survey10");
    
  } else if (runGM>140000) { // 2011 LHC11a pass1 data
    AliOADBContainer emcalgeoCont(Form("emcal2011"));
    emcalgeoCont.InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALlocal2master.root",Form("AliEMCALgeo"));
    mobj=(TObjArray*)emcalgeoCont.GetObject(100,"survey11byS");			
  }
 }

 if(fMisalignSurvey == kSurveybyS){ //take alignment at module level
  if (runGM <= 140000) { //2010 data
    AliOADBContainer emcalgeoCont(Form("emcal2010"));
    emcalgeoCont.InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALlocal2master.root",Form("AliEMCALgeo"));
    mobj=(TObjArray*)emcalgeoCont.GetObject(100,"survey10");
    
  } else if (runGM>140000) { // 2011 LHC11a pass1 data
    AliOADBContainer emcalgeoCont(Form("emcal2011"));
    emcalgeoCont.InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALlocal2master.root",Form("AliEMCALgeo"));
    mobj=(TObjArray*)emcalgeoCont.GetObject(100,"survey11byM");			
  }
 }

 for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
   fEMCALMatrix[mod] = (TGeoHMatrix*) mobj->At(mod);
   fEMCALGeo->SetMisalMatrix(fEMCALMatrix[mod],mod); 
   fEMCALMatrix[mod]->Print();
 }
  return kTRUE;
}

//_____________________________________________________
Int_t AliEMCALTenderSupply::InitBadChannels()
{
  // Initialising bad channel maps
  AliESDEvent *event = fTender->GetEvent();
  if (!event) 
    return 0;
  
  if (fDebugLevel>0) 
    AliInfo("Initialising Bad channel map");
  
  if (fFiducial){
    fEMCALRecoUtils->SetNumberOfCellsFromEMCALBorder(fNCellsFromEMCALBorder);
    fEMCALRecoUtils->SwitchOnNoFiducialBorderInEMCALEta0();
  }
  
  fEMCALRecoUtils->SwitchOnBadChannelsRemoval();
  if (fRecalDistToBadChannels) 
    fEMCALRecoUtils->SwitchOnDistToBadChannelRecalculation();
  
  Int_t runBC = event->GetRunNumber();
  
  AliOADBContainer *contBC=new AliOADBContainer("");
  if(fBasePath!=""){ //if fBasePath specified in the ->SetBasePath()
    if (fDebugLevel>0) AliInfo(Form("Loading Bad Channels OADB from given path %s",fBasePath.Data()));
    TFile *fbad=new TFile(Form("%s/EMCALBadChannels.root",fBasePath.Data()),"read");
    if(fbad->IsZombie()){
      AliFatal(Form("EMCALBadChannels.root was not found in the path provided: %s, aborting",fBasePath.Data()));
      return 0;
    }  
    if(fbad) delete fbad;
    contBC->InitFromFile(Form("%s/EMCALBadChannels.root",fBasePath.Data()),"AliEMCALBadChannels");    
  }
    else { // Else choose the one in the $ALICE_ROOT directory
      if (fDebugLevel>0) AliInfo("Loading Bad Channels OADB from $ALICE_ROOT/OADB/EMCAL");
      TFile *fbad=new TFile("$ALICE_ROOT/OADB/EMCAL/EMCALBadChannels.root","read");
      if(fbad->IsZombie()){
	AliFatal("$ALICE_ROOT/OADB/EMCAL/EMCALBadChannels.root was not found, aborting");
	return 0;
      }  
      if(fbad) delete fbad;
      contBC->InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALBadChannels.root","AliEMCALBadChannels"); 
    }
  
  TObjArray *arrayBC=(TObjArray*)contBC->GetObject(runBC);
  if(!arrayBC){
    AliError(Form("No external hot channel set for run number: %d", runBC));
    return 2; 
  }
  
  TObjArray *arrayBCpass=(TObjArray*)arrayBC->FindObject(fFilepass); // Here, it looks for a specific pass
  if(!arrayBCpass){
    AliError(Form("No external hot channel set for: %d -%s", runBC,fFilepass.Data()));
    return 2; 
  }

  if(fDebugLevel>0) arrayBCpass->Print();

  Int_t sms = fEMCALGeo->GetEMCGeometry()->GetNumberOfSuperModules();
  for (Int_t i=0; i<sms; ++i) {
    TH2I *h = fEMCALRecoUtils->GetEMCALChannelStatusMap(i);
    if (h)
      delete h;
    h=(TH2I*)arrayBCpass->FindObject(Form("EMCALBadChannelMap_Mod%d",i));

    if (!h) {
      AliError(Form("Can not get EMCALBadChannelMap_Mod%d",i));
      continue;
    }
    h->SetDirectory(0);
    fEMCALRecoUtils->SetEMCALChannelStatusMap(i,h);
  }
  
  return 1;  
}


//_____________________________________________________
Int_t AliEMCALTenderSupply::InitRecalib()
{
  // Initialising recalibration factors.
  
  AliESDEvent *event = fTender->GetEvent();
  if (!event) 
    return 0;
  
  if (fDebugLevel>0) 
    AliInfo("Initialising recalibration factors");
  
  fEMCALRecoUtils->SwitchOnRecalibration();
  Int_t runRC = event->GetRunNumber();
      
  AliOADBContainer *contRF=new AliOADBContainer("");
  if(fBasePath!="") { //if fBasePath specified in the ->SetBasePath()
    if (fDebugLevel>0)	AliInfo(Form("Loading Recalib OADB from given path %s",fBasePath.Data()));
    TFile *fRecalib= new TFile(Form("%s/EMCALRecalib.root",fBasePath.Data()),"read");
    if (fRecalib->IsZombie()) {
      AliFatal(Form("EMCALRecalib.root not found in %s, aborting",fBasePath.Data()));
      return 0;
    }
    if(fRecalib) delete fRecalib;
    contRF->InitFromFile(Form("%s/EMCALRecalib.root",fBasePath.Data()),"AliEMCALRecalib");
  }
    else{ // Else choose the one in the $ALICE_ROOT directory
      if (fDebugLevel>0)	AliInfo("Loading Recalib OADB from $ALICE_ROOT/OADB/EMCAL");
      TFile *fRecalib= new TFile("$ALICE_ROOT/OADB/EMCAL/EMCALRecalib.root","read");
      if (fRecalib->IsZombie()) {
	AliFatal("$ALICE_ROOT/OADB/EMCAL/EMCALRecalib.root was not found, aborting");
	return 0;
      }
      if(fRecalib) delete fRecalib;
      contRF->InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALRecalib.root","AliEMCALRecalib");     
    }

  TObjArray *recal=(TObjArray*)contRF->GetObject(runRC); //GetObject(int runnumber)
  if(!recal){
    AliError(Form("No Objects for run: %d",runRC));
    return 2;
  } 

  TObjArray *recalpass=(TObjArray*)recal->FindObject(fFilepass);
  if(!recalpass){
    AliError(Form("No Objects for run: %d - %s",runRC,fFilepass.Data()));
    return 2;
  }

  TObjArray *recalib=(TObjArray*)recalpass->FindObject("Recalib");
  if(!recalib){
    AliError(Form("No Recalib histos found for  %d - %s",runRC,fFilepass.Data())); 
    return 2;
  }

  if(fDebugLevel>0) recalib->Print();

  Int_t sms = fEMCALGeo->GetEMCGeometry()->GetNumberOfSuperModules();
  for (Int_t i=0; i<sms; ++i) {
    TH2F *h = fEMCALRecoUtils->GetEMCALChannelRecalibrationFactors(i);
    if (h)
      delete h;
    h = (TH2F*)recalib->FindObject(Form("EMCALRecalFactors_SM%d",i));
    if (!h) {
      AliError(Form("Could not load EMCALRecalFactors_SM%d",i));
      continue;
    }
    h->SetDirectory(0);
    fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(i,h);
  }
  return 1;
}

//_____________________________________________________
void AliEMCALTenderSupply::RecalibrateCells()
{
  // Recalibrate cells.
  
  AliESDCaloCells *cells = fTender->GetEvent()->GetEMCALCells();
  
  Int_t nEMCcell = cells->GetNumberOfCells();
  for(Int_t icell=0; icell<nEMCcell; ++icell) {
    Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
    fEMCALGeo->GetCellIndex(cells->GetCellNumber(icell),imod,iTower,iIphi,iIeta);
    fEMCALGeo->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);	
    Double_t calibFactor = fEMCALRecoUtils->GetEMCALChannelRecalibrationFactor(imod,ieta,iphi);
    cells->SetCell(icell,cells->GetCellNumber(icell),cells->GetAmplitude(icell)*calibFactor,cells->GetTime(icell));	
  }	
}


//_____________________________________________________
void AliEMCALTenderSupply::UpdateCells()
{
  //Remove bad cells from the cell list 
  //This is required for later reclusterization

  AliESDCaloCells *cells = fTender->GetEvent()->GetEMCALCells();
  
  Int_t nEMCcell = cells->GetNumberOfCells();
  for(Int_t icell=0; icell<nEMCcell; ++icell) 
    {
      Int_t cellId = cells->GetCellNumber(icell);
      Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
      fEMCALGeo->GetCellIndex(cellId,imod,iTower,iIphi,iIeta);
      fEMCALGeo->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);
      if(fEMCALRecoUtils->GetEMCALChannelStatus(imod,ieta,iphi))
	{
	  cells->SetCell(icell,cellId,(-1)*cells->GetAmplitude(icell),cells->GetCellTime(cellId));
	}
    }
}


//_____________________________________________________
void AliEMCALTenderSupply::InitRecParam()
{
  if (fDebugLevel>0) 
    AliInfo("Initialize the recParam");
  fRecParam = new AliEMCALRecParam;
}


//_____________________________________________________
Bool_t AliEMCALTenderSupply::InitClusterization()
{
  // Initialing clusterization/unfolding algorithm and set all the needed parameters.
  
  AliESDEvent *event=fTender->GetEvent();
  if (!event) 
    return kFALSE;
  
  if (fDebugLevel>0) 
    AliInfo(Form("Initialising reclustering parameters: Clusterizer-%d",fRecParam->GetClusterizerFlag()));
  
  //---setup clusterizer
  delete fClusterizer;
  if     (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv1)
    fClusterizer = new AliEMCALClusterizerv1 (fEMCALGeo);
  else if(fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv2) 
    fClusterizer = new AliEMCALClusterizerv2(fEMCALGeo);
  else if(fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerNxN) {
    AliEMCALClusterizerNxN *clusterizer = new AliEMCALClusterizerNxN(fEMCALGeo);
    clusterizer->SetNRowDiff(fRecParam->GetNRowDiff());
    clusterizer->SetNColDiff(fRecParam->GetNColDiff());
    fClusterizer = clusterizer;
  } else {
    AliFatal(Form("Clusterizer < %d > not available", fRecParam->GetClusterizerFlag()));
    return kFALSE;
  }
  
  // Set the clustering parameters
  fClusterizer->SetECAClusteringThreshold( fRecParam->GetClusteringThreshold() );
  fClusterizer->SetECALogWeight          ( fRecParam->GetW0()                  );
  fClusterizer->SetMinECut               ( fRecParam->GetMinECut()             );    
  fClusterizer->SetUnfolding             ( fRecParam->GetUnfold()              );
  fClusterizer->SetECALocalMaxCut        ( fRecParam->GetLocMaxCut()           );
  fClusterizer->SetTimeCut               ( fRecParam->GetTimeCut()             );
  fClusterizer->SetTimeMin               ( fRecParam->GetTimeMin()             );
  fClusterizer->SetTimeMax               ( fRecParam->GetTimeMax()             );
  fClusterizer->SetInputCalibrated       ( kTRUE                               );
  fClusterizer->SetJustClusters          ( kTRUE                               );  
  
  // In case of unfolding after clusterization is requested, set the corresponding parameters
  if (fRecParam->GetUnfold()) {
    for (Int_t i = 0; i < 8; ++i) {
      fClusterizer->SetSSPars(i, fRecParam->GetSSPars(i));
    }
    for (Int_t i = 0; i < 3; ++i) {
      fClusterizer->SetPar5  (i, fRecParam->GetPar5(i));
      fClusterizer->SetPar6  (i, fRecParam->GetPar6(i));
    }
    fClusterizer->InitClusterUnfolding();
  }
  
  fClusterizer->SetDigitsArr(fDigitsArr);
  fClusterizer->SetOutput(0);
  fClusterArr = const_cast<TObjArray *>(fClusterizer->GetRecPoints());
  return kTRUE;
}

//_____________________________________________________
void AliEMCALTenderSupply::FillDigitsArray()
{
  // Fill digits from cells to a TClonesArray.
  
  AliESDEvent *event = fTender->GetEvent();
  if (!event)
    return;
  
  fDigitsArr->Clear("C");
  AliESDCaloCells *cells = event->GetEMCALCells();
  Int_t ncells = cells->GetNumberOfCells();
  for (Int_t icell = 0, idigit = 0; icell < ncells; ++icell) {
    Double_t cellAmplitude=0, cellTime=0;
    Short_t cellNumber=0;
    if (cells->GetCell(icell, cellNumber, cellAmplitude, cellTime) != kTRUE)
      break;
    AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->New(idigit));
    digit->SetId(cellNumber);
    digit->SetTime(cellTime);
    digit->SetTimeR(cellTime);
    digit->SetIndexInList(idigit);
    digit->SetType(AliEMCALDigit::kHG);
    digit->SetAmplitude(cellAmplitude);
    idigit++;
  }
}

//_____________________________________________________
void AliEMCALTenderSupply::Clusterize()
{
  // Clusterize.
  
  fClusterizer->Digits2Clusters("");
}

//_____________________________________________________
void AliEMCALTenderSupply::UpdateClusters()
{
  // Update ESD cluster list.
  
  AliESDEvent *event = fTender->GetEvent();
  if (!event)
    return;
  
  TClonesArray *clus = dynamic_cast<TClonesArray*>(event->FindListObject("caloClusters"));
  if (!clus) 
    clus = dynamic_cast<TClonesArray*>(event->FindListObject("CaloClusters"));
  if (!clus) {
    AliError(" Null pointer to calo clusters array, returning");
    return;
  }
  
  Int_t nents = clus->GetEntriesFast();
  for (Int_t i=0; i < nents; ++i) {
    AliESDCaloCluster *c = static_cast<AliESDCaloCluster*>(clus->At(i));
    if (!c)
      continue;
    if (c->IsEMCAL())
      delete clus->RemoveAt(i);
  }
  clus->Compress();
  RecPoints2Clusters(clus);
}

//_____________________________________________________
void AliEMCALTenderSupply::RecPoints2Clusters(TClonesArray *clus)
{
  // Convert AliEMCALRecoPoints to AliESDCaloClusters.
  // Cluster energy, global position, cells and their amplitude fractions are restored.
  
  AliESDEvent *event = fTender->GetEvent();
  if (!event)
    return;
  
  Int_t ncls = fClusterArr->GetEntriesFast();
  for(Int_t i=0, nout=clus->GetEntriesFast(); i < ncls; ++i) {
    AliEMCALRecPoint *recpoint = static_cast<AliEMCALRecPoint*>(fClusterArr->At(i));
    
    Int_t ncells_true = 0;
    const Int_t ncells = recpoint->GetMultiplicity();
    UShort_t   absIds[ncells];  
    Double32_t ratios[ncells];
    Int_t *dlist = recpoint->GetDigitsList();
    Float_t *elist = recpoint->GetEnergiesList();
    for (Int_t c = 0; c < ncells; ++c) {
      AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->At(dlist[c]));
      absIds[ncells_true] = digit->GetId();
      ratios[ncells_true] = elist[c]/digit->GetAmplitude();
      if (ratios[ncells_true] < 0.001) 
        continue;
      ++ncells_true;
    }
    
    if (ncells_true < 1) {
      AliWarning("Skipping cluster with no cells");
      continue;
    }
    
    // calculate new cluster position
    TVector3 gpos;
    recpoint->GetGlobalPosition(gpos);
    Float_t g[3];
    gpos.GetXYZ(g);
    
    AliESDCaloCluster *c = static_cast<AliESDCaloCluster*>(clus->New(nout++));
    c->SetType(AliVCluster::kEMCALClusterv1);
    c->SetE(recpoint->GetEnergy());
    c->SetPosition(g);
    c->SetNCells(ncells_true);
    c->SetDispersion(recpoint->GetDispersion());
    c->SetEmcCpvDistance(-1);            //not yet implemented
    c->SetChi2(-1);                      //not yet implemented
    c->SetTOF(recpoint->GetTime()) ;     //time-of-flight
    c->SetNExMax(recpoint->GetNExMax()); //number of local maxima
    Float_t elipAxis[2];
    recpoint->GetElipsAxis(elipAxis);
    c->SetM02(elipAxis[0]*elipAxis[0]) ;
    c->SetM20(elipAxis[1]*elipAxis[1]) ;
    AliESDCaloCluster *cesd = static_cast<AliESDCaloCluster*>(c);
    cesd->SetCellsAbsId(absIds);
    cesd->SetCellsAmplitudeFraction(ratios);
  }
}

//_____________________________________________________
void AliEMCALTenderSupply::GetPass()
{
  // Get passx from filename.
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  fInputTree = mgr->GetTree();
  
  if (!fInputTree) {
    AliError("Pointer to tree = 0, returning");
    return;
  }
  
  fInputFile = fInputTree->GetCurrentFile();
  if (!fInputFile) {
    AliError("Null pointer input file, returning");
    return;
  }
  
  TString fname(fInputFile->GetName());
  if     (fname.Contains("pass1")) fFilepass = TString("pass1");
  else if(fname.Contains("pass2")) fFilepass = TString("pass2");
  else if(fname.Contains("pass3")) fFilepass = TString("pass3");
  else {
    AliError(Form("Pass number string not found: %s", fname.Data()));
    return;            
  }
}

