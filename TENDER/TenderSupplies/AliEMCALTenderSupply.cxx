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

#include <TClonesArray.h>
#include <TFile.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TInterpreter.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TTree.h>
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliDataFile.h"
#include "AliEMCALAfterBurnerUF.h"
#include "AliEMCALClusterizer.h"
#include "AliEMCALClusterizerNxN.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALClusterizerv2.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecParam.h"
#include "AliEMCALRecParam.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliMagF.h"
#include "AliOADBContainer.h"
#include "AliTender.h"
#include "AliEMCALTenderSupply.h"

ClassImp(AliEMCALTenderSupply)

AliEMCALTenderSupply::AliEMCALTenderSupply() :
  AliTenderSupply()
  ,fTask(0)
  ,fRun(0)
  ,fEMCALGeo(0x0)
  ,fEMCALGeoName("")
  ,fEMCALRecoUtils(0)
  ,fConfigName("")
  ,fDebugLevel(0)
  ,fNonLinearFunc(-1) 
  ,fNonLinearThreshold(-1)
  ,fReCalibCluster(kFALSE)
  ,fUpdateCell(kFALSE)  
  ,fCalibrateEnergy(kFALSE)
  ,fCalibrateTime(kFALSE)
  ,fCalibrateTimeParamAvailable(kFALSE)
  ,fCalibrateTimeL1Phase(kFALSE)
  ,fDoNonLinearity(kFALSE)
  ,fBadCellRemove(kFALSE)
  ,fLoad1DBadChMap(kFALSE)
  ,fLoad1DRecalibFactors(kFALSE)
  ,fRejectExoticCells(kFALSE)
  ,fRejectExoticClusters(kFALSE)
  ,fClusterBadChannelCheck(kFALSE)
  ,fRecalClusPos(kFALSE)
  ,fFiducial(kFALSE) 
  ,fNCellsFromEMCALBorder(-1)
  ,fRecalDistToBadChannels(kFALSE)
  ,fRecalShowerShape(kFALSE)
  ,fInputTree(0)
  ,fInputFile(0)
  ,fGetPassFromFileName(kTRUE)
  ,fFilepass(0) 
  ,fMass(-1)
  ,fStep(-1)
  ,fCutEtaPhiSum(kTRUE)
  ,fCutEtaPhiSeparate(kFALSE)
  ,fRcut(-1)
  ,fEtacut(-1)
  ,fPhicut(-1)
  ,fBasePath("")
  ,fCustomBC("")
  ,fCustomTempCalibSM("")
  ,fCustomTempCalibParams("")
  ,fReClusterize(kFALSE)
  ,fClusterizer(0)
  ,fGeomMatrixSet(kFALSE)
  ,fLoadGeomMatrices(kFALSE)
  ,fRecParam(0x0)
  ,fDoTrackMatch(kFALSE)
  ,fDoUpdateOnly(kFALSE)
  ,fUnfolder(0)
  ,fDigitsArr(0)
  ,fClusterArr(0)
  ,fMisalignSurvey(kdefault)  
  ,fExoticCellFraction(-1)
  ,fExoticCellDiffTime(-1)
  ,fExoticCellMinAmplitude(-1)
  ,fDoMergedBCs(kFALSE)
  ,fSetCellMCLabelFromCluster(0)
  ,fSetCellMCLabelFromEdepFrac(0)
  ,fTempClusterArr(0)
  ,fRemapMCLabelForAODs(0)
  ,fUseAutomaticRecalib(1)
  ,fUseAutomaticRunDepRecalib(1)
  ,fUseNewRunDepTempCalib(0)
  ,fUseAutomaticTimeCalib(1)
  ,fUseAutomaticRecParam(1)
{
  // Default constructor.

  for(Int_t i = 0; i < AliEMCALGeoParams::fgkEMCALModules; i++) fEMCALMatrix[i] = 0 ;
  for(Int_t j = 0; j < fgkTotalCellNumber;                 j++) 
  { fOrgClusterCellId[j] =-1; fCellLabels[j] =-1 ; }
}

//_____________________________________________________
AliEMCALTenderSupply::AliEMCALTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender)
  ,fTask(0)
  ,fRun(0)
  ,fEMCALGeo(0x0)
  ,fEMCALGeoName("")
  ,fEMCALRecoUtils(0)
  ,fConfigName("")
  ,fDebugLevel(0)
  ,fNonLinearFunc(-1) 
  ,fNonLinearThreshold(-1)        
  ,fReCalibCluster(kFALSE)  
  ,fUpdateCell(kFALSE)  
  ,fCalibrateEnergy(kFALSE)
  ,fCalibrateTime(kFALSE)
  ,fCalibrateTimeParamAvailable(kFALSE)
  ,fCalibrateTimeL1Phase(kFALSE)
  ,fDoNonLinearity(kFALSE)
  ,fBadCellRemove(kFALSE)
  ,fLoad1DBadChMap(kFALSE)
  ,fLoad1DRecalibFactors(kFALSE)
  ,fRejectExoticCells(kFALSE)
  ,fRejectExoticClusters(kFALSE)
  ,fClusterBadChannelCheck(kFALSE)
  ,fRecalClusPos(kFALSE)
  ,fFiducial(kFALSE) 
  ,fNCellsFromEMCALBorder(-1)  
  ,fRecalDistToBadChannels(kFALSE)  
  ,fRecalShowerShape(kFALSE)
  ,fInputTree(0)  
  ,fInputFile(0)
  ,fGetPassFromFileName(kTRUE)
  ,fFilepass("") 
  ,fMass(-1)
  ,fStep(-1)
  ,fCutEtaPhiSum(kTRUE)
  ,fCutEtaPhiSeparate(kFALSE)
  ,fRcut(-1)  
  ,fEtacut(-1)  
  ,fPhicut(-1)  
  ,fBasePath("")
  ,fCustomBC("")
  ,fCustomTempCalibSM("")
  ,fCustomTempCalibParams("")
  ,fReClusterize(kFALSE)
  ,fClusterizer(0)
  ,fGeomMatrixSet(kFALSE)
  ,fLoadGeomMatrices(kFALSE)
  ,fRecParam(0x0)
  ,fDoTrackMatch(kFALSE)
  ,fDoUpdateOnly(kFALSE)
  ,fUnfolder(0)
  ,fDigitsArr(0)
  ,fClusterArr(0)
  ,fMisalignSurvey(kdefault)  
  ,fExoticCellFraction(-1)
  ,fExoticCellDiffTime(-1)
  ,fExoticCellMinAmplitude(-1)
  ,fDoMergedBCs(kFALSE)
  ,fSetCellMCLabelFromCluster(0)
  ,fSetCellMCLabelFromEdepFrac(0)
  ,fTempClusterArr(0)
  ,fRemapMCLabelForAODs(0)
  ,fUseAutomaticRecalib(1)
  ,fUseAutomaticRunDepRecalib(1)
  ,fUseNewRunDepTempCalib(0)
  ,fUseAutomaticTimeCalib(1)
  ,fUseAutomaticRecParam(1)
{
  // Named constructor
  
  for(Int_t i = 0; i < AliEMCALGeoParams::fgkEMCALModules; i++) fEMCALMatrix[i] = 0 ;
  for(Int_t j = 0; j < fgkTotalCellNumber;                 j++) 
  { fOrgClusterCellId[j] =-1; fCellLabels[j] =-1 ; }
}

//_____________________________________________________
AliEMCALTenderSupply::AliEMCALTenderSupply(const char *name, AliAnalysisTaskSE *task) :
  AliTenderSupply(name)
  ,fTask(task)
  ,fRun(0)
  ,fEMCALGeo(0x0)
  ,fEMCALGeoName("")
  ,fEMCALRecoUtils(0)
  ,fConfigName("")
  ,fDebugLevel(0)
  ,fNonLinearFunc(-1) 
  ,fNonLinearThreshold(-1)        
  ,fReCalibCluster(kFALSE)  
  ,fUpdateCell(kFALSE)  
  ,fCalibrateEnergy(kFALSE)
  ,fCalibrateTime(kFALSE)
  ,fCalibrateTimeParamAvailable(kFALSE)
  ,fCalibrateTimeL1Phase(kFALSE)
  ,fDoNonLinearity(kFALSE)
  ,fBadCellRemove(kFALSE)
  ,fRejectExoticCells(kFALSE)
  ,fRejectExoticClusters(kFALSE)
  ,fClusterBadChannelCheck(kFALSE)
  ,fRecalClusPos(kFALSE)
  ,fFiducial(kFALSE) 
  ,fNCellsFromEMCALBorder(-1)  
  ,fRecalDistToBadChannels(kFALSE)  
  ,fRecalShowerShape(kFALSE)
  ,fInputTree(0)  
  ,fInputFile(0)
  ,fGetPassFromFileName(kTRUE)
  ,fFilepass("") 
  ,fMass(-1)
  ,fStep(-1)
  ,fCutEtaPhiSum(kTRUE)
  ,fCutEtaPhiSeparate(kFALSE)
  ,fRcut(-1)  
  ,fEtacut(-1)  
  ,fPhicut(-1)  
  ,fBasePath("")
  ,fCustomBC("")
  ,fCustomTempCalibSM("")
  ,fCustomTempCalibParams("")
  ,fReClusterize(kFALSE)
  ,fClusterizer(0)
  ,fGeomMatrixSet(kFALSE)
  ,fLoadGeomMatrices(kFALSE)
  ,fRecParam(0x0)
  ,fDoTrackMatch(kFALSE)
  ,fDoUpdateOnly(kFALSE)
  ,fUnfolder(0)
  ,fDigitsArr(0)
  ,fClusterArr(0)
  ,fMisalignSurvey(kdefault)  
  ,fExoticCellFraction(-1)
  ,fExoticCellDiffTime(-1)
  ,fExoticCellMinAmplitude(-1)
  ,fDoMergedBCs(kFALSE)
  ,fSetCellMCLabelFromCluster(0)
  ,fSetCellMCLabelFromEdepFrac(0)
  ,fTempClusterArr(0)
  ,fRemapMCLabelForAODs(0)
  ,fUseAutomaticRecalib(1)
  ,fUseAutomaticRunDepRecalib(1)
  ,fUseNewRunDepTempCalib(0)
  ,fUseAutomaticTimeCalib(1)
  ,fUseAutomaticRecParam(1)
{
  // Named constructor.
  
  for(Int_t i = 0; i < AliEMCALGeoParams::fgkEMCALModules; i++) fEMCALMatrix[i] = 0 ;
  for(Int_t j = 0; j < fgkTotalCellNumber;                 j++) 
  { fOrgClusterCellId[j] =-1; fCellLabels[j] =-1 ; }
}

//_____________________________________________________
AliEMCALTenderSupply::~AliEMCALTenderSupply()
{
  //Destructor

  if (!AliAnalysisManager::GetAnalysisManager())  return;  

  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) 
  {
    delete fEMCALRecoUtils;
    delete fRecParam;
    delete fUnfolder;
    
    if (!fClusterizer) 
    {
      if (fDigitsArr) 
      { 
     	fDigitsArr->Clear("C");
      	delete fDigitsArr; 
      }
    } 
    else 
    {
      delete fClusterizer;
      fDigitsArr = 0;
    }
  }
}

//_____________________________________________________
void AliEMCALTenderSupply::SetDefaults()
{
  // Set default settings.

  SwitchOnReclustering();
  SwitchOnTrackMatch();
}

//_____________________________________________________
Bool_t AliEMCALTenderSupply::RunChanged() const
{
  // Get run number.

  return (fTender && fTender->RunChanged()) || (fTask && fRun != fTask->InputEvent()->GetRunNumber()); 
}

//_____________________________________________________
void AliEMCALTenderSupply::Init()
{
  // Initialise EMCAL tender.
  
  if (fDebugLevel>0) {
    AliWarning("Init EMCAL Tender supply");
    AliWarning("======================================================================================");
    AliWarning("--------------------------------------------------------------------------------------");
    AliWarning("============  The EMCal Tender is no longer supported for development!  ==============");
    AliWarning("=========  It is recommended to use the new EMCal Correction Framework!  =============");
    AliWarning("===  See http://alidoc.cern.ch/AliPhysics/master/_r_e_a_d_m_eemc_corrections.html  ===");
    AliWarning("--------------------------------------------------------------------------------------");
    AliWarning("======================================================================================");
  }
  
  if (fConfigName.Length()>0 && gROOT->LoadMacro(fConfigName) >=0) {
    AliDebug(1, Form("Loading settings from macro %s", fConfigName.Data()));
    AliEMCALTenderSupply *tender = (AliEMCALTenderSupply*)gInterpreter->ProcessLine("ConfigEMCALTenderSupply()");
    fDebugLevel             = tender->fDebugLevel;
    fEMCALGeoName           = tender->fEMCALGeoName; 
    fEMCALRecoUtils         = tender->fEMCALRecoUtils; 
    fConfigName             = tender->fConfigName;
    fNonLinearFunc          = tender->fNonLinearFunc;
    fNonLinearThreshold     = tender->fNonLinearThreshold;
    fReCalibCluster         = tender->fReCalibCluster;
    fUpdateCell             = tender->fUpdateCell;
    fRecalClusPos           = tender->fRecalClusPos;
    fCalibrateEnergy        = tender->fCalibrateEnergy;
    fCalibrateTime          = tender->fCalibrateTime;
    fCalibrateTimeParamAvailable = tender->fCalibrateTimeParamAvailable;
    fCalibrateTimeL1Phase   = tender->fCalibrateTimeL1Phase;
    fFiducial               = tender->fFiducial;
    fNCellsFromEMCALBorder  = tender->fNCellsFromEMCALBorder;
    fRecalDistToBadChannels = tender->fRecalDistToBadChannels;    
    fRecalShowerShape       = tender->fRecalShowerShape;
    fClusterBadChannelCheck = tender->fClusterBadChannelCheck;
    fBadCellRemove          = tender->fBadCellRemove;
    fRejectExoticCells      = tender->fRejectExoticCells;
    fRejectExoticClusters   = tender->fRejectExoticClusters;
    fMass                   = tender->fMass;
    fStep                   = tender->fStep;
    fCutEtaPhiSum           = tender->fCutEtaPhiSum;
    fCutEtaPhiSeparate      = tender->fCutEtaPhiSeparate;
    fRcut                   = tender->fRcut;
    fEtacut                 = tender->fEtacut;
    fPhicut                 = tender->fPhicut;
    fReClusterize           = tender->fReClusterize;
    fLoadGeomMatrices       = tender->fLoadGeomMatrices;
    fRecParam               = tender->fRecParam;
    fDoNonLinearity         = tender->fDoNonLinearity;
    fDoTrackMatch           = tender->fDoTrackMatch;
    fDoUpdateOnly           = tender->fDoUpdateOnly;
    fMisalignSurvey         = tender->fMisalignSurvey;
    fExoticCellFraction     = tender->fExoticCellFraction;
    fExoticCellDiffTime     = tender->fExoticCellDiffTime;
    fExoticCellMinAmplitude = tender->fExoticCellMinAmplitude;

    for(Int_t i = 0; i < AliEMCALGeoParams::fgkEMCALModules; i++) 
      fEMCALMatrix[i] = tender->fEMCALMatrix[i] ;
  }
  
  if (fDebugLevel>0){
    AliInfo("Emcal Tender settings:");
    AliInfo("------------ Switches --------------------------"); 
    AliInfo(Form("BadCellRemove : %d", fBadCellRemove)); 
    AliInfo(Form("ExoticCellRemove : %d", fRejectExoticCells)); 
    AliInfo(Form("CalibrateEnergy : %d", fCalibrateEnergy)); 
    AliInfo(Form("CalibrateTime : %d", fCalibrateTime)); 
    AliInfo(Form("CalibrateTimeL1Phase : %d", fCalibrateTimeL1Phase));
    AliInfo(Form("UpdateCell : %d", fUpdateCell)); 
    AliInfo(Form("DoUpdateOnly : %d", fDoUpdateOnly)); 
    AliInfo(Form("Reclustering : %d", fReClusterize)); 
    AliInfo(Form("ClusterBadChannelCheck : %d", fClusterBadChannelCheck)); 
    AliInfo(Form("ClusterExoticChannelCheck : %d", fRejectExoticClusters)); 
    AliInfo(Form("CellFiducialRegion : %d", fFiducial)); 
    AliInfo(Form("ReCalibrateCluster : %d", fReCalibCluster)); 
    AliInfo(Form("RecalculateClusPos : %d", fRecalClusPos)); 
    AliInfo(Form("RecalShowerShape : %d", fRecalShowerShape)); 
    AliInfo(Form("NonLinearityCorrection : %d", fDoNonLinearity)); 
    AliInfo(Form("RecalDistBadChannel : %d", fRecalDistToBadChannels)); 
    AliInfo(Form("TrackMatch : %d", fDoTrackMatch)); 
    AliInfo("------------ Variables -------------------------"); 
    AliInfo(Form("DebugLevel : %d", fDebugLevel)); 
    AliInfo(Form("BasePath : %s", fBasePath.Data())); 
    AliInfo(Form("CustomBCFile : %s", fCustomBC.Data()));
    AliInfo(Form("CustomSMTempFile : %s", fCustomTempCalibSM.Data()));
    AliInfo(Form("CustomTempParamFile : %s", fCustomTempCalibParams.Data()));
    AliInfo(Form("ConfigFileName : %s", fConfigName.Data())); 
    AliInfo(Form("EMCALGeometryName : %s", fEMCALGeoName.Data())); 
    AliInfo(Form("NonLinearityFunction : %d", fNonLinearFunc)); 
    AliInfo(Form("NonLinearityThreshold : %d", fNonLinearThreshold)); 
    AliInfo(Form("MisalignmentMatrixSurvey : %d", fMisalignSurvey)); 
    AliInfo(Form("NumberOfCellsFromEMCALBorder : %d", fNCellsFromEMCALBorder)); 
    AliInfo(Form("RCut : %f", fRcut)); 
    AliInfo(Form("Mass : %f", fMass)); 
    AliInfo(Form("Step : %f", fStep)); 
    AliInfo(Form("EtaCut : %f", fEtacut)); 
    AliInfo(Form("PhiCut : %f", fPhicut)); 
    AliInfo(Form("ExoticCellFraction : %f", fExoticCellFraction)); 
    AliInfo(Form("ExoticCellDiffTime : %f", fExoticCellDiffTime)); 
    AliInfo(Form("ExoticCellMinAmplitude : %f", fExoticCellMinAmplitude)); 
    AliInfo("============================================================="); 
  }

  // init reco utils
  
  if (!fEMCALRecoUtils)
    fEMCALRecoUtils  = new AliEMCALRecoUtils;

  // init geometry if requested
  if (fEMCALGeoName.Length()>0) 
    fEMCALGeo = AliEMCALGeometry::GetInstance(fEMCALGeoName) ;

  // Use one histogram for all BCs
  if (fDoMergedBCs)
    fEMCALRecoUtils->SetUseOneHistForAllBCs(fDoMergedBCs);

  // digits array
  fDigitsArr       = new TClonesArray("AliEMCALDigit",1000);

  // initialising non-linearity parameters
  if (fNonLinearThreshold != -1)
    fEMCALRecoUtils->SetNonLinearityThreshold(fNonLinearThreshold);
  if (fNonLinearFunc != -1)
    fEMCALRecoUtils->SetNonLinearityFunction(fNonLinearFunc);

  // missalignment function
  fEMCALRecoUtils->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);

  // fiducial cut
  // do not do the eta0 fiducial cut
  if (fNCellsFromEMCALBorder != -1)
    fEMCALRecoUtils->SetNumberOfCellsFromEMCALBorder(fNCellsFromEMCALBorder);
  fEMCALRecoUtils->SwitchOnNoFiducialBorderInEMCALEta0();
    
  // exotic cell rejection
  if (fExoticCellFraction != -1)
    fEMCALRecoUtils->SetExoticCellFractionCut(fExoticCellFraction);
  if (fExoticCellDiffTime != -1)
    fEMCALRecoUtils->SetExoticCellDiffTimeCut(fExoticCellDiffTime);
  if (fExoticCellMinAmplitude != -1)
    fEMCALRecoUtils->SetExoticCellMinAmplitudeCut(fExoticCellMinAmplitude);

  // setting track matching parameters ... mass, step size and residual cut
  if (fMass != -1)
    fEMCALRecoUtils->SetMass(fMass);
  if (fStep != -1)
    fEMCALRecoUtils->SetStep(fStep);
  
  // spatial cut based on separate eta/phi or common processing
  if (fCutEtaPhiSum) { 
    fEMCALRecoUtils->SwitchOnCutEtaPhiSum(); 
    if (fRcut != -1)
      fEMCALRecoUtils->SetCutR(fRcut);
  } else if (fCutEtaPhiSeparate) {
    fEMCALRecoUtils->SwitchOnCutEtaPhiSeparate();
    if (fEtacut != -1)
      fEMCALRecoUtils->SetCutEta(fEtacut);
    if (fPhicut != -1)
      fEMCALRecoUtils->SetCutPhi(fPhicut);
  }
}

//_____________________________________________________
AliVEvent* AliEMCALTenderSupply::GetEvent()
{
  // Return the event pointer.
  
  if (fTender) {
    return fTender->GetEvent();
  } else if (fTask) {
    return fTask->InputEvent();
  }
  
  return 0;
}

//_____________________________________________________
AliMCEvent* AliEMCALTenderSupply::GetMCEvent()
{
  // Return the event pointer.
  
  if (fTender) {
    return fTender->MCEvent();
  } else if (fTask) {
    return fTask->MCEvent();
  }
  
  return 0;
}

//_____________________________________________________
void AliEMCALTenderSupply::ProcessEvent()
{
  // Event loop.
  
  AliVEvent *event = GetEvent();

  if (!event) {
    AliError("Event ptr = 0, returning");
    return;
  }
  
  // Initialising parameters once per run number
  
  if (RunChanged()) 
  { 
    fRun = event->GetRunNumber();
    AliWarning(Form("Run changed, initializing parameters for %d", fRun));
    if (dynamic_cast<AliAODEvent*>(event)) {
      AliWarning("============================================================="); 
      AliWarning("===  Running on AOD is not equivalent to running on ESD!  ===");
      AliWarning("============================================================="); 
    }

    // init geometry if not already done
    if (fEMCALGeoName.Length()==0) 
    {
      fEMCALGeo = AliEMCALGeometry::GetInstanceFromRunNumber(fRun);
      
      if (!fEMCALGeo) 
      {
        AliFatal(Form("Can not create geometry: %s", fEMCALGeoName.Data()));
        return;
      }
    } 

    // get pass
    if (fGetPassFromFileName)
      GetPass();

    // define what recalib parameters are needed for various switches
    // this is based on implementation in AliEMCALRecoUtils
    Bool_t needRecoParam   = fReClusterize;
    Bool_t needBadChannels = fBadCellRemove   | fClusterBadChannelCheck | fRecalDistToBadChannels | fReClusterize;
    Bool_t needRecalib     = fCalibrateEnergy | fReClusterize;
    Bool_t needTimecalib   = fCalibrateTime   | fReClusterize;
    Bool_t needMisalign    = fRecalClusPos    | fReClusterize;
    Bool_t needClusterizer = fReClusterize;
    
    if ( fRun > 209121 && fCalibrateTime ) fCalibrateTimeL1Phase = kTRUE;
    Bool_t needTimecalibL1Phase = (fCalibrateTime | fReClusterize) & fCalibrateTimeL1Phase;

    // init bad channels
    if (needBadChannels) {
      Int_t fInitBC = InitBadChannels();
      if (fInitBC==0)
        AliError("InitBadChannels returned false, returning");
      if (fInitBC==1)
        AliWarning("InitBadChannels OK");
      if (fInitBC>1)
        AliWarning(Form("No external hot channel set: %d - %s", event->GetRunNumber(), fFilepass.Data()));
    }

    // init recalibration factors
    if (needRecalib) 
    {
      if(fUseAutomaticRecalib)
      {
        Int_t fInitRecalib = InitRecalib();
        if (fInitRecalib==0)
          AliError("InitRecalib returned false, returning");
        if (fInitRecalib==1)
          AliWarning("InitRecalib OK");
        if (fInitRecalib>1)
          AliWarning(Form("No recalibration available: %d - %s", event->GetRunNumber(), fFilepass.Data()));
      }
      
      if(fUseAutomaticRunDepRecalib)
      {
        Int_t fInitRunDepRecalib = InitRunDepRecalib();
        if (fInitRunDepRecalib==0)
          AliError("InitrunDepRecalib returned false, returning");
        if (fInitRunDepRecalib==1)
          AliWarning("InitRecalib OK");
        if (fInitRunDepRecalib>1)
          AliWarning(Form("No Temperature recalibration available: %d - %s", event->GetRunNumber(), fFilepass.Data()));
      }
    }
    
    // init time calibration
    if (needTimecalib && fUseAutomaticTimeCalib) {
      Int_t initTC = InitTimeCalibration();
      if (!initTC) 
        AliError("InitTimeCalibration returned false, returning");
      if (initTC==1) {
        fCalibrateTimeParamAvailable = kTRUE;
        AliWarning("InitTimeCalib OK");
      }
      if (initTC > 1)
        AliWarning(Form("No external time calibration available: %d - %s", event->GetRunNumber(), fFilepass.Data()));
    }

    // init time calibration with L1 phase
    if (needTimecalibL1Phase && fUseAutomaticTimeCalib) {
      Int_t initTCL1Phase = InitTimeCalibrationL1Phase();
      if (!initTCL1Phase) 
	AliError("InitTimeCalibrationL1Phase returned false, returning");
      if (initTCL1Phase==1) {
	fCalibrateTimeParamAvailable = kTRUE;
	AliWarning("InitTimeCalibL1Phase OK");
      }
      if (initTCL1Phase > 1)
	AliWarning(Form("No external time calibration L1 phase available: %d - %s", event->GetRunNumber(), fFilepass.Data()));
    }

    // init misalignment matrix
    if (needMisalign) {
      if (!InitMisalignMatrix())
        AliError("InitMisalignmentMatrix returned false, returning");
      else
        AliWarning("InitMisalignMatrix OK");
    }
    
    // initiate reco params with some defaults
    // will not overwrite, if those have been provided by user
    if (needRecoParam && fUseAutomaticRecParam) {
      Int_t initRC = InitRecParam();
      
      if (initRC == 0)
        AliInfo("Defaults reco params loaded.");
      if (initRC > 1)
        AliWarning("User defined reco params.");
    }
    
    // init clusterizer
    if (needClusterizer) {
      if (!InitClusterization()) 
        AliError("InitClusterization returned false, returning");
      else
        AliWarning("InitClusterization OK");
    }
    
    if (fDebugLevel>1) 
      fEMCALRecoUtils->Print("");
  }
  
  // disable implied switches -------------------------------------------------
  // AliEMCALRecoUtils or clusterizer functions alredy include some
  // recalibration so based on those implied callibration te switches
  // here are disabled to avoid duplication
    
  // clusterizer does cluster energy recalibration, position recomputation
  // and shower shape
  if (fReClusterize) {
    fReCalibCluster   = kFALSE;
    fRecalClusPos     = kFALSE;
    fRecalShowerShape = kFALSE;
  }
  
  // CONFIGURE THE RECO UTILS -------------------------------------------------
  // configure the reco utils
  // this option does energy recalibration
  if (fCalibrateEnergy)
    fEMCALRecoUtils->SwitchOnRecalibration();
  else
    fEMCALRecoUtils->SwitchOffRecalibration();
  
  // allows time calibration
  if (fCalibrateTime)
    fEMCALRecoUtils->SwitchOnTimeRecalibration();
  else
    fEMCALRecoUtils->SwitchOffTimeRecalibration();

  // allows time calibration with L1 phase
  if (fCalibrateTimeL1Phase)
    fEMCALRecoUtils->SwitchOnL1PhaseInTimeRecalibration();
  else
    fEMCALRecoUtils->SwitchOffL1PhaseInTimeRecalibration();

  // allows to zero bad cells
  if (fBadCellRemove)
    fEMCALRecoUtils->SwitchOnBadChannelsRemoval();
  else
    fEMCALRecoUtils->SwitchOffBadChannelsRemoval();
  
  // distance to bad channel recalibration
  if (fRecalDistToBadChannels)
    fEMCALRecoUtils->SwitchOnDistToBadChannelRecalculation();
  else
    fEMCALRecoUtils->SwitchOffDistToBadChannelRecalculation();

  // exclude exotic cells
  if (fRejectExoticCells)
    fEMCALRecoUtils->SwitchOnRejectExoticCell();
  else
    fEMCALRecoUtils->SwitchOffRejectExoticCell();
  
  // exclude clusters with exotic cells
  if (fRejectExoticClusters)
    fEMCALRecoUtils->SwitchOnRejectExoticCluster();
  else
    fEMCALRecoUtils->SwitchOffRejectExoticCluster();

  // TODO: not implemented switches
  // SwitchOnClusterEnergySmearing
  // SwitchOnRunDepCorrection

  // START PROCESSING ---------------------------------------------------------
  // Test if cells present
  AliVCaloCells *cells= event->GetEMCALCells();
  if (cells->GetNumberOfCells()<=0) 
  {
    if (fDebugLevel>1) 
      AliWarning(Form("Number of EMCAL cells = %d, returning", cells->GetNumberOfCells()));
    return;
  }
  
  if (fDebugLevel>2)
    AliInfo(Form("Re-calibrate cluster %d\n",fReCalibCluster));

  // mark the cells not recalibrated in case of selected
  // time, energy recalibration or bad channel removal
  if (fCalibrateEnergy || fCalibrateTime || fBadCellRemove || fCalibrateTimeL1Phase)
    fEMCALRecoUtils->ResetCellsCalibrated();
  
 // CELL RECALIBRATION -------------------------------------------------------
  // cell objects will be updated
  // the cell calibrations are also processed locally any time those are needed
  // in case that the cell objects are not to be updated here for later use
  if (fUpdateCell) {
    // do the update
    // includes exotic cell check in the UpdateCells function - is not provided
    // by the reco utils
    UpdateCells();

    // switch off recalibrations so those are not done multiple times
    // this is just for safety, the recalibrated flag of cell object
    // should not allow for farther processing anyways
    fEMCALRecoUtils->SwitchOffRecalibration();
    fEMCALRecoUtils->SwitchOffTimeRecalibration();  
    fEMCALRecoUtils->SwitchOffL1PhaseInTimeRecalibration();

    if (fDoUpdateOnly)
      return;
  }

  // RECLUSTERIZATION ---------------------------------------------------------
  if (fReClusterize)
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

  // PROCESS SINGLE CLUSTER RECALIBRATION -------------------------------------
  // now go through clusters one by one and process separate correction
  // as those were defined or not
  Int_t nclusters = clusArr->GetEntriesFast();
  for (Int_t icluster=0; icluster < nclusters; ++icluster) 
  { 
    AliVCluster *clust = static_cast<AliVCluster*>(clusArr->At(icluster));
    if (!clust) 
      continue;
    if  (!clust->IsEMCAL()) 
      continue;

    // REMOVE CLUSTERS WITH BAD CELLS -----------------------------
    if (fClusterBadChannelCheck) {
      // careful, the the ClusterContainsBadChannel is dependent on
      // SwitchOnBadChannelsRemoval, switching it ON automatically
      // and returning to original value after processing
      Bool_t badRemoval = fEMCALRecoUtils->IsBadChannelsRemovalSwitchedOn();
      fEMCALRecoUtils->SwitchOnBadChannelsRemoval();
      
      Bool_t badResult = fEMCALRecoUtils->ClusterContainsBadChannel(fEMCALGeo, clust->GetCellsAbsId(), clust->GetNCells());

      // switch the bad channels removal back
      if (!badRemoval)
        fEMCALRecoUtils->SwitchOffBadChannelsRemoval();
      
      if (badResult)
      {
        delete clusArr->RemoveAt(icluster);
        continue; //TODO is it really needed to remove it? Or should we flag it?
      }
    }
    
    // REMOVE EXOTIC CLUSTERS -------------------------------------
    // does process local cell recalibration energy and time without replacing
    // the global cell values, in case of no cell recalib done yet
    if (fRejectExoticClusters)
    {
      // careful, the IsExoticCluster is dependent on
      // SwitchOnRejectExoticCell, switching it ON automatically
      // and returning to original value after processing
      Bool_t exRemoval = fEMCALRecoUtils->IsRejectExoticCell();
      fEMCALRecoUtils->SwitchOnRejectExoticCell();

      // get bunch crossing
      Int_t bunchCrossNo = event->GetBunchCrossNumber();

      Bool_t exResult = fEMCALRecoUtils->IsExoticCluster(clust, cells, bunchCrossNo);

      // switch the exotic channels removal back
      if (!exRemoval)
        fEMCALRecoUtils->SwitchOffRejectExoticCell();
      
      if (exResult) {
        delete clusArr->RemoveAt(icluster);
        continue; //TODO is it really needed to remove it? Or should we flag it?
      }
    }
    
    // FIDUCIAL CUT -----------------------------------------------
    if (fFiducial) {
      // depends on SetNumberOfCellsFromEMCALBorder
      // SwitchOnNoFiducialBorderInEMCALEta0
      if (!fEMCALRecoUtils->CheckCellFiducialRegion(fEMCALGeo, clust, cells)){
        delete clusArr->RemoveAt(icluster);
        continue; //TODO it would be nice to store the distance
      }
    }
    
    // CLUSTER ENERGY ---------------------------------------------
    // does process local cell recalibration energy and time without replacing
    // the global cell values, in case of no cell recalib done yet
    if (fReCalibCluster) {
      fEMCALRecoUtils->RecalibrateClusterEnergy(fEMCALGeo, clust, cells);
      if (clust->E() < 1e-9) {
        delete clusArr->RemoveAt(icluster);
        continue;
      }
    }
    
    // CLUSTER POSITION -------------------------------------------
    // does process local cell energy recalibration, if enabled and cells
    // not calibrated yet
    if (fRecalClusPos) 
      fEMCALRecoUtils->RecalculateClusterPosition(fEMCALGeo, cells, clust);
    
    // SHOWER SHAPE -----------------------------------------------
    if (fRecalShowerShape)
      fEMCALRecoUtils->RecalculateClusterShowerShapeParameters(fEMCALGeo, cells, clust);  

    // NONLINEARITY -----------------------------------------------
    if (fDoNonLinearity) {
      Float_t correctedEnergy = fEMCALRecoUtils->CorrectClusterEnergyLinearity(clust);
      clust->SetE(correctedEnergy);
    }

    // DISTANCE TO BAD CHANNELS -----------------------------------
    if (fRecalDistToBadChannels)
      fEMCALRecoUtils->RecalculateClusterDistanceToBadChannel(fEMCALGeo, cells, clust);  
  }

  clusArr->Compress();

  if (!fDoTrackMatch)
    return;

  // TRACK MATCHING -----------------------------------------------------------
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    event->InitMagneticField();
  }

  fEMCALRecoUtils->FindMatches(event,0x0,fEMCALGeo,GetMCEvent());
  fEMCALRecoUtils->SetClusterMatchedToTrack(event);
  fEMCALRecoUtils->SetTracksMatchedToCluster(event);
}

//_____________________________________________________
Bool_t AliEMCALTenderSupply::InitMisalignMatrix()
{
  // Initialising misalignment matrices

  AliVEvent *event = GetEvent();
  
  if (!event) 
    return kFALSE;
  
  if (fGeomMatrixSet) 
  {
    AliInfo("Misalignment matrix already set");  
    return kTRUE;
  }
  
  if (fDebugLevel>0) 
    AliInfo("Initialising misalignment matrix");  
  
  if (fLoadGeomMatrices) 
  {
    for(Int_t mod=0; mod < fEMCALGeo->GetNumberOfSuperModules(); ++mod)
    {
      if (fEMCALMatrix[mod])
      {
        if (fDebugLevel > 2) fEMCALMatrix[mod]->Print();
        
        fEMCALGeo->SetMisalMatrix(fEMCALMatrix[mod],mod);  
      }
      else if(gGeoManager)
      {
        AliWarning(Form("Set matrix for SM %d from gGeoManager\n",mod));
        fEMCALGeo->SetMisalMatrix(fEMCALGeo->GetMatrixForSuperModuleFromGeoManager(mod),mod) ;
      }
      else
      {
        AliError(Form("EMCal geometry matrix for SM %d is not available!",mod));
      }
    }
    fGeomMatrixSet = kTRUE;
    return kTRUE;
  }
  
  Int_t runGM = event->GetRunNumber();
  TObjArray *mobj = 0;

  if (fMisalignSurvey == kdefault)
  { //take default alignment corresponding to run no
    AliOADBContainer emcalgeoCont(Form("emcal"));
    emcalgeoCont.InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALlocal2master.root").data(),Form("AliEMCALgeo"));
    mobj=(TObjArray*)emcalgeoCont.GetObject(runGM,"EmcalMatrices");
  }
  
  if (fMisalignSurvey == kSurveybyS)
  { //take alignment at sector level
    if (runGM <= 140000) { //2010 data
      AliOADBContainer emcalgeoCont(Form("emcal2010"));
      emcalgeoCont.InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALlocal2master.root").data(),Form("AliEMCALgeo"));
      mobj=(TObjArray*)emcalgeoCont.GetObject(100,"survey10");
    }
    else if (runGM>140000)
    { // 2011 LHC11a pass1 data
      AliOADBContainer emcalgeoCont(Form("emcal2011"));
      emcalgeoCont.InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALlocal2master.root").data(),Form("AliEMCALgeo"));
      mobj=(TObjArray*)emcalgeoCont.GetObject(100,"survey11byS");
    }
  }

  if (fMisalignSurvey == kSurveybyM)
  { //take alignment at module level
    if (runGM <= 140000) { //2010 data
      AliOADBContainer emcalgeoCont(Form("emcal2010"));
      emcalgeoCont.InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALlocal2master.root").data(),Form("AliEMCALgeo"));
      mobj=(TObjArray*)emcalgeoCont.GetObject(100,"survey10");
    }
    else if (runGM>140000)
    { // 2011 LHC11a pass1 data
      AliOADBContainer emcalgeoCont(Form("emcal2011"));
      emcalgeoCont.InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALlocal2master.root").data(),Form("AliEMCALgeo"));
      mobj=(TObjArray*)emcalgeoCont.GetObject(100,"survey11byM");
    }
  }

  if (!mobj)
  {
    AliFatal("Geometry matrix array not found");
    return kFALSE;
  }
  
  for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
  {
    fEMCALMatrix[mod] = (TGeoHMatrix*) mobj->At(mod);
    
    if ( fEMCALMatrix[mod] )
    {
      if (fDebugLevel > 2) fEMCALMatrix[mod]->Print();

      fEMCALGeo->SetMisalMatrix(fEMCALMatrix[mod],mod);       
    }
    else
    {
      AliWarning(Form("EMCal geometry matrix for SM %d is not available!",mod));
    }
  }
  
  return kTRUE;
}

//_____________________________________________________
Int_t AliEMCALTenderSupply::InitBadChannels()
{
  // Initialising bad channel maps

  AliVEvent *event = GetEvent();

  if (!event) 
    return 0;
  
  if (fDebugLevel>0) 
    AliInfo("Initialising Bad channel map");
  
  // init default maps first
  if (!fEMCALRecoUtils->GetEMCALBadChannelStatusMapArray())
    fEMCALRecoUtils->InitEMCALBadChannelStatusMap1D() ;
  
  Int_t runBC = event->GetRunNumber();
  
  AliOADBContainer *contBC = new AliOADBContainer("");
  if (fBasePath!="")
  { //if fBasePath specified in the ->SetBasePath()
    if (fDebugLevel>0) AliInfo(Form("Loading Bad Channels OADB from given path %s",fBasePath.Data()));
    
    TFile *fbad=new TFile(Form("%s/EMCALBadChannels%s.root",fBasePath.Data(), fLoad1DBadChMap ? "_1D" : ""),"read");
    if (!fbad || fbad->IsZombie())
    {
      AliFatal(Form("EMCALBadChannels%s.root was not found in the path provided: %s", fLoad1DBadChMap ? "_1D" : "", fBasePath.Data()));
      return 0;
    }  
    
    if (fbad) delete fbad;
    
    contBC->InitFromFile(Form("%s/EMCALBadChannels%s.root",fBasePath.Data(), fLoad1DBadChMap ? "_1D" : ""),"AliEMCALBadChannels");
  }
  else if (fCustomBC!="")
  { //if fCustomBC specified in the ->SetCustomBC()
    if (fDebugLevel>0) AliInfo(Form("Loading Bad Channels OADB from given path %s",fCustomBC.Data()));
    
    TFile *fbad=new TFile(Form("%s",fCustomBC.Data()),"read");
    if (!fbad || fbad->IsZombie())
    {
      AliFatal(Form("Input file was not found in the path provided: %s",fCustomBC.Data()));
      return 0;
    }
    
    if (fbad) delete fbad;
    
    contBC->InitFromFile(Form("%s",fCustomBC.Data()),"AliEMCALBadChannels");
  }
  else
  { // Else choose the one in the $ALICE_PHYSICS directory
    if (fDebugLevel>0) AliInfo("Loading Bad Channels OADB from /OADB/EMCAL");
    
    TFile *fbad=new TFile(AliDataFile::GetFileNameOADB(Form("EMCAL/EMCALBadChannels%s.root", fLoad1DBadChMap ? "_1D" : "")).data(),"read");
    if (!fbad || fbad->IsZombie())
    {
      AliFatal(Form("OADB/EMCAL/EMCALBadChannels%s.root was not found", fLoad1DBadChMap ? "_1D" : ""));
      return 0;
    }
      
    if (fbad) delete fbad;
    
    contBC->InitFromFile(AliDataFile::GetFileNameOADB(Form("EMCAL/EMCALBadChannels%s.root", fLoad1DBadChMap ? "_1D" : "")).data(),"AliEMCALBadChannels");
  }
  
  TObjArray *arrayBC=(TObjArray*)contBC->GetObject(runBC);
  if (!arrayBC)
  {
    AliError(Form("No external hot channel set for run number: %d", runBC));
    delete contBC;
    return 2;
  }
  if(fLoad1DBadChMap){
    TH1C *h = fEMCALRecoUtils->GetEMCALChannelStatusMap1D();
    if (h)
      delete h;
    h=(TH1C*)arrayBC->FindObject("EMCALBadChannelMap");

    if (!h)
    {
      AliError("Can not get EMCALBadChannelMap");
    }
    h->SetDirectory(0);
    fEMCALRecoUtils->SetEMCALChannelStatusMap1D(h);
  }else{
    Int_t sms = fEMCALGeo->GetEMCGeometry()->GetNumberOfSuperModules();
    for (Int_t i=0; i<sms; ++i)
    {
      TH2I *h = fEMCALRecoUtils->GetEMCALChannelStatusMap(i);
      if (h)
        delete h;
      h=(TH2I*)arrayBC->FindObject(Form("EMCALBadChannelMap_Mod%d",i));

      if (!h)
      {
        AliError(Form("Can not get EMCALBadChannelMap_Mod%d",i));
        continue;
      }
      h->SetDirectory(0);
      fEMCALRecoUtils->SetEMCALChannelStatusMap(i,h);
    }
  }

  delete contBC;
  
  return 1;
}

//_____________________________________________________
Int_t AliEMCALTenderSupply::InitRecalib()
{
  // Initialising recalibration factors.
  
  AliVEvent *event = GetEvent();

  if (!event)
    return 0;
  
  if (fDebugLevel>0)
    AliInfo("Initialising recalibration factors");
  
  // init default maps first
  if (!fEMCALRecoUtils->GetEMCALRecalibrationFactorsArray())
    fEMCALRecoUtils->InitEMCALRecalibrationFactors() ;

  Int_t runRC = event->GetRunNumber();
      
  AliOADBContainer *contRF=new AliOADBContainer("");
  if (fBasePath!="")
  { //if fBasePath specified in the ->SetBasePath()
    if (fDebugLevel>0)  AliInfo(Form("Loading Recalib OADB from given path %s",fBasePath.Data()));
    
    TFile *fRecalib= new TFile(Form("%s/EMCALRecalib%s.root",fBasePath.Data(), fLoad1DRecalibFactors ? "_1D" : ""),"read");
    if (!fRecalib || fRecalib->IsZombie())
    {
      AliFatal(Form("EMCALRecalib%s.root not found in %s", fLoad1DRecalibFactors ? "_1D" : "" ,fBasePath.Data()));
      return 0;
    }
    
    if (fRecalib) delete fRecalib;
    
    contRF->InitFromFile(Form("%s/EMCALRecalib%s.root",fBasePath.Data(), fLoad1DRecalibFactors ? "_1D" : ""),"AliEMCALRecalib");
  }
  else
  { // Else choose the one in the $ALICE_PHYSICS directory
    if (fDebugLevel>0)  AliInfo("Loading Recalib OADB from OADB/EMCAL");
    
    TFile *fRecalib= new TFile(AliDataFile::GetFileNameOADB(Form("EMCAL/EMCALRecalib%s.root", fLoad1DRecalibFactors ? "_1D" : "")).data(),"read");
    if (!fRecalib || fRecalib->IsZombie())
    {
      AliFatal(Form("OADB/EMCAL/EMCALRecalib%s.root was not found", fLoad1DRecalibFactors ? "_1D" : ""));
      return 0;
    }
    
    if (fRecalib) delete fRecalib;
      
    contRF->InitFromFile(AliDataFile::GetFileNameOADB(Form("EMCAL/EMCALRecalib%s.root", fLoad1DRecalibFactors ? "_1D" : "")).data(),"AliEMCALRecalib");
  }

  TObjArray *recal=(TObjArray*)contRF->GetObject(runRC);
  if (!recal)
  {
    AliError(Form("No Objects for run: %d",runRC));
    delete contRF;
    return 2;
  } 

  TObjArray *recalpass=(TObjArray*)recal->FindObject(fFilepass);
  if (!recalpass)
  {
    AliError(Form("No Objects for run: %d - %s",runRC,fFilepass.Data()));
    delete contRF;
    return 2;
  }

  TObjArray *recalib=(TObjArray*)recalpass->FindObject("Recalib");
  if (!recalib)
  {
    AliError(Form("No Recalib histos found for  %d - %s",runRC,fFilepass.Data())); 
    delete contRF;
    return 2;
  }

  if (fDebugLevel>0) recalib->Print();



  if(fLoad1DRecalibFactors){
    TH1S *h = fEMCALRecoUtils->GetEMCALChannelRecalibrationFactors1D();
    if (h)
      delete h;
    h=(TH1S*)recalib->FindObject("EMCALRecalFactors");

    if (!h)
    {
      AliError("Can not get EMCALRecalFactors");
    }
    h->SetDirectory(0);
    fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors1D(h);
  }else{
    Int_t sms = fEMCALGeo->GetEMCGeometry()->GetNumberOfSuperModules();
    for (Int_t i=0; i<sms; ++i) 
    {
      TH2F *h = fEMCALRecoUtils->GetEMCALChannelRecalibrationFactors(i);
      if (h)
        delete h;
      h = (TH2F*)recalib->FindObject(Form("EMCALRecalFactors_SM%d",i));
      if (!h) 
      {
        AliError(Form("Could not load EMCALRecalFactors_SM%d",i));
        continue;
      }
      h->SetDirectory(0);
      fEMCALRecoUtils->SetEMCALChannelRecalibrationFactors(i,h);
    }
  }
  
  delete contRF;
  
  return 1;
}

//_____________________________________________________
Int_t AliEMCALTenderSupply::InitRunDepRecalib()
{
  // Initialising recalibration factors.
  
  AliVEvent *event = GetEvent();
  
  if (!event) 
    return 0;

  // MF Disable temperature calibration for run2 for the moment
  // A memory leak in the AliOADBContainer was observed when the 
  // container is deleted. As consequence when processing multiple
  // runs the correction task leaks ~100 MB per run included in the 
  // job, where 70% of the memory leak comes from the temperature
  // calibration. Untill a fix is provided the temperature calibration
  // has to be disabled for run2 runs. Disabling has to be done before
  // opening the OADB container
  // 
  // For more information see https://alice.its.cern.ch/jira/browse/EMCAL-135
  if(fUseNewRunDepTempCalib){

    if (fDebugLevel>0)
      AliInfo("Initialising Run2 temperature recalibration factors");

    // init default maps first
    if (!fEMCALRecoUtils->GetEMCALRecalibrationFactorsArray())
      fEMCALRecoUtils->InitEMCALRecalibrationFactors() ;

    Int_t runRC = event->GetRunNumber();

    AliOADBContainer *contTemperature = new AliOADBContainer("");
    AliOADBContainer *contParams      = new AliOADBContainer("");
    if (fBasePath!="")
    { //if fBasePath specified in the ->SetBasePath()
      if (fDebugLevel>0)  AliInfo(Form("Loading Recalib OADB from given path %s",fBasePath.Data()));

      TFile *fRunDepTemperature       = new TFile(Form("%s/EMCALTemperatureCalibSM.root",fBasePath.Data()),"read");
      if (!fRunDepTemperature || fRunDepTemperature->IsZombie()) {
        AliFatal(Form("EMCALTemperatureCalibSM.root not found in %s",fBasePath.Data()));
        return 0;
      }
      if (fRunDepTemperature) delete fRunDepTemperature;

      TFile *fTemperatureCalibParam   = new TFile(Form("%s/EMCALTemperatureCalibParam.root",fBasePath.Data()),"read");
      if (!fTemperatureCalibParam || fTemperatureCalibParam->IsZombie()) {
        AliFatal(Form("EMCALTemperatureCalibParam.root not found in %s",fBasePath.Data()));
        return 0;
      }

      if (fTemperatureCalibParam) delete fTemperatureCalibParam;

      contTemperature->InitFromFile(Form("%s/EMCALTemperatureCalibSM.root",fBasePath.Data()),"AliEMCALTemperatureCalibSM");
      contParams->InitFromFile(Form("%s/EMCALTemperatureCalibParam.root",fBasePath.Data()),"AliEMCALTemperatureCalibParam");
    }
    else if (fCustomTempCalibSM!="" && fCustomTempCalibParams!="" )
    { //if fBasePath specified in the ->SetBasePath()
      if (fDebugLevel>0)  AliInfo(Form("Loading custom recalib OADB from given paths %s and %s",fCustomTempCalibSM.Data(),fCustomTempCalibParams.Data()));

      TFile *fRunDepTemperature       = new TFile(Form("%s",fCustomTempCalibSM.Data()),"read");
      if (!fRunDepTemperature || fRunDepTemperature->IsZombie()) {
        AliFatal(Form("SM-wise temperature file %s not found",fCustomTempCalibSM.Data()));
        return 0;
      }
      if (fRunDepTemperature) delete fRunDepTemperature;

      TFile *fTemperatureCalibParam   = new TFile(Form("%s",fCustomTempCalibParams.Data()),"read");
      if (!fTemperatureCalibParam || fTemperatureCalibParam->IsZombie()) {
        AliFatal(Form("Temp calibration params file %s not found",fCustomTempCalibParams.Data()));
        return 0;
      }

      if (fTemperatureCalibParam) delete fTemperatureCalibParam;

      contTemperature->InitFromFile(Form("%s",fCustomTempCalibSM.Data()),"AliEMCALTemperatureCalibSM");
      contParams->InitFromFile(Form("%s",fCustomTempCalibParams.Data()),"AliEMCALTemperatureCalibParam");
    }
    else
    { // Else choose the one in the $ALICE_PHYSICS directory
      if (fDebugLevel>0)  AliInfo("Loading Recalib OADB from OADB/EMCAL");
      TFile *fRunDepTemperature= new TFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTemperatureCalibSM.root").data(),"read");
      if (!fRunDepTemperature || fRunDepTemperature->IsZombie()) {
        AliFatal("OADB/EMCAL/EMCALTemperatureCalibSM.root was not found");
        return 0;
      }
      if (fRunDepTemperature) delete fRunDepTemperature;

      TFile *fTemperatureCalibParam= new TFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTemperatureCalibParam.root").data(),"read");
      if (!fTemperatureCalibParam || fTemperatureCalibParam->IsZombie()) {
        AliFatal("OADB/EMCAL/EMCALTemperatureCalibParam.root was not found");
        return 0;
      }
      if (fTemperatureCalibParam) delete fTemperatureCalibParam;


      contTemperature->InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTemperatureCalibSM.root").data(),"AliEMCALTemperatureCalibSM");
      contParams->InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTemperatureCalibParam.root").data(),"AliEMCALTemperatureCalibParam");
    }


    TObjArray *arrayParams=(TObjArray*)contParams->GetObject(runRC);
    if (!arrayParams)
    {
      AliError(Form("No external temperature calibration set for run number: %d", runRC));
      delete contParams;
      delete contTemperature;
      return 0;
    }
    TH1D *hRundepTemp = (TH1D*)contTemperature->GetObject(runRC);
    TH1F *hSlopeParam = (TH1F*)arrayParams->FindObject("hParamSlope");
    TH1F *hA0Param    = (TH1F*)arrayParams->FindObject("hParamA0");

    if (!hRundepTemp || !hSlopeParam || !hA0Param)
    {
      AliError(Form("Histogram missing for temperature calibration for run number: %d", runRC));
      delete contParams;
      delete contTemperature;
      return 0;
    }

    Int_t nSM = fEMCALGeo->GetEMCGeometry()->GetNumberOfSuperModules();
    for (Int_t ism=0; ism<nSM; ++ism)
    {
      Double_t temperature = hRundepTemp->GetBinContent(ism+1);
      for (Int_t icol=0; icol<48; ++icol)
      {
        for (Int_t irow=0; irow<24; ++irow)
        {
          Int_t absID     = fEMCALGeo->GetAbsCellIdFromCellIndexes(ism, irow, icol); // original calibration factor
          Float_t factor  = fEMCALRecoUtils->GetEMCALChannelRecalibrationFactor(ism,icol,irow);
          Float_t slope   = 0;
          Float_t offset  = 0;
          slope           = hSlopeParam->GetBinContent(absID+1);
          offset          = hA0Param->GetBinContent(absID+1);
          // Correction is the inverse of the calculated factor
          if(slope || offset)
            factor *= 1 / (offset + (slope * temperature) ); // correction dependent on T
          fEMCALRecoUtils->SetEMCALChannelRecalibrationFactor(ism,icol,irow,factor);
        } // columns
      } // rows
    } // SM loop

    delete contParams;
    delete contTemperature;

    return 1;

    // Run1 treatment with old calibration
  } else {
    if(event->GetRunNumber() > 197692){
      AliInfo("Temperature calibration could not be loaded. Please use useNewRWTempCalib = kTRUE in the AddTaskEMCALTender for Run2 data!");
      return 0;
    }

    if (fDebugLevel>0)
      AliInfo("Initialising Run1 recalibration factors");

    // init default maps first
    if (!fEMCALRecoUtils->GetEMCALRecalibrationFactorsArray())
      fEMCALRecoUtils->InitEMCALRecalibrationFactors() ;

    Int_t runRC = event->GetRunNumber();

    AliOADBContainer *contRF=new AliOADBContainer("");
    if (fBasePath!="")
    { //if fBasePath specified in the ->SetBasePath()
      if (fDebugLevel>0)  AliInfo(Form("Loading Recalib OADB from given path %s",fBasePath.Data()));

      TFile *fRunDepRecalib= new TFile(Form("%s/EMCALTemperatureCorrCalib.root",fBasePath.Data()),"read");
      if (!fRunDepRecalib || fRunDepRecalib->IsZombie())
      {
        AliFatal(Form("EMCALTemperatureCorrCalib.root not found in %s",fBasePath.Data()));
        return 0;
      }

      if (fRunDepRecalib) delete fRunDepRecalib;

      contRF->InitFromFile(Form("%s/EMCALTemperatureCorrCalib.root",fBasePath.Data()),"AliEMCALRunDepTempCalibCorrections");
    }
    else
    { // Else choose the one in the $ALICE_PHYSICS directory
      if (fDebugLevel>0)  AliInfo("Loading Recalib OADB from OADB/EMCAL");

      TFile *fRunDepRecalib= new TFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTemperatureCorrCalib.root").data(),"read");
      if (!fRunDepRecalib || fRunDepRecalib->IsZombie())
      {
        AliFatal("OADB/EMCAL/EMCALTemperatureCorrCalib.root was not found");
        return 0;
      }

      if (fRunDepRecalib) delete fRunDepRecalib;

      contRF->InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTemperatureCorrCalib.root").data(),"AliEMCALRunDepTempCalibCorrections");
    }

    TH1S *rundeprecal=(TH1S*)contRF->GetObject(runRC);

    if (!rundeprecal)
    {
      AliWarning(Form("No TemperatureCorrCalib Objects for run: %d",runRC));
      // let's get the closest runnumber instead then..
      Int_t lower = 0;
      Int_t ic = 0;
      Int_t maxEntry = contRF->GetNumberOfEntries();

      while ((ic < maxEntry) && (contRF->UpperLimit(ic) < runRC)) {
        lower = ic;
        ic++;
      }

      Int_t closest = lower;
      if ((ic<maxEntry) &&
    (contRF->LowerLimit(ic)-runRC) < (runRC - contRF->UpperLimit(lower))) {
    closest = ic;
      }

      AliWarning(Form("TemperatureCorrCalib Objects found closest id %d from run: %d", closest, contRF->LowerLimit(closest)));
      rundeprecal = (TH1S*) contRF->GetObjectByIndex(closest);
    }

    Int_t nSM = fEMCALGeo->GetEMCGeometry()->GetNumberOfSuperModules();
    Int_t nbins = rundeprecal->GetNbinsX();

    // Avoid use of Run1 param for Run2
    if(nSM > 12 && nbins < 12288)
    {
      AliError(Form("Total SM is %d but T corrections available for %d channels, skip Init of T recalibration factors",nSM,nbins));

      delete contRF;

      return 2;
    }

    if (fDebugLevel>0) rundeprecal->Print();

    for (Int_t ism=0; ism<nSM; ++ism)
    {
      for (Int_t icol=0; icol<48; ++icol)
      {
        for (Int_t irow=0; irow<24; ++irow)
        {
          Float_t factor = fEMCALRecoUtils->GetEMCALChannelRecalibrationFactor(ism,icol,irow);

          Int_t absID = fEMCALGeo->GetAbsCellIdFromCellIndexes(ism, irow, icol); // original calibration factor
          factor *= rundeprecal->GetBinContent(absID) / 10000. ; // correction dependent on T

          fEMCALRecoUtils->SetEMCALChannelRecalibrationFactor(ism,icol,irow,factor);
        } // columns
      } // rows
    } // SM loop

    delete contRF;

    return 1;
  }
}


//_____________________________________________________
Int_t AliEMCALTenderSupply::InitTimeCalibration()
{
  // Initialising bad channel maps
  AliVEvent *event = GetEvent();

  if (!event) 
    return 0;
  
  if (fDebugLevel>0) 
    AliInfo("Initialising time calibration map");
  
  // init default maps first
  if (!fEMCALRecoUtils->GetEMCALTimeRecalibrationFactorsArray())
    fEMCALRecoUtils->InitEMCALTimeRecalibrationFactors() ;

  Int_t runBC = event->GetRunNumber();
  
  AliOADBContainer *contBC = new AliOADBContainer("");
  if (fBasePath!="")
  { //if fBasePath specified in the ->SetBasePath()
    if (fDebugLevel>0) AliInfo(Form("Loading time calibration OADB from given path %s",fBasePath.Data()));
    
    TFile *fbad=new TFile(Form("%s/EMCALTimeCalib.root",fBasePath.Data()),"read");
    if (!fbad || fbad->IsZombie())
    {
      AliFatal(Form("EMCALTimeCalib.root was not found in the path provided: %s",fBasePath.Data()));
      return 0;
    }  
    
    if (fbad) delete fbad;
    
    contBC->InitFromFile(Form("%s/EMCALTimeCalib.root",fBasePath.Data()),"AliEMCALTimeCalib");    
  } 
  else 
  { // Else choose the one in the $ALICE_PHYSICS directory
    if (fDebugLevel>0) AliInfo("Loading time calibration OADB from OADB/EMCAL");
    
    TFile *fbad=new TFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTimeCalib.root").data(),"read");
    if (!fbad || fbad->IsZombie())
    {
      AliFatal("OADB/EMCAL/EMCALTimeCalib.root was not found");
      return 0;
    }  
      
    if (fbad) delete fbad;
    
    contBC->InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTimeCalib.root").data(),"AliEMCALTimeCalib");
  }
  
  TObjArray *arrayBC=(TObjArray*)contBC->GetObject(runBC);
  if (!arrayBC)
  {
    AliError(Form("No external time calibration set for run number: %d", runBC));
    delete contBC;
    return 2; 
  }
  
  // Here, it looks for a specific pass
  TString pass = fFilepass;
  if (fFilepass=="spc_calo") pass ="pass1";
  TObjArray *arrayBCpass=(TObjArray*)arrayBC->FindObject(pass);
  if (!arrayBCpass)
  {
    AliError(Form("No external time calibration set for: %d -%s", runBC,pass.Data()));
    delete contBC;
    return 2; 
  }

  if (fDebugLevel>0) arrayBCpass->Print();

  if(!fDoMergedBCs){
    for(Int_t i = 0; i < 4; i++)
    {
      TH1F *h = (TH1F*)fEMCALRecoUtils->GetEMCALChannelTimeRecalibrationFactors(i);
      if (h)
        delete h;
    
      h = (TH1F*)arrayBCpass->FindObject(Form("hAllTimeAvBC%d",i));

      if (!h)
      {
        AliError(Form("Can not get hAllTimeAvBC%d",i));
        continue;
      }
   
      // Shift parameters for bc0 and bc1 in this pass
      if ( fFilepass=="spc_calo" && (i==0 || i==1) ) 
      {
        for(Int_t icell = 0; icell < h->GetNbinsX(); icell++) 
          h->SetBinContent(icell,h->GetBinContent(icell)-100);
      }
    
      h->SetDirectory(0);
      fEMCALRecoUtils->SetEMCALChannelTimeRecalibrationFactors(i,h);
    }
  }else{

    TH1S *h = (TH1S*)fEMCALRecoUtils->GetEMCALChannelTimeRecalibrationFactors(0);
    if (h)
      delete h;
    
    h = (TH1S*)arrayBCpass->FindObject("hAllTimeAv"); //only HG cells

    if (!h)
      AliError("Can not get hAllTimeAv");
    
    h->SetDirectory(0);
    fEMCALRecoUtils->SetEMCALChannelTimeRecalibrationFactors(0,h);

  }
  
  delete contBC;
  
  return 1;  
}

//_____________________________________________________
Int_t AliEMCALTenderSupply::InitTimeCalibrationL1Phase()
{
  // Initialising run-by-run L1 phase in time calibration maps
  AliVEvent *event = GetEvent();

  if (!event) 
    return 0;
  
  if (fDebugLevel>0) 
    AliInfo("Initialising run-by-run L1 phase in time calibration map");
  
  // init default maps first
  if (!fEMCALRecoUtils->GetEMCALL1PhaseInTimeRecalibrationArray())
    fEMCALRecoUtils->InitEMCALL1PhaseInTimeRecalibration() ;

  Int_t runBC = event->GetRunNumber();
  
  AliOADBContainer *contBC = new AliOADBContainer("");
  if (fBasePath!="")
  { //if fBasePath specified in the ->SetBasePath()
    if (fDebugLevel>0) AliInfo(Form("Loading time calibration OADB from given path %s",fBasePath.Data()));
    
    TFile *timeFile=new TFile(Form("%s/EMCALTimeL1PhaseCalib.root",fBasePath.Data()),"read");
    if (!timeFile || timeFile->IsZombie())
    {
      AliFatal(Form("EMCALTimeL1PhaseCalib.root was not found in the path provided: %s",fBasePath.Data()));
      return 0;
    }  
    
    if (timeFile) delete timeFile;
    
    contBC->InitFromFile(Form("%s/EMCALTimeL1PhaseCalib.root",fBasePath.Data()),"AliEMCALTimeL1PhaseCalib");
  }
  else
  { // Else choose the one in the $ALICE_PHYSICS directory
    if (fDebugLevel>0) AliInfo("Loading L1 phase in time calibration OADB from OADB/EMCAL");
    
    TFile *timeFile=new TFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTimeL1PhaseCalib.root").data(),"read");
    if (!timeFile || timeFile->IsZombie())
    {
      AliFatal("OADB/EMCAL/EMCALTimeL1PhaseCalib.root was not found");
      return 0;
    }
      
    if (timeFile) delete timeFile;
    
    contBC->InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTimeL1PhaseCalib.root").data(),"AliEMCALTimeL1PhaseCalib");
  }
  
  TObjArray *arrayBC=(TObjArray*)contBC->GetObject(runBC);
  if (!arrayBC)
  {
    AliError(Form("No external L1 phase in time calibration set for run number: %d", runBC));
    delete contBC;
    return 2; 
  }
  
  // Only 1 L1 phase correction possible, except special cases
  TString pass = "pass1";
  
  if ( fFilepass=="muon_calo_pass1" && fRun > 209121 && fRun < 244284 ) 
    pass = "pass0";//period LHC15a-m

  TObjArray *arrayBCpass=(TObjArray*)arrayBC->FindObject(pass);
  if (!arrayBCpass)
  {
    AliError(Form("No external L1 phase in time calibration set for: %d -%s", runBC,pass.Data()));
    delete contBC;
    return 2; 
  }

  if (fDebugLevel>0) arrayBCpass->Print();


  TH1C *h = fEMCALRecoUtils->GetEMCALL1PhaseInTimeRecalibrationForAllSM();
  if (h) delete h;
    
  h = (TH1C*)arrayBCpass->FindObject(Form("h%d",runBC));
    
  if (!h) {
    AliFatal(Form("There is no calibration histogram h%d for this run",runBC));
  }
  h->SetDirectory(0);
  fEMCALRecoUtils->SetEMCALL1PhaseInTimeRecalibrationForAllSM(h);
  
  delete contBC;
  
  return 1;  
}

//_____________________________________________________
void AliEMCALTenderSupply::UpdateCells()
{
  //Remove bad cells from the cell list
  //Recalibrate energy and time cells 
  //This is required for later reclusterization

  AliVEvent *event = GetEvent();

  if (!event) return ;
  
  AliVCaloCells *cells = event->GetEMCALCells();
  Int_t bunchCrossNo = event->GetBunchCrossNumber();

  fEMCALRecoUtils->RecalibrateCells(cells, bunchCrossNo); 
  
  // remove exotic cells - loop through cells and zero the exotic ones
  // just like with bad cell rejection in reco utils (inside RecalibrateCells)
  if (fRejectExoticCells)
  {
    Short_t  absId  =-1;
    Double_t ecell = 0;
    Double_t tcell = 0;
    Double_t efrac = 0;
    Int_t  mclabel = -1;
    Bool_t   isExot = kFALSE;
  
    // loop through cells
    Int_t nEMcell  = cells->GetNumberOfCells() ;  
    for (Int_t iCell = 0; iCell < nEMcell; iCell++) 
    { 
      cells->GetCell(iCell, absId, ecell, tcell, mclabel, efrac);
    
      isExot = fEMCALRecoUtils->IsExoticCell(absId, cells, bunchCrossNo); 
      // zero if exotic
      if (isExot)
        cells->SetCell(iCell, absId, 0.0, -1.0, mclabel, efrac);
    } // cell loop
  } // reject exotic cells

  cells->Sort();
}

//_____________________________________________________
TString AliEMCALTenderSupply::GetBeamType()
{
  // Get beam type : pp-AA-pA
  // ESDs have it directly, AODs get it from hardcoded run number ranges
  
  AliVEvent *event = GetEvent();

  if (!event) { 
    AliError("Couldn't retrieve event!");
    return "";
  }

  TString beamType;

  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event);
  if (esd) {
    const AliESDRun *run = esd->GetESDRun();
    beamType = run->GetBeamType();
  }
  else
  {
    Int_t runNumber = event->GetRunNumber();
    if ((runNumber >= 136851 && runNumber <= 139517)  // LHC10h
  || (runNumber >= 166529 && runNumber <= 170593))  // LHC11h
    {
      beamType = "A-A";
    }
    else 
    {
      beamType = "p-p";
    }
  }

  return beamType;    
}

//_____________________________________________________
Int_t AliEMCALTenderSupply::InitRecParam()
{
  // Init reco params if not yet exist (probably shipped by the user already)

  if (fRecParam != 0)
    return 2;

  TString beamType = GetBeamType();

  // set some default reco params
  fRecParam = new AliEMCALRecParam();
  fRecParam->SetClusteringThreshold(0.100);
  fRecParam->SetMinECut(0.050);
  
  if (!fCalibrateTimeParamAvailable) {
    fRecParam->SetTimeCut(250*1.e-9);
    fRecParam->SetTimeMin(425*1.e-9);
    fRecParam->SetTimeMax(825*1.e-9);
  } else {
    fRecParam->SetTimeCut(100*1.e-9);
    fRecParam->SetTimeMin(-50*1.e-9);
    fRecParam->SetTimeMax(50*1.e-9);
  }
  
  if (beamType == "A-A") {
    fRecParam->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv2);
  } else {
    fRecParam->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv1);
  }

  return 0;
}

//_____________________________________________________
Bool_t AliEMCALTenderSupply::InitClusterization()
{
  // Initialing clusterization/unfolding algorithm and set all the needed parameters.
  
  AliVEvent *event = GetEvent();

  if (!event) 
    return kFALSE;

  if (fDebugLevel>0) 
    AliInfo(Form("Initialising reclustering parameters: Clusterizer type: %d",fRecParam->GetClusterizerFlag()));
  
  //---setup clusterizer
  if (fClusterizer) {
    // avoid deleting fDigitsArr
    fClusterizer->SetDigitsArr(0);
    delete fClusterizer;
    fClusterizer = 0;
  }
  if     (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv1)
    fClusterizer = new AliEMCALClusterizerv1 (fEMCALGeo);
  else if (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv2) 
    fClusterizer = new AliEMCALClusterizerv2(fEMCALGeo);
  else if (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerNxN) 
  {
    AliEMCALClusterizerNxN *clusterizer = new AliEMCALClusterizerNxN(fEMCALGeo);
    clusterizer->SetNRowDiff(fRecParam->GetNRowDiff());
    clusterizer->SetNColDiff(fRecParam->GetNColDiff());
    fClusterizer = clusterizer;
  } 
  else 
  {
    AliFatal(Form("Clusterizer < %d > not available", fRecParam->GetClusterizerFlag()));
    return kFALSE;
  }
  
  // Set the clustering parameters
  fClusterizer->SetECAClusteringThreshold(fRecParam->GetClusteringThreshold());
  fClusterizer->SetECALogWeight          (fRecParam->GetW0()                 );
  fClusterizer->SetMinECut               (fRecParam->GetMinECut()            );    
  fClusterizer->SetUnfolding             (fRecParam->GetUnfold()             );
  fClusterizer->SetECALocalMaxCut        (fRecParam->GetLocMaxCut()          );
  fClusterizer->SetTimeCut               (fRecParam->GetTimeCut()            );
  fClusterizer->SetTimeMin               (fRecParam->GetTimeMin()            );
  fClusterizer->SetTimeMax               (fRecParam->GetTimeMax()            );
  fClusterizer->SetInputCalibrated       (kTRUE                              );
  fClusterizer->SetJustClusters          (kTRUE                              );  
  
  // In case of unfolding after clusterization is requested, set the corresponding parameters
  if (fRecParam->GetUnfold()) 
  {
    for (Int_t i = 0; i < 8; ++i) 
    {
      fClusterizer->SetSSPars(i, fRecParam->GetSSPars(i));
    }
    for (Int_t i = 0; i < 3; ++i)
    {
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
  
  AliVEvent *event = GetEvent();

 if (!event)
    return;
    
  // In case of MC productions done before aliroot tag v5-02-Rev09
  // assing the cluster label to all the cells belonging to this cluster
  // very rough
  if (fSetCellMCLabelFromCluster || fSetCellMCLabelFromEdepFrac)
  {
    for (Int_t i = 0; i < fgkTotalCellNumber; i++)
    {
      fCellLabels      [i] =-1 ;
      fOrgClusterCellId[i] =-1 ;
    }
    
    Int_t nClusters = event->GetNumberOfCaloClusters();
    for (Int_t i = 0; i < nClusters; i++)
    {
      AliVCluster *clus =  event->GetCaloCluster(i);
      
      if (!clus) continue;
      
      if (!clus->IsEMCAL()) continue ;
      
      Int_t      label = clus->GetLabel();
      UShort_t * index = clus->GetCellsAbsId() ;
      
      for(Int_t icell=0; icell < clus->GetNCells(); icell++)
      {
        if(!fSetCellMCLabelFromEdepFrac) 
          fCellLabels[index[icell]] = label;
        
        fOrgClusterCellId[index[icell]] = i ; // index of the original cluster
      } // cell in cluster loop
    }// cluster loop
  }

  fDigitsArr->Clear("C");
  AliVCaloCells *cells = event->GetEMCALCells();
  Int_t ncells = cells->GetNumberOfCells();
  for (Int_t icell = 0, idigit = 0; icell < ncells; ++icell) 
  {
    Double_t amp=0, cellTime=0, efrac = 0;
    Float_t cellAmplitude=0;
    Short_t  cellNumber=0;
    Int_t mcLabel=-1;

    if (cells->GetCell(icell, cellNumber, amp, cellTime, mcLabel, efrac) != kTRUE)
      break;

    cellAmplitude = amp; // compilation problem

    // Do not add if energy already too low (some cells set to 0 if bad channels)
    if (cellAmplitude < fRecParam->GetMinECut())
      continue;

    // If requested, do not include exotic cells
   if (fEMCALRecoUtils->IsExoticCell(cellNumber,cells,event->GetBunchCrossNumber())) 
      continue;
    
    if(!fSetCellMCLabelFromEdepFrac)
    {
      if      (fSetCellMCLabelFromCluster) mcLabel = fCellLabels[cellNumber];
      else if (fRemapMCLabelForAODs      ) RemapMCLabelForAODs(mcLabel);
    }
    
    //
    // New way to set the cell MC labels, 
    // valid only for MC productions with aliroot > v5-07-21
    //
    TArrayI labeArr(0);
    TArrayF eDepArr(0);
    Int_t nLabels = 0;
    
    if(fSetCellMCLabelFromEdepFrac && fOrgClusterCellId[cellNumber] >= 0) // index can be negative if noisy cell that did not form cluster 
    {
      mcLabel = -1;

      fCellLabels[cellNumber] = idigit;
      
      Int_t iclus = fOrgClusterCellId[cellNumber];
      
      if(iclus < 0)
      {
        AliInfo("Negative original cluster index, skip \n");
        continue;
      }
      
      AliVCluster *clus = event->GetCaloCluster(iclus);
      
      
      for(Int_t icluscell=0; icluscell < clus->GetNCells(); icluscell++ )
      {
        if(cellNumber != clus->GetCellAbsId(icluscell)) continue ;
        
        // Select the MC label contributing, only if enough energy deposition

        fEMCALRecoUtils->RecalculateCellLabelsRemoveAddedGenerator(cellNumber, clus, GetMCEvent(), cellAmplitude, labeArr, eDepArr);
        nLabels = labeArr.GetSize();
      }

      if(cellAmplitude <= 0.01) continue ; // accept if > 10 MeV
    } // recalculate MC labels 

    
    if (mcLabel > 0 && efrac < 1e-6) efrac = 1;
    
    AliEMCALDigit* digit = new((*fDigitsArr)[idigit]) AliEMCALDigit(mcLabel, mcLabel, cellNumber,
                                                                    cellAmplitude, (Float_t)cellTime,
                                                                    AliEMCALDigit::kHG,idigit, 0, 0, efrac*cellAmplitude);
    
    if(nLabels > 0)
    {
      digit->SetListOfParents(nLabels,labeArr.GetArray(),eDepArr.GetArray());
    }
    
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
  
  AliVEvent *event = GetEvent();

  if (!event)
    return;
  
  TClonesArray *clus = dynamic_cast<TClonesArray*>(event->FindListObject("caloClusters"));
  if (!clus) 
    clus = dynamic_cast<TClonesArray*>(event->FindListObject("CaloClusters"));
  if (!clus) 
  {
    AliError(" Null pointer to calo clusters array, returning");
    return;
  }
    
  // Before destroying the orignal list, assign to the rec points the MC labels
  // of the original clusters, if requested
  if (fSetCellMCLabelFromCluster == 2) 
    SetClustersMCLabelFromOriginalClusters() ;

  Int_t nents = clus->GetEntriesFast();
  for (Int_t i=0; i < nents; ++i) 
  {
    AliVCluster *c = dynamic_cast<AliVCluster*>(clus->At(i));
    if (!c)
      continue;
    if (c->IsEMCAL())
    {
      delete clus->RemoveAt(i);
    }
  }
  
  clus->Compress();
  
  RecPoints2Clusters(clus);
}

//_____________________________________________________
void AliEMCALTenderSupply::RecPoints2Clusters(TClonesArray *clus)
{
  // Convert AliEMCALRecoPoints to AliESDCaloClusters/AliAODCaloClusters.
  // Cluster energy, global position, cells and their amplitude fractions are restored.
  
  AliVEvent *event = GetEvent();

  if (!event)
    return;

  Int_t ncls = fClusterArr->GetEntriesFast();
  for(Int_t i=0, nout=clus->GetEntriesFast(); i < ncls; ++i) 
  {
    AliEMCALRecPoint *recpoint = static_cast<AliEMCALRecPoint*>(fClusterArr->At(i));
    
    Int_t ncellsTrue = 0;
    const Int_t ncells = recpoint->GetMultiplicity();
    UShort_t   absIds[ncells];  
    Double32_t ratios[ncells];
    Int_t   *dlist = recpoint->GetDigitsList();
    Float_t *elist = recpoint->GetEnergiesList();
    for (Int_t c = 0; c < ncells; ++c) 
    {
      AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->At(dlist[c]));
      absIds[ncellsTrue] = digit->GetId();
      ratios[ncellsTrue] = elist[c]/digit->GetAmplitude();
      if (ratios[ncellsTrue] < 0.001) 
        continue;
      ++ncellsTrue;
    }
    
    if (ncellsTrue < 1) 
    {
      AliWarning("Skipping cluster with no cells");
      continue;
    }
    
    // calculate new cluster position
    TVector3 gpos;
    recpoint->GetGlobalPosition(gpos);
    Float_t g[3];
    gpos.GetXYZ(g);
    
    AliVCluster *c = static_cast<AliVCluster*>(clus->New(nout++));
    c->SetID(nout-1); 
    c->SetType(AliVCluster::kEMCALClusterv1);
    c->SetE(recpoint->GetEnergy());
    c->SetPosition(g);
    c->SetNCells(ncellsTrue);
    c->SetDispersion(recpoint->GetDispersion());
    c->SetEmcCpvDistance(-1);            //not yet implemented
    c->SetChi2(-1);                      //not yet implemented
    c->SetTOF(recpoint->GetTime()) ;     //time-of-flight
    c->SetNExMax(recpoint->GetNExMax()); //number of local maxima
    Float_t elipAxis[2];
    recpoint->GetElipsAxis(elipAxis);
    c->SetM02(elipAxis[0]*elipAxis[0]) ;
    c->SetM20(elipAxis[1]*elipAxis[1]) ;
    c->SetCellsAbsId(absIds);
    c->SetCellsAmplitudeFraction(ratios);

    //
    // MC labels
    //
    Int_t    parentMult   = 0;
    Int_t   *parentList   = recpoint->GetParents(parentMult);
    Float_t *parentListDE = recpoint->GetParentsDE();  // deposited energy

    if (parentMult > 0)
    {
      c->SetLabel(parentList, parentMult);
      c->SetClusterMCEdepFractionFromEdepArray(parentListDE);
    }
    
    //
    // Set the cell energy deposition fraction map:
    //
    if( parentMult > 0 && fSetCellMCLabelFromEdepFrac )
    {
      UInt_t * mcEdepFracPerCell = new UInt_t[ncellsTrue];
      
      // Get the digit that originated this cell cluster
      
      for(Int_t icell = 0; icell < ncellsTrue ; icell++) 
      {
        Int_t   idigit  = fCellLabels[absIds[icell]];
        
        const AliEMCALDigit * dig = (const AliEMCALDigit*)fDigitsArr->At(idigit);
        
        // Find the 4 MC labels that contributed to the cluster and their 
        // deposited energy in the current digit
        
        mcEdepFracPerCell[icell] = 0; // init
        
        Int_t  nparents   = dig->GetNiparent();
        if ( nparents > 0 ) 
        {
          Int_t   digLabel   =-1 ; 
          Float_t edep       = 0 ;
          Float_t edepTot    = 0 ;
          Float_t mcEDepFrac[4] = {0,0,0,0};
          
          // all parents in digit
          for ( Int_t jndex = 0 ; jndex < nparents ; jndex++ ) 
          { 
            digLabel = dig->GetIparent (jndex+1);
            edep     = dig->GetDEParent(jndex+1);
            edepTot += edep;
            
            if       ( digLabel == parentList[0] ) mcEDepFrac[0] = edep; 
            else  if ( digLabel == parentList[1] ) mcEDepFrac[1] = edep;
            else  if ( digLabel == parentList[2] ) mcEDepFrac[2] = edep;
            else  if ( digLabel == parentList[3] ) mcEDepFrac[3] = edep;
          } // all prarents in digit
          
          // Divide energy deposit by total deposited energy
          // Do this only when deposited energy is significant, use 10 MeV although 50 MeV should be expected
          if(edepTot > 0.01) 
          {
            mcEdepFracPerCell[icell] = c->PackMCEdepFraction(mcEDepFrac);
          }
        } // at least one parent label in digit
      } // cell in cluster loop
      
      c->SetCellsMCEdepFractionMap(mcEdepFracPerCell);
      
      delete [] mcEdepFracPerCell;
      
    } // at least one parent in cluster, do the cell primary packing
  }
}

//___________________________________________________________
void AliEMCALTenderSupply::RemapMCLabelForAODs(Int_t & label)
{
  // MC label for Cells not remapped after ESD filtering, do it here.
  
  if (label < 0) return;
  
  AliAODEvent  * evt = dynamic_cast<AliAODEvent*> (GetEvent()) ;
  if (!evt) return ;
  
  TClonesArray * arr = dynamic_cast<TClonesArray*>(evt->FindListObject("mcparticles")) ;
  if (!arr) return ;
  
  if (label < arr->GetEntriesFast())
  {
    AliAODMCParticle * particle = dynamic_cast<AliAODMCParticle *>(arr->At(label));
    if (!particle) return ;
    
    if (label == particle->Label()) return ; // label already OK
  }
  
  // loop on the particles list and check if there is one with the same label
  for (Int_t ind = 0; ind < arr->GetEntriesFast(); ind++)
  {
    AliAODMCParticle * particle = dynamic_cast<AliAODMCParticle *>(arr->At(ind));
    if (!particle) continue ;
    
    if (label == particle->Label())
    {
      label = ind;
      return;
    }
  }
  
  label = -1;
}

//_____________________________________________________________________________________________
void AliEMCALTenderSupply::SetClustersMCLabelFromOriginalClusters()
{
  // Get the original clusters that contribute to the new rec point cluster,
  // assign the labels of such clusters to the new rec point cluster.
  // Only approximatedly valid  when output are V1 clusters, or input/output clusterizer
  // are the same handle with care
  // Copy from same method in AliAnalysisTaskEMCALClusterize, but here modify the recpoint and
  // not the output calocluster
  
  Int_t ncls = fClusterArr->GetEntriesFast();
  for(Int_t irp=0; irp < ncls; ++irp)
  {
    TArrayI clArray(300) ; //Weird if more than a few clusters are in the origin ...
    clArray.Reset();
    Int_t nClu = 0;
    Int_t nLabTotOrg = 0;
    Float_t emax = -1;
    Int_t idMax = -1;
    
    AliEMCALRecPoint *clus = static_cast<AliEMCALRecPoint*>(fClusterArr->At(irp));
    
    //Find the clusters that originally had the cells
    const Int_t ncells = clus->GetMultiplicity();
    Int_t *digList     = clus->GetDigitsList();
    
    for (Int_t iLoopCell = 0 ; iLoopCell < ncells ; iLoopCell++)
    {
      AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->At(digList[iLoopCell]));
      Int_t idCell = digit->GetId();

      if (idCell>=0)
      {
        Int_t idCluster = fOrgClusterCellId[idCell];
        Bool_t set = kTRUE;
        for (Int_t icl =0; icl < nClu; icl++)
        {
          if (((Int_t)clArray.GetAt(icl))==-1) continue;
          if (idCluster == ((Int_t)clArray.GetAt(icl))) set = kFALSE;
        }
        if (set && idCluster >= 0)
        {
          clArray.SetAt(idCluster,nClu++);
          nLabTotOrg+=(GetEvent()->GetCaloCluster(idCluster))->GetNLabels();
                    
          //Search highest E cluster
          AliVCluster * clOrg = GetEvent()->GetCaloCluster(idCluster);
          if (emax < clOrg->E())
          {
            emax  = clOrg->E();
            idMax = idCluster;
          }
        }
      }
    }// cell loop
        
    // Put the first in the list the cluster with highest energy
    if (idMax != ((Int_t)clArray.GetAt(0))) // Max not at first position
    {
      Int_t maxIndex = -1;
      Int_t firstCluster = ((Int_t)clArray.GetAt(0));
      for (Int_t iLoopCluster = 0 ; iLoopCluster < nClu ; iLoopCluster++)
      {
        if (idMax == ((Int_t)clArray.GetAt(iLoopCluster))) maxIndex = iLoopCluster;
      }
      
      if (firstCluster >=0 && idMax >=0)
      {
        clArray.SetAt(idMax,0);
        clArray.SetAt(firstCluster,maxIndex);
      }
    }
    
    // Get the labels list in the original clusters, assign all to the new cluster
    TArrayI clMCArray(nLabTotOrg) ;
    clMCArray.Reset();
    
    Int_t nLabTot = 0;
    for (Int_t iLoopCluster = 0 ; iLoopCluster < nClu ; iLoopCluster++)
    {
      Int_t idCluster = (Int_t) clArray.GetAt(iLoopCluster);
      AliVCluster * clOrg = GetEvent()->GetCaloCluster(idCluster);
      Int_t nLab = clOrg->GetNLabels();
      
      for (Int_t iLab = 0 ; iLab < nLab ; iLab++)
      {
        Int_t lab = clOrg->GetLabelAt(iLab) ;
        if (lab>=0)
        {
          Bool_t set = kTRUE;
          for(Int_t iLabTot =0; iLabTot < nLabTot; iLabTot++)
          {
            if (lab == ((Int_t)clMCArray.GetAt(iLabTot))) set = kFALSE;
          }
          if (set) clMCArray.SetAt(lab,nLabTot++);
        }
      }
    }// cluster loop
    
    // Set the final list of labels to rec point
    
    Int_t *labels = new Int_t[nLabTot];
    for(Int_t il = 0; il < nLabTot; il++) labels[il] = clMCArray.GetArray()[il];
    clus->SetParents(nLabTot,labels);
    
  } // rec point array
}

//_____________________________________________________
void AliEMCALTenderSupply::GetPass()
{
  // Get passx from filename
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  fInputTree = mgr->GetTree();
  
  if (!fInputTree) 
  {
    AliError("Pointer to tree = 0, returning");
    return;
  }
  
  fInputFile = fInputTree->GetCurrentFile();
  if (!fInputFile) {
    AliError("Null pointer input file, returning");
    return;
  }
  
  TString fname(fInputFile->GetName());
  if      (fname.Contains("pass1")) fFilepass = TString("pass1");
  else if (fname.Contains("pass2")) fFilepass = TString("pass2");
  else if (fname.Contains("pass3")) fFilepass = TString("pass3");
  else if (fname.Contains("pass4")) fFilepass = TString("pass4");
  else if (fname.Contains("pass5")) fFilepass = TString("pass5");
  else if (fname.Contains("LHC11c") &&
           fname.Contains("spc_calo")) fFilepass = TString("spc_calo");
  else if (fname.Contains("calo") || fname.Contains("high_lumi"))

  {
    //printf("AliEMCALTenderSupply::GetPass() - Path contains <calo> or <high-lumi>, set as <pass1>\n");
    fFilepass = TString("pass1");
  }
  else if (fname.Contains("LHC14a1a"))
  {
    fCalibrateEnergy     = kTRUE;
    fUseAutomaticRecalib = kTRUE;
    
    AliInfo("Energy calibration activated for this MC production!");
    
    fFilepass = TString("LHC14a1a");
  }
  else
  {
    AliError(Form("Pass number string not found: %s", fname.Data()));
    return;            
  }
}
