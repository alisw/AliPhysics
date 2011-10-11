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
/* $Id$ */

//-----------------------------------------------------------------------------
/// \class AliMUONReconstructor
///
/// Implementation of AliReconstructor for MUON subsystem.
///
/// The clustering mode and the associated parameters can be changed through the
/// AliMUONRecoParam object set in the reconstruction macro or read from the CDB
/// (see methods in AliMUONRecoParam.h file for details)
///
/// Valid modes are :
///
/// SIMPLEFIT : use the AliMUONClusterFinderSimpleFit clusterizer
///
/// SIMPLEFITV3 : SIMPLEFIT with preclustering=PRECLUSTERV3
///
/// MLEM : use AliMUONClusterFinderMLEM and AliMUONPreClusterFinder for preclustering (default)
/// MLEMV2 : MLEM with preclustering=PRECLUSTERV2
/// MLEMV3 : MLEM with preclustering=PRECLUSTERV3
///
/// PRECLUSTER : use only AliMUONPreClusterFinder. Only for debug as
/// the produced clusters do not have a position, hence the tracking will not
/// work
/// PRECLUSTERV2 : another version of the preclustering
/// PRECLUSTERV3 : yet another version of the preclustering
///
/// COG : use AliMUONClusterFinderCOG clusterizer. Not really a production
/// option either, as center-of-gravity is generally not a good estimate
/// of the cluster position...
///
/// PEAKCOG : COG cluster finder around local maxima
/// PEAKFIT : fit around local maxima with up to 3 peaks, COG otherwise
///
/// NOCLUSTERING : bypass completely the clustering stage
///
/// ------
///
/// The behavior of the MUON reconstruction can also be changed, besides
/// the usual methods found in AliReconstruction (e.g. to disable tracking)
/// by using AliReconstruction::SetOption("MUON",options)
/// where options should be a space separated string.
///
/// Valid options are :
///
/// SAVEDIGITS : if you want to save in the TreeD the *calibrated* digits
///     that are used for the clustering
///
/// DIGITSTOREV1 : use the V1 implementation of the digitstore 
/// DIGITSTOREV2R : use the V2R implementation of the digitstore 
///
/// NOLOCALRECONSTRUCTION : for debug, to disable local reconstruction (and hence
/// "recover" old behavior)
///
/// TRIGGERDISABLE : disable the treatment of MUON trigger
///
/// ENABLEERRORLOGGING : enable error logging (this activates also warnings in RawStreamTracker)
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONReconstructor.h"

#include "AliMUONCalibrationData.h"
#include "AliMUONClusterFinderCOG.h"
#include "AliMUONClusterFinderMLEM.h"
#include "AliMUONClusterFinderSimpleFit.h"
#include "AliMUONClusterFinderPeakCOG.h"
#include "AliMUONClusterFinderPeakFit.h"
#include "AliMUONClusterStoreV1.h"
#include "AliMUONClusterStoreV2.h"
#include "AliMUONConstants.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONDigitMaker.h"
#include "AliMUONDigitStoreV1.h"
#include "AliMUONDigitStoreV2R.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONPreClusterFinder.h"
#include "AliMUONPreClusterFinderV2.h"
#include "AliMUONPreClusterFinderV3.h"
#include "AliMUONSimpleClusterServer.h"
#include "AliMUONTracker.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONTriggerStoreV1.h"
#include "AliMUONVClusterFinder.h"
#include "AliMUONVClusterServer.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTriggerElectronics.h"
#include "AliMUONTriggerUtilities.h"

#include "AliMpArea.h"
#include "AliMpCDB.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpSegmentation.h"

#include "AliRawReader.h"
#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliRunInfo.h"

#include <Riostream.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TTree.h>

/// \cond CLASSIMP
ClassImp(AliMUONReconstructor)
/// \endcond 

//_____________________________________________________________________________
AliMUONReconstructor::AliMUONReconstructor() : 
AliReconstructor(),
fDigitMaker(0x0),
fTransformer(new AliMUONGeometryTransformer()),
fDigitStore(0x0),
fTriggerCircuit(0x0),
fCalibrationData(0x0),
fDigitCalibrator(0x0),
fTriggerStore(0x0),
fTrackStore(0x0),
fClusterStore(0x0),
fTriggerProcessor(0x0),
fTriggerUtilities(0x0),
fClusterServers(),
fTrackers()
{
  /// normal ctor

  AliDebug(1,"");
  
  // Unload mapping objects
  // if they have been loaded from the obsolete OCDB mapping objects

  if ( AliMpDDLStore::Instance(false) ) {
    AliCDBManager::Instance()->UnloadFromCache("MUON/Calib/DDLStore");
    delete AliMpDDLStore::Instance();
  }  

  if ( AliMpSegmentation::Instance(false) ) { 
    AliCDBManager::Instance()->UnloadFromCache("MUON/Calib/Mapping");
    delete AliMpSegmentation::Instance();
  }  

  // Load mapping
  if ( ! AliMpCDB::LoadDDLStore() ) {
    AliFatal("Could not access mapping from OCDB !");
  }
  
  // Load geometry data
  fTransformer->LoadGeometryData();
 
  fClusterServers.SetOwner(kTRUE);
  fTrackers.SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUONReconstructor::~AliMUONReconstructor()
{
  /// dtor

  AliDebug(1,"");

  delete fDigitCalibrator;

  delete fDigitMaker;
  delete fDigitStore;
  delete fTransformer;
  delete fTriggerCircuit;
  delete fTriggerStore;
  delete fTrackStore;
  delete fClusterStore;
  delete fTriggerProcessor;
  delete fTriggerUtilities;

  delete AliMpSegmentation::Instance(false);
  delete AliMpDDLStore::Instance(false);  

  delete fCalibrationData;  
}

//_____________________________________________________________________________
void
AliMUONReconstructor::Calibrate(AliMUONVDigitStore& digitStore) const
{
  /// Calibrate the digitStore
  if (!fDigitCalibrator)
  {
    CreateCalibrator();
  }
  AliCodeTimerAuto(Form("%s::Calibrate(AliMUONVDigitStore*)",fDigitCalibrator->ClassName()),0)
  fDigitCalibrator->Calibrate(digitStore);      
}

//_____________________________________________________________________________
void
AliMUONReconstructor::ConvertDigits(AliRawReader* rawReader, 
                                    AliMUONVDigitStore* digitStore,
                                    AliMUONVTriggerStore* triggerStore) const
{
  /// Convert raw data into digit and trigger stores
  CreateDigitMaker();

  // Skip reconstruction if event is Calibration
  if ( GetRecoParam()->GetEventSpecie() == AliRecoParam::kCalib ) {
    digitStore->Clear(); // Remove possible digits from previous event
    triggerStore->Clear(); // Remove possible triggers from previous event
    AliInfo("Calibration event: do not convert digits");
    return;
  }
  
  AliCodeTimerStart(Form("%s::Raw2Digits(AliRawReader*,AliMUONVDigitStore*,AliMUONVTriggerStore*)",
                    fDigitMaker->ClassName()))
  fDigitMaker->Raw2Digits(rawReader,digitStore,triggerStore);
  AliCodeTimerStop(Form("%s::Raw2Digits(AliRawReader*,AliMUONVDigitStore*,AliMUONVTriggerStore*)",
                         fDigitMaker->ClassName()))
  Calibrate(*digitStore);
}

//_____________________________________________________________________________
void 
AliMUONReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{
   /// convert raw data into a digit tree
  AliCodeTimerAuto("",0)

  Bool_t alone = ( TriggerStore() == 0 );
  
  Bool_t ok = DigitStore()->Connect(*digitsTree,alone);
  if ( TriggerStore() ) 
  {
    ok = ok && TriggerStore()->Connect(*digitsTree,kFALSE);
  }
  
  if (!ok)
  {
    AliError("Could not make branches on TreeD");
  }
  else
  {
    ConvertDigits(rawReader,DigitStore(),TriggerStore());
    AliCodeTimerStart("Fill digits")
    digitsTree->Fill();
    AliCodeTimerStop("Fill digits")
    DigitStore()->Clear();
  }
}

//_____________________________________________________________________________
void
AliMUONReconstructor::CreateDigitMaker() const
{
  /// Create (and create if necessary) the digit maker
  if (fDigitMaker) return;

  AliCodeTimerAuto("",0)

  TString option = GetOption();
  
  Bool_t enableErrorLogging = kFALSE;

  if (option.Contains("ENABLEERRORLOGGING"))
  {
    enableErrorLogging = kTRUE;
  }

  fDigitMaker = new AliMUONDigitMaker(enableErrorLogging);
  option.ToUpper();

  // Always make trigger digits
  // (needed when calculating trigger chamber efficiency)
  fDigitMaker->SetMakeTriggerDigits(kTRUE);

  if ( GetRecoParam()->TryRecover() )
  {
    fDigitMaker->SetTryRecover(kTRUE);
  }
  else
  {
    fDigitMaker->SetTryRecover(kFALSE);    
  }
}

//_____________________________________________________________________________
void 
AliMUONReconstructor::CreateTriggerCircuit() const
{
  /// Return (and create if necessary) the trigger circuit object
  if (fTriggerCircuit) return;

  AliCodeTimerAuto("",0)

  fTriggerCircuit = new AliMUONTriggerCircuit(fTransformer);

}

//_____________________________________________________________________________
void 
AliMUONReconstructor::CreateTriggerUtilities() const
{
  /// Return (and create if necessary) the trigger utilities object
  if ( fTriggerUtilities ) return;
  
  AliCodeTimerAuto("",0)
  
  if ( ! fCalibrationData ) CreateCalibrationData();
  
  fTriggerUtilities = new AliMUONTriggerUtilities(fCalibrationData);
}

//_____________________________________________________________________________
AliTracker* 
AliMUONReconstructor::CreateTracker() const
{
  /// Create the MUONTracker object
  
  CreateTriggerCircuit();
  CreateTriggerUtilities();
  
  const AliMUONRecoParam* rp = GetRecoParam();
  
  Int_t es = rp->GetEventSpecie();
  
  AliTracker* tracker = static_cast<AliTracker*>(fTrackers.At(es));
  
  if (!tracker ) 
  {
    if ( ! rp->CombineClusterTrackReco() )
    {
      tracker = new AliMUONTracker(rp,
                                   0x0,
                                   *DigitStore(),
                                   fTransformer,
                                   fTriggerCircuit,
                                   fTriggerUtilities);

      AliInfo(Form("Created tracker %p for recoparam of type %s es=%d",
                   tracker,
                   AliRecoParam::GetEventSpecieName(AliRecoParam::Convert(rp->GetEventSpecie())),es));
    }
    else
    {
      
      tracker = new AliMUONTracker(rp,
                                   CreateClusterServer(*rp),
                                   *DigitStore(),
                                   fTransformer,
                                   fTriggerCircuit,
                                   fTriggerUtilities);

      AliInfo(Form("Created (combined) tracker %p for recoparam of type %s es=%d",
                   tracker,
                   AliRecoParam::GetEventSpecieName(AliRecoParam::Convert(rp->GetEventSpecie())),es));      
    }
    
    fTrackers.AddAtAndExpand(tracker,es);
  }
  
  
  return tracker;
}

//_____________________________________________________________________________
AliMUONVClusterFinder*
AliMUONReconstructor::CreateClusterFinder(const char* clusterFinderType)
{
  /// Create a given cluster finder instance
  
  AliCodeTimerAutoGeneral("",0)

  AliMUONVClusterFinder* clusterFinder(0x0);
  
  TString opt(clusterFinderType);
  opt.ToUpper();
  
  if ( strstr(opt,"PRECLUSTERV2") )
  {
    clusterFinder = new AliMUONPreClusterFinderV2;
  }    
  else if ( strstr(opt,"PRECLUSTERV3") )
  {
    clusterFinder = new AliMUONPreClusterFinderV3;
  }  
  else if ( strstr(opt,"PRECLUSTER") )
  {
    clusterFinder = new AliMUONPreClusterFinder;
  }  
  else if ( strstr(opt,"PEAKCOG") )
  {
    clusterFinder = new AliMUONClusterFinderPeakCOG(kFALSE,new AliMUONPreClusterFinder);
  }
  else if ( strstr(opt,"PEAKFIT") )
  {
    clusterFinder = new AliMUONClusterFinderPeakFit(kFALSE,new AliMUONPreClusterFinder);
  }
  else if ( strstr(opt,"COG") )
  {
    clusterFinder = new AliMUONClusterFinderCOG(new AliMUONPreClusterFinder);
  }  
  else if ( strstr(opt,"SIMPLEFITV3") )
  {
    clusterFinder = new AliMUONClusterFinderSimpleFit(new AliMUONClusterFinderCOG(new AliMUONPreClusterFinderV3));
  }
  else if ( strstr(opt,"SIMPLEFIT") )
  {
    clusterFinder = new AliMUONClusterFinderSimpleFit(new AliMUONClusterFinderCOG(new AliMUONPreClusterFinder));
  }
  else if ( strstr(opt,"MLEM:DRAW") )
  {
    clusterFinder = new AliMUONClusterFinderMLEM(kTRUE,new AliMUONPreClusterFinder);
  }
  else if ( strstr(opt,"MLEMV3") )
  {
    clusterFinder = new AliMUONClusterFinderMLEM(kFALSE,new AliMUONPreClusterFinderV3);
  } 
  else if ( strstr(opt,"MLEMV2") )
  {
    clusterFinder = new AliMUONClusterFinderMLEM(kFALSE,new AliMUONPreClusterFinderV2);
  } 
  else if ( strstr(opt,"MLEM") )
  {
    clusterFinder = new AliMUONClusterFinderMLEM(kFALSE,new AliMUONPreClusterFinder);
  } 
  else
  {
    AliErrorClass(Form("clustering mode \"%s\" does not exist",opt.Data()));
    return 0x0;
  }
  
  return clusterFinder;
}

//_____________________________________________________________________________
AliMUONVClusterServer*
AliMUONReconstructor::CreateClusterServer(const AliMUONRecoParam& rp) const
{
  /// Create cluster server
  
  AliCodeTimerAuto("",0);
  
  AliMUONVClusterServer* clusterServer = static_cast<AliMUONVClusterServer*>(fClusterServers.At(rp.GetEventSpecie()));
  
  if (!clusterServer )
  {
    AliMUONVClusterFinder* clusterFinder = CreateClusterFinder(rp.GetClusteringMode());
    
    if ( !clusterFinder ) return 0x0;
    
    clusterServer = new AliMUONSimpleClusterServer(clusterFinder,*fTransformer);

    AliInfo(Form("Created AliMUONSimpleClusterServer (%p) for specie %d with clustering = %s",
                 clusterServer,rp.GetEventSpecie(),clusterFinder->ClassName()));
    
    fClusterServers.AddAtAndExpand(clusterServer,rp.GetEventSpecie());
  }
  
  return clusterServer;
}

//_____________________________________________________________________________
void
AliMUONReconstructor::CreateCalibrationData() const
{
  /// Create the calibrator
  
  AliCodeTimerAuto("",0);
  
  Int_t runNumber = AliCDBManager::Instance()->GetRun();

  fCalibrationData = new AliMUONCalibrationData(runNumber);
  if ( !fCalibrationData->IsValid() )
  {
    AliError("Could not retrieve calibrations !");
    delete fCalibrationData;
    fCalibrationData = 0x0;
    return;
  }    
  
  // It is now time to check whether we have everything to proceed.
  // What we need depends on whether both tracker and trigger
  // are in the readout chain, and what specific "bad channel policy"
  // we use
  
  Bool_t kTracker(kFALSE);
  Bool_t kTrigger(kFALSE);
  
  const AliRunInfo* runInfo = GetRunInfo();
  if (!runInfo)
  {
    AliError("Could not get runinfo ?");
  }
  else
  {
    TString detectors(runInfo->GetActiveDetectors());
    if (detectors.Contains("MUONTRK")) kTracker=kTRUE;
    if (detectors.Contains("MUONTRG")) kTrigger=kTRUE;
  }
   
  AliInfo(Form("Run with MUON TRIGGER : %s and MUON TRACKER : %s",
               kTrigger ? "YES":"NO" , 
               kTracker ? "YES":"NO"));
  
  if ( kTracker ) 
  {
    // Check that we get all the calibrations we'll need
    if ( !fCalibrationData->Pedestals() ||
        !fCalibrationData->Gains() )
    {
      AliFatal(Form("Could not access all required calibration data (PED %p GAIN %p)",
                    fCalibrationData->Pedestals(),fCalibrationData->Gains()));      
    }
    
    if ( !fCalibrationData->HV() ) 
    {
      // Special treatment of HV. We only break if the values
      // are not there *AND* we cut on them.
      UInt_t mask = GetRecoParam()->PadGoodnessMask();
      TString smask(AliMUONPadStatusMaker::AsCondition(mask));
      if ( smask.Contains("HV") )
      {
        AliFatal("Could not access all required calibration data (HV)");      
      }
    } 
  }
}

//_____________________________________________________________________________
void
AliMUONReconstructor::CreateCalibrator() const
{
  /// Create the calibrator

  AliCodeTimerAuto("",0);

  if ( ! fCalibrationData )
    CreateCalibrationData();

  AliInfo("Calibration will occur.");

  TString opt(GetOption());
  opt.ToUpper();
  
  if ( strstr(opt,"NOSTATUSMAP") )
  {
    AliWarning("NOSTATUSMAP is obsolete");
  }

  fDigitCalibrator = new AliMUONDigitCalibrator(*fCalibrationData,GetRecoParam());
}

//_____________________________________________________________________________
void
AliMUONReconstructor::ResponseRemovingChambers(AliMUONVTriggerStore* triggerStore) const
{
  /// Update trigger information with informatins obtained after
  /// re-calculation of trigger response
  AliCodeTimerAuto("",0);

  if ( ! fCalibrationData )
    CreateCalibrationData();

  if ( ! fTriggerProcessor )
    fTriggerProcessor = new AliMUONTriggerElectronics(fCalibrationData);

  fTriggerProcessor->ResponseRemovingChambers(*triggerStore);
}

//_____________________________________________________________________________
AliMUONVDigitStore*
AliMUONReconstructor::DigitStore() const
{
  /// Return (and create if necessary) the digit container
  if (!fDigitStore) 
  {
    TString sopt(GetOption());
    sopt.ToUpper();
    
    AliInfo(Form("Options=%s",sopt.Data()));
    
    if ( sopt.Contains("DIGITSTOREV1") )
    {
      fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV1");
    }
    else if ( sopt.Contains("DIGITSTOREV2R") ) 
    {
      fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV2R");
    }
    else if ( sopt.Contains("DIGITSTOREV2S") ) 
    {
      fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV2S");
    }
    
    if (!fDigitStore) fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV2R");
    
    AliInfo(Form("Will use %s to store digits during reconstruction",fDigitStore->ClassName()));
  }
  return fDigitStore;
}

//_____________________________________________________________________________
void
AliMUONReconstructor::FillTreeR(AliMUONVTriggerStore* triggerStore,
                                TTree& clustersTree) const
{
  /// Write the trigger and cluster information into TreeR

  AliCodeTimerAuto("",0)

  Bool_t ok(kFALSE);
  Bool_t alone(kTRUE); // is trigger the only info in TreeR ?
  
   const AliMUONRecoParam* rp = GetRecoParam();
  
  if ( ! rp->CombineClusterTrackReco() )
  {
    alone = kFALSE; // we'll get both tracker and trigger information in TreeR
  }
  
  if ( triggerStore ) 
  {
    ResponseRemovingChambers(triggerStore);
    ok = triggerStore->Connect(clustersTree,alone);
    if (!ok)
    {
      AliError("Could not create triggerStore branches in TreeR");
    }
  }

  if ( !alone )
  {
    if (!fClusterStore)
    {
      fClusterStore = new AliMUONClusterStoreV2;
    }
    
    AliMUONVClusterServer* clusterServer = CreateClusterServer(*rp);
    
    TIter next(DigitStore()->CreateIterator());
    clusterServer->UseDigits(next,DigitStore());

    AliMpArea area;
    
    AliDebug(1,Form("Doing full clusterization in local reconstruction using %s ",clusterServer->ClassName()));
    
    for ( Int_t i = 0; i < AliMpConstants::NofTrackingChambers(); ++i ) 
    {
      if (rp->UseChamber(i))
      {
        if ( ( i == 6 || i == 7 )  && rp->BypassSt4() ) continue;
        if ( ( i == 8 || i == 9 )  && rp->BypassSt5() ) continue;
        
        clusterServer->Clusterize(i,*fClusterStore,area,rp);
      }
    }
    
    Bool_t cok = fClusterStore->Connect(clustersTree,alone);
    
    if (!cok) AliError("Could not connect clusterStore to clusterTree");
    
    AliDebug(1,Form("Number of clusters found = %d",fClusterStore->GetSize()));
    
    StdoutToAliDebug(1,fClusterStore->Print());
  }
         
  if (ok) // at least one type of branches created successfully
  {
    clustersTree.Fill();
  }
  
  if (fClusterStore) fClusterStore->Clear();
}

//_____________________________________________________________________________
const AliMUONRecoParam* AliMUONReconstructor::GetRecoParam()
{ 
  /// Get the recoparam from reconstruction
  return dynamic_cast<const AliMUONRecoParam*>(AliReconstructor::GetRecoParam(AliReconstruction::GetDetIndex("MUON"))); 
}


//_____________________________________________________________________________
Bool_t 
AliMUONReconstructor::HasDigitConversion() const
{
  /// We *do* have digit conversion, but we might advertise it only 
  /// if we want to save the digits.
  
  TString opt(GetOption());
  opt.ToUpper();
  if ( opt.Contains("SAVEDIGITS" ) && !opt.Contains("NOLOCALRECONSTRUCTION") )
  {
    return kTRUE;
  }
  else
  {
    return kFALSE;
  }
}

//_____________________________________________________________________________
void 
AliMUONReconstructor::Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const
{
  /// This method is called by AliReconstruction if HasLocalReconstruction()==kTRUE AND
  /// HasDigitConversion()==kFALSE
  
  if ( !clustersTree ) 
  {
    AliError("clustersTree is 0x0 !");
    return;
  }
  
  ConvertDigits(rawReader,DigitStore(),TriggerStore());

  FillTreeR(TriggerStore(),*clustersTree);
}

//_____________________________________________________________________________
void 
AliMUONReconstructor::Reconstruct(TTree* digitsTree, TTree* clustersTree) const
{
  /// This method is called by AliReconstruction if HasLocalReconstruction()==kTRUE
  /// AND HasDigitConversion()==kTRUE
  
  AliCodeTimerAuto("",0)
  
  AliDebug(1,"");
  
  if (!digitsTree || !clustersTree) 
  {
    AliError(Form("Tree is null : digitsTree=%p clustersTree=%p",
                  digitsTree,clustersTree));
    return;
  }

  if (!fDigitStore)
  {
    fDigitStore = AliMUONVDigitStore::Create(*digitsTree);
    if (!fDigitStore)
    {
      AliError(Form("Could not get DigitStore from %s",digitsTree->GetName()));
    }
    else
    {
      AliInfo(Form("Created %s from %s",fDigitStore->ClassName(),digitsTree->GetName()));
    }
  }
  if (!fTriggerStore)
  {
    fTriggerStore = AliMUONVTriggerStore::Create(*digitsTree);
    if (!fTriggerStore)
    {
      AliError(Form("Could not get TriggerStore from %s",digitsTree->GetName()));
    }
    else
    {
      AliInfo(Form("Created %s from %s",fTriggerStore->ClassName(),digitsTree->GetName()));
    }
  }
  
  if (!fTriggerStore && !fDigitStore)
  {
    AliError("No store at all. Nothing to do.");
    return;
  }
  
  // insure we start with empty stores
  if ( fDigitStore ) 
  {
    fDigitStore->Clear(); 
    Bool_t alone = ( fTriggerStore ? kFALSE : kTRUE );
    Bool_t ok = fDigitStore->Connect(*digitsTree,alone);
    if (!ok)
    {
      AliError("Could not connect digitStore to digitsTree");
      return;
    }
  }
  if ( fTriggerStore ) 
  {
    fTriggerStore->Clear();
    Bool_t alone = ( fDigitStore ? kFALSE : kTRUE );
    Bool_t ok = fTriggerStore->Connect(*digitsTree,alone);
    if (!ok)
    {
      AliError("Could not connect triggerStore to digitsTree");
      return;
    }
  }
  
  digitsTree->GetEvent(0);
  
  if ( fDigitStore ) 
  {
    // Insure we got calibrated digits (if we reconstruct from pure simulated,
    // i.e. w/o going through raw data, this will be the case)
    TIter next(fDigitStore->CreateIterator());
    AliMUONVDigit* digit = static_cast<AliMUONVDigit*>(next());
    if (digit && !digit->IsCalibrated())
    {
      Calibrate(*fDigitStore);
    }
  }
    
  FillTreeR(fTriggerStore,*clustersTree);
}

//_____________________________________________________________________________
AliMUONVTriggerStore*
AliMUONReconstructor::TriggerStore() const
{
  /// Return (and create if necessary and allowed) the trigger container
  TString sopt(GetOption());
  sopt.ToUpper();
  
  if (sopt.Contains("TRIGGERDISABLE"))
  {
    delete fTriggerStore;
    fTriggerStore = 0x0;
  }
  else
  {
    if (!fTriggerStore)
    {
      fTriggerStore = new AliMUONTriggerStoreV1;
    }
  }
  return fTriggerStore;
}
