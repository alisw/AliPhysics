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

// --- Root ---
#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TObjArray.h>
#include <TString.h>
#include <TTree.h>

// --- AliRoot ---
#include "AliAODCaloCluster.h"
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALAfterBurnerUF.h"
#include "AliEMCALCalibData.h"
#include "AliEMCALClusterizerNxN.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALClusterizerv2.h"
#include "AliEMCALClusterizerFixedWindow.h"
#include "AliEMCALFixedWindowClusterInfo.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecParam.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"

#include "AliAnalysisTaskEMCALClusterizeFast.h"

ClassImp(AliAnalysisTaskEMCALClusterizeFast)

//________________________________________________________________________
AliAnalysisTaskEMCALClusterizeFast::AliAnalysisTaskEMCALClusterizeFast() 
  : AliAnalysisTaskSE(), 
    fRun(-1),
    fDigitsArr(0),       
    fClusterArr(0),       
    fRecParam(0),
    fClusterizer(0),
    fUnfolder(0),
    fJustUnfold(kFALSE),
    fGeomName(),
    fGeomMatrixSet(kFALSE), 
    fLoadGeomMatrices(kFALSE),
    fOCDBpath(),
    fCalibData(0),
    fPedestalData(0),
    fOutputAODBranch(0),
    fOutputAODBrName(),
    fRecoUtils(0),
    fLoadCalib(0),
    fLoadPed(0),
    fAttachClusters(0),
    fRecalibOnly(0),
    fSubBackground(0),
    fCreatePattern(0),
    fOverwrite(1),
    fNewClusterArrayName("newCaloClusters"),
    fNPhi(4),
    fNEta(4),
    fShiftPhi(2),
    fShiftEta(2),
    fTRUShift(0),
    fStoreAdditionalInformation(0)
{ 
  // Constructor
}

//________________________________________________________________________
AliAnalysisTaskEMCALClusterizeFast::AliAnalysisTaskEMCALClusterizeFast(const char *name) 
  : AliAnalysisTaskSE(name), 
    fRun(-1),
    fDigitsArr(0),       
    fClusterArr(0),       
    fRecParam(new AliEMCALRecParam),
    fClusterizer(0),
    fUnfolder(0),
    fJustUnfold(kFALSE),
    fGeomName("EMCAL_FIRSTYEARV1"),
    fGeomMatrixSet(kFALSE), 
    fLoadGeomMatrices(kFALSE),
    fOCDBpath(),
    fCalibData(0),
    fPedestalData(0),
    fOutputAODBranch(0),
    fOutputAODBrName(),
    fRecoUtils(0),
    fLoadCalib(0),
    fLoadPed(0),
    fAttachClusters(0),
    fRecalibOnly(0),
    fSubBackground(0),
    fCreatePattern(0),
    fOverwrite(1),
    fNewClusterArrayName("newCaloClusters"),
    fNPhi(4),
    fNEta(4),
    fShiftPhi(2),
    fShiftEta(2),
    fTRUShift(0),
    fStoreAdditionalInformation(0)
{ 
  // Constructor

  fBranchNames     = "ESD:AliESDHeader.,AliESDRun.,EMCALCells. AOD:header,emcalCells";
  for(Int_t i = 0; i < 12; ++i) 
    fGeomMatrix[i] = 0;
}

//________________________________________________________________________
AliAnalysisTaskEMCALClusterizeFast::~AliAnalysisTaskEMCALClusterizeFast()
{
  // Destructor.

  delete fDigitsArr; 
  delete fClusterizer;
  delete fUnfolder;   
  delete fRecoUtils;
}

//-------------------------------------------------------------------
void AliAnalysisTaskEMCALClusterizeFast::UserCreateOutputObjects()
{
  // Create output objects.

  if (!fOutputAODBrName.IsNull()) {
    fOutputAODBranch = new TClonesArray("AliAODCaloCluster", 0);
    fOutputAODBranch->SetName(fOutputAODBrName);
    AddAODBranch("TClonesArray", &fOutputAODBranch);
    AliInfo(Form("Created Branch: %s",fOutputAODBrName.Data()));
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::UserExec(Option_t *) 
{
  // Main loop, called for each event

  // remove the contents of output list set in the previous event 
  if (fOutputAODBranch)
    fOutputAODBranch->Clear("C");

  AliESDEvent *esdevent = dynamic_cast<AliESDEvent*>(InputEvent());
  AliAODEvent *aodevent = dynamic_cast<AliAODEvent*>(InputEvent());

  if (!esdevent&&!aodevent) {
    Error("UserExec","Event not available");
    return;
  }

  LoadBranches();

  UInt_t offtrigger = 0;
  if (esdevent) {
    UInt_t mask1 = esdevent->GetESDRun()->GetDetectorsInDAQ();
    UInt_t mask2 = esdevent->GetESDRun()->GetDetectorsInReco();
    Bool_t desc1 = (mask1 >> 18) & 0x1;
    Bool_t desc2 = (mask2 >> 18) & 0x1;
    if (desc1==0 || desc2==0) { //AliDAQ::OfflineModuleName(180=="EMCAL"
      AliError(Form("EMCAL not in DAQ/RECO: %u (%u)/%u (%u)", 
		    mask1, esdevent->GetESDRun()->GetDetectorsInReco(),
		    mask2, esdevent->GetESDRun()->GetDetectorsInDAQ()));
      return;
    }
    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    offtrigger = ((AliInputEventHandler*)(am->GetInputEventHandler()))->IsEventSelected();
  } else if (aodevent) {
    offtrigger =  aodevent->GetHeader()->GetOfflineTrigger();
  }
  if (offtrigger & AliVEvent::kFastOnly) {
    AliWarning(Form("EMCAL not in fast only partition"));
    return;
  }
  
  Init();

  if (fJustUnfold) {
    AliWarning("Unfolding not implemented");
    return;
  }

  FillDigitsArray();

  if (fRecalibOnly) {
    UpdateCells();
    return; // not requested to run clusterizer
  }

  Clusterize();
  UpdateCells();
  UpdateClusters();
  
  if (fStoreAdditionalInformation)
    StoreAdditionalInformation();

  if (fOutputAODBranch)
    RecPoints2Clusters(fOutputAODBranch);
}

//________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::StoreAdditionalInformation()
{
  if (fClusterizer->ClassName() != TString("AliEMCALClusterizerFixedWindow"))
    return;
  
  TString addInfoName(fNewClusterArrayName);
  addInfoName.Append("_AbsIds");
  
  AliEMCALFixedWindowClusterInfo *clusInfo = dynamic_cast<AliEMCALFixedWindowClusterInfo*>(InputEvent()->FindListObject(addInfoName));
  
  if(!clusInfo)
  {
    AliEMCALClusterizerFixedWindow *clusterizer = dynamic_cast<AliEMCALClusterizerFixedWindow*> (fClusterizer);
    if (!clusterizer)
      return;
    clusInfo = clusterizer->GetClustersInfo();
    if (!clusInfo)
      return;
    clusInfo->SetName(addInfoName);
    InputEvent()->AddObject(clusInfo);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::Clusterize()
{
  // Clusterize

  if (fSubBackground) {
    fClusterizer->SetInputCalibrated(kTRUE);   
    fClusterizer->SetCalibrationParameters(0);
  }

  fClusterizer->Digits2Clusters("");
  if (fSubBackground) {
    if (fCalibData) {
      fClusterizer->SetInputCalibrated(kFALSE);   
      fClusterizer->SetCalibrationParameters(fCalibData);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::FillDigitsArray()
{
  if (fCreatePattern)
  {
    
    AliEMCALGeometry *fGeom = AliEMCALGeometry::GetInstance(fGeomName);
    
    fDigitsArr->Clear("C");
    Int_t maxd = fGeom->GetNCells() / 4;
    
    for (Int_t idigit = 0; idigit < maxd; idigit++)
    {
      if (idigit % 24 == 12) idigit += 12;
      AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->New(idigit));
      digit->SetId(idigit * 4);
      digit->SetTime(600);
      digit->SetTimeR(600);
      digit->SetIndexInList(idigit);
      digit->SetType(AliEMCALDigit::kHG);
      digit->SetAmplitude(0.1);
      
    }
    
  }
  else
	{
		
    // Fill digits from cells.
    
    fDigitsArr->Clear("C");
    AliVCaloCells *cells = InputEvent()->GetEMCALCells();
    Double_t avgE = 0; // for background subtraction
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
      if (fRecalibOnly||fSubBackground) {
        Float_t energy = cellAmplitude;
        Float_t time    = cellTime;
        fClusterizer->Calibrate(energy,time,cellNumber);
        digit->SetAmplitude(energy);
        avgE += energy;
      } else {
        digit->SetAmplitude(cellAmplitude);
      }
      idigit++;
    }
    
    if (fSubBackground) {
      avgE /= AliEMCALGeometry::GetInstance(fGeomName)->GetNumberOfSuperModules()*48*24;
      Int_t ndigis = fDigitsArr->GetEntries();
      for (Int_t i = 0; i < ndigis; ++i) {
        AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->At(i));
        Double_t energy = digit->GetAmplitude() - avgE;
        if (energy<=0.001) {
          digit->SetAmplitude(0);
        } else {
          digit->SetAmplitude(energy);
        }
      }
    }
	}
}

//________________________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::RecPoints2Clusters(TClonesArray *clus)
{
  // Cluster energy, global position, cells and their amplitude fractions are restored.

  Bool_t esdobjects = 0;
  if (strcmp(clus->GetClass()->GetName(),"AliESDCaloCluster")==0)
    esdobjects = 1;

  AliVCaloCells *cells = InputEvent()->GetEMCALCells();
  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance(fGeomName);
  
  Int_t Ncls = fClusterArr->GetEntriesFast();
  for(Int_t i=0, nout=clus->GetEntries(); i < Ncls; ++i) {
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
    
    AliVCluster *c = static_cast<AliVCluster*>(clus->New(nout++));
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
    if (fRecoUtils && fRecoUtils->IsBadChannelsRemovalSwitchedOn()) {
      fRecoUtils->RecalculateClusterDistanceToBadChannel(geom, cells, c);
    } else {
      if (fPedestalData) 
        recpoint->EvalDistanceToBadChannels(fPedestalData);
      c->SetDistanceToBadChannel(recpoint->GetDistanceToBadTower()); 
    }

    if (esdobjects) {
      AliESDCaloCluster *cesd = static_cast<AliESDCaloCluster*>(c);
      cesd->SetCellsAbsId(absIds);
      cesd->SetCellsAmplitudeFraction(ratios);
    } else {
      AliAODCaloCluster *caod = static_cast<AliAODCaloCluster*>(c);
      caod->SetCellsAbsId(absIds);
      caod->SetCellsAmplitudeFraction(ratios);
    }
  }
 
  AliESDEvent *esdevent = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esdevent)
    return;
  if (!fRecoUtils)
    return;

  AliAnalysisManager::GetAnalysisManager()->LoadBranch("Tracks");
  fRecoUtils->FindMatches(esdevent,clus);
  
  if (!esdobjects) {
    Int_t Nclus = clus->GetEntries();
    for(Int_t i=0; i < Nclus; ++i) {
      AliAODCaloCluster *c = static_cast<AliAODCaloCluster*>(clus->At(i));
      Int_t trackIndex = fRecoUtils->GetMatchedTrackIndex(i);
      if(trackIndex >= 0) {
        Float_t dR, dZ;
        fRecoUtils->GetMatchedResiduals(i,dR, dZ);
        c->AddTrackMatched(0x0); //esdevent->GetTrack(trackIndex));
        c->SetTrackDistance(dR,dZ); // not implemented
        c->SetEmcCpvDistance(dR);
        c->SetChi2(dZ);
        if(DebugLevel() > 1) 
          AliInfo(Form("Matched Track index %d to new cluster %d\n",trackIndex,i));
      }
    }
  } else {
    Int_t Nclus = clus->GetEntries();
    for(Int_t i=0; i < Nclus; ++i) {
      AliESDCaloCluster *c = static_cast<AliESDCaloCluster*>(clus->At(i));
      Int_t trackIndex = fRecoUtils->GetMatchedTrackIndex(i);
      if(trackIndex >= 0) {
        Float_t dR, dZ;
        fRecoUtils->GetMatchedResiduals(i,dR, dZ);
        c->SetTrackDistance(dR,dZ);
        c->SetEmcCpvDistance(dR); //to be consistent with AODs
        c->SetChi2(dZ);           //to be consistent with AODs
        TArrayI tm(1,&trackIndex);
        c->AddTracksMatched(tm);
        if(DebugLevel() > 1) 
          AliInfo(Form("Matched Track index %d to new cluster %d\n",trackIndex,i));
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::UpdateCells()
{
  // Update cells in case re-calibration was done.

  if (!fCalibData&&!fSubBackground)
    return;

  AliVCaloCells *cells = InputEvent()->GetEMCALCells();
  Int_t ncells = cells->GetNumberOfCells();
  Int_t ndigis = fDigitsArr->GetEntries();

  if (ncells!=ndigis) {
    cells->DeleteContainer();
    cells->CreateContainer(ndigis);
  }
  for (Int_t idigit = 0; idigit < ndigis; ++idigit) {
    AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->At(idigit));
    Double_t cellAmplitude = digit->GetCalibAmp();
    Short_t cellNumber = digit->GetId();
    Double_t cellTime = digit->GetTime();
    cells->SetCell(idigit, cellNumber, cellAmplitude, cellTime);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::UpdateClusters()
{
  // Update cells in case re-calibration was done.
  
  if (!fAttachClusters)
    return;
  
  TClonesArray *clus;
  
  if (fOverwrite)
  {
    clus = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("caloClusters"));
    if (!clus)
      clus = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("CaloClusters"));
    if(!clus)
      return;
    
    Int_t nents = clus->GetEntries();
    for (Int_t i=0;i<nents;++i) {
      AliVCluster *c = static_cast<AliVCluster*>(clus->At(i));
      if (!c)
        continue;
      if (c->IsEMCAL())
        delete clus->RemoveAt(i);
    }
  }
  else
  {
    clus = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fNewClusterArrayName));
    if(!clus)
    {
      clus = new TClonesArray("AliESDCaloCluster");
      clus->SetName(fNewClusterArrayName);
      InputEvent()->AddObject(clus);
    }
    else 
    {
      
      Int_t nents = clus->GetEntries();
      for (Int_t i=0;i<nents;++i) {
        AliVCluster *c = static_cast<AliVCluster*>(clus->At(i));
        if (!c)
          continue;
        
        delete clus->RemoveAt(i);
      }
      clus->Compress();
    }
  }
  
  RecPoints2Clusters(clus);
}

//________________________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::Init()
{
  //Select clusterization/unfolding algorithm and set all the needed parameters

  AliVEvent * event = InputEvent();
  if (!event) {
    AliWarning("Event not available!!!");
    return;
  }

  if (event->GetRunNumber()==fRun)
    return;
  fRun = event->GetRunNumber();

  if (fJustUnfold){
    // init the unfolding afterburner 
    delete fUnfolder;
    fUnfolder = new AliEMCALAfterBurnerUF(fRecParam->GetW0(),fRecParam->GetLocMaxCut());
    return;
  }

  AliEMCALGeometry *geometry = AliEMCALGeometry::GetInstance(fGeomName);
  if (!geometry) {
    AliFatal("Geometry not available!!!");
    return;
  }

  if (!fGeomMatrixSet) {
    if (fLoadGeomMatrices) {
      for(Int_t mod=0; mod < geometry->GetNumberOfSuperModules(); ++mod) {
        if(fGeomMatrix[mod]){
          if(DebugLevel() > 2) 
            fGeomMatrix[mod]->Print();
          geometry->SetMisalMatrix(fGeomMatrix[mod],mod);  
        }
      }
    } else { // get matrix from file (work around bug in aliroot)
      for(Int_t mod=0; mod < geometry->GetEMCGeometry()->GetNumberOfSuperModules(); ++mod) {
        const TGeoHMatrix *gm = 0;
        AliESDEvent *esdevent = dynamic_cast<AliESDEvent*>(event);
        if (esdevent) {
          gm = esdevent->GetEMCALMatrix(mod);
        } else {
          AliAODHeader *aodheader = dynamic_cast<AliAODHeader*>(event->GetHeader());
          if (aodheader) {
            gm = aodheader->GetEMCALMatrix(mod);
          }
        }
        if (gm) {
          if(DebugLevel() > 2) 
            gm->Print();
          geometry->SetMisalMatrix(gm,mod);
        }
      }
    }
    fGeomMatrixSet=kTRUE;
  }
  
  // setup digit array if needed
  if (!fDigitsArr) {
    fDigitsArr = new TClonesArray("AliEMCALDigit", 1000);
    fDigitsArr->SetOwner(1);
  }

  // then setup clusterizer
  delete fClusterizer;
  if     (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv1)
    fClusterizer = new AliEMCALClusterizerv1(geometry);
  else if(fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerNxN) {
   AliEMCALClusterizerNxN *clusterizer = new AliEMCALClusterizerNxN(geometry);
   clusterizer->SetNRowDiff(fRecParam->GetNRowDiff());
   clusterizer->SetNColDiff(fRecParam->GetNColDiff());
    fClusterizer = clusterizer;
  } else if(fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv2) 
    fClusterizer = new AliEMCALClusterizerv2(geometry);
  else if(fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerFW){
    AliEMCALClusterizerFixedWindow *clusterizer = new AliEMCALClusterizerFixedWindow(geometry);
    clusterizer->SetNphi(fNPhi);
    clusterizer->SetNeta(fNEta);
    clusterizer->SetShiftPhi(fShiftPhi);
    clusterizer->SetShiftEta(fShiftEta);
    clusterizer->SetTRUshift(fTRUShift);
    fClusterizer = clusterizer;
  }
  else{
    AliFatal(Form("Clusterizer < %d > not available", fRecParam->GetClusterizerFlag()));
  }
  fClusterizer->InitParameters(fRecParam);

  if ((!fCalibData&&fLoadCalib) || (!fPedestalData&&fLoadPed)) {
    AliCDBManager *cdb = AliCDBManager::Instance();
    if (!cdb->IsDefaultStorageSet() && !fOCDBpath.IsNull())
      cdb->SetDefaultStorage(fOCDBpath);
    if (fRun!=cdb->GetRun())
      cdb->SetRun(fRun);
  }
  if (!fCalibData&&fLoadCalib&&fRun>0) {
    AliCDBEntry *entry = static_cast<AliCDBEntry*>(AliCDBManager::Instance()->Get("EMCAL/Calib/Data"));
    if (entry) 
      fCalibData =  static_cast<AliEMCALCalibData*>(entry->GetObject());
    if (!fCalibData)
      AliFatal("Calibration parameters not found in CDB!");
  }
  if (!fPedestalData&&fLoadPed&&fRun>0) {
    AliCDBEntry *entry = static_cast<AliCDBEntry*>(AliCDBManager::Instance()->Get("EMCAL/Calib/Pedestals"));
    if (entry) 
      fPedestalData =  static_cast<AliCaloCalibPedestal*>(entry->GetObject());
  }
  if (fCalibData) {
    fClusterizer->SetInputCalibrated(kFALSE);   
    fClusterizer->SetCalibrationParameters(fCalibData);
  } else {
    fClusterizer->SetInputCalibrated(kTRUE);   
  }
  fClusterizer->SetCaloCalibPedestal(fPedestalData);
  fClusterizer->SetJustClusters(kTRUE);
  fClusterizer->SetDigitsArr(fDigitsArr);
  fClusterizer->SetOutput(0);
  fClusterArr = const_cast<TObjArray *>(fClusterizer->GetRecPoints());
}
