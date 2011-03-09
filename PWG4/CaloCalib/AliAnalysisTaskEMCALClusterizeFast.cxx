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
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecParam.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDEvent.h"
#include "AliLog.h"

#include "AliAnalysisTaskEMCALClusterizeFast.h"

ClassImp(AliAnalysisTaskEMCALClusterizeFast)

//________________________________________________________________________
AliAnalysisTaskEMCALClusterizeFast::AliAnalysisTaskEMCALClusterizeFast() 
  : AliAnalysisTaskSE(), 
    fRun(0),
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
    fLoadPed(0)
{ 
  // Constructor
}

//________________________________________________________________________
AliAnalysisTaskEMCALClusterizeFast::AliAnalysisTaskEMCALClusterizeFast(const char *name) 
  : AliAnalysisTaskSE(name), 
    fRun(0),
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
    fLoadPed(0)
{ 
  // Constructor

  fBranchNames     = "ESD:AliESDHeader.,EMCALCells. AOD:header,emcalCells";
  for(Int_t i = 0; i < 10; ++i) 
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
  
  Init();

  if (fJustUnfold) {
    AliWarning("Unfolding not implemented");
  } else {
    if (esdevent) 
      FillDigitsArray(esdevent);
    else 
      FillDigitsArray(aodevent);
    fClusterizer->Digits2Clusters("");
    if (esdevent && fRecoUtils)
      fRecoUtils->FindMatches(esdevent,fClusterArr);
    if (fOutputAODBranch) {
      RecPoints2Clusters();
    }
    if (esdevent) {
      UpdateCells(esdevent);
    } else {
      UpdateCells(aodevent);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::UpdateCells(AliAODEvent *event)
{
  // Update cells in case re-calibration was done.

  if (!fCalibData)
    return;

  AliAODCaloCells *cells = event->GetEMCALCells();
  Int_t ncells = cells->GetNumberOfCells();
  for (Int_t icell = 0, idigit = 0; icell < ncells; ++icell) {
    Double_t cellAmplitude=0, cellTime=0;
    Short_t cellNumber=0;
    if (cells->GetCell(icell, cellNumber, cellAmplitude, cellTime) != kTRUE)
      break;
    AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->At(idigit));
    cellAmplitude = digit->GetCalibAmp();
    cells->SetCell(icell, cellNumber, cellAmplitude, cellTime);
    idigit++;
  }
  cells->Sort();
}

//________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::UpdateCells(AliESDEvent *event)
{
  // Update cells in case re-calibration was done.

  if (!fCalibData)
    return;

  AliESDCaloCells *cells = event->GetEMCALCells();
  Int_t ncells = cells->GetNumberOfCells();
  for (Int_t icell = 0, idigit = 0; icell < ncells; ++icell) {
    Double_t cellAmplitude=0, cellTime=0;
    Short_t cellNumber=0;
    if (cells->GetCell(icell, cellNumber, cellAmplitude, cellTime) != kTRUE)
      break;
    AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->At(idigit));
    cellAmplitude = digit->GetCalibAmp();
    cells->SetCell(icell, cellNumber, cellAmplitude, cellTime);
    idigit++;
  }
  cells->Sort();
}

//________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::FillDigitsArray(AliAODEvent *event)
{
  // Fill digits from cells.

  fDigitsArr->Clear("C");
  AliAODCaloCells *cells = event->GetEMCALCells();
  Int_t ncells = cells->GetNumberOfCells();
  if (ncells>fDigitsArr->GetSize())
    fDigitsArr->Expand(2*ncells);
  for (Int_t icell = 0, idigit = 0; icell < ncells; ++icell) {
    Double_t cellAmplitude=0, cellTime=0;
    Short_t cellNumber=0;
    if (cells->GetCell(icell, cellNumber, cellAmplitude, cellTime) != kTRUE)
      break;
    AliEMCALDigit *digit = (AliEMCALDigit*) fDigitsArr->New(idigit);
    digit->SetId(cellNumber);
    digit->SetAmplitude(cellAmplitude);
    digit->SetTime(cellTime);
    digit->SetTimeR(cellTime);
    digit->SetIndexInList(idigit);
    digit->SetType(AliEMCALDigit::kHG);
    idigit++;
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::FillDigitsArray(AliESDEvent *event)
{
  // Fill digits from cells.

  fDigitsArr->Clear("C");
  AliESDCaloCells *cells = event->GetEMCALCells();
  Int_t ncells = cells->GetNumberOfCells();
  if (ncells>fDigitsArr->GetSize())
    fDigitsArr->Expand(2*ncells);
  for (Int_t icell = 0, idigit = 0; icell < ncells; ++icell) {
    Double_t cellAmplitude=0, cellTime=0;
    Short_t cellNumber=0;
    if (cells->GetCell(icell, cellNumber, cellAmplitude, cellTime) != kTRUE)
      break;
    AliEMCALDigit *digit = (AliEMCALDigit*) fDigitsArr->New(idigit);
    digit->SetId(cellNumber);
    digit->SetAmplitude(cellAmplitude);
    digit->SetTime(cellTime);
    digit->SetTimeR(cellTime);
    digit->SetIndexInList(idigit);
    digit->SetType(AliEMCALDigit::kHG);
    idigit++;
  }
}

//________________________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::RecPoints2Clusters()
{
  // Cluster energy, global position, cells and their amplitude fractions are restored.

  AliESDEvent *esdevent = dynamic_cast<AliESDEvent*>(InputEvent());

  Int_t Ncls = fClusterArr->GetEntriesFast();
  for(Int_t i=0, nout=0; i < Ncls; ++i) {
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
      if (ratios[ncells_true] > 0.001) 
        ++ncells_true;
    }
    
    if (ncells_true < 1) {
      AliWarning("Skipping cluster with no cells");
      continue;
    }
    
    // calculate new cluster position
    TVector3 gpos;
    Float_t g[3];

    recpoint->EvalGlobalPosition(fRecParam->GetW0(), fDigitsArr);
    recpoint->GetGlobalPosition(gpos);
    gpos.GetXYZ(g);
    
    AliAODCaloCluster *clus = static_cast<AliAODCaloCluster*>(fOutputAODBranch->New(nout++));
    clus->SetType(AliVCluster::kEMCALClusterv1);
    clus->SetE(recpoint->GetEnergy());
    clus->SetPosition(g);
    clus->SetNCells(ncells_true);
    clus->SetCellsAbsId(absIds);
    clus->SetCellsAmplitudeFraction(ratios);
    clus->SetDispersion(recpoint->GetDispersion());
    clus->SetChi2(-1);                      //not yet implemented
    clus->SetTOF(recpoint->GetTime()) ;     //time-of-flight
    clus->SetNExMax(recpoint->GetNExMax()); //number of local maxima
    Float_t elipAxis[2];
    recpoint->GetElipsAxis(elipAxis);
    clus->SetM02(elipAxis[0]*elipAxis[0]) ;
    clus->SetM20(elipAxis[1]*elipAxis[1]) ;
    clus->SetDistToBadChannel(recpoint->GetDistanceToBadTower()); 
    if (esdevent && fRecoUtils) {
      Int_t trackIndex = fRecoUtils->GetMatchedTrackIndex(i);
      if(trackIndex >= 0) {
        clus->AddTrackMatched(esdevent->GetTrack(trackIndex));
        if(DebugLevel() > 1) 
          AliInfo(Form("Matched Track index %d to new cluster %d\n",trackIndex,i));
      }
    }
  }
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

  AliCDBManager *cdb = AliCDBManager::Instance();
  if (!cdb->IsDefaultStorageSet() && !fOCDBpath.IsNull())
    cdb->SetDefaultStorage(fOCDBpath);
  if (fRun!=cdb->GetRun())
    cdb->SetRun(fRun);

  AliEMCALGeometry *geometry = AliEMCALGeometry::GetInstance(fGeomName);
  if (!geometry) {
    AliFatal("Geometry not available!!!");
    return;
  }

  if (!fGeomMatrixSet) {
    if (fLoadGeomMatrices) {
      for(Int_t mod=0; mod < (geometry->GetEMCGeometry())->GetNumberOfSuperModules(); ++mod) {
        if(fGeomMatrix[mod]){
          if(DebugLevel() > 2) 
            fGeomMatrix[mod]->Print();
          geometry->SetMisalMatrix(fGeomMatrix[mod],mod);  
        }
      }
    } else {
      for(Int_t mod=0; mod < geometry->GetEMCGeometry()->GetNumberOfSuperModules(); ++mod) {
        if(event->GetEMCALMatrix(mod)) {
          if(DebugLevel() > 2) 
            event->GetEMCALMatrix(mod)->Print();
          geometry->SetMisalMatrix(event->GetEMCALMatrix(mod),mod);
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
  else if(fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerNxN) 
    fClusterizer = new AliEMCALClusterizerNxN(geometry);
  else if(fRecParam->GetClusterizerFlag() > AliEMCALRecParam::kClusterizerNxN) {
   AliEMCALClusterizerNxN *clusterizer = new AliEMCALClusterizerNxN(geometry);
   clusterizer->SetNRowDiff(2);
   clusterizer->SetNColDiff(2);
   fClusterizer = clusterizer;
  } else{
    AliFatal(Form("Clusterizer < %d > not available", fRecParam->GetClusterizerFlag()));
  }
  fClusterizer->InitParameters(fRecParam);
  if (!fCalibData&&fLoadCalib) {
    AliCDBEntry *entry = static_cast<AliCDBEntry*>(AliCDBManager::Instance()->Get("EMCAL/Calib/Data"));
    if (entry) 
      fCalibData =  static_cast<AliEMCALCalibData*>(entry->GetObject());
    if (!fCalibData)
      AliFatal("Calibration parameters not found in CDB!");
  }
  if (!fPedestalData&&fLoadPed) {
    AliCDBEntry *entry = static_cast<AliCDBEntry*>(AliCDBManager::Instance()->Get("EMCAL/Calib/Pedestals"));
    if (entry) 
      fPedestalData =  static_cast<AliCaloCalibPedestal*>(entry->GetObject());
  }
  if (fCalibData) {
    fClusterizer->SetInputCalibrated(kFALSE);   
    fClusterizer->SetCalibrationParameters(fCalibData);
    fClusterizer->SetCaloCalibPedestal(fPedestalData);
  } else {
    fClusterizer->SetInputCalibrated(kTRUE);   
  }
  fClusterizer->SetJustClusters(kTRUE);
  fClusterizer->SetDigitsArr(fDigitsArr);
  fClusterizer->SetOutput(0);
  fClusterArr = const_cast<TObjArray *>(fClusterizer->GetRecPoints());
}
