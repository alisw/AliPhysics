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

//_________________________________________________________________________
// This analysis provides a new list of clusters to be used in other analysis
//
// Author: Constantin Loizides, Salvatore Aiola
//         Adapted from analysis class from Deepa Thomas
//
// $Id$
//_________________________________________________________________________

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
    fAttachClusters(1),
    fRecalibOnly(0),
    fSubBackground(0),
    fOverwrite(0),
    fNewClusterArrayName("newCaloClusters"),
    fNPhi(4),
    fNEta(4),
    fShiftPhi(2),
    fShiftEta(2),
    fTRUShift(0),
    fInputCellType(kFEEData),
    fTrackName()
{ 
  // Constructor

  for(Int_t i = 0; i < 12; ++i) 
    fGeomMatrix[i] = 0;
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
    fGeomName(""),
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
    fAttachClusters(1),
    fRecalibOnly(0),
    fSubBackground(0),
    fOverwrite(0),
    fNewClusterArrayName("newCaloClusters"),
    fNPhi(4),
    fNEta(4),
    fShiftPhi(2),
    fShiftEta(2),
    fTRUShift(0),
    fInputCellType(kFEEData),
    fTrackName()
{ 
  // Constructor

  fBranchNames     = "ESD:AliESDHeader.,AliESDRun.,EMCALCells.,EMCALTrigger. AOD:header,emcalCells";
  for(Int_t i = 0; i < 12; ++i) 
    fGeomMatrix[i] = 0;
}

//________________________________________________________________________
AliAnalysisTaskEMCALClusterizeFast::~AliAnalysisTaskEMCALClusterizeFast()
{
  // Destructor.

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

  if (esdevent) {
    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    Bool_t mcmode = 0;
    if (am->GetMCtruthEventHandler())
      mcmode = 1;
    if (!mcmode) {
      if (offtrigger & AliVEvent::kFastOnly) {
        AliWarning(Form("EMCAL not in fast only partition"));
        return;
      }
    }
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

  if (fOutputAODBranch)
    RecPoints2Clusters(fOutputAODBranch);
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
  // Fill digits array

  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();

  fDigitsArr->Clear("C");

  switch (fInputCellType) {

  case kFEEData :
    {
      AliVCaloCells *cells = InputEvent()->GetEMCALCells();
      Double_t avgE        = 0; // for background subtraction
      const Int_t ncells   = cells->GetNumberOfCells();
      for (Int_t icell = 0, idigit = 0; icell < ncells; ++icell) {
	Double_t cellAmplitude=0, cellTime=0, cellEFrac = 0;
	Short_t  cellNumber=0;
	Int_t cellMCLabel=-1;
	if (cells->GetCell(icell, cellNumber, cellAmplitude, cellTime, cellMCLabel, cellEFrac) != kTRUE)
        break;
	AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->New(idigit));
	digit->SetId(cellNumber);
	digit->SetTime(cellTime);
	digit->SetTimeR(cellTime);
	digit->SetIndexInList(idigit);
	digit->SetType(AliEMCALDigit::kHG);
	if (fRecalibOnly||fSubBackground) {
	  Float_t energy = cellAmplitude;
	  Float_t time   = cellTime;
	  fClusterizer->Calibrate(energy,time,cellNumber);
	  digit->SetAmplitude(energy);
        avgE += energy;
	} else {
	  digit->SetAmplitude(cellAmplitude);
	}
	idigit++;
      }
      
      fDigitsArr->Sort();

      if (fSubBackground) {
	avgE /= geom->GetNumberOfSuperModules()*48*24;
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
    break;
    
  case kPattern :    
    {
      // Fill digits from a pattern
      Int_t maxd = geom->GetNCells() / 4;
      for (Int_t idigit = 0; idigit < maxd; idigit++){
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
    break;

  case kL0FastORs    : 
  case kL0FastORsTC  :
  case kL1FastORs    :
    {
      // Fill digits from FastORs
      
      AliVCaloTrigger *triggers = InputEvent()->GetCaloTrigger("EMCAL");
      
      if (!triggers || !(triggers->GetEntries() > 0))
	return;
      
      Int_t idigit = 0;
      triggers->Reset();
      
      while ((triggers->Next())) {
	Float_t L0Amplitude = 0;
	triggers->GetAmplitude(L0Amplitude);
	
	if (L0Amplitude <= 0 && fInputCellType != kL1FastORs)
	  continue;

	Int_t L1Amplitude = 0;
	triggers->GetL1TimeSum(L1Amplitude);
	
	if (L1Amplitude <= 0 && fInputCellType == kL1FastORs)
	  continue;
      
	Int_t triggerTime = 0;
	Int_t ntimes = 0;
	triggers->GetNL0Times(ntimes);
	
	if (ntimes < 1 && fInputCellType == kL0FastORsTC) 
	  continue;
	
	if (ntimes > 0) {
	  Int_t trgtimes[25];
	  triggers->GetL0Times(trgtimes);
	  triggerTime = trgtimes[0];
	}

	Int_t triggerCol = 0, triggerRow = 0;
	triggers->GetPosition(triggerCol, triggerRow);
	
	Int_t find = -1;
	geom->GetAbsFastORIndexFromPositionInEMCAL(triggerCol, triggerRow, find);
	
	if (find < 0)
	  continue;
      
	Int_t cidx[4] = {-1};
	Bool_t ret = geom->GetCellIndexFromFastORIndex(find, cidx);
	
	if (!ret)
	  continue;

	Float_t triggerAmplitude = 0;
	
	if (fInputCellType == kL1FastORs) {
	  triggerAmplitude = 0.25 * L1Amplitude;  // it will add 4 cells for 1 amplitude
	}
	else {
	  triggerAmplitude = L0Amplitude;      // 10 bit truncated, so it is already divided by 4
	}
	
	for (Int_t idxpos = 0; idxpos < 4; idxpos++) {
	  Int_t triggerNumber = cidx[idxpos];
	  AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->New(idigit));
	  digit->SetId(triggerNumber);
	  digit->SetTime(triggerTime);
	  digit->SetTimeR(triggerTime);
	  digit->SetIndexInList(idigit);
	  digit->SetType(AliEMCALDigit::kHG);
	  digit->SetAmplitude(triggerAmplitude);
	  idigit++;
	} 
      }
    }
    break;
  }
}

//________________________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::RecPoints2Clusters(TClonesArray *clus)
{
  // Cluster energy, global position, cells and their amplitude fractions are restored.

  AliVCaloCells *cells = InputEvent()->GetEMCALCells();
  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
  AliVEvent *event = InputEvent();
  if (!event) {
    AliError("Cannot get the event");
    return;
  } 

  // tracks array for track/cluster matching
  TClonesArray *tarr = 0;
  if (!fTrackName.IsNull()) {
    tarr = dynamic_cast<TClonesArray*>(event->FindListObject(fTrackName));
    if (!tarr) {
      AliError(Form("Cannot get tracks named %s", fTrackName.Data()));
    }
  }
  
  const Int_t Ncls = fClusterArr->GetEntries();
  AliDebug(1, Form("total no of clusters %d", Ncls)); 
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
    
    AliDebug(1, Form("energy %f", recpoint->GetEnergy()));
    
    AliVCluster *c = static_cast<AliVCluster*>(clus->New(nout++));
    c->SetType(AliVCluster::kEMCALClusterv1);
    c->SetE(recpoint->GetEnergy());
    c->SetPosition(g);
    c->SetNCells(ncells_true);
    c->SetCellsAbsId(absIds);
    c->SetCellsAmplitudeFraction(ratios);
    c->SetID(recpoint->GetUniqueID());
    c->SetDispersion(recpoint->GetDispersion());
    c->SetEmcCpvDistance(-1);
    c->SetChi2(-1);
    c->SetTOF(recpoint->GetTime()) ;     //time-of-flight
    c->SetNExMax(recpoint->GetNExMax()); //number of local maxima
    Float_t elipAxis[2];
    recpoint->GetElipsAxis(elipAxis);
    c->SetM02(elipAxis[0]*elipAxis[0]);
    c->SetM20(elipAxis[1]*elipAxis[1]);
    if (fPedestalData) {
      c->SetDistanceToBadChannel(recpoint->GetDistanceToBadTower()); 
    } else {
      if (fRecoUtils && fRecoUtils->IsBadChannelsRemovalSwitchedOn()) {
        fRecoUtils->RecalculateClusterDistanceToBadChannel(geom, cells, c);
      } 
    }

    if (tarr) {
      Double_t dEtaMin  = 1e9;
      Double_t dPhiMin  = 1e9;
      Int_t    imin     = -1;
      Double_t ceta     = gpos.Eta();
      Double_t cphi     = gpos.Phi();
      const Int_t ntrks = tarr->GetEntries();
      for(Int_t t = 0; t<ntrks; ++t) {
        AliVTrack *track = static_cast<AliVTrack*>(tarr->At(t));
        if (!track)
          continue;
        const AliExternalTrackParam *outp = track->GetOuterParam();
        if (!outp)
          continue;
        Double_t trkPos[3] = {0.,0.,0.};
        if (!outp->GetXYZ(trkPos)) 
          continue;
        TVector3 vec(trkPos[0],trkPos[1],trkPos[2]);
        Double_t veta = vec.Eta();
        Double_t vphi = vec.Phi();
        if(vphi<0)
          vphi += 2*TMath::Pi();
        if (TMath::Abs(veta)>0.75 || (vphi<70*TMath::DegToRad()) || (vphi>190*TMath::DegToRad()))
          continue;
        Double_t dR = vec.DeltaR(gpos);
        if(dR > 25) 
          continue;
        Float_t tmpEta=0, tmpPhi=0;
        if (0) {
          AliExternalTrackParam trkParTemp(*outp); // retrieve the starting point every time before the extrapolation
          Bool_t ret = fRecoUtils->ExtrapolateTrackToCluster(&trkParTemp, c, fRecoUtils->GetMass(), fRecoUtils->GetStep(), tmpEta, tmpPhi); 
          if (!ret)
            continue;
        } else {
          tmpEta = ceta - veta;
          tmpPhi = cphi - vphi;
        }
        if (TMath::Abs(tmpEta)<TMath::Abs(dEtaMin) && TMath::Abs(tmpPhi)<TMath::Abs(dPhiMin)) {
          dEtaMin = tmpEta;
          dPhiMin = tmpPhi;
          imin = t;
        }
      }
      c->SetEmcCpvDistance(imin);
      c->SetTrackDistance(dPhiMin, dEtaMin);
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
  const Int_t   ncells = cells->GetNumberOfCells();
  const Int_t   ndigis = fDigitsArr->GetEntries();
  if (ncells!=ndigis) {
    cells->DeleteContainer();
    cells->CreateContainer(ndigis);
  }
  for (Int_t idigit = 0; idigit < ndigis; ++idigit) {
    AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->At(idigit));
    Double_t cellAmplitude = digit->GetCalibAmp();
    Short_t cellNumber     = digit->GetId();
    Double_t cellTime      = digit->GetTime();
    cells->SetCell(idigit, cellNumber, cellAmplitude, cellTime);
  }
}

//________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::UpdateClusters()
{
  // Update cells in case re-calibration was done.
  
  if (!fAttachClusters)
    return;
  
  TClonesArray *clus = 0;
  
  if (fOverwrite) {
    clus = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("caloClusters"));
    if (!clus)
      clus = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("CaloClusters"));
    if (!clus)
      return;
    
    const Int_t nents = clus->GetEntries();
    for (Int_t i=0;i<nents;++i) {
      AliVCluster *c = static_cast<AliVCluster*>(clus->At(i));
      if (!c)
        continue;
      if (c->IsEMCAL())
        delete clus->RemoveAt(i);
    }
  } else {
    clus = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fNewClusterArrayName));
    if (!clus) {
      clus = new TClonesArray("AliESDCaloCluster");
      clus->SetName(fNewClusterArrayName);
      InputEvent()->AddObject(clus);
    } else {
      clus->Delete();
    }
  }
  
  RecPoints2Clusters(clus);
}

//________________________________________________________________________________________
void AliAnalysisTaskEMCALClusterizeFast::Init()
{
  // Select clusterization/unfolding algorithm and set all the needed parameters.

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
    fUnfolder = new AliEMCALAfterBurnerUF(fRecParam->GetW0(),fRecParam->GetLocMaxCut(),fRecParam->GetMinECut());
    return;
  }

  AliEMCALGeometry *geometry = 0;
  if (fGeomName.Length()>0) 
    geometry = AliEMCALGeometry::GetInstance(fGeomName);
  else
    geometry = AliEMCALGeometry::GetInstance();
  if (!geometry) {
    AliFatal("Geometry not available!!!");
    return;
  }

  if (!fGeomMatrixSet) {
    if (fLoadGeomMatrices) {
      for(Int_t mod=0; mod < geometry->GetNumberOfSuperModules(); ++mod) {
        if (fGeomMatrix[mod]){
          if (DebugLevel() > 2) 
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
          if (DebugLevel() > 2) 
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
  else if (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerNxN) {
   AliEMCALClusterizerNxN *clusterizer = new AliEMCALClusterizerNxN(geometry);
   clusterizer->SetNRowDiff(fRecParam->GetNRowDiff());
   clusterizer->SetNColDiff(fRecParam->GetNColDiff());
    fClusterizer = clusterizer;
  } else if (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv2) 
    fClusterizer = new AliEMCALClusterizerv2(geometry);
  else if (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerFW){
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
