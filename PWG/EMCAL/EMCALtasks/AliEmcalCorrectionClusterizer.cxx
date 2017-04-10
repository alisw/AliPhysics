// AliEmcalCorrectionClusterizer
//
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

// --- Root ---
#include <TObjArray.h>
#include <TArrayI.h>
#include <TStopwatch.h>

// --- AliRoot ---
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
#include "AliInputEventHandler.h"
#include "AliClusterContainer.h"
#include "AliAODMCParticle.h"
#include "AliEMCALRecoUtils.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisManager.h"

#include "AliEmcalCorrectionClusterizer.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionClusterizer);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionClusterizer> AliEmcalCorrectionClusterizer::reg("AliEmcalCorrectionClusterizer");

const std::map <std::string, AliEmcalCorrectionClusterizer::EmbeddedCellEnergyType> AliEmcalCorrectionClusterizer::fgkEmbeddedCellEnergyTypeMap = {
  {"kNonEmbedded", EmbeddedCellEnergyType::kNonEmbedded },
  {"kEmbeddedDataMCOnly", EmbeddedCellEnergyType::kEmbeddedDataMCOnly },
  {"kEmbeddedDataExcludeMC", EmbeddedCellEnergyType::kEmbeddedDataExcludeMC }
};

const std::map <std::string, AliEMCALRecParam::AliEMCALClusterizerFlag> AliEmcalCorrectionClusterizer::fgkClusterizerTypeMap = {
  {"kClusterizerv1", AliEMCALRecParam::kClusterizerv1 },
  {"kClusterizerNxN", AliEMCALRecParam::kClusterizerNxN },
  {"kClusterizerv2", AliEMCALRecParam::kClusterizerv2 },
  {"kClusterizerFW", AliEMCALRecParam::kClusterizerFW }
};

/**
 * Default constructor
 */
AliEmcalCorrectionClusterizer::AliEmcalCorrectionClusterizer() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionClusterizer"),
  fDigitsArr(0),
  fClusterArr(0),
  fRecParam(new AliEMCALRecParam),
  fClusterizer(0),
  fUnfolder(0),
  fJustUnfold(kFALSE),
  fGeomName(),
  fGeomMatrixSet(kFALSE),
  fLoadGeomMatrices(kFALSE),
  fOCDBpath(),
  fCalibData(0),
  fPedestalData(0),
  fLoadCalib(kFALSE),
  fLoadPed(kFALSE),
  fSubBackground(kFALSE),
  fNPhi(4),
  fNEta(4),
  fShiftPhi(2),
  fShiftEta(2),
  fTRUShift(0),
  fEmbeddedCellEnergyType(kNonEmbedded),
  fTestPatternInput(kFALSE),
  fSetCellMCLabelFromCluster(0),
  fSetCellMCLabelFromEdepFrac(0),
  fRemapMCLabelForAODs(0),
  fCaloClusters(0),
  fEsd(0),
  fAod(0),
  fRecalDistToBadChannels(kFALSE),
  fRecalShowerShape(kFALSE)
{
  for(Int_t i = 0; i < AliEMCALGeoParams::fgkEMCALModules; i++) fGeomMatrix[i] = 0 ;
  for(Int_t j = 0; j < fgkTotalCellNumber;                 j++)
  { fOrgClusterCellId[j] =-1; fCellLabels[j] =-1 ; }
}

/**
 * Destructor
 */
AliEmcalCorrectionClusterizer::~AliEmcalCorrectionClusterizer()
{
  delete fClusterizer;
  delete fUnfolder;
  delete fRecParam;
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionClusterizer::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();

  std::string clusterizerTypeStr = "";
  GetProperty("clusterizer", clusterizerTypeStr);
  UInt_t clusterizerType = fgkClusterizerTypeMap.at(clusterizerTypeStr);
  Double_t cellE  = 0.05;
  GetProperty("cellE", cellE);
  Double_t seedE  = 0.1;
  GetProperty("seedE", seedE);
  Float_t timeMin = -1;             //minimum time of physical signal in a cell/digit (s) (in run macro, -50e-6)
  GetProperty("cellTimeMin", timeMin);
  Float_t timeMax = +1;             //maximum time of physical signal in a cell/digit (s) (in run macro, 50e-6)
  GetProperty("cellTimeMax", timeMax);
  Float_t timeCut =  1;             //maximum time difference between the digits inside EMC cluster (s) (in run macro, 1e-6)
  GetProperty("clusterTimeLength", timeCut);
  Float_t w0 = 4.5;
  GetProperty("w0", w0);
  GetProperty("recalDistToBadChannels", fRecalDistToBadChannels);
  GetProperty("recalShowerShape", fRecalShowerShape);
  GetProperty("remapMcAod", fRemapMCLabelForAODs);
  Bool_t enableFracEMCRecalc = kFALSE;
  GetProperty("enableFracEMCRecalc", enableFracEMCRecalc);
  GetProperty("setCellMCLabelFromCluster", fSetCellMCLabelFromCluster);
  Float_t diffEAggregation = 0.;
  GetProperty("diffEAggregation", diffEAggregation);
  GetProperty("useTestPatternForInput", fTestPatternInput);
  
  Int_t removeNMCGenerators = 0;
  GetProperty("removeNMCGenerators", removeNMCGenerators);
  Bool_t enableMCGenRemovTrack = 1;
  GetProperty("enableMCGenRemovTrack", enableMCGenRemovTrack);
  std::string removeMcGen1 = "";
  GetProperty("removeMCGen1", removeMcGen1);
  TString removeMCGen1 = removeMcGen1.c_str();
  std::string removeMcGen2 = "";
  GetProperty("removeMCGen2", removeMcGen2);
  TString removeMCGen2 = removeMcGen2.c_str();

  fRecParam->SetClusterizerFlag(clusterizerType);
  fRecParam->SetMinECut(cellE);
  fRecParam->SetClusteringThreshold(seedE);
  fRecParam->SetW0(w0);
  fRecParam->SetTimeMin(timeMin);
  fRecParam->SetTimeMax(timeMax);
  fRecParam->SetTimeCut(timeCut);
  fRecParam->SetLocMaxCut(diffEAggregation);      // Set minimum energy difference to start new cluster
  
  if (clusterizerType == AliEMCALRecParam::kClusterizerNxN)
    fRecParam->SetNxM(1,1); // -> (1,1) means 3x3!
  
  // init reco utils
  if (!fRecoUtils)
    fRecoUtils = new AliEMCALRecoUtils;
  
  if(enableFracEMCRecalc){
    fSetCellMCLabelFromEdepFrac = kTRUE;
    fSetCellMCLabelFromCluster  = 0;
    
    printf("Enable frac MC recalc, remove generators %d \n",removeNMCGenerators);
    if(removeNMCGenerators > 0)
    {
      printf("\t gen1 <%s>, gen2 <%s>, remove tracks %d\n",removeMCGen1.Data(),removeMCGen2.Data(),enableMCGenRemovTrack);
      fRecoUtils->SetNumberOfMCGeneratorsToAccept(removeNMCGenerators) ;
      fRecoUtils->SetNameOfMCGeneratorsToAccept(0,removeMCGen1);
      fRecoUtils->SetNameOfMCGeneratorsToAccept(1,removeMCGen2);
      
      if(!enableMCGenRemovTrack) fRecoUtils->SwitchOffMCGeneratorToAcceptForTrackMatching();
    }
  }
  
  std::string embeddedCellEnergyTypeStr = "";
  GetProperty("embeddedCellEnergyType", embeddedCellEnergyTypeStr);
  fEmbeddedCellEnergyType = fgkEmbeddedCellEnergyTypeMap.at(embeddedCellEnergyTypeStr);
  //Printf("embeddedCellEnergyType: %d",fEmbeddedCellEnergyType);

  // Only support one cluster container for the clusterizer!
  if (fClusterCollArray.GetEntries() > 1) {
    AliFatal("Passed more than one cluster container to the clusterizer, but the clusterizer only supports one cluster container!");
  }
  
  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionClusterizer::UserCreateOutputObjects()
{   
  AliEmcalCorrectionComponent::UserCreateOutputObjects();

  if (fCreateHisto){
    fHistCPUTime = new TH1F("hCPUTime","hCPUTime;CPU Time (ms)", 2000, 0, 1000);
    fOutput->Add(fHistCPUTime);

    fHistRealTime = new TH1F("hRealTime","hRealTime;Real Time (ms)", 2000, 0, 1000);
    fOutput->Add(fHistRealTime);

    fTimer = new TStopwatch();
  }
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionClusterizer::Run()
{
  // Time the event loop if histogram creation is enabled
  if (fCreateHisto) {
    fTimer->Start(kTRUE);
  }

  AliEmcalCorrectionComponent::Run();
  
  fEsd = dynamic_cast<AliESDEvent*>(fEvent);
  fAod = dynamic_cast<AliAODEvent*>(fEvent);

  // Only support one cluster container in the clusterizer!
  AliClusterContainer * clusCont = GetClusterContainer(0);
  if (!clusCont) {
    AliFatal("Could not retrieve cluster container!");
  }
  
  fCaloClusters = clusCont->GetArray();

  // If cells are empty, clear clusters and return
  if (fCaloCells->GetNumberOfCells()<=0)
  {
    AliWarning(Form("Number of EMCAL cells = %d, returning", fCaloCells->GetNumberOfCells()));
    ClearEMCalClusters();
    return kFALSE;
  }
  
  UInt_t offtrigger = 0;
  if (fEsd) {
    UInt_t mask1 = fEsd->GetESDRun()->GetDetectorsInDAQ();
    UInt_t mask2 = fEsd->GetESDRun()->GetDetectorsInReco();
    Bool_t desc1 = (mask1 >> 18) & 0x1;
    Bool_t desc2 = (mask2 >> 18) & 0x1;
    if (desc1==0 || desc2==0) { //AliDAQ::OfflineModuleName(180=="EMCAL"
      AliError(Form("EMCAL not in DAQ/RECO: %u (%u)/%u (%u)",
                    mask1, fEsd->GetESDRun()->GetDetectorsInReco(),
                    mask2, fEsd->GetESDRun()->GetDetectorsInDAQ()));
      return kFALSE;
    }
    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    offtrigger = ((AliInputEventHandler*)(am->GetInputEventHandler()))->IsEventSelected();
  } else {
    offtrigger =  ((AliVAODHeader*)fAod->GetHeader())->GetOfflineTrigger();
  }
  
  if (!fMCEvent) {
    if (offtrigger & AliVEvent::kFastOnly) {
      AliError(Form("EMCAL not in fast only partition"));
      return kFALSE;
    }
  }
  
  Init();
  
  if (fJustUnfold) {
    AliWarning("Unfolding not implemented");
    return kTRUE;
  }
  
  FillDigitsArray();
  
  Clusterize();
  
  UpdateClusters();
  
  CalibrateClusters();

  if (fCreateHisto) {
    fTimer->Stop();
    // Ensure that it is stored in ms
    fHistCPUTime->Fill(fTimer->CpuTime() * 1000.);
    fHistRealTime->Fill(fTimer->RealTime() * 1000.);
  }
  
  return kTRUE;
}

/**
 * Clusterize digits into clusters.
 */
void AliEmcalCorrectionClusterizer::Clusterize()
{
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

/**
 * Fill digits array from input cell collection.
 */
void AliEmcalCorrectionClusterizer::FillDigitsArray()
{
  fDigitsArr->Clear("C");
  if (fTestPatternInput) {
    // Fill digits from a pattern
    Int_t maxd = fGeom->GetNCells() / 4;
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
  else
  {
    // In case of MC productions done before aliroot tag v5-02-Rev09
    // passing the cluster label to all the cells belonging to this cluster
    // very rough
    // Copied and simplified from AliEMCALTenderSupply
    if (fSetCellMCLabelFromCluster || fSetCellMCLabelFromEdepFrac)
    {
      for (Int_t i = 0; i < fgkTotalCellNumber; i++)
      {
        fCellLabels      [i] =-1 ;
        fOrgClusterCellId[i] =-1 ;
      }
      
      Int_t nClusters = fEvent->GetNumberOfCaloClusters();
      for (Int_t i = 0; i < nClusters; i++)
      {
        AliVCluster *clus =  fEvent->GetCaloCluster(i);
        
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
      } // cluster loop
    }

    Double_t avgE        = 0; // for background subtraction
    const Int_t ncells   = fCaloCells->GetNumberOfCells();
    for (Int_t icell = 0, idigit = 0; icell < ncells; ++icell)
    {
      Float_t cellAmplitude=0;
      Double_t cellTime=0, amp = 0, cellEFrac = 0;
      Short_t  cellNumber=0;
      Int_t cellMCLabel=-1;
      if (fCaloCells->GetCell(icell, cellNumber, amp, cellTime, cellMCLabel, cellEFrac) != kTRUE)
        break;

      cellAmplitude = amp; // compilation problem

      //if (fSetCellMCLabelFromCluster) cellMCLabel = fCellLabels[cellNumber];
      if(!fSetCellMCLabelFromEdepFrac)
      {
        if      (fSetCellMCLabelFromCluster) cellMCLabel = fCellLabels[cellNumber];
        else if (fRemapMCLabelForAODs      ) RemapMCLabelForAODs(cellMCLabel);
      }

      if (cellMCLabel > 0 && cellEFrac < 1e-6)
        cellEFrac = 1;

      if (cellAmplitude < 1e-6 || cellNumber < 0)
        continue;

      if (fEmbeddedCellEnergyType == kEmbeddedDataMCOnly) {
        if (cellMCLabel <= 0)
          continue;
        else {
          cellAmplitude *= cellEFrac;
          cellEFrac = 1;
        }
      }
      else if (fEmbeddedCellEnergyType == kEmbeddedDataExcludeMC) {
        if (cellMCLabel > 0)
          continue;
        else {
          cellAmplitude *= 1 - cellEFrac;
          cellEFrac = 0;
        }
      }

      // New way to set the cell MC labels,
      // valid only for MC productions with aliroot > v5-07-21
      //
      TArrayI labeArr(0);
      TArrayF eDepArr(0);
      Int_t nLabels = 0;

      if(fSetCellMCLabelFromEdepFrac && fOrgClusterCellId[cellNumber] >= 0) // index can be negative if noisy cell that did not form cluster
      {
        cellMCLabel = -1;
        
        fCellLabels[cellNumber] = idigit;

        Int_t iclus = fOrgClusterCellId[cellNumber];

        if(iclus < 0)
        {
          AliInfo("Negative original cluster index, skip \n");
          continue;
        }

        AliVCluster* clus = fEvent->GetCaloCluster(iclus);

        for(Int_t icluscell=0; icluscell < clus->GetNCells(); icluscell++ )
        {
          if(cellNumber != clus->GetCellAbsId(icluscell)) continue ;

          // Select the MC label contributing, only if enough energy deposition
          fRecoUtils->RecalculateCellLabelsRemoveAddedGenerator(cellNumber, clus, fMCEvent, cellAmplitude, labeArr, eDepArr);
          nLabels = labeArr.GetSize();
        }
      }

      AliEMCALDigit *digit = new((*fDigitsArr)[idigit]) AliEMCALDigit(cellMCLabel, cellMCLabel, cellNumber,
                                                                  cellAmplitude, (Float_t)cellTime,
                                                                  AliEMCALDigit::kHG,idigit, 0, 0, cellEFrac*cellAmplitude);

      if(nLabels > 0)
      {
        digit->SetListOfParents(nLabels,labeArr.GetArray(),eDepArr.GetArray());
      }

      if (fSubBackground)
      {
        Float_t energy = cellAmplitude;
        Float_t time   = cellTime;
        fClusterizer->Calibrate(energy,time,cellNumber);
        digit->SetAmplitude(energy);
        avgE += energy;
      }

      idigit++;
    }

    if (fSubBackground) {
      avgE /= fGeom->GetNumberOfSuperModules()*48*24;
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

/**
 * Convert AliEMCALRecoPoints to AliESDCaloClusters/AliAODCaloClusters.
 * Cluster energy, global position, cells and their amplitude fractions are restored.
 */
void AliEmcalCorrectionClusterizer::RecPoints2Clusters(TClonesArray *clus)
{
  const Int_t Ncls = fClusterArr->GetEntries();
  AliDebug(1, Form("total no of clusters %d", Ncls));
  
  for(Int_t i=0, nout=clus->GetEntries(); i < Ncls; ++i)
  {
    AliEMCALRecPoint *recpoint = static_cast<AliEMCALRecPoint*>(fClusterArr->At(i));
    
    Int_t ncellsTrue = 0;
    const Int_t ncells = recpoint->GetMultiplicity();
    UShort_t   absIds[ncells];
    Double32_t ratios[ncells];
    Int_t   *dlist = recpoint->GetDigitsList();
    Float_t *elist = recpoint->GetEnergiesList();
    Double_t mcEnergy = 0;
    
    for (Int_t c = 0; c < ncells; ++c)
    {
      AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->At(dlist[c]));
      absIds[ncellsTrue] = digit->GetId();
      ratios[ncellsTrue] = elist[c]/digit->GetAmplitude();
      
      if (digit->GetIparent(1) > 0)
        mcEnergy += digit->GetDEParent(1)/recpoint->GetEnergy();
      
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
    
    AliDebug(1, Form("energy %f", recpoint->GetEnergy()));
    
    AliVCluster *c = static_cast<AliVCluster*>(clus->New(nout++));
    c->SetType(AliVCluster::kEMCALClusterv1);
    c->SetE(recpoint->GetEnergy());
    c->SetPosition(g);
    c->SetNCells(ncellsTrue);
    c->SetCellsAbsId(absIds);
    c->SetCellsAmplitudeFraction(ratios);
    c->SetID(nout-1);
    c->SetDispersion(recpoint->GetDispersion());
    c->SetEmcCpvDistance(-1);
    c->SetChi2(-1);
    c->SetTOF(recpoint->GetTime()) ;     //time-of-flight
    c->SetNExMax(recpoint->GetNExMax()); //number of local maxima
    Float_t elipAxis[2];
    recpoint->GetElipsAxis(elipAxis);
    c->SetM02(elipAxis[0]*elipAxis[0]);
    c->SetM20(elipAxis[1]*elipAxis[1]);
    c->SetMCEnergyFraction(mcEnergy);
    
    //
    // MC labels
    //
    Int_t    parentMult   = 0;
    Int_t   *parentList   = recpoint->GetParents(parentMult);
    Float_t *parentListDE = recpoint->GetParentsDE();  // deposited energy
    
    c->SetLabel(parentList, parentMult);
    if(fSetCellMCLabelFromEdepFrac) {
      c->SetClusterMCEdepFractionFromEdepArray(parentListDE);
    }
    
    //
    // Set the cell energy deposition fraction map:
    //
    if( parentMult > 0 && fSetCellMCLabelFromEdepFrac )
    {
      UInt_t * mcEdepFracPerCell = new UInt_t[ncellsTrue];
      
      // Get the digit that originated this cell cluster
      //AliVCaloCells* cells = InputEvent()->GetEMCALCells();
      
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

/**
 * Clear the old clusters and fill the new clusters.
 */
void AliEmcalCorrectionClusterizer::UpdateClusters()
{
  // Before destroying the orignal list, assign to the rec points the MC labels
  // of the original clusters, if requested
  if (fSetCellMCLabelFromCluster == 2)
    SetClustersMCLabelFromOriginalClusters() ;
  
  ClearEMCalClusters();
  
  fCaloClusters->Compress();
  
  RecPoints2Clusters(fCaloClusters);
}

/**
 * Go through clusters one by one and process separate correction.
 */
void AliEmcalCorrectionClusterizer::CalibrateClusters()
{
  Int_t nclusters = fCaloClusters->GetEntriesFast();
  for (Int_t icluster=0; icluster < nclusters; ++icluster) {
    AliVCluster *clust = static_cast<AliVCluster*>(fCaloClusters->At(icluster));
    if (!clust)
      continue;
    
    // SHOWER SHAPE -----------------------------------------------
    if (fRecalShowerShape)
      fRecoUtils->RecalculateClusterShowerShapeParameters(fGeom, fCaloCells, clust);
    
    // DISTANCE TO BAD CHANNELS -----------------------------------
    if (fRecalDistToBadChannels)
      fRecoUtils->RecalculateClusterDistanceToBadChannel(fGeom, fCaloCells, clust);
    
  }
  
  fCaloClusters->Compress();
}

/**
 * MC label for Cells not remapped after ESD filtering -- do it here.
 */
void AliEmcalCorrectionClusterizer::RemapMCLabelForAODs(Int_t & label)
{
  if (label < 0) return;
  
  TClonesArray * arr = dynamic_cast<TClonesArray*>(fAod->FindListObject("mcparticles")) ;
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

/**
 * Get the original clusters that contribute to the new rec point cluster,
 * assign the labels of such clusters to the new rec point cluster.
 * Only approximatedly valid  when output are V1 clusters, or input/output clusterizer
 * are the same handle with care
 * Copy from same method in AliAnalysisTaskEMCALClusterize, but here modify the recpoint and
 * not the output calocluster
 */
void AliEmcalCorrectionClusterizer::SetClustersMCLabelFromOriginalClusters()
{
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
          nLabTotOrg+=(fEvent->GetCaloCluster(idCluster))->GetNLabels();
          
          //Search highest E cluster
          AliVCluster * clOrg = fEvent->GetCaloCluster(idCluster);
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
      AliVCluster * clOrg = fEvent->GetCaloCluster(idCluster);
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

/**
 * Initialize the clusterizer.
 */
void AliEmcalCorrectionClusterizer::Init()
{
  // Select clusterization/unfolding algorithm and set all the needed parameters.
  
  // If distBC option enabled, init and fill bad channel map
  if (fRecalDistToBadChannels)
    fRecoUtils->SwitchOnDistToBadChannelRecalculation();
  else
    fRecoUtils->SwitchOffDistToBadChannelRecalculation();
  
  CheckIfRunChanged();
  
  if (fJustUnfold){
    // init the unfolding afterburner
    delete fUnfolder;
    fUnfolder = new AliEMCALAfterBurnerUF(fRecParam->GetW0(),fRecParam->GetLocMaxCut(),fRecParam->GetMinECut());
    return;
  }
  
  if (!fGeomMatrixSet) {
    if (fLoadGeomMatrices) {
      for(Int_t mod=0; mod < fGeom->GetNumberOfSuperModules(); ++mod) {
        if (fGeomMatrix[mod]){
          if (AliAnalysisManager::GetAnalysisManager()->GetDebugLevel() > 2)
            fGeomMatrix[mod]->Print();
        fGeom->SetMisalMatrix(fGeomMatrix[mod],mod);
        }
      }
    } else { // get matrix from file (work around bug in aliroot)
      for(Int_t mod=0; mod < fGeom->GetEMCGeometry()->GetNumberOfSuperModules(); ++mod) {
        const TGeoHMatrix *gm = 0;
        if (fEsd) {
          gm = fEsd->GetEMCALMatrix(mod);
        } else {
          AliAODHeader *aodheader = dynamic_cast<AliAODHeader*>(fAod->GetHeader());
          if(!aodheader) AliFatal("Not a standard AOD");
          if (aodheader) {
            gm = aodheader->GetEMCALMatrix(mod);
          }
        }
        if (gm) {
          if (AliAnalysisManager::GetAnalysisManager()->GetDebugLevel() > 2)
            gm->Print();
        fGeom->SetMisalMatrix(gm,mod);
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
  if (fClusterizer) {
    // avoid to delete digits array
    fClusterizer->SetDigitsArr(0);
    delete fClusterizer;
  }
  if (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv1)
    fClusterizer = new AliEMCALClusterizerv1(fGeom);
  else if (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerNxN) {
    AliEMCALClusterizerNxN *clusterizer = new AliEMCALClusterizerNxN(fGeom);
    clusterizer->SetNRowDiff(fRecParam->GetNRowDiff()); //MV: already done in AliEMCALClusterizer::InitParameters
    clusterizer->SetNColDiff(fRecParam->GetNColDiff()); //MV: already done in AliEMCALClusterizer::InitParameters
    fClusterizer = clusterizer;
  }
  else if (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv2)
    fClusterizer = new AliEMCALClusterizerv2(fGeom);
  else if (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerFW) {
    AliEMCALClusterizerFixedWindow *clusterizer = new AliEMCALClusterizerFixedWindow(fGeom);
    clusterizer->SetNphi(fNPhi);
    clusterizer->SetNeta(fNEta);
    clusterizer->SetShiftPhi(fShiftPhi);
    clusterizer->SetShiftEta(fShiftEta);
    clusterizer->SetTRUshift(fTRUShift);
    fClusterizer = clusterizer;
  }
  else {
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

/**
 * This function is called if the run changes (it inherits from the base component),
 * to load a new time calibration and fill relevant variables.
 */
Bool_t AliEmcalCorrectionClusterizer::CheckIfRunChanged()
{
  Bool_t runChanged = AliEmcalCorrectionComponent::CheckIfRunChanged();
  
  if (runChanged && fRecalDistToBadChannels) {
    // init bad channels
    Int_t fInitBC = InitBadChannels();
    if (fInitBC==0) {
      AliError("InitBadChannels returned false, returning");
    }
    if (fInitBC==1) {
      AliWarning("InitBadChannels OK");
    }
    if (fInitBC>1) {
      AliWarning(Form("No external hot channel set: %d - %s", fEvent->GetRunNumber(), fFilepass.Data()));
    }
  }
  return runChanged;
}

/**
 * Clear the EMCal clusters from the cluster TClonesArray.
 */
void AliEmcalCorrectionClusterizer::ClearEMCalClusters()
{
  const Int_t nents = fCaloClusters->GetEntries();
  for (Int_t i=0;i<nents;++i) {
    AliVCluster *c = static_cast<AliVCluster*>(fCaloClusters->At(i));
    if (!c) {
      continue;
    }
    if (c->IsEMCAL()) {
      delete fCaloClusters->RemoveAt(i);
    }
  }
}
