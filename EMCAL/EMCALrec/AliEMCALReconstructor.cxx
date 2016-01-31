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

// --- ROOT system ---
#include <TClonesArray.h>
#include "TGeoManager.h"
#include "TGeoMatrix.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliEMCALReconstructor.h"

#include "AliCodeTimer.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALCalibData.h"
#include "AliEMCALCalibTime.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliESDtrack.h"
#include "AliEMCALLoader.h"
#include "AliEMCALRawUtils.h"
#include "AliEMCALDigit.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALClusterizerv2.h"
#include "AliEMCALClusterizerNxN.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALPID.h"
#include "AliEMCALRecoUtils.h"
#include "AliRawReader.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliEMCALGeometry.h"
#include "AliEMCAL.h"
#include "AliESDVZERO.h"
#include "AliCDBManager.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliEMCALTriggerData.h"
#include "AliEMCALTriggerElectronics.h"
#include "AliEMCALTriggerDCSConfigDB.h"
#include "AliEMCALTriggerDCSConfig.h"
#include "AliEMCALTriggerData.h"
#include "AliEMCALTriggerRawDigit.h"
#include "AliEMCALTriggerPatch.h"
#include "AliEMCALTriggerTypes.h"

/// \cond CLASSIMP
ClassImp(AliEMCALReconstructor) 
/// \endcond
  
const AliEMCALRecParam*     AliEMCALReconstructor::fgkRecParam        = 0;   // EMCAL rec. parameters
AliEMCALRawUtils*           AliEMCALReconstructor::fgRawUtils         = 0;   // EMCAL raw utilities class
AliEMCALClusterizer*        AliEMCALReconstructor::fgClusterizer      = 0;   // EMCAL clusterizer class
TClonesArray*               AliEMCALReconstructor::fgDigitsArr        = 0;   // list of digits, to be used multiple times
TObjArray*                  AliEMCALReconstructor::fgClustersArr      = 0;   // list of clusters, to be used multiple times
TClonesArray*               AliEMCALReconstructor::fgTriggerDigits    = 0;   // list of trigger digits, to be used multiple times
AliEMCALTriggerElectronics *AliEMCALReconstructor::fgTriggerProcessor = 0x0;
TClonesArray               *AliEMCALReconstructor::fgTriggerData      = 0x0;

///
/// Constructor
/// Init all arrays and stuff.
///
//____________________________________________________________________________
AliEMCALReconstructor::AliEMCALReconstructor() 
  : fGeom(0),fCalibData(0),fCalibTime(0),fPedestalData(0), fMatches(0x0)
{
  fgRawUtils = new AliEMCALRawUtils;
  
  AliCDBManager* man = AliCDBManager::Instance();
  
  //----------------------------------
  // Get the geometry, 3 posibilities:
  //  * To make sure we match with the geometry in a simulation file,
  //    let's try to get it first.  
  //  * If not, check the run number assigned for this chunk and set the 
  //    geometry depending on the run number
  //  * If not, take the default geometry
  
  AliRunLoader *rl = AliRunLoader::Instance();
  if (rl->GetAliRun())
  {
    AliEMCAL * emcal = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"));
    if(emcal) fGeom = emcal->GetGeometry();
  }
  
  if(!fGeom)
  {
    Int_t runNumber = man->GetRun();
    fGeom =  AliEMCALGeometry::GetInstanceFromRunNumber(runNumber);
  }
  
  if(!fGeom) 
  {
    AliWarning(Form("Using default geometry in reconstruction!!!"));
    fGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
  }
  
  if ( !fGeom ) AliFatal(Form("Could not get geometry!"));
  else          AliInfo (Form("Geometry name: <<%s>>",fGeom->GetName())); 
  
  //---------------------------
  // Get energy calibration parameters	
  if(!fCalibData)
  {
      AliCDBEntry *entry = (AliCDBEntry*)  man->Get("EMCAL/Calib/Data");
      if (entry) fCalibData =  (AliEMCALCalibData*) entry->GetObject();
  }
  
  if(!fCalibData)
    AliFatal("Energy Calibration parameters not found in CDB!");
  
  
  //---------------------------
  // Get time calibration parameters if requested	
  if(!fCalibTime)
  {
    AliCDBEntry *entry = (AliCDBEntry*)  man->Get("EMCAL/Calib/Time");
    if (entry) fCalibTime =  (AliEMCALCalibTime*) entry->GetObject();
  }
  
  if(!fCalibTime)
    AliFatal("Time Calibration parameters not found in CDB!");
  
  //------------------
  // Get bad channels	
  if(!fPedestalData)
  {
      AliCDBEntry *entry = (AliCDBEntry*) man->Get("EMCAL/Calib/Pedestals");
      if (entry) fPedestalData =  (AliCaloCalibPedestal*) entry->GetObject();
  }
  
  if(!fPedestalData)
    AliFatal("Dead map not found in CDB!");
  
  //----------------------------------------------------
  // Get trigger parameters and init other trigger stuff
  AliEMCALTriggerDCSConfigDB* dcsConfigDB = AliEMCALTriggerDCSConfigDB::Instance();
  
  const AliEMCALTriggerDCSConfig* dcsConfig = dcsConfigDB->GetTriggerDCSConfig();
  
  if (!dcsConfig) 
    AliFatal("No Trigger DCS Configuration from OCDB!");
  
  fgTriggerProcessor = new AliEMCALTriggerElectronics( dcsConfig );
  
  int dsize = (fGeom->GetTriggerMappingVersion() == 2) ? 2 : 1;
  fgTriggerData = new TClonesArray("AliEMCALTriggerData",dsize);
  for (int i=0;i<dsize;i++) {
    new((*fgTriggerData)[i]) AliEMCALTriggerData();
  }

  //-----------------------------
  // Init temporary list of digits
  fgDigitsArr     = new TClonesArray("AliEMCALDigit",1000);
  fgClustersArr   = new TObjArray(1000);

  const int kNTRU = fGeom->GetNTotalTRU();
  fgTriggerDigits = new TClonesArray("AliEMCALTriggerRawDigit", kNTRU * 96);	
	
  //--------------------------
  // Init Track matching array
  fMatches = new TList();
  fMatches->SetOwner(kTRUE);
} 

///
/// Destructor
///
//____________________________________________________________________________
AliEMCALReconstructor::~AliEMCALReconstructor()
{
  ////RS  if(fGeom)              delete fGeom;
  
  //No need to delete, recovered from OCDB
  //if(fCalibData)         delete fCalibData;
  //if(fCalibTime)         delete fCalibTime;
  //if(fPedestalData)      delete fPedestalData;
  
  if(fgDigitsArr) fgDigitsArr->Clear("C");
  delete fgDigitsArr; 
  fgDigitsArr = 0;
  
  if(fgClustersArr) fgClustersArr->Clear();
  delete fgClustersArr; 
  fgClustersArr = 0;
  
  if(fgTriggerDigits) fgTriggerDigits->Clear();
  delete fgTriggerDigits; 
  fgTriggerDigits = 0;
  
  delete fgRawUtils;
  fgRawUtils = 0;
  delete fgClusterizer;
  fgClusterizer = 0;
  
  delete fgTriggerProcessor;
  fgTriggerProcessor = 0;
  
  if(fMatches) { fMatches->Delete(); delete fMatches; fMatches = 0;}
  
  AliCodeTimer::Instance()->Print();
} 

///
/// Init the fgClusterizer with geometry and calibration pointers, avoid doing it twice.                          
///
//____________________________________________________________________________                                  
void AliEMCALReconstructor::InitClusterizer() const
{
  Int_t clusterizerType = -1;
  Int_t eventType = -1;
  if(GetRecParam()) 
  {
    clusterizerType = GetRecParam()->GetClusterizerFlag();
    eventType       = GetRecParam()->GetEventSpecie();
  }
  else
  {
    AliCDBEntry *entry = (AliCDBEntry*)
      AliCDBManager::Instance()->Get("EMCAL/Calib/RecoParam");
  
    // Get The reco param for the default event specie                                                           
    if (entry) 
    {
      AliEMCALRecParam *recParam  = (AliEMCALRecParam*)((TObjArray *) entry->GetObject())->At(0);
      if(recParam) clusterizerType = recParam->GetClusterizerFlag(); 
    }
  }
  
  // Check if clusterizer previously set corresponds to what is needed for this event type                       
  if(fgClusterizer)
  {
    if(eventType!=AliRecoParam::kCalib)
    {
      //printf("ReCreate clusterizer? Clusterizer set <%d>, Clusterizer in use <%s>\n",
      //     clusterizerType, fgClusterizer->Version());
      
      if     (clusterizerType == AliEMCALRecParam::kClusterizerv1  && !strcmp(fgClusterizer->Version(),"clu-v1"))  return;
      
      else if(clusterizerType == AliEMCALRecParam::kClusterizerNxN && !strcmp(fgClusterizer->Version(),"clu-NxN")) return;
      
      else if(clusterizerType == AliEMCALRecParam::kClusterizerv2  && !strcmp(fgClusterizer->Version(),"clu-v2"))  return;
      
      //Need to create new clusterizer, the one set previously is not the correct one     
      delete fgClusterizer;
    }
    else return;
  }
    
  if      (clusterizerType  == AliEMCALRecParam::kClusterizerv1)
    {
      fgClusterizer = new AliEMCALClusterizerv1 (fGeom, fCalibData,fCalibTime,fPedestalData);
    }
  else if (clusterizerType  == AliEMCALRecParam::kClusterizerNxN)
    {
      fgClusterizer = new AliEMCALClusterizerNxN(fGeom, fCalibData,fCalibTime,fPedestalData);
    }
  else if (clusterizerType  == AliEMCALRecParam::kClusterizerv2)
  {
    fgClusterizer = new AliEMCALClusterizerv2   (fGeom, fCalibData,fCalibTime,fPedestalData);
  }
  else 
  {
    AliFatal(Form("Unknown clusterizer %d ", clusterizerType));
  }
}

///
/// Method called by AliReconstruction; 
/// Only the clusterization is performed, the rest of the reconstruction is done 
/// in FillESD because the track matching needs access to the AliESD object to 
/// retrieve the tracks reconstructed by the global tracking.
/// Works on the current event.
///
//____________________________________________________________________________
void AliEMCALReconstructor::Reconstruct(TTree* digitsTree, TTree* clustersTree) const
{
  AliCodeTimerAuto("",0)
    
  // Get input digits and put them in fgDigitsArr, clear the list before 
  ReadDigitsArrayFromTree(digitsTree);
  
  InitClusterizer();
  
  fgClusterizer->InitParameters();
  fgClusterizer->SetOutput(clustersTree);
  
  // Skip clusterization of LED events
  if (GetRecParam()->GetEventSpecie()!=AliRecoParam::kCalib)
  {
    if(fgDigitsArr && fgDigitsArr->GetEntries()) 
    {
      fgClusterizer->SetInput(digitsTree);
      
      //fgClusterizer->Digits2Clusters("deb all") ; //For debugging
      fgClusterizer->Digits2Clusters("");
      
      fgClusterizer->Clear();
      
    }//digits array exists and has somethind
  }//not a LED event
  
  clustersTree->Fill();	

  // Deleting the recpoints at the end of the reconstruction call
  fgClusterizer->DeleteRecPoints();
}

///
/// Conversion from raw data to EMCAL digits.
/// Works on a single-event basis.
///
//____________________________________________________________________________
void AliEMCALReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{
  rawReader->Reset() ; 
  
  for (int  i=0;i<fgTriggerData->GetEntriesFast();i++) 
  {
    ((AliEMCALTriggerData*)fgTriggerData->At(i))->SetMode(1);  
  }

  if(fgDigitsArr) fgDigitsArr->Clear("C");
  
  const int kNTRU = fGeom->GetNTotalTRU();
  TClonesArray *digitsTrg = new TClonesArray("AliEMCALTriggerRawDigit", kNTRU * 96);
  
  Int_t bufsize = 32000;
  digitsTree->Branch("EMCAL", &fgDigitsArr, bufsize);
  digitsTree->Branch("EMTRG", &digitsTrg, bufsize);
  
  // Skip calibration events do the rest
  Bool_t doFit = kTRUE;
  if ( !(GetRecParam()->FitLEDEvents()) && GetRecParam()->GetEventSpecie()==AliRecoParam::kCalib) doFit = kFALSE;
  
  if (doFit)
  {
    // must be done here because, in constructor, option is not yet known
    fgRawUtils->SetOption(GetOption());
    
    //  fgRawUtils->SetRawFormatHighLowGainFactor(GetRecParam()->GetHighLowGainFactor());
    
    //   fgRawUtils->SetRawFormatOrder(GetRecParam()->GetOrderParameter());
    //    fgRawUtils->SetRawFormatTau(GetRecParam()->GetTau());
    fgRawUtils->SetNoiseThreshold(GetRecParam()->GetNoiseThreshold());
    fgRawUtils->SetNPedSamples(GetRecParam()->GetNPedSamples());
    fgRawUtils->SetRemoveBadChannels(GetRecParam()->GetRemoveBadChannels());
    if (!fgRawUtils->GetFittingAlgorithm()) fgRawUtils->SetFittingAlgorithm(GetRecParam()->GetFittingAlgorithm());
    fgRawUtils->SetFALTROUsage(GetRecParam()->UseFALTRO());
    fgRawUtils->SetL1PhaseUsage(GetRecParam()->UseL1Phase());
    //  fgRawUtils->SetFALTROUsage(0);
 
    //fgRawUtils->SetTimeMin(GetRecParam()->GetTimeMin());
    //fgRawUtils->SetTimeMax(GetRecParam()->GetTimeMax());
    
    //  fgRawUtils->SetTimeMin(-99999 );
    //  fgRawUtils->SetTimeMax( 99999 );
    
    fgRawUtils->Raw2Digits(rawReader,fgDigitsArr,fPedestalData,digitsTrg,fgTriggerData);
    
  }//skip calibration event
  else
  {
    AliDebug(1," Calibration Event, skip!");
  }
  
  digitsTree->Fill();
  digitsTrg->Delete();
  delete digitsTrg;
  
}

///
/// Create AliESDCaloCells, AliESDCaloClusters and AliESDCaloTrigger entries 
/// from AliEMCALRecPoints, AliEMCALDigits. 
/// Called by AliReconstruct after Reconstruct() and global tracking and 
/// vertexing and V0. 
/// Works on the current event
///
//____________________________________________________________________________
void AliEMCALReconstructor::FillESD(TTree* digitsTree, TTree* clustersTree, 
				    AliESDEvent* esd) const
{
  //########################################
  // Trigger
  //########################################
  
  static int saveOnce[2] = {0,0};
  
  Int_t v0M[2] = {0, 0};
  
  AliESDVZERO* esdV0 = esd->GetVZEROData();
  
  if (esdV0) 
  {
    v0M[0] = esdV0->GetTriggerChargeC();
    v0M[1] = esdV0->GetTriggerChargeA();
  }
  else 
  {
    AliWarning("No V0 ESD! Run trigger processor w/ null V0 charges");
  }
  
  if (fgTriggerDigits && fgTriggerDigits->GetEntriesFast()) fgTriggerDigits->Delete();
  
  TBranch *branchtrg = digitsTree->GetBranch("EMTRG");
  
  if (!branchtrg) 
  { 
    AliError("Can't get the branch with the EMCAL trigger digits!");
    return;
  }
  
  branchtrg->SetAddress(&fgTriggerDigits);
  branchtrg->GetEntry(0);
  
  // Note: fgTriggerProcessor reset done at the end of this method
  //   fgTriggerProcessor->Digits2Trigger(fgTriggerDigits, v0M, fTriggerData);
  
  // Fill ESD
  AliESDCaloTrigger* trgESD = esd->GetCaloTrigger("EMCAL");
  
  if (trgESD) 
  {
    trgESD->Allocate(fgTriggerDigits->GetEntriesFast());
    int trgBitWord=0;
    for (Int_t i = 0; i < fgTriggerDigits->GetEntriesFast(); i++) {
      AliEMCALTriggerRawDigit* rdig = (AliEMCALTriggerRawDigit*)fgTriggerDigits->At(i);
      if (AliDebugLevel() > 999) rdig->Print("");
      
      Int_t px, py, trgBitWord=0;
      if (fGeom->GetPositionInEMCALFromAbsFastORIndex(rdig->GetId(), px, py)) {
        Int_t a = -1, t = -1, times[10]; 
        
        rdig->GetMaximum(a, t);
        rdig->GetL0Times(times);
        
        trgESD->Add(px, py, a, t, times, rdig->GetNL0Times(), rdig->GetL1TimeSum(), rdig->GetL1SubRegion(), rdig->GetTriggerBits());
        trgBitWord |= rdig->GetTriggerBits();
        if (rdig->GetNL0Times()) trgBitWord |= 1 << 5; // Overwrite L0 id from STU payload
      }
    }
    
    trgESD->SetTriggerBitWord(trgBitWord);
    
    for (int i=0;i<fgTriggerData->GetEntriesFast();i++) 
    {
      for (int j=0;j<2;j++) {
        trgESD->SetL1Threshold(i,2 * j    , ((AliEMCALTriggerData*)fgTriggerData->At(i))->GetL1JetThreshold(  j));
        trgESD->SetL1Threshold(i,2 * j + 1, ((AliEMCALTriggerData*)fgTriggerData->At(i))->GetL1GammaThreshold(j));
      }
      trgESD->SetL1FrameMask(i,((AliEMCALTriggerData*)fgTriggerData->At(i))->GetL1FrameMask());
    
      Int_t v0[2];
      ((AliEMCALTriggerData*)fgTriggerData->At(i))->GetL1V0(v0);
      trgESD->SetL1V0(i,v0);        

      trgESD->SetMedian(i,((AliEMCALTriggerData*)fgTriggerData->At(i))->GetMedian());
      if (!saveOnce[i] && ((AliEMCALTriggerData*)fgTriggerData->At(i))->GetL1DataDecoded()) 
      {
        int type[19] = {0};
        ((AliEMCALTriggerData*)fgTriggerData->At(i))->GetL1TriggerType(type);
        
        esd->SetCaloTriggerType(i,type);
        
        saveOnce[i] = 1;
      }
    }
  }
  
  // Resetting
  for (int i=0;i<fgTriggerData->GetEntriesFast();i++) 
  {
    ((AliEMCALTriggerData*)fgTriggerData->At(i))->Reset();
  }
  
  //########################################
  //##############Fill CaloCells############
  //########################################
  
  // Get input digits and put them in fgDigitsArr, clear the list before 
  ReadDigitsArrayFromTree(digitsTree);
  
  Int_t nDigits = fgDigitsArr->GetEntries(), idignew = 0 ;
  AliDebug(1,Form("%d digits",nDigits));
  
  AliESDCaloCells &emcCells = *(esd->GetEMCALCells());
  emcCells.CreateContainer(nDigits);
  emcCells.SetType(AliVCaloCells::kEMCALCell);
  
  Float_t energy = 0;
  Float_t time   = 0;
  Int_t   mapDigitAndCellIndex[20000] = { 0 } ; // needed to pack cell MC labels in cluster and not loop all digits
  
  for (Int_t idig = 0 ; idig < nDigits ; idig++) 
  {
    const AliEMCALDigit * dig = (const AliEMCALDigit*)fgDigitsArr->At(idig);
    time   = dig->GetTime();      // Time already calibrated in clusterizer
    energy = dig->GetAmplitude(); // energy calibrated in clusterizer
    
    if (energy <= 0 ) continue ;
    
    fgClusterizer->Calibrate(energy,time,dig->GetId()); // Digits already calibrated in clusterizers
      
    if (energy <= 0 ) continue ; // Digits tagged as bad (dead, hot, not alive) are set to 0 in calibrate, remove them
       
    // Only for MC
    // Get the label of the primary particle that generated the cell
    // Assign the particle that deposited more energy
    Int_t   nparents   = dig->GetNiparent();
    Int_t   digLabel   =-1 ;
    Float_t edep       =-1.;
    if ( nparents > 0 )
    {
      for ( Int_t jndex = 0 ; jndex < nparents ; jndex++ ) 
      { 
         if(edep >= dig->GetDEParent(jndex+1)) continue ;
            
         digLabel = dig->GetIparent (jndex+1);
         edep     = dig->GetDEParent(jndex+1);          
      } // all primaries in digit      
    } // select primary label
        
    Bool_t highGain = kFALSE;
    if( dig->GetType() == AliEMCALDigit::kHG ) highGain = kTRUE;
        
    emcCells.SetCell(idignew,dig->GetId(),energy, time,digLabel,0.,highGain);
        
    mapDigitAndCellIndex[idignew] = idig; // needed to pack cell MC labels in cluster 
        
    idignew++;
  }
  
  emcCells.SetNumberOfCells(idignew);
  emcCells.Sort();
  
  //------------------------------------------------------------
  //-----------------CLUSTERS-----------------------------------
  //------------------------------------------------------------
  clustersTree->SetBranchStatus("*",0); //disable all branches
  clustersTree->SetBranchStatus("EMCALECARP",1); //Enable only the branch we need
  if(fgClustersArr) fgClustersArr->Clear();
  TBranch *branch = clustersTree->GetBranch("EMCALECARP");
  branch->SetAddress(&fgClustersArr);
  branch->GetEntry(0);
  //clustersTree->GetEvent(0);
  
  Int_t nClusters = fgClustersArr->GetEntries(),  nClustersNew=0;
  AliDebug(1,Form("%d clusters",nClusters));

  
  //########################################
  //##############Fill CaloClusters#########
  //########################################
  for (Int_t iClust = 0 ; iClust < nClusters ; iClust++) 
  {
    const AliEMCALRecPoint * clust = (const AliEMCALRecPoint*)fgClustersArr->At(iClust);
    
    if(!clust) continue;
    
    // clust->Print(); //For debugging
    
    // Get information from EMCAL reconstruction points
    Float_t xyz[3];
    TVector3 gpos;
    clust->GetGlobalPosition(gpos);
    for (Int_t ixyz=0; ixyz<3; ixyz++)
      xyz[ixyz] = gpos[ixyz];
    Float_t elipAxis[2];
    clust->GetElipsAxis(elipAxis);
    
    //Create digits lists
    Int_t cellMult = clust->GetMultiplicity();

    Float_t *amplFloat = clust->GetEnergiesList();
    Int_t   *digitInts = clust->GetAbsId();
    TArrayS absIdList(cellMult);
    TArrayD fracList(cellMult);

    Int_t newCellMult = 0;
    for (Int_t iCell=0; iCell<cellMult; iCell++) 
    {
      if (amplFloat[iCell] > 0) 
      {
        absIdList[newCellMult] = (UShort_t)(digitInts[iCell]);
        
        //Calculate Fraction
        if(emcCells.GetCellAmplitude(digitInts[iCell])>0 && GetRecParam()->GetUnfold())
          fracList[newCellMult] = amplFloat[iCell]/(emcCells.GetCellAmplitude(digitInts[iCell]));//get cell calibration value 
        else
          fracList[newCellMult] = 0; 

        newCellMult++;
      }
    }
    
    absIdList.Set(newCellMult);
    fracList.Set(newCellMult);
    
    // Accept cluster if it has some digit
    if(newCellMult > 0) 
    { 
      nClustersNew++;
      
      // Fill the ESDCaloCluster
      AliESDCaloCluster * ec = new AliESDCaloCluster() ;
      ec->SetType(AliVCluster::kEMCALClusterv1);
      ec->SetPosition(xyz);
      ec->SetE(clust->GetEnergy());
      ec->SetTOF(clust->GetTime()) ; //time-of-flight

      ec->SetChi2(-1); //not yet implemented

      // Distance to the nearest bad tower
      ec->SetDistanceToBadChannel(clust->GetDistanceToBadTower()); 
      
      //
      // Shower shape related
      //
      ec->SetNCells(newCellMult);
      ec->SetDispersion(clust->GetDispersion());
      ec->SetM02(elipAxis[0]*elipAxis[0]) ;
      ec->SetM20(elipAxis[1]*elipAxis[1]) ;
      ec->SetNExMax(clust->GetNExMax()); //number of local maxima
  
      //
      // Primaries
      //
      Int_t    parentMult   = 0;
      Int_t   *parentList   =  clust->GetParents(parentMult); // label index
      Float_t *parentListDE =  clust->GetParentsDE();         // deposited energy

      ec->SetLabel(parentList,parentMult);    
     
      ec->SetClusterMCEdepFractionFromEdepArray(parentListDE);
     
      //
      // Change type of cells list from short to ushort
      // Pack the at maximum 4 different primary contributions on cells into a single int.
      // (Profit from the same loop on cells in cluster)
      //
      UShort_t *newAbsIdList  = new UShort_t[newCellMult];
      Double_t *newFracList   = new Double_t[newCellMult];
      UInt_t   *cellEDepFrac  = new UInt_t  [newCellMult];

      for(Int_t i = 0; i < newCellMult ; i++) 
      {
        // int to short
        newAbsIdList[i] = absIdList[i];
        newFracList [i] = fracList[i];
        cellEDepFrac[i] = 0;
        
        // *Pack cells parent energy deposit*
        if(parentMult > 0)
        {          
          // Get the digit that originated this cell cluster
          Int_t   cellPos = emcCells.GetCellPosition (newAbsIdList[i]);
          Int_t   idigit  = mapDigitAndCellIndex[cellPos];
          const AliEMCALDigit * dig = (const AliEMCALDigit*)fgDigitsArr->At(idigit);
          
          // Find the 4 MC labels that contributed to the cluster and their 
          // deposited energy in the current digit
          Int_t  nparents   = dig->GetNiparent();
          if ( nparents > 0 ) 
          {
            Int_t   digLabel   =-1 ; 
            Float_t edep       = 0 ;
            Float_t edepTot    = 0 ;
            Float_t mcEDepFrac[4] = {0,0,0,0};
            
            for ( Int_t jndex = 0 ; jndex < nparents ; jndex++ ) 
            { // all primaries in digit
              digLabel = dig->GetIparent (jndex+1);
              edep     = dig->GetDEParent(jndex+1);
              edepTot += edep;
              
              if       ( digLabel == parentList[0] ) mcEDepFrac[0] = edep; 
              else  if ( digLabel == parentList[1] ) mcEDepFrac[1] = edep;
              else  if ( digLabel == parentList[2] ) mcEDepFrac[2] = edep;
              else  if ( digLabel == parentList[3] ) mcEDepFrac[3] = edep;
            } // all primaries in digit
            
            // Divide energy deposit by total deposited energy
            // Do not take the stored digit energy since it is smeared and with noise
            // One could go back to the sdigit, but it is simpler to take the added 
            // deposited energy of all primaries.
            // Do this only when deposited energy is significant, use 10 MeV although 50 MeV should be expected
            if(edepTot > 0.01) 
            {
              for(Int_t idep = 0; idep < 4; idep++) 
                mcEDepFrac[idep] /= edepTot;
              
              cellEDepFrac[i] = ec->PackMCEdepFraction(mcEDepFrac);
            }
          } // at least one parent in digit
        } // at least one parent in cluster, do the cell primary packing
      } // cell cluster loop 
      
      ec->SetCellsAbsId(newAbsIdList);
      ec->SetCellsAmplitudeFraction(newFracList);
      ec->SetCellsMCEdepFractionMap(cellEDepFrac);
      
      //
      // Track matching
      //
      fMatches->Clear();
      Int_t nTracks = esd->GetNumberOfTracks();
      for (Int_t itrack = 0; itrack < nTracks; itrack++)
      {
        AliESDtrack * track = esd->GetTrack(itrack) ; // retrieve track
        if(track->GetEMCALcluster()==iClust)
        {
          Float_t dEta=-999, dPhi=-999;
          Bool_t isMatch =  CalculateResidual(track, ec, dEta, dPhi);
          if(!isMatch) 
          {
            // AliDebug(10, "Not good");
            continue;
          }
          
          AliEMCALMatch *match = new AliEMCALMatch();
          match->SetIndexT(itrack);
          match->SetDistance(TMath::Sqrt(dEta*dEta+dPhi*dPhi));
          match->SetdEta(dEta);
          match->SetdPhi(dPhi);
          fMatches->Add(match);
        }
      } 
      
      fMatches->Sort(kSortAscending); //Sort matched tracks from closest to furthest
      Int_t nMatch = fMatches->GetEntries();
      TArrayI arrayTrackMatched(nMatch);
      for(Int_t imatch=0; imatch<nMatch; imatch++)
      {
        AliEMCALMatch *match = (AliEMCALMatch*)fMatches->At(imatch);
        arrayTrackMatched[imatch] = match->GetIndexT();
        if(imatch==0)
        {
          ec->SetTrackDistance(match->GetdPhi(), match->GetdEta());
        }
      }
      
      ec->AddTracksMatched(arrayTrackMatched);
      
      //add the cluster to the esd object
      esd->AddCaloCluster(ec);

      delete ec;
      delete [] newAbsIdList ;
      delete [] newFracList ;
      delete [] cellEDepFrac ;
    }
  } // cycle on clusters

  //
  // Reset the index of matched cluster for tracks
  // to the one in CaloCluster array
  //
  Int_t ncls = esd->GetNumberOfCaloClusters();
  for(Int_t icl=0; icl<ncls; icl++)
  {
    AliESDCaloCluster *cluster = esd->GetCaloCluster(icl);
    if(!cluster || !cluster->IsEMCAL()) continue;
    TArrayI *trackIndex = cluster->GetTracksMatched();
    for(Int_t itr=0; itr<trackIndex->GetSize(); itr++)
    {
      AliESDtrack *track = esd->GetTrack(trackIndex->At(itr));
      track->SetEMCALcluster(cluster->GetID());
    }
  }
  
  // Fill ESDCaloCluster with PID weights
  AliEMCALPID *pid = new AliEMCALPID;
  //pid->SetPrintInfo(kTRUE);
  pid->SetReconstructor(kTRUE);
  pid->RunPID(esd);
  delete pid;
  
  // Store EMCAL misalignment matrixes
  FillMisalMatrixes(esd) ;
  
}

///
/// Store EMCAL matrixes in ESD Header.
///
//==================================================================================
void AliEMCALReconstructor::FillMisalMatrixes(AliESDEvent* esd)const
{  
  // Check, if matrixes was already stored
  for(Int_t sm = 0 ; sm < fGeom->GetNumberOfSuperModules(); sm++)
  {
    if(esd->GetEMCALMatrix(sm)!=0)
      return ;
  }
  
  // Create and store matrixes
  if(!gGeoManager)
  {
    AliError("Can not store misal. matrixes: no gGeoManager! \n") ;
    return ;
  }
  
  // Note, that owner of copied matrixes will be header
  const Int_t bufsize = 255;
  char path[bufsize] ;
  TGeoHMatrix * m = 0x0;
  Int_t tmpType = -1;
  Int_t SMOrder = 0;
  TString SMName;
  for(Int_t sm = 0; sm < fGeom->GetNumberOfSuperModules(); sm++)
  {
    if(fGeom->GetSMType(sm) == AliEMCALGeometry::kEMCAL_Standard )      SMName = "SMOD";
    else if(fGeom->GetSMType(sm) == AliEMCALGeometry::kEMCAL_Half )     SMName = "SM10";
    else if(fGeom->GetSMType(sm) == AliEMCALGeometry::kEMCAL_3rd )      SMName = "SM3rd";
    else if( fGeom->GetSMType(sm) == AliEMCALGeometry::kDCAL_Standard ) SMName = "DCSM";
    else if( fGeom->GetSMType(sm) == AliEMCALGeometry::kDCAL_Ext )      SMName = "DCEXT";
    else AliError("Unkown SM Type!!");

    if(fGeom->GetSMType(sm) == tmpType) 
    {
      SMOrder++;
    } 
    else 
    {
      tmpType = fGeom->GetSMType(sm);
      SMOrder = 1;
    }
    snprintf(path,bufsize,"/ALIC_1/XEN1_1/%s_%d", SMName.Data(), SMOrder) ;

    if (gGeoManager->CheckPath(path))
    {
      gGeoManager->cd(path);
      m = gGeoManager->GetCurrentMatrix() ;
      //			printf("================================================= \n");
      //			printf("AliEMCALReconstructor::FixMisalMatrixes(), sm %d, \n",sm);
      //			m->Print("");
      esd->SetEMCALMatrix(new TGeoHMatrix(*m),sm) ;
      //			printf("================================================= \n");
    }
    else
    {
      esd->SetEMCALMatrix(NULL,sm) ;
    }
  }
}

///
/// Read the digits from the input tree
/// See AliEMCALClusterizer::SetInput(TTree *digitsTree);    
///
//__________________________________________________________________________
void AliEMCALReconstructor::ReadDigitsArrayFromTree(TTree *digitsTree) const
{
  // Clear previous digits in the list
  if(fgDigitsArr)
  { 
    fgDigitsArr->Clear("C");
  }
  else
  {
    // It should not happen, but just in case ...
    fgDigitsArr = new TClonesArray("AliEMCALDigit",100); 
  }
  
  // Read the digits from the input tree
  TBranch *branch = digitsTree->GetBranch("EMCAL");
  if (!branch)
  { 
    AliError("can't get the branch with the EMCAL digits !");
    return;
  }  
  
  branch->SetAddress(&fgDigitsArr);
  branch->GetEntry(0);
}

///
/// Calculate the residual between track and cluster.
///
/// If the esdFriend is available, use the TPCOuter point as the starting 
/// point of extrapolation, otherwise use the TPCInner point.
///
//==================================================================================
Bool_t AliEMCALReconstructor::CalculateResidual(AliESDtrack *track, AliESDCaloCluster *cluster, Float_t &dEta, Float_t &dPhi)const
{
  dEta = -999, dPhi = -999;
  Bool_t ITSTrackSA = 0;

  AliExternalTrackParam *trkParam = 0;
  
  const AliESDfriendTrack*  friendTrack = track->GetFriendTrack();
  if(friendTrack && friendTrack->GetTPCOut())
    trkParam = const_cast<AliExternalTrackParam*>(friendTrack->GetTPCOut());
  else if(track->GetInnerParam())
    trkParam = const_cast<AliExternalTrackParam*>(track->GetInnerParam());
  else
  {
    trkParam = new AliExternalTrackParam(*track); //If there is ITSSa track 
    ITSTrackSA = 1;	
  }
  if(!trkParam) return kFALSE;
  
  AliExternalTrackParam trkParamTmp (*trkParam);
  if(!AliEMCALRecoUtils::ExtrapolateTrackToCluster(&trkParamTmp, cluster, track->GetMass(kTRUE), GetRecParam()->GetExtrapolateStep(), dEta, dPhi))
  {
    if(ITSTrackSA) delete trkParam;
    return kFALSE;
  }

  if(ITSTrackSA) delete trkParam;
  return kTRUE;
}

///
/// Default constructor.
///
//==================================================================================
AliEMCALReconstructor::AliEMCALMatch::AliEMCALMatch() 
  : TObject(),  
    fIndexT(-1), 
    fDistance(-999.),
    fdEta(-999.),
    fdPhi(-999.)
{
}

///
/// Copy constructor.
///
//==================================================================================
AliEMCALReconstructor::AliEMCALMatch::AliEMCALMatch(const AliEMCALMatch& copy)
  : TObject(),
    fIndexT(copy.fIndexT),
    fDistance(copy.fDistance),
    fdEta(copy.fdEta),
    fdPhi(copy.fdPhi)
{
}

///
/// Assignment operator; use copy ctor.
///
//_____________________________________________________________________
AliEMCALReconstructor::AliEMCALMatch& AliEMCALReconstructor::AliEMCALMatch::AliEMCALMatch::operator = (const AliEMCALMatch &source)
{ 
  if (&source == this) return *this;

  new (this) AliEMCALMatch(source);
  return *this;
}

///
/// Compare wrt the residual.
///
//==================================================================================
Int_t AliEMCALReconstructor::AliEMCALMatch::Compare(const TObject *obj) const 
{	
  AliEMCALReconstructor::AliEMCALMatch *that = (AliEMCALReconstructor::AliEMCALMatch*)obj;
	
  Double_t thisDist = fDistance;//fDistance;
  Double_t thatDist = that->fDistance;//that->GetDistance();
	
  if (thisDist > thatDist) return 1;
  else if (thisDist < thatDist) return -1;
  return 0;
}
