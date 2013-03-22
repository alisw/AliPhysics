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
// Author: Gustavo Conesa Balbastre,
//         Adapted from analysis class from Deepa Thomas
//
// $Id$
//_________________________________________________________________________

// --- Root ---
#include <TString.h>
#include <TRefArray.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TROOT.h>
#include <TInterpreter.h>
#include <TFile.h>

// --- AliRoot Analysis Steering
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliGeomManager.h"
#include "AliVCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliLog.h"
#include "AliVEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliOADBContainer.h"
#include "AliAODMCParticle.h"

// --- EMCAL
#include "AliEMCALAfterBurnerUF.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALClusterizerNxN.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALClusterizerv2.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALDigit.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALCalibData.h"

#include "AliAnalysisTaskEMCALClusterize.h"

ClassImp(AliAnalysisTaskEMCALClusterize)

//______________________________________________________________________________
AliAnalysisTaskEMCALClusterize::AliAnalysisTaskEMCALClusterize(const char *name) 
: AliAnalysisTaskSE(name)
, fEvent(0)
, fGeom(0),               fGeomName("") 
, fGeomMatrixSet(kFALSE), fLoadGeomMatrices(kFALSE)
, fCalibData(0),          fPedestalData(0)
, fOCDBpath(""),          fAccessOCDB(kFALSE)
, fDigitsArr(0),          fClusterArr(0),             fCaloClusterArr(0)
, fRecParam(0),           fClusterizer(0)
, fUnfolder(0),           fJustUnfold(kFALSE) 
, fOutputAODBranch(0),    fOutputAODBranchName(""),   fOutputAODBranchSet(0)
, fFillAODFile(kFALSE),   fFillAODHeader(0)
, fFillAODCaloCells(0),   fRun(-1)
, fRecoUtils(0),          fConfigName("")
, fOrgClusterCellId()
, fCellLabels(),          fCellSecondLabels(),        fCellTime()
, fCellMatchdEta(),       fCellMatchdPhi()
, fRecalibrateWithClusterTime(0)
, fMaxEvent(0),           fDoTrackMatching(kFALSE)
, fSelectCell(kFALSE),    fSelectCellMinE(0),         fSelectCellMinFrac(0)
, fRemoveLEDEvents(kTRUE),fRemoveExoticEvents(kFALSE)
, fImportGeometryFromFile(kFALSE), fImportGeometryFilePath("") 
, fOADBSet(kFALSE),       fAccessOADB(kTRUE),         fOADBFilePath("")
, fCentralityClass(""),   fSelectEMCALEvent(0)
, fEMCALEnergyCut(0.),    fEMCALNcellsCut (0)
, fSetCellMCLabelFromCluster(0)
, fRemapMCLabelForAODs(0)
{
  // Constructor
  
  for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;
  
  ResetArrays();
  
  fCentralityBin[0] = fCentralityBin[1]=-1;
  
}

//______________________________________________________________
AliAnalysisTaskEMCALClusterize::AliAnalysisTaskEMCALClusterize() 
: AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskEMCALClusterize")
, fEvent(0)
, fGeom(0),                 fGeomName("") 
, fGeomMatrixSet(kFALSE),   fLoadGeomMatrices(kFALSE)
, fCalibData(0),            fPedestalData(0)
, fOCDBpath(""),            fAccessOCDB(kFALSE)
, fDigitsArr(0),            fClusterArr(0),             fCaloClusterArr(0)
, fRecParam(0),             fClusterizer(0)
, fUnfolder(0),             fJustUnfold(kFALSE) 
, fOutputAODBranch(0),      fOutputAODBranchName(""),   fOutputAODBranchSet(0)
, fFillAODFile(kFALSE),     fFillAODHeader(0)
, fFillAODCaloCells(0),     fRun(-1)
, fRecoUtils(0),            fConfigName("")
, fOrgClusterCellId()
, fCellLabels(),            fCellSecondLabels(),        fCellTime()
, fCellMatchdEta(),         fCellMatchdPhi()
, fRecalibrateWithClusterTime(0)
, fMaxEvent(0),             fDoTrackMatching(kFALSE)
, fSelectCell(kFALSE),      fSelectCellMinE(0),         fSelectCellMinFrac(0)
, fRemoveLEDEvents(kTRUE),  fRemoveExoticEvents(kFALSE)
, fImportGeometryFromFile(kFALSE), fImportGeometryFilePath("")
, fOADBSet(kFALSE),         fAccessOADB(kTRUE),        fOADBFilePath("")
, fCentralityClass(""),     fSelectEMCALEvent(0)
, fEMCALEnergyCut(0.),      fEMCALNcellsCut (0)
, fSetCellMCLabelFromCluster(0)
, fRemapMCLabelForAODs(0)
{
  // Constructor
  
  for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;
  
  ResetArrays();
  
  fCentralityBin[0] = fCentralityBin[1]=-1;

}


//_______________________________________________________________
AliAnalysisTaskEMCALClusterize::~AliAnalysisTaskEMCALClusterize()
{
  //dtor 
  
  if (fDigitsArr)
  {
    fDigitsArr->Clear("C");
    delete fDigitsArr; 
  }
  
  if (fClusterArr)
  {
    fClusterArr->Delete();
    delete fClusterArr;
  }
  
  if (fCaloClusterArr)
  {
    fCaloClusterArr->Delete();
    delete fCaloClusterArr; 
  }
  
  if(fClusterizer) delete fClusterizer;
  if(fUnfolder)    delete fUnfolder;   
  if(fRecoUtils)   delete fRecoUtils;
  
}

//_______________________________________________________
Bool_t AliAnalysisTaskEMCALClusterize::AcceptEventEMCAL()
{
  // Accept event given there is a EMCAL cluster with enough energy, and not noisy, exotic
  
  if(!fSelectEMCALEvent)   return kTRUE; // accept
  
  if(fEMCALEnergyCut <= 0) return kTRUE; // accept
  
  Int_t           nCluster = InputEvent() -> GetNumberOfCaloClusters();
  AliVCaloCells * caloCell = InputEvent() -> GetEMCALCells();
  Int_t           bc       = InputEvent() -> GetBunchCrossNumber();

  for(Int_t icalo = 0; icalo < nCluster; icalo++)
  {
    AliVCluster *clus = (AliVCluster*) (InputEvent()->GetCaloCluster(icalo));
    
    if( ( clus->IsEMCAL() ) && ( clus->GetNCells() > fEMCALNcellsCut ) && ( clus->E() > fEMCALEnergyCut ) &&
       fRecoUtils->IsGoodCluster(clus,fGeom,caloCell,bc))
    {
      
      if (fDebug > 0)
        printf("AliAnalysisTaskEMCALClusterize::AcceptEventEMCAL() - Accept :  E %2.2f > %2.2f, nCells %d > %d \n",
                             clus->E(), fEMCALEnergyCut, clus->GetNCells(), fEMCALNcellsCut);
      
      return kTRUE;
    }
    
  }// loop
  
  if (fDebug > 0)
    printf("AliAnalysisTaskEMCALClusterize::AcceptEventEMCAL() - Reject \n");
  
  return kFALSE;
  
}  

//_______________________________________________
void AliAnalysisTaskEMCALClusterize::AccessOADB()
{
  // Set the AODB calibration, bad channels etc. parameters at least once
  // alignment matrices from OADB done in SetGeometryMatrices
  
  //Set it only once
  if(fOADBSet) return ; 
  
  Int_t   runnumber = InputEvent()->GetRunNumber() ;
  TString pass      = GetPass();
  
  printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Get AODB parameters from EMCAL in %s for run %d, and <%s> \n",fOADBFilePath.Data(),runnumber,pass.Data());
  
  Int_t nSM = fGeom->GetNumberOfSuperModules();
  
  // Bad map
  if(fRecoUtils->IsBadChannelsRemovalSwitchedOn())
  {
    AliOADBContainer *contBC=new AliOADBContainer("");
    contBC->InitFromFile(Form("%s/EMCALBadChannels.root",fOADBFilePath.Data()),"AliEMCALBadChannels"); 
    
    TObjArray *arrayBC=(TObjArray*)contBC->GetObject(runnumber);
    
    if(arrayBC)
    {
        fRecoUtils->SwitchOnDistToBadChannelRecalculation();
        printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Remove EMCAL bad cells \n");
        
        for (Int_t i=0; i<nSM; ++i) 
        {
          TH2I *hbm = fRecoUtils->GetEMCALChannelStatusMap(i);
          
          if (hbm)
            delete hbm;
          
          hbm=(TH2I*)arrayBC->FindObject(Form("EMCALBadChannelMap_Mod%d",i));
          
          if (!hbm) 
          {
            AliError(Form("Can not get EMCALBadChannelMap_Mod%d",i));
            continue;
          }
          
          hbm->SetDirectory(0);
          fRecoUtils->SetEMCALChannelStatusMap(i,hbm);
          
        } // loop
    } else printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Do NOT remove EMCAL bad channels\n"); // run array
  }  // Remove bad
  
  // Energy Recalibration
  if(fRecoUtils->IsRecalibrationOn())
  {
    AliOADBContainer *contRF=new AliOADBContainer("");
    
    contRF->InitFromFile(Form("%s/EMCALRecalib.root",fOADBFilePath.Data()),"AliEMCALRecalib");
    
    TObjArray *recal=(TObjArray*)contRF->GetObject(runnumber); 
    
    if(recal)
    {
      TObjArray *recalpass=(TObjArray*)recal->FindObject(pass);
      
      if(recalpass)
      {
        TObjArray *recalib=(TObjArray*)recalpass->FindObject("Recalib");
        
        if(recalib)
        {
          printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Recalibrate EMCAL \n");
          for (Int_t i=0; i<nSM; ++i) 
          {
            TH2F *h = fRecoUtils->GetEMCALChannelRecalibrationFactors(i);
            
            if (h)
              delete h;
            
            h = (TH2F*)recalib->FindObject(Form("EMCALRecalFactors_SM%d",i));
            
            if (!h) 
            {
              AliError(Form("Could not load EMCALRecalFactors_SM%d",i));
              continue;
            }
            
            h->SetDirectory(0);
            
            fRecoUtils->SetEMCALChannelRecalibrationFactors(i,h);
          } // SM loop
        }else printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Do NOT recalibrate EMCAL, no params object array \n"); // array ok
      }else printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Do NOT recalibrate EMCAL, no params for pass\n"); // array pass ok
    }else printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Do NOT recalibrate EMCAL, no params for run\n");  // run number array ok
        
  } // Recalibration on
  
  // Energy Recalibration, apply on top of previous calibration factors
  if(fRecoUtils->IsRunDepRecalibrationOn())
  {
    AliOADBContainer *contRFTD=new AliOADBContainer("");
    
    contRFTD->InitFromFile(Form("%s/EMCALTemperatureCorrCalib.root",fOADBFilePath.Data()),"AliEMCALRunDepTempCalibCorrections");
    
    TH1S *htd=(TH1S*)contRFTD->GetObject(runnumber); 
    
    //If it did not exist for this run, get closes one
    if (!htd)
    {
      AliWarning(Form("No TemperatureCorrCalib Objects for run: %d",runnumber));
      // let's get the closest runnumber instead then..
      Int_t lower = 0;
      Int_t ic = 0;
      Int_t maxEntry = contRFTD->GetNumberOfEntries();
      
      while ( (ic < maxEntry) && (contRFTD->UpperLimit(ic) < runnumber) ) {
        lower = ic;
        ic++;
      }
      
      Int_t closest = lower;
      if ( (ic<maxEntry) &&
          (contRFTD->LowerLimit(ic)-runnumber) < (runnumber - contRFTD->UpperLimit(lower)) ) {
        closest = ic;
      }
      
      AliWarning(Form("TemperatureCorrCalib Objects found closest id %d from run: %d", closest, contRFTD->LowerLimit(closest)));
      htd = (TH1S*) contRFTD->GetObjectByIndex(closest);
    } 

    if(htd)
    {
      printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Recalibrate (Temperature) EMCAL \n");
      
      for (Int_t ism=0; ism<nSM; ++ism) 
      {        
        for (Int_t icol=0; icol<48; ++icol) 
        {        
          for (Int_t irow=0; irow<24; ++irow) 
          {
            Float_t factor = fRecoUtils->GetEMCALChannelRecalibrationFactor(ism,icol,irow);
            
            Int_t absID = fGeom->GetAbsCellIdFromCellIndexes(ism, irow, icol); // original calibration factor
            factor *= htd->GetBinContent(absID) / 10000. ; // correction dependent on T
            //printf("\t ism %d, icol %d, irow %d,absID %d, corrA %2.3f, corrB %2.3f, corrAB %2.3f\n",ism, icol, irow, absID, 
            //      GetEMCALChannelRecalibrationFactor(ism,icol,irow) , htd->GetBinContent(absID) / 10000., factor);
            fRecoUtils->SetEMCALChannelRecalibrationFactor(ism,icol,irow,factor);
          } // columns
        } // rows 
      } // SM loop
    }else printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Do NOT recalibrate EMCAL with T variations, no params TH1 \n"); 
  } // Run by Run T calibration
  
  // Time Recalibration
  if(fRecoUtils->IsTimeRecalibrationOn())
  {
    AliOADBContainer *contTRF=new AliOADBContainer("");
    
    contTRF->InitFromFile(Form("%s/EMCALTimeCalib.root",fOADBFilePath.Data()),"AliEMCALTimeCalib");
    
    TObjArray *trecal=(TObjArray*)contTRF->GetObject(runnumber); 
    
    if(trecal)
    {
      TObjArray *trecalpass=(TObjArray*)trecal->FindObject(pass);

      if(trecalpass)
      {
        printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Time Recalibrate EMCAL \n");
        for (Int_t ibc = 0; ibc < 4; ++ibc) 
        {
          TH1F *h = fRecoUtils->GetEMCALChannelTimeRecalibrationFactors(ibc);
          
          if (h)
            delete h;
          
          h = (TH1F*)trecalpass->FindObject(Form("hAllTimeAvBC%d",ibc));
          
          if (!h) 
          {
            AliError(Form("Could not load hAllTimeAvBC%d",ibc));
            continue;
          }
          
          h->SetDirectory(0);
          
          fRecoUtils->SetEMCALChannelTimeRecalibrationFactors(ibc,h);
        } // bunch crossing loop
      }else printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Do NOT recalibrate time EMCAL, no params for pass\n"); // array pass ok
    }else printf("AliAnalysisTaskEMCALClusterize::SetOADBParameters() - Do NOT recalibrate time EMCAL, no params for run\n");  // run number array ok
    
  } // Time recalibration on    
  
  // Parameters already set once, so do not it again
  fOADBSet = kTRUE;
  
}  

//_________________________________________________
Bool_t AliAnalysisTaskEMCALClusterize::AccessOCDB()
{
  //Access to OCDB stuff
  
  fEvent = InputEvent();
  if (!fEvent)
  {
    Warning("AccessOCDB","Event not available!!!");
    return kFALSE;
  }
  
  if (fEvent->GetRunNumber()==fRun)
    return kTRUE;
  fRun = fEvent->GetRunNumber();
  
  if(DebugLevel() > 1 )
    printf("AliAnalysisTaksEMCALClusterize::AccessOCDB() - Begin");
  
  AliCDBManager *cdb = AliCDBManager::Instance();
  
  
  if (fOCDBpath.Length())
  {
    cdb->SetDefaultStorage(fOCDBpath.Data());
    printf("AliAnalysisTaksEMCALClusterize::AccessOCDB() - Default storage %s",fOCDBpath.Data());
  }
  
  cdb->SetRun(fEvent->GetRunNumber());
  
  //
  // EMCAL from RAW OCDB
  if (fOCDBpath.Contains("alien:"))
  {
    cdb->SetSpecificStorage("EMCAL/Calib/Data","alien://Folder=/alice/data/2010/OCDB");
    cdb->SetSpecificStorage("EMCAL/Calib/Pedestals","alien://Folder=/alice/data/2010/OCDB");
  }
  
  TString path = cdb->GetDefaultStorage()->GetBaseFolder();
  
  // init parameters:
  
  //Get calibration parameters	
  if(!fCalibData)
  {
    AliCDBEntry *entry = (AliCDBEntry*) 
    AliCDBManager::Instance()->Get("EMCAL/Calib/Data");
    
    if (entry) fCalibData =  (AliEMCALCalibData*) entry->GetObject();
  }
  
  if(!fCalibData)
    AliFatal("Calibration parameters not found in CDB!");
  
  //Get calibration parameters	
  if(!fPedestalData)
  {
    AliCDBEntry *entry = (AliCDBEntry*) 
    AliCDBManager::Instance()->Get("EMCAL/Calib/Pedestals");
    
    if (entry) fPedestalData =  (AliCaloCalibPedestal*) entry->GetObject();
  }
  
  if(!fPedestalData)
    AliFatal("Dead map not found in CDB!");
  
  return kTRUE;
}

//_____________________________________________________
void AliAnalysisTaskEMCALClusterize::CheckAndGetEvent()
{
  // Get the input event, it can depend in embedded events what you want to get
  // Also check if the quality of the event is good if not reject it
  
  fEvent = 0x0;
  
  //Process events if there is a high energy cluster
  if(!AcceptEventEMCAL())  return ; 
  
  AliAODInputHandler* aodIH = dynamic_cast<AliAODInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  Int_t eventN = Entry();
  if(aodIH) eventN = aodIH->GetReadEntry(); 
  
  if (eventN > fMaxEvent) 
    return ;
  
  //printf("Clusterizer --- Event %d-- \n",eventN);
  
  //Check if input event are embedded events
  //If so, take output event
  if (aodIH && aodIH->GetMergeEvents()) 
  {
    fEvent  = AODEvent();
    
    if(!aodIH->GetMergeEMCALCells()) 
      AliFatal("Events merged but not EMCAL cells, check analysis settings!");
    
    if(DebugLevel() > 1)
    {
      printf("AliAnalysisTaksEMCALClusterize::UserExec() - Use embedded events\n");
      
      printf("\t InputEvent  N Clusters %d, N Cells %d\n",InputEvent()->GetNumberOfCaloClusters(),
             InputEvent()->GetEMCALCells()->GetNumberOfCells());
      
      printf("\t MergedEvent  N Clusters %d, N Cells %d\n",aodIH->GetEventToMerge()->GetNumberOfCaloClusters(),
             aodIH->GetEventToMerge()->GetEMCALCells()->GetNumberOfCells());
      
      for (Int_t icl=0; icl < aodIH->GetEventToMerge()->GetNumberOfCaloClusters(); icl++) 
      {
        AliAODCaloCluster *sigCluster = aodIH->GetEventToMerge()->GetCaloCluster(icl);
        if(sigCluster->IsEMCAL()) printf("\t \t Signal cluster: i %d, E  %f\n",icl,sigCluster->E());
      }
      
      printf("\t OutputEvent N Clusters %d, N Cells %d\n", AODEvent()->GetNumberOfCaloClusters(),
             AODEvent()->GetEMCALCells()->GetNumberOfCells());
    }
  }
  else 
  {
    fEvent =  InputEvent();
    if(fFillAODCaloCells) FillAODCaloCells();   
    if(fFillAODHeader)    FillAODHeader();
  }
  
  if (!fEvent) 
  {
    Error("UserExec","Event not available");
    return ;
  }
  
  //-------------------------------------------------------------------------------------
  // Reject events if LED was firing, use only for LHC11a data 
  // Reject event if triggered by exotic cell and remove exotic cells if not triggered
  //-------------------------------------------------------------------------------------
  
  if( IsLEDEvent( InputEvent()->GetRunNumber() ) ) { fEvent = 0x0 ; return ; }
  
  if( IsExoticEvent() )                            { fEvent = 0x0 ; return ; }
  
  //-------------------------------------------------------------------------------------
  // Set the cluster array in the event (output or input)
  //-------------------------------------------------------------------------------------
  
  if     ( fFillAODFile ) 
  {
    //Magic line to write events to AOD filem put after event rejection
    AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
  }
  else if( !fOutputAODBranchSet )
  {
    // Create array and put it in the input event, if output AOD not selected, only once
    InputEvent()->AddObject(fOutputAODBranch);
    fOutputAODBranchSet = kTRUE;
    printf("AliAnalysisTaskEMCALClusterize::UserExec() - Add AOD branch <%s> to input event\n",fOutputAODBranchName.Data());
  }
  
}

//____________________________________________________
void AliAnalysisTaskEMCALClusterize::ClusterizeCells()
{ 
  // Recluster calocells, transform them into digits, 
  // feed the clusterizer with them and get new list of clusters
  
  //In case of MC, first loop on the clusters and fill MC label to array  
  Int_t nClusters     = fEvent->GetNumberOfCaloClusters();
  Int_t nClustersOrg  = 0;
  
  AliAODInputHandler* aodIH = dynamic_cast<AliAODInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if(aodIH && aodIH->GetEventToMerge())  //Embedding
    nClusters = aodIH->GetEventToMerge()->GetNumberOfCaloClusters(); //Get clusters directly from embedded signal

  ResetArrays();
   
  for (Int_t i = 0; i < nClusters; i++)
  {
    AliVCluster *clus = 0;
    if(aodIH && aodIH->GetEventToMerge()) //Embedding
      clus = aodIH->GetEventToMerge()->GetCaloCluster(i); //Get clusters directly from embedded signal
    else      
      clus = fEvent->GetCaloCluster(i);
    
    if(!clus) return;
    
    if(clus->IsEMCAL())
    {
      Int_t label = clus->GetLabel();
      Int_t label2 = -1 ;
      //printf("Org cluster E %f, Time  %e, Index = %d, ID %d, MC label %d\n", clus->E(), clus->GetTOF(),i, clus->GetID(),label );
      //printf("Original list of labels from old cluster : \n");
      //for(Int_t imc = 0; imc < clus->GetNLabels(); imc++) printf("\t Label %d\n",clus->GetLabelAt(imc));
      
      if (clus->GetNLabels()>=2) label2 = clus->GetLabelAt(1) ;
      UShort_t * index    = clus->GetCellsAbsId() ;
      for(Int_t icell=0; icell < clus->GetNCells(); icell++ )
      {
        //printf("\t cell %d, MC label %d\n",index[icell],fEvent->GetEMCALCells()->GetCellMCLabel(index[icell]));
        fOrgClusterCellId[index[icell]] = i;
        fCellLabels[index[icell]]       = label;
        fCellSecondLabels[index[icell]] = label2;
        fCellTime[index[icell]]         = clus->GetTOF();
        fCellMatchdEta[index[icell]]    = clus->GetTrackDz();
        fCellMatchdPhi[index[icell]]    = clus->GetTrackDx();
      }
      nClustersOrg++;
    }
    // printf("\n");
  } 
  
  // Transform CaloCells into Digits
  
  Int_t    idigit =  0;
  Int_t    id     = -1;
  Float_t  amp    = -1; 
  Double_t time   = -1; 
  
  AliVCaloCells *cells = fEvent->GetEMCALCells();
  
  Int_t bc = InputEvent()->GetBunchCrossNumber();

  for (Int_t icell = 0; icell < cells->GetNumberOfCells(); icell++)
  {
    // Get cell values, recalibrate and not include bad channels found in analysis, nor cells with too low energy, nor exotic cell
    id = cells->GetCellNumber(icell);
    Bool_t accept = fRecoUtils->AcceptCalibrateCell(id,bc,amp,time,cells);
    
    // Do not include cells with too low energy, nor exotic cell
    if( amp  < fRecParam->GetMinECut() ||
        time > fRecParam->GetTimeMax() ||
        time < fRecParam->GetTimeMin()    ) accept = kFALSE;
    
    // In case of old AOD analysis cell time is -1 s, approximate replacing by time of the cluster the digit belongs.
    if (fRecalibrateWithClusterTime)
    { 
      time = fCellTime[id];
      //printf("cell %d time org %f - ",id, time*1.e9);
      fRecoUtils->RecalibrateCellTime(id,bc,time);
      //printf("recal %f\n",time*1.e9);
    }
    
    //Exotic?
    if (accept && fRecoUtils->IsExoticCell(id,cells,bc))
        accept = kFALSE;
    
    if( !accept )
    {
      if( DebugLevel() > 2 )
        printf("AliAnalysisTaksEMCALClusterize::ClusterizeCells() - Remove channel absId %d, index %d of %d, amp %f, time %f\n",
               id,icell, cells->GetNumberOfCells(), amp, time*1.e9);
      continue;
    }
    
    Int_t mcLabel = cells->GetMCLabel(icell);
    //printf("AliAnalysisTaksEMCALClusterize::ClusterizeCells() - cell %d, mc label %d\n",id,mcLabel);

    //if(fCellLabels[id]!=mcLabel)printf("mcLabel %d - %d\n",mcLabel,fCellLabels[id]);
    if     ( fSetCellMCLabelFromCluster == 1 ) mcLabel = fCellLabels[id]; // Older aliroot MC productions
    else if( fSetCellMCLabelFromCluster == 0 && fRemapMCLabelForAODs) RemapMCLabelForAODs(mcLabel);
    else mcLabel = -1; // found later
    
    //printf("\t new label %d\n",mcLabel);
    
    // Create the digit, put a fake primary deposited energy to trick the clusterizer
    // when checking the most likely primary
    
    Float_t efrac = cells->GetEFraction(icell);
    
    //When checking the MC of digits, give weight to cells with embedded signal
    if (mcLabel > 0 && efrac < 1.e-6) efrac = 1;
    
    //printf("******* Cell %d, id %d, e %f,  fraction %f, MC label %d, used MC label %d\n",icell,id,amp,cells->GetEFraction(icell),cells->GetMCLabel(icell),mcLabel);
    
    new((*fDigitsArr)[idigit]) AliEMCALDigit( mcLabel, mcLabel, id, amp, time,AliEMCALDigit::kHG,idigit, 0, 0, amp*efrac);
    // Last parameter should be MC deposited energy, since it is not available, add just the cell amplitude so that
    // we give more weight to the MC label of the cell with highest energy in the cluster
        
    idigit++;
  }
  
  fDigitsArr->Sort();
  
  //-------------------------------------------------------------------------------------
  //Do the clusterization
  //-------------------------------------------------------------------------------------        
  
  fClusterizer->Digits2Clusters("");
  
  //-------------------------------------------------------------------------------------
  //Transform the recpoints into AliVClusters
  //-------------------------------------------------------------------------------------
  
  RecPoints2Clusters();
  
  if(!fCaloClusterArr)
  { 
    printf("AliAnalysisTaksEMCALClusterize::UserExec() - No array with CaloClusters, input RecPoints entries %d\n",fClusterArr->GetEntriesFast());
    return;    
  }
  
  if( DebugLevel() > 0 )
  {
    printf("AliAnalysisTaksEMCALClusterize::ClusterizeCells() - N clusters: before recluster %d, after recluster %d\n",nClustersOrg, fCaloClusterArr->GetEntriesFast());
    
    if(fCaloClusterArr->GetEntriesFast() != fClusterArr->GetEntriesFast())
    {
      printf("\t Some RecRoints not transformed into CaloClusters (clusterizer %d, unfold %d): Input entries %d - Output entries %d - %d (not fast)\n",
             fRecParam->GetClusterizerFlag(),fRecParam->GetUnfold(),
             fClusterArr->GetEntriesFast(), fCaloClusterArr->GetEntriesFast(), fCaloClusterArr->GetEntries());
    }
  }
}

//_____________________________________________________
void AliAnalysisTaskEMCALClusterize::ClusterUnfolding()
{
  // Take the event clusters and unfold them
  
  AliVCaloCells *cells   = fEvent->GetEMCALCells();
  Double_t cellAmplitude = 0;
  Double_t cellTime      = 0;
  Short_t  cellNumber    = 0;
  Int_t    cellMCLabel   = 0;
  Double_t cellEFrac     = 0;
  Int_t    nClustersOrg  = 0;
  
  // Fill the array with the EMCAL clusters, copy them
  for (Int_t i = 0; i < fEvent->GetNumberOfCaloClusters(); i++)
  {
    AliVCluster *clus = fEvent->GetCaloCluster(i);
    if(clus->IsEMCAL())
    {        
      //recalibrate/remove bad channels/etc if requested
      if(fRecoUtils->ClusterContainsBadChannel(fGeom,clus->GetCellsAbsId(), clus->GetNCells()))
      {
        continue;
      } 
      
      if(fRecoUtils->IsRecalibrationOn())
      {
        //Calibrate cluster
        fRecoUtils->RecalibrateClusterEnergy(fGeom, clus, cells);
        
        //CalibrateCells
        for (Int_t icell = 0; icell < cells->GetNumberOfCells(); icell++)
        {
          if (cells->GetCell(icell, cellNumber, cellAmplitude, cellTime, cellMCLabel, cellEFrac) != kTRUE)
            break;
          
          Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
          fGeom->GetCellIndex(cellNumber,imod,iTower,iIphi,iIeta); 
          fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);	
          
          //Do not include bad channels found in analysis?
          if( fRecoUtils->IsBadChannelsRemovalSwitchedOn() && 
              fRecoUtils->GetEMCALChannelStatus(imod, ieta, iphi))
            continue;
          
          cells->SetCell(icell, cellNumber, cellAmplitude*fRecoUtils->GetEMCALChannelRecalibrationFactor(imod,ieta,iphi),cellTime);
          
        }// cells loop            
      }// recalibrate
      
      //Cast to ESD or AOD, needed to create the cluster array
      AliESDCaloCluster * esdCluster = dynamic_cast<AliESDCaloCluster*> (clus);
      AliAODCaloCluster * aodCluster = dynamic_cast<AliAODCaloCluster*> (clus);
      
      if     (esdCluster)
      {
        fCaloClusterArr->Add( new AliESDCaloCluster(*esdCluster) );   
      }//ESD
      else if(aodCluster)
      {
        fCaloClusterArr->Add( new AliAODCaloCluster(*aodCluster) );   
      }//AOD
      else 
        Warning("UserExec()"," - Wrong CaloCluster type?");
      
      nClustersOrg++;
    }
  }
  
  //Do the unfolding
  fUnfolder->UnfoldClusters(fCaloClusterArr, cells);
  
  //CLEAN-UP
  fUnfolder->Clear();
  
}

//_____________________________________________________
void AliAnalysisTaskEMCALClusterize::FillAODCaloCells() 
{
  // Put calo cells in standard branch  
  AliVCaloCells &eventEMcells = *(fEvent->GetEMCALCells());
  Int_t nEMcell = eventEMcells.GetNumberOfCells() ;
  
  AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
  aodEMcells.CreateContainer(nEMcell);
  aodEMcells.SetType(AliVCaloCells::kEMCALCell);
  Double_t calibFactor = 1.;   
  for (Int_t iCell = 0; iCell < nEMcell; iCell++) 
  { 
    Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
    fGeom->GetCellIndex(eventEMcells.GetCellNumber(iCell),imod,iTower,iIphi,iIeta); 
    fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);	
    
    if(fRecoUtils->IsRecalibrationOn())
    { 
      calibFactor = fRecoUtils->GetEMCALChannelRecalibrationFactor(imod,ieta,iphi);
    }
    
    if(!fRecoUtils->GetEMCALChannelStatus(imod, ieta, iphi))
    { //Channel is not declared as bad
      aodEMcells.SetCell(iCell,eventEMcells.GetCellNumber(iCell),eventEMcells.GetAmplitude(iCell)*calibFactor,
                         eventEMcells.GetTime(iCell),eventEMcells.GetMCLabel(iCell),eventEMcells.GetEFraction(iCell));
    }
    else 
    {
      aodEMcells.SetCell(iCell,eventEMcells.GetCellNumber(iCell),0,-1,-1,0);
    }
  }
  aodEMcells.Sort();
  
}

//__________________________________________________
void AliAnalysisTaskEMCALClusterize::FillAODHeader() 
{
  //Put event header information in standard AOD branch
  
  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*> (fEvent);
  AliAODEvent* aodevent = dynamic_cast<AliAODEvent*> (fEvent);
  
  Double_t pos[3]   ;
  Double_t covVtx[6];
  for (Int_t i = 0; i < 6; i++)  covVtx[i] = 0.;
  
  AliAODHeader* header = AODEvent()->GetHeader();
  header->SetRunNumber(fEvent->GetRunNumber());
  
  if(esdevent)
  {
    TTree* tree = fInputHandler->GetTree();
    if (tree) 
    {
      TFile* file = tree->GetCurrentFile();
      if (file) header->SetESDFileName(file->GetName());
    }
  }
  else if (aodevent) header->SetESDFileName(aodevent->GetHeader()->GetESDFileName());
  
  header->SetBunchCrossNumber(fEvent->GetBunchCrossNumber());
  header->SetOrbitNumber(fEvent->GetOrbitNumber());
  header->SetPeriodNumber(fEvent->GetPeriodNumber());
  header->SetEventType(fEvent->GetEventType());
  
  //Centrality
  if(fEvent->GetCentrality())
    header->SetCentrality(new AliCentrality(*(fEvent->GetCentrality())));
  else
    header->SetCentrality(0);
  
  //Trigger  
  header->SetOfflineTrigger(fInputHandler->IsEventSelected()); // propagate the decision of the physics selection
  if      (esdevent) header->SetFiredTriggerClasses(esdevent->GetFiredTriggerClasses());
  else if (aodevent) header->SetFiredTriggerClasses(aodevent->GetFiredTriggerClasses());
  
  header->SetTriggerMask(fEvent->GetTriggerMask()); 
  header->SetTriggerCluster(fEvent->GetTriggerCluster());
  
  if      (esdevent)
  {
    header->SetL0TriggerInputs(esdevent->GetHeader()->GetL0TriggerInputs());    
    header->SetL1TriggerInputs(esdevent->GetHeader()->GetL1TriggerInputs());    
    header->SetL2TriggerInputs(esdevent->GetHeader()->GetL2TriggerInputs());    
  }
  else if (aodevent)
  {
    header->SetL0TriggerInputs(aodevent->GetHeader()->GetL0TriggerInputs());    
    header->SetL1TriggerInputs(aodevent->GetHeader()->GetL1TriggerInputs());    
    header->SetL2TriggerInputs(aodevent->GetHeader()->GetL2TriggerInputs());    
  }
  
  header->SetMagneticField(fEvent->GetMagneticField());
  //header->SetMuonMagFieldScale(esdevent->GetCurrentDip()/6000.); 
  
  header->SetZDCN1Energy(fEvent->GetZDCN1Energy());
  header->SetZDCP1Energy(fEvent->GetZDCP1Energy());
  header->SetZDCN2Energy(fEvent->GetZDCN2Energy());
  header->SetZDCP2Energy(fEvent->GetZDCP2Energy());
  header->SetZDCEMEnergy(fEvent->GetZDCEMEnergy(0),fEvent->GetZDCEMEnergy(1));
  
  Float_t diamxy[2]={fEvent->GetDiamondX(),fEvent->GetDiamondY()};
  Float_t diamcov[3];
  fEvent->GetDiamondCovXY(diamcov);
  header->SetDiamond(diamxy,diamcov);
  if      (esdevent) header->SetDiamondZ(esdevent->GetDiamondZ(),esdevent->GetSigma2DiamondZ());
  else if (aodevent) header->SetDiamondZ(aodevent->GetDiamondZ(),aodevent->GetSigma2DiamondZ());
  
  //
  Int_t nVertices = 1 ;/* = prim. vtx*/;
  Int_t nCaloClus = fEvent->GetNumberOfCaloClusters();
  
  AODEvent()->ResetStd(0, nVertices, 0, 0, 0, nCaloClus, 0, 0);
  
  // Access to the AOD container of vertices
  TClonesArray &vertices = *(AODEvent()->GetVertices());
  Int_t jVertices=0;
  
  // Add primary vertex. The primary tracks will be defined
  // after the loops on the composite objects (V0, cascades, kinks)
  fEvent->GetPrimaryVertex()->GetXYZ(pos);
  Float_t chi = 0;
  if      (esdevent)
  {
    esdevent->GetPrimaryVertex()->GetCovMatrix(covVtx);
    chi = esdevent->GetPrimaryVertex()->GetChi2toNDF();
  }
  else if (aodevent)
  {
    aodevent->GetPrimaryVertex()->GetCovMatrix(covVtx);
    chi = aodevent->GetPrimaryVertex()->GetChi2perNDF();//Different from ESD?
  }
  
  AliAODVertex * primary = new(vertices[jVertices++])
  AliAODVertex(pos, covVtx, chi, NULL, -1, AliAODVertex::kPrimary);
  primary->SetName(fEvent->GetPrimaryVertex()->GetName());
  primary->SetTitle(fEvent->GetPrimaryVertex()->GetTitle());
  
}

//___________________________________________________________
void AliAnalysisTaskEMCALClusterize::FillCaloClusterInEvent()
{
  // Get the CaloClusters array, do some final calculations 
  // and put the clusters in the output or input event 
  // as a separate branch 
  
  //First recalculate track-matching for the new clusters
  if(fDoTrackMatching) 
  {
    fRecoUtils->FindMatches(fEvent,fCaloClusterArr,fGeom);
  }
  //Put the new clusters in the AOD list
  
  Int_t kNumberOfCaloClusters   = fCaloClusterArr->GetEntriesFast();

  for(Int_t i = 0; i < kNumberOfCaloClusters; i++)
  {
    AliAODCaloCluster *newCluster = (AliAODCaloCluster *) fCaloClusterArr->At(i);
    
    newCluster->SetID(i);
    
    // Correct cluster energy non linearity    
    newCluster->SetE(fRecoUtils->CorrectClusterEnergyLinearity(newCluster));
    
    //Add matched track
    if(fDoTrackMatching)
    {
      Int_t trackIndex = fRecoUtils->GetMatchedTrackIndex(i);
      if(trackIndex >= 0)
      {
        newCluster->AddTrackMatched(fEvent->GetTrack(trackIndex));
        if(DebugLevel() > 1) 
          printf("AliAnalysisTaksEMCALClusterize::UserExec() - Matched Track index %d to new cluster %d \n",trackIndex,i);
      }
      
      Float_t dR = 999., dZ = 999.;
      fRecoUtils->GetMatchedResiduals(newCluster->GetID(),dR,dZ);
      newCluster->SetTrackDistance(dR,dZ);
      
    }
    else 
    {// Assign previously assigned matched track in reco, very very rough
      Int_t absId0 = newCluster->GetCellsAbsId()[0]; // Assign match of first cell in cluster
      newCluster->SetTrackDistance(fCellMatchdPhi[absId0],fCellMatchdEta[absId0]);
    }
    
    //printf("New cluster E %f, Time  %e, Id = ", newCluster->E(), newCluster->GetTOF() );
    //for(Int_t icell=0; icell < newCluster->GetNCells(); icell++ ) printf(" %d,", newCluster->GetCellsAbsId() [icell] );
    //printf("\n");
    
    // Calculate distance to bad channel for new cluster. Make sure you give the list of bad channels.
    fRecoUtils->RecalculateClusterDistanceToBadChannel(fGeom, fEvent->GetEMCALCells(), newCluster);
    
    new((*fOutputAODBranch)[i])  AliAODCaloCluster(*newCluster);
    
    if(DebugLevel() > 1 )
      printf("AliAnalysisTaksEMCALClusterize::UserExec() - New cluster %d of %d, energy %f, mc label %d \n",newCluster->GetID(), kNumberOfCaloClusters, newCluster->E(), newCluster->GetLabel());
    
  } // cluster loop
  
  fOutputAODBranch->Expand(kNumberOfCaloClusters); // resize TObjArray to 'remove' slots
  
  
}

//_______________________________________________
TString AliAnalysisTaskEMCALClusterize::GetPass()
{
  // Get passx from filename.
  
  if (!AliAnalysisManager::GetAnalysisManager()->GetTree()) 
  {
    AliError("AliAnalysisTaskEMCALClusterize::GetPass() - Pointer to tree = 0, returning null\n");
    return TString("");
  }
  
  if (!AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()) 
  {
    AliError("AliAnalysisTaskEMCALClusterize::GetPass() - Null pointer input file, returning null\n");
    return TString("");
  }
  
  TString pass(AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()->GetName());
  if      (pass.Contains("ass1")) return TString("pass1");
  else if (pass.Contains("ass2")) return TString("pass2");
  else if (pass.Contains("ass3")) return TString("pass3");
  else if (pass.Contains("ass4")) return TString("pass4");
  else if (pass.Contains("ass5")) return TString("pass5");
  
  // No condition fullfilled, give a default value
  printf("AliAnalysisTaskEMCALClusterize::GetPass() - Pass number string not found \n");
  return TString("");            
  
}

//_________________________________________
void AliAnalysisTaskEMCALClusterize::Init()
{
  //Init analysis with configuration macro if available
  
  fOADBSet           = kFALSE;
  if(fOADBFilePath == "") fOADBFilePath = "$ALICE_ROOT/OADB/EMCAL" ;          
  
  fBranchNames       = "ESD:AliESDHeader.,EMCALCells.";
  
  if(!fRecParam)     fRecParam  = new AliEMCALRecParam;
  if(!fRecoUtils)    fRecoUtils = new AliEMCALRecoUtils();  
  
  if(fMaxEvent          <= 0) fMaxEvent          = 1000000000;
  if(fSelectCellMinE    <= 0) fSelectCellMinE    = 0.005;     
  if(fSelectCellMinFrac <= 0) fSelectCellMinFrac = 0.001;
  
  //Centrality
  if(fCentralityClass  == "") fCentralityClass  = "V0M";
  
  if (fOCDBpath            == "") fOCDBpath            = "raw://" ;
  if (fOutputAODBranchName == "") fOutputAODBranchName = "newEMCALClusters" ;
  
  if(gROOT->LoadMacro(fConfigName) >=0)
  {
    printf("AliAnalysisTaksEMCALClusterize::Init() - Configure analysis with %s\n",fConfigName.Data());
    AliAnalysisTaskEMCALClusterize *clus = (AliAnalysisTaskEMCALClusterize*)gInterpreter->ProcessLine("ConfigEMCALClusterize()");
    fGeomName         = clus->fGeomName; 
    fLoadGeomMatrices = clus->fLoadGeomMatrices;
    fOCDBpath         = clus->fOCDBpath;   
    fAccessOCDB       = clus->fAccessOCDB;
    fRecParam         = clus->fRecParam;
    fJustUnfold       = clus->fJustUnfold;
    fFillAODFile      = clus->fFillAODFile;
    fRecoUtils        = clus->fRecoUtils; 
    fConfigName       = clus->fConfigName;
    fMaxEvent         = clus->fMaxEvent;
    fDoTrackMatching  = clus->fDoTrackMatching;
    fOutputAODBranchName = clus->fOutputAODBranchName;
    for(Int_t i = 0; i < 12; i++) fGeomMatrix[i] = clus->fGeomMatrix[i] ;
    fCentralityClass  = clus->fCentralityClass;
    fCentralityBin[0] = clus->fCentralityBin[0];
    fCentralityBin[1] = clus->fCentralityBin[1];
  }
  
  // Init geometry, I do not like much to do it like this ...
  if(fImportGeometryFromFile && !gGeoManager) 
  {
    if (fImportGeometryFilePath == "") fImportGeometryFilePath = "$ALICE_ROOT/OADB/EMCAL/geometry_2011.root" ; // "$ALICE_ROOT/EVE/alice-data/default_geo.root"
    printf("AliAnalysisTaskEMCALClusterize::Init() - Import %s\n",fImportGeometryFilePath.Data());
    TGeoManager::Import(fImportGeometryFilePath) ; 
  }
  
}  

//_______________________________________________________
void AliAnalysisTaskEMCALClusterize::InitClusterization()
{
  //Select clusterization/unfolding algorithm and set all the needed parameters
  
  if (fJustUnfold)
  {
    // init the unfolding afterburner 
    delete fUnfolder;
    fUnfolder =  new AliEMCALAfterBurnerUF(fRecParam->GetW0(),fRecParam->GetLocMaxCut(),fRecParam->GetMinECut());
    return;
  }
  
  //First init the clusterizer
  delete fClusterizer;
  if     (fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv1)
    fClusterizer = new AliEMCALClusterizerv1 (fGeom, fCalibData, fPedestalData);
  else if(fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv2) 
    fClusterizer = new AliEMCALClusterizerv2(fGeom, fCalibData, fPedestalData);
  else if(fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerNxN)
  { 
    fClusterizer = new AliEMCALClusterizerNxN(fGeom, fCalibData, fPedestalData);
    fClusterizer->SetNRowDiff(fRecParam->GetNRowDiff());
    fClusterizer->SetNColDiff(fRecParam->GetNColDiff());
  } 
  else 
  {
    AliFatal(Form("Clusterizer < %d > not available", fRecParam->GetClusterizerFlag()));
  }
  
  //Now set the parameters
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
  
  // Initialize the cluster rec points and digits arrays and get them.
  fClusterizer->SetOutput(0);
  fClusterArr = const_cast<TObjArray *>(fClusterizer->GetRecPoints());
  fDigitsArr  = fClusterizer->GetDigits();
  
  //In case of unfolding after clusterization is requested, set the corresponding parameters
  if(fRecParam->GetUnfold())
  {
    Int_t i=0;
    for (i = 0; i < 8; i++) 
    {
      fClusterizer->SetSSPars(i, fRecParam->GetSSPars(i));
    }//end of loop over parameters
    
    for (i = 0; i < 3; i++) 
    {
      fClusterizer->SetPar5  (i, fRecParam->GetPar5(i));
      fClusterizer->SetPar6  (i, fRecParam->GetPar6(i));
    }//end of loop over parameters
    
    fClusterizer->InitClusterUnfolding();
    
  }// to unfold
}

//_________________________________________________
void AliAnalysisTaskEMCALClusterize::InitGeometry()
{
  // Init geometry and set the geometry matrix, for the first event, skip the rest
  // Also set once the run dependent calibrations
  
  if(fGeomMatrixSet) return;
  
  Int_t runnumber = InputEvent()->GetRunNumber() ;
  if (!fGeom)
  {
    if(fGeomName=="")
    {
      if     (runnumber < 140000) fGeomName = "EMCAL_FIRSTYEARV1";
      else if(runnumber < 171000) fGeomName = "EMCAL_COMPLETEV1";
      else                        fGeomName = "EMCAL_COMPLETE12SMV1";  
      printf("AliAnalysisTaskEMCALClusterize::InitGeometry() - Set EMCAL geometry name to <%s> for run %d\n",fGeomName.Data(),runnumber);
    }
    
		fGeom = AliEMCALGeometry::GetInstance(fGeomName);
    
		if(fDebug > 0)
    {
			printf("AliAnalysisTaskEMCALClusterize::InitGeometry(run=%d)",runnumber);
			if (!gGeoManager) printf(" - Careful!, gGeoManager not loaded, load misalign matrices");
			printf("\n");
		}
	} // geometry pointer did not exist before
  
  if(fLoadGeomMatrices)
  {
    printf("AliAnalysisTaskEMCALClusterize::InitGeometry() - Load user defined EMCAL geometry matrices\n");
    
    // OADB if available
    AliOADBContainer emcGeoMat("AliEMCALgeo");
    emcGeoMat.InitFromFile(Form("%s/EMCALlocal2master.root",fOADBFilePath.Data()),"AliEMCALgeo");
    TObjArray *matEMCAL=(TObjArray*)emcGeoMat.GetObject(runnumber,"EmcalMatrices");
    
    
    for(Int_t mod=0; mod < (fGeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
    {
      
      if (!fGeomMatrix[mod]) // Get it from OADB
      {
        if(fDebug > 1 ) 
          printf("AliAnalysisTaskEMCALClusterize::InitGeometry() - EMCAL matrices SM %d, %p\n",
                 mod,((TGeoHMatrix*) matEMCAL->At(mod)));
        //((TGeoHMatrix*) matEMCAL->At(mod))->Print();
        
        fGeomMatrix[mod] = (TGeoHMatrix*) matEMCAL->At(mod) ;
      }        
      
      if(fGeomMatrix[mod])
      {
        if(DebugLevel() > 1) 
          fGeomMatrix[mod]->Print();
        
        fGeom->SetMisalMatrix(fGeomMatrix[mod],mod) ;  
      }
      
      fGeomMatrixSet=kTRUE;
      
    }//SM loop
  }//Load matrices
  else if(!gGeoManager)
  {
    printf("AliAnalysisTaksEMCALClusterize::InitGeometry() - Get geo matrices from data");
    //Still not implemented in AOD, just a workaround to be able to work at least with ESDs	
    if(!strcmp(fEvent->GetName(),"AliAODEvent")) 
    {
      if(DebugLevel() > 1) 
        Warning("UserExec","Use ideal geometry, values geometry matrix not kept in AODs.");
    }//AOD
    else 
    {	
      if(DebugLevel() > 1) 
        printf("AliAnalysisTaksEMCALClusterize::InitGeometry() - AliAnalysisTaskEMCALClusterize Load Misaligned matrices.");
      
      AliESDEvent* esd = dynamic_cast<AliESDEvent*>(fEvent) ;
      
      if(!esd)
      {
        Error("InitGeometry"," - This event does not contain ESDs?");
        if(fFillAODFile) AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kFALSE);
        return;
      }
      
      for(Int_t mod=0; mod < (fGeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
      {
        if(DebugLevel() > 1) 
          esd->GetEMCALMatrix(mod)->Print();
        
        if(esd->GetEMCALMatrix(mod)) fGeom->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
        
      } 
      
      fGeomMatrixSet=kTRUE;
      
    }//ESD
  }//Load matrices from Data 
  
}

//____________________________________________________
Bool_t AliAnalysisTaskEMCALClusterize::IsExoticEvent()
{
  
  // Check if event is exotic, get an exotic cell and compare with triggered patch
  // If there is a match remove event ... to be completed, filled with something provisional
  
  if(!fRemoveExoticEvents) return kFALSE;
  
  // Loop on cells
  AliVCaloCells * cells = fEvent->GetEMCALCells();
  Float_t totCellE = 0;
  Int_t bc = InputEvent()->GetBunchCrossNumber();
  for(Int_t icell = 0; icell < cells->GetNumberOfCells(); icell++)
  {
    
    Float_t  ecell = 0 ;
    Double_t tcell = 0 ;
    
    Int_t absID   = cells->GetCellNumber(icell);
    Bool_t accept = fRecoUtils->AcceptCalibrateCell(absID,bc,ecell,tcell,cells);
    if(accept && !fRecoUtils->IsExoticCell(absID,cells,bc)) totCellE += ecell;
  }
  
  //  TString triggerclasses = "";
  //  if(esdevent) triggerclasses = esdevent             ->GetFiredTriggerClasses();
  //  else         triggerclasses = ((AliAODEvent*)event)->GetFiredTriggerClasses();
  //    //  
  //    printf("AliAnalysisTaskEMCALClusterize - reject event %d with cluster  - reject event with ncells in SM3 %d and SM4 %d\n",(Int_t)Entry(),ncellsSM3, ncellsSM4);
  //    if(fFillAODFile) AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kFALSE);;
  //    return;
  //  
  
  //printf("TotE cell %f\n",totCellE);
  if(totCellE < 1) return kTRUE;
  
  return kFALSE;
  
} 

//________________________________________________________________
Bool_t AliAnalysisTaskEMCALClusterize::IsLEDEvent(const Int_t run)
{
  //Check if event is LED
  
  if(!fRemoveLEDEvents) return kFALSE;
  
  //check events of LHC11a period
  if(run < 146858 || run > 146860) return kFALSE ;
  
  // Count number of cells with energy larger than 0.1 in SM3, cut on this number
  Int_t ncellsSM3 = 0;
  AliVCaloCells * cells = fEvent->GetEMCALCells();
  for(Int_t icell = 0; icell < cells->GetNumberOfCells(); icell++)
  {
    if(cells->GetAmplitude(icell) > 0.1 && cells->GetCellNumber(icell)/(24*48)==3) ncellsSM3++;
  }
  
  TString triggerclasses = "";
  
  AliESDEvent *esdevent = dynamic_cast<AliESDEvent*>(fEvent);
  if(esdevent) triggerclasses = esdevent              ->GetFiredTriggerClasses();
  else         triggerclasses = ((AliAODEvent*)fEvent)->GetFiredTriggerClasses();
  
  Int_t ncellcut = 21;
  if(triggerclasses.Contains("EMC")) ncellcut = 35;
  
  if( ncellsSM3 >= ncellcut)
  {
    printf("AliAnalysisTaksEMCALClusterize::IsLEDEvent() - reject event %d with ncells in SM3 %d\n",(Int_t)Entry(),ncellsSM3);
    if(fFillAODFile) AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kFALSE);;
    return kTRUE;
  }
  
  return kFALSE;
  
} 

//_______________________________________________________
void AliAnalysisTaskEMCALClusterize::RecPoints2Clusters()
{
  // Restore clusters from recPoints
  // Cluster energy, global position, cells and their amplitude 
  // fractions are restored
  
  Int_t j = 0;
  for(Int_t i = 0; i < fClusterArr->GetEntriesFast(); i++)
  {
    AliEMCALRecPoint *recPoint = (AliEMCALRecPoint*) fClusterArr->At(i);
    
    const Int_t ncells = recPoint->GetMultiplicity();
    Int_t ncellsTrue = 0;
    
    if(recPoint->GetEnergy() < fRecParam->GetClusteringThreshold()) continue;
    
    // cells and their amplitude fractions
    UShort_t   absIds[ncells];  
    Double32_t ratios[ncells];
    
    //For later check embedding
    AliAODInputHandler* aodIH = dynamic_cast<AliAODInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    
    Float_t clusterE = 0; 
    for (Int_t c = 0; c < ncells; c++) 
    {
      AliEMCALDigit *digit = (AliEMCALDigit*) fDigitsArr->At(recPoint->GetDigitsList()[c]);
      
      absIds[ncellsTrue] = digit->GetId();
      ratios[ncellsTrue] = recPoint->GetEnergiesList()[c]/digit->GetAmplitude();
            
      // In case of unfolding, remove digits with unfolded energy too low      
      if(fSelectCell)
      {
        if     (recPoint->GetEnergiesList()[c] < fSelectCellMinE || ratios[ncellsTrue] < fSelectCellMinFrac)  
        {
          if(DebugLevel() > 1)
          {
            printf("AliAnalysisTaksEMCALClusterize::RecPoints2Clusters() - Too small energy in cell of cluster: cluster cell %f, digit %f\n",
                   recPoint->GetEnergiesList()[c],digit->GetAmplitude());
          }
          
          continue;
          
        } // if cuts
      }// Select cells
      
      //Recalculate cluster energy and number of cluster cells in case any of the cells was rejected
      clusterE  +=recPoint->GetEnergiesList()[c];
      
      // In case of embedding, fill ratio with amount of signal, 
      if (aodIH && aodIH->GetMergeEvents()) 
      {
        //AliVCaloCells* inEMCALCells = InputEvent()->GetEMCALCells();
        AliVCaloCells* meEMCALCells = aodIH->GetEventToMerge()->GetEMCALCells();
        AliVCaloCells* ouEMCALCells = AODEvent()->GetEMCALCells();
        
        Float_t sigAmplitude = meEMCALCells->GetCellAmplitude(absIds[ncellsTrue]);
        //Float_t bkgAmplitude = inEMCALCells->GetCellAmplitude(absIds[ncellsTrue]);
        Float_t sumAmplitude = ouEMCALCells->GetCellAmplitude(absIds[ncellsTrue]);
        //printf("\t AbsID %d, amplitude : bkg %f, sigAmplitude %f, summed %f - %f\n",absIds[ncellsTrue], bkgAmplitude, sigAmplitude, sumAmplitude, digit->GetAmplitude());
        
        if(sumAmplitude > 0) ratios[ncellsTrue] = sigAmplitude/sumAmplitude;
        //printf("\t \t ratio %f\n",ratios[ncellsTrue]);
        
      }//Embedding
      
      ncellsTrue++;
      
    }// cluster cell loop
    
    if (ncellsTrue < 1) 
    {
      if (DebugLevel() > 1) 
        printf("AliAnalysisTaskEMCALClusterize::RecPoints2Clusters() - Skipping cluster with no cells avobe threshold E = %f, ncells %d\n",
               recPoint->GetEnergy(), ncells);
      continue;
    }
    
    //if(ncellsTrue != ncells) printf("Old E %f, ncells %d; New E %f, ncells %d\n",recPoint->GetEnergy(),ncells,clusterE,ncellsTrue);
    
    if(clusterE <  fRecParam->GetClusteringThreshold()) 
    {
      if (DebugLevel()>1)
        printf("AliAnalysisTaskEMCALClusterize::RecPoints2Clusters() - Remove cluster with energy below seed threshold %f\n",clusterE);
      continue;
    }
    
    TVector3 gpos;
    Float_t g[3];
    
    // calculate new cluster position
    
    recPoint->EvalGlobalPosition(fRecParam->GetW0(), fDigitsArr);
    recPoint->GetGlobalPosition(gpos);
    gpos.GetXYZ(g);
    
    // create a new cluster
    
    (*fCaloClusterArr)[j] = new AliAODCaloCluster() ;
    AliAODCaloCluster *clus = dynamic_cast<AliAODCaloCluster *>( fCaloClusterArr->At(j) ) ;
    j++;
    clus->SetType(AliVCluster::kEMCALClusterv1);
    clus->SetE(clusterE);
    clus->SetPosition(g);
    clus->SetNCells(ncellsTrue);
    clus->SetCellsAbsId(absIds);
    clus->SetCellsAmplitudeFraction(ratios);
    clus->SetChi2(-1); //not yet implemented
    clus->SetTOF(recPoint->GetTime()) ; //time-of-flight
    clus->SetNExMax(recPoint->GetNExMax()); //number of local maxima
    clus->SetDistanceToBadChannel(recPoint->GetDistanceToBadTower()); 
    
    if(ncells == ncellsTrue)
    {
      Float_t elipAxis[2];
      recPoint->GetElipsAxis(elipAxis);
      clus->SetM02(elipAxis[0]*elipAxis[0]) ;
      clus->SetM20(elipAxis[1]*elipAxis[1]) ;
      clus->SetDispersion(recPoint->GetDispersion());
    }
    else if(fSelectCell)
    {
      // In case some cells rejected, in unfolding case, recalculate
      // shower shape parameters and position
      if(DebugLevel() > 1) 
        printf("AliAnalysisTaskEMCALClusterize::RecPoints2Clusters() - Cells removed from cluster (ncells %d, ncellsTrue %d), recalculate Shower Shape\n",ncells,ncellsTrue);
      
      AliVCaloCells* cells = 0x0; 
      if (aodIH && aodIH->GetMergeEvents()) cells = AODEvent()  ->GetEMCALCells();
      else                                  cells = InputEvent()->GetEMCALCells();
      
      fRecoUtils->RecalculateClusterShowerShapeParameters(fGeom,cells,clus);
      fRecoUtils->RecalculateClusterPID(clus);
      fRecoUtils->RecalculateClusterPosition(fGeom,cells,clus); 
      
    }
    
    // MC

    if     ( fSetCellMCLabelFromCluster == 1 ) SetClustersMCLabelFrom2SelectedLabels(recPoint,clus) ;
    else if( fSetCellMCLabelFromCluster == 2 ) SetClustersMCLabelFromOriginalClusters(clus) ;
    else
    {
      // Normal case, trust what the clusterizer has found
      Int_t  parentMult = 0;
      Int_t *parentList = recPoint->GetParents(parentMult);
      clus->SetLabel(parentList, parentMult);
//      printf("Label list : ");
//      for(Int_t ilabel = 0; ilabel < parentMult; ilabel++ ) printf(" %d ",parentList[ilabel]);
//      printf("\n");
    }
    
  } // recPoints loop
  
}

//___________________________________________________________________________
void AliAnalysisTaskEMCALClusterize::RemapMCLabelForAODs(Int_t & label)
{
  // MC label for Cells not remapped after ESD filtering, do it here.

  if(label < 0) return ;
  
  AliAODEvent  * evt = dynamic_cast<AliAODEvent*> (fEvent) ;
  if(!evt) return ;
  
  TClonesArray * arr = dynamic_cast<TClonesArray*>(evt->FindListObject("mcparticles")) ;
  if(!arr) return ;
    
  if(label < arr->GetEntriesFast())
  {
    AliAODMCParticle * particle = dynamic_cast<AliAODMCParticle *>(arr->At(label));
    if(!particle) return ;
        
    if(label == particle->Label()) return ; // label already OK
    //else printf("AliAnalysisTaskEMCALClusterize::RemapMCLabelForAODs() - Label  %d - AOD stack %d \n",label, particle->Label());
  }
  //else printf("AliAnalysisTaskEMCALClusterize::RemapMCLabelForAODs() - Label  %d > AOD labels %d \n",label, arr->GetEntriesFast());
  
  // loop on the particles list and check if there is one with the same label
  for(Int_t ind = 0; ind < arr->GetEntriesFast(); ind++ )
  {
    AliAODMCParticle * particle = dynamic_cast<AliAODMCParticle *>(arr->At(ind));
    if(!particle) continue ;

    if(label == particle->Label())
    {
      label = ind;
      //printf("AliAnalysisTaskEMCALClusterize::RemapMCLabelForAODs() - New Label Index  %d \n",label);
      return;
    }
  }
  
  label = -1;
  
  //printf("AliAnalysisTaskEMCALClusterize::RemapMCLabelForAODs() - Label not found set to -1 \n");
 
}

//________________________________________________
void AliAnalysisTaskEMCALClusterize::ResetArrays()
{
  // Reset arrays containing information for all possible cells
  for(Int_t j = 0; j < 12672; j++)
  {
    fOrgClusterCellId[j] = -1;
    fCellLabels[j]       = -1;
    fCellSecondLabels[j] = -1;
    fCellTime[j]         =  0.;
    fCellMatchdEta[j]    = -999;
    fCellMatchdPhi[j]    = -999;
  }
}

//_____________________________________________________________________________________________________
void AliAnalysisTaskEMCALClusterize::SetClustersMCLabelFrom2SelectedLabels(AliEMCALRecPoint  * recPoint,
                                                                           AliAODCaloCluster * clus)
{
  // Set the cluster MC label, the digizer was filled with most likely MC label for all cells in original cluster
  // Now check the second most likely MC label and add it to the new cluster
  
  Int_t  parentMult = 0;
  Int_t *parentList = recPoint->GetParents(parentMult);
  clus->SetLabel(parentList, parentMult);
  
  //Write the second major contributor to each MC cluster.
  Int_t iNewLabel ;
  for ( Int_t iLoopCell = 0 ; iLoopCell < clus->GetNCells() ; iLoopCell++ )
  {
    
    Int_t idCell = clus->GetCellAbsId(iLoopCell) ;
    if(idCell>=0)
    {
      iNewLabel = 1 ; //iNewLabel makes sure we  don't write twice the same label.
      for ( UInt_t iLoopLabels = 0 ; iLoopLabels < clus->GetNLabels() ; iLoopLabels++ )
      {
        if ( fCellSecondLabels[idCell] == -1 )  iNewLabel = 0;  // -1 is never a good second label.
          if ( fCellSecondLabels[idCell] == clus->GetLabelAt(iLoopLabels) )  iNewLabel = 0;
            }
      if (iNewLabel == 1)
      {
        Int_t * newLabelArray = new Int_t[clus->GetNLabels()+1] ;
        for ( UInt_t iLoopNewLabels = 0 ; iLoopNewLabels < clus->GetNLabels() ; iLoopNewLabels++ )
        {
          newLabelArray[iLoopNewLabels] = clus->GetLabelAt(iLoopNewLabels) ;
        }
        
        newLabelArray[clus->GetNLabels()] = fCellSecondLabels[idCell] ;
        clus->SetLabel(newLabelArray,clus->GetNLabels()+1) ;
        delete [] newLabelArray;
      }
    }//positive cell number
  }
}

//___________________________________________________________________________________________________
void AliAnalysisTaskEMCALClusterize::SetClustersMCLabelFromOriginalClusters(AliAODCaloCluster * clus)
{
  // Get the original clusters that contribute to the new cluster, assign the labels of such clusters
  // to the new cluster.
  // Only approximatedly valid  when output are V1 clusters, handle with care
    
  TArrayI clArray(300) ; //Weird if more than a few clusters are in the origin ...
  clArray.Reset();
  Int_t nClu = 0;
  Int_t nLabTotOrg = 0;
  Float_t emax = -1;
  Int_t idMax = -1;
  
  AliVEvent * event = fEvent;
  
  //In case of embedding the MC signal is in the added event
  AliAODInputHandler* aodIH = dynamic_cast<AliAODInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if(aodIH && aodIH->GetEventToMerge())  //Embedding
    event = aodIH->GetEventToMerge(); //Get clusters directly from embedded signal

  
  //Find the clusters that originally had the cells
  for ( Int_t iLoopCell = 0 ; iLoopCell < clus->GetNCells() ; iLoopCell++ )
  {
    Int_t idCell = clus->GetCellAbsId(iLoopCell) ;
    
    if(idCell>=0)
    {
      Int_t idCluster = fOrgClusterCellId[idCell];
      
      Bool_t set = kTRUE;
      for(Int_t icl =0; icl < nClu; icl++)
      {
        if(((Int_t)clArray.GetAt(icl))==-1 ) continue;
        if( idCluster == ((Int_t)clArray.GetAt(icl)) ) set = kFALSE;
        //   printf("\t \t icell %d  Cluster in array %d, IdCluster %d, in array %d, set %d\n",
        //          iLoopCell,                 icl,  idCluster,((Int_t)clArray.GetAt(icl)),set);
      }
      if( set && idCluster >= 0)
      {
        clArray.SetAt(idCluster,nClu++);
        //printf("******** idCluster %d \n",idCluster);
        nLabTotOrg+=(event->GetCaloCluster(idCluster))->GetNLabels();

        //printf("Cluster in array %d, IdCluster %d\n",nClu-1,  idCluster);

        //Search highest E cluster
        AliVCluster * clOrg = event->GetCaloCluster(idCluster);
        //printf("\t E %f\n",clOrg->E());
        if(emax < clOrg->E())
        {
          emax  = clOrg->E();
          idMax = idCluster;
        }        
      }
    }
  }// cell loop

  if(nClu==0 || nLabTotOrg == 0)
  {
    //if(clus->E() > 0.25) printf("AliAnalysisTaskEMCALClusterize::SetClustersMCLabelFromOriginalClusters() - Check: N org clusters %d, n tot labels %d, cluster E %f,  n cells %d\n",nClu,nLabTotOrg,clus->E(), clus->GetNCells());
    //for(Int_t icell = 0; icell < clus->GetNCells(); icell++) printf("\t cell %d",clus->GetCellsAbsId()[icell]);
    //printf("\n");
  }
  
  // Put the first in the list the cluster with highest energy
  if(idMax != ((Int_t)clArray.GetAt(0))) // Max not at first position
  {
    Int_t maxIndex = -1;
    Int_t firstCluster = ((Int_t)clArray.GetAt(0));
    for ( Int_t iLoopCluster = 0 ; iLoopCluster < nClu ; iLoopCluster++ )
    {
      if(idMax == ((Int_t)clArray.GetAt(iLoopCluster))) maxIndex = iLoopCluster;
    }
    
    if(firstCluster >=0 && idMax >=0)
    {
      clArray.SetAt(idMax,0);
      clArray.SetAt(firstCluster,maxIndex);
    }
  }
  
  // Get the labels list in the original clusters, assign all to the new cluster
  TArrayI clMCArray(nLabTotOrg) ;
  clMCArray.Reset();
  
  Int_t nLabTot = 0;
  for ( Int_t iLoopCluster = 0 ; iLoopCluster < nClu ; iLoopCluster++ )
  {
    Int_t idCluster = clArray.GetAt(iLoopCluster);
    //printf("New Cluster in Array %d,  idCluster %d \n",iLoopCluster,idCluster);
    AliVCluster * clOrg = event->GetCaloCluster(idCluster);
    Int_t nLab = clOrg->GetNLabels();
    
    for ( Int_t iLab = 0 ; iLab < nLab ; iLab++ )
    {
      Int_t lab = clOrg->GetLabelAt(iLab) ;
      if(lab>=0)
      {
        Bool_t set = kTRUE;
        //printf("\t \t Set Label %d \n", lab);
        for(Int_t iLabTot =0; iLabTot < nLabTot; iLabTot++)
        {
          if( lab == ((Int_t)clMCArray.GetAt(iLabTot)) ) set = kFALSE;
          //printf("iLoopCluster %d, Label ID in Org Cluster %d,label %d  Label ID in array %d, label in array %d, set %d\n",
          //       iLoopCluster,                           iLab,     lab,           iLabTot,  ((Int_t)clMCArray.GetAt(iLabTot)),set);
        }
        if( set ) clMCArray.SetAt(lab,nLabTot++);
      }
    }
  }// cluster loop
  
  // Set the final list of labels
  
  clus->SetLabel(clMCArray.GetArray(), nLabTot);
  
//  printf("Final list of labels for new cluster : \n");
//  for(Int_t ice = 0; ice < clus->GetNCells() ; ice++)
//  {
//    printf("\t Cell %d ",clus->GetCellsAbsId()[ice]);
//    Int_t label = InputEvent()->GetEMCALCells()->GetCellMCLabel(clus->GetCellsAbsId()[ice]);
//    printf(" org %d ",label);
//    RemapMCLabelForAODs(label);
//    printf(" new %d \n",label);
//  }
//  for(Int_t imc = 0; imc < clus->GetNLabels(); imc++) printf("\t Label %d\n",clus->GetLabelAt(imc));
}


//____________________________________________________________
void AliAnalysisTaskEMCALClusterize::UserCreateOutputObjects()
{
  // Init geometry, create list of output clusters
  

  fOutputAODBranch = new TClonesArray("AliAODCaloCluster", 0);

  if(fOutputAODBranchName.Length()==0)
  {
    fOutputAODBranchName = "newEMCALClustersArray";
    printf("Cluster branch name not set, set it to newEMCALClustersArray \n");
  }
  
  fOutputAODBranch->SetName(fOutputAODBranchName);
  
  if( fFillAODFile )
  {
    //fOutputAODBranch = new TClonesArray("AliAODCaloCluster", 0);

    
    //fOutputAODBranch->SetOwner(kFALSE);
    
    AddAODBranch("TClonesArray", &fOutputAODBranch);
  }
  
}

//_______________________________________________________
void AliAnalysisTaskEMCALClusterize::UserExec(Option_t *) 
{
  // Do clusterization event by event, execute different steps
  // 1) Do some checks on the kind of events (ESD, AOD) or if some filtering is needed, initializations
  //    load things and clear arrays
  // 2) Clusterize a) just unfolding existing clusters (fJustUnfold)
  //               b) recluster cells
  //                   +  convert cells into digits (calibrating them)
  //                   +  recluster digits into recPoints with the chosen clusterizer (V1, V1+Unfold,V2, NxN) 
  //                      with methods in AliEMCALClusterizer
  //                   +  transform recPoints into CaloClusters
  // 3) Do some final calculations in the found clusters (track-matching) and put them in an AOD branch
  
  //-------
  // Step 1
  

  //Remove the contents of AOD branch output list set in the previous event 
  fOutputAODBranch->Clear("C");

  LoadBranches();
  
  //Check if there is a centrality value, PbPb analysis, and if a centrality bin selection is requested
  //If we need a centrality bin, we select only those events in the corresponding bin.
  if( GetCentrality() && fCentralityBin[0] >= 0 && fCentralityBin[1] >= 0 )
  {
    Float_t cen = GetEventCentrality();
    if(cen > fCentralityBin[1] || cen < fCentralityBin[0]) return ; //reject events out of bin.
  }  
    
  // intermediate array with new clusters : init the array only once or clear from previous event
  if(!fCaloClusterArr) fCaloClusterArr    = new TObjArray(10000);
  else                 fCaloClusterArr->Delete();//Clear("C"); it leaks?

  InitGeometry(); // only once, must be done before OADB, geo OADB accessed here
  
  //Get the event, do some checks and settings
  CheckAndGetEvent() ;
  
  if (!fEvent) 
  {
    if(DebugLevel() > 0 ) printf("AliAnalysisTaksEMCALClusterize::UserExec() - Skip Event %d", (Int_t) Entry());
    return ;
  }

  //Init pointers, geometry, clusterizer, ocdb, aodb
  
  
  if(fAccessOCDB) AccessOCDB();
  if(fAccessOADB) AccessOADB(); // only once
  
  InitClusterization();
  
  //-------
  // Step 2
  
  // Make clusters
  if (fJustUnfold) ClusterUnfolding();
  else             ClusterizeCells() ;
  
  //-------
  // Step 3
  
  FillCaloClusterInEvent();
  
}      

