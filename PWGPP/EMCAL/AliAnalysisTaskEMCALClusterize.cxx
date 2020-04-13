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
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliDataFile.h"

// --- EMCAL
#include "AliEMCALAfterBurnerUF.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALClusterizerNxN.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALClusterizerv2.h"
#include "AliEMCALClusterizerv3.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALDigit.h"

#include "AliAnalysisTaskEMCALClusterize.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEMCALClusterize) ;
/// \endcond

//______________________________________________________________________________
// Constructor.
//______________________________________________________________________________
AliAnalysisTaskEMCALClusterize::AliAnalysisTaskEMCALClusterize(const char *name)
: AliAnalysisTaskSE(name)
, fEvent(0)
, fGeom(0),               fGeomName("") 
, fGeomMatrixSet(kFALSE), fLoadGeomMatrices(kFALSE)
, fOCDBpath(""),          fAccessOCDB(kFALSE)
, fDigitsArr(0),          fClusterArr(0)             
, fCaloClusterArr(0),     fCaloCells(0)
, fRecParam(0),           fClusterizer(0)
, fUnfolder(0),           fJustUnfold(kFALSE) 
, fOutputAODBranch(0),    fOutputAODBranchName("")   
, fOutputAODCells (0),    fOutputAODCellsName (""),   fInputCaloCellsName ("")    
, fOutputAODBranchSet(0)
, fFillAODFile(kFALSE),   fFillAODHeader(0)
, fFillAODCaloCells(0),   fRun(-1)
, fRecoUtils(0),          fConfigName("")
, fOrgClusterCellId()
, fCellLabels(),          fCellSecondLabels(),        fCellTime()
, fCellMatchdEta(),       fCellMatchdPhi()
, fRecalibrateWithClusterTime(0)
, fMaxEvent(0),           fMinEvent(0)
, fDoTrackMatching(0),    fUpdateCell(0)
, fSelectCell(kFALSE),    fSelectCellMinE(0),         fSelectCellMinFrac(0)
, fRejectBelowThreshold(kFALSE)
, fRemoveLEDEvents(kTRUE),fRemoveExoticEvents(kFALSE)
, fImportGeometryFromFile(kTRUE), fImportGeometryFilePath("")
, fOADBSet(kFALSE),       fAccessOADB(kTRUE),         fOADBFilePath("")
, fConstantTimeShift(0)
, fCentralityClass(""),   fUseAliCentrality(0),       fSelectEMCALEvent(0)
, fEMCALEnergyCut(0.),    fEMCALNcellsCut (0)
, fSetCellMCLabelFromCluster(0)
, fSetCellMCLabelFromEdepFrac(0)
, fRemapMCLabelForAODs(0)
, fInputFromFilter(0) 
, fTCardCorrEmulation(0), fTCardCorrClusEnerConserv(0)
, fRandom(0),             fRandomizeTCard(1)
, fTCardCorrMinAmp(0.01), fTCardCorrMinInduced(0) 
, fTCardCorrMaxInducedLowE(0), fTCardCorrMaxInduced(100)
, fPrintOnce(0)
, fDoMergedBCs(0x0)
, fLoad1DRecalibFactors(0)

{
  for(Int_t i = 0; i < 22;    i++)  
  {
    fGeomMatrix[i] = 0;
    fTCardCorrInduceEnerProb   [i] = 0;
    fTCardCorrInduceEnerFracMax[i] = 100;
    fTCardCorrInduceEnerFracMin[i] =-100;

    for(Int_t j = 0; j < 4 ; j++)
    {
      fTCardCorrInduceEner         [j][i] =  0 ;   
      fTCardCorrInduceEnerFrac     [j][i] =  0 ;   
      fTCardCorrInduceEnerFracP1   [j][i] =  0 ;   
      fTCardCorrInduceEnerFracWidth[j][i] =  0 ;   
    }
  }
  
  ResetArrays();
  
  fCentralityBin[0] = fCentralityBin[1]=-1;
}

//______________________________________________________________
/// Constructor.
//______________________________________________________________
AliAnalysisTaskEMCALClusterize::AliAnalysisTaskEMCALClusterize()
: AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskEMCALClusterize")
, fEvent(0)
, fGeom(0),                 fGeomName("") 
, fGeomMatrixSet(kFALSE),   fLoadGeomMatrices(kFALSE)
, fOCDBpath(""),            fAccessOCDB(kFALSE)
, fDigitsArr(0),            fClusterArr(0)             
, fCaloClusterArr(0),       fCaloCells(0)
, fRecParam(0),             fClusterizer(0)
, fUnfolder(0),             fJustUnfold(kFALSE) 
, fOutputAODBranch(0),      fOutputAODBranchName("")
, fOutputAODCells (0),      fOutputAODCellsName (""),  fInputCaloCellsName ("")  
, fOutputAODBranchSet(0)
, fFillAODFile(kFALSE),     fFillAODHeader(0)
, fFillAODCaloCells(0),     fRun(-1)
, fRecoUtils(0),            fConfigName("")
, fOrgClusterCellId()
, fCellLabels(),            fCellSecondLabels(),        fCellTime()
, fCellMatchdEta(),         fCellMatchdPhi()
, fRecalibrateWithClusterTime(0)
, fMaxEvent(0),             fMinEvent(0)             
, fDoTrackMatching(kFALSE), fUpdateCell(0)
, fSelectCell(kFALSE),      fSelectCellMinE(0),         fSelectCellMinFrac(0)
, fRejectBelowThreshold(kFALSE)
, fRemoveLEDEvents(kTRUE),  fRemoveExoticEvents(kFALSE)
, fImportGeometryFromFile(kTRUE), fImportGeometryFilePath("")
, fOADBSet(kFALSE),         fAccessOADB(kTRUE),        fOADBFilePath("")
, fConstantTimeShift(0)
, fCentralityClass(""),     fUseAliCentrality(0),      fSelectEMCALEvent(0)
, fEMCALEnergyCut(0.),      fEMCALNcellsCut (0)
, fSetCellMCLabelFromCluster(0)
, fSetCellMCLabelFromEdepFrac(0)
, fRemapMCLabelForAODs(0)
, fInputFromFilter(0)
, fTCardCorrEmulation(0),   fTCardCorrClusEnerConserv(0)
, fRandom(0),               fRandomizeTCard(1)
, fTCardCorrMinAmp(0.01),   fTCardCorrMinInduced(0)
, fTCardCorrMaxInducedLowE(0), fTCardCorrMaxInduced(100)
, fPrintOnce(0)
, fDoMergedBCs(0x0)
, fLoad1DRecalibFactors(0)
{
  for(Int_t i = 0; i < 22;    i++)  
  {
    fGeomMatrix[i] = 0;
    fTCardCorrInduceEnerProb   [i] = 0;
    fTCardCorrInduceEnerFracMax[i] = 100;
    fTCardCorrInduceEnerFracMin[i] =-100;

    for(Int_t j = 0; j < 4 ; j++)
    {
      fTCardCorrInduceEner         [j][i] =  0 ;   
      fTCardCorrInduceEnerFrac     [j][i] =  0 ;   
      fTCardCorrInduceEnerFracP1   [j][i] =  0 ;   
      fTCardCorrInduceEnerFracWidth[j][i] =  0 ;   
    }
  }
  
  ResetArrays();
  
  fCentralityBin[0] = fCentralityBin[1]=-1;
}

//_______________________________________________________________
/// Destructor.
//_______________________________________________________________
AliAnalysisTaskEMCALClusterize::~AliAnalysisTaskEMCALClusterize()
{
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
/// Reject cell if acceptance criteria not passed: 
///   * correct cell number
///   * is it bad channel 
///
/// \param absID: absolute cell ID number
/// \param badmap: consider the bad channel map, on by default, not needed for cross-Talk
///
/// \return bool quality of cell, exists or not 
//_______________________________________________________________________________
Bool_t AliAnalysisTaskEMCALClusterize::AcceptCell( Int_t absID, Bool_t badmap ) 
{  
  if ( absID < 0 || absID >= 24*48*fGeom->GetNumberOfSuperModules() ) 
    return kFALSE;
  
  Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1, status = 0; 
  if (!fGeom->GetCellIndex(absID,imod,iTower,iIphi,iIeta)) 
    return kFALSE; 
  
  // Do not include bad channels found in analysis,
  if ( badmap )
  {
    fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);  
    
    if ( fRecoUtils->IsBadChannelsRemovalSwitchedOn() && 
         fRecoUtils->GetEMCALChannelStatus(imod, ieta, iphi, status) ) 
      return kFALSE;
  }
  
  return kTRUE;
}

//_______________________________________________________
/// \return True if there is in the event an EMCal cluster with enough energy and with good quality.
/// Accept event given there is a EMCAL cluster with enough energy,
/// number of cells and not noisy, exotic.
//_______________________________________________________
Bool_t AliAnalysisTaskEMCALClusterize::AcceptEventEMCAL()
{
  if(!fSelectEMCALEvent)   return kTRUE; // accept
  
  if(fEMCALEnergyCut <= 0) return kTRUE; // accept
  
  Int_t           nCluster = fEvent -> GetNumberOfCaloClusters();
  Int_t           bc       = fEvent -> GetBunchCrossNumber();
  
  for(Int_t icalo = 0; icalo < nCluster; icalo++)
  {
    AliVCluster *clus = (AliVCluster*) (fEvent->GetCaloCluster(icalo));
    
    if( ( clus->IsEMCAL() ) && ( clus->GetNCells() > fEMCALNcellsCut ) && ( clus->E() > fEMCALEnergyCut ) &&
       fRecoUtils->IsGoodCluster(clus,fGeom,fCaloCells,bc))
    {
      
      AliDebug(1, Form("Accept :  E %2.2f > %2.2f, nCells %d > %d",
                       clus->E(), fEMCALEnergyCut, clus->GetNCells(), fEMCALNcellsCut));
      
      return kTRUE;
    }
    
  }// loop
  
  AliDebug(1,"Reject");
  
  return kFALSE;
}

//_______________________________________________
/// Set the AODB calibration, bad channels etc. parameters at least once
/// alignment matrices from OADB done in SetGeometryMatrices.
//_______________________________________________
void AliAnalysisTaskEMCALClusterize::AccessOADB()
{
  // Set it only once, unless run changed
  if ( fOADBSet ) return ; 
    
  TString pass      = GetPass();
  
  AliInfo(Form("Get AODB parameters from EMCAL in %s for run %d, and <%s>",fOADBFilePath.Data(),fRun,pass.Data()));
  
  Int_t nSM = fGeom->GetNumberOfSuperModules();
  
  // Bad map
  if(fRecoUtils->IsBadChannelsRemovalSwitchedOn())
  {
    AliOADBContainer *contBC=new AliOADBContainer("");
    if(fOADBFilePath!="")
      contBC->InitFromFile(Form("%s/EMCALBadChannels.root",fOADBFilePath.Data()),"AliEMCALBadChannels");
    else
      contBC->InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALBadChannels.root").data(),"AliEMCALBadChannels");

    TObjArray *arrayBC=(TObjArray*)contBC->GetObject(fRun);
    
    if(arrayBC)
    {
        fRecoUtils->SwitchOnDistToBadChannelRecalculation();
        AliInfo("Remove EMCAL bad cells");
        
        for (Int_t i=0; i < nSM; ++i)
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
    } else AliInfo("Do NOT remove EMCAL bad channels"); // run array
    
    delete contBC;
  }  // Remove bad
  
  // Energy Recalibration
  if(fRecoUtils->IsRecalibrationOn())
  {
    AliOADBContainer *contRF=new AliOADBContainer("");
    
    if(fOADBFilePath!="")
      contRF->InitFromFile(Form("%s/EMCALRecalib%s.root",fOADBFilePath.Data(), fLoad1DRecalibFactors ? "_1D" : ""),"AliEMCALRecalib");
    else
      contRF->InitFromFile(AliDataFile::GetFileNameOADB(Form("EMCAL/EMCALRecalib%s.root", fLoad1DRecalibFactors ? "_1D" : "")).data(),"AliEMCALRecalib");

    TObjArray *recal=(TObjArray*)contRF->GetObject(fRun); 
    
    if(recal)
    {
      TObjArray *recalpass=(TObjArray*)recal->FindObject(pass);
      
      if(recalpass)
      {
        TObjArray *recalib=(TObjArray*)recalpass->FindObject("Recalib");
        
        if(recalib)
        {
          AliInfo("Recalibrate EMCAL");

          if(fLoad1DRecalibFactors){
            TH1S *h = fRecoUtils->GetEMCALChannelRecalibrationFactors1D();
            if (h)
              delete h;
            h=(TH1S*)recalib->FindObject("EMCALRecalFactors");
              
            if (!h)
            {
              AliError("Can not get EMCALRecalFactors");
            }
            h->SetDirectory(0);
            fRecoUtils->SetEMCALChannelRecalibrationFactors1D(h);
          }else{
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
          }
        } else AliInfo("Do NOT recalibrate EMCAL, no params object array"); // array ok
      } else AliInfo("Do NOT recalibrate EMCAL, no params for pass"); // array pass ok
    } else AliInfo("Do NOT recalibrate EMCAL, no params for run");  // run number array ok
        
    delete contRF;
  } // Recalibration on
  
  // Energy Recalibration, apply on top of previous calibration factors
  if ( fRun > 200000 )
  {
    AliInfo(Form("Switch off Temperature corrections for Run %d (remember to remove when Run2 T corrections available!)",fRun)); 
    fRecoUtils->SwitchOffRunDepCorrection();
  }
  
  if ( fRecoUtils->IsRunDepRecalibrationOn() )
  {
    AliOADBContainer *contRFTD=new AliOADBContainer("");
    
    if(fOADBFilePath!="")
      contRFTD->InitFromFile(Form("%s/EMCALTemperatureCorrCalib.root",fOADBFilePath.Data()),"AliEMCALRunDepTempCalibCorrections");
    else
      contRFTD->InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTemperatureCorrCalib.root").data(),"AliEMCALRunDepTempCalibCorrections");

    TH1S *htd=(TH1S*)contRFTD->GetObject(fRun); 
    
    //If it did not exist for this run, get closes one
    if (!htd)
    {
      AliWarning(Form("No TemperatureCorrCalib Objects for run: %d",fRun));
      // let's get the closest fRun instead then..
      Int_t lower = 0;
      Int_t ic = 0;
      Int_t maxEntry = contRFTD->GetNumberOfEntries();
      
      while ( (ic < maxEntry) && (contRFTD->UpperLimit(ic) < fRun) ) {
        lower = ic;
        ic++;
      }
      
      Int_t closest = lower;
      if ( (ic<maxEntry) &&
          (contRFTD->LowerLimit(ic)-fRun) < (fRun - contRFTD->UpperLimit(lower)) ) {
        closest = ic;
      }
      
      AliWarning(Form("TemperatureCorrCalib Objects found closest id %d from run: %d", closest, contRFTD->LowerLimit(closest)));
      htd = (TH1S*) contRFTD->GetObjectByIndex(closest);
    } 

    if(htd)
    {
      AliInfo("Recalibrate (Temperature) EMCAL");
      
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
    } else AliInfo("Do NOT recalibrate EMCAL with T variations, no params TH1");
    
    delete contRFTD;
  } // Run by Run T calibration
  
  // Time Recalibration
  if(fRecoUtils->IsTimeRecalibrationOn())
  {
    AliOADBContainer *contTRF=new AliOADBContainer("");
    
    if(fOADBFilePath!="")
      contTRF->InitFromFile(Form("%s/EMCALTimeCalib.root",fOADBFilePath.Data()),"AliEMCALTimeCalib");
    else
      contTRF->InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTimeCalib.root").data(),"AliEMCALTimeCalib");

    TObjArray *trecal=(TObjArray*)contTRF->GetObject(fRun);
    
    if(trecal)
    {
      // pass number should be pass1 except on Run1 and special cases
      TString passM = pass;
      if ( pass=="spc_calo"   ) passM = "pass3";
      if ( fRun > 209121 ) passM = "pass1"; // run2 periods
      if ( pass == "muon_calo_pass1" && fRun > 209121 && fRun < 244284 )
        passM = "pass0";//period LHC15a-m

      TObjArray *trecalpass=(TObjArray*)trecal->FindObject(passM);

        if(trecalpass)
        {
          AliInfo("Time Recalibrate EMCAL");

          if(fDoMergedBCs){

            TH1S *h = (TH1S*)fRecoUtils->GetEMCALChannelTimeRecalibrationFactors(0);
            
            if (h)
              delete h;
          
            h = (TH1S*)trecalpass->FindObject("hAllTimeAv");// High Gain only
          
            if (!h) 
              AliError("Could not load hAllTimeAv");
            
            h->SetDirectory(0);

            fRecoUtils->SetEMCALChannelTimeRecalibrationFactors(0,h);

          }else{
            for (Int_t ibc = 0; ibc < 4; ++ibc) 
            {
              TH1F *h = (TH1F*)fRecoUtils->GetEMCALChannelTimeRecalibrationFactors(ibc);
          
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
          }
      } else AliInfo("Do NOT recalibrate time EMCAL, no params for pass"); // array pass ok
    } else AliInfo("Do NOT recalibrate time EMCAL, no params for run");  // run number array ok
    
    delete contTRF;
  } // Time recalibration on    
    
  // L1 Phase Time Recalibration
  if( fRecoUtils->IsL1PhaseInTimeRecalibrationOn() )
  {  
    // init default maps first
    if (!fRecoUtils->GetEMCALL1PhaseInTimeRecalibrationArray())
      fRecoUtils->InitEMCALL1PhaseInTimeRecalibration() ;
    
    AliOADBContainer *contBC = new AliOADBContainer("");
        
    TFile *timeFile;
    if(fOADBFilePath!="")
      timeFile = new TFile(Form("%s/EMCALTimeL1PhaseCalib.root",fOADBFilePath.Data()),"read");
    else
      timeFile = new TFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTimeL1PhaseCalib.root").data(),"read");

    if (!timeFile || timeFile->IsZombie())
    {
      AliFatal(Form("EMCALTimeL1PhaseCalib.root was not found in the path provided: %s",fOADBFilePath.Data()));
      return ;
    }
    
    if (timeFile) delete timeFile;
    
    if(fOADBFilePath!="")
      contBC->InitFromFile(Form("%s/EMCALTimeL1PhaseCalib.root",fOADBFilePath.Data()),"AliEMCALTimeL1PhaseCalib");
    else
      contBC->InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALTimeL1PhaseCalib.root").data(),"AliEMCALTimeL1PhaseCalib");

    TObjArray *arrayBC=(TObjArray*)contBC->GetObject(fRun);
    if (!arrayBC)
    {
      AliError(Form("No external L1 phase in time calibration set for run number: %d", fRun));
      fRecoUtils->SwitchOffL1PhaseInTimeRecalibration();
    }
    else
    {
      // Only 1 L1 phase correction possible, except special cases
      TString pass2 =  "pass1"; 
      
      if ( pass=="muon_calo_pass1" && fRun > 209121 && fRun < 244284 ) 
        pass2 = "pass0"; // period LHC15a-m

      TObjArray *arrayBCpass=(TObjArray*)arrayBC->FindObject(pass2);
      if (!arrayBCpass)
      {
        AliError(Form("No external L1 phase in time calibration set for: %d -%s", fRun,pass2.Data()));
        fRecoUtils->SwitchOffL1PhaseInTimeRecalibration();
      }
      else AliInfo("Recalibrate L1 Phase time");
            
      if(arrayBCpass)
      {
        if ( DebugLevel()>0 ) arrayBCpass->Print();
        
        TH1C *h = fRecoUtils->GetEMCALL1PhaseInTimeRecalibrationForAllSM(0);
        if (h) delete h;
        
        h = (TH1C*)arrayBCpass->FindObject(Form("h%d",fRun));
        
        if (!h) 
        {
          AliFatal(Form("There is no calibration histogram h%d for this run",fRun));
          return;
        }
        
        h->SetDirectory(0);
        fRecoUtils->SetEMCALL1PhaseInTimeRecalibrationForAllSM(h,0);

	//Now special case for PAR runs
	fRecoUtils->SwitchOffParRun();
	//access tree from OADB file
	TTree *tGID = (TTree*)arrayBCpass->FindObject(Form("h%d_GID",fRun));
	if(tGID){//check whether present = PAR run
	  fRecoUtils->SwitchOnParRun();
	  //access tree branch with PARs
	  ULong64_t parGlobalBCs;
	  tGID->SetBranchAddress("GID",&parGlobalBCs);
	  //set number of PARs in run
	  Short_t nPars = (Short_t) tGID->GetEntries();
	  fRecoUtils->SetNPars((Short_t)nPars);
	  //set global ID for each PAR
	  for (Short_t iParNumber = 0; iParNumber < nPars; ++iParNumber) {
	    tGID->GetEntry(iParNumber);
	    fRecoUtils->SetGlobalIDPar(parGlobalBCs,iParNumber);
	  }//loop over entries  
	  
	  //access GlobalID hiostograms for each PAR
	  for(Short_t iParNumber=1; iParNumber< fRecoUtils->GetNPars()+1;iParNumber++){
	    TH1C *hPar = (TH1C*)arrayBCpass->FindObject( Form("h%d_%llu",fRun,fRecoUtils->GetGlobalIDPar(iParNumber-1) ) );
	    if (!hPar) AliError( Form("Could not load h%d_%llu",fRun,fRecoUtils->GetGlobalIDPar(iParNumber-1) ) );
	    hPar->SetDirectory(0);
	    fRecoUtils->SetEMCALL1PhaseInTimeRecalibrationForAllSM(hPar,iParNumber);
	  }//loop over PARs
	}//end if tGID present  
      }
    }

    delete contBC;
  }   // L1 Phase Time Recalibration
    
  // Parameters already set once, so do not it again, unless run changes
  fOADBSet = kTRUE;
}  

//_________________________________________________
/// Access to OCDB stuff, avoid. 
/// Not sure it works anymore.
//_________________________________________________
void AliAnalysisTaskEMCALClusterize::AccessOCDB()
{
  // Set once per run
  if ( fOADBSet ) return;
  
  AliDebug(1,"Begin");
  
  AliCDBManager *cdb = AliCDBManager::Instance();
  
  if (fOCDBpath.Length())
  {
    cdb->SetDefaultStorage(fOCDBpath.Data());
    AliInfo(Form("Default storage %s",fOCDBpath.Data()));
  }
  
  cdb->SetRun(fRun);
  
  //
  // EMCAL from RAW OCDB
  if (fOCDBpath.Contains("alien:"))
  {
    cdb->SetSpecificStorage("EMCAL/Calib/Data","raw://");
    cdb->SetSpecificStorage("EMCAL/Calib/Time","raw://");
    cdb->SetSpecificStorage("EMCAL/Calib/Pedestals","raw://");
  }
  
  TString path = cdb->GetDefaultStorage()->GetBaseFolder();
  
  fOADBSet = kTRUE;
  
  return ;
}

//_____________________________________________________
/// Add to the digits the found induced energies in  MakeCellTCardCorrelation()
/// to new cells that before had no signal if new signal is larger than 10 MeV.
/// It is MC, but no MC label is assigned, just -2 to signal a new created energy
//_____________________________________________________
void AliAnalysisTaskEMCALClusterize::AddNewTCardInducedCellsToDigit() 
{
  // Get the number of stored digits, to assign the index of the new ones
  Int_t idigit = fDigitsArr->GetEntriesFast();
  
  Double_t fixTime =  615.*1e-9;
  
  for(Int_t j = 0; j < fgkNEMCalCells; j++)
  {
    // Newly created?
    if ( !fTCardCorrCellsNew[j] ) continue;
    
    // Accept only if at least 10 MeV
    if (  fTCardCorrCellsEner[j] < 0.01 ) continue; 
    
    // Check if it was not masked
    Float_t  amp  = 0;
    Double_t time = 0;
    if ( !fRecoUtils->AcceptCalibrateCell(j,0,amp,time,fCaloCells) ) continue;
    
//    printf("add new digit absId %d, accept? %d, digit %d, induced amp %2.2f\n",
//           j,fRecoUtils->AcceptCalibrateCell(j,0,amp,time,fCaloCells),idigit,fTCardCorrCellsEner[j]);
    
    // Now add the cell to the digits list
    //AliEMCALDigit* digit = 
    new((*fDigitsArr)[idigit]) AliEMCALDigit( -2, -2, j, fTCardCorrCellsEner[j], fixTime,AliEMCALDigit::kHG,idigit, 0, 0, 0);
    
    //digit->SetListOfParents(0,0x0,0x0); // not needed
        
    idigit++;
  }// loop on all possible cells
}


//_____________________________________________________
/// Get the input event, it can depend in embedded events what you want to get.
/// Also check if the quality of the event is good (has it EMCal clusters,
/// is the event triggered by LED or exotic ), if not reject it.
/// If we add more than clusters, add also Header and CaloCells here.
//_____________________________________________________
void AliAnalysisTaskEMCALClusterize::CheckAndGetEvent()
{
  fEvent = 0x0;
    
  Int_t eventN = Entry();
  
  AliAODInputHandler* aodIH = dynamic_cast<AliAODInputHandler*>
  ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  
  // Entry() does not work for AODs
  if ( eventN <= 0 && aodIH)
    eventN = aodIH->GetReadEntry(); 
  
  if ( eventN > fMaxEvent || eventN < fMinEvent ) 
    return ;
  
  //printf("AliAnalysisTaskEMCALClusterize::CheckAndGetEvent() - Event %d - Entry %d - (First,Last)=(%d,%d) \n", 
  //       eventN, (Int_t) Entry(), fMinEvent, fMaxEvent);
  
  // Check if input event are embedded events
  // If so, take output event
  if (aodIH && aodIH->GetMergeEvents())
  {
    fEvent  = AODEvent();
    
    if(!aodIH->GetMergeEMCALCells())
      AliFatal("Events merged but not EMCAL cells, check analysis settings!");
    
    AliDebug(1,"Use embedded events");
    
    AliDebug(1,Form("\t InputEvent  N Clusters %d, N Cells %d",
                    InputEvent()->GetNumberOfCaloClusters(),InputEvent()->GetEMCALCells()->GetNumberOfCells()));
    
    AliDebug(1,Form("\t MergedEvent  N Clusters %d, N Cells %d",
                    aodIH->GetEventToMerge()->GetNumberOfCaloClusters(), aodIH->GetEventToMerge()->GetEMCALCells()->GetNumberOfCells()));
    
//    if(DebugLevel() > 1)
//    {
//      for (Int_t icl=0; icl < aodIH->GetEventToMerge()->GetNumberOfCaloClusters(); icl++)
//      {
//        AliAODCaloCluster *sigCluster = aodIH->GetEventToMerge()->GetCaloCluster(icl);
//        if(sigCluster->IsEMCAL()) AliInfo(Form("\t \t Signal cluster: i %d, E  %f",icl,sigCluster->E()));
//      }
//    }
    
    AliDebug(1,Form("\t OutputEvent N Clusters %d, N Cells %d",
                    AODEvent()->GetNumberOfCaloClusters(), AODEvent()->GetEMCALCells()->GetNumberOfCells()));
  }
  else if(fInputFromFilter)
  {
    //printf("Get Input From Filtered AOD\n");
    fEvent =  AODEvent();
  }
  else 
  {
    fEvent =  InputEvent();
    if(fFillAODCaloCells) FillAODCaloCells();   
    if(fFillAODHeader)    FillAODHeader();
  }
  
  if (!fEvent) 
  {
    AliError("Event not available");
    return ;
  }

  //check global id for event in case of PAR in the run
  Short_t currentParIndex = 0;
  if(fRecoUtils->IsParRun()) {
    ULong64_t globalEventID = (ULong64_t)fEvent->GetBunchCrossNumber() + (ULong64_t)fEvent->GetOrbitNumber() * (ULong64_t)3564 + (ULong64_t)fEvent->GetPeriodNumber() * (ULong64_t)59793994260;
    for(Short_t ipar=0;ipar< fRecoUtils->GetNPars();ipar++){
      if(globalEventID >= fRecoUtils->GetGlobalIDPar(ipar)) {
	currentParIndex++;
      }
    }
  }
  fRecoUtils->SetCurrentParNumber(currentParIndex);
  
  //Recover the pointer to CaloCells container
  if ( fInputCaloCellsName.Length() == 0 ) 
    fCaloCells = fEvent->GetEMCALCells();
  else      
  {
    fCaloCells = (AliVCaloCells*) fEvent->FindListObject(fInputCaloCellsName);
    if ( !fCaloCells ) 
      AliWarning(Form("CaloCells branch <%s> not found use STD!",fInputCaloCellsName.Data()));
    else 
      fCaloCells = fEvent->GetEMCALCells();
  }
  
  //Process events if there is a high energy cluster
  if(!AcceptEventEMCAL())  { fEvent = 0x0 ; return ; }
  
  //-------------------------------------------------------------------------------------
  // Reject events if LED was firing, use only for LHC11a data 
  // Reject event if triggered by exotic cell and remove exotic cells if not triggered
  //-------------------------------------------------------------------------------------
  
  if( IsLEDEvent( fRun ) ) { fEvent = 0x0 ; return ; }
  
  if( IsExoticEvent() )    { fEvent = 0x0 ; return ; }
  
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
    // Create array of clusters/cells and put it in the input event, if output AOD not selected, only once
    InputEvent()->AddObject(fOutputAODBranch);
    AliInfo(Form("Add AOD clusters branch <%s> to input event",fOutputAODBranchName.Data()));
    
    if ( fOutputAODBranchName.Length() > 0 )
    {
      InputEvent()->AddObject(fOutputAODCells);
      AliInfo(Form("Add AOD cells branch <%s> to input event",fOutputAODCellsName.Data()));
    }
    
    fOutputAODBranchSet = kTRUE;
  }
}

//____________________________________________________
/// Recluster calorimeter cells, transform them into digits,
/// feed the clusterizer with them and get new list of clusters.
/// In case of MC, first loop on the clusters and fill MC label to array.
/// Filter the cells not being exotic, bad and recalibrate them before clusterizing.
//____________________________________________________
void AliAnalysisTaskEMCALClusterize::ClusterizeCells()
{
  Int_t nClusters     = fEvent->GetNumberOfCaloClusters();
  Int_t nClustersOrg  = 0;
  
  AliAODInputHandler* aodIH = dynamic_cast<AliAODInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if(aodIH && aodIH->GetEventToMerge())  //Embedding
    nClusters = aodIH->GetEventToMerge()->GetNumberOfCaloClusters(); //Get clusters directly from embedded signal

  ResetArrays();

  // Loop on original clusters, get MC labels, cluster time (OLD AODs), 
  // or track matching residuals (if matching is not requested)
  if ( fSetCellMCLabelFromEdepFrac || fSetCellMCLabelFromCluster || 
       fRecalibrateWithClusterTime || !fDoTrackMatching)
  {
    for (Int_t i = 0; i < nClusters; i++)
    {
      AliVCluster *clus = 0;
      if(aodIH && aodIH->GetEventToMerge()) //Embedding
        clus = aodIH->GetEventToMerge()->GetCaloCluster(i); //Get clusters directly from embedded signal
      else      
        clus = fEvent->GetCaloCluster(i);
      
      if ( !clus ) continue;
            
      if ( !clus || !clus->IsEMCAL() ) continue;
      
      nClustersOrg++;

      Int_t label = clus->GetLabel();
      Int_t label2 = -1 ;
      if (clus->GetNLabels() >=2 ) label2 = clus->GetLabelAt(1) ;
      
      AliDebug(2, Form("recover original cluster %d info: Id %d, E %2.3f, N cells %d, TOF %3.2f, N labels %d, label %d;",
                       i,clus->GetID(), clus->E(),clus->GetNCells(),clus->GetTOF()*1e9, clus->GetNLabels(), label) );

      if ( fSetCellMCLabelFromEdepFrac && fDebug > 1 )
      {
        for(UInt_t imc = 0; imc < clus->GetNLabels(); imc++) 
        {
          printf("\t mc %d) Label %d, E dep frac %1.3f; ",
                 imc, clus->GetLabelAt(imc),clus->GetClusterMCEdepFraction(imc));
        }
        printf("\n");
      }
      
      UShort_t * index    = clus->GetCellsAbsId() ;
      for(Int_t icell=0; icell < clus->GetNCells(); icell++ )
      {
        fOrgClusterCellId[index[icell]] = i;
        fCellTime[index[icell]]         = clus->GetTOF();
        fCellMatchdEta[index[icell]]    = clus->GetTrackDz();
        fCellMatchdPhi[index[icell]]    = clus->GetTrackDx();
        
        if ( !fSetCellMCLabelFromEdepFrac )
        {
          fCellLabels[index[icell]]       = label;
          fCellSecondLabels[index[icell]] = label2;
        }
        
        AliDebug(2, Form("\t : cell %d Id %d, clus %d, time %2.3e, MatchEta %2.3f, MatchPhi %2.3f; 1st label %d, 2nd label %d",
                         icell,index[icell],fOrgClusterCellId[index[icell]],
                         fCellTime[index[icell]], fCellMatchdEta[index[icell]] ,  fCellMatchdPhi[index[icell]],
                         fCellLabels[index[icell]],fCellSecondLabels[index[icell]] ) );
      } // cell in cluster loop
    } // cluster loop 
    AliDebug(2, Form("N original cluster %d",nClustersOrg) );
  } // use cluster org info
  
  // Do here induced cell energy assignation by T-Card correlation emulation, ONLY MC
  if(fTCardCorrEmulation) MakeCellTCardCorrelation();

  //  
  // Transform CaloCells into Digits
  //
  Int_t    idigit =  0;
  Int_t    id     = -1;
  Float_t  amp    = -1; 
  Double_t time   = -1; 
  
  Int_t bc = InputEvent()->GetBunchCrossNumber();

  for (Int_t icell = 0; icell < fCaloCells->GetNumberOfCells(); icell++)
  {
    // Get cell values, recalibrate and not include bad channels found in analysis, nor cells with too low energy, nor exotic cell
    id = fCaloCells->GetCellNumber(icell);
    Bool_t accept = fRecoUtils->AcceptCalibrateCell(id,bc,amp,time,fCaloCells);
    time-=fConstantTimeShift*1e-9; // only in case of simulations done before 2015

    // Do not include cells with too low energy, nor exotic cell
    // Comment out since it removes some cells that could be accepted by the clusterizer, not clear why.
    // To get inline with what is done in the EMCal correction framework
//    if( amp  < fRecParam->GetMinECut() ||
//        time > fRecParam->GetTimeMax() ||
//        time < fRecParam->GetTimeMin()    ) accept = kFALSE;
    
    // In case of old AOD analysis cell time is -1 s, approximate replacing by time of the cluster the digit belongs.
    if (fRecalibrateWithClusterTime)
    { 
      time = fCellTime[id];
      //printf("cell %d time org %f - ",id, time*1.e9);
      fRecoUtils->RecalibrateCellTime(id,bc,time);
      //printf("recal %f\n",time*1.e9);
    }
    
    //Exotic?
    if (accept && fRecoUtils->IsExoticCell(id,fCaloCells,bc))
        accept = kFALSE;
    
    if( !accept )
    {
      AliDebug(2,Form("Remove channel absId %d, index %d of %d, amp %f, time %f",
               id,icell, fCaloCells->GetNumberOfCells(), amp, time*1.e9));
      continue;
    }
    
    //
    // MC
    //
    
    // Old way
    Int_t mcLabel = fCaloCells->GetMCLabel(icell);
    Float_t eDep = amp;
    
    //printf("--- Selected cell %d--- amp %f\n",id,amp);

    // New way
    TArrayI labeArr(0);
    TArrayF eDepArr(0);
    Int_t nLabels = 0;

    if ( !fSetCellMCLabelFromEdepFrac )
    {
      // Old way to recover/set the cell MC label
      // Only possibility for old Run1 productions

      if     ( fSetCellMCLabelFromCluster == 1 ) mcLabel = fCellLabels[id]; // Older aliroot MC productions
     
      else if( fSetCellMCLabelFromCluster == 0 && 
               fRemapMCLabelForAODs)             RemapMCLabelForAODs(mcLabel);
      
      else                                       mcLabel = -1; // found later
      
      // Last parameter of the digit object should be MC deposited energy, 
      // since it is not available in aliroot before year 2016, add just the cell amplitude so that
      // we give more weight to the MC label of the cell with highest energy in the cluster
            
      Float_t efrac = fCaloCells->GetEFraction(icell);
      
      // When checking the MC of digits, give weight to cells with embedded signal
      if (mcLabel > 0 && efrac < 1.e-6) efrac = 1;
      
      eDep *= efrac ; 
    }
    else if ( fOrgClusterCellId[id] >= 0 ) // fSetCellMCLabelFromEdepFrac = true
    {
      // New way, valid only for MC productions with aliroot > v5-07-21
      mcLabel = -1;
      
      // Map the digit to cell index for later to calculate the cell MC energy deposition map
      fCellLabels[id] = idigit; 

      AliVCluster *clus = 0;
      Int_t iclus = fOrgClusterCellId[id];
      
      AliDebug(1, Form("EdepFrac use for : absId %d, idigit %d, iclus %d, amp %2.3f",id,idigit,iclus,amp) );
      
      if(aodIH && aodIH->GetEventToMerge()) //Embedding
        clus = aodIH->GetEventToMerge()->GetCaloCluster(iclus); //Get clusters directly from embedded signal
      else      
        clus = fEvent->GetCaloCluster(iclus);
      
      fRecoUtils->RecalculateCellLabelsRemoveAddedGenerator(id, clus, MCEvent(), amp, labeArr, eDepArr);
      nLabels = labeArr.GetSize();
      
      AliDebug(1, Form("N labels after EdepFrac info use: %d, amp %2.3f (check changes)",nLabels,amp) );

    } // cell MC label, new
    
    // Apply here the found induced energies
    if ( fTCardCorrEmulation )
    {
//      if( TMath::Abs(fTCardCorrCellsEner[id]) > 0.001 )
//        printf("add energy to digit %d, absId %d: amp %2.2f + %2.2f\n",idigit,id,amp,fTCardCorrCellsEner[id]);
      amp+=fTCardCorrCellsEner[id];
    }

    //
    // Create the digit
    //
    if(amp <= 0.01) continue ; // accept if > 10 MeV
    
    AliDebug(5,Form("*** Add digit *** digit %d, AbsId %d, amp %2.3f, time %3.2f, nlabels %d, label %d",
                    idigit, id, amp, time*1.e9,nLabels,mcLabel));
    
    AliEMCALDigit* digit = new((*fDigitsArr)[idigit]) 
      AliEMCALDigit( mcLabel, mcLabel, id, amp, time,AliEMCALDigit::kHG,idigit, 0, 0, eDep);
    
    if(nLabels > 0)
      digit->SetListOfParents(nLabels,labeArr.GetArray(),eDepArr.GetArray());
    
    idigit++;
    
  }
  
  fDigitsArr->Sort();
  
  // Call after Sort, if not it screws up digits index order in clusterization
  if ( fTCardCorrEmulation ) AddNewTCardInducedCellsToDigit();


  //-------------------------------------------------------------------------------------
  // Do the clusterization
  //-------------------------------------------------------------------------------------        
  
  fClusterizer->Digits2Clusters("");
  
  //-------------------------------------------------------------------------------------
  // Transform the recpoints into AliVClusters
  //-------------------------------------------------------------------------------------

  RecPoints2Clusters();
  
  if(!fCaloClusterArr)
  { 
    AliWarning(Form("No array with CaloClusters, input RecPoints entries %d",fClusterArr->GetEntriesFast()));
    return;    
  }
  
  AliDebug(1,Form("N clusters: before recluster %d, after recluster %d, recpoints %d",
                  nClustersOrg, fCaloClusterArr->GetEntriesFast(),fClusterArr->GetEntriesFast()));
  
//    if(fCaloClusterArr->GetEntriesFast() != fClusterArr->GetEntriesFast())
//    {
//      AliInfo("\t Some RecRoints not transformed into CaloClusters (clusterizer %d, unfold %d): Input entries %d - Output entries %d - %d (not fast)\n",
//             fRecParam->GetClusterizerFlag(),fRecParam->GetUnfold(),
//             fClusterArr->GetEntriesFast(), fCaloClusterArr->GetEntriesFast(), fCaloClusterArr->GetEntries());
//    }
}


//_____________________________________________________
/// Take the event clusters and unfold them.
//_____________________________________________________
void AliAnalysisTaskEMCALClusterize::ClusterUnfolding()
{
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
        fRecoUtils->RecalibrateClusterEnergy(fGeom, clus, fCaloCells);
        
        //CalibrateCells
        for (Int_t icell = 0; icell < fCaloCells->GetNumberOfCells(); icell++)
        {
          if (fCaloCells->GetCell(icell, cellNumber, cellAmplitude, cellTime, cellMCLabel, cellEFrac) != kTRUE)
            break;
          
          Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1, status = 0; 
          fGeom->GetCellIndex(cellNumber,imod,iTower,iIphi,iIeta); 
          fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);	
          
          //Do not include bad channels found in analysis?
          if( fRecoUtils->IsBadChannelsRemovalSwitchedOn() && 
              fRecoUtils->GetEMCALChannelStatus(imod, ieta, iphi, status))
            continue;
          
          fCaloCells->SetCell(icell, cellNumber, cellAmplitude*fRecoUtils->GetEMCALChannelRecalibrationFactor(imod,ieta,iphi),cellTime);
          
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
        AliWarning("Wrong CaloCluster type?");
      
      nClustersOrg++;
    }
  }
  
  //Do the unfolding
  fUnfolder->UnfoldClusters(fCaloClusterArr, fCaloCells);
  
  //CLEAN-UP
  fUnfolder->Clear();
}

//_______________________________________________________________
/// Configure fRecoUtils with some standard arguments for common analysis configurations
///
/// The input parameters:
/// \param reco: pointer to object to initialize in this macro.
/// \param bMC: Bool, indicates if data is MC.
/// \param bExotic: Bool, indicates if exotic clusters are removed.
/// \param bNonLin: Bool, indicates if non linearity correction is applied on clusters.
/// \param bRecalE: Bool, indicates if energy recalibration is applied.
/// \param bBad: Bool, indicates if bad channels/clusters are removed.
/// \param bRecalT: Bool, indicates if time is calibrated.
/// \param debug: int debug level, print info on settings in the macro
///
//_______________________________________________________________
void AliAnalysisTaskEMCALClusterize::ConfigureEMCALRecoUtils
(Bool_t  bMC    , Bool_t  bExotic, Bool_t  bNonLin,  
 Bool_t  bRecalE, Bool_t  bBad   , Bool_t  bRecalT, Int_t   debug)
{
  if ( debug > 0 ) printf("ConfigureEMCALRecoUtils() - **** Start ***\n");
  
  // Init
  if(!fRecoUtils) fRecoUtils = new AliEMCALRecoUtils ;
  
  // Exotic cells removal
  
  if(bExotic)
  {
    if ( debug > 0 ) printf("ConfigureEMCALRecoUtils() - Remove exotics in EMCAL\n");
    fRecoUtils->SwitchOnRejectExoticCell() ;
    fRecoUtils->SwitchOnRejectExoticCluster(); 
    
    //  fRecoUtils->SetExoticCellDiffTimeCut(50);     // If |t cell max - t cell in cross| > 50 do not add its energy, avoid 
    fRecoUtils->SetExoticCellFractionCut(0.97);   // 1-Ecross/Ecell > 0.97 -> out
    fRecoUtils->SetExoticCellMinAmplitudeCut(4.); // 4 GeV    
  }  
  
  // Recalibration factors
  
  if(bRecalE && ! bMC)
  {
    if ( debug > 0 ) printf("ConfigureEMCALRecoUtils() - Switch on energy recalibration in EMCAL\n");
    fRecoUtils->SwitchOnRecalibration();
    fRecoUtils->SwitchOnRunDepCorrection();    
  } 
  
  // Use one histogram for all BCs
  if (fDoMergedBCs)
    fRecoUtils->SetUseOneHistForAllBCs(fDoMergedBCs);

  // Remove EMCAL hot channels 
  
  if(bBad)
  {
    if ( debug > 0 ) printf("ConfigureEMCALRecoUtils() - Switch on bad channels removal in EMCAL\n");
    fRecoUtils->SwitchOnBadChannelsRemoval();
    fRecoUtils->SwitchOnDistToBadChannelRecalculation();
  }
  
  // *** Time recalibration settings ***
  
  if(bRecalT && ! bMC)
  {
    if ( debug > 0 ) printf("ConfigureEMCALRecoUtils() - Switch on time recalibration in EMCAL\n");
    fRecoUtils->SwitchOnTimeRecalibration();
    fRecoUtils->SwitchOnL1PhaseInTimeRecalibration() ;
  }
  
  // Recalculate position with method
  
  fRecoUtils->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);   
  
  // Non linearity
  
  if( bNonLin ) 
  { 
    if(!bMC)
    {
      if ( debug > 0 ) printf("ConfigureEMCALRecoUtils() xxx SET Non linearity correction kBeamTestCorrected xxx\n");
      fRecoUtils->SetNonLinearityFunction(AliEMCALRecoUtils::kBeamTestCorrectedv3);
    }
    else
    {       
      if ( debug > 0 ) printf("ConfigureEMCALRecoUtils() xxx SET Non linearity correction kPi0MCv3 xxx\n");
      fRecoUtils->SetNonLinearityFunction(AliEMCALRecoUtils::kPi0MCv3);
    }
  }
  else 
  {
    if ( debug > 0 ) printf("ConfigureEMCALRecoUtils() xxx DON'T SET Non linearity correction xxx\n");
    fRecoUtils->SetNonLinearityFunction(AliEMCALRecoUtils::kNoCorrection);
  }
  
}

//_____________________________________________________
/// Put calo cells in standard branch.
//_____________________________________________________
void AliAnalysisTaskEMCALClusterize::FillAODCaloCells()
{
  AliVCaloCells &eventEMcells = *(fEvent->GetEMCALCells());
  Int_t nEMcell = eventEMcells.GetNumberOfCells() ;
  
  AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
  aodEMcells.CreateContainer(nEMcell);
  aodEMcells.SetType(AliVCaloCells::kEMCALCell);
  Double_t calibFactor = 1.;   
  Int_t    status = 0;
  for (Int_t iCell = 0; iCell < nEMcell; iCell++) 
  { 
    Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
    fGeom->GetCellIndex(eventEMcells.GetCellNumber(iCell),imod,iTower,iIphi,iIeta); 
    fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);	
    
    if(fRecoUtils->IsRecalibrationOn())
    { 
      calibFactor = fRecoUtils->GetEMCALChannelRecalibrationFactor(imod,ieta,iphi);
    }
    
    if(!fRecoUtils->GetEMCALChannelStatus(imod, ieta, iphi, status))
    { 
      // Channel is not declared as bad
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
/// Put event header information in standard AOD branch.
//__________________________________________________
void AliAnalysisTaskEMCALClusterize::FillAODHeader()
{
  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*> (fEvent);
  AliAODEvent* aodevent = dynamic_cast<AliAODEvent*> (fEvent);
  
  Double_t pos[3]   ;
  Double_t covVtx[6];
  for (Int_t i = 0; i < 6; i++)  covVtx[i] = 0.;
  
  AliAODHeader* header = dynamic_cast<AliAODHeader*>(AODEvent()->GetHeader());
  if(!header) AliFatal("Not a standard AOD");
  header->SetRunNumber(fRun);
  
  if(esdevent)
  {
    TTree* tree = fInputHandler->GetTree();
    if (tree) 
    {
      TFile* file = tree->GetCurrentFile();
      if (file) header->SetESDFileName(file->GetName());
    }
  }
  else if (aodevent) {
    AliAODHeader * aodheader = dynamic_cast<AliAODHeader*>(aodevent->GetHeader());
    if(!aodheader) AliFatal("Not a standard AOD");
    header->SetESDFileName(aodheader->GetESDFileName());
  }
  
  header->SetBunchCrossNumber(fEvent->GetBunchCrossNumber());
  header->SetOrbitNumber(fEvent->GetOrbitNumber());
  header->SetPeriodNumber(fEvent->GetPeriodNumber());
  header->SetEventType(fEvent->GetEventType());
  
  // Centrality
  if(fUseAliCentrality)
  {
    if(fEvent->GetCentrality())
      header->SetCentrality(new AliCentrality(*(fEvent->GetCentrality())));
    else
      header->SetCentrality(0);
  }
  
  // Trigger
  header->SetOfflineTrigger(fInputHandler->IsEventSelected()); // propagate the decision of the physics selection

  header->SetFiredTriggerClasses(fEvent->GetFiredTriggerClasses());
  
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
  
  Float_t diamxy[2]={(Float_t)fEvent->GetDiamondX(),(Float_t)fEvent->GetDiamondY()};
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
/// Get the CaloClusters array, do some final calculations
/// and put the clusters in the output or input event
/// as a separate branch.
//___________________________________________________________
void AliAnalysisTaskEMCALClusterize::FillCaloClusterInEvent()
{
  // First recalculate track-matching for the new clusters
  if(fDoTrackMatching) 
  {
    fRecoUtils->FindMatches(fEvent,fCaloClusterArr,fGeom,MCEvent());
  }
    
  // Put the new clusters in the AOD list
  
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
        AliDebug(2,Form("Matched Track index %d to new cluster %d",trackIndex,i));
      }
      
      Float_t dR = 999., dZ = 999.;
      fRecoUtils->GetMatchedResiduals(newCluster->GetID(),dZ,dR);
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
    //
    //printf("New cluster %d) ID = %d, E %2.2f, Time  %3.0f,  N Cells %d, N MC labels %d, main MC label %d, all MC labels:\n", 
    //       i, newCluster->GetID(), newCluster->E(), newCluster->GetTOF()*1e9, newCluster->GetNCells(), newCluster->GetNLabels(), newCluster->GetLabel());
    //
    //for(Int_t imc = 0; imc < newCluster->GetNLabels(); imc++) 
    //  printf("%d) Label %d, E dep frac %0.2f; ",
    //         imc,newCluster->GetLabelAt(imc),newCluster->GetClusterMCEdepFraction(imc));
    //if(newCluster->GetNLabels() > 0) printf("\n");
    
    // Calculate distance to bad channel for new cluster. Make sure you give the list of bad channels.
    fRecoUtils->RecalculateClusterDistanceToBadChannel(fGeom, fCaloCells, newCluster);
    
    new((*fOutputAODBranch)[i])  AliAODCaloCluster(*newCluster);
    
    AliDebug(2,Form("New cluster %d of %d, energy %f, mc label %d",
                    newCluster->GetID(), kNumberOfCaloClusters, newCluster->E(), newCluster->GetLabel()));
    
  } // cluster loop
  
  fOutputAODBranch->Expand(kNumberOfCaloClusters); // resize TObjArray to 'remove' slots
}


//________________________________________________________________
/// Get centrality/multiplicity percentile
//________________________________________________________________
Float_t AliAnalysisTaskEMCALClusterize::GetEventCentrality() const                          
{ 
  if(fUseAliCentrality)
  {
    if(GetCentrality()) return GetCentrality()->GetCentralityPercentile(fCentralityClass) ;
    else                return -1.       ; 
  }
  else
  {
    if(GetMultSelCen()) return GetMultSelCen()->GetMultiplicityPercentile(fCentralityClass,kTRUE) ;
    else                return -1.       ; 
  }
  
  return -1.;
}

//_______________________________________________
/// Get or guess pass number/string from path of filename.
//_______________________________________________
TString AliAnalysisTaskEMCALClusterize::GetPass()
{
  if (!AliAnalysisManager::GetAnalysisManager()->GetTree()) 
  {
    AliError("AliAnalysisTaskEMCALClusterize::GetPass() - Pointer to tree = 0, returning null");
    return TString("");
  }
  
  if (!AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()) 
  {
    AliError("AliAnalysisTaskEMCALClusterize::GetPass() - Null pointer input file, returning null");
    return TString("");
  }
  
  TString pass(AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()->GetName());
  if      (pass.Contains("ass1")) return TString("pass1");
  else if (pass.Contains("ass2")) return TString("pass2");
  else if (pass.Contains("ass3")) return TString("pass3");
  else if (pass.Contains("ass4")) return TString("pass4");
  else if (pass.Contains("ass5")) return TString("pass5");
  else if (pass.Contains("LHC11c") && pass.Contains("spc_calo") ) return TString("spc_calo");
  else if (pass.Contains("calo") || pass.Contains("high_lumi"))
  {
    AliInfo("Path contains <calo> or <high-lumi>, set as <pass1>");
    return TString("pass1");
  }
  else if (pass.Contains("LHC14a1a"))
  {
    AliInfo("Enable EMCal energy calibration for this MC production!!");
    
    fRecoUtils->SwitchOnRecalibration();

    return TString("LHC14a1a");
  }
  
  // No condition fullfilled, give a default value
  AliInfo("Pass number string not found");
    
  return TString("");            
}

//_________________________________________
/// Init analysis with configuration macro if available.
/// Init other parameters, pointers if not done before with default settings.
//_________________________________________
void AliAnalysisTaskEMCALClusterize::Init()
{
  if(fDebug >=0) (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(this->ClassName(),fDebug);

  fOADBSet           = kFALSE;
  
  fBranchNames       = "ESD:AliESDHeader.,EMCALCells.";
  
  if(!fRecParam)     fRecParam  = new AliEMCALRecParam;
  if(!fRecoUtils)    fRecoUtils = new AliEMCALRecoUtils();  

  // Use one histogram for all BCs
  if (fDoMergedBCs)
    fRecoUtils->SetUseOneHistForAllBCs(fDoMergedBCs);
  
  if(fMaxEvent          <= 0) fMaxEvent          = 1000000000;
  if(fSelectCellMinE    <= 0) fSelectCellMinE    = 0.005;     
  if(fSelectCellMinFrac <= 0) fSelectCellMinFrac = 0.001;
  fRejectBelowThreshold = kFALSE;

  // Centrality
  if(fCentralityClass  == "") fCentralityClass  = "V0M";
  
  if (fOCDBpath            == "") fOCDBpath            = "raw://" ;
  if (fOutputAODBranchName == "") fOutputAODBranchName = "newEMCALClusters" ;
  
  if(gROOT->LoadMacro(fConfigName) >=0)
  {
    AliInfo(Form("Configure analysis with %s",fConfigName.Data()));
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
    fMinEvent         = clus->fMinEvent;
    fDoTrackMatching  = clus->fDoTrackMatching;
    fUpdateCell       = clus->fUpdateCell;
    fOutputAODBranchName = clus->fOutputAODBranchName;
    for(Int_t i = 0; i < 22; i++) fGeomMatrix[i] = clus->fGeomMatrix[i] ;
    fCentralityClass  = clus->fCentralityClass;
    fCentralityBin[0] = clus->fCentralityBin[0];
    fCentralityBin[1] = clus->fCentralityBin[1];
    fUseAliCentrality = clus->fUseAliCentrality;
  }
}  

//_______________________________________________________
/// Select clusterization/unfolding algorithm and
/// set all the needed parameters.
//_______________________________________________________
void AliAnalysisTaskEMCALClusterize::InitClusterization()
{
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
    fClusterizer = new AliEMCALClusterizerv1 (fGeom);
  else if(fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv2) 
    fClusterizer = new AliEMCALClusterizerv2(fGeom);
  else if(fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerv3) 
    fClusterizer = new AliEMCALClusterizerv3(fGeom);
  else if(fRecParam->GetClusterizerFlag() == AliEMCALRecParam::kClusterizerNxN)
  { 
    fClusterizer = new AliEMCALClusterizerNxN(fGeom);
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
  fClusterizer->SetTimeCalibration       ( fRecParam->IsTimeCalibrationOn()    );  
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
   
    fClusterizer->SetRejectBelowThreshold(fRejectBelowThreshold);//here we set option of unfolding: split or reject energy
    fClusterizer->InitClusterUnfolding();
    
  }// to unfold
}

//________________________________________________________________
/// Init geometry and set the geometry matrix,
/// for the first event, skip the rest.
/// Even if run number changes, geom only changes from year to year so first is enough.
//________________________________________________________________
void AliAnalysisTaskEMCALClusterize::InitGeometry()
{
  if(fGeomMatrixSet) return;
  
  if (!fGeom)
  {
    if(fGeomName=="")
    {
      fGeom = AliEMCALGeometry::GetInstanceFromRunNumber(fRun);
      AliInfo(Form("Get EMCAL geometry name <%s> for run %d",fGeom->GetName(),fRun));
    }
    else
    {
      fGeom = AliEMCALGeometry::GetInstance(fGeomName);
      AliInfo(Form("Set EMCAL geometry name to <%s>",fGeomName.Data()));
    }
    
    // Init geometry, I do not like much to do it like this ...
    if(fImportGeometryFromFile && !gGeoManager)
    {
      if(fImportGeometryFilePath=="") // If not specified, set location depending on run number
      {
        // "$ALICE_ROOT/EVE/alice-data/default_geo.root"
        if      (fRun <  140000) fImportGeometryFilePath = AliDataFile::GetFileNameOADB("EMCAL/geometry_2010.root").data();
        else if (fRun <  171000) fImportGeometryFilePath = AliDataFile::GetFileNameOADB("EMCAL/geometry_2011.root").data();
        else if (fRun <  198000) fImportGeometryFilePath = AliDataFile::GetFileNameOADB("EMCAL/geometry_2012.root").data(); // 2012-2013
        else                     fImportGeometryFilePath = AliDataFile::GetFileNameOADB("EMCAL/geometry_2015.root").data(); // >=2015
      }
      
      AliInfo(Form("Import %s",fImportGeometryFilePath.Data()));
      
      TGeoManager::Import(fImportGeometryFilePath) ;
    }

    AliDebug(1,Form("Init for run=%d",fRun));
    if (!gGeoManager) AliDebug(1,"Careful!, gGeoManager not loaded, load misalign matrices");
	} // geometry pointer did not exist before
  
  if(fLoadGeomMatrices)
  {
    AliInfo("Load user defined EMCAL geometry matrices");
    
    // OADB if available
    AliOADBContainer emcGeoMat("AliEMCALgeo");
    if(fOADBFilePath!="")
      emcGeoMat.InitFromFile(Form("%s/EMCALlocal2master.root",fOADBFilePath.Data()),"AliEMCALgeo");
    else
      emcGeoMat.InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALlocal2master.root").data(),"AliEMCALgeo");

    TObjArray *matEMCAL=(TObjArray*)emcGeoMat.GetObject(fRun,"EmcalMatrices");
    
    for(Int_t mod=0; mod < (fGeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
    {
      
      if (!fGeomMatrix[mod]) // Get it from OADB
      {
        AliDebug(2,Form("EMCAL matrices SM %d, %p",mod,((TGeoHMatrix*) matEMCAL->At(mod))));
        //((TGeoHMatrix*) matEMCAL->At(mod))->Print();
        
        fGeomMatrix[mod] = (TGeoHMatrix*) matEMCAL->At(mod) ;
      }        
      
      if(fGeomMatrix[mod])
      {
        if(DebugLevel() > 1) 
          fGeomMatrix[mod]->Print();
        
        fGeom->SetMisalMatrix(fGeomMatrix[mod],mod) ;  
      }
      else if(gGeoManager)
      {
        AliWarning(Form("Set matrix for SM %d from gGeoManager",mod));
        fGeom->SetMisalMatrix(fGeom->GetMatrixForSuperModuleFromGeoManager(mod),mod) ;
      }
      else
      {
        AliError(Form("Alignment atrix for SM %d is not available",mod));
      }
      
      fGeomMatrixSet=kTRUE;
      
    }//SM loop
  }//Load matrices
  else if(!gGeoManager)
  {
    AliInfo("Get geo matrices from data");
    //Still not implemented in AOD, just a workaround to be able to work at least with ESDs	
    if(!strcmp(fEvent->GetName(),"AliAODEvent")) 
    {
      AliWarning("Use ideal geometry, values geometry matrix not kept in AODs");
    }//AOD
    else 
    {	
      AliDebug(1,"Load Misaligned matrices");
      
      AliESDEvent* esd = dynamic_cast<AliESDEvent*>(fEvent) ;
      
      if(!esd)
      {
        AliError("This event does not contain ESDs?");
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
      
    } // ESD
  } // Load matrices from Data
}

//____________________________________________________
/// Check if event is exotic, get an exotic cell and
/// compare with triggered patch
/// If there is a match remove event ... to be completed,
/// filled with something provisional
//____________________________________________________
Bool_t AliAnalysisTaskEMCALClusterize::IsExoticEvent()
{
  if(!fRemoveExoticEvents) return kFALSE;
  
  // Loop on cells
    
  Float_t totCellE = 0;
  Int_t bc = InputEvent()->GetBunchCrossNumber();
    
  for(Int_t icell = 0; icell < fCaloCells->GetNumberOfCells(); icell++)
  {
    Float_t  ecell = 0 ;
    Double_t tcell = 0 ;
    
    Int_t absID   = fCaloCells->GetCellNumber(icell);
    Bool_t accept = fRecoUtils->AcceptCalibrateCell(absID,bc,ecell,tcell,fCaloCells);
    tcell-=fConstantTimeShift*1e-9;// Only for MC simulations done before 2015

    if(accept && !fRecoUtils->IsExoticCell(absID,fCaloCells,bc)) totCellE += ecell;
  }
  
  //  TString triggerclasses = event->GetFiredTriggerClasses();
  //    printf("AliAnalysisTaskEMCALClusterize - reject event %d with cluster  - reject event with ncells in SM3 %d and SM4 %d\n",(Int_t)Entry(),ncellsSM3, ncellsSM4);
  //    if(fFillAODFile) AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kFALSE);;
  //    return;
  //  
  
  //printf("TotE cell %f\n",totCellE);
  if(totCellE < 1) return kTRUE;
  
  return kFALSE;
} 

//________________________________________________________________
/// Check if event is LED, is so remove it. Affected LHC11a runs.
//________________________________________________________________
Bool_t AliAnalysisTaskEMCALClusterize::IsLEDEvent(const Int_t run)
{
  if(!fRemoveLEDEvents) return kFALSE;
  
  //check events of LHC11a period
  if(run < 146858 || run > 146860) return kFALSE ;
  
  // Count number of cells with energy larger than 0.1 in SM3, cut on this number
  Int_t ncellsSM3 = 0;
  for(Int_t icell = 0; icell < fCaloCells->GetNumberOfCells(); icell++)
  {
    if ( fCaloCells->GetAmplitude (icell) > 0.1 && 
         fCaloCells->GetCellNumber(icell)/(24*48)==3 ) ncellsSM3++;
  }
  
  TString triggerclasses = fEvent->GetFiredTriggerClasses();
  
  Int_t ncellcut = 21;
  if(triggerclasses.Contains("EMC")) ncellcut = 35;
  
  if( ncellsSM3 >= ncellcut)
  {
    AliInfo(Form("Reject event %d with ncells in SM3 %d",(Int_t)Entry(),ncellsSM3));
    if(fFillAODFile) AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kFALSE);;
    return kTRUE;
  }
  
  return kFALSE;
}

//_______________________________________________________
/// Calculate the induced energy in a cell belonging to the
/// same T-Card as the reference cell.
/// Used in MakeCellTCardCorrelation()
/// \param absId Id number of cell in same T-Card as reference cell
/// \param absIdRef Id number of reference cell
/// \param sm Supermodule number of cell 
/// \param ampRef Amplitude of the reference cell
/// \param cellCase Type of cell with respect reference cell 0: up or down, 1: up or down on the diagonal, 2: left or right, 3: 2nd row up/down both left/right
//_______________________________________________________
void AliAnalysisTaskEMCALClusterize::CalculateInducedEnergyInTCardCell
(Int_t absId, Int_t absIdRef, Int_t sm, Float_t ampRef, Int_t cellCase) 
{
  // Check that the cell exists
  if( !AcceptCell(absId,0) ) return ; 
  
  // Get the fraction
  Float_t frac = fTCardCorrInduceEnerFrac[cellCase][sm] + ampRef * fTCardCorrInduceEnerFracP1[cellCase][sm];

  // Use an absolute minimum and maximum fraction if calculated one is out of range
  if ( frac < fTCardCorrInduceEnerFracMin[sm] ) frac = fTCardCorrInduceEnerFracMin[sm];
  if ( frac > fTCardCorrInduceEnerFracMax[sm] ) frac = fTCardCorrInduceEnerFracMax[sm];   
  
  AliDebug(1,Form("\t fraction %2.3f",frac));
  
  // Randomize the induced fraction, if requested
  if ( fRandomizeTCard )
  {
    frac = fRandom.Gaus(frac, fTCardCorrInduceEnerFracWidth[cellCase][sm]);
    
    AliDebug(1,Form("\t randomized fraction %2.3f",frac));
  }
  
  // If too small or negative, do nothing else
  if ( frac < 0.0001 ) return;
  
  // Calculate induced energy
  Float_t inducedE = fTCardCorrInduceEner[cellCase][sm] + ampRef * frac;

  // Check if we induce too much energy, in such case use a constant value
  if ( fTCardCorrMaxInduced < inducedE ) inducedE = fTCardCorrMaxInduced;
  
  AliDebug(1,Form("\t induced E %2.3f",inducedE));
  
  // Add the induced energy, check if cell existed
  // Check that the induced+amp is large enough to avoid extra linearity effects
  // typically of the order of the clusterization cell energy cut
  // But if it is below 1 ADC, typically 10 MeV, also do it, to match Beam test linearity
  Float_t amp = fCaloCells->GetCellAmplitude(absId) ;
  if ( (amp+inducedE) > fTCardCorrMinInduced || inducedE < fTCardCorrMaxInducedLowE )
  {
    fTCardCorrCellsEner[absId] += inducedE;
    
    // If original energy of cell was null, create new one 
    if ( amp < 0.01 ) fTCardCorrCellsNew[absId] = kTRUE;
  }
  else return ;
  
  // Subtract the added energy to main cell, if energy conservation is requested
  if ( fTCardCorrClusEnerConserv )
     fTCardCorrCellsEner[absIdRef] -= inducedE;
}

//_______________________________________________________
/// Recover each cell amplitude and absId and induce energy 
/// in cells around of the same T-Card, up to second row
//_______________________________________________________
void AliAnalysisTaskEMCALClusterize::MakeCellTCardCorrelation()
{
  Int_t    id     = -1;
  Float_t  amp    = -1; 
  
  // Loop on all cells with signal
  for (Int_t icell = 0; icell < fCaloCells->GetNumberOfCells(); icell++)
  {
    id  = fCaloCells->GetCellNumber(icell);
    amp = fCaloCells->GetAmplitude (icell); // fCaloCells->GetCellAmplitude(id);
    
    if ( amp <= fTCardCorrMinAmp ) continue ;
    
    //
    // First get the SM, col-row of this tower
    Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
    fGeom->GetCellIndex(id,imod,iTower,iIphi,iIeta); 
    fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);  

    //
    // Determine randomly if we want to create a correlation for this cell, 
    // depending the SM number of the cell
    if ( fTCardCorrInduceEnerProb[imod] < 1 )
    {  
      Float_t rand = fRandom.Uniform(0, 1);
      
      if ( rand > fTCardCorrInduceEnerProb[imod] )
      {
        AliDebug(1,Form("Do not difuse E of cell %d, sm %d, amp %2.2f: SM fraction %2.2f > %2.2f",
                        id,imod,amp,fTCardCorrInduceEnerProb[imod],rand));
        continue;
      }
    }
    
    AliDebug(1,Form("Reference cell absId %d, iEta %d, iPhi %d, sm %d, amp %2.2f",id,ieta,iphi,imod,amp));

    //
    // Get the absId of the cells in the cross and same T-Card
    Int_t absIDup    = -1;
    Int_t absIDdo    = -1;
    Int_t absIDlr    = -1;
    Int_t absIDuplr  = -1;
    Int_t absIDdolr  = -1;
    Int_t absIDup2   = -1;
    Int_t absIDup2lr = -1;
    Int_t absIDdo2   = -1;
    Int_t absIDdo2lr = -1;
    
    // Only 2 columns in the T-Card, +1 for even and -1 for odd with respect reference cell
    Int_t colShift = 0;
    if (  (ieta%2) && ieta <= AliEMCALGeoParams::fgkEMCALCols-1 ) colShift = -1; 
    if ( !(ieta%2) && ieta >= 0 )                                 colShift = +1;               
    
    absIDlr = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta+colShift); 
    
    // Check up / down cells from reference cell not out of SM and in same T-Card
    if (  iphi < AliEMCALGeoParams::fgkEMCALRows-1 ) 
    {
      absIDup   = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta);
      absIDuplr = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta+colShift); 
    }
    
    if (  iphi > 0 ) 
    {
      absIDdo   = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta);
      absIDdolr = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta+colShift); 
    }
    
    // Check 2 up / 2 down cells from reference cell not out of SM
    if (  iphi < AliEMCALGeoParams::fgkEMCALRows-2 ) 
    {
      absIDup2   = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+2, ieta);
      absIDup2lr = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+2, ieta+colShift); 
    }
    
    if (  iphi > 1 )   
    {
      absIDdo2   = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-2, ieta);
      absIDdo2lr = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-2, ieta+colShift); 
    }
    
    // In same T-Card?
    if ( TMath::FloorNint(iphi/8) != TMath::FloorNint((iphi+1)/8) ) { absIDup  = -1 ; absIDuplr  = -1 ; }
    if ( TMath::FloorNint(iphi/8) != TMath::FloorNint((iphi-1)/8) ) { absIDdo  = -1 ; absIDdolr  = -1 ; }
    if ( TMath::FloorNint(iphi/8) != TMath::FloorNint((iphi+2)/8) ) { absIDup2 = -1 ; absIDup2lr = -1 ; }
    if ( TMath::FloorNint(iphi/8) != TMath::FloorNint((iphi-2)/8) ) { absIDdo2 = -1 ; absIDdo2lr = -1 ; }
    
    // Calculate induced energy to T-Card cells
    
    AliDebug(1,Form("cell up %d:"  ,absIDup));
    CalculateInducedEnergyInTCardCell(absIDup   , id, imod, amp, 0);    
    AliDebug(1,Form("cell down %d:",absIDdo));
    CalculateInducedEnergyInTCardCell(absIDdo   , id, imod, amp, 0);
 
    AliDebug(1,Form("cell up left-right %d:"  ,absIDuplr));
    CalculateInducedEnergyInTCardCell(absIDuplr , id, imod, amp, 1);    
    AliDebug(1,Form("cell down left-right %d:",absIDdolr));
    CalculateInducedEnergyInTCardCell(absIDdolr , id, imod, amp, 1);

    AliDebug(1,Form("cell left-right %d:",absIDlr));
    CalculateInducedEnergyInTCardCell(absIDlr   , id, imod, amp, 2);

    AliDebug(1,Form("cell up 2nd row %d:"  ,absIDup2));
    CalculateInducedEnergyInTCardCell(absIDup2  , id, imod, amp, 3);    
    AliDebug(1,Form("cell down 2nd row %d:",absIDdo2));
    CalculateInducedEnergyInTCardCell(absIDdo2  , id, imod, amp, 3);
    
    AliDebug(1,Form("cell up left-right 2nd row %d:"  ,absIDup2lr));
    CalculateInducedEnergyInTCardCell(absIDup2lr, id, imod, amp, 3);    
    AliDebug(1,Form("cell down left-right 2nd row %d:",absIDdo2lr));
    CalculateInducedEnergyInTCardCell(absIDdo2lr, id, imod, amp, 3);
    
  } // cell loop
  
}

//_______________________________________________________
/// Print clusterization task parameters.
//_______________________________________________________
void AliAnalysisTaskEMCALClusterize::PrintParam()
{
  AliInfo(Form("Geometry: name <%s>, matrix set <%d>, load matrix <%d>, import geo <%d> from path <%s>",
               fGeomName.Data(), fGeomMatrixSet, fLoadGeomMatrices, fImportGeometryFromFile, fImportGeometryFilePath.Data()));
  
  if ( fAccessOCDB ) AliInfo(Form("OCDB path name <%s>", fOCDBpath.Data()));
  if ( fAccessOADB ) AliInfo(Form("OADB path name <%s>", fOADBFilePath.Data()));
 
  if ( fInputCaloCellsName.Length() > 0 ) 
    AliInfo(Form("Input CaloCells <%s>", fInputCaloCellsName.Data()));
  
  AliInfo(Form("Just Unfold clusters <%d>, new clusters list name <%s>, new cells name <%s>", 
               fJustUnfold, fOutputAODBranchName.Data(), fOutputAODCellsName.Data()));
  
  if ( fFillAODFile ) AliInfo(Form("Fill new AOD file with: header <%d>, cells <%d>",fFillAODHeader,fFillAODCaloCells));
  
  AliInfo(Form("Use cell time for cluster <%d>, Apply constant time shift <%2.2f>, Do track-matching <%d>, Update cells <%d>, Input from ESD filter <%d>",
               fRecalibrateWithClusterTime, fConstantTimeShift, fDoTrackMatching, fUpdateCell, fInputFromFilter));
  
  AliInfo(Form("Reject events out of range: %d < N event < %d, LED <%d>, exotics <%d>", 
               fMinEvent, fMaxEvent, fRemoveLEDEvents, fRemoveExoticEvents));
  
  if (fCentralityBin[0] != -1 && fCentralityBin[1] != -1 ) 
    AliInfo(Form("Centrality bin [%2.2f,%2.2f], class <%s>, use AliCentrality? <%d>", 
                 fCentralityBin[0], fCentralityBin[1], fCentralityClass.Data(), fUseAliCentrality));
  
  if ( fSelectEMCALEvent ) 
    AliInfo(Form("Select events with signal in EMCal: E min <%2.2f>, n cell min <%d>", fEMCALEnergyCut, fEMCALNcellsCut));
  
  AliInfo(Form("MC label from cluster <%d>, Use EdepFrac <%d>, remap AODs <%d>",
               fSetCellMCLabelFromCluster, fSetCellMCLabelFromEdepFrac, fRemapMCLabelForAODs));
}

//_______________________________________________________
/// Print parameters for T-Card correlation emulation.
//_______________________________________________________
void AliAnalysisTaskEMCALClusterize::PrintTCardParam()
{
  if(!fTCardCorrEmulation)
  {
    AliInfo("T-Card emulation not activated");
    return;
  }
  
  AliInfo(Form("T-Card emulation activated, energy conservation <%d>, randomize E <%d>, induced energy parameters:",
               fTCardCorrClusEnerConserv,fRandomizeTCard));
  
  AliInfo(Form("T-Card emulation super-modules fraction: Min cell E %2.1f MeV; induced Min E %2.1f MeV; Max at low E %2.1f MeV; Max E %2.2f GeV",
               fTCardCorrMinAmp*1000,fTCardCorrMinInduced*1000,fTCardCorrMaxInducedLowE*1000,fTCardCorrMaxInduced));
  
  for(Int_t ism = 0; ism < 22; ism++)
  {
    printf("\t sm %d, fraction %2.3f, E frac abs min %2.3e max %2.3e \n",
           ism, fTCardCorrInduceEnerProb[ism],fTCardCorrInduceEnerFracMin[ism],fTCardCorrInduceEnerFracMax[ism]);
    
    for(Int_t icell = 0; icell < 4; icell++)
    {
      printf("\t \t cell type %d, c %2.4e, p0 %2.4e, p1 %2.4e, sigma %2.4e \n",
             icell,fTCardCorrInduceEner[icell][ism],fTCardCorrInduceEnerFrac[icell][ism],
             fTCardCorrInduceEnerFracP1[icell][ism],fTCardCorrInduceEnerFracWidth[icell][ism]);     
    }
  }
}

//_______________________________________________________
/// Restore clusters from recPoints.
/// Cluster energy, global position, cells and their amplitude
/// fractions are restored.
//_______________________________________________________
void AliAnalysisTaskEMCALClusterize::RecPoints2Clusters()
{
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
      
      if ( !fRecParam->GetUnfold() && (ratios[ncellsTrue] > 1 || ratios[ncellsTrue] < 1)  ) 
        AliWarning(Form("recpoint cell E %2.3f but digit E %2.3f and no unfolding", recPoint->GetEnergiesList()[c], digit->GetAmplitude()));
      
      // In case of unfolding, remove digits with unfolded energy too low      
      if(fSelectCell)
      {
        if     (recPoint->GetEnergiesList()[c] < fSelectCellMinE || ratios[ncellsTrue] < fSelectCellMinFrac)  
        {
         
          AliDebug(2,Form("Too small energy in cell of cluster: cluster cell %f, digit %f",
                          recPoint->GetEnergiesList()[c],digit->GetAmplitude()));
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
        AliDebug(2,Form("Skipping cluster with no cells avobe threshold E = %f, ncells %d",
                        recPoint->GetEnergy(), ncells));
      continue;
    }
    
    //if(ncellsTrue != ncells) printf("Old E %f, ncells %d; New E %f, ncells %d\n",recPoint->GetEnergy(),ncells,clusterE,ncellsTrue);
    
    if(clusterE <  fRecParam->GetClusteringThreshold()) 
    {
      AliDebug(2,Form("Remove cluster with energy below seed threshold %f",clusterE));
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
      AliDebug(2,Form("Cells removed from cluster (ncells %d, ncellsTrue %d), recalculate Shower Shape",ncells,ncellsTrue));
      
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
      //
      // Normal case, trust what the clusterizer has found
      //
      Int_t    parentMult = 0;
      Int_t   *parentList   = recPoint->GetParents(parentMult);
      clus->SetLabel(parentList, parentMult);
      
      //
      // Set the cell energy deposition fraction map:
      //
      if( parentMult > 0 && fSetCellMCLabelFromEdepFrac )
      {
        Float_t *parentListDE = recPoint->GetParentsDE();         // deposited energy
        clus->SetClusterMCEdepFractionFromEdepArray(parentListDE);
        
        UInt_t * mcEdepFracPerCell = new UInt_t[ncellsTrue];
        
        // Get the digit that originated this cell cluster
//        AliVCaloCells* cells = 0x0; 
//        if (aodIH && aodIH->GetMergeEvents()) cells = AODEvent()  ->GetEMCALCells();
//        else                                  cells = InputEvent()->GetEMCALCells();
        
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
              mcEdepFracPerCell[icell] = clus->PackMCEdepFraction(mcEDepFrac);
            }
          } // at least one parent label in digit
        } // cell in cluster loop
        
        clus->SetCellsMCEdepFractionMap(mcEdepFracPerCell);
        
        delete [] mcEdepFracPerCell;
        
      } // at least one parent in cluster, do the cell primary packing
    } /// Set the MC labels, normal procedure in reconstruction
  } // recPoints loop
}

//_____________________________________________________________________
/// MC label for Cells not remapped after ESD filtering,
/// it happened in old productions, do it here.
//_____________________________________________________________________
void AliAnalysisTaskEMCALClusterize::RemapMCLabelForAODs(Int_t & label)
{
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
/// Reset arrays containing information for all possible cells.
//________________________________________________
void AliAnalysisTaskEMCALClusterize::ResetArrays()
{
  for(Int_t j = 0; j < fgkNEMCalCells; j++)
  {
    fOrgClusterCellId[j] = -1;
    fCellLabels[j]       = -1;
    fCellSecondLabels[j] = -1;
    fCellTime[j]         =  0.;
    fCellMatchdEta[j]    = -999;
    fCellMatchdPhi[j]    = -999;
    
    fTCardCorrCellsEner[j] = 0.; 
    fTCardCorrCellsNew [j] = kFALSE; 
   }
}

//_____________________________________________________________________________________________________
/// Set the cluster MC label, the digizer was filled
/// with the most likely MC label for all cells in original cluster.
/// Now check the second most likely MC label and add it to the new cluster.
//_____________________________________________________________________________________________________
void AliAnalysisTaskEMCALClusterize::SetClustersMCLabelFrom2SelectedLabels(AliEMCALRecPoint  * recPoint,
                                                                           AliAODCaloCluster * clus)
{
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
    } // positive cell number
  }
}

//___________________________________________________________________________________________________
/// Get the original clusters that contribute to the new cluster, assign the labels of such clusters
/// to the new cluster. Only approximatedly valid when input and output are V1 clusters, or input are V2
/// clusters and output are any other type of clusters. Handle with care.
//___________________________________________________________________________________________________
void AliAnalysisTaskEMCALClusterize::SetClustersMCLabelFromOriginalClusters(AliAODCaloCluster * clus)
{
  TArrayI clArray(300) ; //Weird if more than a few clusters are in the origin ...
  clArray.Reset();
  Int_t nClu = 0;
  Int_t nLabTotOrg = 0;
  Float_t emax = -1;
  Int_t idMax = -1;
  
  AliVEvent * event = fEvent;
  
  // In case of embedding the MC signal is in the added event
  AliAODInputHandler* aodIH = dynamic_cast<AliAODInputHandler*>((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if(aodIH && aodIH->GetEventToMerge())  //Embedding
    event = aodIH->GetEventToMerge(); //Get clusters directly from embedded signal

  
  // Find the clusters that originally had the cells
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
    Int_t idCluster = (Int_t) clArray.GetAt(iLoopCluster);
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
// Init geometry, create list of output clusters and cells.
//____________________________________________________________
void AliAnalysisTaskEMCALClusterize::UserCreateOutputObjects()
{
  // Clusters
  fOutputAODBranch = new TClonesArray("AliAODCaloCluster", 0);

  if(fOutputAODBranchName.Length()==0)
  {
    fOutputAODBranchName = "newEMCALClustersArray";
    AliInfo("Cluster branch name not set, set it to newEMCALClustersArray");
  }
  
  fOutputAODBranch->SetName(fOutputAODBranchName);
  
  if( fFillAODFile )
  {
    //fOutputAODBranch = new TClonesArray("AliAODCaloCluster", 0);
    
    //fOutputAODBranch->SetOwner(kFALSE);
    
    AddAODBranch("TClonesArray", &fOutputAODBranch);
  }
  
  // Cells
  if ( fOutputAODBranchName.Length() > 0 )
  {
    fOutputAODCells = new AliAODCaloCells(fOutputAODCellsName,fOutputAODCellsName,AliAODCaloCells::kEMCALCell); 
  
    if( fFillAODFile ) AddAODBranch("AliAODCaloCells", &fOutputAODCells);
  }
}

//_______________________________________________________________________
/// Update or create CaloCells container if calibration or some changes were applied.
/// Delete previouly existing content in the container.
//________________________________________________________________________
void AliAnalysisTaskEMCALClusterize::UpdateCells()
{  
  // Update cells only in case re-calibration was done 
  // or bad map applied or additional T-Card cells added.
  if(!fRecoUtils->IsBadChannelsRemovalSwitchedOn() && 
     !fRecoUtils->IsRecalibrationOn()              && 
     !fRecoUtils->IsRunDepRecalibrationOn()        && 
     !fRecoUtils->IsTimeRecalibrationOn()          && 
     !fRecoUtils->IsL1PhaseInTimeRecalibrationOn() &&
     !fTCardCorrEmulation                            ) return;
  
  const Int_t   ncells = fCaloCells->GetNumberOfCells();
  const Int_t   ndigis = fDigitsArr->GetEntries();

  if ( fOutputAODCellsName.Length() > 0 ) 
  {
    fOutputAODCells->DeleteContainer();
    fOutputAODCells->CreateContainer(ndigis);
  }
  else if ( ncells != ndigis ) // update case 
  {
    fCaloCells->DeleteContainer();
    fCaloCells->CreateContainer(ndigis);
  }
  
  for (Int_t idigit = 0; idigit < ndigis; ++idigit) 
  {
    AliEMCALDigit *digit = static_cast<AliEMCALDigit*>(fDigitsArr->At(idigit));
    
    Double_t cellAmplitude  = digit->GetAmplitude();
    Short_t  cellNumber     = digit->GetId();
    Double_t cellTime       = digit->GetTime();
    
    Bool_t highGain = kFALSE;
    if( digit->GetType() == AliEMCALDigit::kHG ) highGain = kTRUE;
    
    // Only for MC
    // Get the label of the primary particle that generated the cell
    // Assign the particle that deposited more energy
    Int_t   nparents  = digit->GetNiparent();
    Int_t   cellMcEDepFrac =-1 ;
    Float_t cellMcLabel    =-1.;
    if ( nparents > 0 )
    {
      for ( Int_t jndex = 0 ; jndex < nparents ; jndex++ ) 
      { 
        if(cellMcEDepFrac >= digit->GetDEParent(jndex+1)) continue ;
        
        cellMcLabel   = digit->GetIparent (jndex+1);
        cellMcEDepFrac= digit->GetDEParent(jndex+1);          
      } // all primaries in digit      
    } // select primary label
    
    if ( cellMcEDepFrac < 0 ) cellMcEDepFrac = 0.;
      
    if ( fUpdateCell ) 
      fCaloCells     ->SetCell(idigit, cellNumber, cellAmplitude, cellTime, cellMcLabel, cellMcEDepFrac, highGain);
    else
      fOutputAODCells->SetCell(idigit, cellNumber, cellAmplitude, cellTime, cellMcLabel, cellMcEDepFrac, highGain);
  }
  
   if ( ncells != ndigis )
   {
     if ( fUpdateCell ) 
       fCaloCells     ->Sort();
     else
       fOutputAODCells->Sort();
   }
}

//_______________________________________________________
/// Do clusterization event by event, execute different steps
///  * 1) Do some checks on the kind of events (ESD, AOD) or if some filtering is needed, initializations
///    load things and clear arrays
///  * 2) Clusterize
///     * a) just unfolding existing clusters (fJustUnfold)
///     * b) recluster cells
///        * +  convert cells into digits (calibrating them)
///        * +  recluster digits into recPoints with the chosen clusterizer (V1, V1+Unfold,V2, NxN) with methods in AliEMCALClusterizer
///        * +  transform recPoints into CaloClusters
///  * 3) Do some final calculations in the found clusters (track-matching) and put them in an AOD branch
//_______________________________________________________
void AliAnalysisTaskEMCALClusterize::UserExec(Option_t *)
{
  //-------
  // Step 1

  // Remove the contents of AOD branch output list set in the previous event
  fOutputAODBranch->Clear("C");
  
  LoadBranches();
  
  // Check if there is a centrality value, PbPb analysis, and if a centrality bin selection is requested
  // If we need a centrality bin, we select only those events in the corresponding bin.
  if( fCentralityBin[0] >= 0 && fCentralityBin[1] >= 0 )
  {
    Float_t cen = GetEventCentrality();
    if(cen > fCentralityBin[1] || cen < fCentralityBin[0]) return ; //reject events out of bin.
  }  
    
  // Intermediate array with new clusters : init the array only once or clear from previous event
  if(!fCaloClusterArr) fCaloClusterArr    = new TObjArray(10000);
  else                 fCaloClusterArr->Delete();//Clear("C"); it leaks?

  
  // In case of analysis in multiple runs, check the OADB again
  if ( InputEvent()->GetRunNumber() != fRun )  
  {
    fRun = InputEvent()->GetRunNumber();
    
    fOADBSet = kFALSE; // recover the OADB for this run
    
    AliInfo(Form("Set run to %d",fRun));
  }
  
  InitGeometry(); // only once, must be done before OADB, geo OADB accessed here

  //access OADB, before first event because of PAR runs
  if(fAccessOADB) AccessOADB(); // only once

  // Get the event, do some checks and settings
  CheckAndGetEvent() ;
  
  if (!fEvent) 
  {
    AliDebug(1,Form("Skip Event %d", (Int_t) Entry()));
    return ;
  }

  // Init pointers, geometry, clusterizer, ocdb, aodb
  
  if(fAccessOCDB) AccessOCDB();
  //if(fAccessOADB) AccessOADB(); // only once
  
  InitClusterization();
  
  // Print once the analysis parameters
  if ( fDebug > 0 || !fPrintOnce )
  {
    //fRecParam->Print("reco"); // AliInfo not printed ...
    
    PrintParam();
    
    PrintTCardParam();
    
    fPrintOnce = kTRUE;
  }
  
  //-------
  // Step 2
  
  // Make clusters
  if (fJustUnfold) ClusterUnfolding();
  else             ClusterizeCells() ;
  
  //-------
  // Step 3
  
  FillCaloClusterInEvent();
  
  if ( fUpdateCell || fOutputAODCellsName.Length() > 0 ) 
    UpdateCells();
}


