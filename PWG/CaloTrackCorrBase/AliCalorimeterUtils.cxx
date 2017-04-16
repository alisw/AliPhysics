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
#include <TGeoManager.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <TFile.h>
#include <TParticle.h>
#include <AliMCEvent.h>

// --- ANALYSIS system ---
#include "AliCalorimeterUtils.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAODPWG4Particle.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliOADBContainer.h"
#include "AliAnalysisManager.h"
#include "AliAODMCParticle.h"
#include "AliLog.h"

// --- Detector ---
#include "AliEMCALGeometry.h"
#include "AliPHOSGeoUtils.h"

#include "AliFiducialCut.h" // Needed for detector flag enum kEMCAL, kPHOS

/// \cond CLASSIMP
ClassImp(AliCalorimeterUtils) ;
/// \endcond

  
//____________________________________________
/// Constructor. Initialize parameters.
//____________________________________________
AliCalorimeterUtils::AliCalorimeterUtils() : 
TObject(), fDebug(0), 
fEMCALGeoName(""),
fPHOSGeoName (""), 
fEMCALGeo(0x0),                   fPHOSGeo(0x0), 
fEMCALGeoMatrixSet(kFALSE),       fPHOSGeoMatrixSet(kFALSE), 
fLoadEMCALMatrices(kFALSE),       fLoadPHOSMatrices(kFALSE),
fRemoveBadChannels(kFALSE),       fPHOSBadChannelMap(0x0), 
fNCellsFromPHOSBorder(0),
fNMaskCellColumns(0),             fMaskCellColumns(0x0),
fRecalibration(kFALSE),           fRunDependentCorrection(kFALSE),
fPHOSRecalibrationFactors(),      fEMCALRecoUtils(new AliEMCALRecoUtils),
fRecalculatePosition(kFALSE),     fCorrectELinearity(kFALSE),
fRecalculateMatching(kFALSE),
fCutR(20),                        fCutZ(20),
fCutEta(20),                      fCutPhi(20),
fLocMaxCutE(0),                   fLocMaxCutEDiff(0),
fPlotCluster(0),                  fOADBSet(kFALSE),
fOADBForEMCAL(kFALSE),            fOADBForPHOS(kFALSE),
fOADBFilePathEMCAL(""),           fOADBFilePathPHOS(""),
fImportGeometryFromFile(0),       fImportGeometryFilePath(""),
fNSuperModulesUsed(0),            fRunNumber(0),
fMCECellClusFracCorrOn(0),        fMCECellClusFracCorrParam()
{
  InitParameters();
  for(Int_t i = 0; i < 22; i++) fEMCALMatrix[i] = 0 ;
  for(Int_t i = 0; i < 5 ; i++) fPHOSMatrix [i] = 0 ;
}

//_________________________________________
// Destructor.
//_________________________________________
AliCalorimeterUtils::~AliCalorimeterUtils()
{
  //if(fPHOSGeo)  delete fPHOSGeo  ;
  if(fEMCALGeo) delete fEMCALGeo ;
  
  if(fPHOSBadChannelMap) 
{ 
    fPHOSBadChannelMap->Clear();
    delete  fPHOSBadChannelMap;
  }
	
  if(fPHOSRecalibrationFactors) 
{ 
    fPHOSRecalibrationFactors->Clear();
    delete  fPHOSRecalibrationFactors;
  }
	
  if(fEMCALRecoUtils)   delete fEMCALRecoUtils ;
  if(fNMaskCellColumns) delete [] fMaskCellColumns;
}

//____________________________________________________
/// Set the AODB calibration, bad channels etc. parameters 
/// at least once.
//____________________________________________________
void AliCalorimeterUtils::AccessOADB(AliVEvent* event)
{
  // Set it only once
  if(fOADBSet) return ; 
  
  if(fRunNumber <= 0) fRunNumber = event->GetRunNumber() ; // take the run number from the event itself
  TString pass      = GetPass();
  
  // EMCAL
  if(fOADBForEMCAL)
  {
    AliInfo(Form("Get AODB parameters from EMCAL in %s for run %d, and <%s>",fOADBFilePathEMCAL.Data(),fRunNumber,pass.Data()));
    
    Int_t nSM = fEMCALGeo->GetNumberOfSuperModules();
    
    // Bad map
    if(fRemoveBadChannels)
    {
      AliOADBContainer *contBC=new AliOADBContainer("");
      contBC->InitFromFile(Form("%s/EMCALBadChannels.root",fOADBFilePathEMCAL.Data()),"AliEMCALBadChannels"); 
      
      TObjArray *arrayBC=(TObjArray*)contBC->GetObject(fRunNumber);
      
      if(arrayBC)
      {
        SwitchOnDistToBadChannelRecalculation();
        AliInfo("Remove EMCAL bad cells");
        
        for (Int_t i=0; i<nSM; ++i) 
        {
          TH2I *hbm = GetEMCALChannelStatusMap(i);
          
          if (hbm)
            delete hbm;
          
          hbm=(TH2I*)arrayBC->FindObject(Form("EMCALBadChannelMap_Mod%d",i));
          
          if (!hbm) 
          {
            AliError(Form("Can not get EMCALBadChannelMap_Mod%d",i));
            continue;
          }
          
          hbm->SetDirectory(0);
          SetEMCALChannelStatusMap(i,hbm);
          
        } // loop
      } else AliInfo("Do NOT remove EMCAL bad channels\n"); // run array
      
      delete contBC;
    }  // Remove bad
    
    // Energy Recalibration
    if(fRecalibration)
    {
      AliOADBContainer *contRF=new AliOADBContainer("");
      
      contRF->InitFromFile(Form("%s/EMCALRecalib.root",fOADBFilePathEMCAL.Data()),"AliEMCALRecalib");
      
      TObjArray *recal=(TObjArray*)contRF->GetObject(fRunNumber); 
      
      if(recal)
      {
        TObjArray *recalpass=(TObjArray*)recal->FindObject(pass);
        
        if(recalpass)
        {
          TObjArray *recalib=(TObjArray*)recalpass->FindObject("Recalib");
          
          if(recalib)
          {
            AliInfo("Recalibrate EMCAL");
            for (Int_t i=0; i < nSM; ++i)
            {
              TH2F *h = GetEMCALChannelRecalibrationFactors(i);
              
              if (h)
                delete h;
              
              h = (TH2F*)recalib->FindObject(Form("EMCALRecalFactors_SM%d",i));
              
              if (!h) 
              {
                AliError(Form("Could not load EMCALRecalFactors_SM%d",i));
                continue;
              }
              
              h->SetDirectory(0);
              
              SetEMCALChannelRecalibrationFactors(i,h);
            } // SM loop
          } else AliInfo("Do NOT recalibrate EMCAL, no params object array"); // array ok
        } else AliInfo("Do NOT recalibrate EMCAL, no params for pass"); // array pass ok
      } else AliInfo("Do NOT recalibrate EMCAL, no params for run");  // run number array ok
      
      delete contRF;
      // once set, apply run dependent corrections if requested
      //fEMCALRecoUtils->SetRunDependentCorrections(fRunNumber);
            
    } // Recalibration on
    
    // Energy Recalibration, apply on top of previous calibration factors
    if(fRunDependentCorrection)
    {
      AliOADBContainer *contRFTD=new AliOADBContainer("");
      
      contRFTD->InitFromFile(Form("%s/EMCALTemperatureCorrCalib.root",fOADBFilePathEMCAL.Data()),"AliEMCALRunDepTempCalibCorrections");
      
      TH1S *htd=(TH1S*)contRFTD->GetObject(fRunNumber); 
      
      //If it did not exist for this run, get closes one
      if (!htd)
      {
        AliWarning(Form("No TemperatureCorrCalib Objects for run: %d",fRunNumber));
        // let's get the closest fRunNumber instead then..
        Int_t lower = 0;
        Int_t ic = 0;
        Int_t maxEntry = contRFTD->GetNumberOfEntries();
        
        while ( (ic < maxEntry) && (contRFTD->UpperLimit(ic) < fRunNumber) ) 
        {
          lower = ic;
          ic++;
        }
        
        Int_t closest = lower;
        if ( (ic<maxEntry) &&
            (contRFTD->LowerLimit(ic)-fRunNumber) < (fRunNumber - contRFTD->UpperLimit(lower)) ) 
        {
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
              Float_t factor = GetEMCALChannelRecalibrationFactor(ism,icol,irow);
              
              Int_t absID = fEMCALGeo->GetAbsCellIdFromCellIndexes(ism, irow, icol); // original calibration factor
              factor *= htd->GetBinContent(absID) / 10000. ; // correction dependent on T
        
              //printf("\t ism %d, icol %d, irow %d,absID %d, corrA %2.3f, corrB %2.3f, corrAB %2.3f\n",ism, icol, irow, absID, 
              //      GetEMCALChannelRecalibrationFactor(ism,icol,irow) , htd->GetBinContent(absID) / 10000., factor);
        
              SetEMCALChannelRecalibrationFactor(ism,icol,irow,factor);
            } // columns
          } // rows 
        } // SM loop
      } else AliInfo("Do NOT recalibrate EMCAL with T variations, no params TH1");
      
      delete contRFTD;
    } // Run by Run T calibration    
    
    // Time Recalibration
    if(fEMCALRecoUtils->IsTimeRecalibrationOn())
    {
      AliOADBContainer *contTRF=new AliOADBContainer("");
      
      contTRF->InitFromFile(Form("%s/EMCALTimeCalib.root",fOADBFilePathEMCAL.Data()),"AliEMCALTimeCalib");
      
      TObjArray *trecal=(TObjArray*)contTRF->GetObject(fRunNumber); 
      
      if(trecal)
      {
        TString passM = pass;
        if(pass=="spc_calo") passM = "pass3";
        TObjArray *trecalpass=(TObjArray*)trecal->FindObject(passM);
        
        if(trecalpass)
        {
          AliInfo("Time Recalibrate EMCAL");
          for (Int_t ibc = 0; ibc < 4; ++ibc) 
          {
            TH1F *h = GetEMCALChannelTimeRecalibrationFactors(ibc);
            
            if (h)
              delete h;
            
            h = (TH1F*)trecalpass->FindObject(Form("hAllTimeAvBC%d",ibc));
            
            if (!h) 
            {
              AliError(Form("Could not load hAllTimeAvBC%d",ibc));
              continue;
            }
            
            h->SetDirectory(0);
            
            SetEMCALChannelTimeRecalibrationFactors(ibc,h);
          } // bunch crossing loop
        } else AliInfo("Do NOT recalibrate time EMCAL, no params for pass"); // array pass ok
      } else AliInfo("Do NOT recalibrate time EMCAL, no params for run");  // run number array ok
      
      delete contTRF;
    } // Time Recalibration on    
    
    // Time L1 phase racalibration    
    if(fEMCALRecoUtils->IsL1PhaseInTimeRecalibrationOn()) 
    {
      AliOADBContainer *contTRF=new AliOADBContainer("");
      contTRF->InitFromFile(Form("%s/EMCALTimeL1PhaseCalib.root",fOADBFilePathEMCAL.Data()),"AliEMCALTimeL1PhaseCalib");
      
      TObjArray *trecal=(TObjArray*)contTRF->GetObject(fRunNumber); 
      if(!trecal) 
      {
        AliError(Form("L1 phase time recal: No params for run %d. Default used.",fRunNumber));  // run number array ok
        trecal=(TObjArray*)contTRF->GetObject(0); // Try default object
      }
      
      if(trecal)
      {
        // Only 1 L1 phase correction possible, except special cases
        TString passM = "pass1";
        
        if ( pass == "muon_calo_pass1" && fRunNumber > 209121 && fRunNumber < 244284 ) 
          passM = "pass0";//period LHC15a-m

        TObjArray *trecalpass=(TObjArray*)trecal->FindObject(passM);
        if(!trecalpass) 
        {
          AliInfo(Form("L1 phase time recal: No params for run %d and pass %s, try default", fRunNumber, passM.Data())); 
          
          trecal->Delete();
          
          trecal=(TObjArray*)contTRF->GetObject(0);
          
          if(trecal)          
            trecalpass=(TObjArray*)trecal->FindObject("pass1");
          
          AliInfo("Time L1 phase Recalibrate EMCAL");
        }
        
        if(trecalpass)
        {
          TH1C *h =GetEMCALL1PhaseInTimeRecalibrationForAllSM();
          
          if (h) delete h;
          
          h = (TH1C*)trecalpass->FindObject(Form("h%d",fRunNumber));
          
          if (!h) AliError(Form("Could not load h%d",fRunNumber));
          
          h->SetDirectory(0);
          
          SetEMCALL1PhaseInTimeRecalibrationForAllSM(h);
        }
        else 
        {       
          AliError("Do not calibrate L1 phase time");
          fEMCALRecoUtils->SwitchOffL1PhaseInTimeRecalibration();
        }
      }
      else 
      {       
        AliError("Do not calibrate L1 phase time");
        fEMCALRecoUtils->SwitchOffL1PhaseInTimeRecalibration();
      }
      
      delete contTRF;
    }//End of Time L1 phase racalibration 
    
  }// EMCAL
  
  // PHOS
  if(fOADBForPHOS)
  {
    AliInfo(Form("Get AODB parameters from PHOS in %s for run %d, and <%s>",fOADBFilePathPHOS.Data(),fRunNumber,pass.Data()));
    
    // Bad map
    if(fRemoveBadChannels)
    {
      AliOADBContainer badmapContainer(Form("phosBadMap"));
      TString fileName="$ALICE_PHYSICS/OADB/PHOS/PHOSBadMaps.root";
      badmapContainer.InitFromFile(Form("%s/PHOSBadMaps.root",fOADBFilePathPHOS.Data()),"phosBadMap");
      
      //Use a fixed run number from year 2010, this year not available yet.
      TObjArray *maps = (TObjArray*)badmapContainer.GetObject(139000,"phosBadMap");
      if(!maps)
      {
        AliInfo(Form("Can not read PHOS bad map for run %d",fRunNumber)) ;
      }
      else
      {
        AliInfo(Form("Setting PHOS bad map with name %s",maps->GetName())) ;
       
        for(Int_t mod = 1; mod < 5; mod++)
        {
          TH2I *hbmPH = GetPHOSChannelStatusMap(mod);
          
          if(hbmPH) 
            delete hbmPH ;  
          
          hbmPH = (TH2I*)maps->At(mod);
          
          if(hbmPH) hbmPH->SetDirectory(0);
          
          SetPHOSChannelStatusMap(mod-1,hbmPH);
          
        } // modules loop  
      } // maps exist
    } // Remove bad channels
  } // PHOS
  
  // Parameters already set once, so do not it again
  fOADBSet = kTRUE;
}  

//_____________________________________________________________
/// Set the calorimeters transformation, alignmnet matrices 
/// and init geometry at least once.
//_____________________________________________________________
void AliCalorimeterUtils::AccessGeometry(AliVEvent* inputEvent) 
{  
  // First init the geometry, a priory not done before
  if(fRunNumber <=0 ) fRunNumber = inputEvent->GetRunNumber() ;
  
  InitPHOSGeometry ();
  InitEMCALGeometry();
  
  //Get the EMCAL transformation geometry matrices from ESD 
  if(!fEMCALGeoMatrixSet && fEMCALGeo)
  {
    if(fLoadEMCALMatrices)
    {
      AliInfo("Load user defined EMCAL geometry matrices");
      
      // OADB if available
      AliOADBContainer emcGeoMat("AliEMCALgeo");
      emcGeoMat.InitFromFile(Form("%s/EMCALlocal2master.root",fOADBFilePathEMCAL.Data()),"AliEMCALgeo");
      TObjArray *matEMCAL=(TObjArray*)emcGeoMat.GetObject(fRunNumber,"EmcalMatrices");
      
      for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
      {
        if (!fEMCALMatrix[mod]) // Get it from OADB
        {
          AliDebug(1,Form("EMCAL matrices SM %d, %p", mod,((TGeoHMatrix*) matEMCAL->At(mod))));
          //((TGeoHMatrix*) matEMCAL->At(mod))->Print();
          
          fEMCALMatrix[mod] = (TGeoHMatrix*) matEMCAL->At(mod) ;
        }
        
        if(fEMCALMatrix[mod])
        {
          if(fDebug > 1) 
            fEMCALMatrix[mod]->Print();
          
          fEMCALGeo->SetMisalMatrix(fEMCALMatrix[mod],mod) ;  
        }
        else if(gGeoManager)
        {
          AliWarning(Form("Set matrix for SM %d from gGeoManager",mod));
          fEMCALGeo->SetMisalMatrix(fEMCALGeo->GetMatrixForSuperModuleFromGeoManager(mod),mod) ;
        }
        else
        {
          AliError(Form("Alignment atrix for SM %d is not available",mod));
        }
      } // SM loop
      
      fEMCALGeoMatrixSet = kTRUE;//At least one, so good

    } // Load matrices
    else if (!gGeoManager) 
    { 
      AliDebug(1,"Load EMCAL misalignment matrices");
      if(!strcmp(inputEvent->GetName(),"AliESDEvent"))  
      {
        for(Int_t mod = 0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
        { 
          //printf("Load ESD matrix %d, %p\n",mod,((AliESDEvent*)inputEvent)->GetEMCALMatrix(mod));
          if(((AliESDEvent*)inputEvent)->GetEMCALMatrix(mod)) 
          {
            fEMCALGeo->SetMisalMatrix(((AliESDEvent*)inputEvent)->GetEMCALMatrix(mod),mod) ;
          }
        }// loop over super modules	
        
        fEMCALGeoMatrixSet = kTRUE;//At least one, so good
        
      }//ESD as input
      else 
      {
        AliDebug(1,"Setting of EMCAL transformation matrixes for AODs not implemented yet. \n Import geometry.root file");
      }//AOD as input
    }//Get matrix from data
    else if(gGeoManager)
    {
      fEMCALGeoMatrixSet = kTRUE;
    }
  }//EMCAL geo && no geoManager
  
	//Get the PHOS transformation geometry matrices from ESD 
  if(!fPHOSGeoMatrixSet && fPHOSGeo)
  {
    if(fLoadPHOSMatrices)
    {
      AliInfo("Load user defined PHOS geometry matrices");
      
      // OADB if available
      AliOADBContainer geomContainer("phosGeo");
      geomContainer.InitFromFile(Form("%s/PHOSGeometry.root",fOADBFilePathPHOS.Data()),"PHOSRotationMatrixes");
      TObjArray *matPHOS = (TObjArray*)geomContainer.GetObject(139000,"PHOSRotationMatrixes");    
      
      for(Int_t mod = 0 ; mod < 5 ; mod++)
      {
        if (!fPHOSMatrix[mod]) // Get it from OADB
        {
          AliDebug(1,Form("PHOS matrices module %d, %p",mod,((TGeoHMatrix*) matPHOS->At(mod))));
          //((TGeoHMatrix*) matPHOS->At(mod))->Print();
          
          fPHOSMatrix[mod] = (TGeoHMatrix*) matPHOS->At(mod) ;
        }
        
        // Set it, if it exists
        if(fPHOSMatrix[mod])
        {
          if(fDebug > 1 ) 
            fPHOSMatrix[mod]->Print();
          
          fPHOSGeo->SetMisalMatrix(fPHOSMatrix[mod],mod) ;  
        }      
      }// SM loop
      
      fPHOSGeoMatrixSet = kTRUE;//At least one, so good
      
    }//Load matrices
    else if (!gGeoManager) 
    { 
      AliDebug(1,"Load PHOS misalignment matrices.");
      if(!strcmp(inputEvent->GetName(),"AliESDEvent"))  
      {
        for(Int_t mod = 0; mod < 5; mod++)
        { 
	  if( ((AliESDEvent*)inputEvent)->GetPHOSMatrix(mod)) 
          {
	    //printf("PHOS: mod %d, matrix %p\n",mod, ((AliESDEvent*)inputEvent)->GetPHOSMatrix(mod));
	    fPHOSGeo->SetMisalMatrix( ((AliESDEvent*)inputEvent)->GetPHOSMatrix(mod),mod) ;
          }
        } // loop over modules

     fPHOSGeoMatrixSet  = kTRUE; //At least one so good
    } // ESD as input
    else 
    {
      AliDebug(1,"Setting of EMCAL transformation matrixes for AODs not implemented yet. \n Import geometry.root file");
    } // AOD as input
   } // get matrix from data
  else if(gGeoManager)
  {
    fPHOSGeoMatrixSet = kTRUE;
  }
 }//PHOS geo	and  geoManager was not set
}

//______________________________________________________________________________________
/// Decide if two cells are neighbours
/// A neighbour is defined as being two cells which share a side or corner.
//______________________________________________________________________________________
Bool_t AliCalorimeterUtils::AreNeighbours(Int_t calo, Int_t absId1, Int_t absId2 ) const
{
  Bool_t areNeighbours = kFALSE ;
  
  Int_t iRCU1 = -1, irow1 = -1, icol1 = -1;
  Int_t iRCU2 = -1, irow2 = -1, icol2 = -1;
  
  Int_t rowdiff =  0, coldiff =  0;
  
  Int_t nSupMod1 = GetModuleNumberCellIndexes(absId1, calo, icol1, irow1, iRCU1); 
  Int_t nSupMod2 = GetModuleNumberCellIndexes(absId2, calo, icol2, irow2, iRCU2); 
  
  if(calo==AliFiducialCut::kEMCAL && nSupMod1!=nSupMod2)
  {
    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2-1
    // C Side impair SM, nSupMod%2=1; A side pair SM nSupMod%2=0
    if(nSupMod1%2) icol1+=AliEMCALGeoParams::fgkEMCALCols;
    else           icol2+=AliEMCALGeoParams::fgkEMCALCols;    
  }
  
  rowdiff = TMath::Abs( irow1 - irow2 ) ;  
  coldiff = TMath::Abs( icol1 - icol2 ) ;  
  
  //if (( coldiff <= 1 )  && ( rowdiff <= 1 ) && (coldiff + rowdiff > 0))
  if ((coldiff + rowdiff == 1 ))
    areNeighbours = kTRUE ;
  
  return areNeighbours;
}

//_________________________________________________________________________________________
/// Checks if all of the cells in the cluster belongs to the same SuperModule.
/// EMCal (not DCal, except 1/3 SM and not PHOS) can share clusters at eta=0.
//_________________________________________________________________________________________
Bool_t AliCalorimeterUtils::IsClusterSharedByTwoSuperModules(const AliEMCALGeometry * geom,
                                                             AliVCluster* cluster)
{  
  Int_t    iSupMod = -1;
  Int_t    iSM0    = -1;
  Int_t    iTower  = -1;
  Int_t    iIphi   = -1;
  Int_t    iIeta   = -1;
  Int_t    iphi    = -1;
  Int_t    ieta    = -1;
  
  for(Int_t iDigit = 0; iDigit < cluster->GetNCells(); iDigit++)
  {
    // Get from the absid the supermodule, tower and eta/phi numbers
    geom->GetCellIndex(cluster->GetCellAbsId(iDigit),iSupMod,iTower,iIphi,iIeta);
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi,iIeta, iphi,ieta);
    
    // Check if there are cells of different SM
    if     (iDigit  == 0   ) iSM0 = iSupMod;
    else if(iSupMod != iSM0) 
    {
      if(iSupMod > 11 && iSupMod < 18) 
        AliWarning(Form("Cluster shared in 2 DCal: SM%d, SM%d??",iSupMod,iSM0));
      
      return kTRUE;
    }
  }
  
  return kFALSE;
}

//______________________________________________________________________________
/// Given the list of AbsId of the cluster, get the maximum cell and 
/// check if there are fNCellsFromBorder from the calorimeter border.
//______________________________________________________________________________
Bool_t AliCalorimeterUtils::CheckCellFiducialRegion(AliVCluster* cluster, 
                                                    AliVCaloCells* cells) const 
{
  //If the distance to the border is 0 or negative just exit accept all clusters
  
  if ( cells->GetType()==AliVCaloCells::kEMCALCell && fEMCALRecoUtils->GetNumberOfCellsFromEMCALBorder() <= 0 ) return kTRUE;
	
  if ( cells->GetType()==AliVCaloCells::kPHOSCell  && fNCellsFromPHOSBorder  <= 0 ) return kTRUE;
  
  Int_t absIdMax = -1;
  Float_t ampMax = -1;
	
  for(Int_t i = 0; i < cluster->GetNCells() ; i++)
  {
    Int_t absId = cluster->GetCellAbsId(i) ;
    Float_t amp	= cells->GetCellAmplitude(absId);

    if(amp > ampMax)
    {
      ampMax   = amp;
      absIdMax = absId;
    }
  }
	
  AliDebug(1,Form("Cluster Max AbsId %d, Cell Energy %2.2f, Cluster Energy %2.2f",
                  absIdMax, ampMax, cluster->E()));
	
  if ( absIdMax == -1 ) return kFALSE;
	
  // Check if the cell is close to the borders:
  Bool_t okrow = kFALSE;
  Bool_t okcol = kFALSE;
  
  if ( cells->GetType() == AliVCaloCells::kEMCALCell )
  {
    // It should be the same as AliEMCALRecoUtils::CheckCellFiducialRegion()
    // Why not calling it?
    
    Int_t iTower = -1, iIphi = -1, iIeta = -1, iphi = -1, ieta = -1, iSM = -1; 
    
    fEMCALGeo->GetCellIndex(absIdMax,iSM,iTower,iIphi,iIeta); 
		
    fEMCALGeo->GetCellPhiEtaIndexInSModule(iSM,iTower,iIphi, iIeta,iphi,ieta);
		
    if(iSM < 0 || iphi < 0 || ieta < 0 )
    {
      AliFatal(Form("Negative value for super module: %d, or cell ieta: %d, or cell iphi: %d, check EMCAL geometry name",iSM,ieta,iphi));
    }
    
    // Check rows/phi
    Int_t nborder = fEMCALRecoUtils->GetNumberOfCellsFromEMCALBorder();
    
    if ( iSM < 10 || (iSM > 11 && iSM < 18) ) // Full EMCal (SM0-9) and DCal 2/3 (SM12-17)
    {
      if(iphi >= nborder && iphi < 24-nborder) okrow = kTRUE; 
    }
    else // 1/3 SMs (SM10-11, SM18-19)
    {
      if(iphi >= nborder && iphi <  8-nborder) okrow = kTRUE; 
    }
		
    // Check columns/eta
    
    // Remove all borders if IsEMCALNoBorderAtEta0 or DCal SMs(12-17)
    if(!fEMCALRecoUtils->IsEMCALNoBorderAtEta0() || (iSM > 11 && iSM < 18))
    {
      if(ieta  > nborder && ieta < 48-nborder) okcol =kTRUE; 
    }
    else // Do not remove borders close at eta = 0 for Full EMCal SMs and 1/3 EMCal
    {
      if ( iSM%2 == 0 )
      {
        if(ieta >= nborder)     okcol = kTRUE;	
      }
      else 
      {
        if(ieta <  48-nborder)  okcol = kTRUE;	
      }
    } // eta 0 not checked
    
    AliDebug(1,Form("EMCAL Cluster in %d cells fiducial volume: ieta %d, iphi %d, SM %d ? ok row %d, ok column %d",
                    nborder, ieta, iphi, iSM,okrow,okcol));

  }//EMCAL
  else if ( cells->GetType() == AliVCaloCells::kPHOSCell )
  {
    Int_t relId[4];
    Int_t irow = -1, icol = -1;
    fPHOSGeo->AbsToRelNumbering(absIdMax,relId);
    
    if (relId[1] != 0 ) return kFALSE; // skip CPV only PHOS

    irow = relId[2];
    icol = relId[3];
    //imod = relId[0]-1;
	
    if(irow >= fNCellsFromPHOSBorder && irow < 64-fNCellsFromPHOSBorder) okrow =kTRUE;
		if(icol >= fNCellsFromPHOSBorder && icol < 56-fNCellsFromPHOSBorder) okcol =kTRUE; 

    AliDebug(1,Form("PHOS Cluster in %d cells fiducial volume: ieta %d, iphi %d, SM %d ? ok row %d, ok column %d",
                    fNCellsFromPHOSBorder, icol, irow, relId[0]-1,okrow,okcol));
    }//PHOS
	
    if (okcol && okrow) return kTRUE; 
    else                return kFALSE;
}	

//________________________________________________________________________________________________________
/// Check that in the cluster cells, there is no bad channel of those stored 
/// in fEMCALBadChannelMap or fPHOSBadChannelMap
//________________________________________________________________________________________________________
Bool_t AliCalorimeterUtils::ClusterContainsBadChannel(Int_t calorimeter, UShort_t* cellList, Int_t nCells)
{	
  if (!fRemoveBadChannels) return kFALSE;
  
  //printf("fEMCALBadChannelMap %p, fPHOSBadChannelMap %p \n",fEMCALBadChannelMap,fPHOSBadChannelMap);
  if(calorimeter == AliFiducialCut::kEMCAL && !fEMCALRecoUtils->GetEMCALChannelStatusMap(0)) return kFALSE;
  if(calorimeter == AliFiducialCut::kPHOS  && !fPHOSBadChannelMap)  return kFALSE;
  
  Int_t icol = -1;
  Int_t irow = -1;
  Int_t imod = -1;
  for(Int_t iCell = 0; iCell<nCells; iCell++)
  {
    // Get the column and row
    if ( calorimeter == AliFiducialCut::kEMCAL )
    {
      return fEMCALRecoUtils->ClusterContainsBadChannel((AliEMCALGeometry*)fEMCALGeo,cellList,nCells);
    }
    else if ( calorimeter == AliFiducialCut::kPHOS )
    {
      Int_t    relId[4];
      fPHOSGeo->AbsToRelNumbering(cellList[iCell],relId);
			
      if (relId[1] != 0 ) return kTRUE; // skip CPV only PHOS
      
      irow = relId[2];
      icol = relId[3];
      imod = relId[0]-1;
      
      if ( fPHOSBadChannelMap->GetEntries() <= imod ) continue;
      
      //printf("PHOS bad channels imod %d, icol %d, irow %d\n",imod, irow, icol);
      if ( GetPHOSChannelStatus(imod, irow, icol) ) return kTRUE;
    }
    else return kFALSE;
  } // cell cluster loop
  
  return kFALSE;
}

//_______________________________________________________________
/// Correct cluster energy non linearity.
//_______________________________________________________________
void AliCalorimeterUtils::CorrectClusterEnergy(AliVCluster *clus)
{  
  clus->SetE(fEMCALRecoUtils->CorrectClusterEnergyLinearity(clus));
}

//_______________________________________________________________
/// Select EMCal SM regions, depending on its location in a SM, 
/// behind frames, close to borders, etc.
/// Current regions are valid for EMCal, rethink for DCal and 1/3 SMs
///
/// \param clus: cluster, access to highest energy tower
/// \param cells: list of cells, needed to find highest energy tower
/// \param regEta: eta sub-region index
/// \param regPhi: phi sub-region index
///
/// \return integer with location
//______________________________________________________________________________
void AliCalorimeterUtils::GetEMCALSubregion(AliVCluster   * clus, AliVCaloCells * cells,
                                            Int_t & regEta, Int_t & regPhi) const
{
  regEta = regPhi = -1 ;

  if(!clus->IsEMCAL()) return ;
  
  Int_t icol = -1, irow = -1, iRCU = -1;
  Float_t clusterFraction = 0;
  
  Int_t absId = GetMaxEnergyCell(cells,clus,clusterFraction);
  
  Int_t sm    = GetModuleNumberCellIndexes(absId,AliFiducialCut::kEMCAL,icol,irow,iRCU);
  
  // Shift by 48 to for impair SM
  if( sm%2 == 1) icol+=AliEMCALGeoParams::fgkEMCALCols;
    
  // Avoid borders
  if(icol < 2 || icol > 93 || irow < 2 || irow > 21) return;
  
  //
  // Eta regions
  //
  
  // Region 0: center of SM ~0.18<|eta|<0.55
  if      ( icol >   9 && icol <  34 ) regEta = 0;
  else if ( icol >  62 && icol <  87 ) regEta = 0;
  
  // Region 3: frames ~0.1<|eta|<~0.22 ~0.51<|eta|<0.62
  
  else if ( icol <=  9 && icol >=  5 )  regEta =  3;
  else if ( icol <= 38 && icol >= 34 )  regEta =  3;
  else if ( icol <= 62 && icol >= 58 )  regEta =  3;
  else if ( icol <= 91 && icol >= 87 )  regEta =  3;

  // Region 1: |eta| < ~0.15 
  
  else if ( icol <  58  && icol >  38 ) regEta =  1 ;
  
  // Region 2: |eta| > ~0.6

  else                                  regEta =  2 ;
  
  //
  // Phi regions
  //
    
  if      ( irow >=  2 && irow <=  5 ) regPhi = 0; // External
  else if ( irow >= 18 && irow <= 21 ) regPhi = 0; // External
  else if ( irow >=  6 && irow <=  9 ) regPhi = 1; // Mid
  else if ( irow >= 14 && irow <= 17 ) regPhi = 1; // Mid
  else                                 regPhi = 2; //10-13 Central
  
}

//________________________________________________________________________________________
/// Fill array with 4 possibly correlated channels absId
///
///  \param absId: Reference absId cell
///  \param absIdCorr: List of cells correlated to absId, absId is included
///  \return true if 4 channels found
///
//________________________________________________________________________________________
Bool_t  AliCalorimeterUtils::GetFECCorrelatedCellAbsId(Int_t absId, Int_t absIdCorr[4]) const 
{
  // Get SM number
  Int_t sm = fEMCALGeo->GetSuperModuleNumber(absId);

  // Get first absId of SM
  Int_t absIdSMMin = fEMCALGeo->GetAbsCellIdFromCellIndexes(sm,0,1); // for absIds, first is in col 1, module/tower ordering map ...

  // Get reference n and correlated 3
  for(Int_t k = 0; k < 4; k++ )
  {
    for(Int_t p = 0; p < 72; p++ )
    {
     Int_t  n = absIdSMMin + 2*k + 16 *p;
     
      if ( absId == n   || absId == n+1 || 
           absId == n+8 || absId == n+9   ) 
      {
        absIdCorr[0] = n   ;
        absIdCorr[1] = n+1 ;
        absIdCorr[2] = n+8 ;
        absIdCorr[3] = n+9 ;
        
        //printf("n=%d, n+1=%d, n+8=%d, n+9=%d\n",
        //       absIdCorr[0],absIdCorr[1],absIdCorr[2],absIdCorr[3]);
        
        return kTRUE;
      }
    }
  }
  
  // Not found;
  absIdCorr[0] = -1 ;
  absIdCorr[1] = -1 ;
  absIdCorr[2] = -1 ;
  absIdCorr[3] = -1 ;
  
  return kFALSE;
}

//________________________________________________________________________________________
/// Check if 2 cells belong to the same TCard
///
///  \param absId1: Reference absId cell
///  \param absId2: Cross checked cell absId
///  \param rowDiff: Distance in rows
///  \param colDiff: Distance in columns
///  \return true if belong to same TCard
///
//________________________________________________________________________________________
Bool_t  AliCalorimeterUtils::IsAbsIDsFromTCard(Int_t absId1, Int_t absId2, 
                                               Int_t & rowDiff, Int_t & colDiff) const
{  
  rowDiff = -100;
  colDiff = -100;
  
  if(absId1 == absId2) return kFALSE;
  
  // Check if in same SM, if not for sure not same TCard
  Int_t sm1 = fEMCALGeo->GetSuperModuleNumber(absId1);
  Int_t sm2 = fEMCALGeo->GetSuperModuleNumber(absId2);
  if ( sm1 != sm2 ) return kFALSE ;
  
  // Get the column and row of each absId
  Int_t iTower = -1, iIphi = -1, iIeta = -1;

  Int_t col1, row1;
  fEMCALGeo->GetCellIndex(absId1,sm1,iTower,iIphi,iIeta);
  fEMCALGeo->GetCellPhiEtaIndexInSModule(sm1,iTower,iIphi, iIeta,row1,col1);
  
  Int_t col2, row2;
  fEMCALGeo->GetCellIndex(absId2,sm2,iTower,iIphi,iIeta);
  fEMCALGeo->GetCellPhiEtaIndexInSModule(sm2,iTower,iIphi, iIeta,row2,col2);
  
  Int_t row0 = Int_t(row1-row1%8);
  Int_t col0 = Int_t(col1-col1%2);
  
  Int_t rowDiff0 = row2-row0;
  Int_t colDiff0 = col2-col0;
  
  rowDiff = row1-row2;
  colDiff = col1-col2;
  
  // TCard is made by 2x8 towers
  if ( colDiff0 >=0 && colDiff0 < 2 && rowDiff0 >=0 && rowDiff0 < 8 ) 
  {
 
//    printf("\t absId (%d,%d), sm %d; col (%d,%d), colDiff %d; row (%d,%d),rowDiff %d\n",
//           absId1 , absId2, sm1, 
//           col1, col2, colDiff, 
//           row1, row2, rowDiff);
    return kTRUE ;
  }
  else
    return kFALSE;
}


//________________________________________________________________________________________
/// For a given CaloCluster, it gets the absId of the cell with maximum energy deposit.
//________________________________________________________________________________________
Int_t  AliCalorimeterUtils::GetMaxEnergyCell(AliVCaloCells * cells,
                                             AliVCluster   * clu, 
                                             Float_t & clusterFraction) const 
{  
  if( !clu || !cells )
  {
    AliInfo("Cluster or cells pointer is null!");
    return -1;
  }
  
  Double_t eMax        =-1.;
  Double_t eTot        = 0.;
  Double_t eCell       =-1.;
  Float_t  fraction    = 1.;
  Float_t  recalFactor = 1.;
  Int_t    cellAbsId   =-1 , absId =-1 ;
  Int_t    iSupMod     =-1 , ieta  =-1 , iphi = -1, iRCU = -1;
  
  Int_t             calo = AliFiducialCut::kEMCAL;
  if(clu->IsPHOS()) calo = AliFiducialCut::kPHOS ;
  
  for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++) 
  {
    cellAbsId = clu->GetCellAbsId(iDig);
    
    fraction  = clu->GetCellAmplitudeFraction(iDig);
    if(fraction < 1e-4) fraction = 1.; // in case unfolding is off
    
    iSupMod = GetModuleNumberCellIndexes(cellAbsId, calo, ieta, iphi, iRCU);
    
    if(IsRecalibrationOn())
    {
      if(calo == AliFiducialCut::kEMCAL) 
        recalFactor = GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
      else               
        recalFactor = GetPHOSChannelRecalibrationFactor (iSupMod,iphi,ieta);
    }
    
    eCell  = cells->GetCellAmplitude(cellAbsId)*fraction*recalFactor;
    
    if(eCell > eMax)  
    { 
      eMax  = eCell; 
      absId = cellAbsId;
    }
    
    eTot+=eCell;
    
  }// cell loop
  
  if(eTot > 0.1)
    clusterFraction = (eTot-eMax)/eTot; //Do not use cluster energy in case it was corrected for non linearity.
  else 
    clusterFraction =1.;
  
  return absId;
}

//___________________________________________________________________________________
/// Return the matched track to the cluster given its index.
/// It is usually just the first match, default value.
/// Since it is different for ESDs and AODs here it is a wrap method to do it.
//___________________________________________________________________________________
AliVTrack * AliCalorimeterUtils::GetMatchedTrack(AliVCluster* cluster,
                                                 AliVEvent* event, Int_t index) const
{  
  AliVTrack *track = 0x0;
  
  //
  // EMCAL case only when matching is recalculated
  //
  if(cluster->IsEMCAL() && IsRecalculationOfClusterTrackMatchingOn())
  {
    Int_t trackIndex = fEMCALRecoUtils->GetMatchedTrackIndex(cluster->GetID());
    //printf("track index %d, cluster ID %d \n ",trackIndex,cluster->GetID());
    
    if(trackIndex < 0 )
      AliInfo(Form("Wrong track index %d, from recalculation", trackIndex));
    else 
      track = dynamic_cast<AliVTrack*> (event->GetTrack(trackIndex));

    return track ;
  }   
  
  //
  // Normal case, get info from ESD or AOD
  //
  
  // No tracks matched
  if( cluster->GetNTracksMatched() < 1 ) return 0x0;
  
  // At least one match
  Int_t iTrack = 0; // only one match for AODs with index 0.
    
  // ESDs
  if(!strcmp("AliESDCaloCluster",Form("%s",cluster->ClassName())))
  {
    if( index >= 0 ) iTrack = index;
    else             iTrack = ((AliESDCaloCluster*)cluster)->GetTracksMatched()->At(0); //cluster->GetTrackMatchedIndex();
              
    track = dynamic_cast<AliVTrack*> ( event->GetTrack(iTrack) );    
  }
  else // AODs
  {        
    if( index > 0 ) iTrack = index;

    track = dynamic_cast<AliVTrack*>( cluster->GetTrackMatched(iTrack) );
  }
  
  return track ;
}

///
/// Correction factor for cell energy in cluster to temptatively match Data and MC.
/// Not used.
///
Float_t AliCalorimeterUtils::GetMCECellClusFracCorrection(Float_t eCell, Float_t eCluster) const
{
  if( eCluster <= 0 || eCluster < eCell )
  {
    AliWarning(Form("Bad values eCell=%f, eCluster %f",eCell,eCluster));
    return 1;
  }
  
  Float_t frac       = eCell / eCluster;
  
  Float_t correction = fMCECellClusFracCorrParam[0] +
                       TMath::Exp( frac*fMCECellClusFracCorrParam[2]+fMCECellClusFracCorrParam[1] ) +
                       fMCECellClusFracCorrParam[3]/TMath::Sqrt(frac);
  
//  printf("AliCalorimeterUtils::GetMCECellClusFracCorrection(eCell=%f, eCluster %f, frac %f) = %f\n",eCell, eCluster, frac, correction);
//  printf("\t %2.2f + TMath::Exp( %2.3f*%2.2f + %2.2f ) + %2.2f/TMath::Sqrt(%2.3f)) = %f\n",
//         fMCECellClusFracCorrParam[0],frac,fMCECellClusFracCorrParam[2],fMCECellClusFracCorrParam[1],fMCECellClusFracCorrParam[3], frac, correction);

  return correction;
}

//_____________________________________________________________________________________________________
/// Get the EMCAL/PHOS module number that corresponds to this particle.
//_____________________________________________________________________________________________________
Int_t AliCalorimeterUtils::GetModuleNumber(AliAODPWG4Particle * particle, AliVEvent * inputEvent) const
{	
  Int_t absId = -1;
  
  if(particle->GetDetectorTag()==AliFiducialCut::kEMCAL)
  {
    fEMCALGeo->GetAbsCellIdFromEtaPhi(particle->Eta(),particle->Phi(), absId);
    
    AliDebug(2,Form("EMCAL: cluster eta %f, phi %f, absid %d, SuperModule %d",
                    particle->Eta(), particle->Phi()*TMath::RadToDeg(),absId, fEMCALGeo->GetSuperModuleNumber(absId)));
    
    return fEMCALGeo->GetSuperModuleNumber(absId) ;
  } // EMCAL
  else if ( particle->GetDetectorTag() == AliFiducialCut::kPHOS )
  {
    // In case we use the MC reader, the input are TParticles,
    // in this case use the corresponing method in PHOS Geometry to get the particle.
    if(strcmp(inputEvent->ClassName(), "AliMCEvent") == 0 )
    {
      Int_t mod =-1;
      Double_t z = 0., x=0.;
      TParticle* primary = 0x0;
      AliStack * stack = ((AliMCEvent*)inputEvent)->Stack();
      
      if(stack)
      {
        primary = stack->Particle(particle->GetCaloLabel(0));
      }
      else
      {
        AliFatal("Stack not available, stop!");
      }
      
      if(primary)
      {
        fPHOSGeo->ImpactOnEmc(primary,mod,z,x) ;
      }
      else
      {
        AliFatal("Primary not available, stop!");
      }
      return mod;
    }
    // Input are ESDs or AODs, get the PHOS module number like this.
    else
    {
      //FIXME
      //AliVCluster *cluster = inputEvent->GetCaloCluster(particle->GetCaloLabel(0));
      //return GetModuleNumber(cluster);
      //MEFIX
      return -1;
    }
  } // PHOS
	
  return -1;
}

//_____________________________________________________________________
/// Get the EMCAL/PHOS module number that corresponds to this cluster.
//_____________________________________________________________________
Int_t AliCalorimeterUtils::GetModuleNumber(AliVCluster * cluster) const
{  
  if(!cluster)
  {
    AliDebug(1,"AliCalorimeterUtils::GetModuleNumber() - NUL Cluster, please check!!!");
    
    return -1;
  }
  
  if ( cluster->GetNCells() <= 0 ) return -1;
  
  Int_t absId = cluster->GetCellAbsId(0);
  
  if ( absId < 0 ) return -1;
  
  if( cluster->IsEMCAL() )
  {
    AliDebug(2,Form("EMCAL absid %d, SuperModule %d",absId, fEMCALGeo->GetSuperModuleNumber(absId)));
    
    return fEMCALGeo->GetSuperModuleNumber(absId) ;
  } // EMCAL
  else if ( cluster->IsPHOS() )
  {
    Int_t    relId[4];
    fPHOSGeo->AbsToRelNumbering(absId,relId);
    
    if (relId[1] != 0 ) return -1; // skip CPV only PHOS
    
    AliDebug(2,Form("PHOS absid %d Module %d",absId, relId[0]-1));
    
    return relId[0]-1;
  } // PHOS
	
  return -1;
}

//___________________________________________________________________________________________________
/// Get the EMCAL/PHOS module, columns, row and RCU/DDL number that corresponds to this absId.
//___________________________________________________________________________________________________
Int_t AliCalorimeterUtils::GetModuleNumberCellIndexes(Int_t absId, Int_t calo,
                                                      Int_t & icol, Int_t & irow, Int_t & iRCU) const
{
  Int_t imod = -1;
  
  if ( absId < 0 ) return -1 ;
  
  if ( calo == AliFiducialCut::kEMCAL )
  {
    Int_t iTower = -1, iIphi = -1, iIeta = -1;
    fEMCALGeo->GetCellIndex(absId,imod,iTower,iIphi,iIeta);
    fEMCALGeo->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,irow,icol);
    
    if(imod < 0 || irow < 0 || icol < 0 )
    {
      AliFatal(Form("Negative value for super module: %d, or cell icol: %d, or cell irow: %d, check EMCAL geometry name",imod,icol,irow));
    }
    
    // In case of DCal C side, shift columns to match offline/online numbering
    // when calculating the DDL for Run2
    // See AliEMCALRawUtils::Digits2Raw and Raw2Digits.
    Int_t ico2 = icol ;
    if ( imod == 13 || imod == 15 || imod == 17 ) ico2 += 16; 

    // RCU / DDL
    if(imod < 10 || (imod > 11 && imod < 18)) // (EMCAL Full || DCAL 2/3)
    {
      // RCU0 / DDL0
      if      ( 0 <= irow && irow <  8 ) iRCU = 0; // first cable row
      else if ( 8 <= irow && irow < 16 &&  
                0 <= ico2 && ico2 < 24 ) iRCU = 0; // first half;
      //second cable row
      
      // RCU1 / DDL1
      else if (  8 <= irow && irow < 16 && 
                24 <= ico2 && ico2 < 48 ) iRCU = 1; // second half;
      //second cable row
      else if ( 16 <= irow && irow < 24 ) iRCU = 1; // third cable row
      
      if ( imod%2 == 1 ) iRCU = 1 - iRCU; // swap for odd=C side, to allow us to cable both sides the same
    }
    else
    {
      // 1/3 SM have one single SRU, just assign RCU/DDL 0
      iRCU = 0 ;
    }
    
    if ( iRCU < 0 )
      AliFatal(Form("Wrong EMCAL RCU number = %d", iRCU));
    
    return imod ;
  } // EMCAL
  else // PHOS
  {
    Int_t relId[4];
    fPHOSGeo->AbsToRelNumbering(absId,relId);
    
    if (relId[1] != 0 ) return -1; // skip CPV only PHOS
    
    irow = relId[2];
    icol = relId[3];
    imod = relId[0]-1;
    iRCU= (Int_t)(relId[2]-1)/16 ;
    
    //Int_t iBranch= (Int_t)(relid[3]-1)/28 ; //0 to 1
    
    if ( iRCU >= 4 )
      AliFatal(Form("Wrong PHOS RCU number = %d", iRCU));
    
    return imod;
  } // PHOS
	
  return -1;
}

//___________________________________________________________________________________________________
/// Same as GetModuleCellIndexes, but add an additional shift in col/row to have continuous 
/// cell distribution from supermodule to supermodule.
//___________________________________________________________________________________________________
Int_t AliCalorimeterUtils::GetModuleNumberCellIndexesAbsCaloMap(Int_t absId    , Int_t calo,
                                                                Int_t & icol   , Int_t & irow   , Int_t & iRCU, 
                                                                Int_t & icolAbs, Int_t & irowAbs) const
{
  Int_t imod = GetModuleNumberCellIndexes(absId, calo, icol, irow,iRCU);

  icolAbs = icol;
  irowAbs = irow;
  
  //
  // PHOS
  //
  if(calo == AliFiducialCut::kPHOS) 
  {    
    irowAbs = irow + 64 * imod;

    return imod;
  }
  //
  // EMCal/DCal
  //
  else
  {
    //
    // Shift collumns in even SM
    Int_t shiftEta = 48;
    
    // Shift collumn even more due to smaller acceptance of DCal collumns
    if ( imod >  11 && imod < 18) shiftEta+=48/3;
    
    icolAbs = (imod % 2) ? icol + shiftEta : icol;	
    
    //
    // Shift rows per sector
    irowAbs = irow + 24 * Int_t(imod / 2); 
    
    // Shift row less due to smaller acceptance of SM 10 and 11 to count DCal rows
    if ( imod >  11 && imod < 20) irowAbs -= (2*24 / 3);
    
    return imod ;
  }
}


//___________________________________________________________________________________________
///  Find the number of local maxima in cluster.
//___________________________________________________________________________________________
Int_t AliCalorimeterUtils::GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells) 
{  
  const Int_t   nc = cluster->GetNCells();
  
  Int_t   absIdList[nc]; 
  Float_t maxEList[nc]; 
  
  Int_t nMax = GetNumberOfLocalMaxima(cluster, cells, absIdList, maxEList);
  
  return nMax;
}

//___________________________________________________________________________________________
/// Find the number of local maxima in cluster.
//___________________________________________________________________________________________
Int_t AliCalorimeterUtils::GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells,
                                                  Int_t *absIdList,     Float_t *maxEList) 
{
  Int_t iDigitN = 0 ;
  Int_t iDigit  = 0 ;
  Int_t absId1 = -1 ;
  Int_t absId2 = -1 ;
  const Int_t nCells = cluster->GetNCells();
  
  Float_t eCluster = RecalibrateClusterEnergy(cluster, cells);// recalculate cluster energy, avoid non lin correction.

  Float_t simuTotWeight = 0;
  if(fMCECellClusFracCorrOn)
  {
    simuTotWeight = RecalibrateClusterEnergyWeightCell(cluster, cells,eCluster);// same but apply a weight
    simuTotWeight/= eCluster;
  }
  
  Int_t calorimeter = AliFiducialCut::kEMCAL;
  if(!cluster->IsEMCAL()) calorimeter = AliFiducialCut::kPHOS;
  
  //printf("cluster : ncells %d \n",nCells);
  
  Float_t emax  = 0;
  Int_t   idmax =-1;
  for(iDigit = 0; iDigit < nCells ; iDigit++)
  {
    absIdList[iDigit] = cluster->GetCellsAbsId()[iDigit]  ; 
    Float_t en = cells->GetCellAmplitude(absIdList[iDigit]);
    RecalibrateCellAmplitude(en,calorimeter,absIdList[iDigit]);
    
    if(fMCECellClusFracCorrOn)
      en*=GetMCECellClusFracCorrection(en,eCluster)/simuTotWeight;
    
    if( en > emax )
    {
      emax  = en ;
      idmax = absIdList[iDigit] ;
    }
    //Int_t icol = -1, irow = -1, iRCU = -1;
    //Int_t sm = GetModuleNumberCellIndexes(absIdList[iDigit], calorimeter, icol, irow, iRCU) ;
    //printf("\t cell %d, id %d, sm %d, col %d, row %d, e %f\n", iDigit, absIdList[iDigit], sm, icol, irow, en );
  }
  
  for(iDigit = 0 ; iDigit < nCells; iDigit++) 
  {   
    if( absIdList[iDigit] >= 0 ) 
    {
      absId1 = cluster->GetCellsAbsId()[iDigit];
      
      Float_t en1 = cells->GetCellAmplitude(absId1);
      RecalibrateCellAmplitude(en1,calorimeter,absId1);  
      
      if(fMCECellClusFracCorrOn)
        en1*=GetMCECellClusFracCorrection(en1,eCluster)/simuTotWeight;
      
      //printf("%d : absIDi %d, E %f\n",iDigit, absId1,en1);
      
      for(iDigitN = 0; iDigitN < nCells; iDigitN++) 
      {	
        absId2 = cluster->GetCellsAbsId()[iDigitN] ;
        
        if(absId2==-1 || absId2==absId1) continue;
        
        //printf("\t %d : absIDj %d\n",iDigitN, absId2);
        
        Float_t en2 = cells->GetCellAmplitude(absId2);
        RecalibrateCellAmplitude(en2,calorimeter,absId2);
        
        if(fMCECellClusFracCorrOn)
          en2*=GetMCECellClusFracCorrection(en2,eCluster)/simuTotWeight;

        //printf("\t %d : absIDj %d, E %f\n",iDigitN, absId2,en2);
        
        if ( AreNeighbours(calorimeter, absId1, absId2) ) 
        {
          // printf("\t \t Neighbours \n");
          if ( en1 > en2 ) 
          {    
            absIdList[iDigitN] = -1 ;
            //printf("\t \t indexN %d not local max\n",iDigitN);
            // but may be digit too is not local max ?
            if(en1 < en2 + fLocMaxCutEDiff) {
              //printf("\t \t index %d not local max cause locMaxCutEDiff\n",iDigit);
              absIdList[iDigit] = -1 ;
            }
          }
          else 
          {
            absIdList[iDigit] = -1 ;
            //printf("\t \t index %d not local max\n",iDigitN);
            // but may be digitN too is not local max ?
            if(en1 > en2 - fLocMaxCutEDiff) 
            {
              absIdList[iDigitN] = -1 ; 
              //printf("\t \t indexN %d not local max cause locMaxCutEDiff\n",iDigit);
            }
          } 
        } // if Are neighbours
        //else printf("\t \t NOT Neighbours \n");
      } // while digitN
    } // slot not empty
  } // while digit
  
  iDigitN = 0 ;
  for(iDigit = 0; iDigit < nCells; iDigit++) 
  { 
    if( absIdList[iDigit] >= 0 )
    {
      absIdList[iDigitN] = absIdList[iDigit] ;
      
      Float_t en = cells->GetCellAmplitude(absIdList[iDigit]);
      RecalibrateCellAmplitude(en,calorimeter,absIdList[iDigit]);
      
      if(fMCECellClusFracCorrOn)
        en*=GetMCECellClusFracCorrection(en,eCluster)/simuTotWeight;
      
      if(en < fLocMaxCutE) continue; // Maxima only with seed energy at least
      
      maxEList[iDigitN] = en ;
      
      //printf("Local max %d, id %d, en %f\n", iDigit,absIdList[iDigitN],en);
      iDigitN++ ; 
    }
  }
  
  if ( iDigitN == 0 )
  {
    AliDebug(1,Form("No local maxima found, assign highest energy cell as maxima, id %d, en cell %2.2f, en cluster %2.2f",
                    idmax,emax,cluster->E()));
    iDigitN      = 1     ;
    maxEList[0]  = emax  ;
    absIdList[0] = idmax ; 
  }
  
  
  AliDebug(1,Form("In cluster E %2.2f (wth non lin. %2.2f), M02 %2.2f, M20 %2.2f, N maxima %d",
                  cluster->E(),eCluster, cluster->GetM02(),cluster->GetM20(), iDigitN));
  
//  if(fDebug > 1) for(Int_t imax = 0; imax < iDigitN; imax++)
//  {
//    printf(" \t i %d, absId %d, Ecell %f\n",imax,absIdList[imax],maxEList[imax]);
//  }
  
  return iDigitN ;
}

//____________________________________
/// Get passx from filename.
//____________________________________
TString AliCalorimeterUtils::GetPass()
{    
  if (!AliAnalysisManager::GetAnalysisManager()->GetTree()) 
  {
    AliError("AliCalorimeterUtils::GetPass() - Pointer to tree = 0, returning null");
    return TString("");
  }
  
  if (!AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()) 
  {
    AliError("AliCalorimeterUtils::GetPass() - Null pointer input file, returning null");
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
    AliInfo("Check that Energy calibration was enabled for this MC production in the tender, clusterizer or here!!");
                        
    return TString("LHC14a1a"); 
  }
  
  // No condition fullfilled, give a default value
  AliInfo("Pass number string not found");
  return TString("");            
}

//________________________________________
/// Initialize the parameters of the analysis.
//________________________________________
void AliCalorimeterUtils::InitParameters()
{
  fEMCALGeoName         = "";
  fPHOSGeoName          = "";
  
  fEMCALGeoMatrixSet    = kFALSE;
  fPHOSGeoMatrixSet     = kFALSE;
  
  fRemoveBadChannels    = kFALSE;
  
  fNCellsFromPHOSBorder = 0;
  
  fLocMaxCutE           = 0.1 ;
  fLocMaxCutEDiff       = 0.0 ;

  //  fMaskCellColumns = new Int_t[fNMaskCellColumns];
  //  fMaskCellColumns[0] = 6 ;  fMaskCellColumns[1] = 7 ;  fMaskCellColumns[2] = 8 ; 
  //  fMaskCellColumns[3] = 35;  fMaskCellColumns[4] = 36;  fMaskCellColumns[5] = 37; 
  //  fMaskCellColumns[6] = 12+AliEMCALGeoParams::fgkEMCALCols; fMaskCellColumns[7] = 13+AliEMCALGeoParams::fgkEMCALCols;
  //  fMaskCellColumns[8] = 40+AliEMCALGeoParams::fgkEMCALCols; fMaskCellColumns[9] = 41+AliEMCALGeoParams::fgkEMCALCols; 
  //  fMaskCellColumns[10]= 42+AliEMCALGeoParams::fgkEMCALCols; 
  
  fOADBSet      = kFALSE;
  fOADBForEMCAL = kTRUE ;            
  fOADBForPHOS  = kFALSE;
  
  fOADBFilePathEMCAL = "$ALICE_PHYSICS/OADB/EMCAL" ;
  fOADBFilePathPHOS  = "$ALICE_PHYSICS/OADB/PHOS"  ;
  
  fImportGeometryFromFile = kTRUE;
  fImportGeometryFilePath = "";
 
  fNSuperModulesUsed = 22;
  
  fMCECellClusFracCorrParam[0] = 0.78;
  fMCECellClusFracCorrParam[1] =-1.8;
  fMCECellClusFracCorrParam[2] =-6.3;
  fMCECellClusFracCorrParam[3] = 0.014;
}

//_____________________________________________________
/// Init PHOS bad channels map
//_____________________________________________________
void AliCalorimeterUtils::InitPHOSBadChannelStatusMap()
{
  AliDebug(1,"Init bad channel map");
  
  // In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fPHOSBadChannelMap = new TObjArray(5);	
  for (int i = 0; i < 5; i++)fPHOSBadChannelMap->Add(new TH2I(Form("PHOS_BadMap_mod%d",i),
                                                              Form("PHOS_BadMap_mod%d",i), 
                                                              64, 0, 64, 56, 0, 56));
  
  fPHOSBadChannelMap->SetOwner(kTRUE);
  fPHOSBadChannelMap->Compress();
  
  //In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);		
}

//______________________________________________________
/// Init PHOS recalibration factors
//______________________________________________________
void AliCalorimeterUtils::InitPHOSRecalibrationFactors()
{
  AliDebug(1,"Init recalibration map");

  // In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fPHOSRecalibrationFactors = new TObjArray(5);
  for (int i = 0; i < 5; i++)
  {
    fPHOSRecalibrationFactors->Add(new TH2F(Form("PHOSRecalFactors_Mod%d",i),
                                            Form("PHOSRecalFactors_Mod%d",i), 
                                            64, 0, 64, 56, 0, 56));
  }
  
  // Init the histograms with 1
  for (Int_t m = 0; m < 5; m++) 
  {
    for (Int_t i = 0; i < 56; i++) 
    {
      for (Int_t j = 0; j < 64; j++) 
      {
        SetPHOSChannelRecalibrationFactor(m,j,i,1.);
      }
    }
  }
  
  fPHOSRecalibrationFactors->SetOwner(kTRUE);
  fPHOSRecalibrationFactors->Compress();
	
  // In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);		
}

///
/// Initialize EMCAL geometry if it did not exist previously.
///
void AliCalorimeterUtils::InitEMCALGeometry()
{  
  if (fEMCALGeo) return;
  
  AliDebug(1,Form(" for run=%d",fRunNumber));
  
  if(fEMCALGeoName=="")
  {
    fEMCALGeo = AliEMCALGeometry::GetInstanceFromRunNumber(fRunNumber);
    AliInfo(Form("Get EMCAL geometry name to <%s> for run %d",fEMCALGeo->GetName(),fRunNumber));
  }
  else
  {
    fEMCALGeo = AliEMCALGeometry::GetInstance(fEMCALGeoName);
    AliInfo(Form("Set EMCAL geometry name to <%s>",fEMCALGeoName.Data()));
  }

  // Init geometry, I do not like much to do it like this ...
  if(fImportGeometryFromFile && !gGeoManager)
  {
    if(fImportGeometryFilePath=="") // If not specified, set location depending on run number
    {
      // "$ALICE_ROOT/EVE/alice-data/default_geo.root"
      if     (fRunNumber <  140000) fImportGeometryFilePath = "$ALICE_PHYSICS/OADB/EMCAL/geometry_2010.root";
      else if(fRunNumber <  171000) fImportGeometryFilePath = "$ALICE_PHYSICS/OADB/EMCAL/geometry_2011.root";
      else if(fRunNumber <  198000) fImportGeometryFilePath = "$ALICE_PHYSICS/OADB/EMCAL/geometry_2012.root"; // 2012-2013
      else                          fImportGeometryFilePath = "$ALICE_PHYSICS/OADB/EMCAL/geometry_2015.root"; // >= 2015
    }
    
    AliInfo(Form("Import %s",fImportGeometryFilePath.Data()));
    
    TGeoManager::Import(fImportGeometryFilePath) ; // default need file "geometry.root" in local dir!!!!
  }
  else if (!gGeoManager) AliInfo("Careful!, gGeoManager not loaded, load misalign matrices");
}

///
/// Initialize PHOS geometry if it did not exist previously.
///
void AliCalorimeterUtils::InitPHOSGeometry()
{ 
  if (fPHOSGeo) return;
  
  AliDebug(1,Form(" for run=%d",fRunNumber));
  
  if(fPHOSGeoName=="") fPHOSGeoName = "PHOSgeo";
  
  fPHOSGeo = new AliPHOSGeoUtils(fPHOSGeoName);
  
  //if (!gGeoManager) AliInfo("Careful!, gGeoManager not loaded, load misalign matrices");
}

//______________________________________________________________________________________________
// Check that a MC ESD is in the calorimeter acceptance
//______________________________________________________________________________________________
Bool_t AliCalorimeterUtils::IsMCParticleInCalorimeterAcceptance(Int_t calo, TParticle* particle)
{  
  if(!particle || (calo!=AliFiducialCut::kEMCAL && calo!=AliFiducialCut::kPHOS)) return kFALSE ;
    
  if( (!IsPHOSGeoMatrixSet () && calo == AliFiducialCut::kPHOS ) ||
      (!IsEMCALGeoMatrixSet() && calo == AliFiducialCut::kEMCAL)   )
  {
    AliFatal(Form("Careful Geo Matrix for calo <%d> is not set, use AliFidutialCut instead",calo));
    return kFALSE ;
  }

  if(calo == AliFiducialCut::kPHOS )
  {
    Int_t mod = 0 ;
    Double_t x = 0, z = 0 ;
    return GetPHOSGeometry()->ImpactOnEmc( particle, mod, z, x);
  }
  else if(calo == AliFiducialCut::kEMCAL)
  {
    Int_t absID = 0 ;
    Bool_t ok = GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(particle->Eta(),particle->Phi(),absID);
    if(ok)
    {
      Int_t icol = -1, irow = -1, iRCU = -1;
      Int_t nModule = GetModuleNumberCellIndexes(absID,calo, icol, irow, iRCU);
      Int_t status  = GetEMCALChannelStatus(nModule,icol,irow);
      if(status > 0) ok = kFALSE;
    }

    return ok ;
  }
  
  return kFALSE ;
}

//______________________________________________________________________________________________________
/// Check that a MC AOD is in the calorimeter acceptance.
//______________________________________________________________________________________________________
Bool_t AliCalorimeterUtils::IsMCParticleInCalorimeterAcceptance(Int_t calo, AliAODMCParticle* particle)
{  
  if(!particle || (calo!=AliFiducialCut::kEMCAL && calo!=AliFiducialCut::kPHOS)) return kFALSE ;
  
  if( (!IsPHOSGeoMatrixSet () && calo == AliFiducialCut::kPHOS ) ||
      (!IsEMCALGeoMatrixSet() && calo == AliFiducialCut::kEMCAL)   )
  {
    AliFatal(Form("Careful Geo Matrix for calo <%d> is not set, use AliFidutialCut instead",calo));
    return kFALSE ;
  }

  Float_t phi = particle->Phi();
  if(phi < 0) phi+=TMath::TwoPi();
  
  if(calo == AliFiducialCut::kPHOS )
  {
    Int_t mod = 0 ;
    Double_t x = 0, z = 0 ;
    Double_t vtx[]={ particle->Xv(), particle->Yv(), particle->Zv() } ;
    return GetPHOSGeometry()->ImpactOnEmc(vtx, particle->Theta(), phi, mod, z, x) ;
  }
  else if(calo == AliFiducialCut::kEMCAL)
  {
    Int_t absID = 0 ;
    Bool_t ok = GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(particle->Eta(),phi,absID);
    if(ok)
    {
      Int_t icol = -1, irow = -1, iRCU = -1;
      Int_t nModule = GetModuleNumberCellIndexes(absID,calo, icol, irow, iRCU);
      Int_t status  = GetEMCALChannelStatus(nModule,icol,irow);
      if(status > 0) ok = kFALSE;
    }

    return ok ;
  }
  
  return kFALSE ;
}

//_____________________________________________________________________________________________________
// Check that a TLorentzVector is in the calorimeter acceptance, give the cell number where it hit.
//_____________________________________________________________________________________________________
Bool_t AliCalorimeterUtils::IsMCParticleInCalorimeterAcceptance(Int_t calo, Float_t eta, Float_t theta,
                                                                Float_t phiOrg, Int_t & absID)
{  
  if(calo!=AliFiducialCut::kEMCAL && calo!=AliFiducialCut::kPHOS) return kFALSE ;
  
  if( (!IsPHOSGeoMatrixSet () && calo == AliFiducialCut::kPHOS ) ||
      (!IsEMCALGeoMatrixSet() && calo == AliFiducialCut::kEMCAL)   )
  {
    AliFatal(Form("Careful Geo Matrix for calo <%d> is not set, use AliFidutialCut instead",calo));
    return kFALSE ;
  }

  Float_t phi = phiOrg;
  if(phi < 0) phi+=TMath::TwoPi();

  if(calo == AliFiducialCut::kPHOS )
  {
    Int_t mod = 0 ;
    Double_t x = 0, z = 0 ;
    Double_t vtx[]={0,0,0} ;
    return GetPHOSGeometry()->ImpactOnEmc(vtx, theta, phi, mod, z, x) ;
  }
  else if(calo == AliFiducialCut::kEMCAL)
  {
    Bool_t ok = GetEMCALGeometry()->GetAbsCellIdFromEtaPhi(eta,phi,absID);
    if(ok)
    {
      Int_t icol = -1, irow = -1, iRCU = -1;
      Int_t nModule = GetModuleNumberCellIndexes(absID,calo, icol, irow, iRCU);
      Int_t status  = GetEMCALChannelStatus(nModule,icol,irow);
      if(status > 0) ok = kFALSE;
    }

    return ok ;
  }
  
  return kFALSE ;
}

//_______________________________________________________________________
/// Check if cell is in one of the regions where we have significant amount 
/// of material in front. Only EMCAL.
//_______________________________________________________________________
Bool_t AliCalorimeterUtils::MaskFrameCluster(Int_t iSM, Int_t ieta) const
{
  Int_t icol = ieta;
  if ( iSM%2 ) icol+=48; // Impair SM, shift index [0-47] to [48-96]
  
  if (fNMaskCellColumns && fMaskCellColumns) 
  {
    for (Int_t imask = 0; imask < fNMaskCellColumns; imask++) 
    {
      if(icol==fMaskCellColumns[imask]) return kTRUE;
    }
  }
  
  return kFALSE;
}

//_________________________________________________________
/// Print some relevant parameters set for the analysis
//_________________________________________________________
void AliCalorimeterUtils::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  printf("***** Print: %s %s ******\n", GetName(), GetTitle() ) ;
  printf("Remove Clusters with bad channels? %d\n",fRemoveBadChannels);
  printf("Remove Clusters with max cell at less than %d cells from EMCAL border and %d cells from PHOS border\n",
         fEMCALRecoUtils->GetNumberOfCellsFromEMCALBorder(), fNCellsFromPHOSBorder);
  if(fEMCALRecoUtils->IsEMCALNoBorderAtEta0()) printf("Do not remove EMCAL clusters at Eta = 0\n");
  printf("Recalibrate Clusters? %d, run by run  %d\n",fRecalibration,fRunDependentCorrection);
  printf("Recalculate Clusters Position? %d\n",fRecalculatePosition);
  printf("Recalculate Clusters Energy? %d\n",fCorrectELinearity);
  printf("Matching criteria: dR < %2.2f[cm], dZ < %2.2f[cm]\n",fCutR,fCutZ);
  
  printf("Recalibrate time? %d, With L1 phase run by run? %d\n",IsTimeRecalibrationOn(),IsL1PhaseInTimeRecalibrationOn());

  printf("Loc. Max. E > %2.2f\n",       fLocMaxCutE);
  printf("Loc. Max. E Diff > %2.2f\n",  fLocMaxCutEDiff);
  
  printf("    \n") ;
} 

//_____________________________________________________________________________________________
/// Recalculate cell energy if recalibration factor.
//_____________________________________________________________________________________________
void AliCalorimeterUtils::RecalibrateCellAmplitude(Float_t & amp, Int_t calo, Int_t id) const
{  
  Int_t icol     = -1; Int_t irow     = -1; Int_t iRCU     = -1;
  Int_t nModule  = GetModuleNumberCellIndexes(id,calo, icol, irow, iRCU);
  
  if (IsRecalibrationOn()) 
  {
    if(calo == AliFiducialCut::kPHOS)
    {
      amp *= GetPHOSChannelRecalibrationFactor(nModule,icol,irow);
    }
    else		                   
    {
      amp *= GetEMCALChannelRecalibrationFactor(nModule,icol,irow);
    }
  }
}

//____________________________________________________________________________________________________
/// Recalculate time if time recalibration available for EMCAL not ready for PHOS.
//____________________________________________________________________________________________________
void AliCalorimeterUtils::RecalibrateCellTime(Double_t & time, Int_t calo, Int_t id, Int_t bc) const
{  
  if ( calo == AliFiducialCut::kEMCAL && GetEMCALRecoUtils()->IsTimeRecalibrationOn() ) 
  {
    GetEMCALRecoUtils()->RecalibrateCellTime(id,bc,time);
  }
}


//____________________________________________________________________________________________________
/// Recalculate time L1 phase shift if time recalibration available for EMCAL.
//____________________________________________________________________________________________________
void AliCalorimeterUtils::RecalibrateCellTimeL1Phase(Double_t & time, Int_t calo, Int_t iSM, Int_t bunchCrossNumber) const
{  
  if ( calo == AliFiducialCut::kEMCAL && GetEMCALRecoUtils()->IsL1PhaseInTimeRecalibrationOn() ) 
  {
    GetEMCALRecoUtils()->RecalibrateCellTimeL1Phase(iSM, bunchCrossNumber, time);
  }
}


//__________________________________________________________________________
/// Recalibrate the cluster energy, considering the recalibration map and the energy of the cells that compose the cluster.
//__________________________________________________________________________
Float_t AliCalorimeterUtils::RecalibrateClusterEnergy(AliVCluster * cluster, 
                                                      AliVCaloCells * cells)
{  
  // Initialize some used variables
  Float_t frac  = 0., energy = 0.;  
  
  if(cells) 
  {
    //Get the cluster number of cells and list of absId, check what kind of cluster do we have.
    
    UShort_t * index    = cluster->GetCellsAbsId() ;
    Double_t * fraction = cluster->GetCellsAmplitudeFraction() ;
    
    Int_t ncells     = cluster->GetNCells();	
    
    Int_t calo = AliFiducialCut::kEMCAL;
    if(cluster->IsPHOS()) calo = AliFiducialCut::kPHOS ;
    
    // Loop on the cells, get the cell amplitude and recalibration factor, multiply and and to the new energy
    for(Int_t icell = 0; icell < ncells; icell++)
    {      
      Int_t absId = index[icell];
      
      frac =  fraction[icell];
      if(frac < 1e-3) frac = 1; //in case of EMCAL, this is set as 0, not used.
      
      Float_t amp = cells->GetCellAmplitude(absId);
      RecalibrateCellAmplitude(amp,calo, absId);
      
      AliDebug(2,Form("Recalibrate cell: calo <%d>, cell fraction %f, cell energy: before cal %f; after cal %f",
                      calo,frac,cells->GetCellAmplitude(absId),amp));
      
      energy += amp*frac;
    }
    
    AliDebug(1,Form("Energy before %f, after %f",cluster->E(),energy));
    
  } // cells available
  else
  {
    AliFatal("Cells pointer does not exist!");
  }
  
  return energy;
}

//_______________________________________________________________________________________________________
/// Recalibrate the cluster energy, considering the recalibration map and the energy of the cells that compose the cluster.
/// Also consider reweighting of cells energy.
//_______________________________________________________________________________________________________
Float_t AliCalorimeterUtils::RecalibrateClusterEnergyWeightCell(AliVCluster * cluster,
                                                                AliVCaloCells * cells, Float_t energyOrg)
{
  //Initialize some used variables
  Float_t frac  = 0., energy = 0.;
  
  if(cells)
  {
    // Get the cluster number of cells and list of absId, check what kind of cluster do we have.
    
    UShort_t * index    = cluster->GetCellsAbsId() ;
    Double_t * fraction = cluster->GetCellsAmplitudeFraction() ;
    
    Int_t ncells     = cluster->GetNCells();
    
    Int_t calo = AliFiducialCut::kEMCAL;
    if(cluster->IsPHOS()) calo = AliFiducialCut::kPHOS ;
    
    // Loop on the cells, get the cell amplitude and recalibration factor, multiply and and to the new energy
    for(Int_t icell = 0; icell < ncells; icell++)
    {
      Int_t absId = index[icell];
      
      frac =  fraction[icell];
      if(frac < 1e-3) frac = 1; //in case of EMCAL, this is set as 0, not used.
      
      Float_t amp = cells->GetCellAmplitude(absId);
      RecalibrateCellAmplitude(amp,calo, absId);
      
      amp*=GetMCECellClusFracCorrection(amp,energyOrg);
      
      AliDebug(2,Form("Recalibrate cell: calo <%d>, cell fraction %f, cell energy %f",
                      calo,frac,cells->GetCellAmplitude(absId)));
      
      energy += amp*frac;
    }
    
    AliDebug(1,Form("Energy before %f, after %f",cluster->E(),energy));
  } // cells available
  else
  {
    AliFatal("Cells pointer does not exist!");
  }
  
  return energy;
}

//__________________________________________________________________________________________
/// Recalculate EMCAL cluster position.
/// The cluster new position is already modified.
//__________________________________________________________________________________________
void AliCalorimeterUtils::RecalculateClusterPosition(AliVCaloCells* cells, AliVCluster* clu)
{
  fEMCALRecoUtils->RecalculateClusterPosition((AliEMCALGeometry*)fEMCALGeo, cells,clu);
}

//________________________________________________________________________________
/// Recalculate track matching and set the new residuals in the cluster.
///
/// \param event: pointer to input event
/// \param clusterArray: list of clusters
/// \param mc: access to MC event
///
//________________________________________________________________________________
void AliCalorimeterUtils::RecalculateClusterTrackMatching(AliVEvent * event, 
                                                          TObjArray* clusterArray,
                                                          AliMCEvent* mc) 
{   
  if (!fRecalculateMatching) return ; 
  
  fEMCALRecoUtils->FindMatches(event,clusterArray,fEMCALGeo,mc) ;
  
  Float_t dZ  = 2000;
  Float_t dR  = 2000;

  Int_t nClusters = event->GetNumberOfCaloClusters();
  if(clusterArray) nClusters = clusterArray->GetEntriesFast();
  
  AliVCluster * clus = 0;

  for (Int_t iclus =  0; iclus < nClusters ; iclus++)
  {
    if  ( clusterArray ) clus = (AliVCluster*) clusterArray->At(iclus) ;
    else                 clus = event->GetCaloCluster(iclus) ;
    
    if (!clus->IsEMCAL()) continue ;
    
    //
    // Put track residuals in cluster
    //
    fEMCALRecoUtils->GetMatchedResiduals(clus->GetID(),dZ,dR);
    
    if ( TMath::Abs(clus->GetTrackDx()) < 500 )
      AliDebug(2,Form("Residuals (Old, New): z (%2.4f,%2.4f), x (%2.4f,%2.4f)\n",
                      clus->GetTrackDz(),dZ,clus->GetTrackDx(),dR));
    
    clus->SetTrackDistance(dR,dZ);
    
    //
    // Remove old matches in cluster
    //
    if(clus->GetNTracksMatched() > 0)
    {
      if(!strcmp("AliESDCaloCluster",Form("%s",clus->ClassName())))
      {
        TArrayI arrayTrackMatched(0);
        ((AliESDCaloCluster*)clus)->AddTracksMatched(arrayTrackMatched);
      }
      else
      {
        for(Int_t iTrack = 0; iTrack < clus->GetNTracksMatched(); iTrack++)
        {
          ((AliAODCaloCluster*)clus)->RemoveTrackMatched((TObject*)((AliAODCaloCluster*)clus)->GetTrackMatched(iTrack));
        }
      }
    }
    
    //
    // Now put first track index in cluster. 
    //
    Int_t trackIndex = fEMCALRecoUtils->GetMatchedTrackIndex(iclus);
    if ( trackIndex >= 0 )
    {
      if(!strcmp("AliESDCaloCluster",Form("%s",clus->ClassName())))
      {
        TArrayI arrayTrackMatched(1);
        arrayTrackMatched[0] = trackIndex;
        ((AliESDCaloCluster*)clus)->AddTracksMatched(arrayTrackMatched);
      }
      else
      {
        ((AliAODCaloCluster*)clus)->AddTrackMatched((TObject*)event->GetTrack(trackIndex));
      }
    } 
    
  } // cluster loop
}

//___________________________________________________________________________
/// Split energy of cluster between the 2 local maxima, sum energy on 3x3, and if the 2 
/// maxima are too close and have common cells, split the energy between the 2.
/// \param absId1: index of highest energy cell in the cluster.
/// \param absId2: index of second highest energy cell in the cluster.
/// \param cluster: original cluster pointer.
/// \param cells: list of cells.
/// \param cluster1: output sub-cluster.
/// \param cluster2: output sub-cluster.
/// \param nMax: Number of local maxima of original cluster.
/// \param eventNumber: Event number needed for debugging and plotting.
/// 
/// Posibility to plot the clusters and sub-clusters.
//___________________________________________________________________________
void AliCalorimeterUtils::SplitEnergy(Int_t absId1, Int_t absId2,
                                      AliVCluster* cluster,
                                      AliVCaloCells* cells,
                                      //Float_t & e1, Float_t & e2,
                                      AliAODCaloCluster* cluster1,
                                      AliAODCaloCluster* cluster2,
                                      Int_t nMax, Int_t eventNumber)
{  
  TH2F* hClusterMap    = 0 ;
  TH2F* hClusterLocMax = 0 ;
  TH2F* hCluster1      = 0 ;
  TH2F* hCluster2      = 0 ;
  
  if(fPlotCluster)
  {
    hClusterMap    = new TH2F("hClusterMap","Cluster Map",48,0,48,24,0,24);
    hClusterLocMax = new TH2F("hClusterLocMax","Cluster 2 highest local maxima",48,0,48,24,0,24);
    hCluster1      = new TH2F("hCluster1","Cluster 1",48,0,48,24,0,24);
    hCluster2      = new TH2F("hCluster2","Cluster 2",48,0,48,24,0,24);
    hClusterMap    ->SetXTitle("column");
    hClusterMap    ->SetYTitle("row");
    hClusterLocMax ->SetXTitle("column");
    hClusterLocMax ->SetYTitle("row");
    hCluster1      ->SetXTitle("column");
    hCluster1      ->SetYTitle("row");
    hCluster2      ->SetXTitle("column");
    hCluster2      ->SetYTitle("row");
  }
  
  Int_t calorimeter = AliFiducialCut::kEMCAL;
  if(cluster->IsPHOS())
  {
    calorimeter = AliFiducialCut::kPHOS;
    AliWarning("Not supported for PHOS yet");
    return;
  }
  
  const Int_t ncells  = cluster->GetNCells();  
  Int_t absIdList[ncells]; 
  
  Float_t e1 = 0,  e2   = 0 ;
  Int_t icol = -1, irow = -1, iRCU = -1, sm = -1;  
  Float_t eCluster = 0;
  Float_t minCol = 100, minRow = 100, maxCol = -1, maxRow = -1; 
  for(Int_t iDigit  = 0; iDigit < ncells; iDigit++ ) 
  {
    absIdList[iDigit] = cluster->GetCellsAbsId()[iDigit];
    
    Float_t ec = cells->GetCellAmplitude(absIdList[iDigit]);
    RecalibrateCellAmplitude(ec,calorimeter, absIdList[iDigit]);
    eCluster+=ec;
    
    if(fPlotCluster) 
    {
      //printf("iDigit %d, absId %d, Ecell %f\n",iDigit,absIdList[iDigit], cells->GetCellAmplitude(absIdList[iDigit]));
      sm = GetModuleNumberCellIndexes(absIdList[iDigit], calorimeter, icol, irow, iRCU) ;
      if(sm > -1 && sm  < 12) // just to avoid compilation warning
      {
        if(icol > maxCol) maxCol = icol;
        if(icol < minCol) minCol = icol;
        if(irow > maxRow) maxRow = irow;
        if(irow < minRow) minRow = irow;
        hClusterMap->Fill(icol,irow,ec);
      }
    }
  }
  
  // Init counters and variables
  Int_t ncells1 = 1 ;
  UShort_t absIdList1[9] ;  
  Double_t fracList1 [9] ;  
  absIdList1[0] = absId1 ;
  fracList1 [0] = 1. ;
  
  Float_t ecell1 = cells->GetCellAmplitude(absId1);
  RecalibrateCellAmplitude(ecell1, calorimeter, absId1);
  e1 =  ecell1;  
  
  Int_t ncells2 = 1 ;
  UShort_t absIdList2[9] ;  
  Double_t fracList2 [9] ; 
  absIdList2[0] = absId2 ;
  fracList2 [0] = 1. ;
  
  Float_t ecell2 = cells->GetCellAmplitude(absId2);
  RecalibrateCellAmplitude(ecell2, calorimeter, absId2);
  e2 =  ecell2;  
  
  if(fPlotCluster)
  {
    Int_t icol1 = -1, irow1 = -1, icol2 = -1, irow2 = -1;
    sm = GetModuleNumberCellIndexes(absId1, calorimeter, icol1, irow1, iRCU) ;
    hClusterLocMax->Fill(icol1,irow1,ecell1);
    sm = GetModuleNumberCellIndexes(absId2, calorimeter, icol2, irow2, iRCU) ;
    hClusterLocMax->Fill(icol2,irow2,ecell2);
  }
  
  // Very rough way to share the cluster energy
  Float_t eRemain = (eCluster-ecell1-ecell2)/2;
  Float_t shareFraction1 = ecell1/eCluster+eRemain/eCluster;
  Float_t shareFraction2 = ecell2/eCluster+eRemain/eCluster;
  
  for(Int_t iDigit = 0; iDigit < ncells; iDigit++)
  {
    Int_t absId = absIdList[iDigit];
    
    if ( absId==absId1 || absId==absId2 || absId < 0 ) continue;
    
    Float_t ecell = cells->GetCellAmplitude(absId);
    RecalibrateCellAmplitude(ecell, calorimeter, absId);
    
    if(AreNeighbours(calorimeter, absId1,absId ))
    { 
      absIdList1[ncells1]= absId;
      
      if(AreNeighbours(calorimeter, absId2,absId ))
      { 
        fracList1[ncells1] = shareFraction1; 
        e1 += ecell*shareFraction1;
      }
      else 
      {
        fracList1[ncells1] = 1.; 
        e1 += ecell;
      }
      
      ncells1++;
      
    } // neigbour to cell1
    
    if(AreNeighbours(calorimeter, absId2,absId ))
    { 
      absIdList2[ncells2]= absId;
      
      if(AreNeighbours(calorimeter, absId1,absId ))
      { 
        fracList2[ncells2] = shareFraction2; 
        e2 += ecell*shareFraction2;
      }
      else
      { 
        fracList2[ncells2] = 1.; 
        e2 += ecell;
      }
      
      ncells2++;
      
    } // neigbour to cell2  
  }
  
  AliDebug(1,Form("N Local Max %d, Cluster energy  = %f, Ecell1 = %f, Ecell2 = %f, Enew1 = %f, Enew2 = %f, Remain %f, \n ncells %d, ncells1 %d, ncells2 %d, f1 %f, f2  %f, sum f12 = %f",
                  nMax, eCluster,ecell1,ecell2,e1,e2,eCluster-e1-e2,ncells,ncells1,ncells2,shareFraction1,shareFraction2,shareFraction1+shareFraction2));
           
  cluster1->SetE(e1);
  cluster2->SetE(e2);  
  
  cluster1->SetNCells(ncells1);
  cluster2->SetNCells(ncells2);  
  
  cluster1->SetCellsAbsId(absIdList1);
  cluster2->SetCellsAbsId(absIdList2);
  
  cluster1->SetCellsAmplitudeFraction(fracList1);
  cluster2->SetCellsAmplitudeFraction(fracList2);
  
  // Correct linearity
  CorrectClusterEnergy(cluster1) ;
  CorrectClusterEnergy(cluster2) ;
  
  if(calorimeter==AliFiducialCut::kEMCAL)
  {
    GetEMCALRecoUtils()->RecalculateClusterPosition(GetEMCALGeometry(), cells, cluster1);
    GetEMCALRecoUtils()->RecalculateClusterPosition(GetEMCALGeometry(), cells, cluster2);
  }
  
  if(fPlotCluster)
  {
    //printf("Cells of cluster1: ");
    for(Int_t iDigit  = 0; iDigit < ncells1; iDigit++ ) 
    {
      //printf(" %d ",absIdList1[iDigit]);
      
      sm = GetModuleNumberCellIndexes(absIdList1[iDigit], calorimeter, icol, irow, iRCU) ;
      
      Float_t ecell = cells->GetCellAmplitude(absIdList1[iDigit]);
      RecalibrateCellAmplitude(ecell, calorimeter, absIdList1[iDigit]);
      
      if( AreNeighbours(calorimeter, absId2,absIdList1[iDigit]) && absId1!=absIdList1[iDigit])
        hCluster1->Fill(icol,irow,ecell*shareFraction1);
      else 
        hCluster1->Fill(icol,irow,ecell);
    }
    
    //printf(" \n ");
    //printf("Cells of cluster2: ");
    
    for(Int_t iDigit  = 0; iDigit < ncells2; iDigit++ ) 
    {
      //printf(" %d ",absIdList2[iDigit]);
      
      sm = GetModuleNumberCellIndexes(absIdList2[iDigit], calorimeter, icol, irow, iRCU) ;
      
      Float_t ecell = cells->GetCellAmplitude(absIdList2[iDigit]);
      RecalibrateCellAmplitude(ecell, calorimeter, absIdList2[iDigit]);
      
      if( AreNeighbours(calorimeter, absId1,absIdList2[iDigit])  && absId2!=absIdList2[iDigit])
        hCluster2->Fill(icol,irow,ecell*shareFraction2);
      else
        hCluster2->Fill(icol,irow,ecell); 
    }
    //printf(" \n ");
    
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadLeftMargin(0.1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(000000);
    
    if(maxCol-minCol > maxRow-minRow)
    {
      maxRow+= maxCol-minCol;
    }
    else 
    {
      maxCol+= maxRow-minRow;
    }
    
    TCanvas  * c= new TCanvas("canvas", "canvas", 4000, 4000) ;
    c->Divide(2,2);  
    c->cd(1);
    gPad->SetGridy();
    gPad->SetGridx();
    gPad->SetLogz();
    hClusterMap    ->SetAxisRange(minCol, maxCol,"X");
    hClusterMap    ->SetAxisRange(minRow, maxRow,"Y");
    hClusterMap    ->Draw("colz TEXT");
    c->cd(2);
    gPad->SetGridy();
    gPad->SetGridx();
    gPad->SetLogz();
    hClusterLocMax ->SetAxisRange(minCol, maxCol,"X");
    hClusterLocMax ->SetAxisRange(minRow, maxRow,"Y");
    hClusterLocMax ->Draw("colz TEXT");
    c->cd(3);
    gPad->SetGridy();
    gPad->SetGridx();
    gPad->SetLogz();
    hCluster1      ->SetAxisRange(minCol, maxCol,"X");
    hCluster1      ->SetAxisRange(minRow, maxRow,"Y");
    hCluster1      ->Draw("colz TEXT");
    c->cd(4);
    gPad->SetGridy();
    gPad->SetGridx();
    gPad->SetLogz();
    hCluster2      ->SetAxisRange(minCol, maxCol,"X");
    hCluster2      ->SetAxisRange(minRow, maxRow,"Y");
    hCluster2      ->Draw("colz TEXT");
    
    if(eCluster > 6 )c->Print(Form("clusterFigures/Event%d_E%1.0f_nMax%d_NCell1_%d_NCell2_%d.eps",
                                   eventNumber,cluster->E(),nMax,ncells1,ncells2));
    
    delete c;
    delete hClusterMap;
    delete hClusterLocMax;
    delete hCluster1;
    delete hCluster2;
  }
}

