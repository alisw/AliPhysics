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
// Class utility for Calorimeter specific selection methods                ///
//
//
//
//-- Author: Gustavo Conesa (LPSC-Grenoble) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include <TGeoManager.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <TFile.h>

// --- ANALYSIS system ---
#include "AliCalorimeterUtils.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAODPWG4Particle.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliMixedEvent.h"
#include "AliAODCaloCluster.h"
#include "AliOADBContainer.h"
#include "AliAnalysisManager.h"

// --- Detector ---
#include "AliEMCALGeometry.h"
#include "AliPHOSGeoUtils.h"

ClassImp(AliCalorimeterUtils)
  
  
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
    fRecalibration(kFALSE),           fRunDependentCorrection(kFALSE), fPHOSRecalibrationFactors(),
    fEMCALRecoUtils(new AliEMCALRecoUtils),
    fRecalculatePosition(kFALSE),     fCorrectELinearity(kFALSE),
    fRecalculateMatching(kFALSE),
    fCutR(20),                        fCutZ(20),
    fCutEta(20),                      fCutPhi(20),
    fLocMaxCutE(0),                   fLocMaxCutEDiff(0),
    fPlotCluster(0),                  fOADBSet(kFALSE),
    fOADBForEMCAL(kFALSE),            fOADBForPHOS(kFALSE),
    fOADBFilePathEMCAL(""),           fOADBFilePathPHOS("")
{
  //Ctor
  
  //Initialize parameters
  InitParameters();
  for(Int_t i = 0; i < 12; i++) fEMCALMatrix[i] = 0 ;
  for(Int_t i = 0; i < 5 ; i++) fPHOSMatrix [i] = 0 ;
  
}

//_________________________________________
AliCalorimeterUtils::~AliCalorimeterUtils()
{
  //Dtor
  
  //if(fPHOSGeo)  delete fPHOSGeo  ;
  if(fEMCALGeo) delete fEMCALGeo ;
  
  if(fPHOSBadChannelMap) { 
    fPHOSBadChannelMap->Clear();
    delete  fPHOSBadChannelMap;
  }
	
  if(fPHOSRecalibrationFactors) { 
    fPHOSRecalibrationFactors->Clear();
    delete  fPHOSRecalibrationFactors;
  }
	
  if(fEMCALRecoUtils)   delete fEMCALRecoUtils ;
  if(fNMaskCellColumns) delete [] fMaskCellColumns;
  
}

//____________________________________________________
void AliCalorimeterUtils::AccessOADB(AliVEvent* event)
{
  // Set the AODB calibration, bad channels etc. parameters at least once
  // alignment matrices from OADB done in SetGeometryMatrices
  
  //Set it only once
  if(fOADBSet) return ; 
  
  Int_t   runnumber = event->GetRunNumber() ;
  TString pass      = GetPass();
  
  // EMCAL
  if(fOADBForEMCAL)
  {
    printf("AliCalorimeterUtils::SetOADBParameters() - Get AODB parameters from EMCAL in %s for run %d, and <%s> \n",fOADBFilePathEMCAL.Data(),runnumber,pass.Data());
    
    Int_t nSM = fEMCALGeo->GetNumberOfSuperModules();
    
    // Bad map
    if(fRemoveBadChannels)
    {
      AliOADBContainer *contBC=new AliOADBContainer("");
      contBC->InitFromFile(Form("%s/EMCALBadChannels.root",fOADBFilePathEMCAL.Data()),"AliEMCALBadChannels"); 
      
      TObjArray *arrayBC=(TObjArray*)contBC->GetObject(runnumber);
      
      if(arrayBC)
      {
        SwitchOnDistToBadChannelRecalculation();
        printf("AliCalorimeterUtils::SetOADBParameters() - Remove EMCAL bad cells \n");
        
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
      } else printf("AliCalorimeterUtils::SetOADBParameters() - Do NOT remove EMCAL bad channels\n"); // run array
    }  // Remove bad
    
    // Energy Recalibration
    if(fRecalibration)
    {
      AliOADBContainer *contRF=new AliOADBContainer("");
      
      contRF->InitFromFile(Form("%s/EMCALRecalib.root",fOADBFilePathEMCAL.Data()),"AliEMCALRecalib");
      
      TObjArray *recal=(TObjArray*)contRF->GetObject(runnumber); 
      
      if(recal)
      {
        TObjArray *recalpass=(TObjArray*)recal->FindObject(pass);
        
        if(recalpass)
        {
          TObjArray *recalib=(TObjArray*)recalpass->FindObject("Recalib");
          
          if(recalib)
          {
            printf("AliCalorimeterUtils::SetOADBParameters() - Recalibrate EMCAL \n");
            for (Int_t i=0; i<nSM; ++i) 
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
          }else printf("AliCalorimeterUtils::SetOADBParameters() - Do NOT recalibrate EMCAL, no params object array\n"); // array ok
        }else printf("AliCalorimeterUtils::SetOADBParameters() - Do NOT recalibrate EMCAL, no params for pass\n"); // array pass ok
      }else printf("AliCalorimeterUtils::SetOADBParameters() - Do NOT recalibrate EMCAL, no params for run\n");  // run number array ok
      
      // once set, apply run dependent corrections if requested
      //fEMCALRecoUtils->SetRunDependentCorrections(runnumber);
            
    } // Recalibration on
    
    // Energy Recalibration, apply on top of previous calibration factors
    if(fRunDependentCorrection)
    {
      AliOADBContainer *contRFTD=new AliOADBContainer("");
      
      contRFTD->InitFromFile(Form("%s/EMCALTemperatureCorrCalib.root",fOADBFilePathEMCAL.Data()),"AliEMCALRunDepTempCalibCorrections");
      
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
        printf("AliCalorimeterUtils::SetOADBParameters() - Recalibrate (Temperature) EMCAL \n");
        
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
      }else printf("AliCalorimeterUtils::SetOADBParameters() - Do NOT recalibrate EMCAL with T variations, no params TH1 \n"); 
    } // Run by Run T calibration    
    
    // Time Recalibration
    if(fEMCALRecoUtils->IsTimeRecalibrationOn())
    {
      AliOADBContainer *contTRF=new AliOADBContainer("");
      
      contTRF->InitFromFile(Form("%s/EMCALTimeCalib.root",fOADBFilePathEMCAL.Data()),"AliEMCALTimeCalib");
      
      TObjArray *trecal=(TObjArray*)contTRF->GetObject(runnumber); 
      
      if(trecal)
      {
        TString passTmp = pass;
        if(pass!="pass1" && pass!="pass2") passTmp = "pass2"; // TEMPORARY FIX FOR LHC11a analysis
        
        TObjArray *trecalpass=(TObjArray*)trecal->FindObject(passTmp);
        
        if(trecalpass)
        {
          printf("AliCalorimeterUtils::SetOADBParameters() - Time Recalibrate EMCAL \n");
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
        }else printf("AliCalorimeterUtils::SetOADBParameters() - Do NOT recalibrate time EMCAL, no params for pass\n"); // array pass ok
      }else printf("AliCalorimeterUtils::SetOADBParameters() - Do NOT recalibrate time EMCAL, no params for run\n");  // run number array ok
      
    } // Recalibration on    
    
  }// EMCAL
  
  // PHOS
  if(fOADBForPHOS)
  {
    printf("AliCalorimeterUtils::SetOADBParameters() - Get AODB parameters from PHOS in %s for run %d, and <%s> \n",fOADBFilePathPHOS.Data(),runnumber,pass.Data());
    
    // Bad map
    if(fRemoveBadChannels)
    {
      AliOADBContainer badmapContainer(Form("phosBadMap"));
      TString fileName="$ALICE_ROOT/OADB/PHOS/PHOSBadMaps.root";
      badmapContainer.InitFromFile(Form("%s/PHOSBadMaps.root",fOADBFilePathPHOS.Data()),"phosBadMap");
      
      //Use a fixed run number from year 2010, this year not available yet.
      TObjArray *maps = (TObjArray*)badmapContainer.GetObject(139000,"phosBadMap");
      if(!maps)
      {
        printf("AliCalorimeterUtils::SetOADBParameters() - Can not read PHOS bad map for run %d.\n",runnumber) ;    
      }
      else
      {
        printf("AliCalorimeterUtils::SetOADBParameters() - Setting PHOS bad map with name %s \n",maps->GetName()) ;
        for(Int_t mod=1; mod<5;mod++)
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
void AliCalorimeterUtils::AccessGeometry(AliVEvent* inputEvent) 
{
  //Set the calorimeters transformation matrices and init geometry
  
  // First init the geometry, a priory not done before
  Int_t runnumber = inputEvent->GetRunNumber() ;
  InitPHOSGeometry (runnumber);
  InitEMCALGeometry(runnumber);
  
  //Get the EMCAL transformation geometry matrices from ESD 
  if(!fEMCALGeoMatrixSet && fEMCALGeo)
  {
    if(fLoadEMCALMatrices)
    {
      printf("AliCalorimeterUtils::AccessGeometry() - Load user defined EMCAL geometry matrices\n");
      
      // OADB if available
      AliOADBContainer emcGeoMat("AliEMCALgeo");
      emcGeoMat.InitFromFile(Form("%s/EMCALlocal2master.root",fOADBFilePathEMCAL.Data()),"AliEMCALgeo");
      TObjArray *matEMCAL=(TObjArray*)emcGeoMat.GetObject(runnumber,"EmcalMatrices");
      
      for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
      {
        if (!fEMCALMatrix[mod]) // Get it from OADB
        {
          if(fDebug > 1 ) 
            printf("AliCalorimeterUtils::AccessGeometry() - EMCAL matrices SM %d, %p\n",
                   mod,((TGeoHMatrix*) matEMCAL->At(mod)));
          //((TGeoHMatrix*) matEMCAL->At(mod))->Print();
          
          fEMCALMatrix[mod] = (TGeoHMatrix*) matEMCAL->At(mod) ;
        }
        
        if(fEMCALMatrix[mod])
        {
          if(fDebug > 1) 
            fEMCALMatrix[mod]->Print();
          
          fEMCALGeo->SetMisalMatrix(fEMCALMatrix[mod],mod) ;  
        }
      }//SM loop
      
      fEMCALGeoMatrixSet = kTRUE;//At least one, so good
      
    }//Load matrices
    else if (!gGeoManager) 
    { 
      if(fDebug > 1) 
        printf(" AliCalorimeterUtils::AccessGeometry() - Load EMCAL misalignment matrices. \n");
      if(!strcmp(inputEvent->GetName(),"AliESDEvent"))  
      {
        for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
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
        if(fDebug > 1)
          printf("AliCalorimeterUtils::SetGeometryTransformationMatrices() - Setting of EMCAL transformation matrixes for AODs not implemented yet. \n Import geometry.root file\n");
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
      printf("AliCalorimeterUtils::SetGeometryTransformationMatrices() - Load user defined PHOS geometry matrices\n");
      
      // OADB if available
      AliOADBContainer geomContainer("phosGeo");
      geomContainer.InitFromFile(Form("%s/PHOSGeometry.root",fOADBFilePathPHOS.Data()),"PHOSRotationMatrixes");
      TObjArray *matPHOS = (TObjArray*)geomContainer.GetObject(139000,"PHOSRotationMatrixes");    
      
      for(Int_t mod = 0 ; mod < 5 ; mod++)
      {
        if (!fPHOSMatrix[mod]) // Get it from OADB
        {
          if(fDebug > 1 ) 
            printf("AliCalorimeterUtils::SetGeometryTransformationMatrices() - PHOS matrices module %d, %p\n",
                   mod,((TGeoHMatrix*) matPHOS->At(mod)));
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
      if(fDebug > 1) 
        printf(" AliCalorimeterUtils::SetGeometryTransformationMatrices() - Load PHOS misalignment matrices. \n");
			if(!strcmp(inputEvent->GetName(),"AliESDEvent"))  
      {
				for(Int_t mod = 0; mod < 5; mod++)
        { 
					if( ((AliESDEvent*)inputEvent)->GetPHOSMatrix(mod)) 
          {
						//printf("PHOS: mod %d, matrix %p\n",mod, ((AliESDEvent*)inputEvent)->GetPHOSMatrix(mod));
						fPHOSGeo->SetMisalMatrix( ((AliESDEvent*)inputEvent)->GetPHOSMatrix(mod),mod) ;
					}
				}// loop over modules
        fPHOSGeoMatrixSet  = kTRUE; //At least one so good
			}//ESD as input
			else 
      {
				if(fDebug > 1) 
					printf("AliCalorimeterUtils::SetGeometryTransformationMatrices() - Setting of EMCAL transformation matrixes for AODs not implemented yet. \n Import geometry.root file\n");
      }//AOD as input
    }// get matrix from data
    else if(gGeoManager)
    {
      fPHOSGeoMatrixSet = kTRUE;
    }
	}//PHOS geo	and  geoManager was not set
  
}

//______________________________________________________________________________________
Bool_t AliCalorimeterUtils::AreNeighbours(const TString calo, 
                                          const Int_t absId1, const Int_t absId2 ) const
{
  // Tells if (true) or not (false) two cells are neighbours
  // A neighbour is defined as being two cells which share a side or corner
	
  Bool_t areNeighbours = kFALSE ;
  
  Int_t iRCU1 = -1, irow1 = -1, icol1 = -1;
  Int_t iRCU2 = -1, irow2 = -1, icol2 = -1;
  
  Int_t rowdiff =  0, coldiff =  0;
  
  Int_t nSupMod1 = GetModuleNumberCellIndexes(absId1, calo, icol1, irow1, iRCU1); 
  Int_t nSupMod2 = GetModuleNumberCellIndexes(absId2, calo, icol2, irow2, iRCU2); 
  
  if(calo=="EMCAL" && nSupMod1!=nSupMod2)
  {
    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2-1
    // C Side impair SM, nSupMod%2=1; A side pair SM nSupMod%2=0
    if(nSupMod1%2) icol1+=AliEMCALGeoParams::fgkEMCALCols;
    else           icol2+=AliEMCALGeoParams::fgkEMCALCols;    
	}
  
  rowdiff = TMath::Abs( irow1 - irow2 ) ;  
  coldiff = TMath::Abs( icol1 - icol2 ) ;  
  
  if (( coldiff <= 1 )  && ( rowdiff <= 1 ) && (coldiff + rowdiff > 0)) 
    areNeighbours = kTRUE ;
  
  return areNeighbours;
  
}


//_____________________________________________________________________________________
Bool_t AliCalorimeterUtils::CheckCellFiducialRegion(AliVCluster* cluster, 
                                                    AliVCaloCells* cells, 
                                                    AliVEvent * event, Int_t iev) const 
{
  
	// Given the list of AbsId of the cluster, get the maximum cell and 
	// check if there are fNCellsFromBorder from the calorimeter border
	
  //If the distance to the border is 0 or negative just exit accept all clusters
	if(cells->GetType()==AliVCaloCells::kEMCALCell && fEMCALRecoUtils->GetNumberOfCellsFromEMCALBorder() <= 0 ) return kTRUE;
	if(cells->GetType()==AliVCaloCells::kPHOSCell  && fNCellsFromPHOSBorder  <= 0 ) return kTRUE;
  
  Int_t absIdMax	= -1;
	Float_t ampMax  = -1;
	
  AliMixedEvent * mixEvent = dynamic_cast<AliMixedEvent*> (event);
  Int_t nMixedEvents = 0 ; 
  Int_t * cellsCumul = NULL ;
  Int_t numberOfCells = 0 ;  
  if (mixEvent){
    nMixedEvents = mixEvent->GetNumberOfEvents() ; 
    if (cells->GetType()==AliVCaloCells::kEMCALCell) {
      cellsCumul =  mixEvent->GetEMCALCellsCumul() ; 
      numberOfCells = mixEvent->GetNumberOfEMCALCells() ;
    } 
    
    else if (cells->GetType()==AliVCaloCells::kPHOSCell) {
      cellsCumul =  mixEvent->GetPHOSCellsCumul() ; 
      numberOfCells = mixEvent->GetNumberOfPHOSCells() ;
    } 
    
    if(cellsCumul){
      
      Int_t startCell = cellsCumul[iev] ; 
      Int_t endCell   = (iev+1 < nMixedEvents)?cellsCumul[iev+1]:numberOfCells;
      //Find cells with maximum amplitude
      for(Int_t i = 0; i < cluster->GetNCells() ; i++){
        Int_t absId = cluster->GetCellAbsId(i) ;
        for (Int_t j = startCell; j < endCell ;  j++) {
          Short_t cellNumber;
          Int_t mclabel; 
          Double_t amp, time, efrac; 
          cells->GetCell(j, cellNumber, amp, time,mclabel,efrac) ; 
          if (absId == cellNumber) {
            if(amp > ampMax){
              ampMax   = amp;
              absIdMax = absId;
            }        
          }
        }
      }//loop on cluster cells
    }// cells cumul available
    else {
      printf("AliCalorimeterUtils::CheckCellFiducialRegion() - CellsCumul is NULL!!!\n");
      abort();
    }
  } else {//Normal SE Events
    for(Int_t i = 0; i < cluster->GetNCells() ; i++){
      Int_t absId = cluster->GetCellAbsId(i) ;
      Float_t amp	= cells->GetCellAmplitude(absId);
      if(amp > ampMax){
        ampMax   = amp;
        absIdMax = absId;
      }
    }
  }
	
	if(fDebug > 1)
		printf("AliCalorimeterUtils::CheckCellFiducialRegion() - Cluster Max AbsId %d, Cell Energy %2.2f, Cluster Energy %2.2f\n", 
           absIdMax, ampMax, cluster->E());
	
	if(absIdMax==-1) return kFALSE;
	
	//Check if the cell is close to the borders:
	Bool_t okrow = kFALSE;
	Bool_t okcol = kFALSE;
  
	if(cells->GetType()==AliVCaloCells::kEMCALCell){
    
		Int_t iTower = -1, iIphi = -1, iIeta = -1, iphi = -1, ieta = -1, iSM = -1; 
		fEMCALGeo->GetCellIndex(absIdMax,iSM,iTower,iIphi,iIeta); 
		fEMCALGeo->GetCellPhiEtaIndexInSModule(iSM,iTower,iIphi, iIeta,iphi,ieta);
		if(iSM < 0 || iphi < 0 || ieta < 0 ) {
      Fatal("CheckCellFidutialRegion","Negative value for super module: %d, or cell ieta: %d, or cell iphi: %d, check EMCAL geometry name\n",iSM,ieta,iphi);
    }
    
		//Check rows/phi
    Int_t nborder = fEMCALRecoUtils->GetNumberOfCellsFromEMCALBorder();
		if(iSM < 10)
    {
			if(iphi >= nborder && iphi < 24-nborder) okrow =kTRUE; 
    }
		else
    {
      if(fEMCALGeoName.Contains("12SM")) // 1/3 SM
      {
        if(iphi >= nborder && iphi < 8-nborder) okrow =kTRUE; 
      }
      else // 1/2 SM
      {
        if(iphi >= nborder && iphi <12-nborder) okrow =kTRUE; 
      }
		}
		
		//Check columns/eta
		if(!fEMCALRecoUtils->IsEMCALNoBorderAtEta0())
    {
			if(ieta  > nborder && ieta < 48-nborder) okcol =kTRUE; 
		}
		else
    {
			if(iSM%2==0)
      {
				if(ieta >= nborder)     okcol = kTRUE;	
			}
			else 
      {
				if(ieta <  48-nborder)  okcol = kTRUE;	
			}
		}//eta 0 not checked
    
		if(fDebug > 1)
		{
			printf("AliCalorimeterUtils::CheckCellFiducialRegion() - EMCAL Cluster in %d cells fiducial volume: ieta %d, iphi %d, SM %d ?",
             nborder, ieta, iphi, iSM);
			if (okcol && okrow ) printf(" YES \n");
			else  printf(" NO: column ok? %d, row ok? %d \n",okcol,okrow);
		}
	}//EMCAL
	else if(cells->GetType()==AliVCaloCells::kPHOSCell){
		Int_t relId[4];
		Int_t irow = -1, icol = -1;
		fPHOSGeo->AbsToRelNumbering(absIdMax,relId);
		irow = relId[2];
		icol = relId[3];
		//imod = relId[0]-1;
		if(irow >= fNCellsFromPHOSBorder && irow < 64-fNCellsFromPHOSBorder) okrow =kTRUE; 
		if(icol >= fNCellsFromPHOSBorder && icol < 56-fNCellsFromPHOSBorder) okcol =kTRUE; 
		if(fDebug > 1)
		{
			printf("AliCalorimeterUtils::CheckCellFiducialRegion() - PHOS Cluster in %d cells fiducial volume: icol %d, irow %d, Module %d?",
             fNCellsFromPHOSBorder, icol, irow, relId[0]-1);
			if (okcol && okrow ) printf(" YES \n");
			else  printf(" NO: column ok? %d, row ok? %d \n",okcol,okrow);
		}
	}//PHOS
	
	if (okcol && okrow) return kTRUE; 
	else                return kFALSE;
	
}	

//_________________________________________________________________________________________________________
Bool_t AliCalorimeterUtils::ClusterContainsBadChannel(TString calorimeter,UShort_t* cellList, Int_t nCells)
{
	// Check that in the cluster cells, there is no bad channel of those stored 
	// in fEMCALBadChannelMap or fPHOSBadChannelMap
	
	if (!fRemoveBadChannels) return kFALSE;
	//printf("fEMCALBadChannelMap %p, fPHOSBadChannelMap %p \n",fEMCALBadChannelMap,fPHOSBadChannelMap);
	if(calorimeter == "EMCAL" && !fEMCALRecoUtils->GetEMCALChannelStatusMap(0)) return kFALSE;
	if(calorimeter == "PHOS"  && !fPHOSBadChannelMap)  return kFALSE;
  
	Int_t icol = -1;
	Int_t irow = -1;
	Int_t imod = -1;
	for(Int_t iCell = 0; iCell<nCells; iCell++){
    
		//Get the column and row
		if(calorimeter == "EMCAL"){
      return fEMCALRecoUtils->ClusterContainsBadChannel((AliEMCALGeometry*)fEMCALGeo,cellList,nCells);
		}
		else if(calorimeter=="PHOS"){
			Int_t    relId[4];
			fPHOSGeo->AbsToRelNumbering(cellList[iCell],relId);
			irow = relId[2];
			icol = relId[3];
			imod = relId[0]-1;
			if(fPHOSBadChannelMap->GetEntries() <= imod)continue;
      //printf("PHOS bad channels imod %d, icol %d, irow %d\n",imod, irow, icol);
			if(GetPHOSChannelStatus(imod, irow, icol)) return kTRUE;
		}
		else return kFALSE;
		
	}// cell cluster loop
  
	return kFALSE;
  
}

//_______________________________________________________________
void AliCalorimeterUtils::CorrectClusterEnergy(AliVCluster *clus)
{
  // Correct cluster energy non linearity
  
  clus->SetE(fEMCALRecoUtils->CorrectClusterEnergyLinearity(clus));

}

//________________________________________________________________________________________
Int_t  AliCalorimeterUtils::GetMaxEnergyCell(AliVCaloCells* cells, const AliVCluster* clu, 
                                             Float_t & clusterFraction) const 
{
  
  //For a given CaloCluster gets the absId of the cell 
  //with maximum energy deposit.
  
  if( !clu || !cells ){
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
  
  TString           calo = "EMCAL";
  if(clu->IsPHOS()) calo = "PHOS";
  
  for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++) {
    
    cellAbsId = clu->GetCellAbsId(iDig);
    
    fraction  = clu->GetCellAmplitudeFraction(iDig);
    if(fraction < 1e-4) fraction = 1.; // in case unfolding is off
    
    iSupMod = GetModuleNumberCellIndexes(cellAbsId, calo, ieta, iphi, iRCU);
    
    if(IsRecalibrationOn()) {
      if(calo=="EMCAL") recalFactor = GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
      else              recalFactor = GetPHOSChannelRecalibrationFactor (iSupMod,iphi,ieta);
    }
    
    eCell  = cells->GetCellAmplitude(cellAbsId)*fraction*recalFactor;
    
    if(eCell > eMax)  { 
      eMax  = eCell; 
      absId = cellAbsId;
    }
    
    eTot+=eCell;
    
  }// cell loop
  
  clusterFraction = (eTot-eMax)/eTot; //Do not use cluster energy in case it was corrected for non linearity.
  
  return absId;
  
}

//__________________________________________________________________________
AliVTrack * AliCalorimeterUtils::GetMatchedTrack(const AliVCluster* cluster, 
                                                 const AliVEvent* event, 
                                                 const Int_t index) const
{
  // Get the matched track given its index, usually just the first match
  // Since it is different for ESDs and AODs here it is a wrap method to do it
  
  AliVTrack *track = 0;
  
  // EMCAL case only when matching is recalculated
  if(cluster->IsEMCAL() && IsRecalculationOfClusterTrackMatchingOn())
  {
    Int_t trackIndex = fEMCALRecoUtils->GetMatchedTrackIndex(cluster->GetID());
    //printf("track index %d, cluster ID %d \n ",trackIndex,cluster->GetID());
    
    if(trackIndex < 0 )
    { 
      printf("AliCalorimeterUtils::GetMatchedTrack() - Wrong track index %d, from recalculation\n", trackIndex);
    }
    else 
    {
      track = dynamic_cast<AliVTrack*> (event->GetTrack(trackIndex));
    }

    return track ;
    
  }   
  
  // Normal case, get info from ESD or AOD
  // ESDs
  if(!strcmp("AliESDCaloCluster",Form("%s",cluster->ClassName())))
  {
    Int_t iESDtrack = cluster->GetTrackMatchedIndex();
    
    if(iESDtrack < 0 )
    { 
      printf("AliCalorimeterUtils::GetMatchedTrack() - Wrong track index %d\n", index);
      return 0x0;
    }
    
    track = dynamic_cast<AliVTrack*> (event->GetTrack(iESDtrack));
    
  }
  else // AODs
  {
    if(cluster->GetNTracksMatched() > 0 )
      track = dynamic_cast<AliVTrack*>(cluster->GetTrackMatched(index));
  }
  
  return track ;
  
}

//_____________________________________________________________________________________________________
Int_t AliCalorimeterUtils::GetModuleNumber(AliAODPWG4Particle * particle, AliVEvent * inputEvent) const
{
	//Get the EMCAL/PHOS module number that corresponds to this particle
	
	Int_t absId = -1;
	if(particle->GetDetector()=="EMCAL"){
		fEMCALGeo->GetAbsCellIdFromEtaPhi(particle->Eta(),particle->Phi(), absId);
		if(fDebug > 2) 
		  printf("AliCalorimeterUtils::GetModuleNumber(PWG4AOD) - EMCAL: cluster eta %f, phi %f, absid %d, SuperModule %d\n",
             particle->Eta(), particle->Phi()*TMath::RadToDeg(),absId, fEMCALGeo->GetSuperModuleNumber(absId));
		return fEMCALGeo->GetSuperModuleNumber(absId) ;
	}//EMCAL
	else if(particle->GetDetector()=="PHOS")
  {
    // In case we use the MC reader, the input are TParticles, 
    // in this case use the corresponing method in PHOS Geometry to get the particle.
    if(strcmp(inputEvent->ClassName(), "AliMCEvent") == 0 )
    {
      Int_t mod =-1;
      Double_t z = 0., x=0.;
      TParticle* primary = 0x0;
      AliStack * stack = ((AliMCEvent*)inputEvent)->Stack();
      if(stack) {
        primary = stack->Particle(particle->GetCaloLabel(0));
      }
      else {
        Fatal("GetModuleNumber(PWG4AOD)", "Stack not available, stop!");
      }
      
      if(primary){
        fPHOSGeo->ImpactOnEmc(primary,mod,z,x) ;
      }
      else{
        Fatal("GetModuleNumber(PWG4AOD)", "Primary not available, stop!");
      }
      return mod;
    }
    // Input are ESDs or AODs, get the PHOS module number like this.
    else{
      //FIXME
      //AliVCluster *cluster = inputEvent->GetCaloCluster(particle->GetCaloLabel(0));
      //return GetModuleNumber(cluster);
      //MEFIX
      return -1;
    }
	}//PHOS
	
	return -1;
}

//_____________________________________________________________________
Int_t AliCalorimeterUtils::GetModuleNumber(AliVCluster * cluster) const
{
	//Get the EMCAL/PHOS module number that corresponds to this cluster
	TLorentzVector lv;
	Double_t v[]={0.,0.,0.}; //not necessary to pass the real vertex.
  if(!cluster)
  {
    if(fDebug > 1) printf("AliCalorimeterUtils::GetModuleNumber() - NUL Cluster, please check!!!");
    return -1;
  }
  
	cluster->GetMomentum(lv,v);
	Float_t phi = lv.Phi();
	if(phi < 0) phi+=TMath::TwoPi();	
	Int_t absId = -1;
	if(cluster->IsEMCAL()){
		fEMCALGeo->GetAbsCellIdFromEtaPhi(lv.Eta(),phi, absId);
		if(fDebug > 2) 
		  printf("AliCalorimeterUtils::GetModuleNumber() - EMCAL: cluster eta %f, phi %f, absid %d, SuperModule %d\n",
             lv.Eta(), phi*TMath::RadToDeg(),absId, fEMCALGeo->GetSuperModuleNumber(absId));
		return fEMCALGeo->GetSuperModuleNumber(absId) ;
	}//EMCAL
	else if(cluster->IsPHOS()) 
  {
		Int_t    relId[4];
		if ( cluster->GetNCells() > 0) 
    {
			absId = cluster->GetCellAbsId(0);
			if(fDebug > 2) 
				printf("AliCalorimeterUtils::GetModuleNumber() - PHOS: cluster eta %f, phi %f, e %f, absId %d\n",
               lv.Eta(), phi*TMath::RadToDeg(), lv.E(), absId);
		}
		else return -1;
		
		if ( absId >= 0) 
    {
			fPHOSGeo->AbsToRelNumbering(absId,relId);
			if(fDebug > 2) 
			  printf("AliCalorimeterUtils::GetModuleNumber() - PHOS: Module %d\n",relId[0]-1);
			return relId[0]-1;
		}
		else return -1;
	}//PHOS
	
	return -1;
}

//___________________________________________________________________________________________________
Int_t AliCalorimeterUtils::GetModuleNumberCellIndexes(const Int_t absId, const TString calo, 
                                                      Int_t & icol, Int_t & irow, Int_t & iRCU) const
{
	//Get the EMCAL/PHOS module, columns, row and RCU number that corresponds to this absId
	Int_t imod = -1;
	if ( absId >= 0) 
  {
		if(calo=="EMCAL")
    {
			Int_t iTower = -1, iIphi = -1, iIeta = -1; 
			fEMCALGeo->GetCellIndex(absId,imod,iTower,iIphi,iIeta); 
			fEMCALGeo->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,irow,icol);
      if(imod < 0 || irow < 0 || icol < 0 ) 
      {
        Fatal("GetModuleNumberCellIndexes()","Negative value for super module: %d, or cell icol: %d, or cell irow: %d, check EMCAL geometry name\n",imod,icol,irow);
      }
      
			//RCU0
      if(imod < 10 )
      {
        if      (0<=irow&&irow<8)                       iRCU=0; // first cable row
        else if (8<=irow&&irow<16 &&  0<=icol&&icol<24) iRCU=0; // first half; 
        //second cable row
        //RCU1
        else if (8<=irow&&irow<16 && 24<=icol&&icol<48) iRCU=1; // second half; 
        //second cable row
        else if (16<=irow&&irow<24)                     iRCU=1; // third cable row
        
        if (imod%2==1) iRCU = 1 - iRCU; // swap for odd=C side, to allow us to cable both sides the same
      }
      else 
      {
        // Last 2 SM have one single SRU, just assign RCU 0
        iRCU = 0 ;
      }

			if (iRCU<0) 
      {
				Fatal("GetModuleNumberCellIndexes()","Wrong EMCAL RCU number = %d\n", iRCU);
			}			
			
			return imod ;
      
		}//EMCAL
		else //PHOS
    {
			Int_t    relId[4];
			fPHOSGeo->AbsToRelNumbering(absId,relId);
			irow = relId[2];
			icol = relId[3];
			imod = relId[0]-1;
			iRCU= (Int_t)(relId[2]-1)/16 ;
			//Int_t iBranch= (Int_t)(relid[3]-1)/28 ; //0 to 1
			if (iRCU >= 4) 
      {
				Fatal("GetModuleNumberCellIndexes()","Wrong PHOS RCU number = %d\n", iRCU);
			}			
			return imod;
		}//PHOS	
	}
	
	return -1;
}

//___________________________________________________________________________________________
Int_t AliCalorimeterUtils::GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells) 
{
  // Find local maxima in cluster
  
  const Int_t   nc = cluster->GetNCells();
  Int_t   absIdList[nc]; 
  Float_t maxEList[nc]; 
  
  Int_t nMax = GetNumberOfLocalMaxima(cluster, cells, absIdList, maxEList);
  
  return nMax;
  
}

//___________________________________________________________________________________________
Int_t AliCalorimeterUtils::GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells,
                                                  Int_t *absIdList,     Float_t *maxEList) 
{
  // Find local maxima in cluster
  
  Int_t iDigitN = 0 ;
  Int_t iDigit  = 0 ;
  Int_t absId1 = -1 ;
  Int_t absId2 = -1 ;
  const Int_t nCells = cluster->GetNCells();
  
  TString calorimeter = "EMCAL";
  if(!cluster->IsEMCAL()) calorimeter = "PHOS";
  
  //printf("cluster : ncells %d \n",nCells);
  
  Float_t emax  = 0;
  Int_t   idmax =-1;
  for(iDigit = 0; iDigit < nCells ; iDigit++)
  {
    absIdList[iDigit] = cluster->GetCellsAbsId()[iDigit]  ; 
    Float_t en = cells->GetCellAmplitude(absIdList[iDigit]);
    RecalibrateCellAmplitude(en,calorimeter,absIdList[iDigit]);  
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
    if(absIdList[iDigit]>=0) 
    {
      absId1 = cluster->GetCellsAbsId()[iDigit];
      
      Float_t en1 = cells->GetCellAmplitude(absId1);
      RecalibrateCellAmplitude(en1,calorimeter,absId1);  
      
      //printf("%d : absIDi %d, E %f\n",iDigit, absId1,en1);
      
      for(iDigitN = 0; iDigitN < nCells; iDigitN++) 
      {	
        absId2 = cluster->GetCellsAbsId()[iDigitN] ;
        
        if(absId2==-1 || absId2==absId1) continue;
        
        //printf("\t %d : absIDj %d\n",iDigitN, absId2);
        
        Float_t en2 = cells->GetCellAmplitude(absId2);
        RecalibrateCellAmplitude(en2,calorimeter,absId2);
        
        //printf("\t %d : absIDj %d, E %f\n",iDigitN, absId2,en2);
        
        if ( AreNeighbours(calorimeter, absId1, absId2) ) 
        {
          // printf("\t \t Neighbours \n");
          if (en1 > en2 ) 
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
    if(absIdList[iDigit]>=0 )
    {
      absIdList[iDigitN] = absIdList[iDigit] ;
      Float_t en = cells->GetCellAmplitude(absIdList[iDigit]);
      RecalibrateCellAmplitude(en,calorimeter,absIdList[iDigit]);  
      if(en < fLocMaxCutE) continue; // Maxima only with seed energy at least
      maxEList[iDigitN] = en ;
      //printf("Local max %d, id %d, en %f\n", iDigit,absIdList[iDigitN],en);
      iDigitN++ ; 
    }
  }
  
  if(iDigitN == 0)
  {
    if(fDebug > 0) 
      printf("AliCalorimeterUtils::GetNumberOfLocalMaxima() - No local maxima found, assign highest energy cell as maxima, id %d, en cell %2.2f, en cluster %2.2f\n",
             idmax,emax,cluster->E());
    iDigitN      = 1     ;
    maxEList[0]  = emax  ;
    absIdList[0] = idmax ; 
  }
  
  if(fDebug > 0) 
  {    
    printf("AliCalorimeterUtils::GetNumberOfLocalMaxima() - In cluster E %2.2f, M02 %2.2f, M20 %2.2f, N maxima %d \n", 
           cluster->E(),cluster->GetM02(),cluster->GetM20(), iDigitN);
  
    if(fDebug > 1) for(Int_t imax = 0; imax < iDigitN; imax++) 
    {
      printf(" \t i %d, absId %d, Ecell %f\n",imax,absIdList[imax],maxEList[imax]);
    }
  }
  
  return iDigitN ;
  
}

//____________________________________
TString AliCalorimeterUtils::GetPass()
{
  // Get passx from filename.
    
  if (!AliAnalysisManager::GetAnalysisManager()->GetTree()) 
  {
    AliError("AliCalorimeterUtils::GetPass() - Pointer to tree = 0, returning null\n");
    return TString("");
  }
  
  if (!AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()) 
  {
    AliError("AliCalorimeterUtils::GetPass() - Null pointer input file, returning null\n");
    return TString("");
  }
  
  TString pass(AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()->GetName());
  if      (pass.Contains("ass1")) return TString("pass1");
  else if (pass.Contains("ass2")) return TString("pass2");
  else if (pass.Contains("ass3")) return TString("pass3");
  else if (pass.Contains("ass4")) return TString("pass4");
  else if (pass.Contains("ass5")) return TString("pass5");

  // No condition fullfilled, give a default value
  printf("AliCalorimeterUtils::GetPass() - Pass number string not found \n");
  return TString("");            
  
}

//________________________________________
void AliCalorimeterUtils::InitParameters()
{
  //Initialize the parameters of the analysis.
  
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
  
  fOADBFilePathEMCAL = "$ALICE_ROOT/OADB/EMCAL" ;          
  fOADBFilePathPHOS  = "$ALICE_ROOT/OADB/PHOS"  ;     
  
}


//_____________________________________________________
void AliCalorimeterUtils::InitPHOSBadChannelStatusMap()
{
  //Init PHOS bad channels map
  if(fDebug > 0 )printf("AliCalorimeterUtils::InitPHOSBadChannelStatusMap()\n");
  //In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fPHOSBadChannelMap = new TObjArray(5);	
  for (int i = 0; i < 5; i++)fPHOSBadChannelMap->Add(new TH2I(Form("PHOS_BadMap_mod%d",i),Form("PHOS_BadMap_mod%d",i), 64, 0, 64, 56, 0, 56));
  
  fPHOSBadChannelMap->SetOwner(kTRUE);
  fPHOSBadChannelMap->Compress();
  
  //In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);		
}

//______________________________________________________
void AliCalorimeterUtils::InitPHOSRecalibrationFactors()
{
	//Init EMCAL recalibration factors
	if(fDebug > 0 )printf("AliCalorimeterUtils::InitPHOSRecalibrationFactors()\n");
	//In order to avoid rewriting the same histograms
	Bool_t oldStatus = TH1::AddDirectoryStatus();
	TH1::AddDirectory(kFALSE);
  
	fPHOSRecalibrationFactors = new TObjArray(5);
	for (int i = 0; i < 5; i++)fPHOSRecalibrationFactors->Add(new TH2F(Form("PHOSRecalFactors_Mod%d",i),Form("PHOSRecalFactors_Mod%d",i), 64, 0, 64, 56, 0, 56));
	//Init the histograms with 1
	for (Int_t m = 0; m < 5; m++) {
		for (Int_t i = 0; i < 56; i++) {
			for (Int_t j = 0; j < 64; j++) {
				SetPHOSChannelRecalibrationFactor(m,j,i,1.);
			}
		}
	}
	fPHOSRecalibrationFactors->SetOwner(kTRUE);
	fPHOSRecalibrationFactors->Compress();
	
	//In order to avoid rewriting the same histograms
	TH1::AddDirectory(oldStatus);		
}


//__________________________________________________________
void AliCalorimeterUtils::InitEMCALGeometry(Int_t runnumber)
{
	//Initialize EMCAL geometry if it did not exist previously
  
	if (!fEMCALGeo)
  {
    if(fEMCALGeoName=="")
    {
      if     (runnumber <  140000 && 
              runnumber >= 100000)   fEMCALGeoName = "EMCAL_FIRSTYEARV1";
      else if(runnumber >= 140000 &&
              runnumber <  171000)   fEMCALGeoName = "EMCAL_COMPLETEV1";
      else                           fEMCALGeoName = "EMCAL_COMPLETE12SMV1";  
      printf("AliCalorimeterUtils::InitEMCALGeometry() - Set EMCAL geometry name to <%s> for run %d\n",fEMCALGeoName.Data(),runnumber);
    }
    
		fEMCALGeo = AliEMCALGeometry::GetInstance(fEMCALGeoName);
    
		if(fDebug > 0)
    {
			printf("AliCalorimeterUtils::InitEMCALGeometry(run=%d)",runnumber);
			if (!gGeoManager) printf(" - Careful!, gGeoManager not loaded, load misalign matrices");
			printf("\n");
		}
	}
}

//_________________________________________________________
void AliCalorimeterUtils::InitPHOSGeometry(Int_t runnumber)
{
	//Initialize PHOS geometry if it did not exist previously
  
	if (!fPHOSGeo)
  {
    if(fPHOSGeoName=="") fPHOSGeoName = "PHOSgeo";
      
		fPHOSGeo = new AliPHOSGeoUtils(fPHOSGeoName); 
    
		if(fDebug > 0)
    {
			printf("AliCalorimeterUtils::InitPHOSGeometry(run=%d)",runnumber);
			if (!gGeoManager) printf(" - Careful!, gGeoManager not loaded, load misalign matrices");
			printf("\n");
		}	
	}	
}

//__________________________________________________________________
Bool_t AliCalorimeterUtils::MaskFrameCluster(const Int_t iSM,  
                                             const Int_t ieta) const 
{
  //Check if cell is in one of the regions where we have significant amount 
  //of material in front. Only EMCAL
  
  Int_t icol = ieta;
  if(iSM%2) icol+=48; // Impair SM, shift index [0-47] to [48-96]
  
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
void AliCalorimeterUtils::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
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
  
  printf("Loc. Max. E > %2.2f\n",       fLocMaxCutE);
  printf("Loc. Max. E Diff > %2.2f\n",  fLocMaxCutEDiff);
  
  printf("    \n") ;
} 

//__________________________________________________________________________________________
void AliCalorimeterUtils::RecalibrateCellAmplitude(Float_t & amp,
                                                   const TString calo, const Int_t id) const
{
  //Recaculate cell energy if recalibration factor
  
  Int_t icol     = -1; Int_t irow     = -1; Int_t iRCU     = -1;
  Int_t nModule  = GetModuleNumberCellIndexes(id,calo, icol, irow, iRCU);
  
  if (IsRecalibrationOn()) 
  {
    if(calo == "PHOS") 
    {
      amp *= GetPHOSChannelRecalibrationFactor(nModule,icol,irow);
    }
    else		                   
    {
      amp *= GetEMCALChannelRecalibrationFactor(nModule,icol,irow);
    }
  }
}

//_________________________________________________________________________________
void AliCalorimeterUtils::RecalibrateCellTime(Double_t & time, 
                                              const TString calo, 
                                              const Int_t id, const Int_t bc) const
{
  // Recalculate time if time recalibration available for EMCAL
  // not ready for PHOS
  
  if(calo == "EMCAL" && GetEMCALRecoUtils()->IsTimeRecalibrationOn()) 
  {
    GetEMCALRecoUtils()->RecalibrateCellTime(id,bc,time);
  }
  
}

//__________________________________________________________________________
Float_t AliCalorimeterUtils::RecalibrateClusterEnergy(AliVCluster * cluster, 
                                                      AliVCaloCells * cells)
{
	// Recalibrate the cluster energy, considering the recalibration map and the energy of the cells that compose the cluster.
  
  //Initialize some used variables
	Float_t frac  = 0., energy = 0.;  
  
	if(cells) 
  {
    //Get the cluster number of cells and list of absId, check what kind of cluster do we have.
    
    UShort_t * index    = cluster->GetCellsAbsId() ;
    Double_t * fraction = cluster->GetCellsAmplitudeFraction() ;
    
    Int_t ncells     = cluster->GetNCells();	
    
    TString calo     = "EMCAL";
    if(cluster->IsPHOS()) 
      calo = "PHOS";
    
    //Loop on the cells, get the cell amplitude and recalibration factor, multiply and and to the new energy
    for(Int_t icell = 0; icell < ncells; icell++){
      
      Int_t absId = index[icell];
      
      frac =  fraction[icell];
      if(frac < 1e-3) frac = 1; //in case of EMCAL, this is set as 0, not used.
      
      Float_t amp = cells->GetCellAmplitude(absId);
      RecalibrateCellAmplitude(amp,calo, absId);
      
      if(fDebug>2)
        printf("AliCalorimeterUtils::RecalibrateClusterEnergy() - recalibrate cell: %s, cell fraction %f, cell energy %f\n", 
               calo.Data(),frac,cells->GetCellAmplitude(absId));
      
      energy += amp*frac;
    }
    
    if(fDebug>1)
      printf("AliCalorimeterUtils::RecalibrateClusterEnergy() - Energy before %f, after %f\n",cluster->E(),energy);
    
	}// cells available
  else
  {
    Fatal("RecalibrateClusterEnergy()","Cells pointer does not exist!");
  }
  
	return energy;
}

//__________________________________________________________________________________________
void AliCalorimeterUtils::RecalculateClusterPosition(AliVCaloCells* cells, AliVCluster* clu)
{
  
  //Recalculate EMCAL cluster position
  
  fEMCALRecoUtils->RecalculateClusterPosition((AliEMCALGeometry*)fEMCALGeo, cells,clu);
  
}

//________________________________________________________________________________
void AliCalorimeterUtils::RecalculateClusterTrackMatching(AliVEvent * event, 
                                                          TObjArray* clusterArray) 
{ 
  //Recalculate track matching
  
  if (fRecalculateMatching) 
  {
    fEMCALRecoUtils->FindMatches(event,clusterArray,fEMCALGeo)   ; 
    //AliESDEvent* esdevent = dynamic_cast<AliESDEvent*> (event);
    //if(esdevent){
    //  fEMCALRecoUtils->SetClusterMatchedToTrack(esdevent)  ;
    //  fEMCALRecoUtils->SetTracksMatchedToCluster(esdevent) ; 
    //}
  }
}

//___________________________________________________________________________
void AliCalorimeterUtils::SplitEnergy(const Int_t absId1, const Int_t absId2,
                                      AliVCluster* cluster, 
                                      AliVCaloCells* cells,
                                      //Float_t & e1, Float_t & e2,
                                      AliAODCaloCluster* cluster1,
                                      AliAODCaloCluster* cluster2,
                                      const Int_t nMax,
                                      const Int_t eventNumber)
{
  
  // Split energy of cluster between the 2 local maxima, sum energy on 3x3, and if the 2 
  // maxima are too close and have common cells, split the energy between the 2
  
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
  
  TString calorimeter = "EMCAL";
  if(cluster->IsPHOS())
  {
    calorimeter="PHOS";
    printf("AliCalorimeterUtils::SplitEnerg() Not supported for PHOS yet \n");
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
      if(icol > maxCol) maxCol = icol;
      if(icol < minCol) minCol = icol;
      if(irow > maxRow) maxRow = irow;
      if(irow < minRow) minRow = irow;
      hClusterMap->Fill(icol,irow,ec);
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
    
    if(absId==absId1 || absId==absId2 || absId < 0) continue;
    
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
  
  if(GetDebug() > 1) printf("AliCalorimeterUtils::SplitEnergy() - n Local Max %d, Cluster energy  = %f, Ecell1 = %f, Ecell2 = %f, Enew1 = %f, Enew2 = %f, Remain %f, \n ncells %d, ncells1 %d, ncells2 %d, f1 %f, f2  %f, sum f12 = %f \n",
                            nMax, eCluster,ecell1,ecell2,e1,e2,eCluster-e1-e2,ncells,ncells1,ncells2,shareFraction1,shareFraction2,shareFraction1+shareFraction2);
  
  cluster1->SetE(e1);
  cluster2->SetE(e2);  
  
  cluster1->SetNCells(ncells1);
  cluster2->SetNCells(ncells2);  
  
  cluster1->SetCellsAbsId(absIdList1);
  cluster2->SetCellsAbsId(absIdList2);
  
  cluster1->SetCellsAmplitudeFraction(fracList1);
  cluster2->SetCellsAmplitudeFraction(fracList2);
  
  //Correct linearity
  CorrectClusterEnergy(cluster1) ;
  CorrectClusterEnergy(cluster2) ;
  
  if(calorimeter=="EMCAL")
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
      
      if( AreNeighbours(calorimeter, absId2,absIdList1[iDigit]) )
        hCluster1->Fill(icol,irow,cells->GetCellAmplitude(absIdList1[iDigit])*shareFraction1);
      else 
        hCluster1->Fill(icol,irow,cells->GetCellAmplitude(absIdList1[iDigit]));
    }
    
    //printf(" \n ");
    //printf("Cells of cluster2: ");
    
    for(Int_t iDigit  = 0; iDigit < ncells2; iDigit++ ) 
    {
      //printf(" %d ",absIdList2[iDigit]);
      
      sm = GetModuleNumberCellIndexes(absIdList2[iDigit], calorimeter, icol, irow, iRCU) ;
      if( AreNeighbours(calorimeter, absId1,absIdList2[iDigit]) )
        hCluster2->Fill(icol,irow,cells->GetCellAmplitude(absIdList2[iDigit])*shareFraction2);
      else
        hCluster2->Fill(icol,irow,cells->GetCellAmplitude(absIdList2[iDigit]));
      
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
    hClusterMap    ->SetAxisRange(minCol, maxCol,"X");
    hClusterMap    ->SetAxisRange(minRow, maxRow,"Y");
    hClusterMap    ->Draw("colz");
    c->cd(2);
    gPad->SetGridy();
    gPad->SetGridx();
    hClusterLocMax ->SetAxisRange(minCol, maxCol,"X");
    hClusterLocMax ->SetAxisRange(minRow, maxRow,"Y");
    hClusterLocMax ->Draw("colz");
    c->cd(3);
    gPad->SetGridy();
    gPad->SetGridx();
    hCluster1      ->SetAxisRange(minCol, maxCol,"X");
    hCluster1      ->SetAxisRange(minRow, maxRow,"Y");
    hCluster1      ->Draw("colz");
    c->cd(4);
    gPad->SetGridy();
    gPad->SetGridx();
    hCluster2      ->SetAxisRange(minCol, maxCol,"X");
    hCluster2      ->SetAxisRange(minRow, maxRow,"Y");
    hCluster2      ->Draw("colz");
    
    if(eCluster > 6 )c->Print(Form("clusterFigures/Event%d_E%1.0f_nMax%d_NCell1_%d_NCell2_%d.eps",
                                   eventNumber,cluster->E(),nMax,ncells1,ncells2));
    
    delete c;
    delete hClusterMap;
    delete hClusterLocMax;
    delete hCluster1;
    delete hCluster2;
  }
}

