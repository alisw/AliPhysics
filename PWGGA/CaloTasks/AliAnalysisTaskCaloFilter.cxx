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

/* $Id: AliAnalysisTaskCaloFilter.cxx $ */

//////////////////////////////////////////////////////////
// Filter the ESDCaloClusters and ESDCaloCells of EMCAL,
// PHOS or both, creating the corresponing AODCaloClusters
// and AODCaloCells.
// Keep also the AODHeader information and the vertex.
// Needed for calorimeter calibration.
// Copy of AliAnalysisTaskESDfilter.
// Author: Gustavo Conesa Balbastre (INFN - Frascati)
//////////////////////////////////////////////////////////

//Root
#include "TGeoManager.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TInterpreter.h"

//STEER
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliLog.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliVEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDtrackCuts.h"
#include "AliTriggerAnalysis.h"

//EMCAL
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"

#include "AliAnalysisTaskCaloFilter.h"

ClassImp(AliAnalysisTaskCaloFilter)
  
////////////////////////////////////////////////////////////////////////

AliAnalysisTaskCaloFilter::AliAnalysisTaskCaloFilter():
  AliAnalysisTaskSE("CaloFilterTask"), //fCuts(0x0),
  fCaloFilter(0), fCorrect(kFALSE), 
  fEMCALGeo(0x0),fEMCALGeoName("EMCAL_FIRSTYEARV1"), 
  fEMCALRecoUtils(new AliEMCALRecoUtils),
  fESDtrackCuts(0), fTriggerAnalysis (new AliTriggerAnalysis), fTrackMultEtaCut(0.8),
  fLoadEMCALMatrices(kFALSE), //fLoadPHOSMatrices(kFALSE),
  fGeoMatrixSet(kFALSE), fEventNtuple(0),fConfigName(""),fFillAODFile(kTRUE)
{
  // Default constructor
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
  for(Int_t i = 0; i < 10; i++) fEMCALMatrix[i] = 0 ;
  //for(Int_t i = 0; i < 5 ; i++) fPHOSMatrix[i]  = 0 ;

  DefineOutput(1, TNtuple::Class());
  
}

//__________________________________________________
AliAnalysisTaskCaloFilter::AliAnalysisTaskCaloFilter(const char* name):
  AliAnalysisTaskSE(name), //fCuts(0x0),
  fCaloFilter(0), fCorrect(kFALSE),
  fEMCALGeo(0x0),fEMCALGeoName("EMCAL_FIRSTYEARV1"), 
  fEMCALRecoUtils(new AliEMCALRecoUtils),
  fESDtrackCuts(0), fTriggerAnalysis (new AliTriggerAnalysis), fTrackMultEtaCut(0.8),
  fLoadEMCALMatrices(kFALSE), //fLoadPHOSMatrices(kFALSE),
  fGeoMatrixSet(kFALSE), fEventNtuple(0),fConfigName(""),fFillAODFile(kTRUE)
{
  // Constructor
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
  for(Int_t i = 0; i < 10; i++) fEMCALMatrix[i] = 0 ;
  //for(Int_t i = 0; i < 5 ; i++) fPHOSMatrix[i]  = 0 ;

  DefineOutput(1, TNtuple::Class());

}

//__________________________________________________
AliAnalysisTaskCaloFilter::~AliAnalysisTaskCaloFilter()
{
  //Destructor.
	
  if(fEMCALGeo)        delete fEMCALGeo;	
  if(fEMCALRecoUtils)  delete fEMCALRecoUtils;
  if(fESDtrackCuts)    delete fESDtrackCuts;
  if(fTriggerAnalysis) delete fTriggerAnalysis;

  if(fEventNtuple)     delete fEventNtuple;
  
}

//-------------------------------------------------------------------
void AliAnalysisTaskCaloFilter::Init()
{
  //Init analysis with configuration macro if available
  
  if(gROOT->LoadMacro(fConfigName) >=0){
    printf("Configure analysis with %s\n",fConfigName.Data());
    AliAnalysisTaskCaloFilter *filter = (AliAnalysisTaskCaloFilter*)gInterpreter->ProcessLine("ConfigCaloFilter()");
    fEMCALGeoName      = filter->fEMCALGeoName; 
    fLoadEMCALMatrices = filter->fLoadEMCALMatrices;
    fFillAODFile       = filter->fFillAODFile;
    fEMCALRecoUtils    = filter->fEMCALRecoUtils; 
    fConfigName        = filter->fConfigName;
    fCaloFilter        = filter->fCaloFilter;
    fCorrect           = filter->fCorrect;
    fTrackMultEtaCut   = filter->fTrackMultEtaCut;
    fESDtrackCuts      = filter->fESDtrackCuts;
    for(Int_t i = 0; i < 10; i++) fEMCALMatrix[i] = filter->fEMCALMatrix[i] ;
  }
} 

//__________________________________________________
void AliAnalysisTaskCaloFilter::UserCreateOutputObjects()
{
  // Init geometry 
	
  fEMCALGeo =  AliEMCALGeometry::GetInstance(fEMCALGeoName) ;	
  
  OpenFile(1);
  
  fEventNtuple = new TNtuple("EventSelection","EventSelection", "bPileup:bGoodVertex:bV0AND:trackMult");
 
  PostData(1, fEventNtuple);

}  

//__________________________________________________
void AliAnalysisTaskCaloFilter::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event
  // Copy input ESD or AOD header, vertex, CaloClusters and CaloCells to output AOD

  if (fDebug > 0)  
    printf("CaloFilter: Analysing event # %d\n", (Int_t)Entry());
  
  AliVEvent* event = InputEvent();
  if(!event) {
    printf("AliAnalysisTaskCaloFilter::UserExec - This event does not contain Input?");
    return;
  }

  //Magic line to write events to file
  AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(fFillAODFile);
   
  // cast event, depending on input we will have to use one or the other type of event
  AliAODEvent* aodevent = dynamic_cast<AliAODEvent*> (event);  
  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*> (event);
  
  //-------------------------------------------
  //Event selection parameters
  //-------------------------------------------
  //Is it a pileup event?
  Bool_t bPileup = event->IsPileupFromSPD(3, 0.8, 3., 2., 5.); //Default values, if not it does not compile
  //Bool_t bPileup = event->IsPileupFromSPD(); 
  //if(bPileup) return kFALSE;
  Int_t  trackMult    = 0;
  Bool_t bV0AND       = kFALSE;
  Bool_t bGoodVertex  = kFALSE;
  if(esdevent){
    //Get track multiplicity
    Int_t nTracks   = InputEvent()->GetNumberOfTracks() ;
    for (Int_t itrack =  0; itrack <  nTracks; itrack++) {////////////// track loop
      AliVTrack * track = (AliVTrack*)InputEvent()->GetTrack(itrack) ; // retrieve track from esd
      if(!fESDtrackCuts->AcceptTrack((AliESDtrack*)track)) continue;
      //Count the tracks in eta < 0.8
      if(TMath::Abs(track->Eta())< fTrackMultEtaCut) trackMult++;
    }  
    //V0AND?   
    bV0AND = fTriggerAnalysis->IsOfflineTriggerFired(esdevent, AliTriggerAnalysis::kV0AND);
    //if(!bV0AND) return kFALSE;
    //Well reconstructed vertex
    bGoodVertex = CheckForPrimaryVertex();
    //if(!bGoodVertex) return kFALSE;
    
  }//ESDs
  
  if(fDebug > 0)
    printf("AliAnalysisTaskCaloFilter::UserExec() - PileUp %d, Good Vertex %d, v0AND %d, Track Mult in |eta| < %2.1f = %d\n",
           bPileup,bGoodVertex,bV0AND, fTrackMultEtaCut, trackMult);
  
  //Put bools with event selection parameters in a TNtuple
  //Int_t eventSelection[] = {bPileup,bGoodVertex,bV0AND,trackMult};
  fEventNtuple->Fill(bPileup,bGoodVertex,bV0AND,trackMult);
  
  //--------------------------------------------------------------------
  //Set in AOD General Event Parameters, vertex, runnumber, trigger, etc 
  //-------------------------------------------------------------------
  
  // set arrays and pointers
  Float_t  posF[3]  ;
  Double_t pos[3]   ;
  Double_t covVtx[6];
  for (Int_t i = 0; i < 6; i++)  covVtx[i] = 0.;
      
  AliAODHeader* header = AODEvent()->GetHeader();
  
  header->SetRunNumber(event->GetRunNumber());
  
  if(esdevent){
    TTree* tree = fInputHandler->GetTree();
    if (tree) {
      TFile* file = tree->GetCurrentFile();
      if (file) header->SetESDFileName(file->GetName());
    }
  }
  else if (aodevent) header->SetESDFileName(aodevent->GetHeader()->GetESDFileName());
  
  header->SetBunchCrossNumber(event->GetBunchCrossNumber());
  header->SetOrbitNumber(event->GetOrbitNumber());
  header->SetPeriodNumber(event->GetPeriodNumber());
  header->SetEventType(event->GetEventType());
  
  //Centrality
  if(event->GetCentrality()){
    header->SetCentrality(new AliCentrality(*(event->GetCentrality())));
  }
  else{
    header->SetCentrality(0);
  }
  
  //Trigger  
  header->SetOfflineTrigger(fInputHandler->IsEventSelected()); // propagate the decision of the physics selection
  if      (esdevent) header->SetFiredTriggerClasses(esdevent->GetFiredTriggerClasses());
  else if (aodevent) header->SetFiredTriggerClasses(aodevent->GetFiredTriggerClasses());
  header->SetTriggerMask(event->GetTriggerMask()); 
  header->SetTriggerCluster(event->GetTriggerCluster());
  if(esdevent){
    header->SetL0TriggerInputs(esdevent->GetHeader()->GetL0TriggerInputs());    
    header->SetL1TriggerInputs(esdevent->GetHeader()->GetL1TriggerInputs());    
    header->SetL2TriggerInputs(esdevent->GetHeader()->GetL2TriggerInputs());    
  }
  else if (aodevent){
    header->SetL0TriggerInputs(aodevent->GetHeader()->GetL0TriggerInputs());    
    header->SetL1TriggerInputs(aodevent->GetHeader()->GetL1TriggerInputs());    
    header->SetL2TriggerInputs(aodevent->GetHeader()->GetL2TriggerInputs());    
  }

  header->SetMagneticField(event->GetMagneticField());
  //header->SetMuonMagFieldScale(esdevent->GetCurrentDip()/6000.); 

  header->SetZDCN1Energy(event->GetZDCN1Energy());
  header->SetZDCP1Energy(event->GetZDCP1Energy());
  header->SetZDCN2Energy(event->GetZDCN2Energy());
  header->SetZDCP2Energy(event->GetZDCP2Energy());
  header->SetZDCEMEnergy(event->GetZDCEMEnergy(0),event->GetZDCEMEnergy(1));
  
  Float_t diamxy[2]={event->GetDiamondX(),event->GetDiamondY()};
  Float_t diamcov[3];
  event->GetDiamondCovXY(diamcov);
  header->SetDiamond(diamxy,diamcov);
  if      (esdevent) header->SetDiamondZ(esdevent->GetDiamondZ(),esdevent->GetSigma2DiamondZ());
  else if (aodevent) header->SetDiamondZ(aodevent->GetDiamondZ(),aodevent->GetSigma2DiamondZ());
  
  //
  //
  Int_t nVertices = 1 ;/* = prim. vtx*/;
  Int_t nCaloClus = event->GetNumberOfCaloClusters();
  
  AODEvent()->ResetStd(0, nVertices, 0, 0, 0, nCaloClus, 0, 0);
  
  // Access to the AOD container of vertices
  TClonesArray &vertices = *(AODEvent()->GetVertices());
  Int_t jVertices=0;
  
  // Add primary vertex. The primary tracks will be defined
  // after the loops on the composite objects (V0, cascades, kinks)
  event->GetPrimaryVertex()->GetXYZ(pos);
  Float_t chi = 0;
  if      (esdevent){
    esdevent->GetPrimaryVertex()->GetCovMatrix(covVtx);
    chi = esdevent->GetPrimaryVertex()->GetChi2toNDF();
  }
  else if (aodevent){
    aodevent->GetPrimaryVertex()->GetCovMatrix(covVtx);
    chi = aodevent->GetPrimaryVertex()->GetChi2perNDF();//Different from ESD?
  }
  
  AliAODVertex * primary = new(vertices[jVertices++])
    AliAODVertex(pos, covVtx, chi, NULL, -1, AliAODVertex::kPrimary);
  primary->SetName(event->GetPrimaryVertex()->GetName());
  primary->SetTitle(event->GetPrimaryVertex()->GetTitle());
  
  // Access to the AOD container of clusters
  TClonesArray &caloClusters = *(AODEvent()->GetCaloClusters());
  Int_t jClusters=0;
  
  //-------------------------------------------
  //Do Corrections in EMCAL 
  //-------------------------------------------
  //If EMCAL, and requested, correct energy, position ...
  //Need to do this in a separate loop before filling the ESDs because of the track matching recalculations
  if(fCorrect && (fCaloFilter==kEMCAL || fCaloFilter==kBoth) ) {
    
    if(!fGeoMatrixSet){
      if(fLoadEMCALMatrices){
        printf("AliAnalysisTaskCaloFilter::UserExec() - Load user defined EMCAL geometry matrices\n");
        for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
          if(fEMCALMatrix[mod]){
            if(DebugLevel() > 1) 
              fEMCALMatrix[mod]->Print();
            fEMCALGeo->SetMisalMatrix(fEMCALMatrix[mod],mod) ;  
          }
          fGeoMatrixSet=kTRUE;
        }//SM loop
      }//Load matrices
      else if(!gGeoManager){
        printf("AliAnalysisTaskCaloFilter::UserExec() - Get geo matrices from data\n");
        //Still not implemented in AOD, just a workaround to be able to work at least with ESDs	
        if(!strcmp(event->GetName(),"AliAODEvent")) {
          if(DebugLevel() > 1) 
            printf("AliAnalysisTaskCaloFilter Use ideal geometry, values geometry matrix not kept in AODs.\n");
        }//AOD
        else {	
          if(DebugLevel() > 1) printf("AliAnalysisTaskCaloFilter Load Misaligned matrices. \n");
          AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event) ;
          if(!esd) {
            printf("AliAnalysisTaskCaloFilter::UserExec() - This event does not contain ESDs?");
            return;
          }
          for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
            //if(DebugLevel() > 1) 
            esd->GetEMCALMatrix(mod)->Print();
            if(esd->GetEMCALMatrix(mod)) fEMCALGeo->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
          } 
          fGeoMatrixSet=kTRUE;
        }//ESD
      }//Load matrices from Data 
      
      //Recover time dependent corrections, put then in recalibration histograms. Do it once
      fEMCALRecoUtils->SetRunDependentCorrections(InputEvent()->GetRunNumber());

    }//first event
    
    
    //Cluster Loop
    for (Int_t iClust=0; iClust<nCaloClus; ++iClust) {
      
      AliVCluster * cluster = event->GetCaloCluster(iClust);
      if(cluster->IsPHOS()) continue ;

      Float_t position[]={0,0,0};
      if(DebugLevel() > 2)
        printf("Check cluster %d for bad channels and close to border\n",cluster->GetID());
      if(fEMCALRecoUtils->ClusterContainsBadChannel(fEMCALGeo,cluster->GetCellsAbsId(), cluster->GetNCells())) continue;	
      //      if(!fEMCALRecoUtils->CheckCellFiducialRegion(fEMCALGeo, cluster, event->GetEMCALCells())) {
      //        printf("Finally reject\n");
      //        continue;
      //      }
      if(DebugLevel() > 2)
      { 
        printf("Filter, before  : i %d, E %f, dispersion %f, m02 %f, m20 %f, distToBad %f\n",iClust,cluster->E(),
               cluster->GetDispersion(),cluster->GetM02(),cluster->GetM20(), cluster->GetDistanceToBadChannel());
        cluster->GetPosition(position);
        printf("Filter, before  : i %d, x %f, y %f, z %f\n",cluster->GetID(), position[0], position[1], position[2]);
      }
      
      //Recalculate distance to bad channels, if new list of bad channels provided
      fEMCALRecoUtils->RecalculateClusterDistanceToBadChannel(fEMCALGeo, event->GetEMCALCells(), cluster);

      if(fEMCALRecoUtils->IsRecalibrationOn())	{
        fEMCALRecoUtils->RecalibrateClusterEnergy(fEMCALGeo, cluster, event->GetEMCALCells());
        fEMCALRecoUtils->RecalculateClusterShowerShapeParameters(fEMCALGeo, event->GetEMCALCells(),cluster);
        fEMCALRecoUtils->RecalculateClusterPID(cluster);
      }
            
      fEMCALRecoUtils->RecalculateClusterPosition(fEMCALGeo, event->GetEMCALCells(),cluster);
      
      if(DebugLevel() > 2)
      { 
        printf("Filter, after   : i %d, E %f, dispersion %f, m02 %f, m20 %f, distToBad %f\n",cluster->GetID(),cluster->E(),
               cluster->GetDispersion(),cluster->GetM02(),cluster->GetM20(), cluster->GetDistanceToBadChannel());
        cluster->GetPosition(position);
        printf("Filter, after   : i %d, x %f, y %f, z %f\n",cluster->GetID(), position[0], position[1], position[2]);
      }    
      
      cluster->SetE(fEMCALRecoUtils->CorrectClusterEnergyLinearity(cluster));
      
    }
    //Recalculate track-matching
    fEMCALRecoUtils->FindMatches(event,0,fEMCALGeo);
    
  } // corrections in EMCAL
  
  //-------------------------------------------
  // Now loop on clusters to fill AODs
  //-------------------------------------------
  for (Int_t iClust=0; iClust<nCaloClus; ++iClust) {
    
    AliVCluster * cluster = event->GetCaloCluster(iClust);
    
    //Check which calorimeter information we want to keep.
    
    if(fCaloFilter!=kBoth){
      if     (fCaloFilter==kPHOS  && cluster->IsEMCAL()) continue ;
      else if(fCaloFilter==kEMCAL && cluster->IsPHOS())  continue ;
    }  
    
    //Temporary trick, FIXME
    Float_t dR = cluster->GetTrackDx();
    Float_t dZ = cluster->GetTrackDz();
    if(DebugLevel() > 2)
      printf("Original residuals : dZ %f, dR %f\n ",dZ, dR);
    //--------------------------------------------------------------
    //If EMCAL and corrections done, get the new matching parameters, do not copy noisy clusters
    if(cluster->IsEMCAL() && fCorrect && esdevent){
      if(DebugLevel() > 2)
        printf("Check cluster %d for bad channels and close to border\n",cluster->GetID());
      if(fEMCALRecoUtils->ClusterContainsBadChannel(fEMCALGeo,cluster->GetCellsAbsId(), cluster->GetNCells())) continue;	
      
      //      if(!fEMCALRecoUtils->CheckCellFiducialRegion(fEMCALGeo, cluster, event->GetEMCALCells())) {
      //        printf("Finally reject\n");
      //        continue;
      //      }

      if(fEMCALRecoUtils->IsExoticCluster(cluster, InputEvent()->GetEMCALCells(),InputEvent()->GetBunchCrossNumber())) continue;	

      fEMCALRecoUtils->GetMatchedResiduals(cluster->GetID(),dR,dZ);
      cluster->SetTrackDistance(dR,dZ);
    }
    
    if(DebugLevel() > 2){
      if(cluster->IsEMCAL()) printf("EMCAL Track-Cluster Residuals : dZ %f, dR %f\n ",dZ, dR);
      if(cluster->IsPHOS())  printf("PHOS  Track-Cluster Residuals : dZ %f, dR %f\n ",dZ, dR);
    }
    
    //--------------------------------------------------------------

    //Now fill AODs
    
    Int_t id       = cluster->GetID();
    Float_t energy = cluster->E();
    cluster->GetPosition(posF);
    
    AliAODCaloCluster *caloCluster = new(caloClusters[jClusters++]) 
      AliAODCaloCluster(id,
			cluster->GetNLabels(),
			cluster->GetLabels(),
			energy,
			posF,
			NULL,
			cluster->GetType());
    
    caloCluster->SetChi2(dZ);//Temporary trick, FIXME
    caloCluster->SetCaloCluster(cluster->GetDistanceToBadChannel(),
				cluster->GetDispersion(),
				cluster->GetM20(), cluster->GetM02(),
				dR,  //Temporary trick, FIXME
				cluster->GetNExMax(),cluster->GetTOF()) ;
    
    caloCluster->SetPIDFromESD(cluster->GetPID());
    caloCluster->SetNCells(cluster->GetNCells());
    caloCluster->SetCellsAbsId(cluster->GetCellsAbsId());
    caloCluster->SetCellsAmplitudeFraction(cluster->GetCellsAmplitudeFraction());
    
    if(DebugLevel() > 2)
    { 
      printf("Filter, aod     : i %d, E %f, dispersion %f, m02 %f, m20 %f\n",caloCluster->GetID(),caloCluster->E(),
             caloCluster->GetDispersion(),caloCluster->GetM02(),caloCluster->GetM20());
      caloCluster->GetPosition(posF);
      printf("Filter, aod     : i %d, x %f, y %f, z %f\n",caloCluster->GetID(), posF[0], posF[1], posF[2]);
    }    
    
    //Matched tracks, just to know if there was any match, the track pointer is useless.
    //Temporary trick, FIXME
    if(esdevent){
      if(TMath::Abs(dR) < 990 && TMath::Abs(dZ) < 990) { //Default value in PHOS 999, in EMCAL 1024, why?
        if(DebugLevel() > 2) 
          printf("*** Cluster Track-Matched *** dR %f, dZ %f\n",caloCluster->GetEmcCpvDistance(),caloCluster->Chi2());
        caloCluster->AddTrackMatched(new AliAODTrack); 
      }// fill the array with one entry to signal a possible match
      //TArrayI* matchedT = 	((AliESDCaloCluster*)cluster)->GetTracksMatched();
      //if (InputEvent()->GetNumberOfTracks() > 0 && matchedT && cluster->GetTrackMatchedIndex() >= 0) {	
      //  for (Int_t im = 0; im < matchedT->GetSize(); im++) {
      //    Int_t iESDtrack = matchedT->At(im);;
      //    if ((AliVTrack*)InputEvent()->GetTrack(iESDtrack) != 0) {
      //      caloCluster->AddTrackMatched((AliVTrack*)InputEvent()->GetTrack(iESDtrack));
      //    }
      //  }
      //}// There is at least a match with a track
    }
  } 
  caloClusters.Expand(jClusters); // resize TObjArray to 'remove' slots for pseudo clusters	 
  // end of loop on calo clusters
  
  // fill EMCAL cell info
  if ((fCaloFilter==kBoth ||  fCaloFilter==kEMCAL) && event->GetEMCALCells()) { // protection against missing ESD information
    AliVCaloCells &eventEMcells = *(event->GetEMCALCells());
    Int_t nEMcell = eventEMcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
    aodEMcells.CreateContainer(nEMcell);
    aodEMcells.SetType(AliVCaloCells::kEMCALCell);
    Double_t calibFactor = 1.;   
    for (Int_t iCell = 0; iCell < nEMcell; iCell++) { 
      Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
      fEMCALGeo->GetCellIndex(eventEMcells.GetCellNumber(iCell),imod,iTower,iIphi,iIeta); 
      fEMCALGeo->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);	
      
      if(fCorrect && fEMCALRecoUtils->IsRecalibrationOn()){ 
        calibFactor = fEMCALRecoUtils->GetEMCALChannelRecalibrationFactor(imod,ieta,iphi);
      }
      
      if(!fEMCALRecoUtils->GetEMCALChannelStatus(imod, ieta, iphi)){ //Channel is not declared as bad
        aodEMcells.SetCell(iCell,eventEMcells.GetCellNumber(iCell),eventEMcells.GetAmplitude(iCell)*calibFactor);
        //printf("GOOD channel\n");
      }
      else {
        aodEMcells.SetCell(iCell,eventEMcells.GetCellNumber(iCell),0);
        //printf("BAD channel\n");
      }
    }
    aodEMcells.Sort();
  }
  
  // fill PHOS cell info
  if ((fCaloFilter==kBoth ||  fCaloFilter==kPHOS) && event->GetPHOSCells()) { // protection against missing ESD information
    AliVCaloCells &eventPHcells = *(event->GetPHOSCells());
    Int_t nPHcell = eventPHcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodPHcells = *(AODEvent()->GetPHOSCells());
    aodPHcells.CreateContainer(nPHcell);
    aodPHcells.SetType(AliVCaloCells::kPHOSCell);
    for (Int_t iCell = 0; iCell < nPHcell; iCell++) {      
      aodPHcells.SetCell(iCell,eventPHcells.GetCellNumber(iCell),eventPHcells.GetAmplitude(iCell));
    }
    aodPHcells.Sort();
  }
  
   PostData(1, fEventNtuple);
  
  //return;
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskCaloFilter::CheckForPrimaryVertex(){
  //Check if the vertex was well reconstructed, copy from V0Reader of conversion group
  //It only works for ESDs
  
  AliESDEvent * event = dynamic_cast<AliESDEvent*> (InputEvent());
  if(!event) return kFALSE;
  
  if(event->GetPrimaryVertexTracks()->GetNContributors() > 0) {
    return kTRUE;
  }
  
  if(event->GetPrimaryVertexTracks()->GetNContributors() < 1) {
    // SPD vertex
    if(event->GetPrimaryVertexSPD()->GetNContributors() > 0) {
      //cout<<"spd vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
      return kTRUE;
      
    }
    if(event->GetPrimaryVertexSPD()->GetNContributors() < 1) {
      //      cout<<"bad vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
      return kFALSE;
    }
  }
  return kFALSE;

}


//_____________________________________________________
void AliAnalysisTaskCaloFilter::PrintInfo(){

  //Print settings

  printf("TASK: AnalysisCaloFilter \n");
  printf("\t Not only filter, correct Clusters? %d\n",fCorrect);
  printf("\t Calorimeter Filtering Option     ? %d\n",fCaloFilter);
  //printf("\t Use handmade geo matrices?   EMCAL %d, PHOS %d\n",fLoadEMCALMatrices, fLoadPHOSMatrices);
  printf("\t Use handmade geo matrices?   EMCAL %d, PHOS 0\n",fLoadEMCALMatrices);
}

//_____________________________________________________
//void AliAnalysisTaskCaloFilter::LocalInit()
//{
//	// Local Initialization
//	
//	// Create cuts/param objects and publish to slot
//	const Int_t buffersize = 255;
//	char onePar[buffersize] ;
//	fCuts = new TList();
//  
//	snprintf(onePar,buffersize, "Calorimeter Filtering Option %d", fCaloFilter) ;
//	fCuts->Add(new TObjString(onePar));
//	snprintf(onePar,buffersize, "Not only filter but correct? %d cells;", fCorrect) ;
//	fCuts->Add(new TObjString(onePar));
//  
//	// Post Data
//	PostData(2, fCuts);
//	
//}


//__________________________________________________
void AliAnalysisTaskCaloFilter::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
    if (fDebug > 1) printf("AnalysisCaloFilter: Terminate() \n");
}

