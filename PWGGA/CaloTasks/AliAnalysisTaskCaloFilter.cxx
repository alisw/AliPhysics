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

//////////////////////////////////////////////////////////
// Filter the ESDCaloClusters and ESDCaloCells of EMCAL,
// PHOS or both, creating the corresponing AODCaloClusters
// and AODCaloCells.
// Keep also the AODHeader information and the vertex.
// Keep tracks if requested
// Copy of AliAnalysisTaskESDfilter.
// Author: Gustavo Conesa Balbastre (INFN - Frascati)
//////////////////////////////////////////////////////////

//Root
#include "TGeoManager.h"
#include "TFile.h"
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
fLoadEMCALMatrices(kFALSE), //fLoadPHOSMatrices(kFALSE),
fGeoMatrixSet(kFALSE),
fConfigName(""),fFillAODFile(kTRUE), fFillTracks(kFALSE),
fEnergyCut(0.), fNcellsCut (0), fVzCut(100.)
{
  // Default constructor
  
  for(Int_t i = 0; i < 12; i++) fEMCALMatrix[i] = 0 ;
  //for(Int_t i = 0; i < 5 ; i++) fPHOSMatrix[i]  = 0 ;
  
}

//__________________________________________________
AliAnalysisTaskCaloFilter::AliAnalysisTaskCaloFilter(const char* name):
AliAnalysisTaskSE(name), //fCuts(0x0),
fCaloFilter(0), fCorrect(kFALSE),
fEMCALGeo(0x0),fEMCALGeoName("EMCAL_FIRSTYEARV1"), 
fEMCALRecoUtils(new AliEMCALRecoUtils),
fLoadEMCALMatrices(kFALSE), //fLoadPHOSMatrices(kFALSE),
fGeoMatrixSet(kFALSE),
fConfigName(""),fFillAODFile(kTRUE),fFillTracks(kFALSE),
fEnergyCut(0.), fNcellsCut (0), fVzCut(100.)
{
  // Constructor
  
  for(Int_t i = 0; i < 12; i++) fEMCALMatrix[i] = 0 ;
  //for(Int_t i = 0; i < 5 ; i++) fPHOSMatrix[i]  = 0 ;
  
}

//__________________________________________________
AliAnalysisTaskCaloFilter::~AliAnalysisTaskCaloFilter()
{
  //Destructor.
	
  if(fEMCALGeo)        delete fEMCALGeo;	
  if(fEMCALRecoUtils)  delete fEMCALRecoUtils;
  
}

//-------------------------------------------------------------------
void AliAnalysisTaskCaloFilter::Init()
{
  //Init analysis with configuration macro if available
  
  if(gROOT->LoadMacro(fConfigName) >=0)
  {
    printf("Configure analysis with %s\n",fConfigName.Data());
    
    AliAnalysisTaskCaloFilter *filter = (AliAnalysisTaskCaloFilter*)gInterpreter->ProcessLine("ConfigCaloFilter()");
    
    fEMCALGeoName      = filter->fEMCALGeoName; 
    fLoadEMCALMatrices = filter->fLoadEMCALMatrices;
    fFillAODFile       = filter->fFillAODFile;
    fFillTracks        = filter->fFillTracks;
    fEMCALRecoUtils    = filter->fEMCALRecoUtils; 
    fConfigName        = filter->fConfigName;
    fCaloFilter        = filter->fCaloFilter;
    fCorrect           = filter->fCorrect;
    fEnergyCut         = filter->fEnergyCut;
    fNcellsCut         = filter->fNcellsCut;
    fVzCut             = filter->fVzCut;
    
    for(Int_t i = 0; i < 12; i++) fEMCALMatrix[i] = filter->fEMCALMatrix[i] ;
  }
} 

//__________________________________________________
void AliAnalysisTaskCaloFilter::UserCreateOutputObjects()
{
  // Init geometry 
	
  fEMCALGeo =  AliEMCALGeometry::GetInstance(fEMCALGeoName) ;	
  
}  

//____________________________________________________________
void AliAnalysisTaskCaloFilter::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event
  // Copy input ESD or AOD header, vertex, CaloClusters and CaloCells to output AOD
  
  if (fDebug > 0)  
    printf("CaloFilter: Analysing event # %d\n", (Int_t)Entry());
  
  AliVEvent* event = InputEvent();
  if(!event) 
  {
    printf("AliAnalysisTaskCaloFilter::UserExec - This event does not contain Input?");
    return;
  }
  
  // Select the event
  
  if(!AcceptEventVertex()) return;
  if(!AcceptEventEMCAL())  return;
  
  //Magic line to write events to file
  AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(fFillAODFile);
  
  Int_t nVertices = event->GetNumberOfV0s()+event->GetNumberOfCascades();;
  Int_t nCaloClus = event->GetNumberOfCaloClusters();
  Int_t nTracks   = event->GetNumberOfTracks();
  
  AODEvent()->ResetStd(nTracks, nVertices, 0, 0, 0, nCaloClus, 0, 0);
  
  //
  FillAODHeader();
  
  //
  FillAODVertices();
  
  //
  FillAODTracks();
  
  //
  CorrectionsInEMCAL();
  
  //
  FillAODCaloClusters();
  
  //
  FillAODCaloCells();
  
}

//___________________________________________________
Bool_t AliAnalysisTaskCaloFilter::AcceptEventVertex()
{
  // Accept event with good vertex
  
  Double_t v[3];
  InputEvent()->GetPrimaryVertex()->GetXYZ(v) ;
  
  if(TMath::Abs(v[2]) > fVzCut) return kFALSE ;
  
  return CheckForPrimaryVertex();
}  

//__________________________________________________
Bool_t AliAnalysisTaskCaloFilter::AcceptEventEMCAL()
{
  // Accept event given there is a cluster with enough energy
  
  if(fCaloFilter!=kEMCAL) return kTRUE; // accept all
  
  Int_t           nCluster = InputEvent() -> GetNumberOfCaloClusters();
  AliVCaloCells * caloCell = InputEvent() -> GetEMCALCells();
  Int_t           bc       = InputEvent() -> GetBunchCrossNumber();
  
  for(Int_t icalo = 0; icalo < nCluster; icalo++)
  {
    AliESDCaloCluster *clus = (AliESDCaloCluster*) (InputEvent()->GetCaloCluster(icalo));
    
    if( ( clus->IsEMCAL() ) && ( clus->GetNCells() > fNcellsCut ) && ( clus->E() > fEnergyCut ) &&
       fEMCALRecoUtils->IsGoodCluster(clus,fEMCALGeo,caloCell,bc))
    { 
      
      // printf("Accept event %d, E %2.2f > %2.2f, nCells %d > %d \n", 
      //        (Int_t) Entry(),clus->E(), fEnergyCut, clus->GetNCells(), fNcellsCut);
      
      return kTRUE;
    }
    
  }// loop
  
  return kFALSE;
  
}  

//________________________________________________
void AliAnalysisTaskCaloFilter::FillAODCaloCells()
{  
  // Fill EMCAL/PHOS cell info
  
  AliVEvent * event = InputEvent();
  
  // EMCAL
  if ((fCaloFilter==kBoth ||  fCaloFilter==kEMCAL) && event->GetEMCALCells()) 
  { // protection against missing ESD information
    AliVCaloCells &eventEMcells = *(event->GetEMCALCells());
    Int_t nEMcell = eventEMcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
    aodEMcells.CreateContainer(nEMcell);
    aodEMcells.SetType(AliVCaloCells::kEMCALCell);
    Double_t calibFactor = 1.;   
    for (Int_t iCell = 0; iCell < nEMcell; iCell++) 
    { 
      Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
      fEMCALGeo->GetCellIndex(eventEMcells.GetCellNumber(iCell),imod,iTower,iIphi,iIeta); 
      fEMCALGeo->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);	
      
      if(fCorrect && fEMCALRecoUtils->IsRecalibrationOn())
      { 
        calibFactor = fEMCALRecoUtils->GetEMCALChannelRecalibrationFactor(imod,ieta,iphi);
      }
      
      if(!fEMCALRecoUtils->GetEMCALChannelStatus(imod, ieta, iphi))
      { //Channel is not declared as bad
        aodEMcells.SetCell(iCell,eventEMcells.GetCellNumber(iCell),eventEMcells.GetAmplitude(iCell)*calibFactor,
                           eventEMcells.GetTime(iCell),eventEMcells.GetMCLabel(iCell),eventEMcells.GetEFraction(iCell));
        //printf("GOOD channel\n");
      }
      else 
      {
        aodEMcells.SetCell(iCell,eventEMcells.GetCellNumber(iCell),0,-1,-1,0);
        //printf("BAD channel\n");
      }
    }
    aodEMcells.Sort();
  }
  
  // PHOS 
  if ((fCaloFilter==kBoth ||  fCaloFilter==kPHOS) && event->GetPHOSCells()) 
  { // protection against missing ESD information
    AliVCaloCells &eventPHcells = *(event->GetPHOSCells());
    Int_t nPHcell = eventPHcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodPHcells = *(AODEvent()->GetPHOSCells());
    aodPHcells.CreateContainer(nPHcell);
    aodPHcells.SetType(AliVCaloCells::kPHOSCell);
    
    for (Int_t iCell = 0; iCell < nPHcell; iCell++) 
    {      
      aodPHcells.SetCell(iCell,eventPHcells.GetCellNumber(iCell),eventPHcells.GetAmplitude(iCell),
                         eventPHcells.GetTime(iCell),eventPHcells.GetMCLabel(iCell),eventPHcells.GetEFraction(iCell));
    }
    
    aodPHcells.Sort();
  }
}


//___________________________________________________
void AliAnalysisTaskCaloFilter::FillAODCaloClusters()
{
  // Fill the AOD with caloclusters
  
  AliVEvent * event = InputEvent();
  
  // Access to the AOD container of clusters
  
  TClonesArray &caloClusters = *(AODEvent()->GetCaloClusters());
  Int_t jClusters=0;
  Float_t  posF[3]  ;
  
  Int_t nCaloClus = event->GetNumberOfCaloClusters();
  for (Int_t iClust=0; iClust<nCaloClus; ++iClust) 
  {
    
    AliVCluster * cluster = event->GetCaloCluster(iClust);
    
    //Check which calorimeter information we want to keep.
    
    if(fCaloFilter!=kBoth)
    {
      if     (fCaloFilter==kPHOS  && cluster->IsEMCAL()) continue ;
      else if(fCaloFilter==kEMCAL && cluster->IsPHOS())  continue ;
    }  
    
    // Get original residuals, in case of previous recalculation, reset them 
    Float_t dR = cluster->GetTrackDx();
    Float_t dZ = cluster->GetTrackDz();
    
    if(DebugLevel() > 2)
      printf("Original residuals : dZ %f, dR %f\n ",dZ, dR);
    
    //--------------------------------------------------------------
    //If EMCAL and corrections done, get the new matching parameters, do not copy noisy clusters
    if(cluster->IsEMCAL() && fCorrect)
    {
      if(DebugLevel() > 2)
        printf("Check cluster %d for bad channels and close to border\n",cluster->GetID());
      
      if(fEMCALRecoUtils->ClusterContainsBadChannel(fEMCALGeo,cluster->GetCellsAbsId(), cluster->GetNCells())) continue;	
      
      if(fEMCALRecoUtils->IsExoticCluster(cluster, InputEvent()->GetEMCALCells(),InputEvent()->GetBunchCrossNumber())) continue;	
      
      fEMCALRecoUtils->GetMatchedResiduals(cluster->GetID(),dR,dZ);
      cluster->SetTrackDistance(dR,dZ);
    }
    
    if(DebugLevel() > 2)
    {
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
    
    caloCluster->SetCaloCluster(cluster->GetDistanceToBadChannel(),
                                cluster->GetDispersion(),
                                cluster->GetM20(), cluster->GetM02(),
                                -1,  
                                cluster->GetNExMax(),cluster->GetTOF()) ;
    
    caloCluster->SetPIDFromESD(cluster->GetPID());
    caloCluster->SetNCells(cluster->GetNCells());
    caloCluster->SetCellsAbsId(cluster->GetCellsAbsId());
    caloCluster->SetCellsAmplitudeFraction(cluster->GetCellsAmplitudeFraction());
    caloCluster->SetTrackDistance(dR, dZ);
    
    if(DebugLevel() > 2)
    { 
      printf("Filter, aod     : i %d, E %f, dispersion %f, m02 %f, m20 %f\n",caloCluster->GetID(),caloCluster->E(),
             caloCluster->GetDispersion(),caloCluster->GetM02(),caloCluster->GetM20());
      caloCluster->GetPosition(posF);
      printf("Filter, aod     : i %d, x %f, y %f, z %f\n",caloCluster->GetID(), posF[0], posF[1], posF[2]);
    }    
    
    //Matched tracks, just to know if there was any match, the track pointer is useless if tracks not stored
    if(TMath::Abs(dR) < 990 && TMath::Abs(dZ) < 990) 
    { //Default value in PHOS 999, in EMCAL 1024, why?
      caloCluster->AddTrackMatched(new AliAODTrack); 
    }
    // TO DO, in case Tracks available, think how to put the matched track in AOD
  }
  
  caloClusters.Expand(jClusters); // resize TObjArray to 'remove' slots  
  
}


//______________________________________________
void AliAnalysisTaskCaloFilter::FillAODHeader()
{
  // AOD header copy
  
  AliVEvent* event = InputEvent();
  AliAODEvent* aodevent = dynamic_cast<AliAODEvent*> (event);  
  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*> (event);
  
  AliAODHeader* header = AODEvent()->GetHeader();
  
  // Copy from AODs
  if(aodevent)
  {
    *header = *(aodevent->GetHeader());
    return;
  }
  
  // Copy from ESDs
  
  header->SetRunNumber(event->GetRunNumber());
  
  TTree* tree = fInputHandler->GetTree();
  if (tree) 
  {
    TFile* file = tree->GetCurrentFile();
    if (file) header->SetESDFileName(file->GetName());
  }
  
  header->SetBunchCrossNumber(event->GetBunchCrossNumber());
  header->SetOrbitNumber(event->GetOrbitNumber());
  header->SetPeriodNumber(event->GetPeriodNumber());
  header->SetEventType(event->GetEventType());
  
  //Centrality
  if(event->GetCentrality())
  {
    header->SetCentrality(new AliCentrality(*(event->GetCentrality())));
  }
  else{
    header->SetCentrality(0);
  }
  
  //Trigger  
  header->SetOfflineTrigger(fInputHandler->IsEventSelected()); // propagate the decision of the physics selection
  header->SetFiredTriggerClasses(esdevent->GetFiredTriggerClasses());
  header->SetTriggerMask(event->GetTriggerMask()); 
  header->SetTriggerCluster(event->GetTriggerCluster());
  header->SetL0TriggerInputs(esdevent->GetHeader()->GetL0TriggerInputs());    
  header->SetL1TriggerInputs(esdevent->GetHeader()->GetL1TriggerInputs());    
  header->SetL2TriggerInputs(esdevent->GetHeader()->GetL2TriggerInputs());    
  
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
  header->SetDiamondZ(esdevent->GetDiamondZ(),esdevent->GetSigma2DiamondZ());
  
}

//_____________________________________________
void AliAnalysisTaskCaloFilter::FillAODTracks()
{
  // AOD track copy
  
  if(!fFillTracks) return;
  
  AliVEvent* event = InputEvent();
  AliAODEvent* aodevent = dynamic_cast<AliAODEvent*> (event);  
  //AliESDEvent* esdevent = dynamic_cast<AliESDEvent*> (event);
  
  // Copy from AODs
  if(aodevent)
  {
    TClonesArray* inTracks = aodevent  ->GetTracks();
    TClonesArray* ouTracks = AODEvent()->GetTracks();
    
		new (ouTracks) TClonesArray(*inTracks);
    
    return;
  }
  
}  

//_______________________________________________
void AliAnalysisTaskCaloFilter::FillAODVertices()
{
  // Copy vertices
  
  AliVEvent* event = InputEvent();
  AliAODEvent* aodevent = dynamic_cast<AliAODEvent*> (event);  
  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*> (event);
  
  // set arrays and pointers
  Double_t pos[3]   ;
  Double_t covVtx[6];
  for (Int_t i = 0; i < 6; i++)  covVtx[i] = 0.;
  
  // Copy from AODs
  if(aodevent)
  {
    TClonesArray* inVertices = aodevent  ->GetVertices();
    TClonesArray* ouVertices = AODEvent()->GetVertices();
    
		new (ouVertices) TClonesArray(*inVertices);
    return;
  }
  
  // Copy from ESDs
  
  // Access to the AOD container of vertices
  Int_t jVertices=0;
  TClonesArray &vertices = *(AODEvent()->GetVertices());
  
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
  
}

//__________________________________________________
void AliAnalysisTaskCaloFilter::CorrectionsInEMCAL()
{
  //If EMCAL, and requested, correct energy, position ...
  
  //Need to do this in a separate loop before filling the ESDs because of the track matching recalculations
  
  if(fCorrect && (fCaloFilter==kEMCAL || fCaloFilter==kBoth) ) 
  {
    if(!fGeoMatrixSet)
    {
      if(fLoadEMCALMatrices)
      {
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
      else if(!gGeoManager)
      {
        printf("AliAnalysisTaskCaloFilter::UserExec() - Get geo matrices from data\n");
        //Still not implemented in AOD, just a workaround to be able to work at least with ESDs	
        if(!strcmp(InputEvent()->GetName(),"AliAODEvent")) 
        {
          if(DebugLevel() > 1) 
            printf("AliAnalysisTaskCaloFilter Use ideal geometry, values geometry matrix not kept in AODs.\n");
        }//AOD
        else 
        {	
          if(DebugLevel() > 1) printf("AliAnalysisTaskCaloFilter Load Misaligned matrices. \n");
          AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent()) ;
          if(!esd) 
          {
            printf("AliAnalysisTaskCaloFilter::UserExec() - This event does not contain ESDs?");
            return;
          }
          for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
          {
            //if(DebugLevel() > 1) 
            esd->GetEMCALMatrix(mod)->Print();
            if(esd->GetEMCALMatrix(mod)) fEMCALGeo->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
          } 
          fGeoMatrixSet=kTRUE;
        }//ESD
      }//Load matrices from Data 
      
    }//first event
    
    //Cluster Loop
    Int_t nCaloClus = InputEvent()->GetNumberOfCaloClusters();
    
    for (Int_t iClust=0; iClust<nCaloClus; ++iClust) 
    {
      AliVCluster * cluster = InputEvent()->GetCaloCluster(iClust);
      
      if(cluster->IsPHOS()) continue ;
      
      Float_t position[]={0,0,0};
      if(DebugLevel() > 2)
        printf("Check cluster %d for bad channels and close to border\n",cluster->GetID());
      if(fEMCALRecoUtils->ClusterContainsBadChannel(fEMCALGeo,cluster->GetCellsAbsId(), cluster->GetNCells())) continue;	
      
      if(DebugLevel() > 2)
      { 
        printf("Filter, before  : i %d, E %f, dispersion %f, m02 %f, m20 %f, distToBad %f\n",iClust,cluster->E(),
               cluster->GetDispersion(),cluster->GetM02(),cluster->GetM20(), cluster->GetDistanceToBadChannel());
        cluster->GetPosition(position);
        printf("Filter, before  : i %d, x %f, y %f, z %f\n",cluster->GetID(), position[0], position[1], position[2]);
      }
      
      //Recalculate distance to bad channels, if new list of bad channels provided
      fEMCALRecoUtils->RecalculateClusterDistanceToBadChannel(fEMCALGeo, InputEvent()->GetEMCALCells(), cluster);
      
      if(fEMCALRecoUtils->IsRecalibrationOn())	
      {
        fEMCALRecoUtils->RecalibrateClusterEnergy(fEMCALGeo, cluster, InputEvent()->GetEMCALCells());
        fEMCALRecoUtils->RecalculateClusterShowerShapeParameters(fEMCALGeo, InputEvent()->GetEMCALCells(),cluster);
        fEMCALRecoUtils->RecalculateClusterPID(cluster);
      }
      
      fEMCALRecoUtils->RecalculateClusterPosition(fEMCALGeo, InputEvent()->GetEMCALCells(),cluster);
      
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
    fEMCALRecoUtils->FindMatches(InputEvent(),0,fEMCALGeo);
    
  } // corrections in EMCAL
}

//_______________________________________________________
Bool_t AliAnalysisTaskCaloFilter::CheckForPrimaryVertex()
{
  //Check if the vertex was well reconstructed, copy from V0Reader of conversion group
  //It only works for ESDs
  
  AliESDEvent * event = dynamic_cast<AliESDEvent*> (InputEvent());
  
  // AODs
  if(!event) 
  {
    // Check that the vertex is not (0,0,0)
    Double_t v[3];
    InputEvent()->GetPrimaryVertex()->GetXYZ(v) ;
    
    if(TMath::Abs(v[2]) < 1e-6 && 
       TMath::Abs(v[1]) < 1e-6 && 
       TMath::Abs(v[0]) < 1e-6 ) return kFALSE ;
    
    return kTRUE;
  }
  
  // ESDs
  if(event->GetPrimaryVertexTracks()->GetNContributors() > 0)
  {
    return kTRUE;
  }
  
  if(event->GetPrimaryVertexTracks()->GetNContributors() < 1) 
  {
    // SPD vertex
    if(event->GetPrimaryVertexSPD()->GetNContributors() > 0) 
    {
      //cout<<"spd vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
      return kTRUE;
      
    }
    if(event->GetPrimaryVertexSPD()->GetNContributors() < 1) 
    {
      //      cout<<"bad vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
      return kFALSE;
    }
  }
  return kFALSE;
  
}

//_________________________________________
void AliAnalysisTaskCaloFilter::PrintInfo()
{
  //Print settings
  
  printf("TASK: AnalysisCaloFilter \n");
  printf("\t Not only filter, correct Clusters? %d\n",fCorrect);
  printf("\t Calorimeter Filtering Option     ? %d\n",fCaloFilter);
  //printf("\t Use handmade geo matrices?   EMCAL %d, PHOS %d\n",fLoadEMCALMatrices, fLoadPHOSMatrices);
  printf("\t Use handmade geo matrices?   EMCAL %d, PHOS 0\n",fLoadEMCALMatrices);
  printf("\t Fill AOD file? %d\n", fFillAODFile);
  printf("\t Fill Tracks? %d\n"  , fFillTracks);
  printf("Event Selection : EMCAL min E %f, EMCAL NCells %d, Vertex %f\n",fEnergyCut,fNcellsCut,fVzCut);
}

