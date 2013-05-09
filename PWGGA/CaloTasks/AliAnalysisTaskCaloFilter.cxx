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
// Keep tracks, v0s, VZERO if requested
// Select events containing a cluster or track avobe a given threshold
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
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

//EMCAL
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"

#include "AliAnalysisTaskCaloFilter.h"

ClassImp(AliAnalysisTaskCaloFilter)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskCaloFilter::AliAnalysisTaskCaloFilter():
AliAnalysisTaskSE("CaloFilterTask"), 
fCaloFilter(0),           fEventSelection(), 
fAcceptAllMBEvent(kFALSE),fMBTriggerMask(AliVEvent::kMB), 
fCorrect(kFALSE), 
fEMCALGeo(0x0),           fEMCALGeoName("EMCAL_COMPLETE12SMV1"), 
fEMCALRecoUtils(new AliEMCALRecoUtils),
fLoadEMCALMatrices(kFALSE), //fLoadPHOSMatrices(kFALSE),
fGeoMatrixSet(kFALSE),
fConfigName(""),          fFillAODFile(kTRUE), 
fFillMCParticles(kFALSE),
fFillTracks(kFALSE),      fFillHybridTracks(kFALSE),
fFillAllVertices(kFALSE), fFillv0s(kFALSE),  
fFillVZERO(kFALSE),
fEMCALEnergyCut(0.),      fEMCALNcellsCut (0),
fPHOSEnergyCut(0.),       fPHOSNcellsCut (0), 
fTrackPtCut(-1),
fVzCut(100.),             fCheckEventVertex(kTRUE),
fEvent(0x0),              
fESDEvent(0x0),           fAODEvent(0x0)
{
  // Default constructor
  
  fEventSelection[0] = kFALSE;
  fEventSelection[1] = kFALSE;
  fEventSelection[2] = kFALSE;
  
  for(Int_t i = 0; i < 12; i++) fEMCALMatrix[i] = 0 ;
  //for(Int_t i = 0; i < 5 ; i++) fPHOSMatrix[i]  = 0 ;
  
}

//__________________________________________________
AliAnalysisTaskCaloFilter::AliAnalysisTaskCaloFilter(const char* name):
AliAnalysisTaskSE(name), 
fCaloFilter(0),           fEventSelection(), 
fAcceptAllMBEvent(kFALSE),fMBTriggerMask(AliVEvent::kMB), 
fCorrect(kFALSE),
fEMCALGeo(0x0),           fEMCALGeoName("EMCAL_COMPLETE12SMV1"), 
fEMCALRecoUtils(new AliEMCALRecoUtils),
fLoadEMCALMatrices(kFALSE), //fLoadPHOSMatrices(kFALSE),
fGeoMatrixSet(kFALSE),
fConfigName(""),          fFillAODFile(kTRUE),
fFillMCParticles(kFALSE),
fFillTracks(kFALSE),      fFillHybridTracks(kFALSE),
fFillAllVertices(kFALSE), fFillv0s(kFALSE),
fFillVZERO(kFALSE),
fEMCALEnergyCut(0.),      fEMCALNcellsCut(0), 
fPHOSEnergyCut(0.),       fPHOSNcellsCut(0), 
fTrackPtCut(-1),
fVzCut(100.),             fCheckEventVertex(kTRUE),
fEvent(0x0),              
fESDEvent(0x0),           fAODEvent(0x0)
{
  // Constructor
  
  fEventSelection[0] = kFALSE;
  fEventSelection[1] = kFALSE;
  fEventSelection[2] = kFALSE;  
  
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

//_____________________________________________
Bool_t AliAnalysisTaskCaloFilter::AcceptEvent()
{
  // Define conditions to accept the event to be filtered
  
  if(!AcceptEventVertex()) return kFALSE;
  
  Bool_t eventSel = kFALSE;
  
  Bool_t isMB = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fMBTriggerMask);
  
  if     ( isMB && fAcceptAllMBEvent )                eventSel = kTRUE; // accept any MB event
  
  else if( fEventSelection[0] && AcceptEventEMCAL() ) eventSel = kTRUE; // accept event depending on EMCAL activity
  
  else if( fEventSelection[1] && AcceptEventPHOS () ) eventSel = kTRUE; // accept event depending on PHOS  activity
  
  else if( fEventSelection[2] && AcceptEventTrack() ) eventSel = kTRUE; // accept event depending on Track activity
    
  return eventSel ;
  
}

//__________________________________________________
Bool_t AliAnalysisTaskCaloFilter::AcceptEventEMCAL()
{
  // Accept event given there is a EMCAL cluster with enough energy, and not noisy, exotic
  
  if(fCaloFilter==kPHOS)   return kTRUE; // accept 
  
  if(fEMCALEnergyCut <= 0) return kTRUE; // accept
  
  Int_t           nCluster = InputEvent() -> GetNumberOfCaloClusters();
  AliVCaloCells * caloCell = InputEvent() -> GetEMCALCells();
  Int_t           bc       = InputEvent() -> GetBunchCrossNumber();
  
  for(Int_t icalo = 0; icalo < nCluster; icalo++)
  {
    AliVCluster *clus = (AliVCluster*) (InputEvent()->GetCaloCluster(icalo));
    
    if( ( clus->IsEMCAL() ) && ( clus->GetNCells() > fEMCALNcellsCut ) && ( clus->E() > fEMCALEnergyCut ) &&
       fEMCALRecoUtils->IsGoodCluster(clus,fEMCALGeo,caloCell,bc))
    { 
      
      if (fDebug > 0) printf("AliAnalysisTaskCaloFilter::AcceptEventEMCAL() - Accept :  E %2.2f > %2.2f, nCells %d > %d \n",
                             clus->E(), fEMCALEnergyCut, clus->GetNCells(), fEMCALNcellsCut);
      
      return kTRUE;
    }
    
  }// loop
  
  if (fDebug > 0)  printf("AliAnalysisTaskCaloFilter::AcceptEventEMCAL() - Reject \n");

  //printf("Fired %s\n",((AliESDEvent*)InputEvent())->GetFiredTriggerClasses().Data());

  return kFALSE;
  
}  

//_________________________________________________
Bool_t AliAnalysisTaskCaloFilter::AcceptEventPHOS()
{
  // Accept event given there is a PHOS cluster with enough energy and not noisy/exotic
  
  if(fCaloFilter==kEMCAL) return kTRUE; // accept 
  
  if(fPHOSEnergyCut <= 0) return kTRUE; // accept
  
  Int_t nCluster = InputEvent() -> GetNumberOfCaloClusters();
  
  for(Int_t icalo = 0; icalo < nCluster; icalo++)
  {
    AliVCluster *clus = (AliVCluster*) (InputEvent()->GetCaloCluster(icalo));
    
    if( ( clus->IsPHOS() ) && ( clus->GetNCells() > fPHOSNcellsCut ) && ( clus->E() > fPHOSEnergyCut ))
    { 
      
      if (fDebug > 0) printf("AliAnalysisTaskCaloFilter::AcceptEventPHOS() - Accept :  E %2.2f > %2.2f, nCells %d > %d \n", 
                             clus->E(), fPHOSEnergyCut, clus->GetNCells(), fPHOSNcellsCut);
      
      return kTRUE;
    }
    
  }// loop
  
  if (fDebug > 0)  printf("AliAnalysisTaskCaloFilter::AcceptEventPHOS() - Reject \n");
  
  return kFALSE;
  
}  

//__________________________________________________
Bool_t AliAnalysisTaskCaloFilter::AcceptEventTrack()
{
  // Accept event if there is a track avobe a certain pT
  
  if(fTrackPtCut <= 0) return kTRUE; // accept
  
  Double_t pTrack[3] = {0,0,0};
  
  for (Int_t nTrack = 0; nTrack < fEvent->GetNumberOfTracks(); ++nTrack) 
  {
    AliVTrack *track = (AliVTrack*) fEvent->GetTrack(nTrack);
    
    // Select only hybrid tracks?
    if(fAODEvent && fFillHybridTracks && !((AliAODTrack*)track)->IsHybridGlobalConstrainedGlobal()) continue;
    
    track->GetPxPyPz(pTrack) ;
    
    TLorentzVector momentum(pTrack[0],pTrack[1],pTrack[2],0);
    
    if(momentum.Pt() > fTrackPtCut) 
    {
      if (fDebug > 0) printf("AliAnalysisTaskCaloFilter::AcceptEventTrack() - Accept :  pT %2.2f > %2.2f \n", 
                             momentum.Pt(), fTrackPtCut);

      return kTRUE;
    }
    
  } 
  
  if (fDebug > 0)  printf("AliAnalysisTaskCaloFilter::AcceptEventTrack() - Reject \n");
  
  return kFALSE;
  
}  

//___________________________________________________
Bool_t AliAnalysisTaskCaloFilter::AcceptEventVertex()
{
  // Accept event with good vertex
  
  Double_t v[3];
  InputEvent()->GetPrimaryVertex()->GetXYZ(v) ;
  
  if(TMath::Abs(v[2]) > fVzCut) 
  {
    if (fDebug > 0) printf("AliAnalysisTaskCaloFilter::AcceptEventVertex() - Vz Reject : vz %2.2f > %2.2f\n",v[2],fVzCut);
    
    return kFALSE ;
  }
  
  return CheckForPrimaryVertex();
}  

//_______________________________________________________
Bool_t AliAnalysisTaskCaloFilter::CheckForPrimaryVertex()
{
  //Check if the vertex was well reconstructed, copy from v0Reader of conversion group
  //It only works for ESDs
  
  if(!fCheckEventVertex) return kTRUE;

  // AODs
  if(!fESDEvent) 
  {
    // Check that the vertex is not (0,0,0)
    Double_t v[3];
    InputEvent()->GetPrimaryVertex()->GetXYZ(v) ;
    
    if(TMath::Abs(v[2]) < 1e-6 && 
       TMath::Abs(v[1]) < 1e-6 && 
       TMath::Abs(v[0]) < 1e-6 ) 
    {
      if (fDebug > 0)  printf("AliAnalysisTaskCaloFilter::CheckForPrimaryVertex() - Reject v(0,0,0) \n");
      
      return kFALSE ;
    }
    
    return kTRUE;
  }
  
  // ESDs
  if(fESDEvent->GetPrimaryVertexTracks()->GetNContributors() > 0)
  {
    return kTRUE;
  }
  
  if(fESDEvent->GetPrimaryVertexTracks()->GetNContributors() < 1) 
  {
    // SPD vertex
    if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors() > 0) 
    {
      //cout<<"spd vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
      return kTRUE;
      
    }
    if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors() < 1) 
    {
      //      cout<<"bad vertex type::"<< fESDEvent->GetPrimaryVertex()->GetName() << endl;
      if (fDebug > 0)  printf("AliAnalysisTaskCaloFilter::CheckForPrimaryVertex() - Reject, GetPrimaryVertexSPD()->GetNContributors() < 1 \n");
      
      return kFALSE;
    }
  }
  
  if (fDebug > 0)  printf("AliAnalysisTaskCaloFilter::CheckForPrimaryVertex() - Reject, GetPrimaryVertexTracks()->GetNContributors() > 1 \n");
  
  return kFALSE;
  
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

//________________________________________________
void AliAnalysisTaskCaloFilter::FillAODCaloCells()
{  
  // Fill EMCAL/PHOS cell info
  
  // EMCAL
  if ((fCaloFilter==kBoth ||  fCaloFilter==kEMCAL) && fEvent->GetEMCALCells()) 
  { // protection against missing ESD information
    AliVCaloCells &eventEMcells = *(fEvent->GetEMCALCells());
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
  if ((fCaloFilter==kBoth ||  fCaloFilter==kPHOS) && fEvent->GetPHOSCells()) 
  { // protection against missing ESD information
    AliVCaloCells &eventPHcells = *(fEvent->GetPHOSCells());
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
  
  // Access to the AOD container of clusters
  
  TClonesArray &caloClusters = *(AODEvent()->GetCaloClusters());
  Int_t jClusters=0;
  Float_t  posF[3]  ;
  
  Int_t nCaloClus = fEvent->GetNumberOfCaloClusters();
  for (Int_t iClust=0; iClust<nCaloClus; ++iClust) 
  {
    
    AliVCluster * cluster = fEvent->GetCaloCluster(iClust);
    
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

//__________________________________________________
void AliAnalysisTaskCaloFilter::FillAODCaloTrigger()
{
  // AOD CaloTrigger copy
  
  if( !AODEvent() || !fAODEvent ) return;
  
  AliAODCaloTrigger* triggerEM   = AODEvent()->GetCaloTrigger("EMCAL");
  AliAODCaloTrigger* triggerPH   = AODEvent()->GetCaloTrigger("PHOS");
  
  // Copy from AODs
  
  AliAODCaloTrigger* inTriggerEM = fAODEvent ->GetCaloTrigger("EMCAL");
  AliAODCaloTrigger* inTriggerPH = fAODEvent ->GetCaloTrigger("PHOS");
  
  if(inTriggerPH && (fCaloFilter==kBoth || fCaloFilter==kPHOS))  *triggerPH = *inTriggerPH;
  
  if(inTriggerEM && (fCaloFilter==kBoth || fCaloFilter==kEMCAL)) *triggerEM = *inTriggerEM;  
}  

//______________________________________________
void AliAnalysisTaskCaloFilter::FillAODHeader()
{
  // AOD header copy
  
  AliAODHeader* header = AODEvent()->GetHeader();
  
  // Copy from AODs
  if(fAODEvent)
  {
    *header = *(fAODEvent->GetHeader());
    return;
  }
  
  if(!fESDEvent) return;
  
  // Copy from ESDs
  
  header->SetRunNumber(fEvent->GetRunNumber());
  
  TTree* tree = fInputHandler->GetTree();
  if (tree) 
  {
    TFile* file = tree->GetCurrentFile();
    if (file) header->SetESDFileName(file->GetName());
  }
  
  header->SetBunchCrossNumber(fEvent->GetBunchCrossNumber());
  header->SetOrbitNumber(fEvent->GetOrbitNumber());
  header->SetPeriodNumber(fEvent->GetPeriodNumber());
  header->SetEventType(fEvent->GetEventType());
  
  //Centrality
  if(fEvent->GetCentrality())
  {
    header->SetCentrality(new AliCentrality(*(fEvent->GetCentrality())));
  }
  else
  {
    header->SetCentrality(0);
  }
  
  //Trigger  
  header->SetOfflineTrigger(fInputHandler->IsEventSelected()); // propagate the decision of the physics selection
  header->SetFiredTriggerClasses(fESDEvent->GetFiredTriggerClasses());
  header->SetTriggerMask(fEvent->GetTriggerMask()); 
  header->SetTriggerCluster(fEvent->GetTriggerCluster());
  header->SetL0TriggerInputs(fESDEvent->GetHeader()->GetL0TriggerInputs());    
  header->SetL1TriggerInputs(fESDEvent->GetHeader()->GetL1TriggerInputs());    
  header->SetL2TriggerInputs(fESDEvent->GetHeader()->GetL2TriggerInputs());    
  
  header->SetMagneticField(fEvent->GetMagneticField());
  header->SetMuonMagFieldScale(fESDEvent->GetCurrentDip()/6000.); 
  
  header->SetZDCN1Energy(fEvent->GetZDCN1Energy());
  header->SetZDCP1Energy(fEvent->GetZDCP1Energy());
  header->SetZDCN2Energy(fEvent->GetZDCN2Energy());
  header->SetZDCP2Energy(fEvent->GetZDCP2Energy());
  header->SetZDCEMEnergy(fEvent->GetZDCEMEnergy(0),fEvent->GetZDCEMEnergy(1));
  
  Float_t diamxy[2]={fEvent->GetDiamondX(),fEvent->GetDiamondY()};
  Float_t diamcov[3];
  fEvent->GetDiamondCovXY(diamcov);
  header->SetDiamond(diamxy,diamcov);
  header->SetDiamondZ(fESDEvent->GetDiamondZ(),fESDEvent->GetSigma2DiamondZ());
  
}


//__________________________________________________
void AliAnalysisTaskCaloFilter::FillAODMCParticles()
{
  // Copy MC particles
  
  if(!fFillMCParticles) return;
  
  TClonesArray* inMCParticles = (TClonesArray*) (fAODEvent  ->FindListObject("mcparticles"));
  TClonesArray* ouMCParticles = (TClonesArray*) ( AODEvent()->FindListObject("mcparticles"));
  
  if( inMCParticles &&  ouMCParticles ) new (ouMCParticles) TClonesArray(*inMCParticles);
    
}  

//_____________________________________________
void AliAnalysisTaskCaloFilter::FillAODTracks()
{
  // AOD track copy
  
  if(!fFillTracks) return;
  
  AliAODTrack* aodTrack(0x0);
  
  Double_t pos[3]   = { 0. };      
  Double_t covTr[21]= { 0. };
  Double_t pid[10]  = { 0. };  
  Double_t p[3]     = { 0. };
    
  // Copy from AODs
  if(fAODEvent)
  {
    //TClonesArray* inTracks = fAODEvent  ->GetTracks();
    TClonesArray* ouTracks = AODEvent()->GetTracks();
		//new (ouTracks) TClonesArray(*inTracks);
    
    //printf("N tracks %d\n",fAODEvent->GetNumberOfTracks());
    Int_t nCopyTrack = 0;
    for (Int_t nTrack = 0; nTrack < fAODEvent->GetNumberOfTracks(); ++nTrack) 
    {
      AliAODTrack *track = fAODEvent->GetTrack(nTrack);
      
      // Select only hybrid tracks?
      if(fFillHybridTracks && !track->IsHybridGlobalConstrainedGlobal()) continue;
      
      // Remove PID object to save space
      //track->SetDetPID(0x0);
      
      //new((*ouTracks)[nCopyTrack++]) AliAODTrack(*track);
      
      track->GetPxPyPz(p);
      Bool_t isDCA  = track->GetPosition(pos);
      track->GetCovMatrix(covTr);
      track->GetPID(pid);
      
      AliAODVertex* primVertex = (AliAODVertex*) AODEvent()->GetVertices()->At(0); // primary vertex, copied previously!!!

      aodTrack = new((*ouTracks)[nCopyTrack++]) AliAODTrack(
                                                            track->GetID(),
                                                            track->GetLabel(),
                                                            p,
                                                            kTRUE,
                                                            pos,
                                                            isDCA,
                                                            covTr, 
                                                            track->Charge(),
                                                            track->GetITSClusterMap(), 
                                                            pid,
                                                            primVertex,
                                                            track->GetUsedForVtxFit(),
                                                            track->GetUsedForPrimVtxFit(),
                                                            (AliAODTrack::AODTrk_t) track->GetType(), 
                                                            track->GetFilterMap(),
                                                            track->Chi2perNDF());
      
      
      aodTrack->SetIsHybridGlobalConstrainedGlobal(track->IsHybridGlobalConstrainedGlobal());   
      aodTrack->SetIsHybridTPCConstrainedGlobal   (track->IsHybridTPCConstrainedGlobal());    
      aodTrack->SetIsGlobalConstrained            (track->IsGlobalConstrained());  
      aodTrack->SetIsTPCConstrained               (track->IsTPCConstrained());    
      
      aodTrack->SetTPCFitMap    (track->GetTPCFitMap());
      aodTrack->SetTPCClusterMap(track->GetTPCClusterMap());
      aodTrack->SetTPCSharedMap (track->GetTPCSharedMap());
      
      aodTrack->SetChi2MatchTrigger(track->GetChi2MatchTrigger());
      
      // set the DCA values to the AOD track
    
      aodTrack->SetPxPyPzAtDCA(track->PxAtDCA(),track->PyAtDCA(),track->PzAtDCA());
      aodTrack->SetXYAtDCA    (track->XAtDCA() ,track->YAtDCA());
      
      aodTrack->SetFlags      (track->GetFlags());
      aodTrack->SetTPCPointsF (track->GetTPCNclsF());
      
      // Calo 
      
      if(track->IsEMCAL()) aodTrack->SetEMCALcluster(track->GetEMCALcluster());
      if(track->IsPHOS())  aodTrack->SetPHOScluster (track->GetPHOScluster());
      aodTrack->SetTrackPhiEtaOnEMCal( track->GetTrackPhiOnEMCal(),  track->GetTrackPhiOnEMCal() );
          
    } 
    
    //printf("Final N tracks %d\n",nCopyTrack);
    
    return;
  }
  
}  

//_________________________________________
void AliAnalysisTaskCaloFilter::FillAODv0s()
{
  // Copy v0s (use if you know what you do, use quite a lot of memory)
  
  if(!fFillv0s) return;
  
  // Copy from AODs
  if(fAODEvent)
  {
    TClonesArray* inv0 = fAODEvent  ->GetV0s();
    TClonesArray* ouv0 =  AODEvent()->GetV0s();
    
		//new (ouv0s) TClonesArray(*inv0s);
    
    Int_t allv0s = inv0->GetEntriesFast();
    
    for (Int_t nv0s = 0; nv0s < allv0s; ++nv0s) 
    {
      AliAODv0 *v0 = (AliAODv0*)inv0->At(nv0s);
      
      new((*ouv0)[nv0s]) AliAODv0(*v0);
    } 
    
    return;
  }
  
}  

//____________________________________________
void AliAnalysisTaskCaloFilter::FillAODVZERO()
{
  // Copy VZERO
  
  if(!fFillVZERO) return;
  
  AliAODVZERO* vzeroData = AODEvent()->GetVZEROData();
  
  if(fESDEvent) *vzeroData = *(fESDEvent->GetVZEROData());
  else          *vzeroData = *(fAODEvent->GetVZEROData());
  
}  

//_______________________________________________
void AliAnalysisTaskCaloFilter::FillAODVertices()
{
  // Copy vertices
  
  // set arrays and pointers
  Double_t pos[3]   ;
  Double_t covVtx[6];
  for (Int_t i = 0; i < 6; i++)  covVtx[i] = 0.;
  
  // Copy from AODs
  if(fAODEvent)
  {
    TClonesArray* inVertices = fAODEvent  ->GetVertices();
    TClonesArray* ouVertices = AODEvent()->GetVertices();
    
		//new (ouVertices) TClonesArray(*inVertices);
    
    //Keep only the first 3 vertices if not requested
    Int_t allVertices = inVertices->GetEntriesFast();
    
    //printf("n Vertices %d\n",allVertices);
    
    if(!fFillAllVertices) 
    {
      if(allVertices >  3) allVertices = 3;
    }
    
    //printf("Final n Vertices %d\n",allVertices);
    
    for (Int_t nVertices = 0; nVertices < allVertices; ++nVertices) 
    {
      AliAODVertex *vertex = (AliAODVertex*)inVertices->At(nVertices);
      
      new((*ouVertices)[nVertices]) AliAODVertex(*vertex);
    } 
    
    return;
  }
  
  if(!fESDEvent) return;
  
  // Copy from ESDs
  
  // Access to the AOD container of vertices
  Int_t jVertices=0;
  TClonesArray &vertices = *(AODEvent()->GetVertices());
  
  // Add primary vertex. The primary tracks will be defined
  // after the loops on the composite objects (v0, cascades, kinks)
  fEvent   ->GetPrimaryVertex()->GetXYZ(pos);
  fESDEvent->GetPrimaryVertex()->GetCovMatrix(covVtx);
  Float_t chi = fESDEvent->GetPrimaryVertex()->GetChi2toNDF();
  
  AliAODVertex * primary = new(vertices[jVertices++])
  AliAODVertex(pos, covVtx, chi, NULL, -1, AliAODVertex::kPrimary);
  primary->SetName(fEvent->GetPrimaryVertex()->GetName());
  primary->SetTitle(fEvent->GetPrimaryVertex()->GetTitle());
  
}

//____________________________________
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
    fFillHybridTracks  = filter->fFillHybridTracks;
    fFillv0s           = filter->fFillv0s;
    fFillVZERO         = filter->fFillVZERO;
    fFillAllVertices   = filter->fFillAllVertices;
    fEMCALRecoUtils    = filter->fEMCALRecoUtils; 
    fConfigName        = filter->fConfigName;
    fCaloFilter        = filter->fCaloFilter;
    fEventSelection[0] = filter->fEventSelection[0];
    fEventSelection[1] = filter->fEventSelection[1];
    fEventSelection[2] = filter->fEventSelection[2];
    fAcceptAllMBEvent  = filter->fAcceptAllMBEvent;
    fMBTriggerMask     = filter->fMBTriggerMask;
    fCorrect           = filter->fCorrect;
    fEMCALEnergyCut    = filter->fEMCALEnergyCut;
    fEMCALNcellsCut    = filter->fEMCALNcellsCut;
    fPHOSEnergyCut     = filter->fPHOSEnergyCut;
    fPHOSNcellsCut     = filter->fPHOSNcellsCut;
    fTrackPtCut        = filter->fTrackPtCut;
    fVzCut             = filter->fVzCut;
    fCheckEventVertex  = filter->fCheckEventVertex;

    for(Int_t i = 0; i < 12; i++) fEMCALMatrix[i] = filter->fEMCALMatrix[i] ;
  }
} 

//_________________________________________
void AliAnalysisTaskCaloFilter::PrintInfo()
{
  //Print settings
  
  printf("AnalysisCaloFilter::PrintInfo() \n");
  
  printf("\t Not only filter, correct Clusters? %d\n",fCorrect);
  printf("\t Calorimeter Filtering Option     ? %d\n",fCaloFilter);
  
  //printf("\t Use handmade geo matrices?   EMCAL %d, PHOS %d\n",fLoadEMCALMatrices, fLoadPHOSMatrices);
  printf("\t Use handmade geo matrices?   EMCAL %d, PHOS 0\n",fLoadEMCALMatrices);
  
  printf("\t Fill: AOD file? %d Tracks? %d; all Vertex? %d; v0s? %d; VZERO ? %d\n", 
         fFillAODFile,fFillTracks,fFillAllVertices, fFillv0s, fFillVZERO);
  
  printf("\t Event Selection based : EMCAL?  %d, PHOS?  %d Tracks? %d - Accept all MB with mask %d? %d\n",
         fEventSelection[0],fEventSelection[1],fEventSelection[2],fMBTriggerMask, fAcceptAllMBEvent);
  
  printf("\t \t EMCAL E > %2.2f, EMCAL nCells >= %d, PHOS E > %2.2f, PHOS nCells >= %d, Track pT > %2.2f, |vz| < %2.2f\n",
         fEMCALEnergyCut,fEMCALNcellsCut,fPHOSEnergyCut,fPHOSNcellsCut, fTrackPtCut,fVzCut);
}

//_______________________________________________________
void AliAnalysisTaskCaloFilter::UserCreateOutputObjects()
{
  // Init geometry 
	
  fEMCALGeo =  AliEMCALGeometry::GetInstance(fEMCALGeoName) ;	
  
  if(fFillMCParticles)
  {
    TClonesArray * aodMCParticles = new TClonesArray("AliAODMCParticle",500);
		aodMCParticles->SetName("mcparticles");
		((AliAODHandler*)AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler())->AddBranch("TClonesArray", &aodMCParticles);
  }
  
}  

//____________________________________________________________
void AliAnalysisTaskCaloFilter::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event
  // Copy input ESD or AOD header, vertex, CaloClusters and CaloCells to output AOD
  
  if (fDebug > 0)  
    printf("CaloFilter: Analysing event # %d\n", (Int_t)Entry());
  
  fEvent    = InputEvent();
  fAODEvent = dynamic_cast<AliAODEvent*> (fEvent);  
  fESDEvent = dynamic_cast<AliESDEvent*> (fEvent);
  
  if(!fEvent) 
  {
    printf("AliAnalysisTaskCaloFilter::UserExec - This event does not contain Input?");
    return;
  }
  
  // printf("Start processing : %s\n",fAODEvent->GetFiredTriggerClasses().Data());
  
  // Select the event
  
  if(!AcceptEvent()) return ;
  
  //Magic line to write events to file
  
  AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(fFillAODFile);
  
  // Reset output AOD
  
  Int_t nVertices = 0;
  if(fFillv0s) nVertices = fEvent->GetNumberOfV0s();
  Int_t nCaloClus = fEvent->GetNumberOfCaloClusters();
  Int_t nTracks   = fEvent->GetNumberOfTracks();
  
  AODEvent()->ResetStd(nTracks, nVertices, 0, 0, 0, nCaloClus, 0, 0);
  
  // Copy
  
  FillAODHeader();
  
  // 
  FillAODv0s();
  
  //
  FillAODVertices(); // Do it before the track filtering to have the reference to the vertex
  
  // 
  FillAODVZERO();
  
  //
  FillAODTracks();
  
  //
  CorrectionsInEMCAL();
  
  //
  FillAODCaloClusters();
  
  //
  FillAODCaloCells();
  
  //
  FillAODCaloTrigger();
  
  // 
  FillAODMCParticles();
  
  //printf("Filtered event, end processing : %s\n",fAODEvent->GetFiredTriggerClasses().Data());
  
}

