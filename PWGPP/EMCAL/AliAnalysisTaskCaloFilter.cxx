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

// Root
#include "TGeoManager.h"
#include "TFile.h"
#include "TROOT.h"
#include "TInterpreter.h"

// STEER
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliLog.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliVEventHandler.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

// EMCAL
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"

#include "AliAnalysisTaskCaloFilter.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskCaloFilter) ;
/// \endcond

//_____________________________________________________
/// Default constructor.
//_____________________________________________________
AliAnalysisTaskCaloFilter::AliAnalysisTaskCaloFilter():
AliAnalysisTaskSE("CaloFilterTask"), 
fCaloFilter(0),           fEventSelection(), 
fAcceptAllMBEvent(kFALSE),fMBTriggerMask(AliVEvent::kMB), 
fCorrect(kFALSE), 
fEMCALGeo(0x0),           fEMCALGeoName("EMCAL_COMPLETE12SMV1_DCAL_8SM"),
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
  fEventSelection[0] = kFALSE;
  fEventSelection[1] = kFALSE;
  fEventSelection[2] = kFALSE;
  
  for(Int_t i = 0; i < 22; i++) fEMCALMatrix[i] = 0 ;
  //for(Int_t i = 0; i < 5 ; i++) fPHOSMatrix[i]  = 0 ;
}

//_____________________________________________________________________
/// Constructor.
//_____________________________________________________________________
AliAnalysisTaskCaloFilter::AliAnalysisTaskCaloFilter(const char* name):
AliAnalysisTaskSE(name), 
fCaloFilter(0),           fEventSelection(), 
fAcceptAllMBEvent(kFALSE),fMBTriggerMask(AliVEvent::kMB), 
fCorrect(kFALSE),
fEMCALGeo(0x0),           fEMCALGeoName("EMCAL_COMPLETE12SMV1_DCAL_8SM"),
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
  fEventSelection[0] = kFALSE;
  fEventSelection[1] = kFALSE;
  fEventSelection[2] = kFALSE;  
  
  for(Int_t i = 0; i < 22; i++) fEMCALMatrix[i] = 0 ;
  //for(Int_t i = 0; i < 5 ; i++) fPHOSMatrix[i]  = 0 ;
}

//_____________________________________________________
/// Destructor.
//_____________________________________________________
AliAnalysisTaskCaloFilter::~AliAnalysisTaskCaloFilter()
{
  if(fEMCALGeo)        delete fEMCALGeo;	
 
  if(fEMCALRecoUtils)  delete fEMCALRecoUtils;
}

//_____________________________________________
/// \return True if the event is accepted
///  Criteria to accept the event:
///  * Vertex selection
///  * If all Min Bias events are selected
///  * Events with activity in EMCal
///  * Events with activity in PHOS
///  * Events with tracks
//_____________________________________________
Bool_t AliAnalysisTaskCaloFilter::AcceptEvent()
{
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
/// \return True if there is signal in EMCal
/// Accept event given there is a EMCAL cluster with
/// enough energy, cells and not noisy, exotic.
//__________________________________________________
Bool_t AliAnalysisTaskCaloFilter::AcceptEventEMCAL()
{
  if( fCaloFilter == kPHOS ) return kTRUE; // accept
  
  if( fEMCALEnergyCut <= 0 ) return kTRUE; // accept
  
  Int_t           nCluster = InputEvent() -> GetNumberOfCaloClusters();
  AliVCaloCells * caloCell = InputEvent() -> GetEMCALCells();
  Int_t           bc       = InputEvent() -> GetBunchCrossNumber();
  
  for(Int_t icalo = 0; icalo < nCluster; icalo++)
  {
    AliVCluster *clus = (AliVCluster*) (InputEvent()->GetCaloCluster(icalo));
    
    if( ( clus->IsEMCAL() ) && ( clus->GetNCells() > fEMCALNcellsCut ) && ( clus->E() > fEMCALEnergyCut ) &&
       fEMCALRecoUtils->IsGoodCluster(clus,fEMCALGeo,caloCell,bc))
    {
      
      AliDebug(1,Form("Accept :  E %2.2f > %2.2f, nCells %d > %d",
                      clus->E(), fEMCALEnergyCut, clus->GetNCells(), fEMCALNcellsCut));
      
      return kTRUE;
    }
  }// loop
  
  AliDebug(1,"Reject");

  //printf("Fired %s\n",((AliESDEvent*)InputEvent())->GetFiredTriggerClasses().Data());

  return kFALSE;
}  

//__________________________________________________
/// \return True if there is signal in PHOS
/// Accept event given there is a PHOS cluster with
/// enough energy and cells.
//__________________________________________________
Bool_t AliAnalysisTaskCaloFilter::AcceptEventPHOS()
{
  if( fCaloFilter == kEMCAL ) return kTRUE; // accept
  
  if( fPHOSEnergyCut <= 0   ) return kTRUE; // accept
  
  Int_t nCluster = InputEvent() -> GetNumberOfCaloClusters();
  
  for(Int_t icalo = 0; icalo < nCluster; icalo++)
  {
    AliVCluster *clus = (AliVCluster*) (InputEvent()->GetCaloCluster(icalo));
    
    if( ( clus->IsPHOS() ) && ( clus->GetNCells() > fPHOSNcellsCut ) && ( clus->E() > fPHOSEnergyCut ))
    { 
      
      AliDebug(1,Form("Accept :  E %2.2f > %2.2f, nCells %d > %d",
                      clus->E(), fPHOSEnergyCut, clus->GetNCells(), fPHOSNcellsCut));
      
      return kTRUE;
    }
  }// loop
  
  AliDebug(1,"Reject");
  
  return kFALSE;
}  

//__________________________________________________
/// \return True if there is signal in TPC
/// Accept event given there is a track with
/// enough pT. If requested, select only hybrid tracks.
//__________________________________________________
Bool_t AliAnalysisTaskCaloFilter::AcceptEventTrack()
{
  // Accept event if there is a track avobe a certain pT
  
  if( fTrackPtCut <= 0 ) return kTRUE; // accept
  
  //Double_t pTrack[3] = {0,0,0};
  
  for (Int_t nTrack = 0; nTrack < fEvent->GetNumberOfTracks(); ++nTrack) 
  {
    AliVTrack *track = (AliVTrack*) fEvent->GetTrack(nTrack);
    
    // Select only hybrid tracks?
    if(fAODEvent && fFillHybridTracks && !((AliAODTrack*)track)->IsHybridGlobalConstrainedGlobal()) continue ;
    
    //track->GetPxPyPz(pTrack) ;
    
    //TLorentzVector momentum(pTrack[0],pTrack[1],pTrack[2],0);
    
    //if( momentum.Pt() > fTrackPtCut )
    if( track->Pt() > fTrackPtCut )
    {
      AliDebug(1,Form("Accept :  pT %2.2f > %2.2f",track->Pt(), fTrackPtCut));

      return kTRUE;
    }
  } // loop
  
  AliDebug(1,"Reject");
  
  return kFALSE;
}  

//___________________________________________________
/// \return True if event with good vertex
//___________________________________________________
Bool_t AliAnalysisTaskCaloFilter::AcceptEventVertex()
{
  Double_t v[3];
  InputEvent()->GetPrimaryVertex()->GetXYZ(v) ;
  
  if(TMath::Abs(v[2]) > fVzCut) 
  {
    AliDebug(1,Form("Vz Reject : vz %2.2f > %2.2f",v[2],fVzCut));
    
    return kFALSE ;
  }
  
  return CheckForPrimaryVertex();
}  

//_______________________________________________________
/// Check if the vertex was well reconstructed, copy from
/// v0Reader of conversion group. Call corresponding selection
/// for ESDs and AODs.
//_______________________________________________________
Bool_t AliAnalysisTaskCaloFilter::CheckForPrimaryVertex()
{
  if(!fCheckEventVertex) return kTRUE;

//  // Check that the vertex is not (0,0,0)
//  Double_t v[3];
//  InputEvent()->GetPrimaryVertex()->GetXYZ(v) ;
//    
//  if(TMath::Abs(v[2]) < 1e-6 &&
//     TMath::Abs(v[1]) < 1e-6 &&
//     TMath::Abs(v[0]) < 1e-6 )
//  {
//    AliDebug(1,"Reject v(0,0,0)");
//        
//    return kFALSE ;
//  }
    
  if ( fESDEvent ) return CheckForPrimaryVertexInESDs();
  else             return CheckForPrimaryVertexInAODs();
}

//_____________________________________________________________
/// Check if the ESDs vertex was well reconstructed,
/// copy from v0Reader of conversion group.
//_____________________________________________________________
Bool_t AliAnalysisTaskCaloFilter::CheckForPrimaryVertexInESDs()
{
  if(fESDEvent->GetPrimaryVertexTracks()->GetNContributors() > 0)
    return kTRUE;

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
      AliDebug(1,"Reject, GetPrimaryVertexSPD()->GetNContributors() < 1");

      return kFALSE;
    }
  }
    
  AliDebug(1,"Reject, GetPrimaryVertexTracks()->GetNContributors() > 1");
    
  return kFALSE;
}

//_____________________________________________________________
/// Check if the AODs vertex was well reconstructed,
/// copy from v0Reader of conversion group.
//_____________________________________________________________
Bool_t AliAnalysisTaskCaloFilter::CheckForPrimaryVertexInAODs()
{
  if (fAODEvent->GetPrimaryVertex() != NULL)
  {
    if(fAODEvent->GetPrimaryVertex()->GetNContributors() > 0) return kTRUE;
  }
    
  if(fAODEvent->GetPrimaryVertexSPD() != NULL)
  {
    if(fAODEvent->GetPrimaryVertexSPD()->GetNContributors() > 0)
    {
      return kTRUE;
    }
    else
    {
      AliWarning(Form("Number of contributors from bad vertex type:: %s",
                      fAODEvent->GetPrimaryVertex()->GetName()));
      return kFALSE;
    }
  }
    
  AliDebug(1,"Reject, GetPrimaryVertexTracks()->GetNContributors() > 1");
    
  return kFALSE;
}

//__________________________________________________
/// If filtering EMCal, correct energy, position ...
/// Need to do this in a separate loop before filling
/// the output CaloClusters because of the track-matching recalculations
//__________________________________________________
void AliAnalysisTaskCaloFilter::CorrectionsInEMCAL()
{
  if(fCorrect && (fCaloFilter==kEMCAL || fCaloFilter==kBoth) )
  {
    if(!fGeoMatrixSet)
    {
      if(fLoadEMCALMatrices)
      {
        AliInfo("Load user defined EMCAL geometry matrices");
        for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
        {
          if(fEMCALMatrix[mod])
          {
            if(DebugLevel() > 1) 
              fEMCALMatrix[mod]->Print();
            fEMCALGeo->SetMisalMatrix(fEMCALMatrix[mod],mod) ;  
          }
          fGeoMatrixSet=kTRUE;
        }//SM loop
      }//Load matrices
      else if(!gGeoManager)
      {
        AliInfo("Get geo matrices from data");
        //Still not implemented in AOD, just a workaround to be able to work at least with ESDs	
        if(!strcmp(InputEvent()->GetName(),"AliAODEvent")) 
        {
          AliDebug(1,"Use ideal geometry, values geometry matrix not kept in AODs");
        }//AOD
        else 
        {	
          AliDebug(1,"Load Misaligned matrices");
          AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent()) ;
          if(!esd) 
          {
            AliInfo("This event does not contain ESDs?");
            return;
          }
          for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
          {
            if(DebugLevel() > 1)
              esd->GetEMCALMatrix(mod)->Print();
            if(esd->GetEMCALMatrix(mod)) fEMCALGeo->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
          } 
          fGeoMatrixSet=kTRUE;
        }//ESD
      }//Load matrices from Data 
      
    }//first event
    
    // Cluster Loop
    Int_t nCaloClus = InputEvent()->GetNumberOfCaloClusters();
    
    for (Int_t iClust=0; iClust<nCaloClus; ++iClust) 
    {
      AliVCluster * cluster = InputEvent()->GetCaloCluster(iClust);
      
      if(cluster->IsPHOS()) continue ;
      
      Float_t position[]={0,0,0};
      
      AliDebug(1,Form("Check cluster %d for bad channels and close to border",cluster->GetID()));
      
      if(fEMCALRecoUtils->ClusterContainsBadChannel(fEMCALGeo,cluster->GetCellsAbsId(), cluster->GetNCells())) continue;	
      
      AliDebug(2,Form("Filter, before  : i %d, E %f, dispersion %f, m02 %f, m20 %f, distToBad %f",iClust,cluster->E(),
                      cluster->GetDispersion(),cluster->GetM02(),cluster->GetM20(), cluster->GetDistanceToBadChannel()));
      cluster->GetPosition(position);
      AliDebug(2,Form("Filter, before  : i %d, x %f, y %f, z %f",cluster->GetID(), position[0], position[1], position[2]));
      
      //Recalculate distance to bad channels, if new list of bad channels provided
      fEMCALRecoUtils->RecalculateClusterDistanceToBadChannel(fEMCALGeo, InputEvent()->GetEMCALCells(), cluster);
      
      if(fEMCALRecoUtils->IsRecalibrationOn())	
      {
        fEMCALRecoUtils->RecalibrateClusterEnergy(fEMCALGeo, cluster, InputEvent()->GetEMCALCells());
        fEMCALRecoUtils->RecalculateClusterShowerShapeParameters(fEMCALGeo, InputEvent()->GetEMCALCells(),cluster);
        fEMCALRecoUtils->RecalculateClusterPID(cluster);
      }
      
      fEMCALRecoUtils->RecalculateClusterPosition(fEMCALGeo, InputEvent()->GetEMCALCells(),cluster);
      
      AliDebug(2,Form("Filter, after   : i %d, E %f, dispersion %f, m02 %f, m20 %f, distToBad %f",cluster->GetID(),cluster->E(),
                      cluster->GetDispersion(),cluster->GetM02(),cluster->GetM20(), cluster->GetDistanceToBadChannel()));
      cluster->GetPosition(position);
      AliDebug(1,Form("Filter, after   : i %d, x %f, y %f, z %f",cluster->GetID(), position[0], position[1], position[2]));
      
      cluster->SetE(fEMCALRecoUtils->CorrectClusterEnergyLinearity(cluster));
      
    }
    
    //Recalculate track-matching
    fEMCALRecoUtils->FindMatches(InputEvent(),0,fEMCALGeo);
    
  } // corrections in EMCAL
}

//________________________________________________
/// Fill EMCAL/PHOS cell info
//________________________________________________
void AliAnalysisTaskCaloFilter::FillAODCaloCells()
{
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
/// Fill the output AOD with input CaloClusters (ESD or AOD)
/// Access to the output AOD container of clusters.
//___________________________________________________
void AliAnalysisTaskCaloFilter::FillAODCaloClusters()
{
  TClonesArray &caloClusters = *(AODEvent()->GetCaloClusters());
  Int_t jClusters=0;
  Float_t  posF[3]  ;
  
  Int_t nCaloClus = fEvent->GetNumberOfCaloClusters();
  for (Int_t iClust=0; iClust<nCaloClus; ++iClust) 
  {
    AliVCluster * cluster = fEvent->GetCaloCluster(iClust);
    
    // Check which calorimeter information we want to keep.
    
    if(fCaloFilter!=kBoth)
    {
      if     (fCaloFilter==kPHOS  && cluster->IsEMCAL()) continue ;
      else if(fCaloFilter==kEMCAL && cluster->IsPHOS())  continue ;
    }  
    
    // Get original residuals, in case of previous recalculation, reset them 
    Float_t dR = cluster->GetTrackDx();
    Float_t dZ = cluster->GetTrackDz();
    
    AliDebug(2,Form("Original residuals : dZ %f, dR %f",dZ, dR));
    
    //--------------------------------------------------------------
    // If EMCAL and corrections done, get the new matching parameters, do not copy noisy clusters
    if(cluster->IsEMCAL() && fCorrect)
    {
      AliDebug(2,Form("Check cluster %d for bad channels and close to border",cluster->GetID()));
      
      if(fEMCALRecoUtils->ClusterContainsBadChannel(fEMCALGeo,cluster->GetCellsAbsId(), cluster->GetNCells())) continue;	
      
      if(fEMCALRecoUtils->IsExoticCluster(cluster, InputEvent()->GetEMCALCells(),InputEvent()->GetBunchCrossNumber())) continue;	
      
      fEMCALRecoUtils->GetMatchedResiduals(cluster->GetID(),dR,dZ);
      cluster->SetTrackDistance(dR,dZ);
    }
    
    AliDebug(2,Form("EMCAL? %d, PHOS? %d Track-Cluster Residuals : dZ %f, dR %f",cluster->IsEMCAL(), cluster->IsPHOS(),dZ, dR));
    
    //--------------------------------------------------------------
    
    // Now fill AODs
    
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
    
    AliDebug(2,Form("Filter, aod     : i %d, E %f, dispersion %f, m02 %f, m20 %f",caloCluster->GetID(),caloCluster->E(),
                    caloCluster->GetDispersion(),caloCluster->GetM02(),caloCluster->GetM20()));
    caloCluster->GetPosition(posF);
    AliDebug(2,Form("Filter, aod     : i %d, x %f, y %f, z %f",caloCluster->GetID(), posF[0], posF[1], posF[2]));
    
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
/// AliAODCaloTrigger direct copy.
//__________________________________________________
void AliAnalysisTaskCaloFilter::FillAODCaloTrigger()
{
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
/// AOD Header copy.
//______________________________________________
void AliAnalysisTaskCaloFilter::FillAODHeader()
{
  AliAODHeader* header = dynamic_cast<AliAODHeader*>(AODEvent()->GetHeader());
  if(!header)
  {
    AliFatal("Not a standard AOD");
    return; // not needed but coverity complains
  }
  
  // Copy from AODs
  if(fAODEvent)
  {
    *header = *((AliAODHeader*)fAODEvent->GetHeader());
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
  
  Float_t diamxy[2]={(Float_t)fEvent->GetDiamondX(),(Float_t)fEvent->GetDiamondY()};
  Float_t diamcov[3];
  fEvent->GetDiamondCovXY(diamcov);
  header->SetDiamond(diamxy,diamcov);
  header->SetDiamondZ(fESDEvent->GetDiamondZ(),fESDEvent->GetSigma2DiamondZ());
}

//__________________________________________________
/// Copy AOD MC particles.
//__________________________________________________
void AliAnalysisTaskCaloFilter::FillAODMCParticles()
{
  if(!fFillMCParticles) return;
  
  TClonesArray* inMCParticles = (TClonesArray*) (fAODEvent  ->FindListObject("mcparticles"));
  TClonesArray* ouMCParticles = (TClonesArray*) ( AODEvent()->FindListObject("mcparticles"));
  
  if( inMCParticles &&  ouMCParticles ) new (ouMCParticles) TClonesArray(*inMCParticles);
}

//_____________________________________________
/// Copy AOD track.
//_____________________________________________
void AliAnalysisTaskCaloFilter::FillAODTracks()
{
  if(!fFillTracks) return;
  
  AliAODTrack* aodTrack(0x0);
  
  Double_t pos[3]   = { 0. };      
  Double_t covTr[21]= { 0. };
  //Double_t pid[10]  = { 0. };  
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
      AliAODTrack *track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(nTrack));
      if(!track) AliFatal("Not a standard AOD");
      
      // Select only hybrid tracks?
      if(fFillHybridTracks && !track->IsHybridGlobalConstrainedGlobal()) continue;
      
      // Remove PID object to save space
      //track->SetDetPID(0x0);
      
      //new((*ouTracks)[nCopyTrack++]) AliAODTrack(*track);
      
      track->GetPxPyPz(p);
      Bool_t isDCA  = track->GetPosition(pos);
      track->GetCovMatrix(covTr);
      //track->GetPID(pid);
      
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
                                                            // pid,
                                                            primVertex,
                                                            track->GetUsedForVtxFit(),
                                                            track->GetUsedForPrimVtxFit(),
                                                            (AliAODTrack::AODTrk_t) track->GetType(), 
                                                            track->GetFilterMap());
      
      
      aodTrack->SetPIDForTracking(track->GetPIDForTracking());
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
      aodTrack->SetTrackPhiEtaPtOnEMCal( track->GetTrackPhiOnEMCal(),  track->GetTrackPhiOnEMCal(),  track->GetTrackPtOnEMCal() );
          
    } 
    
    //printf("Final N tracks %d\n",nCopyTrack);
    
    return;
  }
}

//_________________________________________
/// Copy v0s (use if you know what you do,
/// it consumes quite a lot of memory).
//_________________________________________
void AliAnalysisTaskCaloFilter::FillAODv0s()
{
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
/// Copy VZERO.
//____________________________________________
void AliAnalysisTaskCaloFilter::FillAODVZERO()
{
  if(!fFillVZERO) return;
  
  AliAODVZERO* vzeroData = AODEvent()->GetVZEROData();
  
  if(fESDEvent) *vzeroData = *(fESDEvent->GetVZEROData());
  else          *vzeroData = *(fAODEvent->GetVZEROData());
}  

//_______________________________________________
/// Copy vertices.
//_______________________________________________
void AliAnalysisTaskCaloFilter::FillAODVertices()
{
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
/// Init analysis with configuration macro, if available.
/// Example of configuration macro in: macros/ConfigCaloFilter.C
/// But it is a legacy from first configuration macros,
/// the suggested macro is: macros/AddTaskCaloFilter.C
//____________________________________
void AliAnalysisTaskCaloFilter::Init()
{
  if(gROOT->LoadMacro(fConfigName) >=0)
  {
    AliInfo(Form("Configure analysis with %s",fConfigName.Data()));
    
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

    for(Int_t i = 0; i < 22; i++) fEMCALMatrix[i] = filter->fEMCALMatrix[i] ;
  }
} 

//_________________________________________
/// Print settings.
//_________________________________________
void AliAnalysisTaskCaloFilter::PrintInfo()
{
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
/// Init EMCal geometry and create the AOD MC particles branch.
//_______________________________________________________
void AliAnalysisTaskCaloFilter::UserCreateOutputObjects()
{
  fEMCALGeo =  AliEMCALGeometry::GetInstance(fEMCALGeoName) ;	
  
  if(fFillMCParticles)
  {
    TClonesArray * aodMCParticles = new TClonesArray("AliAODMCParticle",500);
		aodMCParticles->SetName("mcparticles");
		((AliAODHandler*)AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler())->AddBranch("TClonesArray", &aodMCParticles);
  }
}  

//____________________________________________________________
/// Execute analysis for current event.
/// Copy input ESD or AOD header, vertex, CaloClusters, CaloCells,
/// Tracks, vertex, V0 and MC Particles to output AOD.
//____________________________________________________________
void AliAnalysisTaskCaloFilter::UserExec(Option_t */*option*/)
{
  AliDebug(1,Form("Analysing event # %d", (Int_t)Entry()));
  
  fEvent    = InputEvent();
  fAODEvent = dynamic_cast<AliAODEvent*> (fEvent);  
  fESDEvent = dynamic_cast<AliESDEvent*> (fEvent);
  
  if(!fEvent) 
  {
    AliInfo("This event does not contain Input?");
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

