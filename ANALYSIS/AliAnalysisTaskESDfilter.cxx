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

/* $Id$ */
 
#include <TChain.h>
#include <TTree.h>
#include <TList.h>
#include <TArrayI.h>
#include <TRandom.h>
#include <TParticle.h>
#include <TFile.h>

#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDRun.h"
#include "AliStack.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisFilter.h"
#include "AliESDMuonTrack.h"
#include "AliESDVertex.h"
#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliESDv0.h"
#include "AliESDkink.h"
#include "AliESDcascade.h"
#include "AliESDPmdTrack.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliMultiplicity.h"
#include "AliLog.h"
#include "AliCodeTimer.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "Riostream.h"

ClassImp(AliAnalysisTaskESDfilter)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskESDfilter::AliAnalysisTaskESDfilter():
    AliAnalysisTaskSE(),
    fTrackFilter(0x0),
    fKinkFilter(0x0),
    fV0Filter(0x0),
    fCascadeFilter(0x0),
    fHighPthreshold(0),
    fPtshape(0x0),
    fEnableFillAOD(kTRUE),
    fUsedTrack(0x0),
    fUsedKink(0x0),
    fUsedV0(0x0),
    fAODTrackRefs(0x0),
    fAODV0VtxRefs(0x0),
    fAODV0Refs(0x0),
    fMChandler(0x0),
    fNumberOfTracks(0),
    fNumberOfPositiveTracks(0),
    fNumberOfV0s(0),
    fNumberOfVertices(0),
    fNumberOfCascades(0),
    fNumberOfKinks(0),
    fOldESDformat(kFALSE),
    fPrimaryVertex(0x0),
  fTPCConstrainedFilterMask(0),
  fHybridFilterMaskTPCCG(0),
  fWriteHybridTPCCOnly(kFALSE),
  fGlobalConstrainedFilterMask(0),
  fHybridFilterMaskGCG(0),
  fWriteHybridGCOnly(kFALSE),
    fIsVZEROEnabled(kTRUE),
    fIsZDCEnabled(kTRUE),
    fAreCascadesEnabled(kTRUE),
    fAreV0sEnabled(kTRUE),
    fAreKinksEnabled(kTRUE),
    fAreTracksEnabled(kTRUE),
    fArePmdClustersEnabled(kTRUE),
    fAreCaloClustersEnabled(kTRUE),
    fAreEMCALCellsEnabled(kTRUE),
    fArePHOSCellsEnabled(kTRUE),
    fAreTrackletsEnabled(kTRUE),
    fESDpid(0x0),
    fIsPidOwner(kFALSE),
    fTimeZeroType(AliESDpid::kTOF_T0),
    fTPCaloneTrackCuts(0)
{
  // Default constructor
}

//______________________________________________________________________________
AliAnalysisTaskESDfilter::AliAnalysisTaskESDfilter(const char* name):
    AliAnalysisTaskSE(name),
    fTrackFilter(0x0),
    fKinkFilter(0x0),
    fV0Filter(0x0),
    fCascadeFilter(0x0),
    fHighPthreshold(0),
    fPtshape(0x0),
    fEnableFillAOD(kTRUE),
    fUsedTrack(0x0),
    fUsedKink(0x0),
    fUsedV0(0x0),
    fAODTrackRefs(0x0),
    fAODV0VtxRefs(0x0),
    fAODV0Refs(0x0),
    fMChandler(0x0),
    fNumberOfTracks(0),
    fNumberOfPositiveTracks(0),
    fNumberOfV0s(0),
    fNumberOfVertices(0),
    fNumberOfCascades(0),
    fNumberOfKinks(0),
    fOldESDformat(kFALSE),
    fPrimaryVertex(0x0),
  fTPCConstrainedFilterMask(0),
  fHybridFilterMaskTPCCG(0),
  fWriteHybridTPCCOnly(kFALSE),
  fGlobalConstrainedFilterMask(0),
  fHybridFilterMaskGCG(0),
  fWriteHybridGCOnly(kFALSE),
    fIsVZEROEnabled(kTRUE),
    fIsZDCEnabled(kTRUE),
    fAreCascadesEnabled(kTRUE),
    fAreV0sEnabled(kTRUE),
    fAreKinksEnabled(kTRUE),
    fAreTracksEnabled(kTRUE),
    fArePmdClustersEnabled(kTRUE),
    fAreCaloClustersEnabled(kTRUE),
    fAreEMCALCellsEnabled(kTRUE),
    fArePHOSCellsEnabled(kTRUE),
    fAreTrackletsEnabled(kTRUE),
    fESDpid(0x0),
    fIsPidOwner(kFALSE),
    fTimeZeroType(AliESDpid::kTOF_T0),
    fTPCaloneTrackCuts(0)
{
  // Constructor
}
AliAnalysisTaskESDfilter::~AliAnalysisTaskESDfilter(){
    if(fIsPidOwner) delete fESDpid;
}
//______________________________________________________________________________
void AliAnalysisTaskESDfilter::UserCreateOutputObjects()
{
  //
  // Create Output Objects conenct filter to outputtree
  // 
  if(OutputTree())
  {
    OutputTree()->GetUserInfo()->Add(fTrackFilter);
  }
  else
  {
    AliError("No OutputTree() for adding the track filter");
  }
  fTPCaloneTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::Init()
{
  // Initialization
  if (fDebug > 1) AliInfo("Init() \n");
  // Call configuration file
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::PrintTask(Option_t *option, Int_t indent) const
{
  AliInfo("");
  
  AliAnalysisTaskSE::PrintTask(option,indent);
  
  TString spaces(' ',indent+3);
  
  cout << spaces.Data() << Form("Cascades     are %s",fAreCascadesEnabled ? "ENABLED":"DISABLED") << endl;
  cout << spaces.Data() << Form("V0s          are %s",fAreV0sEnabled ? "ENABLED":"DISABLED") << endl;
  cout << spaces.Data() << Form("Kinks        are %s",fAreKinksEnabled ? "ENABLED":"DISABLED") << endl;
  cout << spaces.Data() << Form("Tracks       are %s",fAreTracksEnabled ? "ENABLED":"DISABLED") << endl;
  cout << spaces.Data() << Form("PmdClusters  are %s",fArePmdClustersEnabled ? "ENABLED":"DISABLED") << endl;
  cout << spaces.Data() << Form("CaloClusters are %s",fAreCaloClustersEnabled ? "ENABLED":"DISABLED") << endl;
  cout << spaces.Data() << Form("EMCAL cells  are %s",fAreEMCALCellsEnabled ? "ENABLED":"DISABLED") << endl;
  cout << spaces.Data() << Form("Tracklets    are %s",fAreTrackletsEnabled ? "ENABLED":"DISABLED") << endl;  
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::UserExec(Option_t */*option*/)
{
// Execute analysis for current event
//
					    
  Long64_t ientry = Entry();
  
  if (fDebug > 0) {
      printf("Filter: Analysing event # %5d\n", (Int_t) ientry);
      if (fHighPthreshold == 0) AliInfo("detector PID signals are stored in each track");
      if (!fPtshape) AliInfo("detector PID signals are not stored below the pt threshold");
  }
  // Filters must explicitely enable AOD filling in their UserExec (AG)
  if (!AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()) AliFatal("Cannot run ESD filter without an output event handler");
  if(fEnableFillAOD) {
     AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
     AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillExtension(kTRUE);
  }   
  ConvertESDtoAOD();
}

//______________________________________________________________________________
TClonesArray& AliAnalysisTaskESDfilter::Cascades()
{
  return *(AODEvent()->GetCascades());
}

//______________________________________________________________________________
TClonesArray& AliAnalysisTaskESDfilter::Tracks()
{
  return *(AODEvent()->GetTracks());
}

//______________________________________________________________________________
TClonesArray& AliAnalysisTaskESDfilter::V0s()
{
  return *(AODEvent()->GetV0s());
}

//______________________________________________________________________________
TClonesArray& AliAnalysisTaskESDfilter::Vertices()
{
  return *(AODEvent()->GetVertices());
}

//______________________________________________________________________________
AliAODHeader* AliAnalysisTaskESDfilter::ConvertHeader(const AliESDEvent& esd)
{
  AliCodeTimerAuto("",0);
  
  AliAODHeader* header = AODEvent()->GetHeader();
  
  header->SetRunNumber(esd.GetRunNumber());
  header->SetOfflineTrigger(fInputHandler->IsEventSelected()); // propagate the decision of the physics selection

  TTree* tree = fInputHandler->GetTree();
  if (tree) {
    TFile* file = tree->GetCurrentFile();
    if (file) header->SetESDFileName(file->GetName());
  }
  
  if (fOldESDformat) {
    header->SetBunchCrossNumber(0);
    header->SetOrbitNumber(0);
    header->SetPeriodNumber(0);
    header->SetEventType(0);
    header->SetMuonMagFieldScale(-999.);
    header->SetCentrality(0);       
    header->SetEventplane(0);
  } else {
    header->SetBunchCrossNumber(esd.GetBunchCrossNumber());
    header->SetOrbitNumber(esd.GetOrbitNumber());
    header->SetPeriodNumber(esd.GetPeriodNumber());
    header->SetEventType(esd.GetEventType());
    
    header->SetEventNumberESDFile(esd.GetHeader()->GetEventNumberInFile());
    if(const_cast<AliESDEvent&>(esd).GetCentrality()){
      header->SetCentrality(const_cast<AliESDEvent&>(esd).GetCentrality());
    }
    else{
      header->SetCentrality(0);
    }
    if(const_cast<AliESDEvent&>(esd).GetEventplane()){
      header->SetEventplane(const_cast<AliESDEvent&>(esd).GetEventplane());
    }
    else{
      header->SetEventplane(0);
    }
  }
  
  // Trigger
  header->SetFiredTriggerClasses(esd.GetFiredTriggerClasses());
  header->SetTriggerMask(esd.GetTriggerMask()); 
  header->SetTriggerCluster(esd.GetTriggerCluster());
  header->SetL0TriggerInputs(esd.GetHeader()->GetL0TriggerInputs());    
  header->SetL1TriggerInputs(esd.GetHeader()->GetL1TriggerInputs());    
  header->SetL2TriggerInputs(esd.GetHeader()->GetL2TriggerInputs());    
  
  header->SetMagneticField(esd.GetMagneticField());
  header->SetMuonMagFieldScale(esd.GetCurrentDip()/6000.);
  header->SetZDCN1Energy(esd.GetZDCN1Energy());
  header->SetZDCP1Energy(esd.GetZDCP1Energy());
  header->SetZDCN2Energy(esd.GetZDCN2Energy());
  header->SetZDCP2Energy(esd.GetZDCP2Energy());
  header->SetZDCEMEnergy(esd.GetZDCEMEnergy(0),esd.GetZDCEMEnergy(1));
  
  // ITS Cluster Multiplicty
  const AliMultiplicity *mult = esd.GetMultiplicity();
  for (Int_t ilay = 0; ilay < 6; ilay++) header->SetITSClusters(ilay, mult->GetNumberOfITSClusters(ilay));
  
  // TPC only Reference Multiplicty
  Int_t refMult  = fTPCaloneTrackCuts ? (Short_t)fTPCaloneTrackCuts->GetReferenceMultiplicity(&esd, kTRUE) : -1;
  header->SetTPConlyRefMultiplicity(refMult);
  
  //
  Float_t diamxy[2]={esd.GetDiamondX(),esd.GetDiamondY()};
  Float_t diamcov[3]; 
  esd.GetDiamondCovXY(diamcov);
  header->SetDiamond(diamxy,diamcov);
  header->SetDiamondZ(esd.GetDiamondZ(),esd.GetSigma2DiamondZ());
  
  return header;
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::ConvertCascades(const AliESDEvent& esd) 
{
  // Convert the cascades part of the ESD.
  // Return the number of cascades
 
  AliCodeTimerAuto("",0);
  
  // Create vertices starting from the most complex objects
  Double_t chi2 = 0.;
  
  const AliESDVertex* vtx = esd.GetPrimaryVertex();
  Double_t pos[3] = { 0. };
  Double_t covVtx[6] = { 0. };
  Double_t momBach[3]={0.};
  Double_t covTr[21]={0.};
  Double_t pid[10]={0.};
  AliAODPid* detpid(0x0);
  AliAODVertex* vV0FromCascade(0x0);
  AliAODv0* aodV0(0x0);
  AliAODcascade* aodCascade(0x0);
  AliAODTrack* aodTrack(0x0);
  Double_t momPos[3]={0.};
  Double_t momNeg[3] = { 0. };
  Double_t momPosAtV0vtx[3]={0.};
  Double_t momNegAtV0vtx[3]={0.};

  TClonesArray& verticesArray = Vertices();
  TClonesArray& tracksArray = Tracks();
  TClonesArray& cascadesArray = Cascades();
  
  // Cascades (Modified by A.Maire - February 2009)
  for (Int_t nCascade = 0; nCascade < esd.GetNumberOfCascades(); ++nCascade) {
    
    // 0- Preparation
    //
    AliESDcascade *esdCascade = esd.GetCascade(nCascade);
		Int_t  idxPosFromV0Dghter  = esdCascade->GetPindex();
		Int_t  idxNegFromV0Dghter  = esdCascade->GetNindex();
		Int_t  idxBachFromCascade  = esdCascade->GetBindex();
    
    AliESDtrack  *esdCascadePos  = esd.GetTrack( idxPosFromV0Dghter);
    AliESDtrack  *esdCascadeNeg  = esd.GetTrack( idxNegFromV0Dghter);
    AliESDtrack  *esdCascadeBach = esd.GetTrack( idxBachFromCascade);
    
    // Identification of the V0 within the esdCascade (via both daughter track indices)
    AliESDv0 * currentV0   = 0x0;
    Int_t      idxV0FromCascade = -1;
    
    for (Int_t iV0=0; iV0<esd.GetNumberOfV0s(); ++iV0) {
      
      currentV0 = esd.GetV0(iV0);
      Int_t posCurrentV0 = currentV0->GetPindex();
      Int_t negCurrentV0 = currentV0->GetNindex();
      
      if (posCurrentV0==idxPosFromV0Dghter && negCurrentV0==idxNegFromV0Dghter) {
        idxV0FromCascade = iV0;
        break;
      }
    }
    
    if(idxV0FromCascade < 0){
      printf("Cascade - no matching for the V0 (index V0 = -1) ! Skip ... \n");
      continue;
    }// a priori, useless check, but safer ... in case of pb with tracks "out of bounds"
    
    AliESDv0 *esdV0FromCascade   = esd.GetV0(idxV0FromCascade);
        
    // 1 - Cascade selection 
    
    //	AliESDVertex *esdPrimVtx = new AliESDVertex(*(esd.GetPrimaryVertex()));
    // 	TList cascadeObjects;
    // 	cascadeObjects.AddAt(esdV0FromCascade, 0);
    // 	cascadeObjects.AddAt(esdCascadePos,    1);
    // 	cascadeObjects.AddAt(esdCascadeNeg,    2);
    // 	cascadeObjects.AddAt(esdCascade,       3);
    // 	cascadeObjects.AddAt(esdCascadeBach,   4);
    // 	cascadeObjects.AddAt(esdPrimVtx,       5);
    // 
    // 	UInt_t selectCascade = 0;
    // 	if (fCascadeFilter) {
    // 	  // selectCascade = fCascadeFilter->IsSelected(&cascadeObjects); 
    // 	  	// FIXME AliESDCascadeCuts to be implemented ...
    // 
    // 		// Here we may encounter a moot point at the V0 level 
    // 		// between the cascade selections and the V0 ones :
    // 		// the V0 selected along with the cascade (secondary V0) may 
    // 		// usually be removed from the dedicated V0 selections (prim V0) ...
    // 		// -> To be discussed !
    // 
    // 	  // this is a little awkward but otherwise the 
    // 	  // list wants to access the pointer (delete it) 
    // 	  // again when going out of scope
    // 	  delete cascadeObjects.RemoveAt(5); // esdPrimVtx created via copy construct
    // 	  esdPrimVtx = 0;
    // 	  if (!selectCascade) 
    // 	    continue;
    // 	}
    // 	else{
    // 	  delete cascadeObjects.RemoveAt(5); // esdPrimVtx created via copy construct
    // 	  esdPrimVtx = 0;
    // 	}
    
    // 2 - Add the cascade vertex
    
    esdCascade->GetXYZcascade(pos[0], pos[1], pos[2]);
    esdCascade->GetPosCovXi(covVtx);
    chi2 = esdCascade->GetChi2Xi(); 
    
    AliAODVertex *vCascade = new(verticesArray[fNumberOfVertices++]) AliAODVertex( pos,
                                                                     covVtx,
                                                                     chi2, // FIXME = Chi2/NDF will be needed
                                                                     fPrimaryVertex,
                                                                     nCascade, // id
                                                                     AliAODVertex::kCascade);
    fPrimaryVertex->AddDaughter(vCascade);
    
//    if (fDebug > 2) {
//      printf("---- Cascade / Cascade Vertex (AOD) : \n");
//      vCascade->Print();
//    }
    
    if(esd.GetTOFHeader() && fIsPidOwner) fESDpid->SetTOFResponse(const_cast<AliESDEvent*>(&esd), (AliESDpid::EStartTimeType_t)fTimeZeroType); //in case of AOD production strating form LHC10e without Tender. 


    // 3 - Add the bachelor track from the cascade
    
    if (!fUsedTrack[idxBachFromCascade]) {
      
      esdCascadeBach->GetPxPyPz(momBach);
      esdCascadeBach->GetXYZ(pos);
      esdCascadeBach->GetCovarianceXYZPxPyPz(covTr);
      esdCascadeBach->GetESDpid(pid);
      
	    fUsedTrack[idxBachFromCascade] = kTRUE;
	    UInt_t selectInfo = 0;
	    if (fTrackFilter) selectInfo = fTrackFilter->IsSelected(esdCascadeBach);
	    if (fMChandler) fMChandler->SelectParticle(esdCascadeBach->GetLabel());
	    aodTrack = new(tracksArray[fNumberOfTracks++]) AliAODTrack(esdCascadeBach->GetID(),
                                                                           esdCascadeBach->GetLabel(), 
                                                                           momBach, 
                                                                           kTRUE,
                                                                           pos,
                                                                           kFALSE, // Why kFALSE for "isDCA" ? FIXME
                                                                           covTr, 
                                                                           (Short_t)esdCascadeBach->GetSign(),
                                                                           esdCascadeBach->GetITSClusterMap(), 
                                                                           pid,
                                                                           vCascade,
                                                                           kTRUE,  // usedForVtxFit = kFALSE ? FIXME
                                                                           vtx->UsesTrack(esdCascadeBach->GetID()),
                                                                           AliAODTrack::kSecondary,
                                                                           selectInfo);
	    aodTrack->SetTPCClusterMap(esdCascadeBach->GetTPCClusterMap());
	    aodTrack->SetTPCSharedMap (esdCascadeBach->GetTPCSharedMap());
	    aodTrack->SetChi2perNDF(Chi2perNDF(esdCascadeBach));
	    aodTrack->SetTPCPointsF(esdCascadeBach->GetTPCNclsF());
	    fAODTrackRefs->AddAt(aodTrack,idxBachFromCascade);
	    
	    if (esdCascadeBach->GetSign() > 0) ++fNumberOfPositiveTracks;
	    aodTrack->ConvertAliPIDtoAODPID();
	    aodTrack->SetFlags(esdCascadeBach->GetStatus());
      SetAODPID(esdCascadeBach,aodTrack,detpid);
    }
    else {
	    aodTrack = static_cast<AliAODTrack*>( fAODTrackRefs->At(idxBachFromCascade) );
    }
    
    vCascade->AddDaughter(aodTrack);
    
//    if (fDebug > 4) {
//      printf("---- Cascade / bach dghter : \n");
//      aodTrack->Print();
//    }
    
    
    // 4 - Add the V0 from the cascade. 
    // = V0vtx + both pos and neg daughter tracks + the aodV0 itself
    //
    
    if ( !fUsedV0[idxV0FromCascade] ) {
      // 4.A - if VO structure hasn't been created yet
      
      // 4.A.1 - Create the V0 vertex of the cascade
      
      esdV0FromCascade->GetXYZ(pos[0], pos[1], pos[2]);
      esdV0FromCascade->GetPosCov(covVtx);
      chi2 = esdV0FromCascade->GetChi2V0();  // = chi2/NDF since NDF = 2*2-3 ?
			
      vV0FromCascade = new(verticesArray[fNumberOfVertices++]) AliAODVertex(pos,
                                                                         covVtx,
                                                                         chi2,
                                                                         vCascade,
                                                                         idxV0FromCascade, //id of ESDv0
                                                                         AliAODVertex::kV0);
      // Note:
      //    one V0 can be used by several cascades.
      // So, one AOD V0 vtx can have several parent vtx.
      // This is not directly allowed by AliAODvertex.
      // Setting the parent vtx (here = param "vCascade") doesn't lead to a crash
      // but to a problem of consistency within AODEvent.
      // -> See below paragraph 4.B, for the proposed treatment of such a case.
      
      // Add the vV0FromCascade to the aodVOVtxRefs
      fAODV0VtxRefs->AddAt(vV0FromCascade,idxV0FromCascade);
      
      
      // 4.A.2 - Add the positive tracks from the V0
      
      esdCascadePos->GetPxPyPz(momPos);
      esdCascadePos->GetXYZ(pos);
      esdCascadePos->GetCovarianceXYZPxPyPz(covTr);
      esdCascadePos->GetESDpid(pid);
      
      
      if (!fUsedTrack[idxPosFromV0Dghter]) {
        fUsedTrack[idxPosFromV0Dghter] = kTRUE;
        
        UInt_t selectInfo = 0;
        if (fTrackFilter) selectInfo = fTrackFilter->IsSelected(esdCascadePos);
        if(fMChandler) fMChandler->SelectParticle(esdCascadePos->GetLabel());
        aodTrack = new(tracksArray[fNumberOfTracks++]) 
        AliAODTrack(  esdCascadePos->GetID(),
                    esdCascadePos->GetLabel(), 
                    momPos, 
                    kTRUE,
                    pos,
                    kFALSE, // Why kFALSE for "isDCA" ? FIXME
                    covTr, 
                    (Short_t)esdCascadePos->GetSign(),
                    esdCascadePos->GetITSClusterMap(), 
                    pid,
                    vV0FromCascade,
                    kTRUE,  // usedForVtxFit = kFALSE ? FIXME
                    vtx->UsesTrack(esdCascadePos->GetID()),
                    AliAODTrack::kSecondary,
                    selectInfo);
        aodTrack->SetTPCClusterMap(esdCascadePos->GetTPCClusterMap());
        aodTrack->SetTPCSharedMap (esdCascadePos->GetTPCSharedMap());
        aodTrack->SetChi2perNDF(Chi2perNDF(esdCascadePos));
	aodTrack->SetTPCPointsF(esdCascadePos->GetTPCNclsF());
        fAODTrackRefs->AddAt(aodTrack,idxPosFromV0Dghter);
        
        if (esdCascadePos->GetSign() > 0) ++fNumberOfPositiveTracks;
        aodTrack->ConvertAliPIDtoAODPID();
        aodTrack->SetFlags(esdCascadePos->GetStatus());
        SetAODPID(esdCascadePos,aodTrack,detpid);
      }
      else {
        aodTrack = static_cast<AliAODTrack*>(fAODTrackRefs->At(idxPosFromV0Dghter));
      }
      vV0FromCascade->AddDaughter(aodTrack);
      
      
      // 4.A.3 - Add the negative tracks from the V0
      
      esdCascadeNeg->GetPxPyPz(momNeg);
      esdCascadeNeg->GetXYZ(pos);
      esdCascadeNeg->GetCovarianceXYZPxPyPz(covTr);
      esdCascadeNeg->GetESDpid(pid);
      
      
      if (!fUsedTrack[idxNegFromV0Dghter]) {
        fUsedTrack[idxNegFromV0Dghter] = kTRUE;
        
        UInt_t selectInfo = 0;
        if (fTrackFilter) selectInfo = fTrackFilter->IsSelected(esdCascadeNeg);
        if(fMChandler)fMChandler->SelectParticle(esdCascadeNeg->GetLabel());
        aodTrack = new(tracksArray[fNumberOfTracks++]) AliAODTrack(  esdCascadeNeg->GetID(),
                                                      esdCascadeNeg->GetLabel(),
                                                      momNeg,
                                                      kTRUE,
                                                      pos,
                                                      kFALSE, // Why kFALSE for "isDCA" ? FIXME
                                                      covTr, 
                                                      (Short_t)esdCascadeNeg->GetSign(),
                                                      esdCascadeNeg->GetITSClusterMap(), 
                                                      pid,
                                                      vV0FromCascade,
                                                      kTRUE,  // usedForVtxFit = kFALSE ? FIXME
                                                      vtx->UsesTrack(esdCascadeNeg->GetID()),
                                                      AliAODTrack::kSecondary,
                                                      selectInfo);
        aodTrack->SetTPCClusterMap(esdCascadeNeg->GetTPCClusterMap());
        aodTrack->SetTPCSharedMap (esdCascadeNeg->GetTPCSharedMap());
        aodTrack->SetChi2perNDF(Chi2perNDF(esdCascadeNeg));
	aodTrack->SetTPCPointsF(esdCascadeNeg->GetTPCNclsF());
        fAODTrackRefs->AddAt(aodTrack,idxNegFromV0Dghter);
        
        if (esdCascadeNeg->GetSign() > 0) ++fNumberOfPositiveTracks;
        aodTrack->ConvertAliPIDtoAODPID();
        aodTrack->SetFlags(esdCascadeNeg->GetStatus());
        SetAODPID(esdCascadeNeg,aodTrack,detpid);
      }
      else {
        aodTrack = static_cast<AliAODTrack*>(fAODTrackRefs->At(idxNegFromV0Dghter));
      }
      
      vV0FromCascade->AddDaughter(aodTrack);
      
			
      // 4.A.4 - Add the V0 from cascade to the V0 array
      
      Double_t  dcaV0Daughters      = esdV0FromCascade->GetDcaV0Daughters();
      Double_t  dcaV0ToPrimVertex   = esdV0FromCascade->GetD( esd.GetPrimaryVertex()->GetX(),
                                                             esd.GetPrimaryVertex()->GetY(),
                                                             esd.GetPrimaryVertex()->GetZ() );
      esdV0FromCascade->GetPPxPyPz( momPosAtV0vtx[0],momPosAtV0vtx[1],momPosAtV0vtx[2] ); 
      esdV0FromCascade->GetNPxPyPz( momNegAtV0vtx[0],momNegAtV0vtx[1],momNegAtV0vtx[2] ); 
      
      Double_t dcaDaughterToPrimVertex[2] = { 999., 999.}; // ..[0] = DCA in (x,y) for Pos and ..[1] = Neg
      dcaDaughterToPrimVertex[0] = TMath::Abs(esdCascadePos->GetD(	esd.GetPrimaryVertex()->GetX(),
                                                                  esd.GetPrimaryVertex()->GetY(),
                                                                  esd.GetMagneticField())        );
      dcaDaughterToPrimVertex[1] = TMath::Abs(esdCascadeNeg->GetD(	esd.GetPrimaryVertex()->GetX(),
                                                                  esd.GetPrimaryVertex()->GetY(),
                                                                  esd.GetMagneticField())        );
      
      aodV0 = new(V0s()[fNumberOfV0s++]) AliAODv0( vV0FromCascade, 
                                                  dcaV0Daughters,
                                                  dcaV0ToPrimVertex, 
                                                  momPosAtV0vtx, 
                                                  momNegAtV0vtx, 
                                                  dcaDaughterToPrimVertex); 
      // set the aod v0 on-the-fly status
      aodV0->SetOnFlyStatus(esdV0FromCascade->GetOnFlyStatus());
      
      // Add the aodV0 to the aodVORefs
      fAODV0Refs->AddAt(aodV0,idxV0FromCascade);
      
      fUsedV0[idxV0FromCascade] = kTRUE;
      
    } else { 
      // 4.B - if V0 structure already used
      
      // Note :
      //    one V0 can be used by several cascades (frequent in PbPb evts) : 
      // same V0 which used but attached to different bachelor tracks
      // -> aodVORefs and fAODV0VtxRefs are needed.
      // Goal : avoid a redundancy of the info in "Vertices" and "v0s" clones array.
      
      vV0FromCascade = static_cast<AliAODVertex*>( fAODV0VtxRefs->At(idxV0FromCascade) );
      aodV0          = static_cast<AliAODv0*>    ( fAODV0Refs   ->At(idxV0FromCascade) );
      
      // - Treatment of the parent for such a "re-used" V0 :
      // Insert the cascade that reuses the V0 vertex in the lineage chain
      // Before : vV0 -> vCascade1 -> vPrimary
      //  - Hyp : cascade2 uses the same V0 as cascade1
      //  After :  vV0 -> vCascade2 -> vCascade1 -> vPrimary
      
      AliAODVertex *vCascadePreviousParent = static_cast<AliAODVertex*> (vV0FromCascade->GetParent());
			vV0FromCascade->SetParent(vCascade);
			vCascade      ->SetParent(vCascadePreviousParent);
      
//      if(fDebug > 2)	
//        printf("---- Cascade / Lineage insertion\n"
//               "Parent of V0 vtx                 = Cascade vtx %p\n"
//               "Parent of the cascade vtx        = Cascade vtx %p\n"
//               "Parent of the parent cascade vtx = Cascade vtx %p\n", 
//               static_cast<void*> (vV0FromCascade->GetParent()),
//               static_cast<void*> (vCascade->GetParent()),
//               static_cast<void*> (vCascadePreviousParent->GetParent()) );
      
    }// end if V0 structure already used
    
//    if (fDebug > 2) {
//      printf("---- Cascade / V0 vertex: \n");
//      vV0FromCascade->Print();
//    }
//    
//    if (fDebug > 4) {
//      printf("---- Cascade / pos dghter : \n");
//			aodTrack->Print();
//      printf("---- Cascade / neg dghter : \n");
//			aodTrack->Print();
//      printf("---- Cascade / aodV0 : \n");
//			aodV0->Print();
//    }
    
    // In any case (used V0 or not), add the V0 vertex to the cascade one.
    vCascade->AddDaughter(vV0FromCascade);	
    
		
    // 5 - Add the primary track of the cascade (if any)
    
    
    // 6 - Add the cascade to the AOD array of cascades
    
    Double_t dcaBachToPrimVertexXY = TMath::Abs(esdCascadeBach->GetD(esd.GetPrimaryVertex()->GetX(),
                                                                     esd.GetPrimaryVertex()->GetY(),
                                                                     esd.GetMagneticField())        );
    
    Double_t momBachAtCascadeVtx[3]={0.};

    esdCascade->GetBPxPyPz(momBachAtCascadeVtx[0], momBachAtCascadeVtx[1], momBachAtCascadeVtx[2]);
    
    aodCascade = new(cascadesArray[fNumberOfCascades++]) AliAODcascade(  vCascade,
                                                          esdCascade->Charge(),
                                                          esdCascade->GetDcaXiDaughters(),
                                                          -999.,
                                                          // DCAXiToPrimVtx -> needs to be calculated   ----|
                                                          // doesn't exist at ESD level;
                                                          // See AODcascade::DcaXiToPrimVertex(Double, Double, Double)
                                                          dcaBachToPrimVertexXY,
                                                          momBachAtCascadeVtx,
                                                          *aodV0);
    
    if (fDebug > 10) {
      printf("---- Cascade / AOD cascade : \n\n");
      aodCascade->PrintXi(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    }
    
  } // end of the loop on cascades
  
  Cascades().Expand(fNumberOfCascades);
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::ConvertV0s(const AliESDEvent& esd)
{
  // Access to the AOD container of V0s
  
  AliCodeTimerAuto("",0);

  //
  // V0s
  //
  
  Double_t pos[3] = { 0. };      
  Double_t chi2(0.0);
  Double_t covVtx[6] = { 0. };
  Double_t momPos[3]={0.};
  Double_t covTr[21]={0.};
  Double_t pid[10]={0.};
  AliAODTrack* aodTrack(0x0);
  AliAODPid* detpid(0x0);
  Double_t momNeg[3]={0.};
  Double_t momPosAtV0vtx[3]={0.};
  Double_t momNegAtV0vtx[3]={0.};
    
  for (Int_t nV0 = 0; nV0 < esd.GetNumberOfV0s(); ++nV0) 
  {
    if (fUsedV0[nV0]) continue; // skip if already added to the AOD
    
    AliESDv0 *v0 = esd.GetV0(nV0);
    Int_t posFromV0 = v0->GetPindex();
    Int_t negFromV0 = v0->GetNindex();
    
    // V0 selection 
    //
    AliESDVertex *esdVtx   = new AliESDVertex(*(esd.GetPrimaryVertex()));
    AliESDtrack  *esdV0Pos = esd.GetTrack(posFromV0);
    AliESDtrack  *esdV0Neg = esd.GetTrack(negFromV0);
    TList v0objects;
    v0objects.AddAt(v0,                      0);
    v0objects.AddAt(esdV0Pos,                1);
    v0objects.AddAt(esdV0Neg,                2);
    v0objects.AddAt(esdVtx,                  3);
    UInt_t selectV0 = 0;
    if (fV0Filter) {
      selectV0 = fV0Filter->IsSelected(&v0objects);
      // this is a little awkward but otherwise the 
      // list wants to access the pointer (delete it) 
      // again when going out of scope
      delete v0objects.RemoveAt(3); // esdVtx created via copy construct
      esdVtx = 0;
      if (!selectV0) 
        continue;
    }
    else{
      delete v0objects.RemoveAt(3); // esdVtx created via copy construct
      esdVtx = 0;
    }
    
    v0->GetXYZ(pos[0], pos[1], pos[2]);
    
    if (!fOldESDformat) {
	    chi2 = v0->GetChi2V0(); // = chi2/NDF since NDF = 2*2-3
	    v0->GetPosCov(covVtx);
    } else {
	    chi2 = -999.;
	    for (Int_t i = 0; i < 6; i++)  covVtx[i] = 0.;
    }
    
    
    AliAODVertex * vV0 = 
	  new(Vertices()[fNumberOfVertices++]) AliAODVertex(pos,
                                                      covVtx,
                                                      chi2,
                                                      fPrimaryVertex,
                                                      nV0,
                                                      AliAODVertex::kV0);
    fPrimaryVertex->AddDaughter(vV0);
    
    
    // Add the positive tracks from the V0
    

    esdV0Pos->GetPxPyPz(momPos);
    esdV0Pos->GetXYZ(pos);
    esdV0Pos->GetCovarianceXYZPxPyPz(covTr);
    esdV0Pos->GetESDpid(pid);
    
    const AliESDVertex *vtx = esd.GetPrimaryVertex();

    if (!fUsedTrack[posFromV0]) {
	    fUsedTrack[posFromV0] = kTRUE;
	    UInt_t selectInfo = 0;
	    if (fTrackFilter) selectInfo = fTrackFilter->IsSelected(esdV0Pos);
	    if(fMChandler)fMChandler->SelectParticle(esdV0Pos->GetLabel());
	    aodTrack = new(Tracks()[fNumberOfTracks++]) AliAODTrack(esdV0Pos->GetID(),
                                                    esdV0Pos->GetLabel(), 
                                                    momPos, 
                                                    kTRUE,
                                                    pos,
                                                    kFALSE,
                                                    covTr, 
                                                    (Short_t)esdV0Pos->GetSign(),
                                                    esdV0Pos->GetITSClusterMap(), 
                                                    pid,
                                                    vV0,
                                                    kTRUE,  // check if this is right
                                                    vtx->UsesTrack(esdV0Pos->GetID()),
                                                    AliAODTrack::kSecondary,
                                                    selectInfo);
	    aodTrack->SetTPCClusterMap(esdV0Pos->GetTPCClusterMap());
	    aodTrack->SetTPCSharedMap (esdV0Pos->GetTPCSharedMap());
	    aodTrack->SetChi2perNDF(Chi2perNDF(esdV0Pos));
	    aodTrack->SetTPCPointsF(esdV0Pos->GetTPCNclsF());
	    fAODTrackRefs->AddAt(aodTrack,posFromV0);
	    //	    if (fDebug > 0) printf("-------------------Bo: pos track from original pt %.3f \n",aodTrack->Pt());
	    if (esdV0Pos->GetSign() > 0) ++fNumberOfPositiveTracks;
	    aodTrack->ConvertAliPIDtoAODPID();
	    aodTrack->SetFlags(esdV0Pos->GetStatus());
      SetAODPID(esdV0Pos,aodTrack,detpid);
    }
    else {
	    aodTrack = static_cast<AliAODTrack*>(fAODTrackRefs->At(posFromV0));
	    //	    if (fDebug > 0) printf("-------------------Bo pos track from refArray pt %.3f \n",aodTrack->Pt());
    }
    vV0->AddDaughter(aodTrack);
    
    // Add the negative tracks from the V0
    
    esdV0Neg->GetPxPyPz(momNeg);
    esdV0Neg->GetXYZ(pos);
    esdV0Neg->GetCovarianceXYZPxPyPz(covTr);
    esdV0Neg->GetESDpid(pid);
    
    if (!fUsedTrack[negFromV0]) {
	    fUsedTrack[negFromV0] = kTRUE;
	    UInt_t selectInfo = 0;
	    if (fTrackFilter) selectInfo = fTrackFilter->IsSelected(esdV0Neg);
	    if(fMChandler)fMChandler->SelectParticle(esdV0Neg->GetLabel());
	    aodTrack = new(Tracks()[fNumberOfTracks++]) AliAODTrack(esdV0Neg->GetID(),
                                                    esdV0Neg->GetLabel(),
                                                    momNeg,
                                                    kTRUE,
                                                    pos,
                                                    kFALSE,
                                                    covTr, 
                                                    (Short_t)esdV0Neg->GetSign(),
                                                    esdV0Neg->GetITSClusterMap(), 
                                                    pid,
                                                    vV0,
                                                    kTRUE,  // check if this is right
                                                    vtx->UsesTrack(esdV0Neg->GetID()),
                                                    AliAODTrack::kSecondary,
                                                    selectInfo);
	    aodTrack->SetTPCClusterMap(esdV0Neg->GetTPCClusterMap());
	    aodTrack->SetTPCSharedMap (esdV0Neg->GetTPCSharedMap());
	    aodTrack->SetChi2perNDF(Chi2perNDF(esdV0Neg));
	    aodTrack->SetTPCPointsF(esdV0Neg->GetTPCNclsF());
	    
	    fAODTrackRefs->AddAt(aodTrack,negFromV0);
	    //	    if (fDebug > 0) printf("-------------------Bo: neg track from original pt %.3f \n",aodTrack->Pt());
	    if (esdV0Neg->GetSign() > 0) ++fNumberOfPositiveTracks;
	    aodTrack->ConvertAliPIDtoAODPID();
	    aodTrack->SetFlags(esdV0Neg->GetStatus());
      SetAODPID(esdV0Neg,aodTrack,detpid);
    }
    else {
	    aodTrack = static_cast<AliAODTrack*>(fAODTrackRefs->At(negFromV0));
	    //	    if (fDebug > 0) printf("-------------------Bo neg track from refArray pt %.3f \n",aodTrack->Pt());
    }
    vV0->AddDaughter(aodTrack);
    
    
    // Add the V0 the V0 array as well
    
    Double_t  dcaV0Daughters      = v0->GetDcaV0Daughters();
    Double_t  dcaV0ToPrimVertex   = v0->GetD(esd.GetPrimaryVertex()->GetX(),
                                             esd.GetPrimaryVertex()->GetY(),
                                             esd.GetPrimaryVertex()->GetZ());
    
    v0->GetPPxPyPz(momPosAtV0vtx[0],momPosAtV0vtx[1],momPosAtV0vtx[2]); 
    v0->GetNPxPyPz(momNegAtV0vtx[0],momNegAtV0vtx[1],momNegAtV0vtx[2]); 
    
    Double_t dcaDaughterToPrimVertex[2] = { 999., 999.}; // ..[0] = DCA in (x,y) for Pos and ..[1] = Neg
    dcaDaughterToPrimVertex[0] = TMath::Abs(esdV0Pos->GetD(  esd.GetPrimaryVertex()->GetX(),
                                                           esd.GetPrimaryVertex()->GetY(),
                                                           esd.GetMagneticField()) );
    dcaDaughterToPrimVertex[1] = TMath::Abs(esdV0Neg->GetD(  esd.GetPrimaryVertex()->GetX(),
                                                           esd.GetPrimaryVertex()->GetY(),
                                                           esd.GetMagneticField()) );
    
    AliAODv0* aodV0 = new(V0s()[fNumberOfV0s++]) AliAODv0(vV0, 
                                                dcaV0Daughters,
                                                dcaV0ToPrimVertex,
                                                momPosAtV0vtx,
                                                momNegAtV0vtx,
                                                dcaDaughterToPrimVertex);
    
    // set the aod v0 on-the-fly status
    aodV0->SetOnFlyStatus(v0->GetOnFlyStatus());
  }//End of loop on V0s 
  
  V0s().Expand(fNumberOfV0s);	 
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::ConvertTPCOnlyTracks(const AliESDEvent& esd)
{
  // Convert TPC only tracks
  // Here we have wo hybrid appraoch to remove fakes
  // ******* ITSTPC ********
  // Uses a cut on the ITS properties to select global tracks
  // which are than marked as HybdridITSTPC for the remainder 
  // the TPC only tracks are flagged as HybridITSTPConly. 
  // Note, in order not to get fakes back in the TPC cuts, one needs 
  // two "ITS" cuts one tight (1) (to throw out fakes) and one lose (2) (to NOT flag the trakcs in the TPC only)
  // using cut number (3)
  // so fHybridFilterMask == (1)|(2) fTPCFilterMask = (3), Usercode needs to slect with mask = (1)|(3) and track->IsHybridITSTPC()
  // ******* TPC ********
  // Here, only TPC tracks are flagged that pass the tight ITS cuts and tracks that pass the TPC cuts and NOT the loose ITS cuts
  // the ITS cuts neeed to be added to the filter as extra cuts, since here the selections info is reset in the global and put to the TPC only track

  AliCodeTimerAuto("",0);
  
  // Loop over the tracks and extract and mask out all aod tracks that pass the selections for AODt racks
  for(int it = 0;it < fNumberOfTracks;++it)
  {
    AliAODTrack *tr = (AliAODTrack*)(Tracks().At(it));
    if(!tr)continue;
    UInt_t map = tr->GetFilterMap();
    if(map&fTPCConstrainedFilterMask){
      // we only reset the track select ionfo, no deletion...
      tr->SetFilterMap(map&~fTPCConstrainedFilterMask);
    }
    if(map&fHybridFilterMaskTPCCG){
      // this is one part of the hybrid tracks
      // the others not passing the selection will be TPC only selected below
      tr->SetIsHybridTPCConstrainedGlobal(kTRUE);
    }
  }
  // Loop over the ESD trcks and pick out the tracks passing TPC only cuts
  
  
  const AliESDVertex *vtxSPD = esd.GetPrimaryVertexSPD();
  const AliESDVertex *vtx = esd.GetPrimaryVertex();

  Double_t pos[3] = { 0. };      
  Double_t covTr[21]={0.};
  Double_t pid[10]={0.};  
  Double_t p[3] = { 0. };

  Double_t pDCA[3] = { 0. }; // momentum at DCA
  Double_t rDCA[3] = { 0. }; // position at DCA
  Float_t  dDCA[2] = {0.};    // DCA to the vertex d and z
  Float_t  cDCA[3] = {0.};    // covariance of impact parameters


  AliAODTrack* aodTrack(0x0);
  
  // account for change in pT after the constraint
  Float_t ptMax = 1E10;
  Float_t ptMin = 0;
  for(int i = 0;i<32;i++){
    if(fTPCConstrainedFilterMask&(1<<i)){
      AliESDtrackCuts*cuts = (AliESDtrackCuts*)fTrackFilter->GetCuts()->At(i);
      Float_t tmp1= 0,tmp2 = 0;
      cuts->GetPtRange(tmp1,tmp2);
      if(tmp1>ptMin)ptMin=tmp1;
      if(tmp2<ptMax)ptMax=tmp2;
    }
  } 

  for (Int_t nTrack = 0; nTrack < esd.GetNumberOfTracks(); ++nTrack) 
  {
    AliESDtrack* esdTrack = esd.GetTrack(nTrack); //carefull do not modify it othwise  need to work with a copy 
    
    UInt_t selectInfo = 0;
    Bool_t isHybridITSTPC = false;
    //
    // Track selection
    if (fTrackFilter) {
      selectInfo = fTrackFilter->IsSelected(esdTrack);
    }

    if(!(selectInfo&fHybridFilterMaskTPCCG)){
      // not already selected tracks, use second part of hybrid tracks
      isHybridITSTPC = true;
      // too save space one could only store these...
    }

    selectInfo &= fTPCConstrainedFilterMask;
    if (!selectInfo)continue;
    if (fWriteHybridTPCCOnly&&!isHybridITSTPC)continue; // write only complementary tracks
    // create a tpc only tracl
    AliESDtrack *track = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(&esd),esdTrack->GetID());
    if(!track) continue;
    
    if(track->Pt()>0.)
    {
      // only constrain tracks above threshold
      AliExternalTrackParam exParam;
      // take the B-field from the ESD, no 3D fieldMap available at this point
      Bool_t relate = false;
      relate = track->RelateToVertexTPC(vtxSPD,esd.GetMagneticField(),kVeryBig,&exParam);
      if(!relate){
        delete track;
        continue;
      }
      // fetch the track parameters at the DCA (unconstraint)
      if(track->GetTPCInnerParam()){
	track->GetTPCInnerParam()->GetPxPyPz(pDCA);
	track->GetTPCInnerParam()->GetXYZ(rDCA);
      }
      // get the DCA to the vertex:
      track->GetImpactParametersTPC(dDCA,cDCA);
      // set the constrained parameters to the track
      track->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),exParam.GetCovariance());
    }
    
    track->GetPxPyPz(p);

    Float_t pT = track->Pt();
    if(pT<ptMin||pT>ptMax){
      delete track;
      continue;
    }

    // 


    track->GetXYZ(pos);
    track->GetCovarianceXYZPxPyPz(covTr);
    track->GetESDpid(pid);

    if(fMChandler)fMChandler->SelectParticle(esdTrack->GetLabel());
    aodTrack = new(Tracks()[fNumberOfTracks++]) AliAODTrack((track->GetID()+1)*-1,
                                                            track->GetLabel(),
                                                            p,
                                                            kTRUE,
                                                            pos,
                                                            kFALSE,
                                                            covTr, 
                                                            (Short_t)track->GetSign(),
                                                            track->GetITSClusterMap(), 
                                                            pid,
                                                            fPrimaryVertex,
                                                            kTRUE, // check if this is right
                                                            vtx->UsesTrack(track->GetID()),
                                                            AliAODTrack::kPrimary, 
                                                            selectInfo);
    aodTrack->SetIsHybridTPCConstrainedGlobal(isHybridITSTPC);    
    aodTrack->SetTPCClusterMap(track->GetTPCClusterMap());
    aodTrack->SetTPCSharedMap (track->GetTPCSharedMap());
    aodTrack->SetIsTPCConstrained(kTRUE);    
    Float_t ndf = track->GetTPCNcls()+1 - 5 ;
    if(ndf>0){
      aodTrack->SetChi2perNDF(track->GetConstrainedChi2TPC());
    }
    else{
      aodTrack->SetChi2perNDF(-1);
    }

    // set the DCA values to the AOD track
    aodTrack->SetPxPyPzAtDCA(pDCA[0],pDCA[1],pDCA[2]);
    aodTrack->SetXYAtDCA(rDCA[0],rDCA[1]);
    aodTrack->SetDCA(dDCA[0],dDCA[1]);

    aodTrack->SetFlags(track->GetStatus());
    aodTrack->SetTPCPointsF(track->GetTPCNclsF());

    delete track;
  } // end of loop on tracks
  
}


void AliAnalysisTaskESDfilter::ConvertGlobalConstrainedTracks(const AliESDEvent& esd)
{

  // Here we have the option to store the complement from global constraint information
  // to tracks passing tight cuts (1) in order not to get fakes back in, one needs 
  // two sets of cuts one tight (1) (to throw out fakes) and one lose (2) (fakes/bad tracks would pass (2) but not (1))
  // using cut number (3) selects the tracks that complement (1) e.g. tracks witout ITS refit or cluster requirement


  AliCodeTimerAuto("",0);
  
  // Loop over the tracks and extract and mask out all aod tracks that pass the selections for AODt racks
  for(int it = 0;it < fNumberOfTracks;++it)
  {
    AliAODTrack *tr = (AliAODTrack*)(Tracks().At(it));
    if(!tr)continue;
    UInt_t map = tr->GetFilterMap();
    if(map&fGlobalConstrainedFilterMask){
      // we only reset the track select info, no deletion...
      // mask reset mask in case track is already taken
      tr->SetFilterMap(map&~fGlobalConstrainedFilterMask);
    }
    if(map&fHybridFilterMaskGCG){
      // this is one part of the hybrid tracks
      // the others not passing the selection will be the ones selected below
      tr->SetIsHybridGlobalConstrainedGlobal(kTRUE);
    }
  }
  // Loop over the ESD trcks and pick out the tracks passing the GlobalConstraint cuts
 

  Double_t pos[3] = { 0. };      
  Double_t covTr[21]={0.};
  Double_t pid[10]={0.};  
  Double_t p[3] = { 0. };

  Double_t pDCA[3] = { 0. }; // momentum at DCA
  Double_t rDCA[3] = { 0. }; // position at DCA
  Float_t  dDCA[2] = {0.};    // DCA to the vertex d and z
  Float_t  cDCA[3] = {0.};    // covariance of impact parameters


  AliAODTrack* aodTrack(0x0);
  const AliESDVertex *vtx = esd.GetPrimaryVertex();

  // account for change in pT after the constraint
  Float_t ptMax = 1E10;
  Float_t ptMin = 0;
  for(int i = 0;i<32;i++){
    if(fGlobalConstrainedFilterMask&(1<<i)){
      AliESDtrackCuts*cuts = (AliESDtrackCuts*)fTrackFilter->GetCuts()->At(i);
      Float_t tmp1= 0,tmp2 = 0;
      cuts->GetPtRange(tmp1,tmp2);
      if(tmp1>ptMin)ptMin=tmp1;
      if(tmp2<ptMax)ptMax=tmp2;
    }
  } 



 for (Int_t nTrack = 0; nTrack < esd.GetNumberOfTracks(); ++nTrack) 
  {
    AliESDtrack* esdTrack = esd.GetTrack(nTrack); //carefull do not modify it othwise  need to work with a copy 
    const AliExternalTrackParam * exParamGC = esdTrack->GetConstrainedParam();
    if(!exParamGC)continue;

    UInt_t selectInfo = 0;
    Bool_t isHybridGC = false;

    //
    // Track selection
    if (fTrackFilter) {
      selectInfo = fTrackFilter->IsSelected(esdTrack);
    }


    if(!(selectInfo&fHybridFilterMaskGCG))isHybridGC = true;
    if (fWriteHybridGCOnly&&!isHybridGC)continue; // write only complementary tracks

    selectInfo &= fGlobalConstrainedFilterMask;
    if (!selectInfo)continue;
    // fetch the track parameters at the DCA (unconstrained)
    esdTrack->GetPxPyPz(pDCA);
    esdTrack->GetXYZ(rDCA);
    // get the DCA to the vertex:
    esdTrack->GetImpactParameters(dDCA,cDCA);

    esdTrack->GetConstrainedPxPyPz(p);


    Float_t pT = exParamGC->Pt();
    if(pT<ptMin||pT>ptMax){
      continue;
    }


    esdTrack->GetConstrainedXYZ(pos);
    exParamGC->GetCovarianceXYZPxPyPz(covTr);
    esdTrack->GetESDpid(pid);
    if(fMChandler)fMChandler->SelectParticle(esdTrack->GetLabel());
    aodTrack = new(Tracks()[fNumberOfTracks++]) AliAODTrack((esdTrack->GetID()+1)*-1,
                                                            esdTrack->GetLabel(),
                                                            p,
                                                            kTRUE,
                                                            pos,
                                                            kFALSE,
                                                            covTr, 
                                                            (Short_t)esdTrack->GetSign(),
                                                            esdTrack->GetITSClusterMap(), 
                                                            pid,
                                                            fPrimaryVertex,
                                                            kTRUE, // check if this is right
                                                            vtx->UsesTrack(esdTrack->GetID()),
                                                            AliAODTrack::kPrimary, 
                                                            selectInfo);
    aodTrack->SetIsHybridGlobalConstrainedGlobal(isHybridGC);    
    aodTrack->SetIsGlobalConstrained(kTRUE);    
    aodTrack->SetTPCClusterMap(esdTrack->GetTPCClusterMap());
    aodTrack->SetTPCSharedMap (esdTrack->GetTPCSharedMap());
    Float_t ndf = esdTrack->GetTPCNcls()+1 - 5 ;
    if(ndf>0){
      aodTrack->SetChi2perNDF(esdTrack->GetConstrainedChi2TPC());
    }
    else{
      aodTrack->SetChi2perNDF(-1);
    }

    // set the DCA values to the AOD track
    aodTrack->SetPxPyPzAtDCA(pDCA[0],pDCA[1],pDCA[2]);
    aodTrack->SetXYAtDCA(rDCA[0],rDCA[1]);
    aodTrack->SetDCA(dDCA[0],dDCA[1]);

    aodTrack->SetFlags(esdTrack->GetStatus());
    aodTrack->SetTPCPointsF(esdTrack->GetTPCNclsF());
  } // end of loop on tracks
  
}


//______________________________________________________________________________
void AliAnalysisTaskESDfilter::ConvertTracks(const AliESDEvent& esd)
{
  // Tracks (primary and orphan)

  AliCodeTimerAuto("",0);
  
  AliDebug(1,Form("NUMBER OF ESD TRACKS %5d\n", esd.GetNumberOfTracks()));
  
  const AliESDVertex *vtx = esd.GetPrimaryVertex();
  Double_t p[3] = { 0. };
  Double_t pos[3] = { 0. };
  Double_t covTr[21] = { 0. };
  Double_t pid[10] = { 0. };
  AliAODTrack* aodTrack(0x0);
  AliAODPid* detpid(0x0);
  
  for (Int_t nTrack = 0; nTrack < esd.GetNumberOfTracks(); ++nTrack) 
  {
    if (fUsedTrack[nTrack]) continue;
    
    AliESDtrack *esdTrack = esd.GetTrack(nTrack);
    UInt_t selectInfo = 0;
    //
    // Track selection
    if (fTrackFilter) {
	    selectInfo = fTrackFilter->IsSelected(esdTrack);
	    if (!selectInfo && !vtx->UsesTrack(esdTrack->GetID())) continue;
    }
    
    
    esdTrack->GetPxPyPz(p);
    esdTrack->GetXYZ(pos);
    esdTrack->GetCovarianceXYZPxPyPz(covTr);
    esdTrack->GetESDpid(pid);
    if(fMChandler)fMChandler->SelectParticle(esdTrack->GetLabel());
    fPrimaryVertex->AddDaughter(aodTrack =
                         new(Tracks()[fNumberOfTracks++]) AliAODTrack(esdTrack->GetID(),
                                                            esdTrack->GetLabel(),
                                                            p,
                                                            kTRUE,
                                                            pos,
                                                            kFALSE,
                                                            covTr, 
                                                            (Short_t)esdTrack->GetSign(),
                                                            esdTrack->GetITSClusterMap(), 
                                                            pid,
                                                            fPrimaryVertex,
                                                            kTRUE, // check if this is right
                                                            vtx->UsesTrack(esdTrack->GetID()),
                                                            AliAODTrack::kPrimary, 
                                                            selectInfo)
                         );
    aodTrack->SetTPCClusterMap(esdTrack->GetTPCClusterMap());
    aodTrack->SetTPCSharedMap (esdTrack->GetTPCSharedMap());
    aodTrack->SetChi2perNDF(Chi2perNDF(esdTrack));
    aodTrack->SetTPCPointsF(esdTrack->GetTPCNclsF());
    if(esdTrack->IsEMCAL()) aodTrack->SetEMCALcluster(esdTrack->GetEMCALcluster());
    if(esdTrack->IsPHOS())  aodTrack->SetPHOScluster(esdTrack->GetPHOScluster());

    fAODTrackRefs->AddAt(aodTrack, nTrack);
    
    
    if (esdTrack->GetSign() > 0) ++fNumberOfPositiveTracks;
    aodTrack->SetFlags(esdTrack->GetStatus());
    aodTrack->ConvertAliPIDtoAODPID();
    SetAODPID(esdTrack,aodTrack,detpid);
  } // end of loop on tracks
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::ConvertPmdClusters(const AliESDEvent& esd)
{
  AliCodeTimerAuto("",0);
  Int_t jPmdClusters=0;
  // Access to the AOD container of PMD clusters
  TClonesArray &pmdClusters = *(AODEvent()->GetPmdClusters());
  for (Int_t iPmd = 0; iPmd < esd.GetNumberOfPmdTracks(); ++iPmd) {
    // file pmd clusters, to be revised!
    AliESDPmdTrack *pmdTrack = esd.GetPmdTrack(iPmd);
    Int_t nLabel = 0;
    Int_t *label = 0x0;
    Double_t posPmd[3] = { pmdTrack->GetClusterX(), pmdTrack->GetClusterY(), pmdTrack->GetClusterZ()};
    Double_t pidPmd[13] = { 0.}; // to be revised!
    // type not set!
    // assoc cluster not set
    new(pmdClusters[jPmdClusters++]) AliAODPmdCluster(iPmd, nLabel, label, pmdTrack->GetClusterADC(), posPmd, pidPmd);
  }
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::ConvertCaloClusters(const AliESDEvent& esd)
{
  AliCodeTimerAuto("",0);
  
  // Access to the AOD container of clusters
  TClonesArray &caloClusters = *(AODEvent()->GetCaloClusters());
  Int_t jClusters(0);
  
  for (Int_t iClust=0; iClust<esd.GetNumberOfCaloClusters(); ++iClust) {
    
    AliESDCaloCluster * cluster = esd.GetCaloCluster(iClust);
    
    Int_t  id        = cluster->GetID();
    Int_t  nLabel    = cluster->GetNLabels();
    Int_t *labels    = cluster->GetLabels();
    if(labels){ 
		  for(int i = 0;i < nLabel;++i){
			  if(fMChandler)fMChandler->SelectParticle(labels[i]);
		  }
	  }		
    
    Float_t energy = cluster->E();
    Float_t posF[3] = { 0.};
    cluster->GetPosition(posF);
    
    AliAODCaloCluster *caloCluster = new(caloClusters[jClusters++]) AliAODCaloCluster(id,
                                                                                      nLabel,
                                                                                      labels,
                                                                                      energy,
                                                                                      posF,
                                                                                      NULL,
                                                                                      cluster->GetType(),0);
    
    caloCluster->SetCaloCluster(cluster->GetDistanceToBadChannel(),
                                cluster->GetDispersion(),
                                cluster->GetM20(), cluster->GetM02(),
                                cluster->GetEmcCpvDistance(),  
                                cluster->GetNExMax(),cluster->GetTOF()) ;
    
    caloCluster->SetPIDFromESD(cluster->GetPID());
    caloCluster->SetNCells(cluster->GetNCells());
    caloCluster->SetCellsAbsId(cluster->GetCellsAbsId());
    caloCluster->SetCellsAmplitudeFraction(cluster->GetCellsAmplitudeFraction());
    
    TArrayI* matchedT = 	cluster->GetTracksMatched();
    if (fNumberOfTracks>0 && matchedT && cluster->GetTrackMatchedIndex() >= 0) {	
      for (Int_t im = 0; im < matchedT->GetSize(); im++) {
        Int_t iESDtrack = matchedT->At(im);;
        if (fAODTrackRefs->At(iESDtrack) != 0) {
          caloCluster->AddTrackMatched((AliAODTrack*)fAODTrackRefs->At(iESDtrack));
        }
      }
    }
    
  } 
  caloClusters.Expand(jClusters); // resize TObjArray to 'remove' slots for pseudo clusters	 
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::ConvertEMCALCells(const AliESDEvent& esd)
{
  AliCodeTimerAuto("",0);
  // fill EMCAL cell info
  if (esd.GetEMCALCells()) { // protection against missing ESD information
    AliESDCaloCells &esdEMcells = *(esd.GetEMCALCells());
    Int_t nEMcell = esdEMcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
    aodEMcells.CreateContainer(nEMcell);
    aodEMcells.SetType(AliAODCaloCells::kEMCALCell);
    for (Int_t iCell = 0; iCell < nEMcell; iCell++) {      
      aodEMcells.SetCell(iCell,esdEMcells.GetCellNumber(iCell),esdEMcells.GetAmplitude(iCell));
    }
    aodEMcells.Sort();
  }
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::ConvertPHOSCells(const AliESDEvent& esd)
{
  AliCodeTimerAuto("",0);
  // fill PHOS cell info
  if (esd.GetPHOSCells()) { // protection against missing ESD information
    AliESDCaloCells &esdPHcells = *(esd.GetPHOSCells());
    Int_t nPHcell = esdPHcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodPHcells = *(AODEvent()->GetPHOSCells());
    aodPHcells.CreateContainer(nPHcell);
    aodPHcells.SetType(AliAODCaloCells::kPHOSCell);
    for (Int_t iCell = 0; iCell < nPHcell; iCell++) {      
      aodPHcells.SetCell(iCell,esdPHcells.GetCellNumber(iCell),esdPHcells.GetAmplitude(iCell));
    }
    aodPHcells.Sort();
  }
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::ConvertTracklets(const AliESDEvent& esd)
{
  // tracklets    
  AliCodeTimerAuto("",0);

  AliAODTracklets &SPDTracklets = *(AODEvent()->GetTracklets());
  const AliMultiplicity *mult = esd.GetMultiplicity();
  if (mult) {
    if (mult->GetNumberOfTracklets()>0) {
      SPDTracklets.CreateContainer(mult->GetNumberOfTracklets());
      
      for (Int_t n=0; n<mult->GetNumberOfTracklets(); n++) {
        if(fMChandler){
          fMChandler->SelectParticle(mult->GetLabel(n, 0));
          fMChandler->SelectParticle(mult->GetLabel(n, 1));
        }
        SPDTracklets.SetTracklet(n, mult->GetTheta(n), mult->GetPhi(n), mult->GetDeltaPhi(n), mult->GetLabel(n, 0),mult->GetLabel(n, 1));
      }
    }
  } else {
    //Printf("ERROR: AliMultiplicity could not be retrieved from ESD");
  }
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::ConvertKinks(const AliESDEvent& esd)
{
  AliCodeTimerAuto("",0);
  
  // Kinks: it is a big mess the access to the information in the kinks
  // The loop is on the tracks in order to find the mother and daugther of each kink
  
  Double_t covTr[21]={0.};
  Double_t pid[10]={0.};
  AliAODPid* detpid(0x0);
  
  fNumberOfKinks = esd.GetNumberOfKinks();

  const AliESDVertex* vtx = esd.GetPrimaryVertex();
  
  for (Int_t iTrack=0; iTrack<esd.GetNumberOfTracks(); ++iTrack) 
  {
    AliESDtrack * esdTrack = esd.GetTrack(iTrack);
    
    Int_t ikink = esdTrack->GetKinkIndex(0);
    
    if (ikink && fNumberOfKinks) {
	    // Negative kink index: mother, positive: daughter
	    
	    // Search for the second track of the kink
	    
	    for (Int_t jTrack = iTrack+1; jTrack<esd.GetNumberOfTracks(); ++jTrack) {
        
        AliESDtrack * esdTrack1 = esd.GetTrack(jTrack);
        
        Int_t jkink = esdTrack1->GetKinkIndex(0);
        
        if ( TMath::Abs(ikink)==TMath::Abs(jkink) ) {
          
          // The two tracks are from the same kink
          
          if (fUsedKink[TMath::Abs(ikink)-1]) continue; // skip used kinks
          
          Int_t imother = -1;
          Int_t idaughter = -1;
          
          if (ikink<0 && jkink>0) {
            
            imother = iTrack;
            idaughter = jTrack;
          }
          else if (ikink>0 && jkink<0) {
            
            imother = jTrack;
            idaughter = iTrack;
          }
          else {
            //			cerr << "Error: Wrong combination of kink indexes: "
            //			     << ikink << " " << jkink << endl;
            continue;
          }
          
          // Add the mother track if it passed primary track selection cuts
          
          AliAODTrack * mother = NULL;
          
          UInt_t selectInfo = 0;
          if (fTrackFilter) {
            selectInfo = fTrackFilter->IsSelected(esd.GetTrack(imother));
            if (!selectInfo) continue;
          }
          
          if (!fUsedTrack[imother]) {
            
            fUsedTrack[imother] = kTRUE;
            
            AliESDtrack *esdTrackM = esd.GetTrack(imother);
            Double_t p[3] = { 0. };
            Double_t pos[3] = { 0. };
            esdTrackM->GetPxPyPz(p);
            esdTrackM->GetXYZ(pos);
            esdTrackM->GetCovarianceXYZPxPyPz(covTr);
            esdTrackM->GetESDpid(pid);
            if(fMChandler)fMChandler->SelectParticle(esdTrackM->GetLabel());
            mother = 
            new(Tracks()[fNumberOfTracks++]) AliAODTrack(esdTrackM->GetID(),
                                               esdTrackM->GetLabel(),
                                               p,
                                               kTRUE,
                                               pos,
                                               kFALSE,
                                               covTr, 
                                               (Short_t)esdTrackM->GetSign(),
                                               esdTrackM->GetITSClusterMap(), 
                                               pid,
                                               fPrimaryVertex,
                                               kTRUE, // check if this is right
                                               vtx->UsesTrack(esdTrack->GetID()),
                                               AliAODTrack::kPrimary,
                                               selectInfo);
            mother->SetTPCClusterMap(esdTrackM->GetTPCClusterMap());
            mother->SetTPCSharedMap (esdTrackM->GetTPCSharedMap());
            mother->SetChi2perNDF(Chi2perNDF(esdTrackM));
            mother->SetTPCPointsF(esdTrackM->GetTPCNclsF());

            fAODTrackRefs->AddAt(mother, imother);
            
            if (esdTrackM->GetSign() > 0) ++fNumberOfPositiveTracks;
            mother->SetFlags(esdTrackM->GetStatus());
            mother->ConvertAliPIDtoAODPID();
            fPrimaryVertex->AddDaughter(mother);
            mother->ConvertAliPIDtoAODPID();
            SetAODPID(esdTrackM,mother,detpid);
          }
          else {
            //			cerr << "Error: event " << esd.GetEventNumberInFile() << " kink " << TMath::Abs(ikink)-1
            //			     << " track " << imother << " has already been used!" << endl;
          }
          
          // Add the kink vertex
          AliESDkink * kink = esd.GetKink(TMath::Abs(ikink)-1);
          
          AliAODVertex * vkink = 
          new(Vertices()[fNumberOfVertices++]) AliAODVertex(kink->GetPosition(),
                                                  NULL,
                                                  0.,
                                                  mother,
                                                  esdTrack->GetID(),  // This is the track ID of the mother's track!
                                                  AliAODVertex::kKink);
          // Add the daughter track
          
          AliAODTrack * daughter = NULL;
          
          if (!fUsedTrack[idaughter]) {
            
            fUsedTrack[idaughter] = kTRUE;
            
            AliESDtrack *esdTrackD = esd.GetTrack(idaughter);
            Double_t p[3] = { 0. };
            Double_t pos[3] = { 0. };

            esdTrackD->GetPxPyPz(p);
            esdTrackD->GetXYZ(pos);
            esdTrackD->GetCovarianceXYZPxPyPz(covTr);
            esdTrackD->GetESDpid(pid);
            selectInfo = 0;
            if (fTrackFilter) selectInfo = fTrackFilter->IsSelected(esdTrackD);
            if(fMChandler)fMChandler->SelectParticle(esdTrackD->GetLabel());
            daughter = 
            new(Tracks()[fNumberOfTracks++]) AliAODTrack(esdTrackD->GetID(),
                                               esdTrackD->GetLabel(),
                                               p,
                                               kTRUE,
                                               pos,
                                               kFALSE,
                                               covTr, 
                                               (Short_t)esdTrackD->GetSign(),
                                               esdTrackD->GetITSClusterMap(), 
                                               pid,
                                               vkink,
                                               kTRUE, // check if this is right
                                               vtx->UsesTrack(esdTrack->GetID()),
                                               AliAODTrack::kSecondary,
                                               selectInfo);
            daughter->SetTPCClusterMap(esdTrackD->GetTPCClusterMap());
            daughter->SetTPCSharedMap (esdTrackD->GetTPCSharedMap());
	    daughter->SetTPCPointsF(esdTrackD->GetTPCNclsF());
            fAODTrackRefs->AddAt(daughter, idaughter);
            
            if (esdTrackD->GetSign() > 0) ++fNumberOfPositiveTracks;
            daughter->SetFlags(esdTrackD->GetStatus());
            daughter->ConvertAliPIDtoAODPID();
            vkink->AddDaughter(daughter);
            daughter->ConvertAliPIDtoAODPID();
            SetAODPID(esdTrackD,daughter,detpid);
          }
          else {
            //			cerr << "Error: event " << esd.GetEventNumberInFile() << " kink " << TMath::Abs(ikink)-1
            //			     << " track " << idaughter << " has already been used!" << endl;
          }
        }
	    }
    }      
  }
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::ConvertPrimaryVertices(const AliESDEvent& esd)
{
  AliCodeTimerAuto("",0);
  
  // Access to the AOD container of vertices
  fNumberOfVertices = 0;
  
  Double_t pos[3] = { 0. };
  Double_t covVtx[6] = { 0. };

  // Add primary vertex. The primary tracks will be defined
  // after the loops on the composite objects (V0, cascades, kinks)
  const AliESDVertex *vtx = esd.GetPrimaryVertex();
  
  vtx->GetXYZ(pos); // position
  vtx->GetCovMatrix(covVtx); //covariance matrix
  
  fPrimaryVertex = new(Vertices()[fNumberOfVertices++])
  AliAODVertex(pos, covVtx, vtx->GetChi2toNDF(), NULL, -1, AliAODVertex::kPrimary);
  fPrimaryVertex->SetName(vtx->GetName());
  fPrimaryVertex->SetTitle(vtx->GetTitle());
  
  TString vtitle = vtx->GetTitle();
  if (!vtitle.Contains("VertexerTracks")) 
    fPrimaryVertex->SetNContributors(vtx->GetNContributors());
  
  if (fDebug > 0) fPrimaryVertex->Print();  
  
  // Add SPD "main" vertex 
  const AliESDVertex *vtxS = esd.GetPrimaryVertexSPD();
  vtxS->GetXYZ(pos); // position
  vtxS->GetCovMatrix(covVtx); //covariance matrix
  AliAODVertex * mVSPD = new(Vertices()[fNumberOfVertices++])
  AliAODVertex(pos, covVtx, vtxS->GetChi2toNDF(), NULL, -1, AliAODVertex::kMainSPD);
  mVSPD->SetName(vtxS->GetName());
  mVSPD->SetTitle(vtxS->GetTitle());
  mVSPD->SetNContributors(vtxS->GetNContributors()); 
  
  // Add SPD pileup vertices
  for(Int_t iV=0; iV<esd.GetNumberOfPileupVerticesSPD(); ++iV)
  {
    const AliESDVertex *vtxP = esd.GetPileupVertexSPD(iV);
    vtxP->GetXYZ(pos); // position
    vtxP->GetCovMatrix(covVtx); //covariance matrix
    AliAODVertex * pVSPD =  new(Vertices()[fNumberOfVertices++])
    AliAODVertex(pos, covVtx, vtxP->GetChi2toNDF(), NULL, -1, AliAODVertex::kPileupSPD);
    pVSPD->SetName(vtxP->GetName());
    pVSPD->SetTitle(vtxP->GetTitle());
    pVSPD->SetNContributors(vtxP->GetNContributors()); 
    pVSPD->SetBC(vtxP->GetBC());
  }
  
  // Add TRK pileup vertices
  for(Int_t iV=0; iV<esd.GetNumberOfPileupVerticesTracks(); ++iV)
  {
    const AliESDVertex *vtxP = esd.GetPileupVertexTracks(iV);
    vtxP->GetXYZ(pos); // position
    vtxP->GetCovMatrix(covVtx); //covariance matrix
    AliAODVertex * pVTRK = new(Vertices()[fNumberOfVertices++])
    AliAODVertex(pos, covVtx, vtxP->GetChi2toNDF(), NULL, -1, AliAODVertex::kPileupTracks);
    pVTRK->SetName(vtxP->GetName());
    pVTRK->SetTitle(vtxP->GetTitle());
    pVTRK->SetNContributors(vtxP->GetNContributors());
    pVTRK->SetBC(vtxP->GetBC());
  }
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::ConvertVZERO(const AliESDEvent& esd)
{
  // Convert VZERO data
  AliAODVZERO* vzeroData = AODEvent()->GetVZEROData();
  *vzeroData = *(esd.GetVZEROData());
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::ConvertZDC(const AliESDEvent& esd)
{
  // Convert ZDC data
  AliESDZDC* esdZDC = esd.GetZDCData();
  
  const Double_t zem1Energy = esdZDC->GetZEM1Energy();
  const Double_t zem2Energy = esdZDC->GetZEM2Energy();
   
  const Double_t *towZNC = esdZDC->GetZNCTowerEnergy();
  const Double_t *towZPC = esdZDC->GetZPCTowerEnergy();
  const Double_t *towZNA = esdZDC->GetZNATowerEnergy();
  const Double_t *towZPA = esdZDC->GetZPATowerEnergy();
  const Double_t *towZNCLG = esdZDC->GetZNCTowerEnergyLR();
  const Double_t *towZNALG = esdZDC->GetZNATowerEnergyLR();
  
  AliAODZDC* zdcAOD = AODEvent()->GetZDCData();

  zdcAOD->SetZEM1Energy(zem1Energy);
  zdcAOD->SetZEM2Energy(zem2Energy);
  zdcAOD->SetZNCTowers(towZNC, towZNCLG);
  zdcAOD->SetZNATowers(towZNA, towZNALG);
  zdcAOD->SetZPCTowers(towZPC);
  zdcAOD->SetZPATowers(towZPA);
  
  zdcAOD->SetZDCParticipants(esdZDC->GetZDCParticipants(), esdZDC->GetZDCPartSideA(), esdZDC->GetZDCPartSideC());
  zdcAOD->SetZDCImpactParameter(esdZDC->GetImpactParameter(), esdZDC->GetImpactParamSideA(), 
  	esdZDC->GetImpactParamSideC());
  zdcAOD->SetZDCTDCSum(esdZDC->GetZNTDCSum(0));	
  zdcAOD->SetZDCTDCDiff(esdZDC->GetZNTDCDiff(0));	

}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::ConvertESDtoAOD() 
{
  // ESD Filter analysis task executed for each event
  
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  
  if(!esd)return;

  AliCodeTimerAuto("",0);
  
  fOldESDformat = ( esd->GetAliESDOld() != 0x0 );
  
  fNumberOfTracks = 0;
  fNumberOfPositiveTracks = 0;
  fNumberOfV0s = 0;
  fNumberOfVertices = 0;
  fNumberOfCascades = 0;
  fNumberOfKinks = 0;
    
  AliAODHeader* header = ConvertHeader(*esd);

  if ( fIsVZEROEnabled ) ConvertVZERO(*esd);
  
  // Fetch Stack for debuggging if available 
  fMChandler=0x0;
  if(MCEvent())
  {
    fMChandler = (AliMCEventHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler()); 
  }
  
  // loop over events and fill them
  // Multiplicity information needed by the header (to be revised!)
  Int_t nTracks    = esd->GetNumberOfTracks();
  for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) esd->GetTrack(iTrack)->SetESDEvent(esd);

  // Update the header

  Int_t nV0s      = esd->GetNumberOfV0s();
  Int_t nCascades = esd->GetNumberOfCascades();
  Int_t nKinks    = esd->GetNumberOfKinks();
  Int_t nVertices = nV0s + nCascades /*V0 wihtin cascade already counted*/+ nKinks + 1 /* = prim. vtx*/;
  Int_t nPileSPDVertices=1+esd->GetNumberOfPileupVerticesSPD(); // also SPD main vertex
  Int_t nPileTrkVertices=esd->GetNumberOfPileupVerticesTracks();
  nVertices+=nPileSPDVertices;
  nVertices+=nPileTrkVertices;
  Int_t nJets     = 0;
  Int_t nCaloClus = esd->GetNumberOfCaloClusters();
  Int_t nFmdClus  = 0;
  Int_t nPmdClus  = esd->GetNumberOfPmdTracks();
    
  AliDebug(1,Form("   NV0=%d  NCASCADES=%d  NKINKS=%d", nV0s, nCascades, nKinks));
       
  AODEvent()->ResetStd(nTracks, nVertices, nV0s, nCascades, nJets, nCaloClus, nFmdClus, nPmdClus);

  if (nV0s > 0) 
  {
    // RefArray to store a mapping between esd V0 number and newly created AOD-Vertex V0
    fAODV0VtxRefs = new TRefArray(nV0s);
    // RefArray to store the mapping between esd V0 number and newly created AOD-V0
    fAODV0Refs = new TRefArray(nV0s); 
    // Array to take into account the V0s already added to the AOD (V0 within cascades)
    fUsedV0 = new Bool_t[nV0s];
    for (Int_t iV0=0; iV0<nV0s; ++iV0) fUsedV0[iV0]=kFALSE;
  }
  
  if (nTracks>0) 
  {
    // RefArray to store the mapping between esd track number and newly created AOD-Track
    
    fAODTrackRefs = new TRefArray(nTracks);

    // Array to take into account the tracks already added to the AOD    
    fUsedTrack = new Bool_t[nTracks];
    for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) fUsedTrack[iTrack]=kFALSE;
  }
  
  // Array to take into account the kinks already added to the AOD
  if (nKinks>0) 
  {
    fUsedKink = new Bool_t[nKinks];
    for (Int_t iKink=0; iKink<nKinks; ++iKink) fUsedKink[iKink]=kFALSE;
  }
    
  ConvertPrimaryVertices(*esd);

  //setting best TOF PID
  AliESDInputHandler* esdH = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (esdH)
      fESDpid = esdH->GetESDpid();

  if (fIsPidOwner && fESDpid){
    delete fESDpid;
    fESDpid = 0;
  }
  if(!fESDpid)
  { //in case of no Tender attached 
    fESDpid = new AliESDpid;
    fIsPidOwner = kTRUE;
  }
  
  if(!esd->GetTOFHeader())
  { //protection in case the pass2 LHC10b,c,d have been processed without tender. 
    Float_t t0spread[10];
    Float_t intrinsicTOFres=100; //ps ok for LHC10b,c,d pass2!! 
    for (Int_t i=0; i<10; i++) t0spread[i] = (TMath::Sqrt(esd->GetSigma2DiamondZ()))/0.03; //0.03 to convert from cm to ps
    fESDpid->GetTOFResponse().SetT0resolution(t0spread);
    fESDpid->GetTOFResponse().SetTimeResolution(intrinsicTOFres);
    
    fESDpid->SetTOFResponse(esd, (AliESDpid::EStartTimeType_t)fTimeZeroType);    
  }
  
  if(esd->GetTOFHeader() && fIsPidOwner) fESDpid->SetTOFResponse(esd, (AliESDpid::EStartTimeType_t)fTimeZeroType); //in case of AOD production strating form LHC10e without Tender. 
  
  if ( fAreCascadesEnabled ) ConvertCascades(*esd);

  if ( fAreV0sEnabled ) ConvertV0s(*esd);
  
  if ( fAreKinksEnabled ) ConvertKinks(*esd);
  
  if ( fAreTracksEnabled ) ConvertTracks(*esd);
  
  // Update number of AOD tracks in header at the end of track loop (M.G.)
  header->SetRefMultiplicity(fNumberOfTracks);
  header->SetRefMultiplicityPos(fNumberOfPositiveTracks);
  header->SetRefMultiplicityNeg(fNumberOfTracks - fNumberOfPositiveTracks);

  if ( fTPCConstrainedFilterMask ) ConvertTPCOnlyTracks(*esd);
  if( fGlobalConstrainedFilterMask) ConvertGlobalConstrainedTracks(*esd);  

  if ( fArePmdClustersEnabled ) ConvertPmdClusters(*esd);
  
  if ( fAreCaloClustersEnabled ) ConvertCaloClusters(*esd);
  
  if ( fAreEMCALCellsEnabled )ConvertEMCALCells(*esd);
  
  if ( fArePHOSCellsEnabled )ConvertPHOSCells(*esd);
  
  if ( fAreTrackletsEnabled ) ConvertTracklets(*esd);
  
  delete fAODTrackRefs; fAODTrackRefs=0x0;
  delete fAODV0VtxRefs; fAODV0VtxRefs=0x0;
  delete fAODV0Refs; fAODV0Refs=0x0;
  
  delete[] fUsedTrack; fUsedTrack=0x0;
  delete[] fUsedV0; fUsedV0=0x0;
  delete[] fUsedKink; fUsedKink=0x0;

  if ( fIsPidOwner){
    delete fESDpid;
    fESDpid = 0x0;
  }


}


//______________________________________________________________________________
void AliAnalysisTaskESDfilter::SetAODPID(AliESDtrack *esdtrack, AliAODTrack *aodtrack, AliAODPid *detpid)
{
  //
  // Setter for the raw PID detector signals
  //

  // Save PID object for candidate electrons
    Bool_t pidSave = kFALSE;
    if (fTrackFilter) {
	Bool_t selectInfo = fTrackFilter->IsSelected((char*) "Electrons");
	if (selectInfo)  pidSave = kTRUE;
    }


    // Tracks passing pt cut 
    if(esdtrack->Pt()>fHighPthreshold) {
	pidSave = kTRUE;
    } else {
	if(fPtshape){
	    if(esdtrack->Pt()> fPtshape->GetXmin()){
		Double_t y = fPtshape->Eval(esdtrack->Pt())/fPtshape->Eval(fHighPthreshold);
		if(gRandom->Rndm(0)<1./y){
		    pidSave = kTRUE;
		}//end rndm
	    }//end if p < pmin
	}//end if p function
    }// end else

    if (pidSave) {
      if(!aodtrack->GetDetPid()){// prevent memory leak when calling SetAODPID twice for the same track
	detpid = new AliAODPid();
	SetDetectorRawSignals(detpid,esdtrack);
	aodtrack->SetDetPID(detpid);
      }
    }
}

//______________________________________________________________________________
void AliAnalysisTaskESDfilter::SetDetectorRawSignals(AliAODPid *aodpid, AliESDtrack *track)
{
//
//assignment of the detector signals (AliXXXesdPID inspired)
//
 if(!track){
 AliInfo("no ESD track found. .....exiting");
 return;
 }
 // TPC momentum
 const AliExternalTrackParam *in=track->GetInnerParam();
 if (in) {
   aodpid->SetTPCmomentum(in->GetP());
 }else{
   aodpid->SetTPCmomentum(-1.);
 }


 aodpid->SetITSsignal(track->GetITSsignal());
 Double_t itsdedx[4]; // dE/dx samples for individual ITS layers
 track->GetITSdEdxSamples(itsdedx);
 aodpid->SetITSdEdxSamples(itsdedx);

 aodpid->SetTPCsignal(track->GetTPCsignal());
 aodpid->SetTPCsignalN(track->GetTPCsignalN());

 //n TRD planes = 6
 Int_t nslices = track->GetNumberOfTRDslices()*6;
 TArrayD trdslices(nslices);
 for(Int_t iSl =0; iSl < track->GetNumberOfTRDslices(); iSl++) {
   for(Int_t iPl =0; iPl<6; iPl++) trdslices[iPl*track->GetNumberOfTRDslices()+iSl] = track->GetTRDslice(iPl,iSl);
 }
 
//TRD momentum
 for(Int_t iPl=0;iPl<6;iPl++){
   Double_t trdmom=track->GetTRDmomentum(iPl);
   aodpid->SetTRDmomentum(iPl,trdmom);
 }

 aodpid->SetTRDsignal(track->GetNumberOfTRDslices()*6,trdslices.GetArray());

 //TRD clusters and tracklets
 aodpid->SetTRDncls(track->GetTRDncls());
 aodpid->SetTRDntrackletsPID(track->GetTRDntrackletsPID());
 
 //TOF PID  
 Double_t times[AliAODPid::kSPECIES]; track->GetIntegratedTimes(times);
 aodpid->SetIntegratedTimes(times);

  Float_t tzeroTrack = fESDpid->GetTOFResponse().GetStartTime(track->P());
  aodpid->SetTOFsignal(track->GetTOFsignal()-tzeroTrack);
  
  Double_t tofRes[5];
  for (Int_t iMass=0; iMass<5; iMass++){
    tofRes[iMass]=(Double_t)fESDpid->GetTOFResponse().GetExpectedSigma(track->P(), times[iMass], AliPID::ParticleMass(iMass));
  }
  aodpid->SetTOFpidResolution(tofRes);

  aodpid->SetHMPIDsignal(track->GetHMPIDsignal());

 //Extrapolate track to EMCAL surface for AOD-level track-cluster matching
 Double_t emcpos[3] = {0.,0.,0.};
 Double_t emcmom[3] = {0.,0.,0.};
 aodpid->SetEMCALPosition(emcpos);
 aodpid->SetEMCALMomentum(emcmom);

 Double_t hmpPid[5] = {0};
 track->GetHMPIDpid(hmpPid);
 aodpid->SetHMPIDprobs(hmpPid);

 AliExternalTrackParam *outerparam = (AliExternalTrackParam*)track->GetOuterParam();
 if(!outerparam) return;

 //To be replaced by call to AliEMCALGeoUtils when the class becomes available
 Bool_t okpos = outerparam->GetXYZ(emcpos);
 Bool_t okmom = outerparam->GetPxPyPz(emcmom);
 if(!(okpos && okmom)) return;

 aodpid->SetEMCALPosition(emcpos);
 aodpid->SetEMCALMomentum(emcmom);

}

Double_t  AliAnalysisTaskESDfilter::Chi2perNDF(AliESDtrack* track)
{
    // Calculate chi2 per ndf for track
    Int_t  nClustersTPC = track->GetTPCNcls();

    if ( nClustersTPC > 5) {
       return (track->GetTPCchi2()/Float_t(nClustersTPC - 5));
    } else {
       return (-1.);
    }
 }


//______________________________________________________________________________
void AliAnalysisTaskESDfilter::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisESDfilter: Terminate() \n");
}

//______________________________________________________________________________
void  AliAnalysisTaskESDfilter::PrintMCInfo(AliStack *pStack,Int_t label){
// Print MC info
  if(!pStack)return;
  label = TMath::Abs(label);
  TParticle *part = pStack->Particle(label);
  Printf("########################");
  Printf("%s:%d %d UniqueID %d PDG %d P %3.3f",(char*)__FILE__,__LINE__,label,part->GetUniqueID(),part->GetPdgCode(),part->P());
  part->Print();
  TParticle* mother = part;
  Int_t imo = part->GetFirstMother();
  Int_t nprim = pStack->GetNprimary();
  //  while((imo >= nprim) && (mother->GetUniqueID() == 4)) {
  while((imo >= nprim)) {
    mother =  pStack->Particle(imo);
    Printf("Mother %s:%d Label %d UniqueID %d PDG %d P %3.3f",(char*)__FILE__,__LINE__,imo,mother->GetUniqueID(),mother->GetPdgCode(),mother->P());
    mother->Print();
    imo =  mother->GetFirstMother();
  }
  Printf("########################");
}
