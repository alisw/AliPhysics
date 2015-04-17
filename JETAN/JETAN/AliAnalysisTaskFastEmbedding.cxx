/**************************************************************************
 *                                                                        *
 * Task for fast embedding                                                *
 * read extra input from AOD                                              *
 *                                                                        *
 **************************************************************************/


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

/* $Id: */

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TRef.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TKey.h>
#include <TGrid.h>


#include "AliAnalysisTaskFastEmbedding.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODJet.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliInputEventHandler.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "AliAODv0.h"
#include "AliPIDResponse.h"
#include "AliAODPid.h"
#include "TVector3.h"

#include "AliLog.h"

ClassImp(AliAnalysisTaskFastEmbedding)

//__________________________________________________________________________
AliAnalysisTaskFastEmbedding::AliAnalysisTaskFastEmbedding()
: AliAnalysisTaskSE()
  ,fESD(0)
  ,fAODout(0)
  ,fAODevent(0)
  ,fAODtree(0)
  ,fAODfile(0)
  ,mcHeader(0)
  ,fFFRadius(0)
  ,fCutPostrackEta(0)
  ,fCutNegtrackEta(0)
  ,fCutEta(0)
  ,fCutV0cosPointAngle(0)
  ,fKinkDaughters(0)
  ,fRequireTPCRefit(0)
  ,fCutArmenteros(0)
  ,fCutV0DecayMin(0)
  ,fCutV0DecayMax(0)
  ,fCutV0totMom(0)
  ,fCutDcaV0Daughters(0)
  ,fCutDcaPosToPrimVertex(0)
  ,fCutDcaNegToPrimVertex(0)
  ,fCutV0RadiusMin(0)
  ,fCutV0RadiusMax(0)
  ,fCutBetheBloch(0)
  ,fK0Type(0) 
  ,fLaType(0) 
  ,fALaType(0)
  ,fListK0s(0)
  ,fListLa(0)
  ,fListALa(0)
  ,fListK0sCone(0)
  ,fListLaCone(0)
  ,fListALaCone(0)
  ,fPIDResponse(0)
  ,rndm(0)
  ,fInputEntries(0)
  ,fAODPathArray(0)
  ,fAODEntriesArray(0)
  ,fAODPath("AliAOD.root")
  ,fAODEntries(-1)
  ,fAODEntriesSum(0)
  ,fAODEntriesMax(0)
  ,fOfflineTrgMask(AliVEvent::kAny)
  ,fMinContribVtx(1)
  ,fVtxZMin(-8.)
  ,fVtxZMax(8.)
  ,fEvtClassMin(0)
  ,fEvtClassMax(4)
  ,fCentMin(0.)
  ,fCentMax(100.)
  ,fNInputTracksMin(0)
  ,fNInputTracksMax(-1)
  ,fTrackBranch("aodExtraTracks")
  ,fMCparticlesBranch("aodExtraMCparticles")
  ,fJetBranch("")
  ,fFileId(-1)
  ,fAODEntry(0)
  ,fCountEvents(-1)
  ,fEmbedMode(0)
  ,fEvtSelecMode(0)
  ,fEvtSelMinJetPt(-1)
  ,fEvtSelMaxJetPt(-1)
  ,fEvtSelMinJetEta(-999.)
  ,fEvtSelMaxJetEta( 999.)
  ,fEvtSelMinJetPhi(0.)
  ,fEvtSelMaxJetPhi(TMath::Pi()*2.)
  ,fEvtSelJetMinLConstPt(0)
  ,fExtraEffPb(1)
  ,fToyMinNbOfTracks(1)
  ,fToyMaxNbOfTracks(1)
  ,fToyMinTrackPt(50.)
  ,fToyMaxTrackPt(50.)
  ,fToyDistributionTrackPt(0.)
  ,fToyMinTrackEta(-.5)
  ,fToyMaxTrackEta(.5)
  ,fToyMinTrackPhi(0.)
  ,fToyMaxTrackPhi(2*TMath::Pi())
  ,fToyFilterMap(0)
  ,fTrackFilterMap(0)
  ,fNPtHard(10)
  ,fPtHard(0)
  ,fPtHardBin(-1)
  ,fAODJets(0x0)
  ,fNevents(0)
  ,fXsection(0)
  ,fAvgTrials(0)
  ,fHistList(0)
  ,fHistEvtSelection(0)
  ,fh1Xsec(0)
  ,fh1Trials(0)
  ,fh1TrialsEvtSel(0)
  ,fh2PtHard(0)
  ,fh2PtHardEvtSel(0)
  ,fh2PtHardTrials(0)
  ,fh1TrackPt(0)
  ,fh2TrackEtaPhi(0)
  ,fh1TrackN(0)
  ,fh1V0Pt(0)
  ,fh2V0EtaPhi(0)
  ,fh1K0Pt(0)
  ,fh2K0EtaPhi(0)
  ,fh2K0sPtJetPtCone(0)
  ,fh1LaPt(0)
  ,fh2LaEtaPhi(0)
  ,fh2LaPtJetPtCone(0)
  ,fh1ALaPt(0)
  ,fh2ALaEtaPhi(0)
  ,fh2ALaPtJetPtCone(0)
  ,fh1JetPt(0)
  ,fh2JetEtaPhi(0)
  ,fh1JetN(0)
  ,fh1MCTrackPt(0)
  ,fh2MCTrackEtaPhi(0)
  ,fh1MCTrackN(0)
  ,fh1AODfile(0)
  ,fh2AODevent(0)
  ,fh1trackPosEta(0)            
  ,fh1trackNegEta(0)          
  ,fh1V0Eta(0)        
  ,fh1CosPointAngle(0)           
  ,fh1DecayLengthV0(0)    
  ,fh2ProperLifetimeK0sVsPtBeforeCut(0)    
  ,fh2ProperLifetimeK0sVsPtAfterCut(0)
  ,fh1V0Radius(0)     
  ,fh1DcaV0Daughters(0)        
  ,fh1DcaPosToPrimVertex(0)   
  ,fh1DcaNegToPrimVertex(0)    
  ,fh2ArmenterosBeforeCuts(0)
  ,fh2ArmenterosAfterCuts(0)
  ,fh2BBLaPos(0)
  ,fh2BBLaNeg(0)
{
  // default constructor
  
}


//__________________________________________________________________________
AliAnalysisTaskFastEmbedding::AliAnalysisTaskFastEmbedding(const char *name)
: AliAnalysisTaskSE(name)
,fESD(0)
,fAODout(0)
,fAODevent(0)
,fAODtree(0)
,fAODfile(0)
,mcHeader(0)
,fFFRadius(0)
,fCutPostrackEta(0)
,fCutNegtrackEta(0)
,fCutEta(0)
,fCutV0cosPointAngle(0)
,fKinkDaughters(0)
,fRequireTPCRefit(0)
,fCutArmenteros(0)
,fCutV0DecayMin(0)
,fCutV0DecayMax(0)
,fCutV0totMom(0)
,fCutDcaV0Daughters(0)
,fCutDcaPosToPrimVertex(0)
,fCutDcaNegToPrimVertex(0)
,fCutV0RadiusMin(0)
,fCutV0RadiusMax(0)
,fCutBetheBloch(0)
,fK0Type(0) 
,fLaType(0) 
,fALaType(0)
,fListK0s(0)
,fListLa(0)
,fListALa(0) 
,fListK0sCone(0)
,fListLaCone(0)
,fListALaCone(0) 
,fPIDResponse(0)
,rndm(0)
,fInputEntries(0)
,fAODPathArray(0)
,fAODEntriesArray(0)
,fAODPath("AliAOD.root")
,fAODEntries(-1)
,fAODEntriesSum(0)
,fAODEntriesMax(0)
,fOfflineTrgMask(AliVEvent::kAny)
,fMinContribVtx(1)
,fVtxZMin(-8.)
,fVtxZMax(8.)
,fEvtClassMin(0)
,fEvtClassMax(4)
,fCentMin(0.)
,fCentMax(100.)
,fNInputTracksMin(0)
,fNInputTracksMax(-1)
,fTrackBranch("aodExtraTracks")
,fMCparticlesBranch("aodExtraMCparticles")
,fJetBranch("")
,fFileId(-1)
,fAODEntry(0)
,fCountEvents(-1)
,fEmbedMode(0)
,fEvtSelecMode(0)
,fEvtSelMinJetPt(-1)
,fEvtSelMaxJetPt(-1)
,fEvtSelMinJetEta(-999.)
,fEvtSelMaxJetEta( 999.)
,fEvtSelMinJetPhi(0.)
,fEvtSelMaxJetPhi(TMath::Pi()*2.)
,fEvtSelJetMinLConstPt(0)
,fExtraEffPb(1)
,fToyMinNbOfTracks(1)
,fToyMaxNbOfTracks(1)
,fToyMinTrackPt(50.)
,fToyMaxTrackPt(50.)
,fToyDistributionTrackPt(0.)
,fToyMinTrackEta(-.5)
,fToyMaxTrackEta(.5)
,fToyMinTrackPhi(0.)
,fToyMaxTrackPhi(2*TMath::Pi())
,fToyFilterMap(0)
,fTrackFilterMap(0)
,fNPtHard(10)
,fPtHard(0)
,fPtHardBin(-1)
,fAODJets(0x0)
,fNevents(0)
,fXsection(0)
,fAvgTrials(0)
,fHistList(0)
,fHistEvtSelection(0)
,fh1Xsec(0)
,fh1Trials(0)
,fh1TrialsEvtSel(0)
,fh2PtHard(0)
,fh2PtHardEvtSel(0)
,fh2PtHardTrials(0)
,fh1TrackPt(0)
,fh2TrackEtaPhi(0)
,fh1TrackN(0)
,fh1V0Pt(0)
,fh2V0EtaPhi(0)
,fh1K0Pt(0)
,fh2K0EtaPhi(0)
,fh2K0sPtJetPtCone(0)
,fh1LaPt(0)
,fh2LaEtaPhi(0)
,fh2LaPtJetPtCone(0)
,fh1ALaPt(0)
,fh2ALaEtaPhi(0)
,fh2ALaPtJetPtCone(0)
,fh1JetPt(0)
,fh2JetEtaPhi(0)
,fh1JetN(0)
,fh1MCTrackPt(0)
,fh2MCTrackEtaPhi(0)
,fh1MCTrackN(0)
,fh1AODfile(0)
,fh2AODevent(0)
,fh1trackPosEta(0)            
,fh1trackNegEta(0)          
,fh1V0Eta(0)           
,fh1CosPointAngle(0)           
,fh1DecayLengthV0(0)    
,fh2ProperLifetimeK0sVsPtBeforeCut(0)    
,fh2ProperLifetimeK0sVsPtAfterCut(0)
,fh1V0Radius(0)     
,fh1DcaV0Daughters(0)        
,fh1DcaPosToPrimVertex(0)   
,fh1DcaNegToPrimVertex(0)    
,fh2ArmenterosBeforeCuts(0)
,fh2ArmenterosAfterCuts(0)
,fh2BBLaPos(0)
,fh2BBLaNeg(0)
{
   // constructor
   DefineOutput(1, TList::Class());
}


//__________________________________________________________________________
AliAnalysisTaskFastEmbedding::AliAnalysisTaskFastEmbedding(const AliAnalysisTaskFastEmbedding &copy)
: AliAnalysisTaskSE()
,fESD(copy.fESD)
,fAODout(copy.fAODout)
,fAODevent(copy.fAODevent)
,fAODtree(copy.fAODtree)
,fAODfile(copy.fAODfile)
,mcHeader(copy.mcHeader)
,fFFRadius(copy.fFFRadius)
,fCutPostrackEta(copy.fCutPostrackEta)
,fCutNegtrackEta(copy.fCutNegtrackEta)
,fCutEta(copy.fCutEta)
,fCutV0cosPointAngle(copy.fCutV0cosPointAngle)
,fKinkDaughters(copy.fKinkDaughters)
,fRequireTPCRefit(copy.fRequireTPCRefit)
,fCutArmenteros(copy.fCutArmenteros)
,fCutV0DecayMin(copy.fCutV0DecayMin)
,fCutV0DecayMax(copy.fCutV0DecayMax)
,fCutV0totMom(copy.fCutV0totMom)
,fCutDcaV0Daughters(copy.fCutDcaV0Daughters)
,fCutDcaPosToPrimVertex(copy.fCutDcaPosToPrimVertex)
,fCutDcaNegToPrimVertex(copy.fCutDcaNegToPrimVertex)
,fCutV0RadiusMin(copy.fCutV0RadiusMin)
,fCutV0RadiusMax(copy.fCutV0RadiusMax)
,fCutBetheBloch(copy.fCutBetheBloch)
,fK0Type(copy.fK0Type)
,fLaType(copy.fLaType)          
,fALaType(copy.fALaType)   
,fListK0s(copy.fListK0s)
,fListLa(copy.fListLa)
,fListALa(copy.fListALa)
,fListK0sCone(copy.fListK0sCone)
,fListLaCone(copy.fListLaCone)
,fListALaCone(copy.fListALaCone)
,fPIDResponse(copy.fPIDResponse)
,rndm(copy.rndm)
,fInputEntries(copy.fInputEntries)
,fAODPathArray(copy.fAODPathArray)
,fAODEntriesArray(copy.fAODEntriesArray)
,fAODPath(copy.fAODPath)
,fAODEntries(copy.fAODEntries)
,fAODEntriesSum(copy.fAODEntriesSum)
,fAODEntriesMax(copy.fAODEntriesMax)
,fOfflineTrgMask(copy.fOfflineTrgMask)
,fMinContribVtx(copy.fMinContribVtx)
,fVtxZMin(copy.fVtxZMin)
,fVtxZMax(copy.fVtxZMax)
,fEvtClassMin(copy.fEvtClassMin)
,fEvtClassMax(copy.fEvtClassMax)
,fCentMin(copy.fCentMin)
,fCentMax(copy.fCentMax)
,fNInputTracksMin(copy.fNInputTracksMin)
,fNInputTracksMax(copy.fNInputTracksMax)
,fTrackBranch(copy.fTrackBranch)
,fMCparticlesBranch(copy.fMCparticlesBranch)
,fJetBranch(copy.fJetBranch)
,fFileId(copy.fFileId)
,fAODEntry(copy.fAODEntry)
,fCountEvents(copy.fCountEvents)
,fEmbedMode(copy.fEmbedMode)
,fEvtSelecMode(copy.fEvtSelecMode)
,fEvtSelMinJetPt(copy.fEvtSelMinJetPt)
,fEvtSelMaxJetPt(copy.fEvtSelMaxJetPt)
,fEvtSelMinJetEta(copy.fEvtSelMinJetEta)
,fEvtSelMaxJetEta(copy.fEvtSelMaxJetEta)
,fEvtSelMinJetPhi(copy.fEvtSelMinJetPhi)
,fEvtSelMaxJetPhi(copy.fEvtSelMaxJetPhi)
,fEvtSelJetMinLConstPt(copy.fEvtSelJetMinLConstPt)
,fExtraEffPb(copy.fExtraEffPb)
,fToyMinNbOfTracks(copy.fToyMinNbOfTracks)
,fToyMaxNbOfTracks(copy.fToyMaxNbOfTracks)
,fToyMinTrackPt(copy.fToyMinTrackPt)
,fToyMaxTrackPt(copy.fToyMaxTrackPt)
,fToyDistributionTrackPt(copy.fToyDistributionTrackPt)
,fToyMinTrackEta(copy.fToyMinTrackEta)
,fToyMaxTrackEta(copy.fToyMaxTrackEta)
,fToyMinTrackPhi(copy.fToyMinTrackPhi)
,fToyMaxTrackPhi(copy.fToyMaxTrackPhi)
,fToyFilterMap(copy.fToyFilterMap)
,fTrackFilterMap(copy.fTrackFilterMap)
,fNPtHard(copy.fNPtHard)
,fPtHard(copy.fPtHard)
,fPtHardBin(copy.fPtHardBin)
,fAODJets(copy.fAODJets)
,fNevents(copy.fNevents)
,fXsection(copy.fXsection)
,fAvgTrials(copy.fAvgTrials)
,fHistList(copy.fHistList)
,fHistEvtSelection(copy.fHistEvtSelection)
,fh1Xsec(copy.fh1Xsec)
,fh1Trials(copy.fh1Trials)
,fh1TrialsEvtSel(copy.fh1TrialsEvtSel)
,fh2PtHard(copy.fh2PtHard)
,fh2PtHardEvtSel(copy.fh2PtHardEvtSel)
,fh2PtHardTrials(copy.fh2PtHardTrials)
,fh1TrackPt(copy.fh1TrackPt)
,fh2TrackEtaPhi(copy.fh2TrackEtaPhi)
,fh1TrackN(copy.fh1TrackN)
,fh1V0Pt(copy.fh1V0Pt)
,fh2V0EtaPhi(copy.fh2V0EtaPhi)
  //,fh1V0N(copy.fh1V0N)
,fh1K0Pt(copy.fh1K0Pt)
,fh2K0EtaPhi(copy.fh2K0EtaPhi)
,fh2K0sPtJetPtCone(copy.fh2K0sPtJetPtCone)
,fh1LaPt(copy.fh1LaPt)
,fh2LaEtaPhi(copy.fh2LaEtaPhi)
,fh2LaPtJetPtCone(copy.fh2LaPtJetPtCone)
,fh1ALaPt(copy.fh1ALaPt)
,fh2ALaEtaPhi(copy.fh2ALaEtaPhi)
,fh2ALaPtJetPtCone(copy.fh2ALaPtJetPtCone)
,fh1JetPt(copy.fh1JetPt)
,fh2JetEtaPhi(copy.fh2JetEtaPhi)
,fh1JetN(copy.fh1JetN)
,fh1MCTrackPt(copy.fh1MCTrackPt)
,fh2MCTrackEtaPhi(copy.fh2MCTrackEtaPhi)
,fh1MCTrackN(copy.fh1MCTrackN)
,fh1AODfile(copy.fh1AODfile)
,fh2AODevent(copy.fh2AODevent)
,fh1trackPosEta(copy.fh1trackPosEta)            
,fh1trackNegEta(copy.fh1trackNegEta)
,fh1V0Eta(copy.fh1V0Eta)              
,fh1CosPointAngle(copy.fh1CosPointAngle)           
,fh1DecayLengthV0(copy.fh1DecayLengthV0)  
,fh2ProperLifetimeK0sVsPtBeforeCut(copy.fh2ProperLifetimeK0sVsPtBeforeCut) 
,fh2ProperLifetimeK0sVsPtAfterCut(copy.fh2ProperLifetimeK0sVsPtAfterCut)   
,fh1V0Radius(copy.fh1V0Radius)          
,fh1DcaV0Daughters(copy.fh1DcaV0Daughters)        
,fh1DcaPosToPrimVertex(copy.fh1DcaPosToPrimVertex)   
,fh1DcaNegToPrimVertex(copy.fh1DcaNegToPrimVertex)    
,fh2ArmenterosBeforeCuts(copy.fh2ArmenterosBeforeCuts)
,fh2ArmenterosAfterCuts(copy.fh2ArmenterosAfterCuts)
,fh2BBLaPos(copy.fh2BBLaPos)
,fh2BBLaNeg(copy.fh2BBLaPos)
{
   // copy constructor
}


//__________________________________________________________________________
AliAnalysisTaskFastEmbedding& AliAnalysisTaskFastEmbedding::operator=(const AliAnalysisTaskFastEmbedding& o)
{
   // assignment

   if(this!=&o){
      AliAnalysisTaskSE::operator=(o);
      fESD               = o.fESD;
      fAODout            = o.fAODout;
      fAODevent          = o.fAODevent;
      fAODtree           = o.fAODtree;
      fAODfile           = o.fAODfile;
      mcHeader           = o.mcHeader;
      fFFRadius          = o.fFFRadius;
      fCutPostrackEta    = o.fCutPostrackEta;
      fCutNegtrackEta    = o.fCutNegtrackEta;  
      fCutEta            = o.fCutEta;
      fCutV0cosPointAngle  = o.fCutV0cosPointAngle;
      fKinkDaughters     = o.fKinkDaughters;
      fRequireTPCRefit   = o.fRequireTPCRefit;
      fCutArmenteros     = o.fCutArmenteros;
      fCutV0DecayMin     = o.fCutV0DecayMin;
      fCutV0DecayMax     = o.fCutV0DecayMax;
      fCutV0totMom       = o.fCutV0totMom;
      fCutDcaV0Daughters = o.fCutDcaV0Daughters;
      fCutDcaPosToPrimVertex = o.fCutDcaPosToPrimVertex;
      fCutDcaNegToPrimVertex = o.fCutDcaNegToPrimVertex;
      fCutV0RadiusMin    = o.fCutV0RadiusMin;
      fCutV0RadiusMax    = o.fCutV0RadiusMax;
      fCutBetheBloch     = o.fCutBetheBloch; 
      fK0Type            = o.fK0Type;
      fLaType            = o.fLaType;
      fALaType           = o.fALaType;
      fListK0s           = o.fListK0s;
      fListLa            = o.fListLa;
      fListALa           = o.fListALa;
      fListK0sCone       = o.fListK0sCone;
      fListLaCone        = o.fListLaCone;
      fListALaCone       = o.fListALaCone;
      fPIDResponse       = o.fPIDResponse;
      rndm               = o.rndm;
      fInputEntries      = o.fInputEntries;
      fAODPathArray      = o.fAODPathArray;
      fAODEntriesArray   = o.fAODEntriesArray;
      fAODPath           = o.fAODPath;
      fAODEntries        = o.fAODEntries;
      fAODEntriesSum     = o.fAODEntriesSum;
      fAODEntriesMax     = o.fAODEntriesMax;
      fOfflineTrgMask    = o.fOfflineTrgMask;
      fMinContribVtx     = o.fMinContribVtx;
      fVtxZMin           = o.fVtxZMin;
      fVtxZMax           = o.fVtxZMax;
      fEvtClassMin       = o.fEvtClassMin;
      fEvtClassMax       = o.fEvtClassMax;
      fCentMin           = o.fCentMin;
      fCentMax           = o.fCentMax;
      fNInputTracksMin   = o.fNInputTracksMin;
      fNInputTracksMax   = o.fNInputTracksMax;
      fTrackBranch       = o.fTrackBranch;
      fMCparticlesBranch = o.fMCparticlesBranch;
      fJetBranch         = o.fJetBranch;
      fFileId            = o.fFileId;
      fAODEntry          = o.fAODEntry;
      fCountEvents       = o.fCountEvents;
      fEmbedMode         = o.fEmbedMode;
      fEvtSelecMode      = o.fEvtSelecMode;
      fEvtSelMinJetPt    = o.fEvtSelMinJetPt;
      fEvtSelMaxJetPt    = o.fEvtSelMaxJetPt;
      fEvtSelMinJetEta   = o.fEvtSelMinJetEta;
      fEvtSelMaxJetEta   = o.fEvtSelMaxJetEta;
      fEvtSelMinJetPhi   = o.fEvtSelMinJetPhi;
      fEvtSelMaxJetPhi   = o.fEvtSelMaxJetPhi;
      fEvtSelJetMinLConstPt = o.fEvtSelJetMinLConstPt;
      fExtraEffPb        = o.fExtraEffPb;
      fToyMinNbOfTracks  = o.fToyMinNbOfTracks;
      fToyMaxNbOfTracks  = o.fToyMaxNbOfTracks;
      fToyMinTrackPt     = o.fToyMinTrackPt;
      fToyMaxTrackPt     = o.fToyMaxTrackPt;
      fToyDistributionTrackPt = o.fToyDistributionTrackPt;
      fToyMinTrackEta    = o.fToyMinTrackEta;
      fToyMaxTrackEta    = o.fToyMaxTrackEta;
      fToyMinTrackPhi    = o.fToyMinTrackPhi;
      fToyMaxTrackPhi    = o.fToyMaxTrackPhi;
      fToyFilterMap      = o.fToyFilterMap;
      fTrackFilterMap    = o.fTrackFilterMap;
      fNPtHard           = o.fNPtHard;
      fPtHard            = o.fPtHard;
      fPtHardBin         = o.fPtHardBin;
      fAODJets           = o.fAODJets;
      fNevents           = o.fNevents;
      fXsection          = o.fXsection;
      fAvgTrials         = o.fAvgTrials;
      fHistList          = o.fHistList;
      fHistEvtSelection  = o.fHistEvtSelection;
      fh1Xsec            = o.fh1Xsec;
      fh1Trials          = o.fh1Trials;
      fh1TrialsEvtSel    = o.fh1TrialsEvtSel;
      fh2PtHard          = o.fh2PtHard;
      fh2PtHardEvtSel    = o.fh2PtHardEvtSel;
      fh2PtHardTrials    = o.fh2PtHardTrials;
      fh1TrackPt         = o.fh1TrackPt;
      fh2TrackEtaPhi     = o.fh2TrackEtaPhi;
      fh1TrackN          = o.fh1TrackN;
      fh1V0Pt            = o.fh1V0Pt;
      fh2V0EtaPhi        = o.fh2V0EtaPhi;
      fh1K0Pt            = o.fh1K0Pt;
      fh2K0EtaPhi        = o.fh2K0EtaPhi;
      fh2K0sPtJetPtCone  = o.fh2K0sPtJetPtCone;
      fh1LaPt            = o.fh1LaPt;
      fh2LaEtaPhi        = o.fh2LaEtaPhi;
      fh2LaPtJetPtCone   = o.fh2LaPtJetPtCone;
      fh1ALaPt           = o.fh1ALaPt;
      fh2ALaEtaPhi       = o.fh2ALaEtaPhi;
      fh2ALaPtJetPtCone  = o.fh2ALaPtJetPtCone;
      fh1JetPt           = o.fh1JetPt;
      fh2JetEtaPhi       = o.fh2JetEtaPhi;
      fh1JetN            = o.fh1JetN;
      fh1MCTrackPt       = o.fh1MCTrackPt;
      fh2MCTrackEtaPhi   = o.fh2MCTrackEtaPhi;
      fh1MCTrackN        = o.fh1MCTrackN;
      fh1AODfile         = o.fh1AODfile;
      fh2AODevent        = o.fh2AODevent;
      fh1trackPosEta     = o.fh1trackPosEta;            
      fh1trackNegEta     = o.fh1trackNegEta;   
      fh1V0Eta           = o.fh1V0Eta;             
      fh1CosPointAngle   = o.fh1CosPointAngle;              
      fh1DecayLengthV0   = o.fh1DecayLengthV0;  
      fh2ProperLifetimeK0sVsPtBeforeCut = o.fh2ProperLifetimeK0sVsPtBeforeCut;
      fh2ProperLifetimeK0sVsPtAfterCut = o.fh2ProperLifetimeK0sVsPtAfterCut;
      fh1V0Radius        = o.fh1V0Radius;         
      fh1DcaV0Daughters  = o.fh1DcaV0Daughters;        
      fh1DcaPosToPrimVertex = o.fh1DcaPosToPrimVertex;   
      fh1DcaNegToPrimVertex = o.fh1DcaNegToPrimVertex;    
      fh2ArmenterosBeforeCuts = o.fh2ArmenterosBeforeCuts;
      fh2ArmenterosAfterCuts = o.fh2ArmenterosAfterCuts;
      fh2BBLaPos         = o.fh2BBLaPos;
      fh2BBLaNeg         = o.fh2BBLaPos;
   }

   return *this;
}


//__________________________________________________________________________
AliAnalysisTaskFastEmbedding::~AliAnalysisTaskFastEmbedding()
{
  // destructor
  if(fListK0s) delete fListK0s;
  if(fListLa) delete fListLa;
  if(fListALa) delete fListALa;
  if(fListK0sCone) delete fListK0sCone;
  if(fListLaCone) delete fListLaCone;
  if(fListALaCone) delete fListALaCone;
  delete rndm;
}


//__________________________________________________________________________
void AliAnalysisTaskFastEmbedding::UserCreateOutputObjects()
{
   // create output objects
   if(fDebug > 1) Printf("AliAnalysisTaskFastEmbedding::UserCreateOutputObjects()");
   AliLog::SetClassDebugLevel("AliAnalysisTaskFastEmbedding", AliLog::kInfo);

   OpenFile(1);
   if(!fHistList) fHistList = new TList();
   fHistList->SetOwner(kTRUE);
   
   fListK0s = new TList();
   fListK0s->SetOwner(kFALSE);

   fListLa = new TList();
   fListLa->SetOwner(kFALSE);

   fListALa = new TList();
   fListALa->SetOwner(kFALSE);

   fListK0sCone = new TList();
   fListK0sCone->SetOwner(kFALSE);

   fListLaCone = new TList();
   fListLaCone->SetOwner(kFALSE);

   fListALaCone = new TList();
   fListALaCone->SetOwner(kFALSE);


   // set seed
   rndm = new TRandom3();
   Int_t id = GetJobID();
   if(id>-1) rndm->SetSeed(id);
   else      rndm->SetSeed();   // a TTUID is generated and used for seed
   AliInfo(Form("TRandom3 seed: %d", rndm->GetSeed()));



   // embed mode with AOD
   if(fEmbedMode==kAODFull || fEmbedMode==kAODJetTracks || fEmbedMode==kAODJet4Mom){

      // open input AOD
      fFileId = OpenAODfile();
      fNevents = 0; // force to open another aod in UserExec()
      if(fFileId<0){
         AliError("");
         PostData(1, fHistList);
         return;
      }
   } //end: embed mode with AOD


   // connect output aod
   // create a new branch for extra tracks
   fAODout = AODEvent();
   if(!fAODout){
      AliError("Output AOD not found.");
      PostData(1, fHistList);
      return;
   }
   if(!fAODout->FindListObject(fTrackBranch.Data()) && strlen(fTrackBranch.Data())){//aodExtraTracks
     AliInfo(Form("Add AOD branch %s", fTrackBranch.Data()));
      TClonesArray *tracks = new TClonesArray("AliAODTrack",0);
      tracks->SetName(fTrackBranch.Data());
      AddAODBranch("TClonesArray", &tracks);
   }
   if(!fAODout->FindListObject("aodExtraV0s")){
      AliInfo("Add AOD branch embedded V0s");
      TClonesArray *extrav0s = new TClonesArray("AliAODv0",0);
      extrav0s->SetName("aodExtraV0s");
      AddAODBranch("TClonesArray", &extrav0s);
      }
   if(!fAODout->FindListObject("aodExtraK0s")){
      AliInfo("Add AOD branch embedded K0s candidates");
      TClonesArray *extraK0s = new TClonesArray("AliAODv0",0);
      extraK0s->SetName("aodExtraK0s");
      AddAODBranch("TClonesArray", &extraK0s);
      }
  if(!fAODout->FindListObject("aodExtraLa")){
      AliInfo("Add AOD branch embedded Lambda candidates");
      TClonesArray *extraLa = new TClonesArray("AliAODv0",0);
      extraLa->SetName("aodExtraLa");
      AddAODBranch("TClonesArray", &extraLa);
      }
  if(!fAODout->FindListObject("aodExtraALa")){
      AliInfo("Add AOD branch embedded Antilambda candidates");
      TClonesArray *extraALa = new TClonesArray("AliAODv0",0);
      extraALa->SetName("aodExtraALa");
      AddAODBranch("TClonesArray", &extraALa);
      }

   //for the v0 particle candidates inside the jets (only in kAODJetTracks mode)
   if(!fAODout->FindListObject("aodExtraK0sCone")){
      AliInfo("Add AOD branch embedded K0s Cone candidates");
      TClonesArray *extraK0sCone = new TClonesArray("AliAODv0",0);
      extraK0sCone->SetName("aodExtraK0sCone");
      AddAODBranch("TClonesArray", &extraK0sCone);
      }
   
   if(!fAODout->FindListObject("aodExtraLaCone")){
     AliInfo("Add AOD branch embedded Lambda Cone candidates");
     TClonesArray *extraLaCone = new TClonesArray("AliAODv0",0);
     extraLaCone->SetName("aodExtraLaCone");
     AddAODBranch("TClonesArray", &extraLaCone);
   }   
   if(!fAODout->FindListObject("aodExtraALaCone")){
     AliInfo("Add AOD branch embedded Antilambda Cone candidates");
    TClonesArray *extraALaCone = new TClonesArray("AliAODv0",0);
    extraALaCone->SetName("aodExtraALaCone");
    AddAODBranch("TClonesArray", &extraALaCone);
   }
   // create new branch for extra mcparticle if available as input
   if(fAODevent && fAODevent->FindListObject("mcparticles") && strlen(fMCparticlesBranch.Data())){
     AliInfo(Form("Add AOD branch %s", fMCparticlesBranch.Data()));
     TClonesArray *mcparticles = new TClonesArray("AliAODMCParticle",0);
     mcparticles->SetName(fMCparticlesBranch.Data());
     AddAODBranch("TClonesArray", &mcparticles);
   }
   


   //qa histograms
   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   
   fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 6, -0.5, 5.5);
   fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
   fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
   fHistEvtSelection->GetXaxis()->SetBinLabel(3,"event selection (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(5,"centrality (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(6,"multiplicity (rejected)");
   
   fh1Xsec         = new TProfile("fh1Xsec","xsec from pyxsec.root;p_{T,hard} bin;<#sigma>",fNPtHard+1,-1.5,fNPtHard-0.5);
   fh1Trials       = new TH1F("fh1Trials","trials (simulation) from pyxsec.root;p_{T,hard} bin;#sum{ntrials}",fNPtHard+1,-1.5,fNPtHard-0.5);
   fh1TrialsEvtSel = new TH1F("fh1TrialsEvtSel","trials (event selection) from pyxsec.root;p_{T,hard} bin;#sum{ntrials}",fNPtHard+1,-1.5,fNPtHard-0.5);
   fh2PtHard       = new TH2F("fh2PtHard","PYTHIA Pt hard;p_{T,hard} bin;p_{T,hard}",fNPtHard+1,-1.5,fNPtHard-0.5,350,-.5,349.5);
   fh2PtHardEvtSel = new TH2F("fh2PtHardEvtSel","PYTHIA Pt hard (event selection);p_{T,hard} bin;p_{T,hard}",fNPtHard+1,-1.5,fNPtHard-0.5,350,-.5,349.5);
   fh2PtHardTrials = new TH2F("fh2PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard} bin;#sum{p_{T,hard}}",fNPtHard+1,-1.5,fNPtHard-0.5,350,-.5,349.5);
   
   fHistList->Add(fHistEvtSelection);
   fHistList->Add(fh1Xsec);
   fHistList->Add(fh1Trials);
   fHistList->Add(fh1TrialsEvtSel);
   fHistList->Add(fh2PtHard);
   fHistList->Add(fh2PtHardEvtSel);
   fHistList->Add(fh2PtHardTrials);

   fh1TrackPt      =  new TH1F("fh1TrackPt","pT of extra tracks;p_{T};entries", 120, 0., 12.);
   fh2TrackEtaPhi  =  new TH2F("fh2TrackEtaPhi","eta-phi distribution of extra tracks;#eta;#phi", 20, -1., 1., 60, 0., 2*TMath::Pi());
   fh1TrackN       =  new TH1F("fh1TrackN", "nb. of extra tracks per event;nb. of tracks;entries",300, 0., 300.);
   fh1V0Pt         =  new TH1F("fh1V0Pt","pT of extra V0s from PYTHIA file;p_{T};entries", 120, 0., 12.);
   fh2V0EtaPhi     =  new TH2F("fh2V0EtaPhi","eta-phi distribution of extra V0s;#eta;#phi", 20, -1., 1., 60, 0., 2*TMath::Pi());
   fh1K0Pt         =  new TH1F("fh1K0Pt","pT of extra associated K0s;p_{T};entries", 120, 0., 12.);
   fh2K0EtaPhi     =  new TH2F("fh2K0EtaPhi","eta-phi distribution of extra V0s;#eta;#phi", 20, -1., 1., 60, 0., 2*TMath::Pi());
   fh2K0sPtJetPtCone  =  new TH2F("fh2K0sPtJetPtCone","v0-pt jet pt distribution of extra K0s;#it{p_{T}^{jet,ch}}};#it{p_{T}}", 20, 0., 100.,120, 0., 12.);
   fh1LaPt         =  new TH1F("fh1LaPt","pT of extra associated La;p_{T};entries", 120, 0., 12.);
   fh2LaEtaPhi     =  new TH2F("fh2LaEtaPhi","eta-phi distribution of extra associated #Lamdba;#eta;#phi", 20, -1., 1., 60, 0., 2*TMath::Pi());
   fh2LaPtJetPtCone  =  new TH2F("fh2LaPtJetPtCone","v0-pt jet pt distribution of extra associated #Lamdba;#it{p_{T}^{jet,ch}}};#it{p_{T}}", 20, 0., 100.,120, 0., 12.);
   fh1ALaPt        =  new TH1F("fh1ALaPt","pT of extra associated #bar{#Lamdba};p_{T};entries", 120, 0., 12.);
   fh2ALaEtaPhi    =  new TH2F("fh2ALaEtaPhi","eta-phi distribution of extra associated #bar{#Lamdba};#eta;#phi", 20, -1., 1., 60, 0., 2*TMath::Pi());
   fh2ALaPtJetPtCone  =  new TH2F("fh2ALaPtJetPtCone","v0-pt jet pt distribution of extra associated #bar{#Lamdba};#it{p_{T}^{jet,ch}}};#it{p_{T}}", 20, 0., 100.,120, 0., 12.);

   fHistList->Add(fh1TrackPt);
   fHistList->Add(fh2TrackEtaPhi);
   fHistList->Add(fh1TrackN);
   fHistList->Add(fh1V0Pt);
   fHistList->Add(fh2V0EtaPhi);
   fHistList->Add(fh1K0Pt);
   fHistList->Add(fh2K0EtaPhi);
   fHistList->Add(fh2K0sPtJetPtCone);
   fHistList->Add(fh1LaPt);
   fHistList->Add(fh2LaEtaPhi);
   fHistList->Add(fh2LaPtJetPtCone);
   fHistList->Add(fh1ALaPt);
   fHistList->Add(fh2ALaEtaPhi);
   fHistList->Add(fh2ALaPtJetPtCone);
   fHistList->Add(fh1V0Eta);         
   fHistList->Add(fh1CosPointAngle);                      
   fHistList->Add(fh1DecayLengthV0); 
   fHistList->Add(fh2ProperLifetimeK0sVsPtBeforeCut);
   fHistList->Add(fh2ProperLifetimeK0sVsPtAfterCut);
   fHistList->Add(fh1V0Radius);     
   fHistList->Add(fh1DcaV0Daughters);        
   fHistList->Add(fh1DcaPosToPrimVertex);   
   fHistList->Add(fh1DcaNegToPrimVertex);    
   fHistList->Add(fh2ArmenterosBeforeCuts);
   fHistList->Add(fh2ArmenterosAfterCuts);
   fHistList->Add(fh2BBLaPos);
   fHistList->Add(fh2BBLaNeg);
  
   if(fEmbedMode==kAODFull || fEmbedMode==kAODJetTracks || fEmbedMode==kAODJet4Mom){
      
      fh1JetPt        =  new TH1F("fh1JetPt", "pT of extra jets;p_{T};entries", 120, 0., 120.);
      fh2JetEtaPhi    =  new TH2F("fh2JetEtaPhi", "eta-phi distribution of extra jets;#eta;#phi", 20, -1., 1., 60, 0., 2*TMath::Pi());
      fh1JetN         =  new TH1F("fh1JetN", "nb. of extra jets per event;nb. of jets;entries",20,0.,20.);
      
      fHistList->Add(fh1JetPt);
      fHistList->Add(fh2JetEtaPhi);
      fHistList->Add(fh1JetN);
   }


   if(fAODevent && fAODevent->FindListObject("mcparticles") && strlen(fMCparticlesBranch.Data())){ 

      fh1MCTrackPt      =  new TH1F("fh1MCTrackPt","pT of MC extra tracks;p_{T};entries", 250, 0., 250.);
      fh2MCTrackEtaPhi  =  new TH2F("fh2MCTrackEtaPhi","eta-phi distribution of MC extra tracks;#eta;#phi", 20, -1., 1., 60, 0., 2*TMath::Pi());
      fh1MCTrackN       =  new TH1F("fh1MCTrackN", "nb. of MC extra tracks per event;nb. of tracks;entries",300, 0., 300.);
      
      fHistList->Add(fh1MCTrackPt);
      fHistList->Add(fh2MCTrackEtaPhi);
      fHistList->Add(fh1MCTrackN);
      
   }
   
   fh1AODfile = new TH1I("fh1AODfile", "overview of opened AOD files from the array", 23, -0.5, 2299.5);
   fh2AODevent = new TH2I("fh1AODevent","selected events;file;event", 25,-.5,2499.5,50,-.5,4999.5);
   fHistList->Add(fh1AODfile);
   fHistList->Add(fh2AODevent);

   // =========== Switch on Sumw2 for all histos ===========
   for (Int_t i=0; i<fHistList->GetEntries(); ++i) {
      TH1 *h1 = dynamic_cast<TH1*>(fHistList->At(i));
      if (h1){
         h1->Sumw2();
         continue;
      }
   }

   TH1::AddDirectory(oldStatus);

   PostData(1, fHistList);
}


//__________________________________________________________________________
void AliAnalysisTaskFastEmbedding::Init()
{
   // Initialization
   if(fDebug > 1) Printf("AliAnalysisTaskFastEmbedding::Init()");

}

//__________________________________________________________________________
void AliAnalysisTaskFastEmbedding::CalculateInvMass(AliAODv0* v0vtx, const Int_t particletype, Double_t& invM, Double_t& trackPt){ 
  
   //particletype:
   // * kaon = 1
   // * lambda = 2
   // * antilambda = 3

     invM    = 0;
     trackPt = 0;
   
     Double_t pp[3]={0,0,0}; //3-momentum positive charged track
     Double_t pm[3]={0,0,0}; //3-momentum negative charged track

     const Double_t massPi = 0.13957018; //better use PDG code at this point
     const Double_t massP  = 0.93827203;
    
     Double_t mass1=0;  
     Double_t mass2=0;
    
     TLorentzVector vector;  //lorentzvector V0 particle
     TLorentzVector fourmom1;//lorentzvector positive daughter
     TLorentzVector fourmom2;//lorentzvector negative daughter
  
   //--------------------------------------------------------------
   
   AliAODTrack *trackPos = (AliAODTrack *) (v0vtx->GetSecondaryVtx()->GetDaughter(0));//index 0 defined as positive charged track in AliESDFilter
  
   if( trackPos->Charge() == 1 ){
    
     pp[0]=v0vtx->MomPosX();
     pp[1]=v0vtx->MomPosY();
     pp[2]=v0vtx->MomPosZ();

     pm[0]=v0vtx->MomNegX();
     pm[1]=v0vtx->MomNegY();
     pm[2]=v0vtx->MomNegZ();
   }

   if( trackPos->Charge() == -1 ){ 
    
     pm[0]=v0vtx->MomPosX();
     pm[1]=v0vtx->MomPosY();
     pm[2]=v0vtx->MomPosZ();

     pp[0]=v0vtx->MomNegX();
     pp[1]=v0vtx->MomNegY();
     pp[2]=v0vtx->MomNegZ();
   }

   if (particletype == kK0){ // case K0s 
     mass1 = massPi;//positive particle
     mass2 = massPi;//negative particle
   } else if (particletype == kLambda){ // case Lambda
     mass1 = massP;//positive particle
     mass2 = massPi;//negative particle
   } else if (particletype == kAntiLambda){ //case AntiLambda
     mass1 = massPi;//positive particle
     mass2 = massP; //negative particle
   }
  
   fourmom1.SetXYZM(pp[0],pp[1],pp[2],mass1);//positive track
   fourmom2.SetXYZM(pm[0],pm[1],pm[2],mass2);//negative track
   vector=fourmom1 + fourmom2;
  
   invM    = vector.M(); 
   trackPt = vector.Pt();

   return;
}


//__________________________________________________________________________
Bool_t AliAnalysisTaskFastEmbedding::AcceptBetheBloch(AliAODv0 *v0, AliPIDResponse *PIDResponse, const Int_t particletype) //dont use for MC Analysis
{ 
  
	Int_t nnum = 1; 
	Int_t pnum = 0;
       
	const AliAODTrack *ntracktest=(AliAODTrack *)v0->GetDaughter(nnum); 
	if(ntracktest->Charge() > 0){nnum = 0; pnum = 1;}
	
	const AliAODTrack *trackNeg=(AliAODTrack *)(v0->GetDaughter(nnum));
	const AliAODTrack *trackPos=(AliAODTrack *)(v0->GetDaughter(pnum));
       
	//Check if both tracks are available
	if (!trackPos || !trackNeg) {
	  Printf("AliAnalysisTaskFastEmbedding::AcceptBetheBloch:: Error:Could not retrieve one of the daughter tracks\n");
	  return kFALSE;
	}
	
	//remove like sign V0s
	if ( trackPos->Charge() == trackNeg->Charge() ){
	  return kFALSE;
	  }  

        Double_t nsig_p = 0; //number of sigmas that positive daughter track has got in TPC pid information
        Double_t nsig_n = 0;

        const AliAODPid *pid_p=trackPos->GetDetPid();  // returns fDetPID, more detailed or detector specific pid information
        const AliAODPid *pid_n=trackNeg->GetDetPid();

	if(!pid_p)return kFALSE;
	if(!pid_n)return kFALSE;
	
        if (pid_p)
	  {
	    if(particletype == 1) //PID cut on positive charged Lambda daughters (only those with pt < 1 GeV/c)
	      {		
		nsig_p=PIDResponse->NumberOfSigmasTPC(trackPos,AliPID::kProton);
		Double_t protonPt = trackPos->Pt();
		if ((TMath::Abs(nsig_p) >= fCutBetheBloch) && (fCutBetheBloch >0) && (protonPt < 1)) return kFALSE;     
	      }	    
	  }
	
        if (pid_n)
	  {
	    if(particletype == 2)
	      {	
		nsig_n=PIDResponse->NumberOfSigmasTPC(trackNeg,AliPID::kProton);
		Double_t antiprotonPt = trackNeg->Pt();
	        if ((TMath::Abs(nsig_n) >= fCutBetheBloch) && (fCutBetheBloch >0) && (antiprotonPt < 1)) return kFALSE;
	      }
	    	      
	  }

        return kTRUE;
}


//___________________________________________________________________
Bool_t AliAnalysisTaskFastEmbedding::IsK0InvMass(const Double_t mass) const
{
  // K0 mass ?
  
  if(0.3 <= mass && mass < 0.7) return kTRUE;

  return kFALSE;
}
//___________________________________________________________________
Bool_t AliAnalysisTaskFastEmbedding::IsLaInvMass(const Double_t mass) const
{
  // La mass ? Use FF histo limits

  
  if(1.05 <= mass && mass < 1.25) return kTRUE;

  return kFALSE;
}
//__________________________________________________________________________
Bool_t AliAnalysisTaskFastEmbedding::DaughterTrackCheck(AliAODv0* v0, Int_t& nnum, Int_t& pnum){
  

  if(v0->GetNDaughters() != 2) return kFALSE;//case v0 has more or less than 2 daughters, avoids seg. break at some AOD files                                 
  // safety check of input parameters
  if(v0 == NULL)
    {
      if(fDebug > 1){std::cout << std::endl
		<< "Warning in AliAnalysisTaskFastEmbedding::DaughterTrackCheck:" << std::endl
			   << "v0 = " << v0 << std::endl;}

      return kFALSE;
    }
  else
    {
      //Daughters track check
      nnum = 1; 
      pnum = 0;
         
      AliAODTrack *ntracktest =(AliAODTrack*)(v0->GetDaughter(nnum));
  
      if(ntracktest == NULL)
	{
	 
	  return kFALSE;
	}

      if(ntracktest->Charge() > 0)
	{
	  nnum = 0; 
	  pnum = 1;
	}
  
      const AliAODTrack *trackNeg=(AliAODTrack *)(v0->GetDaughter(nnum));
      const AliAODTrack *trackPos=(AliAODTrack *)(v0->GetDaughter(pnum));
  
      //Check if both tracks are available
      if (!trackPos || !trackNeg) {

	return kFALSE;
      }
      
      //remove like sign V0s
      if ( trackPos->Charge() == trackNeg->Charge() ){

	return kFALSE;
      }  
  
      return kTRUE;
    }
}

//__________________________________________________________________________
Bool_t AliAnalysisTaskFastEmbedding::ApplyV0Cuts(AliAODv0* v0, const Int_t type, const Int_t particletype, AliAODVertex* primVertex, AliAODEvent* fAOD)
{

  if(type==kTrackUndef) return kFALSE;

  if(!primVertex) return kFALSE;// this is the vertex of the PYTHIA event
    
    Double_t lPrimaryVtxPosition[3];
    Double_t lV0Position[3];
    lPrimaryVtxPosition[0] = primVertex->GetX();
    lPrimaryVtxPosition[1] = primVertex->GetY();
    lPrimaryVtxPosition[2] = primVertex->GetZ();
       
       if(!v0)
	 {
	   std::cout << std::endl
		     << "Warning in AliAnalysisTaskFastEmbedding::ApplyV0Cuts:" << std::endl
		     << "v0 = " << v0 << std::endl;
	   return kFALSE;
	 }
       
       Bool_t isOnFly = v0->GetOnFlyStatus();
       
       if(!isOnFly &&  (type == kOnFly || type == kOnFlyPID || type == kOnFlydEdx || type == kOnFlyPrim)) return kFALSE; 
       if( isOnFly &&  (type == kOffl  || type == kOfflPID  || type == kOffldEdx  || type == kOfflPrim))  return kFALSE; 

       
       Double_t pp[3]={0,0,0}; //3-momentum positive charged track
       Double_t pm[3]={0,0,0}; //3-momentum negative charged track
       Double_t v0mom[3]={0,0,0};
       
       Double_t invM = 0;
       Double_t invMK0s=0;
       Double_t invMLa=0;
       Double_t invMALa=0;
       Double_t trackPt=0;
       Int_t nnum = -1;
       Int_t pnum = -1;
       
       //get MC stack

       TClonesArray *stackmc = 0x0;
       stackmc = dynamic_cast<TClonesArray*>(fAODevent->FindListObject(AliAODMCParticle::StdBranchName())); //get MCAOD branch in data
       if (!stackmc)
	 {
	   Printf("AliAnalysisTaskFastEmbedding::ERROR: stack not available");
	   return kFALSE;
	 }
	 

 
       Bool_t daughtercheck = DaughterTrackCheck(v0, nnum, pnum);
       
       if(daughtercheck == kFALSE)return kFALSE;
       
       const AliAODTrack *trackNeg=(AliAODTrack *)(v0->GetDaughter(nnum));
       const AliAODTrack *trackPos=(AliAODTrack *)(v0->GetDaughter(pnum));
       
       //Double_t trackPt=0;   
       Double_t fV0Radius      = -999;
       Double_t fDcaV0Daughters = v0->DcaV0Daughters();
       Double_t fDcaPosToPrimVertex = v0->DcaPosToPrimVertex();//IP of positive charged daughter
       Double_t fDcaNegToPrimVertex = v0->DcaNegToPrimVertex();//IP of negative charged daughter
       Int_t negDaughterpdg = 0;
       Int_t posDaughterpdg = 0;
       Int_t v0Label = -1;
       Double_t MCPt = 0;
       Bool_t fPhysicalPrimary = kFALSE;//don't use IsPhysicalPrimary() anymore for MC analysis, use instead 2D distance from primary to secondary vertex
       Int_t MCv0PDGCode = 0;     
       Int_t negAssLabel = TMath::Abs(trackNeg->GetLabel());                       //negative (reconstructed) charged track label in MC stack
       Int_t posAssLabel = TMath::Abs(trackPos->GetLabel());                       //positive (reconstructed) charged track label in MC stack

       Double_t mcXv = 0;
       Double_t mcYv = 0;
       Double_t mcZv = 0;
       
       TList *listmc = fAODevent->GetList(); //AliAODEvent* is inherited from AliVEvent*, listmc is pointer to reconstructed event in MC list, member of AliAODEvent
       if(!listmc)return kFALSE;
       
       AliAODMCHeader *header=(AliAODMCHeader*)listmc->FindObject(AliAODMCHeader::StdBranchName());
      
        if (!header)
	 {
	   Printf("AliAnalysisTaskFastEmbedding::ERROR: mc header not available!");
	   return kFALSE;
	 }
	 
       mcXv=header->GetVtxX(); mcYv=header->GetVtxY(); mcZv=header->GetVtxZ();       
       //mc label checks
       if(negAssLabel>=0 && negAssLabel < stackmc->GetEntriesFast() && posAssLabel>=0 && posAssLabel < stackmc->GetEntriesFast()){//safety check if label has valid value of stack	 
	 AliAODMCParticle *mcNegPart = dynamic_cast<AliAODMCParticle*>(stackmc->UncheckedAt(negAssLabel));//fetch the, with one MC truth track associated, negative charged track 
	 v0Label = mcNegPart->GetMother();// mother of negative charged particle is v0, get v0 label here
	 negDaughterpdg = mcNegPart->GetPdgCode();
	 AliAODMCParticle *mcPosPart = dynamic_cast<AliAODMCParticle*>(stackmc->UncheckedAt(posAssLabel));//fetch the, with one MC truth track associated (reconstructed), positive charged track 
	 Int_t v0PosLabel = mcPosPart->GetMother();                                        //get mother label of positive charged track label
	 posDaughterpdg = mcPosPart->GetPdgCode();
	 
	 if(v0Label >= 0 && v0Label < stackmc->GetEntriesFast() && v0Label == v0PosLabel){//first v0 mc label check, then: check if both daughters are stemming from same particle
	   
	   AliAODMCParticle *mcv0 = dynamic_cast<AliAODMCParticle*>(stackmc->UncheckedAt(v0Label));  //fetch MC ass. particle to v0 (mother of the both charged daughter tracks)
	   
	   Float_t fDistPrimaryMax = 0.01; // [cm] max distance of production point to the primary vertex (criterion for choice of MC particles considered as primary)
	   
	   // Get the distance between production point of the MC mother particle and the primary vertex
	 
	   Double_t dx = mcXv-mcv0->Xv();//mc primary vertex - mc particle production vertex 
	   Double_t dy = mcYv-mcv0->Yv();
	   Double_t dz = mcZv-mcv0->Zv();
	   
	   Float_t fDistPrimary = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
	   fPhysicalPrimary = kFALSE;//init
	   
	   fPhysicalPrimary = (fDistPrimary < fDistPrimaryMax);
	   if(!fPhysicalPrimary)return kFALSE;//instead of IsPhysicalPrimary()

	   MCv0PDGCode = mcv0->GetPdgCode();
	   
	   MCPt = mcv0->Pt();
	  
	   if(particletype == kK0){
	     
	     if(TMath::Abs(posDaughterpdg) != 211){return kFALSE;}
	     if(TMath::Abs(negDaughterpdg) != 211){return kFALSE;}
	     
	     if(MCv0PDGCode != 310)  {return kFALSE;}
	   }
	   
	   if(particletype == kLambda){
	     if(MCv0PDGCode != 3122)return kFALSE;//if particle is not Antilambda, v0 is rejected
	     if(posDaughterpdg != 2212)return kFALSE;
	     if(negDaughterpdg != -211)return kFALSE;    //pdg code check for Lambda daughters	    
	   }
	   
	   if(particletype == kAntiLambda){
	     if(MCv0PDGCode != -3122)return kFALSE;
	     if(posDaughterpdg != 211)return kFALSE;
	     if(negDaughterpdg !=-2212)return kFALSE;    //pdg code check for Antilambda daughters	     
	   }
	 }  
       }   
	 
       //calculate InvMass for every V0 particle assumption
       switch(particletype){
       case kK0: 
	 CalculateInvMass(v0, kK0, invM, trackPt); //function to calculate invMass with TLorentzVector class
	 invMK0s=invM;
	 break; 
       case kLambda: 
	 CalculateInvMass(v0, kLambda, invM, trackPt); 
	 invMLa=invM;
	 break;   
       case kAntiLambda: 
	 CalculateInvMass(v0, kAntiLambda, invM, trackPt); 
	 invMALa=invM; 
	 break;
       default: 
	 std::cout<<"***NO VALID PARTICLETYPE***"<<std::endl; 
	 return 0;   
       }
  
       if(fDebug>7){if(!(IsK0InvMass(invMK0s)) && !(IsLaInvMass(invMLa)) && !(IsLaInvMass(invMALa))){std::cout<<"AliAnalysisTaskFastEmbedding::ApplyV0Cuts: invM not in selected mass window "<<std::endl;}}
       
       if(!(IsK0InvMass(invMK0s)) && !(IsLaInvMass(invMLa)) && !(IsLaInvMass(invMALa)))return kFALSE; 

       Double_t PosEta = trackPos->Eta();
       Double_t NegEta = trackNeg->Eta();
      
       //DistOverTotMom_in_2D___________
       
       Float_t fMassK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
       Float_t fMassLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
       
       
       AliAODVertex* primVtx = fAOD->GetPrimaryVertex(); // get the primary vertex
       Double_t dPrimVtxPos[3]; // primary vertex position {x,y,z}
       primVtx->GetXYZ(dPrimVtxPos);
       
       Float_t fPtV0 = TMath::Sqrt(v0->Pt2V0()); // transverse momentum of V0
       Double_t dSecVtxPos[3]; // V0 vertex position {x,y,z}
       v0->GetSecondaryVtx(dSecVtxPos);
       Double_t dDecayPath[3];
       for (Int_t iPos = 0; iPos < 3; iPos++)
	 dDecayPath[iPos] = dSecVtxPos[iPos]-dPrimVtxPos[iPos]; // vector of the V0 path
       Float_t fDecLen2D = TMath::Sqrt(dDecayPath[0]*dDecayPath[0]+dDecayPath[1]*dDecayPath[1]); //transverse path length R
       Float_t fROverPt = fDecLen2D/fPtV0; // R/pT
       
       Float_t fMROverPtK0s = fMassK0s*fROverPt; // m*R/pT
       Float_t fMROverPtLambda = fMassLambda*fROverPt; // m*R/pT
    
       Double_t fEta = -999;
       Double_t fV0cosPointAngle = -999;
       Double_t fV0DecayLength = v0->DecayLengthV0(lPrimaryVtxPosition);
       
       Double_t fV0mom[3];
       
       fV0mom[0]=v0->MomV0X();
       fV0mom[1]=v0->MomV0Y();
       fV0mom[2]=v0->MomV0Z();
       
       Double_t fV0TotalMomentum = TMath::Sqrt(fV0mom[0]*fV0mom[0]+fV0mom[1]*fV0mom[1]+fV0mom[2]*fV0mom[2]);

       const Double_t K0sPDGmass = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass(); 
       const Double_t LambdaPDGmass = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
       
       Double_t fDistOverTotMomK0s = 0;
       Double_t fDistOverTotMomLa = 0;
     
       if(particletype == kK0){
	 
	 fDistOverTotMomK0s = fV0DecayLength * K0sPDGmass;
	 fDistOverTotMomK0s /= (fV0TotalMomentum+1e-10);
       }
       
       if((particletype == kLambda)||(particletype == kAntiLambda)){
	 
	 fDistOverTotMomLa = fV0DecayLength * LambdaPDGmass;
	 fDistOverTotMomLa /= (fV0TotalMomentum+1e-10);
       } 
   
       if(fRequireTPCRefit==kTRUE){//if kTRUE: accept only if daughter track is refitted in TPC!!
	 Bool_t isPosTPCRefit = (trackPos->AliAODTrack::IsOn(AliESDtrack::kTPCrefit));
	 Bool_t isNegTPCRefit = (trackNeg->AliAODTrack::IsOn(AliESDtrack::kTPCrefit));
	 if (!isPosTPCRefit)return kFALSE;
	 if (!isNegTPCRefit)return kFALSE;
       }
       
       if(fKinkDaughters==kFALSE){//if kFALSE: no acceptance of kink daughters
	 AliAODVertex* ProdVtxDaughtersPos = (AliAODVertex*) (trackPos->AliAODTrack::GetProdVertex());
	 Char_t isAcceptKinkDaughtersPos  = ProdVtxDaughtersPos->GetType();
	 if(isAcceptKinkDaughtersPos==AliAODVertex::kKink)return kFALSE;
	 
	 AliAODVertex* ProdVtxDaughtersNeg = (AliAODVertex*) (trackNeg->AliAODTrack::GetProdVertex());
	 Char_t isAcceptKinkDaughtersNeg  = ProdVtxDaughtersNeg->GetType();
	 if(isAcceptKinkDaughtersNeg==AliAODVertex::kKink)return kFALSE;
	 
       }
     
       Double_t avDecayLengthK0s = 2.6844;
       Double_t avDecayLengthLa = 7.89;
     
       fV0cosPointAngle = v0->CosPointingAngle(lPrimaryVtxPosition);
       lV0Position[0]= v0->DecayVertexV0X();  
       lV0Position[1]= v0->DecayVertexV0Y();  
       lV0Position[2]= v0->DecayVertexV0Z();  
       
       fV0Radius  = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
       
       if(particletype == kK0)         {//fRap = v0->RapK0Short();
	 fEta = v0->PseudoRapV0();}
       if(particletype == kLambda)     {//fRap = v0->RapLambda();
	 fEta = v0->PseudoRapV0();}
       if(particletype == kAntiLambda) {//fRap = v0->Y(-3122);
	 fEta = v0->PseudoRapV0();}
      
       //cut on 2D DistOverTransMom: (recommended from Iouri)
       if((particletype == kLambda)||(particletype == kAntiLambda)){if(fMROverPtLambda > (fCutV0DecayMax * avDecayLengthLa))return kFALSE;}//fCutV0DecayMax set to 5 in AddTask macro
       
       //Armenteros Podolanski Plot for K0s      
       Double_t ArmenterosAlpha=-999;  
       Double_t ArmenterosPt=-999;
       Double_t PosPl;
       Double_t NegPl;
 
       if(particletype == kK0){
	 
	 pp[0]=v0->MomPosX();
	 pp[1]=v0->MomPosY();
	 pp[2]=v0->MomPosZ();
	 
	 pm[0]=v0->MomNegX();
	 pm[1]=v0->MomNegY();
	 pm[2]=v0->MomNegZ();
	 
	 v0mom[0]=v0->MomV0X();
	 v0mom[1]=v0->MomV0Y();
	 v0mom[2]=v0->MomV0Z();
	 
	 TVector3 v0Pos(pp[0],pp[1],pp[2]);
	 TVector3 v0Neg(pm[0],pm[1],pm[2]);
	 TVector3 v0totMom(v0mom[0], v0mom[1], v0mom[2]); //vector for tot v0 momentum

	 PosPl = v0Pos.Dot(v0totMom)/v0totMom.Mag();  //transversal momentum of positive charged daughter track
	 NegPl = v0Neg.Dot(v0totMom)/v0totMom.Mag();  //transversal momentum of nergative charged daughter track
	 ArmenterosAlpha = 1.-2./(1+(PosPl/NegPl));  
	 ArmenterosPt= v0->PtArmV0();	 
       }      
       if((TMath::Abs(PosEta)>fCutPostrackEta) || (TMath::Abs(NegEta)>fCutNegtrackEta))return kFALSE;   //daughters pseudorapidity cut
       if(particletype == kK0){
       if (fV0cosPointAngle < fCutV0cosPointAngle)return kFALSE;                                       //cospointangle cut
       }
       if((particletype == kLambda)||(particletype == kAntiLambda)){
	 if (fV0cosPointAngle < 0.995)return kFALSE;                                       //cospointangle cut
       }
       if(TMath::Abs(fEta) > fCutEta) return kFALSE;                                                  //V0 Eta Cut
       if (fDcaV0Daughters > fCutDcaV0Daughters)return kFALSE;
       if ((fDcaPosToPrimVertex < fCutDcaPosToPrimVertex) || (fDcaNegToPrimVertex < fCutDcaNegToPrimVertex))return kFALSE;
       if ((fV0Radius < fCutV0RadiusMin) || (fV0Radius > fCutV0RadiusMax))return kFALSE;
       //cut on 2D DistOverTransMom
       if(particletype == kK0){//the cut on Lambdas you can find above
       if(fMROverPtK0s > (fCutV0DecayMax * avDecayLengthK0s))return kFALSE;
       }   
       return kTRUE;
}
	
//__________________________________________________________________________
void AliAnalysisTaskFastEmbedding::GetTracksInCone(TList* inputlist, TList* outputlist, const AliAODJet* jet,const Double_t radius, Double_t& sumPt)
{
  // fill list of V0 tracks in cone around jet axis  

  if(!inputlist)return;
  if(!jet)return;
  sumPt = 0;
 
  Double_t jetMom[3];
  jet->PxPyPz(jetMom);
  TVector3 jet3mom(jetMom);

  for (Int_t itrack=0; itrack<inputlist->GetSize(); itrack++){//loop over all K0s found in event

    AliVParticle* track = dynamic_cast<AliVParticle*>(inputlist->At(itrack));
    if(!track)continue;
    Double_t trackMom[3];
    track->PxPyPz(trackMom);
    TVector3 track3mom(trackMom);
  
    Double_t dR = jet3mom.DeltaR(track3mom);

    if(dR<radius){//fill all the V0s inside cone into outputlist, radius is return value of GetFFRadius() 
      outputlist->Add(track);
      
      sumPt += track->Pt();      
    }
  }

  outputlist->Sort();
  
}

//__________________________________________________________________________
Bool_t AliAnalysisTaskFastEmbedding::UserNotify()
{
   // User defined Notify(), called once
   if(fDebug > 1) Printf("AliAnalysisTaskFastEmbedding::UserNotify()");

   // get total nb of events in tree (of this subjob)
   AliInputEventHandler* inputHandler = (AliInputEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
   fInputEntries = inputHandler->GetTree()->GetEntriesFast();
   AliInfo(Form("Total nb. of events: %d", fInputEntries));

   return kTRUE;

}

//__________________________________________________________________________
void AliAnalysisTaskFastEmbedding::UserExec(Option_t *)
{
   // Analysis, called once per event
   if(fDebug > 1) Printf("AliAnalysisTaskFastEmbedding::UserExec()");

   if(!fAODout){
      AliError("Need output AOD, but is not connected."); 
      PostData(1, fHistList);
      return;
   }

   fESD=dynamic_cast<AliESDEvent*>(InputEvent());
   if (!fESD) {
      AliError("ESD not available");
      PostData(1, fHistList);
      return;
   }

   // -- event selection --
   fHistEvtSelection->Fill(1); // number of events before event selection

   // physics selection
   AliInputEventHandler* inputHandler = (AliInputEventHandler*)
   ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
   if(!(inputHandler->IsEventSelected() & fOfflineTrgMask)){
      if(fDebug) Printf(" Trigger Selection: event REJECTED ... ");
      fHistEvtSelection->Fill(2);
      PostData(1, fHistList);
      return;
   } 

   fPIDResponse = inputHandler->GetPIDResponse();
   
   if (!fPIDResponse){if(fDebug > 1) Printf("AliAnalysisTaskFastEmbedding::UserExec(): fPIDResponse does not exist!"); return;}
 
   // vertex selection
   AliAODVertex* primVertex = fAODout->GetPrimaryVertex();
   Int_t nTracksPrim = primVertex->GetNContributors();
   if ((nTracksPrim < fMinContribVtx) ||
         (primVertex->GetZ() < fVtxZMin) ||
         (primVertex->GetZ() > fVtxZMax) ){
      if(fDebug) Printf("%s:%d primary vertex z = %f: event REJECTED...",(char*)__FILE__,__LINE__,primVertex->GetZ());
      fHistEvtSelection->Fill(3);
      PostData(1, fHistList);
      return;
   }

   /** takes wrong value, since jet service tasks runs later
   // event class selection (from jet helper task)
   Int_t eventClass = AliAnalysisHelperJetTasks::EventClass();
   if(fDebug) Printf("Event class %d", eventClass);
   if (eventClass < fEvtClassMin || eventClass > fEvtClassMax){
   fHistEvtSelection->Fill(4);
   PostData(1, fHistList);
   return;
   }*/

   // centrality selection
   AliCentrality *cent = 0x0;
   Float_t centValue = 0.; 
   if(fESD) cent = fESD->GetCentrality();
   if(cent) centValue = cent->GetCentralityPercentile("V0M");
   if(fDebug) printf("centrality: %f\n", centValue);
   if (centValue < fCentMin || centValue > fCentMax){
      fHistEvtSelection->Fill(4);
      PostData(1, fHistList);
      return;
   }


   /*  ** not implemented **
   // multiplicity due to input tracks
   Int_t nInputTracks = GetNInputTracks();

   if (nInputTracks < fNInputTracksMin || (fNInputTracksMax > -1 && nInputTracks > fNInputTracksMax)){
   fHistEvtSelection->Fill(5);
   PostData(1, fHistList);
   return;
   }
   */

   fHistEvtSelection->Fill(0); // accepted events  
   // -- end event selection --

   // connect aod out
   TClonesArray *tracks = (TClonesArray*)(fAODout->FindListObject(fTrackBranch.Data()));
   if(!tracks){
      AliError("Extra track branch not found in output.");
      PostData(1, fHistList);
      return;
   }
   tracks->Delete();
   Int_t nAODtracks=0;

   TClonesArray *extrav0s = (TClonesArray*)(fAODout->FindListObject("aodExtraV0s"));
   if(!extrav0s){
      AliError("Extra v0s branch not found in output.");
      PostData(1, fHistList);
      return;
   }
   TClonesArray *extraK0s = (TClonesArray*)(fAODout->FindListObject("aodExtraK0s"));
   if(!extraK0s){
     AliError("Extra K0s branch not found in output.");
     PostData(1, fHistList);
     return;
   }   
   TClonesArray *extraLa = (TClonesArray*)(fAODout->FindListObject("aodExtraLa"));
   if(!extraLa){
     AliError("Extra La branch not found in output.");
     PostData(1, fHistList);
     return;
   }   
   TClonesArray *extraALa = (TClonesArray*)(fAODout->FindListObject("aodExtraALa"));
   if(!extraALa){
     AliError("Extra ALa branch not found in output.");
     PostData(1, fHistList);
     return;
   }

  TClonesArray *extraK0sCone = (TClonesArray*)(fAODout->FindListObject("aodExtraK0sCone"));
   if(!extraK0sCone){
     AliError("Extra K0s Cone branch not found in output.");
     PostData(1, fHistList);
     return;
   }
   
   TClonesArray *extraLaCone = (TClonesArray*)(fAODout->FindListObject("aodExtraLaCone"));
   if(!extraLaCone){
     AliError("Extra La Cone branch not found in output.");
     PostData(1, fHistList);
     return;
   }
   
   TClonesArray *extraALaCone = (TClonesArray*)(fAODout->FindListObject("aodExtraALaCone"));
   if(!extraALaCone){
     AliError("Extra ALa Cone branch not found in output.");
     PostData(1, fHistList);
     return;
   }

   extrav0s->Delete();
   extraK0s->Delete();
   extraLa->Delete();
   extraALa->Delete();
   extraK0sCone->Delete();
   extraLaCone->Delete();
   extraALaCone->Delete();

   Int_t nAODv0s=0;
   Int_t nAODK0s=0;
   Int_t nAODLa=0;
   Int_t nAODALa=0;

   TRef dummy;

   // === embed mode with AOD ===
   if(fEmbedMode==kAODFull || fEmbedMode==kAODJetTracks || fEmbedMode==kAODJet4Mom){
      if(!fAODevent){
         AliError("Need input AOD, but is not connected."); 
         PostData(1, fHistList);
         return;
      }

      //essential part of PYTHIA embedding -> randomly picked PYTHIA event is checked, if good it is attached to PbPb real data event in AliAOD.root output file, if not the next PYTHIA event is randomly picked and checked 

      Bool_t useEntry = kFALSE;
      while(!useEntry){  // protection needed, if no (PYTHIA?) event fulfills requirement

         fAODEntry++; // go to next event 
         fCountEvents++;
         if(fAODEntry>=fNevents) fAODEntry=0; 

         // open new aod file
         if(fCountEvents>=fNevents){ 
            if(!fAODPathArray){
               AliDebug(AliLog::kDebug, "No list of AODs, keep current AOD.");
            } 
            else {
               AliDebug(AliLog::kDebug, "Select new AOD file ...");

               fFileId = OpenAODfile();
               if(fFileId<0) {
                  PostData(1, fHistList);
                  return;
               }
               fh1AODfile->Fill(fFileId);
               if(fAODEntries!=fNevents){
                  AliError("File with list of AODs and file with nb. of entries inconsistent");
                  PostData(1, fHistList);
                  return;
               }
            }
         }

         // get next event
         fAODtree->GetEvent(fAODEntry);

         // get pt hard
         if(mcHeader){
            fPtHard = mcHeader->GetPtHard();
            GetPtHard(kTRUE,fPtHard); // make value available for other tasks (e.g. jet response task)
            AliDebug(AliLog::kDebug,Form("pt hard %.2f", fPtHard));
         }
         else{
            AliDebug(AliLog::kDebug,"No mcHeader");
         }
         fPtHardBin = GetPtHardBin(fPtHard);

         //Printf("pt hard (bin) %.2f (%d), xSection %.2e, trials %.2f",fPtHard, fPtHardBin, fXsection, fAvgTrials); 

         // fill event variables for each event
         fh1Xsec->Fill(fPtHardBin,fXsection);
         fh2PtHard->Fill(fPtHardBin,fPtHard);

         fh1Trials->Fill(fPtHardBin, fAvgTrials);
         fh1TrialsEvtSel->Fill(fPtHardBin);
         fh2PtHardTrials->Fill(fPtHardBin,fPtHard,fAvgTrials);

         // jet pt selection
         if(fEvtSelecMode==kEventsJetPt){
            Int_t nJets = fAODJets->GetEntries();
            //Printf("nb. jets: %d",nJets);
            for(Int_t ij=0; ij<nJets; ++ij){
               AliAODJet *jet = dynamic_cast<AliAODJet*>(fAODJets->At(ij));
               if(!jet) continue;

               if(   (fEvtSelMinJetPt==-1. || jet->Pt()>=fEvtSelMinJetPt)//AZ: selection of PYTHIA jets with minimum 5 GeV/c
                     && (fEvtSelMaxJetPt==-1. || jet->Pt()<=fEvtSelMaxJetPt)
                     && (jet->Eta()>=fEvtSelMinJetEta && jet->Eta()<=fEvtSelMaxJetEta)
                     && (jet->Phi()>=fEvtSelMinJetPhi && jet->Phi()<=fEvtSelMaxJetPhi)
		     && HasMinLConstPt(jet)){

                  useEntry = kTRUE;

		  if(fDebug > 10){
		    for(Int_t i=0; i<jet->GetRefTracks()->GetEntriesFast(); i++){

		    
		      AliVParticle* track = dynamic_cast<AliVParticle*>(jet->GetRefTracks()->At(i));
		      
		      if(!track) continue; 

		      Double_t jetMom[3];
		      jet->PxPyPz(jetMom);
		      TVector3 jet3mom(jetMom);

		      Double_t trackMom[3];
		      track->PxPyPz(trackMom);
		      TVector3 track3mom(trackMom);

		      //Double_t dR = jet3mom.DeltaR(track3mom);
		      
		    } 
		  }
	       }
	    }
	 }
         // no selection
         if(fEvtSelecMode==kEventsAll){
	   useEntry = kTRUE;
         }
	 
      }
      AliDebug(AliLog::kDebug,Form("Use entry %d from extra AOD.", fAODEntry));

      fh2PtHardEvtSel->Fill(fPtHardBin,fPtHard);
      fh2AODevent->Fill(fFileId,fAODEntry);

      TClonesArray *mcpartIN  = (TClonesArray*)(fAODevent->FindListObject("mcparticles"));
      TClonesArray *mcpartOUT = 0x0;
      if(mcpartIN){
         mcpartOUT = (TClonesArray*)(fAODout->FindListObject(fMCparticlesBranch.Data()));
         if(!mcpartOUT){
            AliError("Extra MC particles branch not found in output.");
            PostData(1, fHistList);
            return;
         }
         mcpartOUT->Delete();
      } else {
         AliInfo("No extra MC particles found.");
      }


      if(fEmbedMode==kAODFull || fEmbedMode==kAODJetTracks){ // take all tracks or jet tracks
         // loop over input tracks
         // add to output aod
         Int_t nTracks = 0;
	 Int_t nV0s = 0;
         nV0s = fAODevent->GetNumberOfV0s();
         Int_t nJets = fAODJets->GetEntries();
         Int_t nSelectedJets = 0;
         AliAODJet *leadJet = 0x0; // used for jet tracks only

	 for(Int_t it=0; it<nV0s; ++it){//for all PYTHIA embedding v0 particles
	   AliAODv0 *tmpv0 = 0x0;

	   if(fEmbedMode==kAODFull) tmpv0 = fAODevent->GetV0(it);
	   if(!tmpv0)continue;

	   Double_t rd=rndm->Uniform(0.,1.);
	   if(rd>fExtraEffPb) continue; 

	   new ((*extrav0s)[nAODv0s++]) AliAODv0(*tmpv0); 
	   dummy = (*extrav0s)[nAODv0s-1];

	   fh1V0Pt->Fill(tmpv0->Pt());
	   fh2V0EtaPhi->Fill(tmpv0->Eta(), tmpv0->Phi());	  
	 }

	 for(Int_t it=0; it<nV0s; ++it){//for K0s candidates, apply cuts, fill into inclusive K0s list
           AliAODv0 *tmpv0 = 0x0;

	   if(fEmbedMode==kAODFull) tmpv0 = fAODevent->GetV0(it);
	   if(!tmpv0)continue;

	   Double_t rd=rndm->Uniform(0.,1.);
	   if(rd>fExtraEffPb) continue; 
	   Bool_t IsGoodV0 = kFALSE;
	   IsGoodV0 = ApplyV0Cuts(tmpv0, fK0Type, kK0, primVertex, fAODevent);
          
	   if(IsGoodV0 == kFALSE)continue;
	   new ((*extraK0s)[nAODK0s++]) AliAODv0(*tmpv0); 
	   dummy = (*extraK0s)[nAODK0s-1];
	   fh1K0Pt->Fill(tmpv0->Pt());
	   fh2K0EtaPhi->Fill(tmpv0->Eta(), tmpv0->Phi());
	   fListK0s->Add(tmpv0);
	   
	 }

	 for(Int_t it=0; it<nV0s; ++it){//for Lamdba candidates apply cuts, fill into inclusive La list
	   AliAODv0 *tmpv0 = 0x0;

	   if(fEmbedMode==kAODFull) tmpv0 = fAODevent->GetV0(it);

	   if(!tmpv0)continue;

	   Double_t rd=rndm->Uniform(0.,1.);

	   if(rd>fExtraEffPb) continue; 

	   Bool_t IsGoodV0 = ApplyV0Cuts(tmpv0, fLaType, kLambda, primVertex, fAODevent);
	;
	   if(!IsGoodV0)continue;

	   new ((*extraLa)[nAODLa++]) AliAODv0(*tmpv0); 
	   dummy = (*extraLa)[nAODLa-1];	
           fh1LaPt->Fill(tmpv0->Pt());
	   fh2LaEtaPhi->Fill(tmpv0->Eta(), tmpv0->Phi());
	   fListLa->Add(tmpv0);
	 }

	 for(Int_t it=0; it<nV0s; ++it){//for Antilambda candidates apply cuts, fill into inclusive ALa list
	   AliAODv0 *tmpv0 = 0x0;
 
	   if(fEmbedMode==kAODFull) tmpv0 = fAODevent->GetV0(it);
	   if(!tmpv0)continue;

	   Double_t rd=rndm->Uniform(0.,1.);
	   if(rd>fExtraEffPb) continue; 

	   Bool_t IsGoodV0 = ApplyV0Cuts(tmpv0, fALaType, kAntiLambda, primVertex, fAODevent);
	   if(!IsGoodV0)continue;

	   new ((*extraALa)[nAODALa++]) AliAODv0(*tmpv0); 
	   dummy = (*extraALa)[nAODALa-1];
	   fh1ALaPt->Fill(tmpv0->Pt());
	   fh2ALaEtaPhi->Fill(tmpv0->Eta(), tmpv0->Phi());
	   fListALa->Add(tmpv0);
	 }
	 
         if(fEmbedMode==kAODFull){
            nTracks = fAODevent->GetNumberOfTracks();
	    nV0s = fAODevent->GetNumberOfV0s();
	    

            for(Int_t ij=0; ij<nJets; ++ij){
	      AliAODJet *jet = dynamic_cast<AliAODJet*>(fAODJets->At(ij));
	      if(!jet) continue;
	      if(   ((fEvtSelMinJetPt==-1.) || (jet->Pt()>=fEvtSelMinJetPt))//min jet pt
		    && (fEvtSelMaxJetPt==-1. || jet->Pt()<=fEvtSelMaxJetPt)
		    && (jet->Eta()>=fEvtSelMinJetEta && jet->Eta()<=fEvtSelMaxJetEta)
		    && (jet->Phi()>=fEvtSelMinJetPhi && jet->Eta()<=fEvtSelMaxJetPhi)
		    && (jet->Pt()>=fEvtSelMinJetPt)//min jet pt cut
		    && HasMinLConstPt(jet)){
				  
		fh1JetPt->Fill(jet->Pt());
		fh2JetEtaPhi->Fill(jet->Eta(), jet->Phi());
		nSelectedJets++;
		Double_t radius = GetFFRadius();
		Double_t sumPtK0s = 0;

		if(fListK0s != 0){
		GetTracksInCone(fListK0s, fListK0sCone, jet, radius, sumPtK0s);

		if(fListK0sCone->GetSize() > 0){
		  for(Int_t it=0; it<fListK0sCone->GetSize(); ++it){ // loop K0s in jet cone
		    
		    AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListK0sCone->At(it));
		    if(!v0) continue;                    
		    
		    fh2K0sPtJetPtCone->Fill(jet->Pt(), v0->Pt());
		  }
		}
		
		if(fListK0sCone->GetSize() == 0){ 
		  fh2K0sPtJetPtCone->Fill(-1., -1.);
		  
		}
		}
		
		Double_t sumPtLa = 0;
		if(fListLa != 0){
		GetTracksInCone(fListLa, fListLaCone, jet, radius, sumPtLa);

		if(fListLaCone->GetSize() > 0){
		  for(Int_t it=0; it<fListLaCone->GetSize(); ++it){ // loop La in jet cone
		    
		    AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListLaCone->At(it));
		    if(!v0) continue;                    
		    
		    fh2LaPtJetPtCone->Fill(jet->Pt(), v0->Pt());
		  }}
		
		if(fListLaCone->GetSize() == 0){ 
		  fh2LaPtJetPtCone->Fill(-1., -1.);
		}
		}

		if(fListALa != 0){
		Double_t sumPtALa = 0;
		GetTracksInCone(fListALa, fListALaCone, jet, radius, sumPtALa);
		
		if(fListALaCone->GetSize() > 0){
		  for(Int_t it=0; it<fListALaCone->GetSize(); ++it){ // loop ALa in jet cone
		    
		    AliAODv0* v0 = dynamic_cast<AliAODv0*>(fListALaCone->At(it));
		    if(!v0) continue;                    
		    
		    fh2ALaPtJetPtCone->Fill(jet->Pt(), v0->Pt());
		  }
		}
		
		if(fListALaCone->GetSize() == 0){ 
		  fh2ALaPtJetPtCone->Fill(-1., -1.);
		}
		}
		
	      }
	      
	      if(fListK0sCone != 0){ 
		fListK0sCone->Clear();
	      }
	      
	      if(fListLaCone != 0){ 
		fListLaCone->Clear();
	      }
	      
	      if(fListALaCone != 0){ 
		fListALaCone->Clear();
	      }
            }				   
            fh1JetN->Fill(nSelectedJets);
         }
	 
         if(fEmbedMode==kAODJetTracks){
	   // look for leading jet within selection
	   // get reference tracks for this jet
	   Float_t leadJetPt = 0.;
	   for(Int_t ij=0; ij<nJets; ++ij){
               AliAODJet *jet = dynamic_cast<AliAODJet*>(fAODJets->At(ij));
               if(!jet) continue;
               if(   (fEvtSelMinJetPt==-1. || jet->Pt()>=fEvtSelMinJetPt)
                     && (fEvtSelMaxJetPt==-1. || jet->Pt()<=fEvtSelMaxJetPt)
                     && (jet->Eta()>=fEvtSelMinJetEta && jet->Eta()<=fEvtSelMaxJetEta)
                     && (jet->Phi()>=fEvtSelMinJetPhi && jet->Eta()<=fEvtSelMaxJetPhi)
		     && HasMinLConstPt(jet)){
		       
                  if(jet->Pt()>leadJetPt){
                     leadJetPt = jet->Pt();
                     leadJet = jet;
                  } 
               }
            }
            if(leadJet){
               nTracks = leadJet->GetRefTracks()->GetEntriesFast();
               fh1JetPt->Fill(leadJet->Pt());
               fh2JetEtaPhi->Fill(leadJet->Eta(), leadJet->Pt());
               fh1JetN->Fill(1);				   
            }
         }

         fh1TrackN->Fill((Float_t)nTracks);
	 
         for(Int_t it=0; it<nTracks; ++it){
	   AliAODTrack *tmpTr = 0x0;
	   if(fEmbedMode==kAODFull)      tmpTr = dynamic_cast<AliAODTrack*>(fAODevent->GetTrack(it));
	   if(fEmbedMode==kAODJetTracks) tmpTr = dynamic_cast<AliAODTrack*>(leadJet->GetRefTracks()->At(it));
	   if(!tmpTr) continue;
            Double_t rd=rndm->Uniform(0.,1.);
            if(rd>fExtraEffPb) continue; 
            tmpTr->SetStatus(AliESDtrack::kEmbedded);

            new ((*tracks)[nAODtracks++]) AliAODTrack(*tmpTr); 
            dummy = (*tracks)[nAODtracks-1];

            if(fTrackFilterMap<=0 || tmpTr->TestFilterBit(fTrackFilterMap)){
	      if(tmpTr->Pt()>0.15 && TMath::Abs(tmpTr->Eta())<0.9) fh1TrackPt->Fill(tmpTr->Pt());
	      if(tmpTr->Pt()>0.15 && TMath::Abs(tmpTr->Eta())<0.9) fh2TrackEtaPhi->Fill(tmpTr->Eta(), tmpTr->Phi());

           }
	 }



         if(mcpartIN){
            Int_t nMCpart = mcpartIN->GetEntriesFast();

            Int_t nPhysicalPrimary=0;
            Int_t nAODmcpart=0;
            for(Int_t ip=0; ip<nMCpart; ++ip){
               AliAODMCParticle *tmpPart = (AliAODMCParticle*) mcpartIN->At(ip);
	       if(!tmpPart) continue;
               if(fEmbedMode==kAODJetTracks){
                  // jet track? else continue
                  // not implemented yet
                  continue;
               } 

               if(tmpPart->IsPhysicalPrimary() && tmpPart->Charge()!=0. && tmpPart->Charge()!=-99.  && tmpPart->Pt()>0.){
		 new((*mcpartOUT)[nAODmcpart++]) AliAODMCParticle(*tmpPart);
		 dummy = (*mcpartOUT)[nAODmcpart-1];

		 if(fDebug>10) printf("added track %d with pT=%.2f to extra branch\n",nAODmcpart,tmpPart->Pt());
		 
		 fh1MCTrackPt->Fill(tmpPart->Pt());
		 fh2MCTrackEtaPhi->Fill(tmpPart->Eta(), tmpPart->Phi());
		 nPhysicalPrimary++;

               }
            }
            fh1MCTrackN->Fill((Float_t)nPhysicalPrimary);

         }
      } // end: embed all tracks or jet tracks




      if(fEmbedMode==kAODJet4Mom){

         // loop over jets
         Int_t nJets = fAODJets->GetEntries();
         fh1TrackN->Fill((Float_t)nJets);
         for(Int_t ij=0; ij<nJets; ++ij){
            AliAODJet *jet = dynamic_cast<AliAODJet*>(fAODJets->At(ij));
            if(!jet) continue;
            AliAODTrack *tmpTr = (AliAODTrack*)jet;
            tmpTr->SetFlags(AliESDtrack::kEmbedded);

            new ((*tracks)[nAODtracks++]) AliAODTrack(*tmpTr);
            dummy = (*tracks)[nAODtracks-1]; 

            fh1TrackPt->Fill(tmpTr->Pt());
            fh2TrackEtaPhi->Fill(tmpTr->Eta(), tmpTr->Phi());
         }

      } // end: embed jets as 4-momenta


   } //end: embed mode with AOD


   // === embed mode with toy ===
   if(fEmbedMode==kToyTracks){
      Int_t nT = (Int_t)(rndm->Uniform(fToyMinNbOfTracks, fToyMaxNbOfTracks)+0.5);

      fh1TrackN->Fill((Float_t)nT);

      for(Int_t i=0; i<nT; ++i){

         Double_t pt = -1.;
         if(fToyMinTrackPt!=fToyMaxTrackPt){
            if(fToyDistributionTrackPt==0){
               pt = rndm->Uniform(fToyMinTrackPt, fToyMaxTrackPt);
            } else {
               while(pt<fToyMinTrackPt||pt>fToyMaxTrackPt){
                  pt = rndm->Exp(fToyDistributionTrackPt);   // no power law yet!!
                  pt += fToyMinTrackPt;
               }
            }
         } else {
            pt = fToyMinTrackPt;
         }
         Double_t eta = rndm->Uniform(fToyMinTrackEta,fToyMaxTrackEta);
         Double_t phi = rndm->Uniform(fToyMinTrackPhi,fToyMaxTrackPhi);
         phi = TVector2::Phi_0_2pi(phi);

         if(fDebug) Printf("Add track #%d: pt %.2f, eta %.2f, phi %.2f", i, pt, eta, phi);

         Double_t theta = 2*TMath::ATan(1/TMath::Exp(eta));
         Float_t mom[3];
         mom[0] = pt;
         mom[1] = phi;
         mom[2] = theta;

         Float_t xyz[3];
         xyz[0] = -999.;
         xyz[1] = -999.;
         xyz[2] = -999.;

         AliAODTrack *tmpTr = new AliAODTrack( -999,   // id
         -999,   // label
         mom,    // momentum[3]
         kFALSE, // cartesian
         xyz,    // position
         kFALSE, // DCA
         NULL,   // covMatrix,
         -99,    // charge
         0,      // itsClusMap
	 //         NULL,   // pid 
         NULL,   // prodVertex
         kFALSE, // used for vertex fit
         kFALSE, // used for prim vtx fit
         AliAODTrack::kUndef, // type
         fToyFilterMap,  // select info
         -999.    // chi2 per NDF
         );
         tmpTr->SetFlags(AliESDtrack::kEmbedded);

         new((*tracks)[nAODtracks++]) AliAODTrack(*tmpTr);
         dummy = (*tracks)[nAODtracks-1];

         fh1TrackPt->Fill(pt);
         fh2TrackEtaPhi->Fill(eta,phi);

         delete tmpTr;
      }
   } //end: embed mode with toy



   fListK0s->Clear();
   fListLa->Clear();
   fListALa->Clear();
   fListK0sCone->Clear();
   fListLaCone->Clear();
   fListALaCone->Clear();
   
   

   PostData(1, fHistList);
}


//__________________________________________________________________________
void AliAnalysisTaskFastEmbedding::Terminate(Option_t *)
{
   // terminate
   if(fDebug > 1) Printf("AliAnalysisTaskFastEmbedding::Terminate()");

   if(fAODfile && fAODfile->IsOpen())
   fAODfile->Close();  

}

//__________________________________________________________________________
Int_t AliAnalysisTaskFastEmbedding::GetJobID()
{
   // gives back the alien subjob id, if available, else -1

   Int_t id=-1;

   const char* env = gSystem->Getenv("ALIEN_PROC_ID"); // GRID
   //   if(!env || !strlen(env)) env = gSystem->Getenv("LSB_JOBINDEX"); // GSI

   if(env && strlen(env)){
      id= atoi(env);
      AliInfo(Form("Job index %d", id));
   }
   else{
      AliInfo("Job index not found. Okay if running locally.");
   }

   return id;
}

//__________________________________________________________________________
Int_t AliAnalysisTaskFastEmbedding::SelectAODfile()
{
   // select an AOD file from fAODPathArray

   Int_t nFiles = fAODPathArray->GetEntries();

   // choose randomly file and event
   Int_t n = -1;
   Float_t tmpProp = -1.;
   Int_t pendingEvents = fInputEntries-Entry();
   //Printf("input entries %d, entry %d, pending events %d", fInputEntries, (Int_t)Entry(), pendingEvents);
   if(fAODEntriesArray){
     while(rndm->Rndm()>tmpProp){
       n = rndm->Integer(nFiles);
       fAODEntries = fAODEntriesArray->At(n);
       Int_t tmpEntries = fAODEntries<pendingEvents ? pendingEvents : fAODEntries;
       tmpProp = fAODEntriesMax ? (Float_t)tmpEntries/fAODEntriesMax : 1.;
     }
     fAODEntry = (Int_t)(rndm->Uniform(fAODEntries));
   }
   else {
      AliWarning("Number of entries in extra AODs not set!");
      n = rndm->Integer(nFiles);
      fAODEntry = 0;
   }

   AliInfo(Form("Select AOD file %d", n));
   TObjString *objStr = (TObjString*) fAODPathArray->At(n);
   if(!objStr){
      AliError("Could not get path of aod file.");
      return -1;
   }
   fAODPath = objStr->GetString();

   return n;
}

//__________________________________________________________________________
Int_t AliAnalysisTaskFastEmbedding::OpenAODfile(Int_t trial)
{
   // select and open an AOD file from fAODPathArray
  
   if(fAODPathArray) fFileId = SelectAODfile();
   if(fFileId<0) return -1;

   if (!gGrid) {
     AliInfo("Trying to connect to AliEn ...");
     TGrid::Connect("alien://");
   }

   TDirectory *owd = gDirectory;
   if (fAODfile && fAODfile->IsOpen()) fAODfile->Close();
   fAODfile = TFile::Open(fAODPath.Data(),"TIMEOUT=180");
   owd->cd();

   cout<<" Embedding task open AOD file "<<fAODPath.Data()<<" ID "<<fFileId<<endl;
   
   if(!fAODfile){

      fFileId = -1;
      if(fAODPathArray){
         if(trial<=3){ 
            AliError(Form("Could not open AOD file %s (trial %d), try again ...", fAODPath.Data(), trial));
            fFileId = OpenAODfile(trial+1);
         } else {
            AliError(Form("Could not open AOD file %s, give up ...",fAODPath.Data()));
            return -1;
         }
      } else {
         AliError(Form("Could not open AOD file %s.", fAODPath.Data()));
      }

      return fFileId;
   }

   fAODtree = (TTree*)fAODfile->Get("aodTree");

   if(!fAODtree){
      AliError("AOD tree not found.");
      return -1;
   }


   /*
      fAODtree->SetBranchStatus("*",0);
      fAODtree->SetBranchStatus("header",1);
      fAODtree->SetBranchStatus("tracks*",1);
      fAODtree->SetBranchStatus("jets*",1);
      fAODtree->SetBranchStatus("clusters*",1);
      fAODtree->SetBranchStatus("mcparticles*",1);
      fAODtree->SetBranchStatus("mcHeader*",1);
   */

   delete fAODevent;
   fAODevent = new AliAODEvent();
   fAODevent->ReadFromTree(fAODtree);


   // fetch header, jets, etc. from new file
   fNevents = fAODtree->GetEntries();
   mcHeader = (AliAODMCHeader*)fAODevent->FindListObject(AliAODMCHeader::StdBranchName());

   if(fJetBranch.Length()) fAODJets = dynamic_cast<TClonesArray*>(fAODevent->FindListObject(fJetBranch.Data()));
   else                    fAODJets = fAODevent->GetJets();
   if(!fAODJets){
      AliError("Could not find jets in AOD. Check jet branch when indicated.");
      return -1;
   }

   TFile *curfile = fAODtree->GetCurrentFile();
   if (!curfile) {
      AliError("No current file.");
      return -1;
   }

   Float_t trials = 1.;
   Float_t xsec = 0.;
   PythiaInfoFromFile(curfile->GetName(),xsec,trials);
   fXsection = xsec;

   // construct a poor man average trials 
   Float_t nEntries = (Float_t)fAODtree->GetTree()->GetEntries();
   if(trials>=nEntries && nEntries>0.)fAvgTrials = trials/nEntries;

   if(fFileId>=0){
      AliInfo(Form("Read successfully AOD event from file %d",fFileId));
   }

   fCountEvents=0; // new file, reset counter

   return fFileId;  // file position in AOD path array, if array available
}


//____________________________________________________________________________
Float_t AliAnalysisTaskFastEmbedding::GetPtHard(Bool_t bSet, Float_t newValue){

   // static stored, available for other tasks in train

   static Float_t ptHard = -1.;
   if(bSet) ptHard = newValue;

   return ptHard;
}

//____________________________________________________________________________
Int_t AliAnalysisTaskFastEmbedding::GetPtHardBin(Double_t ptHard){

   // returns pt hard bin (for LHC10e14, LHC11a1x, LHC11a2x simulations)

   const Int_t nBins = 10;
   Double_t binLimits[nBins] = { 5., 11., 21., 36., 57., 84., 117., 156., 200., 249. }; // lower limits

   Int_t bin = -1;
   while(bin<nBins-1 && binLimits[bin+1]<ptHard){
      bin++;
   }

   return bin;
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskFastEmbedding::PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials){
   //
   // get the cross section and the trails either from pyxsec.root or from pysec_hists.root
   // This should provide the path to the AOD/ESD file
   // (taken from PWG4/JetTasks/AliAnalysisHelperJetTasks) 

   TString file(currFile);  
   fXsec = 0;
   fTrials = 1;

   if(file.Contains("root_archive.zip#")){
     Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
     Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
     Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
      file.Replace(pos+1,pos2-pos1,"");
      //      file.Replace(pos+1,20,"");
   }
   else {
      // not an archive take the basename....
      file.ReplaceAll(gSystem->BaseName(file.Data()),"");
   }

   TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to live with open error message from root
   if(!fxsec){
      // next trial fetch the histgram file
      fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
      if(!fxsec){
         // not a severe condition but inciate that we have no information
         AliDebug(AliLog::kDebug,Form("Neither pyxsec.root nor pyxsec_hists.root found in %s",file.Data()));
         return kFALSE;
      }
      else{
         // find the tlist we want to be independtent of the name so use the Tkey
         TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0); 
         if(!key){
            fxsec->Close();
            return kFALSE;
         }
         TList *list = dynamic_cast<TList*>(key->ReadObj());
         if(!list){
            fxsec->Close();
            return kFALSE;
         }
         fXsec = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
         fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
         fxsec->Close();
      }
   } // no tree pyxsec.root
   else {
      TTree *xtree = (TTree*)fxsec->Get("Xsection");
      if(!xtree){
         fxsec->Close();
         return kFALSE;
      }
      UInt_t   ntrials  = 0;
      Double_t  xsection  = 0;
      xtree->SetBranchAddress("xsection",&xsection);
      xtree->SetBranchAddress("ntrials",&ntrials);
      xtree->GetEntry(0);
      fTrials = ntrials;
      fXsec = xsection;
      fxsec->Close();
   }
   return kTRUE;
}

//__________________________________________________________________________
Bool_t AliAnalysisTaskFastEmbedding::HasMinLConstPt(AliAODJet* jet){
  
  if(!(fEvtSelJetMinLConstPt>0)) return kTRUE;

  for(Int_t i=0; i<jet->GetRefTracks()->GetEntriesFast(); i++){
		    
    AliVParticle* track = dynamic_cast<AliVParticle*>(jet->GetRefTracks()->At(i));
    if(! track){
      continue; 
    }
    Double_t pt = track->Pt();		     
    if(pt>fEvtSelJetMinLConstPt) return kTRUE;
  }
 
  return kFALSE;
}


