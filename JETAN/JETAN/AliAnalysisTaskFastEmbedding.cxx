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

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskFastEmbedding)

//__________________________________________________________________________
AliAnalysisTaskFastEmbedding::AliAnalysisTaskFastEmbedding()
: AliAnalysisTaskSE()
  ,fESD(0)
  ,fQAMode(0)
  ,fAODout(0)
  ,fAODevent(0)
  ,fAODeventJets(0)
  ,fAODtree(0)
  ,fAODtreeJets(0)
  ,fAODfile(0)
  ,fAODfileJets(0)
  ,mcHeader(0)
  ,fFFRadius(0)
  ,fK0Type(0) 
  ,fLaType(0) 
  ,fALaType(0)
  ,fListK0s(0)
  ,fListLa(0)
  ,fListALa(0)
  /*  ,fListK0sCone(0)
  ,fListLaCone(0)
  ,fListALaCone(0)*/
  ,fPIDResponse(0)
  ,rndm(0)
  ,fInputEntries(0)
  ,fAODPathArray(0)
  ,fAODEntriesArray(0)
  ,fAODPath("AliAOD.root")
  ,fStrJetFriends("")
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
  ,fEPMode(0)
  ,fEvtSelecMode(0)
  ,fEvtSelMinJetPt(-1)
  ,fEvtSelMaxJetPt(-1)
  ,fEvtSelMinJetEta(-999.)
  ,fEvtSelMaxJetEta( 999.)
  ,fEvtSelMinJetPhi(0.)
  ,fEvtSelMaxJetPhi(TMath::Pi()*2.)
  ,fEvtSelJetMinLConstPt(0)
  ,fExtraEffPb(1)
  ,fDiceMapEff(0)
  ,fPathDiceTrackMap("")
  ,fhEffH1(0)
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
  //,fh2K0sPtJetPtCone(0)
  ,fh1LaPt(0)
  ,fh2LaEtaPhi(0)
  //,fh2LaPtJetPtCone(0)
  ,fh1ALaPt(0)
  ,fh2ALaEtaPhi(0)
  //,fh2ALaPtJetPtCone(0)
  ,fh1JetPt(0)
  ,fh2JetEtaPhi(0)
  ,fh1JetN(0)
  ,fh1MCTrackPt(0)
  ,fh2MCTrackEtaPhi(0)
  ,fh1MCTrackN(0)
  ,fh1AODfile(0)
  ,fh2AODevent(0)       
  ,fh1EP2(0)
  ,fh1EP3(0)

{
  // default constructor
  
}


//__________________________________________________________________________
AliAnalysisTaskFastEmbedding::AliAnalysisTaskFastEmbedding(const char *name)
: AliAnalysisTaskSE(name)
,fESD(0)
,fQAMode(0)
,fAODout(0)
,fAODevent(0)
,fAODeventJets(0)
,fAODtree(0)
,fAODtreeJets(0)
,fAODfile(0)
,fAODfileJets(0)
,mcHeader(0)
,fFFRadius(0)
,fK0Type(0) 
,fLaType(0) 
,fALaType(0)
,fListK0s(0)
,fListLa(0)
,fListALa(0) 
  /*,fListK0sCone(0)
,fListLaCone(0)
,fListALaCone(0)*/ 
,fPIDResponse(0)
,rndm(0)
,fInputEntries(0)
,fAODPathArray(0)
,fAODEntriesArray(0)
,fAODPath("AliAOD.root")
,fStrJetFriends("")
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
,fEPMode(0)
,fEvtSelecMode(0)
,fEvtSelMinJetPt(-1)
,fEvtSelMaxJetPt(-1)
,fEvtSelMinJetEta(-999.)
,fEvtSelMaxJetEta( 999.)
,fEvtSelMinJetPhi(0.)
,fEvtSelMaxJetPhi(TMath::Pi()*2.)
,fEvtSelJetMinLConstPt(0)
,fExtraEffPb(1)
,fDiceMapEff(0)
,fPathDiceTrackMap("")
,fhEffH1(0)
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
  //,fh2K0sPtJetPtCone(0)
,fh1LaPt(0)
,fh2LaEtaPhi(0)
  //,fh2LaPtJetPtCone(0)
,fh1ALaPt(0)
,fh2ALaEtaPhi(0)
  //,fh2ALaPtJetPtCone(0)
,fh1JetPt(0)
,fh2JetEtaPhi(0)
,fh1JetN(0)
,fh1MCTrackPt(0)
,fh2MCTrackEtaPhi(0)
,fh1MCTrackN(0)
,fh1AODfile(0)
,fh2AODevent(0)         
,fh1EP2(0)
,fh1EP3(0)

{
   // constructor
   DefineOutput(1, TList::Class());
}


//__________________________________________________________________________
AliAnalysisTaskFastEmbedding::AliAnalysisTaskFastEmbedding(const AliAnalysisTaskFastEmbedding &copy)
: AliAnalysisTaskSE()
,fESD(copy.fESD)
,fQAMode(copy.fQAMode)
,fAODout(copy.fAODout)
,fAODevent(copy.fAODevent)
,fAODeventJets(copy.fAODeventJets)
,fAODtree(copy.fAODtree)
,fAODtreeJets(copy.fAODtreeJets)
,fAODfile(copy.fAODfile)
,fAODfileJets(copy.fAODfileJets)
,mcHeader(copy.mcHeader)
,fFFRadius(copy.fFFRadius)
,fK0Type(copy.fK0Type)
,fLaType(copy.fLaType)          
,fALaType(copy.fALaType)   
,fListK0s(copy.fListK0s)
,fListLa(copy.fListLa)
,fListALa(copy.fListALa)
  /*,fListK0sCone(copy.fListK0sCone)
,fListLaCone(copy.fListLaCone)
,fListALaCone(copy.fListALaCone)*/
,fPIDResponse(copy.fPIDResponse)
,rndm(copy.rndm)
,fInputEntries(copy.fInputEntries)
,fAODPathArray(copy.fAODPathArray)
,fAODEntriesArray(copy.fAODEntriesArray)
,fAODPath(copy.fAODPath)
,fStrJetFriends(copy.fStrJetFriends)
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
,fEPMode(copy.fEPMode)
,fEvtSelecMode(copy.fEvtSelecMode)
,fEvtSelMinJetPt(copy.fEvtSelMinJetPt)
,fEvtSelMaxJetPt(copy.fEvtSelMaxJetPt)
,fEvtSelMinJetEta(copy.fEvtSelMinJetEta)
,fEvtSelMaxJetEta(copy.fEvtSelMaxJetEta)
,fEvtSelMinJetPhi(copy.fEvtSelMinJetPhi)
,fEvtSelMaxJetPhi(copy.fEvtSelMaxJetPhi)
,fEvtSelJetMinLConstPt(copy.fEvtSelJetMinLConstPt)
,fExtraEffPb(copy.fExtraEffPb)
,fDiceMapEff(copy.fDiceMapEff)
,fPathDiceTrackMap(copy.fPathDiceTrackMap)
,fhEffH1(copy.fhEffH1)
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
  //,fh2K0sPtJetPtCone(copy.fh2K0sPtJetPtCone)
,fh1LaPt(copy.fh1LaPt)
,fh2LaEtaPhi(copy.fh2LaEtaPhi)
  //,fh2LaPtJetPtCone(copy.fh2LaPtJetPtCone)
,fh1ALaPt(copy.fh1ALaPt)
,fh2ALaEtaPhi(copy.fh2ALaEtaPhi)
  //,fh2ALaPtJetPtCone(copy.fh2ALaPtJetPtCone)
,fh1JetPt(copy.fh1JetPt)
,fh2JetEtaPhi(copy.fh2JetEtaPhi)
,fh1JetN(copy.fh1JetN)
,fh1MCTrackPt(copy.fh1MCTrackPt)
,fh2MCTrackEtaPhi(copy.fh2MCTrackEtaPhi)
,fh1MCTrackN(copy.fh1MCTrackN)
,fh1AODfile(copy.fh1AODfile)
,fh2AODevent(copy.fh2AODevent)
,fh1EP2(copy.fh1EP2)
,fh1EP3(copy.fh1EP3)


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
      fQAMode            = o.fQAMode;
      fAODout            = o.fAODout;
      fAODevent          = o.fAODevent;
      fAODeventJets      = o.fAODeventJets;
      fAODtree           = o.fAODtree;
      fAODtreeJets       = o.fAODtreeJets;
      fStrJetFriends     = o.fStrJetFriends; 
      fAODfile           = o.fAODfile;
      fAODfileJets       = o.fAODfileJets;
      mcHeader           = o.mcHeader;
      fFFRadius          = o.fFFRadius;
      fK0Type            = o.fK0Type;
      fLaType            = o.fLaType;
      fALaType           = o.fALaType;
      fListK0s           = o.fListK0s;
      fListLa            = o.fListLa;
      fListALa           = o.fListALa;
      /*fListK0sCone       = o.fListK0sCone;
      fListLaCone        = o.fListLaCone;
      fListALaCone       = o.fListALaCone;*/
      fPIDResponse       = o.fPIDResponse;
      rndm               = o.rndm;
      fInputEntries      = o.fInputEntries;
      fAODPathArray      = o.fAODPathArray;
      fAODEntriesArray   = o.fAODEntriesArray;
      fAODPath           = o.fAODPath;
      fStrJetFriends     = o.fStrJetFriends; 
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
      fEPMode            = o.fEPMode;
      fEvtSelecMode      = o.fEvtSelecMode;
      fEvtSelMinJetPt    = o.fEvtSelMinJetPt;
      fEvtSelMaxJetPt    = o.fEvtSelMaxJetPt;
      fEvtSelMinJetEta   = o.fEvtSelMinJetEta;
      fEvtSelMaxJetEta   = o.fEvtSelMaxJetEta;
      fEvtSelMinJetPhi   = o.fEvtSelMinJetPhi;
      fEvtSelMaxJetPhi   = o.fEvtSelMaxJetPhi;
      fEvtSelJetMinLConstPt = o.fEvtSelJetMinLConstPt;
      fExtraEffPb        = o.fExtraEffPb;
      fDiceMapEff          = o.fDiceMapEff;
      fPathDiceTrackMap       =o.fPathDiceTrackMap;
      fhEffH1                   = o.fhEffH1;
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
      //fh2K0sPtJetPtCone  = o.fh2K0sPtJetPtCone;
      fh1LaPt            = o.fh1LaPt;
      fh2LaEtaPhi        = o.fh2LaEtaPhi;
      //fh2LaPtJetPtCone   = o.fh2LaPtJetPtCone;
      fh1ALaPt           = o.fh1ALaPt;
      fh2ALaEtaPhi       = o.fh2ALaEtaPhi;
      //fh2ALaPtJetPtCone  = o.fh2ALaPtJetPtCone;
      fh1JetPt           = o.fh1JetPt;
      fh2JetEtaPhi       = o.fh2JetEtaPhi;
      fh1JetN            = o.fh1JetN;
      fh1MCTrackPt       = o.fh1MCTrackPt;
      fh2MCTrackEtaPhi   = o.fh2MCTrackEtaPhi;
      fh1MCTrackN        = o.fh1MCTrackN;
      fh1AODfile         = o.fh1AODfile;
      fh2AODevent        = o.fh2AODevent;
      fh1EP2             = o.fh1EP2;
      fh1EP3             = o.fh1EP3;
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
  /* if(fListK0sCone) delete fListK0sCone;
  if(fListLaCone) delete fListLaCone;
  if(fListALaCone) delete fListALaCone;*/
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

   /* fListK0sCone = new TList();
   fListK0sCone->SetOwner(kFALSE);

   fListLaCone = new TList();
   fListLaCone->SetOwner(kFALSE);

   fListALaCone = new TList();
   fListALaCone->SetOwner(kFALSE);*/


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

  /*  //for the v0 particle candidates inside the jets (only in kAODJetTracks mode)
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
    }*/
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
   //fh2K0sPtJetPtCone  =  new TH2F("fh2K0sPtJetPtCone","v0-pt jet pt distribution of extra K0s;#it{p_{T}^{jet,ch}}};#it{p_{T}}", 20, 0., 100.,120, 0., 12.);
   fh1LaPt         =  new TH1F("fh1LaPt","pT of extra associated La;p_{T};entries", 120, 0., 12.);
   fh2LaEtaPhi     =  new TH2F("fh2LaEtaPhi","eta-phi distribution of extra associated #Lamdba;#eta;#phi", 20, -1., 1., 60, 0., 2*TMath::Pi());
   //fh2LaPtJetPtCone  =  new TH2F("fh2LaPtJetPtCone","v0-pt jet pt distribution of extra associated #Lamdba;#it{p_{T}^{jet,ch}}};#it{p_{T}}", 20, 0., 100.,120, 0., 12.);
   fh1ALaPt        =  new TH1F("fh1ALaPt","pT of extra associated #bar{#Lamdba};p_{T};entries", 120, 0., 12.);
   fh2ALaEtaPhi    =  new TH2F("fh2ALaEtaPhi","eta-phi distribution of extra associated #bar{#Lamdba};#eta;#phi", 20, -1., 1., 60, 0., 2*TMath::Pi());
   //fh2ALaPtJetPtCone  =  new TH2F("fh2ALaPtJetPtCone","v0-pt jet pt distribution of extra associated #bar{#Lamdba};#it{p_{T}^{jet,ch}}};#it{p_{T}}", 20, 0., 100.,120, 0., 12.);
   fh1EP2          = new TH1F("fh1EP2","2nd order harmonic event plane distribution ESD events",200,-2.,2.);
   fh1EP3          = new TH1F("fh1EP3","3rd order harmonic event plane distribution ESD events",200,-2.,2.);
  
   fHistList->Add(fh1TrackPt);
   fHistList->Add(fh2TrackEtaPhi);
   fHistList->Add(fh1TrackN);
   fHistList->Add(fh1V0Pt);
   fHistList->Add(fh2V0EtaPhi);
   fHistList->Add(fh1K0Pt);
   fHistList->Add(fh2K0EtaPhi);
   //fHistList->Add(fh2K0sPtJetPtCone);
   fHistList->Add(fh1LaPt);
   fHistList->Add(fh2LaEtaPhi);
   //fHistList->Add(fh2LaPtJetPtCone);
   fHistList->Add(fh1ALaPt);
   fHistList->Add(fh2ALaEtaPhi);
   //fHistList->Add(fh2ALaPtJetPtCone);
   if(fEPMode != 0){
   fHistList->Add(fh1EP2);
   fHistList->Add(fh1EP3);}
  
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
   
   fh1AODfile = new TH1I("fh1AODfile", "overview of opened AOD files from the array", 250, -0.5, 249.5);
   fh2AODevent = new TH2I("fh2AODevent","selected events;file;event", 500,-0.5,499.5,500,-0.5,499.5);
  
   if(fQAMode){ 
   fHistList->Add(fh1AODfile);
   fHistList->Add(fh2AODevent);
   }
   
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
   if(fDiceMapEff==kTRUE) LoadDiceTrackMapRootFile();
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

       Int_t nnum = -1;
       Int_t pnum = -1;
       
       //get MC stack

       TClonesArray *stackmc = 0x0;
       stackmc = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName())); //get MCAOD branch in data
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
       
       TList *listmc = fAOD->GetList(); //AliAODEvent* is inherited from AliVEvent*, listmc is pointer to reconstructed event in MC list, member of AliAODEvent
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

       return kTRUE;

}
/*	
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
  
  }*/

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
 
   // vertex selection, get ESD (PbPb) vertex
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

   //Get event plane

   Double_t qx2(0), qy2(0), qx3(0), qy3(0);
  
   // get the q-vectors from AliEventPlane
   Double_t event_plane_2 = fESD->GetEventplane()->CalculateVZEROEventPlane(fESD, 10, 2, qx2, qy2); // returns V0 2nd order EP (EP phi)
   Double_t event_plane_3 = fESD->GetEventplane()->CalculateVZEROEventPlane(fESD, 10,  3, qx3, qy3);  // return V0 3rd order EP (EP phi)

   if(fDebug>3) std::cout<<"2nd order harmonic EP - qx2: "<<qx2<<" , qy2: "<<qy2<<std::endl;
   if(fDebug>3) std::cout<<"3rd order harmonic EP - qx3: "<<qx3<<" , qy3: "<<qy3<<std::endl;
   if(fDebug>3) std::cout<<"event_plane_2: "<<event_plane_2<<" , event_plane_3: "<<event_plane_3<<std::endl;
 
   //std::cout<<"PbPb vertex z coordinate: "<<primVertex->GetZ()<<std::endl;

   fh1EP2->Fill(event_plane_2);
   fh1EP3->Fill(event_plane_3);

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
   /*
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
   */
   extrav0s->Delete();
   extraK0s->Delete();
   extraLa->Delete();
   extraALa->Delete();
   /* extraK0sCone->Delete();
   extraLaCone->Delete();
   extraALaCone->Delete();
   */
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
	       if(fQAMode){
               fh1AODfile->Fill(fFileId);}
               if(fAODEntries!=fNevents){
                  AliError("File with list of AODs and file with nb. of entries inconsistent");
                  PostData(1, fHistList);
                  return;
               }
            }
         }

         // get next event
         fAODtree->GetEvent(fAODEntry);
         if(fAODtreeJets) fAODtreeJets->GetEvent(fAODEntry);


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

	 //AliAODVertex* pythiaVertex = fAODevent->GetPrimaryVertex();

	 //std::cout<<" pythia vertex z coordinate: "<<pythiaVertex->GetZ()<<std::endl;

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
      if(fQAMode){
      fh2AODevent->Fill(fFileId,fAODEntry);
      }

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
	   AliAODVertex* pythiaVertex = fAODevent->GetPrimaryVertex();


	   IsGoodV0 = ApplyV0Cuts(tmpv0, fK0Type, kK0, pythiaVertex, fAODevent);
          
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
	   AliAODVertex* pythiaVertex = fAODevent->GetPrimaryVertex();

	   Bool_t IsGoodV0 = ApplyV0Cuts(tmpv0, fLaType, kLambda, pythiaVertex, fAODevent);
	
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
	   AliAODVertex* pythiaVertex = fAODevent->GetPrimaryVertex();

	   Bool_t IsGoodV0 = ApplyV0Cuts(tmpv0, fALaType, kAntiLambda, pythiaVertex, fAODevent);
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
		
		Double_t jetPhi = jet->Phi();
		
		fh1JetPt->Fill(jet->Pt());
		fh2JetEtaPhi->Fill(jet->Eta(), jetPhi);
		
		
		if(fEPMode != 0){  //do embedding event-plane-dependent
		  
		  Double_t jetPhiDelta = 0;
		  
		  if(fEPMode == 1){//rotate jet axis in-plane
		    
		    //new (in-plane-rotated) jet phi fulfills the following condition
		    jetPhiDelta = jetPhi-event_plane_2;//absolute difference of event plane and PYTHIA jet phi
		    
		    if(fDebug>3)std::cout<<"In-plane jetPhiDelta: "<<jetPhiDelta<<std::endl;
		    
		  }
		  
		  if(fEPMode == 2){//rotate jet axis out-of-plane
		    //new (out-of-plane-rotated) jet phi fulfills the following condition
		    jetPhiDelta = (jetPhi-(event_plane_2))+(TMath::Pi()*0.5);
		    
		    if(fDebug>3)std::cout<<"Out-of-plane jetPhiDelta: "<<jetPhiDelta<<std::endl;
		    
		    
		  }
		  
		  Int_t nJetTracks = jet->GetRefTracks()->GetEntriesFast();
		  
		  for(Int_t it=0; it<nJetTracks; ++it){
		    AliAODTrack *tmpTr = 0x0;
		    tmpTr = dynamic_cast<AliAODTrack*>(jet->GetRefTracks()->At(it));
		    if(!tmpTr) continue;
		    Double_t trackphi = tmpTr->Phi(); 
		    
		    if(fDebug>3)std::cout<<"jet track phi before rotation: "<<trackphi<<std::endl;
		    
		    trackphi=trackphi-jetPhiDelta;//rotate jet constituent in-/out-of-plane
		    tmpTr->SetPhi(trackphi);//save new phi angle of jet track
		    
		    if(fDebug>3)std::cout<<"jet track phi after rotation: "<<trackphi<<std::endl;
		    
		  }

		}
	      }
	      
	      nSelectedJets++;
	      
	      
		
	      
	      /* Double_t radius = GetFFRadius();
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
		}*/
	      
	      
	      
	      /* if(fListK0sCone != 0){ 
		 fListK0sCone->Clear();
		 }
		 
		 if(fListLaCone != 0){ 
		 fListLaCone->Clear();
		 }
		 
		 if(fListALaCone != 0){ 
		 fListALaCone->Clear();
		 }*/
	      
	      fh1JetN->Fill(nSelectedJets);
	    }
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
             if(fDiceMapEff==kTRUE){
	      Double_t pTtmp=tmpTr->Pt();
	      if(pTtmp>100) pTtmp=100;
              Double_t efffrac=fhEffH1->GetBinContent(fhEffH1->FindBin(pTtmp));
              if(rd>efffrac) continue;}


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

               if(tmpPart->IsPhysicalPrimary() && tmpPart->Charge()!=-99.  && tmpPart->Pt()>0.){
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
   /*  fListK0sCone->Clear();
   fListLaCone->Clear();
   fListALaCone->Clear();
   */  
   

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
   
   if(fStrJetFriends.Length()){
     TString fAODPathJets = fAODPath; 
     fAODPathJets.ReplaceAll("AliAOD.root",fStrJetFriends.Data()); 

     if (fAODfileJets && fAODfileJets->IsOpen()) fAODfileJets->Close();
     fAODfileJets = TFile::Open(fAODPathJets.Data(),"TIMEOUT=180");

     cout<<" Embedding task open AOD jets file "<<fAODPathJets.Data()<<endl;
   }
   owd->cd();

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
   if(fAODfileJets) fAODtreeJets = (TTree*)fAODfileJets->Get("aodTree");



   if(!fAODtree){
      AliError("AOD tree not found.");
      return -1;
   }

  
   if(fStrJetFriends.Length() && !fAODtreeJets){
      AliError("AOD tree jets not found.");
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

   delete fAODeventJets;
   fAODeventJets = new AliAODEvent();
   fAODeventJets->ReadFromTree(fAODtreeJets);


   // fetch header, jets, etc. from new file
   fNevents = fAODtree->GetEntries();
   mcHeader = (AliAODMCHeader*)fAODevent->FindListObject(AliAODMCHeader::StdBranchName());
   
   
   if(fJetBranch.Length()){ 
     fAODJets = dynamic_cast<TClonesArray*>(fAODevent->FindListObject(fJetBranch.Data()));
     if(!fAODJets){
       fAODJets = dynamic_cast<TClonesArray*>(fAODeventJets->FindListObject(fJetBranch.Data()));
     }
   }
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

void AliAnalysisTaskFastEmbedding::LoadDiceTrackMapRootFile() {

   if (!gGrid) {
     AliInfo("Trying to connect to AliEn ...");
     TGrid::Connect("alien://");
   }

  TFile *f = TFile::Open(fPathDiceTrackMap.Data());
  if(!f)return;
   TH1D *hEffMap = (TH1D*)f->Get("hEffPosGlobSt");

   SetEfficiencyMap(hEffMap);

}

void AliAnalysisTaskFastEmbedding:: SetEfficiencyMap(TH1 *h1) {
  fhEffH1 = (TH1*)h1->Clone("fhEffH1");
}
