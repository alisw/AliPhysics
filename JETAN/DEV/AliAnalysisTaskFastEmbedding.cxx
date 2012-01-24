/*************************************************************************
 *                                                                       *
 * Task for fast embedding                                               *
 * read extra input from AOD                                             *
 *                                                                       *
 *************************************************************************/


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
,fh1JetPt(0)
,fh2JetEtaPhi(0)
,fh1JetN(0)
,fh1MCTrackPt(0)
,fh2MCTrackEtaPhi(0)
,fh1MCTrackN(0)
,fh1AODfile(0)
,fh2AODevent(0)
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
,fh1JetPt(0)
,fh2JetEtaPhi(0)
,fh1JetN(0)
,fh1MCTrackPt(0)
,fh2MCTrackEtaPhi(0)
,fh1MCTrackN(0)
,fh1AODfile(0)
,fh2AODevent(0)
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
,fh1JetPt(copy.fh1JetPt)
,fh2JetEtaPhi(copy.fh2JetEtaPhi)
,fh1JetN(copy.fh1JetN)
,fh1MCTrackPt(copy.fh1MCTrackPt)
,fh2MCTrackEtaPhi(copy.fh2MCTrackEtaPhi)
,fh1MCTrackN(copy.fh1MCTrackN)
,fh1AODfile(copy.fh1AODfile)
,fh2AODevent(copy.fh2AODevent)
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
      fh1JetPt           = o.fh1JetPt;
      fh2JetEtaPhi       = o.fh2JetEtaPhi;
      fh1JetN            = o.fh1JetN;
      fh1MCTrackPt       = o.fh1MCTrackPt;
      fh2MCTrackEtaPhi   = o.fh2MCTrackEtaPhi;
      fh1MCTrackN        = o.fh1MCTrackN;
      fh1AODfile         = o.fh1AODfile;
      fh2AODevent        = o.fh2AODevent;
   }

   return *this;
}


//__________________________________________________________________________
AliAnalysisTaskFastEmbedding::~AliAnalysisTaskFastEmbedding()
{
   // destructor
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
   if(!fAODout->FindListObject(fTrackBranch.Data()) && strlen(fTrackBranch.Data())){
      AliInfo(Form("Add AOD branch %s", fTrackBranch.Data()));
      TClonesArray *tracks = new TClonesArray("AliAODTrack",0);
      tracks->SetName(fTrackBranch.Data());
      AddAODBranch("TClonesArray", &tracks);
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

   fh1TrackPt      =  new TH1F("fh1TrackPt","pT of extra tracks;p_{T};entries", 250, 0., 250.);
   fh2TrackEtaPhi  =  new TH2F("fh2TrackEtaPhi","eta-phi distribution of extra tracks;#eta;#phi", 20, -1., 1., 60, 0., 2*TMath::Pi());
   fh1TrackN       =  new TH1F("fh1TrackN", "nb. of extra tracks per event;nb. of tracks;entries",300, 0., 300.);

   fHistList->Add(fh1TrackPt);
   fHistList->Add(fh2TrackEtaPhi);
   fHistList->Add(fh1TrackN);
   
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
   
   fh1AODfile = new TH1I("fh1AODfile", "overview of opened AOD files from the array", 2300, -0.5, 2299.5);
   fh2AODevent = new TH2I("fh1AODevent","selected events;file;event", 2500,-.5,2499.5,5000,-.5,4999.5);
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
Bool_t AliAnalysisTaskFastEmbedding::UserNotify()
{

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

   // vertex selection
   AliAODVertex* primVtx = fAODout->GetPrimaryVertex();
   Int_t nTracksPrim = primVtx->GetNContributors();
   if ((nTracksPrim < fMinContribVtx) ||
         (primVtx->GetZ() < fVtxZMin) ||
         (primVtx->GetZ() > fVtxZMax) ){
      if(fDebug) Printf("%s:%d primary vertex z = %f: event REJECTED...",(char*)__FILE__,__LINE__,primVtx->GetZ());
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

   TRef dummy;

   // === embed mode with AOD ===
   if(fEmbedMode==kAODFull || fEmbedMode==kAODJetTracks || fEmbedMode==kAODJet4Mom){
      if(!fAODevent){
         AliError("Need input AOD, but is not connected."); 
         PostData(1, fHistList);
         return;
      }



      Bool_t useEntry = kFALSE;
      while(!useEntry){  // protection need, if no event fulfills requierment

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

               if(   (fEvtSelMinJetPt==-1. || jet->Pt()>=fEvtSelMinJetPt)
                     && (fEvtSelMaxJetPt==-1. || jet->Pt()<=fEvtSelMaxJetPt)
                     && (jet->Eta()>=fEvtSelMinJetEta && jet->Eta()<=fEvtSelMaxJetEta)
                     && (jet->Phi()>=fEvtSelMinJetPhi && jet->Phi()<=fEvtSelMaxJetPhi)){
                  useEntry = kTRUE;
                  break;
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
         Int_t nJets = fAODJets->GetEntries();
         Int_t nSelectedJets = 0;
         AliAODJet *leadJet = 0x0; // used for jet tracks only

         if(fEmbedMode==kAODFull){
            nTracks = fAODevent->GetNumberOfTracks();

            for(Int_t ij=0; ij<nJets; ++ij){
               AliAODJet *jet = dynamic_cast<AliAODJet*>(fAODJets->At(ij));
               if(!jet) continue;
               if(   (fEvtSelMinJetPt==-1. || jet->Pt()>=fEvtSelMinJetPt)
                     && (fEvtSelMaxJetPt==-1. || jet->Pt()<=fEvtSelMaxJetPt)
                     && (jet->Eta()>=fEvtSelMinJetEta && jet->Eta()<=fEvtSelMaxJetEta)
                     && (jet->Phi()>=fEvtSelMinJetPhi && jet->Eta()<=fEvtSelMaxJetPhi)){

                  fh1JetPt->Fill(jet->Pt());
                  fh2JetEtaPhi->Fill(jet->Eta(), jet->Phi());
                  nSelectedJets++;

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
                     && (jet->Phi()>=fEvtSelMinJetPhi && jet->Eta()<=fEvtSelMaxJetPhi)){
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
            if(fEmbedMode==kAODFull)      tmpTr = fAODevent->GetTrack(it);
            if(fEmbedMode==kAODJetTracks) tmpTr = dynamic_cast<AliAODTrack*>(leadJet->GetRefTracks()->At(it));
            if(!tmpTr) continue;

            tmpTr->SetStatus(AliESDtrack::kEmbedded);

            new ((*tracks)[nAODtracks++]) AliAODTrack(*tmpTr); 
            dummy = (*tracks)[nAODtracks-1];

            if(fTrackFilterMap<=0 || tmpTr->TestFilterBit(fTrackFilterMap)){
               fh1TrackPt->Fill(tmpTr->Pt());
               fh2TrackEtaPhi->Fill(tmpTr->Eta(), tmpTr->Phi());
            }
         }

         if(mcpartIN){
            Int_t nMCpart = mcpartIN->GetEntriesFast();

            Int_t nPhysicalPrimary=0;
            Int_t nAODmcpart=0;
            for(Int_t ip=0; ip<nMCpart; ++ip){
               AliAODMCParticle *tmpPart = (AliAODMCParticle*) mcpartIN->At(ip);

               if(fEmbedMode==kAODJetTracks){
                  // jet track? else continue
                  // not implemented yet
                  continue;
               } 

               new((*mcpartOUT)[nAODmcpart++]) AliAODMCParticle(*tmpPart);
               dummy = (*mcpartOUT)[nAODmcpart-1];

               if(tmpPart->IsPhysicalPrimary() && tmpPart->Charge()!=0. && tmpPart->Charge()!=-99. ){
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
         NULL,   // pid 
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


   PostData(1, fHistList);
}


//__________________________________________________________________________
void AliAnalysisTaskFastEmbedding::Terminate(Option_t *)
{
   if(fDebug > 1) Printf("AliAnalysisTaskFastEmbedding::Terminate()");

   if(fAODfile && fAODfile->IsOpen())
   fAODfile->Close();  

}

//__________________________________________________________________________
Int_t AliAnalysisTaskFastEmbedding::GetJobID()
{
   Int_t id=-1;

   const char* env = gSystem->Getenv("ALIEN_PROC_ID"); // GRID
   //if(!env || !strlen(env)) env = gSystem->Getenv("LSB_JOBINDEX"); // GSI

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

Int_t AliAnalysisTaskFastEmbedding::SelectAODfile(){

   Int_t nFiles = fAODPathArray->GetEntries();

   // choose randomly file and event
   Int_t n = -1;
   Float_t tmpProp = -1.;
   Int_t pendingEvents = fInputEntries-Entry();
   //Printf("input entries %d, entry %d, pending events %d", fInputEntries, (Int_t)Entry(), pendingEvents);
   while(rndm->Rndm()>tmpProp){
      n = rndm->Integer(nFiles);
      fAODEntries = fAODEntriesArray->At(n);
      Int_t tmpEntries = fAODEntries<pendingEvents ? pendingEvents : fAODEntries;
      tmpProp = fAODEntriesMax ? (Float_t)tmpEntries/fAODEntriesMax : 1.;
   }
   fAODEntry = (Int_t)(rndm->Uniform(fAODEntries));

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

Int_t AliAnalysisTaskFastEmbedding::OpenAODfile(Int_t trial){

   if(fAODPathArray) fFileId = SelectAODfile();
   if(fFileId<0) return -1;

   TDirectory *owd = gDirectory;
   if (fAODfile) fAODfile->Close();
   fAODfile = TFile::Open(fAODPath.Data(),"TIMEOUT=180");
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

   static Float_t ptHard = -1.;
   if(bSet) ptHard = newValue;

   return ptHard;
}

//____________________________________________________________________________
Int_t AliAnalysisTaskFastEmbedding::GetPtHardBin(Double_t ptHard){

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
      Ssiz_t pos1 = file.Index("root_archive",12,TString::kExact);
      Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
      file.Replace(pos+1,20,"");
   }
   else {
      // not an archive take the basename....
      file.ReplaceAll(gSystem->BaseName(file.Data()),"");
   }

   TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
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


