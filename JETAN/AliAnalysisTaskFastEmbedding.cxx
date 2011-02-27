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


#include "AliAnalysisTaskFastEmbedding.h"
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODJet.h"
#include "AliAODMCParticle.h"

#include "AliLog.h"

ClassImp(AliAnalysisTaskFastEmbedding)

//__________________________________________________________________________
AliAnalysisTaskFastEmbedding::AliAnalysisTaskFastEmbedding()
    : AliAnalysisTaskSE()
      ,fAODout(0)
      ,fAODevent(0)
      ,fAODtree(0)
      ,fAODfile(0)
      ,rndm(0)
      ,fAODPathArray(0)
      ,fAODPath("AliAOD.root")
      ,fTrackBranch("aodExtraTracks")
      ,fMCparticlesBranch("aodExtraMCparticles")
      ,fJetBranch("")
      ,fEntry(0)
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
      ,fHistList(0)
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


{
    // default constructor

}


//__________________________________________________________________________
AliAnalysisTaskFastEmbedding::AliAnalysisTaskFastEmbedding(const char *name)
    : AliAnalysisTaskSE(name)
      ,fAODout(0)
      ,fAODevent(0)
      ,fAODtree(0)
      ,fAODfile(0)
      ,rndm(0)
      ,fAODPathArray(0)
      ,fAODPath("AliAOD.root")
      ,fTrackBranch("aodExtraTracks")
      ,fMCparticlesBranch("aodExtraMCparticles")
      ,fJetBranch("")
      ,fEntry(0)
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
      ,fHistList(0)
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
{
    // constructor
    DefineOutput(1, TList::Class());
}


//__________________________________________________________________________
AliAnalysisTaskFastEmbedding::AliAnalysisTaskFastEmbedding(const AliAnalysisTaskFastEmbedding &copy)
    : AliAnalysisTaskSE()
      ,fAODout(copy.fAODout)
      ,fAODevent(copy.fAODevent)
      ,fAODtree(copy.fAODtree)
      ,fAODfile(copy.fAODfile)
      ,rndm(copy.rndm)
      ,fAODPathArray(copy.fAODPathArray)
      ,fAODPath(copy.fAODPath)
      ,fTrackBranch(copy.fTrackBranch)
      ,fMCparticlesBranch(copy.fMCparticlesBranch)
      ,fJetBranch(copy.fJetBranch)
      ,fEntry(copy.fEntry)
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
      ,fHistList(copy.fHistList)
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
{
    // copy constructor
}


//__________________________________________________________________________
AliAnalysisTaskFastEmbedding& AliAnalysisTaskFastEmbedding::operator=(const AliAnalysisTaskFastEmbedding& o)
{
    // assignment

    if(this!=&o){
	AliAnalysisTaskSE::operator=(o);
	fAODout            = o.fAODout;
	fAODevent          = o.fAODevent;
	fAODtree           = o.fAODtree;
	fAODfile           = o.fAODfile;
        rndm               = o.rndm;
        fAODPathArray       = o.fAODPathArray;
        fAODPath           = o.fAODPath;
        fTrackBranch       = o.fTrackBranch;
        fMCparticlesBranch = o.fMCparticlesBranch;
        fJetBranch         = o.fJetBranch;
        fEntry             = o.fEntry;
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
        fHistList          = o.fHistList;
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

    // connect output aod
    // create a new branch for extra tracks
    fAODout = AODEvent();
    if(!fAODout){
        AliError("Output AOD not found.");
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

    OpenFile(1);
    if(!fHistList) fHistList = new TList();
    fHistList->SetOwner(kTRUE);

    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    fh1TrackPt      =  new TH1F("fh1TrackPt","pT of extra tracks;p_{T};entries", 250, 0., 250.);
    fh2TrackEtaPhi  =  new TH2F("fh2TrackEtaPhi","eta-phi distribution of extra tracks;#eta;#phi", 20, -1., 1., 60, 0., 2*TMath::Pi());
    fh1TrackN       =  new TH1F("fh1TrackN", "nb. of extra tracks per event;nb. of tracks;entries",300, 0., 300.);

    fHistList->Add(fh1TrackPt);
    fHistList->Add(fh2TrackEtaPhi);
    fHistList->Add(fh1TrackN);
	
	if(fEmbedMode==kAODFull || fEmbedMode==kAODJetTracks || fEmbedMode==kAODJet4Mom){
	
	    fh1JetPt        =  new TH1F("fh1JetPt", "pT of extra jets;p_{T};entries", 250, 0., 250.);
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
	
    fh1AODfile = new TH1I("fh1AODfile", "overview of opened AOD files from the array", 100, 0, 100);
    fHistList->Add(fh1AODfile);

    // =========== Switch on Sumw2 for all histos ===========
    for (Int_t i=0; i<fHistList->GetEntries(); ++i) {
        TH1 *h1 = dynamic_cast<TH1*>(fHistList->At(i));
        if (h1){
          h1->Sumw2();
          continue;
        }
    }

    TH1::AddDirectory(oldStatus);


    // set seed
    rndm = new TRandom3();
    Int_t id = GetJobID();
    if(id>-1) rndm->SetSeed(id);
    else      rndm->SetSeed();   // a TTUID is generated and used for seed
    AliInfo(Form("TRandom3 seed: %d", rndm->GetSeed()));

    // embed mode with AOD
    if(fEmbedMode==kAODFull || fEmbedMode==kAODJetTracks || fEmbedMode==kAODJet4Mom){

       // open input AOD
       Int_t rc = OpenAODfile();
       if(rc<0) return;
       fh1AODfile->Fill(rc);

    } //end: embed mode with AOD

   
    
   PostData(1, fHistList);
}


//__________________________________________________________________________
void AliAnalysisTaskFastEmbedding::Init()
{
    // Initialization
    if(fDebug > 1) Printf("AliAnalysisTaskFastEmbedding::Init()");

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

       // fetch jets
       TClonesArray *aodJets = 0;
       if(fJetBranch.Length()) aodJets = dynamic_cast<TClonesArray*>(fAODevent->FindListObject(fJetBranch.Data()));
       else                    aodJets = fAODevent->GetJets();
       if(!aodJets){
          AliError("Could not find jets in AOD. Check jet branch when indicated.");
          PostData(1, fHistList);
          return;
       }
        

       Int_t nEvents = fAODtree->GetEntries();

       Bool_t useEntry = kFALSE;
       while(!useEntry){  // protection need, if no event fulfills requierment
          if(fEntry>=nEvents){
              fEntry=0;
              if(!fAODPathArray){
                 AliDebug(AliLog::kDebug, "Last event in AOD reached, start from entry 0 again.");
              } 
              else {
                 AliDebug(AliLog::kDebug, "Last event in AOD reached, select new AOD file ...");

                 Int_t rc = OpenAODfile();
                 if(rc<0) {
                    PostData(1, fHistList);
                    return;
                 }
                 fh1AODfile->Fill(rc);

		 // new file => we must use the new jet array
		 if(fJetBranch.Length()) aodJets = dynamic_cast<TClonesArray*>(fAODevent->FindListObject(fJetBranch.Data()));
		 else                    aodJets = fAODevent->GetJets();
		 if(!aodJets){
		   AliError("Could not find jets in AOD. Check jet branch when indicated.");
                   PostData(1, fHistList);
		   return;
		 }
              }
          }
    
          fAODtree->GetEvent(fEntry);

	  // jet pt selection
	  if(fEvtSelecMode==kEventsJetPt){
             Int_t nJets = aodJets->GetEntries();
             for(Int_t ij=0; ij<nJets; ++ij){
                 AliAODJet *jet = dynamic_cast<AliAODJet*>(aodJets->At(ij));
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

	  fEntry++;
       }
       AliInfo(Form("Use entry %d from extra AOD.", fEntry-1));


       TClonesArray *mcpartIN  = (TClonesArray*)(fAODevent->FindListObject("mcparticles"));
       TClonesArray *mcpartOUT = 0x0;
       if(mcpartIN){
	   mcpartOUT = (TClonesArray*)(fAODout->FindListObject(fMCparticlesBranch.Data()));
	   mcpartOUT->Delete();
       } else {
	   AliInfo("No extra MC particles found.");
       }
    

       if(fEmbedMode==kAODFull || fEmbedMode==kAODJetTracks){ // take all tracks or jet tracks
	   // loop over input tracks
	   // add to output aod
	   Int_t nTracks = 0;
	   Int_t nJets = aodJets->GetEntries();
	   Int_t nSelectedJets = 0;
       AliAODJet *leadJet = 0x0; // used for jet tracks only
           
           if(fEmbedMode==kAODFull){
        		   nTracks = fAODevent->GetNumberOfTracks();
				   
			       for(Int_t ij=0; ij<nJets; ++ij){
                       AliAODJet *jet = dynamic_cast<AliAODJet*>(aodJets->At(ij));
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
                  AliAODJet *jet = dynamic_cast<AliAODJet*>(aodJets->At(ij));
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

               fh1TrackPt->Fill(tmpTr->Pt());
               fh2TrackEtaPhi->Fill(tmpTr->Eta(), tmpTr->Phi());
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
	   Int_t nJets = aodJets->GetEntries();
           fh1TrackN->Fill((Float_t)nJets);
	   for(Int_t ij=0; ij<nJets; ++ij){
               AliAODJet *jet = dynamic_cast<AliAODJet*>(aodJets->At(ij));
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
     Int_t n = rndm->Integer(nFiles);
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

    if(trial>3){
        AliError("Could not open AOD files, give up ...");
        return -1;
    }
	
	Int_t rc = 0;
	if(fAODPathArray) rc = SelectAODfile();
	if(rc<0) return -1;

    TDirectory *owd = gDirectory;
    if (fAODfile)
      fAODfile->Close();
    fAODfile = TFile::Open(fAODPath.Data());
    owd->cd();
    if(!fAODfile){
	
       rc = -1;
       if(fAODPathArray){
           AliError(Form("Could not open AOD file %s, try another from the list ...", fAODPath.Data()));
           rc = OpenAODfile(trial+1);
       } else {
	       AliError(Form("Could not open AOD file %s.", fAODPath.Data()));
	   }
        
       return rc;
    }

    fAODtree = (TTree*)fAODfile->Get("aodTree");

    if(!fAODtree){
       AliError("AOD tree not found.");
       return -1;
    }

    delete fAODevent;
    fAODevent = new AliAODEvent();
    fAODevent->ReadFromTree(fAODtree);
    return rc;  // file position in AOD path array, if array available
}

