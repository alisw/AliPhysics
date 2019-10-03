// **************************************
// Task used for jet finding in the Kine Train (generation and analysis on the fly, no detector effects)
// Output is stored in an exchance container
// *******************************************


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

 
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TRefArray.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TF1.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include  "TDatabasePDG.h"
#include <TGrid.h>

#include "AliAnalysisTaskJetClusterKine.h"
#include "AliAnalysisManager.h"
#include "AliJetFinder.h"
#include "AliJetHeader.h"
#include "AliJetReader.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODExtension.h"
#include "AliAODTrack.h"
#include "AliAODJet.h"
#include "AliVParticle.h" //FK//
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h" //FK//
#include "AliGenPythiaEventHeader.h"
#include "AliJetKineReaderHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliAODJetEventBackground.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"
// get info on how fastjet was configured
#include "fastjet/config.h"

using std::vector;

ClassImp(AliAnalysisTaskJetClusterKine)

AliAnalysisTaskJetClusterKine::~AliAnalysisTaskJetClusterKine(){
  //
  // Destructor
  //

  delete fRef;

  if(fTCAJetsOut)fTCAJetsOut->Delete();
  delete fTCAJetsOut;
  
}

//_____________________________________________________________________

AliAnalysisTaskJetClusterKine::AliAnalysisTaskJetClusterKine(): 
  AliAnalysisTaskSE(),
  fMcEvent(0x0),  
  fMcHandler(0x0),
  fRef(new TRefArray),
  fTrackTypeGen(kTrackKineCharged), //Kine charged? 
  fAvgTrials(1),
  fTrackEtaWindow(0.9),    
  fTrackPtCut(0.),							
  fJetOutputMinPt(0.150),
  fMaxTrackPtInJet(100.),
  fVtxZCut(10.0),
  fNonStdBranch(""),
  fOutContainer(kNoOutput), //FF//
  fNonStdFile(""),
  fRparam(1.0), 
  fAlgorithm(fastjet::kt_algorithm),
  fStrategy(fastjet::Best),
  fRecombScheme(fastjet::BIpt_scheme),
  fAreaType(fastjet::active_area), 
  fGhostArea(0.01),
  fActiveAreaRepeats(1),
  fGhostEtamax(1.5),
  fTCAJetsOut(0x0),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1NJetsGen(0x0),
  fh1NConstGen(0x0),
  fh1NConstLeadingGen(0x0),
  fh1PtJetsGenIn(0x0),
  fh1PtJetsLeadingGenIn(0x0),
  fh1PtJetConstGen(0x0),
  fh1PtJetConstLeadingGen(0x0),
  fh1PtTracksGenIn(0x0),
  fh1Nch(0x0),
  fh1Z(0x0), 
  fh2NConstPt(0x0),
  fh2NConstLeadingPt(0x0),
  fh2JetPhiEta(0x0),
  fh2LeadingJetPhiEta(0x0),
  fh2JetEtaPt(0x0),
  fh2LeadingJetEtaPt(0x0),
  fh2TrackEtaPt(0x0),
  fh2JetsLeadingPhiEta(0x0),
  fh2JetsLeadingPhiPt(0x0),
  fh2JetsLeadingPhiPtW(0x0),
  fHistList(0x0)  
{
  //
  // Constructor
  //
}

//_____________________________________________________________________

AliAnalysisTaskJetClusterKine::AliAnalysisTaskJetClusterKine(const char* name):
  AliAnalysisTaskSE(name),
  fMcEvent(0x0),  
  fMcHandler(0x0),
  fRef(new TRefArray),
  fTrackTypeGen(kTrackKineCharged),  //kine charged?
  fAvgTrials(1),
  fTrackEtaWindow(0.9),    
  fTrackPtCut(0.),							
  fJetOutputMinPt(0.150),
  fMaxTrackPtInJet(100.),
  fVtxZCut(10.0),
  fNonStdBranch(""),
  fOutContainer(kNoOutput),//FF//
  fNonStdFile(""),
  fRparam(1.0), 
  fAlgorithm(fastjet::kt_algorithm),
  fStrategy(fastjet::Best),
  fRecombScheme(fastjet::BIpt_scheme),
  fAreaType(fastjet::active_area), 
  fGhostArea(0.01),
  fActiveAreaRepeats(1),
  fGhostEtamax(1.5),
  fTCAJetsOut(0x0),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1NJetsGen(0x0),
  fh1NConstGen(0x0),
  fh1NConstLeadingGen(0x0),
  fh1PtJetsGenIn(0x0),
  fh1PtJetsLeadingGenIn(0x0),
  fh1PtJetConstGen(0x0),
  fh1PtJetConstLeadingGen(0x0),
  fh1PtTracksGenIn(0x0),
  fh1Nch(0x0),
  fh1Z(0x0), 
  fh2NConstPt(0x0),
  fh2NConstLeadingPt(0x0),
  fh2JetPhiEta(0x0),
  fh2LeadingJetPhiEta(0x0),
  fh2JetEtaPt(0x0),
  fh2LeadingJetEtaPt(0x0),
  fh2TrackEtaPt(0x0),
  fh2JetsLeadingPhiEta(0x0),
  fh2JetsLeadingPhiPt(0x0),
  fh2JetsLeadingPhiPtW(0x0),
  fHistList(0x0)
{
  //
  // named ctor
  //
  DefineOutput(1, TList::Class());  
  DefineOutput(2, TClonesArray::Class());  
}


//_____________________________________________________________________

Bool_t AliAnalysisTaskJetClusterKine::Notify(){
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 
  return kTRUE;
}

//_____________________________________________________________________

void AliAnalysisTaskJetClusterKine::UserCreateOutputObjects(){

   //
   // Create the output container
   //

   // Connect the AOD


   if(fDebug > 1) printf("AnalysisTaskJetCluster::UserCreateOutputObjects() \n");


   if(fNonStdBranch.Length()!=0){
      // only create the output branch if we have a name
      // Create a new branch for jets...
      //  -> cleared in the UserExec....
      // here we can also have the case that the brnaches are written to a separate file
      fTCAJetsOut = new TClonesArray("AliAODJet", 0);
      fTCAJetsOut->SetName(fNonStdBranch.Data());
      if(fOutContainer==kAODBranch){   //FF//
         AddAODBranch("TClonesArray",&fTCAJetsOut,fNonStdFile.Data());
      }

    
      //if(fNonStdFile.Length()!=0){
	// 
	// case that we have an AOD extension we need to fetch the jets from the extended output
	// we identify the extension aod event by looking for the branchname
	//AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
	// case that we have an AOD extension we need can fetch the background maybe from the extended output                                                                  
	//fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);
      //}
   }


   if(!fHistList) fHistList = new TList();
   fHistList->SetOwner(kTRUE);
   PostData(1, fHistList); // post data in any case once


   if(fOutContainer==kExchCont){
      fTCAJetsOut->SetOwner();
      PostData(2, fTCAJetsOut); //FF// post data in any case once
   }

   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);


   //
   //  Histogram
    
   const Int_t nBinPt = 100;
   Double_t binLimitsPt[nBinPt+1];
   for(Int_t iPt = 0;iPt <= nBinPt;iPt++){
      if(iPt == 0){
         binLimitsPt[iPt] = 0.0;
      }else {// 1.0
         binLimitsPt[iPt] =  binLimitsPt[iPt-1] + 2.0;
      }
   }
  
   const Int_t nBinPhi = 90;
   Double_t binLimitsPhi[nBinPhi+1];
   for(Int_t iPhi = 0; iPhi<=nBinPhi; iPhi++){
       if(iPhi==0){
          binLimitsPhi[iPhi] = -1.*TMath::Pi();
       }else{
          binLimitsPhi[iPhi] = binLimitsPhi[iPhi-1] + 1/(Float_t)nBinPhi * TMath::Pi()*2;
       }
    }
   
    const Int_t nBinEta = 40;
    Double_t binLimitsEta[nBinEta+1];
    for(Int_t iEta = 0;iEta<=nBinEta;iEta++){
       if(iEta==0){
          binLimitsEta[iEta] = -2.0;
       }else{
          binLimitsEta[iEta] = binLimitsEta[iEta-1] + 0.1;
     }
   }
   
   const int nChMax = 5000;
   
   fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
   fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  
   fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
   fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
   
   
   fh1NJetsGen = new TH1F("fh1NJetsGen","N reconstructed jets",120,-0.5,119.5);
   
   fh1NConstGen = new TH1F("fh1NConstGen","# jet constituents",120,-0.5,119.5);
   fh1NConstLeadingGen = new TH1F("fh1NConstLeadingGen","jet constituents",120,-0.5,119.5);
   
   
   fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",nBinPt,binLimitsPt);
   fh1PtHardNoW = new TH1F("fh1PtHardNoW","PYTHIA Pt hard no weight;p_{T,hard}",nBinPt,binLimitsPt);
   fh1PtHardTrials = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",nBinPt,binLimitsPt);
   
   fh1PtJetsGenIn  = new TH1F("fh1PtJetsGenIn","Gen jets P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
   fh1PtJetsLeadingGenIn = new TH1F("fh1PtJetsLeadingGenIn","Gen jets P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
   fh1PtJetConstGen = new TH1F("fh1PtJetsConstGen","Rec jets constituents P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
   fh1PtJetConstLeadingGen = new TH1F("fh1PtJetsConstLeadingGen","Gen jets constituents P_T;p_{T} (GeV/c)",nBinPt,binLimitsPt);
   fh1PtTracksGenIn  = new TH1F("fh1PtTracksGenIn",Form("Gen tracks P_T #eta < %1.2f ;p_{T} (GeV/c)",fTrackEtaWindow),nBinPt,binLimitsPt);
   fh1Nch = new TH1F("fh1Nch","charged multiplicity; N_{ch}",nChMax,-0.5,nChMax-0.5);
   
   
   fh1Z = new TH1F("fh1Z",";zvtx",100,-25,25);
   
   
 
   fh2NConstPt = new TH2F("fh2NConstPt","Number of constituents ;p_{T} (GeV/c);N",nBinPt,binLimitsPt,50,-0.5,49.5);
   fh2NConstLeadingPt = new TH2F("fh2NConstLeadingPt","Number of constituents ;p_{T} (GeV/c);N",nBinPt,binLimitsPt,50,-0.5,49.5);
 
 
   fh2JetPhiEta  = new TH2F("fh2JetPhiEta","eta vs phi all jets;#phi;#eta",
 			   nBinPhi,0.,2.*TMath::Pi(),nBinEta,binLimitsEta);
   fh2LeadingJetPhiEta  = new TH2F("fh2LeadingJetPhiEta","eta vs phi leading jets;#phi;#eta",
 				  nBinPhi,0.,2.*TMath::Pi(),nBinEta,binLimitsEta);
 
   fh2JetEtaPt  = new TH2F("fh2JetEtaPt","pt vs eta all jets;#eta;p_{T}",
 			  nBinEta,binLimitsEta,nBinPt,binLimitsPt);
   fh2LeadingJetEtaPt  = new TH2F("fh2LeadingJetEtaPt","pT vs eta leading jets;#eta;p_{T}",
 				 nBinEta,binLimitsEta,nBinPt,binLimitsPt);
 
   fh2TrackEtaPt  = new TH2F("fh2TrackEtaPt","pt vs eta all jets;#eta;p_{T}",
 			    nBinEta,binLimitsEta,nBinPt,binLimitsPt);
 
 
 
   fh2JetsLeadingPhiEta = new TH2F("fh2JetsLeadingPhiEta","delta eta vs delta phi to leading jet;#Delta#phi;#Delta#eta",
 				  nBinPhi,binLimitsPhi,nBinEta,binLimitsEta);
   fh2JetsLeadingPhiPt = new TH2F("fh2JetsLeadingPhiPt","leading p_T vs delta phi to leading jet;#Delta#phi;p_{T} (GeV/c)",
 				 nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);

   fh2JetsLeadingPhiPtW = new TH2F("fh2JetsLeadingPhiPtW","leading p_T vs delta phi p_T weigted to leading jet;#Delta#phi;p_{T} (GeV/c)",
		                 nBinPhi,binLimitsPhi,nBinPt,binLimitsPt);


  

   const Int_t saveLevel = 3; // large save level more histos
   if(saveLevel>0){
      fHistList->Add(fh1Xsec);
      fHistList->Add(fh1Trials);
      
      fHistList->Add(fh1NJetsGen);
      fHistList->Add(fh1NConstGen);
      fHistList->Add(fh1NConstLeadingGen);
      fHistList->Add(fh1PtJetsGenIn);
      fHistList->Add(fh1PtJetsLeadingGenIn);
      fHistList->Add(fh1PtTracksGenIn);
      fHistList->Add(fh1PtJetConstGen);
      fHistList->Add(fh1PtJetConstLeadingGen);
      fHistList->Add(fh1Nch);
      fHistList->Add(fh1Z);
      
      fHistList->Add(fh2NConstPt);
      fHistList->Add(fh2NConstLeadingPt);
      fHistList->Add(fh2JetPhiEta);
      fHistList->Add(fh2LeadingJetPhiEta);
      fHistList->Add(fh2JetEtaPt);
      fHistList->Add(fh2LeadingJetEtaPt);
      fHistList->Add(fh2TrackEtaPt);
      fHistList->Add(fh2JetsLeadingPhiEta);
      fHistList->Add(fh2JetsLeadingPhiPt);
      fHistList->Add(fh2JetsLeadingPhiPtW);
   }

   // =========== Switch on Sumw2 for all histos ===========
   for(Int_t i=0; i<fHistList->GetEntries(); ++i){
      TH1 *h1 = dynamic_cast<TH1*>(fHistList->At(i));
      if(h1){
         h1->Sumw2();
         continue;
      }
      THnSparse *hn = dynamic_cast<THnSparse*>(fHistList->At(i));
      if(hn) hn->Sumw2();
   }
   TH1::AddDirectory(oldStatus);
}
//_________________________________________________________________________

void AliAnalysisTaskJetClusterKine::Init(){
   // MC handler
   if(fDebug > 1) printf("AnalysisTaskJetClusterKine::Init() \n");
   fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()); //FK//
}

//_________________________________________________________________________
void AliAnalysisTaskJetClusterKine::LocalInit(){
   // MC handler
   //
   // Initialization
   //

   if(fDebug > 1) printf("AnalysisTaskJetClusterKine::LocalInit() \n");
   Init();
}

//_________________________________________________________________________
void AliAnalysisTaskJetClusterKine::UserExec(Option_t* /*option*/){
   // handle and reset the output jet branch 
  
   if(fDebug > 1) printf("AliAnalysisTaskJetClusterKine::UserExec \n");
   if(fTCAJetsOut) fTCAJetsOut->Delete(); //clean TClonesArray

   //
   // Execute analysis for current event
   //
   Init();
   if(fMcHandler){
      fMcEvent = fMcHandler->MCEvent(); 
   }else{
      if(fDebug > 1) printf("AnalysisTaskJetClusterKine::Handler() fMcHandler=NULL\n");
      PostData(1, fHistList);

      if(fOutContainer==kExchCont){
         PostData(2, fTCAJetsOut); //FF//
      }

      return;
   }
   if(!fMcEvent){
      if(fDebug > 1) printf("AnalysisTaskJetClusterKine::Exec()   fMcEvent=NULL \n");
      PostData(1, fHistList);

      if(fOutContainer==kExchCont){
         PostData(2, fTCAJetsOut); //FF// 
      }

      return;
   }

   const AliVVertex *vtxMC = fMcEvent->GetPrimaryVertex();
   Float_t zVtx = vtxMC->GetZ();
   if(TMath::Abs(zVtx) > fVtxZCut){ //vertex cut
      PostData(1, fHistList);

      if(fOutContainer==kExchCont){
         PostData(2, fTCAJetsOut); //FF// 
      }

      return;
   }

   fh1Z->Fill(zVtx);
 
   fh1Trials->Fill("#sum{ntrials}",1);
   if(fDebug > 10) Printf("%s:%d",(char*)__FILE__,__LINE__);

   // we simply fetch the mc particles as a list of AliVParticles

   TList genParticles; //list of particles  with above min pT and  good eta range

   Int_t nParticles = GetListOfTracks(&genParticles, fTrackTypeGen); //fill the list
   fh1Nch->Fill(nParticles);

   if(fDebug>2) Printf("%s:%d Selected Gen tracks: %d %d",(char*)__FILE__,__LINE__, nParticles, genParticles.GetEntries());

   // find the jets....

   vector<fastjet::PseudoJet> inputParticlesKine; 

   for(int i = 0; i < genParticles.GetEntries(); i++){  //loop over generated particles
      AliVParticle *vp = (AliVParticle*) genParticles.At(i);
      fastjet::PseudoJet jInp(vp->Px(),vp->Py(),vp->Pz(),vp->P());
      jInp.set_user_index(i);
      inputParticlesKine.push_back(jInp);

      if(fTCAJetsOut){
         if(i == 0){
            fRef->Delete(); // make sure to delete before placement new...
            new(fRef) TRefArray(TProcessID::GetProcessWithUID(vp)); //TRefArray does not work with toy model ...
         }
         fRef->Add(vp);  //TRefArray does not work with toy model ...
      }
   }//generated particles

   if(inputParticlesKine.size()==0){ //FK//
      if(fDebug)Printf("%s:%d No input particles found, skipping event",(char*)__FILE__,__LINE__);
      PostData(1, fHistList);

      if(fOutContainer==kExchCont){
         PostData(2, fTCAJetsOut); //FF// 
      }

      return;
   } //FK//
 
   // run fast jet
   // employ setters for these...
   // now create the object that holds info about ghosts                        

   fastjet::GhostedAreaSpec ghostSpec(fGhostEtamax, fActiveAreaRepeats, fGhostArea);
   fastjet::AreaType areaType =   fastjet::active_area;
   fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType,ghostSpec);
   fastjet::JetDefinition jetDef(fAlgorithm, fRparam, fRecombScheme, fStrategy);
   fastjet::ClusterSequenceArea clustSeq(inputParticlesKine, jetDef, areaDef); //FK//
  
   //range where to compute background
   Double_t phiMin = 0, phiMax = 0, rapMin = 0, rapMax = 0;
   phiMin = 0;
   phiMax = 2*TMath::Pi();
   rapMax = fGhostEtamax - fRparam;
   rapMin = - fGhostEtamax + fRparam;
   fastjet::RangeDefinition range(rapMin,rapMax, phiMin, phiMax);
 

   const vector <fastjet::PseudoJet> &inclusiveJets = clustSeq.inclusive_jets();
   const vector <fastjet::PseudoJet> &sortedJets = sorted_by_pt(inclusiveJets);

 
   fh1NJetsGen->Fill(sortedJets.size()); //the number of found jets

   // loop over all jets an fill information, first on is the leading jet

   Int_t nGen     = inclusiveJets.size();
   if(inclusiveJets.size()>0){

      //leading Jet
      AliAODJet leadingJet (sortedJets[0].px(), sortedJets[0].py(), sortedJets[0].pz(), sortedJets[0].E());
      Double_t area = clustSeq.area(sortedJets[0]);
      leadingJet.SetEffArea(area,0);
      Float_t pt = leadingJet.Pt();
      Float_t phi = leadingJet.Phi();
      if(phi<0) phi += TMath::Pi()*2.;    
      Float_t eta = leadingJet.Eta();

      //inclusive jet
      Int_t nAodOutJets = 0;
      AliAODJet *aodOutJet = NULL;


      // correlation of leading jet with tracks
  //FK//    TIterator *particleIter = genParticles.MakeIterator();//FK//
  //FK//    particleIter->Reset();
  //FK//    AliVParticle *tmpParticle = NULL;
   //FK//   while((tmpParticle = (AliVParticle*)(particleIter->Next()))){
   //FK//      Float_t tmpPt = tmpParticle->Pt();
    //FK//     // correlation
    //FK//     Float_t tmpPhi =  tmpParticle->Phi();     
   //FK//      if(tmpPhi<0) tmpPhi+= TMath::Pi()*2.;    
   //FK//      Float_t dPhi = phi - tmpPhi;
     //FK//    if(dPhi > TMath::Pi())    dPhi = dPhi - 2.*TMath::Pi(); //-pi,pi
    //FK//     if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();  //-pi,pi     
     //FK//    fh2TracksLeadingJetPhiPt->Fill(dPhi,pt);
     //FK//    fh2TracksLeadingJetPhiPtW->Fill(dPhi,pt,tmpPt);
    //FK//  }  
   
      TLorentzVector vecareab;
      for(int j = 0; j < nGen;j++){ //loop over inclusive jets
         AliAODJet tmpGenJet (sortedJets[j].px(), sortedJets[j].py(), sortedJets[j].pz(), sortedJets[j].E());
         aodOutJet = NULL;
         Float_t tmpPt = tmpGenJet.Pt(); //incl jet pt 
      
         if((tmpPt > fJetOutputMinPt) && fTCAJetsOut){// cut on the non-background subtracted...
	    aodOutJet =  new ((*fTCAJetsOut)[nAodOutJets++]) AliAODJet(tmpGenJet);
	    aodOutJet->GetRefTracks()->Clear();
	    Double_t area1 = clustSeq.area(sortedJets[j]);
	    aodOutJet->SetEffArea(area1,0);
            fastjet::PseudoJet vecarea=clustSeq.area_4vector(sortedJets[j]);  
            vecareab.SetPxPyPzE(vecarea.px(),vecarea.py(),vecarea.pz(),vecarea.e());     
	    aodOutJet->SetVectorAreaCharged(&vecareab);
         }

         fh1PtJetsGenIn->Fill(tmpPt); //incl jet Pt
         // Fill Spectra with constituentsemacs
         const vector<fastjet::PseudoJet> &constituents = clustSeq.constituents(sortedJets[j]);

         fh1NConstGen->Fill(constituents.size());  //number of constituents
         fh2NConstPt->Fill(tmpPt,constituents.size()); //number of constituents vs jet pt 

         // loop over constiutents and fill spectrum
         AliVParticle *partLead = 0;
         Float_t ptLead = -1;

         for(unsigned int ic = 0; ic < constituents.size(); ic++){
  	    AliVParticle *part = (AliVParticle*) genParticles.At(constituents[ic].user_index());
	    if(!part) continue;
	    fh1PtJetConstGen->Fill(part->Pt()); //pt of constituents

	    if(aodOutJet){
	       aodOutJet->AddTrack(fRef->At(constituents[ic].user_index()));//FK//

	       if(part->Pt()>fMaxTrackPtInJet){
	          aodOutJet->SetTrigger(AliAODJet::kHighTrackPtTriggered);
	       }
	    }

  	    if(part->Pt()>ptLead){
	       ptLead = part->Pt();
	       partLead = part;
	    }

	    if(j==0) fh1PtJetConstLeadingGen->Fill(part->Pt()); //pt of leading jet constituents
         }

         if(partLead){
   	    if(aodOutJet){
	       //set pT of leading constituent of jet
	       aodOutJet->SetPtLeading(partLead->Pt());
  	    }
         }
    
         // correlation
         Float_t tmpPhi =  tmpGenJet.Phi(); //incl jet phi
         Float_t tmpEta =  tmpGenJet.Eta(); //incl jet eta

         if(tmpPhi<0) tmpPhi += TMath::Pi()*2.;        
         fh2JetPhiEta->Fill(tmpGenJet.Phi(),tmpEta);
         fh2JetEtaPt->Fill(tmpEta,tmpPt);

         if(j==0){   //leading jet
	    fh1PtJetsLeadingGenIn->Fill(tmpPt); //leading jet pt
	    fh2LeadingJetPhiEta->Fill(tmpPhi,tmpEta);
	    fh2LeadingJetEtaPt->Fill(tmpEta,tmpPt);
	    fh1NConstLeadingGen->Fill(constituents.size());  //number of constituents in leading jet
	    fh2NConstLeadingPt->Fill(tmpPt,constituents.size());
	    continue;
         }

         Float_t dPhi = phi - tmpPhi;
         if(dPhi > TMath::Pi())     dPhi = dPhi - 2.*TMath::Pi(); //-pi,pi
         if(dPhi<(-1.*TMath::Pi())) dPhi = dPhi + 2.*TMath::Pi(); //-pi.pi     
         Float_t dEta = eta - tmpGenJet.Eta();
         fh2JetsLeadingPhiEta->Fill(dPhi,dEta);
         fh2JetsLeadingPhiPt->Fill(dPhi,pt);
         fh2JetsLeadingPhiPtW->Fill(dPhi,pt,tmpPt);
      }// loop over reconstructed jets
    //delete particleIter;
   }//number of jets>0

 
   // fill track information
   Int_t nTrackOver = genParticles.GetSize();
   // do the same for tracks and jets

   if(nTrackOver>0){
      TIterator *particleIter = genParticles.MakeIterator();
      AliVParticle *tmpGen = 0;
      particleIter->Reset();

      while((tmpGen = (AliVParticle*)(particleIter->Next()))){
         Float_t tmpPt = tmpGen->Pt();
         Float_t tmpEta = tmpGen->Eta();
         fh1PtTracksGenIn->Fill(tmpPt);
         fh2TrackEtaPt->Fill(tmpEta,tmpPt);
      }  
      delete particleIter;
   }





   if(fDebug > 2){
      if(fTCAJetsOut)Printf("%s:%d Rec Jets %d",(char*)__FILE__,__LINE__,fTCAJetsOut->GetEntriesFast());
   }
   PostData(1, fHistList);
 
  if(fOutContainer==kExchCont){
      PostData(2, fTCAJetsOut); //FF// 
   }


}
//____________________________________________________________________________

void AliAnalysisTaskJetClusterKine::Terminate(Option_t* /*option*/){
  //
  // Terminate analysis
  //
  if (fDebug > 1) printf("AnalysisJetCluster: Terminate() \n");

    
}

//_____________________________________________________________________________________

Int_t  AliAnalysisTaskJetClusterKine::GetListOfTracks(TList *list, Int_t type){

  //
  // get list of tracks/particles for different types
  //

   if(fDebug>2) Printf("%s:%d Selecting tracks with %d",(char*)__FILE__,__LINE__,type);

   Int_t iCount = 0; //number of particles 

   if(type ==  kTrackKineAll || type == kTrackKineCharged){
      if(! fMcEvent) return iCount;

      //we want to have alivpartilces so use get track
      for(int it = 0; it < fMcEvent->GetNumberOfTracks(); ++it){
         if(!fMcEvent->IsPhysicalPrimary(it)) continue;

         AliMCParticle* part = (AliMCParticle*)fMcEvent->GetTrack(it);
         if(TMath::Abs(part->Eta()) > fTrackEtaWindow) continue;   //apply eta cut
         if(part->Pt() < fTrackPtCut) continue;  //apply pT cut

         if(type == kTrackKineAll){  // full jets
	    list->Add(part);
	    iCount++;

         }else if(type == kTrackKineCharged){  //charged jets
            if(part->Particle()->GetPDG()->Charge()==0) continue;
	    list->Add(part);
  	    iCount++;
         }
      }
   }
  
   list->Sort();  //?//

   return iCount; //number of particles in the list
}


