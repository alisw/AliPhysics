/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
* Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes hereby granted      *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id: $ */

//_________________________________________________________________________
//
// Class for the electron identification and B-tagging.
// Clusters from EMCAL matched to tracks to id electrons.
// Btagger is run on all electrons, then jets are tagged as well.
// 
//
// -- Author: T.R.P.Aronsson (Yale), M. Heinz (Yale)
//////////////////////////////////////////////////////////////////////////////
  
#include <TH2F.h>
#include <TH3F.h>
#include <TParticle.h>
#include <TNtuple.h>
#include <TClonesArray.h>

#include "AliAnaBtag.h" 
#include "AliCaloTrackReader.h"
#include "AliMCAnalysisUtils.h"
#include "AliVCluster.h"
#include "AliFiducialCut.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliCaloPID.h"
#include "AliAODMCParticle.h"
#include "AliStack.h"
#include "AliExternalTrackParam.h"
#include "AliESDv0.h"
#include "AliESDtrack.h"
#include "AliAODJet.h"
#include "AliAODEvent.h"
#include "AliGenPythiaEventHeader.h"
ClassImp(AliAnaBtag)
  
//____________________________________________________________________________
AliAnaBtag::AliAnaBtag() 
: AliAnaPartCorrBaseClass(),
  fWriteNtuple(0),electrons(0),pairs(0),events(0),fEventNumber(0),fNElec(0),fNElecEv(0),fNPair(0),fDrCut(0.),fPairDcaCut(0.),fDecayLenCut(0.),fImpactCut(0.),
  fAssocPtCut(0.),fMassCut(0.),fSdcaCut(0.),fITSCut(0),
  fNTagTrkCut(0),fIPSigCut(0.),fJetEtaCut(0.3),fJetPhiMin(1.8),fJetPhiMax(2.9),
  fhEmcalElectrons(0),fhTRDElectrons(0),fhTPCElectrons(0),fhEmcalMCE(0),fhTRDMCE(0),fhTPCMCE(0),fhEmcalMCEFake(0),fhTRDMCEFake(0),fhTPCMCEFake(0),fhEmcalMCP(0),fhTRDMCP(0),fhTPCMCP(0),fhEmcalMCK(0),fhTRDMCK(0),fhTPCMCK(0),fhSpecies(0),fhDVM1(0),fhDVM2(0),fhNVTX(0),fhNVTXMC(0),fhJets(0),fhJetsAllEtaPhi(0),fhJetsLeadingBElectronEtaPhi(0),fhDVM1EtaPhi(0),fhBJetElectronDetaDphi(0),fhClusterMap(0),fhClusterEnergy(0),fhTestalle(0),fhResidual(0),fhPairPt(0),fhElectrons(0),fhTracks(0)
{
  //default ctor
  //Initialize parameters
  InitParameters();

}

//____________________________________________________________________________
AliAnaBtag::~AliAnaBtag() 
{
  //dtor

}


//________________________________________________________________________
TList *  AliAnaBtag::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ElectronHistos") ; 

  if(fWriteNtuple){
    electrons = new TNtuple("electrons","Electron Ntuple","electronnumber:eventnumber:electronineventnumber:pairs:start:stop:tpc:trd:emc:mcpdg:parentbit:btag:pt");
    outputContainer->Add(electrons);

    pairs = new TNtuple("pairs","Pair Ntuple","pairnumber:electronnumber:eventnumber:pdca:sdca:minv:pth:massphoton:decaylength");
    outputContainer->Add(pairs);

    events = new TNtuple("events","Event NTuple","event");
    outputContainer->Add(events);
  }
  
  fhEmcalElectrons = new TH1F("fhEmcalElectrons","",400,0,400);
  outputContainer->Add(fhEmcalElectrons);
  
  fhTRDElectrons = new TH1F("fhTRDElectrons","",400,0,400);
  outputContainer->Add(fhTRDElectrons);
  
  fhTPCElectrons = new TH1F("fhTPCElectrons","",400,0,400);
  outputContainer->Add(fhTPCElectrons);
  
  fhNVTX = new TH1F("fhNVTX","",20,0,20);
  outputContainer->Add(fhNVTX);

  fhDVM1 = new TH1F("fhDVM1","",400,0,400);
  outputContainer->Add(fhDVM1);
  
  fhDVM2 = new TH1F("fhDVM2","",400,0,400);
  outputContainer->Add(fhDVM2);
  
  if(IsDataMC()){
    fhEmcalMCE = new TH1F("fhEmcalMCE","",400,0,400);
    outputContainer->Add(fhEmcalMCE);
    
    fhTRDMCE = new TH1F("fhTRDMCE","",400,0,400);
    outputContainer->Add(fhTRDMCE);
    
    fhTPCMCE = new TH1F("fhTPCMCE","",400,0,400);
    outputContainer->Add(fhTPCMCE);
    
    fhEmcalMCEFake = new TH1F("fhEmcalMCEFake","",400,0,400);
    outputContainer->Add(fhEmcalMCEFake);
    
    fhTRDMCEFake = new TH1F("fhTRDMCEFake","",400,0,400);
    outputContainer->Add(fhTRDMCEFake);
    
    fhTPCMCEFake = new TH1F("fhTPCMCEFake","",400,0,400);
    outputContainer->Add(fhTPCMCEFake);
    
    fhEmcalMCP = new TH1F("fhEmcalMCP","",400,0,400);
    outputContainer->Add(fhEmcalMCP);
    
    fhTRDMCP = new TH1F("fhTRDMCP","",400,0,400);
    outputContainer->Add(fhTRDMCP);
    
    fhTPCMCP = new TH1F("fhTPCMCP","",400,0,400);
    outputContainer->Add(fhTPCMCP);
    
    fhEmcalMCK = new TH1F("fhEmcalMCK","",400,0,400);
    outputContainer->Add(fhEmcalMCK);
    
    fhTRDMCK = new TH1F("fhTRDMCK","",400,0,400);
    outputContainer->Add(fhTRDMCK);
    
    fhTPCMCK = new TH1F("fhTPCMCK","",400,0,400);
    outputContainer->Add(fhTPCMCK);
    
    fhSpecies  = new TH1F("fhSpecies","",1000,0,1000);
    outputContainer->Add(fhSpecies);
    
    fhNVTXMC = new TH1F("fhNVTXMC","",20,0,20);
    outputContainer->Add(fhNVTXMC);
  }
  

  fhJets = new TH2F("fhJets","",400,0,400,20,0,20);
  outputContainer->Add(fhJets);
  
  fhJetsAllEtaPhi = new TH2F("fhJetsAllEtaPhi","",100,-2,2,100,-2,8);
  outputContainer->Add(fhJetsAllEtaPhi);
  
  fhJetsLeadingBElectronEtaPhi = new TH2F("fhJetsLeadingBElectronEtaPhi","",100,-5,5,200,-10,10);
  outputContainer->Add(fhJetsLeadingBElectronEtaPhi);
  
  fhDVM1EtaPhi = new TH2F("fhDVM1EtaPhi","",100,-2,2,100,-2,8);
  outputContainer->Add(fhDVM1EtaPhi);
  
  fhBJetElectronDetaDphi = new TH2F("fhBJetElectronDetaDphi","",100,-5,5,200,-10,10);
  outputContainer->Add(fhBJetElectronDetaDphi);
  
  fhClusterMap = new TH2F("fhClusterMap","",100,-2,2,100,-2,8);
  outputContainer->Add(fhClusterMap);
  
  fhClusterEnergy = new TH1F("fhClusterEnergy","",100,0,10);
  outputContainer->Add(fhClusterEnergy);
  
  fhTestalle = new TH1F("fhTestalle","",400,0,400);
  outputContainer->Add(fhTestalle);
  
  fhResidual = new TH1F("fhResidual","",500,0,5);
  outputContainer->Add(fhResidual);

  fhPairPt = new TH1F("fhPairPt","",400,0,400);
  outputContainer->Add(fhPairPt);
  
  fhElectrons = new TH2F("fhElectrons","",200,0,100,20,0,20);
  outputContainer->Add(fhElectrons);
  
  fhTracks = new TH2F("fhTracks","",200,0,100,20,0,20);
  outputContainer->Add(fhTracks);
    
  return outputContainer ; 
}

//____________________________________________________________________________
void AliAnaBtag::Init()
{
  //do some initialization
    printf("Rubbish init step AliAnaBtag::Init()");
}


//____________________________________________________________________________
void AliAnaBtag::InitParameters()
{ 
  //Initialize the parameters of the analysis.
  SetOutputAODClassName("AliAODPWG4Particle");
  SetOutputAODName("PWG4Particle");
  //AddToHistogramsName("");
  //DVM B-tagging
  fDrCut       = 1.0; 
  fPairDcaCut  = 0.02;
  fDecayLenCut = 1.0;
  fImpactCut   = 0.5;
  fAssocPtCut  = 1.0;
  fMassCut     = 1.5;
  fSdcaCut     = 0.1;
  fITSCut      = 4;
  //IPSig B-tagging
  fNTagTrkCut  = 2;
  fIPSigCut    = 3.0;
  //Jet fiducial cuts
  fJetEtaCut = 0.3;
  fJetPhiMin = 1.8;
  fJetPhiMax = 2.9;
}

//__________________________________________________________________
void  AliAnaBtag::MakeAnalysisFillAOD() 
{
  fEventNumber++;
  fNElecEv=0;
  if(fWriteNtuple)
    events->Fill(fEventNumber);
  
  //This reads in tracks, extrapolates to EMCAL, does p/E selectrons, identifies electron candidates
  //After candidates are obtained, btagging and saving into AOD.
  AliStack *stack =0x0;
  //Double_t bfield = 0.;
  if(GetDebug()>0) printf("AliAnaBtag::MakeAnalysisFillAOD() - Write ntuple flag is %d \n",fWriteNtuple);
  //if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) bfield = GetReader()->GetBField();
  TObjArray *cl = GetEMCALClusters();
  
  if(!GetCTSTracks() || GetCTSTracks()->GetEntriesFast() == 0) return ;
  Int_t ntracks = GetCTSTracks()->GetEntriesFast();
  if(GetDebug() > 0)
    printf("AliAnaBtag::MakeAnalysisFillAOD() - In CTS aod entries %d\n", ntracks);
  
  Int_t iCluster = -999;
  Int_t ntot = cl->GetEntriesFast();
  
  //CLUSTER STUFF 
  for(Int_t iclus = 0; iclus < ntot; iclus++) {
    AliVCluster * clus = (AliVCluster*) (cl->At(iclus));
    if(!clus) continue;
    fhClusterEnergy->Fill(clus->E());
    Float_t xclus[3];
    clus->GetPosition(xclus);
    TVector3 cluspos(xclus[0],xclus[1],xclus[2]);
    
    fhClusterMap->Fill(cluspos.Eta(),cluspos.Phi());
  }
  
  
  
  for (Int_t itrk =  0; itrk <  ntracks; itrk++) {////////////// track loop
    iCluster = -999; //start with no match
    AliAODTrack * track = (AliAODTrack*) (GetCTSTracks()->At(itrk)) ;
    if(track->GetLabel()<0){
      if(GetDebug()>0)
        printf("Negative track label, aborting!\n");
      continue;
    }
    Double_t imp[2] = {-999.,-999.}; Double_t cov[3] = {-999.,-999.,-999.};
    Bool_t dcaOkay = GetDCA(track,imp,cov);  //homegrown dca calculation until AOD is fixed
    if(!dcaOkay&&GetDebug()>0) printf("AliAnaBtag::FillAOD - Problem computing DCA to primary vertex for track %d.  Skipping it...\n",itrk);
    fhTracks->Fill(track->Pt(),1);
    
    if(track->Pt()<0)
      continue;
    
    AliAODPid* pid = (AliAODPid*) track->GetDetPid();
    if(pid == 0) {
      if(GetDebug() > 0) printf("AliAnaBtag::MakeAnalysisFillAOD() - No PID object - skipping track %d",itrk);
      continue;
    } 
    else {
      Double_t emcpos[3];
      pid->GetEMCALPosition(emcpos);
      Double_t emcmom[3];
      pid->GetEMCALMomentum(emcmom);
      
      TVector3 pos(emcpos[0],emcpos[1],emcpos[2]);
      TVector3 mom(emcmom[0],emcmom[1],emcmom[2]);
      Double_t tphi = pos.Phi();
      Double_t teta = pos.Eta();
      Double_t tmom = mom.Mag();
      
      Bool_t in = kFALSE;
      if(track->Phi()*180./TMath::Pi() > 80. && track->Phi()*180./TMath::Pi() < 190. &&
         track->Eta() > -0.7 && track->Eta() < 0.7) in = kTRUE;
      
      
      Double_t dEdx = pid->GetTPCsignal();
      Int_t pidProb = track->GetMostProbablePID();
      Bool_t tpcEle = kFALSE; if(dEdx > 70.) tpcEle = kTRUE;
      Bool_t trkEle = kFALSE; if(pidProb == AliAODTrack::kElectron) trkEle = kTRUE;
      Bool_t emcEle = kFALSE;   
      
      
      
      
      ////////////////////////////////////////////////EMCAL//////////////////////
      if(mom.Pt() > 1.0 && in) {
        fhTracks->Fill(track->Pt(),3);
        Double_t res = 999.;
        Double_t pOverE = -999.;
        
        //Track Matching parameters
        double minRes=100.;
        Double_t minR  = 99;
        Double_t minPe =-1;
        Double_t minEp =-1;
        Double_t minMult = -1;
        Double_t minPt   = -1;
        
        for(Int_t iclus = 0; iclus < ntot; iclus++) {
          AliVCluster * clus = (AliVCluster*) (cl->At(iclus));
          if(!clus) continue;
          
          // new optimized from ben. 2010May
          if (clus->GetNCells()       < 2    ) continue;
          if (clus->GetNCells()       > 35   ) continue;
          if (clus->E()               < 0 ) continue;
          if (clus->GetDispersion()   > 1.08    ) continue;
          if (clus->GetM20()          > 0.42  ) continue;
          if (clus->GetM02()          > 0.4  ) continue;
          if (clus->GetM20()          < 0 ) continue;
          if (clus->GetM02()          < 0.06 ) continue;
          
          
          Float_t x[3];
          clus->GetPosition(x);
          TVector3 cluspos(x[0],x[1],x[2]);
          Double_t deta = teta - cluspos.Eta();
          Double_t dphi = tphi - cluspos.Phi();
          if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
          if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
          
          res = sqrt(dphi*dphi + deta*deta);
          if(res<minRes)  minRes=res;	    
          
          if(res < 0.0275) { // { //Optimized from Ben
            iCluster = iclus;
            Double_t energy = clus->E(); 
            if(energy > 0) pOverE = tmom/energy;
            if (res< minR) {
              minR  = res;
              minPe = pOverE;
              minEp = energy/tmom;
              minMult = clus->GetNCells() ;
              minPt = track->Pt();
            }
          } else {
            //unmatched
          }//res cut
          
        }//calo cluster loop
        fhResidual->Fill(minRes);
        
        if(minPe > 0.9 && minPe < 1.08) emcEle = kTRUE;//	if(minPe > fpOverEmin && minPe < fpOverEmax) emcEle = kTRUE;
      	
      }//pt, fiducial selection = EMCAL ////////////////////END EMCAL/////////////////////////
      
      
      
      
      ////////////////////////////////////////////////////Electrons/////////////////////     
      if(emcEle ||tpcEle || trkEle) { //Obsolete (kinda...)
        fhTestalle->Fill(track->Pt());
        
        //B-tagging
        if(GetDebug() > 1) printf("Found Electron - do b-tagging\n");
        Int_t pairs1=0,start=0,stop=0;
        Int_t dvmbtag = GetDVMBtag(track,pairs1,start,stop); //add: get back #pairs, start stop in pair-Ntuple.
        if(GetDebug() > 0) printf("AliAnaBtag::MakeAnalysisFillAOD -  Analyze, got back result from dvm: pair counts: pairs: %d, start %d, stop %d. \n",pairs1,start,stop);
        fhNVTX->Fill(dvmbtag);
        
        if(dvmbtag>0){
          fhDVM1->Fill(track->Pt());
          fhDVM1EtaPhi->Fill(track->Eta(),track->Phi());	    
        }
        if(dvmbtag>1)
          fhDVM2->Fill(track->Pt());
        
        //Create particle to save in AOD///////////// Purpose of this is the AODJets needs to check against this.
        Double_t eMass = 0.511/1000; //mass in GeV
        Double_t eleE = sqrt(track->P()*track->P() + eMass*eMass);
        AliAODPWG4Particle tr = AliAODPWG4Particle(track->Px(),track->Py(),track->Pz(),eleE);
        tr.SetLabel(track->GetLabel());
        tr.SetCaloLabel(iCluster,-1); //sets the indices of the original caloclusters
        tr.SetTrackLabel(track->GetID(),-1); //sets the indices of the original tracks tr.SetTrackLabel(track->GetID(),-1) instead of itrk;
        tr.SetBtag(dvmbtag);
        if(track->Charge() < 0) tr.SetPdg(11); //electron is 11
        else  tr.SetPdg(-11); //positron is -11	
        
        //Set detector flags
        Int_t emcflag=0,tpcflag=0,trdflag=0;
        if(emcEle){
          fhEmcalElectrons->Fill(track->Pt());
          emcflag=1;}
        if(trkEle){
          fhTRDElectrons->Fill(track->Pt());
          trdflag=1;}
        if(tpcEle){
          fhTPCElectrons->Fill(track->Pt());
          tpcflag=1;}
        
        if(emcEle) {//PID determined by EMCAL
          tr.SetDetector("EMCAL");
        } else {
          tr.SetDetector("CTS"); //PID determined by CTS
        }
        
        
        /////////////////////////MC stuff////////////////////////////////////////////////
        Int_t pdg = 0;
        if(IsDataMC()){
          stack=GetMCStack();
          if(!stack) {printf("AliAnaBtag::MakeAnalysisFillHistograms() - Crap, no stack: \n");
          }
          else{
            //Is it really an electron?
            TParticle *partX = stack->Particle(TMath::Abs(track->GetLabel()));
            pdg = partX->GetPdgCode();
            fhSpecies->Fill(TMath::Abs(pdg));
            if(TMath::Abs(pdg)==11){ //Check MC electrons
              if(emcEle)
                fhEmcalMCE->Fill(track->Pt());
              if(trkEle)
                fhTRDMCE->Fill(track->Pt());
              if(tpcEle)
                fhTPCMCE->Fill(track->Pt());
            }else{ //Fake histos!
              if(emcEle)
                fhEmcalMCEFake->Fill(track->Pt());
              if(trkEle)
                fhTRDMCEFake->Fill(track->Pt());
              if(tpcEle)
                fhTPCMCEFake->Fill(track->Pt());
            }
            if(TMath::Abs(pdg)==211){ //Check MC pions
              if(emcEle)
                fhEmcalMCP->Fill(track->Pt());
              if(trkEle)
                fhTRDMCP->Fill(track->Pt());
              if(tpcEle)
                fhTPCMCP->Fill(track->Pt());
            }
            if(TMath::Abs(pdg)==321){ //Check MC Kaons
              if(emcEle)
                fhEmcalMCK->Fill(track->Pt());
              if(trkEle)
                fhTRDMCK->Fill(track->Pt());
              if(tpcEle)
                fhTPCMCK->Fill(track->Pt());
            }
          }
          
          //Take care of where it came from (parent bit)
          tr.SetTag(GetMCAnalysisUtils()->CheckOrigin(tr.GetLabel(),GetReader(),tr.GetInputFileIndex())); //Gets a tag bit which contains info about super grandfather particle. Use (tag&(1<<11)), 11 for direct b, and 9 for B->C
          
          if(tr.GetTag()&(1<<9)||tr.GetTag()&(1<<11)) //MC particle from b-decay
            fhNVTXMC->Fill(dvmbtag);
          
          if(fWriteNtuple){
            fNElec++;
            fNElecEv++;
            electrons->Fill(fNElec,fEventNumber,fNElecEv,pairs1,start,stop,tpcflag,trdflag,emcflag,pdg,tr.GetTag(),tr.GetBtag(),tr.Pt());
            
            
          }
          
          if(GetDebug() > 0) 
            printf("AliAnaBtag::MakeAnalysisFillAOD() - Origin of candidate (parent bit) %d\n",tr.GetTag());
        }//MonteCarlo MC done
        
        AddAODParticle(tr);		
      }//electron
      
    } //pid check
  }//track loop                         
  if(GetDebug() > 1) printf("AliAnaBtag::MakeAnalysisFillAOD()  End fill AODs \n");  
  
}

//__________________________________________________________________
void  AliAnaBtag::MakeAnalysisFillHistograms() 
{
  //Do analysis and fill histograms
  AliStack * stack = 0x0;
  //   TParticle * primary = 0x0;
  if(IsDataMC()) {
    if(GetReader()->ReadStack()){
      stack =  GetMCStack() ;      
      if(!stack)
        printf("AliAnaBtag::MakeAnalysisFillHistograms() *** no stack ***: \n");   
    }
  }// is data and MC
  

  double maxjetEta=-4.;
  double maxjetPhi=-4.;
  
  Int_t njets = 0; 
  if(GetReader()->GetOutputEvent()) njets = (GetReader()->GetOutputEvent())->GetNJets();
  if(njets > 0) {
    if(GetDebug() > 0) printf("AliAnaBtag::MakeAnalysisFillHistograms() - Jet AOD branch has %d jets.  Performing b-jet tag analysis\n",njets);
    
    ///////////////////////////////////Jet loop//////////////////////////////////////////////
    for(Int_t ijet = 0; ijet < njets ; ijet++) {
      AliAODJet * jet = (AliAODJet*)(GetReader()->GetOutputEvent())->GetJet(ijet) ;
      fhJets->Fill(jet->Pt(),1);                                                                                ////////////////FILL
      fhJetsAllEtaPhi->Fill(jet->Eta(),jet->Phi());
      if(jet->Pt() < 0.) continue; //This has to be adjusted depending on pp or AA!
      
      //Geometric EMCAL cut
      if(TMath::Abs(jet->Eta()) > 0.3) continue;
      if(jet->Phi() < 1.8 || jet->Phi() > 2.9) continue; //This is BAD FIXME 
      fhJets->Fill(jet->Pt(),4); //All jets after geometric cut                                                 ////////////////FILL
      
      Bool_t leadJet  = kFALSE;
      if (ijet==0) {leadJet= kTRUE; fhJets->Fill(jet->Pt(),5);} //Leading jets                                   ////////////////FILL
      
      /////////////////////////Track loop in Jet////////////////////////////////////
      Bool_t dvmJet = kFALSE; 
      Bool_t dvmMCJet = kFALSE; 
      TRefArray* rt = jet->GetRefTracks();
      Int_t ntrk = rt->GetEntries();
      
      for(Int_t itrk = 0; itrk < ntrk; itrk++) {
      	AliAODTrack* jetTrack = (AliAODTrack*)jet->GetTrack(itrk);
	
	Int_t trackId = jetTrack->GetID(); //get the index in the reader
	Int_t naod = GetOutputAODBranch()->GetEntriesFast();
	if(GetDebug() > 3) printf("AliAnaBtag::CheckIfBjet() - aod branch entries %d\n", naod);
	
	for(Int_t iaod = 0; iaod < naod ; iaod++){
	  AliAODPWG4Particle* ele =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
	  Int_t electronlabel = ele->GetTrackLabel(0);
	  if(electronlabel != trackId) continue;  //skip to the next one if they don't match
	  if(ele->GetBtag()>0){
	    dvmJet = kTRUE;
	    if(ele->GetTag()&(1<<9) || ele->GetTag()&(1<<11) )
	      dvmMCJet = kTRUE;
	  }
	} //Electron check
	
      }//Track loop of jet tracks
      
      if(dvmJet) fhJets->Fill(jet->Pt(),6);                                                                      ////////////////FILL
      //MC stuff
      if(IsDataMC()) {
        Bool_t isTrueBjet = IsMcBJet(jet->Eta(), jet->Phi());
        //Bool_t isTrueDjet = IsMcDJet(jet->Eta(), jet->Phi());
	if(dvmJet){
	  if(dvmMCJet){
	    fhJets->Fill(jet->Pt(),8);   //True                                                                   ////////////////FILL
	  }else{
	    fhJets->Fill(jet->Pt(),9);   //False                                                                 ////////////////FILL
	  }
	  if(isTrueBjet){ 
	    fhJets->Fill(jet->Pt(),10);     //True                                                                 ////////////////FILL
	  }else{
	    fhJets->Fill(jet->Pt(),11);     //False                                                                ////////////////FILL
	  }
	}
	if(isTrueBjet){ 
	  fhJets->Fill(jet->Pt(),12);     //True                                                                 ////////////////FILL
	}else{
	  fhJets->Fill(jet->Pt(),13);     //False                                                                ////////////////FILL
	}

      }//MC stuff

    }//jet loop
  } //jets exist
  
  //Electron loop, read back electrons, fill histos; mostly photonic shit.
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaBtag::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* ele =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t pdg = ele->GetPdg();
    
    
    if(TMath::Abs(pdg) != AliCaloPID::kElectron) continue; //not necessary..
    
    //MC tag of this electron
    //    Int_t mctag = ele->GetTag();
    
    
    fhElectrons->Fill(ele->Pt(),1); //All electrons
    Bool_t photonic = kFALSE;
    photonic = PhotonicV0(ele->GetTrackLabel(0)); //check against V0s
    if(!photonic) fhElectrons->Fill(ele->Pt(),3); //nonphotonic electrons
    if(photonic) fhElectrons->Fill(ele->Pt(),4);  //photonic electrons
    
    //Fill electron histograms 
    Float_t phiele = ele->Phi();
    Float_t etaele = ele->Eta();
    
    
    if(ele->GetBtag()>0){ // removed bit tag shit
      fhElectrons->Fill(ele->Pt(),5);
      if(!photonic) fhElectrons->Fill(ele->Pt(),6);
      if(photonic) fhElectrons->Fill(ele->Pt(),7);
      fhJetsLeadingBElectronEtaPhi->Fill(maxjetEta,maxjetPhi); 
      double deta=etaele-maxjetEta;
      double dphi=phiele-maxjetPhi;
      
      //double r = sqrt((dphi*dphi)+(deta*deta));
      fhBJetElectronDetaDphi->Fill(deta,dphi);
      
    }
    
  }//electron aod loop
  
}

//__________________________________________________________________
Int_t AliAnaBtag::GetDVMBtag(AliAODTrack * tr, Int_t &pair, Int_t &start, Int_t &stop)
{
  //This method uses the Displaced Vertex between electron-hadron
  //pairs and the primary vertex to determine whether an electron is
  //likely from a B hadron.
  Int_t pairstart = fNPair;
  Int_t pairstop = 0;
  Int_t pairn = 0;
  Int_t ncls1 = 0;
  for(Int_t l = 0; l < 6; l++) if(TESTBIT(tr->GetITSClusterMap(),l)) ncls1++;

  if (ncls1 < fITSCut) return 0;

  Double_t imp[2] = {-999.,-999.}; Double_t cov[3] = {-999.,-999.,-999.};
  Bool_t dcaOkay = GetDCA(tr,imp,cov);  //homegrown dca calculation until AOD is fixed                  
  if(!dcaOkay) {
    printf("AliAnaBtag::GetDVMBtag - Problem computing DCA to primary vertex for track %d",tr->GetID());
    return 0;
  }

  if (TMath::Abs(imp[0])   > fImpactCut ) return 0;
  if (TMath::Abs(imp[1])   > fImpactCut ) return 0;

  Int_t nvtx = 0;
//   Int_t nvtx2 = 0;
//   Int_t nvtx3 = 0;

  for (Int_t k2 =0; k2 < GetCTSTracks()->GetEntriesFast() ; k2++) {

    //loop over assoc
    AliAODTrack* track2 = (AliAODTrack*)GetCTSTracks()->At(k2);
    Int_t id1 = tr->GetID();
    Int_t id2 = track2->GetID();
    if(id1 == id2) continue;

    Int_t ncls2 = 0;
    for(Int_t l = 0; l < 6; l++) if(TESTBIT(track2->GetITSClusterMap(),l)) ncls2++;
    if (ncls2 < fITSCut) continue;


    Double_t sdca=0,pdca=0,minv=0,pth=0,massphoton=0,decaylength=0;

    sdca = ComputeSignDca(tr,track2,minv,pdca,massphoton,decaylength);
    pth=track2->Pt();




    Double_t dphi = tr->Phi() - track2->Phi();
    if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
    Double_t deta = tr->Eta() - track2->Eta();
    Double_t dr = sqrt(deta*deta + dphi*dphi);

    if(dr > fDrCut) continue;
    fNPair++;
    pairn++;
    if(GetDebug() > 0) 
      printf("pairs: %d !!!!!!!!!!! \n",fNPair);
    if(fWriteNtuple){
      pairs->Fill(fNPair,fNElec,fEventNumber,pdca,sdca,minv,pth,massphoton,decaylength);
    }
    pairstop=fNPair;
    //pairs->Fill(1,1,1,1,1,1,1);
    fhPairPt->Fill(pth);
    if(decaylength>1.0) continue;
    if(tr->Pt()<6. && pth>0.4 && minv>1.4 && pdca<0.025 && sdca>0.06 && massphoton>0.1) nvtx++; 
    if(tr->Pt()>6.&&tr->Pt()<10. && pth>0.2 && minv>1.7 && pdca<0.012 && sdca>0.06 && massphoton>0.1) nvtx++;
    if(tr->Pt()>10.&& pth>0.6 && minv>1.5 && pdca<0.14 && sdca>0.04 && massphoton>0.1) nvtx++;



  } //loop over hadrons

  if(GetDebug() > 0) printf("AliAnaBtag::GetDVMBtag - result of btagging: %d \n",nvtx);


  pair=pairn;
  start=pairstart+1;
  stop=pairstop;
  if(GetDebug() > 0) printf("End of DVM, pair counts: pairs: %d, start %d, stop %d. \n",pair,start,stop);
  return nvtx;
}

//__________________________________________________________________
Double_t AliAnaBtag::ComputeSignDca(AliAODTrack *tr, AliAODTrack *tr2 , Double_t &masscut, Double_t &pdcacut, Double_t &massphoton, Double_t &decay)
{
  //Compute the signed dca between two tracks
  //and return the result

  Double_t signDca=-999.;
  if(GetDebug() > 2 ) printf(">>ComputeSdca:: track1 %d, track2 %d, masscut %f \n", tr->GetLabel(), tr2->GetLabel(), masscut);

  //=====Now calculate DCA between both tracks=======  
  Double_t massE = 0.000511;
  Double_t massK = 0.493677;
  Double_t vertex[3] = {-999.,-999.,-999}; //vertex

  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) {
    GetVertex(vertex); //If only one file, get the vertex from there  
  }
  
  TVector3 primV(vertex[0],vertex[1],vertex[2]) ; 
  if(GetDebug()>0) printf(">>ComputeSdca:: primary vertex = %2.2f,%2.2f,%2.2f \n",vertex[0],vertex[1],vertex[2]) ;

  AliExternalTrackParam *param1 = new AliExternalTrackParam(tr);
  AliExternalTrackParam *param2 = new AliExternalTrackParam(tr2);

  Double_t bfield[3];
  param1->GetBxByBz(bfield);
  bfield[0]=0;
  bfield[1]=0;
  bfield[2]=5.;
  Double_t bz = 5.; // = param1->GetBz();
  Double_t xplane1 = 0.; Double_t xplane2 = 0.;
  Double_t pairdca = param1->GetDCA(param2,bz,xplane1,xplane2);

  param1->PropagateToBxByBz(xplane1,bfield);
  param2->PropagateToBxByBz(xplane2,bfield);

  Int_t id1 = 0, id2 = 0;
  AliESDv0 bvertex(*param1,id1,*param2,id2);
  Double_t vx,vy,vz;
  bvertex.GetXYZ(vx,vy,vz);

  Double_t emom[3];
  Double_t hmom[3];
  param1->PxPyPz(emom);
  param2->PxPyPz(hmom);
  TVector3 emomAtB(emom[0],emom[1],emom[2]);
  TVector3 hmomAtB(hmom[0],hmom[1],hmom[2]);
  TVector3 secvtxpt(vx,vy,vz);
  TVector3 decayvector(0,0,0);
  decayvector = secvtxpt - primV; //decay vector from PrimVtx
  Double_t decaylength = decayvector.Mag();
  decay=decaylength;

  if(GetDebug() > 0) {
    printf(">>ComputeSdca:: mom1=%f, mom2=%f \n", emomAtB.Perp(), hmomAtB.Perp() );
    printf(">>ComputeSdca:: pairDCA=%f, length=%f \n", pairdca,decaylength );
  }

  if (emomAtB.Mag()>0 /*&& decaylength < fDecayLenCut*/ ) {
    TVector3 sumMom = emomAtB+hmomAtB;
    Double_t ener1 = sqrt(pow(emomAtB.Mag(),2) + massE*massE);
    Double_t ener2 = sqrt(pow(hmomAtB.Mag(),2) + massK*massK);
    Double_t ener3 = sqrt(pow(hmomAtB.Mag(),2) + massE*massE);
    Double_t mass = sqrt(pow((ener1+ener2),2) - pow(sumMom.Mag(),2));
    Double_t massPhot = sqrt(pow((ener1+ener3),2) - pow(sumMom.Mag(),2));
    Double_t sDca = decayvector.Dot(emomAtB)/emomAtB.Mag();
    pdcacut=pairdca;//
    masscut=mass; // 
    massphoton=massPhot;//                                                                     Send it back!
    signDca = sDca;
    if(GetDebug() > 0) printf("\t>>ComputeSdca:: mass=%f \n", mass);
    if(GetDebug() > 0) printf("\t>>ComputeSdca:: sec vtx-signdca :%f\n",signDca);
  }
  //clean up
  delete param1;
  delete param2;
  return signDca;
}


//__________________________________________________________________
Bool_t AliAnaBtag::PhotonicV0(Int_t id) 
{
  //This method checks to see whether a track that has been flagged as
  //an electron was determined to match to a V0 candidate with
  //invariant mass consistent with photon conversion
  Bool_t itIS = kFALSE;
  Double_t massEta = 0.547;
  Double_t massRho0 = 0.770;
  Double_t massOmega = 0.782;
  Double_t massPhi = 1.020;
  
  //---Get V0s---
  AliAODEvent *aod = (AliAODEvent*) GetReader()->GetInputEvent();
  int nv0s = aod->GetNumberOfV0s();
  for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {
    AliAODv0 *v0 = aod->GetV0(iV0);
    if (!v0) continue;
    double radius = v0->RadiusV0();
    double mass = v0->InvMass2Prongs(0,1,11,11);
    if(GetDebug() > 0) {
      printf("## PhotonicV0() :: v0: %d, radius: %f \n", iV0 , radius );
      printf("## PhotonicV0() :: neg-id: %d, pos-id: %d, THIS id: %d\n", v0->GetNegID(), v0->GetPosID(), id);
      printf("## PhotonicV0() :: Minv(e,e): %f \n", v0->InvMass2Prongs(0,1,11,11) );
    }
    if (mass < 0.100 ||
	(mass > massEta-0.05 || mass < massEta+0.05) ||
	(mass > massRho0-0.05 || mass < massRho0+0.05) ||
	(mass > massOmega-0.05 || mass < massOmega+0.05) ||
	(mass > massPhi-0.05 || mass < massPhi+0.05)) {
      if ( id == v0->GetNegID() || id == v0->GetPosID()) {
	itIS=kTRUE;
	if(GetDebug() > 0) printf("## PhotonicV0() :: It's a conversion electron!!! \n" );
      }
    } }
  return itIS;
}

//__________________________________________________________________
Bool_t AliAnaBtag::GetDCA(const AliAODTrack* track,Double_t impPar[2], Double_t cov[3]) 
{
  //Use the Event vertex and AOD track information to get
  //a real impact parameter for the track
  Double_t maxD = 100000.; //max transverse IP
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) {
    AliVEvent* ve = (AliVEvent*)GetReader()->GetInputEvent();
    AliVVertex *vv = (AliVVertex*)ve->GetPrimaryVertex();
    AliESDtrack esdTrack(track);
    Double_t bfield[3];
    esdTrack.GetBxByBz(bfield);
    bfield[0]=0;
    bfield[1]=0;
    bfield[2]=5.;
    
    Bool_t gotit = esdTrack.PropagateToDCABxByBz(vv,bfield,maxD,impPar,cov);
    return gotit;
  }
  return kFALSE;
}
//__________________________________________________________________
Bool_t AliAnaBtag::CheckIfBjet(const AliAODTrack* track)
{
  Bool_t bjet = kFALSE;
  Int_t trackId = track->GetID(); //get the index in the reader
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 3) printf("AliAnaBtag::CheckIfBjet() - aod branch entries %d\n", naod);

  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* ele =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t electronlabel = ele->GetTrackLabel(0);
    if(electronlabel != trackId) continue;  //skip to the next one if they don't match
    if(ele->GetBtag()>0)
      bjet = kTRUE;
  } 

  return bjet;
} 

//__________________________________________________________________
AliAODMCParticle* AliAnaBtag::GetMCParticle(Int_t ipart) 
{
  //Get the MC particle at position ipart
  
  AliAODMCParticle* aodprimary = 0x0;
  TClonesArray * mcparticles0 = 0x0;
  
  if(GetReader()->ReadAODMCParticles()){
    //Get the list of MC particles                                                                                                                           
    mcparticles0 = GetReader()->GetAODMCParticles(0);
    if(!mcparticles0) {
      if(GetDebug() > 0)printf("AliAnaBtag::MakeAnalysisFillHistograms() -  Standard MCParticles not available!\n");
    }
    else{
      Int_t npart0 = mcparticles0->GetEntriesFast();
      if(ipart < npart0) aodprimary = (AliAODMCParticle*)mcparticles0->At(ipart);
      
      if(!aodprimary) {
        printf("AliAnaBtag::GetMCParticle() *** no primary ***:  label %d \n", ipart);
        return 0x0;
      }
    }
  } else {
    printf("AliAnaBtag::GetMCParticle() - Asked for AliAODMCParticle but we have a stack reader.\n");
  }
  return aodprimary;
}

//__________________________________________________________________
Bool_t  AliAnaBtag::IsMcBJet(Double_t jeta, Double_t jphi)
{
  //Check the jet eta,phi against that of the b-quark
  //to decide whether it is an MC B-jet
  Bool_t bjet=kFALSE;
  AliStack* stack = 0x0;
  
  for(Int_t ipart = 0; ipart < 100; ipart++) {

    Double_t pphi = -999.;
    Double_t peta = -999.;
    Int_t pdg = 0;
    if(GetReader()->ReadStack()) {
      stack = GetMCStack();
      if(!stack) {
	printf("AliAnaBtag::IsMCBJet() *** no stack ***: \n");
	return kFALSE;
      }
      TParticle* primary = stack->Particle(ipart);
      if (!primary) continue;
      pdg = primary->GetPdgCode();
      pphi = primary->Phi();
      peta = primary->Eta();
    } else if(GetReader()->ReadAODMCParticles()) {
      AliAODMCParticle* aodprimary = GetMCParticle(ipart);
      if(!aodprimary) continue;
      pdg = aodprimary->GetPdgCode();
      pphi = aodprimary->Phi();
      peta = aodprimary->Eta();
    }
    if ( TMath::Abs(pdg) != 5) continue;
    Double_t dphi = jphi - pphi;
    Double_t deta = jeta - peta;
    Double_t dr = sqrt(deta*deta + dphi*dphi);
    
    if (dr < 0.2) {
      bjet=kTRUE;
      break;
    }
  }
  return bjet;
}

//__________________________________________________________________
Bool_t  AliAnaBtag::IsMcDJet(Double_t jeta, Double_t jphi)
{
  //Check if this jet is a charm jet
  Bool_t cjet=kFALSE;

  AliStack* stack = 0x0;

  for(Int_t ipart = 0; ipart < 100; ipart++) {
    
    Double_t pphi = -999.;
    Double_t peta = -999.;
    Int_t pdg = 0;
    if(GetReader()->ReadStack()) {
      stack = GetMCStack();
      if(!stack) {
	printf("AliAnaBtag::IsMCDJet() *** no stack ***: \n");
	return kFALSE;
      }
      TParticle* primary = stack->Particle(ipart);
      if (!primary) continue;
      pdg = primary->GetPdgCode();
      pphi = primary->Phi();
      peta = primary->Eta();
    } else if(GetReader()->ReadAODMCParticles()) {
      AliAODMCParticle* aodprimary = GetMCParticle(ipart);
      if(!aodprimary) continue;
      pdg = aodprimary->GetPdgCode();
      pphi = aodprimary->Phi();
      peta = aodprimary->Eta();
    }

    if ( TMath::Abs(pdg) != 4) continue;
    Double_t dphi = jphi - pphi;
    Double_t deta = jeta - peta;
    Double_t dr = sqrt(deta*deta + dphi*dphi);
   
    if (dr < 0.2) {
      cjet=kTRUE;
      break;
    }
  }
  return cjet;
}


//__________________________________________________________________
void  AliAnaBtag::Terminate(TList* outputList)
{
  printf(" AliAnaBtag::Terminate()  *** %s Report: %d outputs/histograms \n", GetName(), outputList->GetEntries()) ;
}
