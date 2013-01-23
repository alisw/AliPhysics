// $Id$

//////////
//Measure Jet-hadron correlations
//Does event Mixing using AliEventPoolManager
/////////

#include "AliAnalysisTaskEmcalJetHMEC.h"

#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include <TClonesArray.h>
#include <TParticle.h>
#include "AliVTrack.h"
#include "TParameter.h"

#include "AliAODEvent.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDCaloCluster.h"
#include "AliESDVertex.h"
#include "AliCentrality.h"
#include "AliAODJet.h"
#include "AliEmcalJet.h"
#include "AliESDtrackCuts.h"

#include "TVector3.h"
#include "AliPicoTrack.h"
#include "AliEventPoolManager.h"

ClassImp(AliAnalysisTaskEmcalJetHMEC)

//________________________________________________________________________
AliAnalysisTaskEmcalJetHMEC::AliAnalysisTaskEmcalJetHMEC() : 
  AliAnalysisTaskSE(),
  fTracksName("tracks"),
  fJetsName("jets"),
  fPhimin(-10), 
  fPhimax(10),
  fEtamin(-0.9), 
  fEtamax(0.9),
  fAreacut(0.0),
  fTrkBias(5),
  fClusBias(5),
  fTrkEta(0.9),
  fDoEventMixing(0),
  fMixingTracks(50000),
  fESD(0), 
  fPoolMgr(0x0), 
  fOutputList(0),
  fHistTrackPt(0),
  fHistCentrality(0), 
  fHistJetEtaPhi(0), 
  fHistJetHEtaPhi(0), 
  fhnMixedEvents(0x0),
  fhnJH(0x0)
{
  // Default Constructor

  for(Int_t ipta=0; ipta<7; ipta++){
    fHistTrackEtaPhi[ipta]=0;
  }


  for(Int_t icent = 0; icent<6; ++icent){
    fHistJetPt[icent]=0;
    fHistJetPtBias[icent]=0;
    fHistLeadJetPt[icent]=0;
    fHistLeadJetPtBias[icent]=0;
    fHistJetPtTT[icent]=0;
    for(Int_t iptjet = 0; iptjet<5; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){	
	fHistJetH[icent][iptjet][ieta]=0;
	fHistJetHBias[icent][iptjet][ieta]=0;
	fHistJetHTT[icent][iptjet][ieta]=0;

      }
    }
  }

}
//________________________________________________________________________
AliAnalysisTaskEmcalJetHMEC::AliAnalysisTaskEmcalJetHMEC(const char *name) : 
  AliAnalysisTaskSE(name),
  fTracksName("tracks"),
  fJetsName("jets"),
  fPhimin(-10), 
  fPhimax(10),
  fEtamin(-0.9), 
  fEtamax(0.9),
  fAreacut(0.0),
  fTrkBias(5),
  fClusBias(5),
  fTrkEta(0.9),
  fDoEventMixing(0),
  fMixingTracks(50000),
  fESD(0), 
  fPoolMgr(0x0), 
  fOutputList(0), 
  fHistTrackPt(0),
  fHistCentrality(0), 
  fHistJetEtaPhi(0), 
  fHistJetHEtaPhi(0),
  fhnMixedEvents(0x0),
  fhnJH(0x0)
{
  // Constructor
  for(Int_t ipta=0; ipta<7; ipta++){
    fHistTrackEtaPhi[ipta]=0;
  }
  for(Int_t icent = 0; icent<6; ++icent){
    fHistJetPt[icent]=0;
    fHistJetPtBias[icent]=0;
    fHistLeadJetPt[icent]=0;
    fHistLeadJetPtBias[icent]=0;
    fHistJetPtTT[icent]=0;
    for(Int_t iptjet = 0; iptjet<5; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){	
	fHistJetH[icent][iptjet][ieta]=0;
	fHistJetHBias[icent][iptjet][ieta]=0;
	fHistJetHTT[icent][iptjet][ieta]=0;
      }
    }
  }


  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());

}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetHMEC::UserCreateOutputObjects()
{
  // Called once

 
  AliVEventHandler* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!handler) {
    AliError("Input handler not available!");
    return;
  }

  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner();


  // Create histograms
  fHistTrackPt = new TH1F("fHistTrackPt", "P_{T} distribution", 1000, 0.0, 100.0);

  fHistCentrality = new TH1F("fHistCentrality","centrality",100,0,100);

  fHistJetEtaPhi = new TH2F("fHistJetEtaPhi","Jet eta-phi",900,-1.8,1.8,720,-3.2,3.2);
  fHistJetHEtaPhi = new TH2F("fHistJetHEtaPhi","Jet-Hadron deta-dphi",900,-1.8,1.8,720,-1.6,4.8);

  TString name;

  for(Int_t ipta=0; ipta<7; ++ipta){
    name = Form("fHistTrackEtaPhi_%i", ipta);
    fHistTrackEtaPhi[ipta] = new TH2F(name,name,400,-1,1,720,0.0,2.0*TMath::Pi());
    fOutputList->Add(fHistTrackEtaPhi[ipta]);

  }
 
  for(Int_t icent = 0; icent<6; ++icent){
    name = Form("fHistJetPt_%i",icent);   
    fHistJetPt[icent] = new TH1F(name,name,200,0,200);
    fOutputList->Add(fHistJetPt[icent]);

    name = Form("fHistJetPtBias_%i",icent);   
    fHistJetPtBias[icent] = new TH1F(name,name,200,0,200);
    fOutputList->Add(fHistJetPtBias[icent]);

    name = Form("fHistLeadJetPt_%i",icent);   
    fHistLeadJetPt[icent] = new TH1F(name,name,200,0,200);
    fOutputList->Add(fHistLeadJetPt[icent]);

    name = Form("fHistLeadJetPtBias_%i",icent);   
    fHistLeadJetPtBias[icent] = new TH1F(name,name,200,0,200);
    fOutputList->Add(fHistLeadJetPtBias[icent]);

    name = Form("fHistJetPtTT_%i",icent);   
    fHistJetPtTT[icent] = new TH1F(name,name,200,0,200);
    fOutputList->Add(fHistJetPtTT[icent]);

    for(Int_t iptjet = 0; iptjet<5; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){	
	name = Form("fHistJetH_%i_%i_%i",icent,iptjet,ieta);   
	fHistJetH[icent][iptjet][ieta]=new TH2F(name,name,72,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
	fOutputList->Add(fHistJetH[icent][iptjet][ieta]);

	name = Form("fHistJetHBias_%i_%i_%i",icent,iptjet,ieta);   
	fHistJetHBias[icent][iptjet][ieta]=new TH2F(name,name,72,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
	fOutputList->Add(fHistJetHBias[icent][iptjet][ieta]);

	name = Form("fHistJetHTT_%i_%i_%i",icent,iptjet,ieta);   
	fHistJetHTT[icent][iptjet][ieta]=new TH2F(name,name,72,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
	fOutputList->Add(fHistJetHTT[icent][iptjet][ieta]);

      }
    }
  }



  UInt_t cifras = 0; // bit coded, see GetDimParams() below 
  cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<7 | 1<<8; 
  fhnJH = NewTHnSparseF("fhnJH", cifras);
  
  fhnJH->Sumw2();

  fOutputList->Add(fhnJH);


  if(fDoEventMixing){    
    cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<7 | 1<<8; 
    fhnMixedEvents = NewTHnSparseF("fhnMixedEvents", cifras);
    
    fhnMixedEvents->Sumw2();
    
    fOutputList->Add(fhnMixedEvents);
    
  }
  

  fOutputList->Add(fHistTrackPt);
  fOutputList->Add(fHistCentrality);
  fOutputList->Add(fHistJetEtaPhi);
  fOutputList->Add(fHistJetHEtaPhi);


  PostData(1, fOutputList);


  //Event Mixing
  Int_t trackDepth = fMixingTracks; 
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
 
  Int_t nZvtxBins  = 5+1+5;
  // bins for second buffer are shifted by 100 cm
  Double_t vertexBins[] = { -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, };
  Double_t* zvtxbin = vertexBins;

  Int_t nCentralityBins  = 100;
  Double_t centralityBins[nCentralityBins];
  for(Int_t ic=0; ic<nCentralityBins; ic++){
    centralityBins[ic]=1.0*ic;
  }

  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centralityBins, nZvtxBins, zvtxbin);

  

}

//________________________________________________________________________

Double_t AliAnalysisTaskEmcalJetHMEC:: RelativePhi(Double_t mphi,Double_t vphi) {


  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  Double_t dphi = vphi-mphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());

  return dphi;//dphi in [-Pi, Pi]                                                                                                    
}


//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHMEC::GetCentBin(Double_t cent) const 
{
  // Get centrality bin.

  Int_t centbin = -1;
  if (cent>=0 && cent<10)
    centbin = 0;
  else if (cent>=10 && cent<20)
    centbin = 1;
  else if (cent>=20 && cent<30)
    centbin = 2;
  else if (cent>=30 && cent<40)
    centbin = 3;
  else if (cent>=40 && cent<50)
    centbin = 4;
  else if (cent>=50 && cent<90)
    centbin = 5;
  return centbin;
}


//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHMEC::GetEtaBin(Double_t eta) const 
{
  // Get eta bin for histos.

  Int_t etabin = -1;
  if (TMath::Abs(eta)<=0.4)
    etabin = 0;
  else if (TMath::Abs(eta)>0.4 && TMath::Abs(eta)<0.8)
    etabin = 1;
  else if (TMath::Abs(eta)>=0.8)
    etabin = 2;
  return etabin;
}
//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHMEC::GetpTjetBin(Double_t pt) const 
{
  // Get jet pt  bin for histos.

  Int_t ptbin = -1;
  if (pt>=15 && pt<20)
    ptbin = 0;
  else if (pt>=20 && pt<25)
    ptbin = 1;
  else if (pt>=25 && pt<30)
    ptbin = 2;
  else if (pt>=30 && pt<60)
    ptbin = 3;
  else if (pt>=60)
    ptbin = 4;


  return ptbin;
}


//________________________________________________________________________
void AliAnalysisTaskEmcalJetHMEC::UserExec(Option_t *) 
{


  // Main loop called for each event
 // esd or aod mode
  Bool_t esdMode = kTRUE;
  if (dynamic_cast<AliAODEvent*>(InputEvent()))
    esdMode = kFALSE;


  if (esdMode) {
    // optimization in case autobranch loading is off
    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    if (fTracksName == "Tracks")
      am->LoadBranch("Tracks");
  }
 

  //get centrality
  TList *list = InputEvent()->GetList();
  AliCentrality *centrality = InputEvent()->GetCentrality() ;
  Double_t fcent=-1; 
  if(centrality)
    fcent = centrality->GetCentralityPercentile("V0M");
  else
    fcent=99;//probably pp data

  if (fcent<0) {
    AliError(Form("Centrality negative: %f", fcent));
    return;
  }

  Double_t fvertex[3]={0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(fvertex);
  Double_t zVtx=fvertex[2];

  if(fabs(zVtx)>10.0)
    return;

  fHistCentrality->Fill(fcent);
  Int_t centbin = GetCentBin(fcent);

  if(centbin<0)
    return;
    
  TClonesArray *jets = 0;
  TClonesArray *tracks = 0;

  tracks = dynamic_cast<TClonesArray*>(list->FindObject(fTracksName));
  if (!tracks) {
    AliError(Form("Pointer to tracks %s == 0", fTracksName.Data() ));
    return;
  }
  const Int_t Ntracks=tracks->GetEntries();

  jets= dynamic_cast<TClonesArray*>(list->FindObject(fJetsName));
  if (!jets) {
    AliError(Form("Pointer to tracks %s == 0", fJetsName.Data() ));
    return;
  }

  const Int_t Njets = jets->GetEntries();
 
  //Leticia's loop to find hardest track

  Int_t iTT=-1;
  Double_t ptmax=-10;

  for (Int_t iTracks = 0; iTracks < Ntracks; iTracks++) 
    {
      AliVTrack* track = static_cast<AliVTrack*>(tracks->At(iTracks));
      if (!track) {
	printf("ERROR: Could not receive track %d\n", iTracks);
	continue;
      }
      
      if(TMath::Abs(track->Eta())>0.9) continue;
      if(track->Pt()<0.15)continue;
      //iCount++;
      if(track->Pt()>ptmax){
	ptmax=track->Pt();
	iTT=iTracks;
      }
    }

  
  Int_t ijethi=-1;

  Double_t highestjetpt=0.0;

  Int_t passedTTcut=0;

  for (Int_t ijet = 0; ijet < Njets; ijet++)
    {
 
      AliEmcalJet *jet = static_cast<AliEmcalJet*>(jets->At(ijet));
      
      if (!jet)
	continue;
      
      if(!AcceptJet(jet))
	continue;

      Double_t jetPt = jet->Pt();

      if(highestjetpt<jetPt){
	ijethi=ijet;
	highestjetpt=jetPt;
      }

    }


    
  for (Int_t ijet = 0; ijet < Njets; ijet++){
    
    AliEmcalJet *jet = static_cast<AliEmcalJet*>(jets->At(ijet));

    if(!AcceptJet(jet))
      continue;

        
    Double_t jetphi = jet->Phi();
    Double_t jetPt = jet->Pt();
    Double_t jeteta=jet->Eta();

    Double_t leadjet=0;
    if (ijet==ijethi) leadjet=1;


    fHistJetPt[centbin]->Fill(jet->Pt());
    fHistLeadJetPt[centbin]->Fill(jet->Pt());
    
    if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fClusBias)){
      fHistJetPtBias[centbin]->Fill(jet->Pt());
      fHistLeadJetPtBias[centbin]->Fill(jet->Pt());
    }

      fHistJetEtaPhi->Fill(jet->Eta(),jetphi);


      if(iTT>0){
	AliVTrack* TT = static_cast<AliVTrack*>(tracks->At(iTT));
	if(TMath::Abs(jetphi-TT->Phi()-TMath::Pi())<0.6) passedTTcut=1;
	else passedTTcut=0;
      }

      if(passedTTcut)
	fHistJetPtTT[centbin]->Fill(jet->Pt());
   
  
      Int_t iptjet=-1;
      iptjet=GetpTjetBin(jetPt);
      if(iptjet<0) continue;


   if (highestjetpt>15) {

  
    for (Int_t iTracks = 0; iTracks < Ntracks; iTracks++) 
      {
	AliVTrack* track = static_cast<AliVTrack*>(tracks->At(iTracks));
	if (!track) {
	  printf("ERROR: Could not receive track %d\n", iTracks);
	  continue;
	}
	
	  if(TMath::Abs(track->Eta())>fTrkEta) continue;

	  fHistTrackPt->Fill(track->Pt());
	  	  
	  if (track->Pt()<0.15)
	    continue;
	  
	  Double_t trackphi = track->Phi();
	  if (trackphi > TMath::Pi())
	    trackphi = trackphi-2*TMath::Pi();

	  Double_t tracketa=track->Eta();
	  Double_t trackpt=track->Pt();
	  Double_t deta=tracketa-jeteta;
	  Int_t ieta=GetEtaBin(deta);
	  if (ieta<0) {
	    AliError(Form("Eta Bin negative: %f", deta));
	    continue;
	  }
      

	  //Jet pt, track pt, dPhi,deta,fcent
	  Double_t dphijh = RelativePhi(jetphi,trackphi);
	  if (dphijh < -0.5*TMath::Pi())
	    dphijh+= 2*TMath::Pi();
	  if (dphijh > 1.5*TMath::Pi()) dphijh-=2.*TMath::Pi();


      
	  fHistJetH[centbin][iptjet][ieta]->Fill(dphijh,track->Pt());
	  fHistJetHEtaPhi->Fill(deta,dphijh);

	  Double_t dR=sqrt(deta*deta+dphijh*dphijh);

	  if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fClusBias)){
	    fHistJetHBias[centbin][iptjet][ieta]->Fill(dphijh,trackpt);


	    Double_t triggerEntries[8] = {fcent,jetPt,trackpt,dR,deta,dphijh,0.0,leadjet};                      
	    fhnJH->Fill(triggerEntries);
	  }

	  if(passedTTcut)
	    fHistJetHTT[centbin][iptjet][ieta]->Fill(dphijh,trackpt);


	} //track loop
  }//jet pt cut

  }//jet loop
  

  //Prepare to do event mixing

  // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
  TObjArray* tracksClone = CloneAndReduceTrackList(tracks);
  //delete tracks;


  if(fDoEventMixing>0){
    
    // event mixing
    
    // 1. First get an event pool corresponding in mult (cent) and
    //    zvertex to the current event. Once initialized, the pool
    //    should contain nMix (reduced) events. This routine does not
    //    pre-scan the chain. The first several events of every chain
    //    will be skipped until the needed pools are filled to the
    //    specified depth. If the pool categories are not too rare, this
    //    should not be a problem. If they are rare, you could lose
    //    statistics.

    // 2. Collect the whole pool's content of tracks into one TObjArray
    //    (bgTracks), which is effectively a single background super-event.

    // 3. The reduced and bgTracks arrays must both be passed into
    //    FillCorrelations(). Also nMix should be passed in, so a weight
    //    of 1./nMix can be applied.





    AliEventPool* pool = fPoolMgr->GetEventPool(fcent, zVtx);
    
    if (!pool){
      AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fcent, zVtx));
      return;
    }

    //check for a trigger jet
        

    if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 || pool->GetCurrentNEvents() >= 5) 
      {
	
	for (Int_t ijet = 0; ijet < Njets; ijet++){

	  Double_t leadjet=0;
	  if (ijet==ijethi) leadjet=1;
	  
	  AliEmcalJet *jet = static_cast<AliEmcalJet*>(jets->At(ijet));

	  if(!AcceptJet(jet))
	    continue;

	  Double_t jetPt = jet->Pt();	
	  Double_t jetphi = jet->Phi();
	  Double_t jeteta=jet->Eta();
	  
	  
	  Int_t nMix = pool->GetCurrentNEvents();
	  
	  //Fill for biased jet triggers only
	  if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fClusBias)){

	    // Fill mixed-event histos here  
	    for (Int_t jMix=0; jMix<nMix; jMix++) 
	      {
		TObjArray* bgTracks = pool->GetEvent(jMix);
		const Int_t Nbgtrks = bgTracks->GetEntries();
		for(Int_t ibg=0; ibg<Nbgtrks; ibg++){
		  AliPicoTrack *part = static_cast<AliPicoTrack*>(bgTracks->At(ibg));         
		  if(!part) continue;
  
		  Double_t DEta = part->Eta()-jeteta;
		  Double_t DPhi = RelativePhi(jetphi,part->Phi());

		  Double_t DR=TMath::Sqrt(DPhi*DPhi+DEta*DEta);
		  if(DPhi<-0.5*TMath::Pi()) DPhi+=2.*TMath::Pi();
		  if(DPhi>3./2.*TMath::Pi()) DPhi-=2.*TMath::Pi();
		  Double_t triggerEntries[8] = {fcent,jetPt,part->Pt(),DR,DEta,DPhi,0.0,leadjet};                      
		  fhnMixedEvents->Fill(triggerEntries,1./nMix);
		  
		  
		}
	      }
	  }
	}
    }

    //update pool if jet in event or not
    pool->UpdatePool(tracksClone);
    
  }

  
  
   
  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskEmcalJetHMEC::Terminate(Option_t *) 
{
  //just terminate

}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHMEC::AcceptJet(AliEmcalJet *jet) 
{
  //applies all jet cuts except pt
  float jetphi = jet->Phi();
  if (jetphi>TMath::Pi())
    jetphi = jetphi-2*TMath::Pi();
  
  if ((jet->Phi()<fPhimin)||(jet->Phi()>fPhimax))
    return 0;
  if ((jet->Eta()<fEtamin)||(jet->Eta()>fEtamax))
    return 0;
  if (jet->Area()<fAreacut)
    return 0;
  //prevents 0 area jets from sneaking by when area cut == 0
  if (jet->Area()==0)
    return 0;
  
  //exclude jets with extremely high pt tracks which are likely misreconstructed
  if(jet->MaxTrackPt()>100)
    return 0;

  //passed all above cuts
  return 1;
  
}

//________________________________________________________________________

THnSparse* AliAnalysisTaskEmcalJetHMEC::NewTHnSparseF(const char* name, UInt_t entries)
{
   // generate new THnSparseF, axes are defined in GetDimParams()

   Int_t count = 0;
   UInt_t tmp = entries;
   while(tmp!=0){
      count++;
      tmp = tmp &~ -tmp;  // clear lowest bit
   }

   TString hnTitle(name);
   const Int_t dim = count;
   Int_t nbins[dim];
   Double_t xmin[dim];
   Double_t xmax[dim];

   Int_t i=0;
   Int_t c=0;
   while(c<dim && i<32){
      if(entries&(1<<i)){
      
         TString label("");
         GetDimParams(i, label, nbins[c], xmin[c], xmax[c]);
         hnTitle += Form(";%s",label.Data());
         c++;
      }
      
      i++;
   }
   hnTitle += ";";

   return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
}

void AliAnalysisTaskEmcalJetHMEC::GetDimParams(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   // stores label and binning of axis for THnSparse

   const Double_t pi = TMath::Pi();
   
   switch(iEntry){
      
   case 0:
      label = "V0 centrality (%)";
     
         nbins = 10;
         xmin = 0.;
         xmax = 100.;
         break;
      
      
   case 1:
      label = "corrected jet pt";
         nbins = 20;
         xmin = 0.;
         xmax = 200.;
          break;
      
      
   case 2:
      label = "track pT";
     
         nbins = 100;
         xmin = 0.;
         xmax = 10;
         break;
      
      
    case 3:
      label = "deltaR";
      nbins = 10;
      xmin = 0.;
      xmax = 5.0;
      break;



   case 4:
      label = "deltaEta";
      nbins = 24;
      xmin = -1.2;
      xmax = 1.2;
      break;


  case 5:
      label = "deltaPhi";
      nbins = 72;
      xmin = -0.5*pi;
      xmax = 1.5*pi;
      break;   
   
      
        
    case 6:
      label = "leading track";
      nbins = 13;
      xmin = 0;
      xmax = 50;
      break;
           
     case 7:
    
      label = "trigger track";
      nbins =10;
      xmin = 0;
      xmax = 50;
      break;

    case 8:
      label = "leading jet";
      nbins = 3;
      xmin = -0.5;
      xmax = 2.5;
      break;
  

   }

}


//_________________________________________________
// From CF event mixing code PhiCorrelations
TObjArray* AliAnalysisTaskEmcalJetHMEC::CloneAndReduceTrackList(TObjArray* tracks)
{
  // clones a track list by using AliPicoTrack which uses much less memory (used for event mixing)
  
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  
  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
  {
    AliVParticle* particle = (AliVParticle*) tracks->At(i);
    if(TMath::Abs(particle->Eta())>fTrkEta) continue;
    if(particle->Pt()<0.15)continue;

    Double_t trackpt=particle->Pt();

    Int_t hadbin=-1;
    if(trackpt<0.5) hadbin=0;
    else if(trackpt<1) hadbin=1;
    else if(trackpt<2) hadbin=2;
    else if(trackpt<3) hadbin=3;
    else if(trackpt<5) hadbin=4;
    else if(trackpt<8) hadbin=5;
    else if(trackpt<20) hadbin=6;

    if(hadbin>-1) fHistTrackEtaPhi[hadbin]->Fill(particle->Eta(),particle->Phi());


    tracksClone->Add(new AliPicoTrack(particle->Pt(), particle->Eta(), particle->Phi(), particle->Charge(), 0, 0, 0, 0));
  }
  
  return tracksClone;
}




