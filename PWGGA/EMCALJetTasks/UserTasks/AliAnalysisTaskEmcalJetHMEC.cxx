// $Id: $

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
  fESD(0), 
  fPoolMgr(0x0), 
  fOutputList(0),
  fHistTrackPt(0),
  fHistCentrality(0), 
  fHistJetEtaPhi(0), 
  fHistTrackEtaPhi(0), 
  fHistJetHEtaPhi(0) 
{
  // Default Constructor

  for(Int_t icent = 0; icent<6; ++icent){
    fHistJetPt[icent]=0;
    fHistJetPtBias[icent]=0;
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
  fDoEventMixing(0),
  fMixingTracks(50000),
  fESD(0), 
  fPoolMgr(0x0), 
  fOutputList(0), 
  fHistTrackPt(0),
  fHistCentrality(0), 
  fHistJetEtaPhi(0), 
  fHistTrackEtaPhi(0), 
  fHistJetHEtaPhi(0),
  fNevents(0),
  fTindex(0),
  fTrigBufferIndex(0),
  fCountAgain(0), 
  fhnMixedEvents(0x0)
{
  // Constructor
  for(Int_t icent = 0; icent<6; ++icent){
    fHistJetPt[icent]=0;
    fHistJetPtBias[icent]=0;
    fHistJetPtTT[icent]=0;
    for(Int_t iptjet = 0; iptjet<5; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){	
	fHistJetH[icent][iptjet][ieta]=0;
	fHistJetHBias[icent][iptjet][ieta]=0;
	fHistJetHTT[icent][iptjet][ieta]=0;
      }
    }
  }

    for(Int_t i=0; i<10; i++) {
       for(Int_t j=0; j<6; j++) {
	    fTrigBuffer[i][j]=0;
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

  fHistJetEtaPhi = new TH2F("fHistJetEtaPhi","Jet eta-phi",900,-1.8,1.8,640,-3.2,3.2);
  fHistTrackEtaPhi = new TH2F("fHistTrackEtaPhi","Track eta-phi",900,-1.8,1.8,640,-3.2,3.2);
  fHistJetHEtaPhi = new TH2F("fHistJetHEtaPhi","Jet-Hadron deta-dphi",900,-1.8,1.8,640,-1.6,4.8);

  char name[200];

 
  for(Int_t icent = 0; icent<6; ++icent){
    sprintf(name,"fHistJetPt_%i",icent);   
    fHistJetPt[icent] = new TH1F(name,name,200,0,200);
    fOutputList->Add(fHistJetPt[icent]);

    sprintf(name,"fHistJetPtBias_%i",icent);   
    fHistJetPtBias[icent] = new TH1F(name,name,200,0,200);
    fOutputList->Add(fHistJetPtBias[icent]);

    sprintf(name,"fHistJetPtTT_%i",icent);   
    fHistJetPtTT[icent] = new TH1F(name,name,200,0,200);
    fOutputList->Add(fHistJetPtTT[icent]);

    for(Int_t iptjet = 0; iptjet<5; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){	
	sprintf(name,"fHistJetH_%i_%i_%i",icent,iptjet,ieta);   
	fHistJetH[icent][iptjet][ieta]=new TH2F(name,name,64,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
	fOutputList->Add(fHistJetH[icent][iptjet][ieta]);

	sprintf(name,"fHistJetHBias_%i_%i_%i",icent,iptjet,ieta);   
	fHistJetHBias[icent][iptjet][ieta]=new TH2F(name,name,64,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
	fOutputList->Add(fHistJetHBias[icent][iptjet][ieta]);

	sprintf(name,"fHistJetHTT_%i_%i_%i",icent,iptjet,ieta);   
	fHistJetHTT[icent][iptjet][ieta]=new TH2F(name,name,64,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
	fOutputList->Add(fHistJetHTT[icent][iptjet][ieta]);

      }
    }
  }

  if(fDoEventMixing){    
     UInt_t cifras = 0; // bit coded, see GetDimParams() below 
     cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<7; 
     fhnMixedEvents = NewTHnSparseF("fhnMixedEvents", cifras);
     }



  fOutputList->Add(fHistTrackPt);
  fOutputList->Add(fHistCentrality);
  fOutputList->Add(fHistJetEtaPhi);
  fOutputList->Add(fHistTrackEtaPhi);
  fOutputList->Add(fHistJetHEtaPhi);
   fOutputList->Add(fhnMixedEvents);


  PostData(1, fOutputList);


  //Event Mixing
  Int_t trackDepth = fMixingTracks; 
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
 
  Int_t nZvtxBins  = 7+1+7;
  // bins for second buffer are shifted by 100 cm
  Double_t vertexBins[] = { -7, -5, -3, -1, 1, 3, 5, 7, 93, 95, 97, 99, 101, 103, 105, 107 };
  Double_t* zvtxbin = vertexBins;

  Int_t nCentralityBins  = 100;
  Double_t centralityBins[nCentralityBins];
  for(Int_t ic=0; ic<nCentralityBins; ic++){
    centralityBins[ic]=1.0*ic;
  }
  //Double_t* centbin = centralityBins;

  //cout << "filling centrality bins" <<endl;
  //Int_t nCentralityBins  = fHistCentrality->GetNbinsX();
  //Double_t* centralityBins = (Double_t*)fHistCentrality->GetXaxis()->GetXbins()->GetArray();

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
  if (pt>=10 && pt<15)
    ptbin = 0;
  else if (pt>=15 && pt<20)
    ptbin = 1;
  else if (pt>=20 && pt<25)
    ptbin = 2;
  else if (pt>=25 && pt<30)
    ptbin = 3;
  else if (pt>=30)
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


  fHistCentrality->Fill(fcent);
  Int_t centbin = GetCentBin(fcent);

  TClonesArray *jets = 0;
  TClonesArray *tracks = 0;

  tracks = dynamic_cast<TClonesArray*>(list->FindObject(fTracksName));
  if (!tracks) {
    AliError(Form("Pointer to tracks %s == 0", fTracksName.Data() ));
    return;
  }
  const Int_t Ntracks=tracks->GetEntries();

  jets= dynamic_cast<TClonesArray*>(list->FindObject(fJetsName));
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
      
      //pt,eta,phi,centrality
      float jetphi = jet->Phi();
      if (jetphi>TMath::Pi())
	jetphi = jetphi-2*TMath::Pi();
      
      if ((jet->Phi()<fPhimin)||(jet->Phi()>fPhimax))
	continue;
      if ((jet->Eta()<fEtamin)||(jet->Eta()>fEtamax))
	continue;
      //fHistAreavsRawPt[centbin]->Fill(jet->Pt(),jet->Area());
      if (jet->Area()<fAreacut)
	continue;
      //prevents 0 area jets from sneaking by when area cut == 0
      if (jet->Area()==0)
	continue;

      Double_t jetPt = jet->Pt();

      if(highestjetpt<jetPt){
	ijethi=ijet;
	highestjetpt=jetPt;
      }

    }


    
  //Only look at the highest pT jet in the event

  if(ijethi>-1){
    AliEmcalJet *jet = static_cast<AliEmcalJet*>(jets->At(ijethi)); 
    
    Double_t jetphi = jet->Phi();
    Double_t jetPt = jet->Pt();
    Double_t jeteta=jet->Eta();

      fHistJetPt[centbin]->Fill(jet->Pt());

      if ((jet->MaxTrackPt()>6) || (jet->MaxClusterPt()>6))
	fHistJetPtBias[centbin]->Fill(jet->Pt());


      fHistJetEtaPhi->Fill(jet->Eta(),jetphi);

      //fHistDeltaPtvsArea->Fill(jetPt,jet->Area());

      if(iTT>0){
	AliVTrack* TT = static_cast<AliVTrack*>(tracks->At(iTT));
	if(TMath::Abs(jetphi-TT->Phi()-TMath::Pi())<0.6) passedTTcut=1;
	else passedTTcut=0;
      }

      if(passedTTcut)
	fHistJetPtTT[centbin]->Fill(jet->Pt());
   

  if (highestjetpt>10) {
   
    for (Int_t iTracks = 0; iTracks < Ntracks; iTracks++) 
      {
	AliVTrack* track = static_cast<AliVTrack*>(tracks->At(iTracks));
	if (!track) {
	  printf("ERROR: Could not receive track %d\n", iTracks);
	  continue;
	}
	
	  if(TMath::Abs(track->Eta())>0.9) continue;

	  fHistTrackPt->Fill(track->Pt());
	  	  
	  if (track->Pt()<0.15)
	    continue;
	  
	  Double_t trackphi = track->Phi();
	  if (trackphi > TMath::Pi())
	    trackphi = trackphi-2*TMath::Pi();

	  Double_t tracketa=track->Eta();
      
	  Double_t deta=tracketa-jeteta;
	  Int_t ieta=GetEtaBin(deta);

	  //Jet pt, track pt, dPhi,deta,fcent
	  Double_t dphijh = RelativePhi(jetphi,trackphi);
	  if (dphijh < -0.5*TMath::Pi())
	    dphijh= dphijh+ 2*TMath::Pi();

	  Int_t iptjet=-1;
	  iptjet=GetpTjetBin(jetPt);

	  fHistJetH[centbin][iptjet][ieta]->Fill(dphijh,track->Pt());
	  fHistJetHEtaPhi->Fill(deta,dphijh);
	  fHistTrackEtaPhi->Fill(tracketa,trackphi);
	  if ((jet->MaxTrackPt()>6) || (jet->MaxClusterPt()>6))
	    fHistJetHBias[centbin][iptjet][ieta]->Fill(dphijh,track->Pt());
	  
	  if(passedTTcut)
	  fHistJetHTT[centbin][iptjet][ieta]->Fill(dphijh,track->Pt());


	} //track loop
  }//jet pt cut

  
  // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
  TObjArray* tracksClone = CloneAndReduceTrackList(tracks);
  //delete tracks;

  Double_t fvertex[3]={0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(fvertex);
  Double_t zVtx=fvertex[3];
  
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

    pool->PrintInfo();

    
    if (!pool)
      AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fcent, zVtx));


    if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 || pool->GetCurrentNEvents() >= 5) 
    {
      
      Int_t nMix = pool->GetCurrentNEvents();
      
    
      // Fill mixed-event histos here  
      for (Int_t jMix=0; jMix<nMix; jMix++) 
      {
	TObjArray* bgTracks = pool->GetEvent(jMix);
	const Int_t Nbgtrks = bgTracks->GetEntries();
	for(Int_t ibg=0; ibg<Nbgtrks; ibg++){
	  AliPicoTrack *part = static_cast<AliPicoTrack*>(bgTracks->At(ibg));         
	  if(!part) continue;

	  Double_t DPhi = jetphi - part->Phi();
	  Double_t DEta = jeteta - part->Eta();
	  Double_t DR=TMath::Sqrt(DPhi*DPhi+DEta*DEta);
	  if(DPhi<-0.5*TMath::Pi()) DPhi+=2.*TMath::Pi();
	  if(DPhi>3./2.*TMath::Pi()) DPhi-=2.*TMath::Pi();
	  Double_t triggerEntries[7] = {fcent,jetPt,part->Pt(),DR,DEta,DPhi,0.0};                      
	  fhnMixedEvents->Fill(triggerEntries,1./nMix);
	}
 
      }
    }
    pool->UpdatePool(tracksClone);
  
  }



  }

  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskEmcalJetHMEC::Terminate(Option_t *) 
{
  //just terminate

}

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
      nbins = 15;
      xmin = 0.;
      xmax = 1.5;
      break;



   case 4:
      label = "deltaEta";
      nbins = 8;
      xmin = -1.6;
      xmax = 1.6;
      break;


  case 5:
      label = "deltaPhi";
      nbins = 64;
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


  

   }

}


//_________________________________________________
// From CF event mixing code PhiCorrelations
TObjArray* AliAnalysisTaskEmcalJetHMEC::CloneAndReduceTrackList(TObjArray* tracks)
{
  // clones a track list by using AliDPhiBasicParticle which uses much less memory (used for event mixing)
  
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  
  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
  {
    AliVParticle* particle = (AliVParticle*) tracks->At(i);
    tracksClone->Add(new AliPicoTrack(particle->Pt(), particle->Eta(), particle->Phi(), particle->Charge(), 0, 0, 0, 0));
  }
  
  return tracksClone;
}



