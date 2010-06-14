#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TRandom.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"

#include "AliLog.h"
#include "AliESDVertex.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliMultiplicity.h"

#include "AliAnalysisTaskHardSoft.h"
#include "AliExternalTrackParam.h"
#include "AliTrackReference.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenDPMjetEventHeader.h"

// Analysis Task for Hard and Soft event characteristics

// Author: E.Sicking

ClassImp(AliAnalysisTaskHardSoft)

//________________________________________________________________________
  AliAnalysisTaskHardSoft::AliAnalysisTaskHardSoft(const char *name) 
    : AliAnalysisTaskSE(name) 
    ,fUseMc(false)
    ,fUseNeutral(false)
    ,fRadiusCut(0.7)
    ,fTriggerPtCut(0.7)
    ,fAssociatePtCut(0.4)
    ,fMap(0x0)
    ,fMapLeading(0x0)
    ,fCuts(0)
    ,fFieldOn(kTRUE)
    ,fHists(0)
    
{
  for(Int_t i = 0;i< 2;i++){
    fPt[i]                  = 0;
    fEta[i]                 = 0;
    fPhi[i]                 = 0;
    fEtaPhi[i]              = 0;
    fEtaPhiLeading[i]       = 0;
    fNch[i]                 = 0;
    fPtLeading[i]           = 0;
    fPtLeadingNch[i]        = 0;
    fPtSumNch[i]            = 0;
    fPtAvNch[i]             = 0;
    fDPhiLeading[i]         = 0;
    fRadiusLeading[i]       = 0;
    fDPhiLeadingR[i]        = 0;
    fRadiusLeadingR[i]      = 0;
    fDPhiLeadingRS[i]       = 0;
    fRadiusLeadingRS[i]     = 0;
    fNchAssInR[i]           = 0;
    fTrigger[i]             = 0;

    for(Int_t j=0;j<100;j++){
      fDPhiLeadingNchBin[i][j]   = 0;
    }

    for(Int_t j=0;j<2;j++){
      fNchHardSoft[i][j]   = 0;
      fPtHardSoft[i][j]   = 0;
      fPtAvHardSoftNch[i][j]   = 0;
    }
  }
  DefineOutput(1,  TList::Class()); 
}


//________________________________________________________________________
void AliAnalysisTaskHardSoft::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fHists = new TList();

 

  TString labels[2];
  labels[0]="MC";
  labels[1]="ESD";

  for(Int_t i=0;i<2;i++){
    fPt[i]                 = new TH1F(Form("fPt%s",labels[i].Data()),
				      Form("fPt%s",labels[i].Data()) ,  
				      500, 0., 50);
    fEta[i]                = new TH1F (Form("fEta%s",labels[i].Data()),
				       Form("fEta%s",labels[i].Data()) ,  
				       100, -1., 1);
    fPhi[i]                = new TH1F(Form("fPhi%s",labels[i].Data()),
				      Form("fPhi%s",labels[i].Data()) ,  
				      360, 0.,2*TMath::Pi());
    fEtaPhi[i]             = new TH2F(Form("fEtaPhi%s",labels[i].Data()),
				      Form("fEtaPhi%s",labels[i].Data()) ,  
				      100,-1.,1.,
				      360, 0.,2*TMath::Pi());
    fEtaPhiLeading[i]      = new TH2F(Form("fEtaPhiLeading%s",labels[i].Data()),
				      Form("fEtaPhiLeading%s",labels[i].Data()) ,  
				      100,-1.,1.,
				      360, 0.,2*TMath::Pi());
    fNch[i]                = new TH1F(Form("fNch%s",labels[i].Data()),
				      Form("fNch%s",labels[i].Data()) ,  
				      250, -0.5, 249.5);
    fPtLeading[i]          = new TH1F(Form("fPtLeading%s",labels[i].Data()),
				      Form("fPtLeading%s",labels[i].Data()) ,  
				      500, 0., 50);
    fPtLeadingNch[i]       = new TProfile(Form("fPtLeadingNch%s",labels[i].Data()),
					  Form("fPtLeadingNch%s",labels[i].Data()) ,  
					  250, -0.5, 249.5);
    fPtSumNch[i]           = new TProfile(Form("fPtSumNch%s",labels[i].Data()),
					  Form("fPtSumNch%s",labels[i].Data()) ,  
					  250, -0.5, 249.5);
    fPtAvNch[i]            = new TProfile(Form("fPtAvNch%s",labels[i].Data()),
					  Form("fPtAvNch%s",labels[i].Data()) ,  
					  250, -0.5, 249.5);
    fDPhiLeading[i]        = new TH1F(Form("fDPhiLeading%s",labels[i].Data()),
				      Form("fDPhiLeading%s",labels[i].Data()) ,  
				      180, 0., TMath::Pi());
    fRadiusLeading[i]      = new TH1F(Form("fRadiusLeading%s",labels[i].Data()),
				      Form("fRadiusLeading%s",labels[i].Data()) ,  
				      180, 0., 2*TMath::Pi());
    fDPhiLeadingR[i]       = new TH1F(Form("fDPhiLeadingR%s",labels[i].Data()),
				      Form("fDPhiLeadingR%s",labels[i].Data()) ,  
				      180, 0., TMath::Pi());
    fRadiusLeadingR[i]     = new TH1F(Form("fRadiusLeadingR%s",labels[i].Data()),
				      Form("fRadiusLeadingR%s",labels[i].Data()) ,  
				      180, 0., 2*TMath::Pi());
    fDPhiLeadingRS[i]      = new TH1F(Form("fDPhiLeadingRS%s",labels[i].Data()),
				      Form("fDPhiLeadingRS%s",labels[i].Data()) ,  
				      180, 0., TMath::Pi());
    fRadiusLeadingRS[i]    = new TH1F(Form("fRadiusLeadingRS%s",labels[i].Data()),
				      Form("fRadiusLeadingRS%s",labels[i].Data()) ,  
				      180, 0., 2*TMath::Pi());
    fNchAssInR[i]          = new TProfile(Form("fNchAssInR%s",labels[i].Data()),
					  Form("fNchAssInR%s",labels[i].Data()) ,  
					  250, -0.5, 249.5);
    fTrigger[i]            = new TH1F(Form("fTrigger%s",labels[i].Data()),
				      Form("fTrigger%s",labels[i].Data()) ,  
				      250, -0.5, 249.5);
    for(Int_t j=0;j<100;j++){
      fDPhiLeadingNchBin[i][j]     = new TH1F(Form("fDPhiLeadingNchBin%s%02d",labels[i].Data(),j),
					      Form("fDPhiLeadingNchBin%s%02d",labels[i].Data(),j) ,  
					      180, 0., TMath::Pi());
    }
    
    for(Int_t j=0;j<2;j++){
      fNchHardSoft[i][j]     = new TH1F(Form("fNchHardSoft%s%d",labels[i].Data(),j),
					Form("fNchHardSoft%s%d",labels[i].Data(),j) ,  
					250, -0.5, 249.5);
      fPtHardSoft[i][j]     = new TH1F(Form("fPtHardSoft%s%d",labels[i].Data(),j),
				       Form("fPtHardSoft%s%d",labels[i].Data(),j) ,  
				       500, 0., 50.0);
      fPtAvHardSoftNch[i][j]     = new TProfile(Form("fPtAvHardSoftNch%s%d",labels[i].Data(),j),
						Form("fPtAvHardSoftNch%s%d",labels[i].Data(),j) ,  
						250, -0.5, 249.5);
    }
    
  }
  
 
  fHists->SetOwner();
  for(Int_t i=0;i<2;i++){
    fHists->Add(fPt[i]);
    fHists->Add(fEta[i]);
    fHists->Add(fPhi[i]);
    fHists->Add(fEtaPhi[i]);
    fHists->Add(fEtaPhiLeading[i]);
    fHists->Add(fNch[i]);
    fHists->Add(fPtLeading[i]);
    fHists->Add(fPtLeadingNch[i]);
    fHists->Add(fPtSumNch[i]);
    fHists->Add(fPtAvNch[i]);
    fHists->Add(fDPhiLeading[i]);
    fHists->Add(fRadiusLeading[i]);
    fHists->Add(fDPhiLeadingR[i]);
    fHists->Add(fRadiusLeadingR[i]);
    fHists->Add(fDPhiLeadingRS[i]);
    fHists->Add(fRadiusLeadingRS[i]);
    fHists->Add(fNchAssInR[i]);
    fHists->Add(fTrigger[i]);

    for(Int_t j=0;j<100;j++){
      fHists->Add(fDPhiLeadingNchBin[i][j]);
    }

    for(Int_t j=0;j<2;j++){
      fHists->Add(fNchHardSoft[i][j]);
      fHists->Add(fPtHardSoft[i][j]);
      fHists->Add(fPtAvHardSoftNch[i][j]);
    }
  }

 
  TH1::AddDirectory(oldStatus);
}

//__________________________________________________________

void AliAnalysisTaskHardSoft::UserExec(Option_t *) 
{
  Int_t nentries[2]         = {0,0};  // tracks in event (before selection)
  Int_t nTracksAll[2]       = {0,0};  // accepted tracks in event
  Int_t nTracksAssociate[2] = {0,0};  // accepted tracks in event within R=fRadiusCut around leading track
  
  Float_t pt[2]    = {0.,0.}; // pt
  Float_t eta[2]   = {0.,0.}; // eta
  Float_t phi[2]   = {0.,0.}; // phi
  Float_t ptSum[2] = {0.,0.}; // pt sum
  Float_t ptAv[2]  = {0.,0.}; // average pt


  Float_t ptLeading[2]     = {0.,0.}; // pt leading
  Float_t etaLeading[2]    = {0.,0.}; // eta leading
  Float_t phiLeading[2]    = {0.,0.}; // phi leading

  Float_t ptOthers[2]     = {0.,0.}; // pt others // for second track loop around leading track
  Float_t etaOthers[2]    = {0.,0.}; // eta others
  Float_t phiOthers[2]    = {0.,0.}; // phi others

  Double_t etaLeadingRandom[2]   = {0.,0.}; // eta leading for random particle position (from fMapLeading)
  Double_t phiLeadingRandom[2]   = {0.,0.}; // phi leading  "
  Double_t etaOthersRandom[2]    = {0.,0.}; // eta others  for random particle position (from fMap)
  Double_t phiOthersRandom[2]    = {0.,0.}; // phi others   "

  Int_t fEventType[2]={1,1}; //hard=0, soft=1 -> set to hard if trigger and associate particle 
                             //are within R=0.7 

  //get event and vertex cut :(MC and ESD)
  //---------------------
  //MC
  //---------------------
  AliStack *stack = 0x0; // needed for MC studies
  Float_t vzMC=0.;       // MC vertex position in z
  if(fUseMc==true){
    stack = MCEvent()->Stack();
    AliGenEventHeader*  header = MCEvent()->GenEventHeader();
    TArrayF mcV;
    header->PrimaryVertex(mcV);
    vzMC = mcV[2];
  }
  //---------------------
  //ESD
  //---------------------
  AliVEvent* event = InputEvent();
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    return;
  }
  if(Entry()==0){
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
    if(esd)Printf("We are reading from ESD");
  }
  const AliVVertex* vertex = event->GetPrimaryVertex();
  Float_t vz = vertex->GetZ();
  //---------------------




  //Selection of particle(mcesd=0) or esdtracks(mcesd=1)
  for(Int_t mcesd=0;mcesd<2;mcesd++){
    if(mcesd==0){
      if(fUseMc==false) //MC part can be switched off for real data by function SetUserMc(kFALSE)
	continue;
    }

    // vertex cut and number of particles/tracks per event
    //---------------------------------------
    if(mcesd==0){//mc particles
      if (TMath::Abs(vzMC) > 10.) return;
      nentries[mcesd]=MCEvent()->GetNumberOfTracks();
    }
    else{// esd tracks
      if (TMath::Abs(vz) > 10.) return;
      nentries[mcesd]=event->GetNumberOfTracks();
    }//---------------------------------------
   



    // arrays for leading track determination (done with TMath::Sort of Array)
    Float_t * ptArray = new Float_t[nentries[mcesd]];

    //array of track properties
    Float_t * etaArray = new Float_t[nentries[mcesd]];
    Float_t * phiArray = new Float_t[nentries[mcesd]];

    Int_t *pindex  = new Int_t[nentries[mcesd]];
    for (Int_t i = 0; i < nentries[mcesd]; i++) {
      ptArray[i]=0.;
      etaArray[i]=0.;
      phiArray[i]=0.;
      pindex[i]=0;
    }
  
      

   
    //first track loop
    for (Int_t iTrack = 0; iTrack < nentries[mcesd]; iTrack++) {

      pt[mcesd]  = 0.;
      eta[mcesd] = 0.;
      phi[mcesd] = 0;
         
      //get properties of mc particle
      if(mcesd==0){//mc
	AliMCParticle* mcP = (AliMCParticle*) MCEvent()->GetTrack(iTrack);
	// Primaries only
	if (!(stack->IsPhysicalPrimary(mcP->Label()))) continue;
	if(!fUseNeutral)if (mcP->Particle()->GetPDG()->Charge() == 0) continue;
	//same cuts as on ESDtracks
	if(TMath::Abs(mcP->Eta())>0.9)continue;
	if(mcP->Pt()<0.2)continue;
	if(mcP->Pt()>200)continue;

	pt[mcesd]  = mcP->Pt();
	eta[mcesd] = mcP->Eta();
	phi[mcesd] = mcP->Phi();
      } 

      //get properties of esdtracks
      else{//esd
	AliVParticle *track = event->GetTrack(iTrack);
	if (!track) {
	  Printf("ERROR: Could not receive track %d", iTrack);
	  continue;
	}
	AliESDtrack *esdtrack =  dynamic_cast<AliESDtrack*>(track);
	if(!esdtrack)continue;
	if (!fCuts->AcceptTrack(esdtrack)) continue;
 
	pt[mcesd]=esdtrack->Pt();
	eta[mcesd]=esdtrack->Eta();
	phi[mcesd]=esdtrack->Phi();
      }
    
      ptArray[nTracksAll[mcesd]]    = pt[mcesd];  // fill pt array
      etaArray[nTracksAll[mcesd]]   = eta[mcesd]; // fill eta array
      phiArray[nTracksAll[mcesd]++] = phi[mcesd]; // count tracks and fill phi array

      fPt[mcesd]     -> Fill(pt[mcesd]);
      fEta[mcesd]    -> Fill(eta[mcesd]);
      fPhi[mcesd]    -> Fill(phi[mcesd]);
      fEtaPhi[mcesd] -> Fill(eta[mcesd],phi[mcesd]);
      

    }//end first track loop

    fNch[mcesd]    -> Fill(nTracksAll[mcesd]);
    

    //find leading pt tracks
    if(nentries[mcesd]>0) TMath::Sort(nentries[mcesd], ptArray, pindex, kTRUE);  
    //for(Int_t i=0;i<nTracksAll[mcesd];i++){//show just the filled entries, skip empty ones.
    //     printf("%d:  pt = %f, number %i \n",mcesd, ptArray[pindex[i]],i);
    //}
    
    
    if(nTracksAll[mcesd]>0){
      fPtLeading[mcesd]->Fill(ptArray[pindex[0]]);  //first entry in array is highest
      fPtLeadingNch[mcesd]->Fill(nTracksAll[mcesd],ptArray[pindex[0]]);

      for(Int_t i=0;i<nTracksAll[mcesd];i++){      //calculate ptsum
	ptSum[mcesd]+=ptArray[pindex[i]];
      }
      ptAv[mcesd]=ptSum[mcesd]/nTracksAll[mcesd];  //calculate <pt>

      fPtSumNch[mcesd]->Fill(nTracksAll[mcesd],ptSum[mcesd]);
      fPtAvNch[mcesd]->Fill(nTracksAll[mcesd],ptAv[mcesd]);
    }
    

    if(nTracksAll[mcesd]>1){ // require at least two tracks (leading and prob. accosicates)
      
      //Leading track properties
      ptLeading[mcesd]  = ptArray[pindex[0]];
      etaLeading[mcesd] = etaArray[pindex[0]];
      phiLeading[mcesd] = phiArray[pindex[0]];

      fEtaPhiLeading[mcesd] -> Fill(etaLeading[mcesd],phiLeading[mcesd]);

      fMapLeading->GetRandom2(etaLeadingRandom[mcesd],phiLeadingRandom[mcesd]);
           
      //second track loop for event propoerties around leading tracks with pt>triggerPtCut
      //loop only over already accepted tracks except leading track 
      if(ptLeading[mcesd]>fTriggerPtCut){
	
	fTrigger[mcesd]->Fill(nTracksAll[mcesd]); // how often is there a trigger particle at a certain Nch bin
	
	for (Int_t iTrack = 1; iTrack < nTracksAll[mcesd]; iTrack++) { // starting at second highest track

	  ptOthers[mcesd]   = ptArray[pindex[iTrack]];
	  etaOthers[mcesd]  = etaArray[pindex[iTrack]];
	  phiOthers[mcesd]  = phiArray[pindex[iTrack]];

	  fMap->GetRandom2(etaOthersRandom[mcesd],phiOthersRandom[mcesd]);
	  
	  if(ptOthers[mcesd]>fAssociatePtCut){ // only tracks which fullfill associate pt cut

	    //1. real data

	    Float_t dPhi=TMath::Abs(phiOthers[mcesd]-phiLeading[mcesd]);
	    if(dPhi>TMath::Pi())      dPhi=2*TMath::Pi()-dPhi;
	    Float_t dEta=etaOthers[mcesd]-etaLeading[mcesd];
	    
	    Float_t radius=TMath::Sqrt(dPhi*dPhi+dEta*dEta);
	    fRadiusLeading[mcesd]->Fill(radius);
	    fDPhiLeading[mcesd]->Fill(dPhi);
	    if(nTracksAll[mcesd]<100)fDPhiLeadingNchBin[mcesd][nTracksAll[mcesd]]->Fill(dPhi);
	    
	    if(radius<fRadiusCut){
	      fEventType[mcesd]=0;
	      nTracksAssociate[mcesd]++;
	    }
	    
	    //2. random position
	    Float_t dPhiR=TMath::Abs(phiOthersRandom[mcesd]-phiLeadingRandom[mcesd]);
	    if(dPhiR>TMath::Pi())      dPhiR=dPhiR-2*TMath::Pi();
	    Float_t dEtaR=etaOthersRandom[mcesd]-etaLeadingRandom[mcesd];
	    
	    Float_t radiusR=TMath::Sqrt(dPhiR*dPhiR+dEtaR*dEtaR);
	    fRadiusLeadingR[mcesd]->Fill(radiusR);
	    fDPhiLeadingR[mcesd]->Fill(dPhiR);


	    
	    //3. random position of leading particle
	    Float_t dPhiRS=TMath::Abs(phiOthers[mcesd]-phiLeadingRandom[mcesd]);
	    if(dPhiRS>TMath::Pi())      dPhiRS=dPhiRS-2*TMath::Pi();
	    Float_t dEtaRS=etaOthers[mcesd]-etaLeadingRandom[mcesd];
	    
	    Float_t radiusRS=TMath::Sqrt(dPhiRS*dPhiRS+dEtaRS*dEtaRS);
	    fRadiusLeadingRS[mcesd]->Fill(radiusRS);
	    fDPhiLeadingRS[mcesd]->Fill(dPhiRS);
	  }
	  

	  
	}
	//fill histogram with number of tracks (pt>fAssociatePtCut) around leading track
	fNchAssInR[mcesd]->Fill(nTracksAll[mcesd],nTracksAssociate[mcesd]);
	
      }
    }
    
    fNchHardSoft[mcesd][fEventType[mcesd]]->Fill(nTracksAll[mcesd]);

    for(Int_t i=0;i<nTracksAll[mcesd];i++){
      fPtHardSoft[mcesd][fEventType[mcesd]]->Fill(ptArray[i]);
    }
    
    fPtAvHardSoftNch[mcesd][fEventType[mcesd]]->Fill(nTracksAll[mcesd],ptAv[mcesd]);

   
  }//double loop over mcP and ESD


 
  // Post output data.
  PostData(1, fHists);
}      





//________________________________________________________________________
void AliAnalysisTaskHardSoft::Terminate(Option_t *) 
{


}  





