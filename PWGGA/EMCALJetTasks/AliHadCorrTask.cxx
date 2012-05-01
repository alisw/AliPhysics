#include <TClonesArray.h>
#include <TParticle.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include <TLorentzVector.h>

#include "AliCentrality.h"
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliFJWrapper.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliPicoTrack.h"
#include "AliVEventHandler.h"

#include "AliHadCorrTask.h"

ClassImp(AliHadCorrTask)

//________________________________________________________________________
AliHadCorrTask::AliHadCorrTask(const char *name) : 
  AliAnalysisTaskSE("AliHadCorrTask"),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fHadCorr(0),
  fJets(0),
  fMinPt(0.15),
  fOutCaloName("CaloClustersOut"),
  fOutClusters(0),
  fOutputList(0),
  fHistNclusvsCent(0),
  fHistNclusMatchvsCent(0),
  fHistEbefore(0),
  fHistEafter(0),
  fHistEoPCent(0),
  fHistNMatchCent(0)
{

 for(Int_t i=0; i<4; i++){
    for(Int_t j=0; j<5; j++){
      fHistMatchEtaPhi[i][j]=0x0;
    }
    fHistMatchEvsP[i]=0x0;
    fHistMatchdRvsEP[i]=0x0;
  }

  // Standard constructor.
  if (!name)
    return;

  SetName(name);
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

AliHadCorrTask::AliHadCorrTask() : 
  AliAnalysisTaskSE("AliHadCorrTask"),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fHadCorr(0),
  fJets(0),
  fMinPt(0.15),
  fOutCaloName("CaloClustersOut"),
  fOutClusters(0),
  fOutputList(0),
  fHistNclusvsCent(0),
  fHistNclusMatchvsCent(0),
  fHistEbefore(0),
  fHistEafter(0),
  fHistEoPCent(0),
  fHistNMatchCent(0)
{
  // Standard constructor.

  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";
}


//________________________________________________________________________
AliHadCorrTask::~AliHadCorrTask()
{
  // Destructor
}

//________________________________________________________________________
void AliHadCorrTask::UserCreateOutputObjects()
{

  AliVEventHandler* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  
  if (!handler) {
    AliError("Input handler not available!");
    return;
  }
 
  if (handler->InheritsFrom("AliESDInputHandler")) {
    fOutClusters = new TClonesArray("AliESDCaloCluster");
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    fOutClusters = new TClonesArray("AliAODCaloCluster");
  }
  else {
    AliError("Input handler not recognized!");
    return;
  }

  fOutClusters->SetName(fOutCaloName);
 
  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner();

  char name[200];

   for(int icent=0; icent<4; icent++){
    for(int ipt=0; ipt<5; ipt++){
      sprintf(name,"fHistMatchEtaPhi_%i_%i",icent,ipt);
      fHistMatchEtaPhi[icent][ipt] = new TH2F(name,name,400,-0.2,0.2,1600,-0.8,0.8);
      fOutputList->Add(fHistMatchEtaPhi[icent][ipt]);
    }

    sprintf(name,"fHistMatchEvsP_%i",icent);
    fHistMatchEvsP[icent]=new TH2F(name,name,400,0.,200.,1000,0.,10.);
    sprintf(name,"fHistMatchdRvsEP_%i",icent);
    fHistMatchdRvsEP[icent]=new TH2F(name,name,1000,0.,1.,1000,0.,10.);
    
    fOutputList->Add(fHistMatchEvsP[icent]);
    fOutputList->Add(fHistMatchdRvsEP[icent]);
   }


  fHistNclusvsCent=new TH1F("Nclusvscent","NclusVsCent",100,0,100);
  fHistNclusMatchvsCent=new TH1F("NclusMatchvscent","NclusMatchVsCent",100,0,100);
  fHistEbefore=new TH1F("Ebefore","Ebefore",100,0,100);
  fHistEafter=new TH1F("Eafter","Eafter",100,0,100);
  fHistEoPCent = new TH2F("EoPCent","EoPCent",100,0.0,100.,1000,0.0,10.);
  fHistNMatchCent = new TH2F("NMatchesCent","NMatchesCent",100,0.0,100.,101,-0.5,100.5);

  fOutputList->Add(fHistNclusMatchvsCent);
  fOutputList->Add(fHistNclusvsCent);
  fOutputList->Add(fHistEbefore);
  fOutputList->Add(fHistEafter);
  fOutputList->Add(fHistEoPCent);
  fOutputList->Add(fHistNMatchCent);

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliHadCorrTask::UserExec(Option_t *) 
{

  if (!(InputEvent()->FindListObject(fOutCaloName)))
    InputEvent()->AddObject(fOutClusters);
  else fOutClusters->Delete();


  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  TClonesArray *tracks = 0;
  TClonesArray *clus   = 0;
  TList *l = InputEvent()->GetList();

  AliCentrality *centrality = 0;
  centrality = dynamic_cast<AliCentrality*>(l->FindObject("Centrality"));
  Float_t fCent=-1; 
  if(centrality)
    fCent = centrality->GetCentralityPercentile("V0M");
  else
    fCent=99;//probably pp data

  if(fCent<0) return;

  if (fTracksName == "Tracks")
    am->LoadBranch("Tracks");
  tracks = dynamic_cast<TClonesArray*>(l->FindObject(fTracksName));
  if (!tracks) {
    AliError(Form("Pointer to tracks %s == 0", fTracksName.Data() ));
    return;
  }
  const Int_t Ntrks = tracks->GetEntries();
  
  if (fCaloName == "CaloClusters")
    am->LoadBranch("CaloClusters");
  clus = dynamic_cast<TClonesArray*>(l->FindObject(fCaloName));
  if (!clus) {
    AliError(Form("Pointer to clus %s == 0", fCaloName.Data() ));
    return;
  }

  Int_t centbin=-1;
  if(fCent>=0 && fCent<10) centbin=0;
  else if(fCent>=10 && fCent<30) centbin=1;
  else if(fCent>=30 && fCent<50) centbin=2;
  else if(fCent>=50 && fCent<=100) centbin=3;
 
  if (clus) {
    Double_t vertex[3] = {0, 0, 0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
    const Int_t Nclus = clus->GetEntries();
    for (Int_t iClus = 0, clusCount=0; iClus < Nclus; ++iClus) {
      AliVCluster *c = dynamic_cast<AliVCluster*>(clus->At(iClus));
      if (!c)
        continue;
      if (!c->IsEMCAL())
        continue;
      TLorentzVector nPart;

      c->SetEmcCpvDistance(-1);
      c->SetTrackDistance(999,999);
      c->GetMomentum(nPart, vertex);
      Double_t dEtaMin  = 1e9;
      Double_t dPhiMin  = 1e9;
      Int_t    imin     = -1;
      Double_t totaltrkpt =0.0;
      Int_t Nmatches = 0;
      Double_t energy = nPart.P();
      if(energy<fMinPt) continue;
      for(Int_t t = 0; t<Ntrks; ++t) {
        AliVTrack *track = dynamic_cast<AliVTrack*>(tracks->At(t));
	if (! track)
	  continue;
        Double_t etadiff=999;
        Double_t phidiff=999;
        AliPicoTrack::GetEtaPhiDiff(track,c,phidiff,etadiff);
        Double_t dR = TMath::Sqrt(etadiff*etadiff+phidiff*phidiff);
        Double_t dRmin = TMath::Sqrt(dEtaMin*dEtaMin+dPhiMin*dPhiMin);
        if(dR > 25) 
          continue;
	if(dR<dRmin){
          dEtaMin = etadiff;
          dPhiMin = phidiff;
          imin = t;
	}
	  if(fHadCorr>1) {
	    Double_t fpt=track->P();
	    Int_t ptbin=-1;
	    if(fpt>=0 && fpt<0.5) ptbin=0;
	    else if(fpt>=0.5 && fpt<2.) ptbin=1;
	    else if(fpt>=2. && fpt<3.) ptbin=2;
	    else if(fpt>=3. && fpt<5.) ptbin=3;
	    else if(fpt>=5.) ptbin=4;
	    if(fpt>0) fHistMatchEtaPhi[centbin][ptbin]->Fill(etadiff,phidiff);
	    if(track->P()>0) fHistMatchdRvsEP[centbin]->Fill(dR,energy/track->P());
	  }

	  if (TMath::Abs(phidiff)<0.05 && TMath::Abs(etadiff)<0.025) { // pp cuts!!!
	    Nmatches++;
	    totaltrkpt=totaltrkpt+track->P();
	  }
        
      }
      fHistNMatchCent->Fill(fCent,Nmatches);
      c->SetEmcCpvDistance(imin);
      c->SetTrackDistance(dPhiMin, dEtaMin);

      fHistNclusvsCent->Fill(fCent);
      if(Nmatches>0) fHistNclusMatchvsCent->Fill(fCent);
      fHistEbefore->Fill(fCent,energy);
      if (fHadCorr>0) {
	//to subtract only the closest track set fHadCor to a %
	//to subtract all tracks within the cut set fHadCor to %+1
	if(fHadCorr>1){
	  if (totaltrkpt>0){
	    double EoP=energy/totaltrkpt;
	    fHistEoPCent->Fill(fCent,EoP);
	    fHistMatchEvsP[centbin]->Fill(energy,EoP);
	  }
	  energy -= (fHadCorr-1)*totaltrkpt;
	  if (energy<0)
	    continue;
	}
	else{
	  if (imin>=0) {
	    dPhiMin = c->GetTrackDx();
	    dEtaMin = c->GetTrackDz();
	    Double_t dR=TMath::Sqrt(dEtaMin*dEtaMin+dPhiMin*dPhiMin);
	    
	    AliVTrack *t = dynamic_cast<AliVTrack*>(tracks->At(imin));
	    if (t) {
	      if (t->Pt()<fMinPt)
		continue;
		
	      Double_t fpt=t->P();
	      Int_t ptbin=-1;
	      if(fpt>=0 && fpt<1.) ptbin=0;
	      else if(fpt>=1. && fpt<2.) ptbin=1;
	      else if(fpt>=2. && fpt<3.) ptbin=2;
	      else if(fpt>=3. && fpt<5.) ptbin=3;
	      else if(fpt>=5.) ptbin=4;
	      
	      fHistMatchEtaPhi[centbin][ptbin]->Fill(dEtaMin,dPhiMin);
	      

	      if (t->P()>0){
		fHistMatchEvsP[centbin]->Fill(energy,energy/t->P());
		fHistEoPCent->Fill(fCent,energy/t->P());
		fHistMatchdRvsEP[centbin]->Fill(dR,energy/t->P());
	      }
	      if (TMath::Abs(dPhiMin)<0.05 && TMath::Abs(dEtaMin)<0.025) { // pp cuts!!!
		energy -= fHadCorr*t->P();
	      }
	      if (energy<0)
		continue;
	    }
	  }
	}
	fHistEafter->Fill(fCent,energy);

      }//end had correction if

      if (energy>0){//Output new corrected clusters
	AliVCluster *oc;
	if (c->InheritsFrom("AliESDCaloCluster")) {
	  oc = new ((*fOutClusters)[clusCount]) AliESDCaloCluster(*(dynamic_cast<AliESDCaloCluster*>(c)));
	}
	else if (c->InheritsFrom("AliAODCaloCluster")) {
	  oc = new ((*fOutClusters)[clusCount]) AliAODCaloCluster(*(dynamic_cast<AliAODCaloCluster*>(c)));
	}
	else {
	  AliError("Cluster type not recognized (nor ESD neither AOD)!");
	  continue;
	}
	oc->SetE(energy);
	clusCount++;
      }
    }
  }
}


//________________________________________________________________________
void AliHadCorrTask::Terminate(Option_t *) 
{
  
}

//________________________________________________________________________
