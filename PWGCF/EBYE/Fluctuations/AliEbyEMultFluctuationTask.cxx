
#include "TChain.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TArrayD.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODHeader.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliAODpidUtil.h"
#include "AliVTrack.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliCentrality.h"
#include "AliEbyEMultFluctuationTask.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliCentrality.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAODpidUtil.h"
#include "TRandom.h"
#include "AliTPCPIDResponse.h"
#include "AliTRDPIDResponse.h"


ClassImp(AliEbyEMultFluctuationTask)

//________________________________________________________________________
AliEbyEMultFluctuationTask::AliEbyEMultFluctuationTask(const char *name) 
: AliAnalysisTaskSE(name),fAOD(0), fAODVertex(0),fHistNchPt(0),fHistNchEta(0),fHistNchEtaCent(0),fHistNchPhi(0),fHistDCAxy(0),fHistDCAz(0),fHistnclus(0),fHistchi2ndf(0),fHistchi2ndfvscs(0),fHistVz(0),fHistMultV0A(0),fHistMultV0C(0),fHistMultV0total(0),My_ntuple(0),fOutputList(0),fCentralityEstimator("V0M"),fCentralityBins20(kFALSE),fCentralityCounter(0),fEventCounter(0),histcounter(0)


{
  // Constructor
  for(Int_t ibin=0;ibin<91;ibin++)
		{
			fMult[ibin]=NULL;
		}
		for(Int_t jbin=0;jbin<46;jbin++)
		{
		fMultTwo[jbin]=NULL;
		}
		for(Int_t kbin=0;kbin<15;kbin++)
		{
		fMultFive[kbin]=NULL;
		}
	
	    
  DefineInput(0, TChain::Class());
  
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
Bool_t AliEbyEMultFluctuationTask :: SelectEvent(AliAODVertex* vertex)
{

  if(vertex){ 
  
    if(vertex->GetNContributors() < 0) return kFALSE;
    
    
	  Double_t lvx = vertex->GetX();
	  Double_t lvy = vertex->GetY();
	  Double_t lvz = vertex->GetZ();
	  
	  fEventCounter->Fill(3);
	  
	  
	  if(TMath::Abs(lvx) > 0.3)  return kFALSE;
	  if(TMath::Abs(lvy) > 0.3) return kFALSE;
	  if(TMath::Abs(lvz) > 10) return kFALSE;
	  if(vertex->GetType()==AliAODVertex::kPileupSPD) return kFALSE;  
	  
	  fEventCounter->Fill(5);
	  fHistVz->Fill(lvz);		
	  return kTRUE;
  }
  else{
    return kFALSE;
  }
}
//_____________________________________________________________________________
Int_t AliEbyEMultFluctuationTask :: SelectTrack(AliAODTrack* track)
{
 
	Double_t eta = track->Eta();

     if(TMath :: Abs(eta)>0.8) return 0;

 if(!track->TestFilterBit(128)) return 0;

 Float_t dxy, dz;
		  
  dxy = track->DCA();
 dz  = track->ZAtDCA();
		    if(TMath ::Abs(dxy) >2.4 || TMath ::Abs(dz)>3.0) return 0;
		  // cout<<dxy<<dz<<endl;


  Double_t nclus = track->GetTPCClusterInfo(2,1);
  if(nclus<80) return 0;

//chi2cut

 Double_t chi2ndf = track->Chi2perNDF();
  if(chi2ndf>4) return 0;

return 1;
 
  
}

//________________________________________________________________________________________________
void AliEbyEMultFluctuationTask::UserCreateOutputObjects()
{
  fOutputList = new TList();


	
  My_ntuple = new TNtuple("My_ntuple", "My_ntuple", "Mult:Mult1:Mult2:nCentrality:fV0total:nBin2:nBin5:spdmult0:spdmult1:tracklets:run");
    fOutputList->Add(My_ntuple);

fOutputList->SetOwner(kTRUE);
fEventCounter = new TH1D("fEventCounter","EventCounter", 10, 0.5,10.5);
fOutputList->Add(fEventCounter);


histcounter = new TH1D("histcounter","histcounter", 10, 0.5,10.5);
fOutputList->Add(histcounter);


TString BinLebel;
fCentralityCounter = new TH1D("fEventCounter","EventCounter", 100, 0.5,100.5);
TAxis *xaxis = fCentralityCounter->GetXaxis();
 for(Int_t ibin=0;ibin<100;ibin++)
{
BinLebel="cent";BinLebel+=ibin;
xaxis->SetBinLabel(ibin,BinLebel.Data());

}
fOutputList->Add(fCentralityCounter);

fHistMultV0A = new TH1F("fHistMultV0A","Mult-VOA",22000,0,22000);
fHistMultV0C = new TH1F("fHistMultV0C","Mult-VOC",22000,0,22000);
fHistMultV0total = new TH1F("fHistMultV0total","Mult-VOTotal",22000,0,22000);


 fHistNchEta = new TH1D("fHistNchEta", "#eta distribution of charged tracks(for hybrid)", 200, -0.8, 0.8);
fHistNchEta->GetXaxis()->SetTitle("#eta");
fHistNchEta->GetYaxis()->SetTitle("dNch/d#eta");


fHistNchEtaCent = new TH1D("fHistNchEtaCent", "Central #eta distribution of charged tracks(for hybrid)", 200, -0.8, 0.8);
fHistNchEtaCent->GetXaxis()->SetTitle("#eta");
fHistNchEtaCent->GetYaxis()->SetTitle("dNch/d#eta(Central)");



fHistNchPhi = new TH1D("fHistNchPhi", "#phi distribution of charged tracks(for hybrid)", 200, -10, 10);
fHistNchPhi->GetXaxis()->SetTitle("#phi");
fHistNchPhi->GetYaxis()->SetTitle("dNch/d#phi");

fHistNchPt = new TH1D("fHistNchPt", "P_{T} distribution(charged tracks for hybrid)", 150, 0.1, 3.1);
fHistNchPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
fHistNchPt->GetYaxis()->SetTitle("dNch/dP_{T}(c/GeV)");


 fHistchi2ndfvscs = new TH2D("fHistchi2ndfvscs","Chi2/ndf from Gaussian fits of Multiplicity Distribution",95,0,95,60,0,6);
 fHistchi2ndfvscs->SetXTitle("%cs");
 fHistchi2ndfvscs->SetYTitle("Chi2/ndf");
 
 
fHistnclus=new TH1D("no of clusters","no of tracks after TPCcluster cut",200,0,200);
fHistchi2ndf= new TH1D("chi_square/ndf","no of tracks after chi_square/ndf cut",10,0,10);

 fHistVz = new TH1D("fHistVz","Primary vertex distribution - z coordinate",90,-15,15);
 fHistVz->SetXTitle("Vz");
 fHistVz->SetYTitle("Entries");


  fOutputList->Add(fHistMultV0A);
  fOutputList->Add(fHistMultV0C);
  fOutputList->Add(fHistMultV0total);
  fOutputList->Add(fHistnclus);
  fOutputList->Add(fHistchi2ndf); 
  fOutputList->Add(fHistNchPt);
  fOutputList->Add(fHistNchEta);
  fOutputList->Add(fHistNchEtaCent);

  fOutputList->Add(fHistNchPhi);
  fOutputList->Add(fHistVz);
  fOutputList->Add(fHistchi2ndfvscs);

 

TString hname;
TString htitle;
  
  for(Int_t i = 0; i < 25; i++)  {
    hname  = "fMult"; hname+=i;
    htitle = "Multiplicity in Cent Bin "; htitle+=i;
    fMult[i] = new TH1F(hname.Data(),htitle.Data(),5000,0.5,5000.5);
    fOutputList->Add(fMult[i]);
  }
  
  for(Int_t i = 25; i < 50; i++) {
   hname  = "fMult"; hname+=i;
    htitle = "Multiplicity in Cent Bin "; htitle+=i;
    fMult[i] = new TH1F(hname.Data(),htitle.Data(),3000,0.5,3000.5);
    fOutputList->Add(fMult[i]);
  }
  
  for(Int_t i = 50; i < 91; i++) {
    hname  = "fMult"; hname+=i;
    htitle = "Multiplicity in Cent Bin "; htitle+=i;
    fMult[i] = new TH1F(hname.Data(),htitle.Data(),2000,0.5,2000.5);
    fOutputList->Add(fMult[i]);
  }

//****************for 2% centrality bin**************************//
for(Int_t i = 0; i < 12; i++)  {
    hname  = "fMultTwo"; hname+=i;
    htitle = "Multiplicity in Cent Bin "; htitle+=i;
    fMultTwo[i] = new TH1F(hname.Data(),htitle.Data(),5000,0.5,5000.5);
    fOutputList->Add(fMultTwo[i]);
  }
  
  for(Int_t i = 12; i < 25; i++) {
   hname  = "fMultTwo"; hname+=i;
    htitle = "Multiplicity in Cent Bin "; htitle+=i;
    fMultTwo[i] = new TH1F(hname.Data(),htitle.Data(),3000,0.5,3000.5);
    fOutputList->Add(fMultTwo[i]);
  }
  
  for(Int_t i = 25; i <= 45; i++) {
    hname  = "fMultTwo"; hname+=i;
    htitle = "Multiplicity in Cent Bin "; htitle+=i;
    fMultTwo[i] = new TH1F(hname.Data(),htitle.Data(),2000,0.5,2000.5);
    fOutputList->Add(fMultTwo[i]);
  }

  
  
//*****************************for 5 % centrality bin********************//
for(Int_t i = 0; i < 5; i++)  {
    hname  = "fMultFive"; hname+=i;
    htitle = "Multiplicity in Cent Bin "; htitle+=i;
    fMultFive[i] = new TH1F(hname.Data(),htitle.Data(),5000,0.5,5000.5);
    fOutputList->Add(fMultFive[i]);
  }
  
  for(Int_t i = 5; i < 10; i++) {
   hname  = "fMultFive"; hname+=i;
    htitle = "Multiplicity in Cent Bin "; htitle+=i;
    fMultFive[i] = new TH1F(hname.Data(),htitle.Data(),3000,0.5,3000.5);
    fOutputList->Add(fMultFive[i]);
  }
  
  for(Int_t i = 10; i <= 14; i++) {
    hname  = "fMultFive"; hname+=i;
    htitle = "Multiplicity in Cent Bin "; htitle+=i;
    fMultFive[i] = new TH1F(hname.Data(),htitle.Data(),2000,0.5,2000.5);
    fOutputList->Add(fMultFive[i]);
  }

	
	}

//________________________________________________________________________
void AliEbyEMultFluctuationTask::UserExec(Option_t *) 
{  
  


// Post output data.
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    printf("ERROR: fAOD not available\n");
    return;

  }

Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
if(!isSelected) return;
 fEventCounter->Fill(1);
 //Counter++;


// Get the event vertex.
	fAODVertex = fAOD->GetPrimaryVertex();

	if (!fAODVertex) {	
	  return;
	}
	
	if(TMath::Abs((fAODVertex->GetZ())-((fAOD->GetPrimaryVertexSPD())->GetZ())) > 0.5) return;


    if (!SelectEvent(fAODVertex)) return;
	
    //***************VZERO-AMPLITUDE*************************//	

    	AliAODVZERO *aodV0 = fAOD->GetVZEROData();
	Float_t fV0Amult = aodV0->GetMTotV0A();
	Float_t fV0Cmult = aodV0->GetMTotV0C();      
       	Float_t fV0total = fV0Amult + fV0Cmult;

	//*****************SPD-CLUSTERS*************************//


	
		AliAODHeader *fHeader = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
		if(!fHeader) AliFatal("Not a standard AOD");
		Int_t spdmult0 = fHeader->GetNumberOfITSClusters(0);
		Int_t spdmult1 = fHeader->GetNumberOfITSClusters(1);
		Int_t run = fHeader->GetRunNumber();
		
		AliAODTracklets *fTracklets = fAOD->GetTracklets(); 
		Int_t tracklets = fTracklets->GetNumberOfTracklets();
       

		//********************CENTRALITY******************************************//	

    AliCentrality *centrality= fAOD->GetCentrality();
    Double_t cent=centrality->GetCentralityPercentile("V0M");
    //    cout<<"cent= "<<cent<<endl;	

 if(cent < 0 || cent >= 91) return;

Int_t nCentrality = (Int_t)cent;
	
Int_t nBin2= (Int_t)( nCentrality/2.0);
	
Int_t nBin5 = (Int_t)( nCentrality/5.0);
 if(nBin5==0){

histcounter->Fill(1);
 }
	

	
fCentralityCounter->Fill(nCentrality);

fEventCounter->Fill(7);

 
Double_t Mult=0.0;
	Double_t Mult1=0.0;
	Double_t Mult2=0.0;

 //printf("There are %d tracks in this event\n", fAOD->GetNumberOfTracks());

 
  for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
  AliAODTrack*  track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
  if(!track) AliFatal("Not a standard AOD");
    if (!track) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }




// Find the track classification.
Int_t tracktype = SelectTrack(track);


 if (tracktype==0) {
   continue;          
		}

 if(track->GetType()== AliAODTrack::kPrimary)
{
  Short_t Charge = track->Charge();
	Double_t eta = track->Eta();
 if(TMath :: Abs(Charge)==1 && nBin5==0){
		       fHistNchEtaCent->Fill(eta);
		       
		     }
		
	Double_t pt = track->Pt();
	if( pt < 0.3 || pt > 2.0 ) continue;
  Double_t nclus = track->GetTPCClusterInfo(2,1);
  Double_t chi2ndf = track->Chi2perNDF();
  fHistnclus->Fill(nclus);
  fHistchi2ndf->Fill(chi2ndf);
  fHistchi2ndfvscs->Fill(cent,chi2ndf);

		//Count the charge particles-----
		  if(TMath :: Abs(Charge)==1) Mult +=1;
		  if(Charge==1) Mult1 +=1;
		  	if(Charge==-1) Mult2 +=1;


	Double_t phi = track->Phi();
		     if(TMath :: Abs(Charge)==1) fHistNchPt->Fill(pt);
		    
		     
		  if(TMath :: Abs(Charge)==1) fHistNchEta->Fill(eta);
		   if(TMath :: Abs(Charge)==1) fHistNchPhi->Fill(phi);
		  

 }


  } //track loop 


  //****************Filling Tree*************************//
  
  My_ntuple->Fill(Mult,Mult1,Mult2,nCentrality,fV0total,nBin2,nBin5,spdmult0,spdmult1,tracklets,run);


fMult[nCentrality]->Fill(Mult);

if(nBin2<=45) {
fMultTwo[nBin2]->Fill(Mult);
 }
if(nBin5<=14){
fMultFive[nBin5]->Fill(Mult);
 }

 fHistMultV0A->Fill(fV0Amult);
 fHistMultV0C->Fill(fV0Cmult);
 fHistMultV0total->Fill(fV0total);

  PostData(1, fOutputList);
}
//________________________________________________________________________
void AliEbyEMultFluctuationTask::Terminate(Option_t *) 
{
  

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
 
}
