#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObjArray.h"
#include "TString.h"
#include "TParticle.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"

#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODpidUtil.h"
#include "AliAnalysisUtils.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliInputEventHandler.h"

#include "AliAnalysisTaskParticleEfficiency.h"


ClassImp(AliAnalysisTaskParticleEfficiency)
//ClassImp(AliAnalysisTaskParticleEfficiency)

using std::cout;
using std::endl;

//_______________________________________________________

AliAnalysisTaskParticleEfficiency::AliAnalysisTaskParticleEfficiency(const Char_t *partName) :
AliAnalysisTaskSE(partName), centrality(0), fHistoList(0),  fHistEv(0), fpidResponse(0)
{

  for(Int_t i = 0; i < MULTBINS*PARTTYPES; i++)  {
    fGeneratedMCPrimaries[i] = NULL;
    fMCPrimariesThatAreReconstructed[i] = NULL;
    fReconstructedAfterCuts[i] = NULL;
    fReconstructedNotPrimaries[i] = NULL;
    fReconstructedPrimaries[i] = NULL;
    fContamination[i] = NULL;
  }
  for ( Int_t i = 0; i < 11; i++) { 
    fHistQA[i] = NULL;
    if(i<3) fHistQA2D[i] = NULL;
  }
  
  /* init track cuts */
  //fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
  /*fTrackCuts =  AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    if( !fTrackCuts ) return;
    fTrackCuts->SetMinNClustersTPC(70);*/
  
  
  /* create output */
  fHistoList = new TList();
  fHistoList->SetOwner(kTRUE);

  DefineInput(0, TChain::Class());
  // DefineOutput(0, TTree::Class()); 
  DefineOutput(1, TList::Class());

 
}

//_______________________________________________________

AliAnalysisTaskParticleEfficiency::~AliAnalysisTaskParticleEfficiency()
{
  /* if(centrality) delete centrality;
s128     if(fHistoList) delete fHistoList;
     if(vertex) delete vertex;
     if(vtxSPD) delete vtxSPD;*/
}

//_______________________________________________________

void AliAnalysisTaskParticleEfficiency::UserCreateOutputObjects()
{
   
  TString hname1, hname2, hname3, hname4, hname5;
  
  TString htitle1, htitle2, htitle3, htitle4,htitle5;
  
  TString hname1M, hname2M, hname3M, hname4M, hname5M, hname;
  
  TString htitle1M, htitle2M, htitle3M, htitle4M, htitle5M, htitle;

  TString parttypename = "None";

  for(Int_t j = 0; j < PARTTYPES; j++)  {
    if (j==0) parttypename="All";
    else if (j==1) parttypename="Pion";
    else if (j==2) parttypename="Kaon";
    else if (j==3) parttypename="Proton";

    for(Int_t i = 0; i < MULTBINS; i++)  {
      hname1  = "hGeneratedMCPrimariesEffM"; hname1+=i; hname1+=parttypename;
      htitle1 = "Kinematic level eta_pT (prim only) M"; htitle1+=i; htitle1+=parttypename;
      fGeneratedMCPrimaries[i*PARTTYPES+j] = new TH2F(hname1.Data(),htitle1.Data(),50, -1.5, 1.5,1000,0.,10.0);
    
      hname3  = "hMCPrimariesThatAreReconstructedM"; hname3+=i; hname3+=parttypename;
      htitle3 = "Reconstructed level eta_pT (prim only) M"; htitle3+=i; htitle3+=parttypename;
      fMCPrimariesThatAreReconstructed[i*PARTTYPES+j] = new TH2F(hname3.Data(),htitle3.Data(),50, -1.5, 1.5,1000,0.,10.0);
    
      hname2  = "hHistoReconstructedAfterCutsM"; hname2+=i; hname2+=parttypename;
      htitle2 = "Total Reconstructed tracks M "; htitle2+=i; htitle2+=parttypename;
      fReconstructedAfterCuts[i*PARTTYPES+j] = new TH2F(hname2.Data(),htitle2.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname4  = "hHistoReconstructedNotPrimariesM"; hname4+=i; hname4+=parttypename;
      htitle4 = "Reconstructed level eta_pT (not primaries) M"; htitle4+=i; htitle4+=parttypename;
      fReconstructedNotPrimaries[i*PARTTYPES+j] = new TH2F(hname4.Data(),htitle4.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname4  = "hHistoReconstructedPrimariesM"; hname4+=i; hname4+=parttypename;
      htitle4 = "Reconstructed level eta_pT (primaries) M"; htitle4+=i; htitle4+=parttypename;
      fReconstructedPrimaries[i*PARTTYPES+j] = new TH2F(hname4.Data(),htitle4.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname5  = "hContaminationM"; hname5+=i; hname5+=parttypename;
      htitle5 = "Contamination M"; htitle5+=i; htitle5+=parttypename;
      fContamination[i*PARTTYPES+j] = new TH2F(hname5.Data(),htitle5.Data(),6000, -3000, 3000.,1000,0.,10.0);
 

      fReconstructedAfterCuts[i*PARTTYPES+j]->Sumw2();
      fReconstructedNotPrimaries[i*PARTTYPES+j]->Sumw2();
      fReconstructedPrimaries[i*PARTTYPES+j]->Sumw2();
    }
    hname  = "pidTPCdEdx";  hname+=parttypename;
    htitle = parttypename + " TPC dEdx vs. momentum";
    fHistQAPID[0][j] = new TH2F(hname, htitle, 100, 0.0, 5.0, 250, 0.0, 500.0);
    hname  = "pidTOFTime";  hname+=parttypename;
    htitle = parttypename + " TOF Time vs. momentum";
    fHistQAPID[1][j] = new TH2F(hname, htitle, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    hname  = "pidTOFNSigma";  hname+=parttypename;
    htitle = parttypename + " TOF NSigma vs. momentum";
    fHistQAPID[2][j] = new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hname  = "pidTPCNSigma";  hname+=parttypename;
    htitle = parttypename + " TPC NSigma vs. momentum";
    fHistQAPID[3][j] = new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hname  = "pidTPCTOFNSigma";  hname+=parttypename;
    htitle = parttypename + " TPC vs TOF NSigma";
    fHistQAPID[4][j] = new TH2F(hname,htitle, 200, -10.0, 10.0, 200, -10.0, 10.0);

  }

  fHistEv = new TH1F("fHistEv", "Multiplicity", 100, 0, 100);
  fHistoList->Add(fHistEv);

  for(Int_t i = 0; i < MULTBINS; i++)  {
    hname = "fHistEventCutsM";
    hname+= i;
    
    fHistEvCuts[i] = new TH1F(hname,Form("Event Cuts M%d",i) , 4, 0, 5);
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(1,"All");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(2,"NoVertex");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(3,"PileUp");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(4,"z-vertex>10");
    fHistoList->Add(fHistEvCuts[i]);


    hname  = "hMisidentificationM"; hname+=i; 
    htitle = "Misidentification Fraction M"; htitle+=i; 
    fMisidentification[i] = new TH2F(hname.Data(),htitle.Data(), 3, 0.5, 3.5, 4 , 0, 4);
    fMisidentification[i]->GetXaxis()->SetBinLabel(1,"Pions, MC");
    fMisidentification[i]->GetXaxis()->SetBinLabel(2,"Kaons, MC");
    fMisidentification[i]->GetXaxis()->SetBinLabel(3,"Protons, MC");
    fMisidentification[i]->GetYaxis()->SetBinLabel(1,"Pions, Data");
    fMisidentification[i]->GetYaxis()->SetBinLabel(2,"Kaons, Data");
    fMisidentification[i]->GetYaxis()->SetBinLabel(3,"Protons, Data");
    fMisidentification[i]->GetYaxis()->SetBinLabel(4,"Other, Data");
    fHistoList->Add(fMisidentification[i]);

  }

  fHistQA[0] = new TH1F("fHistVtx", "Z vertex distribution", 100, -15., 15.);
  fHistQA[1] = new TH1F("fHistnTpcCluster", "n TPC Cluster", 100, 0., 200.);
  fHistQA[2] = new TH1F("fHistnTpcClusterF", "n TPC Cluster findable", 100, 0., 200.);
  fHistQA[3] = new TH1F("dcaHistDcaXY1D", "DCA XY", 210, -2.1, 2.1);
  fHistQA[4] = new TH1F("dcaHistDcaZ1D", "DCA Z", 210, -2.1, 2.1);
  fHistQA[5] = new TH1F("fHistChi2Tpc", "Chi2 TPC", 100, 0., 8.);
  fHistQA[6] = new TH1F("fHistpT", "pT distribution",1000,0.,10.0);
  fHistQA[7] = new TH1F("fHistPhi", "Phi distribution" , 100, -TMath::Pi(), TMath::Pi());
  fHistQA[8] = new TH1F("fHistEta", "Eta distribution" , 100, -2, 2);
 
  fHistQA[9] = new TH1F("fHistEventCuts", "Event Cuts" , 4, 0, 5);
  fHistQA[9]->GetXaxis()->SetBinLabel(1,"All");
  fHistQA[9]->GetXaxis()->SetBinLabel(2,"NoVertex");
  fHistQA[9]->GetXaxis()->SetBinLabel(3,"PileUp");
  fHistQA[9]->GetXaxis()->SetBinLabel(4,"z-vertex>10");


  fHistQA[10] = new TH1F("fHistTrackCuts", "Track Cuts" , 7, 0.5, 7.5);
  fHistQA[10]->GetXaxis()->SetBinLabel(1,"AllTracksInEvents");
  fHistQA[10]->GetXaxis()->SetBinLabel(2,"GetTrack");
  fHistQA[10]->GetXaxis()->SetBinLabel(3,"Filter bit");
  fHistQA[10]->GetXaxis()->SetBinLabel(4,"Eta");
  fHistQA[10]->GetXaxis()->SetBinLabel(5,"Pt");
  fHistQA[10]->GetXaxis()->SetBinLabel(6,"DCA");
  fHistQA[10]->GetXaxis()->SetBinLabel(7,"Electron Rejection");

  fHistQA2D[0] = new TH2F("dcaHistDcaXY","DCA XY",50, 0, 5,210, -2.1, 2.1);
  fHistQA2D[1] = new TH2F("dcaHistDcaZ","DCA Z", 50, 0, 5, 210, -2.1, 2.1);
  fHistQA2D[2] = new TH2F("fPhiEta","Eta-Phi",100, -2, 2, 100, -TMath::Pi(), TMath::Pi());



  for ( Int_t i = 0; i < 11; i++)
    {
      fHistoList->Add(fHistQA[i]);
      if(i<3) fHistoList->Add(fHistQA2D[i]);
      if(i<5) {
	for(Int_t j = 0 ; j<PARTTYPES; j++)
	  fHistoList->Add(fHistQAPID[i][j]);
      }
    }
    
  for (Int_t i = 0; i < MULTBINS*PARTTYPES; i++){
    fHistoList->Add(fGeneratedMCPrimaries[i]);
    fHistoList->Add(fMCPrimariesThatAreReconstructed[i]);
    fHistoList->Add(fReconstructedAfterCuts[i]);
    fHistoList->Add(fReconstructedNotPrimaries[i]);
    fHistoList->Add(fReconstructedPrimaries[i]);
    fHistoList->Add(fContamination[i]);

  }
 
  //********** PID ****************

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fpidResponse = inputHandler->GetPIDResponse();
  cout<<"*******"<< fpidResponse<<endl;

  // ************************

  PostData(1, fHistoList);
}


//_____________________________________________________________________

bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{

    if (mom > 0.5) {
        if (TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < 3)
            return true;
	}
    else {
        if (TMath::Abs(nsigmaTPCPi) < 3)
            return true;
    }
/*
  if(mom<0.65)
    {
      if(nsigmaTOFPi<-999.)
	{
	  if(mom<0.35 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
	  else if(mom<0.5 && mom>=0.35 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
	  else if(mom>=0.5 && TMath::Abs(nsigmaTPCPi)<2.0) return true;
	  else return false;
	}
      else if(TMath::Abs(nsigmaTOFPi)<3.0 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
    }
  else if(nsigmaTOFPi<-999.)
    {
      return false;
    }
  else if(mom<1.5 && TMath::Abs(nsigmaTOFPi)<3.0 && TMath::Abs(nsigmaTPCPi)<5.0) return true;
  else if(mom>=1.5 && TMath::Abs(nsigmaTOFPi)<2.0 && TMath::Abs(nsigmaTPCPi)<5.0) return true;
*/
  return false;
}

bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{
    if (mom > 0.5) {
        if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < 3)
            return true;
	}
    else {
        if (TMath::Abs(nsigmaTPCK) < 3)
            return true;
    }

 /* if(mom<0.4)
    {
      if(nsigmaTOFK<-999.)
	{
	  if(TMath::Abs(nsigmaTPCK)<2.0) return true;
	}
      else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
    }
  else if(mom>=0.4 && mom<=0.6)
    {
      if(nsigmaTOFK<-999.)
	{
	  if(TMath::Abs(nsigmaTPCK)<2.0) return true;
	}
      else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
    }
  else if(nsigmaTOFK<-999.)
    {
      return false;
    }
  else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
*/
  return false;
}

bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP)
{

    if (mom > 0.5) {
        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < 3)
            return true;
	}
    else {
        if (TMath::Abs(nsigmaTPCP) < 3)
            return true;
    }

  // if (mom > 0.8 && mom < 2.5) {
  //   if ( TMath::Abs(nsigmaTPCP) < 3.0 && TMath::Abs(nsigmaTOFP) < 3.0)
  //     return true;
  // }
  // else if (mom > 2.5) {
  //   if ( TMath::Abs(nsigmaTPCP) < 3.0 && TMath::Abs(nsigmaTOFP) < 2.0)
  //     return true;
  // }
  // else {
  //   if (TMath::Abs(nsigmaTPCP) < 3.0)
  //     return true;
  // }
  

  return false;
}


bool IsElectron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP)
{
  if(TMath::Abs(nsigmaTPCE)<3 && TMath::Abs(nsigmaTPCPi)>3 && TMath::Abs(nsigmaTPCK)>3 && TMath::Abs(nsigmaTPCP)>3)
      return true;
   else
     return false;
}

//_______________________________________________________

void AliAnalysisTaskParticleEfficiency::UserExec(Option_t *)
{

  /***Get Event****/
  //AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!aodEvent) return;
  AliAODHeader *fAODheader = dynamic_cast<AliAODHeader*>(aodEvent->GetHeader());
  if(!fAODheader) AliFatal("Not a standard AOD");
  Double_t mult = fAODheader->GetRefMultiplicity();
// AliCentrality* alicent= aodEvent->GetCentrality(); //in PbPb and pPb
//  Double_t mult = alicent->GetCentralityPercentile("V0A"); //in pPb
//  Double_t mult = alicent->GetCentralityPercentile("V0A"); //in PbPb
  fHistEv->Fill(mult); 
 

  //******************
  // load MC array
  // arrayMC =  (TClonesArray*)aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  //  if(!arrayMC) {
  //  printf("AliAnalysisTaskParticleEficiency::UserExec: MC particles branch not found!\n");
  //return;
  //}
    
  // load MC header
  //mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  //if(!mcHeader) {
  //printf("AliAnalysisTaskSEDplusCorrelations::UserExec: MC header branch not found!\n");
  //return;
  // }
  //*********************   
 


  // EVENT SELECTION ********************

  fHistQA[9]->Fill(1);

  // collision candidate 
  // if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB)) return;


  //****** Multiplicity selection *********
  Int_t fcent = -999;  
  //if(mult >= 0 && mult <=20)  fcent = 0;
  //else if(mult >= 20 && mult <=39) fcent = 1;
  //else if(mult >= 40 && mult <=59) fcent = 2;
  //else if(mult >= 60 && mult <=90) fcent = 3;
  //else if(mult >= 99990 && mult <=99936) fcent = 4;
  //else if(mult >= 999937 && mult <=99944) fcent = 5;
  //else if(mult >= 999945 && mult <=99957) fcent = 6;
  //else if(mult >= 999958 && mult <=99149) fcent = 6;
  //else fcent = 7;
  //if (fcent == 7) return;

  if(mult >= 2&& mult <=150)  fcent = 0;
  else if(mult >= 2 && mult <=19) fcent = 1;
  else if(mult >= 20 && mult <=49) fcent = 2;
  else if(mult >= 50 && mult <=150) fcent = 3;
  else return;
 
  if(fcent==0)fHistEvCuts[0]->Fill(1);
  else if(fcent==1)fHistEvCuts[1]->Fill(1);
  else if(fcent==2)fHistEvCuts[2]->Fill(1);
  else if(fcent==3)fHistEvCuts[3]->Fill(1);

  //"ESDs/pass2/AOD049/*AliAOD.root");
  const AliAODVertex* vertex =(AliAODVertex*) aodEvent->GetPrimaryVertex();
  if (!vertex || vertex->GetNContributors()<=0) return;

  fHistQA[9]->Fill(2);
 if(fcent==0)fHistEvCuts[0]->Fill(2);
  else if(fcent==1)fHistEvCuts[1]->Fill(2);
  else if(fcent==2)fHistEvCuts[2]->Fill(2);
  else if(fcent==3)fHistEvCuts[3]->Fill(2);

  AliAnalysisUtils *anaUtil=new AliAnalysisUtils();
    
  Bool_t fpA2013 = kFALSE;
  Bool_t fMVPlp = kFALSE;
  Bool_t fisPileUp = kTRUE;
  Int_t fMinPlpContribMV = 0;
  Int_t fMinPlpContribSPD = 3;

  if(fpA2013)
    if(anaUtil->IsVertexSelected2013pA(aodEvent)==kFALSE) return;
 
  if(fMVPlp) anaUtil->SetUseMVPlpSelection(kTRUE);
  else anaUtil->SetUseMVPlpSelection(kFALSE);
 
  if(fMinPlpContribMV) anaUtil->SetMinPlpContribMV(fMinPlpContribMV);
  if(fMinPlpContribSPD) anaUtil->SetMinPlpContribSPD(fMinPlpContribSPD);

  if(fisPileUp)
    if(anaUtil->IsPileUpEvent(aodEvent)) return;

  delete anaUtil;   

 fHistQA[9]->Fill(3);
  if(fcent==0)fHistEvCuts[0]->Fill(3);
  else if(fcent==1)fHistEvCuts[1]->Fill(3);
  else if(fcent==2)fHistEvCuts[2]->Fill(3);
  else if(fcent==3)fHistEvCuts[3]->Fill(3);

  //TString vtxTtl = vertex->GetTitle();
  //if (!vtxTtl.Contains("VertexerTracks")) return;
  Float_t zvtx = vertex->GetZ();
  if (TMath::Abs(zvtx) > 10) return;
  fHistQA[0]->Fill(zvtx);
  fHistQA[9]->Fill(4);
  if(fcent==0)fHistEvCuts[0]->Fill(4);
  else if(fcent==1)fHistEvCuts[1]->Fill(4);
  else if(fcent==2)fHistEvCuts[2]->Fill(4);
  else if(fcent==3)fHistEvCuts[3]->Fill(4);

 //**** getting MC array ******
  TClonesArray  *arrayMC;

  arrayMC = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));


  //get stack 
  //AliStack *mcStack = mcEvent->Stack();
  //if (!mcStack) return;
  //***********************


  // old vertex selection 
  /*const AliESDVertex *vertex = esdEvent->GetPrimaryVertex();
    if (vertex->GetNContributors() < 1) return;

    //z-vertex cut
    if (TMath::Abs(vertex->GetZ()) > 10.) return;
  

    const AliESDVertex *vtxSPD = esdEvent->GetPrimaryVertexSPD();*/
  // Double_t zVertex = vtxSPD->GetZ();

  //cout << "Event  Z vtx ==========> " << vertex->GetZ() <<endl;
  // centrality selection 
  //AliCentrality *centrality = aodEvent->GetHeader()->GetCentralityP();
  //if (centrality->GetQuality() != 0) return;
  //Double_t cent = centrality->GetCentralityPercentileUnchecked("V0M");
  //if(cent < 0 || cent > 100.) return;


//copying pid information for FB 128
  int labels[20000];
  for (int il=0; il<20000; il++) labels[il] = -1;

  // looking for global tracks and saving their numbers to copy from them PID information to TPC-only tracks in the main loop over tracks
  for (int i=0;i<aodEvent->GetNumberOfTracks();i++) {
    const AliAODTrack *aodtrack=dynamic_cast<const AliAODTrack*>(aodEvent->GetTrack(i));
    if(!aodtrack) AliFatal("Not a standard AOD");
    if (!aodtrack->TestFilterBit(128)) {
      if(aodtrack->GetID() < 0) continue;
      labels[aodtrack->GetID()] = i;
    }
  }

  //RECONSTRUCTED TRACKS 

  TObjArray recoParticleArray[PARTTYPES];

  fHistQA[10]->Fill(1,aodEvent->GetNumberOfTracks());
  //loop over AOD tracks 
  for (Int_t iTracks = 0; iTracks < aodEvent->GetNumberOfTracks(); iTracks++) {
    //get track 
    
    //AliESDtrack* track = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(esdEvent),iTracks);
    AliAODTrack *track = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTracks));
    if(!track) AliFatal("Not a standard AOD"); 
    if (!track)continue;
    fHistQA[10]->Fill(2);

    //UInt_t filterBit = (1 << (0));
    UInt_t filterBit = 96;
    if(!track->TestFilterBit(filterBit))continue;	

    //if(!track->IsHybridGlobalConstrainedGlobal())continue;
    //if((track->IsHybridGlobalConstrainedGlobal())==false)continue;
    // if(!track->IsHybridTPCConstrainedGlobal())continue;	
    // if(!track->IsTPCConstrained())continue;	
    //if(!track->IsGlobalConstrained())continue;
    //if((track->TestFilterMask(AliAODTrack::kTrkTPCOnly)==false))continue;//cut0_BIT(0)
  
    //   if((track->IsHybridGlobalConstrainedGlobal())==false)
    //  continue;//def_BIT(272)

    //if((track->TestFilterMask(AliAODTrack::kTrkGlobal)==false))continue;//cut1_BIT(5)

    fHistQA[10]->Fill(3);
     
    if(track->Eta() < -0.8 || track->Eta() > 0.8)
      continue; 
    fHistQA[10]->Fill(4);

    if (track->Pt() < 0.2 || track->Pt() > 20)
      continue;
    fHistQA[10]->Fill(5);

    //single track cuts
    // if(track->Chi2perNDF() > 4.0) continue;
    // if(track->GetTPCNcls() < 70) continue;

    //DCA
    
    Double_t DCAXY;
    Double_t DCAZ;
    //  if(filterBit==(1 << (7))){
    DCAXY = -TMath::Abs(track->DCA());
    DCAZ = -TMath::Abs(track->ZAtDCA());
 
      if(!(DCAXY==-999 || DCAZ==-999)){
	//if(TMath::Abs(DCAXY) > 0.0182 + 0.035*TMath::Power(track->Pt(), -1.01)) continue; //XY, Pt dep
	//no DCA cut
	//if(TMath::Abs(DCAXY) > 1000.0) {continue;} //XY
	//if(TMath::Abs(DCAZ) > 1000.0) {continue;} //Z
      }
    else {
      // code from Michael and Prabhat from AliAnalysisTaskDptDptCorrelations
      // const AliAODVertex* vertex = (AliAODVertex*) aodEvent->GetPrimaryVertex(); (already defined above)
      float vertexX  = -999.;
      float vertexY  = -999.;
      float vertexZ  = -999.;

      if(vertex) {
	Double32_t fCov[6];
	vertex->GetCovarianceMatrix(fCov);
	if(vertex->GetNContributors() > 0) {
	  if(fCov[5] != 0) {
	    vertexX = vertex->GetX();
	    vertexY = vertex->GetY();
	    vertexZ = vertex->GetZ();

	  }
	}
      }

      Double_t pos[3];
      track->GetXYZ(pos);

      Double_t DCAX = pos[0] - vertexX;
      Double_t DCAY = pos[1] - vertexY;
      DCAZ = pos[2] - vertexZ;
      DCAXY = TMath::Sqrt((DCAX*DCAX) + (DCAY*DCAY));

      //if(TMath::Abs(DCAXY) > 0.0182 + 0.035*TMath::Power(track->Pt(), -1.01)) continue; //XY, Pt dep
      //if(TMath::Abs(impactD) > 0.44 + 0.07*TMath::Power(tPt, -1.94)) continue; //XY, Pt dep
      //no DCA cut
      //if(TMath::Abs(DCAXY) > 1000.0) continue;
      //if(TMath::Abs(DCAZ) > 1000.0) continue;
    }

    fHistQA[10]->Fill(6);

    AliAODTrack* aodtrackpid;

    //for FB 128 - tpc only tracks
    if(filterBit==(1 << (7))) {
      aodtrackpid = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(labels[-1-aodEvent->GetTrack(iTracks)->GetID()]));
      if(!aodtrackpid) AliFatal("Not a standard AOD");
    }
    else
      aodtrackpid = track;


 //PANOS--------------------------
    
    AliMCEvent* mcEvent = MCEvent();
    if (!mcEvent) {
      AliError("ERROR: Could not retrieve MC event");
      return;
    }

    Int_t nMCParticles = mcEvent->GetNumberOfTracks();
    Int_t labelp = TMath::Abs(track->GetLabel());
    if(labelp > nMCParticles) continue;

    AliAODMCParticle *AODmcTrack = (AliAODMCParticle*) mcEvent->GetTrack(labelp);

    /*
    //contamination from secondaries without electron rejection
    if (!AODmcTrack->IsPhysicalPrimary()) {
      fReconstructedNotPrimaries[PARTTYPES*fcent]->Fill(track->Eta(), track->Pt());
    }
    else{
      fReconstructedPrimaries[PARTTYPES*fcent]->Fill(track->Eta(), track->Pt());
    }
    */
    //-------------------------
   
    //Electron rejection
    double nSigmaTPCPi = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kPion);
    double nSigmaTPCK = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kKaon);
    double nSigmaTPCP = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kProton);
    double nSigmaTPCe = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kElectron);
    if(IsElectron(nSigmaTPCe,nSigmaTPCPi,nSigmaTPCK,nSigmaTPCP))
      continue;
    fHistQA[10]->Fill(7);
    
    fHistQA[1]->Fill(track->GetTPCClusterInfo(2,1)); 
    //fHistQA[2]->Fill(track->GetTPCNclsF());
     fHistQA[3]->Fill(DCAXY);
     fHistQA[4]->Fill(DCAZ);
    Float_t chi2Tpc = track->Chi2perNDF();
    fHistQA[5]->Fill(chi2Tpc);
    fHistQA[6]->Fill(track->Pt());

    float px=track->Px(); float py=track->Py();  float ph=atan2(py,px); //track->Phi()
    float tPt = track->Pt();

    fHistQA[7]->Fill(ph);
    fHistQA[8]->Fill(track->Eta());
    fHistQA2D[2]->Fill(track->Eta(),ph);

    fHistQA2D[0]->Fill(tPt,DCAXY);
    fHistQA2D[1]->Fill(tPt,DCAZ);


 //PANOS
    //contamination from secondaries with electron rejection
    // if (!AODmcTrack->IsPhysicalPrimary()) {
    //   fReconstructedNotPrimaries[PARTTYPES*fcent]->Fill(track->Eta(), track->Pt());
    // }
    // else{
    //   fReconstructedPrimaries[PARTTYPES*fcent]->Fill(track->Eta(), track->Pt());
    // }
    

    //PID monitors
    double nSigmaTOFPi = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kPion);
    double nSigmaTOFK = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kKaon);
    double nSigmaTOFP = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kProton);

    float tdEdx = aodtrackpid->GetTPCsignal();
    float tTofSig = aodtrackpid->GetTOFsignal();
    double pidTime[5]; aodtrackpid->GetIntegratedTimes(pidTime);


    fHistQAPID[0][0]->Fill(tPt,tdEdx);
    fHistQAPID[1][0]->Fill(tPt,tTofSig-pidTime[2]);//pion
    fHistQAPID[2][0]->Fill(tPt,nSigmaTOFPi);
    fHistQAPID[3][0]->Fill(tPt,nSigmaTPCPi);
    fHistQAPID[4][0]->Fill(nSigmaTPCPi,nSigmaTOFPi);

    if (IsPionNSigma(track->P(),nSigmaTPCPi, nSigmaTOFPi)){
      fHistQAPID[0][1]->Fill(tPt,tdEdx);
      fHistQAPID[1][1]->Fill(tPt,tTofSig-pidTime[2]);//pion
      fHistQAPID[2][1]->Fill(tPt,nSigmaTOFPi);
      fHistQAPID[3][1]->Fill(tPt,nSigmaTPCPi);
      fHistQAPID[4][1]->Fill(nSigmaTPCPi,nSigmaTOFPi);
    }
    if (IsKaonNSigma(track->P(),nSigmaTPCK, nSigmaTOFK)){
      fHistQAPID[0][2]->Fill(tPt,tdEdx);
      fHistQAPID[1][2]->Fill(tPt,tTofSig-pidTime[3]);//kaon
      fHistQAPID[2][2]->Fill(tPt,nSigmaTOFK);
      fHistQAPID[3][2]->Fill(tPt,nSigmaTPCK);
      fHistQAPID[4][2]->Fill(nSigmaTPCK,nSigmaTOFK);
    }
    if (IsProtonNSigma(track->P(),nSigmaTPCP, nSigmaTOFP)){
      fHistQAPID[0][3]->Fill(tPt,tdEdx);
      fHistQAPID[1][3]->Fill(tPt,tTofSig-pidTime[4]);//proton
      fHistQAPID[2][3]->Fill(tPt,nSigmaTOFP);
      fHistQAPID[3][3]->Fill(tPt,nSigmaTPCP);
      fHistQAPID[4][3]->Fill(nSigmaTPCP,nSigmaTOFP);
    }

    fReconstructedAfterCuts[fcent]->Fill(track->Eta(), track->Pt());//Fills hist. for all reconstructed particles after cuts

    if(!arrayMC){
      continue;
    }
    //get coresponding MC particle 
    Int_t label = TMath::Abs(track->GetLabel());
    AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(label);

   //getting no. of tracks for each particle species after all the cuts:

    //********* PID - pions ********
     if (IsPionNSigma(track->P(),nSigmaTPCPi, nSigmaTOFPi)){
       fReconstructedAfterCuts[PARTTYPES*fcent+1]->Fill(track->Eta(), track->Pt());
       if (!MCtrk) continue;
       recoParticleArray[1].Add(MCtrk);
       }
       //Fills for all identified pions found after cuts (reconstructed) - numerator for Efficiency

     //********* PID - kaons ********
     if (IsKaonNSigma(track->P(),nSigmaTPCK, nSigmaTOFK)){
       fReconstructedAfterCuts[PARTTYPES*fcent+2]->Fill(track->Eta(), track->Pt());
       if (!MCtrk) continue;
       recoParticleArray[2].Add(MCtrk);
       }
       //Fills for all identified kaons found after cuts (reconstructed) - numerator for Efficiency

    //********* PID - protons ********
     if (IsProtonNSigma(track->P(),nSigmaTPCP, nSigmaTOFP)){
       fReconstructedAfterCuts[PARTTYPES*fcent+3]->Fill(track->Eta(), track->Pt());
       if (!MCtrk) continue;
       recoParticleArray[3].Add(MCtrk);
       }

      //Fills for all identified protos found after cuts (reconstructed) - numerator for Efficiency
   //******************************

     //get coresponding MC particle 
     // Int_t label = TMath::Abs(track->GetLabel()); //moved up
     // if(!label) cout<<"no label"<<endl;
     //if(label) cout<<"label = "<<label<<endl;
       
    //AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(label); //moved up
    if (!MCtrk) continue;
    if(MCtrk->Charge()==0){cout<<"!!!"<<endl; continue;}
    recoParticleArray[0].Add(MCtrk);


    //Fills histogram for particles that are contamination from secondaries:
    if (!AODmcTrack->IsPhysicalPrimary()) {
      fReconstructedNotPrimaries[PARTTYPES*fcent]->Fill(track->Eta(), track->Pt());
    }
    else{
      fReconstructedPrimaries[PARTTYPES*fcent]->Fill(track->Eta(), track->Pt());
    }
 
    int PDGcode = MCtrk->GetPdgCode();

   //And secondaries for different particle species:
    if (!AODmcTrack->IsPhysicalPrimary() && (IsPionNSigma(track->P(),nSigmaTPCPi, nSigmaTOFPi) && PDGcode==211)) { //secondaries in pions
      fReconstructedNotPrimaries[PARTTYPES*fcent+1]->Fill(track->Eta(), track->Pt());
    }
    else if(AODmcTrack->IsPhysicalPrimary() && (IsPionNSigma(track->P(),nSigmaTPCPi, nSigmaTOFPi) && PDGcode==211)) {
      fReconstructedPrimaries[PARTTYPES*fcent+1]->Fill(track->Eta(), track->Pt());
    }
    
    if (!AODmcTrack->IsPhysicalPrimary() && (IsKaonNSigma(track->P(),nSigmaTPCK, nSigmaTOFK) && PDGcode==321)) { //secondaries in kaons
      fReconstructedNotPrimaries[PARTTYPES*fcent+2]->Fill(track->Eta(), track->Pt());
    }
    else if(AODmcTrack->IsPhysicalPrimary() && (IsKaonNSigma(track->P(),nSigmaTPCK, nSigmaTOFK) && PDGcode==321)) {
      fReconstructedPrimaries[PARTTYPES*fcent+2]->Fill(track->Eta(), track->Pt());
      }

    if (!AODmcTrack->IsPhysicalPrimary() && (IsProtonNSigma(track->P(),nSigmaTPCP, nSigmaTOFP) && PDGcode==2212)) { //secondaries in protons
      fReconstructedNotPrimaries[PARTTYPES*fcent+3]->Fill(track->Eta(), track->Pt());
    } 
    else if(AODmcTrack->IsPhysicalPrimary() && (IsProtonNSigma(track->P(),nSigmaTPCP, nSigmaTOFP) && PDGcode==2212)) {
      fReconstructedPrimaries[PARTTYPES*fcent+3]->Fill(track->Eta(), track->Pt());
    } 

    //Misidentification fraction
    if(PDGcode==211)
      {
	if(IsPionNSigma(track->P(),nSigmaTPCPi, nSigmaTOFPi))
	  fMisidentification[fcent]-> Fill(1,0.5);
	if(IsKaonNSigma(track->P(),nSigmaTPCK, nSigmaTOFK))
	  fMisidentification[fcent]-> Fill(1,1.5);
	if(IsProtonNSigma(track->P(),nSigmaTPCP, nSigmaTOFP))
	  fMisidentification[fcent]-> Fill(1,2.5);
	if(!IsPionNSigma(track->P(),nSigmaTPCPi, nSigmaTOFPi) && !IsKaonNSigma(track->P(),nSigmaTPCK, nSigmaTOFK) && !IsProtonNSigma(track->P(),nSigmaTPCP, nSigmaTOFP))
	  fMisidentification[fcent]-> Fill(1,3.5);


      }
    else if(PDGcode==321)
      {
	if(IsPionNSigma(track->P(),nSigmaTPCPi, nSigmaTOFPi))
	  fMisidentification[fcent]-> Fill(2,0.5);
	if(IsKaonNSigma(track->P(),nSigmaTPCK, nSigmaTOFK))
	  fMisidentification[fcent]-> Fill(2,1.5);
	if(IsProtonNSigma(track->P(),nSigmaTPCP, nSigmaTOFP))
	  fMisidentification[fcent]-> Fill(2,2.5);
	if(!IsPionNSigma(track->P(),nSigmaTPCPi, nSigmaTOFPi) && !IsKaonNSigma(track->P(),nSigmaTPCK, nSigmaTOFK) && !IsProtonNSigma(track->P(),nSigmaTPCP, nSigmaTOFP))
	  fMisidentification[fcent]-> Fill(2,3.5);


      }
    else if(PDGcode == 2212)
      {
	if(IsPionNSigma(track->P(),nSigmaTPCPi, nSigmaTOFPi))
	  fMisidentification[fcent]-> Fill(3,0.5);
	if(IsKaonNSigma(track->P(),nSigmaTPCK, nSigmaTOFK))
	  fMisidentification[fcent]-> Fill(3,1.5);
	if(IsProtonNSigma(track->P(),nSigmaTPCP, nSigmaTOFP))
	  {
	  fMisidentification[fcent]-> Fill(3,2.5);
	  }
	if(!IsPionNSigma(track->P(),nSigmaTPCPi, nSigmaTOFPi) && !IsKaonNSigma(track->P(),nSigmaTPCK, nSigmaTOFK) && !IsProtonNSigma(track->P(),nSigmaTPCP, nSigmaTOFP))
	  fMisidentification[fcent]-> Fill(3,3.5);
      }


    //Contaminations: "how many pions are in the kaons sample"? etc.
    //Do not use for corrections: using those values will be dependant on i.e. Pi/K ratio in MC
    //Use misidentification fraction instead
    if(IsPionNSigma(track->P(),nSigmaTPCPi, nSigmaTOFPi))
      {
	fContamination[PARTTYPES*fcent+1]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for pions
      }
    if(IsKaonNSigma(track->P(),nSigmaTPCK, nSigmaTOFK))
      {
	fContamination[PARTTYPES*fcent+2]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for kaons
      }
    if(IsProtonNSigma(track->P(),nSigmaTPCP, nSigmaTOFP))
      {
	fContamination[PARTTYPES*fcent+3]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for protons
      }
    


  } 


  // MONTECARLO PARTICLES 
  if(!arrayMC){
    AliError("Array of MC particles not found");
    return;
  }
  // loop over MC stack 
  for (Int_t ipart = 0; ipart < arrayMC->GetEntries(); ipart++) {
    //cout<<"Entered MC loop"<<endl;
    
    AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(ipart);

    if (!MCtrk) continue;
    //cout<<"particle obtained"<<endl;
    
    Int_t PDGcode = TMath::Abs(MCtrk->GetPdgCode()); 

    
    if(MCtrk->Charge() == 0) continue;
     
    //*** PID - check if pion ***
    //if(PDGcode!=211) continue; //(PDGcode==11 || PDGcode==321 || PDGcode==2212 || PDGcode==13)

      if(MCtrk->Eta() < -0.8 || MCtrk->Eta() > 0.8){
	continue; }
	
      if (MCtrk->Pt() < 0.2 || MCtrk->Pt() > 20){
	continue;}
      	
      // check physical primary 
      if (MCtrk->IsPhysicalPrimary()){

	// Filling histograms for MC truth particles
	fGeneratedMCPrimaries[fcent*PARTTYPES]->Fill(MCtrk->Eta(), MCtrk->Pt());
	if(PDGcode==211)
	  fGeneratedMCPrimaries[fcent*PARTTYPES+1]->Fill(MCtrk->Eta(), MCtrk->Pt());
	else if(PDGcode==321)
	  fGeneratedMCPrimaries[fcent*PARTTYPES+2]->Fill(MCtrk->Eta(), MCtrk->Pt());
	else if(PDGcode==2212)
	  fGeneratedMCPrimaries[fcent*PARTTYPES+3]->Fill(MCtrk->Eta(), MCtrk->Pt());

	  //Filling data from MC truth particles only for particles that were reconstruced
	if (recoParticleArray[0].Contains(MCtrk)){ //All
	  fMCPrimariesThatAreReconstructed[fcent*PARTTYPES]->Fill(MCtrk->Eta(), MCtrk->Pt());
	}
	if (recoParticleArray[1].Contains(MCtrk)){ //Pions
	  if(PDGcode==211)
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+1]->Fill(MCtrk->Eta(), MCtrk->Pt());
	}
	if (recoParticleArray[2].Contains(MCtrk)){ //Kaons
	  if(PDGcode==321)
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+2]->Fill(MCtrk->Eta(), MCtrk->Pt());
	}
	if (recoParticleArray[3].Contains(MCtrk)){ //Protons
	  if(PDGcode==2212)
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+3]->Fill(MCtrk->Eta(), MCtrk->Pt());
	}

      }
    
  }
  PostData(1, fHistoList);
}
//-----------------------------------------------------------------

//void AliAnalysisTaskParticleEfficiency::Terminate(Option_t *) 
//{}
