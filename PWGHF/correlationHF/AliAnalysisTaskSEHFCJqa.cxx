/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// Class AliAnalysisTaskSEHFCJqa
// AliAnalysisTaskSE for the extraction of the fraction of prompt charm
// using the charm hadron impact parameter to the primary vertex
//
// Author: Andrea Rossi, andrea.rossi@pd.infn.it
/////////////////////////////////////////////////////////////


#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TAxis.h>
#include <THnSparse.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TROOT.h>
#include "AliAODEvent.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecay.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliESDtrack.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliAODMCParticle.h"
#include "AliAODPid.h"
#include "AliTPCPIDResponse.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSEHFCJqa.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliNormalizationCounter.h"

class TCanvas;
class TTree;
class TChain;
class AliAnalysisTaskSE;


ClassImp(AliAnalysisTaskSEHFCJqa)

  AliAnalysisTaskSEHFCJqa::AliAnalysisTaskSEHFCJqa() 
: AliAnalysisTaskSE(),
  fReadMC(),
  ffilterbit(),
  fKeepTrackNegID(),
  fpidResp(),
  fCuts(),
  fhEventCounter(),
  fhImpParResolITSsel(),
  fhImpParResolITSselGoodTracks(),
  fhSparseFilterMask(),
  fhSparseFilterMaskTrackAcc(),
  fhSparseFilterMaskImpPar(),
  fhSparseEoverPeleTPC(),
  fhSparseShowShapeEleTPC(),
  fhnSigmaTPCTOFEle(),
  fhnSigmaTPCTOFPion(),
  fhnSigmaTPCTOFKaon(),
  fhnSigmaTPCTOFProton(),
  fhTrackEMCal(),
  fSparseRecoJets(),
  fLoadJet(),
  fJetArrayString(),
  fListTrackAndPID(),
  fListJets()
{// standard constructor
  
}



//________________________________________________________________________
AliAnalysisTaskSEHFCJqa::AliAnalysisTaskSEHFCJqa(const char *name) 
  : AliAnalysisTaskSE(name),
    fReadMC(),
    ffilterbit(AliAODTrack::kTrkGlobalNoDCA),
    fKeepTrackNegID(),
    fpidResp(),
    fCuts(),
    fhEventCounter(),
    fhImpParResolITSsel(),
    fhImpParResolITSselGoodTracks(),
    fhSparseFilterMask(),
    fhSparseFilterMaskTrackAcc(),
    fhSparseFilterMaskImpPar(),
    fhSparseEoverPeleTPC(),
    fhSparseShowShapeEleTPC(),
    fhnSigmaTPCTOFEle(),
    fhnSigmaTPCTOFPion(),
    fhnSigmaTPCTOFKaon(),
    fhnSigmaTPCTOFProton(),
    fhTrackEMCal(),
    fSparseRecoJets(),
    fLoadJet(),
    fJetArrayString(),
    fListTrackAndPID(),
    fListJets()
{ // default constructor
  
  DefineOutput(1, TH1F::Class()); // counter
  DefineOutput(2, TList::Class()); // single track properties list and PID
  DefineOutput(3, TList::Class()); // jet properties list

}

//________________________________________________________________________
AliAnalysisTaskSEHFCJqa::~AliAnalysisTaskSEHFCJqa(){
  // destructor
  delete fCuts;
  delete fpidResp;
  delete fhEventCounter;
  delete fhSparseFilterMask;
  delete   fhSparseFilterMaskTrackAcc;
  delete fhSparseFilterMaskImpPar;
  delete fhSparseEoverPeleTPC;
  delete fhSparseShowShapeEleTPC;
  delete fhnSigmaTPCTOFEle;
  delete fhnSigmaTPCTOFPion;
  delete fhnSigmaTPCTOFKaon;
  delete fhnSigmaTPCTOFProton;
  delete fhTrackEMCal;
  delete fSparseRecoJets;
  delete fhImpParResolITSsel;
  delete fhImpParResolITSselGoodTracks;
  delete fListTrackAndPID;
  delete fListJets;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFCJqa::Init()
{
  // Initialization

}


//________________________________________________________________________
void AliAnalysisTaskSEHFCJqa::UserCreateOutputObjects(){


  //##########  DEFINE THE TLISTS ##################
  fListTrackAndPID=new TList();
  fListTrackAndPID->SetOwner();
  fListTrackAndPID->SetName("fListTrackAndPID");

  fListJets=new TList();
  fListJets->SetOwner();
  fListJets->SetName("fListJets");

  
  //  ########### DEFINE THE EVENT COUNTER ############
  fhEventCounter=new TH1F("fhEventCounter","Counter of event selected",20,-0.5,19.5);
  fhEventCounter->GetXaxis()->SetBinLabel(1,"Events analyzed");
  fhEventCounter->GetXaxis()->SetBinLabel(2,"Event selected");
  fhEventCounter->GetXaxis()->SetBinLabel(3,"Jet array present");
  fhEventCounter->GetXaxis()->SetBinLabel(4,"Vtx Track Ncont");

  
  
  Int_t nbinsRecoJets[8]={50,20,20,20,5,5,10,60};
  Double_t binlowRecoJets[8]={5.,-1.,-TMath::Pi(),0.99,0.,-0.5,0,0.};
  Double_t binupRecoJets[8]={55.,1.,TMath::Pi(),20.99,4.99,4.5,2.,60.};
  fSparseRecoJets=new THnSparseF("fSparseRecoJets","fSparseRecoJets;jetpt;eta;phi;ntrks;nEle;parton;partContr;ptPart;",8,nbinsRecoJets,binlowRecoJets,binupRecoJets);

  fListJets->Add(fSparseRecoJets);

  // Num axes: 10 filter bits + ID + TPCrefit,ITSrefit,bothTPCandITSrefit + kAny,kFirst,kBoth+ 20Nclust TPC+ 20 NTPC crossed padRows + 20 DCA
  Int_t nbinsFilterMask[16]={2,2,2,2,2,2,2,2,2,2,2,3,3,20,20,70};
  Double_t binlowFilterMask[16]={-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,0,0,-3.5};
  Double_t binupFilterMask[16]={1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,2.5,2.5,160,160,3.5};
  fhSparseFilterMask=new THnSparseF("fhSparseFilterMask","fhSparseFilterMask;fb0;fb1;fb2;fb3;fb4;fb5;fb6;fb7;fb8;fb9;trackID;refitting;SPD;NTPCclust;NTPCcrossRows;DCA",16,nbinsFilterMask,binlowFilterMask,binupFilterMask);
  fListTrackAndPID->Add(fhSparseFilterMask);


// Num axes: 10 filter bits + ID*isSelected 5 + kAny,kFirst,kBoth + phi+ eta +pt 
  Int_t nbinsFilterMaskTrackAcc[15]={2,2,2,2,2,2,2,2,2,2,5,3,36,30,30};
  Double_t binlowFilterMaskTrackAcc[15]={-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-2.5,-0.5,0.,-1.5,0.};
  Double_t binupFilterMaskTrackAcc[15]={1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,2.5,2.5,TMath::Pi()*2.,1.5,15.};
  fhSparseFilterMaskTrackAcc=new THnSparseF("fhSparseFilterMaskTrackAcc","fhSparseFilterMaskTrackAcc;fb0;fb1;fb2;fb3;fb4;fb5;fb6;fb7;fb8;fb9;trackIDandsel;SPD;#phi (rad);#eta;p_{T} (GeV/c)",15,nbinsFilterMaskTrackAcc,binlowFilterMaskTrackAcc,binupFilterMaskTrackAcc);
  fListTrackAndPID->Add(fhSparseFilterMaskTrackAcc);  


// Num axes: ID*isSelected 5 + kAny,kFirst,kBoth + phi+ eta +pt + imp par
  Int_t nbinsFilterMaskImpPar[6]={5,3,36,30,30,200};
  Double_t binlowFilterMaskImpPar[6]={-2.5,-0.5,0.,-1.5,0.,-300.};
  Double_t binupFilterMaskImpPar[6]={2.5,2.5,TMath::Pi()*2.,1.5,15.,300.};
  fhSparseFilterMaskImpPar=new THnSparseF("fhSparseFilterMaskImpPar","fhSparseFilterMaskImpPar;fb0;fb1;fb2;fb3;fb4;fb5;fb6;fb7;fb8;fb9;trackIDandsel;SPD;#phi (rad);#eta;p_{T} (GeV/c; imp par (#mum)",6,nbinsFilterMaskImpPar,binlowFilterMaskImpPar,binupFilterMaskImpPar);
  fListTrackAndPID->Add(fhSparseFilterMaskImpPar);  

  
fhImpParResolITSsel=new TH3F("fhImpParResolITSsel","fhImpParResolITSsel;p_{T} (GeV/c);imp. par (#mum);ITS clust",50.,0.,10.,200,-800.,800.,38,-0.5,37.5);
  // the convention for ITS clust axis: number between 0 and 48, first digit (units) = number of clust in ITS, second digits (x10): SPD status: 0 = none, 1=kFirst, 2=kSecond; 3= kBoth 
  fListTrackAndPID->Add(fhImpParResolITSsel);

  fhImpParResolITSselGoodTracks=new TH3F("fhImpParResolITSselGoodTracks","fhImpParResolITSselGoodTracks;p_{T} (GeV/c);imp. par (#mum);ITS clust",50.,0.,10.,200,-800.,800.,38,-0.5,37.5);
  // the convention for ITS clust axis: number between 0 and 48, first digit (units) = number of clust in ITS, second digits (x10): SPD status: 0 = none, 1=kFirst, 2=kSecond; 3= kBoth 
  fListTrackAndPID->Add(fhImpParResolITSselGoodTracks);


  // PID PLOTS
  fhnSigmaTPCTOFEle=new TH3F("fhnSigmaTPCTOFEle","fhnSigmaTPCTOFEle;p (GeV/c);n#sigma_{TPC}^{ele};n#sigma_{TOF}^{ele}",120,0.,30,64,-8.,8.,64,-8.,8.);
  fhnSigmaTPCTOFPion=new TH3F("fhnSigmaTPCTOFPion","fhnSigmaTPCTOFPion;p (GeV/c);n#sigma_{TPC}^{ele};n#sigma_{TOF}^{ele}",120,0.,30,64,-8.,8.,64,-8.,8.);
  fhnSigmaTPCTOFKaon=new TH3F("fhnSigmaTPCTOFKaon","fhnSigmaTPCTOFKaon;p (GeV/c);n#sigma_{TPC}^{ele};n#sigma_{TOF}^{ele}",120,0.,30,64,-8.,8.,64,-8.,8.);
  fhnSigmaTPCTOFProton=new TH3F("fhnSigmaTPCTOFProton","fhnSigmaTPCTOFProton;p (GeV/c);n#sigma_{TPC}^{ele};n#sigma_{TOF}^{ele}",120,0.,30,64,-8.,8.,64,-8.,8.);

  fListTrackAndPID->Add(fhnSigmaTPCTOFEle);
  fListTrackAndPID->Add(fhnSigmaTPCTOFPion);
  fListTrackAndPID->Add(fhnSigmaTPCTOFKaon);
  fListTrackAndPID->Add(fhnSigmaTPCTOFProton);


  /// EMCAL PLOTS
  //study of NsigmaTPC, e/p, p
  Int_t nbinsEoP[4]={202,400,100,400};
  Double_t binlowEoP[4]= {-1., -20.,-20., -1};
  Double_t binupEoP[4]= {100., 20, 20.,9};
  fhSparseEoverPeleTPC = new THnSparseF("fhSparseEoP", "fhSparseEoP; p;nsigmatpc;nsigmaElePIDresp;E/p;",4, nbinsEoP, binlowEoP, binupEoP);
  fListTrackAndPID->Add(fhSparseEoverPeleTPC);
  //fSparseEoverPallHadr = new THnSparseF("fSparseEoP", "fSparseEoP; p;nsigmatpc;E/p;",3, nbinsEoP, binlowEoP, binupEoP);

  Int_t nbinsEmShSh[7]={35,120,100,50,50,50,50};
  Double_t binlowEmShSh[7]= {5., -20., -1,0,0,0,0};
  Double_t binupEmShSh[7]= {40., 10, 9,1,1,2,50};
  fhSparseShowShapeEleTPC = new THnSparseF("fhSparseShowShapeEleTPC", "fhSparseShowShapeEleTPC; pt;nsigmatpc;E/p;M02;M20;disp;Nclust",7, nbinsEmShSh, binlowEmShSh, binupEmShSh);
  fListTrackAndPID->Add(fhSparseShowShapeEleTPC);

  
  
  Int_t nbinsTrEM[6]={124,20,124,25,20,20};
  Double_t binlowTrEM[6]={-1,-1.,-1.0,0.,0.};
  Double_t binupTrEM[6]={31,1.,31.,5.,1.,1.};

  fhTrackEMCal=new THnSparseF("fhTrackEMCal","fhTrackEMCal;p(GeV/c);eta;clusterE(GeV/c);E/p;trkDistX;trkDistZ",6,nbinsTrEM,binlowTrEM,binupTrEM);
  
  

  //fhPhotonicEle=new TH3F("fhPhotonicEle","fhPhotonicEle; pele;pgamma;mass",20,0.,20.,20,0.,20.,50,0.,5);
  //fhPhotonicEleLS=new TH3F("fhPhotonicEleLS","fhPhotonicEleLS; pele;pgamma;mass",20,0.,20.,20,0.,20.,50,0.,5);


  PostData(1,fhEventCounter);
  PostData(2,fListTrackAndPID);
  PostData(3,fListJets);
  
}



//________________________________________________________________________
void AliAnalysisTaskSEHFCJqa::UserExec(Option_t */*option*/){



// Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if(aod){
    fhEventCounter->Fill(0);

    if(!fCuts->IsEventSelected(aod)){
      PostData(1,fhEventCounter);
      return;
    }
    fhEventCounter->Fill(1);
  }

  TClonesArray *arrayJets=0x0;
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    if (!aod) {
      Printf("ERROR: aod not available");
      return;
    }
    else {
      fhEventCounter->Fill(0);

      if(!fCuts->IsEventSelected(aod)){
	PostData(1,fhEventCounter);
	return;
      }
      fhEventCounter->Fill(1);
    }

    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());

    if(fLoadJet>=1&&aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.Jets.root");
      AliAODEvent* aodFromExt = ext->GetAOD();
      
      // load jet array
      arrayJets = (TClonesArray*)aodFromExt->GetList()->FindObject(fJetArrayString.Data());
      if(!arrayJets) {
	Printf("AliAnalysisTaskSEHFCJqa::UserExec: desired branch %s not found!",fJetArrayString.Data());
	PostData(1,fhEventCounter);
	return;

      }
        
    }
  } else {
    // load jet array
    if(fLoadJet>=1){
      arrayJets = (TClonesArray*)aod->GetList()->FindObject(fJetArrayString.Data());
      if(!arrayJets) {
	Printf("AliAnalysisTaskSEHFCJqa::UserExec: desired Jet branch %s not found!",fJetArrayString.Data());
	PostData(1,fhEventCounter);
	return;
      }
      
    }
  }

  if(fLoadJet!=0 && !arrayJets) {
    printf("AliAnalysisTaskSEHFCJqa::UserExec: desired jet input branch not found!\n");
    PostData(1,fhEventCounter);
    return;
  }

  fhEventCounter->Fill(2);
  
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  TString primTitle = vtx1->GetTitle();
    if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) { 
    fhEventCounter->Fill(3);

  }
  else {
    PostData(1,fhEventCounter);
    return;
    
  }
    // Convert primary vertex to esd vertex (needed for later usage of RelateToVertex)
    Double_t pos[3],cov[6];
    vtx1->GetXYZ(pos);
    vtx1->GetCovarianceMatrix(cov);
    const AliESDVertex vESD(pos,cov,100.,100);
    Double_t magfield=aod->GetMagneticField();

  TClonesArray *arrayMC=0x0;
  AliAODMCHeader *aodmcHeader=0x0;
  Double_t vtxTrue[3];
 
  
  if(fReadMC){
    // load MC particles
    arrayMC = 
      (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      Printf("AliAnalysisTaskSEHFCJqa::UserExec: MC particles branch not found!\n");
      return;
    }
    // load MC header
    aodmcHeader = 
      (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!aodmcHeader) {
      Printf("AliAnalysisTaskSEHFCJqa::UserExec: MC header branch not found!\n");
      return;
    }
    // MC primary vertex
    aodmcHeader->GetVertex(vtxTrue);

  }

  // Starting the fun part
  SetupPIDresponse();// this sets the pid reponse to fPPIDRespons; could also get from the cut object (AliAODPidHF::GetPidResponse)

  // Looping over aod tracks
  for(Int_t j=0;j<aod->GetNumberOfTracks();j++){

    AliAODTrack *aodtrack=dynamic_cast<AliAODTrack*>(aod->GetTrack(j));
    if(!aodtrack) AliFatal("Not a standard AOD");
    // CHECK FILTER MAPS
    if(!FillTrackHistosAndSelectTrack(aodtrack,&vESD,magfield))continue;
    //    if(j%100==0)  
  
    Double_t p=aodtrack->P();
    Double_t eta=aodtrack->Eta();
    // START PID: TPC

    Double_t nsigmaEleTPC=fpidResp->NumberOfSigmasTPC(aodtrack, AliPID::kElectron);
    Double_t nsigmaPionTPC=fpidResp->NumberOfSigmasTPC(aodtrack, AliPID::kPion);
    Double_t nsigmaKaonTPC=fpidResp->NumberOfSigmasTPC(aodtrack, AliPID::kKaon);
    Double_t nsigmaProtonTPC=fpidResp->NumberOfSigmasTPC(aodtrack, AliPID::kProton);


    // TOF
    Double_t nsigmaEleTOF=fpidResp->NumberOfSigmasTOF(aodtrack, AliPID::kElectron);
    Double_t nsigmaPionTOF=fpidResp->NumberOfSigmasTOF(aodtrack, AliPID::kPion);
    Double_t nsigmaKaonTOF=fpidResp->NumberOfSigmasTOF(aodtrack, AliPID::kKaon);
    Double_t nsigmaProtonTOF=fpidResp->NumberOfSigmasTOF(aodtrack, AliPID::kProton);


    fhnSigmaTPCTOFEle->Fill(p,nsigmaEleTPC,nsigmaEleTOF);
    fhnSigmaTPCTOFPion->Fill(p,nsigmaPionTPC,nsigmaPionTOF);
    fhnSigmaTPCTOFKaon->Fill(p,nsigmaKaonTPC,nsigmaKaonTOF);
    fhnSigmaTPCTOFProton->Fill(p,nsigmaProtonTPC,nsigmaProtonTOF);
    //    if(j%100==0)


    // NOW EMCAL
    // CHECK WHETHER THERE IS A EMCAL CLUSTER
    Int_t nClsId = aodtrack->GetEMCALcluster();
    if(nClsId >0) {
      AliVCluster *cluster = aod->GetCaloCluster(nClsId);
      Double_t clsE = cluster->E();
      if(!cluster->IsEMCAL())       continue;
      
      // Check whether it is an electron candidate: do not reject here hadrons to allow QA checks of (E/p,nsigmaTPC,pt)
      Double_t nEoverP = clsE/p;
      Double_t eOverPpidResp;
      Double_t showerShape[4];
      Double_t nsigmaEleEMCal=fpidResp->NumberOfSigmasEMCAL(aodtrack,AliPID::kElectron,eOverPpidResp,showerShape);
      Double_t poix[4]={p, nsigmaEleTPC, nsigmaEleEMCal, nEoverP};
      fhSparseEoverPeleTPC->Fill(poix);
      
      
      Double_t emcTrackDx=cluster->GetTrackDx();
      Double_t emcTrackDz=cluster->GetTrackDz();
      Double_t pointET[6]={p,eta,clsE,nEoverP,emcTrackDx,emcTrackDz};
      fhTrackEMCal->Fill(pointET);    
      
      Double_t pointEmShSh[7]={aodtrack->Pt(), static_cast<Double_t>(nsigmaEleTPC), static_cast<Double_t>(nEoverP),cluster->GetM02(),cluster->GetM20(),cluster->GetDispersion(), static_cast<Double_t>(cluster->GetNCells())};

      fhSparseShowShapeEleTPC->Fill(pointEmShSh);

    }
  }
  
  
  // NOW LOOP OVER JETS
  
  Int_t nJets=arrayJets->GetEntries();//calcolo numero di jet nell'event
  AliAODJet *jet;
  for(Int_t jetcand=0;jetcand<nJets;jetcand++){
    //    if(jetcand%100==0)

    jet=(AliAODJet*)arrayJets->UncheckedAt(jetcand);
       
    Double_t contribution=0,ptpart=-1;
    Int_t partonnat=0;
    if(fReadMC){
      
      AliAODMCParticle *parton=IsMCJet(arrayMC,jet,contribution);
      if(parton){
	Int_t pdg=TMath::Abs(parton->PdgCode());
	//printf("pdg parton: %d \n",pdg);
	if(pdg==21)partonnat=1;
	else if(pdg<4)partonnat=2;
	else if(pdg==4)partonnat=3;
	else if(pdg==5)partonnat=4;
	ptpart=parton->Pt();
      }
      
    }
    
    
    FillJetRecoHisto(jet,partonnat,contribution,ptpart);
  }

  

  PostData(1,fhEventCounter);
  PostData(2,fListTrackAndPID);
  PostData(3,fListJets);

  return;
}



void AliAnalysisTaskSEHFCJqa::SetupPIDresponse(){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler(); 
  fpidResp=inputHandler->GetPIDResponse();
  if(!fpidResp)AliFatal("No PID response could be set");
}

//_______________________________________________________________
Bool_t AliAnalysisTaskSEHFCJqa::FillTrackHistosAndSelectTrack(AliAODTrack *aodtr, const AliESDVertex *primary, const Double_t magfield){
  
  Bool_t retval=kTRUE;
  // THnSparse for filter bits
  Double_t point[16]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,999.,-999.,-999.,-999.,-999.,-999.};
  Double_t pointAcc[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,999.,-999.,-999.,-999.,-999.};
  Double_t pointImpPar[6]={999.,-999.,-999.,-999.,-999.,-999.};

  for(Int_t j=0;j<10;j++){
    if(aodtr->TestFilterBit(TMath::Power(2,j))){
      point[j]=1;      
      pointAcc[j]=1;
    }
  }
  // check ID
  Int_t trid=aodtr->GetID();
  if(aodtr->GetID()>0)point[10]=1.;
  else if(aodtr->GetID()>0)point[10]=0.;
  if(aodtr->GetID()==0)point[10]=-1.;

  Float_t iparxy,iparz;
  Float_t pt=aodtr->Pt();
  Int_t clustITS=aodtr->GetITSNcls();// for histo
  AliESDtrack esdtrack(aodtr);
  // needed to calculate impact parameters
  

  // check refit status
  Int_t refit=-1;
  ULong64_t status = esdtrack.GetStatus();
  if(status&AliESDtrack::kTPCrefit)refit+=1;
  if(status&AliESDtrack::kITSrefit)refit+=2;
  point[11]=refit;
  // CHECK SPD
  point[12]=-1;
  if(aodtr->HasPointOnITSLayer(0)){
    clustITS+=10;
    point[12]+=1;
  }
  if(aodtr->HasPointOnITSLayer(1)){
    point[12]+=2;
    clustITS+=20;
  }
  point[13]=aodtr->GetTPCNcls();
  point[14]=aodtr->GetTPCNCrossedRows();
  esdtrack.RelateToVertex(primary,magfield,4.);// CHECK THIS : I put 4.. usually we set it to 3 
  esdtrack.GetImpactParameters(iparxy,iparz);
  point[15]=iparxy;



  
  fhSparseFilterMask->Fill(point);
  if(!aodtr->TestBit(ffilterbit))retval =kFALSE;

  if(aodtr->GetID()<0&&!fKeepTrackNegID)retval = kFALSE;
  if(retval)  fhImpParResolITSsel->Fill(pt,iparxy*10000.,clustITS);

  AliESDtrackCuts *cuts=fCuts->GetTrackCuts();
  if(!cuts->IsSelected(&esdtrack)) retval = kFALSE;

  if(fCuts->GetUseKinkRejection()){
    AliAODVertex *maybeKink=aodtr->GetProdVertex();
    if(maybeKink->GetType()==AliAODVertex::kKink) retval=kFALSE;
  }

  if(retval){
    pointAcc[10]=1*trid;
  }
  else {
    if(trid!=0)
    pointAcc[10]=2*trid;
    else pointAcc[10]=-3;
  }

  pointAcc[11]=point[12];
  pointAcc[12]=aodtr->Phi();
  if(pointAcc[12]<0.)pointAcc[12]+=2.*TMath::Pi();    
  pointAcc[13]=aodtr->Eta();
  pointAcc[14]=pt;
  fhSparseFilterMaskTrackAcc->Fill(pointAcc);

  pointImpPar[0]=pointAcc[10];
  pointImpPar[1]=pointAcc[11];
  pointImpPar[2]=pointAcc[12];
  pointImpPar[3]=pointAcc[13];
  pointImpPar[4]=pointAcc[14];
  pointImpPar[5]=iparxy;
  fhSparseFilterMaskImpPar->Fill(pointImpPar);
  
  if(retval)  {
    fhImpParResolITSselGoodTracks->Fill(pt,iparxy*10000.,clustITS);  
  }
  
  return retval;
  
}


//---------------------------------------------------------------
AliAODMCParticle* AliAnalysisTaskSEHFCJqa::IsMCJet(TClonesArray *arrayMC,const AliAODJet *jet, Double_t &contribution){// assignment of parton ID to jet
  // method by L. Feldkamp
  std::vector< int >           idx;
  std::vector< int >           idx2;
  std::vector< double >     weight;

  int counter =0;
  int num = jet->GetRefTracks()->GetEntries();
  
  for(int k=0;k<num;++k){
    
    AliAODTrack * track = (AliAODTrack*)(jet->GetRefTracks()->At(k));
    if (track->GetLabel() >=0){
      AliAODMCParticle* part =  (AliAODMCParticle*)  arrayMC->At(track->GetLabel());
      if(!part)continue;

      int label =0 ;
      AliAODMCParticle* motherParton=GetMCPartonOrigin(arrayMC,part, label);
      if (!motherParton) Printf("no mother");
      else {
	counter++;
	idx.push_back(label);                       //! Label  of Mother
	idx2.push_back(label);        
	weight.push_back(track->Pt());  //! Weight : P_t trak /  P_t jet ... the latter used at the end
      }
    }///END LOOP OVER REFERENCED TRACKS   
  }
  //! Remove duplicate indices for counting
  sort( idx2.begin(), idx2.end() );
  idx2.erase( unique( idx2.begin(), idx2.end() ), idx2.end() );
  if (idx2.size() == 0) return 0x0;
  Double_t* arrayOfWeights = new Double_t[(UInt_t)idx2.size()];
  if(!arrayOfWeights){
    return 0x0;
  }
  for(UInt_t ii=0;ii<(UInt_t)idx2.size();ii++)arrayOfWeights[ii]=0;

  for (unsigned int idxloop =0 ;idxloop<idx2.size();idxloop++){
    for (unsigned int z=0; z< idx.size() ; ++z){
      int     a = idx.at(z);
      double w = weight.at(z);
      if(a == idx2.at(idxloop))    arrayOfWeights[idxloop] += w;
    }
  }
  
  int winner = -1;
  double c=-1.;
  for (unsigned int idxloop =0 ;idxloop<idx2.size();idxloop++){
    if(c < arrayOfWeights[idxloop]){
      winner =idxloop; 
      c=arrayOfWeights[idxloop];
    }
  }
  
  AliAODMCParticle *parton = 0x0;
  if(winner>0){
    parton=(AliAODMCParticle*)arrayMC->At(idx.at(winner));
    contribution = arrayOfWeights[winner]/jet->Pt();
  }
  else {
  
    if(arrayOfWeights)    delete [] arrayOfWeights;
    return 0x0;
  }
  if(arrayOfWeights)    delete [] arrayOfWeights;

  return parton;
  
  
}

//---------------------------------------------------------------
AliAODMCParticle *AliAnalysisTaskSEHFCJqa::GetMCPartonOrigin(TClonesArray *arrayMC,AliAODMCParticle *p, Int_t &idx)
{  //Follows chain of track mothers until q/g or idx = -1	
  AliAODMCParticle *p2=0x0;
  Int_t mlabel = TMath::Abs(p->GetMother()) ; 
  Double_t pz=0.;
  while(mlabel > 1){
    p2 = (AliAODMCParticle*)arrayMC->At(mlabel);
    pz=TMath::Abs(p2->Pz());
    //printf("Mother label %d, pdg %d, pz %f\n",mlabel,p2->PdgCode(),pz);
    if( p2->PdgCode() == 21 || (p2->PdgCode() != 0 && abs(p2->PdgCode()) <6) )
      {
	idx = mlabel; 
	return p2;
      }
    mlabel = TMath::Abs(p2->GetMother()); 
  }
  idx=-1;
  return 0x0;

} 

//_____________________________________________________________
void AliAnalysisTaskSEHFCJqa::FillJetRecoHisto(const AliAODJet *jet,Int_t partonnat,Double_t contribution,Double_t ptpart){
//FIll sparse with reco jets properties
  Double_t point[8]={jet->Pt(),jet->Eta(),jet->Phi()-TMath::Pi(), static_cast<Double_t>(jet->GetRefTracks()->GetEntriesFast()),0, static_cast<Double_t>(partonnat),contribution,ptpart};
  fSparseRecoJets->Fill(point);
}


//_______________________________________________________________
//void AliAnalysisTaskSEHFCJqa::FillTrackHistosPID(AliAODTrack *aodtr){
//
//}

//_______________________________________________________________
void AliAnalysisTaskSEHFCJqa::Terminate(const Option_t*){
  //TERMINATE METHOD: NOTHING TO DO



}

