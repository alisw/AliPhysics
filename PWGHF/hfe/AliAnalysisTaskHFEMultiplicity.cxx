
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


//////////////////////////////////////////////
//    Task for Measurement of Heavy Flavour //
//    electron as a function of charged     //	
//    particle multiplicity                 //
//    Author: Preeti Dhankher               //
//////////////////////////////////////////////


#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"
#include "TVector3.h"
#include <TRandom3.h>
#include "TProfile.h"
#include "TGeoManager.h"

#include "stdio.h"
#include "iostream"
#include "fstream"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAnalysisTaskHFEMultiplicity.h"
#include "TGeoGlobalMagField.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"
#include "AliEventPoolManager.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliESDpid.h"
#include "AliAODPid.h"

#include "AliMultSelection.h"
#include "AliMultiplicity.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliAODTracklets.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliESDCaloTrigger.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "AliCentrality.h"
#include "AliMagF.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include "AliESDVZERO.h"
#include "AliAODVZERO.h"
#include "AliESDUtils.h"
#include "AliNormalizationCounter.h"
#include "AliAnalysisVertexingHF.h"
#include "AliVertexingHFUtils.h"

#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"

class AliAnalysisTaskHFEMultiplicity;    
using namespace std;            
ClassImp(AliAnalysisTaskHFEMultiplicity)

AliAnalysisTaskHFEMultiplicity::AliAnalysisTaskHFEMultiplicity() : AliAnalysisTaskSE(), 
 
  fAOD(0),
  fMCarray(0),
  fMCparticle(0),
  fNevents(0),
  fTenderClusterName("caloClusters"),
  fTenderTrackName("tracks"),
  fOutputList(0),
  fListProfiles(0),
  fTracks_tender(0),
  fCaloClusters_tender(0),
  fUseTender(kTRUE),
  fEMCEG1(kFALSE),
  fEMCEG2(kFALSE),
  fDCalDG1(kFALSE),
  fDCalDG2(kFALSE),
  fFlagClsTypeEMC(kTRUE),
  fFlagClsTypeDCAL(kTRUE),
  fClusPhi(0),
  fClusEta(0),
  fClusEtaPhi(0x0),								 
  fClusE(0),								 
  fNCells(0),								 
  fClusT(0),								
  fCellE(0),
  fCellT(0),							 
  fpidResponse(0),
  fVtxZ(0),
  fVtxX(0),
  fVtxY(0),
  fTPCdEdx(0x0),
  fTPCnsigma(0x0),
  fTrkPt(0),
  fTrketa(0),
  fTrkphi(0),
  fTrkMatchTrkPt(0),
  fTrkMatchTrketa(0),
  fTrkMatchTrkphi(0),
  fTrkMatchClusetaphi(0x0),						 
  fEMCTrkMatchcluster(0x0),
  fReadMC(kFALSE),
 
  fRejectPUFromSPD(kFALSE),
  fRefMult(61.26),
  gRandom(new TRandom3(0)),
  fInvmassLS(0),
  fInvmassULS(0),
  fInvmassLSPt(0),
  fInvmassULSPt(0),
  fULSElecPt(0),
  fLSElecPt(0),

  fSparseElectron(0),
  fvalueElectron(0),
  fSparseLSElectron(0),
  fSparseULSElectron(0),
  fvaluePHElectron(0),
  fSparseMulti(0),
  fvalueMulti(0)
  
{   fvalueElectron = new Double_t[9];
  fvaluePHElectron = new Double_t[5];
  fvalueMulti = new Double_t[11];
  for(Int_t i=0; i<2; i++) fMultEstimatorAvg[i]=0;

}

//_____________________________________________________________________________
AliAnalysisTaskHFEMultiplicity::AliAnalysisTaskHFEMultiplicity(const char* name) : AliAnalysisTaskSE(name),
  
  fAOD(0),
  fMCarray(0),
  fMCparticle(0),
  fNevents(0),
  fTenderClusterName("caloClusters"),
  fTenderTrackName("tracks"),
  fOutputList(0),
  fListProfiles(0),
  fTracks_tender(0),
  fCaloClusters_tender(0),
  fUseTender(kTRUE),
  fEMCEG1(kFALSE),
  fEMCEG2(kFALSE),
  fDCalDG1(kFALSE),
  fDCalDG2(kFALSE),
  fFlagClsTypeEMC(kTRUE),
  fFlagClsTypeDCAL(kTRUE),
  fClusPhi(0),
  fClusEta(0),
  fClusEtaPhi(0x0),								 
  fClusE(0),								 
  fNCells(0),								 
  fClusT(0),								
  fCellE(0),
  fCellT(0),							 
  fpidResponse(0),
  fVtxZ(0),
  fVtxX(0),
  fVtxY(0),
  fTPCdEdx(0x0),
  fTPCnsigma(0x0),
  fTrkPt(0),
  fTrketa(0),
  fTrkphi(0),
  fTrkMatchTrkPt(0),
  fTrkMatchTrketa(0),
  fTrkMatchTrkphi(0),
  fTrkMatchClusetaphi(0x0),						 
  fEMCTrkMatchcluster(0x0),
  fReadMC(kFALSE),
 
  fRejectPUFromSPD(kFALSE),
  fRefMult(61.26),
  gRandom(new TRandom3(0)),
  fInvmassLS(0),
  fInvmassULS(0),
  fInvmassLSPt(0),
  fInvmassULSPt(0),
  fULSElecPt(0),
  fLSElecPt(0),

  fSparseElectron(0),
  fvalueElectron(0),
  fSparseLSElectron(0),
  fSparseULSElectron(0),
  fvaluePHElectron(0),
  fSparseMulti(0),
  fvalueMulti(0)

{
  // constructor
  fvalueElectron = new Double_t[9];
  fvaluePHElectron = new Double_t[5];
  fvalueMulti = new Double_t[11];
  for(Int_t i=0; i<2; i++) fMultEstimatorAvg[i]=0;
  DefineInput(0, TChain::Class());   
  DefineOutput(1, TList::Class());
 // DefineOutput(2, TList::Class());    
}


//_____________________________________________________________________________
AliAnalysisTaskHFEMultiplicity::~AliAnalysisTaskHFEMultiplicity()
{
  // destructor
  if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    delete fSparseElectron;
    delete fSparseLSElectron;
    delete fSparseULSElectron;
    delete []fvalueElectron;
    delete []fvaluePHElectron;
    delete fSparseMulti;
    delete []fvalueMulti;
    delete fTracks_tender;
    delete fCaloClusters_tender;
  }

  for(Int_t i=0; i<2; i++) {
    if (fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i];
  }
  delete fListProfiles;
  delete gRandom;
}
//_____________________________________________________________________________
void AliAnalysisTaskHFEMultiplicity::Init()
{
  // Initialization
  	
	
  if(fDebug > 1) printf("AliAnalysisTaskHFEMultiplicity::Init() \n");

  
  fListProfiles = new TList();
  fListProfiles->SetOwner();
  TString period[1];
  Int_t nProfiles=1;
  period[0]="LHC16s";
  period[1]="LHC16r"; 
   
    
  for(Int_t i=0; i<nProfiles; i++){
    if(fMultEstimatorAvg[i]){
      TProfile* hprof=new TProfile(*fMultEstimatorAvg[i]);
      hprof->SetName(Form("ProfileTrkVsZvtx%s\n",period[i].Data()));
      fListProfiles->Add(hprof);
    }
  }
    
  //  PostData(2,fListProfiles);
    
   
  return;
}

//_____________________________________________________________________________

void AliAnalysisTaskHFEMultiplicity::UserCreateOutputObjects()
{
 
  printf("\nseed = %u\n\n", gRandom->GetSeed()); 
 
  fOutputList = new TList();          
  fOutputList->SetOwner(kTRUE);       
 
  
  fNevents 		= new TH1F ("fNevents","Number of events",4,-0.5,3.5);
  fNevents->GetYaxis()->SetTitle("counts");
  fNevents->GetXaxis()->SetBinLabel(1,"All");
  fNevents->GetXaxis()->SetBinLabel(2,"With >2 Trks");
  fNevents->GetXaxis()->SetBinLabel(3,"Vtx_{z}<10cm");
  fNevents->GetXaxis()->SetBinLabel(4,"Vtx_{z}<10cm with Trigger");
	
  
  
  fClusPhi    		= new TH1F("fClusPhi", "Cluster Phi distribution; #phi ; counts",100,0.,7);
  fClusEta    		= new TH1F("fClusEta", "Cluster Eta distribution; #eta ; counts",50,-2,2);
  fClusEtaPhi		= new TH2F( "fClusEtaPhi","Cluster Eta Phi distribution; #eta ; #phi",50,-2,2,100,0.,7);
  fClusT     		= new TH1F( "fClusT","Cluster time distribution ; Time(ns) ; counts",500,-1000,1000);
  fNCells   		= new TH1F("fNCells","ncells distribution ; cell counts ; cluster counts", 50,-10,40);
  fClusE   		= new TH1F("fClusE","Cluster Energy ; Energy(GeV); counts",200,0,100);
  fCellE    		= new TH1F("EnergyCell","Cell Energy ; Energy(GeV) ; counts",200,0,100);
  fCellT   		= new TH1F("fCellT","cell time distribution ; Time(ns) ;counts",500,-1000,1000);
  fVtxZ 		= new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
  fVtxY 		= new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",1000,-50,50);
  fVtxX 		= new TH1F("fVtxX","X vertex position;Vtx_{x};counts",1000,-50,50);
  fTPCdEdx 		= new TH2F("fTPCdEdx","All Track dE/dx distribution;p (GeV/c);dE/dx",200,0,20,500,0,160);
  fTPCnsigma 		= new TH2F("fTPCnsig","All Track TPC Nsigma distribution;p (GeV/c);#sigma_{TPC-dE/dx}",1000,0,50,200,-10,10);
  fTrkPt 		= new TH1F("fTrkPt","p_{T} distribution of all tracks;p_{T} (GeV/c);counts",1000,0,100);
  fTrketa 		= new TH1F("fTrketa","All Track #eta distribution;#eta;counts",100,-1.5,1.5);
  fTrkphi 		= new TH1F("fTrkphi","All Track #phi distribution;#phi;counts",100,0,2*3.141);
  fTrkMatchTrkPt	= new TH1F("fTrkMatchTrkPt","p_{T} distribution of tracks with cluster;p_{T} (GeV/c);counts",1000,0,100);
  fTrkMatchTrketa	= new TH1F("fTrkMatchTrketa","#eta distribution of tracks matched to Cluster;#eta;counts",100,-1.5,1.5);
  fTrkMatchTrkphi	= new TH1F("fTrkMatchTrkphi","#phi distribution of tracks matched to Cluster;#phi;counts",100,0,2*3.141);
  fTrkMatchClusetaphi	= new TH2F( "fTrkMatchClusetaphi","#eta#phi distribution of Clusters matched to tracks;#eta;#phi",50,-2,2,100,0.,7);
  fEMCTrkMatchcluster   = new TH2F("fEMCTrkMatchcluster","Distance of EMCAL cluster to its closest track Method 1",100,-0.3,0.3,100,-0.3,0.3);
  fInvmassLS 		= new TH1F("fInvmassLS", "Inv mass of LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
  fInvmassULS 		= new TH1F("fInvmassULS", "Inv mass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
  fInvmassLSPt 		= new TH2F("fInvmassLSPt", "Inv mass of LS (e,e) vs pT; p_{T}(GeV/c); mass(GeV/c^2); counts;", 500,0,50,1000,0,1.0);
  fInvmassULSPt		= new TH2F("fInvmassULSPt", "Inv mass of ULS (e,e) vs pT; p_{T}(GeV/c); mass(GeV/c^2); counts;", 500,0,50,1000,0,1.0);

  fULSElecPt  		= new TH1F("fULSElecPt","p_{T} distribution of ULS electrons;p_{T} (GeV/c);counts",500,0,50);
  fLSElecPt 		= new TH1F("fLSElecPt","p_{T} distribution of LS electrons;p_{T} (GeV/c);counts",500,0,50);
  


  

  
  Int_t bins[9]		=      	{280, 160, 100, 100, 100, 10, 10, 200,1000};
  Double_t xmin[9]	=	{  2,  -8,   0,   0,   0, 0, 0, 0 ,-50};
  Double_t xmax[9]	=	{  30,   8,   2,   2,  2, 100, 100, 100, 50};

  fSparseElectron 	= new THnSparseD ("Electron","Electron;pT;nSigma;E/P;m02;m20;V0M;SPDTracklets;Cluster Energy;zvtx;",9 ,bins,xmin,xmax);
 

  Int_t binsls[5]	=      	{280, 160, 10, 10, 200};
  Double_t xminls[5]	=	{  2, -8, 0, 0, 0};
  Double_t xmaxls[5]	=	{  30,  8, 100, 100, 100};

  fSparseLSElectron 	= new THnSparseD ("LSElectron","LSElectron;pT;nSigma;V0M;SPDTracklets;Cluster Energy;",5 ,binsls,xminls,xmaxls);
  fSparseULSElectron 	= new THnSparseD ("ULSElectron","ULSElectron;pT;nSigma;V0M;SPDTracklets;Cluster Energy;",5 ,binsls,xminls,xmaxls);
  
  Int_t binsm[11]	=      	{ 200, 200, 200, 200,1000,1000,1000,1000,1000,1000,1000};
  Double_t xminm[11]	=	{     0, 0, 0, 0,-50,0,0,0,0,0,0};
  Double_t xmaxm[11]	=	{   100, 100, 100, 100,50,2000,1000,2000,1000,1000,2000};
  fSparseMulti 		= new THnSparseD ("Multiplicity","Multiplicity;V0M_percertile;V0A_percertile;V0C_percertile;SPDTracklets_percertile;zvtx;V0M_class;SPDTracklets_class;V0M_data;SPDTracklets_data;Corrected_SPDTracklets;Corrected_V0M",11,binsm,xminm,xmaxm);
    
    

    

 
  fTrkMatchTrkPt->Sumw2();
  fSparseElectron->Sumw2();
  fSparseLSElectron->Sumw2();
  fSparseULSElectron->Sumw2();
  fSparseMulti->Sumw2();

 
  fOutputList->Add(fNevents);
  fOutputList->Add(fClusPhi);
  fOutputList->Add(fClusEta);
  fOutputList->Add(fClusEtaPhi);
  fOutputList->Add(fClusT);
  fOutputList->Add(fNCells);
  fOutputList->Add(fClusE);
  fOutputList->Add(fCellE);
  fOutputList->Add(fCellT);
  fOutputList->Add(fVtxZ);
  fOutputList->Add(fVtxY);
  fOutputList->Add(fVtxX);
  fOutputList->Add(fTPCdEdx);
  fOutputList->Add(fTPCnsigma);
  fOutputList->Add(fTrkPt);
  fOutputList->Add(fTrketa);
  fOutputList->Add(fTrkphi);
  fOutputList->Add(fTrkMatchTrkPt);
  fOutputList->Add(fTrkMatchTrketa);
  fOutputList->Add(fTrkMatchTrkphi);
  fOutputList->Add(fTrkMatchClusetaphi);
  fOutputList->Add(fEMCTrkMatchcluster);

  fOutputList->Add(fInvmassLS);
  fOutputList->Add(fInvmassULS);
  fOutputList->Add(fInvmassLSPt);
  fOutputList->Add(fInvmassULSPt);
  fOutputList->Add(fULSElecPt);
  fOutputList->Add(fLSElecPt);

  fOutputList->Add(fSparseElectron);
  fOutputList->Add(fSparseLSElectron);
  fOutputList->Add(fSparseULSElectron);
  fOutputList->Add(fSparseMulti);
             
  PostData(1,fOutputList);  
 // PostData(2,fListProfiles);         

}
//_____________________________________________________________________________
void AliAnalysisTaskHFEMultiplicity::UserExec(Option_t *)
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
		
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    
  if(!fAOD) return;

  if(fReadMC){fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));} //MC information

  if(!PassEventSelect(fAOD)) return;
			

  if(fUseTender){
    fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTenderTrackName));
    fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTenderClusterName));
  }
	
 
  //--------------Multiplicity from Multselection class-------------------

  Float_t lPercentiles[4];
  Int_t fMultV0A = -999, fMultV0C = -999, fMultV0M = -999, fMultSPDTracklets = -999;
  Double_t fMultSPDTracklets_Percentile = -999;
  TString lNames[4] = {"V0M", "V0A", "V0C", "SPDTracklets"};
  for(Int_t iEst=0; iEst<4; iEst++) lPercentiles[iEst] = 200;	
  AliMultSelection *MultSelection = 0x0;
  MultSelection= (AliMultSelection*)fAOD->FindListObject("MultSelection");
  if( MultSelection ){

    fMultV0A =  MultSelection->GetEstimator("V0A")->GetValue();
    fMultV0C =  MultSelection->GetEstimator("V0C")->GetValue();
    fMultV0M =  MultSelection->GetEstimator("V0M")->GetValue();
    fMultSPDTracklets =  MultSelection->GetEstimator("SPDTracklets")->GetValue();
    for(Int_t iEst=0; iEst<4; iEst++)
      lPercentiles[iEst] = MultSelection->GetMultiplicityPercentile(lNames[iEst].Data());


  }
  else{
    AliInfo("Didn't find MultSelection!"); 
  }

  //---------------------multiplicity ------------------------------
  Double_t Zvertex1 = -100;						    
  const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
  Zvertex1 =pVtx->GetZ();
 
  if (fRejectPUFromSPD && fAOD->IsPileupFromSPDInMultBins()) return; // pile-up cut
  //--------------------vertex selection cuts-----------------------
  AliAODVertex* vtxSPD = fAOD->GetPrimaryVertexSPD();
	
  if (!vtxSPD || vtxSPD->GetNContributors() < 1) return;
  Double_t cov[6]={0};
  vtxSPD->GetCovarianceMatrix(cov);
  if (TMath::Sqrt(cov[5]) > 0.25) return;

  //----------V0M Multiplicity------------------
  AliAODVZERO *vzeroAOD = dynamic_cast<AliAODVZERO *>( dynamic_cast<AliAODEvent *>(fAOD)->GetVZEROData());
  Int_t V0AMult = static_cast<Int_t>(vzeroAOD->GetMTotV0A());
  Int_t V0CMult = static_cast<Int_t>(vzeroAOD->GetMTotV0C());
  Int_t V0Mult=V0AMult+V0CMult;
  
  //------------SPDTracklets--------------------
  Int_t nTracklets = 0;
  Int_t nAcc = 0;
  Double_t etaRange = 1.0;

  AliAODTracklets *tracklets = static_cast<const AliAODEvent*>(fAOD)->GetTracklets();
  nTracklets = tracklets->GetNumberOfTracklets();
  for (Int_t nn = 0; nn < nTracklets; nn++) {
    Double_t theta = tracklets->GetTheta(nn);
    Double_t eta = -TMath::Log(TMath::Tan(theta/2.0));
    if (TMath::Abs(eta) < etaRange) nAcc++;
  }
  // Data driven multiplicity z-vertex correction
  

  Int_t countMult = nAcc;

  //SPDTracklets correction
  Double_t correctednAcc   = nAcc;
  Double_t countCorr       = countMult;
  TProfile* estimatorAvg = GetEstimatorHistogram(fAOD);
  if(estimatorAvg){
    correctednAcc=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,nAcc,Zvertex1,fRefMult));
    countCorr=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,countMult,Zvertex1,fRefMult));
  } 

  cout<< " corrected tracklets spd = "<< correctednAcc << " incorrected tracklets spd = "<<nAcc<< endl;
  //V0M Correction
  Int_t vzeroMultACorr=V0AMult, vzeroMultCCorr=V0CMult, vzeroMultCorr=V0Mult;
  vzeroMultACorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(V0AMult,Zvertex1));
  vzeroMultCCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0C(V0CMult,Zvertex1));
  vzeroMultCorr = vzeroMultACorr + vzeroMultCCorr;
  cout<< " corrected tracklets V0m = "<< vzeroMultCorr << " incorrected tracklets V0m= "<<V0Mult<< endl;

  fvalueMulti[0] = lPercentiles[0];
  fvalueMulti[1] = lPercentiles[1];
  fvalueMulti[2] = lPercentiles[2];
  fvalueMulti[3] = lPercentiles[3];
  fvalueMulti[4] = Zvertex1;
  fvalueMulti[5] = fMultV0M;
  fvalueMulti[6] = fMultSPDTracklets;
  fvalueMulti[7] = V0Mult;
  fvalueMulti[8] = nAcc;
  fvalueMulti[9] = correctednAcc;
  fvalueMulti[10] = vzeroMultCorr;
  

  fSparseMulti->Fill(fvalueMulti);    



  //-------------------selecting trigger for calorimeter( EMCAL + DCAL )
  TString firedTrigger;
  TString TriggerEG1("EG1");
  TString TriggerEG2("EG2");
  TString TriggerDG1("DG1");
  TString TriggerDG2("DG2");
    
  if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();
    
  Bool_t EG1tr = kFALSE;
  Bool_t EG2tr = kFALSE;
  if(firedTrigger.Contains(TriggerEG1))EG1tr = kTRUE;
  if(firedTrigger.Contains(TriggerEG2))EG2tr = kTRUE;
    
  if(fEMCEG1){if(!firedTrigger.Contains(TriggerEG1))return;}
  if(fEMCEG2){if(!firedTrigger.Contains(TriggerEG2))return;}
  if(fDCalDG1){if(!firedTrigger.Contains(TriggerDG1))return;}
  if(fDCalDG2){if(!firedTrigger.Contains(TriggerDG2))return;}

  fNevents->Fill(3);
        
 
  fpidResponse = fInputHandler->GetPIDResponse();
  if(!fpidResponse) return;

  


  //-----------cluster information---------------------------------------------------------------------------		
 
  ClusterInfo();
	  

  //--------------------Track information------------------------	

 
  Int_t ntracks = -999;
  if(!fUseTender)ntracks = fAOD->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries(); 
		
		
  for (Int_t iTracks=0; iTracks< ntracks; iTracks++) {
    AliAODTrack* track = 0x0;
    if(!fUseTender) track = (AliAODTrack*)fAOD->GetTrack(iTracks);
    if(fUseTender) track =  dynamic_cast<AliAODTrack*>(fTracks_tender->At(iTracks));
    if(!track) continue;
    if(!Passtrackcuts(track)) continue;
   
			
    Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999,dEdx = -999, nsigma = -999.0;
		
    nsigma=fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron); 
 
    dEdx = track->GetTPCsignal();
    
    TrkPhi = track->Phi();
    fTrkphi->Fill(TrkPhi);
    TrkPt = track->Pt();
    fTrkPt->Fill(TrkPt);
			
    TrkEta = track->Eta();
    fTrketa->Fill(TrkEta);
    TrkP = track->P();
				
    fTPCdEdx->Fill(TrkP,dEdx);
    fTPCnsigma->Fill(TrkP,nsigma);

    if(TrkPt<3.0) continue ;

    Int_t EMCalIndex = -1;
    Double_t emcphi = -999, emceta = -999;
    Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;
    EMCalIndex = track -> GetEMCALcluster();
    AliAODCaloCluster *clustMatch=0x0;
    
    if(!fUseTender) if(EMCalIndex >= 0) clustMatch = (AliAODCaloCluster*)fAOD->GetCaloCluster(EMCalIndex) ; 
    if(fUseTender) if(EMCalIndex >= 0)clustMatch = dynamic_cast<AliAODCaloCluster*>(fCaloClusters_tender->At(EMCalIndex));
    if(clustMatch && clustMatch->IsEMCAL())
      {
	Double_t fPhiDiff = -999, fEtaDiff = -999;
        GetTrkClsEtaPhiDiff(track, clustMatch, fPhiDiff, fEtaDiff);
	fEMCTrkMatchcluster->Fill(fPhiDiff,fEtaDiff);
	if(TMath::Abs(fPhiDiff) > 0.05 || TMath::Abs(fEtaDiff)> 0.05) continue;
	
	Float_t EMCalpos[3];
	clustMatch -> GetPosition(EMCalpos);
	TVector3 clustpos(EMCalpos[0],EMCalpos[1],EMCalpos[2]);
	emcphi = clustpos.Phi();
	if(emcphi < 0) emcphi = emcphi + (2*TMath::Pi());
	emceta = clustpos.Eta();
              
	if(emcphi > 1.39 && emcphi < 3.265) fClsTypeEMC = kTRUE; 
	if(emcphi > 4.53 && emcphi < 5.708) fClsTypeDCAL = kTRUE;  
	
	if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
	  if(!fClsTypeEMC) continue; //selecting only EMCAL clusters

	if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
	  if(!fClsTypeDCAL) continue; //selecting only DCAL clusters
            
	fTrkMatchTrkPt->Fill(TrkPt);
        fTrkMatchTrketa->Fill(TrkEta);
	fTrkMatchTrkphi->Fill(TrkPhi);
        fTrkMatchClusetaphi->Fill(emceta,emcphi);

	Double_t Etrkmatch = -999.0, Eoptrk = -999.0 , M02trkmatch = -999.0, M20trkmatch = -999.0;
      
	Etrkmatch = clustMatch->E();
	Eoptrk = Etrkmatch/TrkP ;
	M02trkmatch = clustMatch->GetM02();
	M20trkmatch = clustMatch->GetM20();
 

     
        fvalueElectron[0] = TrkPt;
        fvalueElectron[1] = nsigma;
	fvalueElectron[2] = Eoptrk;
	fvalueElectron[3] = M02trkmatch;
	fvalueElectron[4] = M20trkmatch;
	fvalueElectron[5] = vzeroMultCorr; //V0M, Multiplicity information
	fvalueElectron[6] = correctednAcc; //SPD Tracklets
	fvalueElectron[7] = Etrkmatch;  //cluster energy after matching
	fvalueElectron[8] = Zvertex1;					
	fSparseElectron->Fill(fvalueElectron);   //Electron information sparse         
	
	if(nsigma < -1.  || nsigma > 3.) continue;
	if(M20trkmatch < 0.02 || M20trkmatch> 0.35) continue;
	if(Eoptrk < 0.8 || Eoptrk > 1.2) continue;

	fvaluePHElectron[0] = TrkPt;
	fvaluePHElectron[1] = nsigma;
	fvaluePHElectron[2] = vzeroMultCorr; //V0M, Multiplicity information
	fvaluePHElectron[3] = correctednAcc; //SPD Tracklets
	fvaluePHElectron[4] = Etrkmatch;  //cluster energy after matching
	
	Bool_t fFlagPhotonicElec = kFALSE, fFlagElecLS=kFALSE;
	SelectNonHFElectron(iTracks,track,fFlagPhotonicElec,fFlagElecLS);					    

      }
			
  }
	
  // cout<<"cluster to track match= "<< actrk<<"track to cluster match = " <<btrkc<<endl;	
  PostData(1, fOutputList);   
		
}


//______________________________________________________________________
Bool_t AliAnalysisTaskHFEMultiplicity::Passtrackcuts(AliAODTrack *atrack)
{ 
  Double_t d0z0[2]={-999,-999}, cov[3];
  Double_t DCAxyCut = 2.4, DCAzCut = 3.2;
  Double_t dEdx =-999;
  Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999;

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();

  //kink daughters
  Int_t numberofvertices = 100;
  numberofvertices = fAOD->GetNumberOfVertices();
  Double_t listofmotherkink[numberofvertices];
  Int_t numberofmotherkink = 0;
  for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
    AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
    if(!aodvertex) continue;
    if(aodvertex->GetType()==AliAODVertex::kKink) {
      AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
      if(!mother) continue;
      Int_t idmother = mother->GetID();
      listofmotherkink[numberofmotherkink] = idmother;
      numberofmotherkink++;
    }
  }
  if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE; //minimum cuts- filter bit 4
  //reject kink
  Bool_t kinkmotherpass = kTRUE;
  for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
    if(atrack->GetID() == listofmotherkink[kinkmother]) {
      kinkmotherpass = kFALSE;
      continue;
    }
  }
  if(!kinkmotherpass) return kFALSE;

  //other cuts
  if(atrack->GetTPCNcls() < 80.) return kFALSE;
  if(atrack->GetITSNcls() < 3.) return kFALSE;
  if((!(atrack->GetStatus()&AliAODTrack::kITSrefit)|| (!(atrack->GetStatus()&AliAODTrack::kTPCrefit)))) return kFALSE;
  if(!(atrack->HasPointOnITSLayer(0) || atrack->HasPointOnITSLayer(1))) return kFALSE;

  if(atrack->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov))
    if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) return kFALSE;



}

//______________________________________________________________________
Bool_t AliAnalysisTaskHFEMultiplicity::PassEventSelect(AliAODEvent *fAOD)
{
  Int_t ntracks = -999;
  fNevents->Fill(0);
  Double_t Zvertex=-100, Xvertex=-100, Yvertex=-100;

  const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
  
  Double_t NContV = pVtx->GetNContributors();

  if(NContV<2) return kFALSE;
  fNevents->Fill(1);
 
  Zvertex =pVtx->GetZ();
  Yvertex =pVtx->GetY();
  Xvertex =pVtx->GetX();
  
  if(TMath::Abs(Zvertex)>10.0) return kFALSE;
  fNevents->Fill(2);
  
  
  fVtxZ->Fill(Zvertex);
  fVtxX->Fill(Xvertex);
  fVtxY->Fill(Yvertex);

  return kTRUE;
  
} 

//______________________________________________________________________
void AliAnalysisTaskHFEMultiplicity::ClusterInfo() 
{
  Double_t cluphi = -999.0;
  Double_t clueta =-999.0 ;
  Int_t ncells = -999.0;
  Float_t energy = -999.0;
  Float_t clut = -999.0;
  Double_t  energycell = -999.0;
  Double_t CellId =0;
  Int_t Nclust = -999;               

 
  if(!fUseTender) Nclust = fAOD->GetNumberOfCaloClusters(); 
  if(fUseTender) Nclust = fCaloClusters_tender->GetEntries();

  Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;
  for ( Int_t index = 0; index < Nclust ; index++ ) {
    AliAODCaloCluster * clu =0x0;
    if(!fUseTender) clu  = (AliAODCaloCluster*)fAOD->GetCaloCluster(index) ; 
    if(fUseTender) clu = dynamic_cast<AliAODCaloCluster*>(fCaloClusters_tender->At(index));
    if(!clu) continue;
	  
    fClsTypeEMC = kFALSE; fClsTypeDCAL = kFALSE;
    if (clu->IsEMCAL()){
	
      AliAODCaloCells &cells = *(fAOD->GetEMCALCells());
	
      Float_t  x[3]; // cluster pos
      Double_t V[3];
      Double_t clup,Etrans,b;
      fAOD->GetVertex()->GetXYZ(V);
      TLorentzVector p;
      clu->GetMomentum(p,V);
      Etrans = p.Et();
      b = clu->GetM02();


      clu->GetPosition(x);
      TVector3 clustposi(x[0],x[1],x[2]);
   

      cluphi = clustposi.Phi();
      clueta = clustposi.Eta();
      if(cluphi < 0) cluphi = cluphi+(2*TMath::Pi());
      if(cluphi > 1.39 && cluphi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
      if(cluphi > 4.53 && cluphi < 5.708) fClsTypeDCAL = kTRUE;//DCAL  : 260 < phi < 327

      if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
	if(!fClsTypeEMC) continue; //selecting only EMCAL clusters

      if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
	if(!fClsTypeDCAL) continue; //selecting only DCAL clusters

      clut = clu->GetTOF()*1e9 ;
      energy = clu->E();
      ncells= clu->GetNCells();
		

      //-----------Cell information
 
      /*    UShort_t * C    = clu->GetCellsAbsId() ;
	    Double_t *fraction = clu->GetCellsAmplitudeFraction() ;
	    for(Int_t i = 0; i < ncells ; i++){
	    Int_t absId       = C[i]; 
	    Double_t ampFract = fraction[i];
	    Double_t Ecell    = cells.GetCellAmplitude(absId);
	    Double_t Tcell 	  = cells.GetCellTime(absId)*1e9;
	    fCellE->Fill(Ecell); 
	    fCellT->Fill(Tcell);
	
	    }*/

	
      fClusPhi->Fill(cluphi);
      fClusEtaPhi->Fill(clueta,cluphi);
      fClusEta->Fill(clueta);
      fNCells->Fill(ncells);
      fClusE->Fill(energy);
      fClusT->Fill(clut);

    
    }	    


		
  }
	
}

//______________________________________________________________________
void AliAnalysisTaskHFEMultiplicity::SelectNonHFElectron(Int_t itrack, AliAODTrack *track, Bool_t &fFlagPhotonicElec, Bool_t &fFlagElecLS)
{
  //Photonic electron selection

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
  Double_t d0z0[2]={-999,-999}, cov[3];
  Double_t DCAxyCut = 2.4, DCAzCut = 3.2;

  Bool_t flagPhotonicElec = kFALSE, flagLSElec = kFALSE;
  Double_t ptAsso=-999., nsigmaAsso=-999.0;
  Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;


  Int_t ntracks = -999;
  if(!fUseTender)ntracks = fAOD->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries();

  for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
    if(jTracks==itrack) continue;


    AliAODTrack* atrackAsso = 0x0;
    if(!fUseTender) atrackAsso = (AliAODTrack*)fAOD->GetTrack(jTracks);
    if(fUseTender) atrackAsso =  dynamic_cast<AliAODTrack*>(fTracks_tender->At(jTracks));
    if(!atrackAsso) continue;

    if(!atrackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
    if(atrackAsso->GetTPCNcls() < 80.) continue;
  
    nsigmaAsso = fpidResponse->NumberOfSigmasTPC(atrackAsso, AliPID::kElectron);
    ptAsso = atrackAsso->Pt();
    Int_t chargeAsso = atrackAsso->Charge();
    Int_t charge = track->Charge();

    if(ptAsso < 0.1) continue; //partner electron pt 0.3
    if(atrackAsso->Eta()<-0.9 || atrackAsso->Eta()>0.9) continue;
    if(nsigmaAsso < -3 || nsigmaAsso > 3) continue;

    if(atrackAsso->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov))
      if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;

    Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
    if(charge>0) fPDGe1 = -11;
    if(chargeAsso>0) fPDGe2 = -11;

    fFlagLS=kFALSE; fFlagULS=kFALSE;
    if(charge == chargeAsso) fFlagLS = kTRUE;
    if(charge != chargeAsso) fFlagULS = kTRUE;

    AliKFParticle::SetField(fAOD->GetMagneticField());

    AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
    AliKFParticle ge2 = AliKFParticle(*atrackAsso, fPDGe2);
    AliKFParticle recg(ge1, ge2);

    if(recg.GetNDF()<1) continue;
    Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
    if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;

    Double_t mass=-999., width = -999.;
    Int_t MassCorrect;
    MassCorrect = recg.GetMass(mass,width);

    if(fFlagLS && track->Pt()>1) fInvmassLS->Fill(mass);
    if(fFlagULS && track->Pt()>1) fInvmassULS->Fill(mass);
    if(fFlagLS) fInvmassLSPt->Fill(track->Pt(),mass);
    if(fFlagULS) fInvmassULSPt->Fill(track->Pt(),mass);

    if(fFlagLS && mass<0.14) {
      fLSElecPt->Fill(track->Pt());
      fSparseLSElectron->Fill(fvaluePHElectron);
    }
    if(fFlagULS && mass<0.14){
      fULSElecPt->Fill(track->Pt());
      fSparseULSElectron->Fill(fvaluePHElectron); 
    }

    if(mass<0.14 && fFlagULS && !flagPhotonicElec){
      flagPhotonicElec = kTRUE;
    }
    if(mass<0.14 && fFlagULS && !flagPhotonicElec){
      flagLSElec = kTRUE;
    }
  }
  fFlagPhotonicElec = flagPhotonicElec;
  fFlagElecLS = flagLSElec;

}
//____________________________________________________________________________
TProfile* AliAnalysisTaskHFEMultiplicity::GetEstimatorHistogram(const AliAODEvent* fAOD)
{
    
  Int_t runNo  = fAOD->GetRunNumber();
  Int_t period = -1; 
   
        
  if (runNo>266437 && runNo<267110) period = 0;
  if (runNo>265744 && runNo<266318) period = 1;
  if (period < 0 || period > 1) return 0;
    
  cout<<"period ="<<period<<endl;
    
  return fMultEstimatorAvg[period];
}
//____________________________________________________________________________
void AliAnalysisTaskHFEMultiplicity::GetTrkClsEtaPhiDiff(AliAODTrack *t, AliAODCaloCluster *v, Double_t &phidiff, Double_t &etadiff)
{
  // Calculate phi and eta difference between a track and a cluster. The position of the track is obtained on the EMCAL surface

  phidiff = 999;
  etadiff = 999;

  if (!t||!v) return;

  Double_t veta = t->GetTrackEtaOnEMCal();
  Double_t vphi = t->GetTrackPhiOnEMCal();

  Float_t pos[3] = {0};
  v->GetPosition(pos);
  TVector3 cpos(pos);
  Double_t ceta     = cpos.Eta();
  Double_t cphi     = cpos.Phi();
  etadiff=veta-ceta;
  phidiff=TVector2::Phi_mpi_pi(vphi-cphi);
}
//______________________________________________________________________________
void AliAnalysisTaskHFEMultiplicity::Terminate(Option_t *)
{
  
}
//_________________________________________________________________________
