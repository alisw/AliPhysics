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
#include "AliESDtrackCuts.h"
#include "AliCentralitySelectionTask.h"
#include "AliMultSelection.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliESDCaloTrigger.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"
#include "AliCentrality.h"
#include "AliMagF.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "TProfile.h"
#include "AliESDVZERO.h"
#include "AliAODVZERO.h"
#include "TVector3.h"
#include "TRandom2.h"

class AliAnalysisTaskHFEMultiplicity;    
using namespace std;            
ClassImp(AliAnalysisTaskHFEMultiplicity)

AliAnalysisTaskHFEMultiplicity::AliAnalysisTaskHFEMultiplicity() : AliAnalysisTaskSE(), 
  fAOD(0),
  fNevents(0),
  fTenderClusterName("caloClusters"),
  fTenderTrackName("tracks"),
  fOutputList(0),
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
  fClusMatchTrkPt(0),
  fClusMatchTrketa(0),
  fClusMatchTrkphi(0),
  fMatchClusphi(0), 
  fMatchCluseta(0),
  fMatchClusetaphi(0x0),
  fMatchClusEnergy(0),							 
  fCentrOrMult(-1),
  fHistCent(0),								 
  fEMCTrkMatch(0x0),
  fSparseElectron(0),
  fvalueElectron(0),
  fSparseMulti(0),
  fvalueMulti(0)
  
{   fvalueElectron = new Double_t[8];
  fvalueMulti = new Double_t[4];
  
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}

//_____________________________________________________________________________
AliAnalysisTaskHFEMultiplicity::AliAnalysisTaskHFEMultiplicity(const char* name) : AliAnalysisTaskSE(name),

  fAOD(0),
  fNevents(0),
  fTenderClusterName("caloClusters"),
  fTenderTrackName("tracks"),
  fOutputList(0),
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
  fClusMatchTrkPt(0),
  fClusMatchTrketa(0),
  fClusMatchTrkphi(0),
  fMatchClusphi(0), 
  fMatchCluseta(0),
  fMatchClusetaphi(0x0),
  fMatchClusEnergy(0),															                                                          
  fCentrOrMult(-1),
  fHistCent(0),								 
  fEMCTrkMatch(0x0),
  fSparseElectron(0),
  fvalueElectron(0),
  fSparseMulti(0),
  fvalueMulti(0)								 
{
  // constructor
  fvalueElectron = new Double_t[8];
  fvalueMulti = new Double_t[4];
  DefineInput(0, TChain::Class());   
  DefineOutput(1, TList::Class());    
}


//_____________________________________________________________________________
AliAnalysisTaskHFEMultiplicity::~AliAnalysisTaskHFEMultiplicity()
{
  // destructor
  if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    delete fSparseElectron;
    delete []fvalueElectron;
    delete fSparseMulti;
    delete []fvalueMulti;
    delete fTracks_tender;
    delete fCaloClusters_tender;
  }

}
//_____________________________________________________________________________

void AliAnalysisTaskHFEMultiplicity::UserCreateOutputObjects()
{
    
  fOutputList = new TList();          
  fOutputList->SetOwner(kTRUE);       
  
  // example of a histogram
  
  fNevents 		= new TH1F ("fNevents","Number of events",4,-0.5,3.5);
  fNevents->GetYaxis()->SetTitle("counts");
  fNevents->GetXaxis()->SetBinLabel(1,"All");
  fNevents->GetXaxis()->SetBinLabel(2,"With >2 Trks");
  fNevents->GetXaxis()->SetBinLabel(3,"Vtx_{z}<10cm");
  fNevents->GetXaxis()->SetBinLabel(4,"Vtx_{z}<10cm with Trigger");
	
  
  fHistCent		= new TH1F("fHistCent", "centrality distribution ; centrality(%) ; counts", 100, 0, 100);
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
  fClusMatchTrkPt 	= new TH1F("fClusMatchTrkPt","p_{T} distribution of tracks with cluster;p_{T} (GeV/c);counts",1000,0,100);
  fClusMatchTrketa 	= new TH1F("fClusMatchTrketa","#eta distribution of tracks matched to Cluster;#eta;counts",100,-1.5,1.5);
  fClusMatchTrkphi 	= new TH1F("fClusMatchTrkphi","#phi distribution of tracks matched to Cluster;#phi;counts",100,0,2*3.141);
  fMatchClusphi  	= new TH1F("fMatchClusphi", "#phi distribution of Clusters matched to track;#phi;counts",100,0.,7);
  fMatchCluseta  	= new TH1F("fMatchCluseta", "#eta distribution of Clusters matched to tracks;#eta;counts",50,-2,2);
  fMatchClusetaphi	= new TH2F( "fMatchClusetaphi","#eta#phi distribution of Clusters matched to tracks;#eta;#phi",50,-2,2,100,0.,7);
  fMatchClusEnergy	= new TH1F("fMatchClusEnergy","Cluster Energy after matching to tracks",200,0,100);
  fEMCTrkMatch 		= new TH2F("fEMCTrkMatch","Distance of EMCAL cluster to its closest track Method 1",100,-0.3,0.3,100,-0.3,0.3);
  
  Int_t bins[8]		=      	{280, 160, 100, 100, 100, 10, 10, 200};
  Double_t xmin[8]	=	{  2,  -8,   0,   0,   0, 0, 0, 0 };
  Double_t xmax[8]	=	{  30,   8,   2,   2,  2, 100, 100, 100};
  fSparseElectron 	= new THnSparseD ("Electron","Electron;pT;nSigma;E/P;m02;m20;V0M;SPDTracklets;Cluster Energy;",8 ,bins,xmin,xmax);
  
  Int_t binsm[4]	=      	{ 10, 10, 10, 10};
  Double_t xminm[4]	=	{     0, 0, 0, 0};
  Double_t xmaxm[4]	=	{   100, 100, 100, 100};
  fSparseMulti 		= new THnSparseD ("Multiplicity","Multiplicity;V0M;V0A;V0C;SPDTracklets;",4,binsm,xminm,xmaxm);
    
    

    
  fHistCent->Sumw2();
  fClusMatchTrkPt->Sumw2();
  fSparseElectron->Sumw2();
  fSparseMulti->Sumw2();


  fOutputList->Add(fNevents);
  fOutputList->Add(fHistCent);
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
  fOutputList->Add(fClusMatchTrkPt);
  fOutputList->Add(fClusMatchTrketa);
  fOutputList->Add(fClusMatchTrkphi);
  fOutputList->Add(fMatchClusphi); 
  fOutputList->Add(fMatchCluseta);
  fOutputList->Add(fMatchClusetaphi);
  fOutputList->Add(fMatchClusEnergy);
  fOutputList->Add(fEMCTrkMatch);
  fOutputList->Add(fSparseElectron);
  fOutputList->Add(fSparseMulti);
                            
  PostData(1, fOutputList);           

}
//_____________________________________________________________________________
void AliAnalysisTaskHFEMultiplicity::UserExec(Option_t *)
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
		
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    
  if(!fAOD) return;
  if(!PassEventSelect(fAOD)) return;
			

  if(fUseTender){
    fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTenderTrackName));
    fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTenderClusterName));
  }
	
 
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

 Float_t lPercentiles[4];
  TString lNames[4] = {"V0M", "V0A", "V0C", "SPDTracklets"};
  for(Int_t iEst=0; iEst<4; iEst++) lPercentiles[iEst] = 200;	
  AliMultSelection *MultSelection = 0x0;
  MultSelection= (AliMultSelection*)fAOD->FindListObject("MultSelection");
  if( MultSelection ){
    for(Int_t iEst=0; iEst<4; iEst++)
      lPercentiles[iEst] = MultSelection->GetMultiplicityPercentile(lNames[iEst].Data());
  }
  else{
    AliInfo("Didn't find MultSelection!"); 
  }



  fvalueMulti[0] = lPercentiles[0];
  fvalueMulti[1] = lPercentiles[1];
  fvalueMulti[2] = lPercentiles[2];
  fvalueMulti[3] = lPercentiles[3];

  fSparseMulti->Fill(fvalueMulti);            
						    
  fHistCent->Fill(lPercentiles[0] );
 

  fpidResponse = fInputHandler->GetPIDResponse();
  if(!fpidResponse) return;


  //-----------cluster information---------------------------------------------------------------------------		
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
 
      UShort_t * C    = clu->GetCellsAbsId() ;
      Double_t *fraction = clu->GetCellsAmplitudeFraction() ;
	

      for(Int_t i = 0; i < ncells ; i++){
	Int_t absId       = C[i]; 
	Double_t ampFract = fraction[i];
	Double_t Ecell    = cells.GetCellAmplitude(absId);
	Double_t Tcell 	  = cells.GetCellTime(absId)*1e9;
	  
	
	fCellE->Fill(Ecell); 
	fCellT->Fill(Tcell);
	
      }

	
      fClusPhi->Fill(cluphi);
      fClusEtaPhi->Fill(clueta,cluphi);
      fClusEta->Fill(clueta);
      fNCells->Fill(ncells);
      fClusE->Fill(energy);
      fClusT->Fill(clut);

      //--------------------------------------------cluster track matching--------------------------------------------------------------

	
		    
      AliAODTrack *atrack=0x0;
      if(clu->GetNTracksMatched()>0){

	atrack=dynamic_cast<AliAODTrack*>(clu->GetTrackMatched(0)); //closest track to cluster
							 
	if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue ;
	Double_t mTrkPhi=-999, mTrkPt=-999, mTrkEta=-999, mTrkP = -999,mdEdx = -999, mnsigma = -999.0;
		
	mnsigma=fpidResponse->NumberOfSigmasTPC(atrack, AliPID::kElectron);
						
	//--------------matched track properties						 
	mTrkPt=atrack->Pt();
	mTrkEta=atrack->Eta();
	mTrkPhi=atrack->Phi();
	mTrkP = atrack->P();
							  
	if(mTrkPt<3.0) continue ;

	fClusMatchTrkPt->Fill(mTrkPt);
	fEMCTrkMatch->Fill(clu->GetTrackDx(),clu->GetTrackDz());
	fClusMatchTrketa->Fill(mTrkEta);
	fClusMatchTrkphi->Fill(mTrkPhi);
	//----------------matched cluster properties
	Double_t clustMatchE = clu->E();
	fMatchClusEnergy->Fill(clustMatchE);
	fMatchCluseta->Fill(clueta);
        fMatchClusphi->Fill(cluphi);
	fMatchClusetaphi->Fill(clueta,cluphi);
                                                
		
	Double_t Ematch = -999.0, Eop = -999.0 , M02match = -999.0, M20match = -999.0;
	Ematch = clu->E();
	Eop = Ematch/mTrkP ;
	M02match = clu->GetM02();
	M20match = clu->GetM20();
//-------------Electron information							   
	fvalueElectron[0] = mTrkPt;
	fvalueElectron[1] = mnsigma;
	fvalueElectron[2] = Eop;
	fvalueElectron[3] = M02match;
	fvalueElectron[4] = M20match;
	fvalueElectron[5] = lPercentiles[0]; //V0M, Multiplicity information
	fvalueElectron[6] = lPercentiles[3]; 
	fvalueElectron[7] = clustMatchE; 
								
	fSparseElectron->Fill(fvalueElectron);   //Electron information sparse         
						    
		   
		   						          								
      }
		           

    }	    


		
  }
		


	  

  //--------------------Track information------------------------	

  const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();	
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
  Int_t ntracks = -999;
  if(!fUseTender)ntracks = fAOD->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries(); 
		
		
  for (Int_t iTracks=0; iTracks< ntracks; iTracks++) {
    AliAODTrack* track = 0x0;
    if(!fUseTender) track = (AliAODTrack*)fAOD->GetTrack(iTracks);
    if(fUseTender) track =  dynamic_cast<AliAODTrack*>(fTracks_tender->At(iTracks));
		  
    if(!track) continue;
		  
		  
		  
    Bool_t kinkmotherpass = kTRUE;
    for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
      if(track->GetID() == listofmotherkink[kinkmother]) {
	kinkmotherpass = kFALSE;
	continue;
      }
    }
		  
    if(!kinkmotherpass) continue;
		  
    //other cuts
    Double_t d0z0[2]={-999,-999}, cov[3];
    Double_t DCAxyCut = 2.4, DCAzCut = 3.2;
    if(!Passtrackcuts(track)) continue;
    if(track->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov))
      if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;
			
 			
    if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue ;
			
    Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999,dEdx = -999, nsigma = -999.0;
		
    nsigma=fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron); 
    cout<<"nsigma value"<<nsigma<<endl;
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
			
 }
		
		
  PostData(1, fOutputList);   
		
}

//____________________________________________________________________________________

//__________________________________________________________________________________

//______________________________________________________________________
Bool_t AliAnalysisTaskHFEMultiplicity::Passtrackcuts(AliAODTrack *atrack)
{ 

  
  if(atrack->GetTPCNcls() < 80) return kFALSE;
  if(atrack->GetITSNcls() < 3) return kFALSE;
  if((!(atrack->GetStatus()&AliAODTrack::kITSrefit)|| (!(atrack->GetStatus()&AliAODTrack::kTPCrefit)))) return kFALSE;
  if(!(atrack->HasPointOnITSLayer(0) || atrack->HasPointOnITSLayer(1))) return kFALSE;

  
  

}
//______________________________________________________________________
Bool_t AliAnalysisTaskHFEMultiplicity::PassEventSelect(AliAODEvent *fAOD)
{
  Int_t ntracks = -999;
  ntracks= fAOD->GetNumberOfTracks();
  cout<<"number of tracks="<<ntracks<<endl;
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

//____________________________________________________________________________

void AliAnalysisTaskHFEMultiplicity::Terminate(Option_t *)
{
  
}
//_____________________________________________________________________________
