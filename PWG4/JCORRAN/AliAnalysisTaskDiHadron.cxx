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
//2- and 3-particle trigger particle correlation analysis
//Author: Jason Glyndwr Ulery, ulery@uni-frankfurt.de
//version: 3.4,  last revised: 2010/08/15

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFormula.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TVector3.h"
#include "TMath.h"


#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "TParticle.h"


//#include "AliHeader.h"
//#include "AliGenEventHeader.h"



#include "AliAnalysisTaskDiHadron.h"


ClassImp(AliAnalysisTaskDiHadron)

//----------------------------------------
AliAnalysisTaskDiHadron::AliAnalysisTaskDiHadron(const char *name):
AliAnalysisTask(name,""), fESD(0), fMC(0), fOutput(0),fMinClustersTPC(0),fMinClusterRatio(0),fMaxTPCchi2(0),fMinClustersITS(0),fEtaCut(0),fTrigEtaCut(0),fNearPhiCut(0),fXECut(0),fMaxDCA(0),fMaxDCAXY(0),fMaxDCAZ(0),fDCA2D(0),fTPCRefit(0),fITSRefit(0),fSPDCut(0),fMinPtAssoc(0),fMaxPtAssoc(0),fVzCut(0),fEfficiencyCorr(0),fDEBUG(0),fnBinPhi(0),fnBinEta(0),fnBinPhiEtaPhi(0),fnBinPhiEtaEta(0),fnBinPhi3(0),fnBinEta3(0),fPi(3.1415926535898),fdPhiMin(0),fdPhiMax(0),fNTPtBins(0),fNMix(0),fNCentBins(0),fNAPtBins(0),fNAPt3Bins(0),fNVertexBins(0),fNXEBins(0),fNIDs(0),fEffFitPt(0),fNFitLowParam(0),fNFitHighParam(0),fMCHistos(0),fFitLow(NULL),fFitHigh(NULL),fFitLowParam(NULL),fFitHighParam(NULL),fPtTrigArray(NULL),fPtAssocArray(NULL),fPtAssoc3Array1(NULL),fPtAssoc3Array2(NULL),fCentArrayMin(NULL),fCentArrayMax(NULL),fXEArray(NULL),fTrigIDArray(NULL),ftPhi(NULL),ftEta(NULL),ftPt(NULL),ftCharge(NULL),ftEff(NULL),ftPtAssoc3(NULL),ftNPtAssoc3(NULL)
  
  {
 
  //IO Slots
  DefineInput(0, TChain::Class());
  DefineOutput(0,TList::Class());
  
  
  for(int c=0;c<fNCentBins;c++){
    for(int v=0;v<fNVertexBins;v++){
      for(int jmc=0;jmc<2;jmc++){
	fMixPointer[c][v][jmc]=-1;
	fMixEnd[c][v][jmc]=-1;
	for(int ievts=0;ievts<fNMix;ievts++){
	  fMPt[ievts][c][v][jmc]=NULL;
	  fMPhi[ievts][c][v][jmc]=NULL;
	  fMEta[ievts][c][v][jmc]=NULL;
	  for(int dd=0;dd<10;dd++)fMPtAssoc3[ievts][c][v][jmc][dd]=NULL;
	  fMNPtAssoc3[ievts][c][v][jmc]=NULL;
	  fMixTrack[ievts][c][v][jmc]=0;
	}
      }
    }
  }
 
  }
//--------------------------------------
void AliAnalysisTaskDiHadron::SetCuts(Int_t MinClustersTPC,  Float_t MinClusterRatio, Float_t MaxTPCchi2, Int_t MinClustersITS, Float_t EtaCut, Float_t TrigEtaCut, Float_t NearPhiCut, Float_t XECut, Float_t MaxDCA, Float_t MaxDCAXY, Float_t MaxDCAZ, Int_t DCA2D, Int_t TPCRefit, Int_t ITSRefit, Int_t SPDCut, Float_t MinPtAssoc, Float_t MaxPtAssoc, Float_t VzCut, Int_t NIDs, const char * TrigIDArray){
//Sets the varibles for track and event cuts
  fMinClustersTPC=MinClustersTPC;
  fMinClusterRatio=MinClusterRatio;
  fMaxTPCchi2=MaxTPCchi2;
  fMinClustersITS=MinClustersITS;
  fEtaCut=EtaCut;
  fTrigEtaCut=TrigEtaCut;
  fNearPhiCut=NearPhiCut;
  fXECut=XECut;
  fMaxDCA=MaxDCA;
  fMaxDCAXY=MaxDCAXY;
  fMaxDCAZ=MaxDCAZ;
  fDCA2D=DCA2D;
  fTPCRefit=TPCRefit;
  fITSRefit=ITSRefit;
  fSPDCut=SPDCut;
  fMinPtAssoc=MinPtAssoc;
  fMaxPtAssoc=MaxPtAssoc;
  fVzCut=VzCut;
  fNIDs=NIDs;
  fTrigIDArray=(char*)TrigIDArray;
}
//--------------------------------------------------------
void AliAnalysisTaskDiHadron::SetOptions(Int_t EfficiencyCorr, Int_t ffDEBUG,  Int_t MCHistos){
//Sets some options
  fEfficiencyCorr=EfficiencyCorr;
  fDEBUG=ffDEBUG;
  fMCHistos=MCHistos;
}
//------------------------------------------------------
void AliAnalysisTaskDiHadron::SetBins(Int_t nBinPhi, Int_t nBinEta, Int_t nBinPhiEtaPhi, Int_t nBinPhiEtaEta, Int_t nBinPhi3, Int_t nBinEta3,Float_t dPhiMin, Float_t dPhiMax, Int_t NTPtBins, Int_t NMixBins, Int_t NCentBins,Int_t NAPtBins, Int_t NAPt3Bins, Int_t NVertexBins, Int_t NXEBins,Float_t *PtTrigArray, Float_t *PtAssocArray,Float_t *PtAssoc3Array1, Float_t *PtAssoc3Array2, Int_t *CentArrayMin, Int_t *CentArrayMax, Float_t *XEArray){
//sets up the histogram binning
  fnBinPhi=nBinPhi;
  fnBinEta=nBinEta;
  fnBinPhiEtaPhi=nBinPhiEtaPhi;
  fnBinPhiEtaEta=nBinPhiEtaEta;
  fnBinPhi3=nBinPhi3;
  fnBinEta3=nBinEta3;
  fdPhiMin=dPhiMin;
  fdPhiMax=dPhiMax;
  fNTPtBins=NTPtBins;
  fNMix=NMixBins;
  fNCentBins=NCentBins;
  fNAPtBins=NAPtBins;
  fNAPt3Bins=NAPt3Bins;
  fNVertexBins=NVertexBins;
  fNXEBins=NXEBins;
  fPtTrigArray=new Float_t [fNTPtBins];
  for(int i=0;i<fNTPtBins;i++)fPtTrigArray[i]=PtTrigArray[i];
  fPtAssocArray=new Float_t [fNAPtBins];
  for(int i=0;i<fNAPtBins;i++)fPtAssocArray[i]=PtAssocArray[i];
  fPtAssoc3Array1=new Float_t [fNAPt3Bins];
  for(int i=0;i<fNAPt3Bins;i++)fPtAssoc3Array1[i]=PtAssoc3Array1[i];
  fPtAssoc3Array2=new Float_t [fNAPt3Bins];
  for(int i=0;i<fNAPt3Bins;i++)fPtAssoc3Array2[i]=PtAssoc3Array2[i];
  fCentArrayMin=new Int_t [fNCentBins];
  for(int i=0;i<NCentBins;i++)fCentArrayMin[i]=CentArrayMin[i];
  fCentArrayMax=new Int_t [fNCentBins];
  for(int i=0;i<NCentBins;i++)fCentArrayMax[i]=CentArrayMax[i];
  fXEArray=new Float_t [fNXEBins];
  for(int i=0;i<fNXEBins;i++)fXEArray[i]=XEArray[i];
 for(int i=0;i<=fNVertexBins;i++)fVertexArray[i]=(2.*i/fNVertexBins-1)*fVzCut;
}
//-------------------------------------------------------
void AliAnalysisTaskDiHadron::SetEfficiencies(Float_t EffFitPt, const TF1 *FitLow, const TF1 *FitHigh, Int_t NFitLowParam, Int_t NFitHighParam, Float_t *FitLowParam, Float_t *FitHighParam){
//Sets up the efficiency corrections
  fEffFitPt=EffFitPt;
  fFitLow=(TF1*)FitLow;
  fFitHigh=(TF1*)FitHigh;
  fNFitLowParam=NFitLowParam;
  fNFitHighParam=NFitHighParam;
  fFitLowParam=new Float_t [fNFitLowParam*fNCentBins];
  for(int i=0;i<fNFitLowParam*fNCentBins;i++)fFitLowParam[i]=FitLowParam[i];
  fFitHighParam=new Float_t [fNFitHighParam*fNCentBins];
  for(int i=0;i<fNFitHighParam*fNCentBins;i++)fFitHighParam[i]=FitHighParam[i];
}

//-----------------------------------------------------------
void AliAnalysisTaskDiHadron::ConnectInputData(Option_t *){
  //Connect to ESD
  if(fDEBUG)Printf("Connecting");
   TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
   if (!tree&&fDEBUG) {Printf("ERROR: Could not read chain from input slot 0");} 
  else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if (!esdH){
      if(fDEBUG)Printf("ERROR: Could not get ESDInputHandler");
    }
    else{ 
      fESD = esdH->GetEvent();
    }
    
    //MC Data handler (so one can calcualte eff)
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    if(mcH)fMC=mcH->MCEvent();
  }
  if(fDEBUG)Printf("Connected");
}
  
//---------------------------------------------------------  
void AliAnalysisTaskDiHadron::CreateOutputObjects(){
  //Creates the histograms and list
  if(fDEBUG)Printf("Output");
  fOutput=new TList();
  fOutput->SetName(GetName());
  char histname[100];
  char histtitle[200];
  int nptbins=fNAPtBins;
  int lptbins=0;
  const char *cmc1[2]={"","_MC"};
  const char *cmc2[2]={""," MC"};
  const char *sign1[3]={"","_LS","_ULS"};
  const char *sign2[3]={""," Like-Sign"," Unlike-Sign"};
  const char *sign31[4]={"","_LS","_ULT","_ULA"};
  const char *sign32[4]={""," Like-Sign"," Trigger-Diff"," Assoc-Diff"};
  Float_t etaEdge=fEtaCut+fTrigEtaCut;
  Float_t phiArray[fnBinPhi+1];
  Float_t etaArray[fnBinEta+1];
  Float_t phiEtaArrayPhi[fnBinPhiEtaPhi+1];
  Float_t phiEtaArrayEta[fnBinPhiEtaEta+1];
  for(int iphi=0;iphi<=fnBinPhi;iphi++){
    phiArray[iphi]=fdPhiMin+iphi*2*fPi/fnBinPhi;
  }
  for(int ieta=0;ieta<=fnBinEta;ieta++){
    etaArray[ieta]=-etaEdge+ieta*2*etaEdge/fnBinEta;
  }
  for(int iphi=0;iphi<=fnBinPhiEtaPhi;iphi++){
    phiEtaArrayPhi[iphi]=fdPhiMin+iphi*2*fPi/fnBinPhiEtaPhi;
  }
  for(int ieta=0;ieta<=fnBinPhiEtaEta;ieta++){
    phiEtaArrayEta[ieta]=-etaEdge+ieta*2*etaEdge/fnBinPhiEtaEta;
  }
  for(int imc=0;imc<=1;imc++){//MC loop
    if(imc==1&&!fMCHistos) continue;
    //Create the histograms
    Int_t buffersize = 256;
    snprintf(histname,buffersize,"fHistMult%s",cmc1[imc]);
    snprintf(histtitle,buffersize,"Multiplicity%s",cmc2[imc]);
    fHistMult[imc]=new TH1F(histname,histtitle,2000,-0.5,1999.5);
    fHistMult[imc]->Sumw2();
    fHistMult[imc]->GetXaxis()->SetTitle("Number of tracks");
    fHistMult[imc]->GetYaxis()->SetTitle("Counts");
    fOutput->Add(fHistMult[imc]);
    
    for(int imult=0;imult<fNCentBins;imult++){//loop for multiplicity bins
      
      //Histograms that are independent of the trigger
      sprintf(histname,"fHistPt_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"P_{T} Distribution of Tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPt[imult][imc]=new TH1F(histname,histtitle,nptbins,fPtAssocArray);
      fHistPt[imult][imc]->Sumw2();
      fHistPt[imult][imc]->GetXaxis()->SetTitle("p_{T}");
      fHistPt[imult][imc]->GetYaxis()->SetTitle("Counts");
      fOutput->Add(fHistPt[imult][imc]);

 //Histograms that are independent of the trigger
      sprintf(histname,"fHistPtEff_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"P_{T} Distribution of Tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPtEff[imult][imc]=new TH1F(histname,histtitle,1000,0,100);
      fHistPtEff[imult][imc]->Sumw2();
      fHistPtEff[imult][imc]->GetXaxis()->SetTitle("p_{T}");
      fHistPtEff[imult][imc]->GetYaxis()->SetTitle("Counts");
      fOutput->Add(fHistPtEff[imult][imc]);
      
      sprintf(histname,"fHistPhi_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"#phi Distribution of Tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPhi[imult][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,nptbins,fPtAssocArray);
      fHistPhi[imult][imc]->Sumw2();
      fHistPhi[imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistPhi[imult][imc]->GetYaxis()->SetTitle("P_{T}");
      fOutput->Add(fHistPhi[imult][imc]);
      
      sprintf(histname,"fHistPhiPt_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"P_{T} weighted #phi Distribution of tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPhiPt[imult][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,nptbins,fPtAssocArray);
      fHistPhiPt[imult][imc]->Sumw2();
      fHistPhiPt[imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistPhiPt[imult][imc]->GetYaxis()->SetTitle("P_{T}");
      fOutput->Add(fHistPhiPt[imult][imc]);
      
      sprintf(histname,"fHistEta_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"#eta Distribution of Tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistEta[imult][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,nptbins,fPtAssocArray);
      fHistEta[imult][imc]->Sumw2();
      fHistEta[imult][imc]->GetXaxis()->SetTitle("#eta");
      fHistEta[imult][imc]->GetYaxis()->SetTitle("P_{T}");
      fOutput->Add(fHistEta[imult][imc]);
      
      sprintf(histname,"fHistEtaPt_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"P_{T} weighted #eta Distribution of tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistEtaPt[imult][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,nptbins,fPtAssocArray);
      fHistEtaPt[imult][imc]->Sumw2();
      fHistEtaPt[imult][imc]->GetXaxis()->SetTitle("#eta");
      fHistEtaPt[imult][imc]->GetYaxis()->SetTitle("P_{T}");
      fOutput->Add(fHistEtaPt[imult][imc]);
      
      sprintf(histname,"fHistNEvents_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"Number of Events and Number Passing Cuts %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistNEvents[imult][imc]=new TH1F(histname,histtitle,2,-0.5,1.5);
      fHistNEvents[imult][imc]->Sumw2();
      fHistNEvents[imult][imc]->GetXaxis()->SetTitle("Events,Passing Cuts");
      fHistNEvents[imult][imc]->GetYaxis()->SetTitle("Number of Events");
      fOutput->Add(fHistNEvents[imult][imc]);
      
      sprintf(histname,"fHistNTrigger_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"Number of Triggers %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistNTrigger[imult][imc]=new TH1F(histname,histtitle,fNTPtBins,-0.5,fNTPtBins-0.5);
      fHistNTrigger[imult][imc]->Sumw2();
      fHistNTrigger[imult][imc]->GetXaxis()->SetTitle("Trigger Number");
      fHistNTrigger[imult][imc]->GetYaxis()->SetTitle("Number of Triggers");
      fOutput->Add(fHistNTrigger[imult][imc]);
      
      sprintf(histname,"fHistNTriggerPt_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"P_{T} Weighted Number of Triggers %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistNTriggerPt[imult][imc]=new TH1F(histname,histtitle,fNTPtBins,-0.5,fNTPtBins-0.5);
      fHistNTriggerPt[imult][imc]->Sumw2();
      fHistNTriggerPt[imult][imc]->GetXaxis()->SetTitle("Trigger Number");
      fHistNTriggerPt[imult][imc]->GetYaxis()->SetTitle("Number of Triggers");
      fOutput->Add(fHistNTriggerPt[imult][imc]);
      
      sprintf(histname,"fHistNMix_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"Number of Mixed Events %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistNMix[imult][imc]=new TH1F(histname,histtitle,fNTPtBins,-0.5,fNTPtBins-0.5);
      fHistNMix[imult][imc]->Sumw2();
      fHistNMix[imult][imc]->GetXaxis()->SetTitle("Trigger Number");
      fHistNMix[imult][imc]->GetYaxis()->SetTitle("Number of Mixed Events");
      fOutput->Add(fHistNMix[imult][imc]);
      
      sprintf(histname,"fHistPhiEta_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"#phi-#eta distribution of tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPhiEta[imult][imc]=new TH3F(histname, histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,nptbins,fPtAssocArray);
      fHistPhiEta[imult][imc]->Sumw2();
      fHistPhiEta[imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistPhiEta[imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistPhiEta[imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistPhiEta[imult][imc]);
      
      sprintf(histname,"fHistPhiEtaPt_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"Pt Weighted #phi-#eta distribution of tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPhiEtaPt[imult][imc]=new TH3F(histname, histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,nptbins,fPtAssocArray);
      fHistPhiEtaPt[imult][imc]->Sumw2();
      fHistPhiEtaPt[imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistPhiEtaPt[imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistPhiEtaPt[imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistPhiEtaPt[imult][imc]);
      
      //if(fDEBUG)Printf("OutPut2");
      //Histograms with a trigger dependence
      
      for(int i=0;i<fNTPtBins;i++){
	for(int j=1;j<fNAPtBins;j++){
	  if(fabs(fPtTrigArray[i]-fPtAssocArray[j])<1E-5)lptbins=j;
	}
	//if(fDEBUG)Printf("Loop: %d Pt %3.2f",i,fPtTrigArray[i]/fPtBinWidth);
	
	//Ones with no centrality binning
	if(imult==0){
	  sprintf(histname,"fHistMultTrig_P%d%s",i,cmc1[imc]);
	  sprintf(histtitle,"Distrubition of number of tracks in triggered events with %3.1f<p_{T}^{Trig}<%3.1f%s",fPtTrigArray[i],fPtTrigArray[i+1],cmc2[imc]);
	  fHistMultTrig[i][imc]=new TH1F(histname,histtitle,2000,0,2000);
	  fHistMultTrig[i][imc]->Sumw2();
	  fHistMultTrig[i][imc]->GetXaxis()->SetTitle("Number of Tracks");
	  fHistMultTrig[i][imc]->GetYaxis()->SetTitle("Counts");
	  fOutput->Add(fHistMultTrig[i][imc]);
	}
	sprintf(histname,"fHistPtTrig_P%d_C%d%s",i,imult,cmc1[imc]);
	sprintf(histtitle,"P_{T} distribution in triggered events with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistPtTrig[i][imult][imc]=new TH1F(histname,histtitle,nptbins,fPtAssocArray);
	fHistPtTrig[i][imult][imc]->Sumw2();
	fHistPtTrig[i][imult][imc]->GetXaxis()->SetTitle("p_{T}");
	fHistPtTrig[i][imult][imc]->GetYaxis()->SetTitle("Counts");
	fOutput->Add(fHistPtTrig[i][imult][imc]);
	
	sprintf(histname,"fHistPhiTrig_P%d_C%d%s",i,imult,cmc1[imc]);
	sprintf(histtitle,"Phi Distribution of triggered events with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistPhiTrig[i][imult][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	fHistPhiTrig[i][imult][imc]->Sumw2();
	fHistPhiTrig[i][imult][imc]->GetXaxis()->SetTitle("#phi");
	fHistPhiTrig[i][imult][imc]->GetYaxis()->SetTitle("p_{T}");
	fOutput->Add(fHistPhiTrig[i][imult][imc]);
	
	sprintf(histname,"fHistPhiTrigPt_P%d_C%d%s",i,imult,cmc1[imc]);
	sprintf(histtitle,"P_{T} Weighted Phi Distribution of triggered events with %3.1f<p_{T}^{Trig}<%3.1f%s",fPtTrigArray[i],fPtTrigArray[i+1],cmc2[imc]);
	fHistPhiTrigPt[i][imult][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	fHistPhiTrigPt[i][imult][imc]->Sumw2();
	fHistPhiTrigPt[i][imult][imc]->GetXaxis()->SetTitle("#phi");
	fHistPhiTrigPt[i][imult][imc]->GetYaxis()->SetTitle("p_{T}");
	fOutput->Add(fHistPhiTrigPt[i][imult][imc]);
	
	sprintf(histname,"fHistEtaTrig_P%d_C%d%s",i,imult,cmc1[imc]);
	sprintf(histtitle,"Eta Distribution of triggered events with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistEtaTrig[i][imult][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	fHistEtaTrig[i][imult][imc]->Sumw2();
	fHistEtaTrig[i][imult][imc]->GetXaxis()->SetTitle("#eta");
	fHistEtaTrig[i][imult][imc]->GetYaxis()->SetTitle("p_{T}");
	fOutput->Add(fHistEtaTrig[i][imult][imc]);
	
	sprintf(histname,"fHistEtaTrigPt_P%d_C%d%s",i,imult,cmc1[imc]);
	sprintf(histtitle,"P_{T} Weighted Eta Distribution of triggered events with %3.1f<p_{T}^{Trig}<%3.1f%s",fPtTrigArray[i],fPtTrigArray[i+1],cmc2[imc]);
	fHistEtaTrigPt[i][imult][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	fHistEtaTrigPt[i][imult][imc]->Sumw2();
	fHistEtaTrigPt[i][imult][imc]->GetXaxis()->SetTitle("#eta");
	fHistEtaTrigPt[i][imult][imc]->GetYaxis()->SetTitle("p_{T}");
	fOutput->Add(fHistEtaTrigPt[i][imult][imc]);
	
	sprintf(histname,"fHistPhiEtaTrig_P%d_C%d%s",i,imult,cmc1[imc]);
	sprintf(histtitle,"#phi-#eta distribution in triggered events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistPhiEtaTrig[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
	fHistPhiEtaTrig[i][imult][imc]->Sumw2();
	fHistPhiEtaTrig[i][imult][imc]->GetXaxis()->SetTitle("#phi");
	fHistPhiEtaTrig[i][imult][imc]->GetYaxis()->SetTitle("#eta");
	fHistPhiEtaTrig[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
	fOutput->Add(fHistPhiEtaTrig[i][imult][imc]);

	sprintf(histname,"fHistXEN_P%d_C%d%s",i,imult,cmc1[imc]);
	sprintf(histtitle,"Near-Side X_{E} distribution for %3.1f<p_{T}^{Lead}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistXEN[i][imult][imc]=new TH1F(histname,histtitle,fNXEBins,fXEArray);
	fHistXEN[i][imult][imc]->Sumw2();
	fHistXEN[i][imult][imc]->GetXaxis()->SetTitle("X_{E}");
	fOutput->Add(fHistXEN[i][imult][imc]);

	sprintf(histname,"fHistXENMixed_P%d_C%d%s",i,imult,cmc1[imc]);
	sprintf(histtitle,"Mixed Near-Side X_{E} distribution for %3.1f<p_{T}^{Lead}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistXENMix[i][imult][imc]=new TH1F(histname,histtitle,fNXEBins,fXEArray);
	fHistXENMix[i][imult][imc]->Sumw2();
	fHistXENMix[i][imult][imc]->GetXaxis()->SetTitle("X_{E}");
	fOutput->Add(fHistXENMix[i][imult][imc]);
	
	sprintf(histname,"fHistXEA_P%d_C%d%s",i,imult,cmc1[imc]);
	sprintf(histtitle,"Away-Side X_{E} distribution for %3.1f<p_{T}^{Lead}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistXEA[i][imult][imc]=new TH1F(histname,histtitle,fNXEBins,fXEArray);
	fHistXEA[i][imult][imc]->Sumw2();
	fHistXEA[i][imult][imc]->GetXaxis()->SetTitle("X_{E}");
	fOutput->Add(fHistXEA[i][imult][imc]);

	sprintf(histname,"fHistXEAMixed_P%d_C%d%s",i,imult,cmc1[imc]);
	sprintf(histtitle,"Mixed Away-Side X_{E} distribution for %3.1f<p_{T}^{Lead}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistXEAMix[i][imult][imc]=new TH1F(histname,histtitle,fNXEBins,fXEArray);
	fHistXEAMix[i][imult][imc]->Sumw2();
	fHistXEAMix[i][imult][imc]->GetXaxis()->SetTitle("X_{E}");
	fOutput->Add(fHistXEAMix[i][imult][imc]);

	//signloop
	for(int isign=0;isign<3;isign++){
	  sprintf(histname,"fHistDeltaPhi_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"#Delta#phi Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhi[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhi[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhi[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhi[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhi[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaPhiPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"P_{T} Weighted #Delta#phi Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiPt[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaPhiMix_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"#Delta#phi Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiMix[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiMix[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiMix[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiMix[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiMix[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaPhiMixPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"P_{T} Weighted #Delta#phi Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiMixPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiMixPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiMixPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiMixPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiMixPt[i][imult][isign][imc]);
	  
	  //etaNear
	  sprintf(histname,"fHistDeltaEtaN_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Near-Side #Delta#eta Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaN[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaN[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaN[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaN[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaN[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaEtaNPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Near-Side P_{T} Weighted #Delta#eta Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNPt[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaEtaNMix_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Near-Side #Delta#eta Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNMix[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNMix[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNMix[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNMix[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNMix[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaEtaNMixPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Near-Side P_{T} Weighted #Delta#eta Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNMixPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNMixPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNMixPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNMixPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNMixPt[i][imult][isign][imc]);
	  
	  //Away Eta
	  sprintf(histname,"fHistDeltaEtaA_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Away-Side #Delta#eta Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaA[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaA[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaA[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaA[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaA[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaEtaAPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Away-Side P_{T} Weighted #Delta#eta Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAPt[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaEtaAMix_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Away-Side #Delta#eta Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAMix[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAMix[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAMix[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAMix[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAMix[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaEtaAMixPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Away-Side P_{T} Weighted #Delta#eta Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAMixPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAMixPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAMixPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAMixPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAMixPt[i][imult][isign][imc]);


      //====
	}//end isignloop
      sprintf(histname,"fHistDeltaPhiEta_P%d_C%d%s",i,imult,cmc1[imc]);
      sprintf(histtitle,"#Delta#phi-#Delta#eta %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEta[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEta[i][imult][imc]->Sumw2();
      fHistDeltaPhiEta[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEta[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEta[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEta[i][imult][imc]);
      
      sprintf(histname,"fHistDeltaPhiEtaMix_P%d_C%d%s",i,imult,cmc1[imc]);
      sprintf(histtitle,"#Delta#phi-#Delta#eta from Mixed Events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaMix[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEtaMix[i][imult][imc]->Sumw2();
      fHistDeltaPhiEtaMix[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEtaMix[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEtaMix[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEtaMix[i][imult][imc]);
      
      sprintf(histname,"fHistPhiEtaTrigPt_P%d_C%d%s",i,imult,cmc1[imc]);
      sprintf(histtitle,"P_{T}-Weighted #phi-#eta distribution in triggered events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPhiEtaTrigPt[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistPhiEtaTrigPt[i][imult][imc]->Sumw2();
      fHistPhiEtaTrigPt[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistPhiEtaTrigPt[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistPhiEtaTrigPt[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistPhiEtaTrigPt[i][imult][imc]);
    
      sprintf(histname,"fHistDeltaPhiEtaPt_P%d_C%d%s",i,imult,cmc1[imc]);
      sprintf(histtitle,"P_{T}-Weighted #Delta#phi-#Delta#eta %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaPt[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEtaPt[i][imult][imc]->Sumw2();
      fHistDeltaPhiEtaPt[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEtaPt[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEtaPt[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEtaPt[i][imult][imc]);
      
      sprintf(histname,"fHistDeltaPhiEtaMixPt_P%d_C%d%s",i,imult,cmc1[imc]);
      sprintf(histtitle,"P_{T}-Weighted #Delta#phi-#Delta#eta from Mixed Events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaMixPt[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEtaMixPt[i][imult][imc]->Sumw2();
      fHistDeltaPhiEtaMixPt[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEtaMixPt[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEtaMixPt[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEtaMixPt[i][imult][imc]);

      //Three-Particle Histograms
      for(int ipt=0;ipt<fNAPt3Bins;ipt++){
	for(int iSign=0;iSign<4;iSign++){
	  sprintf(histname,"fHistDeltaPhiPhi_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  sprintf(histtitle,"Raw #Delta#phi-#Delta#phi %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaPhiPhi[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinPhi3,fdPhiMin,fdPhiMax,fnBinPhi3,fdPhiMin,fdPhiMax);
	  fHistDeltaPhiPhi[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaPhiPhi[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#phi_{1}");
	  fHistDeltaPhiPhi[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#phi_{2}");
	  fOutput->Add(fHistDeltaPhiPhi[i][ipt][imult][iSign][imc]);
	  
	  sprintf(histname,"fHistDeltaPhiPhiMix_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  sprintf(histtitle,"Mixed #Delta#phi-#Delta#phi %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaPhiPhiMix[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinPhi3,fdPhiMin,fdPhiMax,fnBinPhi3,fdPhiMin,fdPhiMax);
	  fHistDeltaPhiPhiMix[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaPhiPhiMix[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#phi_{1}");
	  fHistDeltaPhiPhiMix[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#phi_{2}");
	  fOutput->Add(fHistDeltaPhiPhiMix[i][ipt][imult][iSign][imc]);

	  sprintf(histname,"fHistDeltaPhiPhiSS_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  sprintf(histtitle,"Soft-Soft #Delta#phi-#Delta#phi %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaPhiPhiSS[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinPhi3,fdPhiMin,fdPhiMax,fnBinPhi3,fdPhiMin,fdPhiMax);
	  fHistDeltaPhiPhiSS[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaPhiPhiSS[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#phi_{1}");
	  fHistDeltaPhiPhiSS[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#phi_{2}");
	  fOutput->Add(fHistDeltaPhiPhiSS[i][ipt][imult][iSign][imc]);

	  sprintf(histname,"fHistDeltaEtaEta_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  sprintf(histtitle,"Raw #Delta#eta-#Delta#eta %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaEtaEta[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinEta3,-etaEdge,etaEdge,fnBinEta3,-etaEdge,etaEdge);
	  fHistDeltaEtaEta[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaEtaEta[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#eta_{1}");
	  fHistDeltaEtaEta[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#eta_{2}");
	  fOutput->Add(fHistDeltaEtaEta[i][ipt][imult][iSign][imc]);

sprintf(histname,"fHistDeltaEtaEtaMix_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  sprintf(histtitle,"Mixed #Delta#eta-#Delta#eta %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinEta3,-etaEdge,etaEdge,fnBinEta3,-etaEdge,etaEdge);
	  fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#eta_{1}");
	  fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#eta_{2}");
	  fOutput->Add(fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]);

sprintf(histname,"fHistDeltaEtaEtaSS_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  sprintf(histtitle,"Soft-Soft #Delta#eta-#Delta#eta %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaEtaEtaSS[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinEta3,-etaEdge,etaEdge,fnBinEta3,-etaEdge,etaEdge);
	  fHistDeltaEtaEtaSS[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaEtaEtaSS[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#eta_{1}");
	  fHistDeltaEtaEtaSS[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#eta_{2}");
	  fOutput->Add(fHistDeltaEtaEtaSS[i][ipt][imult][iSign][imc]);

	}//iSign
      }//associated pt (ipt)
      }//pt loop (i) 
    }//centrality loop (imult)
  }//imc
  if(fDEBUG)Printf("OutPut Created");
}//CreateOutputObjects    

Int_t AliAnalysisTaskDiHadron::CheckVertex(const AliESDEvent *rESD){
  //checks whether the vertex passes cuts
  Int_t rGood=-1;
  Float_t vtx[3];
  vtx[0]=rESD->GetPrimaryVertex()->GetX();
  vtx[1]=rESD->GetPrimaryVertex()->GetY();
  vtx[2]=rESD->GetPrimaryVertex()->GetZ();
  if((vtx[0]*vtx[0]+vtx[1]*vtx[1])<9) rGood=0; //vertex out of beam pipe
  if(fabs(vtx[2])<fVzCut)rGood=0;//Vertex Z cut
  if(fDEBUG)Printf("vtxZ %f",vtx[2]);
  for(int i=0;i<fNVertexBins;i++){
    if(vtx[2]>fVertexArray[i]&&vtx[2]<=fVertexArray[i+1]&&rGood==0)rGood=i;
  }
  return rGood;
}

Int_t AliAnalysisTaskDiHadron::CheckTrigger(const AliESDEvent *rESD){
  //checks whether the trigger passes cuts
  Int_t rGood=0;
  TString trigID=rESD->GetFiredTriggerClasses();
  int count=0;
  char trigID2[50];
  int stop=0;//in as a safety

  for(int i=0;i<fNIDs;i++){
    if(stop==1)continue;
    for(int j=0;j<50;j++){
      if(fTrigIDArray[count]==',')trigID2[j]='\0';
      else if(fTrigIDArray[count]=='\0'){trigID2[j]='\0';stop=1;}
      else trigID2[j]=fTrigIDArray[count];
      count++;
      if(trigID2[j]=='\0') break;
      }
      if(trigID.Contains(trigID2)) rGood=1;
  }
    return rGood;
}

Int_t AliAnalysisTaskDiHadron::TrackCuts(const AliESDEvent *rESD, Float_t *rPt, Float_t *rEta, Float_t *rPhi, Short_t *rCharge, Float_t *rEff, Int_t **rPtAssoc3, Int_t *rNPtAssoc3, Int_t *rGoodTracks){
    //fills arrays with all of the tracks passing cuts
  rGoodTracks[0]=0;
  Int_t lead=0;
  Float_t leadPt=0;
  Int_t rTrack=fESD->GetNumberOfTracks();
  Float_t sPt, sEta, sPhi, sChi, sb[2], sbCov[3];
  Int_t sNcls, sNclsF, sITScls;
  Short_t sCharge;
  for(int iTrack=0;iTrack<rTrack;iTrack++){
    AliESDtrack *eSDtrack=rESD->GetTrack(iTrack);
    const AliExternalTrackParam *conTrack = eSDtrack->GetConstrainedParam();
    if(!conTrack)continue;
    sPt=conTrack->Pt();
    //if(fDEBUG)Printf("Pt%f",rPt);
    sEta=conTrack->Eta();
    sPhi=conTrack->Phi();
    sCharge=conTrack->Charge();
    if(sPhi<fdPhiMin)sPhi+=2*fPi;
    if(sPhi>fdPhiMax)sPhi-=2*fPi;
    if(sPt<fMinPtAssoc||sPt>fMaxPtAssoc)continue;//set Pt range
    if(fabs(sEta)>fEtaCut)continue;//set Eta Range
    if(!sCharge)continue;
    sNcls=eSDtrack->GetTPCNcls();
    //if(fDEBUG)Printf("NCLS%d",sNcls);
    if(sNcls<fMinClustersTPC)continue;
    sNclsF=eSDtrack->GetTPCNclsF();
    if((1.0*sNcls/sNclsF)<fMinClusterRatio)continue;//Clusters fit/ Possible
    sChi=(eSDtrack->GetTPCchi2())/sNcls;
    if(sChi>fMaxTPCchi2)continue;
    sITScls=eSDtrack->GetNcls(0);
    if(sITScls<fMinClustersITS)continue;
    eSDtrack->GetImpactParameters(sb,sbCov);
    if(!fDCA2D&&(sb[0]*sb[0]+sb[1]*sb[1])>(fMaxDCA*fMaxDCA))continue;//DCA cut
    if(fDCA2D==1&&(sb[0]*sb[0]/fMaxDCAXY/fMaxDCAXY+sb[1]*sb[1]/fMaxDCAZ/fMaxDCAZ)>1)continue;
    if(fDCA2D==2&&(0.35+0.42*std::pow(double(sPt),-0.9))<(sb[0]*sb[0]))continue;
    if(eSDtrack->GetKinkIndex(0)>0)continue;//removes kinked tracks
    if(!eSDtrack->GetStatus()&&AliESDtrack::kTPCrefit&&fTPCRefit)continue;//refit in TPC
    if((fITSRefit==1||(fITSRefit==2&&sPt>5))&&!eSDtrack->GetStatus()&&AliESDtrack::kITSrefit)continue;//refit of its tracks either for none,all, or >5 GeV/c
    if(fSPDCut&&!eSDtrack->HasPointOnITSLayer(0)&&!eSDtrack->HasPointOnITSLayer(1))continue;
    rPt[rGoodTracks[0]]=sPt;
    rEta[rGoodTracks[0]]=sEta;
    rPhi[rGoodTracks[0]]=sPhi;
    rCharge[rGoodTracks[0]]=sCharge;
    if(fEfficiencyCorr){
    if(sPt<fEffFitPt)rEff[rGoodTracks[0]]=1./fFitLow->Eval(sPt);
    else rEff[rGoodTracks[0]]=1./fFitHigh->Eval(sPt);
    }
    else rEff[rGoodTracks[0]]=1;
    if(sPt>leadPt)lead=rGoodTracks[0];
    //rPtAssoc3[rGoodTracks[0]]=new Int_t [10];
    rNPtAssoc3[rGoodTracks[0]]=0;
    for(int apt3=0;apt3<fNAPt3Bins;apt3++){
      if(sPt<fPtAssoc3Array2[apt3]&&sPt>=fPtAssoc3Array1[apt3]){
	rPtAssoc3[rGoodTracks[0]][rNPtAssoc3[rGoodTracks[0]]]=apt3;
	rNPtAssoc3[rGoodTracks[0]]++;
      }
    }

    rGoodTracks[0]++;
    
  }
  return lead;
}

Int_t AliAnalysisTaskDiHadron::TrackCutsMC(AliMCEvent *rMC, Float_t *rPt, Float_t *rEta, Float_t *rPhi, Short_t *rCharge, Float_t *rEff, Int_t **rPtAssoc3, Int_t *rNPtAssoc3, Int_t *rGoodTracks){
//Fills Arrays of MC particles
  rGoodTracks[1]=0;
  AliStack *rStack=rMC->Stack();
  Int_t rTrack=rStack->GetNtrack();
  Float_t sPt, sEta, sPhi;
  Short_t sCharge;
  Int_t lead=0;
  Float_t leadPt=0;
  for(int iTrack=0;iTrack<rTrack;iTrack++){
    TParticle *rParticle=rStack->Particle(iTrack);
    sPt=rParticle->Pt();
    //if(fDEBUG)Printf("MCPt%f",rPt);
    sEta=rParticle->Eta();
    sPhi=rParticle->Phi();
    sCharge=rMC->GetTrack(iTrack)->Charge();
    if(sPhi<fdPhiMin)sPhi+=2*fPi;
    if(sPhi>fdPhiMax)sPhi-=2*fPi;
    if(sPt<fMinPtAssoc||sPt>fMaxPtAssoc)continue;//set Pt range
    if(fabs(sEta)>fEtaCut)continue;//set Eta Range
    if(!rStack->IsPhysicalPrimary(iTrack))continue;//primary particles only
    if(!sCharge)continue;//only charged particles kept
    rPt[rGoodTracks[1]]=sPt;
    rEta[rGoodTracks[1]]=sEta;
    rPhi[rGoodTracks[1]]=sPhi;
    rCharge[rGoodTracks[1]]=sCharge;
    rEff[rGoodTracks[1]]=1;
    if(sPt>leadPt)lead=rGoodTracks[1];
    rNPtAssoc3[rGoodTracks[1]]=0;
    for(int apt3=0;apt3<fNAPt3Bins;apt3++){
      if(sPt<fPtAssoc3Array2[apt3]&&sPt>=fPtAssoc3Array1[apt3]){
	rPtAssoc3[rGoodTracks[1]][rNPtAssoc3[rGoodTracks[1]]]=apt3;
	rNPtAssoc3[rGoodTracks[1]]++;
      }
    }
    rGoodTracks[1]++;
  }
  return lead;
}
//------------------------------------------------------------
void AliAnalysisTaskDiHadron::Exec(Option_t *)
{ 
  //Main executable
  if(fDEBUG)Printf("Exec");

  const int nTPtBins=fNTPtBins;
  const int nCentBins=fNCentBins;
  for(int ievent=0;ievent<=1;ievent++){

  if(!fESD&&ievent==0){
    if(fDEBUG)Printf("Error: fESD not found");
    break;
  }
  if(!fMC&&ievent==1){
    break;
  }
  if(ievent==1&&!fMCHistos)break;//break out if MC event and we don't have fill of those set
  //Secondary check
  if(ievent==0){
    if(fESD->GetNumberOfTracks()<=0){
      if(fDEBUG)Printf("Error: no tracks");
      break;
    }
  }
  //The previous check doesn't seem to work as a fMC is bad not NULL
  if(ievent==1){
    if(fMC->GetNumberOfTracks()<=0){
      if(fDEBUG)Printf("<=0 MCTracks");
      break;
    }
  }
  
  //Check for Trigger only on real data
  if(!fMC){
    if(!CheckTrigger(fESD)) break;
  }
  //I'll only cut on the reconstructed vertex since these are the events that will be used
  int vertexBin;
  vertexBin=CheckVertex(fESD);
  //else vertexBin=CheckVertex(fMC);
  if(vertexBin<0)break;

  Int_t nGoodTracks[2]={0,0}, nTriggers[nTPtBins][nCentBins][2];
  Int_t nTrack;
  if(!ievent)nTrack=fESD->GetNumberOfTracks();
  else nTrack=fMC->Stack()->GetNtrack();
 
  Float_t tdPhi, tdEta, tXE;
  Float_t tdPhi2, tdEta2;
  ftPhi=new Float_t [nTrack];
  ftEta=new Float_t [nTrack];
  ftPt=new Float_t [nTrack];
  ftCharge=new Short_t [nTrack];
  ftEff=new Float_t [nTrack];
  ftPtAssoc3=new Int_t *[nTrack];
  for(int i=0;i<nTrack;i++){
    ftPtAssoc3[i]=new Int_t [10];
  }
  ftNPtAssoc3=new Int_t [nTrack];
  Short_t sign;

  //will need to do something exta for the effieiency once it comes from embedding
  for(int i=0;i<fNTPtBins;i++){
    for(int c=0;c<fNCentBins;c++){
    nTriggers[i][c][ievent]=0;
    }
  }
  Int_t tMult=fESD->GetMultiplicity()->GetNumberOfTracklets();//I think this is the correct multiplicity to use
  if(fDEBUG)Printf("Mult%d",tMult);
  
  //Decide what multiplicy bins are filled with this event, note set to max of 4 as I didn't think more then 2 overlapping bins likely at one time, easliy changed
  Int_t multArray[4]={0,0,0,0};
  Int_t maxArray=0;
  for(int imult=0;imult<fNCentBins;imult++){
    if(tMult>=fCentArrayMin[imult]&&tMult<fCentArrayMax[imult]){
      multArray[maxArray]=imult;
      maxArray++;
    }
  }
  if(fDEBUG)Printf("maxArray%d",maxArray);
  //Set Efficiency for the centrality bin (lowest bin used in array if multiple overlap)
  for(int ipar=0;ipar<fNFitLowParam;ipar++){
    fFitLow->SetParameter(ipar,fFitLowParam[multArray[0]*fNCentBins+ipar]);
  }
   for(int ipar=0;ipar<fNFitHighParam;ipar++){
    fFitHigh->SetParameter(ipar,fFitHighParam[multArray[0]*fNCentBins+ipar]);
  }
  fHistMult[ievent]->Fill(tMult);
  for(int c=0;c<maxArray;c++){fHistNEvents[multArray[c]][ievent]->Fill(0);}//count the number of events used
  Int_t leadPart;
  
  //returns arrays filled up to nGoodTracks with tracks passing cuts
  if(!ievent)leadPart=TrackCuts(fESD,ftPt,ftEta,ftPhi,ftCharge,ftEff,ftPtAssoc3,ftNPtAssoc3,nGoodTracks);
  else leadPart=TrackCutsMC(fMC,ftPt,ftEta,ftPhi,ftCharge,ftEff,ftPtAssoc3,ftNPtAssoc3,nGoodTracks);
  int nearEta=0,NearXE=0;
  int nearEta2=0;
  if(fDEBUG)Printf("Track Loop");
  for(int iTrack=0;iTrack<nGoodTracks[ievent];iTrack++){
    if(fDEBUG)Printf("Track%d Pt%f",iTrack,ftPt[iTrack]);
    //if(ftPhi[iTrack]<fdPhiMin)ftPhi[iTrack]+=2*fPi;
    //if(ftPhi[iTrack]>fdPhiMax)ftPhi[iTrack]-=2*fPi;
    for(int c=0;c<maxArray;c++){
      fHistPt[multArray[c]][ievent]->Fill(ftPt[iTrack],ftEff[iTrack]);
    fHistPtEff[multArray[c]][ievent]->Fill(ftPt[iTrack]);
    fHistPhi[multArray[c]][ievent]->Fill(ftPhi[iTrack],ftPt[iTrack],ftEff[iTrack]);
    fHistPhiPt[multArray[c]][ievent]->Fill(ftPhi[iTrack],ftPt[iTrack],ftPt[iTrack]*ftEff[iTrack]);
    fHistEta[multArray[c]][ievent]->Fill(ftEta[iTrack],ftPt[iTrack],ftEff[iTrack]);
    fHistEtaPt[multArray[c]][ievent]->Fill(ftEta[iTrack],ftPt[iTrack],ftPt[iTrack]*ftEff[iTrack]);
    fHistPhiEta[multArray[c]][ievent]->Fill(ftPhi[iTrack],ftEta[iTrack],ftPt[iTrack],ftEff[iTrack]);
    fHistPhiEtaPt[multArray[c]][ievent]->Fill(ftPhi[iTrack],ftEta[iTrack],ftPt[iTrack],ftPt[iTrack]*ftEff[iTrack]);
    }
    for(int i=0;i<fNTPtBins;i++){
      if(ftPt[iTrack]>fPtTrigArray[i]&&ftPt[iTrack]<fPtTrigArray[i+1]&&fabs(ftEta[iTrack])<fTrigEtaCut){
	if(fDEBUG)Printf("In %fpt%f",fPtTrigArray[i],fPtTrigArray[i+1]);
	fHistMultTrig[i][ievent]->Fill(tMult);
	for(int c=0;c<maxArray;c++){
	  nTriggers[i][multArray[c]][ievent]++;
	  fHistNTrigger[multArray[c]][ievent]->Fill(i);
	  fHistNTriggerPt[multArray[c]][ievent]->Fill(i,ftPt[iTrack]);
	}
	if(fDEBUG)Printf("Assiciated Particle Loop");
	for(int iTrack2=0;iTrack2<nGoodTracks[ievent];iTrack2++){
	  if(iTrack==iTrack2) continue;
	  if(ftPt[iTrack2]>ftPt[iTrack])continue;
	  tdPhi=ftPhi[iTrack]-ftPhi[iTrack2];
	  if(tdPhi<-fPi)tdPhi+=2*fPi;
	  if(tdPhi>fPi)tdPhi-=2*fPi;
	  if(fabs(tdPhi)<fNearPhiCut)nearEta=1;
	  else nearEta=0;
	  if(fabs(tdPhi)<fXECut)NearXE=1;
	  else NearXE=0;
	  if(fabs(tdPhi)<(fPi/2))tdEta=ftEta[iTrack]-ftEta[iTrack2];
	  else tdEta=ftEta[iTrack]+ftEta[iTrack2];
	  if(tdPhi<fdPhiMin)tdPhi+=2*fPi;
	  if(tdPhi>fdPhiMax)tdPhi-=2*fPi;
	  if((ftCharge[iTrack]<0&&ftCharge[iTrack2]<0)||(ftCharge[iTrack]>0&&ftCharge[iTrack2]>0))sign=1;
	  else sign=2;
	  if(fDEBUG) Printf("dPhi %f  dEta %f",tdPhi,tdEta);
	  for(int c=0;c<maxArray;c++){//loop over multiplicity bins
	    fHistPtTrig[i][multArray[c]][ievent]->Fill(ftPt[iTrack2],ftEff[iTrack2]);
	    fHistPhiTrig[i][multArray[c]][ievent]->Fill(ftPhi[iTrack2],ftPt[iTrack2],ftEff[iTrack2]);
	    fHistPhiTrigPt[i][multArray[c]][ievent]->Fill(ftPhi[iTrack2],ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]);
	    fHistEtaTrig[i][multArray[c]][ievent]->Fill(ftEta[iTrack2],ftPt[iTrack2],ftEff[iTrack2]);
	    fHistEtaTrigPt[i][multArray[c]][ievent]->Fill(ftEta[iTrack2],ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]);

	    fHistPhiEtaTrig[i][multArray[c]][ievent]->Fill(ftPhi[iTrack2],ftEta[iTrack2],ftPt[iTrack2],ftEff[iTrack2]);
	    fHistPhiEtaTrigPt[i][multArray[c]][ievent]->Fill(ftPhi[iTrack2],ftEta[iTrack2],ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]);
	    fHistDeltaPhi[i][multArray[c]][0][ievent]->Fill(tdPhi,ftPt[iTrack2],ftEff[iTrack2]);
	    fHistDeltaPhiPt[i][multArray[c]][0][ievent]->Fill(tdPhi,ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]);
	    fHistDeltaPhi[i][multArray[c]][sign][ievent]->Fill(tdPhi,ftPt[iTrack2],ftEff[iTrack2]);
	    fHistDeltaPhiPt[i][multArray[c]][sign][ievent]->Fill(tdPhi,ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]);

	    if(nearEta){
	    fHistDeltaEtaN[i][multArray[c]][0][ievent]->Fill(tdEta,ftPt[iTrack2],ftEff[iTrack2]);
	    fHistDeltaEtaNPt[i][multArray[c]][0][ievent]->Fill(tdEta,ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]);
	    fHistDeltaEtaN[i][multArray[c]][sign][ievent]->Fill(tdEta,ftPt[iTrack2],ftEff[iTrack2]);
	    fHistDeltaEtaNPt[i][multArray[c]][sign][ievent]->Fill(tdEta,ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]);
	    }
	    else{
	    fHistDeltaEtaA[i][multArray[c]][0][ievent]->Fill(tdEta,ftPt[iTrack2],ftEff[iTrack2]);
	    fHistDeltaEtaAPt[i][multArray[c]][0][ievent]->Fill(tdEta,ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]);
	    fHistDeltaEtaA[i][multArray[c]][sign][ievent]->Fill(tdEta,ftPt[iTrack2],ftEff[iTrack2]);
	    fHistDeltaEtaAPt[i][multArray[c]][sign][ievent]->Fill(tdEta,ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]);
	    }
	    fHistDeltaPhiEta[i][multArray[c]][ievent]->Fill(tdPhi,tdEta,ftPt[iTrack2],ftEff[iTrack2]);
	    fHistDeltaPhiEtaPt[i][multArray[c]][ievent]->Fill(tdPhi,tdEta,ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]);
	    
	    //only fill these if trigger particle is the leading particle
	    if(iTrack==leadPart){
	      if(NearXE){
		tXE=ftPt[iTrack2]*cos(tdPhi)/ftPt[iTrack];
		fHistXEN[i][multArray[c]][ievent]->Fill(tXE,ftEff[iTrack2]);
	      }
	      else{
		tXE=ftPt[iTrack2]*cos(tdPhi+fPi)/ftPt[iTrack];
		fHistXEA[i][multArray[c]][ievent]->Fill(tXE,ftEff[iTrack2]);
	      }
	    }

	  }//Centrality loop (c)

	  //3-particle Correlations
	  for(int iTrack3=0;iTrack3<nGoodTracks[ievent];iTrack3++){
	    if(iTrack2==iTrack3)continue;
	    if(ftPt[iTrack3]>ftPt[iTrack])continue;
	    tdPhi2=ftPhi[iTrack]-ftPhi[iTrack3];
	    if(tdPhi2<-fPi)tdPhi2+=2*fPi;
	    if(tdPhi2>fPi)tdPhi2-=2*fPi;
	    if(fabs(tdPhi2)<fNearPhiCut&&nearEta==1)nearEta2=1;
	    else nearEta2=0;
	    //if(fabs(tdPhi)<fXECut)NearXE=1;
	    //else NearXE=0;
	    if(fabs(tdPhi2)<(fPi/2))tdEta2=ftEta[iTrack]-ftEta[iTrack3];
	    else tdEta2=ftEta[iTrack]+ftEta[iTrack3];
	    if(tdPhi2<fdPhiMin)tdPhi2+=2*fPi;
	    if(tdPhi2>fdPhiMax)tdPhi2-=2*fPi;
	    // if((ftCharge[iTrack]<0&&ftCharge[iTrack2]<0)||(ftCharge[iTrack]>0&&ftCharge[iTrack2]>0))sign=1;
	    if((ftCharge[iTrack]<0&&ftCharge[iTrack2]<0&&ftCharge[iTrack3]<0)||(ftCharge[iTrack]>0&&ftCharge[iTrack2]>0&&ftCharge[iTrack3]>0))sign=1;
	    else if((ftCharge[iTrack3]<0&&ftCharge[iTrack2]<0)||(ftCharge[iTrack3]>0&&ftCharge[iTrack2]>0))sign=2;
	    else sign=3;
	    for(int e=0;e<ftNPtAssoc3[iTrack2];e++){//check associated pT bin
	      for(int f=0;f<ftNPtAssoc3[iTrack3];f++){
		if(ftPtAssoc3[iTrack2][e]==ftPtAssoc3[iTrack3][f]){
		  for(int c=0;c<maxArray;c++){//loop over multiplicity bins
		    fHistDeltaPhiPhi[i][ftPtAssoc3[iTrack2][e]][multArray[c]][0][ievent]->Fill(tdPhi,tdPhi2,ftEff[iTrack2]*ftEff[iTrack3]);
		    fHistDeltaPhiPhi[i][ftPtAssoc3[iTrack2][e]][multArray[c]][sign][ievent]->Fill(tdPhi2,tdPhi,ftEff[iTrack2]*ftEff[iTrack3]);
		 

		    if(nearEta2){
		      fHistDeltaEtaEta[i][ftPtAssoc3[iTrack2][e]][multArray[c]][0][ievent]->Fill(tdEta,tdEta2,ftEff[iTrack2]*ftEff[iTrack3]);
		      fHistDeltaEtaEta[i][ftPtAssoc3[iTrack2][e]][multArray[c]][sign][ievent]->Fill(tdEta,tdEta2,ftEff[iTrack2]*ftEff[iTrack3]);
		    }
		  }//multiplicity loop (c)
		}
	      }
	    }//track checking loops
	  }//iTrack3
	}//iTrack2 (associated track loop)
	
	if(fDEBUG)Printf("Mixed Event Loop");
	for(int c=0;c<maxArray;c++){
	  int d=multArray[c];//Centrality bin we are in
	  if(fMixEnd[d][vertexBin][ievent]>=0){//check if there are any mixed events for this bin
	    for(int imix=0;imix<=fMixEnd[d][vertexBin][ievent];imix++){//loop over the stored mixed events
	      fHistNMix[d][ievent]->Fill(i);
	      for(int iTrack2=0;iTrack2<fMixTrack[imix][d][vertexBin][ievent];iTrack2++){
		if(ftPt[iTrack]<fMPt[imix][d][vertexBin][ievent][iTrack2])continue;
		tdPhi=ftPhi[iTrack]-fMPhi[imix][d][vertexBin][ievent][iTrack2];
		if(tdPhi<-fPi)tdPhi+=2*fPi;
		if(tdPhi>fPi)tdPhi-=2*fPi;
		if(fabs(tdPhi)<fNearPhiCut)nearEta=1;
		else nearEta=0;
		if(fabs(tdPhi)<fXECut)NearXE=1;
		else NearXE=0;
		if(fabs(tdPhi)<(fPi/2))tdEta=ftEta[iTrack]-fMEta[imix][d][vertexBin][ievent][iTrack2];
		else tdEta=ftEta[iTrack]+fMEta[imix][d][vertexBin][ievent][iTrack2];
		if(tdPhi<fdPhiMin)tdPhi+=2*fPi;	
		if(tdPhi>fdPhiMax)tdPhi-=2*fPi;
		if((ftCharge[iTrack]<0&&fMCharge[imix][d][vertexBin][ievent][iTrack2]<0)||(ftCharge[iTrack]>0&&fMCharge[imix][d][vertexBin][ievent][iTrack2]>0))sign=1;
		else sign=2;

		fHistDeltaPhiMix[i][d][0][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]);
		fHistDeltaPhiMixPt[i][d][0][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMPt[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack2]);
		fHistDeltaPhiMix[i][d][sign][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]);
		fHistDeltaPhiMixPt[i][d][sign][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMPt[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack2]);
		if(nearEta){
		  fHistDeltaEtaNMix[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]);
		fHistDeltaEtaNMixPt[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMPt[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack2]);
		fHistDeltaEtaNMix[i][d][sign][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]);
		fHistDeltaEtaNMixPt[i][d][sign][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMPt[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack2]);
		}
		else{
		  fHistDeltaEtaAMix[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]);
		fHistDeltaEtaAMixPt[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMPt[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack2]);
		fHistDeltaEtaAMix[i][d][sign][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]);
		fHistDeltaEtaAMixPt[i][d][sign][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMPt[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack2]);
		}	

		fHistDeltaPhiEtaMix[i][d][ievent]->Fill(tdPhi,tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]);
		fHistDeltaPhiEtaMixPt[i][d][ievent]->Fill(tdPhi,tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMPt[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack2]);

		if(iTrack==leadPart){
		  if(NearXE){
		    tXE=fMPt[imix][d][vertexBin][ievent][iTrack2]*cos(tdPhi)/ftPt[iTrack];
		    fHistXENMix[i][d][ievent]->Fill(tXE,fMEff[imix][d][vertexBin][ievent][iTrack2]);
		  }
		  else{
		    tXE=fMPt[imix][d][vertexBin][ievent][iTrack2]*cos(tdPhi+fPi)/ftPt[iTrack];
		    fHistXEAMix[i][multArray[c]][ievent]->Fill(tXE,fMEff[imix][d][vertexBin][ievent][iTrack2]);
		  }
		}
		//3-particle correlation soft-soft term (both associated from the same event)
		for(int iTrack3=0;iTrack3<fMixTrack[imix][d][vertexBin][ievent];iTrack3++){
		  if(iTrack3==iTrack2)continue;
		  if(ftPt[iTrack]<fMPt[imix][d][vertexBin][ievent][iTrack3])continue;
		tdPhi2=ftPhi[iTrack]-fMPhi[imix][d][vertexBin][ievent][iTrack3];
		if(tdPhi2<-fPi)tdPhi2+=2*fPi;
		if(tdPhi2>fPi)tdPhi2-=2*fPi;
		if(fabs(tdPhi2)<fNearPhiCut&&nearEta)nearEta2=1;
		else nearEta2=0;
		if(fabs(tdPhi2)<(fPi/2))tdEta2=ftEta[iTrack]-fMEta[imix][d][vertexBin][ievent][iTrack3];
		else tdEta2=ftEta[iTrack]+fMEta[imix][d][vertexBin][ievent][iTrack3];
		if(tdPhi2<fdPhiMin)tdPhi2+=2*fPi;	
		if(tdPhi2>fdPhiMax)tdPhi2-=2*fPi;
		//if((ftCharge[iTrack]<0&&fMCharge[imix][d][vertexBin][ievent][iTrack2]<0)||(ftCharge[iTrack]>0&&fMCharge[imix][d][vertexBin][ievent][iTrack2]>0))sign=1;
		//else sign=2;
		if((ftCharge[iTrack]<0&&fMCharge[imix][d][vertexBin][ievent][iTrack2]<0&&fMCharge[imix][d][vertexBin][ievent][iTrack3]<0)||(ftCharge[iTrack]>0&&fMCharge[imix][d][vertexBin][ievent][iTrack2]>0&&fMCharge[imix][d][vertexBin][ievent][iTrack3]>0))sign=1;
		else if((fMCharge[imix][d][vertexBin][ievent][iTrack3]<0&&fMCharge[imix][d][vertexBin][ievent][iTrack2]<0)||(fMCharge[imix][d][vertexBin][ievent][iTrack3]>0&&fMCharge[imix][d][vertexBin][ievent][iTrack2]>0))sign=2;
		else sign=3;
		for(int e=0;e<fMNPtAssoc3[imix][d][vertexBin][ievent][iTrack2];e++){//check associated pT bin
		  for(int f=0;f<fMNPtAssoc3[imix][d][vertexBin][ievent][iTrack3];f++){
		    if(fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]==fMPtAssoc3[imix][d][vertexBin][ievent][f][iTrack3]){
		      fHistDeltaPhiPhiSS[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack3]);
		      fHistDeltaPhiPhiSS[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack3]); 

		      if(nearEta2){
			fHistDeltaEtaEtaSS[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack3]);
			fHistDeltaEtaEtaSS[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack3]);
		      }//near-side
		    }
		  }
		}//associated pt bin
		}//iTrack3

		//3-particle mixed event (associated from different events)
		//for(int imix2=0;imix2<=fMixEnd[d][vertexBin][ievent];imix2++){//loop over the stored mixed events
		//if(imix2==imix)continue;
		  int imix2=imix+1;
		  if(imix2>=fMixEnd[d][vertexBin][ievent])imix2=0;
		  if(imix2==imix)continue;//will kill it when there is only 1 mixed event (remember to scale by 1 less then the number of mixed events in others (number of mixed-mixed is 2*(mixed-1)))
		  for(int iTrack3=0;iTrack3<fMixTrack[imix2][d][vertexBin][ievent];iTrack3++){
		    if(ftPt[iTrack]<fMPt[imix2][d][vertexBin][ievent][iTrack3])continue;
		    tdPhi2=ftPhi[iTrack]-fMPhi[imix2][d][vertexBin][ievent][iTrack3];
		    if(tdPhi2<-fPi)tdPhi2+=2*fPi;
		    if(tdPhi2>fPi)tdPhi2-=2*fPi;
		    if(fabs(tdPhi2)<fNearPhiCut&&nearEta)nearEta2=1;
		    else nearEta2=0;
		    if(fabs(tdPhi2)<(fPi/2))tdEta2=ftEta[iTrack]-fMEta[imix2][d][vertexBin][ievent][iTrack3];
		    else tdEta2=ftEta[iTrack]+fMEta[imix2][d][vertexBin][ievent][iTrack3];
		    if(tdPhi2<fdPhiMin)tdPhi2+=2*fPi;	
		    if(tdPhi2>fdPhiMax)tdPhi2-=2*fPi;
		    //if((ftCharge[iTrack]<0&&fMCharge[imix2][d][vertexBin][ievent][iTrack2]<0)||(ftCharge[iTrack]>0&&fMCharge[imix2][d][vertexBin][ievent][iTrack2]>0))sign=1;
		    //else sign=2;
		    if((ftCharge[iTrack]<0&&fMCharge[imix][d][vertexBin][ievent][iTrack2]<0&&fMCharge[imix2][d][vertexBin][ievent][iTrack3]<0)||(ftCharge[iTrack]>0&&fMCharge[imix][d][vertexBin][ievent][iTrack2]>0&&fMCharge[imix2][d][vertexBin][ievent][iTrack3]>0))sign=1;
		    else if((fMCharge[imix2][d][vertexBin][ievent][iTrack3]<0&&fMCharge[imix][d][vertexBin][ievent][iTrack2]<0)||(fMCharge[imix2][d][vertexBin][ievent][iTrack3]>0&&fMCharge[imix][d][vertexBin][ievent][iTrack2]>0))sign=2;
		    else sign=3;
		    for(int e=0;e<fMNPtAssoc3[imix][d][vertexBin][ievent][iTrack2];e++){//check associated pT bin
		      for(int f=0;f<fMNPtAssoc3[imix2][d][vertexBin][ievent][iTrack3];f++){
			if(fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]==fMPtAssoc3[imix2][d][vertexBin][ievent][f][iTrack3]){
			  fHistDeltaPhiPhiMix[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]);
			  fHistDeltaPhiPhiMix[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]); 
			  fHistDeltaPhiPhiMix[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]);//free factor of 2 in statistics
			  fHistDeltaPhiPhiMix[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]); 
		      if(nearEta2){
			fHistDeltaEtaEtaMix[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]);
			fHistDeltaEtaEtaMix[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]);
			fHistDeltaEtaEtaMix[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta2,tdEta,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]);
			fHistDeltaEtaEtaMix[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdEta2,tdEta,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]);
		      }//near-side
			}
		      }
		    }//associated pt bin
		  }//iTrack3

	      }//iTrack2
	      }//imix
	    }//fMixEnd
	  }//Centrality bins (c)
	}//pt trig cuts
      }//i Pt Trig
    }//itrack    
  
  //now store this event for mixing (using these dynamic arrays instead of the other fix to save memory)
    if(fDEBUG)Printf("Store Event For Mixing");
  for(int c=0;c<maxArray;c++){//loops over centrality bins
    int d=multArray[c];//too many nested arrays looked confusing d=which centrality bin
    if(fMixEnd[d][vertexBin][ievent]<(fNMix-1))fMixEnd[d][vertexBin][ievent]++;
    if(fMixPointer[d][vertexBin][ievent]<(fNMix-1))fMixPointer[d][vertexBin][ievent]++;
    else fMixPointer[d][vertexBin][ievent]=0;
    int e=fMixPointer[d][vertexBin][ievent];//nested arrays (e is event number in pool)
    delete [] fMPt[e][d][vertexBin][ievent];
    delete [] fMPhi[e][d][vertexBin][ievent];
    delete [] fMEta[e][d][vertexBin][ievent];
    delete [] fMCharge[e][d][vertexBin][ievent];
    delete [] fMEff[e][d][vertexBin][ievent];
    delete [] fMNPtAssoc3[e][d][vertexBin][ievent];
    for(int jj=0;jj<10;jj++){
    delete [] fMPtAssoc3[e][d][vertexBin][ievent][jj];
    }
    fMPt[e][d][vertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
    fMPhi[e][d][vertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
    fMEta[e][d][vertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
    fMCharge[e][d][vertexBin][ievent]=new Short_t [nGoodTracks[ievent]];
    fMEff[e][d][vertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
    fMixTrack[e][d][vertexBin][ievent]=nGoodTracks[ievent];
    fMNPtAssoc3[e][d][vertexBin][ievent]=new Short_t [nGoodTracks[ievent]];
    for(int jj=0;jj<10;jj++){
    fMPtAssoc3[e][d][vertexBin][ievent][jj]=new Short_t [nGoodTracks[ievent]];
    }

    for(int iTrack=0;iTrack<nGoodTracks[ievent];iTrack++){
      fMPt[e][d][vertexBin][ievent][iTrack]=ftPt[iTrack];
      fMPhi[e][d][vertexBin][ievent][iTrack]=ftPhi[iTrack];
      fMEta[e][d][vertexBin][ievent][iTrack]=ftEta[iTrack];
      fMCharge[e][d][vertexBin][ievent][iTrack]=ftCharge[iTrack];
      fMEff[e][d][vertexBin][ievent][iTrack]=ftEff[iTrack];
      fMNPtAssoc3[e][d][vertexBin][ievent][iTrack]=ftNPtAssoc3[iTrack];
      // fMPtAssoc3[e][d][vertexBin][ievent][iTrack]=new Int_t [ftNPtAssoc3[iTrack]];
      for(int jj=0;jj<ftNPtAssoc3[iTrack];jj++){
	//if(fDEBUG) Printf("%d",ftPtAssoc3[iTrack][jj]);
	fMPtAssoc3[e][d][vertexBin][ievent][jj][iTrack]=ftPtAssoc3[iTrack][jj];
      }
      
    }//iTracks
  }//Centrality (c)
  }//ievent
  //track=0;
  //track2=0;

  PostData(0, fOutput);
  //get rid of these arrays from memory
  delete [] ftPhi;
  delete [] ftEta;
  delete [] ftPt;
  delete [] ftCharge;
  delete [] ftEff;
  delete [] ftNPtAssoc3;
  delete [] ftPtAssoc3;
  ftPhi=NULL;
  ftEta=NULL;
  ftPt=NULL;
  ftCharge=NULL;
  ftEff=NULL;
  ftNPtAssoc3=NULL;
  ftPtAssoc3=NULL;

}//Exec

//---------------------------------------------------
void AliAnalysisTaskDiHadron::Terminate(Option_t *){
  //Terminates the code, frees up memory
  for(int ii=0;ii<fNMix;ii++){
    for(int cc=0;cc<fNCentBins;cc++){
      for(int vtx=0;vtx<fNVertexBins;vtx++){
	for(int jj=0;jj<2;jj++){
	  delete [] fMPt[ii][cc][vtx][jj];
	  delete [] fMPhi[ii][cc][vtx][jj];
	  delete [] fMEta[ii][cc][vtx][jj];
	  delete [] fMCharge[ii][cc][vtx][jj];
	  delete [] fMEff[ii][cc][vtx][jj];
	  delete [] fMNPtAssoc3[ii][cc][vtx][jj];
	  for(int qq=0;qq<10;qq++){
	    delete [] fMPtAssoc3[ii][cc][vtx][jj][qq];
	    fMPtAssoc3[ii][cc][vtx][jj][qq]=NULL;
	  }
	  fMPt[ii][cc][vtx][jj]=NULL;
	  fMPhi[ii][cc][vtx][jj]=NULL;
	  fMEta[ii][cc][vtx][jj]=NULL;
	  fMCharge[ii][cc][vtx][jj]=NULL;
	  fMEff[ii][cc][vtx][jj]=NULL;
	  fMNPtAssoc3[ii][cc][vtx][jj]=NULL;

	}
      }
    }
  }
  delete [] fFitLowParam;
  delete [] fFitHighParam;
  delete [] fPtTrigArray;
  delete [] fPtAssocArray;
  delete [] fPtAssoc3Array1;
  delete [] fPtAssoc3Array2;
  delete [] fCentArrayMin;
  delete [] fCentArrayMax;
  delete [] fXEArray;
  fFitLowParam=NULL;
  fFitHighParam=NULL;
  fPtTrigArray=NULL;
  fPtAssocArray=NULL;
  fPtAssoc3Array1=NULL;
  fPtAssoc3Array2=NULL;
  fCentArrayMin=NULL;
  fCentArrayMax=NULL;
  Printf("Terminate AliAnalysisTaskDiHadron");
}
