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
AliAnalysisTask(name,""), fESD(0), fMC(0), fOutput(0),fMinClustersTPC(0),fMinClusterRatio(0),fMaxTPCchi2(0),fMinClustersITS(0),fEtaCut(0),fTrigEtaCut(0),fNearPhiCut(0),fXECut(0),fMaxDCA(0),fMaxDCAXY(0),fMaxDCAZ(0),fDCA2D(0),fTPCRefit(0),fITSRefit(0),fSPDCut(0),fMinPtAssoc(0),fMaxPtAssoc(0),fVzCut(0),fEfficiencyCorr(0),DEBUG(0),fnBinPhi(0),fnBinEta(0),fnBinPhiEtaPhi(0),fnBinPhiEtaEta(0),fnBinPhi3(0),fnBinEta3(0),fPi(3.1415926535898),fdPhiMin(0),fdPhiMax(0),fNTPtBins(0),fNMix(0),fNCentBins(0),fNAPtBins(0),fNAPt3Bins(0),fNVertexBins(0),fNXEBins(0),fNIDs(0),fEffFitPt(0),fNFitLowParam(0),fNFitHighParam(0),fMCHistos(0),fNFitLow(0),fNFitHigh(0),fFitLow(NULL),fFitHigh(NULL),fFitLowParam(NULL),fFitHighParam(NULL),fPtTrigArray(NULL),fPtAssocArray(NULL),fPtAssoc3Array1(NULL),fPtAssoc3Array2(NULL),fCentArrayMin(NULL),fCentArrayMax(NULL),fXEArray(NULL),fTrigIDArray(NULL),tPhi(NULL),tEta(NULL),tPt(NULL),tCharge(NULL),tEff(NULL),tPtAssoc3(NULL),tNPtAssoc3(NULL)
  
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
  // tPhi=NULL;
  //  tEta=NULL;
  // tPt=NULL;

  //fPi=3.1415926535898;


  }
//--------------------------------------
void AliAnalysisTaskDiHadron::SetCuts(Int_t MinClustersTPC,  Float_t MinClusterRatio, Float_t MaxTPCchi2, Int_t MinClustersITS, Float_t EtaCut, Float_t TrigEtaCut, Float_t NearPhiCut, Float_t XECut, Float_t MaxDCA, Float_t MaxDCAXY, Float_t MaxDCAZ, Int_t DCA2D, Int_t TPCRefit, Int_t ITSRefit, Int_t SPDCut, Float_t MinPtAssoc, Float_t MaxPtAssoc, Float_t VzCut, Int_t NIDs, char **TrigIDArray){
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
   fTrigIDArray=(char**)TrigIDArray;

   //  Printf("TPC%d R%2.1f Chi%2.1f ITS%d Eta%2.1f Near%2.1f NearX%2.1f DCA%2.1f DCAR%2.1f DCAZ%2.1f DCAM%d TPCFit%d ITSFit%d SPD%d MinPt%2.1f MaxPt%2.1f Vz%2.1f Trig_%s",fMinClustersTPC,fMinClusterRatio,fMaxTPCchi2,fMinClustersITS,fEtaCut,fNearPhiCut,fXECut,fMaxDCA,fMaxDCAXY,fMaxDCAZ,fDCA2D,fTPCRefit,fITSRefit,fSPDCut,fMinPtAssoc,fMaxPtAssoc,fVzCut,fTrigIDArray[0]);
}
//--------------------------------------------------------
void AliAnalysisTaskDiHadron::SetOptions(Int_t EfficiencyCorr, Int_t fDEBUG,  Int_t MCHistos){
  fEfficiencyCorr=EfficiencyCorr;
  DEBUG=fDEBUG;
  fMCHistos=MCHistos;
  //Printf("Eff%d DEBUG%d MCHistos%d",fEfficiencyCorr,DEBUG,fMCHistos);
}
//------------------------------------------------------
void AliAnalysisTaskDiHadron::SetBins(Int_t nBinPhi, Int_t nBinEta, Int_t nBinPhiEtaPhi, Int_t nBinPhiEtaEta, Int_t nBinPhi3, Int_t nBinEta3,Float_t dPhiMin, Float_t dPhiMax, Int_t NTPtBins, Int_t NMixBins, Int_t NCentBins,Int_t NAPtBins, Int_t NAPt3Bins, Int_t NVertexBins, Int_t NXEBins,Float_t *PtTrigArray, Float_t *PtAssocArray,Float_t *PtAssoc3Array1, Float_t *PtAssoc3Array2, Int_t *CentArrayMin, Int_t *CentArrayMax, Float_t *XEArray){
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
  fPtTrigArray=(Float_t*)PtTrigArray;
  fPtAssocArray=(Float_t*)PtAssocArray;
  fPtAssoc3Array1=(Float_t*)PtAssoc3Array1;
  fPtAssoc3Array2=(Float_t*)PtAssoc3Array2;
  fCentArrayMin=(Int_t*)CentArrayMin;
  fCentArrayMax=(Int_t*)CentArrayMax;
  fXEArray=(Float_t*)XEArray;
 for(int i=0;i<=fNVertexBins;i++)fVertexArray[i]=(2.*i/fNVertexBins-1)*fVzCut;
 //Printf("BinPhi%d BinEta%d BinPhiEta%d %d Bin3Phi%d Bin3Eta%d VertexArray%2.1f %2.1f",fnBinPhi,fnBinEta,fnBinPhiEtaPhi,fnBinPhiEtaEta,fnBinPhi3,fnBinEta3,fVertexArray[0],fVertexArray[fNVertexBins]);
}
//-------------------------------------------------------
void AliAnalysisTaskDiHadron::SetEfficiencies(Float_t EffFitPt, TF1 *FitLow, TF1 *FitHigh, Int_t NFitLowParam, Int_t NFitHighParam, Float_t *FitLowParam, Float_t *FitHighParam){
  fEffFitPt=EffFitPt;
  fFitLow=(TF1*)FitLow;
  fFitHigh=(TF1*)FitHigh;
  fNFitLowParam=NFitLowParam;
  fNFitHighParam=NFitHighParam;
  // Printf("NCent%d NFitLow%d NFitHigh%d",NCent,fNFitLowParam,fNFitHighParam);
  //  Printf("FitLowParm%2.4f",FitLowParam[0]);
  fFitLowParam=(Float_t*)FitLowParam;
  fFitHighParam=(Float_t*)FitHighParam;
  //Printf("fEffFitPt%2.1f FitParam%2.1f",fEffFitPt,fFitLowParam[0]);
}

//-----------------------------------------------------------
void AliAnalysisTaskDiHadron::ConnectInputData(Option_t *){
  //Connect to ESD
  if(DEBUG)Printf("Connecting");
   TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
   if (!tree&&DEBUG) {Printf("ERROR: Could not read chain from input slot 0");} 
  else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!esdH&&DEBUG) Printf("ERROR: Could not get ESDInputHandler");
    else fESD = esdH->GetEvent();
    
    //MC Data handler (so one can calcualte eff)
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    if(mcH)fMC=mcH->MCEvent();
  }
  if(DEBUG)Printf("Connected");
}
  
//---------------------------------------------------------  
void AliAnalysisTaskDiHadron::CreateOutputObjects(){
  if(DEBUG)Printf("Output");
  //Creates the histograms and list
  fOutput=new TList();
  fOutput->SetName(GetName());
  char histname[100];
  char histtitle[200];
  int nptbins=fNAPtBins;
  int lptbins=0;
  char *cmc1[2]={"","_MC"};
  char *cmc2[2]={""," MC"};
  char *sign1[3]={"","_LS","_ULS"};
  char *sign2[3]={""," Like-Sign"," Unlike-Sign"};
  char *sign31[4]={"","_LS","_ULT","_ULA"};
  char *sign32[4]={""," Like-Sign"," Trigger-Diff"," Assoc-Diff"};
  Float_t EtaEdge=2*fEtaCut;
  Float_t PhiArray[fnBinPhi+1];
  Float_t EtaArray[fnBinEta+1];
  Float_t PhiEtaArrayPhi[fnBinPhiEtaPhi+1];
  Float_t PhiEtaArrayEta[fnBinPhiEtaEta+1];
  for(int iphi=0;iphi<=fnBinPhi;iphi++){
    PhiArray[iphi]=fdPhiMin+iphi*2*fPi/fnBinPhi;
  }
  for(int ieta=0;ieta<=fnBinEta;ieta++){
    EtaArray[ieta]=-EtaEdge+ieta*2*EtaEdge/fnBinEta;
  }
  for(int iphi=0;iphi<=fnBinPhiEtaPhi;iphi++){
    PhiEtaArrayPhi[iphi]=fdPhiMin+iphi*2*fPi/fnBinPhiEtaPhi;
  }
  for(int ieta=0;ieta<=fnBinPhiEtaEta;ieta++){
    PhiEtaArrayEta[ieta]=-EtaEdge+ieta*2*EtaEdge/fnBinPhiEtaEta;
  }
  for(int imc=0;imc<=1;imc++){//MC loop
    if(imc==1&&!fMCHistos) continue;
    //Create the histograms
    sprintf(histname,"fHistMult%s",cmc1[imc]);
    sprintf(histtitle,"Multiplicity%s",cmc2[imc]);
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
      fHistPhi[imult][imc]=new TH2F(histname,histtitle,fnBinPhi,PhiArray,nptbins,fPtAssocArray);
      fHistPhi[imult][imc]->Sumw2();
      fHistPhi[imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistPhi[imult][imc]->GetYaxis()->SetTitle("P_{T}");
      fOutput->Add(fHistPhi[imult][imc]);
      
      sprintf(histname,"fHistPhiPt_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"P_{T} weighted #phi Distribution of tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPhiPt[imult][imc]=new TH2F(histname,histtitle,fnBinPhi,PhiArray,nptbins,fPtAssocArray);
      fHistPhiPt[imult][imc]->Sumw2();
      fHistPhiPt[imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistPhiPt[imult][imc]->GetYaxis()->SetTitle("P_{T}");
      fOutput->Add(fHistPhiPt[imult][imc]);
      
      sprintf(histname,"fHistEta_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"#eta Distribution of Tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistEta[imult][imc]=new TH2F(histname,histtitle,fnBinEta,EtaArray,nptbins,fPtAssocArray);
      fHistEta[imult][imc]->Sumw2();
      fHistEta[imult][imc]->GetXaxis()->SetTitle("#eta");
      fHistEta[imult][imc]->GetYaxis()->SetTitle("P_{T}");
      fOutput->Add(fHistEta[imult][imc]);
      
      sprintf(histname,"fHistEtaPt_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"P_{T} weighted #eta Distribution of tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistEtaPt[imult][imc]=new TH2F(histname,histtitle,fnBinEta,EtaArray,nptbins,fPtAssocArray);
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
      fHistPhiEta[imult][imc]=new TH3F(histname, histtitle,fnBinPhiEtaPhi,PhiEtaArrayPhi,fnBinPhiEtaEta,PhiEtaArrayEta,nptbins,fPtAssocArray);
      fHistPhiEta[imult][imc]->Sumw2();
      fHistPhiEta[imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistPhiEta[imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistPhiEta[imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistPhiEta[imult][imc]);
      
      sprintf(histname,"fHistPhiEtaPt_C%d%s",imult,cmc1[imc]);
      sprintf(histtitle,"Pt Weighted #phi-#eta distribution of tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPhiEtaPt[imult][imc]=new TH3F(histname, histtitle,fnBinPhiEtaPhi,PhiEtaArrayPhi,fnBinPhiEtaEta,PhiEtaArrayEta,nptbins,fPtAssocArray);
      fHistPhiEtaPt[imult][imc]->Sumw2();
      fHistPhiEtaPt[imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistPhiEtaPt[imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistPhiEtaPt[imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistPhiEtaPt[imult][imc]);
      
      //if(DEBUG)Printf("OutPut2");
      //Histograms with a trigger dependence
      
      for(int i=0;i<fNTPtBins;i++){
	for(int j=1;j<fNAPtBins;j++){
	  if(fPtTrigArray[i]==fPtAssocArray[j])lptbins=j;
	}
	//if(DEBUG)Printf("Loop: %d Pt %3.2f",i,fPtTrigArray[i]/fPtBinWidth);
	
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
	fHistPhiTrig[i][imult][imc]=new TH2F(histname,histtitle,fnBinPhi,PhiArray,lptbins,fPtAssocArray);
	fHistPhiTrig[i][imult][imc]->Sumw2();
	fHistPhiTrig[i][imult][imc]->GetXaxis()->SetTitle("#phi");
	fHistPhiTrig[i][imult][imc]->GetYaxis()->SetTitle("p_{T}");
	fOutput->Add(fHistPhiTrig[i][imult][imc]);
	
	sprintf(histname,"fHistPhiTrigPt_P%d_C%d%s",i,imult,cmc1[imc]);
	sprintf(histtitle,"P_{T} Weighted Phi Distribution of triggered events with %3.1f<p_{T}^{Trig}<%3.1f%s",fPtTrigArray[i],fPtTrigArray[i+1],cmc2[imc]);
	fHistPhiTrigPt[i][imult][imc]=new TH2F(histname,histtitle,fnBinPhi,PhiArray,lptbins,fPtAssocArray);
	fHistPhiTrigPt[i][imult][imc]->Sumw2();
	fHistPhiTrigPt[i][imult][imc]->GetXaxis()->SetTitle("#phi");
	fHistPhiTrigPt[i][imult][imc]->GetYaxis()->SetTitle("p_{T}");
	fOutput->Add(fHistPhiTrigPt[i][imult][imc]);
	
	sprintf(histname,"fHistEtaTrig_P%d_C%d%s",i,imult,cmc1[imc]);
	sprintf(histtitle,"Eta Distribution of triggered events with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistEtaTrig[i][imult][imc]=new TH2F(histname,histtitle,fnBinEta,EtaArray,lptbins,fPtAssocArray);
	fHistEtaTrig[i][imult][imc]->Sumw2();
	fHistEtaTrig[i][imult][imc]->GetXaxis()->SetTitle("#eta");
	fHistEtaTrig[i][imult][imc]->GetYaxis()->SetTitle("p_{T}");
	fOutput->Add(fHistEtaTrig[i][imult][imc]);
	
	sprintf(histname,"fHistEtaTrigPt_P%d_C%d%s",i,imult,cmc1[imc]);
	sprintf(histtitle,"P_{T} Weighted Eta Distribution of triggered events with %3.1f<p_{T}^{Trig}<%3.1f%s",fPtTrigArray[i],fPtTrigArray[i+1],cmc2[imc]);
	fHistEtaTrigPt[i][imult][imc]=new TH2F(histname,histtitle,fnBinEta,EtaArray,lptbins,fPtAssocArray);
	fHistEtaTrigPt[i][imult][imc]->Sumw2();
	fHistEtaTrigPt[i][imult][imc]->GetXaxis()->SetTitle("#eta");
	fHistEtaTrigPt[i][imult][imc]->GetYaxis()->SetTitle("p_{T}");
	fOutput->Add(fHistEtaTrigPt[i][imult][imc]);
	
	sprintf(histname,"fHistPhiEtaTrig_P%d_C%d%s",i,imult,cmc1[imc]);
	sprintf(histtitle,"#phi-#eta distribution in triggered events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistPhiEtaTrig[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,PhiEtaArrayPhi,fnBinPhiEtaEta,PhiEtaArrayEta,lptbins,fPtAssocArray);
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
	  fHistDeltaPhi[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,PhiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhi[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhi[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhi[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhi[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaPhiPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"P_{T} Weighted #Delta#phi Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,PhiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiPt[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaPhiMix_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"#Delta#phi Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiMix[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,PhiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiMix[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiMix[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiMix[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiMix[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaPhiMixPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"P_{T} Weighted #Delta#phi Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiMixPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,PhiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiMixPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiMixPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiMixPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiMixPt[i][imult][isign][imc]);
	  
	  //etaNear
	  sprintf(histname,"fHistDeltaEtaN_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Near-Side #Delta#eta Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaN[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,EtaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaN[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaN[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaN[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaN[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaEtaNPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Near-Side P_{T} Weighted #Delta#eta Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,EtaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNPt[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaEtaNMix_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Near-Side #Delta#eta Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNMix[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,EtaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNMix[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNMix[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNMix[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNMix[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaEtaNMixPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Near-Side P_{T} Weighted #Delta#eta Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNMixPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,EtaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNMixPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNMixPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNMixPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNMixPt[i][imult][isign][imc]);
	  
	  //Away Eta
	  sprintf(histname,"fHistDeltaEtaA_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Away-Side #Delta#eta Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaA[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,EtaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaA[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaA[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaA[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaA[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaEtaAPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Away-Side P_{T} Weighted #Delta#eta Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,EtaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAPt[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaEtaAMix_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Away-Side #Delta#eta Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAMix[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,EtaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAMix[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAMix[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAMix[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAMix[i][imult][isign][imc]);
	  
	  sprintf(histname,"fHistDeltaEtaAMixPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  sprintf(histtitle,"Away-Side P_{T} Weighted #Delta#eta Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAMixPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,EtaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAMixPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAMixPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAMixPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAMixPt[i][imult][isign][imc]);


      //====
	}//end isignloop
      sprintf(histname,"fHistDeltaPhiEta_P%d_C%d%s",i,imult,cmc1[imc]);
      sprintf(histtitle,"#Delta#phi-#Delta#eta %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEta[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,PhiEtaArrayPhi,fnBinPhiEtaEta,PhiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEta[i][imult][imc]->Sumw2();
      fHistDeltaPhiEta[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEta[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEta[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEta[i][imult][imc]);
      
      sprintf(histname,"fHistDeltaPhiEtaMix_P%d_C%d%s",i,imult,cmc1[imc]);
      sprintf(histtitle,"#Delta#phi-#Delta#eta from Mixed Events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaMix[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,PhiEtaArrayPhi,fnBinPhiEtaEta,PhiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEtaMix[i][imult][imc]->Sumw2();
      fHistDeltaPhiEtaMix[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEtaMix[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEtaMix[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEtaMix[i][imult][imc]);
      
      sprintf(histname,"fHistPhiEtaTrigPt_P%d_C%d%s",i,imult,cmc1[imc]);
      sprintf(histtitle,"P_{T}-Weighted #phi-#eta distribution in triggered events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPhiEtaTrigPt[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,PhiEtaArrayPhi,fnBinPhiEtaEta,PhiEtaArrayEta,lptbins,fPtAssocArray);
      fHistPhiEtaTrigPt[i][imult][imc]->Sumw2();
      fHistPhiEtaTrigPt[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistPhiEtaTrigPt[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistPhiEtaTrigPt[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistPhiEtaTrigPt[i][imult][imc]);
    
      sprintf(histname,"fHistDeltaPhiEtaPt_P%d_C%d%s",i,imult,cmc1[imc]);
      sprintf(histtitle,"P_{T}-Weighted #Delta#phi-#Delta#eta %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaPt[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,PhiEtaArrayPhi,fnBinPhiEtaEta,PhiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEtaPt[i][imult][imc]->Sumw2();
      fHistDeltaPhiEtaPt[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEtaPt[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEtaPt[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEtaPt[i][imult][imc]);
      
      sprintf(histname,"fHistDeltaPhiEtaMixPt_P%d_C%d%s",i,imult,cmc1[imc]);
      sprintf(histtitle,"P_{T}-Weighted #Delta#phi-#Delta#eta from Mixed Events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaMixPt[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,PhiEtaArrayPhi,fnBinPhiEtaEta,PhiEtaArrayEta,lptbins,fPtAssocArray);
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
	  fHistDeltaEtaEta[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinEta3,-EtaEdge,EtaEdge,fnBinEta3,-EtaEdge,EtaEdge);
	  fHistDeltaEtaEta[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaEtaEta[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#eta_{1}");
	  fHistDeltaEtaEta[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#eta_{2}");
	  fOutput->Add(fHistDeltaEtaEta[i][ipt][imult][iSign][imc]);

sprintf(histname,"fHistDeltaEtaEtaMix_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  sprintf(histtitle,"Mixed #Delta#eta-#Delta#eta %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinEta3,-EtaEdge,EtaEdge,fnBinEta3,-EtaEdge,EtaEdge);
	  fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#eta_{1}");
	  fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#eta_{2}");
	  fOutput->Add(fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]);

sprintf(histname,"fHistDeltaEtaEtaSS_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  sprintf(histtitle,"Soft-Soft #Delta#eta-#Delta#eta %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaEtaEtaSS[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinEta3,-EtaEdge,EtaEdge,fnBinEta3,-EtaEdge,EtaEdge);
	  fHistDeltaEtaEtaSS[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaEtaEtaSS[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#eta_{1}");
	  fHistDeltaEtaEtaSS[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#eta_{2}");
	  fOutput->Add(fHistDeltaEtaEtaSS[i][ipt][imult][iSign][imc]);

	}//iSign
      }//associated pt (ipt)
      }//pt loop (i) 
    }//centrality loop (imult)
  }//imc
  if(DEBUG)Printf("OutPut Created");
}//CreateOutputObjects    

Int_t AliAnalysisTaskDiHadron::CheckVertex(AliESDEvent *rESD){
  Int_t rGood=-1;
  Float_t Vtx[3];
  Vtx[0]=rESD->GetPrimaryVertex()->GetX();
  Vtx[1]=rESD->GetPrimaryVertex()->GetY();
  Vtx[2]=rESD->GetPrimaryVertex()->GetZ();
  if((Vtx[0]*Vtx[0]+Vtx[1]*Vtx[1])<9) rGood=0; //vertex out of beam pipe
  if(fabs(Vtx[2])<fVzCut)rGood=0;//Vertex Z cut
  if(DEBUG)Printf("VtxZ %f",Vtx[2]);
  for(int i=0;i<fNVertexBins;i++){
    if(Vtx[2]>fVertexArray[i]&&Vtx[2]<=fVertexArray[i+1]&&rGood==0)rGood=i;
  }
  return rGood;
}

Int_t AliAnalysisTaskDiHadron::CheckTrigger(AliESDEvent *rESD){
  Int_t rGood=0;
  TString trigID=rESD->GetFiredTriggerClasses();
  for(int i=0;i<fNIDs;i++){
    if(trigID.Contains(fTrigIDArray[i])) rGood=1;
  }
    return rGood;
}

Int_t AliAnalysisTaskDiHadron::TrackCuts(AliESDEvent *rESD, Float_t *rPt, Float_t *rEta, Float_t *rPhi, Short_t *rCharge, Float_t *rEff, Int_t **rPtAssoc3, Int_t *rNPtAssoc3, Int_t *rGoodTracks){
    //fills arrays with all of the tracks passing cuts
  rGoodTracks[0]=0;
  Int_t Lead=0;
  Float_t LeadPt=0;
  Int_t rTrack=fESD->GetNumberOfTracks();
  Float_t sPt, sEta, sPhi, sChi, sb[2], sbCov[3];
  Int_t sNcls, sNclsF, sITScls;
  Short_t sCharge;
  for(int iTrack=0;iTrack<rTrack;iTrack++){
    AliESDtrack *ESDtrack=rESD->GetTrack(iTrack);
    const AliExternalTrackParam *ConTrack = ESDtrack->GetConstrainedParam();
    if(!ConTrack)continue;
    sPt=ConTrack->Pt();
    //if(DEBUG)Printf("Pt%f",rPt);
    sEta=ConTrack->Eta();
    sPhi=ConTrack->Phi();
    sCharge=ConTrack->Charge();
    if(sPhi<fdPhiMin)sPhi+=2*fPi;
    if(sPhi>fdPhiMax)sPhi-=2*fPi;
    if(sPt<fMinPtAssoc||sPt>fMaxPtAssoc)continue;//set Pt range
    if(fabs(sEta)>fEtaCut)continue;//set Eta Range
    if(!sCharge)continue;
    sNcls=ESDtrack->GetTPCNcls();
    //if(DEBUG)Printf("NCLS%d",sNcls);
    if(sNcls<fMinClustersTPC)continue;
    sNclsF=ESDtrack->GetTPCNclsF();
    if((1.0*sNcls/sNclsF)<fMinClusterRatio)continue;//Clusters fit/ Possible
    sChi=(ESDtrack->GetTPCchi2())/sNcls;
    if(sChi>fMaxTPCchi2)continue;
    sITScls=ESDtrack->GetNcls(0);
    if(sITScls<fMinClustersITS)continue;
    ESDtrack->GetImpactParameters(sb,sbCov);
    if(!fDCA2D&&(sb[0]*sb[0]+sb[1]*sb[1])>(fMaxDCA*fMaxDCA))continue;//DCA cut
    if(fDCA2D==1&&(sb[0]*sb[0]/fMaxDCAXY/fMaxDCAXY+sb[1]*sb[1]/fMaxDCAZ/fMaxDCAZ)>1)continue;
    if(fDCA2D==2&&(0.35+0.42*std::pow(double(sPt),-0.9))<(sb[0]*sb[0]))continue;
    if(ESDtrack->GetKinkIndex(0)>0)continue;//removes kinked tracks
    if(!ESDtrack->GetStatus()&AliESDtrack::kTPCrefit&&fTPCRefit)continue;//refit in TPC
    if((fITSRefit==1||(fITSRefit==2&&sPt>5))&&!ESDtrack->GetStatus()&AliESDtrack::kITSrefit)continue;//refit of its tracks either for none,all, or >5 GeV/c
    if(fSPDCut&&!ESDtrack->HasPointOnITSLayer(0)&&!ESDtrack->HasPointOnITSLayer(1))continue;
    rPt[rGoodTracks[0]]=sPt;
    rEta[rGoodTracks[0]]=sEta;
    rPhi[rGoodTracks[0]]=sPhi;
    rCharge[rGoodTracks[0]]=sCharge;
    if(fEfficiencyCorr){
    if(sPt<fEffFitPt)rEff[rGoodTracks[0]]=1./fFitLow->Eval(sPt);
    else rEff[rGoodTracks[0]]=1./fFitHigh->Eval(sPt);
    }
    else rEff[rGoodTracks[0]]=1;
    if(sPt>LeadPt)Lead=rGoodTracks[0];
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
  return Lead;
}

Int_t AliAnalysisTaskDiHadron::TrackCutsMC(AliMCEvent *rMC, Float_t *rPt, Float_t *rEta, Float_t *rPhi, Short_t *rCharge, Float_t *rEff, Int_t **rPtAssoc3, Int_t *rNPtAssoc3, Int_t *rGoodTracks){//Fills Arrays of MC particles
  rGoodTracks[1]=0;
  AliStack *rStack=rMC->Stack();
  Int_t rTrack=rStack->GetNtrack();
  Float_t sPt, sEta, sPhi;
  Short_t sCharge;
  Int_t Lead=0;
  Float_t LeadPt=0;
  for(int iTrack=0;iTrack<rTrack;iTrack++){
    TParticle *rParticle=rStack->Particle(iTrack);
    sPt=rParticle->Pt();
    //if(DEBUG)Printf("MCPt%f",rPt);
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
    if(sPt>LeadPt)Lead=rGoodTracks[1];
    rNPtAssoc3[rGoodTracks[1]]=0;
    for(int apt3=0;apt3<fNAPt3Bins;apt3++){
      if(sPt<fPtAssoc3Array2[apt3]&&sPt>=fPtAssoc3Array1[apt3]){
	rPtAssoc3[rGoodTracks[1]][rNPtAssoc3[rGoodTracks[1]]]=apt3;
	rNPtAssoc3[rGoodTracks[1]]++;
      }
    }
    rGoodTracks[1]++;
  }
  return Lead;
}
//------------------------------------------------------------
void AliAnalysisTaskDiHadron::Exec(Option_t *)
{ 
  if(DEBUG)Printf("Exec");

  const int NTPtBins=fNTPtBins;
  const int NCentBins=fNCentBins;
  for(int ievent=0;ievent<=1;ievent++){

  if(!fESD&&ievent==0){
    if(DEBUG)Printf("Error: fESD not found");
    break;
  }
  if(!fMC&&ievent==1){
    break;
  }
  if(ievent==1&&!fMCHistos)break;//break out if MC event and we don't have fill of those set
  //Secondary check
  if(ievent==0){
    if(fESD->GetNumberOfTracks()<=0){
      if(DEBUG)Printf("Error: no tracks");
      break;
    }
  }
  //The previous check doesn't seem to work as a fMC is bad not NULL
  if(ievent==1){
    if(fMC->GetNumberOfTracks()<=0){
      if(DEBUG)Printf("<=0 MCTracks");
      break;
    }
  }
  
  //Check for Trigger only on real data
  if(!fMC){
    if(!CheckTrigger(fESD)) break;
  }
  //I'll only cut on the reconstructed vertex since these are the events that will be used
  int VertexBin;
  VertexBin=CheckVertex(fESD);
  //else VertexBin=CheckVertex(fMC);
  if(VertexBin<0)break;

  Int_t nGoodTracks[2]={0,0}, nTriggers[NTPtBins][NCentBins][2];
  Int_t nTrack;
  if(!ievent)nTrack=fESD->GetNumberOfTracks();
  else nTrack=fMC->Stack()->GetNtrack();
 
  Float_t tdPhi, tdEta, tXE;
  Float_t tdPhi2, tdEta2;
  tPhi=new Float_t [nTrack];
  tEta=new Float_t [nTrack];
  tPt=new Float_t [nTrack];
  tCharge=new Short_t [nTrack];
  tEff=new Float_t [nTrack];
  tPtAssoc3=new Int_t *[nTrack];
  for(int i=0;i<nTrack;i++){
    tPtAssoc3[i]=new Int_t [10];
  }
  tNPtAssoc3=new Int_t [nTrack];
  Short_t Sign;

  //will need to do something exta for the effieiency once it comes from embedding
  for(int i=0;i<fNTPtBins;i++){
    for(int c=0;c<fNCentBins;c++){
    nTriggers[i][c][ievent]=0;
    }
  }
  Int_t tMult=fESD->GetMultiplicity()->GetNumberOfTracklets();//I think this is the correct multiplicity to use
  if(DEBUG)Printf("Mult%d",tMult);
  
  //Decide what multiplicy bins are filled with this event, note set to max of 4 as I didn't think more then 2 overlapping bins likely at one time, easliy changed
  Int_t MultArray[4]={0,0,0,0};
  Int_t MaxArray=0;
  for(int imult=0;imult<fNCentBins;imult++){
    if(tMult>=fCentArrayMin[imult]&&tMult<fCentArrayMax[imult]){
      MultArray[MaxArray]=imult;
      MaxArray++;
    }
  }
  if(DEBUG)Printf("MaxArray%d",MaxArray);
  //Set Efficiency for the centrality bin (lowest bin used in array if multiple overlap)
  for(int ipar=0;ipar<fNFitLowParam;ipar++){
    fFitLow->SetParameter(ipar,fFitLowParam[MultArray[0]*fNCentBins+ipar]);
  }
   for(int ipar=0;ipar<fNFitHighParam;ipar++){
    fFitHigh->SetParameter(ipar,fFitHighParam[MultArray[0]*fNCentBins+ipar]);
  }
  fHistMult[ievent]->Fill(tMult);
  for(int c=0;c<MaxArray;c++){fHistNEvents[MultArray[c]][ievent]->Fill(0);}//count the number of events used
  Int_t LeadPart;
  
  //returns arrays filled up to nGoodTracks with tracks passing cuts
  if(!ievent)LeadPart=TrackCuts(fESD,tPt,tEta,tPhi,tCharge,tEff,tPtAssoc3,tNPtAssoc3,nGoodTracks);
  else LeadPart=TrackCutsMC(fMC,tPt,tEta,tPhi,tCharge,tEff,tPtAssoc3,tNPtAssoc3,nGoodTracks);
  int NearEta=0,NearXE=0;
  int NearEta2=0;
  if(DEBUG)Printf("Track Loop");
  for(int iTrack=0;iTrack<nGoodTracks[ievent];iTrack++){
    if(DEBUG)Printf("Track%d Pt%f",iTrack,tPt[iTrack]);
    //if(tPhi[iTrack]<fdPhiMin)tPhi[iTrack]+=2*fPi;
    //if(tPhi[iTrack]>fdPhiMax)tPhi[iTrack]-=2*fPi;
    for(int c=0;c<MaxArray;c++){
      fHistPt[MultArray[c]][ievent]->Fill(tPt[iTrack],tEff[iTrack]);
    fHistPtEff[MultArray[c]][ievent]->Fill(tPt[iTrack]);
    fHistPhi[MultArray[c]][ievent]->Fill(tPhi[iTrack],tPt[iTrack],tEff[iTrack]);
    fHistPhiPt[MultArray[c]][ievent]->Fill(tPhi[iTrack],tPt[iTrack],tPt[iTrack]*tEff[iTrack]);
    fHistEta[MultArray[c]][ievent]->Fill(tEta[iTrack],tPt[iTrack],tEff[iTrack]);
    fHistEtaPt[MultArray[c]][ievent]->Fill(tEta[iTrack],tPt[iTrack],tPt[iTrack]*tEff[iTrack]);
    fHistPhiEta[MultArray[c]][ievent]->Fill(tPhi[iTrack],tEta[iTrack],tPt[iTrack],tEff[iTrack]);
    fHistPhiEtaPt[MultArray[c]][ievent]->Fill(tPhi[iTrack],tEta[iTrack],tPt[iTrack],tPt[iTrack]*tEff[iTrack]);
    }
    for(int i=0;i<fNTPtBins;i++){
      if(tPt[iTrack]>fPtTrigArray[i]&&tPt[iTrack]<fPtTrigArray[i+1]&&fabs(tEta[iTrack])<fTrigEtaCut){
	if(DEBUG)Printf("In %fpt%f",fPtTrigArray[i],fPtTrigArray[i+1]);
	fHistMultTrig[i][ievent]->Fill(tMult);
	for(int c=0;c<MaxArray;c++){
	  nTriggers[i][MultArray[c]][ievent]++;
	  fHistNTrigger[MultArray[c]][ievent]->Fill(i);
	  fHistNTriggerPt[MultArray[c]][ievent]->Fill(i,tPt[iTrack]);
	}
	if(DEBUG)Printf("Assiciated Particle Loop");
	for(int iTrack2=0;iTrack2<nGoodTracks[ievent];iTrack2++){
	  if(iTrack==iTrack2) continue;
	  if(tPt[iTrack2]>tPt[iTrack])continue;
	  tdPhi=tPhi[iTrack]-tPhi[iTrack2];
	  if(tdPhi<-fPi)tdPhi+=2*fPi;
	  if(tdPhi>fPi)tdPhi-=2*fPi;
	  if(fabs(tdPhi)<fNearPhiCut)NearEta=1;
	  else NearEta=0;
	  if(fabs(tdPhi)<fXECut)NearXE=1;
	  else NearXE=0;
	  if(fabs(tdPhi)<(fPi/2))tdEta=tEta[iTrack]-tEta[iTrack2];
	  else tdEta=tEta[iTrack]+tEta[iTrack2];
	  if(tdPhi<fdPhiMin)tdPhi+=2*fPi;
	  if(tdPhi>fdPhiMax)tdPhi-=2*fPi;
	  if((tCharge[iTrack]<0&&tCharge[iTrack2]<0)||(tCharge[iTrack]>0&&tCharge[iTrack2]>0))Sign=1;
	  else Sign=2;
	  if(DEBUG) Printf("dPhi %f  dEta %f",tdPhi,tdEta);
	  for(int c=0;c<MaxArray;c++){//loop over multiplicity bins
	    fHistPtTrig[i][MultArray[c]][ievent]->Fill(tPt[iTrack2],tEff[iTrack2]);
	    fHistPhiTrig[i][MultArray[c]][ievent]->Fill(tPhi[iTrack2],tPt[iTrack2],tEff[iTrack2]);
	    fHistPhiTrigPt[i][MultArray[c]][ievent]->Fill(tPhi[iTrack2],tPt[iTrack2],tPt[iTrack2]*tEff[iTrack2]);
	    fHistEtaTrig[i][MultArray[c]][ievent]->Fill(tEta[iTrack2],tPt[iTrack2],tEff[iTrack2]);
	    fHistEtaTrigPt[i][MultArray[c]][ievent]->Fill(tEta[iTrack2],tPt[iTrack2],tPt[iTrack2]*tEff[iTrack2]);

	    fHistPhiEtaTrig[i][MultArray[c]][ievent]->Fill(tPhi[iTrack2],tEta[iTrack2],tPt[iTrack2],tEff[iTrack2]);
	    fHistPhiEtaTrigPt[i][MultArray[c]][ievent]->Fill(tPhi[iTrack2],tEta[iTrack2],tPt[iTrack2],tPt[iTrack2]*tEff[iTrack2]);
	    fHistDeltaPhi[i][MultArray[c]][0][ievent]->Fill(tdPhi,tPt[iTrack2],tEff[iTrack2]);
	    fHistDeltaPhiPt[i][MultArray[c]][0][ievent]->Fill(tdPhi,tPt[iTrack2],tPt[iTrack2]*tEff[iTrack2]);
	    fHistDeltaPhi[i][MultArray[c]][Sign][ievent]->Fill(tdPhi,tPt[iTrack2],tEff[iTrack2]);
	    fHistDeltaPhiPt[i][MultArray[c]][Sign][ievent]->Fill(tdPhi,tPt[iTrack2],tPt[iTrack2]*tEff[iTrack2]);

	    if(NearEta){
	    fHistDeltaEtaN[i][MultArray[c]][0][ievent]->Fill(tdEta,tPt[iTrack2],tEff[iTrack2]);
	    fHistDeltaEtaNPt[i][MultArray[c]][0][ievent]->Fill(tdEta,tPt[iTrack2],tPt[iTrack2]*tEff[iTrack2]);
	    fHistDeltaEtaN[i][MultArray[c]][Sign][ievent]->Fill(tdEta,tPt[iTrack2],tEff[iTrack2]);
	    fHistDeltaEtaNPt[i][MultArray[c]][Sign][ievent]->Fill(tdEta,tPt[iTrack2],tPt[iTrack2]*tEff[iTrack2]);
	    }
	    else{
	    fHistDeltaEtaA[i][MultArray[c]][0][ievent]->Fill(tdEta,tPt[iTrack2],tEff[iTrack2]);
	    fHistDeltaEtaAPt[i][MultArray[c]][0][ievent]->Fill(tdEta,tPt[iTrack2],tPt[iTrack2]*tEff[iTrack2]);
	    fHistDeltaEtaA[i][MultArray[c]][Sign][ievent]->Fill(tdEta,tPt[iTrack2],tEff[iTrack2]);
	    fHistDeltaEtaAPt[i][MultArray[c]][Sign][ievent]->Fill(tdEta,tPt[iTrack2],tPt[iTrack2]*tEff[iTrack2]);
	    }
	    fHistDeltaPhiEta[i][MultArray[c]][ievent]->Fill(tdPhi,tdEta,tPt[iTrack2],tEff[iTrack2]);
	    fHistDeltaPhiEtaPt[i][MultArray[c]][ievent]->Fill(tdPhi,tdEta,tPt[iTrack2],tPt[iTrack2]*tEff[iTrack2]);
	    
	    //only fill these if trigger particle is the leading particle
	    if(iTrack==LeadPart){
	      if(NearXE){
		tXE=tPt[iTrack2]*cos(tdPhi)/tPt[iTrack];
		fHistXEN[i][MultArray[c]][ievent]->Fill(tXE,tEff[iTrack2]);
	      }
	      else{
		tXE=tPt[iTrack2]*cos(tdPhi+fPi)/tPt[iTrack];
		fHistXEA[i][MultArray[c]][ievent]->Fill(tXE,tEff[iTrack2]);
	      }
	    }

	  }//Centrality loop (c)

	  //3-particle Correlations
	  for(int iTrack3=0;iTrack3<nGoodTracks[ievent];iTrack3++){
	    if(iTrack2==iTrack3)continue;
	    if(tPt[iTrack3]>tPt[iTrack])continue;
	    tdPhi2=tPhi[iTrack]-tPhi[iTrack3];
	    if(tdPhi2<-fPi)tdPhi2+=2*fPi;
	    if(tdPhi2>fPi)tdPhi2-=2*fPi;
	    if(fabs(tdPhi2)<fNearPhiCut&&NearEta==1)NearEta2=1;
	    else NearEta2=0;
	    //if(fabs(tdPhi)<fXECut)NearXE=1;
	    //else NearXE=0;
	    if(fabs(tdPhi2)<(fPi/2))tdEta2=tEta[iTrack]-tEta[iTrack3];
	    else tdEta2=tEta[iTrack]+tEta[iTrack3];
	    if(tdPhi2<fdPhiMin)tdPhi2+=2*fPi;
	    if(tdPhi2>fdPhiMax)tdPhi2-=2*fPi;
	    // if((tCharge[iTrack]<0&&tCharge[iTrack2]<0)||(tCharge[iTrack]>0&&tCharge[iTrack2]>0))Sign=1;
	    if((tCharge[iTrack]<0&&tCharge[iTrack2]<0&&tCharge[iTrack3]<0)||(tCharge[iTrack]>0&&tCharge[iTrack2]>0&&tCharge[iTrack3]>0))Sign=1;
	    else if((tCharge[iTrack3]<0&&tCharge[iTrack2]<0)||(tCharge[iTrack3]>0&&tCharge[iTrack2]>0))Sign=2;
	    else Sign=3;
	    for(int e=0;e<tNPtAssoc3[iTrack2];e++){//check associated pT bin
	      for(int f=0;f<tNPtAssoc3[iTrack3];f++){
		if(tPtAssoc3[iTrack2][e]==tPtAssoc3[iTrack3][f]){
		  for(int c=0;c<MaxArray;c++){//loop over multiplicity bins
		    fHistDeltaPhiPhi[i][tPtAssoc3[iTrack2][e]][MultArray[c]][0][ievent]->Fill(tdPhi,tdPhi2,tEff[iTrack2]*tEff[iTrack3]);
		    fHistDeltaPhiPhi[i][tPtAssoc3[iTrack2][e]][MultArray[c]][Sign][ievent]->Fill(tdPhi2,tdPhi,tEff[iTrack2]*tEff[iTrack3]);
		 

		    if(NearEta2){
		      fHistDeltaEtaEta[i][tPtAssoc3[iTrack2][e]][MultArray[c]][0][ievent]->Fill(tdEta,tdEta2,tEff[iTrack2]*tEff[iTrack3]);
		      fHistDeltaEtaEta[i][tPtAssoc3[iTrack2][e]][MultArray[c]][Sign][ievent]->Fill(tdEta,tdEta2,tEff[iTrack2]*tEff[iTrack3]);
		    }
		  }//multiplicity loop (c)
		}
	      }
	    }//track checking loops
	  }//iTrack3
	}//iTrack2 (associated track loop)
	
	if(DEBUG)Printf("Mixed Event Loop");
	for(int c=0;c<MaxArray;c++){
	  int d=MultArray[c];//Centrality bin we are in
	  if(fMixEnd[d][VertexBin][ievent]>=0){//check if there are any mixed events for this bin
	    for(int imix=0;imix<=fMixEnd[d][VertexBin][ievent];imix++){//loop over the stored mixed events
	      fHistNMix[d][ievent]->Fill(i);
	      for(int iTrack2=0;iTrack2<fMixTrack[imix][d][VertexBin][ievent];iTrack2++){
		if(tPt[iTrack]<fMPt[imix][d][VertexBin][ievent][iTrack2])continue;
		tdPhi=tPhi[iTrack]-fMPhi[imix][d][VertexBin][ievent][iTrack2];
		if(tdPhi<-fPi)tdPhi+=2*fPi;
		if(tdPhi>fPi)tdPhi-=2*fPi;
		if(fabs(tdPhi)<fNearPhiCut)NearEta=1;
		else NearEta=0;
		if(fabs(tdPhi)<fXECut)NearXE=1;
		else NearXE=0;
		if(fabs(tdPhi)<(fPi/2))tdEta=tEta[iTrack]-fMEta[imix][d][VertexBin][ievent][iTrack2];
		else tdEta=tEta[iTrack]+fMEta[imix][d][VertexBin][ievent][iTrack2];
		if(tdPhi<fdPhiMin)tdPhi+=2*fPi;	
		if(tdPhi>fdPhiMax)tdPhi-=2*fPi;
		if((tCharge[iTrack]<0&&fMCharge[imix][d][VertexBin][ievent][iTrack2]<0)||(tCharge[iTrack]>0&&fMCharge[imix][d][VertexBin][ievent][iTrack2]>0))Sign=1;
		else Sign=2;

		fHistDeltaPhiMix[i][d][0][ievent]->Fill(tdPhi,fMPt[imix][d][VertexBin][ievent][iTrack2],fMEff[imix][d][VertexBin][ievent][iTrack2]);
		fHistDeltaPhiMixPt[i][d][0][ievent]->Fill(tdPhi,fMPt[imix][d][VertexBin][ievent][iTrack2],fMPt[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix][d][VertexBin][ievent][iTrack2]);
		fHistDeltaPhiMix[i][d][Sign][ievent]->Fill(tdPhi,fMPt[imix][d][VertexBin][ievent][iTrack2],fMEff[imix][d][VertexBin][ievent][iTrack2]);
		fHistDeltaPhiMixPt[i][d][Sign][ievent]->Fill(tdEta,fMPt[imix][d][VertexBin][ievent][iTrack2],fMPt[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix][d][VertexBin][ievent][iTrack2]);
		if(NearEta){
		  fHistDeltaEtaNMix[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][VertexBin][ievent][iTrack2],fMEff[imix][d][VertexBin][ievent][iTrack2]);
		fHistDeltaEtaNMixPt[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][VertexBin][ievent][iTrack2],fMPt[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix][d][VertexBin][ievent][iTrack2]);
		fHistDeltaEtaNMix[i][d][Sign][ievent]->Fill(tdEta,fMPt[imix][d][VertexBin][ievent][iTrack2],fMEff[imix][d][VertexBin][ievent][iTrack2]);
		fHistDeltaEtaNMixPt[i][d][Sign][ievent]->Fill(tdEta,fMPt[imix][d][VertexBin][ievent][iTrack2],fMPt[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix][d][VertexBin][ievent][iTrack2]);
		}
		else{
		  fHistDeltaEtaAMix[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][VertexBin][ievent][iTrack2],fMEff[imix][d][VertexBin][ievent][iTrack2]);
		fHistDeltaEtaAMixPt[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][VertexBin][ievent][iTrack2],fMPt[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix][d][VertexBin][ievent][iTrack2]);
		fHistDeltaEtaAMix[i][d][Sign][ievent]->Fill(tdEta,fMPt[imix][d][VertexBin][ievent][iTrack2],fMEff[imix][d][VertexBin][ievent][iTrack2]);
		fHistDeltaEtaAMixPt[i][d][Sign][ievent]->Fill(tdEta,fMPt[imix][d][VertexBin][ievent][iTrack2],fMPt[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix][d][VertexBin][ievent][iTrack2]);
		}	

		fHistDeltaPhiEtaMix[i][d][ievent]->Fill(tdPhi,tdEta,fMPt[imix][d][VertexBin][ievent][iTrack2],fMEff[imix][d][VertexBin][ievent][iTrack2]);
		fHistDeltaPhiEtaMixPt[i][d][ievent]->Fill(tdPhi,tdEta,fMPt[imix][d][VertexBin][ievent][iTrack2],fMPt[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix][d][VertexBin][ievent][iTrack2]);

		if(iTrack==LeadPart){
		  if(NearXE){
		    tXE=fMPt[imix][d][VertexBin][ievent][iTrack2]*cos(tdPhi)/tPt[iTrack];
		    fHistXENMix[i][d][ievent]->Fill(tXE,fMEff[imix][d][VertexBin][ievent][iTrack2]);
		  }
		  else{
		    tXE=fMPt[imix][d][VertexBin][ievent][iTrack2]*cos(tdPhi+fPi)/tPt[iTrack];
		    fHistXEAMix[i][MultArray[c]][ievent]->Fill(tXE,fMEff[imix][d][VertexBin][ievent][iTrack2]);
		  }
		}
		//3-particle correlation soft-soft term (both associated from the same event)
		for(int iTrack3=0;iTrack3<fMixTrack[imix][d][VertexBin][ievent];iTrack3++){
		  if(iTrack3==iTrack2)continue;
		  if(tPt[iTrack]<fMPt[imix][d][VertexBin][ievent][iTrack3])continue;
		tdPhi2=tPhi[iTrack]-fMPhi[imix][d][VertexBin][ievent][iTrack3];
		if(tdPhi2<-fPi)tdPhi2+=2*fPi;
		if(tdPhi2>fPi)tdPhi2-=2*fPi;
		if(fabs(tdPhi2)<fNearPhiCut&&NearEta)NearEta2=1;
		else NearEta2=0;
		if(fabs(tdPhi2)<(fPi/2))tdEta2=tEta[iTrack]-fMEta[imix][d][VertexBin][ievent][iTrack3];
		else tdEta2=tEta[iTrack]+fMEta[imix][d][VertexBin][ievent][iTrack3];
		if(tdPhi2<fdPhiMin)tdPhi2+=2*fPi;	
		if(tdPhi2>fdPhiMax)tdPhi2-=2*fPi;
		//if((tCharge[iTrack]<0&&fMCharge[imix][d][VertexBin][ievent][iTrack2]<0)||(tCharge[iTrack]>0&&fMCharge[imix][d][VertexBin][ievent][iTrack2]>0))Sign=1;
		//else Sign=2;
		if((tCharge[iTrack]<0&&fMCharge[imix][d][VertexBin][ievent][iTrack2]<0&&fMCharge[imix][d][VertexBin][ievent][iTrack3]<0)||(tCharge[iTrack]>0&&fMCharge[imix][d][VertexBin][ievent][iTrack2]>0&&fMCharge[imix][d][VertexBin][ievent][iTrack3]>0))Sign=1;
		else if((fMCharge[imix][d][VertexBin][ievent][iTrack3]<0&&fMCharge[imix][d][VertexBin][ievent][iTrack2]<0)||(fMCharge[imix][d][VertexBin][ievent][iTrack3]>0&&fMCharge[imix][d][VertexBin][ievent][iTrack2]>0))Sign=2;
		else Sign=3;
		for(int e=0;e<fMNPtAssoc3[imix][d][VertexBin][ievent][iTrack2];e++){//check associated pT bin
		  for(int f=0;f<fMNPtAssoc3[imix][d][VertexBin][ievent][iTrack3];f++){
		    if(fMPtAssoc3[imix][d][VertexBin][ievent][e][iTrack2]==fMPtAssoc3[imix][d][VertexBin][ievent][f][iTrack3]){
		      fHistDeltaPhiPhiSS[i][fMPtAssoc3[imix][d][VertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix][d][VertexBin][ievent][iTrack3]);
		      fHistDeltaPhiPhiSS[i][fMPtAssoc3[imix][d][VertexBin][ievent][e][iTrack2]][d][Sign][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix][d][VertexBin][ievent][iTrack3]); 

		      if(NearEta2){
			fHistDeltaEtaEtaSS[i][fMPtAssoc3[imix][d][VertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix][d][VertexBin][ievent][iTrack3]);
			fHistDeltaEtaEtaSS[i][fMPtAssoc3[imix][d][VertexBin][ievent][e][iTrack2]][d][Sign][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix][d][VertexBin][ievent][iTrack3]);
		      }//near-side
		    }
		  }
		}//associated pt bin
		}//iTrack3

		//3-particle mixed event (associated from different events)
		//for(int imix2=0;imix2<=fMixEnd[d][VertexBin][ievent];imix2++){//loop over the stored mixed events
		//if(imix2==imix)continue;
		  int imix2=imix+1;
		  if(imix2>=fMixEnd[d][VertexBin][ievent])imix2=0;
		  if(imix2==imix)continue;//will kill it when there is only 1 mixed event (remember to scale by 1 less then the number of mixed events in others (number of mixed-mixed is 2*(mixed-1)))
		  for(int iTrack3=0;iTrack3<fMixTrack[imix2][d][VertexBin][ievent];iTrack3++){
		    if(tPt[iTrack]<fMPt[imix2][d][VertexBin][ievent][iTrack3])continue;
		    tdPhi2=tPhi[iTrack]-fMPhi[imix2][d][VertexBin][ievent][iTrack3];
		    if(tdPhi2<-fPi)tdPhi2+=2*fPi;
		    if(tdPhi2>fPi)tdPhi2-=2*fPi;
		    if(fabs(tdPhi2)<fNearPhiCut&&NearEta)NearEta2=1;
		    else NearEta2=0;
		    if(fabs(tdPhi2)<(fPi/2))tdEta2=tEta[iTrack]-fMEta[imix2][d][VertexBin][ievent][iTrack3];
		    else tdEta2=tEta[iTrack]+fMEta[imix2][d][VertexBin][ievent][iTrack3];
		    if(tdPhi2<fdPhiMin)tdPhi2+=2*fPi;	
		    if(tdPhi2>fdPhiMax)tdPhi2-=2*fPi;
		    //if((tCharge[iTrack]<0&&fMCharge[imix2][d][VertexBin][ievent][iTrack2]<0)||(tCharge[iTrack]>0&&fMCharge[imix2][d][VertexBin][ievent][iTrack2]>0))Sign=1;
		    //else Sign=2;
		    if((tCharge[iTrack]<0&&fMCharge[imix][d][VertexBin][ievent][iTrack2]<0&&fMCharge[imix2][d][VertexBin][ievent][iTrack3]<0)||(tCharge[iTrack]>0&&fMCharge[imix][d][VertexBin][ievent][iTrack2]>0&&fMCharge[imix2][d][VertexBin][ievent][iTrack3]>0))Sign=1;
		    else if((fMCharge[imix2][d][VertexBin][ievent][iTrack3]<0&&fMCharge[imix][d][VertexBin][ievent][iTrack2]<0)||(fMCharge[imix2][d][VertexBin][ievent][iTrack3]>0&&fMCharge[imix][d][VertexBin][ievent][iTrack2]>0))Sign=2;
		    else Sign=3;
		    for(int e=0;e<fMNPtAssoc3[imix][d][VertexBin][ievent][iTrack2];e++){//check associated pT bin
		      for(int f=0;f<fMNPtAssoc3[imix2][d][VertexBin][ievent][iTrack3];f++){
			if(fMPtAssoc3[imix][d][VertexBin][ievent][e][iTrack2]==fMPtAssoc3[imix2][d][VertexBin][ievent][f][iTrack3]){
			  fHistDeltaPhiPhiMix[i][fMPtAssoc3[imix][d][VertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix2][d][VertexBin][ievent][iTrack3]);
			  fHistDeltaPhiPhiMix[i][fMPtAssoc3[imix][d][VertexBin][ievent][e][iTrack2]][d][Sign][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix2][d][VertexBin][ievent][iTrack3]); 
			  fHistDeltaPhiPhiMix[i][fMPtAssoc3[imix][d][VertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix2][d][VertexBin][ievent][iTrack3]);//free factor of 2 in statistics
			  fHistDeltaPhiPhiMix[i][fMPtAssoc3[imix][d][VertexBin][ievent][e][iTrack2]][d][Sign][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix2][d][VertexBin][ievent][iTrack3]); 
		      if(NearEta2){
			fHistDeltaEtaEtaMix[i][fMPtAssoc3[imix][d][VertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix2][d][VertexBin][ievent][iTrack3]);
			fHistDeltaEtaEtaMix[i][fMPtAssoc3[imix][d][VertexBin][ievent][e][iTrack2]][d][Sign][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix2][d][VertexBin][ievent][iTrack3]);
			fHistDeltaEtaEtaMix[i][fMPtAssoc3[imix][d][VertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta2,tdEta,fMEff[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix2][d][VertexBin][ievent][iTrack3]);
			fHistDeltaEtaEtaMix[i][fMPtAssoc3[imix][d][VertexBin][ievent][e][iTrack2]][d][Sign][ievent]->Fill(tdEta2,tdEta,fMEff[imix][d][VertexBin][ievent][iTrack2]*fMEff[imix2][d][VertexBin][ievent][iTrack3]);
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
    if(DEBUG)Printf("Store Event For Mixing");
  for(int c=0;c<MaxArray;c++){//loops over centrality bins
    int d=MultArray[c];//too many nested arrays looked confusing d=which centrality bin
    if(fMixEnd[d][VertexBin][ievent]<(fNMix-1))fMixEnd[d][VertexBin][ievent]++;
    if(fMixPointer[d][VertexBin][ievent]<(fNMix-1))fMixPointer[d][VertexBin][ievent]++;
    else fMixPointer[d][VertexBin][ievent]=0;
    int e=fMixPointer[d][VertexBin][ievent];//nested arrays (e is event number in pool)
    delete [] fMPt[e][d][VertexBin][ievent];
    delete [] fMPhi[e][d][VertexBin][ievent];
    delete [] fMEta[e][d][VertexBin][ievent];
    delete [] fMCharge[e][d][VertexBin][ievent];
    delete [] fMEff[e][d][VertexBin][ievent];
    delete [] fMNPtAssoc3[e][d][VertexBin][ievent];
    for(int jj=0;jj<10;jj++){
    delete [] fMPtAssoc3[e][d][VertexBin][ievent][jj];
    }
    fMPt[e][d][VertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
    fMPhi[e][d][VertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
    fMEta[e][d][VertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
    fMCharge[e][d][VertexBin][ievent]=new Short_t [nGoodTracks[ievent]];
    fMEff[e][d][VertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
    fMixTrack[e][d][VertexBin][ievent]=nGoodTracks[ievent];
    fMNPtAssoc3[e][d][VertexBin][ievent]=new Short_t [nGoodTracks[ievent]];
    for(int jj=0;jj<10;jj++){
    fMPtAssoc3[e][d][VertexBin][ievent][jj]=new Short_t [nGoodTracks[ievent]];
    }

    for(int iTrack=0;iTrack<nGoodTracks[ievent];iTrack++){
      fMPt[e][d][VertexBin][ievent][iTrack]=tPt[iTrack];
      fMPhi[e][d][VertexBin][ievent][iTrack]=tPhi[iTrack];
      fMEta[e][d][VertexBin][ievent][iTrack]=tEta[iTrack];
      fMCharge[e][d][VertexBin][ievent][iTrack]=tCharge[iTrack];
      fMEff[e][d][VertexBin][ievent][iTrack]=tEff[iTrack];
      fMNPtAssoc3[e][d][VertexBin][ievent][iTrack]=tNPtAssoc3[iTrack];
      // fMPtAssoc3[e][d][VertexBin][ievent][iTrack]=new Int_t [tNPtAssoc3[iTrack]];
      for(int jj=0;jj<tNPtAssoc3[iTrack];jj++){
	//if(DEBUG) Printf("%d",tPtAssoc3[iTrack][jj]);
	fMPtAssoc3[e][d][VertexBin][ievent][jj][iTrack]=tPtAssoc3[iTrack][jj];
      }
      
    }//iTracks
  }//Centrality (c)
  }//ievent
  //track=0;
  //track2=0;

  PostData(0, fOutput);
  //get rid of these arrays from memory
  delete [] tPhi;
  delete [] tEta;
  delete [] tPt;
  delete [] tCharge;
  delete [] tEff;
  delete [] tNPtAssoc3;
  delete [] tPtAssoc3;
  tPhi=NULL;
  tEta=NULL;
  tPt=NULL;
  tCharge=NULL;
  tEff=NULL;
  tNPtAssoc3=NULL;
  tPtAssoc3=NULL;

}//Exec

//---------------------------------------------------
void AliAnalysisTaskDiHadron::Terminate(Option_t *){
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
  Printf("Terminate AliAnalysisTaskDiHadron");
}
