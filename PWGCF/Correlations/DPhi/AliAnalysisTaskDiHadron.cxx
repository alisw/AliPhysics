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
#include "TROOT.h"
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
#include "TRandom3.h"
#include "TSystem.h"


#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliESDVZERO.h"
#include "TParticle.h"
#include "AliCentrality.h"

//#include "AliHeader.h"
//#include "AliGenEventHeader.h"



#include "AliAnalysisTaskDiHadron.h"


ClassImp(AliAnalysisTaskDiHadron)

//----------------------------------------
AliAnalysisTaskDiHadron::AliAnalysisTaskDiHadron(const char *name):
AliAnalysisTaskSE(name), fESD(0), fAOD(0), fMC(0), fOutput(0),fMinClustersTPC(0),fMinClusterRatio(0),fMaxTPCchi2(0),fMinClustersITS(0),fEtaCut(0),fTrigEtaCut(0),fNearPhiCut(0),fXECut(0),fMaxDCA(0),fMaxDCAXY(0),fMaxDCAZ(0),fDCA2D(0),fTPCRefit(0),fITSRefit(0),fSPDCut(0),fMinPtAssoc(0),fMaxPtAssoc(0),fVzCut(0),fAODData(0),fEfficiencyCorr(0),fDEBUG(0),fnBinPhi(0),fnBinEta(0),fnBinPhiEtaPhi(0),fnBinPhiEtaEta(0),fnBinPhi3(0),fnBinEta3(0),fPi(3.1415926535898),fdPhiMin(0),fdPhiMax(0),fNTPtBins(0),fNMix(0),fNCentBins(0),fCentPercent(0),fNAPtBins(0),fNAPt3Bins(0),fNVertexBins(0),fNXEBins(0),fNIDs(0),fEffFitPt(0),fNFitLowParam(0),fNFitHighParam(0),fV2FitPt(0),fV3FitPt(0),fV4FitPt(0),fNFitLowParamV2(0),fNFitHighParamV2(0),fNFitLowParamV3(0),fNFitHighParamV3(0),fNFitLowParamV4(0),fNFitHighParamV4(0),fMCHistos(0),fFitLow(NULL),fFitHigh(NULL),fFitLowParam(NULL),fFitHighParam(NULL),fFitLowV2(NULL),fFitHighV2(NULL),fFitLowParamV2(NULL),fFitHighParamV2(NULL),fFitLowV3(NULL),fFitHighV3(NULL),fFitLowParamV3(NULL),fFitHighParamV3(NULL),fFitLowV4(NULL),fFitHighV4(NULL),fFitLowParamV4(NULL),fFitHighParamV4(NULL),fPtTrigArray(NULL),fPtAssocArray(NULL),fPtAssoc3Array1(NULL),fPtAssoc3Array2(NULL),fCentArrayMin(NULL),fCentArrayMax(NULL),fXEArray(NULL),fTrigIDArray(NULL),fSimulate(0),fSimNBgPart(0),fSimNJetPart(0),fSimNJet(0),fSimNEvents(0),ftPhi(NULL),ftEta(NULL),ftPt(NULL),ftCharge(NULL),ftEff(NULL),ftV2(NULL),ftV3(NULL),ftV4(NULL),ftPtAssoc3(NULL),ftNPtAssoc3(NULL)
  
  {
    //TRandom *gRandom=new TRandom3();
  
  //IO Slots
  DefineInput(0, TChain::Class());
  DefineOutput(0,TList::Class());
  
  
  for(int c=0;c<fNCentBins;c++){
    for(int v=0;v<fNVertexBins;v++){
      for(int jmc=0;jmc<2;jmc++){
	fMixPointer[c][v][jmc]=0;
	fMixEnd[c][v][jmc]=0;
	for(int ievts=0;ievts<fNMix;ievts++){
	  fMPt[ievts][c][v][jmc]=NULL;
	  fMPhi[ievts][c][v][jmc]=NULL;
	  fMEta[ievts][c][v][jmc]=NULL;
	  fMCharge[ievts][c][v][jmc]=NULL;
	  fMEff[ievts][c][v][jmc]=NULL;
	  fMV2[ievts][c][v][jmc]=NULL;
	  fMV4[ievts][c][v][jmc]=NULL;
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
  if(fDEBUG)Printf("Setting Cuts");
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
void AliAnalysisTaskDiHadron::SetOptions(Int_t AODData, Int_t EfficiencyCorr, Int_t ffDEBUG,  Int_t MCHistos){
//Sets some options
  if(fDEBUG) Printf("Setting Options");
  fAODData=AODData;
  fEfficiencyCorr=EfficiencyCorr;
  fDEBUG=ffDEBUG;
  fMCHistos=MCHistos;
}
//------------------------------------------------------
void AliAnalysisTaskDiHadron::SetBins(Int_t nBinPhi, Int_t nBinEta, Int_t nBinPhiEtaPhi, Int_t nBinPhiEtaEta, Int_t nBinPhi3, Int_t nBinEta3,Float_t dPhiMin, Float_t dPhiMax, Int_t NTPtBins, Int_t NMixBins, Int_t NCentBins, Int_t CentPercent, Int_t NAPtBins, Int_t NAPt3Bins, Int_t NVertexBins, Int_t NXEBins,Float_t *PtTrigArray, Float_t *PtAssocArray,Float_t *PtAssoc3Array1, Float_t *PtAssoc3Array2, Int_t *CentArrayMin, Int_t *CentArrayMax, Float_t *XEArray){
//sets up the histogram binning
  if(fDEBUG)Printf("Setting Binning");
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
  fCentPercent=CentPercent;
  fNAPtBins=NAPtBins;
  fNAPt3Bins=NAPt3Bins;
  fNVertexBins=NVertexBins;
  fNXEBins=NXEBins;
  fPtTrigArray=new Float_t [fNTPtBins];
  for(int i=0;i<fNTPtBins;i++)fPtTrigArray[i]=PtTrigArray[i];
  fPtAssocArray=new Float_t [fNAPtBins];
  for(int i=0;i<=fNAPtBins;i++)fPtAssocArray[i]=PtAssocArray[i];
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
  if(fDEBUG)Printf("Setting Efficiencies");
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
void AliAnalysisTaskDiHadron::SetFlow(Float_t V2FitPt, Float_t V3FitPt, Float_t V4FitPt, const TF1 *FitLowV2, const TF1 *FitHighV2, const TF1 *FitLowV3, const TF1 *FitHighV3, const TF1 *FitLowV4, const TF1 *FitHighV4, Int_t NFitLowParamV2, Int_t NFitHighParamV2, Int_t NFitLowParamV3, Int_t NFitHighParamV3, Int_t NFitLowParamV4, Int_t NFitHighParamV4, Float_t *FitLowParamV2, Float_t *FitHighParamV2, Float_t *FitLowParamV3, Float_t *FitHighParamV3, Float_t *FitLowParamV4, Float_t *FitHighParamV4){
  if(fDEBUG)Printf("Setting Flow");
  fV2FitPt=V2FitPt;
  fV3FitPt=V3FitPt;
  fV4FitPt=V4FitPt;
  fFitLowV2=(TF1*)FitLowV2;
  fFitHighV2=(TF1*)FitHighV2;
  fFitLowV3=(TF1*)FitLowV3;
  fFitHighV3=(TF1*)FitHighV3;
  fFitLowV4=(TF1*)FitLowV4;
  fFitHighV4=(TF1*)FitHighV4;
  fNFitLowParamV2=NFitLowParamV2;
  fNFitHighParamV2=NFitHighParamV2;
  fNFitLowParamV3=NFitLowParamV3;
  fNFitHighParamV3=NFitHighParamV3;
  fNFitLowParamV4=NFitLowParamV4;
  fNFitHighParamV4=NFitHighParamV4;

  fFitLowParamV2=new Float_t [fNFitLowParamV2*fNCentBins];
  for(int i=0;i<fNFitLowParamV2*fNCentBins;i++)fFitLowParamV2[i]=FitLowParamV2[i];
  fFitHighParamV2=new Float_t [fNFitHighParamV2*fNCentBins];
  for(int i=0;i<fNFitHighParamV2*fNCentBins;i++)fFitHighParamV2[i]=FitHighParamV2[i];

  fFitLowParamV3=new Float_t [fNFitLowParamV3*fNCentBins];
  for(int i=0;i<fNFitLowParamV3*fNCentBins;i++)fFitLowParamV3[i]=FitLowParamV3[i];
  fFitHighParamV3=new Float_t [fNFitHighParamV3*fNCentBins];
  for(int i=0;i<fNFitHighParamV3*fNCentBins;i++)fFitHighParamV3[i]=FitHighParamV3[i];

  fFitLowParamV4=new Float_t [fNFitLowParamV4*fNCentBins];
  for(int i=0;i<fNFitLowParamV4*fNCentBins;i++)fFitLowParamV4[i]=FitLowParamV4[i];
  fFitHighParamV4=new Float_t [fNFitHighParamV4*fNCentBins];
  for(int i=0;i<fNFitHighParamV4*fNCentBins;i++)fFitHighParamV4[i]=FitHighParamV4[i];
  if(fDEBUG)Printf("FlowSet");
  }
//----------------------------------------------------------
void AliAnalysisTaskDiHadron::SetSimulation(Int_t Simulate, Float_t SimNBgPart, Float_t SimNJetPart, Float_t SimNJet, Int_t SimNEvents){
  fSimulate=Simulate;
  fSimNBgPart=SimNBgPart;
  fSimNJetPart=SimNJetPart;
  fSimNJet=SimNJet;
  fSimNEvents=SimNEvents;

}
  
//-----------------------------------------------------------
void AliAnalysisTaskDiHadron::ConnectInputData(Option_t *){
  //Connect to ESD
  if(fDEBUG)Printf("Connecting");
   TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
   if (!tree&&fDEBUG) {Printf("ERROR: Could not read chain from input slot 0");} 
  else {
    if(!fAODData){
AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!esdH){if(fDEBUG) Printf("ERROR: Could not get ESDInputHandler");}
    else fESD = esdH->GetEvent();
    }
    else{
        AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
     if (!aodH) {
       Printf("ERROR: Could not get AODInputHandler");
     }
     else{
       fAOD = aodH->GetEvent();
     }
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
  char histname[300];
  char histtitle[300];
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
  Int_t BufferSize=256;
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
    snprintf(histname,BufferSize,"fHistMult%s",cmc1[imc]);
    snprintf(histtitle,BufferSize,"Multiplicity%s",cmc2[imc]);
    fHistMult[imc]=new TH1F(histname,histtitle,2000,-0.5,1999.5);
    fHistMult[imc]->Sumw2();
    fHistMult[imc]->GetXaxis()->SetTitle("Number of tracks");
    fHistMult[imc]->GetYaxis()->SetTitle("Counts");
    fOutput->Add(fHistMult[imc]);
    
    for(int imult=0;imult<fNCentBins;imult++){//loop for multiplicity bins
      
      //Histograms that are independent of the trigger
      snprintf(histname,BufferSize,"fHistPt_C%d%s",imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"P_{T} Distribution of Tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPt[imult][imc]=new TH1F(histname,histtitle,nptbins,fPtAssocArray);
      fHistPt[imult][imc]->Sumw2();
      fHistPt[imult][imc]->GetXaxis()->SetTitle("p_{T}");
      fHistPt[imult][imc]->GetYaxis()->SetTitle("Counts");
      fOutput->Add(fHistPt[imult][imc]);

 //Histograms that are independent of the trigger
      snprintf(histname,BufferSize,"fHistPtEff_C%d%s",imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"P_{T} Distribution of Tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPtEff[imult][imc]=new TH1F(histname,histtitle,1000,0,100);
      fHistPtEff[imult][imc]->Sumw2();
      fHistPtEff[imult][imc]->GetXaxis()->SetTitle("p_{T}");
      fHistPtEff[imult][imc]->GetYaxis()->SetTitle("Counts");
      fOutput->Add(fHistPtEff[imult][imc]);
      
      snprintf(histname,BufferSize,"fHistPhi_C%d%s",imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"#phi Distribution of Tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPhi[imult][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,nptbins,fPtAssocArray);
      fHistPhi[imult][imc]->Sumw2();
      fHistPhi[imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistPhi[imult][imc]->GetYaxis()->SetTitle("P_{T}");
      fOutput->Add(fHistPhi[imult][imc]);
      
      snprintf(histname,BufferSize,"fHistPhiPt_C%d%s",imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"P_{T} weighted #phi Distribution of tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPhiPt[imult][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,nptbins,fPtAssocArray);
      fHistPhiPt[imult][imc]->Sumw2();
      fHistPhiPt[imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistPhiPt[imult][imc]->GetYaxis()->SetTitle("P_{T}");
      fOutput->Add(fHistPhiPt[imult][imc]);
      
      snprintf(histname,BufferSize,"fHistEta_C%d%s",imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"#eta Distribution of Tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistEta[imult][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,nptbins,fPtAssocArray);
      fHistEta[imult][imc]->Sumw2();
      fHistEta[imult][imc]->GetXaxis()->SetTitle("#eta");
      fHistEta[imult][imc]->GetYaxis()->SetTitle("P_{T}");
      fOutput->Add(fHistEta[imult][imc]);
      
      snprintf(histname,BufferSize,"fHistEtaPt_C%d%s",imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"P_{T} weighted #eta Distribution of tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistEtaPt[imult][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,nptbins,fPtAssocArray);
      fHistEtaPt[imult][imc]->Sumw2();
      fHistEtaPt[imult][imc]->GetXaxis()->SetTitle("#eta");
      fHistEtaPt[imult][imc]->GetYaxis()->SetTitle("P_{T}");
      fOutput->Add(fHistEtaPt[imult][imc]);
      
      snprintf(histname,BufferSize,"fHistNEvents_C%d%s",imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"Number of Events and Number Passing Cuts %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistNEvents[imult][imc]=new TH1F(histname,histtitle,2,-0.5,1.5);
      fHistNEvents[imult][imc]->Sumw2();
      fHistNEvents[imult][imc]->GetXaxis()->SetTitle("Events,Passing Cuts");
      fHistNEvents[imult][imc]->GetYaxis()->SetTitle("Number of Events");
      fOutput->Add(fHistNEvents[imult][imc]);
      
      snprintf(histname,BufferSize,"fHistNTrigger_C%d%s",imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"Number of Triggers %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistNTrigger[imult][imc]=new TH1F(histname,histtitle,fNTPtBins,-0.5,fNTPtBins-0.5);
      fHistNTrigger[imult][imc]->Sumw2();
      fHistNTrigger[imult][imc]->GetXaxis()->SetTitle("Trigger Number");
      fHistNTrigger[imult][imc]->GetYaxis()->SetTitle("Number of Triggers");
      fOutput->Add(fHistNTrigger[imult][imc]);
      
      snprintf(histname,BufferSize,"fHistNTriggerPt_C%d%s",imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"P_{T} Weighted Number of Triggers %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistNTriggerPt[imult][imc]=new TH1F(histname,histtitle,fNTPtBins,-0.5,fNTPtBins-0.5);
      fHistNTriggerPt[imult][imc]->Sumw2();
      fHistNTriggerPt[imult][imc]->GetXaxis()->SetTitle("Trigger Number");
      fHistNTriggerPt[imult][imc]->GetYaxis()->SetTitle("Number of Triggers");
      fOutput->Add(fHistNTriggerPt[imult][imc]);
      
      snprintf(histname,BufferSize,"fHistNMix_C%d%s",imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"Number of Mixed Events %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistNMix[imult][imc]=new TH1F(histname,histtitle,fNTPtBins,-0.5,fNTPtBins-0.5);
      fHistNMix[imult][imc]->Sumw2();
      fHistNMix[imult][imc]->GetXaxis()->SetTitle("Trigger Number");
      fHistNMix[imult][imc]->GetYaxis()->SetTitle("Number of Mixed Events");
      fOutput->Add(fHistNMix[imult][imc]);
      
      snprintf(histname,BufferSize,"fHistPhiEta_C%d%s",imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"#phi-#eta distribution of tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPhiEta[imult][imc]=new TH3F(histname, histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,nptbins,fPtAssocArray);
      fHistPhiEta[imult][imc]->Sumw2();
      fHistPhiEta[imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistPhiEta[imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistPhiEta[imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistPhiEta[imult][imc]);
      
      snprintf(histname,BufferSize,"fHistPhiEtaPt_C%d%s",imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"Pt Weighted #phi-#eta distribution of tracks %dMult%d%s",fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
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
	  snprintf(histname,BufferSize,"fHistMultTrig_P%d%s",i,cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Distrubition of number of tracks in triggered events with %3.1f<p_{T}^{Trig}<%3.1f%s",fPtTrigArray[i],fPtTrigArray[i+1],cmc2[imc]);
	  fHistMultTrig[i][imc]=new TH1F(histname,histtitle,2000,0,2000);
	  fHistMultTrig[i][imc]->Sumw2();
	  fHistMultTrig[i][imc]->GetXaxis()->SetTitle("Number of Tracks");
	  fHistMultTrig[i][imc]->GetYaxis()->SetTitle("Counts");
	  fOutput->Add(fHistMultTrig[i][imc]);
	}
	snprintf(histname,BufferSize,"fHistPtTrig_P%d_C%d%s",i,imult,cmc1[imc]);
	snprintf(histtitle,BufferSize,"P_{T} distribution in triggered events with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistPtTrig[i][imult][imc]=new TH1F(histname,histtitle,nptbins,fPtAssocArray);
	fHistPtTrig[i][imult][imc]->Sumw2();
	fHistPtTrig[i][imult][imc]->GetXaxis()->SetTitle("p_{T}");
	fHistPtTrig[i][imult][imc]->GetYaxis()->SetTitle("Counts");
	fOutput->Add(fHistPtTrig[i][imult][imc]);
	
	snprintf(histname,BufferSize,"fHistPhiTrig_P%d_C%d%s",i,imult,cmc1[imc]);
	snprintf(histtitle,BufferSize,"Phi Distribution of triggered events with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistPhiTrig[i][imult][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	fHistPhiTrig[i][imult][imc]->Sumw2();
	fHistPhiTrig[i][imult][imc]->GetXaxis()->SetTitle("#phi");
	fHistPhiTrig[i][imult][imc]->GetYaxis()->SetTitle("p_{T}");
	fOutput->Add(fHistPhiTrig[i][imult][imc]);
	
	snprintf(histname,BufferSize,"fHistPhiTrigPt_P%d_C%d%s",i,imult,cmc1[imc]);
	snprintf(histtitle,BufferSize,"P_{T} Weighted Phi Distribution of triggered events with %3.1f<p_{T}^{Trig}<%3.1f%s",fPtTrigArray[i],fPtTrigArray[i+1],cmc2[imc]);
	fHistPhiTrigPt[i][imult][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	fHistPhiTrigPt[i][imult][imc]->Sumw2();
	fHistPhiTrigPt[i][imult][imc]->GetXaxis()->SetTitle("#phi");
	fHistPhiTrigPt[i][imult][imc]->GetYaxis()->SetTitle("p_{T}");
	fOutput->Add(fHistPhiTrigPt[i][imult][imc]);
	
	snprintf(histname,BufferSize,"fHistEtaTrig_P%d_C%d%s",i,imult,cmc1[imc]);
	snprintf(histtitle,BufferSize,"Eta Distribution of triggered events with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistEtaTrig[i][imult][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	fHistEtaTrig[i][imult][imc]->Sumw2();
	fHistEtaTrig[i][imult][imc]->GetXaxis()->SetTitle("#eta");
	fHistEtaTrig[i][imult][imc]->GetYaxis()->SetTitle("p_{T}");
	fOutput->Add(fHistEtaTrig[i][imult][imc]);
	
	snprintf(histname,BufferSize,"fHistEtaTrigPt_P%d_C%d%s",i,imult,cmc1[imc]);
	snprintf(histtitle,BufferSize,"P_{T} Weighted Eta Distribution of triggered events with %3.1f<p_{T}^{Trig}<%3.1f%s",fPtTrigArray[i],fPtTrigArray[i+1],cmc2[imc]);
	fHistEtaTrigPt[i][imult][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	fHistEtaTrigPt[i][imult][imc]->Sumw2();
	fHistEtaTrigPt[i][imult][imc]->GetXaxis()->SetTitle("#eta");
	fHistEtaTrigPt[i][imult][imc]->GetYaxis()->SetTitle("p_{T}");
	fOutput->Add(fHistEtaTrigPt[i][imult][imc]);
	
	snprintf(histname,BufferSize,"fHistPhiEtaTrig_P%d_C%d%s",i,imult,cmc1[imc]);
	snprintf(histtitle,BufferSize,"#phi-#eta distribution in triggered events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistPhiEtaTrig[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
	fHistPhiEtaTrig[i][imult][imc]->Sumw2();
	fHistPhiEtaTrig[i][imult][imc]->GetXaxis()->SetTitle("#phi");
	fHistPhiEtaTrig[i][imult][imc]->GetYaxis()->SetTitle("#eta");
	fHistPhiEtaTrig[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
	fOutput->Add(fHistPhiEtaTrig[i][imult][imc]);

	snprintf(histname,BufferSize,"fHistXEN_P%d_C%d%s",i,imult,cmc1[imc]);
	snprintf(histtitle,BufferSize,"Near-Side X_{E} distribution for %3.1f<p_{T}^{Lead}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistXEN[i][imult][imc]=new TH1F(histname,histtitle,fNXEBins,fXEArray);
	fHistXEN[i][imult][imc]->Sumw2();
	fHistXEN[i][imult][imc]->GetXaxis()->SetTitle("X_{E}");
	fOutput->Add(fHistXEN[i][imult][imc]);

	snprintf(histname,BufferSize,"fHistXENMixed_P%d_C%d%s",i,imult,cmc1[imc]);
	snprintf(histtitle,BufferSize,"Mixed Near-Side X_{E} distribution for %3.1f<p_{T}^{Lead}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistXENMix[i][imult][imc]=new TH1F(histname,histtitle,fNXEBins,fXEArray);
	fHistXENMix[i][imult][imc]->Sumw2();
	fHistXENMix[i][imult][imc]->GetXaxis()->SetTitle("X_{E}");
	fOutput->Add(fHistXENMix[i][imult][imc]);
	
	snprintf(histname,BufferSize,"fHistXEA_P%d_C%d%s",i,imult,cmc1[imc]);
	snprintf(histtitle,BufferSize,"Away-Side X_{E} distribution for %3.1f<p_{T}^{Lead}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistXEA[i][imult][imc]=new TH1F(histname,histtitle,fNXEBins,fXEArray);
	fHistXEA[i][imult][imc]->Sumw2();
	fHistXEA[i][imult][imc]->GetXaxis()->SetTitle("X_{E}");
	fOutput->Add(fHistXEA[i][imult][imc]);

	snprintf(histname,BufferSize,"fHistXEAMixed_P%d_C%d%s",i,imult,cmc1[imc]);
	snprintf(histtitle,BufferSize,"Mixed Away-Side X_{E} distribution for %3.1f<p_{T}^{Lead}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
	fHistXEAMix[i][imult][imc]=new TH1F(histname,histtitle,fNXEBins,fXEArray);
	fHistXEAMix[i][imult][imc]->Sumw2();
	fHistXEAMix[i][imult][imc]->GetXaxis()->SetTitle("X_{E}");
	fOutput->Add(fHistXEAMix[i][imult][imc]);

	//signloop
	for(int isign=0;isign<3;isign++){
	  snprintf(histname,BufferSize,"fHistDeltaPhi_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"#Delta#phi Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhi[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhi[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhi[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhi[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhi[i][imult][isign][imc]);
	  
	  snprintf(histname,BufferSize,"fHistDeltaPhiPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"P_{T} Weighted #Delta#phi Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiPt[i][imult][isign][imc]);
	  
	  snprintf(histname,BufferSize,"fHistDeltaPhiMix_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"#Delta#phi Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiMix[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiMix[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiMix[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiMix[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiMix[i][imult][isign][imc]);

	    snprintf(histname,BufferSize,"fHistDeltaPhiMixV2_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"#Delta#phi Mixed Event V2 Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiMixV2[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiMixV2[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiMixV2[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiMixV2[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiMixV2[i][imult][isign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaPhiMixV3_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"#Delta#phi Mixed Event V3 Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiMixV3[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiMixV3[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiMixV3[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiMixV3[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiMixV3[i][imult][isign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaPhiMixV4_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"#Delta#phi Mixed Event V4 Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiMixV4[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiMixV4[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiMixV4[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiMixV4[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiMixV4[i][imult][isign][imc]);
	  
	  snprintf(histname,BufferSize,"fHistDeltaPhiMixPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"P_{T} Weighted #Delta#phi Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiMixPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiMixPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiMixPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiMixPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiMixPt[i][imult][isign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaPhiMixV2Pt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"P_{T} Weighted #Delta#phi Mixed Event V2 Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiMixV2Pt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiMixV2Pt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiMixV2Pt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiMixV2Pt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiMixV2Pt[i][imult][isign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaPhiMixV3Pt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"P_{T} Weighted #Delta#phi Mixed Event V3 Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiMixV3Pt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiMixV3Pt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiMixV3Pt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiMixV3Pt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiMixV3Pt[i][imult][isign][imc]);
	  
	  snprintf(histname,BufferSize,"fHistDeltaPhiMixV4Pt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"P_{T} Weighted #Delta#phi Mixed Event V4 Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaPhiMixV4Pt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinPhi,phiArray,lptbins,fPtAssocArray);
	  fHistDeltaPhiMixV4Pt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaPhiMixV4Pt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#phi");
	  fHistDeltaPhiMixV4Pt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaPhiMixV4Pt[i][imult][isign][imc]);


	  //etaNear
	  snprintf(histname,BufferSize,"fHistDeltaEtaN_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Near-Side #Delta#eta Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaN[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaN[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaN[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaN[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaN[i][imult][isign][imc]);
	  
	  snprintf(histname,BufferSize,"fHistDeltaEtaNPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Near-Side P_{T} Weighted #Delta#eta Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNPt[i][imult][isign][imc]);
	  
	  snprintf(histname,BufferSize,"fHistDeltaEtaNMix_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Near-Side #Delta#eta Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNMix[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNMix[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNMix[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNMix[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNMix[i][imult][isign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaEtaNMixV2_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Near-Side #Delta#eta Mixed Event V2 Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNMixV2[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNMixV2[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNMixV2[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNMixV2[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNMixV2[i][imult][isign][imc]);

  snprintf(histname,BufferSize,"fHistDeltaEtaNMixV3_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Near-Side #Delta#eta Mixed Event V2 Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNMixV3[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNMixV3[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNMixV3[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNMixV3[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNMixV3[i][imult][isign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaEtaNMixV4_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Near-Side #Delta#eta Mixed Event V4 Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNMixV4[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNMixV4[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNMixV4[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNMixV4[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNMixV4[i][imult][isign][imc]);

	  
	  snprintf(histname,BufferSize,"fHistDeltaEtaNMixPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Near-Side P_{T} Weighted #Delta#eta Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNMixPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNMixPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNMixPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNMixPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNMixPt[i][imult][isign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaEtaNMixV2Pt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Near-Side P_{T} Weighted #Delta#eta Mixed Event V2 Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNMixV2Pt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNMixV2Pt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNMixV2Pt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNMixV2Pt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNMixV2Pt[i][imult][isign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaEtaNMixV3Pt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Near-Side P_{T} Weighted #Delta#eta Mixed Event V2 Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNMixV3Pt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNMixV3Pt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNMixV3Pt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNMixV3Pt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNMixV3Pt[i][imult][isign][imc]);

	   snprintf(histname,BufferSize,"fHistDeltaEtaNMixV4Pt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Near-Side P_{T} Weighted #Delta#eta V4 Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaNMixV4Pt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaNMixV4Pt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaNMixV4Pt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaNMixV4Pt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaNMixV4Pt[i][imult][isign][imc]);

	  //Away Eta
	  snprintf(histname,BufferSize,"fHistDeltaEtaA_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Away-Side #Delta#eta Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaA[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaA[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaA[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaA[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaA[i][imult][isign][imc]);
	  
	  snprintf(histname,BufferSize,"fHistDeltaEtaAPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Away-Side P_{T} Weighted #Delta#eta Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAPt[i][imult][isign][imc]);
	  
	  snprintf(histname,BufferSize,"fHistDeltaEtaAMix_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Away-Side #Delta#eta Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAMix[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAMix[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAMix[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAMix[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAMix[i][imult][isign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaEtaAMixV2_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Away-Side #Delta#eta Mixed Event V2 Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAMixV2[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAMixV2[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAMixV2[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAMixV2[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAMixV2[i][imult][isign][imc]);


	  snprintf(histname,BufferSize,"fHistDeltaEtaAMixV3_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Away-Side #Delta#eta Mixed Event V3 Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAMixV3[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAMixV3[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAMixV3[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAMixV3[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAMixV3[i][imult][isign][imc]);


	    snprintf(histname,BufferSize,"fHistDeltaEtaAMixV4_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Away-Side #Delta#eta Mixed Event V4 Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAMixV4[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAMixV4[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAMixV4[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAMixV4[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAMixV4[i][imult][isign][imc]);
	  
	  snprintf(histname,BufferSize,"fHistDeltaEtaAMixPt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Away-Side P_{T} Weighted #Delta#eta Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAMixPt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAMixPt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAMixPt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAMixPt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAMixPt[i][imult][isign][imc]);

	   snprintf(histname,BufferSize,"fHistDeltaEtaAMixV2Pt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Away-Side P_{T} Weighted V2 #Delta#eta Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAMixV2Pt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAMixV2Pt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAMixV2Pt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAMixV2Pt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAMixV2Pt[i][imult][isign][imc]);

	   snprintf(histname,BufferSize,"fHistDeltaEtaAMixV3Pt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Away-Side P_{T} Weighted V3 #Delta#eta Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAMixV3Pt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAMixV3Pt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAMixV3Pt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAMixV3Pt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAMixV3Pt[i][imult][isign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaEtaAMixV4Pt_P%d_C%d%s%s",i,imult,sign1[isign],cmc1[imc]);
	  snprintf(histtitle,BufferSize,"Away-Side P_{T} Weighted #Delta#eta V4 Mixed Event Distribution with %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],sign2[isign],cmc2[imc]);
	  fHistDeltaEtaAMixV4Pt[i][imult][isign][imc]=new TH2F(histname,histtitle,fnBinEta,etaArray,lptbins,fPtAssocArray);
	  fHistDeltaEtaAMixV4Pt[i][imult][isign][imc]->Sumw2();
	  fHistDeltaEtaAMixV4Pt[i][imult][isign][imc]->GetXaxis()->SetTitle("#Delta#eta");
	  fHistDeltaEtaAMixV4Pt[i][imult][isign][imc]->GetYaxis()->SetTitle("p_{T}");
	  fOutput->Add(fHistDeltaEtaAMixV4Pt[i][imult][isign][imc]);


      //====
	}//end isignloop
      snprintf(histname,BufferSize,"fHistDeltaPhiEta_P%d_C%d%s",i,imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"#Delta#phi-#Delta#eta %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEta[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEta[i][imult][imc]->Sumw2();
      fHistDeltaPhiEta[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEta[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEta[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEta[i][imult][imc]);
      
      
      snprintf(histname,BufferSize,"fHistDeltaPhiEtaMix_P%d_C%d%s",i,imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"#Delta#phi-#Delta#eta from Mixed Events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaMix[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEtaMix[i][imult][imc]->Sumw2();
      fHistDeltaPhiEtaMix[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEtaMix[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEtaMix[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEtaMix[i][imult][imc]);

      snprintf(histname,BufferSize,"fHistDeltaPhiEtaMixV2_P%d_C%d%s",i,imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"#Delta#phi-#Delta#eta V2 from Mixed Events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaMixV2[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEtaMixV2[i][imult][imc]->Sumw2();
      fHistDeltaPhiEtaMixV2[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEtaMixV2[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEtaMixV2[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEtaMixV2[i][imult][imc]);

 snprintf(histname,BufferSize,"fHistDeltaPhiEtaMixV3_P%d_C%d%s",i,imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"#Delta#phi-#Delta#eta V3 from Mixed Events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaMixV3[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEtaMixV3[i][imult][imc]->Sumw2();
      fHistDeltaPhiEtaMixV3[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEtaMixV3[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEtaMixV3[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEtaMixV3[i][imult][imc]);

	snprintf(histname,BufferSize,"fHistDeltaPhiEtaMixV4_P%d_C%d%s",i,imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"#Delta#phi-#Delta#eta V4 from Mixed Events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaMixV4[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEtaMixV4[i][imult][imc]->Sumw2();
      fHistDeltaPhiEtaMixV4[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEtaMixV4[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEtaMixV4[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEtaMixV4[i][imult][imc]);
      
      snprintf(histname,BufferSize,"fHistPhiEtaTrigPt_P%d_C%d%s",i,imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"P_{T}-Weighted #phi-#eta distribution in triggered events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistPhiEtaTrigPt[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistPhiEtaTrigPt[i][imult][imc]->Sumw2();
      fHistPhiEtaTrigPt[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistPhiEtaTrigPt[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistPhiEtaTrigPt[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistPhiEtaTrigPt[i][imult][imc]);
    
      snprintf(histname,BufferSize,"fHistDeltaPhiEtaPt_P%d_C%d%s",i,imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"P_{T}-Weighted #Delta#phi-#Delta#eta %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaPt[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEtaPt[i][imult][imc]->Sumw2();
      fHistDeltaPhiEtaPt[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEtaPt[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEtaPt[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEtaPt[i][imult][imc]);
      
      snprintf(histname,BufferSize,"fHistDeltaPhiEtaMixPt_P%d_C%d%s",i,imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"P_{T}-Weighted #Delta#phi-#Delta#eta from Mixed Events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaMixPt[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEtaMixPt[i][imult][imc]->Sumw2();
      fHistDeltaPhiEtaMixPt[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEtaMixPt[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEtaMixPt[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEtaMixPt[i][imult][imc]);

      snprintf(histname,BufferSize,"fHistDeltaPhiEtaMixV2Pt_P%d_C%d%s",i,imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"P_{T}-Weighted #Delta#phi-#Delta#eta V2 from Mixed Events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaMixV2Pt[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEtaMixV2Pt[i][imult][imc]->Sumw2();
      fHistDeltaPhiEtaMixV2Pt[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEtaMixV2Pt[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEtaMixV2Pt[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEtaMixV2Pt[i][imult][imc]);

  snprintf(histname,BufferSize,"fHistDeltaPhiEtaMixV3Pt_P%d_C%d%s",i,imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"P_{T}-Weighted #Delta#phi-#Delta#eta V3 from Mixed Events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaMixV3Pt[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEtaMixV3Pt[i][imult][imc]->Sumw2();
      fHistDeltaPhiEtaMixV3Pt[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEtaMixV3Pt[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEtaMixV3Pt[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEtaMixV3Pt[i][imult][imc]);

      snprintf(histname,BufferSize,"fHistDeltaPhiEtaMixV4Pt_P%d_C%d%s",i,imult,cmc1[imc]);
      snprintf(histtitle,BufferSize,"P_{T}-Weighted #Delta#phi-#Delta#eta V4 from Mixed Events %3.1f<p_{T}^{Trig}<%3.1f %dMult%d%s",fPtTrigArray[i],fPtTrigArray[i+1],fCentArrayMin[imult],fCentArrayMax[imult],cmc2[imc]);
      fHistDeltaPhiEtaMixV4Pt[i][imult][imc]=new TH3F(histname,histtitle,fnBinPhiEtaPhi,phiEtaArrayPhi,fnBinPhiEtaEta,phiEtaArrayEta,lptbins,fPtAssocArray);
      fHistDeltaPhiEtaMixV4Pt[i][imult][imc]->Sumw2();
      fHistDeltaPhiEtaMixV4Pt[i][imult][imc]->GetXaxis()->SetTitle("#phi");
      fHistDeltaPhiEtaMixV4Pt[i][imult][imc]->GetYaxis()->SetTitle("#eta");
      fHistDeltaPhiEtaMixV4Pt[i][imult][imc]->GetZaxis()->SetTitle("p_{T}");
      fOutput->Add(fHistDeltaPhiEtaMixV4Pt[i][imult][imc]);

      //Three-Particle Histograms
      for(int ipt=0;ipt<fNAPt3Bins;ipt++){
	for(int iSign=0;iSign<4;iSign++){
	  snprintf(histname,BufferSize,"fHistDeltaPhiPhi_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Raw #Delta#phi-#Delta#phi %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaPhiPhi[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinPhi3,fdPhiMin,fdPhiMax,fnBinPhi3,fdPhiMin,fdPhiMax);
	  fHistDeltaPhiPhi[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaPhiPhi[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#phi_{1}");
	  fHistDeltaPhiPhi[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#phi_{2}");
	  fOutput->Add(fHistDeltaPhiPhi[i][ipt][imult][iSign][imc]);
	  
	  snprintf(histname,BufferSize,"fHistDeltaPhiPhiMix_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Mixed #Delta#phi-#Delta#phi %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaPhiPhiMix[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinPhi3,fdPhiMin,fdPhiMax,fnBinPhi3,fdPhiMin,fdPhiMax);
	  fHistDeltaPhiPhiMix[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaPhiPhiMix[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#phi_{1}");
	  fHistDeltaPhiPhiMix[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#phi_{2}");
	  fOutput->Add(fHistDeltaPhiPhiMix[i][ipt][imult][iSign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaPhiPhiMixV2_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Mixed #Delta#phi-#Delta#phi V2 %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaPhiPhiMixV2[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinPhi3,fdPhiMin,fdPhiMax,fnBinPhi3,fdPhiMin,fdPhiMax);
	  fHistDeltaPhiPhiMixV2[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaPhiPhiMixV2[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#phi_{1}");
	  fHistDeltaPhiPhiMixV2[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#phi_{2}");
	  fOutput->Add(fHistDeltaPhiPhiMixV2[i][ipt][imult][iSign][imc]);

  snprintf(histname,BufferSize,"fHistDeltaPhiPhiMixV3_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Mixed #Delta#phi-#Delta#phi V3 %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaPhiPhiMixV3[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinPhi3,fdPhiMin,fdPhiMax,fnBinPhi3,fdPhiMin,fdPhiMax);
	  fHistDeltaPhiPhiMixV3[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaPhiPhiMixV3[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#phi_{1}");
	  fHistDeltaPhiPhiMixV3[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#phi_{2}");
	  fOutput->Add(fHistDeltaPhiPhiMixV3[i][ipt][imult][iSign][imc]);

	    snprintf(histname,BufferSize,"fHistDeltaPhiPhiMixV4_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Mixed #Delta#phi-#Delta#phi V4 %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaPhiPhiMixV4[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinPhi3,fdPhiMin,fdPhiMax,fnBinPhi3,fdPhiMin,fdPhiMax);
	  fHistDeltaPhiPhiMixV4[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaPhiPhiMixV4[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#phi_{1}");
	  fHistDeltaPhiPhiMixV4[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#phi_{2}");
	  fOutput->Add(fHistDeltaPhiPhiMixV4[i][ipt][imult][iSign][imc]);

	    snprintf(histname,BufferSize,"fHistDeltaPhiPhiMixV2V2V4_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Mixed #Delta#phi-#Delta#phi %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaPhiPhiMixV2V2V4[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinPhi3,fdPhiMin,fdPhiMax,fnBinPhi3,fdPhiMin,fdPhiMax);
	  fHistDeltaPhiPhiMixV2V2V4[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaPhiPhiMixV2V2V4[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#phi_{1}");
	  fHistDeltaPhiPhiMixV2V2V4[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#phi_{2}");
	  fOutput->Add(fHistDeltaPhiPhiMixV2V2V4[i][ipt][imult][iSign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaPhiPhiSS_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Soft-Soft #Delta#phi-#Delta#phi %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaPhiPhiSS[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinPhi3,fdPhiMin,fdPhiMax,fnBinPhi3,fdPhiMin,fdPhiMax);
	  fHistDeltaPhiPhiSS[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaPhiPhiSS[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#phi_{1}");
	  fHistDeltaPhiPhiSS[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#phi_{2}");
	  fOutput->Add(fHistDeltaPhiPhiSS[i][ipt][imult][iSign][imc]);

	   snprintf(histname,BufferSize,"fHistDeltaPhiPhiSSV2_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Soft-Soft V2  #Delta#phi-#Delta#phi %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaPhiPhiSSV2[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinPhi3,fdPhiMin,fdPhiMax,fnBinPhi3,fdPhiMin,fdPhiMax);
	  fHistDeltaPhiPhiSSV2[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaPhiPhiSSV2[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#phi_{1}");
	  fHistDeltaPhiPhiSSV2[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#phi_{2}");
	  fOutput->Add(fHistDeltaPhiPhiSSV2[i][ipt][imult][iSign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaPhiPhiSSV3_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Soft-Soft V3 #Delta#phi-#Delta#phi %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaPhiPhiSSV3[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinPhi3,fdPhiMin,fdPhiMax,fnBinPhi3,fdPhiMin,fdPhiMax);
	  fHistDeltaPhiPhiSSV3[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaPhiPhiSSV3[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#phi_{1}");
	  fHistDeltaPhiPhiSSV3[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#phi_{2}");
	  fOutput->Add(fHistDeltaPhiPhiSSV3[i][ipt][imult][iSign][imc]);

	   snprintf(histname,BufferSize,"fHistDeltaPhiPhiSSV4_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Soft-Soft V4 #Delta#phi-#Delta#phi %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaPhiPhiSSV4[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinPhi3,fdPhiMin,fdPhiMax,fnBinPhi3,fdPhiMin,fdPhiMax);
	  fHistDeltaPhiPhiSSV4[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaPhiPhiSSV4[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#phi_{1}");
	  fHistDeltaPhiPhiSSV4[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#phi_{2}");
	  fOutput->Add(fHistDeltaPhiPhiSSV4[i][ipt][imult][iSign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaEtaEta_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Raw #Delta#eta-#Delta#eta %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaEtaEta[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinEta3,-etaEdge,etaEdge,fnBinEta3,-etaEdge,etaEdge);
	  fHistDeltaEtaEta[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaEtaEta[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#eta_{1}");
	  fHistDeltaEtaEta[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#eta_{2}");
	  fOutput->Add(fHistDeltaEtaEta[i][ipt][imult][iSign][imc]);

snprintf(histname,BufferSize,"fHistDeltaEtaEtaMix_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Mixed #Delta#eta-#Delta#eta %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinEta3,-etaEdge,etaEdge,fnBinEta3,-etaEdge,etaEdge);
	  fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#eta_{1}");
	  fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#eta_{2}");
	  fOutput->Add(fHistDeltaEtaEtaMix[i][ipt][imult][iSign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaEtaEtaMixV2_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Mixed #Delta#eta-#Delta#eta V2 %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaEtaEtaMixV2[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinEta3,-etaEdge,etaEdge,fnBinEta3,-etaEdge,etaEdge);
	  fHistDeltaEtaEtaMixV2[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaEtaEtaMixV2[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#eta_{1}");
	  fHistDeltaEtaEtaMixV2[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#eta_{2}");
	  fOutput->Add(fHistDeltaEtaEtaMixV2[i][ipt][imult][iSign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaEtaEtaMixV3_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Mixed #Delta#eta-#Delta#eta V3 %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaEtaEtaMixV3[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinEta3,-etaEdge,etaEdge,fnBinEta3,-etaEdge,etaEdge);
	  fHistDeltaEtaEtaMixV3[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaEtaEtaMixV3[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#eta_{1}");
	  fHistDeltaEtaEtaMixV3[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#eta_{2}");
	  fOutput->Add(fHistDeltaEtaEtaMixV3[i][ipt][imult][iSign][imc]);


	  snprintf(histname,BufferSize,"fHistDeltaEtaEtaMixV4_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Mixed #Delta#eta-#Delta#eta V4 %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaEtaEtaMixV4[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinEta3,-etaEdge,etaEdge,fnBinEta3,-etaEdge,etaEdge);
	  fHistDeltaEtaEtaMixV4[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaEtaEtaMixV4[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#eta_{1}");
	  fHistDeltaEtaEtaMixV4[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#eta_{2}");
	  fOutput->Add(fHistDeltaEtaEtaMixV4[i][ipt][imult][iSign][imc]);

	  snprintf(histname,BufferSize,"fHistDeltaEtaEtaMixV2V2V4_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Mixed #Delta#eta-#Delta#eta V2V2V4 %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
	  fHistDeltaEtaEtaMixV2V2V4[i][ipt][imult][iSign][imc]=new TH2F(histname,histtitle,fnBinEta3,-etaEdge,etaEdge,fnBinEta3,-etaEdge,etaEdge);
	  fHistDeltaEtaEtaMixV2V2V4[i][ipt][imult][iSign][imc]->Sumw2();
	  fHistDeltaEtaEtaMixV2V2V4[i][ipt][imult][iSign][imc]->GetXaxis()->SetTitle("#Delta#eta_{1}");
	  fHistDeltaEtaEtaMixV2V2V4[i][ipt][imult][iSign][imc]->GetYaxis()->SetTitle("#Delta#eta_{2}");
	  fOutput->Add(fHistDeltaEtaEtaMixV2V2V4[i][ipt][imult][iSign][imc]);

snprintf(histname,BufferSize,"fHistDeltaEtaEtaSS_P%dp%d_C%d%s%s",i,ipt,imult,cmc1[imc],sign31[iSign]);
	  snprintf(histtitle,BufferSize,"Soft-Soft #Delta#eta-#Delta#eta %3.1f<p_{T}^{Trig}<%3.1f %3.2f<p_{T}^{Assoc}<%3.2f %dMult%d%s%s",fPtTrigArray[i],fPtTrigArray[i+1],fPtAssoc3Array1[ipt],fPtAssoc3Array2[ipt],fCentArrayMin[imult],fCentArrayMax[imult],sign32[iSign],cmc2[imc]);
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
/////////////////////////////////
Int_t AliAnalysisTaskDiHadron::CheckVertex(const AliESDEvent *rESD){
  //checks whether the vertex passes cuts
  Int_t rGood=-1;
  Float_t vtx[3];
  vtx[0]=rESD->GetPrimaryVertex()->GetX();
  vtx[1]=rESD->GetPrimaryVertex()->GetY();
  vtx[2]=rESD->GetPrimaryVertex()->GetZ();
  if((vtx[0]*vtx[0]+vtx[1]*vtx[1])<9&&fabs(vtx[2])<fVzCut) rGood=0; //vertex out of beam pipe
  //if(fabs(vtx[2])<fVzCut)rGood=0;//Vertex Z cut
  if(fDEBUG)Printf("vtxZ %f",vtx[2]);
  for(int i=0;i<fNVertexBins;i++){
    if(vtx[2]>fVertexArray[i]&&vtx[2]<=fVertexArray[i+1]&&rGood==0)rGood=i;
  }
  return rGood;
}
///////////////////////////
Int_t AliAnalysisTaskDiHadron::CheckVertexAOD(const AliAODEvent *rAOD){
  //checks whether the vertex passes cuts
  Int_t rGood=-1;
  Float_t vtx[3];
  vtx[0]=rAOD->GetPrimaryVertex()->GetX();
  vtx[1]=rAOD->GetPrimaryVertex()->GetY();
  vtx[2]=rAOD->GetPrimaryVertex()->GetZ();
  if((vtx[0]*vtx[0]+vtx[1]*vtx[1])<9&&fabs(vtx[2])<fVzCut) rGood=0; //vertex out of beam pipe
  //if(fabs(vtx[2])<fVzCut)rGood=0;//Vertex Z cut
  if(fDEBUG)Printf("vtxZ %f",vtx[2]);
  for(int i=0;i<fNVertexBins;i++){
    if(vtx[2]>fVertexArray[i]&&vtx[2]<=fVertexArray[i+1]&&rGood==0)rGood=i;
  }
  return rGood;
}
///////////////////////////////
Int_t AliAnalysisTaskDiHadron::CheckTrigger(const AliESDEvent *rESD){
  //checks whether the trigger passes cuts
  if(fDEBUG)Printf("Checking Trigger");
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
///////////////////////////////////////////
Int_t AliAnalysisTaskDiHadron::CheckTriggerAOD(const AliAODEvent *rAOD){
  //checks whether the trigger passes cuts
  if(fDEBUG)Printf("Checking Trigger");
  Int_t rGood=0;
  TString trigID=rAOD->GetFiredTriggerClasses();
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
/////////////////////////////////////////////

Int_t AliAnalysisTaskDiHadron::TrackCuts(const AliESDEvent *rESD, Float_t *rPt, Float_t *rEta, Float_t *rPhi, Short_t *rCharge, Float_t *rEff, Float_t *rV2, Float_t *rV3, Float_t *rV4, Int_t **rPtAssoc3, Int_t *rNPtAssoc3, Int_t *rGoodTracks){
  if(fDEBUG) Printf("Selecting Tracks");
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
    if(fDEBUG) Printf("Pt%2.2f Eta%2.2f Phi%2.2f ", sPt,sEta,sPhi);
    if(sPhi<fdPhiMin)sPhi+=2*fPi;
    if(sPhi>fdPhiMax)sPhi-=2*fPi;
    if(sPt<fMinPtAssoc||sPt>fMaxPtAssoc)continue;//set Pt range
    if(fabs(sEta)>fEtaCut)continue;//set Eta Range
    if(!sCharge)continue;
    sNcls=eSDtrack->GetTPCNcls();
    //if(fDEBUG)Printf("NCLS%d",sNcls);
    if(sNcls<fMinClustersTPC)continue;
    sNclsF=eSDtrack->GetTPCnclsS();
    if((1-1.0*sNclsF/sNcls)<fMinClusterRatio)continue;//Clusters fit/ total
    sChi=(eSDtrack->GetTPCchi2())/sNcls;
    if(sChi>fMaxTPCchi2)continue;
    sITScls=eSDtrack->GetNcls(0);
    if(sITScls<fMinClustersITS)continue;
    eSDtrack->GetImpactParameters(sb,sbCov);
    if(fDEBUG)Printf("dca %2.2f %2.2f",sb[0],sb[1]);
    if(!fDCA2D&&(sb[0]*sb[0]+sb[1]*sb[1])>(fMaxDCA*fMaxDCA))continue;//DCA cut
    if(fDCA2D==1&&(sb[0]*sb[0]/fMaxDCAXY/fMaxDCAXY+sb[1]*sb[1]/fMaxDCAZ/fMaxDCAZ)>1)continue;
    if(fDCA2D==2&&(0.35+0.42*std::pow(double(sPt),-0.9))<(sb[0]*sb[0]))continue;
    if(eSDtrack->GetKinkIndex(0)>0)continue;//removes kinked tracks
    if(!eSDtrack->GetStatus()&AliESDtrack::kTPCrefit&&fTPCRefit)continue;//refit in TPC
    if((fITSRefit==1||(fITSRefit==2&&sPt>5))&&!eSDtrack->GetStatus()&AliESDtrack::kITSrefit)continue;//refit of its tracks either for none,all, or >5 GeV/c
    if(fDEBUG)Printf("SPD %d %d ", eSDtrack->HasPointOnITSLayer(0), eSDtrack->HasPointOnITSLayer(1));
    if(fSPDCut&&!eSDtrack->HasPointOnITSLayer(0)&&!eSDtrack->HasPointOnITSLayer(1))continue;
    if(fDEBUG)Printf("Pass \n");
    rPt[rGoodTracks[0]]=sPt;
    rEta[rGoodTracks[0]]=sEta;
    rPhi[rGoodTracks[0]]=sPhi;
    rCharge[rGoodTracks[0]]=sCharge;
    if(fEfficiencyCorr){
    if(sPt<fEffFitPt)rEff[rGoodTracks[0]]=1./fFitLow->Eval(sPt);
    else rEff[rGoodTracks[0]]=1./fFitHigh->Eval(sPt);
    }
    else rEff[rGoodTracks[0]]=1;
    if(rEff[rGoodTracks[0]]!=rEff[rGoodTracks[0]]||rEff[rGoodTracks[0]]>1E8||rEff[rGoodTracks[0]]<-1E8){
      Printf("Efficiency Error %f %f",rEff[rGoodTracks[0]],rPt[rGoodTracks[0]]);
      continue;
    }
    if(sPt>leadPt)lead=rGoodTracks[0];
    if(sPt<fV2FitPt)rV2[rGoodTracks[0]]=fFitLowV2->Eval(sPt);
    else rV2[rGoodTracks[0]]=fFitHighV2->Eval(sPt);
    if(sPt<fV3FitPt)rV3[rGoodTracks[0]]=fFitLowV3->Eval(sPt);
    else rV3[rGoodTracks[0]]=fFitHighV3->Eval(sPt);
    if(sPt<fV4FitPt)rV4[rGoodTracks[0]]=fFitLowV4->Eval(sPt);
    else rV4[rGoodTracks[0]]=fFitHighV4->Eval(sPt);
    if(rV2[rGoodTracks[0]]!=rV2[rGoodTracks[0]]||rV2[rGoodTracks[0]]>1E8||rV2[rGoodTracks[0]]<-1E8){
      Printf("V2 Error %f %f",rV2[rGoodTracks[0]],rPt[rGoodTracks[0]]);
      continue;
    }
      if(rV4[rGoodTracks[0]]!=rV4[rGoodTracks[0]]||rV4[rGoodTracks[0]]>1E8||rV4[rGoodTracks[0]]<-1E8){
	Printf("V4 Error %f %f",rV4[rGoodTracks[0]],rPt[rGoodTracks[0]]);
      continue;
    }
      
      //Printf("V2 %2.2f V4 %2.4f 1.15V2V2 %2.4f",rV2[rGoodTracks[0]],rV4[rGoodTracks[0]],1.15*pow(rV2[rGoodTracks[0]],2));
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
/////////////////////////////////////////////

Int_t AliAnalysisTaskDiHadron::TrackCutsAOD(const AliAODEvent *rAOD, Float_t *rPt, Float_t *rEta, Float_t *rPhi, Short_t *rCharge, Float_t *rEff, Float_t *rV2, Float_t *rV3, Float_t *rV4, Int_t **rPtAssoc3, Int_t *rNPtAssoc3, Int_t *rGoodTracks){
  if(fDEBUG) Printf("Selecting Tracks");
    //fills arrays with all of the tracks passing cuts
  rGoodTracks[0]=0;
  Int_t lead=0;
  Float_t leadPt=0;
  Int_t rTrack=fAOD->GetNTracks();
  Float_t sPt, sEta, sPhi, sChi, sb[2];
  Int_t sNcls, sNclsF, sITScls;
  Short_t sCharge;
  for(int iTrack=0;iTrack<rTrack;iTrack++){
    AliAODTrack *aodTrack=rAOD->GetTrack(iTrack);
    sPt=aodTrack->Pt();
    sEta=aodTrack->Eta();
    sPhi=aodTrack->Phi();
    sCharge=aodTrack->Charge();
    if(fDEBUG) Printf("Pt%2.2f Eta%2.2f Phi%2.2f ", sPt,sEta,sPhi);
    if(sPhi<fdPhiMin)sPhi+=2*fPi;
    if(sPhi>fdPhiMax)sPhi-=2*fPi;
    if(sPt<fMinPtAssoc||sPt>fMaxPtAssoc)continue;//set Pt range
    if(fabs(sEta)>fEtaCut)continue;//set Eta Range
    if(!sCharge)continue;
    sNcls=aodTrack->GetTPCNcls();
    if(sNcls<fMinClustersTPC)continue;
    sNclsF=aodTrack->GetTPCSharedMap().CountBits();
    if((1-1.0*sNclsF/sNcls)<fMinClusterRatio)continue;//Clusters shared/ total;
    sChi=aodTrack->Chi2perNDF();
    if(sChi>fMaxTPCchi2)continue;
    sITScls=aodTrack->GetNcls(0);
    if(sITScls<fMinClustersITS)continue;
    sb[0]=aodTrack->DCA();
    sb[1]=aodTrack->ZAtDCA();
    if(fDEBUG)Printf("dca %2.2f %2.2f",sb[0],sb[1]);
    if(!fDCA2D&&(sb[0]*sb[0]+sb[1]*sb[1])>(fMaxDCA*fMaxDCA))continue;//DCA cut
    if(fDCA2D==1&&(sb[0]*sb[0]/fMaxDCAXY/fMaxDCAXY+sb[1]*sb[1]/fMaxDCAZ/fMaxDCAZ)>1)continue;
    if(fDCA2D==2&&(0.35+0.42*std::pow(double(sPt),-0.9))<(sb[0])*sb[0])continue;
    //if(eSDtrack->GetKinkIndex(0)>0)continue;//removes kinked tracks
    if(!aodTrack->IsPrimaryCandidate())continue;//I assume this removes kinks
    //if(!aodTrack->GetStatus()&AliAODTrack::kTPCrefit&&fTPCRefit)continue;//refit in TPC
    //if((fITSRefit==1||(fITSRefit==2&&sPt>5))&&!aodTrack->GetStatus()&AliAODTrack::kITSrefit)continue;//refit of its tracks either for none,all, or >5 GeV/c
    if(fDEBUG)Printf("SPD %d %d ", aodTrack->HasPointOnITSLayer(0), aodTrack->HasPointOnITSLayer(1));
    if(fSPDCut&&!aodTrack->HasPointOnITSLayer(0)&&!aodTrack->HasPointOnITSLayer(1))continue;
    if(fDEBUG)Printf("Pass \n");
    rPt[rGoodTracks[0]]=sPt;
    rEta[rGoodTracks[0]]=sEta;
    rPhi[rGoodTracks[0]]=sPhi;
    rCharge[rGoodTracks[0]]=sCharge;
    if(fEfficiencyCorr){
    if(sPt<fEffFitPt)rEff[rGoodTracks[0]]=1./fFitLow->Eval(sPt);
    else rEff[rGoodTracks[0]]=1./fFitHigh->Eval(sPt);
    }
    else rEff[rGoodTracks[0]]=1;
    if(rEff[rGoodTracks[0]]!=rEff[rGoodTracks[0]]||rEff[rGoodTracks[0]]>1E8||rEff[rGoodTracks[0]]<-1E8){
      Printf("Efficiency Error %f %f",rEff[rGoodTracks[0]],rPt[rGoodTracks[0]]);
      continue;
    }
    if(sPt>leadPt)lead=rGoodTracks[0];
    if(sPt<fV2FitPt)rV2[rGoodTracks[0]]=fFitLowV2->Eval(sPt);
    else rV2[rGoodTracks[0]]=fFitHighV2->Eval(sPt);
    if(sPt<fV3FitPt)rV3[rGoodTracks[0]]=fFitLowV3->Eval(sPt);
    else rV3[rGoodTracks[0]]=fFitHighV3->Eval(sPt);
    if(sPt<fV4FitPt)rV4[rGoodTracks[0]]=fFitLowV4->Eval(sPt);
    else rV4[rGoodTracks[0]]=fFitHighV4->Eval(sPt);
    if(rV2[rGoodTracks[0]]!=rV2[rGoodTracks[0]]||rV2[rGoodTracks[0]]>1E8||rV2[rGoodTracks[0]]<-1E8){
      Printf("V2 Error %f %f",rV2[rGoodTracks[0]],rPt[rGoodTracks[0]]);
      continue;
    }
      if(rV4[rGoodTracks[0]]!=rV4[rGoodTracks[0]]||rV4[rGoodTracks[0]]>1E8||rV4[rGoodTracks[0]]<-1E8){
	Printf("V4 Error %f %f",rV4[rGoodTracks[0]],rPt[rGoodTracks[0]]);
      continue;
    }
      
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
///////////////////////////////////////////////////////

Int_t AliAnalysisTaskDiHadron::TrackCutsMC(AliMCEvent *rMC, Float_t *rPt, Float_t *rEta, Float_t *rPhi, Short_t *rCharge, Float_t *rEff, Float_t *rV2, Float_t *rV3, Float_t *rV4, Int_t **rPtAssoc3, Int_t *rNPtAssoc3, Int_t *rGoodTracks){
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
    rV2[rGoodTracks[1]]=0;
    rV3[rGoodTracks[1]]=0;
    rV4[rGoodTracks[1]]=0;
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
//---------------------------------------------------------
Int_t AliAnalysisTaskDiHadron::TrackCutsSim(Float_t *rPt, Float_t *rEta, Float_t *rPhi, Short_t *rCharge, Float_t *rEff, Float_t *rV2, Float_t *rV3, Float_t *rV4, Int_t **rPtAssoc3, Int_t *rNPtAssoc3, Int_t *rGoodTracks){

  float v2=fFitHighV2->Eval(3.);
  float v3=fFitHighV3->Eval(3.);
  float v4=fFitHighV4->Eval(3.);
   Int_t lead=0;
  rGoodTracks[0]=0;

  TF1 *SimFlow=new TF1("SimFlow","1+2*[0]*cos(2*x)+2*[1]*cos(3*(x-[2]))+2*[3]*cos(4*x)",-TMath::Pi(),TMath::Pi());
  Float_t V3PlaneAngle=2*TMath::Pi()*gRandom->Rndm();
  SimFlow->SetParameters(v2,v3,V3PlaneAngle,v4);
 // SimFlow->SetParameters(0,0,0,0);
  TF1 *SimNear=new TF1("SimNear","exp(-0.5*x*x/[0]/[0])",-TMath::Pi(),TMath::Pi());
  SimNear->SetParameter(0,0.3);
  TF1 *SimAway=new TF1("SimAway","exp(-0.5*x*x/[0]/[0])",-TMath::Pi(),TMath::Pi());
 SimAway->SetParameter(0,0.5);//0.5 for DiJet 0.3 for cone
  
  TF1 *SimAway2=new TF1("SimAway2","exp(-0.5*(x-[1])*(x-[1])/[0]/[0])+exp(-0.5*(x+[1])*(x+[1])/[0]/[0])",-TMath::Pi(),TMath::Pi());
  SimAway2->SetParameter(0,0.3);//0.5 for DiJet 0.3 for cone
  SimAway2->SetParameter(1,1.4);//Cone Angle

  // TF1 *AwayProb=new TF1("AwayProb","[0]");
   TF1 *AwayProb=new TF1("AwayProb","[0]*cos(2*x)*cos(2*x)");
   //AwayProb->SetParameter(0,0.5);
  AwayProb->SetParameter(0,1);
  Int_t AwayDeflected=0;
  

  // TF1 *SimAway=new TF1("SimAway","0.035+2*[0]*[0]*cos(2*x)+2*[1]*[1]*cos(3*x)+2*[2]*[2]*cos(4*x)",1.0,TMath::Pi()*2-1.0);
  // SimAway->SetParameters(v2,v3,v4);
  // TF1 *SimAway2=new TF1("SimAway2","0.035+2*[0]*[0]*cos(2*x)+2*[1]*[1]*cos(3*x)+2*[2]*[2]*cos(4*x)",1.05,TMath::Pi()*2-1.05);
  //SimAway->SetParameters(v2,v3,v4);
 
  Float_t RPAngle=2*TMath::Pi()*gRandom->Rndm();
  Float_t TrigAngle;
  Float_t sPt,sPhi;
  Int_t InAccpt;
  Int_t AccptPercent=4;//1 over this is % in aceptance on away-side 
  Int_t AwaySidePM=0;

  Int_t AwaySide1=1;
  //Use SimAway1 or 2
  // if(gRandom->Rndm()<AwayProb->Eval(RPAngle))AwaySide1=1;
  //else AwaySide1=2;

  for(int i=0;i<gRandom->Poisson(fSimNBgPart);i++){
    sPt=1.5;
    rPt[rGoodTracks[0]]=sPt;
    rEta[rGoodTracks[0]]=0;
    sPhi=SimFlow->GetRandom()+RPAngle;
    if(sPhi<fdPhiMin)sPhi+=2*fPi;
    if(sPhi>fdPhiMax)sPhi-=2*fPi;
    rPhi[rGoodTracks[0]]=sPhi;
    rCharge[rGoodTracks[0]]=1;
    rEff[rGoodTracks[0]]=1;
    rV2[rGoodTracks[0]]=v2;
    rV3[rGoodTracks[0]]=v3;
    rV4[rGoodTracks[0]]=v4;
    rNPtAssoc3[rGoodTracks[0]]=0;
    for(int apt3=0;apt3<fNAPt3Bins;apt3++){
      if(sPt<fPtAssoc3Array2[apt3]&&sPt>=fPtAssoc3Array1[apt3]){
	rPtAssoc3[rGoodTracks[0]][rNPtAssoc3[rGoodTracks[0]]]=apt3;
	rNPtAssoc3[rGoodTracks[0]]++;
      }
    } 
    rGoodTracks[0]++;
  }
  for(int i=0;i<gRandom->Poisson(fSimNJet);i++){
    TrigAngle=SimFlow->GetRandom()+RPAngle;
     if(gRandom->Rndm()<AwayProb->Eval(TrigAngle-RPAngle))AwaySide1=1;
     else AwaySide1=2;
    sPhi=TrigAngle;
    if(sPhi<fdPhiMin)sPhi+=2*fPi;
    if(sPhi>fdPhiMax)sPhi-=2*fPi;
    sPt=3.1;
    rPt[rGoodTracks[0]]=sPt;
    rEta[rGoodTracks[0]]=0;
    rPhi[rGoodTracks[0]]=sPhi;
    rCharge[rGoodTracks[0]]=1;
    rEff[rGoodTracks[0]]=1;
    rV2[rGoodTracks[0]]=v2;
    rV3[rGoodTracks[0]]=v3;
    rV4[rGoodTracks[0]]=v4;
    rNPtAssoc3[rGoodTracks[0]]=0;
    lead=rGoodTracks[0];
    rGoodTracks[0]++;
    
    for(int k=0;k<gRandom->Poisson(fSimNJetPart);k++){
      sPhi=SimNear->GetRandom()+TrigAngle;
      if(sPhi<fdPhiMin)sPhi+=2*fPi;
      if(sPhi>fdPhiMax)sPhi-=2*fPi;
      sPt=1.5;
      rPt[rGoodTracks[0]]=sPt;
      rEta[rGoodTracks[0]]=0;
      rPhi[rGoodTracks[0]]=sPhi;
      rCharge[rGoodTracks[0]]=1;
      rEff[rGoodTracks[0]]=1;
      rV2[rGoodTracks[0]]=v2;
      rV3[rGoodTracks[0]]=v3;
      rV4[rGoodTracks[0]]=v4;
      rNPtAssoc3[rGoodTracks[0]]=0;
      for(int apt3=0;apt3<fNAPt3Bins;apt3++){
	if(sPt<fPtAssoc3Array2[apt3]&&sPt>=fPtAssoc3Array1[apt3]){
	  rPtAssoc3[rGoodTracks[0]][rNPtAssoc3[rGoodTracks[0]]]=apt3;
	  rNPtAssoc3[rGoodTracks[0]]++;
	}
      } 
      rGoodTracks[0]++;
    }

    if(gRandom->Rndm()<1./AccptPercent)InAccpt=1;
    else InAccpt=0;
    if(gRandom->Rndm()<0.5)AwaySidePM=0;
    else AwaySidePM=1;
    for(int k=0;k<gRandom->Poisson(InAccpt*AccptPercent*fSimNJetPart);k++){
      //sPhi=SimAway->GetRandom()+TrigAngle+TMath::Pi();
      if(AwaySide1==1)sPhi=SimAway->GetRandom();
      else(sPhi=SimAway2->GetRandom());
      if(AwayDeflected){
	if(sPhi>0&&AwaySidePM)sPhi=-sPhi;
	else if(sPhi<0&&!AwaySidePM)sPhi=-sPhi;
      }
	sPhi+=TrigAngle+fPi;
	if(sPhi<fdPhiMin)sPhi+=2*fPi;
	if(sPhi>fdPhiMax)sPhi-=2*fPi;
	sPt=1.5;
	rPt[rGoodTracks[0]]=sPt;
	rEta[rGoodTracks[0]]=0;
	rPhi[rGoodTracks[0]]=sPhi;
	rCharge[rGoodTracks[0]]=1;
	rEff[rGoodTracks[0]]=1;
	rV2[rGoodTracks[0]]=v2;
	rV3[rGoodTracks[0]]=v3;
	rV4[rGoodTracks[0]]=v4;
	rNPtAssoc3[rGoodTracks[0]]=0;
	for(int apt3=0;apt3<fNAPt3Bins;apt3++){
	  if(sPt<fPtAssoc3Array2[apt3]&&sPt>=fPtAssoc3Array1[apt3]){
	    rPtAssoc3[rGoodTracks[0]][rNPtAssoc3[rGoodTracks[0]]]=apt3;
	    rNPtAssoc3[rGoodTracks[0]]++;
	  }
	}
	rGoodTracks[0]++;
    }


  }//njet

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

    if((!fESD&&!fAOD)&&ievent==0){
      if(fDEBUG)Printf("Error: fESD not found");
      break;
    }
    if(!fMC&&ievent==1){
      break;
    }
    // Printf("fSimulate %d",fSimulate);
    if(fSimulate==1&&ievent==1) break;
    if(ievent==1&&!fMCHistos)break;//break out if MC event and we don't have fill of those set
    //Secondary check
    if(ievent==0){
      if(!fAODData){
	if(fESD->GetNumberOfTracks()<=0){
	  if(fDEBUG)Printf("Error: no tracks");
	  break;
	}
      }
      else{
	if(fAOD->GetNTracks()<=0){
	  if(fDEBUG)Printf("Error: no tracks");
	  break;
	}
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
      if(!fAODData){
	if(!CheckTrigger(fESD)) break;
      }
      else{
	if(!CheckTriggerAOD(fAOD)) break;
      }
    }
  
    //I'll only cut on the reconstructed vertex since these are the events that will be used
    int vertexBin;
    if(!fAODData)vertexBin=CheckVertex(fESD);
    else vertexBin=CheckVertexAOD(fAOD);
    if(vertexBin<0)break;

 
    
    Int_t nGoodTracks[2]={0,0}, nTriggers[nTPtBins][nCentBins][2];
    Int_t nTrack;
    if(!ievent){
      if(!fAODData)nTrack=fESD->GetNumberOfTracks();
      else nTrack=fAOD->GetNumberOfTracks();
    }
    else nTrack=fMC->Stack()->GetNtrack();
    if(fSimulate)nTrack=10*(fSimNBgPart+10*fSimNJetPart);
    gRandom->SetSeed(time(0)+gSystem->GetPid());
    Float_t tdPhi, tdEta, tXE;
    Float_t tdPhi2, tdEta2;
    Float_t V2_T1, V2_T2, V3_T1, V3_T2, V4_T1, V4_T2, V42_T12, V42_1T2, V42_2T1, V2V2V4, V2, V3, V4;
    Float_t V2_A, V3_A, V4_A;
    ftPhi=new Float_t [nTrack];
    ftEta=new Float_t [nTrack];
    ftPt=new Float_t [nTrack];
    ftCharge=new Short_t [nTrack];
    ftEff=new Float_t [nTrack];
    ftV2=new Float_t [nTrack];
    ftV3=new Float_t [nTrack];
    ftV4=new Float_t [nTrack];
    ftPtAssoc3=new Int_t *[nTrack];
    for(int i=0;i<nTrack;i++){
      ftPtAssoc3[i]=new Int_t [10];
    }
    ftNPtAssoc3=new Int_t [nTrack];
    Short_t sign;
    //trigger particle arrays
    for(int i=0;i<fNTPtBins;i++){
      for(int c=0;c<fNCentBins;c++){
	nTriggers[i][c][ievent]=0;
      }
    }
    //Int_t tMult=fESD->GetMultiplicity()->GetNumberOfTracklets();//I think this is the correct multiplicity to use
    
    //AliESDVZERO* esdV0 = fESD->GetVZEROData();
    Float_t tMult=0;
    if(!fAODData){
      if(fCentPercent) tMult=fESD->GetCentrality()->GetCentralityPercentile("V0M");
      else tMult=fESD->GetVZEROData()->GetMTotV0A()+fESD->GetVZEROData()->GetMTotV0C();
    }
    else{
      AliAODHeader *tHeader=fAOD->GetHeader();
      tMult=tHeader->GetCentrality();
    }

    if(fDEBUG)Printf("Mult/Cent%6.1f",tMult);
  
    //Decide what multiplicy bins are filled with this event, note set to max of 4 overlapping bins as I didn't think more then 2 overlapping bins likely at one time, easliy changed
    Int_t multArray[4]={0,0,0,0};
    Int_t maxArray=0;
    
    for(int imult=0;imult<fNCentBins;imult++){
      if(tMult>=fCentArrayMin[imult]&&tMult<fCentArrayMax[imult]){
	multArray[maxArray]=imult;
	maxArray++;
      }
    }
    if(maxArray==0)break;
    //Printf("maxArray%d Mult%1.2f",maxArray,tMult);
    if(fDEBUG)Printf("maxArray%d",maxArray);
    //Set Efficiency and flow for the centrality bin (lowest bin used in array if multiple overlap)
    for(int ipar=0;ipar<fNFitLowParam;ipar++){
      fFitLow->SetParameter(ipar,fFitLowParam[multArray[0]*fNFitLowParam+ipar]);
    }
    for(int ipar=0;ipar<fNFitHighParam;ipar++){
      fFitHigh->SetParameter(ipar,fFitHighParam[multArray[0]*fNFitHighParam+ipar]);
    }
    for(int ipar=0;ipar<fNFitLowParamV2;ipar++){
      fFitLowV2->SetParameter(ipar,fFitLowParamV2[multArray[0]*fNFitLowParamV2+ipar]);
    }
    for(int ipar=0;ipar<fNFitHighParamV2;ipar++){
      fFitHighV2->SetParameter(ipar,fFitHighParamV2[multArray[0]*fNFitHighParamV2+ipar]);
    }
    for(int ipar=0;ipar<fNFitLowParamV3;ipar++){
      fFitLowV3->SetParameter(ipar,fFitLowParamV3[multArray[0]*fNFitLowParamV3+ipar]);
    }
    for(int ipar=0;ipar<fNFitHighParamV3;ipar++){
      fFitHighV3->SetParameter(ipar,fFitHighParamV3[multArray[0]*fNFitHighParamV2+ipar]);
    }
    for(int ipar=0;ipar<fNFitLowParamV4;ipar++){
      fFitLowV4->SetParameter(ipar,fFitLowParamV4[multArray[0]*fNFitLowParamV4+ipar]);
    }
    for(int ipar=0;ipar<fNFitHighParamV4;ipar++){
      fFitHighV4->SetParameter(ipar,fFitHighParamV4[multArray[0]*fNFitHighParamV4+ipar]);
    }
    fHistMult[ievent]->Fill(tMult);
    for(int c=0;c<maxArray;c++){fHistNEvents[multArray[c]][ievent]->Fill(0);}//count the number of events used
    Int_t leadPart=-1;
    
    //returns arrays filled up to nGoodTracks with tracks passing cuts
    for(int nSimEvents=0;nSimEvents<=(fSimulate*fSimNEvents);nSimEvents++){//only 1 loop if not simulation
      //Printf("nSimEvents %d",nSimEvents);
      if(fSimulate)leadPart=TrackCutsSim(ftPt,ftEta,ftPhi,ftCharge,ftEff,ftV2,ftV3,ftV4,ftPtAssoc3,ftNPtAssoc3,nGoodTracks);
      else if(!ievent){
	if(!fAODData)leadPart=TrackCuts(fESD,ftPt,ftEta,ftPhi,ftCharge,ftEff,ftV2,ftV3,ftV4,ftPtAssoc3,ftNPtAssoc3,nGoodTracks);
	else leadPart=TrackCutsAOD(fAOD,ftPt,ftEta,ftPhi,ftCharge,ftEff,ftV2,ftV3,ftV4,ftPtAssoc3,ftNPtAssoc3,nGoodTracks);
      }
      else leadPart=TrackCutsMC(fMC,ftPt,ftEta,ftPhi,ftCharge,ftEff,ftV2,ftV3,ftV4,ftPtAssoc3,ftNPtAssoc3,nGoodTracks);
      //Printf("nGoodTracks %d",nGoodTracks[0]);
      int nearEta=0,NearXE=0;
      int nearEta2=0;
 
      if(fDEBUG)Printf("Track Loop");
      for(int iTrack=0;iTrack<nGoodTracks[ievent];iTrack++){
	if(fDEBUG)Printf("Track%d Pt%f",iTrack,ftPt[iTrack]);
	//if(ftPhi[iTrack]<fdPhiMin)ftPhi[iTrack]+=2*fPi;
	//if(ftPhi[iTrack]>fdPhiMax)ftPhi[iTrack]-=2*fPi;
	for(int c=0;c<maxArray;c++){
	  // Printf("c%d mult%d",c,multArray[c]);
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
	  if(ftPt[iTrack]>fPtTrigArray[i]&&ftPt[iTrack]<=fPtTrigArray[i+1]&&fabs(ftEta[iTrack])<fTrigEtaCut){
	    if(fDEBUG)Printf("In %fpt%f",fPtTrigArray[i],fPtTrigArray[i+1]);
	    fHistMultTrig[i][ievent]->Fill(tMult);
	    for(int c=0;c<maxArray;c++){
	      nTriggers[i][multArray[c]][ievent]++;
	      fHistNTrigger[multArray[c]][ievent]->Fill(i,ftEff[iTrack]);
	      fHistNTriggerPt[multArray[c]][ievent]->Fill(i,ftPt[iTrack]*ftEff[iTrack]);
	    }
	  
	    if(fDEBUG)Printf("Assiciated Particle Loop");
	    //	Printf("GoodTracks %d Cent%d",nGoodTracks[ievent],multArray[0]);
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
		fHistPtTrig[i][multArray[c]][ievent]->Fill(ftPt[iTrack2],ftEff[iTrack2]*ftEff[iTrack]);
		fHistPhiTrig[i][multArray[c]][ievent]->Fill(ftPhi[iTrack2],ftPt[iTrack2],ftEff[iTrack2]*ftEff[iTrack]);
		fHistPhiTrigPt[i][multArray[c]][ievent]->Fill(ftPhi[iTrack2],ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]*ftEff[iTrack]);
		fHistEtaTrig[i][multArray[c]][ievent]->Fill(ftEta[iTrack2],ftPt[iTrack2],ftEff[iTrack2]*ftEff[iTrack]);
		fHistEtaTrigPt[i][multArray[c]][ievent]->Fill(ftEta[iTrack2],ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]*ftEff[iTrack]);
	      
		fHistPhiEtaTrig[i][multArray[c]][ievent]->Fill(ftPhi[iTrack2],ftEta[iTrack2],ftPt[iTrack2],ftEff[iTrack2]*ftEff[iTrack]);
		fHistPhiEtaTrigPt[i][multArray[c]][ievent]->Fill(ftPhi[iTrack2],ftEta[iTrack2],ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]*ftEff[iTrack]);
		fHistDeltaPhi[i][multArray[c]][0][ievent]->Fill(tdPhi,ftPt[iTrack2],ftEff[iTrack2]*ftEff[iTrack]);
		fHistDeltaPhiPt[i][multArray[c]][0][ievent]->Fill(tdPhi,ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]*ftEff[iTrack]);
		fHistDeltaPhi[i][multArray[c]][sign][ievent]->Fill(tdPhi,ftPt[iTrack2],ftEff[iTrack2]*ftEff[iTrack]);
		fHistDeltaPhiPt[i][multArray[c]][sign][ievent]->Fill(tdPhi,ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]*ftEff[iTrack]);
	      
		if(nearEta){
		  fHistDeltaEtaN[i][multArray[c]][0][ievent]->Fill(tdEta,ftPt[iTrack2],ftEff[iTrack2]*ftEff[iTrack]);
		  fHistDeltaEtaNPt[i][multArray[c]][0][ievent]->Fill(tdEta,ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]*ftEff[iTrack]);
		  fHistDeltaEtaN[i][multArray[c]][sign][ievent]->Fill(tdEta,ftPt[iTrack2],ftEff[iTrack2]*ftEff[iTrack]);
		  fHistDeltaEtaNPt[i][multArray[c]][sign][ievent]->Fill(tdEta,ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]*ftEff[iTrack]);
		}
		else{
		  fHistDeltaEtaA[i][multArray[c]][0][ievent]->Fill(tdEta,ftPt[iTrack2],ftEff[iTrack2]*ftEff[iTrack]);
		  fHistDeltaEtaAPt[i][multArray[c]][0][ievent]->Fill(tdEta,ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]*ftEff[iTrack]);
		  fHistDeltaEtaA[i][multArray[c]][sign][ievent]->Fill(tdEta,ftPt[iTrack2],ftEff[iTrack2]*ftEff[iTrack]);
		  fHistDeltaEtaAPt[i][multArray[c]][sign][ievent]->Fill(tdEta,ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]*ftEff[iTrack]);
		}
		fHistDeltaPhiEta[i][multArray[c]][ievent]->Fill(tdPhi,tdEta,ftPt[iTrack2],ftEff[iTrack2]*ftEff[iTrack]);
		fHistDeltaPhiEtaPt[i][multArray[c]][ievent]->Fill(tdPhi,tdEta,ftPt[iTrack2],ftPt[iTrack2]*ftEff[iTrack2]*ftEff[iTrack]);
	      
		//only fill these if trigger particle is the leading particle
		if(iTrack==leadPart){
		  if(NearXE){
		    tXE=ftPt[iTrack2]*cos(tdPhi)/ftPt[iTrack];
		    fHistXEN[i][multArray[c]][ievent]->Fill(tXE,ftEff[iTrack2]*ftEff[iTrack]);
		  }
		  else{
		    tXE=ftPt[iTrack2]*cos(tdPhi+fPi)/ftPt[iTrack];
		    fHistXEA[i][multArray[c]][ievent]->Fill(tXE,ftEff[iTrack2]*ftEff[iTrack]);
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
			fHistDeltaPhiPhi[i][ftPtAssoc3[iTrack2][e]][multArray[c]][0][ievent]->Fill(tdPhi,tdPhi2,ftEff[iTrack2]*ftEff[iTrack3]*ftEff[iTrack]);
			fHistDeltaPhiPhi[i][ftPtAssoc3[iTrack2][e]][multArray[c]][sign][ievent]->Fill(tdPhi2,tdPhi,ftEff[iTrack2]*ftEff[iTrack3]*ftEff[iTrack]);
		      
		      
			if(nearEta2){
			  fHistDeltaEtaEta[i][ftPtAssoc3[iTrack2][e]][multArray[c]][0][ievent]->Fill(tdEta,tdEta2,ftEff[iTrack2]*ftEff[iTrack3]*ftEff[iTrack]);
			  fHistDeltaEtaEta[i][ftPtAssoc3[iTrack2][e]][multArray[c]][sign][ievent]->Fill(tdEta,tdEta2,ftEff[iTrack2]*ftEff[iTrack3]*ftEff[iTrack]);
			}
		      }//multiplicity loop (c)
		    }
		  }
		}//track checking loops
	      }//iTrack3
	    }//iTrack2 (associated track loop)
	  
	    if(fDEBUG)Printf("Mixed Event Loop");
	    for(int c=0;c<maxArray;c++){
	      //Printf("c%d mult%d",c,multArray[c]);
	      int d=multArray[c];//Centrality bin we are in
	      if(fMixEnd[d][vertexBin][ievent]>2){//check if there are any mixed events for this bin, require 2 for soft-soft mixing
		for(int imix=0;imix<fMixEnd[d][vertexBin][ievent];imix++){//loop over the stored mixed events
		  fHistNMix[d][ievent]->Fill(i);
		  //Printf("GoodTracksMixed %d Cent%d fMixEnd%d",fMixTrack[imix][d][vertexBin][ievent],d,fMixEnd[d][vertexBin][ievent]);
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
		    V2_T1=2*fMV2[imix][d][vertexBin][ievent][iTrack2]*ftV2[iTrack]*cos(2*tdPhi);
		    V3_T1=2*fMV3[imix][d][vertexBin][ievent][iTrack2]*ftV3[iTrack]*cos(3*tdPhi);
		    V4_T1=2*fMV4[imix][d][vertexBin][ievent][iTrack2]*ftV4[iTrack]*cos(4*tdPhi);
		  
		    fHistDeltaPhiMix[i][d][0][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]);
		    fHistDeltaPhiMixV2[i][d][0][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V2_T1);
		    fHistDeltaPhiMixV3[i][d][0][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V3_T1);
		    fHistDeltaPhiMixV4[i][d][0][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V4_T1);
		    fHistDeltaPhiMixPt[i][d][0][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMPt[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack2]);
		    fHistDeltaPhiMixV2Pt[i][d][0][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V2_T1*fMPt[imix][d][vertexBin][ievent][iTrack2]);
		    fHistDeltaPhiMixV3Pt[i][d][0][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V3_T1*fMPt[imix][d][vertexBin][ievent][iTrack2]);
		    fHistDeltaPhiMixV4Pt[i][d][0][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V4_T1*fMPt[imix][d][vertexBin][ievent][iTrack2]);
		    fHistDeltaPhiMix[i][d][sign][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]);
		    fHistDeltaPhiMixV2[i][d][sign][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V2_T1);
		    fHistDeltaPhiMixV3[i][d][sign][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V3_T1);
		    fHistDeltaPhiMixV4[i][d][sign][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V4_T1);
		    fHistDeltaPhiMixPt[i][d][sign][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMPt[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack2]);
		    fHistDeltaPhiMixV2Pt[i][d][sign][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V2_T1*fMPt[imix][d][vertexBin][ievent][iTrack2]);
		    fHistDeltaPhiMixV3Pt[i][d][sign][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V3_T1*fMPt[imix][d][vertexBin][ievent][iTrack2]);
		    fHistDeltaPhiMixV4Pt[i][d][sign][ievent]->Fill(tdPhi,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V4_T1*fMPt[imix][d][vertexBin][ievent][iTrack2]);
		  
		    if(nearEta){
		      fHistDeltaEtaNMix[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]);
		      fHistDeltaEtaNMixV2[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V2_T1);
		      fHistDeltaEtaNMixV3[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V3_T1);
		      fHistDeltaEtaNMixV4[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V4_T1);
		      fHistDeltaEtaNMixPt[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMPt[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack2]);
		      fHistDeltaEtaNMix[i][d][sign][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]);
		      fHistDeltaEtaNMixPt[i][d][sign][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMPt[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack2]);
		      fHistDeltaEtaNMixV2Pt[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V2_T1*fMPt[imix][d][vertexBin][ievent][iTrack2]);
		      fHistDeltaEtaNMixV3Pt[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V3_T1*fMPt[imix][d][vertexBin][ievent][iTrack2]);
		      fHistDeltaEtaNMixV4Pt[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V4_T1*fMPt[imix][d][vertexBin][ievent][iTrack2]);
		    }
		    else{
		      fHistDeltaEtaAMix[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]);
		      fHistDeltaEtaAMixV2[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V2_T1);
		      fHistDeltaEtaAMixV3[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V3_T1);
		      fHistDeltaEtaAMixV4[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V4_T1);
		      fHistDeltaEtaAMixPt[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMPt[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack2]);
		      fHistDeltaEtaAMix[i][d][sign][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]);
		      fHistDeltaEtaAMixPt[i][d][sign][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMPt[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix][d][vertexBin][ievent][iTrack2]);
		      fHistDeltaEtaAMixV2Pt[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V2_T1*fMPt[imix][d][vertexBin][ievent][iTrack2]);
		      fHistDeltaEtaAMixV3Pt[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V3_T1*fMPt[imix][d][vertexBin][ievent][iTrack2]);
		      fHistDeltaEtaAMixV4Pt[i][d][0][ievent]->Fill(tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V4_T1*fMPt[imix][d][vertexBin][ievent][iTrack2]);
		    }	
		  
		    fHistDeltaPhiEtaMix[i][d][ievent]->Fill(tdPhi,tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]);
		    fHistDeltaPhiEtaMixV2[i][d][ievent]->Fill(tdPhi,tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V2_T1);
		    fHistDeltaPhiEtaMixV3[i][d][ievent]->Fill(tdPhi,tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V3_T1);
		    fHistDeltaPhiEtaMixV4[i][d][ievent]->Fill(tdPhi,tdEta,fMPt[imix][d][vertexBin][ievent][iTrack2],fMEff[imix][d][vertexBin][ievent][iTrack2]*V4_T1);
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
		    int imix2=imix+1;
		    if(imix2>=fMixEnd[d][vertexBin][ievent])imix2=0;
		    if(imix2==imix)continue;
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
		      V2_T2=2*ftV2[iTrack]*fMV2[imix2][d][vertexBin][ievent][iTrack3]*cos(2*tdPhi2);
		      V3_T2=2*ftV3[iTrack]*fMV3[imix2][d][vertexBin][ievent][iTrack3]*cos(3*tdPhi2);
		      V4_T2=2*ftV4[iTrack]*fMV4[imix2][d][vertexBin][ievent][iTrack3]*cos(4*tdPhi2);
		      V42_T12=2*ftV4[iTrack]*fMV2[imix2][d][vertexBin][ievent][iTrack3]*fMV2[imix][d][vertexBin][ievent][iTrack2]*cos(2*tdPhi+2*tdPhi2);
		      V42_1T2=2*fMV4[imix][d][vertexBin][ievent][iTrack2]*ftV2[iTrack]*fMV2[imix2][d][vertexBin][ievent][iTrack3]*cos(4*tdPhi-2*tdPhi2);
		      V42_2T1=2*fMV4[imix2][d][vertexBin][ievent][iTrack3]*ftV2[iTrack]*fMV2[imix][d][vertexBin][ievent][iTrack2]*cos(4*tdPhi2-2*tdPhi);
		      V2_A=2*fMV2[imix2][d][vertexBin][ievent][iTrack3]*fMV2[imix][d][vertexBin][ievent][iTrack2]*cos(2*(tdPhi-tdPhi2));
		      V3_A=2*fMV3[imix2][d][vertexBin][ievent][iTrack3]*fMV3[imix][d][vertexBin][ievent][iTrack2]*cos(3*(tdPhi-tdPhi2));
		      V4_A=2*fMV4[imix2][d][vertexBin][ievent][iTrack3]*fMV4[imix][d][vertexBin][ievent][iTrack2]*cos(4*(tdPhi-tdPhi2)); 
		    
		      V2=V2_T1+V2_T2;
		      V3=V3_T1+V3_T2;
		      V4=V4_T1+V4_T2;
		      V2V2V4=V42_T12+V42_1T2+V42_2T1;
		    
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
			    //v2
			    fHistDeltaPhiPhiMixV2[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2);
			    fHistDeltaPhiPhiMixV2[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2); 
			    fHistDeltaPhiPhiMixV2[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2);//free factor of 2 in statistics
			    fHistDeltaPhiPhiMixV2[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2);
			    fHistDeltaPhiPhiSSV2[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2_A);
			    fHistDeltaPhiPhiSSV2[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2_A);
			    //v3
			    fHistDeltaPhiPhiMixV3[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V3);
			    fHistDeltaPhiPhiMixV3[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V3); 
			    fHistDeltaPhiPhiMixV3[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V3);//free factor of 2 in statistics
			    fHistDeltaPhiPhiMixV3[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V3);
			    fHistDeltaPhiPhiSSV3[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V3_A);
			    fHistDeltaPhiPhiSSV3[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V3_A);
			    //v4
			    fHistDeltaPhiPhiMixV4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V4);
			    fHistDeltaPhiPhiMixV4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V4); 
			    fHistDeltaPhiPhiMixV4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V4);//free factor of 2 in statistics
			    fHistDeltaPhiPhiMixV4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V4);
			    fHistDeltaPhiPhiSSV4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V4_A);
			    fHistDeltaPhiPhiSSV4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V4_A);
			    //v2v2v4
			    fHistDeltaPhiPhiMixV2V2V4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2V2V4);
			    fHistDeltaPhiPhiMixV2V2V4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdPhi,tdPhi2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2V2V4); 
			    fHistDeltaPhiPhiMixV2V2V4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2V2V4);//free factor of 2 in statistics
			    fHistDeltaPhiPhiMixV2V2V4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdPhi2,tdPhi,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2V2V4);
			  
			    if(nearEta2){
			      fHistDeltaEtaEtaMix[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]);
			      fHistDeltaEtaEtaMix[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]);
			      fHistDeltaEtaEtaMix[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta2,tdEta,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]);
			      fHistDeltaEtaEtaMix[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdEta2,tdEta,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]);
			      //v2
			      fHistDeltaEtaEtaMixV2[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2);
			      fHistDeltaEtaEtaMixV2[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2);
			      fHistDeltaEtaEtaMixV2[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta2,tdEta,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2);
			      fHistDeltaEtaEtaMixV2[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdEta2,tdEta,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2);
			      //v3
			      fHistDeltaEtaEtaMixV3[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V3);
			      fHistDeltaEtaEtaMixV3[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V3);
			      fHistDeltaEtaEtaMixV3[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta2,tdEta,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V3);
			      fHistDeltaEtaEtaMixV3[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdEta2,tdEta,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V3);
			      //v4
			      fHistDeltaEtaEtaMixV4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V4);
			      fHistDeltaEtaEtaMixV4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V4);
			      fHistDeltaEtaEtaMixV4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta2,tdEta,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V4);
			      fHistDeltaEtaEtaMixV4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdEta2,tdEta,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V4);
			      //v2v2v4	
			      fHistDeltaEtaEtaMixV2V2V4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2V2V4);
			      fHistDeltaEtaEtaMixV2V2V4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdEta,tdEta2,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2V2V4);
			      fHistDeltaEtaEtaMixV2V2V4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][0][ievent]->Fill(tdEta2,tdEta,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2V2V4);
			      fHistDeltaEtaEtaMixV2V2V4[i][fMPtAssoc3[imix][d][vertexBin][ievent][e][iTrack2]][d][sign][ievent]->Fill(tdEta2,tdEta,fMEff[imix][d][vertexBin][ievent][iTrack2]*fMEff[imix2][d][vertexBin][ievent][iTrack3]*V2V2V4);
			
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
  
      //now store this event for mixing (using these dynamic arrays to save memory)
      if(fDEBUG)Printf("Store Event For Mixing");
      for(int c=0;c<maxArray;c++){//loops over centrality bins
	int d=multArray[c];//too many nested arrays looked confusing d=which centrality bin
	if(fMixEnd[d][vertexBin][ievent]<=fNMix)fMixEnd[d][vertexBin][ievent]++;
	if(fMixPointer[d][vertexBin][ievent]<(fNMix-1)&&fMixEnd[d][vertexBin][ievent]!=1)fMixPointer[d][vertexBin][ievent]++;
	else fMixPointer[d][vertexBin][ievent]=0;
	int e=fMixPointer[d][vertexBin][ievent];//nested arrays (e is event number in pool)
	delete [] fMPt[e][d][vertexBin][ievent];
	delete [] fMPhi[e][d][vertexBin][ievent];
	delete [] fMEta[e][d][vertexBin][ievent];
	delete [] fMCharge[e][d][vertexBin][ievent];
	delete [] fMEff[e][d][vertexBin][ievent];
	delete [] fMV2[e][d][vertexBin][ievent];
	delete [] fMV3[e][d][vertexBin][ievent];
	delete [] fMV4[e][d][vertexBin][ievent];
	delete [] fMNPtAssoc3[e][d][vertexBin][ievent];
	for(int jj=0;jj<10;jj++){
	  delete [] fMPtAssoc3[e][d][vertexBin][ievent][jj];
	}
	fMPt[e][d][vertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
	fMPhi[e][d][vertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
	fMEta[e][d][vertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
	fMCharge[e][d][vertexBin][ievent]=new Short_t [nGoodTracks[ievent]];
	fMEff[e][d][vertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
	fMV2[e][d][vertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
	fMV3[e][d][vertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
	fMV4[e][d][vertexBin][ievent]=new Float_t [nGoodTracks[ievent]];
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
	  fMV2[e][d][vertexBin][ievent][iTrack]=ftV2[iTrack];
	  fMV3[e][d][vertexBin][ievent][iTrack]=ftV3[iTrack];
	  fMV4[e][d][vertexBin][ievent][iTrack]=ftV4[iTrack];
	  fMNPtAssoc3[e][d][vertexBin][ievent][iTrack]=ftNPtAssoc3[iTrack];
	  for(int jj=0;jj<ftNPtAssoc3[iTrack];jj++){
	    fMPtAssoc3[e][d][vertexBin][ievent][jj][iTrack]=ftPtAssoc3[iTrack][jj];
	  }
	}//iTracks
      }//Centrality (c)
    } //sim
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
  delete [] ftV2;
  delete [] ftV3;
  delete [] ftV4;
  delete [] ftNPtAssoc3;
  delete [] ftPtAssoc3;
  ftPhi=NULL;
  ftEta=NULL;
  ftPt=NULL;
  ftCharge=NULL;
  ftEff=NULL;
  ftV2=NULL;
  ftV3=NULL;
  ftV4=NULL;
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
	  delete [] fMV2[ii][cc][vtx][jj];
	  delete [] fMV3[ii][cc][vtx][jj];
	  delete [] fMV4[ii][cc][vtx][jj];
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
	  fMV2[ii][cc][vtx][jj]=NULL;
	  fMV3[ii][cc][vtx][jj]=NULL;
	  fMV4[ii][cc][vtx][jj]=NULL;
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
