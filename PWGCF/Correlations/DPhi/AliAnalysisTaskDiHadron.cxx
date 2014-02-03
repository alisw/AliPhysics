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
AliAnalysisTaskSE(name), fESD(0), fAOD(0), fMC(0), fOutput(0),fMinClustersTPC(0),fMinClusterRatio(0),fMaxTPCchi2(0),fMinClustersITS(0),fEtaCut(0),fTrigEtaCut(0),fNearPhiCut(0),fXECut(0),fMaxDCA(0),fMaxDCAXY(0),fMaxDCAZ(0),fDCA2D(0),fTPCRefit(0),fITSRefit(0),fSPDCut(0),fMinPtAssoc(0),fMaxPtAssoc(0),fVzCut(0),fAODData(0),fEfficiencyCorr(0),fDEBUG(0),fnBinPhi(0),fnBinEta(0),fnBinPhiEtaPhi(0),fnBinPhiEtaEta(0),fnBinPhi3(0),fnBinEta3(0),fPi(3.1415926535898),fdPhiMin(0),fdPhiMax(0),fNTPtBins(0),fNMix(0),fNCentBins(0),fCentPercent(0),fNAPtBins(0),fNAPt3Bins(0),fNVertexBins(0),fNXEBins(0),fNIDs(0),fEffFitPt(0),fNFitLowParam(0),fNFitHighParam(0),fV2FitPt(0),fV3FitPt(0),fV4FitPt(0),fNFitLowParamV2(0),fNFitHighParamV2(0),fNFitLowParamV3(0),fNFitHighParamV3(0),fNFitLowParamV4(0),fNFitHighParamV4(0),fMCHistos(0),fFitLow(NULL),fFitHigh(NULL),fFitLowParam(NULL),fFitHighParam(NULL),fFitLowV2(NULL),fFitHighV2(NULL),fFitLowParamV2(NULL),fFitHighParamV2(NULL),fFitLowV3(NULL),fFitHighV3(NULL),fFitLowParamV3(NULL),fFitHighParamV3(NULL),fFitLowV4(NULL),fFitHighV4(NULL),fFitLowParamV4(NULL),fFitHighParamV4(NULL),fPtTrigArray(NULL),fPtAssocArray(NULL),fPtAssoc3Array1(NULL),fPtAssoc3Array2(NULL),fCentArrayMin(NULL),fCentArrayMax(NULL),fXEArray(NULL),fTrigIDArray(NULL),fSimulate(0),fSimNBgPart(0),fSimNJetPart(0),fSimNJet(0),fSimNEvents(0),fSimAwayDeflected(0),fSimPsi2(0),fSimPsi3(0),fSimPsi4(0),fSimFlowMark(0),ftPhi(NULL),ftEta(NULL),ftPt(NULL),ftCharge(NULL),ftEff(NULL),ftV2(NULL),ftV3(NULL),ftV4(NULL),ftPtAssoc3(NULL),ftNPtAssoc3(NULL)
  
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
void AliAnalysisTaskDiHadron::SetSimulation(Int_t Simulate, Float_t SimNBgPart, Float_t SimNJetPart, Float_t SimNJet, Int_t SimNEvents, Int_t SimAwayDeflected){
  fSimulate=Simulate;
  fSimNBgPart=SimNBgPart;
  fSimNJetPart=SimNJetPart;
  fSimNJet=SimNJet;
  fSimNEvents=SimNEvents;
  fSimAwayDeflected=SimAwayDeflected;
  fSimFlowMark=0;
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
  char histname[512]={0};
  char histtitle[512]={0};
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
  Int_t BufferSize=sizeof(histname);
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
    if((!eSDtrack->GetStatus()&AliESDtrack::kTPCrefit)&&fTPCRefit)continue;//refit in TPC
    if((fITSRefit==1||(fITSRefit==2&&sPt>5))&&(!eSDtrack->GetStatus()&AliESDtrack::kITSrefit))continue;//refit of its tracks either for none,all, or >5 GeV/c
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

void AliAnalysisTaskDiHadron::CalcFlow(Float_t *rPt, Float_t *rEta, Float_t *rPhi, Int_t *rGoodTracks, Int_t LeadPart){
  float psin2=0, psin3=0, psin4=0, pcos2=0, pcos3=0, pcos4=0;
  float cos2=0, cos3=0, cos4=0, sum=0;
  float rNAssoc=0, rNTrig=0;
  for(int i=0;i<rGoodTracks[0];i++){
    if(i==LeadPart)continue;
    if(rPt[i]>1&&rPt[i]<2)rNAssoc++;
    if(rPt[i]>3&&rPt[i]<4)rNTrig++;
    if(fabs(rEta[i]-rEta[LeadPart])<0.5) continue;
    if(rPt[i]<2){
      psin2+=rPt[i]*sin(2*rPhi[i]);
      pcos2+=rPt[i]*cos(2*rPhi[i]);
      psin3+=rPt[i]*sin(3*rPhi[i]);
      pcos3+=rPt[i]*cos(3*rPhi[i]);
      psin4+=rPt[i]*sin(4*rPhi[i]);
      pcos4+=rPt[i]*cos(4*rPhi[i]);
   
    }
    if(rPt[i]>3&&rPt[i]<4){
      for(int j=0;j<rGoodTracks[0];j++){
	if(rPt[j]<2&&rPt[j]>1){
	  cos2+=cos(2*rPhi[i]-2*rPhi[j]);
	  cos3+=cos(3*rPhi[i]-3*rPhi[j]);
	  cos4+=cos(4*rPhi[i]-4*rPhi[j]);
	  sum++;
	}//assoc pt
      }//assoc loop
    }//trigger pt
  }//loop

  if(sum==0)sum=1E-10;
  fSimPsi2=atan(psin2/pcos2)/2;
  fSimPsi3=atan(psin3/pcos3)/3;
  fSimPsi4=atan(psin4/pcos4)/4;
  fFitHighV2->SetParameter(0,sqrt(fabs(cos2)/sum));
  fFitHighV3->SetParameter(0,sqrt(fabs(cos3)/sum));
  fFitHighV4->SetParameter(0,sqrt(fabs(cos4)/sum));
  fSimNBgPart=rNAssoc;
  fSimNJet=rNTrig;

}


//---------------------------------------------------------
Int_t AliAnalysisTaskDiHadron::TrackCutsSim(Float_t *rPt, Float_t *rEta, Float_t *rPhi, Short_t *rCharge, Float_t *rEff, Float_t *rV2, Float_t *rV3, Float_t *rV4, Int_t **rPtAssoc3, Int_t *rNPtAssoc3, Int_t *rGoodTracks){

  Float_t FlowValues[500][8]={{0.0035, 0.0136, 0.0105, 0.0010, 1.1650, 1.0527, 0.1205, 0.8645},
			    {0.0030, 0.0176, 0.0105, 0.0014, 0.3562, 0.6850, 0.0579, 0.8537},
			    {0.0047, 0.0153, 0.0062, 0.0047, 1.6184, 0.3231, 1.0531, 0.1074},
			    {0.0095, 0.0201, 0.0096, 0.0036, 1.1438, 0.7005, 1.4932, 0.6805},
			    {0.0147, 0.0187, 0.0090, 0.0041, 0.3040, 0.7596, 1.2148, 0.3205},
			    {0.0085, 0.0145, 0.0058, 0.0006, 2.8200, 0.6537, 1.4212, 0.4427},
			    {0.0027, 0.0098, 0.0067, 0.0010, 0.7530, 1.7227, 0.8277, 0.2671},
			    {0.0211, 0.0141, 0.0040, 0.0032, 2.8158, 1.2640, 1.5320, 0.8592},
			    {0.0214, 0.0034, 0.0082, 0.0028, 0.9113, 1.3545, 0.8183, 0.4192},
			    {0.0101, 0.0178, 0.0093, 0.0023, 0.1782, 0.4587, 1.3575, 0.4537},
			    {0.0077, 0.0030, 0.0107, 0.0018, 3.0208, 1.2572, 1.3529, 1.0689},
			    {0.0156, 0.0217, 0.0069, 0.0011, 2.7169, 0.9773, 0.2563, 0.5131},
			    {0.0139, 0.0060, 0.0009, 0.0009, 1.3356, 1.0727, 0.8756, 0.0257},
			    {0.0187, 0.0250, 0.0067, 0.0040, 1.2627, 1.9085, 0.7025, 1.0722},
			    {0.0052, 0.0196, 0.0104, 0.0050, 0.3652, 0.8782, 0.8404, 1.0192},
			    {0.0079, 0.0095, 0.0061, 0.0037, 0.3816, 0.6331, 0.1517, 0.4630},
			    {0.0060, 0.0169, 0.0042, 0.0012, 3.0543, 0.7154, 0.1500, 0.8865},
			    {0.0022, 0.0220, 0.0080, 0.0031, 1.5914, 0.4051, 0.5805, 0.5541},
			    {0.0010, 0.0055, 0.0184, 0.0044, 1.3026, 1.9182, 1.5478, 0.2334},
			    {0.0128, 0.0089, 0.0105, 0.0007, 1.6267, 0.0414, 1.3179, 0.8428},
			    {0.0078, 0.0258, 0.0082, 0.0025, 1.3202, 1.4525, 0.8713, 0.1477},
			    {0.0105, 0.0163, 0.0098, 0.0042, 0.0085, 0.0798, 0.3661, 0.6749},
			    {0.0098, 0.0186, 0.0117, 0.0018, 0.6094, 0.7213, 1.2876, 0.4772},
			    {0.0074, 0.0071, 0.0036, 0.0046, 0.7562, 0.5789, 0.2272, 1.0092},
			    {0.0111, 0.0040, 0.0094, 0.0026, 1.3796, 0.5633, 1.3583, 0.5369},
			    {0.0151, 0.0224, 0.0057, 0.0031, 0.1665, 0.9981, 0.1623, 1.1031},
			    {0.0110, 0.0314, 0.0071, 0.0009, 1.0761, 1.2344, 0.3076, 0.0696},
			    {0.0123, 0.0118, 0.0037, 0.0046, 1.7767, 1.3157, 1.5148, 1.1963},
			    {0.0081, 0.0068, 0.0087, 0.0023, 0.6168, 0.3762, 0.0063, 0.3572},
			    {0.0099, 0.0113, 0.0087, 0.0050, 0.5912, 2.0366, 1.3408, 1.2208},
			    {0.0135, 0.0052, 0.0091, 0.0031, 0.8365, 0.9733, 1.2951, 0.9458},
			    {0.0052, 0.0256, 0.0041, 0.0028, 1.6871, 0.7141, 0.2864, 0.4447},
			    {0.0026, 0.0109, 0.0046, 0.0040, 1.8569, 0.8634, 0.3792, 0.3380},
			    {0.0048, 0.0032, 0.0100, 0.0008, 0.2143, 1.1711, 0.4115, 1.0100},
			    {0.0190, 0.0082, 0.0011, 0.0017, 2.3545, 0.5168, 0.5337, 0.0927},
			    {0.0066, 0.0051, 0.0105, 0.0023, 1.4765, 1.0413, 0.6329, 0.5209},
			    {0.0071, 0.0041, 0.0138, 0.0011, 2.6051, 1.2165, 0.3467, 0.8268},
			    {0.0112, 0.0213, 0.0082, 0.0030, 2.4043, 0.1020, 0.2713, 1.2024},
			    {0.0088, 0.0135, 0.0061, 0.0043, 2.8336, 0.9174, 0.4703, 0.8523},
			    {0.0120, 0.0212, 0.0066, 0.0052, 1.8337, 0.7224, 0.0820, 0.9440},
			    {0.0032, 0.0065, 0.0035, 0.0012, 1.2218, 1.4102, 0.8424, 0.1840},
			    {0.0049, 0.0066, 0.0050, 0.0046, 0.5469, 1.0646, 1.3841, 0.6618},
			    {0.0122, 0.0079, 0.0027, 0.0017, 2.0599, 1.2457, 0.1777, 0.3349},
			    {0.0117, 0.0206, 0.0035, 0.0019, 1.1937, 1.6757, 1.5304, 0.2341},
			    {0.0073, 0.0180, 0.0071, 0.0017, 0.9287, 1.7990, 1.3380, 0.2897},
			    {0.0138, 0.0134, 0.0076, 0.0020, 1.4069, 0.7750, 0.7245, 0.1373},
			    {0.0032, 0.0214, 0.0068, 0.0036, 2.7958, 0.1488, 0.0752, 0.3032},
			    {0.0097, 0.0206, 0.0044, 0.0035, 2.9542, 0.6978, 0.3074, 0.8538},
			    {0.0014, 0.0176, 0.0128, 0.0045, 0.4687, 1.8134, 0.3071, 0.0202},
			    {0.0194, 0.0068, 0.0086, 0.0031, 0.0459, 1.0151, 0.8554, 0.9847},
			    {0.0136, 0.0121, 0.0086, 0.0007, 1.2042, 0.1901, 1.0211, 0.9832},
			    {0.0092, 0.0110, 0.0041, 0.0051, 2.3322, 1.6906, 0.6865, 0.8827},
			    {0.0099, 0.0176, 0.0113, 0.0016, 1.9565, 1.2077, 1.4654, 0.4928},
			    {0.0077, 0.0132, 0.0081, 0.0038, 2.2578, 1.3617, 0.8156, 1.0693},
			    {0.0053, 0.0342, 0.0238, 0.0037, 2.7501, 1.6721, 0.7829, 0.2254},
			    {0.0099, 0.0152, 0.0068, 0.0026, 1.2651, 1.4901, 1.3494, 1.2040},
			    {0.0047, 0.0084, 0.0092, 0.0024, 2.3898, 2.0839, 1.0090, 0.6047},
			    {0.0224, 0.0082, 0.0065, 0.0015, 1.7053, 0.5598, 0.6944, 1.2247},
			    {0.0245, 0.0209, 0.0024, 0.0013, 1.3785, 0.6178, 1.3101, 0.7325},
			    {0.0097, 0.0053, 0.0060, 0.0025, 1.8877, 1.5127, 0.6411, 0.0542},
			    {0.0107, 0.0120, 0.0108, 0.0021, 0.4713, 1.0744, 1.0183, 0.1119},
			    {0.0148, 0.0040, 0.0092, 0.0014, 1.8979, 1.8480, 1.0089, 1.1614},
			    {0.0132, 0.0082, 0.0091, 0.0023, 2.6065, 1.5960, 1.3907, 0.8666},
			    {0.0070, 0.0155, 0.0102, 0.0025, 0.5482, 1.9525, 0.5050, 0.8651},
			    {0.0147, 0.0086, 0.0093, 0.0015, 2.8907, 0.4051, 0.3219, 0.1822},
			    {0.0108, 0.0207, 0.0028, 0.0037, 1.0302, 1.5407, 0.5012, 0.7318},
			    {0.0126, 0.0047, 0.0079, 0.0028, 2.8788, 1.2861, 0.8888, 0.5333},
			    {0.0055, 0.0158, 0.0085, 0.0026, 2.2349, 1.9289, 1.0859, 0.0347},
			    {0.0109, 0.0086, 0.0062, 0.0038, 0.6779, 0.7527, 0.2708, 0.2838},
			    {0.0114, 0.0083, 0.0107, 0.0007, 1.5915, 1.1871, 1.5686, 0.5562},
			    {0.0040, 0.0102, 0.0074, 0.0024, 1.7429, 0.0685, 0.8020, 0.2097},
			    {0.0101, 0.0100, 0.0035, 0.0033, 1.2561, 0.5735, 0.9764, 0.4423},
			    {0.0145, 0.0160, 0.0092, 0.0057, 0.0503, 0.5607, 0.8943, 0.4559},
			    {0.0139, 0.0126, 0.0187, 0.0027, 0.7941, 0.7881, 1.3813, 0.1958},
			    {0.0211, 0.0038, 0.0025, 0.0028, 2.2224, 1.9979, 0.9060, 0.0495},
			    {0.0140, 0.0169, 0.0036, 0.0003, 0.5556, 1.6196, 0.0197, 0.2456},
			    {0.0133, 0.0334, 0.0097, 0.0004, 0.4103, 1.9339, 0.3124, 0.8308},
			    {0.0043, 0.0021, 0.0061, 0.0013, 2.5553, 1.8116, 0.7002, 0.2462},
			    {0.0071, 0.0183, 0.0083, 0.0016, 1.6071, 0.3276, 0.4626, 0.6513},
			    {0.0087, 0.0188, 0.0111, 0.0051, 2.4206, 0.1839, 1.4459, 0.1188},
			    {0.0118, 0.0019, 0.0077, 0.0039, 2.3648, 0.3677, 0.2608, 0.2471},
			    {0.0055, 0.0119, 0.0057, 0.0028, 2.7920, 2.0417, 1.2721, 1.1066},
			    {0.0021, 0.0112, 0.0055, 0.0045, 1.3043, 0.2018, 0.8319, 0.6204},
			    {0.0035, 0.0261, 0.0105, 0.0019, 1.2034, 2.0184, 0.0761, 0.0667},
			    {0.0077, 0.0135, 0.0101, 0.0039, 2.0438, 0.3174, 1.0312, 0.1953},
			    {0.0272, 0.0054, 0.0047, 0.0029, 2.8039, 1.4388, 1.4157, 0.7614},
			    {0.0049, 0.0085, 0.0018, 0.0028, 1.5170, 1.9910, 0.1080, 0.7917},
			    {0.0061, 0.0117, 0.0062, 0.0033, 0.2065, 0.6670, 1.2914, 1.1443},
			    {0.0196, 0.0069, 0.0041, 0.0008, 2.6236, 1.2682, 0.1914, 1.1040},
			    {0.0131, 0.0189, 0.0070, 0.0011, 0.1012, 0.0506, 0.3475, 0.3875},
			    {0.0080, 0.0058, 0.0086, 0.0025, 0.9333, 0.4077, 1.4301, 0.3130},
			    {0.0161, 0.0028, 0.0082, 0.0039, 1.3052, 0.2207, 0.6978, 0.7626},
			    {0.0151, 0.0092, 0.0095, 0.0019, 1.4551, 0.5386, 0.9054, 0.4826},
			    {0.0155, 0.0170, 0.0086, 0.0007, 0.9909, 1.9847, 0.2312, 1.1555},
			    {0.0060, 0.0073, 0.0063, 0.0030, 1.9213, 0.2065, 1.1014, 0.8468},
			    {0.0107, 0.0114, 0.0101, 0.0007, 0.1207, 0.2177, 1.4986, 0.7271},
			    {0.0252, 0.0091, 0.0045, 0.0027, 2.8788, 0.3174, 0.9813, 0.9926},
			    {0.0013, 0.0132, 0.0092, 0.0035, 0.4808, 0.4746, 0.6482, 0.0875},
			    {0.0046, 0.0115, 0.0075, 0.0031, 0.5798, 1.5405, 0.6005, 1.0992},
			    {0.0119, 0.0071, 0.0136, 0.0020, 2.3262, 2.0037, 0.1128, 0.2491},
			    {0.0094, 0.0087, 0.0083, 0.0017, 0.0152, 0.3550, 0.6614, 0.4152},
			    {0.0140, 0.0128, 0.0068, 0.0030, 2.9511, 1.9677, 1.5402, 0.4670},
			    {0.0086, 0.0119, 0.0073, 0.0017, 1.6429, 0.9199, 0.7814, 0.0791},
			    {0.0072, 0.0222, 0.0125, 0.0010, 1.9273, 1.6578, 1.4650, 0.6725},
			    {0.0226, 0.0080, 0.0130, 0.0058, 1.0463, 1.0896, 1.5551, 1.0986},
			    {0.0085, 0.0050, 0.0062, 0.0036, 0.6272, 1.1132, 1.2206, 0.9006},
			    {0.0140, 0.0056, 0.0165, 0.0017, 3.0113, 0.6747, 0.7976, 0.9445},
			    {0.0124, 0.0049, 0.0126, 0.0008, 1.7943, 0.4494, 0.8758, 0.1415},
			    {0.0239, 0.0105, 0.0035, 0.0014, 2.3162, 0.3501, 0.3520, 0.9412},
			    {0.0082, 0.0165, 0.0023, 0.0029, 0.1396, 1.9216, 1.3603, 0.7348},
			    {0.0201, 0.0114, 0.0001, 0.0018, 0.3037, 0.0699, 1.4453, 1.2160},
			    {0.0076, 0.0096, 0.0097, 0.0024, 0.0014, 0.4569, 1.1910, 1.1057},
			    {0.0118, 0.0109, 0.0046, 0.0020, 0.9412, 0.7662, 0.0690, 1.2182},
			    {0.0070, 0.0141, 0.0025, 0.0021, 1.1799, 0.2160, 0.0447, 0.8131},
			    {0.0075, 0.0136, 0.0065, 0.0002, 3.1297, 0.8345, 0.7822, 0.4496},
			    {0.0083, 0.0136, 0.0135, 0.0018, 1.3088, 1.7144, 0.8048, 0.7971},
			    {0.0063, 0.0029, 0.0106, 0.0020, 0.1659, 0.1256, 0.1526, 1.1587},
			    {0.0057, 0.0153, 0.0117, 0.0022, 2.7526, 1.2744, 1.5398, 0.0102},
			    {0.0125, 0.0117, 0.0084, 0.0006, 3.1014, 1.2814, 0.5067, 0.0922},
			    {0.0023, 0.0087, 0.0032, 0.0058, 0.1023, 1.2995, 0.7316, 0.6699},
			    {0.0072, 0.0085, 0.0055, 0.0013, 0.9979, 1.8470, 1.1499, 0.7783},
			    {0.0162, 0.0124, 0.0104, 0.0030, 1.0651, 1.7943, 1.1829, 0.6595},
			    {0.0016, 0.0219, 0.0049, 0.0018, 3.0881, 1.3501, 1.3021, 0.2507},
			    {0.0126, 0.0229, 0.0070, 0.0028, 2.8194, 1.4202, 0.2444, 0.6877},
			    {0.0026, 0.0146, 0.0107, 0.0002, 1.3380, 1.1348, 0.1109, 0.8790},
			    {0.0159, 0.0135, 0.0078, 0.0013, 2.3884, 0.2951, 0.7011, 0.7658},
			    {0.0130, 0.0120, 0.0040, 0.0026, 2.4187, 2.0936, 1.4007, 0.6005},
			    {0.0128, 0.0069, 0.0101, 0.0044, 0.2427, 1.4806, 0.0431, 1.1931},
			    {0.0112, 0.0112, 0.0052, 0.0025, 2.4205, 0.3927, 0.6790, 0.8408},
			    {0.0083, 0.0293, 0.0096, 0.0012, 1.3709, 0.3885, 1.1137, 0.9000},
			    {0.0104, 0.0055, 0.0104, 0.0027, 2.9709, 1.8674, 0.7436, 0.2211},
			    {0.0246, 0.0067, 0.0020, 0.0015, 0.3625, 1.2125, 1.2845, 0.2084},
			    {0.0157, 0.0089, 0.0103, 0.0032, 0.7275, 0.6670, 0.1099, 0.0024},
			    {0.0042, 0.0072, 0.0100, 0.0042, 0.2025, 0.5046, 0.9835, 0.7001},
			    {0.0112, 0.0171, 0.0030, 0.0021, 0.4649, 1.0000, 1.0488, 0.9570},
			    {0.0160, 0.0208, 0.0105, 0.0052, 0.6183, 1.4269, 1.2784, 0.6166},
			    {0.0128, 0.0144, 0.0041, 0.0007, 0.3341, 1.4410, 1.5192, 0.0196},
			    {0.0074, 0.0058, 0.0054, 0.0007, 0.6399, 0.5058, 0.4173, 0.4141},
			    {0.0151, 0.0118, 0.0122, 0.0027, 1.8413, 1.1375, 0.2133, 0.9026},
			    {0.0095, 0.0139, 0.0091, 0.0018, 0.7901, 1.5670, 0.6291, 0.7537},
			    {0.0051, 0.0198, 0.0104, 0.0019, 0.4317, 1.0279, 0.6844, 0.9169},
			    {0.0116, 0.0064, 0.0075, 0.0035, 2.9749, 0.4117, 1.4933, 0.6747},
			    {0.0194, 0.0220, 0.0017, 0.0061, 1.2260, 0.4513, 1.3605, 0.1055},
			    {0.0034, 0.0116, 0.0053, 0.0044, 1.0303, 0.4084, 0.8530, 0.1180},
			    {0.0058, 0.0114, 0.0045, 0.0023, 2.6296, 1.1444, 0.1301, 0.7712},
			    {0.0084, 0.0024, 0.0155, 0.0059, 0.3500, 0.8980, 1.3498, 0.1263},
			    {0.0125, 0.0207, 0.0079, 0.0024, 2.4658, 1.5593, 0.8358, 1.1297},
			    {0.0067, 0.0057, 0.0075, 0.0012, 2.4197, 1.2181, 0.9845, 0.4672},
			    {0.0050, 0.0187, 0.0092, 0.0026, 1.8256, 1.6191, 1.3539, 0.7405},
			    {0.0091, 0.0060, 0.0048, 0.0010, 1.8467, 1.5122, 1.3381, 0.2104},
			    {0.0121, 0.0229, 0.0009, 0.0005, 0.6522, 1.0287, 0.3377, 1.1998},
			    {0.0063, 0.0114, 0.0037, 0.0023, 2.0782, 2.0215, 0.8559, 1.1946},
			    {0.0139, 0.0134, 0.0030, 0.0016, 0.5237, 0.3628, 1.4648, 0.7004},
			    {0.0155, 0.0100, 0.0082, 0.0024, 1.6575, 0.7365, 1.2304, 0.7334},
			    {0.0117, 0.0169, 0.0040, 0.0019, 0.3560, 0.6076, 0.8880, 0.6829},
			    {0.0121, 0.0172, 0.0058, 0.0038, 2.9573, 1.5347, 0.2898, 0.5248},
			    {0.0133, 0.0067, 0.0017, 0.0048, 0.9285, 0.1553, 1.1599, 0.1877},
			    {0.0103, 0.0122, 0.0028, 0.0028, 0.4975, 0.0563, 0.0072, 0.4217},
			    {0.0180, 0.0236, 0.0049, 0.0048, 2.8175, 1.1828, 1.3087, 0.2473},
			    {0.0138, 0.0145, 0.0042, 0.0018, 1.1974, 0.4486, 1.2734, 0.6685},
			    {0.0134, 0.0026, 0.0018, 0.0039, 2.0075, 2.0158, 0.9690, 0.6126},
			    {0.0346, 0.0185, 0.0110, 0.0039, 1.1108, 1.8899, 0.4688, 0.8630},
			    {0.0192, 0.0072, 0.0076, 0.0033, 1.8109, 0.6485, 0.0930, 0.0311},
			    {0.0127, 0.0026, 0.0090, 0.0013, 2.3166, 1.4705, 1.0740, 0.7471},
			    {0.0115, 0.0036, 0.0115, 0.0036, 1.1312, 1.2609, 0.5996, 0.9439},
			    {0.0085, 0.0106, 0.0066, 0.0025, 0.0152, 1.5191, 0.0413, 0.2578},
			    {0.0102, 0.0067, 0.0030, 0.0024, 2.1989, 1.8232, 0.2486, 1.2449},
			    {0.0138, 0.0103, 0.0085, 0.0029, 0.5054, 0.9836, 1.0786, 0.7059},
			    {0.0123, 0.0086, 0.0076, 0.0019, 1.8668, 1.6994, 0.9880, 1.1505},
			    {0.0018, 0.0209, 0.0062, 0.0020, 2.2329, 0.1766, 0.0475, 0.8768},
			    {0.0062, 0.0156, 0.0039, 0.0055, 0.1620, 0.7111, 0.0465, 0.8902},
			    {0.0129, 0.0118, 0.0124, 0.0004, 1.0527, 2.0856, 1.4895, 0.9802},
			    {0.0130, 0.0224, 0.0067, 0.0023, 1.3647, 1.9240, 0.9774, 0.2057},
			    {0.0095, 0.0120, 0.0131, 0.0033, 1.3549, 0.0578, 1.0622, 0.1492},
			    {0.0181, 0.0095, 0.0048, 0.0018, 0.3657, 1.2464, 0.3167, 0.4021},
			    {0.0101, 0.0144, 0.0037, 0.0026, 0.4662, 1.8786, 0.2706, 0.0059},
			    {0.0070, 0.0282, 0.0065, 0.0019, 0.7439, 0.7537, 0.3298, 0.9921},
			    {0.0072, 0.0062, 0.0052, 0.0013, 0.1273, 2.0858, 1.3163, 0.7964},
			    {0.0126, 0.0085, 0.0094, 0.0022, 1.2513, 1.4222, 0.2437, 0.0859},
			    {0.0110, 0.0094, 0.0039, 0.0021, 1.1705, 1.2577, 0.6227, 1.0916},
			    {0.0147, 0.0115, 0.0034, 0.0048, 2.0546, 1.3980, 0.6977, 0.2892},
			    {0.0071, 0.0057, 0.0022, 0.0022, 0.4330, 0.3782, 0.3191, 0.4790},
			    {0.0028, 0.0083, 0.0053, 0.0010, 2.2178, 0.7593, 0.2022, 0.8853},
			    {0.0163, 0.0146, 0.0025, 0.0067, 0.6594, 1.8242, 0.5305, 0.3463},
			    {0.0131, 0.0159, 0.0064, 0.0011, 0.2721, 0.4152, 1.2204, 0.2509},
			    {0.0017, 0.0102, 0.0114, 0.0033, 1.7444, 0.7768, 0.5442, 0.2794},
			    {0.0064, 0.0136, 0.0053, 0.0016, 0.0440, 0.7431, 1.2696, 1.2395},
			    {0.0052, 0.0019, 0.0068, 0.0025, 1.7642, 0.1807, 0.5063, 0.0677},
			    {0.0053, 0.0123, 0.0099, 0.0016, 2.4677, 0.6623, 0.0884, 0.2277},
			    {0.0021, 0.0086, 0.0155, 0.0024, 0.4650, 1.9816, 1.3812, 1.0381},
			    {0.0207, 0.0054, 0.0099, 0.0018, 1.7223, 1.0553, 1.3898, 0.3571},
			    {0.0150, 0.0114, 0.0081, 0.0007, 0.2526, 0.1555, 0.1727, 0.8274},
			    {0.0153, 0.0051, 0.0085, 0.0030, 2.8053, 0.6042, 1.2212, 1.2452},
			    {0.0114, 0.0114, 0.0056, 0.0054, 2.8354, 2.0479, 1.5487, 1.1835},
			    {0.0089, 0.0099, 0.0060, 0.0028, 1.5864, 1.1997, 1.1121, 0.4190},
			    {0.0123, 0.0216, 0.0035, 0.0026, 0.9786, 0.7528, 0.7907, 0.5562},
			    {0.0032, 0.0109, 0.0024, 0.0028, 2.2915, 0.9272, 0.6995, 1.1891},
			    {0.0091, 0.0270, 0.0068, 0.0024, 2.8772, 0.8494, 1.1551, 0.3609},
			    {0.0129, 0.0144, 0.0049, 0.0031, 0.1821, 0.0452, 1.5706, 1.2319},
			    {0.0127, 0.0097, 0.0046, 0.0043, 0.2283, 1.6986, 0.8899, 1.1361},
			    {0.0021, 0.0287, 0.0132, 0.0002, 2.7721, 0.3207, 0.4076, 0.2857},
			    {0.0115, 0.0110, 0.0034, 0.0010, 0.5521, 1.4541, 1.0152, 0.6732},
			    {0.0282, 0.0102, 0.0123, 0.0023, 0.6063, 1.4107, 1.4128, 0.8899},
			    {0.0208, 0.0216, 0.0059, 0.0015, 2.6100, 0.4442, 0.3081, 0.8974},
			    {0.0084, 0.0071, 0.0019, 0.0037, 2.6944, 1.1982, 0.9367, 0.9872},
			    {0.0201, 0.0043, 0.0089, 0.0039, 2.1109, 0.3362, 1.0562, 1.0923},
			    {0.0119, 0.0056, 0.0107, 0.0060, 3.1011, 1.8405, 0.1200, 0.5129},
			    {0.0096, 0.0089, 0.0069, 0.0046, 3.0687, 1.2951, 0.4996, 0.9441},
			    {0.0266, 0.0038, 0.0123, 0.0024, 2.1246, 0.2021, 1.2239, 0.6545},
			    {0.0063, 0.0122, 0.0156, 0.0022, 2.2736, 0.4226, 1.0690, 0.2122},
			    {0.0271, 0.0119, 0.0026, 0.0003, 1.6497, 1.2688, 1.1221, 0.9900},
			    {0.0100, 0.0132, 0.0046, 0.0004, 0.9965, 0.9618, 1.1305, 0.6669},
			    {0.0056, 0.0023, 0.0111, 0.0020, 1.9226, 1.5838, 0.9000, 0.1033},
			    {0.0037, 0.0070, 0.0058, 0.0023, 2.0981, 2.0770, 1.1765, 0.5544},
			    {0.0239, 0.0168, 0.0095, 0.0013, 1.6507, 0.8755, 0.6177, 0.1818},
			    {0.0044, 0.0168, 0.0081, 0.0060, 1.8151, 0.2016, 1.0153, 0.0109},
			    {0.0157, 0.0155, 0.0038, 0.0047, 0.4066, 1.9201, 1.0704, 0.3928},
			    {0.0096, 0.0121, 0.0100, 0.0022, 2.1289, 2.0646, 0.4981, 0.7227},
			    {0.0179, 0.0150, 0.0003, 0.0060, 1.9955, 1.4489, 1.2964, 0.4596},
			    {0.0099, 0.0091, 0.0013, 0.0050, 2.3420, 1.2983, 0.8489, 0.6094},
			    {0.0085, 0.0083, 0.0137, 0.0060, 2.7641, 0.5335, 0.7058, 0.8528},
			    {0.0072, 0.0128, 0.0070, 0.0017, 2.5533, 1.1833, 1.1426, 0.3015},
			    {0.0178, 0.0059, 0.0123, 0.0014, 2.3903, 0.8255, 1.4489, 1.0062},
			    {0.0177, 0.0097, 0.0117, 0.0023, 0.3326, 1.1289, 0.0181, 0.5852},
			    {0.0092, 0.0097, 0.0024, 0.0018, 1.4998, 1.6096, 0.9054, 0.1004},
			    {0.0217, 0.0111, 0.0016, 0.0028, 0.9820, 1.4327, 1.4867, 0.4863},
			    {0.0047, 0.0080, 0.0100, 0.0045, 2.8094, 0.2580, 1.2854, 0.2201},
			    {0.0033, 0.0169, 0.0125, 0.0023, 2.1266, 1.0499, 0.7047, 1.0886},
			    {0.0114, 0.0103, 0.0087, 0.0033, 1.3277, 1.6925, 1.1173, 0.9879},
			    {0.0088, 0.0067, 0.0067, 0.0023, 1.1789, 1.5199, 0.4483, 0.1783},
			    {0.0066, 0.0240, 0.0061, 0.0035, 1.0477, 0.2353, 0.1057, 1.0293},
			    {0.0209, 0.0076, 0.0039, 0.0071, 1.3516, 1.0116, 1.1254, 0.4666},
			    {0.0144, 0.0172, 0.0150, 0.0019, 1.4880, 0.6142, 0.1007, 1.1708},
			    {0.0120, 0.0177, 0.0064, 0.0020, 0.3051, 0.7682, 1.2798, 0.4748},
			    {0.0127, 0.0049, 0.0064, 0.0048, 2.9249, 0.1780, 0.6767, 1.2302},
			    {0.0078, 0.0048, 0.0054, 0.0026, 3.0607, 0.7733, 0.0750, 0.6874},
			    {0.0121, 0.0143, 0.0085, 0.0030, 0.4726, 2.0846, 0.4344, 0.0027},
			    {0.0137, 0.0083, 0.0079, 0.0026, 0.1603, 0.7542, 0.2484, 0.2484},
			    {0.0165, 0.0115, 0.0092, 0.0030, 2.0919, 1.8560, 0.8101, 0.3244},
			    {0.0023, 0.0108, 0.0079, 0.0049, 2.9305, 1.8473, 1.3830, 0.5775},
			    {0.0122, 0.0186, 0.0128, 0.0029, 2.1956, 1.3914, 0.5269, 1.0969},
			    {0.0033, 0.0097, 0.0074, 0.0017, 2.1390, 0.1310, 0.5053, 1.1234},
			    {0.0104, 0.0012, 0.0130, 0.0013, 2.3245, 0.6883, 0.9553, 0.5956},
			    {0.0083, 0.0048, 0.0048, 0.0018, 2.9457, 0.6589, 1.2839, 0.7465},
			    {0.0112, 0.0148, 0.0051, 0.0040, 2.1775, 1.7562, 1.3573, 0.5653},
			    {0.0120, 0.0096, 0.0163, 0.0024, 1.3642, 0.5761, 0.0222, 0.4604},
			    {0.0166, 0.0117, 0.0084, 0.0009, 2.0378, 0.0877, 0.3153, 0.1669},
			    {0.0184, 0.0181, 0.0105, 0.0061, 3.0140, 2.0353, 1.3437, 1.2154},
			    {0.0156, 0.0103, 0.0081, 0.0007, 1.7236, 1.7334, 1.1225, 0.8141},
			    {0.0123, 0.0111, 0.0086, 0.0034, 2.8389, 1.5805, 0.0658, 0.3745},
			    {0.0122, 0.0093, 0.0069, 0.0030, 1.3305, 0.6095, 1.3544, 1.0004},
			    {0.0164, 0.0089, 0.0013, 0.0031, 1.6222, 0.5727, 0.1916, 1.2456},
			    {0.0143, 0.0036, 0.0112, 0.0014, 2.9718, 1.6541, 0.3906, 0.6607},
			    {0.0181, 0.0081, 0.0015, 0.0035, 2.5648, 0.3615, 0.3317, 0.4304},
			    {0.0182, 0.0134, 0.0042, 0.0028, 1.2501, 2.0225, 0.1822, 0.2684},
			    {0.0196, 0.0046, 0.0091, 0.0029, 2.2566, 1.1413, 0.2947, 0.4836},
			    {0.0194, 0.0120, 0.0080, 0.0026, 0.6228, 0.4399, 0.1849, 1.2185},
			    {0.0171, 0.0050, 0.0052, 0.0040, 2.7070, 0.7457, 1.0050, 0.4375},
			    {0.0223, 0.0155, 0.0091, 0.0014, 2.2510, 0.0226, 0.5278, 0.9606},
			    {0.0189, 0.0077, 0.0197, 0.0033, 0.9264, 0.8394, 1.2462, 0.6076},
			    {0.0089, 0.0181, 0.0032, 0.0054, 0.8931, 1.4755, 0.1916, 0.8407},
			    {0.0155, 0.0077, 0.0061, 0.0012, 0.1528, 0.4265, 0.4886, 0.0510},
			    {0.0088, 0.0099, 0.0114, 0.0012, 1.9023, 0.1775, 0.0271, 1.1395},
			    {0.0109, 0.0106, 0.0028, 0.0027, 0.5548, 0.0210, 1.3142, 0.5683},
			    {0.0082, 0.0041, 0.0044, 0.0043, 2.7331, 1.9869, 1.1365, 0.6074},
			    {0.0075, 0.0139, 0.0085, 0.0015, 0.2715, 0.4933, 0.9131, 1.2118},
			    {0.0105, 0.0125, 0.0073, 0.0048, 1.4184, 1.6238, 0.9015, 0.9230},
			    {0.0277, 0.0081, 0.0058, 0.0033, 0.1277, 1.5198, 1.4977, 0.0958},
			    {0.0114, 0.0293, 0.0091, 0.0027, 1.8356, 1.9889, 0.8336, 0.0298},
			    {0.0147, 0.0161, 0.0042, 0.0074, 0.1122, 1.8607, 0.6541, 0.0888},
			    {0.0206, 0.0058, 0.0086, 0.0009, 2.8307, 2.0617, 0.3877, 0.2164},
			    {0.0165, 0.0173, 0.0105, 0.0026, 0.9997, 1.2328, 1.4245, 0.6943},
			    {0.0044, 0.0182, 0.0091, 0.0014, 3.0848, 1.3357, 0.5952, 0.8780},
			    {0.0194, 0.0087, 0.0056, 0.0020, 3.1364, 0.1080, 1.1027, 0.1262},
			    {0.0123, 0.0043, 0.0095, 0.0033, 0.5469, 1.9568, 0.8445, 0.4266},
			    {0.0150, 0.0127, 0.0049, 0.0043, 1.6956, 1.3601, 0.8500, 1.1742},
			    {0.0115, 0.0062, 0.0047, 0.0006, 2.1551, 0.6965, 0.1356, 0.7501},
			    {0.0090, 0.0063, 0.0063, 0.0028, 1.7330, 0.2039, 0.8943, 1.2285},
			    {0.0140, 0.0106, 0.0039, 0.0009, 1.1902, 0.8273, 1.4574, 0.3725},
			    {0.0198, 0.0057, 0.0025, 0.0016, 0.0369, 0.5238, 0.9717, 0.2579},
			    {0.0086, 0.0060, 0.0057, 0.0003, 0.1255, 1.0270, 1.2220, 0.5405},
			    {0.0231, 0.0214, 0.0044, 0.0017, 0.7016, 1.1596, 1.0214, 0.8180},
			    {0.0152, 0.0100, 0.0117, 0.0005, 0.2866, 0.2131, 1.2717, 0.5844},
			    {0.0100, 0.0083, 0.0037, 0.0015, 3.1330, 1.8672, 0.8620, 0.1803},
			    {0.0138, 0.0096, 0.0064, 0.0034, 2.1978, 0.0754, 1.2565, 0.2397},
			    {0.0199, 0.0184, 0.0131, 0.0027, 1.3119, 2.0688, 1.0833, 0.9499},
			    {0.0041, 0.0233, 0.0067, 0.0054, 1.9096, 2.0399, 1.0795, 1.0596},
			    {0.0247, 0.0088, 0.0063, 0.0010, 2.9687, 1.0976, 0.2781, 0.4587},
			    {0.0108, 0.0041, 0.0111, 0.0028, 0.4263, 1.9856, 1.0337, 1.0204},
			    {0.0115, 0.0118, 0.0055, 0.0038, 2.4529, 0.9822, 0.5233, 1.1745},
			    {0.0107, 0.0163, 0.0067, 0.0024, 1.3844, 1.9843, 0.2476, 0.0320},
			    {0.0051, 0.0059, 0.0061, 0.0005, 1.8910, 1.3375, 0.7820, 0.0796},
			    {0.0024, 0.0083, 0.0101, 0.0022, 1.5223, 1.0696, 0.9698, 0.7865},
			    {0.0216, 0.0123, 0.0095, 0.0015, 2.6832, 0.6256, 1.1340, 0.9306},
			    {0.0065, 0.0182, 0.0075, 0.0017, 0.9990, 1.2047, 0.7111, 1.1854},
			    {0.0170, 0.0068, 0.0106, 0.0010, 0.4692, 0.4148, 1.2138, 0.7078},
			    {0.0043, 0.0183, 0.0089, 0.0007, 0.8870, 0.0637, 0.1933, 0.5632},
			    {0.0123, 0.0020, 0.0018, 0.0016, 2.0266, 0.4336, 0.8154, 1.2074},
			    {0.0138, 0.0064, 0.0043, 0.0024, 2.4468, 0.7047, 0.9627, 0.2181},
			    {0.0178, 0.0016, 0.0098, 0.0026, 1.9066, 1.7179, 0.8600, 0.3462},
			    {0.0146, 0.0043, 0.0039, 0.0019, 0.7890, 0.7445, 1.2881, 0.0924},
			    {0.0102, 0.0103, 0.0132, 0.0040, 3.0502, 1.7786, 1.4099, 0.0962},
			    {0.0114, 0.0123, 0.0019, 0.0019, 0.3048, 0.7776, 0.1184, 1.2314},
			    {0.0135, 0.0062, 0.0101, 0.0008, 0.5751, 0.5802, 1.4956, 1.2402},
			    {0.0104, 0.0134, 0.0129, 0.0010, 0.2119, 0.0125, 1.1365, 0.8467},
			    {0.0134, 0.0050, 0.0081, 0.0016, 1.6876, 1.8856, 0.4082, 0.1917},
			    {0.0100, 0.0196, 0.0054, 0.0030, 0.1195, 0.3899, 0.1431, 0.4342},
			    {0.0061, 0.0077, 0.0007, 0.0039, 0.6148, 1.0659, 0.8027, 1.2052},
			    {0.0048, 0.0191, 0.0069, 0.0033, 0.8793, 0.4803, 0.3874, 0.7770},
			    {0.0057, 0.0136, 0.0015, 0.0025, 0.5037, 0.1824, 0.6932, 0.8595},
			    {0.0091, 0.0177, 0.0092, 0.0042, 2.4779, 1.3925, 0.1751, 0.1453},
			    {0.0065, 0.0299, 0.0017, 0.0008, 1.4181, 0.2696, 0.0579, 0.7625},
			    {0.0376, 0.0202, 0.0033, 0.0016, 0.7216, 1.7748, 0.9500, 0.5592},
			    {0.0097, 0.0109, 0.0113, 0.0033, 2.5866, 1.9141, 0.9019, 0.0441},
			    {0.0128, 0.0008, 0.0060, 0.0046, 0.4965, 0.6755, 0.2703, 0.1433},
			    {0.0056, 0.0097, 0.0155, 0.0017, 1.0428, 1.1469, 0.8274, 1.1227},
			    {0.0029, 0.0133, 0.0109, 0.0049, 0.3445, 1.9402, 0.0050, 0.4723},
			    {0.0098, 0.0128, 0.0057, 0.0024, 1.8701, 0.9815, 1.3405, 1.0143},
			    {0.0111, 0.0250, 0.0034, 0.0010, 2.0078, 0.2191, 0.1272, 0.8923},
			    {0.0197, 0.0067, 0.0107, 0.0032, 0.1235, 0.5695, 0.7542, 1.1683},
			    {0.0181, 0.0207, 0.0094, 0.0072, 2.6084, 1.3190, 0.7297, 0.6604},
			    {0.0129, 0.0115, 0.0079, 0.0069, 3.0786, 1.1767, 1.5533, 0.7349},
			    {0.0122, 0.0238, 0.0129, 0.0022, 2.6770, 0.8137, 1.2751, 1.0020},
			    {0.0070, 0.0142, 0.0113, 0.0011, 2.6926, 1.5312, 0.9978, 0.8892},
			    {0.0050, 0.0141, 0.0180, 0.0033, 0.1019, 1.4817, 0.4370, 1.0886},
			    {0.0101, 0.0145, 0.0032, 0.0017, 1.8140, 1.2059, 0.2978, 1.1480},
			    {0.0175, 0.0151, 0.0098, 0.0019, 0.3334, 0.0315, 1.0528, 0.7383},
			    {0.0084, 0.0301, 0.0058, 0.0013, 0.2123, 1.2699, 0.5463, 0.5370},
			    {0.0119, 0.0071, 0.0121, 0.0016, 1.3785, 0.5978, 1.1371, 0.4401},
			    {0.0140, 0.0171, 0.0112, 0.0041, 2.2485, 1.4067, 0.2384, 0.4286},
			    {0.0195, 0.0166, 0.0110, 0.0038, 2.2849, 0.4291, 0.3410, 0.1126},
			    {0.0198, 0.0193, 0.0033, 0.0058, 2.2928, 0.0158, 0.0165, 0.7798},
			    {0.0050, 0.0189, 0.0036, 0.0021, 0.0582, 0.8407, 1.1294, 0.9280},
			    {0.0070, 0.0092, 0.0036, 0.0013, 1.1529, 1.9635, 1.0651, 0.0392},
			    {0.0027, 0.0058, 0.0051, 0.0048, 0.7714, 0.4794, 0.6610, 0.9076},
			    {0.0262, 0.0160, 0.0036, 0.0021, 2.4938, 1.8068, 0.5645, 0.0898},
			    {0.0126, 0.0157, 0.0105, 0.0020, 1.8734, 0.0085, 1.2675, 0.1232},
			    {0.0122, 0.0124, 0.0018, 0.0032, 1.2492, 1.3196, 0.6273, 0.3287},
			    {0.0131, 0.0114, 0.0049, 0.0014, 2.1214, 1.4683, 0.2798, 0.1446},
			    {0.0098, 0.0173, 0.0043, 0.0010, 0.8372, 0.9289, 0.3664, 1.1893},
			    {0.0056, 0.0096, 0.0006, 0.0015, 1.1941, 0.2447, 0.6154, 0.1081},
			    {0.0268, 0.0128, 0.0087, 0.0026, 1.3308, 0.6560, 0.1689, 0.7455},
			    {0.0195, 0.0086, 0.0133, 0.0039, 2.1138, 1.5816, 0.4573, 0.9342},
			    {0.0046, 0.0104, 0.0042, 0.0057, 0.3878, 0.5130, 1.4724, 0.0190},
			    {0.0033, 0.0215, 0.0134, 0.0027, 2.2140, 0.8516, 1.0710, 0.9896},
			    {0.0080, 0.0049, 0.0063, 0.0008, 2.6538, 1.8995, 1.2771, 0.6742},
			    {0.0111, 0.0076, 0.0082, 0.0039, 1.1687, 1.7572, 1.1559, 0.9726},
			    {0.0055, 0.0043, 0.0068, 0.0041, 2.1161, 1.4606, 1.1673, 1.2388},
			    {0.0151, 0.0276, 0.0042, 0.0026, 0.3085, 0.6697, 0.6069, 1.1600},
			    {0.0328, 0.0099, 0.0118, 0.0011, 1.9957, 0.2306, 0.9909, 0.3344},
			    {0.0191, 0.0089, 0.0038, 0.0043, 2.3634, 1.9017, 0.6815, 0.0040},
			    {0.0106, 0.0200, 0.0047, 0.0017, 2.5114, 0.5668, 0.6796, 1.1615},
			    {0.0051, 0.0082, 0.0032, 0.0057, 1.5984, 0.9255, 0.8304, 1.2412},
			    {0.0132, 0.0195, 0.0104, 0.0039, 0.1513, 1.6562, 0.4538, 0.9850},
			    {0.0126, 0.0129, 0.0054, 0.0020, 1.9958, 1.2345, 0.0817, 0.5748},
			    {0.0180, 0.0114, 0.0025, 0.0019, 0.3068, 1.5582, 0.7094, 0.7932},
			    {0.0014, 0.0179, 0.0113, 0.0016, 0.6548, 0.0357, 0.3491, 0.3729},
			    {0.0133, 0.0081, 0.0108, 0.0010, 1.2681, 0.0784, 0.8892, 1.1707},
			    {0.0094, 0.0149, 0.0071, 0.0025, 1.9509, 1.4684, 0.2666, 0.2388},
			    {0.0019, 0.0126, 0.0041, 0.0055, 2.5134, 1.7356, 0.9304, 0.1254},
			    {0.0212, 0.0222, 0.0111, 0.0037, 1.8128, 0.7586, 1.5264, 0.5764},
			    {0.0018, 0.0232, 0.0032, 0.0021, 1.8452, 1.7118, 0.4042, 0.3207},
			    {0.0116, 0.0129, 0.0032, 0.0031, 1.7869, 1.8972, 1.0492, 0.4099},
			    {0.0204, 0.0061, 0.0045, 0.0015, 1.5150, 0.1407, 0.8372, 0.2237},
			    {0.0110, 0.0303, 0.0076, 0.0019, 0.4509, 1.1437, 0.7510, 1.2134},
			    {0.0144, 0.0179, 0.0081, 0.0015, 0.1253, 1.4838, 0.6993, 1.2482},
			    {0.0032, 0.0323, 0.0085, 0.0012, 1.2771, 2.0519, 0.7691, 0.5044},
			    {0.0172, 0.0266, 0.0050, 0.0012, 0.0311, 2.0403, 1.4857, 0.0537},
			    {0.0233, 0.0184, 0.0080, 0.0044, 2.5444, 1.9774, 0.2290, 1.2328},
			    {0.0149, 0.0082, 0.0046, 0.0018, 1.1703, 2.0565, 1.2349, 0.3740},
			    {0.0241, 0.0233, 0.0117, 0.0035, 2.0348, 1.0137, 0.7176, 0.7898},
			    {0.0120, 0.0068, 0.0069, 0.0023, 2.2351, 1.9502, 0.7107, 0.8255},
			    {0.0180, 0.0071, 0.0061, 0.0029, 3.0668, 0.5662, 0.4843, 0.3449},
			    {0.0046, 0.0176, 0.0146, 0.0005, 2.5892, 1.6997, 0.5470, 0.5327},
			    {0.0217, 0.0279, 0.0095, 0.0017, 3.1101, 0.1986, 0.7152, 0.8793},
			    {0.0076, 0.0038, 0.0048, 0.0013, 2.5997, 1.8645, 0.6158, 1.0507},
			    {0.0114, 0.0122, 0.0160, 0.0053, 3.0017, 0.2995, 0.6446, 0.0752},
			    {0.0117, 0.0118, 0.0054, 0.0031, 0.5203, 0.2113, 0.4019, 0.2653},
			    {0.0170, 0.0176, 0.0109, 0.0014, 2.9938, 0.4597, 0.9603, 0.5441},
			    {0.0032, 0.0166, 0.0104, 0.0018, 2.2736, 0.1183, 1.5353, 0.9630},
			    {0.0140, 0.0026, 0.0123, 0.0038, 0.9841, 0.5333, 0.1106, 0.9818},
			    {0.0186, 0.0095, 0.0154, 0.0044, 1.0710, 1.1085, 0.5786, 0.4476},
			    {0.0036, 0.0074, 0.0065, 0.0054, 0.1991, 0.4954, 1.4186, 1.2296},
			    {0.0171, 0.0106, 0.0057, 0.0015, 1.3980, 1.9400, 0.5324, 0.2100},
			    {0.0017, 0.0108, 0.0061, 0.0068, 0.2364, 1.9401, 0.6312, 0.6766},
			    {0.0072, 0.0100, 0.0064, 0.0027, 2.6298, 1.9040, 1.5324, 0.4255},
			    {0.0049, 0.0078, 0.0108, 0.0019, 2.8918, 0.5057, 1.0193, 0.0276},
			    {0.0134, 0.0150, 0.0094, 0.0052, 2.7824, 1.4790, 1.0799, 0.1170},
			    {0.0059, 0.0174, 0.0037, 0.0020, 1.1473, 1.5172, 1.4863, 0.7157},
			    {0.0104, 0.0147, 0.0050, 0.0034, 2.8548, 1.6373, 1.1067, 0.1058},
			    {0.0078, 0.0141, 0.0156, 0.0013, 0.3839, 2.0424, 1.2524, 1.2236},
			    {0.0168, 0.0060, 0.0140, 0.0059, 0.3661, 0.1142, 1.0174, 0.6876},
			    {0.0075, 0.0126, 0.0113, 0.0044, 1.0701, 1.7237, 0.3929, 0.5876},
			    {0.0105, 0.0052, 0.0051, 0.0033, 2.8024, 1.7995, 0.9513, 0.0596},
			    {0.0313, 0.0098, 0.0116, 0.0041, 0.0316, 2.0876, 0.9225, 0.2930},
			    {0.0095, 0.0102, 0.0045, 0.0034, 1.0195, 0.1796, 0.3268, 0.1556},
			    {0.0100, 0.0119, 0.0015, 0.0008, 1.3102, 1.9364, 1.4832, 0.9270},
			    {0.0073, 0.0040, 0.0015, 0.0031, 1.0762, 0.4398, 0.4016, 0.7630},
			    {0.0059, 0.0223, 0.0041, 0.0047, 1.2775, 0.9970, 1.1847, 0.4077},
			    {0.0288, 0.0175, 0.0099, 0.0036, 1.2865, 0.8372, 0.3624, 0.3249},
			    {0.0038, 0.0037, 0.0070, 0.0042, 0.2823, 1.3240, 0.6374, 1.0790},
			    {0.0046, 0.0203, 0.0088, 0.0013, 0.0241, 0.5188, 1.1429, 0.6752},
			    {0.0085, 0.0225, 0.0054, 0.0019, 1.2467, 0.1586, 0.5029, 0.9428},
			    {0.0200, 0.0110, 0.0077, 0.0034, 1.7224, 1.0048, 1.4214, 1.0760},
			    {0.0107, 0.0110, 0.0062, 0.0036, 1.3003, 1.8216, 1.4260, 0.8645},
			    {0.0098, 0.0108, 0.0033, 0.0040, 2.3488, 1.6367, 1.2396, 1.1924},
			    {0.0100, 0.0072, 0.0106, 0.0017, 0.9517, 1.7579, 1.1148, 0.1767},
			    {0.0070, 0.0092, 0.0032, 0.0028, 0.3840, 0.0812, 0.2991, 0.2998},
			    {0.0174, 0.0177, 0.0100, 0.0008, 0.9373, 0.4938, 0.3044, 0.3579},
			    {0.0092, 0.0134, 0.0080, 0.0009, 1.2092, 0.4288, 0.5479, 0.9964},
			    {0.0037, 0.0164, 0.0051, 0.0016, 2.1923, 0.8158, 0.0667, 1.0408},
			    {0.0050, 0.0137, 0.0139, 0.0019, 0.0359, 0.0012, 1.0014, 0.5098},
			    {0.0055, 0.0043, 0.0076, 0.0029, 2.3821, 1.3232, 0.4201, 1.2520},
			    {0.0077, 0.0066, 0.0113, 0.0029, 2.8059, 2.0480, 0.2036, 0.0844},
			    {0.0209, 0.0207, 0.0049, 0.0017, 0.0976, 1.6299, 1.1144, 1.0213},
			    {0.0100, 0.0192, 0.0047, 0.0010, 1.0883, 0.9764, 0.2155, 1.1940},
			    {0.0038, 0.0075, 0.0128, 0.0014, 0.6894, 0.8497, 1.5570, 0.2692},
			    {0.0110, 0.0042, 0.0027, 0.0056, 1.7175, 1.7685, 0.6648, 1.1603},
			    {0.0133, 0.0046, 0.0028, 0.0042, 1.2470, 1.3693, 0.4291, 0.1843},
			    {0.0084, 0.0207, 0.0192, 0.0010, 1.7049, 0.6654, 0.1892, 0.9108},
			    {0.0114, 0.0192, 0.0111, 0.0035, 0.1111, 0.2245, 0.6405, 0.4412},
			    {0.0050, 0.0130, 0.0058, 0.0041, 2.9179, 1.3901, 0.5407, 0.4174},
			    {0.0180, 0.0028, 0.0052, 0.0020, 0.7291, 1.9041, 1.0980, 0.4874},
			    {0.0177, 0.0050, 0.0144, 0.0021, 1.2765, 1.9586, 1.5685, 0.4539},
			    {0.0135, 0.0226, 0.0091, 0.0029, 2.7607, 1.1069, 1.4427, 1.2103},
			    {0.0076, 0.0059, 0.0095, 0.0032, 2.3900, 0.3554, 0.8471, 1.1871},
			    {0.0065, 0.0125, 0.0113, 0.0019, 0.4172, 0.9931, 0.7135, 0.3500},
			    {0.0100, 0.0019, 0.0052, 0.0018, 0.7083, 1.9296, 1.2251, 0.3835},
			    {0.0175, 0.0027, 0.0041, 0.0032, 0.3757, 1.3088, 0.7385, 0.0136},
			    {0.0131, 0.0012, 0.0037, 0.0018, 2.7581, 0.4501, 0.8392, 0.4057},
			    {0.0162, 0.0055, 0.0063, 0.0008, 0.9309, 1.9077, 0.6246, 0.6869},
			    {0.0092, 0.0066, 0.0061, 0.0027, 0.4270, 0.5592, 0.5984, 0.3664},
			    {0.0093, 0.0155, 0.0096, 0.0012, 0.0888, 0.0138, 0.9761, 0.9387},
			    {0.0232, 0.0030, 0.0012, 0.0018, 1.0103, 1.3105, 1.1543, 0.4580},
			    {0.0091, 0.0138, 0.0117, 0.0010, 2.0275, 0.1862, 1.0082, 0.7152},
			    {0.0046, 0.0146, 0.0058, 0.0006, 2.7164, 0.4973, 0.3711, 0.8929},
			    {0.0124, 0.0039, 0.0064, 0.0029, 0.0734, 0.5001, 0.1700, 0.7015},
			    {0.0024, 0.0064, 0.0024, 0.0009, 0.9627, 1.4399, 0.3649, 0.2004},
			    {0.0037, 0.0133, 0.0097, 0.0031, 0.2014, 1.1774, 1.2642, 0.6371},
			    {0.0344, 0.0101, 0.0085, 0.0027, 1.4060, 0.1155, 0.2356, 0.0542},
			    {0.0073, 0.0013, 0.0048, 0.0004, 2.8539, 1.5881, 0.3362, 0.6969},
			    {0.0128, 0.0161, 0.0071, 0.0053, 2.3949, 0.1941, 0.1206, 0.1499},
			    {0.0135, 0.0178, 0.0122, 0.0061, 2.3731, 0.0056, 1.3851, 0.3842},
			    {0.0097, 0.0247, 0.0013, 0.0018, 0.8335, 1.3187, 1.3054, 0.6448},
			    {0.0045, 0.0117, 0.0104, 0.0024, 2.5110, 1.1067, 0.9396, 0.5052},
			    {0.0072, 0.0117, 0.0074, 0.0039, 0.7208, 0.6408, 1.1329, 0.5908},
			    {0.0093, 0.0234, 0.0133, 0.0038, 2.3006, 0.2541, 1.1034, 0.2324},
			    {0.0109, 0.0205, 0.0155, 0.0067, 1.5636, 2.0396, 0.8667, 0.0112},
			    {0.0090, 0.0172, 0.0041, 0.0018, 0.2388, 1.4317, 0.2854, 1.1759},
			    {0.0192, 0.0101, 0.0088, 0.0010, 2.1545, 0.7999, 1.5177, 0.3061},
			    {0.0105, 0.0095, 0.0007, 0.0018, 1.3741, 1.2807, 0.6832, 0.7749},
			    {0.0064, 0.0121, 0.0078, 0.0034, 0.5180, 1.6747, 1.3113, 0.2402},
			    {0.0039, 0.0031, 0.0129, 0.0026, 0.7717, 1.4115, 0.7525, 0.0632},
			    {0.0099, 0.0237, 0.0029, 0.0057, 1.5999, 1.6389, 0.6722, 0.7793},
			    {0.0031, 0.0119, 0.0068, 0.0051, 1.4428, 1.6590, 0.6001, 0.3696},
			    {0.0016, 0.0177, 0.0055, 0.0007, 1.3757, 0.5399, 1.4051, 0.3330},
			    {0.0132, 0.0062, 0.0092, 0.0042, 1.8616, 1.7299, 1.2276, 0.3553},
			    {0.0102, 0.0149, 0.0062, 0.0026, 3.0159, 0.9631, 0.4066, 0.7634},
			    {0.0202, 0.0032, 0.0028, 0.0011, 3.0440, 1.6730, 0.0659, 0.4190},
			    {0.0040, 0.0199, 0.0072, 0.0049, 2.8376, 1.1899, 0.6880, 0.3716},
			    {0.0072, 0.0105, 0.0133, 0.0076, 2.2572, 0.5076, 0.6271, 0.1691},
			    {0.0101, 0.0109, 0.0068, 0.0055, 0.2971, 0.5915, 1.1698, 0.0515},
			    {0.0170, 0.0167, 0.0075, 0.0015, 0.3839, 0.9573, 0.6347, 0.6510},
			    {0.0117, 0.0061, 0.0030, 0.0009, 2.6093, 0.4259, 0.4273, 0.8861},
			    {0.0042, 0.0022, 0.0094, 0.0020, 1.3445, 0.5421, 0.6545, 0.9748},
			    {0.0178, 0.0110, 0.0057, 0.0070, 0.9030, 1.1982, 0.9715, 1.1535},
			    {0.0074, 0.0145, 0.0118, 0.0015, 3.0026, 1.8525, 0.9991, 0.5912},
			    {0.0235, 0.0126, 0.0106, 0.0017, 1.5768, 0.9108, 0.5716, 0.2727},
			    {0.0167, 0.0148, 0.0116, 0.0017, 1.0353, 1.7619, 0.3951, 0.2947},
			    {0.0131, 0.0139, 0.0109, 0.0044, 2.2321, 0.1821, 1.5191, 0.3270},
			    {0.0144, 0.0137, 0.0050, 0.0013, 3.1220, 0.4640, 1.2267, 1.1019},
			    {0.0214, 0.0026, 0.0054, 0.0024, 2.3978, 1.7252, 0.1594, 0.6990},
			    {0.0141, 0.0074, 0.0069, 0.0013, 1.0683, 0.0628, 0.9449, 0.0837},
			    {0.0116, 0.0226, 0.0047, 0.0007, 2.1564, 1.8269, 0.8714, 1.2005},
			    {0.0151, 0.0241, 0.0097, 0.0034, 0.2904, 0.6657, 0.5823, 1.0473},
			    {0.0158, 0.0054, 0.0013, 0.0023, 2.0291, 0.1255, 0.2251, 0.1777},
			    {0.0095, 0.0194, 0.0045, 0.0010, 0.1391, 0.6952, 0.7655, 0.3528},
			    {0.0323, 0.0077, 0.0093, 0.0025, 2.5059, 0.8143, 0.9753, 0.4365},
			    {0.0088, 0.0196, 0.0029, 0.0056, 1.0531, 1.2262, 0.7902, 0.5696},
			    {0.0058, 0.0103, 0.0016, 0.0027, 1.3645, 1.0926, 0.6892, 0.1396},
			    {0.0048, 0.0235, 0.0062, 0.0052, 0.1360, 2.0105, 1.4386, 0.0478},
			    {0.0180, 0.0097, 0.0112, 0.0019, 2.5605, 1.0203, 0.3819, 1.1079},
			    {0.0189, 0.0054, 0.0045, 0.0030, 2.4685, 1.2378, 1.3858, 1.0875},
			    {0.0201, 0.0094, 0.0029, 0.0010, 0.9477, 1.4278, 0.8179, 1.0957},
			    {0.0079, 0.0127, 0.0083, 0.0043, 1.3163, 1.3207, 0.2192, 0.0845},
			    {0.0144, 0.0169, 0.0012, 0.0010, 2.0262, 1.9272, 1.4234, 1.2528},
			    {0.0083, 0.0047, 0.0096, 0.0035, 3.1353, 1.9535, 0.7535, 0.9366},
			    {0.0129, 0.0134, 0.0091, 0.0025, 0.5241, 0.5672, 1.3602, 0.9937},
			    {0.0129, 0.0165, 0.0059, 0.0034, 2.0923, 1.1969, 1.0181, 0.3233},
			    {0.0218, 0.0118, 0.0010, 0.0037, 1.5175, 0.3109, 0.3593, 0.8776},
			    {0.0020, 0.0075, 0.0032, 0.0020, 1.2901, 0.8482, 0.6895, 0.2847},
			    {0.0136, 0.0145, 0.0114, 0.0019, 3.0816, 1.1904, 1.2484, 0.1264},
			    {0.0126, 0.0114, 0.0034, 0.0060, 2.5748, 1.5052, 1.2225, 0.1405},
			    {0.0163, 0.0186, 0.0033, 0.0023, 1.9710, 1.2308, 1.1154, 1.0442},
			    {0.0088, 0.0158, 0.0052, 0.0011, 0.7200, 0.2429, 0.3846, 0.5091},
			    {0.0083, 0.0065, 0.0064, 0.0042, 1.7709, 1.4932, 1.4849, 0.9873},
			    {0.0158, 0.0158, 0.0097, 0.0006, 1.3281, 1.2261, 1.5095, 1.0002},
			    {0.0147, 0.0122, 0.0062, 0.0027, 2.6676, 2.0461, 0.4477, 1.1191},
			    {0.0119, 0.0136, 0.0100, 0.0039, 2.1697, 1.7528, 1.0626, 0.1126},
			    {0.0290, 0.0204, 0.0063, 0.0020, 1.2412, 1.8353, 0.1385, 1.2566}}; 

  float v2=fFitHighV2->Eval(3.);
  float v3=fFitHighV3->Eval(3.);
  float v4=fFitHighV4->Eval(3.);
  float v5=0;
  float Psi5=0;

  if(fSimulate==3){
    fSimFlowMark=int(0.5+500*gRandom->Rndm());
    int aa=fSimFlowMark;
    v2=3*FlowValues[aa][0];
    v3=3*FlowValues[aa][1];
    v4=3*FlowValues[aa][2];
    v5=3*FlowValues[aa][3];
    fSimPsi2=FlowValues[aa][4];
    fSimPsi3=FlowValues[aa][5];
    fSimPsi4=FlowValues[aa][6];
    Psi5=FlowValues[aa][7];
   
    if(fSimFlowMark>=500)fSimFlowMark=0;
  }
  if(fSimulate==1){
    fSimPsi2=2*TMath::Pi()*gRandom->Rndm();
    fSimPsi3=2*TMath::Pi()*gRandom->Rndm();
    fSimPsi4=fSimPsi2;
  }
  
  float randAng=2*TMath::Pi()*gRandom->Rndm();
  fSimPsi2+=randAng;
  fSimPsi3+=randAng;
  fSimPsi4+=randAng;
  Psi5+=randAng;
  
  Int_t lead=0;
  rGoodTracks[0]=0;
  TF1 *SimFlow=new TF1("SimFlow","1+2*[0]*cos(2*(x-[1]))+2*[2]*cos(3*(x-[3]))+2*[4]*cos(4*(x-[5]))+2*[6]*cos(5*(x-[6]))",-TMath::Pi(),TMath::Pi());
  
  SimFlow->SetParameters(v2,fSimPsi2,v3,fSimPsi3,v4,fSimPsi4,v5,Psi5);
  // SimFlow->SetParameters(0,0,0,0);
  TF1 *SimNear=new TF1("SimNear","exp(-0.5*x*x/[0]/[0])",-TMath::Pi(),TMath::Pi());
  SimNear->SetParameter(0,0.3);
  TF1 *SimAway=new TF1("SimAway","exp(-0.5*x*x/[0]/[0])",-TMath::Pi(),TMath::Pi());
  SimAway->SetParameter(0,0.5);//0.5 for DiJet 0.3 for cone

  //used fixed v2,v3,v4 for mixing
  v2=0.07,v3=0.05,v4=0.03;
  TF1 *SimAway2=new TF1("SimAway2","exp(-0.5*(x-[1])*(x-[1])/[0]/[0])+exp(-0.5*(x+[1])*(x+[1])/[0]/[0])",-TMath::Pi(),TMath::Pi());
  SimAway2->SetParameter(0,0.3);//0.5 for DiJet 0.3 for cone
  SimAway2->SetParameter(1,1.4);//Cone Angle

  TF1 *AwayProb=new TF1("AwayProb","[0]");
  // TF1 *AwayProb=new TF1("AwayProb","[0]*cos(2*x)*cos(2*x)");
  AwayProb->SetParameter(0,0.5);
  // AwayProb->SetParameter(0,1);

  // TF1 *SimAway=new TF1("SimAway","0.035+2*[0]*[0]*cos(2*x)+2*[1]*[1]*cos(3*x)+2*[2]*[2]*cos(4*x)",1.0,TMath::Pi()*2-1.0);
  // SimAway->SetParameters(v2,v3,v4);
  // TF1 *SimAway2=new TF1("SimAway2","0.035+2*[0]*[0]*cos(2*x)+2*[1]*[1]*cos(3*x)+2*[2]*[2]*cos(4*x)",1.05,TMath::Pi()*2-1.05);
  //SimAway->SetParameters(v2,v3,v4);
 
  
  Float_t TrigAngle;
  Float_t sPt,sPhi;
  Int_t InAccpt;
  Int_t AccptPercent=4;//1 over this is % in aceptance on away-side 
  Int_t AwaySidePM=0;

  Int_t AwaySide1=1;
  //Use SimAway1 or 2
  // if(gRandom->Rndm()<AwayProb->Eval(RPAngle))AwaySide1=1;
  //else AwaySide1=2;
  if(fSimNJet<1)fSimNJet=1;
  for(int i=0;i<gRandom->Poisson(fSimNBgPart);i++){
    sPt=1.5;
    rPt[rGoodTracks[0]]=sPt;
    rEta[rGoodTracks[0]]=0;
    sPhi=SimFlow->GetRandom();
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
    TrigAngle=SimFlow->GetRandom();
    if(gRandom->Rndm()<AwayProb->Eval(TrigAngle-fSimPsi2))AwaySide1=1;
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
      if(fSimAwayDeflected){
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
    /* removed this cross check b/c not coverity complient
    if((!fESD&&!fAOD)&&ievent==0){
      if(fDEBUG)Printf("Error: fESD not found");
      break;
    }
    */
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
    if(fSimulate==1||fSimulate==3)nTrack=int(10*(fSimNBgPart+10*fSimNJetPart));
    if(fSimulate==2){
       if(!fAODData)nTrack=fESD->GetNumberOfTracks();
      else nTrack=fAOD->GetNumberOfTracks();
    }
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
    if(ievent)leadPart=TrackCutsMC(fMC,ftPt,ftEta,ftPhi,ftCharge,ftEff,ftV2,ftV3,ftV4,ftPtAssoc3,ftNPtAssoc3,nGoodTracks);
    else if(fSimulate!=1||fSimulate!=3){
      if(!fAODData)leadPart=TrackCuts(fESD,ftPt,ftEta,ftPhi,ftCharge,ftEff,ftV2,ftV3,ftV4,ftPtAssoc3,ftNPtAssoc3,nGoodTracks);
      else leadPart=TrackCutsAOD(fAOD,ftPt,ftEta,ftPhi,ftCharge,ftEff,ftV2,ftV3,ftV4,ftPtAssoc3,ftNPtAssoc3,nGoodTracks);
      if(fSimulate==2)CalcFlow(ftPt,ftEta,ftPhi,nGoodTracks,leadPart);
    }
    
    for(int nSimEvents=0;nSimEvents<=(fSimulate*fSimNEvents);nSimEvents++){//only 1 loop if not simulation
      //Printf("nSimEvents %d",nSimEvents);
      if(fSimulate)leadPart=TrackCutsSim(ftPt,ftEta,ftPhi,ftCharge,ftEff,ftV2,ftV3,ftV4,ftPtAssoc3,ftNPtAssoc3,nGoodTracks);
      
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
