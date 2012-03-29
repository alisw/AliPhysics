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

//==============================================================================
// AliHMPIDPIDEfficiencyV0AnalysysTask - Class representing a analysis tool to evaluate the PID efficiency for the HMPID detector using the V0's 
// 
//
// HMPID PID Efficiency using the V0
//
// Francesco Barile
//
//==============================================================================
//
// By means of AliHMPIDAnalysisTask.C macro it is possible to use this class
// to perform the analysis on local data, on data on alien using local machine
// and on CAF.

#ifndef AliHMPIDPIDEfficiencyV0AnalysisTASK_CXX
#define AliHMPIDPIDEfficiencyV0AnalysisTASK_CXX


#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCEvent.h"

#include "AliESDtrackCuts.h"

#include "AliESDv0.h"
#include "AliESDv0Cuts.h"

#include "AliPID.h"
#include "AliLog.h"
#include "AliHMPIDPIDEfficiencyV0AnalysisTask.h"

ClassImp(AliHMPIDPIDEfficiencyV0AnalysisTask)

//__________________________________________________________________________
AliHMPIDPIDEfficiencyV0AnalysisTask::AliHMPIDPIDEfficiencyV0AnalysisTask() :
  fESD(0x0),fMC(0x0),fUseMC(kTRUE),
  fHmpHistList(0x0),
  fHmpNevents(0x0),
  fThetavsMom(0x0),
  fmass(0x0),
  fpdg(0x0),
  farme(0x0),
  farmeb(0x0),
  hangle(0x0),
  massap(0x0),
  massan(0x0),
  fmassaHMPID(0x0),
  hdiff(0x0),
  hdiffn(0x0),
  fTree(0x0)
{
  //
  //Default ctor
  //
  for (Int_t i=0; i<90; i++) fVar[i]=0;
}

//___________________________________________________________________________
AliHMPIDPIDEfficiencyV0AnalysisTask::AliHMPIDPIDEfficiencyV0AnalysisTask(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fESD(0x0), fMC(0x0), fUseMC(kTRUE),
  fHmpHistList(0x0),
  fHmpNevents(0x0),
  fThetavsMom(0x0),
  fmass(0x0),
  fpdg(0x0),
  farme(0x0),
  farmeb(0x0),
  hangle(0x0),
  massap(0x0),
  massan(0x0),
  fmassaHMPID(0x0),
  hdiff(0x0),
  hdiffn(0x0),  
  fTree(0x0)
  {
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  for (Int_t i=0; i<90; i++) fVar[i]=0;
  DefineOutput(1,TList::Class());
  DefineOutput(2,TTree::Class());
}

//___________________________________________________________________________
AliHMPIDPIDEfficiencyV0AnalysisTask& AliHMPIDPIDEfficiencyV0AnalysisTask::operator=(const AliHMPIDPIDEfficiencyV0AnalysisTask& c)
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c);
    fESD             = c.fESD;
    fMC              = c.fMC;
    fUseMC           = c.fUseMC;
    fHmpHistList     = c.fHmpHistList;
    fHmpNevents      = c.fHmpNevents;
    fThetavsMom      = c.fThetavsMom;
    fmass            = c.fmass;
    fpdg             = c.fpdg;
    farme            = c.farme;
    farmeb           = c.farmeb;
    hangle           = c.hangle;
    massap           = c.massap;
    massan           = c.massan;
    fmassaHMPID      = c.fmassaHMPID;
    hdiff            = c.hdiff;
    hdiffn           = c.hdiffn;
    fTree            = c.fTree;
    for(Int_t i=0; i<90; i++) fVar[i]=c.fVar[i];
 }
  return *this;
}

//___________________________________________________________________________
AliHMPIDPIDEfficiencyV0AnalysisTask::AliHMPIDPIDEfficiencyV0AnalysisTask(const AliHMPIDPIDEfficiencyV0AnalysisTask& c) :
  AliAnalysisTaskSE(c),
  fESD(c.fESD),fMC(c.fMC),fUseMC(c.fUseMC),
  fHmpHistList(c.fHmpHistList),
  fHmpNevents(c.fHmpNevents),
  fThetavsMom(c.fThetavsMom),
  fmass(c.fmass),
  fpdg(c.fpdg),
  farme(c.farme),
  farmeb(c.farmeb),
  hangle(c.hangle),
  massap(c.massap),
  massan(c.massan),
  fmassaHMPID(c.fmassaHMPID),
  hdiff(c.hdiff),
  hdiffn(c.hdiffn),
  fTree(c.fTree)  
{
  //
  // Copy Constructor
  //
  for (Int_t i=0; i<90; i++) fVar[i]=c.fVar[i];
}
 
//___________________________________________________________________________
AliHMPIDPIDEfficiencyV0AnalysisTask::~AliHMPIDPIDEfficiencyV0AnalysisTask() {
  //
  //destructor
  //
  Info("~AliHMPIDPIDEfficiencyV0AnalysisTask","Calling Destructor");
  if (fHmpHistList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fHmpHistList;
}

//___________________________________________________________________________
void AliHMPIDPIDEfficiencyV0AnalysisTask::ConnectInputData(Option_t *)
{
  // Connect ESD here

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    AliDebug(2,Form("ERROR: Could not get ESDInputHandler"));
  } else
    fESD = esdH->GetEvent();

  if (fUseMC){
    // Connect MC
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!mcH) {
      AliDebug(2,Form("ERROR: Could not get MCEventHandler"));
      fUseMC = kFALSE;
    } else
      fMC = mcH->MCEvent();
      if (!fMC) AliDebug(2,Form("ERROR: Could not get MCEvent"));
  }
}

//___________________________________________________________________________
void AliHMPIDPIDEfficiencyV0AnalysisTask::UserExec(Option_t *)
{
  
  fHmpNevents->Fill(0);

  // Variables:
  Double_t  lPrimaryVtxPosition[3];
  Double_t lV0Position[3];
  
  Double_t lDcaPosToPrimVertex = 0;
  Double_t lDcaNegToPrimVertex = 0;
  //Double_t lDcaV0Daughters     = 0;
  
  Double_t lMagneticField = 999;

  
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  AliESDtrackCuts *myTracksCuts = NULL;
  if (!vertex) return;
  if (!(vertex->GetStatus())) return;  
  vertex->GetXYZ(lPrimaryVtxPosition);
  
  lMagneticField = fESD->GetMagneticField();

  myTracksCuts = new AliESDtrackCuts();
  myTracksCuts->SetRequireTPCRefit(kTRUE);
    
  if(!vertex || !vertex->GetStatus() || vertex->GetNContributors()<1) {
    
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    if(!vertex) return;
    if(!vertex->GetStatus()) return;
    if(vertex->GetNContributors()<1) return; // no good vertex, skip event
    
  }
  
  Double_t vtxPos[3] = {999, 999, 999};
  if(vertex) vertex->GetXYZ(vtxPos);
  
  if(vertex) {
    lPrimaryVtxPosition[0] = vertex->GetX();
    lPrimaryVtxPosition[1] = vertex->GetY();
    lPrimaryVtxPosition[2] = vertex->GetZ();
             
             }
  
  fHmpNevents->Fill(1);
  
  
  
  // Daughters' momentum:
  Double_t  lMomPos[3] = {999,999,999};
  Double_t  lMomNeg[3] = {999,999,999};

  // Daughters    
  AliESDtrack *pTrack = NULL;
  AliESDtrack *nTrack = NULL;
  
  AliStack* pStack = 0;
  if (fUseMC){
    pStack = fMC->Stack();
  }

  Double_t r[3]; 
  Double_t rout[3];
  
  //
  // Main loop function, executed on Event basis
  //
  //...........................................................
  //...........................................................
  //...........................................................
  
  // A N A L I S I   V 0
  
  //...........................................................
  //...........................................................
  //...........................................................
     
  if(fESD->GetNumberOfV0s()>0) { 
  
  for(Int_t iV0 = 0; iV0<fESD->GetNumberOfV0s(); iV0++) {
    
     AliESDv0 *v0 = fESD->GetV0(iV0);
     if(!v0) continue;
    
     //..................................................... 
     
     if(v0->GetOnFlyStatus()) continue;
     
     //.....................................................
     
     Int_t pIndex = TMath::Abs(v0->GetPindex()); 
     Int_t nIndex = TMath::Abs(v0->GetNindex()); 
     
     //Double_t mass=v0->GetEffMass();
     
     // negative..................................
     
     Int_t qN, nphN;
     Float_t xN, yN;
     Float_t xpcN, ypcN, thN, phN;
     
     // positive...................................
     
     Int_t qP, nphP;
     Float_t xP, yP;
     Float_t xpcP, ypcP, thP, phP;
     
     //.............................................
     
     
     //cout<<"xP="<<xP<<"yP="<<yP<<"qP="<<qP<<"nphP="<<nphP<<endl;
         
      AliESDtrack *pTrackTest = fESD->GetTrack(pIndex);
      AliESDtrack *nTrackTest = fESD->GetTrack(nIndex);
     
      if (!pTrackTest || !nTrackTest) { Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n"); continue;  }
      
      // Remove like-sign
      if(pTrackTest->GetSign() == nTrackTest->GetSign()){ continue;  } 

      // IMPONGO CHE LE FIGLIE SIANO UNA POSITIVA E L'ALTRA NEGATIVA....
      //if( pTrack->GetSign() <0 || nTrack->GetSign() > 0) continue;
     
      if( pTrackTest->GetSign() ==1){
	
	pTrack = fESD->GetTrack(pIndex);
	nTrack = fESD->GetTrack(nIndex);

	// Daughters' momentum;
	v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
	v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);
	
	}
           
      if( pTrackTest->GetSign() ==-1){
	
	pTrack = fESD->GetTrack(nIndex);
	nTrack = fESD->GetTrack(pIndex);

        // Daughters' momentum;
	v0->GetPPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);
	v0->GetNPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);

	}

      // DCA between daughter and Primary Vertex:
      if (pTrack) lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lMagneticField) );
      if (nTrack) lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lMagneticField) );
      
      // Quality tracks cuts:
      //if ( !(myTracksCuts->IsSelected(pTrack)) || !(myTracksCuts->IsSelected(nTrack)) ) continue;
 
     v0->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lPrimaryVtxPosition[2]);
     v0->GetXYZ(lV0Position[0], lV0Position[1], lV0Position[2]);
     
     //---------------------------------------------------------------------------------------------------------------------------------------------------------
     //  TAGLIO PER RIDURRE IL PESO DEI .ROOT CREATO.................
     
     if(v0->Theta() < 0.873  ||  v0->Theta() > 2.268  )  continue;
     //cout<<"Phi................................"<<  v0->Phi()<<endl;
     if(v0->Phi() > 1.3962  &&  v0->Phi() < 5.9341 ) continue;
     //cout<<"*****************************Phi After................................"<<  v0->Phi()<<endl;
     
     //---------------------------------------------------------------------------------------------------------------------------------------------------------
     
        pTrack->GetHMPIDmip(xP,yP,qP,nphP);
        pTrack->GetHMPIDtrk(xpcP,ypcP,thP,phP);
      
        nTrack->GetHMPIDmip(xN,yN,qN,nphN);
        nTrack->GetHMPIDtrk(xpcN,ypcN,thN,phN);
     
     
	      
     Float_t bp[2];
     Float_t bpCov[3];
     pTrack->GetImpactParameters(bp,bpCov);    

     Float_t bn[2];
     Float_t bnCov[3];
     nTrack->GetImpactParameters(bn,bnCov);    

     
     fVar[0] = pTrack->GetHMPIDcluIdx()/1000000;
     fVar[1] = nTrack->GetHMPIDcluIdx()/1000000;
     
     //fVar[1] = pHmp3;
     fVar[2] = (Float_t)pTrack->P();
     fVar[3] = (Float_t)nTrack->P();
     
     fVar[4] = xpcP;
     fVar[5] = ypcP;
     fVar[6] = xP;
     fVar[7] = yP;
     fVar[8] = (Float_t)pTrack->GetHMPIDsignal();
     fVar[9] = (Float_t)nTrack->GetHMPIDsignal();
     
     fVar[10] = qP;
     fVar[11] = thP;
     fVar[12] = phP;
     
     fVar[13] = (Float_t)pTrack->GetSign();
     fVar[14] = (Float_t)nTrack->GetSign();
     
     fVar[15] = (Float_t)nphP;
     fVar[16] = (Float_t)pTrack->GetNcls(1);
     fVar[17] = (Float_t)nTrack->GetNcls(1);
     
     //fVar[14] = (Float_t)probs[0];
     //fVar[15] = (Float_t)probs[1];
     //fVar[16] = (Float_t)probs[2];
     //fVar[17] = (Float_t)probs[3];
     //fVar[18] = (Float_t)probs[4];
     fVar[18] = (Float_t)pTrack->GetTOFsignal();
     fVar[19] = (Float_t)nTrack->GetTOFsignal();
     
     fVar[20] = (Float_t)pTrack->GetKinkIndex(0);
     fVar[21] = (Float_t)nTrack->GetKinkIndex(0);
     
     fVar[22] = (Float_t)pTrack->Xv();
     fVar[23] = (Float_t)nTrack->Xv();
     
     fVar[24] = (Float_t)pTrack->Yv();
     fVar[25] = (Float_t)nTrack->Yv();
     
     fVar[26] = (Float_t)pTrack->Zv();
     fVar[27] = (Float_t)nTrack->Zv();
     
     fVar[28] = (Float_t)pTrack->GetTPCchi2();
     fVar[29] = (Float_t)nTrack->GetTPCchi2();
     
     // fVar[25] = b[0];
     // fVar[26] = b[1];
     fVar[30] = pTrack->GetHMPIDcluIdx()%1000000/1000;
     fVar[31] = nTrack->GetHMPIDcluIdx()%1000000/1000;
     
     fVar[32] = vtxPos[0];
     fVar[33] = vtxPos[1];
     fVar[34] = vtxPos[2];
//     fVar[31] = (Float_t)ITSrefit;
//     fVar[32] = (Float_t)TPCrefit;
     fVar[35] = (Float_t)pTrack->Eta();
     fVar[36] = (Float_t)nTrack->Eta();
     
     fVar[37] = (Float_t)r[0];
     fVar[38] = (Float_t)r[1];
     fVar[39] = (Float_t)r[2];
     fVar[40] = (Float_t)rout[0];           
     fVar[41] = (Float_t)rout[1];
     fVar[42] = (Float_t)rout[2];
     fVar[43] = pTrack->GetMass();
     fVar[44] = nTrack->GetMass();
     
     fVar[45] = v0->GetDcaV0Daughters(); // OK
     fVar[46] = v0->GetV0CosineOfPointingAngle();
     
     fVar[47] = lMomPos[0];//pTrack->Px();
     fVar[48] = lMomPos[1];//pTrack->Py();
     fVar[49] = lMomPos[2];//pTrack->Pz();
     fVar[50] = lMomNeg[0];//nTrack->Px();
     fVar[51] = lMomNeg[1];//nTrack->Py();
     fVar[52] = lMomNeg[2];//nTrack->Pz();
   
     fVar[53] = v0->P(); // impulso della V0;
     fVar[54] = v0->GetEffMass();
     
     fVar[55] = v0->Xv();
     fVar[56] = v0->Yv();
     fVar[57] = v0->Zv();
     

     fVar[58] = bp[0];
     fVar[59] = bp[1];
     
     fVar[60] = bn[0];
     fVar[61] = bn[1];
     
     fVar[62] = lDcaPosToPrimVertex;
     fVar[63] = lDcaNegToPrimVertex;
     
     fVar[64] = v0->Eta(); // pseudorapidity
     
     fVar[65] = v0->GetChi2V0();
     fVar[66] = v0->GetOnFlyStatus();
     
     fVar[67] = lPrimaryVtxPosition[0];
     fVar[68] = lPrimaryVtxPosition[1];
     fVar[69] = lPrimaryVtxPosition[2];
         
     fVar[70] = v0->GetV0CosineOfPointingAngle(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1], lPrimaryVtxPosition[2]); // new cos pointing angle
     
     
     fVar[71] = lV0Position[0];
     fVar[72] = lV0Position[1];
     fVar[73] = lV0Position[2];
     
     fVar[74] = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);  // v0 radius;
     fVar[75] = TMath::Sqrt(TMath::Power(lV0Position[0] - lPrimaryVtxPosition[0],2) +
		                   TMath::Power(lV0Position[1] - lPrimaryVtxPosition[1],2) +
		                   TMath::Power(lV0Position[2] - lPrimaryVtxPosition[2],2 ));  // v0 decay lenght;
     /*
     fVar[76] = xpcN;
     fVar[77] = ypcN;
     fVar[78] = xN;
     fVar[79] = yN;
     fVar[80] = qN;
     fVar[81] = thN;
     fVar[82] = phN;
     fVar[83] = (Float_t)nphN;
     */
     fVar[84] = v0->Theta();//*TMath::RadToDeg();
     fVar[85] = v0->Phi();//*TMath::RadToDeg();
     
     
     fVar[76] = xpcN;
     fVar[77] = ypcN;
     fVar[78] = xN;
     fVar[79] = yN;
     fVar[80] = qN;
     fVar[81] = thN;
     fVar[82] = phN;
     fVar[83] = (Float_t)nphN;
     
     fTree->Fill();
 
     
     
     
    }  // endl loop sulle V0;
  }
  //fVar[0] = track->GetHMPIDcluIdx()/1000000;
  //fTree->Fill();
  
  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(1,fHmpHistList);
  PostData(2,fTree);
}

//___________________________________________________________________________
void AliHMPIDPIDEfficiencyV0AnalysisTask::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  Info("Terminate"," ");

  if (!fUseMC) return;

  fHmpHistList = dynamic_cast<TList*> (GetOutputData(1));

  if (!fHmpHistList) {
    AliError("Histogram List is not available");
    return;
  }


  AliAnalysisTaskSE::Terminate();

}

//___________________________________________________________________________
void AliHMPIDPIDEfficiencyV0AnalysisTask::UserCreateOutputObjects() {
  //
  //HERE ONE CAN CREATE OUTPUT OBJECTS
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //

  //slot #1
//   OpenFile(1);
   fHmpHistList = new TList();
   fHmpHistList->SetOwner();

   fHmpNevents = new TH1F("fHmpNevents","Number of events",2,0,2);
   fHmpHistList->Add(fHmpNevents);
   
   fThetavsMom = new TH2F("fThetavsMom","Theta Cherenkov vs momentum",250,0,5,250,0,1);
   fHmpHistList->Add(fThetavsMom);
   
   fmass = new TH1D("fmass","mass",400,0,10);
   fHmpHistList->Add(fmass);
   
   massap = new TH1D("massap","massap",400,0,10);
   fHmpHistList->Add(massap);
   
   massan = new TH1D("massan","massan",400,0,10);
   fHmpHistList->Add(massan);
   
   hdiff = new TH1D("hdiff","hdiff",2000,-100,100);
   fHmpHistList->Add(hdiff);
   
   hdiffn = new TH1D("hdiffn","hdiffn",2000,-100,100);
   fHmpHistList->Add(hdiffn);
   
   fmassaHMPID = new TH1D("fmassaHMPID","fmassaHMPID",400,0,10);
   fHmpHistList->Add(fmassaHMPID);
   
   fpdg = new TH1D("fpdg","pdg",5000,0,5000);
   fHmpHistList->Add(fpdg);
   
   farme = new TH2F("farme","p ARMENTEROS vs alfa", 300, -10.5, 10.5, 300, 0, 10.3);
   fHmpHistList->Add(farme);
   
   farmeb = new TH2F("farmeb","p ARMENTEROS vs alfa", 300, -10.5, 10.5, 300, 0, 10.3);
   fHmpHistList->Add(farmeb);
   
   hangle = new TH1D("hangle","pointing angle",1000,-10,10);
   fHmpHistList->Add(hangle);
   
   
   fTree = new TTree("Tree","Tree with data");
   /*
   fTree->Branch("Chamber",&fVar[0]);
   fTree->Branch("pHmp3",&fVar[1]);
   fTree->Branch("P",&fVar[2]);
   fTree->Branch("h",&fVar[3]);
   fTree->Branch("Ypc",&fVar[4]);
   fTree->Branch("X",&fVar[5]);
   fTree->Branch("Y",&fVar[6]);
   fTree->Branch("HMPIDsignal",&fVar[7]);
   fTree->Branch("Charge",&fVar[8]);
   fTree->Branch("Theta",&fVar[9]);
   fTree->Branch("Phi",&fVar[10]);
   
   fTree->Branch("pSign",&fVar[11]);
   fTree->Branch("nSign",&fVar[52]);
   
   fTree->Branch("NumPhotons",&fVar[12]);
   fTree->Branch("NumTPCclust",&fVar[13]);
   fTree->Branch("Prob0",&fVar[14]);
   fTree->Branch("Prob1",&fVar[15]);
   fTree->Branch("Prob2",&fVar[16]);
   fTree->Branch("Prob3",&fVar[17]);
   fTree->Branch("Prob4",&fVar[18]);
   fTree->Branch("TOFsignal",&fVar[19]);
   fTree->Branch("KinkIndex",&fVar[20]);
   fTree->Branch("Xv",&fVar[21]);
   fTree->Branch("Yv",&fVar[22]);
   fTree->Branch("Zv",&fVar[23]);
   fTree->Branch("TPCchi2",&fVar[24]);
   fTree->Branch("b0",&fVar[25]);
   fTree->Branch("b1",&fVar[26]);
   fTree->Branch("ClustSize",&fVar[27]);
   fTree->Branch("PrimVertexX",&fVar[28]);
   fTree->Branch("PrimVertexY",&fVar[29]);
   fTree->Branch("PrimVertexZ",&fVar[30]);
   fTree->Branch("ITSrefit",&fVar[31]);
   fTree->Branch("TPCrefit",&fVar[32]);   
   fTree->Branch("Eta",&fVar[33]);
   fTree->Branch("xTrack",&fVar[34]);
   fTree->Branch("yTrack",&fVar[35]);
   fTree->Branch("zTrack",&fVar[36]);
   fTree->Branch("xOuterTrack",&fVar[37]);
   fTree->Branch("yOuterTrack",&fVar[38]);
   fTree->Branch("zOuterTrack",&fVar[39]);
   fTree->Branch("massaHMPID",&fVar[40]);
   fTree->Branch("massa_pos",&fVar[41]);
   fTree->Branch("massa_neg",&fVar[42]);
   fTree->Branch("alfa",&fVar[43]);
   fTree->Branch("cosPos",&fVar[44]);
   fTree->Branch("cosNeg",&fVar[45]);
   fTree->Branch("pmadre",&fVar[46]);
   fTree->Branch("distance",&fVar[47]);
   fTree->Branch("pARMENTEROSa",&fVar[48]);
   fTree->Branch("pARMENTEROSb",&fVar[49]);
   fTree->Branch("pointingangle",&fVar[50]);
   fTree->Branch("dcafiglie",&fVar[51]);
   fTree->Branch("alfa2",&fVar[52]);
   fTree->Branch("qt",&fVar[53]);
   fTree->Branch("pMadreSomma",&fVar[54]);
   
   fTree->Branch("pxp",&fVar[55]);
   fTree->Branch("pyp",&fVar[56]);
   fTree->Branch("pzp",&fVar[57]);
   fTree->Branch("pxn",&fVar[58]);
   fTree->Branch("pyn",&fVar[59]);
   fTree->Branch("pzn",&fVar[60]);
   */
   
       
   //==========================================================================================    
       
       
   fTree->Branch("Chamberp",&fVar[0]);
   fTree->Branch("Chambern",&fVar[1]);
   fTree->Branch("Pp",&fVar[2]);
   fTree->Branch("Pn",&fVar[3]);
   fTree->Branch("HMPIDsignalp",&fVar[8]);
   fTree->Branch("HMPIDsignaln",&fVar[9]);
   fTree->Branch("GetNclsp",&fVar[16]);
   fTree->Branch("GetNclsn",&fVar[17]);
   fTree->Branch("TOFp",&fVar[18]);
   fTree->Branch("TOFn",&fVar[19]);
   fTree->Branch("kinkp",&fVar[20]);
   fTree->Branch("kinkn",&fVar[21]);
   fTree->Branch("Xvp",&fVar[22]);
   fTree->Branch("Xvn",&fVar[23]);
   fTree->Branch("Yvp",&fVar[24]);
   fTree->Branch("Yvn",&fVar[25]);
   fTree->Branch("Zvp",&fVar[26]);
   fTree->Branch("Zvn",&fVar[27]);
   
   fTree->Branch("TPCp",&fVar[28]);
   fTree->Branch("TPCn",&fVar[29]);
   fTree->Branch("CHp",&fVar[30]);
   fTree->Branch("CHn",&fVar[31]);
   fTree->Branch("vtxpos",&fVar[32]);
   fTree->Branch("vtypos",&fVar[33]);
   fTree->Branch("vtzpos",&fVar[34]);
   fTree->Branch("Etap",&fVar[35]);
   fTree->Branch("Etan",&fVar[36]);
   fTree->Branch("R0",&fVar[37]);
   fTree->Branch("R1",&fVar[38]);
   fTree->Branch("R2",&fVar[39]);
   fTree->Branch("ROUT0",&fVar[40]);
   fTree->Branch("ROUT1",&fVar[41]);
   fTree->Branch("ROUT2",&fVar[42]);
   fTree->Branch("massp",&fVar[43]);
   fTree->Branch("massn",&fVar[44]);
   fTree->Branch("dcafiglie",&fVar[45]);
   fTree->Branch("pointingangle",&fVar[46]);
   fTree->Branch("pxp",&fVar[47]);
   fTree->Branch("pyp",&fVar[48]);
   fTree->Branch("pzp",&fVar[49]);
   fTree->Branch("pxn",&fVar[50]);
   fTree->Branch("pyn",&fVar[51]);
   fTree->Branch("pzn",&fVar[52]);
   fTree->Branch("Pmadre",&fVar[53]);
   fTree->Branch("MassEffic",&fVar[54]);
   
   fTree->Branch("v0x",&fVar[55]);
   fTree->Branch("v0y",&fVar[56]);
   fTree->Branch("v0z",&fVar[57]);
   
   fTree->Branch("bpos0",&fVar[58]);
   fTree->Branch("bpos1",&fVar[59]);
   fTree->Branch("bneg0",&fVar[60]);
   fTree->Branch("bneg1",&fVar[61]);
   
   fTree->Branch("DcaPosToPrimVertex",&fVar[62]);
   fTree->Branch("DcaNegToPrimVertex",&fVar[63]);
   
   fTree->Branch("etaV0",&fVar[64]);
   fTree->Branch("Chi2V0",&fVar[65]);
   fTree->Branch("ONFLY",&fVar[66]);
   
   fTree->Branch("PrimVertex0",&fVar[67]);
   fTree->Branch("PrimVertex1",&fVar[68]);
   fTree->Branch("PrimVertex2",&fVar[69]);
   
   fTree->Branch("NewCosPointingAngle",&fVar[70]);
   
   fTree->Branch("V0position0",&fVar[71]);
   fTree->Branch("V0position1",&fVar[72]);
   fTree->Branch("V0position2",&fVar[73]);
     
   fTree->Branch("V0radius",&fVar[74]);
   fTree->Branch("V0decayLenght",&fVar[75]);
   // 8 POSITIVE.......................
   fTree->Branch("XpcP",&fVar[4]);
   fTree->Branch("YpcP",&fVar[5]);
   fTree->Branch("XP"  ,&fVar[6]);
   fTree->Branch("YP"  ,&fVar[7]);
   fTree->Branch("QP"  ,&fVar[10]);
   fTree->Branch("THP" ,&fVar[11]);
   fTree->Branch("PHP" ,&fVar[12]);
   fTree->Branch("NPHP",&fVar[15]);
   // 8 NEGATIVE.......................
   
   
   
   //==========================================================================================    
   
   /*
   fTree->Branch("XpcN  ",&fVar[76]);
   fTree->Branch("YpcN  ",&fVar[77]);
   fTree->Branch("XN    ",&fVar[78]);
   fTree->Branch("YN    ",&fVar[79]);
   fTree->Branch("QN    ",&fVar[80]);
   fTree->Branch("THN   ",&fVar[81]);
   fTree->Branch("PHN   ",&fVar[82]);
   fTree->Branch("NPHN  ",&fVar[83]);
   */
   fTree->Branch("XpcN",&fVar[76]);
   fTree->Branch("YpcN",&fVar[77]);
   fTree->Branch("XN",&fVar[78]);
   fTree->Branch("YN",&fVar[79]);
   fTree->Branch("QN",&fVar[80]);
   fTree->Branch("THN",&fVar[81]);
   fTree->Branch("PHN",&fVar[82]);
   fTree->Branch("NPHN",&fVar[83]);
   
   
   fTree->Branch("Theta_VO",&fVar[84]);
   fTree->Branch("Phi_VO",&fVar[85]);
  
   
          
   PostData(1,fHmpHistList);
   PostData(2,fTree);
}

//____________________________________________________________________________________________________________________________________
Bool_t AliHMPIDPIDEfficiencyV0AnalysisTask::Equal(Double_t x, Double_t y, Double_t tolerance)
{
 return abs(x - y) <= tolerance ;
}
   
#endif
