
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

//-----------------------------------------------------------------
//         AliSpectraAODEventCuts class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TSpline.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDVZERO.h"
#include "AliAODVZERO.h"
#include "AliSpectraAODEventCuts.h"
#include "AliSpectraAODTrackCuts.h"
#include <iostream>

using namespace std;

ClassImp(AliSpectraAODEventCuts)

AliSpectraAODEventCuts::AliSpectraAODEventCuts(const char *name) : 
  TNamed(name, "AOD Event Cuts"),
  fAOD(0),
  fSelectBit(AliVEvent::kMB),
  fCentralityMethod("V0M"),
  fTrackBits(1),
  fIsMC(0),
  fIsLHC10h(1),
  fTrackCuts(0),
  fIsSelected(0),
  fCentralityCutMin(0.),
  fCentralityCutMax(999),
  fQVectorCutMin(-999.),
  fQVectorCutMax(999.),
  fVertexCutMin(-10.),
  fVertexCutMax(10.),
  fMultiplicityCutMin(-999.),
  fMultiplicityCutMax(99999.),
  fqV0C(-999.),
  fqV0A(-999.),
  fqV0Cx(-999.),
  fqV0Ax(-999.),
  fqV0Cy(-999.),
  fqV0Ay(-999.),
  fPsiV0C(-999.),
  fPsiV0A(-999.),
  fCent(-999.),
  fOutput(0),
  fCalib(0),
  fRun(-1),
  fMultV0(0),
  fV0Cpol1(-1),
  fV0Cpol2(-1),
  fV0Cpol3(-1),
  fV0Cpol4(-1),
  fV0Apol1(-1),
  fV0Apol2(-1),
  fV0Apol3(-1),
  fV0Apol4(-1),
  fQvecIntList(0),
  fQvecIntegral(0), 
  fSplineArrayV0A(0),
  fSplineArrayV0C(0)
{
  // Constructor
  fOutput=new TList();
  fOutput->SetOwner();
  fOutput->SetName("fOutput");
  
  fCalib=new TList();
  fCalib->SetOwner();
  fCalib->SetName("fCalib");
  
  fQvecIntList=new TList();
  fQvecIntList->SetOwner();
  fQvecIntList->SetName("fQvecIntList");
  
  TH1I *fHistoCuts = new TH1I("fHistoCuts", "Event Cuts", kNVtxCuts, -0.5, kNVtxCuts - 0.5);
  TH1F *fHistoVtxBefSel = new TH1F("fHistoVtxBefSel", "Vtx distr before event selection;z (cm)",500,-15,15);
  TH1F *fHistoVtxAftSel = new TH1F("fHistoVtxAftSel", "Vtx distr after event selection;z (cm)",500,-15,15);
  TH1F *fHistoEtaBefSel = new TH1F("fHistoEtaBefSel", "Eta distr before event selection;eta",500,-2,2);
  TH1F *fHistoEtaAftSel = new TH1F("fHistoEtaAftSel", "Eta distr after event selection;eta",500,-2,2);
  TH1F *fHistoNChAftSel = new TH1F("fHistoNChAftSel", "NCh distr after event selection;Nch",2000,-0.5,1999.5);
  TH2F *fHistoQVector = new TH2F("fHistoQVector", "QVector with VZERO distribution;centrality;Q vector from EP task",20,0,100,100,0,10);
  TH2F *fHistoEP = new TH2F("fHistoEP", "EP with VZERO distribution;centrality;Psi_{EP} from EP task",20,0,100,100,-2,2);
  TH2F *fPsiACor = new TH2F("fPsiACor", "EP with VZERO A distribution;centrality;Psi_{EP} VZERO-A",20,0,100,100,-2,2);
  TH2F *fPsiCCor = new TH2F("fPsiCCor", "EP with VZERO C distribution;centrality;Psi_{EP} VZERO-C",20,0,100,100,-2,2);
  TH2F *fQVecACor = new TH2F("fQVecACor", "QVec VZERO A;centrality;Qvector VZERO-A",20,0,100,100,0,10);
  TH2F *fQVecCCor = new TH2F("fQVecCCor", "QVec VZERO C;centrality;Qvector VZERO-C",20,0,100,100,0,10);
  TH2F *fV0M = new TH2F("fV0M", "V0 Multiplicity, before correction;V0 sector",64,-.5,63.5,500,0,1000);
  TH2F *fV0MCor = new TH2F("fV0MCor", "V0 Multiplicity, after correction;V0 sector",64,-.5,63.5,500,0,1000);
  
  fSplineArrayV0A = new TObjArray();
  fSplineArrayV0A->SetOwner();
  fSplineArrayV0C = new TObjArray();
  fSplineArrayV0C->SetOwner();
  
  fOutput->Add(fHistoCuts);
  fOutput->Add(fHistoVtxBefSel);
  fOutput->Add(fHistoVtxAftSel);
  fOutput->Add(fHistoEtaBefSel);
  fOutput->Add(fHistoEtaAftSel);
  fOutput->Add(fHistoNChAftSel);
  fOutput->Add(fHistoQVector);
  fOutput->Add(fHistoEP);
  fOutput->Add(fPsiACor);
  fOutput->Add(fPsiCCor);
  fOutput->Add(fQVecACor);
  fOutput->Add(fQVecCCor);
  fOutput->Add(fV0M);
  fOutput->Add(fV0MCor);
  
  for (Int_t i = 0; i<10; i++){
    fMeanQxa2[i] = -1;
    fMeanQya2[i] = -1;
    fMeanQxc2[i] = -1;
    fMeanQyc2[i] = -1;   
  }
}

//______________________________________________________
Bool_t AliSpectraAODEventCuts::IsSelected(AliAODEvent * aod,AliSpectraAODTrackCuts     *trackcuts)
{
  // Returns true if Event Cuts are selected and applied
  fAOD = aod;
  fTrackCuts = trackcuts; // FIXME: if track cuts is 0, do not set (use the pre-set member). Do we need to pass this here at all??
  // FIXME: all those references by name are slow.
  ((TH1I*)fOutput->FindObject("fHistoCuts"))->Fill(kProcessedEvents);
  Bool_t IsPhysSel = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fSelectBit);
  if(!IsPhysSel)return IsPhysSel;
  ((TH1I*)fOutput->FindObject("fHistoCuts"))->Fill(kPhysSelEvents);
  //loop on tracks, before event selection, filling QA histos
  AliAODVertex * vertex = fAOD->GetPrimaryVertex();//FIXME vertex is recreated
  if(vertex)((TH1F*)fOutput->FindObject("fHistoVtxBefSel"))->Fill(vertex->GetZ());
  fIsSelected =kFALSE;
  if(CheckVtxRange() && CheckCentralityCut() && CheckMultiplicityCut()){ //selection on vertex and Centrality
    fIsSelected=CheckQVectorCut(); // QVector is calculated only if the centrality and vertex are correct (performance)
  }
  if(fIsSelected){
    ((TH1I*)fOutput->FindObject("fHistoCuts"))->Fill(kAcceptedEvents);
    if(vertex)((TH1F*)fOutput->FindObject("fHistoVtxAftSel"))->Fill(vertex->GetZ());
  }
  Int_t Nch=0;
  for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
    AliAODTrack* track = fAOD->GetTrack(iTracks);
    if (!fTrackCuts->IsSelected(track,kFALSE)) continue;
    ((TH1F*)fOutput->FindObject("fHistoEtaBefSel"))->Fill(track->Eta());
    if(fIsSelected){
      ((TH1F*)fOutput->FindObject("fHistoEtaAftSel"))->Fill(track->Eta());
      Nch++;
    }
  }
  if(fIsSelected)((TH1F*)fOutput->FindObject("fHistoNChAftSel"))->Fill(Nch);
  return fIsSelected;
}

//______________________________________________________
Bool_t AliSpectraAODEventCuts::CheckVtxRange()
{
  // reject events outside of range
  AliAODVertex * vertex = fAOD->GetPrimaryVertex();
  //when moving to 2011 wìone has to add a cut using SPD vertex.
  //The point is that for events with |z|>20 the vertexer tracks is not working (only 2011!). One has to put a safety cut using SPD vertex large e.g. 15cm
  if (!vertex)
    {
      ((TH1I*)fOutput->FindObject("fHistoCuts"))->Fill(kVtxNoEvent);
      return kFALSE;
    }
  if (vertex->GetZ() > fVertexCutMin && vertex->GetZ() < fVertexCutMax)
    {
      return kTRUE;
    }
  ((TH1I*)fOutput->FindObject("fHistoCuts"))->Fill(kVtxRange);
  return kFALSE;
}

//______________________________________________________
Bool_t AliSpectraAODEventCuts::CheckCentralityCut()
{
  // Check centrality cut
  fCent=-999.;
  fCent=fAOD->GetCentrality()->GetCentralityPercentile(fCentralityMethod.Data());
  if ( (fCent <= fCentralityCutMax)  &&  (fCent >= fCentralityCutMin) )  return kTRUE;   
  ((TH1I*)fOutput->FindObject("fHistoCuts"))->Fill(kVtxCentral);
  return kFALSE;
}

//______________________________________________________
Bool_t AliSpectraAODEventCuts::CheckMultiplicityCut()
{
  // Check multiplicity cut
  // FIXME: why this is not tracket in the track stats histos?
  Int_t Ncharged=0;
  for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++){
    AliAODTrack* track = fAOD->GetTrack(iTracks);
    if (!fTrackCuts->IsSelected(track,kFALSE)) continue;
    Ncharged++;
  }
  if(Ncharged>fMultiplicityCutMin && Ncharged<fMultiplicityCutMax)return kTRUE;
  
  return kFALSE;
}

//______________________________________________________
Bool_t AliSpectraAODEventCuts::CheckQVectorCut()
{ 
  Double_t qxEPVZERO = -999., qyEPVZERO = -999.;
  Double_t qVZERO = -999.;
  Double_t psi=-999.;
  
  if(fIsLHC10h){
    qVZERO=CalculateQVectorLHC10h();
    psi=fPsiV0A;
  }else{
    psi=fAOD->GetEventplane()->CalculateVZEROEventPlane(fAOD,10,2,qxEPVZERO,qyEPVZERO);//FIXME we can a flag for 2010 and 2011
    qVZERO= TMath::Sqrt(qxEPVZERO*qxEPVZERO + qyEPVZERO*qyEPVZERO);
  }
  
  //cut on q vector
  if(qVZERO<fQVectorCutMin || qVZERO>fQVectorCutMax)return kFALSE;
  Double_t cent=fAOD->GetCentrality()->GetCentralityPercentile(fCentralityMethod.Data());
  ((TH2F*)fOutput->FindObject("fHistoQVector"))->Fill(cent,qVZERO);
  ((TH2F*)fOutput->FindObject("fHistoEP"))->Fill(cent,psi);
  ((TH1I*)fOutput->FindObject("fHistoCuts"))->Fill(kQVector);
  

  return kTRUE;
}

//______________________________________________________
Double_t AliSpectraAODEventCuts::CalculateQVectorLHC10h(){
  
  Int_t run = fAOD->GetRunNumber();
  if(run != fRun){
    // Load the calibrations run dependent
    if(OpenInfoCalbration(run))fRun=run;
    else{
      fqV0C=-999.;
      fqV0A=-999.;
      return -999.;
    }
  }
  
  //V0 info    
  Double_t Qxa2 = 0, Qya2 = 0;
  Double_t Qxc2 = 0, Qyc2 = 0;
  Double_t sumMc = 0, sumMa = 0;
  
  AliAODVZERO* aodV0 = fAOD->GetVZEROData();
  
  for (Int_t iv0 = 0; iv0 < 64; iv0++) {
    
    Double_t phiV0 = TMath::PiOver4()*(0.5 + iv0 % 8);
    
    Float_t multv0 = aodV0->GetMultiplicity(iv0);
    ((TH2F*)fOutput->FindObject("fV0M"))->Fill(iv0,multv0);
    
    if (iv0 < 32){
      
      Double_t multCorC = -10;
      
      if (iv0 < 8)
	multCorC = multv0*fV0Cpol1/fMultV0->GetBinContent(iv0+1);
      else if (iv0 >= 8 && iv0 < 16)
	multCorC = multv0*fV0Cpol2/fMultV0->GetBinContent(iv0+1);
      else if (iv0 >= 16 && iv0 < 24)
	multCorC = multv0*fV0Cpol3/fMultV0->GetBinContent(iv0+1);
      else if (iv0 >= 24 && iv0 < 32)
	multCorC = multv0*fV0Cpol4/fMultV0->GetBinContent(iv0+1);
      
      if (multCorC < 0){
	cout<<"Problem with multiplicity in V0C"<<endl;
	fqV0C=-999.;
	fqV0A=-999.;
	return -999.;
      }
      
      Qxc2 += TMath::Cos(2*phiV0) * multCorC;
      Qyc2 += TMath::Sin(2*phiV0) * multCorC;
      
      sumMc += multCorC;
      ((TH2F*)fOutput->FindObject("fV0MCor"))->Fill(iv0,multCorC);
      
    } else {
      
      Double_t multCorA = -10;
      
      if (iv0 >= 32 && iv0 < 40)
	multCorA = multv0*fV0Apol1/fMultV0->GetBinContent(iv0+1);
      else if (iv0 >= 40 && iv0 < 48)
	multCorA = multv0*fV0Apol2/fMultV0->GetBinContent(iv0+1);
      else if (iv0 >= 48 && iv0 < 56)
	multCorA = multv0*fV0Apol3/fMultV0->GetBinContent(iv0+1);
      else if (iv0 >= 56 && iv0 < 64)
	multCorA = multv0*fV0Apol4/fMultV0->GetBinContent(iv0+1);
      
      if (multCorA < 0){
	cout<<"Problem with multiplicity in V0A"<<endl;
	fqV0C=-999.;
	fqV0A=-999.;
	return -999.;
      }
      
      Qxa2 += TMath::Cos(2*phiV0) * multCorA;
      Qya2 += TMath::Sin(2*phiV0) * multCorA;
      
      sumMa += multCorA;
      ((TH2F*)fOutput->FindObject("fV0MCor"))->Fill(iv0,multCorA);
    }
  }
  
  Short_t centrV0  = GetCentrCode(fAOD);
  
  Double_t Qxamean2 = 0.;
  Double_t Qyamean2 = 0.;
  Double_t Qxcmean2 = 0.;
  Double_t Qycmean2 = 0.;
  
  if(centrV0!=-1){
    Qxamean2 = fMeanQxa2[centrV0];
    Qyamean2 = fMeanQya2[centrV0];
    Qxcmean2 = fMeanQxc2[centrV0];
    Qycmean2 = fMeanQyc2[centrV0];
  }
  
  Double_t QxaCor2 = Qxa2 - Qxamean2*sumMa;
  Double_t QyaCor2 = Qya2 - Qyamean2*sumMa;
  Double_t QxcCor2 = Qxc2 - Qxcmean2*sumMc;
  Double_t QycCor2 = Qyc2 - Qycmean2*sumMc;
  
  fPsiV0A = TMath::ATan2(QyaCor2, QxaCor2)/2.;
  fPsiV0C = TMath::ATan2(QycCor2, QxcCor2)/2.;
  
  ((TH2F*)fOutput->FindObject("fPsiACor"))->Fill((Float_t)fAOD->GetCentrality()->GetCentralityPercentile("V0M"), fPsiV0A);
  ((TH2F*)fOutput->FindObject("fPsiCCor"))->Fill((Float_t)fAOD->GetCentrality()->GetCentralityPercentile("V0M"), fPsiV0C);
  
  fqV0A = TMath::Sqrt((QxaCor2*QxaCor2 + QyaCor2*QyaCor2)/sumMa);
  fqV0C = TMath::Sqrt((QxcCor2*QxcCor2 + QycCor2*QycCor2)/sumMc);
  fqV0Ax = QxaCor2*TMath::Sqrt(1./sumMa);
  fqV0Cx = QxcCor2*TMath::Sqrt(1./sumMc);
  fqV0Ay = QyaCor2*TMath::Sqrt(1./sumMa);
  fqV0Cy = QycCor2*TMath::Sqrt(1./sumMc);
  
  ((TH2F*)fOutput->FindObject("fQVecACor"))->Fill((Float_t)fAOD->GetCentrality()->GetCentralityPercentile("V0M"), fqV0A);
  ((TH2F*)fOutput->FindObject("fQVecCCor"))->Fill((Float_t)fAOD->GetCentrality()->GetCentralityPercentile("V0M"), fqV0C);
  
  return fqV0A; //FIXME we have to combine VZERO-A and C
}

//______________________________________________________
Short_t AliSpectraAODEventCuts::GetCentrCode(AliVEvent* ev)
{
  
  Short_t centrCode = -1;
  
  AliCentrality* centrality = 0;
  AliAODEvent* aod = (AliAODEvent*)ev;  
  centrality = aod->GetHeader()->GetCentralityP();
  
  Float_t centV0 = centrality->GetCentralityPercentile("V0M");
  Float_t centTrk = centrality->GetCentralityPercentile("TRK"); 
  
  
  if (TMath::Abs(centV0 - centTrk) < 5.0 && centV0 <= 80 && centV0 > 0){  
    
     if ((centV0 > 0) && (centV0 <= 5.0))
       centrCode = 0; 
     else if ((centV0 > 5.0) && (centV0 <= 10.0))
       centrCode = 1;
     else if ((centV0 > 10.0) && (centV0 <= 20.0))
       centrCode = 2;
     else if ((centV0 > 20.0) && (centV0 <= 30.0))
       centrCode = 3;   
     else if ((centV0 > 30.0) && (centV0 <= 40.0))
       centrCode = 4; 
     else if ((centV0 > 40.0) && (centV0 <= 50.0))
       centrCode = 5;  
     else if ((centV0 > 50.0) && (centV0 <= 60.0))
       centrCode = 6;
     else if ((centV0 > 60.0) && (centV0 <= 70.0))
       centrCode = 7;         			       
     else if ((centV0 > 70.0) && (centV0 <= 80.0))
       centrCode = 8;  
  }
  
  return centrCode;
  
}

//______________________________________________________
void AliSpectraAODEventCuts::PrintCuts()
{
  // print info about event cuts
  cout << "Event Stats" << endl;
  cout << " > Trigger Selection: " << fSelectBit << endl;
  cout << " > Centrality estimator: " << fCentralityMethod << endl;
  cout << " > Number of accepted events: " << NumberOfEvents() << endl;
  cout << " > Number of processed events: " << NumberOfProcessedEvents() << endl;
  cout << " > Number of PhysSel events: " <<  NumberOfPhysSelEvents()  << endl;
  cout << " > Vertex out of range: " << ((TH1I*)fOutput->FindObject("fHistoCuts"))->GetBinContent(kVtxRange + 1) << endl;
  cout << " > Events cut by centrality: " << ((TH1I*)fOutput->FindObject("fHistoCuts"))->GetBinContent(kVtxCentral + 1) << endl;
  cout << " > Events without vertex: " << ((TH1I*)fOutput->FindObject("fHistoCuts"))->GetBinContent(kVtxNoEvent + 1) << endl;
  cout << " > QVector cut: " << ((TH1I*)fOutput->FindObject("fHistoCuts"))->GetBinContent(kQVector + 1) << endl;
  cout << " > Track type used for the QVector calculation: " << fTrackBits << endl;
  cout << " > QRange: [" << fQVectorCutMin <<"," <<fQVectorCutMax<<"]"<< endl;
  cout << " > Vertex: [" << fVertexCutMin <<"," <<fVertexCutMax<<"]"<< endl;
  cout << " > Multiplicity: [" << fMultiplicityCutMin <<"," <<fMultiplicityCutMax<<"]"<< endl;
  cout << " > Centrality: [" << fCentralityCutMin <<"," <<fCentralityCutMax<<"]"<< endl;
}

//_____________________________________________________________________________
Bool_t AliSpectraAODEventCuts::OpenInfoCalbration(Int_t run)
{
  
  AliOADBContainer* cont = (AliOADBContainer*) fCalib->FindObject("hMultV0BefCorr");
  if(!cont){
    printf("OADB object hMultV0BefCorr is not available in the file\n");
    return 0;	
  }
  
  if(!(cont->GetObject(run))){
    printf("OADB object hMultV0BefCorr is not available for run %i\n",run);
    return 0;	
  }
  fMultV0 = ((TH2F*) cont->GetObject(run))->ProfileX();
  
  TF1* fpolc1 = new TF1("fpolc1","pol0", 0, 7); 
  fMultV0->Fit(fpolc1, "RN");
  fV0Cpol1 = fpolc1->GetParameter(0);
  
  TF1* fpolc2 = new TF1("fpolc2","pol0", 8, 15); 
  fMultV0->Fit(fpolc2, "RN");
  fV0Cpol2 = fpolc2->GetParameter(0);
  
  TF1* fpolc3 = new TF1("fpolc3","pol0", 16, 23); 
  fMultV0->Fit(fpolc3, "RN");
  fV0Cpol3 = fpolc3->GetParameter(0);
  
  TF1* fpolc4 = new TF1("fpolc4","pol0", 24, 31); 
  fMultV0->Fit(fpolc4, "RN");
  fV0Cpol4 = fpolc4->GetParameter(0);
  
  TF1* fpola1 = new TF1("fpola1","pol0", 32, 39);
  fMultV0->Fit(fpola1, "RN");
  fV0Apol1 = fpola1->GetParameter(0);

  TF1* fpola2 = new TF1("fpola2","pol0", 40, 47);
  fMultV0->Fit(fpola2, "RN");
  fV0Apol2 = fpola2->GetParameter(0);
  
  TF1* fpola3 = new TF1("fpola3","pol0", 48, 55);
  fMultV0->Fit(fpola3, "RN");
  fV0Apol3 = fpola3->GetParameter(0);
  
  TF1* fpola4 = new TF1("fpola4","pol0", 56, 63);
  fMultV0->Fit(fpola4, "RN");
  fV0Apol4 = fpola4->GetParameter(0);

  for(Int_t i=0; i < 10; i++){

    char nameQxa2[100];
    snprintf(nameQxa2,100, "hQxa2m_%i", i);

    char nameQya2[100];
    snprintf(nameQya2,100, "hQya2m_%i", i);

    char nameQxc2[100];
    snprintf(nameQxc2,100, "hQxc2m_%i", i);

    char nameQyc2[100];
    snprintf(nameQyc2,100, "hQyc2m_%i", i);
	
    AliOADBContainer* contQxa2 = (AliOADBContainer*) fCalib->FindObject(nameQxa2);
    if(!contQxa2){
      printf("OADB object %s is not available in the file\n", nameQxa2);
      return 0;	
    }
	
    if(!(contQxa2->GetObject(run))){
      printf("OADB object %s is not available for run %i\n", nameQxa2, run);
      return 0;	
    }

    fMeanQxa2[i] = ((TH1F*) contQxa2->GetObject(run))->GetMean();

 
    AliOADBContainer* contQya2 = (AliOADBContainer*) fCalib->FindObject(nameQya2);
    if(!contQya2){
      printf("OADB object %s is not available in the file\n", nameQya2);
      return 0;	
    }
	
    if(!(contQya2->GetObject(run))){
      printf("OADB object %s is not available for run %i\n", nameQya2, run);
      return 0;	
    }

    fMeanQya2[i] = ((TH1F*) contQya2->GetObject(run))->GetMean();


    AliOADBContainer* contQxc2 = (AliOADBContainer*) fCalib->FindObject(nameQxc2);
    if(!contQxc2){
      printf("OADB object %s is not available in the file\n", nameQxc2);
      return 0;	
    }
	
    if(!(contQxc2->GetObject(run))){
      printf("OADB object %s is not available for run %i\n", nameQxc2, run);
      return 0;	
    }

    fMeanQxc2[i] = ((TH1F*) contQxc2->GetObject(run))->GetMean();
 

    AliOADBContainer* contQyc2 = (AliOADBContainer*) fCalib->FindObject(nameQyc2);
    if(!contQyc2){
      printf("OADB object %s is not available in the file\n", nameQyc2);
      return 0;	
    }
	
    if(!(contQyc2->GetObject(run))){
      printf("OADB object %s is not available for run %i\n", nameQyc2, run);
      return 0;	
    }

    fMeanQyc2[i] = ((TH1F*) contQyc2->GetObject(run))->GetMean();

  }
  return 1;
}

//______________________________________________________


Long64_t AliSpectraAODEventCuts::Merge(TCollection* list)
{
  // Merge a list of AliSpectraAODEventCuts objects with this.
  // Returns the number of merged objects (including this).
  
  AliInfo("Merging");
  
  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;
  
  TIterator* iter = list->MakeIterator();
  TObject* obj;
  
  // collections of all histograms
  TList collections;
  
  Int_t count = 0;

  while ((obj = iter->Next())) {
    AliSpectraAODEventCuts* entry = dynamic_cast<AliSpectraAODEventCuts*> (obj);
    if (entry == 0) 
      continue;

    TList * l = entry->GetOutputList();      
    collections.Add(l);
    count++;
  }
  
  fOutput->Merge(&collections);
  
  delete iter;

  return count+1;
}

//______________________________________________________
Double_t AliSpectraAODEventCuts::GetQvecPercentile(Int_t v0side){
 
  //check Qvec and Centrality consistency
  if(fCent>90.) return -999.;
  if(v0side==0/*V0A*/){ if(fqV0A== -999.) return -999.; }
  if(v0side==1/*V0C*/){ if(fqV0C== -999.) return -999.; }
  
  fQvecIntegral = 0x0;

  if(v0side==0/*V0A*/){ fQvecIntegral = (TH2D*)fQvecIntList->FindObject("VZEROA"); }
  if(v0side==1/*V0C*/){ fQvecIntegral = (TH2D*)fQvecIntList->FindObject("VZEROC"); }

  Double_t ic = (Int_t)fCent; //fQvecIntegral: 1% centrality bin
  
  TH1D *h1D = (TH1D*)fQvecIntegral->ProjectionY("h1D",ic+1,ic+1);
  
  TSpline *spline = 0x0;
  
  if(v0side==0/*V0A*/){
    if( CheckSplineArray(fSplineArrayV0A) ) {
      spline = (TSpline*)fSplineArrayV0A->At(ic);
      //cout<<"Qvec V0A - ic: "<<ic<<" - Found TSpline..."<<endl;
    } else {
      spline = new TSpline3(h1D,"sp3");
      fSplineArrayV0A->AddAtAndExpand(spline,ic);
      //cout<<"Qvec V0A - ic: "<<ic<<" - TSpline not found. Creating..."<<endl;
    }
  }
  else if(v0side==1/*V0C*/){
    if( CheckSplineArray(fSplineArrayV0C) ) {
      spline = (TSpline*)fSplineArrayV0C->At(ic);
    } else {
      spline = new TSpline3(h1D,"sp3");
      fSplineArrayV0C->AddAtAndExpand(spline,ic);
    }
  }
  
  Double_t percentile =  -999.;
  if(v0side==0/*V0A*/){ percentile = 100*spline->Eval(fqV0A); }
  if(v0side==1/*V0C*/){ percentile = 100*spline->Eval(fqV0C); }
  
  if(percentile>100. || percentile<0.) return -999.;

  return percentile;
}

//______________________________________________________
Bool_t AliSpectraAODEventCuts::CheckSplineArray(TObjArray * splarr){
  //check TSpline array consistency
  
  Int_t n = (Int_t)fCent;
  
  if(splarr->GetSize()<n) return kFALSE;
  
  if(!splarr->At(n)) return kFALSE;
    
  return kTRUE;
  
}
