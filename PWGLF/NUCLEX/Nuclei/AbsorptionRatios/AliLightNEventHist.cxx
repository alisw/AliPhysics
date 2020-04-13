/*
 * AliLightNEventHist.cxx
 *
 *  Created on: Nov 22,2017
 *      Author: gu74req
 */

#include "AliLightNEventHist.h"
ClassImp(AliLightNEventHist)
AliLightNEventHist::AliLightNEventHist()
{
    fEventCutList=new TList();
    fEventCutList->SetName("Event Cuts");
    fEventCutList->SetOwner();
    
    fEvtCounter=new TH1F("EventCounter","Event Counter",10,0,10);
    fEvtCounter->GetXaxis()->SetBinLabel(1,"Events");
    fEvtCounter->GetXaxis()->SetBinLabel(2,"AliEventCuts");
    fEvtCounter->GetXaxis()->SetBinLabel(3,"Phys. Sel.");
    fEvtCounter->GetXaxis()->SetBinLabel(4,"nContrib");
    fEvtCounter->GetXaxis()->SetBinLabel(5,"zVtx");
    fEvtCounter->GetXaxis()->SetBinLabel(6,"PileUp");
    fEvtCounter->GetXaxis()->SetBinLabel(7,"SPD Mult Cleanup");
    fEvtCounter->GetXaxis()->SetBinLabel(8,"V0A Mult Cleanup");
    fEvtCounter->GetXaxis()->SetBinLabel(9,"V0C Mult Cleanup");
    fEvtCounter->GetXaxis()->SetBinLabel(10,"RefMult08 Cleanup");
    fEventCutList->Add(fEvtCounter);
    
    fCutConfig=new TProfile("CutConfig","Cut Config",20,0,20);
    fCutConfig->GetXaxis()->SetBinLabel(1,"Min Contrib");
    fCutConfig->GetXaxis()->SetBinLabel(2,"CutZvtx");
    fCutConfig->GetXaxis()->SetBinLabel(3,"Min Zvtx");
    fCutConfig->GetXaxis()->SetBinLabel(4,"Max ZVtx");
    fCutConfig->GetXaxis()->SetBinLabel(5,"PileUp Rejection");
    fCutConfig->GetXaxis()->SetBinLabel(6,"MV PileUp Rejection");
    fCutConfig->GetXaxis()->SetBinLabel(7,"SPD Mult");
    fCutConfig->GetXaxis()->SetBinLabel(8,"V0A Mult");
    fCutConfig->GetXaxis()->SetBinLabel(9,"V0C Mult");
    fCutConfig->GetXaxis()->SetBinLabel(10,"RefMult08");
    fCutConfig->GetXaxis()->SetBinLabel(11,"AliEvtCuts");
    fEventCutList->Add(fCutConfig);
    
    fV0Mpercentile=new TH1F("fV0Mpercentile","fV0Mpercentile",300,0,300);
    fV0Mpercentile->Sumw2();
    fV0Mpercentile->GetXaxis()->SetTitle("V0M percentile");
    fEventCutList->Add(fV0Mpercentile);
    
    fV0MpercentileHM=new TH1F("fV0MpercentileHM","fV0MpercentileHM",1000,0,10);
    fV0MpercentileHM->Sumw2();
    fV0MpercentileHM->GetXaxis()->SetTitle("V0M percentileHM");
    fEventCutList->Add(fV0MpercentileHM);
    
    TString sName[2]={"before","after"};
    
    for(int i=0;i<2;++i){
        
        fEvtCutQA[i]=new TList();
        fEvtCutQA[i]->SetName(sName[i].Data());
        fEvtCutQA[i]->SetOwner();
        
        fEventCutList->Add(fEvtCutQA[i]);
        
        TString nEvtNContName=Form("nContributors_%s",sName[i].Data());
        fEvtNCont[i]=new TH1F(nEvtNContName.Data(),nEvtNContName.Data(),350.,-0.5,349.5);
        fEvtNCont[i]->Sumw2();
        fEvtNCont[i]->GetXaxis()->SetTitle("Number of Contributors");
        fEvtCutQA[i]->Add(fEvtNCont[i]);
        
        TString EvtVtxXName=Form("VtxX_%s",sName[i].Data());
        fEvtVtxX[i]=new TH1F(EvtVtxXName.Data(),EvtVtxXName.Data(),50,-2.5 ,2.5);
        fEvtVtxX[i]->Sumw2();
        fEvtVtxX[i]->GetXaxis()->SetTitle("DCA_{x}");
        fEvtCutQA[i]->Add(fEvtVtxX[i]);
        
        TString EvtVtxYName=Form("VtxY_%s",sName[i].Data());
        fEvtVtxY[i]=new TH1F(EvtVtxYName.Data(),EvtVtxYName.Data(),50,-2.5 ,2.5);
        fEvtVtxY[i]->Sumw2();
        fEvtVtxY[i]->GetXaxis()->SetTitle("DCA_{y}");
        fEvtCutQA[i]->Add(fEvtVtxY[i]);
        
        TString EvtVtxZName=Form("VtxZ_%s",sName[i].Data());
        fEvtVtxZ[i]=new TH1F(EvtVtxZName.Data(),EvtVtxZName.Data(),300,-15. ,15.);
        fEvtVtxZ[i]->Sumw2();
        fEvtVtxZ[i]->GetXaxis()->SetTitle("DCA_{z}");
        fEvtCutQA[i]->Add(fEvtVtxZ[i]);
        
        TString MultNameSPD=Form("MultiplicitySPD_%s",sName[i].Data());
        fMultDistSPD[i]=new TH1F(MultNameSPD.Data(),MultNameSPD.Data(),600,0.,600.);
        fMultDistSPD[i]->Sumw2();
        fMultDistSPD[i]->GetXaxis()->SetTitle("Multiplicity (SPD)");
        fEvtCutQA[i]->Add(fMultDistSPD[i]);
        
        TString MultNameV0A=Form("MultiplicityV0A_%s",sName[i].Data());
        fMultDistV0A[i]=new TH1F(MultNameV0A.Data(),MultNameV0A.Data(),600,0.,600.);
        fMultDistV0A[i]->Sumw2();
        fMultDistV0A[i]->GetXaxis()->SetTitle("Multiplicity (V0A)");
        fEvtCutQA[i]->Add(fMultDistV0A[i]);
        
        TString MultNameV0C=Form("MultiplicityV0C_%s",sName[i].Data());
        fMultDistV0C[i]=new TH1F(MultNameV0C.Data(),MultNameV0C.Data(),600,0.,600.);
        fMultDistV0C[i]->Sumw2();
        fMultDistV0C[i]->GetXaxis()->SetTitle("Multiplicity (V0C)");
        fEvtCutQA[i]->Add(fMultDistV0C[i]);
        
        TString MultNameRefMult08=Form("MultiplicityRef08_%s",sName[i].Data());
        fMultDistRef08[i]=new TH1F(
                                   MultNameRefMult08.Data(),MultNameRefMult08.Data(),600,0.,600.);
        fMultDistRef08[i]->Sumw2();
        fMultDistRef08[i]->GetXaxis()->SetTitle("Multiplicity (RefMult08)");
        fEvtCutQA[i]->Add(fMultDistRef08[i]);
    }
}

AliLightNEventHist::~AliLightNEventHist() {
    
}

