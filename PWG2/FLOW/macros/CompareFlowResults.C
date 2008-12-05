void CompareFlowResults()
{
 gSystem->AddIncludePath("-I$ALICE_ROOT/include");
 gSystem->AddIncludePath("-I$ROOTSYS/include");
 
 //load needed libraries:
 gSystem->Load("libTree.so");
 gSystem->Load("libANALYSIS.so");
 gSystem->Load("libPWG2flow.so");
 cerr<<"libPWG2flow.so loaded ..."<<endl;
 cout<<endl;

 //removing the title and stat. box from all histograms:
 gStyle->SetOptTitle(0);
 gStyle->SetOptStat(0);
 
 //open the output files:
 TString inputFileNameMCEP = "outputMCEPanalysisESD.root";
 TFile* file_MCEP = NULL;
 file_MCEP = TFile::Open(inputFileNameMCEP.Data(), "READ"); 
 
 TString inputFileNameLYZ1 = "outputLYZ1analysisESD_firstrun.root";
 TFile* file_LYZ1 = NULL;
 file_LYZ1 = TFile::Open(inputFileNameLYZ1.Data(), "READ"); 
 
 /*
 TString inputFileNameSP = "outputSPanalysisESD.root";
 TFile* file_SP = NULL;
 file_SP = TFile::Open(inputFileNameSP.Data(), "READ"); 
 */
 
 TString inputFileNameFQD = "outputFQDanalysisESD.root";
 TFile* file_FQD = NULL;
 file_FQD = TFile::Open(inputFileNameFQD.Data(), "READ"); 
 
 TString inputFileNameGFC = "outputGFCanalysisESD.root";
 TFile* file_GFC = NULL;
 file_GFC = TFile::Open(inputFileNameGFC.Data(), "READ"); 
 
 TString inputFileNameQC = "outputQCanalysisESD.root";
 TFile* file_QC = NULL;
 file_QC = TFile::Open(inputFileNameQC.Data(), "READ"); 

 //==================================================================================
 //                              INTEGRATED FLOW
 //==================================================================================
 
 //booking the histogram for the integrated flow results from all methods
 TH1D* intFlow = new TH1D("intFlow","Integrated Flow",11,0,11);  

 intFlow->SetLabelSize(0.044);
 intFlow->SetMarkerStyle(25);
 (intFlow->GetXaxis())->SetBinLabel(1,"v_{n}{MC}");
 (intFlow->GetXaxis())->SetBinLabel(2,"v_{n}{LYZ}");
 (intFlow->GetXaxis())->SetBinLabel(3,"v_{n}{FQD}");
 (intFlow->GetXaxis())->SetBinLabel(4,"v_{n}{2}");
 (intFlow->GetXaxis())->SetBinLabel(5,"v_{n}{4}");
 (intFlow->GetXaxis())->SetBinLabel(6,"v_{n}{6}");
 (intFlow->GetXaxis())->SetBinLabel(7,"v_{n}{8}");
 (intFlow->GetXaxis())->SetBinLabel(8,"v_{n}^{Q}{2}");
 (intFlow->GetXaxis())->SetBinLabel(9,"v_{n}^{Q}{4}");
 (intFlow->GetXaxis())->SetBinLabel(10,"v_{n}^{Q}{6}");
 (intFlow->GetXaxis())->SetBinLabel(11,"v_{n}^{Q}{8}"); 
 
 //MCEP = Monte Carlo Event Plane
 TList *pListMCEP = NULL;
 AliFlowCommonHistResults *mcep = NULL; 
 if(file_MCEP) 
 {
  file_MCEP->GetObject("cobjMCEP",pListMCEP); 
  if(pListMCEP) 
  {
   mcep = dynamic_cast<AliFlowCommonHistResults*> (pListMCEP->FindObject("AliFlowCommonHistResultsMCEP"));
   if(mcep)
   {
    intFlow->SetBinContent(1,(mcep->GetHistIntFlow())->GetBinContent(1));   
    intFlow->SetBinError(1,(mcep->GetHistIntFlow())->GetBinError(1));
   }
  }
 }
 
 //LYZ1 = Lee-Yang Zeros (1st run)
 TList *pListLYZ1 = NULL;
 AliFlowCommonHistResults *lyz1 = NULL; 
 if(file_LYZ1) 
 {
  file_LYZ1->GetObject("cobjLYZ1",pListLYZ1); 
  if(pListLYZ1) 
  {
   lyz1 = dynamic_cast<AliFlowCommonHistResults*> (pListLYZ1->FindObject("AliFlowCommonHistResultsLYZ"));
   if(lyz1)
   {
    intFlow->SetBinContent(2,(lyz1->GetHistIntFlow())->GetBinContent(1));   
    intFlow->SetBinError(2,(lyz1->GetHistIntFlow())->GetBinError(1));
   }
  }
 }
 
 //FQD = Fitting q-distribution
 TList *pListFQD = NULL;
 AliFlowCommonHistResults *fqd = NULL; 
 if(file_FQD) 
 {
  file_FQD->GetObject("cobjFQD",pListFQD); 
  if(pListFQD) 
  {
   fqd = dynamic_cast<AliFlowCommonHistResults*> (pListFQD->FindObject("AliFlowCommonHistResultsFQD"));
   if(fqd)
   {
    intFlow->SetBinContent(3,(fqd->GetHistIntFlow())->GetBinContent(1));   
    intFlow->SetBinError(3,(fqd->GetHistIntFlow())->GetBinError(1));
   }
  }
 }
 
 //GFC = Generating Function Cumulants
 TList *pListGFC = NULL;
 AliFlowCommonHistResults *gfc2 = NULL; 
 AliFlowCommonHistResults *gfc4 = NULL; 
 AliFlowCommonHistResults *gfc6 = NULL; 
 AliFlowCommonHistResults *gfc8 = NULL; 
 if(file_GFC) 
 {
  file_GFC->GetObject("cobjGFC",pListGFC);
  if(pListGFC) 
  {
   gfc2 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults2ndOrderGFC"));
   if(gfc2) 
   {
    intFlow->SetBinContent(4,(gfc2->GetHistIntFlow())->GetBinContent(1)); 
    intFlow->SetBinError(4,(gfc2->GetHistIntFlow())->GetBinError(1));
   }
   gfc4 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults4thOrderGFC"));
   if(gfc4) 
   {
    intFlow->SetBinContent(5,(gfc4->GetHistIntFlow())->GetBinContent(1));
    intFlow->SetBinError(5,(gfc4->GetHistIntFlow())->GetBinError(1));
   }
   gfc6 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults6thOrderGFC"));
   if(gfc6) 
   {
    intFlow->SetBinContent(6,(gfc6->GetHistIntFlow())->GetBinContent(1)); 
    intFlow->SetBinError(6,(gfc6->GetHistIntFlow())->GetBinError(1));
   }
   gfc8 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults8thOrderGFC"));
   if(gfc8) 
   {
    intFlow->SetBinContent(7,(gfc8->GetHistIntFlow())->GetBinContent(1));
    intFlow->SetBinError(7,(gfc8->GetHistIntFlow())->GetBinError(1));
   }
  }
 }
 
 //QC = Q-cumulants
 TList *pListQC = NULL;
 AliFlowCommonHist *qcCommonHist = NULL;//to be improved
 AliFlowCommonHistResults *qc2 = NULL; 
 AliFlowCommonHistResults *qc4 = NULL; 
 AliFlowCommonHistResults *qc6 = NULL; 
 AliFlowCommonHistResults *qc8 = NULL; 
 if(file_QC) 
 {
  file_QC->GetObject("cobjQC",pListQC);
  if(pListQC) 
  {
   qcCommonHist = dynamic_cast<AliFlowCommonHist*> (pListQC->FindObject("AliFlowCommonHistQC"));//to be improved
   qc2 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults2ndOrderQC"));
   if(qc2) 
   {
    intFlow->SetBinContent(8,(qc2->GetHistIntFlow())->GetBinContent(1)); 
    intFlow->SetBinError(8,(qc2->GetHistIntFlow())->GetBinError(1));
   }
   qc4 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults4thOrderQC"));
   if(qc4) 
   {
    intFlow->SetBinContent(9,(qc4->GetHistIntFlow())->GetBinContent(1));
    intFlow->SetBinError(9,(qc4->GetHistIntFlow())->GetBinError(1));
   }
   qc6 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults6thOrderQC"));
   if(qc6) 
   {
    intFlow->SetBinContent(10,(qc6->GetHistIntFlow())->GetBinContent(1)); 
    intFlow->SetBinError(10,(qc6->GetHistIntFlow())->GetBinError(1));
   }
   qc8 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults8thOrderQC"));
   if(qc8) 
   {
    intFlow->SetBinContent(11,(qc8->GetHistIntFlow())->GetBinContent(1));
    intFlow->SetBinError(11,(qc8->GetHistIntFlow())->GetBinError(1));
   }
  }
 }
 
 /*
 //SP = Scalar Product
 TList *pListSP = NULL;
 AliFlowCommonHistResults *sp = NULL; 
 if(file_SP) 
 {
  file_SP->GetObject("cobjSP",pListSP);
  if(pListSP) 
  {
   sp = dynamic_cast<AliFlowCommonHistResults*> (pListSP->FindObject("AliFlowCommonHistResultsSP"));
  }
 }
 */
 
 //cosmetics: Monte Carlo error bands for integrated flow
 if(intFlow)
 {
  Double_t valueMC = intFlow->GetBinContent(1);
  Double_t errorMC = intFlow->GetBinError(1);  
  Int_t nPts       = intFlow->GetNbinsX();     
     
  TGraph* pGraphIF = new TGraph(nPts); 
  
  pGraphIF->SetPoint(1,0,valueMC+errorMC);
  pGraphIF->SetPoint(2,nPts+1,valueMC+errorMC);
  pGraphIF->SetPoint(3,nPts+1,valueMC-errorMC);
  pGraphIF->SetPoint(4,0,valueMC-errorMC);
  
  pGraphIF->SetFillStyle(3044);
  pGraphIF->SetFillColor(kRed-4);
 }
            
 //cosmetics: legend     
 TString *avM = new TString("AvM = ");//to be improved
 TString *nEvts = new TString("nEvts = ");//to be improved
 TLegend* legendIntFlow = new TLegend(0.15,0.15,0.44,0.35);
 legendIntFlow->SetTextFont(72);
 legendIntFlow->SetTextSize(0.03);
 if(qcCommonHist)
 {
  (*avM)+=(qcCommonHist->GetHistMultInt())->GetMean();//to be improved
  (*nEvts)+=(qcCommonHist->GetHistMultInt())->GetEntries();//to be improved
  legendIntFlow->AddEntry(qcCommonHist->GetHistMultInt(),avM->Data(),"");
  legendIntFlow->AddEntry(qcCommonHist->GetHistMultInt(),nEvts->Data(),"");
 }
 
 //drawing everything for integrated flow:
 TCanvas* intFlowCanvas = new TCanvas("Integrated Flow","Integrated Flow",1000,600);
 
 intFlow->Draw();      
 pGraphIF->Draw("LFSAME");                    
 legendIntFlow->Draw("");  
 
 //==================================================================================   
 



 //==================================================================================
 //                            DIFFERENTIAL FLOW
 //==================================================================================

 TCanvas* diffFlowCanvas = new TCanvas("Differential Flow","Differential Flow",1000,600);
 
 TLegend* legendDiffFlow = new TLegend(0.15,0.15,0.3,0.35);
 legendDiffFlow->SetTextFont(72);
 legendDiffFlow->SetTextSize(0.03);
 
 //to be improved
 TLegend* legendDiffFlow2 = new TLegend(0.15,0.60,0.40,0.80);
 legendDiffFlow2->SetTextFont(72);
 legendDiffFlow2->SetTextSize(0.03);
 if(qcCommonHist)
 {
  legendDiffFlow2->AddEntry(qcCommonHist->GetHistMultInt(),avM->Data(),"");
  legendDiffFlow2->AddEntry(qcCommonHist->GetHistMultInt(),nEvts->Data(),"");
 }
 
 //GFC = Generating Function Cumulants
 if(file_GFC)
 {
  if(gfc2)
  {
   (gfc2->GetHistDiffFlow())->SetMarkerColor(28);
   (gfc2->GetHistDiffFlow())->SetMarkerStyle(20);
   (gfc2->GetHistDiffFlow())->SetTitle("Differential Flow");
   (gfc2->GetHistDiffFlow())->SetXTitle("p_{t} [GeV]");
   (gfc2->GetHistDiffFlow())->SetYTitle("v_{n}");
   (gfc2->GetHistDiffFlow())->Draw(""); 
   legendDiffFlow->AddEntry(gfc2->GetHistDiffFlow(),"v_{n}{2}","p");
  }
  if(gfc4)
  { 
   (gfc4->GetHistDiffFlow())->SetMarkerColor(28);
   (gfc4->GetHistDiffFlow())->SetMarkerStyle(21);
   (gfc4->GetHistDiffFlow())->Draw("E1PSAME");
   legendDiffFlow->AddEntry(gfc4->GetHistDiffFlow(),"v_{n}{4}","p"); 
   //(gfc6->GetHistDiffFlow())->Draw("SAME"); 
   //(gfc8->GetHistDiffFlow())->Draw("SAME"); 
  }
 }

 //QC = Q-cumulants
 if(file_QC)
 {
  if(qc2)
  {
   (qc2->GetHistDiffFlow())->SetMarkerColor(44);
   (qc2->GetHistDiffFlow())->SetMarkerStyle(20);
   (qc2->GetHistDiffFlow())->Draw("E1PSAME");
   legendDiffFlow->AddEntry(qc2->GetHistDiffFlow(),"v_{n}^{Q}{2}","p"); 
  }
  if(qc4)
  {
   (qc4->GetHistDiffFlow())->SetMarkerColor(44);
   (qc4->GetHistDiffFlow())->SetMarkerStyle(21);
   (qc4->GetHistDiffFlow())->Draw("E1PSAME"); 
   legendDiffFlow->AddEntry(qc4->GetHistDiffFlow(),"v_{n}^{Q}{4}","p");
  }
 } 
 
 //MCEP = Monte Carlo Event Plane
 if(file_MCEP)
 {
  if(mcep)
  {
   (mcep->GetHistDiffFlow())->Scale(0.01);//to be improved
   (mcep->GetHistDiffFlow())->SetMarkerColor(2);
   (mcep->GetHistDiffFlow())->SetMarkerStyle(20);
   (mcep->GetHistDiffFlow())->SetFillStyle(3044);
   (mcep->GetHistDiffFlow())->SetFillColor(kRed-4);
   //(mcep->GetHistDiffFlow())->Draw("E1PSAME"); 
   legendDiffFlow->AddEntry(mcep->GetHistDiffFlow(),"v_{n}{MC}","f");
  } 
 } 

 /*
 //SP = Scalar Product
 if(file_SP)
 {
  if(sp)
  {
   (sp->GetHistDiffFlow())->SetMarkerColor(3);
   (sp->GetHistDiffFlow())->SetMarkerStyle(28);
   (sp->GetHistDiffFlow())->Draw("E1PSAME"); 
   legend1->AddEntry(sp->GetHistDiffFlow(),"v_{n}{SP}","p");
  } 
 } 
 */
 
 //cosmetics: Monte Carlo error bands for differential flow
 if(mcep)
 {
  Int_t nPtsDF       = (mcep->GetHistDiffFlow())->GetNbinsX();
  Double_t binWidth  = (mcep->GetHistDiffFlow())->GetBinWidth(1);//assuming that all bins have the same width
       
  TGraph* pGraphDF = new TGraph(2*nPts);
  
  for(Int_t i=1;i<nPtsDF+1;i++)
  {
   Double_t valueMC = (mcep->GetHistDiffFlow())->GetBinContent(i);
   Double_t errorMC = (mcep->GetHistDiffFlow())->GetBinError(i);       
   pGraphDF->SetPoint(i,(i-0.5)*binWidth,valueMC+errorMC);
  }
     
  for(Int_t i=nPtsDF;i<2*nPtsDF;i++)
  {
   Double_t valueMC = (mcep->GetHistDiffFlow())->GetBinContent(2*nPtsDF-i);
   Double_t errorMC = (mcep->GetHistDiffFlow())->GetBinError(2*nPtsDF-i);       
   pGraphDF->SetPoint(i,(2*nPtsDF-i-0.5)*binWidth,valueMC-errorMC); 
  }

  pGraphDF->SetFillStyle(3044);
  pGraphDF->SetFillColor(kRed-4);
 }            
                     
 legendDiffFlow->Draw("");   
 legendDiffFlow2->Draw(""); 
 pGraphDF->Draw("LFSAME");        
                                                                                            
 //================================================================================== 
 
 
 
 //==========================================================
/*
 TGraph* fillAreaBetweenFunctions(TF1* fUp, TF1*fDown, double nfmin, 
				 double nfmax, int npf, int fillStyle, 
				 int fillColor){
  
  fDown->SetLineColor(fillColor);
  fDown->SetLineWidth(0.3);
  fUp->SetLineColor(fillColor);
  fUp->SetLineWidth(0.3);
  
  TArrayD xf(2*npf+1), yf(2*npf+1);
  
  //make and closed area and fill it
  Double_t xfmin = nfmin; Double_t xfmax = nfmax;
  Double_t dxf = (xfmax-xfmin)/(npf-1);
  for (Int_t i=0;i<npf;i++) {
    xf[i] = xfmin + dxf*i;
    yf[i] = fDown->Eval(xf[i]);
    xf[npf+i] = xfmax - dxf*i;
    yf[npf+i] = fUp->Eval(xf[npf+i]);
  }
  xf[2*npf] = xf[0]; yf[2*npf] = yf[0];
  TGraph *grf = new TGraph(2*npf+1);
  for (int i=0; i<2*npf+1; i++) grf->SetPoint(i,xf[i],yf[i]);
  grf->SetFillStyle(fillStyle);
  grf->SetFillColor(fillColor);
  return grf;
  
  
}
*/
//=====================================================================================
 
 
 
 
				 
  
}