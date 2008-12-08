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
 
 //choosing the style and color of mesh for MC error bands:
 Int_t meshStyle = 3044;//see documentation of TAttFill
 Int_t meshColor = kRed-4;
 
 //type of analysis was: ESD, AOD, MC, ESDMC0, ESDMC1
 const TString type = "ESD";
 
 //open the output files:
 TString inputFileNameMCEP = "outputMCEPanalysis";
 TFile* file_MCEP = NULL;
 file_MCEP = TFile::Open(((inputFileNameMCEP.Append(type)).Append(".root")).Data(), "READ"); 
 
 TString inputFileNameLYZ1 = "outputLYZ1analysis";
 TFile* file_LYZ1 = NULL;
 file_LYZ1 = TFile::Open(((inputFileNameLYZ1.Append(type)).Append("_firstrun.root")).Data(), "READ"); 
 
 /*
 TString inputFileNameSP = "outputSPanalysis";
 TFile* file_SP = NULL;
 file_SP = TFile::Open(((inputFileNameSP.Append(type)).Append(".root")).Data(), "READ"); 
 */
 
 TString inputFileNameFQD = "outputFQDanalysis";
 TFile* file_FQD = NULL;
 file_FQD = TFile::Open(((inputFileNameFQD.Append(type)).Append(".root")).Data(), "READ"); 
 
 TString inputFileNameGFC = "outputGFCanalysis";
 TFile* file_GFC = NULL;
 file_GFC = TFile::Open(((inputFileNameGFC.Append(type)).Append(".root")).Data(), "READ"); 
 
 TString inputFileNameQC = "outputQCanalysis";
 TFile* file_QC = NULL;
 file_QC = TFile::Open(((inputFileNameQC.Append(type)).Append(".root")).Data(), "READ"); 

 //==================================================================================
 //                              INTEGRATED FLOW
 //==================================================================================
 
 //booking the histogram for the integrated flow results from all methods
 TH1D* intFlow = new TH1D("intFlow","Integrated Flow",11,0,11);  

 intFlow->SetLabelSize(0.044);
 intFlow->SetMarkerStyle(21);
 intFlow->SetMarkerColor(kRed-4);
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
 AliFlowCommonHistResults *mcepCommonHistRes = NULL; 
 if(file_MCEP) 
 {
  file_MCEP->GetObject("cobjMCEP",pListMCEP); 
  if(pListMCEP) 
  {
   mcepCommonHist    = dynamic_cast<AliFlowCommonHist*> (pListMCEP->FindObject("AliFlowCommonHistMCEP"));
   mcepCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pListMCEP->FindObject("AliFlowCommonHistResultsMCEP"));
   if(mcepCommonHistRes)
   {
    intFlow->SetBinContent(1,(mcepCommonHistRes->GetHistIntFlow())->GetBinContent(1));   
    intFlow->SetBinError(1,(mcepCommonHistRes->GetHistIntFlow())->GetBinError(1));
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
 AliFlowCommonHist *qcCommonHist2 = NULL, *qcCommonHist4 = NULL;
 AliFlowCommonHistResults *qcCommonHistRes2 = NULL; 
 AliFlowCommonHistResults *qcCommonHistRes4 = NULL; 
 AliFlowCommonHistResults *qcCommonHistRes6 = NULL; 
 AliFlowCommonHistResults *qcCommonHistRes8 = NULL; 
 if(file_QC) 
 {
  file_QC->GetObject("cobjQC",pListQC);
  if(pListQC) 
  {
   qcCommonHist2 = dynamic_cast<AliFlowCommonHist*> (pListQC->FindObject("AliFlowCommonHistQC"));//to be improved
   qcCommonHistRes2 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults2ndOrderQC"));
   if(qcCommonHistRes2) 
   {
    intFlow->SetBinContent(8,(qcCommonHistRes2->GetHistIntFlow())->GetBinContent(1)); 
    intFlow->SetBinError(8,(qcCommonHistRes2->GetHistIntFlow())->GetBinError(1));
   }
   qcCommonHist4 = dynamic_cast<AliFlowCommonHist*> (pListQC->FindObject("AliFlowCommonHistQC"));//to be improved
   qcCommonHistRes4 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults4thOrderQC"));
   if(qcCommonHistRes4) 
   {
    intFlow->SetBinContent(9,(qcCommonHistRes4->GetHistIntFlow())->GetBinContent(1));
    intFlow->SetBinError(9,(qcCommonHistRes4->GetHistIntFlow())->GetBinError(1));
   }
   qcCommonHistRes6 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults6thOrderQC"));
   if(qcCommonHistRes6) 
   {
    intFlow->SetBinContent(10,(qcCommonHistRes6->GetHistIntFlow())->GetBinContent(1)); 
    intFlow->SetBinError(10,(qcCommonHistRes6->GetHistIntFlow())->GetBinError(1));
   }
   qcCommonHistRes8 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults8thOrderQC"));
   if(qcCommonHistRes8) 
   {
    intFlow->SetBinContent(11,(qcCommonHistRes8->GetHistIntFlow())->GetBinContent(1));
    intFlow->SetBinError(11,(qcCommonHistRes8->GetHistIntFlow())->GetBinError(1));
   }
  }
 }
 
 //----------------------------------------------------------------------------------
 //cosmetics: mesh for MC error bands (integrated flow)
 TGraph* pGraphIF = NULL;
 if(intFlow && mcepCommonHistRes)
 {
  Double_t valueMC = intFlow->GetBinContent(1);
  Double_t errorMC = intFlow->GetBinError(1);  
  Int_t nPts       = intFlow->GetNbinsX();     
     
  pGraphIF = new TGraph(nPts); 
  
  pGraphIF->SetPoint(1,0,valueMC+errorMC);
  pGraphIF->SetPoint(2,nPts+1,valueMC+errorMC);
  pGraphIF->SetPoint(3,nPts+1,valueMC-errorMC);
  pGraphIF->SetPoint(4,0,valueMC-errorMC);
  pGraphIF->SetPoint(5,0,valueMC+errorMC);
  
  pGraphIF->SetFillStyle(meshStyle);
  pGraphIF->SetFillColor(meshColor);
 }
  
 /*                     
 //cosmetics: legend     
 TString *avM = new TString("M = ");//to be improved
 TString *nEvts = new TString("N = ");//to be improved
 TLegend *legendIntFlow = new TLegend(0.15,0.15,0.44,0.35);
 legendIntFlow->SetTextFont(72);
 legendIntFlow->SetTextSize(0.03);
 if(qcCommonHist2)
 {
  (*avM)+=(qcCommonHist->GetHistMultInt())->GetMean();//to be improved
  (*N)+=(Long_t)(qcCommonHist->GetHistMultInt())->GetEntries();//to be improved
  legendIntFlow->AddEntry(qcCommonHist->GetHistMultInt(),avM->Data(),"");
  legendIntFlow->AddEntry(qcCommonHist->GetHistMultInt(),N->Data(),"");
 }
 */
 
 /*
 //cosmetics: legend (integrated flow)    
 TString *entryMC  = new TString("MC ......... ");
 TString *entryQC2 = new TString("QC{2} .... ");
 TString *entryQC4 = new TString("QC{4} .... ");
 
 TLegend *legendIntFlow = new TLegend(0.15,0.15,0.44,0.35); 
 legendIntFlow->SetTextFont(72);
 legendIntFlow->SetTextSize(0.03);
 
 
 Double_t avMultMC=0., avMultQC2=0., avMultQC4=0.;
 Long_t nEvtsMC=0, nEvtsQC2=0, nEvtsQC4=0;
 if(mcepCommonHist)
 {
  avMultMC = (mcepCommonHist->GetHistMultInt())->GetMean();
  nEvtsMC  = (mcepCommonHist->GetHistMultInt())->GetEntries();
  entryMC->Append("M = ");
  (*entryMC)+=(Long_t)avMultMC;
  entryMC->Append(", N = ");
  (*entryMC)+=(Long_t)nEvtsMC;
  legendIntFlow->AddEntry(mcepCommonHist->GetHistMultInt(),entryMC->Data(),"f");
 }
 if(qcCommonHist2)
 {
  avMultQC2 = (qcCommonHist2->GetHistMultInt())->GetMean();
  nEvtsQC2  = (qcCommonHist2->GetHistMultInt())->GetEntries();
  entryQC2->Append("M = ");
  (*entryQC2)+=(Long_t)avMultQC2;
  entryQC2->Append(", N = ");
  (*entryQC2)+=(Long_t)nEvtsQC2;
  legendIntFlow->AddEntry(qcCommonHist2->GetHistMultInt(),entryQC2->Data(),"f");
 }
 if(qcCommonHist4)
 {
  avMultQC4 = (qcCommonHist4->GetHistMultInt())->GetMean();
  nEvtsQC4  = (qcCommonHist4->GetHistMultInt())->GetEntries();
  entryQC4->Append("M = ");
  (*entryQC4)+=(Long_t)avMultQC4;
  entryQC4->Append(", N = ");
  (*entryQC4)+=(Long_t)nEvtsQC4;
  legendIntFlow->AddEntry(qcCommonHist4->GetHistMultInt(),entryQC4->Data(),"f");
 }
 */ 
 
 
 
 //cosmetics: text (integrated flow)    
 TString *entryMC  = new TString("MC ......... ");
 TString *entryQC2 = new TString("QC{2} .... ");
 TString *entryQC4 = new TString("QC{4} .... ");
 
 TPaveText *text = new TPaveText(0.15,0.15,0.44,0.35,"NDC");
 text->SetTextFont(72);
 text->SetTextSize(0.03);
 
 Double_t avMultMC=0., avMultQC2=0., avMultQC4=0.;
 Long_t nEvtsMC=0, nEvtsQC2=0, nEvtsQC4=0;
 if(mcepCommonHist)
 {
  avMultMC = (mcepCommonHist->GetHistMultInt())->GetMean();
  nEvtsMC  = (mcepCommonHist->GetHistMultInt())->GetEntries();
  entryMC->Append("M = ");
  (*entryMC)+=(Long_t)avMultMC;
  entryMC->Append(", N = ");
  (*entryMC)+=(Long_t)nEvtsMC;
  text->AddText(entryMC->Data());
 }
 if(qcCommonHist2)
 {
  avMultQC2 = (qcCommonHist2->GetHistMultInt())->GetMean();
  nEvtsQC2  = (qcCommonHist2->GetHistMultInt())->GetEntries();
  entryQC2->Append("M = ");
  (*entryQC2)+=(Long_t)avMultQC2;
  entryQC2->Append(", N = ");
  (*entryQC2)+=(Long_t)nEvtsQC2;
  text->AddText(entryQC2->Data());
 }
 if(qcCommonHist4)
 {
  avMultQC4 = (qcCommonHist4->GetHistMultInt())->GetMean();
  nEvtsQC4  = (qcCommonHist4->GetHistMultInt())->GetEntries();
  entryQC4->Append("M = ");
  (*entryQC4)+=(Long_t)avMultQC4;
  entryQC4->Append(", N = ");
  (*entryQC4)+=(Long_t)nEvtsQC4;
  text->AddText(entryQC4->Data());
 }
 
 //----------------------------------------------------------------------------------
 
 //drawing everything for integrated flow:
 TCanvas* intFlowCanvas = new TCanvas("Integrated Flow","Integrated Flow",1000,600);
 
 if(intFlow) intFlow->Draw("E1");      
 if(pGraphIF) pGraphIF->Draw("LFSAME");                    
 //if(legendIntFlow) legendIntFlow->Draw("");   
 text->Draw("SAME");
 
 //==================================================================================   
 



 //==================================================================================
 //                            DIFFERENTIAL FLOW
 //==================================================================================

 TCanvas* diffFlowCanvas = new TCanvas("Differential Flow","Differential Flow",1000,600);
 
 Int_t iNbinsPt  = AliFlowCommonConstants::GetNbinsPt();
 Double_t dPtMin = AliFlowCommonConstants::GetPtMin();
 Double_t dPtMax = AliFlowCommonConstants::GetPtMax();
 
 TH1D *styleHist = new TH1D("styleHist","styleHist",iNbinsPt,dPtMin,dPtMax);
 styleHist->SetTitle("Differential Flow");
 styleHist->SetXTitle("p_{t} [GeV]");
 styleHist->SetYTitle("v_{n}");
 styleHist->Draw();
 
 TString *entryDiffMC  = new TString("v_{n}{MC} ");
 TString *entryDiffQC2 = new TString("v_{n}^{Q}{2} ");
 TString *entryDiffQC4 = new TString("v_{n}^{Q}{4} ");
 
 TLegend* legendDiffFlow = new TLegend(0.15,0.15,0.44,0.35);
 legendDiffFlow->SetTextFont(72);
 legendDiffFlow->SetTextSize(0.03);
 
 //MCEP = Monte Carlo Event Plane
 if(file_MCEP)
 {
  if(mcepCommonHistRes)
  {
   (mcepCommonHistRes->GetHistDiffFlow())->Scale(0.01);//to be improved
   (mcepCommonHistRes->GetHistDiffFlow())->SetMarkerColor(2);
   (mcepCommonHistRes->GetHistDiffFlow())->SetMarkerStyle(20);
   (mcepCommonHistRes->GetHistDiffFlow())->SetFillStyle(meshStyle);
   (mcepCommonHistRes->GetHistDiffFlow())->SetFillColor(meshColor);
   entryDiffMC->Append("(M = ");
   (*entryDiffMC)+=(Long_t)avMultMC;
   entryDiffMC->Append(", N = ");
   (*entryDiffMC)+=(Long_t)nEvtsMC; 
   entryDiffMC->Append(")");
   //(mcepCommonHistRes->GetHistDiffFlow())->Draw("E1PSAME"); 
   legendDiffFlow->AddEntry(mcepCommonHistRes->GetHistDiffFlow(),entryDiffMC->Data(),"f");
  } 
 } 
 
 //GFC = Generating Function Cumulants
 if(file_GFC)
 {
  if(gfc2)
  {
   (gfc2->GetHistDiffFlow())->SetMarkerColor(28);
   (gfc2->GetHistDiffFlow())->SetMarkerStyle(20);
   (gfc2->GetHistDiffFlow())->Draw("E1PSAME"); 
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
  if(qcCommonHistRes2)
  {
   (qcCommonHistRes2->GetHistDiffFlow())->SetMarkerColor(44);
   (qcCommonHistRes2->GetHistDiffFlow())->SetMarkerStyle(20);
   (qcCommonHistRes2->GetHistDiffFlow())->Draw("E1PSAME");
   entryDiffQC2->Append("(M = ");
   (*entryDiffQC2)+=(Long_t)avMultQC2;
   entryDiffQC2->Append(", N = ");
   (*entryDiffQC2)+=(Long_t)nEvtsQC2; 
   entryDiffQC2->Append(")");
   legendDiffFlow->AddEntry(qcCommonHistRes2->GetHistDiffFlow(),entryDiffQC2->Data(),"p"); 
  }
  if(qcCommonHistRes4)
  {
   (qcCommonHistRes4->GetHistDiffFlow())->SetMarkerColor(44);
   (qcCommonHistRes4->GetHistDiffFlow())->SetMarkerStyle(21);
   (qcCommonHistRes4->GetHistDiffFlow())->Draw("E1PSAME"); 
   entryDiffQC4->Append("(M = ");
   (*entryDiffQC4)+=(Long_t)avMultQC4;
   entryDiffQC4->Append(", N = ");
   (*entryDiffQC4)+=(Long_t)nEvtsQC4; 
   entryDiffQC4->Append(")");
   legendDiffFlow->AddEntry(qcCommonHistRes4->GetHistDiffFlow(),entryDiffQC4->Data(),"p");
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
 
 
 
 
 
 
 
 
 

 //cosmetics: additional legend (differential flow) 
 
 
 /*
 
 //to be improved
 TLegend* legendDiffFlow2 = new TLegend(0.15,0.60,0.40,0.80);
 legendDiffFlow2->SetTextFont(72);
 legendDiffFlow2->SetTextSize(0.03);
 if(qcCommonHist2)
 {
  legendDiffFlow2->AddEntry(qcCommonHist2->GetHistMultInt(),avM->Data(),"");
  legendDiffFlow2->AddEntry(qcCommonHist2->GetHistMultInt(),nEvts->Data(),"");
 }
 if(qcCommonHist2)
 {
  legendDiffFlow2->AddEntry(qcCommonHist2->GetHistMultInt(),avM->Data(),"");
  legendDiffFlow2->AddEntry(qcCommonHist2->GetHistMultInt(),nEvts->Data(),"");
 }
 if(qcCommonHist2)
 {
  legendDiffFlow2->AddEntry(qcCommonHist2->GetHistMultInt(),avM->Data(),"");
  legendDiffFlow2->AddEntry(qcCommonHist2->GetHistMultInt(),nEvts->Data(),"");
 }
 
 */
 
 //cosmetics: Monte Carlo error bands for differential flow
 TGraph* pGraphDF = NULL;
 if(mcepCommonHistRes)
 {
  Int_t nPtsDF       = (mcepCommonHistRes->GetHistDiffFlow())->GetNbinsX();
  Double_t binWidth  = (mcepCommonHistRes->GetHistDiffFlow())->GetBinWidth(1);//assuming that all bins have the same width
       
  pGraphDF = new TGraph(2*nPts);
  
  for(Int_t i=1;i<nPtsDF+1;i++)
  {
   Double_t valueMC = (mcepCommonHistRes->GetHistDiffFlow())->GetBinContent(i);
   Double_t errorMC = (mcepCommonHistRes->GetHistDiffFlow())->GetBinError(i);       
   pGraphDF->SetPoint(i,(i-0.5)*binWidth,valueMC+errorMC);
  }
     
  for(Int_t i=nPtsDF;i<2*nPtsDF;i++)
  {
   Double_t valueMC = (mcepCommonHistRes->GetHistDiffFlow())->GetBinContent(2*nPtsDF-i);
   Double_t errorMC = (mcepCommonHistRes->GetHistDiffFlow())->GetBinError(2*nPtsDF-i);       
   pGraphDF->SetPoint(i,(2*nPtsDF-i-0.5)*binWidth,valueMC-errorMC); 
  }

  pGraphDF->SetFillStyle(meshStyle);
  pGraphDF->SetFillColor(meshColor);
 }            
                     
 if(legendDiffFlow) legendDiffFlow->Draw("");   
 if(pGraphDF) pGraphDF->Draw("LFSAME");        
                                                                                            
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
