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
  

 //==================================================================================
 //                         accessing output files
 //==================================================================================
 //type of analysis was: ESD, AOD, MC, ESDMC0, ESDMC1
 const TString type = "ESD";
 
 //open the output files:
 TString inputFileNameMCEP = "outputMCEPanalysis";
 TFile* fileMCEP = NULL;
 fileMCEP = TFile::Open(((inputFileNameMCEP.Append(type)).Append(".root")).Data(), "READ"); 
 
 TString inputFileNameLYZ1 = "outputLYZ1analysis";
 TFile* fileLYZ1 = NULL;
 fileLYZ1 = TFile::Open(((inputFileNameLYZ1.Append(type)).Append("_firstrun.root")).Data(), "READ"); 

 TString inputFileNameLYZ2 = "outputLYZ2analysis";
 TFile* fileLYZ2 = NULL;
 fileLYZ2 = TFile::Open(((inputFileNameLYZ2.Append(type)).Append("_secondrun.root")).Data(), "READ"); 
 
 TString inputFileNameLYZEP = "outputLYZEPanalysis";
 TFile* fileLYZEP = NULL;
 fileLYZEP = TFile::Open(((inputFileNameLYZEP.Append(type)).Append(".root")).Data(), "READ");
 
 TString inputFileNameFQD = "outputFQDanalysis";
 TFile* fileFQD = NULL;
 fileFQD = TFile::Open(((inputFileNameFQD.Append(type)).Append(".root")).Data(), "READ"); 
 
 TString inputFileNameGFC = "outputGFCanalysis";
 TFile* fileGFC = NULL;
 fileGFC = TFile::Open(((inputFileNameGFC.Append(type)).Append(".root")).Data(), "READ"); 
 
 TString inputFileNameQC = "outputQCanalysis";
 TFile* fileQC = NULL;
 fileQC = TFile::Open(((inputFileNameQC.Append(type)).Append(".root")).Data(), "READ"); 
 //==================================================================================
 
 //==================================================================================
 //                                 cosmetics
 //==================================================================================
 //removing the title and stat. box from all histograms:
 gStyle->SetOptTitle(0);
 gStyle->SetOptStat(0);
 
 //choosing the style and color of mesh for MC error bands:
 Int_t meshStyle = 1001;
 Int_t meshColor = kRed-10;

 //marker style and color (int. flow) 
 Int_t markerStyle = 21;
 Int_t markerColor = kRed-3;
 //==================================================================================

   
      
         
               
 //==================================================================================
 //                              INTEGRATED FLOW
 //==================================================================================
 //the number of different methods:
 const Int_t nMethods=12;

 //booking the histogram for the integrated flow results from all methods:
 TH1D* intFlowAll = new TH1D("intFlowAll","Integrated Flow",nMethods,0,nMethods);      
 //intFlowAll->SetLabelSize(0.036,"X");
 //intFlowAll->SetLabelSize(0.036,"Y");
 intFlowAll->SetMarkerStyle(markerStyle);
 intFlowAll->SetMarkerColor(markerColor);
 (intFlowAll->GetXaxis())->SetBinLabel(1,"v_{2}{MC}");
 (intFlowAll->GetXaxis())->SetBinLabel(2,"v_{2}{2,GFC}");
 (intFlowAll->GetXaxis())->SetBinLabel(3,"v_{2}{2,QC}");
 (intFlowAll->GetXaxis())->SetBinLabel(4,"v_{2}{4,GFC}");
 (intFlowAll->GetXaxis())->SetBinLabel(5,"v_{2}{4,QC}");
 (intFlowAll->GetXaxis())->SetBinLabel(6,"v_{2}{6,GFC}");
 (intFlowAll->GetXaxis())->SetBinLabel(7,"v_{2}{6,QC}");
 (intFlowAll->GetXaxis())->SetBinLabel(8,"v_{2}{8,GFC}");
 (intFlowAll->GetXaxis())->SetBinLabel(9,"v_{2}{8,QC}");
 (intFlowAll->GetXaxis())->SetBinLabel(10,"v_{2}{FQD}");
 (intFlowAll->GetXaxis())->SetBinLabel(11,"v_{2}{LYZ}");
 (intFlowAll->GetXaxis())->SetBinLabel(12,"v_{2}{LYZEP}");
 
 //booking the graph to store flow values and errors from all methods:  
 Double_t x[nMethods] = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5};
 Double_t xError[nMethods] = {0.};
 Double_t flowValue[nMethods] = {0.};
 Double_t flowError[nMethods] = {0.};
 
 //accessing the results for each method:
 //MCEP = Monte Carlo Event Plane
 TList *pListMCEP = NULL;
 AliFlowCommonHist *mcepCommonHist = NULL;
 AliFlowCommonHistResults *mcepCommonHistRes = NULL; 
 if(fileMCEP) 
 {
  fileMCEP->GetObject("cobjMCEP",pListMCEP); 
  if(pListMCEP) 
  {
   mcepCommonHist    = dynamic_cast<AliFlowCommonHist*> (pListMCEP->FindObject("AliFlowCommonHistMCEP"));
   mcepCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pListMCEP->FindObject("AliFlowCommonHistResultsMCEP"));
   if(mcepCommonHistRes)
   {
    flowValue[0] = (mcepCommonHistRes->GetHistIntFlow())->GetBinContent(1);
    flowError[0] = (mcepCommonHistRes->GetHistIntFlow())->GetBinError(1);
   }
  }
 }
 
 //LYZ1 = Lee-Yang Zeros (1st run)
 TList *pListLYZ1 = NULL;
 AliFlowCommonHist *lyz1CommonHist = NULL;
 AliFlowCommonHistResults *lyz1CommonHistRes = NULL; 
 if(fileLYZ1) 
 {
  fileLYZ1->GetObject("cobjLYZ1",pListLYZ1); 
  if(pListLYZ1) 
  {
   lyz1CommonHist = dynamic_cast<AliFlowCommonHist*> (pListLYZ1->FindObject("AliFlowCommonHistLYZ1"));
   lyz1CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pListLYZ1->FindObject("AliFlowCommonHistResultsLYZ1"));
   if(lyz1CommonHistRes)
   {
    flowValue[10] = (lyz1CommonHistRes->GetHistIntFlow())->GetBinContent(1);
    flowError[10] = (lyz1CommonHistRes->GetHistIntFlow())->GetBinError(1);
   }
  }
 }
 
 //LYZ2 = Lee-Yang Zeros (2nd run) (needed only for differential flow)
 TList *pListLYZ2 = NULL;
 AliFlowCommonHist *lyz2CommonHist = NULL;
 AliFlowCommonHistResults *lyz2CommonHistRes = NULL; 
 if(fileLYZ2) 
 {
  fileLYZ2->GetObject("cobjLYZ2",pListLYZ2); 
  if(pListLYZ2) 
  {
   lyz2CommonHist = dynamic_cast<AliFlowCommonHist*> (pListLYZ2->FindObject("AliFlowCommonHistLYZ2"));
   lyz2CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pListLYZ2->FindObject("AliFlowCommonHistResultsLYZ2"));
  }
 }

 //LYZEP = Lee-Yang Zeros Event Plane
 TList *pListLYZEP = NULL;
 AliFlowCommonHist *lyzepCommonHist = NULL;
 AliFlowCommonHistResults *lyzepCommonHistRes = NULL; 
 if(fileLYZEP) 
 {
  fileLYZEP->GetObject("cobjLYZEP",pListLYZEP); 
  if(pListLYZEP) 
  {
   lyzepCommonHist = dynamic_cast<AliFlowCommonHist*> (pListLYZEP->FindObject("AliFlowCommonHistLYZEP"));
   lyzepCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pListLYZEP->FindObject("AliFlowCommonHistResultsLYZEP"));
   if(lyzepCommonHistRes)
   {
    flowValue[11] = (lyzepCommonHistRes->GetHistIntFlow())->GetBinContent(1);
    //flowError[11] = (lyzepCommonHistRes->GetHistIntFlow())->GetBinError(1);
   }
  }
 }
 
 //FQD = Fitting q-distribution
 TList *pListFQD = NULL;
 AliFlowCommonHist *fqdCommonHist = NULL;
 AliFlowCommonHistResults *fqdCommonHistRes = NULL; 
 if(fileFQD) 
 {
  fileFQD->GetObject("cobjFQD",pListFQD); 
  if(pListFQD) 
  {
   fqdCommonHist = dynamic_cast<AliFlowCommonHist*> (pListFQD->FindObject("AliFlowCommonHistFQD"));
   fqdCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pListFQD->FindObject("AliFlowCommonHistResultsFQD"));
   if(fqdCommonHistRes)
   {
    flowValue[9] = (fqdCommonHistRes->GetHistIntFlow())->GetBinContent(1);
    flowError[9] = (fqdCommonHistRes->GetHistIntFlow())->GetBinError(1);
   }
  }
 }
 
 //GFC = Generating Function Cumulants
 TList *pListGFC = NULL;
 AliFlowCommonHist *gfcCommonHist = NULL;
 AliFlowCommonHistResults *gfcCommonHistRes2 = NULL; 
 AliFlowCommonHistResults *gfcCommonHistRes4 = NULL; 
 AliFlowCommonHistResults *gfcCommonHistRes6 = NULL; 
 AliFlowCommonHistResults *gfcCommonHistRes8 = NULL; 
 if(fileGFC) 
 {
  fileGFC->GetObject("cobjGFC",pListGFC);
  if(pListGFC) 
  {
   gfcCommonHist = dynamic_cast<AliFlowCommonHist*> (pListGFC->FindObject("AliFlowCommonHistGFC"));
   gfcCommonHistRes2 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults2ndOrderGFC"));
   if(gfcCommonHistRes2) 
   {
    flowValue[1] = (gfcCommonHistRes2->GetHistIntFlow())->GetBinContent(1);
    flowError[1] = (gfcCommonHistRes2->GetHistIntFlow())->GetBinError(1);
   }
   gfcCommonHistRes4 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults4thOrderGFC"));
   if(gfcCommonHistRes4) 
   {
    flowValue[3] = (gfcCommonHistRes4->GetHistIntFlow())->GetBinContent(1);
    flowError[3] = (gfcCommonHistRes4->GetHistIntFlow())->GetBinError(1);
   }
   gfcCommonHistRes6 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults6thOrderGFC"));
   if(gfcCommonHistRes6) 
   {
    flowValue[5] = (gfcCommonHistRes6->GetHistIntFlow())->GetBinContent(1);
    flowError[5] = (gfcCommonHistRes6->GetHistIntFlow())->GetBinError(1);
   }
   gfcCommonHistRes8 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults8thOrderGFC"));
   if(gfcCommonHistRes8) 
   {
    flowValue[7] = (gfcCommonHistRes8->GetHistIntFlow())->GetBinContent(1);
    flowError[7] = (gfcCommonHistRes8->GetHistIntFlow())->GetBinError(1);
   }
  }
 }
 
 //QC = Q-cumulants
 TList *pListQC = NULL;
 AliFlowCommonHist *qcCommonHist2 = NULL; 
 AliFlowCommonHist *qcCommonHist4 = NULL; 
 AliFlowCommonHist *qcCommonHist6 = NULL; 
 AliFlowCommonHist *qcCommonHist8 = NULL; 
 AliFlowCommonHistResults *qcCommonHistRes2 = NULL; 
 AliFlowCommonHistResults *qcCommonHistRes4 = NULL; 
 AliFlowCommonHistResults *qcCommonHistRes6 = NULL; 
 AliFlowCommonHistResults *qcCommonHistRes8 = NULL; 
 
 if(fileQC) 
 {
  fileQC->GetObject("cobjQC",pListQC);
  if(pListQC) 
  {
   qcCommonHist2 = dynamic_cast<AliFlowCommonHist*> (pListQC->FindObject("AliFlowCommonHist2ndOrderQC"));
   qcCommonHistRes2 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults2ndOrderQC"));
   if(qcCommonHistRes2) 
   {
    flowValue[2] = (qcCommonHistRes2->GetHistIntFlow())->GetBinContent(1);
    //flowError[2] = (qcCommonHistRes2->GetHistIntFlow())->GetBinError(1);
   }
   qcCommonHist4 = dynamic_cast<AliFlowCommonHist*> (pListQC->FindObject("AliFlowCommonHist4thOrderQC"));
   qcCommonHistRes4 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults4thOrderQC"));
   if(qcCommonHistRes4) 
   {
    flowValue[4] = (qcCommonHistRes4->GetHistIntFlow())->GetBinContent(1);
    //flowError[4] = (qcCommonHistRes4->GetHistIntFlow())->GetBinError(1);
   }
   qcCommonHist6 = dynamic_cast<AliFlowCommonHist*> (pListQC->FindObject("AliFlowCommonHist6thOrderQC"));
   qcCommonHistRes6 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults6thOrderQC"));
   if(qcCommonHistRes6) 
   {
    flowValue[6] = (qcCommonHistRes6->GetHistIntFlow())->GetBinContent(1);
    //flowError[6] = (qcCommonHistRes6->GetHistIntFlow())->GetBinError(1);
   }
   qcCommonHist8 = dynamic_cast<AliFlowCommonHist*> (pListQC->FindObject("AliFlowCommonHist8thOrderQC"));
   qcCommonHistRes8 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults8thOrderQC"));
   if(qcCommonHistRes8) 
   {
    flowValue[8] = (qcCommonHistRes8->GetHistIntFlow())->GetBinContent(1);
    //flowError[8] = (qcCommonHistRes8->GetHistIntFlow())->GetBinError(1);
   }
  }
 }        

 Double_t dMax=flowValue[0]+flowError[0];
 Double_t dMin=flowValue[0]-flowError[0];

 for(Int_t i=1;i<nMethods;i++)
 {
  if(!(flowValue[i]==0. && flowError[i]==0.)) 
  {
   if(dMax<flowValue[i]+flowError[i]) dMax=flowValue[i]+flowError[i];
   if(dMin>flowValue[i]-flowError[i]) dMin=flowValue[i]-flowError[i];
  } 
 }  
                                                                                                                                                                                    
 TGraph* flowResults = new TGraphErrors(nMethods, x, flowValue, xError, flowError);
 
 flowResults->SetMarkerStyle(markerStyle);
 flowResults->SetMarkerColor(markerColor);
 //-----------------------------------------------------------------------------------
  
 //----------------------------------------------------------------------------------
 //cosmetics: mesh for MC error bands (integrated flow)
 TGraph* pMesh = NULL;
 if(intFlowAll && mcepCommonHistRes)
 {
  //Double_t valueMC = intFlowAll->GetBinContent(1);
  //Double_t errorMC = intFlowAll->GetBinError(1);  
  //Int_t nPts       = intFlowAll->GetNbinsX();     

  Double_t valueMC = flowValue[0];
  Double_t errorMC = flowError[0];  
  Int_t nPts       = nMethods;     
       
  pMesh = new TGraph(nPts); 
  
  pMesh->SetPoint(1,0,valueMC+errorMC);
  pMesh->SetPoint(2,nPts+1,valueMC+errorMC);
  pMesh->SetPoint(3,nPts+1,valueMC-errorMC);
  pMesh->SetPoint(4,0,valueMC-errorMC);
  pMesh->SetPoint(5,0,valueMC+errorMC);
  
  pMesh->SetFillStyle(meshStyle);
  pMesh->SetFillColor(meshColor);
 }
 //---------------------------------------------------------------------------------- 
 
 
 //----------------------------------------------------------------------------------
 //cosmetics: text (integrated flow) 
 //default text:
 TPaveText *textDefault = new TPaveText(0.05,0.77,0.95,0.90,"NDC");
 textDefault->SetTextFont(72);
 textDefault->SetTextSize(0.08);
 //textDefault->SetLineColor(kFALSE);
 //textDefault->SetShadowColor(kFALSE);

 TString *entryDefaultAvM = new TString("Average Multiplicity");
 TString *entryDefaultAnd = new TString("and"); 
 TString *entryDefaultNumOfEvts = new TString("Number of Events:");

 textDefault->AddText(entryDefaultAvM->Data());
 textDefault->AddText(entryDefaultAnd->Data());
 textDefault->AddText(entryDefaultNumOfEvts->Data());
 
 //results:
 TPaveText *textResults = new TPaveText(0.05,0.12,0.95,0.70,"NDC");
 textResults->SetTextFont(72);
 textResults->SetTextSize(0.06);
 //textResults->SetLineColor(kFALSE);
 //textResults->SetShadowColor(kFALSE);
       
 TString *entryMC    = new TString("MC ...... ");
 TString *entryGFC   = new TString("GFC ..... "); 
 TString *entryQC2   = new TString("QC{2} ... ");
 TString *entryQC4   = new TString("QC{4} ... ");
 TString *entryQC6   = new TString("QC{6} ... ");
 TString *entryQC8   = new TString("QC{8} ... ");
 TString *entryFQD   = new TString("FQD ..... "); 
 TString *entryLYZ1  = new TString("LYZ ..... "); 
 TString *entryLYZEP = new TString("LYZEP ... "); 

 Double_t avMultMC=0.;
 Long_t nEvtsMC=0;

 Double_t avMultGFC=0.;
 Long_t nEvtsGFC=0;
  
 Double_t avMultQC2=0., avMultQC4=0., avMultQC6=0., avMultQC8=0.;
 Long_t nEvtsQC2=0, nEvtsQC4=0, nEvtsQC6=0, nEvtsQC8=0;

 Double_t avMultFQD=0.;
 Long_t nEvtsFQD=0;

 Double_t avMultLYZ1=0.;
 Long_t nEvtsLYZ1=0;
 
 Double_t avMultLYZEP=0.;
 Long_t nEvtsLYZEP=0;
 
 //MC:  
 if(mcepCommonHist)
 {
  avMultMC = (mcepCommonHist->GetHistMultInt())->GetMean();
  nEvtsMC  = (mcepCommonHist->GetHistMultInt())->GetEntries();
 }
 if(entryMC)
 {   
  entryMC->Append("M = ");
  (*entryMC)+=(Long_t)avMultMC;
  entryMC->Append(", N = ");
  (*entryMC)+=(Long_t)nEvtsMC;
 }
 
 //GFC:
 if(gfcCommonHist)
 {
  avMultGFC = (gfcCommonHist->GetHistMultInt())->GetMean();
  nEvtsGFC  = (gfcCommonHist->GetHistMultInt())->GetEntries();
 }
 if(entryGFC)
 { 
  entryGFC->Append("M = ");
  (*entryGFC)+=(Long_t)avMultGFC;
  entryGFC->Append(", N = ");
  (*entryGFC)+=(Long_t)nEvtsGFC;
 }
 
 //QC:
 if(qcCommonHist2)
 {
  avMultQC2 = (qcCommonHist2->GetHistMultInt())->GetMean();
  nEvtsQC2  = (qcCommonHist2->GetHistMultInt())->GetEntries();
 } 
 if(entryQC2)
 { 
  entryQC2->Append("M = ");
  (*entryQC2)+=(Long_t)avMultQC2;
  entryQC2->Append(", N = ");
  (*entryQC2)+=(Long_t)nEvtsQC2;
 }
 if(qcCommonHist4)
 {
  avMultQC4 = (qcCommonHist4->GetHistMultInt())->GetMean();
  nEvtsQC4  = (qcCommonHist4->GetHistMultInt())->GetEntries();
 }
 if(entryQC4)
 {
  entryQC4->Append("M = ");
  (*entryQC4)+=(Long_t)avMultQC4;
  entryQC4->Append(", N = ");
  (*entryQC4)+=(Long_t)nEvtsQC4;
 }
 if(qcCommonHist6)
 {
  avMultQC6 = (qcCommonHist6->GetHistMultInt())->GetMean();
  nEvtsQC6  = (qcCommonHist6->GetHistMultInt())->GetEntries();
 }
 if(entryQC6)
 {  
  entryQC6->Append("M = ");
  (*entryQC6)+=(Long_t)avMultQC6;
  entryQC6->Append(", N = ");
  (*entryQC6)+=(Long_t)nEvtsQC6;
 }
 if(qcCommonHist8)
 {
  avMultQC8 = (qcCommonHist8->GetHistMultInt())->GetMean();
  nEvtsQC8  = (qcCommonHist8->GetHistMultInt())->GetEntries();
 }
 if(entryQC8)
 {
  entryQC8->Append("M = ");
  (*entryQC8)+=(Long_t)avMultQC8;
  entryQC8->Append(", N = ");
  (*entryQC8)+=(Long_t)nEvtsQC8;
 }
 
 //FQD:
 if(fqdCommonHist)
 {
  avMultFQD = (fqdCommonHist->GetHistMultInt())->GetMean();
  nEvtsFQD  = (fqdCommonHist->GetHistMultInt())->GetEntries();
 } 
 if(entryFQD)
 {
  entryFQD->Append("M = ");
  (*entryFQD)+=(Long_t)avMultFQD;
  entryFQD->Append(", N = ");
  (*entryFQD)+=(Long_t)nEvtsFQD;
 }
 
 //LYZ1:
 if(lyz1CommonHist)
 {
  avMultLYZ1 = (lyz1CommonHist->GetHistMultInt())->GetMean();
  nEvtsLYZ1  = (lyz1CommonHist->GetHistMultInt())->GetEntries();
 }
 if(entryLYZ1) 
 {
  entryLYZ1->Append("M = ");
  (*entryLYZ1)+=(Long_t)avMultLYZ1;
  entryLYZ1->Append(", N = ");
  (*entryLYZ1)+=(Long_t)nEvtsLYZ1;
 }
 
 //LYZEP:
 if(lyzepCommonHist)
 {
  avMultcLYZEP = (lyzepCommonHist->GetHistMultInt())->GetMean();
  nEvtsLYZEP  = (lyzepCommonHist->GetHistMultInt())->GetEntries();
 }
 if(entryLYZEP)
 {
  entryLYZEP->Append("M = ");
  (*entryLYZEP)+=(Long_t)avMultLYZEP;
  entryLYZEP->Append(", N = ");
  (*entryLYZEP)+=(Long_t)nEvtsLYZEP;
 }
 
 if(textResults)
 {
  textResults->AddText(entryMC->Data());
  textResults->AddText(entryGFC->Data());
  textResults->AddText(entryQC2->Data());
  textResults->AddText(entryQC4->Data());
  textResults->AddText(entryQC6->Data());
  textResults->AddText(entryQC8->Data());
  textResults->AddText(entryFQD->Data());
  textResults->AddText(entryLYZ1->Data());
  textResults->AddText(entryLYZEP->Data());
 }
 //----------------------------------------------------------------------------------
 
 //----------------------------------------------------------------------------------
 //final drawing for integrated flow:
 TCanvas* intFlowAllCanvas = new TCanvas("Integrated Flow","Integrated Flow",1000,600);
 
 intFlowAllCanvas->Divide(2,1);
 
 //1st pad is for plot:
 (intFlowAllCanvas->cd(1))->SetPad(0.0,0.0,0.75,1.0);
 
 if(intFlowAll)
 {
  if(dMin>0. && dMax>0.)
  {
   (intFlowAll->GetYaxis())->SetRangeUser(0.9744*dMin,1.0144*dMax);
  } else if(dMin<0. && dMax>0.)
    {
     if(!(-1.*dMin<4.*dMax))
     {  
      (intFlowAll->GetYaxis())->SetRangeUser(1.0266*dMin,1.0144*dMax);
     } else {(intFlowAll->GetYaxis())->SetRangeUser(1.1266*dMin,1.0144*dMax);}  
    } else if(dMin<0. && dMax<0.)
      {
       (intFlowAll->GetYaxis())->SetRangeUser(1.0266*dMin,0.9866*dMax);      
      }
  intFlowAll->Draw("E1");
 }
                        
 if(pMesh) pMesh->Draw("LFSAME");
  
 if(flowResults) flowResults->Draw("PSAME");

 //2nd pad is for legend:
 (intFlowAllCanvas->cd(2))->SetPad(0.75,0.0,1.0,1.0);
 
 if(textDefault)
 {
  textDefault->Draw();
  textResults->Draw();
 }
 //----------------------------------------------------------------------------------
 //==================================================================================   
 



 //==================================================================================
 //                            DIFFERENTIAL FLOW
 //==================================================================================
 Int_t iNbinsPt  = AliFlowCommonConstants::GetNbinsPt();
 Double_t dPtMin = AliFlowCommonConstants::GetPtMin();
 Double_t dPtMax = AliFlowCommonConstants::GetPtMax();
 
 //----------------------------------------------------------------------------------
 //cosmetics: the style histogram:
 TH1D *styleHist = new TH1D("styleHist","styleHist",iNbinsPt,dPtMin,dPtMax);
 styleHist->SetTitle("Differential Flow");
 styleHist->SetXTitle("p_{t} [GeV]");
 styleHist->SetYTitle("v_{n}");
 //----------------------------------------------------------------------------------
 
 //----------------------------------------------------------------------------------
 //cosmetics: Monte Carlo error bands for differential flow
 TGraph* pMeshDiffFlow = NULL;
 if(mcepCommonHistRes)
 {
  Int_t nPtsDF = (mcepCommonHistRes->GetHistDiffFlow())->GetNbinsX();
  Double_t binWidth = (mcepCommonHistRes->GetHistDiffFlow())->GetBinWidth(1);//assuming that all bins have the same width
       
  pMeshDiffFlow = new TGraph(2*nPts+1);
  
  Double_t valueMC=0., errorMC=0.;
  for(Int_t i=1;i<nPtsDF+1;i++)
  {
   valueMC = (mcepCommonHistRes->GetHistDiffFlow())->GetBinContent(i);
   errorMC = (mcepCommonHistRes->GetHistDiffFlow())->GetBinError(i);       
   pMeshDiffFlow->SetPoint(i,(i-0.5)*binWidth,valueMC+errorMC);
  }    
  for(Int_t i=nPtsDF+1;i<2*nPtsDF+1;i++)
  {
   valueMC = (mcepCommonHistRes->GetHistDiffFlow())->GetBinContent(2*nPtsDF+1-i);
   errorMC = (mcepCommonHistRes->GetHistDiffFlow())->GetBinError(2*nPtsDF+1-i);       
   pMeshDiffFlow->SetPoint(i,(2*nPtsDF-i+0.5)*binWidth,valueMC-errorMC); 
  }
  pMeshDiffFlow->SetPoint(2*nPtsDF+1,0.5*binWidth,valueMC+errorMC); 
  pMeshDiffFlow->SetFillStyle(meshStyle);
  pMeshDiffFlow->SetFillColor(meshColor);
 } 
 //----------------------------------------------------------------------------------
  
 //MCEP = Monte Carlo Event Plane
 Double_t avMultDiffFlowMC=0.;
 Double_t nEvtsDiffFlowMC=0;
 if(fileMCEP)
 {
  if(mcepCommonHistRes)
  {
   (mcepCommonHistRes->GetHistDiffFlow())->SetMarkerColor(2);
   (mcepCommonHistRes->GetHistDiffFlow())->SetMarkerStyle(20);
  } 
  if(mcepCommonHist)
  {
   avMultDiffFlowMC = (mcepCommonHist->GetHistMultDiff())->GetMean();
   nEvtsDiffFlowMC  = (mcepCommonHist->GetHistMultDiff())->GetEntries();
  } 
 } 

 //GFC = Generating Function Cumulants
 Double_t avMultDiffFlowGFC=0.;
 Double_t nEvtsDiffFlowGFC=0.;
 if(fileGFC)
 {
  if(gfcCommonHistRes2)
  {
   (gfcCommonHistRes2->GetHistDiffFlow())->SetMarkerColor(kViolet+3);
   (gfcCommonHistRes2->GetHistDiffFlow())->SetMarkerStyle(20);
  }
  if(gfcCommonHistRes4)
  { 
   (gfcCommonHistRes4->GetHistDiffFlow())->SetMarkerColor(kViolet-6);
   (gfcCommonHistRes4->GetHistDiffFlow())->SetMarkerStyle(21);
  }
  if(gfcCommonHist)
  {
   avMultDiffFlowGFC = (gfcCommonHist->GetHistMultDiff())->GetMean();
   nEvtsDiffFlowGFC  = (gfcCommonHist->GetHistMultDiff())->GetEntries();
  } 
 }
  
 //QC = Q-cumulants
 Double_t avMultDiffFlowQC2=0., avMultDiffFlowQC4=0.;
 Double_t nEvtsDiffFlowQC2=0., nEvtsDiffFlowQC4=0.;
 if(fileQC)
 {
  //QC{2}
  if(qcCommonHistRes2)
  {
   (qcCommonHistRes2->GetHistDiffFlow())->SetMarkerColor(kOrange+3);
   (qcCommonHistRes2->GetHistDiffFlow())->SetMarkerStyle(20);
  }
  if(qcCommonHist2)
  {
   avMultDiffFlowQC2 = (qcCommonHist2->GetHistMultDiff())->GetMean();
   nEvtsDiffFlowQC2  = (qcCommonHist2->GetHistMultDiff())->GetEntries();
  }
  //QC{4}
  if(qcCommonHistRes4)
  {
   (qcCommonHistRes4->GetHistDiffFlow())->SetMarkerColor(kOrange-6);
   (qcCommonHistRes4->GetHistDiffFlow())->SetMarkerStyle(21);
  }
  if(qcCommonHist4)
  {
   avMultDiffFlowQC4 = (qcCommonHist4->GetHistMultDiff())->GetMean();
   nEvtsDiffFlowQC4  = (qcCommonHist4->GetHistMultDiff())->GetEntries();
  }
 } 

 //LYZ2 = Lee-Yang Zeros (2nd run)
 Double_t avMultDiffFlowLYZ2=0.;
 Double_t nEvtsDiffFlowLYZ2=0;
 if(fileLYZ2)
 {
  if(lyz2CommonHistRes)
  {
   (lyz2CommonHistRes->GetHistDiffFlow())->Scale(0.01);//to be improved
   (lyz2CommonHistRes->GetHistDiffFlow())->SetMarkerColor(kGreen+3);
   (lyz2CommonHistRes->GetHistDiffFlow())->SetMarkerStyle(22);
  } 
  if(lyz2CommonHist)
  {
   avMultDiffFlowLYZ2 = (lyz2CommonHist->GetHistMultDiff())->GetMean();
   nEvtsDiffFlowLYZ2  = (lyz2CommonHist->GetHistMultDiff())->GetEntries();
  } 
 } 

 //----------------------------------------------------------------------------------
 //final drawing for differential flow:
 TCanvas* diffFlowAllCanvas = new TCanvas("Differential Flow","Differential Flow",1000,600);
 
 diffFlowAllCanvas->Divide(2,1);
 
 //1st pad is for plot:
 (diffFlowAllCanvas->cd(1))->SetPad(0.0,0.0,0.75,1.0);
 
 if(styleHist)
 {
  styleHist->Draw();
 }
 if(pMeshDiffFlow)
 {
  pMeshDiffFlow->Draw("LFSAME");
 }
 //MC 
 if(mcepCommonHistRes)
 { 
  //(mcepCommonHistRes->GetHistDiffFlow())->Draw("E1PSAME");
 }
 //GFC
 if(gfcCommonHistRes2)
 { 
  (gfcCommonHistRes2->GetHistDiffFlow())->Draw("E1PSAME"); 
 } 
 if(gfcCommonHistRes4)
 { 
  (gfcCommonHistRes4->GetHistDiffFlow())->Draw("E1PSAME"); 
 } 
 //QC
 if(qcCommonHistRes2)
 { 
  (qcCommonHistRes2->GetHistDiffFlow())->Draw("E1PSAME");
 }
 if(qcCommonHistRes4)
 { 
  (qcCommonHistRes4->GetHistDiffFlow())->Draw("E1PSAME");
 }
 //LYZ2
 if(lyz2CommonHistRes)
 { 
  (lyz2CommonHistRes->GetHistDiffFlow())->Draw("E1PSAME");
 }
 
 //2nd pad is for legend:
 (diffFlowAllCanvas->cd(2))->SetPad(0.75,0.0,1.0,1.0);

 TLegend* legendDiffFlow = new TLegend(0.02,0.25,0.97,0.75);
 legendDiffFlow->SetTextFont(72);
 legendDiffFlow->SetTextSize(0.06);

 //legend's entries:
 TString *entryDiffMC   = new TString("MC ....... ");
 TString *entryDiffGFC2 = new TString("GFC{2} ... ");
 TString *entryDiffGFC4 = new TString("GFC{4} ... ");
 TString *entryDiffQC2  = new TString("QC{2} .... ");
 TString *entryDiffQC4  = new TString("QC{4} .... ");
 TString *entryDiffLYZ2 = new TString("LYZ ...... ");
 
 //MC
 if(mcepCommonHistRes)
 {
  (mcepCommonHistRes->GetHistDiffFlow())->SetFillStyle(meshStyle);
  (mcepCommonHistRes->GetHistDiffFlow())->SetFillColor(meshColor);
  entryDiffMC->Append("M = ");
  (*entryDiffMC)+=(Long_t)avMultDiffFlowMC;
  entryDiffMC->Append(", N = ");
  (*entryDiffMC)+=(Long_t)nEvtsDiffFlowMC; 
  legendDiffFlow->AddEntry(mcepCommonHistRes->GetHistDiffFlow(),entryDiffMC->Data(),"f");
 }

 //GFC
 if(gfcCommonHistRes2)
 {
  entryDiffGFC2->Append("M = ");
  (*entryDiffGFC2)+=(Long_t)avMultDiffFlowGFC;
  entryDiffGFC2->Append(", N = ");
  (*entryDiffGFC2)+=(Long_t)nEvtsDiffFlowGFC; 
  legendDiffFlow->AddEntry(gfcCommonHistRes2->GetHistDiffFlow(),entryDiffGFC2->Data(),"p");
 }
 if(gfcCommonHistRes4)
 {
  entryDiffGFC4->Append("M = ");
  (*entryDiffGFC4)+=(Long_t)avMultDiffFlowGFC;
  entryDiffGFC4->Append(", N = ");
  (*entryDiffGFC4)+=(Long_t)nEvtsDiffFlowGFC; 
  legendDiffFlow->AddEntry(gfcCommonHistRes4->GetHistDiffFlow(),entryDiffGFC4->Data(),"p");
 }
  
 //QC
 if(qcCommonHistRes2)
 {
  entryDiffQC2->Append("M = ");
  (*entryDiffQC2)+=(Long_t)avMultDiffFlowQC2;
  entryDiffQC2->Append(", N = ");
  (*entryDiffQC2)+=(Long_t)nEvtsDiffFlowQC2; 
  legendDiffFlow->AddEntry(qcCommonHistRes2->GetHistDiffFlow(),entryDiffQC2->Data(),"p");
 }
 if(qcCommonHistRes4)
 {
  entryDiffQC4->Append("M = ");
  (*entryDiffQC4)+=(Long_t)avMultDiffFlowQC4;
  entryDiffQC4->Append(", N = ");
  (*entryDiffQC4)+=(Long_t)nEvtsDiffFlowQC4; 
  legendDiffFlow->AddEntry(qcCommonHistRes4->GetHistDiffFlow(),entryDiffQC4->Data(),"p");
 }
 
 //LYZ
 if(lyz2CommonHistRes)
 {
  entryDiffLYZ2->Append("M = ");
  (*entryDiffLYZ2)+=(Long_t)avMultDiffFlowLYZ2;
  entryDiffLYZ2->Append(", N = ");
  (*entryDiffLYZ2)+=(Long_t)nEvtsDiffFlowLYZ2; 
  legendDiffFlow->AddEntry(lyz2CommonHistRes->GetHistDiffFlow(),entryDiffLYZ2->Data(),"p");
 }

 //drawing finally the legend in the 2nd pad:     
 if(legendDiffFlow)
 {
  legendDiffFlow->Draw();
 }
 
 //=====================================================================================

}
