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
 const TString type = "";
 
 //open the output files:
 TString inputFileNameMCEP = "outputMCEPanalysis";
 TFile* fileMCEP = NULL;
 fileMCEP = TFile::Open(((inputFileNameMCEP.Append(type)).Append(".root")).Data(), "READ"); 
 
 TString inputFileNameLYZ1 = "outputLYZ1analysis";
 TFile* fileLYZ1 = NULL;
 fileLYZ1 = TFile::Open(((inputFileNameLYZ1.Append(type)).Append(".root")).Data(), "READ"); 

 TString inputFileNameLYZ2 = "outputLYZ2analysis";
 TFile* fileLYZ2 = NULL;
 fileLYZ2 = TFile::Open(((inputFileNameLYZ2.Append(type)).Append(".root")).Data(), "READ"); 
 
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
 Double_t flowValue[nMethods] = {0.};//to be removed
 Double_t flowError[nMethods] = {0.};//to be removed
 Double_t flowValueRP[nMethods] = {0.};
 Double_t flowErrorRP[nMethods] = {0.};
 Double_t flowValuePOI[nMethods] = {0.};
 Double_t flowErrorPOI[nMethods] = {0.};
 
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
    flowValue[0] = (mcepCommonHistRes->GetHistIntFlow())->GetBinContent(1);//to be removed
    flowError[0] = (mcepCommonHistRes->GetHistIntFlow())->GetBinError(1);//to be removed
    flowValueRP[0] = (mcepCommonHistRes->GetHistIntFlowRP())->GetBinContent(1);
    flowErrorRP[0] = (mcepCommonHistRes->GetHistIntFlowRP())->GetBinError(1);
    flowValuePOI[0] = (mcepCommonHistRes->GetHistIntFlowPOI())->GetBinContent(1);
    flowErrorPOI[0] = (mcepCommonHistRes->GetHistIntFlowPOI())->GetBinError(1);
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
    flowValue[10] = (lyz1CommonHistRes->GetHistIntFlow())->GetBinContent(1);//to be removed
    flowError[10] = (lyz1CommonHistRes->GetHistIntFlow())->GetBinError(1);//to be removed
    flowValueRP[10] = (lyz1CommonHistRes->GetHistIntFlowRP())->GetBinContent(1);
    flowErrorRP[10] = (lyz1CommonHistRes->GetHistIntFlowRP())->GetBinError(1);
    flowValuePOI[10] = (lyz1CommonHistRes->GetHistIntFlowPOI())->GetBinContent(1);
    flowErrorPOI[10] = (lyz1CommonHistRes->GetHistIntFlowPOI())->GetBinError(1);
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
    flowValue[11] = (lyzepCommonHistRes->GetHistIntFlow())->GetBinContent(1);//to be removed
    //flowError[11] = (lyzepCommonHistRes->GetHistIntFlow())->GetBinError(1);//to be removed
    flowValueRP[11] = (lyzepCommonHistRes->GetHistIntFlowRP())->GetBinContent(1);
    //flowErrorRP[11] = (lyzepCommonHistRes->GetHistIntFlowRP())->GetBinError(1);
    flowValuePOI[11] = (lyzepCommonHistRes->GetHistIntFlowPOI())->GetBinContent(1);
    //flowErrorPOI[11] = (lyzepCommonHistRes->GetHistIntFlowPOI())->GetBinError(1);
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
    flowValue[9] = (fqdCommonHistRes->GetHistIntFlow())->GetBinContent(1);//to be removed
    flowError[9] = (fqdCommonHistRes->GetHistIntFlow())->GetBinError(1);//to be removed
    flowValueRP[9] = (fqdCommonHistRes->GetHistIntFlowRP())->GetBinContent(1);
    flowErrorRP[9] = (fqdCommonHistRes->GetHistIntFlowRP())->GetBinError(1);
    flowValuePOI[9] = (fqdCommonHistRes->GetHistIntFlowPOI())->GetBinContent(1);
    flowErrorPOI[9] = (fqdCommonHistRes->GetHistIntFlowPOI())->GetBinError(1);
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
    flowValue[1] = (gfcCommonHistRes2->GetHistIntFlow())->GetBinContent(1);//to be removed
    flowError[1] = (gfcCommonHistRes2->GetHistIntFlow())->GetBinError(1);//to be removed
    flowValueRP[1] = (gfcCommonHistRes2->GetHistIntFlowRP())->GetBinContent(1);
    flowErrorRP[1] = (gfcCommonHistRes2->GetHistIntFlowRP())->GetBinError(1);
    flowValuePOI[1] = (gfcCommonHistRes2->GetHistIntFlowPOI())->GetBinContent(1);
    flowErrorPOI[1] = (gfcCommonHistRes2->GetHistIntFlowPOI())->GetBinError(1);
   }
   gfcCommonHistRes4 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults4thOrderGFC"));
   if(gfcCommonHistRes4) 
   {
    flowValue[3] = (gfcCommonHistRes4->GetHistIntFlow())->GetBinContent(1);//to be removed
    flowError[3] = (gfcCommonHistRes4->GetHistIntFlow())->GetBinError(1);//to be removed
    flowValueRP[3] = (gfcCommonHistRes4->GetHistIntFlowRP())->GetBinContent(1);
    flowErrorRP[3] = (gfcCommonHistRes4->GetHistIntFlowRP())->GetBinError(1);
    flowValuePOI[3] = (gfcCommonHistRes4->GetHistIntFlowPOI())->GetBinContent(1);
    flowErrorPOI[3] = (gfcCommonHistRes4->GetHistIntFlowPOI())->GetBinError(1);
   }
   gfcCommonHistRes6 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults6thOrderGFC"));
   if(gfcCommonHistRes6) 
   {
    flowValue[5] = (gfcCommonHistRes6->GetHistIntFlow())->GetBinContent(1);//to be removed
    flowError[5] = (gfcCommonHistRes6->GetHistIntFlow())->GetBinError(1);//to be removed
    flowValueRP[5] = (gfcCommonHistRes6->GetHistIntFlowRP())->GetBinContent(1);
    flowErrorRP[5] = (gfcCommonHistRes6->GetHistIntFlowRP())->GetBinError(1);
    flowValuePOI[5] = (gfcCommonHistRes6->GetHistIntFlowPOI())->GetBinContent(1);
    flowErrorPOI[5] = (gfcCommonHistRes6->GetHistIntFlowPOI())->GetBinError(1);
   }
   gfcCommonHistRes8 = dynamic_cast<AliFlowCommonHistResults*> (pListGFC->FindObject("AliFlowCommonHistResults8thOrderGFC"));
   if(gfcCommonHistRes8) 
   {
    flowValue[7] = (gfcCommonHistRes8->GetHistIntFlow())->GetBinContent(1);//to be removed
    flowError[7] = (gfcCommonHistRes8->GetHistIntFlow())->GetBinError(1);//to be removed
    flowValueRP[7] = (gfcCommonHistRes8->GetHistIntFlowRP())->GetBinContent(1);
    flowErrorRP[7] = (gfcCommonHistRes8->GetHistIntFlowRP())->GetBinError(1);
    flowValuePOI[7] = (gfcCommonHistRes8->GetHistIntFlowPOI())->GetBinContent(1);
    flowErrorPOI[7] = (gfcCommonHistRes8->GetHistIntFlowPOI())->GetBinError(1);
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
    flowValue[2] = (qcCommonHistRes2->GetHistIntFlow())->GetBinContent(1);//to be removed
    //flowError[2] = (qcCommonHistRes2->GetHistIntFlow())->GetBinError(1);//to be removed
    flowValueRP[2] = (qcCommonHistRes2->GetHistIntFlowRP())->GetBinContent(1);
    //flowErrorRP[2] = (qcCommonHistRes2->GetHistIntFlowRP())->GetBinError(1);
    flowValuePOI[2] = (qcCommonHistRes2->GetHistIntFlowPOI())->GetBinContent(1);
    //flowErrorPOI[2] = (qcCommonHistRes2->GetHistIntFlowPOI())->GetBinError(1);
   }
   qcCommonHist4 = dynamic_cast<AliFlowCommonHist*> (pListQC->FindObject("AliFlowCommonHist4thOrderQC"));
   qcCommonHistRes4 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults4thOrderQC"));
   if(qcCommonHistRes4) 
   {
    flowValue[4] = (qcCommonHistRes4->GetHistIntFlow())->GetBinContent(1);//to be removed
    //flowError[4] = (qcCommonHistRes4->GetHistIntFlow())->GetBinError(1);//to be removed
    flowValueRP[4] = (qcCommonHistRes4->GetHistIntFlowRP())->GetBinContent(1);
    //flowErrorRP[4] = (qcCommonHistRes4->GetHistIntFlowRP())->GetBinError(1);
    flowValuePOI[4] = (qcCommonHistRes4->GetHistIntFlowPOI())->GetBinContent(1);
    //flowErrorPOI[4] = (qcCommonHistRes4->GetHistIntFlowPOI())->GetBinError(1);
   }
   qcCommonHist6 = dynamic_cast<AliFlowCommonHist*> (pListQC->FindObject("AliFlowCommonHist6thOrderQC"));
   qcCommonHistRes6 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults6thOrderQC"));
   if(qcCommonHistRes6) 
   {
    flowValue[6] = (qcCommonHistRes6->GetHistIntFlow())->GetBinContent(1);//to be removed
    //flowError[6] = (qcCommonHistRes6->GetHistIntFlow())->GetBinError(1);//to be removed
    flowValueRP[6] = (qcCommonHistRes6->GetHistIntFlowRP())->GetBinContent(1);
    //flowErrorRP[6] = (qcCommonHistRes6->GetHistIntFlowRP())->GetBinError(1);
    flowValuePOI[6] = (qcCommonHistRes6->GetHistIntFlowPOI())->GetBinContent(1);
    //flowErrorPOI[6] = (qcCommonHistRes6->GetHistIntFlowPOI())->GetBinError(1);
   }
   qcCommonHist8 = dynamic_cast<AliFlowCommonHist*> (pListQC->FindObject("AliFlowCommonHist8thOrderQC"));
   qcCommonHistRes8 = dynamic_cast<AliFlowCommonHistResults*> (pListQC->FindObject("AliFlowCommonHistResults8thOrderQC"));
   if(qcCommonHistRes8) 
   {
    flowValue[8] = (qcCommonHistRes8->GetHistIntFlow())->GetBinContent(1);//to be removed
    //flowError[8] = (qcCommonHistRes8->GetHistIntFlow())->GetBinError(1);//to be removed
    flowValueRP[8] = (qcCommonHistRes8->GetHistIntFlowRP())->GetBinContent(1);
    //flowErrorRP[8] = (qcCommonHistRes8->GetHistIntFlowRP())->GetBinError(1);
    flowValuePOI[8] = (qcCommonHistRes8->GetHistIntFlowPOI())->GetBinContent(1);
    //flowErrorPOI[8] = (qcCommonHistRes8->GetHistIntFlowPOI())->GetBinError(1);
   }
  }
 }        
 
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //              !!!!  to be removed  !!!!
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
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
 //RP:
 Double_t dMaxRP=flowValueRP[0]+flowErrorRP[0];
 Double_t dMinRP=flowValueRP[0]-flowErrorRP[0];

 for(Int_t i=1;i<nMethods;i++)
 {
  if(!(flowValueRP[i]==0. && flowErrorRP[i]==0.)) 
  {
   if(dMaxRP<flowValueRP[i]+flowErrorRP[i]) dMaxRP=flowValueRP[i]+flowErrorRP[i];
   if(dMinRP>flowValueRP[i]-flowErrorRP[i]) dMinRP=flowValueRP[i]-flowErrorRP[i];
  } 
 }  

 //POI:
 Double_t dMaxPOI=flowValuePOI[0]+flowErrorPOI[0];
 Double_t dMinPOI=flowValuePOI[0]-flowErrorPOI[0];

 for(Int_t i=1;i<nMethods;i++)
 {
  if(!(flowValuePOI[i]==0. && flowErrorPOI[i]==0.)) 
  {
   if(dMaxPOI<flowValuePOI[i]+flowErrorPOI[i]) dMaxPOI=flowValuePOI[i]+flowErrorPOI[i];
   if(dMinPOI>flowValuePOI[i]-flowErrorPOI[i]) dMinPOI=flowValuePOI[i]-flowErrorPOI[i];
  } 
 }  
 
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
 //                    !!!!  to be removed  !!!!
 TGraph* flowResults = new TGraphErrors(nMethods, x, flowValue, xError, flowError);
 
 flowResults->SetMarkerStyle(markerStyle);
 flowResults->SetMarkerColor(markerColor);
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
 //RP:
 TGraph* flowResultsRP = new TGraphErrors(nMethods, x, flowValueRP, xError, flowErrorRP);
 
 flowResultsRP->SetMarkerStyle(markerStyle);
 flowResultsRP->SetMarkerColor(markerColor);
 
 //POI:
 TGraph* flowResultsPOI = new TGraphErrors(nMethods, x, flowValuePOI, xError, flowErrorPOI);
 
 flowResultsPOI->SetMarkerStyle(markerStyle);
 flowResultsPOI->SetMarkerColor(markerColor);
 
 //-----------------------------------------------------------------------------------
  
 //----------------------------------------------------------------------------------
 //cosmetics: mesh for MC error bands (integrated flow)
 TGraph* pMesh = NULL;//to be removed
 TGraph* pMeshRP = NULL;
 TGraph* pMeshPOI = NULL;
 
 if(intFlowAll && mcepCommonHistRes)
 {
  //Double_t valueMC = intFlowAll->GetBinContent(1);
  //Double_t errorMC = intFlowAll->GetBinError(1);  
  //Int_t nPts       = intFlowAll->GetNbinsX();     

  Int_t nPts       = nMethods;
  Double_t valueMC = flowValue[0];//to be removed
  Double_t errorMC = flowError[0];//to be removed  
  Double_t valueMCRP = flowValueRP[0];
  Double_t errorMCRP = flowErrorRP[0]; 
  Double_t valueMCPOI = flowValuePOI[0];
  Double_t errorMCPOI = flowErrorPOI[0]; 
              
  pMesh = new TGraph(nPts);//to be removed
  pMeshRP = new TGraph(nPts);
  pMeshPOI = new TGraph(nPts);
  
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  //             !!!! to be removed !!!!
  pMesh->SetPoint(1,0,valueMC+errorMC);
  pMesh->SetPoint(2,nPts+1,valueMC+errorMC);
  pMesh->SetPoint(3,nPts+1,valueMC-errorMC);
  pMesh->SetPoint(4,0,valueMC-errorMC);
  pMesh->SetPoint(5,0,valueMC+errorMC);
  
  pMesh->SetFillStyle(meshStyle);
  pMesh->SetFillColor(meshColor);
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  //RP:
  pMeshRP->SetPoint(1,0,valueMCRP+errorMCRP);
  pMeshRP->SetPoint(2,nPts+1,valueMCRP+errorMCRP);
  pMeshRP->SetPoint(3,nPts+1,valueMCRP-errorMCRP);
  pMeshRP->SetPoint(4,0,valueMCRP-errorMCRP);
  pMeshRP->SetPoint(5,0,valueMCRP+errorMCRP);
  
  pMeshRP->SetFillStyle(meshStyle);
  pMeshRP->SetFillColor(meshColor);

  //POI:
  pMeshPOI->SetPoint(1,0,valueMCPOI+errorMCPOI);
  pMeshPOI->SetPoint(2,nPts+1,valueMCPOI+errorMCPOI);
  pMeshPOI->SetPoint(3,nPts+1,valueMCPOI-errorMCPOI);
  pMeshPOI->SetPoint(4,0,valueMCPOI-errorMCPOI);
  pMeshPOI->SetPoint(5,0,valueMCPOI+errorMCPOI);
  
  pMeshPOI->SetFillStyle(meshStyle);
  pMeshPOI->SetFillColor(meshColor);     
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
 TPaveText *textResults = new TPaveText(0.05,0.12,0.95,0.70,"NDC");//to be removed
 textResults->SetTextFont(72);//to be removed
 textResults->SetTextSize(0.06);//to be removed
 //textResults->SetLineColor(kFALSE);
 //textResults->SetShadowColor(kFALSE);

 //results (RP):
 TPaveText *textResultsRP = new TPaveText(0.05,0.12,0.95,0.70,"NDC");
 textResultsRP->SetTextFont(72);
 textResultsRP->SetTextSize(0.06);

 //results (POI):
 TPaveText *textResultsPOI = new TPaveText(0.05,0.12,0.95,0.70,"NDC");
 textResultsPOI->SetTextFont(72);
 textResultsPOI->SetTextSize(0.06);
 
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //      !!!! to be removed !!!!!!                                      
 TString *entryMC    = new TString("MC ...... ");
 TString *entryGFC   = new TString("GFC ..... "); 
 TString *entryQC2   = new TString("QC{2} ... ");
 TString *entryQC4   = new TString("QC{4} ... ");
 TString *entryQC6   = new TString("QC{6} ... ");
 TString *entryQC8   = new TString("QC{8} ... ");
 TString *entryFQD   = new TString("FQD ..... "); 
 TString *entryLYZ1  = new TString("LYZ ..... "); 
 TString *entryLYZEP = new TString("LYZEP ... "); 
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
 //RP: 
 TString *entryMCRP    = new TString("MC ...... ");
 TString *entryGFCRP   = new TString("GFC ..... "); 
 TString *entryQC2RP   = new TString("QC{2} ... ");
 TString *entryQC4RP   = new TString("QC{4} ... ");
 TString *entryQC6RP   = new TString("QC{6} ... ");
 TString *entryQC8RP   = new TString("QC{8} ... ");
 TString *entryFQDRP   = new TString("FQD ..... "); 
 TString *entryLYZ1RP  = new TString("LYZ ..... "); 
 TString *entryLYZEPRP = new TString("LYZEP ... "); 

 //POI: 
 TString *entryMCPOI    = new TString("MC ...... ");
 TString *entryGFCPOI   = new TString("GFC ..... "); 
 TString *entryQC2POI   = new TString("QC{2} ... ");
 TString *entryQC4POI   = new TString("QC{4} ... ");
 TString *entryQC6POI   = new TString("QC{6} ... ");
 TString *entryQC8POI   = new TString("QC{8} ... ");
 TString *entryFQDPOI   = new TString("FQD ..... "); 
 TString *entryLYZ1POI  = new TString("LYZ ..... "); 
 TString *entryLYZEPPOI = new TString("LYZEP ... "); 
 
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //                  !!!! to be removed !!!!
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
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 //RP:
 Double_t avMultMCRP=0.;
 Long_t nEvtsMCRP=0;

 Double_t avMultGFCRP=0.;
 Long_t nEvtsGFCRP=0;
  
 Double_t avMultQC2RP=0., avMultQC4RP=0., avMultQC6RP=0., avMultQC8RP=0.;
 Long_t nEvtsQC2RP=0, nEvtsQC4RP=0, nEvtsQC6RP=0, nEvtsQC8RP=0;

 Double_t avMultFQDRP=0.;
 Long_t nEvtsFQDRP=0;

 Double_t avMultLYZ1RP=0.;
 Long_t nEvtsLYZ1RP=0;
 
 Double_t avMultLYZEPRP=0.;
 Long_t nEvtsLYZEPRP=0;
 
 //POI:
 Double_t avMultMCPOI=0.;
 Long_t nEvtsMCPOI=0;

 Double_t avMultGFCPOI=0.;
 Long_t nEvtsGFCPOI=0;
  
 Double_t avMultQC2POI=0., avMultQC4POI=0., avMultQC6POI=0., avMultQC8POI=0.;
 Long_t nEvtsQC2POI=0, nEvtsQC4POI=0, nEvtsQC6POI=0, nEvtsQC8POI=0;

 Double_t avMultFQDPOI=0.;
 Long_t nEvtsFQDPOI=0;

 Double_t avMultLYZ1POI=0.;
 Long_t nEvtsLYZ1POI=0;
 
 Double_t avMultLYZEPPOI=0.;
 Long_t nEvtsLYZEPPOI=0;
 
 //MC:  
 if(mcepCommonHist)
 {
  avMultMC = (mcepCommonHist->GetHistMultInt())->GetMean();//to be removed
  nEvtsMC  = (mcepCommonHist->GetHistMultInt())->GetEntries();//to be removed
  avMultMCRP = (mcepCommonHist->GetHistMultInt())->GetMean();
  nEvtsMCRP  = (mcepCommonHist->GetHistMultInt())->GetEntries();
  avMultMCPOI = (mcepCommonHist->GetHistMultDiff())->GetMean();
  nEvtsMCPOI  = (mcepCommonHist->GetHistMultDiff())->GetEntries();
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //   !!!! to be removed !!!!
 if(entryMC)
 {   
  entryMC->Append("M = ");
  (*entryMC)+=(Long_t)avMultMC;
  entryMC->Append(", N = ");
  (*entryMC)+=(Long_t)nEvtsMC;
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 if(entryMCRP)
 {   
  entryMCRP->Append("M = ");
  (*entryMCRP)+=(Long_t)avMultMCRP;
  entryMCRP->Append(", N = ");
  (*entryMCRP)+=(Long_t)nEvtsMCRP;
 }
 if(entryMCPOI)
 {   
  entryMCPOI->Append("M = ");
  (*entryMCPOI)+=(Long_t)avMultMCPOI;
  entryMCPOI->Append(", N = ");
  (*entryMCPOI)+=(Long_t)nEvtsMCPOI;
 }
  
 //GFC:
 if(gfcCommonHist)
 {
  avMultGFC = (gfcCommonHist->GetHistMultInt())->GetMean();//to be removed
  nEvtsGFC  = (gfcCommonHist->GetHistMultInt())->GetEntries();//to be removed
  avMultGFCRP = (gfcCommonHist->GetHistMultInt())->GetMean();
  nEvtsGFCRP  = (gfcCommonHist->GetHistMultInt())->GetEntries();
  avMultGFCPOI = (gfcCommonHist->GetHistMultDiff())->GetMean();
  nEvtsGFCPOI  = (gfcCommonHist->GetHistMultDiff())->GetEntries();
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //   !!!! to be removed !!!! 
 if(entryGFC)
 { 
  entryGFC->Append("M = ");
  (*entryGFC)+=(Long_t)avMultGFC;
  entryGFC->Append(", N = ");
  (*entryGFC)+=(Long_t)nEvtsGFC;
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 if(entryGFCRP)
 { 
  entryGFCRP->Append("M = ");
  (*entryGFCRP)+=(Long_t)avMultGFCRP;
  entryGFCRP->Append(", N = ");
  (*entryGFCRP)+=(Long_t)nEvtsGFCRP;
 }
 if(entryGFCPOI)
 { 
  entryGFCPOI->Append("M = ");
  (*entryGFCPOI)+=(Long_t)avMultGFCPOI;
  entryGFCPOI->Append(", N = ");
  (*entryGFCPOI)+=(Long_t)nEvtsGFCPOI;
 }
 
 //QC:
 if(qcCommonHist2)
 {
  avMultQC2 = (qcCommonHist2->GetHistMultInt())->GetMean();//to be removed
  nEvtsQC2  = (qcCommonHist2->GetHistMultInt())->GetEntries();//to be removed
  avMultQC2RP = (qcCommonHist2->GetHistMultInt())->GetMean();
  nEvtsQC2RP  = (qcCommonHist2->GetHistMultInt())->GetEntries();
  avMultQC2POI = (qcCommonHist2->GetHistMultDiff())->GetMean();
  nEvtsQC2POI  = (qcCommonHist2->GetHistMultDiff())->GetEntries();
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //   !!!! to be removed !!!!   
 if(entryQC2)
 { 
  entryQC2->Append("M = ");
  (*entryQC2)+=(Long_t)avMultQC2;
  entryQC2->Append(", N = ");
  (*entryQC2)+=(Long_t)nEvtsQC2;
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 if(entryQC2RP)
 { 
  entryQC2RP->Append("M = ");
  (*entryQC2RP)+=(Long_t)avMultQC2RP;
  entryQC2RP->Append(", N = ");
  (*entryQC2RP)+=(Long_t)nEvtsQC2RP;
 } 
 if(entryQC2POI)
 { 
  entryQC2POI->Append("M = ");
  (*entryQC2POI)+=(Long_t)avMultQC2POI;
  entryQC2POI->Append(", N = ");
  (*entryQC2POI)+=(Long_t)nEvtsQC2POI;
 } 

 if(qcCommonHist4)
 {
  avMultQC4 = (qcCommonHist4->GetHistMultInt())->GetMean();//to be removed
  nEvtsQC4  = (qcCommonHist4->GetHistMultInt())->GetEntries();//to be removed
  avMultQC4RP = (qcCommonHist4->GetHistMultInt())->GetMean();
  nEvtsQC4RP  = (qcCommonHist4->GetHistMultInt())->GetEntries();
  avMultQC4POI = (qcCommonHist4->GetHistMultDiff())->GetMean();
  nEvtsQC4POI  = (qcCommonHist4->GetHistMultDiff())->GetEntries();
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //   !!!! to be removed !!!!  
 if(entryQC4)
 {
  entryQC4->Append("M = ");
  (*entryQC4)+=(Long_t)avMultQC4;
  entryQC4->Append(", N = ");
  (*entryQC4)+=(Long_t)nEvtsQC4;
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 if(entryQC4RP)
 {
  entryQC4RP->Append("M = ");
  (*entryQC4RP)+=(Long_t)avMultQC4RP;
  entryQC4RP->Append(", N = ");
  (*entryQC4RP)+=(Long_t)nEvtsQC4RP;
 }
 if(entryQC4POI)
 {
  entryQC4POI->Append("M = ");
  (*entryQC4POI)+=(Long_t)avMultQC4POI;
  entryQC4POI->Append(", N = ");
  (*entryQC4POI)+=(Long_t)nEvtsQC4POI;
 }
   
 if(qcCommonHist6)
 {
  avMultQC6 = (qcCommonHist6->GetHistMultInt())->GetMean();//to be removed
  nEvtsQC6  = (qcCommonHist6->GetHistMultInt())->GetEntries();//to be removed
  avMultQC6RP = (qcCommonHist6->GetHistMultInt())->GetMean();
  nEvtsQC6RP  = (qcCommonHist6->GetHistMultInt())->GetEntries();
  avMultQC6POI = (qcCommonHist6->GetHistMultDiff())->GetMean();
  nEvtsQC6POI  = (qcCommonHist6->GetHistMultDiff())->GetEntries();
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //   !!!! to be removed !!!!   
 if(entryQC6)
 {  
  entryQC6->Append("M = ");
  (*entryQC6)+=(Long_t)avMultQC6;
  entryQC6->Append(", N = ");
  (*entryQC6)+=(Long_t)nEvtsQC6;
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 if(entryQC6RP)
 {  
  entryQC6RP->Append("M = ");
  (*entryQC6RP)+=(Long_t)avMultQC6RP;
  entryQC6RP->Append(", N = ");
  (*entryQC6RP)+=(Long_t)nEvtsQC6RP;
 }
 if(entryQC6POI)
 {  
  entryQC6POI->Append("M = ");
  (*entryQC6POI)+=(Long_t)avMultQC6POI;
  entryQC6POI->Append(", N = ");
  (*entryQC6POI)+=(Long_t)nEvtsQC6POI;
 }
   
 if(qcCommonHist8)
 {
  avMultQC8 = (qcCommonHist8->GetHistMultInt())->GetMean();//to be removed
  nEvtsQC8  = (qcCommonHist8->GetHistMultInt())->GetEntries();//to be removed
  avMultQC8RP = (qcCommonHist8->GetHistMultInt())->GetMean();
  nEvtsQC8RP  = (qcCommonHist8->GetHistMultInt())->GetEntries();
  avMultQC8POI = (qcCommonHist8->GetHistMultDiff())->GetMean();
  nEvtsQC8POI  = (qcCommonHist8->GetHistMultDiff())->GetEntries();    
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //   !!!! to be removed !!!!    
 if(entryQC8)
 {
  entryQC8->Append("M = ");
  (*entryQC8)+=(Long_t)avMultQC8;
  entryQC8->Append(", N = ");
  (*entryQC8)+=(Long_t)nEvtsQC8;
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 if(entryQC8RP)
 {
  entryQC8RP->Append("M = ");
  (*entryQC8RP)+=(Long_t)avMultQC8RP;
  entryQC8RP->Append(", N = ");
  (*entryQC8RP)+=(Long_t)nEvtsQC8RP;
 }
 if(entryQC8POI)
 {
  entryQC8POI->Append("M = ");
  (*entryQC8POI)+=(Long_t)avMultQC8POI;
  entryQC8POI->Append(", N = ");
  (*entryQC8POI)+=(Long_t)nEvtsQC8POI;
 }
  
 //FQD:
 if(fqdCommonHist)
 {
  avMultFQD = (fqdCommonHist->GetHistMultInt())->GetMean();//to be removed
  nEvtsFQD  = (fqdCommonHist->GetHistMultInt())->GetEntries();//to be removed
  avMultFQDRP = (fqdCommonHist->GetHistMultInt())->GetMean();
  nEvtsFQDRP  = (fqdCommonHist->GetHistMultInt())->GetEntries();
  avMultFQDPOI = (fqdCommonHist->GetHistMultDiff())->GetMean();
  nEvtsFQDPOI  = (fqdCommonHist->GetHistMultDiff())->GetEntries();
 } 
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //   !!!! to be removed !!!!  
 if(entryFQD)
 {
  entryFQD->Append("M = ");
  (*entryFQD)+=(Long_t)avMultFQD;
  entryFQD->Append(", N = ");
  (*entryFQD)+=(Long_t)nEvtsFQD;
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 if(entryFQDRP)
 {
  entryFQDRP->Append("M = ");
  (*entryFQDRP)+=(Long_t)avMultFQDRP;
  entryFQDRP->Append(", N = ");
  (*entryFQDRP)+=(Long_t)nEvtsFQDRP;
 }
 if(entryFQDPOI)
 {
  entryFQDPOI->Append("M = ");
  (*entryFQDPOI)+=(Long_t)avMultFQDPOI;
  entryFQDPOI->Append(", N = ");
  (*entryFQDPOI)+=(Long_t)nEvtsFQDPOI;
 }  
  
 //LYZ1:
 if(lyz1CommonHist)
 {
  avMultLYZ1 = (lyz1CommonHist->GetHistMultInt())->GetMean();//to be removed
  nEvtsLYZ1  = (lyz1CommonHist->GetHistMultInt())->GetEntries();//to be removed
  avMultLYZ1RP = (lyz1CommonHist->GetHistMultInt())->GetMean();
  nEvtsLYZ1RP  = (lyz1CommonHist->GetHistMultInt())->GetEntries();
  avMultLYZ1POI = (lyz1CommonHist->GetHistMultDiff())->GetMean();
  nEvtsLYZ1POI  = (lyz1CommonHist->GetHistMultDiff())->GetEntries();
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //   !!!! to be removed !!!!   
 if(entryLYZ1) 
 {
  entryLYZ1->Append("M = ");
  (*entryLYZ1)+=(Long_t)avMultLYZ1;
  entryLYZ1->Append(", N = ");
  (*entryLYZ1)+=(Long_t)nEvtsLYZ1;
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
 if(entryLYZ1RP) 
 {
  entryLYZ1RP->Append("M = ");
  (*entryLYZ1RP)+=(Long_t)avMultLYZ1RP;
  entryLYZ1RP->Append(", N = ");
  (*entryLYZ1RP)+=(Long_t)nEvtsLYZ1RP;
 }
 if(entryLYZ1POI) 
 {
  entryLYZ1POI->Append("M = ");
  (*entryLYZ1POI)+=(Long_t)avMultLYZ1POI;
  entryLYZ1POI->Append(", N = ");
  (*entryLYZ1POI)+=(Long_t)nEvtsLYZ1POI;
 }
 
 //LYZEP:
 if(lyzepCommonHist)
 {
  avMultcLYZEP = (lyzepCommonHist->GetHistMultInt())->GetMean();//to be removed
  nEvtsLYZEP  = (lyzepCommonHist->GetHistMultInt())->GetEntries();//to be removed
  avMultcLYZEPRP = (lyzepCommonHist->GetHistMultInt())->GetMean();
  nEvtsLYZEPRP  = (lyzepCommonHist->GetHistMultInt())->GetEntries();
  avMultcLYZEPPOI = (lyzepCommonHist->GetHistMultDiff())->GetMean();
  nEvtsLYZEPPOI  = (lyzepCommonHist->GetHistMultDiff())->GetEntries();    
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //   !!!! to be removed !!!!   
 if(entryLYZEP)
 {
  entryLYZEP->Append("M = ");
  (*entryLYZEP)+=(Long_t)avMultLYZEP;
  entryLYZEP->Append(", N = ");
  (*entryLYZEP)+=(Long_t)nEvtsLYZEP;
 }
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 if(entryLYZEPRP)
 {
  entryLYZEPRP->Append("M = ");
  (*entryLYZEPRP)+=(Long_t)avMultLYZEPRP;
  entryLYZEPRP->Append(", N = ");
  (*entryLYZEPRP)+=(Long_t)nEvtsLYZEPRP;
 }
 if(entryLYZEPPOI)
 {
  entryLYZEPPOI->Append("M = ");
  (*entryLYZEPPOI)+=(Long_t)avMultLYZEPPOI;
  entryLYZEPPOI->Append(", N = ");
  (*entryLYZEPPOI)+=(Long_t)nEvtsLYZEPPOI;
 }

 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //   !!!! to be removed !!!!   
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
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //RP:
 if(textResultsRP)
 {
  textResultsRP->AddText(entryMCRP->Data());
  textResultsRP->AddText(entryGFCRP->Data());
  textResultsRP->AddText(entryQC2RP->Data());
  textResultsRP->AddText(entryQC4RP->Data());
  textResultsRP->AddText(entryQC6RP->Data());
  textResultsRP->AddText(entryQC8RP->Data());
  textResultsRP->AddText(entryFQDRP->Data());
  textResultsRP->AddText(entryLYZ1RP->Data());
  textResultsRP->AddText(entryLYZEPRP->Data());
 }
 //POI:
 if(textResultsPOI)
 {
  textResultsPOI->AddText(entryMCPOI->Data());
  textResultsPOI->AddText(entryGFCPOI->Data());
  textResultsPOI->AddText(entryQC2POI->Data());
  textResultsPOI->AddText(entryQC4POI->Data());
  textResultsPOI->AddText(entryQC6POI->Data());
  textResultsPOI->AddText(entryQC8POI->Data());
  textResultsPOI->AddText(entryFQDPOI->Data());
  textResultsPOI->AddText(entryLYZ1POI->Data());
  textResultsPOI->AddText(entryLYZEPPOI->Data());
 }
 //----------------------------------------------------------------------------------
 
 
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //                           !!!! to be removed !!!!
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
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 //----------------------------------------------------------------------------------
 //final drawing for integrated flow of RP (i.e. of particles used to determine the reaction plane):
 TCanvas* intFlowAllCanvasRP = new TCanvas("Integrated Flow RP","Integrated Flow RP",1000,600);
 
 intFlowAllCanvasRP->Divide(2,1);
 
 //1st pad is for plot:
 (intFlowAllCanvasRP->cd(1))->SetPad(0.0,0.0,0.75,1.0);
 
 if(intFlowAll)
 {
  if(dMinRP>0. && dMaxRP>0.)
  {
   (intFlowAll->GetYaxis())->SetRangeUser(0.9744*dMinRP,1.0144*dMaxRP);
  } else if(dMinRP<0. && dMaxRP>0.)
    {
     if(!(-1.*dMinRP<4.*dMaxRP))
     {  
      (intFlowAll->GetYaxis())->SetRangeUser(1.0266*dMinRP,1.0144*dMaxRP);
     } else {(intFlowAll->GetYaxis())->SetRangeUser(1.1266*dMinRP,1.0144*dMaxRP);}  
    } else if(dMinRP<0. && dMaxRP<0.)
      {
       (intFlowAll->GetYaxis())->SetRangeUser(1.0266*dMinRP,0.9866*dMaxRP);      
      }
  intFlowAll->Draw("E1");
 }
                            
 if(pMeshRP) pMeshRP->Draw("LFSAME");
  
 if(flowResultsRP) flowResultsRP->Draw("PSAME");

 //2nd pad is for legend:
 (intFlowAllCanvasRP->cd(2))->SetPad(0.75,0.0,1.0,1.0);
 
 if(textDefault) textDefault->Draw();

 if(textResultsRP) textResultsRP->Draw();
 //----------------------------------------------------------------------------------
 
 //----------------------------------------------------------------------------------
 //final drawing for integrated flow of POI (i.e. of particles of interest):
 TCanvas* intFlowAllCanvasPOI = new TCanvas("Integrated Flow POI","Integrated Flow POI",1000,600);
 
 intFlowAllCanvasPOI->Divide(2,1);
 
 //1st pad is for plot:
 (intFlowAllCanvasPOI->cd(1))->SetPad(0.0,0.0,0.75,1.0);
 
 if(intFlowAll)
 {
  if(dMinPOI>0. && dMaxPOI>0.)
  {
   (intFlowAll->GetYaxis())->SetRangeUser(0.9744*dMinPOI,1.0144*dMaxPOI);
  } else if(dMinPOI<0. && dMaxPOI>0.)
    {
     if(!(-1.*dMinPOI<4.*dMaxPOI))
     {  
      (intFlowAll->GetYaxis())->SetRangeUser(1.0266*dMinPOI,1.0144*dMaxPOI);
     } else {(intFlowAll->GetYaxis())->SetRangeUser(1.1266*dMinPOI,1.0144*dMaxPOI);}  
    } else if(dMinPOI<0. && dMaxPOI<0.)
      {
       (intFlowAll->GetYaxis())->SetRangeUser(1.0266*dMinPOI,0.9866*dMaxPOI);      
      }
  intFlowAll->Draw("E1");
 }
                            
 if(pMeshPOI) pMeshPOI->Draw("LFSAME");
  
 if(flowResultsPOI) flowResultsPOI->Draw("PSAME");

 //2nd pad is for legend:
 (intFlowAllCanvasPOI->cd(2))->SetPad(0.75,0.0,1.0,1.0);
 
 if(textDefault) textDefault->Draw();

 if(textResultsPOI) textResultsPOI->Draw();
 //----------------------------------------------------------------------------------
 
 //==================================================================================   
 



 //==================================================================================
 //                            DIFFERENTIAL FLOW
 //==================================================================================
 Int_t iNbinsPt  = AliFlowCommonConstants::GetNbinsPt();
 Double_t dPtMin = AliFlowCommonConstants::GetPtMin();
 Double_t dPtMax = AliFlowCommonConstants::GetPtMax();
 
 Int_t iNbinsEta  = AliFlowCommonConstants::GetNbinsEta();
 Double_t dEtaMin = AliFlowCommonConstants::GetEtaMin();
 Double_t dEtaMax = AliFlowCommonConstants::GetEtaMax();
 
 //----------------------------------------------------------------------------------
 //cosmetics: the style histogram for differential flow (pt):
 TH1D *styleHistPt = new TH1D("styleHistPt","styleHistPt",iNbinsPt,dPtMin,dPtMax);
 styleHistPt->SetTitle("Differential Flow");
 styleHistPt->SetXTitle("p_{t} [GeV]");
 styleHistPt->SetYTitle("v_{n}");
 
 //cosmetics: the style histogram for differential flow (eta):
 TH1D *styleHistEta = new TH1D("styleHistEta","styleHistEta",iNbinsEta,dEtaMin,dEtaMax);
 styleHistEta->SetTitle("Differential Flow");
 styleHistEta->SetXTitle("#eta");
 styleHistEta->SetYTitle("v_{n}");
 //----------------------------------------------------------------------------------
 
 
 
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 //               !!!! to be removed !!!!
 //----------------------------------------------------------------------------------
 //cosmetics: Monte Carlo error bands for differential flow
 TGraph* pMeshDiffFlow = NULL;
 if(mcepCommonHistRes)
 {
  Int_t nPtsDF = (mcepCommonHistRes->GetHistDiffFlow())->GetNbinsX();
  Double_t binWidth = (mcepCommonHistRes->GetHistDiffFlow())->GetBinWidth(1);//assuming that all bins have the same width
       
  pMeshDiffFlow = new TGraph(2*nPtsDF+1);
  
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
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


 //----------------------------------------------------------------------------------
 //RP:
 //cosmetics: Monte Carlo error bands for differential flow (Pt)
 TGraph* pMeshDiffFlowPtRP = NULL;
 if(mcepCommonHistRes)
 {
  Int_t nPtsDFPtRP = (mcepCommonHistRes->GetHistDiffFlowPtRP())->GetNbinsX();
  Double_t binWidthPtRP = (mcepCommonHistRes->GetHistDiffFlowPtRP())->GetBinWidth(1);//assuming that all bins have the same width
       
  pMeshDiffFlowPtRP = new TGraph(2*nPtsDFPtRP+1);
  
  Double_t valueMCPtRP=0., errorMCPtRP=0.;
  for(Int_t i=1;i<nPtsDFPtRP+1;i++)
  {
   valueMCPtRP = (mcepCommonHistRes->GetHistDiffFlowPtRP())->GetBinContent(i);
   errorMCPtRP = (mcepCommonHistRes->GetHistDiffFlowPtRP())->GetBinError(i);       
   pMeshDiffFlowPtRP->SetPoint(i,(i-0.5)*binWidthPtRP,valueMCPtRP+errorMCPtRP);
  }    
  for(Int_t i=nPtsDFPtRP+1;i<2*nPtsDFPtRP+1;i++)
  {
   valueMCPtRP = (mcepCommonHistRes->GetHistDiffFlowPtRP())->GetBinContent(2*nPtsDFPtRP+1-i);
   errorMCPtRP = (mcepCommonHistRes->GetHistDiffFlowPtRP())->GetBinError(2*nPtsDFPtRP+1-i);       
   pMeshDiffFlowPtRP->SetPoint(i,(2*nPtsDFPtRP-i+0.5)*binWidthPtRP,valueMCPtRP-errorMCPtRP); 
  }
  pMeshDiffFlowPtRP->SetPoint(2*nPtsDFPtRP+1,0.5*binWidthPtRP,valueMCPtRP+errorMCPtRP); 
  pMeshDiffFlowPtRP->SetFillStyle(meshStyle);
  pMeshDiffFlowPtRP->SetFillColor(meshColor);
 }
 
 //cosmetics: Monte Carlo error bands for differential flow (Eta)
 TGraph* pMeshDiffFlowEtaRP = NULL;
 if(mcepCommonHistRes)
 {
  Int_t nPtsDFEtaRP = (mcepCommonHistRes->GetHistDiffFlowEtaRP())->GetNbinsX();
  Double_t binWidthEtaRP = (mcepCommonHistRes->GetHistDiffFlowEtaRP())->GetBinWidth(1);//assuming that all bins have the same width
       
  pMeshDiffFlowEtaRP = new TGraph(2*nPtsDFEtaRP+1);
  
  Double_t valueMCEtaRP=0., errorMCEtaRP=0.;
  for(Int_t i=1;i<nPtsDFEtaRP+1;i++)
  {
   valueMCEtaRP = (mcepCommonHistRes->GetHistDiffFlowEtaRP())->GetBinContent(i);
   errorMCEtaRP = (mcepCommonHistRes->GetHistDiffFlowEtaRP())->GetBinError(i);       
   pMeshDiffFlowEtaRP->SetPoint(i,(i-0.5)*binWidthEtaRP,valueMCEtaRP+errorMCEtaRP);
  }    
  for(Int_t i=nPtsDFEtaRP+1;i<2*nPtsDFEtaRP+1;i++)
  {
   valueMCEtaRP = (mcepCommonHistRes->GetHistDiffFlowEtaRP())->GetBinContent(2*nPtsDFEtaRP+1-i);
   errorMCEtaRP = (mcepCommonHistRes->GetHistDiffFlowEtaRP())->GetBinError(2*nPtsDFEtaRP+1-i);       
   pMeshDiffFlowEtaRP->SetPoint(i,(2*nPtsDFEtaRP-i+0.5)*binWidthEtaRP,valueMCEtaRP-errorMCEtaRP); 
  }
  pMeshDiffFlowEtaRP->SetPoint(2*nPtsDFEtaRP+1,0.5*binWidthEtaRP,valueMCEtaRP+errorMCEtaRP); 
  pMeshDiffFlowEtaRP->SetFillStyle(meshStyle);
  pMeshDiffFlowEtaRP->SetFillColor(meshColor);
 } 
 //----------------------------------------------------------------------------------


 //----------------------------------------------------------------------------------
 //POI:
 //cosmetics: Monte Carlo error bands for differential flow (Pt)
 TGraph* pMeshDiffFlowPtPOI = NULL;
 if(mcepCommonHistRes)
 {
  Int_t nPtsDFPtPOI = (mcepCommonHistRes->GetHistDiffFlowPtPOI())->GetNbinsX();
  Double_t binWidthPtPOI = (mcepCommonHistRes->GetHistDiffFlowPtPOI())->GetBinWidth(1);//assuming that all bins have the same width
       
  pMeshDiffFlowPtPOI = new TGraph(2*nPtsDFPtPOI+1);
  
  Double_t valueMCPtPOI=0., errorMCPtPOI=0.;
  for(Int_t i=1;i<nPtsDFPtPOI+1;i++)
  {
   valueMCPtPOI = (mcepCommonHistRes->GetHistDiffFlowPtPOI())->GetBinContent(i);
   errorMCPtPOI = (mcepCommonHistRes->GetHistDiffFlowPtPOI())->GetBinError(i);       
   pMeshDiffFlowPtPOI->SetPoint(i,(i-0.5)*binWidthPtPOI,valueMCPtPOI+errorMCPtPOI);
  }    
  for(Int_t i=nPtsDFPtPOI+1;i<2*nPtsDFPtPOI+1;i++)
  {
   valueMCPtPOI = (mcepCommonHistRes->GetHistDiffFlowPtPOI())->GetBinContent(2*nPtsDFPtPOI+1-i);
   errorMCPtPOI = (mcepCommonHistRes->GetHistDiffFlowPtPOI())->GetBinError(2*nPtsDFPtPOI+1-i);       
   pMeshDiffFlowPtPOI->SetPoint(i,(2*nPtsDFPtPOI-i+0.5)*binWidthPtPOI,valueMCPtPOI-errorMCPtPOI); 
  }
  pMeshDiffFlowPtPOI->SetPoint(2*nPtsDFPtPOI+1,0.5*binWidthPtPOI,valueMCPtPOI+errorMCPtPOI); 
  pMeshDiffFlowPtPOI->SetFillStyle(meshStyle);
  pMeshDiffFlowPtPOI->SetFillColor(meshColor);
 }
 
 //cosmetics: Monte Carlo error bands for differential flow (Eta)
 TGraph* pMeshDiffFlowEtaPOI = NULL;
 if(mcepCommonHistRes)
 {
  Int_t nPtsDFEtaPOI = (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->GetNbinsX();
  Double_t binWidthEtaPOI = (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->GetBinWidth(1);//assuming that all bins have the same width
       
  pMeshDiffFlowEtaPOI = new TGraph(2*nPtsDFEtaPOI+1);
  
  Double_t valueMCEtaPOI=0., errorMCEtaPOI=0.;
  for(Int_t i=1;i<nPtsDFEtaPOI+1;i++)
  {
   valueMCEtaPOI = (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->GetBinContent(i);
   errorMCEtaPOI = (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->GetBinError(i);       
   pMeshDiffFlowEtaPOI->SetPoint(i,(i-0.5)*binWidthEtaPOI,valueMCEtaPOI+errorMCEtaPOI);
  }    
  for(Int_t i=nPtsDFEtaPOI+1;i<2*nPtsDFEtaPOI+1;i++)
  {
   valueMCEtaPOI = (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->GetBinContent(2*nPtsDFEtaPOI+1-i);
   errorMCEtaPOI = (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->GetBinError(2*nPtsDFEtaPOI+1-i);       
   pMeshDiffFlowEtaPOI->SetPoint(i,(2*nPtsDFEtaPOI-i+0.5)*binWidthEtaPOI,valueMCEtaPOI-errorMCEtaPOI); 
  }
  pMeshDiffFlowEtaPOI->SetPoint(2*nPtsDFEtaPOI+1,0.5*binWidthEtaPOI,valueMCEtaPOI+errorMCEtaPOI); 
  pMeshDiffFlowEtaPOI->SetFillStyle(meshStyle);
  pMeshDiffFlowEtaPOI->SetFillColor(meshColor);
 } 
 //----------------------------------------------------------------------------------
   
 //MCEP = Monte Carlo Event Plane
 Double_t avMultDiffFlowMC=0.;//to be removed
 Double_t nEvtsDiffFlowMC=0;//to be removed
 Double_t avMultDiffFlowMCRP=0.;
 Double_t nEvtsDiffFlowMCRP=0;
 Double_t avMultDiffFlowMCPOI=0.;
 Double_t nEvtsDiffFlowMCPOI=0;
 if(fileMCEP)
 {
  if(mcepCommonHistRes)
  {
   (mcepCommonHistRes->GetHistDiffFlow())->SetMarkerColor(2);//to be removed
   (mcepCommonHistRes->GetHistDiffFlow())->SetMarkerStyle(20);//to be removed
   (mcepCommonHistRes->GetHistDiffFlowPtRP())->SetMarkerColor(2);
   (mcepCommonHistRes->GetHistDiffFlowPtRP())->SetMarkerStyle(20);
   (mcepCommonHistRes->GetHistDiffFlowEtaRP())->SetMarkerColor(2);
   (mcepCommonHistRes->GetHistDiffFlowEtaRP())->SetMarkerStyle(20);
   (mcepCommonHistRes->GetHistDiffFlowPtPOI())->SetMarkerColor(2);
   (mcepCommonHistRes->GetHistDiffFlowPtPOI())->SetMarkerStyle(20);
   (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->SetMarkerColor(2);
   (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->SetMarkerStyle(20);
  } 
  if(mcepCommonHist)
  {
   avMultDiffFlowMC = (mcepCommonHist->GetHistMultDiff())->GetMean();//to be removed
   nEvtsDiffFlowMC  = (mcepCommonHist->GetHistMultDiff())->GetEntries();//to be removed
   avMultDiffFlowMCRP = (mcepCommonHist->GetHistMultInt())->GetMean();
   nEvtsDiffFlowMCRP  = (mcepCommonHist->GetHistMultInt())->GetEntries();
   avMultDiffFlowMCPOI = (mcepCommonHist->GetHistMultDiff())->GetMean();
   nEvtsDiffFlowMCPOI  = (mcepCommonHist->GetHistMultDiff())->GetEntries();      
  } 
 } 

 //GFC = Generating Function Cumulants
 Double_t avMultDiffFlowGFC=0.;//to be removed
 Double_t nEvtsDiffFlowGFC=0.;//to be removed
 Double_t avMultDiffFlowGFCRP=0.;
 Double_t nEvtsDiffFlowGFCRP=0.;
 Double_t avMultDiffFlowGFCPOI=0.;
 Double_t nEvtsDiffFlowGFCPOI=0.; 
 if(fileGFC)
 {
  if(gfcCommonHistRes2)
  {
   (gfcCommonHistRes2->GetHistDiffFlow())->SetMarkerColor(kViolet+3);//to be removed
   (gfcCommonHistRes2->GetHistDiffFlow())->SetMarkerStyle(20);//to be removed
   (gfcCommonHistRes2->GetHistDiffFlowPtRP())->SetMarkerColor(kViolet+3);
   (gfcCommonHistRes2->GetHistDiffFlowPtRP())->SetMarkerStyle(20);
   (gfcCommonHistRes2->GetHistDiffFlowEtaRP())->SetMarkerColor(kViolet+3);
   (gfcCommonHistRes2->GetHistDiffFlowEtaRP())->SetMarkerStyle(20);
   (gfcCommonHistRes2->GetHistDiffFlowPtPOI())->SetMarkerColor(kViolet+3);
   (gfcCommonHistRes2->GetHistDiffFlowPtPOI())->SetMarkerStyle(20);
   (gfcCommonHistRes2->GetHistDiffFlowEtaPOI())->SetMarkerColor(kViolet+3);
   (gfcCommonHistRes2->GetHistDiffFlowEtaPOI())->SetMarkerStyle(20);
  }
  if(gfcCommonHistRes4)
  { 
   (gfcCommonHistRes4->GetHistDiffFlow())->SetMarkerColor(kViolet-6);//to be removed
   (gfcCommonHistRes4->GetHistDiffFlow())->SetMarkerStyle(21);//to be removed
   (gfcCommonHistRes4->GetHistDiffFlowPtRP())->SetMarkerColor(kViolet-6);
   (gfcCommonHistRes4->GetHistDiffFlowPtRP())->SetMarkerStyle(21);
   (gfcCommonHistRes4->GetHistDiffFlowEtaRP())->SetMarkerColor(kViolet-6);
   (gfcCommonHistRes4->GetHistDiffFlowEtaRP())->SetMarkerStyle(21);
   (gfcCommonHistRes4->GetHistDiffFlowPtPOI())->SetMarkerColor(kViolet-6);
   (gfcCommonHistRes4->GetHistDiffFlowPtPOI())->SetMarkerStyle(21);
   (gfcCommonHistRes4->GetHistDiffFlowEtaPOI())->SetMarkerColor(kViolet-6);
   (gfcCommonHistRes4->GetHistDiffFlowEtaPOI())->SetMarkerStyle(21);         
  }
  if(gfcCommonHistRes6)
  { 
   (gfcCommonHistRes6->GetHistDiffFlow())->SetMarkerColor(kViolet+5);//to be removed
   (gfcCommonHistRes6->GetHistDiffFlow())->SetMarkerStyle(24);//to be removed
   (gfcCommonHistRes6->GetHistDiffFlowPtRP())->SetMarkerColor(kViolet+5);
   (gfcCommonHistRes6->GetHistDiffFlowPtRP())->SetMarkerStyle(24);
   (gfcCommonHistRes6->GetHistDiffFlowEtaRP())->SetMarkerColor(kViolet+5);
   (gfcCommonHistRes6->GetHistDiffFlowEtaRP())->SetMarkerStyle(24);
   (gfcCommonHistRes6->GetHistDiffFlowPtPOI())->SetMarkerColor(kViolet+5);
   (gfcCommonHistRes6->GetHistDiffFlowPtPOI())->SetMarkerStyle(24);
   (gfcCommonHistRes6->GetHistDiffFlowEtaPOI())->SetMarkerColor(kViolet+5);
   (gfcCommonHistRes6->GetHistDiffFlowEtaPOI())->SetMarkerStyle(24);
  }
  if(gfcCommonHistRes8)
  { 
   (gfcCommonHistRes8->GetHistDiffFlow())->SetMarkerColor(kViolet-5);//to be removed
   (gfcCommonHistRes8->GetHistDiffFlow())->SetMarkerStyle(25);//to be removed
   (gfcCommonHistRes8->GetHistDiffFlowPtRP())->SetMarkerColor(kViolet-5);
   (gfcCommonHistRes8->GetHistDiffFlowPtRP())->SetMarkerStyle(25);
   (gfcCommonHistRes8->GetHistDiffFlowEtaRP())->SetMarkerColor(kViolet-5);
   (gfcCommonHistRes8->GetHistDiffFlowEtaRP())->SetMarkerStyle(25);
   (gfcCommonHistRes8->GetHistDiffFlowPtPOI())->SetMarkerColor(kViolet-5);
   (gfcCommonHistRes8->GetHistDiffFlowPtPOI())->SetMarkerStyle(25);
   (gfcCommonHistRes8->GetHistDiffFlowEtaPOI())->SetMarkerColor(kViolet-5);
   (gfcCommonHistRes8->GetHistDiffFlowEtaPOI())->SetMarkerStyle(25);
  }
  if(gfcCommonHist)
  {
   avMultDiffFlowGFC = (gfcCommonHist->GetHistMultDiff())->GetMean();//to be removed
   nEvtsDiffFlowGFC  = (gfcCommonHist->GetHistMultDiff())->GetEntries();//to be removed
   avMultDiffFlowGFCRP = (gfcCommonHist->GetHistMultInt())->GetMean();   
   nEvtsDiffFlowGFCRP  = (gfcCommonHist->GetHistMultInt())->GetEntries();
   avMultDiffFlowGFCPOI = (gfcCommonHist->GetHistMultDiff())->GetMean();
   nEvtsDiffFlowGFCPOI  = (gfcCommonHist->GetHistMultDiff())->GetEntries();   
  } 
 }
  
 //QC = Q-cumulants
 Double_t avMultDiffFlowQC2=0., avMultDiffFlowQC4=0.;//to be removed
 Double_t nEvtsDiffFlowQC2=0., nEvtsDiffFlowQC4=0.;//to be removed
 Double_t avMultDiffFlowQC2RP=0.;
 Double_t nEvtsDiffFlowQC2RP=0.;
 Double_t avMultDiffFlowQC2POI=0.;
 Double_t nEvtsDiffFlowQC2POI=0.;
 Double_t avMultDiffFlowQC4RP=0.;
 Double_t nEvtsDiffFlowQC4RP=0.;
 Double_t avMultDiffFlowQC4POI=0.;
 Double_t nEvtsDiffFlowQC4POI=0.;

 if(fileQC)
 {
  //QC{2}
  if(qcCommonHistRes2)
  {
   (qcCommonHistRes2->GetHistDiffFlow())->SetMarkerColor(kOrange+3);//to be removed
   (qcCommonHistRes2->GetHistDiffFlow())->SetMarkerStyle(20);//to be removed
   (qcCommonHistRes2->GetHistDiffFlowPtRP())->SetMarkerColor(kOrange+3);
   (qcCommonHistRes2->GetHistDiffFlowPtRP())->SetMarkerStyle(20);
   (qcCommonHistRes2->GetHistDiffFlowEtaRP())->SetMarkerColor(kOrange+3);
   (qcCommonHistRes2->GetHistDiffFlowEtaRP())->SetMarkerStyle(20);
   (qcCommonHistRes2->GetHistDiffFlowPtPOI())->SetMarkerColor(kOrange+3);
   (qcCommonHistRes2->GetHistDiffFlowPtPOI())->SetMarkerStyle(20);
   (qcCommonHistRes2->GetHistDiffFlowEtaPOI())->SetMarkerColor(kOrange+3);
   (qcCommonHistRes2->GetHistDiffFlowEtaPOI())->SetMarkerStyle(20);
  }
  if(qcCommonHist2)
  {
   avMultDiffFlowQC2 = (qcCommonHist2->GetHistMultDiff())->GetMean();//to be removed
   nEvtsDiffFlowQC2  = (qcCommonHist2->GetHistMultDiff())->GetEntries();//to be removed
   avMultDiffFlowQC2RP = (qcCommonHist2->GetHistMultInt())->GetMean();
   nEvtsDiffFlowQC2RP  = (qcCommonHist2->GetHistMultInt())->GetEntries();
   avMultDiffFlowQC2POI = (qcCommonHist2->GetHistMultDiff())->GetMean();
   nEvtsDiffFlowQC2POI  = (qcCommonHist2->GetHistMultDiff())->GetEntries();
  }
  //QC{4}
  if(qcCommonHistRes4)
  {
   (qcCommonHistRes4->GetHistDiffFlow())->SetMarkerColor(kOrange-6);//to be removed
   (qcCommonHistRes4->GetHistDiffFlow())->SetMarkerStyle(21);//to be removed
   (qcCommonHistRes4->GetHistDiffFlowPtRP())->SetMarkerColor(kOrange-6);
   (qcCommonHistRes4->GetHistDiffFlowPtRP())->SetMarkerStyle(21);
   (qcCommonHistRes4->GetHistDiffFlowEtaRP())->SetMarkerColor(kOrange-6);
   (qcCommonHistRes4->GetHistDiffFlowEtaRP())->SetMarkerStyle(21);
   (qcCommonHistRes4->GetHistDiffFlowPtPOI())->SetMarkerColor(kOrange-6);
   (qcCommonHistRes4->GetHistDiffFlowPtPOI())->SetMarkerStyle(21);
   (qcCommonHistRes4->GetHistDiffFlowEtaPOI())->SetMarkerColor(kOrange-6);
   (qcCommonHistRes4->GetHistDiffFlowEtaPOI())->SetMarkerStyle(21);
  }
  if(qcCommonHist4)
  {
   avMultDiffFlowQC4 = (qcCommonHist4->GetHistMultDiff())->GetMean();//to be removed
   nEvtsDiffFlowQC4  = (qcCommonHist4->GetHistMultDiff())->GetEntries();//to be removed
   avMultDiffFlowQC4RP = (qcCommonHist4->GetHistMultInt())->GetMean();
   nEvtsDiffFlowQC4RP  = (qcCommonHist4->GetHistMultInt())->GetEntries();
   avMultDiffFlowQC4POI = (qcCommonHist4->GetHistMultDiff())->GetMean();
   nEvtsDiffFlowQC4POI  = (qcCommonHist4->GetHistMultDiff())->GetEntries();
  }
 } 

 //LYZ2 = Lee-Yang Zeros (2nd run)
 Double_t avMultDiffFlowLYZ2=0.;//to be removed
 Double_t nEvtsDiffFlowLYZ2=0;//to be removed
 Double_t avMultDiffFlowLYZ2RP=0.;
 Double_t nEvtsDiffFlowLYZ2RP=0;
 Double_t avMultDiffFlowLYZ2POI=0.;
 Double_t nEvtsDiffFlowLYZ2POI=0;
 if(fileLYZ2)
 {
  if(lyz2CommonHistRes)
  {
   (lyz2CommonHistRes->GetHistDiffFlow())->SetMarkerColor(kGreen+3);//to be removed
   (lyz2CommonHistRes->GetHistDiffFlow())->SetMarkerStyle(22);//to be removed
   (lyz2CommonHistRes->GetHistDiffFlowPtRP())->SetMarkerColor(kGreen+3);
   (lyz2CommonHistRes->GetHistDiffFlowPtRP())->SetMarkerStyle(22);
   (lyz2CommonHistRes->GetHistDiffFlowEtaRP())->SetMarkerColor(kGreen+3);
   (lyz2CommonHistRes->GetHistDiffFlowEtaRP())->SetMarkerStyle(22);
   (lyz2CommonHistRes->GetHistDiffFlowPtPOI())->SetMarkerColor(kGreen+3);
   (lyz2CommonHistRes->GetHistDiffFlowPtPOI())->SetMarkerStyle(22);
   (lyz2CommonHistRes->GetHistDiffFlowEtaPOI())->SetMarkerColor(kGreen+3);
   (lyz2CommonHistRes->GetHistDiffFlowEtaPOI())->SetMarkerStyle(22);
  } 
  if(lyz2CommonHist)
  {
   avMultDiffFlowLYZ2 = (lyz2CommonHist->GetHistMultDiff())->GetMean();//to be removed
   nEvtsDiffFlowLYZ2  = (lyz2CommonHist->GetHistMultDiff())->GetEntries();//to be removed
   avMultDiffFlowLYZ2RP = (lyz2CommonHist->GetHistMultInt())->GetMean();
   nEvtsDiffFlowLYZ2RP  = (lyz2CommonHist->GetHistMultInt())->GetEntries();
   avMultDiffFlowLYZ2POI = (lyz2CommonHist->GetHistMultDiff())->GetMean();
   nEvtsDiffFlowLYZ2POI  = (lyz2CommonHist->GetHistMultDiff())->GetEntries();
  } 
 } 

 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 // !!!! to be removed !!!!
 //----------------------------------------------------------------------------------
 //final drawing for differential flow:
 TCanvas* diffFlowAllCanvas = new TCanvas("Differential Flow","Differential Flow",1000,600);
 
 diffFlowAllCanvas->Divide(2,1);
 
 //1st pad is for plot:
 (diffFlowAllCanvas->cd(1))->SetPad(0.0,0.0,0.75,1.0);
  
 if(styleHistPt)
 {
  styleHistPt->Draw();
 }

 if(pMeshDiffFlow)
 {
  pMeshDiffFlow->Draw("LFSAME");
 }
 //MC 
 if(mcepCommonHistRes)
 { 
  (mcepCommonHistRes->GetHistDiffFlow())->Draw("E1PSAME");
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
 if(gfcCommonHistRes6)
 { 
  (gfcCommonHistRes6->GetHistDiffFlow())->Draw("E1PSAME"); 
 } 
 if(gfcCommonHistRes8)
 { 
  (gfcCommonHistRes8->GetHistDiffFlow())->Draw("E1PSAME"); 
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
 TString *entryDiffGFC6 = new TString("GFC{6} ... ");
 TString *entryDiffGFC8 = new TString("GFC{8} ... "); 
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
 if(gfcCommonHistRes6)
 {
  entryDiffGFC6->Append("M = ");
  (*entryDiffGFC6)+=(Long_t)avMultDiffFlowGFC;
  entryDiffGFC6->Append(", N = ");
  (*entryDiffGFC6)+=(Long_t)nEvtsDiffFlowGFC; 
  legendDiffFlow->AddEntry(gfcCommonHistRes6->GetHistDiffFlow(),entryDiffGFC6->Data(),"p");
 } 
 if(gfcCommonHistRes8)
 {
  entryDiffGFC8->Append("M = ");
  (*entryDiffGFC8)+=(Long_t)avMultDiffFlowGFC;
  entryDiffGFC8->Append(", N = ");
  (*entryDiffGFC8)+=(Long_t)nEvtsDiffFlowGFC; 
  legendDiffFlow->AddEntry(gfcCommonHistRes8->GetHistDiffFlow(),entryDiffGFC8->Data(),"p");
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
  legendDiffFlow->SetMargin(0.15);
  legendDiffFlow->Draw();
 }
 //----------------------------------------------------------------------------------
 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 
 //----------------------------------------------------------------------------------
 //final drawing for differential flow (Pt, RP):
 TCanvas* diffFlowPtAllCanvasRP = new TCanvas("Differential Flow (Pt) of RP","Differential Flow (Pt) of RP ",1000,600);
 
 diffFlowPtAllCanvasRP->Divide(2,1);
 
 //1st pad is for plot:
 (diffFlowPtAllCanvasRP->cd(1))->SetPad(0.0,0.0,0.75,1.0);
 
 if(styleHistPt)
 {
  styleHistPt->Draw();
 }
 if(pMeshDiffFlowPtRP)
 {
  pMeshDiffFlowPtRP->Draw("LFSAME");
 }
 
 //MC 
 if(mcepCommonHistRes)
 { 
  (mcepCommonHistRes->GetHistDiffFlowPtRP())->Draw("E1PSAME");
 }
 //GFC
 if(gfcCommonHistRes2)
 { 
  (gfcCommonHistRes2->GetHistDiffFlowPtRP())->Draw("E1PSAME"); 
 } 
 if(gfcCommonHistRes4)
 { 
  (gfcCommonHistRes4->GetHistDiffFlowPtRP())->Draw("E1PSAME"); 
 } 
 if(gfcCommonHistRes6)
 { 
  (gfcCommonHistRes6->GetHistDiffFlowPtRP())->Draw("E1PSAME"); 
 } 
 if(gfcCommonHistRes8)
 { 
  (gfcCommonHistRes8->GetHistDiffFlowPtRP())->Draw("E1PSAME"); 
 }    
 //QC
 if(qcCommonHistRes2)
 { 
  (qcCommonHistRes2->GetHistDiffFlowPtRP())->Draw("E1PSAME");
 }
 if(qcCommonHistRes4)
 { 
  (qcCommonHistRes4->GetHistDiffFlowPtRP())->Draw("E1PSAME");
 }
 //LYZ2
 if(lyz2CommonHistRes)
 { 
  (lyz2CommonHistRes->GetHistDiffFlowPtRP())->Draw("E1PSAME");
 }
 
 //2nd pad is for legend:
 (diffFlowPtAllCanvasRP->cd(2))->SetPad(0.75,0.0,1.0,1.0);

 TLegend* legendDiffFlowPtRP = new TLegend(0.02,0.25,0.97,0.75);
 legendDiffFlowPtRP->SetTextFont(72);
 legendDiffFlowPtRP->SetTextSize(0.06);
 
 //legend's entries:
 TString *entryDiffMCPtRP   = new TString("MC ....... ");
 TString *entryDiffGFC2PtRP = new TString("GFC{2} ... ");
 TString *entryDiffGFC4PtRP = new TString("GFC{4} ... ");
 TString *entryDiffGFC6PtRP = new TString("GFC{6} ... ");
 TString *entryDiffGFC8PtRP = new TString("GFC{8} ... "); 
 TString *entryDiffQC2PtRP  = new TString("QC{2} .... ");
 TString *entryDiffQC4PtRP  = new TString("QC{4} .... ");
 TString *entryDiffLYZ2PtRP = new TString("LYZ ...... ");
 
 //MC
 if(mcepCommonHistRes)
 {
  (mcepCommonHistRes->GetHistDiffFlow())->SetFillStyle(meshStyle);
  (mcepCommonHistRes->GetHistDiffFlow())->SetFillColor(meshColor);
  entryDiffMCPtRP->Append("M = ");
  (*entryDiffMCPtRP)+=(Long_t)avMultDiffFlowMCRP;
  entryDiffMCPtRP->Append(", N = ");
  (*entryDiffMCPtRP)+=(Long_t)nEvtsDiffFlowMCRP; 
  legendDiffFlowPtRP->AddEntry(mcepCommonHistRes->GetHistDiffFlowPtRP(),entryDiffMCPtRP->Data(),"f");
 }

 //GFC
 if(gfcCommonHistRes2)
 {
  entryDiffGFC2PtRP->Append("M = ");
  (*entryDiffGFC2PtRP)+=(Long_t)avMultDiffFlowGFCRP;
  entryDiffGFC2PtRP->Append(", N = ");
  (*entryDiffGFC2PtRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
  legendDiffFlowPtRP->AddEntry(gfcCommonHistRes2->GetHistDiffFlowPtRP(),entryDiffGFC2PtRP->Data(),"p");
 }
 if(gfcCommonHistRes4)
 {
  entryDiffGFC4PtRP->Append("M = ");
  (*entryDiffGFC4PtRP)+=(Long_t)avMultDiffFlowGFCRP;
  entryDiffGFC4PtRP->Append(", N = ");
  (*entryDiffGFC4PtRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
  legendDiffFlowPtRP->AddEntry(gfcCommonHistRes4->GetHistDiffFlowPtRP(),entryDiffGFC4PtRP->Data(),"p");
 }
 if(gfcCommonHistRes6)
 {
  entryDiffGFC6PtRP->Append("M = ");
  (*entryDiffGFC6PtRP)+=(Long_t)avMultDiffFlowGFCRP;
  entryDiffGFC6PtRP->Append(", N = ");
  (*entryDiffGFC6PtRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
  legendDiffFlowPtRP->AddEntry(gfcCommonHistRes6->GetHistDiffFlowPtRP(),entryDiffGFC6PtRP->Data(),"p");
 } 
 if(gfcCommonHistRes8)
 {
  entryDiffGFC8PtRP->Append("M = ");
  (*entryDiffGFC8PtRP)+=(Long_t)avMultDiffFlowGFCRP;
  entryDiffGFC8PtRP->Append(", N = ");
  (*entryDiffGFC8PtRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
  legendDiffFlowPtRP->AddEntry(gfcCommonHistRes8->GetHistDiffFlowPtRP(),entryDiffGFC8PtRP->Data(),"p");
 }  
 
 //QC
 if(qcCommonHistRes2)
 {
  entryDiffQC2PtRP->Append("M = ");
  (*entryDiffQC2PtRP)+=(Long_t)avMultDiffFlowQC2RP;
  entryDiffQC2PtRP->Append(", N = ");
  (*entryDiffQC2PtRP)+=(Long_t)nEvtsDiffFlowQC2RP; 
  legendDiffFlowPtRP->AddEntry(qcCommonHistRes2->GetHistDiffFlowPtRP(),entryDiffQC2PtRP->Data(),"p");
 }
 if(qcCommonHistRes4)
 {
  entryDiffQC4PtRP->Append("M = ");
  (*entryDiffQC4PtRP)+=(Long_t)avMultDiffFlowQC4RP;
  entryDiffQC4PtRP->Append(", N = ");
  (*entryDiffQC4PtRP)+=(Long_t)nEvtsDiffFlowQC4RP; 
  legendDiffFlowPtRP->AddEntry(qcCommonHistRes4->GetHistDiffFlowPtRP(),entryDiffQC4PtRP->Data(),"p");
 }
 
 //LYZ
 if(lyz2CommonHistRes)
 {
  entryDiffLYZ2PtRP->Append("M = ");
  (*entryDiffLYZ2PtRP)+=(Long_t)avMultDiffFlowLYZ2RP;
  entryDiffLYZ2PtRP->Append(", N = ");
  (*entryDiffLYZ2PtRP)+=(Long_t)nEvtsDiffFlowLYZ2RP; 
  legendDiffFlowPtRP->AddEntry(lyz2CommonHistRes->GetHistDiffFlowPtRP(),entryDiffLYZ2PtRP->Data(),"p");
 }

 //drawing finally the legend in the 2nd pad:     
 if(legendDiffFlowPtRP)
 {
  legendDiffFlowPtRP->SetMargin(0.15);
  legendDiffFlowPtRP->Draw();
 }

 //----------------------------------------------------------------------------------
 
 
 //----------------------------------------------------------------------------------
 //final drawing for differential flow (Eta, RP):
 TCanvas* diffFlowEtaAllCanvasRP = new TCanvas("Differential Flow (Eta) of RP","Differential Flow (Eta) of RP ",1000,600);
 
 diffFlowEtaAllCanvasRP->Divide(2,1);
 
 //1st pad is for plot:
 (diffFlowEtaAllCanvasRP->cd(1))->SetPad(0.0,0.0,0.75,1.0);
 
 if(styleHistEta)
 {
  styleHistEta->Draw();
 }
 if(pMeshDiffFlowEtaRP)
 {
  pMeshDiffFlowEtaRP->Draw("LFSAME");
 }
 
 //MC 
 if(mcepCommonHistRes)
 { 
  (mcepCommonHistRes->GetHistDiffFlowEtaRP())->Draw("E1PSAME");
 }
 //GFC
 if(gfcCommonHistRes2)
 { 
  (gfcCommonHistRes2->GetHistDiffFlowEtaRP())->Draw("E1PSAME"); 
 } 
 if(gfcCommonHistRes4)
 { 
  (gfcCommonHistRes4->GetHistDiffFlowEtaRP())->Draw("E1PSAME"); 
 } 
 if(gfcCommonHistRes6)
 { 
  (gfcCommonHistRes6->GetHistDiffFlowEtaRP())->Draw("E1PSAME"); 
 } 
 if(gfcCommonHistRes8)
 { 
  (gfcCommonHistRes8->GetHistDiffFlowEtaRP())->Draw("E1PSAME"); 
 }    
 //QC
 if(qcCommonHistRes2)
 { 
  (qcCommonHistRes2->GetHistDiffFlowEtaRP())->Draw("E1PSAME");
 }
 if(qcCommonHistRes4)
 { 
  (qcCommonHistRes4->GetHistDiffFlowEtaRP())->Draw("E1PSAME");
 }
 //LYZ2
 if(lyz2CommonHistRes)
 { 
  (lyz2CommonHistRes->GetHistDiffFlowEtaRP())->Draw("E1PSAME");
 }
 
 //2nd pad is for legend:
 (diffFlowEtaAllCanvasRP->cd(2))->SetPad(0.75,0.0,1.0,1.0);

 TLegend* legendDiffFlowEtaRP = new TLegend(0.02,0.25,0.97,0.75);
 legendDiffFlowEtaRP->SetTextFont(72);
 legendDiffFlowEtaRP->SetTextSize(0.06);
 
 //legend's entries:
 TString *entryDiffMCEtaRP   = new TString("MC ....... ");
 TString *entryDiffGFC2EtaRP = new TString("GFC{2} ... ");
 TString *entryDiffGFC4EtaRP = new TString("GFC{4} ... ");
 TString *entryDiffGFC6EtaRP = new TString("GFC{6} ... ");
 TString *entryDiffGFC8EtaRP = new TString("GFC{8} ... "); 
 TString *entryDiffQC2EtaRP  = new TString("QC{2} .... ");
 TString *entryDiffQC4EtaRP  = new TString("QC{4} .... ");
 TString *entryDiffLYZ2EtaRP = new TString("LYZ ...... ");
 
 //MC
 if(mcepCommonHistRes)
 {
  (mcepCommonHistRes->GetHistDiffFlow())->SetFillStyle(meshStyle);
  (mcepCommonHistRes->GetHistDiffFlow())->SetFillColor(meshColor);
  entryDiffMCEtaRP->Append("M = ");
  (*entryDiffMCEtaRP)+=(Long_t)avMultDiffFlowMCRP;
  entryDiffMCEtaRP->Append(", N = ");
  (*entryDiffMCEtaRP)+=(Long_t)nEvtsDiffFlowMCRP; 
  legendDiffFlowEtaRP->AddEntry(mcepCommonHistRes->GetHistDiffFlowEtaRP(),entryDiffMCEtaRP->Data(),"f");
 }

 //GFC
 if(gfcCommonHistRes2)
 {
  entryDiffGFC2EtaRP->Append("M = ");
  (*entryDiffGFC2EtaRP)+=(Long_t)avMultDiffFlowGFCRP;
  entryDiffGFC2EtaRP->Append(", N = ");
  (*entryDiffGFC2EtaRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
  legendDiffFlowEtaRP->AddEntry(gfcCommonHistRes2->GetHistDiffFlowEtaRP(),entryDiffGFC2EtaRP->Data(),"p");
 }
 if(gfcCommonHistRes4)
 {
  entryDiffGFC4EtaRP->Append("M = ");
  (*entryDiffGFC4EtaRP)+=(Long_t)avMultDiffFlowGFCRP;
  entryDiffGFC4EtaRP->Append(", N = ");
  (*entryDiffGFC4EtaRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
  legendDiffFlowEtaRP->AddEntry(gfcCommonHistRes4->GetHistDiffFlowEtaRP(),entryDiffGFC4EtaRP->Data(),"p");
 }
 if(gfcCommonHistRes6)
 {
  entryDiffGFC6EtaRP->Append("M = ");
  (*entryDiffGFC6EtaRP)+=(Long_t)avMultDiffFlowGFCRP;
  entryDiffGFC6EtaRP->Append(", N = ");
  (*entryDiffGFC6EtaRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
  legendDiffFlowEtaRP->AddEntry(gfcCommonHistRes6->GetHistDiffFlowEtaRP(),entryDiffGFC6EtaRP->Data(),"p");
 } 
 if(gfcCommonHistRes8)
 {
  entryDiffGFC8EtaRP->Append("M = ");
  (*entryDiffGFC8EtaRP)+=(Long_t)avMultDiffFlowGFCRP;
  entryDiffGFC8EtaRP->Append(", N = ");
  (*entryDiffGFC8EtaRP)+=(Long_t)nEvtsDiffFlowGFCRP; 
  legendDiffFlowEtaRP->AddEntry(gfcCommonHistRes8->GetHistDiffFlowEtaRP(),entryDiffGFC8EtaRP->Data(),"p");
 }  
 
 //QC
 if(qcCommonHistRes2)
 {
  entryDiffQC2EtaRP->Append("M = ");
  (*entryDiffQC2EtaRP)+=(Long_t)avMultDiffFlowQC2RP;
  entryDiffQC2EtaRP->Append(", N = ");
  (*entryDiffQC2EtaRP)+=(Long_t)nEvtsDiffFlowQC2RP; 
  legendDiffFlowEtaRP->AddEntry(qcCommonHistRes2->GetHistDiffFlowEtaRP(),entryDiffQC2EtaRP->Data(),"p");
 }
 if(qcCommonHistRes4)
 {
  entryDiffQC4EtaRP->Append("M = ");
  (*entryDiffQC4EtaRP)+=(Long_t)avMultDiffFlowQC4RP;
  entryDiffQC4EtaRP->Append(", N = ");
  (*entryDiffQC4EtaRP)+=(Long_t)nEvtsDiffFlowQC4RP; 
  legendDiffFlowEtaRP->AddEntry(qcCommonHistRes4->GetHistDiffFlowEtaRP(),entryDiffQC4EtaRP->Data(),"p");
 }
 
 //LYZ
 if(lyz2CommonHistRes)
 {
  entryDiffLYZ2EtaRP->Append("M = ");
  (*entryDiffLYZ2EtaRP)+=(Long_t)avMultDiffFlowLYZ2RP;
  entryDiffLYZ2EtaRP->Append(", N = ");
  (*entryDiffLYZ2EtaRP)+=(Long_t)nEvtsDiffFlowLYZ2RP; 
  legendDiffFlowEtaRP->AddEntry(lyz2CommonHistRes->GetHistDiffFlowEtaRP(),entryDiffLYZ2EtaRP->Data(),"p");
 }

 //drawing finally the legend in the 2nd pad:     
 if(legendDiffFlowEtaRP)
 {
  legendDiffFlowEtaRP->SetMargin(0.15);
  legendDiffFlowEtaRP->Draw();
 }

 //----------------------------------------------------------------------------------

 //----------------------------------------------------------------------------------
 //final drawing for differential flow (Pt, POI):
 TCanvas* diffFlowPtAllCanvasPOI = new TCanvas("Differential Flow (Pt) of POI","Differential Flow (Pt) of POI ",1000,600);
 
 diffFlowPtAllCanvasPOI->Divide(2,1);
 
 //1st pad is for plot:
 (diffFlowPtAllCanvasPOI->cd(1))->SetPad(0.0,0.0,0.75,1.0);
 
 if(styleHistPt)
 {
  styleHistPt->Draw();
 }
 if(pMeshDiffFlowPtPOI)
 {
  pMeshDiffFlowPtPOI->Draw("LFSAME");
 }
 
 //MC 
 if(mcepCommonHistRes)
 { 
  (mcepCommonHistRes->GetHistDiffFlowPtPOI())->Draw("E1PSAME");
 }
 //GFC
 if(gfcCommonHistRes2)
 { 
  (gfcCommonHistRes2->GetHistDiffFlowPtPOI())->Draw("E1PSAME"); 
 } 
 if(gfcCommonHistRes4)
 { 
  (gfcCommonHistRes4->GetHistDiffFlowPtPOI())->Draw("E1PSAME"); 
 } 
 if(gfcCommonHistRes6)
 { 
  (gfcCommonHistRes6->GetHistDiffFlowPtPOI())->Draw("E1PSAME"); 
 } 
 if(gfcCommonHistRes8)
 { 
  (gfcCommonHistRes8->GetHistDiffFlowPtPOI())->Draw("E1PSAME"); 
 }    
 //QC
 if(qcCommonHistRes2)
 { 
  (qcCommonHistRes2->GetHistDiffFlowPtPOI())->Draw("E1PSAME");
 }
 if(qcCommonHistRes4)
 { 
  (qcCommonHistRes4->GetHistDiffFlowPtPOI())->Draw("E1PSAME");
 }
 //LYZ2
 if(lyz2CommonHistRes)
 { 
  (lyz2CommonHistRes->GetHistDiffFlowPtPOI())->Draw("E1PSAME");
 }
 
 //2nd pad is for legend:
 (diffFlowPtAllCanvasPOI->cd(2))->SetPad(0.75,0.0,1.0,1.0);

 TLegend* legendDiffFlowPtPOI = new TLegend(0.02,0.25,0.97,0.75);
 legendDiffFlowPtPOI->SetTextFont(72);
 legendDiffFlowPtPOI->SetTextSize(0.06);
 
 //legend's entries:
 TString *entryDiffMCPtPOI   = new TString("MC ....... ");
 TString *entryDiffGFC2PtPOI = new TString("GFC{2} ... ");
 TString *entryDiffGFC4PtPOI = new TString("GFC{4} ... ");
 TString *entryDiffGFC6PtPOI = new TString("GFC{6} ... ");
 TString *entryDiffGFC8PtPOI = new TString("GFC{8} ... "); 
 TString *entryDiffQC2PtPOI  = new TString("QC{2} .... ");
 TString *entryDiffQC4PtPOI  = new TString("QC{4} .... ");
 TString *entryDiffLYZ2PtPOI = new TString("LYZ ...... ");
 
 //MC
 if(mcepCommonHistRes)
 {
  (mcepCommonHistRes->GetHistDiffFlow())->SetFillStyle(meshStyle);
  (mcepCommonHistRes->GetHistDiffFlow())->SetFillColor(meshColor);
  entryDiffMCPtPOI->Append("M = ");
  (*entryDiffMCPtPOI)+=(Long_t)avMultDiffFlowMCPOI;
  entryDiffMCPtPOI->Append(", N = ");
  (*entryDiffMCPtPOI)+=(Long_t)nEvtsDiffFlowMCPOI; 
  legendDiffFlowPtPOI->AddEntry(mcepCommonHistRes->GetHistDiffFlowPtPOI(),entryDiffMCPtPOI->Data(),"f");
 }

 //GFC
 if(gfcCommonHistRes2)
 {
  entryDiffGFC2PtPOI->Append("M = ");
  (*entryDiffGFC2PtPOI)+=(Long_t)avMultDiffFlowGFCPOI;
  entryDiffGFC2PtPOI->Append(", N = ");
  (*entryDiffGFC2PtPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
  legendDiffFlowPtPOI->AddEntry(gfcCommonHistRes2->GetHistDiffFlowPtPOI(),entryDiffGFC2PtPOI->Data(),"p");
 }
 if(gfcCommonHistRes4)
 {
  entryDiffGFC4PtPOI->Append("M = ");
  (*entryDiffGFC4PtPOI)+=(Long_t)avMultDiffFlowGFCPOI;
  entryDiffGFC4PtPOI->Append(", N = ");
  (*entryDiffGFC4PtPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
  legendDiffFlowPtPOI->AddEntry(gfcCommonHistRes4->GetHistDiffFlowPtPOI(),entryDiffGFC4PtPOI->Data(),"p");
 }
 if(gfcCommonHistRes6)
 {
  entryDiffGFC6PtPOI->Append("M = ");
  (*entryDiffGFC6PtPOI)+=(Long_t)avMultDiffFlowGFCPOI;
  entryDiffGFC6PtPOI->Append(", N = ");
  (*entryDiffGFC6PtPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
  legendDiffFlowPtPOI->AddEntry(gfcCommonHistRes6->GetHistDiffFlowPtPOI(),entryDiffGFC6PtPOI->Data(),"p");
 } 
 if(gfcCommonHistRes8)
 {
  entryDiffGFC8PtPOI->Append("M = ");
  (*entryDiffGFC8PtPOI)+=(Long_t)avMultDiffFlowGFCPOI;
  entryDiffGFC8PtPOI->Append(", N = ");
  (*entryDiffGFC8PtPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
  legendDiffFlowPtPOI->AddEntry(gfcCommonHistRes8->GetHistDiffFlowPtPOI(),entryDiffGFC8PtPOI->Data(),"p");
 }  
 
 //QC
 if(qcCommonHistRes2)
 {
  entryDiffQC2PtPOI->Append("M = ");
  (*entryDiffQC2PtPOI)+=(Long_t)avMultDiffFlowQC2POI;
  entryDiffQC2PtPOI->Append(", N = ");
  (*entryDiffQC2PtPOI)+=(Long_t)nEvtsDiffFlowQC2POI; 
  legendDiffFlowPtPOI->AddEntry(qcCommonHistRes2->GetHistDiffFlowPtPOI(),entryDiffQC2PtPOI->Data(),"p");
 }
 if(qcCommonHistRes4)
 {
  entryDiffQC4PtPOI->Append("M = ");
  (*entryDiffQC4PtPOI)+=(Long_t)avMultDiffFlowQC4POI;
  entryDiffQC4PtPOI->Append(", N = ");
  (*entryDiffQC4PtPOI)+=(Long_t)nEvtsDiffFlowQC4POI; 
  legendDiffFlowPtPOI->AddEntry(qcCommonHistRes4->GetHistDiffFlowPtPOI(),entryDiffQC4PtPOI->Data(),"p");
 }
 
 //LYZ
 if(lyz2CommonHistRes)
 {
  entryDiffLYZ2PtPOI->Append("M = ");
  (*entryDiffLYZ2PtPOI)+=(Long_t)avMultDiffFlowLYZ2POI;
  entryDiffLYZ2PtPOI->Append(", N = ");
  (*entryDiffLYZ2PtPOI)+=(Long_t)nEvtsDiffFlowLYZ2POI; 
  legendDiffFlowPtPOI->AddEntry(lyz2CommonHistRes->GetHistDiffFlowPtPOI(),entryDiffLYZ2PtPOI->Data(),"p");
 }

 //drawing finally the legend in the 2nd pad:     
 if(legendDiffFlowPtPOI)
 {
  legendDiffFlowPtPOI->SetMargin(0.15);
  legendDiffFlowPtPOI->Draw();
 }

 //----------------------------------------------------------------------------------
 

 //----------------------------------------------------------------------------------
 //final drawing for differential flow (Eta, POI):
 TCanvas* diffFlowEtaAllCanvasPOI = new TCanvas("Differential Flow (Eta) of POI","Differential Flow (Eta) of POI ",1000,600);
 
 diffFlowEtaAllCanvasPOI->Divide(2,1);
 
 //1st pad is for plot:
 (diffFlowEtaAllCanvasPOI->cd(1))->SetPad(0.0,0.0,0.75,1.0);
 
 if(styleHistEta)
 {
  styleHistEta->Draw();
 }
 if(pMeshDiffFlowEtaPOI)
 {
  pMeshDiffFlowEtaPOI->Draw("LFSAME");
 }
 
 //MC 
 if(mcepCommonHistRes)
 { 
  (mcepCommonHistRes->GetHistDiffFlowEtaPOI())->Draw("E1PSAME");
 }
 //GFC
 if(gfcCommonHistRes2)
 { 
  (gfcCommonHistRes2->GetHistDiffFlowEtaPOI())->Draw("E1PSAME"); 
 } 
 if(gfcCommonHistRes4)
 { 
  (gfcCommonHistRes4->GetHistDiffFlowEtaPOI())->Draw("E1PSAME"); 
 } 
 if(gfcCommonHistRes6)
 { 
  (gfcCommonHistRes6->GetHistDiffFlowEtaPOI())->Draw("E1PSAME"); 
 } 
 if(gfcCommonHistRes8)
 { 
  (gfcCommonHistRes8->GetHistDiffFlowEtaPOI())->Draw("E1PSAME"); 
 }    
 //QC
 if(qcCommonHistRes2)
 { 
  (qcCommonHistRes2->GetHistDiffFlowEtaPOI())->Draw("E1PSAME");
 }
 if(qcCommonHistRes4)
 { 
  (qcCommonHistRes4->GetHistDiffFlowEtaPOI())->Draw("E1PSAME");
 }
 //LYZ2
 if(lyz2CommonHistRes)
 { 
  (lyz2CommonHistRes->GetHistDiffFlowEtaPOI())->Draw("E1PSAME");
 }
 
 //2nd pad is for legend:
 (diffFlowEtaAllCanvasPOI->cd(2))->SetPad(0.75,0.0,1.0,1.0);

 TLegend* legendDiffFlowEtaPOI = new TLegend(0.02,0.25,0.97,0.75);
 legendDiffFlowEtaPOI->SetTextFont(72);
 legendDiffFlowEtaPOI->SetTextSize(0.06);
 
 //legend's entries:
 TString *entryDiffMCEtaPOI   = new TString("MC ....... ");
 TString *entryDiffGFC2EtaPOI = new TString("GFC{2} ... ");
 TString *entryDiffGFC4EtaPOI = new TString("GFC{4} ... ");
 TString *entryDiffGFC6EtaPOI = new TString("GFC{6} ... ");
 TString *entryDiffGFC8EtaPOI = new TString("GFC{8} ... "); 
 TString *entryDiffQC2EtaPOI  = new TString("QC{2} .... ");
 TString *entryDiffQC4EtaPOI  = new TString("QC{4} .... ");
 TString *entryDiffLYZ2EtaPOI = new TString("LYZ ...... ");
 
 //MC
 if(mcepCommonHistRes)
 {
  (mcepCommonHistRes->GetHistDiffFlow())->SetFillStyle(meshStyle);
  (mcepCommonHistRes->GetHistDiffFlow())->SetFillColor(meshColor);
  entryDiffMCEtaPOI->Append("M = ");
  (*entryDiffMCEtaPOI)+=(Long_t)avMultDiffFlowMCPOI;
  entryDiffMCEtaPOI->Append(", N = ");
  (*entryDiffMCEtaPOI)+=(Long_t)nEvtsDiffFlowMCPOI; 
  legendDiffFlowEtaPOI->AddEntry(mcepCommonHistRes->GetHistDiffFlowEtaPOI(),entryDiffMCEtaPOI->Data(),"f");
 }

 //GFC
 if(gfcCommonHistRes2)
 {
  entryDiffGFC2EtaPOI->Append("M = ");
  (*entryDiffGFC2EtaPOI)+=(Long_t)avMultDiffFlowGFCPOI;
  entryDiffGFC2EtaPOI->Append(", N = ");
  (*entryDiffGFC2EtaPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
  legendDiffFlowEtaPOI->AddEntry(gfcCommonHistRes2->GetHistDiffFlowEtaPOI(),entryDiffGFC2EtaPOI->Data(),"p");
 }
 if(gfcCommonHistRes4)
 {
  entryDiffGFC4EtaPOI->Append("M = ");
  (*entryDiffGFC4EtaPOI)+=(Long_t)avMultDiffFlowGFCPOI;
  entryDiffGFC4EtaPOI->Append(", N = ");
  (*entryDiffGFC4EtaPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
  legendDiffFlowEtaPOI->AddEntry(gfcCommonHistRes4->GetHistDiffFlowEtaPOI(),entryDiffGFC4EtaPOI->Data(),"p");
 }
 if(gfcCommonHistRes6)
 {
  entryDiffGFC6EtaPOI->Append("M = ");
  (*entryDiffGFC6EtaPOI)+=(Long_t)avMultDiffFlowGFCPOI;
  entryDiffGFC6EtaPOI->Append(", N = ");
  (*entryDiffGFC6EtaPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
  legendDiffFlowEtaPOI->AddEntry(gfcCommonHistRes6->GetHistDiffFlowEtaPOI(),entryDiffGFC6EtaPOI->Data(),"p");
 } 
 if(gfcCommonHistRes8)
 {
  entryDiffGFC8EtaPOI->Append("M = ");
  (*entryDiffGFC8EtaPOI)+=(Long_t)avMultDiffFlowGFCPOI;
  entryDiffGFC8EtaPOI->Append(", N = ");
  (*entryDiffGFC8EtaPOI)+=(Long_t)nEvtsDiffFlowGFCPOI; 
  legendDiffFlowEtaPOI->AddEntry(gfcCommonHistRes8->GetHistDiffFlowEtaPOI(),entryDiffGFC8EtaPOI->Data(),"p");
 }  
 
 //QC
 if(qcCommonHistRes2)
 {
  entryDiffQC2EtaPOI->Append("M = ");
  (*entryDiffQC2EtaPOI)+=(Long_t)avMultDiffFlowQC2POI;
  entryDiffQC2EtaPOI->Append(", N = ");
  (*entryDiffQC2EtaPOI)+=(Long_t)nEvtsDiffFlowQC2POI; 
  legendDiffFlowEtaPOI->AddEntry(qcCommonHistRes2->GetHistDiffFlowEtaPOI(),entryDiffQC2EtaPOI->Data(),"p");
 }
 if(qcCommonHistRes4)
 {
  entryDiffQC4EtaPOI->Append("M = ");
  (*entryDiffQC4EtaPOI)+=(Long_t)avMultDiffFlowQC4POI;
  entryDiffQC4EtaPOI->Append(", N = ");
  (*entryDiffQC4EtaPOI)+=(Long_t)nEvtsDiffFlowQC4POI; 
  legendDiffFlowEtaPOI->AddEntry(qcCommonHistRes4->GetHistDiffFlowEtaPOI(),entryDiffQC4EtaPOI->Data(),"p");
 }
 
 //LYZ
 if(lyz2CommonHistRes)
 {
  entryDiffLYZ2EtaPOI->Append("M = ");
  (*entryDiffLYZ2EtaPOI)+=(Long_t)avMultDiffFlowLYZ2POI;
  entryDiffLYZ2EtaPOI->Append(", N = ");
  (*entryDiffLYZ2EtaPOI)+=(Long_t)nEvtsDiffFlowLYZ2POI; 
  legendDiffFlowEtaPOI->AddEntry(lyz2CommonHistRes->GetHistDiffFlowEtaPOI(),entryDiffLYZ2EtaPOI->Data(),"p");
 }

 //drawing finally the legend in the 2nd pad:     
 if(legendDiffFlowEtaPOI)
 {
  legendDiffFlowEtaPOI->SetMargin(0.15);
  legendDiffFlowEtaPOI->Draw();
 }

 //----------------------------------------------------------------------------------


 //=====================================================================================

}
