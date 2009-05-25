enum libModes {mLocal,mLocalSource};
//mLocal: Analyze data on your computer using aliroot
//mLocalSource: Analyze data on your computer using root + source files

void showSpread(const Int_t nRuns=10, TString type="",Int_t mode=mLocal)
{ 
 // load needed libraries:                       
 LoadSpreadLibraries(mode);  
 
 // standard magic:
 TString execDir(gSystem->pwd());  
 TSystemDirectory* baseDir = new TSystemDirectory(".",execDir.Data());          
 TList* dirList = baseDir->GetListOfFiles();
 Int_t nDirs = dirList->GetEntries();
 gSystem->cd(execDir);
 
 // arrays to store estimates of each method from different runs:
 // MCEP:
 Double_t mcepValueNONAME[nRuns] = {0.}; 
 Double_t mcepMaxValueNONAME = 0.;       
 Double_t mcepMinValueNONAME = 1000.;   
 
 // SP:
 Double_t spValueNONAME[nRuns] = {0.}; 
 Double_t spMaxValueNONAME = 0.;       
 Double_t spMinValueNONAME = 1000.;   
 
 // GFC{2}
 Double_t gfc2ValueNONAME[nRuns] = {0.}; 
 Double_t gfc2MaxValueNONAME = 0.;       
 Double_t gfc2MinValueNONAME = 1000.;       

 // GFC{4}
 Double_t gfc4ValueNONAME[nRuns] = {0.}; 
 Double_t gfc4MaxValueNONAME = 0.;       
 Double_t gfc4MinValueNONAME = 1000.;       

 // GFC{6}
 Double_t gfc6ValueNONAME[nRuns] = {0.}; 
 Double_t gfc6MaxValueNONAME = 0.;       
 Double_t gfc6MinValueNONAME = 1000.;       

 // GFC{8}
 Double_t gfc8ValueNONAME[nRuns] = {0.}; 
 Double_t gfc8MaxValueNONAME = 0.;       
 Double_t gfc8MinValueNONAME = 1000.;       
 
 // QC{2}
 Double_t qc2ValueNONAME[nRuns] = {0.}; 
 Double_t qc2MaxValueNONAME = 0.;       
 Double_t qc2MinValueNONAME = 1000.;       

 // QC{4}
 Double_t qc4ValueNONAME[nRuns] = {0.}; 
 Double_t qc4MaxValueNONAME = 0.;       
 Double_t qc4MinValueNONAME = 1000.;       

 // QC{6}
 Double_t qc6ValueNONAME[nRuns] = {0.}; 
 Double_t qc6MaxValueNONAME = 0.;       
 Double_t qc6MinValueNONAME = 1000.;       

 // QC{8}
 Double_t qc8ValueNONAME[nRuns] = {0.}; 
 Double_t qc8MaxValueNONAME = 0.;       
 Double_t qc8MinValueNONAME = 1000.; 
 
 // FQD:
 Double_t fqdValueNONAME[nRuns] = {0.}; 
 Double_t fqdMaxValueNONAME = 0.;       
 Double_t fqdMinValueNONAME = 1000.;   
 
 // LYZ1:
 Double_t lyz1ValueNONAME[nRuns] = {0.}; 
 Double_t lyz1MaxValueNONAME = 0.;       
 Double_t lyz1MinValueNONAME = 1000.;
 
 // LYZEP:
 Double_t lyzepValueNONAME[nRuns] = {0.}; 
 Double_t lyzepMaxValueNONAME = 0.;       
 Double_t lyzepMinValueNONAME = 1000.;            
             
 Int_t counter = 0;
  
 for(Int_t iDir=0;iDir<nDirs;++iDir)
 {
  TSystemFile* presentDir = (TSystemFile*)dirList->At(iDir);
  if(!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || 
     strcmp(presentDir->GetName(), "..") == 0) continue; 
          
  if (counter >= nRuns) break;       
                       
  TString presentDirName(gSystem->pwd()); 
  presentDirName += "/";
  presentDirName += presentDir->GetName();
  presentDirName += "/";
   
  // accessing the output .root files for each method:
  // MCEP:     
  TString fileNameMCEP = presentDirName;   
  fileNameMCEP+="outputMCEPanalysis.root";
  TFile *fileMCEP = TFile::Open(fileNameMCEP.Data(), "READ");      
  TList *listMCEP = NULL;
  AliFlowCommonHistResults *mcepCommonHistRes = NULL; 
  
  if(fileMCEP) 
  {
   fileMCEP->GetObject("cobjMCEP",listMCEP); 
   if(listMCEP) 
   {
    mcepCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listMCEP->FindObject("AliFlowCommonHistResultsMCEP")); 
   }  
  }
  
  if(mcepCommonHistRes && mcepCommonHistRes->GetHistIntFlow())
  {
   mcepValueNONAME[counter] = (mcepCommonHistRes->GetHistIntFlow())->GetBinContent(1);
   if(mcepValueNONAME[counter]>0.) // modify this condition for RPs and POIs !!!
   {
    if(mcepMaxValueNONAME < mcepValueNONAME[counter]) mcepMaxValueNONAME = mcepValueNONAME[counter]; 
    if(mcepMinValueNONAME > mcepValueNONAME[counter]) mcepMinValueNONAME = mcepValueNONAME[counter]; 
   } 
  }
  
  // SP:     
  TString fileNameSP = presentDirName;   
  fileNameSP+="outputSPanalysis.root";
  TFile *fileSP = TFile::Open(fileNameSP.Data(), "READ");      
  TList *listSP = NULL;
  AliFlowCommonHistResults *spCommonHistRes = NULL; 
  
  if(fileSP) 
  {
   fileSP->GetObject("cobjSP",listSP); 
   if(listSP) 
   {
    spCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listSP->FindObject("AliFlowCommonHistResultsSP")); 
   }  
  }
  
  if(spCommonHistRes && spCommonHistRes->GetHistIntFlow())
  {
   spValueNONAME[counter] = (spCommonHistRes->GetHistIntFlow())->GetBinContent(1);
   if(spValueNONAME[counter]>0.) // modify this condition for RPs and POIs !!!
   {
    if(spMaxValueNONAME < spValueNONAME[counter]) spMaxValueNONAME = spValueNONAME[counter]; 
    if(spMinValueNONAME > spValueNONAME[counter]) spMinValueNONAME = spValueNONAME[counter]; 
   } 
  }
  
  // GFC:     
  TString fileNameGFC = presentDirName;   
  fileNameGFC+="outputGFCanalysis.root";
  TFile *fileGFC = TFile::Open(fileNameGFC.Data(), "READ");      
  TList *listGFC = NULL;
  AliFlowCommonHistResults *gfc2CommonHistRes = NULL; 
  AliFlowCommonHistResults *gfc4CommonHistRes = NULL; 
  AliFlowCommonHistResults *gfc6CommonHistRes = NULL; 
  AliFlowCommonHistResults *gfc8CommonHistRes = NULL; 
 
  if(fileGFC) 
  {
   fileGFC->GetObject("cobjGFC",listGFC); 
   if(listGFC) 
   {
    gfc2CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults2ndOrderGFC")); 
    gfc4CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults4thOrderGFC")); 
    gfc6CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults6thOrderGFC")); 
    gfc8CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults8thOrderGFC")); 
   }  
  }
    
  if(gfc2CommonHistRes && gfc2CommonHistRes->GetHistIntFlow())
  {
   gfc2ValueNONAME[counter] = (gfc2CommonHistRes->GetHistIntFlow())->GetBinContent(1);
   if(gfc2ValueNONAME[counter]>0.) // modify this condition for RPs and POIs !!!
   {
    if(gfc2MaxValueNONAME < gfc2ValueNONAME[counter]) gfc2MaxValueNONAME = gfc2ValueNONAME[counter]; 
    if(gfc2MinValueNONAME > gfc2ValueNONAME[counter]) gfc2MinValueNONAME = gfc2ValueNONAME[counter]; 
   } 
  }
  
  if(gfc4CommonHistRes && gfc4CommonHistRes->GetHistIntFlow())
  {
   gfc4ValueNONAME[counter] = (gfc4CommonHistRes->GetHistIntFlow())->GetBinContent(1);
   if(gfc4ValueNONAME[counter]>0.) // modify this condition for RPs and POIs !!!
   {
    if(gfc4MaxValueNONAME < gfc4ValueNONAME[counter]) gfc4MaxValueNONAME = gfc4ValueNONAME[counter]; 
    if(gfc4MinValueNONAME > gfc4ValueNONAME[counter]) gfc4MinValueNONAME = gfc4ValueNONAME[counter]; 
   } 
  }
   
  if(gfc6CommonHistRes && gfc6CommonHistRes->GetHistIntFlow())
  {
   gfc6ValueNONAME[counter] = (gfc6CommonHistRes->GetHistIntFlow())->GetBinContent(1);
   if(gfc6ValueNONAME[counter]>0.) // modify this condition for RPs and POIs !!!
   {
    if(gfc6MaxValueNONAME < gfc6ValueNONAME[counter]) gfc6MaxValueNONAME = gfc6ValueNONAME[counter]; 
    if(gfc6MinValueNONAME > gfc6ValueNONAME[counter]) gfc6MinValueNONAME = gfc6ValueNONAME[counter]; 
   } 
  }
  
  if(gfc8CommonHistRes && gfc8CommonHistRes->GetHistIntFlow())
  {
   gfc8ValueNONAME[counter] = (gfc8CommonHistRes->GetHistIntFlow())->GetBinContent(1);
   if(gfc8ValueNONAME[counter]>0.) // modify this condition for RPs and POIs !!!
   {
    if(gfc8MaxValueNONAME < gfc8ValueNONAME[counter]) gfc8MaxValueNONAME = gfc8ValueNONAME[counter]; 
    if(gfc8MinValueNONAME > gfc8ValueNONAME[counter]) gfc8MinValueNONAME = gfc8ValueNONAME[counter]; 
   } 
  } 
                             
  // QC:     
  TString fileNameQC = presentDirName;   
  fileNameQC+="outputQCanalysis.root";
  TFile *fileQC = TFile::Open(fileNameQC.Data(), "READ");      
  TList *listQC = NULL;
  AliFlowCommonHistResults *qc2CommonHistRes = NULL; 
  AliFlowCommonHistResults *qc4CommonHistRes = NULL; 
  AliFlowCommonHistResults *qc6CommonHistRes = NULL; 
  AliFlowCommonHistResults *qc8CommonHistRes = NULL; 
 
  if(fileQC) 
  {
   fileQC->GetObject("cobjQC",listQC); 
   if(listQC) 
   {
    qc2CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults2ndOrderQC")); 
    qc4CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults4thOrderQC")); 
    qc6CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults6thOrderQC")); 
    qc8CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults8thOrderQC")); 
   }  
  }
    
  if(qc2CommonHistRes && qc2CommonHistRes->GetHistIntFlow())
  {
   qc2ValueNONAME[counter] = (qc2CommonHistRes->GetHistIntFlow())->GetBinContent(1);
   if(qc2ValueNONAME[counter]>0.) // modify this condition for RPs and POIs !!!
   {
    if(qc2MaxValueNONAME < qc2ValueNONAME[counter]) qc2MaxValueNONAME = qc2ValueNONAME[counter]; 
    if(qc2MinValueNONAME > qc2ValueNONAME[counter]) qc2MinValueNONAME = qc2ValueNONAME[counter]; 
   } 
  }
  
  if(qc4CommonHistRes && qc4CommonHistRes->GetHistIntFlow())
  {
   qc4ValueNONAME[counter] = (qc4CommonHistRes->GetHistIntFlow())->GetBinContent(1);
   if(qc4ValueNONAME[counter]>0.) // modify this condition for RPs and POIs !!!
   {
    if(qc4MaxValueNONAME < qc4ValueNONAME[counter]) qc4MaxValueNONAME = qc4ValueNONAME[counter]; 
    if(qc4MinValueNONAME > qc4ValueNONAME[counter]) qc4MinValueNONAME = qc4ValueNONAME[counter]; 
   } 
  }
   
  if(qc6CommonHistRes && qc6CommonHistRes->GetHistIntFlow())
  {
   qc6ValueNONAME[counter] = (qc6CommonHistRes->GetHistIntFlow())->GetBinContent(1);
   if(qc6ValueNONAME[counter]>0.) // modify this condition for RPs and POIs !!!
   {
    if(qc6MaxValueNONAME < qc6ValueNONAME[counter]) qc6MaxValueNONAME = qc6ValueNONAME[counter]; 
    if(qc6MinValueNONAME > qc6ValueNONAME[counter]) qc6MinValueNONAME = qc6ValueNONAME[counter]; 
   } 
  }
  
  if(qc8CommonHistRes && qc8CommonHistRes->GetHistIntFlow())
  {
   qc8ValueNONAME[counter] = (qc8CommonHistRes->GetHistIntFlow())->GetBinContent(1);
   if(qc8ValueNONAME[counter]>0.) // modify this condition for RPs and POIs !!!
   {
    if(qc8MaxValueNONAME < qc8ValueNONAME[counter]) qc8MaxValueNONAME = qc8ValueNONAME[counter]; 
    if(qc8MinValueNONAME > qc8ValueNONAME[counter]) qc8MinValueNONAME = qc8ValueNONAME[counter]; 
   } 
  } 
  
  // FQD:     
  TString fileNameFQD = presentDirName;   
  fileNameFQD+="outputFQDanalysis.root";
  TFile *fileFQD = TFile::Open(fileNameFQD.Data(), "READ");      
  TList *listFQD = NULL;
  AliFlowCommonHistResults *fqdCommonHistRes = NULL; 
  
  if(fileFQD) 
  {
   fileFQD->GetObject("cobjFQD",listFQD); 
   if(listFQD) 
   {
    fqdCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listFQD->FindObject("AliFlowCommonHistResultsFQD")); 
   }  
  }
  
  if(fqdCommonHistRes && fqdCommonHistRes->GetHistIntFlow())
  {
   fqdValueNONAME[counter] = (fqdCommonHistRes->GetHistIntFlow())->GetBinContent(1);
   if(fqdValueNONAME[counter]>0.) // modify this condition for RPs and POIs !!!
   {
    if(fqdMaxValueNONAME < fqdValueNONAME[counter]) fqdMaxValueNONAME = fqdValueNONAME[counter]; 
    if(fqdMinValueNONAME > fqdValueNONAME[counter]) fqdMinValueNONAME = fqdValueNONAME[counter]; 
   } 
  }
  
  // LYZ1:     
  TString fileNameLYZ1 = presentDirName;   
  fileNameLYZ1+="outputLYZ1analysis.root";
  TFile *fileLYZ1 = TFile::Open(fileNameLYZ1.Data(), "READ");      
  TList *listLYZ1 = NULL;
  AliFlowCommonHistResults *lyz1CommonHistRes = NULL; 
  
  if(fileLYZ1) 
  {
   fileLYZ1->GetObject("cobjLYZ1",listLYZ1); 
   if(listLYZ1) 
   {
    lyz1CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listLYZ1->FindObject("AliFlowCommonHistResultsLYZ1")); 
   }  
  }
  
  if(lyz1CommonHistRes && lyz1CommonHistRes->GetHistIntFlow())
  {
   lyz1ValueNONAME[counter] = (lyz1CommonHistRes->GetHistIntFlow())->GetBinContent(1);
   if(lyz1ValueNONAME[counter]>0.) // modify this condition for RPs and POIs !!!
   {
    if(lyz1MaxValueNONAME < lyz1ValueNONAME[counter]) lyz1MaxValueNONAME = lyz1ValueNONAME[counter]; 
    if(lyz1MinValueNONAME > lyz1ValueNONAME[counter]) lyz1MinValueNONAME = lyz1ValueNONAME[counter]; 
   } 
  }
  
  // LYZEP:     
  TString fileNameLYZEP = presentDirName;   
  fileNameLYZEP+="outputLYZEPanalysis.root";
  TFile *fileLYZEP = TFile::Open(fileNameLYZEP.Data(), "READ");      
  TList *listLYZEP = NULL;
  AliFlowCommonHistResults *lyzepCommonHistRes = NULL; 
  
  if(fileLYZEP) 
  {
   fileLYZEP->GetObject("cobjLYZEP",listLYZEP); 
   if(listLYZEP) 
   {
    lyzepCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listLYZEP->FindObject("AliFlowCommonHistResultsLYZEP")); 
   }  
  }
  
  if(lyzepCommonHistRes && lyzepCommonHistRes->GetHistIntFlow())
  {
   lyzepValueNONAME[counter] = (lyzepCommonHistRes->GetHistIntFlow())->GetBinContent(1);
   if(lyzepValueNONAME[counter]>0.) // modify this condition for RPs and POIs !!!
   {
    if(lyzepMaxValueNONAME < lyzepValueNONAME[counter]) lyzepMaxValueNONAME = lyzepValueNONAME[counter]; 
    if(lyzepMinValueNONAME > lyzepValueNONAME[counter]) lyzepMinValueNONAME = lyzepValueNONAME[counter]; 
   } 
  }
  
  counter++;
  
 } // end of for(Int_t iDir=0;iDir<nDirs;++iDir)
 
 
 

 
 
 // removing the title and stat. box from all histograms:
 // gStyle->SetOptTitle(0);
 gStyle->SetOptStat(0);  
 
 // box width:
 const Double_t boxWidth = 0.25;
  
 // the number of different methods:
 const Int_t nMethods = 13;
 
 // the number of runs:
 const Int_t nPoints = counter;
   
 // booking the style histogram for the integrated flow results from all methods for NONAME, RPs and POIs:
 // assigning bin numbers to methods:
 Int_t binMCEP = 1; 
 Int_t binSP = 2;
 Int_t binGFC2 = 3; 
 Int_t binGFC4 = 5; 
 Int_t binGFC6 = 7; 
 Int_t binGFC8 = 9; 
 Int_t binQC2 = 4; 
 Int_t binQC4 = 6; 
 Int_t binQC6 = 8; 
 Int_t binQC8 = 10; 
 Int_t binFQD = 11; 
 Int_t binLYZ1 = 12; 
 Int_t binLYZEP = 13; 
 
 TH1D* intFlowAll = new TH1D("intFlowAll","Integrated Flow",nMethods,0,nMethods);      
 //intFlowAll->SetLabelSize(0.036,"X");
 //intFlowAll->SetLabelSize(0.036,"Y");
 //intFlowAll->SetMarkerStyle(markerStyle);
 //intFlowAll->SetMarkerColor(markerColor);
 (intFlowAll->GetYaxis())->SetTickLength(0.01);
 (intFlowAll->GetXaxis())->SetBinLabel(binMCEP,"v_{2}{MC}");
 (intFlowAll->GetXaxis())->SetBinLabel(binSP,"v_{2}{SP}");
 (intFlowAll->GetXaxis())->SetBinLabel(binGFC2,"v_{2}{2,GFC}");
 (intFlowAll->GetXaxis())->SetBinLabel(binQC2,"v_{2}{2,QC}");
 (intFlowAll->GetXaxis())->SetBinLabel(binGFC4,"v_{2}{4,GFC}");
 (intFlowAll->GetXaxis())->SetBinLabel(binQC4,"v_{2}{4,QC}");
 (intFlowAll->GetXaxis())->SetBinLabel(binGFC6,"v_{2}{6,GFC}");
 (intFlowAll->GetXaxis())->SetBinLabel(binQC6,"v_{2}{6,QC}");
 (intFlowAll->GetXaxis())->SetBinLabel(binGFC8,"v_{2}{8,GFC}");
 (intFlowAll->GetXaxis())->SetBinLabel(binQC8,"v_{2}{8,QC}");
 (intFlowAll->GetXaxis())->SetBinLabel(binFQD,"v_{2}{FQD}");
 (intFlowAll->GetXaxis())->SetBinLabel(binLYZ1,"v_{2}{LYZ}");
 (intFlowAll->GetXaxis())->SetBinLabel(binLYZEP,"v_{2}{LYZEP}");
  
 
 
 //=============================================================================
 //                             NONAME 
 //=============================================================================
   
 TH1D *intFlowNONAME = new TH1D(*intFlowAll); 
 
 Double_t mcepNONAME[nPoints] = {0.};  
 Double_t spNONAME[nPoints] = {0.}; 
 Double_t gfc2NONAME[nPoints] = {0.};
 Double_t gfc4NONAME[nPoints] = {0.};
 Double_t gfc6NONAME[nPoints] = {0.};
 Double_t gfc8NONAME[nPoints] = {0.};
 Double_t qc2NONAME[nPoints] = {0.};
 Double_t qc4NONAME[nPoints] = {0.};
 Double_t qc6NONAME[nPoints] = {0.};
 Double_t qc8NONAME[nPoints] = {0.};
 Double_t fqdNONAME[nPoints] = {0.}; 
 Double_t lyz1NONAME[nPoints] = {0.}; 
 Double_t lyzepNONAME[nPoints] = {0.}; 
 
 for(Int_t i=0;i<nPoints;i++)
 {
  mcepNONAME[i]=binMCEP-0.5;
  spNONAME[i]=binSP-0.5;
  gfc2NONAME[i]=binGFC2-0.5;
  gfc4NONAME[i]=binGFC4-0.5;
  gfc6NONAME[i]=binGFC6-0.5;
  gfc8NONAME[i]=binGFC8-0.5;
  qc2NONAME[i]=binQC2-0.5;
  qc4NONAME[i]=binQC4-0.5;
  qc6NONAME[i]=binQC6-0.5;
  qc8NONAME[i]=binQC8-0.5;
  fqdNONAME[i]=binFQD-0.5;
  lyz1NONAME[i]=binLYZ1-0.5;
  lyzepNONAME[i]=binLYZEP-0.5;
 }
 
 // MCEP:
 TGraph* mcepTGraphNONAME = new TGraph(nPoints, mcepNONAME, mcepValueNONAME);
 mcepTGraphNONAME->SetMarkerStyle(21);
 mcepTGraphNONAME->SetMarkerColor(kBlack); 
 mcepTGraphNONAME->SetMarkerSize(0.75); 
 
 TGraph *mcepBoxNONAME = new TGraph(5);
 mcepBoxNONAME->SetPoint(1,(binMCEP-0.5)-boxWidth,mcepMinValueNONAME);
 mcepBoxNONAME->SetPoint(2,(binMCEP-0.5)+boxWidth,mcepMinValueNONAME);
 mcepBoxNONAME->SetPoint(3,(binMCEP-0.5)+boxWidth,mcepMaxValueNONAME);
 mcepBoxNONAME->SetPoint(4,(binMCEP-0.5)-boxWidth,mcepMaxValueNONAME);
 mcepBoxNONAME->SetPoint(5,(binMCEP-0.5)-boxWidth,mcepMinValueNONAME);    
 mcepBoxNONAME->SetFillStyle(1001);
 mcepBoxNONAME->SetFillColor(kGray);
  
 // SP:
 TGraph* spTGraphNONAME = new TGraph(nPoints, spNONAME, spValueNONAME);
 spTGraphNONAME->SetMarkerStyle(21);
 spTGraphNONAME->SetMarkerColor(kViolet+3); 
 spTGraphNONAME->SetMarkerSize(0.75); 
 
 TGraph *spBoxNONAME = new TGraph(5);
 spBoxNONAME->SetPoint(1,(binSP-0.5)-boxWidth,spMinValueNONAME);
 spBoxNONAME->SetPoint(2,(binSP-0.5)+boxWidth,spMinValueNONAME);
 spBoxNONAME->SetPoint(3,(binSP-0.5)+boxWidth,spMaxValueNONAME);
 spBoxNONAME->SetPoint(4,(binSP-0.5)-boxWidth,spMaxValueNONAME);
 spBoxNONAME->SetPoint(5,(binSP-0.5)-boxWidth,spMinValueNONAME);    
 spBoxNONAME->SetFillStyle(1001);
 spBoxNONAME->SetFillColor(kViolet-9);
 
 // GFC{2}:
 TGraph* gfc2TGraphNONAME = new TGraph(nPoints, gfc2NONAME, gfc2ValueNONAME);
 gfc2TGraphNONAME->SetMarkerStyle(21);
 gfc2TGraphNONAME->SetMarkerColor(kBlue);
 gfc2TGraphNONAME->SetMarkerSize(0.75); 
 
 TGraph *gfc2BoxNONAME = new TGraph(5);
 gfc2BoxNONAME->SetPoint(1,(binGFC2-0.5)-boxWidth,gfc2MinValueNONAME);
 gfc2BoxNONAME->SetPoint(2,(binGFC2-0.5)+boxWidth,gfc2MinValueNONAME);
 gfc2BoxNONAME->SetPoint(3,(binGFC2-0.5)+boxWidth,gfc2MaxValueNONAME);
 gfc2BoxNONAME->SetPoint(4,(binGFC2-0.5)-boxWidth,gfc2MaxValueNONAME);
 gfc2BoxNONAME->SetPoint(5,(binGFC2-0.5)-boxWidth,gfc2MinValueNONAME);    
 gfc2BoxNONAME->SetFillStyle(1001);
 gfc2BoxNONAME->SetFillColor(kBlue-10);
 
 // GFC{4}:
 TGraph* gfc4TGraphNONAME = new TGraph(nPoints, gfc4NONAME, gfc4ValueNONAME);
 gfc4TGraphNONAME->SetMarkerStyle(21);
 gfc4TGraphNONAME->SetMarkerColor(kBlue); 
 gfc4TGraphNONAME->SetMarkerSize(0.75); 
 
 TGraph *gfc4BoxNONAME = new TGraph(5);
 gfc4BoxNONAME->SetPoint(1,(binGFC4-0.5)-boxWidth,gfc4MinValueNONAME);
 gfc4BoxNONAME->SetPoint(2,(binGFC4-0.5)+boxWidth,gfc4MinValueNONAME);
 gfc4BoxNONAME->SetPoint(3,(binGFC4-0.5)+boxWidth,gfc4MaxValueNONAME);
 gfc4BoxNONAME->SetPoint(4,(binGFC4-0.5)-boxWidth,gfc4MaxValueNONAME);
 gfc4BoxNONAME->SetPoint(5,(binGFC4-0.5)-boxWidth,gfc4MinValueNONAME);    
 gfc4BoxNONAME->SetFillStyle(1001);
 gfc4BoxNONAME->SetFillColor(kBlue-10);

 // GFC{6}:
 TGraph* gfc6TGraphNONAME = new TGraph(nPoints, gfc6NONAME, gfc6ValueNONAME);
 gfc6TGraphNONAME->SetMarkerStyle(21);
 gfc6TGraphNONAME->SetMarkerColor(kBlue); 
 gfc6TGraphNONAME->SetMarkerSize(0.75); 
 
 TGraph *gfc6BoxNONAME = new TGraph(5);
 gfc6BoxNONAME->SetPoint(1,(binGFC6-0.5)-boxWidth,gfc6MinValueNONAME);
 gfc6BoxNONAME->SetPoint(2,(binGFC6-0.5)+boxWidth,gfc6MinValueNONAME);
 gfc6BoxNONAME->SetPoint(3,(binGFC6-0.5)+boxWidth,gfc6MaxValueNONAME);
 gfc6BoxNONAME->SetPoint(4,(binGFC6-0.5)-boxWidth,gfc6MaxValueNONAME);
 gfc6BoxNONAME->SetPoint(5,(binGFC6-0.5)-boxWidth,gfc6MinValueNONAME);    
 gfc6BoxNONAME->SetFillStyle(1001);
 gfc6BoxNONAME->SetFillColor(kBlue-10);

 // GFC{8}:
 TGraph* gfc8TGraphNONAME = new TGraph(nPoints, gfc8NONAME, gfc8ValueNONAME);
 gfc8TGraphNONAME->SetMarkerStyle(21);
 gfc8TGraphNONAME->SetMarkerColor(kBlue); 
 gfc8TGraphNONAME->SetMarkerSize(0.75); 
 
 TGraph *gfc8BoxNONAME = new TGraph(5);
 gfc8BoxNONAME->SetPoint(1,(binGFC8-0.5)-boxWidth,gfc8MinValueNONAME);
 gfc8BoxNONAME->SetPoint(2,(binGFC8-0.5)+boxWidth,gfc8MinValueNONAME);
 gfc8BoxNONAME->SetPoint(3,(binGFC8-0.5)+boxWidth,gfc8MaxValueNONAME);
 gfc8BoxNONAME->SetPoint(4,(binGFC8-0.5)-boxWidth,gfc8MaxValueNONAME);
 gfc8BoxNONAME->SetPoint(5,(binGFC8-0.5)-boxWidth,gfc8MinValueNONAME);    
 gfc8BoxNONAME->SetFillStyle(1001);
 gfc8BoxNONAME->SetFillColor(kBlue-10);
 
 // QC{2}:
 TGraph* qc2TGraphNONAME = new TGraph(nPoints, qc2NONAME, qc2ValueNONAME);
 qc2TGraphNONAME->SetMarkerStyle(21);
 qc2TGraphNONAME->SetMarkerColor(kRed); 
 qc2TGraphNONAME->SetMarkerSize(0.75); 
 
 TGraph *qc2BoxNONAME = new TGraph(5);
 qc2BoxNONAME->SetPoint(1,(binQC2-0.5)-boxWidth,qc2MinValueNONAME);
 qc2BoxNONAME->SetPoint(2,(binQC2-0.5)+boxWidth,qc2MinValueNONAME);
 qc2BoxNONAME->SetPoint(3,(binQC2-0.5)+boxWidth,qc2MaxValueNONAME);
 qc2BoxNONAME->SetPoint(4,(binQC2-0.5)-boxWidth,qc2MaxValueNONAME);
 qc2BoxNONAME->SetPoint(5,(binQC2-0.5)-boxWidth,qc2MinValueNONAME);    
 qc2BoxNONAME->SetFillStyle(1001);
 qc2BoxNONAME->SetFillColor(kRed-10);
 
 // QC{4}:
 TGraph* qc4TGraphNONAME = new TGraph(nPoints, qc4NONAME, qc4ValueNONAME);
 qc4TGraphNONAME->SetMarkerStyle(21);
 qc4TGraphNONAME->SetMarkerColor(kRed); 
 qc4TGraphNONAME->SetMarkerSize(0.75); 
 
 TGraph *qc4BoxNONAME = new TGraph(5);
 qc4BoxNONAME->SetPoint(1,(binQC4-0.5)-boxWidth,qc4MinValueNONAME);
 qc4BoxNONAME->SetPoint(2,(binQC4-0.5)+boxWidth,qc4MinValueNONAME);
 qc4BoxNONAME->SetPoint(3,(binQC4-0.5)+boxWidth,qc4MaxValueNONAME);
 qc4BoxNONAME->SetPoint(4,(binQC4-0.5)-boxWidth,qc4MaxValueNONAME);
 qc4BoxNONAME->SetPoint(5,(binQC4-0.5)-boxWidth,qc4MinValueNONAME);    
 qc4BoxNONAME->SetFillStyle(1001);
 qc4BoxNONAME->SetFillColor(kRed-10);

 // QC{6}:
 TGraph* qc6TGraphNONAME = new TGraph(nPoints, qc6NONAME, qc6ValueNONAME);
 qc6TGraphNONAME->SetMarkerStyle(21);
 qc6TGraphNONAME->SetMarkerColor(kRed); 
 qc6TGraphNONAME->SetMarkerSize(0.75); 
 
 TGraph *qc6BoxNONAME = new TGraph(5);
 qc6BoxNONAME->SetPoint(1,(binQC6-0.5)-boxWidth,qc6MinValueNONAME);
 qc6BoxNONAME->SetPoint(2,(binQC6-0.5)+boxWidth,qc6MinValueNONAME);
 qc6BoxNONAME->SetPoint(3,(binQC6-0.5)+boxWidth,qc6MaxValueNONAME);
 qc6BoxNONAME->SetPoint(4,(binQC6-0.5)-boxWidth,qc6MaxValueNONAME);
 qc6BoxNONAME->SetPoint(5,(binQC6-0.5)-boxWidth,qc6MinValueNONAME);    
 qc6BoxNONAME->SetFillStyle(1001);
 qc6BoxNONAME->SetFillColor(kRed-10);

 // QC{8}:
 TGraph* qc8TGraphNONAME = new TGraph(nPoints, qc8NONAME, qc8ValueNONAME);
 qc8TGraphNONAME->SetMarkerStyle(21);
 qc8TGraphNONAME->SetMarkerColor(kRed); 
 qc8TGraphNONAME->SetMarkerSize(0.75); 
 
 TGraph *qc8BoxNONAME = new TGraph(5);
 qc8BoxNONAME->SetPoint(1,(binQC8-0.5)-boxWidth,qc8MinValueNONAME);
 qc8BoxNONAME->SetPoint(2,(binQC8-0.5)+boxWidth,qc8MinValueNONAME);
 qc8BoxNONAME->SetPoint(3,(binQC8-0.5)+boxWidth,qc8MaxValueNONAME);
 qc8BoxNONAME->SetPoint(4,(binQC8-0.5)-boxWidth,qc8MaxValueNONAME);
 qc8BoxNONAME->SetPoint(5,(binQC8-0.5)-boxWidth,qc8MinValueNONAME);    
 qc8BoxNONAME->SetFillStyle(1001);
 qc8BoxNONAME->SetFillColor(kRed-10);
 
 // FQD:
 TGraph* fqdTGraphNONAME = new TGraph(nPoints, fqdNONAME, fqdValueNONAME);
 fqdTGraphNONAME->SetMarkerStyle(21);
 fqdTGraphNONAME->SetMarkerColor(kOrange+7); 
 fqdTGraphNONAME->SetMarkerSize(0.75); 
 
 TGraph *fqdBoxNONAME = new TGraph(5);
 fqdBoxNONAME->SetPoint(1,(binFQD-0.5)-boxWidth,fqdMinValueNONAME);
 fqdBoxNONAME->SetPoint(2,(binFQD-0.5)+boxWidth,fqdMinValueNONAME);
 fqdBoxNONAME->SetPoint(3,(binFQD-0.5)+boxWidth,fqdMaxValueNONAME);
 fqdBoxNONAME->SetPoint(4,(binFQD-0.5)-boxWidth,fqdMaxValueNONAME);
 fqdBoxNONAME->SetPoint(5,(binFQD-0.5)-boxWidth,fqdMinValueNONAME);    
 fqdBoxNONAME->SetFillStyle(1001);
 fqdBoxNONAME->SetFillColor(kOrange-9);
 
 // LYZ1:
 TGraph* lyz1TGraphNONAME = new TGraph(nPoints, lyz1NONAME, lyz1ValueNONAME);
 lyz1TGraphNONAME->SetMarkerStyle(21);
 lyz1TGraphNONAME->SetMarkerColor(kYellow+3); 
 lyz1TGraphNONAME->SetMarkerSize(0.75); 
  
 TGraph *lyz1BoxNONAME = new TGraph(5);
 lyz1BoxNONAME->SetPoint(1,(binLYZ1-0.5)-boxWidth,lyz1MinValueNONAME);
 lyz1BoxNONAME->SetPoint(2,(binLYZ1-0.5)+boxWidth,lyz1MinValueNONAME);
 lyz1BoxNONAME->SetPoint(3,(binLYZ1-0.5)+boxWidth,lyz1MaxValueNONAME);
 lyz1BoxNONAME->SetPoint(4,(binLYZ1-0.5)-boxWidth,lyz1MaxValueNONAME);
 lyz1BoxNONAME->SetPoint(5,(binLYZ1-0.5)-boxWidth,lyz1MinValueNONAME);    
 lyz1BoxNONAME->SetFillStyle(1001);
 lyz1BoxNONAME->SetFillColor(kYellow-8);
   
 // LYZEP:
 TGraph* lyzepTGraphNONAME = new TGraph(nPoints, lyzepNONAME, lyzepValueNONAME);
 lyzepTGraphNONAME->SetMarkerStyle(21);
 lyzepTGraphNONAME->SetMarkerColor(kGreen+3); 
 lyzepTGraphNONAME->SetMarkerSize(0.75); 
 
 TGraph *lyzepBoxNONAME = new TGraph(5);
 lyzepBoxNONAME->SetPoint(1,(binLYZEP-0.5)-boxWidth,lyzepMinValueNONAME);
 lyzepBoxNONAME->SetPoint(2,(binLYZEP-0.5)+boxWidth,lyzepMinValueNONAME);
 lyzepBoxNONAME->SetPoint(3,(binLYZEP-0.5)+boxWidth,lyzepMaxValueNONAME);
 lyzepBoxNONAME->SetPoint(4,(binLYZEP-0.5)-boxWidth,lyzepMaxValueNONAME);
 lyzepBoxNONAME->SetPoint(5,(binLYZEP-0.5)-boxWidth,lyzepMinValueNONAME);    
 lyzepBoxNONAME->SetFillStyle(1001);
 lyzepBoxNONAME->SetFillColor(kGreen-9);
   
 TCanvas* intFlowCanvasNONAME = new TCanvas("Integrated Flow NONAME","Integrated Flow NONAME",1000,600);
 
 intFlowCanvasNONAME->Divide(2,1);
 
 // 1st pad is for plot:
 (intFlowCanvasNONAME->cd(1))->SetPad(0.0,0.0,0.75,1.0);
 
 if(intFlowNONAME) 
 {
  TString intFlowNameNONAME("Superimposing ");
  intFlowNameNONAME+=nRuns;
  intFlowNameNONAME+=" independent runs";
  intFlowNONAME->SetTitle(intFlowNameNONAME.Data());
  (intFlowNONAME->GetYaxis())->SetRangeUser(0,0.20);
  intFlowNONAME->Draw();
 }
 
 // MCEP 
 if(mcepBoxNONAME && mcepMinValueNONAME < 1000.) mcepBoxNONAME->Draw("LFSAME");
 if(mcepTGraphNONAME) mcepTGraphNONAME->Draw("PSAME");
 
 // SP 
 if(spBoxNONAME && spMinValueNONAME < 1000.) spBoxNONAME->Draw("LFSAME");
 if(spTGraphNONAME) spTGraphNONAME->Draw("PSAME");
 
 // GFC{2} 
 if(gfc2BoxNONAME && gfc2MinValueNONAME < 1000.) gfc2BoxNONAME->Draw("LFSAME");
 if(gfc2TGraphNONAME) gfc2TGraphNONAME->Draw("PSAME");
  
 // GFC{4} 
 if(gfc4BoxNONAME && gfc4MinValueNONAME < 1000.) gfc4BoxNONAME->Draw("LFSAME");
 if(gfc4TGraphNONAME) gfc4TGraphNONAME->Draw("PSAME");

 // GFC{6} 
 if(gfc6BoxNONAME && gfc6MinValueNONAME < 1000.) gfc6BoxNONAME->Draw("LFSAME");
 if(gfc6TGraphNONAME) gfc6TGraphNONAME->Draw("PSAME");

 // GFC{8} 
 if(gfc8BoxNONAME && gfc8MinValueNONAME < 1000.) gfc8BoxNONAME->Draw("LFSAME");
 if(gfc8TGraphNONAME) gfc8TGraphNONAME->Draw("PSAME");
 
 // QC{2} 
 if(qc2BoxNONAME && qc2MinValueNONAME < 1000.) qc2BoxNONAME->Draw("LFSAME");
 if(qc2TGraphNONAME) qc2TGraphNONAME->Draw("PSAME");
  
 // QC{4} 
 if(qc4BoxNONAME && qc4MinValueNONAME < 1000.) qc4BoxNONAME->Draw("LFSAME");
 if(qc4TGraphNONAME) qc4TGraphNONAME->Draw("PSAME");

 // QC{6} 
 if(qc6BoxNONAME && qc6MinValueNONAME < 1000.) qc6BoxNONAME->Draw("LFSAME");
 if(qc6TGraphNONAME) qc6TGraphNONAME->Draw("PSAME");

 // QC{8} 
 if(qc8BoxNONAME && qc8MinValueNONAME < 1000.) qc8BoxNONAME->Draw("LFSAME");
 if(qc8TGraphNONAME) qc8TGraphNONAME->Draw("PSAME");
 
 // FQD 
 if(fqdBoxNONAME && fqdMinValueNONAME < 1000.) fqdBoxNONAME->Draw("LFSAME");
 if(fqdTGraphNONAME) fqdTGraphNONAME->Draw("PSAME");
 
 // LYZ1 
 if(lyz1BoxNONAME && lyz1MinValueNONAME < 1000.) lyz1BoxNONAME->Draw("LFSAME");
 if(lyz1TGraphNONAME) lyz1TGraphNONAME->Draw("PSAME");
 
 // LYZEP 
 if(lyzepBoxNONAME && lyzepMinValueNONAME < 1000.) lyzepBoxNONAME->Draw("LFSAME");
 if(lyzepTGraphNONAME) lyzepTGraphNONAME->Draw("PSAME");
 
 // 2nd pad is for legend:   
 (intFlowCanvasNONAME->cd(2))->SetPad(0.75,0.0,1.0,1.0);
 
 // count real estimates:
 Int_t mcepCountRealNONAME = 0;
 Int_t spCountRealNONAME = 0;
 Int_t gfc2CountRealNONAME = 0;
 Int_t gfc4CountRealNONAME = 0;
 Int_t gfc6CountRealNONAME = 0;
 Int_t gfc8CountRealNONAME = 0;
 Int_t qc2CountRealNONAME = 0;
 Int_t qc4CountRealNONAME = 0;
 Int_t qc6CountRealNONAME = 0;
 Int_t qc8CountRealNONAME = 0;
 Int_t fqdCountRealNONAME = 0;
 Int_t lyz1CountRealNONAME = 0;
 Int_t lyzepCountRealNONAME = 0;
 for(Int_t i=0;i<nRuns;i++)
 {
  if(mcepValueNONAME[i]>0.) mcepCountRealNONAME++; 
  if(spValueNONAME[i]>0.) spCountRealNONAME++; 
  if(gfc2ValueNONAME[i]>0.) gfc2CountRealNONAME++; 
  if(gfc4ValueNONAME[i]>0.) gfc4CountRealNONAME++; 
  if(gfc6ValueNONAME[i]>0.) gfc6CountRealNONAME++; 
  if(gfc8ValueNONAME[i]>0.) gfc8CountRealNONAME++; 
  if(qc2ValueNONAME[i]>0.) qc2CountRealNONAME++; 
  if(qc4ValueNONAME[i]>0.) qc4CountRealNONAME++; 
  if(qc6ValueNONAME[i]>0.) qc6CountRealNONAME++; 
  if(qc8ValueNONAME[i]>0.) qc8CountRealNONAME++;
  if(fqdValueNONAME[i]>0.) fqdCountRealNONAME++; 
  if(lyz1ValueNONAME[i]>0.) lyz1CountRealNONAME++; 
  if(lyzepValueNONAME[i]>0.) lyzepCountRealNONAME++; 
 }

 TPaveText *textDefaultNONAME = new TPaveText(0.05,0.67,0.95,0.90,"NDC");
 textDefaultNONAME->SetTextFont(72);
 textDefaultNONAME->SetTextSize(0.08);
  
 TString *entryDefaultRealNONAME = new TString("Real estimates");
 //TString *entryDefaultOutOfNONAME = new TString("out of"); 
 TString *entryDefaultTotalNumberNONAME = new TString(" out of total number of ");
 TString *entryDefaultTotalIndNONAME = new TString("independent");
 TString *entryDefaultTotalSimNONAME = new TString("simulations:");
  
 textDefaultNONAME->AddText(entryDefaultRealNONAME->Data());
 //textDefaultNONAME->AddText(entryDefaultOutOfNONAME->Data());
 textDefaultNONAME->AddText(entryDefaultTotalNumberNONAME->Data());
 textDefaultNONAME->AddText(entryDefaultTotalIndNONAME->Data());
 textDefaultNONAME->AddText(entryDefaultTotalSimNONAME->Data());
 
 // results:
 TPaveText *textResultsNONAME = new TPaveText(0.05,0.12,0.95,0.60,"NDC");
 textResultsNONAME->SetTextFont(72);
 textResultsNONAME->SetTextSize(0.06);
 
 // entries:
 TString *entryIntFlowMCNONAME    = new TString("MC ................ ");
 TString *entryIntFlowSPNONAME    = new TString("SP ................ ");  
 TString *entryIntFlowGFC2NONAME  = new TString("GFC{2} ........ ");
 TString *entryIntFlowGFC4NONAME  = new TString("GFC{4} ........ ");
 TString *entryIntFlowGFC6NONAME  = new TString("GFC{6} ........ ");
 TString *entryIntFlowGFC8NONAME  = new TString("GFC{8} ........ "); 
 TString *entryIntFlowQC2NONAME   = new TString("QC{2} .......... ");
 TString *entryIntFlowQC4NONAME   = new TString("QC{4} .......... ");
 TString *entryIntFlowQC6NONAME   = new TString("QC{6} .......... ");
 TString *entryIntFlowQC8NONAME   = new TString("QC{8} .......... ");
 TString *entryIntFlowFQDNONAME   = new TString("FQD ............. ");
 TString *entryIntFlowLYZ1NONAME  = new TString("LYZ ............. ");
 TString *entryIntFlowLYZEPNONAME = new TString("LYZEP ........ ");

 if(entryIntFlowMCNONAME)
 { 
  (*entryIntFlowMCNONAME)+=(Long_t)mcepCountRealNONAME;
  entryIntFlowMCNONAME->Append(" out of ");
  (*entryIntFlowMCNONAME)+=(Long_t)counter;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowMCNONAME->Data());
 }
 
 if(entryIntFlowSPNONAME)
 { 
  (*entryIntFlowSPNONAME)+=(Long_t)spCountRealNONAME;
  entryIntFlowSPNONAME->Append(" out of ");
  (*entryIntFlowSPNONAME)+=(Long_t)counter;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowSPNONAME->Data());
 }
 
 if(entryIntFlowGFC2NONAME)
 { 
  (*entryIntFlowGFC2NONAME)+=(Long_t)gfc2CountRealNONAME;
  entryIntFlowGFC2NONAME->Append(" out of ");
  (*entryIntFlowGFC2NONAME)+=(Long_t)counter;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowGFC2NONAME->Data());
 }
 
 if(entryIntFlowGFC4NONAME)
 { 
  (*entryIntFlowGFC4NONAME)+=(Long_t)gfc4CountRealNONAME;
  entryIntFlowGFC4NONAME->Append(" out of ");
  (*entryIntFlowGFC4NONAME)+=(Long_t)counter;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowGFC4NONAME->Data());
 }

 if(entryIntFlowGFC6NONAME)
 { 
  (*entryIntFlowGFC6NONAME)+=(Long_t)gfc6CountRealNONAME;
  entryIntFlowGFC6NONAME->Append(" out of ");
  (*entryIntFlowGFC6NONAME)+=(Long_t)counter;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowGFC6NONAME->Data());
 }

 if(entryIntFlowGFC8NONAME)
 { 
  (*entryIntFlowGFC8NONAME)+=(Long_t)gfc8CountRealNONAME;
  entryIntFlowGFC8NONAME->Append(" out of ");
  (*entryIntFlowGFC8NONAME)+=(Long_t)counter;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowGFC8NONAME->Data());
 }
 
 if(entryIntFlowQC2NONAME)
 { 
  (*entryIntFlowQC2NONAME)+=(Long_t)qc2CountRealNONAME;
  entryIntFlowQC2NONAME->Append(" out of ");
  (*entryIntFlowQC2NONAME)+=(Long_t)counter;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowQC2NONAME->Data());
 }
 
 if(entryIntFlowQC4NONAME)
 { 
  (*entryIntFlowQC4NONAME)+=(Long_t)qc4CountRealNONAME;
  entryIntFlowQC4NONAME->Append(" out of ");
  (*entryIntFlowQC4NONAME)+=(Long_t)counter;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowQC4NONAME->Data());
 }

 if(entryIntFlowQC2NONAME)
 { 
  (*entryIntFlowQC6NONAME)+=(Long_t)qc6CountRealNONAME;
  entryIntFlowQC6NONAME->Append(" out of ");
  (*entryIntFlowQC6NONAME)+=(Long_t)counter;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowQC6NONAME->Data());
 }

 if(entryIntFlowQC8NONAME)
 { 
  (*entryIntFlowQC8NONAME)+=(Long_t)qc8CountRealNONAME;
  entryIntFlowQC8NONAME->Append(" out of ");
  (*entryIntFlowQC8NONAME)+=(Long_t)counter;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowQC8NONAME->Data());
 }
 
 if(entryIntFlowFQDNONAME)
 { 
  (*entryIntFlowFQDNONAME)+=(Long_t)fqdCountRealNONAME;
  entryIntFlowFQDNONAME->Append(" out of ");
  (*entryIntFlowFQDNONAME)+=(Long_t)counter;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowFQDNONAME->Data());
 }
 
 if(entryIntFlowLYZ1NONAME)
 { 
  (*entryIntFlowLYZ1NONAME)+=(Long_t)lyz1CountRealNONAME;
  entryIntFlowLYZ1NONAME->Append(" out of ");
  (*entryIntFlowLYZ1NONAME)+=(Long_t)counter;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowLYZ1NONAME->Data());
 }
 
 if(entryIntFlowLYZEPNONAME)
 { 
  (*entryIntFlowLYZEPNONAME)+=(Long_t)lyzepCountRealNONAME;
  entryIntFlowLYZEPNONAME->Append(" out of ");
  (*entryIntFlowLYZEPNONAME)+=(Long_t)counter;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowLYZEPNONAME->Data());
 }
 
 if(textDefaultNONAME) textDefaultNONAME->Draw();
 if(textResultsNONAME) textResultsNONAME->Draw();
 








// ===============================================================================================
//                  calculate results from the merged ouput for each method
// ===============================================================================================

 TString pwd;
 
 (intFlowCanvasNONAME->cd(1))->SetPad(0.0,0.0,0.75,1.0); // to be improved

 // MCEP:
 TString mergedFileNameMCEP("outputMCEPanalysis.root");
 TFile* fileMCEP = NULL;
 TList *listHistosMCEP = NULL;
 AliFlowAnalysisWithMCEventPlane* mcep = new AliFlowAnalysisWithMCEventPlane();
 AliFlowCommonHistResults *mcepCommonHistRes = NULL; 
 TGraph *meanValueMCEP = new TGraph(1);
 pwd = gSystem->pwd();
 
 if(gSystem->AccessPathName(((pwd+="/")+=mergedFileNameMCEP).Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file for MCEP !!!!"<<endl;
 } else
   { 
    fileMCEP = TFile::Open(mergedFileNameMCEP.Data(),"UPDATE");
    if(fileMCEP) fileMCEP->GetObject("cobjMCEP",listHistosMCEP);
    if(listHistosMCEP) mcep->GetOutputHistograms(listHistosMCEP);
    mcep->Finish(); 
    mcepCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listHistosMCEP->FindObject("AliFlowCommonHistResultsMCEP")); 
    if(mcepCommonHistRes) 
    {
     meanValueMCEP->SetPoint(1,binMCEP-0.5,(mcepCommonHistRes->GetHistIntFlow())->GetBinContent(1));
     meanValueMCEP->SetMarkerStyle(25);
     meanValueMCEP->SetMarkerColor(kBlack); 
     meanValueMCEP->SetMarkerSize(1.25); 
     meanValueMCEP->Draw("PSAME");
    } 
    fileMCEP->Close();
    fileMCEP = TFile::Open(mergedFileNameMCEP.Data(),"RECREATE");
    listHistosMCEP->SetName("cobjMCEP");
    listHistosMCEP->Write(listHistosMCEP->GetName(), TObject::kSingleKey);
   }
   
 // SP:
 TString mergedFileNameSP("outputSPanalysis.root");
 TFile* fileSP = NULL;
 TList *listHistosSP = NULL;
 AliFlowAnalysisWithScalarProduct* sp = new AliFlowAnalysisWithScalarProduct();
 AliFlowCommonHistResults *spCommonHistRes = NULL; 
 TGraph *meanValueSP = new TGraph(1);
 pwd = gSystem->pwd();
 
 if(gSystem->AccessPathName(((pwd+="/")+=mergedFileNameSP).Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file for SP !!!!"<<endl;
 } else
   { 
    fileSP = TFile::Open(mergedFileNameSP.Data(),"UPDATE");
    if(fileSP) fileSP->GetObject("cobjSP",listHistosSP);
    if(listHistosSP) sp->GetOutputHistograms(listHistosSP);
    sp->Finish(); 
    spCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listHistosSP->FindObject("AliFlowCommonHistResultsSP")); 
    if(spCommonHistRes) 
    {
     meanValueSP->SetPoint(1,binSP-0.5,(spCommonHistRes->GetHistIntFlow())->GetBinContent(1));
     meanValueSP->SetMarkerStyle(25);
     meanValueSP->SetMarkerColor(kViolet+3); 
     meanValueSP->SetMarkerSize(1.25); 
     meanValueSP->Draw("PSAME");
    } 
    fileSP->Close();
    fileSP = TFile::Open(mergedFileNameSP.Data(),"RECREATE");
    listHistosSP->SetName("cobjSP");
    listHistosSP->Write(listHistosSP->GetName(), TObject::kSingleKey);
   }
 
 // GFC:
 TString mergedFileNameGFC("outputGFCanalysis.root");
 TFile* fileGFC = NULL;
 TList *listHistosGFC = NULL;
 AliFlowAnalysisWithCumulants* gfc = new AliFlowAnalysisWithCumulants(); 
 AliFlowCommonHistResults *gfc2CommonHistRes = NULL; 
 AliFlowCommonHistResults *gfc4CommonHistRes = NULL; 
 AliFlowCommonHistResults *gfc6CommonHistRes = NULL; 
 AliFlowCommonHistResults *gfc8CommonHistRes = NULL; 
 TGraph *meanValueGFC2 = new TGraph(1);
 TGraph *meanValueGFC4 = new TGraph(1);
 TGraph *meanValueGFC6 = new TGraph(1);
 TGraph *meanValueGFC8 = new TGraph(1);
 pwd = gSystem->pwd();
 
 if(gSystem->AccessPathName(((pwd+="/")+=mergedFileNameGFC).Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file for GFC !!!!"<<endl;
 } else
   { 
    fileGFC = TFile::Open(mergedFileNameGFC.Data(),"UPDATE");
    if(fileGFC) fileGFC->GetObject("cobjGFC",listHistosGFC);
    if(listHistosGFC) gfc->GetOutputHistograms(listHistosGFC);
    gfc->Finish(); 
    gfc2CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults2ndOrderGFC")); 
    gfc4CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults4thOrderGFC")); 
    gfc6CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults6thOrderGFC")); 
    gfc8CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults8thOrderGFC")); 
    if(gfc2CommonHistRes) 
    {
     meanValueGFC2->SetPoint(1,binGFC2-0.5,(gfc2CommonHistRes->GetHistIntFlow())->GetBinContent(1));
     meanValueGFC2->SetMarkerStyle(25);
     meanValueGFC2->SetMarkerColor(kBlue); 
     meanValueGFC2->SetMarkerSize(1.25); 
     meanValueGFC2->Draw("PSAME");
    } 
    if(gfc4CommonHistRes) 
    {
     meanValueGFC4->SetPoint(1,binGFC4-0.5,(gfc4CommonHistRes->GetHistIntFlow())->GetBinContent(1));
     meanValueGFC4->SetMarkerStyle(25);
     meanValueGFC4->SetMarkerColor(kBlue); 
     meanValueGFC4->SetMarkerSize(1.25); 
     meanValueGFC4->Draw("PSAME");
    } 
    if(gfc6CommonHistRes) 
    {
     meanValueGFC6->SetPoint(1,binGFC6-0.5,(gfc6CommonHistRes->GetHistIntFlow())->GetBinContent(1));
     meanValueGFC6->SetMarkerStyle(25);
     meanValueGFC6->SetMarkerColor(kBlue); 
     meanValueGFC6->SetMarkerSize(1.25); 
     meanValueGFC6->Draw("PSAME");
    } 
    if(gfc8CommonHistRes) 
    {
     meanValueGFC8->SetPoint(1,binGFC8-0.5,(gfc8CommonHistRes->GetHistIntFlow())->GetBinContent(1));
     meanValueGFC8->SetMarkerStyle(25);
     meanValueGFC8->SetMarkerColor(kBlue); 
     meanValueGFC8->SetMarkerSize(1.25); 
     meanValueGFC8->Draw("PSAME");
    } 
    fileGFC->Close();
    fileGFC = TFile::Open(mergedFileNameGFC.Data(),"RECREATE");
    listHistosGFC->SetName("cobjGFC");
    listHistosGFC->Write(listHistosGFC->GetName(), TObject::kSingleKey);
   }
  
 // QC:
 TString mergedFileNameQC("outputQCanalysis.root");
 TFile* fileQC = NULL;
 TList *listHistosQC = NULL;
 AliFlowAnalysisWithQCumulants* qc = new AliFlowAnalysisWithQCumulants(); 
 AliFlowCommonHistResults *qc2CommonHistRes = NULL; 
 AliFlowCommonHistResults *qc4CommonHistRes = NULL; 
 AliFlowCommonHistResults *qc6CommonHistRes = NULL; 
 AliFlowCommonHistResults *qc8CommonHistRes = NULL; 
 TGraph *meanValueQC2 = new TGraph(1);
 TGraph *meanValueQC4 = new TGraph(1);
 TGraph *meanValueQC6 = new TGraph(1);
 TGraph *meanValueQC8 = new TGraph(1);
 pwd = gSystem->pwd();
 
 if(gSystem->AccessPathName(((pwd+="/")+=mergedFileNameQC).Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file for QC !!!!"<<endl;
 } else
   { 
    fileQC = TFile::Open(mergedFileNameQC.Data(),"UPDATE");
    if(fileQC) fileQC->GetObject("cobjQC",listHistosQC);
    if(listHistosQC) qc->GetOutputHistograms(listHistosQC);
    qc->Finish(); 
    qc2CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults2ndOrderQC")); 
    qc4CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults4thOrderQC")); 
    qc6CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults6thOrderQC")); 
    qc8CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults8thOrderQC")); 
    if(qc2CommonHistRes) 
    {
     meanValueQC2->SetPoint(1,binQC2-0.5,(qc2CommonHistRes->GetHistIntFlow())->GetBinContent(1));
     meanValueQC2->SetMarkerStyle(25);
     meanValueQC2->SetMarkerColor(kRed); 
     meanValueQC2->SetMarkerSize(1.25); 
     meanValueQC2->Draw("PSAME");
    } 
    if(qc4CommonHistRes) 
    {
     meanValueQC4->SetPoint(1,binQC4-0.5,(qc4CommonHistRes->GetHistIntFlow())->GetBinContent(1));
     meanValueQC4->SetMarkerStyle(25);
     meanValueQC4->SetMarkerColor(kRed); 
     meanValueQC4->SetMarkerSize(1.25); 
     meanValueQC4->Draw("PSAME");
    } 
    if(qc6CommonHistRes) 
    {
     meanValueQC6->SetPoint(1,binQC6-0.5,(qc6CommonHistRes->GetHistIntFlow())->GetBinContent(1));
     meanValueQC6->SetMarkerStyle(25);
     meanValueQC6->SetMarkerColor(kRed); 
     meanValueQC6->SetMarkerSize(1.25); 
     meanValueQC6->Draw("PSAME");
    } 
    if(qc8CommonHistRes) 
    {
     meanValueQC8->SetPoint(1,binQC8-0.5,(qc8CommonHistRes->GetHistIntFlow())->GetBinContent(1));
     meanValueQC8->SetMarkerStyle(25);
     meanValueQC8->SetMarkerColor(kRed); 
     meanValueQC8->SetMarkerSize(1.25); 
     meanValueQC8->Draw("PSAME");
    } 
    fileQC->Close();
    fileQC = TFile::Open(mergedFileNameQC.Data(),"RECREATE");
    listHistosQC->SetName("cobjQC");
    listHistosQC->Write(listHistosGFC->GetName(), TObject::kSingleKey);
   }
  
 // FQD:
 TString mergedFileNameFQD("outputFQDanalysis.root");
 TFile* fileFQD = NULL;
 TList *listHistosFQD = NULL;
 AliFittingQDistribution* fqd = new AliFittingQDistribution();
 AliFlowCommonHistResults *fqdCommonHistRes = NULL; 
 TGraph *meanValueFQD = new TGraph(1);
 pwd = gSystem->pwd();
 
 if(gSystem->AccessPathName(((pwd+="/")+=mergedFileNameFQD).Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file for FQD !!!!"<<endl;
 } else
   { 
    fileFQD = TFile::Open(mergedFileNameFQD.Data(),"UPDATE");
    if(fileFQD) fileFQD->GetObject("cobjFQD",listHistosFQD);
    if(listHistosFQD) fqd->GetOutputHistograms(listHistosFQD);
    fqd->Finish(); 
    fqdCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listHistosFQD->FindObject("AliFlowCommonHistResultsFQD")); 
    if(fqdCommonHistRes) 
    {
     meanValueFQD->SetPoint(1,binFQD-0.5,(fqdCommonHistRes->GetHistIntFlow())->GetBinContent(1));
     meanValueFQD->SetMarkerStyle(25);
     meanValueFQD->SetMarkerColor(kOrange+7); 
     meanValueFQD->SetMarkerSize(1.25); 
     meanValueFQD->Draw("PSAME");
    } 
    fileFQD->Close();
    fileFQD = TFile::Open(mergedFileNameFQD.Data(),"RECREATE");
    listHistosFQD->SetName("cobjFQD");
    listHistosFQD->Write(listHistosFQD->GetName(), TObject::kSingleKey);
   }
  
 /*   
 // LYZEP:
 TString mergedFileNameLYZEP("outputLYZEPanalysis.root");
 TFile* fileLYZEP = NULL;
 TList *listHistosLYZEP = NULL;
 AliFlowAnalysisWithLYZEventPlane* LYZEP = new AliFlowAnalysisWithLYZEventPlane();
 pwd = gSystem->pwd();
 
 if(gSystem->AccessPathName(((pwd+="/")+=mergedFileNameLYZEP).Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file for LYZEP !!!!"<<endl;
 } else
   { 
    fileLYZEP = TFile::Open(mergedFileNameLYZEP.Data(),"UPDATE");
    if(fileLYZEP) fileLYZEP->GetObject("cobjLYZEP",listHistosLYZEP);
    if(listHistosLYZEP) LYZEP->GetOutputHistograms(listHistosLYZEP);
    LYZEP->Finish(); 
    fileLYZEP->Close();
    fileLYZEP = TFile::Open(mergedFileNameLYZEP.Data(),"RECREATE");
    listHistosLYZEP->SetName("cobjLYZEP");
    listHistosLYZEP->Write(listHistosLYZEP->GetName(), TObject::kSingleKey);
   }
 */  
  
}



















void LoadSpreadLibraries(const libModes mode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libXMLIO.so");
  gSystem->Load("libPhysics.so");
  
  //----------------------------------------------------------
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  if (mode==mLocal) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------

  //==================================================================================  
  //load needed libraries:
  gSystem->AddIncludePath("-I$ROOTSYS/include");
  gSystem->Load("libTree.so");

  // for AliRoot
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libPWG2flowCommon.so");
  cerr<<"libPWG2flowCommon.so loaded ..."<<endl;
  
  }
  
  else if (mode==mLocalSource) {
 
    // In root inline compile
  
    // Constants  
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonConstants.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZConstants.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowCumuConstants.cxx+");
    
    // Flow event
    gROOT->LoadMacro("AliFlowCommon/AliFlowVector.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimple.cxx+");    
    gROOT->LoadMacro("AliFlowCommon/AliFlowEventSimple.cxx+");
    
    // Cuts
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimpleCuts.cxx+");    
    
    // Output histosgrams
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHist.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHistResults.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist1.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist2.cxx+");
       
    cout << "finished loading macros!" << endl;  
    
  }  
  
  
}


