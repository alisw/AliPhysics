//type of analysis can be: ESD, AOD, MC, ESDMC0, ESDMC1
//const TString type = "ESD"; 

enum libModes {mLocal,mLocalSource};
//mLocal: Analyze data on your computer using aliroot
//mLocalSource: Analyze data on your computer using root + source files

//void compareFlowResults(TString type="",Int_t mode=mLocalSource)
void showSpread(const Int_t nRuns=10, TString type="",Int_t mode=mLocal)
{ 
 // load needed libraries:                       
 LoadPlotLibraries(mode);  
 
 // standard magic:
 TString execDir(gSystem->pwd());  
 TSystemDirectory* baseDir = new TSystemDirectory(".",execDir.Data());          
 TList* dirList = baseDir->GetListOfFiles();
 Int_t nDirs = dirList->GetEntries();
 gSystem->cd(execDir);
 
 // array to store estimates of each method from different runs:
 
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
 
 
 

 
 
 //removing the title and stat. box from all histograms:
 //gStyle->SetOptTitle(0);
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
  intFlowNONAME->SetTitle("Integrated Flow NONAME");
  (intFlowNONAME->GetYaxis())->SetRangeUser(0,0.20);
  intFlowNONAME->Draw();
 }
 
 // MCEP 
 if(mcepBoxNONAME) mcepBoxNONAME->Draw("LFSAME");
 if(mcepTGraphNONAME) mcepTGraphNONAME->Draw("PSAME");
 
 // SP 
 if(spBoxNONAME) spBoxNONAME->Draw("LFSAME");
 if(spTGraphNONAME) spTGraphNONAME->Draw("PSAME");
 
 // GFC{2} 
 if(gfc2BoxNONAME) gfc2BoxNONAME->Draw("LFSAME");
 if(gfc2TGraphNONAME) gfc2TGraphNONAME->Draw("PSAME");
  
 // GFC{4} 
 if(gfc4BoxNONAME) gfc4BoxNONAME->Draw("LFSAME");
 if(gfc4TGraphNONAME) gfc4TGraphNONAME->Draw("PSAME");

 // GFC{6} 
 if(gfc6BoxNONAME) gfc6BoxNONAME->Draw("LFSAME");
 if(gfc6TGraphNONAME) gfc6TGraphNONAME->Draw("PSAME");

 // GFC{8} 
 if(gfc8BoxNONAME) gfc8BoxNONAME->Draw("LFSAME");
 if(gfc8TGraphNONAME) gfc8TGraphNONAME->Draw("PSAME");
 
 // QC{2} 
 if(qc2BoxNONAME) qc2BoxNONAME->Draw("LFSAME");
 if(qc2TGraphNONAME) qc2TGraphNONAME->Draw("PSAME");
  
 // QC{4} 
 if(qc4BoxNONAME) qc4BoxNONAME->Draw("LFSAME");
 if(qc4TGraphNONAME) qc4TGraphNONAME->Draw("PSAME");

 // QC{6} 
 if(qc6BoxNONAME) qc6BoxNONAME->Draw("LFSAME");
 if(qc6TGraphNONAME) qc6TGraphNONAME->Draw("PSAME");

 // QC{8} 
 if(qc8BoxNONAME) qc8BoxNONAME->Draw("LFSAME");
 if(qc8TGraphNONAME) qc8TGraphNONAME->Draw("PSAME");
 
 // FQD 
 if(fqdBoxNONAME) fqdBoxNONAME->Draw("LFSAME");
 if(fqdTGraphNONAME) fqdTGraphNONAME->Draw("PSAME");
 
 // LYZ1 
 if(lyz1BoxNONAME) lyz1BoxNONAME->Draw("LFSAME");
 if(lyz1TGraphNONAME) lyz1TGraphNONAME->Draw("PSAME");
 
 // LYZEP 
 //if(lyzepBoxNONAME) lyzepBoxNONAME->Draw("LFSAME");
 //if(lyzepTGraphNONAME) lyzepTGraphNONAME->Draw("PSAME");
  
 (intFlowCanvasNONAME->cd(2))->SetPad(0.75,0.0,1.0,1.0);

}


void LoadPlotLibraries(const libModes mode) {
  
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


