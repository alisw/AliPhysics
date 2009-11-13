enum libModes {mLocal,mLocalSource};

void showSpread(TString type="", Int_t mode=mLocal)
{
 // type:  type of analysis can be ESD, AOD, MC, ESDMC0, ESDMC1
 //        (if type="" output files are from MC simulation (default))
 // mode:  if mode=mLocal analyze data on your computer using aliroot
 //        if mode=mLocalSource analyze data on your computer using root + source files
 
 Bool_t plotLegend = kTRUE; // plot legend on the plot
 
 // load needed libraries:                       
 LoadSpreadLibraries(mode);  

 // access the path of current directory:
 TString pwd(gSystem->pwd());
 pwd+="/"; 
 
 // standard magic:
 TString execDir(gSystem->pwd());  
 TSystemDirectory* baseDir = new TSystemDirectory(".",execDir.Data());          
 TList* dirList = baseDir->GetListOfFiles();
 Int_t nDirs = dirList->GetEntries();
 gSystem->cd(execDir);
 
 Int_t countDirectories = 0; // count all directories except the ones with merged results
 Int_t countMergedSubsets = 0; // count only directories with merged results (name is subset*)
 if(dirList)
 {
  for(Int_t i=0;i<nDirs;i++)
  {
   if((dirList->At(i))->InheritsFrom("TSystemDirectory") &&  
      !(strcmp((dirList->At(i))->GetName(),".") == 0 ) &&  
      !(strcmp((dirList->At(i))->GetName(),"..") == 0 )) 
   {
    TSystemFile *presentFile = (TSystemFile*)dirList->At(i);
    TString presentFileName(gSystem->pwd());
    (presentFileName+="/")+=presentFile->GetName(); 
    if(presentFileName.Contains("subset"))
    {
     countMergedSubsets++;
    } else 
      {
       countDirectories++;
      } 
   } 
  } // end of for(Int_t i=0;i<nDirs;i++)
 } // end of if(dirList)  
 
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
 Int_t binLYZ1SUM = 12; 
 Int_t binLYZ1PROD = 13;
 Int_t binLYZEP = 14;
 
 // option A: one subdirectory gives one estimate for each method (merged subsets ignored): 
 const Int_t nEstimates = countDirectories;
 
 // option B: one merged subset give one estimate for each method
 const Int_t nSubsets = countMergedSubsets;

 if(nSubsets == 0)
 {
  // arrays to store estimates and errors of each method from different small statistics runs:
  // MCEP:
  Double_t mcepValueNONAME[nEstimates] = {0.}; 
  Double_t mcepErrorNONAME[nEstimates] = {0.}; 
  Double_t mcepMaxValueNONAME = 0.; // to be improved      
  Double_t mcepMinValueNONAME = 1000.; // to be improved
 
  // SP:
  Double_t spValueNONAME[nEstimates] = {0.}; 
  Double_t spErrorNONAME[nEstimates] = {0.}; 
  Double_t spMaxValueNONAME = 0.; // to be improved      
  Double_t spMinValueNONAME = 1000.; // to be improved
 
  // GFC{2}:
  Double_t gfc2ValueNONAME[nEstimates] = {0.}; 
  Double_t gfc2ErrorNONAME[nEstimates] = {0.}; 
  Double_t gfc2MaxValueNONAME = 0.; // to be improved      
  Double_t gfc2MinValueNONAME = 1000.; // to be improved
 
  // GFC{4}:
  Double_t gfc4ValueNONAME[nEstimates] = {0.}; 
  Double_t gfc4ErrorNONAME[nEstimates] = {0.}; 
  Double_t gfc4MaxValueNONAME = 0.; // to be improved      
  Double_t gfc4MinValueNONAME = 1000.; // to be improved
    
  // GFC{6}:
  Double_t gfc6ValueNONAME[nEstimates] = {0.}; 
  Double_t gfc6ErrorNONAME[nEstimates] = {0.}; 
  Double_t gfc6MaxValueNONAME = 0.; // to be improved      
  Double_t gfc6MinValueNONAME = 1000.; // to be improved
 
  // GFC{8}:
  Double_t gfc8ValueNONAME[nEstimates] = {0.}; 
  Double_t gfc8ErrorNONAME[nEstimates] = {0.}; 
  Double_t gfc8MaxValueNONAME = 0.; // to be improved      
  Double_t gfc8MinValueNONAME = 1000.; // to be improved
 
  // QC{2}:
  Double_t qc2ValueNONAME[nEstimates] = {0.}; 
  Double_t qc2ErrorNONAME[nEstimates] = {0.}; 
  Double_t qc2MaxValueNONAME = 0.; // to be improved      
  Double_t qc2MinValueNONAME = 1000.; // to be improved
 
  // QC{4}:
  Double_t qc4ValueNONAME[nEstimates] = {0.}; 
  Double_t qc4ErrorNONAME[nEstimates] = {0.}; 
  Double_t qc4MaxValueNONAME = 0.; // to be improved      
  Double_t qc4MinValueNONAME = 1000.; // to be improved
    
  // QC{6}:
  Double_t qc6ValueNONAME[nEstimates] = {0.}; 
  Double_t qc6ErrorNONAME[nEstimates] = {0.}; 
  Double_t qc6MaxValueNONAME = 0.; // to be improved      
  Double_t qc6MinValueNONAME = 1000.; // to be improved

  // QC{8}:
  Double_t qc8ValueNONAME[nEstimates] = {0.}; 
  Double_t qc8ErrorNONAME[nEstimates] = {0.}; 
  Double_t qc8MaxValueNONAME = 0.; // to be improved      
  Double_t qc8MinValueNONAME = 1000.; // to be improved
 
  // FQD:
  Double_t fqdValueNONAME[nEstimates] = {0.}; 
  Double_t fqdErrorNONAME[nEstimates] = {0.}; 
  Double_t fqdMaxValueNONAME = 0.; // to be improved      
  Double_t fqdMinValueNONAME = 1000.; // to be improved
       
  // LYZ1SUM:
  Double_t lyz1sumValueNONAME[nEstimates] = {0.}; 
  Double_t lyz1sumErrorNONAME[nEstimates] = {0.}; 
  Double_t lyz1sumMaxValueNONAME = 0.; // to be improved      
  Double_t lyz1sumMinValueNONAME = 1000.; // to be improved
 
  // LYZ1PROD:
  Double_t lyz1prodValueNONAME[nEstimates] = {0.}; 
  Double_t lyz1prodErrorNONAME[nEstimates] = {0.}; 
  Double_t lyz1prodMaxValueNONAME = 0.; // to be improved      
  Double_t lyz1prodMinValueNONAME = 1000.; // to be improved
             
  // LYZEP:
  Double_t lyzepValueNONAME[nEstimates] = {0.}; 
  Double_t lyzepErrorNONAME[nEstimates] = {0.}; 
  Double_t lyzepMaxValueNONAME = 0.; // to be improved      
  Double_t lyzepMinValueNONAME = 1000.; // to be improved        
             
  Int_t countFilesMCEP = 0;
  Int_t countFilesSP = 0;
  Int_t countFilesGFC = 0;
  Int_t countFilesQC = 0;
  Int_t countFilesFQD = 0;
  Int_t countFilesLYZ1SUM = 0;
  Int_t countFilesLYZ1PROD = 0;
  Int_t countFilesLYZEP = 0;
  
  for(Int_t iDir=0;iDir<nDirs;++iDir)
  { 
   TSystemFile* presentDir = (TSystemFile*)dirList->At(iDir);
   if(!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || 
     strcmp(presentDir->GetName(), "..") == 0) continue; 
                                     
   TString presentDirName(gSystem->pwd()); 
   presentDirName += "/";
   presentDirName += presentDir->GetName();
   presentDirName += "/";
    
   // accessing the small statistics output .root files for each method:
   // MCEP:     
   TString fileNameMCEP = presentDirName;   
   fileNameMCEP+="outputMCEPanalysis";
   (fileNameMCEP+=type.Data())+=".root";
   TFile *fileMCEP = TFile::Open(fileNameMCEP.Data(), "READ");      
   if(fileMCEP) 
   {
    TList *listMCEP = NULL;
    AliFlowCommonHistResults *mcepCommonHistRes = NULL; 
    fileMCEP->GetObject("cobjMCEP",listMCEP); 
    fileMCEP->Close();
    if(listMCEP) 
    {
     mcepCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listMCEP->FindObject("AliFlowCommonHistResultsMCEP")); 
     if(mcepCommonHistRes && mcepCommonHistRes->GetHistIntFlow())
     {
      mcepValueNONAME[countFilesMCEP] = (mcepCommonHistRes->GetHistIntFlow())->GetBinContent(1);
      mcepErrorNONAME[countFilesMCEP] = (mcepCommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(mcepValueNONAME[countFilesMCEP]>=0.) // to be improved 
      {
       if(mcepMaxValueNONAME < mcepValueNONAME[countFilesMCEP]) mcepMaxValueNONAME = mcepValueNONAME[countFilesMCEP]; 
       if(mcepMinValueNONAME > mcepValueNONAME[countFilesMCEP]) mcepMinValueNONAME = mcepValueNONAME[countFilesMCEP]; 
      } 
     }
     countFilesMCEP++;
    } // end of if(listMCEP)
   } // end of if(fileMCEP) 
  
   // SP:     
   TString fileNameSP = presentDirName;   
   fileNameSP+="outputSPanalysis";
   (fileNameSP+=type.Data())+=".root";
   TFile *fileSP = TFile::Open(fileNameSP.Data(), "READ");      
   if(fileSP) 
   {
    TList *listSP = NULL;
    AliFlowCommonHistResults *spCommonHistRes = NULL; 
    fileSP->GetObject("cobjSP",listSP); 
    fileSP->Close();
    if(listSP) 
    {
     spCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listSP->FindObject("AliFlowCommonHistResultsSP")); 
     if(spCommonHistRes && spCommonHistRes->GetHistIntFlow())
     {
      spValueNONAME[countFilesSP] = (spCommonHistRes->GetHistIntFlow())->GetBinContent(1);
      spErrorNONAME[countFilesSP] = (spCommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(spValueNONAME[countFilesSP]>=0.) // to be improved 
      {
       if(spMaxValueNONAME < spValueNONAME[countFilesSP]) spMaxValueNONAME = spValueNONAME[countFilesSP]; 
       if(spMinValueNONAME > spValueNONAME[countFilesSP]) spMinValueNONAME = spValueNONAME[countFilesSP]; 
      } 
     }
     countFilesSP++;
    } // end of if(listSP)
   } // end of if(fileSP) 
   
   // GFC:     
   TString fileNameGFC = presentDirName;   
   fileNameGFC+="outputGFCanalysis";
   (fileNameGFC+=type.Data())+=".root";
   TFile *fileGFC = TFile::Open(fileNameGFC.Data(), "READ");      
   if(fileGFC) 
   {
    TList *listGFC = NULL;
    AliFlowCommonHistResults *gfc2CommonHistRes = NULL; 
    AliFlowCommonHistResults *gfc4CommonHistRes = NULL; 
    AliFlowCommonHistResults *gfc6CommonHistRes = NULL; 
    AliFlowCommonHistResults *gfc8CommonHistRes = NULL; 
    fileGFC->GetObject("cobjGFC",listGFC); 
    fileGFC->Close();
    if(listGFC) 
    {
     gfc2CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults2ndOrderGFC"));
     gfc4CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults4thOrderGFC"));
     gfc6CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults6thOrderGFC"));
     gfc8CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults8thOrderGFC")); 
     if(gfc2CommonHistRes && gfc2CommonHistRes->GetHistIntFlow())
     {
      gfc2ValueNONAME[countFilesGFC] = (gfc2CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      gfc2ErrorNONAME[countFilesGFC] = (gfc2CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(gfc2ValueNONAME[countFilesGFC]>=0.) // to be improved 
      {
       if(gfc2MaxValueNONAME < gfc2ValueNONAME[countFilesGFC]) gfc2MaxValueNONAME = gfc2ValueNONAME[countFilesGFC]; 
       if(gfc2MinValueNONAME > gfc2ValueNONAME[countFilesGFC]) gfc2MinValueNONAME = gfc2ValueNONAME[countFilesGFC]; 
      } 
     }
     if(gfc4CommonHistRes && gfc4CommonHistRes->GetHistIntFlow())
     {
      gfc4ValueNONAME[countFilesGFC] = (gfc4CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      gfc4ErrorNONAME[countFilesGFC] = (gfc4CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(gfc4ValueNONAME[countFilesGFC]>=0.) // to be improved 
      {
       if(gfc4MaxValueNONAME < gfc4ValueNONAME[countFilesGFC]) gfc4MaxValueNONAME = gfc4ValueNONAME[countFilesGFC]; 
       if(gfc4MinValueNONAME > gfc4ValueNONAME[countFilesGFC]) gfc4MinValueNONAME = gfc4ValueNONAME[countFilesGFC]; 
      } 
     }
     if(gfc6CommonHistRes && gfc6CommonHistRes->GetHistIntFlow())
     {
      gfc6ValueNONAME[countFilesGFC] = (gfc6CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      gfc6ErrorNONAME[countFilesGFC] = (gfc6CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(gfc6ValueNONAME[countFilesGFC]>=0.) // to be improved 
      {
       if(gfc6MaxValueNONAME < gfc6ValueNONAME[countFilesGFC]) gfc6MaxValueNONAME = gfc6ValueNONAME[countFilesGFC]; 
       if(gfc6MinValueNONAME > gfc6ValueNONAME[countFilesGFC]) gfc6MinValueNONAME = gfc6ValueNONAME[countFilesGFC]; 
      } 
     }
     if(gfc8CommonHistRes && gfc8CommonHistRes->GetHistIntFlow())
     {
      gfc8ValueNONAME[countFilesGFC] = (gfc8CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      gfc8ErrorNONAME[countFilesGFC] = (gfc8CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(gfc8ValueNONAME[countFilesGFC]>=0.) // to be improved 
      {
       if(gfc8MaxValueNONAME < gfc8ValueNONAME[countFilesGFC]) gfc8MaxValueNONAME = gfc8ValueNONAME[countFilesGFC]; 
       if(gfc8MinValueNONAME > gfc8ValueNONAME[countFilesGFC]) gfc8MinValueNONAME = gfc8ValueNONAME[countFilesGFC]; 
      } 
     }
     countFilesGFC++;
    } // end of if(listGFC)
   } // end of if(fileGFC) 
   
   // QC:     
   TString fileNameQC = presentDirName;   
   fileNameQC+="outputQCanalysis";
   (fileNameQC+=type.Data())+=".root";
   TFile *fileQC = TFile::Open(fileNameQC.Data(), "READ");      
   if(fileQC) 
   {
    TList *listQC = NULL;
    AliFlowCommonHistResults *qc2CommonHistRes = NULL; 
    AliFlowCommonHistResults *qc4CommonHistRes = NULL;     
    AliFlowCommonHistResults *qc6CommonHistRes = NULL; 
    AliFlowCommonHistResults *qc8CommonHistRes = NULL; 
    fileQC->GetObject("cobjQC",listQC); 
    fileQC->Close();
    if(listQC) 
    {
     qc2CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults2ndOrderQC"));
     qc4CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults4thOrderQC"));
     qc6CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults6thOrderQC"));
     qc8CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults8thOrderQC")); 
     if(qc2CommonHistRes && qc2CommonHistRes->GetHistIntFlow())
     {
      qc2ValueNONAME[countFilesQC] = (qc2CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      qc2ErrorNONAME[countFilesQC] = (qc2CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(qc2ValueNONAME[countFilesQC]>=0.) // to be improved 
      {
       if(qc2MaxValueNONAME < qc2ValueNONAME[countFilesQC]) qc2MaxValueNONAME = qc2ValueNONAME[countFilesQC]; 
       if(qc2MinValueNONAME > qc2ValueNONAME[countFilesQC]) qc2MinValueNONAME = qc2ValueNONAME[countFilesQC]; 
      } 
     }
     if(qc4CommonHistRes && qc4CommonHistRes->GetHistIntFlow())
     {
      qc4ValueNONAME[countFilesQC] = (qc4CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      qc4ErrorNONAME[countFilesQC] = (qc4CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(qc4ValueNONAME[countFilesQC]>=0.) // to be improved 
      {
       if(qc4MaxValueNONAME < qc4ValueNONAME[countFilesQC]) qc4MaxValueNONAME = qc4ValueNONAME[countFilesQC]; 
       if(qc4MinValueNONAME > qc4ValueNONAME[countFilesQC]) qc4MinValueNONAME = qc4ValueNONAME[countFilesQC]; 
      } 
     }
     if(qc6CommonHistRes && qc6CommonHistRes->GetHistIntFlow())
     {
      qc6ValueNONAME[countFilesQC] = (qc6CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      qc6ErrorNONAME[countFilesQC] = (qc6CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(qc6ValueNONAME[countFilesQC]>=0.) // to be improved 
      {
       if(qc6MaxValueNONAME < qc6ValueNONAME[countFilesQC]) qc6MaxValueNONAME = qc6ValueNONAME[countFilesQC]; 
       if(qc6MinValueNONAME > qc6ValueNONAME[countFilesQC]) qc6MinValueNONAME = qc6ValueNONAME[countFilesQC]; 
      } 
     }
     if(qc8CommonHistRes && qc8CommonHistRes->GetHistIntFlow())
     {
      qc8ValueNONAME[countFilesQC] = (qc8CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      qc8ErrorNONAME[countFilesQC] = (qc8CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(qc8ValueNONAME[countFilesQC]>=0.) // to be improved 
      {
       if(qc8MaxValueNONAME < qc8ValueNONAME[countFilesQC]) qc8MaxValueNONAME = qc8ValueNONAME[countFilesQC]; 
       if(qc8MinValueNONAME > qc8ValueNONAME[countFilesQC]) qc8MinValueNONAME = qc8ValueNONAME[countFilesQC]; 
      } 
     }
     countFilesQC++;
    } // end of if(listQC)
   } // end of if(fileQC) 
  
   // FQD:     
   TString fileNameFQD = presentDirName;   
   fileNameFQD+="outputFQDanalysis";
   (fileNameFQD+=type.Data())+=".root";
   TFile *fileFQD = TFile::Open(fileNameFQD.Data(), "READ");      
   if(fileFQD) 
   {
    TList *listFQD = NULL;
    AliFlowCommonHistResults *fqdCommonHistRes = NULL; 
    fileFQD->GetObject("cobjFQD",listFQD); 
    fileFQD->Close();
    if(listFQD) 
    {
     fqdCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listFQD->FindObject("AliFlowCommonHistResultsFQD")); 
     if(fqdCommonHistRes && fqdCommonHistRes->GetHistIntFlow())
     {
      fqdValueNONAME[countFilesFQD] = (fqdCommonHistRes->GetHistIntFlow())->GetBinContent(1);
      fqdErrorNONAME[countFilesFQD] = (fqdCommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(fqdValueNONAME[countFilesFQD]>=0.) // to be improved 
      {
       if(fqdMaxValueNONAME < fqdValueNONAME[countFilesFQD]) fqdMaxValueNONAME = fqdValueNONAME[countFilesFQD]; 
       if(fqdMinValueNONAME > fqdValueNONAME[countFilesFQD]) fqdMinValueNONAME = fqdValueNONAME[countFilesFQD]; 
      } 
     }
     countFilesFQD++;
    } // end of if(listFQD)
   } // end of if(fileFQD)   
  
   // LYZ1SUM:     
   TString fileNameLYZ1SUM = presentDirName;   
   fileNameLYZ1SUM+="outputLYZ1SUManalysis";
   (fileNameLYZ1SUM+=type.Data())+=".root";
   TFile *fileLYZ1SUM = TFile::Open(fileNameLYZ1SUM.Data(), "READ");      
   if(fileLYZ1SUM) 
   {
    TList *listLYZ1SUM = NULL;
    AliFlowCommonHistResults *lyz1sumCommonHistRes = NULL; 
    fileLYZ1SUM->GetObject("cobjLYZ1SUM",listLYZ1SUM); 
    fileLYZ1SUM->Close();
    if(listLYZ1SUM) 
    {
     lyz1sumCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listLYZ1SUM->FindObject("AliFlowCommonHistResultsLYZ1SUM")); 
     if(lyz1sumCommonHistRes && lyz1sumCommonHistRes->GetHistIntFlow())
     {
      lyz1sumValueNONAME[countFilesLYZ1SUM] = (lyz1sumCommonHistRes->GetHistIntFlow())->GetBinContent(1);
      lyz1sumErrorNONAME[countFilesLYZ1SUM] = (lyz1sumCommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(lyz1sumValueNONAME[countFilesLYZ1SUM]>=0.) // to be improved 
      {
       if(lyz1sumMaxValueNONAME < lyz1sumValueNONAME[countFilesLYZ1SUM]) lyz1sumMaxValueNONAME = lyz1sumValueNONAME[countFilesLYZ1SUM]; 
       if(lyz1sumMinValueNONAME > lyz1sumValueNONAME[countFilesLYZ1SUM]) lyz1sumMinValueNONAME = lyz1sumValueNONAME[countFilesLYZ1SUM]; 
      } 
     }
     countFilesLYZ1SUM++;
    } // end of if(listLYZ1SUM)
   } // end of if(fileLYZ1SUM)   
  
   // LYZ1PROD:     
   TString fileNameLYZ1PROD = presentDirName;   
   fileNameLYZ1PROD+="outputLYZ1PRODanalysis";
   (fileNameLYZ1PROD+=type.Data())+=".root";
   TFile *fileLYZ1PROD = TFile::Open(fileNameLYZ1PROD.Data(), "READ");      
   if(fileLYZ1PROD) 
   {
    TList *listLYZ1PROD = NULL;
    AliFlowCommonHistResults *lyz1prodCommonHistRes = NULL; 
    fileLYZ1PROD->GetObject("cobjLYZ1PROD",listLYZ1PROD); 
    fileLYZ1PROD->Close();
    if(listLYZ1PROD) 
    {
     lyz1prodCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listLYZ1PROD->FindObject("AliFlowCommonHistResultsLYZ1PROD")); 
     if(lyz1prodCommonHistRes && lyz1prodCommonHistRes->GetHistIntFlow())
     {
      lyz1prodValueNONAME[countFilesLYZ1PROD] = (lyz1prodCommonHistRes->GetHistIntFlow())->GetBinContent(1); 
      lyz1prodErrorNONAME[countFilesLYZ1PROD] = (lyz1prodCommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(lyz1prodValueNONAME[countFilesLYZ1PROD]>=0.) // to be improved 
      {
       if(lyz1prodMaxValueNONAME < lyz1prodValueNONAME[countFilesLYZ1PROD]) lyz1prodMaxValueNONAME = lyz1prodValueNONAME[countFilesLYZ1PROD]; 
       if(lyz1prodMinValueNONAME > lyz1prodValueNONAME[countFilesLYZ1PROD]) lyz1prodMinValueNONAME = lyz1prodValueNONAME[countFilesLYZ1PROD]; 
      } 
     }
     countFilesLYZ1PROD++;
    } // end of if(listLYZ1PROD)
   } // end of if(fileLYZ1PROD)   
   
  } // end of for(Int_t iDir=0;iDir<nDirs;++iDir)
 } else if(nSubsets > 0)
 {
  // arrays to store estimates and errors of each method from merged subsets:
  // MCEP:
  Double_t mcepValueNONAME[nSubsets] = {0.}; 
  Double_t mcepErrorNONAME[nSubsets] = {0.}; 
  Double_t mcepMaxValueNONAME = 0.; // to be improved      
  Double_t mcepMinValueNONAME = 1000.; // to be improved
 
  // SP:
  Double_t spValueNONAME[nSubsets] = {0.}; 
  Double_t spErrorNONAME[nSubsets] = {0.}; 
  Double_t spMaxValueNONAME = 0.; // to be improved      
  Double_t spMinValueNONAME = 1000.; // to be improved
 
  // GFC{2}:
  Double_t gfc2ValueNONAME[nSubsets] = {0.}; 
  Double_t gfc2ErrorNONAME[nSubsets] = {0.}; 
  Double_t gfc2MaxValueNONAME = 0.; // to be improved      
  Double_t gfc2MinValueNONAME = 1000.; // to be improved
 
  // GFC{4}:
  Double_t gfc4ValueNONAME[nSubsets] = {0.}; 
  Double_t gfc4ErrorNONAME[nSubsets] = {0.}; 
  Double_t gfc4MaxValueNONAME = 0.; // to be improved      
  Double_t gfc4MinValueNONAME = 1000.; // to be improved
    
  // GFC{6}:
  Double_t gfc6ValueNONAME[nSubsets] = {0.}; 
  Double_t gfc6ErrorNONAME[nSubsets] = {0.}; 
  Double_t gfc6MaxValueNONAME = 0.; // to be improved      
  Double_t gfc6MinValueNONAME = 1000.; // to be improved
 
  // GFC{8}:
  Double_t gfc8ValueNONAME[nSubsets] = {0.}; 
  Double_t gfc8ErrorNONAME[nSubsets] = {0.}; 
  Double_t gfc8MaxValueNONAME = 0.; // to be improved      
  Double_t gfc8MinValueNONAME = 1000.; // to be improved
 
  // QC{2}:
  Double_t qc2ValueNONAME[nSubsets] = {0.}; 
  Double_t qc2ErrorNONAME[nSubsets] = {0.}; 
  Double_t qc2MaxValueNONAME = 0.; // to be improved      
  Double_t qc2MinValueNONAME = 1000.; // to be improved
 
  // QC{4}:
  Double_t qc4ValueNONAME[nSubsets] = {0.}; 
  Double_t qc4ErrorNONAME[nSubsets] = {0.}; 
  Double_t qc4MaxValueNONAME = 0.; // to be improved      
  Double_t qc4MinValueNONAME = 1000.; // to be improved
    
  // QC{6}:
  Double_t qc6ValueNONAME[nSubsets] = {0.}; 
  Double_t qc6ErrorNONAME[nSubsets] = {0.}; 
  Double_t qc6MaxValueNONAME = 0.; // to be improved      
  Double_t qc6MinValueNONAME = 1000.; // to be improved

  // QC{8}:
  Double_t qc8ValueNONAME[nSubsets] = {0.}; 
  Double_t qc8ErrorNONAME[nSubsets] = {0.}; 
  Double_t qc8MaxValueNONAME = 0.; // to be improved      
  Double_t qc8MinValueNONAME = 1000.; // to be improved
 
  // FQD:
  Double_t fqdValueNONAME[nSubsets] = {0.}; 
  Double_t fqdErrorNONAME[nSubsets] = {0.}; 
  Double_t fqdMaxValueNONAME = 0.; // to be improved      
  Double_t fqdMinValueNONAME = 1000.; // to be improved
       
  // LYZ1SUM:
  Double_t lyz1sumValueNONAME[nSubsets] = {0.}; 
  Double_t lyz1sumErrorNONAME[nSubsets] = {0.}; 
  Double_t lyz1sumMaxValueNONAME = 0.; // to be improved      
  Double_t lyz1sumMinValueNONAME = 1000.; // to be improved
 
  // LYZ1PROD:
  Double_t lyz1prodValueNONAME[nSubsets] = {0.}; 
  Double_t lyz1prodErrorNONAME[nSubsets] = {0.}; 
  Double_t lyz1prodMaxValueNONAME = 0.; // to be improved      
  Double_t lyz1prodMinValueNONAME = 1000.; // to be improved
             
  // LYZEP:
  Double_t lyzepValueNONAME[nSubsets] = {0.}; 
  Double_t lyzepErrorNONAME[nSubsets] = {0.}; 
  Double_t lyzepMaxValueNONAME = 0.; // to be improved      
  Double_t lyzepMinValueNONAME = 1000.; // to be improved        
             
  Int_t countFilesMCEP = 0;
  Int_t countFilesSP = 0;
  Int_t countFilesGFC = 0;
  Int_t countFilesQC = 0;
  Int_t countFilesFQD = 0;
  Int_t countFilesLYZ1SUM = 0;
  Int_t countFilesLYZ1PROD = 0;          
  Int_t countFilesLYZEP = 0; 
 
  for(Int_t iDir=0;iDir<nDirs;iDir++)
  { 
   TSystemFile* presentDir = (TSystemFile*)dirList->At(iDir);
   if(!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || 
     strcmp(presentDir->GetName(), "..") == 0) continue; 
                                                         
   TString presentDirName(gSystem->pwd()); 
   presentDirName += "/";
   presentDirName += presentDir->GetName();
   presentDirName += "/";
   
   if(!(presentDirName.Contains("subset"))) continue;
  
   // accessing the output *.root files from merged subsets for each method:
   // MCEP:     
   TString fileNameMCEP = gSystem->pwd();
   fileNameMCEP+="/";
   fileNameMCEP+=presentDir->GetName();
   fileNameMCEP+="/";
   fileNameMCEP+="outputMCEPanalysis";
   (fileNameMCEP+=type.Data())+=".root";
   TFile *fileMCEP = TFile::Open(fileNameMCEP.Data(), "READ");      
   if(fileMCEP) 
   {
    TList *listMCEP = NULL;
    AliFlowCommonHistResults *mcepCommonHistRes = NULL; 
    fileMCEP->GetObject("cobjMCEP",listMCEP); 
    fileMCEP->Close();
    if(listMCEP) 
    {
     mcepCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listMCEP->FindObject("AliFlowCommonHistResultsMCEP")); 
     if(mcepCommonHistRes && mcepCommonHistRes->GetHistIntFlow())
     {
      mcepValueNONAME[countFilesMCEP] = (mcepCommonHistRes->GetHistIntFlow())->GetBinContent(1);
      mcepErrorNONAME[countFilesMCEP] = (mcepCommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(mcepValueNONAME[countFilesMCEP]>=0.) // to be improved 
      {
       if(mcepMaxValueNONAME < mcepValueNONAME[countFilesMCEP]) mcepMaxValueNONAME = mcepValueNONAME[countFilesMCEP]; 
       if(mcepMinValueNONAME > mcepValueNONAME[countFilesMCEP]) mcepMinValueNONAME = mcepValueNONAME[countFilesMCEP]; 
      } 
       countFilesMCEP++;
     }
    } // end of if(listMCEP)
   } // end of if(fileMCEP) 
  
   // SP:     
   TString fileNameSP = gSystem->pwd();
   fileNameSP+="/";
   fileNameSP+=presentDir->GetName();
   fileNameSP+="/";
   fileNameSP+="outputSPanalysis";
   (fileNameSP+=type.Data())+=".root";
   TFile *fileSP = TFile::Open(fileNameSP.Data(), "READ");      
   if(fileSP) 
   {
    TList *listSP = NULL;
    AliFlowCommonHistResults *spCommonHistRes = NULL; 
    fileSP->GetObject("cobjSP",listSP); 
    fileSP->Close();
    if(listSP) 
    {
     spCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listSP->FindObject("AliFlowCommonHistResultsSP")); 
     if(spCommonHistRes && spCommonHistRes->GetHistIntFlow())
     {
      spValueNONAME[countFilesSP] = (spCommonHistRes->GetHistIntFlow())->GetBinContent(1);
      spErrorNONAME[countFilesSP] = (spCommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(spValueNONAME[countFilesSP]>=0.) // to be improved 
      {
       if(spMaxValueNONAME < spValueNONAME[countFilesSP]) spMaxValueNONAME = spValueNONAME[countFilesSP]; 
       if(spMinValueNONAME > spValueNONAME[countFilesSP]) spMinValueNONAME = spValueNONAME[countFilesSP]; 
      } 
      countFilesSP++;     
     }
    } // end of if(listSP)
   } // end of if(fileSP) 
   
   // GFC:     
   TString fileNameGFC = gSystem->pwd();
   fileNameGFC+="/";
   fileNameGFC+=presentDir->GetName();
   fileNameGFC+="/";
   fileNameGFC+="outputGFCanalysis";
   (fileNameGFC+=type.Data())+=".root";
   TFile *fileGFC = TFile::Open(fileNameGFC.Data(), "READ");      
   if(fileGFC) 
   {
    TList *listGFC = NULL;
    AliFlowCommonHistResults *gfc2CommonHistRes = NULL; 
    AliFlowCommonHistResults *gfc4CommonHistRes = NULL; 
    AliFlowCommonHistResults *gfc6CommonHistRes = NULL; 
    AliFlowCommonHistResults *gfc8CommonHistRes = NULL; 
    fileGFC->GetObject("cobjGFC",listGFC); 
    fileGFC->Close();
    if(listGFC) 
    {
     gfc2CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults2ndOrderGFC"));
     gfc4CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults4thOrderGFC"));
     gfc6CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults6thOrderGFC"));
     gfc8CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listGFC->FindObject("AliFlowCommonHistResults8thOrderGFC")); 
     if(gfc2CommonHistRes && gfc2CommonHistRes->GetHistIntFlow())
     {
      gfc2ValueNONAME[countFilesGFC] = (gfc2CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      gfc2ErrorNONAME[countFilesGFC] = (gfc2CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(gfc2ValueNONAME[countFilesGFC]>=0.) // to be improved 
      {
       if(gfc2MaxValueNONAME < gfc2ValueNONAME[countFilesGFC]) gfc2MaxValueNONAME = gfc2ValueNONAME[countFilesGFC]; 
       if(gfc2MinValueNONAME > gfc2ValueNONAME[countFilesGFC]) gfc2MinValueNONAME = gfc2ValueNONAME[countFilesGFC]; 
      } 
     }
     if(gfc4CommonHistRes && gfc4CommonHistRes->GetHistIntFlow())
     {
      gfc4ValueNONAME[countFilesGFC] = (gfc4CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      gfc4ErrorNONAME[countFilesGFC] = (gfc4CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(gfc4ValueNONAME[countFilesGFC]>=0.) // to be improved 
      {
       if(gfc4MaxValueNONAME < gfc4ValueNONAME[countFilesGFC]) gfc4MaxValueNONAME = gfc4ValueNONAME[countFilesGFC]; 
       if(gfc4MinValueNONAME > gfc4ValueNONAME[countFilesGFC]) gfc4MinValueNONAME = gfc4ValueNONAME[countFilesGFC]; 
      } 
     }
     if(gfc6CommonHistRes && gfc6CommonHistRes->GetHistIntFlow())
     {
      gfc6ValueNONAME[countFilesGFC] = (gfc6CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      gfc6ErrorNONAME[countFilesGFC] = (gfc6CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(gfc6ValueNONAME[countFilesGFC]>=0.) // to be improved 
      {
       if(gfc6MaxValueNONAME < gfc6ValueNONAME[countFilesGFC]) gfc6MaxValueNONAME = gfc6ValueNONAME[countFilesGFC]; 
       if(gfc6MinValueNONAME > gfc6ValueNONAME[countFilesGFC]) gfc6MinValueNONAME = gfc6ValueNONAME[countFilesGFC]; 
      } 
     }
     if(gfc8CommonHistRes && gfc8CommonHistRes->GetHistIntFlow())
     {
      gfc8ValueNONAME[countFilesGFC] = (gfc8CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      gfc8ErrorNONAME[countFilesGFC] = (gfc8CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(gfc8ValueNONAME[countFilesGFC]>=0.) // to be improved 
      {
       if(gfc8MaxValueNONAME < gfc8ValueNONAME[countFilesGFC]) gfc8MaxValueNONAME = gfc8ValueNONAME[countFilesGFC]; 
       if(gfc8MinValueNONAME > gfc8ValueNONAME[countFilesGFC]) gfc8MinValueNONAME = gfc8ValueNONAME[countFilesGFC]; 
      }
      countFilesGFC++; 
     }
    } // end of if(listGFC)
   } // end of if(fileGFC) 
   
   // QC:     
   TString fileNameQC = gSystem->pwd();
   fileNameQC+="/";
   fileNameQC+=presentDir->GetName();
   fileNameQC+="/";
   fileNameQC+="outputQCanalysis";
   (fileNameQC+=type.Data())+=".root";
   TFile *fileQC = TFile::Open(fileNameQC.Data(), "READ");      
   if(fileQC) 
   {
    TList *listQC = NULL;
    AliFlowCommonHistResults *qc2CommonHistRes = NULL; 
    AliFlowCommonHistResults *qc4CommonHistRes = NULL;     
    AliFlowCommonHistResults *qc6CommonHistRes = NULL; 
    AliFlowCommonHistResults *qc8CommonHistRes = NULL; 
    fileQC->GetObject("cobjQC",listQC); 
    fileQC->Close();
    if(listQC) 
    {
     qc2CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults2ndOrderQC"));
     qc4CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults4thOrderQC"));
     qc6CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults6thOrderQC"));
     qc8CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listQC->FindObject("AliFlowCommonHistResults8thOrderQC")); 
     if(qc2CommonHistRes && qc2CommonHistRes->GetHistIntFlow())
     {
      qc2ValueNONAME[countFilesQC] = (qc2CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      qc2ErrorNONAME[countFilesQC] = (qc2CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(qc2ValueNONAME[countFilesQC]>=0.) // to be improved 
      {
       if(qc2MaxValueNONAME < qc2ValueNONAME[countFilesQC]) qc2MaxValueNONAME = qc2ValueNONAME[countFilesQC]; 
       if(qc2MinValueNONAME > qc2ValueNONAME[countFilesQC]) qc2MinValueNONAME = qc2ValueNONAME[countFilesQC]; 
      } 
     }
     if(qc4CommonHistRes && qc4CommonHistRes->GetHistIntFlow())
     {
      qc4ValueNONAME[countFilesQC] = (qc4CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      qc4ErrorNONAME[countFilesQC] = (qc4CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(qc4ValueNONAME[countFilesQC]>=0.) // to be improved 
      {
       if(qc4MaxValueNONAME < qc4ValueNONAME[countFilesQC]) qc4MaxValueNONAME = qc4ValueNONAME[countFilesQC]; 
       if(qc4MinValueNONAME > qc4ValueNONAME[countFilesQC]) qc4MinValueNONAME = qc4ValueNONAME[countFilesQC]; 
      } 
     }
     if(qc6CommonHistRes && qc6CommonHistRes->GetHistIntFlow())
     {
      qc6ValueNONAME[countFilesQC] = (qc6CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      qc6ErrorNONAME[countFilesQC] = (qc6CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(qc6ValueNONAME[countFilesQC]>=0.) // to be improved 
      {
       if(qc6MaxValueNONAME < qc6ValueNONAME[countFilesQC]) qc6MaxValueNONAME = qc6ValueNONAME[countFilesQC]; 
       if(qc6MinValueNONAME > qc6ValueNONAME[countFilesQC]) qc6MinValueNONAME = qc6ValueNONAME[countFilesQC]; 
      } 
     }
     if(qc8CommonHistRes && qc8CommonHistRes->GetHistIntFlow())
     {
      qc8ValueNONAME[countFilesQC] = (qc8CommonHistRes->GetHistIntFlow())->GetBinContent(1);
      qc8ErrorNONAME[countFilesQC] = (qc8CommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(qc8ValueNONAME[countFilesQC]>=0.) // to be improved 
      {
       if(qc8MaxValueNONAME < qc8ValueNONAME[countFilesQC]) qc8MaxValueNONAME = qc8ValueNONAME[countFilesQC]; 
       if(qc8MinValueNONAME > qc8ValueNONAME[countFilesQC]) qc8MinValueNONAME = qc8ValueNONAME[countFilesQC]; 
      } 
      countFilesQC++;
     }
    } // end of if(listQC)
   } // end of if(fileQC) 
  
   // FQD:     
   TString fileNameFQD = gSystem->pwd();
   fileNameFQD+="/";
   fileNameFQD+=presentDir->GetName();
   fileNameFQD+="/";
   fileNameFQD+="outputFQDanalysis";
   (fileNameFQD+=type.Data())+=".root";
   TFile *fileFQD = TFile::Open(fileNameFQD.Data(), "READ");      
   if(fileFQD) 
   {
    TList *listFQD = NULL;
    AliFlowCommonHistResults *fqdCommonHistRes = NULL; 
    fileFQD->GetObject("cobjFQD",listFQD); 
    fileFQD->Close();
    if(listFQD) 
    {
     fqdCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listFQD->FindObject("AliFlowCommonHistResultsFQD")); 
     if(fqdCommonHistRes && fqdCommonHistRes->GetHistIntFlow())
     {
      fqdValueNONAME[countFilesFQD] = (fqdCommonHistRes->GetHistIntFlow())->GetBinContent(1);
      fqdErrorNONAME[countFilesFQD] = (fqdCommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(fqdValueNONAME[countFilesFQD]>=0.) // to be improved 
      {
       if(fqdMaxValueNONAME < fqdValueNONAME[countFilesFQD]) fqdMaxValueNONAME = fqdValueNONAME[countFilesFQD]; 
       if(fqdMinValueNONAME > fqdValueNONAME[countFilesFQD]) fqdMinValueNONAME = fqdValueNONAME[countFilesFQD]; 
      } 
      countFilesFQD++;
     }     
    } // end of if(listFQD)
   } // end of if(fileFQD)   
  
   // LYZ1SUM:     
   TString fileNameLYZ1SUM = gSystem->pwd();
   fileNameLYZ1SUM+="/";
   fileNameLYZ1SUM+=presentDir->GetName();
   fileNameLYZ1SUM+="/";
   fileNameLYZ1SUM+="outputLYZ1SUManalysis";
   (fileNameLYZ1SUM+=type.Data())+=".root";
   TFile *fileLYZ1SUM = TFile::Open(fileNameLYZ1SUM.Data(), "READ");      
   if(fileLYZ1SUM) 
   {
    TList *listLYZ1SUM = NULL;
    AliFlowCommonHistResults *lyz1sumCommonHistRes = NULL; 
    fileLYZ1SUM->GetObject("cobjLYZ1SUM",listLYZ1SUM); 
    fileLYZ1SUM->Close();
    if(listLYZ1SUM) 
    {
     lyz1sumCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listLYZ1SUM->FindObject("AliFlowCommonHistResultsLYZ1SUM")); 
     if(lyz1sumCommonHistRes && lyz1sumCommonHistRes->GetHistIntFlow())
     {
      lyz1sumValueNONAME[countFilesLYZ1SUM] = (lyz1sumCommonHistRes->GetHistIntFlow())->GetBinContent(1);
      lyz1sumErrorNONAME[countFilesLYZ1SUM] = (lyz1sumCommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(lyz1sumValueNONAME[countFilesLYZ1SUM]>=0.) // to be improved 
      {
       if(lyz1sumMaxValueNONAME < lyz1sumValueNONAME[countFilesLYZ1SUM]) lyz1sumMaxValueNONAME = lyz1sumValueNONAME[countFilesLYZ1SUM]; 
       if(lyz1sumMinValueNONAME > lyz1sumValueNONAME[countFilesLYZ1SUM]) lyz1sumMinValueNONAME = lyz1sumValueNONAME[countFilesLYZ1SUM]; 
      } 
      countFilesLYZ1SUM++;
     }
    } // end of if(listLYZ1SUM)
   } // end of if(fileLYZ1SUM)   
  
   // LYZ1PROD:     
   TString fileNameLYZ1PROD = gSystem->pwd();
   fileNameLYZ1PROD+="/";
   fileNameLYZ1PROD+=presentDir->GetName();
   fileNameLYZ1PROD+="/";
   fileNameLYZ1PROD+="outputLYZ1PRODanalysis";
   (fileNameLYZ1PROD+=type.Data())+=".root";
   TFile *fileLYZ1PROD = TFile::Open(fileNameLYZ1PROD.Data(), "READ");      
   if(fileLYZ1PROD) 
   {
    TList *listLYZ1PROD = NULL;
    AliFlowCommonHistResults *lyz1prodCommonHistRes = NULL; 
    fileLYZ1PROD->GetObject("cobjLYZ1PROD",listLYZ1PROD); 
    fileLYZ1PROD->Close();
    if(listLYZ1PROD) 
    {
     lyz1prodCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> (listLYZ1PROD->FindObject("AliFlowCommonHistResultsLYZ1PROD")); 
     if(lyz1prodCommonHistRes && lyz1prodCommonHistRes->GetHistIntFlow())
     {
      lyz1prodValueNONAME[countFilesLYZ1PROD] = (lyz1prodCommonHistRes->GetHistIntFlow())->GetBinContent(1); 
      lyz1prodErrorNONAME[countFilesLYZ1PROD] = (lyz1prodCommonHistRes->GetHistIntFlow())->GetBinError(1);
      if(lyz1prodValueNONAME[countFilesLYZ1PROD]>=0.) // to be improved 
      {
       if(lyz1prodMaxValueNONAME < lyz1prodValueNONAME[countFilesLYZ1PROD]) lyz1prodMaxValueNONAME = lyz1prodValueNONAME[countFilesLYZ1PROD]; 
       if(lyz1prodMinValueNONAME > lyz1prodValueNONAME[countFilesLYZ1PROD]) lyz1prodMinValueNONAME = lyz1prodValueNONAME[countFilesLYZ1PROD]; 
      } 
      countFilesLYZ1PROD++;
     }
    } // end of if(listLYZ1PROD)
   } // end of if(fileLYZ1PROD)    
  } // end of for(Int_t iDir=0;iDir<nDirs;++iDir)  
 } // end of else if(nSubsets > 0)
 
 
 // accessing the large statistics merged output .root files for each method:
 // MCEP:
 TString mergedOutputFileNameMCEP(pwd.Data());
 ((mergedOutputFileNameMCEP+="outputMCEPanalysis")+=type.Data())+=".root";
 TFile *mergedOutputFileMCEP = NULL;
 TList *mergedOutputListMCEP = NULL;
 Double_t mergedValueMCEP = 0.;
 Double_t mergedErrorMCEP = 0.; 
 Double_t mergedAvMMCEP = 0.; // average multiplicity 
 Double_t mergedNMCEP = 0.;  // number of events
 if(gSystem->AccessPathName(mergedOutputFileNameMCEP.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<mergedOutputFileNameMCEP.Data()<<endl;
 } else 
   {     
    mergedOutputFileMCEP = TFile::Open(mergedOutputFileNameMCEP.Data(),"READ");      
    if(mergedOutputFileMCEP) 
    { 
     mergedOutputFileMCEP->GetObject("cobjMCEP",mergedOutputListMCEP); 
     mergedOutputFileMCEP->Close();
     if(mergedOutputListMCEP) 
     {
      AliFlowCommonHistResults *mcepCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> 
                                                    (mergedOutputListMCEP->FindObject("AliFlowCommonHistResultsMCEP")); 
      if(mcepCommonHistRes && mcepCommonHistRes->GetHistIntFlow())
      {
       mergedValueMCEP = (mcepCommonHistRes->GetHistIntFlow())->GetBinContent(1);
       mergedErrorMCEP = (mcepCommonHistRes->GetHistIntFlow())->GetBinError(1);
      } 
      
      AliFlowCommonHist *mcepCommonHist = dynamic_cast<AliFlowCommonHist*> 
                                                    (mergedOutputListMCEP->FindObject("AliFlowCommonHistMCEP")); 
      if(mcepCommonHist && mcepCommonHist->GetHistMultRP())
      {
       mergedAvMMCEP = (mcepCommonHist->GetHistMultRP())->GetMean();
       mergedNMCEP = (mcepCommonHist->GetHistMultRP())->GetEntries();
      }                                                                                                                                              
     }
    } // end of if(mergedOutputFileMCEP)
   } // end of else 
   
 // SP:
 TString mergedOutputFileNameSP(pwd.Data());
 ((mergedOutputFileNameSP+="outputSPanalysis")+=type.Data())+=".root";
 TFile *mergedOutputFileSP = NULL;
 TList *mergedOutputListSP = NULL;
 Double_t mergedValueSP = 0.;
 Double_t mergedErrorSP = 0.; 
 if(gSystem->AccessPathName(mergedOutputFileNameSP.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<mergedOutputFileNameSP.Data()<<endl;
 } else 
   {     
    mergedOutputFileSP = TFile::Open(mergedOutputFileNameSP.Data(),"READ");      
    if(mergedOutputFileSP) 
    { 
     mergedOutputFileSP->GetObject("cobjSP",mergedOutputListSP); 
     mergedOutputFileSP->Close();
     if(mergedOutputListSP) 
     {
      AliFlowCommonHistResults *spCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> 
                                                    (mergedOutputListSP->FindObject("AliFlowCommonHistResultsSP")); 
      if(spCommonHistRes && spCommonHistRes->GetHistIntFlow())
      {
       mergedValueSP = (spCommonHistRes->GetHistIntFlow())->GetBinContent(1);
       mergedErrorSP = (spCommonHistRes->GetHistIntFlow())->GetBinError(1);
      } 
     }
    } // end of if(mergedOutputFileSP)
   } // end of else 
 
 // GFC:
 TString mergedOutputFileNameGFC(pwd.Data());
 ((mergedOutputFileNameGFC+="outputGFCanalysis")+=type.Data())+=".root";
 TFile *mergedOutputFileGFC = NULL;
 TList *mergedOutputListGFC = NULL;
 Double_t mergedValueGFC2 = 0.;
 Double_t mergedErrorGFC2 = 0.; 
 Double_t mergedValueGFC4 = 0.;
 Double_t mergedErrorGFC4 = 0.; 
 Double_t mergedValueGFC6 = 0.;
 Double_t mergedErrorGFC6 = 0.; 
 Double_t mergedValueGFC8 = 0.;
 Double_t mergedErrorGFC8 = 0.; 
 if(gSystem->AccessPathName(mergedOutputFileNameGFC.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<mergedOutputFileNameGFC.Data()<<endl;
 } else 
   {     
    mergedOutputFileGFC = TFile::Open(mergedOutputFileNameGFC.Data(),"READ");      
    if(mergedOutputFileGFC) 
    { 
     mergedOutputFileGFC->GetObject("cobjGFC",mergedOutputListGFC); 
     mergedOutputFileGFC->Close();
     if(mergedOutputListGFC) 
     {
      AliFlowCommonHistResults *gfc2CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> 
                                                    (mergedOutputListGFC->FindObject("AliFlowCommonHistResults2ndOrderGFC")); 
      AliFlowCommonHistResults *gfc4CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> 
                                                    (mergedOutputListGFC->FindObject("AliFlowCommonHistResults4thOrderGFC")); 
      AliFlowCommonHistResults *gfc6CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> 
                                                    (mergedOutputListGFC->FindObject("AliFlowCommonHistResults6thOrderGFC")); 
      AliFlowCommonHistResults *gfc8CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> 
                                                    (mergedOutputListGFC->FindObject("AliFlowCommonHistResults8thOrderGFC")); 
      if(gfc2CommonHistRes && gfc2CommonHistRes->GetHistIntFlow())
      {
       mergedValueGFC2 = (gfc2CommonHistRes->GetHistIntFlow())->GetBinContent(1);
       mergedErrorGFC2 = (gfc2CommonHistRes->GetHistIntFlow())->GetBinError(1);
      } 
      if(gfc4CommonHistRes && gfc4CommonHistRes->GetHistIntFlow())
      {
       mergedValueGFC4 = (gfc4CommonHistRes->GetHistIntFlow())->GetBinContent(1);
       mergedErrorGFC4 = (gfc4CommonHistRes->GetHistIntFlow())->GetBinError(1);
      }
      if(gfc6CommonHistRes && gfc6CommonHistRes->GetHistIntFlow())
      {
       mergedValueGFC6 = (gfc6CommonHistRes->GetHistIntFlow())->GetBinContent(1);
       mergedErrorGFC6 = (gfc6CommonHistRes->GetHistIntFlow())->GetBinError(1);
      } 
      if(gfc8CommonHistRes && gfc8CommonHistRes->GetHistIntFlow())
      {
       mergedValueGFC8 = (gfc8CommonHistRes->GetHistIntFlow())->GetBinContent(1);
       mergedErrorGFC8 = (gfc8CommonHistRes->GetHistIntFlow())->GetBinError(1);
      }  
     }
    } // end of if(mergedOutputFileGFC)
   } // end of else  

 // QC:
 TString mergedOutputFileNameQC(pwd.Data());
 ((mergedOutputFileNameQC+="outputQCanalysis")+=type.Data())+=".root";
 TFile *mergedOutputFileQC = NULL;
 TList *mergedOutputListQC = NULL;
 Double_t mergedValueQC2 = 0.;
 Double_t mergedErrorQC2 = 0.; 
 Double_t mergedValueQC4 = 0.;
 Double_t mergedErrorQC4 = 0.; 
 Double_t mergedValueQC6 = 0.;
 Double_t mergedErrorQC6 = 0.; 
 Double_t mergedValueQC8 = 0.;
 Double_t mergedErrorQC8 = 0.; 
 if(gSystem->AccessPathName(mergedOutputFileNameQC.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<mergedOutputFileNameQC.Data()<<endl;
 } else 
   {     
    mergedOutputFileQC = TFile::Open(mergedOutputFileNameQC.Data(),"READ");      
    if(mergedOutputFileQC) 
    { 
     mergedOutputFileQC->GetObject("cobjQC",mergedOutputListQC); 
     mergedOutputFileQC->Close();
     if(mergedOutputListQC) 
     {
      AliFlowCommonHistResults *qc2CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> 
                                                    (mergedOutputListQC->FindObject("AliFlowCommonHistResults2ndOrderQC")); 
      AliFlowCommonHistResults *qc4CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> 
                                                    (mergedOutputListQC->FindObject("AliFlowCommonHistResults4thOrderQC")); 
      AliFlowCommonHistResults *qc6CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> 
                                                    (mergedOutputListQC->FindObject("AliFlowCommonHistResults6thOrderQC")); 
      AliFlowCommonHistResults *qc8CommonHistRes = dynamic_cast<AliFlowCommonHistResults*> 
                                                    (mergedOutputListQC->FindObject("AliFlowCommonHistResults8thOrderQC")); 
      if(qc2CommonHistRes && qc2CommonHistRes->GetHistIntFlow())
      {
       mergedValueQC2 = (qc2CommonHistRes->GetHistIntFlow())->GetBinContent(1);
       mergedErrorQC2 = (qc2CommonHistRes->GetHistIntFlow())->GetBinError(1);
      } 
      if(qc4CommonHistRes && qc4CommonHistRes->GetHistIntFlow())
      {
       mergedValueQC4 = (qc4CommonHistRes->GetHistIntFlow())->GetBinContent(1);
       mergedErrorQC4 = (qc4CommonHistRes->GetHistIntFlow())->GetBinError(1);
      }
      if(qc6CommonHistRes && qc6CommonHistRes->GetHistIntFlow())
      {
       mergedValueQC6 = (qc6CommonHistRes->GetHistIntFlow())->GetBinContent(1);
       mergedErrorQC6 = (qc6CommonHistRes->GetHistIntFlow())->GetBinError(1);
      } 
      if(qc8CommonHistRes && qc8CommonHistRes->GetHistIntFlow())
      {
       mergedValueQC8 = (qc8CommonHistRes->GetHistIntFlow())->GetBinContent(1);
       mergedErrorQC8 = (qc8CommonHistRes->GetHistIntFlow())->GetBinError(1);
      }  
     }
    } // end of if(mergedOutputFileGFC)
   } // end of else  

 // FQD:
 TString mergedOutputFileNameFQD(pwd.Data());
 ((mergedOutputFileNameFQD+="outputFQDanalysis")+=type.Data())+=".root";
 TFile *mergedOutputFileFQD = NULL;
 TList *mergedOutputListFQD = NULL;
 Double_t mergedValueFQD = 0.;
 Double_t mergedErrorFQD = 0.; 
 if(gSystem->AccessPathName(mergedOutputFileNameFQD.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<mergedOutputFileNameFQD.Data()<<endl;
 } else 
   {     
    mergedOutputFileFQD = TFile::Open(mergedOutputFileNameFQD.Data(),"READ");      
    if(mergedOutputFileFQD) 
    { 
     mergedOutputFileFQD->GetObject("cobjFQD",mergedOutputListFQD); 
     mergedOutputFileFQD->Close();
     if(mergedOutputListFQD) 
     {
      AliFlowCommonHistResults *fqdCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> 
                                                    (mergedOutputListFQD->FindObject("AliFlowCommonHistResultsFQD")); 
      if(fqdCommonHistRes && fqdCommonHistRes->GetHistIntFlow())
      {
       mergedValueFQD = (fqdCommonHistRes->GetHistIntFlow())->GetBinContent(1);
       mergedErrorFQD = (fqdCommonHistRes->GetHistIntFlow())->GetBinError(1);
      } 
     }
    } // end of if(mergedOutputFileFQD)
   } // end of else  
   
 // LYZ1SUM:
 TString mergedOutputFileNameLYZ1SUM(pwd.Data());
 ((mergedOutputFileNameLYZ1SUM+="outputLYZ1SUManalysis")+=type.Data())+=".root";
 TFile *mergedOutputFileLYZ1SUM = NULL;
 TList *mergedOutputListLYZ1SUM = NULL;
 Double_t mergedValueLYZ1SUM = 0.;
 Double_t mergedErrorLYZ1SUM = 0.; 
 if(gSystem->AccessPathName(mergedOutputFileNameLYZ1SUM.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<mergedOutputFileNameLYZ1SUM.Data()<<endl;
 } else 
   {     
    mergedOutputFileLYZ1SUM = TFile::Open(mergedOutputFileNameLYZ1SUM.Data(),"READ");      
    if(mergedOutputFileLYZ1SUM) 
    { 
     mergedOutputFileLYZ1SUM->GetObject("cobjLYZ1SUM",mergedOutputListLYZ1SUM); 
     mergedOutputFileLYZ1SUM->Close();
     if(mergedOutputListLYZ1SUM) 
     {
      AliFlowCommonHistResults *lyz1sumCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> 
                                                    (mergedOutputListLYZ1SUM->FindObject("AliFlowCommonHistResultsLYZ1SUM")); 
      if(lyz1sumCommonHistRes && lyz1sumCommonHistRes->GetHistIntFlow())
      {
       mergedValueLYZ1SUM = (lyz1sumCommonHistRes->GetHistIntFlow())->GetBinContent(1);
       mergedErrorLYZ1SUM = (lyz1sumCommonHistRes->GetHistIntFlow())->GetBinError(1);
      } 
     }
    } // end of if(mergedOutputFileLYZ1SUM)
   } // end of else  
   
 // LYZ1PROD:
 TString mergedOutputFileNameLYZ1PROD(pwd.Data());
 ((mergedOutputFileNameLYZ1PROD+="outputLYZ1PRODanalysis")+=type.Data())+=".root";
 TFile *mergedOutputFileLYZ1PROD = NULL;
 TList *mergedOutputListLYZ1PROD = NULL;
 Double_t mergedValueLYZ1PROD = 0.;
 Double_t mergedErrorLYZ1PROD = 0.; 
 if(gSystem->AccessPathName(mergedOutputFileNameLYZ1PROD.Data(),kFileExists))
 {
  cout<<"WARNING: You do not have a merged output file "<<mergedOutputFileNameLYZ1PROD.Data()<<endl;
 } else 
   {     
    mergedOutputFileLYZ1PROD = TFile::Open(mergedOutputFileNameLYZ1PROD.Data(),"READ");      
    if(mergedOutputFileLYZ1PROD) 
    { 
     mergedOutputFileLYZ1PROD->GetObject("cobjLYZ1PROD",mergedOutputListLYZ1PROD); 
     mergedOutputFileLYZ1PROD->Close();
     if(mergedOutputListLYZ1PROD) 
     {
      AliFlowCommonHistResults *lyz1prodCommonHistRes = dynamic_cast<AliFlowCommonHistResults*> 
                                                    (mergedOutputListLYZ1PROD->FindObject("AliFlowCommonHistResultsLYZ1PROD")); 
      if(lyz1prodCommonHistRes && lyz1prodCommonHistRes->GetHistIntFlow())
      {
       mergedValueLYZ1PROD = (lyz1prodCommonHistRes->GetHistIntFlow())->GetBinContent(1);
       mergedErrorLYZ1PROD = (lyz1prodCommonHistRes->GetHistIntFlow())->GetBinError(1);
      } 
     }
    } // end of if(mergedOutputFileLYZ1PROD)
   } // end of else    
                   
 // removing the title and stat. box from all histograms:
 // gStyle->SetOptTitle(0);
 gStyle->SetOptStat(0);  
 
 // box width:
 const Double_t boxWidth = 0.25;
  
 // the number of different methods:
 const Int_t nMethods = 14;
 
 // the number of small statistics runs:
 const Int_t nPointsMCEP = countFilesMCEP;
 const Int_t nPointsSP = countFilesSP;
 const Int_t nPointsGFC = countFilesGFC;
 const Int_t nPointsQC = countFilesQC;
 const Int_t nPointsFQD = countFilesFQD;
 const Int_t nPointsLYZ1SUM = countFilesLYZ1SUM;   
 const Int_t nPointsLYZ1PROD = countFilesLYZ1PROD; 
 const Int_t nPointsLYZEP = countFilesLYZEP;    
 
 // booking the style histogram for the integrated flow results from all methods for NONAME, RPs and POIs:
 TH1D* intFlowAll = new TH1D("intFlowAll","Integrated Flow",nMethods,0,nMethods);      
 if(!plotLegend) intFlowAll->SetLabelSize(0.044,"X");
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
 (intFlowAll->GetXaxis())->SetBinLabel(binLYZ1SUM,"v_{2}{LYZ,sum}");
 (intFlowAll->GetXaxis())->SetBinLabel(binLYZ1PROD,"v_{2}{LYZ,prod}");
 (intFlowAll->GetXaxis())->SetBinLabel(binLYZEP,"v_{2}{LYZEP}");
 
 //=============================================================================
 //                             NONAME histogram
 //=============================================================================
   
 TH1D *intFlowNONAME = new TH1D(*intFlowAll); 
 
 Double_t mcepNONAME[nPointsMCEP];  
 Double_t spNONAME[nPointsSP]; 
 Double_t gfc2NONAME[nPointsGFC];
 Double_t gfc4NONAME[nPointsGFC];
 Double_t gfc6NONAME[nPointsGFC];
 Double_t gfc8NONAME[nPointsGFC];
 Double_t qc2NONAME[nPointsQC];
 Double_t qc4NONAME[nPointsQC];
 Double_t qc6NONAME[nPointsQC];
 Double_t qc8NONAME[nPointsQC];
 Double_t fqdNONAME[nPointsFQD]; 
 Double_t lyz1sumNONAME[nPointsLYZ1SUM]; 
 Double_t lyz1prodNONAME[nPointsLYZ1PROD]; 
 Double_t lyzepNONAME[nPointsLYZEP];
 
 for(Int_t i=0;i<nPointsMCEP;i++)
 {
  mcepNONAME[i]=binMCEP-0.5;
 } 
 for(Int_t i=0;i<nPointsSP;i++)
 {
  spNONAME[i]=binSP-0.5;
 } 
 for(Int_t i=0;i<nPointsGFC;i++)
 {
  gfc2NONAME[i]=binGFC2-0.5;
  gfc4NONAME[i]=binGFC4-0.5;
  gfc6NONAME[i]=binGFC6-0.5;
  gfc8NONAME[i]=binGFC8-0.5;
 } 
 for(Int_t i=0;i<nPointsQC;i++)
 {
  qc2NONAME[i]=binQC2-0.5;
  qc4NONAME[i]=binQC4-0.5;
  qc6NONAME[i]=binQC6-0.5;
  qc8NONAME[i]=binQC8-0.5;
 }  
 for(Int_t i=0;i<nPointsFQD;i++)
 {
  fqdNONAME[i]=binFQD-0.5;
 } 
 for(Int_t i=0;i<nPointsLYZ1SUM;i++)
 {
  lyz1sumNONAME[i]=binLYZ1SUM-0.5;
 }  
 for(Int_t i=0;i<nPointsLYZ1PROD;i++)
 {
  lyz1prodNONAME[i]=binLYZ1PROD-0.5;
 }  
 for(Int_t i=0;i<nPointsLYZEP;i++)
 {
  lyzepNONAME[i]=binLYZEP-0.5;
 }  
 
 // MCEP:
 TGraphErrors *mcepMeanNONAME = new TGraphErrors(1);    
 mcepMeanNONAME->SetPoint(0,binMCEP-0.5,mergedValueMCEP);
 mcepMeanNONAME->SetPointError(0,0,mergedErrorMCEP);   
 mcepMeanNONAME->SetMarkerStyle(25);
 mcepMeanNONAME->SetMarkerColor(kBlack); 
 mcepMeanNONAME->SetLineColor(kBlack);
 mcepMeanNONAME->SetMarkerSize(1.25); 
 
 TGraph *mcepTGraphNONAME = NULL;
 if(nPointsMCEP>0)
 {
  mcepTGraphNONAME = new TGraph(nPointsMCEP, mcepNONAME, mcepValueNONAME);
  mcepTGraphNONAME->SetMarkerStyle(21);
  mcepTGraphNONAME->SetMarkerColor(kGray+1); 
  mcepTGraphNONAME->SetMarkerSize(0.75); 
 }
 
 TGraph *mcepBoxNONAME = new TGraph(5);
 mcepBoxNONAME->SetPoint(0,(binMCEP-0.5)-boxWidth,mcepMinValueNONAME);
 mcepBoxNONAME->SetPoint(1,(binMCEP-0.5)+boxWidth,mcepMinValueNONAME);
 mcepBoxNONAME->SetPoint(2,(binMCEP-0.5)+boxWidth,mcepMaxValueNONAME);
 mcepBoxNONAME->SetPoint(3,(binMCEP-0.5)-boxWidth,mcepMaxValueNONAME);
 mcepBoxNONAME->SetPoint(4,(binMCEP-0.5)-boxWidth,mcepMinValueNONAME);    
 mcepBoxNONAME->SetFillStyle(1001);
 mcepBoxNONAME->SetFillColor(kGray);
 
 // SP:
 TGraphErrors *spMeanNONAME = new TGraphErrors(1);    
 spMeanNONAME->SetPoint(0,binSP-0.5,mergedValueSP);  
 spMeanNONAME->SetPointError(0,0,mergedErrorSP);  
 spMeanNONAME->SetMarkerStyle(25);
 spMeanNONAME->SetMarkerColor(kBlack); 
 spMeanNONAME->SetLineColor(kBlack);
 spMeanNONAME->SetMarkerSize(1.25); 
 
 TGraph *spTGraphNONAME = NULL;
 if(nPointsSP>0)
 {
  spTGraphNONAME = new TGraph(nPointsSP, spNONAME, spValueNONAME);
  spTGraphNONAME->SetMarkerStyle(21);
  spTGraphNONAME->SetMarkerColor(kViolet-8); 
  spTGraphNONAME->SetMarkerSize(0.75); 
 }
 
 TGraph *spBoxNONAME = new TGraph(5);
 spBoxNONAME->SetPoint(0,(binSP-0.5)-boxWidth,spMinValueNONAME);
 spBoxNONAME->SetPoint(1,(binSP-0.5)+boxWidth,spMinValueNONAME);
 spBoxNONAME->SetPoint(2,(binSP-0.5)+boxWidth,spMaxValueNONAME);
 spBoxNONAME->SetPoint(3,(binSP-0.5)-boxWidth,spMaxValueNONAME);
 spBoxNONAME->SetPoint(4,(binSP-0.5)-boxWidth,spMinValueNONAME);    
 spBoxNONAME->SetFillStyle(1001);
 spBoxNONAME->SetFillColor(kViolet-9);
 
 // GFC{2}:
 TGraphErrors *gfc2MeanNONAME = new TGraphErrors(1);    
 gfc2MeanNONAME->SetPoint(0,binGFC2-0.5,mergedValueGFC2);  
 gfc2MeanNONAME->SetPointError(0,0,mergedErrorGFC2);  
 gfc2MeanNONAME->SetMarkerStyle(25);
 gfc2MeanNONAME->SetMarkerColor(kBlack); 
 gfc2MeanNONAME->SetLineColor(kBlack);
 gfc2MeanNONAME->SetMarkerSize(1.25); 
 
 
 TGraph *gfc2TGraphNONAME = NULL;
 if(nPointsGFC>0)
 {
  gfc2TGraphNONAME = new TGraph(nPointsGFC, gfc2NONAME, gfc2ValueNONAME);
  gfc2TGraphNONAME->SetMarkerStyle(21);
  gfc2TGraphNONAME->SetMarkerColor(kBlue-9); 
  gfc2TGraphNONAME->SetMarkerSize(0.75); 
 }
 
 TGraph *gfc2BoxNONAME = new TGraph(5);
 gfc2BoxNONAME->SetPoint(0,(binGFC2-0.5)-boxWidth,gfc2MinValueNONAME);
 gfc2BoxNONAME->SetPoint(1,(binGFC2-0.5)+boxWidth,gfc2MinValueNONAME);
 gfc2BoxNONAME->SetPoint(2,(binGFC2-0.5)+boxWidth,gfc2MaxValueNONAME);
 gfc2BoxNONAME->SetPoint(3,(binGFC2-0.5)-boxWidth,gfc2MaxValueNONAME);
 gfc2BoxNONAME->SetPoint(4,(binGFC2-0.5)-boxWidth,gfc2MinValueNONAME);    
 gfc2BoxNONAME->SetFillStyle(1001);
 gfc2BoxNONAME->SetFillColor(kBlue-10);
 
 // GFC{4}:
 TGraphErrors *gfc4MeanNONAME = new TGraphErrors(1);    
 gfc4MeanNONAME->SetPoint(0,binGFC4-0.5,mergedValueGFC4);  
 gfc4MeanNONAME->SetPointError(0,0,mergedErrorGFC4);  
 gfc4MeanNONAME->SetMarkerStyle(25);
 gfc4MeanNONAME->SetMarkerColor(kBlack); 
 gfc4MeanNONAME->SetLineColor(kBlack);
 gfc4MeanNONAME->SetMarkerSize(1.25); 
 
 TGraph *gfc4TGraphNONAME = NULL;
 if(nPointsGFC>0)
 { 
  gfc4TGraphNONAME = new TGraph(nPointsGFC, gfc4NONAME, gfc4ValueNONAME);
  gfc4TGraphNONAME->SetMarkerStyle(21);
  gfc4TGraphNONAME->SetMarkerColor(kBlue-9); 
  gfc4TGraphNONAME->SetMarkerSize(0.75); 
 }
 
 TGraph *gfc4BoxNONAME = new TGraph(5);
 gfc4BoxNONAME->SetPoint(0,(binGFC4-0.5)-boxWidth,gfc4MinValueNONAME);
 gfc4BoxNONAME->SetPoint(1,(binGFC4-0.5)+boxWidth,gfc4MinValueNONAME);
 gfc4BoxNONAME->SetPoint(2,(binGFC4-0.5)+boxWidth,gfc4MaxValueNONAME);
 gfc4BoxNONAME->SetPoint(3,(binGFC4-0.5)-boxWidth,gfc4MaxValueNONAME);
 gfc4BoxNONAME->SetPoint(4,(binGFC4-0.5)-boxWidth,gfc4MinValueNONAME);    
 gfc4BoxNONAME->SetFillStyle(1001);
 gfc4BoxNONAME->SetFillColor(kBlue-10);
 
 // GFC{6}:
 TGraphErrors *gfc6MeanNONAME = new TGraphErrors(1);    
 gfc6MeanNONAME->SetPoint(0,binGFC6-0.5,mergedValueGFC6);  
 gfc6MeanNONAME->SetPointError(0,0,mergedErrorGFC6);  
 gfc6MeanNONAME->SetMarkerStyle(25);
 gfc6MeanNONAME->SetMarkerColor(kBlack); 
 gfc6MeanNONAME->SetLineColor(kBlack);
 gfc6MeanNONAME->SetMarkerSize(1.25); 
 
 TGraph *gfc6TGraphNONAME = NULL;
 if(nPointsGFC>0)
 { 
  gfc6TGraphNONAME = new TGraph(nPointsGFC, gfc6NONAME, gfc6ValueNONAME);
  gfc6TGraphNONAME->SetMarkerStyle(21);
  gfc6TGraphNONAME->SetMarkerColor(kBlue-9); 
  gfc6TGraphNONAME->SetMarkerSize(0.75); 
 }
 
 TGraph *gfc6BoxNONAME = new TGraph(5);
 gfc6BoxNONAME->SetPoint(0,(binGFC6-0.5)-boxWidth,gfc6MinValueNONAME);
 gfc6BoxNONAME->SetPoint(1,(binGFC6-0.5)+boxWidth,gfc6MinValueNONAME);
 gfc6BoxNONAME->SetPoint(2,(binGFC6-0.5)+boxWidth,gfc6MaxValueNONAME);
 gfc6BoxNONAME->SetPoint(3,(binGFC6-0.5)-boxWidth,gfc6MaxValueNONAME);
 gfc6BoxNONAME->SetPoint(4,(binGFC6-0.5)-boxWidth,gfc6MinValueNONAME);    
 gfc6BoxNONAME->SetFillStyle(1001);
 gfc6BoxNONAME->SetFillColor(kBlue-10);
 
 // GFC{8}:
 TGraphErrors *gfc8MeanNONAME = new TGraphErrors(1);    
 gfc8MeanNONAME->SetPoint(0,binGFC8-0.5,mergedValueGFC8);  
 gfc8MeanNONAME->SetPointError(0,0,mergedErrorGFC8);  
 gfc8MeanNONAME->SetMarkerStyle(25);
 gfc8MeanNONAME->SetMarkerColor(kBlack); 
 gfc8MeanNONAME->SetLineColor(kBlack);
 gfc8MeanNONAME->SetMarkerSize(1.25); 
 
 TGraph *gfc8TGraphNONAME = NULL;
 if(nPointsGFC>0)
 { 
  gfc8TGraphNONAME = new TGraph(nPointsGFC, gfc8NONAME, gfc8ValueNONAME);
  gfc8TGraphNONAME->SetMarkerStyle(21);
  gfc8TGraphNONAME->SetMarkerColor(kBlue-9); 
  gfc8TGraphNONAME->SetMarkerSize(0.75); 
 }
 
 TGraph *gfc8BoxNONAME = new TGraph(5);
 gfc8BoxNONAME->SetPoint(0,(binGFC8-0.5)-boxWidth,gfc8MinValueNONAME);
 gfc8BoxNONAME->SetPoint(1,(binGFC8-0.5)+boxWidth,gfc8MinValueNONAME);
 gfc8BoxNONAME->SetPoint(2,(binGFC8-0.5)+boxWidth,gfc8MaxValueNONAME);
 gfc8BoxNONAME->SetPoint(3,(binGFC8-0.5)-boxWidth,gfc8MaxValueNONAME);
 gfc8BoxNONAME->SetPoint(4,(binGFC8-0.5)-boxWidth,gfc8MinValueNONAME);    
 gfc8BoxNONAME->SetFillStyle(1001);
 gfc8BoxNONAME->SetFillColor(kBlue-10);
 
 // QC{2}:
 TGraphErrors *qc2MeanNONAME = new TGraphErrors(1);    
 qc2MeanNONAME->SetPoint(0,binQC2-0.5,mergedValueQC2);  
 qc2MeanNONAME->SetPointError(0,0,mergedErrorQC2);  
 qc2MeanNONAME->SetMarkerStyle(25);
 qc2MeanNONAME->SetMarkerColor(kBlack); 
 qc2MeanNONAME->SetLineColor(kBlack);
 qc2MeanNONAME->SetMarkerSize(1.25); 
 
 TGraph *qc2TGraphNONAME = NULL;
 if(nPointsQC>0)
 {  
  qc2TGraphNONAME = new TGraph(nPointsQC, qc2NONAME, qc2ValueNONAME);
  qc2TGraphNONAME->SetMarkerStyle(21);
  qc2TGraphNONAME->SetMarkerColor(kRed-7); 
  qc2TGraphNONAME->SetMarkerSize(0.75); 
 }
 
 TGraph *qc2BoxNONAME = new TGraph(5);
 qc2BoxNONAME->SetPoint(0,(binQC2-0.5)-boxWidth,qc2MinValueNONAME);
 qc2BoxNONAME->SetPoint(1,(binQC2-0.5)+boxWidth,qc2MinValueNONAME);
 qc2BoxNONAME->SetPoint(2,(binQC2-0.5)+boxWidth,qc2MaxValueNONAME);
 qc2BoxNONAME->SetPoint(3,(binQC2-0.5)-boxWidth,qc2MaxValueNONAME);
 qc2BoxNONAME->SetPoint(4,(binQC2-0.5)-boxWidth,qc2MinValueNONAME);    
 qc2BoxNONAME->SetFillStyle(1001);
 qc2BoxNONAME->SetFillColor(kRed-10);
 
 // QC{4}:
 TGraphErrors *qc4MeanNONAME = new TGraphErrors(1);    
 qc4MeanNONAME->SetPoint(0,binQC4-0.5,mergedValueQC4);  
 qc4MeanNONAME->SetPointError(0,0,mergedErrorQC4);  
 qc4MeanNONAME->SetMarkerStyle(25);
 qc4MeanNONAME->SetMarkerColor(kBlack); 
 qc4MeanNONAME->SetLineColor(kBlack);
 qc4MeanNONAME->SetMarkerSize(1.25); 
 
 TGraph *qc4TGraphNONAME = NULL;
 if(nPointsQC>0)
 {   
  qc4TGraphNONAME = new TGraph(nPointsQC, qc4NONAME, qc4ValueNONAME);
  qc4TGraphNONAME->SetMarkerStyle(21);
  qc4TGraphNONAME->SetMarkerColor(kRed-7); 
  qc4TGraphNONAME->SetMarkerSize(0.75); 
 }
  
 TGraph *qc4BoxNONAME = new TGraph(5);
 qc4BoxNONAME->SetPoint(0,(binQC4-0.5)-boxWidth,qc4MinValueNONAME);
 qc4BoxNONAME->SetPoint(1,(binQC4-0.5)+boxWidth,qc4MinValueNONAME);
 qc4BoxNONAME->SetPoint(2,(binQC4-0.5)+boxWidth,qc4MaxValueNONAME);
 qc4BoxNONAME->SetPoint(3,(binQC4-0.5)-boxWidth,qc4MaxValueNONAME);
 qc4BoxNONAME->SetPoint(4,(binQC4-0.5)-boxWidth,qc4MinValueNONAME);    
 qc4BoxNONAME->SetFillStyle(1001);
 qc4BoxNONAME->SetFillColor(kRed-10);
 
 // QC{6}:
 TGraphErrors *qc6MeanNONAME = new TGraphErrors(1);    
 qc6MeanNONAME->SetPoint(0,binQC6-0.5,mergedValueQC6);  
 qc6MeanNONAME->SetPointError(0,0,mergedErrorQC6);  
 qc6MeanNONAME->SetMarkerStyle(25);
 qc6MeanNONAME->SetMarkerColor(kBlack); 
 qc6MeanNONAME->SetLineColor(kBlack);
 qc6MeanNONAME->SetMarkerSize(1.25); 
 
 TGraph *qc6TGraphNONAME = NULL;
 if(nPointsQC>0)
 {   
  qc6TGraphNONAME = new TGraph(nPointsQC, qc6NONAME, qc6ValueNONAME);
  qc6TGraphNONAME->SetMarkerStyle(21);
  qc6TGraphNONAME->SetMarkerColor(kRed-7);
  qc6TGraphNONAME->SetMarkerSize(0.75); 
 }
 
 TGraph *qc6BoxNONAME = new TGraph(5);
 qc6BoxNONAME->SetPoint(0,(binQC6-0.5)-boxWidth,qc6MinValueNONAME);
 qc6BoxNONAME->SetPoint(1,(binQC6-0.5)+boxWidth,qc6MinValueNONAME);
 qc6BoxNONAME->SetPoint(2,(binQC6-0.5)+boxWidth,qc6MaxValueNONAME);
 qc6BoxNONAME->SetPoint(3,(binQC6-0.5)-boxWidth,qc6MaxValueNONAME);
 qc6BoxNONAME->SetPoint(4,(binQC6-0.5)-boxWidth,qc6MinValueNONAME);    
 qc6BoxNONAME->SetFillStyle(1001);
 qc6BoxNONAME->SetFillColor(kRed-10);
 
 // QC{8}:
 TGraphErrors *qc8MeanNONAME = new TGraphErrors(1);    
 qc8MeanNONAME->SetPoint(0,binQC8-0.5,mergedValueQC8);  
 qc8MeanNONAME->SetPointError(0,0,mergedErrorQC8);  
 qc8MeanNONAME->SetMarkerStyle(25);
 qc8MeanNONAME->SetMarkerColor(kBlack);
 qc8MeanNONAME->SetLineColor(kBlack); 
 qc8MeanNONAME->SetMarkerSize(1.25); 
 
 TGraph *qc8TGraphNONAME = NULL;
 if(nPointsQC>0)
 {   
  qc8TGraphNONAME = new TGraph(nPointsQC, qc8NONAME, qc8ValueNONAME);
  qc8TGraphNONAME->SetMarkerStyle(21);
  qc8TGraphNONAME->SetMarkerColor(kRed-7); 
  qc8TGraphNONAME->SetMarkerSize(0.75); 
 }
  
 TGraph *qc8BoxNONAME = new TGraph(5);
 qc8BoxNONAME->SetPoint(0,(binQC8-0.5)-boxWidth,qc8MinValueNONAME);
 qc8BoxNONAME->SetPoint(1,(binQC8-0.5)+boxWidth,qc8MinValueNONAME);
 qc8BoxNONAME->SetPoint(2,(binQC8-0.5)+boxWidth,qc8MaxValueNONAME);
 qc8BoxNONAME->SetPoint(3,(binQC8-0.5)-boxWidth,qc8MaxValueNONAME);
 qc8BoxNONAME->SetPoint(4,(binQC8-0.5)-boxWidth,qc8MinValueNONAME);    
 qc8BoxNONAME->SetFillStyle(1001);
 qc8BoxNONAME->SetFillColor(kRed-10);
 
 // FQD:
 TGraphErrors *fqdMeanNONAME = new TGraphErrors(1);    
 fqdMeanNONAME->SetPoint(0,binFQD-0.5,mergedValueFQD);  
 fqdMeanNONAME->SetPointError(0,0,mergedErrorFQD);  
 fqdMeanNONAME->SetMarkerStyle(25);
 fqdMeanNONAME->SetMarkerColor(kBlack); 
 fqdMeanNONAME->SetLineColor(kBlack);
 fqdMeanNONAME->SetMarkerSize(1.25); 
 
 TGraph *fqdTGraphNONAME = NULL;
 if(nPointsFQD>0)
 {  
  fqdTGraphNONAME = new TGraph(nPointsFQD, fqdNONAME, fqdValueNONAME);
  fqdTGraphNONAME->SetMarkerStyle(21);
  fqdTGraphNONAME->SetMarkerColor(kOrange-8); 
  fqdTGraphNONAME->SetMarkerSize(0.75); 
 }
 
 TGraph *fqdBoxNONAME = new TGraph(5);
 fqdBoxNONAME->SetPoint(0,(binFQD-0.5)-boxWidth,fqdMinValueNONAME);
 fqdBoxNONAME->SetPoint(1,(binFQD-0.5)+boxWidth,fqdMinValueNONAME);
 fqdBoxNONAME->SetPoint(2,(binFQD-0.5)+boxWidth,fqdMaxValueNONAME);
 fqdBoxNONAME->SetPoint(3,(binFQD-0.5)-boxWidth,fqdMaxValueNONAME);
 fqdBoxNONAME->SetPoint(4,(binFQD-0.5)-boxWidth,fqdMinValueNONAME);    
 fqdBoxNONAME->SetFillStyle(1001);
 fqdBoxNONAME->SetFillColor(kOrange-9);
 
 // LYZ1SUM:
 TGraphErrors *lyz1sumMeanNONAME = new TGraphErrors(1);    
 lyz1sumMeanNONAME->SetPoint(0,binLYZ1SUM-0.5,mergedValueLYZ1SUM);  
 lyz1sumMeanNONAME->SetPointError(0,0,mergedErrorLYZ1SUM);  
 lyz1sumMeanNONAME->SetMarkerStyle(25);
 lyz1sumMeanNONAME->SetMarkerColor(kBlack); 
 lyz1sumMeanNONAME->SetLineColor(kBlack);
 lyz1sumMeanNONAME->SetMarkerSize(1.25); 
 
 TGraph *lyz1sumTGraphNONAME = NULL;
 if(nPointsLYZ1SUM>0)
 { 
  lyz1sumTGraphNONAME = new TGraph(nPointsLYZ1SUM, lyz1sumNONAME, lyz1sumValueNONAME);
  lyz1sumTGraphNONAME->SetMarkerStyle(21);
  lyz1sumTGraphNONAME->SetMarkerColor(kYellow-5); 
  lyz1sumTGraphNONAME->SetMarkerSize(0.75); 
 }
 
 TGraph *lyz1sumBoxNONAME = new TGraph(5);
 lyz1sumBoxNONAME->SetPoint(0,(binLYZ1SUM-0.5)-boxWidth,lyz1sumMinValueNONAME);
 lyz1sumBoxNONAME->SetPoint(1,(binLYZ1SUM-0.5)+boxWidth,lyz1sumMinValueNONAME);
 lyz1sumBoxNONAME->SetPoint(2,(binLYZ1SUM-0.5)+boxWidth,lyz1sumMaxValueNONAME);
 lyz1sumBoxNONAME->SetPoint(3,(binLYZ1SUM-0.5)-boxWidth,lyz1sumMaxValueNONAME);
 lyz1sumBoxNONAME->SetPoint(4,(binLYZ1SUM-0.5)-boxWidth,lyz1sumMinValueNONAME);    
 lyz1sumBoxNONAME->SetFillStyle(1001);
 lyz1sumBoxNONAME->SetFillColor(kYellow-8);   
 
 // LYZ1PROD:
 TGraphErrors *lyz1prodMeanNONAME = new TGraphErrors(1);    
 lyz1prodMeanNONAME->SetPoint(0,binLYZ1PROD-0.5,mergedValueLYZ1PROD);  
 lyz1prodMeanNONAME->SetPointError(0,0,mergedErrorLYZ1PROD);  
 lyz1prodMeanNONAME->SetMarkerStyle(25);
 lyz1prodMeanNONAME->SetMarkerColor(kBlack); 
 lyz1prodMeanNONAME->SetLineColor(kBlack);
 lyz1prodMeanNONAME->SetMarkerSize(1.25); 
 
 TGraph *lyz1prodTGraphNONAME = NULL;
 if(nPointsLYZ1PROD>0)
 {  
  lyz1prodTGraphNONAME = new TGraph(nPointsLYZ1PROD, lyz1prodNONAME, lyz1prodValueNONAME);
  lyz1prodTGraphNONAME->SetMarkerStyle(21);
  lyz1prodTGraphNONAME->SetMarkerColor(kYellow-2); 
  lyz1prodTGraphNONAME->SetMarkerSize(0.75); 
 }
 
 TGraph *lyz1prodBoxNONAME = new TGraph(5);
 lyz1prodBoxNONAME->SetPoint(0,(binLYZ1PROD-0.5)-boxWidth,lyz1prodMinValueNONAME);
 lyz1prodBoxNONAME->SetPoint(1,(binLYZ1PROD-0.5)+boxWidth,lyz1prodMinValueNONAME);
 lyz1prodBoxNONAME->SetPoint(2,(binLYZ1PROD-0.5)+boxWidth,lyz1prodMaxValueNONAME);
 lyz1prodBoxNONAME->SetPoint(3,(binLYZ1PROD-0.5)-boxWidth,lyz1prodMaxValueNONAME);
 lyz1prodBoxNONAME->SetPoint(4,(binLYZ1PROD-0.5)-boxWidth,lyz1prodMinValueNONAME);    
 lyz1prodBoxNONAME->SetFillStyle(1001);
 lyz1prodBoxNONAME->SetFillColor(kYellow-6);
 
 TCanvas* intFlowCanvasNONAME = new TCanvas("Integrated Flow NONAME","Integrated Flow NONAME",1000,600);
 
 if(plotLegend)
 {
  intFlowCanvasNONAME->Divide(2,1);
  
  // 1st pad is for plot:
  (intFlowCanvasNONAME->cd(1))->SetPad(0.0,0.0,0.75,1.0);
 }
 
 if(intFlowNONAME) 
 {
  TString intFlowNameNONAME("Superimposing ");
  if(nPointsMCEP != 0) // to be improved
  {
   intFlowNameNONAME+=nPointsMCEP;
  } else if (nPointsSP != 0)  
    {
     intFlowNameNONAME+=nPointsSP;
    } else if (nPointsGFC != 0)  
      {
       intFlowNameNONAME+=nPointsGFC;
      } else if (nPointsQC != 0)  
        {
         intFlowNameNONAME+=nPointsQC;
        } else if (nPointsFQD != 0)  
          {
           intFlowNameNONAME+=nPointsFQD;
          } else if (nPointsLYZ1SUM != 0)  
            {
             intFlowNameNONAME+=nPointsLYZ1SUM;
            } else if (nPointsLYZ1PROD != 0)  
              {
               intFlowNameNONAME+=nPointsLYZ1PROD;
              } else if (nPointsLYZEP != 0)  
                {
                 intFlowNameNONAME+=nPointsLYZEP;
                }
       
  intFlowNameNONAME+=" independent runs";
  intFlowNONAME->SetTitle(intFlowNameNONAME.Data());
  (intFlowNONAME->GetYaxis())->SetRangeUser(0,0.20);
  intFlowNONAME->Draw();
 }
 
 // MCEP: 
 if(mcepBoxNONAME && mcepMinValueNONAME < 1000.) mcepBoxNONAME->Draw("LFSAME");
 if(mcepTGraphNONAME) mcepTGraphNONAME->Draw("PSAME");
 if(mcepMeanNONAME) mcepMeanNONAME->Draw("PSAME");
 
 // SP: 
 if(spBoxNONAME && spMinValueNONAME < 1000.) spBoxNONAME->Draw("LFSAME");
 if(spTGraphNONAME) spTGraphNONAME->Draw("PSAME");
 if(spMeanNONAME) spMeanNONAME->Draw("PSAME");
 
 // GFC{2}: 
 if(gfc2BoxNONAME && gfc2MinValueNONAME < 1000.) gfc2BoxNONAME->Draw("LFSAME");
 if(gfc2TGraphNONAME) gfc2TGraphNONAME->Draw("PSAME");
 if(gfc2MeanNONAME) gfc2MeanNONAME->Draw("PSAME");
 
 // GFC{4}: 
 if(gfc4BoxNONAME && gfc4MinValueNONAME < 1000.) gfc4BoxNONAME->Draw("LFSAME");
 if(gfc4TGraphNONAME) gfc4TGraphNONAME->Draw("PSAME");
 if(gfc4MeanNONAME) gfc4MeanNONAME->Draw("PSAME");
 
 // GFC{6}: 
 if(gfc6BoxNONAME && gfc6MinValueNONAME < 1000.) gfc6BoxNONAME->Draw("LFSAME");
 if(gfc6TGraphNONAME) gfc6TGraphNONAME->Draw("PSAME");
 if(gfc6MeanNONAME) gfc6MeanNONAME->Draw("PSAME");
 
 // GFC{8}: 
 if(gfc8BoxNONAME && gfc8MinValueNONAME < 1000.) gfc8BoxNONAME->Draw("LFSAME");
 if(gfc8TGraphNONAME) gfc8TGraphNONAME->Draw("PSAME");
 if(gfc8MeanNONAME) gfc8MeanNONAME->Draw("PSAME");
 
 // QC{2}: 
 if(qc2BoxNONAME && qc2MinValueNONAME < 1000.) qc2BoxNONAME->Draw("LFSAME");
 if(qc2TGraphNONAME) qc2TGraphNONAME->Draw("PSAME");
 if(qc2MeanNONAME) qc2MeanNONAME->Draw("PSAME");
 
 // QC{4}: 
 if(qc4BoxNONAME && qc4MinValueNONAME < 1000.) qc4BoxNONAME->Draw("LFSAME");
 if(qc4TGraphNONAME) qc4TGraphNONAME->Draw("PSAME");
 if(qc4MeanNONAME) qc4MeanNONAME->Draw("PSAME");
 
 // QC{6}: 
 if(qc6BoxNONAME && qc6MinValueNONAME < 1000.) qc6BoxNONAME->Draw("LFSAME");
 if(qc6TGraphNONAME) qc6TGraphNONAME->Draw("PSAME");
 if(qc6MeanNONAME) qc6MeanNONAME->Draw("PSAME");
 
 // QC{8}: 
 if(qc8BoxNONAME && qc8MinValueNONAME < 1000.) qc8BoxNONAME->Draw("LFSAME");
 if(qc8TGraphNONAME) qc8TGraphNONAME->Draw("PSAME");
 if(qc8MeanNONAME) qc8MeanNONAME->Draw("PSAME");
 
 // FQD: 
 if(fqdBoxNONAME && fqdMinValueNONAME < 1000.) fqdBoxNONAME->Draw("LFSAME");
 if(fqdTGraphNONAME) fqdTGraphNONAME->Draw("PSAME");
 if(fqdMeanNONAME) fqdMeanNONAME->Draw("PSAME");

 // LYZ1SUM: 
 if(lyz1sumBoxNONAME && lyz1sumMinValueNONAME < 1000.) lyz1sumBoxNONAME->Draw("LFSAME");
 if(lyz1sumTGraphNONAME) lyz1sumTGraphNONAME->Draw("PSAME");
 if(lyz1sumMeanNONAME) lyz1sumMeanNONAME->Draw("PSAME"); 
 
 // LYZ1PROD: 
 if(lyz1prodBoxNONAME && lyz1prodMinValueNONAME < 1000.) lyz1prodBoxNONAME->Draw("LFSAME");
 if(lyz1prodTGraphNONAME) lyz1prodTGraphNONAME->Draw("PSAME");
 if(lyz1prodMeanNONAME) lyz1prodMeanNONAME->Draw("PSAME"); 
 
 if(plotLegend)
 {
  // 2nd pad is for legend:   
  (intFlowCanvasNONAME->cd(2))->SetPad(0.75,0.0,1.0,1.0);
 }

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
 Int_t lyz1sumCountRealNONAME = 0;
 Int_t lyz1prodCountRealNONAME = 0;
 Int_t lyzepCountRealNONAME = 0;
 
 Int_t N = 0;
 if(nSubsets > 0)
 {
  N = nSubsets;
 } else
   {
    N = nEstimates;
   }
   
 for(Int_t i=0;i<N;i++)
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
  if(lyz1sumValueNONAME[i]>0.) lyz1sumCountRealNONAME++; 
  if(lyz1prodValueNONAME[i]>0.) lyz1prodCountRealNONAME++;
  if(lyzepValueNONAME[i]>0.) lyzepCountRealNONAME++;
 }

 TPaveText *textDefaultNONAME = new TPaveText(0.05,0.71,0.95,0.95,"NDC");
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
 
 // average multiplicity and number of events: (to be improved)
 TPaveText *textAvMNevtsNONAME = new TPaveText(0.05,0.05,0.95,0.14,"NDC");
 textAvMNevtsNONAME->SetTextFont(72);
 textAvMNevtsNONAME->SetTextSize(0.06);
 
 TString *entryNevtsNONAME = new TString("MC .... nEvts = ");
 TString *entryAvMNONAME   = new TString("MC .... AvM = ");
 
 (*entryNevtsNONAME)+=(Long_t)mergedNMCEP;
 (*entryAvMNONAME)+=(Long_t)mergedAvMMCEP;
 
 textAvMNevtsNONAME->AddText(entryNevtsNONAME->Data());
 textAvMNevtsNONAME->AddText(entryAvMNONAME->Data());
 
 // results:
 TPaveText *textResultsNONAME = new TPaveText(0.05,0.17,0.95,0.68,"NDC");
 textResultsNONAME->SetTextFont(72);
 textResultsNONAME->SetTextSize(0.06);
 
 // entries:
 TString *entryIntFlowMCNONAME       = new TString("MC ................ ");
 TString *entryIntFlowSPNONAME       = new TString("SP ................ ");  
 TString *entryIntFlowGFC2NONAME     = new TString("GFC{2} ........ ");
 TString *entryIntFlowGFC4NONAME     = new TString("GFC{4} ........ ");
 TString *entryIntFlowGFC6NONAME     = new TString("GFC{6} ........ ");
 TString *entryIntFlowGFC8NONAME     = new TString("GFC{8} ........ "); 
 TString *entryIntFlowQC2NONAME      = new TString("QC{2} .......... ");
 TString *entryIntFlowQC4NONAME      = new TString("QC{4} .......... ");
 TString *entryIntFlowQC6NONAME      = new TString("QC{6} .......... ");
 TString *entryIntFlowQC8NONAME      = new TString("QC{8} .......... ");
 TString *entryIntFlowFQDNONAME      = new TString("FQD ............. ");
 TString *entryIntFlowLYZ1SUMNONAME  = new TString("LYZ{sum} ......... ");
 TString *entryIntFlowLYZ1PRODNONAME = new TString("LYZ{prod} ........ ");
 TString *entryIntFlowLYZEPNONAME    = new TString("LYZEP ........ ");

 if(entryIntFlowMCNONAME)
 { 
  (*entryIntFlowMCNONAME)+=(Long_t)mcepCountRealNONAME;
  entryIntFlowMCNONAME->Append(" out of ");
  (*entryIntFlowMCNONAME)+=(Long_t)nPointsMCEP;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowMCNONAME->Data());
 }
 
 if(entryIntFlowSPNONAME)
 { 
  (*entryIntFlowSPNONAME)+=(Long_t)spCountRealNONAME;
  entryIntFlowSPNONAME->Append(" out of ");
  (*entryIntFlowSPNONAME)+=(Long_t)nPointsSP;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowSPNONAME->Data());
 }
 
 if(entryIntFlowGFC2NONAME)
 { 
  (*entryIntFlowGFC2NONAME)+=(Long_t)gfc2CountRealNONAME;
  entryIntFlowGFC2NONAME->Append(" out of ");
  (*entryIntFlowGFC2NONAME)+=(Long_t)nPointsGFC;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowGFC2NONAME->Data());
 }
 
 if(entryIntFlowGFC4NONAME)
 { 
  (*entryIntFlowGFC4NONAME)+=(Long_t)gfc4CountRealNONAME;
  entryIntFlowGFC4NONAME->Append(" out of ");
  (*entryIntFlowGFC4NONAME)+=(Long_t)nPointsGFC;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowGFC4NONAME->Data());
 }

 if(entryIntFlowGFC6NONAME)
 { 
  (*entryIntFlowGFC6NONAME)+=(Long_t)gfc6CountRealNONAME;
  entryIntFlowGFC6NONAME->Append(" out of ");
  (*entryIntFlowGFC6NONAME)+=(Long_t)nPointsGFC;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowGFC6NONAME->Data());
 }

 if(entryIntFlowGFC8NONAME)
 { 
  (*entryIntFlowGFC8NONAME)+=(Long_t)gfc8CountRealNONAME;
  entryIntFlowGFC8NONAME->Append(" out of ");
  (*entryIntFlowGFC8NONAME)+=(Long_t)nPointsGFC;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowGFC8NONAME->Data());
 }
 
 if(entryIntFlowQC2NONAME)
 { 
  (*entryIntFlowQC2NONAME)+=(Long_t)qc2CountRealNONAME;
  entryIntFlowQC2NONAME->Append(" out of ");
  (*entryIntFlowQC2NONAME)+=(Long_t)nPointsQC;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowQC2NONAME->Data());
 }
 
 if(entryIntFlowQC4NONAME)
 { 
  (*entryIntFlowQC4NONAME)+=(Long_t)qc4CountRealNONAME;
  entryIntFlowQC4NONAME->Append(" out of ");
  (*entryIntFlowQC4NONAME)+=(Long_t)nPointsQC;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowQC4NONAME->Data());
 }

 if(entryIntFlowQC2NONAME)
 { 
  (*entryIntFlowQC6NONAME)+=(Long_t)qc6CountRealNONAME;
  entryIntFlowQC6NONAME->Append(" out of ");
  (*entryIntFlowQC6NONAME)+=(Long_t)nPointsQC;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowQC6NONAME->Data());
 }

 if(entryIntFlowQC8NONAME)
 { 
  (*entryIntFlowQC8NONAME)+=(Long_t)qc8CountRealNONAME;
  entryIntFlowQC8NONAME->Append(" out of ");
  (*entryIntFlowQC8NONAME)+=(Long_t)nPointsQC;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowQC8NONAME->Data());
 }
 
 if(entryIntFlowFQDNONAME)
 { 
  (*entryIntFlowFQDNONAME)+=(Long_t)fqdCountRealNONAME;
  entryIntFlowFQDNONAME->Append(" out of ");
  (*entryIntFlowFQDNONAME)+=(Long_t)nPointsFQD;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowFQDNONAME->Data());
 }
 
 if(entryIntFlowLYZ1SUMNONAME)
 { 
  (*entryIntFlowLYZ1SUMNONAME)+=(Long_t)lyz1sumCountRealNONAME;
  entryIntFlowLYZ1SUMNONAME->Append(" out of ");
  (*entryIntFlowLYZ1SUMNONAME)+=(Long_t)nPointsLYZ1SUM;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowLYZ1SUMNONAME->Data());
 }
 
 if(entryIntFlowLYZ1PRODNONAME)
 { 
  (*entryIntFlowLYZ1PRODNONAME)+=(Long_t)lyz1prodCountRealNONAME;
  entryIntFlowLYZ1PRODNONAME->Append(" out of ");
  (*entryIntFlowLYZ1PRODNONAME)+=(Long_t)nPointsLYZ1PROD;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowLYZ1PRODNONAME->Data());
 }
 
 if(entryIntFlowLYZEPNONAME)
 { 
  (*entryIntFlowLYZEPNONAME)+=(Long_t)lyzepCountRealNONAME;
  entryIntFlowLYZEPNONAME->Append(" out of ");
  (*entryIntFlowLYZEPNONAME)+=(Long_t)nPointsLYZEP;
  if(textResultsNONAME)textResultsNONAME->AddText(entryIntFlowLYZEPNONAME->Data());
 }
 
 if(plotLegend)
 {
  if(textDefaultNONAME) textDefaultNONAME->Draw();
  if(textResultsNONAME) textResultsNONAME->Draw();
  if(textAvMNevtsNONAME) textAvMNevtsNONAME->Draw();
 }
}

void LoadSpreadLibraries(const libModes mode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  
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
  gSystem->Load("libTree");

  // for AliRoot
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG2flowCommon");
  cerr<<"libPWG2flowCommon loaded ..."<<endl;
  
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


