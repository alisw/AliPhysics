// add comment
 
// Set how many output analysis files in total you want to access:
const Int_t nFiles = 1;
 
// Set how many of those output analysis files you want to represent with a mesh (usually used to represent results of simulations):
const Int_t nSim = 0;
 
// Set paths of all output analysis files (first the ones to be represented with mesh (simulations), then the ones to be represented with markers (real data))
TString files[nFiles] = {"data/LHC10bcd"};

// Set analysis types for all output analysis files (can be "ESD","AOD","MC",""):
TString type[nFiles] = {"ESD"};
 
// Set mesh color:
Int_t meshColor[nSim] = {};

// Set marker styles:
Int_t markerStyle[nFiles-nSim] = {kStar};

// Set marker colors:
Int_t markerColor[nFiles-nSim] = {kBlack};

// Set legend entries:
TString legendEntry[nFiles] = {"LHC10bcd"};
 
// Set if you want to rebin the histograms into wider multiplicity bins:
Bool_t rebin = kTRUE;
Int_t nMergedBins = 10; // set how many original multiplicity bins will be merged into 1 new one
 
// Set if you whish to plot cumulants versus <reference multiplicity> (by default they are plotted versus # of RPs):
Bool_t plotCumulantsVsReferenceMultiplicity = kFALSE;

// Set flow values whose theoretical contribution to cumulants will be shown on the plots with the straight coloured lines: 
Bool_t showTheoreticalLines = kTRUE;
const Int_t nFlowValues = 1;
Double_t v[nFlowValues] = {0.1,};
Int_t lineColor[nFlowValues] = {kRed}; 

// If the statistical error of 6th and 8th order cumulant is huge you may prefer not to show them:
Bool_t plotOnly2ndAnd4thOrderCumulant = kFALSE;

// For comparison sake show also GFC results with dotted line:
Bool_t showAlsoGFCResults = kFALSE;
Int_t gfcLineStyle = 3;

// Set method names which calculate cumulants vs multiplicity:
const Int_t nMethods = 2;
TString method[nMethods] = {"QC","GFC"}; 

TFile *commonOutputFiles[nFiles] = {NULL}; // common output files "AnalysisResults.root"
TList *lists[nFiles][nMethods] = {{NULL}}; // lists cobj<method> holding objects with results for each method
TH1D *cumulantsVsM[nFiles][nMethods][4] = {{{NULL}}}; // histograms with results for cumulants vs multiplicity (4 stands for 4 cumulant orders)
TGraph *lines[nFlowValues][4] = {{NULL}}; // lines denoting theoretical flow contribution to cumulants
TProfile *refMultVsNoOfRPs[nFiles] = {NULL}; // <reference multipicity> versus # of RPs

// Ranges for plots:
Double_t xMin[4]={0.};
Double_t xMax[4]={0.};
Double_t yMin[4]={0.};
Double_t yMax[4]={0.};

enum libModes {mLocal,mLocalSource};

void plotCumulants(Int_t analysisMode=mLocal)
{
 // analysisMode: if analysisMode = mLocal -> analyze data on your computer using aliroot
 //               if analysisMode = mLocalSource -> analyze data on your computer using root + source files 
    
 // Load needed libraries:
 LoadLibrariesPC(analysisMode);
 
 // Access all common output files:
 TString commonOutputFileName = "AnalysisResults.root"; 
 AccessCommonOutputFiles(commonOutputFileName);
 
 // Get from common output files the lists holding histograms for each method:
 GetLists();
 
 // Get histograms with results for cumulants vs multiplicity:
 GetHistograms();
 
 // Determine ranges for plots:
 DetermineMinMax();
 
 // Make lines which indicate theoretical contributions of flow to cumulants:
 Lines();
 
 // Print number of events and average multiplicities for each common output file:
 Print();  
 
 // Make plots:
 Plot(); 
 
 // Global settings which will affect all plots:
 GlobalSettings();
 
} // end of void plotCumulants(Int_t analysisMode=mLocal) 
 
// =====================================================================================

void Plot()
{
 // Make all plots.
 
 TCanvas *c = NULL;
 Int_t coMax = 0;
 if(!plotOnly2ndAnd4thOrderCumulant)
 {
  c = new TCanvas("c","cumulants");
  c->Divide(2,2);
  coMax = 4; 
 } else 
   {
    c = new TCanvas("c","cumulants",1200,500);
    c->Divide(2,1);  
    coMax = 2; 
   } 
   
 TLegend *legend = new TLegend(0.1,0.7,0.33,0.9);
 legend->SetFillStyle(0);
 //legend->SetHeader("     minClustersTpcRP");

 TString qcFlag[4] = {"QC{2}","QC{4}","QC{6}","QC{8}"};
 
 for(Int_t co=0;co<coMax;co++) // cumulant order
 {
  c->cd(co+1);
  StyleHist(qcFlag[co].Data(),co)->Draw(); 
  // simulations:
  for(Int_t s=0;s<nSim;s++) 
  {
   TGraph *errorMesh = GetErrorMesh(cumulantsVsM[s][0][co]);
   if(errorMesh)
   {
    errorMesh->SetFillColor(meshColor[s]);
    errorMesh->Draw("lfsame");
   }
   if(co==0){legend->AddEntry(errorMesh,legendEntry[s].Data(),"f");}
  } // end of if(Int_t s=0;s<nSim;s++) 
  // data:
  for(Int_t f=nSim;f<nFiles;f++)
  {
   // Theoretical lines:
   if(showTheoreticalLines)
   {
    if(f==nSim) // plot them only once
    {
     for(Int_t fv=0;fv<nFlowValues;fv++)
     { 
      lines[fv][co]->Draw("lsame");
      if(co==0){legend->AddEntry(lines[fv][co],Form("v_{2} = %g",v[fv]),"l");}  
     } 
    }
   } 
   // QC results:
   if(cumulantsVsM[f][0][co])
   {
    cumulantsVsM[f][0][co]->Draw("e1same"); 
    cumulantsVsM[f][0][co]->SetMarkerStyle(markerStyle[f-nSim]);  
    cumulantsVsM[f][0][co]->SetMarkerColor(markerColor[f-nSim]); 
    if(co==0)
    {
     if(showAlsoGFCResults)
     {
      legend->AddEntry(cumulantsVsM[f][0][co],Form("%s (QC)",legendEntry[f].Data()),"p");
     } else
       {
        legend->AddEntry(cumulantsVsM[f][0][co],legendEntry[f].Data(),"p");     
       }
    }
   }
   // GFC results:
   if(showAlsoGFCResults && cumulantsVsM[f][1][co])
   {
    cumulantsVsM[f][1][co]->Draw("lsame");  
    cumulantsVsM[f][1][co]->SetLineStyle(gfcLineStyle);    
    if(co==0){legend->AddEntry(cumulantsVsM[f][1][co],Form("%s (GFC)",legendEntry[f].Data()),"l");}
   }
  } 
  // Draw legend:
  if(co==0){legend->Draw("same");}
 } // end of for(Int_t co=0;co<4;co++) // cumulant order
 
 // Plot also <reference multiplicity> vs # of RPs:
 if(plotCumulantsVsReferenceMultiplicity)
 {
  TCanvas *cRefMultVsNoOfRPs = new TCanvas("cRefMultVsNoOfRPs","#LTreference multiplicity#GT vs # of RPs",1200,600);
  cRefMultVsNoOfRPs->Divide(nFiles,1);
  for(Int_t f=0;f<nFiles;f++)
  {
   cRefMultVsNoOfRPs->cd(f+1);
   if(refMultVsNoOfRPs[f])
   {
    refMultVsNoOfRPs[f]->SetTitle(legendEntry[f].Data());
    refMultVsNoOfRPs[f]->GetXaxis()->SetRangeUser(0,refMultVsNoOfRPs[f]->FindLastBinAbove());
    refMultVsNoOfRPs[f]->Draw();   
   } else
     {
      cout<<endl;
      cout<<"WARNING: refMultVsNoOfRPs[f] is NULL in Plot(), f = "<<f<<" !!!!"<<endl;
      cout<<endl;
     }
  }   
 }   
 
} // end of void Plot()
 
// =====================================================================================

void Lines()
{
 // Make lines denoting theoretical contribution of flow to cumulants.
 
 for(Int_t fv=0;fv<nFlowValues;fv++)
 {
  lines[fv][0] = new TGraph(2);
  lines[fv][0]->SetPoint(0,xMin[0],pow(v[fv],2));  
  lines[fv][0]->SetPoint(1,xMax[0]+0.5,pow(v[fv],2));  
  lines[fv][0]->SetLineColor(lineColor[fv]);
  lines[fv][1] = new TGraph(2);
  lines[fv][1]->SetPoint(0,xMin[1],-pow(v[fv],4));  
  lines[fv][1]->SetPoint(1,xMax[1]+0.5,-pow(v[fv],4));  
  lines[fv][1]->SetLineColor(lineColor[fv]); 
  lines[fv][2] = new TGraph(2);
  lines[fv][2]->SetPoint(0,xMin[2],4.*pow(v[fv],6));  
  lines[fv][2]->SetPoint(1,xMax[2]+0.5,4.*pow(v[fv],6));  
  lines[fv][2]->SetLineColor(lineColor[fv]);  
  lines[fv][3] = new TGraph(2);
  lines[fv][3]->SetPoint(0,xMin[3],-33.*pow(v[fv],8));  
  lines[fv][3]->SetPoint(1,xMax[3]+0.5,-33.*pow(v[fv],8));  
  lines[fv][3]->SetLineColor(lineColor[fv]);
 }

} // end of void Lines()

// =====================================================================================

void Print()
{
 // Print number of events and average multiplicities for each common output file.
 
 cout<<endl;
 cout<<"Accessed files:"<<endl;
 cout<<endl;
 for(Int_t f=0;f<nFiles;f++)
 {
  cout<<commonOutputFiles[f]->GetName()<<endl;
  for(Int_t m=0;m<nMethods;m++)
  {
   AliFlowCommonHist *commonHist = NULL;
   if(lists[f][m])
   {
    commonHist = dynamic_cast<AliFlowCommonHist*> (lists[f][m]->FindObject(Form("AliFlowCommonHist%s",method[m].Data())));
   }
   Double_t nEvts = 0.;
   Double_t AvM = 0.;
   if(commonHist && commonHist->GetHistMultRP())
   {
    nEvts = commonHist->GetHistMultRP()->GetEntries();
    AvM = commonHist->GetHistMultRP()->GetMean();
   }
   if(!(strcmp(method[m].Data(),"QC")))
   {
    cout<<Form("%s:",method[m].Data())<<"  <M> = "<<AvM<<", N = "<<nEvts<<endl;
   }
   if(!(strcmp(method[m].Data(),"GFC")) && showAlsoGFCResults)
   {
    cout<<Form("%s:",method[m].Data())<<"  <M> = "<<AvM<<", N = "<<nEvts<<endl;
   }
  } // end of for(Int_t m=0;m<nMethods;m++)
  cout<<endl;
 } // end of for(Int_t f=0;f<nFiles;f++) 

} // end of void Print()

// =====================================================================================

void DetermineMinMax()
{
 // Determine ranges for plots.
 
 for(Int_t co=0;co<4;co++)
 {
  xMin[co] = 0.; yMin[co] = 44.;
  xMax[co] = -440000.; yMax[co] = -44.;
 }
 
 Double_t tfc[nFlowValues][4] = {{0.}}; // theoretical flow contribution
 for(Int_t fv=0;fv<nFlowValues;fv++)
 {
  tfc[fv][0] = pow(v[fv],2);
  tfc[fv][1] = -pow(v[fv],4);
  tfc[fv][2] = 4.*pow(v[fv],6);
  tfc[fv][3] = -33.*pow(v[fv],8);
 }
  
 for(Int_t f=0;f<nFiles;f++)
 {
  for(Int_t m=0;m<nMethods;m++)
  { 
   for(Int_t co=0;co<4;co++)
   { 
    if(cumulantsVsM[f][m][co]) 
    {
     Int_t nBins = cumulantsVsM[f][m][co]->GetXaxis()->GetNbins();
     for(Int_t b=1;b<=nBins;b++)
     {
      Double_t result = cumulantsVsM[f][m][co]->GetBinContent(b);
      Double_t error = cumulantsVsM[f][m][co]->GetBinError(b);
      if(TMath::Abs(result)>1.e-44 && TMath::Abs(error)>1.e-44)
      {
       // y-axis:
       if(yMin[co] > result-error){yMin[co] = result-error;} // min value
       if(yMax[co] < result+error) {yMax[co] = result+error;} // max value    
       // x-axis:
       xMax[co] = cumulantsVsM[f][m][co]->GetBinLowEdge(b+1); 
      }
     } // end of for(Int_t b=1;b<=cumulantsVsM[f][m][co]->GetXaxis()->GetNbins();b++) 
     // theoretical contributions:
     if(showTheoreticalLines)
     {
      for(Int_t fv=0;fv<nFlowValues;fv++)
      {
       if(yMin[co] > tfc[fv][co]) {yMin[co] = tfc[fv][co];} // min value
       if(yMax[co] < tfc[fv][co]) {yMax[co] = tfc[fv][co];} // max value      
      } // end of for(Int_t fv=0;fv<nFlowValues;fv++)
     } // end of if(showTheoreticalLines) 
    } // end of if(cumulantsVsM[f][m][co])
   } // end of for(Int_t co=0;co<4;co++)
  } // end of for(Int_t m=0;m<nMethods;m++)
 } // end of for(Int_t f=0;f<nFiles;f++)
 
} // end of void DetermineMinMax()

// =====================================================================================

void GetHistograms()
{
 // Get histograms with results for cumulants vs multiplicity.
 
 if(plotCumulantsVsReferenceMultiplicity) {GetProfileForRefMultVsNoOfRPs();}
  
 TString qcFlag[4] = {"QC{2}","QC{4}","QC{6}","QC{8}"};
 TString gfcFlag[4] = {"GFC{2}","GFC{4}","GFC{6}","GFC{8}"};
 for(Int_t f=0;f<nFiles;f++)
 {
  for(Int_t m=0;m<nMethods;m++)
  { 
   TList *temp = NULL;
   if(!(strcmp(method[m].Data(),"QC")) && lists[f][m])
   {
    temp = dynamic_cast<TList*> (lists[f][m]->FindObject("Integrated Flow"));
    if(temp) {temp = dynamic_cast<TList*> (temp->FindObject("Results"));}
    if(temp) 
    {
     for(Int_t co=0;co<4;co++)
     {
      cumulantsVsM[f][m][co] = dynamic_cast<TH1D*> (temp->FindObject(Form("fIntFlowQcumulantsVsM, %s",qcFlag[co].Data())));
      if(rebin && cumulantsVsM[f][m][co])
      {
       cumulantsVsM[f][m][co] = Rebin(cumulantsVsM[f][m][co]);
      }
      if(plotCumulantsVsReferenceMultiplicity && cumulantsVsM[f][m][co])
      {
       Map(cumulantsVsM[f][m][co],f);
      }    
     } 
    } 
   } // end of if(!(strcmp(method[m].Data(),"QC")))
   else if(!(strcmp(method[m].Data(),"GFC")) && lists[f][m])
   {
    temp = dynamic_cast<TList*> (lists[f][m]->FindObject("Reference Flow"));
    if(temp) {temp = dynamic_cast<TList*> (temp->FindObject("Results"));}
    if(temp) 
    {
     for(Int_t co=0;co<4;co++)
     {
      cumulantsVsM[f][m][co] = dynamic_cast<TH1D*> (temp->FindObject(Form("fReferenceFlowCumulantsVsM, %s",gfcFlag[co].Data())));
      if(plotCumulantsVsReferenceMultiplicity && cumulantsVsM[f][m][co])
      {
       Map(cumulantsVsM[f][m][co],f);
      }
     } // end of for(Int_t co=0;co<4;co++)
    } // end of if(temp)
   } // end of else if(!(strcmp(method[m].Data(),"QC")))
  } // end of  for(Int_t m=0;m<nMethods;m++)
 } // end of for(Int_t f=0;f<nFiles;f++)

 Int_t counter = 0;
 for(Int_t f=0;f<nFiles;f++)
 {
  for(Int_t m=0;m<nMethods;m++)
  { 
   for(Int_t co=0;co<4;co++)
   {
    if(cumulantsVsM[f][m][co]){counter++;}
   }
  }
 }
 
 if(counter == 0)
 {
  cout<<endl;
  cout<<" WARNING: Couldn't access a single histogram with results vs multiplicity !!!!"<<endl;
  cout<<"          Did you enable this calculation before running over data?"<<endl;
  cout<<endl;
  exit(0);
 }   

} // end of void GetHistograms()

// =====================================================================================

void GetProfileForRefMultVsNoOfRPs()
{
 // Get profile holding <reference multiplicity> versus # of RPs.
 
 for(Int_t f=0;f<nFiles;f++)
 {
  AliFlowCommonHist *commonHist = NULL;
  if(lists[f][0])
  {
   commonHist = dynamic_cast<AliFlowCommonHist*> (lists[f][0]->FindObject("AliFlowCommonHistQC"));
  } else
    {
     cout<<endl;
     cout<<"WARNING: lists[f][0] is NULL in GetProfileForRefMultVsNoOfRPs(), f = "<<f<<" !!!!"<<endl;
     cout<<endl;    
    } 
  if(commonHist && commonHist->GetRefMultVsNoOfRPs())
  {
   refMultVsNoOfRPs[f] = commonHist->GetRefMultVsNoOfRPs();
  } else
    {
     cout<<endl;
     cout<<"WARNING: commonHist && commonHist->GetRefMultVsNoOfRPs() is NULL in GetProfileForRefMultVsNoOfRPs() !!!!"<<endl;
     cout<<endl;
    }
 }   
 
} // end of void GetProfileForRefMultVsNoOfRPs()

// =====================================================================================

void Map(TH1D *hist, Int_t f)
{
 // Map cumulant versus # of RPs into cumulants versus <reference multiplicity>.
 
 TH1D *temp = NULL;
 if(!refMultVsNoOfRPs[f])
 {
  cout<<endl;
  cout<<"WARNING: refMultVsNoOfRPs[f] is NULL in Map(...), f = "<<f<<" !!!!"<<endl;
  cout<<endl;
 }
 
 if(hist)
 {
  temp = (TH1D*) hist->Clone();
  Int_t nBins = hist->GetNbinsX();
  for(Int_t b=1;b<=nBins;b++)
  {
   hist->SetBinContent((Int_t)refMultVsNoOfRPs[f]->GetBinContent(b),temp->GetBinContent(b));
   hist->SetBinError((Int_t)refMultVsNoOfRPs[f]->GetBinContent(b),temp->GetBinError(b));
  }
 }
     
} // end of Map()

// =====================================================================================

TH1D* Rebin(TH1D *hist)
{
 // Rebin original histograms.
 
 if(nMergedBins == 0)
 {
  cout<<endl;
  cout<<" WARNING: nMergedBins == 0 !!!!"<<endl;
  cout<<endl;
  exit(0);
 } 
 if(!hist)
 {
  cout<<endl;
  cout<<" WARNING: hist is NULL in Rebin() !!!!"<<endl;
  cout<<endl;
  exit(0); 
 } 
 
 Int_t nBinsOld = hist->GetXaxis()->GetNbins(); 
 if(nBinsOld == 0){cout<<" WARNING: nBinsOld == 0 !!!!"<<endl;exit(0);}
 Double_t xMinOld = hist->GetXaxis()->GetXmin(); 
 Double_t xMaxOld = hist->GetXaxis()->GetXmax(); 
 Double_t binWidthOld = (xMaxOld-xMinOld)/nBinsOld;
 Int_t nBinsNew = TMath::Floor(nBinsOld/nMergedBins);
 Double_t xMinNew = xMinOld;
 Double_t xMaxNew = xMinOld + nBinsNew*nMergedBins*binWidthOld;
  
 TH1D *temp = new TH1D("","",nBinsNew,xMinNew,xMaxNew); // rebinned histogram 
 Int_t binNew = 1;
 Double_t value = 0.;
 Double_t error = 0.;
 Double_t dSum1 = 0.; // sum value_i/(error_i)^2
 Double_t dSum2 = 0.; // sum 1/(error_i)^2
 for(Int_t b=1;b<=nBinsOld;b++)
 {
  value = hist->GetBinContent(b);  
  error = hist->GetBinError(b);  
  if(error>0.)
  {
   dSum1+=value/(error*error);
   dSum2+=1./(error*error);
  }
  if(b%nMergedBins == 0)
  {
   if(dSum2>0.)
   {
    temp->SetBinContent(binNew,dSum1/dSum2);
    temp->SetBinError(binNew,pow(1./dSum2,0.5));
   }
   binNew++;
   dSum1 = 0.;
   dSum2 = 0.;
  } // end of if(b%nMergedBins == 0)
 } // end of for(Int_t b=1;b<=nBinsOld;b++)
  
 return temp;
      
} // end of Rebin()

// =====================================================================================

TGraphErrors* GetGraphErrors(Int_t bin, Int_t nFiles, TH1D** qc)
{
 TGraphErrors *ge = new TGraphErrors(nFiles);
 for(Int_t f=0;f<nFiles;f++)
 {
  ge->SetPoint(f,f+0.5,qc[f]->GetBinContent(bin+1));
  ge->SetPointError(f,0,qc[f]->GetBinError(bin+1));
 }

 return ge;

} // end of TGraphErrors* GetGraphErrors(Int_t bin, Int_t nFiles, TH1D** qc)

// =====================================================================================

TH1D* StyleHist(TString yAxisTitle, Int_t co) 
{
 // Style histogram.
 
 TH1D *styleHist = new TH1D(Form("%d",co),"",(Int_t)xMax[co],0,xMax[co]);
 // y-axis:
 styleHist->GetYaxis()->SetRangeUser(yMin[co],yMax[co]);
 
 if(plotCumulantsVsReferenceMultiplicity)
 {
  styleHist->GetXaxis()->SetTitle("#LTreference multiplicity#GT");     
 } else
   {
    styleHist->GetXaxis()->SetTitle("# of RPs");         
   }
 
 styleHist->GetYaxis()->SetTitle(yAxisTitle.Data());
 
 return styleHist;

} // end of TH1D* StyleHist(TString yAxisTitle, Int_t co)  
  
// ===========================================================================================

TGraph* GetErrorMesh(TH1D *hist)
{
 // Error mesh.
 
 TGraph* errorMesh = NULL;
 if(hist)
 {
  Int_t nBins = hist->GetNbinsX();
  Double_t binWidth = hist->GetBinWidth(1); // assuming that all bins have the same width
  // Counting the non-empty bins: 
  Int_t nNonEmptyBins=0;
  for(Int_t i=1;i<nBins+1;i++)
  {
   if(!(hist)->GetBinError(i)==0.0))
   {
    nNonEmptyBins++;
   }
  }   
  errorMesh = new TGraph(2*nNonEmptyBins+1); 
  Double_t value=0.,error=0.;
  Int_t count=1;
  Double_t xFirst=0.,yUpFirst=0.; // needed to close up the mesh
  for(Int_t i=1;i<nBins+1;i++)
  {
   // Setting up the upper limit of the mesh:
   value = hist->GetBinContent(i);
   error = hist->GetBinError(i);   
   if(!(error==0.0))
   {    
    errorMesh->SetPoint(count++,(i-0.5)*binWidth,value+error);
    if(xFirst==0.)
    {
     xFirst=(i-0.5)*binWidth;
     yUpFirst=value+error;
    }
   } 
  }   
  for(Int_t i=nBins+1;i<2*nBins+1;i++)
  {
   // Setting up the lower limit of the mesh:
   value = hist->GetBinContent(2*nBins+1-i);
   error = hist->GetBinError(2*nBins+1-i); 
   if(!(error==0.0))
   {      
    errorMesh->SetPoint(count++,(2*nBins-i+0.5)*binWidth,value-error);
   }  
  }
  // Closing the mesh area:
  errorMesh->SetPoint(2*nNonEmptyBins+1,xFirst,yUpFirst);   
 } // end if(hist)
 
 errorMesh->SetFillStyle(1001);
 
 return errorMesh;
 
} // end of TGraph* GetErrorMesh(TH1D *hist)

// ===========================================================================================

void GlobalSettings()
{
 // Settings which will affect all plots.
 
 gROOT->SetStyle("Plain"); // default color is white instead of gray
 gStyle->SetOptStat(0); // remove stat. box from all histos
 TGaxis::SetMaxDigits(4); 
 
} // end of void GlobalSettings()

// ===========================================================================================

void GetLists() 
{
 // Get from common output files the lists holding histograms for each method.

 TString fileName[nFiles][nMethods]; 
 TDirectoryFile *dirFile[nFiles][nMethods] = {{NULL}}; 
 TString listName[nFiles][nMethods]; 
 for(Int_t f=0;f<nFiles;f++)
 { 
  for(Int_t i=0;i<nMethods;i++)
  {
   // Form a file name for each method:
   fileName[f][i]+="output";
   fileName[f][i]+=method[i].Data();
   fileName[f][i]+="analysis";
   fileName[f][i]+=type[f].Data();
   // Access this file:
   if(commonOutputFiles[f]){dirFile[f][i] = (TDirectoryFile*)commonOutputFiles[f]->FindObjectAny(fileName[f][i].Data());}
   // Form a list name for each method:
   listName[f][i]+="cobj";
   listName[f][i]+=method[i].Data();
   // Access this lists:
   if(dirFile[f][i])
   {
    dirFile[f][i]->GetObject(listName[f][i].Data(),lists[f][i]); 
    if(!lists[f][i])
    {
     cout<<endl;
     cout<<" WARNING: Couldn't find a list "<<listName[f][i].Data()<<" in "<<commonOutputFiles[f]->GetName()<<" !!!!"<<endl;
     cout<<"          Did you use method "<<method[i].Data()<<" in the analysis which produced this file?"<<endl;
     cout<<endl;
    }
   } else 
     {
      cout<<endl;
      cout<<" WARNING: Couldn't find a file "<<fileName[f][i].Data()<<".root in "<<commonOutputFiles[f]->GetName()<<" !!!!"<<endl;exit(0);
      cout<<endl;
     }
  } // end of for(Int_t i=0;i<nMethods;i++)   
 } // end of for(Int_t f=0;f<nFiles;f++)

} // end of void GetLists() 

// ===========================================================================================

void AccessCommonOutputFiles(TString commonOutputFileName)
{
 // Access all output files.
 for(Int_t f=0;f<nFiles;f++)
 { 
  if(!(gSystem->AccessPathName(Form("%s/%s/%s",gSystem->pwd(),files[f].Data(),commonOutputFileName.Data()),kFileExists)))
  {
   commonOutputFiles[f] = TFile::Open(Form("%s/%s/%s",gSystem->pwd(),files[f].Data(),commonOutputFileName.Data()),"READ");
  } else
    { 
     cout<<endl;
     cout<<"WARNING: Couldn't find the file "<<Form("%s/%s/%s",gSystem->pwd(),files[f].Data(),commonOutputFileName.Data())<<" !!!!"<<endl;
     cout<<endl;
     exit(0);
    }
 } // end of for(Int_t f=0;f<nFiles;f++)
 
} // void AccessCommonOutputFiles(TString commonOutputFileName);

// ===========================================================================================

void LoadLibrariesPC(const libModes analysisMode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  //gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  
  //----------------------------------------------------------
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  if (analysisMode==mLocal) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------
    
    //==================================================================================  
    //load needed libraries:
    gSystem->AddIncludePath("-I$ROOTSYS/include");
    //gSystem->Load("libTree");
    
    // for AliRoot
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    //gSystem->Load("libANALYSIS");
    gSystem->Load("libPWG2flowCommon");
    //cerr<<"libPWG2flowCommon loaded ..."<<endl;
    
  }
  
  else if (analysisMode==mLocalSource) {
    
    // In root inline compile
   
    // Constants  
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonConstants.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZConstants.cxx+");
    
    // Flow event
    gROOT->LoadMacro("AliFlowCommon/AliFlowVector.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimple.cxx+");    
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimpleCuts.cxx+");    
    gROOT->LoadMacro("AliFlowCommon/AliFlowEventSimple.cxx+");
   
    // Output histosgrams
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHist.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHistResults.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist1.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist2.cxx+");
    
    cout << "finished loading macros!" << endl;  
    
  }  
  
} // end of void LoadLibraries(const libModes analysisMode) 
