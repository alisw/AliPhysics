// Macro plotCumulants.C is used to represent true correlations, a.k.a. cumulants,
// as a function of reference multiplicity. Different order cumulants provide
// independent estimates of the flow harmonics. In the ideal case when only flow
// correlations are present we have: 
//         c{2} = v^2, c{4} = -v^4, c{6} = 4v^6 and c{8} = -33v^8.
// In practice cumulants can are calculated via Q-vectors (QC) or via formalism
// of generating functions (GFC). 

// Set how many output analysis files in total you want to access:
const Int_t nFiles = 2;
 
// Set how many of those output analysis files you want to represent with a mesh (usually used to represent results of simulations):
const Int_t nSim = 0;

// Set paths of all output analysis files (first the ones to be represented with mesh (simulations), then the ones to be represented with markers (real data))
TString files[nFiles] = {"trackletsCorrectedNUA","default"}; // subdirectory names holding <commonOutputFileName>.root
//TString files[nFiles] = {"outputCentrality0","outputCentrality1"}; // file names (use case: centrality train)

// Set analysis types for all output analysis files (can be "ESD","AOD","MC","","MK", ....):
TString type[nFiles] = {"ESD","ESD"};
//TString type[nFiles] = {"MK","MK"};
 
// Set mesh color:
Int_t meshColor[nSim] = {};

// Set marker styles:
Int_t markerStyle[nFiles-nSim] = {kFullSquare,kOpenSquare};

// Set marker colors:
Int_t markerColor[nFiles-nSim] = {kBlack,kRed};

// Set legend entries:
TString legendEntry[nFiles] = {"",""};
 
// Set if you want to rebin the histograms into wider multiplicity bins (set for each cumulant order separately):
Bool_t rebin = kFALSE;
Int_t nMergedBins[4] = {10,10,10,10}; // set how many original multiplicity bins will be merged into 1 new one 
 
// Set if you whish to plot cumulants versus <reference multiplicity> (by default they are plotted versus # of RPs):
Bool_t plotCumulantsVsReferenceMultiplicity = kFALSE;
Bool_t showReferenceMultiplicityVsNoOfRPs = kFALSE; 

// Set flow values whose theoretical contribution to cumulants will be shown on the plots with the straight coloured lines: 
Bool_t showTheoreticalLines = kFALSE;
const Int_t nFlowValues = 1;
Double_t v[nFlowValues] = {0.05};
Int_t lineColor[nFlowValues] = {kRed}; 

// If the statistical error of 6th and 8th order cumulant is huge you may prefer not to show them:
Bool_t plotOnly2ndAnd4thOrderCumulant = kFALSE;

// Set if you want independent canvas with results for reference flow vs multiplicity:
Bool_t showRefFlowVsM = kTRUE;
Int_t refFlowVsMMarkerStyle[nFiles] = {kFullSquare,kOpenSquare}; // marker style is different for different file
Int_t refFlowVsMMarkerColor[4] = {kBlack,kRed,kBlue,kGreen+2}; // marker color is different for different cumulant order
Int_t refFlowVsMeshOrder = -1; // set results of which order will be plotted as mesh in all pads [0=2nd, 1=4th, 2=6th, 3=8th, -1=do not show]
Int_t refFlowVsMMeshColor = kRed-10; // mesh color for above specified order
Double_t refFlowVsMxRange[2] = {1.,11000.}; // x range on the plots for reference multiplicity vs M
Double_t refFlowVsMyRange[2] = {0.0,0.194}; // x range on the plots for reference multiplicity vs M

// Set if you want to show theoretical curves for the toy model:
Bool_t showToyModel = kFALSE;
const Int_t nToyModels = 2; // number of toy models with different parameters k, vn, v2n
Double_t k[nToyModels] = {2.,4.};
Double_t vn[nToyModels] = {0.25,0.25};
Double_t v2n[nToyModels] = {0.0,0.0};
Int_t lineColorToyModel[nToyModels] = {kBlack,kRed};

// For comparison sake show also GFC results with dotted line:
Bool_t showAlsoGFCResults = kFALSE;
Int_t gfcLineStyle = 3;

// For comparison sake show also Monte Carlo QC results with coloured mesh:
Bool_t showAlsoMCEPResults = kFALSE; 
Bool_t showOnlyMCEPResults = kFALSE;
Int_t mcepMeshColor[nSim] = {};

// Set the naming convention:
Bool_t bLookInSubdirectories = kFALSE; // if kTRUE: Look in subdirectories <files[0]>, <files[1]>, ..., for output files <commonOutputFileName>
                                      // if kFALSE: Look in <pwd> for files <files[0]>.root, <files[1]>.root, ....
TString commonOutputFileName = "AnalysisResults.root"; 

// Set method names which calculate cumulants vs multiplicity (do not touch these settings unless you are looking for a trouble):
const Int_t nMethods = 3;
TString method[nMethods] = {"QC","GFC","MCEP"}; 

TFile *commonOutputFiles[nFiles] = {NULL}; // common output files "AnalysisResults.root"
TList *lists[nFiles][nMethods] = {{NULL}}; // lists cobj<method> holding objects with results for each method
TH1D *cumulantsVsM[nFiles][nMethods][4] = {{{NULL}}}; // histograms with results for cumulants vs multiplicity (4 stands for 4 cumulant orders)
TH1D *refFlowVsM[nFiles][nMethods][4] = {{{NULL}}}; // histograms with results for reference flow vs multiplicity (4 stands for 4 cumulant orders)
TGraph *lines[nFlowValues][4] = {{NULL}}; // lines denoting theoretical flow contribution to cumulants
TProfile *refMultVsNoOfRPs[nFiles] = {NULL}; // <reference multipicity> versus # of RPs
TF1 *toyModel[2][nToyModels] = {{NULL}}; // [cumulant order][number of toy models]

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
 
 // Cross-check user settings:
 CrossCheckSettings();
    
 // Load needed libraries:
 LoadLibrariesPC(analysisMode);
 
 // Global settings which will affect all plots:
 GlobalSettings();
  
 // Access all common output files:
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
 
 // Make plots for cumulants vs multiplicity: // to be improved
 PlotCumulantsVsM(); 

 // Make plots for reference flow vs multiplicity:  
 if(showRefFlowVsM){PlotRefFlowVsM();}

} // end of void plotCumulants(Int_t analysisMode=mLocal)  

// =====================================================================================

void PlotRefFlowVsM()
{
 // Make plots for reference flow vs multiplicity.
 
 // Calculate reference flow from cumulants vs M:
 CalculateReferenceFlowFromCumulantsVsM();
 
 TCanvas *cRefFlowVsM = new TCanvas("cRefFlowVsM","Reference Flow");
 cRefFlowVsM->Divide(2,2);

 TLegend *lRefFlowVsM = new TLegend(0.1,0.7,0.33,0.9);
 lRefFlowVsM->SetFillStyle(0);
 lRefFlowVsM->SetHeader("     Therminator (20-30%)");

 TString refFlowVsMFlag[4] = {"v_{2}{2,QC}","v_{2}{4,QC}","v_{2}{6,QC}","v_{2}{8,QC}"};

 for(Int_t f=0;f<nFiles;f++)
 {
  for(Int_t m=0;m<1+(Int_t)(showAlsoGFCResults)+(Int_t)(showAlsoMCEPResults);m++)
  {
   for(Int_t co=0;co<4;co++) // cumulant order
   {
    cRefFlowVsM->cd(co+1);
    // Style histogram:
    if(f==0) // superimpose histograms from other files on top of this one
    {
     TH1D *styleHist = (TH1D*)StyleHist(refFlowVsMFlag[co].Data(),co)->Clone();
     if(styleHist)
     {
      styleHist->GetXaxis()->SetRangeUser(refFlowVsMxRange[0],refFlowVsMxRange[1]);
      styleHist->GetYaxis()->SetRangeUser(refFlowVsMyRange[0],refFlowVsMyRange[1]);
      styleHist->Draw();
     }
    } // end of if(f==0)    
    // Plot first the mesh for reference flow of specified order in all pads:
    if(refFlowVsMeshOrder != -1)
    {
     TGraph *rfMeshOutOf = GetErrorMesh(refFlowVsM[f][m][refFlowVsMeshOrder]);
     if(rfMeshOutOf)
     {
      rfMeshOutOf->SetFillColor(refFlowVsMMeshColor);
      rfMeshOutOf->Draw("lfsame"); 
      if(refFlowVsMeshOrder==co){lRefFlowVsM->AddEntry(rfMeshOutOf,Form("v_{2}{%i,QC} stat. error",2*(refFlowVsMeshOrder+1)),"f");}
     }
    } // end of if(refFlowVsMeshOrder != -1)
    // Plot reference flow vs M:
    refFlowVsM[f][m][co]->SetMarkerStyle(refFlowVsMMarkerStyle[f]);
    refFlowVsM[f][m][co]->SetMarkerColor(refFlowVsMMarkerColor[co]);
    refFlowVsM[f][m][co]->Draw("e1same");
    lRefFlowVsM->AddEntry(refFlowVsM[f][m][co],Form("v_{2}{%i,QC}, %s",2*(co+1),files[f].Data()),"p"); 
   } // end of for(Int_t co=0;co<4;co++) // cumulant order
  } // end of for(Int_t m=0;m<nMethods;m++)
 } // end of for(Int_t f=0;f<nFiles;f++)   
 
 // Draw a common legend in the 1st pad:
 cRefFlowVsM->cd(1);
 lRefFlowVsM->Draw("same");

} // end of void PlotRefFlowVsM()

// =====================================================================================

void PlotCumulantsVsM()
{
 // Make plots for cumulants vs multiplicity: // to be improved
 
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
  TH1D *styleHist = (TH1D*) StyleHist(qcFlag[co].Data(),co)->Clone(); 
  if(styleHist){styleHist->Draw();}
  if(co==0)
  {
   if(showToyModel)
   {
    for(Int_t ntm=0;ntm<nToyModels;ntm++)
    {
     TF1 *tm = ToyModel(co,k[ntm],vn[ntm],v2n[ntm]);
     if(tm)
     {
      tm->SetLineColor(lineColorToyModel[ntm]);
      tm->Draw("same");
     }
    } // for(Int_t ntm=0;ntm<nToyModels;ntm++)
   } // end of if(showToyModel)
  } // end of if(co==0)
  if(co==1)
  {    
   if(showToyModel)
   {
    for(Int_t ntm=0;ntm<nToyModels;ntm++)
    {
     TF1 *tm = ToyModel(co,k[ntm],vn[ntm],v2n[ntm]);
     if(tm)
     {
      tm->SetLineColor(lineColorToyModel[ntm]);
      tm->Draw("same");
     }
    } // for(Int_t ntm=0;ntm<nToyModels;ntm++)
   } // end of if(showToyModel)   
  } // end of if(co==1)
  
  // simulations:
  for(Int_t s=0;s<nSim;s++) 
  {
   // Monte Carlo QC:
   if(showAlsoMCEPResults)
   {
    TGraph *mcepQC = GetErrorMesh(cumulantsVsM[s][2][co]);    
    if(mcepQC)
    {
     mcepQC->SetFillColor(mcepMeshColor[s]);
     mcepQC->Draw("lfsame");
    }
    if(co==0){legend->AddEntry(mcepQC,Form("%s (MC)",legendEntry[s].Data()),"f");}   
   } // end of if(showAlsoMCEPResults)  
   // QC:
   TGraph *errorMesh = GetErrorMesh(cumulantsVsM[s][0][co]);
   if(errorMesh)
   {
    errorMesh->SetFillColor(meshColor[s]);
    if(!(showOnlyMCEPResults && showAlsoMCEPResults)){errorMesh->Draw("lfsame");}
   }
   if(co==0 && !(showOnlyMCEPResults && showAlsoMCEPResults)){legend->AddEntry(errorMesh,legendEntry[s].Data(),"f");}  
  } // end of if(Int_t s=0;s<nSim;s++) 
  // Theoretical lines:
  if(showTheoreticalLines)
  {
   for(Int_t fv=0;fv<nFlowValues;fv++)
   { 
    lines[fv][co]->Draw("lsame");
    if(co==0){legend->AddEntry(lines[fv][co],Form("v_{2} = %g",v[fv]),"l");}  
   } 
  } 
  // data:
  for(Int_t f=nSim;f<nFiles;f++)
  {
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
 if(plotCumulantsVsReferenceMultiplicity && showReferenceMultiplicityVsNoOfRPs)
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
 
} // end of void PlotCumulantsVsM()
 
// =====================================================================================

void CrossCheckSettings()
{
 // Cross-check user settings in this method.
 
 if(showAlsoGFCResults && rebin)
 {
  cout<<endl;
  cout<<" WARNING: Rebinning in M not supported for GFC yet !!!!"<<endl;
  cout<<endl;
  exit(0);
 }

 if(showAlsoGFCResults && plotCumulantsVsReferenceMultiplicity)
 {
  cout<<endl;
  cout<<" WARNING: Showing GFC versus <reference multiplicity> not supported yet !!!!"<<endl;
  cout<<endl;
  exit(0);
 }

} // end of void CrossCheckSettings()

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

TF1* ToyModel(Int_t co, Double_t k, Double_t vn, Double_t v2n)
{
 // Make theoretical curves for the toy model.
 
 TF1 *function = NULL;
 
 if(co==0) 
 {
  function = new TF1("function","([1]*[1]*(x-0.5-[0])+[0]-1)/((x-0.5)-1)",1.5,1000);
  function->SetParameter(0,k);
  function->SetParameter(1,vn);
 } 
 else if(co==1)
 {
  function = new TF1("function","-1.*([1]*[1]*[1]*[1]*([0]-(x-0.5))*(12.*[0]-6.*[0]*[0]-12.*(x-0.5)-5.*[0]*(x-0.5)+6.*[0]*[0]*(x-0.5)+9.*(x-0.5)*(x-0.5)-3.*[0]*(x-0.5)*(x-0.5)-(x-0.5)*(x-0.5)*(x-0.5))+[1]*[1]*4.*([0]-1.)*(x-0.5-[0])*([0]*(x-0.5-1.)-2.*(x-0.5-2.))-[1]*[1]*[2]*2.*(x-0.5-1.)*([0]-1.)*(x-0.5-[0])*(x-0.5-2.*[0])-[2]*[2]*([0]-1.)*([0]-1.)*((x-0.5)-[0])*(x-0.5-1.)+([0]-1.)*(-6.+9*[0]-[0]*[0]+2.*(x-0.5)-5.*[0]*(x-0.5)+[0]*[0]*(x-0.5)))/((x-0.5-1.)*(x-0.5-1.)*(x-0.5-2.)*(x-0.5-3.))",3.5,100);         
  function->SetParameter(0,k);
  function->SetParameter(1,vn);
  function->SetParameter(2,v2n); 
 }
 
 return function;
 
} // end of TF1* ToyModel(Int_t co)

// =====================================================================================

void CalculateReferenceFlowFromCumulantsVsM()
{
 // Calculate reference flow from cumulants vs M:
 
 for(Int_t f=0;f<nFiles;f++)
 {
  for(Int_t m=0;m<1+(Int_t)(showAlsoGFCResults)+(Int_t)(showAlsoMCEPResults);m++)
  {
   for(Int_t co=0;co<4;co++) // cumulant order
   {
    if(cumulantsVsM[f][m][co])
    {
     refFlowVsM[f][m][co] = (TH1D*)cumulantsVsM[f][m][co]->Clone(Form("%i,%i,%i",f,m,co));
     if(!refFlowVsM[f][m][co]){cout<<" WARNING: "<<Form("refFlowVsM[%i][%i][%i]",f,m,co)<<" is NULL in PlotRefFlowVsM() !!!!"<<exit(0)<<endl;}
     refFlowVsM[f][m][co]->Reset(); // to have empty histogram with the same binning as the one with cumulants
     Int_t nBins = refFlowVsM[f][m][co]->GetNbinsX();  
     for(Int_t b=1;b<=nBins;b++)
     {
      Double_t qcVsM = cumulantsVsM[f][m][co]->GetBinContent(b); // QC vs M 
      Double_t qcErrorVsM = cumulantsVsM[f][m][co]->GetBinError(b); // QC vs M stat. error  
      Double_t vVsM = 0.; // reference flow vs M  
      Double_t vErrorVsM = 0.; // reference flow vs M stat. error  
      if(co==0) // 2nd order
      {
       if(qcVsM>0.)
       {
        vVsM = pow(qcVsM,1./2.);
        vErrorVsM = (1./2.)*pow(qcVsM,-1./2.)*qcErrorVsM;
       }      
      } // end of if(co==0) 2nd order
      else if(co==1) // 4th order
      {
       if(qcVsM<0.)
       {
        vVsM = pow(-1.*qcVsM,1./4.);
        vErrorVsM = (1./4.)*pow(-qcVsM,-3./4.)*qcErrorVsM;
       } 
      } // end of if(co==1) 4th order
      else if(co==2) // 6th order
      {
       if(qcVsM>0.)
       {
        vVsM = pow((1./4.)*qcVsM,1./6.);
        vErrorVsM = (1./6.)*pow(2.,-1./3.)*pow(qcVsM,-5./6.)*qcErrorVsM;
       }
      } // end of if(co==2) 6th order
      else if(co==3) // 8th order
      {
       if(qcVsM<0.)
       {
        vVsM = pow((-1./33.)*qcVsM,1./8.);
        vErrorVsM = (1./8.)*pow(33.,-1./8.)*pow(-qcVsM,-7./8.)*qcErrorVsM;
       }
      } // end of if(co==3) 8th order     
      // Store final results and statisticial errror for reference flow:
      refFlowVsM[f][m][co]->SetBinContent(b,vVsM);
      refFlowVsM[f][m][co]->SetBinError(b,vErrorVsM);
     } // end of for(Int_t b=1;b<=nBins;b++)
    } else
      {
       cout<<endl;
       cout<<" WARNING: "<<Form("cumulantsVsM[%i][%i][%i]",f,m,co)<<" is NULL in CalculateReferenceFlowFromCumulantsVsM() !!!!"<<endl;
       cout<<endl;
      } 
   } // end of for(Int_t co=0;co<4;co++) // cumulant order
  } // end of for(Int_t m=0;m<nMethods;m++)
 } // end of for(Int_t f=0;f<nFiles;f++)   

} // end of void CalculateReferenceFlowFromCumulantsVsM()

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
  for(Int_t m=0;m<1+(Int_t)(showAlsoGFCResults)+(Int_t)(showAlsoMCEPResults);m++)
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
       if(yMax[co] < result+error){yMax[co] = result+error;} // max value    
       // x-axis:
       if(xMax[co] < cumulantsVsM[f][m][co]->GetBinLowEdge(b+1)){xMax[co] = cumulantsVsM[f][m][co]->GetBinLowEdge(b+1);}
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
 // Get all histograms and profiles.
 
 // a) Get profiles holding results for <refMult> vs number of Reference Particles (RPs); 
 // b) Get histograms holding results for cumulants and reference flow vs multiplicity.
 
 // a) Get profiles holding results for <refMult> vs number of Reference Particles (RPs): 
 if(plotCumulantsVsReferenceMultiplicity){GetProfileRefMultVsNoOfRPs();}
  
 // b) Get histograms holding results for cumulants and reference flow vs multiplicity:
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
      // Cumulants vs multiplicity:
      cumulantsVsM[f][m][co] = dynamic_cast<TH1D*> (temp->FindObject(Form("fIntFlowQcumulantsVsM, %s",qcFlag[co].Data())));
      if(plotCumulantsVsReferenceMultiplicity && cumulantsVsM[f][m][co])
      {
       cumulantsVsM[f][m][co] = Map(cumulantsVsM[f][m][co],f);
      }    
      if(rebin && cumulantsVsM[f][m][co])
      {
       cumulantsVsM[f][m][co] = Rebin(cumulantsVsM[f][m][co],co);
      }
     } // end of for(Int_t co=0;co<4;co++)
    } // end of if(temp) 
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
      if(showAlsoGFCResults && !cumulantsVsM[f][m][co])
      {
       cout<<endl;
       cout<<Form(" WARNING: Couldn't access histogram fReferenceFlowCumulantsVsM, %s in the ",gfcFlag[co].Data())<<endl;
       cout<<"          file "<<commonOutputFiles[f]->GetName()<<" !!!!"<<endl;
       cout<<"          Did you enable calculation of GFC vs M in the analysis which produced this file?"<<endl;
       cout<<endl;  
      }
      if(plotCumulantsVsReferenceMultiplicity && cumulantsVsM[f][m][co])
      {
       cumulantsVsM[f][m][co] = Map(cumulantsVsM[f][m][co],f);
      }
      if(rebin && cumulantsVsM[f][m][co])
      {
       cumulantsVsM[f][m][co] = Rebin(cumulantsVsM[f][m][co],co);
      }       
     } // end of for(Int_t co=0;co<4;co++)
    } // end of if(temp)
   } // end of else if(!(strcmp(method[m].Data(),"GFC")))
   else if(!(strcmp(method[m].Data(),"MCEP")) && lists[f][m])
   {
    TProfile *mcepVsM = dynamic_cast<TProfile*> (lists[f][m]->FindObject("FlowPro_VsM_MCEP"));
    if(!mcepVsM)
    {
     cout<<endl;
     cout<<" WARNING: Couldn't access TProfile FlowPro_VsM_MCEP in the file "<<commonOutputFiles[f]->GetName()<<" !!!!"<<endl;
     cout<<endl;
     return;
    } 
    for(Int_t co=0;co<4;co++)
    {
     cumulantsVsM[f][m][co] = new TH1D("",Form("Monte Carlo QC{%i} #font[72]{vs} M",2*(co+1)),
                                       mcepVsM->GetNbinsX(),mcepVsM->GetBinLowEdge(1),mcepVsM->GetBinLowEdge(1+mcepVsM->GetNbinsX()));
    } // end of for(Int_t co=0;co<4;co++)
    for(Int_t b=1;b<=mcepVsM->GetNbinsX();b++)
    {
     Double_t v = mcepVsM->GetBinContent(b);
     Double_t vError = mcepVsM->GetBinError(b);
     if(TMath::Abs(v)<1.e-44 || TMath::Abs(vError)<1.e-44){continue;}
     // Monte Carlo QC{2}:     
     Double_t qc2 = pow(v,2.);
     Double_t qc2Error = 2.*TMath::Abs(v)*vError; // error propagation for f(x) = x^2
     cumulantsVsM[f][m][0]->SetBinContent(b,qc2);
     cumulantsVsM[f][m][0]->SetBinError(b,qc2Error);
     // Monte Carlo QC{4}:     
     Double_t qc4 = -pow(v,4.);
     Double_t qc4Error = 4.*TMath::Abs(pow(v,3.))*vError; // error propagation for f(x) = -x^4
     cumulantsVsM[f][m][1]->SetBinContent(b,qc4);
     cumulantsVsM[f][m][1]->SetBinError(b,qc4Error);    
     // Monte Carlo QC{6}:     
     Double_t qc6 = 4.*pow(v,6.);
     Double_t qc6Error = 24.*TMath::Abs(pow(v,5.))*vError; // error propagation for f(x) = 4x^6
     cumulantsVsM[f][m][2]->SetBinContent(b,qc6);
     cumulantsVsM[f][m][2]->SetBinError(b,qc6Error);
     // Monte Carlo QC{8}:     
     Double_t qc8 = -33.*pow(v,8.);
     Double_t qc8Error = 264.*TMath::Abs(pow(v,7.))*vError; // error propagation for f(x) = -33x^8
     cumulantsVsM[f][m][3]->SetBinContent(b,qc8);
     cumulantsVsM[f][m][3]->SetBinError(b,qc8Error);
    } // end of for(Int_t b=1;b<=mcepVsM->GetNbinsX();b++)
    for(Int_t co=0;co<4;co++)
    {
     if(plotCumulantsVsReferenceMultiplicity && cumulantsVsM[f][m][co])
     {
      cumulantsVsM[f][m][co] = Map(cumulantsVsM[f][m][co],f);
     }     
     if(rebin && cumulantsVsM[f][m][co])
     {
      cumulantsVsM[f][m][co] = Rebin(cumulantsVsM[f][m][co],co);
     }    
    } // end of for(Int_t co=0;co<4;co++)
   } // end of else if(!(strcmp(method[m].Data(),"MCEP")))
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

void GetProfileRefMultVsNoOfRPs()
{
 // Get profiles holding results for <reference multiplicity> vs number of Reference Particles (RPs).
 
 // Set here from which method's output file this profile will be accessed:
 TString methodName = "QC"; // Alternatives are GFC and MCEP
 Int_t i = -1; 
 for(Int_t m=0;m<nMethods;m++)
 {
  if(method[m] == methodName){i=m;}
 }
 if(i==-1)
 {
  cout<<endl;
  cout<<" WARNING: Unknown method name in GetProfileRefMultVsNoOfRPs() !!!!"<<endl;
  cout<<"          Try something else for TString methodName in GetProfileRefMultVsNoOfRPs()."<<endl;
  cout<<endl;exit(0);    
 }
  
 for(Int_t f=0;f<nFiles;f++)
 {
  AliFlowCommonHist *commonHist = NULL;
  if(lists[f][i])
  {
   commonHist = dynamic_cast<AliFlowCommonHist*> (lists[f][i]->FindObject(Form("AliFlowCommonHist%s",method[i].Data())));
  } else
    {
     cout<<endl;
     cout<<Form(" WARNING: lists[%i][%i] is NULL in GetProfileRefMultVsNoOfRPs() !!!!",f,i)<<endl;
     cout<<Form("          Did you use method %s in the analysis which has produced output",method[i].Data())<<endl;
     cout<<Form("          file %s?",commonOutputFiles[f]->GetName())<<endl;
    } 
  if(commonHist && commonHist->GetRefMultVsNoOfRPs())
  {
   refMultVsNoOfRPs[f] = commonHist->GetRefMultVsNoOfRPs();
  } else
    {
     cout<<endl;
     cout<<" WARNING: commonHist && commonHist->GetRefMultVsNoOfRPs() is NULL in GetProfileRefMultVsNoOfRPs() !!!!"<<endl;
     cout<<endl;
    }
 } // end of for(Int_t f=0;f<nFiles;f++)
  
 return;  
 
} // end of void GetProfileRefMultVsNoOfRPs()

// =====================================================================================

TH1D* Map(TH1D *hist, Int_t f)
{
 // Map cumulant versus # of RPs into cumulants versus <reference multiplicity>.
 
 if(!refMultVsNoOfRPs[f])
 {
  cout<<endl;
  cout<<Form(" WARNING: refMultVsNoOfRPs[%i] is NULL in Map(...) !!!!",f)<<endl;
  cout<<endl;
 }
 
 Int_t rpMinBin = refMultVsNoOfRPs[f]->FindFirstBinAbove(); // FindFirstBinAbove(Double_t threshold = 0, Int_t axis = 1) 
 Int_t rpMaxBin = refMultVsNoOfRPs[f]->FindLastBinAbove(); // FindLastBinAbove(Double_t threshold = 0, Int_t axis = 1) 
 Int_t rmMinBin = 440000;; 
 Int_t rmMaxBin = (Int_t)TMath::Floor(refMultVsNoOfRPs[f]->GetMaximum()); 
 for(Int_t rpBin=rpMinBin;rpBin<=rpMaxBin;rpBin++) // non-empty # of RPs bins
 {
  if(refMultVsNoOfRPs[f]->GetBinContent(rpBin)>0. && refMultVsNoOfRPs[f]->GetBinContent(rpBin)<rmMinBin)
  {
   rmMinBin = (Int_t)TMath::Floor(refMultVsNoOfRPs[f]->GetBinContent(rpBin));
  }
 } // end of for(Int_t rpBin=rpMinBin;rpBin<=rpMaxBin;rpBin++) // non-empty # of RPs bins 
  
 if(hist)
 {
  temp = (TH1D*) hist->Clone();
  temp->Reset();
  for(Int_t rmBin=rmMinBin;rmBin<=rmMaxBin;rmBin++) // reference multiplicity bins
  { 
   Double_t value = 0.;
   Double_t error = 0.;
   Double_t dSum1 = 0.; // sum value_i/(error_i)^2
   Double_t dSum2 = 0.; // sum 1/(error_i)^2
   for(Int_t rpBin=rpMinBin;rpBin<=rpMaxBin;rpBin++) // # of RPs bins
   {
    if((Int_t)TMath::Floor(refMultVsNoOfRPs[f]->GetBinContent(rpBin)) >= temp->GetBinLowEdge(rmBin) &&
       (Int_t)TMath::Floor(refMultVsNoOfRPs[f]->GetBinContent(rpBin)) < temp->GetBinLowEdge(rmBin+1))
    {        
     value = hist->GetBinContent(rpBin);  
     error = hist->GetBinError(rpBin);  
     if(error>0.)
     {
      dSum1+=value/(error*error);
      dSum2+=1./(error*error);
     }
    }
   } // end of for(Int_t rpBin=1;rpBin<=nBins;rpBin++) // # of RPs bins
   if(dSum2>0.)
   {
    temp->SetBinContent(rmBin,dSum1/dSum2);
    temp->SetBinError(rmBin,pow(1./dSum2,0.5));
   } 
  } // end of for(Int_t rmBin=1;rmBin<=nBins;rmBin++) // reference multiplicity bins
 } // end of if(hist)
     
 return temp;
     
} // end of Map()

// =====================================================================================

TH1D* Rebin(TH1D *hist, Int_t co)
{
 // Rebin original histograms.
 
 if(nMergedBins[co] == 0)
 {
  cout<<endl;
  cout<<Form(" WARNING: nMergedBins[%i] == 0 !!!!",co)<<endl;
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
 Int_t nBinsNew = TMath::Floor(nBinsOld/nMergedBins[co]);
 Double_t xMinNew = xMinOld;
 Double_t xMaxNew = xMinOld + nBinsNew*nMergedBins[co]*binWidthOld;
  
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
  if(b%nMergedBins[co] == 0)
  {
   if(dSum2>0.)
   {
    temp->SetBinContent(binNew,dSum1/dSum2);
    temp->SetBinError(binNew,pow(1./dSum2,0.5));
   }
   binNew++;
   dSum1 = 0.;
   dSum2 = 0.;
  } // end of if(b%nMergedBins[co] == 0)
 } // end of for(Int_t b=1;b<=nBinsOld;b++)
  
 return temp;
      
} // end of TH1D* Rebin(TH1D *hist, Int_t co)

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
 
 TH1D *styleHist = new TH1D(yAxisTitle.Data(),"",(Int_t)xMax[co],0,xMax[co]);
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
 
 if(hist)
 {
  Int_t nBins = hist->GetNbinsX();
  Double_t value = 0.;
  Double_t error = 0.;
  // Counting the non-empty bins: 
  Int_t nNonEmptyBins = 0;
  for(Int_t i=1;i<nBins+1;i++)
  {
   value = hist->GetBinContent(i);
   error = hist->GetBinError(i);
   if(TMath::Abs(value)>0.0 && error>0.0)
   {
    nNonEmptyBins++;
   }
  } // end of for(Int_t i=1;i<nBins+1;i++)  
  // Error mesh:
  TGraph *errorMesh = new TGraph(2*nNonEmptyBins+1); 
  Int_t count = 0;
  Double_t binCenter = 0.;
  for(Int_t i=1;i<nBins+1;i++)
  {
   value = hist->GetBinContent(i);
   error = hist->GetBinError(i);
   // Setting up the the mesh:
   if(TMath::Abs(value)>0.0 && error>0.0)
   {    
    binCenter = hist->GetBinCenter(i);   
    errorMesh->SetPoint(count,binCenter,value+error);
    errorMesh->SetPoint(2*nNonEmptyBins-(count+1),binCenter,value-error);
    count++;
   }
  } // end of for(Int_t i=1;i<nBins+1;i++)
  // Closing the mesh area:
  Double_t xStart = 0.;
  Double_t yStart = 0.;
  errorMesh->GetPoint(0,xStart,yStart);
  errorMesh->SetPoint(2*nNonEmptyBins,xStart,yStart);   
 } // end if(hist)
 
 return errorMesh;
 
} // end of TGraph* GetErrorMesh(TH1D *hist)

// ===========================================================================================

void GlobalSettings()
{
 // Settings which will affect all plots.
 
 gROOT->SetStyle("Plain"); // default color is white instead of gray
 gStyle->SetOptStat(0); // remove stat. box from all histos
 TGaxis::SetMaxDigits(4); // prefer exp notation for 5 and more significant digits
 
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
   if(dirFile[f][i])
   {
    TList* listTemp = dirFile[f][i]->GetListOfKeys();
    if(listTemp && listTemp->GetEntries() == 1)
    {
     listName[f][i] = listTemp->At(0)->GetName(); // to be improved - implemented better
     dirFile[f][i]->GetObject(listName[f][i].Data(),lists[f][i]);
    } else
      {
       cout<<endl;
       cout<<" WARNING: Couldn't find a list "<<listName[f][i].Data()<<" in "<<commonOutputFiles[f]->GetName()<<" !!!!"<<endl;
       cout<<"          Did you use method "<<method[i].Data()<<" in the analysis which produced this file?"<<endl;
      }
   } else 
     {
      cout<<" WARNING: Couldn't find a file "<<fileName[f][i].Data()<<".root in "<<commonOutputFiles[f]->GetName()<<" !!!!"<<endl;
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
  if(bLookInSubdirectories)
  {
   if(!(gSystem->AccessPathName(Form("%s/%s/%s",gSystem->pwd(),files[f].Data(),commonOutputFileName.Data()),kFileExists)))
   {
    commonOutputFiles[f] = TFile::Open(Form("%s/%s/%s",gSystem->pwd(),files[f].Data(),commonOutputFileName.Data()),"READ");
   } else
     { 
      cout<<endl;
      cout<<" WARNING: Couldn't find the file "<<Form("%s/%s/%s",gSystem->pwd(),files[f].Data(),commonOutputFileName.Data())<<" !!!!"<<endl;
      cout<<"          Did you specify correctly all paths in TString files[nFiles]?"<<endl;
      cout<<"          Or you are misusing 'Bool_t bLookInSubdirectories'? "<<endl;
      cout<<endl;
      exit(0);
     }
  } else // to if(bLookInSubdirectories)
    {
     if(!(gSystem->AccessPathName(Form("%s/%s%s",gSystem->pwd(),files[f].Data(),".root"),kFileExists)))
     {
      commonOutputFiles[f] = TFile::Open(Form("%s/%s%s",gSystem->pwd(),files[f].Data(),".root"),"READ");
     } else
       { 
        cout<<endl;
        cout<<" WARNING: Couldn't find the file "<<Form("%s/%s%s",gSystem->pwd(),files[f].Data(),".root")<<" !!!!"<<endl;
        cout<<"          Did you specify correctly all file names in TString files[nFiles]?"<<endl;
        cout<<endl;
        exit(0);
       }
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
    gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
    //gSystem->Load("libANALYSIS");
    gSystem->Load("libPWGflowBase");
    //cerr<<"libPWGflowBase loaded ..."<<endl;
    
  }
  
  else if (analysisMode==mLocalSource) {
    
    // In root inline compile
   
    // Constants  
    gROOT->LoadMacro("Base/AliFlowCommonConstants.cxx+");
    gROOT->LoadMacro("Base/AliFlowLYZConstants.cxx+");
    
    // Flow event
    gROOT->LoadMacro("Base/AliFlowVector.cxx+"); 
    gROOT->LoadMacro("Base/AliFlowTrackSimple.cxx+");    
    gROOT->LoadMacro("Base/AliFlowTrackSimpleCuts.cxx+");    
    gROOT->LoadMacro("Base/AliFlowEventSimple.cxx+");
   
    // Output histosgrams
    gROOT->LoadMacro("Base/AliFlowCommonHist.cxx+");
    gROOT->LoadMacro("Base/AliFlowCommonHistResults.cxx+");
    gROOT->LoadMacro("Base/AliFlowLYZHist1.cxx+");
    gROOT->LoadMacro("Base/AliFlowLYZHist2.cxx+");
    
    cout << "finished loading macros!" << endl;  
    
  }  
  
} // end of void LoadLibraries(const libModes analysisMode) 
