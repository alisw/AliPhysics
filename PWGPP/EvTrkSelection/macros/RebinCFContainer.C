#include "TDirectoryFile.h"
#include "THnBase.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TFile.h"
#include "AliCFGridSparse.h"
#include "AliCFContainer.h"
#include "AliCFEffGrid.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "Riostream.h"

//
// Macro to project (and rebin) the output of CFSingleTrackEfficiencyTask
//
// Parameters:
//   input file name
//   variable to project and rebin
//   steps to project (1:GenAcc/LimAcc, 2: RecCuts/GenAcc, 3: RecPID/GenAcc, 4: Rec/GenAcc, 5: RecAcc/GenAcc)
//   container directory suffix
//
const Double_t etamin = -0.9;
const Double_t etamax =  0.9;
const Double_t ptmin = 0.0;
const Double_t ptmax = 24.0;
const Double_t phimin = -2*TMath::Pi();
const Double_t phimax = 2*TMath::Pi();
const Double_t thetamin = 0;
const Double_t thetamax = TMath::Pi();
const Double_t zvtxmin = -10.0;
const Double_t zvtxmax =  10.0;
const Float_t multmin = -1;
const Float_t multmax = 103;//102; (adding +-1 to the limit to include possible over(under)fow events in the projection)
const Float_t centmin = -1;
const Float_t centmax = 101;//100; (adding +-1 to the limit to include possible over(under)fow events in the projection)

void RebinCFContainer(const char *infile="AnalysisResults.root",Int_t rebinVar=0, Int_t myEff=3, const char * name="Nch"){

  gSystem->SetIncludePath("-I. -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/RAW -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_PHYSICS/PWGPP -g");

  // gROOT->LoadMacro("AliSingleTrackEffCuts.cxx++g");
  // gROOT->LoadMacro("AliCFSingleTrackEfficiencyTask.cxx++g");

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetOptTitle(0);
  
  
  TFile* file = TFile::Open(infile,"read");
  
  TDirectoryFile *d = 0;
  AliCFContainer *dataIni = 0;
  d = (TDirectoryFile*)file->Get(Form("PWGPP_CFSingleTrack"));
  if(!d) {
    cout<<" no directory "<<endl;
    return;
  }
  dataIni = (AliCFContainer*) (d->Get(Form("container%s",name)));
  if(!dataIni){
    cout <<" no container"<<endl;
  }
  //
  // Do an slice of the container to make sure to remove outliers
  //
  const UInt_t ipt = 0;
  const UInt_t ieta  = 1;
  const UInt_t iphi = 2;
  const UInt_t itheta = 3;
  const UInt_t izvtx = 4;
  const UInt_t imult = 5;
  const UInt_t icent = 6;
  Int_t nvars = 7;
  Int_t* ivarSlice = new Int_t[nvars];
  ivarSlice[0]= ipt;    ivarSlice[1] = ieta;   ivarSlice[2] = iphi;  ivarSlice[3] = itheta;
  ivarSlice[4]= izvtx;  ivarSlice[5] = imult;  ivarSlice[6] = icent;
  Double_t *mins = new Double_t[nvars];
  Double_t *maxs = new Double_t[nvars];
  mins[ipt] = ptmin;        maxs[ipt] = ptmax;
  mins[ieta] = etamin;      maxs[ieta] = etamax;
  mins[iphi] = phimin;      maxs[iphi] = phimax;
  mins[itheta] = thetamin;  maxs[itheta] = thetamax;
  mins[izvtx] = zvtxmin;    maxs[izvtx] = zvtxmax;
  mins[imult] = multmin;    maxs[imult] = multmax;
  mins[icent] = centmin;    maxs[icent] = centmax;
  AliCFContainer *data = (AliCFContainer*)dataIni->MakeSlice(nvars,ivarSlice,mins,maxs);
  cout<< "  ... slice done"<<endl;
  cout<< "  the new container has "<< data->GetNStep() << " steps"<<endl;
  

  //
  // *********** NUMERATOR
  //
  AliCFGridSparse* gridSparse = 0;		
  if( myEff == 1 )  gridSparse = (AliCFGridSparse*)data->GetGrid(1); // GenAcc
  else if( myEff == 2 )  gridSparse = (AliCFGridSparse*)data->GetGrid(6); // Rec (acc + cuts) draw reco properties
  else if( myEff == 3 )  gridSparse = (AliCFGridSparse*)data->GetGrid(7); // RecPID
  else if( myEff == 4 )  gridSparse = (AliCFGridSparse*)data->GetGrid(3); // Rec (no cuts)
  else if( myEff == 5 )  gridSparse = (AliCFGridSparse*)data->GetGrid(4); // RecAcc
  else if( myEff == 6 )  gridSparse = (AliCFGridSparse*)data->GetGrid(5); // Rec (acc + cuts) draw kine properties

  THnSparse* numData = (THnSparse*)gridSparse->GetGrid();
  
  // method 2: defining a new THnSparse changing only the axis of interest
  Printf("Method 2 ");
  THnSparse* newnumData = (THnSparse*)numData->Clone("numNew");
  newnumData->Reset();

  Int_t nLimits=0;
  Double_t* newLimits =0;
  TString varname="";
  if(rebinVar==0) {
    varname="pt";
    nLimits = 17;
    newLimits = new Double_t[nLimits+1];
    newLimits[0]=0;
    newLimits[1]=0.25;
    newLimits[2]=0.5;
    newLimits[3]=0.75;
    newLimits[4]=1;
    newLimits[5]=1.5;
    newLimits[6]=2;
    newLimits[7]=2.5;
    newLimits[8]=3;
    newLimits[9]=4;
    newLimits[10]=5;
    newLimits[11]=6;
    newLimits[12]=7;
    newLimits[13]=8;
    newLimits[14]=10;
    newLimits[15]=12;
    newLimits[16]=14;
    newLimits[17]=16;
  } else if (rebinVar==1) {
    varname="eta";
    nLimits = 9;
    newLimits = new Double_t[nLimits+1];
    newLimits[0]=-0.9;
    newLimits[1]=-0.7;
    newLimits[2]=-0.5;
    newLimits[3]=-0.3;
    newLimits[4]=-0.1;
    newLimits[5]=0.1;
    newLimits[6]=0.3;
    newLimits[7]=0.5;
    newLimits[8]=0.7;
    newLimits[9]=0.9;
  } else if (rebinVar==2) {
    varname="phi";
    nLimits = 9;
    newLimits = new Double_t[nLimits+1];
    for(Int_t i=0; i<=nLimits; i++) newLimits[i]=(Double_t)phimin + (phimax-phimin)/nLimits*(Double_t)i ;
    
  }
  const Int_t nnewLimits = nLimits;
  
  
  TAxis* axis = (TAxis*)newnumData->GetAxis(rebinVar);
  axis->Set(nnewLimits,newLimits);
  newnumData->SetBinEdges(rebinVar,newLimits);
  
  newnumData->RebinnedAdd(numData, 1);
  
  // checking the bin contents
  TH1D* h1 = (TH1D*)numData->Projection(rebinVar);
  Float_t sum = 0;
  Float_t sumnew = 0;
  Printf("Original THnSparse");
  for (Int_t ibin = 1; ibin<=h1->GetNbinsX(); ibin++){
    Printf("ibin = %d, low edge = %f, content = %f",ibin,h1->GetBinLowEdge(ibin),h1->GetBinContent(ibin));
    sum+=h1->GetBinContent(ibin);
  }

  Printf("THnSparse changing only one axis");
  TH1D* h2num = (TH1D*)newnumData->Projection(rebinVar);
  for (Int_t ibin = 1; ibin<=h2num->GetNbinsX(); ibin++){
    Printf("ibin = %d, low edge = %f, content = %f",ibin,h2num->GetBinLowEdge(ibin),h2num->GetBinContent(ibin));
    sumnew+=h2num->GetBinContent(ibin);
  }
  
  Printf("sum = %f, sumnew = %f",sum, sumnew);

  //  
  // *********** DENOMINATOR RECPID
  //
  Printf("DENOMINATOR");
  
  AliCFGridSparse* gridSparse2 = 0;
  if( myEff == 1 )  gridSparse2 = (AliCFGridSparse*)data->GetGrid(0); // LimAcc
  else gridSparse2 = (AliCFGridSparse*)data->GetGrid(1);              // GenAcc

  THnSparse* denData = (THnSparse*)gridSparse2->GetGrid();
  
  // method 2: defining a new THnSparse changing only the axis of interest
  Printf("Method 2 ");
  THnSparse* newdenData = (THnSparse*)denData->Clone("denNew");
  newdenData->Reset();
  
  TAxis* axis2 = (TAxis*)newdenData->GetAxis(rebinVar);
  axis2->Set(nnewLimits,newLimits);
  newdenData->SetBinEdges(rebinVar,newLimits);
  
  newdenData->RebinnedAdd(denData, 1);
  
  // checking the bin contents
  TH1D* h1d = (TH1D*)denData->Projection(rebinVar);
  sum = 0;
  sumnew = 0;
  Printf("Original THnSparse");
  for (Int_t ibin = 1; ibin<=h1d->GetNbinsX(); ibin++){
    Printf("ibin = %d, low edge = %f, content = %f",ibin,h1d->GetBinLowEdge(ibin),h1d->GetBinContent(ibin));
    sum+=h1d->GetBinContent(ibin);
  }
  
  Printf("THnSparse changing only one axis");
  TH1D* h2 = (TH1D*)newdenData->Projection(rebinVar);
  for (Int_t ibin = 1; ibin<=h2->GetNbinsX(); ibin++){
    Printf("ibin = %d, low edge = %f, content = %f",ibin,h2->GetBinLowEdge(ibin),h2->GetBinContent(ibin));
    sumnew+=h2->GetBinContent(ibin);
  }
  
  Printf("sum = %f, sumnew = %f",sum, sumnew);
  
  TH1D* heff = (TH1D*)h2num->Clone("heff");
  heff->Divide(h2num, h2,1,1,"B");
  //  heff->GetXaxis()->SetRangeUser(0,23.5);
  TFile* fout = NULL;
  if( myEff == 1 ) {
    fout = new TFile(Form("efficiency_STE_GenAcc_over_LimAcc_rebinned_%d_%s.root",rebinVar,name),"RECREATE");
    heff->Write(Form("h_%s_effCFLimAccGenAcc",varname.Data()));
    h2num->Write(Form("h_%s_LimAcc",varname.Data()));
    h2->Write(Form("h_%s_GenAcc",varname.Data()));
  }
  else if( myEff == 2 ) {
    fout = new TFile(Form("efficiency_STE_RecPPR_over_GenAcc_rebinned_%d_%s.root",rebinVar,name),"RECREATE");
    heff->Write(Form("h_%s_effCFRecPPRGenAcc",varname.Data()));
    h2num->Write(Form("h_%s_RecPPR",varname.Data()));
    h2->Write(Form("h_%s_GenAcc",varname.Data()));
  }
  else if( myEff == 3 ) {
    fout = new TFile(Form("efficiency_STE_RecPID_over_GenAcc_rebinned_%d_%s.root",rebinVar,name),"RECREATE");
    heff->Write(Form("h_%s_effCFRecPIDGenAcc",varname.Data()));
    h2num->Write(Form("h_%s_RecPID",varname.Data()));
    h2->Write(Form("h_%s_GenAcc",varname.Data()));
  }
  else if( myEff == 4 ) {
    fout = new TFile(Form("efficiency_STE_Rec_over_GenAcc_rebinned_%d_%s.root",rebinVar,name),"RECREATE");
    heff->Write(Form("h_%s_effCFRecGenAcc",varname.Data()));
    h2num->Write(Form("h_%s_Rec",varname.Data()));
    h2->Write(Form("h_%s_GenAcc",varname.Data()));
  }
  else if( myEff == 5 ) {
    fout = new TFile(Form("efficiency_STE_RecAcc_over_GenAcc_rebinned_%d_%s.root",rebinVar,name),"RECREATE");
    heff->Write(Form("h_%s_effCFRecAccGenAcc",varname.Data()));
    h2num->Write(Form("h_%s_RecAcc",varname.Data()));
    h2->Write(Form("h_%s_GenAcc",varname.Data()));
  }
  else if( myEff == 6 ) {
    fout = new TFile(Form("efficiency_STE_RecPPRKine_over_GenAcc_rebinned_%d_%s.root",rebinVar,name),"RECREATE");
    heff->Write(Form("h_%s_effCFRecPPRKineGenAcc",varname.Data()));
    h2num->Write(Form("h_%s_RecPPRKine",varname.Data()));
    h2->Write(Form("h_%s_GenAcc",varname.Data()));
  }
  fout->Close();

}
