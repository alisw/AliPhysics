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

void RebinCFContainer(const char *infile="AnalysisResults.root",Int_t rebinVar=0, Int_t myEff=3, const char * name="Nch"){

  gSystem->SetIncludePath("-I. -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/RAW -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWGPP -g");

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
  AliCFContainer *data = 0;
  d = (TDirectoryFile*)file->Get(Form("PWGPP_CFSingleTrack%s",name));
  if(!d) {
    cout<<" no directory "<<endl;
    return;
  }
  data = (AliCFContainer*) (d->Get("container"));
  if(!data){
    cout <<" no container"<<endl;
  }

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
  if(rebinVar==0) {
    nLimits = 15;
    newLimits = new Double_t[nLimits+1];
    newLimits[0]=0;
    newLimits[1]=0.5;
    newLimits[2]=1;
    newLimits[3]=1.5;
    newLimits[4]=2;
    newLimits[5]=2.5;
    newLimits[6]=3;
    newLimits[7]=4;
    newLimits[8]=5;
    newLimits[9]=6;
    newLimits[10]=7;
    newLimits[11]=8;
    newLimits[12]=10;
    newLimits[13]=12;
    newLimits[14]=14;
    newLimits[15]=16;
  } else if (rebinVar==1) {
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
    heff->Write("hpteffCFLimAccGenAcc");
    h2num->Write("hptLimAcc");
    h2->Write("hptGenAcc");
  }
  else if( myEff == 2 ) {
    fout = new TFile(Form("efficiency_STE_RecPPR_over_GenAcc_rebinned_%d_%s.root",rebinVar,name),"RECREATE");
    heff->Write("hpteffCFRecPPRGenAcc");
    h2num->Write("hptRecPPR");
    h2->Write("hptGenAcc");
  }
  else if( myEff == 3 ) {
    fout = new TFile(Form("efficiency_STE_RecPID_over_GenAcc_rebinned_%d_%s.root",rebinVar,name),"RECREATE");
    heff->Write("hpteffCFRecPIDGenAcc");
    h2num->Write("hptRecPID");
    h2->Write("hptGenAcc");
  }
  else if( myEff == 4 ) {
    fout = new TFile(Form("efficiency_STE_Rec_over_GenAcc_rebinned_%d_%s.root",rebinVar,name),"RECREATE");
    heff->Write("hpteffCFRecGenAcc");
    h2num->Write("hptRec");
    h2->Write("hptGenAcc");
  }
  else if( myEff == 5 ) {
    fout = new TFile(Form("efficiency_STE_RecAcc_over_GenAcc_rebinned_%d_%s.root",rebinVar,name),"RECREATE");
    heff->Write("hpteffCFRecAccGenAcc");
    h2num->Write("hptRecAcc");
    h2->Write("hptGenAcc");
  }
  else if( myEff == 6 ) {
    fout = new TFile(Form("efficiency_STE_RecPPRKine_over_GenAcc_rebinned_%d_%s.root",rebinVar,name),"RECREATE");
    heff->Write("hpteffCFRecPPRKineGenAcc");
    h2num->Write("hptRecPPRKine");
    h2->Write("hptGenAcc");
  }
  fout->Close();

}
