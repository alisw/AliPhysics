/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


/*
 .L $ALICE_ROOT/../src/STAT/AliTreeTrending.cxx+

*/

#include "TStatToolkit.h"
#include "Riostream.h"
#include <iostream>
#include "TSystem.h"
#include "TNamed.h"
#include "TFile.h"
#include "TTree.h"
#include "TPRegexp.h"
#include "TFriendElement.h"
#include "AliExternalInfo.h"
#include "TTreeFormula.h"
#include "TTreeFormulaManager.h"
#include "AliTreeTrending.h"
#include "TEntryList.h"
#include "AliLog.h"
#include "TStyle.h"
#include "TROOT.h"
#include "AliTreeTrending.h"


ClassImp(AliTreeTrending)

// constants
const Int_t   kMaxCanvasWidth=3000;
const Int_t   kMinCanvasWidth=1000;
const Int_t   kStepCanvasWidth=500;



AliTreeTrending::AliTreeTrending():
  TNamed(),
  fTree(NULL),
  fUserDescription(NULL),
  fLatexDescription(NULL),
  fStatusGraph(NULL)    
{
}
//_____________________________________________________________________________
AliTreeTrending::AliTreeTrending(const char *name, const char *title):
  TNamed(name, title),
  fTree(NULL),
  fUserDescription(NULL),
  fLatexDescription(NULL),
  fStatusGraph(NULL) 
{
  //
  //
}


TObjArray * AliTreeTrending::GetTexDescription(TLatex *latex){
  //
  // Get latex version of description
  //
  TObjArray * description= new TObjArray();
  TString sTimestamp  = TString::Format("Creation time:%s",gSystem->GetFromPipe("date").Data());
  TString sUser=TString::Format("User:%s",gSystem->GetFromPipe("echo $USER").Data());  
  TString sAlirootVer;
  TString sAliphysicsVer;
  if (gSystem->GetFromPipe("echo $ALICE_VER") == "master" || gSystem->GetFromPipe("echo $ALICE_VER") == ""){
    sAlirootVer = "AliRoot: " + gSystem->GetFromPipe("wdir=`pwd`; cd $ALICE_ROOT/../src; git describe; cd $wdir;");
  }else {
    sAlirootVer = "AliRoot: " + gSystem->GetFromPipe("echo $ALIROOT_VERSION");
  }
  if (gSystem->GetFromPipe("echo $ALIPHYSICS_VER") == "master" || gSystem->GetFromPipe("echo $ALIPHYSICS_VER") == ""){
    sAliphysicsVer = "AliPhysics: " + gSystem->GetFromPipe("wdir=`pwd`; cd $ALICE_PHYSICS/../src; git describe; cd $wdir;");
  }else {
    sAliphysicsVer = "AliPhysics: " + gSystem->GetFromPipe("echo $ALIPHYSICS_VERSION");
  }
  //
  description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY(),sTimestamp.Data()));
  description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY()-1.1*latex->GetTextSize(),sUser.Data()));
  description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY()-4.4*latex->GetTextSize(),sAlirootVer.Data()));
  description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY()-5.5*latex->GetTextSize(),sAliphysicsVer.Data()));
  if (fUserDescription!=NULL){
    Int_t entries=fUserDescription->GetEntries();
    for (Int_t i=0; i<entries; i++){
      description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY()-((6+i)*1.1)*latex->GetTextSize(),
					       TString::Format("%s: %s",fUserDescription->At(i)->GetName(), fUserDescription->At(i)->GetTitle()).Data()));
    }
  }
  return description;
}

void AliTreeTrending::AddUserDescription(TNamed * description){
  //
  // Add user description 
  // AliTreeTrending is owner
  // 
  if (fUserDescription==NULL) fUserDescription = new TObjArray;
  fUserDescription->AddLast(description);
}


Bool_t  AliTreeTrending::InitSummaryTrending(TString statusDescription[3], Float_t descriptionSize){
  //
  // Init drawing for the <detector> QA
  // Detector specific qaConfig() has to be called before invoking this function
  //    0.) Make descriptor 
  //    1.) Make default canvas - addopt canvas width to the number of entries to draw
  //    3.) compute detector  status graphs

  //
  //   0.) Make descriptor 
  //
  if (fTree==NULL) {
    AliError("Input tree not defined");
    return 0;
  }
  if (fTree->GetAlias("tagID")==NULL && fTree->GetBranch("tagID")==NULL) {
    AliError("tagID is not defined");
    return 0;
  }
  TLatex *latex= new TLatex;
  latex->SetX(0.11);
  latex->SetY(0.8);
  latex->SetTextSize(descriptionSize);
  fLatexDescription = GetTexDescription(latex);
  //
  //   1.) Make default canvas - addopt canvas width to the number of entries to draw
  //
  TGraphErrors *gr = (TGraphErrors*) TStatToolkit::MakeGraphSparse(fTree,"tagID:tagID","");
  Int_t numberOfTags = gr->GetN();
  cout<<"number of graph entries: "<<numberOfTags<<endl;
  double SpaceForLegend = 0.;
  Int_t canvas_width  = SpaceForLegend + (numberOfTags+5)*30;
  Int_t canvas_height = 600; 
  if ( canvas_width>kMaxCanvasWidth)   canvas_width=kMaxCanvasWidth;
  if ( canvas_width<kMinCanvasWidth)   canvas_width=kMinCanvasWidth;

  fWorkingCanvas = new TCanvas("fWorkingCanvas","fWorkingCanvas",canvas_width,canvas_height);
  fWorkingCanvas->SetGrid(3);
  fWorkingCanvas->cd();
  gPad->SetTicks(1,2);
  // fWorkingCanvas->SetRightMargin(SpaceForLegend/canvas_width);
  double leftlegend  = 1 - 180./fWorkingCanvas->GetWw();
  double rightlegend = 1 - 10./fWorkingCanvas->GetWw();
  //
  // 2.) process config file qaConfig.C to initialize status aliases (outliers etc.), status bar criteria, status lines, ...
  //
  TString sStatusbarVars  = statusDescription[0];
  TString sStatusbarNames = statusDescription[1];
  TString sCriteria       = statusDescription[2];
  cout << "sStatusbarVars = " << sStatusbarVars.Data() << endl;
  cout << "sCriteria      = " << sCriteria.Data() << endl;
  //
  // 3.) compute detector status graphs
  //
  TObjArray* oaStatusbarVars = sStatusbarVars.Tokenize(";");
  TObjArray* oaStatusbarNames = sStatusbarNames.Tokenize(";");
  fStatusGraph = new TObjArray();
  int igr=0;
  for (Int_t vari=oaStatusbarVars->GetEntriesFast()-1; vari>=0; vari--){ // invert the order of the status graphs
    TString sVar = Form("%s:tagID", oaStatusbarVars->At(vari)->GetName()); //e.g. -> dcar:run
    fStatusGraph->Add( TStatToolkit::MakeStatusMultGr(fTree, sVar.Data(),  "", sCriteria.Data(), igr) );
    TString sYtitle = oaStatusbarNames->At(vari)->GetName(); // set better name for y axis of statuspad
    ((TMultiGraph*) fStatusGraph->At(igr))->SetTitle(sYtitle.Data());
    igr++;
  }
  return kTRUE;
}

void  AliTreeTrending::SetDefaultStyle(){
  //
  // ??? Who is owner of drawing style ??? 
  //
  //
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetLabelSize(0.04,"x");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}



void AliTreeTrending::AppendStatusPad(Float_t padratio, Float_t bottomMargin, Float_t rightMargin){
  //
  // Add status pad and latex description to the canvas (can be pad?)
  // status pad and latex description using data mbers of trendingDraw class
  //
  // 1.) Split canvase
  // 2.) Draw status bar in lower part of canvas
  //
  if (fLatexDescription==NULL){
    AliError("LatexDescription not initialized. Run AliTreeTrending::InitSummaryTrending");
    return;
  }
  TCanvas *c1 = fWorkingCanvas;
  c1->cd();
  fLatexDescription->DrawClone();
  //
  // 1.) Split canvas
  //
  TCanvas* c1_clone = (TCanvas*) c1->Clone("c1_clone");  
  c1->Clear();
  // produce new pads
  c1->cd();
  TPad* pad1 = new TPad("pad1", "pad1", 0., padratio, 1., 1.); 
  pad1->Draw();
  pad1->SetNumber(1); // so it can be called via "c1->cd(1);"
  c1->cd();
  TPad* pad2 = new TPad("pad2", "pad2", 0., 0., 1., padratio);
  pad2->Draw();
  pad2->SetNumber(2);
  // draw original canvas into first pad
  c1->cd(1);
  c1_clone->DrawClonePad();
  pad1->SetBottomMargin(0.001);
  pad1->SetRightMargin(rightMargin);
  // set up second pad
  c1->cd(2);
  pad2->SetGrid(3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(bottomMargin); // for the long x-axis labels (runnumbers)
  pad2->SetRightMargin(rightMargin);  
  //
  //
  const Int_t nvars = fStatusGraph->GetEntriesFast();
  TGraph* grAxis = (TGraph*) ((TMultiGraph*) fStatusGraph->At(0))->GetListOfGraphs()->At(0);
  grAxis->SetMaximum(0.5*nvars+0.5);
  grAxis->SetMinimum(0);
  grAxis->GetYaxis()->SetLabelSize(0);
  grAxis->GetYaxis()->SetTitle("");
  grAxis->SetTitle("");
  Int_t entries = grAxis->GetN();
  grAxis->GetXaxis()->SetLabelSize(5.7*TMath::Min(TMath::Max(5./entries,0.01),0.03));
  grAxis->GetXaxis()->LabelsOption("v");
  grAxis->Draw("ap");
  //
  // draw multigraphs & names of status variables on the y axis
  for (Int_t i=0; i<nvars; i++){
    ((TMultiGraph*) fStatusGraph->At(i))->Draw("p");
    TLatex* ylabel = new TLatex(-0.1, 0.5*i+0.5, ((TMultiGraph*) fStatusGraph->At(i))->GetTitle());
    ylabel->SetTextAlign(32); //hor:right & vert:centered
    ylabel->SetTextSize(0.025/gPad->GetHNDC());
    ylabel->Draw();
  } 
  c1->cd(1)->Draw();
}

