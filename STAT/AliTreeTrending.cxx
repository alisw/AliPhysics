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

///\ingroup STAT
///\class AliTreeTrending
///\brief AliTreeTrending class for the visualization of the QA trending/alarms
///\author Marian Ivanov
/*!
 Generalization of the original TPC QA trending visualization code
 Example usage in the $AliPhysics_SRC/PWGPP/QA
 Related JIRA task - https://alice.its.cern.ch/jira/browse/ATO-361
 // TODO -automatic resizing of pad to addapt to number of entries shown
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
#include "TLegend.h"
#include "AliTreeTrending.h"
#include <stdexcept>      // std::invalid_argument
#include "AliDrawStyle.h"

using std::cout;
using std::cerr;
using std::endl;



ClassImp(AliTreeTrending)

using std::cout;
using std::endl;

// constants
const Int_t   kMaxCanvasWidth=3000;
const Int_t   kMinCanvasWidth=1000;
const Int_t   kStepCanvasWidth=500;



AliTreeTrending::AliTreeTrending():
  TNamed(),
  fTree(NULL),
  fUserDescription(NULL),
  fLatexDescription(NULL),
  fStatusGraphM(NULL),
  fCurrentCssStyle(""),
  fReport(NULL) {
}
//_____________________________________________________________________________
AliTreeTrending::AliTreeTrending(const char *name, const char *title, const char * cssStyle):
  TNamed(name, title),
  fTree(NULL),
  fUserDescription(NULL),
  fLatexDescription(NULL),
  fStatusGraphM(NULL),
  fCurrentCssStyle(cssStyle),
  fReport(NULL) {
  //
  fReport=new TFile("report.root","recreate");
  // Register default style
  if (fCurrentCssStyle.Length()>0){
     if (AliDrawStyle::GetCssStyle(fCurrentCssStyle.Data())==NULL){
       ::Error("AliTreeTrending::AliTreeTrending", "Requested CSS Style %s not registered", fCurrentCssStyle.Data());
       fCurrentCssStyle="";
     }
  }
  if (fCurrentCssStyle.Length()==0){
    if (AliDrawStyle::GetCssStyle("AliTreeTrending") == NULL) {
      AliDrawStyle::RegisterCssStyle("AliTreeTrending", AliDrawStyle::ReadCSSFile("$AliRoot_SRC/STAT/test/AliTreeTrending.css", 0));
    }
    fCurrentCssStyle="AliTreeTrending";
  }
}

/// Change the style for drawing
/// cssStyle should be registered before  using AliDrawStyle::SetSccStyle();
/// \param cssStyle
/// \return
Bool_t AliTreeTrending::SetCssStyle(const char * cssStyle){
  if (cssStyle){
     if (AliDrawStyle::GetCssStyle(cssStyle)==NULL){
       ::Error("AliTreeTrending::SetCssStyle", "Requested CSS Style %s not registered", cssStyle);
     }else{
       fCurrentCssStyle = cssStyle;
     }
  }
}


///
/// \param latex
/// \return
TObjArray * AliTreeTrending::GetTexDescription(TLatex *latex){
  //
  // Get latex version of description
  //
  TObjArray * description= new TObjArray();
  TString sTimestamp  = TString::Format("Creation time:%s",gSystem->GetFromPipe("date").Data());
  TString sUser=TString::Format("User:%s",gSystem->GetFromPipe("echo $USER").Data());  
  TString sAlirootVer;
  TString sAliPhysicsVer;
  if (gSystem->GetFromPipe("echo $ALICE_VER") == "master" || gSystem->GetFromPipe("echo $ALICE_VER") == ""){
    sAlirootVer = "AliRoot: " + gSystem->GetFromPipe("wdir=`pwd`; cd $AliRoot_SRC/; git describe; cd $wdir;");
  }else {
    sAlirootVer = "AliRoot: " + gSystem->GetFromPipe("echo $ALIROOT_VERSION");
  }
  if (gSystem->GetFromPipe("echo $ALIPHYSICS_VER") == "master" || gSystem->GetFromPipe("echo $ALIPHYSICS_VER") == ""){
    sAliPhysicsVer = "AliPhysics: " + gSystem->GetFromPipe("wdir=`pwd`; cd $AliPhysics_SRC/; git describe; cd $wdir;");
  }else {
    sAliPhysicsVer = "AliPhysics: " + gSystem->GetFromPipe("echo $ALIPHYSICS_VERSION");
  }
  //
  description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY(),sTimestamp.Data()));
  description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY()-1.1*latex->GetTextSize(),sUser.Data()));
  description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY()-4.4*latex->GetTextSize(),sAlirootVer.Data()));
  description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY()-5.5*latex->GetTextSize(),sAliPhysicsVer.Data()));
  if (fUserDescription!=NULL){
    Int_t entries=fUserDescription->GetEntries();
    for (Int_t i=0; i<entries; i++){
      description->AddLast(latex->DrawLatexNDC(latex->GetX(), latex->GetY()-((6+i)*1.1)*latex->GetTextSize(),
					       TString::Format("%s: %s",fUserDescription->At(i)->GetName(), fUserDescription->At(i)->GetTitle()).Data()));
    }
  }
  return description;
}

/// Add user description to the report
/// \param description
void AliTreeTrending::AddUserDescription(TNamed * description){
  //
  // Add user description 
  // AliTreeTrending is owner
  // 
  if (fUserDescription==NULL) fUserDescription = new TObjArray;
  fUserDescription->AddLast(description);
}

/// ### Initialize drawing for <detector> QA
///  \param statusDescription       - contains status bar description array [variables,names,criteria]
///  \param descriptionSize         - Text size of description
///  \param cutString               - string "selection"
/// #### Example usage:
/// \code
/// Example usage of the AliTreeTrending::InitSummaryTrending
/// trendingDraw->InitSummaryTrending(statusString,0.015,"defaultCut")
/// here "defaultCut" refers to an alias set like this:   treeMC->SetAlias("defaultCut","run==TPC.Anchor.run");
/// \endcode
Bool_t  AliTreeTrending::InitSummaryTrending(TString statusDescription[3], Float_t descriptionSize, TString cutString){
  //
  // Init drawing for the <detector> QA
  // Detector specific qaConfig() has to be called before invoking this function
  //    0.) Make descriptor 
  //    1.) Make default canvas - addapt canvas width to the number of entries to draw
  //    2.) initialize status aliases (outliers etc.), status bar criteria, status lines, ...
  //    3.) compute detector  status graphs

  //   0.) Make descriptor
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
  //   1.) Make default canvas - addapt canvas width to the number of entries to draw
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
  double leftLegend  = 1 - 180./fWorkingCanvas->GetWw();
  double rightLegend = 1 - 10./fWorkingCanvas->GetWw();

  // 2.) initialize status aliases (outliers etc.), status bar criteria, status lines, ...
  TString sStatusBarVars  = statusDescription[0];
  TString sStatusBarNames = statusDescription[1];
  TString sCriteria       = statusDescription[2];
  cout << "sStatusBarVars = " << sStatusBarVars.Data() << endl;
  cout << "sCriteria      = " << sCriteria.Data() << endl;

  // 3.) compute detector status graphs
  fStatusGraphM = MakeMultiGraphStatus(fTree,"statusBar", sStatusBarVars+";tagID",sStatusBarNames,cutString.Data(),sCriteria,kTRUE);
  fStatusGraphM->GetYaxis()->SetTitle("");
  fStatusGraphM->GetHistogram()->SetTitle("");
  return kTRUE;
}

/// TODO - use CSS style
void  AliTreeTrending::SetDefaultStyle(){
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetLabelSize(0.04,"x");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

///
/// \param padRatio
/// \param bottomMargin
/// \param rightMargin
void AliTreeTrending::AppendStatusPad(Float_t padRatio, Float_t bottomMargin, Float_t rightMargin){
  //
  // Add status pad and latex description to the canvas (can be pad?)
  // status pad and latex description using data members of trendingDraw class
  //
  // 1.) Split canvas
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
  TPad* pad1 = new TPad("pad1", "pad1", 0., padRatio, 1., 1.);
  pad1->Draw();
  pad1->SetNumber(1); // so it can be called via "c1->cd(1);"
  c1->cd();
  TPad* pad2 = new TPad("pad2", "pad2", 0., 0., 1., padRatio);
  pad2->Draw();
  pad2->SetNumber(2);
  // draw original canvas into first pad
  c1->cd(1);
  c1_clone->DrawClonePad();

  pad1->SetName("Top.class(statusPad)");
  AliDrawStyle::TPadApplyStyle(fCurrentCssStyle.Data(),pad1);

//  pad1->SetBottomMargin(0.001);
//  pad1->SetRightMargin(rightMargin);
  // set up second pad
  c1->cd(2);
//  pad2->SetGrid(3);
//  pad2->SetTopMargin(0);
//  pad2->SetBottomMargin(bottomMargin); // for the long x-axis labels (run numbers)
//  pad2->SetRightMargin(rightMargin);
  pad2->SetName("Bottom.class(statusPad)");
  AliDrawStyle::TPadApplyStyle(fCurrentCssStyle.Data(),pad2); //"testStyle" to be changed to
  
  const Int_t nVars = fStatusGraphM->GetYaxis()->GetNbins();
  TGraph* grAxis = (TGraph*) fStatusGraphM->GetListOfGraphs()->At(0);
  Int_t entries = grAxis->GetN();
  fStatusGraphM->GetXaxis()->SetLabelSize(5.7*TMath::Min(TMath::Max(5./entries,0.01),0.03));
  fStatusGraphM->GetYaxis()->SetLabelSize(0.025/gPad->GetHNDC());
  fStatusGraphM->Draw("ap");
  c1->cd(1)->Draw();
}


/// \param fTree          - input tree
/// \param mgrName        - name of the output TMultiGraph // TODO - to use together CSS
/// \param expression     - expression to draw
/// \param varTitle       - title description
/// \param cutString      - selection string
/// \param sCriteria      - criteria to draw
/// \return               - return status MultiGraph
/*!
### Example usage
\code
  TString expression="dcarAP0;dcarAP1;dcarCP0;dcarCP1;dcar_posA_0;dcar_posA_1;dcar_posA_2;dcar_posC_0;dcar_posC_1;dcar_posC_2;";
  TString varTitle="dcarAP0;dcarAP1;dcarCP0;dcarCP1;dcar_posA_0;dcar_posA_1;dcar_posA_2;dcar_posC_0;dcar_posC_1;dcar_posC_2;";
  expression+="TPC.Anchor.dcarAP0;TPC.Anchor.dcarAP1;TPC.Anchor.dcarCP0;TPC.Anchor.dcarCP1;TPC.Anchor.dcar_posA_0;TPC.Anchor.dcar_posA_1;TPC.Anchor.dcar_posA_2;TPC.Anchor.dcar_posC_0;TPC.Anchor.dcar_posC_1;TPC.Anchor.dcar_posC_2;run";
  varTitle+="TPC.Anchor.dcarAP0;TPC.Anchor.dcarAP1;TPC.Anchor.dcarCP0;TPC.Anchor.dcarCP1;TPC.Anchor.dcar_posA_0;TPC.Anchor.dcar_posA_1;TPC.Anchor.dcar_posA_2;TPC.Anchor.dcar_posC_0;TPC.Anchor.dcar_posC_1;TPC.Anchor.dcar_posC_2";
  TString sCriteria("1-2*(varname_Warning):(statisticOK):(varname_Warning):(varname_Outlier):(varname_PhysAcc)"); // or to just show veto : (varname_PhysAcc&&varname_Warning)
  cutString="defaultCut";
  grTest0=AliTreeTrending::MakeMultiGraphStatus(treeMC, "DCA Status", expression, varTitle, cutString,sCriteria);
  TStatToolkit::DrawMultiGraph(grTest0,"ap");
  trendingDraw->AppendStatusPad(0.3, 0.4, 0.05);

\endcode
 */
TMultiGraph * AliTreeTrending::MakeMultiGraphStatus(TTree *fTree, TString mgrName, TString expression, TString varTitle, TCut cutString, TString sCriteria, Bool_t setAxis) {
  // 1. Make TMultiGraph for each group   
  TObjArray *oaStatusBarVars = expression.Tokenize(";:");
  TObjArray *oaStatusBarNames = varTitle.Tokenize(";:");

  if(oaStatusBarVars->GetEntries()<1) ::Error("AliTreeTrending::MakeMultiGraphStatus", "Empty variable string provided");
  if(oaStatusBarVars->GetEntries()-1!=oaStatusBarNames->GetEntries()) { 
    ::Error("AliTreeTrending::MakeMultiGraphStatus", "Variable\tNames");
    for(int i=0; i<oaStatusBarVars->GetEntries(); i++){        
      TString varVars = TString(oaStatusBarVars->At(i)!=NULL?oaStatusBarVars->At(i)->GetName():"undef.");
      TString varNames = TString(oaStatusBarNames->At(i)!=NULL?oaStatusBarNames->At(i)->GetName():"undef.");
      ::Error("AliTreeTrending::MakeMultiGraphStatus", "%s\t%s",varVars.Data(),varNames.Data());
    }
    return 0;
  }
  Int_t nVars = oaStatusBarVars->GetEntriesFast() - 1;
  TObjArray *graphArray = new TObjArray(nVars);
  for (Int_t iVar = 0; iVar < nVars; iVar++) {
    TString sVar = Form("%s:%s", oaStatusBarVars->At(iVar)->GetName(),
                        oaStatusBarVars->At(nVars)->GetName()); //e.g. -> dcar:run
    TMultiGraph *multGr = TStatToolkit::MakeStatusMultGr(fTree, sVar.Data(), cutString, sCriteria.Data(), 2*iVar);
    if (multGr) {
      graphArray->AddAt(multGr, iVar);
      ((TMultiGraph *) graphArray->At(iVar))->SetTitle(oaStatusBarNames->At(iVar)->GetName());
      for (Int_t igr = 0; igr < multGr->GetListOfGraphs()->GetEntries(); igr++) {  // Code names and codes to the object names
        TGraph *cgr = (TGraph *) multGr->GetListOfGraphs()->At(igr);
        //cgr->SetName(TString((oaStatusBarNames->At(iVar)->GetName())).Data());
      }
    } else {
      ::Error("AliTreeTrending::MakeMultiGraphStatus", "TStatToolkit::MakeStatusMultGr() returned with error -> next");
      continue;
    }
  }
  // 2.) Set Y ranges and Y labels for each graph
  Double_t *yBins = new Double_t[nVars+1];
  for(Int_t i=0; i<=nVars;i++) yBins[i] = Double_t(i);
  TMultiGraph *mgrCombined=new TMultiGraph(mgrName,"Status");
  for (Int_t i=0; i<nVars; i++) {
    TMultiGraph * mgr=((TMultiGraph*) graphArray->At(i));
    if (mgr==NULL) continue;
    for (Int_t igr=0; igr<mgr->GetListOfGraphs()->GetEntries(); igr++) {
      TGraph *cgr = (TGraph *) mgr->GetListOfGraphs()->At(igr);
      if (cgr==NULL){
        ::Error("AliTreeTrending::MakeMultiGraphStatus","Graph %d of multi-graph %s is 0", igr, mgr->GetName());
        continue;
      }
      cgr->SetTitle(""); cgr->GetYaxis()->SetTitle("");
      cgr->GetYaxis()->Set(nVars, yBins);
      cgr->GetYaxis()->SetRangeUser(0,nVars);
      for (Int_t jgr=0; jgr<nVars; jgr++)  cgr->GetYaxis()->SetBinLabel(jgr+1, graphArray->At(jgr)->GetTitle());
      mgrCombined->Add(cgr);
    }
  }
  mgrCombined->SetMinimum(0); mgrCombined->SetMaximum(nVars);
  // TStatToolkit::DrawMultiGraph(mgrCombined,"ap");
  if (setAxis) {      //TODO  - to get axis graph has to be drawn - is there other option ?
    TGraph *gr0=(TGraph *) mgrCombined->GetListOfGraphs()->At(0);
    mgrCombined->Draw("ap");
    mgrCombined->GetXaxis()->Set(gr0->GetXaxis()->GetNbins(), gr0->GetXaxis()->GetXbins()->GetArray());
    for (Int_t jgr = 0; jgr < gr0->GetXaxis()->GetNbins(); jgr++)
      mgrCombined->GetXaxis()->SetBinLabel(jgr + 1, gr0->GetXaxis()->GetBinLabel(jgr+1));
    mgrCombined->GetYaxis()->Set(nVars, yBins);
    mgrCombined->GetYaxis()->SetRangeUser(0, nVars);
    for (Int_t jgr = 0; jgr < nVars; jgr++)
      mgrCombined->GetYaxis()->SetBinLabel(jgr + 1, graphArray->At(jgr)->GetTitle());
    mgrCombined->GetYaxis()->SetTitle("");
    mgrCombined->Draw("ap");
  }
  delete oaStatusBarVars;
  delete oaStatusBarNames;
  delete graphArray;
  return mgrCombined;
}

/// #### MakePlot
/// * make a MultiGraph using array of TTree expression queries
/// * format legend
///   * crate legend box
///   * append variable.<legend> metadata into legend
/// * save resulting figure into file
/// TODO - use css style for default variables
/// \param outputDir      - directory to save figure
/// \param figureName     - name of the output figure file
/// \param LegendTitle    - legend  title
/// \param legendPos      - boundary box of legend
/// \param groupName      - name of the plot group (TO BE used in combination with the CSS )
/// \param expr           - semicolon separated list of graphs (TFormula) expression  to draw
/// \param cut            - selection
/// \param markers        - marker style specification (AliDrawStyle::GetXXXStyle())  //TODO - finish  CSS style
/// \param colors         - marker style specification (AliDrawStyle::GetXXXStyle())  //TODO - finish  CSS style
/// \param drawSparse     - kTRUE  - SparseGraph created, kFALSE - normal graph is produced
/// \param markerSize     - size of marker  //TODO - finish  CSS style
/// \param sigmaRange     - sigma range for automatic // TO use CSS style
/// \param comp           - ????
void AliTreeTrending::MakePlot(const char* outputDir, const char *figureName, const char *LegendTitle, std::vector<Double_t>& legendPos, const char *groupName, const char* expr, const char * cut, const char * markers, const char *colors, Bool_t drawSparse, Float_t markerSize, Float_t sigmaRange, Bool_t comp) {
  TMultiGraph *mGraph=0;
  fWorkingCanvas->Clear();
  TLegend *legend = new TLegend(legendPos[0],legendPos[1],legendPos[2],legendPos[3],LegendTitle);
  legend->SetBorderSize(0);
  mGraph = TStatToolkit::MakeMultGraph(fTree,groupName,expr,cut,markers,colors,drawSparse,markerSize,sigmaRange,legend,comp);
  
  if (mGraph) for(Int_t it=0; it<mGraph->GetListOfGraphs()->GetSize(); it++){
    TGraph* graph = (TGraph*) mGraph->GetListOfGraphs()->At(it);
    if (groupName!=NULL && strlen(groupName)>0){    // example  groupName=".class(multiGraphPair).{marker_style:25,21,22,23;marker_color:1,2,4,5;}"
      graph->SetName(TString::Format("graph[%d].%s",it,groupName).Data());
    }else {
      graph->SetName(TString::Format("graph[%d].class(multiGraphPair)", it).Data());
    }
	  // TODO add group name if exist
    AliDrawStyle::TGraphApplyStyle(fCurrentCssStyle.Data(),graph);
  }
  

  if(!mGraph){
    ::Error("AliTreeTrending::MakePlot","No plot returned -> dummy plot! ");
    ::Error("AliTreeTrending::MakePlot","Invalid query expression. Try tree->Draw(\"%s\",\"\t%s\")",expr, cut);
  }
  else {
    mGraph->SetName(groupName);
    if (drawSparse) TStatToolkit::RebinSparseMultiGraph(mGraph,(TGraph*)fStatusGraphM->GetListOfGraphs()->At(0));
    TStatToolkit::DrawMultiGraph(mGraph,"alp");
    AppendStatusPad(0.3, 0.4, 0.05);
    legend->SetFillStyle(0);
    legend->Draw();
  }
  if(outputDir!=0){
    fWorkingCanvas->SaveAs(TString(outputDir)+"/"+TString(figureName));
    fWorkingCanvas->Print(TString(outputDir)+"/report.pdf");
    if (fReport) {fReport->cd(); fWorkingCanvas->Write(figureName); fReport->Flush();}
  }
}

/// ### AppendBand to existing plots
/// * draw line bands
/// * optionally save figure into file
/// TO BE used for indication of dead bands
/// TODO  - add css style
/// \param outputDir     -  directory to save figure - Optional
/// \param figureName    -  output figure name - Optional
/// \param expr          -  semicolon separate list of lines (TTreeFormula)
/// \param cut           -  selection criteria
/// \param lineStyle     -  line style     //TODO use CSS style
/// \param colors        -  line colors    //TODO use CSS style
/// \param drawSparse    -  indication draw sparse
/// \param sigmaRange    -  Not used
/// \param comp          -  to remove
void AliTreeTrending::AppendBand(const char* outputDir, const char *figureName, const char* expr, const char * cut, const char * lineStyle, const char *colors, Bool_t drawSparse, Float_t sigmaRange, Bool_t comp) {
  TMultiGraph *graph=0;
  graph = TStatToolkit::MakeMultGraph(fTree,"",expr,cut,lineStyle,colors,drawSparse,0,sigmaRange,0,comp);
  if(!graph){
    ::Error("MakePlot","No plot returned -> dummy plot!");
  }
  else {
    if (drawSparse) TStatToolkit::RebinSparseMultiGraph(graph,(TGraph*)fStatusGraphM->GetListOfGraphs()->At(0));
    TStatToolkit::DrawMultiGraph(graph,"l");
  }
  if(outputDir!=0){
    fWorkingCanvas->SaveAs(TString(outputDir)+"/"+TString(figureName));
    fWorkingCanvas->Print(TString(outputDir)+"/report.pdf");
    if (fReport) {fReport->cd();fWorkingCanvas->Write(figureName);}
  }
  static Int_t counter=0;
  AliSysInfo::AddStamp(expr,2,counter++);
}

/// #### MakeStatusPlot
/// * make status graphs
/// * optionally save figure
/// \param outputDir      - directory to save file
/// \param figureName     - figure name
/// \param expression     - semicolon separated list of the status variables (variables to be evaluated)
/// \param varTitle       - titles to show on y axis
/// \param cutString      - selection string - cut
/// \param sCriteria      - criteria to draw ()
/// \return               - return status TMultiGraph
/*!
#### Example usage:
* assuming trendingDraw was properly initialized)
\code
  TString expression="dcarAP0;dcarAP1;dcarCP0;dcarCP1;dcar_posA_0;dcar_posA_1;dcar_posA_2;dcar_posC_0;dcar_posC_1;dcar_posC_2;";
  TString varTitle="dcarAP0;dcarAP1;dcarCP0;dcarCP1;dcar_posA_0;dcar_posA_1;dcar_posA_2;dcar_posC_0;dcar_posC_1;dcar_posC_2;";
  expression+="TPC.Anchor.dcarAP0;TPC.Anchor.dcarAP1;TPC.Anchor.dcarCP0;TPC.Anchor.dcarCP1;TPC.Anchor.dcar_posA_0;TPC.Anchor.dcar_posA_1;TPC.Anchor.dcar_posA_2;TPC.Anchor.dcar_posC_0;TPC.Anchor.dcar_posC_1;TPC.Anchor.dcar_posC_2;run";
  varTitle+="TPC.Anchor.dcarAP0;TPC.Anchor.dcarAP1;TPC.Anchor.dcarCP0;TPC.Anchor.dcarCP1;TPC.Anchor.dcar_posA_0;TPC.Anchor.dcar_posA_1;TPC.Anchor.dcar_posA_2;TPC.Anchor.dcar_posC_0;TPC.Anchor.dcar_posC_1;TPC.Anchor.dcar_posC_2";
  cutString="defaultCut";
  trendingDraw->MakeStatusPlot("./", "dcaStatus.png", expression, varTitle, cutString,sCriteria);
\endcode
*/
void AliTreeTrending::MakeStatusPlot(const char *outputDir, const char *figureName,  TString expression, TString varTitle,  TCut cutString, TString sCriteria, TString friendName) {
  fWorkingCanvas->Clear();
  TMultiGraph *multiGraph = NULL;
  if (friendName.Length() > 0 && fTree->GetFriend(friendName)){
    multiGraph =  AliTreeTrending::MakeMultiGraphStatus(fTree->GetFriend(friendName), "", expression, varTitle, cutString, sCriteria,kTRUE);
  }else{
    multiGraph =  AliTreeTrending::MakeMultiGraphStatus(fTree, "", expression, varTitle, cutString, sCriteria,kTRUE);
  }

  if(multiGraph==0) {
    ::Error("AliTreeTrending::MakeStatusPlot", "MakeMultiGraphStatus returned NULL pointer");
    return;
  }
  TStatToolkit::RebinSparseMultiGraph(multiGraph,(TGraph*)fStatusGraphM->GetListOfGraphs()->At(0)); // rebin to the fStatusBar
  TStatToolkit::DrawMultiGraph(multiGraph, "ap");
  //multiGraph->Draw("ap");
  AppendStatusPad(0.3, 0.4, 0.05);
  if (outputDir != NULL) {
    fWorkingCanvas->SaveAs(TString(outputDir) + "/" + TString(figureName));
    fWorkingCanvas->Print(TString(outputDir) + "/report.pdf");
    if (fReport) {fReport->cd();fWorkingCanvas->Write(figureName);}
  }
  static Int_t counter=0;
  AliSysInfo::AddStamp(expression.Data(),3,counter++);
}



///
/// \param figureName       - figure name
/// \param refVariable      - reference  variable
/// \param bandNamePrefix   - band prefix
/// \param selection        - selection
/// \param groupName        - group name (used to specify class or properties)
void AliTreeTrending::AppendDefaultBands(const char *outputDir, const char *figureName,  const char * refVariable,const char * bandNamePrefix, const char * selection, const char* groupName /*other drawing variable  TString style*/ ){
  TMultiGraph *mGraph=0;
  const char* aType[3]={"_WarningBand","_OutlierBand","_PhysAccBand"}; //yellow 400 ,red 632 ,green 416/
  TString expr;
  for(Int_t itype=0; itype<3; itype++){
    expr = refVariable+TString("+")+bandNamePrefix+TString(aType[itype])+TString(";")+refVariable+TString("-")+bandNamePrefix+TString(aType[itype])+TString(":run");
    mGraph = TStatToolkit::MakeMultGraph(fTree,"",expr,selection, "figTemplateTRDPair", "figTemplateTRDPair",kTRUE,0,6,0,kTRUE);
    if(!mGraph){
      ::Error("MakePlot","No plot returned -> dummy plot!");
    }
    else {
      if (kTRUE) TStatToolkit::RebinSparseMultiGraph(mGraph,(TGraph*)fStatusGraphM->GetListOfGraphs()->At(0));
      for(Int_t it=0; it<mGraph->GetListOfGraphs()->GetSize(); it++){
        TGraph* band = (TGraph*) mGraph->GetListOfGraphs()->At(it);
        if (groupName && strlen(groupName)>1){    // example  groupName=".class(multiGraphPair).{marker_style:25,21,22,23;marker_color:1,2,4,5;}"
          band->SetName(TString::Format("graph[%d]%s",itype,groupName).Data());
        }else {
          band->SetName(TString::Format("graph[%d].class(deadBand)", itype).Data());
        }
        // TODO add group name if exist
        AliDrawStyle::TGraphApplyStyle(fCurrentCssStyle,band); //"testStyle" to be changed to fCurrenCSSStyle
      }
      TStatToolkit::DrawMultiGraph(mGraph,"l");
    }
  }
  if(outputDir!=0){
    fWorkingCanvas->SaveAs(TString(outputDir)+"/"+TString(figureName));
    fWorkingCanvas->Print(TString(outputDir)+"/report.pdf");
    if (fReport) {fReport->cd();fWorkingCanvas->Write(figureName);}
  }
  static Int_t counter=0;
  AliSysInfo::AddStamp(expr,2,counter++);
}

/// Recursive function to decompose the status string into sub-contribution
/// \param tree            - input tree
/// \param currentString   - current status string
/// \param statusVar       - resulting status variable
/// \param statusTitle     - resulting status variable
/// \param suffix          - suffix to parse status string "_Warning"
/// \param counter         - counter for debug purposes
/*!
\code
 TString currentString="dcar_Warning";
 TString statusVar="", statusTitle="", bitMaskAlias="";
 TPRegexp suffix("_Warning$");
 Int_t counter=0;
 AliTreeTrending::DecomposeStatusAlias(treeMC, currentString,statusVar,statusTitle,suffix,counter, bitMaskAlias);
 statusVar+="run";
 trendingDraw->MakeStatusPlot("./", "dcarStatusMC.png", statusVar, statusTitle, "defaultCut",sCriteria);

 trendingDraw->MakeStatusPlot("./", "dcarStatusAnchor.png", statusVar, statusTitle, "defaultCut",sCriteria,"TPC.Anchor");
 \endcode
*/
void AliTreeTrending::DecomposeStatusAlias(TTree* tree, TString currentString, TString &statusVar, TString &statusTitle, TPRegexp &suffix, Int_t &counter, TString &bitMaskAlias){
  //
  if (tree==NULL) throw std::invalid_argument("invalid tree argument");
  counter++;
  TString toAdd=currentString;
  suffix.Substitute(toAdd,"");
  statusVar+=toAdd;
  statusVar+=";";
  if (counter>1) bitMaskAlias+="+";
  bitMaskAlias+=toAdd;
  bitMaskAlias+="_#";
  bitMaskAlias+=TString::Format("*(1<<%d)",counter);
  TString title=toAdd;
  if (TStatToolkit::GetMetadata(tree,(toAdd+".Title").Data())) title=TStatToolkit::GetMetadata(tree,toAdd+".Title")->GetTitle();
  statusTitle+=title;
  statusTitle+=";";
  TString content=tree->GetAlias(currentString.Data());
  TObjArray  * array = content.Tokenize("|&");
  printf("%d\t%s\t%s\n",counter,toAdd.Data(), content.Data());
  for (Int_t i=0; i<array->GetEntries(); i++){
    TString cString=array->At(i)->GetName();
    cString.ReplaceAll("(","");
    cString.ReplaceAll(")","");
    if (suffix.Match(cString,"")>0) {
      DecomposeStatusAlias(tree, cString, statusVar, statusTitle,suffix,counter,bitMaskAlias);
    }
  }
}


/// Print to the output string object names fulfilling criteria sRegExp
/// \param array            - input TCollection
/// \param regExp           - regular expression to select
/// \param separator        -
/// \return                 - string
TString AliTreeTrending::ArrayNameToString(TCollection *array, TString sRegExp, TString separator){
  TString output="";
  TPRegexp regExp(sRegExp);
  TIter next(array);
  TObject *object;
  while ( ( object =next())) {
    if (regExp.Match(object->GetName())) {
      output += object->GetName();
      output += separator;
    }
  }
  return output;
}


/// AddHtmlLink - helper function to build html page for JSROOT
/// \param pFile   - pointer to the output html text file
/// \param title   - Title of the section
/// \param prefix  - prefix Query
/// \param items   - comma separated list of figures to include
/*!
    TString prefix="http://localhost:90/data/jsroot/index.htm?file=http://localhost:90/data/alice-tpc-notes/JIRA/ATO-83/test/report.root"
    items=AliTreeTrending::ArrayNameToString(fReport->GetListOfKeys(),".*tatus.*",",");
    AliTreeTrending::AddJSROOTHtmlLink(pFile,"Status", prefix, Form("items=[%s]",items.Data()));
*/
void AliTreeTrending::AddJSROOTHtmlLink(FILE * pFile, TString title,  TString prefix, TString items){
  fprintf(pFile, "<b>%s</b>",title.Data());
  fprintf(pFile, "<a href=\"%s&%s&layout=tabs\">[all,</a>",prefix.Data(),items.Data());
  fprintf(pFile, "<a href=\"%s&%s&layout=collapsible&nobrowser\">col,</a>",prefix.Data(),items.Data());
  fprintf(pFile, "<a href=\"%s&%s&layout=tabs&nobrowser\">tabs,</a>",prefix.Data(),items.Data());
  fprintf(pFile, "<a href=\"%s&%s&layout=flex&nobrowser\">flex]</a>",prefix.Data(),items.Data());
  fprintf(pFile, "<br>");
}


