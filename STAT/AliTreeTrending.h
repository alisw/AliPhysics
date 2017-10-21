#ifndef ALITREETRENDING_H
#define ALITREETRENDING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */
///\ingroup STAT
///\class AliTreeTrending
///\brief AliTreeTrending class for the visualization of the QA trending/alarms
///\author Marian Ivanov
/*!
 Generalization of the original TPC QA trending visualization code
 Example usage in the $AliPhysics_SRC/PWGPP/QA
 Related JIRA task - https://alice.its.cern.ch/jira/browse/ATO-361
*/

class TTree;
class TObjArray;
class TNamed;
class TStyle;
class TMultiGraph;
class TLegend;

class AliTreeTrending : public TNamed {
public:
  AliTreeTrending();
  AliTreeTrending(const char *name, const char *title);
  ~AliTreeTrending(){;}
  void SetTree(TTree * tree) {fTree=tree;};
  void AddUserDescription(TNamed *description);  //
  Bool_t  InitSummaryTrending(TString statusDescription[3], Float_t descriptionSize=0.015, TString cutString="");
  void SetDefaultStyle();  // own style to be used - not yet
  void AppendStatusPad(Float_t padRatio, Float_t bottomMargin, Float_t rightMargin);
  //
  void MakePlot(const char* outputDir, const char *figureName, const char *LegendTitle, std::vector<Double_t>& legendPos, const char *groupName, const char* expr, const char * cut, const char * markers, const char *colors, Bool_t drawSparse, Float_t markerSize, Float_t sigmaRange, Bool_t comp);
  void AppendBand(const char* outputDir, const char *figureName, const char* expr, const char * cut, const char * lineStyle, const char *colors, Bool_t drawSparse, Float_t sigmaRange, Bool_t comp) ;
  void MakeStatusPlot(const char* outputDir, const char *figureName, TString expression, TString varTitle, TCut cutString, TString sCriteria);
  static TMultiGraph * MakeMultiGraphStatus(TTree *fTree, TString mgrName, TString expression, TString varTitle, TCut cutString, TString sCriteria, Bool_t setAxis=kFALSE);
  // TODO void MakeHtml(char *htmlName, char *varList)
public:
  TObjArray *  GetTexDescription(TLatex *latex);  /// Currently only standard variables
  TTree     *  fTree;              /// working tree with friends
  TObjArray *  fUserDescription;   /// user defined description
  TObjArray *  fLatexDescription;  /// description of process  (time/aliroot/root/user)
  TCanvas   *  fWorkingCanvas;     /// default canvas
  TMultiGraph *fStatusGraphM;      /// status graph
  TStyle    *  fDrawStyle;         /// TreeTrending owner of drawing style - currently gStyle used
  ClassDef(AliTreeTrending,2)
};


#endif

