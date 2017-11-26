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
  AliTreeTrending(const char *name, const char *title, const char *cssStyle="");
  ~AliTreeTrending(){;}
  void SetTree(TTree * tree) {fTree=tree;};
  void AddUserDescription(TNamed *description);  //
  Bool_t  InitSummaryTrending(TString statusDescription[3], Float_t descriptionSize=0.015, TString cutString="");
  void SetDefaultStyle();  // own style to be used - not yet
  Bool_t SetCssStyle(const char *cssStyle="");
  void AppendStatusPad(Float_t padRatio, Float_t bottomMargin, Float_t rightMargin);
  //
  void MakePlot(const char* outputDir, const char *figureName, const char *LegendTitle, std::vector<Double_t>& legendPos, const char *groupName, const char* expr, const char * cut, const char * markers, const char *colors, Bool_t drawSparse, Float_t markerSize, Float_t sigmaRange, Bool_t comp);
  void AppendBand(const char* outputDir, const char *figureName, const char* expr, const char * cut, const char * lineStyle, const char *colors, Bool_t drawSparse, Float_t sigmaRange, Bool_t comp) ;
  void AppendDefaultBands(const char *outputDir, const char *figureName,const char * refVariable,const char * bandNamePrefix, const char * selection, const char* groupName);
  void MakeStatusPlot(const char* outputDir, const char *figureName, TString expression, TString varTitle, TCut cutString, TString sCriteria, TString friendName="");
  static TMultiGraph * MakeMultiGraphStatus(TTree *fTree, TString mgrName, TString expression, TString varTitle, TCut cutString, TString sCriteria, Bool_t setAxis=kFALSE);
  static void  DecomposeStatusAlias(TTree* tree, TString currentString, TString &statusVar, TString &statusTitle, TPRegexp &suffix, Int_t &counter, TString &maskAlias);
  // TODO void MakeHtml(char *htmlName, char *varList)
  // JSROOT export related helper functions
  static TString ArrayNameToString(TCollection *array, TString sRegExp, TString separator);
  static void AddJSROOTHtmlLink(FILE * pFile, TString title,  TString prefix, TString items);
public:
  TObjArray *  GetTexDescription(TLatex *latex);  /// Currently only standard variables
  TTree     *  fTree;              /// working tree with friends
  TObjArray *  fUserDescription;   /// user defined description
  TObjArray *  fLatexDescription;  /// description of process  (time/aliroot/root/user)
  TCanvas   *  fWorkingCanvas;     /// default canvas
  TMultiGraph *fStatusGraphM;      /// status graph
  TStyle    *  fDrawStyle;         /// TreeTrending owner of drawing style - currently gStyle used  // TODO make  OBSOLETE
  TString      fCurrentCssStyle;   /// Name of the  CSS style to be used - if empty string current CSS style used
  TFile     *  fReport;            /// root report file
  ClassDef(AliTreeTrending,2)
};


#endif

