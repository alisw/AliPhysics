#ifndef AliTPCCalibViewerGUItime_H
#define AliTPCCalibViewerGUItime_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliTPCCalibViewerGUItime
///
/// \brief GUI for displaying calibration entires over time
///  Calibration Trees are created using the macro TPC/CalibMacros/CalibEnv.C
///  used for the calibration monitor


#ifndef ROOT_TGButton
#include "TGWidget.h"
#endif
#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif

#include <TGComboBox.h>
#include <TString.h>
#include <TVectorT.h>

class TGCompositeFrame;
class TRootEmbeddedCanvas;
class TGTextButton;
class TGListBox;
class TGRadioButton;
class TGGroupFrame;
class TGLabel;
class TGTabElement;
class TGTextEntry;

class TFile;
class TTree;
class TChain;
class TGraph;
class TObjArray;

class TMap;

class AliTPCCalibViewerGUI;
class AliTPCConfigParser;



class AliTPCCalibViewerGUItime : public TGCompositeFrame {
public:
  AliTPCCalibViewerGUItime(const TGWindow *p, UInt_t w, UInt_t h);
  virtual ~AliTPCCalibViewerGUItime();
  
  static TObjArray* ShowGUI(const char* fileName = 0, const char* treeName="dcs");             // initialize and show GUI standalone
  static TObjArray* ShowGUI(TChain* chain);                                         // initialize and show GUI standalone
  
  void DrawGUI(const TGWindow */*p*/, UInt_t w, UInt_t h);
  
  void UseFile(const char* fileName, const char* treeName);
  void UseChain(TChain *const chain);
  void UseConfigFile(const char* file="");
  void Reload(Int_t first=1);
  void AddReferenceTree(const char* treeFileName, const char* refName="R");

  
  void SetCalibViewerGUI(AliTPCCalibViewerGUI *const gui) {fCalibViewerGUI=gui;}
  void SetCalibViewerGUItab(TGTabElement *const tab) {fCalibViewerGUItab=tab;}
  void SetCacheDir(const char* cachedir) {fOutputCacheDir=cachedir;}
  void SetConfigFileName(const char* file) {fConfigFile=file;}
  
  const TString GetDrawString();
  const TString GetDrawOptionString();
  const char* GetCustomDrawString() const {return fComboCustomDraw->GetTextEntry()?fComboCustomDraw->GetTextEntry()->GetText():"";}
  void GetCutString(TString &cutStr);
  TChain* GetChain() const {return fTree;}
  //
  TGTextEntry* GetDrawEntry() {return fComboCustomDraw->GetTextEntry();}
  TGTextEntry* GetCutsEntry() {return fComboCustomCuts->GetTextEntry();}
  TGTextEntry* GetDrawOptEntry() {return fComboAddDrawOpt->GetTextEntry();}
  //Slots
  void DoDraw();
  void DoDumpRuns();
  void DoCustomDraw();
  void DoCustomCutsDraw();
  void DoParLimitChange();
  void DoNewSelection();                    // decides whether to redraw if user makes another selection
  void DoChangeSelectionList() {Reload(0);}
  void HandleButtonsDrawSel(Int_t id = -1);              
  void MouseMove(Int_t event, Int_t x, Int_t y, TObject */*selected*/);
  void DoNewSelectionAliases();
  void DoAddAlias();
  void DoDelAlias();
  void UpdateAliasList();
  TCanvas * GetCanvas(){ return fCanvMain->GetCanvas();}
 private:
  TFile*  fFile;                          ///< file that keeps the tree
  TChain*  fTree;                         ///< internal tree
  AliTPCCalibViewerGUI *fCalibViewerGUI;  ///< calib viewer gui used to display verbose information for one run
  TGTabElement *fCalibViewerGUItab;       ///< tab the calib view gui redies in
  TH1*    fCurrentHist;                   ///< histogram currently drawn in main canvas
  TGraph* fCurrentGraph;                  ///< current graph
  Int_t   fCurrentRunDetails;             ///< run number for wich details are currently shown
  TString fOutputCacheDir;                ///< output cache diretory for AliTPCCalibViewerGUI trees, created on the fly
  TString fDrawString;                    ///< current draw string
  TString fConfigFile;                    ///< configuration file keeping active branches and branch descriptions
  AliTPCConfigParser *fConfigParser;      ///< configuration parser
  Bool_t  fIsCustomDraw;                  ///< if custom draw string is selected
  TVectorD fRunNumbers;                   ///< run numbers of current selection
  TVectorD fTimeStamps;                   ///< timr stamps of current selection
  TVectorD fValuesX;                      ///< values of current selection
  TVectorD fValuesY;                      ///< values of current selection
  //
  Bool_t  fNoGraph;                       ///< Whether to create a graph
  Long64_t fGraphLimitEntries;            ///< limit in number of entries in the chain for producing a graph
  //
  TMap *fMapRefTrees;                      ///< map of reference trees for the CalibViewer
  //GUI elements
  //main canvas Top part, bottom part
  TGCompositeFrame    *fContTopBottom;      ///< container for all GUI elements, vertical divided
  //top left, centre, right
  TGCompositeFrame    *fContLCR;            ///< container for all GUI elements, horizontal divided
  //content left
  TGCompositeFrame    *fContLeft;           ///< container for GUI elements on left side
  TGGroupFrame        *fContDrawOpt;        ///< Subgroup for draw selection
  TGCheckButton       *fChkDrawOptSame;     ///< draw option 'same'
  TGCheckButton       *fChkDrawOptSparse;   ///< draw option 'sparse'
  TGComboBox          *fComboAddDrawOpt;    ///< additional draw options combo box
  TGGroupFrame        *fContDrawSel;        ///< Subgroup for draw selection
  TGCompositeFrame    *fContDrawSelSubRunTime; ///< Radio button subframe
  TGRadioButton       *fRadioXhist;         ///< Radio button x-variable: show only 1D distribution
  TGRadioButton       *fRadioXrun;          ///< Radio button x-variable: run
  TGRadioButton       *fRadioXtime;         ///< Radio button x-variable: time
  TGListBox           *fListVariables;      ///< listbox with possible variables
  TGComboBox          *fComboRunType;       ///< run type selection box
  TGLabel             *fLblRunType;         ///< run type label
  TGNumberEntry       *fNmbPar;             ///< parameter number
  TGLabel             *fLblPar;             ///< parameter name
  TGListBox           *fListCalibType;      ///< calibration type selection box
  TGGroupFrame        *fContCalibType;      ///< calibration type label
  //content centre
  TGCompositeFrame    *fContCenter;         ///< container for GUI elements at the center
  TRootEmbeddedCanvas *fCanvMain;           ///< main drawing canvas
  //content right 
  TGCompositeFrame    *fContRight;          ///< container for GUI elements on right side
  TGGroupFrame        *fContValues;         ///< container to keep data point information
  TGLabel             *fLblRunNumber;       ///< run number label
  TGLabel             *fLblRunTime;         ///< time stamp label
  TGLabel             *fLblValueX;          ///< value label
  TGLabel             *fLblValueY;          ///< value label
  TGLabel             *fLblRunNumberVal;    ///< run number of the data point hoovered
  TGLabel             *fLblRunTimeVal;      ///< time stamp of the data point hoovered
  TGLabel             *fLblValueXVal;       ///< value of the data point hoovered
  TGLabel             *fLblValueYVal;       ///< value of the data point hoovered
  TGTextButton        *fBtnDumpRuns;        ///< draw button
  TGGroupFrame        *fContAliases;         ///< container to keep data point information
  TGListBox           *fListAliases;        ///< list of aliases
  //content bottom
  TGCompositeFrame    *fContCustom;         ///< container for custom draw command GUI elements
  TGCompositeFrame    *fContCustomCuts;     ///< container for custom cut options GUI elements
  TGLabel             *fLblCustomDraw;      ///< label for custom draw string
  TGLabel             *fLblCustomCuts;      ///< label for custom cuts string
  TGComboBox          *fComboCustomDraw;    ///< combo box custom draw string
  TGComboBox          *fComboCustomCuts;    ///< combo box custom cuts string
  //
  TObjArray *fTrashBox;                   ///< graphics objects to be deleted (histograms, graphs,...)
  
  enum { kRadioXhist=10, kRadioXrun=11, kRadioXtime=12 };
  enum { kBranchOnOff=0, kBranchTitle=1, kCalibType=2, kParamNames=3 };
  
  void UpdateParLimits();
  void UpdateParName();
  void SetGuiTree(Int_t run);
  void FillRunTypes();
  void FillCalibTypes();
  void SetInitialValues();
  void CheckDrawGraph();
  Bool_t CheckChain();
  void UpdateValueArrays(Bool_t withGraph, const Double_t *xArr);
  void SubstituteUnderscores(TString &s);
  void GetHistogramTitle(TString &title);
  void AdjustYRange();
private:
  AliTPCCalibViewerGUItime(const AliTPCCalibViewerGUItime &v);
  AliTPCCalibViewerGUItime &operator = (const AliTPCCalibViewerGUItime &v);         // assignment operator
  
  ClassDef(AliTPCCalibViewerGUItime, 0)
    
};

////////////////////////////////////////////////////////////////////////
//
//   GUI Alias frame
//
////////////////////////////////////////////////////////////////////////

class AliTPCCalibViewerGUItimeAddAliasFrame : public TObject {
public:
  AliTPCCalibViewerGUItimeAddAliasFrame(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
             UInt_t options, AliTPCCalibViewerGUItime *gui, TString strAlias="");
  virtual ~AliTPCCalibViewerGUItimeAddAliasFrame();
  
   // slots
  void DoOK();
  void DoCancel();
  

private:
  TGTransientFrame    *fMain;           ///< Main frame
  TGTextEntry         *fTxt1, *fTxt2;   ///< text input

  AliTPCCalibViewerGUItime *fGUI;       ///< pointer to mother process

  AliTPCCalibViewerGUItimeAddAliasFrame(const AliTPCCalibViewerGUItimeAddAliasFrame &r);
  AliTPCCalibViewerGUItimeAddAliasFrame &operator = (const AliTPCCalibViewerGUItimeAddAliasFrame &r);

  /// \cond CLASSIMP
  ClassDef(AliTPCCalibViewerGUItimeAddAliasFrame,0)
  /// \endcond
};

#endif
