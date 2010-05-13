#ifndef ALICALIBVIEWERGUITIME_H
#define ALICALIBVIEWERGUITIME_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliCalibViewerGUItime.h,v */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  GUI for displaying calibration entries over run/time                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#ifndef ROOT_TGButton
#include <TGWidget.h>
#endif
#ifndef ROOT_TGFrame
#include <TGFrame.h>
#endif
#include <TString.h>
#include <TVectorT.h>

class TGCompositeFrame;
class TRootEmbeddedCanvas;
class TGTextButton;
class TGListBox;
class TGRadioButton;
class TGComboBox;
class TGGroupFrame;
class TGLabel;
class TGTabElement;
class TGWindow;

class TFile;
class TTree;
class TChain;
class TGraph;
class TObjArray;

class TMap;

class AliBaseCalibViewerGUI;

class AliCalibViewerGUItime : public TGCompositeFrame {
public:
  AliCalibViewerGUItime(const TGWindow *p, UInt_t w, UInt_t h, const Char_t* det = "TRD");
  virtual ~AliCalibViewerGUItime();

  static TObjArray* ShowGUI(TChain* chain);   // launch the time trending GUI and load a chain
  static TObjArray* ShowGUI();        // launch the empty time trending GUI

  void DrawGUI(const TGWindow */*p*/, UInt_t w, UInt_t h);
  
  void UseFile(const char* fileName, const char* treeName);
  void UseChain(TChain *chain = 0);
  void Reload(Int_t first=1);
  void AddReferenceTree(const char* treeFileName, const char* refName="R");
  
  void SetCalibViewerGUI(AliBaseCalibViewerGUI *gui) {fCalibViewerGUI=gui;}
  void SetCalibViewerGUItab(TGTabElement *tab) {fCalibViewerGUItab=tab;}
  void SetCacheDir(const char* cachedir) {fOutputCacheDir=cachedir;}
  
  const char* GetDrawString();
  const char* GetDrawOption();
  void GetCutString(TString &cutStr);
  TChain* GetChain() const {return fTree;}
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
  void SavePicture();
  void HandleLoadRunTextEntry();
  void HandleLoadRunButtons();
  
 private:
  TString fDetector;                       // TPC/TRD

  TFile*  fFile;                          //file that keeps the tree
  TChain*  fTree;                         //internal tree
  AliBaseCalibViewerGUI *fCalibViewerGUI;  //calib viewer gui used to display verbose information for one run
  TGTabElement *fCalibViewerGUItab;       //tab the calib view gui redies in
  TH1*    fCurrentHist;                   //histogram currently drawn in main canvas
  TGraph* fCurrentGraph;                  //current graph
  Int_t   fCurrentRunDetails;             //run number for wich details are currently shown
  TString fOutputCacheDir;                //output cache diretory for AliTPCCalibViewerGUI trees, created on the fly
  TString fDrawString;                    //current draw string
  Bool_t  fIsCustomDraw;                  //if custom draw string is selected
  TVectorD fRunNumbers;                   //run numbers of current selection
  TVectorD fTimeStamps;                   //timr stamps of current selection
  TVectorD fValuesX;                      //values of current selection
  TVectorD fValuesY;                      //values of current selection
  //
  Bool_t  fNoGraph;                       //Whether to create a graph
  Long64_t fGraphLimitEntries;            //limit in number of entries in the chain for producing a graph
  //
  TMap *fMapRefTrees;                      //map of reference trees for the CalibViewer
  //GUI elements
  //main canvas Top part, bottom part
  TGCompositeFrame    *fContTopBottom;      // container for all GUI elements, vertical divided
  //top left, centre, right
  TGCompositeFrame    *fContLCR;            // container for all GUI elements, horizontal divided
  //content left
  TGCompositeFrame    *fContLeft;           // container for GUI elements on left side
  TGGroupFrame        *fContDrawOpt;        // Subgroup for draw selection
  TGCheckButton       *fChkDrawOptSame;     // draw option 'same'
  TGComboBox          *fComboAddDrawOpt;    // additional draw options combo box
  TGGroupFrame        *fContDrawSel;        // Subgroup for draw selection
  TGCompositeFrame    *fContDrawSelSubRunTime; //Radio button subframe
  TGRadioButton       *fRadioXhist;         // Radio button x-variable: show only 1D distribution
  TGRadioButton       *fRadioXrun;          // Radio button x-variable: run
  TGRadioButton       *fRadioXtime;         // Radio button x-variable: time
  TGListBox           *fListVariables;      // listbox with possible variables
  TGComboBox          *fComboRunType;       // run type selection box
  TGLabel             *fLblRunType;         // run type label
  TGNumberEntry       *fNmbPar;             // parameter number
  TGLabel             *fLblPar;             // parameter name
  TGListBox           *fListCalibType;      // calibration type selection box
  TGGroupFrame        *fContCalibType;      // calibration type label
  //content centre
  TGCompositeFrame    *fContCenter;         // container for GUI elements at the center
  TRootEmbeddedCanvas *fCanvMain;           // main drawing canvas
  //content right 
  TGCompositeFrame    *fContRight;          // container for GUI elements on right side
  TGGroupFrame        *fContValues;         // container to keep data point information
  TGLabel             *fLblRunNumber;       // run number label
  TGLabel             *fLblRunTime;         // time stamp label
  TGLabel             *fLblValueX;          // value label
  TGLabel             *fLblValueY;          // value label
  TGLabel             *fLblRunNumberVal;    // run number of the data point hoovered
  TGLabel             *fLblRunTimeVal;      // time stamp of the data point hoovered
  TGLabel             *fLblValueXVal;       // value of the data point hoovered
  TGLabel             *fLblValueYVal;       // value of the data point hoovered
  TGTextButton        *fBtnDumpRuns;        // draw button
  TGTextButton        *fBtnSave;            // Save button

  TGGroupFrame        *fContLoad;           // Load file container
  TGCompositeFrame    *fContFilename;       // filename container
  TGCompositeFrame    *fContConfigFile;     // config file container
  TGCompositeFrame    *fContTreeName;       // tree name container
  TGLabel             *fLblTreeName;        // tree name label
  TGTextEntry         *fTxtFilename;        // filename text entry
  TGTextEntry         *fTxtConfigFile;      // config file text entry
  TGTextEntry         *fTxtTreeName;        // tree name text entry
  TGButton            *fBtnLoadFile;        // load file button

  //content bottom
  TGCompositeFrame    *fContCustom;         // container for custom draw command GUI elements
  TGCompositeFrame    *fContCustomCuts;     // container for custom cut options GUI elements
  TGLabel             *fLblCustomDraw;      // label for custom draw string
  TGLabel             *fLblCustomCuts;      // label for custom cuts string
  TGComboBox          *fComboCustomDraw;    // combo box custom draw string
  TGComboBox          *fComboCustomCuts;    // combo box custom cuts string
  //
  TObjArray *fTrashBox;                   //graphics objects to be deleted (histograms, graphs,...)
  
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
  void UpdateValueArrays(Bool_t withGraph);
  const char* SubstituteUnderscores(const char* in);
  void GetHistogramTitle(TString &title);
  void AdjustYRange();
private:
  AliCalibViewerGUItime(const AliCalibViewerGUItime &v);
  AliCalibViewerGUItime &operator = (const AliCalibViewerGUItime &v);         // assignment operator
  
  ClassDef(AliCalibViewerGUItime, 1)
    
};

#endif
