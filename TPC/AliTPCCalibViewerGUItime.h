#ifndef AliTPCCalibViewerGUItime_H
#define AliTPCCalibViewerGUItime_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCCalibViewerGUI.h,v */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  GUI for displaying calibration entires over time                         //
//  Calibration Trees are created using the macro TPC/CalibMacros/CalibEnv.C //
//  used for the calibration monitor                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#ifndef ROOT_TGButton
#include "TGWidget.h"
#endif
#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif

// #include <TGButton.h>
// #include <TGListBox.h>
// #include <TGComboBox.h>
// #include <TGNumberEntry.h>
// #include <TRootEmbeddedCanvas.h>
// #include <TGSplitter.h>
// #include <TGButtonGroup.h>
// #include <TGLabel.h>
// #include <TGTab.h>
#include <TString.h>
#include <TVectorT.h>

class TGCompositeFrame;
class TRootEmbeddedCanvas;
class TGTextButton;
class TGListBox;
class TGRadioButton;
class TGGroupFrame;
class TGLabel;
class TGComboBox;
class TGTabElement;

class TFile;
class TTree;
class TGraph;

class AliTPCCalibViewerGUI;
class AliTPCConfigParser;



class AliTPCCalibViewerGUItime : public TGCompositeFrame {
public:
  AliTPCCalibViewerGUItime(const TGWindow *p, UInt_t w, UInt_t h);
  virtual ~AliTPCCalibViewerGUItime();
  
  static TObjArray* ShowGUI(const char* fileName = 0);             // initialize and show GUI standalone

  void DrawGUI(const TGWindow */*p*/, UInt_t w, UInt_t h);
  
  void UseFile(const char* fileName);
  void Reload(Int_t first=1);

  
  void SetCalibViewerGUI(AliTPCCalibViewerGUI *gui) {fCalibViewerGUI=gui;}
  void SetCalibViewerGUItab(TGTabElement *tab) {fCalibViewerGUItab=tab;}
  void SetConfigFile(const char* file) {fConfigFile=file;}
  const char* GetDrawString();
  const char* GetCutString();
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
  
private:
  TFile*  fFile;                          //file that keeps the tree
  TTree*  fTree;                          //internal tree
  AliTPCCalibViewerGUI *fCalibViewerGUI;  //calib viewer gui used to display verbose information for one run
  TGTabElement *fCalibViewerGUItab;       //tab the calib view gui redies in
  TH1*    fCurrentHist;                   //histogram currently drawn in main canvas
  TGraph* fCurrentGraph;                  //current graph
  Int_t   fCurrentRunDetails;             //run number for wich details are currently shown
  TString fOutputCacheDir;                //output cache diretory for AliTPCCalibViewerGUI trees, created on the fly
  TString fDrawString;                    //current draw string
  TString fConfigFile;                    //configuration file keeping active branches and branch descriptions
  AliTPCConfigParser *fConfigParser;      //configuration parser
  Bool_t  fIsCustomDraw;                  //if custom draw string is selected
  TVectorD fRunNumbers;                   //run numbers of current selection
  TVectorD fTimeStamps;                   //timr stamps of current selection
  TVectorD fValuesX;                      //values of current selection
  TVectorD fValuesY;                      //values of current selection
  
  //GUI elements
  //main canvas Top part, bottom part
  TGCompositeFrame    *fContTopBottom;      // container for all GUI elements, vertical divided
  //top left, centre, right
  TGCompositeFrame    *fContLCR;            // container for all GUI elements, horizontal divided
  //content left
  TGCompositeFrame    *fContLeft;           // container for GUI elements on left side
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
  //content bottom
  TGCompositeFrame    *fContCustom;         // container for custom draw command GUI elements
  TGCompositeFrame    *fContCustomCuts;     // container for custom cut options GUI elements
  TGLabel             *fLblCustomDraw;      // label for custom draw string
  TGLabel             *fLblCustomCuts;      // label for custom cuts string
  TGComboBox          *fComboCustomDraw;    // combo box custom draw string
  TGComboBox          *fComboCustomCuts;    // combo box custom cuts string
  
  enum { kRadioXhist=10, kRadioXrun=11, kRadioXtime=12 };
  enum { kBranchOnOff=0, kBranchTitle=1, kCalibType=2, kParamNames=3 };
  
  void UpdateParLimits();
  void UpdateParName();
  void SetGuiTree(Int_t run);
  void FillRunTypes();
  void FillCalibTypes();
  void SetInitialValues();
  const char* SubstituteUnderscores(const char* in);
  
  AliTPCCalibViewerGUItime(const AliTPCCalibViewerGUItime &v);
  AliTPCCalibViewerGUItime &operator = (const AliTPCCalibViewerGUItime &v);         // assignment operator
  
  ClassDef(AliTPCCalibViewerGUItime, 0)
    
};

#endif
