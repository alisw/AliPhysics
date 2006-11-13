#ifndef ALIANALYSISGUI_H
#define ALIANALYSISGUI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliAnalysisGUIFrame
//   AliAnalysisGUI class that describes the overall analysis GUI
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliAnalysisGUI                                //
//                                                                      //
//              Implementation fo the analysis GUI.                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TGDockableFrame.h>
#include <TGrid.h>

class TGToolBar;
class TGTab;
class TGCanvas;
class TGStatusBar;
class TGPicture;
class TGIcon;
class TGListTreeItem;
class TGPopupMenu;
class TGMenuBar;
class TGHorizontal3DLine;

//GUI
class AliFileListFrame;
class AliLoginFrame;
class AliAlienBrowser;
class AliTagFrame;
class AliPackageFrame;

#include "AliTagAnalysisFrame.h"
#include "AliSelectorFrame.h"

enum ECommandIdentifiers {
   kMFILELOGIN,
   kMFILEOPEN,
   kMFILESAVEAS,
   kMFILETAG,
   kMFILEEXIT
};

//___________________________________________________________________________
class AliAnalysisGUI : public TGMainFrame {
  
 public:
  AliAnalysisGUI(const TGWindow *p, UInt_t w, UInt_t h);
  ~AliAnalysisGUI();
  
  //___________________________________________________________________________
  void   CloseWindow();
  Bool_t LogIn(const char * server, const char* username="");
  Bool_t IsConnected() const {return fIsConnected;}
  
  // slot
  void   HandleMenu(Int_t id);
  void   HandleToolBar(Int_t id);
  void   OnDoubleClick(TGListTreeItem* item, Int_t btn);
  
  //___________________________________________________________________________
 private:
  AliAnalysisGUI(const AliAnalysisGUI&);
  AliAnalysisGUI& operator= (const AliAnalysisGUI&);
  
  // private methods
  void AddMenuBar();
  void AddToolBar();
  void AddStatusBar();
  void ChangeRightLogo(const char *name);
  
  TGHorizontalFrame   *fHFrame1; //horizontal frame
  TGVerticalFrame     *fVFrame1, *fVFrame2; //verticla frames
  TGDockableFrame     *fMenuDock; //main menu
  TGPopupMenu         *fMenuFile; //main popup menu
  TGMenuBar           *fMenuBar; //menu bar
  TGToolBar           *fToolBar; //the button tool bar
  TGTab               *fTab; //tab objects
  TGLayoutHints       *fMenuBarLayout, *fMenuBarItemLayout; //layout
  TGHorizontal3DLine  *fH3DLine; //3d line
  TGCanvas            *fCanvas2; //canvas
  TGStatusBar         *fStatusBar; //status bar
  
  AliAlienBrowser     *fAliEnBrowser; //the catalog browser
  AliFileListFrame    *fFileListFrame; //the file list tab
  AliLoginFrame       *fLogInFrame; //the login frame
  AliTagFrame         *fTagFrame; //the tag frame
  AliTagAnalysisFrame *fTagAnalysisFrame; //the event tag tab
  AliPackageFrame     *fPackageFrame; //the package tab
  AliSelectorFrame    *fSelectorFrame; //the selector tab
  
  const TGPicture     *fIcon; //picture
  TGPicture           *fRightIconPicture; //picture
  TGIcon              *fRightIcon; //icon
  Bool_t               fIsConnected; //alien connected
  TGrid               *fAlien; //api pointer
   
  ClassDef(AliAnalysisGUI, 0); // AliAnalysisGUI
};

#endif
