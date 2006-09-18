#ifndef ALIMONITORCLIENT_H
#define ALIMONITORCLIENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TGFrame.h>
#include <RQ_OBJECT.h>
#include "AliMonitorDialog.h"

class TGListTreeItem;
class TGLayoutHints;
class TGPopupMenu;
class TGMenuBar;
class TGHorizontal3DLine;
class TGToolBar;
class TGNumberEntry;
class TGPicture;
class TGButton;
class TGCanvas;
class TGListTree;
class TGVSplitter;
class TGHSplitter;
class TRootEmbeddedCanvas;
class TGTextView;
class TGLabel;
class TGTextEntry;
class TSocket;
class TFileHandler;
class TFolder;
class TTimer;
class AliMonitorHisto;


class AliMonitorClient : public TGMainFrame {

RQ_OBJECT("AliMonitorClient")

public:
  AliMonitorClient();
  virtual ~AliMonitorClient();

  void               CloseWindow();

  void               OnNewData();

  void               OnMenuActivated(Int_t id);

  void               OnEventNumberChanged();
  void               OnEventButtonPressed();
  void               OnSumButtonPressed();
  void               OnRunButtonPressed();
  void               OnLoopButtonClicked();
  void               OnPreviousButtonClicked();
  void               OnNextButtonClicked();
  void               OnCopyButtonClicked();
  void               OnSaveButtonClicked();
  void               OnPrintButtonClicked();

  void               OnTreeClicked(TGListTreeItem* item, Int_t btn);
  void               OnTreeReturnPressed(TGListTreeItem* item);

  void               OnLoopTimer();

private:
  AliMonitorClient(const AliMonitorClient& client);   // Not implemented
  AliMonitorClient& operator = (const AliMonitorClient& client); // Not implemented

  TFolder*           CreateTopFolder() const;
  AliMonitorHisto*   GetHisto(const char* folderName, const char* histoName);
  TGListTreeItem*    GetItem(TGListTreeItem* base, const char* folderName, 
			     const char* histoName, Bool_t create);

  Bool_t             ConnectToServer();
  void               DisconnectFromServer();
  Bool_t             OpenFile();

  void               ViewToolBar(Bool_t visible);
  void               ViewTree(Bool_t visible);
  void               ViewDescription(Bool_t visible);
  void               ViewReference(Bool_t visible);
  void               ViewStatistics(Bool_t visible);

  Bool_t             AddFavorite();
  Bool_t             DeleteFavorite();
  Bool_t             LoadFavorites(Bool_t dialog = kTRUE);
  Bool_t             SaveFavorites();
  Bool_t             SaveFavoritesAs();

  Bool_t             LoadReference(Bool_t dialog = kTRUE);
  void               SetReference();
  Bool_t             TakeCurrentReference();
  void               TakeReferenceHisto(const char* folderName,
					AliMonitorHisto* histo);
  void               TakeReferenceFolder(TGListTreeItem* item);
  Bool_t             SaveReference();
  Bool_t             SaveReferenceAs();

  void               LoadSettings();
  void               SaveSettings();

  void               StopLoop();
  Bool_t             GetBaseItem();
  Bool_t             GoToNextItem();
  Bool_t             GoToPreviousItem();
  void               UpdateItem(Bool_t highlight);

  Bool_t             CheckForNewData();
  void               ClearItems(TGListTreeItem* base) const;
  void               CleanUpTree(TGListTreeItem* base);
  void               UpdateTree();
  void               UpdateFavoritesTree();
  void               UpdateComparisonTree();
  void               UpdateDescription();
  void               UpdateHisto();
  void               UpdateAll();

  TGLayoutHints*     fMenuBarLayout;           // layout of the menu bar
  TGLayoutHints*     fMenuBarItemLayout;       // layout of the menu items
  TGLayoutHints*     fMenuBarHelpLayout;       // layout of the help menu
  TGPopupMenu*       fMenuFile;                // the file menu
  TGPopupMenu*       fMenuView;                // the view menu
  TGPopupMenu*       fMenuFavorites;           // the favorites menu
  TGPopupMenu*       fMenuReference;           // the reference menu
  TGPopupMenu*       fMenuOptions;             // the options menu
  TGPopupMenu*       fMenuHelp;                // the help menu
  TGMenuBar*         fMenuBar;                 // the menu bar

  TGLayoutHints*     fToolBarLayout;           // layout of the tool bar
  TGHorizontal3DLine* fToolBarSep;             // separation of the tool bar
  TGToolBar*         fToolBar;                 // the tool bar
  TGLayoutHints*     fEventNumberLayout;       // layout of the event number
  TGNumberEntry*     fEventNumber;             // the event number
  TGButton*          fEventButton;             // the button for one event
  TGButton*          fSumButton;               // the button for the sum of events
  TGButton*          fRunButton;               // the button for a run
  TGButton*          fLoopButton;              // the botton for the loop
  const TGPicture*   fLoopOnPicture;           // the picture for running loop
  const TGPicture*   fLoopOffPicture;          // the picture for stoped loop
  TGButton*          fPreviousButton;          // the button for previous histo
  TGButton*          fNextButton;              // the button for next histo
  TGButton*          fCopyButton;              // the button for copy histo
  TGButton*          fSaveButton;              // the button for save histo
  TGButton*          fPrintButton;             // the button for print histo

  TGLayoutHints*     fBottomLayout;            // layout at bottom
  TGLayoutHints*     fLeftLayout;              // layout at left
  TGLayoutHints*     fExpandLayout;            // expanded layout

  TGVerticalFrame*   fVerticalFrame;           // frame for tree/histo and description
  TGHorizontalFrame* fHorizontalFrame;         // frame for tree and histo

  TGCompositeFrame*  fTreeFrame;               // frame for tree
  TGCanvas*          fTreeCanvas;              // canvas for tree
  TGListTree*        fTree;                    // tree with histos
  const TGPicture*   fHistoPicture;            // picture for histo item
  TGListTreeItem*    fAllItem;                 // top item for all histos
  TGListTreeItem*    fFavoritesItem;           // top item for favorites
  TGListTreeItem*    fComparisonItem;          // top item for comparison

  TGVSplitter*       fTreeSplitter;            // splitter for tree and histo

  TGCompositeFrame*  fDrawFrame;               // frame for histo
  TRootEmbeddedCanvas* fDrawCanvas;            // canvas for histo

  TGHSplitter*       fDescriptionSplitter;     // splitter tree/histo and description

  TGCompositeFrame*  fDescriptionFrame;        // frame for description
  TGTextView*        fDescription;             // description text

  TString            fServerName;              // name of the monitor server
  TSocket*           fSocket;                  // socket to the monitor server
  TFileHandler*      fSocketHandler;           // handler for fSocket

  TFolder*           fFolder;                  // folder with histos

  TGListTreeItem*    fCurrentItem;             // current tree item
  TGListTreeItem*    fBaseItem;                // base item of current item
  TTimer*            fLoopTimer;               // timer for loop over histos
  Int_t              fLoopInterval;            // loop interval

  TString            fFavoritesFileName;       // file name of favorites

  TString            fReferenceFileName;       // file name with reference histos
  TFolder*           fReference;               // folder with reference histos

  TString            fPrintCommand;            // print command

  static const char* fgSettingsFileName;       // file name of settings


  class AliMonitorStringDlg : public AliMonitorDialog {

  public:
    AliMonitorStringDlg(TString& string, TGFrame* main, const char* title,
			const char* label);
    virtual ~AliMonitorStringDlg();

    virtual void       OnOkClicked();

  private:
    AliMonitorStringDlg(const AliMonitorStringDlg& dlg);
    AliMonitorStringDlg& operator = (const AliMonitorStringDlg& /*dlg*/);

    TGLayoutHints*     fStringLayout;    // layout of the text entry
    TGLabel*           fStringLabel;     // label for the text entry
    TGTextEntry*       fStringEntry;     // the text enty

    TString&           fString;          // result
  };


  class AliMonitorNumberDlg : public AliMonitorDialog {

  public:
    AliMonitorNumberDlg(Float_t& value, TGFrame* main, const char* title,
			const char* label, Float_t min);
    virtual ~AliMonitorNumberDlg();

    virtual void       OnOkClicked();

  private:
    AliMonitorNumberDlg(const AliMonitorNumberDlg& dlg);
    AliMonitorNumberDlg& operator = (const AliMonitorNumberDlg& /*dlg*/);

    TGLayoutHints*     fNumberLayout;    // layout of the number entry
    TGLabel*           fNumberLabel;     // label for the number entry
    TGNumberEntry*     fNumberEntry;     // the number entry

    Float_t&           fNumber;          // result
  };


  ClassDef(AliMonitorClient, 0)   // class for receiving and displaying monitor histograms
};
 

#endif









