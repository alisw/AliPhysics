#ifndef ALIMONITORCLIENT_H
#define ALIMONITORCLIENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TGFrame.h>
#include <TGMenu.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TGToolBar.h>
#include <TG3DLine.h>
#include <TGNumberEntry.h>
#include <TGCanvas.h>
#include <TGSplitter.h>
#include <TGListTree.h>
#include <TRootEmbeddedCanvas.h>
#include <TGTextView.h>
#include <RQ_OBJECT.h>
#include <TFolder.h>
#include <TSocket.h>
#include <TTimer.h>
#include <TFile.h>
#include "AliMonitorHisto.h"


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
  TFolder*           CreateTopFolder();
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
  void               ClearItems(TGListTreeItem* base);
  void               CleanUpTree(TGListTreeItem* base);
  void               UpdateTree();
  void               UpdateFavoritesTree();
  void               UpdateComparisonTree();
  void               UpdateDescription();
  void               UpdateHisto();
  void               UpdateAll();

  TGLayoutHints*     fMenuBarLayout;
  TGLayoutHints*     fMenuBarItemLayout;
  TGLayoutHints*     fMenuBarHelpLayout;
  TGPopupMenu*       fMenuFile;
  TGPopupMenu*       fMenuView;
  TGPopupMenu*       fMenuFavorites;
  TGPopupMenu*       fMenuReference;
  TGPopupMenu*       fMenuOptions;
  TGPopupMenu*       fMenuHelp;
  TGMenuBar*         fMenuBar;

  TGLayoutHints*     fToolBarLayout;
  TGHorizontal3DLine* fToolBarSep;
  TGToolBar*         fToolBar;
  TGLayoutHints*     fEventNumberLayout;
  TGNumberEntry*     fEventNumber;
  TGButton*          fEventButton;
  TGButton*          fSumButton;
  TGButton*          fRunButton;
  TGButton*          fLoopButton;
  const TGPicture*   fLoopOnPicture;
  const TGPicture*   fLoopOffPicture;
  TGButton*          fPreviousButton;
  TGButton*          fNextButton;
  TGButton*          fCopyButton;
  TGButton*          fSaveButton;
  TGButton*          fPrintButton;

  TGLayoutHints*     fBottomLayout;
  TGLayoutHints*     fLeftLayout;
  TGLayoutHints*     fExpandLayout;

  TGVerticalFrame*   fVerticalFrame;
  TGHorizontalFrame* fHorizontalFrame;

  TGCompositeFrame*  fTreeFrame;
  TGCanvas*          fTreeCanvas;
  TGListTree*        fTree;
  const TGPicture*   fHistoPicture;
  TGListTreeItem*    fAllItem;
  TGListTreeItem*    fFavoritesItem;
  TGListTreeItem*    fComparisonItem;

  TGVSplitter*       fTreeSplitter;

  TGCompositeFrame*  fDrawFrame;
  TRootEmbeddedCanvas* fDrawCanvas;

  TGHSplitter*       fDescriptionSplitter;

  TGCompositeFrame*  fDescriptionFrame;
  TGTextView*        fDescription;

  TString            fServerName;
  TSocket*           fSocket;
  TFileHandler*      fSocketHandler;

  TFolder*           fFolder;

  TGListTreeItem*    fCurrentItem;
  TGListTreeItem*    fBaseItem;
  TTimer*            fLoopTimer;
  Int_t              fLoopInterval;

  TString            fFavoritesFileName;

  TString            fReferenceFileName;
  TFolder*           fReference;

  TString            fPrintCommand;

  static const char* fgSettingsFileName;

  ClassDef(AliMonitorClient, 0)   // class for receiving and displaying monitor histograms
};
 

#endif









