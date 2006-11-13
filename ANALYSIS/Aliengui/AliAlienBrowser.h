#ifndef ALIALIENBROWSER_H
#define ALIALIENBROWSER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliAlienBrowser
//   AliAlienBrowser class that describes the alien browser of the GUI
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliAlienBrowser                               //
//                                                                      //
//                      AliEn browser of the GUI.                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <TGCanvas.h>

class TGFrame;
class TGListTree;
class TGListTreeItem;

enum EBrowseType{
   kGridBrowse,
   kLocalBrowse
};

//___________________________________________________________________________
class AliAlienBrowser : public TGCanvas {
 public:
  AliAlienBrowser(const TGWindow* p, UInt_t w, UInt_t h, TGFrame* frame, const char* objectToConnect, EBrowseType type);
  ~AliAlienBrowser();
  
  void        AddItem(TGListTreeItem* parent, const char* txt);
  void        OnDoubleClick(TGListTreeItem* item, Int_t btn);
  const char* GetPath();
  void        GotoDir(const char* dir);

  EBrowseType GetBrowseType() const {return fBrowseType;}

  //___________________________________________________________________________
 private:   
  AliAlienBrowser(const AliAlienBrowser&); // copy ctor
  AliAlienBrowser& operator= (const AliAlienBrowser&); // assignment operator
  
  TString     DirName(TGListTreeItem* item) const;
  const char* GetPath(TGListTreeItem *item) const;
  void        Refresh() const;
  
  void        OnDoubleClickGrid(TGListTreeItem* item, Int_t btn);
  void        OnDoubleClickLocal(TGListTreeItem* item, Int_t btn);
  
  
  TGFrame     *fFrame; //main browser frame     
  TGListTree  *fListTree; //tree structure
  
  EBrowseType  fBrowseType; // whether is for Local or Grid browsing
  
  ClassDef(AliAlienBrowser, 0) // AliAlienBrowser
};

#endif
