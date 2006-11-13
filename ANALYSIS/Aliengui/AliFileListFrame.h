#ifndef ALIFILELISTFRAME_H
#define ALIFILELISTFRAME_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliFileListFrame
//   AliFileListFrame class that describes the file list frame of the GUI
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliFileListFrame                              //
//                                                                      //
//                      File list frame of the GUI.                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TGFrame.h>
class TGCanvas;

class TGHorizontalFrame;
class TGVerticalFrame;
class TGHorizontalFrame;
class TGLabel;
class TGTextEntry;
class TGNumberEntry;
class TGTextEntry;
class TGButton;
class TGCompositeFrame;
class TGCanvas;
class TGTableLayout;
class TGLVEntry;
  
class TObjArray;

//___________________________________________________________________________
class AliFileListFrame : public TGCompositeFrame {
 public:
  AliFileListFrame(const TGWindow *main, UInt_t w, UInt_t h);
  ~AliFileListFrame();
  
  void        SetQueryPath(const char* path);
  const char* GetQueryPath();
  const char* GetQueryPattern();
  
  // slots
  void OnDoubleClick(TGLVEntry* entry, Int_t btn);
  void RunQuery();
  
  //___________________________________________________________________________
 private:
  AliFileListFrame(const AliFileListFrame&); // cp ctor
  AliFileListFrame& operator= (AliFileListFrame&); // op= 

  // private methods
  TGHorizontalFrame   *fHFrame1; //horiz. frame 
  TGVerticalFrame     *fVFrame1, *fVFrame2; //vertical frames
  
  // VFrame1's widgets
  TGHorizontalFrame   *fHFrame2;  //vertical frame
  TGLabel             *fLabel1, *fLabel2, *fLabel3; //labels
  TGTextEntry         *fTextQueryPath; //text box
  TGNumberEntry       *fNumMaxResults; //results
  TGTextEntry         *fTextQueryPattern; //query pattern
  TGButton            *fButtonRun; //run button
  TGCompositeFrame    *fContents; //frame
  TGCanvas            *fCanvas; //canvas
  TGTableLayout       *fTableLayout; //layout
  
  TObjArray           *fTags; // File "detail" mode tags 
  
  void DisplayObject(const TString& fname,const TString& name) const;
  void BuildQueryPathFrame();

  AliFileListFrame & operator=(const AliFileListFrame & ) {return *this;}

  ClassDef(AliFileListFrame, 0) // AliFileListFrame  
};

#endif
