#ifndef ALISELECTORFRAME_H
#define ALISELECTORFRAME_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliSelectorFrame
//   AliSelectorFrame class that describes the selector frame of the GUI
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliSelectorFrame                              //
//                                                                      //
//                      Selector tab of the GUI.                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TGFrame.h>

class TGCanvas;
class TGVerticalFrame;
class TGTextEntry;
class TGButton;
class TGLabel;

class AliAnalysisGUI;

//___________________________________________________________________________
class AliSelectorFrame : public TGHorizontalFrame {

 public:
  AliSelectorFrame(const TGWindow *main, UInt_t w, UInt_t h, AliAnalysisGUI*, AliTagAnalysisFrame*);
  AliSelectorFrame(const AliSelectorFrame& fSelectorFrame);

  //___________________________________________________________________________
  //slots
  void OnSelect();
  void OnRun();

  //___________________________________________________________________________
 private:
  TGVerticalFrame     *fVFrame1, *fVFrame2; //vertical frames
  TGLabel             *fLabel1; //macro label
  TGTextEntry         *fTextSelector; //selector text box
  TGButton            *fButtonSelect; //select button
  TGButton            *fButtonRun; //run button
  
  AliAnalysisGUI      *fAliAnalysisGUI; //analysis gui pointer
  AliTagAnalysisFrame *fTagAnalysisFrame; //tag frame pointer
  
  AliSelectorFrame & operator=(const AliSelectorFrame & ) {return *this;}

  ClassDef(AliSelectorFrame, 0); // SelectorFrame
};


#endif
