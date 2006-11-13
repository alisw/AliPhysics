#ifndef ALIPACKAGEFRAME_H
#define ALIPACKAGEFRAME_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliPackageFrame
//   AliPackageFrame class that describes the package frame of the GUI
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliPackageFrame                               //
//                                                                      //
//                      Package tab of the GUI.                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TGFrame.h>

class TGLabel;
class TGVerticalFrame;
class TGTextEntry;
class TGButton;

class AliAnalysisGUI;

//___________________________________________________________________________
class AliPackageFrame : public TGHorizontalFrame {

 public:
  AliPackageFrame(const TGWindow *main, UInt_t w, UInt_t h, AliAnalysisGUI*);
   
  //___________________________________________________________________________
  //slots
  void OnBuild();
  void OnSelect();

  //___________________________________________________________________________
 private:
  AliPackageFrame(const AliPackageFrame&);
  AliPackageFrame& operator= (const AliPackageFrame&);
  
  void CreatePARFile(const char* parfile);
  
  TGVerticalFrame     *fVFrame1, *fVFrame2; //vertical frames 
  TGLabel             *fLabel1;  //package label  
  TGTextEntry         *fTextPackage; //package text box
  TGButton            *fButtonSelect; //select button
  TGButton            *fButtonBuild; //run button
  
  AliAnalysisGUI      *fAliAnalysisGUI; //analysis gui pointer
  
  ClassDef(AliPackageFrame, 0); // AliPackageFrame
};


#endif
