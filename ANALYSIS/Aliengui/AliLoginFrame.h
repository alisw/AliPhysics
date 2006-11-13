#ifndef ALILOGINFRAME_H
#define ALILOGINFRAME_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliLoginFrame
//   AliPackageFrame class that describes the login frame of the GUI
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliLoginFrame                                 //
//                                                                      //
//                      Login frame of the GUI.                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TGFrame.h>

class TGLabel;
class TGButton;
class TGTextEntry;
class TGVerticalFrame;
class TGHorizontalFrame;

//___________________________________________________________________________
class AliLoginFrame : public TGTransientFrame {
  
 public:
  AliLoginFrame(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h, UInt_t options = kVerticalFrame);
  ~AliLoginFrame();
  
  void DoLogIn();
  void DoCancel();
  
  //___________________________________________________________________________
 private:
  AliLoginFrame(const AliLoginFrame&); // cp ctor
  AliLoginFrame& operator= (const AliLoginFrame&); // op=
  
  TGLabel            *fLabel1, *fLabel2; //labels
  TGTextEntry        *fTextServer, *fTextUsername; //server - username text box
  TGButton           *fButtonLogIn, *fButtonCancel; //login & cancel buttons
  TGVerticalFrame    *fVFrame1; //vertical frame
  TGHorizontalFrame  *fHFrame1, *fHFrame2, *fHFrame3; //horizontal frames
  
  ClassDef(AliLoginFrame, 0); // LogIn Frame
};

#endif
