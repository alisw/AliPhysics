#ifndef ALITAGFRAME_H
#define ALITAGFRAME_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliTagFrame
//   AliTagFrame class that describes the event tag frame of the GUI
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                           AliTagFrame                                //
//                                                                      //
//                      Event tag tab of the GUI.                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <TGFrame.h>

class TGCanvas;
class TGVerticalFrame;
class TGButton;
class TGNumberEntryField;

enum ETagRangeType {
   kRangeMin,
   kRangeMax,
   kRangeMinMax
};

//___________________________________________________________________________
class AliTagFrame : public TGTransientFrame {
public:
  AliTagFrame(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h, UInt_t options, const char* text, Int_t tagId, ETagRangeType range);
  ~AliTagFrame();

  Int_t GetRangeMin() const {return fMin;}
  Int_t GetRangeMax() const {return fMax;}

  void  OnClicked();
   
//___________________________________________________________________________
private:
  AliTagFrame(const AliTagFrame&); // copy ctor
  AliTagFrame& operator= (const AliTagFrame&); // assignment op

  // methods to build the GUI
   void  CreateTagName(const char* name);
   void  CreateTagRange(TGVerticalFrame* frame, TGNumberEntryField*& entry, const char *name);
   void  CreateTagButton();

   void (AliTagFrame::*fTagCutMethods [3]) (void); //tag fields
   
   Int_t              fMin;   // min range
   Int_t              fMax;   // max range
   ETagRangeType      fRange; // range type

   TGNumberEntryField *fEntry1, *fEntry2; //range entry fields

   TGButton           *fButton; //button
   TGVerticalFrame    *fVFrame1, *fVFrame2; //vertical frames

   ClassDef(AliTagFrame, 0) // Tag Frame
};

#endif
