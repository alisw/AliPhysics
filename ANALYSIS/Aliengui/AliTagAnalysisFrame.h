#ifndef ALITAGANALYSISFRAME_H
#define ALITAGANALYSISFRAME_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliTagAnalysisFrame
//   AliTagAnalysisFrame class that describes the event tag frame of the GUI
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliTagAnalysisFrame                           //
//                                                                      //
//                      Event tag tab of the GUI.                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TGFrame.h>

class TGListTreeItem;
class TGLabel;
class TGTextEntry;
class TGVerticalFrame;
class TGButton;
class TGGroupFrame;
class TGComboBox;
class TGListBox;
class TGTransientFrame;

class TGridResult;
class TChain;

class AliRunTagCuts;
class AliEventTagCuts;
class AliTagAnalysis;

class AliAlienBrowser;
class AliTagFrame;
class AliAnalysisGUI;

//___________________________________________________________________________
class AliTagAnalysisFrame : public TGMainFrame {
public:
  AliTagAnalysisFrame(const TGWindow *main, UInt_t w, UInt_t h, AliAnalysisGUI* a=0);
  ~AliTagAnalysisFrame();
  
  //___________________________________________________________________________
  // slots
  void LocalBrowse();
  void GridBrowse();
  void OnDoubleClick(TGListTreeItem* item, Int_t btn);
  void OnOKButton();
  void InsertTagCutsRangeLocal();
  void InsertTagCutsRangeGrid();
  void RunLocal();
  void RunGrid();
  void ProcessSelector(const char* selectorfile);

  //___________________________________________________________________________
 private:
  AliTagAnalysisFrame(const AliTagAnalysisFrame&); // cp ctor
  AliTagAnalysisFrame& operator= (const AliTagAnalysisFrame&); // op=
   
  // private methods
  void BuildLocalGroup (TGCompositeFrame* frame);
  void BuildGridGroup  (TGCompositeFrame* frame);
  void InsertTagCutsRange(Int_t id);
  void AddResult (const char* line);

  const Int_t fkNumberOfTags; //event tags

  TGVerticalFrame     *fVFrame1, *fVFrame2; //vertical frames
  TGGroupFrame        *fGroup1, *fGroup2, *fGroup3; //group of frames
   
  AliAnalysisGUI      *fAliAnalysisGUI; //analysis gui pointer
  AliTagFrame         *fTagFrame; //tag frame pointer
  AliAlienBrowser     *fAliEnBrowser; //alien browser pointer

  // local 
  TGLabel             *fLocalLabel1; //label - local tags
  TGTextEntry         *fLocalPath; //text box - local tags
  TGButton            *fLocalButton, *fButtonInsert, *fButtonRun; //buttons
  TGComboBox          *fComboEventTagCut; //combo box
     
  // Grid
  TGLabel             *fGridLabel1; //label - grid tags
  TGTextEntry         *fGridPath; //text box - grid tags
  TGButton            *fGridButton, *fButtonInsert2, *fButtonRun2; //buttons
  TGComboBox          *fComboEventTagCut2; //combo box

  TGridResult         *fTagResult; //grid result
  TChain              *fAnalysisChain; //tchain object

  TGListBox           *fListBox; //list box
  TGTransientFrame    *fBrowser; //frame
  TGButton            *fBrowserButton; //browse button

  // AliRoot Tag cut analysis 
  AliTagAnalysis      *fAliTagAnalysis; //alitaganalysis object
  AliRunTagCuts       *fAliRunCuts; //run cuts object
  AliEventTagCuts     *fAliEventCuts; //event cuts object
  
  const char          **fEventTagCutsName; //event tag names

  ClassDef(AliTagAnalysisFrame, 0); // Tag Analysis Frame
};

#endif
