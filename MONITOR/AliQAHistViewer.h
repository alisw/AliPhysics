/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
///////////////////////////////////////////////////////////////////////////////
//
//     QA histogram viewer
//     (see AliQAHistNavigator.cxx for details)
//     Origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

#ifndef ALIQAHISTVIEWER_H
#define ALIQAHISTVIEWER_H

#include <TGFrame.h>

class TRootEmbeddedCanvas;
class AliQAHistNavigator;
class TGComboBox;
class TGCheckButton;
class TGWindow;

class AliQAHistViewer : public TGMainFrame {

private:
    TRootEmbeddedCanvas  *fEcan; //embedded canvas
    AliQAHistNavigator   *fQANavigator; //the navigator engine
    TGComboBox           *fFileListBox; //drop down menu
    TGComboBox           *fDetectorListBox; //drop down menun
    TGComboBox           *fLevelListBox; //drop down menu
    TGComboBox           *fHistListBox; //drop down menu
    TGCheckButton        *fExpertMode;// expertmode
    Bool_t               fIsEmbedded; //whether the window is embedded somewhere
    AliQAHistViewer(const AliQAHistViewer&);            // Not implemented
    AliQAHistViewer& operator=(const AliQAHistViewer&); // Not implemented
   
public:
    AliQAHistViewer(const TGWindow *p, UInt_t w=500, UInt_t h=500, Bool_t embed=kFALSE);
    virtual ~AliQAHistViewer();
    void DoExit();
    void DoDrawNext();
    void DoDrawPrev();
    void DoSetFile(Int_t s);
    void DoSetDetector(Int_t s);
    void DoSetLevel(Int_t s);
    void DoSetHistogram(Int_t s);
    void DoSetExpertMode(Bool_t mode=kTRUE);
    void FillComboBoxWithListEntries( TGComboBox* box, const TList* list );
    void UpdateAllPathComboBoxes();
   
    ClassDef(AliQAHistViewer, 999)
};

#endif
