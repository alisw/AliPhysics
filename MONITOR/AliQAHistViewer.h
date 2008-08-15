/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//
//     (see AliQAHistNavigator.cxx for details)
//
//     Origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

#ifndef ALIQAHISTVIEWER_H
#define ALIQAHISTVIEWER_H

#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TGStatusBar.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TString.h>
#include <TGComboBox.h>
#include "AliQAHistNavigator.h"

class AliQAHistViewer : public TGMainFrame {

private:
    TRootEmbeddedCanvas  *fEcan; //embedded canvas
    AliQAHistNavigator   *fQANavigator; //the navigator engine
    TGComboBox           *fFileListBox; //drop down menu
    TGComboBox           *fDetectorListBox; //drop down menun
    TGComboBox           *fLevelListBox; //drop down menu
    TGComboBox           *fHistListBox; //drop down menu
    Bool_t               fIsEmbedded; //whether the window is embedded somewhere
   
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
    void FillComboBoxWithListEntries( TGComboBox* box, const TList* list );
    void UpdateAllPathComboBoxes();
   
    ClassDef(AliQAHistViewer, 999)
};

#endif
