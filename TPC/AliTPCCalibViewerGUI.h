#ifndef ALITPCCALIBVIEWERGUI
#define ALITPCCALIBVIEWERGUI

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCCalibViewerGUI.h,v */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  GUI for the AliTPCCalibViewer                                            //
//  used for the calibration monitor                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TGButton
#include "TGWidget.h"
#endif
#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif

#include <TGButton.h>
#include <TGListBox.h>
#include <TGComboBox.h>
#include <TGNumberEntry.h>
#include <TRootEmbeddedCanvas.h>
#include <TGSplitter.h>
#include <TGButtonGroup.h>
#include <TGLabel.h>
#include <TGTab.h>


#include <iostream>
#include "AliTPCCalibViewer.h"


class AliTPCCalibViewerGUI : public TGCompositeFrame {
protected:
   AliTPCCalibViewer   *fViewer;             // CalibViewer object used for drawing

   TGCompositeFrame    *fContTopBottom;      // container for all GUI elements, vertical divided
   TGCompositeFrame    *fContLCR;            // container for all GUI elements, horizontal divided
   TGCompositeFrame    *fContLeft;           // container for GUI elements on left side
   TGCompositeFrame    *fContRight;          // container for GUI elements on right side
   TGCompositeFrame    *fContCenter;         // container for GUI elements at the center
   TGCompositeFrame    *fContPlotOpt;        // container for plot options GUI elements
   TGCompositeFrame    *fContDrawOpt;        // container for draw options GUI elements
   TGCompositeFrame    *fContDrawOptSub1D2D; // container for 1D and 2D radio-button
   TGCompositeFrame    *fContNormalized;     // container for normalization options GUI elements
   TGCompositeFrame    *fContCustom;         // container for custom draw command GUI elements
   TGCompositeFrame    *fContCuts;           // container for cut options GUI elements
   TGCompositeFrame    *fContSector;         // container for sector GUI elements
   TGCompositeFrame    *fContAddCuts;        // container for additional cut command GUI elements
   TGCompositeFrame    *fContFit;            // container for fit GUI elements
   TGCompositeFrame    *fContAddFit;         // container for additional fit GUI elements
   TGCompositeFrame    *fContScaling;        // container for scaling GUI elements
   TGCompositeFrame    *fContSetMax;         // container for SetMaximum elements
   TGCompositeFrame    *fContSetMin;         // container for SetMinimum elements
   TGListBox           *fListVariables;      // listbox with possible variables
   TGTextButton        *fBtnDraw;            // draw button
   TGTextButton        *fBtnFit;             // fit button
   TGTextButton        *fBtnAddFitFunction;  // button to add fit function to normalization
   TGTextButton        *fBtnGetMinMax;       // GetMinMax-button
   TRootEmbeddedCanvas *fCanvMain;           // main drawing canvas
   TGRadioButton       *fRadioRaw;           // raw radio button
   TGRadioButton       *fRadioNormalized;    // normalized radio button
   TGRadioButton       *fRadioPredefined;    // predefined plot radio button
   TGRadioButton       *fRadioCustom;        // custom radio button
   TGRadioButton       *fRadio1D;            // 1D radio button
   TGRadioButton       *fRadio2D;            // 2D radio button
   TGRadioButton       *fRadioTPC;           // TPC radio button
   TGRadioButton       *fRadioSideA;         // side A radio button
   TGRadioButton       *fRadioSideC;         // side C radio button
   TGRadioButton       *fRadioSector;        // sector radio button
   TGCheckButton       *fChkAuto;            // automatic redraw checkbox
   TGComboBox          *fComboMethod;        // normalization methods dropdown box
   TGListBox           *fListNormalization;  // listbox with possible normalization variables
   TGComboBox          *fComboCustom;        // combo box for custom draw commands
   TGNumberEntry       *fNmbSector;          // number entry box for specifying the sector
   TGLabel             *fLblSector;          // label that shows the active sector
   TGCheckButton       *fChkAddCuts;         // additional cuts check box
   TGComboBox          *fComboAddCuts;       // additional cuts combo box
   TGComboBox          *fComboCustomFit;     // custom fit combo box
   TGCheckButton       *fChkSetMax;          // Set maximum check box
   TGCheckButton       *fChkSetMin;          // Set maximum check box
   TGCheckButton       *fChkGetMinMaxAuto;   // Get Min & Max automatically from plot
   TGTextEntry         *fTxtSetMax;          // custom maximum text box
   TGTextEntry         *fTxtSetMin;          // custom minimum text box
   

   void Initialize(char* fileName);          // initializes the GUI with default settings and opens tree for drawing
   
public:
   AliTPCCalibViewerGUI(const TGWindow *p, UInt_t w, UInt_t h, char* fileName);  // constructor; fileName specifies the ROOT tree used for drawing
   AliTPCCalibViewerGUI(const AliTPCCalibViewerGUI &c);                          // copy constructor
   AliTPCCalibViewerGUI &operator = (const AliTPCCalibViewerGUI &param);         // assignment operator

   virtual ~AliTPCCalibViewerGUI();
   // virtual void CloseWindow();

   void HandleButtons(Int_t id = -1);        // handles mutual radio button exclusions
   void DoNewSelection();                    // decides whether to redraw if user makes another selection
   void DoDraw();                            // main method for drawing according to user selection
   void DoFit();                             // main method for fitting
   void GetMinMax();                         // Read current Min & Max from the plot and set it to fTxtSetMin & fTxtSetMax
   void ChangeSector();                      // function that is called, when the number of the sector is changed
   void AddFitFunction();                    // adds the last fit function to the normalization list
   static void ShowGUI(const char* fileName); //initialize and show GUI for presentation
   ClassDef(AliTPCCalibViewerGUI, 0)
};

#endif
