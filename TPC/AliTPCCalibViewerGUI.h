#ifndef ALITPCCALIBVIEWERGUI_H
#define ALITPCCALIBVIEWERGUI_H

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
class AliTPCCalibViewer;


// class TGListBox;
// class TGNumberEntry;
// class TGSplitter;
// class TGTab;
// class TGWidget; // ???
// class TGLabel;
// class TGButtonGroup;
// class TGComboBox;
// class TRootEmbeddedCanvas;
// class TGButton;
// class TGRadioButton;
// class TGCheckButton;
// class TGTextEntry;
       
       
class AliTPCCalibViewerGUI : public TGCompositeFrame {
   
public:
   AliTPCCalibViewerGUI(const TGWindow *p, UInt_t w, UInt_t h, char* fileName);  // constructor; fileName specifies the ROOT tree used for drawing
   AliTPCCalibViewerGUI(const AliTPCCalibViewerGUI &c);                          // copy constructor
   AliTPCCalibViewerGUI &operator = (const AliTPCCalibViewerGUI &param);         // assignment operator

   virtual ~AliTPCCalibViewerGUI();
   // virtual void CloseWindow();

   void HandleButtonsGeneral(Int_t id = -1); // handles mutual radio button exclusions for general Tab
   void HandleButtons1D(Int_t id = -1);      // handles mutual radio button exclusions for 1D Tab
   void HandleButtons2D(Int_t id = -1);      // handles mutual radio button exclusions for 2D Tab
   void HandleButtonsStat(Int_t id = -1);    // handles statistic check boxes 
   void HandleButtonsRight(Int_t id = -1);   // handles mutual radio button exclusions for right side
   void DoNewSelection();                    // decides whether to redraw if user makes another selection
   void DoDraw();                            // main method for drawing according to user selection
   void DoFit();                             // main method for fitting
   void GetMinMax();                         // Read current Min & Max from the plot and set it to fTxtSetMin & fTxtSetMax
   void ChangeSector();                      // function that is called, when the number of the sector is changed
   void AddFitFunction() const;              // adds the last fit function to the normalization list
   static void ShowGUI(const char* fileName); //initialize and show GUI for presentation
   
   
protected:
   AliTPCCalibViewer   *fViewer;             // CalibViewer object used for drawing

   TGCompositeFrame    *fContTopBottom;      // container for all GUI elements, vertical divided
   TGCompositeFrame    *fContLCR;            // container for all GUI elements, horizontal divided
   TGCompositeFrame    *fContLeft;           // container for GUI elements on left side
   TGTab               *ftabLeft;            // Tabs on the left side for plot options
   TGCompositeFrame    *ftabLeft0;           // Tab 0 on the left side for general plot options
   TGCompositeFrame    *ftabLeft1;           // Tab 1 on the left side for 1D plot options
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
   TGCompositeFrame    *fContAddDrawOpt;     // additional draw options container
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
   TGComboBox          *fComboAddDrawOpt;    // additional draw options combo box
   TGCheckButton       *fChkAuto;            // automatic redraw checkbox
   TGComboBox          *fComboMethod;        // normalization methods dropdown box
   TGListBox           *fListNormalization;  // listbox with possible normalization variables
   TGComboBox          *fComboCustom;        // combo box for custom draw commands
   TGCheckButton       *fChkAddDrawOpt;      // additional draw options check box
   TGNumberEntry       *fNmbSector;          // number entry box for specifying the sector
   TGLabel             *fLblSector;          // label that shows the active sector
   TGCheckButton       *fChkCutZero;         // cut zeros check box
   TGCheckButton       *fChkAddCuts;         // additional cuts check box
   TGComboBox          *fComboAddCuts;       // additional cuts combo box
   TGComboBox          *fComboCustomFit;     // custom fit combo box
   TGCheckButton       *fChkSetMax;          // Set maximum check box
   TGCheckButton       *fChkSetMin;          // Set maximum check box
   TGCheckButton       *fChkGetMinMaxAuto;   // Get Min & Max automatically from plot
   TGTextEntry         *fTxtSetMax;          // custom maximum text box
   TGTextEntry         *fTxtSetMin;          // custom minimum text box
   TGGroupFrame        *fContDrawOpt1D;      // container in tabLeft1 
   TGCompositeFrame    *fcontDrawOpt1DSubLR; // container in tabLeft1 to divide L/R
   TGCompositeFrame    *fContDrawOpt1DSubNSC; // container in tabLeft1 for following radio buttons 
   TGRadioButton       *fRadioNorm;          // radio button for normal 1D drawing
   TGRadioButton       *fRadioSigma;         // radio button for sigma 1D drawing
   TGTextEntry         *fTxtSigmas;          // text box to specify sigmas
   TGCompositeFrame    *fContCumuLR;         // container in tabLeft1 for two colums for cumulative and integrative
   TGCompositeFrame    *fContCumLeft;        // container in tabLeft1 for cumulative, left
   TGCompositeFrame    *fContCumRight;       // container in tabLeft1 for cumulative, right
   TGLabel             *fLblSigmaMax;        // label to indicate sigmaMax
   TGTextEntry         *fTxtSigmaMax;        // text box to specify sigmaMax
   TGRadioButton       *fRadioCumulative;    // radio button for cumulative 1D drawing
   TGCheckButton       *fCheckCumulativePM;  // checkbox for plus/minus cumulative 1D drawing
   TGRadioButton       *fRadioIntegrate;     // radio button for integral 1D drawing
   TGCompositeFrame    *fContDrawOpt1DSubMML; // container in tabLeft1 for following check boxes
   TGCheckButton       *fChkMean;            // checkbox to plot mean
   TGCheckButton       *fChkMedian;          // checkbox to plot median
   TGCheckButton       *fChkLTM;             // checkbox to plot LTM
   TGGroupFrame        *fContStatOpt;        // container for statistic options in tabLeft1 
   TGCheckButton       *fChkStatName;        // checkbox to display histogram name in statistic legend
   TGCheckButton       *fChkStatEntries;     // checkbox to display entries in statistic legend
   TGCompositeFrame    *fContStatMean;       // container for mean and its error in stat opt
   TGCheckButton       *fChkStatMean;        // checkbox to display mean in statistic legend
   TGCheckButton       *fChkStatMeanPM;      // checkbox to display mean error in statistic legend
   TGCompositeFrame    *fContStatRMS;        // container for RMS and its error in stat opt
   TGCheckButton       *fChkStatRMS;         // checkbox to display RMS in statistic legend
   TGCheckButton       *fChkStatRMSPM;       // checkbox to display RMS error in statistic legend
   TGCheckButton       *fChkStatUnderflow;   // checkbox to display underflow error in statistic legend
   TGCheckButton       *fChkStatOverflow;    // checkbox to display overflow error in statistic legend
   TGCheckButton       *fChkStatIntegral;    // checkbox to display integral in statistic legend
   TGCompositeFrame    *fContStatSkew;       // container for skewness and its error in stat opt
   TGCheckButton       *fChkStatSkewness;    // checkbox to display skewness in statistic legend
   TGCheckButton       *fChkStatSkewnessPM;  // checkbox to display skewness error in statistic legend
   TGCompositeFrame    *fContStatKurt;       // container for kurtosis and its error in stat opt
   TGCheckButton       *fChkStatKurtosis;    // checkbox to display kurtosis in statistic legend
   TGCheckButton       *fChkStatKurtosisPM;  // checkbox to display kurtosis error in statistic legend
   
   void Initialize(char* fileName);          // initializes the GUI with default settings and opens tree for drawing
   
   ClassDef(AliTPCCalibViewerGUI, 0)
};

#endif
