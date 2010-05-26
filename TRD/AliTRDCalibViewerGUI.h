#ifndef ALITRDCALIBVIEWERGUI_H
#define ALITRDCALIBVIEWERGUI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalibViewerGUI.h 34418 2009-08-26 15:47:50Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Class which implements AliBaseCalibViewerGUI for the TRD                 //
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
#include <TString.h>
#include "AliBaseCalibViewer.h"
#include "AliTRDCalibViewer.h"
#include "AliBaseCalibViewerGUI.h"
class TROOTt;
       
class AliTRDCalibViewerGUI : public AliBaseCalibViewerGUI {
   
public:
   AliTRDCalibViewerGUI(const TGWindow *p, UInt_t w, UInt_t h, char* fileName);  // constructor; fileName specifies the ROOT tree used for drawing
   AliTRDCalibViewerGUI(const AliTRDCalibViewerGUI &c);                          // copy constructor
   AliTRDCalibViewerGUI &operator = (const AliTRDCalibViewerGUI &param);         // assignment operator

   virtual ~AliTRDCalibViewerGUI();
   static void ShowGUI();                                                      // show the GUI
   static void ShowGUI(const Char_t* treeFile, const Char_t* treeName="TRDcalibDetails");  // show the GUI
   static void ShowGUIwithTrending();          // show a GUI with 2 tabs: 1 tab for time/run dependent observables
                                               //                        1 tab for run details
   
   virtual void DrawGUI(const TGWindow *p, UInt_t w, UInt_t h);         // to be called by the costructor, here the windows is drawn
   virtual void Initialize(const char* fileName, const char* treeName = "TRDcalibDetails"); // initializes the GUI with default settings and opens tree for drawing
   virtual void Initialize(AliBaseCalibViewer *viewer);                  // initializes the GUI with default settings and opens tree for drawing
   virtual void Reload(){Initialize(fViewer);}                          // reload the viewr after it has been changed, e.g. added a new referenceTree, ...
   virtual void Reset();
   virtual Bool_t CreateDetailsTree(Int_t run, const Char_t* outFile, const Char_t* ocdbStorage="nothing");   // create a tree with pad level info for a given run

   virtual TString* GetDrawString();                                    // create the draw string out of selection
   virtual TString* GetCutString();                                     // create the cut string out of selection
   virtual TString* GetSectorString();                                  // create the sector string out of selection
 //  virtual AliBaseCalibViewer* GetViewer() {return fViewer;}             // returns the internal AliTPCCalibViewer object, which does the work
   virtual void DoDraw();                            // main method for drawing according to user selection
   virtual void MouseMove(Int_t event, Int_t x, Int_t y, TObject *selected); 
   void SetTree();
   void HandleFilesystem();

protected:      
   // TRD specific widgets; these are added to the ones painted by the base class
   TGCompositeFrame    *fContLayer;               // container for layer GUI elements   -
   TGCompositeFrame    *fContSector;              // container for sector GUI elements  -
   TGCompositeFrame    *fContStack;               // container for stack GUI elements    -
   TGLabel             *fLblLayer;                // Layer label                          -
   TGLabel             *fLblSector;               // Sector label                        -
   TGLabel             *fLblStack;                // Stack label                         -
   TGNumberEntry       *fNmbLayer;                // Layer number entry                   -
   TGNumberEntry       *fNmbSector;               // Sector number entry                 -
   TGNumberEntry       *fNmbStack;                // Stack number entry                  -
   TGCompositeFrame    *fContLoad;                // 'Load Run' frame container          -
   TGCompositeFrame    *fContRun;                 // run frame                           -
   TGLabel             *fLblRun;                  // run label                           -
   TGNumberEntry       *fNmbRun;                  // run entry                           -
   TGCompositeFrame    *fContStorage;             // OCDB storage frame                   -
   TGLabel             *fLblStorage;              // OCDB label                           -
   TGTextEntry         *fTxtStorage;              // OCDB text entry                      -
   TGCompositeFrame    *fContVersion;	          // version frame                        -
   TGLabel             *fLblVersion;              // version label                         -
   TGNumberEntry       *fNmbVersion;              // version number entry                 -
   TGCompositeFrame    *fContSubVersion;	  // sub-version frame                  -
   TGLabel             *fLblSubVersion;           // sub-version label                  -
   TGNumberEntry       *fNmbSubVersion;           // sub-version number entry           -
   TGCompositeFrame    *fContChecks;              // check buttons container            -
   TGCheckButton       *fChkCalibs;               // check calibration                  -
   TGCheckButton       *fChkDCS;                  // check DCS                        -
   TGCheckButton       *fChkAlign;                // check alignment                   -
   TGTextButton        *fBtnLoad;                 // load button                        -
   // frame for loading an array of AliTRDCalPad objects
   TGCompositeFrame    *fContLoadCalibObjects;    //  container for loading calib objects
   TGCompositeFrame    *fContCalibInput;          //  container for the input file
   TGLabel             *fLblCalibInputFilename;   //   label for the input file
   TGTextEntry         *fTxtCalibInputFilename;   //   text entry for the input file
   TGCompositeFrame    *fContCalibOutput;         //   container for output file
   TGLabel             *fLblCalibOutputFilename;  //   label for output file
   TGTextEntry         *fTxtCalibOutputFilename;  //   text entry for output file
   TGTextButton        *fBtnLoadCalibObjects;     //   load button
   
   ClassDef(AliTRDCalibViewerGUI, 1)
};

#endif
