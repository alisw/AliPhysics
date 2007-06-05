#ifndef ALIMUONTRIGGERGUI_H
#define ALIMUONTRIGGERGUI_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup evaluation
/// \class AliMUONTriggerGUI
/// \brief Trigger GUI utility class

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Graphical User Interface utility class for the MUON trigger          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TObjArray.h>

class TString;
class TGMainFrame;
class TGTransientFrame;
class TGImageMap;
class TGTextEntry;
class TGTextBuffer;
class TRootEmbeddedCanvas;
class TParticle;
class TGTableLayout;

class AliStack;
class AliLoader;
class AliRunLoader;
class AliCDBManager;
class AliMUONCalibrationData;
class AliMUONTriggerGUIboard;
class AliMUONTriggerGUIdimap;
class AliMUONData;
class AliMUONTriggerElectronics;

class AliMUONTriggerGUI : public TObject
{

public:

  AliMUONTriggerGUI();
  virtual ~AliMUONTriggerGUI() { 
    /// main gui destructor 
  };
  
  AliMUONTriggerGUI (const AliMUONTriggerGUI& board);
  AliMUONTriggerGUI& operator=(const AliMUONTriggerGUI& board);

  void OpenBoard(Int_t id);
  void HandleMenu(Int_t id);

  /* Do functions */

  void DoRunApply();
  void DoRunCancel();
  void DoControlClose();
  void DoCircuitCancel();
  void DoCircuitOpen();
  void DoErrorOK();
  void DoNextEvent();
  void DoPreviousEvent();
  void DoSkipToEvent();
  void DoErrorGUI(const Char_t *wt);

  /* Close functions */

  void CloseWindow();
  void CloseRunInput() const;
  void CloseError() const;
  void CloseControl() const;
  void CloseCircuit() const;

private:
  
  enum { kNBoards = 234, kNMT = 4 };    ///< nr of boards, nr of chambers

  enum EMenuIdentifiers {
    
    kMFILEEXIT,
    kMFILERUN,
    kMFILECNTRL,
    
    kMMAPDIGITS,
    kMRESETDIGITS,

    kMCIRCUITOPEN,

    kMTRIGGERDSET

  };                                    ///< gui menu identifiers

  enum {
    kGood = 0x0001, kWithProblems = 0x0002, kNotWorking = 0x0004, kUnknown = 0x0008
  };                                ///< working status flags

  TGMainFrame      *fMain;          ///< The main frame
  TGImageMap       *fImageMap;      ///< The image map of the main frame
  TGTextBuffer     *fTxtBuffer1;    ///< Path to the data (galice.root)
  TGTextBuffer     *fTxtBuffer2;    ///< Current event number
  TGTextBuffer     *fTxtCircuit;    ///< Circuit to open

  TGTransientFrame *fRunInput;      ///< Run input window
  TGTransientFrame *fError;         ///< Error window
  TGTransientFrame *fControl;       ///< Run control window
  TGTransientFrame *fCircuit;       ///< Circuit window

  TGTextEntry      *fSkipToEventTxt;///< Control field shows current event number

  TString          *fFileName;      ///< Full galice file name
  TString          *fPath;          ///< Path string to galice
  TString          *fEvString;      ///< Event number string

  Int_t             fChamber;       ///< Current MT chamber
  Int_t             fEvent;         ///< Current event number
  Int_t             fEventsPerRun;  ///< Number of events per file (run)

  AliLoader        *fLoader;        ///< The MUON loader
  AliRunLoader     *fRunLoader;     ///< The run loader

  AliCDBManager    *fCDBManager;    ///< Calibration DB manager
  AliMUONCalibrationData *fCalibrationData;   ///< Calibration data for MUON

  Bool_t            fBoardsInit;    ///< Control the InitBoards only once

  AliMUONTriggerGUIdimap *fDiMap;   ///< Digits map

  AliMUONTriggerElectronics *fTriggerProcessor;   ///< The GUI trigger processor
  AliMUONData      *fMUONData;      ///< The MUON data manager

  TObjArray *fBoards;               ///< The array of trigger boards
  TObjArray *Boards() {
    if(!fBoards) fBoards = new TObjArray(kNBoards); return fBoards;
  };                                ///< Access the array of trigger boards
  AliMUONTriggerGUIboard *GetBoard(Int_t id) const;

  virtual void Init();
  virtual void InitBoards();

  void  SetStripBoxes(AliMUONTriggerGUIboard *board);

  ClassDef(AliMUONTriggerGUI,1)      // Main GUI class for the MUON trigger

};

#endif
