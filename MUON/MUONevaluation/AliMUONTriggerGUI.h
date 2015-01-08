#ifndef ALIMUONTRIGGERGUI_H
#define ALIMUONTRIGGERGUI_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup evaluation
/// \class AliMUONTriggerGUI
/// \brief Graphical User Interface utility class for the MUON trigger detector
//  Author Bogdan Vulpescu, LPC Clermont-Ferrand


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
class AliRawReader;
class AliMUONCalibrationData;
class AliMUONTriggerGUIboard;
class AliMUONTriggerGUIdimap;
class AliMUONTriggerElectronics;
class AliMUONTriggerCrateStore;
class AliMUONMCDataInterface;
class AliMUONDigitStoreV2R;
class AliMUONDigitStoreV1;
class AliMUONTriggerStoreV1;

class AliMUONTriggerGUI : public TObject
{

public:

  AliMUONTriggerGUI(Int_t runNumber = 0);
  virtual ~AliMUONTriggerGUI() { 
    /// main gui destructor 
  };
  
  void OpenBoard(Int_t id);
  void HandleMenu(Int_t id);

  void CreateDigitStore();
  void PrintDigitStore() const;
  void ClearDigitStore();

  void CreateTriggerStore();
  void PrintTriggerStore() const;
  void ClearTriggerStore();

  void WriteTriggerRawData();

  void FET(Int_t onoff);
  void FETboard(Int_t ib, Int_t amp);

  /* Do functions */

  void DoRunGalApply();
  void DoRunRawApply();
  void DoRunCancel();
  void DoControlClose();
  void DoCircuitCancel();
  void DoCircuitOpen();
  void DoErrorOK();
  void DoNextEvent();
  void DoPreviousEvent();
  void DoSkipToEvent();
  void DoErrorGUI(const Char_t *wt);
  void DoFETRegRun(Int_t onoff);
  void DoFETRegOnCancel();
  void DoFETRegOffCancel();
  void DoFETRegOnRun();
  void DoFETRegOffRun();

  /* Close functions */

  void CloseWindow();
  void CloseRunInput() const;
  void CloseError() const;
  void CloseControl() const;
  void CloseCircuit() const;
  void CloseFETRegOn() const;
  void CloseFETRegOff() const;

private:
    
  AliMUONTriggerGUI (const AliMUONTriggerGUI& board); ///< copy constructor
  AliMUONTriggerGUI& operator=(const AliMUONTriggerGUI& board);///< assignment operator
  
private:
  
  /// nr of boards, nr of chambers
  enum { kNBoards = 234, kNMT = 4 }; 

  /// working status flags
  enum EMenuIdentifiers {
    
    kMFILEEXIT,
    kMFILERUN,
    kMFILECNTRL,
    
    kMMAPDIGITS,
    kMRESETDIGITS,

    kMCIRCUITOPEN,

    kMFETON,
    kMFETOFF,
    kMFETREGON,
    kMFETREGOFF,

    kMDSTORE,
    kMTSTORE,
    kMDSTORECL,
    kMTSTORECL,
    kMDSTOREP,
    kMTSTOREP,

    kMTRAWDATA

  };

  enum {
    kGood = 0x0001, kWithProblems = 0x0002, kNotWorking = 0x0004, kUnknown = 0x0008
  };                                ///< working status flags

  TGMainFrame      *fMain;          ///< The main frame
  TGImageMap       *fImageMap;      ///< The image map of the main frame
  TGTextBuffer     *fTxtBuffer1;    ///< Path to the data (galice.root)
  TGTextBuffer     *fTxtBuffer2;    ///< Current event number
  TGTextBuffer     *fTxtCircuit;    ///< Circuit to open
  TGTextBuffer     *fTxtFETRegOn;   ///< Regional crate to FET ON
  TGTextBuffer     *fTxtFETRegOff;  ///< Regional crate to FET OFF

  TGTransientFrame *fRunInput;      ///< Run input window
  TGTransientFrame *fError;         ///< Error window
  TGTransientFrame *fControl;       ///< Run control window
  TGTransientFrame *fCircuit;       ///< Circuit window
  TGTransientFrame *fFETRegOn;      ///< FET ON for a regional card
  TGTransientFrame *fFETRegOff;     ///< FET OFF for a regional card

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
  AliMUONTriggerCrateStore *fCrateManager;    ///< trigger boards manager
  AliMUONMCDataInterface *fMCDataInterface;   ///< MC data interface

  Bool_t            fBoardsInit;    ///< Control the InitBoards only once
  Bool_t            fControlOn;     ///< If the control frame is open

  AliMUONTriggerGUIdimap *fDiMap;   ///< Digits map

  AliMUONTriggerElectronics *fTriggerProcessor;   ///< The GUI trigger processor

  TObjArray *fBoards;               ///< The array of trigger boards

  AliMUONDigitStoreV2R *fDigitStore;    ///< GUI digit store (DSET)
  AliMUONTriggerStoreV1* fTriggerStore; ///< Trigger store with GUI digit store
  Bool_t fTStoreOn;                   ///< True if DSET trigger store has data
  Bool_t fRUNRAW;                     ///< True if run with raw data (root)
  AliRawReader *fRawReader;           ///< Reader for raw data input
  Int_t fCurrentRawEvent;             ///< Current event for raw data input
  AliMUONDigitStoreV1 *fRawDigitStore; ///< Digit store from raw data
  AliMUONTriggerStoreV1 *fRawTriggerStore; ///< Trigger store from raw data

  /// Access the array of trigger boards
  TObjArray *Boards() {
    if(!fBoards) fBoards = new TObjArray(kNBoards); return fBoards;
  };                                
  AliMUONTriggerGUIboard *GetBoard(Int_t id) const;

  virtual void Init();
  virtual void InitBoards();

  void  SetStripBoxes(AliMUONTriggerGUIboard *board);

  ClassDef(AliMUONTriggerGUI,2)      // Main GUI class for the MUON trigger

};

#endif
