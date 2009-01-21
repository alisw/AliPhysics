#ifndef ALIMUONTRIGGERGUIBDMAP_H
#define ALIMUONTRIGGERGUIBDMAP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup evaluation
/// \class AliMUONTriggerGUIbdmap
/// \brief Trigger GUI utility class: single board map of the strips/digits
//  Author Bogdan Vulpescu, LPC Clermont-Ferrand

#include <TGFrame.h>

class TCanvas;
class TGCheckButton;
class TGTextEdit;
class TPolyLine;
class TBox;
class TPaveText;
class TObjArray;
class TH1F;
class TLatex;
class TGTableLayout;
class TGLabel;

class AliMUONTriggerGUIboard;
class AliMUONTriggerGUI;
class AliMUONTriggerCircuit;
class AliMUONTriggerCrateStore;
class AliMUONMCDataInterface;
class AliMUONDigitStoreV1;
class AliMUONTriggerStoreV1;
class AliMUONCalibrationData;
class AliLoader;

class AliMUONTriggerGUIbdmap : public TGFrame
{

public:

  AliMUONTriggerGUIbdmap(const TGWindow *p, const TGWindow *mainWindow, UInt_t w, UInt_t h);
  virtual ~AliMUONTriggerGUIbdmap();
  
  /// set the name of the board gui window
  void SetName(const Char_t *name)         { fMain->SetWindowName(name); };
  /// set the board associated to this instance
  void SetBoard(AliMUONTriggerGUIboard *b) { fBoard = b; };  
  /// set the board associated to this instance, from boards array and id
  void SetBoard(TObjArray *boards, Int_t id) { 
    fBoards = boards;
    fBoard  = (AliMUONTriggerGUIboard*)boards->UncheckedAt(id); }
  /// set the current muon loader
  void SetLoader(AliLoader *loader)        { fLoader = loader; };
  /// set the MC data interface
  void SetMCDataInterface(AliMUONMCDataInterface *mc) { fMCDataInterface = mc; };
  /// set the digit store from raw data
  void SetRawDigitStore(AliMUONDigitStoreV1 *ds) { fRawDigitStore = ds; };
  /// set the trigger store from raw data
  void SetRawTriggerStore(AliMUONTriggerStoreV1 *ts) { fRawTriggerStore = ts; };

  /// set the trigger boards manager
  void SetCrateManager(AliMUONTriggerCrateStore *crates) { fCrateManager = crates; };

  void Show();

  void DrawStrips(Bool_t bx, Bool_t by);
  void DrawDigits(Bool_t bx, Bool_t by);
  void DrawClear();
  void EditStrips(Int_t event, Int_t x, Int_t y, TObject *sel);

  void Init();
  void HandleButtons(Int_t id = -1);
  void HandleEditButton();
  void CloseWindow() const;
  void DoClose();
  void DoDigits();
  void ResetDigits();
  void LocalTriggerInfo();

private:
  /// Not implemented
  AliMUONTriggerGUIbdmap (const AliMUONTriggerGUIbdmap& bdmap);
  /// Not implemented
  AliMUONTriggerGUIbdmap& operator=(const AliMUONTriggerGUIbdmap& bdmap);


  enum { kNBoards = 234, kNMT = 4, kNS = 16 };  ///< constants

  TGTransientFrame     *fMain;             ///< Main board frame
  TCanvas              *fCanvas[kNMT];     ///< MT canvases
  TGTextEdit           *fLocTrigE;         ///< Window local trigger info

  AliMUONTriggerGUIboard  *fBoard;           ///< Current board object
  AliLoader               *fLoader;          ///< The MUON loader
  AliMUONMCDataInterface  *fMCDataInterface; ///< MC data interface
  AliMUONDigitStoreV1     *fRawDigitStore;   ///< Raw data digit store
  AliMUONTriggerStoreV1   *fRawTriggerStore; ///< Raw data trigger store

  TGCheckButton        *fXStrips;          ///< Draw x-strips and digits
  TGCheckButton        *fYStrips;          ///< Draw y-strips and digits
  TGCheckButton        *fEditStrips;       ///< Set/unset the strips

  TPolyLine            *fXDigPL[kNMT][kNS];     ///< X-strip polyline 
  TPolyLine            *fYDigPL[kNMT][kNS];     ///< Y-strip polyline
  TBox                 *fXDigBox[kNMT][kNS];    ///< X-digit box
  TBox                 *fYDigBox[kNMT][kNS];    ///< Y-digit box
  TPaveText            *fXLabelL[kNMT][kNS];    ///< X-strip labels left
  TPaveText            *fXLabelR[kNMT][kNS];    ///< X-strip labels right
  TPaveText            *fYLabelL[kNMT][kNS];    ///< Y-strip labels left
  TPaveText            *fYLabelR[kNMT][kNS];    ///< Y-strip labels right

  Float_t               fXWidth[kNMT];     ///< Board x-width
  Float_t               fYWidth[kNMT];     ///< Board y-width
  Float_t               fXCenter[kNMT];    ///< Board x-center
  Float_t               fYCenter[kNMT];    ///< Board y-center
  
  Bool_t                fXOn;              ///< x-strips/digits on canvas ?
  Bool_t                fYOn;              ///< y-strips/digits on canvas ?
  Bool_t                fLabelX;           ///< x-labels exist
  Bool_t                fLabelY;           ///< y-labels exist
  Bool_t                fIsEditable;       ///< allows set/unset the strips

  UInt_t                fCanvasSize;       ///< Size of the canvas
  Int_t                 fNStripX;          ///< Number of x-strips on board
  Int_t                 fNStripY;          ///< Number of y-strips on board

  TObjArray            *fBoards;           ///< Array with all boards

  AliMUONCalibrationData *fCalibrationData;  ///< Pointer to calibration data
  AliMUONTriggerCrateStore *fCrateManager;   ///< trigger boards manager

  ClassDef(AliMUONTriggerGUIbdmap,2)       // board gui class

};

#endif
