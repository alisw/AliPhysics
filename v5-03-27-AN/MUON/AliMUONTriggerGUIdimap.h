#ifndef ALIMUONTRIGGERGUIDIMAP_H
#define ALIMUONTRIGGERGUIDIMAP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup evaluation
/// \class AliMUONTriggerGUIdimap
/// \brief Trigger GUI utility class: digits maps of the trigger chambers
//  Author Bogdan Vulpescu, LPC Clermont-Ferrand

#include <TGFrame.h>

class AliLoader;
class AliMUONMCDataInterface;
class AliMUONDigitStoreV1;
class TGTransientFrame;
class TObjArray;
class TRootEmbeddedCanvas;
class TPave;
class TH1F;
class TGWindow;
class TPad;

class AliMUONTriggerGUIdimap : public TGFrame
{
    
public:

 AliMUONTriggerGUIdimap(TObjArray *boards, 
			const TGWindow *p, const TGWindow *main, 
			UInt_t w, UInt_t h);
 virtual ~AliMUONTriggerGUIdimap();
 
 /// set the current muon loader
 void SetLoader(AliLoader * const loader) { fLoader = loader; };
 /// set the MC data interface
 void SetMCDataInterface(AliMUONMCDataInterface * const mc) { fMCDataInterface = mc; };
 /// set the digit store from raw data
 void SetRawDigitStore(AliMUONDigitStoreV1 * const ds) { fRawDigitStore = ds; };

 /// return info if the map is open
 Bool_t IsOn() const { return fIsOn; };
 void DoClose();
 void DoUpdate();
 void DoTab(Int_t id) const;
 void DoReset();
 void CloseWindow();
 void DrawMaps(Int_t chamber);
 void SelectBoard(Int_t ib);
 void DrawAllMaps();

private:
 /// Not implemented  
 AliMUONTriggerGUIdimap (const AliMUONTriggerGUIdimap& dimap);
 /// Not implemented  
 AliMUONTriggerGUIdimap& operator=(const AliMUONTriggerGUIdimap& dimap);
 
private:

  enum { kNBoards = 234 };        ///< number of boards
  enum { kGood = 0x0001, kWithProblems = 0x0002, kNotWorking = 0x0004, kUnknown = 0x0008 };                      ///< working status flags
  enum { kNSide = 2, kNCol = 7, kNLine = 9, kNMT = 4, kNBoardType = 3 }; ///< other constants
  
  TGTransientFrame    *fMain;     ///< Main frame

  AliLoader   *fLoader;           ///< The MUON loader
  AliMUONMCDataInterface *fMCDataInterface;  ///< MC data interface
  AliMUONDigitStoreV1 *fRawDigitStore;       ///< Raw data digit store

  TRootEmbeddedCanvas *fEc[kNMT]; ///< Canvases for drawing the digits

  TPave *fPaveBoard[kNMT][kNBoards];    ///< Drawing of the board
  TObjArray   *fBoards;           ///< Array of boards

  Bool_t fIsOn;                   ///< True if the map is open

  ClassDef(AliMUONTriggerGUIdimap,2) //Trigger GUI utility class: digits maps of the trigger chambers
   
};

#endif
