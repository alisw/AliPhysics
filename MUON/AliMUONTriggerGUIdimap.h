#ifndef ALIMUONTRIGGERGUIDIMAP_H
#define ALIMUONTRIGGERGUIDIMAP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup evaluation
/// \class AliMUONTriggerGUIdimap
/// \brief Trigger GUI utility class: digits maps of the trigger chambers

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Graphical User Interface utility class for the MUON trigger          //
// - digits maps of the trigger chambers                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <RQ_OBJECT.h>

class AliLoader;
class AliMUONData;
class TGTransientFrame;
class TObjArray;
class TRootEmbeddedCanvas;
class TPave;
class TH1F;
class TGWindow;
class TPad;

class AliMUONTriggerGUIdimap : public TObject
{

RQ_OBJECT("AliMUONTriggerGUIdimap")
    
public:

 AliMUONTriggerGUIdimap(AliLoader *loader,
			TObjArray *boards, 
			const TGWindow *p, const TGWindow *main, 
			UInt_t w, UInt_t h);
 virtual ~AliMUONTriggerGUIdimap();
 
 AliMUONTriggerGUIdimap (const AliMUONTriggerGUIdimap& dimap);
 AliMUONTriggerGUIdimap& operator=(const AliMUONTriggerGUIdimap& dimap);

 /// set the current muon loader
 void SetLoader(AliLoader *loader) { fLoader = loader; };
 /// true if the gui is active
 Bool_t IsOn() const { return fIsOn; };

 void DoClose();
 void DoUpdate();
 void DoTab(Int_t id) const;
 void DoReset();
 void CloseWindow();
 void DrawMaps(Int_t chamber);
 void SelectBoard(Int_t ib);

private:

  enum { kNBoards = 234 };        ///< number of boards
  enum { kGood = 0x0001, kWithProblems = 0x0002, kNotWorking = 0x0004, kUnknown = 0x0008 };                      ///< working status flags
  enum { kNSide = 2, kNCol = 7, kNLine = 9, kNMT = 4, kNBoardType = 3 }; ///< other constants
  
  TGTransientFrame    *fMain;     ///< Main frame

  AliLoader   *fLoader;           ///< The MUON loader
  AliMUONData *fMUONData;         ///< MUON data manager

  TRootEmbeddedCanvas *fEc[kNMT]; ///< Canvases for drawing the digits

  TPave *fPaveBoard[kNMT][kNBoards];    ///< Drawing of the board
  TObjArray   *fBoards;           ///< Array of boards

  Bool_t fIsOn;                   ///< True if the map is open

  ClassDef(AliMUONTriggerGUIdimap,1)
   
};

#endif
