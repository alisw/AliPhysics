#ifndef ALIMUONPAINTERINTERFACEHELPER_H
#define ALIMUONPAINTERINTERFACEHELPER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterInterfaceHelper
/// \brief Helper class to ease building a GUI with button groups...
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class TGWindow;
class TGButton;
class TGButtonGroup;
class TString;

class AliMUONPainterInterfaceHelper : public TObject
{
public:
  AliMUONPainterInterfaceHelper();
  virtual ~AliMUONPainterInterfaceHelper();

  static void AddRadioButton(TGButtonGroup& bg,
                             const TString& name,
                             void* userData=0x0,
                             Bool_t select=kFALSE);

  static void AddCheckButton(TGButtonGroup& bg,
                             const TString& name,
                             void* userData=0x0,
                             Bool_t select=kFALSE);
  
  /// Id of first button in a group
  static Int_t ButtonStartingId() { return 1; }

  static void ClearButtons(TGButtonGroup& bg);
  
  using TObject::Copy;
  
  static void Copy(const TGButtonGroup& src, TGButtonGroup& dest);
    
  using TObject::Dump;
  
  static void Dump(const TGButtonGroup& bg);
    
  static TGButton* FindButtonByName(const TGButtonGroup& bg, const TString& name);

  static TGButton* FindButtonByUserData(const TGButtonGroup& bg, const void* userData);

  static TGButton* FindDownButton(const TGButtonGroup& bg);
  
  static void SetBackgroundColor(const char* resourceBaseName, TGWindow& window);

  static void SetState(TGButtonGroup& bg, Bool_t state);

  static void Select(TGButtonGroup& bg, const TString& buttonName, Bool_t emit=kFALSE);
  
  static void Unselect(TGButtonGroup& bg, const TString& buttonName, Bool_t emit=kFALSE);
    
  static void RemoveButton(TGButtonGroup& bg, const TGButton* button);
  
  ClassDef(AliMUONPainterInterfaceHelper,1) // Helper class for TGButtonGroup manipulation
};

#endif
