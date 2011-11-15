#ifndef AliOADBPhysicsSelection_H
#define AliOADBPhysicsSelection_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     OADB interface for the physics selection
//     Author: Michele Floris, CERN
//    
// This class contains the parameters needed to configure the physics
// selection. The online trigger class must be associated to an
// offline trigger mask, as defined in
// AliVEvent::EOfflineTriggerTypes.  The offline condition, even for
// the same trigger mask, can be different for every different online
// trigger, so it must be given as a parameter in the Add...Class
// methods.  It's up to the user to set the appropriate offline
// conditions when filling the class.
// -------------------------------------------------------------------------

#include <TNamed.h>
#include "AliVEvent.h"
#include "TList.h"
#include "TMap.h"

class TObjArray;
class TArrayI;


#define NTRIGGERBITS   32
#define NTRIGGERLOGICS 64

class AliOADBPhysicsSelection : public TNamed {

 public :
  AliOADBPhysicsSelection();
  AliOADBPhysicsSelection(const char* name);
  virtual ~AliOADBPhysicsSelection();
  AliOADBPhysicsSelection(const AliOADBPhysicsSelection& cont); 
  AliOADBPhysicsSelection& operator=(const AliOADBPhysicsSelection& cont);
  void Init();
  // Getters
  // These take a mask of trigger bits, because you can activate several triggers at once
  // TList * GetCollTrigClass(AliVEvent::EOfflineTriggerTypes trig) const; 
  // TList * GetBGTrigClass(AliVEvent::EOfflineTriggerTypes trig) const;    
  // TList * GetTrigClass(AliVEvent::EOfflineTriggerTypes trig, TList ** listClasses) const ; 

  TList * GetCollTrigClass(UInt_t triggerBit) { return fCollTrigClasses[triggerBit];}
  TList * GetBGTrigClass(UInt_t triggerBit)   { return fBGTrigClasses[triggerBit];  }
  const TString GetBeamSide (const char * trigger) ;
  TMap * Debug() { return fBeamSide; }
  // Thess take a single trigger bit, as the HW/offline conditions are mapped 1 <-> 1 to a single EOfflineTriggerTypes bit
  const TString  GetHardwareTrigger(UInt_t triggerLogic) const { return triggerLogic >= NTRIGGERLOGICS ? "" : fHardwareTrigger[triggerLogic].String(); }
  const TString  GetOfflineTrigger (UInt_t triggerLogic) const { return triggerLogic >= NTRIGGERLOGICS ? "" : fOfflineTrigger [triggerLogic].String(); }
  UInt_t GetNTriggerBits()  const { return fNtriggerBits; }
  // Setters 
  void AddCollisionTriggerClass(AliVEvent::EOfflineTriggerTypes triggerMask, const char* className,const char * beamSide, UInt_t triggerLogic);
  void AddBGTriggerClass       (AliVEvent::EOfflineTriggerTypes triggerMask, const char* className,const char * beamSide, UInt_t triggerLogic);


  void SetHardwareTrigger      (UInt_t triggerLogic, const char * trigger)  { fHardwareTrigger [triggerLogic].SetString(trigger);     }
  void SetOfflineTrigger       (UInt_t triggerLogic, const char * trigger)  { fOfflineTrigger  [triggerLogic].SetString(trigger);     }
  void SetBeamSide             (const char * className, const char * side);
  // MISC
  virtual Bool_t	IsFolder() const { return kTRUE; }
  void Browse(TBrowser *b);
  virtual void	Print(Option_t* option = "") const;

  static UInt_t GetActiveBit(UInt_t mask) ;

protected:
  void CleanKey(TString & str) ;
  const char* ExpandTriggerString(const char* className);
  
 private :
  
  UInt_t fNtriggerBits; // Size of the arrays below. Initialized using NTRIGGERBITS
  UInt_t fNtriggerLogics; // number of possible trigger logics, initialized using NTRIGGERLOGICS. To be increased if needed.

  TList ** fCollTrigClasses  ; //[fNtriggerBits] Array of collision trigger classes, corresponding to the different trigger bits defined in AliVEvent
  TList ** fBGTrigClasses    ; //[fNtriggerBits] Array of background trigger classes, corresponding to the different trigger bits defined in AliVEvent
  TObjString * fHardwareTrigger ; //[fNtriggerLogics] Array of online trigger condition, corresponding to the different trigger logics set in Add...TriggerClass
  TObjString * fOfflineTrigger  ; //[fNtriggerLogics] Array of offline trigger condition, corresponding to the different trigger logics set in Add...TriggerClass
  TMap * fBeamSide;             // Map from the trigger classname to the beam side ("B", "A", "C", "E", "AC")
  ClassDef(AliOADBPhysicsSelection, 1);
};

#endif
