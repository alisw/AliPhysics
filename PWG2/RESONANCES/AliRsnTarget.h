#ifndef ALIRSNTARGET_H
#define ALIRSNTARGET_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
////////////////////////////////////////////////////////////////////////////////
//
//  General class for cut/computation target objects (daughter/pair/event).
//
////////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

#include "AliRsnEvent.h"

class AliRsnDaughter;
class AliRsnMother;

class AliRsnTarget : public TNamed {
public:

   enum ETargetType {
      kDaughter,
      kMother,
      kEvent,
      kTargetTypes
   };

   AliRsnTarget() : fTargetType(kTargetTypes), fDaughter(0x0), fMother(0x0), fEvent(0x0) { /*nothing*/ }
   AliRsnTarget(const char *name, ETargetType type = kTargetTypes) : TNamed(name, ""), fTargetType(type), fDaughter(0x0), fMother(0x0), fEvent(0x0) { /*nothing*/ }
   AliRsnTarget(const AliRsnTarget& copy) : TNamed(copy), fTargetType(copy.fTargetType), fDaughter(0x0), fMother(0x0), fEvent(0x0) { /*nothing*/ }
   AliRsnTarget& operator=(const AliRsnTarget& copy) { TNamed::operator=(copy); fTargetType = copy.fTargetType; return (*this); }
   virtual ~AliRsnTarget() { /*nothing*/ }

   Bool_t           IsAllNull()                       {return (!fDaughter && !fMother && !fEvent);}
   Bool_t           IsTarget(ETargetType targetType)  {return (fTargetType == targetType);}
   Bool_t           TargetOK(TObject *object);
   
   void             SetTargetType(ETargetType type)   {fTargetType = type;}
   ETargetType      GetTargetType() const             {return fTargetType;}
   Char_t           GetTargetTypeChar() const;
   const char*      GetTargetTypeName() const;
   
   AliRsnDaughter*  GetTargetDaughter()               {return fDaughter;}
   AliRsnMother*    GetTargetMother()                 {return fMother;}
   AliRsnEvent*     GetTargetEvent()                  {return fEvent;}

protected:

   ETargetType     fTargetType;  //  target type selected for this object
   AliRsnDaughter *fDaughter;    //! internal pointer to which any checked object is cast if it matches expected type
   AliRsnMother   *fMother;      //! internal pointer to which any checked object is cast if it matches expected type
   AliRsnEvent    *fEvent;       //! internal pointer to which any checked object is cast if it matches expected type

   // ROOT dictionary
   ClassDef(AliRsnTarget, 2)
};

typedef AliRsnTarget::ETargetType RSNTARGET;

#endif
