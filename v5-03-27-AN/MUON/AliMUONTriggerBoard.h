#ifndef ALIMUONTRIGGERBOARD_H
#define ALIMUONTRIGGERBOARD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup trigger
/// \class AliMUONTriggerBoard
/// \brief Trigger board base class
///
//  Author Rachid Guernane (LPCCFd)

#include <TNamed.h>

class AliMUONTriggerBoard : public TNamed
{
   public:
      AliMUONTriggerBoard();
      AliMUONTriggerBoard(const char *name, Int_t islot);
      virtual ~AliMUONTriggerBoard();

      /// virtual method for derivated classes
      virtual void Response() = 0;

      /// virtual method for derivated classes
      virtual void Reset() = 0;

      /// virtual method for derivated classes
      virtual void Scan(Option_t *option) const = 0;

      /// virtual method for derivated classes
      virtual void Resp(Option_t *option) const = 0;

      /// Return response
      virtual UShort_t GetResponse() const {return fResponse;}

      AliMUONTriggerBoard(const AliMUONTriggerBoard &rhs);
      AliMUONTriggerBoard& operator=(const AliMUONTriggerBoard &rhs);

   protected:
      Int_t fSlot;                ///< SLOT NUMBER IN CRATE

      UShort_t fResponse;         ///< RESPONSE

   ClassDef(AliMUONTriggerBoard,1) //Trigger board base class
};
#endif

