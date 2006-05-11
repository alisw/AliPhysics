#ifndef ALIMUONTRIGGERBOARD_H
#define ALIMUONTRIGGERBOARD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup sim
/// \class AliMUONTriggerBoard
/// \brief TRIGGER BOARD BASE CLASS
///
/// \author Rachid Guernane (LPCCFd)

#include <TNamed.h>

class AliMUONTriggerBoard : public TNamed
{
   public:
      AliMUONTriggerBoard();
      AliMUONTriggerBoard(const AliMUONTriggerBoard &entry);
      AliMUONTriggerBoard(const char *name, Int_t islot);
      virtual ~AliMUONTriggerBoard() {}

      AliMUONTriggerBoard& operator=(const AliMUONTriggerBoard &rhs);

      virtual void Response() = 0;

      virtual void Reset() = 0;

      virtual void Scan(Option_t *option) const = 0;

      virtual void Resp(Option_t *option) const = 0;

      virtual UShort_t GetResponse() const {return fResponse;}

      virtual void Mask(Int_t index, UShort_t mask) = 0;
      
   protected:

      Int_t fSlot;                // SLOT NUMBER IN CRATE

      UShort_t fResponse;         // RESPONSE

      void Copy(TObject&) const;

   private:



   ClassDef(AliMUONTriggerBoard,1)
};
#endif

