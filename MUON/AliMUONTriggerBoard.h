#ifndef ALIMUONTRIGGERBOARD_H
#define ALIMUONTRIGGERBOARD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup sim
/// \class AliMUONTriggerBoard
/// \brief TRIGGER BOARD BASE CLASS
///
//  Author Rachid Guernane (LPCCFd)

#include <TNamed.h>

class AliMUONTriggerBoard : public TNamed
{
   public:
      AliMUONTriggerBoard();
      AliMUONTriggerBoard(const char *name, Int_t islot);
      virtual ~AliMUONTriggerBoard();

      /// \todo add comment 
      virtual void Response() = 0;

      /// \todo add comment 
      virtual void Reset() = 0;

      /// \todo add comment 
      virtual void Scan(Option_t *option) const = 0;

      /// \todo add comment 
      virtual void Resp(Option_t *option) const = 0;

      /// Return response
      virtual UShort_t GetResponse() const {return fResponse;}

      /// \todo add comment 
      virtual void Mask(Int_t index, UShort_t mask) = 0;
      
   protected:
      /// Not implemented
      AliMUONTriggerBoard(const AliMUONTriggerBoard &entry);
      /// Not implemented
      AliMUONTriggerBoard& operator=(const AliMUONTriggerBoard &rhs);

      Int_t fSlot;                ///< SLOT NUMBER IN CRATE

      UShort_t fResponse;         ///< RESPONSE

   ClassDef(AliMUONTriggerBoard,1)
};
#endif

