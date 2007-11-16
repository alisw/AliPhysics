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

   protected:
      Int_t fSlot;                ///< SLOT NUMBER IN CRATE

      UShort_t fResponse;         ///< RESPONSE

   private:
      /// Not implemented
      AliMUONTriggerBoard(const AliMUONTriggerBoard &entry);
      /// Not implemented
      AliMUONTriggerBoard& operator=(const AliMUONTriggerBoard &rhs);

   ClassDef(AliMUONTriggerBoard,1) //Trigger board base class
};
#endif

