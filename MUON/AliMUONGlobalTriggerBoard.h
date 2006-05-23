#ifndef ALIMUONGLOBALTRIGGERBOARD_H
#define ALIMUONGLOBALTRIGGERBOARD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup sim
/// \class AliMUONGlobalTriggerBoard
/// \brief GLOBAL TRIGGER
///
/// \author Rachid Guernane (LPCCFd)

#include "AliMUONTriggerBoard.h"

class AliMUONGlobalTriggerBoard : public AliMUONTriggerBoard
{
   public:

      AliMUONGlobalTriggerBoard();  
      AliMUONGlobalTriggerBoard(const char *name, Int_t a);
      virtual ~AliMUONGlobalTriggerBoard() {;}
    
      virtual void     SetRegionalResponse(UShort_t resp[16]) {for (Int_t i=0; i<16; i++) fRegionalResponse[i] = resp[i];}

      virtual void     Response();

      virtual UShort_t Algo(UShort_t i, UShort_t j, char *thres);

      virtual void     Reset() {for (Int_t i=0; i<16; i++) fRegionalResponse[i] = 0;}

      virtual void     Scan(Option_t *option) const;

      virtual void     Resp(Option_t*) const {}

      void Mask(Int_t index, UShort_t mask);
      
   private:

      UShort_t fRegionalResponse[16]; ///< REGIONAL BOARD RESPONSES
      UShort_t fMask[16];             ///< MASK

   ClassDef(AliMUONGlobalTriggerBoard,1) 
};
#endif
