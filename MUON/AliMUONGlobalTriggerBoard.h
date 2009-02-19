#ifndef ALIMUONGLOBALTRIGGERBOARD_H
#define ALIMUONGLOBALTRIGGERBOARD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup trigger
/// \class AliMUONGlobalTriggerBoard
/// \brief Global trigger board
///
//  Author: Rachid Guernane (LPCCFd)

#include "AliMUONTriggerBoard.h"

class AliMUONGlobalTriggerBoard : public AliMUONTriggerBoard
{
   public:

      AliMUONGlobalTriggerBoard();  
      AliMUONGlobalTriggerBoard(const char *name, Int_t a);
      virtual ~AliMUONGlobalTriggerBoard();
                       
                       /// Set regional board responses
      virtual void     SetRegionalResponse(UShort_t resp[16]) {for (Int_t i=0; i<16; i++) fRegionalResponse[i] = resp[i];}

      virtual void     Response();

                       /// response of the algorithm
      virtual UShort_t Algo(UShort_t i, UShort_t j, const char *thres);
                      
                       /// Reset regional board responses
      virtual void     Reset() {for (Int_t i=0; i<16; i++) fRegionalResponse[i] = 0;}

                       /// scan response of regional boards
      virtual void     Scan(Option_t *option) const;

                       /// Dummy implementation
      virtual void     Resp(Option_t*) const {}

                       /// Set mask for global input (from regional boards)
      void             Mask(Int_t index, UInt_t mask);

                       /// Build the 4 words (32bits) global input
      void             BuildGlobalInput();
                       /// Apply masks to global input
      void             MaskGlobalInput();
                       /// Global input 4 words (32bits) from regional responses
      UInt_t*          GetGlobalInput() { return fGlobalInput; }
      
   private:

      UShort_t fRegionalResponse[16]; ///< Regional board responses
      UInt_t   fGlobalInput[4];       ///< Global input 
      UInt_t   fMask[4];              ///< Mask for the global input

   ClassDef(AliMUONGlobalTriggerBoard,2) //Global trigger board
};
#endif
