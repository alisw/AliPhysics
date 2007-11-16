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
      virtual UShort_t Algo(UShort_t i, UShort_t j, char *thres);
                      
                       /// Reset regional board responses
      virtual void     Reset() {for (Int_t i=0; i<16; i++) fRegionalResponse[i] = 0;}

                       /// scan response of regional boards
      virtual void     Scan(Option_t *option) const;

                       /// \todo add comment
      virtual void     Resp(Option_t*) const {}

                      /// Set mask (disable) for regional boards
      void Mask(Int_t index, UShort_t mask);
      
   private:

      UShort_t fRegionalResponse[16]; ///< Regional board responses
      UShort_t fMask[2];              ///< Mask

   ClassDef(AliMUONGlobalTriggerBoard,1) //Global trigger board
};
#endif
