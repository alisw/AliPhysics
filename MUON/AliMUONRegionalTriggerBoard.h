#ifndef ALIMUONREGIONALTRIGGERBOARD_H
#define ALIMUONREGIONALTRIGGERBOARD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup trigger
/// \class AliMUONRegionalTriggerBoard
/// \brief Regional trigger - real HW algorithm is implemented
///
//  Author: Rachid Guernane (LPCCFd)

#include "AliMUONTriggerBoard.h"

class AliMUONRegionalTriggerBoard : public AliMUONTriggerBoard
{
   public: 
      AliMUONRegionalTriggerBoard();  
      AliMUONRegionalTriggerBoard(const char *name, Int_t a);
      virtual ~AliMUONRegionalTriggerBoard();
    
      /// Reset Local trigger inputs
      virtual void Reset() {for (Int_t i=0; i<16; i++) fLocalResponse[i] = 0;}

      virtual void Scan(Option_t *option) const;

      /// Dummy implementation
      virtual void Resp(Option_t*) const {}

      virtual void Response();

      /// Set Local trigger inputs
      virtual void SetLocalResponse(UShort_t val[16]) {for (Int_t i=0;i<16;i++) fLocalResponse[i] = val[i];}

      /// response of the algorithm
      virtual UShort_t Algo(UShort_t i, UShort_t j, char *thres, Int_t level);

      /// set local boards enable
      void Mask(UShort_t mask);
      
   private:
      UShort_t fLocalResponse[16]; ///< Local trigger inputs
      UShort_t fMask;              ///< Entry mask
      
      ClassDef(AliMUONRegionalTriggerBoard,1) // Regional trigger - real HW algorithm is implemented
};
#endif



