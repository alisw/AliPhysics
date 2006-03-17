#ifndef ALIMUONREGIONALTRIGGERRBOARD_H
#define ALIMUONREGIONALTRIGGERRBOARD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*-- Author: Rachid Guernane (LPCCFd)

#include "AliMUONTriggerBoard.h"

class AliMUONRegionalTriggerBoard : public AliMUONTriggerBoard
{
   public: 
      AliMUONRegionalTriggerBoard();  
      AliMUONRegionalTriggerBoard(const char *name, Int_t a);
      virtual ~AliMUONRegionalTriggerBoard() {;}
    
      virtual void Reset() {for (Int_t i=0; i<16; i++) fLocalResponse[i] = 0;}

      virtual void Scan(Option_t *option);

      virtual void Resp(Option_t*) {}

      virtual void Response();

      virtual void SetLocalResponse(UShort_t val[16]) {for (Int_t i=0;i<16;i++) fLocalResponse[i] = val[i];}

      virtual UShort_t Algo(UShort_t i, UShort_t j, char *thres, Int_t level);

      void Mask(Int_t index, UShort_t mask);
      
   private:
      UShort_t fLocalResponse[16];
      UShort_t fMask[16];
      
      ClassDef(AliMUONRegionalTriggerBoard,1)
};
#endif



