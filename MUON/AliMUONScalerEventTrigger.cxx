/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliMUONScalerEventTrigger.h"

ClassImp(AliMUONScalerEventTrigger)
  const Int_t AliMUONScalerEventTrigger::fgkLocalScalerLength  = 45;
  const Int_t AliMUONScalerEventTrigger::fgkRegScalerLength    = 11;
  const Int_t AliMUONScalerEventTrigger::fgkGlobalScalerLength = 10;
  const Int_t AliMUONScalerEventTrigger::fgkDarcScalerLength   =  6;

//___________________________________________
AliMUONScalerEventTrigger::AliMUONScalerEventTrigger()
 : fLocalL0(0),   
   fLocalHold(0), 
   fLocalClk(0),   
   fLocalLPtNTrig(0), 
   fLocalHPtNTrig(0), 
   fLocalLPtRTrig(0), 
   fLocalHPtRTrig(0), 
   fLocalLPtLTrig(0), 
   fLocalHPtLTrig(0), 
   fLocalLPtSTrig(0), 
   fLocalHPtSTrig(0), 
   fLocalEOS(0),         
   fLocalReset(0),       

   fRegL0(0), 
   fRegClk(0),
   fRegHold(0),      

   fGlobalL0(0), 
   fGlobalClk(0),
   fGlobalHold(0),      
   fGlobalSpare(0),     

   fDarcL0R(0),
   fDarcL0U(0),
   fDarcL0P(0),
   fDarcL0S(0),
   fDarcClk(0),
   fDarcHold(0)
{
  for (Int_t i = 0; i < 8*4; i++)
    fLocalScaler[i] = 0;

  for (Int_t i = 0; i < 8; i++)
    fRegScaler[i] = 0;

  for (Int_t i = 0; i < 6; i++)
    fGlobalScaler[i] = 0;

}
