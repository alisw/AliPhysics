/**************************************************************************
 * Copyright(c) 2008-2010, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

/////////////////////////////////////////////////////////////////
// Author: A. Mastroserio                                      //
// This class is used as a container for the FO scan           //
//                                                             //
/////////////////////////////////////////////////////////////////


#include <TObjArray.h>
#include "AliITSOnlineSPDfoChipConfig.h"
#include "AliITSOnlineSPDfoChip.h"

ClassImp(AliITSOnlineSPDfoChip)
//-----------------------------------------------
AliITSOnlineSPDfoChip::AliITSOnlineSPDfoChip():
TObject(),
fActiveHS(0),
fChipId(0),
fNumDACparams(0),
fNumChipConfigs(0),
fDACparams(0x0),
fChipConfigArray(0x0)
{
// default constructor
} 

//-------------------------------------------------
AliITSOnlineSPDfoChip::AliITSOnlineSPDfoChip(Short_t nparams):
TObject(),
fActiveHS(0),
fChipId(0),
fNumDACparams(nparams),
fNumChipConfigs(0),
fDACparams(0x0),
fChipConfigArray(0x0)
{
 //
 // constructor 
 //
  fDACparams = new Short_t[fNumDACparams];
  fChipConfigArray = new TObjArray(); //starting with no objects, then adding...
}
//--------------------------------------------------

AliITSOnlineSPDfoChip::AliITSOnlineSPDfoChip(const AliITSOnlineSPDfoChip &c):
TObject(c),
fActiveHS(c.fActiveHS),
fChipId(c.fChipId),
fNumDACparams(c.fNumDACparams),
fNumChipConfigs(c.fNumChipConfigs),
fDACparams(0x0),
fChipConfigArray(0x0)
{
  //
  //copy constructor 
  // fChipConfigArray is not copied. This method is private
  //

  for(Int_t iPar =0; iPar < fNumDACparams; iPar++) fDACparams[iPar] = c.fDACparams[iPar]; 
}    
//--------------------------------------------------
AliITSOnlineSPDfoChip::~AliITSOnlineSPDfoChip()
{
 // dctor
  
   fChipConfigArray->SetOwner(kTRUE);
   fChipConfigArray->Clear();
   delete fChipConfigArray;
   delete [] fDACparams;
}
//--------------------------------------------------
void AliITSOnlineSPDfoChip::AddMeasurement(AliITSOnlineSPDfoChipConfig *confinfo)
{
  //
  //setter for the Chip global measurement
  //
  
  fChipConfigArray->AddLast(confinfo);
  fNumChipConfigs++;
}
//_____________________________________________________
void AliITSOnlineSPDfoChip::PrintInfo()
{
  //
  // prints container content
  //
  printf(" \n ActiveHS %d   ChipId %d",fActiveHS,fChipId);
  printf(" DAC parameters : %d    %d    %d    %d \n",fDACparams[0],fDACparams[1],fDACparams[2],fDACparams[3]); 
  for(Int_t i=0; i< fNumChipConfigs; i++) {
    if(!(AliITSOnlineSPDfoChipConfig*)(GetChipConfigInfo()->At(i)) ) printf("AliISOnlineSPDfoChipConfig pointer null");
    else ((AliITSOnlineSPDfoChipConfig*)(GetChipConfigInfo()->At(i)))->PrintInfo();   
  }
}
