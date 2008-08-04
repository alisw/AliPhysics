/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                         *
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


#include "AliITSGainSSDv2.h"
#include <cstring>

//////////////////////////////////////////////////////
// Author: Enrico Fragiacomo
// Date: 23/08/2007
// Modified: 08/07/2008
//                                                  //
//////////////////////////////////////////////////////

//const Int_t AliITSGainSSD::fgkDefaultNModulesSSD = 1698;
//const Int_t AliITSGainSSD::fgkDefaultNStripsSSD = 768;

ClassImp(AliITSGainSSDv2)
  
//______________________________________________________________________
  AliITSGainSSDv2::AliITSGainSSDv2() 
{
    // Default Constructor  
    for(Int_t i=0; i<2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD; i++) 
      fGain[i]=0;    
}

//______________________________________________________________________
AliITSGainSSDv2::AliITSGainSSDv2(const AliITSGainSSDv2 &source): 
  TObject(source)
{
    // copy Constructor
  memcpy(fGain,source.fGain,
	 2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD*sizeof(UShort_t));
}

//______________________________________________________________________
AliITSGainSSDv2::~AliITSGainSSDv2(){
    // destructor

}

//______________________________________________________________________
AliITSGainSSDv2& AliITSGainSSDv2::operator=(const AliITSGainSSDv2 &source) {
 // ass. op.
    if (this == &source)return *this;

    memcpy(fGain,source.fGain,
	 2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD*sizeof(UShort_t));
    
    return *this;
}
