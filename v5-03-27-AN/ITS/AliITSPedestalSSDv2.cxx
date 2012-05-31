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


#include "AliITSPedestalSSDv2.h"
#include <cstring>

//////////////////////////////////////////////////////
// Author: Enrico Fragiacomo
// Date: 23/08/2007
// Modified: 08/07/2008
//                                                  //
//////////////////////////////////////////////////////

//const Int_t AliITSPedestalSSD::fgkDefaultNModulesSSD = 1698;
//const Int_t AliITSPedestalSSD::fgkDefaultNStripsSSD = 768;

ClassImp(AliITSPedestalSSDv2)
  
//______________________________________________________________________
  AliITSPedestalSSDv2::AliITSPedestalSSDv2()
    // Default Constructor
      //: fPedestal(new Float_t[2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD]) 
{ 
    for(Int_t i=0; i<2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD; i++) 
      fPedestal[i]=0;    
  }

//______________________________________________________________________
AliITSPedestalSSDv2::AliITSPedestalSSDv2(const AliITSPedestalSSDv2 &source): 
  TObject(source)  
				    //,  fPedestal(new Float_t[2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD])
{
    // copy Constructor
  memcpy(fPedestal,source.fPedestal,
	 2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD*sizeof(Float_t));
}

//______________________________________________________________________
AliITSPedestalSSDv2::~AliITSPedestalSSDv2(){
    // destructor

}

//______________________________________________________________________
AliITSPedestalSSDv2& AliITSPedestalSSDv2::operator=(const AliITSPedestalSSDv2 &source) {
 // ass. op.
    if (this == &source)return *this;

    memcpy(fPedestal,source.fPedestal,
	 2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD*sizeof(Float_t));
    
    return *this;
}
