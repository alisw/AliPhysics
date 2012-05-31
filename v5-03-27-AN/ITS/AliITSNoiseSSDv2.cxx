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


#include "AliITSNoiseSSDv2.h"
#include <cstring>

//////////////////////////////////////////////////////
// Author: Enrico Fragiacomo
// Date: 23/08/2007
// Modified: 08/07/2008
//                                                  //
//////////////////////////////////////////////////////

//const Int_t AliITSNoiseSSD::fgkDefaultNModulesSSD = 1698;
//const Int_t AliITSNoiseSSD::fgkDefaultNStripsSSD = 768;

ClassImp(AliITSNoiseSSDv2)
  
//______________________________________________________________________
  AliITSNoiseSSDv2::AliITSNoiseSSDv2() {
  // Default Constructor
  for(Int_t i=0; i<2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD; i++) 
    fNois[i]=0;    
}

//______________________________________________________________________
AliITSNoiseSSDv2::AliITSNoiseSSDv2(const AliITSNoiseSSDv2 &source): 
  TObject(source)  
{
    // copy Constructor
  memcpy(fNois,source.fNois,
	 2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD*sizeof(Float_t));
}

//______________________________________________________________________
AliITSNoiseSSDv2::~AliITSNoiseSSDv2(){
    // destructor

}

//______________________________________________________________________
AliITSNoiseSSDv2& AliITSNoiseSSDv2::operator=(const AliITSNoiseSSDv2 &source) {
 // ass. op.
    if (this == &source)
      return *this;

    memcpy(fNois,source.fNois,
	 2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD*sizeof(Float_t));
    
    return *this;
}
