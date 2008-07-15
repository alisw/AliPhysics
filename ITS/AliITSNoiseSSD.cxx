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


#include "AliITSNoiseSSD.h"
#include <cstring>

//////////////////////////////////////////////////////
// Author: Enrico Fragiacomo
// Date: 23/08/2007
// Modified: 08/07/2008
//                                                  //
//////////////////////////////////////////////////////

//const Int_t AliITSNoiseSSD::fgkDefaultNModulesSSD = 1698;
//const Int_t AliITSNoiseSSD::fgkDefaultNStripsSSD = 768;

ClassImp(AliITSNoiseSSD)
  
//______________________________________________________________________
  AliITSNoiseSSD::AliITSNoiseSSD() {
  // Default Constructor
  for(Int_t i=0; i<2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD; i++) 
    fNois[i]=0;    
}

//______________________________________________________________________
AliITSNoiseSSD::AliITSNoiseSSD(const AliITSNoiseSSD &source): 
  TObject(source)  
{
    // copy Constructor
  memcpy(fNois,source.fNois,
	 2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD*sizeof(Float_t));
}

//______________________________________________________________________
AliITSNoiseSSD::~AliITSNoiseSSD(){
    // destructor

}

//______________________________________________________________________
AliITSNoiseSSD& AliITSNoiseSSD::operator=(const AliITSNoiseSSD &source) {
 // ass. op.
    if (this == &source)
      return *this;

    memcpy(fNois,source.fNois,
	 2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD*sizeof(Float_t));
    
    return *this;
}
