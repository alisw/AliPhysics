/**************************************************************************
 * Copyright(c) 2000-2004, ALICE Experiment at CERN, All rights reserved. *
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
 
#include "AliITSRawClusterSSD.h"

////////////////////////////////////////////////////
//  Cluster classes for set:ITS                   //
//  Raw Clusters for SSD                          //
//                                                //
////////////////////////////////////////////////////

ClassImp(AliITSRawClusterSSD)

//______________________________________________________________________
AliITSRawClusterSSD::AliITSRawClusterSSD():
fMultiplicityN(0),
fQErr(0),
fSignalP(0),
fSignalN(0),
fStatus(-1),
fNtracks(0)
{
  // Default constructor
}
//______________________________________________________________________
AliITSRawClusterSSD::AliITSRawClusterSSD(Float_t Prob,Int_t Sp,Int_t Sn):
fMultiplicityN(Sn),
fQErr(0),
fSignalP(0),
fSignalN(0),
fStatus(-1),
fNtracks(0) {  
    // constructor

    Prob = 0.0; // added to remove unused variable warning.
    //fProbability   = Prob;
    fMultiplicity  = Sp;
}
