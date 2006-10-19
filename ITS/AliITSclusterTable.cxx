/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////////
// Class used to simplify some operations with clusters.                 //
// -Function FillArray fills an array wich contains, for each            //
//  ITS module, an array with the indices of all the clusters detected   //
//  by the module. The indices correspond to the cluster indices in class// 
//  AliITSlayer of AliITStrackerV2.                                      //
//  This function is used in AliITStrackerSA::FindTracks.                // 
// -Function FillArrayLabel fills an array wich contains, for each       //
//  particle label, and for each layer, the information on clusters:     //
//  0 if there is no cluster, 1 if there is a cluster with this label.   //
//  This function is used to define trackable tracks.                    //
// -Function FillArrayCoorAngles fills arrays wich contains, for each    //
//  layer, the global coordinates, errors on x,y,z and angles lambda     //
//  and phi for each cluster                                             //
///////////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include "AliITSclusterTable.h"
ClassImp(AliITSclusterTable)

//_______________________________________________________________
AliITSclusterTable::AliITSclusterTable():TObject(),
fOrInd(0),
fX(0),
fY(0),
fZ(0),
fSx(0),
fSy(0),
fSz(0),
fPhi(0),
fLam(0){
  //Default constructor

}
//_______________________________________________________________
AliITSclusterTable::AliITSclusterTable(Float_t x, Float_t y, Float_t z, Float_t sx, Float_t sy, Float_t sz, Double_t phi, Double_t lambda, Int_t index):
fOrInd(index),
fX(x),
fY(y),
fZ(z),
fSx(sx),
fSy(sy),
fSz(sz),
fPhi(phi),
fLam(lambda){
  //Default constructor

}


