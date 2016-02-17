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
////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Base Class used to find                                                //
// the reconstructed points for ITS Upgrade                               // 
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include "AliITSUPixelModule.h"


AliITSUPixelModule::AliITSUPixelModule():
  TObject(),
  fCharge(0),
  fModule(0),
  fCol(0),
  fRow(0)
{
  //
  // Default constructor
  //
  for(Int_t i=0;i<kMaxLab;i++)fLabels[i]=-1;
}
//_______________________________________________
AliITSUPixelModule::AliITSUPixelModule( UShort_t mod, UInt_t col, UInt_t row, UInt_t charge, Int_t lab[kMaxLab]):
  TObject(),
  fCharge(charge),
  fModule(mod),
  fCol(col),
  fRow(row)
{
  //
  // Constructor
  //
  for(Int_t i=0;i<kMaxLab;i++)fLabels[i]=lab[i];
}
//____________________________________________
void AliITSUPixelModule::SetLabels(Int_t lab[kMaxLab]){
  // Setter for cluster labels 
  for(Int_t i=0;i<kMaxLab;i++)fLabels[i]=lab[i];
}
//_______________________________________________
void AliITSUPixelModule::PrintInfo(){
  //
  // printout method for debugging
  // 
  printf(" module %d col %i row %i charge %i \n -> labels %i, %i, %i, ", fModule, fCol, fRow,fCharge, fLabels[0], fLabels[1], fLabels[2]);
  for (Int_t i=3;i<kMaxLab;i++) printf("%i, ",fLabels[i]);
  printf("\n");
				       
}



