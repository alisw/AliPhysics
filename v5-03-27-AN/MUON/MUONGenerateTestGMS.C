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

/* $Id: */

/// \ingroup macros
/// \file MUONGenerateTestGMS.C
/// \brief Macro to generate ad hoc GMS alignment matrices in the agreed format
///
/// TClonesArray is saved in the Root file with a key "GMSarray"
/// containing TGeoHMatrix objects with TObject::fUniqueID equal to the geometry
/// module Ids
///
/// \author I. Hrivnacova, IPN Orsay

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpConstants.h"

#include <TFile.h>
#include <TGeoMatrix.h>
#include <TClonesArray.h>

#endif

void MUONGenerateTestGMS(Bool_t print = kFALSE)
{
/// \param print option to switch on printing the processed matrices

  TFile f("data/GMS.root", "RECREATE");
  TClonesArray* array = new TClonesArray("TGeoHMatrix",100);
  
  for (Int_t i=0; i<AliMpConstants::NofGeomModules(); i++) {
    TGeoHMatrix* m = new((*array)[i]) TGeoHMatrix(); 
    m->SetUniqueID(i);
    /// rotate by small angle
    m->RotateX(i*0.01);
    
    if (print) m->Print();
  }
  
  f.WriteObject(array,"GMSarray");
  f.Close();
}  
