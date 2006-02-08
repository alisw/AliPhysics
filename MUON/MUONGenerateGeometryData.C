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

// $Id$
//
// Macro for generating the geometry data files:
// (volpaths.dat, transform.dat, svmap.dat).
// To be run from aliroot:
// .x MUONGenerateGeometryData.C
//
// The generated files do not replace the existing ones
// but have different names (with extension ".out").
//
//  Author: I. Hrivnacova, IPN Orsay

void MUONGenerateGeometryData(Bool_t volpaths = true,
                              Bool_t transforms = true, 
                              Bool_t svmaps = true)
{
  // Initialize
  gAlice->Init("$ALICE_ROOT/MUON/Config.C");
  cout << "Init done " << endl;

  // Get MUON detector
  AliMUON* muon = (AliMUON*)gAlice->GetModule("MUON");
  if (!muon) {
    AliFatal("MUON detector not defined.");
    return 0;
  }  

  // Get geometry builder
  AliMUONGeometryBuilder* builder = muon ->GetGeometryBuilder();
  
  if (volpaths) {
    cout << "Generating volpath file ..." << endl;
    builder->GetTransformer()->WriteVolumePaths("volpath.dat.out");
  }  

  if (transforms) {
    cout << "Generating transformation file ..." << endl;
    builder->GetTransformer()->WriteTransformations("transform.dat.out");
  }  

  if (svmaps) {
    cout << "Generating svmaps file ..." << endl;
    builder->WriteSVMaps();
  }  
}  
