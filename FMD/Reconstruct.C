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

/* $Id$ */
/** @file    Reconstruct.C
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 14:19:56 2006
    @brief   Script to do reconstruction 
*/
// Script to do test the FMD digitization class. 

/** Do reconstruction */
void 
Reconstruct()
{
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT");
  AliLog::SetModuleDebugLevel("FMD", 2);
  AliReconstruction rec;   
  AliCDBEntry* align = cdb->Get("FMD/Align/Data");
  if (align) {
    TClonesArray* array = dynamic_cast<TClonesArray*>(align->GetObject());
    if (array) {
      Int_t nAlign = array->GetEntries();
      for (Int_t i = 0; i < nAlign; i++) {
	AliAlignObjAngles* a = static_cast<AliAlignObjAngles*>(array->At(i));
	if (!a->ApplyToGeometry()) {
	  Warning("ApplyAlignement", "Failed to apply alignment to %s", 
		  a->GetVolPath());
	}
      }
    }
  }
  rec.SetRunLocalReconstruction("FMD");
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunTracking(""); 
  rec.SetFillESD("FMD"); 
  rec.SetInput("./");
  rec.Run(); 
}

//
// EOF
//
