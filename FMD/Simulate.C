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
/** @file    Simulate.C
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 14:20:24 2006
    @brief   Script to do simulation 
*/
/** Script to do test the FMD digitization class.  
 */
void
Simulate()
{
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT");
  AliSimulation sim;
  AliCDBEntry* align = cdb->Get("FMD/Align/Data");
  if (align) {
    TClonesArray* array = dynamic_cast<TClonesArray*>(align->GetObject());
    if (array) sim.SetAlignObjArray(array);
  }
  AliLog::SetModuleDebugLevel("FMD", 2);
  sim.SetConfigFile("$(ALICE_ROOT)/FMD/Config.C");
  // sim.SetMakeSDigits("FMD");
  sim.SetMakeDigits("FMD"); 
  sim.SetWriteRawData("FMD"); 
  // sim.SetMakeDigitsFromHits("FMD"); 
  TStopwatch w; 
  w.Start(); 
  sim.Run(1);  
  w.Stop(); 
  w.Print(); 
}

//
// EOF
//
