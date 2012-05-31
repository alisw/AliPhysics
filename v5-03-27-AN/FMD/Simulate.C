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
Simulate(Int_t n=1)
{
  AliSimulation sim;
  sim.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  sim.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));
  sim.SetConfigFile("./Config.C");
  sim.SetMakeSDigits("FMD");
  sim.SetMakeDigits("FMD"); 
  sim.SetWriteRawData("FMD", "raw.root"); 
  sim.SetRunQA(":");

  AliLog::SetModuleDebugLevel("FMD", 2);

  TStopwatch w; 
  w.Start(); 
  sim.Run(n);  
  w.Stop(); 
  w.Print(); 
}

//
// EOF
//
