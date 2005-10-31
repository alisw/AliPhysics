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

// Script to do test the FMD digitization class. 
void
Simulate()
{
 AliSimulation sim;
 sim.SetConfigFile("$(ALICE_ROOT)/FMD/Config.C");
 // sim.SetMakeSDigits("FMD");
 sim.SetMakeDigits("FMD");
 sim.SetWriteRawData("FMD");
 // sim.SetMakeDigitsFromHits("FMD");
 sim.Run(1); 
}

//
// EOF
//
