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
Reconstruct(bool raw=false)
{
  // To reconstruct raw data from FDR-I, please enable below lines: 
  // AliFMDParameters::Instance()->UseRcuTrailer(false);
  // AliFMDParameters::Instance()->UseCompleteHeader(false);

  AliReconstruction rec;   
  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  rec.SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
  rec.SetRunLocalReconstruction("FMD");
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunReconstruction("FMD");
  rec.SetRunTracking(""); 
  rec.SetFillESD("FMD"); 
  rec.SetRunQA(":");
  if (raw) rec.SetInput("raw.root");

  AliLog::SetModuleDebugLevel("FMD", 2);

  TStopwatch timer;
  timer.Start();
  rec.Run();
  timer.Stop();
  timer.Print();

}

//
// EOF
//
