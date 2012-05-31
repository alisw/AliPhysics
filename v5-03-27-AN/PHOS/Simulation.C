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
//_________________________________________________________________________
// Macros performing the full simulation chain 
// Use Case : 
//          root> .L Simulation.C++
//          root> sim(2, "GSD", "PHOS EMCAL") --> generates 2 events (Config.C)
//                                                and produces Digits and SDigits  
//                                                for PHOS and EMCAL 
//          root> sim(2, "GSDM", "PHOS EMCAL") --> same as before but before making Digits
//                                                 the 2 signal events are merged with one 
//                                                 background event in ./bgrd/galice.root
// author  : Yves Schutz (CERN/SUBATECH)
// February 2004
//_________________________________________________________________________
#include "AliSimulation.h"
#include "TString.h"
#include "Riostream.h"
#include "TStopwatch.h"

void simu(Int_t nevents=1, TString opt="GSD", TString name="all", TString cdb="") 
{
  AliSimulation sim ; 
  // Generation and simulation
  if ( !opt.Contains("G") )
    sim.SetRunGeneration(kFALSE) ;
  // Making SDigits 
  if ( !opt.Contains("S") )
    sim.SetMakeSDigits("") ; 
  else {
   sim.SetMakeSDigits(name.Data()) ;
    if (name.Contains("ITS") || name.Contains("all"))
     sim.SetMakeDigitsFromHits("ITS ") ; 
    if (name.Contains("TPC") || name.Contains("all"))
     sim.SetMakeDigitsFromHits("TPC") ;
  }  
  // Making Digits 
  if ( !opt.Contains("D") )
    sim.SetMakeDigits("") ; 
  else 
    sim.SetMakeDigits(name.Data()) ;    
  //Merging
  if ( opt.Contains("M") )
    sim.MergeWith("bgrd/galice.root") ;  
  //Region of Interest
  sim.SetRegionOfInterest(kFALSE) ; 
  if ( opt.Contains("W") )
     sim.SetWriteRawData(name.Data(),"raw.root",kTRUE) ;
  if ( !cdb.IsNull() ) 
     sim.SetDefaultStorage(cdb.Data()) ; 

  TStopwatch timer;
  timer.Start();
  sim.Run(nevents) ;    
  timer.Stop();
  timer.Print();

}
