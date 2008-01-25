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

/// \ingroup graphics
/// \file mchview.cxx
/// \brief Tracker visualization program
///
/// \author Laurent Aphecetche, Subatech


#include "AliMUONMchViewApplication.h"
#include "AliMUONPainterHelper.h"
#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliLog.h"
#include <TROOT.h>
#include <TStyle.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <Riostream.h>

//______________________________________________________________________________
Int_t Usage()
{
  /// Printout available options of the program
  cout << "mchview " << endl;
  cout << "  --version : shows the current version of the program" << endl;
  return -1;
}

//______________________________________________________________________________
int main(int argc, char** argv)
{
  /// Main function for the program
  TObjArray args;
  
  for ( int i = 1; i < argc; ++i ) 
  {
    args.Add(new TObjString(argv[i]));
  }
  
  Int_t nok(0);
  
  for ( Int_t i = 0; i <= args.GetLast(); ++i ) 
  {
    TString a(static_cast<TObjString*>(args.At(i))->String());
    if ( a == "--version" ) 
    {
      cout << "mchview Version " << AliMUONMchViewApplication::Version() << " ($Id$)" << endl;
      ++nok;
      return 0;
    }
  }
  
  if ( nok < args.GetLast() )
  {
    return Usage();
  }
  
  AliWarningGeneral("main","Remove default storage and run number from here...");
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetRun(0);
 
  gROOT->SetStyle("Plain");  
  gStyle->SetPalette(1);
  Int_t n = gStyle->GetNumberOfColors();
  Int_t* colors = new Int_t[n+2];
  for ( Int_t i = 1; i <= n; ++i )
  {
    colors[i] = gStyle->GetColorPalette(i-1);
  }
  colors[0] = 0;
  colors[n+1] = 1;
  gStyle->SetPalette(n+2,colors);
  delete[] colors;
  
  TRint* theApp = new AliMUONMchViewApplication("mchview", &argc, argv, 0.7, 0.9);
   
  AliCodeTimer::Instance()->Print();

  // --- Start the event loop ---
  theApp->Run(kTRUE);

  AliMUONPainterHelper::Instance()->Save();
}
