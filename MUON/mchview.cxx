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

#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUONPainterHelper.h"
#include <Riostream.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TStyle.h>

#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpDDLStore.h"
#include "AliMpManuStore.h"

//______________________________________________________________________________
Int_t Usage()
{
  /// Printout available options of the program
  cout << "mchview " << endl;
  cout << "  --version : shows the current version of the program" << endl;
  cout << "  --use filename.root : reuse a previously saved (from this program) root file. Several --use can be used ;-)" << endl;
  cout << "  --geometry #x#+#+# : manually specify the geometry of the window, ala X11..., e.g. --geometry 1280x900+1600+0 will" << endl;
  cout << "    get a window of size 1280x900, located at (1600,0) from the top-left of the (multihead) display " << endl;
  cout << "  --asciimapping : load mapping from ASCII files instead of OCDB (for debug and experts only...)" << endl;
  cout << "  --de detElemId : start by displaying the given detection element instead of the default view (which is all the chambers)" << endl;
  cout << "  --chamber chamberId (from 1 to 10) : start by displaying the given chamber instead of the default view (which is all the chambers)" << endl;
  cout << "  --ocdb ocdbPath : read the mapping from the given OCDB" << endl;
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
  
  TObjArray filesToOpen;
  Bool_t isGeometryFixed(kFALSE);
  Int_t gix(0),giy(0);
  Int_t gox(0),goy(0);
  Bool_t ASCIImapping(kFALSE);
  TString defaultOCDB("local://$ALICE_ROOT/OCDB");
  
  for ( Int_t i = 0; i <= args.GetLast(); ++i ) 
  {
    TString a(static_cast<TObjString*>(args.At(i))->String());
    if ( a == "--version" ) 
    {
      cout << "mchview Version " << AliMUONMchViewApplication::Version() << " ($Id$)" << endl;
      ++nok;
      return 0;
    }
    if ( a == "--use" && i < args.GetLast() )
    {
      filesToOpen.Add(args.At(i+1));
      ++i;
      nok += 2;
    }
    else if ( a == "--geometry" )
    {
      isGeometryFixed = kTRUE;
      TString g(static_cast<TObjString*>(args.At(i+1))->String());
      sscanf(g.Data(),"%10dx%10d+%10d+%10d",&gix,&giy,&gox,&goy);
      nok += 2;
      ++i;
    }
    else if ( a == "--asciimapping" )
    {
      ++nok;
      ASCIImapping = kTRUE;
    }
    else if ( a == "--de" || a == "--chamber" )
    {
      // do nothing. Let AliMUONMchViewApplication handle that one. (and the next one as well).
      nok += 2;
      ++i;      
    }
    else if ( a == "--ocdb" )
    {
      defaultOCDB = static_cast<TObjString*>(args.At(i+1))->String();
      cout << "Using default storage  = " << defaultOCDB.Data() << endl;
      nok += 2;
      ++i;
    }
    else
    {
      return Usage();
    }
  }
  
  if ( nok < args.GetLast() )
  {
    return Usage();
  }
  
  AliCDBManager::Instance()->SetDefaultStorage(defaultOCDB.Data());
  AliCDBManager::Instance()->SetRun(0);
 
  if ( ASCIImapping ) 
  {
    AliMpDataProcessor mp;
    {
      AliMpDataMap* datamap = mp.CreateDataMap("data");
      AliMpDataStreams dataStreams(datamap);
      AliMpDDLStore::ReadData(dataStreams);
    }
    {
      AliMpDataMap* datamap = mp.CreateDataMap("data_run");
      AliMpDataStreams dataStreams(datamap);
      AliMpManuStore::ReadData(dataStreams);
    }
    
    AliCDBManager::Instance()->SetSpecificStorage("MUON/Calib/Neighbours","local://$ALICE_ROOT/OCDB");

  }
  
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

  AliMUONMchViewApplication* theApp(0x0);

  if ( isGeometryFixed )
  {
    theApp = new AliMUONMchViewApplication("mchview", &argc, argv,gix,giy,gox,goy);
  }
  else
  {
    theApp = new AliMUONMchViewApplication("mchview",&argc,argv);
  }
   
  TIter next(&filesToOpen);
  TObjString* s;
  while ( ( s = static_cast<TObjString*>(next()) ) )
  {
    theApp->Open(s->String().Data());
  }
  
  // --- Start the event loop ---
  theApp->Run(kTRUE);

  delete AliMUONPainterHelper::Instance(); // important to trigger the saving of the env. file
}
