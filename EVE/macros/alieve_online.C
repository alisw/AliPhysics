/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void alieve_online()
{
  // List of macros to be executed
  gROOT->Macro("its_raw.C");
  gROOT->Macro("its_clusters.C");
}
