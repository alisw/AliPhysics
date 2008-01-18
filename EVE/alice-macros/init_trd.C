// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/
void init_trd()
{
  TString macdir("$(REVESYS)/alice-macros");
  gSystem->ExpandPathName(macdir);
  gROOT->GetListOfBrowsables()->Add
    (new TSystemDirectory(macdir.Data(), macdir.Data()));
  TEveUtil::AssertMacro("region_marker.C");

  AliEveTRDLoaderManager *trd=new AliEveTRDLoaderManager("TRD manager", "Loader manager for TRD data monitoring");
  gEve->AddElement(trd);
  gEve->AddToListTree(trd, kTRUE);
}
