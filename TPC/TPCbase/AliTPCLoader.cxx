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

#include "AliTPCLoader.h"
#include "AliLog.h"

const TString AliTPCLoader::fgkDefaultHitsFileName      = "TPC.Hits.root";
const TString AliTPCLoader::fgkDefaultSDigitsFileName   = "TPC.SDigits.root";
const TString AliTPCLoader::fgkDefaultDigitsFileName    = "TPC.Digits.root";
const TString AliTPCLoader::fgkDefaultRecPointsFileName = "TPC.RecPoints.root";
const TString AliTPCLoader::fgkDefaultTracksFileName    = "TPC.Tracks.root";


ClassImp(AliTPCLoader)
AliTPCLoader::AliTPCLoader()
 {
 }
/*****************************************************************************/ 
AliTPCLoader::AliTPCLoader(const Char_t *name,const Char_t *topfoldername)
 :AliLoader(name,topfoldername)
{
  AliDebug(1,Form("Name = %s; topfolder = %s",name,topfoldername));
}
/*****************************************************************************/ 

AliTPCLoader::AliTPCLoader(const Char_t *name,TFolder *topfolder)
 :AliLoader(name,topfolder)
 {
 }

