//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSISOLATEDMIPTRIGGER_H
#define ALIHLTPHOSISOLATEDMIPTRIGGER_H


/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSMipCluster.h"
#include "AliHLTPHOSConstants.h"
#include "Rtypes.h"
#include "AliHLTPHOSMipClusterManager.h"
#include "TH1.h"

#define MAXFILES 1000

using namespace PhosHLTConst;

class  AliHLTPHOSDigit;

class  AliHLTPHOSIsolatedMipTrigger
{
public:
  AliHLTPHOSIsolatedMipTrigger();
  virtual ~AliHLTPHOSIsolatedMipTrigger();
  void SetFilePath(char *path);
  int  Analyze();
  void SetMipLowCut(Float_t lCut);
  void SetMipHighCut(Float_t hCut); 
  int  ScanFileList(); //scans the file filelist.txt, the files to analyze
  void Init();
  
  void FillClusterHistograms();
 
  TH1F *fClusterHistograms[N_ZROWS_MOD][N_XCOLUMNS_MOD];
  void WriteHistograms(); 


private:
  // char **fFileList;
  char fFileList[MAXFILES][256];
  int fNFiles;
  int fModuleId;

  bool IsMipCandidate(AliHLTPHOSDigit *digit);
  AliHLTPHOSMipClusterManager *fClusterManager;
  //  char fCurrentFilePath[256];
  Float_t fMipLowCut;
  Float_t fMipHighCut;
};

#endif

