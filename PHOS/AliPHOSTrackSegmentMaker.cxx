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
// Algorithm Base class to construct PHOS track segments
// Associates EMC and PPSD clusters
// Unfolds the EMC cluster   
//*-- 
//*-- Author: Dmitri Peressounko (RRC Ki & SUBATECH)


// --- ROOT system ---
#include "TFile.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliPHOSTrackSegmentMaker.h"

ClassImp( AliPHOSTrackSegmentMaker) 


//____________________________________________________________________________
  AliPHOSTrackSegmentMaker:: AliPHOSTrackSegmentMaker() : TTask("","")
{
  // ctor
  fSplitFile= 0 ; 

}

//____________________________________________________________________________
AliPHOSTrackSegmentMaker::AliPHOSTrackSegmentMaker(const char * headerFile, const char * name, const Bool_t toSplit): TTask(name, headerFile)
{
  // ctor
  fSplitFile= 0 ; 
  fToSplit  = toSplit ;
}

//____________________________________________________________________________
  AliPHOSTrackSegmentMaker:: AliPHOSTrackSegmentMaker(const AliPHOSTrackSegmentMaker& ts) : 
    TTask(ts.GetName(), ts.GetTitle())
{
  // ctor
  fSplitFile = new TFile( (ts.fSplitFile)->GetName(), "new") ; 
  fToSplit   = ts.fToSplit ; 
}

//____________________________________________________________________________
AliPHOSTrackSegmentMaker::~AliPHOSTrackSegmentMaker()
{
   
      fSplitFile = 0 ;
}

