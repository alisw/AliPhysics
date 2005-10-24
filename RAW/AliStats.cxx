// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliStats                                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TH1.h>
#include <TFile.h>

#include "AliStats.h"


ClassImp(AliStats)


//______________________________________________________________________________
AliStats::AliStats(const char *filename, Int_t compmode, Bool_t filter)
{
   // Create statistics object.

   fEvents     = 0;
   fFirstRun   = 0;
   fFirstEvent = 0;
   fLastRun    = 0;
   fLastEvent  = 0;
   fChunk      = -0.5;
   fFileName   = filename;
   fCompMode   = compmode;
   fFilter     = filter;
   fRTHist     = 0;
}

//______________________________________________________________________________
AliStats::AliStats(const AliStats &rhs): TObject(rhs)
{
   // AliStats copy constructor.

   operator=(rhs);
}

//______________________________________________________________________________
AliStats::~AliStats()
{
   // Cleanup stats object.

   delete fRTHist;
}

//______________________________________________________________________________
AliStats &AliStats::operator=(const AliStats &rhs)
{
   // AliStats assignment operator.

   if (this != &rhs) {
      TObject::operator=(rhs);
      fEvents     = rhs.fEvents;
      fFirstRun   = rhs.fFirstRun;
      fFirstEvent = rhs.fFirstEvent;
      fLastRun    = rhs.fLastRun;
      fLastEvent  = rhs.fLastEvent;
      fBegin      = rhs.fBegin;
      fEnd        = rhs.fEnd;
      fFileName   = rhs.fFileName;
      fFileSize   = rhs.fFileSize;
      fCompFactor = rhs.fCompFactor;
      fCompMode   = rhs.fCompMode;
      fFilter     = rhs.fFilter;
      fRTHist     = rhs.fRTHist ? (TH1F*) rhs.fRTHist->Clone() : 0;
      fChunk      = rhs.fChunk;
   }
   return *this;
}

//______________________________________________________________________________
void AliStats::Fill(Float_t time)
{
   // Fill histogram. This histogram shows the (hopefully constant) time
   // it takes to fill the ROOT DB.
   // Expects to be called 100 times for each file.

   if (!fRTHist) {
      fRTHist = new TH1F("rtime","Real-time to write data chunk", 100, 0, 100);
      fRTHist->SetDirectory(0);
   }

   fRTHist->Fill(fChunk, time);
   fChunk += 1.0;
}
