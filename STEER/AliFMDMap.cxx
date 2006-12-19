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

//____________________________________________________________________
//                                                                          
// Base class for caches of per-strip information.
// This is used to index a strip. 
// Data stored depends on derived class. 
// This class provides some common infra-structure.
// Derived classes sould define Reset, and operator(). 
//
#include "AliFMDMap.h"		// ALIFMDMAP_H
#include "AliLog.h"
//#include <TClass.h>
//#include <TBuffer.h>
#include <TFile.h>
#include <TList.h>
#include <TStreamerInfo.h>

//____________________________________________________________________
ClassImp(AliFMDMap)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDMap::AliFMDMap(UShort_t maxDet, 
		     UShort_t maxRing, 
		     UShort_t maxSec, 
		     UShort_t maxStr)
  : fMaxDetectors(maxDet), 
    fMaxRings(maxRing), 
    fMaxSectors(maxSec), 
    fMaxStrips(maxStr)
{
  // Construct a map
  //
  // Parameters:
  //     maxDet       Maximum # of detectors
  //     maxRinf      Maximum # of rings
  //     maxSec       Maximum # of sectors
  //     maxStr       Maximum # of strips
  SetBit(kNeedUShort, kFALSE);
}

//____________________________________________________________________
void
AliFMDMap::CheckNeedUShort(TFile* file) 
{
  if (!file) return;
  TObject* o = file->GetStreamerInfoList()->FindObject(ClassName());
  if (!o) return;
  TStreamerInfo* info = static_cast<TStreamerInfo*>(o);
  if (info->GetClassVersion() == 2) SetBit(kNeedUShort);
}
//____________________________________________________________________
Int_t 
AliFMDMap::CheckIndex(UShort_t det, Char_t ring, UShort_t sec, UShort_t str) const
{
  // Check that the index supplied is OK.   Returns true index, or -1
  // on error. 
  if (det < 1) return -1;
  UShort_t ringi = (ring == 'I' ||  ring == 'i' ? 0 : 1);
  Int_t idx = 
    (str + fMaxStrips * (sec + fMaxSectors * (ringi + fMaxRings * (det-1))));
  if (TestBit(kNeedUShort)) idx = UShort_t(idx);
  if (idx < 0 || idx >= fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips) 
    return -1;
  return idx;
}

    
//____________________________________________________________________
Int_t 
AliFMDMap::CalcIndex(UShort_t det, Char_t ring, UShort_t sec, UShort_t str) const
{
  // Calculate index into storage from arguments. 
  // 
  // Parameters: 
  //     det       Detector #
  //     ring      Ring ID
  //     sec       Sector # 
  //     str       Strip # 
  //
  // Returns appropriate index into storage 
  //
  Int_t idx = CheckIndex(det, ring, sec, str);
  if (idx < 0) {
    UShort_t ringi = (ring == 'I' ||  ring == 'i' ? 0 : 1);
    AliFatal(Form("Index (%d,'%c',%d,%d) out of bounds, "
		  "in particular the %s index ", 
		  det, ring, sec, str, 
		  (det > fMaxDetectors ? "Detector" : 
		   (ringi >= fMaxRings ? "Ring" : 
		    (sec >= fMaxSectors ? "Sector" : "Strip")))));
    return 0;
  }
  return idx;
}

#if 0
//___________________________________________________________________
void AliFMDMap::Streamer(TBuffer &R__b)
{
  // Stream an object of class AliFMDMap.
  // This is overridden so that we can know the version of the object
  // that we are reading in.  In this way, we can fix problems that
  // might occur in the class. 
  if (R__b.IsReading()) {
    // read the class version from the buffer
    UInt_t R__s, R__c;
    Version_t version = R__b.ReadVersion(&R__s, &R__c, this->Class());
    TFile *file = (TFile*)R__b.GetParent();
    if (file && file->GetVersion() < 30000) version = -1; 
    AliFMDMap::Class()->ReadBuffer(R__b, this, version, R__s, R__c);
    if (version == 2) SetBit(kNeedUShort);
  } else {
    AliFMDMap::Class()->WriteBuffer(R__b, this);
  }
}
#endif

//___________________________________________________________________
//
// EOF
//
