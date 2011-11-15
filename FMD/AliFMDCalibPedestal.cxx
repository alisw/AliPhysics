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
/** @file    AliFMDCalibPedestal.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:30:36 2006
    @brief   Per strip pedestal calibration 
    @ingroup FMD_base
*/
//____________________________________________________________________
//                                                                          
// This class stores a pedestal and pedestal width for each strip in
// the FMD detectors. 
// The values are stored as floats, since they may be results from a
// fit. 
// Need to make algorithm that makes this data
//
#include "AliFMDCalibPedestal.h"	// ALIFMDCALIBPEDESTAL_H
#include <iostream>
#include <TString.h>
#include <AliLog.h>
#include "AliFMDDebug.h"
#include "AliFMDBoolMap.h"

//____________________________________________________________________
ClassImp(AliFMDCalibPedestal)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDCalibPedestal::AliFMDCalibPedestal()
  : fValue(0), // nDet == 0 mean 51200 entries 
    fWidth(0)  // nDet == 0 mean 51200 entries
{
  // CTOR 
  fValue.Reset(-1.);
  fWidth.Reset(-1.);
}

//____________________________________________________________________
AliFMDCalibPedestal::AliFMDCalibPedestal(const AliFMDCalibPedestal& o)
  : TObject(o), 
    fValue(o.fValue), 
    fWidth(o.fWidth)
{
  // Copy Ctor 
}

//____________________________________________________________________
AliFMDCalibPedestal&
AliFMDCalibPedestal::operator=(const AliFMDCalibPedestal& o)
{
  // Assignment operator 
  if (&o == this) return *this; 
  fValue = o.fValue;
  fWidth = o.fWidth;
  return (*this);
}

//____________________________________________________________________
void
AliFMDCalibPedestal::Set(UShort_t det, Char_t ring, UShort_t sec, 
			 UShort_t str, Float_t ped, Float_t pedW)
{
  // set value and width for a strip 
  if (fValue.CheckIndex(det, ring, sec, str) < 0) return;
  fValue(det, ring, sec, str) = ped;
  fWidth(det, ring, sec, str) = pedW;
}

//____________________________________________________________________
Float_t
AliFMDCalibPedestal::Value(UShort_t det, Char_t ring, UShort_t sec, 
			   UShort_t str)
{
  // Get pedestal value for a strip 
  return fValue(det, ring, sec, str);
}

//____________________________________________________________________
Float_t
AliFMDCalibPedestal::Width(UShort_t det, Char_t ring, UShort_t sec, 
			   UShort_t str)
{
  // Get pedestal width for a strip 
  return fWidth(det, ring, sec, str);
}

//____________________________________________________________________
namespace {
  struct MakeDead : public AliFMDMap::ForOne
  {
    MakeDead(AliFMDBoolMap* dead, Float_t max) 
      : fDead(dead), fMax(max), fCount(0)
    {}
    MakeDead(const MakeDead& other) 
      : AliFMDMap::ForOne(other),
        fDead(other.fDead), fMax(other.fMax), fCount(other.fCount)
    {}
    MakeDead& operator=(const MakeDead& other) 
    { 
      if (&other == this) return *this; 
      fDead   = other.fDead;
      fMax    = other.fMax;
      fCount  = other.fCount;
      return *this;
    }
    Bool_t operator()(UShort_t d, Char_t r, UShort_t s, UShort_t t, Float_t v)
    {
      AliDebugGeneral("AliFMDCalibPedestal::MakeDeadMap", 100, 
		      Form("Checking if noise of FMD%d%c[%2d,%3d]=%f "
			      "is larger than %f", d, r, s, t, v, fMax));
      if (v > fMax) {
	fDead->operator()(d,r,s,t) = kTRUE;
	fCount++;
      }
      return kTRUE;
    }
    Bool_t operator()(UShort_t, Char_t, UShort_t, UShort_t, Int_t) 
    { return kFALSE; }
    Bool_t operator()(UShort_t, Char_t, UShort_t, UShort_t, UShort_t)
    { return kFALSE; }
    Bool_t operator()(UShort_t, Char_t, UShort_t, UShort_t, Bool_t)
    { return kFALSE; }
    AliFMDBoolMap* fDead;
    Float_t        fMax;
    Int_t          fCount;
  };
}

//____________________________________________________________________
AliFMDBoolMap*
AliFMDCalibPedestal::MakeDeadMap(Float_t maxW, AliFMDBoolMap* dead) const
{
  // 
  // Make a dead map based on the noise of the channels.  If the noise
  // of a paraticular channel is larger than @a maxW, then the channel
  // is marked as dead. 
  //
  // If the argument @a dead is non-null, then the map passed is
  // modified.  That is, channels marked as dead in the map will
  // remain marked.   Channels that meat the criterion (noise larger
  // than @a maxW) will in addition be marked as dead. 
  //
  // If the argument @a dead is null, then a new map is created and a
  // pointer to this will be returned. 
  // 
  // Parameters:
  //    maxW Maximum value of noise for a channel before it is
  // marked as dead. 
  //    dead If non-null, then modify this map. 
  // 
  // Return:
  //    A pointer to possibly newly allocated dead map. 
  //
 if (!dead) { 
    dead = new AliFMDBoolMap(0,0,0,0);
    dead->Reset(kFALSE);
  }
  MakeDead dm(dead, maxW);
  fWidth.ForEach(dm);
  AliFMDDebug(1, ("Found a total of %d dead channels", dm.fCount));
  return dead;
}

//____________________________________________________________________
Bool_t
AliFMDCalibPedestal::ReadFromFile(std::istream& in)
{
  //
  // Read information from file and set values
  // 
  // Parameters:
  //    inFile inputFile
  //
  TString header;
  header.ReadLine(in);
  header.ToLower();
  if(!header.Contains("pedestals")) {
    AliError("File does not contain pedestals!");
    return kFALSE;
  }
    
  // Read columns line
  int lineno = 2;
  header.ReadLine(in);
    
  // Loop until EOF
  while(in.peek()!=EOF) {
    if(in.bad()) { 
      AliError(Form("Bad read at line %d in input", lineno));
      break;
    }
    UShort_t det, sec, strip;
    Char_t ring;
    Float_t ped, noise, mu, sigma, chi2ndf;
    Char_t c[8];
	  
    in >> det      >> c[0] 
       >> ring     >> c[1]
       >> sec      >> c[2]
       >> strip    >> c[3]
       >> ped      >> c[4]
       >> noise    >> c[5]
       >> mu       >> c[6]
       >> sigma    >> c[7]
       >> chi2ndf;
    lineno++;
      
    Set(det,ring,sec,strip,ped,noise);
  }
  return kTRUE;
}
//____________________________________________________________________
//
// EOF
//
