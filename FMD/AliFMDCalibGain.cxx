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
/** @file    AliFMDCalibGain.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:30:02 2006
    @brief   Per strip gain calibration 
*/
//____________________________________________________________________
//                                                                          
// Gain value and width for each strip in the FMD. 
// Foo 
// Bar 
// Baz
// Gnus
//
#include "AliFMDCalibGain.h"	// ALIFMDCALIBGAIN_H
#include <iostream>
#include <TString.h>
#include <AliLog.h>

//____________________________________________________________________
ClassImp(AliFMDCalibGain)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDCalibGain::AliFMDCalibGain()
  : fValue(0), // nDet == 0 mean 51200 slots
    fThreshold(-1.)
{
  // CTOR
  fValue.Reset(-1.);
  fThreshold = -1.;
}

//____________________________________________________________________
AliFMDCalibGain::AliFMDCalibGain(const AliFMDCalibGain& o)
  : TObject(o), 
    fValue(o.fValue), 
    fThreshold(o.fThreshold)
{
  // Copy CTOR 
}

//____________________________________________________________________
AliFMDCalibGain&
AliFMDCalibGain::operator=(const AliFMDCalibGain& o)
{
  // Assignment operator 
  fValue     = o.fValue;
  fThreshold = o.fThreshold;
  return (*this);
}

//____________________________________________________________________
void
AliFMDCalibGain::Set(UShort_t det, Char_t ring, UShort_t sec, 
		     UShort_t str, Float_t val)
{
  // Set the value for a strip 
  if (fValue.CheckIndex(det, ring, sec, str) < 0) return;
  fValue(det, ring, sec, str) = val;
}

//____________________________________________________________________
Float_t
AliFMDCalibGain::Value(UShort_t det, Char_t ring, UShort_t sec, 
		       UShort_t str)
{
  // Get the value for a strip 
  return fValue(det, ring, sec, str);
}

//____________________________________________________________________
Bool_t
AliFMDCalibGain::ReadFromFile(std::istream& in)
{
  //Get header (how long is it ?)
  TString header;
  header.ReadLine(in);
  header.ToLower();
  if(!header.Contains("gains")) {
    AliError("File does not contain gains!");
    return kFALSE;;
  }

  // Read column headers
  header.ReadLine(in);
  
  int lineno  = 2;
  // Read until EOF 
  while(in.peek()!=EOF) {
    if(in.bad()) { 
      AliError(Form("Bad read at line %d of input", lineno));
      break;
    }
    UShort_t det, sec, strip;
    Char_t ring;
    
    Float_t gain,error,  chi2ndf;
    Char_t c[6];
    
    in >> det      >> c[0] 
       >> ring     >> c[1]
       >> sec      >> c[2]
       >> strip    >> c[3]
       >> gain     >> c[4]
       >> error    >> c[5]
       >> chi2ndf;
    lineno++;
    Set(det,ring,sec,strip,gain);
  }
  return kTRUE;
}
//____________________________________________________________________
//
// EOF
//
