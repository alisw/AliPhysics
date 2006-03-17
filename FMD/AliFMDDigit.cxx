/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////
//
//  Digits classes for the FMD                
//
//  Digits consists of
//   - Detector #
//   - Ring ID                                             
//   - Sector #     
//   - Strip #
//   - ADC count in this channel                                  
//
//  Digits consists of
//   - Detector #
//   - Ring ID                                             
//   - Sector #     
//   - Strip #
//   - Total energy deposited in the strip
//   - ADC count in this channel                                  
//
// As the Digits and SDigits have so much in common, the classes
// AliFMDDigit and AliFMDSDigit are implemented via a base
// class AliFMDBaseDigit.
///
//              +-----------------+
//              | AliFMDBaseDigit |
//              +-----------------+
//                      ^
//                      |
//                +------------+
//                |            |
//      +-------------+ +--------------+
//      | AliFMDDigit |	| AliFMDSDigit |
//      +-------------+	+--------------+
//
// (Note, that I'd really would have liked to implement AliFMDHit as a
// derived class from some base class - say AliFMDStrip, and the Digit
// classes would (eventually) have derived from that as well.
// However, ROOT doesn't do well with multiple inheritance, so I chose
// not to anyway).
//
// Latest changes by Christian Holm Christensen
//
//////////////////////////////////////////////////////////////////////

#include "AliFMDDigit.h"	// ALIFMDDIGIT_H
#include "Riostream.h"		// ROOT_Riostream
#include <TString.h>

//====================================================================
ClassImp(AliFMDBaseDigit)
#if 0
  ; // This is here to keep Emacs from indenting the next line
#endif

//____________________________________________________________________
AliFMDBaseDigit::AliFMDBaseDigit()
  : fDetector(0), 
    fRing('\0'), 
    fSector(0), 
    fStrip(0)
{}

//____________________________________________________________________
AliFMDBaseDigit::AliFMDBaseDigit(UShort_t detector, 
			 Char_t   ring, 
			 UShort_t sector, 
			 UShort_t strip)
  : fDetector(detector), 
    fRing(ring), 
    fSector(sector), 
    fStrip(strip)
{
  //
  // Creates a base data digit object
  //
  // Parameters 
  //
  //    detector  Detector # (1, 2, or 3)                      
  //    ring	  Ring ID ('I' or 'O')
  //    sector	  Sector # (For inner/outer rings: 0-19/0-39)
  //    strip	  Strip # (For inner/outer rings: 0-511/0-255)
}

//____________________________________________________________________
void
AliFMDBaseDigit::Print(Option_t* /* option*/) const 
{
  // Print digit to standard out 
  cout << ClassName() << ": FMD" << fDetector << fRing << "[" 
       << setw(3) << fSector << ","
       << setw(3) << fStrip << "]" 
       << flush;
}

//____________________________________________________________________
const char*
AliFMDBaseDigit::GetName() const 
{ 
  static TString n;
  n = Form("FMD%d%c[%2d,%3d]", fDetector, fRing, fSector, fStrip);
  return n.Data();
}

//====================================================================
ClassImp(AliFMDDigit)

//____________________________________________________________________
AliFMDDigit::AliFMDDigit()
  : fCount1(0),
    fCount2(-1),
    fCount3(-1)
{}

//____________________________________________________________________
AliFMDDigit::AliFMDDigit(UShort_t detector, 
			 Char_t   ring, 
			 UShort_t sector, 
			 UShort_t strip, 
			 UShort_t count1,
			 Short_t  count2, 
			 Short_t  count3)
  : AliFMDBaseDigit(detector, ring, sector, strip), 
    fCount1(count1),
    fCount2(count2),
    fCount3(count3)
{
  //
  // Creates a real data digit object
  //
  // Parameters 
  //
  //    detector  Detector # (1, 2, or 3)                      
  //    ring	  Ring ID ('I' or 'O')
  //    sector	  Sector # (For inner/outer rings: 0-19/0-39)
  //    strip	  Strip # (For inner/outer rings: 0-511/0-255)
  //    count1    ADC count (a 10-bit word)
  //    count2    ADC count (a 10-bit word) -1 if not used
  //    count3    ADC count (a 10-bit word) -1 if not used
}

//____________________________________________________________________
const char*
AliFMDDigit::GetTitle() const 
{ 
  static TString t;
  t = Form("ADC: %d", Counts());
  return t.Data();
}

//____________________________________________________________________
void
AliFMDDigit::Print(Option_t* /* option*/) const 
{
  // Print digit to standard out 
  AliFMDBaseDigit::Print();
  cout << "\t" 
       << fCount1 << " (+ " << fCount2 << " + " << fCount2 << ") = " 
       << Counts() << endl;
}

//====================================================================
ClassImp(AliFMDSDigit)

//____________________________________________________________________
AliFMDSDigit::AliFMDSDigit()
  : fEdep(0), 
    fCount1(0),
    fCount2(-1),
    fCount3(-1)
{}

//____________________________________________________________________
AliFMDSDigit::AliFMDSDigit(UShort_t detector, 
			   Char_t   ring, 
			   UShort_t sector, 
			   UShort_t strip, 
			   Float_t  edep,
			   UShort_t count1,
			   Short_t  count2, 
			   Short_t  count3)
  : AliFMDBaseDigit(detector, ring, sector, strip), 
    fEdep(edep),
    fCount1(count1),
    fCount2(count2),
    fCount3(count3)
{
  //
  // Creates a real data digit object
  //
  // Parameters 
  //
  //    detector  Detector # (1, 2, or 3)                      
  //    ring	  Ring ID ('I' or 'O')
  //    sector	  Sector # (For inner/outer rings: 0-19/0-39)
  //    strip	  Strip # (For inner/outer rings: 0-511/0-255)
  //    edep      Total energy deposited 
  //    count1    ADC count (a 10-bit word)
  //    count2    ADC count (a 10-bit word) -1 if not used
  //    count3    ADC count (a 10-bit word) -1 if not used
}

//____________________________________________________________________
void
AliFMDSDigit::Print(Option_t* /* option*/) const 
{
  // Print digit to standard out 
  AliFMDBaseDigit::Print();
  cout << "\t" << fEdep << " -> "
       << fCount1 << " (+ " << fCount2 << " + " << fCount2 << ") = " 
       << Counts() << endl;
}

//____________________________________________________________________
//
// EOF
//
