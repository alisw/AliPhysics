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
/** @file    AliFMDSDigit.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:37:41 2006
    @brief   Digits for the FMD 
    @ingroup FMD_base
*/
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

#include "AliFMDSDigit.h"	// ALIFMDDIGIT_H
#include "Riostream.h"		// ROOT_Riostream
#include <TString.h>

//====================================================================
ClassImp(AliFMDSDigit)
#if 0
; // Here to make Emacs happy
#endif
//____________________________________________________________________
AliFMDSDigit::AliFMDSDigit()
  : fEdep(0), 
    fCount1(0),
    fCount2(-1),
    fCount3(-1), 
    fCount4(-1), 
    fNParticles(0),
    fNPrimaries(0)
    // ,     fLabels(0)
{
  // cTOR 
}

//____________________________________________________________________
AliFMDSDigit::AliFMDSDigit(UShort_t       detector, 
			   Char_t         ring, 
			   UShort_t       sector, 
			   UShort_t       strip, 
			   Float_t        edep,
			   UShort_t       count1,
			   Short_t        count2, 
			   Short_t        count3, 
			   Short_t        count4,
			   UShort_t       npart,
			   UShort_t       nprim,
			   const Int_t*   refs)
  : AliFMDBaseDigit(detector, ring, sector, strip), 
    fEdep(edep),
    fCount1(count1),
    fCount2(count2),
    fCount3(count3),
    fCount4(count4),
    fNParticles(npart), 
    fNPrimaries(nprim)
    // , fLabels(refs)
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
  if (!refs) return;
  for (Int_t i = 0; i < npart; i++) AddTrack(refs[i]);
}

//____________________________________________________________________
void
AliFMDSDigit::Print(Option_t* option) const 
{
  // Print digit to standard out 
  AliFMDBaseDigit::Print();
  std::cout << "\t" 
	    << std::setw(10) << fEdep << " -> "
	    << std::setw(4) << fCount1 << " (" 
	    << std::setw(4) << fCount2 << "," 
	    << std::setw(4) << fCount3 << "," 
	    << std::setw(4) << fCount4 << ") = " 
	    << std::setw(4) << Counts() << std::flush;

  TString opt(option);
  if (opt.Contains("p", TString::kIgnoreCase)) 
    std::cout << " [" 
	      << std::setw(2) << fNPrimaries << "/" 
	      << std::setw(2) << fNParticles << "]"
	      << std::flush;
  if (opt.Contains("l", TString::kIgnoreCase)) {
    std::cout << " ";
    for (Int_t i = 0; i < GetNTrack(); i++) 
      std::cout << (i == 0 ? "" : ",") << std::setw(5) << fTracks[i];
#if 0
    for (Int_t i = 0; i < fLabels.fN; i++) 
      std::cout << (i == 0 ? "" : ",") << fLabels.fArray[i];
#endif
  }
  std::cout << std::endl;
}

//____________________________________________________________________
//
// EOF
//
