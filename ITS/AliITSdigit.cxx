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

/*
$Log$
*/

//
#include "AliITS.h"
#include "AliITSdigit.h"
#include "AliRun.h"
////////////////////////////////////////////////////////////////////////
// Version: 0
// Written by Rene Brun, Federico Carminati, and Roberto Barbera
// Minor modifications made and documented by Bjorn S. Nilsen
// July 11 1999
//
// The default ITS digit structure. This should either be replaced
// or added on to later with the proper digit structure defined for
// each detector type. See the proposed Digit structure defined by
// Bjorn S. Nilsen for an example.
//Begin_Html
/*
<img src="picts/ITS/AliITShit_Class_Diagram.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This show the relasionships between the ITS digit class' and the rest of Aliroot.
</font>
<pre>
*/
//End_Html
//_____________________________________________________________________________
ClassImp(AliITSdigit)

AliITSdigit::AliITSdigit(Int_t *tracks, Int_t *digits):
  AliDigit(tracks){
  //
  // Create ITS digit
  //     The creator for the AliITSdigit class. This routine fills the
  // AliITSdigit data members from the array digits. The array of track
  // numbers are passed to the AliDigit creator. The order of the elements
  // in the digits array are fEvent = digits[0], fLayer = digits[1],
  // fLadder = digits[2], fDet = digits[3], and fNoverl = digits[4].
  // Therefore the array digits is expected to be at least 5 elements long.
  //
  fEvent      = digits[0];
  fLayer      = digits[1];
  fLadder     = digits[2];
  fDet        = digits[3];
  fNoverl     = digits[4];
}
