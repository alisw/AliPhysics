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

////////////////////////////////////////////////////////////////////////
//
// AliTOFConstants class
//
// This class serves to group constants needed by TOF detector in 1
// easily accessible place. All constants are public const static data 
// members. The class is never instatiated.
//
// Note: only a few constants are in the first version of this class,
//       more should be added by TOF developpers
//
// Author: Jiri Chudoba (CERN), F. Pierella
//
////////////////////////////////////////////////////////////////////////

#include "AliTOFConstants.h"

const Int_t AliTOFConstants::fgkNStripA = 15;
const Int_t AliTOFConstants::fgkNStripB = 19;
const Int_t AliTOFConstants::fgkNStripC = 20;
const Int_t AliTOFConstants::fgkNpadX   = 48;
const Int_t AliTOFConstants::fgkNpadZ   =  2;
const Int_t AliTOFConstants::fgkPadXSector =
      (fgkNStripA + 2*fgkNStripB + 2*fgkNStripC)*fgkNpadX*fgkNpadZ;
const Int_t AliTOFConstants::fgkNSectors   =  18;
const Int_t AliTOFConstants::fgkNPlates    =  5;

const Float_t AliTOFConstants::fgkrmin     = 370.;
const Float_t AliTOFConstants::fgkrmax     = 399.;
const Int_t AliTOFConstants::fgkmaxtoftree = 5;
const Int_t AliTOFConstants::fgkmaxNstrip  = 20;
const Int_t AliTOFConstants::fgkPadXStrip  = fgkNpadX*fgkNpadZ;
const Float_t AliTOFConstants::fgkzlenA    = 106.0;
const Float_t AliTOFConstants::fgkzlenB    = 141.0;
const Float_t AliTOFConstants::fgkzlenC    = 177.5;
const Float_t AliTOFConstants::fgkXPad     = 2.5;
const Float_t AliTOFConstants::fgkZPad     = 3.5;
const Float_t AliTOFConstants::fgkMaxhZtof = 371.5;
const Float_t AliTOFConstants::fgkSigmaForTail1= 2.;
const Float_t AliTOFConstants::fgkSigmaForTail2=0.5;
const Int_t AliTOFConstants::fgkTimeDiff   =  25000;
const Float_t AliTOFConstants::fgkSpeedOfLight = 0.299792458;
const Float_t AliTOFConstants::fgkPionMass     = 0.13957;
const Float_t AliTOFConstants::fgkKaonMass     = 0.49368;
const Float_t AliTOFConstants::fgkProtonMass   = 0.93827;
const Float_t AliTOFConstants::fgkElectronMass = 0.00051;
const Float_t AliTOFConstants::fgkMuonMass     = 0.10566;
ClassImp(AliTOFConstants)

