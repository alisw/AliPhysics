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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Container for PID information                                        //
//                                                                      //
// Authors:                                                             //
//   Prashant Shukla <shukla@pi0.physi.uni-heidelberg.de>               //
//   Alex Bercuci <a.bercuci@gsi.de>                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TROOT.h>

#include "AliLog.h"
#include "AliPID.h"
#include "AliESD.h"
#include "AliESDtrack.h"

#include "AliTRDCalPID.h"
#include "AliTRDcalibDB.h"

ClassImp(AliTRDCalPID)

Char_t* AliTRDCalPID::fPartName[AliPID::kSPECIES] = { "electron", "muon", "pion", "kaon", "proton"};
Char_t* AliTRDCalPID::fPartSymb[AliPID::kSPECIES] = { "EL", "MU", "PI", "KA", "PR"};
Float_t AliTRDCalPID::fTrackMomentum[kNMom]       = {  0.6,  0.8,  1.0,  1.5,  2.0
                                                   ,   3.0,  4.0,  5.0,  6.0,  8.0
                                                   ,  10.0};

//_________________________________________________________________________
AliTRDCalPID::AliTRDCalPID()
  :TNamed("pid", "PID for TRD")
  ,fModel(0x0)
{
  //
  //  The Default constructor
  //

}

//_____________________________________________________________________________
AliTRDCalPID::AliTRDCalPID(const Text_t *name, const Text_t *title) 
  :TNamed(name,title)
  ,fModel(0x0)
{
  //
  //  The main constructor
  //  

}

//_________________________________________________________________________
AliTRDCalPID::~AliTRDCalPID()
{
  //
  // Destructor
  //
  
  if (fModel) {
    delete fModel;
  }

}
