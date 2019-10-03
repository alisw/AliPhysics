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

//--------------------------------------------------
// Filling of CalTrkEvent objects in the reader task
//
// Author: magali.estienne@subatech.in2p3.fr
//         alexandre.shabetai@cern.ch
//-------------------------------------------------

#include "AliJetFillCalTrkEvent.h"

// --- ROOT system ---
class TSystem;
class TLorentzVector;
class TVector3;
class TGeoManager;

// --- AliRoot header files ---
class AliJetFinder;
class AliJetReader;
class AliJetReader;
class AliJetReaderHeader;
class AliJetCalTrkEvent;

ClassImp(AliJetFillCalTrkEvent)

/////////////////////////////////////////////////////////////////////

AliJetFillCalTrkEvent::AliJetFillCalTrkEvent() :
    fOpt(0),
    fDebug(0),
    fReaderHeader(0x0),
    fCalTrkEvent(0x0),
    fGeom(0x0)
{
  // constructor
}

//-----------------------------------------------------------------------
AliJetFillCalTrkEvent::AliJetFillCalTrkEvent(const AliJetFillCalTrkEvent& cpfrom) :
    fOpt(0),
    fDebug(0),
    fReaderHeader(0x0),
    fCalTrkEvent(0x0),
    fGeom(0x0)
{
    // Copy constructor
    fOpt              = cpfrom.fOpt;
    fDebug            = cpfrom.fDebug;
    fReaderHeader     = cpfrom.fReaderHeader;
    fCalTrkEvent     = cpfrom.fCalTrkEvent;
    fGeom             = cpfrom.fGeom;
}

//-----------------------------------------------------------------------
AliJetFillCalTrkEvent& AliJetFillCalTrkEvent::operator=(const AliJetFillCalTrkEvent& rhs)
{
   // Assignment operator
  if (this != &rhs) {
    fOpt              = rhs.fOpt;
    fDebug            = rhs.fDebug;
    fReaderHeader     = rhs.fReaderHeader;
    fCalTrkEvent     = rhs.fCalTrkEvent;
    fGeom             = rhs.fGeom;
  }
  
  return *this;
    
}

//-----------------------------------------------------------------------
AliJetFillCalTrkEvent::~AliJetFillCalTrkEvent()
{
  // destructor
}

//-----------------------------------------------------------------------
Float_t  AliJetFillCalTrkEvent::EtaToTheta(Float_t arg)
{
  // Transform eta to theta
  return 2.*atan(exp(-arg));
}







