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


// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2F.h> 
#include <TLorentzVector.h>

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAChecker.h"

#include "AliMUONQADataMakerSim.h"

//-----------------------------------------------------------------------------
/// \class AliMUONQADataMakerSim
///
/// MUON base class for quality assurance data (histo) maker
///
/// \author C. Finck

/// \cond CLASSIMP
ClassImp(AliMUONQADataMakerSim)
/// \endcond
           
//____________________________________________________________________________ 
AliMUONQADataMakerSim::AliMUONQADataMakerSim() : 
    AliQADataMakerSim(AliQA::GetDetName(AliQA::kMUON), "MUON Quality Assurance Data Maker")
{
    /// ctor
}

//____________________________________________________________________________ 
AliMUONQADataMakerSim::AliMUONQADataMakerSim(const AliMUONQADataMakerSim& qadm) :
    AliQADataMakerSim()
{
    ///copy ctor 
    SetName((const char*)qadm.GetName()) ; 
    SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliMUONQADataMakerSim& AliMUONQADataMakerSim::operator = (const AliMUONQADataMakerSim& qadm )
{
    /// Equal operator.
    this->~AliMUONQADataMakerSim();
    new(this) AliMUONQADataMakerSim(qadm);
    return *this;
}

//__________________________________________________________________
AliMUONQADataMakerSim::~AliMUONQADataMakerSim()
{
    /// dtor
}

//____________________________________________________________________________ 
void AliMUONQADataMakerSim::EndOfDetectorCycle(AliQA::TASKINDEX task, TObjArray* list)
{
    ///Detector specific actions at end of cycle
    // do the QA checking
    AliQAChecker::Instance()->Run(AliQA::kMUON, task, list) ;  
}


//____________________________________________________________________________ 
void AliMUONQADataMakerSim::StartOfDetectorCycle()
{
    /// Detector specific actions at start of cycle
  
}
