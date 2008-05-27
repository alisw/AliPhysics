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

// $Id$

#include "AliQAChecker.h"
#include "AliMUONQAChecker.h"

//-----------------------------------------------------------------------------
/// \class AliMUONQAChecker
///
/// MUON base class for quality assurance checker
///
/// \author C. Finck
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONQAChecker)
/// \endcond

//__________________________________________________________________
AliMUONQAChecker::AliMUONQAChecker() : 
    AliQACheckerBase("MUON","MUON Quality Assurance Data Maker") 
{
/// ctor

}          

//__________________________________________________________________
AliMUONQAChecker::~AliMUONQAChecker() 
{
/// dtor
}

//__________________________________________________________________
AliMUONQAChecker::AliMUONQAChecker(const AliMUONQAChecker& qac) : 
    AliQACheckerBase(qac.GetName(), qac.GetTitle()) 
{
/// copy ctor 

}   

//__________________________________________________________________
AliMUONQAChecker& AliMUONQAChecker::operator = (const AliMUONQAChecker& /*qac*/ )
{
    /// Equal operator.
    return *this;
}

