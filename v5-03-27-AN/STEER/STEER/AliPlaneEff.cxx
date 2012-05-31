/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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
///////////////////////////////////////////////////////////////////////////
//  Virtual base Plane Efficiency class           
//  Specific detector implementation is done in  
//  AliITSPlaneEff                           
//  AliTRDPlaneEff                               
//  AliTOFPlaneEff                               
//
//  Author: G.E. Bruno 
//          giuseppe.bruno@ba.infn.it
//
///////////////////////////////////////////////////////////////////////////


#include <TMath.h>
#include "AliPlaneEff.h"
#include "AliLog.h"
//#include "AliCDBManager.h"
//#include "AliCDBStorage.h"

ClassImp(AliPlaneEff)
//______________________________________________________________________
AliPlaneEff::AliPlaneEff(): TObject()/*,
fRunNumber(0), 
fCDBUri(""),
fInitCDBCalled(kFALSE)*/
{
    // Default constructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    a default constructed AliPlaneEff class
 //InitCDB();
}
//______________________________________________________________________
AliPlaneEff::AliPlaneEff(const AliPlaneEff &s) : TObject(s)/*,
fRunNumber(s.fRunNumber),
fCDBUri(s.fCDBUri),
fInitCDBCalled(s.fInitCDBCalled)*/
{
    //     Copy Constructor
    // Inputs:
    //    const AliPlaneEff &s  simulation class to copy from
    // Outputs:
    //    none.
    // Return:
    //    a standard constructed AliPlaneEff class with values the same
    //    as that of s.

}
//_________________________________________________________________________
AliPlaneEff&  AliPlaneEff::operator=(const AliPlaneEff &source){
    //    Assignment operator
    // Inputs:
    //    const AliPlaneEff &source  simulation class to copy from
    // Outputs:
    //    none.
    // Return:
    //    a standard constructed AliPlaneEff class with values the same
    //    as that of s.
    if(this != &source){
       source.Copy(*this);
    }
    return *this;
}
//_____________________________________________________________
void AliPlaneEff::Copy(TObject &/*obj*/) const {
  // copy this to obj
/*  ((AliPlaneEff& ) obj).fRunNumber		= fRunNumber;
  ((AliPlaneEff& ) obj).fCDBUri		= fCDBUri;
  ((AliPlaneEff& ) obj).fInitCDBCalled	= fInitCDBCalled;
*/
}
//_________________________________________________________________________
