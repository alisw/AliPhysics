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

//////////////////////////////////////////////////////
//  Calibration class for set:ITS                   //
//  Specific subdetector implementation is done in  //
//  AliITSCalibrationSPD                            //
//  AliITSCalibrationSDD                            //
//  AliITSCalibrationSSD                            //
//////////////////////////////////////////////////////

#include <Riostream.h>
#include <TMath.h>

#include "AliITSCalibration.h"
#include "AliLog.h"

ClassImp(AliITSCalibration)

//______________________________________________________________________
AliITSCalibration::AliITSCalibration():
TObject(),
fDataType(),
fT(300.)
{
    // Default Constructor (300 microns and 80 volts)

}

//______________________________________________________________________
AliITSCalibration::AliITSCalibration(const AliITSCalibration &ob):
TObject(ob),
fDataType(ob.fDataType),
fT(ob.fT)

{
  // Copy constructor

}
//----------------------------------------------------------------------
void AliITSCalibration::Print(ostream *os) const {
  // Standard output format for this class.
  // Inputs:
    *os << fT << " ";
  //    printf("%-10.6e  %-10.6e %-10.6e %-10.6e \n",fdv,fN,fT,fGeVcharge);
    return;
}
//----------------------------------------------------------------------
void AliITSCalibration::Read(istream *is) {
  // Standard input format for this class.
  // Inputs:
  //    ostream *is  Pointer to the output stream
  // Outputs:
  //    none:
  // Return:
  //    none.

    *is >> fT;
    return;
}
//----------------------------------------------------------------------

ostream &operator<<(ostream &os,AliITSCalibration &p){
  // Standard output streaming function.
  // Inputs:
  //    ostream *os  Pointer to the output stream
  // Outputs:
  //    none:
  // Return:
  //    none.

    p.Print(&os);
    return os;
}

//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSCalibration &r){
  // Standard input streaming function.
  // Inputs:
  //    ostream *os  Pointer to the output stream
  // Outputs:
  //    none:
  // Return:
  //    none.

    r.Read(&is);
    return is;
}
//----------------------------------------------------------------------
