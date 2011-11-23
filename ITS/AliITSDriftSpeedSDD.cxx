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

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class for SDD drift speed               //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSDriftSpeedSDD.h"
#include "AliLog.h"

const Float_t AliITSDriftSpeedSDD::fgkDriftSpeedDefault = 7.3;

ClassImp(AliITSDriftSpeedSDD)
//______________________________________________________________________
AliITSDriftSpeedSDD::AliITSDriftSpeedSDD():
TObject(),
fEvNum(0),
fTimestamp(0),
fPolDeg(0){
  // default constructor
  fDriftSpeedParam[0]=fgkDriftSpeedDefault;
  for(Int_t i=1; i<fgkMaxPolDeg+1; i++) fDriftSpeedParam[i]=0;
}
//______________________________________________________________________
AliITSDriftSpeedSDD::AliITSDriftSpeedSDD(Int_t ev, UInt_t timest, Int_t deg, Double_t *coeff):
TObject(),
fEvNum(ev),
fTimestamp(timest),
fPolDeg(deg){
  // standard constructor
  if(deg>fgkMaxPolDeg){
    fPolDeg=fgkMaxPolDeg;
    AliWarning(Form("Polynomial degree must be <%d. Drift speed parametrization limited to %dth degree poly.\n",fgkMaxPolDeg,fgkMaxPolDeg));
  }
  for(Int_t i=0; i<fPolDeg+1; i++) fDriftSpeedParam[i]=coeff[i];
  for(Int_t i=fPolDeg+1; i<fgkMaxPolDeg+1; i++) fDriftSpeedParam[i]=0;
}
//______________________________________________________________________
AliITSDriftSpeedSDD::AliITSDriftSpeedSDD(const AliITSDriftSpeedSDD& drSpeed):
TObject(),
fEvNum(drSpeed.fEvNum),
fTimestamp(drSpeed.fTimestamp),
fPolDeg(drSpeed.fPolDeg)
{
  // copy constructor
  for(Int_t i=0; i<fgkMaxPolDeg+1; i++) fDriftSpeedParam[i]=drSpeed.GetDriftSpeedParameter(i);
  
}
//_____________________________________________________________________________
AliITSDriftSpeedSDD& AliITSDriftSpeedSDD::operator=(const AliITSDriftSpeedSDD &drSpeed){
  // Assignment operator
 if(this==&drSpeed) return *this;
  ((TObject *)this)->operator=(drSpeed);
  fEvNum = drSpeed.fEvNum;
  fTimestamp = drSpeed.fTimestamp;
  fPolDeg = drSpeed.fPolDeg;
  for(Int_t i=0; i<fgkMaxPolDeg+1; i++) fDriftSpeedParam[i]=drSpeed.GetDriftSpeedParameter(i);
}

//______________________________________________________________________
void AliITSDriftSpeedSDD::PrintDriftSpeedParameters() const {
  // printout drift speed parametrization
  printf("Injector event #%d at time %d\n",fEvNum,fTimestamp);
  printf("Coefficients of %d degree poly fit:\n",fPolDeg);
  for(Int_t i=0; i<fgkMaxPolDeg+1; i++) printf("par[%d]=%G\n",i,fDriftSpeedParam[i]);
}
