/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                         *
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

#include <TMath.h>

#include "AliITSresponseSSD.h"
#include "AliITSgeom.h"

ClassImp(AliITSresponseSSD)	
//----------------------------------------------------------
AliITSresponseSSD::AliITSresponseSSD()
{
  // constructor
  SetDiffCoeff();
  SetNoiseParam();
  SetDataType();
  SetSigmaSpread();
  SetParamOptions();
  SetNDetParam();
  fDetPar = new Float_t[fNPar];
  if (fNPar==6) {
    fDetPar[0]=10.;
    fDetPar[1]=5.;
    fDetPar[2]=0.02;
    fDetPar[3]=0.02;
    fDetPar[4]=0.02;
    fDetPar[5]=0.03;
  }
  
  
}

//----------------------------------------------------------
AliITSresponseSSD::~AliITSresponseSSD()
{
  // destructor
  delete [] fDetPar;
  delete fDetPar;
  
}

//__________________________________________________________________________
AliITSresponseSSD::AliITSresponseSSD(const AliITSresponseSSD &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fNPar = source.fNPar;
  this->fDetPar = source.fDetPar;
  this->fNoiseP = source.fNoiseP;
  this->fNoiseN = source.fNoiseN;
  this->fSigmaP = source.fSigmaP;
  this->fSigmaN = source.fSigmaN;
  this->fDiffCoeff = source.fDiffCoeff;
  this->fOption1 = source.fOption1;
  this->fOption2 = source.fOption2;
  this->fDataType = source.fDataType;
  return;
}

//_________________________________________________________________________
AliITSresponseSSD& 
  AliITSresponseSSD::operator=(const AliITSresponseSSD &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fNPar = source.fNPar;
  this->fDetPar = source.fDetPar;
  this->fNoiseP = source.fNoiseP;
  this->fNoiseN = source.fNoiseN;
  this->fSigmaP = source.fSigmaP;
  this->fSigmaN = source.fSigmaN;
  this->fDiffCoeff = source.fDiffCoeff;
  this->fOption1 = source.fOption1;
  this->fOption2 = source.fOption2;
  this->fDataType = source.fDataType;
  return *this;
}

//----------------------------------------------------------
void AliITSresponseSSD::SetDetParam(Float_t  *par)
{
  // set det param
  Int_t i;
  for(i=0; i<fNPar; i++) {
    fDetPar[i]=par[i];
    //printf("\n CompressPar %d %d \n",i,fCPar[i]);
    
  } 
}
void AliITSresponseSSD::GetDetParam(Float_t  *par)
{
  // get det param
  Int_t i;  
  for(i=0; i<fNPar; i++) {
    par[i]=fDetPar[i];
  }
}
