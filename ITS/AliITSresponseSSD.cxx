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


#include "AliITSresponseSSD.h"
//////////////////////////////////////////////////
//  Response class for set:ITS                      //
//  Specific subdetector implementation             //
//  for silicon strips detectors                    //
//                                                  //
//                                                  //
//////////////////////////////////////////////////////

const Double_t AliITSresponseSSD::fgkDiffCoeffDefault = 0.;
const TString AliITSresponseSSD::fgkOption1Default = "";
const TString AliITSresponseSSD::fgkOption2Default = "";
const Double_t AliITSresponseSSD::fgkNoiseNDefault = 625.;
const Double_t AliITSresponseSSD::fgkNoisePDefault = 420.;
const Int_t AliITSresponseSSD::fgkNParDefault = 6;
const Double_t AliITSresponseSSD::fgkSigmaPDefault = 3.;
const Double_t AliITSresponseSSD::fgkSigmaNDefault = 2.;

ClassImp(AliITSresponseSSD)

//______________________________________________________________________
AliITSresponseSSD::AliITSresponseSSD(){
    // Default Constructor

    fDetPar = 0;
    fNPar   = 0;
    fNoiseP = 0;
    fNoiseN = 0;
    fSigmaP = 0;
    fSigmaN = 0;
    fDiffCoeff = 0;
    fADCpereV  = 0;
    SetParamOptions(fgkOption1Default.Data(),fgkOption2Default.Data());
    SetNoiseParam(fgkNoisePDefault,fgkNoiseNDefault);
}
//______________________________________________________________________
AliITSresponseSSD::AliITSresponseSSD(const char *dataType){
    // constructor

    SetDiffCoeff(fgkDiffCoeffDefault,0.);
    SetNoiseParam(fgkNoisePDefault,fgkNoiseNDefault);
    SetDataType(dataType);
    SetSigmaSpread(fgkSigmaPDefault,fgkSigmaNDefault);
    SetParamOptions(fgkOption1Default.Data(),fgkOption2Default.Data());
    SetNDetParam(fgkNParDefault);   // Sets fNPar=6 by default.
    SetADCpereV();
    fDetPar = new Double_t[fNPar];
    if (fNPar==6) {
	fDetPar[0]=10.;
	fDetPar[1]=5.;
	fDetPar[2]=0.02;
	fDetPar[3]=0.02;
	fDetPar[4]=0.02;
	fDetPar[5]=0.03;
    } // end if
}
//______________________________________________________________________
AliITSresponseSSD::~AliITSresponseSSD(){
    // destructor

    delete [] fDetPar;
}
//______________________________________________________________________
AliITSresponseSSD& AliITSresponseSSD::operator=(const AliITSresponseSSD &src) {
    // = operator.

    if(&src == this) return *this;

    this->fNPar      = src.fNPar;
    for(Int_t i=0;i<this->fNPar;i++) this->fDetPar[i] = src.fDetPar[i];
    this->fNoiseP    = src.fNoiseP;
    this->fNoiseN    = src.fNoiseN;
    this->fSigmaP    = src.fSigmaP;
    this->fSigmaN    = src.fSigmaN;
    this->fDiffCoeff = src.fDiffCoeff;
    this->fADCpereV  = src.fADCpereV;
    this->fOption1   = src.fOption1;
    this->fOption2   = src.fOption2;
    this->fDataType  = src.fDataType;

    return *this;
}
//_________________________________________________________________________
AliITSresponseSSD::AliITSresponseSSD(const AliITSresponseSSD &src) :
    AliITSresponse(src) {
    // copy constructor

    *this = src;
}
//______________________________________________________________________
void AliITSresponseSSD::SetDetParam(Double_t  *par){
    // set det param
    Int_t i;

    for (i=0; i<fNPar; i++) {
	fDetPar[i]=par[i];
	//printf("\n CompressPar %d %d \n",i,fCPar[i]);    
    } // end for i
}
//______________________________________________________________________
void AliITSresponseSSD::GetDetParam(Double_t  *par) const {
    // get det param
    Int_t i;

    for (i=0; i<fNPar; i++) {
	par[i]=fDetPar[i];
    } // end for i
}
