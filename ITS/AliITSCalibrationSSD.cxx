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


#include "AliITSCalibrationSSD.h"
//////////////////////////////////////////////////////
//  Calibration class for set:ITS                   //
//  Specific subdetector implementation             //
//  for silicon strips detectors                    //
//                                                  //
//                                                  //
//////////////////////////////////////////////////////

const Double_t AliITSCalibrationSSD::fgkNoiseNDefault = 625.;
const Double_t AliITSCalibrationSSD::fgkNoisePDefault = 420.;
const Int_t AliITSCalibrationSSD::fgkNParDefault = 6;
const Double_t AliITSCalibrationSSD::fgkSigmaPDefault = 3.;
const Double_t AliITSCalibrationSSD::fgkSigmaNDefault = 2.;

ClassImp(AliITSCalibrationSSD)

//______________________________________________________________________
AliITSCalibrationSSD::AliITSCalibrationSSD(){
    // Default Constructor

    fDetPar = 0;
    fNPar   = 0;
    fNoiseP = 0;
    fNoiseN = 0;
    fSigmaP = 0;
    fSigmaN = 0;
    SetNoiseParam(fgkNoisePDefault,fgkNoiseNDefault);
}
//______________________________________________________________________
AliITSCalibrationSSD::AliITSCalibrationSSD(const char *dataType){
    // constructor

    SetNoiseParam(fgkNoisePDefault,fgkNoiseNDefault);
    SetDataType(dataType);
    SetSigmaSpread(fgkSigmaPDefault,fgkSigmaNDefault);
    SetNDetParam(fgkNParDefault);   // Sets fNPar=6 by default.
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
AliITSCalibrationSSD::~AliITSCalibrationSSD(){
    // destructor

    delete [] fDetPar;
}
//______________________________________________________________________
AliITSCalibrationSSD& AliITSCalibrationSSD::operator=(const AliITSCalibrationSSD &src) {
    // = operator.

    if(&src == this) return *this;

    this->fNPar      = src.fNPar;
    for(Int_t i=0;i<this->fNPar;i++) this->fDetPar[i] = src.fDetPar[i];
    this->fNoiseP    = src.fNoiseP;
    this->fNoiseN    = src.fNoiseN;
    this->fSigmaP    = src.fSigmaP;
    this->fSigmaN    = src.fSigmaN;
    this->fDataType  = src.fDataType;

    return *this;
}
//_________________________________________________________________________
AliITSCalibrationSSD::AliITSCalibrationSSD(const AliITSCalibrationSSD &src) :
    AliITSCalibration(src) {
    // copy constructor

    *this = src;
}
//______________________________________________________________________
void AliITSCalibrationSSD::SetDetParam(Double_t  *par){
    // set det param
    Int_t i;

    for (i=0; i<fNPar; i++) {
	fDetPar[i]=par[i];
	//printf("\n CompressPar %d %d \n",i,fCPar[i]);    
    } // end for i
}
//______________________________________________________________________
void AliITSCalibrationSSD::GetDetParam(Double_t  *par) const {
    // get det param
    Int_t i;

    for (i=0; i<fNPar; i++) {
	par[i]=fDetPar[i];
    } // end for i
}
