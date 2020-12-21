/**************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
#include <iostream>
#include <TString.h>
#include <AliConversionCutHandler.h>

AliConversionCutHandler::AliConversionCutHandler(Int_t nMax) :
fValidCuts(true),
fNCuts(0),
fNMaxCuts(nMax),
fEventCutArray(new TString[fNMaxCuts]),
fPhotonCutArray(new TString[fNMaxCuts]),
fMesonCutArray(new TString[fNMaxCuts]),
fClusterCutArray(new TString[fNMaxCuts])
{
  for(Int_t i=0; i<fNMaxCuts; i++) {
    fEventCutArray[i] = "";
    fPhotonCutArray[i] = "";
    fMesonCutArray[i] = "";
    fClusterCutArray[i] = "";
  }
}

AliConversionCutHandler::AliConversionCutHandler(const AliConversionCutHandler &ref):
fValidCuts(ref.fValidCuts),
fNCuts(ref.fNCuts),
fNMaxCuts(ref.fNMaxCuts),
fEventCutArray(new TString[fNMaxCuts]),
fPhotonCutArray(new TString[fNMaxCuts]),
fMesonCutArray(new TString[fNMaxCuts]),
fClusterCutArray(new TString[fNMaxCuts])
{
  for(Int_t i=0; i<fNMaxCuts; i++) {
	fEventCutArray[i] = ref.fEventCutArray[i];
	fPhotonCutArray[i] = ref.fPhotonCutArray[i];
	fMesonCutArray[i] = ref.fMesonCutArray[i];
	fClusterCutArray[i] = ref.fClusterCutArray[i];
  }
}

AliConversionCutHandler &AliConversionCutHandler::operator=(const AliConversionCutHandler &ref) {
  if(this != &ref){
    delete[] fEventCutArray;
    delete[] fPhotonCutArray;
    delete[] fMesonCutArray;
    delete[] fClusterCutArray;

    fValidCuts = ref.fValidCuts;
    fNCuts = ref.fNCuts;
    fNMaxCuts = ref.fNMaxCuts;
    fEventCutArray = new TString[fNMaxCuts];
    fPhotonCutArray = new TString[fNMaxCuts];
    fMesonCutArray = new TString[fNMaxCuts];
    fClusterCutArray = new TString[fNMaxCuts];
    for(Int_t i=0; i<fNMaxCuts; i++) {
      for(Int_t i=0; i<fNMaxCuts; i++) {
      fEventCutArray[i] = ref.fEventCutArray[i];
      fPhotonCutArray[i] = ref.fPhotonCutArray[i];
      fMesonCutArray[i] = ref.fMesonCutArray[i];
      fClusterCutArray[i] = ref.fClusterCutArray[i];
      }
    }
  }
  return *this;
}

AliConversionCutHandler::~AliConversionCutHandler() {
  delete[] fEventCutArray;
  delete[] fPhotonCutArray;
  delete[] fMesonCutArray;
  delete[] fClusterCutArray;
}

void AliConversionCutHandler::AddCut(TString eventCut, TString photonCut, TString mesonCut){
  if(fNCuts>=fNMaxCuts) {
    std::cout << "ERROR in CutHandlerConv: Exceeded maximum number of cuts!" << std::endl;
    fValidCuts = false;
    return;
  }
  if(eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 ) {  // Markus: These numbers should in principle be handled by the cut classes themselves
    std::cout << "ERROR in CutHandlerConv: Incorrect length of cut string!" << std::endl;
    fValidCuts = false;
    return;
  }
  fEventCutArray[fNCuts]=eventCut;
  fPhotonCutArray[fNCuts]=photonCut;
  fMesonCutArray[fNCuts]=mesonCut;
  fClusterCutArray[fNCuts]="";
  fNCuts++;
  return;
}

void AliConversionCutHandler::AddCut(TString eventCut, TString photonCut, TString mesonCut, TString clusterCut){
  if(fNCuts>=fNMaxCuts) {
    std::cout << "ERROR in CutHandlerConv: Exceeded maximum number of cuts!" << std::endl;
    fValidCuts = false;
    return;
  }
  if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 || clusterCut.Length()!=19 ){
    std::cout << "ERROR in CutHandlerConv: Incorrect length of cut string!" << std::endl;
    fValidCuts = false; return;
  }
  fEventCutArray[fNCuts]=eventCut;
  fPhotonCutArray[fNCuts]=photonCut;
  fMesonCutArray[fNCuts]=mesonCut;
  fClusterCutArray[fNCuts]=clusterCut;
  fNCuts++;
  return;
}

TString AliConversionCutHandler::GetEventCut(Int_t i) const {
  if(fValidCuts&&i<fNMaxCuts&&i>=0)
    return fEventCutArray[i];
  else{
    std::cout << "ERROR in CutHandlerConv: GetEventCut wrong index i" << std::endl;
    return "";
  }
}

TString AliConversionCutHandler::GetPhotonCut(Int_t i) const {
  if(fValidCuts&&i<fNMaxCuts&&i>=0)
    return fPhotonCutArray[i];
  else {
    std::cout << "ERROR in CutHandlerConv: GetPhotonCut wrong index i" << std::endl;
    return "";
  }
}

TString AliConversionCutHandler::GetMesonCut(Int_t i) const {
  if(fValidCuts&&i<fNMaxCuts&&i>=0)
    return fMesonCutArray[i];
  else {
    std::cout << "ERROR in CutHandlerConv: GetMesonCut wrong index i" << std::endl;
    return "";
  }
}

TString AliConversionCutHandler::GetClusterCut(Int_t i) const {
  if(fValidCuts&&i<fNMaxCuts&&i>=0)
    return fClusterCutArray[i];
  else {
    std::cout << "ERROR in CutHandlerConv: GetClusterCut wrong index i" << std::endl;
    return "";
  }
}
