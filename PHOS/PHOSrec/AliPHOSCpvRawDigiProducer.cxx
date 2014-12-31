/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
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
#include "TClonesArray.h"

// --- AliRoot header files ---
#include "AliPHOSCpvRawDigiProducer.h"
#include "AliPHOSDigit.h"
#include "AliPHOSCpvRawStream.h"
#include "AliLog.h"
#include<iostream>
using namespace std;

ClassImp(AliPHOSCpvRawDigiProducer);

//--------------------------------------------------------------------------------------
AliPHOSCpvRawDigiProducer::AliPHOSCpvRawDigiProducer():
  TObject(),
  fTurbo(kFALSE),
  fCpvMinE(10.),
  fRawStream(0),
  fhErrors(0),
  fPedFilesRLoaded(kFALSE)
{
  CreateErrHist();
  // create a 2d array to store the pedestals                                        
  for (Int_t iDDL=0;iDDL<AliPHOSCpvParam::kNDDL;iDDL++){
    ped[0][iDDL] = new Int_t *[AliPHOSCpvParam::kPadPcX];
    ped[1][iDDL] = new Int_t *[AliPHOSCpvParam::kPadPcX];
    for(Int_t ix=0; ix<AliPHOSCpvParam::kPadPcX; ix++) {
      ped[0][iDDL][ix] = new Int_t [AliPHOSCpvParam::kPadPcY];
      ped[1][iDDL][ix] = new Int_t [AliPHOSCpvParam::kPadPcY];
    }
  }
}
//-------------------------------------------------------------------------------------
AliPHOSCpvRawDigiProducer::AliPHOSCpvRawDigiProducer(AliRawReader *& rawReader):
  TObject(),
  fTurbo(kFALSE),
  fCpvMinE(10.),
  fRawStream(0),
  fhErrors(0),
  fPedFilesRLoaded(kFALSE)
{
  fRawStream = new AliPHOSCpvRawStream(rawReader);
  fRawStream->SetTurbo(fTurbo);
  CreateErrHist();
  // create a 2d array to store the pedestals                                               
  for (Int_t iDDL=0;iDDL<AliPHOSCpvParam::kNDDL;iDDL++) {
    ped[0][iDDL] = new Int_t *[AliPHOSCpvParam::kPadPcX];
    ped[1][iDDL] = new Int_t *[AliPHOSCpvParam::kPadPcX];
    for(Int_t ix=0; ix<AliPHOSCpvParam::kPadPcX; ix++) {
      ped[0][iDDL][ix] = new Int_t [AliPHOSCpvParam::kPadPcY];
      ped[1][iDDL][ix] = new Int_t [AliPHOSCpvParam::kPadPcY];
    }
  }
}
//--------------------------------------------------------------------------------------
AliPHOSCpvRawDigiProducer::~AliPHOSCpvRawDigiProducer()
{
  if(fRawStream) delete fRawStream;
  if(fhErrors) delete fhErrors; 
  for(Int_t iDDL = 0;iDDL<AliPHOSCpvParam::kNDDL;iDDL++) {
    for(Int_t ix=0; ix<AliPHOSCpvParam::kPadPcX; ix++) { 
      delete [] ped[0][iDDL][ix];
      delete [] ped[1][iDDL][ix];
    }
    delete [] ped[0][iDDL];
    delete [] ped[1][iDDL];
  }
}
//--------------------------------------------------------------------------------------
Bool_t AliPHOSCpvRawDigiProducer::LoadPedFiles() {
  // read pedestals from file                                             
  for(Int_t iDDL = 0;iDDL<AliPHOSCpvParam::kNDDL;iDDL++)
    for(Int_t iCC=0; iCC<AliPHOSCpvParam::kNRows; iCC++) {
      FILE * pedFile;
      pedFile = fopen(Form("thr%d_%02d.dat",iDDL,iCC),"r");
      if(!pedFile) {
	Printf("AliPHOSCpvRawDigiProducer::LoadPedFiles: Error, file thr%d_%02d.dat could not be open",iDDL,iCC);
	continue;
	//return kFALSE;
      }
      Int_t i3g = 0, iPad = 0;
      Int_t lineCnt = 0;
      while(!feof(pedFile)) {
	Int_t abs = AliPHOSCpvParam::Abs(iDDL,iCC,i3g,iPad);
	if(iPad<48&&i3g<10){
	  if(AliPHOSCpvParam::A2DDL(abs)!=iDDL)
	    cout<<"AliPHOSCpvRawDigiProducer::LoadPedFiles(): wrong connection table! abs = "
		<<abs<<", DDL = "<<iDDL<<", A2DDL = "<<AliPHOSCpvParam::A2DDL(abs)<<endl;
	  if(AliPHOSCpvParam::A2CC(abs)!=iCC)
	    cout<<"AliPHOSCpvRawDigiProducer::LoadPedFiles(): wrong connection table! abs = "
		<<abs<<", CC = "<< iCC <<", A2CC = "<<AliPHOSCpvParam::A2CC(abs)<<endl;
	  if(AliPHOSCpvParam::A23G(abs)!=i3g)
	    cout<<"AliPHOSCpvRawDigiProducer::LoadPedFiles(): wrong connection table! abs = "
		<<abs<<", 3G = "<< i3g <<", A23G = "<<AliPHOSCpvParam::A23G(abs)<<endl;
	  if(AliPHOSCpvParam::A2Pad(abs)!=iPad)
	    cout<<"AliPHOSCpvRawDigiProducer::LoadPedFiles(): wrong connection table! abs = "
		<<abs<<", Pad = "<< iPad <<", A2Pad = "<<AliPHOSCpvParam::A2Pad(abs)<<endl;
	}
	Int_t thr;
	fscanf(pedFile,"%x",&thr);
	if(AliPHOSCpvParam::IsValidAbs(abs)) {
	  Int_t s = thr & 0x1ff;
	  Int_t p = thr >> 9;
	  ped[0][iDDL][AliPHOSCpvParam::A2X(abs)][AliPHOSCpvParam::A2Y(abs)] = p-s;
	  ped[1][iDDL][AliPHOSCpvParam::A2X(abs)][AliPHOSCpvParam::A2Y(abs)] = s;
	  int testAbs = AliPHOSCpvParam::XY2A(iDDL,AliPHOSCpvParam::A2X(abs),AliPHOSCpvParam::A2Y(abs));
	  if(abs!=testAbs)
	    cout<<"AliPHOSCpvRawDigiProducer::LoadPedFiles(): wrong connection table! abs = "
		<<abs<<", testAbs = "<<testAbs<<endl;
	  //Printf("ped[%d][%d] = %d, pad = %d, abs = %d",AliPHOSCpvParam::A2X(abs),AliPHOSCpvParam::A2Y(abs), p + s, iPad, abs);
	}
	iPad++;
	if(iPad == 64) {iPad = 0; i3g++;}
	lineCnt++;
      }
      if(lineCnt < AliPHOSCpvParam::kN3GAdd * 64) return kFALSE;
      fclose(pedFile);
    }
  fPedFilesRLoaded = kTRUE;
  return kTRUE;
}
//--------------------------------------------------------------------------------------
void AliPHOSCpvRawDigiProducer::SetTurbo(Bool_t turbo) 
{
  fTurbo = turbo;
  if(fRawStream) fRawStream->SetTurbo(fTurbo);
}
//--------------------------------------------------------------------------------------
Bool_t AliPHOSCpvRawDigiProducer::LoadNewEvent(AliRawReader *& rawReader)
{
  if(fRawStream) delete fRawStream;
  fRawStream = new AliPHOSCpvRawStream(rawReader);
  if(fRawStream) {
    fRawStream->SetTurbo(fTurbo);
    return kTRUE;
  }
  fhErrors->Fill(0);
  return kFALSE;
}
//--------------------------------------------------------------------------------------
void AliPHOSCpvRawDigiProducer::MakeDigits(TClonesArray *& digits) const
{
  // returns histogram of error types

  if(digits) 
    digits->Clear();
  digits = new TClonesArray("AliPHOSDigit", AliPHOSCpvParam::kNDDL * AliPHOSCpvParam::kNRows * AliPHOSCpvParam::kN3GAdd * AliPHOSCpvParam::kNPadAdd);
  Int_t iDigit = 0;
  while(fRawStream->Next()) {
    for(Int_t iPad=0;iPad<fRawStream->GetNPads();iPad++) {
      Int_t charge = fRawStream->GetChargeArray()[iPad];
      Int_t aPad   = fRawStream->GetPadArray()[iPad];
      //cout<<"AliPHOSCpvRawDigiProducer::MakeDigits(): I've got pad "<<aPad<< "with amplitude "<<charge<<endl;
      if(fPedFilesRLoaded) {
	Int_t ix   = AliPHOSCpvParam::A2X(aPad);
	Int_t iy   = AliPHOSCpvParam::A2Y(aPad);
	Int_t iddl = AliPHOSCpvParam::A2DDL(aPad);
	if (charge>ped[0][iddl][ix][iy]+ped[1][iddl][ix][iy]){
	  charge -=ped[0][iddl][ix][iy];
	}
	else charge=0;
	if(charge < fCpvMinE) charge = 0;
      }
      // if(charge) new((*digits)[iDigit++]) AliPHOSDigit(AliPHOSCpvParam::A2X(aPad),AliPHOSCpvParam::A2Y(aPad),charge);
      // Check what is aPad! YK 31.12.2014
      if(charge) new((*digits)[iDigit++]) AliPHOSDigit(-1,aPad,charge,0);
      
    }
  } // while(fRawStream->Next())
  //cout<<"AliPHOSCpvRawDigiProducer::MakeDigits(): I've created "<<iDigit<<" digits."<<endl;
  // fill histogram of errors
  for(Int_t iDDL=0; iDDL<AliPHOSCpvParam::kNDDL; iDDL++) {
    Int_t nErrors = AliPHOSCpvRawStream::GetNErrors();
    for(Int_t iType=0; iType<nErrors; iType++) { // iType - type of error
      fhErrors -> Fill(iType+1,fRawStream -> GetErrors(iDDL,iType));
    }
  }
}
//--------------------------------------------------------------------------------------
void AliPHOSCpvRawDigiProducer::CreateErrHist()
{
  Int_t nErrors = AliPHOSCpvRawStream::GetNErrors();
  const char * errNames[nErrors];
  for(Int_t i=0; i<nErrors; i++) {
    errNames[i] = AliPHOSCpvRawStream::GetErrName(i);
  }
  fhErrors = new TH1I("errorTypes","Errors occured during processing",nErrors+1,0,nErrors+1);
  TAxis* x = fhErrors->GetXaxis();
  x->SetBinLabel(1, "Can't get event");
  for(Int_t i=0; i<nErrors; i++) {
    x->SetBinLabel(i+2,errNames[i]);
  }

}
//--------------------------------------------------------------------------------------
