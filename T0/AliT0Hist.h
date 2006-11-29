#ifndef ALIT0HIST_H
#define ALIT0HIST_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TH1F.h"
#include "TH2F.h"
#include "TArrayI.h"
#include "AliRunLoader.h"
#include "/home/alla/AliRoot/verynew/RAW/AliRawReader.h"

class AliT0Hist: public TObject
{
 public:
  AliT0Hist();
  virtual ~AliT0Hist();
  void  FillHist(AliRunLoader* runLoader,
		 AliRawReader* rawReader) const;
   TH1F * hTimeR() {return fhTimeR;} 
   TH1F * hTimeL() {return fhTimeL;}
   TH1F *hADCR()   {return fhADCL;}
   TH1F *hADCL() {return fhADCR;} 
   TH1F * hBestTimeR() {return fhBestTimeR;} 
   TH1F * hBestTimeL() {return fhBestTimeL;}
   TH1F * hTimeDet() {return fhTimeDet;} 
   TH1F * hADCDet() {return fhADCDet;}
   TH1F * hTimeDiff() {return fhTimeDiff;} 
   TH1F * hMeanTime() {return fhMeanTime;}
   TH1F * hT0detL() {return   fhT0detL;}
   TH1F * hT0detR() {return   fhT0detR;}
   TH2F * hTimevsADCR() {return fhTimevsADCR;}
   TH2F * hTimevsADCL() {return fhTimevsADCL;}
   
 private:
   TH1F *fhTimeR;
   TH1F *fhTimeL;
   TH1F *fhADCL;
   TH1F *fhADCR;
   TH1F *fhBestTimeR;
   TH1F *fhBestTimeL;
   TH1F *fhADCDet;
   TH1F *fhTimeDet;
   TH1F *fhMeanTime;
   TH1F *fhTimeDiff;
   TH1F * fhT0detL;
   TH1F * fhT0detR;
   TH2F  *fhTimevsADCR;
   TH2F  *fhTimevsADCL;
   
   ClassDef(AliT0Hist, 0)   // class for the T0 reconstruction

};

#endif
