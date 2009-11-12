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

///////////////////////////////////////////////////////////////////////
//                                                                   //
//  Produces the data needed to calculate the ZDC quality assurance. //
//  QA objects are 1 & 2 Dimensional histograms.                     //
//  author: C. Oppedisano                                            //
//                                                                   //
///////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH2F.h>
#include <TProfile.h>
#include <Riostream.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAChecker.h"
#include "AliZDCReco.h"
#include "AliRawReader.h"
#include "AliZDCQADataMakerRec.h"
#include "AliZDCRawStream.h"
#include "AliZDCDigit.h"
#include "AliESDZDC.h"
#include "AliESDEvent.h"

ClassImp(AliZDCQADataMakerRec)
           
//____________________________________________________________________________ 
  AliZDCQADataMakerRec::AliZDCQADataMakerRec() : 
  AliQADataMakerRec(AliQAv1::GetDetName(AliQAv1::kZDC), "ZDC Quality Assurance Data Maker"), 
  fDigit(0)
{
  // ctor
}

//____________________________________________________________________________ 
AliZDCQADataMakerRec::AliZDCQADataMakerRec(const AliZDCQADataMakerRec& qadm) :
  AliQADataMakerRec(),      
  fDigit(0)
{
  //copy ctor 
  SetName((const char*)qadm.GetName()); 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliZDCQADataMakerRec& AliZDCQADataMakerRec::operator = (const AliZDCQADataMakerRec& qadm )
{
  // Equal operator.
  this->~AliZDCQADataMakerRec();
  new(this) AliZDCQADataMakerRec(qadm);
  return *this;
}

//____________________________________________________________________________ 
void AliZDCQADataMakerRec::InitDigits()
{
  // create Digits histograms in Digits subdir
  //
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  // ------------------- HIGH GAIN CHAIN ---------------------------
  TH1F * hDigitPMCZNC = new TH1F("hDigitPMCZNC", "Digit common ZNC PMT;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPMCZNA = new TH1F("hDigitPMCZNA", "Digit common ZNA PMT;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPMCZPC = new TH1F("hDigitPMCZPC", "Digit common ZPC PMT;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPMCZPA = new TH1F("hDigitPMCZPA", "Digit common ZPA PMT;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPMZEM1 = new TH1F("hDigitPMZEM1", "Digit ZEM1 PMT;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPMZEM2 = new TH1F("hDigitPMZEM2", "Digit ZEM2 PMT;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2DigitsList(hDigitPMCZNC, 0 , !expert, image);
  Add2DigitsList(hDigitPMCZNA, 1 , !expert, image);
  Add2DigitsList(hDigitPMCZPC, 2, !expert, image);
  Add2DigitsList(hDigitPMCZPA, 3, !expert, image);
  Add2DigitsList(hDigitPMZEM1, 4, !expert, image);
  Add2DigitsList(hDigitPMZEM2, 5, !expert, image);
  // 
  // ------------------- LOW GAIN CHAIN ---------------------------
  TH1F * hDigitPMCZNClg = new TH1F("hDigitPMCZNClg", "Digit common lg ZNC PMT;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigitPMCZNAlg = new TH1F("hDigitPMCZNAlg", "Digit common lg ZNA PMT;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigitPMCZPClg = new TH1F("hDigitPMCZPClg", "Digit common lg ZPC PMT;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigitPMCZPAlg = new TH1F("hDigitPMCZPAlg", "Digit common lg ZPA PMT;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigitPMZEM1lg = new TH1F("hDigitPMZEM1lg", "Digit ZEM1 PMT lg;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigitPMZEM2lg = new TH1F("hDigitPMZEM2lg", "Digit ZEM2 PMT lg;Amplitude [ADC counts];Counts",100, 0., 4000.);
  Add2DigitsList(hDigitPMCZNClg, 6, expert, !image);
  Add2DigitsList(hDigitPMCZNAlg, 7, expert, !image);
  Add2DigitsList(hDigitPMCZPClg, 8, expert, !image);
  Add2DigitsList(hDigitPMCZPAlg, 9, expert, !image);
  Add2DigitsList(hDigitPMZEM1lg, 10, !expert, image);
  Add2DigitsList(hDigitPMZEM2lg, 11, !expert, image);
  //
  // ------------------- HIGH GAIN CHAIN ---------------------------
  TH1F * hDigitPM1ZNC = new TH1F("hDigitPM1ZNC", "Digit PM1 ZNC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM2ZNC = new TH1F("hDigitPM2ZNC", "Digit PM2 ZNC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM3ZNC = new TH1F("hDigitPM3ZNC", "Digit PM3 ZNC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM4ZNC = new TH1F("hDigitPM4ZNC", "Digit PM4 ZNC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2DigitsList(hDigitPM1ZNC, 12, !expert, image);
  Add2DigitsList(hDigitPM2ZNC, 13, !expert, image);
  Add2DigitsList(hDigitPM3ZNC, 14, !expert, image);
  Add2DigitsList(hDigitPM4ZNC, 15, !expert, image);
  // 
  TH1F * hDigitPM1ZPC = new TH1F("hDigitPM1ZPC", "Digit PM1 ZPC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM2ZPC = new TH1F("hDigitPM2ZPC", "Digit PM2 ZPC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM3ZPC = new TH1F("hDigitPM3ZPC", "Digit PM3 ZPC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM4ZPC = new TH1F("hDigitPM4ZPC", "Digit PM4 ZPC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2DigitsList(hDigitPM1ZPC, 16, !expert, image);
  Add2DigitsList(hDigitPM2ZPC, 17, !expert, image);
  Add2DigitsList(hDigitPM3ZPC, 18, !expert, image);
  Add2DigitsList(hDigitPM4ZPC, 19, !expert, image);
  //
  TH1F * hDigitPM1ZNA = new TH1F("hDigitPM1ZNA", "Digit PM1 ZNA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM2ZNA = new TH1F("hDigitPM2ZNA", "Digit PM2 ZNA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM3ZNA = new TH1F("hDigitPM3ZNA", "Digit PM3 ZNA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM4ZNA = new TH1F("hDigitPM4ZNA", "Digit PM4 ZNA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2DigitsList(hDigitPM1ZNA, 20, !expert, image);
  Add2DigitsList(hDigitPM2ZNA, 21, !expert, image);
  Add2DigitsList(hDigitPM3ZNA, 22, !expert, image);
  Add2DigitsList(hDigitPM4ZNA, 23, !expert, image);
  // 
  TH1F * hDigitPM1ZPA = new TH1F("hDigitPM1ZPA", "Digit PM1 ZPA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM2ZPA = new TH1F("hDigitPM2ZPA", "Digit PM2 ZPA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM3ZPA = new TH1F("hDigitPM3ZPA", "Digit PM3 ZPA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM4ZPA = new TH1F("hDigitPM4ZPA", "Digit PM4 ZPA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2DigitsList(hDigitPM1ZPA, 24, !expert, image);
  Add2DigitsList(hDigitPM2ZPA, 25, !expert, image);
  Add2DigitsList(hDigitPM3ZPA, 26, !expert, image);
  Add2DigitsList(hDigitPM4ZPA, 27, !expert, image);
  //
  // ------------------- LOW GAIN CHAIN ---------------------------
  TH1F * hDigitPM1ZNClg = new TH1F("hDigitPM1ZNClg", "Digit PM1 ZNC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM2ZNClg = new TH1F("hDigitPM2ZNClg", "Digit PM2 ZNC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM3ZNClg = new TH1F("hDigitPM3ZNClg", "Digit PM3 ZNC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM4ZNClg = new TH1F("hDigitPM4ZNClg", "Digit PM4 ZNC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2DigitsList(hDigitPM1ZNClg, 28, !expert, image);
  Add2DigitsList(hDigitPM2ZNClg, 29, !expert, image);
  Add2DigitsList(hDigitPM3ZNClg, 30, !expert, image);
  Add2DigitsList(hDigitPM4ZNClg, 31, !expert, image);
  // 
  TH1F * hDigitPM1ZPClg = new TH1F("hDigitPM1ZPClg", "Digit PM1 ZPC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM2ZPClg = new TH1F("hDigitPM2ZPClg", "Digit PM2 ZPC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM3ZPClg = new TH1F("hDigitPM3ZPClg", "Digit PM3 ZPC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM4ZPClg = new TH1F("hDigitPM4ZPClg", "Digit PM4 ZPC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2DigitsList(hDigitPM1ZPClg, 32, !expert, image);
  Add2DigitsList(hDigitPM2ZPClg, 33, !expert, image);
  Add2DigitsList(hDigitPM3ZPClg, 34, !expert, image);
  Add2DigitsList(hDigitPM4ZPClg, 35, !expert, image);
  //
  TH1F * hDigitPM1ZNAlg = new TH1F("hDigitPM1ZNAlg", "Digit PM1 ZNA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM2ZNAlg = new TH1F("hDigitPM2ZNAlg", "Digit PM2 ZNA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM3ZNAlg = new TH1F("hDigitPM3ZNAlg", "Digit PM3 ZNA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM4ZNAlg = new TH1F("hDigitPM4ZNAlg", "Digit PM4 ZNA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2DigitsList(hDigitPM1ZNAlg, 36, !expert, image);
  Add2DigitsList(hDigitPM2ZNAlg, 37, !expert, image);
  Add2DigitsList(hDigitPM3ZNAlg, 38, !expert, image);
  Add2DigitsList(hDigitPM4ZNAlg, 39, !expert, image);
  // 
  TH1F * hDigitPM1ZPAlg = new TH1F("hDigitPM1ZPAlg", "Digit PM1 ZPA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM2ZPAlg = new TH1F("hDigitPM2ZPAlg", "Digit PM2 ZPA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM3ZPAlg = new TH1F("hDigitPM3ZPAlg", "Digit PM3 ZPA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hDigitPM4ZPAlg = new TH1F("hDigitPM4ZPAlg", "Digit PM4 ZPA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2DigitsList(hDigitPM1ZPAlg, 40, !expert, image);
  Add2DigitsList(hDigitPM2ZPAlg, 41, !expert, image);
  Add2DigitsList(hDigitPM3ZPAlg, 42, !expert, image);
  Add2DigitsList(hDigitPM4ZPAlg, 43, !expert, image);
  
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::InitRecPoints()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert = kTRUE ; 
  const Bool_t image  = kTRUE ; 
  //
  // ------------------- HIGH GAIN CHAIN ---------------------------
  TH1F * hRecZNCTot = new TH1F("hRecZNCTot", "Rec signal in ZNC;Amplitude [ADC counts];Counts", 100, 0., 6000.);
  TH1F * hRecZNATot = new TH1F("hRecZNATot", "Rec signal in ZNA;Amplitude [ADC counts];Counts", 100, 0., 6000.);
  TH1F * hRecZPCTot = new TH1F("hRecZPCTot", "Rec signal in ZPC;Amplitude [ADC counts];Counts", 100, 0., 10000.);
  TH1F * hRecZPATot = new TH1F("hRecZPATot", "Rec signal in ZPA;Amplitude [ADC counts];Counts", 100, 0., 10000.);
  Add2RecPointsList(hRecZNCTot, 0, expert, !image);
  Add2RecPointsList(hRecZNATot, 1, expert, !image);
  Add2RecPointsList(hRecZPCTot, 2, expert, !image);
  Add2RecPointsList(hRecZPATot, 3, expert, !image);
  //
  TH1F * hRecSumQZNC = new TH1F("hRecSumQZNC", "Rec summed 4 ZNC quadrants;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hRecSumQZNA = new TH1F("hRecSumQZNA", "Rec summed 4 ZNA quadrants;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hRecSumQZPC = new TH1F("hRecSumQZPC", "Rec summed 4 ZPC quadrants;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hRecSumQZPA = new TH1F("hRecSumQZPA", "Rec summed 4 ZPA quadrants;Amplitude [ADC counts];Counts",100, 0., 4000.);
  Add2RecPointsList(hRecSumQZNC, 4, expert, !image);
  Add2RecPointsList(hRecSumQZNA, 5, expert, !image);
  Add2RecPointsList(hRecSumQZPC, 6, expert, !image);
  Add2RecPointsList(hRecSumQZPA, 7, expert, !image);
  //
  TH1F * hRecPMCZNC = new TH1F("hRecPMCZNC", "Rec common ZNC PMT;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hRecPMCZNA = new TH1F("hRecPMCZNA", "Rec common ZNA PMT;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hRecPMCZPC = new TH1F("hRecPMCZPC", "Rec common ZPC PMT;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hRecPMCZPA = new TH1F("hRecPMCZPA", "Rec common ZPA PMT;Amplitude [ADC counts];Counts",100, 0., 4000.);
  Add2RecPointsList(hRecPMCZNC, 8 , !expert, image);
  Add2RecPointsList(hRecPMCZNA, 9 , !expert, image);
  Add2RecPointsList(hRecPMCZPC, 10, !expert, image);
  Add2RecPointsList(hRecPMCZPA, 11, !expert, image); 
}


//____________________________________________________________________________
void AliZDCQADataMakerRec::InitRaws()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  //
  // ------------------- HIGH GAIN CHAIN ---------------------------
  TH1F * hRawPMCZNC = new TH1F("hRawPMCZNC", "Raw common ZNC PMT;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPMCZNA = new TH1F("hRawPMCZNA", "Raw common ZNA PMT;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPMCZPC = new TH1F("hRawPMCZPC", "Raw common ZPC PMT;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPMCZPA = new TH1F("hRawPMCZPA", "Raw common ZPA PMT;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPMZEM1 = new TH1F("hRawPMZEM1", "Raw ZEM1 PMT;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPMZEM2 = new TH1F("hRawPMZEM2", "Raw ZEM2 PMT;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2RawsList(hRawPMCZNC, 0 , !expert, image);
  Add2RawsList(hRawPMCZNA, 1 , !expert, image);
  Add2RawsList(hRawPMCZPC, 2, !expert, image);
  Add2RawsList(hRawPMCZPA, 3, !expert, image);
  Add2RawsList(hRawPMZEM1, 4, !expert, image);
  Add2RawsList(hRawPMZEM2, 5, !expert, image);
  // 
  // ------------------- LOW GAIN CHAIN ---------------------------
  TH1F * hRawPMCZNClg = new TH1F("hRawPMCZNClg", "Raw common lg ZNC PMT;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hRawPMCZNAlg = new TH1F("hRawPMCZNAlg", "Raw common lg ZNA PMT;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hRawPMCZPClg = new TH1F("hRawPMCZPClg", "Raw common lg ZPC PMT;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hRawPMCZPAlg = new TH1F("hRawPMCZPAlg", "Raw common lg ZPA PMT;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hRawPMZEM1lg = new TH1F("hRawPMZEM1lg", "Raw ZEM1 PMT lg;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hRawPMZEM2lg = new TH1F("hRawPMZEM2lg", "Raw ZEM2 PMT lg;Amplitude [ADC counts];Counts",100, 0., 4000.);
  Add2RawsList(hRawPMCZNClg, 6, expert, !image);
  Add2RawsList(hRawPMCZNAlg, 7, expert, !image);
  Add2RawsList(hRawPMCZPClg, 8, expert, !image);
  Add2RawsList(hRawPMCZPAlg, 9, expert, !image);
  Add2RawsList(hRawPMZEM1lg, 10, !expert, image);
  Add2RawsList(hRawPMZEM2lg, 11, !expert, image);
  //
  // ------------------- HIGH GAIN CHAIN ---------------------------
  TH1F * hRawPM1ZNC = new TH1F("hRawPM1ZNC", "Raw PM1 ZNC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM2ZNC = new TH1F("hRawPM2ZNC", "Raw PM2 ZNC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM3ZNC = new TH1F("hRawPM3ZNC", "Raw PM3 ZNC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM4ZNC = new TH1F("hRawPM4ZNC", "Raw PM4 ZNC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2RawsList(hRawPM1ZNC, 12, !expert, image);
  Add2RawsList(hRawPM2ZNC, 13, !expert, image);
  Add2RawsList(hRawPM3ZNC, 14, !expert, image);
  Add2RawsList(hRawPM4ZNC, 15, !expert, image);
  // 
  TH1F * hRawPM1ZPC = new TH1F("hRawPM1ZPC", "Raw PM1 ZPC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM2ZPC = new TH1F("hRawPM2ZPC", "Raw PM2 ZPC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM3ZPC = new TH1F("hRawPM3ZPC", "Raw PM3 ZPC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM4ZPC = new TH1F("hRawPM4ZPC", "Raw PM4 ZPC;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2RawsList(hRawPM1ZPC, 16, !expert, image);
  Add2RawsList(hRawPM2ZPC, 17, !expert, image);
  Add2RawsList(hRawPM3ZPC, 18, !expert, image);
  Add2RawsList(hRawPM4ZPC, 19, !expert, image);
  //
  TH1F * hRawPM1ZNA = new TH1F("hRawPM1ZNA", "Raw PM1 ZNA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM2ZNA = new TH1F("hRawPM2ZNA", "Raw PM2 ZNA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM3ZNA = new TH1F("hRawPM3ZNA", "Raw PM3 ZNA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM4ZNA = new TH1F("hRawPM4ZNA", "Raw PM4 ZNA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2RawsList(hRawPM1ZNA, 20, !expert, image);
  Add2RawsList(hRawPM2ZNA, 21, !expert, image);
  Add2RawsList(hRawPM3ZNA, 22, !expert, image);
  Add2RawsList(hRawPM4ZNA, 23, !expert, image);
  // 
  TH1F * hRawPM1ZPA = new TH1F("hRawPM1ZPA", "Raw PM1 ZPA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM2ZPA = new TH1F("hRawPM2ZPA", "Raw PM2 ZPA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM3ZPA = new TH1F("hRawPM3ZPA", "Raw PM3 ZPA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM4ZPA = new TH1F("hRawPM4ZPA", "Raw PM4 ZPA;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2RawsList(hRawPM1ZPA, 24, !expert, image);
  Add2RawsList(hRawPM2ZPA, 25, !expert, image);
  Add2RawsList(hRawPM3ZPA, 26, !expert, image);
  Add2RawsList(hRawPM4ZPA, 27, !expert, image);
  //
  // ------------------- LOW GAIN CHAIN ---------------------------
  TH1F * hRawPM1ZNClg = new TH1F("hRawPM1ZNClg", "Raw PM1 ZNC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM2ZNClg = new TH1F("hRawPM2ZNClg", "Raw PM2 ZNC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM3ZNClg = new TH1F("hRawPM3ZNClg", "Raw PM3 ZNC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM4ZNClg = new TH1F("hRawPM4ZNClg", "Raw PM4 ZNC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2RawsList(hRawPM1ZNClg, 28, !expert, image);
  Add2RawsList(hRawPM2ZNClg, 29, !expert, image);
  Add2RawsList(hRawPM3ZNClg, 30, !expert, image);
  Add2RawsList(hRawPM4ZNClg, 31, !expert, image);
  // 
  TH1F * hRawPM1ZPClg = new TH1F("hRawPM1ZPClg", "Raw PM1 ZPC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM2ZPClg = new TH1F("hRawPM2ZPClg", "Raw PM2 ZPC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM3ZPClg = new TH1F("hRawPM3ZPClg", "Raw PM3 ZPC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM4ZPClg = new TH1F("hRawPM4ZPClg", "Raw PM4 ZPC lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2RawsList(hRawPM1ZPClg, 32, !expert, image);
  Add2RawsList(hRawPM2ZPClg, 33, !expert, image);
  Add2RawsList(hRawPM3ZPClg, 34, !expert, image);
  Add2RawsList(hRawPM4ZPClg, 35, !expert, image);
  //
  TH1F * hRawPM1ZNAlg = new TH1F("hRawPM1ZNAlg", "Raw PM1 ZNA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM2ZNAlg = new TH1F("hRawPM2ZNAlg", "Raw PM2 ZNA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM3ZNAlg = new TH1F("hRawPM3ZNAlg", "Raw PM3 ZNA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM4ZNAlg = new TH1F("hRawPM4ZNAlg", "Raw PM4 ZNA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2RawsList(hRawPM1ZNAlg, 36, !expert, image);
  Add2RawsList(hRawPM2ZNAlg, 37, !expert, image);
  Add2RawsList(hRawPM3ZNAlg, 38, !expert, image);
  Add2RawsList(hRawPM4ZNAlg, 39, !expert, image);
  // 
  TH1F * hRawPM1ZPAlg = new TH1F("hRawPM1ZPAlg", "Raw PM1 ZPA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM2ZPAlg = new TH1F("hRawPM2ZPAlg", "Raw PM2 ZPA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM3ZPAlg = new TH1F("hRawPM3ZPAlg", "Raw PM3 ZPA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  TH1F * hRawPM4ZPAlg = new TH1F("hRawPM4ZPAlg", "Raw PM4 ZPA lg;Amplitude [ADC counts];Counts",100, 0., 1000.);
  Add2RawsList(hRawPM1ZPAlg, 40, !expert, image);
  Add2RawsList(hRawPM2ZPAlg, 41, !expert, image);
  Add2RawsList(hRawPM3ZPAlg, 42, !expert, image);
  Add2RawsList(hRawPM4ZPAlg, 43, !expert, image);
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::InitESDs()
{
  //Booking ESDs histograms
  //
  const Bool_t expert = kTRUE ; 
  const Bool_t image  = kTRUE ; 
  
  TH2F * hZNC  = new TH2F("hZNC", "Centroid in ZNC", 100, -5.,5.,100,-5.,5.);
  TH2F * hZNA  = new TH2F("hZNA", "Centroid in ZNA", 100, -5.,5.,100,-5.,5.);
  Add2ESDsList(hZNC, 0, !expert, image);
  Add2ESDsList(hZNA, 1, !expert, image);
  //
  // ------------------- HIGH GAIN CHAIN ---------------------------
  TH1F * hESDZNCTot = new TH1F("hESDZNCTot", "Energy in ZNC", 100, 0., 6000.);
  TH1F * hESDZNATot = new TH1F("hESDZNATot", "Energy in ZNA", 100, 0., 6000.);
  TH1F * hESDZPCTot = new TH1F("hESDZPCTot", "Energy in ZPC", 100, 0., 10000.);
  TH1F * hESDZPATot = new TH1F("hESDZPATot", "Energy in ZPA", 100, 0., 10000.);
  Add2ESDsList(hESDZNCTot, 2, expert, !image);
  Add2ESDsList(hESDZNATot, 3, expert, !image);
  Add2ESDsList(hESDZPCTot, 4, expert, !image);
  Add2ESDsList(hESDZPATot, 5, expert, !image);
  //
  TH1F * hESDSumQZNC = new TH1F("hESDSumQZNC", "Sum of 4 ZNC energy",100, 0., 4000.);
  TH1F * hESDSumQZNA = new TH1F("hESDSumQZNA", "Sum of 4 ZNA energy",100, 0., 4000.);
  TH1F * hESDSumQZPC = new TH1F("hESDSumQZPC", "Sum of 4 ZPC energy",100, 0., 4000.);
  TH1F * hESDSumQZPA = new TH1F("hESDSumQZPA", "Sum of 4 ZPA energy",100, 0., 4000.);
  Add2ESDsList(hESDSumQZNC, 6, expert, !image);
  Add2ESDsList(hESDSumQZNA, 7, expert, !image);
  Add2ESDsList(hESDSumQZPC, 8, expert, !image);
  Add2ESDsList(hESDSumQZPA, 9, expert, !image);
  //
  TH1F * hESDPMCZNC = new TH1F("hESDPMCZNC", "Energy in common ZNC PMT",100, 0., 4000.);
  TH1F * hESDPMCZNA = new TH1F("hESDPMCZNA", "Energy in common ZNA PMT",100, 0., 4000.);
  TH1F * hESDPMCZPC = new TH1F("hESDPMCZPC", "Energy in common ZPC PMT",100, 0., 4000.);
  TH1F * hESDPMCZPA = new TH1F("hESDPMCZPA", "Energy in common ZPA PMT",100, 0., 4000.);
  Add2ESDsList(hESDPMCZNC, 10, !expert, image);
  Add2ESDsList(hESDPMCZNA, 11, !expert, image);
  Add2ESDsList(hESDPMCZPC, 12, !expert, image);
  Add2ESDsList(hESDPMCZPA, 13, !expert, image);
  // 
  // ------------------- LOW GAIN CHAIN ---------------------------
  TH1F * hESDZNCTotlg = new TH1F("hESDZNCTotlg", "ESD lg signal in ZNC", 100, 0., 6000.);
  TH1F * hESDZNATotlg = new TH1F("hESDZNATotlg", "ESD lg signal in ZNA", 100, 0., 6000.);
  TH1F * hESDZPCTotlg = new TH1F("hESDZPCTotlg", "ESD lg signal in ZPC", 100, 0., 10000.);
  TH1F * hESDZPATotlg = new TH1F("hESDZPATotlg", "ESD lg signal in ZPA", 100, 0., 10000.);
  Add2ESDsList(hESDZNCTotlg, expert, !image);
  Add2ESDsList(hESDZNATotlg, expert, !image);
  Add2ESDsList(hESDZPCTotlg, expert, !image);
  Add2ESDsList(hESDZPATotlg, expert, !image);
  //
  TH1F * hESDSumQZNClg = new TH1F("hESDSumQZNClg", "Sum of 4 lg ZNC sectors",100, 0., 4000.);
  TH1F * hESDSumQZNAlg = new TH1F("hESDSumQZNAlg", "Sum of 4 lg ZNA sectors",100, 0., 4000.);
  TH1F * hESDSumQZPClg = new TH1F("hESDSumQZPClg", "Sum of 4 lg ZPC sectors",100, 0., 4000.);
  TH1F * hESDSumQZPAlg = new TH1F("hESDSumQZPAlg", "Sum of 4 lg ZPA sectors",100, 0., 4000.);
  Add2ESDsList(hESDSumQZNClg, 18, expert, !image);
  Add2ESDsList(hESDSumQZNAlg, 19, expert, !image);
  Add2ESDsList(hESDSumQZPClg, 20, expert, !image);
  Add2ESDsList(hESDSumQZPAlg, 21, expert, !image);
  //
  TH1F * hESDPMCZNClg = new TH1F("hESDPMCZNClg", "Signal in common ZNC lg PMT",100, 0., 4000.);
  TH1F * hESDPMCZNAlg = new TH1F("hESDPMCZNAlg", "Signal in common ZNA lg PMT",100, 0., 4000.);
  TH1F * hESDPMCZPClg = new TH1F("hESDPMCZPClg", "Signal in common ZPC lg PMT",100, 0., 4000.);
  TH1F * hESDPMCZPAlg = new TH1F("hESDPMCZPAlg", "Signal in common ZPA lg PMT",100, 0., 4000.);
  Add2ESDsList(hESDPMCZNClg, 22, expert, !image);
  Add2ESDsList(hESDPMCZNAlg, 23, expert, !image);
  Add2ESDsList(hESDPMCZPClg, 24, expert, !image);
  Add2ESDsList(hESDPMCZPAlg, 25, expert, !image);
}

//___________________________________________________________________________
void AliZDCQADataMakerRec::MakeDigits(TTree *digitTree)
{
  // Check id histograms already created for this Event Specie
  if(!GetDigitsData(0)) InitDigits() ;
  
  TBranch * branch = digitTree->GetBranch("ZDC");
  if(!branch){
    AliError("ZDC branch in Digit Tree not found"); 
    return;
  } 
  
  
  branch->SetAddress(&fDigit);
  
  Long64_t ndig = digitTree->GetEntries();

  for(Int_t i=0; i<ndig; i++){
    digitTree->GetEntry(i);
    if(fDigit->GetSector(0)==1 && fDigit->GetSector(1)!=5){
      if(fDigit->GetSector(1)==0){
	      GetDigitsData(0)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(6)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==1){
	      GetDigitsData(12)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(28)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==2){
	      GetDigitsData(13)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(29)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==3){
	      GetDigitsData(14)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(30)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==4){
	      GetDigitsData(15)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(31)->Fill(fDigit->GetADCValue(1));
      }
    }
    else if(fDigit->GetSector(0)==2){
      if(fDigit->GetSector(1)==0){
	      GetDigitsData(2)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(8)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==1){
	      GetDigitsData(16)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(32)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==2){
	      GetDigitsData(17)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(33)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==3){
	      GetDigitsData(18)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(34)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==4){
	      GetDigitsData(19)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(35)->Fill(fDigit->GetADCValue(1));
      }
    }
    else if(fDigit->GetSector(0)==3){
      if(fDigit->GetSector(1)==1){
	      GetDigitsData(4)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(10)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==2){
	      GetDigitsData(5)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(11)->Fill(fDigit->GetADCValue(1));
      }
    }
    else if(fDigit->GetSector(0)==4 && fDigit->GetSector(1)!=5){
      if(fDigit->GetSector(1)==0){
	      GetDigitsData(1)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(7)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==1){
	      GetDigitsData(20)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(36)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==2){
	      GetDigitsData(21)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(37)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==3){
	      GetDigitsData(22)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(38)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==4){
	      GetDigitsData(23)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(39)->Fill(fDigit->GetADCValue(1));
      }
    }
    else if(fDigit->GetSector(0)==5){
      if(fDigit->GetSector(1)==0){
	      GetDigitsData(3)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(9)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==1){
	      GetDigitsData(24)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(40)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==2){
	      GetDigitsData(25)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(41)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==3){
	      GetDigitsData(26)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(42)->Fill(fDigit->GetADCValue(1));
      }
      else if(fDigit->GetSector(1)==4){
	      GetDigitsData(27)->Fill(fDigit->GetADCValue(0));
	      GetDigitsData(43)->Fill(fDigit->GetADCValue(1));
      }
    }
  }
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::MakeRecPoints(TTree * clustersTree)
{
  // Filling QA histos from RecPoints

  TBranch *branch = clustersTree->GetBranch("ZDC");
  if(!branch){ 
    AliError("Can't get the ZDC branch for rec points!");
    return;
  }
  
  if(!GetRecPointsData(0)) InitRecPoints() ;

  Float_t sum_ZNC=0., sum_ZNA=0., sum_ZPC=0., sum_ZPA=0.;
  Float_t sumQ_ZNC=0., sumQ_ZNA=0., sumQ_ZPC=0., sumQ_ZPA=0.;

  AliZDCReco reco;
  AliZDCReco* preco = &reco;
  clustersTree->SetBranchAddress("ZDC", &preco);

  clustersTree->GetEntry(0);
  for(Int_t i=0; i<5; i++){
    sum_ZNC += reco.GetZN1HREnTow(i);
    sum_ZPC += reco.GetZN2HREnTow(i);
    sum_ZNA += reco.GetZP1HREnTow(i);
    sum_ZPA += reco.GetZP2HREnTow(i);
    if(i==0){
      GetRecPointsData(8)->Fill(reco.GetZN1HREnTow(i));
      GetRecPointsData(9)->Fill(reco.GetZN2HREnTow(i));
      GetRecPointsData(10)->Fill(reco.GetZP1HREnTow(i));
      GetRecPointsData(11)->Fill(reco.GetZP2HREnTow(i));
    }
    else{
      sumQ_ZNC += reco.GetZN1HREnTow(i);
      sumQ_ZPC += reco.GetZN2HREnTow(i);
      sumQ_ZNA += reco.GetZP1HREnTow(i);
      sumQ_ZPA += reco.GetZP2HREnTow(i);
    }
  }
  
  GetRecPointsData(0)->Fill(sum_ZNC);
  GetRecPointsData(1)->Fill(sum_ZPC);
  GetRecPointsData(2)->Fill(sum_ZNA);
  GetRecPointsData(3)->Fill(sum_ZPA);
  //
  GetRecPointsData(4)->Fill(sumQ_ZNC);
  GetRecPointsData(5)->Fill(sumQ_ZPC);
  GetRecPointsData(6)->Fill(sumQ_ZNA);
  GetRecPointsData(7)->Fill(sumQ_ZPA);
  
}  

//____________________________________________________________________________
void AliZDCQADataMakerRec::MakeRaws(AliRawReader *rawReader)
{
  // Filling Raws QA histos
  //
  // Check if histograms already created for this Event Specie
  if(!GetRawsData(0)) InitRaws();
  
  AliZDCRawStream stream(rawReader);
  while(stream.Next()){
    if(stream.IsADCDataWord() && 
     (stream.GetADCModule()==0 || stream.GetADCModule()==1)){
     
       if(stream.GetSector(0)==1 && stream.GetSector(1)!=5){ 
         if(stream.GetSector(1)==0){
	   if(stream.GetADCGain()==0) GetRawsData(0)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(6)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==1){
	   if(stream.GetADCGain()==0) GetRawsData(12)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(28)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==2){
	   if(stream.GetADCGain()==0) GetRawsData(13)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(29)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==3){
	   if(stream.GetADCGain()==0) GetRawsData(14)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(30)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==4){
	   if(stream.GetADCGain()==0) GetRawsData(15)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(31)->Fill(stream.GetADCValue());
	 }
       }
       else if(stream.GetSector(0)==2){ 
         if(stream.GetSector(1)==0){
	   if(stream.GetADCGain()==0) GetRawsData(2)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(8)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==1){
	   if(stream.GetADCGain()==0) GetRawsData(16)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(32)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==2){
	   if(stream.GetADCGain()==0) GetRawsData(17)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(33)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==3){
	   if(stream.GetADCGain()==0) GetRawsData(18)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(34)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==4){
	   if(stream.GetADCGain()==0) GetRawsData(19)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(35)->Fill(stream.GetADCValue());
	 }
       }
       else if(stream.GetSector(0)==3){
         if(stream.GetSector(1)==1){ // ZEM1
	   if(stream.GetADCGain()==0) GetRawsData(4)->Fill(stream.GetADCValue());
	   else GetRawsData(10)->Fill(stream.GetADCValue());
	 }
         else if(stream.GetSector(1)==2){ // ZEM2
	   if(stream.GetADCGain()==0) GetRawsData(5)->Fill(stream.GetADCValue());
	   else GetRawsData(11)->Fill(stream.GetADCValue());
	 }
       }
       else if(stream.GetSector(0)==4 && stream.GetSector(1)!=5){
         if(stream.GetSector(1)==0){
	   if(stream.GetADCGain()==0) GetRawsData(1)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(7)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==1){
	   if(stream.GetADCGain()==0) GetRawsData(20)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(36)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==2){
	   if(stream.GetADCGain()==0) GetRawsData(21)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(37)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==3){
	   if(stream.GetADCGain()==0) GetRawsData(22)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(38)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==4){
	   if(stream.GetADCGain()==0) GetRawsData(23)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(39)->Fill(stream.GetADCValue());
	 }
       }
       else if(stream.GetSector(0)==5){
         if(stream.GetSector(1)==0){
	   if(stream.GetADCGain()==0) GetRawsData(3)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(9)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==1){
	   if(stream.GetADCGain()==0) GetRawsData(24)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(40)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==2){
	   if(stream.GetADCGain()==0) GetRawsData(25)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(41)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==3){
	   if(stream.GetADCGain()==0) GetRawsData(26)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(42)->Fill(stream.GetADCValue());
	 }
	 else if(stream.GetSector(1)==4){
	   if(stream.GetADCGain()==0) GetRawsData(27)->Fill(stream.GetADCValue());
	   else if(stream.GetADCGain()==1) GetRawsData(43)->Fill(stream.GetADCValue());
	 }
       }
    }
    
  }

//   stream.Delete();
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs
  //
  
  // Check id histograms already created for this Event Specie
  if(!GetESDsData(0))
    InitESDs() ;

  AliESDZDC * zdcESD =  esd->GetESDZDC();
  //
  TString beamType = esd->GetBeamType();
  const Double_t *centr_ZNC, *centr_ZNA;
  if(((beamType.CompareTo("pp"))==0) || ((beamType.CompareTo("p-p"))==0)
     ||((beamType.CompareTo("PP"))==0) || ((beamType.CompareTo("P-P"))==0)){
    Float_t beamEne = esd->GetBeamEnergy();
    centr_ZNC = zdcESD->GetZNCCentroidInPbPb(beamEne);
    centr_ZNA = zdcESD->GetZNACentroidInPbPb(beamEne);
  }
  else if((beamType.CompareTo("A-A")) == 0){
    centr_ZNC = zdcESD->GetZNCCentroidInpp();
    centr_ZNA = zdcESD->GetZNACentroidInpp();
  }
  GetESDsData(0)->Fill(centr_ZNC[0], centr_ZNC[1]);
  GetESDsData(1)->Fill(centr_ZNA[0], centr_ZNA[1]);

  //
  GetESDsData(2)->Fill(esd->GetZDCN1Energy());
  GetESDsData(3)->Fill(esd->GetZDCN2Energy());
  GetESDsData(4)->Fill(esd->GetZDCP1Energy());
  GetESDsData(5)->Fill(esd->GetZDCP2Energy());
  //
  Double_t sumQZNC=0., sumQZPC=0., sumQZNA=0., sumQZPA=0.;
  Double_t sumQZNC_lg=0., sumQZPC_lg=0., sumQZNA_lg=0., sumQZPA_lg=0.;
  //
  const Double_t *towZNC, *towZPC, *towZNA, *towZPA;
  const Double_t *towZNC_lg, *towZPC_lg, *towZNA_lg, *towZPA_lg;
  //
  towZNC = zdcESD->GetZN1TowerEnergy();
  towZPC = zdcESD->GetZP1TowerEnergy();
  towZNA = zdcESD->GetZN2TowerEnergy();
  towZPA = zdcESD->GetZP2TowerEnergy();
  //
  towZNC_lg = zdcESD->GetZN1TowerEnergyLR();
  towZPC_lg = zdcESD->GetZP1TowerEnergyLR();
  towZNA_lg = zdcESD->GetZN2TowerEnergyLR();
  towZPA_lg = zdcESD->GetZP2TowerEnergyLR();
  //
  for(Int_t i=0; i<5; i++){
     if(i==0){
       GetESDsData(10)->Fill(towZNC[i]);
       GetESDsData(11)->Fill(towZNA[i]);
       GetESDsData(12)->Fill(towZPC[i]);
       GetESDsData(13)->Fill(towZPA[i]);
       //
       GetESDsData(22)->Fill(towZNC_lg[i]);
       GetESDsData(23)->Fill(towZNA_lg[i]);
       GetESDsData(24)->Fill(towZPC_lg[i]);
       GetESDsData(25)->Fill(towZPA_lg[i]);
     }
     else{
       sumQZNC += towZNC[i];
       sumQZPC += towZPC[i];
       sumQZNA += towZNA[i];
       sumQZPA += towZPA[i];
       //
       sumQZNC_lg += towZNC_lg[i];
       sumQZPC_lg += towZPC_lg[i];
       sumQZNA_lg += towZNA_lg[i];
       sumQZPA_lg += towZPA_lg[i];
     }
  }
  GetESDsData(6)->Fill(sumQZNC);
  GetESDsData(7)->Fill(sumQZNA);
  GetESDsData(8)->Fill(sumQZPC);
  GetESDsData(9)->Fill(sumQZPA);
  //
  GetESDsData(18)->Fill(sumQZNC_lg);
  GetESDsData(19)->Fill(sumQZNA_lg);
  GetESDsData(20)->Fill(sumQZPC_lg);
  GetESDsData(21)->Fill(sumQZPA_lg);
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}

//____________________________________________________________________________ 
void AliZDCQADataMakerRec::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kZDC, task, list) ;  
}

