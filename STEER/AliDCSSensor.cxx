/**************************************************************************
 * Copyright(c) 2006-07, ALICE Experiment at CERN, All rights reserved. *
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


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Class describing TPC temperature sensors (including pointers to graphs/fits//
// Authors: Marian Ivanov, Haavard Helstrup and Martin Siska                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


// Running instructions:
/*
  TClonesArray * arr = AliDCSSensor::ReadList("TempSensor.txt");
  TFile f("TempSensors.root","RECREATE");
  TTree * tree = new TTree("TempSensor", "TempSensor");
  tree->Branch("Temp",&arr);
  tree->Fill();
  tree->Write();
  
 */
//


#include "AliDCSSensor.h"
ClassImp(AliDCSSensor)


AliDCSSensor::AliDCSSensor():
  fId(),
  fIdDCS(0),
  fStartTime(0),
  fGraph(0),
  fFit(0),
  fX(0),
  fY(0),
  fZ(0)
{
  //
  //  Standard constructor
  //
}

AliDCSSensor::AliDCSSensor(const AliDCSSensor& source) :
   TNamed(source),
   fId(source.fId),
   fIdDCS(source.fIdDCS),
   fStartTime(source.fStartTime),
   fGraph(source.fGraph),
   fFit(source.fFit),
   fX(source.fX),
   fY(source.fY),
   fZ(source.fZ)
//
//  Copy constructor
//
{ }

AliDCSSensor& AliDCSSensor::operator=(const AliDCSSensor& source){
//
// assignment operator
//
  if (&source == this) return *this;
  new (this) AliDCSSensor(source);
  
  return *this;  
}

//_____________________________________________________________________________
Double_t AliDCSSensor::GetValue(UInt_t timeSec) 
{
 //
 // Get temperature value for actual sensor
 //  timeSec given as offset from start-of-run measured in seconds
 //
 Double_t timeHrs = timeSec/3600.0;
 return  fFit->Eval(timeHrs,0);
}
//_____________________________________________________________________________
Double_t AliDCSSensor::GetValue(TTimeStamp time) 
{
 // Get temperature value for actual sensor
 //  time given as absolute TTimeStamp
 //
 Double_t timeHrs = (time.GetSec() - fStartTime)/3600.0;
 return fFit->Eval(timeHrs,0);
}


