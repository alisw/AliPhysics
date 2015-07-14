#ifndef ALITPCSAMPAEMULATOR_H
#define ALITPCSAMPAEMULATOR_H


/**    @file AliTPCSAMPAEmulator.h
       @brief This the header File for the SAMPA class
 
       author: marian.ivanov@cern.ch
               mesut.arslandok@cern.ch  
*/


///////////////////////////////////////////////////////////////////////////////
//                        Class AliTPCSAMPAEmulator                          //
//  Class for emulation of the ALTRO chip (Altro digital Chain) in C++       //
///////////////////////////////////////////////////////////////////////////////

#include "TSystem.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <AliRawReader.h>
#include <AliTPCRawStreamV3.h>


using namespace std;

class AliTPCSAMPAEmulator : public TNamed {
public:
  AliTPCSAMPAEmulator();
  ~AliTPCSAMPAEmulator();
  Bool_t  DigitalFilterFloat(Int_t npoints, Double_t *dataArray,  Double_t &baseline);
  Bool_t  ZeroSuppression(Int_t npoints, Double_t *dataArray, Double_t threshold);
  //
  Bool_t  BC3SlopeFilterFloat(Int_t npoints, Double_t *dataArray,  Double_t &baseline);
  static Bool_t  BC3SlopeFilterFloat(Int_t npoints, Double_t *dataArray, Double_t slopeDown, Double_t slopeUp, Double_t round, Double_t &baseline);
  void SetBC3Parameters(Double_t slopeDown, Double_t slopeUp, Double_t round);
  void SetBC3DiffCutMI(Double_t BC3DiffCutMI){fBC3DiffCutMI=BC3DiffCutMI;}
  Bool_t  BC3SlopeFilterMI(Int_t npoints, Double_t *dataArray,  Double_t &baseline);
  static Bool_t  BC3SlopeFilterMI(Int_t npoints, Double_t *dataArray, Double_t slopeDown, Double_t slopeUp, Double_t round, Double_t &baseline, Double_t diffCutMI);
  //
  Bool_t  MovingAverageFilter(Int_t npoints, Double_t *dataArray, Double_t &baseline);
  static Bool_t  MovingAverageFilter(Int_t npoints, Double_t *dataArray, Double_t length, Double_t skipDiff,  Bool_t onlyMinima, Double_t &baseline); // local maxim need treatment
  void SetMAFMIParameters(Double_t  MAFMIKernelWidth,  Double_t  MAFMIDiffCut, Bool_t onlyMinima);
private:
  AliTPCSAMPAEmulator(const AliTPCSAMPAEmulator &sig);
  AliTPCSAMPAEmulator& operator = (const  AliTPCSAMPAEmulator &source);
public:
  Int_t fDigitFilterType;     //
  static Int_t fgBaselineExportType;  // swith to export either corrected signal or the perdestal estimator itself  - needed for the performance ssudies   
  //
  // BC3 parameters
  //
  Double_t  fBC3SlopeDown;  // BC3 slope down parameter
  Double_t  fBC3SlopeUp;    // BC3 slope up   parameter
  Double_t  fBC3Round;      // Rounding error of BC3 filter
  Double_t  fBC3DiffCutMI;  // BC3 cut on the signal difference - for MI implementation

  //
  // Moving average filter parameters (MI) implementation
  //
  Double_t  fMAFMIKernelWidth;   // kernel width for MAF filtering 
  Double_t  fMAFMIDiffCut;       // cut on the diff to skip "signal"
  Double_t  fMAFMIOnlyMinima;         // use only local minima



  ClassDef(AliTPCSAMPAEmulator,1);
};
#endif
