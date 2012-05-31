#ifndef ALIEMCALSMCALIBCOSMICRESULT_H
#define ALIEMCALSMCALIBCOSMICRESULT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for EMCAL cosmic calibration results                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


/*
  Objects of this class read files with info on
  1) gain balanced voltage values per APD/tower
  2) MIP peak values per tower
  3) LED peak values per tower
  4) LED reference values per strip

  Other values (arrays of smaller dimensions) are assigned via set/get calls:
  a) SM id a la US1, EU1 etc., or W1,Y1 or G1 or whatever you prefer, set via ctor 
  b) run numbers used
  c) start of run, end of run 
  d) temperatures (min, max, mean)

*/



#include <TObject.h>

class TString;


class AliEMCALSMCalibCosmicResult : public TObject {
public:

  AliEMCALSMCalibCosmicResult(); // default ctor
  AliEMCALSMCalibCosmicResult(const TString &smId); // ctor
  
  virtual ~AliEMCALSMCalibCosmicResult(); // dtor

  // simple setters and getters: b)-d) above
  // first setters
  void SetRunNumber(int ipart, int i)  { fRunNumber[ipart]=i; }; // runnumber used for the datataking of a certain part 'ipart' of the SM
  void SetStartTime(int ipart, int i)  { fStartTime[ipart]=i; }; // start of run timestamp for the datataking of 'ipart'
  void SetEndTime(int ipart, int i)    { fEndTime[ipart]=i;   }; // stop/end of run timestamp for the datataking of 'ipart'
  void SetMinTemp(int ipart, double d) { fMinTemp[ipart]=d;   }; // minimum temperature for the datataking of 'ipart'
  void SetMaxTemp(int ipart, double d) { fMaxTemp[ipart]=d;   }; // maximum temperature for the datataking of 'ipart'
  void SetMeanTemp(int ipart, double d){ fMeanTemp[ipart]=d;  }; // mean (not necessary (max+min)/2) temperature for the datataking of 'ipart'
 
  void InitArrays();                     //init all arrays.
  void ReadSMPart(int ipart);
 
  // then getters
  int GetRunNumber(int ipart)   const { return fRunNumber[ipart]; }; // runnumber used for the datataking of a certain part 'ipart' of the SM
  int GetStartTime(int ipart)   const { return fStartTime[ipart]; }; // start of run timestamp for the datataking of 'ipart'
  int GetEndTime(int ipart)     const { return fEndTime[ipart];   }; // stop/end of run timestamp for the datataking of 'ipart'
  double GetMinTemp(int ipart)  const { return fMinTemp[ipart];   }; // minimum temperature for the datataking of 'ipart'
  double GetMaxTemp(int ipart)  const { return fMaxTemp[ipart];   }; // maximum temperature for the datataking of 'ipart'
  double GetMeanTemp(int ipart) const { return fMeanTemp[ipart];   }; // mean (not necessary (max+min)/2) temperature for the datataking of 'ipart'


  // readers: To interface with Yale's output files 
  int  ReadAPDVoltageValues(int ipart, const TString &txtFile); 
  void ReadLEDRefADCValues(int ipart); 
  int  ReadMIPPeakADCValues(int th); 
  int  ReadLEDPeakADCValues(int th); 
  void ReadTempSensors(int ipart, int run); 
  
  // and some methods to access the info stored in the larger arrays
  double GetAPDVoltage(int icol, int irow) const { return fAPDVoltage[icol][irow]; }; 
  double GetMIPPeakADC(int icol, int irow) const { return fMIPPeakADC[icol][irow]; }; 
  double GetLEDPeakADC(int icol, int irow) const { return fLEDPeakADC[icol][irow]; }; 
  double GetLEDRefADC(int icol)            const { return fLEDRefADC[icol]; }; 

  // print functions
  void PrintAPD();
  void PrintMIP();
  void PrintLED();
  void PrintLEDref();
  void PrintTemps();
 
private:
  AliEMCALSMCalibCosmicResult(const AliEMCALSMCalibCosmicResult &);
  AliEMCALSMCalibCosmicResult &operator = (const AliEMCALSMCalibCosmicResult &);

protected:

  TString fSMId; // assigned in ctor
    
  TTree *fTree;  //tree containing results from MIP and LED peak fitting  
  TH2D  *fhref;	 // histograms of LED reference data
  
  float fMip;   // mean ADC of MIP for a given tower
  float fLed;   // mean ADC of LED for a given tower	
  int fCol;     // column index of tower
  int fRow;	// row index of tower
    
  static const int fgkEmCalRows = 24; // number of rows per module for EMCAL
  static const int fgkEmCalCols = 48; // number of columns per module for EMCAL
  static const int fgkEmCalStrips = 24; // number of strips per module for EMCAL (fgkEmCalCols/2)
  static const int fgkCalibParts = 3; // we do the calibrations in 3 different parts; 8 strips or 16 columns at a time

  // small dim array values; assigned via Set methods above
  int fRunNumber[fgkCalibParts];  //runnumber of calibration part
  int fStartTime[fgkCalibParts];  //start time of calibration part
  int fEndTime[fgkCalibParts];    //end time of calibration part

  double fMinTemp[fgkCalibParts];  // min temperature of calibration part
  double fMaxTemp[fgkCalibParts];  // max temperature of calibration part
  double fMeanTemp[fgkCalibParts]; // mean temperature of calibration part

  // larger arrays; assigned via ReadXValues methods above
  double fAPDVoltage[fgkEmCalCols][fgkEmCalRows];  // Adjusted APD bias voltage
  double fMIPPeakADC[fgkEmCalCols][fgkEmCalRows];  // mean of MIP peak
  double fLEDPeakADC[fgkEmCalCols][fgkEmCalRows];  // mean of LED peak
  double fLEDRefADC[fgkEmCalCols];                 // LED reference 

  ClassDef(AliEMCALSMCalibCosmicResult, 1) //SMCalibCosmicResult data reader
};

#endif
