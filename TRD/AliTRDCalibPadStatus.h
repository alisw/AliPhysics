#ifndef ALITRDCALIBPADSTATUS_H
#define ALITRDCALIBPADSTATUS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for online calibration                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class TObjArray;
class TH2F;

class AliRawReader;

class AliTRDCalDet;
class AliTRDCalPad;
class AliTRDCalROC;
class AliTRDCalPadStatus;
class AliTRDrawStreamBase;
class AliTRDgeometry;

class AliTRDrawFastStream;
class AliTRDdigitsManager;
class AliTRDSignalIndex;

struct eventHeaderStruct;

class AliTRDCalibPadStatus : public TObject {

public:

  AliTRDCalibPadStatus();
  AliTRDCalibPadStatus(const AliTRDCalibPadStatus &ped);
  virtual ~AliTRDCalibPadStatus();

  AliTRDCalibPadStatus& operator = (const  AliTRDCalibPadStatus &source);

  Int_t ProcessEvent(AliTRDrawStreamBase *rawStream, Bool_t nocheck = kFALSE);
  Int_t ProcessEvent(AliRawReader    *rawReader, Bool_t nocheck = kFALSE);
  Int_t ProcessEvent(const eventHeaderStruct   *event, Bool_t nocheck = kFALSE);
  Int_t ProcessEvent2(AliRawReader    *rawReader);
  Int_t ProcessEvent3(AliRawReader    *rawReader);
 
  void  Destroy();
  Int_t UpdateHisto(const Int_t idet, const Int_t iRow, const Int_t iCol,
		    const Int_t signal, const Int_t crowMax, const Int_t ccold, const Int_t icMcm);

  Int_t UpdateHisto2(const Int_t idet, const Int_t iRow, const Int_t iCol,
		     const Int_t signal, const Int_t crowMax, const Int_t ccold, const Int_t icMcm, const Int_t icRob);

  void AnalyseHisto();
  AliTRDCalPadStatus *CreateCalPadStatus();
  AliTRDCalPad *CreateCalPad();
  AliTRDCalDet *CreateCalDet() const;

  void SetCalRocMean(AliTRDCalROC *mean, Int_t det);
  void SetCalRocRMS(AliTRDCalROC *rms, Int_t det);  

  void SetCalRocMeand(AliTRDCalROC *mean, Int_t det);
  void SetCalRocRMSd(AliTRDCalROC *rms, Int_t det);  

  //
  AliTRDCalROC* GetCalRocMean(Int_t det, Bool_t force=kFALSE);    // get calibration object
  AliTRDCalROC* GetCalRocRMS(Int_t det, Bool_t force=kFALSE);     // get calibration object

  AliTRDCalROC* GetCalRocMeand(Int_t det, Bool_t force=kFALSE);   // get calibration object
  AliTRDCalROC* GetCalRocRMSd(Int_t det, Bool_t force=kFALSE);    // get calibration object

  TH2F* GetHisto  (Int_t det, Bool_t force=kFALSE);              // get refernce histogram
  
  void  DumpToFile(const Char_t *filename, const Char_t *dir="", Bool_t append=kFALSE);
  //
  Int_t   GetAdcMin()       const { return fAdcMin;       }
  Int_t   GetAdcMax()       const { return fAdcMax;       }

  void    SetRangeAdc (Int_t aMin, Int_t aMax){ fAdcMin=aMin; fAdcMax=aMax; }  // Set adc range 


  Bool_t TestEventHisto(Int_t nevent, Int_t sm, Int_t ch);  //test the fast approach to fill histograms  

 private:

  // Geometry
  AliTRDgeometry  *fGeo;            //! The TRD geometry

  Int_t fAdcMin;                    //  min adc channel of pedestal value
  Int_t fAdcMax;                    //  max adc channel of pedestal value
  Int_t fDetector;                  //  Current detector
  Int_t fNumberOfTimeBins;          //  Current number of time bins
     
  TObjArray fCalRocArrayMean;       //  Array of AliTRDCalROC class for signal width calibration
  TObjArray fCalRocArrayRMS;        //  Array of AliTRDCalROC class for mean width calibration

  TObjArray fCalRocArrayMeand;      //  Array of AliTRDCalROC class for signal width calibration doubled
  TObjArray fCalRocArrayRMSd;       //  Array of AliTRDCalROC class for mean width calibration doubled

  TObjArray fHistoArray;            //  Array of histos for mean width calibration
  
 
  AliTRDCalROC* GetCalRoc(Int_t det, TObjArray* arr, Bool_t force);
 
  TH2F* GetHisto(Int_t det, TObjArray *arr,
		 Int_t nbinsY, Float_t ymin, Float_t ymax,
		 const Char_t *type, Bool_t force);

  // Some basic geometry function
  virtual Int_t    GetLayer(Int_t d) const;
  virtual Int_t    GetStack(Int_t d) const;
  virtual Int_t    GetSector(Int_t d) const;

  ClassDef(AliTRDCalibPadStatus,1)

};
#endif
