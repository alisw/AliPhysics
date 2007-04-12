#ifndef ALITPCCALIBPEDESTAL_H
#define ALITPCCALIBPEDESTAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */





class TObjArray;
class TH2F;
class TH2S;
class TH1S;
class TH1F;
class TH1D;
class TF1;
class TTreeSRedirector;
class AliTPCROC;
class AliRawReader;


class AliTPCCalibPedestal : public TObject {

public:
  AliTPCCalibPedestal();
  virtual ~AliTPCCalibPedestal();
  
  Bool_t ProcessEvent(AliRawReader *rawReader); 
  Int_t Update(const Int_t isector, const Int_t iRow, const Int_t iPad,
	       const Int_t iTimeBin, const Float_t signal);
  void Analyse();
  //
  AliTPCCalROC* GetCalRocPedestal (Int_t sector, Bool_t force=kFALSE);  //get calibration object - sector
  AliTPCCalROC* GetCalRocRMS(Int_t sector, Bool_t force=kFALSE);        //get calibration object - sector
  const TObjArray* GetCalPadPedestal (){return &fCalRocArrayPedestal;}//get calibration object
  const TObjArray* GetCalPadRMS(){return &fCalRocArrayRMS;}           //get calibration object
  
  TH2S* GetHistoPedestal  (Int_t sector, Bool_t force=kFALSE);          //get refernce histogram
  void DumpToFile(const Char_t *filename, const Char_t *dir="", const Bool_t append=kFALSE);
  //
  Short_t GetDebugLevel(){ return fDebugLevel; }
  void    SetDebugLevel(Short_t debug=1){ fDebugLevel = debug;}


  Bool_t TestEvent();  //test the fast approach to fill histogram  - used for test purposes

private:
  Int_t fFirstTimeBin;              //  First Time bin needed for analysis
  Int_t fLastTimeBin;               //  Last Time bin needed for analysis
  
  Int_t fAdcMin;                    //  min adc channel of pedestal value
  Int_t fAdcMax;                    //  max adc channel of pedestal value
  
  AliTPCROC *fROC;                  //! ROC information
  
  TObjArray fCalRocArrayPedestal;   //  Array of AliTPCCalROC class for Time0 calibration
  TObjArray fCalRocArrayRMS;        //  Array of AliTPCCalROC class for signal width calibration
  
  TObjArray fHistoPedestalArray;    //  Calibration histograms for Pedestal distribution
  
  TTreeSRedirector *fDebugStreamer;  //! debug streamer
  
  Short_t fDebugLevel;
  //! debugging
  
  TH2S* GetHisto(Int_t sector, TObjArray *arr,
		 Int_t nbinsY, Float_t ymin, Float_t ymax,
		 Char_t *type, Bool_t force);
    
  AliTPCCalROC* GetCalRoc(Int_t sector, TObjArray* arr, Bool_t force);

public:


  ClassDef(AliTPCCalibPedestal,1)
};



#endif

