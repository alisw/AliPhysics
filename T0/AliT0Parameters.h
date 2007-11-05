#ifndef ALIT0PARAMETERS_H
#define ALIT0PARAMETERS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */

//____________________________________________________________________
//
//  Singleton class to handle various parameters of the
//  T0 - T0
//  Should get data fromm Conditions DB.
//
# include <TNamed.h>
# include <TF1.h>
# include <TMap.h>
# include <TGraph.h>
#include <TObjArray.h>
class AliT0CalibData;
class AliCDBEntry;

class AliT0Parameters : public TNamed
{
public:
  static AliT0Parameters* Instance();
  
  AliT0Parameters();
  virtual ~AliT0Parameters() {};
 
  void Init();  
  void InitIfOnline();

  // Set various `Fixed' parameters 
  void SetPh2Mip(Int_t r=300)          { fPh2Mip = r; }
  void SetmV2Mip(Int_t r=50)          { fmV2Mip = r; }
  void SetChannelWidth(Int_t s=24.4)   { fChannelWidth = s;}
  void SetmV2channel(Int_t size=320) { fmV2Channel = size; }
  void SetQTmin(Int_t qt=13) {fQTmin = qt;}
  void SetQTmax(Int_t qt=125) {fQTmax = qt;}
  void SetZposition( Float_t valueC=69.7, Float_t valueA=373) {
    fT0zPosition[0]=valueC, fT0zPosition[1]=valueA;}
  void SetPMTeff(Int_t ipmt);  

  void SetTimeDelayTVD(Float_t r=150)   { fTimeDelayTVD = r; };
  Float_t GetTimeDelayTVD()   { return fTimeDelayTVD; }

 
  // Get `Fixed' various parameters
  Int_t GetPh2Mip()          const { return fPh2Mip; }
  Int_t GetmV2Mip()          const { return fmV2Mip; }
  Int_t GetChannelWidth()     const { return fChannelWidth; }
  Int_t GetmV2channel()     const { return fmV2Channel; }
  Int_t GetQTmin() const {return fQTmin;}
  Int_t GetQTmax() const {return fQTmax;}
  Double_t GetZposition(Int_t i) const {return fT0zPosition[i];}
  Double_t GetZPosition(const char* symname) ;
  Double_t GetZPositionShift(const char* symname);


  TGraph *  GetPMTeff(Int_t ipmt) const  
  {return (TGraph*)fPMTeff.At(ipmt);}
  Float_t GetpmtEFF(Int_t ipmt, Float_t lambda) const
  {return((TGraph*)fPMTeff.At(ipmt))->Eval(lambda);} 


  TGraph *  GetAmpLEDRec(Int_t ipmt) const;  
 
  TGraph *GetWalk(Int_t ipmt )  const;
  Float_t  GetWalkVal(Int_t ipmt, Float_t mv ) const ;
   
  Float_t GetTimeDelayCFD(Int_t ipmt);
  Float_t GetTimeDelayDA(Int_t ipmt);

  //  void SetMeanT0(Int_t mean=500) { fMeanT0 = mean; };
  Int_t GetMeanT0 (); //{return fMeanT0;};

    TMap * GetMapLookup();
  Int_t GetChannel(Int_t trm,  Int_t tdc, Int_t chain, Int_t channel);
  Int_t GetNumberOfTRMs();
  void SetNumberOfTRMs(Int_t ntrms=2) {fNumberOfTRMs = ntrms;}

 protected:
  static AliT0Parameters* fgInstance; // Static singleton instance
  
  Bool_t    fIsInit;                // Whether we've been initialised
  Float_t   fT0zPosition[2] ;  // z-position of the two T0s
  Int_t     fPh2Mip;            // # photoelectrons per MIP in radiator
  Int_t     fmV2Mip;            // # mV per MIP in radiator
  Int_t     fChannelWidth;          // channel width in ns   
  Int_t     fmV2Channel;     // ADC mv  2  channel # (200000ps/(25*25).
  Int_t     fQTmin;                 //min  time for QTC
  Int_t     fQTmax;                 //max  time fro QTC 
  TObjArray fAmpLEDRec;  // array of amlitude vs LED-CFD (simulation & reconstruction)
  TObjArray fPMTeff; //array PMT registration efficiency
  TObjArray fWalk; //array time-amplitude walk
  
  Float_t   fTimeDelayDA;  //  sum time delay for LED channel
  Float_t   fTimeDelayCFD;  // sum time delay for CFD channel
  Float_t   fTimeDelayTVD;  //time delay for TVD (vertex trigger channel)
  Int_t     fMeanT0; //mean of T0distribution with vertex=0;
  
  TMap      fLookUp;           //lookup table
  Int_t     fNumberOfTRMs;    // number of TRMs in setup
  

  static AliT0CalibData * fgCalibData; // singleton for Calibration data
  static AliT0CalibData * fgLookUp; // singleton for Calibration data
  static AliT0CalibData * fgSlewCorr; // singleton for Calibration data
  
  AliCDBEntry*   fCalibentry ;  // pointer to T0 calibration object
  AliCDBEntry*   fLookUpentry ;  // pointer to T0 lokkup table
  AliCDBEntry*   fSlewCorr ;  // pointer to slewing correction

 private:
  AliT0Parameters(const  AliT0Parameters&);
  AliT0Parameters& operator=(const AliT0Parameters&);
  
  ClassDef(AliT0Parameters,4)
 
};

#endif
//____________________________________________________________________

