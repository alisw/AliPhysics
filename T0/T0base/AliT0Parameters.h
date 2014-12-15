#ifndef ALIT0PARAMETERS_H
#define ALIT0PARAMETERS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Alla Maevskaya INR RAS alla@inr.ru
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
class AliT0CalibWalk;
class AliT0CalibTimeEq;
class AliT0CalibLatency;

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
  void SetChannelWidth(Float_t s=24.4)   { fChannelWidth = s;}
  void SetmV2channel(Int_t size=320) { fmV2Channel = size; }
  void SetQTmin(Int_t qt=13) {fQTmin = qt;}
  void SetQTmax(Int_t qt=125) {fQTmax = qt;}
  void SetZposition( Float_t valueC=69.7, Float_t valueA=373) {
    fT0zPosition[0]=valueC, fT0zPosition[1]=valueA;}
  void SetPMTeff(Int_t ipmt);  

  void SetTimeDelayTVD(Float_t r=150)   { fTimeDelayTVD = r; };
  Float_t GetTimeDelayTVD() const   { return fTimeDelayTVD; }

 
  // Get `Fixed' various parameters
  Int_t GetPh2Mip()          const { return fPh2Mip; }
  Int_t GetmV2Mip()          const { return fmV2Mip; }
  Float_t GetChannelWidth()     const { return fChannelWidth; }
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


  TGraph *GetAmpLEDRec(Int_t ipmt) const;  
  TGraph *GetWalk(Int_t ipmt )  const;
  TGraph *GetQTC(Int_t ipmt) const;
  TGraph *GetAmpLED(Int_t ipmt) const;
   
  Float_t GetTimeDelayCFD(Int_t ipmt);
//  Float_t GetTimeV0(Int_t ipmt = 512) {return  fTimeV0;}
  Float_t GetCFD (Int_t ipmt);
  void SetMeanT0(Float_t mean=512) { fMeanT0 = mean; };
  Float_t GetMeanT0 () {return fMeanT0;};
  void SetMeanVertex(Float_t mean=0) { fMeanVertex = mean; };
   Float_t GetMeanVertex ();

  TMap * GetMapLookup();
  Int_t GetChannel(Int_t trm,  Int_t tdc, Int_t chain, Int_t channel);
  Int_t GetNumberOfTRMs();
  void SetNumberOfTRMs(Int_t ntrms=2) {fNumberOfTRMs = ntrms;}
 
  Float_t GetLatencyHPTDC();
  Float_t GetLatencyL1();/* {return fLatencyL1;} */
  Float_t GetLatencyL1A(); /* {return fLatencyL1A;}*/ 
  Float_t GetLatencyL1C(); /* {return fLatencyL1C;} */
 
  void SetLatencyHPTDC(Float_t lat) {fLatencyHPTDC=lat;} 
  void SetLatencyL1(Float_t lat) {fLatencyL1=lat;} 
  void SetLatencyL1A(Float_t lat) { fLatencyL1A=lat;} 
  void  SetLatencyL1C(Float_t lat) { fLatencyL1C=lat;} 

 protected:
  static AliT0Parameters* fgInstance; // Static singleton instance
  
  Bool_t    fIsInit;                // Whether we've been initialised
  Float_t   fT0zPosition[2] ;  // z-position of the two T0s
  Int_t     fPh2Mip;            // # photoelectrons per MIP in radiator
  Int_t     fmV2Mip;            // # mV per MIP in radiator
  Float_t     fChannelWidth;          // channel width in ns   
  Int_t     fmV2Channel;     // ADC mv  2  channel # (200000ps/(25*25).
  Int_t     fQTmin;                 //min  time for QTC
  Int_t     fQTmax;                 //max  time fro QTC 
  TObjArray fAmpLEDRec;  // array of amlitude vs LED-CFD (simulation & reconstruction)
  TObjArray fPMTeff; //array PMT registration efficiency
  TObjArray fWalk; //array time-amplitude walk
  TObjArray fQTC; //array of TGraphs for QTC vs number of MIPs
  TObjArray fAmpLED; //array of TGraphs for LED-CFD vs number of MIPs

  
  Float_t   fTimeDelayCFD;  // sum time delay for CFD channel
 // Float_t   fTimeV0;  // sum time delay for CFD channel
  Float_t   fTimeDelayTVD;  //time delay for TVD (vertex trigger channel)
  Float_t     fMeanT0; //mean of T0distribution with vertex=0;
  Float_t     fMeanVertex; // mean of vertex distribution;
   
  Float_t  fLatencyHPTDC; // all latencies;
  Float_t  fLatencyL1; // all latencies;
  Float_t  fLatencyL1A; // all latencies;
  Float_t  fLatencyL1C; // all latencies;
    
  TMap      fLookUp;           //lookup table
  Int_t     fNumberOfTRMs;    // number of TRMs in setup
  
  //latency


  static AliT0CalibTimeEq * fgCalibData; // singleton for Calibration data
  static AliT0CalibData * fgLookUp; // singleton for Calibration data
  static AliT0CalibWalk * fgSlewCorr; // singleton for Calibration data
  static AliT0CalibLatency * fgLatency; // singleton for Calibration data
  
  AliCDBEntry*   fCalibentry ;  // pointer to T0 calibration object
  AliCDBEntry*   fLookUpentry ;  // pointer to T0 lokkup table
  AliCDBEntry*   fSlewCorr ;  // pointer to slewing correction
  AliCDBEntry*   fLatency ;  // pointer to latency

 private:
  AliT0Parameters(const  AliT0Parameters&);
  AliT0Parameters& operator=(const AliT0Parameters&);
  
  ClassDef(AliT0Parameters,6)
 
};

#endif
//____________________________________________________________________

