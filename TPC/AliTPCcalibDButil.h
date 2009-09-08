#ifndef AliTPCcalibDButil_H
#define AliTPCcalibDButil_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class providing the calculation of derived quantities (mean,rms,fits,...) //
//       of calibration entries                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliDCSSensorArray;
class AliTPCcalibDB;
class AliTPCCalPad;
class AliTPCmapper;

class AliTPCcalibDButil : public TObject
{
public:
  AliTPCcalibDButil();
  virtual ~AliTPCcalibDButil();

  void UpdateFromCalibDB();
  //data processing functions
  void ProcessCEdata(const char* fitFormula, TVectorD &fitResultsA, TVectorD &fitResultsC, Int_t &noutliersCE);
  void ProcessCEgraphs(TVectorD &vecTEntries, TVectorD &vecTMean, TVectorD &vecTRMS, TVectorD &vecTMedian,
                       TVectorD &vecQEntries, TVectorD &vecQMean, TVectorD &vecQRMS, TVectorD &vecQMedian,
                       Float_t &driftTimeA, Float_t &driftTimeC );
  void ProcessNoiseData(TVectorD &vNoiseMean, TVectorD &vNoiseMeanSenRegions,
                        TVectorD &vNoiseRMS, TVectorD &vNoiseRMSSenRegions,
                        Int_t &nonMaskedZero);
  void ProcessPulser(TVectorD &vMeanTime);
  void ProcessALTROConfig(Int_t &nMasked);
  void ProcessGoofie(TVectorD & vecEntries, TVectorD & vecMedian, TVectorD &vecMean, TVectorD &vecRMS);
  //getter preprocess information
  Int_t GetNPulserOutliers() const {return fNpulserOutliers;}
  Float_t GetMeanAltro(const AliTPCCalROC *roc, const Int_t row, const Int_t pad, AliTPCCalROC *rocOut=0x0);
  AliTPCCalPad *GetPulserOutlierMap() const {return fPulserOutlier;}
  //setters for pad by pad information
  void SetPulserData(AliTPCCalPad *tmean, AliTPCCalPad *trms=0x0, AliTPCCalPad *qmean=0x0)
                {fPulserTmean=tmean; fPulserTrms=trms; fPulserQmean=qmean;}
  void SetCEData(AliTPCCalPad *tmean, AliTPCCalPad *trms=0x0, AliTPCCalPad *qmean=0x0)
                {fCETmean=tmean; fCETrms=trms; fCEQmean=qmean;}
  void SetNoisePedestal(AliTPCCalPad *noise, AliTPCCalPad *pedestal=0x0)
                {fPadNoise=noise; fPedestals=pedestal;}
  void SetALTROData(AliTPCCalPad *masked)
                {fALTROMasked=masked;}
  void SetGoofieArray(AliDCSSensorArray *arr) {fGoofieArray=arr;}
  //creation of derived pad by pad calibration data
  AliTPCCalPad *CreatePadTime0(Int_t model=0);
  //
  void UpdatePulserOutlierMap();
private:
  AliTPCcalibDB *fCalibDB;            //pointer to calibDB object
  AliTPCCalPad  *fPadNoise;           //noise information
  AliTPCCalPad  *fPedestals;          //pedestal information
  AliTPCCalPad  *fPulserTmean;        //pulser mean time information 
  AliTPCCalPad  *fPulserTrms;         //pulser rms time information
  AliTPCCalPad  *fPulserQmean;        //pulser mean q information
  AliTPCCalPad  *fPulserOutlier;      //pulser outlier map
  AliTPCCalPad  *fCETmean;            //central electrode mean time information
  AliTPCCalPad  *fCETrms;             //central electrode rms time information
  AliTPCCalPad  *fCEQmean;            //central electrode mean q information
  AliTPCCalPad  *fALTROMasked;        //ALTRO masked channels information
  //
  AliDCSSensorArray* fGoofieArray;    //Goofie Data
  //
  AliTPCmapper  *fMapper;             //TPC mapping handler
  Int_t fNpulserOutliers;             //number of outliers from Pulser calibration
  
  Float_t fIrocTimeOffset;               //timing offset between IROC and OROC in timebins
  Float_t fCETmaxLimitAbs;               //maximum variation in CE data before pads will be treated as outliers
  Float_t fPulTmaxLimitAbs;              //maximum variation of Pulser Signals (time) before pads will be treated as outliers
  Float_t fPulQmaxLimitAbs;              //maximum variation of Pulser Signals (charge) before pads will be treated as outliers
  Float_t fPulQminLimit;                 //minimum charge value for Pulser Signals before pads will be treated as outliers
  
  AliTPCcalibDButil (const AliTPCcalibDButil& );
  AliTPCcalibDButil& operator= (const AliTPCcalibDButil& );

    
  ClassDef(AliTPCcalibDButil,0)
};


#endif
