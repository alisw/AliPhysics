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
class AliTPCCalibRaw;

class AliTPCcalibDButil : public TObject
{
public:
  AliTPCcalibDButil();
  virtual ~AliTPCcalibDButil();

  void UpdateFromCalibDB();
  //data processing functions
  void ProcessCEdata(const char* fitFormula, TVectorD &fitResultsA, TVectorD &fitResultsC,
                     Int_t &noutliersCE, Double_t & chi2A, Double_t &chi2C, AliTPCCalPad *outCE=0);
  void ProcessCEgraphs(TVectorD &vecTEntries, TVectorD &vecTMean, TVectorD &vecTRMS, TVectorD &vecTMedian,
                       TVectorD &vecQEntries, TVectorD &vecQMean, TVectorD &vecQRMS, TVectorD &vecQMedian,
                       Float_t &driftTimeA, Float_t &driftTimeC );
  void ProcessNoiseData(TVectorD &vNoiseMean, TVectorD &vNoiseMeanSenRegions,
                        TVectorD &vNoiseRMS, TVectorD &vNoiseRMSSenRegions,
                        Int_t &nonMaskedZero);
  void ProcessPulser(TVectorD &vMeanTime);
  void ProcessALTROConfig(Int_t &nMasked);
  void ProcessGoofie(TVectorD & vecEntries, TVectorD & vecMedian, TVectorD &vecMean, TVectorD &vecRMS);
  //processing functions using reference data
  void ProcessPedestalVariations(TVectorF &pedestalDeviations);
  void ProcessNoiseVariations(TVectorF &noiseDeviations);
  void ProcessPulserVariations(TVectorF &pulserQdeviations, Float_t &varQMean, Int_t &npadsOutOneTB, Int_t &npadsOffAdd);
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
  //setters for pad by pad information
  void SetRefFile(const char* filename);
  void SetRefPulserData(AliTPCCalPad *tmean, AliTPCCalPad *trms=0x0, AliTPCCalPad *qmean=0x0)
                {fRefPulserTmean=tmean; fRefPulserTrms=trms; fRefPulserQmean=qmean;}
  void SetRefCEData(AliTPCCalPad *tmean, AliTPCCalPad *trms=0x0, AliTPCCalPad *qmean=0x0)
                {fRefCETmean=tmean; fRefCETrms=trms; fRefCEQmean=qmean;}
  void SetRefNoisePedestal(AliTPCCalPad *noise, AliTPCCalPad *pedestal=0x0)
                {fRefPadNoise=noise; fRefPedestals=pedestal;}
  void SetRefALTROData(AliTPCCalPad *masked)
                {fRefALTROMasked=masked;}
  
  //creation of derived pad by pad calibration data
  AliTPCCalPad *CreatePadTime0(Int_t model, Double_t &gyA, Double_t &gyC, Double_t &chi2A, Double_t &chi2C);
  //
  // create outlyer maps
  //
  AliTPCCalPad *CreateCEOutlyerMap(Int_t &noutliersCE, AliTPCCalPad *ceOut=0, Float_t minSignal=10, Float_t cutTrmsMin=0.9, Float_t cutTrmsMax=1.2, Float_t cutMaxDistT=0.7);
  AliTPCCalPad *CreatePulserOutlyerMap(Int_t &noutliersPulser, AliTPCCalPad *pulserOut=0, Float_t cutTime=3, Float_t cutnRMSQ=5, Float_t cutnRMSrms=5);
  //
  AliTPCCalPad *CreatePadTime0CE(TVectorD &fitResultsA, TVectorD&fitResultsC, Int_t &nOut, Double_t &chi2A, Double_t &chi2C, const char *dumpfile=0);
  //

  void UpdatePulserOutlierMap();
  void UpdateRefPulserOutlierMap();
  void PulserOutlierMap(AliTPCCalPad *pulOut, const AliTPCCalPad *pulT, const AliTPCCalPad *pulQ);
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
  AliTPCCalibRaw *fCalibRaw;          //raw calibration object
  //reference data
  AliTPCCalPad  *fRefPadNoise;           //Reference noise information
  AliTPCCalPad  *fRefPedestals;          //Reference pedestal information
  AliTPCCalPad  *fRefPulserTmean;        //Reference pulser mean time information
  AliTPCCalPad  *fRefPulserTrms;         //Reference pulser rms time information
  AliTPCCalPad  *fRefPulserQmean;        //Reference pulser mean q information
  AliTPCCalPad  *fRefPulserOutlier;      //Reference pulser outlier map
  AliTPCCalPad  *fRefCETmean;            //Reference central electrode mean time information
  AliTPCCalPad  *fRefCETrms;             //Reference central electrode rms time information
  AliTPCCalPad  *fRefCEQmean;            //Reference central electrode mean q information
  AliTPCCalPad  *fRefALTROMasked;        //Reference ALTRO masked channels information
  //
  AliTPCCalibRaw *fRefCalibRaw;          //Reference raw calibration object
  
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
