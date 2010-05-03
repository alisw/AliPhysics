#ifndef ALITPCCALIBDBUTIL_H
#define ALITPCCALIBDBUTIL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class providing the calculation of derived quantities (mean,rms,fits,...) //
//       of calibration entries                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TArrayI.h>
#include <TVectorD.h>
#include <TVectorF.h>

class TGraph;
class TMap;
class AliDCSSensorArray;
class AliTPCcalibDB;
class AliTPCCalPad;
class AliTPCCalROC;
class AliTPCmapper;
class AliTPCCalibRaw;
class AliCDBEntry;
class AliDCSSensor;
class AliDCSSensorArray;
class AliTPCSensorTempArray;
class AliTPCdataQA;
class TGraphErrors;
class TTreeSRedirector;
class AliTPCCalROC;


class AliTPCcalibDButil : public TObject
{
public:
  AliTPCcalibDButil();
  virtual ~AliTPCcalibDButil();

  void UpdateFromCalibDB();
  //data processing functions
  void ProcessCEdata(const char* fitFormula, TVectorD &fitResultsA, TVectorD &fitResultsC,
                     Int_t &noutliersCE, Double_t & chi2A, Double_t &chi2C, AliTPCCalPad * const outCE=0);
  void ProcessCEgraphs(TVectorD &vecTEntries, TVectorD &vecTMean, TVectorD &vecTRMS, TVectorD &vecTMedian,
                       TVectorD &vecQEntries, TVectorD &vecQMean, TVectorD &vecQRMS, TVectorD &vecQMedian,
                       Float_t &driftTimeA, Float_t &driftTimeC );
  void ProcessNoiseData(TVectorD &vNoiseMean, TVectorD &vNoiseMeanSenRegions,
                        TVectorD &vNoiseRMS, TVectorD &vNoiseRMSSenRegions,
                        Int_t &nonMaskedZero, Int_t &nNaN);
  void ProcessPulser(TVectorD &vMeanTime);
  void ProcessALTROConfig(Int_t &nMasked);
  void ProcessGoofie(TVectorD & vecEntries, TVectorD & vecMedian, TVectorD &vecMean, TVectorD &vecRMS);
  
  //processing functions using reference data
  void ProcessPedestalVariations(TVectorF &pedestalDeviations);
  void ProcessNoiseVariations(TVectorF &noiseDeviations);
  void ProcessPulserVariations(TVectorF &pulserQdeviations, Float_t &varQMean, Int_t &npadsOutOneTB, Int_t &npadsOffAdd);
  
  //getter preprocess information
  Int_t GetNPulserOutliers() const {return fNpulserOutliers;}
  Float_t GetMeanAltro(const AliTPCCalROC *roc, const Int_t row, const Int_t pad, AliTPCCalROC * const rocOut=0x0);
  AliTPCCalPad *GetPulserOutlierMap() const {return fPulserOutlier;}

  //getters ref data
  TMap *GetReferenceMap() const {return fRefMap;}
  Int_t GetReferenceRun(const char* type) const;
  const char* GetRefValidity() const {return fRefValidity.Data();}
  
  AliTPCCalPad* GetRefPadNoise() const {return fRefPadNoise;}
  AliTPCCalPad* GetRefPedestals() const {return fRefPedestals;}
  AliTPCCalPad* GetRefPedestalMasked() const {return fRefPedestalMasked;}
  AliTPCCalPad* GetRefPulserTmean() const {return fRefPulserTmean;}
  AliTPCCalPad* GetRefPulserTrms() const {return fRefPulserTrms;}
  AliTPCCalPad* GetRefPulserQmean() const {return fRefPulserQmean;}
  AliTPCCalPad* GetRefPulserOutlier() const {return fRefPulserOutlier;}
  AliTPCCalPad* GetRefPulserMasked() const {return fRefPulserMasked;}
  AliTPCCalPad* GetRefCETmean() const {return fRefCETmean;}
  AliTPCCalPad* GetRefCETrms() const {return fRefCETrms;}
  AliTPCCalPad* GetRefCEQmean() const {return fRefCEQmean;}
  AliTPCCalPad* GetRefCEMasked() const {return fRefCEMasked;}
  AliTPCCalPad* GetRefALTROFPED() const {return fRefALTROFPED;}
  AliTPCCalPad* GetRefALTROZsThr() const {return fRefALTROZsThr;}
  AliTPCCalPad* GetRefALTROAcqStart() const {return fRefALTROAcqStart;}
  AliTPCCalPad* GetRefALTROAcqStop() const {return fRefALTROAcqStop;}
  AliTPCCalPad* GetRefALTROMasked() const {return fRefALTROMasked;}
  
  
  //setters for pad by pad information
  void SetPulserData(AliTPCCalPad * const tmean, AliTPCCalPad * const trms=0x0, AliTPCCalPad * const qmean=0x0)
                {fPulserTmean=tmean; fPulserTrms=trms; fPulserQmean=qmean;}
  void SetCEData(AliTPCCalPad *const tmean, AliTPCCalPad *const trms=0x0, AliTPCCalPad *const qmean=0x0)
                {fCETmean=tmean; fCETrms=trms; fCEQmean=qmean;}
  void SetNoisePedestal(AliTPCCalPad *const noise, AliTPCCalPad *const pedestal=0x0)
                {fPadNoise=noise; fPedestals=pedestal;}
  void SetALTROData(AliTPCCalPad *const masked)
                {fALTROMasked=masked;}
  void SetGoofieArray(AliDCSSensorArray *const arr) {fGoofieArray=arr;}
  
  //setters for pad by pad information
  void SetRefFile(const char* filename);
  void SetReferenceRun(Int_t run=-1);
  void UpdateRefDataFromOCDB();
  void SetRefPulserData(AliTPCCalPad *const tmean, AliTPCCalPad *const trms=0x0, AliTPCCalPad *const qmean=0x0)
                {fRefPulserTmean=tmean; fRefPulserTrms=trms; fRefPulserQmean=qmean;}
  void SetRefCEData(AliTPCCalPad *const tmean, AliTPCCalPad *const trms=0x0, AliTPCCalPad *const qmean=0x0)
                {fRefCETmean=tmean; fRefCETrms=trms; fRefCEQmean=qmean;}
  void SetRefNoisePedestal(AliTPCCalPad *const noise, AliTPCCalPad *const pedestal=0x0)
                {fRefPadNoise=noise; fRefPedestals=pedestal;}
  void SetRefALTROData(AliTPCCalPad *const masked)
                {fRefALTROMasked=masked;}
  
  //creation of derived pad by pad calibration data
  AliTPCCalPad *CreatePadTime0(Int_t model, Double_t &gyA, Double_t &gyC, Double_t &chi2A, Double_t &chi2C);
  //
  // create outlyer maps
  //
  AliTPCCalPad *CreateCEOutlyerMap(Int_t &noutliersCE, AliTPCCalPad * const ceOut=0, Float_t minSignal=10, Float_t cutTrmsMin=0.9, Float_t cutTrmsMax=1.2, Float_t cutMaxDistT=0.7);
  AliTPCCalPad *CreatePulserOutlyerMap(Int_t &noutliersPulser, AliTPCCalPad * const pulserOut=0, Float_t cutTime=3, Float_t cutnRMSQ=5, Float_t cutnRMSrms=5);
  //
  AliTPCCalPad *CreatePadTime0CE(TVectorD &fitResultsA, TVectorD&fitResultsC, Int_t &nOut, Double_t &chi2A, Double_t &chi2C, const char *dumpfile=0);
  //

  void UpdatePulserOutlierMap();
  void UpdateRefPulserOutlierMap();
  void PulserOutlierMap(AliTPCCalPad *pulOut, const AliTPCCalPad *pulT, const AliTPCCalPad *pulQ);

  const char* GetGUIRefTreeDefaultName();
  
  Bool_t CreateGUIRefTree(const char* filename="");
  //
  // graph tools
  //
  static Double_t GetLaserTime0(Int_t run, Int_t timeStamp, Int_t deltaT, Int_t side);
  static TGraph* FilterGraphMedian(TGraph * const graph, Float_t sigmaCut, Double_t &medianY);
  static TGraph* FilterGraphMedianAbs(TGraph * graph, Float_t cut, Double_t &medianY);
  static TGraphErrors* FilterGraphMedianErr(TGraphErrors * graph, Float_t sigmaCut,Double_t &medianY);
  //
  static void Sort(TGraph *graph);
  static void SmoothGraph(TGraph *graph, Double_t delta);
  static Int_t     GetNearest(TGraph *graph, Double_t xref, Double_t &dx, Double_t &y);
  static Double_t EvalGraphConst(TGraph * const graph, Double_t xref);
  //
  // Filter sensors
  //
  static Float_t FilterSensor(AliDCSSensor * sensor, Double_t ymin, Double_t ymax, Double_t maxdy, Double_t sigmaCut); 
  //
  // Filter AliRelAlignmentKalman - Alignment/Drift velocity
  //
  static TMatrixD* MakeStatRelKalman(TObjArray * const array, Float_t minFraction, Int_t minStat, Float_t maxvd);
  static TObjArray *SmoothRelKalman(TObjArray * const array,const TMatrixD & stat, Bool_t direction, Float_t sigmaCut);
  static TObjArray *SmoothRelKalman(TObjArray * const arrayP, TObjArray * const arrayM);
  static void FilterCE(Double_t deltaT=100, Double_t cutAbs=10, Double_t cutSigma=4., TTreeSRedirector * const pcstream=0);
  static void FilterTracks(Int_t run, Double_t cutSigma=20., TTreeSRedirector * const pcstream=0);
  static Float_t FilterTemperature(AliTPCSensorTempArray *tempArray, Double_t ymin=15, Double_t ymax=22, Double_t sigmaCut=5); 

  static   void FilterGoofie(AliDCSSensorArray * goofieArray, Double_t deltaT=2, Double_t cutSigma=4.,  Double_t minVdn=8.5, Double_t maxVdn=9.05, TTreeSRedirector * const pcstream=0);
  static Double_t  GetTriggerOffsetTPC(Int_t run, Int_t timeStamp, Double_t deltaT=86400, Double_t deltaTLaser=3600, Int_t valType=0);
  static Double_t  GetVDriftTPC(Double_t &dist, Int_t run, Int_t timeStamp, Double_t deltaT=86400, Double_t deltaTLaser=3600, Int_t valType=0);
  static Double_t  GetVDriftTPCLaserTracks(Double_t &dist,Int_t run, Int_t timeStamp, Double_t deltaT=43200, Int_t side=2);
  static Double_t  GetVDriftTPCCE(Double_t &dist, Int_t run, Int_t timeStamp, Double_t deltaT=43200, Int_t side=2);
  static Double_t  GetVDriftTPCITS(Double_t &dist, Int_t run, Int_t timeStamp);
  static Double_t  GetTime0TPCITS(Double_t &dist, Int_t run, Int_t timeStamp);
  Int_t MakeRunList(Int_t startRun, Int_t stopRun); // find the list of usable runs
  Int_t FindRunTPC(Int_t    itime, Bool_t debug=kFALSE);
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
  //
  AliTPCdataQA   *fDataQA;            //data qa
  //reference data
  TMap *fRefMap;                        // latest map to reference information
  TMap *fCurrentRefMap;                 // reference data map of entries currently loaded
  TString fRefValidity;                 // validity range of reference data
  //  
  AliTPCCalPad  *fRefPadNoise;           //Reference noise information
  AliTPCCalPad  *fRefPedestals;          //Reference pedestal information
  AliTPCCalPad  *fRefPedestalMasked;     //Reference masked channels in pedestal run
  AliTPCCalPad  *fRefPulserTmean;        //Reference pulser mean time information
  AliTPCCalPad  *fRefPulserTrms;         //Reference pulser rms time information
  AliTPCCalPad  *fRefPulserQmean;        //Reference pulser mean q information
  AliTPCCalPad  *fRefPulserOutlier;      //Reference pulser outlier map
  AliTPCCalPad  *fRefPulserMasked;       //Reference masked channels in pulser run
  AliTPCCalPad  *fRefCETmean;            //Reference central electrode mean time information
  AliTPCCalPad  *fRefCETrms;             //Reference central electrode rms time information
  AliTPCCalPad  *fRefCEQmean;            //Reference central electrode mean q information
  AliTPCCalPad  *fRefCEMasked;           //Reference masked channels in laser run
  AliTPCCalPad  *fRefALTROFPED;          //Reference fixed pedestal value
  AliTPCCalPad  *fRefALTROZsThr;         //Reference zero suppression threshol
  AliTPCCalPad  *fRefALTROAcqStart;      //Reference accquistion start time bin
  AliTPCCalPad  *fRefALTROAcqStop;       //Reference accquistion stop time bin
  AliTPCCalPad  *fRefALTROMasked;        //Reference ALTRO masked channels information
  //
  AliTPCCalibRaw *fRefCalibRaw;          //Reference raw calibration object
  //
  AliTPCdataQA   *fRefDataQA;            //Reference data QA
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

  //
  // helpers to get the run number for given time stamps
  //
  // filters  
  
  TArrayI fRuns;                         // run list with OCDB info
  TArrayI fRunsStart;                    // start time for given run
  TArrayI fRunsStop;                     // stop time for given run
  
  AliTPCcalibDButil (const AliTPCcalibDButil& );
  AliTPCcalibDButil& operator= (const AliTPCcalibDButil& );

  AliTPCCalPad* GetRefCalPad(AliCDBEntry *entry, const char* objName);
  AliTPCCalPad* GetRefCalPad(AliCDBEntry *entry);
  AliTPCCalPad* GetAltroMasked(const char* cdbPath, const char* name);
  Bool_t HasRefChanged(const char *cdbPath);
  Int_t GetCurrentReferenceRun(const char* type) const;
  AliCDBEntry* GetRefEntry(const char* cdbPath);
  
  ClassDef(AliTPCcalibDButil,0)
};


#endif
