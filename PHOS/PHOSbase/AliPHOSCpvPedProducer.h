#ifndef AliPHOSCpvPedProducer_h
#define AliPHOSCpvPedProducer_h

// Authors: Mikhail.Stolpovskiy@cern.ch, Sergey.Evdokimov@cern.ch
// The AliPHOSCpvPedProducer class calculates pedestals using AliPHOSCpvRawStream
// Also writes pedestals to files
// And creates a ROOT file with some histograms.
// this class supposed to be used in Cpv DA programm

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TFile.h>
#include <TString.h>
#include <THnSparse.h>
#include "AliPHOSCpvParam.h"
#include "AliPHOSCpvRawStream.h"
#include <TClonesArray.h>

class TFile;
class AliPHOSCpvPedProducer: public TObject {


public:
  AliPHOSCpvPedProducer(Int_t sigcut = 3);
  virtual ~AliPHOSCpvPedProducer();
  void SetPermanentBadMap(TH2* badMap, int iDDL);
  void   SetSigCut(Int_t sigcut = 3) {fSigCut = sigcut;} //set n. of pedestal distribution sigmas used to create zero suppresion table
  Bool_t LoadNewEvent(AliRawReader *& rawReader); // returns true, if ok
  void   SetTurbo(Bool_t turbo);                  // if turbo==true then do read without error checking
  Bool_t GetTurbo() const {return fTurbo;}

  Bool_t FillPedestal(Int_t pad,Float_t q);   // pad - absolute pad number; q - charge of the pad
  Bool_t FillPedestal();  // analyse event and fill pedestals
  Bool_t CalcPedestal(Int_t iDDL);        //  analyse pedestals when all events processed for indicated DDL

  TH2F* GetPedMeanMap(Int_t iDDL) const {return fPedMeanMap[iDDL];}  //Get the pedestal mean map for a DDL to send to AMORE
  TH2F* GetPedSigMap(Int_t iDDL)  const {return fPedSigMap[iDDL];}   //Get the pedestal sigma map for a DDL to send to AMORE
  TH1F* GetPedMean(Int_t iChFee)  const {return f1DPedMean[iChFee];} //Get the pedestal mean map for a FEE channel to send to AMORE
  TH1F* GetPedSigma(Int_t iChFee) const {return f1DPedSigma[iChFee];}//Get the pedestal Sigma map for a FEE channel to send to AMORE
  void  WritePedFiles(Int_t iDDL) const;                             // write pedestals to load(?) to RCB for indicated DDL
  void WriteAllHistsToFile(const char* name) const;                  // create and write a new root file with hists

  void SetErrorsHist(TH1I * pHist) {fhErrors = new TH1I(*pHist);}    //Set the histogram of errors, taken from AliPHOSCpvRawDigiProdicer
  static Int_t GetMaxThr() {return fMaxThr;}                         // Get maximal threshold

protected:
  void CreateErrHist();    // initialize histogram of errors
  void CreateDDLHistos(Int_t iDDL);  // initialize histograms for pedestal calculation and representation

  //our ddls are 0,2,4,6,8
  TH1F       *fPadAdc[2*AliPHOSCpvParam::kNDDL][AliPHOSCpvParam::kPadPcX][AliPHOSCpvParam::kPadPcY];        //Charge distribution for pads
  Int_t       fSigCut;                         //n. of pedestal distribution sigmas used to create zero suppresion table
  static const Int_t fMaxThr = 511;            //maximal threshold (9 bits all with 1)
  Bool_t fTurbo;           // if true, then read without error checking
  //our ddls are 0,2,4,6,8
  TH2F       *fPedMeanMap[2*AliPHOSCpvParam::kNDDL]; //2D mean pedestal map to export to AMORE
  TH2F       *fPedSigMap [2*AliPHOSCpvParam::kNDDL]; //2D pedestal sigma map to export to AMORE
  TH1F       *f1DPedMean [2*AliPHOSCpvParam::kNDDL]; //1D mean pedestal map to export to AMORE
  TH1F       *f1DPedSigma[2*AliPHOSCpvParam::kNDDL]; //1D pedestal sigma map to export to AMORE
  TH1I       *fhErrors;                        //histogram of errors from AliPHOSCpvRawDigiProducer
  TH2I       *fPermanentBadMap[2*AliPHOSCpvParam::kNDDL];
  AliPHOSCpvRawStream         * fRawStream;       //! Raw data stream
private:
  ClassDef(AliPHOSCpvPedProducer,1);                                                  //Cpv calibration and pedestal class
};
#endif
