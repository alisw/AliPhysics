#ifndef AliPHOSCpvGainCalibDA_h
#define AliPHOSCpvGainCalibDA_h

// Author: Sergey Evdokimov <sevdokim@cern.ch>
// The AliPHOSCpvGainCalibDA class creates 1D amplitude histos for every channel using AliPHOSDigit array
// produced by AliPHOSCPVRawDigiProducer (with pedestal subtraction!!!)
// Also writes these histos to file  
// And creates a ROOT file with some histograms

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TFile.h>
#include <TString.h>
#include <THnSparse.h>
#include "AliPHOSCpvParam.h"
#include "AliPHOSCpvRawDigiProducer.h"
#include <TClonesArray.h>
#include "AliPHOSGeometry.h"

class TFile;
class AliPHOSCpvGainCalibDA: public TObject { 


public:
  AliPHOSCpvGainCalibDA();
  virtual ~AliPHOSCpvGainCalibDA();
  void InitCalibration(TFile *fCalibrSupplyRoot); //run it before analysing data to create calibration coeffs map. input is file where previously created histos are stored.
  void CreateA0Histos(Int_t iDDL);
  Bool_t SetDeadChannelMapFromFile(const char * filename = "CpvBadMap.root");
  void WriteA0HistosToFile(const char * filename=0x0) const;                  // create and write a new CpvCalibrSupply.root file with hists
  Bool_t IsBad(Int_t ddl, Int_t x, Int_t y) {  // returns true, if the cell is bad
    if(ddl<0||ddl>2*AliPHOSCpvParam::kNDDL) return kTRUE;
    if(!fDeadMap[ddl]) return kFALSE;
    if(fDeadMap[ddl] -> GetBinContent(x+1,y+1)) return kTRUE;
    return kFALSE;
  }
  Bool_t IsBad(Int_t abs) {  // returns true, if the cell is bad
    if(!AliPHOSCpvParam::IsValidAbs(abs)) return kFALSE;
    Int_t ddl = AliPHOSCpvParam::A2DDL(abs),
            x = AliPHOSCpvParam::A2X(abs),
            y = AliPHOSCpvParam::A2Y(abs);
    return IsBad(ddl,x,y);
  }
  Bool_t FillAmplA0Histos(TClonesArray *digits);
  TList * GetQAHistos(){return fHistosList;}
  void CreateQAHistos();
  void SetMinClustSize(Int_t a){fMinClustSize=a;}

 protected: 
  Int_t fMinClustSize;//minimum cluster size
  AliPHOSGeometry * fGeom ;         //! PHOS geometry
  TH2I *fDeadMap[2*AliPHOSCpvParam::kNDDL]; //Dead Channel Map  
  TH2I *fEntriesMap[2*AliPHOSCpvParam::kNDDL];//entries map
  TH1F       *fAmplA0Histo[2*AliPHOSCpvParam::kNDDL][AliPHOSCpvParam::kPadPcX][AliPHOSCpvParam::kPadPcY]; //raw amplitudes specrta for every channel   
  TList* fHistosList;
  TH1F* fhClusterMult;
  TH2F* fhClusterShape;
  TH1F* fhA0Value;
  TH2F* fhAmplInClust;
  TH1F* fhTotalClusterAmplitude;
  
  ClassDef(AliPHOSCpvGainCalibDA,1);                                                  //Cpv calibration class        
};
#endif

