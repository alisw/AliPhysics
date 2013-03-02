#ifndef ALITPCCLUSTERERKR_H
#define ALITPCCLUSTERERKR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCclustererKr.h,v 1.8 2008/02/07 16:07:15 matyja Exp $ */

//-------------------------------------------------------
//                    TPC Kr Clusterer Class
//
//   Origin: Adam Matyja, INP PAN, adam.matyja@ifj.edu.pl
//-------------------------------------------------------

#include "TObject.h"

class AliTPCclusterKr;
class AliPadMax;
class AliSimDigits;
class AliTPCv4;
class AliTPCParam;
class AliTPCDigitsArray;
class AliTPCvtpr;
class AliTPCClustersRow;
class TTree;
class TH1F;
class TH2F;

class AliTPCTransform;

//used in raw data finder
class AliTPCROC;
class AliTPCCalPad;
class AliTPCAltroMapping;
class AliTPCcalibDB;
class AliTPCRecoParam;
class AliTPCReconstructor;
class AliRawReader;
class AliTPCCalROC;
class TTreeSRedirector;
class AliRawEventHeaderBase;
//_____________________________________________________________________________
class AliTPCclustererKr: public TObject{
public:
  AliTPCclustererKr();
  AliTPCclustererKr(const AliTPCclustererKr &param);//copy constructor
  AliTPCclustererKr &operator = (const AliTPCclustererKr & param); 
  virtual ~AliTPCclustererKr();

  //finders
  virtual Int_t FinderIO();//for MC
  virtual Int_t FinderIO(AliRawReader* rawReader);//for data
  virtual Int_t FindClusterKrIO();//main routine for finding clusters
  virtual void CleanSector(Int_t sector); // clean isolated digits

  //other
  void GetXY(Int_t sec,Int_t row,Int_t pad,Double_t& xGlob,Double_t& yGlob);//give XY coordinate of the pad

  virtual void SetInput(TTree * tree);  //set input tree with digits    
  virtual void SetOutput(TTree * tree); //set output tree with clusters
  virtual void SetParam(AliTPCParam *param){fParam=param;}//set TPC parameters
  virtual void SetDigArr(AliTPCDigitsArray *digarr){fDigarr=digarr;}//set current array of digits
  virtual void SetRecoParam(AliTPCRecoParam *recoParam=0);//set reconstruction parameters

  virtual void SetTimeStamp(UInt_t timestamp){ fTimeStamp = timestamp; }
  virtual void SetRun(UInt_t run){ fRun = run; }

  //setters for cluster finder parameters
  virtual void SetZeroSup(Int_t v){fZeroSup=v;}//set zero suppresion parameter
  virtual void SetFirstBin(Int_t v){fFirstBin=v;}//set first considered timebin
  virtual void SetLastBin(Int_t v){fLastBin=v;}//set last considered timebin
  virtual void SetMaxNoiseAbs(Float_t v){fMaxNoiseAbs=v;}//set maximal noise value
  virtual void SetMaxNoiseSigma(Float_t v){fMaxNoiseSigma=v;}//set maximal noise sigma

  virtual void SetMinAdc(Int_t v){v<=0?fMinAdc=1:fMinAdc=v;}//set fMinAdc
  virtual void SetMinTimeBins(Int_t v){fMinTimeBins=v;}//set fMinTimeBins
//  virtual void SetMaxPadRange(Int_t v){fMaxPadRange=v;}//set fMaxPadRange
//  virtual void SetMaxRowRange(Int_t v){fMaxRowRange=v;}//set fMaxRowRange
  virtual void SetMaxTimeRange(Int_t v){fMaxTimeRange=v;}//set fMaxTimeRange
  virtual void SetValueToSize(Float_t v){fValueToSize=v;}//set fValueToSize

  virtual void SetMaxPadRangeCm(Double_t v){fMaxPadRangeCm=v;}//set fMaxPadRangeCm
  virtual void SetMaxRowRangeCm(Double_t v){fMaxRowRangeCm=v;}//set fMaxRowRangeCm

  virtual void SetIsolCut(Short_t v){fIsolCut=v;}

  virtual void SetDebugLevel(Int_t debug){fDebugLevel=debug;}
  //debug = 0 to 71 -sector number to  print
  // = 72 - all sectors
  // = 73 - inners
  // = 74 - outers

  virtual void SetHistoRow(TH1F *histo)   {fHistoRow   =histo;}
  virtual void SetHistoPad(TH1F *histo)   {fHistoPad   =histo;}
  virtual void SetHistoTime(TH1F *histo)  {fHistoTime  =histo;}
  virtual void SetHistoRowPad(TH2F *histo){fHistoRowPad=histo;}

  //getters for cluster finder parameters
  Int_t GetZeroSup() const {return fZeroSup;}//get zero suppresion parameter
  Int_t GetFirstBin() const {return fFirstBin;}//get first considered timebin
  Int_t GetLastBin() const {return fLastBin;}//get last considered timebin
  Float_t GetMaxNoiseAbs() const {return fMaxNoiseAbs;}//get maximal noise value
  Float_t GetMaxNoiseSigma() const {return fMaxNoiseSigma;}//get maximal noise sigma

  Int_t GetMinAdc() const {return fMinAdc;}//get fMinAdc
  Int_t GetMinTimeBins() const {return fMinTimeBins;}//get fMinTimeBins
//  Int_t GetMaxPadRange() const {return fMaxPadRange;}//get fMaxPadRange
//  Int_t GetMaxRowRange() const {return fMaxRowRange;}//get fMaxRowRange
  Int_t GetMaxTimeRange() const {return fMaxTimeRange;}//get fMaxTimeRange
  Float_t GetValueToSize() const {return fValueToSize;}//get fValueToSize

  Double_t GetMaxPadRangeCm() const {return fMaxPadRangeCm;}//get fMaxPadRangeCm
  Double_t GetMaxRowRangeCm() const {return fMaxRowRangeCm;}//get fMaxRowRangeCm

  Short_t GetIsolCut() const {return fIsolCut;}

  Int_t GetDebugLevel() const {return fDebugLevel;}
  TH1F * GetHistoRow(){return fHistoRow;} 
  TH1F * GetHistoPad(){return fHistoPad;}
  TH1F * GetHistoTime(){return fHistoTime;}
  TH2F * GetHistoRowPad(){return fHistoRowPad;}

  UInt_t GetTimeStamp() const {return fTimeStamp;}
  UInt_t GetRun() const {return fRun;}

private:
  void MakeClusters(TObjArray * maximaInSector, Int_t iSec, Int_t &clusterCounter);
  Bool_t fRawData; //flag =0 for MC =1 for real data

  TTree * fInput;   //!input  tree with digits - object not owner
  TTreeSRedirector * fOutput;   //!output tree with clusters - object not owner

  AliTPCParam * fParam;//!TPC parameters
  AliTPCDigitsArray *fDigarr;//! pointer to current array if digits

  //only for raw data :)
  const AliTPCRecoParam  * fRecoParam;        //! reconstruction parameters

  //cluster finder parameters
  Int_t fZeroSup;//zero suppresion parameter = 2 def.
  Int_t fFirstBin;//first considered time bin used by cluster finder = 60 def.
  Int_t fLastBin;//last considered time bin used by cluster finder = 950 def.
  Float_t fMaxNoiseAbs;// maximal noise value on pad used in cluster finder = 2 def.
  Float_t fMaxNoiseSigma;// maximal noise sigma on pad used in cluster finder = 3 def.

  Int_t fMinAdc;//minimal value of acd count in each timebin = 3 def.
  Int_t fMinTimeBins;//minimal value of time bins one after each other = 2 def.
//  Int_t fMaxPadRange;//maximal pad range from maximum = 4 def.
//  Int_t fMaxRowRange;//maximal row range from maximum = 3 def.
  Int_t fMaxTimeRange;//maximal time bin range from maximum = 7 def.
  Float_t fValueToSize;//ratio cluster value to cluster size = 3.5 def.

  Double_t fMaxPadRangeCm;//maximal pad range in cm from maximum = 2.5cm def.
  Double_t fMaxRowRangeCm;//maximal row range in cm from maximum = 3.5cm def.

  Short_t fIsolCut;//isolation cut in 3D = 5 def.

  Int_t fDebugLevel;//! debug level variable
  TH1F *fHistoRow;//!debug histo for rows
  TH1F *fHistoPad;//!debug histo for pads
  TH1F *fHistoTime;//!debug histo for timebins
  TH2F *fHistoRowPad;//!debug histo for rows and pads

  UInt_t fTimeStamp; //!time stamp from event header
  UInt_t fRun; //!run from event header
  ClassDef(AliTPCclustererKr,8)  // Time Projection Chamber Kr clusters
};


#endif
