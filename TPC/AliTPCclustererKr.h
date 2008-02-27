#ifndef ALITPCCLUSTERERKR_H
#define ALITPCCLUSTERERKR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCclusterKr.h,v 1.8 2008/02/07 16:07:15 matyja Exp $ */

//-------------------------------------------------------
//                    TPC Kr Cluster Class
//
//   Origin: Adam Matyja, INP PAN, adam.matyja@ifj.edu.pl
//-------------------------------------------------------

#include "AliTPCclusterKr.h"
#include <vector>
#include "TObject.h"
#include "AliPadMax.h"
#include "AliSimDigits.h"
#include "AliTPCv4.h"
#include "AliTPCParam.h"
#include "AliTPCDigitsArray.h"
#include "AliTPCvtpr.h"
#include "AliTPCClustersRow.h"
#include "TTree.h"

//used in raw data finder
#include "AliTPCROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCAltroMapping.h"
#include "AliTPCcalibDB.h"
#include "AliTPCRawStream.h"
#include "AliTPCRecoParam.h"
#include "AliTPCReconstructor.h"
#include "AliRawReader.h"
#include "AliTPCCalROC.h"

//_____________________________________________________________________________
class AliTPCclustererKr: public TObject{
public:
  AliTPCclustererKr();
  AliTPCclustererKr(const AliTPCclustererKr &param);//copy constructor
  AliTPCclustererKr &operator = (const AliTPCclustererKr & param); 
  ~AliTPCclustererKr();

  //finders
  Int_t FinderIO();//for MC
  Int_t FinderIO(AliRawReader* rawReader);//for data
  Int_t FindClusterKrIO();//main routine for finding clusters

  //other
  void GetXY(Short_t sec,Short_t row,Short_t pad,Double_t& xGlob,Double_t& yGlob);//give XY coordinate of the pad

  virtual void SetInput(TTree * tree);  //set input tree with digits    
  virtual void SetOutput(TTree * tree); //set output tree with clusters

  void SetParam(AliTPCParam *param){fParam=param;}//set TPC parameters
  void SetDigArr(AliTPCDigitsArray *digarr){fDigarr=digarr;}//set current array of digits
  void SetRecoParam(AliTPCRecoParam *recoParam=0);//set reconstruction parameters
  virtual void SetOldRCUFormat(Bool_t rcuFormat = kFALSE)
    { fIsOldRCUFormat = rcuFormat; };

  Bool_t fRawData; //flague =0 for MC =1 for real data
  AliTPCClustersRow * fRowCl;  //! current cluster row (used in rootuple fill)

  //setters for cluster finder parameters
  void SetZeroSup(Int_t v){fZeroSup=v;}//set zero suppresion parameter
  void SetFirstBin(Short_t v){fFirstBin=v;}//set first considered timebin
  void SetLastBin(Short_t v){fLastBin=v;}//set last considered timebin
  void SetMaxNoiseAbs(Float_t v){fMaxNoiseAbs=v;}//set maximal noise value
  void SetMaxNoiseSigma(Float_t v){fMaxNoiseSigma=v;}//set maximal noise sigma

  void SetMinAdc(Short_t v){v<=0?fMinAdc=1:fMinAdc=v;}//set fMinAdc
  void SetMinTimeBins(Short_t v){fMinTimeBins=v;}//set fMinTimeBins
//  void SetMaxPadRange(Short_t v){fMaxPadRange=v;}//set fMaxPadRange
//  void SetMaxRowRange(Short_t v){fMaxRowRange=v;}//set fMaxRowRange
  void SetMaxTimeRange(Short_t v){fMaxTimeRange=v;}//set fMaxTimeRange
  void SetValueToSize(Float_t v){fValueToSize=v;}//set fValueToSize

  void SetMaxPadRangeCm(Double_t v){fMaxPadRangeCm=v;}//set fMaxPadRangeCm
  void SetMaxRowRangeCm(Double_t v){fMaxRowRangeCm=v;}//set fMaxRowRangeCm

  //getters for cluster finder parameters
  Int_t GetZeroSup(){return fZeroSup;}//get zero suppresion parameter
  Short_t GetFirstBin(){return fFirstBin;}//get first considered timebin
  Short_t GetLastBin(){return fLastBin;}//get last considered timebin
  Float_t GetMaxNoiseAbs(){return fMaxNoiseAbs;}//get maximal noise value
  Float_t GetMaxNoiseSigma(){return fMaxNoiseSigma;}//get maximal noise sigma

  Short_t GetMinAdc(){return fMinAdc;}//get fMinAdc
  Short_t GetMinTimeBins(){return fMinTimeBins;}//get fMinTimeBins
//  Short_t GetMaxPadRange(){return fMaxPadRange;}//get fMaxPadRange
//  Short_t GetMaxRowRange(){return fMaxRowRange;}//get fMaxRowRange
  Short_t GetMaxTimeRange(){return fMaxTimeRange;}//get fMaxTimeRange
  Float_t GetValueToSize(){return fValueToSize;}//get fValueToSize

  Double_t GetMaxPadRangeCm(){return fMaxPadRangeCm;}//get fMaxPadRangeCm
  Double_t GetMaxRowRangeCm(){return fMaxRowRangeCm;}//get fMaxRowRangeCm

private:
  TTree * fInput;   //!input  tree with digits - object not owner
  TTree * fOutput;   //!output tree with clusters - object not owner
  AliTPCParam * fParam;//!TPC parameters
  AliTPCDigitsArray *fDigarr;//! pointer to current array if digits

  //only for raw data :)
  const AliTPCRecoParam  * fRecoParam;        //! reconstruction parameters
  Bool_t fIsOldRCUFormat; // assume old RCU raw data format

  //cluster finder parameters
  Int_t fZeroSup;//zero suppresion parameter = 2 def.
  Short_t fFirstBin;//first considered time bin used by cluster finder = 60 def.
  Short_t fLastBin;//last considered time bin used by cluster finder = 950 def.
  Float_t fMaxNoiseAbs;// maximal noise value on pad used in cluster finder = 2 def.
  Float_t fMaxNoiseSigma;// maximal noise sigma on pad used in cluster finder = 3 def.

  Short_t fMinAdc;//minimal value of acd count in each timebin = 3 def.
  Short_t fMinTimeBins;//minimal value of time bins one after each other = 2 def.
//  Short_t fMaxPadRange;//maximal pad range from maximum = 4 def.
//  Short_t fMaxRowRange;//maximal row range from maximum = 3 def.
  Short_t fMaxTimeRange;//maximal time bin range from maximum = 7 def.
  Float_t fValueToSize;//ratio cluster value to cluster size = 3.5 def.

  Double_t fMaxPadRangeCm;//maximal pad range in cm from maximum = 2.5cm def.
  Double_t fMaxRowRangeCm;//maximal row range in cm from maximum = 3.5cm def.

  ClassDef(AliTPCclustererKr,2)  // Time Projection Chamber Kr clusters
};


#endif


