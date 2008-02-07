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
  Int_t finderIO();//for MC
  Int_t finderIO(AliRawReader* rawReader);//for data
  Int_t findClusterKrIO();//main routine for finding clusters

  virtual void SetInput(TTree * tree);  //set input tree with digits    
  virtual void SetOutput(TTree * tree); //set output tree with clusters

  void SetParam(AliTPCParam *param){fParam=param;}//set TPC parameters
  void SetDigArr(AliTPCDigitsArray *digarr){fDigarr=digarr;}//set current array of digits
  void SetRecoParam(AliTPCRecoParam *recoParam=0);//set reconstruction parameters
  virtual void SetOldRCUFormat(Bool_t rcuFormat = kFALSE)
    { fIsOldRCUFormat = rcuFormat; };

  Bool_t fRawData; //flague =0 for MC =1 for real data
  AliTPCClustersRow * fRowCl;  //! current cluster row (used in rootuple fill)

private:
  TTree * fInput;   //!input  tree with digits - object not owner
  TTree * fOutput;   //!output tree with clusters - object not owner
  AliTPCParam * fParam;//!TPC parameters
  AliTPCDigitsArray *fDigarr;//! pointer to current array if digits

  //only for raw data :)
  const AliTPCRecoParam  * fRecoParam;        //! reconstruction parameters
  Bool_t fIsOldRCUFormat; // assume old RCU raw data format
  

  ClassDef(AliTPCclustererKr,1)  // Time Projection Chamber Kr clusters
};


#endif


