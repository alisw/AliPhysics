#ifndef ALITPCQADATAMAKERREC_H
#define ALITPCQADATAMAKERREC_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id: $ */

/*
  Based on AliPHOSQADataMaker
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  P. Christiansen, Lund, January 2008
*/


// --- ROOT system ---
#include <TH1.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include <AliQADataMakerRec.h>
#include <AliRawReader.h>
#include <AliTPCAltroMapping.h>

#include <AliTPCdataQA.h>

class AliTPCQADataMakerRec: public AliQADataMakerRec {

public:
  enum HRawsType_t         {kRawsOccupancy=0, kRawsOccupancyVsSector, kRawsNClustersPerEventVsSector, kRawsQVsSector, kRawsQmaxVsSector, kRawsOccupancyVsEvent, kRawsNclustersVsEvent} ; 
  enum HDigitType_t        {kDigitsADC=0} ; 
  enum HRECPOINTsType_t    {KClusters=0, kRatio, kPt} ; 
  enum HESDsType_t         {kQmaxShort=0, kQmaxMedium, kQmaxLong, kQShort, kQMedium, kQLong, kRow} ; 

  AliTPCQADataMakerRec() ;          // ctor
  AliTPCQADataMakerRec(const AliTPCQADataMakerRec& qadm) ;   
  AliTPCQADataMakerRec& operator = (const AliTPCQADataMakerRec& qadm) ;
  virtual ~AliTPCQADataMakerRec(); 
  
  void SetBeautifyOption(Int_t value)  {fBeautifyOption= value;}
  void SetOccHighLimit(Float_t value)  {fOccHighLimit  = value;}
  void SetQmaxLowLimit(Float_t value)  {fQmaxLowLimit  = value;}
  void SetQmaxHighLimit(Float_t value) {fQmaxHighLimit = value;}

  Int_t   GetBeautifyOption() const {return fBeautifyOption;}
  Float_t GetOccHighLimit() const {return fOccHighLimit; }
  Float_t GetQmaxLowLimit() const {return fQmaxLowLimit; }
  Float_t GetQmaxHighLimit() const {return fQmaxHighLimit;}

private:
  virtual void   StartOfDetectorCycle() {}; // empty 
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray** list) ;

  // ESD QA
  virtual void   InitESDs() ; 
  virtual void   MakeESDs(AliESDEvent *esd) ;

  // Raw QA
  virtual void   InitRaws();
  virtual void   MakeRaws(AliRawReader* rawReader);

  // Digits QA
  virtual void   InitDigits();
  virtual void   MakeDigits()  {return;}
  virtual void   MakeDigits(TTree *digTree);
  
  // RecPoints QA
  virtual void   InitRecPoints();
  virtual void   MakeRecPoints(TTree *recTree);
  
  virtual void LoadMaps();

  AliTPCAltroMapping *fMapping[6]; //! Pointers to ALTRO mapping
  AliTPCdataQA** fTPCdataQA;//! TPC calibration object for making raw data QA

  Int_t   fBeautifyOption;//! 0:no beautify, !=0:beautify RAW 
  Float_t fOccHighLimit;  //! high limit for accepting occupancy values
  Float_t fQmaxLowLimit;  //! low limit for accepting Qmax values
  Float_t fQmaxHighLimit; //! high limit for accepting Qmax values
  
  ClassDef(AliTPCQADataMakerRec,1)  // TPC Rec Quality Assurance Data Maker 
};

#endif // ALITPCQADATAMAKERREC_H
