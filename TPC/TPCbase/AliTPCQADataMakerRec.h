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

// --- Standard library ---

// --- AliRoot header files ---
#include <AliQADataMakerRec.h>
#include <AliRawReader.h>
#include <AliTPCAltroMapping.h>

#include <AliTPCdataQA.h>

class AliTPCQADataMakerRec: public AliQADataMakerRec {

public:
  enum HRawsType_t         {kRawsOccupancyVsSector=0, kRawsQVsSector, kRawsQmaxVsSector, kRawsOccupancy2dVsSector} ; 
  enum HDigitType_t        {kDigitsADC=0} ; 
  enum HRECPOINTsType_t    {kClusters=0, kRatio, kPt} ; 
  enum HESDsType_t         {kQmaxShort=0, kQmaxMedium, kQmaxLong, kQShort, kQMedium, kQLong, kRow} ; 

  AliTPCQADataMakerRec() ;          // ctor
  AliTPCQADataMakerRec(const AliTPCQADataMakerRec& qadm) ;   
  AliTPCQADataMakerRec& operator = (const AliTPCQADataMakerRec& qadm) ;
  virtual ~AliTPCQADataMakerRec(); 
  
  Int_t  GetRawFirstTimeBin() const { return fRawFirstTimeBin; }
  Int_t  GetRawLastTimeBin()  const { return fRawLastTimeBin;  }
  
  void  SetRawRangeTime(Int_t tMin, Int_t tMax){ fRawFirstTimeBin=tMin; fRawLastTimeBin=tMax;}
  
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
  AliTPCdataQA*  fTPCdataQA;//! TPC calibration object for making raw data QA
  
  Int_t fRawFirstTimeBin;   //! First Time bin needed for RAW QA
  Int_t fRawLastTimeBin;    //! Last Time bin needed for RAW QA
  
  ClassDef(AliTPCQADataMakerRec,1)  // TPC Rec Quality Assurance Data Maker 
    };

#endif // ALITPCQADATAMAKERREC_H
