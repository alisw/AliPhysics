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
  
  virtual void ResetDetector(AliQAv1::TASKINDEX_t task);

  Int_t  GetRawMaxEvents()      const { return fRawMaxEvents;     }
  Int_t  GetRawEventsPerBin()   const { return fRawEventsPerBin;  }
  Int_t  GetRawFirstTimeBin() const { return fRawFirstTimeBin; }
  Int_t  GetRawLastTimeBin()  const { return fRawLastTimeBin;  }

  void  SetRawMaxEvents   (Int_t value) { fRawMaxEvents = value; }
  void  SetRawEventsPerBin(Int_t value) { fRawEventsPerBin = value; }
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

  TH1F* CreateEventsHistCopy(const TH1F* hist, const Char_t* copyName);

  AliTPCAltroMapping *fMapping[6]; //! Pointers to ALTRO mapping
  AliTPCdataQA** fTPCdataQA;//! TPC calibration object for making raw data QA

  Int_t fRawMaxEvents;      //! Max events for RAW QA event histograms
  Int_t fRawEventsPerBin;   //! Events per bin for RAW QA event histograms
  Int_t fRawFirstTimeBin;   //! First Time bin needed for RAW QA
  Int_t fRawLastTimeBin;    //! Last Time bin needed for RAW QA

  ClassDef(AliTPCQADataMakerRec,1)  // TPC Rec Quality Assurance Data Maker 
};

#endif // ALITPCQADATAMAKERREC_H
