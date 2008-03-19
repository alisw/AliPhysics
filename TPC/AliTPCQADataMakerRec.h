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

#include <AliTPCdataQA.h>

class AliTPCQADataMakerRec: public AliQADataMakerRec {

public:
  AliTPCQADataMakerRec() ;          // ctor
  AliTPCQADataMakerRec(const AliTPCQADataMakerRec& qadm) ;   
  AliTPCQADataMakerRec& operator = (const AliTPCQADataMakerRec& qadm) ;
  virtual ~AliTPCQADataMakerRec() { delete fTPCdataQA; } // dtor
  
private:
  virtual void   StartOfDetectorCycle() {}; // empty 
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX_t, TObjArray* list) ;

  // ESD QA
  virtual void   InitESDs() ; 
  virtual void   MakeESDs(AliESDEvent *esd) ;

  // Raw QA
  virtual void   InitRaws();
  virtual void   MakeRaws(AliRawReader* rawReader);

  // RecPoints QA
  virtual void   InitRecPoints();
  virtual void   MakeRecPoints(TTree *recTree);

  AliTPCdataQA* fTPCdataQA;//! TPC calibration object for making raw data QA

  TH1F* fHistESDclusters;  //! Clusters per ESD track
  TH1F* fHistESDratio;     //! Ratio of clusters to findables
  TH1F* fHistESDpt;        //! Pt spectrum
  
  TH1F* fHistRawsOccupancy;//! Pad occupancy (1 entry per pad)

  TH1F* fHistRecPointsQmaxShort; //! Qmax (short pads)
  TH1F* fHistRecPointsQmaxMedium;//! Qmax (medium pads)
  TH1F* fHistRecPointsQmaxLong;  //! Qmax (long pads)
  TH1F* fHistRecPointsQShort;    //! Q (short pads)
  TH1F* fHistRecPointsQMedium;   //! Q (medium pads)
  TH1F* fHistRecPointsQLong;     //! Q (long pads)
  TH1F* fHistRecPointsRow;       //! Row distribution

  ClassDef(AliTPCQADataMakerRec,1)  // TPC Rec Quality Assurance Data Maker 
};

#endif // ALITPCQADATAMAKERREC_H
