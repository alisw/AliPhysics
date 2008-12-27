#ifndef ALIEMCALQADataMakerRec_H
#define ALIEMCALQADataMakerRec_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.

  Based on PHOS code written by
  Y. Schutz CERN July 2007
*/


// --- ROOT system ---
class TH1F ; 
class TH1I ; 
class TObjArray ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQADataMakerRec.h"

class AliEMCALQADataMakerRec: public AliQADataMakerRec {

public:
  //Histograms for Raw data control
  enum HRawType_t {kNsmodLG,kNsmodHG,kLGtime,kHGtime,
		   kSpecLG,kSpecHG,kNtotLG,kNtotHG,
		   kEtotLG,kEtotHG} ;

  //Histograms for RecPoints  control
  enum HRPType_t {kRecPE,kRecPM,kRecPDigM};

  //Histograms for ESDs  control
  enum HESDType_t {kESDCaloClusE,kESDCaloClusM,kESDCaloCellA,kESDCaloCellM} ;
                 

public:
  AliEMCALQADataMakerRec() ;          // ctor
  AliEMCALQADataMakerRec(const AliEMCALQADataMakerRec& qadm) ;   
  AliEMCALQADataMakerRec& operator = (const AliEMCALQADataMakerRec& qadm) ;
  virtual ~AliEMCALQADataMakerRec() {;} // dtor
  
private:
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX_t, TObjArray ** list) ;
  virtual void   InitESDs() ; 
  virtual void   InitRecPoints() ; 
  virtual void   InitRaws() ; 
  virtual void   MakeESDs(AliESDEvent * esd) ;
  virtual void   MakeRecPoints(TTree * recpoTree) ; 
  virtual void   MakeRaws(AliRawReader* rawReader) ; 
  virtual void   StartOfDetectorCycle() ; 

  ClassDef(AliEMCALQADataMakerRec,1)  // description 

};

#endif // AliEMCALQADataMakerRec_H
