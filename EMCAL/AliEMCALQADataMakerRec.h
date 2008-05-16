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
  enum HRawType_t {kHGsmod1,kHGsmod2,kHGsmod3,kHGsmod4,kHGsmod5,kHGsmod6,
		   kHGsmod7,kHGsmod8,kHGsmod9,kHGsmod10,kHGsmod11,kHGsmod12,
		   kLGsmod1,kLGsmod2,kLGsmod3,kLGsmod4,kLGsmod5,kLGsmod6,
		   kLGsmod7,kLGsmod8,kLGsmod9,kLGsmod10,kLGsmod11,kLGsmod12,
		   kNsmodLG,kNsmodHG,
		   kNtotLG,kNtotHG,kEtotLG,kEtotHG,
		   kLGtime,kHGtime,kSpecLG,kSpecHG} ;
  //Histograms for RecPoints  control
  enum HRPType_t {kRPsmod1,kRPsmod2,kRPsmod3,kRPsmod4,kRPsmod5,kRPsmod6,
		  kRPsmod7,kRPsmod8,kRPsmod9,kRPsmod10,kRPsmod11,kRPsmod12,
		  kRPNtot,kRPEtot,kRPSpec,kRPTime} ;
  //Histograms for ESDs  control
  enum HESDType_t {kESDNtot,kESDEtot,kESDSpec,kESDpid} ;
                 

public:
  AliEMCALQADataMakerRec() ;          // ctor
  AliEMCALQADataMakerRec(const AliEMCALQADataMakerRec& qadm) ;   
  AliEMCALQADataMakerRec& operator = (const AliEMCALQADataMakerRec& qadm) ;
  virtual ~AliEMCALQADataMakerRec() {;} // dtor
  
private:
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX_t, TObjArray * list) ;
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
