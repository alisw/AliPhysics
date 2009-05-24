#ifndef ALIPMDQADataMakerRec_H
#define ALIPMDQADataMakerRec_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  B.K. Nandi
*/


// --- ROOT system ---
class TH1F ; 
class TH1I ; 
class TObjArray ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQADataMakerRec.h"

class AliPMDQADataMakerRec: public AliQADataMakerRec {

 public:
    AliPMDQADataMakerRec() ;          // ctor
    AliPMDQADataMakerRec(const AliPMDQADataMakerRec& qadm) ;   
    AliPMDQADataMakerRec& operator = (const AliPMDQADataMakerRec& qadm) ;
    virtual ~AliPMDQADataMakerRec() {;} // dtor
    
 private:
    
    virtual void   InitRaws() ; 
    virtual void   InitDigits() ; 
    virtual void   InitRecPoints() ; 
    virtual void   InitESDs() ; 
    virtual void   MakeRaws(AliRawReader* rawReader) ; 
    virtual void   MakeDigits(TClonesArray* digits)  ; 
    virtual void   MakeDigits(TTree * recpoTree) ; 
    virtual void   MakeRecPoints(TTree * recpoTree) ; 
    virtual void   MakeESDs(AliESDEvent * esd) ;
    virtual void   StartOfDetectorCycle() ; 
    virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** list) ;
    
    ClassDef(AliPMDQADataMakerRec,1)  // description 

};

#endif // AliPMDQADataMakerRec_H
