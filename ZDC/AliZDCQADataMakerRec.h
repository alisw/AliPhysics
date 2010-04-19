#ifndef AliZDCQADataMakerRec_H
#define AliZDCQADataMakerRec_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////
//  						    //
//  Produces the data needed for quality assurance. //
//  C. Oppedisano Chiara.Oppedisano@to.infn.it      //
//  						    //
//////////////////////////////////////////////////////


#include "AliQADataMakerRec.h"
#include "AliQAv1.h"

class AliCDBManager;
class AliCDBEntry;
class AliCDBStorage;
class AliZDCPedestals;

class AliZDCQADataMakerRec: public AliQADataMakerRec {

public:
  AliZDCQADataMakerRec() ;          // ctor
  AliZDCQADataMakerRec(const AliZDCQADataMakerRec& qadm) ;   
  AliZDCQADataMakerRec& operator = (const AliZDCQADataMakerRec& qadm) ;
  virtual ~AliZDCQADataMakerRec() {;} // dtor
  AliZDCPedestals  *GetPedCalibData() const; 
  
  
private:
  virtual void   InitESDs(); 
  virtual void   InitRecPoints();
  virtual void   InitRaws(); 
  virtual void   MakeRecPoints(TTree * /*recTree*/);
  virtual void   MakeRaws(AliRawReader* rawReader) ; 
  virtual void   MakeESDs(AliESDEvent * esd) ;
  virtual void   StartOfDetectorCycle(); 
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list) ;

  AliZDCPedestals *fPedCalibData; //! pedestal calibration data
 
  ClassDef(AliZDCQADataMakerRec,3)  // description 

};

#endif // AliZDCQADataMakerRec_H
