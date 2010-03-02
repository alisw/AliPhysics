#ifndef ALITOFQADATAMAKERREC_H
#define ALITOFQADATAMAKERREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////
//                                                                // 
//  Produces the data needed to calculate the quality assurance.  //
//    All data must be mergeable objects.                         //
//    S. Arcelli                                                  //
//                                                                //
//    /* last modified by F. Bellini on 25/02/2010 */             // 
////////////////////////////////////////////////////////////////////


#include "AliQADataMakerRec.h"
#include "AliQAv1.h"
#include "AliTOFChannelOnlineStatusArray.h"

class AliCDBManager;
class AliCDBEntry;
class AliCDBStorage;
class AliTOFChannelOnlineStatusArray;

class AliTOFQADataMakerRec: public AliQADataMakerRec {

public:
  AliTOFQADataMakerRec() ;          // ctor
  AliTOFQADataMakerRec(const AliTOFQADataMakerRec& qadm) ;   
  AliTOFQADataMakerRec& operator = (const AliTOFQADataMakerRec& qadm) ;
  AliTOFChannelOnlineStatusArray *GetCalibData() const;
  virtual ~AliTOFQADataMakerRec() {;} // dtor
  
protected: 
  AliTOFChannelOnlineStatusArray * fCalibData;        //! calibration data

private:
  virtual void   InitESDs() ; 
  virtual void   InitRecPoints() ; 
  virtual void   InitRaws() ; 
  virtual void   MakeESDs(AliESDEvent * esd) ;
  virtual void   MakeRecPoints(TTree * recTree) ; 
  virtual void   MakeRaws(AliRawReader* rawReader) ; 
  virtual void   StartOfDetectorCycle() ; 
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list) ;
  virtual void   GetMapIndeces(Int_t *in, Int_t *out) ; 
  virtual void   EnableNoiseFiltering(Bool_t enable){fEnableNoiseFiltering = enable;};
  Bool_t   CheckVolumeID(Int_t *equipmentID); 
  Bool_t   CheckEquipID(Int_t *equipmentID); 

  Bool_t fEnableNoiseFiltering; //the choice is not implemented so far
  ClassDef(AliTOFQADataMakerRec,2)  // description 

};

#endif // ALITOFQADATAMAKERREC_H
