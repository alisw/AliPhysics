#ifndef ALIITSBEAMTEST_H
#define ALIITSBEAMTEST_H

////////////////////////////////////////////////////
//  Base class to define                          //
//  ITS beam test                                 //
//                                                //   
//  Origin: E. Crescio crescio@to.infn.it         //
////////////////////////////////////////////////////

#include "AliITS.h"

typedef enum {kAug04,kNov04} BeamtestPeriod_t;

class AliITSBeamTest : public AliITS {
 
 public:


  AliITSBeamTest();
  AliITSBeamTest(const char* name,const char *title);
  virtual ~AliITSBeamTest();

  virtual void SetNumberOfSPD(Int_t nSPD) {fNspd=nSPD;}
  virtual void SetNumberOfSDD(Int_t nSDD) {fNsdd=nSDD;}
  virtual void SetNumberOfSSD(Int_t nSSD) {fNssd=nSSD;}

  Int_t GetNSPD() const {return fNspd;}
  Int_t GetNSDD() const {return fNsdd;}
  Int_t GetNSSD() const {return fNssd;}

  Int_t GetNumberOfSubDet(const TString& det) const;
    
  
 protected:

  static const Int_t fgkNumberOfSPD; //number of SPD detectors
  static const Int_t fgkNumberOfSDD; //number of SDD detectors
  static const Int_t fgkNumberOfSSD; //number of SSD detectors


  Int_t     fNspd;                    //Number of SPD modules
  Int_t     fNsdd;                    //Number of SDD modules
  Int_t     fNssd;                    //Number of SSD modules
  

  ClassDef(AliITSBeamTest,0)  // An Alice ITS beam test 

 };


#endif

    
