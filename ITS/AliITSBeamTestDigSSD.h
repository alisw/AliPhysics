#ifndef ALIITSBEAMTESTDIGSSD_H
#define ALIITSBEAMTESTDIGSSD_H

////////////////////////////////////////////////////
//  Class to define                               //
//  SSD beam test raw 2 dig conv.                 //
//                                                //
//  Origin: Enrico Fragiacomo                     //
//                                                //
////////////////////////////////////////////////////

#include "AliITSBeamTestDig.h"




class AliITSBeamTestDigSSD: public AliITSBeamTestDig {
 
 public:

 
  AliITSBeamTestDigSSD();
  AliITSBeamTestDigSSD(const Text_t* name, const Text_t* title);
  virtual ~AliITSBeamTestDigSSD();

  virtual void Exec(Option_t* opt);
  
  //  void SetRawReaderDate(AliRawReaderDate* rd) {fReaderDate=rd;}
  // void SetHeader(Bool_t H){fFlagHeader = H;}

  //void SetTree(TTree* treedig) {fTreeD=treedig;}
    
 protected:      
  
   

  Bool_t fFlagHeader;  //flag for the header


 ClassDef(AliITSBeamTestDigSSD,1)  // An Alice SSD beam test run

 };


#endif

    

