#ifndef ALIITSBEAMTESTDIGSDD_H
#define ALIITSBEAMTESTDIGSDD_H

////////////////////////////////////////////////////
//  Class to define                               //
//  SDD beam test raw 2 dig conv.                 //
//  Origin: E. Crescio crescio@to.infn.it         //
//                                                //
////////////////////////////////////////////////////

#include "AliITSBeamTestDig.h"
#include "AliITSBeamTest.h"

class AliRawReaderDate;
class AliITSRawStreamSDD;

class AliITSBeamTestDigSDD: public AliITSBeamTestDig {
 
 public:

 
  AliITSBeamTestDigSDD();
  AliITSBeamTestDigSDD(const Text_t* name, const Text_t* title);
  AliITSBeamTestDigSDD(const AliITSBeamTestDigSDD& bt);
  AliITSBeamTestDigSDD& operator=(AliITSBeamTestDigSDD &bt);

  virtual ~AliITSBeamTestDigSDD();

  void SetBtPeriod(BeamtestPeriod_t per=kNov04) {fBtPer=per;}
  void SetThreshold(Int_t threshold) {fThreshold=threshold;}

  BeamtestPeriod_t GetBtPeriod() const {return fBtPer;}
  Int_t GetThreshold() const {return fThreshold;}
  
  virtual void Exec(Option_t* opt);

  
 
 protected:      

  Int_t  fSDDEvType;                 //SDD event type (real, calibration)
  const UInt_t* fSubEventAttributes; //SDD sub-event attributes 
  BeamtestPeriod_t fBtPer;           //beam test version 
                                     // November 2004 = kNov04
                                     // August 2004 = kAug04
  Int_t fThreshold;                  // Low carlos threshold 
  AliITSRawStreamSDD* fStreamer;     //! SDD streamer

 private: 

  Int_t GetEventType();


  ClassDef(AliITSBeamTestDigSDD,1)  // An Alice SDD beam test digitizer

 };



#endif

    
