#ifndef ALIITSBEAMTESTDIG_H
#define ALIITSBEAMTESTDIG_H

////////////////////////////////////////////////////
//  Class to define                               //
//  beam test raw 2 dig conv.                     //
//  Origin: E. Crescio crescio@to.infn.it         //
//  November 2004                                 //
////////////////////////////////////////////////////

#include "TTask.h"
#include "AliITSEventHeader.h"

class AliITSBeamTest;
class AliRawReaderDate;
class TTree;


class AliITSBeamTestDig: public TTask {
 
 public:

 
  AliITSBeamTestDig();
  AliITSBeamTestDig(const Text_t* name, const Text_t* title);
  AliITSBeamTestDig(const AliITSBeamTestDig& bt);
  AliITSBeamTestDig& operator=(AliITSBeamTestDig &bt);
  virtual ~AliITSBeamTestDig() {}
 
  void SetRawReaderDate(AliRawReaderDate* rd) {fReaderDate=rd;}
  void SetTree(TTree* treedig) {fTreeD=treedig;}
  void SetITSEventHeader(AliITSEventHeader* header){fITSHeader = header;}

  void SetBeamTest(AliITSBeamTest* bt) {fBt=bt;}


 protected:      
  
  AliITSEventHeader* fITSHeader;     // its event header
  AliRawReaderDate* fReaderDate;     // !reader date;
  TTree* fTreeD;                     // tree of digits
 
  AliITSBeamTest* fBt;               // !beam test object

  ClassDef(AliITSBeamTestDig,1)   // its beam test digitization 

 };



#endif

    
