#ifndef ALIGLOBALQADATAMAKER_H
#define ALIGLOBALQADATAMAKER_H

/*
 The class for calculating the global (not detector specific) quality assurance.
 It reuses the following TLists from its base class 
    AliQADataMaker::fRecPointsQAList (for keeping the track residuals)
    AliQADataMaker::fESDsQAList      (for keeping global ESD QA data)
*/

#include "AliQADataMakerRec.h"

class AliESDEvent;

class AliGlobalQADataMaker: public AliQADataMakerRec {
public:
  enum {
    kClr0,kClr1,kClr2,
    kTrk0,kTrk1,kTrk2,kTrk3,kTrk4,kTrk5,kTrk6,
    kV0s0,kV0s1,
    kPid0,kPid1,kPid2,
    kLast
  };
  AliGlobalQADataMaker(const Char_t *name="Global", 
                       const Char_t *title="Global QA data maker"):
	AliQADataMakerRec(name,title) {;}
  AliGlobalQADataMaker(const AliQADataMakerRec& qadm):
	AliQADataMakerRec(qadm) {;}

private:
	void   EndOfDetectorCycle(AliQA::TASKINDEX_t, TObjArray ** list) ;

	void InitRaws(); 
	void InitRecPoints();
  void InitESDs();

	void MakeRaws(AliRawReader* rawReader) ; 
  void MakeESDs(AliESDEvent *event);

  void StartOfDetectorCycle() {;}

  AliGlobalQADataMaker &operator=(const AliGlobalQADataMaker &qadm);

  ClassDef(AliGlobalQADataMaker,1)  // Global QA 
};

#endif
