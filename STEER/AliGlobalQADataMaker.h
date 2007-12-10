#ifndef ALIGLOBALQADATAMAKER_H
#define ALIGLOBALQADATAMAKER_H

/*
 The class for calculating the global (not detector specific) quality assurance.
 It reuses the following TLists from its base class 
    AliQADataMaker::fRecPointsQAList (for keeping the track residuals)
    AliQADataMaker::fESDsQAList      (for keeping global ESD QA data)
*/

#include "AliQADataMaker.h"

class AliESDEvent;

class AliGlobalQADataMaker: public AliQADataMaker {
public:
  AliGlobalQADataMaker(const Char_t *name="Global", 
                       const Char_t *title="Global QA data maker"):
    AliQADataMaker(name,title) {;}
  AliGlobalQADataMaker(const AliQADataMaker& qadm):
    AliQADataMaker(qadm) {;}

  void InitRecPoints();

  void StartOfDetectorCycle() {;}

private:   
  AliGlobalQADataMaker &operator=(const AliGlobalQADataMaker &qadm);

  ClassDef(AliGlobalQADataMaker,1)  // Global QA 
};

#endif
