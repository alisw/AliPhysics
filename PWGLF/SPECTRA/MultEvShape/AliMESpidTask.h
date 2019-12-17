#ifndef ALIMESPIDTASK_H
#define ALIMESPIDTASK_H

/////////////////////////////////////////////////////////////////////////////
//  PID task for Multiplicity and Event Shape group                       //
//  Authors:                                                              //
//    Cristi Andrei <Cristian.Andrei@cern.ch>                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIMESBASETASK_H
#include "AliMESbaseTask.h"
#endif

class AliMESpidTask : public AliMESbaseTask
{
public:
  AliMESpidTask();
  AliMESpidTask(const char *name);
  virtual ~AliMESpidTask();

  Bool_t         BuildQAHistos();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *opt);
  virtual Bool_t PostProcess();
  Double_t       ComputeDeltaPhi(Double_t, Double_t);

private:
  AliMESpidTask(const AliMESpidTask&);
  AliMESpidTask& operator=(const AliMESpidTask&);

  enum outputs_t {slot_MultEst, slot_Response, slot_Miss, slot_AllESD, slot_Gen, slot_NoEvts, slot_PIDQA, slot_DeltaPhi, slot_testTree, slot_DCA, slot_Fakes};

  ClassDef(AliMESpidTask, 4)            // PID task for the Multi Event Shape
};

#endif
