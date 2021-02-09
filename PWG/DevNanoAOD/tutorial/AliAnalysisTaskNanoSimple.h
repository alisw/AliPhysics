#ifndef AliAnalysisTaskNanoSimple_H
#define AliAnalysisTaskNanoSimple_H

#include "AliAnalysisTaskSE.h"

class  AliAnalysisTaskNanoSimple : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskNanoSimple(const char* name="AliAnalysisTaskNanoSimple");
  virtual           ~AliAnalysisTaskNanoSimple();

  // Implementation of interace methods
  virtual     void   UserCreateOutputObjects();
  virtual     void   UserExec(Option_t *option);

private:
  AliAnalysisTaskNanoSimple(const  AliAnalysisTaskNanoSimple &det);
  AliAnalysisTaskNanoSimple&   operator=(const  AliAnalysisTaskNanoSimple &det);

  // Histogram settings
  TList*              fListOfHistos;    //  Output list of containers

  ClassDef(AliAnalysisTaskNanoSimple, 1); // Analysis task for correlation development
};

#endif


