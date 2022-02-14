#ifndef AliAnalysisTaskNanoValidator_H
#define AliAnalysisTaskNanoValidator_H

#include "AliAnalysisTaskSE.h"

class  AliAnalysisTaskNanoValidator : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskNanoValidator(const char* name="AliAnalysisTaskNanoValidator");
  virtual           ~AliAnalysisTaskNanoValidator();

  // Implementation of interace methods
  virtual     void   UserCreateOutputObjects();
  virtual     void   UserExec(Option_t *option);

private:
  AliAnalysisTaskNanoValidator(const  AliAnalysisTaskNanoValidator &det);
  AliAnalysisTaskNanoValidator&   operator=(const  AliAnalysisTaskNanoValidator &det);

  // Histogram settings
  TList*              fListOfHistos;    //  Output list of containers

  ClassDef(AliAnalysisTaskNanoValidator, 1); // Analysis task for correlation development
};

#endif


