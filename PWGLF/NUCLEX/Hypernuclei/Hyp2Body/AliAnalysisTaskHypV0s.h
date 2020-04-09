#ifndef ALIAANALYSISTASKHYPV0S_H
#define ALIAANALYSISTASKHYPV0S_H

#include <string>
#include <vector>
#include "AliVertexerHyperTriton2Body.h"
#include "AliAnalysisTaskSE.h"

class AliPIDResponse;


class AliAnalysisTaskHypV0s : public AliAnalysisTaskSE
{
public:

  AliAnalysisTaskHypV0s(std::string name = "HypV0s");
  virtual ~AliAnalysisTaskHypV0s();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  static AliAnalysisTaskHypV0s *AddTask(TString suffix = "");


  AliVertexerHyperTriton2Body fV0Vertexer; //


private:

  AliInputEventHandler *fInputHandler; //!
  AliPIDResponse *fPIDResponse;        //! PID response object
  AliAnalysisTaskHypV0s(
      const AliAnalysisTaskHypV0s &); // not implemented
  AliAnalysisTaskHypV0s &operator=(
      const AliAnalysisTaskHypV0s &); // not implemented

  ClassDef(AliAnalysisTaskHypV0s, 2);
};

#endif