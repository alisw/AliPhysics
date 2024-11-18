#ifndef ALIAODRELABELINTERFACE_H
#define ALIAODRELABELINTERFACE_H

#include "AliAnalysisTaskSE.h"

// Interface for AOD relabeling that can happen in the V0Reader and RUN3/AliAnalysisTaskAO2Dconverter
class AliAODRelabelInterface : public AliAnalysisTaskSE {
public:
  AliAODRelabelInterface() = default;
  AliAODRelabelInterface(const char* name) : AliAnalysisTaskSE(name) {}
  virtual ~AliAODRelabelInterface() = default;

  // Pure virtual functions for required methods
  virtual bool AreAODsRelabeled() const = 0;
  virtual int IsReaderPerformingRelabeling() const = 0;

  ClassDef(AliAODRelabelInterface, 1);
};

#endif // ALIAODRELABELINTERFACE_H
