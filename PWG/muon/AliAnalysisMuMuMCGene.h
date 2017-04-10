#ifndef ALIANALYSISMUMUMCGENE_H
#define ALIANALYSISMUMUMCGENE_H

/**
 * \class AliAnalysisMuMuMCGene
 * \brief Histogramming of basic distributions of the MC generated event.
 * \author L. Aphecetche (Subatech)
 */

#include "AliAnalysisMuMuBase.h"

#include <vector>
#include <string>
#include <set>

class AliGenEventHeader;
class AliVParticle;

class AliAnalysisMuMuMCGene : public AliAnalysisMuMuBase
{
public:
  AliAnalysisMuMuMCGene();
  virtual ~AliAnalysisMuMuMCGene() {}

  void FillHistosForMCEvent(const char* eventSelection, const char* triggerClassName,
                            const char* centrality);

  virtual void DefineHistogramCollection(const char* eventSelection, const char* triggerClassName,
                                         const char* centrality, Bool_t mix =kFALSE);

  Bool_t SelectAnyTriggerClass(const TString& firedTriggerClasses, TString& acceptedTriggerClasses) const;


private:

  AliGenEventHeader* GetGenEventHeader(const AliVEvent& event) const;

  std::vector<std::string> fParticlesOfInterest;
  std::set<int> fPDGCodeOfInterest;

  Bool_t ParticleOfInterest(const AliVParticle& part) const;

  ClassDef(AliAnalysisMuMuMCGene,1) // implementation of AliAnalysisMuMuBase for MCGene event properties
};

#endif
