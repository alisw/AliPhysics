#ifndef AliAnalysisTaskHMTFMC_cxx
#define AliAnalysisTaskHMTFMC_cxx
#include <iostream>
#include <vector>

class TH1F;
class TH1I;
class TGraphErrors;


enum {kHistINEL,kHistNSD,kHistND,kHistSiD,kNHist};
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskHMTFMC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHMTFMC();
  AliAnalysisTaskHMTFMC(const char *name );
  virtual ~AliAnalysisTaskHMTFMC() {}

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

 private:
  TH1F  *fHistEta;           //dNdEta
  TH1F  *fHistNch;           // Multiplicity distribution
  TH1F  *fHistNchUnweighted; // Multiplicity distribution
  TH1F  *fHistRawMult;       // Raw multiplicity distribution
  TH1D  *fHistIev;           // Event counter
  TList *fMyOut;             // Output list

  std::map<int,int> fPrimaryPDGs; // Vector of primary PDGs // WARNING: this does not work on CAF/GRID (not meargiable)
  std::map<int,int> fMotherPDGs; // Vector of Mother's PDGs

  AliAnalysisTaskHMTFMC(const AliAnalysisTaskHMTFMC&); // not implemented
  AliAnalysisTaskHMTFMC& operator=(const AliAnalysisTaskHMTFMC&); // not implemented

  ClassDef(AliAnalysisTaskHMTFMC, 1); // example of analysis
};

#endif
