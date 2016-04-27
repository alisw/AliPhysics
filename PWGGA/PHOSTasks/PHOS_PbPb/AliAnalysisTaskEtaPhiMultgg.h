#ifndef AliAnalysisTaskEtaPhiMultgg_cxx
#define AliAnalysisTaskEtaPhiMultgg_cxx

// Class to analyze multiparticle Eta-Phi correlations
// Author: D.Peresunko

class THashList ;
class AliPHOSGeometry;
class AliCaloPhoton ;
class AliAODTrack ;
class AliEPFlattener ;
class AliV0ReaderV1 ;
class AliConvEventCuts ;
class AliConversionPhotonCuts ;
class AliAODConversionPhoton ;
class AliEMCALGeometry ;

#include "AliAnalysisTaskEtaPhigg.h"

class AliAnalysisTaskEtaPhiMultgg : public AliAnalysisTaskEtaPhigg {
public:
    
  
  AliAnalysisTaskEtaPhiMultgg(const char *name = "AliAnalysisTaskEtaPhiMultgg");
  virtual ~AliAnalysisTaskEtaPhiMultgg() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
    
  
private:
  AliAnalysisTaskEtaPhiMultgg(const AliAnalysisTaskEtaPhiMultgg&); // not implemented
  AliAnalysisTaskEtaPhiMultgg& operator=(const AliAnalysisTaskEtaPhiMultgg&); // not implemented
  Bool_t IsSameKtBin(Double_t e1, Double_t e2) ;
  
private:
  
  
  ClassDef(AliAnalysisTaskEtaPhiMultgg, 1); // PHOS analysis task
};

#endif
