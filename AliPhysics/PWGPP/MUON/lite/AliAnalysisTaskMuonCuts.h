#ifndef ALIANALYSISTASKMUONCUTS_H
#define ALIANALYSISTASKMUONCUTS_H

/* $Id: AliAnalysisTaskMuonCuts.h 47782 2011-02-24 18:37:31Z martinez $ */ 

// Class for muon pxDCA cuts tuning
// 
//  Author Diego Stocco

#include "AliVAnalysisMuon.h"
#include "TArrayD.h"

class TObjArray;
class AliMergeableCollection;
class TString;
class TAxis;
class AliVParticle;
class AliAODEvent;

class AliAnalysisTaskMuonCuts : public AliVAnalysisMuon {
 public:
  AliAnalysisTaskMuonCuts();
  AliAnalysisTaskMuonCuts(const char *name, const AliMuonTrackCuts& cuts);
  virtual ~AliAnalysisTaskMuonCuts();

  virtual void   Terminate(Option_t *option);

  void MyUserCreateOutputObjects();
  void ProcessEvent(TString physSel, const TObjArray& selectTrigClasses, TString centrality);

  enum {
    kThetaAbs23,
    kThetaAbs310,
    kNthetaAbs
  };

  void SetSigmaCuts(Int_t nSigmaCuts = -1, Double_t* sigmaCuts = 0x0);

 private:

  AliAnalysisTaskMuonCuts(const AliAnalysisTaskMuonCuts&);
  AliAnalysisTaskMuonCuts& operator=(const AliAnalysisTaskMuonCuts&);

  // Histograms to extract average DCA position
  enum {
    kDCAxVsP,      ///< DCA_x vs momentum
    kDCAyVsP,      ///< DCA_y vs momentum
    kPdcaVsP,      ///< p x DCA vs momentum (binning for fit)
    kPDCAVsPCheck, ///< p x DCA vs momentum (check beam gas)
    kDCAVsPCheck,  ///< DCA vs momentum
    kChiProbVsP,   ///< Chi square probability vs momentum
    kSigmaVsPt,    ///< pt distribution for different p x DCA sigma cuts
    kSigmaVsEta,   ///< eta distribution for different p x DCA sigma cuts
    kNhistoTypes   ///< Number of histograms
  };

  TString GetHistoName(Int_t histoTypeIndex, Int_t thetaAbsIndex, Int_t srcIndex);

  TObjArray* fHistoTypeKeys;   ///< Base histogram name
  TObjArray* fThetaAbsKeys;    ///< Name of theta at absorber end
  TArrayD fSigmaCuts;          ///< List of sigma cuts

  ClassDef(AliAnalysisTaskMuonCuts, 1); // Single muon analysis
};

#endif
