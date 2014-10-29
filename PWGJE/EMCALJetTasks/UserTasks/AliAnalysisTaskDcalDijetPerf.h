#ifndef ALIANALYSISTASKDCALDIJETPERF_H
#define ALIANALYSISTASKDCALDIJETPERF_H

// $Id$

class TH1;
class TH2;
class TH3;
class THnSparse;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskDcalDijetPerf : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskDcalDijetPerf();
  AliAnalysisTaskDcalDijetPerf(const char *name);
  virtual ~AliAnalysisTaskDcalDijetPerf();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

 protected:
  Float_t                     RelativePhi(Double_t mphi,Double_t vphi) const;
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;
  
  // General histograms
  TH1                       **fHistTracksPt;            //!Track pt spectrum
  TH2                       **fHistTracksEtaPhi;        //!Track eta phi
  TH1                       **fHistClustersPt;          //!Cluster pt spectrum
  TH1                       **fHistLeadingJetPt;        //!Leading jet pt spectrum
  TH2                       **fHistJetsPhiEta;          //!Phi-Eta distribution of jets
  TH2                       **fHistJetsPtArea;          //!Jet pt vs. area
  TH2                       **fHistJetsPtLeadHad;       //!Jet pt vs. leading hadron
  TH2                       **fHistJetsCorrPtArea;      //!Jet pt - bkg vs. area

  THnSparse                  *fHistJet1;                //!jet collection 1
  THnSparse                  *fHistJet1m;               //!jet collection 1 matched
  THnSparse                  *fHistJet1nm;              //!jet collection 1 unmatched
  THnSparse                  *fHistJet2;                //!jet collection 2
  THnSparse                  *fHistJet1to2;             //!jet collection 1 and 2
  THnSparse                  *fHistDiJet1;              //!Dijet collection 1 and 3
  THnSparse                  *fHistDiJet1m;              //!Dijet collection 1 and 3 matched
  
  AliJetContainer            *fJetsCont;                   //!Jets Jet 1
  AliJetContainer            *fJetsCont2;                  //!Jets Trigger Jer
  AliJetContainer            *fJetsCont3;                  //!Jets DiJet
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters  

 private:
  AliAnalysisTaskDcalDijetPerf(const AliAnalysisTaskDcalDijetPerf&);            // not implemented
  AliAnalysisTaskDcalDijetPerf &operator=(const AliAnalysisTaskDcalDijetPerf&); // not implemented

  ClassDef(AliAnalysisTaskDcalDijetPerf, 3)
};
#endif
