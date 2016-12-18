#ifndef ALIANALYSISTASKEMCALNEUTRALJETS_H_
#define ALIANALYSISTASKEMCALNEUTRALJETS_H_

#include "AliAnalysisTaskEmcalJet.h"
#include <TString.h>

class THistManager;

class AliClusterContainer;
class AliJetContainer;

class AliAnalysisTaskEmcalNeutralJets : public AliAnalysisTaskEmcalJet {
public:
  AliAnalysisTaskEmcalNeutralJets();
  AliAnalysisTaskEmcalNeutralJets(const char *name);
  virtual ~AliAnalysisTaskEmcalNeutralJets();

  void SetRequestTrigger(ULong_t triggers, const TString &triggerstring){
    fTriggerBits = triggers;
    fTriggerString = triggerstring;
  }

  void SetContainers(const TString &name02, const TString &name04, const TString &nameClusters){
    fNameR02jets = name02;
    fNameR04jets = name04;
    fNameClusters = nameClusters;
  }

protected:
  virtual void UserCreateOutputObjects();
  virtual bool IsEventSelected();
  virtual bool Run();
  virtual bool FillHistograms();

private:
  THistManager                              *fHistos;                   //!<! Histogram handler
  UInt_t                                     fTriggerBits;              ///< Trigger Bits
  TString                                    fTriggerString;            ///< Trigger string (distinguish EGA triggers)

  TString                                    fNameR02jets;              ///< Name of the jet container for R=0.2 jets
  TString                                    fNameR04jets;              ///< Name of the jet container for R=0.4 jets
  TString                                    fNameClusters;             ///< Name of the cluster container
  AliJetContainer                           *fR02jets;                  //!<! Link to jet container for R=0.2 jets
  AliJetContainer                           *fR04jets;                  //!<! Link to jet container for R=0.4 jets
  AliClusterContainer                       *fClusters;                 //!<! Link to cluster container;

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalNeutralJets, 1);
  /// \endcond
};

#endif /* ALIANALYSISTASKEMCALNEUTRALJETS_H_ */
