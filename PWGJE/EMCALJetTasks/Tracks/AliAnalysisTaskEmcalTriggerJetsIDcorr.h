#ifndef ALIANALYSISTASKEMCALTRIGGERJETSIDCORR_H
#define ALIANALYSISTASKEMCALTRIGGERJETSIDCORR_H
// Copyright (C) 2017, Copyright Holders of the ALICE Collaboration
// All rights reserved.

#include <AliAnalysisTaskEmcalJet.h>
#include "AliPID.h"
#include <exception>
#include <vector>

class AliJetContainer;
class AliPIDResponse;
class THistManager;
class TString;
class TVector3;

namespace PWGJE {
  
namespace EMCALJetTasks {

struct CorrParticleInfo {
  Double_t fPt;
  Double_t fDR;
  Double_t fMass;
};

class AliAnalysisTaskEmcalTriggerJetsIDcorr: public AliAnalysisTaskEmcalJet {
public:

  class TOFMassException : public std::exception {
  public:
    TOFMassException() : std::exception() {}
    virtual ~TOFMassException() throw() {}
    virtual const char *what() const throw() { return "TOF mass cannot be calculated for particle"; }
  };

  class TPCdEdxException : public std::exception {
  public:
    TPCdEdxException(AliPID::EParticleType type): std::exception(), fParticle(type), fMessage(TString::Format("TPC dE/dx information not available for particle type %s", AliPID::ParticleName(type))) {}
    virtual ~TPCdEdxException() throw() {}
    virtual const char *what() const throw() { return fMessage.Data(); }
    AliPID::EParticleType GetParticleType() const throw() { return fParticle; }

  private:
    AliPID::EParticleType     fParticle;
    TString                   fMessage;
  };

  AliAnalysisTaskEmcalTriggerJetsIDcorr();
  AliAnalysisTaskEmcalTriggerJetsIDcorr(const char *name);
  virtual ~AliAnalysisTaskEmcalTriggerJetsIDcorr();

  static AliAnalysisTaskEmcalTriggerJetsIDcorr *AddTaskEmcalTriggerJetsIDcorr(const char *name);

protected:

  virtual void UserCreateOutputObjects();
  virtual void UserExecOnce();
  virtual bool Run();

  double GetTOFMass(const AliVTrack *const track) const;
  std::vector<AliVTrack *> GetTPCPIDCandidates(AliPID::EParticleType type) const;
  std::vector<CorrParticleInfo> CorrelateCandidatesToJet(const TVector3 &jet, std::vector<AliVTrack *> candidates) const;

private:
  AliAnalysisTaskEmcalTriggerJetsIDcorr(const AliAnalysisTaskEmcalTriggerJetsIDcorr &);
  AliAnalysisTaskEmcalTriggerJetsIDcorr &operator=(const AliAnalysisTaskEmcalTriggerJetsIDcorr &);

  AliJetContainer       *fJetCont;            //!<! Jet container
  AliPIDResponse        *fPIDResponse;        //!<! PID Response handler
  THistManager          *fHistos;             //!<! Histogram handler

  ClassDef(AliAnalysisTaskEmcalTriggerJetsIDcorr, 1)
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKEMCALTRIGGERJETSIDCORR_H */
