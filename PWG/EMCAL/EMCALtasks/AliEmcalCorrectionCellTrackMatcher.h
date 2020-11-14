#ifndef ALIEMCALCORRECTIONCELLTRACKMATCHER_H
#define ALIEMCALCORRECTIONCELLTRACKMATCHER_H

#include "AliEmcalCorrectionComponent.h"

#if !(defined(__CINT__) || defined(__MAKECINT__))
#include "AliEmcalContainerIndexMap.h"
#endif

class TH1;
class TClonesArray;
class AliVParticle;

/**
 * @class AliEmcalCorrectionCellTrackMatching
 * @ingroup EMCALCORRECTIONFW
 * @brief Track matching of charged tracks to cells, subtracting MIP, prior to clusterization.
 *
 * Track matching of charged tracks to cells, subtracting MIP, prior to clusterization.
 *
 * @author Mike Sas, mike.sas@cern.ch, Yale
 * @date Nov 10 2020
 */

class AliEmcalCorrectionCellTrackMatcher : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellTrackMatcher();
  virtual ~AliEmcalCorrectionCellTrackMatcher();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  void ExecOnce();
  Bool_t Run();
  
protected:
  
  Bool_t        IsTrackInEmcalAcceptance(AliVParticle* part, Double_t edges=0.9) const;

  Double_t               fMipE;                           ///< Energy of MIP used to subtract from cell

  #if !(defined(__CINT__) || defined(__MAKECINT__))
  // Handle mapping between index and containers
  AliEmcalContainerIndexMap <AliParticleContainer, AliVParticle> fParticleContainerIndexMap; //!<! Mapping between index and particle containers
  #endif

  TClonesArray *fEmcalTracks;                  //!<!emcal tracks
  TH2          *fCellTrackMatchdEtadPhi;       //!<!dEtadPhi distribution
  TH2          *fCellNTrackMatch;              //!<!NTrackMatch distribution
  TH1          *fCellTrackMatchEbefore;        //!<!Ebefore distribution
  TH1          *fCellTrackMatchEafter;         //!<!Eafter distribution

 private:
  AliEmcalCorrectionCellTrackMatcher(const AliEmcalCorrectionCellTrackMatcher &);               // Not implemented
  AliEmcalCorrectionCellTrackMatcher &operator=(const AliEmcalCorrectionCellTrackMatcher &);    // Not implemented

  // Allows the registration of the class so that it is available to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellTrackMatcher> reg;
  
  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionCellTrackMatcher, 1); // EMCal cell energy variation component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCELLTRACKMATCHER_H */
