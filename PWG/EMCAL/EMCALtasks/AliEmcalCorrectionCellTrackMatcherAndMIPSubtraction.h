#ifndef ALIEMCALCORRECTIONCELLTRACKMATCHERANDMIPSUBTRACTION_H
#define ALIEMCALCORRECTIONCELLTRACKMATCHERANDMIPSUBTRACTION_H

#include "AliEmcalCorrectionComponent.h"

//#if !(defined(__CINT__) || defined(__MAKECINT__))
//#include "AliEmcalContainerIndexMap.h"
//#endif

class TH1;
class TClonesArray;
class AliVParticle;

/**
 * @class AliEmcalCorrectionCellTrackMatching
 * @ingroup EMCALCORRECTIONFW
 * @brief Track matching of charged tracks to cells, subtracting MIP, prior to clusterization.
 *
 * Track matching of charged tracks to cells. Once matched it will subtract the MIP energy from the cell.
 * In case that Ecell<EmipData(MC), then energy of the cell is set to zero.
 * This all is done prior to clusterization.
 *
 * @author Mike Sas, mike.sas@cern.ch, Yale
 * @date Nov 10 2020
 */

class AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction();
  virtual ~AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  void ExecOnce();
  Bool_t Run();

protected:

  Bool_t        IsTrackInEmcalAcceptance(AliVTrack* part, Double_t edges=0.9) const;

  Double_t               fEmipData;            ///< Energy of MIP used to subtract from cell, in case of data
  Double_t               fEmipMC;              ///< Energy of MIP used to subtract from cell, in case of MC

  Bool_t                 fDoPropagation;      ///< Bool to switch propagation on/off

  //#if !(defined(__CINT__) || defined(__MAKECINT__))
  // Handle mapping between index and containers
  //AliEmcalContainerIndexMap <AliParticleContainer, AliVParticle> fParticleContainerIndexMap; //!<! Mapping between index and particle containers
  //#endif

  TH2          *fCellTrackMatchdEtadPhi;       //!<!dEtadPhi distribution
  TH1          *fCellTrackMatchdEtaDiff;       //!<! Difference of estimated track position (eta) on EMCal before and after propagation
  TH1          *fCellTrackMatchdPhiDiff;       //!<! Difference of estimated track position (phi) on EMCal before and after propagation
  TH2          *fCellNTrackMatch;              //!<!NTrackMatch distribution
  TH1          *fCellTrackMatchEbefore;        //!<!Ebefore distribution
  TH1          *fCellTrackMatchEafter;         //!<!Eafter distribution

 private:
  AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction(const AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction &);               // Not implemented
  AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction &operator=(const AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction &);    // Not implemented

  // Allows the registration of the class so that it is available to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionCellTrackMatcherAndMIPSubtraction, 3); // EMCal cell track matcher
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCELLTRACKMATCHERANDMIPSUBTRACTION_H */
