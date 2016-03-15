#ifndef ALIANALYSISTASKEMCALITERATORTEST_H
#define ALIANALYSISTASKEMCALITERATORTEST_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcal.h"
#include <TString.h>

class THistManager;

/**
 * \class AliAnalysisTaskEmcalIteratorTest
 * \brief Unit test for the c++ stl iterators in the EMCAL containers
 * \ingroup EMCALFWTASKS
 * \author Markus Fasel <mfasel@lbl.gov>, Lawrence Berkeley National Laboratory
 * \date March 16, 2016
 *
 * This class serves as basic unit test of the stl iterators in the
 * different EMCAL containers AliClusterContainer and AliParticleContainers.
 * For the particle containers it tests the both the track and the MC particle
 * container, if available.
 *
 * As the unit tests usually return a number based on the outcome of the tests,
 * this number is monitored in a histogram corresponding to the iterator which
 * is tested.
 *
 * In order to run the test, the user needs to prepare an EMCAL train with
 * wagons for the cluster maker, track selection and MC particle maker. Tests
 * are added in the following way:
 *
 * ~~~{.cxx}
 * AliAnalysisTaskEmcalIteratorTest *test = new AliAnalysisTaskEmcalIteratorTest("iteratorTest");
 * // Adding cluster container
 * test->AddClusterContainer("contname");       // contname needs to match the name of the container to be tested
 * test->SetClusterContainerName("contname");
 * // Adding MC particle container
 * test->AddMCParticleContainer("mcontname");
 * test->SetMCParticleContainerName("mccontname");
 * // Adding Track container
 * test->AddTrackContainer("trackcontname");
 * test->SetTrackContainerName("trackcontname");
 * ~~~
 */
class AliAnalysisTaskEmcalIteratorTest : public AliAnalysisTaskEmcal {
public:
  AliAnalysisTaskEmcalIteratorTest();
  AliAnalysisTaskEmcalIteratorTest(const char *name);
  virtual ~AliAnalysisTaskEmcalIteratorTest();

  void SetClusterContainerName(TString name)        { fNameClusterContainer = name; }
  void SetTrackContainerName(TString name)          { fNameTrackContainer = name; }
  void SetMCParticleContainerName(TString name)     { fNameMCParticleContainer = name; }

protected:

  virtual void UserCreateOutputObjects();
  virtual bool Run();

  THistManager                *fHistos;                     //!<!  Histogram manager
  TString                     fNameClusterContainer;        /// Name of the cluster container
  TString                     fNameMCParticleContainer;     /// Name of the MC particle container
  TString                     fNameTrackContainer;          /// Name of the track container

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalIteratorTest, 1);
  /// \endcond
};

#endif /* ALIANALYSISTASKEMCALITERATORTEST_H */
