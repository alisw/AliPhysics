/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <THistManager.h>

#include "AliClusterContainer.h"
#include "AliMCParticleContainer.h"
#include "AliTrackContainer.h"

#include "AliAnalysisTaskEmcalIteratorTest.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalIteratorTest)
/// \endcond

/**
 * Default (I/O) constructor
 */
AliAnalysisTaskEmcalIteratorTest::AliAnalysisTaskEmcalIteratorTest():
  AliAnalysisTaskEmcal(),
  fHistos(NULL),
  fNameClusterContainer(""),
  fNameMCParticleContainer(""),
  fNameTrackContainer("")
{

}

/**
 * Named constructor, initializing the histograms from the AliAnalysisTaskEmcal.
 */
AliAnalysisTaskEmcalIteratorTest::AliAnalysisTaskEmcalIteratorTest(const char *name):
  AliAnalysisTaskEmcal(name, true),
  fHistos(NULL),
  fNameClusterContainer(""),
  fNameMCParticleContainer(""),
  fNameTrackContainer("")
{

}

/**
 * Destructor, deleting histogram container
 */
AliAnalysisTaskEmcalIteratorTest::~AliAnalysisTaskEmcalIteratorTest(){
}

/**
 * Creating histograms monitoring the test results of the unit
 * test for the iterators of the different EMCAL containers.
 */
void AliAnalysisTaskEmcalIteratorTest::UserCreateOutputObjects() {
  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  fHistos = new THistManager("testhistos");
  fHistos->ReleaseOwner();

  const int kNcontNames = 3, kNiterNames = 2;
  TString containernames[kNcontNames] = {"ClusterContainer", "MCParticleContainer", "TrackContainer"},
          iternames[kNiterNames] = {"Accept", "All"};
  for(int icont = 0; icont < kNcontNames; icont++)
    for(int iiter = 0; iiter < kNiterNames; iiter++)
      fHistos->CreateTH1(
          Form("hTest%sIter%s", containernames[icont].Data(), iternames[iiter].Data()),
          Form("Test results for container %s and iterator %s", containernames[icont].Data(), iternames[iiter].Data()),
          3, -0.5, 2.5);

  for(TIter histiter = TIter(fHistos->GetListOfHistograms()).Begin(); histiter != TIter::End(); ++histiter){
    fOutput->Add(*histiter);
  }
  PostData(1, fOutput);
}

/**
 * Running the unit test suites for the three kinds of containers
 * - Cluster container
 * - MC Particle container
 * - Track container
 * and both iterators
 * - all iterator
 * - accept iterator
 * separately. The result of each test is filled in a histogram
 * as the value returned by the test suite (0 - passed, 1 -
 * missing entries, 2 - excess entries). The test is passed
 * in case 100% of the entries of each test category are at 0.
 *
 * @return Always true
 */
bool AliAnalysisTaskEmcalIteratorTest::Run() {
  AliClusterContainer *clustercont = this->GetClusterContainer(fNameClusterContainer.Data());
  int testresult = -1;
  if(clustercont){
    // Run test suite of the cluster container
    testresult = TestClusterContainerIterator(clustercont, 0);
    fHistos->FillTH1("hTestClusterContainerIterAccept", testresult);
    testresult = TestClusterContainerIterator(clustercont, 1);
    fHistos->FillTH1("hTestClusterContainerIterAll", testresult);
  }

  AliMCParticleContainer *mcpcont = this->GetMCParticleContainer(fNameMCParticleContainer.Data());
  if(mcpcont){
    // Run test suite of the particle container
    testresult = TestParticleContainerIterator(mcpcont, 0);
    fHistos->FillTH1("hTestMCParticleContainerIterAccept", testresult);
    testresult = TestParticleContainerIterator(mcpcont, 1);
    fHistos->FillTH1("hTestMCParticleContainerIterAll", testresult);
  }

  AliTrackContainer *trackcont = this->GetTrackContainer(fNameTrackContainer.Data());
  if(mcpcont){
    // Run test suite of the particle container
    testresult = TestParticleContainerIterator(trackcont, 0);
    fHistos->FillTH1("hTestTrackContainerIterAccept", testresult);
    testresult = TestParticleContainerIterator(trackcont, 1);
    fHistos->FillTH1("hTestTrackContainerIterAll", testresult);
  }

  return kTRUE;
}
