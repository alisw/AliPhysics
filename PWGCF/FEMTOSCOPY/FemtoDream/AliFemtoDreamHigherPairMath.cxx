/*
 * AliFemtoDreamHigherPairMath.cxx
 *
 *  Created on: Jul 3, 2019
 *      Author: schmollweger
 */

#include <AliFemtoDreamHigherPairMath.h>
#include "TMath.h"
static const float piHi = TMath::Pi();

AliFemtoDreamHigherPairMath::AliFemtoDreamHigherPairMath()
    : fBField(-99.),
      fDeltaPhiEtaMax(-99.){

}

AliFemtoDreamHigherPairMath::~AliFemtoDreamHigherPairMath() {
  // TODO Auto-generated destructor stub
}

bool AliFemtoDreamHigherPairMath::PassesPairSelection(
    AliFemtoDreamBasePart& part1, AliFemtoDreamBasePart& part2, bool Recalculate) {
  bool outBool = true;
  //Method calculates the average separation between two tracks
  //at different radii within the TPC and rejects pairs which a
  //too low separation
  if(Recalculate) {
    RecalculatePhiStar(part1);
    RecalculatePhiStar(part2);
  }
  unsigned int nDaug1 = part1.GetPhiAtRaidius().size();
  unsigned int nDaug2 = part2.GetPhiAtRaidius().size();
  // if nDaug == 1 => Single Track, else decay
  std::vector<float> eta1 = part1.GetEta();
  std::vector<float> eta2 = part2.GetEta();
  for (unsigned int iDaug1 = 0; iDaug1 < nDaug1 && outBool; ++iDaug1) {
    std::vector<float> PhiAtRad1 = part1.GetPhiAtRaidius().at(iDaug1);
    float etaPar1;
    if (nDaug1 == 1) {
      etaPar1 = eta1.at(0);
    } else {
      etaPar1 = eta1.at(iDaug1 + 1);
    }
    for (unsigned int iDaug2 = 0;
        iDaug2 < part2.GetPhiAtRaidius().size() && outBool; ++iDaug2) {
      std::vector<float> phiAtRad2 = part2.GetPhiAtRaidius().at(iDaug2);
      float etaPar2;
      if (nDaug2 == 1) {
        etaPar2 = eta2.at(0);
      } else {
        etaPar2 = eta2.at(iDaug2 + 1);
      }
      float deta = etaPar1 - etaPar2;
      const int size =
          (PhiAtRad1.size() > phiAtRad2.size()) ?
              phiAtRad2.size() : PhiAtRad1.size();
      for (int iRad = 0; iRad < size; ++iRad) {
        float dphi = PhiAtRad1.at(iRad) - phiAtRad2.at(iRad);
        if (dphi > piHi) {
          dphi += -piHi * 2;
        } else if (dphi < -piHi) {
          dphi += piHi * 2;
        }
        dphi = TVector2::Phi_mpi_pi(dphi);
        if (dphi * dphi + deta * deta < fDeltaPhiEtaMax) {
          outBool = false;
          break;
        }
      }
    }
  }
  return outBool;
}

void AliFemtoDreamHigherPairMath::RecalculatePhiStar(
    AliFemtoDreamBasePart &part) {
  //Only use this for Tracks, this was implemented for the case where particles
  //were added a random phi to obtain an uncorrelated sample. This should be extended
  //to the daughter tracks. For some reason the momentum is stored in a TVector3,
  //and not a std::vector<TVector3>. Should be changed in time, conceptually the
  //difficulty is to make this work in a proper way after the mother,
  //e.g. the Lambda was phi shifted.
  //
  if (fBField < -95) {
    AliWarning("BField was most probably not set! PhiStar Calculation meaningless. \n");
  }
  part.GetPhiAtRaidius().clear();
  float TPCradii[9] = { 85., 105., 125., 145., 165., 185., 205., 225., 245. };
  std::vector<float> tmpVec;
  auto phi0 = part.GetMomentum().Phi();
  float pt = part.GetPt();
  float chg = part.GetCharge().at(0);
  for (int radius = 0; radius < 9; radius++) {
    tmpVec.push_back(
        phi0
            - TMath::ASin(
                0.1 * chg * fBField * 0.3 * TPCradii[radius] * 0.01
                    / (2. * pt)));
  }
  part.SetPhiAtRadius(tmpVec);
}
