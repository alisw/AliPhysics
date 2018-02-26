//
// Base class for Jet-Hadron correlation analysis
//

#include "AliAnalysisTaskEmcalJetHBase.h"

#include <TMath.h>

#include <AliLog.h>

namespace PWGJE {
namespace EMCALJetTasks {

/**
 * Function to calculate angle between jet and EP in the 1st quadrant (0,Pi/2).
 * Adapted from AliAnalysisTaskEmcalJetHadEPpid.
 *
 * @param jetAngle Phi angle of the jet (could be any particle)
 * @param epAngle Event plane angle
 *
 * @return Angle between jet and EP in the 1st quadrant (0,Pi/2)
 */
double AliAnalysisTaskEmcalJetHBase::RelativeEPAngle(double jetAngle, double epAngle)
{
  double dphi = (epAngle - jetAngle);

  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi<-1*TMath::Pi() ) {
    dphi = dphi + 1*TMath::Pi();
  } // this assumes we are doing full jets currently

  if( (dphi>0) && (dphi<1*TMath::Pi()/2) ) {
    // Do nothing! we are in quadrant 1
  } else if ( (dphi>1*TMath::Pi()/2) && (dphi<1*TMath::Pi()) ) {
    dphi = 1*TMath::Pi() - dphi;
  } else if ( (dphi<0) && (dphi>-1*TMath::Pi()/2) ) {
    dphi = std::abs(dphi);
  } else if ( (dphi<-1*TMath::Pi()/2) && (dphi>-1*TMath::Pi()) ) {
    dphi = dphi + 1*TMath::Pi();
  }

  // Warn if we are not in the proper range
  if ( dphi < 0 || dphi > TMath::Pi()/2 ) {
    AliWarningGeneralStream("AliAnalysisTaskEmcalJetHBase") << ": dPHI not in range [0, 0.5*Pi]!\n";
  }

  return dphi;   // dphi in [0, Pi/2]
}

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */
