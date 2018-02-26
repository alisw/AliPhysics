#ifndef ALIANALYSISTASKEMCALJETHBASE_H
#define ALIANALYSISTASKEMCALJETHBASE_H

/**
 * @class AliAnalysisTaskEmcalJetHBase
 * @brief Jet-hadron correlations base class
 *
 * Contains funtionality that is shared between the various classes
 *
 * @author Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
 * @date 23 Feb 2018
 */

namespace PWGJE {
namespace EMCALJetTasks {

class AliAnalysisTaskEmcalJetHBase {
 public:
  static double RelativeEPAngle(double jetAngle, double epAngle);
};

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */

#endif /* AliAnalysisTaskEmcalJetHBase.h */
