#ifndef ALIANALYSISTASKEMCALJETHUTILS_H
#define ALIANALYSISTASKEMCALJETHUTILS_H

/**
 * @class AliAnalysisTaskEmcalJetHUtils
 * @brief Jet-hadron correlations utilities class
 *
 * Contains funtionality that is shared between the various classes
 *
 * @author Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
 * @date 23 Feb 2018
 */

#include <string>
#include <map>

class AliEmcalJet;

namespace PWGJE {
namespace EMCALJetTasks {

class AliAnalysisTaskEmcalJetHUtils {
 public:
  /**
   * @enum ELeadingHadronBiasType_t
   * @brief Determine the jet leading hadron bias type
   */
  enum ELeadingHadronBiasType_t {
    kCharged = 0,   //!<! Charged hadron bias
    kNeutral = 1,   //!<! Neutral hadron bias
    kBoth = 2       //!<! Leading is from the leader of both
  };
  static const std::map<std::string, ELeadingHadronBiasType_t> fgkLeadingHadronBiasMap; //!<! Map from name to leading hadron bias used with the YAML config
  static double GetLeadingHadronPt(AliEmcalJet * jet, ELeadingHadronBiasType_t leadingHadronType);

  static double RelativeEPAngle(double jetAngle, double epAngle);
};

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */

#endif /* AliAnalysisTaskEmcalJetHUtils.h */
