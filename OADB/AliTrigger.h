#ifndef AliTrigger_h
#define AliTrigger_h 1

#include "AliBits.h"
#include "TString.h"
namespace AliTrigger{
  static const Int_t NTRIGGERBITS = 55;
  static const AliBits kMB                 = AliBits( 0); // Minimum bias [INT1 = (V0A | V0C | SPD) in pp] or [PBI = (V0A & V0C & ZDCtime) in PbPb2011]
  static const AliBits kINT7               = AliBits( 1); // Minimum bias INT7 = V0A & V0C
  static const AliBits kMUON               = AliBits( 2); // low-pt single muon, INT1 suite
  static const AliBits kHighMult           = AliBits( 3); // High-multiplicity trigger (threshold defined online)
  static const AliBits kEMC1               = AliBits( 4); // EMCAL L0 trigger, INT1 suite
  static const AliBits kCINT5              = AliBits( 5); // Minimum bias INT5 = V0A | V0C
  static const AliBits kCMUS5              = AliBits( 6); // low-pt single muon, INT5 suite
  static const AliBits kMUSPB              = AliBits( 6); // low-pt single muon, PBI suite
  static const AliBits kMUSH7              = AliBits( 7); // high-pt single muon, INT7 suite
  static const AliBits kMUSHPB             = AliBits( 7); // high-pt single muon, PBI suite
  static const AliBits kMUL7               = AliBits( 8); // low-pt like sign dimuon, INT7 suite
  static const AliBits kMuonLikePB         = AliBits( 8); // low-pt like sign dimuon, PBI suite
  static const AliBits kMUU7               = AliBits( 9); // low-pt unlike sign dimuon, INT7 suite
  static const AliBits kMuonUnlikePB       = AliBits( 9); // low-pt unlike sign dimuon, PBI suite
  static const AliBits kEMC7               = AliBits(10); // EMCAL L0 trigger, INT7 suite
  static const AliBits kMUS7               = AliBits(11); // low-pt single muon, INT7 suite
  static const AliBits kPHI1               = AliBits(12); // PHOS L0 trigger, INT1 suite
  static const AliBits kPHI7               = AliBits(13); // PHOS L0 trigger, INT7 suite
  static const AliBits kPHOSPb             = AliBits(13); // PHOS L0 trigger, PBI suite
  static const AliBits kEMCEJE             = AliBits(14); // EMCAL L1 jet trigger
  static const AliBits kEMCEGA             = AliBits(15); // EMCAL L1 gamma trigger
  static const AliBits kCentral            = AliBits(16); // Central PbPb trigger
  static const AliBits kSemiCentral        = AliBits(17); // Semicentral PbPb trigger
  static const AliBits kDG5                = AliBits(18); // Double gap diffractive
  static const AliBits kZED                = AliBits(19); // ZDC electromagnetic dissociation
  static const AliBits kSPI7               = AliBits(20); // Power interaction trigger, INT7 suite
  static const AliBits kSPI                = AliBits(20); // Power interaction trigger
  static const AliBits kINT8               = AliBits(21); // INT8 = 0TVX = T0-vertex requirement
  static const AliBits kMuonSingleLowPt8   = AliBits(22); // low-pt single muon, INT8 suite
  static const AliBits kMuonSingleHighPt8  = AliBits(23); // high-pt single muon, INT8 suite
  static const AliBits kMuonLikeLowPt8     = AliBits(24); // low-pt like sign dimuon, INT8 suite
  static const AliBits kMuonUnlikeLowPt8   = AliBits(25); // low-pt unlike sign dimuon, INT8 suite
  static const AliBits kMuonUnlikeLowPt0   = AliBits(26); // low-pt unlike sign dimuon, no additional requirement
  static const AliBits kUserDefined        = AliBits(27); // Set when custom trigger classes are set in AliPhysicsSelection
  static const AliBits kTRD                = AliBits(28); // Mixture of TRD triggers
  static const AliBits kUserDefined2       = AliBits(29); // Set when custom trigger classes are set in AliPhysicsSelection
  static const AliBits kFastOnly           = AliBits(30); // The fast cluster fired. Set in addition to another trigger bit
  static const AliBits kUserDefined3       = AliBits(31); // Set when custom trigger classes are set in AliPhysicsSelection
  static const AliBits kAny                = AliBits(NTRIGGERBITS,TString("set")); // Accept any trigger
  static const AliBits kAnyINT             = AliTrigger::kMB | AliTrigger::kINT7 | AliTrigger::kCINT5 | AliTrigger::kINT8 | AliTrigger::kSPI7; // Any interaction (minimum bias) trigger
  static const AliBits kPHI8               = AliBits(32); // PHOS L0 trigger, INT8 suite
  static const AliBits kEmcalL1GammaHigh7  = AliBits(33); // EMCAL L1 gamma trigger, high threshold, INT7 suite
  static const AliBits kEmcalL1GammaLow7   = AliBits(34); // EMCAL L1 gamma trigger, low threshold, INT7 suite
  static const AliBits kEmcalL1JetHigh7    = AliBits(35); // EMCAL L1 jet trigger, high threshold, INT7 suite
  static const AliBits kEmcalL1JetLow7     = AliBits(36); // EMCAL L1 jet trigger, low threshold, INT7 suite
  static const AliBits kEMC8               = AliBits(37); // EMCAL L0 trigger, INT8 suite
  static const AliBits kEmcalL1GammaHigh8  = AliBits(38); // EMCAL L1 gamma trigger, high threshold, INT8 suite
  static const AliBits kEmcalL1GammaLow8   = AliBits(39); // EMCAL L1 gamma trigger, low threshold, INT8 suite
  static const AliBits kEmcalL1JetHigh8    = AliBits(40); // EMCAL L1 jet trigger, high threshold, INT8 suite
  static const AliBits kEmcalL1JetLow8     = AliBits(41); // EMCAL L1 jet trigger, low threshold, INT8 suite
  static const AliBits kINT7HJT            = AliBits(42); // TRD jet trigger, INT7 suite
  static const AliBits kINT7HSE            = AliBits(43); // TRD high-pt electron trigger, INT7 suite
  static const AliBits kINT7HQU            = AliBits(44); // TRD quarkonium trigger, INT7 suite
  static const AliBits kEMC7HQU            = AliBits(45); // TRD quarkonium trigger + EMCAL L0, INT7 suite
  static const AliBits kEMC7HEE            = AliBits(46); // TRD high-pt electron trigger in EMCAL acceptance + EMCAL L0, INT7 suite
  static const AliBits kINT8HJT            = AliBits(47); // TRD jet trigger, INT7 suite
  static const AliBits kINT8HSE            = AliBits(48); // TRD high-pt electron trigger, INT7 suite
  static const AliBits kINT8HQU            = AliBits(49); // TRD quarkonium trigger, INT7 suite
  static const AliBits kEMC8HQU            = AliBits(50); // TRD quarkonium trigger + EMCAL L0, INT7 suite
  static const AliBits kEMC8HEE            = AliBits(51); // TRD high-pt electron trigger in EMCAL acceptance + EMCAL L0, INT7 suite
  static const AliBits kSPI8               = AliBits(52); // Power interaction trigger, INT8 suite
  static const AliBits kSTP                = AliBits(53); // SPD topology trigger (2 hits in layer0 + 2 hits in layer1 + topology)
  static const AliBits kOMU                = AliBits(54); // TOF topology trigger (2 hits back-to-back)
  static const AliBits kCUP7               = AliBits(55); // Central barrel ultra-peripheral trigger (SPD and TOF topology, V0 veto)
}
#endif 

