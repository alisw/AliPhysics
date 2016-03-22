#ifndef ALIMUONTRACKERDATASOURCETYPES_H
#define ALIMUONTRACKERDATASOURCETYPES_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

/// \ingroup graphics
/// \class AliMUONTrackerDataSourceTypes
/// \brief Short names and aliases for data source types recognized by AliMUONTrackerData related classes.
///
/// \author Laurent Aphecetche <laurent.aphecetche@cern.ch>, Subatech
///

class AliMUONTrackerDataSourceTypes : public TObject
{
public:

  static const char* ShortNameForConfig() { return "CONF"; }
  static const char* AliasesForConfig() { return "CONF CONFIG Configuration"; }

  static const char* ShortNameForHV() { return "HV"; }
  static const char* AliasesForHV() { return "HV"; }

  static const char* ShortNameForLV() { return "LV"; }
  static const char* AliasesForLV() { return "LV"; }

  static const char* ShortNameForOccupancy() { return "OCC"; }
  static const char* AliasesForOccupancy() { return "OCC Occupancy"; }

  static const char* ShortNameForPedestals() { return "PED"; }
  static const char* AliasesForPedestals() { return "PED PEDESTAL Pedestals"; }

  static const char* ShortNameForRejectList() { return "RL"; }
  static const char* AliasesForRejectList() { return "RL RejectList"; }

  static const char* ShortNameForStatus() { return "STAT"; }
  static const char* AliasesForStatus() { return "STAT Status"; }

  static const char* ShortNameForStatusMap() { return "STATMAP"; }
  static const char* AliasesForStatusMap() { return "STATMAP StatusMap"; }

  static Bool_t IsConfig(const char* type);
  static Bool_t IsHV(const char* type);
  static Bool_t IsLV(const char* type);
  static Bool_t IsOccupancy(const char* type);
  static Bool_t IsPedestals(const char* type);
  static Bool_t IsRejectList(const char* type);
  static Bool_t IsStatus(const char* type);
  static Bool_t IsStatusMap(const char* type);

private:

  static Bool_t IsInAliasList(const char* type, const char* aliases);

  /// \cond CLASSIMP
  ClassDef(AliMUONTrackerDataSourceTypes,0);
  /// \endcond
};

#endif
