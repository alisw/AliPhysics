#include "AliMUONTrackerDataSourceTypes.h"

#include "TObjString.h"
#include "TObjArray.h"

/// \cond CLASSIMP
ClassImp(AliMUONTrackerDataSourceTypes)
/// \endcond

//_____________________________________________________________________________
/// Whether type is within the aliases list
Bool_t AliMUONTrackerDataSourceTypes::IsInAliasList(const char* type, const char* aliases)
{
  Bool_t rv(kFALSE);
  TObjArray* a = TString(aliases).Tokenize(" ");
  TIter next(a);
  TObjString* str;
  TString stype(type);

  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    if ( !str->String().CompareTo(stype,TString::kIgnoreCase) )
    {
      rv = kTRUE;
      break;
    }
  }
  delete a;
  return rv;
}

//_____________________________________________________________________________
/// Whether type is of Configuration flavour
Bool_t AliMUONTrackerDataSourceTypes::IsConfig(const char* type)
{
  return IsInAliasList(type,AliasesForConfig());
}

//_____________________________________________________________________________
/// Whether type is of HV flavour
Bool_t AliMUONTrackerDataSourceTypes::IsHV(const char* type)
{
  return IsInAliasList(type,AliasesForHV());
}

//_____________________________________________________________________________
/// Whether type is of LV flavour
Bool_t AliMUONTrackerDataSourceTypes::IsLV(const char* type)
{
  return IsInAliasList(type,AliasesForLV());
}
//_____________________________________________________________________________
/// Whether type is of Occupancy flavour
Bool_t AliMUONTrackerDataSourceTypes::IsOccupancy(const char* type)
{
  return IsInAliasList(type,AliasesForOccupancy());
}

//_____________________________________________________________________________
/// Whether type is of Pedestals flavour
Bool_t AliMUONTrackerDataSourceTypes::IsPedestals(const char* type)
{
  return IsInAliasList(type,AliasesForPedestals());
}

//_____________________________________________________________________________
/// Whether type is of RejectList flavour
Bool_t AliMUONTrackerDataSourceTypes::IsRejectList(const char* type)
{
  return IsInAliasList(type,AliasesForRejectList());
}

//_____________________________________________________________________________
/// Whether type is of Status flavour
Bool_t AliMUONTrackerDataSourceTypes::IsStatus(const char* type)
{
  return IsInAliasList(type,AliasesForStatus());
}

//_____________________________________________________________________________
/// Whether type is of StatusMap flavour
Bool_t AliMUONTrackerDataSourceTypes::IsStatusMap(const char* type)
{
  return IsInAliasList(type,AliasesForStatusMap());
}
