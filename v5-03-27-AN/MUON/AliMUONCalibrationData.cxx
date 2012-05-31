/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// $Id$

#include "AliMUONCalibrationData.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliMpDCSNamer.h"
#include "AliMpIntPair.h"
#include "AliMUONGlobalCrateConfig.h"
#include "AliMUONRegionalTriggerConfig.h"
#include "AliMUONRejectList.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerLut.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVStore.h"
#include "AliMUONVStore.h"

#include <Riostream.h>
#include <TClass.h>
#include <TMap.h>
#include <TMath.h>

//-----------------------------------------------------------------------------
/// \class AliMUONCalibrationData
///
/// For the moment, this class stores pedestals, gains, hv (for tracker)
/// and lut, masks and efficiencies (for trigger) that are fetched from the CDB.
///
/// This class is to be considered as a convenience class.
/// Its aim is to ease retrieval of calibration data from the 
/// condition database.
///
/// It acts as a "facade" to a bunch of underlying 
/// containers/calibration classes.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONCalibrationData)
/// \endcond

AliMUONVStore* AliMUONCalibrationData::fgBypassPedestals(0x0);
AliMUONVStore* AliMUONCalibrationData::fgBypassGains(0x0);

namespace  
{
  void MarkForDeletion(Int_t* indices, Int_t first, Int_t last)
  {
    for ( Int_t i = first; i <= last; ++i ) 
    {
      indices[i] = 1;
    }
  }
}

//_____________________________________________________________________________
AliMUONCalibrationData::AliMUONCalibrationData(Int_t runNumber, 
                                               Bool_t deferredInitialization) 
: TObject(), 
fIsValid(kTRUE),
fRunNumber(runNumber), 
fGains(0x0), 
fPedestals(0x0),
fHV(0x0),
fTriggerDCS(0x0),
fLocalTriggerBoardMasks(0x0),
fRegionalTriggerConfig(0x0),
fGlobalTriggerCrateConfig(0x0),
fTriggerLut(0x0),
fTriggerEfficiency(0x0),
fCapacitances(0x0),
fNeighbours(0x0),
fOccupancyMap(0x0),
fRejectList(0x0),
fConfig(0x0)
{
/// Default ctor.

  // If deferredInitialization is false, we read *all* calibrations
  // at once.
  // So when using this class to access only one kind of calibrations (e.g.
  // only pedestals), you should put deferredInitialization to kTRUE, which
  // will instruct this object to fetch the data only when neeeded.

  if ( deferredInitialization == kFALSE )
  {
    Gains();
    Pedestals();
    OccupancyMap();
    RejectList();
    HV();
    TriggerDCS();
    LocalTriggerBoardMasks(0);
    RegionalTriggerConfig();
    GlobalTriggerCrateConfig();
    TriggerLut();
    TriggerEfficiency();
    Capacitances();
    Neighbours();
    Config();
  }
}

//_____________________________________________________________________________
AliMUONCalibrationData::~AliMUONCalibrationData()
{
  /// Destructor. Note that we're the owner of our pointers if the OCDB cache
  /// is not set. Otherwise the cache is supposed to take care of them...
  if (!(AliCDBManager::Instance()->GetCacheFlag())) Reset();
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::Capacitances() const
{
  /// Create (if needed) and return the internal store for capacitances.
  
  if (!fCapacitances)
  {
    fCapacitances = CreateCapacitances(fRunNumber);
  }
  return fCapacitances;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateCapacitances(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create capa store from OCDB for a given run
  
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/Capacitances",startOfValidity));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateGains(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create a new gain store from the OCDB for a given run
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/Gains",startOfValidity));
}

//_____________________________________________________________________________
AliMUONGlobalCrateConfig*
AliMUONCalibrationData::CreateGlobalTriggerCrateConfig(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create the internal store for GlobalTriggerCrateConfig from OCDB
  
  return dynamic_cast<AliMUONGlobalCrateConfig*>(CreateObject(runNumber,"MUON/Calib/GlobalTriggerCrateConfig",startOfValidity));
}


//______________________________________________________________________________
Bool_t AliMUONCalibrationData::CheckHVGroup(TObjArray& values, Int_t first, Int_t last, Double_t& value, Int_t& slope, TString* msg)
{
  // Get the HV of the values between first and last indices
  // return the HV slope  (in Volt per second) and a message
  // Return kFALSE if we must discard the group
  //
  
  if (msg) *msg="";
  
  if ( last < first ) return kFALSE;
  if ( last - first < 2 ) return kFALSE;
  
  Double_t a(0.0);
  Double_t b(0.0);

  Float_t HVSAME(1); // 1 volts

  AliDCSValue* vfirst = static_cast<AliDCSValue*>(values.UncheckedAt(first));
  AliDCSValue* vlast = static_cast<AliDCSValue*>(values.UncheckedAt(last));

  Int_t deltaHV = TMath::Nint(TMath::Abs(vfirst->GetFloat()-vlast->GetFloat()));
  
  if ( deltaHV < HVSAME ) return kFALSE;

  for ( Int_t i = first; i <= last; ++i )
  {
    AliDCSValue* v = static_cast<AliDCSValue*>(values.UncheckedAt(i));

    Double_t y = v->GetFloat() - vfirst->GetFloat();
    Double_t x = v->GetTimeStamp() - vfirst->GetTimeStamp();
  
    a += x*y;
    b += x*x;
  }
  
  value = a/b;
  slope = value > 0 ? 1 : -1;
  value = TMath::Abs(value);
  
  UInt_t deltaTime = vlast->GetTimeStamp() - vfirst->GetTimeStamp();
  
  if (msg)
  {
    if (slope>0) (*msg) = Form("RU%d[%d:%d](%d)",TMath::Nint(value),first,last,deltaTime);
    if (slope<0) (*msg) = Form("RD%d[%d:%d](%d)",TMath::Nint(value),first,last,deltaTime);
    
    if ( TMath::Nint(value) == 0 )
    {
      // this is to protect for the few cases 
      // (see e.g. MchHvLvLeft/Chamber00Left/Quad2Sect0.actual.vMon in run 134497)
      // where we can have *lots* of values (2483 in this example) but that
      // are more or less constant...
      //
      // or simply to remove small ramps
      //
      slope = 0;
      value = (vfirst->GetFloat()+vlast->GetFloat())/2.0;
      *msg = Form("FLUCT%d[%d:%d]",TMath::Nint(value),first,last);
    }
  }
  
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliMUONCalibrationData::PatchHVValues(TObjArray& values,
                                             TString* msg)
{
  /// We do here a little bit of massaging of the HV values, if needed.
  ///
  /// The main point is to "gather" values that are within a given small amount
  /// of time (typically 60 seconds) and infer a slope from those values
  /// slope > 0 means it is a ramp-up, slope < 0 that's a ramp-down
  ///
  /// This is to avoid both the "ramp-down-before-end-of-run" and the
  /// "ramp-up-after-start-of-run" syndroms...
  ///
  /// Return kFALSE is the kind of HV (trouble) case we have here
  /// has not been identified...
  ///
  
  UInt_t DELTATIME(60); // in seconds
  Int_t IENDRU(60); // in seconds
  
  // Start by finding groups of values which are not separated (each) by more than
  // deltaTime
  
  Bool_t gather(kFALSE);
  Int_t ifirst(0);
  Int_t ilast(0);
  TObjArray groups;
  groups.SetOwner(kTRUE);
  
  for ( Int_t i = values.GetLast(); i > 0; --i ) 
  {
    AliDCSValue* vi = static_cast<AliDCSValue*>(values.UncheckedAt(i));
    AliDCSValue* vj = static_cast<AliDCSValue*>(values.UncheckedAt(i-1));

    if ( vi->GetTimeStamp() - vj->GetTimeStamp() < DELTATIME )
    {
      if ( !gather ) 
      {
        gather = kTRUE;
        ifirst = i;    
      }
      ilast=i;
    }
    else
    {
      if ( gather ) 
      {
        ilast=i;
        
        groups.Add(new AliMpIntPair(ilast,ifirst));
      }
      gather = kFALSE;
    }
  }
  
  if (gather)
  {
    groups.Add(new AliMpIntPair(0,ifirst));
  }
                 
  TIter nextGroup(&groups,kIterBackward);
  AliMpIntPair* p;
  TString internalMsg;
  Int_t ngroups(0);

  Int_t nRU(0);
  Int_t nRD(0);
  Int_t nStartRU(0);
  Int_t nEndAndShortRU(0);
  Int_t nEndRD(0);
  Int_t nTripRD(0);
  Int_t nFluct(0);
  
  while ( ( p = static_cast<AliMpIntPair*>(nextGroup()) ) )
  {
    Double_t value;
    Int_t slope;
    
    TString groupMsg;
    
    AliDebugClass(1,Form("group %d:%d",p->GetFirst(),p->GetSecond()));
    
    Bool_t ok = CheckHVGroup(values,p->GetFirst(),p->GetSecond(),value,slope,&groupMsg);
    
    if (!ok) continue;
    
    ++ngroups;
    
    if ( slope > 0 )
    {
      if ( p->GetFirst() == 0 ) 
      {
        // start with a ramp-up
        ++nStartRU;
      }
      else if ( p->GetSecond() == values.GetLast() && TMath::Nint(value) < IENDRU )
      {
        ++nEndAndShortRU;
      }
      else
      {
        // ramp-up in the middle of nowhere...
        ++nRU;
      }
    }
    else if ( slope < 0 )
    {
      if ( p->GetSecond() == values.GetLast() ) 
      {
        // end with a ramp-down
        ++nEndRD;
      }
      else
      {
        // ramp-down in the middle of nowhere
        ++nRD;
      }

      AliDCSValue* d = static_cast<AliDCSValue*>(values.At(p->GetSecond()));
      
      if ( d->GetFloat() < AliMpDCSNamer::TrackerHVOFF() )
      {
        ++nTripRD;
      }
    }
    else
    {
      ++nFluct;
    }
    
    internalMsg += groupMsg;
    internalMsg += " ";    
  }
  
  /*
   
   Once we have "decoded" the groups we try to find out which of 
   the following cases we're facing :
  
   case A = -------- = OK(1)
  
   case B = ----
                \
                \   = OK, once we have removed the ramp-down (2)
  
   case C =    ----- 
              /
             /       = OK, once we have removed the ramp-up (3)
  
   case D =    ----- 
              /     \
             /       \ = OK, once we have removed the ramp-down (2) and the ramp-up (3)
  
   case E = ----
                \
                 \____ = TRIP = BAD (here the ramp-down slope should be bigger than in case C)
  
   case F = ----         
                \      ----- = BAD (trip + ramp-up at end of run)
                 \____/   
  
   case G = fluctuations (within a range defined in CheckHVGroup...)
   
   case H =            
                   /
                  /   = ramp-up right at the end-of-run = OK (4)
            ------
   
   (1) OK means the group is identified correctly, still the value can be below ready...
   (2) ramp-down values will be removed if the ramp is indeed the last values in the serie
       i.e. it's really an end-of-run problem (otherwise it's not case B)
   (3) ramp-up values will be removed if the ramp is indeed the first values in the serie
       i.e. it's really a start-of-run problem (otherwise it's not case C)
   (4) OK if short enough...
   
   Any other case is unknown and we'll :
   a) return kFALSE
   b) assume the channel is OFF.
  
  
  */
  
  AliDebugClass(1,Form("msg=%s ngroupds=%d",internalMsg.Data(),ngroups));
  AliDebugClass(1,Form("nRU %d nRD %d nStartRU %d nEndRD %d nTripRD %d nFluct %d",
                       nRU,nRD,nStartRU,nEndRD,nTripRD,nFluct));
  
  TString hvCase("OTHER");
  int dummy(0),a(-1),b(-1);
  char r[81];
  Int_t nvalues = values.GetSize();  
  Int_t* indices = new Int_t[nvalues];
  memset(indices,0,nvalues*sizeof(Int_t));
         
  AliDCSValue* vfirst = static_cast<AliDCSValue*>(values.UncheckedAt(0));
  AliDCSValue* vlast = static_cast<AliDCSValue*>(values.UncheckedAt(values.GetLast()));

  UInt_t meanTimeStamp = ( vfirst->GetTimeStamp() + vlast->GetTimeStamp() ) / 2;
  
  if ( ngroups == 0 ) 
  {
    hvCase = "A"; 
  }
  else if ( nTripRD > 0 )
  {
    if ( nRU > 0 && nRD > 0 )
    {
      hvCase = "F";
    }
    else
    {
      hvCase = "E";
    }
    internalMsg += "TRIP ";
    MarkForDeletion(indices,0,values.GetLast());
    values.Add(new AliDCSValue(static_cast<Float_t>(0),meanTimeStamp));
  }
  else if ( nStartRU > 0 && nRU == 0 && nRD == 0 && nEndRD == 0 )
  {
    hvCase = "C";
    sscanf(internalMsg.Data(),"RU%10d[%10d:%10d]%80s",&dummy,&a,&b,r);
    MarkForDeletion(indices,a,b);
  }
  else if ( nStartRU > 0 && nEndRD > 0 && nRD == 0 && nRU == 0 )
  {
    hvCase = "D";
    sscanf(internalMsg.Data(),"RU%10d[%10d:%10d]%80s",&dummy,&a,&b,r);    
    MarkForDeletion(indices,a,b-1);
    Int_t i = internalMsg.Index("RD",strlen("RD"),0,TString::kExact);
    sscanf(internalMsg(i,internalMsg.Length()-i).Data(),
           "RD%10d[%10d:%10d]%80s",&dummy,&a,&b,r);    
    MarkForDeletion(indices,a+1,b);
  }
  else if ( nEndRD > 0 && nStartRU == 0 && nRU == 0 && nRD == 0 )
  {
    hvCase = "B";
    Int_t i = internalMsg.Index("RD",strlen("RD"),0,TString::kExact);
    sscanf(internalMsg(i,internalMsg.Length()-i).Data(),
           "RD%10d[%10d:%10d]%80s",&dummy,&a,&b,r);    
    MarkForDeletion(indices,a,b);
  }
  else if ( nFluct > 0 )
  {
    hvCase = "G";
    TObjArray* af = internalMsg.Tokenize(" ");
    TIter next(af);
    TObjString* str;
    while ( ( str = static_cast<TObjString*>(next()) ) )
    {
      TString s(str->String());
      if ( s.BeginsWith("FLUCT") )
      {
        sscanf(s.Data(),"FLUCT%d[%d:%d]",&dummy,&a,&b);
        MarkForDeletion(indices,a,b);
      }
    }
    delete af;
  }
  else if ( nEndAndShortRU > 0 && nStartRU == 0 && nRU == 0 && nRD == 0 && nEndRD == 0 )
  {
    hvCase = "H";
    sscanf(internalMsg.Data(),"RU%10d[%10d:%10d]%80s",&dummy,&a,&b,r);
    MarkForDeletion(indices,a,b);
  }
  else
  {
    // last chance... 
    // here we know it's not a trip, so let's assume everything is OK
    // if first and last value are in the same ballpark

    const Double_t HVFLUCT(20); // volts
    
    if ( TMath::Abs(vfirst->GetFloat() - vlast->GetFloat()) < HVFLUCT )
    {
      hvCase = "Z";
    }
    MarkForDeletion(indices,1,nvalues-1);
  }
  
  for ( Int_t i = 0; i < nvalues; ++i ) 
  {
    if ( indices[i] )
    {
      values.RemoveAt(i);
    }
  }
  
  values.Compress();

  delete[] indices;
  
  if ( !values.GetEntries() )
  {
    AliErrorClass(Form("No value left after patch... Check that !!! initial # of values=%d msg=%s",
                       nvalues,internalMsg.Data()));
    hvCase = "OTHER";
  }

  // take the max of the remaining values
  TIter nextA(&values);
  AliDCSValue* val;
  Float_t maxval(-9999);
  
  while ( ( val = static_cast<AliDCSValue*>(nextA()) ) )
  {
    if ( val->GetFloat() > maxval )
    {
      maxval = val->GetFloat();
    }
  }
  
  values.Clear();
  
  values.Add(new AliDCSValue(maxval,meanTimeStamp));
  
  // once the case is inferred, add a "CASE:%10d",hvCase.Data()
  // to the msg
  // so we can them sum up for all channels and get a summary per run...
  
  internalMsg += Form("CASE:%s",hvCase.Data());
 
  if (msg) *msg = internalMsg.Data();
  
  return hvCase=="OTHER" ? kFALSE : kTRUE;
}

//_____________________________________________________________________________
TMap*
AliMUONCalibrationData::CreateHV(Int_t runNumber, 
                                 Int_t* startOfValidity, 
                                 Bool_t patched,
                                 TList* messages)
{
  /// Create a new HV map from the OCDB for a given run
  TMap* hvMap = dynamic_cast<TMap*>(CreateObject(runNumber,"MUON/Calib/HV",startOfValidity));

  if (!hvMap) return 0x0;
  
  if (patched)
  {
    TIter next(hvMap);
    TObjString* hvChannelName;
    
    while ( ( hvChannelName = static_cast<TObjString*>(next()) ) )
    {
      TString name(hvChannelName->String());
      
      if ( name.Contains("sw") ) continue; // skip switches
      
      TPair* hvPair = static_cast<TPair*>(hvMap->FindObject(name.Data()));
      TObjArray* values = static_cast<TObjArray*>(hvPair->Value());
      if (!values)
      {
        AliErrorClass(Form("Could not get values for alias %s",name.Data()));
      }
      else
      {
        TString msg;
        
        AliDebugClass(1,Form("channel %s",name.Data()));
        Bool_t ok = PatchHVValues(*values,&msg);
        
        if ( messages ) 
        {
          messages->Add(new TObjString(Form("%s:%s",hvChannelName->String().Data(),msg.Data())));
        }
        
        if (!ok)
        {
          AliErrorClass(Form("PatchHVValue was not successfull ! This is serious ! "
                             "You'll have to check the logic for channel %s in run %09d",
                             name.Data(),runNumber));
        }
      }
    }
    
  }
  
  if ( messages ) 
  {
    Int_t a(0),b(0),c(0),d(0),e(0),f(0),g(0),h(0),u(0),z(0);
    TIter next(messages);
    TObjString* msg;
    char hvCase('u');
    
    while ( ( msg = static_cast<TObjString*>(next()) ) )
    {
      Int_t i = msg->String().Index("CASE",strlen("CASE"),0,TString::kExact);
      
      if ( i >= 0 )
      {
        sscanf(msg->String()(i,msg->String().Length()-i).Data(),"CASE:%10c",&hvCase);
      }

      switch (hvCase)
      {
        case 'A': ++a; break;
        case 'B': ++b; break;
        case 'C': ++c; break;
        case 'D': ++d; break;
        case 'E': ++e; break;
        case 'F': ++f; break;
        case 'G': ++g; break;
        case 'H': ++h; break;
        case 'Z': ++z; break;
        default: ++u; break;
      }
    }
    
    messages->Add(new TObjString(Form("SUMMARY : # of cases A(%3d) B(%3d) C(%3d) D(%3d) E(%3d) F(%3d) G(%3d) H(%3d) Z(%3d) OTHER(%3d)",
                                      a,b,c,d,e,f,g,h,z,u)));
  }
  
  return hvMap;
}

//_____________________________________________________________________________
TMap*
AliMUONCalibrationData::CreateTriggerDCS(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create a new Trigger HV and curent map from the OCDB for a given run
  return dynamic_cast<TMap*>(CreateObject(runNumber,"MUON/Calib/TriggerDCS",startOfValidity));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateLocalTriggerBoardMasks(Int_t runNumber, Int_t* startOfValidity)
{
  /// Get the internal store for LocalTriggerBoardMasks from OCDB
  
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/LocalTriggerBoardMasks",startOfValidity));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateNeighbours(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create a neighbour store from the OCDB for a given run
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/Neighbours",startOfValidity));
}

//_____________________________________________________________________________
TObject*
AliMUONCalibrationData::CreateObject(Int_t runNumber, const char* path, Int_t* startOfValidity)
{
  /// Access the CDB for a given path (e.g. MUON/Calib/Pedestals),
  /// and return the corresponding TObject.
  
  AliCodeTimerAutoClass(Form("%09d : %s",runNumber,path),0);
  
  AliCDBManager* man = AliCDBManager::Instance();
  
  AliCDBEntry* entry =  man->Get(path,runNumber);
  
  if (entry)
  {
		if ( startOfValidity ) *startOfValidity = entry->GetId().GetFirstRun();
		
    TObject* object = entry->GetObject();
    if (!(man->GetCacheFlag()))
    {
      entry->SetOwner(kFALSE);
      delete entry;      
    }
//    else
//    {
//      entry->SetOwner(kTRUE); //FIXME : this should be done but is causing problems with RecoParams at the end of the reco : investigate why...
//    }
    return object;
  }
	else
	{
		if ( startOfValidity )  *startOfValidity = AliCDBRunRange::Infinity();
  }
	
  {
    
    AliCodeTimerAutoClass(Form("Failed to get %s for run %09d",path,runNumber),1);

  }
  
  return 0x0;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateOccupancyMap(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create a new occupancy map store from the OCDB for a given run
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/OccupancyMap",startOfValidity));
}

//_____________________________________________________________________________
AliMUONRejectList*
AliMUONCalibrationData::CreateRejectList(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create a new rejectlist store from the OCDB for a given run
  return dynamic_cast<AliMUONRejectList*>(CreateObject(runNumber,"MUON/Calib/RejectList",startOfValidity));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreatePedestals(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create a new pedestal store from the OCDB for a given run
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/Pedestals",startOfValidity));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::CreateConfig(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create a new config store from the OCDB for a given run
  return dynamic_cast<AliMUONVStore*>(CreateObject(runNumber,"MUON/Calib/Config",startOfValidity));
}


//_____________________________________________________________________________
AliMUONRegionalTriggerConfig*
AliMUONCalibrationData::CreateRegionalTriggerConfig(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create the internal store for RegionalTriggerConfig from OCDB
  
  return dynamic_cast<AliMUONRegionalTriggerConfig*>(CreateObject(runNumber,"MUON/Calib/RegionalTriggerConfig",startOfValidity));
}

//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells* 
AliMUONCalibrationData::CreateTriggerEfficiency(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create trigger efficiency object from OCBD
  
  return dynamic_cast<AliMUONTriggerEfficiencyCells*>(CreateObject(runNumber,"MUON/Calib/TriggerEfficiency",startOfValidity));
}

//_____________________________________________________________________________
AliMUONTriggerLut* 
AliMUONCalibrationData::CreateTriggerLut(Int_t runNumber, Int_t* startOfValidity)
{
  /// Create trigger LUT from OCDB
  
  return dynamic_cast<AliMUONTriggerLut*>(CreateObject(runNumber,"MUON/Calib/TriggerLut",startOfValidity));
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::Gains() const
{
  /// Create (if needed) and return the internal store for gains.
  if (fgBypassGains) return fgBypassGains;
  
  if (!fGains)
  {
    fGains = CreateGains(fRunNumber);
  }
  return fGains;
}
 
//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONCalibrationData::Gains(Int_t detElemId, Int_t manuId) const
{
/// Return the gains for a given (detElemId, manuId) pair
/// Note that, unlike the DeadChannel case, if the result is 0x0, that's an
/// error (meaning that we should get gains for all channels).

  AliMUONVStore* gains = Gains();
  if (!gains)
  {
    return 0x0;
  }
  
  return static_cast<AliMUONVCalibParam*>(gains->FindObject(detElemId,manuId));
}

//_____________________________________________________________________________
AliMUONGlobalCrateConfig* 
AliMUONCalibrationData::GlobalTriggerCrateConfig() const
{
  /// Return the config for the global trigger board.
  
  if (!fGlobalTriggerCrateConfig)
  {
    fGlobalTriggerCrateConfig = CreateGlobalTriggerCrateConfig(fRunNumber);
  }
  return fGlobalTriggerCrateConfig;
}


//_____________________________________________________________________________
TMap*
AliMUONCalibrationData::HV(Bool_t patched) const
{
  /// Return the calibration for a given (detElemId, manuId) pair
  
  if (!fHV)
  {
    fHV = CreateHV(fRunNumber,0,patched);
  }
  return fHV;
}

//_____________________________________________________________________________
TMap*
AliMUONCalibrationData::TriggerDCS() const
{
  /// Return the calibration for a given (detElemId, manuId) pair
  
  if (!fTriggerDCS)
  {
    fTriggerDCS = CreateTriggerDCS(fRunNumber);
  }
  return fTriggerDCS;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::Neighbours() const
{
  /// Create (if needed) and return the internal store for neighbours.
  if (!fNeighbours)
  {
    fNeighbours = CreateNeighbours(fRunNumber);
  }
  return fNeighbours;
}

//_____________________________________________________________________________
AliMUONVCalibParam* 
AliMUONCalibrationData::LocalTriggerBoardMasks(Int_t localBoardNumber) const
{
/// Return the masks for a given trigger local board.

  if (!fLocalTriggerBoardMasks)
  {
    fLocalTriggerBoardMasks = CreateLocalTriggerBoardMasks(fRunNumber);
  }

  if ( fLocalTriggerBoardMasks ) 
  {
    AliMUONVCalibParam* ltbm = 
      static_cast<AliMUONVCalibParam*>(fLocalTriggerBoardMasks->FindObject(localBoardNumber));
    if (!ltbm)
    {
      AliError(Form("Could not get mask for localBoardNumber=%d",localBoardNumber));
    }
    return ltbm;  
  }
  return 0x0;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::OccupancyMap() const
{
  /// Get occupancy map
  if (!fOccupancyMap)
  {
    fOccupancyMap = CreateOccupancyMap(fRunNumber);
  }
  return fOccupancyMap;
}

//_____________________________________________________________________________
AliMUONRejectList*
AliMUONCalibrationData::RejectList() const
{
  /// Get reject list
  if (!fRejectList)
  {
    fRejectList = CreateRejectList(fRunNumber);
  }
  return fRejectList;
}

//_____________________________________________________________________________
void
AliMUONCalibrationData::BypassStores(AliMUONVStore* ped, AliMUONVStore* gain)
{
  /// Force the use of those pedestals and gains
  fgBypassPedestals = ped;
  fgBypassGains = gain;
  
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::Pedestals() const
{
  /// Return pedestals
  
  if (fgBypassPedestals) return fgBypassPedestals;
  
  if (!fPedestals)
  {
    fPedestals = CreatePedestals(fRunNumber);
  }
  return fPedestals;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONCalibrationData::Config() const
{
  /// Return config
  
  if (!fConfig)
  {
    fConfig = CreateConfig(fRunNumber);
  }
  return fConfig;
}

//_____________________________________________________________________________
AliMUONVCalibParam*
AliMUONCalibrationData::Pedestals(Int_t detElemId, Int_t manuId) const
{
  /// Return the pedestals for a given (detElemId, manuId) pair.
  /// A return value of 0x0 is considered an error, meaning we should get
  /// pedestals for all channels.
  
  AliMUONVStore* pedestals = Pedestals();
  if (!pedestals) 
  {
    return 0x0;
  }
  
  return static_cast<AliMUONVCalibParam*>(pedestals->FindObject(detElemId,manuId));
}

//_____________________________________________________________________________
void
AliMUONCalibrationData::Print(Option_t*) const
{
  /// A very basic dump of our guts.

  cout << "RunNumber " << RunNumber()
  << " fGains=" << fGains
  << " fPedestals=" << fPedestals
  << " fConfig=" << fConfig
  << " fHV=" << fHV
  << " fTriggerDCS=" << fTriggerDCS
  << " fLocalTriggerBoardMasks=" << fLocalTriggerBoardMasks
  << " fRegionalTriggerConfig=" << fRegionalTriggerConfig
  << " fGlobalTriggerCrateConfig=" << fGlobalTriggerCrateConfig
  << " fTriggerLut=" << fTriggerLut
  << endl;
}


//_____________________________________________________________________________
AliMUONRegionalTriggerConfig* 
AliMUONCalibrationData::RegionalTriggerConfig() const
{
  /// Return the config for the regional trigger board.
  
  if (!fRegionalTriggerConfig)
  {
    fRegionalTriggerConfig = CreateRegionalTriggerConfig(fRunNumber);
    }
  return fRegionalTriggerConfig;
}


//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells*
AliMUONCalibrationData::TriggerEfficiency() const
{
/// Return the trigger efficiency.

  if (!fTriggerEfficiency)
  {
    fTriggerEfficiency = CreateTriggerEfficiency(fRunNumber);
  }
  return fTriggerEfficiency;
}


//_____________________________________________________________________________
AliMUONTriggerLut*
AliMUONCalibrationData::TriggerLut() const
{
/// Return the trigger look up table.

  if (!fTriggerLut)
  {
    fTriggerLut = CreateTriggerLut(fRunNumber);
  }
  return fTriggerLut;
}

//_____________________________________________________________________________
void
AliMUONCalibrationData::Reset()
{
/// Reset all data

  AliCodeTimerAuto("",0);
  
  delete fConfig;
  fConfig = 0x0;
  delete fPedestals;
  fPedestals = 0x0;
  delete fGains;
  fGains = 0x0;
  delete fHV;
  fHV = 0x0;
  delete fTriggerDCS;
  fTriggerDCS = 0x0;
  delete fLocalTriggerBoardMasks;
  fLocalTriggerBoardMasks = 0x0;
  delete fRegionalTriggerConfig;
  fRegionalTriggerConfig = 0x0;
  delete fGlobalTriggerCrateConfig;
  fGlobalTriggerCrateConfig = 0x0;
  
  delete fTriggerLut;
  fTriggerLut = 0x0;
  delete fTriggerEfficiency;
  fTriggerEfficiency = 0x0;
  delete fCapacitances;
  fCapacitances = 0x0;
  delete fNeighbours;
  fNeighbours = 0x0;
}

//_____________________________________________________________________________
void
AliMUONCalibrationData::Check(Int_t runNumber)
{
  /// Self-check to see if we can read all data for a given run 
  /// from the current OCDB...
  
  if ( ! CreateCapacitances(runNumber) )
  {
    AliErrorClass("Could not read capacitances");
  }
  else
  {
    AliInfoClass("Capacitances read OK");
  }

  if ( ! CreateGains(runNumber) ) 
  {
    AliErrorClass("Could not read gains");
  }
  else
  {
    AliInfoClass("Gains read OK");
  }

  if ( ! CreateGlobalTriggerCrateConfig(runNumber) ) 
  {
    AliErrorClass("Could not read Trigger Crate Config");
  }
  else
  {
    AliInfoClass("TriggerBoardMasks read OK");
  }

  if ( !  CreateHV(runNumber) )
  {
    AliErrorClass("Could not read HV");
  }
  else
  {
    AliInfoClass("HV read OK");
  }

  if ( !  CreateTriggerDCS(runNumber) )
  {
    AliErrorClass("Could not read Trigger HV and Currents");
  }
  else
  {
    AliInfoClass("Trigger HV and Currents read OK");
  }

  if ( ! CreateNeighbours(runNumber) )
  {
    AliErrorClass("Could not read Neighbours");
  }
  else
  {
    AliInfoClass("Neighbours read OK");
  }

  if ( !  CreateLocalTriggerBoardMasks(runNumber) )
  {
    AliErrorClass("Could not read LocalTriggerBoardMasks");
  }
  else
  {
    AliInfoClass("LocalTriggerBoardMasks read OK");
  }
  
  if ( ! CreatePedestals(runNumber) )
  {
    AliErrorClass("Could not read pedestals");
  }
  else
  {
    AliInfoClass("Pedestals read OK");
  }

  if ( ! CreateConfig(runNumber) )
  {
    AliErrorClass("Could not read config");
  }
  else
  {
    AliInfoClass("Config read OK");
  }
  
  if ( ! CreateRegionalTriggerConfig(runNumber) )
  {
    AliErrorClass("Could not read RegionalTriggerConfig");
  }
  else
  {
    AliInfoClass("RegionalTriggerBoardMasks read OK");
  }
  
  if ( ! CreateTriggerLut(runNumber) )
  {
    AliErrorClass("Could not read TriggerLut");
  }
  else
  {
    AliInfoClass("TriggerLut read OK");
  }

  if ( ! CreateTriggerEfficiency(runNumber) )
  {
    AliErrorClass("Could not read TriggerEfficiency");
  }
  else    
  {
    AliInfoClass("TriggerEfficiency read OK");
  }
}
