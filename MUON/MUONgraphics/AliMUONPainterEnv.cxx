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

#include "AliMUONPainterEnv.h"

#include <TEnv.h>
#include <TSystem.h>
#include "TObjString.h"
#include "TObjArray.h"
#include "Riostream.h"
#include "AliLog.h"
#include <cassert>
#include <limits>

///\class AliMUONPainterEnv
///
/// A resource file handling class.
///
/// Used to get some things persistent between two sessions of the mchview
/// program.
///
/// Main usage is for data sources.
///
/// Syntax is :
///
/// NumberOfDataSources: n
/// DataSource.%d: ID;NAME;RANGES
///
/// ID = ORIGIN|URI|TYPE
/// URI = depends on source
///
/// where ORIGIN might be OCDB or FILE
/// (see AliMUONPainterDataSourceFrame::CreateOCDBDataSource())
/// (see AliMUONPainterDataSourceFrame::CreateRawDataSource())
/// (see AliMUONPainterDataSourceFrame::CreateACFDataSource())
///
/// URI is a source dependent set of information (separated by |) needed to recreate the source, e.g.
/// local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB|216253
///
/// TYPE is the data source type (e.g. PED, RAW, CALZ, see AliMUONPainterDataSourceTypes)
///
/// NAME is a (unique) shortname for display, e.g. PED216253(alien)
///
/// RANGES is a set of "|"-separated range definition. A range definition indicates,
/// for one of the dimension of the data source, the min and max value to be plotted.
/// Syntax is DIMNAME/XMIN/XMAX
///
/// The full line ID NAME RANGES is called the data source descriptor
///
/// Examples :
///
/// DataSource.0: OCDB|local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB|228901|PED;PED228901(cvmfs);Mean/0/500|Sigma/0/2
///
///
/// Apart from data source descriptions, the following keywords are recognized :
///
/// disableAutoPedCanvas: 0 | 1
///
/// which is used to suppress the automatic creation of summary canvases when opening a pedestal source
///
///
/// defaultRange: TYPE1/DIMNAME1/XMIN1/XMAX1|TYPE1/DIMNAME2/XMIN2/XMAX2
///
/// example of defaultRange :
///
/// defaultRange: PED/Mean/0/500|PED/Sigma/0/2|OCC/occ/0/0.01|HV/mean of HV/0/1650|
///
/// which is used to define default ranges for some dimensions of some source types
///
///
/// \author Laurent Aphecetche, Subatech
///

///\cond CLASSIMP
ClassImp(AliMUONPainterEnv)
///\endcond

const char* AliMUONPainterEnv::fgkNumberOfDataSourcesKey = "NumberOfDataSources";
const char* AliMUONPainterEnv::fgkDataSourceKey = "DataSource.%d";
const char* AliMUONPainterEnv::fgkDefaultRangeKey = "defaultRange";
const char* AliMUONPainterEnv::fgkDisableAutoPedCanvasKey = "disableAutoPedCanvas";
const char* AliMUONPainterEnv::fgkSeparatorWithinRange = "/";
const char* AliMUONPainterEnv::fgkSeparatorWithinPart = "|";
const char* AliMUONPainterEnv::fgkSeparatorBetweenDescriptorParts=";";

//_____________________________________________________________________________
AliMUONPainterEnv::AliMUONPainterEnv(const char* resourceFile)
: fEnv(new TEnv(resourceFile))
{
  /// Ctor
  
  if ( String(fgkDisableAutoPedCanvasKey,"").Length()==0 )
  {
    // no default for that one, assume it should be 1
    Set(fgkDisableAutoPedCanvasKey,1);
  }
  
  if ( String(fgkDefaultRangeKey,"").Length()==0 )
  {
    // no default ranges defined at all, let's make some reasonable guesses
    
    AliInfo("No default ranges defined in the resource file, generating a few...");
    
    SetDefaultRange("PED","Mean",0.0,500.0);
    SetDefaultRange("PED","Sigma",0.0,5.0);
    SetDefaultRange("OCC","occ",0.0,0.01);
    SetDefaultRange("HV","Mean of HV",0,1650);
    
    Save();
  }
}

//_____________________________________________________________________________
AliMUONPainterEnv::~AliMUONPainterEnv()
{
  /// dtor
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::AddDataSource(const char* dataSourceDescriptor)
{
  /// Add (or replace) a data source to the list of data sources we handle
  /// Return the full descriptor of the data source
  
  TString name = Descriptor2Name(dataSourceDescriptor);
  TString sid = Descriptor2ID(dataSourceDescriptor);
  TString uri = ID2URI(sid);

  for ( Int_t i = 0; i < NumberOfDataSources(); ++i )
  {
    TString desc = DataSourceDescriptor(i);
    TString sid = Descriptor2ID(desc);
    TString urii = ID2URI(sid);
    
    if ( name == Descriptor2Name(desc) && uri == urii )
      {
        // already there, just replace
        Set(Form(fgkDataSourceKey,i),dataSourceDescriptor);
        return dataSourceDescriptor;
      }
  }
  
  // new source
  
  Int_t n = NumberOfDataSources();
  
  Set(fgkNumberOfDataSourcesKey,n+1);

  Set(Form(fgkDataSourceKey,n),dataSourceDescriptor);
  
  ForceDataSourceToDefaultRange(name);
  
  Save();
  
  return DataSourceDescriptor(name);
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::DataSourceDescriptor(Int_t index) const
{
  if ( index >= 0 && index < NumberOfDataSources() )
  {
    return String(Form(fgkDataSourceKey,index),"");
  }
  return "";
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::DataSourceDescriptor(const char* dataSourceName) const
{
  Int_t index = GetDataSourceIndex(dataSourceName);
  return DataSourceDescriptor(index);
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::DataSourceID(const char* dataSourceName) const
{
  return Descriptor2ID(DataSourceDescriptor(dataSourceName));
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::DataSourceName(const char* dataSourceName) const
{
  return Descriptor2Name(DataSourceDescriptor(dataSourceName));
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::DataSourceOrigin(const char* dataSourceName) const
{
  return ID2Origin(DataSourceID(dataSourceName));
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::DataSourceRanges(const char* dataSourceName) const
{
  return Descriptor2Ranges(DataSourceDescriptor(dataSourceName));
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::DataSourceType(const char* dataSourceName) const
{
  return ID2Type(DataSourceID(dataSourceName));
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::DataSourceURI(const char* dataSourceName) const
{
  return ID2URI(DataSourceID(dataSourceName));
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::Descriptor2Name(const char* dataSourceDescriptor)
{
  /// From descriptor = ID NAME RANGES return the NAME part

  TString rv;
  
  TObjArray* a = TString(dataSourceDescriptor).Tokenize(fgkSeparatorBetweenDescriptorParts);
  if ( a->GetLast() >= 1 )
  {
    rv = static_cast<TObjString*>(a->At(1))->String();
  }
  delete a;
  return rv;
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::Descriptor2Ranges(const char* dataSourceDescriptor)
{
  /// From descriptor = ID NAME RANGES return the RANGES part

  return TupleLast(dataSourceDescriptor,fgkSeparatorBetweenDescriptorParts);
}


//_____________________________________________________________________________
TString AliMUONPainterEnv::Descriptor2ID(const char* dataSourceDescriptor)
{
  /// From descriptor = ID NAME RANGES return the ID part
  
  return TupleFirst(dataSourceDescriptor,fgkSeparatorBetweenDescriptorParts);
}

//_____________________________________________________________________________
Double_t
AliMUONPainterEnv::Double(const char* resourceName, Double_t defaultValue) const
{
  /// Retrieve the value associated with a given source, as a double
  
  return fEnv->GetValue(resourceName,defaultValue);
}

//__________________________________________________________________________________________________
void AliMUONPainterEnv::ForceDataSourceToDefaultRange(const char* dataSourceName, const char* dimensionName)
{
  TString desc = DataSourceDescriptor(dataSourceName);
  
  if ( desc.Length()==0 )
  {
    AliError(Form("Could not find data source named : %s",dataSourceName));
    return;
  }
  
  TString dsname = Descriptor2Name(desc);
  
  assert(dsname==dataSourceName);
  
  TString dataSourceType = DataSourceType(dataSourceName);

  TString newDescriptor = Descriptor2ID(desc);
  
  newDescriptor += fgkSeparatorBetweenDescriptorParts;
  newDescriptor += dsname;
  
  // at this stage the descriptor is back to "URI NAME"
  // now we add the existing ranges and replace the one corresponding to dim
  // (or all if dim="") by the default range(s) (if any is defined)

  newDescriptor += fgkSeparatorBetweenDescriptorParts;

  if ( strlen(dimensionName) == 0 )
  {
    // use all the default ranges for this data source
    TString defaultRange = GetDefaultRange(dataSourceType,"");
    newDescriptor += defaultRange;
  }
  else
  {
    // replace only the default range for dimensionName
    TString ranges = DataSourceRanges(dataSourceName);
  
    TObjArray* r = ranges.Tokenize(fgkSeparatorWithinPart);
    TString defaultRange = GetDefaultRange(dataSourceType,dimensionName);

    for ( Int_t i = 0; i <= r->GetLast(); ++i )
    {
      TObjString* sr = static_cast<TObjString*>(r->At(i));

      TString val = sr->String();
    
      if ( TupleFirst(sr->String(),fgkSeparatorWithinRange)==dimensionName )
      {
        if ( defaultRange.Length() )
        {
          val = defaultRange;
        }
      }
      newDescriptor += fgkSeparatorWithinPart;
      newDescriptor += val;
    }
    
    delete r;
  }
  
  Int_t index = GetDataSourceIndex(dataSourceName);
  
  Set(Form(fgkDataSourceKey,index),newDescriptor.Data());
}

//_____________________________________________________________________________
Int_t AliMUONPainterEnv::GetDataSourceIndex(const char* dataSourceName) const
{
  for ( Int_t i = 0; i < NumberOfDataSources(); ++i )
  {
    TString desc = DataSourceDescriptor(i);
    
    if ( Descriptor2Name(desc) == dataSourceName )
    {
      return i;
    }
  }
  return -1;
}

//_____________________________________________________________________________
TString
AliMUONPainterEnv::GetDefaultRange(const char* dataSourceType, const char* dimensionName) const
{
  /// Get the default range for a given dimension
  /// if dimensionName="" return all the defined default ranges for this data source type
  
  TString defaultRanges = String(fgkDefaultRangeKey,"");

  TString rv;
  
  TObjArray* r = defaultRanges.Tokenize(fgkSeparatorWithinPart);

  for ( Int_t ir = 0; ir <= r->GetLast(); ++ir )
  {
    TObjString* s = static_cast<TObjString*>(r->At(ir));

    TObjArray* a = s->String().Tokenize(fgkSeparatorWithinRange);
    if ( a->GetLast() != 3 )
    {
      AliError(Form("Malformed %s string : %s",fgkDefaultRangeKey,s->String().Data()));
      continue;
    }

    TString typeName = static_cast<TObjString*>(a->At(0))->String();

    if (typeName==dataSourceType)
    {
      TString dimName = static_cast<TObjString*>(a->At(1))->String();

      if ( dimName == dimensionName || strlen(dimensionName) == 0 )
      {
        for ( Int_t i = 1; i <= 3; ++i )
        {
          rv += static_cast<TObjString*>(a->At(i))->String();
          if ( i < 3 )
          {
            rv += fgkSeparatorWithinRange;
          }
        }
        if ( ir < r->GetLast() )
        {
          rv += fgkSeparatorWithinPart;
        }
      }
    }
    delete a;
  }

  delete r;
  
  return rv;
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::ID2Origin(const char* dataSourceID)
{
  /// From ID = ORIGIN|URI|TYPE return ORIGIN
  return TupleFirst(dataSourceID,fgkSeparatorWithinPart);
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::ID2URI(const char* dataSourceID)
{
  /// From ID = ORIGIN|URI|TYPE return URI
  return TupleMiddle(dataSourceID,fgkSeparatorWithinPart);
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::ID2Type(const char*  dataSourceID)
{
  /// From ID = ORIGIN|URI|TYPE return TYPE
  return TupleLast(dataSourceID,fgkSeparatorWithinPart);
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterEnv::Integer(const char* resourceName, Int_t defaultValue) const
{
  /// Retrieve the value associated with a given source, as an integer

  return fEnv->GetValue(resourceName,defaultValue);
}

//_____________________________________________________________________________
Int_t
AliMUONPainterEnv::NumberOfDataSources() const
{
  return Integer(fgkNumberOfDataSourcesKey);
}

//_____________________________________________________________________________
void
AliMUONPainterEnv::Print(Option_t* opt) const
{
  fEnv->Print(opt);
}

//_____________________________________________________________________________
Bool_t AliMUONPainterEnv::Ranges2DimensionRange(const char* ranges,
                                                const char* dimensionName,
                                                Double_t& xmin,
                                                Double_t& xmax)
{
  xmin = std::numeric_limits<double>::max();
  xmax = std::numeric_limits<double>::min();
  
  TString dims;
  
  TObjArray* a = TString(ranges).Tokenize(fgkSeparatorWithinPart);
  
  for ( Int_t i = 0; i <= a->GetLast(); ++i )
  {
    TString drange = static_cast<TObjString*>(a->At(i))->String();
    // drange should be of the form dimname/xmin/xmax
    
    TString bname = TupleFirst(drange,fgkSeparatorWithinRange);
      
    if ( bname == dimensionName )
    {
      xmin = TupleMiddle(drange,fgkSeparatorWithinRange).Atof();
      xmax = TupleLast(drange,fgkSeparatorWithinRange).Atof();
      break;
    }
  }
  
  return (xmin<xmax);
}

//_____________________________________________________________________________
void
AliMUONPainterEnv::Save()
{
  /// Save the resource file
  fEnv->WriteFile(gSystem->ExpandPathName(Form("$HOME/%s",fEnv->GetRcName())));
}

//_____________________________________________________________________________
void 
AliMUONPainterEnv::Set(const char* resourceName, Int_t value)
{
  /// Set an integer resource
  fEnv->SetValue(resourceName,Form("%d",value));
}

//_____________________________________________________________________________
void 
AliMUONPainterEnv::Set(const char* resourceName, const char* value)
{
  /// Set a string resource
  fEnv->SetValue(resourceName,value);
}

//_____________________________________________________________________________
void 
AliMUONPainterEnv::Set(const char* resourceName, Double_t value)
{
  /// Set a double resource

  fEnv->SetValue(resourceName,Form("%g",value));
}

//_____________________________________________________________________________
void AliMUONPainterEnv::SetDefaultRange(const char* dataSourceType, const char* dimensionName,
                                        Double_t xmin, Double_t xmax)
{
  /// Set the default range of a given dimension for a data source type
  
  TString defaultRanges = String(fgkDefaultRangeKey,"");
  
  TObjArray* a = defaultRanges.Tokenize(fgkSeparatorWithinPart);
  TString newRange;
  TString def = Form("%s%s%s%s%g%s%g",
                     dataSourceType,fgkSeparatorWithinRange,
                     dimensionName,fgkSeparatorWithinRange,
                     xmin,fgkSeparatorWithinRange,
                     xmax);
  Bool_t replace(kFALSE);
  
  for ( Int_t i = 0 ; i<= a->GetLast(); ++i )
  {
    TObjString* s = static_cast<TObjString*>(a->At(i));

    TString dim = TupleSecond(s->String(),fgkSeparatorWithinRange);
    
    if ( dim == dimensionName )
    {
      replace = kTRUE;
      newRange += Form("%s%s%s%s%g%s%g",
                       dataSourceType,fgkSeparatorWithinRange,
                       dimensionName,fgkSeparatorWithinRange,
                       xmin,fgkSeparatorWithinRange,
                       xmax);
    }
    else
    {
      newRange += s->String();
    }
    
    newRange += fgkSeparatorWithinPart;
  }

  if (!replace)
  {
    newRange += def;
  }

  delete a;

  Set(fgkDefaultRangeKey,newRange.Data());
}

//_____________________________________________________________________________
void AliMUONPainterEnv::SetDimensionRange(const char* dataSourceName, const char* dimensionName, Double_t xmin, Double_t xmax)
{
  /// Replace i-th source dimension range by xmin,xmax
  
  TString desc = DataSourceDescriptor(dataSourceName);
  TString name = Descriptor2Name(desc);
  TString oldRanges = Descriptor2Ranges(desc);
  
  TObjArray* a = oldRanges.Tokenize(fgkSeparatorWithinPart);
  
  TString newRanges;
  
  for ( Int_t i = 0; i <= a->GetLast(); ++i )
  {
    TString drange = static_cast<TObjString*>(a->At(i))->String();
    // drange should be of the form dimname/xmin/xmax
    
    TString bname = TupleFirst(drange,fgkSeparatorWithinRange);
    
    if ( bname == dimensionName )
    {
      newRanges += Form("%s%s%g%s%g",bname.Data(),fgkSeparatorWithinRange,
                        xmin,fgkSeparatorWithinRange,xmax);
    }
    else
    {
      newRanges += drange;
    }
    if ( i < a->GetLast() )
    {
      newRanges += fgkSeparatorWithinPart;
    }
  }

  delete a;

  AliInfo(Form("oldRanges=%s desc=%s newRanges=%s",oldRanges.Data(),desc.Data(),newRanges.Data()));
  
  AddDataSource(Form("%s%s%s%s%s",Descriptor2ID(desc).Data(),
                     fgkSeparatorBetweenDescriptorParts,
                     name.Data(),
                     fgkSeparatorBetweenDescriptorParts,
                     newRanges.Data()));
}

//_____________________________________________________________________________
TString
AliMUONPainterEnv::String(const char* resourceName, const char* defaultValue) const
{
  /// Retrieve the value associated with a given source, as a string
  
  return fEnv->GetValue(resourceName,defaultValue);
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::TupleFirst(const TString& tuple, const char* sep)
{
  /// Assuming tuple = aaa;bbb;....;yyy;zzz returns aaa
  Ssiz_t index = tuple.First(sep[0]);
  return tuple(0,index);
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::TupleSecond(const TString& tuple, const char* sep)
{
  /// Assuming tuple = aaa;bbb;....;yyy;zzz returns bbb
  return TupleFirst(TupleMiddle(tuple,sep),sep);
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::TupleMiddle(const TString& tuple, const char* sep)
{
  /// Assuming tuple = aaa;bbb;....;yyy;zzz returns bbb;....;yyy
  Ssiz_t i1 = tuple.First(sep[0]);
  Ssiz_t i2 = tuple.Last(sep[0]);
  return tuple(i1+1,i2-1);
}

//_____________________________________________________________________________
TString AliMUONPainterEnv::TupleLast(const TString& tuple, const char* sep)
{
  /// Assuming tuple = aaa;bbb;....;yyy;zzz returns zzz
  Ssiz_t index = tuple.Last(sep[0]);
  return tuple(index+1,tuple.Length()-index-1);
}

