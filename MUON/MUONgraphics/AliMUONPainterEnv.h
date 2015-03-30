#ifndef ALIMUONPAINTERENV_H
#define ALIMUONPAINTERENV_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterEnv
/// \brief Resource file handling
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class TEnv;

class AliMUONPainterEnv : public TObject
{
public:
  AliMUONPainterEnv(const char* resourceFile=".mchviewrc");
  virtual ~AliMUONPainterEnv();
  
  TString AddDataSource(const char* dataSourceDescriptor);

  TString DataSourceDescriptor(const char* dataSourceName) const;

  TString DataSourceDescriptor(Int_t index) const;

  TString DataSourceName(const char* dataSourceName) const;
  TString DataSourceID(const char* dataSourceName) const;

  TString DataSourceOrigin(const char* dataSourceName) const;
  TString DataSourceRanges(const char* dataSourceName) const;
  TString DataSourceType(const char* dataSourceName) const;
  TString DataSourceURI(const char* dataSourceName) const;
  
  static TString Descriptor2ID(const char* dataSourceDescriptor);
  static TString Descriptor2Name(const char* dataSourceDescriptor);
  static TString Descriptor2Ranges(const char* dataSourceDescriptor);
  
  Double_t Double(const char* resourceName, Double_t defaultValue=0.0) const;
  
  void ForceDataSourceToDefaultRange(const char* dataSourceName, const char* dim="");
  
  Int_t GetDataSourceIndex(const char* dataSourceName) const;

  TString GetDefaultRange(const char* dataSourceType, const char* dimensionName) const;
  
  static TString ID2Origin(const char* dataSourceID);
  static TString ID2URI(const char* dataSourceID);
  static TString ID2Type(const char* dataSourceID);
  
  Int_t Integer(const char* resourceName, Int_t defaultValue=0) const;
  
  Int_t NumberOfDataSources() const;

  void Print(Option_t* opt="") const;
  
  static Bool_t Ranges2DimensionRange(const char* ranges, const char* dimensionName, Double_t& xmin, Double_t& xmax);
  
  void Save();
  
  void Set(const char* resourceName, Int_t value);

  void Set(const char* resourceName, const char* value);

  void Set(const char* resourceName, Double_t value);

  void SetDefaultRange(const char* dataSourceType, const char* dimensionName, Double_t xmin, Double_t xmax);

  void SetDimensionRange(const char* dataSourceName, const char* dimensionName, Double_t xmin, Double_t xmax);

  TString String(const char* resourceName, const char* defaultValue="") const;
  
  static TString TupleFirst(const TString& tuple, const char* sep);
  static TString TupleLast(const TString& tuple, const char* sep);
  static TString TupleMiddle(const TString& tuple, const char* sep);
  static TString TupleSecond(const TString& tuple, const char* sep);

  static const char* SeparatorWithinPart() { return fgkSeparatorWithinPart; }

  static const char* SeparatorBetweenDescriptorParts() { return fgkSeparatorBetweenDescriptorParts; }

private:
  /// Not implemented
  AliMUONPainterEnv(const AliMUONPainterEnv& rhs);
  /// Not implemented
  AliMUONPainterEnv& operator=(const AliMUONPainterEnv& rhs);
  
  TEnv* fEnv; ///< the worker class
  
  static const char* fgkNumberOfDataSourcesKey; ///< key used to store the # of data sources in the resource file
  static const char* fgkDataSourceKey; ///< key used to store the data source URIs in the resource file
  static const char* fgkDefaultRangeKey; ///< key used to specify the default ranges for data representation
  static const char* fgkDisableAutoPedCanvasKey;
  static const char* fgkSeparatorWithinRange;
  static const char* fgkSeparatorWithinPart;
  static const char* fgkSeparatorBetweenDescriptorParts;

  ClassDef(AliMUONPainterEnv,3) // Painter display resource file
};

#endif
