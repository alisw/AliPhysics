#ifndef ALIMUONPAINTERPLOTSELECTOR_H
#define ALIMUONPAINTERPLOTSELECTOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterPlotSelector
/// \brief Widget to pick what to plot for the painters
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TGFrame
#  include "TGFrame.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONPainterMatrix;
class AliMUONVTrackerData;
class TGButtonGroup;
class TMap;

class AliMUONPainterPlotSelector : public TGCompositeFrame
{
public:
  AliMUONPainterPlotSelector(const TGWindow* window, UInt_t w=1, UInt_t h=1);
  virtual ~AliMUONPainterPlotSelector();
  
  void DataSourceWasRegistered(AliMUONVTrackerData* data);  
  
  void DataSourceWasUnregistered(AliMUONVTrackerData* data);

  void DataSourceWasChanged(const char* type, 
                            AliMUONVTrackerData* data,
                            Int_t dataIndex); // *SIGNAL*
  
  void DimensionButtonWasClicked(Int_t id);
  
  void SourceButtonWasClicked(Int_t id);
  
  void TypeButtonWasClicked(Int_t id);
    
  void Update(const AliMUONPainterMatrix& painterMatrix);

  void NumberOfEventsChanged();

private:
  /// Not implemented
  AliMUONPainterPlotSelector(const AliMUONPainterPlotSelector& rhs);
  /// Not implemented
  AliMUONPainterPlotSelector& operator=(const AliMUONPainterPlotSelector& rhs);
  
  void BackupDimensionButtons();
  
  void CreateDimensionButtons(const char* dataSourceName);

  void CreateTypeButtons(const TObjArray& types);
  
  void DataSourceWasChanged();

  void ResetDimensionButtonMap();

  void RestoreDimensionButtons(const char* dataSourceName,
                               Bool_t updateCurrentDimension);
  
  void SetCurrentData(AliMUONVTrackerData* data);
  
  void SetCurrentDimension(Long_t i);
  
  void SetCurrentType(const char* type);
  
  void UpdateDimensionButton();
  
  void UpdateSourceButton();
  
  void UpdateTypeButton();
  
private:
  
  TGButtonGroup* fTypes; ///< types buttons
  TGButtonGroup* fDataSourceNames; ///< data source names buttons
  TGButtonGroup* fDataSourceDimensions; ///< data source dimensions buttons  
  TMap* fDimensionButtonMap; ///< cache for button group  
  TString fCurrentType; ///< current type
  AliMUONVTrackerData* fCurrentData; ///< current data
  Long_t fCurrentDimension; ///< current data index
  static const char* fgkDefaultSourceName; ///< default source name
  
  ClassDef(AliMUONPainterPlotSelector,1) // Widget to select what to plot for painters
};

#endif
