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

/// \class AliMUONPainterPlotSelector
/// 
/// Widget to select which data to plot for painters
/// 
/// \author Laurent Aphecetche
///
/// See AliMUONPainterInterfaceHelper for an important implementation note
/// about our use of TGButtonGroup
///

#include "AliMUONPainterPlotSelector.h"

#include "AliMUONPainterGroup.h"
#include "AliMUONPainterInterfaceHelper.h"
#include "AliMUONPainterMatrix.h"
#include "AliMUONPainterDataRegistry.h"
#include "AliMUONVPainter.h"
#include "AliMUONVTrackerData.h"
#include "AliLog.h"
#include <Riostream.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TObjArray.h>
#include <TObjString.h>

///\cond CLASSIMP
ClassImp(AliMUONPainterPlotSelector)
///\endcond

const char* AliMUONPainterPlotSelector::fgkDefaultSourceName = "none";  

//_____________________________________________________________________________
AliMUONPainterPlotSelector::AliMUONPainterPlotSelector(const TGWindow* window, UInt_t w, UInt_t h)
: TGCompositeFrame(window,w,h,kHorizontalFrame),
fTypes(0x0),
fDataSourceNames(0x0),
fDataSourceDimensions(0x0),
fDimensionButtonMap(new TMap),
fCurrentType(""),
fCurrentData(0x0),
fCurrentDimension(-1)
{
  /// ctor
  fTypes = new TGButtonGroup(this,"Plot");

  fDataSourceNames = new TGButtonGroup(this,"Sources");
      
  AliMUONPainterDataRegistry* reg = AliMUONPainterDataRegistry::Instance();
  
  reg->Connect("DataSourceWasRegistered(AliMUONVTrackerData*)",
               "AliMUONPainterPlotSelector",
               this,
               "DataSourceWasRegistered(AliMUONVTrackerData*)");
  
  reg->Connect("DataSourceWasUnregistered(AliMUONVTrackerData*)",
               "AliMUONPainterPlotSelector",
               this,
               "DataSourceWasUnregistered(AliMUONVTrackerData*)");
    
  AliMUONPainterInterfaceHelper::AddRadioButton(*fDataSourceNames,
                                                fgkDefaultSourceName,
                                                0x0,
                                                kTRUE);
  
  CreateDimensionButtons(fgkDefaultSourceName);
  
  fDataSourceDimensions = new TGButtonGroup(this,0,3,5,0,"Dimensions");
  
  for ( Int_t i = 0; i < reg->NumberOfDataSources(); ++i )
  {
    AliMUONVTrackerData* data = reg->DataSource(i);
    DataSourceWasRegistered(data);
  }
  
  fDataSourceNames->Connect("Clicked(Int_t)","AliMUONPainterPlotSelector",
                            this,
                            "SourceButtonWasClicked(Int_t)");

  AddFrame(fTypes);
  AddFrame(fDataSourceNames);
  AddFrame(fDataSourceDimensions);
}

//_____________________________________________________________________________
AliMUONPainterPlotSelector::~AliMUONPainterPlotSelector()
{
  /// dtor
}

//_____________________________________________________________________________
void
AliMUONPainterPlotSelector::BackupDimensionButtons()
{
  /// Backup the dimension button group
  
  TString name = fDataSourceDimensions->GetTitle();

  if ( name !=  fgkDefaultSourceName )
  {
    TGButtonGroup* group = static_cast<TGButtonGroup*>(fDimensionButtonMap->GetValue(name));
    if (!group) 
    {
      AliError(Form("Did not find group %s",name.Data()));
    }
    else
    {
      AliMUONPainterInterfaceHelper::Copy(*fDataSourceDimensions,*group);
    }
  
  }
  
  fDataSourceDimensions->Disconnect("Clicked(Int_t)",
                                    this,
                                    "DimensionButtonWasClicked(Int_t)");  
}


//_____________________________________________________________________________
void
AliMUONPainterPlotSelector::CreateDimensionButtons(const char* dataSourceName)
{
  /// Create the dimension button group for a given data source
  
  AliMUONVTrackerData* data = AliMUONPainterDataRegistry::Instance()->DataSource(dataSourceName);

  TGButtonGroup* bg = new TGButtonGroup(this,0,3,5,0,dataSourceName);
  
  if ( data ) 
  {
    for ( Int_t i = 0; i < data->NumberOfDimensions(); ++i ) 
    {
      AliMUONPainterInterfaceHelper::AddRadioButton(*bg,
                                                    data->DimensionName(i),
                                                    reinterpret_cast<void*>(i));
    }
  }
  
  fDimensionButtonMap->Add(new TObjString(dataSourceName),bg);
  
  bg->Connect("Clicked(Int_t)","AliMUONPainterPlotSelector",this,
                                "DimensionButtonWasClicked(Int_t)");
}

//_____________________________________________________________________________
void
AliMUONPainterPlotSelector::CreateTypeButtons(const TObjArray& types)
{
  /// Create the type button group
  
  AliMUONPainterInterfaceHelper::ClearButtons(*fTypes);

  TIter nextType(&types);
  TObjString* str;

  while ( ( str = static_cast<TObjString*>(nextType()) ) )
  {
    AliMUONPainterInterfaceHelper::AddRadioButton(*fTypes,str->String());
  }

  fTypes->Connect("Clicked(Int_t)","AliMUONPainterPlotSelector",this,
                  "TypeButtonWasClicked(Int_t)");
  
  fTypes->Show();
}

//_____________________________________________________________________________
void
AliMUONPainterPlotSelector::DataSourceWasChanged()
{
  /// Data source was changed
  DataSourceWasChanged(fCurrentType.Data(),fCurrentData,fCurrentDimension);
}

//_____________________________________________________________________________
void 
AliMUONPainterPlotSelector::DataSourceWasChanged(const char* type, 
                                                 AliMUONVTrackerData* data,
                                                 Int_t dataIndex)
{
  /// Emit a signal to tell data source was changed
  
	UpdateTypeButton();
	
  Long_t param[] = { (Long_t)type, (Long_t)data,
    (Long_t)dataIndex };
  
  Emit("DataSourceWasChanged(const char*,AliMUONVTrackerData*,Int_t)",param);
}

//_____________________________________________________________________________
void
AliMUONPainterPlotSelector::DataSourceWasRegistered(AliMUONVTrackerData* data)
{
  /// A new data source has been registered : add it to the interface
  
  AliMUONPainterInterfaceHelper::AddRadioButton(*fDataSourceNames,
                                                data->GetName(),
                                                data);
  
  data->Connect("NumberOfEventsChanged()",
                "AliMUONPainterPlotSelector",
                this,
                "NumberOfEventsChanged()");
  
  CreateDimensionButtons(data->GetName());
  
  fDataSourceNames->Show();
  
  Layout();
}

//_____________________________________________________________________________
void
AliMUONPainterPlotSelector::NumberOfEventsChanged()
{
  /// Change the tool tip of the corresponding data source button

  // find first the sender of the signal
  
//  AliMUONVTrackerData* data = reinterpret_cast<AliMUONVTrackerData*>(gTQSender);
//  
//  TGButton* button = AliMUONPainterInterfaceHelper::FindButtonByUserData(*fDataSourceNames,data);
//  
//  if (button)
//  {
//    button->SetToolTipText(Form("%d events",data->NumberOfEvents()),250);
//  }
}

//_____________________________________________________________________________
void
AliMUONPainterPlotSelector::DataSourceWasUnregistered(AliMUONVTrackerData* data)
{
  /// A data source has been unregistered : remove it from the interface
  
  TGButton* button = AliMUONPainterInterfaceHelper::FindButtonByUserData(*fDataSourceNames,data);

  TGButton* bd = AliMUONPainterInterfaceHelper::FindDownButton(*fDataSourceNames);

  if ( bd == button ) 
  {
    // selected data source is the one we are removing...
    // revert to "none" before actually removing it.
    SourceButtonWasClicked(1);
  }
  
  AliMUONPainterInterfaceHelper::RemoveButton(*fDataSourceNames,button);

  // do not forget to re-connect things
  fDataSourceNames->Connect("Clicked(Int_t)","AliMUONPainterPlotSelector",
                            this,
                            "SourceButtonWasClicked(Int_t)");
  

  TObject* o = fDimensionButtonMap->Remove(new TObjString(data->GetName()));
  
  if (!o)
  {
    AliError("Remove failed. Please check");    
  }

  fDataSourceNames->Show();
  
  Layout();
}

//_____________________________________________________________________________
void 
AliMUONPainterPlotSelector::DimensionButtonWasClicked(Int_t id)
{
  /// One dim button was clicked
  
  TGTextButton* button = (TGTextButton*)fDataSourceDimensions->GetButton(id);
  
  SetCurrentDimension(reinterpret_cast<Long_t>(button->GetUserData()));
  
  if ( fCurrentDimension >= 0 )
  {
    AliMUONPainterInterfaceHelper::SetState(*fTypes,kTRUE);
    AliMUONPainterInterfaceHelper::Select(*fTypes,fCurrentType.Data(),kFALSE);
    fTypes->Show();
  }
  
  DataSourceWasChanged();
}

//_____________________________________________________________________________
void
AliMUONPainterPlotSelector::ResetDimensionButtonMap()
{
  /// Reset the button group map
  
  TIter next(fDimensionButtonMap);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    TGButtonGroup* bg = static_cast<TGButtonGroup*>(fDimensionButtonMap->GetValue(str->String()));
    AliMUONPainterInterfaceHelper::Unselect(*bg,"*");
  }
}

//_____________________________________________________________________________
void
AliMUONPainterPlotSelector::RestoreDimensionButtons(const char* dataSourceName,
                                                    Bool_t updateCurrentDimension)
{
  /// Restore (i.e. contrary of Backup) a given dimension button group
  
  TGButtonGroup* group = static_cast<TGButtonGroup*>(fDimensionButtonMap->GetValue(dataSourceName));
  
  AliMUONPainterInterfaceHelper::Copy(*group,*fDataSourceDimensions);
  
  fDataSourceDimensions->Connect("Clicked(Int_t)",
                                 "AliMUONPainterPlotSelector",this,
                                 "DimensionButtonWasClicked(Int_t)");    
  
  if ( updateCurrentDimension ) 
  {
    TGButton* button = AliMUONPainterInterfaceHelper::FindDownButton(*fDataSourceDimensions);
    if ( button ) 
    {
      SetCurrentDimension(reinterpret_cast<Long_t>(button->GetUserData()));
    }
    else 
    {
      SetCurrentDimension(-1);
    }
  }
  
  fDataSourceDimensions->Show();
}

//_____________________________________________________________________________
void 
AliMUONPainterPlotSelector::SetCurrentData(AliMUONVTrackerData* data)
{
  /// Set the current data pointer
  fCurrentData = data;
}

//_____________________________________________________________________________
void 
AliMUONPainterPlotSelector::SetCurrentDimension(Long_t i)
{
  /// Set the current dimension
  fCurrentDimension = i;
}

//_____________________________________________________________________________
void 
AliMUONPainterPlotSelector::SetCurrentType(const char* type)
{
  /// Set the current type
  fCurrentType = type;
}

//_____________________________________________________________________________
void
AliMUONPainterPlotSelector::SourceButtonWasClicked(Int_t id)
{
  /// A source button was clicked

  BackupDimensionButtons();
  
  TGButton* button = fDataSourceNames->GetButton(id);
  if ( !button ) 
  {
    AliFatal(Form("Could not get DataSource button id=%d",id));
  }
  
  AliMUONVTrackerData* data = reinterpret_cast<AliMUONVTrackerData*>(button->GetUserData());

  TString name =  data ? data->GetName() : fgkDefaultSourceName;

  RestoreDimensionButtons(name,kTRUE);
  
  if ( data != 0 && 
       fCurrentDimension >= 0 && 
       fCurrentDimension < (Long_t)(data->NumberOfDimensions()) )
  {
    AliMUONPainterInterfaceHelper::SetState(*fTypes,kTRUE);
  }
  else
  {
    AliMUONPainterInterfaceHelper::SetState(*fTypes,kFALSE);
  }
  
  SetCurrentData(data);
  
  UpdateTypeButton();
  
  fDataSourceNames->Show();
  fDataSourceDimensions->Show();
  fTypes->Show();
  
  Resize();
  Layout();
  
  DataSourceWasChanged();
}

//_____________________________________________________________________________
void
AliMUONPainterPlotSelector::TypeButtonWasClicked(Int_t id)
{
  /// A type button was clicked
  TGTextButton* button = (TGTextButton*)fTypes->GetButton(id);
  SetCurrentType(button->GetTitle());
  DataSourceWasChanged();
}

//_____________________________________________________________________________
void AliMUONPainterPlotSelector::Update(const AliMUONPainterMatrix& painterMatrix)
{
  /// Update ourselves from a new painter matrix
  
  SetCurrentType("");
  SetCurrentData(0x0);
  SetCurrentDimension(-1);
  
  AliMUONPainterInterfaceHelper::Unselect(*fDataSourceNames,"*");
  AliMUONPainterInterfaceHelper::Unselect(*fDataSourceDimensions,"*");
  
  ResetDimensionButtonMap();
  
  TObjArray types;
  types.SetOwner(kTRUE);  
  painterMatrix.GetTypes(types);  
  CreateTypeButtons(types);
  
  if ( painterMatrix.Size() > 0 ) 
  {
    AliMUONVPainter* painter = painterMatrix.Painter(painterMatrix.Size()-1);
    
    AliMUONPainterGroup* plotterGroup = painter->PlotterGroup();
    
    if ( plotterGroup )
    {
      // the group have some data to plot, so update our internal pointers
      // to reflect that
      SetCurrentData(plotterGroup->Data());
      SetCurrentDimension(plotterGroup->DataIndex());
      SetCurrentType(plotterGroup->Type());
    }
  }
  
  // the *order* of the 3 following lines is *very* important

  UpdateSourceButton();
  UpdateDimensionButton();
  UpdateTypeButton();
  
  Resize();
  Layout();
}

//_____________________________________________________________________________
void 
AliMUONPainterPlotSelector::UpdateDimensionButton()
{
  /// Update the dim buttons
  
  TGTextButton* button = static_cast<TGTextButton*>
  (AliMUONPainterInterfaceHelper::FindButtonByUserData(*fDataSourceDimensions,
                                                       reinterpret_cast<void*>(fCurrentDimension)));
  
  if ( button ) 
  {
    // set this button as "ON"
    AliMUONPainterInterfaceHelper::Select(*fDataSourceDimensions,button->GetTitle());
  }
  else
  {
    AliMUONPainterInterfaceHelper::Unselect(*fDataSourceDimensions,"*");
  }
  
  fDataSourceDimensions->Show();
}

//_____________________________________________________________________________
void
AliMUONPainterPlotSelector::UpdateSourceButton()
{
  /// Update the source buttons
  
  TGTextButton* button = static_cast<TGTextButton*>
  (AliMUONPainterInterfaceHelper::FindButtonByUserData(*fDataSourceNames,
                                                       fCurrentData));
  
  if ( button ) 
  {
    AliMUONPainterInterfaceHelper::Select(*fDataSourceNames,button->GetTitle());
  
    RestoreDimensionButtons(button->GetTitle(),kFALSE);
  }
  
  fDataSourceNames->Show();
}

//_____________________________________________________________________________
void
AliMUONPainterPlotSelector::UpdateTypeButton()
{
  /// Update the type buttons
	
	if (!fCurrentData)
  {
    AliMUONPainterInterfaceHelper::SetState(*fTypes,kFALSE);
  }
	else
	{
		// disable levels that cannot be handled by this data
		TGTextButton* padButton = static_cast<TGTextButton*>
		(AliMUONPainterInterfaceHelper::FindButtonByName(*fTypes,"PAD"));
		if (padButton) 
		{ 
			padButton->SetEnabled(fCurrentData->IsChannelLevelEnabled());
		}
		TGTextButton* manuButton = static_cast<TGTextButton*>
		(AliMUONPainterInterfaceHelper::FindButtonByName(*fTypes,"MANU"));
		if (manuButton) 
		{ 
			manuButton->SetEnabled(fCurrentData->IsManuLevelEnabled());
		}
		TGTextButton* busPatchButton = static_cast<TGTextButton*>
		(AliMUONPainterInterfaceHelper::FindButtonByName(*fTypes,"BUSPATCH"));
		if (busPatchButton) 
		{ 
			busPatchButton->SetEnabled(fCurrentData->IsBusPatchLevelEnabled());
		}
		TGTextButton* pcbButton = static_cast<TGTextButton*>
		(AliMUONPainterInterfaceHelper::FindButtonByName(*fTypes,"PCB"));
		if (pcbButton) 
		{ 
			pcbButton->SetEnabled(fCurrentData->IsPCBLevelEnabled());
		}
    
	}
	
  TGTextButton* button = static_cast<TGTextButton*>
  (AliMUONPainterInterfaceHelper::FindButtonByName(*fTypes,fCurrentType));

  if ( button ) 
  {
    AliMUONPainterInterfaceHelper::Select(*fTypes,button->GetTitle());
  }
  else
  {
    AliMUONPainterInterfaceHelper::Unselect(*fTypes,"*");
  }

	
  fTypes->Show();
}

