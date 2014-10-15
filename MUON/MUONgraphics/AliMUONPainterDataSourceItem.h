#ifndef ALIMUONPAINTERDATASOURCEITEM_H
#define ALIMUONPAINTERDATASOURCEITEM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterDataSourceItem
/// \brief A widget describing a single data source
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TGFrame
#  include "TGFrame.h"
#endif

class AliMUONVTrackerDataMaker;
class TGLabel;
class TGTextButton;
class AliMUONPainterDataReader;
class TThread;

class AliMUONPainterDataSourceItem : public TGCompositeFrame
{
public:
  AliMUONPainterDataSourceItem(const TGWindow* p, UInt_t w, UInt_t h, 
                               AliMUONVTrackerDataMaker* maker);
  virtual ~AliMUONPainterDataSourceItem();
  
  void Run();
  
  void Stop(); 
  
  void Rewind();
  
  void Remove();
  
  void Update();
  
  void Reset();
  
  /// Return data source reader
  AliMUONVTrackerDataMaker* DataMaker() const { return fDataMaker; }
  
  void EnableRun();
  
  void DisableRun();
  
  void StartRunning(); //*SIGNAL*

  void StopRunning(); //*SIGNAL*

  void Save();
  
  void SaveWithDialog();
  
private:
    
    void Save(const char* filename);
  
  /// Not implemented
  AliMUONPainterDataSourceItem(const AliMUONPainterDataSourceItem& rhs);
  /// Not implemented
  AliMUONPainterDataSourceItem& operator=(const AliMUONPainterDataSourceItem& rhs);
  
  AliMUONVTrackerDataMaker* fDataMaker; ///< data source reader (not owner)  
  TGLabel* fSourceName; ///< the (short) name of the data source
  TGLabel* fSource; ///< the full uri of the data source
  TGLabel* fNumberOfEvents; ///< number of evts this source has seen so far
  TGTextButton* fRun; ///< button to start running over the source
  TGTextButton* fStop; ///< button to stop running over the source
  TGTextButton* fRewind; ///< button to rewind events for the source
  TGTextButton* fRemove; ///< button to remove the source
  TGTextButton* fSave; ///< button to save the source (filename is fixed)
  TGTextButton* fSaveAs; ///< button to save as... 
  TThread* fThread; ///< thread used to actually loop over the data
  Long_t fParams[2]; ///< used in conjunction with fThread
  
  Bool_t fShouldReset; ///< whether we should reset or not...
  
  ClassDef(AliMUONPainterDataSourceItem,3) // Data source widget
};

#endif
