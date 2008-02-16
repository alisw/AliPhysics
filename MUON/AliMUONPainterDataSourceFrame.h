#ifndef ALIMUONPAINTERDATASOURCEFRAME_H
#define ALIMUONPAINTERDATASOURCEFRAME_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterDataSourceFrame
/// \brief Frame to select input data source(s) to be displayed
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TGFrame
#  include <TGFrame.h>
#endif

class AliMUONPainterDataSourceItem;
class AliMUONVTrackerDataMaker;
class TGCheckButton;
class TGComboBox;
class TGGroupFrame;
class TGNumberEntry;
class TGTextEntry;
class TObjArray;

class AliMUONPainterDataSourceFrame : public TGCompositeFrame
{
public:
  AliMUONPainterDataSourceFrame(const TGWindow* p, UInt_t w, UInt_t h);
  virtual ~AliMUONPainterDataSourceFrame();
  
  void CalibrateButtonClicked(); 
  
  void CreateOCDBDataSource();

  void CreateRawDataSource();
  
  void DataReaderWasRegistered(AliMUONVTrackerDataMaker* reader);
  
  void DataReaderWasUnregistered(AliMUONVTrackerDataMaker* reader);

  void HistogramButtonClicked();
  
  void OpenFileDialog();
  
  void OpenRecentSource();

  void StartRunning(AliMUONPainterDataSourceItem* item);

  void StopRunning(AliMUONPainterDataSourceItem* item);

private:
  /// Not implemented
  AliMUONPainterDataSourceFrame(const AliMUONPainterDataSourceFrame& rhs);
  /// Not implemented
  AliMUONPainterDataSourceFrame& operator=(const AliMUONPainterDataSourceFrame& rhs);

  void AddRecentSource(const char* name);

  Bool_t CreateRawDataSource(const TString& uri);
  
  void CreateOCDBDataSource(const TString& uri);

  void CreateOCDBDataSource(const TString& cdbPath, Int_t runNumber, const TString& type);

private:
    
  TGGroupFrame* fRecentSourceSelector; ///< to select recently used sources   
  
  TGGroupFrame* fRawSelector; ///< to select a new raw data source
  TGCompositeFrame* fRawSelector2; ///< idem
  TGCompositeFrame* fRawSelector21; ///< idem
  TGCompositeFrame* fRawSelector22; ///< idem
  TGCompositeFrame* fRawSelector23; ///< idem
  TGCheckButton* fCalibrateNoGain; ///< to trig calibration of raw data (only 0 suppression)
  TGCheckButton* fCalibrateGainConstantCapa; ///< to trig calibration of raw data (0-supp and gain w/ constant capacitance)
  TGCheckButton* fCalibrateGain; ///< to trig calibration of raw data (full blown calibration)
  TGCheckButton* fHistogramButton; ///< to trig histogramming of raw data
  TGNumberEntry* fHistoMin; ///< xmin of histo to make
  TGNumberEntry* fHistoMax; ///< xmax of histo to make
  TGTextEntry* fRawOCDBPath; ///< OCDB path for raw data calibration
  
  TGGroupFrame* fOCDBSelector; ///< to select a new OCDB data source
  TGGroupFrame* fDataReaders; ///< to display currently active data sources  
  TGTextEntry* fFilePath; ///< raw data file path text entry widget
  TGTextEntry* fOCDBPath; ///< OCDB path text entry widget
  TGNumberEntry* fRunSelector; ///< OCDB run number entry widget
  TGComboBox* fOCDBTypes; ///< OCDB type combo box entry widget  
  TGComboBox* fRecentSources; ///< recent sources combo box  
  TObjArray* fItems; ///< list of data readers we handle
  
  static const char* fgkNumberOfDataSourcesKey; ///< key used to store the # of data sources in the resource file
  static const char* fgkDataSourceURIKey; ///< key usde to store the data source URIs in the resource file

  ClassDef(AliMUONPainterDataSourceFrame,2) // Data source selection frame
};

#endif
