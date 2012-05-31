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

class AliMUONAttPainter;
class AliMUONPainterDataSourceItem;
class AliMUONPainterMatrix;
class AliMUONVTrackerDataMaker;
class AliMUONVTrackerData;
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

  void CreateACFDataSource();

  void CreateRawDataSource();
  
  void DataMakerWasRegistered(AliMUONVTrackerDataMaker* reader);
  
  void DataMakerWasUnregistered(const AliMUONVTrackerDataMaker* reader);

  void HistogramButtonClicked();

  void EventRangeButtonClicked();

  void OpenFileDialog();
  
  void OpenFileDialogACF();
  
  void OpenRecentSource();

  void StartRunning();

  void StopRunning();

  static void CreatePedestalCanvases(AliMUONVTrackerData* data,
                                     Double_t pedMin=0, Double_t pedMax=500,
                                     Double_t sigmaMin=0, Double_t sigmaMax=5);
  
  static AliMUONPainterMatrix* CreateFullTracker(AliMUONVTrackerData* data, 
                                                 Int_t dim, 
                                                 Double_t xmin, Double_t xmax, 
                                                 const AliMUONAttPainter& att);

private:
  /// Not implemented
  AliMUONPainterDataSourceFrame(const AliMUONPainterDataSourceFrame& rhs);
  /// Not implemented
  AliMUONPainterDataSourceFrame& operator=(const AliMUONPainterDataSourceFrame& rhs);

  void AddRecentSource(const char* name);

  Bool_t CreateRawDataSource(const TString& uri);
  
  void CreateOCDBDataSource(const TString& uri);

  void CreateOCDBDataSource(const TString& cdbPath, Int_t runNumber, const TString& type);

  void CreateACFDataSource(const TString& uri);

  void CreateACFDataSource(const TString& acfPath, const TString& type);
  
  void RegisterDataSource(AliMUONVTrackerDataMaker* reader, const char* dsName);
  
private:
    
  TGGroupFrame* fRecentSourceSelector; ///< to select recently used sources   
  
  TGGroupFrame* fRawSelector; ///< to select a new raw data source
  TGCompositeFrame* fRawSelector2; ///< idem
  TGCompositeFrame* fRawSelector21; ///< idem
  TGCompositeFrame* fRawSelector22; ///< idem
  TGCompositeFrame* fRawSelector24; ///< idem
  TGCompositeFrame* fRawSelector23; ///< idem
  TGCheckButton* fCalibrateNoGain; ///< to trig calibration of raw data (only 0 suppression)
  TGCheckButton* fCalibrateGainConstantCapa; ///< to trig calibration of raw data (0-supp and gain w/ constant capacitance)
  TGCheckButton* fCalibrateGain; ///< to trig calibration of raw data (full blown calibration)
  TGCheckButton* fCalibrateEmelecGain; ///< to trig calibration of raw data (full blown calibration but with factory gains)
  TGCheckButton* fHistogramButton; ///< to trig histogramming of raw data  
  TGNumberEntry* fHistoMin; ///< xmin of histo to make
  TGNumberEntry* fHistoMax; ///< xmax of histo to make
  TGCheckButton* fEventRangeButton; ///< to trig limitation of event range
  TGNumberEntry* fEventMin; ///< min event number to consider
  TGNumberEntry* fEventMax; ///< max event number to consider  
  TGTextEntry* fRawOCDBPath; ///< OCDB path for raw data calibration
  
  TGGroupFrame* fOCDBSelector; ///< to select a new OCDB data source
  TGGroupFrame* fDataReaders; ///< to display currently active data sources  
  TGTextEntry* fFilePath; ///< raw data file path text entry widget
  TGTextEntry* fOCDBPath; ///< OCDB path text entry widget
  TGNumberEntry* fRunSelector; ///< OCDB run number entry widget
  TGComboBox* fOCDBTypes; ///< OCDB type combo box entry widget  
  TGComboBox* fRecentSources; ///< recent sources combo box  
  TObjArray* fItems; ///< list of data readers we handle
  
  TGGroupFrame* fACFSelector; ///< to select ACF (ASCII calibration files)
  TGTextEntry* fACFPath; ///< path to ASCII calibration file
  TGComboBox* fACFTypes; ///< types of ASCII calibration files 

  static const char* fgkNumberOfDataSourcesKey; ///< key used to store the # of data sources in the resource file
  static const char* fgkDataSourceURIKey; ///< key usde to store the data source URIs in the resource file

  ClassDef(AliMUONPainterDataSourceFrame,5) // Data source selection frame
};

#endif
