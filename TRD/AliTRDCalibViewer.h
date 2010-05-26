#ifndef ALITRDCALIBVIEWER_H
#define ALITRDCALIBVIEWER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalibViewer.h 34418 2009-08-26 15:47:50Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Class which implements AliBaseCalibViewer for the TRD                    //
//  used for the calibration monitor                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include "TFriendElement.h"
#include "TVectorD.h"
#include "AliBaseCalibViewer.h"
class AliTRDCalDet;
class AliTRDCalPad;


class AliTRDCalibViewer : public AliBaseCalibViewer {
 public:
  AliTRDCalibViewer();
  AliTRDCalibViewer(const AliTRDCalibViewer &c);
  AliTRDCalibViewer(TTree* tree);
  AliTRDCalibViewer(const char* fileName, const char* treeName = "TRDcalibDetails");
  AliTRDCalibViewer &operator= (const AliTRDCalibViewer& param);
  virtual ~AliTRDCalibViewer();

  virtual TObjArray* GetListOfVariables(Bool_t printList = kFALSE);
  virtual TObjArray* GetListOfNormalizationVariables(Bool_t printList = kFALSE) const;

  //virtual void GetTimeInfoOCDB(const Char_t* runList, const Char_t* outFile,
                               //Int_t firstRun, Int_t lastRun, UInt_t infoFlags,
                               //const Char_t* ocdbStorage);
  // extract pad level OCDB information for a run list and dump it into a tree
  Bool_t DumpOCDBtoTreeDetails(const Char_t* runListFilename, const Char_t* outFilename,
			       Int_t firstRun, Int_t lastRun, const Char_t* storage,
                               Int_t version = -1, Int_t subVersion = -1,
                               Bool_t getCalibs = kTRUE, Bool_t getDCS = kTRUE, Bool_t getAlign = kTRUE);
  // read AliTRDCalPad objects from a root file and dump a root tree
  void DumpCalibToTree(const Char_t* inFilename, const Char_t* outFilename);
  // extract averages from calibration objects
  void ProcessTRDCalibArray(AliTRDCalDet* chamberCalib, AliTRDCalPad *padCalib,
			    TString parName, 
			    Double_t &runValue, Double_t &runRMS,
			    TVectorD &chamberValues, TVectorD &chamberValuesRMS,
			    TVectorD &superModuleValues, TVectorD &superModuleValuesRMS);
  // extract averages from calibration objects
  void ProcessTRDCalibArray(AliTRDCalPad *padCalib,
                            TVectorD &superModuleValues, TVectorD &superModuleValuesRMS);

  virtual const char* AddAbbreviations(char* c, Bool_t printDrawCommand = kFALSE);
  void GetLayerSectorStack(TString trdString, Int_t& layerNo, Int_t& sectorNo, Int_t& stackNo) const;
  // easy drawing of data, use '~' for abbreviation of '.fElements'
  virtual Int_t EasyDraw(const char* drawCommand, const char* sector, const char* cuts = 0, 
			 const char* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const;   
  // easy drawing of data, use '~' for abbreviation of '.fElements'
  virtual Int_t EasyDraw(const char* drawCommand, Int_t sector, const char* cuts = 0, 
			 const char* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const;   
  // easy drawing of data, use '~' for abbreviation of '.fElements'
  virtual Int_t EasyDraw1D(const char* drawCommand, const char* sector, const char* cuts = 0, 
			   const char* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const;   
  // easy drawing of data, use '~' for abbreviation of '.fElements'
  virtual Int_t EasyDraw1D(const char* drawCommand, Int_t sector, const char* cuts = 0, 
			   const char* drawOptions = 0, Bool_t writeDrawCommand = kFALSE) const;  

  ClassDef(AliTRDCalibViewer,1)    //  TRD calibration viewer class
};

#endif
