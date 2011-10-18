#ifndef ALITOFQADATAMAKERREC_H
#define ALITOFQADATAMAKERREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////
//                                                                // 
//  Produces the data needed to calculate the quality assurance.  //
//    All data must be mergeable objects.                         //
//    S. Arcelli                                                  //
//                                                                //
//    /* last modified by F. Bellini on 25/02/2010 */             // 
////////////////////////////////////////////////////////////////////

#include <TLine.h>
#include "AliQADataMakerRec.h"
#include "AliQAv1.h"
#include "AliTOFcalib.h"
#include "AliTOFTrigger.h"
class AliCDBManager;
class AliCDBEntry;
class AliCDBStorage;
class AliTOFChannelOnlineStatusArray;
class AliTOFDecoderSummaryData;

class AliTOFQADataMakerRec: public AliQADataMakerRec {

public:
  AliTOFQADataMakerRec() ;          // ctor
  AliTOFQADataMakerRec(const AliTOFQADataMakerRec& qadm) ;   
  AliTOFQADataMakerRec& operator = (const AliTOFQADataMakerRec& qadm) ;
  AliTOFChannelOnlineStatusArray *GetCalibData() ;
  virtual ~AliTOFQADataMakerRec(); // dtor
  
  void   GetGeo2LTMIndex(const Int_t * const detind, Int_t *indexLTM);
  void   GetGeo2CTTMIndex(Int_t *detind, Int_t *indexCTTM);
  void   GetCTTMIndex(Int_t *equipid, Int_t *indexCTTM);
  void ReadHistogramRangeFromFile(const Char_t * filename);
  void SetDefaultHistogramRange();
  void SetDefaultMultiHistogramRange();
  void SetDefaultTimeHistogramRange();
  void SetDefaultCutNmaxFiredMacropad();

protected: 
  AliTOFChannelOnlineStatusArray * fCalibData;        //! calibration data
  
private:
  virtual void   InitESDs() ; 
  virtual void   InitRecPoints() ; 
  virtual void   InitRaws() ; 
  virtual void   MakeESDs(AliESDEvent * const esd) ;
  virtual void   MakeRecPoints(TTree * recTree) ; 
  virtual void   MakeRaws(AliRawReader* rawReader) ; 
  virtual void   StartOfDetectorCycle() ; 
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list) ;
  virtual void   GetMapIndeces(const Int_t * const in, Int_t *out) ; 
  	  Int_t  GetStripIndex(const Int_t * const in);
  virtual void   EnableNoiseFiltering(Bool_t enable){fEnableNoiseFiltering = enable;};
  virtual void   EnableDqmShifterOpt(Bool_t enable){ fEnableDqmShifterOpt = enable;};
          Bool_t CheckVolumeID(const Int_t * const equipmentID); 
          Bool_t CheckEquipID( const Int_t * const equipmentID); 
          Bool_t FilterLTMData(const Int_t * const equipmentID) const ; 
          Bool_t FilterSpare(  const Int_t * const equipmentID) const ;
	  Int_t  GetNumberOfFiredMacropad(AliRawReader * rawReader);
  static void SetNbinsMultiplicityHisto(Int_t value){fgNbinsMultiplicity=value; return;}; 
  static void SetMultiplicityHistoRange (Int_t valueMin, Int_t valueMax){fgRangeMinMultiplicity=valueMin; fgRangeMaxMultiplicity=valueMax; return;}
  static void SetNbinsTimeHisto(Int_t value){fgNbinsTime=value; return;};
  static void SetTimeHistoRange (Float_t valueMin, Float_t valueMax){fgRangeMinTime=valueMin; fgRangeMaxTime=valueMax; return;};
  static void SetCutNmaxFiredMacropad(Int_t value){fgCutNmaxFiredMacropad=value;return;};

 	  // void   ResetAllTRMcounters();
	  Bool_t fEnableNoiseFiltering; //the choice is not implemented so far
	  Bool_t fEnableDqmShifterOpt;  // draw option flag to help
					// DQM shifter in the
					// interpretation of the TOF
					// raw QA histograms
	  Bool_t fIsSOC;  //flag for StartOfCycle operations
	  //lines for the DQM GUI
	  TLine* fLineExpTimeMin;
	  TLine* fLineExpTimeMax;
	  TLine* fLineExpTotMin;
	  TLine* fLineExpTotMax;
	  TLine* fLineSMid[17];
	  TLine* fLineLTMid[71];
	  TLine* fLineLTMbitId[22];

	  AliTOFRawStream fTOFRawStream; // AliTOFRawStream variable 
	  AliTOFDecoderSummaryData * fDecoderSummary; //pointer to decoder summary data object
	  Int_t fRunNumber; //run number
	  static Int_t fgNbinsMultiplicity;//number of bins in multiplicity plot
	  static Int_t fgRangeMinMultiplicity;//min range in multiplicity plot
	  static Int_t fgRangeMaxMultiplicity;//max range in multiplicity plot
	  static Int_t fgNbinsTime;//number of bins in time plot
	  static const Float_t fgkNbinsWidthTime;//width of bins in time plot
	  static Float_t fgRangeMinTime;//range min in time plot
	  static Float_t fgRangeMaxTime; //range max in time plot
	  static Int_t fgCutNmaxFiredMacropad; //cut on max number of fired macropad 
	  static const Int_t fgkFiredMacropadLimit; //limit on cut on number of fired macropad 
	  AliTOFcalib fCalib;//calibration object

	  ClassDef(AliTOFQADataMakerRec,8)  // description 	    
};

#endif // ALITOFQADATAMAKERREC_H
