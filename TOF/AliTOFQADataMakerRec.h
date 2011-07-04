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
  AliTOFChannelOnlineStatusArray *GetCalibData() const;
  virtual ~AliTOFQADataMakerRec(); // dtor
  
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
	  // void   ResetAllTRMcounters();
	  Bool_t fEnableNoiseFiltering; //the choice is not implemented so far
	  Bool_t fEnableDqmShifterOpt;  // draw option flag to help
					// DQM shifter in the
					// interpretation of the TOF
					// raw QA histograms
	
	  //	  Int_t  fProcessedRawEventN;   // number of processed rawData events
	  Bool_t fIsSOC;  //flag for StartOfCycle operations
	  //lines for the DQM GUI
	  TLine* fLineExpTimeMin;
	  TLine* fLineExpTimeMax;
	  TLine* fLineExpTotMin;
	  TLine* fLineExpTotMax;
	  TLine* fLineSMid[17];
	  
	  AliTOFRawStream fTOFRawStream; // AliTOFRawStream variable
	  //counters for TRM - usage disabled
	  /* Short_t fTRMNoisyArray[720]; */
	  /* Short_t fTRMHwOkArray[720]; */
	  /* Short_t fTRMEnabledArray[720]; */
	  
	  AliTOFDecoderSummaryData * fDecoderSummary;

	  ClassDef(AliTOFQADataMakerRec,6)  // description 
	    
};

#endif // ALITOFQADATAMAKERREC_H
