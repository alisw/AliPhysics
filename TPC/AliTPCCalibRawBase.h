#ifndef ALITPCCALIBRAWBASE_H
#define ALITPCCALIBRAWBASE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
//                  Raw data processing base class                                     //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>

class AliAltroMapping;
class AliAltroRawStream;
class AliRawReader;
class AliTPCAltroMapping;
class AliTPCRawStreamV3;
class AliTPCRawStream;
class AliTPCROC;
class TTreeSRedirector;
class TCollection;
struct eventHeaderStruct;

class AliTPCCalibRawBase : public TNamed {


public:
  AliTPCCalibRawBase();
  AliTPCCalibRawBase(const AliTPCCalibRawBase &calib);

  AliTPCCalibRawBase& operator = (const  AliTPCCalibRawBase &source);

  virtual ~AliTPCCalibRawBase();
  
  //uses the new decoder which is compatible with the new altro format
  Bool_t ProcessEvent(AliTPCRawStreamV3   * const rawStreamV3);
  Bool_t ProcessEvent(AliRawReader        * const rawReader);
  Bool_t ProcessEvent(eventHeaderStruct   * const event);

  //For using the old decoder use the following functions
  Bool_t ProcessEvent(AliTPCRawStream * const rawStream);
  Bool_t ProcessEventOld(AliRawReader * const rawReader);
  
  virtual Int_t Update(const Int_t /*isector*/, const Int_t /*iRow*/, const Int_t /*iPad*/,
                       const Int_t /*iTimeBin*/, const Float_t /*signal*/) { return 0; }
  virtual void UpdateDDL() {return;}
  virtual void ProcessBunch(const Int_t /*sector*/, const Int_t /*row*/, const Int_t /*pad*/,
                            const Int_t /*length*/, const UInt_t /*startTimeBin*/, const UShort_t* /*signal*/) {return; }
  virtual void Analyse(){ return; }
  
  virtual Long64_t Merge(TCollection * /*list*/) {return 0;}
  void MergeBase(const AliTPCCalibRawBase *calib);
  
  //Setters
  void  SetRangeTime (Int_t firstTimeBin, Int_t lastTimeBin) { fFirstTimeBin=firstTimeBin;   fLastTimeBin=lastTimeBin;  } //Set range in which the signal is expected
  void  SetAltroMapping(AliTPCAltroMapping **mapp) { fMapping = mapp; }
  //
  void SetUseL1Phase(Bool_t useL1Phase=kTRUE) {fUseL1Phase=useL1Phase;}
  //
  void  SetTimeStampEvent(UInt_t timestamp){ fTimeStamp = timestamp; }
  void  SetRunNumber(UInt_t eventnumber){ fRunNumber = eventnumber; }

  //
  Int_t GetFirstTimeBin()   const { return fFirstTimeBin;  }
  Int_t GetLastTimeBin()    const { return fLastTimeBin;   }
  Int_t GetNevents() const { return fNevents; }
  //
  Double_t GetL1Phase()   const {return fAltroL1Phase;}
  Double_t GetL1PhaseTB() const {return fAltroL1PhaseTB;}
  Bool_t   GetUseL1Phase()const {return fUseL1Phase;}
//
  UInt_t GetRunNumber()      const {return fRunNumber;}
  UInt_t GetFirstTimeStamp() const {return fFirstTimeStamp;}
  UInt_t GetLastTimeStamp()  const {return fLastTimeStamp;}
  UInt_t GetTimeStamp()      const {return fTimeStamp;}
  UInt_t GetEventType()      const {return fEventType;}
  //
  AliTPCAltroMapping **GetAltroMapping() { return fMapping; }
  const AliAltroRawStream *GetAltroRawStream() const {return fAltroRawStream;}
  const AliTPCROC *GetTPCROC() const {return fROC;}
  //
  void IncrementNevents(){++fNevents;}
  //
  virtual void DumpToFile(const Char_t *filename, const Char_t *dir="", Bool_t append=kFALSE);
  // debug and debug streamer support
  TTreeSRedirector *GetDebugStreamer();
  void       SetStreamLevel(Int_t streamLevel){fStreamLevel=streamLevel;}
  void       SetDebugLevel(Int_t level) {fDebugLevel = level;}
  Int_t      GetStreamLevel() const {return fStreamLevel;}
  Int_t      GetDebugLevel() const {return fDebugLevel;}

protected:
  Int_t fFirstTimeBin;                //  First Time bin used for analysis
  Int_t fLastTimeBin;                 //  Last Time bin used for analysis
  
  Int_t fNevents;                     //  Number of processed events 
  
  Int_t fDebugLevel;                  //! debug level
  Int_t fStreamLevel;                 //! level of streamer output
  //
  UInt_t fRunNumber;                  // current run number from event header
  UInt_t fFirstTimeStamp;             // First event time stamp
  UInt_t fLastTimeStamp;              // Last event time stamp
  UInt_t fTimeStamp;                  //! time stamp from event header
  UInt_t fEventType;                  //! current event Type from event header
  //
  Double_t fAltroL1Phase;             //! L1 Phase
  Float_t  fAltroL1PhaseTB;           //! L1 Phase in time bins
  Int_t    fCurrRCUId;                //! Current RCU Id
  Int_t    fPrevRCUId;                //! Previous RCU Id
  Int_t    fCurrDDLNum;               //! Current DDL number
  Int_t    fPrevDDLNum;               //! Current DDL number
  Bool_t   fUseL1Phase;               //  use L1 Phase information?
  //
  TTreeSRedirector *fDebugStreamer;   //! debug streamer
  //
  AliAltroRawStream *fAltroRawStream; //! pointer to the altro object
  AliTPCAltroMapping **fMapping;      //! Altro Mapping object

  AliTPCROC *fROC;                    //! ROC information
    
  virtual void EndEvent() {++fNevents; return; } //fNevents should be updated in the derived classes in a proper place
  virtual void ResetEvent(){ return; }           //Reset Event counters
  
  
  ClassDef(AliTPCCalibRawBase,3)      //  Calibration base class for raw data processing
    
};


#endif

