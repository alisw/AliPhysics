#ifndef ALITPCCALIBRAW_H
#define ALITPCCALIBRAW_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////////////////
//                                                                                    //
//                        TPC ALTRO Header analysis                                   //
//                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////

#include <TVectorF.h>
#include <TObjArray.h>
#include <THnSparse.h>

#include "AliTPCCalibRawBase.h"
#include "AliTPCCalPad.h"
#include "AliTPCROC.h"

class TH2C;
class TH1F;
class TMap;
class TGraph;
class TCanvas;

class AliTPCCalibRaw : public AliTPCCalibRawBase {
public:
  AliTPCCalibRaw();
  AliTPCCalibRaw(const TMap *config);
  
  virtual ~AliTPCCalibRaw();

  enum {kNRCU=216};
  
  virtual Int_t Update(const Int_t isector, const Int_t iRow, const Int_t iPad,
                       const Int_t iTimeBin, const Float_t signal);
  virtual void UpdateDDL();
  virtual void EndEvent();
  virtual void ResetEvent();
  virtual void Analyse();

  UInt_t GetNFailL1Phase()                const  {return fNFailL1Phase;}
  UInt_t GetNFailL1PhaseEvents()          const  {return fNFailL1PhaseEvent;}
  Int_t   GetPeakDetectionMinus() const {return fPeakDetMinus;}
  Int_t   GetPeakDetectionPlus()  const {return fPeakDetPlus;}
  
  const TVectorF* GetALTROL1PhaseEvents() const {return &fArrALTROL1Phase;}

  const TVectorF *GetALTROL1PhaseEventsRCU(Int_t rcu) const {return (TVectorF*)fArrALTROL1PhaseEvent.At(rcu);}
  const TVectorF *GetALTROL1PhaseFailEventsRCU(Int_t rcu) const {return (TVectorF*)fArrALTROL1PhaseFailEvent.At(rcu);}

  const TVectorF *GetOccupancyEvent()          const {return &fVOccupancyEvent;}
  const TVectorF *GetOccupancyEventSensitive() const {return &fVOccupancySenEvent;}
  const TVectorF *GetSignalSumEvent()          const {return &fVSignalSumEvent;}
  const TVectorF *GetSignalSumEventSensitive() const {return &fVSignalSumSenEvent;}
  const TVectorF *GetFiredPadsSensitive()      const {return &fVNfiredPadsSenEvent;}
  const TVectorF *GetEventTimeStamps()         const {return &fVTimeStampEvent;}
  UInt_t GetFirstTimeStamp() const {return fFirstTimeStamp;}
  
  void  SetRangePeakDetection(Int_t minus, Int_t plus) { fPeakDetMinus=minus; fPeakDetPlus=plus;}
  //Phase info
  TH2C *MakeHistL1RCUEvents(Int_t type=0);
  TH2C *MakeHistL1RCUEventsIROC(Int_t type=0);
  TH2C *MakeHistL1RCUEventsOROC(Int_t type=0);
  TH1F *MakeHistL1PhaseDist();
  TVectorF *MakeVectL1PhaseDist();
  //Occupancy info
  TGraph*  MakeGraphOccupancy(const Int_t type=0, const Int_t xType=0);
//   TGraph*  MakeGraphNoiseEvents();
  TCanvas* MakeCanvasOccupancy(const Int_t xType=1, Bool_t sen=kFALSE);

  const THnSparseI *GetHnDrift() const {return fHnDrift;}
//   AliTPCCalPad *CreateCalPadL1Mean();
//   AliTPCCalPad *CreateCalPadL1RMS();
  
  void Merge(AliTPCCalibRaw * const sig);
  virtual Long64_t Merge(TCollection * const list);
  
private:
  Int_t   fPeakDetMinus;             //  Consecutive timebins on rising edge to be regarded as a signal
  Int_t   fPeakDetPlus;              //  Consecutive timebins on falling edge to be regarded as a signal
  UInt_t  fNFailL1Phase;             //Number of failures in L1 phase
  UInt_t  fNFailL1PhaseEvent;        //Number of events with L1 phase failures
  //binning dv hist
  UInt_t  fNSecTime;                 //Number of seconds per bin in time
  UInt_t  fNBinsTime;                //Number of bin in time
  //processing information
  Bool_t    fPadProcessed;           //! if last pead has already been filled for the current pad
  Int_t     fCurrentChannel;         //! current channel processed
  Int_t     fCurrentSector;          //! current sector processed
  Int_t     fLastSector;             //! current sector processed
  Int_t     fCurrentRow;             //! current row processed
  Int_t     fCurrentPad;             //! current pad processed
  Int_t     fLastTimeBinProc;        //! last time bin processed
  Int_t     fPeakTimeBin;            //! time bin with local maximum
  Int_t     fLastSignal;             //! last signal processed
  Int_t     fNOkPlus;                //! number of processed time bins fullfilling peak criteria
  Int_t     fNOkMinus;               //! number of processed time bins fullfilling peak criteria
  Int_t     fNanoSec;                //! current nano seconds stamp
//
  //L1 phase stuff
  TVectorF fArrCurrentPhaseDist;       //!Phase distribution of the current event
  TVectorF fArrCurrentPhase;           //!Current phase of all RCUs
  TVectorF fArrFailEventNumber;        //event numbers of failed events;
  TVectorF fArrALTROL1Phase;           //Array of L1 phases on an event bases;
  TObjArray fArrALTROL1PhaseEvent;     //L1 phase for each RCU and event
  TObjArray fArrALTROL1PhaseFailEvent; //L1 failure for each RCU and event
  //drift velocity stuff
  enum {kHnBinsDV=3};
  THnSparseI *fHnDrift;                //Histogram last time bin vs. ROC, Time
  //occupancy
  TVectorF fVOccupancyEvent;           //occupancy per event (number of samples above threshold)
  TVectorF fVSignalSumEvent;           //occupancy per event (sum of all adc values)
  TVectorF fVOccupancySenEvent;        //occupancy per event (number of samples abouve threshold) in sensitive regions
  TVectorF fVSignalSumSenEvent;        //occupancy per event (sum of all adc values) in sensitive regions
  TVectorF fVNfiredPadsSenEvent;       //number of pads with a signal above threshold in sensitive regions
  TVectorF fVTimeStampEvent;           //timestamp for all events
  
  TVectorF *MakeArrL1PhaseRCU(Int_t rcu, Bool_t force=kFALSE);
  TVectorF *MakeArrL1PhaseFailRCU(Int_t rcu, Bool_t force=kFALSE);
  
  Bool_t IsEdgePad(Int_t sector, Int_t row, Int_t pad) const;
  void CreateDVhist();
  
  AliTPCCalibRaw(const AliTPCCalibRaw &calib);
  AliTPCCalibRaw& operator = (const  AliTPCCalibRaw &source);

  ClassDef(AliTPCCalibRaw,4) //  Analysis of the Altro header information
};

//----------------------
// Inline Functions
//----------------------
inline TVectorF *AliTPCCalibRaw::MakeArrL1PhaseRCU(Int_t rcu, Bool_t force)
{
  TVectorF *arr=(TVectorF*)fArrALTROL1PhaseEvent.UncheckedAt(rcu);
  if (!arr && force) {
    arr=new TVectorF(100);
    fArrALTROL1PhaseEvent.AddAt(arr,rcu);
  }
  return arr;
}
//
inline TVectorF *AliTPCCalibRaw::MakeArrL1PhaseFailRCU(Int_t rcu, Bool_t force)
{
  TVectorF *arr=(TVectorF*)fArrALTROL1PhaseFailEvent.UncheckedAt(rcu);
  if (!arr && force) {
    arr=new TVectorF(100);
    fArrALTROL1PhaseFailEvent.AddAt(arr,rcu);
  }
  return arr;
}
//_____________________________________________________________________
inline Bool_t AliTPCCalibRaw::IsEdgePad(Int_t sector, Int_t row, Int_t pad) const
{
  //
  // return true if pad is on the edge of a row
  //
  Int_t edge1   = 0;
  if ( pad == edge1 ) return kTRUE;
  Int_t edge2   = fROC->GetNPads(sector,row)-1;
  if ( pad == edge2 ) return kTRUE;
  
  return kFALSE;
}


#endif
