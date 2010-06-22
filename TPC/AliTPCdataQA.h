#ifndef ALITPCDATAQA_H
#define ALITPCDATAQA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



#include <TH1F.h>
class TProfile;
#include "AliRecoParam.h"

class TArrayF;
class TH2F;
class TTreeSRedirector;
class AliTPCROC;
class AliTPCCalROC;
class AliTPCRawStream;
class AliTPCRawStreamV3;
class AliTPCRawStreamFast;
class AliRawReader;
class AliTPCAltroMapping;
class AliTPCCalPad;
class TMap; 
struct eventHeaderStruct;

class AliTPCdataQA : public TH1F {

public:
  AliTPCdataQA();
  AliTPCdataQA(AliRecoParam::EventSpecie_t es);
  AliTPCdataQA(const AliTPCdataQA &ped);
  AliTPCdataQA(const TMap *config);
  virtual ~AliTPCdataQA();

  AliTPCdataQA& operator = (const  AliTPCdataQA &source);
 void  DumpToFile(const Char_t *filename, const Char_t *dir="", Bool_t append=kFALSE);
 void MakeTree(const char *fname="QApad.root") const;

  //
  Bool_t ProcessEventFast(AliTPCRawStreamFast *const rawStreamFast);
  Bool_t ProcessEventFast(AliRawReader        *const rawReader);
  Bool_t ProcessEvent(AliTPCRawStream *const rawStream);
  Bool_t ProcessEvent(AliTPCRawStreamV3 *const rawStreamV3);
  Bool_t ProcessEvent(AliRawReader    *const rawReader);
  Bool_t ProcessEventOld(AliRawReader    *const rawReader);
  Bool_t ProcessEvent(eventHeaderStruct   *const event);

  void   Analyse();
  //
  //
  void SetPedestal(AliTPCCalPad *const pedestalCal){ fPedestal = pedestalCal;}
  void SetNoise(AliTPCCalPad *const noiseCal){ fNoise = noiseCal;}

  AliTPCCalPad *GetNoThreshold() const { return fNoThreshold;}
  AliTPCCalPad *GetMaxCharge() const { return fMaxCharge;}
  AliTPCCalPad *GetMeanCharge() const { return fMeanCharge;}
  AliTPCCalPad *GetNLocalMaxima() const { return fNLocalMaxima;}
  AliTPCCalPad *GetOverThreshold10() const { return fOverThreshold10;}
  AliTPCCalPad *GetOverThreshold20() const { return fOverThreshold20;}
  AliTPCCalPad *GetOverThreshold30() const { return fOverThreshold30;}
  AliTPCCalPad *GetNTimeBins() const { return fNTimeBins;}
  AliTPCCalPad *GetNPads() const { return fNPads;}
  AliTPCCalPad *GetTimePosition() const { return fTimePosition;}
  TProfile* GetHistQVsTimeSideA()    const {return fHistQVsTimeSideA;}
  TProfile* GetHistQVsTimeSideC()    const {return fHistQVsTimeSideC;}
  TProfile* GetHistQMaxVsTimeSideA() const {return fHistQMaxVsTimeSideA;}
  TProfile* GetHistQMaxVsTimeSideC() const {return fHistQMaxVsTimeSideC;}
  TH1F*     GetHistOccupancyVsEventConst() const {return fHistOccupancyVsEvent;}
  TH1F*     GetHistNclustersVsEventConst() const {return fHistNclustersVsEvent;}
  TH1F*     GetHistOccupancyVsEvent();
  TH1F*     GetHistNclustersVsEvent();

  //
  AliTPCAltroMapping **GetAltroMapping() const { return fMapping; };
  void  SetAltroMapping(AliTPCAltroMapping **mapp) { fMapping = mapp; };
  //
  //
  Int_t  GetFirstTimeBin() const { return fFirstTimeBin; }
  Int_t  GetLastTimeBin()  const { return fLastTimeBin;  }
  Int_t  GetAdcMin()       const { return fAdcMin;       }
  Int_t  GetAdcMax()       const { return fAdcMax;       }
  Int_t  GetEventCounter() const { return fEventCounter; }
  Bool_t GetIsAnalysed()   const { return fIsAnalysed;   }
  Int_t  GetMaxEvents()      const { return fMaxEvents;     }
  Int_t  GetEventsPerBin()   const { return fEventsPerBin;  }
  Int_t  GetSignalCounter()  const { return fSignalCounter; }
  Int_t  GetClusterCounter() const { return fClusterCounter;}

  void  SetRangeTime(Int_t tMin, Int_t tMax){ fFirstTimeBin=tMin; fLastTimeBin=tMax;}  // Set time bin range that is used for the pedestal calibration
  void  SetRangeAdc (Int_t aMin, Int_t aMax){ fAdcMin=aMin; fAdcMax=aMax; }  // Set adc range for the pedestal calibration
  void  SetMaxEvents   (Int_t value) { fMaxEvents = value; }
  void  SetEventsPerBin(Int_t value) { fEventsPerBin = value; }

private:
  Int_t Update(const Int_t iSector, const Int_t iRow, const Int_t iPad,
	       const Int_t iTimeBin, Float_t signal);
  void  FindLocalMaxima(const Int_t iSector);

  void MakeArrays();                // Create arrays for random data acces
  void CleanArrays();               // Clean arrays for random data acces
  void SetExpandDigit(const Int_t iRow, Int_t iPad, Int_t iTimeBin, 
		      const Float_t signal); // Fill arrays with signals
  void GetPadAndTimeBin(Int_t bin, Int_t& iPad, Int_t& iTimeBin); // Get pad and time bin corresponding to the 1d bin
  Float_t GetQ(const Float_t* adcArray, const Int_t time,
	       const Int_t pad, const Int_t maxTimeBins, 
	       Int_t& timeMin,Int_t& timeMax,Int_t& padMin,Int_t& padMax) const;
  void UpdateEventHistograms();

  Int_t fFirstTimeBin;              //  First Time bin needed for analysis
  Int_t fLastTimeBin;               //  Last Time bin needed for analysis
  Int_t fAdcMin;                    //  min adc channel of pedestal value
  Int_t fAdcMax;                    //  max adc channel of pedestal value

  AliTPCAltroMapping **fMapping;    //! Altro Mapping object
  //
  //
  AliTPCCalPad * fPedestal;         //! option to set pedestal cal object
  AliTPCCalPad * fNoise;            //! option to set noise cal object
  AliTPCCalPad * fNLocalMaxima;     // local maximas found
  AliTPCCalPad * fMaxCharge;        // max charge
  AliTPCCalPad * fMeanCharge;       // mean charge
  AliTPCCalPad * fNoThreshold;      // number of digits
  AliTPCCalPad * fNTimeBins;        // timebins width of cluster
  AliTPCCalPad * fNPads;            // pads with of cluster
  AliTPCCalPad * fTimePosition;     // Time position of local maximum
  AliTPCCalPad * fOverThreshold10;  //! local maxima with qMax over threshold
  AliTPCCalPad * fOverThreshold20;  //! local maxima with qMax over threshold
  AliTPCCalPad * fOverThreshold30;  //! local maxima with qMax over threshold

  TProfile* fHistQVsTimeSideA;      // Q vs time (side A)
  TProfile* fHistQVsTimeSideC;	    // Q vs time (side C)
  TProfile* fHistQMaxVsTimeSideA;   // QMax vs time (side A)
  TProfile* fHistQMaxVsTimeSideC;   // QMax vs time (side C)

  TH1F* fHistOccupancyVsEvent;      // Occupancy vs event number (~time)
  TH1F* fHistNclustersVsEvent;      // Nclusters vs event number (~time)

  Int_t   fEventCounter;            // event Counter
  Bool_t  fIsAnalysed;              // Set to true after Analyse has been called

  Int_t fMaxEvents;                 // Max events for event histograms
  Int_t fEventsPerBin;              // Events per bin for event histograms
  Int_t fSignalCounter;             // Signal counter
  Int_t fClusterCounter;            // Cluster counter
  //
  //  Expand buffer
  //
  Float_t** fAllBins;              //! array for digit using random access
  Int_t**   fAllSigBins;           //! array of pointers to the indexes over threshold
  Int_t*    fAllNSigBins;          //! 
  Int_t fRowsMax;                  //!  Maximum number of time bins
  Int_t fPadsMax;                  //!  Maximum number of time bins
  Int_t fTimeBinsMax;              //!  Maximum number of time bins

  ClassDef(AliTPCdataQA, 5)  // Implementation of the TPC Raw QA
};



#endif

