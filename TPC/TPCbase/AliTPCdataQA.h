/// \class AliTPCdataQA

#ifndef ALITPCDATAQA_H
#define ALITPCDATAQA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TBits.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include "AliRecoParam.h"

#include <TArray.h>

class TH2F;
class TTreeSRedirector;
class AliTPCROC;
class AliTPCCalROC;
class AliTPCRawStreamV3;
class AliRawReader;
class AliTPCAltroMapping;
class AliTPCCalPad;
class TMap;
class TObjArray;
struct eventHeaderStruct;

class AliTPCdataQA : public TH1F {

public:
  AliTPCdataQA();
  AliTPCdataQA(const AliTPCdataQA &ped);
  AliTPCdataQA(const TMap *config);
  virtual ~AliTPCdataQA();

  AliTPCdataQA& operator = (const  AliTPCdataQA &source);
 void  DumpToFile(const Char_t *filename, const Char_t *dir="", Bool_t append=kFALSE);
 void MakeTree(const char *fname="QApad.root") const;

  //
  Bool_t ProcessEvent(AliTPCRawStreamV3 *const rawStreamV3);
  Bool_t ProcessEvent(AliRawReader    *const rawReader);
  Bool_t ProcessEvent(eventHeaderStruct   *const event);

  void   Analyse();
  //
  //
  void SetPedestal(AliTPCCalPad *const pedestalCal){ fPedestal = pedestalCal;}
  void SetNoise(AliTPCCalPad *const noiseCal){ fNoise = noiseCal;}

  void SetMinQMax               (Float_t minQmax  ) { fMinQMax                = minQmax; }
  void SetRequireNeighbouringPad(Bool_t  req=kTRUE) { fRequireNeighbouringPad = req;     }

  // DQM methods
  void FillOccupancyProfile();
  void ResetProfiles();


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

  // DQM output
  TProfile* GetHistOccVsSector()  const { return fHistOccVsSector; }
  TProfile2D* GetHistOcc2dVsSector() const { return fHistOcc2dVsSector; }
  TProfile* GetHistQVsSector()    const { return fHistQVsSector; }
  TProfile* GetHistQmaxVsSector() const { return fHistQmaxVsSector; }

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

  Float_t GetMinQMax()                const { return fMinQMax;                }
  Bool_t  GetRequireNeighbouringPad() const { return fRequireNeighbouringPad; }

  // DQM getter
  Bool_t    GetIsDQM() const { return fIsDQM; }

  void  SetRangeTime(Int_t tMin, Int_t tMax){ fFirstTimeBin=tMin; fLastTimeBin=tMax;}  // Set time bin range that is used for the pedestal calibration
  void  SetRangeAdc (Int_t aMin, Int_t aMax){ fAdcMin=aMin; fAdcMax=aMax; }  // Set adc range for the pedestal calibration
  void  SetMaxEvents   (Int_t value) { fMaxEvents = value; }
  void  SetEventsPerBin(Int_t value) { fEventsPerBin = value; }

  // DQM setter
  void  SetIsDQM(Bool_t value) { fIsDQM = value; }

  void ResetData();

  void SetChamberStatus(UInt_t roc, Bool_t status) { fActiveChambers.SetBitNumber(roc,status); }
  Bool_t GetChamberStatus(UInt_t roc) {return fActiveChambers.TestBitNumber(roc);}

  // Merge functionality
  void Merge(AliTPCdataQA * const ce);
  virtual Long64_t Merge(TCollection * const list);


private:
  Int_t Update(const Int_t iSector, const Int_t iRow, const Int_t iPad,
	       const Int_t iTimeBin, Float_t signal,
	       const Int_t iPatch=-1, const Int_t iBranch=-1);
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

  void Init();

  TObjArray *ConfigArrRocs(TObjArray *arr, const Text_t* name);

  Int_t fFirstTimeBin;              ///< First Time bin needed for analysis
  Int_t fLastTimeBin;               ///< Last Time bin needed for analysis
  Int_t fAdcMin;                    ///< min adc channel of pedestal value
  Int_t fAdcMax;                    ///< max adc channel of pedestal value
  Float_t fMinQMax;                 ///< Minimun charge for Maximum ADC in cluster
  Bool_t fRequireNeighbouringPad;   ///< If clusterer should require a neighbouring pad to accept it


  AliTPCAltroMapping **fMapping;    //!<! Altro Mapping object
  //
  //
  AliTPCCalPad * fPedestal;         //!<! option to set pedestal cal object
  AliTPCCalPad * fNoise;            //!<! option to set noise cal object
  AliTPCCalPad * fNLocalMaxima;     ///< local maximas found
  AliTPCCalPad * fMaxCharge;        ///< max charge
  AliTPCCalPad * fMeanCharge;       ///< mean charge
  AliTPCCalPad * fNoThreshold;      ///< number of digits
  AliTPCCalPad * fNTimeBins;        ///< timebins width of cluster
  AliTPCCalPad * fNPads;            ///< pads with of cluster
  AliTPCCalPad * fTimePosition;     ///< Time position of local maximum
  AliTPCCalPad * fOverThreshold10;  //!<! local maxima with qMax over threshold
  AliTPCCalPad * fOverThreshold20;  //!<! local maxima with qMax over threshold
  AliTPCCalPad * fOverThreshold30;  //!<! local maxima with qMax over threshold

  TProfile* fHistQVsTimeSideA;      ///< Q vs time (side A)
  TProfile* fHistQVsTimeSideC;	    ///< Q vs time (side C)
  TProfile* fHistQMaxVsTimeSideA;   ///< QMax vs time (side A)
  TProfile* fHistQMaxVsTimeSideC;   ///< QMax vs time (side C)

  TH1F* fHistOccupancyVsEvent;      ///< Occupancy vs event number (~time)
  TH1F* fHistNclustersVsEvent;      ///< Nclusters vs event number (~time)

  Int_t   fEventCounter;            ///< event Counter
  Bool_t  fIsAnalysed;              ///< Set to true after Analyse has been called

  Int_t fMaxEvents;                 ///< Max events for event histograms
  Int_t fEventsPerBin;              ///< Events per bin for event histograms
  Int_t fSignalCounter;             ///< Signal counter
  Int_t fClusterCounter;            ///< Cluster counter

  TBits fActiveChambers;           ///< configured ROCs

  //
  //  Expand buffer
  //
  Float_t** fAllBins;              //!<! array for digit using random access
  Int_t**   fAllSigBins;           //!<! array of pointers to the indexes over threshold
  Int_t*    fAllNSigBins;          //!<!
  Int_t fRowsMax;                  //!<! Maximum number of time bins
  Int_t fPadsMax;                  //!<! Maximum number of time bins
  Int_t fTimeBinsMax;              //!<! Maximum number of time bins

  // DQM variables
  Bool_t fIsDQM;                   //!<! Is DQM -> Simple output (no 2D!)
  TProfile* fHistOccVsSector;      //!<! Occ vs sector (for DQM only)
  TProfile2D* fHistOcc2dVsSector;  //!<! Occ vs sector 2D (for DQM only)
  TProfile* fHistQVsSector;        //!<! Q vs sector (for DQM only)
  TProfile* fHistQmaxVsSector;     //!<! QMax vs sector (for DQM only)
  TArrayD* fOccVec;                //!<! Occupancy help counter for DQM
  TArrayD* fOccMaxVec;             //!<! Occupancy help normlization for DQM
  TArrayD* fOccVecFine;            //!<! "2D" occupancy help counter for DQM
  TArrayD* fOccMaxVecFine;         //!<! "2D" occupancy help normlization for DQM


  /// \cond CLASSIMP
  ClassDef(AliTPCdataQA, 6)  // Implementation of the TPC Raw QA
  /// \endcond
};



#endif

