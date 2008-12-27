#ifndef ALITPCDATAQA_H
#define ALITPCDATAQA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



#include <TH1F.h>
#include <TObjArray.h>
#include "AliRecoParam.h"

class TArrayF;
class TH2F;
class TTreeSRedirector;
class AliTPCROC;
class AliTPCCalROC;
class AliTPCRawStream;
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
 void MakeTree(const char *fname="QApad.root");

  //
  Bool_t ProcessEventFast(AliTPCRawStreamFast *rawStreamFast);
  Bool_t ProcessEventFast(AliRawReader        *rawReader);
  Bool_t ProcessEvent(AliTPCRawStream *rawStream);
  Bool_t ProcessEvent(AliRawReader    *rawReader);
  Bool_t ProcessEvent(eventHeaderStruct   *event);

  void   Analyse();
  //
  //
  void SetPedestal(AliTPCCalPad *pedestalCal){ fPedestal = pedestalCal;}
  void SetNoise(AliTPCCalPad *noiseCal){ fNoise = noiseCal;}

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

  //
  AliTPCAltroMapping **GetAltroMapping() { return fMapping; };
  void  SetAltroMapping(AliTPCAltroMapping **mapp) { fMapping = mapp; };
  //
  //
  Int_t  GetFirstTimeBin() const { return fFirstTimeBin; }
  Int_t  GetLastTimeBin()  const { return fLastTimeBin;  }
  Int_t  GetAdcMin()       const { return fAdcMin;       }
  Int_t  GetAdcMax()       const { return fAdcMax;       }
  Int_t  GetEventCounter() const { return fEventCounter; }
  Bool_t GetIsAnalysed()   const { return fIsAnalysed;   }
  void  SetRangeTime(Int_t tMin, Int_t tMax){ fFirstTimeBin=tMin; fLastTimeBin=tMax;}  // Set time bin range that is used for the pedestal calibration
  void  SetRangeAdc (Int_t aMin, Int_t aMax){ fAdcMin=aMin; fAdcMax=aMax; }  // Set adc range for the pedestal calibration



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
	       Int_t& timeMin,Int_t& timeMax,Int_t& padMin,Int_t& padMax);

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

  Int_t   fEventCounter;            // event Counter
  Bool_t  fIsAnalysed;              // Set to true after Analyse has been called
  //
  //  Expand buffer
  //
  Float_t** fAllBins;              //! array for digit using random access
  Int_t**   fAllSigBins;           //! array of pointers to the indexes over threshold
  Int_t*    fAllNSigBins;          //! 
  Int_t fRowsMax;                  //!  Maximum number of time bins
  Int_t fPadsMax;                  //!  Maximum number of time bins
  Int_t fTimeBinsMax;              //!  Maximum number of time bins


public:
  ClassDef(AliTPCdataQA, 3)  // Implementation of the TPC pedestal and noise calibration
};



#endif

