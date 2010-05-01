#ifndef ALIEMCALQACHECKER_H
#define ALIEMCALQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/*
  Checks the quality assurance. 
  By comparing with reference data

  Based on PHOS code written by
  Y. Schutz CERN July 2007
  
  Implemented by Yaxian Mao (CERN Oct. 2009)
*/


// --- ROOT system ---
class TFile ; 
class TH1F ; 
class TH1I ; 
class TH1 ;
class TLine ;
class TText ;
class TObjArray;
class TPaveText ;

// --- Standard library ---
// --- AliRoot header files ---
#include "AliQACheckerBase.h"
class AliEMCALLoader ; 

class AliEMCALQAChecker: public AliQACheckerBase {

public:
  //Histograms for Raw data control
  enum HRawType_t { 
    // first normal Low Gain and High Gain info
    kNsmodLG,kNsmodHG,kTimeLG,kTimeHG,
    kSigLG,kSigHG,kNtotLG,kNtotHG,kTowerHG,kTowerLG,
    kPedLG,kPedHG,
    k2DRatioAmp,kRatioDist, kLEDMonRatio, kLEDMonRatioDist,
    // then TRU info
    kNsmodTRU,kTimeTRU,
    kSigTRU,kNtotTRU,
    kPedTRU,
    kNL0TRU, kTimeL0TRU,
    // and also LED Mon info
    kNsmodLGLEDMon,kNsmodHGLEDMon,kTimeLGLEDMon,kTimeHGLEDMon,
    kSigLGLEDMon,kSigHGLEDMon,kNtotLGLEDMon,kNtotHGLEDMon,
    kPedLGLEDMon,kPedHGLEDMon
  } ;
	//Histograms for RecPoints  control
	enum HRPType_t {kRecPE,kRecPM,kRecPDigM};
	
	//Histograms for ESDs  control
	enum HESDType_t {kESDCaloClusE,kESDCaloClusM,kESDCaloCellA,kESDCaloCellM} ;

  AliEMCALQAChecker() ;        // ctor
  AliEMCALQAChecker(const AliEMCALQAChecker& qac);
  AliEMCALQAChecker& operator = (const AliEMCALQAChecker& qac) ;
  virtual ~AliEMCALQAChecker() ; // dtor
	
  virtual void   Init(const AliQAv1::DETECTORINDEX_t det) ; 

protected:
		
	virtual void Check( Double_t * test, AliQAv1::ALITASK_t index, TObjArray ** list,  const AliDetectorRecoParam * /*recoParam*/) ;
	//virtual void SetQA(AliQAv1::ALITASK_t index, Double_t * value) const ;	
	
  void CheckRaws(Double_t* test, TObjArray ** list);
	void CheckRecPoints(Double_t* test, TObjArray** /*list*/){;}
	void CheckESD(Double_t* test, TObjArray** /*list*/){;}
	TH1* GetHisto(TObjArray* list, const char* hname, Int_t specie) const;
	Double_t MarkHisto(TH1& histo, Double_t value) const;
	
	
private:
	//TH1F * htemp; //a tempory histrogram for getting the mean and sigma
	//Double_t fMean; //mean value 
	//Double_t fWidth; //sigma of the distribution
  static const Int_t fknSM = 4;    //! number of current SM
//  TLine **     fLine       ; //! line to distinguish the different SM
//  TLine **     fHref       ; //! Line marking the average value for each SM
  TText **    fTextSM        ; //! Text info for each SM
  TLine *     fLineCol       ; //! line to distinguish the different SM side: A side and C side
  TLine *     fLineRow       ; //! line to distinguish the different SM sector 0 and 1 
  TPaveText * fText          ;  //! Information text for the quality of each SM
  ClassDef(AliEMCALQAChecker,2)  // description 

};

#endif // AliEMCALQAChecker_H
