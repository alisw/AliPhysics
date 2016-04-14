#ifndef ALIEMCALQACHECKER_H
#define ALIEMCALQACHECKER_H

//===============================================================
/// \class AliEMCALQAChecker
/// \brief QA checker class
///
///  Checks the quality assurance. 
///  By comparing with reference data
///
///  Based on PHOS code written by
///  Y. Schutz CERN July 2007
///  
///  For the moment we only implement the checking of raw data QA.
///  The checked for ESD and RecPoints will be implemented later.
///
///  \author Yaxian Mao, <Yaxian.Mao@cern.ch>, CCNU, 2009-2011.
///
//===============================================================

//  --- ROOT system ---
class TFile ; 
class TH1F ; 
class TH1I ; 
class TH1 ;
class TLine ;
class TText ;
class TObjArray;
class TPaveText ;

//  --- AliRoot header files ---
#include "AliQACheckerBase.h"
#include "AliEMCALTriggerMappingV2.h"
class AliEMCALLoader ; 

class AliEMCALQAChecker: public AliQACheckerBase {

public:
  
  // Histograms for Raw data control
  enum HRawType_t { 
    //  first normal Low Gain and High Gain info
    kNsmodLG,kNsmodHG,kTimeLG,kTimeHG,
    kNtotLG,kNtotHG,kSigHG,kSigLG,
    kPedLG,kPedHG,
    k2DRatioAmp,kRatioDist, kLEDMonRatio, kLEDMonRatioDist,
    //  then TRU info
    kNsmodTRU,
    kSigTRU,kNtotTRU,
    kNL0TRU, kTimeL0TRU,
    kNL0FirstTRU, kTimeL0FirstTRU,kNL0TRUSamples,kNL0TRURMS,
    //  and also LED Mon info
    kNsmodLGLEDMon,kNsmodHGLEDMon,kTimeLGLEDMon,kTimeHGLEDMon,
    kSigLGLEDMon,kSigHGLEDMon,kNtotLGLEDMon,kNtotHGLEDMon,
    kPedLGLEDMon,kPedHGLEDMon,
    // and STU info
    kAmpL1, kGL1, kJL1,
    kGL1V0, kJL1V0, kSTUTRU  
  } ;
  
  // Histograms for RecPoints  control
  enum HRPType_t {kRecPE,kRecPM,kRecPDigM};
	
  // Histograms for ESDs  control
  enum HESDType_t {kESDCaloClusE,kESDCaloClusM,kESDCaloCellA,kESDCaloCellM} ;

  AliEMCALQAChecker() ;          // ctor
  
  virtual ~AliEMCALQAChecker() ; // dtor
	
  virtual void   Init(const AliQAv1::DETECTORINDEX_t det) ; 

protected:
		
  virtual void Check( Double_t * test, AliQAv1::ALITASK_t index, TObjArray ** list,  
                     const AliDetectorRecoParam * /*recoParam*/) ;
  //virtual void SetQA(AliQAv1::ALITASK_t index, Double_t * value) const ;	
	
  void CheckRaws(Double_t* test, TObjArray ** list);
  
  void CheckRecPoints(Double_t* /*test*/, TObjArray** /*list*/) const {;}

  void CheckESD(Double_t* /*test*/, TObjArray** /*list*/) const {;}
	
  void CleanListOfFunctions(TList *list);
  
  TH1* GetHisto(TObjArray* list, const char* hname, Int_t specie) const;
  
  Double_t MarkHisto(TH1& histo, Double_t value) const;
	
private:
  
  AliEMCALQAChecker(             const AliEMCALQAChecker& qac);
  AliEMCALQAChecker& operator = (const AliEMCALQAChecker& qac) ;
	
//TH1F * htemp; //a tempory histrogram for getting the mean and sigma
//Double_t fMean; //mean value 
//Double_t fWidth; //sigma of the distribution
  
  static const Int_t fgknSM = 20;       //!<! number of current SM; EMCal + DCal
  static const Int_t fgknSectLines = 9; //!<! lines between SM sectors
  
//TLine **    fLine          ; //!<! line to distinguish the different SM
//TLine **    fHref          ; //!<! Line marking the average value for each SM
  TText **    fTextSM        ; //!<! Text info for each SM
  TLine *     fLineCol       ; //!<! line to distinguish the different SM side: A side and C side
  TLine *     fLineRow[fgknSectLines] ; //!<! line to distinguish the different SM sectors (0-5, 9-12) 
  TPaveText * fText          ; //!<! Information text for the quality of each SM
	TPaveText * fTextL1[3]     ; //!<! Information text for the quality of L1 plots (3 plots in total)
  
  /// \cond CLASSIMP	
  ClassDef(AliEMCALQAChecker,5) ;  
  /// \endcond
};

#endif /// AliEMCALQAChecker_H
