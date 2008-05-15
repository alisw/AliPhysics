#ifndef ALICFEVENTCLASSCUTS_H
#define ALICFEVENTCLASSCUTS_H
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
// Cut on the type of Event Class: for the moment 
// the Trigger (pp running configurations) 
// and requirements on the energy observed in the ZDC are implemented
// The argument of IsSelected member function (passed object) is cast into 
// an AliVEvent, although conditions are meaningful only for AliESD(AOD)Event 
// type objects.
// The class derives from AliCFCutBase
// Author:S.Arcelli Silvia.Arcelli@cern.ch

#include "AliCFCutBase.h"
#include "TBits.h"
class TH1F;
class AliVEvent;
//____________________________________________________________________________
class AliCFEventClassCuts: public AliCFCutBase 
{
 public :
  AliCFEventClassCuts() ;
  AliCFEventClassCuts(Char_t* name, Char_t* title) ;
  AliCFEventClassCuts(const AliCFEventClassCuts& c) ;
  AliCFEventClassCuts& operator=(const AliCFEventClassCuts& c) ;
  ~AliCFEventClassCuts();
  Bool_t IsSelected(TObject* obj);
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
  void Init();
  void GetBitMap(TObject *obj,  TBits *bitmap);
  void AddQAHistograms(TList *list) const;

  //Association to The Trigger bits in the mask. 
  //They correspond to the PP running descriptor as in 
  //STEER/createTriggerDescriptor_pp.C plus five MB Trigger combinations
  enum TriggerType { 
    kVZEROLeft=0,kVZERORight,kVZEROBeamGas,
    kSTARTAL0,kSTARTCL0,
    kITSSPDGFOL0,kITSSPDHMultL0,
    kMUSingleLPtL0,kMUUnlikeLPtL0,kMUUnlikeHPtL0,kMULikeLPtL0,kMULikeHPtL0,
    kMB,kTOFMB,
    kMUSingleMB,kMUUnLikeLPtMB,kMULikeLPtMB,
    kMB1,kMB2,kMB3,kMB4,kMB5
  }; 

  //static checker for trigger bits
  static Bool_t IsTriggered(AliVEvent *ev, TriggerType trigger=kMB1); 

  void   SetTriggerType(TriggerType trigger=kMB1) { fTriggerType.SetBitNumber(trigger,kTRUE);} // Set requested trigger bits
  TBits  GetTriggerType() const { return fTriggerType;} // get Triggers bits which were requested
  void  SetTriggersInAND( Bool_t flag){fTriggerAND=flag;} // request Trigger bits in .AND.

  void  SetZDCN1EnergyCut(Double_t min,Double_t max){fZDCN1EnergyMin=min; fZDCN1EnergyMax=max;} // ZDC energy cuts
  void  SetZDCN2EnergyCut(Double_t min,Double_t max){fZDCN2EnergyMin=min; fZDCN2EnergyMax=max;} // ZDC energy cuts
  void  SetZDCP1EnergyCut(Double_t min,Double_t max){fZDCP1EnergyMin=min; fZDCP1EnergyMax=max;} // ZDC energy cuts
  void  SetZDCP2EnergyCut(Double_t min,Double_t max){fZDCP2EnergyMin=min; fZDCP2EnergyMax=max;} // ZDC energy cuts  
  void  SetZDCEM1EnergyCut(Double_t min,Double_t max){fZDCEM1EnergyMin=min; fZDCEM1EnergyMax=max;} // ZDC energy cuts
  void  SetZDCEM2EnergyCut(Double_t min,Double_t max){fZDCEM2EnergyMin=min; fZDCEM2EnergyMax=max;} // ZDC energy cuts

  Double_t  GetZDCN1EnergyCutMin() const {return fZDCN1EnergyMin;};//ZDC N1 energy min
  Double_t  GetZDCN2EnergyCutMin() const {return fZDCN2EnergyMin;};//ZDC N2 Emin
  Double_t  GetZDCP1EnergyCutMin() const {return fZDCP1EnergyMin;};//ZDC P1 Emin
  Double_t  GetZDCP2EnergyCutMin() const {return fZDCP2EnergyMin;};//ZDC P2 Emin
  Double_t  GetZDCEM1EnergyCutMin() const {return fZDCEM1EnergyMin;};//ZDC EM1 Emin
  Double_t  GetZDCEM2EnergyCutMin() const {return fZDCEM2EnergyMin;};//ZDC EM2 Emin

  Double_t  GetZDCN1EnergyCutMax() const {return fZDCN1EnergyMax;};//ZDC N1 Emax
  Double_t  GetZDCN2EnergyCutMax() const {return fZDCN2EnergyMax;};//ZDC N2 Emax
  Double_t  GetZDCP1EnergyCutMax() const {return fZDCP1EnergyMax;};//ZDC P1 Emax
  Double_t  GetZDCP2EnergyCutMax() const {return fZDCP2EnergyMax;};//ZDC P2 Emax
  Double_t  GetZDCEM1EnergyCutMax() const {return fZDCEM1EnergyMax;};//ZDC EM1 Emax
  Double_t  GetZDCEM2EnergyCutMax() const {return fZDCEM2EnergyMax;};//ZDC EM2 Emax


  // QA histograms
  void FillHistogramsBeforeCuts(TObject* obj) {return FillHistograms(obj,kFALSE);}
  void FillHistogramsAfterCuts(TObject* obj)  {return FillHistograms(obj,kTRUE);}
  // QA histogram setter
  // please use indices from the enumeration below
  void SetHistogramBins(Int_t index, Int_t nbins, Double_t *bins);
  void SetHistogramBins(Int_t index, Int_t nbins, Double_t xmin, Double_t xmax);
  enum{kTrigger=0,
	 kZDCEnergyN1,
	 kZDCEnergyP1,
	 kZDCEnergyN2,
	 kZDCEnergyP2,
	 kZDCEnergyEM1,
	 kZDCEnergyEM2,
	 kNTriggers=17,
	 kNTriggersMB=5,
	 kNCuts=7,
         kNStepQA=2
	 };
 private:
  TBits* SelectionBitMap(TObject* obj);
  static void TriggerBitMap(AliVEvent* ev,TBits *bitmapT);
  void DefineHistograms(); 		// books histograms and TList
  void Initialise();			// sets everything to 0
  void FillHistograms(TObject* obj, Bool_t b);

  TBits fTriggerType ; //The type of trigger to be checked
  Bool_t fTriggerAND; //Flag to ak for .AND of all the requested trigger bits (.or.is default)
  Double_t fZDCN1EnergyMin;  //Min Energy in ZDCN1
  Double_t fZDCP1EnergyMin;  //Min Energy in ZDCP1
  Double_t fZDCN2EnergyMin;  //Min Energy in ZDCN2
  Double_t fZDCP2EnergyMin;  //Min Energy in ZDCP2
  Double_t fZDCEM1EnergyMin; //Min Energy in ZDCEM1
  Double_t fZDCEM2EnergyMin; //Min Energy in ZDCEM2
  Double_t fZDCN1EnergyMax;  //Max Energy in ZDCN1
  Double_t fZDCP1EnergyMax;  //Max Energy in ZDCP1
  Double_t fZDCN2EnergyMax;  //Max Energy in ZDCN2
  Double_t fZDCP2EnergyMax;  //Max Energy in ZDCP2
  Double_t fZDCEM1EnergyMax; //Max Energy in ZDCEM1
  Double_t fZDCEM2EnergyMax; //Max Energy in ZDCEM2

  TBits *fBitMap ; //cut mask

  TH1F* fhQA[kNCuts][kNStepQA];		// QA Histograms
  //QA Histogram parameters
  Int_t fhNBinsTrigger;//size of array of bin limits, Trigger Mask
  Double_t *fhBinLimTrigger;//[fhNBinsTrigger] bin limits, Trigger Mask
  Int_t fhNBinsZDCEnN1;//size of array of bin limits, Energy in ZDC N1
  Double_t *fhBinLimZDCEnN1;//[fhNBinsZDCEnN1] bin limits, Energy in ZDC N1
  Int_t fhNBinsZDCEnP1;//size of array of bin limits, Energy in ZDC P1
  Double_t *fhBinLimZDCEnP1;//[fhNBinsZDCEnP1] bin limits, Energy in ZDC P1
  Int_t fhNBinsZDCEnN2;//size of array of bin limits, Energy in ZDC N2
  Double_t *fhBinLimZDCEnN2;//[fhNBinsZDCEnN2] bin limits, Energy in ZDC N2
  Int_t fhNBinsZDCEnP2;//size of array of bin limits, Energy in ZDC P2
  Double_t *fhBinLimZDCEnP2;//[fhNBinsZDCEnP2] bin limits, Energy in ZDC P2
  Int_t fhNBinsZDCEnEM1;//size of array of bin limits, Energy in ZDC EM1
  Double_t *fhBinLimZDCEnEM1;//[fhNBinsZDCEnEM1] bin limits, Energy in ZDC EM1
  Int_t fhNBinsZDCEnEM2;//size of array of bin limits, Energy in ZDC EM2
  Double_t *fhBinLimZDCEnEM2;//[fhNBinsZDCEnEM1] bin limits, Energy in ZDC EM2
 
  ClassDef(AliCFEventClassCuts,1);
};
#endif
