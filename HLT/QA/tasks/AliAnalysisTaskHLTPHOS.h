// $Id$

#ifndef ALIANALYSISTASKHLTPHOS_H
#define ALIANALYSISTASKHLTPHOS_H

//* This file is property of and copyright by the ALICE HLT Project *
//* ALICE Experiment at CERN, All rights reserved.                  *
//* See cxx source for full Copyright notice                        *

/** @file AliAnalysisTaskHLTTPC.h
    @author Zhongbao Yin, Kalliopi Kanaki
    @date
    @brief An analysis task to compare the offline and HLT esd trees
*/


// forward declarations
class TH1F;
class TH2F;
class TList;
class AliESDEvent;
class AliESDtrack;
class AliESDRun;
class TObjArray;

#include "AliAnalysisTaskHLTCalo.h"

class AliAnalysisTaskHLTPHOS : public AliAnalysisTaskHLTCalo {
 
public:  
  AliAnalysisTaskHLTPHOS(const char *name);
  virtual ~AliAnalysisTaskHLTPHOS() {}
  
private:
  

  TH2F *fHistOnlTrk2PHOS; //! track to PHOS 2,3,4 modules in (eta, phi)
  TH2F *fHistOfflTrk2PHOS; //! 
  TH2F *fHistOfflTrk2PHOSTrig; //!
  TH2F *fHistOfflTrk2PHOSNoTrig; //!

  static const Float_t fgkPhiMin[5];
  static const Float_t fgkPhiMax[5];
  static const Float_t fgkEtaMin;
  static const Float_t fgkEtaMax;
  static const Float_t fgkNormX[5];
  static const Float_t fgkNormY[5];
  static const Float_t fgkInitPosX[5];
  static const Float_t fgkInitPosY[5];

  /** copy constructor */
  AliAnalysisTaskHLTPHOS(const AliAnalysisTaskHLTPHOS&); 
  /** assignment operator */
  AliAnalysisTaskHLTPHOS& operator=(const AliAnalysisTaskHLTPHOS&); 

  Bool_t IsInPHOS(Int_t iMod, AliESDtrack * trk, Float_t b, TVector3& v);
  void CreateSpecificStuff(TList * fOutputList);
  void DoSpecificStuff(AliESDEvent * evESD, AliESDEvent * evHLTESD);

  Bool_t IsThisDetector(AliESDCaloCluster * cluster);
  Int_t GetClusters(AliESDEvent * event, TRefArray * clusters);

  ClassDef(AliAnalysisTaskHLTPHOS, 0);
};

#endif
