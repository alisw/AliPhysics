/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/*************************************************
 * Qvec event                                    *
 *                                               *
 * author: Shi Qiu                               *
 *         (s.qiu@nikhef.nl)                     *
 *************************************************/

#ifndef AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple_H
#define AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple_H

#include <fstream>
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TMath.h"
#include <map>
#include <complex>
#include <cmath>

//==============================================================================================================

class AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple {
 public:
  AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple();
  virtual ~AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple();
  void setRunNum( Int_t runNum ) { this->fRunNum = runNum; }
  Int_t getRunNum() { return fRunNum; }
  void setCentrality( Double_t centrality ) { this->fCentrality = centrality; }
  Double_t getCentrality() { return fCentrality; }
  void setVtxPosX( Double_t vtxPosX ) { this->fVtxPosX = vtxPosX; }
  Double_t getVtxPosX() { return fVtxPosX; }
  void setVtxPosY( Double_t vtxPosY ) { this->fVtxPosY = vtxPosY; }
  Double_t getVtxPosY() { return fVtxPosY; }
  void setVtxPosZ( Double_t vtxPosZ ) { this->fVtxPosZ = vtxPosZ; }
  Double_t getVtxPosZ() { return fVtxPosZ; }
  
  void setOrbitNumber(UInt_t orbit) {this->fOrbitNumber = orbit;}
  UInt_t getOrbitNumber() { return fOrbitNumber; }
  
  void setVZCRe( Double_t VZCRe ) { this->fVZCRe = VZCRe; }
  Double_t getVZCRe() { return fVZCRe; }
  void setVZCIm( Double_t VZCIm ) { this->fVZCIm = VZCIm; }
  Double_t getVZCIm() { return fVZCIm; }
  void setVZCM( Double_t VZCM ) { this->fVZCM = VZCM; }
  Double_t getVZCM() { return fVZCM; }
  
  void setVZARe( Double_t VZARe ) { this->fVZARe = VZARe; }
  Double_t getVZARe() { return fVZARe; }
  void setVZAIm( Double_t VZAIm ) { this->fVZAIm = VZAIm; }
  Double_t getVZAIm() { return fVZAIm; }
  void setVZAM( Double_t VZAM ) { this->fVZAM = VZAM; }
  Double_t getVZAM() { return fVZAM; }
  
  void setTowZNCraw0( Double_t TowZNCraw0 ) { this->fTowZNCraw0 = TowZNCraw0; }
  Double_t getTowZNCraw0() { return fTowZNCraw0; }
  void setTowZNCraw1( Double_t TowZNCraw1 ) { this->fTowZNCraw1 = TowZNCraw1; }
  Double_t getTowZNCraw1() { return fTowZNCraw1; }
  void setTowZNCraw2( Double_t TowZNCraw2 ) { this->fTowZNCraw2 = TowZNCraw2; }
  Double_t getTowZNCraw2() { return fTowZNCraw2; }
  void setTowZNCraw3( Double_t TowZNCraw3 ) { this->fTowZNCraw3 = TowZNCraw3; }
  Double_t getTowZNCraw3() { return fTowZNCraw3; }
  void setTowZNCraw4( Double_t TowZNCraw4 ) { this->fTowZNCraw4 = TowZNCraw4; }
  Double_t getTowZNCraw4() { return fTowZNCraw4; }
  
  void setTowZNAraw0( Double_t TowZNAraw0 ) { this->fTowZNAraw0 = TowZNAraw0; }
  Double_t getTowZNAraw0() { return fTowZNAraw0; }
  void setTowZNAraw1( Double_t TowZNAraw1 ) { this->fTowZNAraw1 = TowZNAraw1; }
  Double_t getTowZNAraw1() { return fTowZNAraw1; }
  void setTowZNAraw2( Double_t TowZNAraw2 ) { this->fTowZNAraw2 = TowZNAraw2; }
  Double_t getTowZNAraw2() { return fTowZNAraw2; }
  void setTowZNAraw3( Double_t TowZNAraw3 ) { this->fTowZNAraw3 = TowZNAraw3; }
  Double_t getTowZNAraw3() { return fTowZNAraw3; }
  void setTowZNAraw4( Double_t TowZNAraw4 ) { this->fTowZNAraw4 = TowZNAraw4; }
  Double_t getTowZNAraw4() { return fTowZNAraw4; }
  
  
  // TPC
  void setRPReTPC( Double_t val ) { this->fRPReTPC = val; }
  Double_t getRPReTPC() { return fRPReTPC; }
  void setRPImTPC( Double_t val ) { this->fRPImTPC = val; }
  Double_t getRPImTPC() { return fRPImTPC; }
  void setRPMultTPC( Double_t val ) { this->fRPMultTPC = val; }
  Double_t getRPMultTPC() { return fRPMultTPC; }
  
  void setRP2ReTPC( Double_t val ) { this->fRP2ReTPC = val; }
  Double_t getRP2ReTPC() { return fRP2ReTPC; }
  void setRP2ImTPC( Double_t val ) { this->fRP2ImTPC = val; }
  Double_t getRP2ImTPC() { return fRP2ImTPC; }
  
  void setPOIPosReTPC( Double_t val ) { this->fPOIPosReTPC = val; }
  Double_t getPOIPosReTPC() { return fPOIPosReTPC; }
  void setPOIPosMult( Double_t val ) { this->fPOIPosMult = val; }
  Double_t getPOIPosMult() { return fPOIPosMult; }
  void setPOIPosImTPC( Double_t val ) { this->fPOIPosImTPC = val; }
  Double_t getPOIPosImTPC() { return fPOIPosImTPC; }
  
  void setPOINegReTPC( Double_t val ) { this->fPOINegReTPC = val; }
  Double_t getPOINegReTPC() { return fPOINegReTPC; }
  void setPOINegMult( Double_t val ) { this->fPOINegMult = val; }
  Double_t getPOINegMult() { return fPOINegMult; }
  void setPOINegImTPC( Double_t val ) { this->fPOINegImTPC = val; }
  Double_t getPOINegImTPC() { return fPOINegImTPC; }

  void setPOIPos2ReTPC( Double_t val ) { this->fPOIPos2ReTPC = val; }
  Double_t getPOIPos2ReTPC() { return fPOIPos2ReTPC; }
  void setPOIPos2ImTPC( Double_t val ) { this->fPOIPos2ImTPC = val; }
  Double_t getPOIPos2ImTPC() { return fPOIPos2ImTPC; }
  
  void setPOINeg2ReTPC( Double_t val ) { this->fPOINeg2ReTPC = val; }
  Double_t getPOINeg2ReTPC() { return fPOINeg2ReTPC; }
  void setPOINeg2ImTPC( Double_t val ) { this->fPOINeg2ImTPC = val; }
  Double_t getPOINeg2ImTPC() { return fPOINeg2ImTPC; }
  
  void set2pCorrelatorCos2PsiDiff2PsiV0CRP( Double_t val ) { this->f2pCorrelatorCos2PsiDiff2PsiV0CRP = val; }
  Double_t get2pCorrelatorCos2PsiDiff2PsiV0CRP() { return f2pCorrelatorCos2PsiDiff2PsiV0CRP; }
  void set2pCorrelatorCos2PsiDiff2PsiV0ARP( Double_t val ) { this->f2pCorrelatorCos2PsiDiff2PsiV0ARP = val; }
  Double_t get2pCorrelatorCos2PsiDiff2PsiV0ARP() { return f2pCorrelatorCos2PsiDiff2PsiV0ARP; }
  
  void set2pCorrelatorCos2PsiDiff2PsiZDCCRP( Double_t val ) { this->f2pCorrelatorCos2PsiDiff2PsiZDCCRP = val; }
  Double_t get2pCorrelatorCos2PsiDiff2PsiZDCCRP() { return f2pCorrelatorCos2PsiDiff2PsiZDCCRP; }
  void set2pCorrelatorCos2PsiDiff2PsiZDCARP( Double_t val ) { this->f2pCorrelatorCos2PsiDiff2PsiZDCARP = val; }
  Double_t get2pCorrelatorCos2PsiDiff2PsiZDCARP() { return f2pCorrelatorCos2PsiDiff2PsiZDCARP; }
  void set2pCorrelatorCos2PsiDiff2PsiZDCCARP( Double_t val ) { this->f2pCorrelatorCos2PsiDiff2PsiZDCCARP = val; }
  Double_t get2pCorrelatorCos2PsiDiff2PsiZDCCARP() { return f2pCorrelatorCos2PsiDiff2PsiZDCCARP; }
  
  void setNITCosPsidiff2PsiV0COS( Double_t val ) { this->fNITCosPsidiff2PsiV0COS = val; }
  Double_t getNITCosPsidiff2PsiV0COS() { return fNITCosPsidiff2PsiV0COS; }
  void setNITSinPsidiff2PsiV0COS( Double_t val ) { this->fNITSinPsidiff2PsiV0COS = val; }
  Double_t getNITSinPsidiff2PsiV0COS() { return fNITSinPsidiff2PsiV0COS; }
  void setNITCosPsidiff2PsiV0AOS( Double_t val ) { this->fNITCosPsidiff2PsiV0AOS = val; }
  Double_t getNITCosPsidiff2PsiV0AOS() { return fNITCosPsidiff2PsiV0AOS; }
  void setNITSinPsidiff2PsiV0AOS( Double_t val ) { this->fNITSinPsidiff2PsiV0AOS = val; }
  Double_t getNITSinPsidiff2PsiV0AOS() { return fNITSinPsidiff2PsiV0AOS; }

  void setNITCosPsidiff2PsiZDCCOS( Double_t val ) { this->fNITCosPsidiff2PsiZDCCOS = val; }
  Double_t getNITCosPsidiff2PsiZDCCOS() { return fNITCosPsidiff2PsiZDCCOS; }
  void setNITSinPsidiff2PsiZDCCOS( Double_t val ) { this->fNITSinPsidiff2PsiZDCCOS = val; }
  Double_t getNITSinPsidiff2PsiZDCCOS() { return fNITSinPsidiff2PsiZDCCOS; }
  void setNITCosPsidiff2PsiZDCAOS( Double_t val ) { this->fNITCosPsidiff2PsiZDCAOS = val; }
  Double_t getNITCosPsidiff2PsiZDCAOS() { return fNITCosPsidiff2PsiZDCAOS; }
  void setNITSinPsidiff2PsiZDCAOS( Double_t val ) { this->fNITSinPsidiff2PsiZDCAOS = val; }
  Double_t getNITSinPsidiff2PsiZDCAOS() { return fNITSinPsidiff2PsiZDCAOS; }
  void setNITCosPsidiff2PsiZDCCAOS( Double_t val ) { this->fNITCosPsidiff2PsiZDCCAOS = val; }
  Double_t getNITCosPsidiff2PsiZDCCAOS() { return fNITCosPsidiff2PsiZDCCAOS; }
  void setNITSinPsidiff2PsiZDCCAOS( Double_t val ) { this->fNITSinPsidiff2PsiZDCCAOS = val; }
  Double_t getNITSinPsidiff2PsiZDCCAOS() { return fNITSinPsidiff2PsiZDCCAOS; }

  void setNITCosPsidiff2PsiV0CPOIPos( Double_t val ) { this->fNITCosPsidiff2PsiV0CPOIPos = val; }
  Double_t getNITCosPsidiff2PsiV0CPOIPos() { return fNITCosPsidiff2PsiV0CPOIPos; }
  void setNITSinPsidiff2PsiV0CPOIPos( Double_t val ) { this->fNITSinPsidiff2PsiV0CPOIPos = val; }
  Double_t getNITSinPsidiff2PsiV0CPOIPos() { return fNITSinPsidiff2PsiV0CPOIPos; }
  void setNITCosPsidiff2PsiV0APOIPos( Double_t val ) { this->fNITCosPsidiff2PsiV0APOIPos = val; }
  Double_t getNITCosPsidiff2PsiV0APOIPos() { return fNITCosPsidiff2PsiV0APOIPos; }
  void setNITSinPsidiff2PsiV0APOIPos( Double_t val ) { this->fNITSinPsidiff2PsiV0APOIPos = val; }
  Double_t getNITSinPsidiff2PsiV0APOIPos() { return fNITSinPsidiff2PsiV0APOIPos; }
  
  void setNITCosPsidiff2PsiZDCCPOIPos( Double_t val ) { this->fNITCosPsidiff2PsiZDCCPOIPos = val; }
  Double_t getNITCosPsidiff2PsiZDCCPOIPos() { return fNITCosPsidiff2PsiZDCCPOIPos; }
  void setNITSinPsidiff2PsiZDCCPOIPos( Double_t val ) { this->fNITSinPsidiff2PsiZDCCPOIPos = val; }
  Double_t getNITSinPsidiff2PsiZDCCPOIPos() { return fNITSinPsidiff2PsiZDCCPOIPos; }
  void setNITCosPsidiff2PsiZDCAPOIPos( Double_t val ) { this->fNITCosPsidiff2PsiZDCAPOIPos = val; }
  Double_t getNITCosPsidiff2PsiZDCAPOIPos() { return fNITCosPsidiff2PsiZDCAPOIPos; }
  void setNITSinPsidiff2PsiZDCAPOIPos( Double_t val ) { this->fNITSinPsidiff2PsiZDCAPOIPos = val; }
  Double_t getNITSinPsidiff2PsiZDCAPOIPos() { return fNITSinPsidiff2PsiZDCAPOIPos; }
  void setNITCosPsidiff2PsiZDCCAPOIPos( Double_t val ) { this->fNITCosPsidiff2PsiZDCCAPOIPos = val; }
  Double_t getNITCosPsidiff2PsiZDCCAPOIPos() { return fNITCosPsidiff2PsiZDCCAPOIPos; }
  void setNITSinPsidiff2PsiZDCCAPOIPos( Double_t val ) { this->fNITSinPsidiff2PsiZDCCAPOIPos = val; }
  Double_t getNITSinPsidiff2PsiZDCCAPOIPos() { return fNITSinPsidiff2PsiZDCCAPOIPos; }

  void setNITCosPsidiff2PsiV0CPOINeg( Double_t val ) { this->fNITCosPsidiff2PsiV0CPOINeg = val; }
  Double_t getNITCosPsidiff2PsiV0CPOINeg() { return fNITCosPsidiff2PsiV0CPOINeg; }
  void setNITSinPsidiff2PsiV0CPOINeg( Double_t val ) { this->fNITSinPsidiff2PsiV0CPOINeg = val; }
  Double_t getNITSinPsidiff2PsiV0CPOINeg() { return fNITSinPsidiff2PsiV0CPOINeg; }
  void setNITCosPsidiff2PsiV0APOINeg( Double_t val ) { this->fNITCosPsidiff2PsiV0APOINeg = val; }
  Double_t getNITCosPsidiff2PsiV0APOINeg() { return fNITCosPsidiff2PsiV0APOINeg; }
  void setNITSinPsidiff2PsiV0APOINeg( Double_t val ) { this->fNITSinPsidiff2PsiV0APOINeg = val; }
  Double_t getNITSinPsidiff2PsiV0APOINeg() { return fNITSinPsidiff2PsiV0APOINeg; }

  void setNITCosPsidiff2PsiZDCCPOINeg( Double_t val ) { this->fNITCosPsidiff2PsiZDCCPOINeg = val; }
  Double_t getNITCosPsidiff2PsiZDCCPOINeg() { return fNITCosPsidiff2PsiZDCCPOINeg; }
  void setNITSinPsidiff2PsiZDCCPOINeg( Double_t val ) { this->fNITSinPsidiff2PsiZDCCPOINeg = val; }
  Double_t getNITSinPsidiff2PsiZDCCPOINeg() { return fNITSinPsidiff2PsiZDCCPOINeg; }
  void setNITCosPsidiff2PsiZDCAPOINeg( Double_t val ) { this->fNITCosPsidiff2PsiZDCAPOINeg = val; }
  Double_t getNITCosPsidiff2PsiZDCAPOINeg() { return fNITCosPsidiff2PsiZDCAPOINeg; }
  void setNITSinPsidiff2PsiZDCAPOINeg( Double_t val ) { this->fNITSinPsidiff2PsiZDCAPOINeg = val; }
  Double_t getNITSinPsidiff2PsiZDCAPOINeg() { return fNITSinPsidiff2PsiZDCAPOINeg; }
  void setNITCosPsidiff2PsiZDCCAPOINeg( Double_t val ) { this->fNITCosPsidiff2PsiZDCCAPOINeg = val; }
  Double_t getNITCosPsidiff2PsiZDCCAPOINeg() { return fNITCosPsidiff2PsiZDCCAPOINeg; }
  void setNITSinPsidiff2PsiZDCCAPOINeg( Double_t val ) { this->fNITSinPsidiff2PsiZDCCAPOINeg = val; }
  Double_t getNITSinPsidiff2PsiZDCCAPOINeg() { return fNITSinPsidiff2PsiZDCCAPOINeg; }

  void set2pCorrelatorCosPsiDiff( Double_t val ) { this->f2pCorrelatorCosPsiDiff = val; }
  Double_t get2pCorrelatorCosPsiDiff() { return f2pCorrelatorCosPsiDiff; }
  void set2pCorrelatorCos2PsiDiff( Double_t val ) { this->f2pCorrelatorCos2PsiDiff = val; }
  Double_t get2pCorrelatorCos2PsiDiff() { return f2pCorrelatorCos2PsiDiff; }
  void set2pCorrelatorRPMult( Double_t val ) { this->f2pCorrelatorRPMult = val; }
  Double_t get2pCorrelatorRPMult() { return f2pCorrelatorRPMult; }

  void set2pCorrelatorCosPsiSumPOIOS( Double_t val ) { this->f2pCorrelatorCosPsiSumPOIOS = val; }
  Double_t get2pCorrelatorCosPsiSumPOIOS() { return f2pCorrelatorCosPsiSumPOIOS; }
  void set2pCorrelatorPOIOSMult( Double_t val ) { this->f2pCorrelatorPOIOSMult = val; }
  Double_t get2pCorrelatorPOIOSMult() { return f2pCorrelatorPOIOSMult; }
  void set2pCorrelatorSinPsiSumPOIOS( Double_t val ) { this->f2pCorrelatorSinPsiSumPOIOS = val; }
  Double_t get2pCorrelatorSinPsiSumPOIOS() { return f2pCorrelatorSinPsiSumPOIOS; }
  void set2pCorrelatorCosPsiDiffPOIOS( Double_t val ) { this->f2pCorrelatorCosPsiDiffPOIOS = val; }
  Double_t get2pCorrelatorCosPsiDiffPOIOS() { return f2pCorrelatorCosPsiDiffPOIOS; }
  void set2pCorrelatorCos2PsiDiffPOIOS( Double_t val ) { this->f2pCorrelatorCos2PsiDiffPOIOS = val; }
  Double_t get2pCorrelatorCos2PsiDiffPOIOS() { return f2pCorrelatorCos2PsiDiffPOIOS; }

  void set2pCorrelatorCosPsiSumPOIPP( Double_t val ) { this->f2pCorrelatorCosPsiSumPOIPP = val; }
  Double_t get2pCorrelatorCosPsiSumPOIPP() { return f2pCorrelatorCosPsiSumPOIPP; }
  void set2pCorrelatorPOIPPMult( Double_t val ) { this->f2pCorrelatorPOIPPMult = val; }
  Double_t get2pCorrelatorPOIPPMult() { return f2pCorrelatorPOIPPMult; }
  void set2pCorrelatorSinPsiSumPOIPP( Double_t val ) { this->f2pCorrelatorSinPsiSumPOIPP = val; }
  Double_t get2pCorrelatorSinPsiSumPOIPP() { return f2pCorrelatorSinPsiSumPOIPP; }
  void set2pCorrelatorCosPsiDiffPOIPP( Double_t val ) { this->f2pCorrelatorCosPsiDiffPOIPP = val; }
  Double_t get2pCorrelatorCosPsiDiffPOIPP() { return f2pCorrelatorCosPsiDiffPOIPP; }
  void set2pCorrelatorCos2PsiDiffPOIPP( Double_t val ) { this->f2pCorrelatorCos2PsiDiffPOIPP = val; }
  Double_t get2pCorrelatorCos2PsiDiffPOIPP() { return f2pCorrelatorCos2PsiDiffPOIPP; }
  
  void set2pCorrelatorCosPsiSumPOINN( Double_t val ) { this->f2pCorrelatorCosPsiSumPOINN = val; }
  Double_t get2pCorrelatorCosPsiSumPOINN() { return f2pCorrelatorCosPsiSumPOINN; }
  void set2pCorrelatorPOINNMult( Double_t val ) { this->f2pCorrelatorPOINNMult = val; }
  Double_t get2pCorrelatorPOINNMult() { return f2pCorrelatorPOINNMult; }
  void set2pCorrelatorSinPsiSumPOINN( Double_t val ) { this->f2pCorrelatorSinPsiSumPOINN = val; }
  Double_t get2pCorrelatorSinPsiSumPOINN() { return f2pCorrelatorSinPsiSumPOINN; }
  void set2pCorrelatorCosPsiDiffPOINN( Double_t val ) { this->f2pCorrelatorCosPsiDiffPOINN = val; }
  Double_t get2pCorrelatorCosPsiDiffPOINN() { return f2pCorrelatorCosPsiDiffPOINN; }
  void set2pCorrelatorCos2PsiDiffPOINN( Double_t val ) { this->f2pCorrelatorCos2PsiDiffPOINN = val; }
  Double_t get2pCorrelatorCos2PsiDiffPOINN() { return f2pCorrelatorCos2PsiDiffPOINN; }
  
 private:
  Int_t fRunNum;
  Double_t fCentrality;
  Double_t fVtxPosX;
  Double_t fVtxPosY;
  Double_t fVtxPosZ;
  
  // period, orbit number, bunch cross, time stamp
  UInt_t fOrbitNumber;
  
  // VZ eta < 0
  Double_t fVZCRe;
  Double_t fVZCIm;
  Double_t fVZCM;
  // VZ eta > 0
  Double_t fVZARe;
  Double_t fVZAIm;
  Double_t fVZAM;
  
  // ZNC each tow energy
  Double_t fTowZNCraw0;
  Double_t fTowZNCraw1;
  Double_t fTowZNCraw2;
  Double_t fTowZNCraw3;
  Double_t fTowZNCraw4;
  
  // ZNA each tow energy
  Double_t fTowZNAraw0;
  Double_t fTowZNAraw1;
  Double_t fTowZNAraw2;
  Double_t fTowZNAraw3;
  Double_t fTowZNAraw4;
  
  // TPC 
  Double_t fRPReTPC; // Qvec TPC cos(phi) 
  Double_t fRPImTPC; // Qvec TPC sin(phi) 
  Double_t fRPMultTPC; 
  
  Double_t fRP2ReTPC; // Qvec TPC cos(2phi)
  Double_t fRP2ImTPC; // Qvec TPC sin(2phi)
  
  Double_t fPOIPosReTPC; // Cos(phi_POIPos)
  Double_t fPOIPosMult;
  Double_t fPOIPosImTPC;

  Double_t fPOINegReTPC; // Cos(phi_POINeg)
  Double_t fPOINegMult;
  Double_t fPOINegImTPC;

  Double_t fPOIPos2ReTPC; // Cos(2phi_POIPos)
  Double_t fPOIPos2ImTPC;

  Double_t fPOINeg2ReTPC; // Cos(2phi_POINeg)
  Double_t fPOINeg2ImTPC;

  Double_t f2pCorrelatorCos2PsiDiff2PsiV0CRP; //<cos(2psi1-2phi_V0C)>
  Double_t f2pCorrelatorCos2PsiDiff2PsiV0ARP; //<cos(2psi1-2phi_V0A)>

  Double_t f2pCorrelatorCos2PsiDiff2PsiZDCCRP; //<cos(2psi1-2phi_ZDCC)>
  Double_t f2pCorrelatorCos2PsiDiff2PsiZDCARP; //<cos(2psi1-2phi_ZDCA)>
  Double_t f2pCorrelatorCos2PsiDiff2PsiZDCCARP; //<cos(2psi1-2phi_ZDCCA)>

	    
  Double_t fNITCosPsidiff2PsiV0COS; // <<cos(psi1-2phi_V0C)>> 
  Double_t fNITSinPsidiff2PsiV0COS; // <<sin(psi1-2phi_V0C)>>
  Double_t fNITCosPsidiff2PsiV0AOS; // <<cos(psi1-2phi_V0A)>>
  Double_t fNITSinPsidiff2PsiV0AOS; // <<sin(psi1-2phi_V0A)>>

  Double_t fNITCosPsidiff2PsiZDCCOS; // <<cos(psi1-2phi_ZDCC)>>
  Double_t fNITSinPsidiff2PsiZDCCOS; // <<sin(psi1-2phi_ZDCC)>>
  Double_t fNITCosPsidiff2PsiZDCAOS; // <<cos(psi1-2phi_ZDCA)>>
  Double_t fNITSinPsidiff2PsiZDCAOS; // <<sin(psi1-2phi_ZDCA)>>
  Double_t fNITCosPsidiff2PsiZDCCAOS; // <<cos(psi1-2phi_ZDCCA)>>
  Double_t fNITSinPsidiff2PsiZDCCAOS; // <<sin(psi1-2phi_ZDCCA)>>


  Double_t fNITCosPsidiff2PsiV0CPOIPos; // <<cos(psi1-2phi_V0C)>> 
  Double_t fNITSinPsidiff2PsiV0CPOIPos; // <<sin(psi1-2phi_V0C)>>
  Double_t fNITCosPsidiff2PsiV0APOIPos; // <<cos(psi1-2phi_V0A)>>
  Double_t fNITSinPsidiff2PsiV0APOIPos; // <<sin(psi1-2phi_V0A)>>

  Double_t fNITCosPsidiff2PsiZDCCPOIPos; // <<cos(psi1-2phi_ZDCC)>>
  Double_t fNITSinPsidiff2PsiZDCCPOIPos; // <<sin(psi1-2phi_ZDCC)>>
  Double_t fNITCosPsidiff2PsiZDCAPOIPos; // <<cos(psi1-2phi_ZDCA)>>
  Double_t fNITSinPsidiff2PsiZDCAPOIPos; // <<sin(psi1-2phi_ZDCA)>>
  Double_t fNITCosPsidiff2PsiZDCCAPOIPos; // <<cos(psi1-2phi_ZDCCA)>>
  Double_t fNITSinPsidiff2PsiZDCCAPOIPos; // <<sin(psi1-2phi_ZDCCA)>>

  Double_t fNITCosPsidiff2PsiV0CPOINeg; // <<cos(psi1-2phi_V0C)>> 
  Double_t fNITSinPsidiff2PsiV0CPOINeg; // <<sin(psi1-2phi_V0C)>>
  Double_t fNITCosPsidiff2PsiV0APOINeg; // <<cos(psi1-2phi_V0A)>>
  Double_t fNITSinPsidiff2PsiV0APOINeg; // <<sin(psi1-2phi_V0A)>>

  Double_t fNITCosPsidiff2PsiZDCCPOINeg; // <<cos(psi1-2phi_ZDCC)>>
  Double_t fNITSinPsidiff2PsiZDCCPOINeg; // <<sin(psi1-2phi_ZDCC)>>
  Double_t fNITCosPsidiff2PsiZDCAPOINeg; // <<cos(psi1-2phi_ZDCA)>>
  Double_t fNITSinPsidiff2PsiZDCAPOINeg; // <<sin(psi1-2phi_ZDCA)>>
  Double_t fNITCosPsidiff2PsiZDCCAPOINeg; // <<cos(psi1-2phi_ZDCCA)>>
  Double_t fNITSinPsidiff2PsiZDCCAPOINeg; // <<sin(psi1-2phi_ZDCCA)>>


  Double_t f2pCorrelatorCosPsiDiff; // <cos(dPsi1-dPsi2)>
  Double_t f2pCorrelatorCos2PsiDiff; // <cos(2(dPsi1-dPsi2))>
  Double_t f2pCorrelatorRPMult;

  Double_t f2pCorrelatorCosPsiSumPOIOS; // <cos(dPsi1+dPsi2)>
  Double_t f2pCorrelatorPOIOSMult;
  Double_t f2pCorrelatorSinPsiSumPOIOS; // <sin(dPsi1+dPsi2)>
  Double_t f2pCorrelatorCosPsiDiffPOIOS; // <cos(dPsi1-dPsi2)>
  Double_t f2pCorrelatorCos2PsiDiffPOIOS; // <cos(2(dPsi1-dPsi2))>

  Double_t f2pCorrelatorCosPsiSumPOIPP; // <cos(dPsi1+dPsi2)>
  Double_t f2pCorrelatorPOIPPMult;
  Double_t f2pCorrelatorSinPsiSumPOIPP; // <sin(dPsi1+dPsi2)>
  Double_t f2pCorrelatorCosPsiDiffPOIPP; // <cos(dPsi1-dPsi2)>
  Double_t f2pCorrelatorCos2PsiDiffPOIPP; // <cos(2(dPsi1-dPsi2))>

  Double_t f2pCorrelatorCosPsiSumPOINN; // <cos(dPsi1+dPsi2)>
  Double_t f2pCorrelatorPOINNMult;
  Double_t f2pCorrelatorSinPsiSumPOINN; // <sin(dPsi1+dPsi2)>
  Double_t f2pCorrelatorCosPsiDiffPOINN; // <cos(dPsi1-dPsi2)>
  Double_t f2pCorrelatorCos2PsiDiffPOINN; // <cos(2(dPsi1-dPsi2))>
  
  ClassDef(AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple, 1);
};

#endif
