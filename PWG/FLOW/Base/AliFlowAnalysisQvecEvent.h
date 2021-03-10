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

#ifndef AliFlowAnalysisQvecEvent_H
#define AliFlowAnalysisQvecEvent_H

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

class AliFlowAnalysisQvecEvent {
 public:
  AliFlowAnalysisQvecEvent();
  virtual ~AliFlowAnalysisQvecEvent();
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
  
  void setZCRe( Double_t ZCRe ) { this->fZCRe = ZCRe; }
  Double_t getZCRe() { return fZCRe; }
  void setZCIm( Double_t ZCIm ) { this->fZCIm = ZCIm; }
  Double_t getZCIm() { return fZCIm; }
  void setZCM( Double_t ZCM ) { this->fZCM = ZCM; }
  Double_t getZCM() { return fZCM; }
  
  void setZARe( Double_t ZARe ) { this->fZARe = ZARe; }
  Double_t getZARe() { return fZARe; }
  void setZAIm( Double_t ZAIm ) { this->fZAIm = ZAIm; }
  Double_t getZAIm() { return fZAIm; }
  void setZAM( Double_t ZAM ) { this->fZAM = ZAM; }
  Double_t getZAM() { return fZAM; }
  
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
  
  void setTPCRePosChPosEta( Double_t TPCRePosChPosEta ) { this->fTPCRePosChPosEta = TPCRePosChPosEta; }
  Double_t getTPCRePosChPosEta() { return fTPCRePosChPosEta; }
  void setTPCImPosChPosEta( Double_t TPCImPosChPosEta ) { this->fTPCImPosChPosEta = TPCImPosChPosEta; }
  Double_t getTPCImPosChPosEta() { return fTPCImPosChPosEta; }
  void setTPCMPosChPosEta( Double_t TPCMPosChPosEta ) { this->fTPCMPosChPosEta = TPCMPosChPosEta; }
  Double_t getTPCMPosChPosEta() { return fTPCMPosChPosEta; }
  void setTPCRePosChNegEta( Double_t TPCRePosChNegEta ) { this->fTPCRePosChNegEta = TPCRePosChNegEta; }
  Double_t getTPCRePosChNegEta() { return fTPCRePosChNegEta; }
  void setTPCImPosChNegEta( Double_t TPCImPosChNegEta ) { this->fTPCImPosChNegEta = TPCImPosChNegEta; }
  Double_t getTPCImPosChNegEta() { return fTPCImPosChNegEta; }
  void setTPCMPosChNegEta( Double_t TPCMPosChNegEta ) { this->fTPCMPosChNegEta = TPCMPosChNegEta; }
  Double_t getTPCMPosChNegEta() { return fTPCMPosChNegEta; }
  void setTPCReNegChPosEta( Double_t TPCReNegChPosEta ) { this->fTPCReNegChPosEta = TPCReNegChPosEta; }
  Double_t getTPCReNegChPosEta() { return fTPCReNegChPosEta; }
  void setTPCImNegChPosEta( Double_t TPCImNegChPosEta ) { this->fTPCImNegChPosEta = TPCImNegChPosEta; }
  Double_t getTPCImNegChPosEta() { return fTPCImNegChPosEta; }
  void setTPCMNegChPosEta( Double_t TPCMNegChPosEta ) { this->fTPCMNegChPosEta = TPCMNegChPosEta; }
  Double_t getTPCMNegChPosEta() { return fTPCMNegChPosEta; }
  void setTPCReNegChNegEta( Double_t TPCReNegChNegEta ) { this->fTPCReNegChNegEta = TPCReNegChNegEta; }
  Double_t getTPCReNegChNegEta() { return fTPCReNegChNegEta; }
  void setTPCImNegChNegEta( Double_t TPCImNegChNegEta ) { this->fTPCImNegChNegEta = TPCImNegChNegEta; }
  Double_t getTPCImNegChNegEta() { return fTPCImNegChNegEta; }
  void setTPCMNegChNegEta( Double_t TPCMNegChNegEta ) { this->fTPCMNegChNegEta = TPCMNegChNegEta; }
  Double_t getTPCMNegChNegEta() { return fTPCMNegChNegEta; }

  void setTPC2RePosChPosEta( Double_t TPC2RePosChPosEta ) { this->fTPC2RePosChPosEta = TPC2RePosChPosEta; }
  Double_t getTPC2RePosChPosEta() { return fTPC2RePosChPosEta; }
  void setTPC2ImPosChPosEta( Double_t TPC2ImPosChPosEta ) { this->fTPC2ImPosChPosEta = TPC2ImPosChPosEta; }
  Double_t getTPC2ImPosChPosEta() { return fTPC2ImPosChPosEta; }
  void setTPC2Re2PosChPosEta( Double_t TPC2Re2PosChPosEta ) { this->fTPC2Re2PosChPosEta = TPC2Re2PosChPosEta; }
  Double_t getTPC2Re2PosChPosEta() { return fTPC2Re2PosChPosEta; }
  void setTPC2Im2PosChPosEta( Double_t TPC2Im2PosChPosEta ) { this->fTPC2Im2PosChPosEta = TPC2Im2PosChPosEta; }
  Double_t getTPC2Im2PosChPosEta() { return fTPC2Im2PosChPosEta; }
  void setTPC2MPosChPosEta( Double_t TPC2MPosChPosEta ) { this->fTPC2MPosChPosEta = TPC2MPosChPosEta; }
  Double_t getTPC2MPosChPosEta() { return fTPC2MPosChPosEta; }
  void setTPC2RePosChNegEta( Double_t TPC2RePosChNegEta ) { this->fTPC2RePosChNegEta = TPC2RePosChNegEta; }
  Double_t getTPC2RePosChNegEta() { return fTPC2RePosChNegEta; }
  void setTPC2ImPosChNegEta( Double_t TPC2ImPosChNegEta ) { this->fTPC2ImPosChNegEta = TPC2ImPosChNegEta; }
  Double_t getTPC2ImPosChNegEta() { return fTPC2ImPosChNegEta; }
  void setTPC2Re2PosChNegEta( Double_t TPC2Re2PosChNegEta ) { this->fTPC2Re2PosChNegEta = TPC2Re2PosChNegEta; }
  Double_t getTPC2Re2PosChNegEta() { return fTPC2Re2PosChNegEta; }
  void setTPC2Im2PosChNegEta( Double_t TPC2Im2PosChNegEta ) { this->fTPC2Im2PosChNegEta = TPC2Im2PosChNegEta; }
  Double_t getTPC2Im2PosChNegEta() { return fTPC2Im2PosChNegEta; }
  void setTPC2MPosChNegEta( Double_t TPC2MPosChNegEta ) { this->fTPC2MPosChNegEta = TPC2MPosChNegEta; }
  Double_t getTPC2MPosChNegEta() { return fTPC2MPosChNegEta; }
  void setTPC2ReNegChPosEta( Double_t TPC2ReNegChPosEta ) { this->fTPC2ReNegChPosEta = TPC2ReNegChPosEta; }
  Double_t getTPC2ReNegChPosEta() { return fTPC2ReNegChPosEta; }
  void setTPC2ImNegChPosEta( Double_t TPC2ImNegChPosEta ) { this->fTPC2ImNegChPosEta = TPC2ImNegChPosEta; }
  Double_t getTPC2ImNegChPosEta() { return fTPC2ImNegChPosEta; }
  void setTPC2Re2NegChPosEta( Double_t TPC2Re2NegChPosEta ) { this->fTPC2Re2NegChPosEta = TPC2Re2NegChPosEta; }
  Double_t getTPC2Re2NegChPosEta() { return fTPC2Re2NegChPosEta; }
  void setTPC2Im2NegChPosEta( Double_t TPC2Im2NegChPosEta ) { this->fTPC2Im2NegChPosEta = TPC2Im2NegChPosEta; }
  Double_t getTPC2Im2NegChPosEta() { return fTPC2Im2NegChPosEta; }
  void setTPC2MNegChPosEta( Double_t TPC2MNegChPosEta ) { this->fTPC2MNegChPosEta = TPC2MNegChPosEta; }
  Double_t getTPC2MNegChPosEta() { return fTPC2MNegChPosEta; }
  void setTPC2ReNegChNegEta( Double_t TPC2ReNegChNegEta ) { this->fTPC2ReNegChNegEta = TPC2ReNegChNegEta; }
  Double_t getTPC2ReNegChNegEta() { return fTPC2ReNegChNegEta; }
  void setTPC2ImNegChNegEta( Double_t TPC2ImNegChNegEta ) { this->fTPC2ImNegChNegEta = TPC2ImNegChNegEta; }
  Double_t getTPC2ImNegChNegEta() { return fTPC2ImNegChNegEta; }
  void setTPC2Re2NegChNegEta( Double_t TPC2Re2NegChNegEta ) { this->fTPC2Re2NegChNegEta = TPC2Re2NegChNegEta; }
  Double_t getTPC2Re2NegChNegEta() { return fTPC2Re2NegChNegEta; }
  void setTPC2Im2NegChNegEta( Double_t TPC2Im2NegChNegEta ) { this->fTPC2Im2NegChNegEta = TPC2Im2NegChNegEta; }
  Double_t getTPC2Im2NegChNegEta() { return fTPC2Im2NegChNegEta; }
  void setTPC2MNegChNegEta( Double_t TPC2MNegChNegEta ) { this->fTPC2MNegChNegEta = TPC2MNegChNegEta; }
  Double_t getTPC2MNegChNegEta() { return fTPC2MNegChNegEta; }
  
 private:
  Int_t fRunNum;
  Double_t fCentrality;
  Double_t fVtxPosX;
  Double_t fVtxPosY;
  Double_t fVtxPosZ;
  // ZDC-C (eta < -8.8)
  Double_t fZCRe;
  Double_t fZCIm;
  Double_t fZCM;
  // ZDC-A (eta > 8.8)
  Double_t fZARe;
  Double_t fZAIm;
  Double_t fZAM;
  // VZ eta < 0
  Double_t fVZCRe;
  Double_t fVZCIm;
  Double_t fVZCM;
  // VZ eta > 0
  Double_t fVZARe;
  Double_t fVZAIm;
  Double_t fVZAM;
  
  // TPC (|eta|<0.8)
  Double_t fTPCRePosChPosEta; // w * cos(theta+) eta+
  Double_t fTPCImPosChPosEta; // w * sin(theta+) eta+
  Double_t fTPCMPosChPosEta;   // w eta+
  Double_t fTPCRePosChNegEta; // w * cos(theta+) eta-
  Double_t fTPCImPosChNegEta; // w * sin(theta+) eta-
  Double_t fTPCMPosChNegEta;    // w eta-
  Double_t fTPCReNegChPosEta; // w * cos(theta-) eta+
  Double_t fTPCImNegChPosEta; // w * sin(theta-) eta+
  Double_t fTPCMNegChPosEta;    // w eta+
  Double_t fTPCReNegChNegEta; // w * cos(theta-) eta-
  Double_t fTPCImNegChNegEta; // w * sin(theta-) eta-
  Double_t fTPCMNegChNegEta;    // w eta-
  
  Double_t fTPC2RePosChPosEta; // w * cos(2theta+) eta+
  Double_t fTPC2ImPosChPosEta; // w * sin(2theta+) eta+
  Double_t fTPC2Re2PosChPosEta; // w^2 * cos(2theta+) eta+
  Double_t fTPC2Im2PosChPosEta; // w^2 * sin(2theta+) eta+
  Double_t fTPC2MPosChPosEta;   // w^2 eta+
  Double_t fTPC2RePosChNegEta; // w * cos(2theta+) eta-
  Double_t fTPC2ImPosChNegEta; // w * sin(2theta+) eta-
  Double_t fTPC2Re2PosChNegEta; // w^2 * cos(2theta+) eta-
  Double_t fTPC2Im2PosChNegEta; // w^2 * sin(2theta+) eta-
  Double_t fTPC2MPosChNegEta;    // w^2 eta-
  Double_t fTPC2ReNegChPosEta; // w * cos(2theta-) eta+
  Double_t fTPC2ImNegChPosEta; // w * sin(2theta-) eta+
  Double_t fTPC2Re2NegChPosEta; // w^2 * cos(2theta-) eta+
  Double_t fTPC2Im2NegChPosEta; // w^2 * sin(2theta-) eta+
  Double_t fTPC2MNegChPosEta;    // w^2 eta+
  Double_t fTPC2ReNegChNegEta; // w * cos(2theta-) eta-
  Double_t fTPC2ImNegChNegEta; // w * sin(2theta-) eta-
  Double_t fTPC2Re2NegChNegEta; // w^2 * cos(2theta-) eta-
  Double_t fTPC2Im2NegChNegEta; // w^2 * sin(2theta-) eta-
  Double_t fTPC2MNegChNegEta;    // w^2 eta-
  
  ClassDef(AliFlowAnalysisQvecEvent, 1);
};

#endif
