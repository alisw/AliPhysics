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

#ifndef AliAnalysisTaskGammaDeltaPIDSaveQvecEvent_H
#define AliAnalysisTaskGammaDeltaPIDSaveQvecEvent_H

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

class AliAnalysisTaskGammaDeltaPIDSaveQvecEvent {
 public:
  AliAnalysisTaskGammaDeltaPIDSaveQvecEvent();
  virtual ~AliAnalysisTaskGammaDeltaPIDSaveQvecEvent();
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
  
  
  // TPC pion
  void setTPCPionRePosChPosEta( Double_t TPCPionRePosChPosEta ) { this->fTPCPionRePosChPosEta = TPCPionRePosChPosEta; }
  Double_t getTPCPionRePosChPosEta() { return fTPCPionRePosChPosEta; }
  void setTPCPionImPosChPosEta( Double_t TPCPionImPosChPosEta ) { this->fTPCPionImPosChPosEta = TPCPionImPosChPosEta; }
  Double_t getTPCPionImPosChPosEta() { return fTPCPionImPosChPosEta; }
  void setTPCPionMPosChPosEta( Double_t TPCPionMPosChPosEta ) { this->fTPCPionMPosChPosEta = TPCPionMPosChPosEta; }
  Double_t getTPCPionMPosChPosEta() { return fTPCPionMPosChPosEta; }
  void setTPCPionRePosChNegEta( Double_t TPCPionRePosChNegEta ) { this->fTPCPionRePosChNegEta = TPCPionRePosChNegEta; }
  Double_t getTPCPionRePosChNegEta() { return fTPCPionRePosChNegEta; }
  void setTPCPionImPosChNegEta( Double_t TPCPionImPosChNegEta ) { this->fTPCPionImPosChNegEta = TPCPionImPosChNegEta; }
  Double_t getTPCPionImPosChNegEta() { return fTPCPionImPosChNegEta; }
  void setTPCPionMPosChNegEta( Double_t TPCPionMPosChNegEta ) { this->fTPCPionMPosChNegEta = TPCPionMPosChNegEta; }
  Double_t getTPCPionMPosChNegEta() { return fTPCPionMPosChNegEta; }
  void setTPCPionReNegChPosEta( Double_t TPCPionReNegChPosEta ) { this->fTPCPionReNegChPosEta = TPCPionReNegChPosEta; }
  Double_t getTPCPionReNegChPosEta() { return fTPCPionReNegChPosEta; }
  void setTPCPionImNegChPosEta( Double_t TPCPionImNegChPosEta ) { this->fTPCPionImNegChPosEta = TPCPionImNegChPosEta; }
  Double_t getTPCPionImNegChPosEta() { return fTPCPionImNegChPosEta; }
  void setTPCPionMNegChPosEta( Double_t TPCPionMNegChPosEta ) { this->fTPCPionMNegChPosEta = TPCPionMNegChPosEta; }
  Double_t getTPCPionMNegChPosEta() { return fTPCPionMNegChPosEta; }
  void setTPCPionReNegChNegEta( Double_t TPCPionReNegChNegEta ) { this->fTPCPionReNegChNegEta = TPCPionReNegChNegEta; }
  Double_t getTPCPionReNegChNegEta() { return fTPCPionReNegChNegEta; }
  void setTPCPionImNegChNegEta( Double_t TPCPionImNegChNegEta ) { this->fTPCPionImNegChNegEta = TPCPionImNegChNegEta; }
  Double_t getTPCPionImNegChNegEta() { return fTPCPionImNegChNegEta; }
  void setTPCPionMNegChNegEta( Double_t TPCPionMNegChNegEta ) { this->fTPCPionMNegChNegEta = TPCPionMNegChNegEta; }
  Double_t getTPCPionMNegChNegEta() { return fTPCPionMNegChNegEta; }

  void setTPCPion2RePosChPosEta( Double_t TPCPion2RePosChPosEta ) { this->fTPCPion2RePosChPosEta = TPCPion2RePosChPosEta; }
  Double_t getTPCPion2RePosChPosEta() { return fTPCPion2RePosChPosEta; }
  void setTPCPion2ImPosChPosEta( Double_t TPCPion2ImPosChPosEta ) { this->fTPCPion2ImPosChPosEta = TPCPion2ImPosChPosEta; }
  Double_t getTPCPion2ImPosChPosEta() { return fTPCPion2ImPosChPosEta; }
  void setTPCPion2Re2PosChPosEta( Double_t TPCPion2Re2PosChPosEta ) { this->fTPCPion2Re2PosChPosEta = TPCPion2Re2PosChPosEta; }
  Double_t getTPCPion2Re2PosChPosEta() { return fTPCPion2Re2PosChPosEta; }
  void setTPCPion2Im2PosChPosEta( Double_t TPCPion2Im2PosChPosEta ) { this->fTPCPion2Im2PosChPosEta = TPCPion2Im2PosChPosEta; }
  Double_t getTPCPion2Im2PosChPosEta() { return fTPCPion2Im2PosChPosEta; }
  void setTPCPion2MPosChPosEta( Double_t TPCPion2MPosChPosEta ) { this->fTPCPion2MPosChPosEta = TPCPion2MPosChPosEta; }
  Double_t getTPCPion2MPosChPosEta() { return fTPCPion2MPosChPosEta; }
  void setTPCPion2RePosChNegEta( Double_t TPCPion2RePosChNegEta ) { this->fTPCPion2RePosChNegEta = TPCPion2RePosChNegEta; }
  Double_t getTPCPion2RePosChNegEta() { return fTPCPion2RePosChNegEta; }
  void setTPCPion2ImPosChNegEta( Double_t TPCPion2ImPosChNegEta ) { this->fTPCPion2ImPosChNegEta = TPCPion2ImPosChNegEta; }
  Double_t getTPCPion2ImPosChNegEta() { return fTPCPion2ImPosChNegEta; }
  void setTPCPion2Re2PosChNegEta( Double_t TPCPion2Re2PosChNegEta ) { this->fTPCPion2Re2PosChNegEta = TPCPion2Re2PosChNegEta; }
  Double_t getTPCPion2Re2PosChNegEta() { return fTPCPion2Re2PosChNegEta; }
  void setTPCPion2Im2PosChNegEta( Double_t TPCPion2Im2PosChNegEta ) { this->fTPCPion2Im2PosChNegEta = TPCPion2Im2PosChNegEta; }
  Double_t getTPCPion2Im2PosChNegEta() { return fTPCPion2Im2PosChNegEta; }
  void setTPCPion2MPosChNegEta( Double_t TPCPion2MPosChNegEta ) { this->fTPCPion2MPosChNegEta = TPCPion2MPosChNegEta; }
  Double_t getTPCPion2MPosChNegEta() { return fTPCPion2MPosChNegEta; }
  void setTPCPion2ReNegChPosEta( Double_t TPCPion2ReNegChPosEta ) { this->fTPCPion2ReNegChPosEta = TPCPion2ReNegChPosEta; }
  Double_t getTPCPion2ReNegChPosEta() { return fTPCPion2ReNegChPosEta; }
  void setTPCPion2ImNegChPosEta( Double_t TPCPion2ImNegChPosEta ) { this->fTPCPion2ImNegChPosEta = TPCPion2ImNegChPosEta; }
  Double_t getTPCPion2ImNegChPosEta() { return fTPCPion2ImNegChPosEta; }
  void setTPCPion2Re2NegChPosEta( Double_t TPCPion2Re2NegChPosEta ) { this->fTPCPion2Re2NegChPosEta = TPCPion2Re2NegChPosEta; }
  Double_t getTPCPion2Re2NegChPosEta() { return fTPCPion2Re2NegChPosEta; }
  void setTPCPion2Im2NegChPosEta( Double_t TPCPion2Im2NegChPosEta ) { this->fTPCPion2Im2NegChPosEta = TPCPion2Im2NegChPosEta; }
  Double_t getTPCPion2Im2NegChPosEta() { return fTPCPion2Im2NegChPosEta; }
  void setTPCPion2MNegChPosEta( Double_t TPCPion2MNegChPosEta ) { this->fTPCPion2MNegChPosEta = TPCPion2MNegChPosEta; }
  Double_t getTPCPion2MNegChPosEta() { return fTPCPion2MNegChPosEta; }
  void setTPCPion2ReNegChNegEta( Double_t TPCPion2ReNegChNegEta ) { this->fTPCPion2ReNegChNegEta = TPCPion2ReNegChNegEta; }
  Double_t getTPCPion2ReNegChNegEta() { return fTPCPion2ReNegChNegEta; }
  void setTPCPion2ImNegChNegEta( Double_t TPCPion2ImNegChNegEta ) { this->fTPCPion2ImNegChNegEta = TPCPion2ImNegChNegEta; }
  Double_t getTPCPion2ImNegChNegEta() { return fTPCPion2ImNegChNegEta; }
  void setTPCPion2Re2NegChNegEta( Double_t TPCPion2Re2NegChNegEta ) { this->fTPCPion2Re2NegChNegEta = TPCPion2Re2NegChNegEta; }
  Double_t getTPCPion2Re2NegChNegEta() { return fTPCPion2Re2NegChNegEta; }
  void setTPCPion2Im2NegChNegEta( Double_t TPCPion2Im2NegChNegEta ) { this->fTPCPion2Im2NegChNegEta = TPCPion2Im2NegChNegEta; }
  Double_t getTPCPion2Im2NegChNegEta() { return fTPCPion2Im2NegChNegEta; }
  void setTPCPion2MNegChNegEta( Double_t TPCPion2MNegChNegEta ) { this->fTPCPion2MNegChNegEta = TPCPion2MNegChNegEta; }
  Double_t getTPCPion2MNegChNegEta() { return fTPCPion2MNegChNegEta; }
  
  
  // TPC kaon
  void setTPCKaonRePosChPosEta( Double_t TPCKaonRePosChPosEta ) { this->fTPCKaonRePosChPosEta = TPCKaonRePosChPosEta; }
  Double_t getTPCKaonRePosChPosEta() { return fTPCKaonRePosChPosEta; }
  void setTPCKaonImPosChPosEta( Double_t TPCKaonImPosChPosEta ) { this->fTPCKaonImPosChPosEta = TPCKaonImPosChPosEta; }
  Double_t getTPCKaonImPosChPosEta() { return fTPCKaonImPosChPosEta; }
  void setTPCKaonMPosChPosEta( Double_t TPCKaonMPosChPosEta ) { this->fTPCKaonMPosChPosEta = TPCKaonMPosChPosEta; }
  Double_t getTPCKaonMPosChPosEta() { return fTPCKaonMPosChPosEta; }
  void setTPCKaonRePosChNegEta( Double_t TPCKaonRePosChNegEta ) { this->fTPCKaonRePosChNegEta = TPCKaonRePosChNegEta; }
  Double_t getTPCKaonRePosChNegEta() { return fTPCKaonRePosChNegEta; }
  void setTPCKaonImPosChNegEta( Double_t TPCKaonImPosChNegEta ) { this->fTPCKaonImPosChNegEta = TPCKaonImPosChNegEta; }
  Double_t getTPCKaonImPosChNegEta() { return fTPCKaonImPosChNegEta; }
  void setTPCKaonMPosChNegEta( Double_t TPCKaonMPosChNegEta ) { this->fTPCKaonMPosChNegEta = TPCKaonMPosChNegEta; }
  Double_t getTPCKaonMPosChNegEta() { return fTPCKaonMPosChNegEta; }
  void setTPCKaonReNegChPosEta( Double_t TPCKaonReNegChPosEta ) { this->fTPCKaonReNegChPosEta = TPCKaonReNegChPosEta; }
  Double_t getTPCKaonReNegChPosEta() { return fTPCKaonReNegChPosEta; }
  void setTPCKaonImNegChPosEta( Double_t TPCKaonImNegChPosEta ) { this->fTPCKaonImNegChPosEta = TPCKaonImNegChPosEta; }
  Double_t getTPCKaonImNegChPosEta() { return fTPCKaonImNegChPosEta; }
  void setTPCKaonMNegChPosEta( Double_t TPCKaonMNegChPosEta ) { this->fTPCKaonMNegChPosEta = TPCKaonMNegChPosEta; }
  Double_t getTPCKaonMNegChPosEta() { return fTPCKaonMNegChPosEta; }
  void setTPCKaonReNegChNegEta( Double_t TPCKaonReNegChNegEta ) { this->fTPCKaonReNegChNegEta = TPCKaonReNegChNegEta; }
  Double_t getTPCKaonReNegChNegEta() { return fTPCKaonReNegChNegEta; }
  void setTPCKaonImNegChNegEta( Double_t TPCKaonImNegChNegEta ) { this->fTPCKaonImNegChNegEta = TPCKaonImNegChNegEta; }
  Double_t getTPCKaonImNegChNegEta() { return fTPCKaonImNegChNegEta; }
  void setTPCKaonMNegChNegEta( Double_t TPCKaonMNegChNegEta ) { this->fTPCKaonMNegChNegEta = TPCKaonMNegChNegEta; }
  Double_t getTPCKaonMNegChNegEta() { return fTPCKaonMNegChNegEta; }

  void setTPCKaon2RePosChPosEta( Double_t TPCKaon2RePosChPosEta ) { this->fTPCKaon2RePosChPosEta = TPCKaon2RePosChPosEta; }
  Double_t getTPCKaon2RePosChPosEta() { return fTPCKaon2RePosChPosEta; }
  void setTPCKaon2ImPosChPosEta( Double_t TPCKaon2ImPosChPosEta ) { this->fTPCKaon2ImPosChPosEta = TPCKaon2ImPosChPosEta; }
  Double_t getTPCKaon2ImPosChPosEta() { return fTPCKaon2ImPosChPosEta; }
  void setTPCKaon2Re2PosChPosEta( Double_t TPCKaon2Re2PosChPosEta ) { this->fTPCKaon2Re2PosChPosEta = TPCKaon2Re2PosChPosEta; }
  Double_t getTPCKaon2Re2PosChPosEta() { return fTPCKaon2Re2PosChPosEta; }
  void setTPCKaon2Im2PosChPosEta( Double_t TPCKaon2Im2PosChPosEta ) { this->fTPCKaon2Im2PosChPosEta = TPCKaon2Im2PosChPosEta; }
  Double_t getTPCKaon2Im2PosChPosEta() { return fTPCKaon2Im2PosChPosEta; }
  void setTPCKaon2MPosChPosEta( Double_t TPCKaon2MPosChPosEta ) { this->fTPCKaon2MPosChPosEta = TPCKaon2MPosChPosEta; }
  Double_t getTPCKaon2MPosChPosEta() { return fTPCKaon2MPosChPosEta; }
  void setTPCKaon2RePosChNegEta( Double_t TPCKaon2RePosChNegEta ) { this->fTPCKaon2RePosChNegEta = TPCKaon2RePosChNegEta; }
  Double_t getTPCKaon2RePosChNegEta() { return fTPCKaon2RePosChNegEta; }
  void setTPCKaon2ImPosChNegEta( Double_t TPCKaon2ImPosChNegEta ) { this->fTPCKaon2ImPosChNegEta = TPCKaon2ImPosChNegEta; }
  Double_t getTPCKaon2ImPosChNegEta() { return fTPCKaon2ImPosChNegEta; }
  void setTPCKaon2Re2PosChNegEta( Double_t TPCKaon2Re2PosChNegEta ) { this->fTPCKaon2Re2PosChNegEta = TPCKaon2Re2PosChNegEta; }
  Double_t getTPCKaon2Re2PosChNegEta() { return fTPCKaon2Re2PosChNegEta; }
  void setTPCKaon2Im2PosChNegEta( Double_t TPCKaon2Im2PosChNegEta ) { this->fTPCKaon2Im2PosChNegEta = TPCKaon2Im2PosChNegEta; }
  Double_t getTPCKaon2Im2PosChNegEta() { return fTPCKaon2Im2PosChNegEta; }
  void setTPCKaon2MPosChNegEta( Double_t TPCKaon2MPosChNegEta ) { this->fTPCKaon2MPosChNegEta = TPCKaon2MPosChNegEta; }
  Double_t getTPCKaon2MPosChNegEta() { return fTPCKaon2MPosChNegEta; }
  void setTPCKaon2ReNegChPosEta( Double_t TPCKaon2ReNegChPosEta ) { this->fTPCKaon2ReNegChPosEta = TPCKaon2ReNegChPosEta; }
  Double_t getTPCKaon2ReNegChPosEta() { return fTPCKaon2ReNegChPosEta; }
  void setTPCKaon2ImNegChPosEta( Double_t TPCKaon2ImNegChPosEta ) { this->fTPCKaon2ImNegChPosEta = TPCKaon2ImNegChPosEta; }
  Double_t getTPCKaon2ImNegChPosEta() { return fTPCKaon2ImNegChPosEta; }
  void setTPCKaon2Re2NegChPosEta( Double_t TPCKaon2Re2NegChPosEta ) { this->fTPCKaon2Re2NegChPosEta = TPCKaon2Re2NegChPosEta; }
  Double_t getTPCKaon2Re2NegChPosEta() { return fTPCKaon2Re2NegChPosEta; }
  void setTPCKaon2Im2NegChPosEta( Double_t TPCKaon2Im2NegChPosEta ) { this->fTPCKaon2Im2NegChPosEta = TPCKaon2Im2NegChPosEta; }
  Double_t getTPCKaon2Im2NegChPosEta() { return fTPCKaon2Im2NegChPosEta; }
  void setTPCKaon2MNegChPosEta( Double_t TPCKaon2MNegChPosEta ) { this->fTPCKaon2MNegChPosEta = TPCKaon2MNegChPosEta; }
  Double_t getTPCKaon2MNegChPosEta() { return fTPCKaon2MNegChPosEta; }
  void setTPCKaon2ReNegChNegEta( Double_t TPCKaon2ReNegChNegEta ) { this->fTPCKaon2ReNegChNegEta = TPCKaon2ReNegChNegEta; }
  Double_t getTPCKaon2ReNegChNegEta() { return fTPCKaon2ReNegChNegEta; }
  void setTPCKaon2ImNegChNegEta( Double_t TPCKaon2ImNegChNegEta ) { this->fTPCKaon2ImNegChNegEta = TPCKaon2ImNegChNegEta; }
  Double_t getTPCKaon2ImNegChNegEta() { return fTPCKaon2ImNegChNegEta; }
  void setTPCKaon2Re2NegChNegEta( Double_t TPCKaon2Re2NegChNegEta ) { this->fTPCKaon2Re2NegChNegEta = TPCKaon2Re2NegChNegEta; }
  Double_t getTPCKaon2Re2NegChNegEta() { return fTPCKaon2Re2NegChNegEta; }
  void setTPCKaon2Im2NegChNegEta( Double_t TPCKaon2Im2NegChNegEta ) { this->fTPCKaon2Im2NegChNegEta = TPCKaon2Im2NegChNegEta; }
  Double_t getTPCKaon2Im2NegChNegEta() { return fTPCKaon2Im2NegChNegEta; }
  void setTPCKaon2MNegChNegEta( Double_t TPCKaon2MNegChNegEta ) { this->fTPCKaon2MNegChNegEta = TPCKaon2MNegChNegEta; }
  Double_t getTPCKaon2MNegChNegEta() { return fTPCKaon2MNegChNegEta; }
  
  
  // TPC proton
  void setTPCProtonRePosChPosEta( Double_t TPCProtonRePosChPosEta ) { this->fTPCProtonRePosChPosEta = TPCProtonRePosChPosEta; }
  Double_t getTPCProtonRePosChPosEta() { return fTPCProtonRePosChPosEta; }
  void setTPCProtonImPosChPosEta( Double_t TPCProtonImPosChPosEta ) { this->fTPCProtonImPosChPosEta = TPCProtonImPosChPosEta; }
  Double_t getTPCProtonImPosChPosEta() { return fTPCProtonImPosChPosEta; }
  void setTPCProtonMPosChPosEta( Double_t TPCProtonMPosChPosEta ) { this->fTPCProtonMPosChPosEta = TPCProtonMPosChPosEta; }
  Double_t getTPCProtonMPosChPosEta() { return fTPCProtonMPosChPosEta; }
  void setTPCProtonRePosChNegEta( Double_t TPCProtonRePosChNegEta ) { this->fTPCProtonRePosChNegEta = TPCProtonRePosChNegEta; }
  Double_t getTPCProtonRePosChNegEta() { return fTPCProtonRePosChNegEta; }
  void setTPCProtonImPosChNegEta( Double_t TPCProtonImPosChNegEta ) { this->fTPCProtonImPosChNegEta = TPCProtonImPosChNegEta; }
  Double_t getTPCProtonImPosChNegEta() { return fTPCProtonImPosChNegEta; }
  void setTPCProtonMPosChNegEta( Double_t TPCProtonMPosChNegEta ) { this->fTPCProtonMPosChNegEta = TPCProtonMPosChNegEta; }
  Double_t getTPCProtonMPosChNegEta() { return fTPCProtonMPosChNegEta; }
  void setTPCProtonReNegChPosEta( Double_t TPCProtonReNegChPosEta ) { this->fTPCProtonReNegChPosEta = TPCProtonReNegChPosEta; }
  Double_t getTPCProtonReNegChPosEta() { return fTPCProtonReNegChPosEta; }
  void setTPCProtonImNegChPosEta( Double_t TPCProtonImNegChPosEta ) { this->fTPCProtonImNegChPosEta = TPCProtonImNegChPosEta; }
  Double_t getTPCProtonImNegChPosEta() { return fTPCProtonImNegChPosEta; }
  void setTPCProtonMNegChPosEta( Double_t TPCProtonMNegChPosEta ) { this->fTPCProtonMNegChPosEta = TPCProtonMNegChPosEta; }
  Double_t getTPCProtonMNegChPosEta() { return fTPCProtonMNegChPosEta; }
  void setTPCProtonReNegChNegEta( Double_t TPCProtonReNegChNegEta ) { this->fTPCProtonReNegChNegEta = TPCProtonReNegChNegEta; }
  Double_t getTPCProtonReNegChNegEta() { return fTPCProtonReNegChNegEta; }
  void setTPCProtonImNegChNegEta( Double_t TPCProtonImNegChNegEta ) { this->fTPCProtonImNegChNegEta = TPCProtonImNegChNegEta; }
  Double_t getTPCProtonImNegChNegEta() { return fTPCProtonImNegChNegEta; }
  void setTPCProtonMNegChNegEta( Double_t TPCProtonMNegChNegEta ) { this->fTPCProtonMNegChNegEta = TPCProtonMNegChNegEta; }
  Double_t getTPCProtonMNegChNegEta() { return fTPCProtonMNegChNegEta; }

  void setTPCProton2RePosChPosEta( Double_t TPCProton2RePosChPosEta ) { this->fTPCProton2RePosChPosEta = TPCProton2RePosChPosEta; }
  Double_t getTPCProton2RePosChPosEta() { return fTPCProton2RePosChPosEta; }
  void setTPCProton2ImPosChPosEta( Double_t TPCProton2ImPosChPosEta ) { this->fTPCProton2ImPosChPosEta = TPCProton2ImPosChPosEta; }
  Double_t getTPCProton2ImPosChPosEta() { return fTPCProton2ImPosChPosEta; }
  void setTPCProton2Re2PosChPosEta( Double_t TPCProton2Re2PosChPosEta ) { this->fTPCProton2Re2PosChPosEta = TPCProton2Re2PosChPosEta; }
  Double_t getTPCProton2Re2PosChPosEta() { return fTPCProton2Re2PosChPosEta; }
  void setTPCProton2Im2PosChPosEta( Double_t TPCProton2Im2PosChPosEta ) { this->fTPCProton2Im2PosChPosEta = TPCProton2Im2PosChPosEta; }
  Double_t getTPCProton2Im2PosChPosEta() { return fTPCProton2Im2PosChPosEta; }
  void setTPCProton2MPosChPosEta( Double_t TPCProton2MPosChPosEta ) { this->fTPCProton2MPosChPosEta = TPCProton2MPosChPosEta; }
  Double_t getTPCProton2MPosChPosEta() { return fTPCProton2MPosChPosEta; }
  void setTPCProton2RePosChNegEta( Double_t TPCProton2RePosChNegEta ) { this->fTPCProton2RePosChNegEta = TPCProton2RePosChNegEta; }
  Double_t getTPCProton2RePosChNegEta() { return fTPCProton2RePosChNegEta; }
  void setTPCProton2ImPosChNegEta( Double_t TPCProton2ImPosChNegEta ) { this->fTPCProton2ImPosChNegEta = TPCProton2ImPosChNegEta; }
  Double_t getTPCProton2ImPosChNegEta() { return fTPCProton2ImPosChNegEta; }
  void setTPCProton2Re2PosChNegEta( Double_t TPCProton2Re2PosChNegEta ) { this->fTPCProton2Re2PosChNegEta = TPCProton2Re2PosChNegEta; }
  Double_t getTPCProton2Re2PosChNegEta() { return fTPCProton2Re2PosChNegEta; }
  void setTPCProton2Im2PosChNegEta( Double_t TPCProton2Im2PosChNegEta ) { this->fTPCProton2Im2PosChNegEta = TPCProton2Im2PosChNegEta; }
  Double_t getTPCProton2Im2PosChNegEta() { return fTPCProton2Im2PosChNegEta; }
  void setTPCProton2MPosChNegEta( Double_t TPCProton2MPosChNegEta ) { this->fTPCProton2MPosChNegEta = TPCProton2MPosChNegEta; }
  Double_t getTPCProton2MPosChNegEta() { return fTPCProton2MPosChNegEta; }
  void setTPCProton2ReNegChPosEta( Double_t TPCProton2ReNegChPosEta ) { this->fTPCProton2ReNegChPosEta = TPCProton2ReNegChPosEta; }
  Double_t getTPCProton2ReNegChPosEta() { return fTPCProton2ReNegChPosEta; }
  void setTPCProton2ImNegChPosEta( Double_t TPCProton2ImNegChPosEta ) { this->fTPCProton2ImNegChPosEta = TPCProton2ImNegChPosEta; }
  Double_t getTPCProton2ImNegChPosEta() { return fTPCProton2ImNegChPosEta; }
  void setTPCProton2Re2NegChPosEta( Double_t TPCProton2Re2NegChPosEta ) { this->fTPCProton2Re2NegChPosEta = TPCProton2Re2NegChPosEta; }
  Double_t getTPCProton2Re2NegChPosEta() { return fTPCProton2Re2NegChPosEta; }
  void setTPCProton2Im2NegChPosEta( Double_t TPCProton2Im2NegChPosEta ) { this->fTPCProton2Im2NegChPosEta = TPCProton2Im2NegChPosEta; }
  Double_t getTPCProton2Im2NegChPosEta() { return fTPCProton2Im2NegChPosEta; }
  void setTPCProton2MNegChPosEta( Double_t TPCProton2MNegChPosEta ) { this->fTPCProton2MNegChPosEta = TPCProton2MNegChPosEta; }
  Double_t getTPCProton2MNegChPosEta() { return fTPCProton2MNegChPosEta; }
  void setTPCProton2ReNegChNegEta( Double_t TPCProton2ReNegChNegEta ) { this->fTPCProton2ReNegChNegEta = TPCProton2ReNegChNegEta; }
  Double_t getTPCProton2ReNegChNegEta() { return fTPCProton2ReNegChNegEta; }
  void setTPCProton2ImNegChNegEta( Double_t TPCProton2ImNegChNegEta ) { this->fTPCProton2ImNegChNegEta = TPCProton2ImNegChNegEta; }
  Double_t getTPCProton2ImNegChNegEta() { return fTPCProton2ImNegChNegEta; }
  void setTPCProton2Re2NegChNegEta( Double_t TPCProton2Re2NegChNegEta ) { this->fTPCProton2Re2NegChNegEta = TPCProton2Re2NegChNegEta; }
  Double_t getTPCProton2Re2NegChNegEta() { return fTPCProton2Re2NegChNegEta; }
  void setTPCProton2Im2NegChNegEta( Double_t TPCProton2Im2NegChNegEta ) { this->fTPCProton2Im2NegChNegEta = TPCProton2Im2NegChNegEta; }
  Double_t getTPCProton2Im2NegChNegEta() { return fTPCProton2Im2NegChNegEta; }
  void setTPCProton2MNegChNegEta( Double_t TPCProton2MNegChNegEta ) { this->fTPCProton2MNegChNegEta = TPCProton2MNegChNegEta; }
  Double_t getTPCProton2MNegChNegEta() { return fTPCProton2MNegChNegEta; }
  
  
  
 private:
  Int_t fRunNum;
  Double_t fCentrality;
  Double_t fVtxPosX;
  Double_t fVtxPosY;
  Double_t fVtxPosZ;
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
  
  // TPC Pion (|eta|<0.8)
  Double_t fTPCPionRePosChPosEta; // w * cos(theta+) eta+
  Double_t fTPCPionImPosChPosEta; // w * sin(theta+) eta+
  Double_t fTPCPionMPosChPosEta;   // w eta+
  Double_t fTPCPionRePosChNegEta; // w * cos(theta+) eta-
  Double_t fTPCPionImPosChNegEta; // w * sin(theta+) eta-
  Double_t fTPCPionMPosChNegEta;    // w eta-
  Double_t fTPCPionReNegChPosEta; // w * cos(theta-) eta+
  Double_t fTPCPionImNegChPosEta; // w * sin(theta-) eta+
  Double_t fTPCPionMNegChPosEta;    // w eta+
  Double_t fTPCPionReNegChNegEta; // w * cos(theta-) eta-
  Double_t fTPCPionImNegChNegEta; // w * sin(theta-) eta-
  Double_t fTPCPionMNegChNegEta;    // w eta-
  
  Double_t fTPCPion2RePosChPosEta; // w * cos(2theta+) eta+
  Double_t fTPCPion2ImPosChPosEta; // w * sin(2theta+) eta+
  Double_t fTPCPion2Re2PosChPosEta; // w^2 * cos(2theta+) eta+
  Double_t fTPCPion2Im2PosChPosEta; // w^2 * sin(2theta+) eta+
  Double_t fTPCPion2MPosChPosEta;   // w^2 eta+
  Double_t fTPCPion2RePosChNegEta; // w * cos(2theta+) eta-
  Double_t fTPCPion2ImPosChNegEta; // w * sin(2theta+) eta-
  Double_t fTPCPion2Re2PosChNegEta; // w^2 * cos(2theta+) eta-
  Double_t fTPCPion2Im2PosChNegEta; // w^2 * sin(2theta+) eta-
  Double_t fTPCPion2MPosChNegEta;    // w^2 eta-
  Double_t fTPCPion2ReNegChPosEta; // w * cos(2theta-) eta+
  Double_t fTPCPion2ImNegChPosEta; // w * sin(2theta-) eta+
  Double_t fTPCPion2Re2NegChPosEta; // w^2 * cos(2theta-) eta+
  Double_t fTPCPion2Im2NegChPosEta; // w^2 * sin(2theta-) eta+
  Double_t fTPCPion2MNegChPosEta;    // w^2 eta+
  Double_t fTPCPion2ReNegChNegEta; // w * cos(2theta-) eta-
  Double_t fTPCPion2ImNegChNegEta; // w * sin(2theta-) eta-
  Double_t fTPCPion2Re2NegChNegEta; // w^2 * cos(2theta-) eta-
  Double_t fTPCPion2Im2NegChNegEta; // w^2 * sin(2theta-) eta-
  Double_t fTPCPion2MNegChNegEta;    // w^2 eta-
  
  // TPC Kaon (|eta|<0.8)
  Double_t fTPCKaonRePosChPosEta; // w * cos(theta+) eta+
  Double_t fTPCKaonImPosChPosEta; // w * sin(theta+) eta+
  Double_t fTPCKaonMPosChPosEta;   // w eta+
  Double_t fTPCKaonRePosChNegEta; // w * cos(theta+) eta-
  Double_t fTPCKaonImPosChNegEta; // w * sin(theta+) eta-
  Double_t fTPCKaonMPosChNegEta;    // w eta-
  Double_t fTPCKaonReNegChPosEta; // w * cos(theta-) eta+
  Double_t fTPCKaonImNegChPosEta; // w * sin(theta-) eta+
  Double_t fTPCKaonMNegChPosEta;    // w eta+
  Double_t fTPCKaonReNegChNegEta; // w * cos(theta-) eta-
  Double_t fTPCKaonImNegChNegEta; // w * sin(theta-) eta-
  Double_t fTPCKaonMNegChNegEta;    // w eta-
  
  Double_t fTPCKaon2RePosChPosEta; // w * cos(2theta+) eta+
  Double_t fTPCKaon2ImPosChPosEta; // w * sin(2theta+) eta+
  Double_t fTPCKaon2Re2PosChPosEta; // w^2 * cos(2theta+) eta+
  Double_t fTPCKaon2Im2PosChPosEta; // w^2 * sin(2theta+) eta+
  Double_t fTPCKaon2MPosChPosEta;   // w^2 eta+
  Double_t fTPCKaon2RePosChNegEta; // w * cos(2theta+) eta-
  Double_t fTPCKaon2ImPosChNegEta; // w * sin(2theta+) eta-
  Double_t fTPCKaon2Re2PosChNegEta; // w^2 * cos(2theta+) eta-
  Double_t fTPCKaon2Im2PosChNegEta; // w^2 * sin(2theta+) eta-
  Double_t fTPCKaon2MPosChNegEta;    // w^2 eta-
  Double_t fTPCKaon2ReNegChPosEta; // w * cos(2theta-) eta+
  Double_t fTPCKaon2ImNegChPosEta; // w * sin(2theta-) eta+
  Double_t fTPCKaon2Re2NegChPosEta; // w^2 * cos(2theta-) eta+
  Double_t fTPCKaon2Im2NegChPosEta; // w^2 * sin(2theta-) eta+
  Double_t fTPCKaon2MNegChPosEta;    // w^2 eta+
  Double_t fTPCKaon2ReNegChNegEta; // w * cos(2theta-) eta-
  Double_t fTPCKaon2ImNegChNegEta; // w * sin(2theta-) eta-
  Double_t fTPCKaon2Re2NegChNegEta; // w^2 * cos(2theta-) eta-
  Double_t fTPCKaon2Im2NegChNegEta; // w^2 * sin(2theta-) eta-
  Double_t fTPCKaon2MNegChNegEta;    // w^2 eta-
  
  // TPC Proton (|eta|<0.8)
  Double_t fTPCProtonRePosChPosEta; // w * cos(theta+) eta+
  Double_t fTPCProtonImPosChPosEta; // w * sin(theta+) eta+
  Double_t fTPCProtonMPosChPosEta;   // w eta+
  Double_t fTPCProtonRePosChNegEta; // w * cos(theta+) eta-
  Double_t fTPCProtonImPosChNegEta; // w * sin(theta+) eta-
  Double_t fTPCProtonMPosChNegEta;    // w eta-
  Double_t fTPCProtonReNegChPosEta; // w * cos(theta-) eta+
  Double_t fTPCProtonImNegChPosEta; // w * sin(theta-) eta+
  Double_t fTPCProtonMNegChPosEta;    // w eta+
  Double_t fTPCProtonReNegChNegEta; // w * cos(theta-) eta-
  Double_t fTPCProtonImNegChNegEta; // w * sin(theta-) eta-
  Double_t fTPCProtonMNegChNegEta;    // w eta-
  
  Double_t fTPCProton2RePosChPosEta; // w * cos(2theta+) eta+
  Double_t fTPCProton2ImPosChPosEta; // w * sin(2theta+) eta+
  Double_t fTPCProton2Re2PosChPosEta; // w^2 * cos(2theta+) eta+
  Double_t fTPCProton2Im2PosChPosEta; // w^2 * sin(2theta+) eta+
  Double_t fTPCProton2MPosChPosEta;   // w^2 eta+
  Double_t fTPCProton2RePosChNegEta; // w * cos(2theta+) eta-
  Double_t fTPCProton2ImPosChNegEta; // w * sin(2theta+) eta-
  Double_t fTPCProton2Re2PosChNegEta; // w^2 * cos(2theta+) eta-
  Double_t fTPCProton2Im2PosChNegEta; // w^2 * sin(2theta+) eta-
  Double_t fTPCProton2MPosChNegEta;    // w^2 eta-
  Double_t fTPCProton2ReNegChPosEta; // w * cos(2theta-) eta+
  Double_t fTPCProton2ImNegChPosEta; // w * sin(2theta-) eta+
  Double_t fTPCProton2Re2NegChPosEta; // w^2 * cos(2theta-) eta+
  Double_t fTPCProton2Im2NegChPosEta; // w^2 * sin(2theta-) eta+
  Double_t fTPCProton2MNegChPosEta;    // w^2 eta+
  Double_t fTPCProton2ReNegChNegEta; // w * cos(2theta-) eta-
  Double_t fTPCProton2ImNegChNegEta; // w * sin(2theta-) eta-
  Double_t fTPCProton2Re2NegChNegEta; // w^2 * cos(2theta-) eta-
  Double_t fTPCProton2Im2NegChNegEta; // w^2 * sin(2theta-) eta-
  Double_t fTPCProton2MNegChNegEta;    // w^2 eta-
  
  ClassDef(AliAnalysisTaskGammaDeltaPIDSaveQvecEvent, 1);
};

#endif
