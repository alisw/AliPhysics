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
  
  void setRawPeriod(UInt_t period) {this->fRawPeriod = period;}
  UInt_t getRawPeriod() { return fRawPeriod; }
  void setRawOrbitNumber24(UInt_t orbit24) {this->fRawOrbitNumber24 = orbit24;}
  UInt_t getRawOrbitNumber24() { return fRawOrbitNumber24; }
  void setOrbitNumber(UInt_t orbit) {this->fOrbitNumber = orbit;}
  UInt_t getOrbitNumber() { return fOrbitNumber; }
  void setBunchCrossNumber(UShort_t bunchCrossNumber) {this->fBunchCrossNumber = bunchCrossNumber;}
  UShort_t getBunchCrossNumber() { return fBunchCrossNumber; }
  void setTimeStamp(UInt_t timeStamp) {this->fTimeStamp = timeStamp;}
  UInt_t getTimeStamp() { return fTimeStamp; }
  
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
  
  void setTowV0Craw0( Double_t TowV0Craw0 ) { this->fTowV0Craw0 = TowV0Craw0; }
  Double_t getTowV0Craw0() { return fTowV0Craw0; }
  void setTowV0Craw1( Double_t TowV0Craw1 ) { this->fTowV0Craw1 = TowV0Craw1; }
  Double_t getTowV0Craw1() { return fTowV0Craw1; }
  void setTowV0Craw2( Double_t TowV0Craw2 ) { this->fTowV0Craw2 = TowV0Craw2; }
  Double_t getTowV0Craw2() { return fTowV0Craw2; }
  void setTowV0Craw3( Double_t TowV0Craw3 ) { this->fTowV0Craw3 = TowV0Craw3; }
  Double_t getTowV0Craw3() { return fTowV0Craw3; }
  void setTowV0Craw4( Double_t TowV0Craw4 ) { this->fTowV0Craw4 = TowV0Craw4; }
  Double_t getTowV0Craw4() { return fTowV0Craw4; }
  void setTowV0Craw5( Double_t TowV0Craw5 ) { this->fTowV0Craw5 = TowV0Craw5; }
  Double_t getTowV0Craw5() { return fTowV0Craw5; }
  void setTowV0Craw6( Double_t TowV0Craw6 ) { this->fTowV0Craw6 = TowV0Craw6; }
  Double_t getTowV0Craw6() { return fTowV0Craw6; }
  void setTowV0Craw7( Double_t TowV0Craw7 ) { this->fTowV0Craw7 = TowV0Craw7; }
  Double_t getTowV0Craw7() { return fTowV0Craw7; }
  void setTowV0Craw8( Double_t TowV0Craw8 ) { this->fTowV0Craw8 = TowV0Craw8; }
  Double_t getTowV0Craw8() { return fTowV0Craw8; }
  void setTowV0Craw9( Double_t TowV0Craw9 ) { this->fTowV0Craw9 = TowV0Craw9; }
  Double_t getTowV0Craw9() { return fTowV0Craw9; }
  void setTowV0Craw10( Double_t TowV0Craw10 ) { this->fTowV0Craw10 = TowV0Craw10; }
  Double_t getTowV0Craw10() { return fTowV0Craw10; }
  void setTowV0Craw11( Double_t TowV0Craw11 ) { this->fTowV0Craw11 = TowV0Craw11; }
  Double_t getTowV0Craw11() { return fTowV0Craw11; }
  void setTowV0Craw12( Double_t TowV0Craw12 ) { this->fTowV0Craw12 = TowV0Craw12; }
  Double_t getTowV0Craw12() { return fTowV0Craw12; }
  void setTowV0Craw13( Double_t TowV0Craw13 ) { this->fTowV0Craw13 = TowV0Craw13; }
  Double_t getTowV0Craw13() { return fTowV0Craw13; }
  void setTowV0Craw14( Double_t TowV0Craw14 ) { this->fTowV0Craw14 = TowV0Craw14; }
  Double_t getTowV0Craw14() { return fTowV0Craw14; }
  void setTowV0Craw15( Double_t TowV0Craw15 ) { this->fTowV0Craw15 = TowV0Craw15; }
  Double_t getTowV0Craw15() { return fTowV0Craw15; }
  void setTowV0Craw16( Double_t TowV0Craw16 ) { this->fTowV0Craw16 = TowV0Craw16; }
  Double_t getTowV0Craw16() { return fTowV0Craw16; }
  void setTowV0Craw17( Double_t TowV0Craw17 ) { this->fTowV0Craw17 = TowV0Craw17; }
  Double_t getTowV0Craw17() { return fTowV0Craw17; }
  void setTowV0Craw18( Double_t TowV0Craw18 ) { this->fTowV0Craw18 = TowV0Craw18; }
  Double_t getTowV0Craw18() { return fTowV0Craw18; }
  void setTowV0Craw19( Double_t TowV0Craw19 ) { this->fTowV0Craw19 = TowV0Craw19; }
  Double_t getTowV0Craw19() { return fTowV0Craw19; }
  void setTowV0Craw20( Double_t TowV0Craw20 ) { this->fTowV0Craw20 = TowV0Craw20; }
  Double_t getTowV0Craw20() { return fTowV0Craw20; }
  void setTowV0Craw21( Double_t TowV0Craw21 ) { this->fTowV0Craw21 = TowV0Craw21; }
  Double_t getTowV0Craw21() { return fTowV0Craw21; }
  void setTowV0Craw22( Double_t TowV0Craw22 ) { this->fTowV0Craw22 = TowV0Craw22; }
  Double_t getTowV0Craw22() { return fTowV0Craw22; }
  void setTowV0Craw23( Double_t TowV0Craw23 ) { this->fTowV0Craw23 = TowV0Craw23; }
  Double_t getTowV0Craw23() { return fTowV0Craw23; }
  void setTowV0Craw24( Double_t TowV0Craw24 ) { this->fTowV0Craw24 = TowV0Craw24; }
  Double_t getTowV0Craw24() { return fTowV0Craw24; }
  void setTowV0Craw25( Double_t TowV0Craw25 ) { this->fTowV0Craw25 = TowV0Craw25; }
  Double_t getTowV0Craw25() { return fTowV0Craw25; }
  void setTowV0Craw26( Double_t TowV0Craw26 ) { this->fTowV0Craw26 = TowV0Craw26; }
  Double_t getTowV0Craw26() { return fTowV0Craw26; }
  void setTowV0Craw27( Double_t TowV0Craw27 ) { this->fTowV0Craw27 = TowV0Craw27; }
  Double_t getTowV0Craw27() { return fTowV0Craw27; }
  void setTowV0Craw28( Double_t TowV0Craw28 ) { this->fTowV0Craw28 = TowV0Craw28; }
  Double_t getTowV0Craw28() { return fTowV0Craw28; }
  void setTowV0Craw29( Double_t TowV0Craw29 ) { this->fTowV0Craw29 = TowV0Craw29; }
  Double_t getTowV0Craw29() { return fTowV0Craw29; }
  void setTowV0Craw30( Double_t TowV0Craw30 ) { this->fTowV0Craw30 = TowV0Craw30; }
  Double_t getTowV0Craw30() { return fTowV0Craw30; }
  void setTowV0Craw31( Double_t TowV0Craw31 ) { this->fTowV0Craw31 = TowV0Craw31; }
  Double_t getTowV0Craw31() { return fTowV0Craw31; }

  void setTowV0Araw0( Double_t TowV0Araw0 ) { this->fTowV0Araw0 = TowV0Araw0; }
  Double_t getTowV0Araw0() { return fTowV0Araw0; }
  void setTowV0Araw1( Double_t TowV0Araw1 ) { this->fTowV0Araw1 = TowV0Araw1; }
  Double_t getTowV0Araw1() { return fTowV0Araw1; }
  void setTowV0Araw2( Double_t TowV0Araw2 ) { this->fTowV0Araw2 = TowV0Araw2; }
  Double_t getTowV0Araw2() { return fTowV0Araw2; }
  void setTowV0Araw3( Double_t TowV0Araw3 ) { this->fTowV0Araw3 = TowV0Araw3; }
  Double_t getTowV0Araw3() { return fTowV0Araw3; }
  void setTowV0Araw4( Double_t TowV0Araw4 ) { this->fTowV0Araw4 = TowV0Araw4; }
  Double_t getTowV0Araw4() { return fTowV0Araw4; }
  void setTowV0Araw5( Double_t TowV0Araw5 ) { this->fTowV0Araw5 = TowV0Araw5; }
  Double_t getTowV0Araw5() { return fTowV0Araw5; }
  void setTowV0Araw6( Double_t TowV0Araw6 ) { this->fTowV0Araw6 = TowV0Araw6; }
  Double_t getTowV0Araw6() { return fTowV0Araw6; }
  void setTowV0Araw7( Double_t TowV0Araw7 ) { this->fTowV0Araw7 = TowV0Araw7; }
  Double_t getTowV0Araw7() { return fTowV0Araw7; }
  void setTowV0Araw8( Double_t TowV0Araw8 ) { this->fTowV0Araw8 = TowV0Araw8; }
  Double_t getTowV0Araw8() { return fTowV0Araw8; }
  void setTowV0Araw9( Double_t TowV0Araw9 ) { this->fTowV0Araw9 = TowV0Araw9; }
  Double_t getTowV0Araw9() { return fTowV0Araw9; }
  void setTowV0Araw10( Double_t TowV0Araw10 ) { this->fTowV0Araw10 = TowV0Araw10; }
  Double_t getTowV0Araw10() { return fTowV0Araw10; }
  void setTowV0Araw11( Double_t TowV0Araw11 ) { this->fTowV0Araw11 = TowV0Araw11; }
  Double_t getTowV0Araw11() { return fTowV0Araw11; }
  void setTowV0Araw12( Double_t TowV0Araw12 ) { this->fTowV0Araw12 = TowV0Araw12; }
  Double_t getTowV0Araw12() { return fTowV0Araw12; }
  void setTowV0Araw13( Double_t TowV0Araw13 ) { this->fTowV0Araw13 = TowV0Araw13; }
  Double_t getTowV0Araw13() { return fTowV0Araw13; }
  void setTowV0Araw14( Double_t TowV0Araw14 ) { this->fTowV0Araw14 = TowV0Araw14; }
  Double_t getTowV0Araw14() { return fTowV0Araw14; }
  void setTowV0Araw15( Double_t TowV0Araw15 ) { this->fTowV0Araw15 = TowV0Araw15; }
  Double_t getTowV0Araw15() { return fTowV0Araw15; }
  void setTowV0Araw16( Double_t TowV0Araw16 ) { this->fTowV0Araw16 = TowV0Araw16; }
  Double_t getTowV0Araw16() { return fTowV0Araw16; }
  void setTowV0Araw17( Double_t TowV0Araw17 ) { this->fTowV0Araw17 = TowV0Araw17; }
  Double_t getTowV0Araw17() { return fTowV0Araw17; }
  void setTowV0Araw18( Double_t TowV0Araw18 ) { this->fTowV0Araw18 = TowV0Araw18; }
  Double_t getTowV0Araw18() { return fTowV0Araw18; }
  void setTowV0Araw19( Double_t TowV0Araw19 ) { this->fTowV0Araw19 = TowV0Araw19; }
  Double_t getTowV0Araw19() { return fTowV0Araw19; }
  void setTowV0Araw20( Double_t TowV0Araw20 ) { this->fTowV0Araw20 = TowV0Araw20; }
  Double_t getTowV0Araw20() { return fTowV0Araw20; }
  void setTowV0Araw21( Double_t TowV0Araw21 ) { this->fTowV0Araw21 = TowV0Araw21; }
  Double_t getTowV0Araw21() { return fTowV0Araw21; }
  void setTowV0Araw22( Double_t TowV0Araw22 ) { this->fTowV0Araw22 = TowV0Araw22; }
  Double_t getTowV0Araw22() { return fTowV0Araw22; }
  void setTowV0Araw23( Double_t TowV0Araw23 ) { this->fTowV0Araw23 = TowV0Araw23; }
  Double_t getTowV0Araw23() { return fTowV0Araw23; }
  void setTowV0Araw24( Double_t TowV0Araw24 ) { this->fTowV0Araw24 = TowV0Araw24; }
  Double_t getTowV0Araw24() { return fTowV0Araw24; }
  void setTowV0Araw25( Double_t TowV0Araw25 ) { this->fTowV0Araw25 = TowV0Araw25; }
  Double_t getTowV0Araw25() { return fTowV0Araw25; }
  void setTowV0Araw26( Double_t TowV0Araw26 ) { this->fTowV0Araw26 = TowV0Araw26; }
  Double_t getTowV0Araw26() { return fTowV0Araw26; }
  void setTowV0Araw27( Double_t TowV0Araw27 ) { this->fTowV0Araw27 = TowV0Araw27; }
  Double_t getTowV0Araw27() { return fTowV0Araw27; }
  void setTowV0Araw28( Double_t TowV0Araw28 ) { this->fTowV0Araw28 = TowV0Araw28; }
  Double_t getTowV0Araw28() { return fTowV0Araw28; }
  void setTowV0Araw29( Double_t TowV0Araw29 ) { this->fTowV0Araw29 = TowV0Araw29; }
  Double_t getTowV0Araw29() { return fTowV0Araw29; }
  void setTowV0Araw30( Double_t TowV0Araw30 ) { this->fTowV0Araw30 = TowV0Araw30; }
  Double_t getTowV0Araw30() { return fTowV0Araw30; }
  void setTowV0Araw31( Double_t TowV0Araw31 ) { this->fTowV0Araw31 = TowV0Araw31; }
  Double_t getTowV0Araw31() { return fTowV0Araw31; }
  
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
  
  void setTPC4Re2PosChPosEta( Double_t TPC4Re2PosChPosEta ) { this->fTPC4Re2PosChPosEta = TPC4Re2PosChPosEta; }
  Double_t getTPC4Re2PosChPosEta() { return fTPC4Re2PosChPosEta; }
  void setTPC4Im2PosChPosEta( Double_t TPC4Im2PosChPosEta ) { this->fTPC4Im2PosChPosEta = TPC4Im2PosChPosEta; }
  Double_t getTPC4Im2PosChPosEta() { return fTPC4Im2PosChPosEta; }
  void setTPC2Re3PosChPosEta( Double_t TPC2Re3PosChPosEta ) { this->fTPC2Re3PosChPosEta = TPC2Re3PosChPosEta; }
  Double_t getTPC2Re3PosChPosEta() { return fTPC2Re3PosChPosEta; }
  void setTPC2Im3PosChPosEta( Double_t TPC2Im3PosChPosEta ) { this->fTPC2Im3PosChPosEta = TPC2Im3PosChPosEta; }
  Double_t getTPC2Im3PosChPosEta() { return fTPC2Im3PosChPosEta; }
  void setTPC0MPosChPosEta( Double_t TPC0MPosChPosEta ) { this->fTPC0MPosChPosEta = TPC0MPosChPosEta; }
  Double_t getTPC0MPosChPosEta() { return fTPC0MPosChPosEta; }
  void setTPC3MPosChPosEta( Double_t TPC3MPosChPosEta ) { this->fTPC3MPosChPosEta = TPC3MPosChPosEta; }
  Double_t getTPC3MPosChPosEta() { return fTPC3MPosChPosEta; }
  void setTPC4MPosChPosEta( Double_t TPC4MPosChPosEta ) { this->fTPC4MPosChPosEta = TPC4MPosChPosEta; }
  Double_t getTPC4MPosChPosEta() { return fTPC4MPosChPosEta; }
  void setTPC4Re2PosChNegEta( Double_t TPC4Re2PosChNegEta ) { this->fTPC4Re2PosChNegEta = TPC4Re2PosChNegEta; }
  Double_t getTPC4Re2PosChNegEta() { return fTPC4Re2PosChNegEta; }
  void setTPC4Im2PosChNegEta( Double_t TPC4Im2PosChNegEta ) { this->fTPC4Im2PosChNegEta = TPC4Im2PosChNegEta; }
  Double_t getTPC4Im2PosChNegEta() { return fTPC4Im2PosChNegEta; }
  void setTPC2Re3PosChNegEta( Double_t TPC2Re3PosChNegEta ) { this->fTPC2Re3PosChNegEta = TPC2Re3PosChNegEta; }
  Double_t getTPC2Re3PosChNegEta() { return fTPC2Re3PosChNegEta; }
  void setTPC2Im3PosChNegEta( Double_t TPC2Im3PosChNegEta ) { this->fTPC2Im3PosChNegEta = TPC2Im3PosChNegEta; }
  Double_t getTPC2Im3PosChNegEta() { return fTPC2Im3PosChNegEta; }
  void setTPC0MPosChNegEta( Double_t TPC0MPosChNegEta ) { this->fTPC0MPosChNegEta = TPC0MPosChNegEta; }
  Double_t getTPC0MPosChNegEta() { return fTPC0MPosChNegEta; }
  void setTPC3MPosChNegEta( Double_t TPC3MPosChNegEta ) { this->fTPC3MPosChNegEta = TPC3MPosChNegEta; }
  Double_t getTPC3MPosChNegEta() { return fTPC3MPosChNegEta; }
  void setTPC4MPosChNegEta( Double_t TPC4MPosChNegEta ) { this->fTPC4MPosChNegEta = TPC4MPosChNegEta; }
  Double_t getTPC4MPosChNegEta() { return fTPC4MPosChNegEta; }
  void setTPC4Re2NegChPosEta( Double_t TPC4Re2NegChPosEta ) { this->fTPC4Re2NegChPosEta = TPC4Re2NegChPosEta; }
  Double_t getTPC4Re2NegChPosEta() { return fTPC4Re2NegChPosEta; }
  void setTPC4Im2NegChPosEta( Double_t TPC4Im2NegChPosEta ) { this->fTPC4Im2NegChPosEta = TPC4Im2NegChPosEta; }
  Double_t getTPC4Im2NegChPosEta() { return fTPC4Im2NegChPosEta; }
  void setTPC2Re3NegChPosEta( Double_t TPC2Re3NegChPosEta ) { this->fTPC2Re3NegChPosEta = TPC2Re3NegChPosEta; }
  Double_t getTPC2Re3NegChPosEta() { return fTPC2Re3NegChPosEta; }
  void setTPC2Im3NegChPosEta( Double_t TPC2Im3NegChPosEta ) { this->fTPC2Im3NegChPosEta = TPC2Im3NegChPosEta; }
  Double_t getTPC2Im3NegChPosEta() { return fTPC2Im3NegChPosEta; }
  void setTPC0MNegChPosEta( Double_t TPC0MNegChPosEta ) { this->fTPC0MNegChPosEta = TPC0MNegChPosEta; }
  Double_t getTPC0MNegChPosEta() { return fTPC0MNegChPosEta; }
  void setTPC3MNegChPosEta( Double_t TPC3MNegChPosEta ) { this->fTPC3MNegChPosEta = TPC3MNegChPosEta; }
  Double_t getTPC3MNegChPosEta() { return fTPC3MNegChPosEta; }
  void setTPC4MNegChPosEta( Double_t TPC4MNegChPosEta ) { this->fTPC4MNegChPosEta = TPC4MNegChPosEta; }
  Double_t getTPC4MNegChPosEta() { return fTPC4MNegChPosEta; }
  void setTPC4Re2NegChNegEta( Double_t TPC4Re2NegChNegEta ) { this->fTPC4Re2NegChNegEta = TPC4Re2NegChNegEta; }
  Double_t getTPC4Re2NegChNegEta() { return fTPC4Re2NegChNegEta; }
  void setTPC4Im2NegChNegEta( Double_t TPC4Im2NegChNegEta ) { this->fTPC4Im2NegChNegEta = TPC4Im2NegChNegEta; }
  Double_t getTPC4Im2NegChNegEta() { return fTPC4Im2NegChNegEta; }
  void setTPC2Re3NegChNegEta( Double_t TPC2Re3NegChNegEta ) { this->fTPC2Re3NegChNegEta = TPC2Re3NegChNegEta; }
  Double_t getTPC2Re3NegChNegEta() { return fTPC2Re3NegChNegEta; }
  void setTPC2Im3NegChNegEta( Double_t TPC2Im3NegChNegEta ) { this->fTPC2Im3NegChNegEta = TPC2Im3NegChNegEta; }
  Double_t getTPC2Im3NegChNegEta() { return fTPC2Im3NegChNegEta; }
  void setTPC0MNegChNegEta( Double_t TPC0MNegChNegEta ) { this->fTPC0MNegChNegEta = TPC0MNegChNegEta; }
  Double_t getTPC0MNegChNegEta() { return fTPC0MNegChNegEta; }
  void setTPC3MNegChNegEta( Double_t TPC3MNegChNegEta ) { this->fTPC3MNegChNegEta = TPC3MNegChNegEta; }
  Double_t getTPC3MNegChNegEta() { return fTPC3MNegChNegEta; }
  void setTPC4MNegChNegEta( Double_t TPC4MNegChNegEta ) { this->fTPC4MNegChNegEta = TPC4MNegChNegEta; }
  Double_t getTPC4MNegChNegEta() { return fTPC4MNegChNegEta; }
  
  
  void setTPCRePosChSubPosEta( Double_t TPCRePosChSubPosEta ) { this->fTPCRePosChSubPosEta = TPCRePosChSubPosEta; }
  Double_t getTPCRePosChSubPosEta() { return fTPCRePosChSubPosEta; }
  void setTPCImPosChSubPosEta( Double_t TPCImPosChSubPosEta ) { this->fTPCImPosChSubPosEta = TPCImPosChSubPosEta; }
  Double_t getTPCImPosChSubPosEta() { return fTPCImPosChSubPosEta; }
  void setTPCMPosChSubPosEta( Double_t TPCMPosChSubPosEta ) { this->fTPCMPosChSubPosEta = TPCMPosChSubPosEta; }
  Double_t getTPCMPosChSubPosEta() { return fTPCMPosChSubPosEta; }
  void setTPCRePosChSubNegEta( Double_t TPCRePosChSubNegEta ) { this->fTPCRePosChSubNegEta = TPCRePosChSubNegEta; }
  Double_t getTPCRePosChSubNegEta() { return fTPCRePosChSubNegEta; }
  void setTPCImPosChSubNegEta( Double_t TPCImPosChSubNegEta ) { this->fTPCImPosChSubNegEta = TPCImPosChSubNegEta; }
  Double_t getTPCImPosChSubNegEta() { return fTPCImPosChSubNegEta; }
  void setTPCMPosChSubNegEta( Double_t TPCMPosChSubNegEta ) { this->fTPCMPosChSubNegEta = TPCMPosChSubNegEta; }
  Double_t getTPCMPosChSubNegEta() { return fTPCMPosChSubNegEta; }
  void setTPCReNegChSubPosEta( Double_t TPCReNegChSubPosEta ) { this->fTPCReNegChSubPosEta = TPCReNegChSubPosEta; }
  Double_t getTPCReNegChSubPosEta() { return fTPCReNegChSubPosEta; }
  void setTPCImNegChSubPosEta( Double_t TPCImNegChSubPosEta ) { this->fTPCImNegChSubPosEta = TPCImNegChSubPosEta; }
  Double_t getTPCImNegChSubPosEta() { return fTPCImNegChSubPosEta; }
  void setTPCMNegChSubPosEta( Double_t TPCMNegChSubPosEta ) { this->fTPCMNegChSubPosEta = TPCMNegChSubPosEta; }
  Double_t getTPCMNegChSubPosEta() { return fTPCMNegChSubPosEta; }
  void setTPCReNegChSubNegEta( Double_t TPCReNegChSubNegEta ) { this->fTPCReNegChSubNegEta = TPCReNegChSubNegEta; }
  Double_t getTPCReNegChSubNegEta() { return fTPCReNegChSubNegEta; }
  void setTPCImNegChSubNegEta( Double_t TPCImNegChSubNegEta ) { this->fTPCImNegChSubNegEta = TPCImNegChSubNegEta; }
  Double_t getTPCImNegChSubNegEta() { return fTPCImNegChSubNegEta; }
  void setTPCMNegChSubNegEta( Double_t TPCMNegChSubNegEta ) { this->fTPCMNegChSubNegEta = TPCMNegChSubNegEta; }
  Double_t getTPCMNegChSubNegEta() { return fTPCMNegChSubNegEta; }

  void setTPC2RePosChSubPosEta( Double_t TPC2RePosChSubPosEta ) { this->fTPC2RePosChSubPosEta = TPC2RePosChSubPosEta; }
  Double_t getTPC2RePosChSubPosEta() { return fTPC2RePosChSubPosEta; }
  void setTPC2ImPosChSubPosEta( Double_t TPC2ImPosChSubPosEta ) { this->fTPC2ImPosChSubPosEta = TPC2ImPosChSubPosEta; }
  Double_t getTPC2ImPosChSubPosEta() { return fTPC2ImPosChSubPosEta; }
  void setTPC2Re2PosChSubPosEta( Double_t TPC2Re2PosChSubPosEta ) { this->fTPC2Re2PosChSubPosEta = TPC2Re2PosChSubPosEta; }
  Double_t getTPC2Re2PosChSubPosEta() { return fTPC2Re2PosChSubPosEta; }
  void setTPC2Im2PosChSubPosEta( Double_t TPC2Im2PosChSubPosEta ) { this->fTPC2Im2PosChSubPosEta = TPC2Im2PosChSubPosEta; }
  Double_t getTPC2Im2PosChSubPosEta() { return fTPC2Im2PosChSubPosEta; }
  void setTPC2MPosChSubPosEta( Double_t TPC2MPosChSubPosEta ) { this->fTPC2MPosChSubPosEta = TPC2MPosChSubPosEta; }
  Double_t getTPC2MPosChSubPosEta() { return fTPC2MPosChSubPosEta; }
  void setTPC2RePosChSubNegEta( Double_t TPC2RePosChSubNegEta ) { this->fTPC2RePosChSubNegEta = TPC2RePosChSubNegEta; }
  Double_t getTPC2RePosChSubNegEta() { return fTPC2RePosChSubNegEta; }
  void setTPC2ImPosChSubNegEta( Double_t TPC2ImPosChSubNegEta ) { this->fTPC2ImPosChSubNegEta = TPC2ImPosChSubNegEta; }
  Double_t getTPC2ImPosChSubNegEta() { return fTPC2ImPosChSubNegEta; }
  void setTPC2Re2PosChSubNegEta( Double_t TPC2Re2PosChSubNegEta ) { this->fTPC2Re2PosChSubNegEta = TPC2Re2PosChSubNegEta; }
  Double_t getTPC2Re2PosChSubNegEta() { return fTPC2Re2PosChSubNegEta; }
  void setTPC2Im2PosChSubNegEta( Double_t TPC2Im2PosChSubNegEta ) { this->fTPC2Im2PosChSubNegEta = TPC2Im2PosChSubNegEta; }
  Double_t getTPC2Im2PosChSubNegEta() { return fTPC2Im2PosChSubNegEta; }
  void setTPC2MPosChSubNegEta( Double_t TPC2MPosChSubNegEta ) { this->fTPC2MPosChSubNegEta = TPC2MPosChSubNegEta; }
  Double_t getTPC2MPosChSubNegEta() { return fTPC2MPosChSubNegEta; }
  void setTPC2ReNegChSubPosEta( Double_t TPC2ReNegChSubPosEta ) { this->fTPC2ReNegChSubPosEta = TPC2ReNegChSubPosEta; }
  Double_t getTPC2ReNegChSubPosEta() { return fTPC2ReNegChSubPosEta; }
  void setTPC2ImNegChSubPosEta( Double_t TPC2ImNegChSubPosEta ) { this->fTPC2ImNegChSubPosEta = TPC2ImNegChSubPosEta; }
  Double_t getTPC2ImNegChSubPosEta() { return fTPC2ImNegChSubPosEta; }
  void setTPC2Re2NegChSubPosEta( Double_t TPC2Re2NegChSubPosEta ) { this->fTPC2Re2NegChSubPosEta = TPC2Re2NegChSubPosEta; }
  Double_t getTPC2Re2NegChSubPosEta() { return fTPC2Re2NegChSubPosEta; }
  void setTPC2Im2NegChSubPosEta( Double_t TPC2Im2NegChSubPosEta ) { this->fTPC2Im2NegChSubPosEta = TPC2Im2NegChSubPosEta; }
  Double_t getTPC2Im2NegChSubPosEta() { return fTPC2Im2NegChSubPosEta; }
  void setTPC2MNegChSubPosEta( Double_t TPC2MNegChSubPosEta ) { this->fTPC2MNegChSubPosEta = TPC2MNegChSubPosEta; }
  Double_t getTPC2MNegChSubPosEta() { return fTPC2MNegChSubPosEta; }
  void setTPC2ReNegChSubNegEta( Double_t TPC2ReNegChSubNegEta ) { this->fTPC2ReNegChSubNegEta = TPC2ReNegChSubNegEta; }
  Double_t getTPC2ReNegChSubNegEta() { return fTPC2ReNegChSubNegEta; }
  void setTPC2ImNegChSubNegEta( Double_t TPC2ImNegChSubNegEta ) { this->fTPC2ImNegChSubNegEta = TPC2ImNegChSubNegEta; }
  Double_t getTPC2ImNegChSubNegEta() { return fTPC2ImNegChSubNegEta; }
  void setTPC2Re2NegChSubNegEta( Double_t TPC2Re2NegChSubNegEta ) { this->fTPC2Re2NegChSubNegEta = TPC2Re2NegChSubNegEta; }
  Double_t getTPC2Re2NegChSubNegEta() { return fTPC2Re2NegChSubNegEta; }
  void setTPC2Im2NegChSubNegEta( Double_t TPC2Im2NegChSubNegEta ) { this->fTPC2Im2NegChSubNegEta = TPC2Im2NegChSubNegEta; }
  Double_t getTPC2Im2NegChSubNegEta() { return fTPC2Im2NegChSubNegEta; }
  void setTPC2MNegChSubNegEta( Double_t TPC2MNegChSubNegEta ) { this->fTPC2MNegChSubNegEta = TPC2MNegChSubNegEta; }
  Double_t getTPC2MNegChSubNegEta() { return fTPC2MNegChSubNegEta; }
  
  void setTPC4Re2PosChSubPosEta( Double_t TPC4Re2PosChSubPosEta ) { this->fTPC4Re2PosChSubPosEta = TPC4Re2PosChSubPosEta; }
  Double_t getTPC4Re2PosChSubPosEta() { return fTPC4Re2PosChSubPosEta; }
  void setTPC4Im2PosChSubPosEta( Double_t TPC4Im2PosChSubPosEta ) { this->fTPC4Im2PosChSubPosEta = TPC4Im2PosChSubPosEta; }
  Double_t getTPC4Im2PosChSubPosEta() { return fTPC4Im2PosChSubPosEta; }
  void setTPC2Re3PosChSubPosEta( Double_t TPC2Re3PosChSubPosEta ) { this->fTPC2Re3PosChSubPosEta = TPC2Re3PosChSubPosEta; }
  Double_t getTPC2Re3PosChSubPosEta() { return fTPC2Re3PosChSubPosEta; }
  void setTPC2Im3PosChSubPosEta( Double_t TPC2Im3PosChSubPosEta ) { this->fTPC2Im3PosChSubPosEta = TPC2Im3PosChSubPosEta; }
  Double_t getTPC2Im3PosChSubPosEta() { return fTPC2Im3PosChSubPosEta; }
  void setTPC0MPosChSubPosEta( Double_t TPC0MPosChSubPosEta ) { this->fTPC0MPosChSubPosEta = TPC0MPosChSubPosEta; }
  Double_t getTPC0MPosChSubPosEta() { return fTPC0MPosChSubPosEta; }
  void setTPC3MPosChSubPosEta( Double_t TPC3MPosChSubPosEta ) { this->fTPC3MPosChSubPosEta = TPC3MPosChSubPosEta; }
  Double_t getTPC3MPosChSubPosEta() { return fTPC3MPosChSubPosEta; }
  void setTPC4MPosChSubPosEta( Double_t TPC4MPosChSubPosEta ) { this->fTPC4MPosChSubPosEta = TPC4MPosChSubPosEta; }
  Double_t getTPC4MPosChSubPosEta() { return fTPC4MPosChSubPosEta; }
  void setTPC4Re2PosChSubNegEta( Double_t TPC4Re2PosChSubNegEta ) { this->fTPC4Re2PosChSubNegEta = TPC4Re2PosChSubNegEta; }
  Double_t getTPC4Re2PosChSubNegEta() { return fTPC4Re2PosChSubNegEta; }
  void setTPC4Im2PosChSubNegEta( Double_t TPC4Im2PosChSubNegEta ) { this->fTPC4Im2PosChSubNegEta = TPC4Im2PosChSubNegEta; }
  Double_t getTPC4Im2PosChSubNegEta() { return fTPC4Im2PosChSubNegEta; }
  void setTPC2Re3PosChSubNegEta( Double_t TPC2Re3PosChSubNegEta ) { this->fTPC2Re3PosChSubNegEta = TPC2Re3PosChSubNegEta; }
  Double_t getTPC2Re3PosChSubNegEta() { return fTPC2Re3PosChSubNegEta; }
  void setTPC2Im3PosChSubNegEta( Double_t TPC2Im3PosChSubNegEta ) { this->fTPC2Im3PosChSubNegEta = TPC2Im3PosChSubNegEta; }
  Double_t getTPC2Im3PosChSubNegEta() { return fTPC2Im3PosChSubNegEta; }
  void setTPC0MPosChSubNegEta( Double_t TPC0MPosChSubNegEta ) { this->fTPC0MPosChSubNegEta = TPC0MPosChSubNegEta; }
  Double_t getTPC0MPosChSubNegEta() { return fTPC0MPosChSubNegEta; }
  void setTPC3MPosChSubNegEta( Double_t TPC3MPosChSubNegEta ) { this->fTPC3MPosChSubNegEta = TPC3MPosChSubNegEta; }
  Double_t getTPC3MPosChSubNegEta() { return fTPC3MPosChSubNegEta; }
  void setTPC4MPosChSubNegEta( Double_t TPC4MPosChSubNegEta ) { this->fTPC4MPosChSubNegEta = TPC4MPosChSubNegEta; }
  Double_t getTPC4MPosChSubNegEta() { return fTPC4MPosChSubNegEta; }
  void setTPC4Re2NegChSubPosEta( Double_t TPC4Re2NegChSubPosEta ) { this->fTPC4Re2NegChSubPosEta = TPC4Re2NegChSubPosEta; }
  Double_t getTPC4Re2NegChSubPosEta() { return fTPC4Re2NegChSubPosEta; }
  void setTPC4Im2NegChSubPosEta( Double_t TPC4Im2NegChSubPosEta ) { this->fTPC4Im2NegChSubPosEta = TPC4Im2NegChSubPosEta; }
  Double_t getTPC4Im2NegChSubPosEta() { return fTPC4Im2NegChSubPosEta; }
  void setTPC2Re3NegChSubPosEta( Double_t TPC2Re3NegChSubPosEta ) { this->fTPC2Re3NegChSubPosEta = TPC2Re3NegChSubPosEta; }
  Double_t getTPC2Re3NegChSubPosEta() { return fTPC2Re3NegChSubPosEta; }
  void setTPC2Im3NegChSubPosEta( Double_t TPC2Im3NegChSubPosEta ) { this->fTPC2Im3NegChSubPosEta = TPC2Im3NegChSubPosEta; }
  Double_t getTPC2Im3NegChSubPosEta() { return fTPC2Im3NegChSubPosEta; }
  void setTPC0MNegChSubPosEta( Double_t TPC0MNegChSubPosEta ) { this->fTPC0MNegChSubPosEta = TPC0MNegChSubPosEta; }
  Double_t getTPC0MNegChSubPosEta() { return fTPC0MNegChSubPosEta; }
  void setTPC3MNegChSubPosEta( Double_t TPC3MNegChSubPosEta ) { this->fTPC3MNegChSubPosEta = TPC3MNegChSubPosEta; }
  Double_t getTPC3MNegChSubPosEta() { return fTPC3MNegChSubPosEta; }
  void setTPC4MNegChSubPosEta( Double_t TPC4MNegChSubPosEta ) { this->fTPC4MNegChSubPosEta = TPC4MNegChSubPosEta; }
  Double_t getTPC4MNegChSubPosEta() { return fTPC4MNegChSubPosEta; }
  void setTPC4Re2NegChSubNegEta( Double_t TPC4Re2NegChSubNegEta ) { this->fTPC4Re2NegChSubNegEta = TPC4Re2NegChSubNegEta; }
  Double_t getTPC4Re2NegChSubNegEta() { return fTPC4Re2NegChSubNegEta; }
  void setTPC4Im2NegChSubNegEta( Double_t TPC4Im2NegChSubNegEta ) { this->fTPC4Im2NegChSubNegEta = TPC4Im2NegChSubNegEta; }
  Double_t getTPC4Im2NegChSubNegEta() { return fTPC4Im2NegChSubNegEta; }
  void setTPC2Re3NegChSubNegEta( Double_t TPC2Re3NegChSubNegEta ) { this->fTPC2Re3NegChSubNegEta = TPC2Re3NegChSubNegEta; }
  Double_t getTPC2Re3NegChSubNegEta() { return fTPC2Re3NegChSubNegEta; }
  void setTPC2Im3NegChSubNegEta( Double_t TPC2Im3NegChSubNegEta ) { this->fTPC2Im3NegChSubNegEta = TPC2Im3NegChSubNegEta; }
  Double_t getTPC2Im3NegChSubNegEta() { return fTPC2Im3NegChSubNegEta; }
  void setTPC0MNegChSubNegEta( Double_t TPC0MNegChSubNegEta ) { this->fTPC0MNegChSubNegEta = TPC0MNegChSubNegEta; }
  Double_t getTPC0MNegChSubNegEta() { return fTPC0MNegChSubNegEta; }
  void setTPC3MNegChSubNegEta( Double_t TPC3MNegChSubNegEta ) { this->fTPC3MNegChSubNegEta = TPC3MNegChSubNegEta; }
  Double_t getTPC3MNegChSubNegEta() { return fTPC3MNegChSubNegEta; }
  void setTPC4MNegChSubNegEta( Double_t TPC4MNegChSubNegEta ) { this->fTPC4MNegChSubNegEta = TPC4MNegChSubNegEta; }
  Double_t getTPC4MNegChSubNegEta() { return fTPC4MNegChSubNegEta; }
  
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
  
  void setTPCPion4Re2PosChPosEta( Double_t TPCPion4Re2PosChPosEta ) { this->fTPCPion4Re2PosChPosEta = TPCPion4Re2PosChPosEta; }
  Double_t getTPCPion4Re2PosChPosEta() { return fTPCPion4Re2PosChPosEta; }
  void setTPCPion4Im2PosChPosEta( Double_t TPCPion4Im2PosChPosEta ) { this->fTPCPion4Im2PosChPosEta = TPCPion4Im2PosChPosEta; }
  Double_t getTPCPion4Im2PosChPosEta() { return fTPCPion4Im2PosChPosEta; }
  void setTPCPion2Re3PosChPosEta( Double_t TPCPion2Re3PosChPosEta ) { this->fTPCPion2Re3PosChPosEta = TPCPion2Re3PosChPosEta; }
  Double_t getTPCPion2Re3PosChPosEta() { return fTPCPion2Re3PosChPosEta; }
  void setTPCPion2Im3PosChPosEta( Double_t TPCPion2Im3PosChPosEta ) { this->fTPCPion2Im3PosChPosEta = TPCPion2Im3PosChPosEta; }
  Double_t getTPCPion2Im3PosChPosEta() { return fTPCPion2Im3PosChPosEta; }
  void setTPCPion0MPosChPosEta( Double_t TPCPion0MPosChPosEta ) { this->fTPCPion0MPosChPosEta = TPCPion0MPosChPosEta; }
  Double_t getTPCPion0MPosChPosEta() { return fTPCPion0MPosChPosEta; }
  void setTPCPion3MPosChPosEta( Double_t TPCPion3MPosChPosEta ) { this->fTPCPion3MPosChPosEta = TPCPion3MPosChPosEta; }
  Double_t getTPCPion3MPosChPosEta() { return fTPCPion3MPosChPosEta; }
  void setTPCPion4MPosChPosEta( Double_t TPCPion4MPosChPosEta ) { this->fTPCPion4MPosChPosEta = TPCPion4MPosChPosEta; }
  Double_t getTPCPion4MPosChPosEta() { return fTPCPion4MPosChPosEta; }
  void setTPCPion4Re2PosChNegEta( Double_t TPCPion4Re2PosChNegEta ) { this->fTPCPion4Re2PosChNegEta = TPCPion4Re2PosChNegEta; }
  Double_t getTPCPion4Re2PosChNegEta() { return fTPCPion4Re2PosChNegEta; }
  void setTPCPion4Im2PosChNegEta( Double_t TPCPion4Im2PosChNegEta ) { this->fTPCPion4Im2PosChNegEta = TPCPion4Im2PosChNegEta; }
  Double_t getTPCPion4Im2PosChNegEta() { return fTPCPion4Im2PosChNegEta; }
  void setTPCPion2Re3PosChNegEta( Double_t TPCPion2Re3PosChNegEta ) { this->fTPCPion2Re3PosChNegEta = TPCPion2Re3PosChNegEta; }
  Double_t getTPCPion2Re3PosChNegEta() { return fTPCPion2Re3PosChNegEta; }
  void setTPCPion2Im3PosChNegEta( Double_t TPCPion2Im3PosChNegEta ) { this->fTPCPion2Im3PosChNegEta = TPCPion2Im3PosChNegEta; }
  Double_t getTPCPion2Im3PosChNegEta() { return fTPCPion2Im3PosChNegEta; }
  void setTPCPion0MPosChNegEta( Double_t TPCPion0MPosChNegEta ) { this->fTPCPion0MPosChNegEta = TPCPion0MPosChNegEta; }
  Double_t getTPCPion0MPosChNegEta() { return fTPCPion0MPosChNegEta; }
  void setTPCPion3MPosChNegEta( Double_t TPCPion3MPosChNegEta ) { this->fTPCPion3MPosChNegEta = TPCPion3MPosChNegEta; }
  Double_t getTPCPion3MPosChNegEta() { return fTPCPion3MPosChNegEta; }
  void setTPCPion4MPosChNegEta( Double_t TPCPion4MPosChNegEta ) { this->fTPCPion4MPosChNegEta = TPCPion4MPosChNegEta; }
  Double_t getTPCPion4MPosChNegEta() { return fTPCPion4MPosChNegEta; }
  void setTPCPion4Re2NegChPosEta( Double_t TPCPion4Re2NegChPosEta ) { this->fTPCPion4Re2NegChPosEta = TPCPion4Re2NegChPosEta; }
  Double_t getTPCPion4Re2NegChPosEta() { return fTPCPion4Re2NegChPosEta; }
  void setTPCPion4Im2NegChPosEta( Double_t TPCPion4Im2NegChPosEta ) { this->fTPCPion4Im2NegChPosEta = TPCPion4Im2NegChPosEta; }
  Double_t getTPCPion4Im2NegChPosEta() { return fTPCPion4Im2NegChPosEta; }
  void setTPCPion2Re3NegChPosEta( Double_t TPCPion2Re3NegChPosEta ) { this->fTPCPion2Re3NegChPosEta = TPCPion2Re3NegChPosEta; }
  Double_t getTPCPion2Re3NegChPosEta() { return fTPCPion2Re3NegChPosEta; }
  void setTPCPion2Im3NegChPosEta( Double_t TPCPion2Im3NegChPosEta ) { this->fTPCPion2Im3NegChPosEta = TPCPion2Im3NegChPosEta; }
  Double_t getTPCPion2Im3NegChPosEta() { return fTPCPion2Im3NegChPosEta; }
  void setTPCPion0MNegChPosEta( Double_t TPCPion0MNegChPosEta ) { this->fTPCPion0MNegChPosEta = TPCPion0MNegChPosEta; }
  Double_t getTPCPion0MNegChPosEta() { return fTPCPion0MNegChPosEta; }
  void setTPCPion3MNegChPosEta( Double_t TPCPion3MNegChPosEta ) { this->fTPCPion3MNegChPosEta = TPCPion3MNegChPosEta; }
  Double_t getTPCPion3MNegChPosEta() { return fTPCPion3MNegChPosEta; }
  void setTPCPion4MNegChPosEta( Double_t TPCPion4MNegChPosEta ) { this->fTPCPion4MNegChPosEta = TPCPion4MNegChPosEta; }
  Double_t getTPCPion4MNegChPosEta() { return fTPCPion4MNegChPosEta; }
  void setTPCPion4Re2NegChNegEta( Double_t TPCPion4Re2NegChNegEta ) { this->fTPCPion4Re2NegChNegEta = TPCPion4Re2NegChNegEta; }
  Double_t getTPCPion4Re2NegChNegEta() { return fTPCPion4Re2NegChNegEta; }
  void setTPCPion4Im2NegChNegEta( Double_t TPCPion4Im2NegChNegEta ) { this->fTPCPion4Im2NegChNegEta = TPCPion4Im2NegChNegEta; }
  Double_t getTPCPion4Im2NegChNegEta() { return fTPCPion4Im2NegChNegEta; }
  void setTPCPion2Re3NegChNegEta( Double_t TPCPion2Re3NegChNegEta ) { this->fTPCPion2Re3NegChNegEta = TPCPion2Re3NegChNegEta; }
  Double_t getTPCPion2Re3NegChNegEta() { return fTPCPion2Re3NegChNegEta; }
  void setTPCPion2Im3NegChNegEta( Double_t TPCPion2Im3NegChNegEta ) { this->fTPCPion2Im3NegChNegEta = TPCPion2Im3NegChNegEta; }
  Double_t getTPCPion2Im3NegChNegEta() { return fTPCPion2Im3NegChNegEta; }
  void setTPCPion0MNegChNegEta( Double_t TPCPion0MNegChNegEta ) { this->fTPCPion0MNegChNegEta = TPCPion0MNegChNegEta; }
  Double_t getTPCPion0MNegChNegEta() { return fTPCPion0MNegChNegEta; }
  void setTPCPion3MNegChNegEta( Double_t TPCPion3MNegChNegEta ) { this->fTPCPion3MNegChNegEta = TPCPion3MNegChNegEta; }
  Double_t getTPCPion3MNegChNegEta() { return fTPCPion3MNegChNegEta; }
  void setTPCPion4MNegChNegEta( Double_t TPCPion4MNegChNegEta ) { this->fTPCPion4MNegChNegEta = TPCPion4MNegChNegEta; }
  Double_t getTPCPion4MNegChNegEta() { return fTPCPion4MNegChNegEta; }
  
  
  void setTPCPionRePosChSubPosEta( Double_t TPCPionRePosChSubPosEta ) { this->fTPCPionRePosChSubPosEta = TPCPionRePosChSubPosEta; }
  Double_t getTPCPionRePosChSubPosEta() { return fTPCPionRePosChSubPosEta; }
  void setTPCPionImPosChSubPosEta( Double_t TPCPionImPosChSubPosEta ) { this->fTPCPionImPosChSubPosEta = TPCPionImPosChSubPosEta; }
  Double_t getTPCPionImPosChSubPosEta() { return fTPCPionImPosChSubPosEta; }
  void setTPCPionMPosChSubPosEta( Double_t TPCPionMPosChSubPosEta ) { this->fTPCPionMPosChSubPosEta = TPCPionMPosChSubPosEta; }
  Double_t getTPCPionMPosChSubPosEta() { return fTPCPionMPosChSubPosEta; }
  void setTPCPionRePosChSubNegEta( Double_t TPCPionRePosChSubNegEta ) { this->fTPCPionRePosChSubNegEta = TPCPionRePosChSubNegEta; }
  Double_t getTPCPionRePosChSubNegEta() { return fTPCPionRePosChSubNegEta; }
  void setTPCPionImPosChSubNegEta( Double_t TPCPionImPosChSubNegEta ) { this->fTPCPionImPosChSubNegEta = TPCPionImPosChSubNegEta; }
  Double_t getTPCPionImPosChSubNegEta() { return fTPCPionImPosChSubNegEta; }
  void setTPCPionMPosChSubNegEta( Double_t TPCPionMPosChSubNegEta ) { this->fTPCPionMPosChSubNegEta = TPCPionMPosChSubNegEta; }
  Double_t getTPCPionMPosChSubNegEta() { return fTPCPionMPosChSubNegEta; }
  void setTPCPionReNegChSubPosEta( Double_t TPCPionReNegChSubPosEta ) { this->fTPCPionReNegChSubPosEta = TPCPionReNegChSubPosEta; }
  Double_t getTPCPionReNegChSubPosEta() { return fTPCPionReNegChSubPosEta; }
  void setTPCPionImNegChSubPosEta( Double_t TPCPionImNegChSubPosEta ) { this->fTPCPionImNegChSubPosEta = TPCPionImNegChSubPosEta; }
  Double_t getTPCPionImNegChSubPosEta() { return fTPCPionImNegChSubPosEta; }
  void setTPCPionMNegChSubPosEta( Double_t TPCPionMNegChSubPosEta ) { this->fTPCPionMNegChSubPosEta = TPCPionMNegChSubPosEta; }
  Double_t getTPCPionMNegChSubPosEta() { return fTPCPionMNegChSubPosEta; }
  void setTPCPionReNegChSubNegEta( Double_t TPCPionReNegChSubNegEta ) { this->fTPCPionReNegChSubNegEta = TPCPionReNegChSubNegEta; }
  Double_t getTPCPionReNegChSubNegEta() { return fTPCPionReNegChSubNegEta; }
  void setTPCPionImNegChSubNegEta( Double_t TPCPionImNegChSubNegEta ) { this->fTPCPionImNegChSubNegEta = TPCPionImNegChSubNegEta; }
  Double_t getTPCPionImNegChSubNegEta() { return fTPCPionImNegChSubNegEta; }
  void setTPCPionMNegChSubNegEta( Double_t TPCPionMNegChSubNegEta ) { this->fTPCPionMNegChSubNegEta = TPCPionMNegChSubNegEta; }
  Double_t getTPCPionMNegChSubNegEta() { return fTPCPionMNegChSubNegEta; }

  void setTPCPion2RePosChSubPosEta( Double_t TPCPion2RePosChSubPosEta ) { this->fTPCPion2RePosChSubPosEta = TPCPion2RePosChSubPosEta; }
  Double_t getTPCPion2RePosChSubPosEta() { return fTPCPion2RePosChSubPosEta; }
  void setTPCPion2ImPosChSubPosEta( Double_t TPCPion2ImPosChSubPosEta ) { this->fTPCPion2ImPosChSubPosEta = TPCPion2ImPosChSubPosEta; }
  Double_t getTPCPion2ImPosChSubPosEta() { return fTPCPion2ImPosChSubPosEta; }
  void setTPCPion2Re2PosChSubPosEta( Double_t TPCPion2Re2PosChSubPosEta ) { this->fTPCPion2Re2PosChSubPosEta = TPCPion2Re2PosChSubPosEta; }
  Double_t getTPCPion2Re2PosChSubPosEta() { return fTPCPion2Re2PosChSubPosEta; }
  void setTPCPion2Im2PosChSubPosEta( Double_t TPCPion2Im2PosChSubPosEta ) { this->fTPCPion2Im2PosChSubPosEta = TPCPion2Im2PosChSubPosEta; }
  Double_t getTPCPion2Im2PosChSubPosEta() { return fTPCPion2Im2PosChSubPosEta; }
  void setTPCPion2MPosChSubPosEta( Double_t TPCPion2MPosChSubPosEta ) { this->fTPCPion2MPosChSubPosEta = TPCPion2MPosChSubPosEta; }
  Double_t getTPCPion2MPosChSubPosEta() { return fTPCPion2MPosChSubPosEta; }
  void setTPCPion2RePosChSubNegEta( Double_t TPCPion2RePosChSubNegEta ) { this->fTPCPion2RePosChSubNegEta = TPCPion2RePosChSubNegEta; }
  Double_t getTPCPion2RePosChSubNegEta() { return fTPCPion2RePosChSubNegEta; }
  void setTPCPion2ImPosChSubNegEta( Double_t TPCPion2ImPosChSubNegEta ) { this->fTPCPion2ImPosChSubNegEta = TPCPion2ImPosChSubNegEta; }
  Double_t getTPCPion2ImPosChSubNegEta() { return fTPCPion2ImPosChSubNegEta; }
  void setTPCPion2Re2PosChSubNegEta( Double_t TPCPion2Re2PosChSubNegEta ) { this->fTPCPion2Re2PosChSubNegEta = TPCPion2Re2PosChSubNegEta; }
  Double_t getTPCPion2Re2PosChSubNegEta() { return fTPCPion2Re2PosChSubNegEta; }
  void setTPCPion2Im2PosChSubNegEta( Double_t TPCPion2Im2PosChSubNegEta ) { this->fTPCPion2Im2PosChSubNegEta = TPCPion2Im2PosChSubNegEta; }
  Double_t getTPCPion2Im2PosChSubNegEta() { return fTPCPion2Im2PosChSubNegEta; }
  void setTPCPion2MPosChSubNegEta( Double_t TPCPion2MPosChSubNegEta ) { this->fTPCPion2MPosChSubNegEta = TPCPion2MPosChSubNegEta; }
  Double_t getTPCPion2MPosChSubNegEta() { return fTPCPion2MPosChSubNegEta; }
  void setTPCPion2ReNegChSubPosEta( Double_t TPCPion2ReNegChSubPosEta ) { this->fTPCPion2ReNegChSubPosEta = TPCPion2ReNegChSubPosEta; }
  Double_t getTPCPion2ReNegChSubPosEta() { return fTPCPion2ReNegChSubPosEta; }
  void setTPCPion2ImNegChSubPosEta( Double_t TPCPion2ImNegChSubPosEta ) { this->fTPCPion2ImNegChSubPosEta = TPCPion2ImNegChSubPosEta; }
  Double_t getTPCPion2ImNegChSubPosEta() { return fTPCPion2ImNegChSubPosEta; }
  void setTPCPion2Re2NegChSubPosEta( Double_t TPCPion2Re2NegChSubPosEta ) { this->fTPCPion2Re2NegChSubPosEta = TPCPion2Re2NegChSubPosEta; }
  Double_t getTPCPion2Re2NegChSubPosEta() { return fTPCPion2Re2NegChSubPosEta; }
  void setTPCPion2Im2NegChSubPosEta( Double_t TPCPion2Im2NegChSubPosEta ) { this->fTPCPion2Im2NegChSubPosEta = TPCPion2Im2NegChSubPosEta; }
  Double_t getTPCPion2Im2NegChSubPosEta() { return fTPCPion2Im2NegChSubPosEta; }
  void setTPCPion2MNegChSubPosEta( Double_t TPCPion2MNegChSubPosEta ) { this->fTPCPion2MNegChSubPosEta = TPCPion2MNegChSubPosEta; }
  Double_t getTPCPion2MNegChSubPosEta() { return fTPCPion2MNegChSubPosEta; }
  void setTPCPion2ReNegChSubNegEta( Double_t TPCPion2ReNegChSubNegEta ) { this->fTPCPion2ReNegChSubNegEta = TPCPion2ReNegChSubNegEta; }
  Double_t getTPCPion2ReNegChSubNegEta() { return fTPCPion2ReNegChSubNegEta; }
  void setTPCPion2ImNegChSubNegEta( Double_t TPCPion2ImNegChSubNegEta ) { this->fTPCPion2ImNegChSubNegEta = TPCPion2ImNegChSubNegEta; }
  Double_t getTPCPion2ImNegChSubNegEta() { return fTPCPion2ImNegChSubNegEta; }
  void setTPCPion2Re2NegChSubNegEta( Double_t TPCPion2Re2NegChSubNegEta ) { this->fTPCPion2Re2NegChSubNegEta = TPCPion2Re2NegChSubNegEta; }
  Double_t getTPCPion2Re2NegChSubNegEta() { return fTPCPion2Re2NegChSubNegEta; }
  void setTPCPion2Im2NegChSubNegEta( Double_t TPCPion2Im2NegChSubNegEta ) { this->fTPCPion2Im2NegChSubNegEta = TPCPion2Im2NegChSubNegEta; }
  Double_t getTPCPion2Im2NegChSubNegEta() { return fTPCPion2Im2NegChSubNegEta; }
  void setTPCPion2MNegChSubNegEta( Double_t TPCPion2MNegChSubNegEta ) { this->fTPCPion2MNegChSubNegEta = TPCPion2MNegChSubNegEta; }
  Double_t getTPCPion2MNegChSubNegEta() { return fTPCPion2MNegChSubNegEta; }
  
  void setTPCPion4Re2PosChSubPosEta( Double_t TPCPion4Re2PosChSubPosEta ) { this->fTPCPion4Re2PosChSubPosEta = TPCPion4Re2PosChSubPosEta; }
  Double_t getTPCPion4Re2PosChSubPosEta() { return fTPCPion4Re2PosChSubPosEta; }
  void setTPCPion4Im2PosChSubPosEta( Double_t TPCPion4Im2PosChSubPosEta ) { this->fTPCPion4Im2PosChSubPosEta = TPCPion4Im2PosChSubPosEta; }
  Double_t getTPCPion4Im2PosChSubPosEta() { return fTPCPion4Im2PosChSubPosEta; }
  void setTPCPion2Re3PosChSubPosEta( Double_t TPCPion2Re3PosChSubPosEta ) { this->fTPCPion2Re3PosChSubPosEta = TPCPion2Re3PosChSubPosEta; }
  Double_t getTPCPion2Re3PosChSubPosEta() { return fTPCPion2Re3PosChSubPosEta; }
  void setTPCPion2Im3PosChSubPosEta( Double_t TPCPion2Im3PosChSubPosEta ) { this->fTPCPion2Im3PosChSubPosEta = TPCPion2Im3PosChSubPosEta; }
  Double_t getTPCPion2Im3PosChSubPosEta() { return fTPCPion2Im3PosChSubPosEta; }
  void setTPCPion0MPosChSubPosEta( Double_t TPCPion0MPosChSubPosEta ) { this->fTPCPion0MPosChSubPosEta = TPCPion0MPosChSubPosEta; }
  Double_t getTPCPion0MPosChSubPosEta() { return fTPCPion0MPosChSubPosEta; }
  void setTPCPion3MPosChSubPosEta( Double_t TPCPion3MPosChSubPosEta ) { this->fTPCPion3MPosChSubPosEta = TPCPion3MPosChSubPosEta; }
  Double_t getTPCPion3MPosChSubPosEta() { return fTPCPion3MPosChSubPosEta; }
  void setTPCPion4MPosChSubPosEta( Double_t TPCPion4MPosChSubPosEta ) { this->fTPCPion4MPosChSubPosEta = TPCPion4MPosChSubPosEta; }
  Double_t getTPCPion4MPosChSubPosEta() { return fTPCPion4MPosChSubPosEta; }
  void setTPCPion4Re2PosChSubNegEta( Double_t TPCPion4Re2PosChSubNegEta ) { this->fTPCPion4Re2PosChSubNegEta = TPCPion4Re2PosChSubNegEta; }
  Double_t getTPCPion4Re2PosChSubNegEta() { return fTPCPion4Re2PosChSubNegEta; }
  void setTPCPion4Im2PosChSubNegEta( Double_t TPCPion4Im2PosChSubNegEta ) { this->fTPCPion4Im2PosChSubNegEta = TPCPion4Im2PosChSubNegEta; }
  Double_t getTPCPion4Im2PosChSubNegEta() { return fTPCPion4Im2PosChSubNegEta; }
  void setTPCPion2Re3PosChSubNegEta( Double_t TPCPion2Re3PosChSubNegEta ) { this->fTPCPion2Re3PosChSubNegEta = TPCPion2Re3PosChSubNegEta; }
  Double_t getTPCPion2Re3PosChSubNegEta() { return fTPCPion2Re3PosChSubNegEta; }
  void setTPCPion2Im3PosChSubNegEta( Double_t TPCPion2Im3PosChSubNegEta ) { this->fTPCPion2Im3PosChSubNegEta = TPCPion2Im3PosChSubNegEta; }
  Double_t getTPCPion2Im3PosChSubNegEta() { return fTPCPion2Im3PosChSubNegEta; }
  void setTPCPion0MPosChSubNegEta( Double_t TPCPion0MPosChSubNegEta ) { this->fTPCPion0MPosChSubNegEta = TPCPion0MPosChSubNegEta; }
  Double_t getTPCPion0MPosChSubNegEta() { return fTPCPion0MPosChSubNegEta; }
  void setTPCPion3MPosChSubNegEta( Double_t TPCPion3MPosChSubNegEta ) { this->fTPCPion3MPosChSubNegEta = TPCPion3MPosChSubNegEta; }
  Double_t getTPCPion3MPosChSubNegEta() { return fTPCPion3MPosChSubNegEta; }
  void setTPCPion4MPosChSubNegEta( Double_t TPCPion4MPosChSubNegEta ) { this->fTPCPion4MPosChSubNegEta = TPCPion4MPosChSubNegEta; }
  Double_t getTPCPion4MPosChSubNegEta() { return fTPCPion4MPosChSubNegEta; }
  void setTPCPion4Re2NegChSubPosEta( Double_t TPCPion4Re2NegChSubPosEta ) { this->fTPCPion4Re2NegChSubPosEta = TPCPion4Re2NegChSubPosEta; }
  Double_t getTPCPion4Re2NegChSubPosEta() { return fTPCPion4Re2NegChSubPosEta; }
  void setTPCPion4Im2NegChSubPosEta( Double_t TPCPion4Im2NegChSubPosEta ) { this->fTPCPion4Im2NegChSubPosEta = TPCPion4Im2NegChSubPosEta; }
  Double_t getTPCPion4Im2NegChSubPosEta() { return fTPCPion4Im2NegChSubPosEta; }
  void setTPCPion2Re3NegChSubPosEta( Double_t TPCPion2Re3NegChSubPosEta ) { this->fTPCPion2Re3NegChSubPosEta = TPCPion2Re3NegChSubPosEta; }
  Double_t getTPCPion2Re3NegChSubPosEta() { return fTPCPion2Re3NegChSubPosEta; }
  void setTPCPion2Im3NegChSubPosEta( Double_t TPCPion2Im3NegChSubPosEta ) { this->fTPCPion2Im3NegChSubPosEta = TPCPion2Im3NegChSubPosEta; }
  Double_t getTPCPion2Im3NegChSubPosEta() { return fTPCPion2Im3NegChSubPosEta; }
  void setTPCPion0MNegChSubPosEta( Double_t TPCPion0MNegChSubPosEta ) { this->fTPCPion0MNegChSubPosEta = TPCPion0MNegChSubPosEta; }
  Double_t getTPCPion0MNegChSubPosEta() { return fTPCPion0MNegChSubPosEta; }
  void setTPCPion3MNegChSubPosEta( Double_t TPCPion3MNegChSubPosEta ) { this->fTPCPion3MNegChSubPosEta = TPCPion3MNegChSubPosEta; }
  Double_t getTPCPion3MNegChSubPosEta() { return fTPCPion3MNegChSubPosEta; }
  void setTPCPion4MNegChSubPosEta( Double_t TPCPion4MNegChSubPosEta ) { this->fTPCPion4MNegChSubPosEta = TPCPion4MNegChSubPosEta; }
  Double_t getTPCPion4MNegChSubPosEta() { return fTPCPion4MNegChSubPosEta; }
  void setTPCPion4Re2NegChSubNegEta( Double_t TPCPion4Re2NegChSubNegEta ) { this->fTPCPion4Re2NegChSubNegEta = TPCPion4Re2NegChSubNegEta; }
  Double_t getTPCPion4Re2NegChSubNegEta() { return fTPCPion4Re2NegChSubNegEta; }
  void setTPCPion4Im2NegChSubNegEta( Double_t TPCPion4Im2NegChSubNegEta ) { this->fTPCPion4Im2NegChSubNegEta = TPCPion4Im2NegChSubNegEta; }
  Double_t getTPCPion4Im2NegChSubNegEta() { return fTPCPion4Im2NegChSubNegEta; }
  void setTPCPion2Re3NegChSubNegEta( Double_t TPCPion2Re3NegChSubNegEta ) { this->fTPCPion2Re3NegChSubNegEta = TPCPion2Re3NegChSubNegEta; }
  Double_t getTPCPion2Re3NegChSubNegEta() { return fTPCPion2Re3NegChSubNegEta; }
  void setTPCPion2Im3NegChSubNegEta( Double_t TPCPion2Im3NegChSubNegEta ) { this->fTPCPion2Im3NegChSubNegEta = TPCPion2Im3NegChSubNegEta; }
  Double_t getTPCPion2Im3NegChSubNegEta() { return fTPCPion2Im3NegChSubNegEta; }
  void setTPCPion0MNegChSubNegEta( Double_t TPCPion0MNegChSubNegEta ) { this->fTPCPion0MNegChSubNegEta = TPCPion0MNegChSubNegEta; }
  Double_t getTPCPion0MNegChSubNegEta() { return fTPCPion0MNegChSubNegEta; }
  void setTPCPion3MNegChSubNegEta( Double_t TPCPion3MNegChSubNegEta ) { this->fTPCPion3MNegChSubNegEta = TPCPion3MNegChSubNegEta; }
  Double_t getTPCPion3MNegChSubNegEta() { return fTPCPion3MNegChSubNegEta; }
  void setTPCPion4MNegChSubNegEta( Double_t TPCPion4MNegChSubNegEta ) { this->fTPCPion4MNegChSubNegEta = TPCPion4MNegChSubNegEta; }
  Double_t getTPCPion4MNegChSubNegEta() { return fTPCPion4MNegChSubNegEta; }
  
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
  
  void setTPCKaon4Re2PosChPosEta( Double_t TPCKaon4Re2PosChPosEta ) { this->fTPCKaon4Re2PosChPosEta = TPCKaon4Re2PosChPosEta; }
  Double_t getTPCKaon4Re2PosChPosEta() { return fTPCKaon4Re2PosChPosEta; }
  void setTPCKaon4Im2PosChPosEta( Double_t TPCKaon4Im2PosChPosEta ) { this->fTPCKaon4Im2PosChPosEta = TPCKaon4Im2PosChPosEta; }
  Double_t getTPCKaon4Im2PosChPosEta() { return fTPCKaon4Im2PosChPosEta; }
  void setTPCKaon2Re3PosChPosEta( Double_t TPCKaon2Re3PosChPosEta ) { this->fTPCKaon2Re3PosChPosEta = TPCKaon2Re3PosChPosEta; }
  Double_t getTPCKaon2Re3PosChPosEta() { return fTPCKaon2Re3PosChPosEta; }
  void setTPCKaon2Im3PosChPosEta( Double_t TPCKaon2Im3PosChPosEta ) { this->fTPCKaon2Im3PosChPosEta = TPCKaon2Im3PosChPosEta; }
  Double_t getTPCKaon2Im3PosChPosEta() { return fTPCKaon2Im3PosChPosEta; }
  void setTPCKaon0MPosChPosEta( Double_t TPCKaon0MPosChPosEta ) { this->fTPCKaon0MPosChPosEta = TPCKaon0MPosChPosEta; }
  Double_t getTPCKaon0MPosChPosEta() { return fTPCKaon0MPosChPosEta; }
  void setTPCKaon3MPosChPosEta( Double_t TPCKaon3MPosChPosEta ) { this->fTPCKaon3MPosChPosEta = TPCKaon3MPosChPosEta; }
  Double_t getTPCKaon3MPosChPosEta() { return fTPCKaon3MPosChPosEta; }
  void setTPCKaon4MPosChPosEta( Double_t TPCKaon4MPosChPosEta ) { this->fTPCKaon4MPosChPosEta = TPCKaon4MPosChPosEta; }
  Double_t getTPCKaon4MPosChPosEta() { return fTPCKaon4MPosChPosEta; }
  void setTPCKaon4Re2PosChNegEta( Double_t TPCKaon4Re2PosChNegEta ) { this->fTPCKaon4Re2PosChNegEta = TPCKaon4Re2PosChNegEta; }
  Double_t getTPCKaon4Re2PosChNegEta() { return fTPCKaon4Re2PosChNegEta; }
  void setTPCKaon4Im2PosChNegEta( Double_t TPCKaon4Im2PosChNegEta ) { this->fTPCKaon4Im2PosChNegEta = TPCKaon4Im2PosChNegEta; }
  Double_t getTPCKaon4Im2PosChNegEta() { return fTPCKaon4Im2PosChNegEta; }
  void setTPCKaon2Re3PosChNegEta( Double_t TPCKaon2Re3PosChNegEta ) { this->fTPCKaon2Re3PosChNegEta = TPCKaon2Re3PosChNegEta; }
  Double_t getTPCKaon2Re3PosChNegEta() { return fTPCKaon2Re3PosChNegEta; }
  void setTPCKaon2Im3PosChNegEta( Double_t TPCKaon2Im3PosChNegEta ) { this->fTPCKaon2Im3PosChNegEta = TPCKaon2Im3PosChNegEta; }
  Double_t getTPCKaon2Im3PosChNegEta() { return fTPCKaon2Im3PosChNegEta; }
  void setTPCKaon0MPosChNegEta( Double_t TPCKaon0MPosChNegEta ) { this->fTPCKaon0MPosChNegEta = TPCKaon0MPosChNegEta; }
  Double_t getTPCKaon0MPosChNegEta() { return fTPCKaon0MPosChNegEta; }
  void setTPCKaon3MPosChNegEta( Double_t TPCKaon3MPosChNegEta ) { this->fTPCKaon3MPosChNegEta = TPCKaon3MPosChNegEta; }
  Double_t getTPCKaon3MPosChNegEta() { return fTPCKaon3MPosChNegEta; }
  void setTPCKaon4MPosChNegEta( Double_t TPCKaon4MPosChNegEta ) { this->fTPCKaon4MPosChNegEta = TPCKaon4MPosChNegEta; }
  Double_t getTPCKaon4MPosChNegEta() { return fTPCKaon4MPosChNegEta; }
  void setTPCKaon4Re2NegChPosEta( Double_t TPCKaon4Re2NegChPosEta ) { this->fTPCKaon4Re2NegChPosEta = TPCKaon4Re2NegChPosEta; }
  Double_t getTPCKaon4Re2NegChPosEta() { return fTPCKaon4Re2NegChPosEta; }
  void setTPCKaon4Im2NegChPosEta( Double_t TPCKaon4Im2NegChPosEta ) { this->fTPCKaon4Im2NegChPosEta = TPCKaon4Im2NegChPosEta; }
  Double_t getTPCKaon4Im2NegChPosEta() { return fTPCKaon4Im2NegChPosEta; }
  void setTPCKaon2Re3NegChPosEta( Double_t TPCKaon2Re3NegChPosEta ) { this->fTPCKaon2Re3NegChPosEta = TPCKaon2Re3NegChPosEta; }
  Double_t getTPCKaon2Re3NegChPosEta() { return fTPCKaon2Re3NegChPosEta; }
  void setTPCKaon2Im3NegChPosEta( Double_t TPCKaon2Im3NegChPosEta ) { this->fTPCKaon2Im3NegChPosEta = TPCKaon2Im3NegChPosEta; }
  Double_t getTPCKaon2Im3NegChPosEta() { return fTPCKaon2Im3NegChPosEta; }
  void setTPCKaon0MNegChPosEta( Double_t TPCKaon0MNegChPosEta ) { this->fTPCKaon0MNegChPosEta = TPCKaon0MNegChPosEta; }
  Double_t getTPCKaon0MNegChPosEta() { return fTPCKaon0MNegChPosEta; }
  void setTPCKaon3MNegChPosEta( Double_t TPCKaon3MNegChPosEta ) { this->fTPCKaon3MNegChPosEta = TPCKaon3MNegChPosEta; }
  Double_t getTPCKaon3MNegChPosEta() { return fTPCKaon3MNegChPosEta; }
  void setTPCKaon4MNegChPosEta( Double_t TPCKaon4MNegChPosEta ) { this->fTPCKaon4MNegChPosEta = TPCKaon4MNegChPosEta; }
  Double_t getTPCKaon4MNegChPosEta() { return fTPCKaon4MNegChPosEta; }
  void setTPCKaon4Re2NegChNegEta( Double_t TPCKaon4Re2NegChNegEta ) { this->fTPCKaon4Re2NegChNegEta = TPCKaon4Re2NegChNegEta; }
  Double_t getTPCKaon4Re2NegChNegEta() { return fTPCKaon4Re2NegChNegEta; }
  void setTPCKaon4Im2NegChNegEta( Double_t TPCKaon4Im2NegChNegEta ) { this->fTPCKaon4Im2NegChNegEta = TPCKaon4Im2NegChNegEta; }
  Double_t getTPCKaon4Im2NegChNegEta() { return fTPCKaon4Im2NegChNegEta; }
  void setTPCKaon2Re3NegChNegEta( Double_t TPCKaon2Re3NegChNegEta ) { this->fTPCKaon2Re3NegChNegEta = TPCKaon2Re3NegChNegEta; }
  Double_t getTPCKaon2Re3NegChNegEta() { return fTPCKaon2Re3NegChNegEta; }
  void setTPCKaon2Im3NegChNegEta( Double_t TPCKaon2Im3NegChNegEta ) { this->fTPCKaon2Im3NegChNegEta = TPCKaon2Im3NegChNegEta; }
  Double_t getTPCKaon2Im3NegChNegEta() { return fTPCKaon2Im3NegChNegEta; }
  void setTPCKaon0MNegChNegEta( Double_t TPCKaon0MNegChNegEta ) { this->fTPCKaon0MNegChNegEta = TPCKaon0MNegChNegEta; }
  Double_t getTPCKaon0MNegChNegEta() { return fTPCKaon0MNegChNegEta; }
  void setTPCKaon3MNegChNegEta( Double_t TPCKaon3MNegChNegEta ) { this->fTPCKaon3MNegChNegEta = TPCKaon3MNegChNegEta; }
  Double_t getTPCKaon3MNegChNegEta() { return fTPCKaon3MNegChNegEta; }
  void setTPCKaon4MNegChNegEta( Double_t TPCKaon4MNegChNegEta ) { this->fTPCKaon4MNegChNegEta = TPCKaon4MNegChNegEta; }
  Double_t getTPCKaon4MNegChNegEta() { return fTPCKaon4MNegChNegEta; }
  
  
  void setTPCKaonRePosChSubPosEta( Double_t TPCKaonRePosChSubPosEta ) { this->fTPCKaonRePosChSubPosEta = TPCKaonRePosChSubPosEta; }
  Double_t getTPCKaonRePosChSubPosEta() { return fTPCKaonRePosChSubPosEta; }
  void setTPCKaonImPosChSubPosEta( Double_t TPCKaonImPosChSubPosEta ) { this->fTPCKaonImPosChSubPosEta = TPCKaonImPosChSubPosEta; }
  Double_t getTPCKaonImPosChSubPosEta() { return fTPCKaonImPosChSubPosEta; }
  void setTPCKaonMPosChSubPosEta( Double_t TPCKaonMPosChSubPosEta ) { this->fTPCKaonMPosChSubPosEta = TPCKaonMPosChSubPosEta; }
  Double_t getTPCKaonMPosChSubPosEta() { return fTPCKaonMPosChSubPosEta; }
  void setTPCKaonRePosChSubNegEta( Double_t TPCKaonRePosChSubNegEta ) { this->fTPCKaonRePosChSubNegEta = TPCKaonRePosChSubNegEta; }
  Double_t getTPCKaonRePosChSubNegEta() { return fTPCKaonRePosChSubNegEta; }
  void setTPCKaonImPosChSubNegEta( Double_t TPCKaonImPosChSubNegEta ) { this->fTPCKaonImPosChSubNegEta = TPCKaonImPosChSubNegEta; }
  Double_t getTPCKaonImPosChSubNegEta() { return fTPCKaonImPosChSubNegEta; }
  void setTPCKaonMPosChSubNegEta( Double_t TPCKaonMPosChSubNegEta ) { this->fTPCKaonMPosChSubNegEta = TPCKaonMPosChSubNegEta; }
  Double_t getTPCKaonMPosChSubNegEta() { return fTPCKaonMPosChSubNegEta; }
  void setTPCKaonReNegChSubPosEta( Double_t TPCKaonReNegChSubPosEta ) { this->fTPCKaonReNegChSubPosEta = TPCKaonReNegChSubPosEta; }
  Double_t getTPCKaonReNegChSubPosEta() { return fTPCKaonReNegChSubPosEta; }
  void setTPCKaonImNegChSubPosEta( Double_t TPCKaonImNegChSubPosEta ) { this->fTPCKaonImNegChSubPosEta = TPCKaonImNegChSubPosEta; }
  Double_t getTPCKaonImNegChSubPosEta() { return fTPCKaonImNegChSubPosEta; }
  void setTPCKaonMNegChSubPosEta( Double_t TPCKaonMNegChSubPosEta ) { this->fTPCKaonMNegChSubPosEta = TPCKaonMNegChSubPosEta; }
  Double_t getTPCKaonMNegChSubPosEta() { return fTPCKaonMNegChSubPosEta; }
  void setTPCKaonReNegChSubNegEta( Double_t TPCKaonReNegChSubNegEta ) { this->fTPCKaonReNegChSubNegEta = TPCKaonReNegChSubNegEta; }
  Double_t getTPCKaonReNegChSubNegEta() { return fTPCKaonReNegChSubNegEta; }
  void setTPCKaonImNegChSubNegEta( Double_t TPCKaonImNegChSubNegEta ) { this->fTPCKaonImNegChSubNegEta = TPCKaonImNegChSubNegEta; }
  Double_t getTPCKaonImNegChSubNegEta() { return fTPCKaonImNegChSubNegEta; }
  void setTPCKaonMNegChSubNegEta( Double_t TPCKaonMNegChSubNegEta ) { this->fTPCKaonMNegChSubNegEta = TPCKaonMNegChSubNegEta; }
  Double_t getTPCKaonMNegChSubNegEta() { return fTPCKaonMNegChSubNegEta; }

  void setTPCKaon2RePosChSubPosEta( Double_t TPCKaon2RePosChSubPosEta ) { this->fTPCKaon2RePosChSubPosEta = TPCKaon2RePosChSubPosEta; }
  Double_t getTPCKaon2RePosChSubPosEta() { return fTPCKaon2RePosChSubPosEta; }
  void setTPCKaon2ImPosChSubPosEta( Double_t TPCKaon2ImPosChSubPosEta ) { this->fTPCKaon2ImPosChSubPosEta = TPCKaon2ImPosChSubPosEta; }
  Double_t getTPCKaon2ImPosChSubPosEta() { return fTPCKaon2ImPosChSubPosEta; }
  void setTPCKaon2Re2PosChSubPosEta( Double_t TPCKaon2Re2PosChSubPosEta ) { this->fTPCKaon2Re2PosChSubPosEta = TPCKaon2Re2PosChSubPosEta; }
  Double_t getTPCKaon2Re2PosChSubPosEta() { return fTPCKaon2Re2PosChSubPosEta; }
  void setTPCKaon2Im2PosChSubPosEta( Double_t TPCKaon2Im2PosChSubPosEta ) { this->fTPCKaon2Im2PosChSubPosEta = TPCKaon2Im2PosChSubPosEta; }
  Double_t getTPCKaon2Im2PosChSubPosEta() { return fTPCKaon2Im2PosChSubPosEta; }
  void setTPCKaon2MPosChSubPosEta( Double_t TPCKaon2MPosChSubPosEta ) { this->fTPCKaon2MPosChSubPosEta = TPCKaon2MPosChSubPosEta; }
  Double_t getTPCKaon2MPosChSubPosEta() { return fTPCKaon2MPosChSubPosEta; }
  void setTPCKaon2RePosChSubNegEta( Double_t TPCKaon2RePosChSubNegEta ) { this->fTPCKaon2RePosChSubNegEta = TPCKaon2RePosChSubNegEta; }
  Double_t getTPCKaon2RePosChSubNegEta() { return fTPCKaon2RePosChSubNegEta; }
  void setTPCKaon2ImPosChSubNegEta( Double_t TPCKaon2ImPosChSubNegEta ) { this->fTPCKaon2ImPosChSubNegEta = TPCKaon2ImPosChSubNegEta; }
  Double_t getTPCKaon2ImPosChSubNegEta() { return fTPCKaon2ImPosChSubNegEta; }
  void setTPCKaon2Re2PosChSubNegEta( Double_t TPCKaon2Re2PosChSubNegEta ) { this->fTPCKaon2Re2PosChSubNegEta = TPCKaon2Re2PosChSubNegEta; }
  Double_t getTPCKaon2Re2PosChSubNegEta() { return fTPCKaon2Re2PosChSubNegEta; }
  void setTPCKaon2Im2PosChSubNegEta( Double_t TPCKaon2Im2PosChSubNegEta ) { this->fTPCKaon2Im2PosChSubNegEta = TPCKaon2Im2PosChSubNegEta; }
  Double_t getTPCKaon2Im2PosChSubNegEta() { return fTPCKaon2Im2PosChSubNegEta; }
  void setTPCKaon2MPosChSubNegEta( Double_t TPCKaon2MPosChSubNegEta ) { this->fTPCKaon2MPosChSubNegEta = TPCKaon2MPosChSubNegEta; }
  Double_t getTPCKaon2MPosChSubNegEta() { return fTPCKaon2MPosChSubNegEta; }
  void setTPCKaon2ReNegChSubPosEta( Double_t TPCKaon2ReNegChSubPosEta ) { this->fTPCKaon2ReNegChSubPosEta = TPCKaon2ReNegChSubPosEta; }
  Double_t getTPCKaon2ReNegChSubPosEta() { return fTPCKaon2ReNegChSubPosEta; }
  void setTPCKaon2ImNegChSubPosEta( Double_t TPCKaon2ImNegChSubPosEta ) { this->fTPCKaon2ImNegChSubPosEta = TPCKaon2ImNegChSubPosEta; }
  Double_t getTPCKaon2ImNegChSubPosEta() { return fTPCKaon2ImNegChSubPosEta; }
  void setTPCKaon2Re2NegChSubPosEta( Double_t TPCKaon2Re2NegChSubPosEta ) { this->fTPCKaon2Re2NegChSubPosEta = TPCKaon2Re2NegChSubPosEta; }
  Double_t getTPCKaon2Re2NegChSubPosEta() { return fTPCKaon2Re2NegChSubPosEta; }
  void setTPCKaon2Im2NegChSubPosEta( Double_t TPCKaon2Im2NegChSubPosEta ) { this->fTPCKaon2Im2NegChSubPosEta = TPCKaon2Im2NegChSubPosEta; }
  Double_t getTPCKaon2Im2NegChSubPosEta() { return fTPCKaon2Im2NegChSubPosEta; }
  void setTPCKaon2MNegChSubPosEta( Double_t TPCKaon2MNegChSubPosEta ) { this->fTPCKaon2MNegChSubPosEta = TPCKaon2MNegChSubPosEta; }
  Double_t getTPCKaon2MNegChSubPosEta() { return fTPCKaon2MNegChSubPosEta; }
  void setTPCKaon2ReNegChSubNegEta( Double_t TPCKaon2ReNegChSubNegEta ) { this->fTPCKaon2ReNegChSubNegEta = TPCKaon2ReNegChSubNegEta; }
  Double_t getTPCKaon2ReNegChSubNegEta() { return fTPCKaon2ReNegChSubNegEta; }
  void setTPCKaon2ImNegChSubNegEta( Double_t TPCKaon2ImNegChSubNegEta ) { this->fTPCKaon2ImNegChSubNegEta = TPCKaon2ImNegChSubNegEta; }
  Double_t getTPCKaon2ImNegChSubNegEta() { return fTPCKaon2ImNegChSubNegEta; }
  void setTPCKaon2Re2NegChSubNegEta( Double_t TPCKaon2Re2NegChSubNegEta ) { this->fTPCKaon2Re2NegChSubNegEta = TPCKaon2Re2NegChSubNegEta; }
  Double_t getTPCKaon2Re2NegChSubNegEta() { return fTPCKaon2Re2NegChSubNegEta; }
  void setTPCKaon2Im2NegChSubNegEta( Double_t TPCKaon2Im2NegChSubNegEta ) { this->fTPCKaon2Im2NegChSubNegEta = TPCKaon2Im2NegChSubNegEta; }
  Double_t getTPCKaon2Im2NegChSubNegEta() { return fTPCKaon2Im2NegChSubNegEta; }
  void setTPCKaon2MNegChSubNegEta( Double_t TPCKaon2MNegChSubNegEta ) { this->fTPCKaon2MNegChSubNegEta = TPCKaon2MNegChSubNegEta; }
  Double_t getTPCKaon2MNegChSubNegEta() { return fTPCKaon2MNegChSubNegEta; }
  
  void setTPCKaon4Re2PosChSubPosEta( Double_t TPCKaon4Re2PosChSubPosEta ) { this->fTPCKaon4Re2PosChSubPosEta = TPCKaon4Re2PosChSubPosEta; }
  Double_t getTPCKaon4Re2PosChSubPosEta() { return fTPCKaon4Re2PosChSubPosEta; }
  void setTPCKaon4Im2PosChSubPosEta( Double_t TPCKaon4Im2PosChSubPosEta ) { this->fTPCKaon4Im2PosChSubPosEta = TPCKaon4Im2PosChSubPosEta; }
  Double_t getTPCKaon4Im2PosChSubPosEta() { return fTPCKaon4Im2PosChSubPosEta; }
  void setTPCKaon2Re3PosChSubPosEta( Double_t TPCKaon2Re3PosChSubPosEta ) { this->fTPCKaon2Re3PosChSubPosEta = TPCKaon2Re3PosChSubPosEta; }
  Double_t getTPCKaon2Re3PosChSubPosEta() { return fTPCKaon2Re3PosChSubPosEta; }
  void setTPCKaon2Im3PosChSubPosEta( Double_t TPCKaon2Im3PosChSubPosEta ) { this->fTPCKaon2Im3PosChSubPosEta = TPCKaon2Im3PosChSubPosEta; }
  Double_t getTPCKaon2Im3PosChSubPosEta() { return fTPCKaon2Im3PosChSubPosEta; }
  void setTPCKaon0MPosChSubPosEta( Double_t TPCKaon0MPosChSubPosEta ) { this->fTPCKaon0MPosChSubPosEta = TPCKaon0MPosChSubPosEta; }
  Double_t getTPCKaon0MPosChSubPosEta() { return fTPCKaon0MPosChSubPosEta; }
  void setTPCKaon3MPosChSubPosEta( Double_t TPCKaon3MPosChSubPosEta ) { this->fTPCKaon3MPosChSubPosEta = TPCKaon3MPosChSubPosEta; }
  Double_t getTPCKaon3MPosChSubPosEta() { return fTPCKaon3MPosChSubPosEta; }
  void setTPCKaon4MPosChSubPosEta( Double_t TPCKaon4MPosChSubPosEta ) { this->fTPCKaon4MPosChSubPosEta = TPCKaon4MPosChSubPosEta; }
  Double_t getTPCKaon4MPosChSubPosEta() { return fTPCKaon4MPosChSubPosEta; }
  void setTPCKaon4Re2PosChSubNegEta( Double_t TPCKaon4Re2PosChSubNegEta ) { this->fTPCKaon4Re2PosChSubNegEta = TPCKaon4Re2PosChSubNegEta; }
  Double_t getTPCKaon4Re2PosChSubNegEta() { return fTPCKaon4Re2PosChSubNegEta; }
  void setTPCKaon4Im2PosChSubNegEta( Double_t TPCKaon4Im2PosChSubNegEta ) { this->fTPCKaon4Im2PosChSubNegEta = TPCKaon4Im2PosChSubNegEta; }
  Double_t getTPCKaon4Im2PosChSubNegEta() { return fTPCKaon4Im2PosChSubNegEta; }
  void setTPCKaon2Re3PosChSubNegEta( Double_t TPCKaon2Re3PosChSubNegEta ) { this->fTPCKaon2Re3PosChSubNegEta = TPCKaon2Re3PosChSubNegEta; }
  Double_t getTPCKaon2Re3PosChSubNegEta() { return fTPCKaon2Re3PosChSubNegEta; }
  void setTPCKaon2Im3PosChSubNegEta( Double_t TPCKaon2Im3PosChSubNegEta ) { this->fTPCKaon2Im3PosChSubNegEta = TPCKaon2Im3PosChSubNegEta; }
  Double_t getTPCKaon2Im3PosChSubNegEta() { return fTPCKaon2Im3PosChSubNegEta; }
  void setTPCKaon0MPosChSubNegEta( Double_t TPCKaon0MPosChSubNegEta ) { this->fTPCKaon0MPosChSubNegEta = TPCKaon0MPosChSubNegEta; }
  Double_t getTPCKaon0MPosChSubNegEta() { return fTPCKaon0MPosChSubNegEta; }
  void setTPCKaon3MPosChSubNegEta( Double_t TPCKaon3MPosChSubNegEta ) { this->fTPCKaon3MPosChSubNegEta = TPCKaon3MPosChSubNegEta; }
  Double_t getTPCKaon3MPosChSubNegEta() { return fTPCKaon3MPosChSubNegEta; }
  void setTPCKaon4MPosChSubNegEta( Double_t TPCKaon4MPosChSubNegEta ) { this->fTPCKaon4MPosChSubNegEta = TPCKaon4MPosChSubNegEta; }
  Double_t getTPCKaon4MPosChSubNegEta() { return fTPCKaon4MPosChSubNegEta; }
  void setTPCKaon4Re2NegChSubPosEta( Double_t TPCKaon4Re2NegChSubPosEta ) { this->fTPCKaon4Re2NegChSubPosEta = TPCKaon4Re2NegChSubPosEta; }
  Double_t getTPCKaon4Re2NegChSubPosEta() { return fTPCKaon4Re2NegChSubPosEta; }
  void setTPCKaon4Im2NegChSubPosEta( Double_t TPCKaon4Im2NegChSubPosEta ) { this->fTPCKaon4Im2NegChSubPosEta = TPCKaon4Im2NegChSubPosEta; }
  Double_t getTPCKaon4Im2NegChSubPosEta() { return fTPCKaon4Im2NegChSubPosEta; }
  void setTPCKaon2Re3NegChSubPosEta( Double_t TPCKaon2Re3NegChSubPosEta ) { this->fTPCKaon2Re3NegChSubPosEta = TPCKaon2Re3NegChSubPosEta; }
  Double_t getTPCKaon2Re3NegChSubPosEta() { return fTPCKaon2Re3NegChSubPosEta; }
  void setTPCKaon2Im3NegChSubPosEta( Double_t TPCKaon2Im3NegChSubPosEta ) { this->fTPCKaon2Im3NegChSubPosEta = TPCKaon2Im3NegChSubPosEta; }
  Double_t getTPCKaon2Im3NegChSubPosEta() { return fTPCKaon2Im3NegChSubPosEta; }
  void setTPCKaon0MNegChSubPosEta( Double_t TPCKaon0MNegChSubPosEta ) { this->fTPCKaon0MNegChSubPosEta = TPCKaon0MNegChSubPosEta; }
  Double_t getTPCKaon0MNegChSubPosEta() { return fTPCKaon0MNegChSubPosEta; }
  void setTPCKaon3MNegChSubPosEta( Double_t TPCKaon3MNegChSubPosEta ) { this->fTPCKaon3MNegChSubPosEta = TPCKaon3MNegChSubPosEta; }
  Double_t getTPCKaon3MNegChSubPosEta() { return fTPCKaon3MNegChSubPosEta; }
  void setTPCKaon4MNegChSubPosEta( Double_t TPCKaon4MNegChSubPosEta ) { this->fTPCKaon4MNegChSubPosEta = TPCKaon4MNegChSubPosEta; }
  Double_t getTPCKaon4MNegChSubPosEta() { return fTPCKaon4MNegChSubPosEta; }
  void setTPCKaon4Re2NegChSubNegEta( Double_t TPCKaon4Re2NegChSubNegEta ) { this->fTPCKaon4Re2NegChSubNegEta = TPCKaon4Re2NegChSubNegEta; }
  Double_t getTPCKaon4Re2NegChSubNegEta() { return fTPCKaon4Re2NegChSubNegEta; }
  void setTPCKaon4Im2NegChSubNegEta( Double_t TPCKaon4Im2NegChSubNegEta ) { this->fTPCKaon4Im2NegChSubNegEta = TPCKaon4Im2NegChSubNegEta; }
  Double_t getTPCKaon4Im2NegChSubNegEta() { return fTPCKaon4Im2NegChSubNegEta; }
  void setTPCKaon2Re3NegChSubNegEta( Double_t TPCKaon2Re3NegChSubNegEta ) { this->fTPCKaon2Re3NegChSubNegEta = TPCKaon2Re3NegChSubNegEta; }
  Double_t getTPCKaon2Re3NegChSubNegEta() { return fTPCKaon2Re3NegChSubNegEta; }
  void setTPCKaon2Im3NegChSubNegEta( Double_t TPCKaon2Im3NegChSubNegEta ) { this->fTPCKaon2Im3NegChSubNegEta = TPCKaon2Im3NegChSubNegEta; }
  Double_t getTPCKaon2Im3NegChSubNegEta() { return fTPCKaon2Im3NegChSubNegEta; }
  void setTPCKaon0MNegChSubNegEta( Double_t TPCKaon0MNegChSubNegEta ) { this->fTPCKaon0MNegChSubNegEta = TPCKaon0MNegChSubNegEta; }
  Double_t getTPCKaon0MNegChSubNegEta() { return fTPCKaon0MNegChSubNegEta; }
  void setTPCKaon3MNegChSubNegEta( Double_t TPCKaon3MNegChSubNegEta ) { this->fTPCKaon3MNegChSubNegEta = TPCKaon3MNegChSubNegEta; }
  Double_t getTPCKaon3MNegChSubNegEta() { return fTPCKaon3MNegChSubNegEta; }
  void setTPCKaon4MNegChSubNegEta( Double_t TPCKaon4MNegChSubNegEta ) { this->fTPCKaon4MNegChSubNegEta = TPCKaon4MNegChSubNegEta; }
  Double_t getTPCKaon4MNegChSubNegEta() { return fTPCKaon4MNegChSubNegEta; }
  
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
  
  void setTPCProton4Re2PosChPosEta( Double_t TPCProton4Re2PosChPosEta ) { this->fTPCProton4Re2PosChPosEta = TPCProton4Re2PosChPosEta; }
  Double_t getTPCProton4Re2PosChPosEta() { return fTPCProton4Re2PosChPosEta; }
  void setTPCProton4Im2PosChPosEta( Double_t TPCProton4Im2PosChPosEta ) { this->fTPCProton4Im2PosChPosEta = TPCProton4Im2PosChPosEta; }
  Double_t getTPCProton4Im2PosChPosEta() { return fTPCProton4Im2PosChPosEta; }
  void setTPCProton2Re3PosChPosEta( Double_t TPCProton2Re3PosChPosEta ) { this->fTPCProton2Re3PosChPosEta = TPCProton2Re3PosChPosEta; }
  Double_t getTPCProton2Re3PosChPosEta() { return fTPCProton2Re3PosChPosEta; }
  void setTPCProton2Im3PosChPosEta( Double_t TPCProton2Im3PosChPosEta ) { this->fTPCProton2Im3PosChPosEta = TPCProton2Im3PosChPosEta; }
  Double_t getTPCProton2Im3PosChPosEta() { return fTPCProton2Im3PosChPosEta; }
  void setTPCProton0MPosChPosEta( Double_t TPCProton0MPosChPosEta ) { this->fTPCProton0MPosChPosEta = TPCProton0MPosChPosEta; }
  Double_t getTPCProton0MPosChPosEta() { return fTPCProton0MPosChPosEta; }
  void setTPCProton3MPosChPosEta( Double_t TPCProton3MPosChPosEta ) { this->fTPCProton3MPosChPosEta = TPCProton3MPosChPosEta; }
  Double_t getTPCProton3MPosChPosEta() { return fTPCProton3MPosChPosEta; }
  void setTPCProton4MPosChPosEta( Double_t TPCProton4MPosChPosEta ) { this->fTPCProton4MPosChPosEta = TPCProton4MPosChPosEta; }
  Double_t getTPCProton4MPosChPosEta() { return fTPCProton4MPosChPosEta; }
  void setTPCProton4Re2PosChNegEta( Double_t TPCProton4Re2PosChNegEta ) { this->fTPCProton4Re2PosChNegEta = TPCProton4Re2PosChNegEta; }
  Double_t getTPCProton4Re2PosChNegEta() { return fTPCProton4Re2PosChNegEta; }
  void setTPCProton4Im2PosChNegEta( Double_t TPCProton4Im2PosChNegEta ) { this->fTPCProton4Im2PosChNegEta = TPCProton4Im2PosChNegEta; }
  Double_t getTPCProton4Im2PosChNegEta() { return fTPCProton4Im2PosChNegEta; }
  void setTPCProton2Re3PosChNegEta( Double_t TPCProton2Re3PosChNegEta ) { this->fTPCProton2Re3PosChNegEta = TPCProton2Re3PosChNegEta; }
  Double_t getTPCProton2Re3PosChNegEta() { return fTPCProton2Re3PosChNegEta; }
  void setTPCProton2Im3PosChNegEta( Double_t TPCProton2Im3PosChNegEta ) { this->fTPCProton2Im3PosChNegEta = TPCProton2Im3PosChNegEta; }
  Double_t getTPCProton2Im3PosChNegEta() { return fTPCProton2Im3PosChNegEta; }
  void setTPCProton0MPosChNegEta( Double_t TPCProton0MPosChNegEta ) { this->fTPCProton0MPosChNegEta = TPCProton0MPosChNegEta; }
  Double_t getTPCProton0MPosChNegEta() { return fTPCProton0MPosChNegEta; }
  void setTPCProton3MPosChNegEta( Double_t TPCProton3MPosChNegEta ) { this->fTPCProton3MPosChNegEta = TPCProton3MPosChNegEta; }
  Double_t getTPCProton3MPosChNegEta() { return fTPCProton3MPosChNegEta; }
  void setTPCProton4MPosChNegEta( Double_t TPCProton4MPosChNegEta ) { this->fTPCProton4MPosChNegEta = TPCProton4MPosChNegEta; }
  Double_t getTPCProton4MPosChNegEta() { return fTPCProton4MPosChNegEta; }
  void setTPCProton4Re2NegChPosEta( Double_t TPCProton4Re2NegChPosEta ) { this->fTPCProton4Re2NegChPosEta = TPCProton4Re2NegChPosEta; }
  Double_t getTPCProton4Re2NegChPosEta() { return fTPCProton4Re2NegChPosEta; }
  void setTPCProton4Im2NegChPosEta( Double_t TPCProton4Im2NegChPosEta ) { this->fTPCProton4Im2NegChPosEta = TPCProton4Im2NegChPosEta; }
  Double_t getTPCProton4Im2NegChPosEta() { return fTPCProton4Im2NegChPosEta; }
  void setTPCProton2Re3NegChPosEta( Double_t TPCProton2Re3NegChPosEta ) { this->fTPCProton2Re3NegChPosEta = TPCProton2Re3NegChPosEta; }
  Double_t getTPCProton2Re3NegChPosEta() { return fTPCProton2Re3NegChPosEta; }
  void setTPCProton2Im3NegChPosEta( Double_t TPCProton2Im3NegChPosEta ) { this->fTPCProton2Im3NegChPosEta = TPCProton2Im3NegChPosEta; }
  Double_t getTPCProton2Im3NegChPosEta() { return fTPCProton2Im3NegChPosEta; }
  void setTPCProton0MNegChPosEta( Double_t TPCProton0MNegChPosEta ) { this->fTPCProton0MNegChPosEta = TPCProton0MNegChPosEta; }
  Double_t getTPCProton0MNegChPosEta() { return fTPCProton0MNegChPosEta; }
  void setTPCProton3MNegChPosEta( Double_t TPCProton3MNegChPosEta ) { this->fTPCProton3MNegChPosEta = TPCProton3MNegChPosEta; }
  Double_t getTPCProton3MNegChPosEta() { return fTPCProton3MNegChPosEta; }
  void setTPCProton4MNegChPosEta( Double_t TPCProton4MNegChPosEta ) { this->fTPCProton4MNegChPosEta = TPCProton4MNegChPosEta; }
  Double_t getTPCProton4MNegChPosEta() { return fTPCProton4MNegChPosEta; }
  void setTPCProton4Re2NegChNegEta( Double_t TPCProton4Re2NegChNegEta ) { this->fTPCProton4Re2NegChNegEta = TPCProton4Re2NegChNegEta; }
  Double_t getTPCProton4Re2NegChNegEta() { return fTPCProton4Re2NegChNegEta; }
  void setTPCProton4Im2NegChNegEta( Double_t TPCProton4Im2NegChNegEta ) { this->fTPCProton4Im2NegChNegEta = TPCProton4Im2NegChNegEta; }
  Double_t getTPCProton4Im2NegChNegEta() { return fTPCProton4Im2NegChNegEta; }
  void setTPCProton2Re3NegChNegEta( Double_t TPCProton2Re3NegChNegEta ) { this->fTPCProton2Re3NegChNegEta = TPCProton2Re3NegChNegEta; }
  Double_t getTPCProton2Re3NegChNegEta() { return fTPCProton2Re3NegChNegEta; }
  void setTPCProton2Im3NegChNegEta( Double_t TPCProton2Im3NegChNegEta ) { this->fTPCProton2Im3NegChNegEta = TPCProton2Im3NegChNegEta; }
  Double_t getTPCProton2Im3NegChNegEta() { return fTPCProton2Im3NegChNegEta; }
  void setTPCProton0MNegChNegEta( Double_t TPCProton0MNegChNegEta ) { this->fTPCProton0MNegChNegEta = TPCProton0MNegChNegEta; }
  Double_t getTPCProton0MNegChNegEta() { return fTPCProton0MNegChNegEta; }
  void setTPCProton3MNegChNegEta( Double_t TPCProton3MNegChNegEta ) { this->fTPCProton3MNegChNegEta = TPCProton3MNegChNegEta; }
  Double_t getTPCProton3MNegChNegEta() { return fTPCProton3MNegChNegEta; }
  void setTPCProton4MNegChNegEta( Double_t TPCProton4MNegChNegEta ) { this->fTPCProton4MNegChNegEta = TPCProton4MNegChNegEta; }
  Double_t getTPCProton4MNegChNegEta() { return fTPCProton4MNegChNegEta; }
  
  
  void setTPCProtonRePosChSubPosEta( Double_t TPCProtonRePosChSubPosEta ) { this->fTPCProtonRePosChSubPosEta = TPCProtonRePosChSubPosEta; }
  Double_t getTPCProtonRePosChSubPosEta() { return fTPCProtonRePosChSubPosEta; }
  void setTPCProtonImPosChSubPosEta( Double_t TPCProtonImPosChSubPosEta ) { this->fTPCProtonImPosChSubPosEta = TPCProtonImPosChSubPosEta; }
  Double_t getTPCProtonImPosChSubPosEta() { return fTPCProtonImPosChSubPosEta; }
  void setTPCProtonMPosChSubPosEta( Double_t TPCProtonMPosChSubPosEta ) { this->fTPCProtonMPosChSubPosEta = TPCProtonMPosChSubPosEta; }
  Double_t getTPCProtonMPosChSubPosEta() { return fTPCProtonMPosChSubPosEta; }
  void setTPCProtonRePosChSubNegEta( Double_t TPCProtonRePosChSubNegEta ) { this->fTPCProtonRePosChSubNegEta = TPCProtonRePosChSubNegEta; }
  Double_t getTPCProtonRePosChSubNegEta() { return fTPCProtonRePosChSubNegEta; }
  void setTPCProtonImPosChSubNegEta( Double_t TPCProtonImPosChSubNegEta ) { this->fTPCProtonImPosChSubNegEta = TPCProtonImPosChSubNegEta; }
  Double_t getTPCProtonImPosChSubNegEta() { return fTPCProtonImPosChSubNegEta; }
  void setTPCProtonMPosChSubNegEta( Double_t TPCProtonMPosChSubNegEta ) { this->fTPCProtonMPosChSubNegEta = TPCProtonMPosChSubNegEta; }
  Double_t getTPCProtonMPosChSubNegEta() { return fTPCProtonMPosChSubNegEta; }
  void setTPCProtonReNegChSubPosEta( Double_t TPCProtonReNegChSubPosEta ) { this->fTPCProtonReNegChSubPosEta = TPCProtonReNegChSubPosEta; }
  Double_t getTPCProtonReNegChSubPosEta() { return fTPCProtonReNegChSubPosEta; }
  void setTPCProtonImNegChSubPosEta( Double_t TPCProtonImNegChSubPosEta ) { this->fTPCProtonImNegChSubPosEta = TPCProtonImNegChSubPosEta; }
  Double_t getTPCProtonImNegChSubPosEta() { return fTPCProtonImNegChSubPosEta; }
  void setTPCProtonMNegChSubPosEta( Double_t TPCProtonMNegChSubPosEta ) { this->fTPCProtonMNegChSubPosEta = TPCProtonMNegChSubPosEta; }
  Double_t getTPCProtonMNegChSubPosEta() { return fTPCProtonMNegChSubPosEta; }
  void setTPCProtonReNegChSubNegEta( Double_t TPCProtonReNegChSubNegEta ) { this->fTPCProtonReNegChSubNegEta = TPCProtonReNegChSubNegEta; }
  Double_t getTPCProtonReNegChSubNegEta() { return fTPCProtonReNegChSubNegEta; }
  void setTPCProtonImNegChSubNegEta( Double_t TPCProtonImNegChSubNegEta ) { this->fTPCProtonImNegChSubNegEta = TPCProtonImNegChSubNegEta; }
  Double_t getTPCProtonImNegChSubNegEta() { return fTPCProtonImNegChSubNegEta; }
  void setTPCProtonMNegChSubNegEta( Double_t TPCProtonMNegChSubNegEta ) { this->fTPCProtonMNegChSubNegEta = TPCProtonMNegChSubNegEta; }
  Double_t getTPCProtonMNegChSubNegEta() { return fTPCProtonMNegChSubNegEta; }

  void setTPCProton2RePosChSubPosEta( Double_t TPCProton2RePosChSubPosEta ) { this->fTPCProton2RePosChSubPosEta = TPCProton2RePosChSubPosEta; }
  Double_t getTPCProton2RePosChSubPosEta() { return fTPCProton2RePosChSubPosEta; }
  void setTPCProton2ImPosChSubPosEta( Double_t TPCProton2ImPosChSubPosEta ) { this->fTPCProton2ImPosChSubPosEta = TPCProton2ImPosChSubPosEta; }
  Double_t getTPCProton2ImPosChSubPosEta() { return fTPCProton2ImPosChSubPosEta; }
  void setTPCProton2Re2PosChSubPosEta( Double_t TPCProton2Re2PosChSubPosEta ) { this->fTPCProton2Re2PosChSubPosEta = TPCProton2Re2PosChSubPosEta; }
  Double_t getTPCProton2Re2PosChSubPosEta() { return fTPCProton2Re2PosChSubPosEta; }
  void setTPCProton2Im2PosChSubPosEta( Double_t TPCProton2Im2PosChSubPosEta ) { this->fTPCProton2Im2PosChSubPosEta = TPCProton2Im2PosChSubPosEta; }
  Double_t getTPCProton2Im2PosChSubPosEta() { return fTPCProton2Im2PosChSubPosEta; }
  void setTPCProton2MPosChSubPosEta( Double_t TPCProton2MPosChSubPosEta ) { this->fTPCProton2MPosChSubPosEta = TPCProton2MPosChSubPosEta; }
  Double_t getTPCProton2MPosChSubPosEta() { return fTPCProton2MPosChSubPosEta; }
  void setTPCProton2RePosChSubNegEta( Double_t TPCProton2RePosChSubNegEta ) { this->fTPCProton2RePosChSubNegEta = TPCProton2RePosChSubNegEta; }
  Double_t getTPCProton2RePosChSubNegEta() { return fTPCProton2RePosChSubNegEta; }
  void setTPCProton2ImPosChSubNegEta( Double_t TPCProton2ImPosChSubNegEta ) { this->fTPCProton2ImPosChSubNegEta = TPCProton2ImPosChSubNegEta; }
  Double_t getTPCProton2ImPosChSubNegEta() { return fTPCProton2ImPosChSubNegEta; }
  void setTPCProton2Re2PosChSubNegEta( Double_t TPCProton2Re2PosChSubNegEta ) { this->fTPCProton2Re2PosChSubNegEta = TPCProton2Re2PosChSubNegEta; }
  Double_t getTPCProton2Re2PosChSubNegEta() { return fTPCProton2Re2PosChSubNegEta; }
  void setTPCProton2Im2PosChSubNegEta( Double_t TPCProton2Im2PosChSubNegEta ) { this->fTPCProton2Im2PosChSubNegEta = TPCProton2Im2PosChSubNegEta; }
  Double_t getTPCProton2Im2PosChSubNegEta() { return fTPCProton2Im2PosChSubNegEta; }
  void setTPCProton2MPosChSubNegEta( Double_t TPCProton2MPosChSubNegEta ) { this->fTPCProton2MPosChSubNegEta = TPCProton2MPosChSubNegEta; }
  Double_t getTPCProton2MPosChSubNegEta() { return fTPCProton2MPosChSubNegEta; }
  void setTPCProton2ReNegChSubPosEta( Double_t TPCProton2ReNegChSubPosEta ) { this->fTPCProton2ReNegChSubPosEta = TPCProton2ReNegChSubPosEta; }
  Double_t getTPCProton2ReNegChSubPosEta() { return fTPCProton2ReNegChSubPosEta; }
  void setTPCProton2ImNegChSubPosEta( Double_t TPCProton2ImNegChSubPosEta ) { this->fTPCProton2ImNegChSubPosEta = TPCProton2ImNegChSubPosEta; }
  Double_t getTPCProton2ImNegChSubPosEta() { return fTPCProton2ImNegChSubPosEta; }
  void setTPCProton2Re2NegChSubPosEta( Double_t TPCProton2Re2NegChSubPosEta ) { this->fTPCProton2Re2NegChSubPosEta = TPCProton2Re2NegChSubPosEta; }
  Double_t getTPCProton2Re2NegChSubPosEta() { return fTPCProton2Re2NegChSubPosEta; }
  void setTPCProton2Im2NegChSubPosEta( Double_t TPCProton2Im2NegChSubPosEta ) { this->fTPCProton2Im2NegChSubPosEta = TPCProton2Im2NegChSubPosEta; }
  Double_t getTPCProton2Im2NegChSubPosEta() { return fTPCProton2Im2NegChSubPosEta; }
  void setTPCProton2MNegChSubPosEta( Double_t TPCProton2MNegChSubPosEta ) { this->fTPCProton2MNegChSubPosEta = TPCProton2MNegChSubPosEta; }
  Double_t getTPCProton2MNegChSubPosEta() { return fTPCProton2MNegChSubPosEta; }
  void setTPCProton2ReNegChSubNegEta( Double_t TPCProton2ReNegChSubNegEta ) { this->fTPCProton2ReNegChSubNegEta = TPCProton2ReNegChSubNegEta; }
  Double_t getTPCProton2ReNegChSubNegEta() { return fTPCProton2ReNegChSubNegEta; }
  void setTPCProton2ImNegChSubNegEta( Double_t TPCProton2ImNegChSubNegEta ) { this->fTPCProton2ImNegChSubNegEta = TPCProton2ImNegChSubNegEta; }
  Double_t getTPCProton2ImNegChSubNegEta() { return fTPCProton2ImNegChSubNegEta; }
  void setTPCProton2Re2NegChSubNegEta( Double_t TPCProton2Re2NegChSubNegEta ) { this->fTPCProton2Re2NegChSubNegEta = TPCProton2Re2NegChSubNegEta; }
  Double_t getTPCProton2Re2NegChSubNegEta() { return fTPCProton2Re2NegChSubNegEta; }
  void setTPCProton2Im2NegChSubNegEta( Double_t TPCProton2Im2NegChSubNegEta ) { this->fTPCProton2Im2NegChSubNegEta = TPCProton2Im2NegChSubNegEta; }
  Double_t getTPCProton2Im2NegChSubNegEta() { return fTPCProton2Im2NegChSubNegEta; }
  void setTPCProton2MNegChSubNegEta( Double_t TPCProton2MNegChSubNegEta ) { this->fTPCProton2MNegChSubNegEta = TPCProton2MNegChSubNegEta; }
  Double_t getTPCProton2MNegChSubNegEta() { return fTPCProton2MNegChSubNegEta; }
  
  void setTPCProton4Re2PosChSubPosEta( Double_t TPCProton4Re2PosChSubPosEta ) { this->fTPCProton4Re2PosChSubPosEta = TPCProton4Re2PosChSubPosEta; }
  Double_t getTPCProton4Re2PosChSubPosEta() { return fTPCProton4Re2PosChSubPosEta; }
  void setTPCProton4Im2PosChSubPosEta( Double_t TPCProton4Im2PosChSubPosEta ) { this->fTPCProton4Im2PosChSubPosEta = TPCProton4Im2PosChSubPosEta; }
  Double_t getTPCProton4Im2PosChSubPosEta() { return fTPCProton4Im2PosChSubPosEta; }
  void setTPCProton2Re3PosChSubPosEta( Double_t TPCProton2Re3PosChSubPosEta ) { this->fTPCProton2Re3PosChSubPosEta = TPCProton2Re3PosChSubPosEta; }
  Double_t getTPCProton2Re3PosChSubPosEta() { return fTPCProton2Re3PosChSubPosEta; }
  void setTPCProton2Im3PosChSubPosEta( Double_t TPCProton2Im3PosChSubPosEta ) { this->fTPCProton2Im3PosChSubPosEta = TPCProton2Im3PosChSubPosEta; }
  Double_t getTPCProton2Im3PosChSubPosEta() { return fTPCProton2Im3PosChSubPosEta; }
  void setTPCProton0MPosChSubPosEta( Double_t TPCProton0MPosChSubPosEta ) { this->fTPCProton0MPosChSubPosEta = TPCProton0MPosChSubPosEta; }
  Double_t getTPCProton0MPosChSubPosEta() { return fTPCProton0MPosChSubPosEta; }
  void setTPCProton3MPosChSubPosEta( Double_t TPCProton3MPosChSubPosEta ) { this->fTPCProton3MPosChSubPosEta = TPCProton3MPosChSubPosEta; }
  Double_t getTPCProton3MPosChSubPosEta() { return fTPCProton3MPosChSubPosEta; }
  void setTPCProton4MPosChSubPosEta( Double_t TPCProton4MPosChSubPosEta ) { this->fTPCProton4MPosChSubPosEta = TPCProton4MPosChSubPosEta; }
  Double_t getTPCProton4MPosChSubPosEta() { return fTPCProton4MPosChSubPosEta; }
  void setTPCProton4Re2PosChSubNegEta( Double_t TPCProton4Re2PosChSubNegEta ) { this->fTPCProton4Re2PosChSubNegEta = TPCProton4Re2PosChSubNegEta; }
  Double_t getTPCProton4Re2PosChSubNegEta() { return fTPCProton4Re2PosChSubNegEta; }
  void setTPCProton4Im2PosChSubNegEta( Double_t TPCProton4Im2PosChSubNegEta ) { this->fTPCProton4Im2PosChSubNegEta = TPCProton4Im2PosChSubNegEta; }
  Double_t getTPCProton4Im2PosChSubNegEta() { return fTPCProton4Im2PosChSubNegEta; }
  void setTPCProton2Re3PosChSubNegEta( Double_t TPCProton2Re3PosChSubNegEta ) { this->fTPCProton2Re3PosChSubNegEta = TPCProton2Re3PosChSubNegEta; }
  Double_t getTPCProton2Re3PosChSubNegEta() { return fTPCProton2Re3PosChSubNegEta; }
  void setTPCProton2Im3PosChSubNegEta( Double_t TPCProton2Im3PosChSubNegEta ) { this->fTPCProton2Im3PosChSubNegEta = TPCProton2Im3PosChSubNegEta; }
  Double_t getTPCProton2Im3PosChSubNegEta() { return fTPCProton2Im3PosChSubNegEta; }
  void setTPCProton0MPosChSubNegEta( Double_t TPCProton0MPosChSubNegEta ) { this->fTPCProton0MPosChSubNegEta = TPCProton0MPosChSubNegEta; }
  Double_t getTPCProton0MPosChSubNegEta() { return fTPCProton0MPosChSubNegEta; }
  void setTPCProton3MPosChSubNegEta( Double_t TPCProton3MPosChSubNegEta ) { this->fTPCProton3MPosChSubNegEta = TPCProton3MPosChSubNegEta; }
  Double_t getTPCProton3MPosChSubNegEta() { return fTPCProton3MPosChSubNegEta; }
  void setTPCProton4MPosChSubNegEta( Double_t TPCProton4MPosChSubNegEta ) { this->fTPCProton4MPosChSubNegEta = TPCProton4MPosChSubNegEta; }
  Double_t getTPCProton4MPosChSubNegEta() { return fTPCProton4MPosChSubNegEta; }
  void setTPCProton4Re2NegChSubPosEta( Double_t TPCProton4Re2NegChSubPosEta ) { this->fTPCProton4Re2NegChSubPosEta = TPCProton4Re2NegChSubPosEta; }
  Double_t getTPCProton4Re2NegChSubPosEta() { return fTPCProton4Re2NegChSubPosEta; }
  void setTPCProton4Im2NegChSubPosEta( Double_t TPCProton4Im2NegChSubPosEta ) { this->fTPCProton4Im2NegChSubPosEta = TPCProton4Im2NegChSubPosEta; }
  Double_t getTPCProton4Im2NegChSubPosEta() { return fTPCProton4Im2NegChSubPosEta; }
  void setTPCProton2Re3NegChSubPosEta( Double_t TPCProton2Re3NegChSubPosEta ) { this->fTPCProton2Re3NegChSubPosEta = TPCProton2Re3NegChSubPosEta; }
  Double_t getTPCProton2Re3NegChSubPosEta() { return fTPCProton2Re3NegChSubPosEta; }
  void setTPCProton2Im3NegChSubPosEta( Double_t TPCProton2Im3NegChSubPosEta ) { this->fTPCProton2Im3NegChSubPosEta = TPCProton2Im3NegChSubPosEta; }
  Double_t getTPCProton2Im3NegChSubPosEta() { return fTPCProton2Im3NegChSubPosEta; }
  void setTPCProton0MNegChSubPosEta( Double_t TPCProton0MNegChSubPosEta ) { this->fTPCProton0MNegChSubPosEta = TPCProton0MNegChSubPosEta; }
  Double_t getTPCProton0MNegChSubPosEta() { return fTPCProton0MNegChSubPosEta; }
  void setTPCProton3MNegChSubPosEta( Double_t TPCProton3MNegChSubPosEta ) { this->fTPCProton3MNegChSubPosEta = TPCProton3MNegChSubPosEta; }
  Double_t getTPCProton3MNegChSubPosEta() { return fTPCProton3MNegChSubPosEta; }
  void setTPCProton4MNegChSubPosEta( Double_t TPCProton4MNegChSubPosEta ) { this->fTPCProton4MNegChSubPosEta = TPCProton4MNegChSubPosEta; }
  Double_t getTPCProton4MNegChSubPosEta() { return fTPCProton4MNegChSubPosEta; }
  void setTPCProton4Re2NegChSubNegEta( Double_t TPCProton4Re2NegChSubNegEta ) { this->fTPCProton4Re2NegChSubNegEta = TPCProton4Re2NegChSubNegEta; }
  Double_t getTPCProton4Re2NegChSubNegEta() { return fTPCProton4Re2NegChSubNegEta; }
  void setTPCProton4Im2NegChSubNegEta( Double_t TPCProton4Im2NegChSubNegEta ) { this->fTPCProton4Im2NegChSubNegEta = TPCProton4Im2NegChSubNegEta; }
  Double_t getTPCProton4Im2NegChSubNegEta() { return fTPCProton4Im2NegChSubNegEta; }
  void setTPCProton2Re3NegChSubNegEta( Double_t TPCProton2Re3NegChSubNegEta ) { this->fTPCProton2Re3NegChSubNegEta = TPCProton2Re3NegChSubNegEta; }
  Double_t getTPCProton2Re3NegChSubNegEta() { return fTPCProton2Re3NegChSubNegEta; }
  void setTPCProton2Im3NegChSubNegEta( Double_t TPCProton2Im3NegChSubNegEta ) { this->fTPCProton2Im3NegChSubNegEta = TPCProton2Im3NegChSubNegEta; }
  Double_t getTPCProton2Im3NegChSubNegEta() { return fTPCProton2Im3NegChSubNegEta; }
  void setTPCProton0MNegChSubNegEta( Double_t TPCProton0MNegChSubNegEta ) { this->fTPCProton0MNegChSubNegEta = TPCProton0MNegChSubNegEta; }
  Double_t getTPCProton0MNegChSubNegEta() { return fTPCProton0MNegChSubNegEta; }
  void setTPCProton3MNegChSubNegEta( Double_t TPCProton3MNegChSubNegEta ) { this->fTPCProton3MNegChSubNegEta = TPCProton3MNegChSubNegEta; }
  Double_t getTPCProton3MNegChSubNegEta() { return fTPCProton3MNegChSubNegEta; }
  void setTPCProton4MNegChSubNegEta( Double_t TPCProton4MNegChSubNegEta ) { this->fTPCProton4MNegChSubNegEta = TPCProton4MNegChSubNegEta; }
  Double_t getTPCProton4MNegChSubNegEta() { return fTPCProton4MNegChSubNegEta; }
  
 private:
  Int_t fRunNum;
  Double_t fCentrality;
  Double_t fVtxPosX;
  Double_t fVtxPosY;
  Double_t fVtxPosZ;
  
  // period, orbit number, bunch cross, time stamp
  UInt_t fRawPeriod;
  UInt_t fRawOrbitNumber24;
  UInt_t fOrbitNumber;
  UShort_t fBunchCrossNumber;
  UInt_t fTimeStamp;
  
  // VZ eta < 0
  Double_t fVZCRe;
  Double_t fVZCIm;
  Double_t fVZCM;
  // VZ eta > 0
  Double_t fVZARe;
  Double_t fVZAIm;
  Double_t fVZAM;
  
  // VZ tow
  Double_t fTowV0Craw0;
  Double_t fTowV0Craw1;
  Double_t fTowV0Craw2;
  Double_t fTowV0Craw3;
  Double_t fTowV0Craw4;
  Double_t fTowV0Craw5;
  Double_t fTowV0Craw6;
  Double_t fTowV0Craw7;
  Double_t fTowV0Craw8;
  Double_t fTowV0Craw9;
  Double_t fTowV0Craw10;
  Double_t fTowV0Craw11;
  Double_t fTowV0Craw12;
  Double_t fTowV0Craw13;
  Double_t fTowV0Craw14;
  Double_t fTowV0Craw15;
  Double_t fTowV0Craw16;
  Double_t fTowV0Craw17;
  Double_t fTowV0Craw18;
  Double_t fTowV0Craw19;
  Double_t fTowV0Craw20;
  Double_t fTowV0Craw21;
  Double_t fTowV0Craw22;
  Double_t fTowV0Craw23;
  Double_t fTowV0Craw24;
  Double_t fTowV0Craw25;
  Double_t fTowV0Craw26;
  Double_t fTowV0Craw27;
  Double_t fTowV0Craw28;
  Double_t fTowV0Craw29;
  Double_t fTowV0Craw30;
  Double_t fTowV0Craw31;

  Double_t fTowV0Araw0;
  Double_t fTowV0Araw1;
  Double_t fTowV0Araw2;
  Double_t fTowV0Araw3;
  Double_t fTowV0Araw4;
  Double_t fTowV0Araw5;
  Double_t fTowV0Araw6;
  Double_t fTowV0Araw7;
  Double_t fTowV0Araw8;
  Double_t fTowV0Araw9;
  Double_t fTowV0Araw10;
  Double_t fTowV0Araw11;
  Double_t fTowV0Araw12;
  Double_t fTowV0Araw13;
  Double_t fTowV0Araw14;
  Double_t fTowV0Araw15;
  Double_t fTowV0Araw16;
  Double_t fTowV0Araw17;
  Double_t fTowV0Araw18;
  Double_t fTowV0Araw19;
  Double_t fTowV0Araw20;
  Double_t fTowV0Araw21;
  Double_t fTowV0Araw22;
  Double_t fTowV0Araw23;
  Double_t fTowV0Araw24;
  Double_t fTowV0Araw25;
  Double_t fTowV0Araw26;
  Double_t fTowV0Araw27;
  Double_t fTowV0Araw28;
  Double_t fTowV0Araw29;
  Double_t fTowV0Araw30;
  Double_t fTowV0Araw31;
  
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
  
  // TPC (Pos Eta 0.1<|eta|<0.8, Neg Eta -0.8<|eta|<-0.1, SubPos Eta 0<|eta|<0.1, SubNeg Eta -0.1<|eta|<0)
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
  
  Double_t fTPC4Re2PosChPosEta; // w^2*cos(4phi+) eta+
  Double_t fTPC4Im2PosChPosEta; // w^2*sin(4phi+) eta+
  Double_t fTPC2Re3PosChPosEta; // w^3*cos(2phi+) eta+
  Double_t fTPC2Im3PosChPosEta; // w^3*sin(2phi+) eta+
  Double_t fTPC0MPosChPosEta;   // w^0 eta+
  Double_t fTPC3MPosChPosEta;   // w^3 eta+
  Double_t fTPC4MPosChPosEta;   // w^4 eta+
  Double_t fTPC4Re2PosChNegEta; // w^2*cos(4phi+) eta-
  Double_t fTPC4Im2PosChNegEta; // w^2*sin(4phi+) eta-
  Double_t fTPC2Re3PosChNegEta; // w^3*cos(2phi+) eta-
  Double_t fTPC2Im3PosChNegEta; // w^3*sin(2phi+) eta-
  Double_t fTPC0MPosChNegEta;   // w^0 eta-
  Double_t fTPC3MPosChNegEta;   // w^3 eta-
  Double_t fTPC4MPosChNegEta;   // w^4 eta-
  Double_t fTPC4Re2NegChPosEta; // w^2*cos(4phi-) eta+
  Double_t fTPC4Im2NegChPosEta; // w^2*sin(4phi-) eta+
  Double_t fTPC2Re3NegChPosEta; // w^3*cos(2phi-) eta+
  Double_t fTPC2Im3NegChPosEta; // w^3*sin(2phi-) eta+
  Double_t fTPC0MNegChPosEta;   // w^0 eta+
  Double_t fTPC3MNegChPosEta;   // w^3 eta+
  Double_t fTPC4MNegChPosEta;   // w^4 eta+
  Double_t fTPC4Re2NegChNegEta; // w^2*cos(4phi-) eta-
  Double_t fTPC4Im2NegChNegEta; // w^2*sin(4phi-) eta-
  Double_t fTPC2Re3NegChNegEta; // w^3*cos(2phi-) eta-
  Double_t fTPC2Im3NegChNegEta; // w^3*sin(2phi-) eta-
  Double_t fTPC0MNegChNegEta;   // w^0 eta-
  Double_t fTPC3MNegChNegEta;   // w^3 eta-
  Double_t fTPC4MNegChNegEta;   // w^4 eta-
  
  Double_t fTPCRePosChSubPosEta; // w * cos(theta+) eta+
  Double_t fTPCImPosChSubPosEta; // w * sin(theta+) eta+
  Double_t fTPCMPosChSubPosEta;   // w eta+
  Double_t fTPCRePosChSubNegEta; // w * cos(theta+) eta-
  Double_t fTPCImPosChSubNegEta; // w * sin(theta+) eta-
  Double_t fTPCMPosChSubNegEta;    // w eta-
  Double_t fTPCReNegChSubPosEta; // w * cos(theta-) eta+
  Double_t fTPCImNegChSubPosEta; // w * sin(theta-) eta+
  Double_t fTPCMNegChSubPosEta;    // w eta+
  Double_t fTPCReNegChSubNegEta; // w * cos(theta-) eta-
  Double_t fTPCImNegChSubNegEta; // w * sin(theta-) eta-
  Double_t fTPCMNegChSubNegEta;    // w eta-
  
  Double_t fTPC2RePosChSubPosEta; // w * cos(2theta+) eta+
  Double_t fTPC2ImPosChSubPosEta; // w * sin(2theta+) eta+
  Double_t fTPC2Re2PosChSubPosEta; // w^2 * cos(2theta+) eta+
  Double_t fTPC2Im2PosChSubPosEta; // w^2 * sin(2theta+) eta+
  Double_t fTPC2MPosChSubPosEta;   // w^2 eta+
  Double_t fTPC2RePosChSubNegEta; // w * cos(2theta+) eta-
  Double_t fTPC2ImPosChSubNegEta; // w * sin(2theta+) eta-
  Double_t fTPC2Re2PosChSubNegEta; // w^2 * cos(2theta+) eta-
  Double_t fTPC2Im2PosChSubNegEta; // w^2 * sin(2theta+) eta-
  Double_t fTPC2MPosChSubNegEta;    // w^2 eta-
  Double_t fTPC2ReNegChSubPosEta; // w * cos(2theta-) eta+
  Double_t fTPC2ImNegChSubPosEta; // w * sin(2theta-) eta+
  Double_t fTPC2Re2NegChSubPosEta; // w^2 * cos(2theta-) eta+
  Double_t fTPC2Im2NegChSubPosEta; // w^2 * sin(2theta-) eta+
  Double_t fTPC2MNegChSubPosEta;    // w^2 eta+
  Double_t fTPC2ReNegChSubNegEta; // w * cos(2theta-) eta-
  Double_t fTPC2ImNegChSubNegEta; // w * sin(2theta-) eta-
  Double_t fTPC2Re2NegChSubNegEta; // w^2 * cos(2theta-) eta-
  Double_t fTPC2Im2NegChSubNegEta; // w^2 * sin(2theta-) eta-
  Double_t fTPC2MNegChSubNegEta;    // w^2 eta-
  
  Double_t fTPC4Re2PosChSubPosEta; // w^2*cos(4phi+) eta+
  Double_t fTPC4Im2PosChSubPosEta; // w^2*sin(4phi+) eta+
  Double_t fTPC2Re3PosChSubPosEta; // w^3*cos(2phi+) eta+
  Double_t fTPC2Im3PosChSubPosEta; // w^3*sin(2phi+) eta+
  Double_t fTPC0MPosChSubPosEta;   // w^0 eta+
  Double_t fTPC3MPosChSubPosEta;   // w^3 eta+
  Double_t fTPC4MPosChSubPosEta;   // w^4 eta+
  Double_t fTPC4Re2PosChSubNegEta; // w^2*cos(4phi+) eta-
  Double_t fTPC4Im2PosChSubNegEta; // w^2*sin(4phi+) eta-
  Double_t fTPC2Re3PosChSubNegEta; // w^3*cos(2phi+) eta-
  Double_t fTPC2Im3PosChSubNegEta; // w^3*sin(2phi+) eta-
  Double_t fTPC0MPosChSubNegEta;   // w^0 eta-
  Double_t fTPC3MPosChSubNegEta;   // w^3 eta-
  Double_t fTPC4MPosChSubNegEta;   // w^4 eta-
  Double_t fTPC4Re2NegChSubPosEta; // w^2*cos(4phi-) eta+
  Double_t fTPC4Im2NegChSubPosEta; // w^2*sin(4phi-) eta+
  Double_t fTPC2Re3NegChSubPosEta; // w^3*cos(2phi-) eta+
  Double_t fTPC2Im3NegChSubPosEta; // w^3*sin(2phi-) eta+
  Double_t fTPC0MNegChSubPosEta;   // w^0 eta+
  Double_t fTPC3MNegChSubPosEta;   // w^3 eta+
  Double_t fTPC4MNegChSubPosEta;   // w^4 eta+
  Double_t fTPC4Re2NegChSubNegEta; // w^2*cos(4phi-) eta-
  Double_t fTPC4Im2NegChSubNegEta; // w^2*sin(4phi-) eta-
  Double_t fTPC2Re3NegChSubNegEta; // w^3*cos(2phi-) eta-
  Double_t fTPC2Im3NegChSubNegEta; // w^3*sin(2phi-) eta-
  Double_t fTPC0MNegChSubNegEta;   // w^0 eta-
  Double_t fTPC3MNegChSubNegEta;   // w^3 eta-
  Double_t fTPC4MNegChSubNegEta;   // w^4 eta-
  
  // TPC Pion (Pos Eta 0.1<|eta|<0.8, Neg Eta -0.8<|eta|<-0.1, SubPos Eta 0<|eta|<0.1, SubNeg Eta -0.1<|eta|<0)
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
  
  Double_t fTPCPion4Re2PosChPosEta; // w^2*cos(4phi+) eta+
  Double_t fTPCPion4Im2PosChPosEta; // w^2*sin(4phi+) eta+
  Double_t fTPCPion2Re3PosChPosEta; // w^3*cos(2phi+) eta+
  Double_t fTPCPion2Im3PosChPosEta; // w^3*sin(2phi+) eta+
  Double_t fTPCPion0MPosChPosEta;   // w^0 eta+
  Double_t fTPCPion3MPosChPosEta;   // w^3 eta+
  Double_t fTPCPion4MPosChPosEta;   // w^4 eta+
  Double_t fTPCPion4Re2PosChNegEta; // w^2*cos(4phi+) eta-
  Double_t fTPCPion4Im2PosChNegEta; // w^2*sin(4phi+) eta-
  Double_t fTPCPion2Re3PosChNegEta; // w^3*cos(2phi+) eta-
  Double_t fTPCPion2Im3PosChNegEta; // w^3*sin(2phi+) eta-
  Double_t fTPCPion0MPosChNegEta;   // w^0 eta-
  Double_t fTPCPion3MPosChNegEta;   // w^3 eta-
  Double_t fTPCPion4MPosChNegEta;   // w^4 eta-
  Double_t fTPCPion4Re2NegChPosEta; // w^2*cos(4phi-) eta+
  Double_t fTPCPion4Im2NegChPosEta; // w^2*sin(4phi-) eta+
  Double_t fTPCPion2Re3NegChPosEta; // w^3*cos(2phi-) eta+
  Double_t fTPCPion2Im3NegChPosEta; // w^3*sin(2phi-) eta+
  Double_t fTPCPion0MNegChPosEta;   // w^0 eta+
  Double_t fTPCPion3MNegChPosEta;   // w^3 eta+
  Double_t fTPCPion4MNegChPosEta;   // w^4 eta+
  Double_t fTPCPion4Re2NegChNegEta; // w^2*cos(4phi-) eta-
  Double_t fTPCPion4Im2NegChNegEta; // w^2*sin(4phi-) eta-
  Double_t fTPCPion2Re3NegChNegEta; // w^3*cos(2phi-) eta-
  Double_t fTPCPion2Im3NegChNegEta; // w^3*sin(2phi-) eta-
  Double_t fTPCPion0MNegChNegEta;   // w^0 eta-
  Double_t fTPCPion3MNegChNegEta;   // w^3 eta-
  Double_t fTPCPion4MNegChNegEta;   // w^4 eta-
  
  Double_t fTPCPionRePosChSubPosEta; // w * cos(theta+) eta+
  Double_t fTPCPionImPosChSubPosEta; // w * sin(theta+) eta+
  Double_t fTPCPionMPosChSubPosEta;   // w eta+
  Double_t fTPCPionRePosChSubNegEta; // w * cos(theta+) eta-
  Double_t fTPCPionImPosChSubNegEta; // w * sin(theta+) eta-
  Double_t fTPCPionMPosChSubNegEta;    // w eta-
  Double_t fTPCPionReNegChSubPosEta; // w * cos(theta-) eta+
  Double_t fTPCPionImNegChSubPosEta; // w * sin(theta-) eta+
  Double_t fTPCPionMNegChSubPosEta;    // w eta+
  Double_t fTPCPionReNegChSubNegEta; // w * cos(theta-) eta-
  Double_t fTPCPionImNegChSubNegEta; // w * sin(theta-) eta-
  Double_t fTPCPionMNegChSubNegEta;    // w eta-
  
  Double_t fTPCPion2RePosChSubPosEta; // w * cos(2theta+) eta+
  Double_t fTPCPion2ImPosChSubPosEta; // w * sin(2theta+) eta+
  Double_t fTPCPion2Re2PosChSubPosEta; // w^2 * cos(2theta+) eta+
  Double_t fTPCPion2Im2PosChSubPosEta; // w^2 * sin(2theta+) eta+
  Double_t fTPCPion2MPosChSubPosEta;   // w^2 eta+
  Double_t fTPCPion2RePosChSubNegEta; // w * cos(2theta+) eta-
  Double_t fTPCPion2ImPosChSubNegEta; // w * sin(2theta+) eta-
  Double_t fTPCPion2Re2PosChSubNegEta; // w^2 * cos(2theta+) eta-
  Double_t fTPCPion2Im2PosChSubNegEta; // w^2 * sin(2theta+) eta-
  Double_t fTPCPion2MPosChSubNegEta;    // w^2 eta-
  Double_t fTPCPion2ReNegChSubPosEta; // w * cos(2theta-) eta+
  Double_t fTPCPion2ImNegChSubPosEta; // w * sin(2theta-) eta+
  Double_t fTPCPion2Re2NegChSubPosEta; // w^2 * cos(2theta-) eta+
  Double_t fTPCPion2Im2NegChSubPosEta; // w^2 * sin(2theta-) eta+
  Double_t fTPCPion2MNegChSubPosEta;    // w^2 eta+
  Double_t fTPCPion2ReNegChSubNegEta; // w * cos(2theta-) eta-
  Double_t fTPCPion2ImNegChSubNegEta; // w * sin(2theta-) eta-
  Double_t fTPCPion2Re2NegChSubNegEta; // w^2 * cos(2theta-) eta-
  Double_t fTPCPion2Im2NegChSubNegEta; // w^2 * sin(2theta-) eta-
  Double_t fTPCPion2MNegChSubNegEta;    // w^2 eta-
  
  Double_t fTPCPion4Re2PosChSubPosEta; // w^2*cos(4phi+) eta+
  Double_t fTPCPion4Im2PosChSubPosEta; // w^2*sin(4phi+) eta+
  Double_t fTPCPion2Re3PosChSubPosEta; // w^3*cos(2phi+) eta+
  Double_t fTPCPion2Im3PosChSubPosEta; // w^3*sin(2phi+) eta+
  Double_t fTPCPion0MPosChSubPosEta;   // w^0 eta+
  Double_t fTPCPion3MPosChSubPosEta;   // w^3 eta+
  Double_t fTPCPion4MPosChSubPosEta;   // w^4 eta+
  Double_t fTPCPion4Re2PosChSubNegEta; // w^2*cos(4phi+) eta-
  Double_t fTPCPion4Im2PosChSubNegEta; // w^2*sin(4phi+) eta-
  Double_t fTPCPion2Re3PosChSubNegEta; // w^3*cos(2phi+) eta-
  Double_t fTPCPion2Im3PosChSubNegEta; // w^3*sin(2phi+) eta-
  Double_t fTPCPion0MPosChSubNegEta;   // w^0 eta-
  Double_t fTPCPion3MPosChSubNegEta;   // w^3 eta-
  Double_t fTPCPion4MPosChSubNegEta;   // w^4 eta-
  Double_t fTPCPion4Re2NegChSubPosEta; // w^2*cos(4phi-) eta+
  Double_t fTPCPion4Im2NegChSubPosEta; // w^2*sin(4phi-) eta+
  Double_t fTPCPion2Re3NegChSubPosEta; // w^3*cos(2phi-) eta+
  Double_t fTPCPion2Im3NegChSubPosEta; // w^3*sin(2phi-) eta+
  Double_t fTPCPion0MNegChSubPosEta;   // w^0 eta+
  Double_t fTPCPion3MNegChSubPosEta;   // w^3 eta+
  Double_t fTPCPion4MNegChSubPosEta;   // w^4 eta+
  Double_t fTPCPion4Re2NegChSubNegEta; // w^2*cos(4phi-) eta-
  Double_t fTPCPion4Im2NegChSubNegEta; // w^2*sin(4phi-) eta-
  Double_t fTPCPion2Re3NegChSubNegEta; // w^3*cos(2phi-) eta-
  Double_t fTPCPion2Im3NegChSubNegEta; // w^3*sin(2phi-) eta-
  Double_t fTPCPion0MNegChSubNegEta;   // w^0 eta-
  Double_t fTPCPion3MNegChSubNegEta;   // w^3 eta-
  Double_t fTPCPion4MNegChSubNegEta;   // w^4 eta-
  
  // TPC Kaon (Pos Eta 0.1<|eta|<0.8, Neg Eta -0.8<|eta|<-0.1, SubPos Eta 0<|eta|<0.1, SubNeg Eta -0.1<|eta|<0)
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
  
  Double_t fTPCKaon4Re2PosChPosEta; // w^2*cos(4phi+) eta+
  Double_t fTPCKaon4Im2PosChPosEta; // w^2*sin(4phi+) eta+
  Double_t fTPCKaon2Re3PosChPosEta; // w^3*cos(2phi+) eta+
  Double_t fTPCKaon2Im3PosChPosEta; // w^3*sin(2phi+) eta+
  Double_t fTPCKaon0MPosChPosEta;   // w^0 eta+
  Double_t fTPCKaon3MPosChPosEta;   // w^3 eta+
  Double_t fTPCKaon4MPosChPosEta;   // w^4 eta+
  Double_t fTPCKaon4Re2PosChNegEta; // w^2*cos(4phi+) eta-
  Double_t fTPCKaon4Im2PosChNegEta; // w^2*sin(4phi+) eta-
  Double_t fTPCKaon2Re3PosChNegEta; // w^3*cos(2phi+) eta-
  Double_t fTPCKaon2Im3PosChNegEta; // w^3*sin(2phi+) eta-
  Double_t fTPCKaon0MPosChNegEta;   // w^0 eta-
  Double_t fTPCKaon3MPosChNegEta;   // w^3 eta-
  Double_t fTPCKaon4MPosChNegEta;   // w^4 eta-
  Double_t fTPCKaon4Re2NegChPosEta; // w^2*cos(4phi-) eta+
  Double_t fTPCKaon4Im2NegChPosEta; // w^2*sin(4phi-) eta+
  Double_t fTPCKaon2Re3NegChPosEta; // w^3*cos(2phi-) eta+
  Double_t fTPCKaon2Im3NegChPosEta; // w^3*sin(2phi-) eta+
  Double_t fTPCKaon0MNegChPosEta;   // w^0 eta+
  Double_t fTPCKaon3MNegChPosEta;   // w^3 eta+
  Double_t fTPCKaon4MNegChPosEta;   // w^4 eta+
  Double_t fTPCKaon4Re2NegChNegEta; // w^2*cos(4phi-) eta-
  Double_t fTPCKaon4Im2NegChNegEta; // w^2*sin(4phi-) eta-
  Double_t fTPCKaon2Re3NegChNegEta; // w^3*cos(2phi-) eta-
  Double_t fTPCKaon2Im3NegChNegEta; // w^3*sin(2phi-) eta-
  Double_t fTPCKaon0MNegChNegEta;   // w^0 eta-
  Double_t fTPCKaon3MNegChNegEta;   // w^3 eta-
  Double_t fTPCKaon4MNegChNegEta;   // w^4 eta-
  
  Double_t fTPCKaonRePosChSubPosEta; // w * cos(theta+) eta+
  Double_t fTPCKaonImPosChSubPosEta; // w * sin(theta+) eta+
  Double_t fTPCKaonMPosChSubPosEta;   // w eta+
  Double_t fTPCKaonRePosChSubNegEta; // w * cos(theta+) eta-
  Double_t fTPCKaonImPosChSubNegEta; // w * sin(theta+) eta-
  Double_t fTPCKaonMPosChSubNegEta;    // w eta-
  Double_t fTPCKaonReNegChSubPosEta; // w * cos(theta-) eta+
  Double_t fTPCKaonImNegChSubPosEta; // w * sin(theta-) eta+
  Double_t fTPCKaonMNegChSubPosEta;    // w eta+
  Double_t fTPCKaonReNegChSubNegEta; // w * cos(theta-) eta-
  Double_t fTPCKaonImNegChSubNegEta; // w * sin(theta-) eta-
  Double_t fTPCKaonMNegChSubNegEta;    // w eta-
  
  Double_t fTPCKaon2RePosChSubPosEta; // w * cos(2theta+) eta+
  Double_t fTPCKaon2ImPosChSubPosEta; // w * sin(2theta+) eta+
  Double_t fTPCKaon2Re2PosChSubPosEta; // w^2 * cos(2theta+) eta+
  Double_t fTPCKaon2Im2PosChSubPosEta; // w^2 * sin(2theta+) eta+
  Double_t fTPCKaon2MPosChSubPosEta;   // w^2 eta+
  Double_t fTPCKaon2RePosChSubNegEta; // w * cos(2theta+) eta-
  Double_t fTPCKaon2ImPosChSubNegEta; // w * sin(2theta+) eta-
  Double_t fTPCKaon2Re2PosChSubNegEta; // w^2 * cos(2theta+) eta-
  Double_t fTPCKaon2Im2PosChSubNegEta; // w^2 * sin(2theta+) eta-
  Double_t fTPCKaon2MPosChSubNegEta;    // w^2 eta-
  Double_t fTPCKaon2ReNegChSubPosEta; // w * cos(2theta-) eta+
  Double_t fTPCKaon2ImNegChSubPosEta; // w * sin(2theta-) eta+
  Double_t fTPCKaon2Re2NegChSubPosEta; // w^2 * cos(2theta-) eta+
  Double_t fTPCKaon2Im2NegChSubPosEta; // w^2 * sin(2theta-) eta+
  Double_t fTPCKaon2MNegChSubPosEta;    // w^2 eta+
  Double_t fTPCKaon2ReNegChSubNegEta; // w * cos(2theta-) eta-
  Double_t fTPCKaon2ImNegChSubNegEta; // w * sin(2theta-) eta-
  Double_t fTPCKaon2Re2NegChSubNegEta; // w^2 * cos(2theta-) eta-
  Double_t fTPCKaon2Im2NegChSubNegEta; // w^2 * sin(2theta-) eta-
  Double_t fTPCKaon2MNegChSubNegEta;    // w^2 eta-
  
  Double_t fTPCKaon4Re2PosChSubPosEta; // w^2*cos(4phi+) eta+
  Double_t fTPCKaon4Im2PosChSubPosEta; // w^2*sin(4phi+) eta+
  Double_t fTPCKaon2Re3PosChSubPosEta; // w^3*cos(2phi+) eta+
  Double_t fTPCKaon2Im3PosChSubPosEta; // w^3*sin(2phi+) eta+
  Double_t fTPCKaon0MPosChSubPosEta;   // w^0 eta+
  Double_t fTPCKaon3MPosChSubPosEta;   // w^3 eta+
  Double_t fTPCKaon4MPosChSubPosEta;   // w^4 eta+
  Double_t fTPCKaon4Re2PosChSubNegEta; // w^2*cos(4phi+) eta-
  Double_t fTPCKaon4Im2PosChSubNegEta; // w^2*sin(4phi+) eta-
  Double_t fTPCKaon2Re3PosChSubNegEta; // w^3*cos(2phi+) eta-
  Double_t fTPCKaon2Im3PosChSubNegEta; // w^3*sin(2phi+) eta-
  Double_t fTPCKaon0MPosChSubNegEta;   // w^0 eta-
  Double_t fTPCKaon3MPosChSubNegEta;   // w^3 eta-
  Double_t fTPCKaon4MPosChSubNegEta;   // w^4 eta-
  Double_t fTPCKaon4Re2NegChSubPosEta; // w^2*cos(4phi-) eta+
  Double_t fTPCKaon4Im2NegChSubPosEta; // w^2*sin(4phi-) eta+
  Double_t fTPCKaon2Re3NegChSubPosEta; // w^3*cos(2phi-) eta+
  Double_t fTPCKaon2Im3NegChSubPosEta; // w^3*sin(2phi-) eta+
  Double_t fTPCKaon0MNegChSubPosEta;   // w^0 eta+
  Double_t fTPCKaon3MNegChSubPosEta;   // w^3 eta+
  Double_t fTPCKaon4MNegChSubPosEta;   // w^4 eta+
  Double_t fTPCKaon4Re2NegChSubNegEta; // w^2*cos(4phi-) eta-
  Double_t fTPCKaon4Im2NegChSubNegEta; // w^2*sin(4phi-) eta-
  Double_t fTPCKaon2Re3NegChSubNegEta; // w^3*cos(2phi-) eta-
  Double_t fTPCKaon2Im3NegChSubNegEta; // w^3*sin(2phi-) eta-
  Double_t fTPCKaon0MNegChSubNegEta;   // w^0 eta-
  Double_t fTPCKaon3MNegChSubNegEta;   // w^3 eta-
  Double_t fTPCKaon4MNegChSubNegEta;   // w^4 eta-
  
  // TPC Proton (0.1<|eta|<0.8)
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
  
  Double_t fTPCProton4Re2PosChPosEta; // w^2*cos(4phi+) eta+
  Double_t fTPCProton4Im2PosChPosEta; // w^2*sin(4phi+) eta+
  Double_t fTPCProton2Re3PosChPosEta; // w^3*cos(2phi+) eta+
  Double_t fTPCProton2Im3PosChPosEta; // w^3*sin(2phi+) eta+
  Double_t fTPCProton0MPosChPosEta;   // w^0 eta+
  Double_t fTPCProton3MPosChPosEta;   // w^3 eta+
  Double_t fTPCProton4MPosChPosEta;   // w^4 eta+
  Double_t fTPCProton4Re2PosChNegEta; // w^2*cos(4phi+) eta-
  Double_t fTPCProton4Im2PosChNegEta; // w^2*sin(4phi+) eta-
  Double_t fTPCProton2Re3PosChNegEta; // w^3*cos(2phi+) eta-
  Double_t fTPCProton2Im3PosChNegEta; // w^3*sin(2phi+) eta-
  Double_t fTPCProton0MPosChNegEta;   // w^0 eta-
  Double_t fTPCProton3MPosChNegEta;   // w^3 eta-
  Double_t fTPCProton4MPosChNegEta;   // w^4 eta-
  Double_t fTPCProton4Re2NegChPosEta; // w^2*cos(4phi-) eta+
  Double_t fTPCProton4Im2NegChPosEta; // w^2*sin(4phi-) eta+
  Double_t fTPCProton2Re3NegChPosEta; // w^3*cos(2phi-) eta+
  Double_t fTPCProton2Im3NegChPosEta; // w^3*sin(2phi-) eta+
  Double_t fTPCProton0MNegChPosEta;   // w^0 eta+
  Double_t fTPCProton3MNegChPosEta;   // w^3 eta+
  Double_t fTPCProton4MNegChPosEta;   // w^4 eta+
  Double_t fTPCProton4Re2NegChNegEta; // w^2*cos(4phi-) eta-
  Double_t fTPCProton4Im2NegChNegEta; // w^2*sin(4phi-) eta-
  Double_t fTPCProton2Re3NegChNegEta; // w^3*cos(2phi-) eta-
  Double_t fTPCProton2Im3NegChNegEta; // w^3*sin(2phi-) eta-
  Double_t fTPCProton0MNegChNegEta;   // w^0 eta-
  Double_t fTPCProton3MNegChNegEta;   // w^3 eta-
  Double_t fTPCProton4MNegChNegEta;   // w^4 eta-
  
  Double_t fTPCProtonRePosChSubPosEta; // w * cos(theta+) eta+
  Double_t fTPCProtonImPosChSubPosEta; // w * sin(theta+) eta+
  Double_t fTPCProtonMPosChSubPosEta;   // w eta+
  Double_t fTPCProtonRePosChSubNegEta; // w * cos(theta+) eta-
  Double_t fTPCProtonImPosChSubNegEta; // w * sin(theta+) eta-
  Double_t fTPCProtonMPosChSubNegEta;    // w eta-
  Double_t fTPCProtonReNegChSubPosEta; // w * cos(theta-) eta+
  Double_t fTPCProtonImNegChSubPosEta; // w * sin(theta-) eta+
  Double_t fTPCProtonMNegChSubPosEta;    // w eta+
  Double_t fTPCProtonReNegChSubNegEta; // w * cos(theta-) eta-
  Double_t fTPCProtonImNegChSubNegEta; // w * sin(theta-) eta-
  Double_t fTPCProtonMNegChSubNegEta;    // w eta-
  
  Double_t fTPCProton2RePosChSubPosEta; // w * cos(2theta+) eta+
  Double_t fTPCProton2ImPosChSubPosEta; // w * sin(2theta+) eta+
  Double_t fTPCProton2Re2PosChSubPosEta; // w^2 * cos(2theta+) eta+
  Double_t fTPCProton2Im2PosChSubPosEta; // w^2 * sin(2theta+) eta+
  Double_t fTPCProton2MPosChSubPosEta;   // w^2 eta+
  Double_t fTPCProton2RePosChSubNegEta; // w * cos(2theta+) eta-
  Double_t fTPCProton2ImPosChSubNegEta; // w * sin(2theta+) eta-
  Double_t fTPCProton2Re2PosChSubNegEta; // w^2 * cos(2theta+) eta-
  Double_t fTPCProton2Im2PosChSubNegEta; // w^2 * sin(2theta+) eta-
  Double_t fTPCProton2MPosChSubNegEta;    // w^2 eta-
  Double_t fTPCProton2ReNegChSubPosEta; // w * cos(2theta-) eta+
  Double_t fTPCProton2ImNegChSubPosEta; // w * sin(2theta-) eta+
  Double_t fTPCProton2Re2NegChSubPosEta; // w^2 * cos(2theta-) eta+
  Double_t fTPCProton2Im2NegChSubPosEta; // w^2 * sin(2theta-) eta+
  Double_t fTPCProton2MNegChSubPosEta;    // w^2 eta+
  Double_t fTPCProton2ReNegChSubNegEta; // w * cos(2theta-) eta-
  Double_t fTPCProton2ImNegChSubNegEta; // w * sin(2theta-) eta-
  Double_t fTPCProton2Re2NegChSubNegEta; // w^2 * cos(2theta-) eta-
  Double_t fTPCProton2Im2NegChSubNegEta; // w^2 * sin(2theta-) eta-
  Double_t fTPCProton2MNegChSubNegEta;    // w^2 eta-
  
  Double_t fTPCProton4Re2PosChSubPosEta; // w^2*cos(4phi+) eta+
  Double_t fTPCProton4Im2PosChSubPosEta; // w^2*sin(4phi+) eta+
  Double_t fTPCProton2Re3PosChSubPosEta; // w^3*cos(2phi+) eta+
  Double_t fTPCProton2Im3PosChSubPosEta; // w^3*sin(2phi+) eta+
  Double_t fTPCProton0MPosChSubPosEta;   // w^0 eta+
  Double_t fTPCProton3MPosChSubPosEta;   // w^3 eta+
  Double_t fTPCProton4MPosChSubPosEta;   // w^4 eta+
  Double_t fTPCProton4Re2PosChSubNegEta; // w^2*cos(4phi+) eta-
  Double_t fTPCProton4Im2PosChSubNegEta; // w^2*sin(4phi+) eta-
  Double_t fTPCProton2Re3PosChSubNegEta; // w^3*cos(2phi+) eta-
  Double_t fTPCProton2Im3PosChSubNegEta; // w^3*sin(2phi+) eta-
  Double_t fTPCProton0MPosChSubNegEta;   // w^0 eta-
  Double_t fTPCProton3MPosChSubNegEta;   // w^3 eta-
  Double_t fTPCProton4MPosChSubNegEta;   // w^4 eta-
  Double_t fTPCProton4Re2NegChSubPosEta; // w^2*cos(4phi-) eta+
  Double_t fTPCProton4Im2NegChSubPosEta; // w^2*sin(4phi-) eta+
  Double_t fTPCProton2Re3NegChSubPosEta; // w^3*cos(2phi-) eta+
  Double_t fTPCProton2Im3NegChSubPosEta; // w^3*sin(2phi-) eta+
  Double_t fTPCProton0MNegChSubPosEta;   // w^0 eta+
  Double_t fTPCProton3MNegChSubPosEta;   // w^3 eta+
  Double_t fTPCProton4MNegChSubPosEta;   // w^4 eta+
  Double_t fTPCProton4Re2NegChSubNegEta; // w^2*cos(4phi-) eta-
  Double_t fTPCProton4Im2NegChSubNegEta; // w^2*sin(4phi-) eta-
  Double_t fTPCProton2Re3NegChSubNegEta; // w^3*cos(2phi-) eta-
  Double_t fTPCProton2Im3NegChSubNegEta; // w^3*sin(2phi-) eta-
  Double_t fTPCProton0MNegChSubNegEta;   // w^0 eta-
  Double_t fTPCProton3MNegChSubNegEta;   // w^3 eta-
  Double_t fTPCProton4MNegChSubNegEta;   // w^4 eta-
  
  ClassDef(AliAnalysisTaskGammaDeltaPIDSaveQvecEvent, 1);
};

#endif
