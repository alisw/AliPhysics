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
  
  ClassDef(AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple, 1);
};

#endif
