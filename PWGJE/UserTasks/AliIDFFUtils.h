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

//Utils for identified fragmentation function (IDFF) analysis
//Author: Xianguo Lu (xianguo.lu@cern.ch)

#ifndef ALIIDFFUTILS_H
#define ALIIDFFUTILS_H

#include "TH2D.h"
#include "THnSparse.h"
#include "AliPIDResponse.h"
#include "AliAODTrack.h"

#include <TTreeStream.h>
  
class AliIDFFUtils
{

 friend class AliAnalysisTaskIDFFTCF;

 public:
  enum TYPEID{
    kNOTACCEPTED = -3, 
    kNOINFO = -2,
    kNOTSELECTED = -1,
    kPROTON=0,
    kPION,
    kKAON,
    kELECTRON
  };

  static THnSparseD* GetTHn(const TString name);
  static void FillTHn(THnSparseD * hh, Double_t jetpt, const AliAODTrack * trackk,  AliAODEvent* aodevt, Int_t tofmode);  

  static Bool_t TPCCutPIDN(const AliAODTrack * track);
  static Bool_t TPCCutMIGeo(const AliAODTrack * track, const AliVEvent* evt, TTreeStream * streamer=0x0);

                                                                        
  static AliPIDResponse * fPid;

 private:
  static Double_t Xmin(){return -1;}
  static Double_t Xmax(){return 2;}
  static Int_t Nx(){return 300;}

  static Int_t PDG2Type(Int_t pdg);
  
  static Int_t TOFType(const AliAODTrack * trackptr, Int_t tofmode);
  static Int_t TPCType(const AliAODTrack * trackptr);

  static Int_t CombineTPCTOF(Int_t ktpc, Int_t ktof);

  //static Bool_t HMPIDAcceptance(const AliAODTrack *track);
  //static Bool_t HMPIDQA(const AliAODTrack *track);
  //static Int_t HMPIDType(const AliAODTrack * track);

};

#endif
