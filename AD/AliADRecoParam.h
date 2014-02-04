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


/**************************************************************************************

        AD's RecoParams Version 1.0

        In this version we only consider:

        ->) The AD's Trigger Mode (Single Muon Trigger or Multi Muon Trigger)
        ->) The AD's Trigger Mask (Same in SMT and MMT)
        
        In Runs PbPb, pp, and cosmics by default we have the same RecoParams.

        From: 
		Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch> @ CERN
        	FCFM, BUAP. Puebla, Mexico
        
	Further comments:

		Arturo Fernandez <afernan@mail.cern.ch>

        March 2nd. 2009 

        NOTE: Please suggest improvements if needed.

**************************************************************************************/


#ifndef ALIADRECOPARAM_H
#define ALIADRECOPARAM_H

#include "AliDetectorRecoParam.h"

class AliADRecoParam : public AliDetectorRecoParam
{
 public: 
  AliADRecoParam();
  AliADRecoParam(const AliADRecoParam &p); // Copy constructor 
  AliADRecoParam& operator=(const AliADRecoParam &p); // Assignment operator
  virtual ~AliADRecoParam();

  virtual void PrintParameters() const;

  //Getters
/* 
  Bool_t   GetAcordeSingleMuonTrigger()    const  {return fAcordeSingleMuonTrigger;}
  Bool_t   GetAcordeMultiMuonTrigger()     const  {return fAcordeMultiMuonTrigger;}
  UInt_t   GetAcordeWord0()    const  {return fAcordeWord0;}
  UInt_t   GetAcordeWord1()    const  {return fAcordeWord1;}
  UInt_t   GetAcordeWord2()    const  {return fAcordeWord2;}
  UInt_t   GetAcordeWord3()    const  {return fAcordeWord3;}


  //Setters

  void   SetAcordeSingleMuonTrigger(Bool_t flag)        {fAcordeSingleMuonTrigger=flag;}
  void   SetAcordeMultiMuonTrigger(Bool_t flag)  {fAcordeMultiMuonTrigger=flag;}
  void   SetAcordeWord0(UInt_t flag)  {fAcordeWord0=flag;}
  void   SetAcordeWord1(UInt_t flag)  {fAcordeWord1=flag;}
  void   SetAcordeWord2(UInt_t flag)  {fAcordeWord2=flag;}
  void   SetAcordeWord3(UInt_t flag)  {fAcordeWord3=flag;}



  static   AliADRecoParam *GetPbPbparam();       // reco param for PbPb.
  static   AliADRecoParam *GetPPparam();         // reco param for PP
  static   AliADRecoParam *GetCosmicMuonParam(); // reco param for cosmic muons
 
private:

  Bool_t fAcordeSingleMuonTrigger; // kTRUE if AD triggered in Singe Muon Mode
  Bool_t fAcordeMultiMuonTrigger; // kTRUE if AD triggered in Multi Muon Mode
  UInt_t fAcordeWord0; // [1..30] Acorde's Modules in Single Muon Trigger
  UInt_t fAcordeWord1; // [31..60] Acorde's Modules in Single Muon Trigger
  UInt_t fAcordeWord2; // [1..30] Acorde's Modules in Multi Muon Trigger
  UInt_t fAcordeWord3; // [31..60] Acorde's Modules in Multi Muon Trigger
*/
  ClassDef(AliADRecoParam, 1)
};
#endif
