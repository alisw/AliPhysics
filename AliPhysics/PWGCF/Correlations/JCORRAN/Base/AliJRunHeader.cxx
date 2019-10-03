/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

// Comment describing what this class does needed!

// $Id: AliJRunHeader.cxx,v 1.1 2008/02/04 13:28:47 rak Exp $

////////////////////////////////////////////////////
//
//  \file AliJRunHeader.cxx
//  \brief
//  \author J. Rak, D.J.Kim, F.Krizek(University of Jyvaskyla)
//  \email: djkim@cc.jyu.fi
//  \version $Revision: 1.1 $
//  \date $Date: 2008/02/04 13:28:47 $
//
// Class encapsulation aliroot run header information
////////////////////////////////////////////////////
#include "AliJRunHeader.h"
#include <iostream>
#include "AliJConst.h"

ClassImp(AliJRunHeader);

//----------------------------
AliJRunHeader::AliJRunHeader():
TNamed("AliJRunHeader", ""),
//==== General Info ====//
fRunNumber(0),
fRunType(""),
fBeamType(""),
fBeamTypeI(-1),
fBeamEnergy(0),
fIsMC(kFALSE),
//==== Production Info ====//
fInputFormat(-1),
fWithoutSDD(kFALSE),
fStoreEventPlaneSource(kFALSE),
fStoreEMCalInfo(kFALSE),
fStoreTPCTrackBitMask(0),
fStoreGCGTrackBitMask(0),
fESDInfo(""),
fRefitESDVertexTracks(false),
//==== Detector Status ====//
fL3MagnetPolarity(0),
fMagneticFieldL3(0),
fCurrentL3(0),
fCurrentDip(0),
fUniformBMap(kFALSE),
//==== Trigger Information ====//
fFiredTriggers(0),
fTriggerMask(0),
fTriggerCluster(0),
fSizeOfTableJCorran(0),
fActiveTriggersJCorran(kRangeTriggerTableJCorran,""),
fActiveTriggersAlice(kRangeTriggerTableAlice,"")
{                     //constructor
}

Int_t AliJRunHeader::GetActiveTriggerBitAlice(TString TriggerName){
  //get trigger bit corresponding to trigger name
  for(UInt_t t=0;t<GetActiveTriggersAlice().size();t++){
    if(TriggerName.Contains( GetActiveTriggersAlice(t) )){ // TODO : not equal?
      return t;
    }
  }
  return -1;
}

//________________________________________________________________________

void AliJRunHeader::PrintOut(){
  //print object
  std::cout<<"RUN         = "<<fRunNumber<<std::endl;
  std::cout<<"RunType     = "<<fRunType<<std::endl;
  std::cout<<"BeamType    = "<<fBeamType<<" => "<<fBeamTypeI<<std::endl;
  std::cout<<"BeamEnergy  = "<<fBeamEnergy<<std::endl;
  std::cout<<"IsMC        = "<<fIsMC<<std::endl;
  std::cout<<"ImputFormat = "<<fInputFormat<<std::endl;
  std::cout<<"Alice trigger table:"<<std::endl;
  for(Int_t i=0;i<kRangeTriggerTableAlice;i++){
    std::cout<<i<<"   "<<GetActiveTriggersAlice(i)<<std::endl;
  }
  std::cout<<"============================="<<std::endl;
  std::cout<<"JCorran trigger table:"<<std::endl;
  for(Int_t i=0;i<kRangeTriggerTableJCorran;i++){
    std::cout<<i<<"   "<<GetActiveTriggersJCorran(i)<<std::endl;
  }
  std::cout<<"============================="<<std::endl;
  std::cout<<"Magnet polarity "<<fL3MagnetPolarity<<std::endl;
  std::cout<<"B "<<fMagneticFieldL3<<std::endl;
  std::cout<<"============================="<<std::endl;

}
