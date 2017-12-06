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

#include "AliEMCALTriggerDCSConfig.h"
#include "AliEMCALTriggerSTUDCSConfig.h"
#include "AliEMCALTriggerTRUDCSConfig.h"

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerDCSConfig) ;
/// \endcond

///
/// Default constructor
//_____________________________________________________________________________
AliEMCALTriggerDCSConfig::AliEMCALTriggerDCSConfig() : TObject()
,fTRUArr(0x0)
,fSTUObj(0x0)
,fSTUDCAL(0x0)
{
  /*
	fTRUArr = new TClonesArray("AliEMCALTriggerTRUDCSConfig",62);
	fSTUObj = new AliEMCALTriggerSTUDCSConfig();
	*/
}

///
/// Destructor
//_____________________________________________________________________________
AliEMCALTriggerDCSConfig::~AliEMCALTriggerDCSConfig()
{	
  delete fTRUArr;
  delete fSTUObj;
  delete fSTUDCAL;
}

bool AliEMCALTriggerDCSConfig::operator==(const AliEMCALTriggerDCSConfig & other) const {
  bool isequal = true;
  if(fSTUObj && other.fSTUObj) {
    if(!(*fSTUObj == *other.fSTUObj)) isequal == false; // both EMCAL STU objects there, but not the same
  } else if((fSTUObj && !other.fSTUObj) || (!fSTUObj && other.fSTUObj)) isequal == false; // one of the two missing

  if(fSTUDCAL && other.fSTUDCAL) {
    if(!(*fSTUDCAL == *other.fSTUDCAL)) isequal == false; // both DCAL STU objects there, but not the same
  } else if((fSTUDCAL && !other.fSTUDCAL) || (!fSTUDCAL && other.fSTUDCAL)) isequal == false; // one of the two missing

  // check TRUs
  if(fTRUArr->GetEntries() != other.fTRUArr->GetEntries()){
    // one object has more TRU entries than the other
    isequal = false;
  } else {
    for(int itru = 0; itru < fTRUArr->GetSize(); itru++) {
      AliEMCALTriggerTRUDCSConfig *thistru = GetTRUDCSConfig(itru),
                                  *othertru = other.GetTRUDCSConfig(itru);
      if(thistru && othertru){
        if(!(*thistru == *othertru)) isequal = false;     
      } else if((thistru && !othertru) || (!thistru && othertru)) isequal = false;    // one of the two missing
    }
  }
  return isequal;
}