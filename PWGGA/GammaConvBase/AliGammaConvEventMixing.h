/**************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Joshua Koenig <joshua.konig@cern.ch>                                        *
 * Version 1.0                                                            *
 *                                                                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do event mixing for events with Jet axis
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////

#ifndef ALIGAMMACONVEVENTMIXING_H
#define ALIGAMMACONVEVENTMIXING_H

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <map>
#include <vector>
#include "AliAODConversionPhoton.h"


struct EventWithJetAxis{
  EventWithJetAxis(std::vector<AliAODConversionPhoton*> vecPhotons, bool isCalo, TVector3 jet) : caloPhotons({}), convPhotons({}), jetAxis() {
    if(isCalo){
      for(const auto & i : vecPhotons){
        caloPhotons.push_back(new AliAODConversionPhoton(*i));
      }
    } else {
      for(const auto & i : vecPhotons){
        convPhotons.push_back(new AliAODConversionPhoton(*i));
      }
    }
    jetAxis = jet;
  }
  EventWithJetAxis(std::vector<AliAODConversionPhoton*> vecCaloPhotons, std::vector<AliAODConversionPhoton*> vecConvPhotons, TVector3 jet) : caloPhotons({}), convPhotons({}), jetAxis(){
    for(const auto & i : vecCaloPhotons){
      caloPhotons.push_back(new AliAODConversionPhoton(*i));
    }
    for(const auto & i : vecConvPhotons){
      convPhotons.push_back(new AliAODConversionPhoton(*i));
    }
    jetAxis = jet;
  }
  ~EventWithJetAxis(){
    for(unsigned int i = 0; i < caloPhotons.size(); ++i){
      if(caloPhotons[i] != nullptr) delete caloPhotons[i];
    }
    for(unsigned int i = 0; i < convPhotons.size(); ++i){
      if(convPhotons[i] != nullptr) delete convPhotons[i];
    }
  }
  AliAODConversionPhoton* getPhoton(int index, bool isCaloPhoton){
    if(isCaloPhoton) return caloPhotons[index];
    return convPhotons[index];
  }
  unsigned int getNPhotons(bool isCaloPhoton){
    if(isCaloPhoton) return caloPhotons.size();
    return convPhotons.size();
  }
  std::vector<AliAODConversionPhoton*> caloPhotons;
  std::vector<AliAODConversionPhoton*> convPhotons;
  TVector3 jetAxis;
};

class EventMixPoolMesonJets{
  public:
  EventMixPoolMesonJets();
  EventMixPoolMesonJets(std::vector<float> vec);
  EventMixPoolMesonJets(const EventMixPoolMesonJets&) = delete;
  EventMixPoolMesonJets(EventMixPoolMesonJets&&) = default;
  ~EventMixPoolMesonJets();

  ///\brief Set new jet momentum classes for mixing pool
  ///\param vec vector with jet momenta
  void SetJetPtClasses(std::vector<float> vec);

  ///\brief Set mixing pool depth
  ///\param n new pool depth
  void SetPoolDepth(unsigned int n)             { poolDepth = n;}

  ///\brief get index of jet momentum for mixing pool
  ///\param jetP momentum of the jet
  ///\return index
  int getJetPIndex(float jetP);

  ///\brief Add new event to the mixing pool. If mixing pool is full, delete first element and then add event
  ///\param ev event that should be added
  ///\param jetP highest jet momentum of this event
  void AddEvent(EventWithJetAxis *ev, float jetP);

  ///\brief get number of events in mixing pool for specified jet momentum
  ///\param jetP momentum of the jet
  ///\return number of events in mixing pool for specified jet momentum
  unsigned int GetNBGEvents(float jetP);
  
  ///\brief get number of photons in the specified event in the mixing pool
  ///\param jetP momentum of the jet in the current event
  ///\param nEvt index of the event in the mixing pool 
  ///\param isCaloPhoton switch if calo photons or conversion photons should be returned
  ///\return number of photons in the specified mixing pool class
  unsigned int GetNGammasInEvt(float jetP, int evt, bool isCaloPhoton);

  ///\brief rotate events such that the jet axes with the highest momentum match.
  ///\param nEvt index of the event in the mixing pool 
  ///\param jetP momentum of the jet in the current event
  ///\param jetAxis Jet axis of the current event
  ///\param isCaloPhoton switch if calo photons or conversion photons should be returned
  ///\return vector of rotated photons that can be used for mixing with th current event
  std::vector<std::unique_ptr<AliAODConversionPhoton>> getPhotonsRotated(int nEvt, float jetP, TVector3 jetAxis, bool isCaloPhoton);
  
  private:
  std::vector<float> vecJetPClasses = {-1, 20, 40, 60, 100, 100000};
  std::vector<std::vector<EventWithJetAxis*>> mixingPool;
  unsigned int poolDepth = 20;

};

#endif
