/************************************************************************************
 * Copyright (C) 2020, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <memory>
#include <iostream>
#include "AliAnalysisManager.h"
#include "AliEmcalMCPartonInfo.h"
#include "AliEmcalMCPartonInfoCreator.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliVEvent.h"

ClassImp(PWG::EMCAL::AliEmcalMCPartonInfoCreator)

using namespace PWG::EMCAL;

AliEmcalMCPartonInfoCreator::AliEmcalMCPartonInfoCreator():
  AliAnalysisTaskSE(),
  fMCPartonInfo(nullptr),
  fInitialized(false),
  fNamePartonInfo("EmcalMCPartonInfo")
{

}

AliEmcalMCPartonInfoCreator::AliEmcalMCPartonInfoCreator(const char *name):
  AliAnalysisTaskSE(name),
  fMCPartonInfo(nullptr),
  fInitialized(false),
  fNamePartonInfo("EmcalMCPartonInfo")
{

}

AliEmcalMCPartonInfoCreator::~AliEmcalMCPartonInfoCreator(){
  if(fMCPartonInfo) delete fMCPartonInfo;
}

void AliEmcalMCPartonInfoCreator::UserExec(Option_t * /*option*/){
  if(!fInitialized) {
    if(!(InputEvent()->FindListObject(fNamePartonInfo.Data()))) {
      fMCPartonInfo = new AliEmcalMCPartonInfo(fNamePartonInfo.Data());
      InputEvent()->AddObject(fMCPartonInfo);
      fInitialized = true;
    } else {
      fInitialized = false;
      AliFatal(Form("%s Container with same name %s already present. Aborting", GetName(), fNamePartonInfo.Data())); 
      return;
    }
  }

  if(!fInitialized) return;
  fMCPartonInfo->Reset(); 

  // Build tree in order to find all direct daughters of the beam particles
  struct ParticleNode {
    AliVParticle *fNodeParticle;
    Int_t fParticleIndex;
    std::vector<std::shared_ptr<ParticleNode>> fDirectDaughters;

    ParticleNode *FindNode(Int_t motherID) {
      if(motherID == fParticleIndex) return this;
          for(auto daughter : fDirectDaughters) {
            auto result = daughter->FindNode(motherID);
            if(result) return result;
          }
      return nullptr;
    }

    bool Insert(AliVParticle *part, int index) {
      auto mothernode = FindNode(part->GetMother());
      if(!mothernode){
        return false;
      }
      auto daughterParticle = new ParticleNode;
      daughterParticle->fNodeParticle = part;
      daughterParticle->fParticleIndex = index;
      mothernode->fDirectDaughters.push_back(std::shared_ptr<ParticleNode>(daughterParticle));
      return true;
    }

    int getLength() {
      int max = 0;
      for(auto daughter : fDirectDaughters){
        auto daughterlength = daughter->getLength();
        if(daughterlength > max) max = daughterlength;
      }
      return max+1;
    }
  };

  // Build tree starting from first beam particle
  ParticleNode mothernode;
  mothernode.fNodeParticle = MCEvent()->GetTrack(0);
  mothernode.fParticleIndex = 0;
  for(auto ipart = 2; ipart < MCEvent()->GetNumberOfTracks(); ipart++) {
    mothernode.Insert(MCEvent()->GetTrack(ipart), ipart);
  }

  // Adding direct daughters to the parton info object
  for(auto partnode : mothernode.fDirectDaughters) {
    auto part = partnode.get()->fNodeParticle;
    if(part->PdgCode() < 7 || part->PdgCode() == 21) fMCPartonInfo->AddDirectParton(part, partnode->fParticleIndex);
  }

}

AliEmcalMCPartonInfoCreator *AliEmcalMCPartonInfoCreator::AddTaskMCPartonInfoCreator(const char *taskname) {
  auto mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    AliErrorGeneralStream("AliEmcalMCPartonInfoCreator::AddTaskMCPartonInfoCreator") << "Could not get analysis manager" << std::endl;
    return nullptr;
  }

  auto task = new AliEmcalMCPartonInfoCreator(taskname);
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  return task;
}
