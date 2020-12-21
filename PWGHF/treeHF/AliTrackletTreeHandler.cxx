/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
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

/**
 * \class AliTrackletTreeHandler
 * \brief Helper class to handle a tree for cut optimisation and MVA analyses, based heavily on AliHFTreeHandler and AliAnalysisTaskDmesonJets
 *
 * \author Nima Zardoshti nima.zardoshti@cern.ch
 *         Benedikt Volkel benedikt.volkel@cern.ch
 *         Cristina Terrevoli cristina.terrevoli@cern.ch
 *         Gian Michele Innocenti ginnocen@cern.ch
 * \date June 27 2019
 */

#include "AliTrackletTreeHandler.h"
#include <iostream>

//________________________________________________________________
/// \cond CLASSIMP
ClassImp(AliTrackletTreeHandler);
/// \endcond

//________________________________________________________________
// Default constructor
AliTrackletTreeHandler::AliTrackletTreeHandler():
  TObject(),
  fTreeTracklet(nullptr),
  fTrackletContainer(nullptr),
  fTrackletEta(-999.),
  fTrackletPhi(-999.),
  fRunNumber(0),
  fEventID(0), 
  fEventIDExt(0),
  fEventIDLong(0)
{
}

//________________________________________________________________
// Destructor
AliTrackletTreeHandler::~AliTrackletTreeHandler()
{
  if(fTreeTracklet) delete fTreeTracklet;
}

/**
 * Create particle TTree, with completely flat structure.
 */
//________________________________________________________________
TTree* AliTrackletTreeHandler::BuildTree(TString name, TString title)
{
  if(fTreeTracklet) {
    delete fTreeTracklet;
    fTreeTracklet=0x0;
  }
  fTreeTracklet = new TTree(name.Data(),title.Data());
  
  // Create branches for each particle variable
  fTreeTracklet->Branch("run_number", &fRunNumber);
  fTreeTracklet->Branch("ev_id",&fEventID);
  fTreeTracklet->Branch("ev_id_ext",&fEventIDExt);
  fTreeTracklet->Branch("ev_id_long",&fEventIDLong);
  fTreeTracklet->Branch("TrackletEta",&fTrackletEta);
  fTreeTracklet->Branch("TrackletPhi",&fTrackletPhi);
  
  return fTreeTracklet;
}

/**
 * Set tree variables and fill them
 */
//________________________________________________________________
void AliTrackletTreeHandler::FillTree(int runNumber, int eventID, int eventID_Ext, Long64_t eventID_Long)
{
  
  fRunNumber = runNumber;
  fEventID = eventID;
  fEventIDExt = eventID_Ext;
  fEventIDLong = eventID_Long;
  Int_t nTr=fTrackletContainer->GetNumberOfTracklets();
  for(Int_t iTr=0; iTr<nTr; iTr++){
    Double_t phi=fTrackletContainer->GetPhi(iTr);
    Double_t theta=fTrackletContainer->GetTheta(iTr);
    Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
    fTrackletEta = eta;
    fTrackletPhi = phi;
    fTreeTracklet->Fill();
  } 
}
