/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Author: Svein Lindal <slindal@fys.uio.no>                      *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliAnalysisTaskGammaJet.cxx
/// @author Svein Lindal
/// @brief  Class used to run conversion gamma/pion - hadron/jet analysis


#include <iostream>
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TObjArray.h"


#include "TH2F.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskGammaJet.h"

#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDInputHandler.h"

#include "AliAODPWG4ParticleCorrelation.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODCaloCluster.h"
#include "AliAODConversionPhoton.h"
#include "AliAODJet.h"

#include "AliAODInputHandler.h"

#include "AliAnaConvIsolation.h"

#include "AliAnaConvCorrPhoton.h"
#include "AliAnaConvCorrPhotonJet.h"
#include "AliAnaConvCorrPionJet.h"
#include "AliAnaConvCorrPion.h"

#include "AliAODTrack.h"

// Gamma - jet correlation analysis task
// Authors: Svein Lindal


using namespace std;

ClassImp(AliAnalysisTaskGammaJet)

//________________________________________________________________________
AliAnalysisTaskGammaJet::AliAnalysisTaskGammaJet() : 
AliAnalysisTaskSE(), 
  fOutputList(NULL), 
  fEtaLimit(0.8),
  fDeltaAODFileName("AliAODGammaConversion.root"),
  fGammaCutString("GammaConv"),
  fPionCutString("GammaConv"),
  fAnaIsolation(NULL), 
  fAnaIsolationArray(NULL), 
  fAnaPionArray(NULL), 
  fAnaPhotonArray(NULL),
  fAnaPhotonJetArray(NULL),
  fAnaPionJetArray(NULL),
  fMinPt(1.0), 
  fMinNTracks(20)
{
  // Dummy Constructor
}

//________________________________________________________________________________
AliAnalysisTaskGammaJet::~AliAnalysisTaskGammaJet() {
  //Destructor
}


//________________________________________________________________________
AliAnalysisTaskGammaJet::AliAnalysisTaskGammaJet(const char *name) : 
  AliAnalysisTaskSE(name), 
  fOutputList(0), 
  fEtaLimit(0.8),
  fDeltaAODFileName("AliAODGammaConversion.root"),
  fGammaCutString("GammaConv"),
  fPionCutString("GammaConv"),
  fAnaIsolation(NULL),
  fAnaIsolationArray(NULL),
  fAnaPionArray(NULL),
  fAnaPhotonArray(NULL),
  fAnaPhotonJetArray(NULL), 
  fAnaPionJetArray(NULL),
  fMinPt(1.0),
  fMinNTracks(20)
{
  // Constructor
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());

  fAnaIsolation = new AliAnaConvIsolation();

  fAnaPionArray = new TObjArray();
  fAnaPionArray->SetOwner(kTRUE);

  fAnaPhotonArray = new TObjArray();
  fAnaPhotonArray->SetOwner(kTRUE);

  fAnaPhotonJetArray = new TObjArray();
  fAnaPhotonJetArray->SetOwner(kTRUE);

  fAnaPionJetArray = new TObjArray();
  fAnaPionJetArray->SetOwner(kTRUE);

  fAnaIsolationArray = new TObjArray();
  fAnaIsolationArray->SetOwner(kTRUE);

  fhTracksMissingPt[0] = NULL;
  fhTracksMissingPt[1] = NULL;



}

//________________________________________________________________________
void AliAnalysisTaskGammaJet::UserExec(Option_t *) 
{
  //Inherited from AliAnalysisTaskSE

  ///Get AOD event
  AliAODEvent * aodEvent = GetAODEvent();
  if(!aodEvent) {
    cout << "No AOD event!!" << endl;
    return;
  }

  TClonesArray * photons = GetConversionGammas(aodEvent);
  if(!photons) {
    cout << "No Conversion gamma!!" << endl;
    return;
  }

  TClonesArray * pions = GetConversionPions(aodEvent);
  if(!pions) {
    cout << "No Conversion pion branch!!" << endl;
    return;
  }

  TClonesArray * tracks = aodEvent->GetTracks();
  if(!tracks) {
    cout << "Can't get tracks!!" << endl;
    return;
  }

  if(!((Entry() )%10000)) {
    AliInfo(Form("%s ----> Processing event # %lld", CurrentFileName(), Entry()));
    AliInfo(Form("%d %d", photons->GetEntriesFast(), pions->GetEntriesFast()));
  }


  if(pions->GetEntriesFast() > photons->GetEntriesFast() ) {
    cout << "WTF!!!!"<<endl;
  }

 // if( (photons->GetEntriesFast() > 0) ) {
 //   cout << "Neext evetn !!!!!!!!!!!!!!!!"  << Entry() << endl;
 // } else {
 //   cout << "____________________________"<<endl;
 // }


  if(aodEvent->GetNumberOfTracks() < 2 && photons->GetEntriesFast() > 0 ) {
    cout << "we have a photon but less than 2 tracks" << endl;
    return;
  }

  if(aodEvent->GetNumberOfTracks() < fMinNTracks) return;

  if( (photons->GetEntriesFast() > 0) ) {
    ProcessConvGamma(photons, pions, tracks);
    ProcessPions(pions, photons, tracks);
  }



  TClonesArray * jets = aodEvent->GetJets();
  if(jets && jets->GetEntriesFast() > 0) {
    for(int i = 0; i < fAnaPhotonJetArray->GetEntriesFast(); i++) {
      AliAnaConvCorrPhotonJet * jetAna = dynamic_cast<AliAnaConvCorrPhotonJet*>(fAnaPhotonJetArray->At(i));
      if(jetAna) {
	for(Int_t iJet = 0; iJet < jets->GetEntriesFast(); iJet++) {
	  AliAODJet * jet = dynamic_cast<AliAODJet*>(jets->At(iJet));
	  if(jet) {
	    jetAna->DoJetAnalysisGamma(jet, photons, pions);
	  }		   
	}
      }
    }
  }
  PostData(1, fOutputList);

}



//________________________________________________________________________
void AliAnalysisTaskGammaJet::UserCreateOutputObjects()
{
  //Create histograms add, to outputlist
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);


  fhTracksMissingPt[0] = new TH2F("hpttracksmissing", "hpttracksmissing",200,  0, 200, 200, 0, 200);
  fOutputList->Add(fhTracksMissingPt[0]);
  fhTracksMissingPt[1] = new TH2F("hpttrackpresent", "hptbothtrackspresent" ,200,  0, 200, 200, 0, 200);
  fOutputList->Add(fhTracksMissingPt[1]);

  fAnaIsolation->CreateHistograms();
  fOutputList->Add(fAnaIsolation->GetHistograms());

  for(int i = 0; i < fAnaIsolationArray->GetEntriesFast(); i++) {
    AliAnaConvIsolation * isoAna = dynamic_cast<AliAnaConvIsolation*>(fAnaIsolationArray->At(i));
    if(isoAna) {
      isoAna->CreateHistograms();
      fOutputList->Add(isoAna->GetHistograms());
    } else {
      AliError("problem getting ana pointer!!!");
      cout << i << endl;
    }
  }

  for(int i = 0; i < fAnaPhotonArray->GetEntriesFast(); i++) {
    AliAnaConvCorrPhoton * photonAna = dynamic_cast<AliAnaConvCorrPhoton*>(fAnaPhotonArray->At(i));
    if(photonAna) {
      photonAna->CreateHistograms();
      fOutputList->Add(photonAna->GetHistograms());
    } else {
      AliError("problem getting ana pointer!!!");
    }
  }


  for(int i = 0; i < fAnaPionArray->GetEntriesFast(); i++) {
    AliAnaConvCorrPion * pionAna = dynamic_cast<AliAnaConvCorrPion*>(fAnaPionArray->At(i));
    if(pionAna) {
      pionAna->CreateHistograms();
      fOutputList->Add(pionAna->GetHistograms());
    } else {
      AliError("problem getting ana pointer!!!");
    }
  }

  for(int i = 0; i < fAnaPhotonJetArray->GetEntriesFast(); i++) {
    AliAnaConvCorrPhotonJet * jetAna = dynamic_cast<AliAnaConvCorrPhotonJet*>(fAnaPhotonJetArray->At(i));
    if(jetAna) {
      jetAna->CreateHistograms();
      fOutputList->Add(jetAna->GetHistograms());
    } else {
      AliError("problem getting jet  ana pointer!!!");
    }
  }

  for(int i = 0; i < fAnaPionJetArray->GetEntriesFast(); i++) {
    AliAnaConvCorrPionJet * jetAna = dynamic_cast<AliAnaConvCorrPionJet*>(fAnaPionJetArray->At(i));
    if(jetAna) {
      jetAna->CreateHistograms();
      fOutputList->Add(jetAna->GetHistograms());
    } else {
      AliError("problem getting jet  ana pointer!!!");
    }
  }

  

  PostData(1, fOutputList);


}





//______________________________________________________________________________________________
void AliAnalysisTaskGammaJet::ProcessCalorimeters( const AliAODEvent * const aodEvent ) {
  //See header file for documentation
  TClonesArray * clusters = aodEvent->GetCaloClusters();
  

  for(int ic = 0; ic < clusters->GetEntriesFast(); ic++) {
    AliAODCaloCluster * cluster = dynamic_cast<AliAODCaloCluster*>(clusters->At(ic));
    if (!cluster) { 
      AliError(Form("Error getting cluster"));
      continue;
    }


    if (cluster->GetNCells() < 6) continue;
    if (cluster->GetEmcCpvDistance() < 15) continue;

    TLorentzVector tlvec;
    
    AliAODVertex * vertex = aodEvent->GetPrimaryVertex();
    Double_t vertexPosition[3];
    vertex->GetXYZ(vertexPosition);
    cluster->GetMomentum(tlvec, vertexPosition);
    //if (tlvec.Pt() < GetMinPt()) continue; 
    
  }
  
}

//______________________________________________________________________________________________
Bool_t AliAnalysisTaskGammaJet::EventIsSynced(const TClonesArray * const tracks, const TClonesArray * const convGamma, const TClonesArray * const pions)  {
 //See header file for documentation

  for (Int_t iPhot = 0; iPhot < convGamma->GetEntriesFast(); iPhot++) {
    AliAODConversionPhoton * photon = dynamic_cast<AliAODConversionPhoton*>(convGamma->At(iPhot));
    if(photon) {
      AliAODTrack * track1 = NULL;
      AliAODTrack * track2 = NULL;
      for(Int_t i = 0; i < tracks->GetEntriesFast(); i++) {
	AliAODTrack * track = dynamic_cast<AliAODTrack*>(tracks->At(i));
	if(track) {
	  if(track->GetID() == photon->GetLabel1()) track1 = track;
	  else if (track->GetID() == photon->GetLabel2()) track2 = track;
	  if(track1 && track2) break;
	}
      }
      
      if(track1 && track2) {
	Float_t totE = track1->E() + track2->E();
	if(TMath::Abs(totE - photon->E()) > 1.0)  {
	  cout << "BALLE BALLE "<<TMath::Abs(totE - photon->P()) << endl;

	  cout <<  track2->Px() << " "
	       <<  track2->Py() << " "
	       <<  track2->Pz() << " "
	       <<  track2->P() << " "
	       <<  track2->E() << endl;

	  cout << track1->Px() << " "
	       << track1->Py() << " "
	       << track1->Pz() << " "
	       << track1->P()  << " "
	       << track1->E()  << endl;


	  cout << track1->Px() + track2->Px() << " "
	       << track1->Py() + track2->Py() << " "
	       << track1->Pz() + track2->Pz() << " "
	       << track1->P() + track2->P() << " "
	       << track1->E() + track2->E() << endl;
	  
	  cout << photon->Px() << " " <<  photon->Py() << " " <<  photon->Pz() << " " << photon->P() << " " <<  photon->E() << endl;
	  return kFALSE;
	}
      } else {
	cout << Entry() << " " << convGamma->GetEntriesFast() << " " << photon->GetLabel1() << " " << photon->GetLabel2() << endl;
	cout << "Could not get both tracks!!! " <<endl;
	return kFALSE;
      }
    }
  }

  if(pions) { 
    //placeholder
  }
  
  //cout  <<"insync"<<endl;
  return kTRUE;
}


//______________________________________________________________________________________________
Bool_t AliAnalysisTaskGammaJet::BothTracksPresent(const AliAODConversionPhoton * const photon, const TClonesArray * const tracks)  const {

  AliAODTrack * track1 = NULL;
  AliAODTrack * track2 = NULL;
  for(Int_t i = 0; i < tracks->GetEntriesFast(); i++) {
    AliAODTrack * track = dynamic_cast<AliAODTrack*>(tracks->At(i));
    if(track) {
      if(track->GetID() == photon->GetLabel1()) track1 = track;
      else if (track->GetID() == photon->GetLabel2()) track2 = track;
      if(track1 && track2) break;
    }
  }
  
  if(track1 && track2) {
    return kTRUE;
  }
  cout << "Could not get both tracks!!! labels "  << photon->GetLabel1() << " " << photon->GetLabel2()  <<endl;
  return kFALSE;
  

}

//______________________________________________________________________________________________
Bool_t AliAnalysisTaskGammaJet::BothGammaPresent(const AliAODConversionPhoton * const pion, const TClonesArray * const photons, const TClonesArray * const tracks)  const {

  AliAODConversionPhoton * photon1 = dynamic_cast<AliAODConversionPhoton*>(photons->At(pion->GetLabel1()));
  AliAODConversionPhoton * photon2 = dynamic_cast<AliAODConversionPhoton*>(photons->At(pion->GetLabel2()));

  if(photon1 && photon2) {
    if( BothTracksPresent(photon1, tracks) &&  BothTracksPresent(photon1, tracks)) {
      return kTRUE;
    }
  } else {
    cout << "can't find both photons "<< endl;
  }
  
  return kFALSE;
}




//___________________________________________________________________________________________
void AliAnalysisTaskGammaJet::ProcessConvGamma( const TClonesArray * convGamma, const TClonesArray * const pions, const TClonesArray * tracks ) {
  //See header file for documentation


  for (Int_t iPhot = 0; iPhot < convGamma->GetEntriesFast(); iPhot++) {
    AliAODConversionPhoton * photon = dynamic_cast<AliAODConversionPhoton*>(convGamma->At(iPhot));
    if(!photon) {
      AliError(Form("ERROR: Could not receive ga %d\n", iPhot));
      continue;
    }
    
    Bool_t btp = BothTracksPresent(photon, tracks);
    fhTracksMissingPt[btp]->Fill(photon->Pt(), tracks->GetEntriesFast());
    if(!btp || photon->Pt() < fMinPt || TMath::Abs(photon->Eta()) > fEtaLimit) {
      continue;
    }
    
    Bool_t leading = kTRUE;
    Bool_t decayPion = IsDecayPion(iPhot, pions);

    for(Int_t i = 0; i < fAnaIsolationArray->GetEntriesFast(); i ++) {
      AliAnaConvIsolation * isoAna = dynamic_cast<AliAnaConvIsolation*>(fAnaIsolationArray->At(i));
      if(isoAna)  isoAna->IsIsolated(photon, tracks, leading);
    }

    Bool_t isolated = fAnaIsolation->IsIsolated(photon, tracks, leading);
    if(leading) {
      for(Int_t i = 0; i < fAnaPhotonArray->GetEntriesFast(); i ++) {
	AliAnaConvCorrPhoton * ana = static_cast<AliAnaConvCorrPhoton*>(fAnaPhotonArray->At(i));
	if(ana) {
	  ana->CorrelateWithHadrons(photon, tracks, isolated, decayPion);
	} 
      }
 
      TClonesArray * jets = GetAODEvent()->GetJets();
      if(jets) {
	for(Int_t i = 0; i < fAnaPhotonJetArray->GetEntriesFast(); i ++) {
	  AliAnaConvCorrPhotonJet * ana = static_cast<AliAnaConvCorrPhotonJet*>(fAnaPhotonJetArray->At(i));
	  if(ana) {
	    ana->CorrelateWithHadrons(photon, jets, isolated);
	  } 
	}
      } else {
	cout << "No jets "<<endl;
      }
    } 
  } // 
}

///______________________________________________________________________________________________
void AliAnalysisTaskGammaJet::ProcessPions( const TClonesArray * const pions, const TClonesArray * const photons, const TClonesArray * const tracks) {
  //See header file for documentation


  for (Int_t iPi = 0; iPi < pions->GetEntriesFast(); iPi++) {
    AliAODConversionPhoton * pion = dynamic_cast<AliAODConversionPhoton*>(pions->At(iPi));
    if(!pion) {
      AliError(Form("ERROR: Could not receive ga %d\n", iPi));
      continue;
    }
    
    if (!BothGammaPresent(pion, photons, tracks) || pion->Pt() < fMinPt || TMath::Abs(pion->Eta()) > fEtaLimit ) {
      return;
    }

    
    Int_t trackLabels[4] = {-1, -1, -1, -1};
    GetPionGrandChildren(pion, photons, trackLabels);
    Bool_t leading = kTRUE;

    for(Int_t i = 0; i < fAnaIsolationArray->GetEntriesFast(); i ++) {
      AliAnaConvIsolation * isoAna = dynamic_cast<AliAnaConvIsolation*>(fAnaIsolationArray->At(i));
      if(isoAna) isoAna->IsIsolated(pion, tracks, 4, trackLabels, leading);
    }

    Bool_t isolated = fAnaIsolation->IsIsolated(pion, tracks, 4, trackLabels, leading);
    if(leading) {
      for(Int_t i = 0; i < fAnaPionArray->GetEntriesFast(); i ++) {
  	AliAnaConvCorrPion * ana = dynamic_cast<AliAnaConvCorrPion*>(fAnaPionArray->At(i));
  	if(ana) ana->CorrelateWithHadrons(pion, tracks, isolated, 4, trackLabels );
      }


      TClonesArray * jets = GetAODEvent()->GetJets();
      if(jets) {
	for(Int_t i = 0; i < fAnaPionJetArray->GetEntriesFast(); i ++) {
	  AliAnaConvCorrPionJet * ana = static_cast<AliAnaConvCorrPionJet*>(fAnaPionJetArray->At(i));
	  if(ana) {
	    ana->CorrelateWithHadrons(pion, jets, isolated);
	  } 
	}
      } else {
	cout << "No jets "<<endl;
      }
      
      
    } 
  } // 
 

  
  // for (Int_t iPhot = 0; iPhot < photons->GetEntriesFast(); iPhot++) {
  //   AliAODConversionPhoton * photon = dynamic_cast<AliAODConversionPhoton*>(photons->At(iPhot));
  //   if(photon) {
  //     for (Int_t iPhot2 = iPhot+1; iPhot2 < photons->GetEntriesFast(); iPhot2++) {
  // 	AliAODConversionPhoton * photon2 = dynamic_cast<AliAODConversionPhoton*>(photons->At(iPhot2));
  // 	if(photon2) {
  // 	  Int_t trackLabels[4] = {photon->GetTrackLabel(0), photon->GetTrackLabel(1), photon2->GetTrackLabel(0), photon2->GetTrackLabel(1)};
  // 	  AliAODConversionPhoton * pion = new AliAODConversionPhoton(photon, photon2);
  // 	  Bool_t leading = kTRUE;
  // 	  Bool_t isolated = fAnaIsolation->IsIsolated(pion, tracks, 4, trackLabels, leading);
  // 	  if(leading) {
  // 	    for(Int_t i = 0; i < fAnaPionArray->GetEntriesFast(); i ++) {
  // 	      AliAnaConvCorrPion * ana = dynamic_cast<AliAnaConvCorrPion*>(fAnaPionArray->At(i));
  // 	      if(ana) ana->CorrelateWithHadrons(pion, tracks, isolated, 4, trackLabels );
  // 	    }
  // 	    delete pion;
  // 	  }
  // 	}
  //     }
  //   }
  // }
}

///_______________________________________________________________________________________________
Bool_t AliAnalysisTaskGammaJet::IsDecayPion(Int_t iPhot, const TClonesArray * const pions) {
  //See header file for documentation

  for(int ip = 0; ip < pions->GetEntriesFast(); ip++) {
    AliAODConversionPhoton * pion = dynamic_cast<AliAODConversionPhoton*>(pions->At(ip));
    if(pion) {
      if(pion->GetLabel1() == iPhot || pion->GetLabel2() == iPhot) 
	
	return kTRUE;
    } else {
      AliError("pion corrupted!");
    }
  }

  return kFALSE;
}


///_______________________________________________________________________________
Bool_t AliAnalysisTaskGammaJet::UserNotify() {
    //See header file for documentation

    AliInfo(Form("%s ----> Processing event # %lld", CurrentFileName(), Entry()));

    // AliAnaConvCorrPhoton * phJetAna = dynamic_cast<AliAnaConvCorrPhoton*>(fAnaPhotonArray->At(0));
    // phJetAna->PrintStatistics();

    return kTRUE;

}

///_______________________________________________________________________________
void AliAnalysisTaskGammaJet::NotifyRun() {
  //See header file for documentation
  
  AliInfo(Form("we have a new run: %d", fCurrentRunNumber));
  cout << Form("we have a new run: %d", fCurrentRunNumber) << endl;
  
}

///_______________________________________________________________________________
void AliAnalysisTaskGammaJet::GetPionGrandChildren(const AliAODConversionPhoton * pion, const TClonesArray * photons, Int_t* trackLabels) {
  ///Get the track labels of the electrons reconstructed as gamma forming the pion
  
  for(Int_t i = 0; i< 2; i++) {
    AliAODConversionPhoton * gamma = dynamic_cast<AliAODConversionPhoton*>(photons->At(pion->GetLabel(i)));
    
    if(gamma) { 
      for(Int_t j = 0; j< 2; j++) {
	trackLabels[ i*2+ j] = gamma->GetLabel(j);
      }
      
    } else {
      cout << "AliAnaConvCorrPion::GetTrackLabels() :: Not good!!!"<<endl;
    }
  }
}




//_____________________________________________________________________
void AliAnalysisTaskGammaJet::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
}

//_____________________________________________________________________
AliAODEvent * AliAnalysisTaskGammaJet::GetAODEvent() {
  //Get the AOD event from whereever it might be
  AliAODEvent * aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  

  if(!aodEvent) {
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!Getting AODEvent();" << endl;
    aodEvent = AODEvent();
  } else {
    //cout << "got aod from input event1" << endl;
  }
  return aodEvent;

}

//_____________________________________________________________________
TClonesArray * AliAnalysisTaskGammaJet::GetConversionGammas(const AliAODEvent * aodEvent) {

  //Get Conversion gamma branch of AOD. First try standard AOD
  TClonesArray * convGamma = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(Form("%s_gamma", fGammaCutString.Data())));
  
  //If it's there, send it back
  if(convGamma)  {
    //cout << "found conv gamma branch in aod event!!!"<<endl;
    return convGamma;
  }else {
    cout << "did NOT !!! find conv gamma branch in aod event!!!"<<endl;
  }

  //If AOD not in standard file have to locate it in delta AOD
  if( !(fDeltaAODFileName.Length() > 0)  ) return NULL;
  
  AliAODHandler * aodHandler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()); 
  if(aodHandler) {
    AliAODExtension * gExt = dynamic_cast<AliAODExtension*>(aodHandler->GetExtensions()->FindObject(fDeltaAODFileName));
    if(gExt) {
      AliAODEvent * gcEvent = gExt->GetAOD();
      return dynamic_cast<TClonesArray*>(gcEvent->FindListObject(Form("%s_gamma", fGammaCutString.Data())));
    }
  }

  cout << "could not find branch " << Form("%s_gamma", fPionCutString.Data()) << endl; 
  return NULL;
}



///_____________________________________________________________________
TClonesArray * AliAnalysisTaskGammaJet::GetConversionPions(const AliAODEvent * aodEvent) {

  //Get Conversion pion branch of AOD. First try standard AOD
  TClonesArray * convGamma = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(Form("%s_Pi0", fPionCutString.Data()) ));
  
  //If it's there, send it back
  if(convGamma)  return convGamma;


  //If AOD not in standard file have to locate it in delta AOD
  if( !(fDeltaAODFileName.Length() > 0)  ) return NULL;
  
  AliAODHandler * aodHandler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()); 
  if(aodHandler) {
    AliAODExtension * gExt = dynamic_cast<AliAODExtension*>(aodHandler->GetExtensions()->FindObject(fDeltaAODFileName));
    if(gExt) {
      AliAODEvent * gcEvent = gExt->GetAOD();
      return dynamic_cast<TClonesArray*>(gcEvent->FindListObject(Form("%s_Pi0", fPionCutString.Data())));
    }
  }  
  return NULL;
}
