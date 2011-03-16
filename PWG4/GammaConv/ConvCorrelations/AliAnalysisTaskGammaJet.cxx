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
#include "AliGammaConversionAODObject.h"
#include "AliAODConversionParticle.h"
#include "AliAODJet.h"

#include "AliAODInputHandler.h"

#include "AliAnaConvIsolation.h"

#include "AliAnaConvCorrPhoton.h"
#include "AliAnaConvCorrPhotonJet.h"
#include "AliAnaConvCorrPion.h"
// Gamma - jet correlation analysis task
// Authors: Svein Lindal


using namespace std;

ClassImp(AliAnalysisTaskGammaJet)

//________________________________________________________________________
AliAnalysisTaskGammaJet::AliAnalysisTaskGammaJet() : 
AliAnalysisTaskSE(), 
  fOutputList(NULL), 
  fDeltaAODFileName("AliAODGammaConversion.root"),
  fConversionCutString("GammaConv"),
  fAnaIsolation(NULL), 
  fAnaIsolationArray(NULL), 
  fAnaPionArray(NULL), 
  fAnaPhotonArray(NULL),
  fAnaPhotonJetArray(NULL),
  fMinPt(3.0), 
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
  fDeltaAODFileName("AliAODGammaConversion.root"),
  fConversionCutString("GammaConv"),
  fAnaIsolation(NULL),
  fAnaIsolationArray(NULL),
  fAnaPionArray(NULL),
  fAnaPhotonArray(NULL),
  fAnaPhotonJetArray(NULL), 
  fMinPt(3.0),
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

  fAnaIsolationArray = new TObjArray();
  fAnaIsolationArray->SetOwner(kTRUE);

}

//________________________________________________________________________
void AliAnalysisTaskGammaJet::UserCreateOutputObjects()
{
  //Create histograms add, to outputlist
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

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
  

  PostData(1, fOutputList);


}

//________________________________________________________________________
void AliAnalysisTaskGammaJet::UserExec(Option_t *) 
{
  //Inherited from AliAnalysisTaskSE

  ///Get AOD event
  AliAODEvent * aodEvent = GetAODEvent();
  if(!aodEvent) {
    AliError("No AOD event!!");
    return;
  }

  TClonesArray * photons = GetConversionGammas(aodEvent);
  if(!photons) {
    AliError("No Conversion gamma!!");
    return;
  }

  TClonesArray * pions = GetConversionPions(aodEvent);
  if(!pions) {
    AliError("No Conversion gamma!!");
    return;
  }

  TClonesArray * tracks = aodEvent->GetTracks();
  if(!tracks) {
    AliError("Can't get tracks!!");
    return;
  }

  if(!((Entry() )%10000)) {
    AliInfo(Form("%s ----> Processing event # %lld", CurrentFileName(), Entry()));
    AliInfo(Form("%d %d", photons->GetEntriesFast(), pions->GetEntriesFast()));
  }


 if(aodEvent->GetNumberOfTracks() < fMinNTracks) return;

  if(photons->GetEntriesFast() > aodEvent->GetNumberOfTracks()) {
    AliError(Form("more conv gamma than tracks, ntracks %d, nconvGamma %d:  ", aodEvent->GetNumberOfTracks(), photons->GetEntriesFast()));
    return;

  } else if(photons->GetEntriesFast() > 0) {
    
    ProcessConvGamma(photons, pions, tracks);
    ProcessPions(pions, photons, tracks);
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

//___________________________________________________________________________________________
void AliAnalysisTaskGammaJet::ProcessConvGamma( const TClonesArray * convGamma, const TClonesArray * const pions, const TClonesArray * tracks ) {
  //See header file for documentation
  

  for (Int_t iPhot = 0; iPhot < convGamma->GetEntriesFast(); iPhot++) {
    Bool_t delP = kFALSE;
    AliAODConversionParticle * photon = dynamic_cast<AliAODConversionParticle*>(convGamma->At(iPhot));
    if(!photon) {
      AliGammaConversionAODObject * aodO = dynamic_cast<AliGammaConversionAODObject*>(convGamma->At(iPhot));
      if (!aodO) {
	AliError(Form("ERROR: Could not receive ga %d\n", iPhot));
	continue;
      }
      
      photon = new AliAODConversionParticle(aodO);
      delP = kTRUE;
    } 
    
    if(photon->Pt() < fMinPt) {
      if(delP) delete photon;
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
    if (delP) delete photon;
  } // 
}

///______________________________________________________________________________________________
void AliAnalysisTaskGammaJet::ProcessPions( const TClonesArray * const pions, const TClonesArray * const photons, const TClonesArray * const tracks) {
  //See header file for documentation

  for (Int_t iPi = 0; iPi < pions->GetEntriesFast(); iPi++) {
    Bool_t delP = kFALSE;
    AliAODConversionParticle * pion = dynamic_cast<AliAODConversionParticle*>(pions->At(iPi));
    if(!pion) {
      AliGammaConversionAODObject * aodO = dynamic_cast<AliGammaConversionAODObject*>(pions->At(iPi));
      if (!aodO) {
	AliError(Form("ERROR: Could not receive ga %d\n", iPi));
	continue;
      }
      
      pion = new AliAODConversionParticle(aodO);
      delP = kTRUE;
    } 


    if(pion->Pt() < fMinPt) {
      if(delP) delete pion;
      continue;
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
    } 
    if (delP) delete pion;
  } // 


}

///_______________________________________________________________________________________________
Bool_t AliAnalysisTaskGammaJet::IsDecayPion(Int_t iPhot, const TClonesArray * const pions) {
  //See header file for documentation

  for(int ip = 0; ip < pions->GetEntriesFast(); ip++) {
    AliAODConversionParticle * pion = dynamic_cast<AliAODConversionParticle*>(pions->At(ip));
    if(pion) {
      if(pion->GetLabel1() == iPhot || pion->GetLabel2() == iPhot) 
	
	return kTRUE;
    } else {
      AliGammaConversionAODObject * aodPion = dynamic_cast<AliGammaConversionAODObject*>(pions->At(ip));
      if(aodPion) {
	if(aodPion->GetLabel1() == iPhot || aodPion->GetLabel2() == iPhot) 
	  return kTRUE;
      } else {
	AliError("pion corrupted!");
      }
    }
  }

  return kFALSE;
}


///_______________________________________________________________________________
Bool_t AliAnalysisTaskGammaJet::UserNotify() {
    //See header file for documentation

    AliInfo(Form("%s ----> Processing event # %lld", CurrentFileName(), Entry()));

    AliAnaConvCorrPhoton * phJetAna = dynamic_cast<AliAnaConvCorrPhoton*>(fAnaPhotonArray->At(0));
    phJetAna->PrintStatistics();

    return kTRUE;

}

///_______________________________________________________________________________
void AliAnalysisTaskGammaJet::NotifyRun() {
  //See header file for documentation
  
  AliInfo(Form("we have a new run: %d", fCurrentRunNumber));
  cout << Form("we have a new run: %d", fCurrentRunNumber) << endl;
  
}

///_______________________________________________________________________________
void AliAnalysisTaskGammaJet::GetPionGrandChildren(const AliAODConversionParticle * pion, const TClonesArray * photons, Int_t* trackLabels) {
  ///Get the track labels of the electrons reconstructed as gamma forming the pion

  for(Int_t i = 0; i< 2; i++) {
    AliAODConversionParticle * gamma = dynamic_cast<AliAODConversionParticle*>(photons->At(pion->GetTrackLabel(i)));

    if(gamma) { 
      for(Int_t j = 0; j< 2; j++) {
	trackLabels[ i*2+ j] = gamma->GetTrackLabel(j);
      }

    } else {
      AliGammaConversionAODObject * aodO = dynamic_cast<AliGammaConversionAODObject*>(photons->At(pion->GetTrackLabel(i)));
      if(aodO) {
	trackLabels[i*2] = aodO->GetLabel1();
	trackLabels[i*2 + 1] = aodO->GetLabel2();
      } else {
	cout << "AliAnaConvCorrPion::GetTrackLabels() :: Not good!!!"<<endl;
      }
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
    aodEvent = AODEvent();
  }
  
  return aodEvent;

}

//_____________________________________________________________________
TClonesArray * AliAnalysisTaskGammaJet::GetConversionGammas(const AliAODEvent * aodEvent) {

  //Get Conversion gamma branch of AOD. First try standard AOD
  TClonesArray * convGamma = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(Form("%s_gamma", fConversionCutString.Data())));
  
  //If it's there, send it back
  if(convGamma)  return convGamma;


  //If AOD not in standard file have to locate it in delta AOD
  if( !(fDeltaAODFileName.Length() > 0)  ) return NULL;
  
  AliAODHandler * aodHandler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()); 
  if(aodHandler) {
    AliAODExtension * gExt = dynamic_cast<AliAODExtension*>(aodHandler->GetExtensions()->FindObject(fDeltaAODFileName));
    if(gExt) {
      AliAODEvent * gcEvent = gExt->GetAOD();
      return dynamic_cast<TClonesArray*>(gcEvent->FindListObject(Form("%s_gamma", fConversionCutString.Data())));
    }
  }

  cout << "could not find branch " << Form("%s_gamma", fConversionCutString.Data()) << endl; 
  return NULL;
}



///_____________________________________________________________________
TClonesArray * AliAnalysisTaskGammaJet::GetConversionPions(const AliAODEvent * aodEvent) {

  //Get Conversion pion branch of AOD. First try standard AOD
  TClonesArray * convGamma = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(Form("%s_Pi0", fConversionCutString.Data()) ));
  
  //If it's there, send it back
  if(convGamma)  return convGamma;


  //If AOD not in standard file have to locate it in delta AOD
  if( !(fDeltaAODFileName.Length() > 0)  ) return NULL;
  
  AliAODHandler * aodHandler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()); 
  if(aodHandler) {
    AliAODExtension * gExt = dynamic_cast<AliAODExtension*>(aodHandler->GetExtensions()->FindObject(fDeltaAODFileName));
    if(gExt) {
      AliAODEvent * gcEvent = gExt->GetAOD();
      return dynamic_cast<TClonesArray*>(gcEvent->FindListObject(Form("%s_Pi0", fConversionCutString.Data())));
    }
  }  
  return NULL;
}
