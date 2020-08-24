#include "AliAnalysisTaskEffMatrix.h"
#include "TList.h"
#include "TClonesArray.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliVParticle.h"
#include "AliAODMCParticle.h"
#include "AliESDtrackCuts.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TString.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCHeader.h"

#include <iostream>

/*
 * Efficiency Matrix Analysis
 * author: Roberto Preghenella (R+)
 * email:  preghenella@bo.infn.it
 *
 */

ClassImp(AliAnalysisTaskEffMatrix)

//________________________________________________________________________


//________________________________________________________________________

AliAnalysisTaskEffMatrix::AliAnalysisTaskEffMatrix(const Char_t *name) :
AliAnalysisTaskSE(name),
  fAODfilterBit(AliAODTrack::kTrkGlobal),
  fESDtrackCuts(NULL),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fPtMin(0.15),
  fPtMax(1.e6),
  fRapidityMin(-0.5),
  fRapidityMax(0.5),
  fMCArray(NULL),
  fOutputList(NULL)
{
    
  /*
   * default constructor
   */

  TDatabasePDG *dbpdg = TDatabasePDG::Instance();
  dbpdg->AddParticle("Lambda1520", "Lambda1520", 1.5195, 0, 0.0156, 0, "Unknown", 3124, -3124);
  dbpdg->AddParticle("Lambda1520_bar", "Lambda1520_bar", 1.5195, 0, 0.0156, 0, "Unknown", -3124, 3124);
    
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________

AliAnalysisTaskEffMatrix::~AliAnalysisTaskEffMatrix()
{
    
  /*
   * default destructor
   */
    
  if (fOutputList) delete fOutputList;
    
}

//________________________________________________________________________

void AliAnalysisTaskEffMatrix::UserCreateOutputObjects()
{
    
  /*
   * user create output objects
   */
    
  Int_t fgDataBins[kNData] = {
    20,
    /*32*/16, 
    /*16*/7, 
    /*180*/18, 
    /*20.*/1, 
    10
  };
  Double_t fgDataMin[kNData] = {
    TMath::Log(/*0.1*/0.15), 
    -0.8, 
    -0.7,
    0.,
    -10., 
    0.
  };
  Double_t fgDataMax[kNData] = {
    TMath::Log(10.), 
    0.8, 
    0.7, 
    2. * TMath::Pi(), 
    10., 
    100.
  };
    
  const Char_t *fgStageName[kNStages] = {
    "generated", "reconstructed", "matched"
  };

  const Int_t fgParticleCode[] = {
    211, -211, /* pions */
    321, -321, /* kaons */
    2212, -2212, /* protons */
    113, /* rho */
    313, -313, /* K*0 */
    333, /* phi */
    3124, -3124 /* lambda* */
  };

  TDatabasePDG *dbpdg = TDatabasePDG::Instance();
    
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
    
  for (Int_t istage = 0; istage < kNStages; istage++) {

    fHistoList[istage] = new TList();
    fHistoList[istage]->SetOwner(kTRUE);
    fHistoList[istage]->SetName(fgStageName[istage]);

    for (Int_t ipart = 0; ipart < sizeof(fgParticleCode) / 4; ipart++) {
	
      TParticlePDG *ppdg = dbpdg->GetParticle(fgParticleCode[ipart]);
      const Char_t *pname = ppdg ? ppdg->GetName() : "unknown";
      THnSparseF *histo = new THnSparseF(Form("hHistoData_%s", pname), "", kNData, fgDataBins, fgDataMin, fgDataMax);
      fHistoList[istage]->Add(histo);

    }

    fOutputList->Add(fHistoList[istage]);
      
  }
    
  PostData(1, fOutputList);
}

//________________________________________________________________________

void AliAnalysisTaskEffMatrix::UserExec(Option_t *)
{

  /*
   * user exec
   */
  
  /* get event */
  AliVEvent *event = InputEvent();
  if (!event) return;
  fMCArray = (TClonesArray *)event->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  AliAODMCHeader* mcHeader  = (AliAODMCHeader*) event->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  
  /* accept event */
  if (!AcceptEvent(event)) return;

  /* useful stuff */
  TObjArray reconstructedArray, matchedArray;
  Double_t data[kNData];
  THnSparseF *histo = NULL;
  TDatabasePDG *dbpdg = TDatabasePDG::Instance();
  TString hname;
    
  /* centrality percentile */
  //  AliCentrality *centrality = event->GetCentrality();
  //  data[kData_centrality] = centrality->GetCentralityPercentile("V0M");

 

  AliMultSelection *mult_selection = (AliMultSelection*)event->GetList()->FindObject("MultSelection");
  if (mult_selection == NULL) {
    return;
  }

  data[kData_centrality] = mult_selection->GetMultiplicityPercentile("V0M");

  
  /* vertex z */
  const AliVVertex *vertex = GetPrimaryVertex(event);
  data[kData_zv] = vertex->GetZ();
 
  /* loop over reconstructed tracks */
  for (Int_t itrk = 0; itrk < event->GetNumberOfTracks(); itrk++) {
    
    /* get track */
    AliVTrack *track = (AliVTrack *)event->GetTrack(itrk);
    if (!track) continue;
    /* accept track */
    if (!AcceptTrack(track)) continue;
    
    /* get particle */
    Int_t label = TMath::Abs(track->GetLabel());
    AliVParticle *particle = (AliVParticle *)fMCArray->At(label);
    if (!particle) continue;
    
    /* check TPC PID */
    if (!HasTPCPID(track)) continue;
    /* add to reconstructed array */
    reconstructedArray.Add(particle);
    
    /* check TOF PID */
    if (!HasTOFPID(track)) continue;
    /* add to matched array */
    matchedArray.Add(particle);
    
  } /* end of loop over reconstructed tracks */
  
    /* loop over generated particles */
  for (Int_t imc = 0; imc < fMCArray->GetEntries(); imc++) {
    
    /* get particle */
    AliVParticle *particle = (AliVParticle *)fMCArray->At(imc);
    if (!particle) continue;

    // check if particle is pileup particle
    
    if ( AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(imc, mcHeader, fMCArray) ) {
           // std::cout << " --- Found pile-up particle, the centrality is: " << mult_selection->GetMultiplicityPercentile("V0M") << " with PDG code: "  << 	(dbpdg->GetParticle(particle->PdgCode()))->GetName() << " Accept Particle: " << AcceptParticle(particle) << std::endl;
      continue;
    }

    
    data[kData_logpt] = particle->Pt() > 0. ? TMath::Log(particle->Pt()) : -100.;
    data[kData_eta] = particle->Eta();
    data[kData_y] = particle->Y();
    data[kData_phi] = particle->Phi();
    Int_t pdg = particle->PdgCode();
    TParticlePDG *ppdg = dbpdg->GetParticle(pdg);
    const Char_t *pname = ppdg ? ppdg->GetName() : "unknown";
    hname = Form("hHistoData_%s", pname);
    
    /* switch PDG code */
    switch (pdg) {
	  
    case 211: case -211: /* pions */
    case 321: case -321: /* kaons */
    case 2212: case -2212: /* protons */
	  
      /* accept particle */
      if (!AcceptParticle(particle)) continue;	  
      /* fill generated stage */
      histo = (THnSparseF *)fHistoList[kStage_generated]->FindObject(hname.Data());
      if (!histo) continue;
      histo->Fill(data);
	  
      /* check if reconstructed */
      if (!reconstructedArray.Contains(particle)) continue;
      /* fill reconstructed stage */
      histo = (THnSparseF *)fHistoList[kStage_reconstructed]->FindObject(hname.Data());
      if (!histo) continue;
      histo->Fill(data);
	  
      /* check if matched */
      if (!matchedArray.Contains(particle)) continue;
      /* fill matched stage */
      histo = (THnSparseF *)fHistoList[kStage_matched]->FindObject(hname.Data());
      if (!histo) continue;
      histo->Fill(data);
	  
      break;
      /* end of pions */
      /* end of kaons */
      /* end of protons */

    case 113: /* rho */
    case 313: case -313: /* K*0 */
    case 333: /* phi */
    case 3124: case -3124: /* lambda* */
      
      Int_t dpdg1 = 0;
      Int_t dpdg2 = 0;
      if (pdg == 113) {dpdg1 = 211; dpdg2 = -211;} /* rho */
      if (pdg == 313) {dpdg1 = 321; dpdg2 = -211;} /* K*0 */
      if (pdg == -313) {dpdg1 = -321; dpdg2 = 211;} /* K*0_bar */
      if (pdg == 333) {dpdg1 = 321; dpdg2 = -321;} /* phi */
      if (pdg == 3124) {dpdg1 = 2212; dpdg2 = -321;} /* lambda* */
      if (pdg == -3124) {dpdg1 = -2212; dpdg2 = 321;} /* lambda* */

      /* accept resonance */
      if (!AcceptResonance(particle, dpdg1, dpdg2)) continue;

#if 0
      /* that's debug */
      if (pdg == 3124 || pdg == -3124) {

	/* compute distance from the primary vertex */
	Double_t pVertexXYZ[3];
	vertex->GetXYZ(pVertexXYZ);
	Double_t rVertexXYZ[3];
	particle->XvYvZv(rVertexXYZ);
	Double_t dVertex = TMath::Sqrt((rVertexXYZ[0]-pVertexXYZ[0])*(rVertexXYZ[0]-pVertexXYZ[0]) +
				       (rVertexXYZ[1]-pVertexXYZ[0])*(rVertexXYZ[1]-pVertexXYZ[1]) +
				       (rVertexXYZ[2]-pVertexXYZ[0])*(rVertexXYZ[2]-pVertexXYZ[2]));
				       
	printf(">>> L = %f \n", dVertex);
	particle->Print();
	Int_t imoth = particle->GetMother();
	if (imoth > 0) {
	  AliVParticle *particle = (AliVParticle *)fMCArray->At(imoth);
	  particle->Print();
	}	
      
#endif     
 
      /* fill generated stage */
      histo = (THnSparseF *)fHistoList[kStage_generated]->FindObject(hname.Data());
      if (!histo) continue;
      histo->Fill(data);
      
      /* get daughters */
      Int_t dl1 = particle->GetDaughterFirst();
      Int_t dl2 = particle->GetDaughterLast();
      AliVParticle *daughter1 = (AliVParticle *)fMCArray->At(dl1);
      AliVParticle *daughter2 = (AliVParticle *)fMCArray->At(dl2);
      
      /* check if reconstructed */
      if (!reconstructedArray.Contains(daughter1) ||
	  !reconstructedArray.Contains(daughter2)) continue;
      /* fill reconstructed stage */
      histo = (THnSparseF *)fHistoList[kStage_reconstructed]->FindObject(hname.Data());
      if (!histo) continue;
      histo->Fill(data);
      
      /* check if matched */
      if (!matchedArray.Contains(daughter1) ||
	  !matchedArray.Contains(daughter2)) continue;
      /* fill matched stage */
      histo = (THnSparseF *)fHistoList[kStage_matched]->FindObject(hname.Data());
      if (!histo) continue;
      histo->Fill(data);
      
      break;
      /* end of phi */

#if 0
      /* all the rest */
    default:
      break;
#endif
      
    } /* end of switch PDG code */
  } /* end of loop over generated particles */
  
  PostData(1, fOutputList);
}

//___________________________________________________________

const AliVVertex *
AliAnalysisTaskEffMatrix::GetPrimaryVertex(AliVEvent *event) const
{
  /*
   * get primary vertex
   */

  /* avoid ESD analysis until fully implemented */
  if (event->InheritsFrom("AliESDEvent")) return NULL;
        
  /* check physics selection */
  //    if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB)) return kFALSE;
    
  /* check primary vertex */
  const AliVVertex *vertex = event->GetPrimaryVertex();
  if (vertex->GetNContributors() < 1) {
    /* get AOD vertex SPD */
    if (event->InheritsFrom("AliAODEvent")) {
      AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
      if (!aodevent) return NULL;
      vertex = aodevent->GetPrimaryVertexSPD();
    }
    /* get ESD vertex SPD */
    else if (event->InheritsFrom("AliESDEvent")) {
      AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
      if (!esdevent) return NULL;
      vertex = esdevent->GetPrimaryVertexSPD();
    }
    else return NULL;
  }
  if (vertex->GetNContributors() < 1) return NULL;

  return vertex;
}

//___________________________________________________________

Bool_t
AliAnalysisTaskEffMatrix::AcceptEvent(AliVEvent *event) const
{
  /*
   * accept event
   */
    
  /* avoid ESD analysis until fully implemented */
  if (event->InheritsFrom("AliESDEvent")) return kFALSE;
        
  /* check physics selection */
  //    if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB)) return kFALSE;
    
  /* check primary vertex */
  const AliVVertex *vertex = GetPrimaryVertex(event);
  if (!vertex) return kFALSE;

  /* check vertex position */
  if (TMath::Abs(vertex->GetZ()) > 10.) return kFALSE;
    
  /* event accepted */
  return kTRUE;
}

//___________________________________________________________

Bool_t
AliAnalysisTaskEffMatrix::AcceptTrack(AliVTrack *track) const
{
  /*
   * accept track
   */
    
  /* check AOD filter bit */
  if (track->InheritsFrom("AliAODTrack")) {
    AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(track);
    if (!aodtrack) return kFALSE;
    if (!aodtrack->TestFilterBit(fAODfilterBit)) return kFALSE;
  }
  /* check ESD track cuts */
  else if (track->InheritsFrom("AliESDtrack")) {
    AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(track);
    if (!esdtrack) return kFALSE;
    if (!fESDtrackCuts->AcceptTrack(esdtrack)) return kFALSE;
  }
  else return kFALSE;
    
  /* check eta range */
  if (track->Eta() < fEtaMin ||
      track->Eta() > fEtaMax) return kFALSE;
  /* check pt range */
  if (track->Pt() < fPtMin ||
      track->Pt() > fPtMax) return kFALSE;
    
  /* track accepted */
  return kTRUE;
}

//___________________________________________________________

Bool_t
AliAnalysisTaskEffMatrix::AcceptParticle(AliVParticle *particle) const
{
  /*
   * accept particle
   */
    
  /* check AOD MC particle */
  if (particle->InheritsFrom("AliAODMCParticle")) {
    AliAODMCParticle *aodparticle = dynamic_cast<AliAODMCParticle *>(particle);
    if (!aodparticle) return kFALSE;
    if (!aodparticle->IsPhysicalPrimary()) return kFALSE;
  }
  else return kFALSE;
  
  /* particle accepted */
  return kTRUE;

}


//___________________________________________________________

Bool_t
AliAnalysisTaskEffMatrix::AcceptResonance(AliVParticle *particle, Int_t dpdg1, Int_t dpdg2) const
{
  /*
   * accept resonance
   */
    
  /* check AOD MC particle */
  if (particle->InheritsFrom("AliAODMCParticle")) {
    AliAODMCParticle *aodparticle = dynamic_cast<AliAODMCParticle *>(particle);
    if (!aodparticle) return kFALSE;
    //        if (!aodparticle->IsPhysicalPrimary()) return kFALSE;
  }
  else return kFALSE;

  /* check rapidity */
  if (particle->Y() < fRapidityMin || particle->Y() > fRapidityMax)
    return kFALSE;

  /* check daughters */
  Int_t dl1 = particle->GetDaughterFirst();
  Int_t dl2 = particle->GetDaughterLast();

  if (dl1 < 0 || dl2 < 0) {
    //    particle->Print();
    return kFALSE;
  }
  
  AliVParticle *daughter1 = (AliVParticle *)fMCArray->At(dl1);
  AliVParticle *daughter2 = (AliVParticle *)fMCArray->At(dl2);

  if ((daughter1->PdgCode() != dpdg1 || daughter2->PdgCode() != dpdg2) &&
      (daughter1->PdgCode() != dpdg2 || daughter2->PdgCode() != dpdg1)) return kFALSE;

  /* resonance accepted */
  return kTRUE;

}


//___________________________________________________________

Bool_t
AliAnalysisTaskEffMatrix::HasTPCPID(AliVTrack *track) const
{
  /*
   * has TPC PID
   */
    
  /* check PID signal */
  if (track->GetTPCsignal() <= 0. ||
      track->GetTPCsignalN() == 0) return kFALSE;
  return kTRUE;
    
}

//___________________________________________________________

Bool_t
AliAnalysisTaskEffMatrix::HasTOFPID(AliVTrack *track) const
{
  /*
   * has TOF PID
   */
    
  /* check TOF matched track */
  if (!(track->GetStatus() & AliESDtrack::kTOFout)||
      !(track->GetStatus() & AliESDtrack::kTIME)) return kFALSE;
  return kTRUE;
    
}
