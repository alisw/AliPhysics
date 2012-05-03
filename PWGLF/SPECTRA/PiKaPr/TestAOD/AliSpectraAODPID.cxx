#include "AliSpectraAODPID.h"
#include "AliAODEvent.h"      
#include "TH1F.h"             
#include "TH2F.h"             
#include "TList.h"            
#include "AliAODTrack.h"      
#include "AliAODMCParticle.h" 
#include "AliPIDResponse.h"   
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliSpectraAODTrackCuts.h"

ClassImp(AliSpectraAODPID)

AliSpectraAODPID::AliSpectraAODPID() : TObject(), fPIDType(kNSigmaTPCTOF), fNSigmaPID(3), fPIDResponse(0) {

}

AliSpectraAODPID::AliSpectraAODPID(AODPIDType_t pidType) : TObject(), fPIDType(pidType), fNSigmaPID(3), fPIDResponse(0) {



}



void AliSpectraAODPID::FillQAHistos(AliSpectraAODHistoManager * hman, AliAODTrack * track, AliSpectraAODTrackCuts * trackCuts) {

  // fill a bunch of QA histos


  // Get PID response object, if needed
  if(!fPIDResponse) {
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
  }
  
  //Response
  AliVParticle *inEvHMain = dynamic_cast<AliVParticle *>(track);

  hman->GetPIDHistogram(kHistPIDTPC)->Fill(track->GetTPCmomentum(), track->GetTPCsignal()*track->Charge()); // PID histo


  Double_t nsigmaTPCkProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton));
  Double_t nsigmaTPCkKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon)); 
  Double_t nsigmaTPCkPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion)); 
  Double_t nsigmaTOFkProton=0,nsigmaTOFkKaon=0,nsigmaTOFkPion=0;

  if(track->Pt()>trackCuts->GetPtTOFMatching()){
    hman->GetPIDHistogram(kHistPIDTOF)->Fill(track->P(),(track->GetTOFsignal()/100)*track->Charge()); // PID histo
    nsigmaTOFkProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kProton));
    nsigmaTOFkKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kKaon)); 
    nsigmaTOFkPion = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kPion)); 
    
    //TOF
    hman->GetPtHistogram(kHistNSigProtonTOF)->Fill(track->P(),fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kProton));
    hman->GetPtHistogram(kHistNSigKaonTOF)->Fill(track->P(),fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kKaon));
    hman->GetPtHistogram(kHistNSigPionTOF)->Fill(track->P(),fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kPion));
    hman->GetPtHistogram(kHistNSigProtonPtTOF)->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kProton));
    hman->GetPtHistogram(kHistNSigKaonPtTOF)->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kKaon));
    hman->GetPtHistogram(kHistNSigPionPtTOF)->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kPion));
  }
  
  Double_t nsigmaTPCTOFkProton = TMath::Sqrt(nsigmaTPCkProton*nsigmaTPCkProton+nsigmaTOFkProton*nsigmaTOFkProton);
  Double_t nsigmaTPCTOFkKaon = TMath::Sqrt(nsigmaTPCkKaon*nsigmaTPCkKaon+nsigmaTOFkKaon*nsigmaTOFkKaon);
  Double_t nsigmaTPCTOFkPion = TMath::Sqrt(nsigmaTPCkPion*nsigmaTPCkPion+nsigmaTOFkPion*nsigmaTOFkPion);
  
  //TPC
  hman->GetPtHistogram(kHistNSigProtonTPC)->Fill(track->P(),fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton));
  hman->GetPtHistogram(kHistNSigKaonTPC)->Fill(track->P(),fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon));
  hman->GetPtHistogram(kHistNSigPionTPC)->Fill(track->P(),fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion));
  hman->GetPtHistogram(kHistNSigProtonPtTPC)->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton));
  hman->GetPtHistogram(kHistNSigKaonPtTPC)->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon));
  hman->GetPtHistogram(kHistNSigPionPtTPC)->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion));
  //TPCTOF
  hman->GetPtHistogram(kHistNSigProtonTPCTOF)->Fill(track->P(),nsigmaTPCTOFkProton);
  hman->GetPtHistogram(kHistNSigKaonTPCTOF)->Fill(track->P(),nsigmaTPCTOFkKaon);
  hman->GetPtHistogram(kHistNSigPionTPCTOF)->Fill(track->P(),nsigmaTPCTOFkPion);
  hman->GetPtHistogram(kHistNSigProtonPtTPCTOF)->Fill(track->Pt(),nsigmaTPCTOFkProton);
  hman->GetPtHistogram(kHistNSigKaonPtTPCTOF)->Fill(track->Pt(),nsigmaTPCTOFkKaon);
  hman->GetPtHistogram(kHistNSigPionPtTPCTOF)->Fill(track->Pt(),nsigmaTPCTOFkPion);

}

Int_t AliSpectraAODPID::GetParticleSpecie(AliAODMCParticle * part) {
  // return PID according to MC truth
  switch(TMath::Abs(part->PdgCode())){
  case 2212:
    return kSpProton;
    break;
  case 321:
    return kSpKaon;
    break;
  case 211:
    return kSpPion;
    break;
  default:
    return kSpUndefined;
  } 
}


Int_t AliSpectraAODPID::GetParticleSpecie(AliAODTrack      * trk, AliSpectraAODTrackCuts * trackCuts) {
  // return PID according to detectors
  
  // Get PID response object, if needed
  if(!fPIDResponse) {
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
  }

  if(!fPIDResponse) {
    AliFatal("Cannot get pid response");
    return 0;
  }


  // Compute nsigma for each hypthesis
  AliVParticle *inEvHMain = dynamic_cast<AliVParticle *>(trk);

  // --- TPC
  Double_t nsigmaTPCkProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton));
  Double_t nsigmaTPCkKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon)); 
  Double_t nsigmaTPCkPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion)); 
  // --- TOF
  Double_t nsigmaTOFkProton=0,nsigmaTOFkKaon=0,nsigmaTOFkPion=0;
  if(trk->Pt()>trackCuts->GetPtTOFMatching()){
    nsigmaTOFkProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kProton));
    nsigmaTOFkKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kKaon)); 
    nsigmaTOFkPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kPion)); 
  }
	  
  // --- combined
  Double_t nsigmaTPCTOFkProton = TMath::Sqrt(nsigmaTPCkProton*nsigmaTPCkProton+nsigmaTOFkProton*nsigmaTOFkProton);
  Double_t nsigmaTPCTOFkKaon   = TMath::Sqrt(nsigmaTPCkKaon*nsigmaTPCkKaon+nsigmaTOFkKaon*nsigmaTOFkKaon);
  Double_t nsigmaTPCTOFkPion   = TMath::Sqrt(nsigmaTPCkPion*nsigmaTPCkPion+nsigmaTOFkPion*nsigmaTOFkPion);


  // select the nsigma to be used for the actual PID
  Double_t nsigmaPion, nsigmaKaon, nsigmaProton;

  switch (fPIDType) {
  case kNSigmaTPC:
    nsigmaProton  =  nsigmaTPCkProton;
    nsigmaKaon	  =  nsigmaTPCkKaon  ;
    nsigmaPion    =  nsigmaTPCkPion  ;
    break;
  case kNSigmaTOF:
    nsigmaProton  =  nsigmaTOFkProton;
    nsigmaKaon	  =  nsigmaTOFkKaon  ;
    nsigmaPion    =  nsigmaTOFkPion  ;
    break;
  case kNSigmaTPCTOF:
    nsigmaProton  =  nsigmaTPCTOFkProton;
    nsigmaKaon	  =  nsigmaTPCTOFkKaon  ;
    nsigmaPion    =  nsigmaTPCTOFkPion  ;
    break;
  }

  // guess the particle based on the smaller nsigma
  if( ( nsigmaKaon   < nsigmaPion ) && ( nsigmaKaon < nsigmaProton ) && (nsigmaKaon   < fNSigmaPID)) return kSpKaon;
  if( ( nsigmaPion   < nsigmaKaon ) && ( nsigmaPion < nsigmaProton ) && (nsigmaPion   < fNSigmaPID)) return kSpPion;
  if( ( nsigmaProton < nsigmaKaon ) && ( nsigmaProton < nsigmaPion ) && (nsigmaProton < fNSigmaPID)) return kSpProton;

  // else, return undefined
  return kSpUndefined;

}
