#include <TH1F.h>
#include <TH2F.h>
#include <TRandom3.h>
#include <TMath.h>

#include "AliParticleContainer.h"
#include "AliAnalysisTaskBackFlucRandomCone.h"
#include "AliRhoParameter.h"

ClassImp(AliAnalysisTaskBackFlucRandomCone)

//________________________________________________________________________________________________
AliAnalysisTaskBackFlucRandomCone::AliAnalysisTaskBackFlucRandomCone() : AliAnalysisTaskEmcalJet("AliAnalysisTaskBackFlucRandomCone", kTRUE),
fEtaMin(-0.5),
fEtaMax(0.5),
fEta(0),
fPhi(0),
fRCone(0.3),
fRnd(0),
fRhoName(""),
fhEtaPhiConeCentre(0),
fNInOut(0),
fhMasspTInCone(0),
fhNConstituents(0),
fhDeltapT(0)
{
   /// default constructor
   SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________________________________
AliAnalysisTaskBackFlucRandomCone::AliAnalysisTaskBackFlucRandomCone(const char *name) : AliAnalysisTaskEmcalJet(name, kTRUE),
fEtaMin(-0.5),
fEtaMax(0.5),
fEta(0),
fPhi(0),
fRCone(0.3),
fRnd(0),
fRhoName(""),
fhEtaPhiConeCentre(0),
fNInOut(0),
fhMasspTInCone(0),
fhNConstituents(0),
fhDeltapT(0)

{
   /// standard constructor
   SetMakeGeneralHistograms(kTRUE);
}
//________________________________________________________________________________________________

void AliAnalysisTaskBackFlucRandomCone::UserCreateOutputObjects(){
   /// Create output
   AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
 
   //init TRandom
   fRnd = TRandom3();
   
   Int_t nBinsM = 100;
   Double_t minM = 0, maxM = 20;
   Int_t nBinspT = 120;
   Double_t minpT = 0, maxpT = 119;
   
   fhEtaPhiConeCentre = new TH2F("fhEtaPhiConeCentre", "Cone centre; #eta ; #varphi", 50, fEtaMin, fEtaMax, 50, 0., TMath::TwoPi());
   fOutput->Add(fhEtaPhiConeCentre);
   
   fNInOut = new TH1F("fNInOut", "In = 1 out = 2 of cone", 2, 0.5, 2.5);
   fOutput->Add(fNInOut);
   
   fhMasspTInCone = new TH2F("fhMasspTInCone", "In Cone mass and pT distribution; M ; #it{p}_{T};", nBinsM, minM, maxM, nBinspT, minpT, maxpT);
   fOutput->Add(fhMasspTInCone);
   
   fhNConstituents = new TH1F("fhNConstituents", "Constituents of the random cone jet; #it{N}_{constituents}", 21, 0., 20.);
   fOutput->Add(fhNConstituents);
   
   fhDeltapT = new TH1F("fhDeltapT", "#Delta #it{p}_{T}; #it{p}_{T} - A#rho", nBinspT, minpT-20, maxpT-20);
   fOutput->Add(fhDeltapT);
   
}
//________________________________________________________________________________________________
Bool_t AliAnalysisTaskBackFlucRandomCone::Run(){
   /// Run code before FillHistograms()
   
   fPhi = fRnd.Uniform(0, TMath::TwoPi());
   fEta = fRnd.Uniform(fEtaMin,fEtaMax);
   
   return kTRUE;
}
//________________________________________________________________________________________________
Bool_t AliAnalysisTaskBackFlucRandomCone::FillHistograms(){
   /// Fill histograms defined
   
   fhEtaPhiConeCentre->Fill(fEta, fPhi);
   
   //loop on tracks in the event
   AliParticleContainer *trackCont = GetParticleContainer(0);
   Double_t pTinCone = 0, massinCone = 0, mass = 139.57;
   Int_t Nconst = 0;
   for(Int_t itr=0; itr < trackCont->GetNAcceptedParticles(); itr++){
      AliVParticle *part = trackCont->GetParticle(itr);
      if ((TMath::Abs(part->Eta() - fEta) > fRCone) || (TMath::Abs(part->Phi() - fPhi) > fRCone)) {
      	 fNInOut->Fill(2);
      	 continue;
      }
      Nconst++;
      fNInOut->Fill(1);
      Double_t ptpart = part->Pt();
      pTinCone += ptpart;
      //massinCone += (mass*mass + ptpart*ptpart );
      
   }
   
   fhMasspTInCone->Fill(0.,pTinCone);
   fhNConstituents->Fill(Nconst);
   Double_t rhovalue = 0;
   AliRhoParameter *rhopar = GetRhoFromEvent(fRhoName); 
   if(rhopar) rhovalue = rhopar->GetVal();
   fhDeltapT->Fill(pTinCone - TMath::Pi()*fRCone*fRCone*rhovalue);
   
   return kTRUE;
}


//________________________________________________________________________________________________

AliAnalysisTaskBackFlucRandomCone::~AliAnalysisTaskBackFlucRandomCone() {
   /// Empty destructor.
}

//________________________________________________________________________________________________

void AliAnalysisTaskBackFlucRandomCone::Terminate(Option_t *) 
{
   /// Called once at the end of the analysis.
}
