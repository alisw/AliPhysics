class TParticle;

#include <Riostream.h>
#include "TFile.h"
#include "THistManager.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TParticlePDG.h"

#include "AliAnalysisTaskMultspec_MCpred.h"

ClassImp(AliAnalysisTaskMultspec_MCpred)

AliAnalysisTaskMultspec_MCpred::AliAnalysisTaskMultspec_MCpred() : AliAnalysisTaskSE(),
fHistos_misc(nullptr)
{
  //default constructor
}

AliAnalysisTaskMultspec_MCpred::AliAnalysisTaskMultspec_MCpred(const char *name, TString lExtraOptions) : AliAnalysisTaskSE(name),
fHistos_misc(nullptr)
{

  //Standard output
  DefineOutput(1, TList::Class()); // Miscellaneous Histograms


}

AliAnalysisTaskMultspec_MCpred::~AliAnalysisTaskMultspec_MCpred()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------
    if (fHistos_misc) {
        delete fHistos_misc;
        fHistos_misc = 0x0;
    }

}

//________________________________________________________________________
void AliAnalysisTaskMultspec_MCpred::UserCreateOutputObjects()
{

  //miscellaneous histograms
  fHistos_misc = new THistManager("fHistos_misc");
  fHistos_misc->CreateTH1("hNevts", "", 1, 0, 1);  // #events histo

  //histograms Gen_distrib VS mult
  fHistos_misc->CreateTH2("hmultspec_K0S_nch" , "", 100, 0, 100, 1000, 0, 1000, "s"); //vs nch
  fHistos_misc->CreateTH2("hmultspec_Lam_nch" , "", 100, 0, 100, 1000, 0, 1000, "s");
  fHistos_misc->CreateTH2("hmultspec_ALam_nch", "", 100, 0, 100, 1000, 0, 1000, "s");
  fHistos_misc->CreateTH2("hmultspec_Xi_nch"  , "", 100, 0, 100, 1000, 0, 1000, "s");
  fHistos_misc->CreateTH2("hmultspec_AXi_nch" , "", 100, 0, 100, 1000, 0, 1000, "s");
  fHistos_misc->CreateTH2("hmultspec_Om_nch"  , "", 100, 0, 100, 1000, 0, 1000, "s");
  fHistos_misc->CreateTH2("hmultspec_AOm_nch" , "", 100, 0, 100, 1000, 0, 1000, "s");
  fHistos_misc->CreateTH2("hmultspec_K0S_V0"  , "", 100, 0, 100, 1000, 0, 1000, "s"); //vs V0
  fHistos_misc->CreateTH2("hmultspec_Lam_V0"  , "", 100, 0, 100, 1000, 0, 1000, "s");
  fHistos_misc->CreateTH2("hmultspec_ALam_V0" , "", 100, 0, 100, 1000, 0, 1000, "s");
  fHistos_misc->CreateTH2("hmultspec_Xi_V0"   , "", 100, 0, 100, 1000, 0, 1000, "s");
  fHistos_misc->CreateTH2("hmultspec_AXi_V0"  , "", 100, 0, 100, 1000, 0, 1000, "s");
  fHistos_misc->CreateTH2("hmultspec_Om_V0"   , "", 100, 0, 100, 1000, 0, 1000, "s");
  fHistos_misc->CreateTH2("hmultspec_AOm_V0"  , "", 100, 0, 100, 1000, 0, 1000, "s");

  // multiplicity-related histos
  fHistos_misc->CreateTH2("hmult_nch_V0", "" , 1000, 0, 1000, 1000, 0, 1000, "s");

  //histograms multiplicity progression
  fHistos_misc->CreateTH1("hmultprogr_K0S_nch" , "", 1000, 0, 1000, "s"); //vs nch
  fHistos_misc->CreateTH1("hmultprogr_Lam_nch" , "", 1000, 0, 1000, "s");
  fHistos_misc->CreateTH1("hmultprogr_ALam_nch", "", 1000, 0, 1000, "s");
  fHistos_misc->CreateTH1("hmultprogr_Xi_nch"  , "", 1000, 0, 1000, "s");
  fHistos_misc->CreateTH1("hmultprogr_AXi_nch" , "", 1000, 0, 1000, "s");
  fHistos_misc->CreateTH1("hmultprogr_Om_nch"  , "", 1000, 0, 1000, "s");
  fHistos_misc->CreateTH1("hmultprogr_AOm_nch" , "", 1000, 0, 1000, "s");
  fHistos_misc->CreateTH1("hmultprogr_K0S_V0" , "", 1000, 0, 1000, "s"); //vs V0
  fHistos_misc->CreateTH1("hmultprogr_Lam_V0" , "", 1000, 0, 1000, "s");
  fHistos_misc->CreateTH1("hmultprogr_ALam_V0", "", 1000, 0, 1000, "s");
  fHistos_misc->CreateTH1("hmultprogr_Xi_V0"  , "", 1000, 0, 1000, "s");
  fHistos_misc->CreateTH1("hmultprogr_AXi_V0" , "", 1000, 0, 1000, "s");
  fHistos_misc->CreateTH1("hmultprogr_Om_V0"  , "", 1000, 0, 1000, "s");
  fHistos_misc->CreateTH1("hmultprogr_AOm_V0" , "", 1000, 0, 1000, "s");

  fHistos_misc->CreateTH1("hmultprogr_pi_V0"  , "", 1000, 0, 1000, "s"); //vs V0
  fHistos_misc->CreateTH1("hmultprogr_k_V0"  , "", 1000, 0, 1000, "s");
  fHistos_misc->CreateTH1("hmultprogr_p_V0" , "", 1000, 0, 1000, "s");

  //Output posting
  PostData(1, fHistos_misc->GetListOfHistograms());

}// end UserCreateOutputObjects

//________________________________________________________________________
void AliAnalysisTaskMultspec_MCpred::UserExec(Option_t *)
{

  //get MC event and fill histo
  AliMCEvent  *lMCev  = MCEvent();
  if (!lMCev) {
    Printf("ERROR: Could not retrieve MC event in file %s\n",fInputHandler->GetTree()->GetCurrentFile()->GetName());
    PostData(1, fHistos_misc->GetListOfHistograms());
    return;
  }

  //
  AliStack *lMCstack = lMCev->Stack();
  if (!lMCstack) {
    Printf("ERROR: Could not retrieve MC stack \n");
    return;
  }

  // take track of number of analyzed events
  fHistos_misc->FillTH1("hNevts", 0.5);

  int nch    = 0;
  int nV0    = 0;
  int n_K0S  = 0;
  int n_Lam  = 0;
  int n_ALam = 0;
  int n_Xi   = 0;
  int n_AXi  = 0;
  int n_Om   = 0;
  int n_AOm  = 0;

  int n_pi   = 0;
  int n_k    = 0;
  int n_p    = 0;

  for (int itrk = 0; itrk < lMCstack->GetNtrack(); itrk++){

    TParticle *lPart = lMCstack->Particle(itrk);
    if(!lPart || !lMCstack->IsPhysicalPrimary(itrk)) continue;
    double ch = ((TParticlePDG*)lPart->GetPDG())->Charge();

    if(!(TMath::Abs(ch)<1e-3)) {
      if(lPart->Eta()>-0.5 && lPart->Eta()<0.5) nch++;
      else if((lPart->Eta()>2.8 && lPart->Eta()<5.1)||(lPart->Eta()>-3.7 && lPart->Eta()<-1.7)) nV0++;
    }

    if(TMath::Abs(lPart->Y())>0.5) continue;

    int pdg = (int)lPart->GetPdgCode();
    if     (pdg==310  ) n_K0S++;
    else if(pdg==3122 ) n_Lam++;
    else if(pdg==-3122) n_ALam++;
    else if(pdg==3312 ) n_Xi++;
    else if(pdg==-3312) n_AXi++;
    else if(pdg==3334 ) n_Om++;
    else if(pdg==-3334) n_AOm++;
    else if(TMath::Abs(pdg)==211)  n_pi++;
    else if(TMath::Abs(pdg)==321)  n_k++;
    else if(TMath::Abs(pdg)==2212) n_p++;

  }

  fHistos_misc->FillTH2("hmult_nch_V0", nch, nV0);

  fHistos_misc->FillTH2("hmultspec_K0S_nch" , n_K0S, nch); //vs nch
  fHistos_misc->FillTH2("hmultspec_Lam_nch" , n_Lam, nch);
  fHistos_misc->FillTH2("hmultspec_ALam_nch", n_ALam,nch);
  fHistos_misc->FillTH2("hmultspec_Xi_nch"  , n_Xi,  nch);
  fHistos_misc->FillTH2("hmultspec_AXi_nch" , n_AXi, nch);
  fHistos_misc->FillTH2("hmultspec_Om_nch"  , n_Om,  nch);
  fHistos_misc->FillTH2("hmultspec_AOm_nch" , n_AOm, nch);
  fHistos_misc->FillTH2("hmultspec_K0S_V0"  , n_K0S, nV0); //vs V0
  fHistos_misc->FillTH2("hmultspec_Lam_V0"  , n_Lam, nV0);
  fHistos_misc->FillTH2("hmultspec_ALam_V0" , n_ALam,nV0);
  fHistos_misc->FillTH2("hmultspec_Xi_V0"   , n_Xi,  nV0);
  fHistos_misc->FillTH2("hmultspec_AXi_V0"  , n_AXi, nV0);
  fHistos_misc->FillTH2("hmultspec_Om_V0"   , n_Om,  nV0);
  fHistos_misc->FillTH2("hmultspec_AOm_V0"  , n_AOm, nV0);

  for(int i=0; i<n_K0S; i++) {
      fHistos_misc->FillTH1("hmultprogr_K0S_nch" ,nch);
      fHistos_misc->FillTH1("hmultprogr_K0S_V0"  ,nV0 );
  }
  for(int i=0; i<n_Lam; i++) {
      fHistos_misc->FillTH1("hmultprogr_Lam_nch" ,nch);
      fHistos_misc->FillTH1("hmultprogr_Lam_V0"  ,nV0 );
  }
  for(int i=0; i<n_ALam; i++) {
      fHistos_misc->FillTH1("hmultprogr_ALam_nch" ,nch);
      fHistos_misc->FillTH1("hmultprogr_ALam_V0"  ,nV0 );
  }
  for(int i=0; i<n_Xi; i++) {
      fHistos_misc->FillTH1("hmultprogr_Xi_nch" ,nch);
      fHistos_misc->FillTH1("hmultprogr_Xi_V0"  ,nV0 );
  }
  for(int i=0; i<n_AXi; i++) {
      fHistos_misc->FillTH1("hmultprogr_AXi_nch" ,nch);
      fHistos_misc->FillTH1("hmultprogr_AXi_V0"  ,nV0 );
  }
  for(int i=0; i<n_Om; i++) {
      fHistos_misc->FillTH1("hmultprogr_Om_nch" ,nch);
      fHistos_misc->FillTH1("hmultprogr_Om_V0"  ,nV0 );
  }
  for(int i=0; i<n_AOm; i++) {
      fHistos_misc->FillTH1("hmultprogr_AOm_nch" ,nch);
      fHistos_misc->FillTH1("hmultprogr_AOm_V0"  ,nV0 );
  }
  for(int i=0; i<n_pi; i++) {
      fHistos_misc->FillTH1("hmultprogr_pi_V0"  ,nV0 );
  }
  for(int i=0; i<n_k; i++) {
      fHistos_misc->FillTH1("hmultprogr_k_V0"  ,nV0 );
  }
  for(int i=0; i<n_p; i++) {
      fHistos_misc->FillTH1("hmultprogr_p_V0"  ,nV0 );
  }

  PostData(1, fHistos_misc->GetListOfHistograms());

}

//________________________________________________________________________
void AliAnalysisTaskMultspec_MCpred::Terminate(Option_t *)
{

}
