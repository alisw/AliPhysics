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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// This task is meant to acquire MC-level predictions for strangeness studies vs multiplicity
// and ZDC Energy.
//
// Modified version of AliAnalysisTaskMCPredictions.cxx from David Dobrigkeit Chinellato
//
// --- Francesca Ercolessi
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class TTree;
class TParticle;
class TVector3;

//class AliMCEventHandler;
//class AliMCEvent;
//class AliStack;

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;

#include <Riostream.h>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
//#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliPPVsMultUtils.h"
#include "AliPWG0Helper.h"
#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskMCPredictionsStrgVsMultVsZDC.h"
////
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"

#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliVertexingHFUtils.h"
//#include "AliPythia8.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMCPredictionsStrgVsMultVsZDC)

AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::AliAnalysisTaskMCPredictionsStrgVsMultVsZDC()
: AliAnalysisTaskSE(),
fkNSpecies(21),
fkSelectINELgtZERO(kTRUE),
fListHist(0),
fHistEventCounter(0),
fHistV0MMult(0),
fHistV0AMult(0),
fHistV0CMult(0),
fHistMult05(0),
fHistMult08(0),
fHistMult08to15(0),
fHistSPDClusters(0),
fHistNMPI(0),
fHistQ2(0),
fHistb(0),
fHistLeadingE(0),
fHistEffEnergy(0)
{
  for(Int_t ih=0; ih<21; ih++){
    fHistPt[ih] = 0x0;
    f2DHistPartSPDV0M[ih] = 0x0;
    f2DHistINELgt0SPDV0M[ih] = 0x0;
    f2DHistLeadingESPDV0M[ih] = 0x0;
    f2DHistEffEnergySPDV0M[ih] = 0x0;
    f2DHistNchSPDV0M[ih] = 0x0;
    f2DHistNMPISPDV0M[ih] = 0x0;
    f2DHistQ2SPDV0M[ih] = 0x0;
    f2DHistbSPDV0M[ih] = 0x0;
    f2DHistAvPtSPDV0M[ih] = 0x0;
  }
}

AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::AliAnalysisTaskMCPredictionsStrgVsMultVsZDC(const char *name, Float_t lCenterOfMassEnergy, Bool_t kDoPythia, Bool_t kDoEPOS)
: AliAnalysisTaskSE(name),
fCenterOfMassEnergy(lCenterOfMassEnergy),
fkNSpecies(21),
fkSelectINELgtZERO(kTRUE),
fkDoPythia(kDoPythia),
fkDoEPOS(kDoEPOS),
fListHist(0),
fHistEventCounter(0),
fHistV0MMult(0),
fHistV0AMult(0),
fHistV0CMult(0),
fHistMult05(0),
fHistMult08(0),
fHistMult08to15(0),
fHistSPDClusters(0),
fHistNMPI(0),
fHistQ2(0),
fHistb(0),
fHistLeadingE(0),
fHistEffEnergy(0)
{
  for(Int_t ih=0; ih<21; ih++){
    fHistPt[ih] = 0x0;
    f2DHistPartSPDV0M[ih] = 0x0;
    f2DHistINELgt0SPDV0M[ih] = 0x0;
    f2DHistLeadingESPDV0M[ih] = 0x0;
    f2DHistEffEnergySPDV0M[ih] = 0x0;
    f2DHistNchSPDV0M[ih] = 0x0;
    f2DHistNchSPDV0M[ih] = 0x0;
    f2DHistNMPISPDV0M[ih] = 0x0;
    f2DHistQ2SPDV0M[ih] = 0x0;
    f2DHistbSPDV0M[ih] = 0x0;
    f2DHistAvPtSPDV0M[ih] = 0x0;
  }
  //
  DefineOutput(1, TList::Class()); 
}


AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::~AliAnalysisTaskMCPredictionsStrgVsMultVsZDC()
{
  //------------------------------------------------
  // DESTRUCTOR
  //------------------------------------------------
  
  if (fListHist) {
    delete fListHist;
    fListHist = 0x0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::UserCreateOutputObjects()
{
  //------------------------------------------------
  // Histograms: Basic Analysis Output
  //------------------------------------------------
  // Create histograms
  fListHist = new TList();
  fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

  TString lPartNames[21] = {
    "PiPlus", "PiMinus", 
    "KaPlus", "KaMinus", 
    "Proton", "AntiProton",
    "K0Short", 
    "Lambda", "AntiLambda",
    "XiMinus", "XiPlus", 
    "OmegaMinus", "OmegaPlus",
    "Phi", 
    "D0", "AntiD0", 
    "DPlus", "DMinus", 
    "Lambdac", "AntiLambdac", 
    "JPsi"
  };
  
  //-----------------------------------------------------------------------------
  if(! fHistEventCounter ) {
    fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",2,0,2);
    fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
    fHistEventCounter->GetXaxis()->SetBinLabel(2, "INEL>0");
    fListHist->Add(fHistEventCounter);
  }

  //-----------------------------------------------------------------------------
  if(! fHistV0MMult ) {
    fHistV0MMult = new TH1D( "fHistV0MMult", ";V0M Mult;Count",1000,0,1000);
    fListHist->Add(fHistV0MMult);
  }

  //-----------------------------------------------------------------------------
  if(! fHistV0AMult ) {
    fHistV0AMult = new TH1D( "fHistV0AMult", ";V0A Mult;Count",1000,0,1000);
    fListHist->Add(fHistV0AMult);
  }

  //-----------------------------------------------------------------------------
  if(! fHistV0CMult ) {
    fHistV0CMult = new TH1D( "fHistV0CMult", ";V0C Mult;Count",1000,0,1000);
    fListHist->Add(fHistV0CMult);
  }

  //-----------------------------------------------------------------------------
  if(! fHistMult05 ) {
    fHistMult05 = new TH1D( "fHistMult05", ";Nch in |#eta|<0.5 ;Count",1000,0,1000);
    fListHist->Add(fHistMult05);
  }

  //-----------------------------------------------------------------------------
  if(! fHistMult08 ) {
    fHistMult08 = new TH1D( "fHistMult08", ";Nch in |#eta|<0.8 ;Count",1000,0,1000);
    fListHist->Add(fHistMult08);
  }

  //-----------------------------------------------------------------------------
  if(! fHistMult08to15 ) {
    fHistMult08to15 = new TH1D( "fHistMult08to15", ";Nch in 0.8<|#eta|<1.5 ;Count",1000,0,1000);
    fListHist->Add(fHistMult08to15);
  }

  //-----------------------------------------------------------------------------
  if(! fHistSPDClusters ) {
    fHistSPDClusters = new TH1D( "fHistSPDClusters", ";SPD Clusters ;Count",1000,0,1000);
    fListHist->Add(fHistSPDClusters);
  }

  //-----------------------------------------------------------------------------
  if(!fHistNMPI && fkDoPythia) {
    fHistNMPI = new TH1D( "fHistNMPI", ";NMPI;Count",50,0,50);
    fListHist->Add(fHistNMPI);
  }

  //-----------------------------------------------------------------------------
  if(!fHistQ2 && fkDoPythia) {
    fHistQ2 = new TH1D( "fHistQ2", "; Q^{2} ;Count",1000,0,1000);
    fListHist->Add(fHistQ2);
  }

  //-----------------------------------------------------------------------------
  if(!fHistb && fkDoEPOS) {
    fHistb = new TH1D( "fHistb", ";b;Count",20,0,20);
    fListHist->Add(fHistb);
  }

  //-----------------------------------------------------------------------------
  if(!fHistLeadingE ) {
    fHistLeadingE = new TH1D( "fHistLeadingE", ";Leading energy (GeV);Count",1300,0,fCenterOfMassEnergy);
    fListHist->Add(fHistLeadingE);
  }

  //-----------------------------------------------------------------------------
  if(!fHistEffEnergy ) {
    fHistEffEnergy = new TH1D( "fHistEffEnergy", ";Effective energy (GeV);Count",1300,0,fCenterOfMassEnergy);
    fListHist->Add(fHistEffEnergy);
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(!fHistPt[ih]) {
      fHistPt[ih] = new TH1D(Form("fHistPt_%s",lPartNames[ih].Data()), Form("Generated %s;p_{T} (GeV/c)",lPartNames[ih].Data()), 250, 0, 25.);
      fListHist->Add(fHistPt[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(!f2DHistPartSPDV0M[ih]) {
      f2DHistPartSPDV0M[ih] = new TH2D(Form("f2DHistPartSPDV0M_%s",lPartNames[ih].Data()), Form("Generated %s;SPD Clusters; V0M multiplicity",lPartNames[ih].Data()), 800, 0, 800., 500, 0, 500.);
      fListHist->Add(f2DHistPartSPDV0M[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(!f2DHistINELgt0SPDV0M[ih]) {
      f2DHistINELgt0SPDV0M[ih] = new TH2D(Form("f2DHistINELgt0SPDV0M_%s",lPartNames[ih].Data()), "INEL>0 events;SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
      fListHist->Add(f2DHistINELgt0SPDV0M[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(!f2DHistLeadingESPDV0M[ih]) {
      f2DHistLeadingESPDV0M[ih] = new TH2D(Form("f2DHistLeadingESPDV0M_%s",lPartNames[ih].Data()), "Leading energy (GeV);SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
      fListHist->Add(f2DHistLeadingESPDV0M[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(!f2DHistEffEnergySPDV0M[ih]) {
      f2DHistEffEnergySPDV0M[ih] = new TH2D(Form("f2DHistEffEnergySPDV0M_%s",lPartNames[ih].Data()), "Effective energy (GeV);SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
      fListHist->Add(f2DHistEffEnergySPDV0M[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(!f2DHistNchSPDV0M[ih]) {
      f2DHistNchSPDV0M[ih] = new TH2D(Form("f2DHistNchSPDV0M_%s",lPartNames[ih].Data()), "Nch (|#eta|<0.5);SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
      fListHist->Add(f2DHistNchSPDV0M[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(!f2DHistNMPISPDV0M[ih] && fkDoPythia) {
      f2DHistNMPISPDV0M[ih] = new TH2D(Form("f2DHistNMPISPDV0M_%s",lPartNames[ih].Data()), "NMPI;SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
      fListHist->Add(f2DHistNMPISPDV0M[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(!f2DHistQ2SPDV0M[ih] && fkDoPythia) {
      f2DHistQ2SPDV0M[ih] = new TH2D(Form("f2DHistQ2SPDV0M_%s",lPartNames[ih].Data()), "Q2;SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
      fListHist->Add(f2DHistQ2SPDV0M[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(!f2DHistbSPDV0M[ih] && fkDoEPOS) {
      f2DHistbSPDV0M[ih] = new TH2D(Form("f2DHistbSPDV0M_%s",lPartNames[ih].Data()), "b;SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
      fListHist->Add(f2DHistbSPDV0M[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  for(Int_t ih=0; ih<fkNSpecies; ih++){
    if(!f2DHistAvPtSPDV0M[ih]) {
      f2DHistAvPtSPDV0M[ih] = new TH2D(Form("f2DHistAvPtSPDV0M_%s",lPartNames[ih].Data()), "#LT Pt #GT (GeV/c);SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
      fListHist->Add(f2DHistAvPtSPDV0M[ih]);
    }
  }

  //List of Histograms: Normal
  PostData(1, fListHist);
  
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::UserExec(Option_t *)
{
  // Main loop --> called for each event  
  AliMCEvent  *lMCevent  = 0x0;
  AliStack    *lMCstack  = 0x0;
  //
  lMCevent = MCEvent();
  if (!lMCevent) {
    Printf("ERROR: Could not retrieve MC event \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }
  //
  lMCstack = lMCevent->Stack();
  if (!lMCstack) {
    Printf("ERROR: Could not retrieve MC stack \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }
  //
  AliGenEventHeader* mcGenH = lMCevent->GenEventHeader();
  
  //Events Processed
  fHistEventCounter->Fill(0.5);

  //Multiplicity
  Long_t lNchEta05         = 0;
  Long_t lNchEta08         = 0;
  Long_t lNchEta08to15     = 0;
  Long_t lSPDCl0           = 0;
  Long_t lSPDCl1           = 0;
  Long_t lNchVZEROA        = 0;
  Long_t lNchVZEROC        = 0;
  Bool_t lEvSel_INELgtZERO = kFALSE;

  //Eff energy
  Float_t fLeadingE = 0.;
  Float_t fEffEnergy = fCenterOfMassEnergy;

  //Utility
  const Float_t c = 299792458; //m/s
  Bool_t lIsPhysicalPrimary = kFALSE;

  //Particle info
  Int_t lPartCounter[21];
  for(int i = 0; i<21; i++){
    lPartCounter[i] = 0;
  }

  TString lPartNames[21] = {
    "PiPlus", "PiMinus", 
    "KaPlus", "KaMinus", 
    "Proton", "AntiProton",
    "K0Short", 
    "Lambda", "AntiLambda",
    "XiMinus", "XiPlus", 
    "OmegaMinus", "OmegaPlus",
    "Phi", 
    "D0", "AntiD0", 
    "DPlus", "DMinus", 
    "Lambdac", "AntiLambdac", 
    "JPsi"
  };
  Int_t lPDGCodes[21] = {
    211, -211, 
    321, -321, 
    2212, -2212,
    310, 
    3122, -3122,
    3312, -3312, 
    3334, -3334,
    333, 
    421, -421, 
    411, -411, 
    4122, -4122, 
    443    
  };  

  Bool_t lCheckIsPhysicalPrimary[21] = {
    kTRUE, kTRUE, 
    kTRUE, kTRUE, 
    kTRUE, kTRUE,
    kTRUE, 
    kTRUE, kTRUE,
    kTRUE, kTRUE, 
    kTRUE, kTRUE,
    kFALSE, 
    kFALSE, kFALSE,
    kFALSE, kFALSE, 
    kFALSE, kFALSE, 
    kFALSE
  };
    
  //----- Loop on Stack ----------------------------------------------------------------
  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
  { // This is the begining of the loop on tracks
    
    TParticle* particleOne = lMCstack->Particle(iCurrentLabelStack);
    if(!particleOne) continue;
    if(!particleOne->GetPDG()) continue;

    Double_t geta = particleOne -> Eta();
    Int_t charge = particleOne->GetPDG()->Charge();
    Double_t partcharge = particleOne->GetPDG()->Charge()/3.;

    //Vtx position of the particle
    Double_t vtxpart = particleOne->R();

    //Vtx position of daughters
    Int_t idau = particleOne->GetFirstDaughter();
    TParticle* firstdau;
    Double_t vtxdau = 9999;
    if(idau>=0 && idau<lMCstack->GetNtrack()) {
      firstdau = lMCstack->Particle(idau);
      vtxdau = firstdau->R();
    }
    
    if( TMath::Abs(partcharge)>0.001 ) {
      if(vtxpart <= 7. && vtxdau > 7.) {
        if( (TMath::Abs(geta) < 1.5928278) ) lSPDCl1++; 
      }
      if(vtxpart <= 4. && vtxdau > 4.) {
        if( (TMath::Abs(geta) < 2.1245920) ) lSPDCl0++; 
      }
    }

    if(! (lMCstack->IsPhysicalPrimary(iCurrentLabelStack)) ) continue;
    if(particleOne->GetFirstDaughter()>0) continue;
    
    if( (charge>0 && TMath::Abs(geta)>7. && TMath::Abs(geta)<8.7) || 
        (charge==0 && TMath::Abs(geta)>8.8)
      ) {
        fEffEnergy -= particleOne -> Energy();
        fLeadingE += particleOne -> Energy();
      } 
    
    if(TMath::Abs(partcharge)<0.001) continue; //now only charged primaries

    //
    if( TMath::Abs(geta) < 0.5 ) lNchEta05++;
    if( TMath::Abs(geta) < 0.8 ) lNchEta08++;
    if( (TMath::Abs(geta) > 0.8) && (TMath::Abs(geta) < 1.5) ) lNchEta08to15++;   
    if( TMath::Abs(geta) < 1.0 ) lEvSel_INELgtZERO = kTRUE;
    if( 2.8 < geta && geta < 5.1 ) lNchVZEROA++;
    if(-3.7 < geta && geta <-1.7 ) lNchVZEROC++;

  }//End of loop on tracks
  //----- End Loop on Stack ------------------------------------------------------------

  //Reject non-INEL>0 if requested
  if( !lEvSel_INELgtZERO && fkSelectINELgtZERO ) return;

  //
  Long_t lSPDClusters = lSPDCl0 + lSPDCl1;

  //Events INEL>0
  fHistEventCounter->Fill(1.5);
  
  //
  Int_t fMC_NMPI = -1;
  Float_t fMC_Q2 = -1;
  Float_t fMC_b = -1;
  //Int_t fMC_NPart = -1;
  //Int_t fMC_NColl = -1;
  
  //PYTHIA
  if ( fkDoPythia ){
    if (mcGenH->InheritsFrom(AliGenPythiaEventHeader::Class())){
      AliGenPythiaEventHeader *fMcPythiaHeader = dynamic_cast <AliGenPythiaEventHeader*> (mcGenH);
      if(fMcPythiaHeader){
        fMC_NMPI = fMcPythiaHeader->GetNMPI();
        fMC_Q2 = fMcPythiaHeader->GetPtHard();
      }
    }
  }
  
  //EPOS 
  if ( fkDoEPOS ){
    if (mcGenH->InheritsFrom(AliGenHepMCEventHeader::Class())) {
      AliGenHepMCEventHeader * lHepMCHeader = dynamic_cast <AliGenHepMCEventHeader*> (mcGenH);    
      if (lHepMCHeader ){
        //fMC_NPart = lHepMCHeader->Npart_proj()+lHepMCHeader->Npart_targ();
        //fMC_NColl = lHepMCHeader->N_Nwounded_collisions() + lHepMCHeader->Nwounded_N_collisions() + lHepMCHeader->Nwounded_Nwounded_collisions();      
        fMC_b = lHepMCHeader->impact_parameter();
      }
    }
  }

  //1D histograms
  if(fHistV0MMult)      fHistV0MMult        -> Fill ( lNchVZEROA+lNchVZEROC );
  if(fHistV0AMult)      fHistV0AMult        -> Fill ( lNchVZEROA );
  if(fHistV0CMult)      fHistV0CMult        -> Fill ( lNchVZEROC );
  if(fHistMult05)       fHistMult05         -> Fill ( lNchEta05 );
  if(fHistMult08)       fHistMult08         -> Fill ( lNchEta08 );
  if(fHistMult08to15)   fHistMult08to15     -> Fill ( lNchEta08to15 );
  if(fHistSPDClusters)  fHistSPDClusters    -> Fill ( lSPDClusters );
  if(fHistNMPI)         fHistNMPI           -> Fill ( fMC_NMPI );
  if(fHistQ2)           fHistQ2             -> Fill ( fMC_Q2 );
  if(fHistb)            fHistb              -> Fill ( fMC_b );
  if(fHistLeadingE)     fHistLeadingE       -> Fill ( fLeadingE );
  if(fHistEffEnergy)    fHistEffEnergy      -> Fill ( fEffEnergy );

  //----- Loop on Stack ----------------------------------------------------------------
  for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
  { // This is the begining of the loop on tracks
    
    TParticle* part = lMCstack->Particle(iCurrentLabelStack);
    if(!part) continue;
    if(!part->GetPDG()) continue;
    lIsPhysicalPrimary = lMCstack->IsPhysicalPrimary(iCurrentLabelStack);
    
    Double_t charge = part -> GetPDG()->Charge();
    Float_t px      = part -> Px();
    Float_t py      = part -> Py();
    Float_t pz      = part -> Pz();
    Float_t pt      = part -> Pt();
    Float_t energy  = part -> Energy();
    Float_t eta     = part -> Eta();
    Int_t pdg       = (Int_t) part -> GetPdgCode();
    Float_t y       = TMath::Log( TMath::Sqrt( ( energy + pz*c )/( energy - pz*c ) ) );
        
    if(TMath::Abs(eta)<0.5){   
      for(Int_t ih=0; ih<21; ih++){ //loop over pdg codes
        if( pdg == lPDGCodes[ih] ) {
          //
          if( lCheckIsPhysicalPrimary[ih] == kTRUE && lIsPhysicalPrimary == kFALSE ) continue;
    
          //Fill histos
          if(fHistPt[ih]) fHistPt[ih] -> Fill( pt );
          
          //Counter for specific parricles in this event
          lPartCounter[ih]++;
        }
      }//end loop over pdg codes        
    }  
  } //end loop on tracks 

  for(Int_t ih=0; ih<21; ih++){ //loop over pdg codes
    if(f2DHistPartSPDV0M[ih])      f2DHistPartSPDV0M[ih]      -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, lPartCounter[ih]);    
    if(f2DHistINELgt0SPDV0M[ih])   f2DHistINELgt0SPDV0M[ih]   -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC );    
    if(f2DHistLeadingESPDV0M[ih])  f2DHistLeadingESPDV0M[ih]  -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fLeadingE );    
    if(f2DHistEffEnergySPDV0M[ih]) f2DHistEffEnergySPDV0M[ih] -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fEffEnergy );    
    if(f2DHistNchSPDV0M[ih])       f2DHistNchSPDV0M[ih]       -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, lNchEta05 );    
    if(f2DHistNMPISPDV0M[ih])      f2DHistNMPISPDV0M[ih]      -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fMC_NMPI );   
    if(f2DHistQ2SPDV0M[ih])        f2DHistQ2SPDV0M[ih]        -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fMC_Q2 );
    if(f2DHistbSPDV0M[ih])         f2DHistbSPDV0M[ih]         -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fMC_b );
    if(f2DHistAvPtSPDV0M[ih])      f2DHistAvPtSPDV0M[ih]      -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fHistPt[ih]->GetMean() );
  }

  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  // Post output data.
  PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  
  TList *cRetrievedList = 0x0;
  cRetrievedList = (TList*)GetOutputData(1);
  if(!cRetrievedList) {
    Printf("ERROR - AliAnalysisTaskMCPredictionsStrgVsMultVsZDC : ouput data container list not available\n");
    return;
  }
  
  fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
  if (!fHistEventCounter) {
    Printf("ERROR - AliAnalysisTaskMCPredictionsStrgVsMultVsZDC : fHistEventCounter not available");
    return;
  }
  
  TCanvas *canCheck = new TCanvas("AliAnalysisTaskMCPredictionsStrgVsMultVsZDC","Event Multiplicity",10,10,510,510);
  canCheck->cd(1)->SetLogy();
  
  fHistEventCounter->SetMarkerStyle(22);
  fHistEventCounter->DrawCopy("E");
}

//______________________________________________________________________
Bool_t AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::IsEPOSLHC() const {
  //Function to check if this is DPMJet
  Bool_t lReturnValue = kFALSE;
  AliMCEvent*  mcEvent = MCEvent();
  if (mcEvent) {
    AliGenEventHeader* mcGenH = mcEvent->GenEventHeader();
    //A bit uncivilized, but hey, if it works...
    TString lHeaderTitle = mcGenH->GetName();
    if (lHeaderTitle.Contains("EPOSLHC")) {
      //This header has "EPOS" in its title!
      lReturnValue = kTRUE;
    }
  }
  return lReturnValue;
}
