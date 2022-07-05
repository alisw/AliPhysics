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
fkNSpecies(22),
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
fHistEffEnergy(0),
f2DHistINELgt0SPDV0M(0),
f2DHistLeadingESPDV0M(0),
f2DHistEffEnergySPDV0M(0),
f2DHistNchSPDV0M(0),
f2DHistNMPISPDV0M(0),
f2DHistQ2SPDV0M(0),
f2DHistbSPDV0M(0),
f2dHistZDCVsLE(0),
f2dHistZDCVsEE(0),
f2dHistZDCVsLEA(0),
f2dHistZDCVsLEC(0),
f2dHistZPVsLP(0),
f2dHistZNVsLN(0),
f2dHistSPDClRecoVsTrue(0),
f2dHistV0MRecoVsTrue(0)
{
  for(Int_t ih=0; ih<22; ih++){
    fHistPt[ih] = 0x0;
    f2DHistPartSPDV0M[ih] = 0x0;
    f2DHistAvPtSPDV0M[ih] = 0x0;
  }
}

AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::AliAnalysisTaskMCPredictionsStrgVsMultVsZDC(const char *name, Float_t lCenterOfMassEnergy, Bool_t kDoPythia, Bool_t kDoEPOS)
: AliAnalysisTaskSE(name),
fCenterOfMassEnergy(lCenterOfMassEnergy),
fkNSpecies(22),
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
fHistEffEnergy(0),
f2DHistINELgt0SPDV0M(0),
f2DHistLeadingESPDV0M(0),
f2DHistEffEnergySPDV0M(0),
f2DHistNchSPDV0M(0),
f2DHistNMPISPDV0M(0),
f2DHistQ2SPDV0M(0),
f2DHistbSPDV0M(0),
f2dHistZDCVsLE(0),
f2dHistZDCVsEE(0),
f2dHistZDCVsLEA(0),
f2dHistZDCVsLEC(0),
f2dHistZPVsLP(0),
f2dHistZNVsLN(0),
f2dHistSPDClRecoVsTrue(0),
f2dHistV0MRecoVsTrue(0)
{
  for(Int_t ih=0; ih<22; ih++){
    fHistPt[ih] = 0x0;
    f2DHistPartSPDV0M[ih] = 0x0;
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

  TString lPartNames[22] = {
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
    "JPsi",
    "Pi0"
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
  if(!f2DHistINELgt0SPDV0M) {
    f2DHistINELgt0SPDV0M = new TH2D( "f2DHistINELgt0SPDV0M", "INEL>0 events;SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
    fListHist->Add(f2DHistINELgt0SPDV0M);
  }
  

  //-----------------------------------------------------------------------------
  if(!f2DHistLeadingESPDV0M) {
    f2DHistLeadingESPDV0M = new TH2D( "f2DHistLeadingESPDV0M", "Leading energy (GeV);SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
    fListHist->Add(f2DHistLeadingESPDV0M);
  }
  

  //-----------------------------------------------------------------------------
  if(!f2DHistEffEnergySPDV0M) {
    f2DHistEffEnergySPDV0M = new TH2D( "f2DHistEffEnergySPDV0M", "Effective energy (GeV);SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
    fListHist->Add(f2DHistEffEnergySPDV0M);
  }
  

  //-----------------------------------------------------------------------------
  if(!f2DHistNchSPDV0M) {
    f2DHistNchSPDV0M = new TH2D( "f2DHistNchSPDV0M", "Nch (|#eta|<0.5);SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
    fListHist->Add(f2DHistNchSPDV0M);
  }
  

  //-----------------------------------------------------------------------------
  if(!f2DHistNMPISPDV0M && fkDoPythia) {
    f2DHistNMPISPDV0M = new TH2D( "f2DHistNMPISPDV0M", "NMPI;SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
    fListHist->Add(f2DHistNMPISPDV0M);
  }
  

  //-----------------------------------------------------------------------------
  if(!f2DHistQ2SPDV0M && fkDoPythia) {
    f2DHistQ2SPDV0M = new TH2D( "f2DHistQ2SPDV0M", "Q2;SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
    fListHist->Add(f2DHistQ2SPDV0M);
  }
  

  //-----------------------------------------------------------------------------
  if(!f2DHistbSPDV0M && fkDoEPOS) {
    f2DHistbSPDV0M = new TH2D( "f2DHistbSPDV0M", "b;SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
    fListHist->Add(f2DHistbSPDV0M);
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
    if(!f2DHistAvPtSPDV0M[ih]) {
      f2DHistAvPtSPDV0M[ih] = new TH2D(Form("f2DHistAvPtSPDV0M_%s",lPartNames[ih].Data()), "#LT Pt #GT (GeV/c);SPD Clusters; V0M multiplicity", 800, 0, 800., 500, 0, 500.);
      fListHist->Add(f2DHistAvPtSPDV0M[ih]);
    }
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZDCVsLE) {
    f2dHistZDCVsLE = new TH2D("f2dHistZDCVsLE", ";ZDC Energy Sum (a.u.);Leading energy (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZDCVsLE);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZDCVsEE) {
    f2dHistZDCVsEE = new TH2D("f2dHistZDCVsEE", ";ZDC Energy Sum (a.u.); Effective energy (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZDCVsEE);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZDCVsLEA) {
    f2dHistZDCVsLEA = new TH2D("f2dHistZDCVsLEA", ";ZDC-A Energy Sum (a.u.);Leading energy A-side (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZDCVsLEA);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZDCVsLEC) {
    f2dHistZDCVsLEC = new TH2D("f2dHistZDCVsLEC", ";ZDC-C Energy Sum (a.u.);Leading energy C-side (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZDCVsLEC);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZPVsLP) {
    f2dHistZPVsLP = new TH2D("f2dHistZPVsLP", ";ZP Energy Sum (a.u.);Leading proton energy (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZPVsLP);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistZNVsLN) {
    f2dHistZNVsLN = new TH2D("f2dHistZNVsLN", ";ZN Energy Sum (a.u.);Leading neutron C-side (GeV);",400,0.,4000.,1300, 0., fCenterOfMassEnergy);
    fListHist->Add(f2dHistZNVsLN);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistSPDClRecoVsTrue) {
    f2dHistSPDClRecoVsTrue = new TH2D("f2dHistSPDClRecoVsTrue", "; SPDClusters centrality (%); SPD Clusters (true)", 100,0.,100., 800, 0, 800.);
    fListHist->Add(f2dHistSPDClRecoVsTrue);
  }

  //-----------------------------------------------------------------------------
  if(!f2dHistV0MRecoVsTrue) {
    f2dHistV0MRecoVsTrue = new TH2D("f2dHistV0MRecoVsTrue", ";V0M centrality (%); V0M Multiplicity (true)",100,0.,100., 500, 0, 500.);
    fListHist->Add(f2dHistV0MRecoVsTrue);
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
  Float_t fLeadingE_Aside = 0.;
  Float_t fLeadingE_Cside = 0.;
  Float_t fLeadingP = 0.;
  Float_t fLeadingN = 0.;
  Float_t fEffEnergy = fCenterOfMassEnergy;

  //Utility
  const Float_t c = 299792458; //m/s
  Bool_t lIsPhysicalPrimary = kFALSE;

  //Particle info
  Int_t lPartCounter[22];
  for(int i = 0; i<22; i++){
    lPartCounter[i] = 0;
  }

  TString lPartNames[22] = {
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
    "JPsi",
    "Pi0"
  };
  Int_t lPDGCodes[22] = {
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
    443,
    111   
  };  

  Bool_t lCheckIsPhysicalPrimary[22] = {
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
    kFALSE,
    kTRUE
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

        if( (charge>0 && TMath::Abs(geta)>7. && TMath::Abs(geta)<8.7) ) fLeadingP += particleOne -> Energy();
        if (charge==0 && TMath::Abs(geta)>8.8) fLeadingN += particleOne -> Energy();

        if((charge>0 && geta>7. && geta<8.7) || (charge==0 && geta>8.8)) {
          fLeadingE_Aside += particleOne -> Energy();
        } else if((charge>0 && geta>-8.7 && geta<-7) || (charge==0 && geta<-8.8)) {
          fLeadingE_Cside += particleOne -> Energy();
        }
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
    Float_t y       = Rapidity(energy, pz);
        
    if(TMath::Abs(y)<0.5){   
      for(Int_t ih=0; ih<22; ih++){ //loop over pdg codes
        if( pdg == lPDGCodes[ih] ) {
          //
          if( lCheckIsPhysicalPrimary[ih] == kTRUE && lIsPhysicalPrimary == kFALSE ) continue;
    
          //Fill histos
          if(fHistPt[ih]) fHistPt[ih] -> Fill( pt );
          if(f2DHistAvPtSPDV0M[ih]) f2DHistAvPtSPDV0M[ih] -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, pt );
          
          //Counter for specific parricles in this event
          lPartCounter[ih]++;
        }
      }//end loop over pdg codes        
    }  
  } //end loop on tracks 

  if(f2DHistINELgt0SPDV0M)   f2DHistINELgt0SPDV0M   -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC );    
  if(f2DHistLeadingESPDV0M)  f2DHistLeadingESPDV0M  -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fLeadingE );    
  if(f2DHistEffEnergySPDV0M) f2DHistEffEnergySPDV0M -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fEffEnergy );    
  if(f2DHistNchSPDV0M)       f2DHistNchSPDV0M       -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, lNchEta05 );    
  if(f2DHistNMPISPDV0M)      f2DHistNMPISPDV0M      -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fMC_NMPI );   
  if(f2DHistQ2SPDV0M)        f2DHistQ2SPDV0M        -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fMC_Q2 );
  if(f2DHistbSPDV0M)         f2DHistbSPDV0M         -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, fMC_b );
    
  for(Int_t ih=0; ih<22; ih++){ //loop over pdg codes
    if(f2DHistPartSPDV0M[ih])   f2DHistPartSPDV0M[ih]  -> Fill( lSPDClusters, lNchVZEROA+lNchVZEROC, lPartCounter[ih]);    
  }


  //Reco information

  AliESDEvent *lESDevent = 0x0;    
  lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!lESDevent) {
      AliWarning("ERROR: lESDevent not available \n");
      return;
  }

  Float_t lV0MPercentile = 300;
  Float_t lSPDClusterspercentile = 300;

  AliMultSelection *MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
  if( !MultSelection) {
    //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
    return;
  } else {
      lV0MPercentile = MultSelection->GetMultiplicityPercentile("V0M");
      lSPDClusterspercentile = MultSelection->GetMultiplicityPercentile("SPDClusters");
  }
    
  fCentrality_V0M = lV0MPercentile;
  fCentrality_SPDClusters = lSPDClusterspercentile; 
    
  // ZDC info ==========================================================
  const Double_t *aZDCN1 = lESDevent->GetESDZDC()->GetZNCTowerEnergy();
  fZNCpp = aZDCN1[0];
  const Double_t *aZDCN2 = lESDevent->GetESDZDC()->GetZNATowerEnergy();
  fZNApp = aZDCN2[0];
  const Double_t *aZDCP1 = lESDevent->GetESDZDC()->GetZPCTowerEnergy();
  fZPCpp = aZDCP1[0];  
  const Double_t *aZDCP2 = lESDevent->GetESDZDC()->GetZPATowerEnergy();
  fZPApp = aZDCP2[0];  

  if(f2dHistZDCVsLE)    f2dHistZDCVsLE         -> Fill ( fZNCpp+fZNApp+fZPCpp+fZPApp , fLeadingE );
  if(f2dHistZDCVsEE)    f2dHistZDCVsEE         -> Fill ( fZNCpp+fZNApp+fZPCpp+fZPApp , fEffEnergy );
  if(f2dHistZDCVsLEA)   f2dHistZDCVsLEA        -> Fill ( fZNApp+fZPApp , fLeadingE_Aside );
  if(f2dHistZDCVsLEC)   f2dHistZDCVsLEC        -> Fill ( fZNCpp+fZPCpp , fLeadingE_Cside );
  if(f2dHistZPVsLP)     f2dHistZPVsLP          -> Fill ( fZPCpp+fZPApp , fLeadingP );
  if(f2dHistZNVsLN)     f2dHistZNVsLN          -> Fill ( fZNCpp+fZNApp , fLeadingN );
  if(f2dHistSPDClRecoVsTrue)    f2dHistSPDClRecoVsTrue    -> Fill ( fCentrality_SPDClusters , lSPDClusters );
  if(f2dHistV0MRecoVsTrue)      f2dHistV0MRecoVsTrue      -> Fill ( fCentrality_V0M , lNchVZEROA+lNchVZEROC );
  
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

//______________________________________________________________________
Double_t AliAnalysisTaskMCPredictionsStrgVsMultVsZDC::Rapidity(Double_t E, Double_t Pz) const
{
  // Local calculation for rapidity
  Double_t ReturnValue = -100;
  if( (E - Pz + 1.e-13) != 0 && (E + Pz) != 0 ) {
    ReturnValue =  0.5*TMath::Log((E + Pz)/( E - Pz + 1.e-13));
  }
  return ReturnValue;
}