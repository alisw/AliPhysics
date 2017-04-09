/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Developers: Friederike Bock                                            *
 * Original Authors: Svein Lindal, Daniel Lohner                          *
 *                                                                        *
 * Version 1.0                                                            *
 *                                                                        *
 *                                                                        *
 * based on: on older version (see aliroot up to v5-04-42-AN)             *
 *           AliV0Reader.cxx                                              *
 *           Authors: Kathrin Koch, Kenneth Aamodt, Ana Marin             *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notic    *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class reconstructing conversion photons from V0s
//---------------------------------------------
////////////////////////////////////////////////

// The AliAODConversionPhotons will return the Position (== ESDTrack->GetID())of the positive (negative) ESDTrack 
// in the Array of ESDTracks stored in the ESDEvent via GetTrackLabelPositive() (GetTrackLabelNagative()).
// Since it is the Position in the Array it is always positive.
// After the conversion ESD --> AOD each AOD track will give you the former Position in the ESDArray via the
// Function AODTrack->GetID(). AODTracks are stored for different TrackParameter (e.g. HybridTracks, TPC only ...). No Standard 
// AODTracks (not copies of AliExternalTrackParam) will be stored with negative values of AODTrack->GetID() (Formula: (ESDTrack->GetID()+1)*-1)
// This leads to the possibility of more than one AOD track with the same negative ID. For that reason all AOD tracks are additionally flaged.
// In this analysis we should therefore only use AODTracks with positive values of AODTrack->GetID().

// If you want to find the AODTrack corresponding to the daugher track of a AliAODConversionPhoton you have to find the AODTrack with
// AODTrack->GetID() == GetTrackLabelPositive() (GetTrackLabelNagative()).

#include <vector>
#include <TGeoGlobalMagField.h>

#include "AliV0ReaderV1.h"
#include "AliKFParticle.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliPID.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "TRandom3.h"
#include "AliGenCocktailEventHeader.h"
#include "TList.h"
#include "AliKFConversionPhoton.h"
#include "AliAODConversionPhoton.h"
#include "AliConversionPhotonBase.h"
#include "TVector.h"
#include "AliKFVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TObjArray.h"

class iostream;

using namespace std;

ClassImp(AliV0ReaderV1)

//________________________________________________________________________
AliV0ReaderV1::AliV0ReaderV1(const char *name) : AliAnalysisTaskSE(name),
  fConversionCuts(NULL),
  fEventCuts(NULL),
  fConversionGammas(NULL),
  fUseImprovedVertex(kTRUE),
  fUseOwnXYZCalculation(kTRUE),
  fUseConstructGamma(kFALSE),
  kUseAODConversionPhoton(kTRUE),
  fCreateAOD(kFALSE),
  fDeltaAODBranchName("GammaConv"),
  fDeltaAODFilename("AliAODGammaConversion.root"),
  fRelabelAODs(kFALSE),
  fPreviousV0ReaderPerformsAODRelabeling(0),
  fEventIsSelected(kFALSE),
  fNumberOfPrimaryTracks(0),
  fPeriodName(""),
  fPtHardBin(0),
  fUseMassToZero(kTRUE),
  fProduceV0findingEffi(kFALSE),
  fProduceImpactParamHistograms(kFALSE),
  fCurrentInvMassPair(0),
  fImprovedPsiPair(3),
  fHistograms(NULL),
  fImpactParamHistograms(NULL),
  fHistoMCGammaPtvsR(NULL),
  fHistoMCGammaPtvsPhi(NULL),
  fHistoMCGammaPtvsEta(NULL),
  fHistoMCGammaRvsPhi(NULL),
  fHistoMCGammaRvsEta(NULL),
  fHistoMCGammaPhivsEta(NULL),
  fHistoRecMCGammaPtvsR(NULL),
  fHistoRecMCGammaPtvsPhi(NULL),
  fHistoRecMCGammaPtvsEta(NULL),
  fHistoRecMCGammaRvsPhi(NULL),
  fHistoRecMCGammaRvsEta(NULL),
  fHistoRecMCGammaPhivsEta(NULL),
  fHistoRecMCGammaMultiPt(NULL),
  fHistoRecMCGammaMultiPtvsEta(NULL),
  fHistoRecMCGammaMultiR(NULL),
  fHistoRecMCGammaMultiPhi(NULL),
  fHistoPosTrackImpactParamZ(NULL),
  fHistoPosTrackImpactParamY(NULL),
  fHistoPosTrackImpactParamX(NULL),
  fHistoPosTrackImpactParamZvsPt(NULL),
  fHistoPosTrackImpactParamYvsPt(NULL),
  fHistoPosTrackImpactParamXvsPt(NULL),
  fHistoNegTrackImpactParamZ(NULL),
  fHistoNegTrackImpactParamY(NULL),
  fHistoNegTrackImpactParamX(NULL),
  fHistoNegTrackImpactParamZvsPt(NULL),
  fHistoNegTrackImpactParamYvsPt(NULL),
  fHistoNegTrackImpactParamXvsPt(NULL),
  fHistoImpactParamZvsR(NULL),
  fHistoImpactParamZvsR2(NULL),
  fHistoPt(NULL),
  fHistoPt2(NULL),
  fHistoDCAzPhoton(NULL),
  fHistoDCAzPhoton2(NULL),
  fHistoR(NULL),
  fHistoRrecalc(NULL),
  fHistoRviaAlpha(NULL),
  fHistoRviaAlphaRecalc(NULL),
  fHistoRdiff(NULL),
  fHistoImpactParameterStudy(NULL),
  fImpactParamTree(NULL),
  fVectorFoundGammas(0),
  fCurrentFileName(""),
  fMCFileChecked(kFALSE)
{
  // Default constructor

  DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliV0ReaderV1::~AliV0ReaderV1()
{
  // default deconstructor

  if(fConversionGammas){
    fConversionGammas->Delete();// Clear Objects
    delete fConversionGammas;
    fConversionGammas=0x0;
  }
}

/*
//________________________________________________________________________
AliV0ReaderV1::AliV0ReaderV1(AliV0ReaderV1 &original) : AliAnalysisTaskSE(original),
fConversionCuts(NULL),
fConversionGammas(NULL),
fUseImprovedVertex(original.fUseImprovedVertex),
fUseOwnXYZCalculation(original.fUseOwnXYZCalculation),
fUseConstructGamma(original.fUseConstructGamma),
kUseAODConversionPhoton(original.kUseAODConversionPhoton),
fCreateAOD(original.fCreateAOD),
fDeltaAODBranchName(original.fDeltaAODBranchName),
fDeltaAODFilename(original.fDeltaAODFilename),
fEventIsSelected(original.fEventIsSelected)
{
// Default constructor

DefineInput(0, TChain::Class());
}

//____________________________________________________________
AliV0ReaderV1 &AliV0ReaderV1::operator=(const AliV0ReaderV1 &ref){
//
// Assignment operator
// Only copies pointers, object is not the owner of the references
//
if(this != &ref){
AliAnalysisTaskSE::operator=(ref);
fUseImprovedVertex=ref.fUseImprovedVertex;
fUseOwnXYZCalculation=ref.fUseOwnXYZCalculation;
fUseConstructGamma=ref.fUseConstructGamma;
kUseAODConversionPhoton=ref.kUseAODConversionPhoton;
fCreateAOD=ref.fCreateAOD;
fDeltaAODBranchName=ref.fDeltaAODBranchName;
fDeltaAODFilename=ref.fDeltaAODFilename;
fEventIsSelected=ref.fEventIsSelected;
}
return *this;
}
*/

//________________________________________________________________________
void AliV0ReaderV1::Init()
{
  // Initialize function to be called once before analysis
  if(fConversionCuts==NULL){
    if(fConversionCuts==NULL)AliError("No Conversion Cut Selection initialized");
  }
  if(fEventCuts==NULL){
    if(fEventCuts==NULL)AliError("No Event Cut Selection initialized");
  }

  if(fCreateAOD){kUseAODConversionPhoton=kTRUE;}

  if(fConversionGammas != NULL){
    delete fConversionGammas;
    fConversionGammas=NULL;
  }

  if(fConversionGammas == NULL){
    if(kUseAODConversionPhoton){
      fConversionGammas = new TClonesArray("AliAODConversionPhoton",100);}
    else{
      fConversionGammas = new TClonesArray("AliKFConversionPhoton",100);}
  }
  fConversionGammas->Delete();//Reset the TClonesArray
}

//________________________________________________________________________
void AliV0ReaderV1::UserCreateOutputObjects()
{
  // Create AODs

  if(fCreateAOD){
    if (fEventCuts){
      fDeltaAODBranchName.Append("_");
      fDeltaAODBranchName.Append(fEventCuts->GetCutNumber());
    }
    if(fConversionCuts){
      fDeltaAODBranchName.Append("_");
      fDeltaAODBranchName.Append(fConversionCuts->GetCutNumber());
      fDeltaAODBranchName.Append("_gamma");
    }
    fConversionGammas->SetName(fDeltaAODBranchName.Data());

    AddAODBranch("TClonesArray", &fConversionGammas, fDeltaAODFilename.Data());
    AliAnalysisManager::GetAnalysisManager()->RegisterExtraFile(fDeltaAODFilename.Data());
  }

  if(fProduceImpactParamHistograms){
    if(fImpactParamHistograms != NULL){
       delete fImpactParamHistograms;
       fImpactParamHistograms = NULL;
    }
    if(fImpactParamHistograms==NULL){
      fImpactParamHistograms = new TList();
      fImpactParamHistograms->SetOwner(kTRUE);
      fImpactParamHistograms->SetName(Form("ImpactParamHistograms_%s_%s",fEventCuts->GetCutNumber().Data(),fConversionCuts->GetCutNumber().Data()));
    }
    fHistoPosTrackImpactParamZ = new TH1F("fHistoPosTrackImpactParamZ","",480,-80,80);
    fHistoPosTrackImpactParamZ->SetXTitle("Z (cm)");
    fImpactParamHistograms->Add(fHistoPosTrackImpactParamZ);

    fHistoPosTrackImpactParamY = new TH1F("fHistoPosTrackImpactParamY","",720,-120,120);
    fHistoPosTrackImpactParamY->SetXTitle("Y (cm)");
    fImpactParamHistograms->Add(fHistoPosTrackImpactParamY);

    fHistoPosTrackImpactParamX = new TH1F("fHistoPosTrackImpactParamX","",30,-3,3);
    fHistoPosTrackImpactParamX->SetXTitle("X (cm)");
    fImpactParamHistograms->Add(fHistoPosTrackImpactParamX);

    fHistoNegTrackImpactParamZ = new TH1F("fHistoNegTrackImpactParamZ","",480,-80,80);
    fHistoNegTrackImpactParamZ->SetXTitle("Z (cm)");
    fImpactParamHistograms->Add(fHistoNegTrackImpactParamZ);

    fHistoNegTrackImpactParamY = new TH1F("fHistoNegTrackImpactParamY","",720,-120,120);
    fHistoNegTrackImpactParamY->SetXTitle("Y (cm)");
    fImpactParamHistograms->Add(fHistoNegTrackImpactParamY);

    fHistoNegTrackImpactParamX = new TH1F("fHistoNegTrackImpactParamX","",30,-3,3);
    fHistoNegTrackImpactParamX->SetXTitle("X (cm)");
    fImpactParamHistograms->Add(fHistoNegTrackImpactParamX);

    fHistoPosTrackImpactParamZvsPt = new TH2F("fHistoPosTrackImpactParamZvsPt","",100,0,10,480,-80,80);
    fHistoPosTrackImpactParamZvsPt->SetYTitle("Z (cm)");
    fHistoPosTrackImpactParamZvsPt->SetXTitle("Pt (GeV)");
    fImpactParamHistograms->Add(fHistoPosTrackImpactParamZvsPt);

    fHistoPosTrackImpactParamYvsPt = new TH2F("fHistoPosTrackImpactParamYvsPt","",100,0,10,720,-120,120);
    fHistoPosTrackImpactParamYvsPt->SetYTitle("Y (cm)");
    fHistoPosTrackImpactParamYvsPt->SetXTitle("Pt (GeV)");
    fImpactParamHistograms->Add(fHistoPosTrackImpactParamYvsPt);
 
    fHistoPosTrackImpactParamXvsPt = new TH2F("fHistoPosTrackImpactParamXvsPt","",100,0,10,30,-3,5);
    fHistoPosTrackImpactParamXvsPt->SetYTitle("X (cm)");
    fHistoPosTrackImpactParamXvsPt->SetXTitle("Pt (GeV)");
    fImpactParamHistograms->Add(fHistoPosTrackImpactParamXvsPt);

    fHistoNegTrackImpactParamZvsPt = new TH2F("fHistoNegTrackImpactParamZvsPt","",100,0,10,480,-80,80);
    fHistoNegTrackImpactParamZvsPt->SetYTitle("Z (cm)");
    fHistoNegTrackImpactParamZvsPt->SetXTitle("Pt (GeV)");
    fImpactParamHistograms->Add(fHistoNegTrackImpactParamZvsPt);

    fHistoNegTrackImpactParamYvsPt = new TH2F("fHistoNegTrackImpactParamYvsPt","",100,0,10,720,-120,120);
    fHistoNegTrackImpactParamYvsPt->SetYTitle("Y (cm)");
    fHistoNegTrackImpactParamYvsPt->SetXTitle("Pt (GeV)");
    fImpactParamHistograms->Add(fHistoNegTrackImpactParamYvsPt);

    fHistoNegTrackImpactParamXvsPt = new TH2F("fHistoNegTrackImpactParamXvsPt","",100,0,10,30,-3,3);
    fHistoNegTrackImpactParamXvsPt->SetYTitle("X (cm)");
    fHistoNegTrackImpactParamXvsPt->SetXTitle("Pt (GeV)");
    fImpactParamHistograms->Add(fHistoNegTrackImpactParamXvsPt);

    fHistoImpactParamZvsR = new TH2F("fHistoImpactParamZvsR","Before cuts",300,-150,150,200,0,200);
    fHistoImpactParamZvsR->SetXTitle("Z (cm)");
    fHistoImpactParamZvsR->SetYTitle("R (cm)");
    fImpactParamHistograms->Add(fHistoImpactParamZvsR);

    fHistoImpactParamZvsR2 = new TH2F("fHistoImpactParamZvsR2","After cuts",300,-150,150,200,0,200);
    fHistoImpactParamZvsR2->SetXTitle("Z (cm)");
    fHistoImpactParamZvsR2->SetYTitle("R (cm)");
    fImpactParamHistograms->Add(fHistoImpactParamZvsR2);

    fHistoPt = new TH1F("fHistoPt","Before all cuts",100,0,10);
    fHistoPt->SetXTitle("Pt (GeV)");
    fImpactParamHistograms->Add(fHistoPt);

    fHistoPt2 = new TH1F("fHistoPt2","After all cuts",100,0,10);
    fHistoPt2->SetXTitle("Pt (GeV)");
    fImpactParamHistograms->Add(fHistoPt2);

    fHistoDCAzPhoton = new TH1F("fHistoDCAzPhoton","Before cuts",20,-2,2);
    fHistoDCAzPhoton->SetXTitle("DCAz photon (cm)");
    fImpactParamHistograms->Add(fHistoDCAzPhoton);

    fHistoDCAzPhoton2 = new TH1F("fHistoDCAzPhoton2","After cuts",20,-2,2);
    fHistoDCAzPhoton2->SetXTitle("DCAz photon (cm)");
    fImpactParamHistograms->Add(fHistoDCAzPhoton2);

    fHistoR = new TH1F("fHistoR","",200,0,200);
    fHistoR->SetXTitle("Conversion radius (cm)");
    fImpactParamHistograms->Add(fHistoR);

    fHistoRrecalc = new TH1F("fHistoRrecalc","",200,0,200);
    fHistoRrecalc->SetXTitle("conversion radius (cm)");
    fImpactParamHistograms->Add(fHistoRrecalc);

    fHistoRviaAlpha = new TH1F("fHistoRviaAlpha","",200,0,200);
    fHistoRviaAlpha->SetXTitle("Conversion radius (cm)");
    fImpactParamHistograms->Add(fHistoRviaAlpha);

    fHistoRviaAlphaRecalc = new TH1F("fHistoRviaAlphaRecalc","",200,0,200);
    fHistoRviaAlphaRecalc->SetXTitle("conversion radius (cm)");
    fImpactParamHistograms->Add(fHistoRviaAlphaRecalc);

    fHistoRdiff = new TH1F("fHistoRdiff","",200,0,200);
    fHistoRdiff->SetXTitle("R_conv - R_cluster conflict (cm)");
    fImpactParamHistograms->Add(fHistoRdiff);

    fHistoImpactParameterStudy = new TH1F("fHistoImpactParameterStudy","",7,-0.5,6.5);
    fHistoImpactParameterStudy->GetXaxis()->SetBinLabel(1,"# V0s");
    fHistoImpactParameterStudy->GetXaxis()->SetBinLabel(2,"two TPC-only tracks");
    fHistoImpactParameterStudy->GetXaxis()->SetBinLabel(3,"Z cut not passed");
    fHistoImpactParameterStudy->GetXaxis()->SetBinLabel(4,"Y cut not passed");
    fHistoImpactParameterStudy->GetXaxis()->SetBinLabel(5,"R>80cm");
    fHistoImpactParameterStudy->GetXaxis()->SetBinLabel(6,"causality cut not p.");
    fHistoImpactParameterStudy->GetXaxis()->SetBinLabel(7,"# removed V0s"); // of all V0s
    fImpactParamHistograms->Add(fHistoImpactParameterStudy);

    fImpactParamTree = new TTree("fImpactParamTree","");
    fImpactParamHistograms->Add(fImpactParamTree);
  }

  if (fProduceV0findingEffi){
    TH1::AddDirectory(kFALSE);
    if(fHistograms != NULL){
      delete fHistograms;
      fHistograms  = NULL;
    }
    if(fHistograms==NULL){
      fHistograms = new TList();
      fHistograms->SetOwner(kTRUE);
      fHistograms->SetName(Form("V0FindingEfficiencyInput_%s_%s",fEventCuts->GetCutNumber().Data(),fConversionCuts->GetCutNumber().Data()));
    }

    fHistoMCGammaPtvsR         = new TH2F("MCconvGamma_Pt_R","MC converted gamma Pt vs R (|eta| < 0.9)",250,0.0,25,400,0,200);
    fHistoMCGammaPtvsR->SetXTitle("p_{MC,T} (GeV/c)");
    fHistoMCGammaPtvsR->SetYTitle("R_{MC,conv} (cm)");
    fHistograms->Add(fHistoMCGammaPtvsR);

    fHistoMCGammaPtvsEta       = new TH2F("MCconvGamma_Pt_Eta","MC converted gamma Pt vs Eta ",250,0.0,25,280,-1.4,1.4);
    fHistoMCGammaPtvsEta->SetXTitle("p_{MC,T} (GeV/c)");
    fHistoMCGammaPtvsEta->SetYTitle("#eta_{MC}");
    fHistograms->Add(fHistoMCGammaPtvsEta);

    fHistoMCGammaPtvsPhi       = new TH2F("MCconvGamma_Pt_Phi","MC converted gamma Pt vs Phi (|eta| < 0.9) ",250,0.0,25,400,0,2*TMath::Pi());
    fHistoMCGammaPtvsPhi->SetXTitle("p_{MC,T} (GeV/c)");
    fHistoMCGammaPtvsPhi->SetYTitle("#varphi_{MC} (rad)");
    fHistograms->Add(fHistoMCGammaPtvsPhi);

    fHistoMCGammaRvsPhi       = new TH2F("MCconvGamma_R_Phi","MC converted gamma R vs Phi (|eta| < 0.9) ",400,0,200,400,0,2*TMath::Pi());
    fHistoMCGammaRvsPhi->SetXTitle("R_{MC,conv} (cm)");
    fHistoMCGammaRvsPhi->SetYTitle("#varphi_{MC} (rad)");
    fHistograms->Add(fHistoMCGammaRvsPhi);

    fHistoMCGammaRvsEta       = new TH2F("MCconvGamma_R_Eta","MC converted gamma R vs Eta ",400,0,200,280,-1.4,1.4);
    fHistoMCGammaRvsEta->SetXTitle("R_{MC,conv} (cm)");
    fHistoMCGammaRvsEta->SetYTitle("#eta_{MC}");
    fHistograms->Add(fHistoMCGammaRvsEta);

    fHistoMCGammaPhivsEta       = new TH2F("MCconvGamma_Phi_Eta","MC converted gamma Phi vs Eta ",400,0,2*TMath::Pi(),280,-1.4,1.4);
    fHistoMCGammaPhivsEta->SetXTitle("#phi_{MC} (rad)");
    fHistoMCGammaPhivsEta->SetYTitle("#eta_{MC}");
    fHistograms->Add(fHistoMCGammaPhivsEta);

    fHistoRecMCGammaPtvsR       = new TH2F("RecMCconvGamma_Pt_R","rec MC converted gamma Pt vs R (|eta| < 0.9)",250,0.0,25,400,0,200);
    fHistoRecMCGammaPtvsR->SetXTitle("p_{MC,T} (GeV/c)");
    fHistoRecMCGammaPtvsR->SetYTitle("R_{MC,conv} (cm)");
    fHistograms->Add(fHistoRecMCGammaPtvsR);

    fHistoRecMCGammaPtvsEta     = new TH2F("RecMCconvGamma_Pt_Eta","rec MC converted gamma Pt vs Eta ",250,0.0,25,280,-1.4,1.4);
    fHistoRecMCGammaPtvsEta->SetXTitle("p_{MC,T} (GeV/c)");
    fHistoRecMCGammaPtvsEta->SetYTitle("#eta_{MC}");
    fHistograms->Add(fHistoRecMCGammaPtvsEta);

    fHistoRecMCGammaPtvsPhi     = new TH2F("RecMCconvGamma_Pt_Phi","rec MC converted gamma Pt vs Phi (|eta| < 0.9) ",250,0.0,25,400,0,2*TMath::Pi());
    fHistoRecMCGammaPtvsPhi->SetXTitle("p_{MC,T} (GeV/c)");
    fHistoRecMCGammaPtvsPhi->SetYTitle("#varphi_{MC} (rad)");
    fHistograms->Add(fHistoRecMCGammaPtvsPhi);

    fHistoRecMCGammaRvsPhi       = new TH2F("RecMCconvGamma_R_Phi","rec MC converted gamma R vs Phi (|eta| < 0.9) ",400,0,200,400,0,2*TMath::Pi());
    fHistoRecMCGammaRvsPhi->SetXTitle("R_{MC,conv} (cm)");
    fHistoRecMCGammaRvsPhi->SetYTitle("#varphi_{MC} (rad)");
    fHistograms->Add(fHistoRecMCGammaRvsPhi);

    fHistoRecMCGammaRvsEta      = new TH2F("RecMCconvGamma_R_Eta","rec MC converted gamma R vs Eta ",400,0,200,280,-1.4,1.4);
    fHistoRecMCGammaRvsEta->SetXTitle("R_{MC,conv} (cm)");
    fHistoRecMCGammaRvsEta->SetYTitle("#eta_{MC}");
    fHistograms->Add(fHistoRecMCGammaRvsEta);

    fHistoRecMCGammaPhivsEta     = new TH2F("RecMCconvGamma_Phi_Eta","rec MC converted gamma Phi vs Eta ",400,0,2*TMath::Pi(),280,-1.4,1.4);
    fHistoRecMCGammaPhivsEta->SetXTitle("#phi_{MC} (rad)");
    fHistoRecMCGammaPhivsEta->SetYTitle("#eta_{MC}");
    fHistograms->Add(fHistoRecMCGammaPhivsEta);

    fHistoRecMCGammaMultiPtvsEta   = new TH2F("RecMCconvGammaMulti_Pt_Eta","rec MC converted gamma (at least double counted) Pt vs Eta ",250,0.0,25,280,-1.4,1.4);
    fHistoRecMCGammaMultiPtvsEta->SetXTitle("p_{MC,T} (GeV/c)");
    fHistoRecMCGammaMultiPtvsEta->SetYTitle("#eta_{MC}");
    fHistograms->Add(fHistoRecMCGammaMultiPtvsEta);

    fHistoRecMCGammaMultiPt      = new TH1F("RecMCconvGammaMulti_Pt","rec MC converted gamma (at least double counted) Pt (|eta| < 0.9)",250,0.0,25);
    fHistoRecMCGammaMultiPt->SetXTitle("p_{MC,T} (GeV/c)");
    fHistograms->Add(fHistoRecMCGammaMultiPt);

    fHistoRecMCGammaMultiR      = new TH1F("RecMCconvGammaMulti_R","rec MC converted gamma (at least double counted) R (|eta| < 0.9)",400,0,200);
    fHistoRecMCGammaMultiR->SetXTitle("R_{MC,conv} (cm)");
    fHistograms->Add(fHistoRecMCGammaMultiR);

    fHistoRecMCGammaMultiPhi    = new TH1F("RecMCconvGammaMulti_Phi","rec MC converted gamma (at least double counted) Phi (|eta| < 0.9)",400,0,2*TMath::Pi());
    fHistoRecMCGammaMultiPhi->SetXTitle("#phi_{MC} (rad)");
    fHistograms->Add(fHistoRecMCGammaMultiPhi);

    fVectorFoundGammas.clear();
  }

}
//________________________________________________________________________
Bool_t AliV0ReaderV1::Notify(){
  
  // set file name to empty & reset check flag
  fCurrentFileName = "";
  fMCFileChecked   = kFALSE;
  
  // obtain file name from analysis manager
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if(man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    if (inputHandler){
      TTree* tree = (TTree*) inputHandler->GetTree();
      TFile* file = (TFile*) tree->GetCurrentFile();
      fCurrentFileName = file->GetName();
    }
  }
  
  // try to extract period name from file name if not given
  if (fPeriodName.CompareTo("") == 0){
    TObjArray *arr = fCurrentFileName.Tokenize("/");
    fPtHardBin = -1;
    for (Int_t i = 0; i < arr->GetEntriesFast();i++ ){
      TObjString* testObjString = (TObjString*)arr->At(i);
      if (testObjString->GetString().BeginsWith("LHC")){
        fPeriodName = testObjString->GetString();
        i = arr->GetEntriesFast();
      }
    }
    if (fPeriodName.CompareTo("")==0){
      TObjArray *arr2 = fCurrentFileName.Tokenize("__");
      for (Int_t i = 0; i < arr->GetEntriesFast();i++ ){
        TObjString* testObjString = (TObjString*)arr2->At(i);
        if (testObjString->GetString().BeginsWith("LHC")){
          fPeriodName = testObjString->GetString();
          i = arr2->GetEntriesFast();
        }
      }
    }
    if (fPeriodName.CompareTo("") != 0 && fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){ 
      fEventCuts->SetPeriodEnum (fPeriodName);
    }  
  } else {
    if(fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      fEventCuts->SetPeriodEnum (fPeriodName);
    }  
  }
  // set pt hard bin from file name
  if (  fEventCuts->GetPeriodEnum() == AliConvEventCuts::kLHC15a3a  || fEventCuts->GetPeriodEnum() == AliConvEventCuts::kLHC15a3a_plus || 
        fEventCuts->GetPeriodEnum() == AliConvEventCuts::kLHC15a3b  ||
        fEventCuts->GetPeriodEnum() == AliConvEventCuts::kLHC15g1a  || fEventCuts->GetPeriodEnum() == AliConvEventCuts::kLHC15g1b || 
        fEventCuts->GetPeriodEnum() == AliConvEventCuts::kLHC16c2  || fEventCuts->GetPeriodEnum() == AliConvEventCuts::kLHC16c2_plus  ||
        fEventCuts->GetPeriodEnum() == AliConvEventCuts::kLHC16c3a  || fEventCuts->GetPeriodEnum() == AliConvEventCuts::kLHC16c3b || 
        fEventCuts->GetPeriodEnum() == AliConvEventCuts::kLHC16c3c){
    TObjArray *arr = fCurrentFileName.Tokenize("/");
    fPtHardBin = -1;
    for (Int_t i = 0; i < arr->GetEntriesFast();i++ ){
      TObjString* testObjString = (TObjString*)arr->At(i);
      if (testObjString->GetString().BeginsWith("LHC")){
        TObjString* testObjString2 = (TObjString*)arr->At(i+1);
        fPtHardBin = testObjString2->GetString().Atoi();
        i = arr->GetEntriesFast();
      }
    }
  }

  if(!fEventCuts->GetDoEtaShift()) return kTRUE; // No Eta Shift requested, continue
  if(fEventCuts->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
    fEventCuts->GetCorrectEtaShiftFromPeriod();
    fEventCuts->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
    return kTRUE;
  } else {
    printf(" Gamma Conversion Reader %s_%s :: Eta Shift Manually Set to %f \n\n",
        (fEventCuts->GetCutNumber()).Data(),(fConversionCuts->GetCutNumber()).Data(),fEventCuts->GetEtaShift());
    fEventCuts->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
  }

  return kTRUE;
}
//________________________________________________________________________
void AliV0ReaderV1::UserExec(Option_t *option){

  AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
  if(esdEvent) {
    if (!TGeoGlobalMagField::Instance()->GetField()) esdEvent->InitMagneticField();
  }

  // Check if correctly initialized
  if(!fConversionGammas)Init();

  // User Exec
  fEventIsSelected=ProcessEvent(fInputEvent,fMCEvent);

  // AOD Output
  FillAODOutput();

}

//________________________________________________________________________
Bool_t AliV0ReaderV1::ProcessEvent(AliVEvent *inputEvent,AliMCEvent *mcEvent)
{

  //Reset the TClonesArray
  fConversionGammas->Delete();

  fInputEvent = inputEvent;
  fMCEvent    = mcEvent;

  if(!fInputEvent){
    AliError("No Input event");
    return kFALSE;
  }
  if(!fEventCuts){AliError("No EventCuts");return kFALSE;}
  if(!fConversionCuts){AliError("No ConversionCuts");return kFALSE;}


  // Count Primary Tracks Event
  CountTracks();

  // Event Cuts
  if(!fEventCuts->EventIsSelected(fInputEvent,fMCEvent)){
    if (fEventCuts->GetEventQuality() == 2 && !fMCFileChecked ){
      cout << "ERROR with MC reading for: "<<  fCurrentFileName.Data() << endl;
      fMCFileChecked                  = kTRUE;
    }  
    return kFALSE;
  }  
  // Set Magnetic Field
  AliKFParticle::SetField(fInputEvent->GetMagneticField());

  if(fInputEvent->IsA()==AliAODEvent::Class() && fProduceV0findingEffi){
    fProduceV0findingEffi = kFALSE;
    AliWarning("V0finding effi cannot be run on AODs ");
  }

  if(fProduceV0findingEffi){
    CreatePureMCHistosForV0FinderEffiESD();
    fVectorFoundGammas.clear();
  }

  if(fInputEvent->IsA()==AliESDEvent::Class()){
    ProcessESDV0s();
  }
  if(fInputEvent->IsA()==AliAODEvent::Class()){
    GetAODConversionGammas();
  }

  return kTRUE;
}
///________________________________________________________________________
void AliV0ReaderV1::FillAODOutput()
{
   // Fill AOD Output with reconstructed Photons

  if(fInputEvent->IsA()==AliESDEvent::Class()){
    ///Make sure delta aod is filled if standard aod is filled (for synchronization when reading aod with standard aod)
    if(fCreateAOD) {
      AliAODHandler * aodhandler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
      if (aodhandler && aodhandler->GetFillAOD()) {
        AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillExtension(kTRUE);
        //PostData(0, fConversionGammas);

      }
    }
  }
}

///________________________________________________________________________
const AliExternalTrackParam *AliV0ReaderV1::GetExternalTrackParam(AliESDv0 *fCurrentV0,Int_t &tracklabel,Int_t charge){

  // Get External Track Parameter with given charge

  if(!(charge==1||charge==-1)){AliError("Charge not defined");return 0x0;}

  if(fConversionCuts->GetV0FinderSameSign()==1){
    if(fCurrentV0){
      if((fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetPindex()))->Charge()!=(fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetNindex()))->Charge())return 0x0;
      if((fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetPindex()))->Charge()==(fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetNindex()))->Charge()){
        if(charge==1){
          tracklabel=fCurrentV0->GetPindex();
          return fCurrentV0->GetParamP();
        }else{
          tracklabel=fCurrentV0->GetNindex();
          return fCurrentV0->GetParamN();
        }
      }
    }
  }else if(fConversionCuts->GetV0FinderSameSign()==2){
    if(fCurrentV0){
        if(charge==1){
          tracklabel=fCurrentV0->GetPindex();
          return fCurrentV0->GetParamP();
        }else{
          tracklabel=fCurrentV0->GetNindex();
          return fCurrentV0->GetParamN();
        }
    }
  }else{
    // Check for sign flip
    if(fCurrentV0){
      if(!fCurrentV0->GetParamN()||!fCurrentV0->GetParamP())return 0x0;
      if(!fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetNindex())||!fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetPindex()))return 0x0;
      if((fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetPindex()))->Charge()==charge){
        tracklabel=fCurrentV0->GetPindex();
        return fCurrentV0->GetParamP();}
      if((fConversionCuts->GetTrack(fInputEvent,fCurrentV0->GetNindex()))->Charge()==charge){
        tracklabel=fCurrentV0->GetNindex();
        return fCurrentV0->GetParamN();}
    }
  }
  return 0x0;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::ProcessESDV0s()
{

  // Process ESD V0s for conversion photon reconstruction
  AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(fInputEvent);

  AliKFConversionPhoton *fCurrentMotherKFCandidate=NULL;

  if(fESDEvent){
    for(Int_t currentV0Index=0;currentV0Index<fESDEvent->GetNumberOfV0s();currentV0Index++){
      AliESDv0 *fCurrentV0=(AliESDv0*)(fESDEvent->GetV0(currentV0Index));
      if(!fCurrentV0){
        printf("Requested V0 does not exist");
        continue;
      }

      fCurrentMotherKFCandidate=ReconstructV0(fCurrentV0,currentV0Index);

      if(fCurrentMotherKFCandidate){
        // Add Gamma to the TClonesArray

        if(kUseAODConversionPhoton){
          new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliAODConversionPhoton(fCurrentMotherKFCandidate);
          AliAODConversionPhoton * currentConversionPhoton = (AliAODConversionPhoton*)(fConversionGammas->At(fConversionGammas->GetEntriesFast()-1));
          currentConversionPhoton->SetMass(fCurrentMotherKFCandidate->M());
          if (fUseMassToZero) currentConversionPhoton->SetMassToZero();
          currentConversionPhoton->SetInvMassPair(fCurrentInvMassPair);
        } else {
          new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliKFConversionPhoton(*fCurrentMotherKFCandidate);
        }

        delete fCurrentMotherKFCandidate;
        fCurrentMotherKFCandidate=NULL;
      }
    }
  }
  return kTRUE;
}

///________________________________________________________________________
AliKFConversionPhoton *AliV0ReaderV1::ReconstructV0(AliESDv0 *fCurrentV0,Int_t currentV0Index)
{
  //   cout << currentV0Index << endl;
  // Reconstruct conversion photon from ESD v0
  fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kPhotonIn);

  //checks if on the fly mode is set
  if(!fConversionCuts->SelectV0Finder(fCurrentV0->GetOnFlyStatus())){
    fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kOnFly);
    return 0x0;
  }

  if (fMCEvent && fProduceV0findingEffi ) FillRecMCHistosForV0FinderEffiESD(fCurrentV0);

  // TrackLabels
  Int_t currentTrackLabels[2]={-1,-1};

  // Get Daughter KF Particles

  const AliExternalTrackParam *fCurrentExternalTrackParamPositive=GetExternalTrackParamP(fCurrentV0,currentTrackLabels[0]);
  //    cout << fCurrentExternalTrackParamPositive << "\t" << currentTrackLabels[0] << endl;
  const AliExternalTrackParam *fCurrentExternalTrackParamNegative=GetExternalTrackParamN(fCurrentV0,currentTrackLabels[1]);
  //    cout << fCurrentExternalTrackParamNegative << "\t" << currentTrackLabels[1] << endl;
  if(!fCurrentExternalTrackParamPositive||!fCurrentExternalTrackParamNegative)return 0x0;

  // Apply some Cuts before Reconstruction

  AliVTrack * posTrack = fConversionCuts->GetTrack(fInputEvent,currentTrackLabels[0]);
  AliVTrack * negTrack = fConversionCuts->GetTrack(fInputEvent,currentTrackLabels[1]);
  if(!negTrack || !posTrack) {
    fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kNoTracks);
    return 0x0;
  }
  // Track Cuts
  if(!fConversionCuts->TracksAreSelected(negTrack, posTrack)){
    fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kTrackCuts);
    return 0x0;
  }

  fConversionCuts->FillV0EtaBeforedEdxCuts(fCurrentV0->Eta());
  if (!fConversionCuts->dEdxCuts(posTrack)) {
    fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kdEdxCuts);
    return 0x0;
  }
  // PID Cuts
  if(!fConversionCuts->dEdxCuts(negTrack)) {
    fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kdEdxCuts);
    return 0x0;
  }
  fConversionCuts->FillV0EtaAfterdEdxCuts(fCurrentV0->Eta());
  // Reconstruct Photon
  AliKFConversionPhoton *fCurrentMotherKF=NULL;
  //    fUseConstructGamma = kFALSE;
  //    cout << "construct gamma " << endl;
  AliKFParticle fCurrentNegativeKFParticle(*(fCurrentExternalTrackParamNegative),11);
  //    cout << fCurrentExternalTrackParamNegative << "\t" << endl;
  AliKFParticle fCurrentPositiveKFParticle(*(fCurrentExternalTrackParamPositive),-11);
  //    cout << fCurrentExternalTrackParamPositive << "\t"  << endl;
  //    cout << currentTrackLabels[0] << "\t" << currentTrackLabels[1] << endl;
  //    cout << "construct gamma " <<fUseConstructGamma << endl;

  // Reconstruct Gamma
  if(fUseConstructGamma){
    fCurrentMotherKF = new AliKFConversionPhoton();
    fCurrentMotherKF->ConstructGamma(fCurrentNegativeKFParticle,fCurrentPositiveKFParticle);
  }else{
    fCurrentMotherKF = new AliKFConversionPhoton(fCurrentNegativeKFParticle,fCurrentPositiveKFParticle);
    fCurrentMotherKF->SetMassConstraint(0,0.0001);
  }

  // Set Track Labels

  fCurrentMotherKF->SetTrackLabels(currentTrackLabels[0],currentTrackLabels[1]);

  // Set V0 index

  fCurrentMotherKF->SetV0Index(currentV0Index);

  //Set MC Label
  if(fMCEvent){

    AliStack *fMCStack= fMCEvent->Stack();

    Int_t labelp=TMath::Abs(fConversionCuts->GetTrack(fInputEvent,fCurrentMotherKF->GetTrackLabelPositive())->GetLabel());
    Int_t labeln=TMath::Abs(fConversionCuts->GetTrack(fInputEvent,fCurrentMotherKF->GetTrackLabelNegative())->GetLabel());

//     cout << "rec: " <<  currentTrackLabels[0] << "\t" << currentTrackLabels[1] << endl;
//     cout << "recProp: " <<  fCurrentMotherKF->GetTrackLabelPositive() << "\t" << fCurrentMotherKF->GetTrackLabelNegative() << endl;
//     cout << "MC: " <<  labeln << "\t" << labelp << endl;

    TParticle *fNegativeMCParticle = 0x0;
    if(labeln>-1) fNegativeMCParticle = fMCStack->Particle(labeln);
    TParticle *fPositiveMCParticle = 0x0;
    if(labelp>-1) fPositiveMCParticle = fMCStack->Particle(labelp);

    if(fPositiveMCParticle&&fNegativeMCParticle){
      fCurrentMotherKF->SetMCLabelPositive(labelp);
      fCurrentMotherKF->SetMCLabelNegative(labeln);
    }
  }

  // Update Vertex (moved for same eta compared to old)
  //      cout << currentV0Index <<" \t before: \t" << fCurrentMotherKF->GetPx() << "\t" << fCurrentMotherKF->GetPy() << "\t" << fCurrentMotherKF->GetPz()  << endl;
  if(fUseImprovedVertex == kTRUE){
    AliKFVertex primaryVertexImproved(*fInputEvent->GetPrimaryVertex());
    //        cout << "Prim Vtx: " << primaryVertexImproved.GetX() << "\t" << primaryVertexImproved.GetY() << "\t" << primaryVertexImproved.GetZ() << endl;
    primaryVertexImproved+=*fCurrentMotherKF;
    fCurrentMotherKF->SetProductionVertex(primaryVertexImproved);
  }
  // SetPsiPair
  Double_t convpos[3]={0,0,0};
  if (fImprovedPsiPair == 0){
    Double_t PsiPair=GetPsiPair(fCurrentV0,fCurrentExternalTrackParamPositive,fCurrentExternalTrackParamNegative, convpos);
    fCurrentMotherKF->SetPsiPair(PsiPair);
  }

  // Recalculate ConversionPoint
  Double_t dca[2]={0,0};
  if(fUseOwnXYZCalculation){
    //    Double_t convpos[3]={0,0,0};
    if(!GetConversionPoint(fCurrentExternalTrackParamPositive,fCurrentExternalTrackParamNegative,convpos,dca)){
      fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kConvPointFail);
      delete fCurrentMotherKF;
      fCurrentMotherKF=NULL;
      return 0x0;
    }

    fCurrentMotherKF->SetConversionPoint(convpos);
  }


  // SetPsiPair
   if (fImprovedPsiPair >= 1){
     // the propagation can be more precise after the precise conversion point calculation
     Double_t PsiPair=GetPsiPair(fCurrentV0,fCurrentExternalTrackParamPositive,fCurrentExternalTrackParamNegative,convpos);
     fCurrentMotherKF->SetPsiPair(PsiPair);
     //cout<<" GetPsiPair::"<<fCurrentMotherKF->GetPsiPair() <<endl;
   }

  if(fCurrentMotherKF->GetNDF() > 0.)
    fCurrentMotherKF->SetChi2perNDF(fCurrentMotherKF->GetChi2()/fCurrentMotherKF->GetNDF());   //->Photon is created before all chi2 relevant changes are performed, set it "by hand"


  // Set Dilepton Mass (moved down for same eta compared to old)
  fCurrentMotherKF->SetMass(fCurrentMotherKF->M());

  // Calculating invariant mass
  Double_t mass=-99.0, mass_width=-99.0, Pt=-99.0, Pt_width=-99.0;
  AliKFParticle fCurrentMotherKFForMass(fCurrentNegativeKFParticle,fCurrentPositiveKFParticle);
  fCurrentMotherKFForMass.GetMass(mass,mass_width);
  fCurrentMotherKFForMass.GetPt(Pt,Pt_width);
  fCurrentInvMassPair=mass;

  // apply possible Kappa cut
  if (!fConversionCuts->KappaCuts(fCurrentMotherKF,fInputEvent)){
    fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kdEdxCuts);
    delete fCurrentMotherKF;
    fCurrentMotherKF=NULL;
    return 0x0;
  }
  
  // Apply Photon Cuts
  if(!fConversionCuts->PhotonCuts(fCurrentMotherKF,fInputEvent)){
    fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kPhotonCuts);
    delete fCurrentMotherKF;
    fCurrentMotherKF=NULL;
    return 0x0;
  }

  //    cout << currentV0Index <<" \t after: \t" <<fCurrentMotherKF->GetPx() << "\t" << fCurrentMotherKF->GetPy() << "\t" << fCurrentMotherKF->GetPz()  << endl;

  if(fProduceImpactParamHistograms) FillImpactParamHistograms(posTrack, negTrack, fCurrentV0, fCurrentMotherKF);

  fConversionCuts->FillPhotonCutIndex(AliConversionPhotonCuts::kPhotonOut);
  return fCurrentMotherKF;
}

///________________________________________________________________________
Double_t AliV0ReaderV1::GetPsiPair(const AliESDv0* v0, const AliExternalTrackParam *positiveparam,const AliExternalTrackParam *negativeparam,const Double_t convpos[3]) const {
  //
  // Angle between daughter momentum plane and plane
  //

  AliExternalTrackParam nt(*negativeparam);
  AliExternalTrackParam pt(*positiveparam);

  Float_t magField = fInputEvent->GetMagneticField();

  Double_t xyz[3] = {0.,0.,0.};
  if (fImprovedPsiPair==0 ) {
    v0->GetXYZ(xyz[0],xyz[1],xyz[2]);
  } else if  (fImprovedPsiPair>=1 ) {
    xyz[0]= convpos[0];
    xyz[1]= convpos[1];
    xyz[2]= convpos[2];
  }
 

  // Double_t pPlus[3]  = {pt.Px(),pt.Py(),pt.Pz()};
  // Double_t pMinus[3] = {nt.Px(),nt.Py(),nt.Pz()};

  // Double_t u[3] = {pPlus[0]+pMinus[0],pPlus[1]+pMinus[1],pPlus[2]+pMinus[2]};
  // Double_t normu = sqrt( (u[0]*u[0]) + (u[1]*u[1]) + (u[2]*u[2]) );

  // u[0] = u[0] / normu;
  // u[1] = u[1] / normu;
  // u[2] = u[2] / normu;

  // Double_t normpPlus = sqrt( (pPlus[0]*pPlus[0]) + (pPlus[1]*pPlus[1]) + (pPlus[2]*pPlus[2]) );
  // Double_t normpMinus = sqrt( (pMinus[0]*pMinus[0]) + (pMinus[1]*pMinus[1]) + (pMinus[2]*pMinus[2]) );

  // pPlus[0] = pPlus[0] / normpPlus;
  // pPlus[1] = pPlus[1] / normpPlus;
  // pPlus[2] = pPlus[2] / normpPlus;

  // pMinus[0] = pMinus[0] / normpMinus;
  // pMinus[1] = pMinus[1] / normpMinus;
  // pMinus[2] = pMinus[2] / normpMinus;

  // Double_t v[3] = {0,0,0}; // pPlus X pMinus
  // v[0] = (pPlus[1]*pMinus[2]) - (pPlus[2]*pMinus[1]);
  // v[1] = (pPlus[2]*pMinus[0]) - (pPlus[0]*pMinus[2]);
  // v[2] = (pPlus[0]*pMinus[1]) - (pPlus[1]*pMinus[0]);

  // Double_t w[3] = {0,0,0}; // u X v
  // w[0] = (u[1]*v[2]) - (u[2]*v[1]);
  // w[1] = (u[2]*v[0]) - (u[0]*v[2]);
  // w[2] = (u[0]*v[1]) - (u[1]*v[0]);

  // Double_t z[3] = {0,0,1};
  // Double_t wc[3] = {0,0,0}; // u X v
  // wc[0] = (u[1]*z[2]) - (u[2]*z[1]);
  // wc[1] = (u[2]*z[0]) - (u[0]*z[2]);
  // wc[2] = (u[0]*z[1]) - (u[1]*z[0]);

  // Double_t PhiV = TMath::ACos((w[0]*wc[0]) + (w[1]*wc[1]) + (w[2]*wc[2]));
  //return TMath::Abs(PhiV);


  // TVector3 pPlus(pt.Px(),pt.Py(),pt.Pz());
  // TVector3 pMinus(nt.Px(),nt.Py(),nt.Pz());

  // TVector3 u = pMinus + pPlus;
  // u = u*(1/u.Mag());

  // TVector3 pHPlus = pPlus*(1/pPlus.Mag());
  // TVector3 pHMinus = pMinus*(1/pMinus.Mag());

  // TVector3 v = pHPlus.Cross(pHMinus);
  // TVector3 w = u.Cross(v);
  // TVector3 z(0,0,1);
  // TVector3 wc = u.Cross(z);

  // Double_t PhiV = w * wc;

  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};

  v0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
  v0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter;

  Double_t deltat = 1.;
  deltat = TMath::ATan(mp[2]/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1])+1.e-13)) - TMath::ATan(mn[2]/(TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1])+1.e-13));//difference of angles of the two daughter tracks with z-axis
  Double_t radiussum = TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]) + 50;//radius to which tracks shall be propagated

  Double_t momPosProp[3] = {0,0,0};
  Double_t momNegProp[3] = {0,0,0};

  Double_t psiPair = 4.;
  // cout<< "Momentum before propagation at radius ::"<< TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1])<< endl;
   // cout<< mn[0]<< " " <<mn[1]<<" "<< mn[2]<<" "<<TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1] +mn[2]*mn[2])<< endl;
   // cout<< mp[0]<< " " <<mp[1]<<" "<< mp[2]<<" "<<TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1] +mp[2]*mp[2])<< endl;
  Double_t pEle,pPos;

  if (fImprovedPsiPair==1 || fImprovedPsiPair==0 ){
    if(nt.PropagateTo(radiussum,magField) == 0) return psiPair; //propagate tracks to the outside -> Better Purity and Efficiency
    if(pt.PropagateTo(radiussum,magField) == 0) return psiPair; //propagate tracks to the outside -> Better Purity and Efficiency

    pt.GetPxPyPz(momPosProp);//Get momentum vectors of tracks after propagation
    nt.GetPxPyPz(momNegProp);
 } else if (fImprovedPsiPair>=2) {
    momPosProp[0] = pt.GetParameterAtRadius(radiussum,magField,3);
    momPosProp[1] = pt.GetParameterAtRadius(radiussum,magField,4);
    momPosProp[2] = pt.GetParameterAtRadius(radiussum,magField,5);

    momNegProp[0] = nt.GetParameterAtRadius(radiussum,magField,3);
    momNegProp[1] = nt.GetParameterAtRadius(radiussum,magField,4);
    momNegProp[2] = nt.GetParameterAtRadius(radiussum,magField,5);
    pEle = TMath::Sqrt(momNegProp[0]*momNegProp[0]+momNegProp[1]*momNegProp[1]+momNegProp[2]*momNegProp[2]);//absolute momentum value of negative daughter

    pPos = TMath::Sqrt(momPosProp[0]*momPosProp[0]+momPosProp[1]*momPosProp[1]+momPosProp[2]*momPosProp[2]);//absolute momentum value of positive daughter

    if ( (pEle==0 || pPos==0) &&  fImprovedPsiPair==3) {
      radiussum = TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]) + 30;
      momPosProp[0] = pt.GetParameterAtRadius(radiussum,magField,3);
      momPosProp[1] = pt.GetParameterAtRadius(radiussum,magField,4);
      momPosProp[2] = pt.GetParameterAtRadius(radiussum,magField,5);

      momNegProp[0] = nt.GetParameterAtRadius(radiussum,magField,3);
      momNegProp[1] = nt.GetParameterAtRadius(radiussum,magField,4);
      momNegProp[2] = nt.GetParameterAtRadius(radiussum,magField,5);

    }
  }

  pEle = TMath::Sqrt(momNegProp[0]*momNegProp[0]+momNegProp[1]*momNegProp[1]+momNegProp[2]*momNegProp[2]);//absolute momentum value of negative daughter
  pPos = TMath::Sqrt(momPosProp[0]*momPosProp[0]+momPosProp[1]*momPosProp[1]+momPosProp[2]*momPosProp[2]);//absolute momentum value of positive daughter

  Double_t scalarproduct =
    momPosProp[0]*momNegProp[0]+momPosProp[1]*momNegProp[1]+momPosProp[2]*momNegProp[2];//scalar product of propagated positive and negative daughters' momenta

  if (pEle==0 || pPos==0) return psiPair;


  Double_t chipair = TMath::ACos(scalarproduct/(pEle*pPos));//Angle between propagated daughter tracks

  //    psiPair =  TMath::Abs(TMath::ASin(deltat/chipair));
  psiPair = TMath::ASin(deltat/chipair);
  return psiPair;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::GetHelixCenter(const AliExternalTrackParam *track,Double_t center[2]){

  // Get Center of the helix track parametrization

  Int_t charge=track->Charge();
  Double_t b=fInputEvent->GetMagneticField();

  Double_t  helix[6];
  track->GetHelixParameters(helix,b);

  Double_t xpos =  helix[5];
  Double_t ypos =  helix[0];
  Double_t radius = TMath::Abs(1./helix[4]);
  Double_t phi = helix[2];

  if(phi < 0){
    phi = phi + 2*TMath::Pi();
  }

  phi -= TMath::Pi()/2.;
  Double_t xpoint =  radius * TMath::Cos(phi);
  Double_t ypoint =  radius * TMath::Sin(phi);

  if(b<0){
    if(charge > 0){
      xpoint = - xpoint;
      ypoint = - ypoint;
    }
  }
  if(b>0){
    if(charge < 0){
      xpoint = - xpoint;
      ypoint = - ypoint;
    }
  }
  center[0] =  xpos + xpoint;
  center[1] =  ypos + ypoint;

  return 1;
}
///________________________________________________________________________
Bool_t AliV0ReaderV1::GetConversionPoint(const AliExternalTrackParam *pparam,const AliExternalTrackParam *nparam,Double_t convpos[3],Double_t dca[2]){

  // Recalculate Conversion Point

  if(!pparam||!nparam)return kFALSE;

  Double_t helixcenterpos[2];
  GetHelixCenter(pparam,helixcenterpos);

  Double_t helixcenterneg[2];
  GetHelixCenter(nparam,helixcenterneg);

  Double_t helixpos[6];
  pparam->GetHelixParameters(helixpos,fInputEvent->GetMagneticField());
  Double_t posradius = TMath::Abs(1./helixpos[4]);

  Double_t helixneg[6];
  nparam->GetHelixParameters(helixneg,fInputEvent->GetMagneticField());
  Double_t negradius = TMath::Abs(1./helixneg[4]);

  // Calculate xy-position

  Double_t xpos = helixcenterpos[0];
  Double_t ypos = helixcenterpos[1];
  Double_t xneg = helixcenterneg[0];
  Double_t yneg = helixcenterneg[1];

  convpos[0] = (xpos*negradius + xneg*posradius)/(negradius+posradius);
  convpos[1] = (ypos*negradius+ yneg*posradius)/(negradius+posradius);


  // Calculate Track XY vertex position

  Double_t deltaXPos = convpos[0] -  xpos;
  Double_t deltaYPos = convpos[1] -  ypos;

  Double_t deltaXNeg = convpos[0] -  xneg;
  Double_t deltaYNeg = convpos[1] -  yneg;

  Double_t alphaPos =  TMath::Pi() + TMath::ATan2(-deltaYPos,-deltaXPos);
  Double_t alphaNeg =  TMath::Pi() + TMath::ATan2(-deltaYNeg,-deltaXNeg);

  Double_t vertexXNeg =  xneg +  TMath::Abs(negradius)*TMath::Cos(alphaNeg);
  Double_t vertexYNeg =  yneg +  TMath::Abs(negradius)*TMath::Sin(alphaNeg);

  Double_t vertexXPos =  xpos +  TMath::Abs(posradius)*TMath::Cos(alphaPos);
  Double_t vertexYPos =  ypos +  TMath::Abs(posradius)*TMath::Sin(alphaPos);

  AliExternalTrackParam p(*pparam);
  AliExternalTrackParam n(*nparam);

  TVector2 vertexPos(vertexXPos,vertexYPos);
  TVector2 vertexNeg(vertexXNeg,vertexYNeg);

  // Convert to local coordinate system
  TVector2 vertexPosRot=vertexPos.Rotate(-p.GetAlpha());
  TVector2 vertexNegRot=vertexNeg.Rotate(-n.GetAlpha());

  // Propagate Track Params to Vertex

  if(!p.PropagateTo(vertexPosRot.X(),fInputEvent->GetMagneticField()))return kFALSE;
  if(!n.PropagateTo(vertexNegRot.X(),fInputEvent->GetMagneticField()))return kFALSE;

  // Check whether propagation was sucessful

  if(TMath::Abs(vertexPos.Mod()-TMath::Sqrt(p.GetX()*p.GetX()+p.GetY()*p.GetY()))>0.01)return kFALSE;
  if(TMath::Abs(vertexNeg.Mod()-TMath::Sqrt(n.GetX()*n.GetX()+n.GetY()*n.GetY()))>0.01)return kFALSE;

  // Calculate z position

  convpos[2] = (p.GetZ()*negradius+n.GetZ()*posradius)/(negradius+posradius);

  // Calculate DCA
  TVector2 vdca=vertexPos-vertexNeg;
  dca[0]=vdca.Mod();
  dca[1]=TMath::Abs(n.GetZ()-p.GetZ());

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliV0ReaderV1::GetAODConversionGammas(){

  // Get reconstructed conversion photons from satellite AOD file

  AliAODEvent *fAODEvent=dynamic_cast<AliAODEvent*>(fInputEvent);

  if(fAODEvent){

    if(fConversionGammas == NULL){
      fConversionGammas = new TClonesArray("AliAODConversionPhoton",100);
    }
    fConversionGammas->Delete();//Reset the TClonesArray

    //Get Gammas from satellite AOD gamma branch

    AliAODConversionPhoton *gamma=0x0;

    TClonesArray *fInputGammas=dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(fDeltaAODBranchName.Data()));

    if(!fInputGammas){
      FindDeltaAODBranchName();
      fInputGammas=dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(fDeltaAODBranchName.Data()));}
    if(!fInputGammas){AliError("No Gamma Satellites found");return kFALSE;}
    // Apply Selection Cuts to Gammas and create local working copy
    if(fInputGammas){
      for(Int_t i=0;i<fInputGammas->GetEntriesFast();i++){
        gamma=dynamic_cast<AliAODConversionPhoton*>(fInputGammas->At(i));
        if(gamma){
        if(fRelabelAODs)RelabelAODPhotonCandidates(gamma);
        if(fConversionCuts->PhotonIsSelected(gamma,fInputEvent)){
          new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliAODConversionPhoton(*gamma);}
        }
      }
    }
  }

  if(fConversionGammas->GetEntries()){return kTRUE;}

  return kFALSE;
}

//________________________________________________________________________
void AliV0ReaderV1::FindDeltaAODBranchName(){

  // Find delta AOD branchname containing reconstructed photons

  TList *list=fInputEvent->GetList();
  for(Int_t ii=0;ii<list->GetEntries();ii++){
    TString name((list->At(ii))->GetName());
    if(name.BeginsWith(fDeltaAODBranchName)&&name.EndsWith("gamma")){
      fDeltaAODBranchName=name;
      AliInfo(Form("Set DeltaAOD BranchName to: %s",fDeltaAODBranchName.Data()));
      return;
    }
  }
}

//________________________________________________________________________
void AliV0ReaderV1::RelabelAODPhotonCandidates(AliAODConversionPhoton *PhotonCandidate){

  if(fPreviousV0ReaderPerformsAODRelabeling == 2) return;
  else if(fPreviousV0ReaderPerformsAODRelabeling == 0){
    printf("Running AODs! Determine if V0Reader '%s' should perform relabeling\n",this->GetName());
    TObjArray* obj = (TObjArray*)AliAnalysisManager::GetAnalysisManager()->GetTasks();
    Int_t iPosition = obj->IndexOf(this);
    Bool_t prevV0ReaderRunningButNotRelabeling = kFALSE;
    for(Int_t i=iPosition-1; i>=0; i--){
     if( (obj->At(i))->IsA() == AliV0ReaderV1::Class()){
       AliV0ReaderV1* tempReader = (AliV0ReaderV1*) obj->At(i);
       if( tempReader->AreAODsRelabeled() && tempReader->IsReaderPerformingRelabeling() == 1){
         fPreviousV0ReaderPerformsAODRelabeling = 2;
         prevV0ReaderRunningButNotRelabeling = kFALSE;
         printf("V0Reader '%s' is running before this V0Reader '%s': do _NOT_ relabel AODs by current reader!\n",tempReader->GetName(),this->GetName());
         break;
       }else prevV0ReaderRunningButNotRelabeling = kTRUE;
     }
    }
    if(prevV0ReaderRunningButNotRelabeling) AliFatal(Form("There are V0Readers before '%s', but none of them is relabeling!",this->GetName()));

    if(fPreviousV0ReaderPerformsAODRelabeling == 2) return;
    else{
      printf("This V0Reader '%s' is first to be processed: do relabel AODs by current reader!\n",this->GetName());
      fPreviousV0ReaderPerformsAODRelabeling = 1;
    }
  }

  if(fPreviousV0ReaderPerformsAODRelabeling != 1) AliFatal(Form("In %s: fPreviousV0ReaderPerformsAODRelabeling = '%i' - while it should be impossible it is something different than '1'!",this->GetName(),fPreviousV0ReaderPerformsAODRelabeling));

  // Relabeling For AOD Event
  // ESDiD -> AODiD
  // MCLabel -> AODMCLabel
  Bool_t AODLabelPos = kFALSE;
  Bool_t AODLabelNeg = kFALSE;

  for(Int_t i = 0; i<fInputEvent->GetNumberOfTracks();i++){
    AliAODTrack *tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
    if(!AODLabelPos){
      if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelPositive() ){
        PhotonCandidate->SetMCLabelPositive(TMath::Abs(tempDaughter->GetLabel()));
        PhotonCandidate->SetLabelPositive(i);
        AODLabelPos = kTRUE;
      }
    }
    if(!AODLabelNeg){
      if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()){
        PhotonCandidate->SetMCLabelNegative(TMath::Abs(tempDaughter->GetLabel()));
        PhotonCandidate->SetLabelNegative(i);
        AODLabelNeg = kTRUE;
      }
    }
    if(AODLabelNeg && AODLabelPos){
      return;
    }
  }
  if(!AODLabelPos || !AODLabelNeg){
    AliError(Form("NO AOD Daughters Found Pos: %i %i Neg: %i %i, setting all labels to -999999",AODLabelPos,PhotonCandidate->GetTrackLabelPositive(),AODLabelNeg,PhotonCandidate->GetTrackLabelNegative()));
    if(!AODLabelNeg){
      PhotonCandidate->SetMCLabelNegative(-999999);
      PhotonCandidate->SetLabelNegative(-999999);
    }
    if(!AODLabelPos){
      PhotonCandidate->SetMCLabelPositive(-999999);
      PhotonCandidate->SetLabelPositive(-999999);
    }
  }

}

//************************************************************************
// This function counts the number of primary tracks in the event, the 
// implementation for ESD and AOD had to be different due to the different
// filters which are already applied on AOD level the so-called filter 
// bits, we tried to replicate the conditions in both but an exact match 
// could not be reached. The cuts are similar to the ones used by the 
// jet-group
//************************************************************************
//________________________________________________________________________
void AliV0ReaderV1::CountTracks(){

  // In Case of MC count only MB tracks
  // if(fMCEvent && fInputEvent->IsA()==AliESDEvent::Class()){
  //    fEventCuts->GetNotRejectedParticles(1,NULL,fMCEvent);
  // }
  // if(fMCEvent && fInputEvent->IsA()==AliAODEvent::Class()){
  //    fEventCuts->GetNotRejectedParticles(1,NULL,fInputEvent);
  // }

  if(fInputEvent->IsA()==AliESDEvent::Class()){
    static AliESDtrackCuts *EsdTrackCuts = 0x0;
    static int prevRun = -1;
    // Using standard function for setting Cuts
    Int_t runNumber = fInputEvent->GetRunNumber();

    if (prevRun!=runNumber) {
      delete EsdTrackCuts;
      EsdTrackCuts = 0;
      prevRun = runNumber;
    }
    if (!EsdTrackCuts) {
      // if LHC11a or earlier or if LHC13g or if LHC12a-i -> use 2010 cuts
      if( (runNumber<=146860) || (runNumber>=197470 && runNumber<=197692) || (runNumber>=172440 && runNumber<=193766) ){
	EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
	// else if run2 data use 2015 PbPb cuts
      }else if (runNumber>=209122){
	EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb();
	// else use 2011 version of track cuts
      }else{
	EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
      }
      EsdTrackCuts->SetMaxDCAToVertexZ(2);
      EsdTrackCuts->SetEtaRange(-0.8, 0.8);
      EsdTrackCuts->SetPtRange(0.15);
    }
    fNumberOfPrimaryTracks = 0;
    for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
      AliESDtrack* curTrack = (AliESDtrack*) fInputEvent->GetTrack(iTracks);
      if(!curTrack) continue;
      if(!EsdTrackCuts->AcceptTrack(curTrack)) continue;
      //if(fMCEvent && !(fEventCuts->IsParticleFromBGEvent(TMath::Abs(curTrack->GetLabel()),fMCEvent->Stack(),fInputEvent))) continue;
      fNumberOfPrimaryTracks++;
    }
  }
  else if(fInputEvent->IsA()==AliAODEvent::Class()){
    fNumberOfPrimaryTracks = 0;
    for(Int_t iTracks = 0; iTracks<fInputEvent->GetNumberOfTracks(); iTracks++){
      AliAODTrack* curTrack = (AliAODTrack*) fInputEvent->GetTrack(iTracks);
      if(curTrack->GetID()<0) continue; // Avoid double counting of tracks
      if(!curTrack->IsHybridGlobalConstrainedGlobal()) continue;
      if(TMath::Abs(curTrack->Eta())>0.8) continue;
      if(curTrack->Pt()<0.15) continue;
      //if(fMCEvent && !(fEventCuts->IsParticleFromBGEvent(TMath::Abs(curTrack->GetLabel()),NULL,fInputEvent))) continue;
      //if(TMath::Abs(curTrack->ZAtDCA())>2) continue; // Only Set For TPCOnly tracks
      fNumberOfPrimaryTracks++;
    }
  }

  return;
}

///________________________________________________________________________
Bool_t AliV0ReaderV1::ParticleIsConvertedPhoton(AliStack *MCStack, TParticle *particle, Double_t etaMax, Double_t rMax, Double_t zMax){
  // MonteCarlo Photon Selection
  if(!MCStack)return kFALSE;

  if (particle->GetPdgCode() == 22){
    // check whether particle is within eta range
    if( TMath::Abs(particle->Eta()) > etaMax ) return kFALSE;
    // check if particle doesn't have a photon as mother
    if(particle->GetMother(0) >-1 && MCStack->Particle(particle->GetMother(0))->GetPdgCode() == 22){
      return kFALSE; // no photon as mothers!
    }
    // looking for conversion gammas (electron + positron from pairbuilding (= 5) )
    TParticle* ePos = NULL;
    TParticle* eNeg = NULL;
    if(particle->GetNDaughters() >= 2){
      for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
        if(daughterIndex<0) continue;
        TParticle *tmpDaughter = MCStack->Particle(daughterIndex);
        if(tmpDaughter->GetUniqueID() == 5){
          if(tmpDaughter->GetPdgCode() == 11){
            eNeg = tmpDaughter;
          } else if(tmpDaughter->GetPdgCode() == -11){
            ePos = tmpDaughter;
          }
        }
      }
    }
    if(ePos == NULL || eNeg == NULL){ // means we do not have two daughters from pair production
      return kFALSE;
    }
    // check if electrons are in correct eta window
    if( TMath::Abs(ePos->Eta()) > etaMax ||
      TMath::Abs(eNeg->Eta()) > etaMax )
      return kFALSE;

    // check if photons have converted in reconstructable range
    if(ePos->R() > rMax){
      return kFALSE; // cuts on distance from collision point
    }
    if(TMath::Abs(ePos->Vz()) > zMax){
      return kFALSE;  // outside material
    }
    if(TMath::Abs(eNeg->Vz()) > zMax){
      return kFALSE;  // outside material
    }


    Double_t lineCutZRSlope = tan(2*atan(exp(-etaMax)));
    Double_t lineCutZValue = 7.;
    if( ePos->R() <= ((TMath::Abs(ePos->Vz()) * lineCutZRSlope) - lineCutZValue)){
      return kFALSE;  // line cut to exclude regions where we do not reconstruct
    }
    if( eNeg->R() <= ((TMath::Abs(eNeg->Vz()) * lineCutZRSlope) - lineCutZValue)){
      return kFALSE; // line cut to exclude regions where we do not reconstruct
    }
    if (ePos->Pt() < 0.05 || eNeg->Pt() < 0.05){
      return kFALSE;
    }

    return kTRUE;
  }
  return kFALSE;
}

///_______________________________________________________________________
void AliV0ReaderV1::CreatePureMCHistosForV0FinderEffiESD(){
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();
   //cout << mcProdVtxX <<"\t" << mcProdVtxY << "\t" << mcProdVtxZ << endl;

  AliStack *fMCStack= fMCEvent->Stack();
  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCStack->GetNtrack(); i++) {
    if (fEventCuts->IsConversionPrimaryESD( fMCStack, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
      // fill primary histogram
      TParticle* particle = (TParticle *)fMCStack->Particle(i);
      if (!particle) continue;
      if (ParticleIsConvertedPhoton(fMCStack, particle, 0.9, 180.,250. )){
        if(particle->GetFirstDaughter()<0) continue;
        TParticle *tmpDaughter = fMCStack->Particle(particle->GetFirstDaughter());
        if (!tmpDaughter) continue;
        fHistoMCGammaPtvsR->Fill(particle->Pt(),tmpDaughter->R());
        fHistoMCGammaPtvsPhi->Fill(particle->Pt(),particle->Phi());
        fHistoMCGammaRvsPhi->Fill(tmpDaughter->R(),particle->Phi());
      }
      if (ParticleIsConvertedPhoton(fMCStack, particle, 1.4, 180.,250. )){
        if(particle->GetFirstDaughter()<0) continue;
        TParticle *tmpDaughter = fMCStack->Particle(particle->GetFirstDaughter());
        if (!tmpDaughter) continue;
        fHistoMCGammaPtvsEta->Fill(particle->Pt(),particle->Eta());
        fHistoMCGammaRvsEta->Fill(tmpDaughter->R(),particle->Eta());
        fHistoMCGammaPhivsEta->Fill(particle->Phi(),particle->Eta());
      }
    }
  }
}

///_______________________________________________________________________
void AliV0ReaderV1::FillRecMCHistosForV0FinderEffiESD( AliESDv0* currentV0){

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();
//   cout << mcProdVtxX <<"\t" << mcProdVtxY << "\t" << mcProdVtxZ << endl;

  Int_t tracklabelPos=currentV0->GetPindex();
  Int_t tracklabelNeg=currentV0->GetNindex();

  AliStack *fMCStack= fMCEvent->Stack();

  Int_t labelp=TMath::Abs(fConversionCuts->GetTrack(fInputEvent,tracklabelPos)->GetLabel());
  Int_t labeln=TMath::Abs(fConversionCuts->GetTrack(fInputEvent,tracklabelNeg)->GetLabel());

  TParticle* negPart = 0x0;
  if(labeln>-1) negPart = (TParticle *)fMCStack->Particle(labeln);
  TParticle* posPart = 0x0;
  if(labelp>-1) posPart = (TParticle *)fMCStack->Particle(labelp);

  if ( negPart == NULL || posPart == NULL ) return;
//   if (!(negPart->GetPdgCode() == 11)) return;
//   if (!(posPart->GetPdgCode() == -11)) return;
  Long_t motherlabelNeg = negPart->GetFirstMother();
  Long_t motherlabelPos = posPart->GetFirstMother();

//   cout << "mother neg " << motherlabelNeg << " mother pos " << motherlabelPos << endl;
  if (motherlabelNeg>-1 && motherlabelNeg == motherlabelPos && negPart->GetFirstMother() != -1){
    if (fEventCuts->IsConversionPrimaryESD( fMCStack, negPart->GetFirstMother(), mcProdVtxX, mcProdVtxY, mcProdVtxZ)){

      TParticle* mother =  (TParticle *)fMCStack->Particle(motherlabelNeg);
      if (mother->GetPdgCode() == 22 ){
        if (!CheckVectorForDoubleCount(fVectorFoundGammas,motherlabelNeg ) ){
          if (ParticleIsConvertedPhoton(fMCStack, mother, 0.9, 180.,250. )){
            fHistoRecMCGammaPtvsR->Fill(mother->Pt(),negPart->R());
            fHistoRecMCGammaPtvsPhi->Fill(mother->Pt(),mother->Phi());
            fHistoRecMCGammaRvsPhi->Fill(negPart->R(),mother->Phi());
          }
          if (ParticleIsConvertedPhoton(fMCStack, mother, 1.4, 180.,250. )){
            fHistoRecMCGammaPtvsEta->Fill(mother->Pt(),mother->Eta());
            fHistoRecMCGammaRvsEta->Fill(negPart->R(),mother->Eta());
            fHistoRecMCGammaPhivsEta->Fill(mother->Phi(),mother->Eta());
          }
//           cout << "new gamma found" << endl;
        } else {
          if (ParticleIsConvertedPhoton(fMCStack, mother, 0.9, 180.,250. )){
            fHistoRecMCGammaMultiPt->Fill(mother->Pt());
            fHistoRecMCGammaMultiPhi->Fill(mother->Phi());
            fHistoRecMCGammaMultiR->Fill(negPart->R());
          }
          if (ParticleIsConvertedPhoton(fMCStack, mother, 1.4, 180.,250. )){
            fHistoRecMCGammaMultiPtvsEta->Fill(mother->Pt(),mother->Eta());
          }
//           cout << "this one I had already: " << motherlabelNeg << endl << "-----------------------"  << endl;
        }
//         cout << "event gammas: " << endl;
//        for(Int_t iGamma=0; iGamma<fVectorFoundGammas.size(); iGamma++){cout << fVectorFoundGammas.at(iGamma) << ", ";}
      }
    }
  }
}

//_________________________________________________________________________________
void AliV0ReaderV1::FillImpactParamHistograms( AliVTrack* pTrack, AliVTrack* nTrack, AliESDv0 *fCurrentV0, AliKFConversionPhoton *fCurrentMotherKF){

  // values of cuts to be introduced in 2016 to reduce ESD size
  Float_t fZmax = 25; //cm
  Float_t fYmax = 25; //cm
  Double_t kTPCMargin = 1.0; //cm
  Int_t kMaxTPCV0Conflicts = 1; //# conflicting clusters tolerated

  //conversion point
  Double_t convX, convY, convZ;
  fCurrentV0->GetXYZ(convX,convY,convZ);
  Double_t convR = TMath::Sqrt(convX*convX+convY*convY);
  //recalculated conversion point
  Double_t convXrecalc=fCurrentMotherKF->GetConversionX(); 
  Double_t convYrecalc=fCurrentMotherKF->GetConversionY(); 
  Double_t convZrecalc=fCurrentMotherKF->GetConversionZ(); 
  Double_t convRrecalc = fCurrentMotherKF->GetConversionRadius();

  //count V0s
  fHistoImpactParameterStudy->Fill(0);

  //count V0s with two TPC-only tracks
  AliESDtrack* positiveTrack = (AliESDtrack*) pTrack;
  AliESDtrack* negativeTrack = (AliESDtrack*) nTrack;
  Bool_t TPConly = kFALSE;
  if (!positiveTrack->IsOn(AliESDtrack::kITSin) && !negativeTrack->IsOn(AliESDtrack::kITSin)){
    fHistoImpactParameterStudy->Fill(1);
    TPConly = kTRUE;
  }

  Bool_t RemovedByZcut = kFALSE;
  Bool_t RemovedByYcut = kFALSE;

  //count V0s which have...
  if(TPConly){
    //not passed z cut:    pos.tr.     or           neg.tr. 
    if(((TMath::Abs(positiveTrack->GetZ()))>fZmax) || ((TMath::Abs(negativeTrack->GetZ()))>fZmax)){
      fHistoImpactParameterStudy->Fill(2);
      RemovedByZcut=kTRUE;  
    }

    //not passed y cut:     pos.tr.     or           neg.tr.  
    if(((TMath::Abs(positiveTrack->GetY()))>fYmax) || ((TMath::Abs(negativeTrack->GetY()))>fYmax)){
      fHistoImpactParameterStudy->Fill(3);
      RemovedByYcut=kTRUE;  
    }
  }

  //causality cut
  Bool_t RemovedByCausality=kFALSE;
  AliESDEvent *fESDEvent=dynamic_cast<AliESDEvent*>(fInputEvent);
  const Float_t rTPC[159]={
   85.225, 85.975, 86.725, 87.475, 88.225, 88.975, 89.725, 90.475, 91.225, 91.975, 92.725, 93.475, 94.225, 94.975, 95.725,
   96.475, 97.225, 97.975, 98.725, 99.475,100.225,100.975,101.725,102.475,103.225,103.975,104.725,105.475,106.225,106.975,
   107.725,108.475,109.225,109.975,110.725,111.475,112.225,112.975,113.725,114.475,115.225,115.975,116.725,117.475,118.225,
   118.975,119.725,120.475,121.225,121.975,122.725,123.475,124.225,124.975,125.725,126.475,127.225,127.975,128.725,129.475,
   130.225,130.975,131.725,135.100,136.100,137.100,138.100,139.100,140.100,141.100,142.100,143.100,144.100,145.100,146.100,
   147.100,148.100,149.100,150.100,151.100,152.100,153.100,154.100,155.100,156.100,157.100,158.100,159.100,160.100,161.100,
   162.100,163.100,164.100,165.100,166.100,167.100,168.100,169.100,170.100,171.100,172.100,173.100,174.100,175.100,176.100,
   177.100,178.100,179.100,180.100,181.100,182.100,183.100,184.100,185.100,186.100,187.100,188.100,189.100,190.100,191.100,
   192.100,193.100,194.100,195.100,196.100,197.100,198.100,199.350,200.850,202.350,203.850,205.350,206.850,208.350,209.850,
   211.350,212.850,214.350,215.850,217.350,218.850,220.350,221.850,223.350,224.850,226.350,227.850,229.350,230.850,232.350,
   233.850,235.350,236.850,238.350,239.850,241.350,242.850,244.350,245.850};

  // fill conversion radius histograms
    fHistoR->Fill(convR);
    fHistoRrecalc->Fill(convRrecalc);
    Double_t alpha = TMath::ATan2(convY,convX);  
    if (alpha<0) alpha += TMath::Pi()*2;
    Int_t sec = alpha/(TMath::Pi()/9);
    alpha = (10.+sec*20.)*TMath::DegToRad();
    Double_t cs = TMath::Cos(alpha);
    Double_t sn = TMath::Sin(alpha);
    Double_t xsV0 = convX*cs - convY*sn; 
    fHistoRviaAlpha->Fill(xsV0);
    Double_t alpha_r = TMath::ATan2(convYrecalc,convXrecalc);  
    if (alpha_r<0) alpha_r += TMath::Pi()*2;
    Int_t sec_r = alpha_r/(TMath::Pi()/9);
    alpha_r = (10.+sec_r*20.)*TMath::DegToRad();
    Double_t cs_r = TMath::Cos(alpha_r);
    Double_t sn_r = TMath::Sin(alpha_r);
    Double_t xsV0_r = convXrecalc*cs_r - convYrecalc*sn_r;     
    fHistoRviaAlphaRecalc->Fill(xsV0_r);

  if (convR > 80) {  // conversion within TPC <-> TPC-only tracks
    fHistoImpactParameterStudy->Fill(4);

    for (Int_t it=2;it--;) {
      Int_t trId = fCurrentV0->GetIndex(it);
      AliESDtrack* tr = fESDEvent->GetTrack(trId);
      const TBits& bits = tr->GetTPCClusterMap();
      Int_t nConflict = 0;
      for (Int_t ic=0;ic<159;ic++) {
        if (rTPC[ic]>(xsV0-kTPCMargin)) break;  
        if (bits.TestBitNumber(ic)){
          nConflict++;
          fHistoRdiff->Fill(xsV0-rTPC[ic]);
        }
        if (nConflict>kMaxTPCV0Conflicts) {
          fHistoImpactParameterStudy->Fill(5);
          RemovedByCausality=kTRUE;
          break;
        }
      }
    }
  }

  //not passed y or z o causality cut:
 Bool_t RemovedByAnyCut=kFALSE;
  if(RemovedByZcut||RemovedByYcut||RemovedByCausality){
      fHistoImpactParameterStudy->Fill(6);
      RemovedByAnyCut=kTRUE;
  }

  //Fill tree for further analysis
  Float_t posZ;
  Float_t posY;
  Float_t posX;
  Float_t posPt;
  Float_t negZ;
  Float_t negY;
  Float_t negX;
  Float_t negPt;
  Float_t R;
  TBranch *Branch_Pt = fImpactParamTree->Branch("posPt",&posPt,"posPt/F");
  TBranch *Branch_Y = fImpactParamTree->Branch("posY",&posY,"posY/F");
  TBranch *Branch_R = fImpactParamTree->Branch("R",&R,"R/F");
  posZ = positiveTrack->GetZ();
  posY = positiveTrack->GetY();
  posX = positiveTrack->GetX();
  posPt = positiveTrack->Pt();
  negZ = negativeTrack->GetZ();
  negY = negativeTrack->GetY();
  negX = negativeTrack->GetX();
  negPt = negativeTrack->Pt();
  R=convRrecalc;
  fImpactParamTree->Fill();

  // fill impact parameter histograms
  fHistoPosTrackImpactParamZ->Fill(posZ); 
  fHistoPosTrackImpactParamY->Fill(posY);
  fHistoPosTrackImpactParamX->Fill(posX);
  fHistoPosTrackImpactParamZvsPt->Fill(posPt, posZ); 
  fHistoPosTrackImpactParamYvsPt->Fill(posPt, posY);
  fHistoPosTrackImpactParamXvsPt->Fill(posPt, posX); 
  fHistoNegTrackImpactParamZ->Fill(negZ);  
  fHistoNegTrackImpactParamY->Fill(negY);
  fHistoNegTrackImpactParamX->Fill(negX);
  fHistoNegTrackImpactParamZvsPt->Fill(negPt, negZ);  
  fHistoNegTrackImpactParamYvsPt->Fill(negPt, negY);
  fHistoNegTrackImpactParamXvsPt->Fill(negPt, negX);

  // fill comparison histos before / after new cuts
  fHistoPt->Fill(positiveTrack->Pt());
  fHistoImpactParamZvsR->Fill(convZrecalc,convRrecalc); 
  //Float_t *DCAzPhoton;
  Float_t DCAzPhoton;
  //const AliVVertex *PrimVertex = fInputEvent->GetPrimaryVertex();
  AliAODConversionPhoton* fCurrentMother=(AliAODConversionPhoton*)fCurrentMotherKF;
  //AliAODConversionPhoton *fCurrentMother=dynamic_cast<AliAODConversionPhoton*>(fCurrentMotherKF);
  //fCurrentMotherKF->GetDistanceOfClossetApproachToPrimVtx(PrimVertex, DCAzPhoton);
  DCAzPhoton = fCurrentMother->GetDCAzToPrimVtx();
  fHistoDCAzPhoton->Fill(DCAzPhoton);
  if(!RemovedByAnyCut) {
    fHistoPt2->Fill(positiveTrack->Pt());
    fHistoImpactParamZvsR2->Fill(convZrecalc,convRrecalc);  
    fHistoDCAzPhoton2->Fill(DCAzPhoton);
  }

  return;
}

//_________________________________________________________________________________
Bool_t AliV0ReaderV1::CheckVectorOnly(vector<Int_t> &vec, Int_t tobechecked)
{
  if(tobechecked > -1)
  {
    vector<Int_t>::iterator it;
    it = find (vec.begin(), vec.end(), tobechecked);
    if (it != vec.end()) return true;
    else return false;
  }
  return false;
}

//_________________________________________________________________________________
Bool_t AliV0ReaderV1::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked)
{
  if(tobechecked > -1)
  {
    vector<Int_t>::iterator it;
    it = find (vec.begin(), vec.end(), tobechecked);
    if (it != vec.end()) return true;
    else{
      vec.push_back(tobechecked);
      return false;
    }
  }
  return false;
}

//________________________________________________________________________
void AliV0ReaderV1::Terminate(Option_t *)
{

}
