/*************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
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

//#include <iostream>
//
#include "AliESDVZERO.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliVTrack.h"

#include "TDatabasePDG.h"
#include "TFile.h"
#include "TParticlePDG.h"
#include "TMath.h"
#include "TVectorDfwd.h"
#include "Math/Vector4D.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMultiplicity.h"
#include "AliDataFile.h"
#include "AliExternalTrackParam.h"

#include "AliAODEvent.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliAODVertex.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"

#include "AliESDZDC.h"
#include "AliPIDResponse.h"
#include "AliESDVertex.h"

#include "AliAnalysisTaskUpc4Prongs.h"
//
// TODO: create startup item based on runAnalysis.C macro or even better make from runAnalysis.C startup item
// TODO: add split for trigger string via ';'
// TODO: dynamic switch from esd to aod based on input data
// TODO: add tests - comparison of results from local running on ESD data

ClassImp(AliAnalysisTaskUpc4Prongs);

AliAnalysisTaskUpc4Prongs::AliAnalysisTaskUpc4Prongs()
    : AliAnalysisTaskSE(), fPIDResponse(0), fTriggerName(0), fRhoTree(0),
    BunchCrossNumber(0), OrbitNumber(0), PeriodNumber(0), RunNum(0), Mass(0),
    Q(-10), Pt(0), Rapidity(0), V0Adecision(0), V0Cdecision(0),
    ADAdecision(0), ADCdecision(0), UBAfired(0), UBCfired(0), VBAfired(0),
    VBCfired(0), ZNAenergy(0), ZNCenergy(0), ZPAenergy(0), ZPCenergy(0),
    VtxContrib(0), SpdVtxContrib(0), VtxChi2(0), VtxNDF(0), nTracklets(0), nTracks(0),
    Phi(0), T_NumberOfSigmaTPCPion(0),
    T_NumberOfSigmaTPCElectron(0), T_NumberOfSigmaITSPion(0), T_NumberOfSigmaITSElectron(0),
    T_TPCsignal(0), T_P(0), T_Eta(0), T_Phi(0), T_Px(0), T_Py(0),
    T_Pz(0), T_Q(0), T_HasPointOnITSLayer0(0), T_HasPointOnITSLayer1(0), 
    T_ITSModuleInner(0), T_ITSModuleOuter(0), T_TPCNCls(0), T_ITSNCls(0), T_Dca0(0),
    T_Dca1(0), T_TPCRefit(0), T_ITSRefit(0), T_Lets_Theta(0), T_Lets_Phi(0), T_ITSSensorNum(0), T_ITSsa(0) {}

AliAnalysisTaskUpc4Prongs::AliAnalysisTaskUpc4Prongs(const char* name)
    : AliAnalysisTaskSE(name), fPIDResponse(0), fTriggerName(0), fRhoTree(0),
    BunchCrossNumber(0), OrbitNumber(0), PeriodNumber(0), RunNum(0), Mass(0),
    Q(-10), Pt(0), Rapidity(0), V0Adecision(0), V0Cdecision(0),
    ADAdecision(0), ADCdecision(0), UBAfired(0), UBCfired(0), VBAfired(0),
    VBCfired(0), ZNAenergy(0), ZNCenergy(0), ZPAenergy(0), ZPCenergy(0),
    VtxContrib(0), SpdVtxContrib(0), VtxChi2(0), VtxNDF(0), nTracklets(0), nTracks(0),
    Phi(0), T_NumberOfSigmaTPCPion(0),
    T_NumberOfSigmaTPCElectron(0), T_NumberOfSigmaITSPion(0), T_NumberOfSigmaITSElectron(0), 
    T_TPCsignal(0), T_P(0), T_Eta(0), T_Phi(0), T_Px(0), T_Py(0),
    T_Pz(0), T_Q(0), T_HasPointOnITSLayer0(0), T_HasPointOnITSLayer1(0),
    T_ITSModuleInner(0), T_ITSModuleOuter(0), T_TPCNCls(0), T_ITSNCls(0), T_Dca0(0),
    T_Dca1(0), T_TPCRefit(0), T_ITSRefit(0), T_Lets_Theta(0), T_Lets_Phi(0), T_ITSSensorNum(0), T_ITSsa(0) {
    Init();
    DefineOutput(1, TTree::Class());
    DefineOutput(2, TList::Class());
}

AliAnalysisTaskUpc4Prongs::~AliAnalysisTaskUpc4Prongs()
{
    if (fRhoTree)     { delete fRhoTree;     fRhoTree     = nullptr; }
    if (fPIDResponse) { delete fPIDResponse; fPIDResponse = nullptr; }
}

void AliAnalysisTaskUpc4Prongs::Init()
{
    for (Int_t i = 0; i < 3; i++)
    {
        Vertex[i]    = -1717;
        SpdVertex[i] = -1717;
        ZDCAtime[i]  = -1717.;
        ZDCCtime[i]  = -1717.;
    }

    ZDCAtime[3] = -1717.;
    ZDCCtime[3] = -1717.;

  /*  PIDTPCPion.reserve(200);
    PIDTPCElectron.reserve(200);
    TPCsignal.reserve(200);
    TrackP.reserve(200);
    TrackEta.reserve(200);
    TrackPhi.reserve(200);
    TrackPx.reserve(200);
    TrackPy.reserve(200);
    TrackPz.reserve(200);
    TrackQ.reserve(200);
    TrackITSModuleInner.reserve(200);
    TrackITSModuleOuter.reserve(200);
    TrackHasPointOnITSLayer0.reserve(200);
    TrackHasPointOnITSLayer1.reserve(200);*/
}

void AliAnalysisTaskUpc4Prongs::UserCreateOutputObjects()
{
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    fRhoTree = new TTree("events", "Selected events for 4proungs analysis");

    fRhoTree->Branch("RunNum",                                   &RunNum,                   "RunNum/I");
    fRhoTree->Branch("PeriodNumber",                             &PeriodNumber,             "PeriodNumber/i");
    fRhoTree->Branch("OrbitNumber",                              &OrbitNumber,              "OrbitNumber/i");
    fRhoTree->Branch("BunchCrossNumber",                         &BunchCrossNumber,         "BunchCrossNumber/s");
    fRhoTree->Branch("Mass",                                     &Mass,                     "Mass/F");
    fRhoTree->Branch("Pt",                                       &Pt,                       "Pt/F");
    fRhoTree->Branch("Q",                                        &Q,                        "Q/S");
    fRhoTree->Branch("Rapidity",                                 &Rapidity,                 "Rapidity/F");
    fRhoTree->Branch("Phi",                                      &Phi,                      "Phi/F");
    fRhoTree->Branch("ZNAenergy",                                &ZNAenergy,                "ZNAenergy/F");
    fRhoTree->Branch("ZNCenergy",                                &ZNCenergy,                "ZNCenergy/F");
    fRhoTree->Branch("ZPAenergy",                                &ZPAenergy,                "ZPAenergy/F");
    fRhoTree->Branch("ZPCenergy",                                &ZPCenergy,                "ZPCenergy/F");
    fRhoTree->Branch("VtxX",                                     &Vertex[0],                "VtxX/F");
    fRhoTree->Branch("VtxY",                                     &Vertex[1],                "VtxY/F");
    fRhoTree->Branch("VtxZ",                                     &Vertex[2],                "VtxZ/F");
    fRhoTree->Branch("VtxContrib",                               &VtxContrib,               "VtxContrib/I");
    fRhoTree->Branch("VtxChi2",                                  &VtxChi2,                  "VtxChi2/F");
    fRhoTree->Branch("VtxNDF",                                   &VtxNDF,                   "VtxNDF/F");
    fRhoTree->Branch("SpdVtxX",                                  &SpdVertex[0],             "SpdVtxX/F");
    fRhoTree->Branch("SpdVtxY",                                  &SpdVertex[1],             "SpdVtxY/F");
    fRhoTree->Branch("SpdVtxZ",                                  &SpdVertex[2],             "SpdVtxZ/F");
    fRhoTree->Branch("SpdVtxContrib",                            &SpdVtxContrib,            "SpdVtxContrib/I");
    fRhoTree->Branch("V0Adecision",                              &V0Adecision,              "V0Adecision/I");
    fRhoTree->Branch("V0Cdecision",                              &V0Cdecision,              "V0Cdecision/I");
    fRhoTree->Branch("ADAdecision",                              &ADAdecision,              "ADAdecision/I");
    fRhoTree->Branch("ADCdecision",                              &ADCdecision,              "ADCdecision/I");
    fRhoTree->Branch("UBAfired",                                 &UBAfired,                 "UBAfired/O");
    fRhoTree->Branch("UBCfired",                                 &UBCfired,                 "UBCfired/O");
    fRhoTree->Branch("VBAfired",                                 &VBAfired,                 "VBAfired/O");
    fRhoTree->Branch("VBCfired",                                 &VBCfired,                 "VBCfired/O");
    fRhoTree->Branch("nTracklets",                               &nTracklets,               "nTracklets/I");
    fRhoTree->Branch("nTracks",                                  &nTracks,                  "nTracks/I");
    //fRhoTree->Branch("ZDCAtime",                               &ZDCAtime,                 "ZDCAtime[4]/F");
    //fRhoTree->Branch("ZDCCtime",                               &ZDCCtime,                 "ZDCCtime[4]/F");
    fRhoTree->Branch("T_NumberOfSigmaITSPion",                   &T_NumberOfSigmaITSPion); //,               "PIDTPCPion[4]/F");
    fRhoTree->Branch("T_NumberOfSigmaITSElectron",               &T_NumberOfSigmaITSElectron); //,           "PIDTPCElectron[4]/F");
    fRhoTree->Branch("T_NumberOfSigmaTPCPion",                   &T_NumberOfSigmaTPCPion); //,               "PIDTPCPion[4]/F");
    fRhoTree->Branch("T_NumberOfSigmaTPCElectron",               &T_NumberOfSigmaTPCElectron); //,           "PIDTPCElectron[4]/F");
    fRhoTree->Branch("TPCsignal",                                &T_TPCsignal); //,                "TPCsignal[4]/I");
    fRhoTree->Branch("T_P",                                      &T_P); //,                   "T_P[4]/F");
    fRhoTree->Branch("T_Eta",                                    &T_Eta); //,                 "T_Eta[4]/F");
    fRhoTree->Branch("T_Phi",                                    &T_Phi); //,                 "T_Phi[4]/F");
    fRhoTree->Branch("T_Px",                                     &T_Px); //,                  "T_Px[4]/F");
    fRhoTree->Branch("T_Py",                                     &T_Py); //,                  "T_Py[4]/F");
    fRhoTree->Branch("T_Pz",                                     &T_Pz); //,                  "T_Pz[4]/F");
    fRhoTree->Branch("T_Q",                                      &T_Q); //,                   "T_Q[4]/S");
    fRhoTree->Branch("T_HasPointOnITSLayer0",                    &T_HasPointOnITSLayer0); //, "T_HasPointOnITSLayer0[4]/O");
    fRhoTree->Branch("T_HasPointOnITSLayer1",                    &T_HasPointOnITSLayer1); //, "T_HasPointOnITSLayer1[4]/O");
    fRhoTree->Branch("T_ITSModuleInner",                         &T_ITSModuleInner); //,      "T_ITSModuleInner[4]/I");
    fRhoTree->Branch("T_ITSModuleOuter",                         &T_ITSModuleOuter); //,      "T_ITSModuleOuter[4]/I");
    fRhoTree->Branch("T_TPCNCls",                                &T_TPCNCls);
    fRhoTree->Branch("T_ITSNCls",                                &T_ITSNCls);
    fRhoTree->Branch("T_Dca0",                                   &T_Dca0);
    fRhoTree->Branch("T_Dca1",                                   &T_Dca1);
    fRhoTree->Branch("T_TPCRefit",                               &T_TPCRefit);
    fRhoTree->Branch("T_ITSRefit",                               &T_ITSRefit);
    fRhoTree->Branch("T_ITSsa",                                  &T_ITSsa);
    fRhoTree->Branch("TLets_Theta",                              &T_Lets_Theta);
    fRhoTree->Branch("TLets_Phi",                                &T_Lets_Phi);
    fRhoTree->Branch("T_ITSSensorNum",                           &T_ITSSensorNum);

    //fRhoTree->Branch("ITSModule",                &ITSModule,                "ITSModule/I");

    PostData(1, fRhoTree);
}

void AliAnalysisTaskUpc4Prongs::UserExec(Option_t*)
{
    AliESDEvent* esd = (AliESDEvent*)InputEvent();
    if (!esd) return;

    TString trigger = esd->GetFiredTriggerClasses();
    if (!trigger.Contains("CCUP9-B")) return;

    T_NumberOfSigmaITSPion.clear();
    T_NumberOfSigmaITSElectron.clear();
    T_NumberOfSigmaTPCPion.clear();
    T_NumberOfSigmaTPCElectron.clear();
    T_TPCsignal.clear();
    T_TPCNCls.clear();
    T_ITSNCls.clear();
    T_P.clear();
    T_Eta.clear();
    T_Phi.clear();
    T_Px.clear();
    T_Py.clear();
    T_Pz.clear();
    T_Dca0.clear();
    T_Dca1.clear();
    T_Q.clear();
    T_TPCRefit.clear();
    T_ITSRefit.clear();
    T_HasPointOnITSLayer0.clear();
    T_HasPointOnITSLayer1.clear();
    T_ITSModuleInner.clear();
    T_ITSModuleOuter.clear();
    T_Lets_Theta.clear();
    T_Lets_Phi.clear();
    T_ITSSensorNum.clear();
    T_ITSsa.clear();

    nTracks = 0;
    Q = 0;

    AliESDVertex* fESDVertex = (AliESDVertex*)esd->GetPrimaryVertex();
    TDatabasePDG* pdgdat = TDatabasePDG::Instance();
    TParticlePDG* partPion = pdgdat->GetParticle(211);
    Double_t pionMass = partPion->Mass();
  
    // event info
    RunNum = esd->GetRunNumber();
    OrbitNumber = esd->GetOrbitNumber();
    PeriodNumber = esd->GetPeriodNumber();
    BunchCrossNumber = esd->GetBunchCrossNumber();

    // VZERO, ZDC, AD
    AliESDVZERO* fV0data = esd->GetVZEROData();
    AliESDZDC* fZDCdata = esd->GetESDZDC();
    AliESDAD* fADdata = esd->GetADData();

    V0Adecision = fV0data->GetV0ADecision();
    V0Cdecision = fV0data->GetV0CDecision();

    if (fADdata)
    {
        ADAdecision = fADdata->GetADADecision();
        ADCdecision = fADdata->GetADCDecision();
    }

    UBAfired = esd->GetHeader()->IsTriggerInputFired("0UBA");
    UBCfired = esd->GetHeader()->IsTriggerInputFired("0UBC");
    VBAfired = esd->GetHeader()->IsTriggerInputFired("0VBA");
    VBCfired = esd->GetHeader()->IsTriggerInputFired("0VBC");

    // ZN energy
    ZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
    ZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];
    ZPAenergy = fZDCdata->GetZPATowerEnergy()[0];
    ZPCenergy = fZDCdata->GetZPCTowerEnergy()[0];

    // neutron ZDC time
    Int_t detChZNA = fZDCdata->GetZNATDCChannel();
    Int_t detChZNC = fZDCdata->GetZNCTDCChannel();

   /* for (Int_t i = 0; i < 4; i++)
    {
        ZDCAtime[i] = fZDCdata->GetZDCTDCCorrected(detChZNA, i);
        ZDCCtime[i] = fZDCdata->GetZDCTDCCorrected(detChZNC, i);
    }*/

    // primary vertex
    VtxContrib = fESDVertex->GetNContributors();
    Vertex[0] = fESDVertex->GetX();
    Vertex[1] = fESDVertex->GetY();
    Vertex[2] = fESDVertex->GetZ();
    VtxChi2 = fESDVertex->GetChi2();
    VtxNDF = fESDVertex->GetNDF();

    // SPD primary vertex
    AliESDVertex* fSPDVertex = (AliESDVertex*)esd->GetPrimaryVertexSPD();
    SpdVtxContrib = fSPDVertex->GetNContributors();
    SpdVertex[0] = fSPDVertex->GetX();
    SpdVertex[1] = fSPDVertex->GetY();
    SpdVertex[2] = fSPDVertex->GetZ();

    // Tracklets
    nTracklets = esd->GetMultiplicity()->GetNumberOfTracklets();

    for (auto i = 0; i < nTracklets; ++i)
    {
        T_Lets_Theta.push_back(esd->GetMultiplicity()->GetTheta(i));
        T_Lets_Phi.push_back(esd->GetMultiplicity()->GetPhi(i));
    }

    std::vector<ROOT::Math::PxPyPzMVector> EventVectors;
    //TrackLV.reserve(200);

    // Track loop - cuts
    for (Int_t i = 0; i < esd->GetNumberOfTracks(); ++i)
    {
        AliESDtrack* trk = esd->GetTrack(i);
        if (!trk) continue;
        if (!trk->HasPointOnITSLayer(0) && !trk->HasPointOnITSLayer(1)) continue;

        Float_t dca[2] = { 0.0,0.0 }; AliExternalTrackParam cParam;
        if (!trk->RelateToVertex(fESDVertex, esd->GetMagneticField(), 300., &cParam)) continue;
        trk->GetImpactParameters(dca[0], dca[1]);
        if (TMath::Abs(dca[1]) > 2) continue;
        Double_t cut_DCAxy = (0.0182 + 0.0350 / TMath::Power(trk->Pt(), 1.01));
        if (TMath::Abs(dca[0]) > cut_DCAxy) continue;
        if (i >= 200) return;

        nTracks++;
        Q += trk->Charge();

        T_ITSModuleInner.push_back(trk->GetITSModuleIndex(0));
        T_ITSModuleOuter.push_back(trk->GetITSModuleIndex(1));

        T_HasPointOnITSLayer0.push_back(trk->HasPointOnITSLayer(0));
        T_HasPointOnITSLayer1.push_back(trk->HasPointOnITSLayer(1));

        T_ITSNCls.push_back(trk->GetNumberOfITSClusters());
        T_TPCNCls.push_back(trk->GetNumberOfTPCClusters());
        T_TPCRefit.push_back(trk->IsOn(AliESDtrack::kTPCrefit));
        T_ITSRefit.push_back(trk->IsOn(AliESDtrack::kITSrefit));
        T_ITSsa.push_back(trk->IsPureITSStandalone());

        // TPC&ITS PID n-sigma
        T_NumberOfSigmaTPCElectron.push_back(fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kElectron));
        T_NumberOfSigmaTPCPion.push_back(fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kPion));
        T_NumberOfSigmaITSPion.push_back(fPIDResponse->NumberOfSigmasITS(trk, AliPID::kPion));
        T_NumberOfSigmaITSElectron.push_back(fPIDResponse->NumberOfSigmasITS(trk, AliPID::kElectron));

        T_Q.push_back(trk->Charge());
        T_TPCsignal.push_back(trk->GetTPCsignal());
        T_P.push_back(trk->P());
        T_Phi.push_back(trk->Phi());
        T_Eta.push_back(trk->Eta());
        T_Px.push_back(trk->Px());
        T_Py.push_back(trk->Py());
        T_Pz.push_back(trk->Pz());

        T_Dca0.push_back(dca[0]);
        T_Dca1.push_back(dca[1]);

        ROOT::Math::PxPyPzMVector trackV(trk->Px(), trk->Py(), trk->Pz(), pionMass);
        EventVectors.push_back(trackV);

    } // end loop over tracks

    if (nTracks < 4) return;

    for (Int_t chipkey = 0; chipkey < 1200; chipkey++) {
        if (esd->GetMultiplicity()->TestFastOrFiredChips(chipkey)) 
            T_ITSSensorNum.push_back(chipkey / 5);
    }

  /*  PIDTPCPion.shrink_to_fit();
    PIDTPCElectron.shrink_to_fit();
    TPCsignal.shrink_to_fit();
    TrackP.shrink_to_fit();
    TrackEta.shrink_to_fit();
    TrackPhi.shrink_to_fit();
    TrackPx.shrink_to_fit();
    TrackPy.shrink_to_fit();
    TrackPz.shrink_to_fit();
    TrackQ.shrink_to_fit();
    TrackITSModuleInner.shrink_to_fit();
    TrackITSModuleOuter.shrink_to_fit();
    TrackHasPointOnITSLayer0.shrink_to_fit();
    TrackHasPointOnITSLayer1.shrink_to_fit();*/
    
    ROOT::Math::PxPyPzMVector sumVector;
    for (const auto& tv : EventVectors)
        sumVector += tv;
        
    Mass = sumVector.M();
    Pt = sumVector.Pt();
    Rapidity = sumVector.Rapidity();
    Phi = sumVector.Phi();

    fRhoTree->Fill();

    PostData(1, fRhoTree);

} // UserExec

