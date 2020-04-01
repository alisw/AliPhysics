/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************************
// \class AliAnalysisTaskSECharmHadronvnTMVA
// \brief task for the analysis of D-meson vn
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
// F. Catalano, fabio.catalano@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// G. Luparello, grazia.luparello@cern.ch
// F. Prino, prino@to.infn.it
// A. Rossi, andrea.rossi@cern.ch
// S. Trogolo, stefano.trogolo@cern.ch
///////////////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TFile.h>
#include <TDatabasePDG.h>
#include <TAxis.h>
#include <TSpline.h>
#include <TGrid.h>

#include <AliLog.h>
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAnalysisVertexingHF.h"

#include "AliMultSelection.h"
#include "AliVertexingHFUtils.h"

#include "AliAnalysisTaskSEHFTenderQnVectors.h"
#include "AliAnalysisTaskSECharmHadronvnTMVA.h"

#include "AliAODPidHF.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSECharmHadronvnTMVA)
/// \endcond
//________________________________________________________________________
AliAnalysisTaskSECharmHadronvnTMVA::AliAnalysisTaskSECharmHadronvnTMVA() :
AliAnalysisTaskSE(),
    fAOD(nullptr),
    fOutput(nullptr),
    fHistNEvents(nullptr),
    fHistCandVsCent(nullptr),
    fHistPercqnVsqnVsCentr(nullptr),
    fHistNtrklVsqnVsCentr(nullptr),
    fHistMassPtPhiqnCentr(nullptr),
    fHistMassPtPhiqnCentr_1(nullptr),
    fHistMassPtPhiqnCentr_2(nullptr),
    fHistMassPtPhiqnCentr_3(nullptr),
    fHistMassPtPhiqnCentr_4(nullptr),
    fHistMassPtPhiqnCentr_5(nullptr),
    fRDCuts(nullptr),
    fTenderTaskName("HFTenderQnVectors"),
    fMinCentr(0.),
    fMaxCentr(100.),
    fAODProtection(1),
    fDaughterTracks(),
    fHarmonic(2),
    fCalibType(AliHFQnVectorHandler::kQnFrameworkCalib),
    fNormMethod(AliHFQnVectorHandler::kQoverM),
    fOADBFileName(""),
    fFlowMethod(kEvShapeEP),
    fEvPlaneDet(kFullV0),
    fSubEvDetA(kPosTPC),
    fSubEvDetB(kNegTPC),
    fqnMeth(kq2TPC),
    fScalProdLimit(0.4),
    fPercentileqn(false),
    fqnSplineFileName(""),
    fLoadedSplines(false),
    fDecChannel(kD0toKpi),
    fLowmasslimit(1.669),
    fUpmasslimit(2.069),
    fNMassBins(200),
    fEtaGapInTPCHalves(0),
    fRemoveDauFromqn(0),
    fListRDHFBDT(0),
    fListBDTNtuple(0),
    fRemoveSoftPion(false),
    fEnableDownsamplqn(false),
    fFracToKeepDownSamplqn(1.1)
{
    // Default constructor
    for (Int_t ih=0; ih<13; ih++) fBDT1Cut[ih] = -2;
    for (Int_t in=0; in<13; in++) fBDT2Cut[in] = -2;
    for(int iHisto=0; iHisto<3; iHisto++) {
        fHistCentrality[iHisto]         = nullptr;
        fHistEPResolVsCentrVsqn[iHisto] = nullptr;
    }
    for(int iDet=0; iDet<6; iDet++) {
        fHistEvPlaneQncorr[iDet]        = nullptr;
        fHistqnVsCentrPercCalib[iDet]   = nullptr;
        fqnSplinesList[iDet]            = nullptr;
    }
}

//________________________________________________________________________
AliAnalysisTaskSECharmHadronvnTMVA::AliAnalysisTaskSECharmHadronvnTMVA(const char *name, AliRDHFCuts *rdCuts, int decaychannel) :
    AliAnalysisTaskSE(name),
    fAOD(nullptr),
    fOutput(nullptr),
    fHistNEvents(nullptr),
    fHistCandVsCent(nullptr),
    fHistPercqnVsqnVsCentr(nullptr),
    fHistNtrklVsqnVsCentr(nullptr),
    fHistMassPtPhiqnCentr(nullptr),
    fHistMassPtPhiqnCentr_1(nullptr),
    fHistMassPtPhiqnCentr_2(nullptr),
    fHistMassPtPhiqnCentr_3(nullptr),
    fHistMassPtPhiqnCentr_4(nullptr),
    fHistMassPtPhiqnCentr_5(nullptr),
    fRDCuts(rdCuts),
    fTenderTaskName("HFTenderQnVectors"),
    fMinCentr(0.),
    fMaxCentr(100.),
    fAODProtection(1),
    fDaughterTracks(),
    fHarmonic(2),
    fCalibType(AliHFQnVectorHandler::kQnFrameworkCalib),
    fNormMethod(AliHFQnVectorHandler::kQoverM),
    fOADBFileName(""),
    fFlowMethod(kEvShapeEP),
    fEvPlaneDet(kFullV0),
    fSubEvDetA(kPosTPC),
    fSubEvDetB(kNegTPC),
    fqnMeth(kq2TPC),
    fScalProdLimit(0.4),
    fPercentileqn(false),
    fqnSplineFileName(""),
    fLoadedSplines(false),
    fDecChannel(decaychannel),
    fLowmasslimit(1.669),
    fUpmasslimit(2.069),
    fNMassBins(200),
    fEtaGapInTPCHalves(0),
    fRemoveDauFromqn(0),
    fRemoveSoftPion(false),
    fEnableDownsamplqn(false),
    fListRDHFBDT(0),
    fListBDTNtuple(0),
    fFracToKeepDownSamplqn(1.1)
{
    // standard constructor
    for (Int_t ih=0; ih<13; ih++) fBDT1Cut[ih] = -2;
    for (Int_t in=0; in<13; in++) fBDT2Cut[in] = -2;
    for(int iHisto=0; iHisto<3; iHisto++) {
        fHistCentrality[iHisto]         = nullptr;
        fHistEPResolVsCentrVsqn[iHisto] = nullptr;
    }
    for(int iDet=0; iDet<6; iDet++) {
        fHistEvPlaneQncorr[iDet]        = nullptr;
        fHistqnVsCentrPercCalib[iDet]   = nullptr;
        fqnSplinesList[iDet]            = nullptr;
    }

    int pdg=421;
    switch(fDecChannel){
        case kDplustoKpipi:
            pdg=411;
        break;
        case kD0toKpi:
            pdg=421;
        break;
        case kDstartoKpipi:
            pdg=413;
        break;
        case kDstoKKpi:
            pdg=431;
        break;
    }
    if(pdg==413) SetMassLimits((float)0.1,(float)0.2);
    else SetMassLimits((float)0.2,pdg); //check range

    // Output slot #1 writes into a TList container
    DefineOutput(1,TList::Class());  //Main output
    // Output slot #3 writes into a AliRDHFCuts container (cuts)
    switch(fDecChannel){
        case kDplustoKpipi:
            DefineOutput(2,AliRDHFCutsDplustoKpipi::Class());  //Cut object for Dplus
        break;
        case kD0toKpi:
            DefineOutput(2,AliRDHFCutsD0toKpi::Class());       //Cut object for D0
        break;
        case kDstartoKpipi:
            DefineOutput(2,AliRDHFCutsDStartoKpipi::Class());  //Cut object for D*
        break;
        case kDstoKKpi:
            DefineOutput(2,AliRDHFCutsDstoKKpi::Class());      //Cut object for Ds
        break;
    }

}

//________________________________________________________________________
AliAnalysisTaskSECharmHadronvnTMVA::~AliAnalysisTaskSECharmHadronvnTMVA()
{
    // Destructor
    if(fOutput && !fOutput->IsOwner()){
        delete fHistNEvents;
        delete fHistCandVsCent;
        delete fHistPercqnVsqnVsCentr;
        delete fHistNtrklVsqnVsCentr;
        delete fHistMassPtPhiqnCentr;
        delete fHistMassPtPhiqnCentr_1;
        delete fHistMassPtPhiqnCentr_2;
        delete fHistMassPtPhiqnCentr_3;
        delete fHistMassPtPhiqnCentr_4;
        delete fHistMassPtPhiqnCentr_5;

        for(int iHisto=0; iHisto<3; iHisto++){
            delete fHistCentrality[iHisto];
            delete fHistEPResolVsCentrVsqn[iHisto];
        }
        for(int iDet=0; iDet<6; iDet++) {
            delete fHistEvPlaneQncorr[iDet];
            delete fHistqnVsCentrPercCalib[iDet];
        }
    }
    delete fOutput;
    delete fRDCuts;

    for(int iDet=0; iDet<6; iDet++) {
        if(fqnSplinesList[iDet] && fLoadedSplines) delete fqnSplinesList[iDet];
    }
    if (fListRDHFBDT) {
        delete fListRDHFBDT;
        fListRDHFBDT = 0;
    }
    if (fListBDTNtuple) {
        delete fListBDTNtuple;
        fListBDTNtuple = 0;
    }
}
//_________________________________________________________________
void  AliAnalysisTaskSECharmHadronvnTMVA::SetMassLimits(float range, int pdg){
    // Set limits for mass spectra plots
    float mass=0;
    int abspdg=TMath::Abs(pdg);
    mass=TDatabasePDG::Instance()->GetParticle(abspdg)->Mass();
    fUpmasslimit = mass+range;
    fLowmasslimit = mass-range;
}

//_________________________________________________________________
void  AliAnalysisTaskSECharmHadronvnTMVA::SetMassLimits(float lowlimit, float uplimit){
    // Set limits for mass spectra plots
    if(uplimit>lowlimit) {
        fUpmasslimit = uplimit;
        fLowmasslimit = lowlimit;
    }
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronvnTMVA::LocalInit()
{
    // Initialization
    fRDCuts->SetMinCentrality(fMinCentr);
    fRDCuts->SetMaxCentrality(fMaxCentr);

    switch(fDecChannel) {
        case kDplustoKpipi:
            {
                AliRDHFCutsDplustoKpipi* copycut=new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fRDCuts)));
                PostData(2,copycut);
            }
        break;
        case kD0toKpi:
            {
                AliRDHFCutsD0toKpi* copycut=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fRDCuts)));
                PostData(2,copycut);
            }
        break;
        case kDstartoKpipi:
            {
                AliRDHFCutsDStartoKpipi* copycut=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fRDCuts)));
                PostData(2,copycut);
            }
        break;
        case kDstoKKpi:
            {
                AliRDHFCutsDstoKKpi* copycut=new AliRDHFCutsDstoKKpi(*(static_cast<AliRDHFCutsDstoKKpi*>(fRDCuts)));
                PostData(2,copycut);
            }
        break;
    }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronvnTMVA::UserCreateOutputObjects()
{
    // Create the output container
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("MainOutput");
   // fListRDHFBDT->SetOwner();
    fHistNEvents = new TH1F("fHistNEvents", "Number of AODs scanned",16,-0.5,15.5);
    fHistNEvents->GetXaxis()->SetBinLabel(1,"nEventsRead");
    fHistNEvents->GetXaxis()->SetBinLabel(2,"nEvents Matched dAOD");
    fHistNEvents->GetXaxis()->SetBinLabel(3,"nEvents Mismatched dAOD");
    fHistNEvents->GetXaxis()->SetBinLabel(4,"nEventsAnal");
    fHistNEvents->GetXaxis()->SetBinLabel(5,"n. passing IsEvSelected");
    fHistNEvents->GetXaxis()->SetBinLabel(6,"n. rejected due to trigger");
    fHistNEvents->GetXaxis()->SetBinLabel(7,"n. rejected due to not reco vertex");
    fHistNEvents->GetXaxis()->SetBinLabel(8,"n. rejected for contr vertex");
    fHistNEvents->GetXaxis()->SetBinLabel(9,"n. rejected for vertex out of accept");
    fHistNEvents->GetXaxis()->SetBinLabel(10,"n. rejected for pileup events");
    fHistNEvents->GetXaxis()->SetBinLabel(11,Form("no. of out %.0f-%.0f%% centrality events",fMinCentr,fMaxCentr));
    fHistNEvents->GetXaxis()->SetBinLabel(12,"n. rejected for bad cent corr");
    fHistNEvents->GetXaxis()->SetBinLabel(13,"non valid TPC EP");
    fHistNEvents->GetXaxis()->SetBinLabel(14,"bad event plane");
    fHistNEvents->GetXaxis()->SetBinLabel(15,"no. of sel. candidates");
    fHistNEvents->GetXaxis()->SetBinLabel(16,"no. cand. out of pt bounds");
    fHistNEvents->GetXaxis()->SetNdivisions(1,false);
    fOutput->Add(fHistNEvents);

    fHistCentrality[0]=new TH1F("hCentr","centrality",10000,0.,100.);
    fHistCentrality[1]=new TH1F("hCentr(selectedCent)","centrality(selectedCent)",10000,0.,100.);
    fHistCentrality[2]=new TH1F("hCentr(OutofCent)","centrality(OutofCent)",10000,0.,100.);
    for(int iHisto=0; iHisto<3; iHisto++){
        fOutput->Add(fHistCentrality[iHisto]);
    }

    const int ncentbins = static_cast<int>(fMaxCentr-fMinCentr);

    fHistCandVsCent=new TH2F("hCandVsCent","number of selected candidates vs. centrality;centrality(%);number of candidates",ncentbins,fMinCentr,fMaxCentr,101,-0.5,100.5);
    fOutput->Add(fHistCandVsCent);

    TString detConfName[6] = {"TPC","TPCPosEta","TPCNegEta","V0","V0A","V0C"};

    TString qnaxisname=Form("#it{q}_{%d}",fHarmonic);
    switch(fqnMeth) {
        case kq2TPC:
            qnaxisname=Form("#it{q}_{%d}^{TPC}",fHarmonic);
        break;
        case kq2PosTPC:
            qnaxisname=Form("#it{q}_{%d}^{TPCPosEta}",fHarmonic);
        break;
        case kq2NegTPC:
            qnaxisname=Form("#it{q}_{%d}^{TPCNegEta}",fHarmonic);
        break;
        case kq2VZERO:
            qnaxisname=Form("#it{q}_{%d}^{V0}",fHarmonic);
        break;
        case kq2VZEROA:
            qnaxisname=Form("#it{q}_{%d}^{V0A}",fHarmonic);
        break;
        case kq2VZEROC:
            qnaxisname=Form("#it{q}_{%d}^{V0C}",fHarmonic);
        break;
    }

    TString qnpercaxisname = qnaxisname + " (%)";
    TString qnaxisnamefill = qnaxisname;
    if(fPercentileqn)
        qnaxisnamefill = qnpercaxisname;

    int nqnbins=1; //single bin if unbiased analysis
    double qnmin = 0.;
    double qnmax = 15.;
    if(fFlowMethod==kEvShapeEP || fFlowMethod==kEvShapeSP || fFlowMethod==kEvShapeEPVsMass) {
        if(!fPercentileqn) {
            nqnbins=300;
        }
        else {
            nqnbins=100;
            qnmin = 0.;
            qnmax = 100.;
        }
    }

    for(int iDet = 0; iDet < 6; iDet++) {
        fHistEvPlaneQncorr[iDet]   = new TH3F(Form("fHistEvPlaneQncorr%sVsqnVsCent",detConfName[iDet].Data()),Form("hEvPlaneQncorr%s;centrality(%%);%s;#psi_{%d}",detConfName[iDet].Data(),qnaxisnamefill.Data(),fHarmonic),ncentbins,fMinCentr,fMaxCentr,nqnbins,qnmin,qnmax,200,0.,TMath::Pi());
        fOutput->Add(fHistEvPlaneQncorr[iDet]);

        // histos for qn vs. centrality with fine binning (for qn percentiles calibration)
        if(fFlowMethod==kEvShapeEP || fFlowMethod==kEvShapeSP || fFlowMethod==kEvShapeEPVsMass) {
            fHistqnVsCentrPercCalib[iDet] = new TH2F(Form("fHistqnVsCentr%s",detConfName[iDet].Data()),Form("#it{q}_{%d}^{%s} vs. centrality;centrality(%%);#it{q}_{%d}^{%s}",fHarmonic,detConfName[iDet].Data(),fHarmonic,detConfName[iDet].Data()),ncentbins,fMinCentr,fMaxCentr,15000,0.,15.);
            fOutput->Add(fHistqnVsCentrPercCalib[iDet]);
        }
    }

    //qn percentile vs. qn vs. centrality
    if(fFlowMethod==kEvShapeEP || fFlowMethod==kEvShapeSP || fFlowMethod==kEvShapeEPVsMass) {
        fHistNtrklVsqnVsCentr = new TH3F("fHistNtrklVsqnVsCentr",Form("#it{N}_{tracklets} vs. %s vs. centrality;centrality (%%);%s;#it{N}_{tracklets}",qnpercaxisname.Data(),qnpercaxisname.Data()),ncentbins,fMinCentr,fMaxCentr,nqnbins,qnmin,qnmax,500,-0.5,4999.5);
        fOutput->Add(fHistNtrklVsqnVsCentr);
        if(fPercentileqn) {
            fHistPercqnVsqnVsCentr = new TH3F("fHistPercqnVsqnVsCentr",Form("%s vs. %s vs. centrality;centrality (%%);%s;%s",qnpercaxisname.Data(),qnaxisname.Data(),qnaxisname.Data(),qnpercaxisname.Data()),ncentbins,fMinCentr,fMaxCentr,300,0,15,nqnbins,qnmin,qnmax);
            fOutput->Add(fHistPercqnVsqnVsCentr);
        }
    }

    //EP / Qn resolutions
    TString detLabels[3][2] = {{"A","B"},{"A","C"},{"B","C"}};
    for(int iResoHisto=0; iResoHisto<3; iResoHisto++) {
        if(fFlowMethod==kEP || fFlowMethod==kEvShapeEP || fFlowMethod==kEPVsMass || fFlowMethod==kEvShapeEPVsMass) {
            fHistEPResolVsCentrVsqn[iResoHisto] = new TH3F(Form("fHistEvPlaneReso%d",iResoHisto+1),Form("Event plane angle Resolution;centrality (%%);%s;cos2(#psi_{%s}-#psi_{%s})",qnaxisnamefill.Data(),detLabels[iResoHisto][0].Data(),detLabels[iResoHisto][1].Data()),ncentbins,fMinCentr,fMaxCentr,nqnbins,qnmin,qnmax,220,-1.1,1.1);
        }
        else if(fFlowMethod==kSP || fFlowMethod==kEvShapeSP) {
            fHistEPResolVsCentrVsqn[iResoHisto] = new TH3F(Form("hScalProdQnVectors%d",iResoHisto+1),Form("Scalar product between Q-vectors;centrality (%%);%s;(Q{%s}Q{%s})",qnaxisnamefill.Data(),detLabels[iResoHisto][0].Data(),detLabels[iResoHisto][1].Data()),ncentbins,fMinCentr,fMaxCentr,nqnbins,qnmin,qnmax,200,-fScalProdLimit*fScalProdLimit,fScalProdLimit*fScalProdLimit);
        }
        fOutput->Add(fHistEPResolVsCentrVsqn[iResoHisto]);
    }

    //Sparse with D-meson candidates and useful azimuthal info
    int nptbins          = 200;
    double ptmin         = 0.;
    double ptmax         = 50.;

    int ndeltaphibins    = 1;
    double mindeltaphi   = -1.;
    double maxdeltaphi   = -1.;

    int nfphibins        = 100;
    double fphimin       = -1.;
    double fphimax       = 1.;

    TString deltaphiname = "";
    TString nfphiname1 = Form("Cos(%d#varphi_{D})",fHarmonic);
    TString nfphiname2 = Form("Sin(%d#varphi_{D})",fHarmonic);
    if(fFlowMethod==kEP || fFlowMethod==kEvShapeEP) {
        ndeltaphibins    = 96;
        mindeltaphi      = 0.;
        maxdeltaphi      = TMath::Pi();
        deltaphiname     = "#Delta#varphi";
    }
    else if(fFlowMethod==kSP || fFlowMethod==kEvShapeSP) {
        ndeltaphibins    = 100;
        mindeltaphi      = -fScalProdLimit*fScalProdLimit;
        maxdeltaphi      = fScalProdLimit*fScalProdLimit;
        deltaphiname     = Form("u_{%d,D}Q_{%d,A}",fHarmonic,fHarmonic);
        fphimin          = -fScalProdLimit*fScalProdLimit;
        fphimax          = fScalProdLimit*fScalProdLimit;
        nfphiname1       = Form("Cos(%d#varphi_{D})Q_{%d,y,A}",fHarmonic,fHarmonic);
        nfphiname2       = Form("Sin(%d#varphi_{D})Q_{%d,y,A}",fHarmonic,fHarmonic);
    }
    else if(fFlowMethod==kEPVsMass || fFlowMethod==kEvShapeEPVsMass) {
        ndeltaphibins    = 100;
        mindeltaphi      = -1.;
        maxdeltaphi      = 1.;
        deltaphiname     = Form("cos(%d#Delta#varphi)",fHarmonic);
    }

    int nphibins         = 18;
    double phimin        = 0.;
    double phimax        = 2*TMath::Pi();

    int nNtrkBins        = 100;
    double Ntrkmin       = 0.;
    double Ntrkmax       = 5000.;

    TString massaxisname = "";
    if(fDecChannel==0)      massaxisname = "#it{M}(K#pi#pi) (GeV/#it{c}^{2})";
    else if(fDecChannel==1) massaxisname = "#it{M}(K#pi) (GeV/#it{c}^{2})";
    else if(fDecChannel==2) massaxisname = "#it{M}(K#pi#pi)-#it{M}(K#pi) (GeV/#it{c}^{2})";
    else if(fDecChannel==3) massaxisname = "#it{M}(KK#pi) (GeV/#it{c}^{2})";

    int naxes=kVarForSparse;

    int nbins[kVarForSparse]     = {fNMassBins, nptbins, ndeltaphibins, nfphibins, nfphibins, nphibins, ncentbins, nNtrkBins, nqnbins};
    double xmin[kVarForSparse]   = {fLowmasslimit, ptmin, mindeltaphi, fphimin, fphimin, phimin, fMinCentr, Ntrkmin, qnmin};
    double xmax[kVarForSparse]   = {fUpmasslimit, ptmax, maxdeltaphi, fphimax, fphimax, phimax, fMaxCentr, Ntrkmax, qnmax};
    TString axTit[kVarForSparse] = {massaxisname, "#it{p}_{T} (GeV/#it{c})", deltaphiname, nfphiname1, nfphiname2, "#varphi_{D}", "Centrality (%)", "#it{N}_{tracklets}", qnaxisnamefill};

    fHistMassPtPhiqnCentr = new THnSparseF("fHistMassPtPhiqnCentr",Form("InvMass vs. #it{p}_{T} vs. %s vs. centr vs. #it{q}_{%d} ",deltaphiname.Data(),fHarmonic),naxes,nbins,xmin,xmax);
    fHistMassPtPhiqnCentr_1 = new THnSparseF("fHistMassPtPhiqnCentr_1",Form("InvMass vs. #it{p}_{T} vs. %s vs. centr vs. #it{q}_{%d} ",deltaphiname.Data(),fHarmonic),naxes,nbins,xmin,xmax);
    fHistMassPtPhiqnCentr_2 = new THnSparseF("fHistMassPtPhiqnCentr_2",Form("InvMass vs. #it{p}_{T} vs. %s vs. centr vs. #it{q}_{%d} ",deltaphiname.Data(),fHarmonic),naxes,nbins,xmin,xmax);
    fHistMassPtPhiqnCentr_3 = new THnSparseF("fHistMassPtPhiqnCentr_3",Form("InvMass vs. #it{p}_{T} vs. %s vs. centr vs. #it{q}_{%d} ",deltaphiname.Data(),fHarmonic),naxes,nbins,xmin,xmax);
    fHistMassPtPhiqnCentr_4 = new THnSparseF("fHistMassPtPhiqnCentr_4",Form("InvMass vs. #it{p}_{T} vs. %s vs. centr vs. #it{q}_{%d} ",deltaphiname.Data(),fHarmonic),naxes,nbins,xmin,xmax);
    fHistMassPtPhiqnCentr_5 = new THnSparseF("fHistMassPtPhiqnCentr_5",Form("InvMass vs. #it{p}_{T} vs. %s vs. centr vs. #it{q}_{%d} ",deltaphiname.Data(),fHarmonic),naxes,nbins,xmin,xmax);
    
    for(int iAx=0; iAx<naxes; iAx++){
        fHistMassPtPhiqnCentr->GetAxis(iAx)->SetTitle(axTit[iAx].Data());
        fHistMassPtPhiqnCentr_1->GetAxis(iAx)->SetTitle(axTit[iAx].Data());
        fHistMassPtPhiqnCentr_2->GetAxis(iAx)->SetTitle(axTit[iAx].Data());
        fHistMassPtPhiqnCentr_3->GetAxis(iAx)->SetTitle(axTit[iAx].Data());
        fHistMassPtPhiqnCentr_4->GetAxis(iAx)->SetTitle(axTit[iAx].Data());
        fHistMassPtPhiqnCentr_5->GetAxis(iAx)->SetTitle(axTit[iAx].Data());
    }
    
    fOutput->Add(fHistMassPtPhiqnCentr);
    fOutput->Add(fHistMassPtPhiqnCentr_1);
    fOutput->Add(fHistMassPtPhiqnCentr_2);
    fOutput->Add(fHistMassPtPhiqnCentr_3);
    fOutput->Add(fHistMassPtPhiqnCentr_4);
    fOutput->Add(fHistMassPtPhiqnCentr_5);

 //   fListBDTNtuple = new TList();fListBDTNtuple->SetOwner(); fListBDTNtuple->SetName("NtupleList");
 //   fListRDHFBDT->SetOwner(); fListBDTNtuple->SetName("BDTList");
    //ML model

    PostData(1,fOutput);
  //  PostData(3,fListBDTNtuple);
    return;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronvnTMVA::UserExec(Option_t */*option*/)
{
    
    // Execute analysis for current event:
    // heavy flavor candidates association to MC truth
    fAOD = dynamic_cast<AliAODEvent*> (InputEvent());
    fHistNEvents->Fill(0);

    //   Protection against the mismatch of candidate TRefs:
    //   Check if AOD and corresponding deltaAOD files contain the same number of events.
    //   In case of discrepancy the event is rejected.
    if(fAODProtection>=0) {
        //   Protection against different number of events in the AOD and deltaAOD
        //   In case of discrepancy the event is rejected.
        int matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
        if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
            // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
            fHistNEvents->Fill(2);
            return;
        }
        fHistNEvents->Fill(1);
    }

    // Post the data already here
    PostData(1,fOutput);

    TClonesArray *arrayProng = nullptr, *arrayD0toKpi = nullptr;
    int absPdgMom=0;
    int nDau=0;
    if(!fAOD && AODEvent() && IsStandardAOD()) {
        // In case there is an AOD handler writing a standard AOD, use the AOD
        // event in memory rather than the input (ESD) event.
        fAOD = dynamic_cast<AliAODEvent*> (AODEvent());
        // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
        // have to taken from the AOD event hold by the AliAODExtension
        AliAODHandler* aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
        if(aodHandler->GetExtensions()) {

            AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
            AliAODEvent *aodFromExt = ext->GetAOD();

            switch(fDecChannel) {
                case kDplustoKpipi:
                    absPdgMom=411;
                    nDau=3;
                    arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
                break;
                case kD0toKpi:
                    absPdgMom=421;
                    nDau=2;
                    arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
                break;
                case kDstartoKpipi:
                    absPdgMom=413;
                    nDau=3;
                    arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
                    arrayD0toKpi=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
                break;
                case kDstoKKpi:
                    absPdgMom=431;
                    nDau=3;
                    arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
                break;
            }
        }
    }
    else if(fAOD){
        switch(fDecChannel) {
            case kDplustoKpipi:
                absPdgMom=411;
                nDau = 3;
                arrayProng=(TClonesArray*)fAOD->GetList()->FindObject("Charm3Prong");
            break;
            case kD0toKpi:
                absPdgMom=421;
                nDau = 2;
                arrayProng=(TClonesArray*)fAOD->GetList()->FindObject("D0toKpi");
            break;
            case kDstartoKpipi:
                absPdgMom=413;
                nDau = 3;
                arrayProng=(TClonesArray*)fAOD->GetList()->FindObject("Dstar");
                arrayD0toKpi=(TClonesArray*)fAOD->GetList()->FindObject("D0toKpi");
            break;
            case kDstoKKpi:
                absPdgMom=431;
                nDau = 3;
                arrayProng=(TClonesArray*)fAOD->GetList()->FindObject("Charm3Prong");
            break;
        }
    }

    if(!fAOD || !arrayProng || (!arrayD0toKpi && fDecChannel==kDstartoKpipi)) {
        AliError("AliAnalysisTaskSECharmHadronvnTMVA::UserExec:Branch not found!\n");
        return;
    }

    // fix for temporary bug in ESDfilter
    // the AODs with null vertex pointer didn't pass the PhysSel
    if(!fAOD->GetPrimaryVertex() || TMath::Abs(fAOD->GetMagneticField())<0.001) return;

    fHistNEvents->Fill(3);

    double evCentr = fRDCuts->GetCentrality(fAOD);
    bool isEvSel = fRDCuts->IsEventSelected(fAOD);
    fHistCentrality[0]->Fill(evCentr);

    if(!isEvSel) {
        if(fRDCuts->IsEventRejectedDueToTrigger()) fHistNEvents->Fill(5);
        if(fRDCuts->IsEventRejectedDueToNotRecoVertex()) fHistNEvents->Fill(6);
        if(fRDCuts->IsEventRejectedDueToVertexContributors()) fHistNEvents->Fill(7);
        if(fRDCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion()) fHistNEvents->Fill(8);
        if(fRDCuts->IsEventRejectedDueToPileup()) fHistNEvents->Fill(9);
        if(fRDCuts->IsEventRejectedDueToCentrality()){
            fHistNEvents->Fill(10);
            fHistCentrality[2]->Fill(evCentr);
        }
        return;
    }

    fHistNEvents->Fill(4);
    fHistCentrality[1]->Fill(evCentr);

    int tracklets = AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(fAOD,-1.,1.);

    //Get Qn-vectors from tender task
    AliHFQnVectorHandler *HFQnVectorHandler = nullptr;
    bool isHandlerFound = false;
    AliAnalysisTaskSEHFTenderQnVectors *HFQnVectorTask = dynamic_cast<AliAnalysisTaskSEHFTenderQnVectors*>(AliAnalysisManager::GetAnalysisManager()->GetTask(fTenderTaskName.Data()));
    if(HFQnVectorTask) {
        HFQnVectorHandler = HFQnVectorTask->GetQnVectorHandler();

        if(fPercentileqn && !fqnSplinesList[0]) {
            for(int iDet=0; iDet<6; iDet++) {
                fqnSplinesList[iDet] = dynamic_cast<TList*>(HFQnVectorTask->GetSplineForqnPercentileList(iDet));
            }
        }
    }

    if(HFQnVectorHandler) {
        isHandlerFound = true;
        if(HFQnVectorHandler->GetHarmonic() != fHarmonic) {
            AliWarning("Harmonic of task and Qn-vector handler not consistent!");
            return;
        }
        if(HFQnVectorHandler->GetCalibrationType() != fCalibType) {
            AliWarning("Calibration strategy of task and Qn-vector handler not consistent!");
            return;
        }
        if(HFQnVectorHandler->GetNormalisationMethod() != fNormMethod) {
            AliWarning("Normalisation method of task and Qn-vector handler not consistent!");
            return;
        }
        if(fCalibType==AliHFQnVectorHandler::kQnCalib && HFQnVectorHandler->GetCalibrationsOADBFileName() != fOADBFileName) {
            AliWarning("OADB file name for calibrations of task and Qn-vector handler not consistent!");
            return;
        }
    }
    else { //create a new handler if not found in tender task
        AliWarning("Qn-vector tender task not found! Create a new one");
        HFQnVectorHandler = new AliHFQnVectorHandler(fCalibType,fNormMethod,fHarmonic,fOADBFileName);
        HFQnVectorHandler->SetAODEvent(fAOD);
        HFQnVectorHandler->ComputeCalibratedQnVectorTPC();
        HFQnVectorHandler->ComputeCalibratedQnVectorV0();
    }

    if(fPercentileqn && !fqnSplinesList[0]) { //probably qn splines not found in tender task, let's load them here
        fLoadedSplines=LoadSplinesForqnPercentile();
    }

    double QnFullTPC[2], QnPosTPC[2], QnNegTPC[2];
    double QnFullV0[2], QnV0A[2], QnV0C[2];
    double MultQnFullTPC = -1., MultQnPosTPC = -1., MultQnNegTPC = -1.;
    double MultQnFullV0 = -1., MultQnV0A = -1., MultQnV0C = -1.;
    double PsinFullTPC = -1., PsinPosTPC = -1., PsinNegTPC = -1.;
    double PsinFullV0 = -1., PsinV0A = -1., PsinV0C = -1.;

    //get the unnormalised Qn-vectors --> normalisation can be done in the task
    HFQnVectorHandler->GetUnNormQnVecTPC(QnFullTPC,QnPosTPC,QnNegTPC);
    HFQnVectorHandler->GetUnNormQnVecV0(QnFullV0,QnV0A,QnV0C);

    HFQnVectorHandler->GetMultQnVecTPC(MultQnFullTPC,MultQnPosTPC,MultQnNegTPC);
    HFQnVectorHandler->GetMultQnVecV0(MultQnFullV0,MultQnV0A,MultQnV0C);

    HFQnVectorHandler->GetEventPlaneAngleTPC(PsinFullTPC,PsinPosTPC,PsinNegTPC);
    HFQnVectorHandler->GetEventPlaneAngleV0(PsinFullV0,PsinV0A,PsinV0C);

    double mainQn[2], SubAQn[2], SubBQn[2], SubCQn[2];
    double mainPsin = -1., SubAPsin = -1., SubBPsin = -1., SubCPsin = -1.;
    double mainMultQn = -1., SubAMultQn = -1., SubBMultQn = -1., SubCMultQn = -1.;
    GetMainQnVectorInfo(mainPsin, mainMultQn, mainQn, SubAPsin, SubAMultQn, SubAQn, SubBPsin, SubBMultQn, SubBQn, HFQnVectorHandler);
    SubCPsin   = mainPsin;
    SubCMultQn = mainMultQn;
    SubCQn[0]  = mainQn[0];
    SubCQn[1]  = mainQn[1];

    int nsubevents = 3;
    if(fEvPlaneDet==fSubEvDetA || fEvPlaneDet==fSubEvDetB)
        nsubevents = 2;

    double qnFullTPC = -1., qnPosTPC = -1., qnNegTPC = -1.;
    double qnFullV0 = -1., qnV0A = -1., qnV0C = -1.;
    double mainqn = -1., mainpercqn = -1.;
    TSpline3* qnspline = nullptr;

    //define vars for candidate loop already here, in case daughter tracks removed from qn
    AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
    AliAODRecoDecayHF *d = nullptr;
    AliAODRecoDecayHF2Prong *dD0 = nullptr;

    int nCand = arrayProng->GetEntriesFast();
    bool alreadyLooped = false;
    vector<int> isSelByPrevLoop;
    vector<int> isBDTByPrevLoop;
    vector<AliAODTrack*> trackstoremove;

    if(fFlowMethod==kEvShapeEP || fFlowMethod==kEvShapeSP || fFlowMethod==kEvShapeEPVsMass) {
        HFQnVectorHandler->GetqnTPC(qnFullTPC,qnPosTPC,qnNegTPC);
        HFQnVectorHandler->GetqnV0(qnFullV0,qnV0A,qnV0C);
        if((fEnableDownsamplqn && fFracToKeepDownSamplqn<1.) || fRemoveDauFromqn==2) {

            //random downsampling of tracks for qn
            if(fEnableDownsamplqn && fFracToKeepDownSamplqn<1.) {
                for(int iTrack=0; iTrack<fAOD->GetNumberOfTracks(); iTrack++) {
                    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
                    double pseudoRand = track->Pt()*1000.-(long)(track->Pt()*1000);
                    if(pseudoRand>fFracToKeepDownSamplqn)
                        trackstoremove.push_back(track);
                }
            }

            //remove all daughter tracks of the analysed event, if enabled
            if(fRemoveDauFromqn==2) {
                alreadyLooped = true;
                for (int iCand = 0; iCand < nCand; iCand++) { //loop already here on candidates
                    d = dynamic_cast<AliAODRecoDecayHF*>(arrayProng->UncheckedAt(iCand));
                 AliAODRecoDecayHF2Prong *dd0 = (AliAODRecoDecayHF2Prong*)arrayProng->UncheckedAt(iCand);
                    if(fDecChannel==kDstartoKpipi) {
                        if(d->GetIsFilled()<1)
                            dD0 = dynamic_cast<AliAODRecoDecayHF2Prong*>(arrayD0toKpi->At(d->GetProngID(1)));
                        else
                            dD0 = (dynamic_cast<AliAODRecoCascadeHF*>(d))->Get2Prong();
                    }
                    double modelPred[2] = {-1., -1.};
                    int isSel = IsCandidateSelected(d, nDau, absPdgMom, vHF, dD0);
                   int isBDT = ProcessBDT(fAOD,dd0,isSel,vHF);
                    isSelByPrevLoop.push_back(isSel);
                    isBDTByPrevLoop.push_back(isBDT);
                    if(!isSel) continue;
                    if(!isBDT) continue;

                    GetDaughterTracksToRemove(d,nDau,trackstoremove);
                }
            }

            double QnFullTPCDownSampl[2], QnPosTPCDownSampl[2], QnNegTPCDownSampl[2];
            double MultQnFullTPCDownSampl = -1., MultQnPosTPCDownSampl = -1., MultQnNegTPCDownSampl = -1.;
            HFQnVectorHandler->RemoveTracksFromQnTPC(trackstoremove, QnFullTPCDownSampl, QnPosTPCDownSampl, QnNegTPCDownSampl, MultQnFullTPCDownSampl, MultQnPosTPCDownSampl, MultQnNegTPCDownSampl, true);
            qnFullTPC = TMath::Sqrt((QnFullTPCDownSampl[0]*QnFullTPCDownSampl[0]+QnFullTPCDownSampl[1]*QnFullTPCDownSampl[1])/MultQnFullTPCDownSampl);
            qnPosTPC = TMath::Sqrt((QnPosTPCDownSampl[0]*QnPosTPCDownSampl[0]+QnPosTPCDownSampl[1]*QnPosTPCDownSampl[1])/MultQnPosTPCDownSampl);
            qnNegTPC = TMath::Sqrt((QnNegTPCDownSampl[0]*QnNegTPCDownSampl[0]+QnNegTPCDownSampl[1]*QnNegTPCDownSampl[1])/MultQnNegTPCDownSampl);
        }

        TAxis* centraxis = fHistqnVsCentrPercCalib[0]->GetXaxis();
        int centrbin = centraxis->FindBin(evCentr);
        double centrbinmin = centraxis->GetBinLowEdge(centrbin);
        double centrbinmax = centrbinmin+centraxis->GetBinWidth(centrbin);

        switch(fqnMeth) {
            case kq2TPC:
                mainqn = qnFullTPC;
                if(fPercentileqn) qnspline = static_cast<TSpline3*>(fqnSplinesList[0]->FindObject(Form("sq2Int_centr_%0.f_%0.f",centrbinmin,centrbinmax)));
            break;
            case kq2PosTPC:
                mainqn = qnPosTPC;
                if(fPercentileqn) qnspline = static_cast<TSpline3*>(fqnSplinesList[1]->FindObject(Form("sq2Int_centr_%0.f_%0.f",centrbinmin,centrbinmax)));
            break;
            case kq2NegTPC:
                mainqn = qnNegTPC;
                if(fPercentileqn) qnspline = static_cast<TSpline3*>(fqnSplinesList[2]->FindObject(Form("sq2Int_centr_%0.f_%0.f",centrbinmin,centrbinmax)));
            break;
            case kq2VZERO:
                mainqn = qnFullV0;
                if(fPercentileqn) qnspline = static_cast<TSpline3*>(fqnSplinesList[3]->FindObject(Form("sq2Int_centr_%0.f_%0.f",centrbinmin,centrbinmax)));
            break;
            case kq2VZEROA:
                mainqn = qnV0A;
                if(fPercentileqn) qnspline = static_cast<TSpline3*>(fqnSplinesList[4]->FindObject(Form("sq2Int_centr_%0.f_%0.f",centrbinmin,centrbinmax)));
            break;
            case kq2VZEROC:
                mainqn = qnV0C;
                if(fPercentileqn) qnspline = static_cast<TSpline3*>(fqnSplinesList[5]->FindObject(Form("sq2Int_centr_%0.f_%0.f",centrbinmin,centrbinmax)));
            break;
        }
        if(fPercentileqn) {
            if(qnspline)
                mainpercqn = qnspline->Eval(mainqn);
            else
                AliWarning("Centrality binning and centrality intervals of qn splines do not match!");
        }
        else {
            mainpercqn = mainqn;
        }
    }
    //Fill event-based histos
    //EP and q2 histos
    fHistEvPlaneQncorr[0]->Fill(evCentr,mainpercqn,PsinFullTPC);
    fHistEvPlaneQncorr[1]->Fill(evCentr,mainpercqn,PsinPosTPC);
    fHistEvPlaneQncorr[2]->Fill(evCentr,mainpercqn,PsinNegTPC);
    fHistEvPlaneQncorr[3]->Fill(evCentr,mainpercqn,PsinFullV0);
    fHistEvPlaneQncorr[4]->Fill(evCentr,mainpercqn,PsinV0A);
    fHistEvPlaneQncorr[5]->Fill(evCentr,mainpercqn,PsinV0C);

    if(fFlowMethod==kEvShapeEP || fFlowMethod==kEvShapeSP || fFlowMethod==kEvShapeEPVsMass) {
        fHistqnVsCentrPercCalib[0]->Fill(evCentr,qnFullTPC);
        fHistqnVsCentrPercCalib[1]->Fill(evCentr,qnPosTPC);
        fHistqnVsCentrPercCalib[2]->Fill(evCentr,qnNegTPC);
        fHistqnVsCentrPercCalib[3]->Fill(evCentr,qnFullV0);
        fHistqnVsCentrPercCalib[4]->Fill(evCentr,qnV0A);
        fHistqnVsCentrPercCalib[5]->Fill(evCentr,qnV0C);

        fHistNtrklVsqnVsCentr->Fill(evCentr,mainpercqn,tracklets);

        if(fPercentileqn) {
            fHistPercqnVsqnVsCentr->Fill(evCentr,mainqn,mainpercqn);
        }
    }

    //EP / Qn resolution histos
    if(fFlowMethod==kEP || fFlowMethod==kEvShapeEP || fFlowMethod==kEPVsMass || fFlowMethod==kEvShapeEPVsMass) {
        fHistEPResolVsCentrVsqn[0]->Fill(evCentr,mainpercqn,TMath::Cos(fHarmonic*GetDeltaPsiSubInRange(SubAPsin,SubBPsin)));
        if(nsubevents==3) {
            fHistEPResolVsCentrVsqn[1]->Fill(evCentr,mainpercqn,TMath::Cos(fHarmonic*GetDeltaPsiSubInRange(SubAPsin,SubCPsin)));
            fHistEPResolVsCentrVsqn[2]->Fill(evCentr,mainpercqn,TMath::Cos(fHarmonic*GetDeltaPsiSubInRange(SubBPsin,SubCPsin)));
        }
    }
    else if(fFlowMethod==kSP || fFlowMethod==kEvShapeSP) {
        fHistEPResolVsCentrVsqn[0]->Fill(evCentr,mainpercqn,(SubAQn[0]*SubBQn[0]+SubAQn[1]*SubBQn[1])/(SubAMultQn * SubBMultQn));
        if(nsubevents==3) {
            fHistEPResolVsCentrVsqn[1]->Fill(evCentr,mainpercqn,(SubAQn[0]*SubCQn[0]+SubAQn[1]*SubCQn[1])/(SubAMultQn * SubCMultQn));
            fHistEPResolVsCentrVsqn[2]->Fill(evCentr,mainpercqn,(SubBQn[0]*SubCQn[0]+SubBQn[1]*SubCQn[1])/(SubCMultQn * SubBMultQn));
        }
    }
    //Loop on D candidates
    int nSelCand = 0;
    for (int iCand = 0; iCand < nCand; iCand++) {
        d = dynamic_cast<AliAODRecoDecayHF*>(arrayProng->UncheckedAt(iCand));
        AliAODRecoDecayHF2Prong *dd0 = (AliAODRecoDecayHF2Prong*)arrayProng->UncheckedAt(iCand);
        if(fDecChannel==kDstartoKpipi) {
            if(d->GetIsFilled()<1)
                dD0 = dynamic_cast<AliAODRecoDecayHF2Prong*>(arrayD0toKpi->At(d->GetProngID(1)));
            else
                dD0 = (dynamic_cast<AliAODRecoCascadeHF*>(d))->Get2Prong();
        }

        int isSelected = 0;
        int isBDT = 0;
        double modelPred[2] = {-1., -1.};
        if(alreadyLooped) {//check already here to avoid application of cuts twice
            isSelected = isSelByPrevLoop[iCand];
            isBDT = isBDTByPrevLoop[iCand];
        }
        else {
            isSelected = IsCandidateSelected(d, nDau, absPdgMom, vHF, dD0);
            isBDT = ProcessBDT(fAOD,dd0,isSelected,vHF);
         //     printf("isBDT =  %d\n",isBDT);
         //   printf("isSelected =  %d\n",isSelected);

        }
        if(!isSelected) continue;
        if(!isBDT) continue;

        nSelCand++;
        fHistNEvents->Fill(14); // candidate selected (including ML selection if used)

        float* invMass = nullptr;
        int nmasses = 0;
        CalculateInvMasses(d,invMass,nmasses);

        double ptD = d->Pt();
        double phiD = d->Phi();

        double deltaphi = GetPhiInRange(phiD-mainPsin);
        double scalprod = (TMath::Cos(fHarmonic*phiD)*mainQn[0]+TMath::Sin(fHarmonic*phiD)*mainQn[1]) / mainMultQn;
        double cosndeltaphi = TMath::Cos(fHarmonic*deltaphi);
        double vnfunc = -999.;
        double phifunc1 = -999.;
        double phifunc2 = -999.;
        if(fFlowMethod==kEP || fFlowMethod==kEvShapeEP) {
            vnfunc   = deltaphi;
            phifunc1 = TMath::Cos(fHarmonic*phiD);
            phifunc2 = TMath::Sin(fHarmonic*phiD);
        }
        else if(fFlowMethod==kSP || fFlowMethod==kEvShapeSP) {
            vnfunc = scalprod;
            phifunc1 = TMath::Cos(fHarmonic*phiD)*mainQn[1] / mainMultQn;
            phifunc2 = TMath::Sin(fHarmonic*phiD)*mainQn[0] / mainMultQn;
        }
        else if(fFlowMethod==kEPVsMass || fFlowMethod==kEvShapeEPVsMass) {
            vnfunc   = cosndeltaphi;
            phifunc1 = TMath::Cos(fHarmonic*phiD);
            phifunc2 = TMath::Sin(fHarmonic*phiD);
        }

        //remove daughter tracks from qn (on top of random downsampling, if enabled)
        double candQnFullTPC[2], candQnPosTPC[2], candQnNegTPC[2];
        double candMultQnFullTPC = -1., candMultQnPosTPC = -1., candMultQnNegTPC = -1.;
        double candqn = -1., candpercqn = -1.;
        if(fRemoveDauFromqn==1 && (fqnMeth==kq2TPC || fqnMeth==kq2PosTPC || fqnMeth==kq2NegTPC) && (fFlowMethod==kEvShapeEP || fFlowMethod==kEvShapeSP || fFlowMethod==kEvShapeEPVsMass)) {

            vector<AliAODTrack*> trackstoremoveDau(trackstoremove);
            GetDaughterTracksToRemove(d,nDau,trackstoremoveDau);
            HFQnVectorHandler->RemoveTracksFromQnTPC(trackstoremoveDau, candQnFullTPC, candQnPosTPC, candQnNegTPC, candMultQnFullTPC, candMultQnPosTPC, candMultQnNegTPC, true);
            switch(fqnMeth) {
                case kq2TPC:
                    candqn = TMath::Sqrt((candQnFullTPC[0]*candQnFullTPC[0]+candQnFullTPC[1]*candQnFullTPC[1])/candMultQnFullTPC);
                break;
                case kq2PosTPC:
                    candqn = TMath::Sqrt((candQnPosTPC[0]*candQnPosTPC[0]+candQnPosTPC[1]*candQnPosTPC[1])/candMultQnPosTPC);
                break;
                case kq2NegTPC:
                    candqn = TMath::Sqrt((candQnNegTPC[0]*candQnNegTPC[0]+candQnNegTPC[1]*candQnNegTPC[1])/candMultQnNegTPC);
                break;
            }
            trackstoremoveDau.clear();
            if(fPercentileqn) {
                if(qnspline)
                    candpercqn = qnspline->Eval(candqn);
                else
                    AliWarning("Centrality binning and centrality intervals of qn splines do not match!");
            }
            else {
                candpercqn = candqn;
            }
        }
        else {
            candpercqn = mainpercqn;
        }

        switch(fDecChannel) {
            case kDplustoKpipi:
            {
                double sparsearray[9] = {invMass[0],ptD,vnfunc,phifunc1,phifunc2,phiD,evCentr,static_cast<double>(tracklets),candpercqn};
                fHistMassPtPhiqnCentr->Fill(sparsearray);
                break;
            }
            case kD0toKpi:
            {
                if(isSelected==1 || isSelected==3) {
                    double sparsearray[9] = {invMass[0],ptD,vnfunc,phifunc1,phifunc2,phiD,evCentr,static_cast<double>(tracklets),candpercqn};
                    // cout<< isBDT;
                    if (((isBDT/1000000)%10)==1)fHistMassPtPhiqnCentr->Fill(sparsearray);
                    if (((isBDT/100000)%10)==1)fHistMassPtPhiqnCentr_1->Fill(sparsearray);
                    if (((isBDT/10000)%10)==1)fHistMassPtPhiqnCentr_2->Fill(sparsearray);
                    if (((isBDT/1000)%10)==1)fHistMassPtPhiqnCentr_3->Fill(sparsearray);
                    if (((isBDT/100)%10)==1)fHistMassPtPhiqnCentr_4->Fill(sparsearray);
                    if (((isBDT/10)%10)==1)fHistMassPtPhiqnCentr_5->Fill(sparsearray);
                    
                }
                if(isSelected==2 || isSelected==3) {
                    // cout<< isBDT;
                    double sparsearray[9] = {invMass[1],ptD,vnfunc,phifunc1,phifunc2,phiD,evCentr,static_cast<double>(tracklets),candpercqn};
                    if (((isBDT/1000000)%10)==1)fHistMassPtPhiqnCentr->Fill(sparsearray);
                    if (((isBDT/100000)%10)==1)fHistMassPtPhiqnCentr_1->Fill(sparsearray);
                    if (((isBDT/10000)%10)==1)fHistMassPtPhiqnCentr_2->Fill(sparsearray);
                    if (((isBDT/1000)%10)==1)fHistMassPtPhiqnCentr_3->Fill(sparsearray);
                    if (((isBDT/100)%10)==1)fHistMassPtPhiqnCentr_4->Fill(sparsearray);
                    if (((isBDT/10)%10)==1)fHistMassPtPhiqnCentr_5->Fill(sparsearray);
                }
                break;
            }
            case kDstartoKpipi:
            {
                double sparsearray[9] = {invMass[0],ptD,vnfunc,phifunc1,phifunc2,phiD,evCentr,static_cast<double>(tracklets),candpercqn};
                fHistMassPtPhiqnCentr->Fill(sparsearray);
                break;
            }
            case kDstoKKpi:
            {
                if(isSelected&4) {
                    double sparsearray[9] = {invMass[0],ptD,vnfunc,phifunc1,phifunc2,phiD,evCentr,static_cast<double>(tracklets),candpercqn};
                    fHistMassPtPhiqnCentr->Fill(sparsearray);
                }
                if(isSelected&8) {
                    double sparsearray[9] = {invMass[1],ptD,vnfunc,phifunc1,phifunc2,phiD,evCentr,static_cast<double>(tracklets),candpercqn};
                    fHistMassPtPhiqnCentr->Fill(sparsearray);
                }
                break;
            }
        }
    }

    fHistCandVsCent->Fill(evCentr,nSelCand);

    trackstoremove.clear();
    isSelByPrevLoop.clear();

    fDaughterTracks.Clear();

    delete vHF;
    vHF = nullptr;
    if(!isHandlerFound) { // if not found in the tender task, allocated memory with new --> to be deleted
        delete HFQnVectorHandler;
        HFQnVectorHandler = nullptr;
    }

    PostData(1,fOutput);
   // PostData(3,fListBDTNtuple);
    return;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronvnTMVA::SetQnVectorDetConf(int detconf)
{
  // Interface method to set some specific cases

    switch (detconf) {
        case kTPC:
            fEvPlaneDet=kPosTPC;
            fSubEvDetA=kPosTPC;
            fSubEvDetB=kNegTPC;
        break;
        case kTPCVZERO:
            fEvPlaneDet=kFullTPC;
            fSubEvDetA=kV0A;
            fSubEvDetB=kV0C;
        break;
        case kVZERO:
            fEvPlaneDet=kFullV0;
            fSubEvDetA=kPosTPC;
            fSubEvDetB=kNegTPC;
        break;
        case kVZEROA:
            fEvPlaneDet=kV0A;
            fSubEvDetA=kV0C;
            fSubEvDetB=kFullTPC;
        break;
        case kVZEROC:
            fEvPlaneDet=kV0C;
            fSubEvDetA=kV0A;
            fSubEvDetB=kFullTPC;
        break;
        case kPosTPCVZERO:
            fEvPlaneDet=kPosTPC;
            fSubEvDetA=kV0A;
            fSubEvDetB=kV0C;
        break;
        case kNegTPCVZERO:
            fEvPlaneDet=kNegTPC;
            fSubEvDetA=kV0A;
            fSubEvDetB=kV0C;
        break;
        default:
            AliWarning("Case not defined -> use full V0");
            fEvPlaneDet=kFullV0;
            fSubEvDetA=kPosTPC;
            fSubEvDetB=kNegTPC;
        break;
    }
}

//________________________________________________________________________
double AliAnalysisTaskSECharmHadronvnTMVA::GetPhiInRange(double phi)
{
    // Sets the phi angle in the range [0,2*pi/harmonic]

    double result = phi;
    while(result < 0) {
        result = result + 2. * TMath::Pi() / fHarmonic;
    }
    while(result > 2.*TMath::Pi() / fHarmonic){
        result = result - 2. * TMath::Pi() / fHarmonic;
    }
    return result;
}

//________________________________________________________________________
double AliAnalysisTaskSECharmHadronvnTMVA::GetDeltaPsiSubInRange(double psi1, double psi2)
{
    // difference of subevents reaction plane angle cannot be bigger than pi / n

    double delta = psi1 - psi2;
    if(TMath::Abs(delta) > TMath::Pi() / fHarmonic) {
        if(delta>0.) delta -= 2.*TMath::Pi() / fHarmonic;
        else delta += 2.*TMath::Pi() / fHarmonic;
    }

    return delta;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronvnTMVA::CalculateInvMasses(AliAODRecoDecayHF* d,float*& masses,int& nmasses)
{
  //Calculates all the possible invariant masses for each candidate
    switch(fDecChannel) {
        case kDplustoKpipi:
            {
                nmasses=1;
                masses=new float[nmasses];
                unsigned int pdgdaughters[3] = {211,321,211};
                masses[0]=d->InvMass(3,pdgdaughters);
            }
        break;
        case kD0toKpi:
            {
                nmasses=2;
                masses=new float[nmasses];
                unsigned int pdgdaughtersD0[2]={211,321};//pi,K
                masses[0]=d->InvMass(2,pdgdaughtersD0); //D0
                unsigned int pdgdaughtersD0bar[2]={321,211};//K,pi
                masses[1]=d->InvMass(2,pdgdaughtersD0bar); //D0bar
            }
        break;
        case kDstartoKpipi:
            {
                nmasses=1;
                masses=new float[nmasses];
                masses[0]=((AliAODRecoCascadeHF*)d)->DeltaInvMass();
            }
        break;
        case kDstoKKpi:
            {
                nmasses=2;
                masses=new float[nmasses];
                unsigned int pdgdaughtersKKpi[3] = {321,321,211};
                unsigned int pdgdaughterspiKK[3] = {211,321,321};
                masses[0]=d->InvMass(3,pdgdaughtersKKpi);
                masses[1]=d->InvMass(3,pdgdaughterspiKK);
            }
        break;
    }
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronvnTMVA::GetMainQnVectorInfo(double &mainPsin, double &mainMultQn, double mainQn[2], double &SubAPsin, double &SubAMultQn, double SubAQn[2], double &SubBPsin, double &SubBMultQn, double SubBQn[2], AliHFQnVectorHandler* HFQnVectorHandler)
{
    double QnFullTPC[2], QnPosTPC[2], QnNegTPC[2];
    double QnFullV0[2], QnV0A[2], QnV0C[2];
    double MultQnFullTPC = -1., MultQnPosTPC = -1., MultQnNegTPC = -1.;
    double MultQnFullV0 = -1., MultQnV0A = -1., MultQnV0C = -1.;
    double PsinFullTPC = -1., PsinPosTPC = -1., PsinNegTPC = -1.;
    double PsinFullV0 = -1., PsinV0A = -1., PsinV0C = -1.;

    //get the unnormalised Qn-vectors --> normalisation can be done in the task
    HFQnVectorHandler->GetUnNormQnVecTPC(QnFullTPC,QnPosTPC,QnNegTPC);
    HFQnVectorHandler->GetUnNormQnVecV0(QnFullV0,QnV0A,QnV0C);

    HFQnVectorHandler->GetMultQnVecTPC(MultQnFullTPC,MultQnPosTPC,MultQnNegTPC);
    HFQnVectorHandler->GetMultQnVecV0(MultQnFullV0,MultQnV0A,MultQnV0C);

    HFQnVectorHandler->GetEventPlaneAngleTPC(PsinFullTPC,PsinPosTPC,PsinNegTPC);
    HFQnVectorHandler->GetEventPlaneAngleV0(PsinFullV0,PsinV0A,PsinV0C);

    if(fEtaGapInTPCHalves>0.) {
        vector<AliAODTrack*> trackstoremove;
        for(int iTrack=0; iTrack<fAOD->GetNumberOfTracks(); iTrack++) {
            AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
            if(TMath::Abs(track->Eta())<fEtaGapInTPCHalves/2)
                trackstoremove.push_back(track);
        }
        HFQnVectorHandler->RemoveTracksFromQnTPC(trackstoremove, QnFullTPC, QnPosTPC, QnNegTPC, MultQnFullTPC, MultQnPosTPC, MultQnNegTPC, true);
        PsinFullTPC = (TMath::Pi()+TMath::ATan2(-QnFullTPC[1],-QnFullTPC[0]))/2;
        PsinPosTPC = (TMath::Pi()+TMath::ATan2(-QnPosTPC[1],-QnPosTPC[0]))/2;
        PsinNegTPC = (TMath::Pi()+TMath::ATan2(-QnNegTPC[1],-QnNegTPC[0]))/2;
    }

    switch(fEvPlaneDet) {
        case kFullTPC:
            mainPsin   = PsinFullTPC;
            mainMultQn = MultQnFullTPC;
            mainQn[0]  = QnFullTPC[0];
            mainQn[1]  = QnFullTPC[1];
        break;
        case kPosTPC:
            mainPsin   = PsinPosTPC;
            mainMultQn = MultQnPosTPC;
            mainQn[0]  = QnPosTPC[0];
            mainQn[1]  = QnPosTPC[1];
        break;
        case kNegTPC:
            mainPsin   = PsinNegTPC;
            mainMultQn = MultQnNegTPC;
            mainQn[0]  = QnNegTPC[0];
            mainQn[1]  = QnNegTPC[1];
        break;
        case kFullV0:
            mainPsin   = PsinFullV0;
            mainMultQn = MultQnFullV0;
            mainQn[0]  = QnFullV0[0];
            mainQn[1]  = QnFullV0[1];
        break;
        case kV0A:
            mainPsin   = PsinV0A;
            mainMultQn = MultQnV0A;
            mainQn[0]  = QnV0A[0];
            mainQn[1]  = QnV0A[1];
        break;
        case kV0C:
            mainPsin   = PsinV0C;
            mainMultQn = MultQnV0C;
            mainQn[0]  = QnV0C[0];
            mainQn[1]  = QnV0C[1];
        break;
    }

    switch(fSubEvDetA) {
        case kFullTPC:
            SubAPsin   = PsinFullTPC;
            SubAMultQn = MultQnFullTPC;
            SubAQn[0]  = QnFullTPC[0];
            SubAQn[1]  = QnFullTPC[1];
        break;
        case kPosTPC:
            SubAPsin   = PsinPosTPC;
            SubAMultQn = MultQnPosTPC;
            SubAQn[0]  = QnPosTPC[0];
            SubAQn[1]  = QnPosTPC[1];
        break;
        case kNegTPC:
            SubAPsin   = PsinNegTPC;
            SubAMultQn = MultQnNegTPC;
            SubAQn[0]  = QnNegTPC[0];
            SubAQn[1]  = QnNegTPC[1];
        break;
        case kFullV0:
            SubAPsin   = PsinFullV0;
            SubAMultQn = MultQnFullV0;
            SubAQn[0]  = QnFullV0[0];
            SubAQn[1]  = QnFullV0[1];
        break;
        case kV0A:
            SubAPsin   = PsinV0A;
            SubAMultQn = MultQnV0A;
            SubAQn[0]  = QnV0A[0];
            SubAQn[1]  = QnV0A[1];
        break;
        case kV0C:
            SubAPsin   = PsinV0C;
            SubAMultQn = MultQnV0C;
            SubAQn[0]  = QnV0C[0];
            SubAQn[1]  = QnV0C[1];
        break;
    }

    switch(fSubEvDetB) {
        case kFullTPC:
            SubBPsin   = PsinFullTPC;
            SubBMultQn = MultQnFullTPC;
            SubBQn[0]  = QnFullTPC[0];
            SubBQn[1]  = QnFullTPC[1];
        break;
        case kPosTPC:
            SubBPsin   = PsinPosTPC;
            SubBMultQn = MultQnPosTPC;
            SubBQn[0]  = QnPosTPC[0];
            SubBQn[1]  = QnPosTPC[1];
        break;
        case kNegTPC:
            SubBPsin   = PsinNegTPC;
            SubBMultQn = MultQnNegTPC;
            SubBQn[0]  = QnNegTPC[0];
            SubBQn[1]  = QnNegTPC[1];
        break;
        case kFullV0:
            SubBPsin   = PsinFullV0;
            SubBMultQn = MultQnFullV0;
            SubBQn[0]  = QnFullV0[0];
            SubBQn[1]  = QnFullV0[1];
        break;
        case kV0A:
            SubBPsin   = PsinV0A;
            SubBMultQn = MultQnV0A;
            SubBQn[0]  = QnV0A[0];
            SubBQn[1]  = QnV0A[1];
        break;
        case kV0C:
            SubBPsin   = PsinV0C;
            SubBMultQn = MultQnV0C;
            SubBQn[0]  = QnV0C[0];
            SubBQn[1]  = QnV0C[1];
        break;
    }
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronvnTMVA::GetDaughterTracksToRemove(AliAODRecoDecayHF* d, int nDau, vector<AliAODTrack*> &trackstoremove)
{
    if(fDecChannel!=kDstartoKpipi) {
        for(int iDau=0; iDau<nDau; iDau++) {
            AliAODTrack* dau = dynamic_cast<AliAODTrack*>(d->GetDaughter(iDau));
            trackstoremove.push_back(dau);
        }
    }
    else {
        AliAODRecoDecayHF2Prong *dauD0 = dynamic_cast<AliAODRecoCascadeHF*>(d)->Get2Prong();
        trackstoremove.push_back(dynamic_cast<AliAODTrack*>(dauD0->GetDaughter(0)));
        trackstoremove.push_back(dynamic_cast<AliAODTrack*>(dauD0->GetDaughter(1)));
        if(fRemoveSoftPion) trackstoremove.push_back(dynamic_cast<AliAODRecoCascadeHF*>(d)->GetBachelor());
    }
}

//________________________________________________________________________
int AliAnalysisTaskSECharmHadronvnTMVA::IsCandidateSelected(AliAODRecoDecayHF *&d, int nDau, int absPdgMom, AliAnalysisVertexingHF *vHF, AliAODRecoDecayHF2Prong *dD0) {
    if(!d || !vHF || (fDecChannel==kDstartoKpipi && !dD0)) return false;
    //Preselection to speed up task
    TObjArray arrDauTracks(nDau);
    AliAODTrack *track = nullptr;
    if(fDecChannel!=kDstartoKpipi) {
        for(int iDau=0; iDau<nDau; iDau++){
            AliAODTrack *track = vHF->GetProng(fAOD,d,iDau);
            arrDauTracks.AddAt(track,iDau);
        }
    }
    else {
        for(int iDau=0; iDau<nDau; iDau++){
            if(iDau == 0)
                track=vHF->GetProng(fAOD,d,iDau); //soft pion
            else
                track=vHF->GetProng(fAOD,dD0,iDau-1); //D0 daughters
            arrDauTracks.AddAt(track,iDau);
        }
    }
    if(!fRDCuts->PreSelect(arrDauTracks)){
        return 0;
    }

    bool isSelBit=true;
    switch(fDecChannel) {
        case kDplustoKpipi:
            isSelBit = d->HasSelectionBit(AliRDHFCuts::kDplusCuts);
            if(!isSelBit || !vHF->FillRecoCand(fAOD,(AliAODRecoDecayHF3Prong*)d)) return 0;
        break;
        case kD0toKpi:
            isSelBit = d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts);
            if(!isSelBit || !vHF->FillRecoCand(fAOD,(AliAODRecoDecayHF2Prong*)d)) return 0;
        break;
        case kDstartoKpipi:
            if(!vHF->FillRecoCasc(fAOD,((AliAODRecoCascadeHF*)d),true)) return 0;
        break;
        case kDstoKKpi:
            isSelBit = d->HasSelectionBit(AliRDHFCuts::kDsCuts);
            if(!isSelBit || !vHF->FillRecoCand(fAOD,(AliAODRecoDecayHF3Prong*)d)) return 0;
        break;
    }

    double ptD = d->Pt();
    double yD = d->Y(absPdgMom);
    int ptbin=fRDCuts->PtBin(ptD);
    if(ptbin<0) {
        fHistNEvents->Fill(15);
        return 0;
    }
    bool isFidAcc = fRDCuts->IsInFiducialAcceptance(ptD,yD);
    if(!isFidAcc) return 0;

    int isSelected = fRDCuts->IsSelected(d,AliRDHFCuts::kAll,fAOD);

    //ML application

    if(!isSelected) return 0;
    if(fDecChannel==kDstoKKpi)
        if(!(isSelected&4) && !(isSelected&8)) return 0;

    return isSelected;
}

//________________________________________________________________________
bool AliAnalysisTaskSECharmHadronvnTMVA::LoadSplinesForqnPercentile()
{
    // load splines from file

    TString listname[6] = {"SplineListq2TPC", "SplineListq2TPCPosEta", "SplineListq2TPCNegEta", "SplineListq2V0", "SplineListq2V0A", "SplineListq2V0C"};

    if (!gGrid) {
        TGrid::Connect("alien://");
    }
    TFile* splinesfile = TFile::Open(fqnSplineFileName.Data());
    if(!splinesfile) {
        AliFatal("File with splines for qn percentiles not found!");
        return false;
    }

    for(int iDet=0; iDet<6; iDet++) {
        fqnSplinesList[iDet] = (TList*)splinesfile->Get(listname[iDet].Data());
        if(!fqnSplinesList[iDet]) {
            AliFatal("TList with splines for qn percentiles not found in the spline file!");
            return false;
        }
        fqnSplinesList[iDet]->SetOwner();
    }
    splinesfile->Close();

    return true;
}

//________________________________________________________________________
int AliAnalysisTaskSECharmHadronvnTMVA::ProcessBDT(AliAODEvent *fAOD, AliAODRecoDecayHF2Prong *dD0, int isSelected, AliAnalysisVertexingHF *vHF){
    if(!vHF||!dD0||!fAOD) return 0;
    if(!isSelected) return 0;
     AliAODTrack *prong2 = vHF->GetProng(fAOD,dD0,0);
     AliAODTrack *prong3 = vHF->GetProng(fAOD,dD0,1);
    Double_t normIP[2];
    Double_t d0Prong[2];
    Double_t ptProng[2];
    
    Double_t invmassD0 = dD0->InvMassD0(); Double_t invmassD0bar = dD0->InvMassD0bar();
    Double_t ptB;
    
    Double_t cosPointingAngle;
    Double_t cosThetaStarD0 = 99;
    Double_t cosThetaStarD0bar = 99;
     Double_t diffIP[2], errdiffIP[2];
    dD0->Getd0MeasMinusExpProng(0,fAOD->GetMagneticField(),diffIP[0],errdiffIP[0]);
    dD0->Getd0MeasMinusExpProng(1,fAOD->GetMagneticField(),diffIP[1],errdiffIP[1]);
   if(!diffIP[0]||!errdiffIP[0])return 0;
    normIP[0]=diffIP[0]/errdiffIP[0];
    normIP[1]=diffIP[1]/errdiffIP[1];
    
    AliAODVertex *secvtx  = dD0->GetSecondaryVtx();
    AliAODVertex *primvtx = dD0->GetPrimaryVtx();
    Double_t err2decaylength = secvtx->Error2DistanceXYToVertex(primvtx);
    Double_t lxy = dD0->AliAODRecoDecay::DecayLengthXY(primvtx);
    Bool_t isusepid = fRDCuts->GetIsUsePID();
    //peng
    d0Prong[0]=dD0->Getd0Prong(0); d0Prong[1]=dD0->Getd0Prong(1); //d0*d0 d0Prong[1]*d0Prong[0]
    cosPointingAngle = dD0->CosPointingAngle();
    ptProng[0]=dD0->Pt2Prong(0); ptProng[1]=dD0->Pt2Prong(1);
    cosThetaStarD0 = dD0->CosThetaStarD0();
    cosThetaStarD0bar = dD0->CosThetaStarD0bar();
    //  if (part->Pt() > 10) return;
    // DCA part->GetDCA();
    Float_t tmp[22];
    tmp[8]= -99;                                // Invariant Mass
    tmp[18] = -99;                                 // ptB (if accessible)
    tmp[19] = 0;                                 // PDGCode
    tmp[20] = -99;                                // Rapidity YD0
    tmp[21] = -99;                                // Azimuthal Phi
    
    tmp[0] = dD0->Pt();                         // ptD
    tmp[1] = normIP[0];                         // Normalized d0N-<d0N> (topo1)
    tmp[2] = normIP[1];                         // Normalized d0P-<d0P> (topo2)
    tmp[3] = lxy;                                 // Decay length
    tmp[4] = lxy/TMath::Sqrt(err2decaylength);     // Normalized decay length
    tmp[9] = d0Prong[0]*d0Prong[1];             // d0N*d0P
    tmp[10] = cosPointingAngle;                 // CosThetaPointing
    tmp[11] = dD0->GetDCA();                     // DCAtracks
    tmp[12] = prong2->Pt();                     // pt1
    tmp[13] = prong3->Pt();                     // pt2
    tmp[14] = dD0->CosPointingAngleXY();         // CosThetaPointingXY
    tmp[15] = d0Prong[0];                         // d01
    tmp[16] = d0Prong[1];                         // d02
    tmp[20] = dD0->YD0();
    tmp[21] = dD0->Phi();

    // PID and Cuts
    if(isusepid)fRDCuts->SetUsePID(kFALSE);// if PID on, switch it off
    Int_t isCuts=fRDCuts->IsSelected(dD0,AliRDHFCuts::kAll,fAOD);
    if(isusepid)fRDCuts->SetUsePID(kTRUE);//if PID was on, switch it on
    Int_t isPid= fRDCuts->IsSelectedPID(dD0);
    tmp[5] = (Double_t)isCuts;
    tmp[6] = (Double_t)isPid;

    //~ std::vector<TString> vnameString; vnameString.resize(23);
    /// "ptD:topo1:topo2:lxy:nlxy:iscut:ispid:type:mass:d0d0:cosp:dca:ptk:ptpi:cospxy:d0k:d0pi:cosstar:ptB:pdgcode:YD0:phi
    //~ vnameString.push_back("ptD");     vnameString.push_back("topo1"); vnameString.push_back("topo2");     vnameString.push_back("lxy");     vnameString.push_back("nlxy");
    //~ vnameString.push_back("iscut"); vnameString.push_back("ispid"); vnameString.push_back("type");         vnameString.push_back("mass");     vnameString.push_back("d0d0");
    //~ vnameString.push_back("cosp");     vnameString.push_back("dca");     vnameString.push_back("ptk");         vnameString.push_back("ptpi");     vnameString.push_back("cospxy");
    //~ vnameString.push_back("d0k");     vnameString.push_back("d0pi");     vnameString.push_back("cosstar");     vnameString.push_back("ptB");     vnameString.push_back("pdgcode");
    //~ vnameString.push_back("YD0");     vnameString.push_back("phi");
  //  float BDTmatrix[13][2]={{0.11,0.13},/* 1<pt<2*/
  //      {0.08,0.13},/* 2<pt<3*/
  //      {0.11,0.15},/* 3<pt<4 */
  //      {0.10,0.16},/* 4<pt<5 */
  //      {0.07,0.16},/* 5<pt<6 */
  //      {0.08,0.12},/* 6<pt<7 */
  //      {0.08,0.08},/* 7<pt<8 */
  //      {0.05,0.08},/* 8<pt<10 */
  //      {0.05,0.06},/* 10<pt<12 */
  //      {0.04,0.09},/* 12<pt<16 */
   //     {0.03,0.05},/* 16<pt<24 */
   //     {0.01,0.04},/* 24<pt<36 */
   //     {0.01,0.04}};/* 36<pt<50 */
    Double_t BDTmatrix_1[13][2]={{-0.01,0.095},/* 1<pt<2*/
        {-0.01,0.16},/* 2<pt<3*/
        {-0.01,0.155},/* 3<pt<4 */
        {-0.01,0.165},/* 4<pt<5 */
        {-0.01,0.145},/* 5<pt<6 */
        {-0.01,0.14},/* 6<pt<7 */
        {-0.01,0.105},/* 7<pt<8 */
        {-0.01,0.11},/* 8<pt<10 */
        {0,0.065},/* 10<pt<12 */
        {-0.01,0.08},/* 12<pt<16 */
        {-0.01,0.095},/* 16<pt<24 */
        {-0.01,0.1},/* 24<pt<36 */
        {-0.01,0.1}};/* 36<pt<50 */
    
    Double_t BDTmatrix_2[13][2]={{0.01,0.095},/* 1<pt<2*/
        {0.01,0.16},/* 2<pt<3*/
        {0.01,0.16},/* 3<pt<4 */
        {0.01,0.16},/* 4<pt<5 */
        {0.01,0.15},/* 5<pt<6 */
        {0.01,0.14},/* 6<pt<7 */
        {0.01,0.11},/* 7<pt<8 */
        {0.01,0.115},/* 8<pt<10 */
        {0.01,0.065},/* 10<pt<12 */
        {0,0.08},/* 12<pt<16 */
        {0,0.095},/* 16<pt<24 */
        {0,0.1},/* 24<pt<36 */
        {0,0.1}};/* 36<pt<50 */
    Double_t BDTmatrix_3[13][2]={{0.03,0.135},/* 1<pt<2*/
        {0.03,0.16},/* 2<pt<3*/
        {0.02,0.16},/* 3<pt<4 */
        {0.02,0.145},/* 4<pt<5 */
        {0.02,0.145},/* 5<pt<6 */
        {0.03,0.12},/* 6<pt<7 */
        {0.02,0.095},/* 7<pt<8 */
        {0.02,0.09},/* 8<pt<10 */
        {0.02,0.05},/* 10<pt<12 */
        {0.015,0.095},/* 12<pt<16 */
        {0.015,0.09},/* 16<pt<24 */
        {0.005,0.08},/* 24<pt<36 */
        {0.005,0.08}};/* 36<pt<50 */
    Double_t BDTmatrix_4[13][2]={{0.05,0.1},/* 1<pt<2*/
        {0.05,0.16},/* 2<pt<3*/
        {0.04,0.16},/* 3<pt<4 */
        {0.04,0.145},/* 4<pt<5 */
        {0.04,0.145},/* 5<pt<6 */
        {0.05,0.115},/* 6<pt<7 */
        {0.04,0.1},/* 7<pt<8 */
        {0.04,0.09},/* 8<pt<10 */
        {0.04,0.065},/* 10<pt<12 */
        {0.025,0.095},/* 12<pt<16 */
        {0.025,0.09},/* 16<pt<24 */
        {0.015,0.1},/* 24<pt<36 */
        {0.015,0.1}};/* 36<pt<50 */
    Double_t BDTmatrix_5[13][2]={{0.07,0.14},/* 1<pt<2*/
        {0.07,0.135},/* 2<pt<3*/
        {0.06,0.16},/* 3<pt<4 */
        {0.06,0.155},/* 4<pt<5 */
        {0.05,0.16},/* 5<pt<6 */
        {0.06,0.12},/* 6<pt<7 */
        {0.05,0.11},/* 7<pt<8 */
        {0.05,0.08},/* 8<pt<10 */
        {0.055,0.055},/* 10<pt<12 */
        {0.035,0.095},/* 12<pt<16 */
        {0.03,0.095},/* 16<pt<24 */
        {0.02,0.1},/* 24<pt<36 */
        {0.02,0.1}};/* 36<pt<50 */
    Double_t BDTmatrix_6[13][2]={{0.085,0.135},/* 1<pt<2*/
        {0.08,0.135},/* 2<pt<3*/
        {0.08,0.16},/* 3<pt<4 */
        {0.08,0.155},/* 4<pt<5 */
        {0.07,0.15},/* 5<pt<6 */
        {0.07,0.12},/* 6<pt<7 */
        {0.07,0.13},/* 7<pt<8 */
        {0.07,0.075},/* 8<pt<10 */
        {0.065,0.055},/* 10<pt<12 */
        {0.045,0.095},/* 12<pt<16 */
        {0.035,0.1},/* 16<pt<24 */
        {0.03,0.1},/* 24<pt<36 */
        {0.03,0.1}};/* 36<pt<50 */



    
       // if(!isSelected) return 0;
  //  printf("isSelected =  %d\n",isSelected);
        int isBDTSelected = 0;
        std::vector<Double_t> BDTClsVar;// BDT cls input
        BDTClsVar.resize(10);
        // Data fill this
        Int_t thisptbin = fRDCuts->PtBin(tmp[0]);
        if(thisptbin<0) return 0;
        Float_t *ptbin = fRDCuts->GetPtBinLimits();
        TString ptstring = Form("_%.0f_%.0f",ptbin[thisptbin],ptbin[thisptbin+1]);
  //  printf("ptstring = _%.0f_%.0f\n",ptbin[thisptbin],ptbin[thisptbin+1]);
    Double_t bdt1resp1 = -3; Double_t bdt2resp1_0 = -3; Double_t bdt2resp1_1 = -3; Double_t bdt2resp1_2 = -3;
    Double_t bdt1resp = -3; Double_t bdt2resp_0 = -3; Double_t bdt2resp_1 = -3; Double_t bdt2resp_2 = -3;
        if(isSelected==1 || isSelected==3){
            tmp[7] = 1; tmp[8] = invmassD0; tmp[17] = cosThetaStarD0;
            if(tmp[8]>2.12||tmp[8]<1.65) return 0;
            
            // Link variables to be used as classifier
            BDTClsVar[0] = tmp[1];     BDTClsVar[1] = tmp[2];     BDTClsVar[2] = tmp[4];     BDTClsVar[3] = tmp[9];     BDTClsVar[4] = tmp[10];
            BDTClsVar[5] = tmp[11]; BDTClsVar[6] = tmp[14]; BDTClsVar[7] = tmp[15]; BDTClsVar[8] = tmp[16]; BDTClsVar[9] = tmp[17];
            
            
            AliRDHFBDT *thisbdt1   = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT1%s",ptstring.Data()));
            AliRDHFBDT *thisbdt2ll = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_0",ptstring.Data()));
            AliRDHFBDT *thisbdt2mm = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_1",ptstring.Data()));
            AliRDHFBDT *thisbdt2hh = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_2",ptstring.Data()));

              //  Float_t bdt1resp1; Float_t bdt2resp1;
            bdt1resp1 = thisbdt1->GetResponse(BDTClsVar);
            bdt2resp1_0 = thisbdt2ll->GetResponse(BDTClsVar);
            bdt2resp1_1 = thisbdt2mm->GetResponse(BDTClsVar);
            bdt2resp1_2 = thisbdt2hh->GetResponse(BDTClsVar);
     //       cout<< bdt1resp1;
      //      cout<< bdt2resp1;
          //     printf("bdt1resp1 =  %lf\n",thisbdt1);
          //   printf("bdt1resp2 =  %lf\n",thisbdt2ll);
                //~ bdt1resp1=0; bdt2resp1=0; bdt2resp[1]=0; bdt2resp[2]=0; bdt2resp[3]=0; bdt2resp[4]=0; bdt2resp[5]=0;
                // BDT Responses ready

            
        }
        if (isSelected>1){
            tmp[7] = 2; tmp[8] = invmassD0bar; tmp[17] = cosThetaStarD0bar;
            
            if(tmp[8]>2.12||tmp[8]<1.65) return 0;
            // Link variables to be used as classifier
            BDTClsVar[0] = tmp[1];     BDTClsVar[1] = tmp[2];     BDTClsVar[2] = tmp[4];     BDTClsVar[3] = tmp[9];     BDTClsVar[4] = tmp[10];
            BDTClsVar[5] = tmp[11]; BDTClsVar[6] = tmp[14]; BDTClsVar[7] = tmp[15]; BDTClsVar[8] = tmp[16]; BDTClsVar[9] = tmp[17];

        // Data application
            AliRDHFBDT *thisbdt1 = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT1%s",ptstring.Data()));
            AliRDHFBDT *thisbdt2ll = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_0",ptstring.Data()));
            AliRDHFBDT *thisbdt2mm = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_1",ptstring.Data()));
            AliRDHFBDT *thisbdt2hh = (AliRDHFBDT*)fListRDHFBDT->FindObject(Form("_BDT2%s_2",ptstring.Data()));
            
          //  thisbdt2ll->Print();
            //   Float_t bdt1resp; Float_t bdt2resp;
            bdt1resp = thisbdt1->GetResponse(BDTClsVar);
            bdt2resp_0 = thisbdt2ll->GetResponse(BDTClsVar);
            bdt2resp_1 = thisbdt2mm->GetResponse(BDTClsVar);
            bdt2resp_2 = thisbdt2hh->GetResponse(BDTClsVar);
         //  cout<< bdt1resp;
         //   cout<< bdt2resp;
          //  printf("bdt1resp1 =  %lf\n",thisbdt1);
           // printf("bdt1resp2 =  %lf\n",thisbdt2ll);
                //~ bdt1resp=0; bdt2resp[0]=0; bdt2resp[1]=0; bdt2resp[2]=0; bdt2resp[3]=0; bdt2resp[4]=0; bdt2resp[5]=0;
                // BDT Responses ready
 
            
        }
  //  printf("thisptbin =  %d\n",thisptbin);
  //  printf("BDTmatrix1 =  %f\n",fBDT1Cut[thisptbin]);
   // printf("BDTmatrix2 =  %f\n",fBDT2Cut[thisptbin]);
    int rep1 = 0; int rep2 = 0; int rep3 = 0; int rep4 = 0; int rep5 = 0; int rep6 = 0;
    if((bdt1resp1>BDTmatrix_1[thisptbin][0] && bdt2resp1_0>BDTmatrix_1[thisptbin][1])||(bdt1resp>BDTmatrix_1[thisptbin][0]&& bdt2resp_0>BDTmatrix_1[thisptbin][1])) rep1 = 1;
    if((bdt1resp1>BDTmatrix_2[thisptbin][0] && bdt2resp1_0>BDTmatrix_2[thisptbin][1])||(bdt1resp>BDTmatrix_2[thisptbin][0]&& bdt2resp_0>BDTmatrix_2[thisptbin][1])) rep2 = 1;
    if((bdt1resp1>BDTmatrix_3[thisptbin][0] && bdt2resp1_1>BDTmatrix_3[thisptbin][1])||(bdt1resp>BDTmatrix_3[thisptbin][0]&& bdt2resp_1>BDTmatrix_3[thisptbin][1])) rep3 = 1;
    if((bdt1resp1>BDTmatrix_4[thisptbin][0] && bdt2resp1_1>BDTmatrix_4[thisptbin][1])||(bdt1resp>BDTmatrix_4[thisptbin][0]&& bdt2resp_1>BDTmatrix_4[thisptbin][1])) rep4 = 1;
    if((bdt1resp1>BDTmatrix_5[thisptbin][0] && bdt2resp1_2>BDTmatrix_5[thisptbin][1])||(bdt1resp>BDTmatrix_5[thisptbin][0]&& bdt2resp_2>BDTmatrix_5[thisptbin][1])) rep5 = 1;
    if((bdt1resp1>BDTmatrix_6[thisptbin][0] && bdt2resp1_2>BDTmatrix_6[thisptbin][1])||(bdt1resp>BDTmatrix_6[thisptbin][0]&& bdt2resp_2>BDTmatrix_6[thisptbin][1])) rep6 = 1;
    isBDTSelected = 10*rep6+100*rep5+1000*rep4+10000*rep3+100000*rep2+1000000*rep1;
    return isBDTSelected;
    
    
}
