/* Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************************
// \class AliAnalysisTaskSECharmHadronvn
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
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsLc2V0bachelor.h"

#include "AliMultSelection.h"
#include "AliVertexingHFUtils.h"

#include "AliAnalysisTaskSEHFTenderQnVectors.h"
#include "AliAnalysisTaskSECharmHadronvn.h"

#include "AliAODPidHF.h"
#include "AliHFMLResponseDstoKKpi.h"
#include "AliHFMLResponseDplustoKpipi.h"

ClassImp(AliAnalysisTaskSECharmHadronvn)

//________________________________________________________________________
AliAnalysisTaskSECharmHadronvn::AliAnalysisTaskSECharmHadronvn() :
AliAnalysisTaskSE(),
    fAOD(nullptr),
    fOutput(nullptr),
    fHistNEvents(nullptr),
    fHistCandVsCent(nullptr),
    fHistPercqnVsqnVsCentr(nullptr),
    fHistNtrklVsqnVsCentr(nullptr),
    fHistMassPtPhiqnCentr(nullptr),
    fRDCuts(nullptr),
    fTenderTaskName("HFTenderQnVectors"),
    fMinCentr(0.),
    fMaxCentr(100.),
    fAODProtection(1),
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
    fRemoveSoftPion(false),
    fEnableDownsamplqn(false),
    fFracToKeepDownSamplqn(1.1),
    fApplyML(false),
    fConfigPath(""),
    fMLResponse(nullptr),
    fNMLBins(300),
    fMLOutputMin(0.85),
    fMLOutputMax(1.)
{
    // Default constructor
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
AliAnalysisTaskSECharmHadronvn::AliAnalysisTaskSECharmHadronvn(const char *name, AliRDHFCuts *rdCuts, int decaychannel) :
    AliAnalysisTaskSE(name),
    fAOD(nullptr),
    fOutput(nullptr),
    fHistNEvents(nullptr),
    fHistCandVsCent(nullptr),
    fHistPercqnVsqnVsCentr(nullptr),
    fHistNtrklVsqnVsCentr(nullptr),
    fHistMassPtPhiqnCentr(nullptr),
    fRDCuts(rdCuts),
    fTenderTaskName("HFTenderQnVectors"),
    fMinCentr(0.),
    fMaxCentr(100.),
    fAODProtection(1),
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
    fFracToKeepDownSamplqn(1.1),
    fApplyML(false),
    fConfigPath(""),
    fMLResponse(nullptr),
    fNMLBins(300),
    fMLOutputMin(0.85),
    fMLOutputMax(1.)
{
    // standard constructor
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
        case kLctopK0S:
            pdg=4122;
        break;
    }
    if(pdg==413) SetMassLimits((float)0.1,(float)0.2);
    else if(pdg==4122) SetMassLimits((float)0.25,pdg); //Default from Lc task
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
        case kLctopK0S:
            DefineOutput(2,AliRDHFCutsLc2V0bachelor::Class()); //Cut object for Lc->pK0S
    }
}

//________________________________________________________________________
AliAnalysisTaskSECharmHadronvn::~AliAnalysisTaskSECharmHadronvn()
{
    // Destructor
    if(fOutput && !fOutput->IsOwner()){
        delete fHistNEvents;
        delete fHistCandVsCent;
        delete fHistPercqnVsqnVsCentr;
        delete fHistNtrklVsqnVsCentr;
        delete fHistMassPtPhiqnCentr;

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

    if(fApplyML && fMLResponse)
        delete fMLResponse;
}

//_________________________________________________________________
void  AliAnalysisTaskSECharmHadronvn::SetMassLimits(float range, int pdg){
    // Set limits for mass spectra plots
    float mass=0;
    int abspdg=TMath::Abs(pdg);
    mass=TDatabasePDG::Instance()->GetParticle(abspdg)->Mass();
    fUpmasslimit = mass+range;
    fLowmasslimit = mass-range;
}

//_________________________________________________________________
void  AliAnalysisTaskSECharmHadronvn::SetMassLimits(float lowlimit, float uplimit){
    // Set limits for mass spectra plots
    if(uplimit>lowlimit) {
        fUpmasslimit = uplimit;
        fLowmasslimit = lowlimit;
    }
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronvn::LocalInit()
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
        case kLctopK0S:
            {
               AliRDHFCutsLc2V0bachelor* copycut = new AliRDHFCutsLc2V0bachelor(*(static_cast<AliRDHFCutsLc2V0bachelor*>(fRDCuts)));
               PostData(2,copycut);
            }
        break;
    }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronvn::UserCreateOutputObjects()
{
    // Create the output container
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("MainOutput");

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
    else if(fDecChannel==4) massaxisname = "#it{M}(pK^{0}_{S}) (GeV/#it{c}^{2})";

    int naxes=kVarForSparse;
    if(!fApplyML)
        naxes--;

    int nbins[kVarForSparse]     = {fNMassBins, nptbins, ndeltaphibins, nfphibins, nfphibins, nphibins, ncentbins, nNtrkBins, nqnbins, fNMLBins};
    double xmin[kVarForSparse]   = {fLowmasslimit, ptmin, mindeltaphi, fphimin, fphimin, phimin, fMinCentr, Ntrkmin, qnmin, fMLOutputMin};
    double xmax[kVarForSparse]   = {fUpmasslimit, ptmax, maxdeltaphi, fphimax, fphimax, phimax, fMaxCentr, Ntrkmax, qnmax, fMLOutputMax};
    TString axTit[kVarForSparse] = {massaxisname, "#it{p}_{T} (GeV/#it{c})", deltaphiname, nfphiname1, nfphiname2, "#varphi_{D}", "Centrality (%)", "#it{N}_{tracklets}", qnaxisnamefill, "ML response"};

    fHistMassPtPhiqnCentr = new THnSparseF("fHistMassPtPhiqnCentr",Form("InvMass vs. #it{p}_{T} vs. %s vs. centr vs. #it{q}_{%d} ",deltaphiname.Data(),fHarmonic),naxes,nbins,xmin,xmax);

    for(int iAx=0; iAx<naxes; iAx++)
        fHistMassPtPhiqnCentr->GetAxis(iAx)->SetTitle(axTit[iAx].Data());

    fOutput->Add(fHistMassPtPhiqnCentr);

    //ML model
    if(fApplyML) {
        switch(fDecChannel) {
            case kDplustoKpipi:
                {
                    fMLResponse = new AliHFMLResponseDplustoKpipi("DplustoKpipiMLResponse", "DplustoKpipiMLResponse", fConfigPath.Data());
                    fMLResponse->MLResponseInit();
                }
            break;
            case kD0toKpi:
                {
                    AliFatal("ML application for D0 meson not yet implemented! Exit");
                }
            break;
            case kDstartoKpipi:
                {
                    AliFatal("ML application for Dstar meson not yet implemented! Exit");
                }
            break;
            case kDstoKKpi:
                {
                    fMLResponse = new AliHFMLResponseDstoKKpi("DstoKKpiMLResponse", "DstoKKpiMLResponse", fConfigPath.Data());
                    fMLResponse->MLResponseInit();
                }
            break;
            case kLctopK0S:
               {
                    AliFatal("ML application for Lc->pK0S not yet implemented! Exit");
               }
            break;
        }
    }

    PostData(1,fOutput);
    return;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronvn::UserExec(Option_t */*option*/)
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

    TClonesArray *arrayProng = nullptr, *arrayD0toKpi = nullptr; *arrayV0daugh = nullptr;
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
                case kLctopK0S:
                    absPdgMom=4122;
                    nDau=3;
                    arrayProng = (TClonesArray*)aodFromExt->GetList()->FindObject("CascadesHF");
                    array
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
            case kLctopK0S:
                absPdgMom=4122;
                nDau=3; //!TODO check if 3 is correct here or if array of V0daughters needed
                arrayProng=(TClonesArray*)fAOD->GetList()->FindObject("CascadesHF");
            break;
        }
    }

    if(!fAOD || !arrayProng || (!arrayD0toKpi && fDecChannel==kDstartoKpipi)) {
        AliError("AliAnalysisTaskSECharmHadronvn::UserExec:Branch not found!\n");
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
    AliAODv0 *dV0 = nullptr;
    int nCand = arrayProng->GetEntriesFast();
    bool alreadyLooped = false;
    vector<int> isSelByPrevLoop;
    vector<double> MLPred[2];
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

                    if(fDecChannel==kDstartoKpipi) {
                        if(d->GetIsFilled()<1)
                            dD0 = dynamic_cast<AliAODRecoDecayHF2Prong*>(arrayD0toKpi->At(d->GetProngID(1)));
                        else
                            dD0 = (dynamic_cast<AliAODRecoCascadeHF*>(d))->Get2Prong();
                    }
                    if(fDecChannel==kLctopK0S) {
                            dV0 = (dynamic_cast<AliAODRecoCascadeHF*>(d))->Get2Prong();
                    }
                    double modelPred[2] = {-1., -1.};
                    int isSel = IsCandidateSelected(d, nDau, absPdgMom, vHF, dD0, dV0, modelPred);
                    isSelByPrevLoop.push_back(isSel);
                    MLPred[0].push_back(modelPred[0]);
                    MLPred[1].push_back(modelPred[1]);
                    if(!isSel) continue;

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

        if(fDecChannel==kDstartoKpipi) {
            if(d->GetIsFilled()<1)
                dD0 = dynamic_cast<AliAODRecoDecayHF2Prong*>(arrayD0toKpi->At(d->GetProngID(1)));
            else
                dD0 = (dynamic_cast<AliAODRecoCascadeHF*>(d))->Get2Prong();
        }

        int isSelected = 0;
        double modelPred[2] = {-1., -1.};
        if(alreadyLooped) {//check already here to avoid application of cuts twice
            isSelected = isSelByPrevLoop[iCand];
            modelPred[0] = MLPred[0][iCand];
            modelPred[1] = MLPred[1][iCand];
        }
        else {
            isSelected = IsCandidateSelected(d, nDau, absPdgMom, vHF, dD0, dV0, modelPred);
        }
        if(!isSelected) continue;

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
                double sparsearray[10] = {invMass[0],ptD,vnfunc,phifunc1,phifunc2,phiD,evCentr,static_cast<double>(tracklets),candpercqn,modelPred[0]};
                fHistMassPtPhiqnCentr->Fill(sparsearray);
                break;
            }
            case kD0toKpi:
            {
                if(isSelected==1 || isSelected==3) {
                    double sparsearray[10] = {invMass[0],ptD,vnfunc,phifunc1,phifunc2,phiD,evCentr,static_cast<double>(tracklets),candpercqn,modelPred[0]};
                    fHistMassPtPhiqnCentr->Fill(sparsearray);
                }
                if(isSelected==2 || isSelected==3) {
                    double sparsearray[10] = {invMass[1],ptD,vnfunc,phifunc1,phifunc2,phiD,evCentr,static_cast<double>(tracklets),candpercqn,modelPred[1]};
                    fHistMassPtPhiqnCentr->Fill(sparsearray);
                }
                break;
            }
            case kDstartoKpipi:
            {
                double sparsearray[10] = {invMass[0],ptD,vnfunc,phifunc1,phifunc2,phiD,evCentr,static_cast<double>(tracklets),candpercqn,modelPred[0]};
                fHistMassPtPhiqnCentr->Fill(sparsearray);
                break;
            }
            case kDstoKKpi:
            {
                if(isSelected&4) {
                    double sparsearray[10] = {invMass[0],ptD,vnfunc,phifunc1,phifunc2,phiD,evCentr,static_cast<double>(tracklets),candpercqn,modelPred[0]};
                    fHistMassPtPhiqnCentr->Fill(sparsearray);
                }
                if(isSelected&8) {
                    double sparsearray[10] = {invMass[1],ptD,vnfunc,phifunc1,phifunc2,phiD,evCentr,static_cast<double>(tracklets),candpercqn,modelPred[1]};
                    fHistMassPtPhiqnCentr->Fill(sparsearray);
                }
                break;
            }
            case kLctopK0S:
            {
                //!TODO check selection value
            }
        }
    }

    fHistCandVsCent->Fill(evCentr,nSelCand);

    trackstoremove.clear();
    isSelByPrevLoop.clear();
    MLPred[0].clear();
    MLPred[1].clear();

    delete vHF;
    vHF = nullptr;
    if(!isHandlerFound) { // if not found in the tender task, allocated memory with new --> to be deleted
        delete HFQnVectorHandler;
        HFQnVectorHandler = nullptr;
    }

    PostData(1,fOutput);
    return;
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronvn::SetQnVectorDetConf(int detconf)
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
double AliAnalysisTaskSECharmHadronvn::GetPhiInRange(double phi)
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
double AliAnalysisTaskSECharmHadronvn::GetDeltaPsiSubInRange(double psi1, double psi2)
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
void AliAnalysisTaskSECharmHadronvn::CalculateInvMasses(AliAODRecoDecayHF* d,float*& masses,int& nmasses)
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
        case kLctopK0S:
            {  //!TODO: check inv mass calculations
                nmasses=1; 
                masses = new float[nmasses];
                masses[0]=((AliAODRecoCascadeHF*)d)->InvMassLctoK0sP();
            }
        break;
    }
}

//________________________________________________________________________
void AliAnalysisTaskSECharmHadronvn::GetMainQnVectorInfo(double &mainPsin, double &mainMultQn, double mainQn[2], double &SubAPsin, double &SubAMultQn, double SubAQn[2], double &SubBPsin, double &SubBMultQn, double SubBQn[2], AliHFQnVectorHandler* HFQnVectorHandler)
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
void AliAnalysisTaskSECharmHadronvn::GetDaughterTracksToRemove(AliAODRecoDecayHF* d, int nDau, vector<AliAODTrack*> &trackstoremove)
{
 
    if (fDecChannel == kDstartoKpipi) {
        AliAODRecoDecayHF2Prong *dauD0 = dynamic_cast<AliAODRecoCascadeHF*>(d)->Get2Prong();
        trackstoremove.push_back(dynamic_cast<AliAODTrack*>(dauD0->GetDaughter(0)));
        trackstoremove.push_back(dynamic_cast<AliAODTrack*>(dauD0->GetDaughter(1)));
        if(fRemoveSoftPion) trackstoremove.push_back(dynamic_cast<AliAODRecoCascadeHF*>(d)->GetBachelor());
    }
    else if (fDecChannel == kLctopK0S) {
        AliAODv0 *dauV0 = dynamic_cast<AliAODRecoCascadeHF*>(d)->Getv0();
        trackstoremove.push_back(dynamic_cast<AliAODTrack*>(dauV0->GetDaughter(0)));
        trackstoremove.push_back(dynamic_cast<AliAODTrack*>(dauV0->GetDaughter(1)));
        trackstoremove.push_back(dynamic_cast<AliAODRecoCascadeHF*>(d)->GetBachelor());  
    }
    else {
        for(int iDau=0; iDau<nDau; iDau++) {
            AliAODTrack* dau = dynamic_cast<AliAODTrack*>(d->GetDaughter(iDau));
            trackstoremove.push_back(dau);
        }
    }
}

//________________________________________________________________________
int AliAnalysisTaskSECharmHadronvn::IsCandidateSelected(AliAODRecoDecayHF *&d, int nDau, int absPdgMom, AliAnalysisVertexingHF *vHF, AliAODRecoDecayHF2Prong *dD0, AliAODv0 *dV0, double modelPred[2]) {

    if(!d || !vHF || (fDecChannel==kDstartoKpipi && !dD0) || (fDecChannel==kLctopK0S && !dV0)) return false;
 
    //Preselection to speed up task
    TObjArray arrDauTracks(nDau);
    AliAODTrack *track = nullptr;
    
    if (fDecChannel == kDstartoKpipi){
        for(int iDau=0; iDau<nDau; iDau++){
            if(iDau == 0)
                track=vHF->GetProng(fAOD,d,iDau); //soft pion
            else
                track=vHF->GetProng(fAOD,dD0,iDau-1); //D0 daughters
            arrDauTracks.AddAt(track,iDau);
        }
    } else if (fDecChannel == kLctopK0S) {
         for(int iDau=0; iDau<nDau; iDau++) {
            if(iDau == 0)
               track = vHF->GetProng(fAOD,d,iDau); //bachelor
            else
               track = vHF->GetProng(fAOD,dV0,iDau-1);   //K0S daughters
            arrDauTracks.AddAt(track,iDau);
         }
    } else {
        for(int iDau=0; iDau<nDau; iDau++){
            AliAODTrack *track = vHF->GetProng(fAOD,d,iDau);
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
        case kLctopK0S:
            isSelBit = d->HasSelectionBit(AliRDHFCuts::kLctoV0Cuts);
            if(!isSelBit || !vHF->FillRecoCasc(fAOD,((AliAODRecoCascadeHF*)d),true)) return 0;
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
    if(fApplyML) {
        AliAODPidHF* pidHF = fRDCuts->GetPidHF();
        bool isMLsel = true;

        switch(fDecChannel) {
            case kDplustoKpipi:
            {
                isMLsel = fMLResponse->IsSelected(modelPred[0], d, fAOD->GetMagneticField(), pidHF);
                if(!isMLsel)
                    isSelected = 0;
                break;
            }
            case kD0toKpi:
            {
                AliWarning("ML application for D0 meson not yet implemented! Exit");
                break;
            }
            case kDstartoKpipi:
            {
                AliWarning("ML application for Dstar meson not yet implemented! Exit");
                break;
            }
            case kDstoKKpi:
            {
                if(isSelected&4) {
                    isMLsel = fMLResponse->IsSelected(modelPred[0], d, fAOD->GetMagneticField(), pidHF, 0);
                    if(!isMLsel)
                        isSelected &= ~4; //switch off bit ok KKpi hypo
                }
                if(isSelected&8) {
                    isMLsel = fMLResponse->IsSelected(modelPred[1], d, fAOD->GetMagneticField(), pidHF, 1);
                    if(!isMLsel)
                        isSelected &= ~8; //switch off bit ok piKK hypo
                }
                break;
            }
           case kLctopK0S:
           {
               AliWarning("ML application for Lc->pK0S not yet implemented! exit");
               break;
           }
        }
    }

    if(!isSelected) return 0;
    if(fDecChannel==kDstoKKpi)
        if(!(isSelected&4) && !(isSelected&8)) return 0;

    return isSelected;
}

//________________________________________________________________________
bool AliAnalysisTaskSECharmHadronvn::LoadSplinesForqnPercentile()
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
