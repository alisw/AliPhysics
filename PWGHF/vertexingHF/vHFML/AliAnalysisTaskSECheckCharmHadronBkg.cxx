/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// \class AliAnalysisTaskSECheckCharmHadronBkg
// \brief Analysis task to study correlated backgrounds in D2H analyses
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch

#include "AliAnalysisTaskSECheckCharmHadronBkg.h"

#include <TClonesArray.h>
#include <TList.h>
#include <TString.h>
#include <TDatabasePDG.h>

#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliVertexingHFUtils.h"
#include "AliHFMLResponseDplustoKpipi.h"
#include "AliAODPidHF.h"

#include <vector>
#include <algorithm>

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSECheckCharmHadronBkg);
/// \endcond

//________________________________________________________________________
AliAnalysisTaskSECheckCharmHadronBkg::AliAnalysisTaskSECheckCharmHadronBkg() : AliAnalysisTaskSE()
{
    /// Standrd constructor
}

//________________________________________________________________________
AliAnalysisTaskSECheckCharmHadronBkg::AliAnalysisTaskSECheckCharmHadronBkg(const char *name, AliRDHFCutsDplustoKpipi *dpluscutsana) : 
AliAnalysisTaskSE(name),
fRDCutsAnalysis(dpluscutsana)
{
    /// Standrd constructor

    // Output slot #1 writes into a TList container
    DefineOutput(1, TList::Class()); //My private output
    // Output slot #2 writes cut to private output
    //  DefineOutput(2,AliRDHFCutsDplustoKpipi::Class());
    DefineOutput(2, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskSECheckCharmHadronBkg::~AliAnalysisTaskSECheckCharmHadronBkg()
{
    /// Destructor
    if (fOutput)
        delete fOutput;
    if (fListCuts)
        delete fListCuts;
    if (fRDCutsAnalysis)
        delete fRDCutsAnalysis;

    if (fApplyML && fMLResponse)
        delete fMLResponse;
}

//__________________________________________
void AliAnalysisTaskSECheckCharmHadronBkg::Init()
{

    /// Initialization
    AliRDHFCutsDplustoKpipi *analysis = new AliRDHFCutsDplustoKpipi(*fRDCutsAnalysis);
    if (analysis)
    {
        analysis->SetName("AnalysisCuts");
        float *PtLims = (float *)analysis->GetPtBinLimits();
        int lastpt = analysis->GetNPtBins();
        fPtMin = PtLims[0];
        fPtMax = PtLims[lastpt];
    }

    double PtRange = fPtMax - fPtMin;
    fNPtBins = PtRange / fPtBinWidth;
    fPtMax = fPtMin + fNPtBins * fPtBinWidth;

    double MassRange = fUpMassLimit - fLowMassLimit;
    fNMassBins = MassRange / fMassBinWidth;
    fUpMassLimit = fLowMassLimit + fNMassBins * fMassBinWidth;

    fListCuts = new TList();
    fListCuts->SetOwner();
    fListCuts->Add(analysis);

    PostData(2, fListCuts);

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSECheckCharmHadronBkg::UserCreateOutputObjects()
{
    /// Create the output container
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("OutputHistos");

    fHistNEvents = new TH1F("fHistNEvents", "number of events ", 10, -0.5, 9.5);
    fHistNEvents->GetXaxis()->SetBinLabel(1, "nEvents read");
    fHistNEvents->GetXaxis()->SetBinLabel(2, "Rejected due to mismatch in trees");
    fHistNEvents->GetXaxis()->SetBinLabel(3, "nEvents with good AOD");
    fHistNEvents->GetXaxis()->SetBinLabel(4, "Rejected due to trigger");
    fHistNEvents->GetXaxis()->SetBinLabel(5, "Rejected due to vertex reco");
    fHistNEvents->GetXaxis()->SetBinLabel(6, "Rejected due to pileup");
    fHistNEvents->GetXaxis()->SetBinLabel(7, "Rejected due to centrality");
    fHistNEvents->GetXaxis()->SetBinLabel(8, "Rejected due to vtxz");
    fHistNEvents->GetXaxis()->SetBinLabel(9, "Rejected due to Physics Sel");
    fHistNEvents->GetXaxis()->SetBinLabel(10, "nEvents accepted");
    fHistNEvents->GetXaxis()->SetNdivisions(1, false);
    fHistNEvents->SetMinimum(0);
    fOutput->Add(fHistNEvents);

    fHistCentrality[0] = new TH2F("hCentrMult", "centrality", 100, 0.5, 30000.5, 40, 0., 100.);
    fHistCentrality[1] = new TH2F("hCentrMult(selectedCent)", "centrality(selectedCent)", 100, 0.5, 30000.5, 40, 0., 100.);
    fHistCentrality[2] = new TH2F("hCentrMult(OutofCent)", "centrality(OutofCent)", 100, 0.5, 30000.5, 40, 0., 100.);
    for (int iHisto = 0; iHisto < 3; iHisto++)
    {
        fHistCentrality[iHisto]->Sumw2();
        fOutput->Add(fHistCentrality[iHisto]);
    }
    CreateSparse();

    //Loading of ML models
    if (fApplyML)
    {
        fMLResponse = new AliHFMLResponseDplustoKpipi("DplustoKpipiMLResponse", "DplustoKpipiMLResponse", fConfigPath.Data());
        fMLResponse->MLResponseInit();
    }

    PostData(1, fOutput);

    return;
}

//________________________________________________________________________
void AliAnalysisTaskSECheckCharmHadronBkg::UserExec(Option_t * /*option*/)
{

    fAOD = dynamic_cast<AliAODEvent *>(InputEvent());
    fHistNEvents->Fill(0); // count event

    if (fAODProtection >= 0)
    {
        //   Protection against different number of events in the AOD and deltaAOD
        //   In case of discrepancy the event is rejected.
        int matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
        if (matchingAODdeltaAODlevel < 0 || (matchingAODdeltaAODlevel == 0 && fAODProtection == 1))
        {
            // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
            fHistNEvents->Fill(1);
            return;
        }
    }

    TClonesArray *array3Prong = 0;
    if (!fAOD && AODEvent() && IsStandardAOD())
    {
        fAOD = dynamic_cast<AliAODEvent *>(AODEvent());
        AliAODHandler *aodHandler = (AliAODHandler *)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
        if (aodHandler->GetExtensions())
        {
            AliAODExtension *ext = (AliAODExtension *)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
            AliAODEvent *aodFromExt = ext->GetAOD();
            array3Prong = (TClonesArray *)aodFromExt->GetList()->FindObject("Charm3Prong");
        }
    }
    else if (fAOD)
        array3Prong = (TClonesArray *)fAOD->GetList()->FindObject("Charm3Prong");

    if (!fAOD || !array3Prong)
    {
        printf("AliAnalysisTaskSECheckCharmHadronBkg::UserExec: Charm3Prong branch not found!\n");
        return;
    }

    // fix for temporary bug in ESDfilter
    // the AODs with null vertex pointer didn't pass the PhysSel
    if (!fAOD->GetPrimaryVertex() || TMath::Abs(fAOD->GetMagneticField()) < 0.001)
        return;
    fHistNEvents->Fill(2);

    //Event selection
    bool isEvSel = fRDCutsAnalysis->IsEventSelected(fAOD);
    int ntracks = fAOD->GetNumberOfTracks();
    double evCentr = 0;
    if (fRDCutsAnalysis->GetUseCentrality() > 0)
        evCentr = fRDCutsAnalysis->GetCentrality(fAOD);
    fHistCentrality[0]->Fill(ntracks, evCentr);
    if (fRDCutsAnalysis->GetWhyRejection() == 5)
        fHistNEvents->Fill(3);
    if (!isEvSel && fRDCutsAnalysis->GetWhyRejection() == 0)
        fHistNEvents->Fill(4);
    if (fRDCutsAnalysis->GetWhyRejection() == 1)
        fHistNEvents->Fill(5);
    if (fRDCutsAnalysis->GetWhyRejection() == 2)
    {
        fHistNEvents->Fill(6);
        fHistCentrality[2]->Fill(ntracks, evCentr);
    }
    if (fRDCutsAnalysis->GetWhyRejection() == 6)
        fHistNEvents->Fill(7);
    if (fRDCutsAnalysis->GetWhyRejection() == 7)
        fHistNEvents->Fill(8);

    TClonesArray *arrayMC = (TClonesArray *)fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if (!arrayMC)
    {
        printf("AliAnalysisTaskSECheckCharmHadronBkg::UserExec: MC particles branch not found!\n");
        return;
    }
    AliAODMCHeader *mcHeader = (AliAODMCHeader *)fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader)
    {
        printf("AliAnalysisTaskSECheckCharmHadronBkg::UserExec: MC header branch not found!\n");
        return;
    }
    int runNumber = fAOD->GetRunNumber();

    if (fAOD->GetTriggerMask() == 0 && (runNumber >= 195344 && runNumber <= 195677))
        return;
    if (fRDCutsAnalysis->GetUseCentrality() > 0 && fRDCutsAnalysis->IsEventSelectedInCentrality(fAOD) != 0)
        return;

    // Post the data already here
    PostData(1, fOutput);
    if (!isEvSel)
        return;

    fHistCentrality[1]->Fill(ntracks, evCentr);
    fHistNEvents->Fill(9);

    // AOD primary vertex
    AliAODVertex *trksvtx = (AliAODVertex *)fAOD->GetPrimaryVertex();
    int n3Prong = array3Prong->GetEntriesFast();

    // vHF object is needed to call the method that refills the missing info of the candidates
    // if they have been deleted in dAOD reconstruction phase
    // in order to reduce the size of the file
    AliAnalysisVertexingHF vHF;
    for (int i3Prong = 0; i3Prong < n3Prong; i3Prong++)
    {
        AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong *)array3Prong->UncheckedAt(i3Prong);
        if (fUseBit && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts))
            continue;
        if (!(vHF.FillRecoCand(fAOD, d)))
            continue;

        int passTopolAndPIDCuts = fRDCutsAnalysis->IsSelected(d, AliRDHFCuts::kAll, fAOD);

        if (!fRDCutsAnalysis->GetIsSelectedCuts())
            continue;
        if (!passTopolAndPIDCuts)
            continue;

        //variables for ML application
        AliAODPidHF *Pid_HF = nullptr;
        std::vector<double> modelPred{};
        bool isMLsel = true;
        if (fApplyML)
        {
            Pid_HF = fRDCutsAnalysis->GetPidHF();
            isMLsel = fMLResponse->IsSelectedMultiClass(modelPred, d, fAOD->GetMagneticField(), Pid_HF);
        }

        if(!isMLsel)
            continue;

        bool unsetvtx = false;
        if (!d->GetOwnPrimaryVtx())
        {
            d->SetOwnPrimaryVtx(trksvtx);
            unsetvtx = true;
        }

        bool recVtx = false;
        AliAODVertex *origownvtx = nullptr;
        if (fRDCutsAnalysis->GetIsPrimaryWithoutDaughters())
        {
            if (d->GetOwnPrimaryVtx())
                origownvtx = new AliAODVertex(*d->GetOwnPrimaryVtx());
            if (fRDCutsAnalysis->RecalcOwnPrimaryVtx(d, fAOD))
                recVtx = true;
            else
                fRDCutsAnalysis->CleanOwnPrimaryVtx(d, fAOD, origownvtx);
        }
        int pdgDgDplustoKpipi[3] = {321, 211, 211};
        int labMother = -1;
        int MotherPDG = -1;
        int daulab[3] = {-1, -1, -1};
        int dauPDG[3] = {-1, -1, -1};
        double isDplusToKpipi = 0;
        int decay = 0;
        double labeldaus = -1;

        //check if the three daughter particles have the same mother
        labMother = d->MatchToMC(411, arrayMC, 3, pdgDgDplustoKpipi);
        if (labMother >= 0)
            isDplusToKpipi = 1;
        labMother = -1;
        for (int iTrk = 0; iTrk < 3; iTrk++)
        {
            AliAODTrack *trk = (AliAODTrack *)d->GetDaughter(iTrk);
            daulab[iTrk] = trk->GetLabel();
            daulab[iTrk] = TMath::Abs(daulab[iTrk]);
        }
        if (daulab[0] >= 0 && daulab[1] >= 0 && daulab[2] >= 0)
        {
            labeldaus = 1;
            labMother = CheckMother(arrayMC, daulab, dauPDG);
        }
        if (labMother >= 0)
        {
            AliAODMCParticle *mother = (AliAODMCParticle *)arrayMC->At(labMother);
            MotherPDG = TMath::Abs(mother->GetPdgCode());
            if (MotherPDG == 411)
            { //D+
                if (MatchToSpecie(dauPDG, 321, 211, 211))
                    decay = kDplusToKpipi;
                else if (MatchToSpecie(dauPDG, 321, 321, 211))
                    decay = kDplusToKKpi;
                else if (MatchToSpecie(dauPDG, 211, 211, 211))
                    decay = kDplusTopipipi;
                else if (MatchToSpecie(dauPDG, 321, 211, 11))
                    decay = kDplusToKpie;
                else if (MatchToSpecie(dauPDG, 321, 211, 13))
                    decay = kDplusToKpimu;
            }
            else if (MotherPDG == 413)
            { //D*+
                if (IsInDecayChain(421, labMother, arrayMC) && MatchToSpecie(dauPDG, 321, 211, 211))
                    decay = kDstarToKpipi;
            }
            else if (MotherPDG == 421)
            { //D0
                if (MatchToSpecie(dauPDG, 211, 211, 211))
                    decay = kDzeroTopipipiX;
                else if (MatchToSpecie(dauPDG, 321, 211, 211))
                    decay = kDzeroToKpipiX;
            }
            else if (MotherPDG == 431)
            { //Ds+
                if (MatchToSpecie(dauPDG, 321, 321, 211))
                    decay = kDsToKKpi;
                else if (MatchToSpecie(dauPDG, -321, 211, 211, false))
                    decay = kDsToKminuspipluspiplus;
                else if (MatchToSpecie(dauPDG, 321, 211, -211, false))
                    decay = kDsToKpluspipluspiminus;
            }
            else if (MotherPDG == 511)
            { //B0
                if (IsInDecayChain(411, labMother, arrayMC))
                    decay = kBzeroToDplusX;
                if (IsInDecayChain(421, labMother, arrayMC))
                    decay = kBzeroToDzeroX;
                if (IsInDecayChain(431, labMother, arrayMC))
                    decay = kBzeroToDsX;
                if (IsInDecayChain(4122, labMother, arrayMC))
                    decay = kBzeroToLambdacX;
                if (IsInDecayChain(443, labMother, arrayMC))
                    decay = kBzeroToJpsiX;
            }
            else if (MotherPDG == 521)
            { //B+
                if (IsInDecayChain(411, labMother, arrayMC))
                    decay = kBplusToDplusX;
                if (IsInDecayChain(421, labMother, arrayMC))
                    decay = kBplusToDzeroX;
                if (IsInDecayChain(431, labMother, arrayMC))
                    decay = kBplusToDsX;
                if (IsInDecayChain(4122, labMother, arrayMC))
                    decay = kBplusToLambdacX;
                if (IsInDecayChain(443, labMother, arrayMC))
                    decay = kBplusToJpsiX;
            }
            else if (MotherPDG == 531)
            { //Bs
                if (IsInDecayChain(431, labMother, arrayMC))
                    decay = kBsToDsX;
            }
            else if (MotherPDG == 4122)
            { //Lambdac
                if (MatchToSpecie(dauPDG, 2212, 321, 211))
                    decay = kLambdacTopKpi;
                if (MatchToSpecie(dauPDG, 2212, 211, 211))
                    decay = kLambdacToppipi;
            }
        }
        double rapid = d->YDplus();
        double pt = d->Pt();
        double invMass = d->InvMassDplus();
        bool isFidAcc = fRDCutsAnalysis->IsInFiducialAcceptance(pt, rapid);
        //fill the THnSparseF
        if (isFidAcc)
        {
            double sumabsdauPDG = TMath::Abs(dauPDG[0]) + TMath::Abs(dauPDG[1]) + TMath::Abs(dauPDG[2]);
            double sparsearray[7] = {invMass, pt, (double)decay, (double)MotherPDG, isDplusToKpipi, sumabsdauPDG, labeldaus};
            fMassVsPtVsOriginSparse->Fill(sparsearray);
        }

        if (recVtx)
            fRDCutsAnalysis->CleanOwnPrimaryVtx(d, fAOD, origownvtx);
        if (unsetvtx)
            d->UnsetOwnPrimaryVtx();
    }

    PostData(1, fOutput);
    PostData(2, fListCuts);
    return;
}

//________________________________________________________________________
void AliAnalysisTaskSECheckCharmHadronBkg::CreateSparse()
{
    /// Sparse for candidate-origin studies

    double PDGmin = -1.5;
    double PDGmax = 30000 - 0.5;
    int nPDGbins = PDGmax - PDGmin;

    double decaymin = -0.5;
    double decaymax = 27.5;
    int ndecaybins = decaymax - decaymin;

    TString axTit[7] = {"#it{M}_{K#pi#pi} (GeV/#it{c}^{2})", "#it{p}_{T} (GeV/#it{c})", "decay", "PDG mother particle", "isDtoKpipi", "sumdauPDG", "labdau"};

    int nbins[7] = {fNMassBins, fNPtBins, ndecaybins, nPDGbins, 2, 5999, 3};
    double xmin[7] = {fLowMassLimit, fPtMin, decaymin, PDGmin, -0.5, -2999.5, -1.5};
    double xmax[7] = {fUpMassLimit, fPtMax, decaymax, PDGmax, 1.5, 2999.5, 1.5};

    fMassVsPtVsOriginSparse = new THnSparseF("fMassVsPtVsOriginSparse", "Mass vs. origin", 7, nbins, xmin, xmax);

    for (int iAx = 0; iAx < 7; iAx++)
        fMassVsPtVsOriginSparse->GetAxis(iAx)->SetTitle(axTit[iAx].Data());
    fOutput->Add(fMassVsPtVsOriginSparse);
}

//________________________________________________________________________
int AliAnalysisTaskSECheckCharmHadronBkg::CheckMother(TClonesArray *mcArray, int daulab[3], int dauPDG[3])
{

    AliAODMCParticle *part1 = (AliAODMCParticle *)mcArray->At(daulab[0]);
    AliAODMCParticle *part2 = (AliAODMCParticle *)mcArray->At(daulab[1]);
    AliAODMCParticle *part3 = (AliAODMCParticle *)mcArray->At(daulab[2]);
    if (!part1 || !part2 || !part3)
        return -1;

    dauPDG[0] = TMath::Abs(part1->GetPdgCode());
    dauPDG[1] = TMath::Abs(part2->GetPdgCode());
    dauPDG[2] = TMath::Abs(part3->GetPdgCode());

    std::vector<int> labmother1;
    std::vector<int> labmother2;
    std::vector<int> labmother3;
    std::vector<int> PDGmother1;
    std::vector<int> PDGmother2;
    std::vector<int> PDGmother3;
    std::vector<int> commonlabel12;
    std::vector<int> commonlabel13;
    std::vector<int> commonlabel123;
    labmother1.push_back(part1->GetMother());
    labmother2.push_back(part2->GetMother());
    labmother3.push_back(part3->GetMother());
    int mothercounter1 = 0;
    int mothercounter2 = 0;
    int mothercounter3 = 0;
    while (labmother1[mothercounter1] >= 0)
    {
        AliAODMCParticle *mother1 = (AliAODMCParticle *)mcArray->At(labmother1[mothercounter1]);
        PDGmother1.push_back(mother1->GetPdgCode());
        if (mother1)
            labmother1.push_back(mother1->GetMother());
        mothercounter1++;
    }
    while (labmother2[mothercounter2] >= 0)
    {
        AliAODMCParticle *mother2 = (AliAODMCParticle *)mcArray->At(labmother2[mothercounter2]);
        PDGmother2.push_back(mother2->GetPdgCode());
        if (mother2)
            labmother2.push_back(mother2->GetMother());
        mothercounter2++;
    }
    while (labmother3[mothercounter3] >= 0)
    {
        AliAODMCParticle *mother3 = (AliAODMCParticle *)mcArray->At(labmother3[mothercounter3]);
        PDGmother3.push_back(mother3->GetPdgCode());
        if (mother3)
            labmother3.push_back(mother3->GetMother());
        mothercounter3++;
    }
    bool commonmother12 = false;
    bool commonmother13 = false;
    bool commonmother123 = false;
    for (size_t iMother1 = 0; iMother1 < labmother1.size() - 1; iMother1++)
    {
        for (size_t iMother2 = 0; iMother2 < labmother2.size() - 1; iMother2++)
        {
            if (labmother1[iMother1] == labmother2[iMother2] && isPDGcodeAcceptable(PDGmother1[iMother1]))
            {
                commonmother12 = true;
                commonlabel12.push_back(labmother1[iMother1]);
            }
        }
    }
    for (size_t iMother1 = 0; iMother1 < labmother1.size() - 1; iMother1++)
    {
        for (size_t iMother3 = 0; iMother3 < labmother3.size() - 1; iMother3++)
        {
            if (labmother1[iMother1] == labmother3[iMother3] && isPDGcodeAcceptable(PDGmother1[iMother1]))
            {
                commonmother13 = true;
                commonlabel13.push_back(labmother1[iMother1]);
            }
        }
    }
    if (commonlabel13.size() != 0 && commonlabel12.size() != 0)
    {
        for (size_t iLab12 = 0; iLab12 < commonlabel12.size(); iLab12++)
        {
            for (size_t iLab13 = 0; iLab13 < commonlabel13.size(); iLab13++)
            {
                if (commonlabel12[iLab12] == commonlabel13[iLab13])
                {
                    commonmother123 = true;
                    commonlabel123.push_back(commonlabel12[iLab12]);
                }
            }
        }
        int firstcommon = *max_element(commonlabel123.begin(), commonlabel123.end());
        if (commonmother123)
            return firstcommon;
    }
    return -1;
}

//________________________________________________________________________
bool AliAnalysisTaskSECheckCharmHadronBkg::isPDGcodeAcceptable(int pdg)
{

    pdg = TMath::Abs(pdg);
    if (pdg <= 8)
        return false; //quarks
    if (pdg >= 21 && pdg <= 37)
        return false; //gauge and higgs bosons
    if (pdg >= 30000)
        return false;
    return true;
}

//________________________________________________________________________
bool AliAnalysisTaskSECheckCharmHadronBkg::MatchToSpecie(int dauPDG[3], int pdg1, int pdg2, int pdg3, bool absPDGcode)
{

    double pdg[3];
    for (int iDau = 0; iDau < 3; iDau++)
    {
        if (absPDGcode)
            pdg[iDau] = TMath::Abs(dauPDG[iDau]);
        else
            pdg[iDau] = dauPDG[iDau];
    }

    if (pdg[0] == pdg1 && pdg[1] == pdg2 && pdg[2] == pdg3)
        return true;
    if (pdg[0] == pdg1 && pdg[2] == pdg2 && pdg[1] == pdg3)
        return true;
    if (pdg[1] == pdg1 && pdg[0] == pdg2 && pdg[2] == pdg3)
        return true;
    if (pdg[1] == pdg1 && pdg[2] == pdg2 && pdg[0] == pdg3)
        return true;
    if (pdg[2] == pdg1 && pdg[0] == pdg2 && pdg[1] == pdg3)
        return true;
    if (pdg[2] == pdg1 && pdg[1] == pdg2 && pdg[0] == pdg3)
        return true;

    //if not absolute PDG code -> check c.c.
    if (pdg[0] == -pdg1 && pdg[1] == -pdg2 && pdg[2] == -pdg3)
        return true;
    if (pdg[0] == -pdg1 && pdg[2] == -pdg2 && pdg[1] == -pdg3)
        return true;
    if (pdg[1] == -pdg1 && pdg[0] == -pdg2 && pdg[2] == -pdg3)
        return true;
    if (pdg[1] == -pdg1 && pdg[2] == -pdg2 && pdg[0] == -pdg3)
        return true;
    if (pdg[2] == -pdg1 && pdg[0] == -pdg2 && pdg[1] == -pdg3)
        return true;
    if (pdg[2] == -pdg1 && pdg[1] == -pdg2 && pdg[0] == -pdg3)
        return true;

    return false;
}

//________________________________________________________________________
bool AliAnalysisTaskSECheckCharmHadronBkg::IsInDecayChain(int abspdg, int motherlab, TClonesArray *mcArray)
{

    for (int iPart = 0; iPart < mcArray->GetEntriesFast(); iPart++)
    {
        AliAODMCParticle *part = (AliAODMCParticle *)mcArray->At(iPart);
        if (TMath::Abs(part->GetPdgCode()) == abspdg)
        {
            int lab = part->GetMother();
            AliAODMCParticle *mother = 0x0;
            if (lab >= 0)
                mother = (AliAODMCParticle *)mcArray->At(lab);
            while (lab != motherlab && lab >= 0)
            {
                lab = mother->GetMother();
                if (lab >= 0)
                    mother = (AliAODMCParticle *)mcArray->At(lab);
            }
            if (lab == motherlab)
                return true;
        }
    }
    return false;
}

//________________________________________________________________________
void AliAnalysisTaskSECheckCharmHadronBkg::Terminate(Option_t * /*option*/)
{
    /// Terminate analysis

    fOutput = dynamic_cast<TList *>(GetOutputData(1));
    if (!fOutput)
    {
        printf("ERROR: fOutput not available\n");
        return;
    }

    fHistNEvents = dynamic_cast<TH1F *>(fOutput->FindObject("fHistNEvents"));
    if (fHistNEvents)
        printf("Number of analyzed events = %d\n", (int)fHistNEvents->GetBinContent(10));
    else
    {
        printf("ERROR: fHistNEvents not available\n");
        return;
    }

    return;
}
