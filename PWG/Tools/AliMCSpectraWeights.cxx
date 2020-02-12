/**
 * @file AliMCSpectraWeights.cxx
 * @brief Class for re-weighting particle abundances in MC simulations
 * @author Patrick Huhn
 * @date 25/10/2019
 */
#include "AliMCEvent.h"
#include "AliMCSpectraWeights.h"
#include "AliStack.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH3F.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include <string>

#if defined(__CLING__)
#include <algorithm>
#include <array>
#include <memory>
#endif

int GetBinFromTH3(TH3F* h, std::array<float, 3> const & _values) {
    return h->FindBin(_values[0], _values[1], _values[2]);
}

void FillTH3WithValue(TH3F* h, std::array<float, 3> _xyzValue, float _value) {
    auto const _iBin = GetBinFromTH3(h, _xyzValue);
    h->SetBinContent(_iBin, _value);
    h->SetBinError(_iBin, 1e-30);
}

/**
 * @brief default constructor
 */
AliMCSpectraWeights::AliMCSpectraWeights()
    : fstCollisionSystem("pp"), fstFileMCSpectra(""), fstFilePublished(""),
      fstSavedObjName("fMCSpectraWeights"), fstSavedListName("dNdPt_ParCor"),
      fstPartTypes(0), fstCentralities(0), fBinsMultCent{0}, fBinsPt{0},
      fHistMCGenPrimTrackParticle(nullptr), fHistDataFractions(nullptr),
      fHistMCFractions(nullptr), fHistMCWeights(nullptr), fMCEvent(nullptr),
      fMultOrCent(0), fNPartTypes(6), fNCentralities(0),
      fbTaskStatus(AliMCSpectraWeights::TaskState::kAllEmpty),
      fFlag(AliMCSpectraWeights::SysFlag::kNominal), fUseMultiplicity(kTRUE) {}

/**
 *  @brief standard way for constuctor
 *  @param[in] collisionSystem string for selecting the used collision system.
 * Supported systems are: pp, pPb, PbPb and XeXe
 *  @param[in] stName string for naming the object. Not used anymore. Included
 * only for backward compability.
 *  @param[in] flag used for systematic variations.
 */
AliMCSpectraWeights::AliMCSpectraWeights(std::string collisionSystem,
                                         std::string stName,
                                         AliMCSpectraWeights::SysFlag flag)
    : fstCollisionSystem("pp"), fstFileMCSpectra(""), fstFilePublished(""),
      fstSavedObjName("fMCSpectraWeights"), fstSavedListName("dNdPt_ParCor"),
      fstPartTypes(), fstCentralities(), fBinsMultCent{0}, fBinsPt{0},
      fHistMCGenPrimTrackParticle(nullptr), fHistDataFractions(nullptr),
      fHistMCFractions(nullptr), fHistMCWeights(nullptr), fMCEvent(nullptr),
      fMultOrCent(0), fNPartTypes(6), fNCentralities(0),
      fbTaskStatus(AliMCSpectraWeights::TaskState::kAllEmpty), fFlag(flag),
      fUseMultiplicity(kTRUE) {
    fstCollisionSystem = collisionSystem;
    std::for_each(
        fstCollisionSystem.begin(), fstCollisionSystem.end(),
        [](char& c) { c = std::tolower(static_cast<unsigned char>(c)); });
    // setting uniform name
    if (fstCollisionSystem == "pp" || fstCollisionSystem == "p-p")
        fstCollisionSystem = "pp";
    else if (fstCollisionSystem == "ppb" || fstCollisionSystem == "p-pb")
        fstCollisionSystem = "ppb";
    else if (fstCollisionSystem == "pbpb" || fstCollisionSystem == "pb-pb")
        fstCollisionSystem = "pbpb";
    else if (fstCollisionSystem == "xexe" || fstCollisionSystem == "xe-xe")
        fstCollisionSystem = "xexe";
    else
        fstCollisionSystem = "pp";

    // set default Binning
    // pT binning
    fBinsPt = {0.0,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,   0.45, 0.5,
               0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85,  0.9,  0.95,
               1.0,  1.1,  1.2,  1.4,  1.6,  1.8,  2.0,   2.2,  2.4,
               2.6,  2.8,  3.0,  3.2,  3.6,  4.0,  5.0,   6.0,  8.0,
               10.0, 13.0, 20.0, 30.0, 50.0, 80.0, 100.0, 200.0};
    // multiplicity binning
    if (fstCollisionSystem.find("pp") != std::string::npos) {
        fNCentralities = 10;
        std::vector<std::string> tmpCent{};
        tmpCent.reserve(fNCentralities);
        for (int icent = 0; icent < fNCentralities; icent++) {
            tmpCent.push_back(Form("%d", icent));
        }
        fstCentralities = tmpCent;
        fBinsMultCent.reserve(51);
        for (int i = 0; i < 51; ++i) {
            fBinsMultCent.push_back(i);
        }
    } else if (fstCollisionSystem.find("ppb") != std::string::npos) {
        fBinsMultCent.reserve(301);
        for (int i = 0; i < 301; ++i) {
            fBinsMultCent.push_back(i);
        }
    } else if (fstCollisionSystem.find("pbpb") != std::string::npos) {
        fBinsMultCent.reserve(201);
        for (int i = 0; i < 201; ++i) {
            fBinsMultCent.push_back(i * 25);
        }
    } else if (fstCollisionSystem.find("xexe") != std::string::npos) {
        fBinsMultCent.reserve(201);
        for (int i = 0; i < 201; ++i) {
            fBinsMultCent.push_back(i * 25);
        }
    } else {
        fBinsMultCent.reserve(101);
        for (int i = 0; i < 101; ++i) {
            fBinsMultCent.push_back(i);
        }
    }
    if (fstCollisionSystem.find("pp") == std::string::npos) {
        std::vector<std::string> cent{"MB",    "c0005", "c0510", "c0010",
                                      "c1020", "c2040", "c4060", "c6080"};
        fNCentralities = 8;
        fstCentralities.clear();
        fstCentralities = cent;
    }
    const std::vector<std::string> tmpPart{"Pion",       "Proton",    "Kaon",
                                           "SigmaMinus", "SigmaPlus", "Rest"};
    fstPartTypes = tmpPart;
    fbTaskStatus = AliMCSpectraWeights::TaskState::kAllEmpty;
    fstFilePublished =
        "alien:///alice/cern.ch/user/p/phuhn/AllPublishedFractions.root";
}

AliMCSpectraWeights::~AliMCSpectraWeights() {
    //    if(fHistMCGenPrimTrackParticle) delete fHistMCGenPrimTrackParticle; //
    //    is put into train output in AnalysisTask if(fHistDataFractions) delete
    //    fHistDataFractions; if(fHistMCFractions) delete fHistMCFractions;
    //    if(fHistMCWeights) delete fHistMCWeights;
    if (fMCEvent)
        fMCEvent = 0;
}

void AliMCSpectraWeights::SetBinsPt(std::vector<double> bins) {
    fBinsPt.clear();
    fBinsPt.reserve(bins.size());
    for (auto& x : bins) {
        fBinsPt.push_back(x);
    }
}

void AliMCSpectraWeights::SetBinsMultCent(std::vector<double> bins) {
    fBinsMultCent.clear();
    fBinsMultCent.reserve(bins.size());
    for (auto& x : bins) {
        fBinsMultCent.push_back(x);
    }
}

/**
 * @brief Initialisation of object
 *
 * To be called after all settings were made. Will create all internal
 * histograms, load the input from previous train outputs if available, load the
 * expert input from alien, calculate weight factors. Set proper internal status
 * depending on what is avaiable
 */
void AliMCSpectraWeights::Init() {
    // Histograms
    AliMCSpectraWeights::InitHistos();

    // Loading output from previous train
    if (fstFileMCSpectra.length() > 5) // *.root
    {
        // printf("AliMCSpectraWeights:: Opening %s \n",
        // fstFileMCSpectra.Data());
        TFile* fInput = TFile::Open(fstFileMCSpectra.c_str());
        if (fInput) {
            if (fInput->GetNkeys() != 1) {
                if (!fInput->Get("dNdPt_ParCor"))
                    std::cerr << "AliMCSpectraWeights::WARNING: more than 1 "
                                 "list in the streamed file; please specify; "
                                 "using 1st list;\n\n";
            } else
                fstSavedListName = fInput->GetListOfKeys()->At(0)->GetName();
            // printf("AliMCSpectraWeights:: Loading %s from list %s\n",
            // fstSavedObjName.Data(), fstSavedListName.Data());
            TList* listMC = (TList*)fInput->Get(fstSavedListName.c_str());
            if (!listMC) {
                std::cerr << "AliMCSpectraWeights::ERROR: could not load list "
                             "in streamed "
                             "file\n";
            } else {
                auto tmpHist =
                    (TH3F*)listMC->FindObject("fHistMCGenPrimTrackParticle");
                if (!tmpHist) {
                    std::cerr << "AliMCSpectraWeights::WARNING: Couln't get "
                                 "fHistMCGenPrimTrackParticle\n";
                    //                    AliMCSpectraWeights* inWeights =
                    //                        (AliMCSpectraWeights*)listMC->FindObject(
                    //                            fstSavedObjName.c_str());
                    //                    if (!inWeights) {
                    //                        std::cerr <<
                    //                        "AliMCSpectraWeights::ERROR:
                    //                        Couln't load "
                    //                                     "object from "
                    //                                     "previous output";
                    //                        return;
                    //                    }
                    //                    if
                    //                    (AliMCSpectraWeights::LoadFromAliMCSpectraWeight(
                    //                            inWeights))
                    //                        fbTaskStatus =
                    //                            AliMCSpectraWeights::TaskState::kMCSpectraObtained;
                    //                    else
                    //                        return;
                }
                fHistMCGenPrimTrackParticle =
                    (TH3F*)tmpHist->Clone("fHistMCGenPrimTrackParticle_prev");
                fHistMCGenPrimTrackParticle->SetDirectory(0);
                if (fHistMCGenPrimTrackParticle->GetEntries() > 0) {
                    std::cerr << "Previous train has mc tracks\n";
                    fbTaskStatus =
                        AliMCSpectraWeights::TaskState::kMCSpectraObtained;
                }
                //                delete tmpHist; // gets deleted with
                //                fInput->Close() delete listMC;
            }
            fInput->Close();
            delete fInput;
        }
    }

    // Loading measured fractions
    if (fbTaskStatus == AliMCSpectraWeights::TaskState::kMCSpectraObtained) {
        AliMCSpectraWeights::LoadMeasuredFractions();
        fbTaskStatus = AliMCSpectraWeights::TaskState::kDataFractionLoaded;
    }

    // Calculating weight factors
    if (fbTaskStatus == AliMCSpectraWeights::TaskState::kDataFractionLoaded) {
        if (AliMCSpectraWeights::CalculateMCWeights()) {
            fbTaskStatus = AliMCSpectraWeights::TaskState::kMCWeightCalculated;
        }
    }

    std::cerr << "AliMCSpectraWeights::INFO: Init finished with status "
              << fbTaskStatus << std::endl;
}

/**
 *  @brief Create all internal histograms
 */
void AliMCSpectraWeights::InitHistos() {
    // Initalizing histograms
    // histogram charged patricles pt:multcent:type
    std::array<float, 7> partArray{};
    std::array<float, 6> partArrayDATA{};
    for (int i = 0; i < 7; ++i) {
        partArray[i] = -0.5 + i;
        if (i < 6)
            partArrayDATA[i] = -0.5 + i;
    }

    fHistMCGenPrimTrackParticle =
        new TH3F("fHistMCGenPrimTrackParticle",
                 "histogram for charged particle composition;#it{p}_{T} "
                 "(GeV/#it{c});multiplicity or centrality;Particle type",
                 static_cast<int>(fBinsPt.size()) - 1,
                 static_cast<float*>(fBinsPt.data()),
                 static_cast<int>(fBinsMultCent.size()) - 1,
                 static_cast<float*>(fBinsMultCent.data()),
                 static_cast<int>(partArray.size()) - 1,
                 static_cast<float*>(partArray.data()));
    fHistMCGenPrimTrackParticle->Sumw2();

    fHistDataFractions =
        new TH3F("fHistDataFractions",
                 "DATA fractions histogram;#it{p}_{T} "
                 "(GeV/#it{c});multiplicity or centrality;Particle type",
                 static_cast<int>(fBinsPt.size()) - 1,
                 static_cast<float*>(fBinsPt.data()),
                 static_cast<int>(fBinsMultCent.size()) - 1,
                 static_cast<float*>(fBinsMultCent.data()),
                 static_cast<int>(partArrayDATA.size()) - 1,
                 static_cast<float*>(partArrayDATA.data()));
    fHistDataFractions->Sumw2();

    fHistMCFractions =
        new TH3F("fHistMCFractions",
                 "MC fractions histogram;#it{p}_{T} (GeV/#it{c});multiplicity "
                 "or centrality;Particle type",
                 static_cast<int>(fBinsPt.size()) - 1,
                 static_cast<float*>(fBinsPt.data()),
                 static_cast<int>(fBinsMultCent.size()) - 1,
                 static_cast<float*>(fBinsMultCent.data()),
                 static_cast<int>(partArray.size()) - 1,
                 static_cast<float*>(partArray.data()));
    fHistMCFractions->Sumw2();

    fHistMCWeights = new TH3F(
        "fHistMCWeights",
        "MC weight histogram for charged particle composition;#it{p}_{T} "
        "(GeV/#it{c});multiplicity or centrality;Particle type",
        static_cast<int>(fBinsPt.size()) - 1,
        static_cast<float*>(fBinsPt.data()),
        static_cast<int>(fBinsMultCent.size()) - 1,
        static_cast<float*>(fBinsMultCent.data()),
        static_cast<int>(partArrayDATA.size()) - 1,
        static_cast<float*>(partArrayDATA.data()));
    fHistMCWeights->Sumw2();
    // printf("AliMCSpectraWeights: init histos successful\n"); // works
}

/**
 * @brief Load measured fractions (expert input) from alien
 */
void AliMCSpectraWeights::LoadMeasuredFractions() {
    // TFile *fMeasuredFile = AliDataFile::OpenOADB(fstFilePublished.Data());
    //  TFile *fMeasuredFile = TFile::Open(fstFilePublished.c_str(), "OPEN");
    auto fMeasuredFile = TFile::Open(fstFilePublished.c_str());
    if (!fMeasuredFile) {
        std::cerr << "AliMCSpectraWeights::Error: Could not load measured "
                     "fractions in "
                  << fstFilePublished << "\n";
        return;
    }

    for (auto& part : fstPartTypes) {
        if (part.find("Rest") != std::string::npos ||
            part.find("rest") != std::string::npos)
            continue; // there is no rest particles in measurement
        int _iPart = GetPartTypeNumber(part);
        for (auto& cent : fstCentralities) {
            // CollisionSystem:ParticleType:CentNumber:Stat/Sys:Function:FunctionVar
            std::string stHistName{fstCollisionSystem};
            stHistName += part;
            stHistName += cent;
            stHistName += "Stat";
            stHistName += AliMCSpectraWeights::GetFunctionFromSysFlag(fFlag);
            stHistName += AliMCSpectraWeights::GetSysVarFromSysFlag(fFlag);
            TH1D* hist = (TH1D*)fMeasuredFile->Get(stHistName.c_str());
            if (!hist) {
                std::cerr << "AliMCSpectraWeights::Error: could not find "
                          << stHistName << "\n";
                continue;
            }

            // hist-> pt:multcent:PartType
            std::array<float, 3> binEntry{0.};
            binEntry[2] = static_cast<float>(_iPart); // particle type
            binEntry[1] = AliMCSpectraWeights::GetMultFromCent(cent);
            if (!fUseMultiplicity)
                binEntry[1] = static_cast<float>(GetCentFromString(cent));

            for (int ipt = 0; ipt < fHistDataFractions->GetNbinsX(); ++ipt) {
                binEntry[0] = fHistDataFractions->GetXaxis()->GetBinCenter(ipt);
                if (binEntry[0] <
                    0.1) // pT cut; measurements start at 0.15 at best
                    continue;

                auto const _FractionValue =
                    hist->GetBinContent(hist->FindBin(binEntry[0]));
                FillTH3WithValue(fHistDataFractions, binEntry, _FractionValue);
            }
        }
    }
    fMeasuredFile->Close();
}

/**
 *  @brief Load the information from a previous train output
 *  @param[in] obj pointer to object from previous train output
 */
bool AliMCSpectraWeights::LoadFromAliMCSpectraWeight(AliMCSpectraWeights* obj) {
    // printf("AliMCSpectraWeights::DEBUG: Loading MC histogram from input
    // object\n");
    if (!obj)
        return false;
    fHistMCGenPrimTrackParticle = (TH3F*)obj->GetHistMCGenPrimTrackParticles();
    if (!fHistMCGenPrimTrackParticle) {
        std::cerr
            << "AliMCSpectraWeights::ERROR: problem with loading from object\n";
        return false;
    }
    fHistMCGenPrimTrackParticle->SetDirectory(0);
    return true;
}

/**
 *  @brief calculate the relative fractions of all particle species in MC
 *
 *  This require the input from a previous train output and is only calculated
 *  if the internal histogram of MC information is filled.
 */
bool AliMCSpectraWeights::CalcMCFractions() {
    if (!fHistMCGenPrimTrackParticle)
        return false;
    std::array<std::array<TH1D*, 10>, 20> _histMCFractions{
    nullptr};//FIXME: WARNING HARD CODED RANGES
    std::array<TH1D*, 20> _h1pTMCAll{nullptr};
    for (int icent = 0; icent < fNCentralities; ++icent) {
        auto const multTuple = AliMCSpectraWeights::GetMultTupleFromCent(icent);
        auto const _multFront = multTuple.front();
        auto const _multBack = multTuple.back();

        auto const _multBin1 =
            fHistMCGenPrimTrackParticle->GetYaxis()->FindBin(_multFront);
        auto const _multBbin2 =
            fHistMCGenPrimTrackParticle->GetYaxis()->FindBin(_multBack);

        for (int ipart = 0; ipart < fNPartTypes;
             ++ipart) { // calculate pt spectra of pi+k+p+Sigma+Rest
            int const _iPart =
                AliMCSpectraWeights::GetPartTypeNumber(fstPartTypes[ipart]);
            int const _iPartBin =
                fHistMCGenPrimTrackParticle->GetZaxis()->FindBin(_iPart);
            
            _histMCFractions[icent][ipart] =
                static_cast<TH1D*>(fHistMCGenPrimTrackParticle->ProjectionX(
                    Form("h1MCFraction_%s_%d", fstPartTypes[ipart].c_str(),
                         icent),
                    _multBin1, _multBbin2, _iPartBin, _iPartBin, "e"));
            if (!_histMCFractions[icent][ipart]) {
                std::cerr << "AliMCSpectraWeights::ERROR could not create "
                             "h1MCFraction\n";
                return false;
            }
            _histMCFractions[icent][ipart]->Scale(1, "width");
            if (!_h1pTMCAll[icent]) {
                _h1pTMCAll[icent] = (TH1D*)_histMCFractions[icent][ipart]->Clone(Form("h1pTMCAll_%d", icent));
            } else {
                _h1pTMCAll[icent]->Add(_histMCFractions[icent][ipart]);
            }
        }
        // all hist calculated
        // ------------------
        // calculate fractions
        for (auto _hist : _histMCFractions[icent]) {
            if (nullptr == _hist) {
                std::cerr << "AliMCSpectraWeights::ERROR could not calculate "
                             "fraction\n";
                continue;
            }

            int const ipart =
                AliMCSpectraWeights::GetPartTypeNumber(_hist->GetName());
            TH1D* h1MCFraction =
                (TH1D*)_hist->Clone(Form("%s_fraction", _hist->GetName()));
            h1MCFraction->Divide(_h1pTMCAll[icent]);
            // Set content of fractions to fHistMCFractions
            for (int ipt = 0; ipt < fHistMCFractions->GetNbinsX(); ++ipt) {
                float const pt =
                    fHistMCFractions->GetXaxis()->GetBinCenter(ipt);
                if (pt < 0)
                    continue;
                // fHistMCFractions : pt-mult-ipart
                std::array<float, 3> binEntry{
                    pt, AliMCSpectraWeights::GetMultFromCent(icent),
                    static_cast<float>(AliMCSpectraWeights::GetPartTypeNumber(
                        fstPartTypes[ipart]))};
                auto const _FractionValue =
                    h1MCFraction->GetBinContent(h1MCFraction->FindBin(pt));
                FillTH3WithValue(fHistMCFractions, binEntry, _FractionValue);
            }
            delete h1MCFraction;
        }
    }
    for (auto& histArray : _histMCFractions) {
        for (auto& hist : histArray) {
            if (hist){
                delete hist;
            hist = nullptr;}
        }
    }
    for (auto& hist : _h1pTMCAll) {
        if(hist) {delete hist;
            hist = nullptr;}
    }
    return true;
}

/**
 *  @brief correct the expert input for remaining particle species
 *
 *  The expert input only included pions, kaons, protons, sigma+ and sigma-
 * (constructed using lambda measurements). All other charged particles which
 * are included in the MC information are not measured, yet (xi, electron, muon,
 * omega-, ...). Therefore the fractions from measurement need to be corrected
 * for not having those particles included. (see analysis note)
 */
bool AliMCSpectraWeights::CorrectFractionsforRest() {
    if (!fHistMCGenPrimTrackParticle || !fHistDataFractions)
        return false;

    for (int icent = 0; icent < fNCentralities; ++icent) {
        auto multTuple = AliMCSpectraWeights::GetMultTupleFromCent(icent);
        auto _multFront = multTuple.front();
        auto _multBack = multTuple.back();

        auto _bin1 =
            fHistMCGenPrimTrackParticle->GetYaxis()->FindBin(_multFront);
        auto _bin2 =
            fHistMCGenPrimTrackParticle->GetYaxis()->FindBin(_multBack);

        auto h1pTMCAll = (TH1D*)fHistMCGenPrimTrackParticle->ProjectionX(
            "h1pTMCAll", _bin1, _bin2, 1,
            fHistMCGenPrimTrackParticle->GetNbinsZ(), "e");
        if (!h1pTMCAll) {
            std::cerr
                << "AliMCSpectraWeights::ERROR could not create h1pTMCAll\n";
            return false;
        }
        auto const _iRestPos = AliMCSpectraWeights::GetPartTypeNumber("Rest");
        auto const _RestBin = fHistMCGenPrimTrackParticle->GetZaxis()->FindBin(_iRestPos);
        auto h1RestCorrFactor = fHistMCGenPrimTrackParticle->ProjectionX(
            Form("h1RestCorrFactor_%d", icent), _bin1, _bin2, 1, _RestBin - 1,
            "e");
        h1RestCorrFactor->Divide(h1pTMCAll);
        for (int ipart = 0; ipart < fNPartTypes; ++ipart) {
            if ("Rest" == fstPartTypes[ipart])
                continue;
            for (int ipt = 0; ipt < fHistDataFractions->GetNbinsX(); ++ipt) {
                float pt = fHistDataFractions->GetXaxis()->GetBinCenter(ipt);
                // fHistMCFractions : pt-mult-ipart
                std::array<float, 3> binEntry{
                    pt,
                    static_cast<float>(
                        AliMCSpectraWeights::GetMultFromCent(icent)),
                    static_cast<float>(AliMCSpectraWeights::GetPartTypeNumber(
                    fstPartTypes[ipart]))};

                int const _iBinFind = fHistDataFractions->FindBin(
                    binEntry[0], binEntry[1], binEntry[2]);

                float value = fHistDataFractions->GetBinContent(_iBinFind);
                auto _RestCorrVal = h1RestCorrFactor->GetBinContent(
                    h1RestCorrFactor->FindBin(pt));
                if (0 == _RestCorrVal)
                    _RestCorrVal = 1;
                
                fHistDataFractions->SetBinContent(_iBinFind,
                                                  value * _RestCorrVal);
            }
        }
        delete h1pTMCAll;
        delete h1RestCorrFactor;
    }
    return true;
}

/**
 *  @brief calculate weight factors using internal information
 *
 *  Input from a previous train output is needed.
 *  Calculate for each particle species and multiplicity interval
 *  the ratio of fractionDATA to fractionMC. This is then used as
 *  weight factor later on.
 */
bool AliMCSpectraWeights::CalculateMCWeights() {
    if (!AliMCSpectraWeights::CalcMCFractions())
        return false;
    if (!AliMCSpectraWeights::CorrectFractionsforRest())
        return false;
    // correction of rest particles not measured in data fractions (see
    // AnalysisNote)

    for (int icent = 0; icent < fNCentralities; icent++) {
        for (int ipart = 0; ipart < fNPartTypes; ipart++) {
            //            if (ipart == GetPartTypeNumber("Rest"))
            //                continue;
            for (int ipt = 0; ipt < static_cast<int>(fBinsPt.size()); ++ipt) {
                float pt = fHistMCWeights->GetXaxis()->GetBinCenter(ipt);
                if (pt < 0)
                    continue;
                std::array<float, 3> binEntry{
                    pt, static_cast<float>(GetMultFromCent(icent)),
                    static_cast<float>(AliMCSpectraWeights::GetPartTypeNumber(
                        fstPartTypes[ipart]))};
                // Double_t binEntryData[3] = {pt, static_cast<Double_t>(icent),
                // static_cast<Double_t>(ipart)};

                auto const _iBinMC = fHistMCFractions->FindBin(
                    binEntry[0], binEntry[1], binEntry[2]);
                auto const _iBinData = fHistDataFractions->FindBin(
                    binEntry[0], binEntry[1], binEntry[2]);

                float const dFractionMC =
                    fHistMCFractions->GetBinContent(_iBinMC);

                float const dFractionData =
                    fHistDataFractions->GetBinContent(_iBinData);

                auto const _iBinWeight = fHistMCWeights->FindBin(
                    binEntry[0], binEntry[1], binEntry[2]);

                if (ipart == GetPartTypeNumber("Rest")) {
                    fHistMCWeights->SetBinContent(_iBinWeight, 1);
                } else if (dFractionMC != 0)
                    fHistMCWeights->SetBinContent(_iBinWeight,
                                                  dFractionData / dFractionMC);
                else
                    fHistMCWeights->SetBinContent(_iBinWeight, 1);

                fHistMCWeights->SetBinError(_iBinWeight, 1e-30);
            }
        }
    }
    return true;
}

/**
 *  @brief count the number of charged particles in the current event
 *
 *  Parse through the internal set AliMCEvent and count the number of particles
 *  according to the cuts and requirements of the analysis of the published
 *  spectra of identified charged particles used in the expert input.
 */
void AliMCSpectraWeights::CountEventMult() {
    fMultOrCent = 0;
    float lowPtCut = 0.05;
    float eta = 0.5;
    //    if (fstCollisionSystem.find("pp") != std::string::npos)
    //        eta = 0.5;
    AliStack* fMCStack = fMCEvent->Stack();
    for (int ipart = 0; ipart < fMCStack->GetNtrack(); ipart++) {
        TParticle* mcGenParticle = fMCStack->Particle(ipart);
        if (!fMCStack->IsPhysicalPrimary(ipart))
            continue; // secondary rejection
        if (TMath::Abs(mcGenParticle->GetPDG()->Charge()) < 0.01)
            continue; // neutral rejection
        float pEta = mcGenParticle->Eta();
        if (TMath::Abs(pEta) > eta)
            continue; // acceptance cut
        if (mcGenParticle->Pt() < lowPtCut)
            continue; // TODO: hard coded low pT cut
        ++fMultOrCent;
    }
}

/**
 *  @brief
 *  @param[in] mcGenParticle
 *  @param[in] eventMultiplicityOrCentrality
 *  @return
 */
float AliMCSpectraWeights::GetMCSpectraWeight(
    TParticle* mcGenParticle, float eventMultiplicityOrCentrality) {
    float weight = 1;
    if (!mcGenParticle->GetPDG()) {
        return 1;
    }
    if (TMath::Abs(mcGenParticle->GetPDG()->Charge()) < 0.01) {
        return 1;
    }
    int particleType = AliMCSpectraWeights::IdentifyMCParticle(mcGenParticle);
    if (particleType == GetPartTypeNumber("Rest")) {
        return 1;
    }
    if (fbTaskStatus == AliMCSpectraWeights::TaskState::kMCWeightCalculated) {
        // rest particles can not be tuned
        float icent =
            AliMCSpectraWeights::GetCentFromMult(eventMultiplicityOrCentrality);
        float pt = mcGenParticle->Pt();
        if (pt < 0.15)
            return 1;
        if (pt >= 20)
            pt = 19.9;
        std::array<float, 3> binEntry{
            pt, static_cast<float>(AliMCSpectraWeights::GetMultFromCent(icent)),
            static_cast<float>(particleType)};

        auto const _iBin =
            fHistMCWeights->FindBin(binEntry[0], binEntry[1], binEntry[2]);

        weight = fHistMCWeights->GetBinContent(_iBin);
        if (weight <= 0)
            weight = 1;
    }
    return weight;
}

/**
 *  @brief
 *  @param[in] mcGenParticle
 *  @param[in] mcEvent
 *  @return
 */
float AliMCSpectraWeights::GetMCSpectraWeight(TParticle* mcGenParticle,
                                              AliMCEvent* mcEvent) {
    if (mcEvent != fMCEvent) {
        fMCEvent = mcEvent;
        AliMCSpectraWeights::CountEventMult();
    }
    return AliMCSpectraWeights::GetMCSpectraWeight(mcGenParticle, fMultOrCent);
}

/**
 *  @brief
 *  @param[in] mcEvent
 *
 */
void AliMCSpectraWeights::FillMCSpectra(AliMCEvent* mcEvent) {
    if (fbTaskStatus >= AliMCSpectraWeights::TaskState::kMCSpectraObtained)
        return;
    if (mcEvent != fMCEvent) {
        fMCEvent = mcEvent;
        AliMCSpectraWeights::CountEventMult();
    }

    AliStack* MCStack = fMCEvent->Stack();
    if (!MCStack) {
        printf("AliMCSpectraWeights::ERROR: fMCStack not available\n");
        return;
    }

    for (int iParticle = 0; iParticle < MCStack->GetNtrack(); ++iParticle) {
        TParticle* mcGenParticle = MCStack->Particle(iParticle);
        if (!mcGenParticle) {
            printf(
                "AliMCSpectraWeights::WARNING: mcGenParticle not available\n");
            continue;
        }
        if (!mcGenParticle->GetPDG())
            continue;
        if (!MCStack->IsPhysicalPrimary(iParticle))
            continue;
        if (TMath::Abs(mcGenParticle->GetPDG()->Charge()) < 0.01)
            continue;

        float partEta = mcGenParticle->Eta();
        float _maxEta = 0.5; // hard coded max eta; in all papers 0.5

        if (TMath::Abs(partEta) > _maxEta)
            continue; // apply same acceptance as in published spectra
        int particleType =
            AliMCSpectraWeights::IdentifyMCParticle(mcGenParticle);
        if (particleType < 0)
            continue;
        //        std::array<float, 3>
        //        binEntry{static_cast<float>(mcGenParticle->Pt()),
        //                                      fMultOrCent,
        //                                      static_cast<float>(particleType)};
        fHistMCGenPrimTrackParticle->Fill(
            static_cast<float>(mcGenParticle->Pt()), fMultOrCent,
            static_cast<float>(particleType));
    }
}

/**
 *  @brief
 *  @param[in] mcParticle
 *  @return
 */
int AliMCSpectraWeights::IdentifyMCParticle(TParticle* mcParticle) {
    int ipdg = TMath::Abs(
        mcParticle
            ->GetPdgCode()); // Abs() because antiparticles are negaitve...
    if (ipdg == 211)
        return GetPartTypeNumber(AliMCSpectraWeights::ParticleType::kPion);
    if (ipdg == 321)
        return GetPartTypeNumber(AliMCSpectraWeights::ParticleType::kKaon);
    if (ipdg == 2212)
        return GetPartTypeNumber(AliMCSpectraWeights::ParticleType::kProtons);
    if (ipdg == 3222)
        return GetPartTypeNumber(AliMCSpectraWeights::ParticleType::kSigmaPlus);
    if (ipdg == 3112)
        return GetPartTypeNumber(
            AliMCSpectraWeights::ParticleType::kSigmaMinus);
    // if(ipdg==3334) return AliMCSpectraWeights::ParticleType::kOmegaMinus;
    // if(ipdg==3312) return AliMCSpectraWeights::ParticleType::kXiMinus;
    // if(ipdg==11) return AliMCSpectraWeights::ParticleType::kElectron;
    // if(ipdg==13) return AliMCSpectraWeights::ParticleType::kMuon;
    // printf("AliMCSpectraWeights:: pdf code of rest particle %d\n", ipdg);
    return GetPartTypeNumber(AliMCSpectraWeights::ParticleType::kRest);
}

/**
 *  @brief
 *  @param[in] CentBin
 *  @return
 */
double AliMCSpectraWeights::GetMultFromCent(int CentBin) const {
    if (fstCollisionSystem.find("pp") != std::string::npos &&
        fstCollisionSystem.find("ppb") == std::string::npos) {
        // for | eta | < 0.5
        switch (CentBin) {
        case 0:
            return 21.3;
        case 1:
            return 16.5;
        case 2:
            return 13.5;
        case 3:
            return 11.5;
        case 4:
            return 10.1;
        case 5:
            return 8.45;
        case 6:
            return 6.72;
        case 7:
            return 5.4;
        case 8:
            return 3.9;
        case 9:
            return 2.26;
        default:
            return -2.0;
        }
    } else if (fstCollisionSystem.find("ppb") !=
               std::string::npos)
    {
        // for | eta | < 0.5
        switch (CentBin) {
        case 0:
            return 45.0;
        case 1:
            return 36.2;
        case 2:
            return 30.5;
        case 3:
            return 23.2;
        case 4:
            return 16.1;
        case 5:
            return 9.8;
        case 6:
            return 4.4;
        default:
            return -2.0;
        }
    } else if (fstCollisionSystem.find("pbpb") != std::string::npos) {
        // for | eta | < 0.5
        //        Centrality 0–5% 5–10% 10–20% 20–30% 30–40% 40–50% 50–60%
        //        60–70% 70–80% dNch /dη 1601 ± 60 1294 ± 49 966±37 649±23
        //        426±15 261±9 149±6 76±4 35±2
        switch (CentBin) {
        case 0:
            return 1601;
        case 1:
            return 1294;
        case 2:
            return 966;
        case 3:
            return 649;
        case 4:
            return 426;
        case 5:
            return 261;
        case 6:
            return 149;
        case 7:
            return 76;
        case 8:
            return 35;
        default:
            return -2.0;
        }
    } else if (fstCollisionSystem.find("xexe") != std::string::npos) {
        // for | eta | < 0.5
        switch (CentBin) {
        case 0:
            return 1167;
        case 1:
            return 939;
        case 2:
            return 706;
        case 3:
            return 478;
        case 4:
            return 315;
        case 5:
            return 198;
        case 6:
            return 118;
        case 7:
            return 64.7;
        case 8:
            return 32;
        default:
            return -2.0;
        }
    }

    return -1;
}

/**
 *  @brief
 *  @param[in] cent
 *  @return
 */
double AliMCSpectraWeights::GetMultFromCent(std::string cent) {
    return GetMultFromCent(GetCentFromString(cent));
}

/**
 *  @brief
 *  @param[in] CentBin
 *  @return
 */
std::vector<float>
AliMCSpectraWeights::GetMultTupleFromCent(int CentBin) const {
    float const currentMult = AliMCSpectraWeights::GetMultFromCent(CentBin);
    float nextMult = 0;
    float previousMult = 0;
    float dMultHigh = 0;
    float dMultLow = 0;
    if (CentBin < fNCentralities - 1)
        nextMult = AliMCSpectraWeights::GetMultFromCent(CentBin + 1);
    if (CentBin > 0)
        previousMult = AliMCSpectraWeights::GetMultFromCent(CentBin - 1);

    if (CentBin == 0) {
        if (fstCollisionSystem.find("pp") != std::string::npos && fstCollisionSystem.find("ppb") == std::string::npos)
            dMultHigh = 49;
        else if (fstCollisionSystem.find("ppb") != std::string::npos)
            dMultHigh = 299;
        else if (fstCollisionSystem.find("pbpb") != std::string::npos)
            dMultHigh = 3499;
    } else
        dMultHigh = (previousMult + currentMult) / 2.;
    if (CentBin == fNCentralities - 1)
        dMultLow = 0.;
    else
        dMultLow = (currentMult + nextMult) / 2.;

    return {dMultLow, dMultHigh};
}

/**
 *  @brief
 *  @param[in] cent
 *  @return
 */
// TODO: implement
int AliMCSpectraWeights::GetCentFromString(std::string cent) {
    for (int i = 0; i < static_cast<int>(fstCentralities.size()); ++i) {
        if (fstCentralities[i].find(cent) != std::string::npos)
            return i;
    }
    return -1;
}

/**
 *  @brief
 *  @param[in] dMult
 *  @return
 */
// TODO: implement
double AliMCSpectraWeights::GetCentFromMult(double dMult) {
    if (fstCollisionSystem.find("pp") != std::string::npos) {
        if (dMult > 18)
            return 0;
        if (dMult > 14.5)
            return 1;
        if (dMult > 12)
            return 2;
        if (dMult > 10.9)
            return 3;
        if (dMult > 9)
            return 4;
        if (dMult > 7)
            return 5;
        if (dMult > 5.7)
            return 6;
        if (dMult > 4.3)
            return 7;
        if (dMult > 3)
            return 8;
        if (dMult > 0)
            return 9;
    } else if (fstCollisionSystem.find("ppb") !=
               std::string::npos) // TODO: include other systems
    {
        //        switch (dMult) {
        //            case 0:
        //                return 45.0;
        //            case 1:
        //                return 36.2;
        //            case 2:
        //                return 30.5;
        //            case 3:
        //                return 23.2;
        //            case 4:
        //                return 16.1;
        //            case 5:
        //                return 9.8;
        //            case 6:
        //                return 4.4;
        //            default:
        //                return -2.0;
        //        }

    } else if (fstCollisionSystem.find("pbpb") != std::string::npos) {

    } else if (fstCollisionSystem.find("xexe") != std::string::npos) {
    }

    return -1;
}

/**
 *  @brief
 *  @param[in] type
 *  @return
 */
int AliMCSpectraWeights::GetPartTypeNumber(
    AliMCSpectraWeights::ParticleType type) {
    switch (type) {
    case AliMCSpectraWeights::ParticleType::kPion:
        return AliMCSpectraWeights::ParticleType::kPion;
    case AliMCSpectraWeights::ParticleType::kProtons:
        return AliMCSpectraWeights::ParticleType::kProtons;
    case AliMCSpectraWeights::ParticleType::kKaon:
        return AliMCSpectraWeights::ParticleType::kKaon;
    case AliMCSpectraWeights::ParticleType::kSigmaMinus:
        return AliMCSpectraWeights::ParticleType::kSigmaMinus;
    case AliMCSpectraWeights::ParticleType::kSigmaPlus:
        return AliMCSpectraWeights::ParticleType::kSigmaPlus;
    case AliMCSpectraWeights::ParticleType::kRest:
        return AliMCSpectraWeights::ParticleType::kRest;
    default:
        return -1;
        break;
    }
}

/**
 *  @brief
 *  @param[in] Particle
 *  @return
 */
int AliMCSpectraWeights::GetPartTypeNumber(std::string Particle) {
    if (Particle.find("Pion") != std::string::npos) {
        return AliMCSpectraWeights::GetPartTypeNumber(
            AliMCSpectraWeights::ParticleType::kPion);
    } else if (Particle.find("Proton") != std::string::npos) {
        return AliMCSpectraWeights::GetPartTypeNumber(
            AliMCSpectraWeights::ParticleType::kProtons);
    } else if (Particle.find("Kaon") != std::string::npos) {
        return AliMCSpectraWeights::GetPartTypeNumber(
            AliMCSpectraWeights::ParticleType::kKaon);
    } else if (Particle.find("SigmaMinus") != std::string::npos) {
        return AliMCSpectraWeights::GetPartTypeNumber(
            AliMCSpectraWeights::ParticleType::kSigmaMinus);
    } else if (Particle.find("SigmaPlus") != std::string::npos) {
        return AliMCSpectraWeights::GetPartTypeNumber(
            AliMCSpectraWeights::ParticleType::kSigmaPlus);
    } else if (Particle.find("Rest") != std::string::npos) {
        return AliMCSpectraWeights::GetPartTypeNumber(
            AliMCSpectraWeights::ParticleType::kRest);
    } else
        return -1;
}

// TODO: implement flags
/**
 *  @brief
 *  @param[in] flag
 *  @return
 */
std::string AliMCSpectraWeights::GetFunctionFromSysFlag(SysFlag flag) {
    switch (flag) {
    case SysFlag::kNominal:
        return "Bylinkin";
    case SysFlag::kPionUp:
        return "";
    case SysFlag::kPionDown:
        return "";
    case SysFlag::kProtonUp:
        return "";
    case SysFlag::kProtonDown:
        return "";
    case SysFlag::kKaonUp:
        return "";
    case SysFlag::kKaonDown:
        return "";
    case SysFlag::kSigmaPlusUp:
        return "";
    case SysFlag::kSigmaPlusDown:
        return "";
    case SysFlag::kSigmaMinusUp:
        return "";
    case SysFlag::kSigmaMinusDown:
        return "";
    case SysFlag::kBylinkinLower:
        return "";
    case SysFlag::kBylinkinUpper:
        return "";
    case SysFlag::kHagedorn:
        return "";
    case SysFlag::kHagedornUpper:
        return "";
    case SysFlag::kHagedornLower:
        return "";
    case SysFlag::kExponential:
        return "";
    case SysFlag::kExponentialUpper:
        return "";
    case SysFlag::kExponentialLower:
        return "";
    case SysFlag::kBlastwave:
        return "";
    case SysFlag::kBlastwaveUpper:
        return "";
    case SysFlag::kBlastwaveLower:
        return "";
    default:
        return "Bylinkin";
    }

    return "Bylinkin";
}
// TODO: implement
std::string AliMCSpectraWeights::GetSysVarFromSysFlag(SysFlag flag) {
    switch (flag) {
    case SysFlag::kNominal:
        return "";
    case SysFlag::kPionUp:
        return "";
    case SysFlag::kPionDown:
        return "";
    case SysFlag::kProtonUp:
        return "";
    case SysFlag::kProtonDown:
        return "";
    case SysFlag::kKaonUp:
        return "";
    case SysFlag::kKaonDown:
        return "";
    case SysFlag::kSigmaPlusUp:
        return "";
    case SysFlag::kSigmaPlusDown:
        return "";
    case SysFlag::kSigmaMinusUp:
        return "";
    case SysFlag::kSigmaMinusDown:
        return "";
    case SysFlag::kBylinkinLower:
        return "";
    case SysFlag::kBylinkinUpper:
        return "";
    case SysFlag::kHagedorn:
        return "";
    case SysFlag::kHagedornUpper:
        return "";
    case SysFlag::kHagedornLower:
        return "";
    case SysFlag::kExponential:
        return "";
    case SysFlag::kExponentialUpper:
        return "";
    case SysFlag::kExponentialLower:
        return "";
    case SysFlag::kBlastwave:
        return "";
    case SysFlag::kBlastwaveUpper:
        return "";
    case SysFlag::kBlastwaveLower:
        return "";
    default:
        return "";
    }

    return "";
}
