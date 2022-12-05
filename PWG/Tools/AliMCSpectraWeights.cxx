/**
 * @file AliMCSpectraWeights.cxx
 * @brief Class for re-weighting particle abundances in MC simulations
 * @author Patrick Huhn
 * @date 25/10/2019
 */
#include "AliAnalysisUtils.h"
#include "AliMCEvent.h"
#include "AliMCSpectraWeights.h"
#include "AliMCParticle.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH3F.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TDirectory.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <math.h>
#include <memory>
#include <string>
#include "AliLog.h"

//AliDebug(3, std::to_string(x).c_str())
#ifdef __AliMCSpectraWeights_DebugPCC__
#define DebugPCC(x) std::cerr << x
#else
#define DebugPCC(x)
#endif

#ifdef __AliMCSpectraWeights_DebugTiming__
#define DebugChrono(x) std::cerr << x
#else
#define DebugChrono(x)
#endif

int const GetBinFromTH3(TH3F* h, std::array<float, 3> const& _values) {
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
: TNamed(), fstCollisionSystem("pp"), fstFileMCSpectra(""),
fstFilePublished(""), fstSavedObjName("fMCSpectraWeights"),
fstSavedListName("AliMCWeightsTask"), fstPartTypes(0),
fstCentralities(0), fBinsMultCent{0}, fBinsPt{0},
fHistMCGenPrimTrackParticle(nullptr), fHistDataFractions(nullptr),
fHistMCFractions(nullptr), fHistMCWeights(nullptr), fHistMCWeightsSysUp(nullptr), fHistMCWeightsSysDown(nullptr), fMCEvent(nullptr),
fMultOrCent(0), fNPartTypes(6), fNCentralities(0),
fbTaskStatus(AliMCSpectraWeights::TaskState::kAllEmpty),
fFlag(AliMCSpectraWeights::SysFlag::kNominal), fUseMultiplicity(kTRUE),
fUseMBFractions(kFALSE), fDoInterpolation(kTRUE) {}

/**
 *  @brief standard way for constuctor
 *  @param[in] collisionSystem string for selecting the used collision system.
 * Supported systems are: pp, pPb, PbPb and XeXe
 *  @param[in] stName string for naming the object. Not used anymore. Included
 * only for backward compability.
 *  @param[in] flag used for systematic variations.
 */
AliMCSpectraWeights::AliMCSpectraWeights(std::string const& collisionSystem,
                                         std::string const& stName,
                                         AliMCSpectraWeights::SysFlag flag)
: TNamed(stName.c_str(), stName.c_str()), fstCollisionSystem("pp"),
fstFileMCSpectra(""), fstFilePublished(""),
fstSavedObjName("fMCSpectraWeights"), fstSavedListName("AliMCWeightsTask"),
fstPartTypes(), fstCentralities(), fBinsMultCent{0}, fBinsPt{0},
fHistMCGenPrimTrackParticle(nullptr), fHistDataFractions(nullptr),
fHistMCFractions(nullptr), fHistMCWeights(nullptr),  fHistMCWeightsSysUp(nullptr), fHistMCWeightsSysDown(nullptr), fMCEvent(nullptr),
fMultOrCent(0), fNPartTypes(6), fNCentralities(0),
fbTaskStatus(AliMCSpectraWeights::TaskState::kAllEmpty), fFlag(flag),
fUseMultiplicity(kTRUE), fUseMBFractions(kFALSE), fDoInterpolation(kTRUE) {
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    fstCollisionSystem = collisionSystem; // copy here; ok
    std::for_each(
                  fstCollisionSystem.begin(), fstCollisionSystem.end(),
                  [](char& c) { c = std::tolower(static_cast<unsigned char>(c)); });

    DebugPCC("AliMCSpectraWeights with DebugPCC info for " +fstCollisionSystem
             + " collisions\n");

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

    // init random generator
    frndGen = TRandom3(0);

    // set default Binning
    // pT binning
    fBinsPt = {0.0,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,   0.45, 0.5,
        0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85,  0.9,  0.95,
        1.0,  1.1,  1.2,  1.4,  1.6,  1.8,  2.0,   2.2,  2.4,
        2.6,  2.8,  3.0,  3.2,  3.6,  4.0,  5.0,   6.0,  8.0,
        10.0, 13.0, 20.0, 30.0, 50.0, 80.0, 100.0, 200.0};
    // multiplicity binning
    if (fstCollisionSystem == "pp") {
        fBinsPt = {0.0, 0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.5,5,5.5,6,6.5,7,8,10,13,20, 30, 50, 80, 100, 200};
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
    } else if (fstCollisionSystem == "ppb") {
        fBinsMultCent.reserve(301);
        for (int i = 0; i < 301; ++i) {
            fBinsMultCent.push_back(i);
        }

        std::vector<std::string> cent{"0",    "1", "2", "3",
            "4", "5", "6"};
        fstCentralities.clear();

        for (auto const& _c : cent) {
            fstCentralities.push_back(_c);
        }
        fNCentralities = 7;
//        DebugPCC("size of fstCentralities " << fNCentralities << "\n");
    } else if (fstCollisionSystem == "pbpb") {
        fBinsMultCent.reserve(201);
        for (int i = 0; i < 201; ++i) {
            fBinsMultCent.push_back(i * 25);
        }
        std::vector<std::string> cent{"MB",    "c0005", "c0510", "c0010",
                                      "c1020", "c2040", "c4060", "c6080"};
        fNCentralities = 8;
        fstCentralities.clear();
        fstCentralities = cent;
    } else if (fstCollisionSystem == "xexe") {
        fBinsMultCent.reserve(201);
        for (int i = 0; i < 201; ++i) {
            fBinsMultCent.push_back(i * 25);
        }
        std::vector<std::string> cent{"MB",    "c0005", "c0510", "c0010",
            "c1020", "c2040", "c4060", "c6080"};
        fNCentralities = 8;
        fstCentralities.clear();
        fstCentralities = cent;
    } else {
        fBinsMultCent.reserve(101);
        for (int i = 0; i < 101; ++i) {
            fBinsMultCent.push_back(i);
        }
        fNCentralities = 10;
    }
    //    if (fstCollisionSystem.find("pp") == std::string::npos) {
    //        std::vector<std::string> cent{"MB",    "c0005", "c0510", "c0010",
    //                                      "c1020", "c2040", "c4060", "c6080"};
    //        fNCentralities = 8;
    //        fstCentralities.clear();
    //        fstCentralities = cent;
    //    }
    fstPartTypes = {"Pion",       "Proton",    "Kaon",
        "SigmaMinus", "SigmaPlus", "Rest", "Lambda"};

    // setup systematics
    fAllSystematicFlags = {
        AliMCSpectraWeights::SysFlag::kNominal,
        AliMCSpectraWeights::SysFlag::kPionUp,
        //                           AliMCSpectraWeights::SysFlag::kPionDown,
        AliMCSpectraWeights::SysFlag::kProtonUp,
        //                           AliMCSpectraWeights::SysFlag::kProtonDown,
        AliMCSpectraWeights::SysFlag::kKaonUp,
        //                           AliMCSpectraWeights::SysFlag::kKaonDown,
        AliMCSpectraWeights::SysFlag::kSigmaPlusUp,
        //                           AliMCSpectraWeights::SysFlag::kSigmaPlusDown,
        AliMCSpectraWeights::SysFlag::kSigmaMinusUp,
        //                           AliMCSpectraWeights::SysFlag::kSigmaMinusDown

    };

    if ("pp" == fstCollisionSystem) {
        fAllSystematicFlags.push_back(
                                      AliMCSpectraWeights::SysFlag::kBylinkinUpper);
        fAllSystematicFlags.push_back(
                                      AliMCSpectraWeights::SysFlag::kBylinkinLower);
        fAllSystematicFlags.push_back(AliMCSpectraWeights::SysFlag::kHagedorn);
//        fAllSystematicFlags.push_back(
//                                      AliMCSpectraWeights::SysFlag::kHagedornUpper);
//        fAllSystematicFlags.push_back(
//                                      AliMCSpectraWeights::SysFlag::kHagedornLower);
    } else if ("ppb" == fstCollisionSystem) {
        fAllSystematicFlags.push_back(
                                      AliMCSpectraWeights::SysFlag::kBylinkinUpper);
        fAllSystematicFlags.push_back(
                                      AliMCSpectraWeights::SysFlag::kBylinkinLower);
        fAllSystematicFlags.push_back(AliMCSpectraWeights::SysFlag::kHagedorn);
//        fAllSystematicFlags.push_back(
//                                      AliMCSpectraWeights::SysFlag::kHagedornUpper);
//        fAllSystematicFlags.push_back(
//                                      AliMCSpectraWeights::SysFlag::kHagedornLower);
    } else if ("pbpb" == fstCollisionSystem || "xexe" == fstCollisionSystem) {
//        fAllSystematicFlags.push_back(
//                                      AliMCSpectraWeights::SysFlag::kBlastwaveUpper);
//        fAllSystematicFlags.push_back(
//                                      AliMCSpectraWeights::SysFlag::kBlastwaveLower);
        fAllSystematicFlags.push_back(AliMCSpectraWeights::SysFlag::kHagedorn);
//        fAllSystematicFlags.push_back(
//                                      AliMCSpectraWeights::SysFlag::kHagedornUpper);
//        fAllSystematicFlags.push_back(
//                                      AliMCSpectraWeights::SysFlag::kHagedornLower);
    }

    fNSysFlags = fAllSystematicFlags.size();
    fbTaskStatus = AliMCSpectraWeights::TaskState::kAllEmpty;
    fstFilePublished =
//    "/home/alidock/ExpertInput/ExpertInput.root";
    "alien:///alice/cern.ch/user/p/phuhn/AllPublishedFractions.root";
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    DebugChrono("Construction took " << duration << " microseconds\n");
#endif
}

AliMCSpectraWeights::~AliMCSpectraWeights() {
    //    if(fHistMCGenPrimTrackParticle) delete fHistMCGenPrimTrackParticle; //
    //    is put into train output in AnalysisTask if(fHistDataFractions) delete
    //    fHistDataFractions; if(fHistMCFractions) delete fHistMCFractions;
    //    if(fHistMCWeights) delete fHistMCWeights;
    if (fMCEvent)
        fMCEvent = 0;
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
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
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
                if (!fInput->Get("AliMCWeightsTask"))
                    std::cerr << "AliMCSpectraWeights::WARNING: more than 1 "
                    "list in the streamed file; please specify; "
                    "using 1st list;\n\n";
            } else
                fstSavedListName = fInput->GetListOfKeys()->At(0)->GetName();
            // printf("AliMCSpectraWeights:: Loading %s from list %s\n",
            // fstSavedObjName.Data(), fstSavedListName.Data());
            auto dir = (TDirectory*)fInput->Get(fstSavedListName.c_str());
            TList* listMC = (TList*)dir->Get(fstSavedListName.c_str());
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
                } else {
                    fHistMCGenPrimTrackParticle = (TH3F*)tmpHist->Clone(
                                                                        "fHistMCGenPrimTrackParticle_prev");
                    fHistMCGenPrimTrackParticle->SetDirectory(0);
                    if (fHistMCGenPrimTrackParticle->GetEntries() > 0) {
                        DebugPCC("Previous train has mc tracks\n");
                        fbTaskStatus =
                        AliMCSpectraWeights::TaskState::kMCSpectraObtained;
                        AliMCSpectraWeights::CalcMCFractions();
                    }
                }
            }
            fInput->Close();
            delete fInput;
        }
    }

    if (fDoSystematics) {
        if (fbTaskStatus ==
            AliMCSpectraWeights::TaskState::kMCSpectraObtained) {
            for (auto const& flag : fAllSystematicFlags) {
                fFlag = flag;
                fbTaskStatus =
                AliMCSpectraWeights::TaskState::kMCSpectraObtained;
                AliMCSpectraWeights::LoadMeasuredFractions();
                fbTaskStatus =
                AliMCSpectraWeights::TaskState::kDataFractionLoaded;
                if (AliMCSpectraWeights::CorrectFractionsforRest() &&
                    AliMCSpectraWeights::CalculateMCWeights()) {
                    fbTaskStatus =
                    AliMCSpectraWeights::TaskState::kMCWeightCalculated;
                }
            }
            // Set final outputs to nominal again;
            fFlag = AliMCSpectraWeights::SysFlag::kNominal;
//            AliMCSpectraWeights::LoadMeasuredFractions();
//            AliMCSpectraWeights::CorrectFractionsforRest();
//            AliMCSpectraWeights::CalculateMCWeights();
            AliMCSpectraWeights::CalculateSystematicUncertainties();
        }
    } else {
        // Loading measured fractions
        if (fbTaskStatus ==
            AliMCSpectraWeights::TaskState::kMCSpectraObtained) {
            AliMCSpectraWeights::LoadMeasuredFractions();
            fbTaskStatus = AliMCSpectraWeights::TaskState::kDataFractionLoaded;
        }
        // Calculating weight factors
        if (fbTaskStatus ==
            AliMCSpectraWeights::TaskState::kDataFractionLoaded) {
            if (AliMCSpectraWeights::CorrectFractionsforRest() &&
                AliMCSpectraWeights::CalculateMCWeights()) {
                fbTaskStatus =
                AliMCSpectraWeights::TaskState::kMCWeightCalculated;
            }
        }
    }

DebugPCC("AliMCSpectraWeights::INFO: Init finished with status " + std::to_string(fbTaskStatus)+"\n");
#ifdef __AliMCSpectraWeights_DebugTiming__
auto t2 = std::chrono::high_resolution_clock::now();
auto duration =
std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
DebugChrono("Init took " << duration << " microseconds\n");
#endif
}

/**
 *  @brief Create all internal histograms
 */
void AliMCSpectraWeights::InitHistos() {
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    DebugPCC("Initializing histograms\n");
    // Initalizing histograms
    // histogram charged patricles pt:multcent:type
    std::array<float, 8> partArray{};
    std::array<float, 6> partArrayDATA{};
    for (int i = 0; i < 8; ++i) {
        partArray[i] = -0.5 + i;
        if (i < 6)
            partArrayDATA[i] = -0.5 + i;
    }

    //    std::cout << "Create multiplicity binning\n";
    std::vector<float> _NchBinning{};
    if (!fUseMBFractions) {
        int iCent = 0;
        for (auto const& cent : fstCentralities) {
            auto const _Nch =
            AliMCSpectraWeights::GetMultTupleFromCent(iCent++);
            auto const _front = _Nch.front();
            auto const _back = _Nch.back();
            if (_NchBinning.size() < 1) {
                if (_front > _back) {
                    _NchBinning.push_back(_front);
                    _NchBinning.push_back(_back);
                } else {
                    _NchBinning.push_back(_back);
                    _NchBinning.push_back(_front);
                }
            } else {
                if (_front > _back) {
                    _NchBinning.push_back(_back);
                } else {
                    _NchBinning.push_back(_front);
                }
            }
        }
    } else {
        _NchBinning.push_back(10000.5);
        _NchBinning.push_back(0.5);
    }
    std::reverse(_NchBinning.begin(), _NchBinning.end());
    if (fUseMBFractions && "pp" == fstCollisionSystem) {
        fNCentralities = 1;
        fstCentralities = {"0"};
        fBinsMultCent = _NchBinning;
    }
#ifdef __AliMCSpectraWeights_DebugPCC__
    DebugPCC("_NchBinning:\n");
    auto const printInt = [](float const& n) { std::cout << "\t" << n; };
    auto const printString = [](std::string const& n) { std::cout << "\t" << n; };
    std::for_each(_NchBinning.begin(), _NchBinning.end(), printInt);
    std::cout << "\n .. done" << std::endl;
    DebugPCC("\n");
    DebugPCC("fstCentralities:\n");
    std::for_each(fstCentralities.begin(), fstCentralities.end(), printString);
    std::cout << "\n .. done" << std::endl;
#endif

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

    fHistDataFractions =
    new TH3F("fHistDataFractions",
             "DATA fractions histogram;#it{p}_{T} "
             "(GeV/#it{c});multiplicity or centrality;Particle type",
             static_cast<int>(fBinsPt.size()) - 1,
             static_cast<float*>(fBinsPt.data()),
             static_cast<int>(_NchBinning.size()) - 1,
             static_cast<float*>(_NchBinning.data()),
             static_cast<int>(partArrayDATA.size()) - 1,
             static_cast<float*>(partArrayDATA.data()));

    fHistMCFractions =
    new TH3F("fHistMCFractions",
             "MC fractions histogram;#it{p}_{T} (GeV/#it{c});multiplicity "
             "or centrality;Particle type",
             static_cast<int>(fBinsPt.size()) - 1,
             static_cast<float*>(fBinsPt.data()),
             static_cast<int>(_NchBinning.size()) - 1,
             static_cast<float*>(_NchBinning.data()),
             static_cast<int>(partArray.size()) - 1,
             static_cast<float*>(partArray.data()));

    fHistMCWeights = new TH3F(
                              "fHistMCWeights",
                              "MC weight histogram for charged particle composition;#it{p}_{T} "
                              "(GeV/#it{c});multiplicity or centrality;Particle type",
                              static_cast<int>(fBinsPt.size()) - 1,
                              static_cast<float*>(fBinsPt.data()),
                              static_cast<int>(_NchBinning.size()) - 1,
                              static_cast<float*>(_NchBinning.data()),
                              static_cast<int>(partArrayDATA.size()) - 1,
                              static_cast<float*>(partArrayDATA.data()));
    DebugPCC("AliMCSpectraWeights: init histos successful\n");

    if (fDoSystematics) {
        DebugPCC("AliMCSpectraWeights: Init sys histos\n");
        for (auto& flag : fAllSystematicFlags) {
            std::string const _histName{"fHistMCWeights" +
                GetFunctionFromSysFlag(flag) +
                GetSysVarFromSysFlag(flag)};
            DebugPCC("\tHist name: "+ _histName + "\n");

            fHistMCWeightsSys[flag] = new TH3F(
                                               _histName.c_str(),
                                               "MC weight histogram for charged particle "
                                               "composition;#it{p}_{T} "
                                               "(GeV/#it{c});multiplicity or centrality;Particle type",
                                               static_cast<int>(fBinsPt.size()) - 1,
                                               static_cast<float*>(fBinsPt.data()),
                                               static_cast<int>(_NchBinning.size()) - 1,
                                               static_cast<float*>(_NchBinning.data()),
                                               static_cast<int>(partArrayDATA.size()) - 1,
                                               static_cast<float*>(partArrayDATA.data()));
        }

        fHistMCWeightsSysUp = new TH3F(
                                  "fHistMCWeightsSysUp",
                                  "MC weight histogram for charged particle composition;#it{p}_{T} "
                                  "(GeV/#it{c});multiplicity or centrality;Particle type",
                                  static_cast<int>(fBinsPt.size()) - 1,
                                  static_cast<float*>(fBinsPt.data()),
                                  static_cast<int>(_NchBinning.size()) - 1,
                                  static_cast<float*>(_NchBinning.data()),
                                  static_cast<int>(partArrayDATA.size()) - 1,
                                  static_cast<float*>(partArrayDATA.data()));

        fHistMCWeightsSysDown = new TH3F(
                                       "fHistMCWeightsSysDown",
                                       "MC weight histogram for charged particle composition;#it{p}_{T} "
                                       "(GeV/#it{c});multiplicity or centrality;Particle type",
                                       static_cast<int>(fBinsPt.size()) - 1,
                                       static_cast<float*>(fBinsPt.data()),
                                       static_cast<int>(_NchBinning.size()) - 1,
                                       static_cast<float*>(_NchBinning.data()),
                                       static_cast<int>(partArrayDATA.size()) - 1,
                                       static_cast<float*>(partArrayDATA.data()));


        DebugPCC("AliMCSpectraWeights: init histos systematics successful\n");
    }
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    DebugChrono("InitHistos took " << duration << " microseconds\n");
#endif
}

/**
 * @brief Load measured fractions (expert input) from alien
 */
void AliMCSpectraWeights::LoadMeasuredFractions() {
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    DebugPCC("Load measured fractions\n");
    auto fMeasuredFile = TFile::Open(fstFilePublished.c_str());
    if (!fMeasuredFile) {
        std::cerr << "AliMCSpectraWeights::Error: Could not load measured "
        "fractions in "
        << fstFilePublished << "\n";
        return;
    }
        fHistDataFractions->Reset(); // clean up for new input
    for (auto& part : fstPartTypes) {
//        DebugPCC("\tPart: " << part << "\n");
        if (part.find("Rest") != std::string::npos ||
            part.find("rest") != std::string::npos ||
            part.find("Lambda") != std::string::npos ||
            part.find("lambda") != std::string::npos){
            DebugPCC("\t\t skipping " + part + "\n");
            continue; // there is no rest particles in measurement
        }

        int const _iPart = GetPartTypeNumber(part);
        int iCent = 0;
        for (auto& cent : fstCentralities) {
//            DebugPCC("\t\tCent: " << cent << "\n");
//             CollisionSystem:ParticleType:CentNumber:Stat/Sys:Function:FunctionVar
            std::string stHistName{fstCollisionSystem};
            stHistName += part;
            stHistName += cent;
            stHistName += AliMCSpectraWeights::GetFunctionFromSysFlag(fFlag);
            stHistName += AliMCSpectraWeights::GetSysVarFromSysFlag(fFlag);
            if (fUseMBFractions) {
                stHistName = fstCollisionSystem + part + std::to_string(iCent) +
                "PublishedFractions";
                if (fstCollisionSystem == "pp")
                    stHistName = fstCollisionSystem + part + "MB";
            }
//            DebugPCC("\t\t\tLoading hist " << stHistName << "\n");
            TH1D* hist = (TH1D*)fMeasuredFile->Get(stHistName.c_str());
            if (!hist) {
                std::cerr << "AliMCSpectraWeights::Error: could not find "
                << stHistName << "\n";
                continue;
            }

            // hist-> pt:multcent:PartType
            std::array<float, 3> binEntry{0.};
            binEntry[2] = static_cast<float>(_iPart); // particle type
            if (!fUseMultiplicity)
                binEntry[1] = static_cast<float>(GetCentFromString(cent));
            else if (fUseMBFractions && "pp" == fstCollisionSystem) {
                binEntry[1] = 1.0; // only MB
            } else
                binEntry[1] = AliMCSpectraWeights::GetMultFromCent(cent);
//            DebugPCC("\t\t\tWriting to fHistDataFractions\n");
//            DebugPCC("\t\t\t part: " << binEntry[2] << "\t mult: " << binEntry[1] << "\n");
            for (int ipt = 0; ipt < fHistDataFractions->GetNbinsX(); ++ipt) {
                binEntry[0] = fHistDataFractions->GetXaxis()->GetBinCenter(ipt);
                if (binEntry[0] < 0.1){ // pT cut; measurements start at 0.15 at best
                    DebugPCC("pt is smaller than 0.1\n");
                    continue;
                }
//                DebugPCC("\t\t pT: " << binEntry[0]<<"\n");
                auto const _FractionValue =
                hist->GetBinContent(hist->FindBin(binEntry[0]));
//                DebugPCC("\t\t\t\t fraction: " << _FractionValue << "\n");
                FillTH3WithValue(fHistDataFractions, binEntry, _FractionValue);
            }
//            delete hist;
            ++iCent;
        }
    }
    DebugPCC("Loading measured fractions finished\n");
    fMeasuredFile->Close();
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    DebugChrono("LoadMeasuredFractions took " << duration << " microseconds\n");
#endif
}

/**
 *  @brief calculate the relative fractions of all particle species in MC
 *
 *  This require the input from a previous train output and is only calculated
 *  if the internal histogram of MC information is filled.
 */
bool AliMCSpectraWeights::CalcMCFractions() {
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    DebugPCC("Calculate MC fractions\n");
    if (!fHistMCGenPrimTrackParticle)
        return false;
    std::array<std::array<TH1D*, 10>, 20> _histMCFractions{
        nullptr}; // FIXME: WARNING HARD CODED RANGES
    std::array<TH1D*, 20> _h1pTMCAll{nullptr};
    fHistMCFractions->Reset(); // clean up for new input
    for (int icent = 0; icent < fNCentralities; ++icent) {
        if (fUseMBFractions && "pp" == fstCollisionSystem && icent > 0)
            continue;
        auto const multTuple = AliMCSpectraWeights::GetMultTupleFromCent(icent);
        auto _multFront = multTuple.front();
        auto _multBack = multTuple.back();

        if (fUseMBFractions && "pp" == fstCollisionSystem) {
            _multFront = 0;
            _multBack = 49.9;
        }

        //        DebugPCC("\t cent: " << icent << " => " << _multFront << " - "
        //                             << _multBack << "\n");

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
            //            DebugPCC("\t\t project spectra of " <<
            //            fstPartTypes[ipart] << "\n");
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
                //                DebugPCC("\t\t Clone spectra to
                //                _h1pTMCAll\n");
                _h1pTMCAll[icent] =
                (TH1D*)_histMCFractions[icent][ipart]->Clone(
                                                             Form("h1pTMCAll_%d", icent));
            } else {
                //                DebugPCC("\t\t Add spectra to _h1pTMCAll\n");
                _h1pTMCAll[icent]->Add(_histMCFractions[icent][ipart]);
            }
        }
        //        DebugPCC("\t spectra and _h1pTMCAll created and filled\n");
        // all hist calculated
        // ------------------
        // calculate fractions
        //        DebugPCC("\t Calculate fractions in MC now\n");
        for (int ipart = 0; ipart < fNPartTypes; ++ipart) {
            if (nullptr == _histMCFractions[icent][ipart]) {
                std::cerr << "AliMCSpectraWeights::ERROR could not calculate "
                "fraction\n";
                continue;
            }

            int const _iPart =
            AliMCSpectraWeights::GetPartTypeNumber(fstPartTypes[ipart]);
            //            DebugPCC("\n\n ipart: " << ipart << "\n\n");
            TH1D* h1MCFraction = (TH1D*)_histMCFractions[icent][ipart]->Clone(
                                                                              Form("%s_fraction", _histMCFractions[icent][ipart]->GetName()));
            h1MCFraction->Divide(_h1pTMCAll[icent]);
            // Set content of fractions to fHistMCFractions
            //            DebugPCC("\t Write MC fraction to
            //            fHistMCFractions\n");
            for (int ipt = 0; ipt < fHistMCFractions->GetNbinsX(); ++ipt) {
                float const pt =
                fHistMCFractions->GetXaxis()->GetBinCenter(ipt);
                if (pt < 0)
                    continue;
                // fHistMCFractions : pt-mult-ipart
                std::array<float, 3> binEntry{
                    pt, AliMCSpectraWeights::GetMultFromCent(icent),
                    static_cast<float>(_iPart)};
                if (fUseMBFractions && "pp" == fstCollisionSystem)
                    binEntry[1] = 1.0;
                auto const _FractionValue =
                h1MCFraction->GetBinContent(h1MCFraction->FindBin(pt));
                //                DebugPCC("\t\t pT: " << binEntry[0] << " mult:
                //                " << binEntry[1]
                //                                     << " part: " <<
                //                                     binEntry[2]
                //                                     << " val: " <<
                //                                     _FractionValue << "\n");
                FillTH3WithValue(fHistMCFractions, binEntry, _FractionValue);
            }
            delete h1MCFraction;
        }
        //        DebugPCC("\t Fractions in MC calculated\n");
    }
    //    DebugPCC("\t delete tmp histos\n");
    for (auto& histArray : _histMCFractions) {
        for (auto& hist : histArray) {
            if (hist) {
                delete hist;
                hist = nullptr;
            }
        }
    }
    for (auto& hist : _h1pTMCAll) {
        if (hist) {
            delete hist;
            hist = nullptr;
        }
    }
    DebugPCC("\t ...worked\n");
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    DebugChrono("CalcMCFractions took " << duration << " microseconds\n");
#endif
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
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    //    DebugPCC("Correct data fractions for not having rest particles
    //    measured\n");
    auto dMultHigh = 49;
    if (fstCollisionSystem.find("ppb") != std::string::npos)
        dMultHigh = 299;
    else if (fstCollisionSystem.find("pbpb") != std::string::npos)
        dMultHigh = 3499;


    auto _bin1 =
    fHistMCGenPrimTrackParticle->GetYaxis()->FindBin(0.);
    auto _bin2 =
    fHistMCGenPrimTrackParticle->GetYaxis()->FindBin(dMultHigh);

    auto h1pTMCAll = (TH1D*)fHistMCGenPrimTrackParticle->ProjectionX(
                                                                     "h1pTMCAll", _bin1, _bin2, 1,
                                                                     fHistMCGenPrimTrackParticle->GetNbinsZ(), "e");
    if (!h1pTMCAll) {
        std::cerr
        << "AliMCSpectraWeights::ERROR could not create h1pTMCAll\n";
        return false;
    }
    //        DebugPCC("\t created MC all spectra\n");
    auto const _iRestPos = AliMCSpectraWeights::GetPartTypeNumber("Rest");
    auto const _RestBin =
    fHistMCGenPrimTrackParticle->GetZaxis()->FindBin(_iRestPos);
    auto h1RestCorrFactor = fHistMCGenPrimTrackParticle->ProjectionX(
                                                                     "h1RestCorrFactor", _bin1, _bin2, 1, _RestBin - 1,
                                                                     "e");
    h1RestCorrFactor->Divide(h1pTMCAll);
    //        DebugPCC("\t calculated correction factor\n");
    for (int icent = 0; icent < fNCentralities; ++icent){
        for (int ipart = 0; ipart < fNPartTypes; ++ipart) {
            //            DebugPCC("\t\t correct " << fstPartTypes[ipart] <<
            //            "\n");
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
                if (fUseMBFractions && "pp" == fstCollisionSystem)
                    binEntry[1] = 1.0;
                //                DebugPCC("\t\t\t pT: " << binEntry[0]
                //                                       << " mult: " <<
                //                                       binEntry[1]
                //                                       << " part: " <<
                //                                       binEntry[2] << "\n");
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
    }
    delete h1pTMCAll;
    delete h1RestCorrFactor;
    DebugPCC("\t ...correction finished\n");
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    DebugChrono("CorrectFractionsforRest took " << duration
                << " microseconds\n");
#endif
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
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    DebugPCC("Calculate weight factors\n");
    // correction of rest particles not measured in data fractions (see
    // AnalysisNote)
    fHistMCWeights->Reset(); // clean up for new input
    for (int icent = 0; icent < fNCentralities; icent++) {
        if (fUseMBFractions && "pp" == fstCollisionSystem && icent > 0)
            continue;
//        DebugPCC("\t cent: " << icent << "\n");
        for (int ipart = 0; ipart < fNPartTypes; ipart++) {
//        DebugPCC("\t\t part: " << fstPartTypes[ipart] <<"\n");
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
                if (fUseMBFractions && "pp" == fstCollisionSystem){
                    binEntry[1] = 1.0;
                }
//                DebugPCC("\t\t\t pT: " << binEntry[0]
//                                       << " mult: " <<
//                                       binEntry[1]
//                                       << " part: " <<
//                                       binEntry[2] << "\n");

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

                float value = 1;
                if (dFractionMC != 0 && fstPartTypes[ipart] != "Rest" && fstPartTypes[ipart] != "Lambda") {
                    value = dFractionData / dFractionMC;
                }
//                DebugPCC("\t\t\t fractionMC: "
//                         << dFractionMC << "\t fractionData: "
//                         << dFractionData
//                         << "\t weight factor: " << value
//                         << "\n");
                if (fDoSystematics) {
                    fHistMCWeightsSys[fFlag]->SetBinContent(_iBinWeight, value);
                    fHistMCWeightsSys[fFlag]->SetBinError(_iBinWeight, 1e-30);
                } else {
                    fHistMCWeights->SetBinContent(_iBinWeight, value);
                    fHistMCWeights->SetBinError(_iBinWeight, 1e-30);
                }
            }
        }
    }

//#ifdef __AliMCSpectraWeights_DebugPCC__
//    DebugPCC("weight factors at pT=1 GeV and mult class 1\n");
//    for (int ipart = 0; ipart < fNPartTypes; ipart++) {
//        std::array<float, 3> binEntry{
//            1., static_cast<float>(GetMultFromCent(1)),
//            static_cast<float>(AliMCSpectraWeights::GetPartTypeNumber(fstPartTypes[ipart]))};
//        auto const _iBinWeight = fHistMCWeights->FindBin(
//                                                         binEntry[0], binEntry[1], binEntry[2]);
//        DebugPCC("part " << fstPartTypes[ipart] << " ");
//        double value = 1;
//        if (fDoSystematics) {
//            value = fHistMCWeightsSys[fFlag]->GetBinContent(_iBinWeight);
//        } else {
//            value = fHistMCWeights->GetBinContent(_iBinWeight);
//        }
//        DebugPCC(value << "\n");
//    }
//#endif

    DebugPCC("... calculation completed\n");
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    DebugChrono("CalculateMCWeights took " << duration << " microseconds\n");
#endif
    return true;
}

bool AliMCSpectraWeights::CalculateSystematicUncertainties(){

    DebugPCC("Calculating average systematic effect\n");
//    DebugPCC("fNCentralities is " << fNCentralities << "\n");
    for (int icent = 0; icent < fNCentralities; ++icent) { // cent loop
        DebugPCC("------------------------------------\n");
//        DebugPCC("\n\n\t cent:" << icent << " < " << _testNCent << "\n");
        auto const _mult = GetMultFromCent(icent);
//        DebugPCC("\t\t cent after:" << icent << "\n");

        for (auto const& part : fstPartTypes) { // particle loop
            DebugPCC("------------------\n");
            DebugPCC("\n\t\t part:" + part + "\n");
            if(part.find("Rest") != std::string::npos ||
               part.find("rest") != std::string::npos ||
               part.find("Lambda") != std::string::npos ||
               part.find("lambda") != std::string::npos){
                DebugPCC("\t\t\t skip " + part + " \n");
                continue; //skip rest
            }
            auto const _part = GetPartTypeNumber(part);

            for (int ipT= 1; ipT<fHistMCWeights->GetNbinsX()+1;++ipT){ // pt loop
                auto const _pT = fHistMCWeights->GetXaxis()->GetBinCenter(ipT);
                auto const _bin = fHistMCWeights->FindBin(_pT, _mult, _part);

                auto const _nomVal = fHistMCWeightsSys[AliMCSpectraWeights::SysFlag::kNominal]->GetBinContent(_bin);

                DebugPCC("\n\t\t\t pT: " + std::to_string(_pT) + " mult: " + std::to_string(_mult) + " part: " + std::to_string(_part) + "\n");
                DebugPCC("\t\t\t nominal weight: " + std::to_string(_nomVal) + "\n");

                if(_nomVal <= 0.){
                    DebugPCC("\t\t\t weight is too small\n");
                    fHistMCWeightsSysUp->SetBinContent(_bin, 1);
                    fHistMCWeightsSysDown->SetBinContent(_bin, 1);
                    continue;
                }

                double _rmsSys = 0; //rms of sum

//                for (auto& _sys : fAllSystematicFlags) {
                for (int iSys = 0; iSys < fNSysFlags; ++iSys){ // loop over all systematic variations
//                    DebugPCC("\t\t\t\t looking at sys " << GetFunctionFromSysFlag(fAllSystematicFlags[iSys])+GetSysVarFromSysFlag(fAllSystematicFlags[iSys]) << "\n");
                    if(!fHistMCWeightsSys.at(fAllSystematicFlags[iSys])) continue;
                    auto const _sysVal = fHistMCWeightsSys[fAllSystematicFlags[iSys]]->GetBinContent(_bin);
                    auto const _relSys = _sysVal/_nomVal;
                    auto const _absDiffSys = 1 - _relSys;
                    _rmsSys+=_absDiffSys*_absDiffSys;
                }
                _rmsSys = TMath::Sqrt(_rmsSys);
                DebugPCC("\t\t\t average sys: " + std::to_string(_rmsSys) + " percent\n");

                auto const _avgSysValUp = _nomVal * (1+_rmsSys);
                auto const _avgSysValDown = _nomVal * (1-_rmsSys);

                DebugPCC("\t\t\t sysUp: " + std::to_string(_avgSysValUp) + "\n");
                DebugPCC("\t\t\t sysDown: " + std::to_string(_avgSysValDown) + "\n");

                fHistMCWeightsSysUp->SetBinContent(_bin, _avgSysValUp);
                fHistMCWeightsSysDown->SetBinContent(_bin, _avgSysValDown);
            }
        }

    }
    DebugPCC("\t... done\n");
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
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
//    DebugPCC("Count event multiplicity ");
    fMultOrCent = 0;
    const Double_t lowPtCut = 0.05;
    const Double_t eta = 0.5;
    //    if (fstCollisionSystem.find("pp") != std::string::npos)
    //        eta = 0.5;
//    DebugPCC("in eta < " << eta << " and pT > " << lowPtCut << "\n");
    if (!fMCEvent)
        return;
    for (int ipart = 0; ipart < fMCEvent->GetNumberOfTracks(); ipart++) {
        AliMCParticle* fMCParticle = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(ipart));
        if (!fMCParticle) {
            std::cerr << "AliMCSpectraWeights::ERROR::noMCParticle\n";
            continue;
        }
        if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(ipart,
                                                                      fMCEvent))
            continue;
        if (!fMCEvent->IsPhysicalPrimary(ipart))
            continue; // secondary rejection
        if (TMath::Abs(fMCParticle->Charge()) < 0.01)
            continue; // neutral rejection
        float pEta = fMCParticle->Eta();
        if (TMath::Abs(pEta) > eta)
            continue; // acceptance cut
        if (fMCParticle->Pt() < lowPtCut)
            continue; // TODO: hard coded low pT cut
        ++fMultOrCent;
    }

//    DebugPCC("... counted " << fMultOrCent << " charged particles\n");
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    DebugChrono("CountEventMult took " + std::to_string(duration) + " microseconds\n");
#endif
}

/**
 *  @brief select randomly a systematic variation
 */
void AliMCSpectraWeights::SelectRndSysFlagForEvent() {
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    int const _SysIndex =
    static_cast<int>(std::round(frndGen.Uniform(-0.5, fNSysFlags - 0.5)));
    fFlag = fAllSystematicFlags[_SysIndex];
    DebugPCC(
             "SysFlag: " + std::to_string(_SysIndex) + " "
             + GetFunctionFromSysFlag(fAllSystematicFlags[_SysIndex])
             + GetSysVarFromSysFlag(fAllSystematicFlags[_SysIndex])
             + "\n");
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    DebugChrono("SelectRndSysFlagForEvent took " << duration
                << " microseconds\n");
#endif
}

/**
 *  @brief
 *  @param[in] mcGenParticle
 *  @param[in] mcEvent
 *  @return
 */
float const AliMCSpectraWeights::GetMCSpectraWeight(Int_t mcGenParticle,
                                                    AliMCEvent* mcEvent) {
    return AliMCSpectraWeights::GetMCSpectraWeightNominal(mcGenParticle);
}

int const AliMCSpectraWeights::CheckAndIdentifyParticle(TParticle* part) {
    if (!part->GetPDG()) {
        DebugPCC("Warning: particle has no PDG; skipped\n");
        return -2;
    }
    if (TMath::Abs(part->GetPDG()->Charge()) < 0.01 && TMath::Abs(part->GetPdgCode()!=3122)) {
        DebugPCC("Warning: particle not charged\n");
        return -3;
    }
    return AliMCSpectraWeights::IdentifyMCParticle(part);
}

std::vector<int> const AliMCSpectraWeights::FindBinEntry(float pt, int const part) {
    fMultClass = AliMCSpectraWeights::GetCentFromMult(fMultOrCent); // nominal cent value
    if (pt < 0.15) {
        DebugPCC("Warning: pt too low; pt = " + std::to_string(pt) + "\n");
        return {-1};
    }
    if (pt >= 20) {
        DebugPCC("Info: pt too high; pt = " + std::to_string(pt) + "; set to 19.9\n");
        pt = 19.9;
    }
    auto const _multTuple = AliMCSpectraWeights::GetMultTupleFromCent(fMultClass);

    std::vector<int> results;
    auto const icent_low    = AliMCSpectraWeights::GetCentFromMult(_multTuple.front());
    auto const icent_high   = AliMCSpectraWeights::GetCentFromMult(_multTuple.back());


    if(fDoInterpolation){
        std::array<float, 3> binEntry_low{
            pt, static_cast<float>(AliMCSpectraWeights::GetMultFromCent(icent_low)),
            static_cast<float>(part)};
        auto const _iBin_low = GetBinFromTH3(
                                             fHistMCWeightsSys[AliMCSpectraWeights::SysFlag::kNominal], binEntry_low);

        std::array<float, 3> binEntry_high{
            pt, static_cast<float>(AliMCSpectraWeights::GetMultFromCent(icent_high)),
            static_cast<float>(part)};
        auto const _iBin_high = GetBinFromTH3(
                                              fHistMCWeightsSys[AliMCSpectraWeights::SysFlag::kNominal], binEntry_high);
        results.push_back(_iBin_low);
        results.push_back(_iBin_high);
    } else {
        std::array<float, 3> binEntry{
            pt, static_cast<float>(AliMCSpectraWeights::GetMultFromCent(fMultClass)),
            static_cast<float>(part)};
        auto const _iBin = GetBinFromTH3(
                                             fHistMCWeightsSys[AliMCSpectraWeights::SysFlag::kNominal], binEntry);
        results.push_back(_iBin);
    }

    DebugPCC("Found bin at low: " + std::to_string(_iBin_low) << "\t high: " << std::to_string(_iBin_high) << "\n");
    return results;
}

float const
AliMCSpectraWeights::GetMCSpectraWeightNominal(TParticle* mcGenParticle) {
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    if (fbTaskStatus < AliMCSpectraWeights::TaskState::kMCWeightCalculated) {
        DebugPCC("Warning: Status not kMCWeightCalculated\n");
        return 1;
    }
    int const particleType =
    AliMCSpectraWeights::CheckAndIdentifyParticle(mcGenParticle);
    if (particleType < 0) {
        DebugPCC("Can't find particle type\n");
        return 1;
    }
    auto const _iBin =
    AliMCSpectraWeights::FindBinEntry(mcGenParticle->Pt(), particleType);
    if (_iBin.size() < 1) {
        DebugPCC("Can't find bin\n");
        return 1;
    }
    float _weight_Interpolated = 1;
    if(fDoInterpolation){
        auto const _MultTuple = AliMCSpectraWeights::GetMultTupleFromCent(fMultClass);
        float const weight1 = fHistMCWeightsSys[AliMCSpectraWeights::SysFlag::kNominal]->GetBinContent(_iBin.front());
        float const weight2 = fHistMCWeightsSys[AliMCSpectraWeights::SysFlag::kNominal]->GetBinContent(_iBin.back());
        _weight_Interpolated = weight1 + (weight2 - weight1)/(_MultTuple.back() - _MultTuple.front()) * (fMultOrCent - _MultTuple.front()); // linear interpolation
    }
    else {
        _weight_Interpolated = fHistMCWeightsSys[AliMCSpectraWeights::SysFlag::kNominal]->GetBinContent(_iBin.front());
    }

    if (_weight_Interpolated <= 0) {
        DebugPCC("ERROR: negative weight " << _weight_Interpolated << "; set to 1\n");
        _weight_Interpolated = 1;
    }
    DebugPCC("GetMCSpectraWeight: nominal");
    DebugPCC(fstPartTypes[particleType] << " ");
    DebugPCC("pT: " + std::to_string(mcGenParticle->Pt()) + " ");
    DebugPCC("weight: " + std::to_string(_weight_Interpolated) + "\n");
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    DebugChrono("GetMCSpectraWeightNominal took " << duration
                << " microseconds\n");
#endif
    return _weight_Interpolated;
}

float const
AliMCSpectraWeights::GetMCSpectraWeightNominal(Int_t mcGenParticle) {
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    if (fbTaskStatus < AliMCSpectraWeights::TaskState::kMCWeightCalculated) {
        DebugPCC("Warning: Status not kMCWeightCalculated\n");
        return 1;
    }
    AliMCParticle* _MCpart = (AliMCParticle*)fMCEvent->GetTrack(mcGenParticle);
    return AliMCSpectraWeights::GetMCSpectraWeightNominal(_MCpart->Particle());
}

float const
AliMCSpectraWeights::GetMCSpectraWeightSystematics(TParticle* mcGenParticle, Int_t SysCase) {
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    if (fbTaskStatus < AliMCSpectraWeights::TaskState::kMCWeightCalculated) {
        DebugPCC("Warning: Status not kMCWeightCalculated\n");
        return 1;
    }
    int const particleType =
    AliMCSpectraWeights::CheckAndIdentifyParticle(mcGenParticle);
    if (particleType < 0) {
        DebugPCC("Can't find particle type\n");
        return 1;
    }
    auto const _iBin =
    AliMCSpectraWeights::FindBinEntry(mcGenParticle->Pt(), particleType);
    if (_iBin.size() < 1) {
        DebugPCC("Can't find bin\n");
        return 1;
    }
    float _weight_Interpolated = 1 ;
    auto const _MultTuple = AliMCSpectraWeights::GetMultTupleFromCent(fMultClass);

    if(fDoInterpolation){
        if(SysCase > 0){
            //weight = fHistMCWeightsSysUp->GetBinContent(_iBin);
            float const weight1 = fHistMCWeightsSysUp->GetBinContent(_iBin.front());
            float const weight2 = fHistMCWeightsSysUp->GetBinContent(_iBin.back());
            _weight_Interpolated = weight1 + (weight2 - weight1)/(_MultTuple.back() - _MultTuple.front()) * (fMultOrCent - _MultTuple.front()); // linear interpolation

        }
        if(SysCase < 0){
            //        weight = fHistMCWeightsSysDown->GetBinContent(_iBin);
            float const weight1 = fHistMCWeightsSysDown->GetBinContent(_iBin.front());
            float const weight2 = fHistMCWeightsSysDown->GetBinContent(_iBin.back());
            _weight_Interpolated = weight1 + (weight2 - weight1)/(_MultTuple.back() - _MultTuple.front()) * (fMultOrCent - _MultTuple.front()); // linear interpolation
        }
    }
    else {
        if(SysCase > 0){
            //weight = fHistMCWeightsSysUp->GetBinContent(_iBin);
            _weight_Interpolated = fHistMCWeightsSysUp->GetBinContent(_iBin.front());
        }
        if(SysCase < 0){
            //        weight = fHistMCWeightsSysDown->GetBinContent(_iBin);
            _weight_Interpolated = fHistMCWeightsSysDown->GetBinContent(_iBin.front());
        }
    }

    if (_weight_Interpolated <= 0) {
        DebugPCC("ERROR: negative weight; set to 1\n");
        _weight_Interpolated = 1;
    }
    DebugPCC("GetMCSpectraWeight: with systematics");
    DebugPCC(fstPartTypes[particleType]+ " ");
    DebugPCC("pT: " + std::to_string(mcGenParticle->Pt()) + " ");
    DebugPCC("weight: " + std::to_string(_weight_Interpolated) + "\n");
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    DebugChrono("GetMCSpectraWeightSystematics took " << duration
                << " microseconds\n");
#endif
    return _weight_Interpolated;
}

float const
AliMCSpectraWeights::GetMCSpectraWeightSystematics(Int_t mcGenParticle, Int_t SysCase) {
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    if (fbTaskStatus < AliMCSpectraWeights::TaskState::kMCWeightCalculated) {
        DebugPCC("Warning: Status not kMCWeightCalculated\n");
        return 1;
    }
    AliMCParticle* _MCpart = (AliMCParticle*)fMCEvent->GetTrack(mcGenParticle);
    return AliMCSpectraWeights::GetMCSpectraWeightSystematics(_MCpart->Particle(), SysCase);
}

float const
AliMCSpectraWeights::GetMCSpectraWeight(TParticle* mcGenParticle,
                   Int_t SysCase){
    if(0==SysCase){
        return AliMCSpectraWeights::GetMCSpectraWeightNominal(mcGenParticle);
    }
    else{
        return AliMCSpectraWeights::GetMCSpectraWeightSystematics(mcGenParticle, SysCase);
    }
}

float const
AliMCSpectraWeights::GetMCSpectraWeight(Int_t mcGenParticle,
                   Int_t SysCase){
    if(0==SysCase){
        return AliMCSpectraWeights::GetMCSpectraWeightNominal(mcGenParticle);
    }
    else{
        return AliMCSpectraWeights::GetMCSpectraWeightSystematics(mcGenParticle, SysCase);
    }
}

/**
 *  @brief identify secondary particle depending on mother pid
 *  @param[in] part
 *  @return -1 for error; 0 from lambda, 1 from kaon0Short, 2 other
 */
int const AliMCSpectraWeights::IdentifySecondaryType(Int_t partLabel){
    DebugPCC("IdentifySecondaryType\n");
    AliMCParticle* part = (AliMCParticle*)fMCEvent->GetTrack(partLabel);
    if(!part){
        std::cerr << "AliMCSpectraWeights::Error: particle pointer zero\n";
        return -1;
    }
    // two possible cases:
    // 1) if a proton or a pi- check if mother is a Lambda (and check that it's a physical primary)
    // 2) K0short case: if a pi+ check if mother is a K0short and assign the weight of K+, the same way as for Lambda's

    auto const absPDG = TMath::Abs(part->PdgCode());
    auto motherPartLabel = fMCEvent->GetLabelOfParticleMother(partLabel);
    if (motherPartLabel<0) {
        std::cerr << "AliMCSpectraWeights::Error: mother label negative\n";
        return -1;
    }
    auto const motherPart =  (AliMCParticle*)fMCEvent->GetTrack(motherPartLabel);
    if(!motherPart){
        std::cerr << "AliMCSpectraWeights::Error: could not find mother\n";
        return -1;
    }
    auto const motherPDG = TMath::Abs(motherPart->PdgCode());
    auto const grandMotherLabel = fMCEvent->GetLabelOfParticleMother(motherPartLabel);
    DebugPCC("\t\t part = " + std::to_string(absPDG) + " \t mother = " + std::to_string(motherPDG)+"\n");

    if ((motherPDG == 3122 || motherPDG == 3222 || motherPDG == 3112 || motherPDG == 3212) && fMCEvent->IsPhysicalPrimary(motherPartLabel)) { //&& motherPart->IsPhysicalPrimary()
        DebugPCC("\t\t mother is lambda case\n");
        return 0;
    }// end if primary lambda

    if ( (motherPDG == 310 || motherPDG == 130 || motherPDG == 311 || motherPDG == 321) && (fMCEvent->IsPhysicalPrimary(motherPartLabel) || fMCEvent->IsPhysicalPrimary(grandMotherLabel)) ) {
        DebugPCC("\t\t mother is from K0short case\n");
        return 1;
    }// end if primary K0short

    if(motherPDG == 211 && fMCEvent->IsPhysicalPrimary(motherPartLabel)){
        DebugPCC("\t\t mother is from pion case\n");
        return 2;
    }
    if(motherPDG == 3122 && TMath::Abs(fMCEvent->MotherOfParticle(motherPartLabel)->PdgCode()) == 3312){
        DebugPCC("\t\t mother is from xi case\n");
        return 3;
    }

#ifdef __AliMCSpectraWeights_DebugPCC__
    std::cerr << "\tclimbing up the tree\n";
    bool hasMother=true;
    while (hasMother) {
        motherPartLabel = fMCEvent->GetLabelOfParticleMother(motherPartLabel);
        auto const _motherPart =  (AliMCParticle*)fMCEvent->GetTrack(motherPartLabel);
        if (!_motherPart) {
            std::cerr << "\t no more mother\n";
            hasMother = false;
            break;
        } else {
            auto const _motherPDG = _motherPart->PdgCode();
            std::cerr << "\t mother: " << std::to_string(_motherPDG);
        }
    }
#endif
    return 4;
}

/**
 *  @brief weight factor dependent of mother particle
 *  @param[in] mcGenParticle
 *  @return weight factor for secondary particle
 */
float const AliMCSpectraWeights::GetWeightForSecondaryParticle(Int_t partLabel, Int_t SysCase){
    DebugPCC("\n\nGetWeightForSecondaryParticle\n");
    float _weight_Interpolated = 1;
    if (fbTaskStatus < AliMCSpectraWeights::TaskState::kMCWeightCalculated) {
        DebugPCC("Warning: Status not kMCWeightCalculated\n");
        return 1;
    }

    auto const _SecondaryID = IdentifySecondaryType(partLabel);
    if(_SecondaryID < 0){
        std::cerr << "AliMCSpectraWeights::Error: secondary ID error\n";
        return 1;
    }
    DebugPCC("\t SecondaryID = " + std::to_string(_SecondaryID)+"\n");
    AliMCParticle* part = (AliMCParticle*)fMCEvent->GetTrack(partLabel);
    if(!part){
        std::cerr << "AliMCSpectraWeights::Error: can not found particle for label " << partLabel << "\n";
        return 1;
    }
    auto motherPart = fMCEvent->MotherOfParticle(partLabel);
//    auto const motherPart =  (AliMCParticle*)fMCEvent->GetTrack(motherPartLabel);
    if(!motherPart){
        std::cerr << "AliMCSpectraWeights::Error: mother is not available\n";
        return 1;
    }
    std::vector<int> _iBin = {};
    switch (_SecondaryID) {
        case 0: // Lambda case
            _iBin = AliMCSpectraWeights::FindBinEntry(motherPart->Pt(), AliMCSpectraWeights::ParticleType::kSigmaPlus);
            break;
        case 1: // Kaon0Short case
            _iBin = AliMCSpectraWeights::FindBinEntry(motherPart->Pt(), AliMCSpectraWeights::ParticleType::kKaon);
            break;
        case 2: // electron from primary pion
            _iBin = AliMCSpectraWeights::FindBinEntry(motherPart->Pt(), AliMCSpectraWeights::ParticleType::kPion);
            break;
        case 3: // xi -> lambda -> proton
            _iBin = AliMCSpectraWeights::FindBinEntry(motherPart->Pt(), AliMCSpectraWeights::ParticleType::kSigmaPlus);
            break;
#ifdef __AliMCSpectraWeights_DebugPCC__
        case 4:
            DebugPCC("GetWeightForSecondaryParticle::case 4\n");
            DebugPCC("\t PID: "+std::to_string(part->PdgCode())+"\t mother: "+std::to_string(motherPart->PdgCode())+"\n");
            break;
#endif
        default:
            break;
    }
    if (_iBin.size() < 1) {
        DebugPCC("\tCan't find bin;\n");
        return 1;
    }

    auto const _MultTuple = AliMCSpectraWeights::GetMultTupleFromCent(fMultClass);

    if(SysCase==0){
//        weight = fHistMCWeightsSys[AliMCSpectraWeights::SysFlag::kNominal]->GetBinContent(_iBin);
        if(fDoInterpolation){
            float const weight1 = fHistMCWeightsSys[AliMCSpectraWeights::SysFlag::kNominal]->GetBinContent(_iBin.front());
            float const weight2 = fHistMCWeightsSys[AliMCSpectraWeights::SysFlag::kNominal]->GetBinContent(_iBin.back());
            _weight_Interpolated = weight1 + (weight2 - weight1)/(_MultTuple.back() - _MultTuple.front()) * (fMultOrCent - _MultTuple.front()); // linear interpolation
        } else {
            _weight_Interpolated = fHistMCWeightsSys[AliMCSpectraWeights::SysFlag::kNominal]->GetBinContent(_iBin.front());
        }
    } else if(SysCase<0){
//        weight = fHistMCWeightsSysDown->GetBinContent(_iBin);
        if(fDoInterpolation){
            float const weight1 = fHistMCWeightsSysDown->GetBinContent(_iBin.front());
            float const weight2 = fHistMCWeightsSysDown->GetBinContent(_iBin.back());
            _weight_Interpolated = weight1 + (weight2 - weight1)/(_MultTuple.back() - _MultTuple.front()) * (fMultOrCent - _MultTuple.front()); // linear interpolation
        } else {
            _weight_Interpolated = fHistMCWeightsSysDown->GetBinContent(_iBin.front());
        }
    } else if(SysCase>0){
//        weight = fHistMCWeightsSysUp->GetBinContent(_iBin);
        if(fDoInterpolation){
            float const weight1 = fHistMCWeightsSysUp->GetBinContent(_iBin.front());
            float const weight2 = fHistMCWeightsSysUp->GetBinContent(_iBin.back());
            _weight_Interpolated = weight1 + (weight2 - weight1)/(_MultTuple.back() - _MultTuple.front()) * (fMultOrCent - _MultTuple.front()); // linear interpolation
        } else {
            _weight_Interpolated = fHistMCWeightsSysUp->GetBinContent(_iBin.front());
        }
    }

    if(_weight_Interpolated < 1e-2){
        DebugPCC("AliMCSpectraWeights::WARNING: weight is too small -> set to 1 \n");
        _weight_Interpolated = 1;
    }

    DebugPCC("\t final weight is " + std::to_string(_weight_Interpolated) + "\n");
    return _weight_Interpolated;
}

void AliMCSpectraWeights::StartNewEvent() {
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    AliMCSpectraWeights::CountEventMult();
//    if (fDoSystematics)
//        AliMCSpectraWeights::SelectRndSysFlagForEvent();
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    DebugChrono("StartNewEvent took " << duration << " microseconds\n");
#endif
}

/**
 *  @brief
 *  @param[in] mcEvent
 *
 */
void AliMCSpectraWeights::FillMCSpectra(AliMCEvent* mcEvent) {
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t1 = std::chrono::high_resolution_clock::now();
#endif
    DebugPCC("FillMCSpectra\n");
    if (fbTaskStatus >= AliMCSpectraWeights::TaskState::kMCSpectraObtained)
        return;
    if (mcEvent != fMCEvent) {
        fMCEvent = mcEvent;
        AliMCSpectraWeights::CountEventMult();
    }

    bool const ispPbCollision = fstCollisionSystem.find("ppb")!=std::string::npos;

    for (int iParticle = 0; iParticle < fMCEvent->GetNumberOfTracks(); ++iParticle) {
        AliMCParticle* fMCParticle = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(iParticle));
        if (!fMCParticle) {
            std::cerr << "AliMCSpectraWeights::ERROR::noMCParticle\n";
            continue;
        }
        if (!fMCEvent->IsPhysicalPrimary(iParticle))
            continue;
        if (TMath::Abs(fMCParticle->Charge()) < 0.01 && TMath::Abs(fMCParticle->PdgCode())!=3122)
            continue;
        if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iParticle, mcEvent)) // is from MC pile-up
            continue;

        float partY = fMCParticle->Y();
        float _maxY = 0.5; // hard coded max eta; in all papers 0.5

        //for p-Pb this is 0 < y_cms < 0.5 with y_traf = -0.465
        //-->
        if(ispPbCollision){
            if( partY < 0.465 || partY > 0.965)
                continue;
        } else { // symmetric collision systems
            if (TMath::Abs(partY) > _maxY)
                continue; // apply same acceptance as in published spectra
        }
        int particleType =
        AliMCSpectraWeights::IdentifyMCParticle(fMCParticle->Particle());
        if (particleType < 0)
            continue;
        //        std::array<float, 3>
        //        binEntry{static_cast<float>(mcGenParticle->Pt()),
        //                                      fMultOrCent,
        //                                      static_cast<float>(particleType)};
        DebugPCC("\t Fill with:\n");
        DebugPCC("\t\tpT: " + std::to_string(fMCParticle->Pt()) + "\t mult: " + std::to_string(fMultOrCent) + "\t particle: " + std::to_string(particleType) + "\n");
        fHistMCGenPrimTrackParticle->Fill(
                                          static_cast<float>(fMCParticle->Pt()), fMultOrCent,
                                          static_cast<float>(particleType));
    }
#ifdef __AliMCSpectraWeights_DebugTiming__
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    DebugChrono("FillMCSpectra took " << duration << " microseconds\n");
#endif
}

/**
 *  @brief
 *  @param[in] mcParticle
 *  @return
 */
int const AliMCSpectraWeights::IdentifyMCParticle(TParticle* mcParticle) {
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
    if(ipdg == 3122)
        return GetPartTypeNumber(
                                 AliMCSpectraWeights::ParticleType::kLambda);
    return GetPartTypeNumber(AliMCSpectraWeights::ParticleType::kRest);
}

/**
 *  @brief
 *  @param[in] CentBin
 *  @return
 */
float const AliMCSpectraWeights::GetMultFromCent(int CentBin) const {
    if (fstCollisionSystem == "pp") {
        // for | eta | < 0.5
        // pp 7 TeV
        /* switch (CentBin) {
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
        }*/
        // pp 13 TeV
        switch (CentBin) {
            case 0:
                return 26.02;
            case 1:
                return 20.02;
            case 2:
                return 16.17;
            case 3:
                return 13.77;
            case 4:
                return 12.04;
            case 5:
                return 10.02;
            case 6:
                return 7.95;
            case 7:
                return 6.32;
            case 8:
                return 4.50;
            case 9:
                return 2.55;
            default:
                return -2.0;
        }
    } else if (fstCollisionSystem == "ppb") {
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
    } else if (fstCollisionSystem == "pbpb" || fstCollisionSystem == "xexe") {
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
    }
//    else if (fstCollisionSystem == "xexe") {
//        // for | eta | < 0.5
//        switch (CentBin) {
//            case 0:
//                return 1167;
//            case 1:
//                return 939;
//            case 2:
//                return 706;
//            case 3:
//                return 478;
//            case 4:
//                return 315;
//            case 5:
//                return 198;
//            case 6:
//                return 118;
//            case 7:
//                return 64.7;
//            case 8:
//                return 32;
//            default:
//                return -2.0;
//        }
//    }

    return -1;
}

/**
 *  @brief
 *  @param[in] cent
 *  @return
 */
float const
AliMCSpectraWeights::GetMultFromCent(std::string const& cent) const {
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
        if (fstCollisionSystem.find("pp") != std::string::npos &&
            fstCollisionSystem.find("ppb") == std::string::npos)
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

    if (dMultLow < dMultHigh) {
        return {dMultLow, dMultHigh};
    }
    return {dMultHigh, dMultLow};
}

/**
 *  @brief
 *  @param[in] cent
 *  @return
 */
int const
AliMCSpectraWeights::GetCentFromString(std::string const& cent) const {
    for (int i = 0; i < static_cast<int>(fstCentralities.size()); ++i) {
        if (fstCentralities[i] == cent)
            return i;
    }
    return -1;
}

/**
 *  @brief
 *  @param[in] dMult
 *  @return
 */
float const AliMCSpectraWeights::GetCentFromMult(float const dMult) const {
    if (fstCollisionSystem.find("pp") != std::string::npos && fstCollisionSystem.find("ppb") == std::string::npos) {
/* pp7TeV
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
*/
// pp13TeV
        if (dMult > 23.02)
            return 0;
        if (dMult > 18.09)
            return 1;
        if (dMult > 14.97)
            return 2;
        if (dMult > 12.91)
            return 3;
        if (dMult > 11.03)
            return 4;
        if (dMult > 8.99)
            return 5;
        if (dMult > 7.14)
            return 6;
        if (dMult > 5.41)
            return 7;
        if (dMult > 3.53)
            return 8;
        if (dMult > 0)
            return 9;
    } else if (fstCollisionSystem.find("ppb") !=
               std::string::npos)
    {
        if (dMult > 40.6)
            return 0;
        if (dMult > 33.35)
            return 1;
        if (dMult > 26.85)
            return 2;
        if (dMult > 19.65)
            return 3;
        if (dMult > 12.9)
            return 4;
        if (dMult > 7.1)
            return 5;
        if (dMult > 0)
            return 6;
    } else if (fstCollisionSystem.find("pbpb") != std::string::npos) {
        if (dMult > 1447)
            return 0;
        if (dMult > 1130)
            return 1;
        if (dMult > 807)
            return 2;
        if (dMult > 537)
            return 3;
        if (dMult > 343)
            return 4;
        if (dMult > 205)
            return 5;
        if (dMult > 112)
            return 6;
        if (dMult > 55)
            return 7;
        if (dMult > 0)
            return 8;
    } else if (fstCollisionSystem.find("xexe") != std::string::npos) {

        if (dMult > 1053)
            return 0;
        if (dMult > 822)
            return 1;
        if (dMult > 592)
            return 2;
        if (dMult > 396)
            return 3;
        if (dMult > 256)
            return 4;
        if (dMult > 158)
            return 5;
        if (dMult > 91.2)
            return 6;
        if (dMult > 48.3)
            return 7;
        if (dMult > 0)
            return 8;
    }

    return -1;
}

/**
 *  @brief
 *  @param[in] type
 *  @return
 */
int const AliMCSpectraWeights::GetPartTypeNumber(
                                                 AliMCSpectraWeights::ParticleType type) const {
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
        case AliMCSpectraWeights::ParticleType::kLambda:
            return AliMCSpectraWeights::ParticleType::kLambda;
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
int const
AliMCSpectraWeights::GetPartTypeNumber(std::string const& Particle) const {
    if (Particle == "Pion") {
        return AliMCSpectraWeights::GetPartTypeNumber(
                                                      AliMCSpectraWeights::ParticleType::kPion);
    } else if (Particle == "Proton") {
        return AliMCSpectraWeights::GetPartTypeNumber(
                                                      AliMCSpectraWeights::ParticleType::kProtons);
    } else if (Particle == "Kaon") {
        return AliMCSpectraWeights::GetPartTypeNumber(
                                                      AliMCSpectraWeights::ParticleType::kKaon);
    } else if (Particle == "SigmaMinus") {
        return AliMCSpectraWeights::GetPartTypeNumber(
                                                      AliMCSpectraWeights::ParticleType::kSigmaMinus);
    } else if (Particle == "SigmaPlus") {
        return AliMCSpectraWeights::GetPartTypeNumber(
                                                      AliMCSpectraWeights::ParticleType::kSigmaPlus);
    } else if (Particle == "Rest") {
        return AliMCSpectraWeights::GetPartTypeNumber(
                                                      AliMCSpectraWeights::ParticleType::kRest);
    } else
        return -1;
}

/**
 *  @brief
 *  @param[in] flag
 *  @return
 */
std::string const
AliMCSpectraWeights::GetFunctionFromSysFlag(SysFlag flag) const {
    std::string _default{"Bylinkin"};
    if (fstCollisionSystem == "pbpb") {
        _default = "Blastwave";
    }

    switch (flag) {
        case SysFlag::kNominal:
            return _default;
        case SysFlag::kBylinkinLower:
            return "BylinkinLower";
        case SysFlag::kBylinkinUpper:
            return "BylinkinUpper";
        case SysFlag::kHagedorn:
            return "Hagedorn";
        case SysFlag::kHagedornUpper:
            return "HagedornUpper";
        case SysFlag::kHagedornLower:
            return "HagedornLower";
        case SysFlag::kExponential:
            return "Exponential";
        case SysFlag::kExponentialUpper:
            return "ExponentialUpper";
        case SysFlag::kExponentialLower:
            return "ExponentialLower";
        case SysFlag::kBlastwave:
            return "Blastwave";
        case SysFlag::kBlastwaveUpper:
            return "BlastwaveUpper";
        case SysFlag::kBlastwaveLower:
            return "BlastwaveLower";
        default:
            return _default;
    }

    return _default;
}
std::string const
AliMCSpectraWeights::GetSysVarFromSysFlag(SysFlag flag) const {
    switch (flag) {
        case SysFlag::kPionUp:
            return "PionUp";
        case SysFlag::kPionDown:
            return "PionDown";
        case SysFlag::kProtonUp:
            return "ProtonUp";
        case SysFlag::kProtonDown:
            return "ProtonDown";
        case SysFlag::kKaonUp:
            return "KaonUp";
        case SysFlag::kKaonDown:
            return "KaonDown";
        case SysFlag::kSigmaPlusUp:
            return "SigmaPlusUp";
        case SysFlag::kSigmaPlusDown:
            return "SigmaPlusDown";
        case SysFlag::kSigmaMinusUp:
            return "SigmaMinusUp";
        case SysFlag::kSigmaMinusDown:
            return "SigmaMinusDown";
        default:
            return "";
    }

    return "";
}

AliMCSpectraWeightsHandler::AliMCSpectraWeightsHandler() : TNamed() {
    fMCSpectraWeight = nullptr;
}

AliMCSpectraWeightsHandler::AliMCSpectraWeightsHandler(
                                                       AliMCSpectraWeights* fMCWeight, const char* name)
: TNamed(name, name) {
    fMCSpectraWeight = fMCWeight;
}
