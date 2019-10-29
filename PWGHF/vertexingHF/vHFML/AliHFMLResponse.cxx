// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponse
// \brief helper class to handle application of ML models trained with python libraries
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <TGrid.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TFile.h>

#include "yaml-cpp/yaml.h"

#include "AliHFMLResponse.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliHFMLResponse);
/// \endcond

//________________________________________________________________
AliHFMLResponse::AliHFMLResponse() : TObject(),
                                     fModels{},
                                     fModelLibraries{},
                                     fModelVarNames{},
                                     fModelPaths{},
                                     fModelOutputCuts{},
                                     fPtBinsModel{},
                                     fPtBinCand(-1),
                                     fVars{}
{
    //
    // Default constructor
    //
}

//________________________________________________________________
AliHFMLResponse::AliHFMLResponse(string configfilename) : TObject(),
                                                          fModels{},
                                                          fModelLibraries{},
                                                          fModelVarNames{},
                                                          fModelPaths{},
                                                          fModelOutputCuts{},
                                                          fPtBinsModel{},
                                                          fPtBinCand(-1),
                                                          fVars{}
{
    //
    // Standard constructor
    //

    if (configfilename != "")
        SetConfigFile(configfilename);
}

//________________________________________________________________
AliHFMLResponse::~AliHFMLResponse()
{
    //
    // Destructor
    //
}

//--------------------------------------------------------------------------
AliHFMLResponse::AliHFMLResponse(const AliHFMLResponse &source) : TObject(source),
                                                                  fModels(source.fModels),
                                                                  fModelLibraries(source.fModelLibraries),
                                                                  fModelVarNames(source.fModelVarNames),
                                                                  fModelPaths(source.fModelPaths),
                                                                  fModelOutputCuts(source.fModelOutputCuts),
                                                                  fPtBinsModel(source.fPtBinsModel),
                                                                  fPtBinCand(source.fPtBinCand),
                                                                  fVars(source.fVars)
{
    //
    // Copy constructor
    //
}

AliHFMLResponse &AliHFMLResponse::operator=(const AliHFMLResponse &source)
{
    //
    // assignment operator
    //
    if (&source == this)
        return *this;

    TObject::operator=(source);

    fModels = source.fModels;
    fModelLibraries = source.fModelLibraries;
    fModelVarNames = source.fModelVarNames;
    fModelPaths = source.fModelPaths;
    fModelOutputCuts = source.fModelOutputCuts;
    fPtBinsModel = source.fPtBinsModel;
    fPtBinCand = source.fPtBinCand;
    fVars = source.fVars;

    return *this;
}

//________________________________________________________________
void AliHFMLResponse::SetConfigFile(const string configfilename)
{
    string configFilePath = GetFile(configfilename);
    YAML::Node configFile = YAML::LoadFile(configFilePath.data());
    if (configFile.IsNull())
        AliFatal("Yaml config file not found! Exit");

    fModelVarNames = configFile["VarNames"].as<vector<string> >();
    fModelPaths = configFile["ModelNames"].as<vector<string> >();
    fModelOutputCuts = configFile["ModelOutputCuts"].as<vector<double> >();
    fPtBinsModel = configFile["PtBins"].as<vector<double> >();
    fModelLibraries = configFile["ModelLibraries"].as<vector<string> >();

    //for consistency check
    unsigned int numModels = configFile["NumModels"].as<unsigned int>();
    unsigned int numVars = configFile["NumVars"].as<unsigned int>();

    if (numModels != fModelPaths.size() || numModels != fModelLibraries.size() || numModels != fModelOutputCuts.size() || numModels != fPtBinsModel.size() - 1)
        AliFatal("Inconsistency found in the number of models loaded from your yaml config file, please check it! Exit");

    if (numVars != fModelVarNames.size())
        AliFatal("Inconsistency found in the number of variables (features) loaded from your yaml config file, please check it! Exit");
}

//________________________________________________________________
double AliHFMLResponse::PredictProbaML(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int masshypo)
{
    fPtBinCand = FindPtBin(cand->Pt());
    SetMapOfVariables(cand, bfield, pidHF, masshypo);
    if (fVars.empty())
    {
        AliWarning("Map of features empty!");
        return -999.;
    }

    vector<double> vecOfFeatures;
    for (auto varname : fModelVarNames)
    {
        if (fVars.find(varname) == fVars.end())
        { //variable name not found in variables map!
            AliWarning("Variable (feature) used for ML training not implemented for model application!");
            return -999.;
        }

        vecOfFeatures.push_back(fVars[varname]);
    }

    if (fPtBinCand > fModels.size()-1){
        AliWarning(Form("Model for pT bin %d not loaded!",fPtBinCand));
        return -999.;
    }

    double modelPred = fModels[fPtBinCand].Predict(&vecOfFeatures[0], vecOfFeatures.size());

    return modelPred;
}

//________________________________________________________________
double AliHFMLResponse::PredictProbaML(double pt, vector<double> variables)
{
    fPtBinCand = FindPtBin(pt);
    if (fPtBinCand > fModels.size()-1){
        AliWarning(Form("Model for pT bin %d not loaded!",fPtBinCand));
        return -999.;
    }

    double modelPred = fModels[fPtBinCand].Predict(&variables[0], variables.size());

    return modelPred;
}

//________________________________________________________________
bool AliHFMLResponse::IsSelectedML(double &prob, AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int masshypo)
{
    prob = PredictProbaML(cand, bfield, pidHF, masshypo);
    if (prob > fModelOutputCuts[fPtBinCand])
        return true;

    return false;
}

//________________________________________________________________
bool AliHFMLResponse::IsSelectedML(double &prob, double pt, vector<double> variables)
{
    prob = PredictProbaML(pt, variables);
    if (prob > fModelOutputCuts[fPtBinCand])
        return true;

    return false;
}

//_________________________________________________________________________
string AliHFMLResponse::GetFile(const string path)
{
    if (path.find("alien:") != string::npos)
    {
        size_t currPos = path.find_last_of("/") + 1;
        string modelName = path.substr(currPos);
        if (gGrid == nullptr)
        {
            TGrid::Connect("alien://");
            if (gGrid == nullptr)
            {
                AliFatal("Connection to GRID not established! Exit");
            }
        }
        string newPath = gSystem->pwd() + string("/") + modelName.data();
        const char *OldRootDir = gDirectory->GetPath();
        bool cpStatus = TFile::Cp(path.data(), newPath.data());
        gDirectory->cd(OldRootDir);
        if (!cpStatus)
        {
            AliFatal("Error in coping file from Alien! Exit");
        }
        return newPath;
    }
    else
    {
        return path;
    }
}

//_________________________________________________________________________
void AliHFMLResponse::InitModels()
{
    for (auto iMod = 0; iMod<fModelPaths.size(); iMod++)
    {
        string localpath = GetFile(fModelPaths[iMod]);
        AliExternalBDT model = AliExternalBDT();
        bool loadmodel = false;
        switch (kLibMap[fModelLibraries[iMod]])
        {
        case kXGBoost:
        {
            loadmodel = model.LoadXGBoostModel(localpath);
            break;
        }
        case kLightGBM:
        {
            loadmodel = model.LoadLightGBMModel(localpath);
            break;
        }
        case kModelLibrary:
        {
            loadmodel = model.LoadModelLibrary(localpath);
            break;
        }
        default:
        {
            loadmodel = model.LoadXGBoostModel(localpath);
            break;
        }
        }
        if (!loadmodel)
            AliFatal("Problem in loading model");
        fModels.push_back(model);
    }
}

//_________________________________________________________________________
int AliHFMLResponse::FindPtBin(double pt)
{
    int bin = TMath::BinarySearch(fPtBinsModel.size() - 1, &fPtBinsModel[0], pt);
    if (bin < 0) //underflow --> equal to min value
        bin = 0;

    return bin;
}

//________________________________________________________________
double AliHFMLResponse::ComputeMaxd0MeasMinusExp(AliAODRecoDecayHF *cand, double bfield)
{
    double dd0max = 0;
    unsigned int fNProngsCand = static_cast<unsigned int>(cand->GetNProngs());
    for (unsigned int iProng = 0; iProng < fNProngsCand; iProng++)
    {
        double d0diff, errd0diff;
        cand->Getd0MeasMinusExpProng(iProng, bfield, d0diff, errd0diff);
        double normdd0 = d0diff / errd0diff;
        if (iProng == 0 || TMath::Abs(normdd0) > TMath::Abs(dd0max))
            dd0max = normdd0;
    }
    return dd0max;
}

//________________________________________________________________
double AliHFMLResponse::CombineNsigmaTPCTOF(double nsigmaTPC, double nsigmaTOF)
{
    if (nsigmaTPC > -998. && nsigmaTOF > -998.)
        return TMath::Sqrt((nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF) / 2);
    else if (nsigmaTPC > -998. && nsigmaTOF < -998.)
        return TMath::Abs(nsigmaTPC);
    else if (nsigmaTPC < -998. && nsigmaTOF > -998.)
        return TMath::Abs(nsigmaTOF);
    else
        return -999.;
}
