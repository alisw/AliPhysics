// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file AliExternalBDT.cxx
/// \author maximiliano.puccio@cern.ch, pietro.fecchio@cern.ch

#include "AliExternalBDT.h"

#include <cassert>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <TSystem.h>

namespace {
  inline bool checkFile (const std::string name) {
    FILE *file = fopen(name.c_str(), "r");
    if (file != NULL) {
      fclose(file);
      return true;
    } else {
      return false;
    }
  }
}

AliExternalBDT::AliExternalBDT(std::string name) :
  fBDTname{name},
  fModel{},
  fModelPath{""},
  fModelName{""},
  fCompiler{},
  fPredictor{},
  fOutSize{0u},
  fNumFeatures{0u}
{
  gSystem->Setenv("TREELITE_BIND_THREADS","0");
}


bool AliExternalBDT::CompileAndLoadModelLibrary() {
  std::string path = GetUniquePath();
  if (checkFile(path + "/main.so")) {
    std::cout << "Library found: " << path.data() << "/main.so . Loading it!" << std::endl;
  } else {
    std::cout << "Starting the model compilation, depending on the model size it can take a while..." << std::endl;
    system((std::string("gcc -c -O1 -fPIC ") + path + "/main.c -o " + path + "/main.o && gcc -shared " + path + \
          "/main.o -o " + path + "/main.so").data());
  }
  return LoadModelLibrary(path + "/main.so");
}

bool AliExternalBDT::CreateModelCode() {
  std::string path = GetUniquePath();
  if (checkFile(path + "/main.c")) {
    std::cout << "Code found: " << path.data() << "/main.c . \
      Remove it or unset/change the AliExternalBDT name to force its regeneration." << std::endl;
  } else {
    const int status_comp = TreeliteCompilerCreate("ast_native", &fCompiler);
    if (status_comp != 0) {
      std::cerr << "Compiler creation failed." << std::endl;
      return false;
    }
    const int status_gen = TreeliteCompilerGenerateCode(fCompiler, fModel, 1, path.data());
    if (status_gen != 0) {
      std::cerr << "Code generation failed." << std::endl;
      return false;
    }
  }
  return true;
}

std::string AliExternalBDT::GetUniquePath() {
  if (fBDTname.empty()) {
    return fModelName + std::to_string((unsigned long)this);
  } else {
    return fBDTname + "_" + fModelName;
  }
}

bool AliExternalBDT::LoadModel(const std::string &path, int type) {
  if (path.empty()) {
    std::cout << "Invalid empty model path string" << std::endl;
    return false;
  }

  // skip steps below in case path points already to precompiled model
  if(path.substr(path.find_last_of(".")+1,path.size()) == "so"){
    std::cout << "User provided a pre-compiled model library. This may not work across all systems. Use with care!" << std::endl;
    return LoadModelLibrary(path);
  }

  fModelPath = path;
  fModelName = fModelPath.substr(fModelPath.find_last_of("\\/")+1,fModelPath.size());
  int status = 0;
  switch (type) {
    case 0:
      status = TreeliteLoadXGBoostModel(fModelPath.data(), &fModel);
      break;
    case 1:
      status = TreeliteLoadLightGBMModel(fModelPath.data(), &fModel);
      break;
    default:
      std::cerr << "Invalid model type" << std::endl;
      return false;
  }
  if (status != 0) {
    std::cerr << "Model loading failed" << std::endl;
    return false;
  }
  if (!CreateModelCode()) return false;
  if (!CompileAndLoadModelLibrary()) return false;
  return true;
}

bool AliExternalBDT::LoadXGBoostModel(std::string path) {
  if (!LoadModel(path, 0)) return false;
  return true;
}

bool AliExternalBDT::LoadLightGBMModel(std::string path) {
  if (!LoadModel(path, 1)) return false;
  return true;
}

bool AliExternalBDT::LoadModelLibrary(std::string path) {
  const int status = TreelitePredictorLoad(path.data(), 1, &fPredictor);

  TreelitePredictorQueryResultSizeSingleInst(fPredictor, &fOutSize);
  TreelitePredictorQueryNumFeature(fPredictor, &fNumFeatures);

  if (status != 0) {
    std::cerr << "Library loading failed" << std::endl;
    return false;
  }
  return true;
}

bool AliExternalBDT::Predict(double *features, int size, std::vector<double> &outputScores, bool useRawScore) {
  std::vector<TreelitePredictorEntry> entries(size);
  for (std::size_t iEntry = 0; iEntry < entries.size(); ++iEntry) {
    entries[iEntry].fvalue = static_cast<float>(features[iEntry]);
  }

  std::vector<float> output(fOutSize);
  int predict = TreelitePredictorPredictInst(fPredictor, entries.data(),
      static_cast<int>(useRawScore), &output[0],
      &fOutSize);
  if(predict<0)
    return false;

  for (std::size_t iEntry = 0; iEntry < fOutSize; ++iEntry) {
    outputScores.push_back(static_cast<double>(output[iEntry]));
  }

  return true;
}
