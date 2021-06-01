#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliAnalysisHelperJetTasks+;
#pragma link C++ class AliBasicParticle+;
#pragma link C++ class AliFigure+;
#pragma link C++ class AliCanvas+;
#pragma link C++ class AliHelperPID+;
#pragma link C++ class AliLatexTable+;
#pragma link C++ class AliNamedArrayI+;
#pragma link C++ class AliNamedString+;
#pragma link C++ class AliPWGFunc+;
#pragma link C++ class AliPWGHistoTools+;
#pragma link C++ typedef AliTHn;
#pragma link C++ typedef AliTHnD;
#pragma link C++ class AliTHnBase+;
#pragma link C++ class AliTHnT<TArrayF, Float_t>+;
#pragma link C++ class AliTHnT<TArrayD, Double_t>+;
#pragma link C++ class THistManager+;
#pragma link C++ class AliJSONReader+;
#pragma link C++ class AliJSONData+;
#pragma link C++ class AliJSONValue+;
#pragma link C++ class AliJSONInt+;
#pragma link C++ class AliJSONFloat+;
#pragma link C++ class AliJSONDouble+;
#pragma link C++ class AliJSONBool+;
#pragma link C++ class AliJSONString+;
#pragma link C++ class AliAnalysisTaskDummy+;
#pragma link C++ class AliTLorentzVector+;
#if ROOT_VERSION_CODE > ROOT_VERSION(6,4,0)
#pragma link C++ class AliMCSpectraWeights+;
#pragma link C++ class AliMCSpectraWeightsHandler+;
#pragma link C++ namespace YAML+;
#pragma link C++ class YAML::Node+;
#endif
#pragma link C++ class std::pair<std::string, std::string>+;
#pragma link C++ class std::vector<std::pair<std::string, std::string> >+;
#pragma link C++ namespace PWG+;
#pragma link C++ namespace PWG::Tools+;
#pragma link C++ class PWG::Tools::AliYAMLConfiguration+;
#pragma link C++ namespace TestTHistManager;
#pragma link C++ class TBinning+;
#pragma link C++ class TCustomBinning+;
#pragma link C++ class TLinearBinning+;
#pragma link C++ class TVariableBinning+;
#pragma link C++ class TestTHistManager::THistManagerTestSuite;
#pragma link C++ function TestTHistManager::TestRunAll();
#pragma link C++ function TestTHistManager::TestRunBuildSimple();
#pragma link C++ function TestTHistManager::TestRunBuildGrouped();
#pragma link C++ function TestTHistManager::TestRunFillSimple();
#pragma link C++ function TestTHistManager::TestRunFillGrouped();
#endif
