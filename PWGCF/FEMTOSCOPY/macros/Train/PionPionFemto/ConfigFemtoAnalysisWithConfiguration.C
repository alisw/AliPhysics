///
/// \file ConfigFemtoAnalysisWithConfiguration.C
///

#if !defined(__CINT__) && !defined(__CLING__)

#include "AliFemtoConfigObject.h"

#endif

#if __cplusplus < 201103L
// #ifndef nullptr
#define nullptr NULL
#endif

AliFemtoManager*
ConfigFemtoAnalysis(const TString& param_str="")
{
  AliFemtoConfigObject cfg = AliFemtoConfigObject::Parse(param_str);
  if (!cfg.is_map()) {
    std::cerr << "\n\n ERROR - Configuration Object must be 'map' type; found " << cfg.name_of_type() << "\n\n";
    return nullptr;
  }

  AliFemtoManager *mgr = new AliFemtoManager();

  // setup eventreader
  AliFemtoConfigObject event_reader_cfg;
  cfg.pop_and_load("eventreader", event_reader_cfg);
  mgr->SetEventReader(AliFemtoAnalysisPionPion::ConstructEventReader(event_reader_cfg));

  if (cfg.has_key("analysis")) {
    std::cout << "Analysis\n";
  }

  return mgr;
}
