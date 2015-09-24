///
/// \file ConfigFemtoAnalysis.C
/// \brief The configuration macro which sets up the pion-lambda analysis
/// \author Andrew Kubera, Ohio State University, andrew.kubera@cern.ch
///

#if !defined(__CINT__) || defined(__MAKECINT_)

#include "AliFemtoAnalysisPionLambda.h"

#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoEventReaderAODChain.h"

#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoCutMonitorV0.h"

#endif

typedef struct macro_params_s MacroParams;

struct macro_params_s {
  int pion_charge;
  float pion_pt_range[2];
  float pion_eta_range[2];
  float pion_DCA_range[2];
  float pion_nsigma[2];

};

MacroParams parse_arguments(const char*);

AliFemtoManager*
ConfigFemtoAnalysis(const TString& param_str = "")
{
  const double PionMass = 0.13956995,
             LambdaMass = 1.115683;

  // MacroParams params = parse_arguments(param_str);

  AliFemtoManager *manager = new AliFemtoManager();

  AliFemtoEventReaderAOD *rdr = new AliFemtoEventReaderAODChain();
  rdr->SetFilterBit(7);
  rdr->SetEPVZERO(kTRUE);
  rdr->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
  rdr->SetCentralityFlattening(kFALSE);
  rdr->SetReadV0(1);
  rdr->SetPrimaryVertexCorrectionTPCPoints(kTRUE);
  manager->SetEventReader(rdr);

  AliFemtoAnalysisPionLambda *analysis = new AliFemtoAnalysisPionLambda("PionLambda");

  analysis->EventCut()->AddCutMonitor(new AliFemtoCutMonitorEventMult("EVPass"),
                                         new AliFemtoCutMonitorEventMult("EVFail"));

  analysis->GetPionCut()->AddCutMonitor(new AliFemtoCutMonitorParticleYPt("PionPass", PionMass),
                                        new AliFemtoCutMonitorParticleYPt("PionFail", PionMass));

  analysis->GetLambdaCut()->AddCutMonitor(new AliFemtoCutMonitorV0("LambdaPass"),
                                          new AliFemtoCutMonitorV0("LambdaFail"));

  manager->AddAnalysis(analysis);

  return manager;
}

MacroParams
default_arguments()
{
  MacroParams params;
  params.pion_charge = -1;
  params.pion_pt_range[0] = 0.2;
  params.pion_pt_range[1] = 2.0;
  params.pion_eta_range[0] = -0.8;
  params.pion_eta_range[1] =  0.8;
  params.pion_DCA_range[0] =  0.5;
  params.pion_DCA_range[1] =  4.0;
  params.pion_nsigma[0] = -0.2;
  params.pion_nsigma[1] = 0.2;

  return params;
}

void
process_argument(const TString& key, const TString& val, MacroParams& params)
{
  if (key == 'pion_pt_range') {
    params.pion_pt_range[0] = val.Atoi();
    params.pion_pt_range[1] = val.Atoi();
  }

}

MacroParams
parse_arguments(const TString& in)
{
  // Get default parameters
  MacroParams params = default_arguments();

  TObjArray* tokens = in.Tokenize(',');
  TIter next_param(tokens);
  TObject *param_obj = NULL;
  while (param_obj = next_param()) {
    const TString& key_val = ((TObjString*)param_obj)->String();
    const size_t equal_position = key_val.First('='),
                      val_start = equal_position + 1,
                       val_stop = key_val.Length() - val_start;

    const TString key = key_val(0, equal_position),
                  val = key_val(val_start, val_stop);

    process_argument(key, val, params);
  }

  return params;
}
