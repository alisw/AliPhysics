PWGJE:EMCALJetTasks::AliAnalysisTaskEmcalQoverPtShift *AddTaskEmcalQoverPtShift(const char *trigger, const char *tag) {
  TString tagstr(tag);
  Double_t sign = trigger[0] == 'M' ? -1. : 1.;
  Double_t absshift = tagstr.Substr(1, tagstr.Length()-1).Atoi() * 1e5;
  Double_t shift = sign * absshift

  std::cout << "Using Q/pt shift " << shift << std::endl;
  return PWG::EMCAL::AliAnalysisTaskEmcalQoverPtShift::AddTaskEmcalQOverPtShift()
}