PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalQoverPtShift *AddTaskEmcalQoverPtShift(const char *trigger, const char *tag) {
  TString tagstr(tag);
  TString valstr = tagstr(1, tagstr.Length()-1).Data();
  Double_t sign = tag[0] == 'M' ? -1. : 1.;
  Double_t absshift = double(valstr.Atoi()) / 1e5;
  Double_t shift = sign * absshift;

  std::cout << "Using Q/pt shift " << std::scientific <<  shift << std::dec << " from tag " << tag << std::endl;
  return PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalQoverPtShift::AddTaskQOverPtShift(trigger, shift);
}
