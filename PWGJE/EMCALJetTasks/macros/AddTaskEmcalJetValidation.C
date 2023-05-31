AliAnalysisTaskEmcalJetValidation* AddTaskEmcalJetValidation(
  TString suffix="",
  UInt_t trigger= AliVEvent::kINT7,
  TString jsonconfigfile="",
  Bool_t readMC=kFALSE
)
{
   return AliAnalysisTaskEmcalJetValidation::AddTask(
    suffix,
    trigger,
    jsonconfigfile,
    readMC
    );
}
