AliAnalysisTaskEmcalJetValidation* AddTaskEmcalJetValidation(
  TString suffix="",
  TString jsonconfigfile="",
  Bool_t readMC=kFALSE
)
{
   return AliAnalysisTaskEmcalJetValidation::AddTask(
    suffix,
    jsonconfigfile,
    readMC
    );
}
