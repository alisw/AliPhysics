#ifndef __CINT__
#endif
Bool_t AddRsnToManagerLambda(TList *listRsn) {
   // This will use AddRsnPairs<Name>.C
   // and for cuts AddRsnDaughterCuts<CutName>.C
   // and <opt> string is passed to AddRsnDaughterCuts<CutName>.C
   // so you can control different cut settings
   // string "<Name>:mon" means that it will add monitoring histograms from cuts
   // Note : for now you have to set gRsnUseMiniPackage = 0 to have mon histograms
   //    listRsn->Add(new TNamed("<Name>:mon","<CutName>:<opt>"));

   // defautl PID - NSigma -> TPC=3 and TOF=3
   listRsn->Add(new TNamed("Lambda","LambdaNsigma"));

// //    quality only
//    listRsn->Add(new TNamed("Lambda","LambdaNsigma:qualityonly"));
// 
// //    this is TPC only
//    listRsn->Add(new TNamed("Lambda","LambdaNsigma:PTPCnsig10_KTPCnsig20"));
// 
// //    this is TOF only
//    listRsn->Add(new TNamed("Lambda","LambdaNsigma:PTOFnsig20_KTOFnsig20"));
// 
// //    this is TPC only for KAON and TOF only for PROTON
//    listRsn->Add(new TNamed("Lambda","LambdaNsigma:PTPCnsig20_KTOFnsig30"));

   return kTRUE;
}
