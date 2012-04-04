#ifndef __CINT__
#endif
Bool_t AddRsnToManagerKStar(TList *listRsn) {
   // This will use AddRsnPairs<Name>.C
   // and for cuts AddRsnDaughterCuts<CutName>.C
   // and <opt> string is passed to AddRsnDaughterCuts<CutName>.C
   // so you can control different cut settings
   // string "<Name>:mon" means that it will add monitoring histograms from cuts
   // Note : for now you have to set gRsnUseMiniPackage = 0 to have mon histograms
   //    listRsn->Add(new TNamed("<Name>:mon","<CutName>:<opt>"));

   // defautl PID - NSigma -> TPC=3 and TOF=3
   listRsn->Add(new TNamed("KStar","KStarNsigma"));

//    // quality only
//    listRsn->Add(new TNamed("KStar","KStarNsigma:qualityonly"));

//    // this is TPC only
//    listRsn->Add(new TNamed("KStar","KStarNsigma:PiTPCnsig10_KTPCnsig20"));
//
//    // this is TOF only
//    listRsn->Add(new TNamed("KStar","KStarNsigma:PiTOFnsig20_KTOFnsig20"));
//
//    // this is TPC only for KAON and TOF only for PI
//    listRsn->Add(new TNamed("KStar","KStarNsigma:PiTPCnsig20_KTOFnsig30"));

//    // Subhash for train
//    listRsn->Add(new TNamed("KStar","KStarNsigma:PiTPCnsig10_KTPCnsig10"));
//    listRsn->Add(new TNamed("KStar","KStarNsigma:PiTPCnsig15_KTPCnsig15"));
//    listRsn->Add(new TNamed("KStar","KStarNsigma:PiTPCnsig20_KTPCnsig20"));
//    listRsn->Add(new TNamed("KStar","KStarNsigma:PiTPCnsig25_KTPCnsig25"));
//    listRsn->Add(new TNamed("KStar","KStarNsigma:PiTPCnsig30_KTPCnsig30"));
//    // END Subhash for train

   return kTRUE;
}
