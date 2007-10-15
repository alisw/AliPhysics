Int_t convertFiles(Int_t runNb, Int_t nFiles)
{
   gROOT->ProcessLine(".L ../digits/AliHLTPHOSAltroConfig.cxx++");
   gROOT->ProcessLine(".L ../digits/AliHLTPHOSDebugRawDigit.cxx++");
   gROOT->ProcessLine(".L ../digits/AliHLTPHOSDigit.cxx++");
   gROOT->ProcessLine(".L digitConverter.C++");
 
   for(Int_t i = 0; i < nFiles; i++)
     {
       digitConverter(runNb, i);
     }
   return 0;
}

