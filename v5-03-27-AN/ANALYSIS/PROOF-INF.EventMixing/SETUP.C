Int_t SETUP()
{
 
   // Load the EventMixing library
   TString dypath = TString::Format("%s:%s", gSystem->WorkingDirectory(), gSystem->GetDynamicPath());
   gSystem->SetDynamicPath(dypath);
   
   if (gSystem->Load("libEventMixing")<0) return -1;

   // Set the include paths
   gROOT->ProcessLine(".include EventMixing");
   gROOT->ProcessLine(".include EventMixing/EventMixing");

   // Set our location, so that other packages can find us
   gSystem->Setenv("EventMixing_INCLUDE", "EventMixing");

   // Set our lib coordinates, so that other packages can link to us
   TString lib = TString::Format("-L%s -lEventMixing", gSystem->WorkingDirectory());
   gSystem->Setenv("EventMixng_LIBS", lib.Data());
   
   return 0;
   
}
