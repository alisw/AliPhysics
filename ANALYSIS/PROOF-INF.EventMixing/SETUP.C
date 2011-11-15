void SETUP()
{
  // Load the ANALYSIS library

   gSystem->SetDynamicPath(Form("%s:%s", gSystem->pwd(), gSystem->GetDynamicPath()));
   gSystem->Load("libEventMixing");

   // Set the include paths
   gROOT->ProcessLine(".include EventMixing");
   gROOT->ProcessLine(".include EventMixing/EventMixing");

   // Set our location, so that other packages can find us
   gSystem->Setenv("EventMixing_INCLUDE", "EventMixing");
}
