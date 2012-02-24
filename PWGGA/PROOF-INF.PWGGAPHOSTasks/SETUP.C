// $Id$

void SETUP()
{
   TString ocwd = gSystem->WorkingDirectory();
   gSystem->Load(ocwd+"/libPWGGAPHOSTasks.so");

   // Set the Include paths
   gSystem->SetIncludePath("-I$ROOTSYS/include -IPWGGAPHOSTasks");
   gROOT->ProcessLine(".include PWGGAPHOSTasks/PHOSTasks");
   gROOT->ProcessLine(".include PWGGAPHOSTasks/PHOSTasks/CaloCellQA");
   gROOT->ProcessLine(".include PWGGAPHOSTasks/PHOSTasks/omega3pi");
   gROOT->ProcessLine(".include PWGGAPHOSTasks/PHOSTasks/PHOS_embedding");
   gROOT->ProcessLine(".include PWGGAPHOSTasks/PHOSTasks/PHOS_PbPbQA");
   gROOT->ProcessLine(".include PWGGAPHOSTasks/PHOSTasks/PHOS_pp_pi0");
   gROOT->ProcessLine(".include PWGGAPHOSTasks/PHOSTasks/PHOS_TriggerQA");

   // Set our location, so that other packages can find us
   gSystem->Setenv("PWGGAPHOSTasks_INCLUDE", "PWGGAPHOSTasks/PHOSTasks");
}
