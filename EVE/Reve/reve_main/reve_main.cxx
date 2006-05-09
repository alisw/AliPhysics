#include <Reve/Reve.h>
#include <Reve/RGTopFrame.h>

#include <TInterpreter.h>
#include <Getline.h>

int main(int argc, char **argv)
{
  Reve::SetupEnvironment();

  gROOT->SetMacroPath(Form(".:%s/macros", gSystem->Getenv("REVESYS")));
  gInterpreter->AddIncludePath(Form("%s/macros", gSystem->Getenv("REVESYS")));

  int r = Reve::RGTopFrame::SpawnGuiAndRun(argc, argv);
  Getlinem(kCleanUp, 0);
  return r;
}
