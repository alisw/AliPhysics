#include <Reve/Reve.h>
#include <Reve/RGTopFrame.h>

#include <Getline.h>

int main(int argc, char **argv)
{
  Reve::SetupEnvironment();

  int r = Reve::RGTopFrame::SpawnGuiAndRun(argc, argv);
  Getlinem(kCleanUp, 0);
  return r;
}
