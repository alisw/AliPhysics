#include <Reve/Reve.h>
#include <Reve/ReveManager.h>

#include <Getline.h>

int main(int argc, char **argv)
{
  Reve::SetupEnvironment();

  int r = Reve::ReveManager::SpawnGuiAndRun(argc, argv);
  Getlinem(kCleanUp, 0);
  return r;
}
