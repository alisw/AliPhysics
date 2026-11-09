#include "../../trains/ParUtilities.C"

// Build and load the analysis packages in the current ROOT session.
Bool_t LoadPars()
{
  const char* pkgs[] = {"STEERBase", "ESD", "AOD", "ANALYSIS",
                        "ANALYSISalice", "PWGLFforward2", 0};
  for (const char** pkg = pkgs; *pkg; ++pkg) {
    if (!ParUtilities::Load(*pkg)) return false;
  }
  return true;
}
