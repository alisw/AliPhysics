// this macro runs analyze.cxx, which takes as input an Ascii starlight output
// file, slight.out, and creates a standard set of histograms, which are stored
// in histograms.root

#include "Analyze.cxx"


void Ana()
{
  Analyze a("slight.out", 20);
  a.doAnalysis();
}
