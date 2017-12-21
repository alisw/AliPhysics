///
/// \file AliFemtoTrioMinvFctn.cxx
/// \author Jeremi Niedziela

#include "AliFemtoTrioMinvFctn.h"

AliFemtoTrioMinvFctn::AliFemtoTrioMinvFctn(const char* name, int nBins, double min, double max)
{
  fRealDistribution = new TH1D(Form("real_m_inv_%s",name),Form("real_m_inv_%s",name),nBins,min,max);
  fMixedDistribution = new TH1D(Form("mixed_m_inv_%s",name),Form("mixed_m_inv_%s",name),nBins,min,max);
  
  fRealDistribution->Sumw2();
  fMixedDistribution->Sumw2();
}

AliFemtoTrioMinvFctn::~AliFemtoTrioMinvFctn()
{
  delete fRealDistribution;
  delete fMixedDistribution;
}

void AliFemtoTrioMinvFctn::AddRealTrio(AliFemtoTrio *trio)
{
  // check if particle trio passes all cuts
  if (fTrioCut && !fTrioCut->Pass(trio)){
    return;
  }
  fRealDistribution->Fill(trio->MInv());
}
void AliFemtoTrioMinvFctn::AddMixedTrio(AliFemtoTrio *trio)
{
  // check if particle trio passes all cuts
  if (fTrioCut && !fTrioCut->Pass(trio)){
    return;
  }
  fMixedDistribution->Fill(trio->MInv());
}

void AliFemtoTrioMinvFctn::Write()
{
  // Write out neccessary objects
  fRealDistribution->Write();
  fMixedDistribution->Write();
}

TList* AliFemtoTrioMinvFctn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *outputList = new TList();
  
  outputList->Add(fRealDistribution);
  outputList->Add(fMixedDistribution);
  
  return outputList;
}

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoTrioMinvFctn);
  /// \endcond
#endif
