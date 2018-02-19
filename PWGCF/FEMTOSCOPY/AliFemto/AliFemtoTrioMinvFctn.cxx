///
/// \file AliFemtoTrioMinvFctn.cxx
/// \author Jeremi Niedziela

#include "AliFemtoTrioMinvFctn.h"

AliFemtoTrioMinvFctn::AliFemtoTrioMinvFctn(const char* name, int nBins, double min, double max,bool doMinv, bool doDalitz) :
fDoMinv(doMinv),
fDoDalitz(doDalitz)
{
  if(fDoMinv){
    fRealDistribution = new TH1D(Form("real_m_inv_%s",name),Form("real_m_inv_%s",name),nBins,min,max);
    fMixedDistribution = new TH1D(Form("mixed_m_inv_%s",name),Form("mixed_m_inv_%s",name),nBins,min,max);
    fRealDistribution->Sumw2();
    fMixedDistribution->Sumw2();
  }
  if(fDoDalitz){
    fDalitzPlot = new TH2D(Form("dalitz_%s",name),Form("dalitz_%s",name),300,0.0,3.0,300,0.0,3.0);
  }
  
}

AliFemtoTrioMinvFctn::~AliFemtoTrioMinvFctn()
{
  if(fRealDistribution) delete fRealDistribution;
  if(fMixedDistribution) delete fMixedDistribution;
  if(fDalitzPlot) delete fDalitzPlot;
}

void AliFemtoTrioMinvFctn::AddRealTrio(AliFemtoTrio *trio)
{
  // check if particle trio passes all cuts
  if (fTrioCut && !fTrioCut->Pass(trio)){
    return;
  }
  if(fDoMinv)   fRealDistribution->Fill(trio->MInv());
  if(fDoDalitz) fDalitzPlot->Fill(pow(trio->MInv12(),2),pow(trio->MInv23(),2));
}
void AliFemtoTrioMinvFctn::AddMixedTrio(AliFemtoTrio *trio)
{
  // check if particle trio passes all cuts
  if (fTrioCut && !fTrioCut->Pass(trio)){
    return;
  }
  if(fDoMinv) fMixedDistribution->Fill(trio->MInv());
}

void AliFemtoTrioMinvFctn::Write()
{
  // Write out neccessary objects
  if(fDoMinv){
    fRealDistribution->Write();
    fMixedDistribution->Write();
  }
  if(fDoDalitz){
    fDalitzPlot->Write();
  }
}

TList* AliFemtoTrioMinvFctn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *outputList = new TList();
  if(fDoMinv){
    outputList->Add(fRealDistribution);
    outputList->Add(fMixedDistribution);
  }
  if(fDoDalitz){
    outputList->Add(fDalitzPlot);
  }
  return outputList;
}

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoTrioMinvFctn);
  /// \endcond
#endif
