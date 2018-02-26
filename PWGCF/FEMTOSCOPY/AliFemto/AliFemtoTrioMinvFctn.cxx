///
/// \file AliFemtoTrioMinvFctn.cxx
/// \author Jeremi Niedziela

#include "AliFemtoTrioMinvFctn.h"

AliFemtoTrioMinvFctn::AliFemtoTrioMinvFctn(const char* name, int nBins, double min, double max,bool doMinv, bool doDalitz, bool doAngles) :
fDoMinv(doMinv),
fDoDalitz(doDalitz),
fDoAngles(doAngles)
{
  if(fDoMinv){
    fRealDistribution = new TH1D(Form("real_m_inv_%s",name),Form("real_m_inv_%s",name),nBins,min,max);
    fMixedDistribution = new TH1D(Form("mixed_m_inv_%s",name),Form("mixed_m_inv_%s",name),nBins,min,max);
    fRealDistribution->Sumw2();
    fMixedDistribution->Sumw2();
  }
  if(fDoDalitz){
    fDalitzPlot12_23 = new TH2D(Form("dalitz_12_23_%s",name),Form("dalitz_12_23_%s",name),nBins,min,max,nBins,min,max);
    fDalitzPlot23_31 = new TH2D(Form("dalitz_23_31_%s",name),Form("dalitz_23_31_%s",name),nBins,min,max,nBins,min,max);
    fDalitzPlot12_31 = new TH2D(Form("dalitz_12_31_%s",name),Form("dalitz_12_31_%s",name),nBins,min,max,nBins,min,max);
  }
  if(fDoAngles){
    fAngle12 = new TH1D(Form("angles_12-3_%s",name),Form("angles_12-3_%s",name),200,-TMath::Pi(),TMath::Pi());
    fAngle23 = new TH1D(Form("angles_23-1_%s",name),Form("angles_23-1_%s",name),200,-TMath::Pi(),TMath::Pi());
    fAngle31 = new TH1D(Form("angles_31-2_%s",name),Form("angles_31-2_%s",name),200,-TMath::Pi(),TMath::Pi());
    
    fAngle12->Sumw2();
    fAngle23->Sumw2();
    fAngle31->Sumw2();
  }
}

AliFemtoTrioMinvFctn::~AliFemtoTrioMinvFctn()
{
  if(fRealDistribution) delete fRealDistribution;
  if(fMixedDistribution) delete fMixedDistribution;
  if(fDalitzPlot12_23) delete fDalitzPlot12_23;
  if(fDalitzPlot23_31) delete fDalitzPlot23_31;
  if(fDalitzPlot12_31) delete fDalitzPlot12_31;
  if(fAngle12) delete fAngle12;
  if(fAngle23) delete fAngle23;
  if(fAngle31) delete fAngle31;
}

void AliFemtoTrioMinvFctn::AddRealTrio(AliFemtoTrio *trio)
{
  // check if particle trio passes all cuts
  if (fTrioCut && !fTrioCut->Pass(trio)){
    return;
  }
  if(fDoMinv)   fRealDistribution->Fill(trio->MInv());
  if(fDoDalitz){
    fDalitzPlot12_23->Fill(pow(trio->MInv12(),2),pow(trio->MInv23(),2));
    fDalitzPlot23_31->Fill(pow(trio->MInv23(),2),pow(trio->MInv31(),2));
    fDalitzPlot12_31->Fill(pow(trio->MInv12(),2),pow(trio->MInv31(),2));
  }
  if(fDoAngles){
    double angle = trio->GetCosTheta12();
    while(angle >  TMath::Pi()) angle -= TMath::Pi();
    while(angle < -TMath::Pi()) angle += TMath::Pi();
    fAngle12->Fill(angle);
    
    angle = trio->GetCosTheta23();
    while(angle >  TMath::Pi()) angle -= TMath::Pi();
    while(angle < -TMath::Pi()) angle += TMath::Pi();
    fAngle23->Fill(angle);
    
    angle = trio->GetCosTheta31();
    while(angle >  TMath::Pi()) angle -= TMath::Pi();
    while(angle < -TMath::Pi()) angle += TMath::Pi();
    fAngle31->Fill(angle);
  }
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
    fDalitzPlot12_23->Write();
    fDalitzPlot23_31->Write();
    fDalitzPlot12_31->Write();
  }
  if(fDoAngles){
    fAngle12->Write();
    fAngle23->Write();
    fAngle31->Write();
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
    outputList->Add(fDalitzPlot12_23);
    outputList->Add(fDalitzPlot23_31);
    outputList->Add(fDalitzPlot12_31);
  }
  if(fDoAngles){
    outputList->Add(fAngle12);
    outputList->Add(fAngle23);
    outputList->Add(fAngle31);
  }
  return outputList;
}

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoTrioMinvFctn);
  /// \endcond
#endif
