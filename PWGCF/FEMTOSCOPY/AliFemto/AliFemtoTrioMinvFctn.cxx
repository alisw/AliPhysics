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
    fAngle12 = new TH1D(Form("angles_12_%s",name),Form("angles_12_%s",name),200,-TMath::Pi(),TMath::Pi());
    fAngle23 = new TH1D(Form("angles_23_%s",name),Form("angles_23_%s",name),200,-TMath::Pi(),TMath::Pi());
    fAngle31 = new TH1D(Form("angles_31_%s",name),Form("angles_31_%s",name),200,-TMath::Pi(),TMath::Pi());
    
    fAngle1 = new TH1D(Form("angles_1_%s",name),Form("angles_1_%s",name),200,-TMath::Pi(),TMath::Pi());
    fAngle2 = new TH1D(Form("angles_2_%s",name),Form("angles_2_%s",name),200,-TMath::Pi(),TMath::Pi());
    fAngle3 = new TH1D(Form("angles_3_%s",name),Form("angles_3_%s",name),200,-TMath::Pi(),TMath::Pi());
    
    fCosAngle12 = new TH1D(Form("cos_angles_12-3_%s",name),Form("cos_angles_12-3_%s",name),200,-1,1);
    fCosAngle23 = new TH1D(Form("cos_angles_23-1_%s",name),Form("cos_angles_23-1_%s",name),200,-1,1);
    fCosAngle31 = new TH1D(Form("cos_angles_31-2_%s",name),Form("cos_angles_31-2_%s",name),200,-1,1);
    
    fCosAngle1 = new TH1D(Form("cos_angles_1_%s",name),Form("cos_angles_1_%s",name),200,-1,1);
    fCosAngle2 = new TH1D(Form("cos_angles_2_%s",name),Form("cos_angles_2_%s",name),200,-1,1);
    fCosAngle3 = new TH1D(Form("cos_angles_3_%s",name),Form("cos_angles_3_%s",name),200,-1,1);
    
    fAngle12->Sumw2();
    fAngle23->Sumw2();
    fAngle31->Sumw2();
    
    fAngle1->Sumw2();
    fAngle2->Sumw2();
    fAngle3->Sumw2();
    
    fCosAngle12->Sumw2();
    fCosAngle23->Sumw2();
    fCosAngle31->Sumw2();
    
    fCosAngle1->Sumw2();
    fCosAngle2->Sumw2();
    fCosAngle3->Sumw2();
  }
}

AliFemtoTrioMinvFctn::AliFemtoTrioMinvFctn(const AliFemtoTrioMinvFctn& aTrioFctn) :
  fDoMinv(aTrioFctn.fDoMinv),
  fDoDalitz(aTrioFctn.fDoDalitz),
  fDoAngles(aTrioFctn.fDoAngles)
{
  // copy constructor
  if(aTrioFctn.fRealDistribution)
    fRealDistribution = new TH1D(*aTrioFctn.fRealDistribution);
  else
    fRealDistribution = 0;
  if(aTrioFctn.fMixedDistribution)
    fMixedDistribution = new TH1D(*aTrioFctn.fMixedDistribution);
  else
    fMixedDistribution = 0;   
       


  
  if(fDoDalitz){
    if(aTrioFctn.fDalitzPlot12_23)
      fDalitzPlot12_23 = new TH2D(*aTrioFctn.fDalitzPlot12_23);
    else
      fDalitzPlot12_23 = 0; 
    if(aTrioFctn.fDalitzPlot23_31)
      fDalitzPlot23_31 = new TH2D(*aTrioFctn.fDalitzPlot23_31);
    else
      fDalitzPlot23_31 = 0;
    if(aTrioFctn.fDalitzPlot12_31)
      fDalitzPlot12_31 = new TH2D(*aTrioFctn.fDalitzPlot12_31);
    else
      fDalitzPlot12_31 = 0;
  }
  if(fDoAngles){
    if(aTrioFctn.fAngle12)
      fAngle12 = new TH1D(*aTrioFctn.fAngle12);
    else
      fAngle12 = 0; 
    if(aTrioFctn.fAngle23)
      fAngle23 = new TH1D(*aTrioFctn.fAngle23);
    else
      fAngle23 = 0; 
    if(aTrioFctn.fAngle31)
      fAngle31 = new TH1D(*aTrioFctn.fAngle31);
    else
      fAngle31 = 0; 

    if(aTrioFctn.fAngle1)
      fAngle1 = new TH1D(*aTrioFctn.fAngle1);
    else
      fAngle1 = 0;
    if(aTrioFctn.fAngle2)
      fAngle2 = new TH1D(*aTrioFctn.fAngle2);
    else
      fAngle1 = 0;
    if(aTrioFctn.fAngle3)
      fAngle3 = new TH1D(*aTrioFctn.fAngle3);
    else
      fAngle3 = 0; 


    if(aTrioFctn.fCosAngle12)
      fCosAngle12 = new TH1D(*aTrioFctn.fCosAngle12);
    else
      fCosAngle12 = 0; 
    if(aTrioFctn.fCosAngle23)
      fCosAngle23 = new TH1D(*aTrioFctn.fCosAngle23);
    else
      fCosAngle23 = 0; 
    if(aTrioFctn.fCosAngle31)
      fCosAngle31 = new TH1D(*aTrioFctn.fCosAngle31);
    else
      fCosAngle31 = 0; 

    if(aTrioFctn.fCosAngle1)
      fCosAngle1 = new TH1D(*aTrioFctn.fCosAngle1);
    else
      fCosAngle1 = 0;
    if(aTrioFctn.fCosAngle2)
      fCosAngle2 = new TH1D(*aTrioFctn.fCosAngle2);
    else
      fCosAngle1 = 0;
    if(aTrioFctn.fCosAngle3)
      fCosAngle3 = new TH1D(*aTrioFctn.fCosAngle3);
    else
      fCosAngle3 = 0; 
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
  if(fAngle1) delete fAngle1;
  if(fAngle2) delete fAngle2;
  if(fAngle3) delete fAngle3;
  if(fCosAngle12) delete fCosAngle12;
  if(fCosAngle23) delete fCosAngle23;
  if(fCosAngle31) delete fCosAngle31;
  if(fCosAngle1) delete fCosAngle1;
  if(fCosAngle2) delete fCosAngle2;
  if(fCosAngle3) delete fCosAngle3;
}


AliFemtoTrioMinvFctn& AliFemtoTrioMinvFctn::operator=(const AliFemtoTrioMinvFctn& aTrioFctn)
{
  // assignment operator
  if (this == &aTrioFctn)
    return *this;

  fDoMinv = aTrioFctn.fDoMinv;
  fDoDalitz = aTrioFctn.fDoDalitz;
  fDoAngles = aTrioFctn.fDoAngles;


  if(aTrioFctn.fRealDistribution)
    fRealDistribution = new TH1D(*aTrioFctn.fRealDistribution);
  else
    fRealDistribution = 0;
  if(aTrioFctn.fMixedDistribution)
    fMixedDistribution = new TH1D(*aTrioFctn.fMixedDistribution);
  else
    fMixedDistribution = 0;   
       


  
  if(fDoDalitz){
    if(aTrioFctn.fDalitzPlot12_23)
      fDalitzPlot12_23 = new TH2D(*aTrioFctn.fDalitzPlot12_23);
    else
      fDalitzPlot12_23 = 0; 
    if(aTrioFctn.fDalitzPlot23_31)
      fDalitzPlot23_31 = new TH2D(*aTrioFctn.fDalitzPlot23_31);
    else
      fDalitzPlot23_31 = 0;
    if(aTrioFctn.fDalitzPlot12_31)
      fDalitzPlot12_31 = new TH2D(*aTrioFctn.fDalitzPlot12_31);
    else
      fDalitzPlot12_31 = 0;
  }
  if(fDoAngles){
    if(aTrioFctn.fAngle12)
      fAngle12 = new TH1D(*aTrioFctn.fAngle12);
    else
      fAngle12 = 0; 
    if(aTrioFctn.fAngle23)
      fAngle23 = new TH1D(*aTrioFctn.fAngle23);
    else
      fAngle23 = 0; 
    if(aTrioFctn.fAngle31)
      fAngle31 = new TH1D(*aTrioFctn.fAngle31);
    else
      fAngle31 = 0; 

    if(aTrioFctn.fAngle1)
      fAngle1 = new TH1D(*aTrioFctn.fAngle1);
    else
      fAngle1 = 0;
    if(aTrioFctn.fAngle2)
      fAngle2 = new TH1D(*aTrioFctn.fAngle2);
    else
      fAngle1 = 0;
    if(aTrioFctn.fAngle3)
      fAngle3 = new TH1D(*aTrioFctn.fAngle3);
    else
      fAngle3 = 0; 


    if(aTrioFctn.fCosAngle12)
      fCosAngle12 = new TH1D(*aTrioFctn.fCosAngle12);
    else
      fCosAngle12 = 0; 
    if(aTrioFctn.fCosAngle23)
      fCosAngle23 = new TH1D(*aTrioFctn.fCosAngle23);
    else
      fCosAngle23 = 0; 
    if(aTrioFctn.fCosAngle31)
      fCosAngle31 = new TH1D(*aTrioFctn.fCosAngle31);
    else
      fCosAngle31 = 0; 

    if(aTrioFctn.fCosAngle1)
      fCosAngle1 = new TH1D(*aTrioFctn.fCosAngle1);
    else
      fCosAngle1 = 0;
    if(aTrioFctn.fCosAngle2)
      fCosAngle2 = new TH1D(*aTrioFctn.fCosAngle2);
    else
      fCosAngle1 = 0;
    if(aTrioFctn.fCosAngle3)
      fCosAngle3 = new TH1D(*aTrioFctn.fCosAngle3);
    else
      fCosAngle3 = 0; 
  
  }

  return *this;
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
    double angle = trio->GetTheta12();
    fAngle12->Fill(angle);
    fCosAngle12->Fill(cos(angle));
    angle = trio->GetTheta3();
    fAngle3->Fill(angle);
    fCosAngle3->Fill(cos(angle));
    
    angle = trio->GetTheta23();
    fAngle23->Fill(angle);
    fCosAngle23->Fill(cos(angle));
    angle = trio->GetTheta1();
    fAngle1->Fill(angle);
    fCosAngle1->Fill(cos(angle));
    
    angle = trio->GetTheta31();
    fAngle31->Fill(angle);
    fCosAngle31->Fill(cos(angle));
    angle = trio->GetTheta2();
    fAngle2->Fill(angle);
    fCosAngle2->Fill(cos(angle));
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
    fAngle1->Write();
    fAngle2->Write();
    fAngle3->Write();
    fCosAngle12->Write();
    fCosAngle23->Write();
    fCosAngle31->Write();
    fCosAngle1->Write();
    fCosAngle2->Write();
    fCosAngle3->Write();
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
    outputList->Add(fAngle1);
    outputList->Add(fAngle2);
    outputList->Add(fAngle3);
    outputList->Add(fCosAngle12);
    outputList->Add(fCosAngle23);
    outputList->Add(fCosAngle31);
    outputList->Add(fCosAngle1);
    outputList->Add(fCosAngle2);
    outputList->Add(fCosAngle3);
  }
  return outputList;
}

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoTrioMinvFctn);
  /// \endcond
#endif
