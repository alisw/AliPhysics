///
/// \file AliFemtoTrioDEtaDhiFctn.cxx
/// \author Lukasz Graczykowski

#include "AliFemtoTrioDEtaDPhiFctn.h"

#define PIH 1.57079632679489656
#define PIT 6.28318530717958623

AliFemtoTrioDEtaDPhiFctn::AliFemtoTrioDEtaDPhiFctn(const char* title, const int& aPhiBins, const int& aEtaBins) :
  fRealDPhi(0),
  fMixedDPhi(0),
  fRealDEta(0),
  fMixedDEta(0),
  fphiL(0),
  fphiT(0),
  fEtaBins(aEtaBins),
  fPhiBins(aPhiBins),
  fTitle(title)
{
  fphiL = (-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;
  fphiT = 2*TMath::Pi()+(-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;
  
  // set up numerator DPhi
  char tTitNumPhi[101] = "NumTrioDPhi";
  strncat(tTitNumPhi,title, 100);
  fRealDPhi = new TH2D(tTitNumPhi,title,aPhiBins,fphiL,fphiT,aPhiBins,fphiL,fphiT);
  // set up denominator DPhi
  char tTitDenPhi[101] = "DenTrioDPhi";
  strncat(tTitDenPhi,title, 100);
  fMixedDPhi = new TH2D(tTitDenPhi,title,aPhiBins,fphiL,fphiT,aPhiBins,fphiL,fphiT);

  
  // set up numerator DEta
  char tTitNumEta[101] = "NumTrioDEta";
  strncat(tTitNumEta,title, 100);
  fRealDEta = new TH2D(tTitNumEta,title,aEtaBins,-2.0,2.0,aEtaBins,-2.0,2.0);
  // set up denominator DEta
  char tTitDenDEta[101] = "DenTrioDEta";
  strncat(tTitDenDEta,title, 100);
  fMixedDEta = new TH2D(tTitDenDEta,title,aEtaBins,-2.0,2.0,aEtaBins,-2.0,2.0);

  // to enable error bar calculation...
  fRealDPhi->Sumw2();
  fMixedDPhi->Sumw2();
  fRealDEta->Sumw2();
  fMixedDEta->Sumw2();
  
}

AliFemtoTrioDEtaDPhiFctn::AliFemtoTrioDEtaDPhiFctn(const AliFemtoTrioDEtaDPhiFctn& aTrioFctn) :
  fRealDPhi(aTrioFctn.fRealDPhi),
  fMixedDPhi(aTrioFctn.fMixedDPhi),
  fRealDEta(aTrioFctn.fRealDEta),
  fMixedDEta(aTrioFctn.fMixedDEta),
  fphiL(aTrioFctn.fphiL),
  fphiT(aTrioFctn.fphiT),
  fEtaBins(aTrioFctn.fEtaBins),
  fPhiBins(aTrioFctn.fPhiBins),
  fTitle(aTrioFctn.fTitle)
{
  // copy constructor
  if (aTrioFctn.fRealDPhi)
    fRealDPhi = new TH2D(*aTrioFctn.fRealDPhi);
  else
    fRealDPhi = 0;
  if (aTrioFctn.fMixedDPhi)
    fMixedDPhi = new TH2D(*aTrioFctn.fMixedDPhi);
  else
    fMixedDPhi = 0;
  if (aTrioFctn.fRealDEta)
    fRealDEta = new TH2D(*aTrioFctn.fRealDEta);
  else
    fRealDEta = 0;
  if (aTrioFctn.fMixedDEta)
    fMixedDPhi = new TH2D(*aTrioFctn.fMixedDEta);
  else
    fMixedDEta = 0;
}


AliFemtoTrioDEtaDPhiFctn::~AliFemtoTrioDEtaDPhiFctn()
{
  if(fRealDPhi) delete fRealDPhi;
  if(fMixedDPhi) delete fMixedDPhi;
  if(fRealDEta) delete fRealDEta;
  if(fMixedDEta) delete fMixedDEta;
  
}


AliFemtoTrioDEtaDPhiFctn& AliFemtoTrioDEtaDPhiFctn::operator=(const AliFemtoTrioDEtaDPhiFctn& aTrioFctn)
{
  // assignment operator
  if (this == &aTrioFctn)
    return *this;

  fEtaBins = aTrioFctn.fEtaBins;
  fPhiBins = aTrioFctn.fPhiBins;

  fTitle = aTrioFctn.fTitle;

  if (aTrioFctn.fRealDPhi)
    fRealDPhi = new TH2D(*aTrioFctn.fRealDPhi);
  else
    fRealDPhi = 0;
  if (aTrioFctn.fMixedDPhi)
    fMixedDPhi = new TH2D(*aTrioFctn.fMixedDPhi);
  else
    fMixedDPhi = 0;
  if (aTrioFctn.fRealDEta)
    fRealDEta = new TH2D(*aTrioFctn.fRealDEta);
  else
    fRealDEta = 0;
  if (aTrioFctn.fMixedDEta)
    fMixedDPhi = new TH2D(*aTrioFctn.fMixedDEta);
  else
    fMixedDEta = 0;



  return *this;
}

void AliFemtoTrioDEtaDPhiFctn::AddRealTrio(AliFemtoTrio *trio)
{
  // check if particle trio passes all cuts
  if (fTrioCut && !fTrioCut->Pass(trio)){
    return;
  }


  
  double phi1 = trio->Track1()->FourMomentum().Phi();
  double phi2 = trio->Track2()->FourMomentum().Phi();
  double phi3 = trio->Track3()->FourMomentum().Phi();
  
  double eta1 = trio->Track1()->FourMomentum().PseudoRapidity();
  double eta2 = trio->Track2()->FourMomentum().PseudoRapidity();
  double eta3 = trio->Track3()->FourMomentum().PseudoRapidity();
  
  
  double dphi12 = phi1 - phi2;
  while (dphi12<fphiL) dphi12+=PIT;
  while (dphi12>fphiT) dphi12-=PIT;
  
  double dphi13 = phi1 - phi3;
  while (dphi13<fphiL) dphi13+=PIT;
  while (dphi13>fphiT) dphi13-=PIT;

  //double dphi23 = phi2 - phi3;
  //while (dphi23<fphiL) dphi23+=PIT;
  //while (dphi23>fphiT) dphi23-=PIT;

  double deta12 = eta1 - eta2;
  double deta13 = eta1 - eta3;
  //double deta23 = eta2 - eta3;

  
  fRealDPhi->Fill(dphi12,dphi13);
  fRealDEta->Fill(deta12,deta13);
  

}
void AliFemtoTrioDEtaDPhiFctn::AddMixedTrio(AliFemtoTrio *trio)
{
  // check if particle trio passes all cuts
  if (fTrioCut && !fTrioCut->Pass(trio)){
    return;
  }

  double phi1 = trio->Track1()->FourMomentum().Phi();
  double phi2 = trio->Track2()->FourMomentum().Phi();
  double phi3 = trio->Track3()->FourMomentum().Phi();
  
  double eta1 = trio->Track1()->FourMomentum().PseudoRapidity();
  double eta2 = trio->Track2()->FourMomentum().PseudoRapidity();
  double eta3 = trio->Track3()->FourMomentum().PseudoRapidity();
  
  double dphi12 = phi1 - phi2;
  while (dphi12<fphiL) dphi12+=PIT;
  while (dphi12>fphiT) dphi12-=PIT;
  
  double dphi13 = phi1 - phi3;
  while (dphi13<fphiL) dphi13+=PIT;
  while (dphi13>fphiT) dphi13-=PIT;
  
  //double dphi23 = phi2 - phi3;
  //while (dphi23<fphiL) dphi23+=PIT;
  //while (dphi23>fphiT) dphi23-=PIT;

  double deta12 = eta1 - eta2;
  double deta13 = eta1 - eta3;
  //double deta23 = eta2 - eta3;

  
  fMixedDPhi->Fill(dphi12,dphi13);
  fMixedDEta->Fill(deta12,deta13);
}

void AliFemtoTrioDEtaDPhiFctn::Write()
{
  // Write out neccessary objects
  fRealDPhi->Write();
  fMixedDPhi->Write();
  fRealDEta->Write();
  fMixedDEta->Write();

}


TList* AliFemtoTrioDEtaDPhiFctn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *outputList = new TList();

  outputList->Add(fRealDPhi);
  outputList->Add(fMixedDPhi);
  outputList->Add(fRealDEta);
  outputList->Add(fMixedDEta);
   
  return outputList;
}

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoTrioDEtaDPhiFctn);
  /// \endcond
#endif
