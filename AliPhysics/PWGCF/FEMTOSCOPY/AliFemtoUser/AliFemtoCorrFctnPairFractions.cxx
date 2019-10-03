////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnPairFractions - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Malgorzata Janik, majanik@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnPairFractions.h"
#include "AliFemtoModelHiddenInfo.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoCorrFctnPairFractions);
  /// \endcond
#endif

#define PIH 1.57079632679489656
#define PIT 6.28318530717958623

//____________________________
AliFemtoCorrFctnPairFractions::AliFemtoCorrFctnPairFractions(const char* title):
AliFemtoCorrFctn(),
  fPairFractions(0),
  fPairFractionsDen(0),
  fphiL(0),
  fphiT(0),
  detadphi(0)
{
  TString  hname  = "hPairFraction"; hname+= title;
  TString  htitle = "Pair Fraction "; htitle+= title;
  fPairFractions = new TH1F(hname.Data(),htitle.Data(), 9, 0, 9);
  fPairFractions->GetXaxis()->SetBinLabel(1,"#pi#pi, MC");
  fPairFractions->GetXaxis()->SetBinLabel(2,"KK, MC");
  fPairFractions->GetXaxis()->SetBinLabel(3,"pp, MC");
  fPairFractions->GetXaxis()->SetBinLabel(4,"#pi K, MC");
  fPairFractions->GetXaxis()->SetBinLabel(5,"#pi p, MC");
  fPairFractions->GetXaxis()->SetBinLabel(6,"Kp, MC");
  fPairFractions->GetXaxis()->SetBinLabel(7,"e+, MC");
  fPairFractions->GetXaxis()->SetBinLabel(8,"#mu+, MC");
  fPairFractions->GetXaxis()->SetBinLabel(9,"Other, MC");

  hname  = "hPairFractionDen"; hname+= title;
  htitle = "Pair Fraction in Mixing "; htitle+= title;
  fPairFractionsDen = new TH1F(hname.Data(),htitle.Data(), 9, 0, 9);
  fPairFractionsDen->GetXaxis()->SetBinLabel(1,"#pi#pi, MC");
  fPairFractionsDen->GetXaxis()->SetBinLabel(2,"KK, MC");
  fPairFractionsDen->GetXaxis()->SetBinLabel(3,"pp, MC");
  fPairFractionsDen->GetXaxis()->SetBinLabel(4,"#pi K, MC");
  fPairFractionsDen->GetXaxis()->SetBinLabel(5,"#pi p, MC");
  fPairFractionsDen->GetXaxis()->SetBinLabel(6,"Kp, MC");
  fPairFractionsDen->GetXaxis()->SetBinLabel(7,"e+, MC");
  fPairFractionsDen->GetXaxis()->SetBinLabel(8,"#mu+, MC");
  fPairFractionsDen->GetXaxis()->SetBinLabel(9,"Other, MC");

  // to enable error bar calculation...
  fPairFractions->Sumw2();
  fPairFractionsDen->Sumw2();
}

//____________________________
AliFemtoCorrFctnPairFractions::AliFemtoCorrFctnPairFractions(const AliFemtoCorrFctnPairFractions& aCorrFctn):
  AliFemtoCorrFctn(),
  fPairFractions(0),
  fPairFractionsDen(0),
  fphiL(aCorrFctn.fphiL),
  fphiT(aCorrFctn.fphiT),
  detadphi(aCorrFctn.detadphi)
{
  // copy constructor
  fPairFractions = new TH1F(*aCorrFctn.fPairFractions);
  fPairFractionsDen = new TH1F(*aCorrFctn.fPairFractionsDen);

 if (detadphi) {
   for (int i = 0; i < 7; i++) {
     fPairFractionsDEtaDPhi[i] = new TH2F(*aCorrFctn.fPairFractionsDEtaDPhi[i]);
     fPairFractionsDenDEtaDPhi[i] = new TH2F(*aCorrFctn.fPairFractionsDenDEtaDPhi[i]);
   }
 }
}
//____________________________
AliFemtoCorrFctnPairFractions::~AliFemtoCorrFctnPairFractions()
{
  // destructor
  delete fPairFractions;
  delete fPairFractionsDen;

  if (detadphi) {
    for (int i = 0; i < 7; i++) {
      delete fPairFractionsDEtaDPhi[i];
      delete fPairFractionsDenDEtaDPhi[i];
    }
  }
}
//_________________________
AliFemtoCorrFctnPairFractions& AliFemtoCorrFctnPairFractions::operator=(const AliFemtoCorrFctnPairFractions& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn) {
    return *this;
  }

  *fPairFractions = *aCorrFctn.fPairFractions;
  *fPairFractionsDen = *aCorrFctn.fPairFractionsDen;

  fphiL = aCorrFctn.fphiL;
  fphiT = aCorrFctn.fphiT;

  if (detadphi && aCorrFctn.detadphi) {
    for (int i = 0; i < 7; i++) {
        *fPairFractionsDEtaDPhi[i] = *aCorrFctn.fPairFractionsDEtaDPhi[i];
        *fPairFractionsDenDEtaDPhi[i] = *aCorrFctn.fPairFractionsDenDEtaDPhi[i];
    }
  } else if (aCorrFctn.detadphi) {
    for (int i = 0; i < 7; i++) {
        fPairFractionsDEtaDPhi[i] = new TH2F(*aCorrFctn.fPairFractionsDEtaDPhi[i]);
        fPairFractionsDenDEtaDPhi[i] = new TH2F(*aCorrFctn.fPairFractionsDenDEtaDPhi[i]);
    }
  } else if (detadphi) {
    for (int i = 0; i < 7; i++) {
       delete fPairFractionsDEtaDPhi[i];
       fPairFractionsDEtaDPhi[i] = nullptr;

       delete fPairFractionsDenDEtaDPhi[i];
       fPairFractionsDenDEtaDPhi[i] = nullptr;
    }
  }

  detadphi = aCorrFctn.detadphi;

  return *this;
}
//_________________________
void AliFemtoCorrFctnPairFractions::Finish()
{
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  // mShareDenominator->Draw();
  // mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnPairFractions::Report()
{
  // create report
  AliFemtoString report = "Pair Fractions Correlation Function Report:\n";
  report += Form("Number of entries in numerator:\t%E\n",fPairFractions->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n",fPairFractions->GetEntries());
  //  stemp += mCoulombWeight->Report();
  return report;
}
//____________________________
void AliFemtoCorrFctnPairFractions::AddRealPair(AliFemtoPair* pair)
{
  // add real (effect) pair

  //Applying pair cuts
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }


  Int_t pdg1=0;
  AliFemtoModelHiddenInfo *info1 = ( AliFemtoModelHiddenInfo *) pair->Track1()->GetHiddenInfo();
  if(info1)pdg1 = info1->GetPDGPid();
  Int_t pdg2=0;
  AliFemtoModelHiddenInfo *info2 = ( AliFemtoModelHiddenInfo *) pair->Track2()->GetHiddenInfo();
  if(info2)pdg2 = info2->GetPDGPid();

  if(abs(pdg1)==211 && abs(pdg2)==211) //pi pi
      fPairFractions->Fill(0.5);
  else if(abs(pdg1)==321 && abs(pdg2)==321)// K K
      fPairFractions->Fill(1.5);
  else if(abs(pdg1)==2212 && abs(pdg2)==2212)// p p
      fPairFractions->Fill(2.5);
  else if((abs(pdg1)==211 && abs(pdg2)==321)||(abs(pdg1)==321 && abs(pdg2)==211))// pi K
      fPairFractions->Fill(3.5);
  else if((abs(pdg1)==211 && abs(pdg2)==2212)||(abs(pdg1)==2212 && abs(pdg2)==211))// pi p
      fPairFractions->Fill(4.5);
  else if((abs(pdg1)==321 && abs(pdg2)==2212)||(abs(pdg1)==2212 && abs(pdg2)==321))//K p
      fPairFractions->Fill(5.5);
  else if(abs(pdg1)==13 || abs(pdg2)==13)//one particle from the pair is electron
      fPairFractions->Fill(6.5);
 else if(abs(pdg1)==11 || abs(pdg2)==11)//one particle from the pair is muon
      fPairFractions->Fill(7.5);
  else //other
    {
      fPairFractions->Fill(8.5);
    }

  if(detadphi){
    // double phi1 = pair->Track1()->Track()->P().Phi();
    // double phi2 = pair->Track2()->Track()->P().Phi();
    // double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
    // double eta2 = pair->Track2()->Track()->P().PseudoRapidity();

    double phi1 = pair->Track1()->FourMomentum().Phi();
    double phi2 = pair->Track2()->FourMomentum().Phi();
    double eta1 = pair->Track1()->FourMomentum().PseudoRapidity();
    double eta2 = pair->Track2()->FourMomentum().PseudoRapidity();

    double dphi = phi1 - phi2;
    while (dphi<fphiL) dphi+=PIT;
    while (dphi>fphiT) dphi-=PIT;

    double deta = eta1 - eta2;

    if(abs(pdg1)==211 && abs(pdg2)==211) //pi pi
      fPairFractionsDEtaDPhi[0]->Fill(dphi,deta);
    else if((abs(pdg1)==211 && abs(pdg2)==321)||(abs(pdg1)==321 && abs(pdg2)==211))// pi K
      fPairFractionsDEtaDPhi[1]->Fill(dphi,deta);
    else if((abs(pdg1)==211 && abs(pdg2)==2212)||(abs(pdg1)==2212 && abs(pdg2)==211))// pi p
      fPairFractionsDEtaDPhi[2]->Fill(dphi,deta);
    else if(abs(pdg1)==321 && abs(pdg2)==321)// K K
      fPairFractionsDEtaDPhi[3]->Fill(dphi,deta);
    else if((abs(pdg1)==321 && abs(pdg2)==2212)||(abs(pdg1)==2212 && abs(pdg2)==321))//K p
      fPairFractionsDEtaDPhi[4]->Fill(dphi,deta);
    else if(abs(pdg1)==2212 && abs(pdg2)==2212)// p p
      fPairFractionsDEtaDPhi[5]->Fill(dphi,deta);
    else //other
	fPairFractionsDEtaDPhi[6]->Fill(dphi,deta);

  }
/*
   double px1 = pair->Track1()->Track()->P().x();
   double py1 = pair->Track1()->Track()->P().y();
   double px2 = pair->Track2()->Track()->P().x();
   double py2 = pair->Track2()->Track()->P().y();
   double pt1 = TMath::Hypot(px1, py1);
   double pt2 = TMath::Hypot(px2, py2);



   double PionMass = 0.13956995;*/

}
//____________________________
void AliFemtoCorrFctnPairFractions::AddMixedPair(AliFemtoPair* pair)
{
  // add mixed (background) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  Int_t pdg1=0;
  AliFemtoModelHiddenInfo *info1 = ( AliFemtoModelHiddenInfo *) pair->Track1()->GetHiddenInfo();
  if(info1)pdg1 = info1->GetPDGPid();
  Int_t pdg2=0;
  AliFemtoModelHiddenInfo *info2 = ( AliFemtoModelHiddenInfo *) pair->Track2()->GetHiddenInfo();
  if(info2)pdg2 = info2->GetPDGPid();

  if(abs(pdg1)==211 && abs(pdg2)==211) //pi pi
    fPairFractionsDen->Fill(0.5);
  else if(abs(pdg1)==321 && abs(pdg2)==321)// K K
    fPairFractionsDen->Fill(1.5);
  else if(abs(pdg1)==2212 && abs(pdg2)==2212)// p p
    fPairFractionsDen->Fill(2.5);
  else if((abs(pdg1)==211 && abs(pdg2)==321)||(abs(pdg1)==321 && abs(pdg2)==211))// pi K
    fPairFractionsDen->Fill(3.5);
  else if((abs(pdg1)==211 && abs(pdg2)==2212)||(abs(pdg1)==2212 && abs(pdg2)==211))// pi p
    fPairFractionsDen->Fill(4.5);
  else if((abs(pdg1)==321 && abs(pdg2)==2212)||(abs(pdg1)==2212 && abs(pdg2)==321))//K p
    fPairFractionsDen->Fill(5.5);
  else if(abs(pdg1)==13 || abs(pdg2)==13)//one particle from the pair is electron
    fPairFractionsDen->Fill(6.5);
  else if(abs(pdg1)==11 || abs(pdg2)==11)//one particle from the pair is muon
    fPairFractionsDen->Fill(7.5);
  else //other
    {
      fPairFractionsDen->Fill(8.5);
    }

  if (detadphi) {
    // double phi1 = pair->Track1()->Track()->P().Phi();
    // double phi2 = pair->Track2()->Track()->P().Phi();
    // double eta1 = pair->Track1()->Track()->P().PseudoRapidity();
    // double eta2 = pair->Track2()->Track()->P().PseudoRapidity();

    double phi1 = pair->Track1()->FourMomentum().Phi();
    double phi2 = pair->Track2()->FourMomentum().Phi();
    double eta1 = pair->Track1()->FourMomentum().PseudoRapidity();
    double eta2 = pair->Track2()->FourMomentum().PseudoRapidity();

    double dphi = phi1 - phi2;
    while (dphi<fphiL) dphi+=PIT;
    while (dphi>fphiT) dphi-=PIT;

    double deta = eta1 - eta2;

    if(abs(pdg1)==211 && abs(pdg2)==211) //pi pi
      fPairFractionsDenDEtaDPhi[0]->Fill(dphi,deta);
    else if((abs(pdg1)==211 && abs(pdg2)==321)||(abs(pdg1)==321 && abs(pdg2)==211))// pi K
      fPairFractionsDenDEtaDPhi[1]->Fill(dphi,deta);
    else if((abs(pdg1)==211 && abs(pdg2)==2212)||(abs(pdg1)==2212 && abs(pdg2)==211))// pi p
      fPairFractionsDenDEtaDPhi[2]->Fill(dphi,deta);
    else if(abs(pdg1)==321 && abs(pdg2)==321)// K K
      fPairFractionsDenDEtaDPhi[3]->Fill(dphi,deta);
    else if((abs(pdg1)==321 && abs(pdg2)==2212)||(abs(pdg1)==2212 && abs(pdg2)==321))//K p
      fPairFractionsDenDEtaDPhi[4]->Fill(dphi,deta);
    else if(abs(pdg1)==2212 && abs(pdg2)==2212)// p p
      fPairFractionsDenDEtaDPhi[5]->Fill(dphi,deta);
    else //other
      fPairFractionsDenDEtaDPhi[6]->Fill(dphi,deta);
  }
}


void AliFemtoCorrFctnPairFractions::WriteHistos()
{
  // Write out result histograms
  fPairFractions->Write();
  fPairFractionsDen->Write();
  if (detadphi) {
    for (int i = 0; i < 7; i++) {
      fPairFractionsDEtaDPhi[i]->Write();
      fPairFractionsDenDEtaDPhi[i]->Write();
    }
  }
}

TList* AliFemtoCorrFctnPairFractions::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fPairFractions);
  tOutputList->Add(fPairFractionsDen);

  if (detadphi) {
    for (int i = 0; i < 7; i++) {
      tOutputList->Add(fPairFractionsDenDEtaDPhi[i]);
      tOutputList->Add(fPairFractionsDEtaDPhi[i]);
    }
  }
  return tOutputList;
}


void AliFemtoCorrFctnPairFractions::SetDoDEtaDPhiMaps(bool dodedp)
{
  detadphi = dodedp;
  if (detadphi) {
    double aPhiBins = 5;
    fphiL = (-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;
    fphiT = 2*TMath::Pi()+(-(int)(aPhiBins/4)+0.5)*2.*TMath::Pi()/aPhiBins;
    const char *title = fPairFractions->GetName();
    TString particle;

    for (int i = 0; i < 7; i++) {
      if(i==0) particle = "PiPi";
      else if(i==1) particle = "PiK";
      else if(i==2) particle = "PiP";
      else if(i==3) particle = "KK";
      else if(i==4) particle = "KP";
      else if(i==5) particle = "PP";
      else if(i==6) particle = "Other";



      TString hname  = "NumDEtaDPhi"; hname+= title; hname+= particle;
      TString htitle = " 2D DEtaDPhi Num "; htitle+= title; htitle+= particle;
      fPairFractionsDEtaDPhi[i] = new TH2F(hname.Data(),htitle.Data(),5,fphiL,fphiT,5,-2.0,2.0);

      hname  = "DenDEtaDPhi"; hname+= title;hname+= particle;
      htitle = " 2D DEtaDPhi Den"; htitle+= title;htitle+= particle;
      fPairFractionsDenDEtaDPhi[i] = new TH2F(hname.Data(),htitle.Data(),5,fphiL,fphiT,5,-2.0,2.0);

      fPairFractionsDEtaDPhi[i]->Sumw2();
      fPairFractionsDenDEtaDPhi[i]->Sumw2();
    }
  }
}
