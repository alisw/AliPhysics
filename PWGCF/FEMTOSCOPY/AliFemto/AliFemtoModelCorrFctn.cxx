////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelCorrFctn - the base class for correlation function which    ///
/// uses the model framework and weight generation                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoModelCorrFctn, 1);
  /// \endcond
#endif

#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include <TH1D.h>
#include <TH2D.h>

//_______________________
AliFemtoModelCorrFctn::AliFemtoModelCorrFctn():
AliFemtoCorrFctn(),
  fManager(0),
  fNumeratorTrue(0),
  fNumeratorFake(0),
  fDenominator(0),
  fNumeratorTrueIdeal(0),
  fNumeratorFakeIdeal(0),
  fDenominatorIdeal(0),
  fQgenQrec(0),
  fKaonPDG(kFALSE)
{
  // Default constructor
  fNumeratorTrue = new TH1D("ModelNumTrue","ModelNumTrue",50,0.0,0.5);
  fNumeratorFake = new TH1D("ModelNumFake","ModelNumFake",50,0.0,0.5);
  fDenominator = new TH1D("ModelDen","ModelDen",50,0.0,0.5);

  fNumeratorTrueIdeal = new TH1D("ModelNumTrueIdeal","ModelNumTrueIdeal",50,0.0,0.5);
  fNumeratorFakeIdeal = new TH1D("ModelNumFakeIdeal","ModelNumFakeIdeal",50,0.0,0.5);
  fDenominatorIdeal = new TH1D("ModelDenIdeal","ModelDenIdeal",50,0.0,0.5);

  fQgenQrec = new TH2D("QgenQrec","QgenQrec",50,0.0,0.5,50,0.0,0.5);

  fNumeratorTrue->Sumw2();
  fNumeratorFake->Sumw2();
  fDenominator->Sumw2();

  fNumeratorTrueIdeal->Sumw2();
  fNumeratorFakeIdeal->Sumw2();
  fDenominatorIdeal->Sumw2();

  fQgenQrec->Sumw2();

}
//_______________________
AliFemtoModelCorrFctn::AliFemtoModelCorrFctn(const char *title, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi):
  AliFemtoCorrFctn(),
  fManager(0),
  fNumeratorTrue(0),
  fNumeratorFake(0),
  fDenominator(0),
  fNumeratorTrueIdeal(0),
  fNumeratorFakeIdeal(0),
  fDenominatorIdeal(0),
  fQgenQrec(0),
  fKaonPDG(kFALSE)
{
  // Normal constructor
  char buf[100];
  snprintf(buf , 100,  "NumTrue%s", title);
  fNumeratorTrue = new TH1D(buf,buf,aNbins,aQinvLo,aQinvHi);
  snprintf(buf , 100,  "NumFake%s", title);
  fNumeratorFake = new TH1D(buf,buf,aNbins,aQinvLo,aQinvHi);
  snprintf(buf , 100,  "Den%s", title);
  fDenominator = new TH1D(buf,buf,aNbins,aQinvLo,aQinvHi);

  snprintf(buf , 100,  "NumTrueIdeal%s", title);
  fNumeratorTrueIdeal = new TH1D(buf,buf,aNbins,aQinvLo,aQinvHi);
  snprintf(buf , 100,  "NumFakeIdeal%s", title);
  fNumeratorFakeIdeal = new TH1D(buf,buf,aNbins,aQinvLo,aQinvHi);
  snprintf(buf , 100,  "DenIdeal%s", title);
  fDenominatorIdeal = new TH1D(buf,buf,aNbins,aQinvLo,aQinvHi);

  snprintf(buf , 100,  "QgenQrec%s", title);
  fQgenQrec = new TH2D(buf,buf,aNbins,aQinvLo,aQinvHi,aNbins,aQinvLo,aQinvHi);

  fNumeratorTrue->Sumw2();
  fNumeratorFake->Sumw2();
  fDenominator->Sumw2();

  fNumeratorTrueIdeal->Sumw2();
  fNumeratorFakeIdeal->Sumw2();
  fDenominatorIdeal->Sumw2();

  fQgenQrec->Sumw2();
}
//_______________________
AliFemtoModelCorrFctn::AliFemtoModelCorrFctn(const AliFemtoModelCorrFctn& aCorrFctn) :
  AliFemtoCorrFctn(),
  fManager(0),
  fNumeratorTrue(0),
  fNumeratorFake(0),
  fDenominator(0),
  fNumeratorTrueIdeal(0),
  fNumeratorFakeIdeal(0),
  fDenominatorIdeal(0),
  fQgenQrec(0),
  fKaonPDG(aCorrFctn.fKaonPDG)
{
  // Copy constructor
  if (aCorrFctn.fNumeratorTrue)
    fNumeratorTrue = new TH1D(*(aCorrFctn.fNumeratorTrue));
  if (aCorrFctn.fNumeratorFake)
    fNumeratorFake = new TH1D(*(aCorrFctn.fNumeratorFake));
  if (aCorrFctn.fDenominator)
    fDenominator = new TH1D(*(aCorrFctn.fDenominator));

  if (aCorrFctn.fNumeratorTrueIdeal)
    fNumeratorTrueIdeal = new TH1D(*(aCorrFctn.fNumeratorTrueIdeal));
  if (aCorrFctn.fNumeratorFakeIdeal)
    fNumeratorFakeIdeal = new TH1D(*(aCorrFctn.fNumeratorFakeIdeal));
  if (aCorrFctn.fDenominatorIdeal)
    fDenominatorIdeal = new TH1D(*(aCorrFctn.fDenominatorIdeal));

  if (aCorrFctn.fQgenQrec)
    fQgenQrec = new TH2D(*(aCorrFctn.fQgenQrec));

  fManager = aCorrFctn.fManager;
}
//_______________________
AliFemtoModelCorrFctn::~AliFemtoModelCorrFctn()
{
  // Destructor
  if (fNumeratorTrue) delete fNumeratorTrue;
  if (fNumeratorFake) delete fNumeratorFake;
  if (fDenominator) delete fDenominator;

  if (fNumeratorTrueIdeal) delete fNumeratorTrueIdeal;
  if (fNumeratorFakeIdeal) delete fNumeratorFakeIdeal;
  if (fDenominatorIdeal) delete fDenominatorIdeal;

  if (fQgenQrec) delete fQgenQrec;

}
//_______________________
AliFemtoModelCorrFctn& AliFemtoModelCorrFctn::operator=(const AliFemtoModelCorrFctn& aCorrFctn)
{
  // Assignment operator
  if (this == &aCorrFctn)
    return *this;

  delete fNumeratorTrue;
  fNumeratorTrue = aCorrFctn.fNumeratorTrue
                 ? new TH1D(*aCorrFctn.fNumeratorTrue)
                 : nullptr;

  delete fNumeratorFake;
  fNumeratorFake = (aCorrFctn.fNumeratorFake)
                 ? new TH1D(*aCorrFctn.fNumeratorFake)
                 : nullptr;

  delete fDenominator;
  fDenominator = (aCorrFctn.fDenominator)
               ? new TH1D(*aCorrFctn.fDenominator)
               : nullptr;

  delete fQgenQrec;
  fQgenQrec = (aCorrFctn.fQgenQrec)
            ? new TH2D(*aCorrFctn.fQgenQrec)
            : nullptr;

  delete fNumeratorTrueIdeal;
  fNumeratorTrueIdeal = (aCorrFctn.fNumeratorTrueIdeal)
                      ? new TH1D(*aCorrFctn.fNumeratorTrueIdeal)
                      : nullptr;

  delete fNumeratorFakeIdeal;
  fNumeratorFakeIdeal = (aCorrFctn.fNumeratorFakeIdeal)
                      ? new TH1D(*aCorrFctn.fNumeratorFakeIdeal)
                      : nullptr;

  delete fDenominatorIdeal;
  fDenominatorIdeal = (aCorrFctn.fDenominatorIdeal)
                    ? new TH1D(*aCorrFctn.fDenominatorIdeal)
                    : nullptr;

  fManager = aCorrFctn.fManager;

  fKaonPDG = aCorrFctn.fKaonPDG;

  return *this;
}
//_______________________
void AliFemtoModelCorrFctn::ConnectToManager(AliFemtoModelManager *aManager)
{
  fManager = aManager;
}

//_______________________
AliFemtoString AliFemtoModelCorrFctn::Report()
{
  // Prepare report
  AliFemtoString tStr = "AliFemtoModelCorrFctn report";

  return tStr;
}

//_______________________
void AliFemtoModelCorrFctn::AddRealPair(AliFemtoPair* aPair)
{
  if(!fKaonPDG) {
    if (fPairCut)
      if (!fPairCut->Pass(aPair)) return;
    // cout<<" AliFemtoModelCorrFcn add real pair "<<endl;
    Double_t weight = fManager->GetWeight(aPair);
    //cout<<" wight "<< weight<<endl;
    //cout<<"Qinv"<<aPair->QInv()<<endl;
    
    fNumeratorTrue->Fill(aPair->QInv(), weight);
    
    Double_t tQinvTrue = GetQinvTrue(aPair);
    
    fNumeratorTrueIdeal->Fill(tQinvTrue, weight);
    
    //cout<<"Qinv true"<<tQinvTrue<<endl;
  }
  //Special MC analysis for K selected by PDG code -->
  else {
    Double_t weight = fManager->GetWeight(aPair);
    fNumeratorTrue->Fill(aPair->QInv(), weight);
    Double_t tQinvTrue = GetQinvTrue(aPair);
    fNumeratorTrueIdeal->Fill(tQinvTrue, weight);
  }
}
//_______________________
void AliFemtoModelCorrFctn::AddMixedPair(AliFemtoPair* aPair)
{
  if(!fKaonPDG) {
    if (fPairCut)
      if (!fPairCut->Pass(aPair)) return;
    Double_t weight = fManager->GetWeight(aPair);
    fNumeratorFake->Fill(aPair->QInv(), weight);
    fDenominator->Fill(aPair->QInv(), 1.0);
    
    Double_t tQinvTrue = GetQinvTrue(aPair);
    
    fNumeratorFakeIdeal->Fill(tQinvTrue, weight);
    fDenominatorIdeal->Fill(tQinvTrue, 1.0);
    
    fQgenQrec->Fill(tQinvTrue,aPair->QInv());
  }
  //Special MC analysis for K selected by PDG code -->
  else {
    Double_t weight = fManager->GetWeight(aPair);
    AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
    AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();
    Double_t pdg1 = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetPDGPid();
    Double_t pdg2 = ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetPDGPid();
    // if((aPair->KT())<0.5)cout<<" Corr Func  pdg1 "<<pdg1<<" pdg2 "<<pdg2<<" qinv "<<aPair->QInv()<< " w "<<weight<<endl;
    fNumeratorFake->Fill(aPair->QInv(), weight);
    fDenominator->Fill(aPair->QInv(), 1.0);
    Double_t tQinvTrue = GetQinvTrue(aPair);
    if(tQinvTrue>0)fNumeratorFakeIdeal->Fill(tQinvTrue, weight);
    if(tQinvTrue>0)fDenominatorIdeal->Fill(tQinvTrue, 1.0);
    if(tQinvTrue>0)fQgenQrec->Fill(tQinvTrue,aPair->QInv());
  }
}

//_______________________
Double_t AliFemtoModelCorrFctn::GetQinvTrue(AliFemtoPair* aPair)
{
  if(!fKaonPDG) {
  AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
  AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();

  AliFemtoLorentzVector fm1;
  AliFemtoThreeVector* temp = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum();
  fm1.SetVect(*temp);
  Double_t am1 = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetMass();
  Double_t am2 = ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetMass();
  double ener = TMath::Sqrt(temp->Mag2()+am1*am1);
  fm1.SetE(ener);

  AliFemtoLorentzVector fm2;
  AliFemtoThreeVector* temp2 =  ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetTrueMomentum();
  fm2.SetVect(*temp2);
  ener = TMath::Sqrt(temp2->Mag2()+am2*am2);
  fm2.SetE(ener);

  //std::cout<<" CFModel mass1 mass2 "<<am1<<" "<<am2<<std::endl;

  AliFemtoLorentzVector tQinvTrueVec = (fm1-fm2);
  Double_t tQinvTrue = -1.* tQinvTrueVec.m();

  return tQinvTrue;
  }
  //Special MC analysis for K selected by PDG code -->
  else {
      AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
  AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();

  AliFemtoLorentzVector fm1;
  AliFemtoThreeVector* temp = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum();
  fm1.SetVect(*temp);
  Double_t am1 = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetMass();
  Double_t am2 = ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetMass();
 
  am1=0.493677;
  am2=0.493677;
 
  //Double_t pdg1 = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetPDGPid();
  //Double_t pdg2 = ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetPDGPid();
 
  double ener = TMath::Sqrt(temp->Mag2()+am1*am1);
  fm1.SetE(ener);

  AliFemtoLorentzVector fm2;
  AliFemtoThreeVector* temp2 =  ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetTrueMomentum();
  fm2.SetVect(*temp2);
  ener = TMath::Sqrt(temp2->Mag2()+am2*am2);
  fm2.SetE(ener);

 
  AliFemtoLorentzVector tQinvTrueVec = (fm1-fm2);
  Double_t tQinvTrue = -1.* tQinvTrueVec.m();
  
 // if(tQinvTrue<0 && am1!=0 && am2!=0)std::cout<<" CFModel Qinv mass1 mass2 "<<aPair->QInv()<<" Qinv_true "<<tQinvTrue<<" "<<am1<<" "<<am2<<" pdg1 "<<pdg1<<" pdg2 "<<pdg2<<std::endl;
 // if(pdg1!=211 || pdg2!=211)std::cout<<" CFModel Qinv mass1 mass2 "<<aPair->QInv()<<" Qinv_true "<<tQinvTrue<<" "<<am1<<" "<<am2<<" pdg1 "<<pdg1<<" pdg2 "<<pdg2<<std::endl;

  if(am1==0 || am2==0)tQinvTrue=-10000;

  return tQinvTrue;
  }
}

//_______________________
void AliFemtoModelCorrFctn::EventBegin(const AliFemtoEvent* /* aEvent */)
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctn::EventEnd(const AliFemtoEvent* /* aEvent */)
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctn::Finish()
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctn::Write()
{
  // Write out data histos

  fQgenQrec->Write();

  fNumeratorTrue->Write();
  fNumeratorFake->Write();
  fDenominator->Write();

  fNumeratorTrueIdeal->Write();
  fNumeratorFakeIdeal->Write();
  fDenominatorIdeal->Write();


}
//_______________________
AliFemtoModelCorrFctn* AliFemtoModelCorrFctn::Clone()
{
  // Create clone
  AliFemtoModelCorrFctn *tCopy = new AliFemtoModelCorrFctn(*this);

  return tCopy;
}
//_________________________
TList* AliFemtoModelCorrFctn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fNumeratorTrue);
  tOutputList->Add(fNumeratorFake);
  tOutputList->Add(fDenominator);

  tOutputList->Add(fNumeratorTrueIdeal);
  tOutputList->Add(fNumeratorFakeIdeal);
  tOutputList->Add(fDenominatorIdeal);
  tOutputList->Add(fQgenQrec);


  return tOutputList;
}
void AliFemtoModelCorrFctn::SetKaonPDG(Bool_t aSetKaonAna)
{
  fKaonPDG = aSetKaonAna;
}
