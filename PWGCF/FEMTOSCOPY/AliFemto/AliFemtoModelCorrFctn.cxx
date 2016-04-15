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
  fQgenQrec(0)
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
  fQgenQrec(0)
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
  fQgenQrec(0)
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

  if (aCorrFctn.fNumeratorTrue)
    fNumeratorTrue = new TH1D(*(aCorrFctn.fNumeratorTrue));
  else
    fNumeratorTrue = 0;
  if (aCorrFctn.fNumeratorFake)
    fNumeratorFake = new TH1D(*(aCorrFctn.fNumeratorFake));
  else
    fNumeratorFake = 0;
  if (aCorrFctn.fDenominator)
    fDenominator = new TH1D(*(aCorrFctn.fDenominator));
  else
    fDenominator = 0;

  if (aCorrFctn.fQgenQrec)
    fQgenQrec = new TH2D(*(aCorrFctn.fQgenQrec));
  else
    fQgenQrec = 0;

  fManager = aCorrFctn.fManager;


  if (aCorrFctn.fNumeratorTrueIdeal)
    fNumeratorTrueIdeal = new TH1D(*(aCorrFctn.fNumeratorTrueIdeal));
  else
    fNumeratorTrueIdeal = 0;
  if (aCorrFctn.fNumeratorFakeIdeal)
    fNumeratorFakeIdeal = new TH1D(*(aCorrFctn.fNumeratorFakeIdeal));
  else
    fNumeratorFake = 0;
  if (aCorrFctn.fDenominatorIdeal)
    fDenominatorIdeal = new TH1D(*(aCorrFctn.fDenominatorIdeal));
  else
    fDenominatorIdeal = 0;

  fManager = aCorrFctn.fManager;

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
 
 // cout<<" AliFemtoModelCorrFcn add real pair "<<endl;
  Double_t weight = fManager->GetWeight(aPair);
  //cout<<" wight "<< weight<<endl;
  //cout<<"Qinv"<<aPair->QInv()<<endl;
   
  fNumeratorTrue->Fill(aPair->QInv(), weight);

  Double_t tQinvTrue = GetQinvTrue(aPair);
   
  fNumeratorTrueIdeal->Fill(tQinvTrue, weight);
  
  //cout<<"Qinv true"<<tQinvTrue<<endl;
  
}
//_______________________
void AliFemtoModelCorrFctn::AddMixedPair(AliFemtoPair* aPair)
{
  Double_t weight = fManager->GetWeight(aPair);
  fNumeratorFake->Fill(aPair->QInv(), weight);
  fDenominator->Fill(aPair->QInv(), 1.0);

  Double_t tQinvTrue = GetQinvTrue(aPair);

  fNumeratorFakeIdeal->Fill(tQinvTrue, weight);
  fDenominatorIdeal->Fill(tQinvTrue, 1.0);

  fQgenQrec->Fill(tQinvTrue,aPair->QInv());
}

//_______________________
Double_t AliFemtoModelCorrFctn::GetQinvTrue(AliFemtoPair* aPair)
{
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
