////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelCorrFctnWithWeights - the base class for correlation function which    ///
/// uses the model framework and weight generation                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoModelCorrFctnWithWeights, 1);
  /// \endcond
#endif

#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelCorrFctnWithWeights.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include <TH1D.h>
#include <TH2D.h>

//_______________________
AliFemtoModelCorrFctnWithWeights::AliFemtoModelCorrFctnWithWeights(TH2D *filter1, TH2D *filter2):
AliFemtoCorrFctn(),
  filterHist1(filter1),
  filterHist2(filter2),
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
AliFemtoModelCorrFctnWithWeights::AliFemtoModelCorrFctnWithWeights(const char *title, TH2D *filter1, TH2D *filter2, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi):
  AliFemtoCorrFctn(),
  filterHist1(filter1),
  filterHist2(filter2),
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
AliFemtoModelCorrFctnWithWeights::AliFemtoModelCorrFctnWithWeights(const AliFemtoModelCorrFctnWithWeights& aCorrFctn) :
  AliFemtoCorrFctn(),
  filterHist1(0),
  filterHist2(0),
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

  if (aCorrFctn.filterHist1)
    filterHist1 = new TH2D(*aCorrFctn.filterHist1);
  else
    throw std::runtime_error("You try to copy from an object that has no filterHist1!");
  if (aCorrFctn.filterHist2)
    filterHist2 = new TH2D(*aCorrFctn.filterHist2);
  else
    throw std::runtime_error("You try to copy from an object that has no filterHist2!");

  fManager = aCorrFctn.fManager;
}
//_______________________
AliFemtoModelCorrFctnWithWeights::~AliFemtoModelCorrFctnWithWeights()
{
  // Destructor
  if (fNumeratorTrue) delete fNumeratorTrue;
  if (fNumeratorFake) delete fNumeratorFake;
  if (fDenominator) delete fDenominator;

  if (fNumeratorTrueIdeal) delete fNumeratorTrueIdeal;
  if (fNumeratorFakeIdeal) delete fNumeratorFakeIdeal;
  if (fDenominatorIdeal) delete fDenominatorIdeal;

  if (fQgenQrec) delete fQgenQrec;

  if(filterHist1) delete filterHist1;
  if(filterHist2) delete filterHist2;

}
//_______________________
AliFemtoModelCorrFctnWithWeights& AliFemtoModelCorrFctnWithWeights::operator=(const AliFemtoModelCorrFctnWithWeights& aCorrFctn)
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

  if (aCorrFctn.filterHist1)
    filterHist1 = new TH2D(*aCorrFctn.filterHist1);
  else
    throw std::runtime_error("You try to copy from an object that has no filterHist1!");
  if (aCorrFctn.filterHist2)
    filterHist2 = new TH2D(*aCorrFctn.filterHist2);
  else
    throw std::runtime_error("You try to copy from an object that has no filterHist2!");


  return *this;
}
//_______________________
void AliFemtoModelCorrFctnWithWeights::ConnectToManager(AliFemtoModelManager *aManager)
{
  fManager = aManager;
}

//_______________________
AliFemtoString AliFemtoModelCorrFctnWithWeights::Report()
{
  // Prepare report
  AliFemtoString tStr = "AliFemtoModelCorrFctnWithWeights report";

  return tStr;
}

//_______________________
void AliFemtoModelCorrFctnWithWeights::AddRealPair(AliFemtoPair* aPair)
{

 // cout<<" AliFemtoModelCorrFcn add real pair "<<endl;
  Double_t weight = fManager->GetWeight(aPair);

  double y1 = aPair->Track1()->FourMomentum().Rapidity();
  double y2 = aPair->Track2()->FourMomentum().Rapidity();
  double px1 = aPair->Track1()->Track()->P().x();
  double py1 = aPair->Track1()->Track()->P().y();
  //double pz1 = aPair->Track1()->Track()->P().z();

  double px2 = aPair->Track2()->Track()->P().x();
  double py2 = aPair->Track2()->Track()->P().y();
  // double pz2 = aPair->Track2()->Track()->P().z();

  double pt1 = TMath::Hypot(px1, py1);
  double pt2 = TMath::Hypot(px2, py2);

  double weight1 = filterHist1->GetBinContent(filterHist1->FindBin(y1, pt1));
  double weight2 = filterHist2->GetBinContent(filterHist2->FindBin(y2, pt2));
  double totalWeight = weight*weight1*weight2;

  fNumeratorTrue->Fill(aPair->QInv(), totalWeight);

  Double_t tQinvTrue = GetQinvTrue(aPair);

  fNumeratorTrueIdeal->Fill(tQinvTrue, totalWeight);

  //cout<<"Qinv true"<<tQinvTrue<<endl;

}
//_______________________
void AliFemtoModelCorrFctnWithWeights::AddMixedPair(AliFemtoPair* aPair)
{
  Double_t weight = fManager->GetWeight(aPair);
  double y1 = aPair->Track1()->FourMomentum().Rapidity();
  double y2 = aPair->Track2()->FourMomentum().Rapidity();
  double px1 = aPair->Track1()->Track()->P().x();
  double py1 = aPair->Track1()->Track()->P().y();
  //double pz1 = aPair->Track1()->Track()->P().z();

  double px2 = aPair->Track2()->Track()->P().x();
  double py2 = aPair->Track2()->Track()->P().y();
  //double pz2 = aPair->Track2()->Track()->P().z();

  double pt1 = TMath::Hypot(px1, py1);
  double pt2 = TMath::Hypot(px2, py2);

  double weight1 = filterHist1->GetBinContent(filterHist1->FindBin(y1, pt1));
  double weight2 = filterHist2->GetBinContent(filterHist2->FindBin(y2, pt2));
  double totalWeight = weight*weight1*weight2;

  const double qinv = aPair->QInv();

  fNumeratorFake->Fill(qinv, totalWeight);
  fDenominator->Fill(qinv, weight1*weight2);

  Double_t tQinvTrue = GetQinvTrue(aPair);

  fNumeratorFakeIdeal->Fill(tQinvTrue, totalWeight);
  fDenominatorIdeal->Fill(tQinvTrue, weight1*weight2);

  fQgenQrec->Fill(tQinvTrue, qinv);
}

//_______________________
Double_t AliFemtoModelCorrFctnWithWeights::GetQinvTrue(AliFemtoPair* aPair)
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
void AliFemtoModelCorrFctnWithWeights::EventBegin(const AliFemtoEvent* /* aEvent */)
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctnWithWeights::EventEnd(const AliFemtoEvent* /* aEvent */)
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctnWithWeights::Finish()
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctnWithWeights::Write()
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
AliFemtoModelCorrFctnWithWeights* AliFemtoModelCorrFctnWithWeights::Clone() const
{
  // Create clone
  AliFemtoModelCorrFctnWithWeights *tCopy = new AliFemtoModelCorrFctnWithWeights(*this);

  return tCopy;
}
//_________________________
TList* AliFemtoModelCorrFctnWithWeights::GetOutputList()
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
