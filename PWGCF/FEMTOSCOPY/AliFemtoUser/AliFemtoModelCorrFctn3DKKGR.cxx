/******************************************************************************/
/*                                                                            */
/*  AliFemtoModelCorrFctn3DKKGR -   numerator and denominator                 */
/*      3D histograms from the true Monte Carlo KK tracks                     */
/*         gleb.romanenko@cern.ch                                             */
/*                                                                            */
/******************************************************************************/
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoModelCorrFctn3DKKGR, 1);
  /// \endcond
#endif

#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelCorrFctn3DKKGR.h"
#include "AliFemtoPair.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoModelManager.h"
#include <TH3D.h>
#include <TH2D.h>

//------------------
//_______________________
AliFemtoModelCorrFctn3DKKGR::AliFemtoModelCorrFctn3DKKGR()
  : AliFemtoModelCorrFctn()
  , fManager(nullptr)
  , fNumeratorTrue(nullptr)
  , fNumeratorFake(nullptr)
  , fDenominator(nullptr)
  , fNumeratorTrueIdeal(nullptr)
  , fNumeratorFakeIdeal(nullptr)
  , fDenominatorIdeal(nullptr)
  , fQgenQrec(nullptr)
  , fKaonPDG(kFALSE)
  , fFillkT(kFALSE)
{
  // Default constructor
  fNumeratorTrue = new TH3D("ModelNumTrue","ModelNumTrue",50,0.0,0.5,50,0.0,0.5,50,0.0,0.5);
  fNumeratorFake = new TH3D("ModelNumFake","ModelNumFake",50,0.0,0.5,50,0.0,0.5,50,0.0,0.5);
  fDenominator = new TH3D("ModelDen","ModelDen",50,0.0,0.5,50,0.0,0.5,50,0.0,0.5);

  fNumeratorTrueIdeal = new TH3D("ModelNumTrueIdeal","ModelNumTrueIdeal",50,0.0,0.5,50,0.0,0.5,50,0.0,0.5);
  fNumeratorFakeIdeal = new TH3D("ModelNumFakeIdeal","ModelNumFakeIdeal",50,0.0,0.5,50,0.0,0.5,50,0.0,0.5);
  fDenominatorIdeal = new TH3D("ModelDenIdeal","ModelDenIdeal",50,0.0,0.5,50,0.0,0.5,50,0.0,0.5);

  fQgenQrec = new TH2D("QgenQrec","QgenQrec",50,0.0,0.5,50,0.0,0.5);

  for (int i=0;i<fNbbPairs;i++){
    auto *name = Form("fkTdists[%i]", i);
    fkTdists[i] = new TH1D(name, name,100,0.0,5.0);
    fkTdists[i]->Sumw2();
  }

  fNumeratorTrue->Sumw2();
  fNumeratorFake->Sumw2();
  fDenominator->Sumw2();

  fNumeratorTrueIdeal->Sumw2();
  fNumeratorFakeIdeal->Sumw2();
  fDenominatorIdeal->Sumw2();

  fQgenQrec->Sumw2();

}
//_______________________
AliFemtoModelCorrFctn3DKKGR::AliFemtoModelCorrFctn3DKKGR(const char *title,
                                             Int_t aNbins,
                                             Double_t aQinvLo,
                                             Double_t aQinvHi)
  : AliFemtoModelCorrFctn()
  , fManager(nullptr)
  , fNumeratorTrue(nullptr)
  , fNumeratorFake(nullptr)
  , fDenominator(nullptr)
  , fNumeratorTrueIdeal(nullptr)
  , fNumeratorFakeIdeal(nullptr)
  , fDenominatorIdeal(nullptr)
  , fQgenQrec(nullptr)
  , fKaonPDG(kFALSE)
  , fFillkT(kFALSE)
{
  // Normal constructor
  char *buf;
  buf = Form("NumTrue%s", title);
  fNumeratorTrue = new TH3D(buf,buf,aNbins,aQinvLo,aQinvHi,aNbins,aQinvLo,aQinvHi,aNbins,aQinvLo,aQinvHi);

  buf = Form("NumFake%s", title);
  fNumeratorFake = new TH3D(buf,buf,aNbins,aQinvLo,aQinvHi,aNbins,aQinvLo,aQinvHi,aNbins,aQinvLo,aQinvHi);

  buf = Form("Den%s", title);
  fDenominator = new TH3D(buf,buf,aNbins,aQinvLo,aQinvHi,aNbins,aQinvLo,aQinvHi,aNbins,aQinvLo,aQinvHi);

  buf = Form("NumTrueIdeal%s", title);
  fNumeratorTrueIdeal = new TH3D(buf,buf,aNbins,aQinvLo,aQinvHi,aNbins,aQinvLo,aQinvHi,aNbins,aQinvLo,aQinvHi);

  buf = Form("NumFakeIdeal%s", title);
  fNumeratorFakeIdeal = new TH3D(buf,buf,aNbins,aQinvLo,aQinvHi,aNbins,aQinvLo,aQinvHi,aNbins,aQinvLo,aQinvHi);

  buf = Form("DenIdeal%s", title);
  fDenominatorIdeal = new TH3D(buf,buf,aNbins,aQinvLo,aQinvHi,aNbins,aQinvLo,aQinvHi,aNbins,aQinvLo,aQinvHi);

  buf = Form("QgenQrec%s", title);
  fQgenQrec = new TH2D(buf,buf,aNbins,aQinvLo,aQinvHi,aNbins,aQinvLo,aQinvHi);
  //test
  //fQgenQrec = new TH2D(buf,buf,aNbins,aQinvLo,aQinvHi,aNbins,-0.05,0.05);

  for (int i=0;i<fNbbPairs;i++){
    const char *name = Form("fkTdists[%i]_%s",i,title);
    fkTdists[i] = new TH1D(name, name, 100, 0.0, 5.0);
    fkTdists[i]->Sumw2();
  }

  fNumeratorTrue->Sumw2();
  fNumeratorFake->Sumw2();
  fDenominator->Sumw2();

  fNumeratorTrueIdeal->Sumw2();
  fNumeratorFakeIdeal->Sumw2();
  fDenominatorIdeal->Sumw2();

  fQgenQrec->Sumw2();
}

//_______________________
AliFemtoModelCorrFctn3DKKGR::AliFemtoModelCorrFctn3DKKGR(const AliFemtoModelCorrFctn3DKKGR& aCorrFctn)
  : AliFemtoModelCorrFctn(aCorrFctn)
  , fManager(aCorrFctn.fManager)
  , fNumeratorTrue(nullptr)
  , fNumeratorFake(nullptr)
  , fDenominator(nullptr)
  , fNumeratorTrueIdeal(nullptr)
  , fNumeratorFakeIdeal(nullptr)
  , fDenominatorIdeal(nullptr)
  , fQgenQrec(nullptr)
  , fKaonPDG(aCorrFctn.fKaonPDG)
  , fFillkT(aCorrFctn.fFillkT)
{
  // Copy constructor
  fNumeratorTrue = new TH3D(*aCorrFctn.fNumeratorTrue);
  fNumeratorFake = new TH3D(*aCorrFctn.fNumeratorFake);
  fDenominator = new TH3D(*aCorrFctn.fDenominator);

  fNumeratorTrueIdeal = new TH3D(*aCorrFctn.fNumeratorTrueIdeal);
  fNumeratorFakeIdeal = new TH3D(*aCorrFctn.fNumeratorFakeIdeal);
  fDenominatorIdeal = new TH3D(*aCorrFctn.fDenominatorIdeal);

  fQgenQrec = new TH2D(*aCorrFctn.fQgenQrec);

  for (int i=0;i<fNbbPairs;i++) {
    fkTdists[i] = new TH1D(*aCorrFctn.fkTdists[i]);
  }
}

//_______________________
AliFemtoModelCorrFctn3DKKGR::~AliFemtoModelCorrFctn3DKKGR()
{
  // Destructor
  delete fNumeratorTrue;
  delete fNumeratorFake;
  delete fDenominator;

  delete fNumeratorTrueIdeal;
  delete fNumeratorFakeIdeal;
  delete fDenominatorIdeal;

  delete fQgenQrec;

  for (int i=0; i<fNbbPairs; i++) {
    delete fkTdists[i];
  }

}
//_______________________
AliFemtoModelCorrFctn3DKKGR& AliFemtoModelCorrFctn3DKKGR::operator=(const AliFemtoModelCorrFctn3DKKGR& aCorrFctn)
{
  // Assignment operator
  if (this == &aCorrFctn) {
    return *this;
  }

  AliFemtoCorrFctn::operator=(aCorrFctn);

  *fNumeratorTrue = *aCorrFctn.fNumeratorTrue;
  *fNumeratorFake = *aCorrFctn.fNumeratorFake;
  *fDenominator = *aCorrFctn.fDenominator;
  *fQgenQrec = *aCorrFctn.fQgenQrec;

  for (int i=0; i<fNbbPairs; i++) {
    *fkTdists[i] = *aCorrFctn.fkTdists[i];
  }

  *fNumeratorTrueIdeal = *aCorrFctn.fNumeratorTrueIdeal;
  *fNumeratorFakeIdeal = *aCorrFctn.fNumeratorFakeIdeal;
  *fDenominatorIdeal = *aCorrFctn.fDenominatorIdeal;

  fManager = aCorrFctn.fManager;
  fKaonPDG = aCorrFctn.fKaonPDG;

  return *this;
}
//_______________________
void AliFemtoModelCorrFctn3DKKGR::ConnectToManager(AliFemtoModelManager *aManager)
{
  fManager = aManager;
}

//_______________________
AliFemtoString AliFemtoModelCorrFctn3DKKGR::Report()
{
  // Prepare report
  AliFemtoString tStr = "AliFemtoModelCorrFctn3DKKGR report";

  return tStr;
}

//_______________________
void AliFemtoModelCorrFctn3DKKGR::AddRealPair(AliFemtoPair* aPair)
{
  if (fPairCut && !fPairCut->Pass(aPair)) {
    return;
  }

  if(!fKaonPDG) {
    // cout<<" AliFemtoModelCorrFcn add real pair "<<endl;
    Double_t weight = fManager->GetWeight(aPair);
    //cout<<" wight "<< weight<<endl;
    //cout<<"Qinv"<<aPair->QInv()<<endl;

    fNumeratorTrue->Fill(aPair->QOutCMS(),aPair->QSideCMS(),aPair->QLongCMS(), weight);  //qOut, qSide, qLong

    //Double_t tQinvTrue = GetQinvTrue(aPair);

    fNumeratorTrueIdeal->Fill(GetQoutTrue(aPair),GetQsideTrue(aPair),GetQlongTrue(aPair), weight);  //qOut_ideal, qSide_ideal, qLong_ideal

    //cout<<"Qinv true"<<tQinvTrue<<endl;
  }
  //Special MC analysis for K selected by PDG code -->
  else {
    Double_t weight = fManager->GetWeight(aPair);
    fNumeratorTrue->Fill(aPair->QOutCMS(),aPair->QSideCMS(),aPair->QLongCMS(), weight);  //qOut, qSide, qLong
    //Double_t tQinvTrue = GetQinvTrue(aPair);
    fNumeratorTrueIdeal->Fill(GetQoutTrue(aPair),GetQsideTrue(aPair),GetQlongTrue(aPair), weight);  //qOut_ideal, qSide_ideal, qLong_ideal
  }
}
//_______________________
void AliFemtoModelCorrFctn3DKKGR::AddMixedPair(AliFemtoPair* aPair)
{
  if (fPairCut && !fPairCut->Pass(aPair)) {
    return;
  }

  const Double_t
    weight = fManager->GetWeight(aPair),
    qinv = aPair->QInv(),
    qinv_ideal = GetQinvTrue(aPair);

  if(!fKaonPDG) {
    fNumeratorFake->Fill(aPair->QOutCMS(),aPair->QSideCMS(),aPair->QLongCMS(), weight);  //qOut, qSide, qLong
    fDenominator->Fill(aPair->QOutCMS(),aPair->QSideCMS(),aPair->QLongCMS(), 1.0);  //qOut, qSide, qLong

    fNumeratorFakeIdeal->Fill(GetQoutTrue(aPair),GetQsideTrue(aPair),GetQlongTrue(aPair), weight);  //qOut_ideal, qSide_ideal, qLong_ideal
    fDenominatorIdeal->Fill(GetQoutTrue(aPair),GetQsideTrue(aPair),GetQlongTrue(aPair), 1.0);  //qOut_ideal, qSide_ideal, qLong_ideal

    if (fFillkT) {
        int pairNumber = GetPairNumber(aPair);

        if(pairNumber >= 0){
            if(fkTdists[pairNumber]){
                fkTdists[pairNumber]->Fill(GetParentsKt(aPair));
            }
        }
    }

    fQgenQrec->Fill(qinv_ideal, qinv);
  }
  //Special MC analysis for K selected by PDG code -->
  else {
    // AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
    // AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();
    // Double_t pdg1 = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetPDGPid();
    // Double_t pdg2 = ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetPDGPid();
    // if((aPair->KT())<0.5)cout<<" Corr Func  pdg1 "<<pdg1<<" pdg2 "<<pdg2<<" qinv "<<aPair->QInv()<< " w "<<weight<<endl;
    fNumeratorFake->Fill(aPair->QOutCMS(),aPair->QSideCMS(),aPair->QLongCMS(), weight);  //qOut, qSide, qLong
    fDenominator->Fill(aPair->QOutCMS(),aPair->QSideCMS(),aPair->QLongCMS(), 1.0);  //qOut, qSide, qLong
    //if(qinv_ideal>0) {
      fNumeratorFakeIdeal->Fill(GetQoutTrue(aPair),GetQsideTrue(aPair),GetQlongTrue(aPair), weight);  //qOut_ideal, qSide_ideal, qLong_ideal
      fDenominatorIdeal->Fill(GetQoutTrue(aPair),GetQsideTrue(aPair),GetQlongTrue(aPair), 1.0);  //qOut_ideal, qSide_ideal, qLong_ideal
      fQgenQrec->Fill(qinv_ideal, qinv);
    //}
    //test
    //if(tQinvTrue>0)fQgenQrec->Fill(tQinvTrue,tQinvTrue-aPair->QInv());
  }
}

void AliFemtoModelCorrFctn3DKKGR::SetSpecificPairCut(AliFemtoPairCut* aCut)
{
  fPairCut = aCut;
}

//_______________________
Double_t AliFemtoModelCorrFctn3DKKGR::GetQinvTrue(AliFemtoPair* aPair)
{
  if(!fKaonPDG) {

      const AliFemtoParticle *first = aPair->Track1(),
                             *second = aPair->Track2();

      if(!first || !second) return -1;

      const AliFemtoModelHiddenInfo *inf1 = (AliFemtoModelHiddenInfo*)first->GetHiddenInfo(),
                                    *inf2 = (AliFemtoModelHiddenInfo*)second->GetHiddenInfo();

      if(!inf1 || !inf2){
//          cout<<"no hidden info"<<endl;
          return -1;
      }

      AliFemtoThreeVector* temp = inf1->GetTrueMomentum();
      Double_t am1 = inf1->GetMass();
      double ener = TMath::Sqrt(temp->Mag2()+am1*am1);
      AliFemtoLorentzVector fm1(ener, *temp);

      AliFemtoThreeVector* temp2 = inf2->GetTrueMomentum();
      const Double_t am2 = inf2->GetMass();
      ener = TMath::Sqrt(temp2->Mag2()+am2*am2);
      AliFemtoLorentzVector fm2(ener, *temp2);

      //std::cout<<" CFModel mass1 mass2 "<<am1<<" "<<am2<<std::endl;

      AliFemtoLorentzVector tQinvTrueVec = (fm1-fm2);
      Double_t tQinvTrue = -1.* tQinvTrueVec.m();

      return tQinvTrue;
  }
  //Special MC analysis for K selected by PDG code -->
  else {
    const AliFemtoTrack *inf1 = aPair->Track1()->Track(),
                        *inf2 = aPair->Track2()->Track();

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

Double_t AliFemtoModelCorrFctn3DKKGR::GetQoutTrue(AliFemtoPair* aPair)
{
  const AliFemtoParticle *first = aPair->Track1(),
                         *second = aPair->Track2();

  if(!first || !second) return -1;

  const AliFemtoModelHiddenInfo *inf1 = (AliFemtoModelHiddenInfo*)first->GetHiddenInfo(),
                                *inf2 = (AliFemtoModelHiddenInfo*)second->GetHiddenInfo();

  if(!inf1 || !inf2){
    //          cout<<"no hidden info"<<endl;
    return -1;
  }

  AliFemtoThreeVector* temp = inf1->GetTrueMomentum();

  AliFemtoThreeVector* temp2 = inf2->GetTrueMomentum();

  double dx = temp->x() - temp2->x();
  double dy = temp->y() - temp2->y();

  double px = temp->x() + temp2->x();
  double py = temp->y() + temp2->y();

  double qTkT = dx*px + dy*py;
  double pT = sqrt(px*px + py*py);
  double qOut;

  if(!pT)
    qOut=0;
  else
    qOut=qTkT/pT;

  if(!qOut)
    qOut=-10000;

  return qOut;
}

Double_t AliFemtoModelCorrFctn3DKKGR::GetQsideTrue(AliFemtoPair* aPair)
{
  const AliFemtoParticle *first = aPair->Track1(),
                         *second = aPair->Track2();

  if(!first || !second) return -1;

  const AliFemtoModelHiddenInfo *inf1 = (AliFemtoModelHiddenInfo*)first->GetHiddenInfo(),
                                *inf2 = (AliFemtoModelHiddenInfo*)second->GetHiddenInfo();

  if(!inf1 || !inf2){
    //          cout<<"no hidden info"<<endl;
    return -1;
  }

  AliFemtoThreeVector* temp = inf1->GetTrueMomentum();

  AliFemtoThreeVector* temp2 = inf2->GetTrueMomentum();

  double x1 = temp->x();
  double x2 = temp2->x();
  double y1 = temp->y();
  double y2 = temp2->y();

  double px = temp->x() + temp2->x();
  double py = temp->y() + temp2->y();

  double qTkT = x2*y1 - x1*y2;
  double pT = sqrt(px*px + py*py);
  double qSide;

  if(!pT)
    qSide=0;
  else
    qSide=(2.0*qTkT)/pT;

  if(!qSide)
    qSide=-10000;

  return qSide;
}

Double_t AliFemtoModelCorrFctn3DKKGR::GetQlongTrue(AliFemtoPair* aPair)
{

  const AliFemtoParticle *first = aPair->Track1(),
                         *second = aPair->Track2();

  if(!first || !second) return -1;

  const AliFemtoModelHiddenInfo *inf1 = (AliFemtoModelHiddenInfo*)first->GetHiddenInfo(),
                                *inf2 = (AliFemtoModelHiddenInfo*)second->GetHiddenInfo();

  if(!inf1 || !inf2){
    //          cout<<"no hidden info"<<endl;
    return -1;
  }

  AliFemtoThreeVector* temp = inf1->GetTrueMomentum();
  Double_t am1 = inf1->GetMass();
  double ener = TMath::Sqrt(temp->Mag2()+am1*am1);
  AliFemtoLorentzVector fm1(ener, *temp);

  AliFemtoThreeVector* temp2 = inf2->GetTrueMomentum();
  const Double_t am2 = inf2->GetMass();
  ener = TMath::Sqrt(temp2->Mag2()+am2*am2);
  AliFemtoLorentzVector fm2(ener, *temp2);

  //std::cout<<" CFModel mass1 mass2 "<<am1<<" "<<am2<<std::endl;

  double dz = fm1.z() - fm2.z();
  double zz = fm1.z() + fm2.z();

  double dt = fm1.t() - fm2.t();
  double tt = fm1.t() + fm2.t();

  double beta = zz/tt;
  double gamma = 1.0/TMath::Sqrt((1.-beta)*(1.+beta));

  double qLong = gamma*(dz - beta*dt);

  if(!qLong)
    qLong=-10000;
  
  return qLong;
}

//_______________________
void AliFemtoModelCorrFctn3DKKGR::EventBegin(const AliFemtoEvent* /* aEvent */)
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctn3DKKGR::EventEnd(const AliFemtoEvent* /* aEvent */)
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctn3DKKGR::Finish()
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctn3DKKGR::Write()
{
  // Write out data histos

  fQgenQrec->Write();

  fNumeratorTrue->Write();
  fNumeratorFake->Write();
  fDenominator->Write();

    if(fFillkT)
    {
        for(int i=0;i<fNbbPairs;i++){
            fkTdists[i]->Write();
        }
    }
  fNumeratorTrueIdeal->Write();
  fNumeratorFakeIdeal->Write();
  fDenominatorIdeal->Write();


}
//_________________________
TList* AliFemtoModelCorrFctn3DKKGR::GetOutputList()
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

    if(fFillkT)
    {
        for(int i=0;i<fNbbPairs;i++){
            tOutputList->Add(fkTdists[i]);
        }
    }
  return tOutputList;
}
void AliFemtoModelCorrFctn3DKKGR::SetKaonPDG(Bool_t aSetKaonAna)
{
  fKaonPDG = aSetKaonAna;
}

double AliFemtoModelCorrFctn3DKKGR::GetParentsKt(AliFemtoPair *pair)
{
    AliFemtoParticle *first = new AliFemtoParticle(*(pair->Track1()));
    AliFemtoParticle *second = new AliFemtoParticle(*(pair->Track2()));

    if(!first)
    {
        if(second) delete second;
        return -1;
    }
    if(!second)
    {
        if(first) delete first;
        return -1;
    }

    AliFemtoModelHiddenInfo *info1 = (AliFemtoModelHiddenInfo*)first->GetHiddenInfo();
    AliFemtoModelHiddenInfo *info2 = (AliFemtoModelHiddenInfo*)second->GetHiddenInfo();

    if(!info1 || !info2)
    {
        if(first) delete first;
        if(second) delete second;
        return -1;
    }
    AliFemtoThreeVector* p1 = info1->GetMotherMomentum();
    AliFemtoThreeVector* p2 = info2->GetMotherMomentum();

    if(!p1 || !p2)
    {
        if(first) delete first;
        if(second) delete second;
        return -1;
    }
    double px = p1->x() + p2->x();
    double py = p1->y() + p2->y();
    double pT = sqrt(px*px + py*py);

    delete first;delete second;

    return pT/2.;
}

int AliFemtoModelCorrFctn3DKKGR::GetPairNumber(AliFemtoPair *pair)
{
    const AliFemtoModelHiddenInfo
        *info1 = (AliFemtoModelHiddenInfo*)pair->Track1()->GetHiddenInfo(),
        *info2 = (AliFemtoModelHiddenInfo*)pair->Track2()->GetHiddenInfo();

    if (!info1 || !info2) {
        return -1;
    }

    int pdg1 = TMath::Abs(info1->GetMotherPdgCode());
    int pdg2 = TMath::Abs(info2->GetMotherPdgCode());

    if(pdg2 < pdg1){
        std::swap(pdg1, pdg2);
    }

    if(pdg1 == 2212 && pdg2 == 2212) return 0; // pp
    if(pdg1 == 2212 && pdg2 == 3122) return 1; // pΛ
    if(pdg1 == 3122 && pdg2 == 3122) return 2; // ΛΛ
    if(pdg1 == 2212 && pdg2 == 3222) return 3; // pΣ+
    if(pdg1 == 3122 && pdg2 == 3222) return 4; // ΛΣ+
    if(pdg1 == 3222 && pdg2 == 3222) return 5; // Σ+Σ+
    if(pdg1 == 2212 && pdg2 == 3312) return 6; // pΞ-
    if(pdg1 == 2212 && pdg2 == 3322) return 7; // pΞ0
    if(pdg1 == 3122 && pdg2 == 3312) return 8; // ΛΞ-
    if(pdg1 == 3122 && pdg2 == 3322) return 9; // ΛΞ0
    if(pdg1 == 3222 && pdg2 == 3322) return 10; // Σ+Ξ0
    if(pdg1 == 3222 && pdg2 == 3312) return 11; // Σ+Ξ-
    if(pdg1 == 3212 && pdg2 == 3122) return 12; // Σ0Λ
    if(pdg1 == 2212 && pdg2 == 3212) return 13; // pΣ0
    if(pdg1 == 3212 && pdg2 == 3222) return 14; // Σ0Σ+
    if(pdg1 == 3322 && pdg2 == 3322) return 15; // Ξ0Ξ0
    if(pdg1 == 3312 && pdg2 == 3322) return 16; // Ξ-Ξ0
    if(pdg1 == 3312 && pdg2 == 3312) return 17; // Ξ-Ξ-
    if(pdg1 == 3212 && pdg2 == 3322) return 18; // Σ0Ξ0
    if(pdg1 == 3212 && pdg2 == 3312) return 19; // Σ0Ξ-
    if(pdg1 == 3212 && pdg2 == 3212) return 20; // Σ0Σ0

    return -1;
}
