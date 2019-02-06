/******************************************************************************/
/*                                                                            */
/*  AliFemtoModelCorrFctnKK -   numerator and denominator                     */
/*      histograms from the true Monte Carlo K+K- tracks                      */
/*         Konstantin.Mikhaylov@cern.ch                                       */
/*                                                                            */
/******************************************************************************/

#ifdef __ROOT__
   ClassImp(AliFemtoModelCorrFctnKK, 1);
#endif

#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelCorrFctnKK.h"
#include "AliFemtoPair.h"
//#include "AliFemtoModelManager.h"
//#include <TH1D.h>
//#include <TH2D.h>

//_______________________
AliFemtoModelCorrFctnKK::AliFemtoModelCorrFctnKK():
  AliFemtoModelCorrFctn(),
  fManager(0),
  fNumeratorTrue(0),
  fNumeratorFake(0),
  fDenominator(0),
  fNumeratorTrueIdeal(0),
  fNumeratorFakeIdeal(0),
  fDenominatorIdeal(0),
  fQgenQrec(0),
  fdP(0),
  fdPt(0),
  fdPtvsPt(0),
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

  fdP = new TH1D("Delta_P","Delta_P", 50,-0.1,0.1);
  fdPt = new TH1D("Delta_Pt","Delta_Pt", 50,-0.1,0.1);
  fdPtvsPt = new TH2D("Delta_PtvsPt","Delta_PtvsPt",50,0.0,2.0,50,-0.1,0.1);

  /*
    for(int i=0;i<fNbbPairs;i++){
        fkTdists[i] = new TH1D(Form("fkTdists[%i]",i),Form("fkTdists[%i]",i),100,0.0,5.0);
        fkTdists[i]->Sumw2();
    }
  */

  fNumeratorTrue->Sumw2();
  fNumeratorFake->Sumw2();
  fDenominator->Sumw2();

  fNumeratorTrueIdeal->Sumw2();
  fNumeratorFakeIdeal->Sumw2();
  fDenominatorIdeal->Sumw2();

  fQgenQrec->Sumw2();

  fdP->Sumw2();
  fdPt->Sumw2();

}
//_______________________
AliFemtoModelCorrFctnKK::AliFemtoModelCorrFctnKK(const char *title, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi):
  AliFemtoModelCorrFctn(),
  fManager(0),
  fNumeratorTrue(0),
  fNumeratorFake(0),
  fDenominator(0),
  fNumeratorTrueIdeal(0),
  fNumeratorFakeIdeal(0),
  fDenominatorIdeal(0),
  fQgenQrec(0),
  fdP(0),
  fdPt(0),
  fdPtvsPt(0),
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
  //test
  //fQgenQrec = new TH2D(buf,buf,aNbins,aQinvLo,aQinvHi,aNbins,-0.05,0.05);

  snprintf(buf , 100,  "DeltaP_%s", title);
  fdP = new TH1D(buf,buf,50,-0.1,0.1);
  snprintf(buf , 100,  "DeltaPT_%s", title);
  fdPt = new TH1D(buf,buf,50,-0.1,0.1);
  snprintf(buf , 100,  "DeltaPT_PT_%s", title);
  fdPtvsPt = new TH2D(buf,buf,50,0.0,2.0,50,-0.1,0.1);

  /*
    for(int i=0;i<fNbbPairs;i++){
        fkTdists[i] = new TH1D(Form("fkTdists[%i]_%s",i,title),Form("fkTdists[%i]_%s",i,title),100,0.0,5.0);
        fkTdists[i]->Sumw2();
    }
  */

  fNumeratorTrue->Sumw2();
  fNumeratorFake->Sumw2();
  fDenominator->Sumw2();

  fNumeratorTrueIdeal->Sumw2();
  fNumeratorFakeIdeal->Sumw2();
  fDenominatorIdeal->Sumw2();

  fQgenQrec->Sumw2();

  fdP->Sumw2();
  fdPt->Sumw2();
  fdPtvsPt->Sumw2();

}
//_______________________
AliFemtoModelCorrFctnKK::AliFemtoModelCorrFctnKK(const AliFemtoModelCorrFctnKK& aCorrFctn) :
  AliFemtoModelCorrFctn(aCorrFctn),
  fManager(0),
  fNumeratorTrue(0),
  fNumeratorFake(0),
  fDenominator(0),
  fNumeratorTrueIdeal(0),
  fNumeratorFakeIdeal(0),
  fDenominatorIdeal(0),
  fQgenQrec(0),
  fdP(0),
  fdPt(0),
  fdPtvsPt(0),
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

  if (aCorrFctn.fdP)
    fdP = new TH1D(*(aCorrFctn.fdP));
  if (aCorrFctn.fdPt)
    fdPt = new TH1D(*(aCorrFctn.fdPt));
   if (aCorrFctn.fdPtvsPt)
    fdPtvsPt = new TH2D(*(aCorrFctn.fdPtvsPt));

  /*
    for(int i=0;i<fNbbPairs;i++){
        if(aCorrFctn.fkTdists[i])
            fkTdists[i] = aCorrFctn.fkTdists[i];
    }
  */

  fManager = aCorrFctn.fManager;
}
//_______________________
AliFemtoModelCorrFctnKK::~AliFemtoModelCorrFctnKK()
{
  // Destructor
  if (fNumeratorTrue) delete fNumeratorTrue;
  if (fNumeratorFake) delete fNumeratorFake;
  if (fDenominator) delete fDenominator;

  if (fNumeratorTrueIdeal) delete fNumeratorTrueIdeal;
  if (fNumeratorFakeIdeal) delete fNumeratorFakeIdeal;
  if (fDenominatorIdeal) delete fDenominatorIdeal;

  if (fQgenQrec) delete fQgenQrec;

  if (fdP) delete fdP;
  if (fdPt) delete fdPt;
  if (fdPtvsPt) delete fdPtvsPt;

  /*
    for(int i=0;i<fNbbPairs;i++){
        if(fkTdists[i]) delete fkTdists[i];
    }
  */

}
//_______________________
AliFemtoModelCorrFctnKK& AliFemtoModelCorrFctnKK::operator=(const AliFemtoModelCorrFctnKK& aCorrFctn)
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

  delete fdP;
  fdP = (aCorrFctn.fdP)
            ? new TH1D(*aCorrFctn.fdP)
            : nullptr;

  delete fdPt;
  fdPt = (aCorrFctn.fdPt)
            ? new TH1D(*aCorrFctn.fdPt)
            : nullptr;

  delete fdPtvsPt;
  fdPtvsPt = (aCorrFctn.fdPtvsPt)
            ? new TH2D(*aCorrFctn.fdPtvsPt)
            : nullptr;

  /*
  for(int i=0;i<fNbbPairs;i++){
        delete fkTdists[i];
        fkTdists[i] = (aCorrFctn.fkTdists[i])
        ? new TH1D(*aCorrFctn.fkTdists[i])
        : nullptr;

    }
  */
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
void AliFemtoModelCorrFctnKK::ConnectToManager(AliFemtoModelManager *aManager)
{
  fManager = aManager;
}

//_______________________
AliFemtoString AliFemtoModelCorrFctnKK::Report()
{
  // Prepare report
  AliFemtoString tStr = "AliFemtoModelCorrFctnKK report";

  return tStr;
}

//_______________________
void AliFemtoModelCorrFctnKK::AddRealPair(AliFemtoPair* aPair)
{
  if (fPairCut && !fPairCut->Pass(aPair)) {
    return;
  }

  if(!fKaonPDG) {
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

    AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
    AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();
    //true and recontructed momntum to get delta_p/p -->
    Double_t pdg1 = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetPDGPid();
    Double_t pdg2 = ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetPDGPid();

    //if(pdg1 != 321 || pdg2 != -321) cout<<"____________ ++++Kaons : pdg1="<<pdg1<<" pdg2="<<pdg2<<endl;

    if(pdg1 == 321 && pdg2 == -321) { //to be sure PID does not change MR

    AliFemtoLorentzVector p_true_1;//true momentum
    AliFemtoThreeVector* temp1 =
      ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum();
    p_true_1.SetVect(*temp1);
    AliFemtoLorentzVector p_true_2;//true momentum
    AliFemtoThreeVector* temp2 =
      ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetTrueMomentum();
    p_true_2.SetVect(*temp2);
    //cout<<"###_____________________ True: Px="<<p_true_1.x()<<endl;
    // //AliFemtoLorentzVector p_rec_1 = aPair->Track1()->FourMomentum();//reconstructed momentum
    //AliFemtoThreeVector p_rec_1 = aPair->Track1()->Track()->P();
    AliFemtoLorentzVector p_rec_1 = aPair->Track1()->FourMomentum();
    AliFemtoLorentzVector p_rec_2 = aPair->Track2()->FourMomentum();
    Double_t P, deltaP, Pt, deltaPt;

    if(fP1x != p_true_1.x()) {
      // cout<<"#+1>_____________________ True: x="<<p_true_1.x()<<"  y="<<p_true_1.y()<<" z="<<p_true_1.z()<<
      //" Rec: x="<<p_rec_1.x()<<"  y="<<p_rec_1.y()<<" z="<<p_rec_1.z()<<endl;
    fP1x = p_true_1.x();
    P = TMath::Sqrt(p_true_1.x()*p_true_1.x()+p_true_1.y()*p_true_1.y()+p_true_1.z()*p_true_1.z());
    Pt = TMath::Sqrt(p_true_1.x()*p_true_1.x()+p_true_1.y()*p_true_1.y());
    //cout<<" P= "<<P<<endl;
    if( P !=0 && Pt !=0 ){
    deltaP=
      TMath::Sqrt(p_true_1.x()*p_true_1.x()+p_true_1.y()*p_true_1.y()+p_true_1.z()*p_true_1.z())-
      TMath::Sqrt(p_rec_1.x()*p_rec_1.x()+p_rec_1.y()*p_rec_1.y()+p_rec_1.z()*p_rec_1.z());
    //cout<<" deltaP= "<<deltaP<<endl;
    deltaPt=
      TMath::Sqrt(p_true_1.x()*p_true_1.x()+p_true_1.y()*p_true_1.y())-
      TMath::Sqrt(p_rec_1.x()*p_rec_1.x()+p_rec_1.y()*p_rec_1.y());
    fdP->Fill(deltaP/P);
    fdPt->Fill(deltaPt/Pt);
    fdPtvsPt->Fill(Pt,deltaPt/Pt);
    }
    }


    if(fP2x != p_true_2.x()) {
      //cout<<"#+2>_____________________ True: x="<<p_true_2.x()<<"  y="<<p_true_2.y()<<" z="<<p_true_2.z()<<
      //" Rec: x="<<p_rec_2.x()<<"  y="<<p_rec_2.y()<<" z="<<p_rec_2.z()<<endl;
    fP2x = p_true_2.x();

    P = TMath::Sqrt(p_true_2.x()*p_true_2.x()+p_true_2.y()*p_true_2.y()+p_true_2.z()*p_true_2.z());
    Pt = TMath::Sqrt(p_true_2.x()*p_true_2.x()+p_true_2.y()*p_true_2.y());
    deltaP=
      TMath::Sqrt(p_true_2.x()*p_true_2.x()+p_true_2.y()*p_true_2.y()+p_true_2.z()*p_true_2.z())-
      TMath::Sqrt(p_rec_2.x()*p_rec_2.x()+p_rec_2.y()*p_rec_2.y()+p_rec_2.z()*p_rec_2.z());
    Pt = TMath::Sqrt(p_true_2.x()*p_true_2.x()+p_true_2.y()*p_true_2.y());
    deltaPt=
      TMath::Sqrt(p_true_2.x()*p_true_2.x()+p_true_2.y()*p_true_2.y())-
      TMath::Sqrt(p_rec_2.x()*p_rec_2.x()+p_rec_2.y()*p_rec_2.y());
    fdP->Fill(deltaP/P);
    fdPt->Fill(deltaPt/Pt);
    fdPtvsPt->Fill(Pt,deltaPt/Pt);
    }

    /*Double_t deltaP=TMath::Sqrt(
				(p_true_2.x()-p_rec_2.x())*(p_true_2.x()-p_rec_2.x())+
				(p_true_2.y()-p_rec_2.y())*(p_true_2.y()-p_rec_2.y())+
				(p_true_2.z()-p_rec_2.z())*(p_true_2.z()-p_rec_2.z())
				);
    Double_t deltaPt=TMath::Sqrt(
				(p_true_2.x()-p_rec_2.x())*(p_true_2.x()-p_rec_2.x())+
				(p_true_2.y()-p_rec_2.y())*(p_true_2.y()-p_rec_2.y())
				);*/
    }
    //true and recontructed momntum to get delta_p/p <--

  }
}
//_______________________
void AliFemtoModelCorrFctnKK::AddMixedPair(AliFemtoPair* aPair)
{
  if (fPairCut && !fPairCut->Pass(aPair)) {
    return;
  }

  if(!fKaonPDG) {
    Double_t weight = fManager->GetWeight(aPair);
    fNumeratorFake->Fill(aPair->QInv(), weight);
    fDenominator->Fill(aPair->QInv(), 1.0);

    Double_t tQinvTrue = GetQinvTrue(aPair);

    fNumeratorFakeIdeal->Fill(tQinvTrue, weight);
    fDenominatorIdeal->Fill(tQinvTrue, 1.0);

    /*
      if(fFillkT)
      {
          int pairNumber = GetPairNumber(aPair);

          if(pairNumber >= 0){
              if(fkTdists[pairNumber]){
                  fkTdists[pairNumber]->Fill(GetParentsKt(aPair));
              }
          }
      }
    */
      fQgenQrec->Fill(tQinvTrue,aPair->QInv());
  }
  //Special MC analysis for K selected by PDG code -->
  else {
    Double_t weight = fManager->GetWeight(aPair);
    AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
    AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();
    Double_t pdg1 = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetPDGPid();
    Double_t pdg2 = ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetPDGPid();
    if(pdg1 == 321 && pdg2 == -321) { //to be sure PID does not change MR

    //true and recontructed momntum to get delta_p/p -->
    //AliFemtoLorentzVector p_true_1;//true momentum
    //AliFemtoThreeVector* temp =
    //  ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetTrueMomentum();
    //p_true_1.SetVect(*temp);
    //cout<<"###_____________________ True: Px="<<p_true_1.x()<<endl;
    // //AliFemtoLorentzVector p_rec_1 = aPair->Track1()->FourMomentum();//reconstructed momentum
    // //AliFemtoThreeVector p_rec_1 = aPair->Track1()->Track()->P();
    //AliFemtoLorentzVector p_rec_1 = aPair->Track1()->FourMomentum();
    //cout<<"###_____________________ True: x="<<p_true_1.x()<<"  y="<<p_true_1.y()<<" z="<<p_true_1.z()<<
    //  " Rec: x="<<p_rec_1.x()<<"  y="<<p_rec_1.y()<<" z="<<p_rec_1.z()<<endl;
    //true and recontructed momntum to get delta_p/p <--
    /*
    //not used?
    Double_t pdg1 = ((AliFemtoModelHiddenInfo*)inf1->GetHiddenInfo())->GetPDGPid();
    Double_t pdg2 = ((AliFemtoModelHiddenInfo*)inf2->GetHiddenInfo())->GetPDGPid();
    */
    // if((aPair->KT())<0.5)cout<<" Corr Func  pdg1 "<<pdg1<<" pdg2 "<<pdg2<<" qinv "<<aPair->QInv()<< " w "<<weight<<endl;
    fNumeratorFake->Fill(aPair->QInv(), weight);
    fDenominator->Fill(aPair->QInv(), 1.0);
    Double_t tQinvTrue = GetQinvTrue(aPair);
    if(tQinvTrue>0)fNumeratorFakeIdeal->Fill(tQinvTrue, weight);
    if(tQinvTrue>0)fDenominatorIdeal->Fill(tQinvTrue, 1.0);
    if(tQinvTrue>0)fQgenQrec->Fill(tQinvTrue,aPair->QInv());
    //test
    //if(tQinvTrue>0)fQgenQrec->Fill(tQinvTrue,tQinvTrue-aPair->QInv());
  }
  }
}

//_______________________
Double_t AliFemtoModelCorrFctnKK::GetQinvTrue(AliFemtoPair* aPair)
{
  if(!fKaonPDG) {

      AliFemtoParticle *first = (AliFemtoParticle*)aPair->Track1();
      AliFemtoParticle *second = (AliFemtoParticle*)aPair->Track2();

      if(!first || !second) return -1;

      AliFemtoModelHiddenInfo *inf1 = (AliFemtoModelHiddenInfo*)first->GetHiddenInfo();
      AliFemtoModelHiddenInfo *inf2 = (AliFemtoModelHiddenInfo*)second->GetHiddenInfo();

      if(!inf1 || !inf2){
//          cout<<"no hidden info"<<endl;
          return -1;
      }

      AliFemtoLorentzVector fm1;
      AliFemtoThreeVector* temp = inf1->GetTrueMomentum();
      fm1.SetVect(*temp);
      Double_t am1 = inf1->GetMass();
      Double_t am2 = inf2->GetMass();
      double ener = TMath::Sqrt(temp->Mag2()+am1*am1);
      fm1.SetE(ener);

      AliFemtoLorentzVector fm2;
      AliFemtoThreeVector* temp2 =  inf2->GetTrueMomentum();
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
void AliFemtoModelCorrFctnKK::EventBegin(const AliFemtoEvent* /* aEvent */)
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctnKK::EventEnd(const AliFemtoEvent* /* aEvent */)
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctnKK::Finish()
{
  /* Do nothing */
}
//_______________________
void AliFemtoModelCorrFctnKK::Write()
{
  // Write out data histos

  fQgenQrec->Write();

  fdP->Write();
  fdPt->Write();
  fdPtvsPt->Write();

  fNumeratorTrue->Write();
  fNumeratorFake->Write();
  fDenominator->Write();

  /*  if(fFillkT)
    {
        for(int i=0;i<fNbbPairs;i++){
            fkTdists[i]->Write();
        }
	} */
  fNumeratorTrueIdeal->Write();
  fNumeratorFakeIdeal->Write();
  fDenominatorIdeal->Write();


}
//_______________________
AliFemtoModelCorrFctn* AliFemtoModelCorrFctnKK::Clone() const
{
  // Create clone
  AliFemtoModelCorrFctnKK *tCopy = new AliFemtoModelCorrFctnKK(*this);

  return tCopy;
}
//_________________________
TList* AliFemtoModelCorrFctnKK::GetOutputList()
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

  tOutputList->Add(fdP);
  tOutputList->Add(fdPt);
  tOutputList->Add(fdPtvsPt);

  /*  if(fFillkT)
    {
        for(int i=0;i<fNbbPairs;i++){
            tOutputList->Add(fkTdists[i]);
        }
	} */
  return tOutputList;
}
void AliFemtoModelCorrFctnKK::SetKaonPDG(Bool_t aSetKaonAna)
{
  fKaonPDG = aSetKaonAna;
}

/*
double AliFemtoModelCorrFctnKK::GetParentsKt(AliFemtoPair *pair)
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
*/

int AliFemtoModelCorrFctnKK::GetPairNumber(AliFemtoPair *pair)
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

    int pdg1 = TMath::Abs(info1->GetMotherPdgCode());
    int pdg2 = TMath::Abs(info2->GetMotherPdgCode());

    int tmp;
    if(pdg2 < pdg1){
        tmp = pdg1;
        pdg1 = pdg2;
        pdg2 = tmp;
    }

    delete first;delete second;

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



