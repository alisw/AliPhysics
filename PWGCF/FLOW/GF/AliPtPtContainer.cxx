/*
Author: Emil Gorm Nielsen
Modified from AliCkContainer by Vytautas
Container to store the terms required to calculate higher order moments of the pT spectrum.
If used, modified, or distributed, please aknowledge the original author of this code.
*/
#include "AliPtPtContainer.h"

Double_t AliPtPtContainer::fFactorial[9] = {1.,1.,2.,6.,24.,120.,720.,5040.,40320.};
Int_t AliPtPtContainer::fSign[9] = {1,-1,1,-1,1,-1,1,-1,1};

using namespace PtPtSpace;
AliPtPtContainer::AliPtPtContainer():
    fCMTermList(0),
    fCorrList(0),
    fSubList(0),
    fSubCMList(0),
    fCentralMomentList(0),
    fCumulantList(0),
    fSubCentralMomentList(0),
    fSubCumulantList(0),
    mpar(0),
    fEventWeight(kEventWeight::kUnity),
    fCorr(),
    fSumw(),
    fCorr1(),
    fSumw1(),
    fCorr2(),
    fSumw2(),
    fcmNum(),
    fcmDen()
{};
AliPtPtContainer::~AliPtPtContainer()
{
    delete fCMTermList;
    delete fCorrList;
    delete fSubList;
    delete fSubCMList;
    delete fCumulantList;
    delete fCentralMomentList;
    delete fSubCentralMomentList;
    delete fSubCumulantList;
};
AliPtPtContainer::AliPtPtContainer(const char* name, const char* title, int nbinsx, double* xbins, int m):
    TNamed(name,title),
    fCMTermList(0),
    fCorrList(0),
    fSubList(0),
    fSubCMList(0),
    fCentralMomentList(0),
    fCumulantList(0),
    fSubCentralMomentList(0),
    fSubCumulantList(0),
    mpar(m),
    fEventWeight(kEventWeight::kUnity),
    fCorr(),
    fSumw(),
    fCorr1(),
    fSumw1(),
    fCorr2(),
    fSumw2(),
    fcmNum(),
    fcmDen()
{
    Initialize(nbinsx,xbins);
};
AliPtPtContainer::AliPtPtContainer(const char* name, const char* title, int nbinsx, double xlow, double xhigh, int m):
    TNamed(name,title),
    fCMTermList(0),
    fCorrList(0),
    fSubList(0),
    fSubCMList(0),
    fCentralMomentList(0),
    fCumulantList(0),
    fSubCentralMomentList(0),
    fSubCumulantList(0),
    mpar(m),
    fEventWeight(kEventWeight::kUnity),
    fCorr(),
    fSumw(),
    fCorr1(),
    fSumw1(),
    fCorr2(),
    fSumw2(),
    fcmNum(),
    fcmDen()
{
    Initialize(nbinsx,xlow,xhigh);
};
void AliPtPtContainer::Initialize(int nbinsx, double* xbins) {
    if(fCMTermList) delete fCMTermList;
    fCMTermList = new TList();
    fCMTermList->SetOwner(kTRUE);
    if(fCorrList) delete fCorrList;
    fCorrList = new TList();
    fCorrList->SetOwner(kTRUE);
    for(int m=0;m<mpar;++m) fCorrList->Add(new AliProfileBS(Form("corr_%ipar",m+1),this->GetTitle(),nbinsx,xbins));
    for(int m=0;m<4;++m) {
        for(int i=0;i<=m;++i)
            fCMTermList->Add(new AliProfileBS(Form("cm%i_Mpt%i",m+1,i),this->GetTitle(),nbinsx,xbins));
    }
    printf("Container %s initialized with m = %i\n",this->GetName(),mpar);
};
void AliPtPtContainer::Initialize(int nbinsx, double xlow, double xhigh) {
    if(fCMTermList) delete fCMTermList;
    fCMTermList = new TList();
    fCMTermList->SetOwner(kTRUE);
    if(fCorrList) delete fCorrList;
    fCorrList = new TList();
    fCorrList->SetOwner(kTRUE);
    for(int m=0;m<mpar;++m) fCorrList->Add(new AliProfileBS(Form("corr_%ipar",m+1),this->GetTitle(),nbinsx,xlow,xhigh));
    for(int m=0;m<4;++m) {
        for(int i=0;i<=m;++i)
            fCMTermList->Add(new AliProfileBS(Form("cm%i_Mpt%i",m+1,i),this->GetTitle(),nbinsx,xlow,xhigh));
    }
    printf("Container %s initialized with m = %i\n",this->GetName(),mpar);
};
void AliPtPtContainer::InitializeSubsamples(const int &nsub) {
    if(!fCMTermList || !fCorrList) return;
    for(int i=0; i<fCMTermList->GetEntries();++i)
            ((AliProfileBS*)fCMTermList->At(i))->InitializeSubsamples(nsub);
    for(int i=0; i<fCorrList->GetEntries();++i)
        ((AliProfileBS*)fCorrList->At(i))->InitializeSubsamples(nsub);
    if(fSubList){
        for(int i=0; i<fSubList->GetEntries();++i)
            ((AliProfileBS*)fSubList->At(i))->InitializeSubsamples(nsub);
        for(int i=0; i<fSubCMList->GetEntries();++i)
            ((AliProfileBS*)fSubCMList->At(i))->InitializeSubsamples(nsub);
    }
    return;
};
void AliPtPtContainer::InitializeSubevent(int nbinsx, double* xbins){
    if(fSubList) delete fSubList;
    fSubList = new TList();
    fSubList->SetOwner(kTRUE);
    for(int sub = 0; sub<2; ++sub){
        for(int m=0;m<mpar;++m) {
            fSubList->Add(new AliProfileBS(Form("corr_sub%i_%ipar",sub+1,m+1),this->GetTitle(),nbinsx,xbins));
        }
    }
    for(int m=2;m<=mpar;++m) {
        for(int k=0; k < m-1; ++k){
            fSubList->Add(new AliProfileBS(Form("corr_%isub1_%isub2_%ipar",m-k-1,k+1,m),this->GetTitle(),nbinsx,xbins));
        }
    }
    if(fSubCMList) delete fSubCMList;
    fSubCMList = new TList();
    fSubCMList->SetOwner(kTRUE);
    for(int sub = 0;sub <2; ++sub){
        for(int m=0;m<4;++m) {
            for(int i=0;i<=m;++i){
                fSubCMList->Add(new AliProfileBS(Form("cm%i_sub%i_Mpt%i",m+1,sub+1,i),this->GetTitle(),nbinsx,xbins));
            }
        }
    }
    for(int m = 2;m<=4;++m) {
        for(int first = 1; first < m; ++first) {
            for(int second = first; second < m; ++second) {
                if(first > second) continue;
                int fourth = m-second;
                for(int third = 1; third < m; ++third) {
                    if(third > fourth) continue;
                    fSubCMList->Add(new AliProfileBS(Form("cm%i_%i%isub1_%i%isub2",m,first,second,third,fourth),this->GetTitle(),nbinsx,xbins));
                }
            }
        }
    }

    printf("Subevents initialised for container %s with m = %i\n",this->GetName(),mpar);
}
void AliPtPtContainer::InitializeSubevent(int nbinsx, double xlow, double xhigh){
    if(fSubList) delete fSubList;
    fSubList = new TList();
    fSubList->SetOwner(kTRUE);
    for(int sub = 0; sub<2; ++sub){
        for(int m=0;m<mpar;++m) {
            fSubList->Add(new AliProfileBS(Form("corr_sub%i_%ipar",sub+1,m+1),this->GetTitle(),nbinsx,xlow,xhigh));
        }
    }
    for(int m=2;m<=mpar;++m) {
        for(int k=0; k < m-1; ++k){
            fSubList->Add(new AliProfileBS(Form("corr_%isub1_%isub2_%ipar",m-k-1,k+1,m),this->GetTitle(),nbinsx,xlow,xhigh));
        }
    }
    if(fSubCMList) delete fSubCMList;
    fSubCMList = new TList();
    fSubCMList->SetOwner(kTRUE);
    for(int sub = 0;sub <2; ++sub){
        for(int m=0;m<4;++m) {
            for(int i=0;i<=m;++i){
                fSubCMList->Add(new AliProfileBS(Form("cm%i_sub%i_Mpt%i",m+1,sub+1,i),this->GetTitle(),nbinsx,xlow,xhigh));
            }
        }
    }
    for(int m = 2;m<=4;++m) {
        for(int first = 1; first < m; ++first) {
            for(int second = first; second < m; ++second) {
                if(first > second) continue;
                int fourth = m-second;
                for(int third = 1; third < m; ++third) {
                    if(third > fourth) continue;
                    fSubCMList->Add(new AliProfileBS(Form("cm%i_%i%isub1_%i%isub2",m,first,second,third,fourth),this->GetTitle(),nbinsx,xlow,xhigh));
                }
            }
        }
    }

    printf("Subevents initialised for container %s with m = %i\n",this->GetName(),mpar);
}
vector<double> AliPtPtContainer::getEventCorrelation(int mOrder) {
  vector<double> outvec = {fCorr[mOrder],fSumw[mOrder]};
  return outvec;
}
void AliPtPtContainer::CalculateCorrelations(const vector<vector<double>> &inarr) {
  fCorr.clear(); fCorr.resize(mpar+1,0); fCorr[0] = 1.0;
  fSumw.clear(); fSumw.resize(mpar+1,0); fSumw[0] = 1.0;
  double sumNum = 0;
  double sumDenum = 0;
  vector<double> valNum;
  vector<double> valDenum;
  for(int m(1); m<=mpar; ++m)
  {
    for(int k(1);k<=m;++k)
    {
      valNum.push_back(fSign[k-1]*fCorr[m-k]*(fFactorial[m-1]/fFactorial[m-k])*inarr[k][k]);
      valDenum.push_back(fSign[k-1]*fSumw[m-k]*(fFactorial[m-1]/fFactorial[m-k])*inarr[k][0]);
    }
    sumNum = OrderedAddition(valNum);
    sumDenum = OrderedAddition(valDenum);
    valNum.clear();
    valDenum.clear();

    fCorr[m] = sumNum;
    fSumw[m] = sumDenum;
  }
  return;
}
void AliPtPtContainer::CalculateSubeventCorrelations(const vector<vector<double>> &insub1,const vector<vector<double>> &insub2 ) {
  fCorr1.clear(); fCorr1.resize(mpar+1,0); fCorr1[0] = 1.0;
  fSumw1.clear(); fSumw1.resize(mpar+1,0); fSumw1[0] = 1.0;
  fCorr2.clear(); fCorr2.resize(mpar+1,0); fCorr2[0] = 1.0;
  fSumw2.clear(); fSumw2.resize(mpar+1,0); fSumw2[0] = 1.0;
  double sumNum1 = 0;
  double sumDenum1 = 0;
  vector<double> valNum1;
  vector<double> valDenum1;
  double sumNum2 = 0;
  double sumDenum2 = 0;
  vector<double> valNum2;
  vector<double> valDenum2;
  for(int m(1); m<=mpar; ++m)
  {
    for(int k(1);k<=m;++k)
    {
      //correlations in subevent 1
      valNum1.push_back(fSign[k-1]*fCorr1[m-k]*(fFactorial[m-1]/fFactorial[m-k])*insub1[k][k]);
      valDenum1.push_back(fSign[k-1]*fSumw1[m-k]*(fFactorial[m-1]/fFactorial[m-k])*insub1[k][0]);
      //correlations in subevent 2
      valNum2.push_back(fSign[k-1]*fCorr2[m-k]*(fFactorial[m-1]/fFactorial[m-k])*insub2[k][k]);
      valDenum2.push_back(fSign[k-1]*fSumw2[m-k]*(fFactorial[m-1]/fFactorial[m-k])*insub2[k][0]);
    }
    sumNum1 = OrderedAddition(valNum1);
    sumDenum1 = OrderedAddition(valDenum1);
    sumNum2 = OrderedAddition(valNum2);
    sumDenum2 = OrderedAddition(valDenum2);
    valNum1.clear();
    valDenum1.clear();
    valNum2.clear();
    valDenum2.clear();
    fCorr1[m] = sumNum1;
    fSumw1[m] = sumDenum1;
    fCorr2[m] = sumNum2;
    fSumw2[m] = sumDenum2;
  }
  return;
}

void AliPtPtContainer::FillProfiles(const double &centmult, const double &rn) {
    for(int m=1;m<=mpar;++m)
    {
        if(fSumw[m]!=0) ((AliProfileBS*)fCorrList->At(m-1))->FillProfile(centmult,fCorr[m]/fSumw[m],(fEventWeight==PtPtSpace::kUnity)?1.0:fSumw[m],rn);
    }
    return;
}
void AliPtPtContainer::FillSubeventProfiles(const double &centmult, const double &rn) {
    //Fill the correlations within subevents, requires that the CalculateSubeventCorrelations have been called with the correct input vectors right before
    for(int m=1;m<=mpar;++m)
    {
        if(fSumw1[m]!=0) ((AliProfileBS*)fSubList->At(m-1))->FillProfile(centmult,fCorr1[m]/fSumw1[m],(fEventWeight==PtPtSpace::kUnity)?1.0:fSumw1[m],rn);
        //subevent 2 profiles offset by mpar positions
        if(fSumw2[m]!=0) ((AliProfileBS*)fSubList->At(mpar+m-1))->FillProfile(centmult,fCorr2[m]/fSumw2[m],(fEventWeight==PtPtSpace::kUnity)?1.0:fSumw2[m],rn);
    }

    //Fill the cross-subevent correlations
    for(int m=2;m<=mpar;++m) {
        for(int k=0; k < m-1; ++k){
            if(fSumw1[m-k-1]!=0 && fSumw2[k+1]!=0) ((AliProfileBS*)fSubList->FindObject(Form("corr_%isub1_%isub2_%ipar",m-k-1,k+1,m)))->FillProfile(centmult,fCorr1[m-k-1]/fSumw1[m-k-1]*fCorr2[k+1]/fSumw2[k+1],(fEventWeight==PtPtSpace::kUnity)?1.0:fSumw1[m-k-1]*fSumw2[k+1],rn);
        }
    }
    return;
}
void AliPtPtContainer::FillCMProfiles(const vector<vector<double>> &inarr, const double &centmult, const double &rn) {
  fcmNum.clear();
  fcmDen.clear();
  if (inarr[0][0] == 0)
    return;
  // 0th order correlation
  fcmDen.push_back(1.);
  fcmNum.push_back(1.);

  fcmDen.push_back(inarr[1][0]);
  fcmDen.push_back(inarr[1][0] * inarr[1][0] - inarr[2][0]);
  fcmDen.push_back(inarr[1][0] * inarr[1][0] * inarr[1][0] - 3 * inarr[2][0] * inarr[1][0] + 2 * inarr[3][0]);
  fcmDen.push_back(inarr[1][0] * inarr[1][0] * inarr[1][0] * inarr[1][0] - 6 * inarr[2][0] * inarr[1][0] * inarr[1][0] + 8 * inarr[1][0] * inarr[3][0] + 3 * inarr[2][0] * inarr[2][0] - 6 * inarr[4][0]);
  if (mpar < 1 || fcmDen[1] == 0)
    return;
  fcmNum.push_back(inarr[1][1] / fcmDen[1]);
  dynamic_cast<AliProfileBS*>(fCMTermList->At(0))->FillProfile(centmult, fcmNum[1], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[1], rn);
  if (mpar < 2 || inarr[2][0] == 0 || fcmDen[2] == 0)
    return;
  fcmNum.push_back(1 / fcmDen[2] * (inarr[1][1] * inarr[1][1] - inarr[2][2]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(1))->FillProfile(centmult, fcmNum[2], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[2], rn);
  fcmNum.push_back(1 / fcmDen[2] * (inarr[1][0] * inarr[1][1] - inarr[2][1]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(2))->FillProfile(centmult, fcmNum[3], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[2], rn);
  if (mpar < 3 || inarr[3][0] == 0 || fcmDen[3] == 0)
    return;
  fcmNum.push_back(1 / fcmDen[3] * (inarr[1][1] * inarr[1][1] * inarr[1][1] - 3 * inarr[2][2] * inarr[1][1] + 2 * inarr[3][3]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(3))->FillProfile(centmult, fcmNum[4], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[3], rn);
  fcmNum.push_back(1 / fcmDen[3] * (inarr[1][1] * inarr[1][1] * inarr[1][0] - 2 * inarr[2][1] * inarr[1][1] + 2 * inarr[3][2] - inarr[2][2] * inarr[1][0]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(4))->FillProfile(centmult, fcmNum[5], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[3], rn);
  fcmNum.push_back(1 / fcmDen[3] * (inarr[1][1] * inarr[1][0] * inarr[1][0] - 2 * inarr[2][1] * inarr[1][0] + 2 * inarr[3][1] - inarr[1][1] * inarr[2][0]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(5))->FillProfile(centmult, fcmNum[6], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[3], rn);
  if (mpar < 4 || inarr[4][0] == 0 || fcmDen[4] == 0)
    return;
  fcmNum.push_back(1 / fcmDen[4] * (inarr[1][1] * inarr[1][1] * inarr[1][1] * inarr[1][1] - 6 * inarr[2][2] * inarr[1][1] * inarr[1][1] + 3 * inarr[2][2] * inarr[2][2] + 8 * inarr[3][3] * inarr[1][1] - 6 * inarr[4][4]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(6))->FillProfile(centmult, fcmNum[7], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[4], rn);
  fcmNum.push_back(1 / fcmDen[4] * (inarr[1][1] * inarr[1][1] * inarr[1][1] * inarr[1][0] - 3 * inarr[2][2] * inarr[1][1] * inarr[1][0] - 3 * inarr[1][1] * inarr[1][1] * inarr[2][1] + 3 * inarr[2][2] * inarr[2][1] + 2 *inarr[3][3] * inarr[1][0] + 6 * inarr[1][1] * inarr[3][2] - 6 * inarr[4][3]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(7))->FillProfile(centmult, fcmNum[8], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[4], rn);
  fcmNum.push_back(1 / fcmDen[4] * (inarr[1][1] * inarr[1][1] * inarr[1][0] * inarr[1][0] - inarr[2][2] * inarr[1][0] * inarr[1][0] - inarr[2][0] * inarr[1][1] * inarr[1][1] + inarr[2][0] * inarr[2][2] - 4 * inarr[2][1] * inarr[1][1] * inarr[1][0] + 4 * inarr[3][2] * inarr[1][0] + 4 * inarr[3][1] * inarr[1][1] + 2 * inarr[2][1] * inarr[2][1] - 6 * inarr[4][2]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(8))->FillProfile(centmult, fcmNum[9], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[4], rn);
  fcmNum.push_back(1 / fcmDen[4] * (inarr[1][1] * inarr[1][0] * inarr[1][0] * inarr[1][0] - 3 * inarr[2][1] * inarr[1][0] * inarr[1][0] - 3 * inarr[1][1] * inarr[2][0] * inarr[1][0] + 3 * inarr[2][1] * inarr[2][0] + 2 * inarr[1][1] * inarr[3][0] + 6 * inarr[3][1] * inarr[1][0] - 6 * inarr[4][1]));
  dynamic_cast<AliProfileBS*>(fCMTermList->At(9))->FillProfile(centmult, fcmNum[10], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen[4], rn);
  return;
}
void AliPtPtContainer::FillCMSubeventProfiles(const vector<vector<double>> &sub1, const vector<vector<double>> &sub2, const double &centmult, const double &rn) {
  vector<Double_t> fcmNum1;
  vector<Double_t> fcmDen1;
  vector<Double_t> fcmNum2;
  vector<Double_t> fcmDen2;
  if (mpar < 1)
    return;
  // 0th order correlation
  fcmDen1.push_back(1.);
  fcmNum1.push_back(1.);
  fcmDen2.push_back(1.);
  fcmNum2.push_back(1.);

  fcmDen1.push_back(sub1[1][0]);
  fcmDen1.push_back(sub1[1][0] * sub1[1][0] - sub1[2][0]);
  fcmDen1.push_back(sub1[1][0] * sub1[1][0] * sub1[1][0] - 3 * sub1[2][0] * sub1[1][0] + 2 * sub1[3][0]);
  fcmDen1.push_back(sub1[1][0] * sub1[1][0] * sub1[1][0] * sub1[1][0] - 6 * sub1[2][0] * sub1[1][0] * sub1[1][0] + 8 * sub1[1][0] * sub1[3][0] + 3 * sub1[2][0] * sub1[2][0] - 6 * sub1[4][0]);

  fcmDen2.push_back(sub2[1][0]);
  fcmDen2.push_back(sub2[1][0] * sub2[1][0] - sub2[2][0]);
  fcmDen2.push_back(sub2[1][0] * sub2[1][0] * sub2[1][0] - 3 * sub2[2][0] * sub2[1][0] + 2 * sub2[3][0]);
  fcmDen2.push_back(sub2[1][0] * sub2[1][0] * sub2[1][0] * sub2[1][0] - 6 * sub2[2][0] * sub2[1][0] * sub2[1][0] + 8 * sub2[1][0] * sub2[3][0] + 3 * sub2[2][0] * sub2[2][0] - 6 * sub2[4][0]);


  if(fcmDen1[1]!=0) {
    fcmNum1.push_back(sub1[1][1] / fcmDen1[1]);
    dynamic_cast<AliProfileBS*>(fSubCMList->At(0))->FillProfile(centmult, fcmNum1[1], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen1[1], rn);
  }
  if(fcmDen2[1]!=0){
    fcmNum2.push_back(sub2[1][1] / fcmDen2[1]);
    dynamic_cast<AliProfileBS*>(fSubCMList->At(10+0))->FillProfile(centmult, fcmNum2[1], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen2[1], rn);
  }

  if (mpar < 2)
    return;
  if(sub1[2][0] != 0 && fcmDen1[2] != 0){
    fcmNum1.push_back(1 / fcmDen1[2] * (sub1[1][1] * sub1[1][1] - sub1[2][2]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(1))->FillProfile(centmult, fcmNum1[2], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen1[2], rn);
    fcmNum1.push_back(1 / fcmDen1[2] * (sub1[1][0] * sub1[1][1] - sub1[2][1]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(2))->FillProfile(centmult, fcmNum1[3], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen1[2], rn);
  }
  if(sub2[2][0] != 0 && fcmDen2[2] != 0){
        fcmNum2.push_back(1 / fcmDen2[2] * (sub2[1][1] * sub2[1][1] - sub2[2][2]));
        dynamic_cast<AliProfileBS*>(fSubCMList->At(10+1))->FillProfile(centmult, fcmNum2[2], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen2[2], rn);
        fcmNum2.push_back(1 / fcmDen2[2] * (sub2[1][0] * sub2[1][1] - sub2[2][1]));
        dynamic_cast<AliProfileBS*>(fSubCMList->At(10+2))->FillProfile(centmult, fcmNum2[3], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen2[2], rn);
  }

  if (mpar < 3)
    return;
  if(sub1[3][0] != 0 && fcmDen1[3] != 0){
    fcmNum1.push_back(1 / fcmDen1[3] * (sub1[1][1] * sub1[1][1] * sub1[1][1] - 3 * sub1[2][2] * sub1[1][1] + 2 * sub1[3][3]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(3))->FillProfile(centmult, fcmNum1[4], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen1[3], rn);
    fcmNum1.push_back(1 / fcmDen1[3] * (sub1[1][1] * sub1[1][1] * sub1[1][0] - 2 * sub1[2][1] * sub1[1][1] + 2 * sub1[3][2] - sub1[2][2] * sub1[1][0]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(4))->FillProfile(centmult, fcmNum1[5], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen1[3], rn);
    fcmNum1.push_back(1 / fcmDen1[3] * (sub1[1][1] * sub1[1][0] * sub1[1][0] - 2 * sub1[2][1] * sub1[1][0] + 2 * sub1[3][1] - sub1[1][1] * sub1[2][0]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(5))->FillProfile(centmult, fcmNum1[6], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen1[3], rn);
  }
  if(sub2[3][0] != 0 && fcmDen2[3] != 0){
    fcmNum2.push_back(1 / fcmDen2[3] * (sub2[1][1] * sub2[1][1] * sub2[1][1] - 3 * sub2[2][2] * sub2[1][1] + 2 * sub2[3][3]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(10+3))->FillProfile(centmult, fcmNum2[4], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen2[3], rn);
    fcmNum2.push_back(1 / fcmDen2[3] * (sub2[1][1] * sub2[1][1] * sub2[1][0] - 2 * sub2[2][1] * sub2[1][1] + 2 * sub2[3][2] - sub2[2][2] * sub2[1][0]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(10+4))->FillProfile(centmult, fcmNum2[5], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen2[3], rn);
    fcmNum2.push_back(1 / fcmDen2[3] * (sub2[1][1] * sub2[1][0] * sub2[1][0] - 2 * sub2[2][1] * sub2[1][0] + 2 * sub2[3][1] - sub2[1][1] * sub2[2][0]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(10+5))->FillProfile(centmult, fcmNum2[6], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen2[3], rn);
  }

  if (mpar < 4)
    return;
  if(sub1[4][0] != 0 && fcmDen1[4] != 0){
    fcmNum1.push_back(1 / fcmDen1[4] * (sub1[1][1] * sub1[1][1] * sub1[1][1] * sub1[1][1] - 6 * sub1[2][2] * sub1[1][1] * sub1[1][1] + 3 * sub1[2][2] * sub1[2][2] + 8 * sub1[3][3] * sub1[1][1] - 6 * sub1[4][4]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(6))->FillProfile(centmult, fcmNum1[7], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen1[4], rn);
    fcmNum1.push_back(1 / fcmDen1[4] * (sub1[1][1] * sub1[1][1] * sub1[1][1] * sub1[1][0] - 3 * sub1[2][2] * sub1[1][1] * sub1[1][0] - 3 * sub1[1][1] * sub1[1][1] * sub1[2][1] + 3 * sub1[2][2] * sub1[2][1] + 6 * sub1[1][1] * sub1[3][2] - 6 * sub1[4][3]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(7))->FillProfile(centmult, fcmNum1[8], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen1[4], rn);
    fcmNum1.push_back(1 / fcmDen1[4] * (sub1[1][1] * sub1[1][1] * sub1[1][0] * sub1[1][0] - sub1[2][2] * sub1[1][0] * sub1[1][0] - sub1[2][0] * sub1[1][1] * sub1[1][1] + sub1[2][0] * sub1[2][2] - 4 * sub1[2][1] * sub1[1][1] * sub1[1][0] + 4 * sub1[3][2] * sub1[1][0] + 4 * sub1[3][1] * sub1[1][1] + 2 * sub1[2][1] * sub1[2][1] - 6 * sub1[4][2]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(8))->FillProfile(centmult, fcmNum1[9], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen1[4], rn);
    fcmNum1.push_back(1 / fcmDen1[4] * (sub1[1][1] * sub1[1][0] * sub1[1][0] * sub1[1][0] - 3 * sub1[2][1] * sub1[1][0] * sub1[1][0] - 3 * sub1[1][1] * sub1[2][0] * sub1[1][0] + 3 * sub1[2][1] * sub1[2][0] + 2 * sub1[1][1] * sub1[3][0] + 6 * sub1[3][1] * sub1[1][0] - 6 * sub1[4][1]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(9))->FillProfile(centmult, fcmNum1[10], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen1[4], rn);
  }
  if(sub2[4][0] != 0 && fcmDen2[4] != 0){
    fcmNum2.push_back(1 / fcmDen2[4] * (sub2[1][1] * sub2[1][1] * sub2[1][1] * sub2[1][1] - 6 * sub2[2][2] * sub2[1][1] * sub2[1][1] + 3 * sub2[2][2] * sub2[2][2] + 8 * sub2[3][3] * sub2[1][1] - 6 * sub2[4][4]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(10+6))->FillProfile(centmult, fcmNum2[7], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen2[4], rn);
    fcmNum2.push_back(1 / fcmDen2[4] * (sub2[1][1] * sub2[1][1] * sub2[1][1] * sub2[1][0] - 3 * sub2[2][2] * sub2[1][1] * sub2[1][0] - 3 * sub2[1][1] * sub2[1][1] * sub2[2][1] + 3 * sub2[2][2] * sub2[2][1] + 6 * sub2[1][1] * sub2[3][2] - 6 * sub2[4][3]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(10+7))->FillProfile(centmult, fcmNum2[8], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen2[4], rn);
    fcmNum2.push_back(1 / fcmDen2[4] * (sub2[1][1] * sub2[1][1] * sub2[1][0] * sub2[1][0] - sub2[2][2] * sub2[1][0] * sub2[1][0] - sub2[2][0] * sub2[1][1] * sub2[1][1] + sub2[2][0] * sub2[2][2] - 4 * sub2[2][1] * sub2[1][1] * sub2[1][0] + 4 * sub2[3][2] * sub2[1][0] + 4 * sub2[3][1] * sub2[1][1] + 2 * sub2[2][1] * sub2[2][1] - 6 * sub2[4][2]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(10+8))->FillProfile(centmult, fcmNum2[9], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen2[4], rn);
    fcmNum2.push_back(1 / fcmDen2[4] * (sub2[1][1] * sub2[1][0] * sub2[1][0] * sub2[1][0] - 3 * sub2[2][1] * sub2[1][0] * sub2[1][0] - 3 * sub2[1][1] * sub2[2][0] * sub2[1][0] + 3 * sub2[2][1] * sub2[2][0] + 2 * sub2[1][1] * sub2[3][0] + 6 * sub2[3][1] * sub2[1][0] - 6 * sub2[4][1]));
    dynamic_cast<AliProfileBS*>(fSubCMList->At(10+9))->FillProfile(centmult, fcmNum2[10], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen2[4], rn);
  }

  //Fill cross terms
    for(int m = 2;m<=4;++m) {
        for(int first = 1; first < m; ++first) {
            for(int second = first; second < m; ++second) {
                if(first > second) continue;
                int fourth = m-second;
                for(int third = 1; third < m; ++third) {
                    if(third > fourth) continue;
                    if(sub1[m][0]!=0 && sub2[m][0]!=0 && fcmDen1[m]*fcmDen2[m]!=0) ((AliProfileBS*)fSubCMList->FindObject(Form("cm%i_%i%isub1_%i%isub2",m,first,second,third,fourth)))->FillProfile(centmult,fcmNum1[second*(second-1)/2+second-first+1]*fcmNum2[fourth*(fourth-1)/2+fourth-third+1], (fEventWeight == kEventWeight::kUnity) ? 1.0 : fcmDen1[m]*fcmDen2[m], rn);
                    //printf("m = %i, first = %i, second = %i, index = %i\n",m, first,second,second*(second-1)/2+second-first+1);
                    //printf("third = %i, fourth = %i, index = %i\n",third,fourth,fourth*(fourth-1)/2+fourth-third+1);
                }
            }
        }
    }
  return;
}
double AliPtPtContainer::OrderedAddition(vector<double> vec) {
  double sum = 0;
  std::sort(vec.begin(), vec.end());

  for(int i = 0; i < vec.size(); i++)
  {
    sum += vec[i];
  }
  return sum;
}
void AliPtPtContainer::RebinMulti(Int_t nbins) {
    if(fCMTermList) for(Int_t i=0;i<fCMTermList->GetEntries();i++) ((AliProfileBS*)fCMTermList->At(i))->RebinMulti(nbins);
    if(fCorrList) for(Int_t i=0;i<fCorrList->GetEntries();i++) ((AliProfileBS*)fCorrList->At(i))->RebinMulti(nbins);
    if(fSubList) for(Int_t i=0;i<fSubList->GetEntries();i++) ((AliProfileBS*)fSubList->At(i))->RebinMulti(nbins);
    if(fSubCMList) for(Int_t i=0;i<fSubCMList->GetEntries();i++) ((AliProfileBS*)fSubCMList->At(i))->RebinMulti(nbins);
    return;
}
void AliPtPtContainer::RebinMulti(Int_t nbins, Double_t *binedges) {
    if(fCMTermList) for(Int_t i=0;i<fCMTermList->GetEntries();i++) ((AliProfileBS*)fCMTermList->At(i))->RebinMulti(nbins,binedges);
    if(fCorrList) for(Int_t i=0;i<fCorrList->GetEntries();i++) ((AliProfileBS*)fCorrList->At(i))->RebinMulti(nbins,binedges);
    if(fSubList) for(Int_t i=0;i<fSubList->GetEntries();i++) ((AliProfileBS*)fSubList->At(i))->RebinMulti(nbins,binedges);
    if(fSubCMList) for(Int_t i=0;i<fSubCMList->GetEntries();i++) ((AliProfileBS*)fSubCMList->At(i))->RebinMulti(nbins,binedges);
    return;
}
TH1* AliPtPtContainer::getCorrHist(int ind, int m, int sub) {
    if(!sub) return ((AliProfileBS*)fCorrList->FindObject(Form("corr_%ipar",m)))->getHist(ind);
    if(sub && fSubList->GetEntries() == 0) { printf("The subevent profiles have not been filled! Returning nullptr\n"); return nullptr; }
    return  ((AliProfileBS*)fSubList->FindObject(Form("corr_sub%i_%ipar",sub,m)))->getHist(ind);
}
TH1* AliPtPtContainer::getCentralMomentHist(int ind, int m) {
  if(!fCentralMomentList) CreateCentralMomentList();
  if(!fCentralMomentList) return 0;
  if(ind+1<fCentralMomentList->GetEntries()) return (TH1*)fCentralMomentList->FindObject(Form("cm%i_%i",m,ind));
  return 0;
}
void AliPtPtContainer::CreateCentralMomentList() {
    if(fCentralMomentList) delete fCentralMomentList;
    fCentralMomentList = new TList();
    fCentralMomentList->SetOwner();
    for(auto m(1); m<=4;++m){
        for(int i=-1;i<((AliProfileBS*)fCMTermList->At(0))->getNSubs();++i) {
            TH1* hMpt = ((AliProfileBS*)fCMTermList->At(0))->getHist(i);
            vector<TH1*> hTs;
            for(int j=0;j<m;++j) {
                ((AliProfileBS*)fCMTermList->FindObject(Form("cm%i_Mpt%i",m,j)))->SetErrorOption("g");
                hTs.push_back(((AliProfileBS*)fCMTermList->FindObject(Form("cm%i_Mpt%i",m,j)))->getHist(i));
            }
            CalculateCentralMomentHists(hTs,i,m,hMpt);
        }
    }
    return;
}
void AliPtPtContainer::CalculateCentralMomentHists(vector<TH1*> inh, int ind, int m, TH1* hMpt) {
    auto binomial = [&](const int n, const int m) { return factorial(n)/(factorial(m)*factorial(n-m)); };
    TH1* reth = (TH1*)inh[0]->Clone(Form("cm%i_%i",m,ind));
    for(auto i(1);i<m;++i){
        TH1* mptPow = raiseHistToPower(hMpt,i);
        int coeff = binomial(m,i);
        coeff *= pow(-1,i);
        inh[i]->Multiply(mptPow);
        inh[i]->Scale(coeff);
        reth->Add(inh[i]);
    }
    TH1* mptLast = raiseHistToPower(hMpt,m);
    reth->Add(mptLast,(m%2)?(-1):1);
    fCentralMomentList->Add(reth);
    return;
}
TH1* AliPtPtContainer::getCumulantHist(int ind, int m) {
    if(!fCumulantList) CreateCumulantList();
    if(!fCumulantList) return 0;
    if(ind+1<fCumulantList->GetEntries()) return (TH1*)fCumulantList->At((ind+1)*mpar+m-1);
}
void AliPtPtContainer::CreateCumulantList() {
    if(fCumulantList) delete fCumulantList;
    fCumulantList = new TList();
    fCumulantList->SetOwner();
    //((AliProfileBS*)fCorrList->At(0))->PresetWeights((AliProfileBS*)fCorrList->At(mpar-1));
    for(int i=-1;i<((AliProfileBS*)fCorrList->At(0))->getNSubs();++i) {
        vector<TH1*> hTs;
        for(int j=0;j<mpar;++j) {
            ((AliProfileBS*)fCorrList->FindObject(Form("corr_%ipar",j+1)))->SetErrorOption("g");
            hTs.push_back(((AliProfileBS*)fCorrList->FindObject(Form("corr_%ipar",j+1)))->getHist(i));
        }
        CalculateCumulantHists(hTs,i);
    }
    //((AliProfileBS*)fCorrList->At(0))->PresetWeights(0);
    return;
}
void AliPtPtContainer::CalculateCumulantHists(vector<TH1*> inh, int ind) {
    auto binomial = [&](const int n, const int m) { return factorial(n)/(factorial(m)*factorial(n-m)); };
    for(int m=1;m<=mpar;++m)
    {
        TH1* reth = (TH1*)inh[m-1]->Clone(Form("reth%i_%i",m,ind));
        //TH1* hWeights = (TH1*)inh[m-1]->Clone(Form("hWeights%i_%i",m,ind));
        for(int k=1;k<m;++k)
        {
            TH1* corrh = (TH1*)inh[m-k-1]->Clone(Form("hcorr%i%i_%i",m,k,ind));
            corrh->Multiply((TH1*)fCumulantList->At((ind+1)*mpar+k-1));
            corrh->Scale(binomial(m-1,k-1));
            reth->Add(corrh,-1);
            delete corrh;
        }
        //for(int i=1;i<=hWeights->GetNbinsX();++i) reth->SetBinError(i,hWeights->GetBinError(i));
        //delete hWeights;
        fCumulantList->Add(reth->Clone(Form("kappa%i_%i",m,ind)));
    }
    return;
}
Long64_t AliPtPtContainer::Merge(TCollection *collist) {
  Long64_t nmerged=0;
  TIter all_PTC(collist);
  AliPtPtContainer *l_PTC = 0;
  while ((l_PTC = ((AliPtPtContainer*) all_PTC()))) {
      TList *t_CMTerm = l_PTC->fCMTermList;
      TList* t_Corr = l_PTC->fCorrList;
      TList* t_Sub = l_PTC->fSubList;
      TList* t_SubCMTerm = l_PTC->fSubCMList;
      TList* t_Cum = l_PTC->fCumulantList;
      TList* t_CM = l_PTC->fCentralMomentList;
      TList* t_SubCum = l_PTC->fSubCumulantList;
      TList* t_SubCM = l_PTC->fSubCentralMomentList;
      if(t_CMTerm) {
        if(!fCMTermList) fCMTermList = (TList*)t_CMTerm->Clone();
        else MergeBSLists(fCMTermList,t_CMTerm);
        nmerged++;
      };
      if(t_Corr) {
          if(!fCorrList) fCorrList = (TList*)t_Corr->Clone();
          else MergeBSLists(fCorrList,t_Corr);
      };
    if(t_Sub) {
          if(!fSubList) fSubList = (TList*)t_Sub->Clone();
          else MergeBSLists(fSubList,t_Sub);
      };
    if(t_SubCMTerm) {
          if(!fSubCMList) fSubCMList = (TList*)t_SubCMTerm->Clone();
          else MergeBSLists(fSubCMList,t_SubCMTerm);
      };
      if(t_Cum) {
          if(!fCumulantList) fCumulantList = (TList*)t_Cum->Clone();
          else MergeBSLists(fCumulantList,t_Cum);
      };
      if(t_CM) {
          if(!fCentralMomentList) fCentralMomentList = (TList*)t_CM->Clone();
          else MergeBSLists(fCentralMomentList,t_CM);
      }
      if(t_SubCum) {
          if(!fSubCumulantList) fSubCumulantList = (TList*)t_SubCum->Clone();
          else MergeBSLists(fSubCumulantList,t_SubCum);
      };
      if(t_SubCM) {
          if(!fSubCentralMomentList) fSubCentralMomentList = (TList*)t_SubCM->Clone();
          else MergeBSLists(fSubCentralMomentList,t_SubCM);
      }
  }
  return nmerged;
}
void AliPtPtContainer::MergeBSLists(TList *source, TList *target) {
  if(source->GetEntries()!=target->GetEntries()) { printf("Number in lists to be merged are not the same, skipping...\n"); return; };
  for(Int_t i=0;i<source->GetEntries();i++) {
    AliProfileBS *l_obj = (AliProfileBS*)source->At(i);
    AliProfileBS *t_obj = (AliProfileBS*)target->At(i);
    l_obj->MergeBS(t_obj);
  };
}
TH1* AliPtPtContainer::raiseHistToPower(TH1* inh, double p)
{
    TH1D* reth = (TH1D*)inh->Clone("reth");
    reth->SetName(Form("power%.2f_%s",p,inh->GetName()));
    for(int i=1;i<=inh->GetNbinsX();i++) {
        if(inh->GetBinContent(i)>=0 || std::floor(p) == p)
        {
            reth->SetBinContent(i,pow(inh->GetBinContent(i),p));
            reth->SetBinError(i,p*pow(reth->GetBinContent(i),p-1)*inh->GetBinError(i));
        }
        else
        {
            reth->SetBinContent(i,-999);
            reth->SetBinError(i,0.000000001);
        }

    }
    return reth;
}
