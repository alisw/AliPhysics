///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoModelWeightGeneratorLednicky : the most advanced weight       //
// generator available. Supports a large number of different pair types  //
// and interaction types. Can calculate pair weights coming from         //
// quantum statistics, coulomb interation and strong interaction ot any  //
// combination of the three, as applicable.                              //
// This class is a wrapper for the fortran code provided by Richard      //
// Lednicky.                                                             //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

//#include "StHbtMaker/ThCorrFctn/AliFemtoModelWeightGeneratorLednicky.h"
#include "AliFemtoModelWeightGeneratorLednicky.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoPair.h"
//#include "StarCallf77.h"
//#include <strstream.h>
//#include <iomanip.h>
//#include <stream>
//#include <iomanip>
#include <sstream>

#ifdef SOLARIS
# ifndef false
typedef int bool;
#define false 0
#define true 1
# endif
#endif

#ifdef WIN32
# ifdef CERNLIB_MSSTDCALL
#  define F77_UCASE
#  define type_of_call _stdcall
#  ifndef CERNLIB_QXCAPT
#    define CERNLIB_QXCAPT
#  endif
# else
#  define F77_LCASE
#  ifndef CERNLIB_QXNO_SC
#    define CERNLIB_QXNO_SC
#  endif
# endif
# define type_of_call  _stdcall
# define DEFCHARD   const char* , const int
# define DEFCHARL
# define PASSCHARD(string) string, strlen(string)
# define PASSCHARL(string)
#else
# define DEFCHARD     const char*
# define DEFCHARL   , const int
# define PASSCHARD(string) string
# define PASSCHARL(string) , strlen(string)
#endif
#ifdef CERNLIB_QXCAPT
#  define F77_NAME(name,NAME) NAME
#else
#  if defined(CERNLIB_QXNO_SC)
#    define F77_NAME(name,NAME) name
#  else
#    define F77_NAME(name,NAME) name##_
#  endif
#endif
#ifndef type_of_call
# define type_of_call
#endif

// --- Prototype of the function used in the weight calculator
//     (in FsiWeightLedinicky.F)
//      SUBROUTINE FSIINI(I_ITEST,I_LL,I_NS,I_ICH,I_IQS,I_ISI,I_I3C)
#define fsiini F77_NAME(fsiini,FSIINI)
extern "C" {void type_of_call F77_NAME(fsiini,FSIINI)(const int &itest,const int &ill, const int &ins, const int &ich, const int &iqs, const int &isi,const int &i3c);}
#define llini F77_NAME(llini,LLINI)
extern "C" {void type_of_call F77_NAME(llini,LLINI)(const int &lll,const int &ns, const int &itest);}

#define fsinucl F77_NAME(fsinucl,FSINUCL)
extern "C" {void type_of_call  F77_NAME(fsinucl,FSINUCL)(const double &mn,const double &cn);}
#define fsimomentum F77_NAME(fsimomentum,FSIMOMENTUM)
extern "C" {void type_of_call F77_NAME(fsimomentum,FSIMOMENTUM)(double &p1,double &p2);}
#define fsiposition F77_NAME(fsiposition,FSIPOSITION)
extern "C" {void type_of_call F77_NAME(fsiposition,FSIPOSITION)(double &x1,double &x2);}
#define fsiw F77_NAME(fsiw,FSIW)
extern "C" {void type_of_call F77_NAME(fsiw,FSIW)(const int &i,double &weif,double &wei,double &wein);}
#define ltran12 F77_NAME(ltran12,LTRAN12)
extern "C" {void type_of_call ltran12_();}

//K+K- model type
#define setkpkmmodel F77_NAME(setkpkmmodel,SETKPKMMODEL)
extern "C" {void type_of_call F77_NAME(setkpkmmodel,SETKPKMMODEL)(const int &i_model,const int &i_PhiOffOn);}

// Test function for Lambda potential
//#define printlam F77_NAME(printlam,PRINTLAM)
//extern "C" {void type_of_call printlam_();}
//there is not PRINTLAM in *.F file

// --- Additional prototyping of some CERN functions (in FsiTool.F)
typedef float   REAL;
typedef struct { REAL re; REAL im; } COMPLEX;
#define cgamma F77_NAME(cgamma,CGAMMA)
extern "C" {COMPLEX type_of_call cgamma_(COMPLEX*);}

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoModelWeightGeneratorLednicky);
  /// \endcond
#endif

AliFemtoModelWeightGeneratorLednicky::AliFemtoModelWeightGeneratorLednicky()
  : AliFemtoModelWeightGenerator()
  , fWei(0)
  , fWein(0)
  , fWeif(0)
  , fWeightDen(0)
  , fItest(0)
  , fIch(1)
  , fIqs(1)
  , fIsi(1)
  , fI3c(0)
  , fNuclMass(1.)
  , fNuclCharge(0.)
  , fSphereApp(false)
  , fT0App(false)
  , fLL(1)
  , fNuclChargeSign(1)
  , fSwap(0)
  , fLLMax(30)
  , fNS(4)
  , fLLName({"all-pairs",
             "neutron neutron",
             "proton proton",
             "neutron proton",
             "alpha alpha",
             "pi+ pi-",
             "pi0 pi0",
             "pi+ pi+",
             "neutron deuteron",
             "proton deuteron",
             "pi+ K-",
             "pi+ K+",
             "pi+ proton",
             "pi- proton",
             "K+ K-",
             "K+ K+",
             "K+ proton",
             "K- proton",
             "deuteron deuteron",
             "deuton alpha",
             "triton triton",
             "triton alpha",
             "K0 K0",
             "K0 K0b",
             "deuteron triton",
             "proton triton",
             "proton alpha",
             "proton lambda",
             "neutron lambda",
             "Lambda lambda",
             "Proton Anti-proton"})
  , fNumProcessPair(nullptr)
  , fNumbNonId(0)
  , fKpKmModel(14)
  , fPhi_OffOn(1)
{
  // default constructor
  fNumProcessPair = new int[fLLMax+1];
  for (int i=1;i<=fLLMax;i++) {
    fNumProcessPair[i] = 0;
  }

  SetPid(211,211);
  FsiInit();
  FsiNucl();
}
//______________________
AliFemtoModelWeightGeneratorLednicky
  ::AliFemtoModelWeightGeneratorLednicky(const AliFemtoModelWeightGeneratorLednicky &aWeight)
  : AliFemtoModelWeightGenerator(aWeight)
  , fWei(aWeight.fWei)
  , fWein(aWeight.fWein)
  , fWeif(aWeight.fWeif)
  , fWeightDen(aWeight.fWeightDen)
  , fItest(aWeight.fItest)
  , fIch(aWeight.fIch)
  , fIqs(aWeight.fIqs)
  , fIsi(aWeight.fIsi)
  , fI3c(aWeight.fI3c)
  , fNuclMass(aWeight.fNuclMass)
  , fNuclCharge(aWeight.fNuclCharge)
  , fSphereApp(aWeight.fSphereApp)
  , fT0App(aWeight.fT0App)
  , fLL(aWeight.fLL)
  , fNuclChargeSign(aWeight.fNuclChargeSign)
  , fSwap(aWeight.fSwap)
  , fLLMax(30)
  , fNS(4)
  , fLLName(aWeight.fLLName.begin(), aWeight.fLLName.end())
  , fNumProcessPair(nullptr)
  , fNumbNonId(aWeight.fNumbNonId)
  , fKpKmModel(aWeight.fKpKmModel)
  , fPhi_OffOn(aWeight.fPhi_OffOn)
{
  fNumProcessPair = new int[fLLMax+1];
  for (int i=1;i<=fLLMax;i++) {
    fNumProcessPair[i] = 0;
  }

  FsiInit();
  FsiNucl();
}

AliFemtoModelWeightGeneratorLednicky&
AliFemtoModelWeightGeneratorLednicky::operator=(const AliFemtoModelWeightGeneratorLednicky& aWeight)
{
  // assignment operator
  if (this == &aWeight) {
    return *this;
  }

  AliFemtoModelWeightGenerator::operator=(aWeight);

  fWei = aWeight.fWei;
  fWein = aWeight.fWein;
  fWeif = aWeight.fWeif;
  fWeightDen = aWeight.fWeightDen;

  fItest = aWeight.fItest;
  fIch = aWeight.fIch;
  fIqs = aWeight.fIqs;
  fIsi = aWeight.fIsi;
  fI3c = aWeight.fI3c;
  fNuclMass = aWeight.fNuclMass;
  fNuclCharge = aWeight.fNuclCharge;
  fSphereApp = aWeight.fSphereApp;
  fT0App = aWeight.fT0App;
  fLL = aWeight.fLL;
  fNS = aWeight.fNS;
  fNuclChargeSign = aWeight.fNuclChargeSign;
  fSwap = aWeight.fSwap;
  fNumbNonId = aWeight.fNumbNonId;
  fLLName = aWeight.fLLName;

  if (fNumProcessPair) {
    delete [] fNumProcessPair;
  }
  fNumProcessPair = new int[fLLMax+1];
  fKpKmModel = aWeight.fKpKmModel;
  fPhi_OffOn = aWeight.fPhi_OffOn;

  for (int i=1;i<=fLLMax;i++) {
    fNumProcessPair[i] = 0;
  }

  FsiInit();
  FsiNucl();

  return *this;
}


double AliFemtoModelWeightGeneratorLednicky::GenerateWeight(AliFemtoPair* aPair)
{
  // Get hidden information pointers
  //AliFemtoModelHiddenInfo *inf1 = (AliFemtoModelHiddenInfo *) aPair->Track1()->HiddenInfo();
  //AliFemtoModelHiddenInfo *inf2 = (AliFemtoModelHiddenInfo *) aPair->Track2()->HiddenInfo();
  {
    double cached_weight = aPair->LookupFemtoWeightCache(this);
    if (!std::isnan(cached_weight)) {
      return cached_weight;
    }
  }

  const AliFemtoTrack &track1 = *aPair->Track1()->Track(),
                      &track2 = *aPair->Track2()->Track();

  auto &hinfo1 = *static_cast<AliFemtoModelHiddenInfo*>(track1.GetHiddenInfo()),
       &hinfo2 = *static_cast<AliFemtoModelHiddenInfo*>(track2.GetHiddenInfo());


  const auto &true_p1 = *hinfo1.GetTrueMomentum(),
             &true_p2 = *hinfo2.GetTrueMomentum();

  // Calculate pair variables
  // Double_t tPx = inf1->GetTrueMomentum()->x()+inf2->GetTrueMomentum()->x();
  Double_t tPx = true_p1.x() + true_p2.x();
  Double_t tPy = true_p1.y() + true_p2.y();
  Double_t tPz = true_p1.z() + true_p2.z();

  Double_t tM1 = hinfo1.GetMass();
  Double_t tM2 = hinfo2.GetMass();

  Double_t tE1 = sqrt(tM1*tM1 + true_p1.Mag2());
  Double_t tE2 = sqrt(tM2*tM2 + true_p2.Mag2());


 // if (tPx==0 && tPy==0 && tPz==0 ) {cout<<" zero true momentum "<<endl; return 0;}


  Double_t tE  = tE1 + tE2;
  Double_t tPt = tPx*tPx + tPy*tPy;
  Double_t tMt = tE*tE - tPz*tPz;//mCVK;
  Double_t tM  = sqrt(tMt - tPt);
  tMt = sqrt(tMt);
  tPt = sqrt(tPt);

  if (tMt==0 || tE==0 || tM==0 || tPt==0 ) {
    std::cout << " weight generator zero tPt || tMt || tM || tPt"
              << tM1 << " " << tM2 << "\n";
    return 0.0;
  }


  Double_t tBetat = tPt/tMt;

  // Boost to LCMS
  Double_t tBeta = tPz/tE;
  Double_t tGamma = tE/tMt;

  Double_t pX = true_p1.x();
  Double_t pY = true_p1.y();
  Double_t pZ = true_p1.z();

  fKStarLong = tGamma * (pZ - tBeta * tE1);
  Double_t tE1L = tGamma * (tE1  - tBeta * pZ);

  // Rotate in transverse plane
  fKStarOut  = ( pX*tPx + pY*tPy)/tPt;
  fKStarSide = (-pX*tPy + pY*tPx)/tPt;

  // Boost to pair cms
  fKStarOut = tMt/tM * (fKStarOut - tPt/tMt * tE1L);

  fKStar = ::sqrt(fKStarOut * fKStarOut
                  + fKStarSide * fKStarSide
                  + fKStarLong * fKStarLong);
  tBetat = tPt/tMt;

  const auto &epoint1 = *hinfo1.GetEmissionPoint(),
             &epoint2 = *hinfo2.GetEmissionPoint();

  Double_t tDX = epoint1.x() - epoint2.x();
  Double_t tDY = epoint1.y() - epoint2.y();
  Double_t tRLong = epoint1.z() - epoint2.z();
  Double_t tDTime = epoint1.t() - epoint2.t();

  Double_t tROut = (tDX*tPx + tDY*tPy)/tPt;
  Double_t tRSide = (-tDX*tPy + tDY*tPx)/tPt;


//cout<<"Weight generator"<<" tDX "<<tDX<<" tDY "<<tDY<<"tRLong "<<tRLong<<endl;

 /*
  cout << "Got points 1 " << epoint1.x() << "  "
                          << epoint1.y() << "  "
                          << epoint1.z() << "  "
                          << epoint1.t() << "\n";

  cout << "Got points 2 " << epoint2.x() << "  "
                          << epoint2.y() << "  "
                          << epoint2.z() << "  "
                          << epoint2.t() << "\n";
 */

  fRStarSide = tRSide;

  fRStarLong = tGamma * (tRLong - tBeta * tDTime);
  Double_t tDTimePairLCMS = tGamma * (tDTime - tBeta * tRLong);

  tBeta = tPt/tMt;
  tGamma = tMt/tM;

  fRStarOut = tGamma * (tROut - tBeta * tDTimePairLCMS);
  fRStar = ::sqrt(fRStarOut * fRStarOut
                  + fRStarSide * fRStarSide
                  + fRStarLong * fRStarLong);

  //cout << "-- weights generator : Got out side " << fRStarOut << " " << fRStarSide << endl;

  const int pdg1 = hinfo1.GetPDGPid(),
            pdg2 = hinfo2.GetPDGPid();

  // Check bad PID
  if (!SetPid(pdg1, pdg2)) {
    fWeightDen = 1.0;
    //    cout<<" bad PID weight generator pdg1 "<<hinfo1.GetPDGPid()<<" pdg2 " << hinfo2.GetPDGPid()<<endl;
    return 1; //non-correlated
  }

  // cout<<" good PID weight generator pdg1 "<<hinfo1.GetPDGPid()<<" pdg2 "<<hinfo2.GetPDGPid()<<endl;

  if (true_p1 == true_p2) {
    fWeightDen = 0.;
    return 0;
  }

  double p1[] = {true_p1.x(), true_p1.y(), true_p1.z()},
         p2[] = {true_p2.x(), true_p2.y(), true_p2.z()};

  if (fSwap) {
    fsimomentum(*p2,*p1);
  } else {
    fsimomentum(*p1,*p2);
  }

  if (epoint1 == epoint2) {
    fWeightDen=0.;
    return 0;
  }

//    if(pdg1==!211||pdg2!=211)cout << "Weight pdg1 pdg2 = " << pdg1<<" "<<pdg2<< endl;
//     cout << "LL:in GetWeight = " << mLL << endl;

  double x1[] = {epoint1.x(), epoint1.y(), epoint1.z(), epoint1.t()},
         x2[] = {epoint2.x(), epoint2.y(), epoint2.z(), epoint2.t()};

  if (fSwap) {
    fsiposition(*x2,*x1);
  } else {
    fsiposition(*x1,*x2);
  }

  FsiSetLL();
  ltran12();
  fsiw(1, fWeif, fWei, fWein);

  //  cout<<" fWeif "<<fWeif<<" fWei "<<fWei<<" fWein "<<fWein<<endl;

  if (fI3c == 0) {
    aPair->AddWeightToCache(this, fWein);
    return fWein;
  }

  fWeightDen = fWeif;

  aPair->AddWeightToCache(this, fWei);
  return fWei;
}


AliFemtoString
AliFemtoModelWeightGeneratorLednicky::Report()
{
  // create report
  ostringstream tStr;
  tStr << "Lednicky afterburner calculation for  Correlation -  Report" << endl;
  tStr << "    Setting : Quantum : " << ((fIqs) ? "On" : "Off");
  tStr << " - Coulbomb : " << ((fIch) ? "On" : "Off") ;
  tStr << " - Strong : " << ((fIsi) ? "On" : "Off");
  tStr << endl;
  tStr << "              3-Body : " << ((fI3c) ? "On"  : "Off") ;
  if (fI3c) tStr << " Mass=" <<  fNuclMass << " - Charge= " << fNuclCharge ;
  tStr << endl;
  tStr << "    " << fNumProcessPair[0] << " Pairs have been Processed :" << endl;
  for (int i=1;i<=fLLMax;i++) {
    if (fNumProcessPair[i]) {
      tStr << "         " << fNumProcessPair[i] << " " << fLLName[i] << endl;
    }
  }
  if (fNumbNonId) {
    tStr << "         "<< fNumbNonId << " Non Identified" << endl;
  }
  AliFemtoString returnThis = tStr.str();
  return returnThis;
}
/*
void AliFemtoModelWeightGeneratorLednicky::FsiInit()
{
  // Initialize weight generation module
   cout << "*******************AliFemtoModelWeightGeneratorLednicky check FsiInit ************" << endl;
   cout <<"mItest dans FsiInit() = " << fItest << endl;
   cout <<"mIch dans FsiInit() = " << fIch << endl;
   cout <<"mIqs dans FsiInit() = " << fIqs << endl;
   cout <<"mIsi dans FsiInit() = " << fIsi << endl;
   cout <<"mI3c dans FsiInit() = " << fI3c << endl;

  fsiini(fItest,fIch,fIqs,fIsi,fI3c);
}
*/
void AliFemtoModelWeightGeneratorLednicky::FsiInit()
{
  // Initialize weight generation module
   cout << "*******************AliFemtoModelWeightGeneratorLednicky check FsiInit ************" << endl;
   /*
C-   LL       1  2  3  4  5   6   7   8  9 10  11  12  13  14 15 16 17
C-   part. 1: n  p  n  a  pi+ pi0 pi+ n  p pi+ pi+ pi+ pi- K+ K+ K+ K-
C-   part. 2: n  p  p  a  pi-  pi0  pi+ d d  K-   K+  p     p   K-  K+ p  p
C   NS=1 y/n: +  +  +  +  +   -   -   -  -  -   -   -   -  -  -  -  -
C----------------------------------------------------------------------
C-   LL       18 19 20 21 22 23  24 25 26 27 28 29 30 31  32  33  34
C-   part. 1: d  d  t  t  K0 K0  d  p  p  p  n  /\ p  pi+ pi- p   p
C-   part. 2: d  a  t  a  K0 K0b t  t  a  /\ /\ /\ pb Xi- Xi- Om- Omb
C   NS=1 y/n: -  -  -  -  -  -   -  -  -  +  +  +  -  -   -   -   -
C----------------------------------------------------------------------
C-   LL       35 
C-   part. 1: K+ 
C-   part. 2: K0b 
C   NS=1 y/n: -  
   */
   if (fPairType == fgkPionPlusPionPlus) fLL = 8;
   if (fPairType == fgkPionPlusPionMinus ) fLL = 6;
   if (fPairType == fgkKaonPlusKaonPlus ) fLL = 15;
   if (fPairType == fgkKaonPlusKaonMinus ) fLL = 14;
   if (fPairType == fgkProtonProton ) fLL = 2;
   if (fPairType == fgkProtonAntiproton ) fLL = 30;
   if (fPairType == fgkPionPlusKaonPlus ) fLL = 11;
   if (fPairType == fgkPionPlusKaonMinus ) fLL = 10;
   if (fPairType == fgkPionPlusProton ) fLL = 12;
   if (fPairType == fgkPionPlusAntiproton ) fLL = 13;
   if (fPairType == fgkKaonPlusProton ) fLL = 16;
   if (fPairType == fgkKaonPlusAntiproton ) fLL = 17;

   cout<<"fPairType: "<<fPairType<<endl;
   cout <<"mItest dans FsiInit() = " << fItest << endl; //ok
   cout <<"mLL dans FsiInit() = " << fLL << endl; //ok
   cout <<"mNS dans FsiInit() = " << fNS << endl; //ok
   cout <<"mIch dans FsiInit() = " << fIch << endl; //ok
   cout <<"mIqs dans FsiInit() = " << fIqs << endl; //ok
   cout <<"mIsi dans FsiInit() = " << fIsi << endl;  //ok
   cout <<"mI3c dans FsiInit() = " << fI3c << endl; //ok


  fsiini(fItest,fLL,fNS,fIch,fIqs,fIsi,fI3c);
}

void AliFemtoModelWeightGeneratorLednicky::FsiSetKpKmModelType()
{
  // initialize K+K- model type
  cout<<"******************* AliFemtoModelWeightGeneratorLednicky check FsiInit initialize K+K- model type with FsiSetKpKmModelType(), type= "<<fKpKmModel<<" PhiOffON= "<<fPhi_OffOn<<" *************"<< endl;
   setkpkmmodel(fKpKmModel,fPhi_OffOn);
   cout<<"-----------------END FsiSetKpKmModelType-------"<<endl;
}

void AliFemtoModelWeightGeneratorLednicky::FsiNucl()
{
  // initialize weight generation taking into account the residual charge
//   cout << "*******************AliFemtoModelWeightGeneratorLednicky check FsiNucl ************" << endl;
//   cout <<"fNuclMass dans FsiNucl() = " << fNuclMass << endl;
//   cout <<"fNuclCharge dans FsiNucl() = " << fNuclCharge << endl;
//   cout <<"fNuclChargeSign dans FsiNucl() = " << fNuclChargeSign << endl;
  fsinucl(fNuclMass,fNuclCharge*fNuclChargeSign);
}

void AliFemtoModelWeightGeneratorLednicky::FsiSetLL()
{
  
  // set internal pair type for the module
  int tNS;
  if (fSphereApp||(fLL>5)) {
    if (fT0App) { tNS=4;}
    else {tNS=2;}
  } else { tNS=1;}
  //cout<<"*********************** AliFemtoModelWeightGeneratorLednicky::FsiSetLL() *********************"<<endl;
  if(fNS_4==4) tNS=4;//K+K- analisys
  //cout <<"fLL dans FsiSetLL() = "<< fLL << endl;
  //cout <<"tNS dans FsiSetLL() = "<< tNS << endl;
  //cout <<"fItest dans FsiSetLL() = "<< fItest << endl;

  //cout <<"fLL dans FsiSetLL() = "<< fLL << endl;
   //cout <<"tNS dans FsiSetLL() = "<< tNS << endl;
  // cout <<"fItest dans FsiSetLL() = "<< fItest << endl;
  llini(fLL,tNS,fItest);
 // cout<<" end of FsiSetLL"<<endl;
}

bool AliFemtoModelWeightGeneratorLednicky::SetPid(const int aPid1,const int aPid2)
{
  // set calculated system for basing on particles' pids
  static const int ksPi0Pid=111;
  static const int ksPionPid=211;
  static const int ksK0Pid=311;
  static const int ksKPid=321;
  static const int ksNeutPid=2112;
  static const int ksProtPid=2212;
  static const int ksLamPid=3122;
  //  static const int sLamLamPid=3122;

//    cout << "in SetPiD Setting PID to " << aPid1 << " " << aPid2 << endl;

  int tPidl,tPidh;
  int tChargeFactor=1;

  if (abs(aPid1)<abs(aPid2)) {
    if (aPid1<0) tChargeFactor=-1;
    tPidl=aPid1*tChargeFactor;
    tPidh=aPid2*tChargeFactor;
    fSwap=false;
  } else {
    if (aPid2<0) tChargeFactor=-1;
    tPidl=aPid2*tChargeFactor;
    tPidh=aPid1*tChargeFactor;
    fSwap=true;
  }
  switch (tPidl) {
  case ksPionPid:
    switch (tPidh) {
    case -ksPionPid:   fLL=5; tChargeFactor*=1 ;break;
    case ksPionPid:    fLL=7; tChargeFactor*=1 ;break;
    case -ksKPid:      fLL=10;tChargeFactor*=1 ;break;
    case ksKPid:       fLL=11;tChargeFactor*=1 ;break;
    case ksProtPid:    fLL=12;tChargeFactor*=1 ;break;
    case -ksProtPid:   fLL=13;tChargeFactor*=-1;break;
    default: fLL=0;
    }
    break;
  case ksProtPid:
    switch (tPidh) {
    case ksProtPid:    fLL=2; tChargeFactor*=1 ;break;
    case ksLamPid:     fLL=27;tChargeFactor*=1 ;break;
    case -ksProtPid:   fLL=30;tChargeFactor*=1 ;break;
    default: fLL=0;
    }
    break;
  case ksKPid:
    switch (tPidh) {
    case -ksKPid:      fLL=14;tChargeFactor*=1 ;break;
    case ksKPid:       fLL=15;tChargeFactor*=1 ;break;
    case ksProtPid:    fLL=16;tChargeFactor*=1 ;break;
    case -ksProtPid:   fLL=17;tChargeFactor*=-1 ;break;
    default: fLL=0;
    }
    break;
  case ksK0Pid:
    switch (tPidh) {
    case ksK0Pid:         fLL=22;tChargeFactor*=1 ;break;
    case -ksK0Pid:        fLL=23;tChargeFactor*=1 ;break;
    default: fLL=0;
    }
    break;
  case ksPi0Pid:
    switch (tPidh) {
    case ksPi0Pid:        fLL=6; tChargeFactor*=1 ;break;
    default: fLL=0;
    }
    break;
  case ksNeutPid:
    switch (tPidh) {
    case ksNeutPid:      fLL=1; tChargeFactor*=1 ;break;
    case ksProtPid:      fLL=3; tChargeFactor*=1 ;break;
    case ksLamPid:       fLL=28;tChargeFactor*=1 ;break;
    default: fLL=0;
    }
    break;                                             //Gael 21May02
  case ksLamPid:                                        //Gael 21May02
    switch (tPidh) {                                   //Gael 21May02
    case ksLamPid:       fLL=29;tChargeFactor*=1 ;break;//Gael 21May02
    default: fLL=0;                                    //Gael 21May02
    }                                                 //Gael 21May02
    break;                                             //Gael 21May02
  default:
    fLL=0;
  }
  if (tChargeFactor!=fNuclChargeSign) {
    fNuclChargeSign=tChargeFactor;
    FsiNucl();
  }
  (fNumProcessPair[0])++;
  if (fLL) {
    (fNumProcessPair[fLL])++;
    return true;
  } else {
    fNumbNonId++;
    return false;
  }
//   cout << "*******************AliFemtoModelWeightGeneratorLednicky check SetPid ************" << endl;
 //  cout << "fLL=="<< fLL << endl;
 //  cout << "fNuclCharge=="<< fNuclCharge << endl;

}

AliFemtoModelWeightGeneratorLednicky::~AliFemtoModelWeightGeneratorLednicky()
{
  if (fNumProcessPair) {
    delete [] fNumProcessPair;
  }
}


void
AliFemtoModelWeightGeneratorLednicky::SetPairType(Int_t aPairType)
{
  // set calculated system basing on the pair type
  fPairType = aPairType;
  if (fPairType == fgkPionPlusPionPlus) SetPid(211,211);
  if (fPairType == fgkPionPlusPionMinus ) SetPid(211, -211);
  if (fPairType == fgkKaonPlusKaonPlus ) SetPid(321, 321);
  if (fPairType == fgkKaonPlusKaonMinus ) SetPid(321, -321);
  if (fPairType == fgkProtonProton ) SetPid(2212, 2212);
  if (fPairType == fgkProtonAntiproton ) SetPid(2212, -2212);
  if (fPairType == fgkPionPlusKaonPlus ) SetPid(211, 321);
  if (fPairType == fgkPionPlusKaonMinus ) SetPid(211, -321);
  if (fPairType == fgkPionPlusProton ) SetPid(211, 2212);
  if (fPairType == fgkPionPlusAntiproton ) SetPid(211, -2212);
  if (fPairType == fgkKaonPlusProton ) SetPid(321, 2212);
  if (fPairType == fgkKaonPlusAntiproton ) SetPid(321, -2212);
}

//_____________________________________________
Int_t    AliFemtoModelWeightGeneratorLednicky::GetPairType() const
{
  // return pair type
  return fPairType;
}

//_____________________________________________
void     AliFemtoModelWeightGeneratorLednicky::SetPairTypeFromPair(AliFemtoPair *aPair)
{
  // set calculated system based on the hidden info in the pair
  AliFemtoModelHiddenInfo *inf1 = ( AliFemtoModelHiddenInfo *) aPair->Track1()->HiddenInfo();
  AliFemtoModelHiddenInfo *inf2 = ( AliFemtoModelHiddenInfo *) aPair->Track2()->HiddenInfo();

  const Int_t ktPid1 = inf1->GetPDGPid();
  const Int_t ktPid2 = inf2->GetPDGPid();

  fPairType = (((ktPid1 ==   211) && (ktPid2 ==   211)) ||
               ((ktPid1 ==  -211) && (ktPid2 ==  -211))) ? fgkPionPlusPionPlus

            : (((ktPid1 ==  -211) && (ktPid2 ==   211)) ||
               ((ktPid1 ==   211) && (ktPid2 ==  -211))) ? fgkPionPlusPionMinus

            : (((ktPid1 ==   321) && (ktPid2 ==   321)) ||
               ((ktPid1 ==  -321) && (ktPid2 ==  -321))) ? fgkKaonPlusKaonPlus

            : (((ktPid1 ==  -321) && (ktPid2 ==   321)) ||
               ((ktPid1 ==   321) && (ktPid2 ==  -321))) ? fgkKaonPlusKaonMinus

            : (((ktPid1 ==  2212) && (ktPid2 ==  2212)) ||
               ((ktPid1 == -2212) && (ktPid2 == -2212))) ? fgkProtonProton

            : (((ktPid1 == -2212) && (ktPid2 ==  2212)) ||
               ((ktPid1 ==  2212) && (ktPid2 == -2212))) ? fgkProtonAntiproton

            : (((ktPid1 ==   211) && (ktPid2 ==   321)) ||
               ((ktPid1 ==  -211) && (ktPid2 ==  -321))) ? fgkPionPlusKaonPlus

            : (((ktPid1 ==  -211) && (ktPid2 ==   321)) ||
               ((ktPid1 ==   211) && (ktPid2 ==  -321))) ? fgkPionPlusKaonMinus

            : (((ktPid1 ==   211) && (ktPid2 ==  2212)) ||
               ((ktPid1 ==  -211) && (ktPid2 == -2212))) ? fgkPionPlusProton

            : (((ktPid1 ==  -211) && (ktPid2 ==  2212)) ||
               ((ktPid1 ==   211) && (ktPid2 == -2212))) ? fgkPionPlusAntiproton

            : (((ktPid1 ==   321) && (ktPid2 ==  2212)) ||
               ((ktPid1 ==  -321) && (ktPid2 == -2212))) ? fgkKaonPlusProton

            : (((ktPid1 ==  -321) && (ktPid2 ==  2212)) ||
               ((ktPid1 ==   321) && (ktPid2 == -2212))) ? fgkKaonPlusAntiproton

            // no change
            : fPairType;

  SetPid(ktPid1, ktPid2);
}

//K+K- model type
void AliFemtoModelWeightGeneratorLednicky::SetKpKmModelType(const int aModelType, const int aPhi_OffOn)
{
  fKpKmModel = aModelType;
  fPhi_OffOn = aPhi_OffOn;
  fNS_4 = 4;
  FsiSetKpKmModelType();
}

void AliFemtoModelWeightGeneratorLednicky::SetNuclCharge(const double aNuclCharge)
  { fNuclCharge = aNuclCharge; FsiNucl(); }
void AliFemtoModelWeightGeneratorLednicky::SetNuclMass(const double aNuclMass)
  { fNuclMass = aNuclMass; FsiNucl(); }

void AliFemtoModelWeightGeneratorLednicky::SetSphere()
  { fSphereApp = true; }
void AliFemtoModelWeightGeneratorLednicky::SetSquare()
  { fSphereApp=false; }
void AliFemtoModelWeightGeneratorLednicky::SetT0ApproxOn()
  { fT0App = true; }
void AliFemtoModelWeightGeneratorLednicky::SetT0ApproxOff()
  { fT0App = false; }

void AliFemtoModelWeightGeneratorLednicky::SetNS(int mNS)
{
  fNS = mNS;
}

void AliFemtoModelWeightGeneratorLednicky::SetDefaultCalcPar()
{
  fItest = 1;
  fIqs = 1;
  fIsi = 1;
  fI3c = 0;
  fIch = 1;
  FsiInit();
  fSphereApp=false;
  fT0App=false;
}

void AliFemtoModelWeightGeneratorLednicky::SetCoulOn()    {fItest=1;fIch=1;FsiInit();}
void AliFemtoModelWeightGeneratorLednicky::SetCoulOff()   {fItest=1;fIch=0;FsiInit();}
void AliFemtoModelWeightGeneratorLednicky::SetQuantumOn() {fItest=1;fIqs=1;FsiInit();}
void AliFemtoModelWeightGeneratorLednicky::SetQuantumOff(){fItest=1;fIqs=0;FsiInit();}
void AliFemtoModelWeightGeneratorLednicky::SetStrongOn()  {fItest=1;fIsi=1;FsiInit();}
void AliFemtoModelWeightGeneratorLednicky::SetStrongOff() {fItest=1;fIsi=0;FsiInit();}
void AliFemtoModelWeightGeneratorLednicky::Set3BodyOn()   {fItest=1;fI3c=1;FsiInit();FsiNucl();}
void AliFemtoModelWeightGeneratorLednicky::Set3BodyOff()  {fItest=1;fI3c=0;FsiInit();fWeightDen=1.;FsiNucl();}

AliFemtoModelWeightGenerator*
AliFemtoModelWeightGeneratorLednicky::Clone() const
{
  AliFemtoModelWeightGenerator* tmp = new AliFemtoModelWeightGeneratorLednicky(*this);
  return tmp;
}
