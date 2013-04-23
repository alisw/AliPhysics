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
#define fsiin F77_NAME(fsiin,FSIIN)
extern "C" {void type_of_call F77_NAME(fsiin,FSIIN)(const int &itest,const int &ich, const int &iqs, const int &isi,const int &i3c);}
#define llini F77_NAME(llini,LLINI)
extern "C" {void type_of_call F77_NAME(llini,LLINI)(const int &lll,const int &ns, const int &itest);}

#define fsinucl F77_NAME(fsinucl,FSINUCL)
extern "C" {void type_of_call  F77_NAME(fsinucl,FSINUCL)(const double &mn,const double &cn);}
#define fsimomentum F77_NAME(fsimomentum,FSIMOMENTUM)
extern "C" {void type_of_call F77_NAME(fsimomentum,FSIMOMENTUM)(double &p1,double &p2);}
#define fsiposition F77_NAME(fsiposition,FSIPOSITION)
extern "C" {void type_of_call F77_NAME(fsiposition,FSIPOSITION)(double &x1,double &x2);}
#define fsiw F77_NAME(fsiw,FSIW)
extern "C" {void type_of_call F77_NAME(fsiw,FSIW)(const int &i,double &weif,
						  double &wei,double &wein);}
#define ltran12 F77_NAME(ltran12,LTRAN12)
extern "C" {void type_of_call ltran12_();}

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
ClassImp(AliFemtoModelWeightGeneratorLednicky)
#endif

AliFemtoModelWeightGeneratorLednicky::AliFemtoModelWeightGeneratorLednicky() : 
  AliFemtoModelWeightGenerator(),
  fWei(0), fWein(0), fWeif(0), fWeightDen(0), 
  fItest(0),fIch(1),fIqs(1),fIsi(1),fI3c(0),
  fNuclMass(1.),fNuclCharge(0.),
  fSphereApp(false),fT0App(false) ,
  fLL(0), fNuclChargeSign(1), fSwap(0), fLLMax(30), fLLName(0), 
  fNumProcessPair(0), fNumbNonId(0)
{
  // default constructor
  fLLName=new char*[fLLMax+1];
  fNumProcessPair=new int[fLLMax+1];
  int i;
  for (i=1;i<=fLLMax;i++) {fLLName[i]=new char[40];fNumProcessPair[i]=0;}
  strncpy( fLLName[1],"neutron neutron",40);
  strncpy( fLLName[2],"proton proton",40);
  strncpy( fLLName[3],"neutron proton",40);
  strncpy( fLLName[4],"alpha alpha",40);
  strncpy( fLLName[5],"pi+ pi-",40);
  strncpy( fLLName[6],"pi0 pi0",40);
  strncpy( fLLName[7],"pi+ pi+",40);
  strncpy( fLLName[8],"neutron deuteron",40);
  strncpy( fLLName[9],"proton deuteron",40);
  strncpy( fLLName[10],"pi+ K-",40);
  strncpy( fLLName[11],"pi+ K+",40);
  strncpy( fLLName[12],"pi+ proton",40);
  strncpy( fLLName[13],"pi- proton",40);
  strncpy( fLLName[14],"K+ K-",40);
  strncpy( fLLName[15],"K+ K+",40);
  strncpy( fLLName[16],"K+ proton",40);
  strncpy( fLLName[17],"K- proton",40);
  strncpy( fLLName[18],"deuteron deuteron",40);
  strncpy( fLLName[19],"deuton alpha",40);
  strncpy( fLLName[20],"triton triton",40);
  strncpy( fLLName[21],"triton alpha",40);
  strncpy( fLLName[22],"K0 K0",40);
  strncpy( fLLName[23],"K0 K0b",40);
  strncpy( fLLName[24],"deuteron triton",40);
  strncpy( fLLName[25],"proton triton",40);
  strncpy( fLLName[26],"proton alpha",40);
  strncpy( fLLName[27],"proton lambda",40);
  strncpy( fLLName[28],"neutron lambda",40);
  strncpy( fLLName[29],"Lambda lambda",40);// gael 21May02
  strncpy( fLLName[30],"Proton Anti-proton",40);// gael 21May02
  FsiInit();
  FsiNucl();
}
//______________________
AliFemtoModelWeightGeneratorLednicky::AliFemtoModelWeightGeneratorLednicky(const AliFemtoModelWeightGeneratorLednicky &aWeight):
  AliFemtoModelWeightGenerator(),
  fWei(0), fWein(0), fWeif(0), fWeightDen(0), 
  fItest(0),fIch(1),fIqs(1),fIsi(1),fI3c(0),
  fNuclMass(1.),fNuclCharge(0.),
  fSphereApp(false),fT0App(false) ,
  fLL(0), fNuclChargeSign(1), fSwap(0), fLLMax(30), fLLName(0), 
  fNumProcessPair(0), fNumbNonId(0)
{
  // copy constructor
  fWei = aWeight.fWei; 
  fWein = aWeight.  fWein;
  fWeif = aWeight. fWeif;
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
  fNuclChargeSign = aWeight.fNuclChargeSign;
  fSwap = aWeight.fSwap;
  fLLName = aWeight.fLLName; 
  fNumProcessPair = aWeight.fNumProcessPair;
  fNumbNonId = aWeight.fNumbNonId;
  fLLName=new char*[fLLMax+1];
  fNumProcessPair=new int[fLLMax+1];
  int i;
  for (i=1;i<=fLLMax;i++) {fLLName[i]=new char[40];fNumProcessPair[i]=0;}
  strncpy( fLLName[1],"neutron neutron",40);
  strncpy( fLLName[2],"proton proton",40);
  strncpy( fLLName[3],"neutron proton",40);
  strncpy( fLLName[4],"alpha alpha",40);
  strncpy( fLLName[5],"pi+ pi-",40);
  strncpy( fLLName[6],"pi0 pi0",40);
  strncpy( fLLName[7],"pi+ pi+",40);
  strncpy( fLLName[8],"neutron deuteron",40);
  strncpy( fLLName[9],"proton deuteron",40);
  strncpy( fLLName[10],"pi+ K-",40);
  strncpy( fLLName[11],"pi+ K+",40);
  strncpy( fLLName[12],"pi+ proton",40);
  strncpy( fLLName[13],"pi- proton",40);
  strncpy( fLLName[14],"K+ K-",40);
  strncpy( fLLName[15],"K+ K+",40);
  strncpy( fLLName[16],"K+ proton",40);
  strncpy( fLLName[17],"K- proton",40);
  strncpy( fLLName[18],"deuteron deuteron",40);
  strncpy( fLLName[19],"deuton alpha",40);
  strncpy( fLLName[20],"triton triton",40);
  strncpy( fLLName[21],"triton alpha",40);
  strncpy( fLLName[22],"K0 K0",40);
  strncpy( fLLName[23],"K0 K0b",40);
  strncpy( fLLName[24],"deuteron triton",40);
  strncpy( fLLName[25],"proton triton",40);
  strncpy( fLLName[26],"proton alpha",40);
  strncpy( fLLName[27],"proton lambda",40);
  strncpy( fLLName[28],"neutron lambda",40);
  strncpy( fLLName[29],"Lambda lambda",40);// gael 21May02
  strncpy( fLLName[30],"Proton Anti-proton",40);// gael 21May02
  FsiInit();
  FsiNucl();
}

AliFemtoModelWeightGeneratorLednicky& AliFemtoModelWeightGeneratorLednicky::operator=(const AliFemtoModelWeightGeneratorLednicky& aWeight)
{
  // assignment operator
  if (this == &aWeight)
    return *this;

  fWei = aWeight.fWei; 
  fWein = aWeight.  fWein;
  fWeif = aWeight. fWeif;
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
  fNuclChargeSign = aWeight.fNuclChargeSign;
  fSwap = aWeight.fSwap;
  //  fLLName = aWeight.fLLName; 
  fNumProcessPair = aWeight.fNumProcessPair;
  fNumbNonId = aWeight.fNumbNonId;
  if (fLLName) free(fLLName);
  fLLName=new char*[fLLMax+1];
  if (fNumProcessPair) free(fNumProcessPair);
  fNumProcessPair=new int[fLLMax+1];
  int i;
  for (i=1;i<=fLLMax;i++) {fLLName[i]=new char[40];fNumProcessPair[i]=0;}
  strncpy( fLLName[1],"neutron neutron",40);
  strncpy( fLLName[2],"proton proton",40);
  strncpy( fLLName[3],"neutron proton",40);
  strncpy( fLLName[4],"alpha alpha",40);
  strncpy( fLLName[5],"pi+ pi-",40);
  strncpy( fLLName[6],"pi0 pi0",40);
  strncpy( fLLName[7],"pi+ pi+",40);
  strncpy( fLLName[8],"neutron deuteron",40);
  strncpy( fLLName[9],"proton deuteron",40);
  strncpy( fLLName[10],"pi+ K-",40);
  strncpy( fLLName[11],"pi+ K+",40);
  strncpy( fLLName[12],"pi+ proton",40);
  strncpy( fLLName[13],"pi- proton",40);
  strncpy( fLLName[14],"K+ K-",40);
  strncpy( fLLName[15],"K+ K+",40);
  strncpy( fLLName[16],"K+ proton",40);
  strncpy( fLLName[17],"K- proton",40);
  strncpy( fLLName[18],"deuteron deuteron",40);
  strncpy( fLLName[19],"deuton alpha",40);
  strncpy( fLLName[20],"triton triton",40);
  strncpy( fLLName[21],"triton alpha",40);
  strncpy( fLLName[22],"K0 K0",40);
  strncpy( fLLName[23],"K0 K0b",40);
  strncpy( fLLName[24],"deuteron triton",40);
  strncpy( fLLName[25],"proton triton",40);
  strncpy( fLLName[26],"proton alpha",40);
  strncpy( fLLName[27],"proton lambda",40);
  strncpy( fLLName[28],"neutron lambda",40);
  strncpy( fLLName[29],"Lambda lambda",40);// gael 21May02
  strncpy( fLLName[30],"Proton Anti-proton",40);// gael 21May02
  FsiInit();
  FsiNucl();
  
  return *this;
}


double AliFemtoModelWeightGeneratorLednicky::GenerateWeight(AliFemtoPair* aPair)
{
  // Get hidden information pointers
  //AliFemtoModelHiddenInfo *inf1 = (AliFemtoModelHiddenInfo *) aPair->Track1()->HiddenInfo();
  //AliFemtoModelHiddenInfo *inf2 = (AliFemtoModelHiddenInfo *) aPair->Track2()->HiddenInfo();
  AliFemtoTrack *inf1 = (AliFemtoTrack *) aPair->Track1()->Track();
  AliFemtoTrack *inf2 = (AliFemtoTrack *) aPair->Track2()->Track();

  // Calculate pair variables
  Double_t tPx = inf1->GetTrueMomentum()->x()+inf2->GetTrueMomentum()->x();
  Double_t tPy = inf1->GetTrueMomentum()->y()+inf2->GetTrueMomentum()->y();
  Double_t tPz = inf1->GetTrueMomentum()->z()+inf2->GetTrueMomentum()->z();
  Double_t tM1 = inf1->GetMass();
  Double_t tM2 = inf2->GetMass();
  Double_t tE1 = sqrt(tM1*tM1 + inf1->GetTrueMomentum()->Mag2());
  Double_t tE2 = sqrt(tM2*tM2 + inf2->GetTrueMomentum()->Mag2());
  Double_t tE  = tE1 + tE2;
  Double_t tPt = tPx*tPx + tPy*tPy;
  Double_t tMt = tE*tE - tPz*tPz;//mCVK;
  Double_t tM  = sqrt(tMt - tPt);
  tMt = sqrt(tMt);
  tPt = sqrt(tPt);
  Double_t tBetat = tPt/tMt;

  // Boost to LCMS
  Double_t tBeta = tPz/tE;
  Double_t tGamma = tE/tMt;	    
  fKStarLong = tGamma * (inf1->GetTrueMomentum()->z() - tBeta * tE1);
  Double_t tE1L = tGamma * (tE1  - tBeta * inf1->GetTrueMomentum()->z());
    
  // Rotate in transverse plane
  fKStarOut  = ( inf1->GetTrueMomentum()->x()*tPx + inf1->GetTrueMomentum()->y()*tPy)/tPt;
  fKStarSide = (-inf1->GetTrueMomentum()->x()*tPy + inf1->GetTrueMomentum()->y()*tPx)/tPt;
      
  // Boost to pair cms
  fKStarOut = tMt/tM * (fKStarOut - tPt/tMt * tE1L);
  
  tBetat = tPt/tMt;
  
  Double_t tDX = inf1->GetEmissionPoint()->x()-inf2->GetEmissionPoint()->x();
  Double_t tDY = inf1->GetEmissionPoint()->y()-inf2->GetEmissionPoint()->y();
  Double_t tRLong = inf1->GetEmissionPoint()->z()-inf2->GetEmissionPoint()->z();
  Double_t tDTime = inf1->GetEmissionPoint()->t()-inf2->GetEmissionPoint()->t();

  Double_t tROut = (tDX*tPx + tDY*tPy)/tPt;
  Double_t tRSide = (-tDX*tPy + tDY*tPx)/tPt;

//   cout << "Got points 1 " << inf1->GetEmissionPoint()->x() << "  " <<  inf1->GetEmissionPoint()->y() << " "  << inf1->GetEmissionPoint()->z() << "  " << inf1->GetEmissionPoint()->t() << endl;

//   cout << "Got points 2 " << inf2->GetEmissionPoint()->x() << "  " << inf2->GetEmissionPoint()->y() << " " << inf2->GetEmissionPoint()->z() << "  " << inf2->GetEmissionPoint()->t() << endl;

  fRStarSide = tRSide;

  fRStarLong = tGamma*(tRLong - tBeta* tDTime);
  Double_t tDTimePairLCMS = tGamma*(tDTime - tBeta* tRLong);

  tBeta = tPt/tMt;
  tGamma = tMt/tM;

  fRStarOut = tGamma*(tROut - tBeta* tDTimePairLCMS);
  fRStar = ::sqrt(fRStarOut*fRStarOut + fRStarSide*fRStarSide +
			   fRStarLong*fRStarLong);
  fKStar = ::sqrt(fKStarOut*fKStarOut + fKStarSide*fKStarSide + fKStarLong*fKStarLong);

//   cout << "Got out side " << fRStarOut << " " << fRStarSide << endl;

  if (!SetPid(inf1->GetPDGPid(),inf2->GetPDGPid())) {
    fWeightDen=1.;
    return 1;    
  } 
  else { // Good Pid
    AliFemtoThreeVector*  p;
    p=(inf1->GetTrueMomentum());
    double p1[]={p->x(),p->y(),p->z()};
    p=(inf2->GetTrueMomentum());
    double p2[]={p->x(),p->y(),p->z()};
    if ((p1[0]==p2[0])&&(p1[1]==p2[1])&&(p1[2]==p2[2])) {
      fWeightDen=0.;
      return 0;  
    } 
    if (fSwap) {
      fsimomentum(*p2,*p1);
    } else {
      fsimomentum(*p1,*p2);
    }
    AliFemtoLorentzVector* tPoint;
    tPoint=(inf1->GetEmissionPoint());
//     cout << "Pid1:dans GetWeight = " << aThPair->GetPid1() << endl;
//     cout << "Pid2:dans GetWeight = " << aThPair->GetPid2() << endl;
//     cout << "LL:in GetWeight = " << mLL << endl;

    double x1[]={tPoint->x(),tPoint->y(),tPoint->z(),tPoint->t()};
    tPoint=(inf2->GetEmissionPoint());
    double x2[]={tPoint->x(),tPoint->y(),tPoint->z(),tPoint->t()};
    if ((x1[0]==x2[0])&&(x1[1]==x2[1])&&(x1[2]==x2[2])&&(x1[3]==x2[3])) {
      fWeightDen=0.;
      return 0;  
    } 
    if (fSwap) {
      fsiposition(*x2,*x1);
    } else {
      fsiposition(*x1,*x2);
    }
    FsiSetLL();
    ltran12();
    fsiw(1,fWeif,fWei,fWein);

    if (fI3c==0) return fWein;
    fWeightDen=fWeif;
    return fWei;
  }
}

AliFemtoString AliFemtoModelWeightGeneratorLednicky::Report() {
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
  int i;
  for(i=1;i<=fLLMax;i++) { 
    if (fNumProcessPair[i])
      tStr << "         " << fNumProcessPair[i] << " " << fLLName[i] << endl;
  }
  if (fNumbNonId)
    tStr << "         "<< fNumbNonId << " Non Identified" << endl;
  AliFemtoString returnThis = tStr.str();
  return returnThis;
}

void AliFemtoModelWeightGeneratorLednicky::FsiInit(){
  // Initialize weight generation module
//   cout << "*******************AliFemtoModelWeightGeneratorLednicky check FsiInit ************" << endl;
//   cout <<"mItest dans FsiInit() = " << fItest << endl;
//   cout <<"mIch dans FsiInit() = " << fIch << endl;
//   cout <<"mIqs dans FsiInit() = " << fIqs << endl;
//   cout <<"mIsi dans FsiInit() = " << fIsi << endl;
//   cout <<"mI3c dans FsiInit() = " << fI3c << endl;
  fsiin(fItest,fIch,fIqs,fIsi,fI3c);
}

void AliFemtoModelWeightGeneratorLednicky::FsiNucl(){
  // initialize weight generation taking into account the residual charge
//   cout << "*******************AliFemtoModelWeightGeneratorLednicky check FsiNucl ************" << endl;
//   cout <<"fNuclMass dans FsiNucl() = " << fNuclMass << endl;
//   cout <<"fNuclCharge dans FsiNucl() = " << fNuclCharge << endl;
//   cout <<"fNuclChargeSign dans FsiNucl() = " << fNuclChargeSign << endl;
  fsinucl(fNuclMass,fNuclCharge*fNuclChargeSign);
}

void AliFemtoModelWeightGeneratorLednicky::FsiSetLL(){
  // set internal pair type for the module
  int tNS;
  if (fSphereApp||(fLL>5)) {
    if (fT0App) { tNS=4;} 
    else {tNS=2;}
  } else { tNS=1;}
   //cout <<"fLL dans FsiSetLL() = "<< fLL << endl;
   //cout <<"tNS dans FsiSetLL() = "<< tNS << endl;
   //cout <<"fItest dans FsiSetLL() = "<< fItest << endl;
  llini(fLL,tNS,fItest);
  //cout<<" end of FsiSetLL"<<endl;
}
         
bool AliFemtoModelWeightGeneratorLednicky::SetPid(const int aPid1,const int aPid2) {
  // set calculated system for basing on particles' pids
  static const int ksPi0Pid=111;
  static const int ksPionPid=211; 
  static const int ksK0Pid=311;
  static const int ksKPid=321;
  static const int ksNeutPid=2112;
  static const int ksProtPid=2212;
  static const int ksLamPid=3122;
  //  static const int sLamLamPid=3122;

   // cout << "Setting PID to " << aPid1 << " " << aPid2 << endl;

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
  default: fLL=0;
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
//   cout << "fLL=="<< fLL << endl;
//   cout << "fNuclCharge=="<< fNuclCharge << endl;

}    
AliFemtoModelWeightGeneratorLednicky::~AliFemtoModelWeightGeneratorLednicky() 
{ 
  if (fLLName) delete [] fLLName;
  if (fNumProcessPair) delete [] fNumProcessPair;
/* no-op */ 
}

//_____________________________________________
void     AliFemtoModelWeightGeneratorLednicky::SetPairType(Int_t aPairType)
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

  if      (((ktPid1 ==   211) && (ktPid2 ==   211)) ||
           ((ktPid1 ==  -211) && (ktPid2 ==  -211)))
    fPairType = fgkPionPlusPionPlus;
  else if (((ktPid1 ==  -211) && (ktPid2 ==   211)) ||
           ((ktPid1 ==   211) && (ktPid2 ==  -211)))
    fPairType = fgkPionPlusPionMinus;
  else if (((ktPid1 ==   321) && (ktPid2 ==   321)) ||
           ((ktPid1 ==  -321) && (ktPid2 ==  -321)))
    fPairType = fgkKaonPlusKaonPlus;
  else if (((ktPid1 ==  -321) && (ktPid2 ==   321)) ||
           ((ktPid1 ==   321) && (ktPid2 ==  -321)))
    fPairType = fgkKaonPlusKaonMinus;
  else if (((ktPid1 ==  2212) && (ktPid2 ==  2212)) ||
           ((ktPid1 == -2212) && (ktPid2 == -2212)))
    fPairType = fgkProtonProton;
  else if (((ktPid1 == -2212) && (ktPid2 ==  2212)) ||
           ((ktPid1 ==  2212) && (ktPid2 == -2212)))
    fPairType = fgkProtonAntiproton;
  else if (((ktPid1 ==   211) && (ktPid2 ==   321)) ||
           ((ktPid1 ==  -211) && (ktPid2 ==  -321)))
    fPairType = fgkPionPlusKaonPlus;
  else if (((ktPid1 ==  -211) && (ktPid2 ==   321)) ||
           ((ktPid1 ==   211) && (ktPid2 ==  -321)))
    fPairType = fgkPionPlusKaonMinus;
  else if (((ktPid1 ==   211) && (ktPid2 ==  2212)) ||
           ((ktPid1 ==  -211) && (ktPid2 == -2212)))
    fPairType = fgkPionPlusProton;
  else if (((ktPid1 ==  -211) && (ktPid2 ==  2212)) ||
           ((ktPid1 ==   211) && (ktPid2 == -2212)))
    fPairType = fgkPionPlusAntiproton;
  else if (((ktPid1 ==   321) && (ktPid2 ==  2212)) ||
           ((ktPid1 ==  -321) && (ktPid2 == -2212)))
    fPairType = fgkKaonPlusProton;
  else if (((ktPid1 ==  -321) && (ktPid2 ==  2212)) ||
           ((ktPid1 ==   321) && (ktPid2 == -2212)))
    fPairType = fgkKaonPlusAntiproton;
  SetPid(ktPid1, ktPid2);
}

void AliFemtoModelWeightGeneratorLednicky::SetNuclCharge(const double aNuclCharge) {fNuclCharge=aNuclCharge;FsiNucl();}
void AliFemtoModelWeightGeneratorLednicky::SetNuclMass(const double aNuclMass){fNuclMass=aNuclMass;FsiNucl();}

void AliFemtoModelWeightGeneratorLednicky::SetSphere(){fSphereApp=true;}
void AliFemtoModelWeightGeneratorLednicky::SetSquare(){fSphereApp=false;}
void AliFemtoModelWeightGeneratorLednicky::SetT0ApproxOn(){ fT0App=true;}
void AliFemtoModelWeightGeneratorLednicky::SetT0ApproxOff(){ fT0App=false;}
void AliFemtoModelWeightGeneratorLednicky::SetDefaultCalcPar(){
  fItest=1;fIqs=1;fIsi=1;fI3c=0;fIch=1;FsiInit();
  fSphereApp=false;fT0App=false;}

void AliFemtoModelWeightGeneratorLednicky::SetCoulOn()    {fItest=1;fIch=1;FsiInit();}
void AliFemtoModelWeightGeneratorLednicky::SetCoulOff()   {fItest=1;fIch=0;FsiInit();}
void AliFemtoModelWeightGeneratorLednicky::SetQuantumOn() {fItest=1;fIqs=1;FsiInit();}
void AliFemtoModelWeightGeneratorLednicky::SetQuantumOff(){fItest=1;fIqs=0;FsiInit();}
void AliFemtoModelWeightGeneratorLednicky::SetStrongOn()  {fItest=1;fIsi=1;FsiInit();}
void AliFemtoModelWeightGeneratorLednicky::SetStrongOff() {fItest=1;fIsi=0;FsiInit();}
void AliFemtoModelWeightGeneratorLednicky::Set3BodyOn()   {fItest=1;fI3c=1;FsiInit();FsiNucl();}
void AliFemtoModelWeightGeneratorLednicky::Set3BodyOff()  {fItest=1;fI3c=0;FsiInit();fWeightDen=1.;FsiNucl();}

Double_t AliFemtoModelWeightGeneratorLednicky::GetKStar() const {return AliFemtoModelWeightGenerator::GetKStar();}
Double_t AliFemtoModelWeightGeneratorLednicky::GetKStarOut() const { return AliFemtoModelWeightGenerator::GetKStarOut(); }
Double_t AliFemtoModelWeightGeneratorLednicky::GetKStarSide() const { return AliFemtoModelWeightGenerator::GetKStarSide(); }
Double_t AliFemtoModelWeightGeneratorLednicky::GetKStarLong() const { return AliFemtoModelWeightGenerator::GetKStarLong(); }
Double_t AliFemtoModelWeightGeneratorLednicky::GetRStar() const { return AliFemtoModelWeightGenerator::GetRStar(); }
Double_t AliFemtoModelWeightGeneratorLednicky::GetRStarOut() const { return AliFemtoModelWeightGenerator::GetRStarOut(); }
Double_t AliFemtoModelWeightGeneratorLednicky::GetRStarSide() const { return AliFemtoModelWeightGenerator::GetRStarSide(); }
Double_t AliFemtoModelWeightGeneratorLednicky::GetRStarLong() const { return AliFemtoModelWeightGenerator::GetRStarLong(); }

AliFemtoModelWeightGenerator* AliFemtoModelWeightGeneratorLednicky::Clone() const {
  AliFemtoModelWeightGenerator* tmp = new AliFemtoModelWeightGeneratorLednicky(*this);
  return tmp;
}
