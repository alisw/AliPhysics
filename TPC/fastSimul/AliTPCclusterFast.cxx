/// \class AliTPCclusterFast
///
///  Code for the fast MC simulation of the TPC response
///
///  Microscopic simulation includes all features from the FULL MC simulation:
///    - primary ionization
///    - secondary inization
///    - diffusion
///    - gas gain with fluctuation
///    - pad and time response fucntion
///
///  Usage:
///    -  fast sensitivity studies - dEdx
///    -    development of new algorithm
///    -  resolution studies
///    -  NEW: double track resolutuion studies - are currently implemented
///
/// Documentation for the reconstruction part   - we will refer to the  CHEP paper - (http://arxiv.org/pdf/physics/0306108.pdf
   
/*
   Example usage in the test mode: 
  .x ~/rootlogon.C
  
  .L $ALICE_ROOT/../src/TPC/fastSimul/AliTPCclusterFast.cxx+
  AliTPCclusterFast::InitFormulas();
  AliTPCclusterFast::fPRF = new TF1("fprf","gausn",-5,5);  //MWPC
  AliTPCclusterFast::fPRF->SetParameters(1,0,0.5);
  // 
  AliTPCclusterFast::fPRF =new TF1("fprf"," GEMPRF(x,0.02)",-3,3); // GEM
  //
  //AliTPCclusterFast::fTRF = new TF1("ftrf","gausn",-5,5);
  AliTPCclusterFast::fTRF = new TF1("gamma4Norm","Gamma4Norm(x)",-3,3); //Gamma4 TRF in bin units (100ns)

  TStopwatch timer;  AliTPCtrackFast::Simul("trackerSimul.root",20,0.6);   timer.Print();

  TFile * ftrack = TFile::Open("trackerSimul.root");
  TTree *tree  = (TTree*)ftrack->Get("simulTrack");
  tree->SetMarkerStyle(25);
  tree->SetMarkerSize(0.6);
  testTree=tree;
  AliTPCclusterFast::SetMetadata(tree);
  AliTPCclusterFast::UnitTest();

*/
/// ~~~
///
///  Modifications to add:
///  0.)   Cluster Unfolding methods (from  (http://arxiv.org/pdf/physics/0306108.pdf) )
///  1.)   Track unfolding methods
///  4.)   Create hardware setup class
///        (fNoise, fGain, fBRounding, fBAddpedestal, ....)
///  5.)   Create arrays of registered hardware setups
///  6.)   Extend on the fly functions to use registered hardware setups, identified by ID.
///        hwMode

#include "TObject.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVectorF.h"
#include "TMatrixF.h"
#include "TH1.h"
#include "AliTPCreco.h"
#include "TClonesArray.h"
#include "TTreeStream.h"
#include "TGrid.h"
#include "TStatToolkit.h"
#include "AliTPCParamSR.h"
#include "TFormulaPrimitive.h"

TTree* testTree=0;

class AliTPCclusterFast: public TObject {
public:
  AliTPCclusterFast();
  void Init();
  
  virtual ~AliTPCclusterFast();
  void SetParam(Float_t mnprim, Float_t diff, Float_t diffL, Int_t padrow, Float_t y, Float_t z, Float_t ky, Float_t kz, Float_t yCenter, Float_t zCenter);
  static void GenerElectrons(AliTPCclusterFast *cl0, AliTPCclusterFast *clm, AliTPCclusterFast *clp);
  void Digitize();
  Double_t GetQtot(Float_t gain,Float_t thr, Float_t noise, Bool_t rounding=kTRUE, Bool_t addPedestal=kTRUE);
  Double_t GetQmax(Float_t gain,Float_t thr, Float_t noise, Bool_t rounding=kTRUE, Bool_t addPedestal=kTRUE);
  Double_t GetQmaxCorr(Float_t rmsy0, Float_t rmsz0);
  Double_t GetQtotCorr(Float_t rmsy0, Float_t rmsz0, Float_t gain, Float_t thr);
  Double_t GetClusterProperties(TVectorF &param, Int_t addOverlap, Float_t gain=0.8, Float_t thr=2, Float_t noise=0.7, Bool_t rounding=kTRUE, Bool_t addPedestal=kTRUE, Int_t skipSample=0);
  Double_t GetCOG(Int_t returnType, Int_t addOverlap, Float_t gain=0.8, Float_t thr=2, Float_t noise=0.7, Bool_t rounding=kTRUE, Bool_t addPedestal=kTRUE,Int_t skipSample=0);
  Double_t GetCOGHit(Int_t returnType);
  Float_t  UnfoldCluster(Float_t & meani, Float_t & meanj,  Float_t & sumu, Float_t & overlap, Int_t addOverlap, Float_t gain=0.8, Float_t thr=2, Float_t noise=0.7, Bool_t rounding=kTRUE, Bool_t addPedestal=kTRUE,Int_t skipSample=0);
  Float_t  GetCOGUnfolded(Int_t returnType, Int_t addOverlap=kTRUE, Float_t gain=0.8, Float_t thr=2, Float_t noise=0.7, Bool_t rounding=kTRUE, Bool_t addPedestal=kTRUE,Int_t skipSample=0);

  Bool_t MakeDigitization(Int_t addOverlap, Float_t gain, Float_t thr, Float_t noise, Bool_t rounding, Bool_t addPedestal,Int_t skipSample);
  //
  Double_t  GetExpectedRMS(Int_t dim);
  //Double_t  GetExpectedResolution(Int_t dim);

  Double_t GetNsec();
  Int_t GetYMaxBin(){return fYMaxBin;}
  Int_t GetZMaxBin(){return fZMaxBin;}
  TMatrixF *GetRawDigits(){return &fRawDigits;}
  Float_t GetRawDigit(Int_t i, Int_t j){return fRawDigits(i,j);};
  Float_t GetDigit(Int_t i, Int_t j){return fDigits(i,j);};
  Float_t GetDigitsRawMax(){return TMath::MaxElement(35,fRawDigits.GetMatrixArray());}
  Float_t GetDigitsMax(){return TMath::MaxElement(35,fDigits.GetMatrixArray());}
  //static void Simul(const char* simul, Int_t npoints);
  static void InitFormulas();
  static Double_t GaussConvolution(Double_t x0, Double_t x1, Double_t k0, Double_t k1, Double_t s0, Double_t s1);
  static Double_t GaussExpConvolution(Double_t x0, Double_t s0,Double_t t1);
  static Double_t GaussGamma4(Double_t x, Double_t s0, Double_t p1);
  static Double_t Gamma4(Double_t x, Double_t p0, Double_t p1); 
  static Double_t Gamma4Norm(Double_t x); 
  static Double_t GEMPRF(Double_t x, Double_t sigma);
  static void SetMetadata(TTree * tree);
  //
  //
  static Bool_t UnitTest();
public: // Cluster parameters expressid in bin units
  //
  Float_t fMNprim;     ///< mean number of primary electrons
  //                   //electrons part input
  Int_t   fNprim;      ///< mean number of primary electrons
  Int_t   fNtot;       ///< total number of  electrons
  Float_t fQtot;       ///< total charge - Gas gain flucuation taken into account
  //
  Float_t fDiff;       ///< diffusion sigma
  Float_t fDiffLong;   ///< diffusion sigma longitudinal direction
  Int_t   fPadRow;        ///< cluster pad row number
  Float_t fY;          ///< y ideal position - center bin
  Float_t fZ;          ///< z ideal position - center bin
  Float_t fYCenterBin; ///< y center bin   
  Float_t fZCenterBin; ///< z center bin   
  Int_t fYMaxBin;    //! y maximum position as position as calculated in the  MakeDigitization   (should be 0 in mean )  
  Int_t fZMaxBin;    //! z maximum position as position as calculated in the  MakeDigitization   (should be 0 in mean )
  
  Float_t fAngleY;     ///< y angle - tan(y)
  Float_t fAngleZ;     ///< z angle - tan z
  AliTPCclusterFast *fOverlapCluster; //
  //
  //
  //                   // electron part simul
  TVectorF fSec;       //!< number of secondary electrons
  TVectorF fPosY;      //!< position y for each electron
  TVectorF fPosZ;      //!< position z for each electron
  TVectorF fGain;      //!< gg for each electron
  //
  TVectorF fStatY;     // stat Y
  TVectorF fStatZ;     // stat Y
  //
  // digitization part
  //
  TMatrixF fDigits;    ///< response matrix (ideal signal without noise, trheshold rounding, pile-up)
  TMatrixF fRawDigits; //!  "on the fly caclualted digits matrix" for current version of setup - gain, noise, thr.
  Float_t fPRFRMS;     /// pad respons function width
  Float_t fTRFRMS;     /// time rsponsefunction width
  static TF1* fPRF;    ///< Pad response
  static TF1* fTRF;    ///< Time response function
  static Float_t fgZSamplingFactor;               // z sample at 10 MHz is 2 time narrower as rphi sample 
  static Int_t fgDebugLevel;         // debug output downscaling  factor
  ClassDef(AliTPCclusterFast,2)  // container for
};


class AliTPCtrackFast: public TObject {
public:
  AliTPCtrackFast();
  void MakeTrack();
  static void Simul(const char* simul, Int_t ntracks, Double_t diff, Bool_t simulOverlap=kTRUE);
  Double_t  CookdEdxNtot(Double_t f0,Float_t f1);
  Double_t  CookdEdxQtot(Double_t f0,Float_t f1);
  Double_t  CookdEdxNtotThr(Double_t f0,Float_t f1, Double_t thr, Int_t dEdxMode);
  Double_t  CookdEdxQtotThr(Double_t f0,Float_t f1, Double_t thr, Int_t dEdxMode);
  //
  Double_t  CookdEdxDtot(Double_t f0,Float_t f1, Float_t gain,Float_t thr, Float_t noise, Bool_t corr, Int_t dEdxMode);
  Double_t  CookdEdxDmax(Double_t f0,Float_t f1,Float_t gain,Float_t thr, Float_t noise, Bool_t corr, Int_t dEdxMode);
  //
  Double_t  CookdEdx(Int_t npoints, Double_t *amp, Double_t f0,Float_t f1, Int_t dEdxMode);
  //
  Float_t fMNprim;     ///< mean number of primary electrons per cm
  Float_t fTY;         ///< track Y at the vertex in (cm)
  Float_t fTZ;         ///< track Z at the vertex in (cm)
  Float_t fTAngleY;     ///< y angle - tan(y) - dy/dx (cm/cm)
  Float_t fTAngleZ;     ///< z angle - tan(z) - dy/dx (cm/cm) 
  Float_t fDiff;       ///< diffusion  in mm/sqrt(m) - nominal is 2.2 mm/sqrt(m)
  Float_t fDiffLong;       ///< diffusion sigma longitudinal direction
  Int_t   fN;          ///< number of clusters simulated
  //  overlap track properties
  Bool_t  fBOverlap;          ///< flag generate overlap track
  Float_t fMNprimOverlap;     ///< mean number of primary electrons for overlap track
  Float_t fDAngleYOverlap;    ///< y angle - tan(y) for overlap track  in cm
  Float_t fDAngleZOverlap;    ///< z angle - tan z for overlap track   in cm
  Float_t fDYOverlap;         ///< delta y position of overlap track at row 0 in dy/dx 
  Float_t fDZOverlap;         ///< delta z position of overlap track at row 0 in dz/dx
  TClonesArray *fCl;          ///< array of clusters
  //
  Bool_t   fInit;      ///< initialization flag

  /// \cond CLASSIMP
  ClassDef(AliTPCtrackFast,2)
  /// \endcond
};

/// \cond CLASSIMP
ClassImp(AliTPCclusterFast)
ClassImp(AliTPCtrackFast)
/// \endcond

TF1 *AliTPCclusterFast::fPRF=0;
TF1 *AliTPCclusterFast::fTRF=0;
Float_t AliTPCclusterFast::fgZSamplingFactor=2;
Int_t AliTPCclusterFast::fgDebugLevel=0;
AliTPCParamSR paramSR;    // parameter class to get TPC properties



AliTPCtrackFast::AliTPCtrackFast():
  TObject(),
  fMNprim(0),
  fTAngleY(0),
  fTAngleZ(0),
  fN(0),
  fBOverlap(kFALSE),          ///< flag generate overlap track
  fMNprimOverlap(0),     ///< mean number of primary electrons for overlap track
  fDAngleYOverlap(0),    ///< y angle - tan(y) for overlap track
  fDAngleZOverlap(0),    ///< z angle - tan z for overlap track
  fDYOverlap(0),         ///< delta y position of overlap track at row 0 
  fDZOverlap(0),         ///< delta z position of overlap track at row 0 
  fCl(0),
  fInit(kFALSE)
{
  ///

}

void AliTPCclusterFast::InitFormulas(){
  //
  // Register formula as a build-in formulas to speed up evaluation
  //
  TFormulaPrimitive::AddFormula(new TFormulaPrimitive("GEMPRF","GEMPRF",AliTPCclusterFast::GEMPRF));
  TFormulaPrimitive::AddFormula(new TFormulaPrimitive("Gamma4","Gamma4",AliTPCclusterFast::Gamma4));
  TFormulaPrimitive::AddFormula(new TFormulaPrimitive("Gamma4Norm","Gamma4Norm",AliTPCclusterFast::Gamma4Norm));
}


void AliTPCtrackFast::MakeTrack(){
  ///
  /// track paramters in absolute units
  /// cluster should be expresse in relative units 
  //    - cluster performnce dependes on relative units normalized to bin size
  //    - to keep back compatibility 
  //
  if (!fCl) fCl = new TClonesArray("AliTPCclusterFast",159);
  //
  // 0.) Init data structure
  //
  for (Int_t irow=0;irow<fN;irow++){
    AliTPCclusterFast * cluster = (AliTPCclusterFast*) fCl->UncheckedAt(irow);
    if (!cluster) cluster =   new ((*fCl)[irow]) AliTPCclusterFast;
    cluster->Init();
  }
  //
  // 1.) Create hits - with crosstalk diffusion
  //
  for (Int_t irow=0;irow<fN;irow++){
    Double_t tX = (irow < paramSR.GetNRow(0)) ?  paramSR.GetPadRowRadiiLow(irow): paramSR.GetPadRowRadiiUp(irow-paramSR.GetNRow(0));
    Double_t tY = fTY+tX*fTAngleY;
    Double_t tZ = fTZ+tX*fTAngleZ;
    Double_t padLength= (irow < paramSR.GetNRow(0)) ? paramSR.GetPadPitchLength(0,irow):paramSR.GetPadPitchLength(36,irow-paramSR.GetNRow(0));
    Double_t padWidth= (irow < paramSR.GetNRow(0)) ? paramSR.GetPadPitchWidth(0):paramSR.GetPadPitchWidth(36);
    Double_t zWidth= paramSR.GetZWidth()*0.5;  // 10 MHz width default

    //
    AliTPCclusterFast * cluster = (AliTPCclusterFast*) fCl->UncheckedAt(irow);
    AliTPCclusterFast * clusterm = (AliTPCclusterFast*) fCl->UncheckedAt(TMath::Max(irow-1,0));
    AliTPCclusterFast * clusterp = (AliTPCclusterFast*) fCl->UncheckedAt(TMath::Min(irow+1,kMaxRow-1));
    if (!cluster) cluster =   new ((*fCl)[irow]) AliTPCclusterFast;
    //
    // Transform cm -> bins and angle to "bin angle"
    // 
    Double_t tYBin=tY/padWidth;
    Double_t tZBin=tZ/zWidth;
    Double_t drift=paramSR.GetZLength()-TMath::Abs(tZ);
    if (drift<0) drift=0;
    Double_t driftFactor= TMath::Sqrt(drift/100.); // diffusion convenrsion - drift length in meters
    Double_t fDiffBin=fDiff*driftFactor/padWidth;
    Double_t fDiffLongBin=fDiffLong*driftFactor/padLength;
    Double_t fAngleYBin=fTAngleY*padLength/padWidth;
    Double_t fAngleZBin=fTAngleZ*padLength/zWidth;
    //
    Float_t  yCenterBin= TMath::Nint(tYBin);
    Float_t  zCenterBin= TMath::Nint(tZBin*0.5)*2;
    Double_t posYBin = tYBin-yCenterBin;
    Double_t posZBin = tZBin-zCenterBin;
   

    cluster->SetParam(fMNprim*padLength,fDiffBin, fDiffLongBin, irow, posYBin,posZBin,fAngleYBin,fAngleZBin,yCenterBin,zCenterBin);   // when is the cluster plus parameters done 
    //
    cluster->GenerElectrons(cluster, clusterm, clusterp);
    if (fBOverlap){   // if we simulate the overlap of tracks
      AliTPCclusterFast * clusterOverlap = cluster->fOverlapCluster;
      if (clusterOverlap==NULL){
	cluster->fOverlapCluster=new AliTPCclusterFast;
	clusterOverlap=	cluster->fOverlapCluster;
      }
      clusterOverlap->Init();
      Double_t posYOverlapBin=posYBin+(fDYOverlap+tX*fDAngleYOverlap)/padWidth;
      Double_t posZOverlapBin=posZBin+(fDZOverlap+tX*fDAngleZOverlap)/padLength;
      clusterOverlap->SetParam(fMNprimOverlap*padLength,fDiffBin, fDiffLongBin, irow,posYOverlapBin,posZOverlapBin,fAngleYBin+fDAngleYOverlap/padWidth,fAngleZBin+fDAngleZOverlap/zWidth,yCenterBin,zCenterBin);      clusterOverlap->GenerElectrons(clusterOverlap, 0, 0);
    }
  }
  //
  // 2.) make digitization
  //
  for (Int_t i=0;i<fN;i++){
    AliTPCclusterFast * cluster = (AliTPCclusterFast*) fCl->UncheckedAt(i);
    cluster->Digitize();
    if (cluster->fOverlapCluster) cluster->fOverlapCluster->Digitize();
  }

}

Double_t  AliTPCtrackFast::CookdEdxNtot(Double_t f0,Float_t f1){
  ///    Double_t  CookdEdxNtot(Double_t f0,Float_t f1);   //  dEdx_{hit}  reconstructed meen number of  electrons

  Double_t amp[159];
  for (Int_t i=0;i<fN;i++){ 
    AliTPCclusterFast * cluster = ( AliTPCclusterFast *)((*fCl)[i]);
    amp[i]=cluster->fNtot;
  }
  return CookdEdx(fN,amp,f0,f1,0);
}

Double_t  AliTPCtrackFast::CookdEdxQtot(Double_t f0,Float_t f1){
  ///     dEdx_{Q} reconstructed mean number of electronsxGain

  Double_t amp[159];
  for (Int_t i=0;i<fN;i++){ 
    AliTPCclusterFast * cluster = ( AliTPCclusterFast *)((*fCl)[i]);
    amp[i]=cluster->fQtot;
  }
  return CookdEdx(fN,amp,f0,f1,0);
}


Double_t  AliTPCtrackFast::CookdEdxNtotThr(Double_t f0,Float_t f1, Double_t thr, Int_t dEdxMode){
  ///   dEdx_{hit}  reconstructed mean number of  electrons
  ///     thr  = threshold in terms of the number of electrons
  ///     dEdxMode = algorithm to deal with trhesold values replacing

  Double_t amp[159];
  Int_t nBellow=0;
  //
  Double_t minAbove=-1;
  for (Int_t i=0;i<fN;i++){ 
    AliTPCclusterFast * cluster = ( AliTPCclusterFast *)((*fCl)[i]);
    Double_t clQ= cluster->fNtot;
    if (clQ<thr) {
      nBellow++;
      continue;
    }
    if (minAbove<0) minAbove=clQ;
    if (minAbove>clQ) minAbove=clQ;
  }
  //
  if (dEdxMode==-1) return Double_t(nBellow)/Double_t(fN);

  for (Int_t i=0;i<fN;i++){ 
    AliTPCclusterFast * cluster = ( AliTPCclusterFast *)((*fCl)[i]);
    Double_t clQ= cluster->fNtot;
    //
    if (dEdxMode==0)  amp[i]=clQ;              // dEdxMode0 - not threshold  - keep default
    //
    //
    if (dEdxMode==1 && clQ>thr) amp[i]=clQ;    // dEdxMode1 - skip if bellow 
    if (dEdxMode==1 && clQ<thr) amp[i]=0;      // dEdxMode1 - skip if bellow 
    //
    //
    if (dEdxMode==2 && clQ>thr) amp[i]=clQ;    // dEdxMode2 - use 0 if below
    if (dEdxMode==2 && clQ<thr) amp[i]=0;      // dEdxMode2 - use 0 if below
    //
    //
    if (dEdxMode==3)  amp[i]=(clQ>thr)?clQ:thr; // dEdxMode3 - use thr if below
    if (dEdxMode==4)  amp[i]=(clQ>thr)?clQ:minAbove; // dEdxMode4 -  use minimal above threshold if bellow thr
  }
  return CookdEdx(fN,amp,f0,f1, dEdxMode);
}



Double_t  AliTPCtrackFast::CookdEdxQtotThr(Double_t f0,Float_t f1, Double_t thr, Int_t dEdxMode){
  ///   dEdx_{Q}  reconstructed mean number of  electrons xgain
  ///     thr  = threshold in terms of the number of electrons
  ///     dEdxMode = algorithm to deal with trhesold values replacing

  //
  Double_t amp[159];
  Int_t nBellow=0;
  //
  Double_t minAbove=-1;
  for (Int_t i=0;i<fN;i++){ 
    AliTPCclusterFast * cluster = ( AliTPCclusterFast *)((*fCl)[i]);
    Double_t clQ= cluster->fQtot;
    if (clQ<thr) {
      nBellow++;
      continue;
    }
    if (minAbove<0) minAbove=clQ;
    if (minAbove>clQ) minAbove=clQ;
  }
  //
  if (dEdxMode==-1) return Double_t(nBellow)/Double_t(fN);

  for (Int_t i=0;i<fN;i++){ 
    AliTPCclusterFast * cluster = ( AliTPCclusterFast *)((*fCl)[i]);
    Double_t clQ= cluster->fQtot;
    //
    if (dEdxMode==0)  amp[i]=clQ;              // dEdxMode0 - not threshold  - keep default
    //
    //
    if (dEdxMode==1 && clQ>thr) amp[i]=clQ;    // dEdxMode1 - skip if bellow 
    if (dEdxMode==1 && clQ<thr) amp[i]=0;      // dEdxMode1 - skip if bellow 
    //
    //
    if (dEdxMode==2 && clQ>thr) amp[i]=clQ;    // dEdxMode2 - use 0 if below
    if (dEdxMode==2 && clQ<thr) amp[i]=0;      // dEdxMode2 - use 0 if below
    //
    //
    if (dEdxMode==3)  amp[i]=(clQ>thr)?clQ:thr; // dEdxMode3 - use thr if below
    if (dEdxMode==4)  amp[i]=(clQ>thr)?clQ:minAbove; // dEdxMode4 -  use minimal above threshold if bellow thr
  }
  return CookdEdx(fN,amp,f0,f1, dEdxMode);
}





Double_t   AliTPCtrackFast::CookdEdxDtot(Double_t f0,Float_t f1, Float_t gain,Float_t thr, Float_t noise, Bool_t doCorr, Int_t dEdxMode){
  /// total charge in the cluster (sum of the pad x time matrix ), hits were digitized before, but additional
  /// actions can be specified by switches  // dEdx_{Qtot}

  Double_t amp[159];
  Double_t minAmp=-1;
  //
  for (Int_t i=0;i<fN;i++){ 
    AliTPCclusterFast * cluster = ( AliTPCclusterFast *)((*fCl)[i]);
    Float_t camp = 0;
    if (dEdxMode==0) camp = cluster->GetQtot(gain,0,noise);
    else
      camp = cluster->GetQtot(gain,thr,noise);
    Float_t corr =  1;
    if (doCorr) corr = cluster->GetQtotCorr(0.5,0.5,gain,thr);
    camp/=corr;
    amp[i]=camp;
    if (camp>0){
      if (minAmp <0) minAmp=camp;
      if (minAmp >camp) minAmp=camp;
    }
  }
  if (dEdxMode==3) for (Int_t i=0;i<fN;i++) if (amp[i]<=0) amp[i]=thr;
  if (dEdxMode==4) for (Int_t i=0;i<fN;i++) if (amp[i]<=0) amp[i]=minAmp;
  return CookdEdx(fN,amp,f0,f1, dEdxMode);
}



Double_t   AliTPCtrackFast::CookdEdxDmax(Double_t f0,Float_t f1, Float_t gain,Float_t thr, Float_t noise, Bool_t doCorr, Int_t dEdxMode){
  /// maximal charge in the cluster (maximal amplitude in the digit matrix), hits were digitized before,
  /// but additional actions can be specified by switches

  Double_t amp[159];
  Double_t minAmp=-1;
  //
  for (Int_t i=0;i<fN;i++){ 
    AliTPCclusterFast * cluster = ( AliTPCclusterFast *)((*fCl)[i]);
    Float_t camp = 0;
    if (dEdxMode==0) camp =  cluster->GetQmax(gain,0,noise);
    else
      camp =  cluster->GetQmax(gain,thr,noise);
    Float_t corr =  1;
    if (doCorr) corr = cluster->GetQmaxCorr(0.5,0.5);
    camp/=corr;
    amp[i]=camp;
    if (camp>0){
      if (minAmp <0) minAmp=camp;
      if (minAmp >camp) minAmp=camp;
    }
  }
  if (dEdxMode==3) for (Int_t i=0;i<fN;i++) if (amp[i]<=0) amp[i]=thr;
  if (dEdxMode==4) for (Int_t i=0;i<fN;i++) if (amp[i]<=0) amp[i]=minAmp;
  return CookdEdx(fN,amp,f0,f1, dEdxMode);
}


Double_t  AliTPCtrackFast::CookdEdx(Int_t npoints, Double_t *amp,Double_t f0,Float_t f1, Int_t dEdxMode){
  /// Calculate truncated mean
  ///   npoints   - number of points in array
  ///   amp       - array with points
  ///   f0-f1     - truncation range
  ///   dEdxMode      - specify handling of the 0 clusters, actual handling - filling of amplitude defiend in algorithm above
  ///      dEdxMode =  0     - accept everything
  ///      dEdxMode =  1     - do not count 0 amplitudes
  ///      dEdxMode =  2     - use 0 amplitude as it is
  ///      dEdxMode =  3     - use amplitude as it is (in above function amp. replace by the thr)
  ///      dEdxMode =  4     - use amplitude as it is (in above function amp. replace by the minimal amplitude)

  //
  // 0. sorted the array of amplitudes
  //
  Int_t index[159];
  TMath::Sort(npoints,amp,index,kFALSE);
  //
  // 1.) Calculate truncated mean from the selected range of the array (ranking statistic )
  //     dependening on the dEdxMode 0 amplitude can be skipped
  Float_t sum0=0, sum1=0,sum2=0;
  Int_t   accepted=0;
  Int_t above=0;
  for (Int_t i=0;i<npoints;i++) if  (amp[index[i]]>0) above++;

  for (Int_t i=0;i<npoints;i++){
    //
    if (dEdxMode==1 && amp[index[i]]==0) {
      continue;
    }
    if (accepted<npoints*f0) continue;
    if (accepted>npoints*f1) continue;
    sum0++;
    sum1+= amp[index[i]];
    sum2+= amp[index[i]];
    accepted++;
  }
  if (dEdxMode==-1) return 1-Double_t(above)/Double_t(npoints);
  if (sum0<=0) return 0;
  return sum1/sum0;
}

void AliTPCtrackFast::Simul(const char* fname, Int_t ntracks, Double_t diffFactor, Bool_t simulOverlap){
  ///
  /// simulation is done in bin size:
  /// values has to be rescaled:
  /// bin size = (R=0.75, 1, 1.5  cm, RPhi=0.4 or 0.6 cm, Z=0.26 cm (10 MHz),0.52 cm (5 MHz))   
  ///
  /// diffusion is ~ 0.22 cm per sqrt(m) ~ 0.35 cm for full drift
  /// Simulation is done in 10 MHz scenario
  ///    5 MHz is obtained skipping the every second sample
  ///

  paramSR.Update(); // intitialize TPC parameters used later in the TOY simulation

  AliTPCtrackFast fast;
  TTreeSRedirector *pcstream = new TTreeSRedirector(fname,"recreate");
  for (Int_t itr=0; itr<ntracks; itr++){
    //
    fast.fMNprim=(13.+100*gRandom->Rndm());
    if (gRandom->Rndm()>0.5) fast.fMNprim=1./(0.00001+gRandom->Rndm()*0.1);

    fast.fDiff =0.22*(1+0.2*(gRandom->Rndm()-0.5));     // diffusion in mm/sqrt(cm) 0.22 cm/sqrt(m) +-10%       
    //    fast.fDiffLong =  fast.fDiff*0.6/1.;
    fast.fDiffLong =  fast.fDiff*diffFactor/1.;   
    //
    fast.fTY=gRandom->Gaus(0,0.01); //  cm vertex spread for primaries
    fast.fTZ=gRandom->Gaus(0,7);    //  7 cm vertex spread in z 
    if (gRandom->Rndm()>0.5) {
      fast.fTY=gRandom->Gaus(0,20);  //  20 cm vertex spread for secondaries - generate for half of statistic
      fast.fTZ=gRandom->Gaus(0,20);  //  20 cm vertex spread for secondaries - generate one half of statistic
    }

    fast.fTAngleY   = 4.0*(gRandom->Rndm()-0.5);
    if (gRandom->Rndm()<0.2)  fast.fTAngleY   = (gRandom->Rndm()-0.5)*TMath::Pi()/9.;    // admixture of high momenta tracks perpendicular to pad row +-10 degrees
    fast.fTAngleZ   = 3.0*(gRandom->Rndm()-0.5);  // track range as acceptabel for TPC
    fast.fN  = 159;

    if (simulOverlap){ //
      fast.fBOverlap=kTRUE;        // flag generate overlap track
      fast.fMNprimOverlap=1./(0.00001+gRandom->Rndm()*0.1);        // mean number of primary electrons for overlap track - flat in 1/Q
      //                                                           // better to get "realistic" q distribution
      fast.fDAngleYOverlap=(gRandom->Rndm()-0.5)*20./fast.fN;      // y angle - tan(y) for overlap track - for full lenght +-10 bins  
      fast.fDAngleZOverlap=(gRandom->Rndm()-0.5)*20./fast.fN;      // z angle - tan z for overlap track
      fast.fDYOverlap=((gRandom->Rndm()-0.5)*5.);                   // delta y position of overlap track at row 0 
      fast.fDZOverlap=((gRandom->Rndm()-0.5)*5.);                   // delta z position of overlap track at row 0 
    }
    fast.MakeTrack();
    if (itr%100==0) printf("%d\n",itr);
    (*pcstream)<<"simulTrack"<<
      "tr.="<<&fast<<
      "\n";
  }
  fast.Write("track");
  delete pcstream;
}



AliTPCclusterFast::AliTPCclusterFast():
  fOverlapCluster(0),
  fPRFRMS(0.5),     // pad respons function width
  fTRFRMS(0.5),      // time responsefunction width
  fRawDigits()
{
  ///
  fDigits.ResizeTo(5,7);
  fRawDigits.ResizeTo(5,7);
}

void AliTPCclusterFast::Init(){
  /// reset all counters

  const Int_t knMax=10000;
  fMNprim=0;     // mean number of primary electrons
  //                   //electrons part input
  fNprim=0;      // mean number of primary electrons
  fNtot=0;       // total number of  electrons
  fQtot=0;       // total charge - Gas gain flucuation taken into account
  fYCenterBin=0;
  fZCenterBin=0;
  fPadRow=0;
  //
  fPosY.ResizeTo(knMax);
  fPosZ.ResizeTo(knMax);
  fGain.ResizeTo(knMax);
  fSec.ResizeTo(knMax);
  fStatY.ResizeTo(3);
  fStatZ.ResizeTo(3);
  for (Int_t i=0; i<knMax; i++){
    fPosY[i]=0;
    fPosZ[i]=0;
    fGain[i]=0;
    fSec[i]=0;
  }
  fDiff=0;       ///< diffusion sigma
  fDiffLong=0;   ///< diffusion sigma longitudinal direction
  fY=0;          ///< y ideal position - center bin
  fZ=0;          ///< z ideal position - center bin
  fYCenterBin=0; ///< y center bin  
  fZCenterBin=0; ///< z center bin  
  fAngleY=0;     ///< y angle - tan(y)
  fAngleZ=0;     ///< z angle - tan z
}



AliTPCclusterFast::~AliTPCclusterFast(){
}


void AliTPCclusterFast::SetParam(Float_t mnprim, Float_t diff,  Float_t diffL, Int_t padrow, Float_t y, Float_t z, Float_t ky, Float_t kz, Float_t yCenter, Float_t zCenter){
  ///

  fMNprim = mnprim; fDiff = diff; fDiffLong=diffL;
  fY=y; fZ=z; 
  fAngleY=ky; fAngleZ=kz;
  fYCenterBin=yCenter;
  fZCenterBin=zCenter;
  fPadRow=padrow;
  
}
Double_t AliTPCclusterFast::GetNsec(){
  /// Generate number of secondary electrons
  /// copy of procedure implemented in geant

  const Double_t FPOT=20.77E-9, EEND=10E-6, EEXPO=2.2; // EEND1=1E-6;
  const Double_t XEXPO=-EEXPO+1, YEXPO=1/XEXPO;
  const Double_t W=20.77E-9;
  Float_t RAN = gRandom->Rndm();
  //Double_t edep = TMath::Power((TMath::Power(FPOT,XEXPO)*(1-RAN)+TMath::Power(EEND,XEXPO)*RAN),YEXPO);
  //edep = TMath::Min(edep, EEND);
  //return TMath::Nint(edep/W);
  return TMath::Nint(TMath::Power((TMath::Power(FPOT,XEXPO)*(1-RAN)+TMath::Power(EEND,XEXPO)*RAN),YEXPO)/W);
}

void AliTPCclusterFast::GenerElectrons(AliTPCclusterFast *cl0, AliTPCclusterFast *clm, AliTPCclusterFast *clp){
  ///
  //
  //
  const Int_t knMax=1000;
  cl0->fNprim = gRandom->Poisson(cl0->fMNprim);  //number of primary electrons
  // cl0->fNtot=0; //total number of electrons
  // cl0->fQtot=0; //total number of electrons after gain multiplification
  //
  Double_t sumQ=0;
  Double_t sumYQ=0;
  Double_t sumZQ=0;
  Double_t sumY2Q=0;
  Double_t sumZ2Q=0;
  //  for (Int_t i=0;i<knMax;i++){ 
  //  cl0->fSec[i]=0;
  //}
  for (Int_t iprim=0; iprim<cl0->fNprim;iprim++){   // loop over primary electrons
    Float_t dN   =  cl0->GetNsec();
    cl0->fSec[iprim]=dN;
    Double_t rc = (gRandom->Rndm()-0.5);             // primary electrons distributed randomly along pad row
    Double_t yc = cl0->fY+rc*cl0->fAngleY;           // primary electorns along stright line trajectory +-0.5 bin (pad-row) 
    Double_t zc = cl0->fZ+rc*cl0->fAngleZ;

    for (Int_t isec=0;isec<=dN;isec++){
      //
      //
      Double_t y = gRandom->Gaus(0,cl0->fDiff)+yc;
      Double_t z = gRandom->Gaus(0,cl0->fDiff*fgZSamplingFactor)+zc;
      Double_t r = gRandom->Gaus(0,cl0->fDiffLong)+rc;
      // choose pad row
      AliTPCclusterFast *cl=cl0;
      if (r<-0.5 &&clm) cl=clm;
      if (r>0.5 &&clp)  cl=clp;
      //
      Double_t gg = -TMath::Log(gRandom->Rndm());
      cl->fPosY[cl->fNtot]=y;
      cl->fPosZ[cl->fNtot]=z;
      cl->fGain[cl->fNtot]=gg;
      cl->fQtot+=gg;
      cl->fNtot++;
      //
      //
      if (cl->fNtot>=knMax) continue;
    }
    if (cl0->fNtot>=knMax) break;
  }
}

void AliTPCclusterFast::Digitize(){
  /// 1. Clear digits

  for (Int_t i=0; i<5;i++)
    for (Int_t j=0; j<7;j++){
      fDigits(i,j)=0;
    }
  //
  // Fill digits
  fStatY*=0;
  fStatZ*=0;
  for (Int_t iel = 0; iel<fNtot; iel++){
    Double_t gg=fGain[iel];
    Double_t y=fPosY[iel];
    Double_t z=fPosZ[iel];
    fStatY[0]+=gg;
    fStatY[1]+=gg*y;
    fStatY[2]+=gg*y*y;
    fStatZ[0]+=gg;
    fStatZ[1]+=gg*z;
    fStatZ[2]+=gg*z*z;
    for (Int_t di=-2; di<=2;di++)
      for (Int_t dj=-3; dj<=3;dj++){	
	Float_t fac = fPRF->Eval(di-y)*fTRF->Eval(dj-z);
	fac*=gg;
	fDigits(2+di,3+dj)+=fac;
      }
  }
  //
  //
  //
}


Double_t AliTPCclusterFast::GetQtot(Float_t gain, Float_t thr, Float_t noise, Bool_t brounding, Bool_t baddPedestal){
  ///

  Float_t sum =0;
  for (Int_t ip=0;ip<5;ip++){
    Float_t pedestal=gRandom->Rndm()-0.5; //pedestal offset different for each pad
    for (Int_t it=0;it<7;it++){
      Float_t amp = gain*fDigits(ip,it)+gRandom->Gaus()*noise;
      if (baddPedestal) amp+=pedestal;
      if (brounding) amp=TMath::Nint(amp);
      if (amp>thr) sum+=amp;
    } 
  }
  return sum;
}

Double_t AliTPCclusterFast::GetQmax(Float_t gain, Float_t thr, Float_t noise, Bool_t brounding, Bool_t baddPedestal){
  ///

  Float_t max =0;
  for (Int_t ip=0;ip<5;ip++){
    Float_t pedestal=gRandom->Rndm()-0.5; //pedestal offset different for each pad
    for (Int_t it=0;it<7;it++){
      Float_t amp = gain*fDigits(ip,it)+gRandom->Gaus()*noise;
      if (baddPedestal) amp+=pedestal;
      if (brounding) amp=TMath::Nint(amp);
      if (amp>max && amp>thr) max=amp;
    } 
  }
  return max;
}



Double_t  AliTPCclusterFast::GetQmaxCorr(Float_t rmsy0, Float_t rmsz0){
  /// Gaus distribution convolueted with rectangular
  /// Gaus width sy and sz is determined by RF width and diffusion
  /// Integral of Q is equal 1
  /// Q max is calculated at position fY,fX

  Double_t sy = TMath::Sqrt(rmsy0*rmsy0+fDiff*fDiff);
  Double_t sz = TMath::Sqrt(rmsz0*rmsz0+fDiff*fDiff); 
  return GaussConvolution(fY,fZ, fAngleY,fAngleZ,sy,sz);
}


Double_t  AliTPCclusterFast::GetQtotCorr(Float_t rmsy0, Float_t rmsz0, Float_t gain, Float_t thr){
  ///  Calculates the fraction of the charge over threshol to total charge
  ///  The response function

  Double_t sy = TMath::Sqrt(rmsy0*rmsy0+fDiff*fDiff);
  Double_t sz = TMath::Sqrt(rmsz0*rmsz0+fDiff*fDiff);
  Double_t sumAll=0,sumThr=0;
  Double_t qtot = GetQtot(gain,thr,0); // sum of signal over threshold
  Double_t qmax = GetQmax(gain,thr,0); // qmax
  //
  Double_t corr =1;
  Double_t qnorm=qtot;
  for (Int_t iter=0;iter<2;iter++){
    for (Int_t iy=-2;iy<=2;iy++)
      for (Int_t iz=-2;iz<=2;iz++){      
	Double_t val = GaussConvolution(fY-iy,fZ-iz, fAngleY,fAngleZ,sy,sz);
	Double_t qlocal =TMath::Nint(qnorm*val);
	if (qlocal>thr) sumThr+=qlocal;
	sumAll+=qlocal;
      }
    if (sumAll>0&&sumThr>0) corr=(sumThr)/sumAll;
    if (sumThr==0) corr=GetQmaxCorr(0.5,0.5);
    //corr = sumThr;
    if (corr>0) qnorm=qtot/corr;
    
  }
  return corr;
}

Double_t AliTPCclusterFast::GetClusterProperties(TVectorF &param, Int_t addOverlap, Float_t gain, Float_t thr, Float_t noise, Bool_t rounding, Bool_t addPedestal, Int_t skipSample){
  //
  // calculate area, COG, RMS and  skeewnes  cluster without any correction
  //
  Float_t sumW=0;
  Float_t sumYW=0;
  Float_t sumZW=0;
  Float_t sumY2W=0;
  Float_t sumZ2W=0;
  Float_t sumY3W=0;
  Float_t sumZ3W=0;
  for (Int_t iy=-2;iy<=2;iy++)
    for (Int_t jz=-2;jz<=2;jz+=1+skipSample){      
      Int_t iz=jz;
      if (skipSample>0) iz=jz/(1+skipSample);      
      Double_t val = gain*fDigits(iy+2,jz+3);
      if (addOverlap && fOverlapCluster){
	val+=gain*fOverlapCluster->fDigits(iy+2,jz+3);
      }      
      val+=noise*gRandom->Gaus();
      if (addPedestal) val+=gRandom->Rndm()-0.5;   // add random pedestal assuming mean bias value 0
      if (val>thr){
	sumW+=val;
	sumYW+=iy*val;
	sumZW+=iz*val;
	sumY2W+=iy*iy*val;
	sumZ2W+=iz*iz*val;
      }
    }
  if (sumW>0) {
    sumYW/=sumW;
    sumZW/=sumW;
    sumY2W/=sumW;
    sumZ2W/=sumW;
    if (skipSample>0){
      sumZW*=2.;
      sumZ2W*=4.;
    }
  }
  param[0]=sumW/gain;                        // total charge
  param[1]=sumYW;                            // mean y
  param[2]=sumZW;                            // mean z
  param[3]=TMath::Sqrt(sumY2W-sumYW*sumYW);  // rms y
  param[4]=TMath::Sqrt(sumZ2W-sumZW*sumZW);  // rms z
  //
  // 3 moment
  //
  for (Int_t iy=-2;iy<=2;iy++)
    for (Int_t iz=-2;iz<=2;iz++){      
      Double_t val = fDigits(iy+2,iz+3);
      if (addOverlap && fOverlapCluster){
	val+=fOverlapCluster->fDigits(iy+2,iz+3);
      }      
      val+=noise;
      if (addPedestal) val+=gRandom->Rndm()-0.5;   // add random pedestal assuming mean bias value 0
      if (val>thr){
	Double_t dy=iy-param[1];
	Double_t dz=iz-param[2];
	sumY3W+=dy*dy*dy*val;
	sumZ3W+=dz*dz*dz*val;
      }
    }
  if (sumW>0){
    param[5]=sumY3W/sumW;
    param[6]=sumZ3W/sumW;  
  }
}

Double_t  AliTPCclusterFast::GetExpectedRMS(Int_t dim){
  //
  // Idealized - 
  // Expected RMS of cluster in case of not threhsold effect and noise equal 0 in bin units
  // see:
  //    http://arxiv.org/pdf/physics/0306108.pdf
  //    formulas 26, 27 
  //          in the code - all variables in is bin size units  - not necessay to normalize to pad length, resp pad width)
  if (dim==1) return TMath::Sqrt(fPRFRMS*fPRFRMS+fDiff*fDiff+fAngleY*fAngleY/12.);
  if (dim==2) return TMath::Sqrt(fTRFRMS*fTRFRMS+fDiff*fDiff*fgZSamplingFactor*fgZSamplingFactor+fAngleZ*fAngleZ/12.);
  return 0;
}


Double_t AliTPCclusterFast::GetCOG(Int_t returnType, Int_t addOverlap, Float_t gain, Float_t thr, Float_t noise, Bool_t rounding, Bool_t addPedestal, Int_t skipSample){
  //
  //
  //
  TVectorF properties(20);
  GetClusterProperties(properties, addOverlap, gain, thr, noise,rounding,addPedestal,skipSample);
  return properties(returnType);
}

Double_t AliTPCclusterFast::GetCOGHit(Int_t returnType){
  //
  // Return COG using hit inforamtion
  //
  if (returnType==0) return fStatY[0];
  Float_t meanY=(fStatY[0]>0)?fStatY[1]/fStatY[0]:0;
  if (returnType==1) return meanY;   // meanY
  Float_t meanZ=(fStatZ[0]>0)?fStatZ[1]/fStatZ[0]:0;
  if (returnType==2) return meanZ;   // mean Z
  
  if (returnType==3) return (fStatY[0]>0) ? TMath::Sqrt(fStatY[2]/fStatY[0]-meanY*meanY):0;
  if (returnType==4) return (fStatZ[0]>0) ? TMath::Sqrt(fStatZ[2]/fStatZ[0]-meanZ*meanZ):0;
  return 0;
}



Double_t  AliTPCclusterFast::GaussConvolution(Double_t x0, Double_t x1, Double_t k0, Double_t k1, Double_t s0, Double_t s1){
  /// 2 D gaus convoluted with angular effect
  /// See in mathematica:
  /// Simplify[Integrate[Exp[-(x0-k0*xd)*(x0-k0*xd)/(2*s0*s0)-(x1-k1*xd)*(x1-k1*xd)/(2*s1*s1)]/(s0*s1),{xd,-1/2,1/2}]]
  ///
  /// TF1 f1("f1","AliTPCclusterFast::GaussConvolution(x,0,1,0,0.1,0.1)",-2,2)
  /// TF2 f2("f2","AliTPCclusterFast::GaussConvolution(x,y,1,1,0.1,0.1)",-2,2,-2,2)

  const Float_t kEpsilon = 0.0001;
  if ((TMath::Abs(k0)+TMath::Abs(k1))<kEpsilon*(s0+s1)){
    // small angular effect
    Double_t val = (TMath::Gaus(x0,0,s0)*TMath::Gaus(x1,0,s1))/(s0*s1*2.*TMath::Pi());
    return val;
  }
  
  Double_t sigma2 = k1*k1*s0*s0+k0*k0*s1*s1;
  Double_t exp0 = TMath::Exp(-(k1*x0-k0*x1)*(k1*x0-k0*x1)/(2*sigma2));
  //
  Double_t sigmaErf =  2*s0*s1*TMath::Sqrt(2*sigma2);			     
  Double_t erf0 = TMath::Erf( (k0*s1*s1*(k0-2*x0)+k1*s0*s0*(k1-2*x1))/sigmaErf);
  Double_t erf1 = TMath::Erf( (k0*s1*s1*(k0+2*x0)+k1*s0*s0*(k1+2*x1))/sigmaErf);
  Double_t norm = 1./TMath::Sqrt(sigma2);
  norm/=2.*TMath::Sqrt(2.*TMath::Pi());
  Double_t val  = norm*exp0*(erf0+erf1);
  return val;

}


Double_t  AliTPCclusterFast::GaussExpConvolution(Double_t x0, Double_t s0,Double_t t1){
 /// 2 D gaus convoluted with exponential
 /// Integral nomalized to 1
 /// See in mathematica:
 /// Simplify[Integrate[Exp[-(x0-x1)*(x0-x1)/(2*s0*s0)]*Exp[-x1*t1],{x1,0,Infinity}]]
 /// TF1 fgexp("fgexp","AliTPCclusterFast::GaussExpConvolution(x,0.5,1)",-2,2)

  Double_t exp1 = (s0*s0*t1-2*x0)*t1/2.;
  exp1 = TMath::Exp(exp1);
  Double_t erf = 1+TMath::Erf((-s0*s0*t1+x0)/(s0*TMath::Sqrt(2.)));
  Double_t val = exp1*erf;
  val *=t1/(2.);
  return val;

}

Double_t  Gamma4(Double_t *x, Double_t *param){
  /// Gamma 4 Time response function of ALTRO

  if (x[0]<0) return 0;
  Double_t g1 = TMath::Exp(-4.*x[0]/param[1]);
  Double_t g2 = TMath::Power(x[0]/param[1],4);
  return param[0]*g1*g2;
}
 

Double_t  AliTPCclusterFast::Gamma4(Double_t x, Double_t p0, Double_t p1){
  /// Gamma 4 Time response function of ALTRO
  /// Gamma 4 Time response function of ALTRO

  if (x<0) return 0;
  Double_t g1 = TMath::Exp(-4.*x/p1);
  Double_t g2 = TMath::Power(x/p1,4);
  return p0*g1*g2;
}

Double_t  AliTPCclusterFast::Gamma4Norm(Double_t x){
  /// Gamma 4 Time response function of ALTRO with hardwired normalization parameters
  ///   X in bin unit=100 ns
  ///    Maximum of the TRF ==1  - sum of all time bins ~ 1.9
  ///
  return AliTPCclusterFast::Gamma4((x+1.98094355865156)*0.1,55,0.160);
}
 




Double_t AliTPCclusterFast::GEMPRF(Double_t x, Double_t sigma){
  //
  // GEM PRF response function aprroximated as integral of gaussian
  // sigma of gaussian is the diffusion of electrons within GEM layers
  // in 0.5 cm of the drift lenght + O(0.100)mm GEM pitch ==> 0.2 mm 
  if (x<0) return (1+TMath::Erf((x+0.5)/sigma))*0.5;
  if (x>0) return (1+TMath::Erf(-(x-0.5)/sigma))*0.5;
}


Double_t  AliTPCclusterFast::GaussGamma4(Double_t x, Double_t s0, Double_t p1){
  /// Gamma 4 Time response function of ALTRO convoluted with Gauss
  /// Simplify[Integrate[Exp[-(x0-x1)*(x0-x1)/(2*s0*s0)]*Exp[-4*x1/p1]*(x/p1)^4/s0,{x1,0,Infinity}]]
  /// TF1 fgg4("fgg4","AliTPCclusterFast::GaussGamma4(x,0.5,0.5)",-2,2)

  Double_t exp1 = (8*s0*s0-4.*p1*x)/(p1*p1);
  exp1 = TMath::Exp(exp1);
  Double_t erf1 = 1+TMath::Erf((-4*s0/p1+x/s0)/TMath::Sqrt(2));
  //  Double_t xp14 = TMath::Power(TMath::Abs((x/p1)),4);
  return exp1*erf1;

 
}
Float_t AliTPCclusterFast::GetCOGUnfolded(Int_t returnType, Int_t addOverlap, Float_t gain, Float_t thr, Float_t noise, Bool_t rounding, Bool_t addPedestal,Int_t skipSample){
  //
  //
  //
  Float_t meani, meanj, sumu,overlap;
  UnfoldCluster(meani, meanj, sumu,overlap,addOverlap, gain,thr,noise,rounding,addPedestal,skipSample);
  if (returnType==-1) return overlap;
  if (returnType==0) return sumu;
  if (returnType==1) return meani;
  if (returnType==2) return meanj;
}


Bool_t  AliTPCclusterFast::MakeDigitization(Int_t addOverlap, Float_t gain, Float_t thr, Float_t noise, Bool_t rounding, Bool_t addPedestal,Int_t skipSample){
  //
  // Conversion ideal digits array to digits 
  // In addition local max bin calculated 
  // 

  //
  // 1.) Make digitization and find the maximal bin 
  //     store position if maximal bin in the array
  //       in case of rounding can happen we have the same values in maxima
  //       to avoid bias we should take one position randomly
  //       this bias is present also in the OFFLINE clusterer - can be seen like assymetry in the peak position
  //       Numericaly effect is significant at Qmax~ 10 it is about 0.025 bin size
 
  Float_t maxValue=0;
  const Float_t kEpsilon=0.000001;
  for (Int_t i=0; i<5; i++){
    for (Int_t j=0; j<7; j++){
      Float_t value = fDigits(i,j);
      if (addOverlap && fOverlapCluster) value+=fOverlapCluster->fDigits(i,j);
      value*=gain;
      Float_t cNoise=gRandom->Gaus();
      value+=noise*cNoise;
      if (addPedestal) value+=gRandom->Rndm()-0.5;
      if (rounding) value=TMath::Nint(value);
      if (value<thr) value=0;      
      if (skipSample==0) {
	fRawDigits(i,j)=value;
	if ((value+kEpsilon*cNoise)>=maxValue  && TMath::Abs(i-2)<=1 &&  TMath::Abs(j-3)<=(1+skipSample)){  // check position of local maxima
	  maxValue=value+kEpsilon*cNoise; // in case of equal values alogn peak chose one randomly
	  fYMaxBin=i-2;
	  fZMaxBin=j-3;
	}
      }
      if (skipSample==1 && (j%2)==1  && TMath::Abs(i-2)<=1 &&  TMath::Abs(j-3)<=(1+skipSample) ) {    // check position of local maxima
	fRawDigits(i,2+j/2)=value;      
	if (value+kEpsilon*cNoise>=maxValue){
	  maxValue=value+kEpsilon*cNoise;  // in case of equal values alogn peak chose one randomly
	  fYMaxBin=i-2;
	  fZMaxBin=2+j/2-3;
	}
      }
    }
  }
  //
  // shift digits  - local maxima should be at predeefinde position y=2 z=3
  //
  if (fYMaxBin!=0 ||  fZMaxBin!=0){  //shift digits  - local maxima should be at predeefinde position y=2 z=3
    TMatrixF tmpDigits(5,7);
    for (Int_t i=0; i<5; i++)
      for (Int_t j=0; j<7; j++)	
	if ((i-fYMaxBin)>=0 && (i-fYMaxBin)<5 && (j-fZMaxBin)>=0 && (j-fZMaxBin)<7) {
	  tmpDigits(i-fYMaxBin,j-fZMaxBin)= fRawDigits(i, j);
	}
    fRawDigits= tmpDigits;          
  }
  
  static Int_t dumpCounter=0;
  dumpCounter++;
  if (fgDebugLevel>0 && (dumpCounter%fgDebugLevel)==0){
    ::Info("AliTPCclusterFast::MakeDigitization","fYMaxBin\t%d fZMaxBin\t%d",fYMaxBin,fZMaxBin);    
  }
  return kTRUE;
  /*
   */
}



Float_t  AliTPCclusterFast::UnfoldCluster(Float_t & meani, Float_t & meanj,  Float_t & sumu, Float_t & overlap, Int_t addOverlap, Float_t gain, Float_t thr, Float_t noise, Bool_t rounding, Bool_t addPedestal,Int_t skipSample){
  //
  // Unforling as used in the  AliTPCclusterer::UnfoldCluster - (described in CHEP paper ...)
  //  - we are not using second local maxima as we assume we know the cluster shape
  //  - in case of absence of this information local maxima usage to be considered
  //  
  Float_t sum3i[7] = {0,0,0,0,0,0,0};
  Float_t sum3j[7] = {0,0,0,0,0,0,0};
  meani=0;
  meanj=0;
  sumu=0;
  overlap=0;
  //
  //
  //
  MakeDigitization(addOverlap, gain, thr, noise,rounding, addPedestal, skipSample);

  
  for (Int_t k =0;k<7;k++) // sequence array
    for (Int_t l = -1; l<=1;l++){  // +- 1 bin loop
      if (k>0&&k<6) sum3i[k]+=fRawDigits(k-1,l+3);  //  i- y direction
      sum3j[k]+=fRawDigits(l+2,k);
    }
  if (sum3i[3]<=3 || sum3j[3]<=3) return 0;


  Float_t mratio[3][3]={{1,1,1},{1,1,1},{1,1,1}};
  //
  //unfold  y 
  Float_t sum3wi    = 0;  //charge minus overlap
  Float_t sum3wio   = 0;  //full charge
  Float_t sum3iw    = 0;  //sum for mean value
  for (Int_t dk=-1;dk<=1;dk++){
    sum3wio+=sum3i[dk+3];
    if (dk==0){
      sum3wi+=sum3i[dk+3];     
    }
    else{
      Float_t ratio =1;
      if (  ( ((sum3i[dk+3]+3)/(sum3i[3]-3))+1 < (sum3i[2*dk+3]-3)/(sum3i[dk+3]+3))||
	    (sum3i[dk+3]<=sum3i[2*dk+3] && sum3i[dk+3]>2 )){
	Float_t xm2 = sum3i[-dk+3];
	Float_t xm1 = sum3i[+3];
	Float_t x1  = sum3i[2*dk+3];
	Float_t x2  = sum3i[3*dk+3]; 	
	Float_t w11   = TMath::Max((Float_t)(4.*xm1-xm2),(Float_t)0.000001);	  
	Float_t w12   = TMath::Max((Float_t)(4 *x1 -x2),(Float_t)0.);
	ratio = w11/(w11+w12);	 
	for (Int_t dl=-1;dl<=1;dl++)
	  mratio[dk+1][dl+1] *= ratio;
      }
      Float_t amp = sum3i[dk+3]*ratio;
      sum3wi+=amp;
      sum3iw+= dk*amp;      
    }
  }
  meani = sum3iw/sum3wi;
  Float_t overlapi = (sum3wio-sum3wi)/sum3wio;
  //unfold  z 
  Float_t sum3wj    = 0;  //charge minus overlap
  Float_t sum3wjo   = 0;  //full charge
  Float_t sum3jw    = 0;  //sum for mean value
  for (Int_t dk=-1;dk<=1;dk++){
    sum3wjo+=sum3j[dk+3];
    if (dk==0){
      sum3wj+=sum3j[dk+3];     
    }
    else{
      Float_t ratio =1;
      if ( ( ((sum3j[dk+3]+3)/(sum3j[3]-3))+1 < (sum3j[2*dk+3]-3)/(sum3j[dk+3]+3)) ||
	   (sum3j[dk+3]<=sum3j[2*dk+3] && sum3j[dk+3]>2)){
	Float_t xm2 = sum3j[-dk+3];
	Float_t xm1 = sum3j[+3];
	Float_t x1  = sum3j[2*dk+3];
	Float_t x2  = sum3j[3*dk+3]; 	
	Float_t w11   = TMath::Max((Float_t)(4.*xm1-xm2),(Float_t)0.000001);	  
	Float_t w12   = TMath::Max((Float_t)(4 *x1 -x2),(Float_t)0.);
	ratio = w11/(w11+w12);	 
	for (Int_t dl=-1;dl<=1;dl++)
	  mratio[dl+1][dk+1] *= ratio;
      }
      Float_t amp = sum3j[dk+3]*ratio;
      sum3wj+=amp;
      sum3jw+= dk*amp;      
    }
  }
  meanj = sum3jw/sum3wj;
  Float_t overlapj = (sum3wjo-sum3wj)/sum3wjo;  
  overlap = Int_t(100*TMath::Max(overlapi,overlapj)+3);  
  sumu = (sum3wj+sum3wi)/2.;
}

Bool_t AliTPCclusterFast::UnitTest(){
  //
  //
  //
  ::Info("AliTPCclusterFast::UnitTest","Test BEGIN");
  //
  // Test1. local Max position
  testTree->SetAlias("IsMaxOK","((fCl.GetRawDigit(2,3)>=fCl.GetRawDigit(1,3))&&(fCl.GetRawDigit(2,3)>=fCl.GetRawDigit(4,3)))!=0");
  testTree->Draw("IsMaxOK==0","fCl.GetCOGUnfolded(1, 0, 1,  0.0 ,  0.0,  0, 0, 0)!=0","goff",1000);
  Double_t mean=TMath::Mean(testTree->GetSelectedRows(),testTree->GetV1());
  if (TMath::Abs(mean)<0.001){
    ::Info("AliTPCclusterFast::UnitTest","MaxTest OK");
  }else{
    ::Error("AliTPCclusterFast::UnitTest","MaxTest FAILED");
  }
}

void AliTPCclusterFast::SetMetadata(TTree* tree){
  //
  //
  //
  tree->SetAlias("deltaYHit","(fCl.GetCOGHit(1)-fY-0)");    // ideal resolution using pefect readout"
  tree->SetAlias("errYHit","sqrt(2.*fCl.fDiff**2/fCl.fQtot+(2*fCl.fAngleY**2)/(12.*sqrt(fCl.fNprim)))"); // resoulution for ideal response detector
  tree->SetAlias("errYMWPC","sqrt(0.2**2+2.*fCl.fDiff**2/fCl.fQtot+(2*fCl.fAngleY**2)/(12.*sqrt(fCl.fNprim)))"); // approximate resolution for the MWPC detector 

  tree->SetAlias("deltaYUnfoldedDefault","(GetCOGUnfolded(1, 1, 0.8,   2,  0.6,  1, 1, 0)+fCl.GetYMaxBin()-fY)");
  tree->SetAlias("deltaYNoOverlapDefault","(fCl.GetCOG(1, 0, 0.8,   2,  0.6,  1, 1, 0)-fY-0)");
  tree->SetAlias("deltaYOverlapDefault","(fCl.GetCOG(1, 1, 0.8,   2,  0.6,  1, 1, 0)-fY-0)");
  tree->SetAlias("padrow","Iteration$");

  tree->SetAlias("IROC","padrow<63");
  tree->SetAlias("OROCMedium","padrow>63&&padrow<127");
  tree->SetAlias("OROCLong","padrow>127");  
  tree->SetAlias("padLength","((IROC)*0.75+(OROCMedium)*1+(OROCLong*1.5))");
  tree->SetAlias("padWidth","((IROC)*0.4+(padrow>63)*0.6)");
  tree->SetAlias("driftLength","250-abs(fZCenterBin*0.25)");
  
  //

  //
  //
  //
  TStatToolkit::AddMetadata(tree,"deltaYHit.AxisTitle","#Delta_{r#phi} (bin)");
  TStatToolkit::AddMetadata(tree,"deltaYHit.Title","(fCl.GetCOGHit(1)-fY-0)");
  TStatToolkit::AddMetadata(tree,"deltaYHit.Legend","#Delta_{rphi} hit COG");
  TStatToolkit::AddMetadata(tree,"deltaYHit.Comment","#Delta_{rphi} for COG estimator using hits ");
  //
  TStatToolkit::AddMetadata(tree,"errYHit.AxisTitle","#sigma_{r#phi} (bin)");
  TStatToolkit::AddMetadata(tree,"errYHit.Title","(fCl.GetCOGHit(1)-fY-0)");
  TStatToolkit::AddMetadata(tree,"errYHit.Legend","#sigma_{rphi} hit COG");
  TStatToolkit::AddMetadata(tree,"rrYHit.Comment","#sigma_{rphi} for COG estimator using hits ");
  //
  //
  TStatToolkit::AddMetadata(tree,"deltaYNoOverlapDefault.AxisTitle","#Delta_{r#phi} (bin)");
  TStatToolkit::AddMetadata(tree,"deltaYNoOverlapDefault.Title","(fCl.GetCOG(1, 0, 0.8,   2,  0.6,  1, 1, 0)-fY)");
  TStatToolkit::AddMetadata(tree,"deltaYNoOverlapDefault.Legend","#Delta_{rphi} no overlap (10 MHz)");
  TStatToolkit::AddMetadata(tree,"deltaYNoOverlapDefault.Comment","#Delta_{rphi} no overlap (10 MHz) \n (fCl.GetCOG(1, 0, 0.8,   2,  0.6,  1, 1, 0)-fY)");
  //
  TStatToolkit::AddMetadata(tree,"deltaYUnfoldedDefault.AxisTitle","#Delta_{r#phi} (bin)");
  TStatToolkit::AddMetadata(tree,"deltaYUnfoldedDefault.Title","(GetCOGUnfolded(1, 1, 0.8,   2,  0.6,  1, 1, 0)+fCl.GetYMaxBin()-fY)");
  TStatToolkit::AddMetadata(tree,"deltaYUnfoldedDefault.Legend","#Delta_{rphi} unfolded  (10 MHz)");
  TStatToolkit::AddMetadata(tree,"deltaYUnfoldedDefault.Comment","#Delta_{rphi} unfolded  (10 MHz). Default digitization parameters. \n Sampling rate 10 MHz.\n Alias: (GetCOGUnfolded(1, 1, 0.8,   2,  0.6,  1, 1, 0)+fCl.GetYMaxBin()-fY)");
  //
  TStatToolkit::AddMetadata(tree,"fY.AxisTitle","r#phi (bin)");
  TStatToolkit::AddMetadata(tree,"fY.Title","fY");
  TStatToolkit::AddMetadata(tree,"fY.Legend","MC ideal r#phi position");
  TStatToolkit::AddMetadata(tree,"fY.Comment","MC ideal r#phi position of cluster");
  //
  TStatToolkit::AddMetadata(tree,"padrow.AxisTitle","pad row (bin)");
  TStatToolkit::AddMetadata(tree,"padrow.Title","padrow");
  TStatToolkit::AddMetadata(tree,"padrow.Legend","pad row (bin)");
  TStatToolkit::AddMetadata(tree,"padrow.Comment","Cluster - pad row position ");
  //
  
}

// Analytical sollution only in 1D - too long expression
// Simplify[Integrate[Exp[-(x0-(x1-k*x2))*(x0-(x1-k*x2))/(2*s0*s0)]*Exp[-(x1*t1-k*x2)],{x2,-1,1}]] 
//
//
// No analytical solution
// 
//Simplify[Integrate[Exp[-(x0-k0*xd)*(x0-k0*xd)/(2*s0*s0)-(x1-xt-k1*xd)*(x1-xt-k1*xd)/(2*s1*s1)]*Exp[-kt*xt]/(s0*s1),{xd,-1/2,1/2},{xt,0,Infinity}]]
