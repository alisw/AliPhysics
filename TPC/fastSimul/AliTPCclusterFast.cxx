/// \class AliTPCclusterFast
///
/// ~~~
/// gSystem->Load("libSTAT");
/// .x ~/NimStyle.C
///
/// .L $ALICE_ROOT/TPC/fastSimul/AliTPCclusterFast.cxx+
///
/// AliTPCclusterFast::fPRF = new TF1("fprf","gausn",-5,5);
/// AliTPCclusterFast::fTRF = new TF1("ftrf","gausn",-5,5);
/// AliTPCclusterFast::fPRF->SetParameters(1,0,0.5);
/// AliTPCclusterFast::fTRF->SetParameters(1,0,0.5);
///
/// AliTPCtrackFast::Simul("trackerSimul.root",100);
/// // AliTPCclusterFast::Simul("cluterSimul.root",20000);
/// ~~~
///
///  Modifications to add:
///  1. modigy mode ==> dEdxMode
///  2. Create hardware setup class
///     (fNoise, fGain, fBRounding, fBAddpedestal, ....)
///  3. Create arrays of registered hardware setups
///  4. Extend on the fly functions to use registered hardware setups, identified by ID.
///      hwMode

#include "TObject.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TClonesArray.h"
#include "TTreeStream.h"
#include "TGrid.h"

class AliTPCclusterFast: public TObject {
public:
  AliTPCclusterFast();
  void Init();
  
  virtual ~AliTPCclusterFast();
  void SetParam(Float_t mnprim, Float_t diff, Float_t diffL, Float_t y, Float_t z, Float_t ky, Float_t kz);
  static void GenerElectrons(AliTPCclusterFast *cl0, AliTPCclusterFast *clm, AliTPCclusterFast *clp);
  void Digitize();
  Double_t GetQtot(Float_t gain,Float_t thr, Float_t noise, Bool_t rounding=kTRUE, Bool_t addPedestal=kTRUE);
  Double_t GetQmax(Float_t gain,Float_t thr, Float_t noise, Bool_t rounding=kTRUE, Bool_t addPedestal=kTRUE);
  Double_t GetQmaxCorr(Float_t rmsy0, Float_t rmsz0);
  Double_t GetQtotCorr(Float_t rmsy0, Float_t rmsz0, Float_t gain, Float_t thr);
  
  Double_t GetNsec();
  //static void Simul(const char* simul, Int_t npoints);
  static Double_t GaussConvolution(Double_t x0, Double_t x1, Double_t k0, Double_t k1, Double_t s0, Double_t s1);
  static Double_t GaussExpConvolution(Double_t x0, Double_t s0,Double_t t1);
  static Double_t GaussGamma4(Double_t x, Double_t s0, Double_t p1);
  static Double_t Gamma4(Double_t x, Double_t p0, Double_t p1);
public:
  Float_t fMNprim;     ///< mean number of primary electrons
  //                   //electrons part input
  Int_t   fNprim;      ///< mean number of primary electrons
  Int_t   fNtot;       ///< total number of  electrons
  Float_t fQtot;       ///< total charge - Gas gain flucuation taken into account
  //
  Float_t fDiff;       ///< diffusion sigma
  Float_t fDiffLong;       ///< diffusion sigma longitudinal direction
  Float_t fY;          ///< y position
  Float_t fZ;          ///< z postion
  Float_t fAngleY;     ///< y angle - tan(y)
  Float_t fAngleZ;     ///< z angle - tan z
  //
  //
  //                   // electron part simul
  TVectorD fSec;       //!< number of secondary electrons
  TVectorD fPosY;      //!< position y for each electron
  TVectorD fPosZ;      //!< position z for each electron
  TVectorD fGain;      //!< gg for each electron
  //
  TVectorD fStatY;     //!< stat Y
  TVectorD fStatZ;     //!< stat Y
  //
  // digitization part
  //
  TMatrixD fDigits;    ///< response matrix
  static TF1* fPRF;    ///< Pad response
  static TF1* fTRF;    ///< Time response function
  ClassDef(AliTPCclusterFast,1)  // container for
};


class AliTPCtrackFast: public TObject {
public:
  AliTPCtrackFast();
  void Add(AliTPCtrackFast &track2);
  void MakeTrack();
  static void Simul(const char* simul, Int_t ntracks, Double_t diff);
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
  Float_t fMNprim;     ///< mean number of primary electrons
  Float_t fAngleY;     ///< y angle - tan(y)
  Float_t fAngleZ;     ///< z angle - tan z
  Float_t fDiff;       ///< diffusion
  Float_t fDiffLong;       ///< diffusion sigma longitudinal direction
  Int_t   fN;          ///< number of clusters
  TClonesArray *fCl;   ///< array of clusters
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


AliTPCtrackFast::AliTPCtrackFast():
  TObject(),
  fMNprim(0),
  fAngleY(0),
  fAngleZ(0),
  fN(0),
  fCl(0),
  fInit(kFALSE)
{
  ///

}

void AliTPCtrackFast::Add(AliTPCtrackFast &track2){
  if (!track2.fInit) return;
  
}



void AliTPCtrackFast::MakeTrack(){
  ///

  if (!fCl) fCl = new TClonesArray("AliTPCclusterFast",160);
  //
  // 0.) Init data structure
  //
  for (Int_t i=0;i<fN;i++){
    AliTPCclusterFast * cluster = (AliTPCclusterFast*) fCl->UncheckedAt(i);
    if (!cluster) cluster =   new ((*fCl)[i]) AliTPCclusterFast;
    cluster->Init();
  }
  //
  // 1.) Create hits - with crosstalk diffusion
  //
  for (Int_t i=0;i<fN;i++){
    Double_t tY = i*fAngleY;
    Double_t tZ = i*fAngleZ;
    AliTPCclusterFast * cluster = (AliTPCclusterFast*) fCl->UncheckedAt(i);
    AliTPCclusterFast * clusterm = (AliTPCclusterFast*) fCl->UncheckedAt(TMath::Max(i-1,0));
    AliTPCclusterFast * clusterp = (AliTPCclusterFast*) fCl->UncheckedAt(TMath::Min(i+1,159));
    if (!cluster) cluster =   new ((*fCl)[i]) AliTPCclusterFast;
    //
    Double_t posY = tY-TMath::Nint(tY);
    Double_t posZ = tZ-TMath::Nint(tZ);
    cluster->SetParam(fMNprim,fDiff, fDiffLong, posY,posZ,fAngleY,fAngleZ); 
    //
    cluster->GenerElectrons(cluster, clusterm, clusterp);
  }
  //
  // 2.) make digitization
  //
  for (Int_t i=0;i<fN;i++){
    AliTPCclusterFast * cluster = (AliTPCclusterFast*) fCl->UncheckedAt(i);
    cluster->Digitize();
  }

}

Double_t  AliTPCtrackFast::CookdEdxNtot(Double_t f0,Float_t f1){
  ///    Double_t  CookdEdxNtot(Double_t f0,Float_t f1);   //  dEdx_{hit}  reconstructed meen number of  electrons

  Double_t amp[160];
  for (Int_t i=0;i<fN;i++){ 
    AliTPCclusterFast * cluster = ( AliTPCclusterFast *)((*fCl)[i]);
    amp[i]=cluster->fNtot;
  }
  return CookdEdx(fN,amp,f0,f1,0);
}

Double_t  AliTPCtrackFast::CookdEdxQtot(Double_t f0,Float_t f1){
  ///     dEdx_{Q} reconstructed mean number of electronsxGain

  Double_t amp[160];
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

  Double_t amp[160];
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
  Double_t amp[160];
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

  Double_t amp[160];
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

  Double_t amp[160];
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
  Int_t index[160];
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

void AliTPCtrackFast::Simul(const char* fname, Int_t ntracks, Double_t diffFactor){
  ///

  AliTPCtrackFast fast;
  TTreeSRedirector *pcstream = new TTreeSRedirector(fname,"recreate");
  for (Int_t itr=0; itr<ntracks; itr++){
    //
    fast.fMNprim=(10.+100*gRandom->Rndm());
    if (gRandom->Rndm()>0.5) fast.fMNprim=1./(0.00001+gRandom->Rndm()*0.1);

    fast.fDiff =0.01 +0.35*gRandom->Rndm();
    //    fast.fDiffLong =  fast.fDiff*0.6/1.;
    fast.fDiffLong =  fast.fDiff*diffFactor/1.;
    //
    fast.fAngleY   = 4.0*(gRandom->Rndm()-0.5);
    fast.fAngleZ   = 4.0*(gRandom->Rndm()-0.5);
    fast.fN  = 160;
    fast.MakeTrack();
    if (itr%100==0) printf("%d\n",itr);
    (*pcstream)<<"simulTrack"<<
      "tr.="<<&fast<<
      "\n";
  }
  fast.Write("track");
  delete pcstream;
}



AliTPCclusterFast::AliTPCclusterFast(){
  ///

  fDigits.ResizeTo(5,7);
}

void AliTPCclusterFast::Init(){
  /// reset all counters

  const Int_t knMax=10000;
  fMNprim=0;     // mean number of primary electrons
  //                   //electrons part input
  fNprim=0;      // mean number of primary electrons
  fNtot=0;       // total number of  electrons
  fQtot=0;       // total charge - Gas gain flucuation taken into account
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
}



AliTPCclusterFast::~AliTPCclusterFast(){
}


void AliTPCclusterFast::SetParam(Float_t mnprim, Float_t diff,  Float_t diffL,Float_t y, Float_t z, Float_t ky, Float_t kz){
  ///

  fMNprim = mnprim; fDiff = diff; fDiffLong=diffL;
  fY=y; fZ=z; 
  fAngleY=ky; fAngleZ=kz;
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
  for (Int_t iprim=0; iprim<cl0->fNprim;iprim++){
    Float_t dN   =  cl0->GetNsec();
    cl0->fSec[iprim]=dN;
    Double_t yc = cl0->fY+(gRandom->Rndm()-0.5)*cl0->fAngleY;
    Double_t zc = cl0->fZ+(gRandom->Rndm()-0.5)*cl0->fAngleZ;
    Double_t rc = (gRandom->Rndm()-0.5);

    for (Int_t isec=0;isec<=dN;isec++){
      //
      //
      Double_t y = gRandom->Gaus(0,cl0->fDiff)+yc;
      Double_t z = gRandom->Gaus(0,cl0->fDiff)+zc;
      Double_t r = gRandom->Gaus(0,cl0->fDiffLong)+rc;
      // choose pad row
      AliTPCclusterFast *cl=cl0;
      if (r<-0.5 &&cl) cl=clm;
      if (r>0.5 &&cl)  cl=clp;
      //
      Double_t gg = -TMath::Log(gRandom->Rndm());
      cl->fPosY[cl->fNtot]=y;
      cl->fPosZ[cl->fNtot]=z;
      cl->fGain[cl->fNtot]=gg;
      cl->fQtot+=gg;
      cl->fNtot++;
      //
      //     cl->sumQ+=gg;
      //       cl->sumYQ+=gg*y;
      //       cl->sumY2Q+=gg*y*y;
      //       cl->sumZQ+=gg*z;
      //       cl->sumZ2Q+=gg*z*z;
      if (cl->fNtot>=knMax) continue;
    }
    if (cl0->fNtot>=knMax) break;
  }

 //  if (sumQ>0){
//     fStatY[0]=sumQ;
//     fStatY[1]=sumYQ/sumQ;
//     fStatY[2]=sumY2Q/sumQ-fStatY[1]*fStatY[1];
//     fStatZ[0]=sumQ;
//     fStatZ[1]=sumZQ/sumQ;
//     fStatZ[2]=sumZ2Q/sumQ-fStatZ[1]*fStatZ[1];
//   }
}

void AliTPCclusterFast::Digitize(){
  /// 1. Clear digits

  for (Int_t i=0; i<5;i++)
    for (Int_t j=0; j<7;j++){
      fDigits(i,j)=0;
    }
  //
  // Fill digits
  for (Int_t iel = 0; iel<fNtot; iel++){
    for (Int_t di=-2; di<=2;di++)
      for (Int_t dj=-3; dj<=3;dj++){
	Float_t fac = fPRF->Eval(di-fPosY[iel])*fTRF->Eval(dj-fPosZ[iel]);
	fac*=fGain[iel];
	fDigits(2+di,3+dj)+=fac;
      }
  }
  //
  //
  //
}



// void AliTPCclusterFast::Simul(const char* fname, Int_t npoints){
//   //
//   // Calc rms
//   //
//   AliTPCclusterFast fast;
//   TTreeSRedirector cstream(fname);
//   for (Int_t icl=0; icl<npoints; icl++){
//     Float_t nprim=(10+20*gRandom->Rndm());
//     Float_t diff =0.01 +0.35*gRandom->Rndm();
//     Float_t posY = gRandom->Rndm()-0.5;
//     Float_t posZ = gRandom->Rndm()-0.5;
//     //
//     Float_t ky   = 4.0*(gRandom->Rndm()-0.5);
//     Float_t kz   = 4.0*(gRandom->Rndm()-0.5);
//     fast.SetParam(nprim,diff,posY,posZ,ky,kz);
//     fast.GenerElectrons();
//     fast.Digitize();
//     if (icl%10000==0) printf("%d\n",icl);
//     cstream<<"simul"<<
//       "s.="<<&fast<<
//       "\n";
//   }
// }


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


Double_t  AliTPCclusterFast::Gamma4(Double_t x, Double_t p0, Double_t p1){
  /// Gamma 4 Time response function of ALTRO

  if (x<0) return 0;
  Double_t g1 = TMath::Exp(-4.*x/p1);
  Double_t g2 = TMath::Power(x/p1,4);
  return p0*g1*g2;
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
 
// Analytical sollution only in 1D - too long expression
// Simplify[Integrate[Exp[-(x0-(x1-k*x2))*(x0-(x1-k*x2))/(2*s0*s0)]*Exp[-(x1*t1-k*x2)],{x2,-1,1}]] 
//
//
// No analytical solution
// 
//Simplify[Integrate[Exp[-(x0-k0*xd)*(x0-k0*xd)/(2*s0*s0)-(x1-xt-k1*xd)*(x1-xt-k1*xd)/(2*s1*s1)]*Exp[-kt*xt]/(s0*s1),{xd,-1/2,1/2},{xt,0,Infinity}]]
