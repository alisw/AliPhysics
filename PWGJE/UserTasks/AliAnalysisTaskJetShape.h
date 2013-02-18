#ifndef ALIANALYSISTASKJETSHAPE_h
#define ALIANALYSISTASKJETSHAPE_h

#include "AliAnalysisTaskSE.h"
#include "AliVEvent.h"
//#include "AliTriggerAnalysis.h"
//#include "THnSparse.h"
#include "TMath.h"
#include "TString.h"
#include "TH1F.h"
#include "TClonesArray.h"

class TH2F;
class TH3F;
class TVector3;

class AliAODExtension;
class TNtuple;
class TTree;
class TList;
class AliAODEvent;
class AliAODJet;
class TArrayD;







class AliAnalysisTaskJetShape : public AliAnalysisTaskSE {

 public:


class AliAnalysisTaskJetShapeTool :public TObject {
public:
AliAnalysisTaskJetShapeTool();
//AnalJetSubStrTool(TList *list, TVector3 v);
AliAnalysisTaskJetShapeTool(TClonesArray *list);
 ~AliAnalysisTaskJetShapeTool();


 void SetVecJ(TVector3 v);
 void SetPtMinTr(Double_t a, Double_t b) {fPtMinTr[0] = a; fPtMinTr[1] = b;}

 void SetNEntries(Int_t n) {fEntries = n;}

 TArrayI GetListJ(Int_t b, Int_t i)  {return fListJ[b][i];}
 TArrayI GetListB1(Int_t b, Int_t i) {return fListB1[b][i];}
 TArrayI GetListB2(Int_t b, Int_t i) {return fListB2[b][i];}

 Int_t GetSizeJ(Int_t b, Int_t i)  {return fListJc[b][i];}
 Int_t GetSizeB1(Int_t b, Int_t i) {return fListB1c[b][i];}
 Int_t GetSizeB2(Int_t b, Int_t i) {return fListB2c[b][i];}

 TVector3 *GetJ(Int_t b,Int_t i, Int_t p)
  { return (TVector3*)GetAt(fListJ[b][i].At(p));  }
 TVector3 *GetB1(Int_t b,Int_t i, Int_t p)
  { return (TVector3*)GetAt(fListB1[b][i].At(p)); }
 TVector3 *GetB2(Int_t b,Int_t i, Int_t p)
  { return (TVector3*)GetAt(fListB2[b][i].At(p)); }

 TVector3 GetVecJ() {return  fvecJ;}
 TVector3 GetVecB1() {return fvecB1;}
 TVector3 GetVecB2() {return fvecB2;}

 Double_t GetR(Int_t b) {return fR[b];}

 TVector3 *GetAt(Int_t i)
 { 
   if(i<0 || i>= fEntries) printf(" TVector3 *GetAt(Int_t i) for i= %d\n",i);
return (TVector3*) fList->At(i); 
}
 
 Double_t GetLocPhi(TVector3 v, Int_t i)
 {
   TVector3 p(*GetAt(i));
   //   GetLocalMom(v, &vtmp);
   //   return CalcdPhi0To2pi(vtmp.Phi());

   Double_t phi =  TMath::ATan2(p.Eta() - v.Eta(), CalcdPhi(p.Phi(), v.Phi()) );

   return CalcdPhi0To2pi(phi);
 }


 Double_t GetLocPhiJ(Int_t b, Int_t i, Int_t index)
 {
   TVector3 v1(*GetAt(fListJ[b][i].At(index)));
   Double_t phi =  TMath::ATan2(v1.Eta() - fvecJ.Eta(), CalcdPhi(v1.Phi(), fvecJ.Phi()) );
   return CalcdPhi0To2pi(phi);
 }


 Double_t GetLocPhiB1(Int_t b, Int_t i, Int_t index) 
 {
   TVector3 v0(*GetAt(fListB1[b][i].At(index)));
   Double_t cV = (fvecB1(0)*fvecJ(0) + fvecB1(1)*fvecJ(1))/fvecJ.Perp2();
   Double_t sV = (fvecB1(1)*fvecJ(0) - fvecB1(0)*fvecJ(1))/fvecJ.Perp2();
   TVector3 v1(v0(0)*cV+v0(1)*sV, v0(1)*cV-v0(0)*sV, v0(2));
   //   GetLocalMom(fvecJ, &v1);
   //   return CalcdPhi0To2pi(v1.Phi());

   Double_t phi =  TMath::ATan2(v1.Eta() - fvecJ.Eta(), CalcdPhi(v1.Phi(), fvecJ.Phi()) );
   return CalcdPhi0To2pi(phi);

 }
 Double_t GetLocPhiB2(Int_t b, Int_t i, Int_t index) 
 {
   TVector3 v0(*GetAt(fListB2[b][i].At(index)));
   Double_t cV = (fvecB2(0)*fvecJ(0) + fvecB2(1)*fvecJ(1))/fvecJ.Perp2();
   Double_t sV = (fvecB2(1)*fvecJ(0) - fvecB2(0)*fvecJ(1))/fvecJ.Perp2();
   TVector3 v1(v0(0)*cV+v0(1)*sV, v0(1)*cV-v0(0)*sV, v0(2));
   //   GetLocalMom(fvecJ, &v1);
   //   return CalcdPhi0To2pi(v1.Phi());
   Double_t phi =  TMath::ATan2(v1.Eta() - fvecJ.Eta(), CalcdPhi(v1.Phi(), fvecJ.Phi()) );
   return CalcdPhi0To2pi(phi);

 }


 void Add(TVector3 v) 
   // { new(fList[fList.GetEntriesFast()]) TVector3(v);}
 { new((*fList)[fEntries]) TVector3(v); fEntries++;}

 void AddToJ(Int_t b, Int_t i, Int_t index) {
   fListJ[b][i].AddAt(index, fListJc[b][i]);
   fListJc[b][i]++;
 }

 void AddToB1(Int_t b, Int_t i, Int_t index) {
   fListB1[b][i].AddAt(index, fListB1c[b][i]);
   fListB1c[b][i]++;
 }

 void AddToB2(Int_t b, Int_t i, Int_t index) {
   fListB2[b][i].AddAt(index, fListB2c[b][i]);
   fListB2c[b][i]++;
 }


 

 TVector3 GetPRecJ(Int_t b, Int_t i) {return fPRecJ[b][i];}
 TVector3 GetPRecInRJ() {return fPRecInRJ;}//


 void Init();
 // Bool_t FindCorrelationAxesAnd(TArrayI list, TVector3 vec, Double_t &Phi, Int_t scenario=0);
 Bool_t FindCorrelationAxes(TArrayI list, TVector3 vec, Double_t &Phi, Int_t scenario=3);
 Bool_t FindCorrelationAxesCorr(TArrayI list, TVector3 vec, Double_t &Phi, Int_t scenario, Double_t R);

 // virtual  void  Print(Option_t* /*option = ""*/) const;
 void  PRINT() const;
 void PRINT(TArrayI a, Int_t n, const char *txt="") const;

 void Clean();

static Double_t CalcdPhi0To2pi(Double_t phi1, Double_t phi2=0.)
{
  Double_t dphi = CalcdPhi(phi1, phi2);
  while(dphi<0) dphi+=TMath::TwoPi();
  while(dphi>TMath::TwoPi()) dphi-=TMath::TwoPi();
  return dphi;
}

static Double_t CalcdPhi(Double_t phi1, Double_t phi2);

static Double_t CalcdPhiSigned(Double_t phi1, Double_t phi2)
{

  Double_t dphi = CalcdPhi(phi1, phi2);
  if(dphi <  TMath::Pi()) return dphi;
  else return  (dphi-  TMath::Pi());

  return -999;
}


private:

 TVector3 fvecJ;
 TVector3 fvecB1;
 TVector3 fvecB2;


 Double_t fRmax;

 static const Int_t fgkbinR= 3; // n bins
 Double_t fR[fgkbinR]; 
 TArrayI fListJ[fgkbinR][2]; //
 TArrayI fListB1[fgkbinR][2];//
 TArrayI fListB2[fgkbinR][2];//

 Int_t fListJc[fgkbinR][2]; //
 Int_t fListB1c[fgkbinR][2];//
 Int_t fListB2c[fgkbinR][2];//

 TVector3 fPRecJ[fgkbinR][2];//
 TVector3 fPRecInRJ;//

 Double_t fPtMinTr[2];


 TClonesArray *fList; //!

 Int_t fEntries;

 Double_t CalcR(TVector3 v1, TVector3 v2);
 void GetLocalMom(TVector3 vecJ, TVector3 *q);


  AliAnalysisTaskJetShapeTool(const AliAnalysisTaskJetShapeTool&);            // not implemented
  AliAnalysisTaskJetShapeTool& operator=(const AliAnalysisTaskJetShapeTool&); // not implemented

 ClassDef(AliAnalysisTaskJetShapeTool, 1)   // tbd

};




class AliAnalysisTaskJetShapeHM :public TObject {

public:
  AliAnalysisTaskJetShapeHM(TString comment = "", Bool_t kMC = kFALSE);
 ~AliAnalysisTaskJetShapeHM();


 void SetPtJetBins(Int_t Nbin=-1, Double_t *array = 0);
 void SetRBins(    Int_t Nbin = 8, Double_t Rmax = 0.4) {fPsiVsRNbin = Nbin; fRmax = Rmax;}
 void SetPhiNbins(Int_t Nbin = 1024) {fPsiVsPhiNbin = Nbin;}
 void SetFilterMask(UInt_t i = 0){fFilterMask = i;}
 void SetEtaTrackMax(Double_t e =0.9) {fEtaTrackMax = e;}
 void SetPtTrackRange(Double_t min = 1., Double_t max = 100.) {fPtTrackMin = min; fPtTrackMax = max;}
 void SetPtJetRange(Double_t min = 1., Double_t max = 200. ) {fPtJmin  =min; fPtJmax = max; }

 void MCprod(Bool_t kMCprod=kTRUE) {fkMCprod=kMCprod; }


 Double_t GetPtJ() {return fPtJ;}
 void InitHistos();
 Bool_t AddEvent(AliAODEvent* aodE,  AliAODJet *jet, Double_t DeltaPtJ=0.);
 void AddToList(TList *list);


private:
 TString fComment;
 Bool_t fkMC;
 Bool_t fkMCprod;
 TH1F *fhEvent;//! 
 TH1F *fhMult; //! 
 TH1F *fhPtJ;//! 
 TH2F *fhPtJvsPtCorr;//! 
 TH3F *fhPsiVsR;//! 
 TH2F *fhPsiVsRPtJ;//! 
 TH2F *fhPhiEtaTrack;//! 

 TH1F *fhPsiVsR_MCtr; //!
 TH2F *fhPsiVsRPtJ_MCtr; //!
 TH2F *fhJetTrPtVsPartPt;//!


 TH1F *fhTMA_JAA[3];//! 
 TH1F *fhTMA_JAp[3];//! 
 TH1F *fhTMA_B1AA[3];//! 
 TH1F *fhTMA_B1Ap[3];//! 
 TH1F *fhTMA_B2AA[3];//! 
 TH1F *fhTMA_B2Ap[3];//! 

    TH2F *fhPtresVsPt[3][2] ;//!
    TH2F *fhPhiresVsPhi[3][2];//!
    TH2F *fhEtaresVsEta[3][2];//!
    TH2F *fhRresVsPt[3][2];//!
    TH2F *fhDCAxy[3][2]; //!
    TH2F *fhDCAz[3][2]; //!
    TH3F *fhTrackPtEtaPhi[3][2];//!

 Int_t     fPtJetNbin;
 TArrayD   fPtJetArray;
 Int_t     fPsiVsRNbin;
 Double_t  fRmax;
 Int_t     fPsiVsPhiNbin;

 UInt_t   fFilterMask;
 Double_t fEtaTrackMax;
 Double_t fPtTrackMin;
 Double_t fPtTrackMax;
 Double_t fPtJmin;
 Double_t fPtJmax;
 Double_t fPtJ;

 TVector3 fJvec;


 Double_t CalcR(TVector3 v1, TVector3 v2);
 Double_t CalcdPhi(Double_t phi1, Double_t phi2);

 TH1F* Hist1D(const char* name, Int_t nBins, Double_t xMin, Double_t xMax,  const char* xLabel = "", Int_t color=1, const char* yLabel="");
 TH1F* Hist1D(const char* name, Int_t nBins, Double_t *xArray,  const char* xLabel = "", Int_t color=1, const char* yLabel="");
 TH2F *Hist2D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, Int_t nBinsy, Double_t yMin, Double_t yMax, const char* xLabel = NULL, const char* yLabel = NULL, Int_t color=1);
 TH2F *Hist2D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, Int_t nBinsy, Double_t *yArray, const char* xLabel = NULL, const char* yLabel = NULL, Int_t color=1, const char* zLabel = NULL);
 TH2F *Hist2D(const char* name, Int_t nBinsx, Double_t *yArrax, Int_t nBinsy, Double_t yMin, Double_t yMax, const char* xLabel = NULL, const char* yLabel = NULL, Int_t color=1, const char* zLabel = NULL);
  TH3F *Hist3D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, Int_t nBinsy, Double_t yMin, Double_t yMax, Int_t nBinsz, Double_t zMin, Double_t zMax, const char* xLabel = NULL, const char* yLabel = NULL, const char* zLabel = NULL, Int_t color=1);
  TH3F *Hist3D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, Int_t nBinsy, Double_t yMin, Double_t yMax, Int_t nBinsz, Double_t *z, const char* xLabel = NULL, const char* yLabel = NULL, const char* zLabel = NULL, Int_t color=1);


  AliAnalysisTaskJetShapeHM(const AliAnalysisTaskJetShapeHM&);            // not implemented
  AliAnalysisTaskJetShapeHM& operator=(const AliAnalysisTaskJetShapeHM&); // not implemented


ClassDef(AliAnalysisTaskJetShapeHM,1)   // tbd
};




  AliAnalysisTaskJetShape(const char *name = "");
  ~AliAnalysisTaskJetShape();

  virtual void SetIsMC(Bool_t ismc=kTRUE) {fkMC=ismc;} 
  virtual void SetCMSE(Double_t s = 7000.) {fCMSE=s;}


  virtual void     SetNonStdFile(TString c){fNonStdFile = c;} 
  virtual void     SetBranchNames(const TString &branch1, const TString &branch2);
  virtual void     SetBackgroundBranch(TString &branch1, TString &branch2) { fBackgroundBranch[0] = branch1;  fBackgroundBranch[1] = branch2;};

  virtual void     SetFilterMask(UInt_t i){fFilterMask = i;}
  virtual void     SetOfflineTrgMask(AliVEvent::EOfflineTriggerTypes mask) { fOfflineTrgMask = mask; }
  virtual void     SetCentMin(Float_t cent) { fCentMin = cent; }
  virtual void     SetCentMax(Float_t cent) { fCentMax = cent; }
  virtual void     SetEvtClassMin(Int_t evtClass) { fEvtClassMin = evtClass; }
  virtual void     SetEvtClassMax(Int_t evtClass) { fEvtClassMax = evtClass; }
  virtual void     SetJetPtCorrMin(Float_t ptJ=20, Float_t ptB=20) { fPtJcorrMin = ptJ; fPtJBcorrMin = ptB;}

  virtual void     SetPbPb(Bool_t a = kTRUE) {fkIsPbPb = a;}
  // vtx
  virtual void     SetVtxZRange(Double_t zmin=-10., Double_t zmax=10.) {fVtxZMin = zmin; fVtxZMax = zmax;}
  virtual void     SetVixMinContrib(Int_t n=1) { fVtxMinContrib = n; }

   //  virtual void     SetAngStructCloseTracks(Int_t yesno){fAngStructCloseTracks=yesno;}


  //  virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);




 private:
  //  TTree *fOutputTree; //!Output tree

  TList       *fOutputList; // Output list
  AliESDEvent *fESD;    // ESD object
  AliAODEvent *fAODIn;    // AOD event
  AliAODEvent *fAODOut;    // AOD event
  AliAODExtension  *fAODExtension; //! where we take the jets from can be input or output AOD
 

  //  AliTriggerAnalysis * fTriggerAnalysis; // trigger analysis object, to get the offline triggers


   Bool_t   fkMC;         //is MC
   Double_t fCMSE;      //cms energy
   UInt_t fRunNb;       //run number
   Bool_t fkIsPhysSel;  //tbd
   TString       fNonStdFile; // delta AOD file
   UInt_t  fFilterMask;            // filter bit for slecected tracks
   AliVEvent::EOfflineTriggerTypes fOfflineTrgMask; // mask of offline triggers to accept
   Float_t fCentMin;	  // lower bound on centrality
   Float_t fCentMax;	  // upper bound on centrality
   Int_t   fEvtClassMin;	  // lower bound on event class
   Int_t   fEvtClassMax;	  // upper bound on event class
   Double_t fPtJcorrMin ;
   Double_t fPtJBcorrMin;
   Double_t fJpPtmin;
   Double_t fJaPtmin;
   Int_t fVtxMinContrib;
   Double_t fVtxZMin;
   Double_t fVtxZMax;
   Bool_t fkIsPbPb;


   TString fBackgroundBranch[2];  //tbd
   TString fJetBranchName[2]; //  name of jet branches to compare


   static const Int_t fgkbinNCent=7; // tbd
  Double_t fgkbinCent[fgkbinNCent+1] ;// {0, 5, 10,  20, 40, 60, 80, 100}; 

  TH1F *fhPtJL ;//!
  TH1F *fhAreaJL ;//!

  //  TClonesArray farray;
    AliAnalysisTaskJetShapeHM *fanalJetSubStrHM;//!
    AliAnalysisTaskJetShapeHM *fanalJetSubStrMCHM;//!

    TH2F *fhPtresVsPt[3] ;//!
    TH2F *fhPhiresVsPhi[3];//!
    TH2F *fhEtaresVsEta[3];//!
    TH2F *fhDCAxy[3]; //!
    TH2F *fhDCAz[3]; //!
    TH3F *fhTrackPtEtaPhi[3];//!

   Bool_t IsGoodEvent();


  Double_t CalcdPhi(Double_t phi1, Double_t phi2);
  Double_t CalcdPhi0To2pi(Double_t phi1, Double_t phi2);

  TH1F *Hist1D(const char* name, Int_t nBins , Double_t xMin, Double_t xMax, const char* xLabel = NULL, Int_t color=1);
  TH2F *Hist2D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, Int_t nBinsy, Double_t yMin, Double_t yMax, const char* xLabel = NULL, const char* yLabel = NULL, Int_t color=1);
  TH3F *Hist3D(const char* name, Int_t nBinsx, Double_t xMin, Double_t xMax, Int_t nBinsy, Double_t yMin, Double_t yMax, Int_t nBinsz, Double_t zMin, Double_t zMax, const char* xLabel = NULL, const char* yLabel = NULL, const char* zLabel = NULL, Int_t color=1);


  //  AliAnalysisTaskJetShape(const AnalysisJetMain&); // not implemented
  //  AliAnalysisTaskJetShape& operator=(const AnalysisJetMain&); // not implemented



  AliAnalysisTaskJetShape(const AliAnalysisTaskJetShape&);            // not implemented
  AliAnalysisTaskJetShape& operator=(const AliAnalysisTaskJetShape&); // not implemented


  ClassDef(AliAnalysisTaskJetShape, 1)
};









#endif
