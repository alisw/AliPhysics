#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliTrackReference.h"
#endif

AliTRDclusterCorrection * gCorrection;

void ReadCorrection(){
  TFile f("TRDcorrection.root");
  gCorrection= (AliTRDclusterCorrection *)f.Get("TRDcorrection");
  if (gCorrection==0){
    printf("Correction not found");
  }
}


class AliTRDExactPoint: public TObject {
  public : 
  AliTRDExactPoint();
  Float_t fTX;    //x in rotated coordinate in the center of time bin
  Float_t fTY;    //y 
  Float_t fTZ;    //z
  Float_t fTAY;   //angle y
  Float_t fTAZ;   //angle z
  Float_t fGx;
  Float_t fGy;
  Float_t fGz;
  //
  void SetReference(AliTrackReference *ref);
  Float_t fTRefAngleY;
  Float_t fRefPos[3];
  Float_t fRefMom[3];
  //
  Int_t   fDetector;      // detector (chamber)
  Int_t   fLocalTimeBin;  // local time bin
  Int_t   fPlane;         // plane (layer)
  Int_t   fSector;       // segment
  Int_t   fPlaneMI;  
  // 
  Float_t fTQ;
  Float_t fTPrim;
  //  
  ClassDef(AliTRDExactPoint,1)
};

class AliTRDCI: public TObject {
  public :
  AliTRDCI(){;}
  virtual ~AliTRDCI(){;}
  //
  AliTRDclusterMI fCl;
  AliTRDExactPoint fEp;
  TParticle fP;
  Char_t fStatus;
  //
  Float_t fDYtilt;
  Float_t fTDistZ;  
  //
  Int_t   fNTracks;
  Float_t fPt;
  Float_t fCharge;
  Bool_t  fIsPrim;
  Float_t fCorrection;
  void Update();
  ClassDef(AliTRDCI,1)
};

class AliTRDClusterErrAnal: public TObject{
public: 
  AliTRDClusterErrAnal(Char_t *chloader  = "galice.root");
  void SetIO(Int_t event);  
  Int_t Analyze(Int_t trackmax);
  void LoadClusters();
  void MakeExactPoints(Int_t trackmax);
  void SortReferences();
  AliTrackReference * FindNearestReference(Int_t lab, Float_t pos[3], Float_t dmax=10.);
public:
  AliRunLoader * fRunLoader;
  AliLoader * fTRDLoader;
  AliTRDparameter *fParam;
  AliTRDgeometry *fGeometry;
  TTree * fHitTree;
  TTree * fClusterTree;
  TTree * fReferenceTree;
  AliTRDv1 * fTRD;
  //
  TTree * fTreeA;
  TFile * fFileA;
  AliTRDtracker *fTracker;
  AliStack *fStack;
  TObjArray fClusters[12][100][18];  //first plane, second time bin
  TObjArray fExactPoints;
  TObjArray *fReferences;

  ClassDef(AliTRDClusterErrAnal,1)
};


class AliTRDClusterErrDraw: public TObject{
public:
  AliTRDclusterCorrection*   MakeCorrection(TTree * tree, Float_t offset);

  static TH1F * ResDyVsAmp(TTree* tree, const char* selection, Float_t t0, Float_t ampmin=10, Float_t ampmax=300);
  static TH1F * ResDyVsRelPos(TTree* tree, const char* selection, Float_t t0, Float_t min=-0.5, Float_t max=0.5);
  static TH1F * ResDyVsAngleY(TTree* tree, const char* selection, Float_t t0, Float_t min=-1., Float_t max=1.);
  static void AliLabelAxes(TH1* histo, const char* xAxisTitle, const char* yAxisTitle);
  static TH1F* CreateEffHisto(TH1F* hGen, TH1F* hRec);
  static TH1F* CreateResHisto(TH2F* hRes2, Bool_t draw = kTRUE, Bool_t drawBinFits = kTRUE, 
		     Bool_t overflowBinFits = kFALSE);
  ClassDef(AliTRDClusterErrDraw,1)
};




AliTRDExactPoint::AliTRDExactPoint()
{
  fTX=fTY=fTZ=fTAZ=fTAY=fGx=fGy=fGz=fTRefAngleY=0;
  fRefPos[0]=fRefPos[1]=fRefPos[2]=fRefMom[0]=fRefMom[1]=fRefMom[2]=0;
  fDetector=fLocalTimeBin=fPlane=fSector=fPlaneMI=0;
  fTQ=fTPrim=0;
}

void AliTRDExactPoint::SetReference(AliTrackReference *ref){
  fRefPos[0] = ref->X();
  fRefPos[1] = ref->Y();
  fRefPos[2] = ref->Z();
  //
  fRefMom[0] = ref->Px();
  fRefMom[1] = ref->Py();
  fRefMom[2] = ref->Pz();
}


void AliTRDCI::Update()
{
  //
  //thanks to root
  fPt = fP.Pt();
  fCharge = fP.GetPDG()->Charge();
  fIsPrim = (fP.GetMother(0)>0)? kFALSE :kTRUE;
  fCorrection = gCorrection->GetCorrection(fEp.fPlane,fCl.fTimeBin,fEp.fTAY);
}

