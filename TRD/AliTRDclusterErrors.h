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
  Int_t  fNClusters;
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
  TTree * fTree;
  AliTRDclusterCorrection*   MakeCorrection(TTree * tree, Float_t offset);
  static TH1F * ResDyVsAmp(TTree* tree, const char* selection, Float_t t0, Float_t ampmin=10, Float_t ampmax=300);
  static TH1F * ResDyVsRelPos(TTree* tree, const char* selection, Float_t t0, Float_t min=-0.5, Float_t max=0.5);
  static TH1F * ResDyVsAngleY(TTree* tree, const char* selection, Float_t t0, Float_t min=-1., Float_t max=1.);
  static void   AliLabelAxes(TH1* histo, const char* xAxisTitle, const char* yAxisTitle);
  static TH1F*  CreateEffHisto(TH1F* hGen, TH1F* hRec);
  static TH1F*  CreateResHisto(TH2F* hRes2, Bool_t draw = kTRUE, Bool_t drawBinFits = kTRUE, 
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


/*
//example seesion

.L AliGenInfo.C+
.L AliTRDclusterErrors.C+
gCorrection  = AliTRDclusterCorrection::GetCorrection();
AliTRDClusterErrAnal ana;
ana.Analyze(10000000)



.L AliGenInfo.C+
.L AliTRDclusterErrors.C+
TFile f("trdclusteranal.root");
TTree* tree = (TTree*)f.Get("trdcl");
AliComparisonDraw comp;
comp->fTree = tree;


tree->SetAlias("shapef","(1.-(0.8+0.06*(6-fEp.fPlane))*(fCl.fSigmaY2/(0.17+0.027*abs(fEp.fTAY))))");
tree->SetAlias("shapes","0.08+0.3/sqrt(fCl.fQ)");
tree->SetAlias("sfactor","shapef/shapes");

tree->SetAlias("shapen","(fCl.fNPads-2.7-0.9*abs(fEp.fTAY))");

tree->SetAlias("gshape","sfactor>-2&&fCl.fNPads<6&&shapen<1");


tree->SetAlias("dy"    , "fEp.fTY-fCl.fY-fDYtilt");
tree->SetAlias("angle","abs(fEp.fTAY)");
TCut cbase("cbase","(abs(fP.fPdgCode)!=11||fPt>0.2)&&fPt>0.1&&angle<2");

tree->SetAlias("erry0","(0.028+0.07*angle)");
tree->SetAlias("erry1","erry0*(0.9+15./fCl.fQ)");
tree->SetAlias("erry2","erry1*(0.8+0.5*max(-sfactor,0))"); 




TH1F his("resy","resy",100,-0.2,0.2);
comp->fTree->Draw("dy:0.028*fEp.fTAY","fStatus==0 && abs(dy)<1.0&&fNTracks<2&&angle<2&&gshape","")
comp->DrawXY("sqrt(fCl.fQ)","dy","fStatus==0"+cbase,"gshape",10,0,20,-0.7,0.7);

comp->DrawXY("angle","dy/erry1","fStatus==0"+cbase,"gshape",10,0,2,-5.7,5.7);
comp->DrawXY("sqrt(cl->fQ)","dy/err1","fStatus==0"+cbase,"gshape",10,2,20,-5.7,5.7);



AliTRDClusterErrDraw::ResDyVsAmp(tree,"abs(dy)<0.4&&abs(fEp.fTAY)<0.2&&fNTracks<2&&fPt>0.1"+cbase,-0.03);
AliTRDClusterErrDraw::ResDyVsRelPos(tree,"abs(dy)<0.4&&abs(fEp.fTAY)<0.2&&fNTracks<2&&fPt>0.2"+cbase,-0.03);
AliTRDClusterErrDraw::ResDyVsAngleY(tree,"abs(fEp.fTY-fCl.fY+fDYtilt)<0.4",0.);

AliTRDclusterCorrection * cor = AliTRDClusterErrDraw::MakeCorrection(tree,0.134);
tree->Draw("sqrt(fCl.fRmsY)","fStatus==0&&abs(fEp.fTY-fCl.fY+fDYtilt)<0.4&&abs(fEp.fTAY)<0.2&&fNTracks<2&&fPt>0.2&&fCl.fNPads==0");

tree->Draw("fEp.fTY-fCl.fY+fDYtilt:fCl.fTimeBin","fStatus==0&&abs(fEp.fTY-fCl.fY+fDYtilt)<0.5&&abs(fEp.fTAY)<0.3&&fEp.fTAY>0.13","prof") 


*/

