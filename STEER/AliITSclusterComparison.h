#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliTrackReference.h"
#endif


class AliITSCI: public TObject {
  public :
  AliITSCI(){;}
  virtual ~AliITSCI(){;}
  Float_t GetExpectedQ();
  //
  AliITSclusterV2 fCl;
  AliTrackReference fRef;
  TParticle fP;
  Char_t fStatus;
  //
  Float_t fBeta2;
  Float_t fVertex;
  Int_t fNClusters;
  Int_t   fLabPos;
  Float_t fExpectedQ;
  Float_t fNormQ;
  Int_t   fPDG;
  Int_t   fNTracks;
  Float_t fPt;
  Float_t fPtotal;
  Float_t fR;
  Float_t fCharge;
  Bool_t  fIsPrim;
  Double_t fDy;
  Double_t fDz;
  Double_t fPoolY;
  Double_t fPoolY2;
  Double_t fPoolZ;
  Double_t fPoolZ2;
  Double_t fErrY;
  Double_t fErrZ;
  Int_t    fErrType;
  Double_t fPhiCl;
  Double_t fRCl;
  Double_t fLx;
  Double_t fLy;
  Double_t fGx;
  Double_t fGy;
  Float_t  fTNy;
  Float_t  fTNz;
  Float_t  fTheta;
  Double_t fPhi;
  Double_t fPhiL;
  Int_t fLayer;
  void Update(AliITStrackerMI*tracker);
  ClassDef(AliITSCI,1)
};

class AliITSClusterErrAnal: public TObject{
public: 
  AliITSClusterErrAnal(Char_t *chloader  = "galice.root");
  void SetIO(Int_t event);  
  Int_t Analyze(Int_t trackmax);
  void LoadClusters();
  void LoadParticles();
  void SignDeltas( TObjArray *ClusterArray, Float_t vz);
  void SortReferences();
  AliITSclusterV2 * FindNearestCluster(AliITSCI * clinfo, Float_t dmax=0.5);
  void GetNTeor(Int_t layer, Float_t theta, Float_t phi, Float_t &ny, Float_t &nz);
  Int_t GetError(Int_t layer, const AliITSclusterV2*cl, Float_t theta, Float_t phi, Float_t expQ, Float_t &erry, Float_t &errz);
  void MakeSeeds(Double_t zv);
public:
  AliRunLoader * fRunLoader;
  AliLoader * fITSLoader;
  TTree * fHitTree;
  TTree * fClusterTree;
  TTree * fReferenceTree;
  AliITS * fITS;
  //
  TTree * fTreeA;
  TTree * fTreeB;
  TFile * fFileA;
  AliITStrackerMI *fTracker;
  AliStack *fStack;
  TObjArray *fClusters;  //array of clusters
  TObjArray fExactPoints;
  TObjArray *fReferences;
  TObjArray *fParticles;
  ClassDef(AliITSClusterErrAnal,1)
};


void AliITSCI::Update(AliITStrackerMI*tracker)
{
  //
  //
  fPt     = fRef.Pt();
  fR      = fRef.R();
  fCharge = fP.GetPDG()->Charge();
  //fIsPrim = (fP.GetMother(0)>0)? kFALSE :kTRUE;
  Double_t vertex = TMath::Sqrt( fP.Vx()*fP.Vx()
				+fP.Vy()*fP.Vy());
  fVertex = vertex; 
  fIsPrim = (vertex<1) ? kTRUE:kFALSE;
  fTheta  = fRef.Pz()/fRef.Pt();
  fPhi    =  TMath::ATan2(fRef.Py(),fRef.Px());  
  fPDG    = fP.GetPdgCode();
  fExpectedQ = GetExpectedQ();
}


Float_t AliITSCI::GetExpectedQ()
{
  //
  //
  //
  //fPtotal = fP.P();
  fPtotal = fRef.P();

  Double_t p2= fPtotal*fPtotal;
  Float_t mass = fP.GetPDG()->Mass();
  Double_t beta2=p2/(p2 + mass*mass);
  fBeta2 = beta2;
  Float_t d=0.003;
  //  Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2)) - beta2)*d;
  //dE+=0.33e-3/(beta2)*d;
  //if (beta2/(1-beta2)>3.5*3.5)
  //  dE=0.153e-3/beta2*(log(3.5*5940)+0.5*log(beta2/(1-beta2)) - beta2)*d;
  Double_t dE = 0.22*0.153e-3*(39.2-55.6*beta2+28.7*beta2*beta2+27.41/beta2)*d;
  
  return dE*10000000.; //normalization  to ADC
  //return dE; //normalization  to ADC
}

/*

.L AliITSclusterComparison.C+  
.L AliGenInfo.C+ 
 
AliITSClusterErrAnal anal;
anal.SetIO(0); 
anal.LoadParticles();
anal.SortReferences();
anal.LoadClusters();


anal.Analyze(30000000);


anal.MakeSeeds(-3.63); >seed1.txt


TFile f("itsclusteranal.root");
TTree * tree= (TTree*)f.Get("itscl")
TTree * tree2= (TTree*)f.Get("Clusters")

AliComparisonDraw comp;
comp.fTree = tree;

TH1F * hpools = new TH1F("pools","pools",100,-6,6); hpools->SetLineColor(4); hpools->SetFillColor(0);
TH1F * hpools0 = new TH1F("pools0","pools0",100,-6,6);  hpools0->SetLineColor(2); hpools0->SetFillColor(0);

TH1F * hsigma = new TH1F("hsigma","hsigma",1000,0,0.3)

TCut cprim("cprim","fCl.fTracks[0]<110000&&fCl.fTracks[1]<110000&&fCl.fTracks[2]<110000");

TCut cprim("cprim","fCl.fTracks[0]<24000&&fCl.fTracks[1]<24000&&fCl.fTracks[2]<24000");
TCut cpt("cpt","sqrt(fRef.fPx**2+fRef.fPy**2)>0.05");
TCut cphi("cphi","abs(fPhiL)<0.9");
TCut ctheta("ctheta","abs(fTheta)<1.1");
TCut cangle = cphi+ctheta;
TCut cpools5("cpools5","abs(fPoolY)<5&&abs(fPoolZ)<5");
TCut cpools52("cpools52","abs(fPoolY2)<5&&abs(fPoolZ2)<5");
TCut cpools62("cpools52","sqrt(fPoolY2**2+fPoolZ2**2)<6");
TCut cgold("cgold","abs(fCl.fNy+fCl.fNz-fTNy-fTNz)<1.5");
TCut cback("cback","fRef.fX*fRef.fPx+fRef.fPy*fRef.fY<0");

comp.fTree->Draw("fPoolY2>>pools","fStatus==1&&abs(fLayer-0.5)<0.6&&fIsPrim"+cpt+cangle+cpools62,"")

comp.DrawXY("fErrType","fPoolZ2","fStatus==1&&abs(fLayer-4.5)<0.6","1",8,99.9,108,-5,5)

comp.DrawXY("fCl.fType","fPoolY2","fStatus==1&&abs(fLayer-4.5)<0.6&&fIsPrim","1",4,102.9,106.1,-5,5)
comp.DrawXY("abs(fCl.fChargeRatio)","fPoolY2","fStatus==1&&abs(fLayer-4.5)<0.6&&fIsPrim","1",4,0.0,1,-5,5)

TProfile prof1("prof1","prof1",30,0.1,1); prof1.SetLineColor(2);
TProfile prof2("prof2","prof2",30,0.1,1); prof2.SetLineColor(3);

comp.fTree->Draw("fNormQ/fExpectedQ:fBeta2>>prof1","fStatus==1&&fCl.fLayer>1&&fNormQ/fExpectedQ<1.5&&fCl.fNz<4&&abs(fPDG)==211&&fIsPrim")
comp.fTree->Draw("fNormQ/fExpectedQ:fBeta2>>prof2","fStatus==1&&fCl.fLayer>1&&fNormQ/fExpectedQ<1.5&&fCl.fNz<4&&abs(fPDG)!=211&&fIsPrim")

prof1->Draw();prof2->Draw("same");


TH1F * hdeltam   = new TH1F("hdeltam" ,"hdeltam",100,-3,10);  hdeltam->SetLineColor(4); hdeltam->SetFillColor(0);
TH1F * hdelta0  =  new TH1F("hdelta0" ,"hdelta0",100,-3,10);  hdelta0->SetLineColor(2); hdelta0->SetFillColor(0);



*/
