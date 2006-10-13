#include "AliITSDetTypeRec.h"
#include "AliITSVertexerIons.h"
#include "AliITSVertexerPPZ.h"
#include "AliESDVertex.h"
#include "AliITSgeom.h"
#include "AliITSLoader.h"
#include "AliITSRecPoint.h"
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>

//////////////////////////////////////////////////////////////////////
// AliITSVertexerIons is a class for full 3D primary vertex         //
// finding optimized for Ion-Ion interactions                       //
//                                                                  // 
//                                                                  //
//                                                                  //
//                                                                  //
// Written by Giuseppe Lo Re and Francesco Riggi                    //
// Giuseppe.Lore@ct.infn.it                                         //
// Franco.Riggi@ct.infn.it                                          //
//                                                                  //
// Release date: Mar 2004                                           //
//                                                                  //
//                                                                  //       
//////////////////////////////////////////////////////////////////////


ClassImp(AliITSVertexerIons)
  
//______________________________________________________________________
  AliITSVertexerIons::AliITSVertexerIons():AliITSVertexer(),
fNpThreshold(0),
fMaxDeltaPhi(0),
fMaxDeltaZ(0){
  // Default Constructor

  //fITS = 0;
  SetNpThreshold();
  SetMaxDeltaPhi();
  SetMaxDeltaZ();
}

//______________________________________________________________________
AliITSVertexerIons::AliITSVertexerIons(TString fn):AliITSVertexer(fn),
fNpThreshold(0),
fMaxDeltaPhi(0),
fMaxDeltaZ(0) {
  // Standard constructor
  
  //fITS = 0;
  SetNpThreshold();
  SetMaxDeltaPhi();
  SetMaxDeltaZ();
}
/*
//______________________________________________________________________
AliITSVertexerIons::AliITSVertexerIons(const AliITSVertexerIons &source):AliITSVertexer(source) {
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  Error("AliITSVertexerIons","Copy constructor not allowed\n");
}
*/
//_________________________________________________________________________
//AliITSVertexerIons& AliITSVertexerIons::operator=(const AliITSVertexerIons &/*source*/) {
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  //Error("= operator","Assignment operator not allowed\n");
  //return *this;
//}


//______________________________________________________________________
AliITSVertexerIons::~AliITSVertexerIons() {
  // Default Destructor
  //fITS = 0;
}

//______________________________________________________________________
AliESDVertex* AliITSVertexerIons::FindVertexForCurrentEvent(Int_t evnumber){ 
// Defines the AliESDVertex for the current event

  fCurrentVertex = 0;

  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliITSLoader* itsloader = (AliITSLoader*)rl->GetLoader("ITSLoader");
  TDirectory * olddir = gDirectory;
  rl->CdGAFile();
  AliITSgeom* g2  = (AliITSgeom*)gDirectory->Get("AliITSgeom");
  olddir->cd(); 

  TTree *tr =  itsloader->TreeR();
  AliITSDetTypeRec detTypeRec;

  detTypeRec.SetTreeAddressR(tr);

  TClonesArray  *recpoints = detTypeRec.RecPoints();
  AliITSRecPoint *pnt;


  Int_t npoints=0;
  Int_t nopoints1=40000;
  Int_t nopoints2=40000;
  Float_t l[3], p[3];

  Double_t *z1, *z2, *y1, *y2, *x1, *x2, *phi1, *phi2, *r1, *r2;
  z1=new Double_t[nopoints1];
  z2=new Double_t[nopoints2];
  y1=new Double_t[nopoints1];
  y2=new Double_t[nopoints2];
  x1=new Double_t[nopoints1];
  x2=new Double_t[nopoints2];
  phi1=new Double_t[nopoints1];
  phi2=new Double_t[nopoints2];
  r1=new Double_t[nopoints1];
  r2=new Double_t[nopoints2];

  Double_t mxpiu = 0;
  Double_t mxmeno = 0;
  Double_t mypiu = 0;
  Double_t mymeno = 0;
  Double_t r=0;

  Int_t np1=0, np2=0;
  for(Int_t i=g2->GetStartSPD();i<=g2->GetLastSPD();i++) {
    detTypeRec.ResetRecPoints();
    tr->GetEvent(i);
    npoints = recpoints->GetEntries();
    for (Int_t ipoint=0;ipoint<npoints;ipoint++) {
      pnt = (AliITSRecPoint*)recpoints->UncheckedAt(ipoint);
      l[0]=pnt->GetDetLocalX();
      l[1]=0;
      l[2]=pnt->GetDetLocalZ();
      g2->LtoG(i, l, p);
      r=TMath::Sqrt(TMath::Power(p[0],2)+TMath::Power(p[1],2));

      if(i<80 && TMath::Abs(p[2])<14.35)  {
		y1[np1]=p[1];
		x1[np1]=p[0];
		z1[np1]=p[2];
		if(p[0]>0) mxpiu++;
		if(p[0]<0) mxmeno++;
		if(p[1]>0) mypiu++;
		if(p[1]<0) mymeno++;
		np1++;
      }
      if(i>=80 && TMath::Abs(p[2])<14.35) {
		y2[np2]=p[1];
		x2[np2]=p[0];
		z2[np2]=p[2];
		np2++;
      }
    }
  }

  if(np1<fNpThreshold) {
    Warning("FindVertexForCurrentEvent","AliITSVertexerIons finder is not reliable for low multiplicity events. Switching to AliITSVertexerPPZ with default parameters...\n");
    Warning("FindVertexForCurrentEvent","N rec points = %d - Threshold is %d",np1,fNpThreshold);
    AliITSVertexerPPZ *dovert = new AliITSVertexerPPZ("default");
    fCurrentVertex =dovert->FindVertexForCurrentEvent(rl->GetEventNumber());
    delete dovert;
    return fCurrentVertex;
  }

  if(!np1 || !np2) {
    Error("FindVertexForCurrentEvent","No points in the pixer layers");
    return fCurrentVertex;
  }  
  if(fDebug>0)   cout << "N. points layer 1 and 2 : " << np1 << " " << np2 << endl;

  Double_t asparx = (mxpiu-mxmeno)/(mxpiu+mxmeno);
  Double_t aspary = (mypiu-mymeno)/(mypiu+mymeno);

  Double_t x0 = 2.86*asparx;
  Double_t y0 = 2.86*aspary;

  if(fDebug>0)   cout << "Rough xy vertex = " << x0 << " " << y0 << endl;  

  for(Int_t i=0;i<np1;i++) {
    x1[i]-=x0;
    y1[i]-=y0;
    PhiFunc(x1[i],y1[i],phi1[i]);
    r1[i]=TMath::Sqrt(x1[i]*x1[i]+y1[i]*y1[i]);
  }
  for(Int_t i=0;i<np2;i++) {
    x2[i]-=x0;
    y2[i]-=y0;
    PhiFunc(x2[i],y2[i],phi2[i]);
    r2[i]=TMath::Sqrt(x2[i]*x2[i]+y2[i]*y2[i]);
  }
   
  Int_t nbinxy=400;
  Int_t nbinz=400;
  TH1F *hxv=new TH1F("hxv","",nbinxy,-5,5);
  TH1F *hyv=new TH1F("hyv","",nbinxy,-5,5);
  TH1F *hzv=new TH1F("hzv","",nbinz,-15,15);

  Double_t dphi;
  for(Int_t j=0; j<np1; j++) {
    for(Int_t k=0; k<np2; k++) {
      dphi=TMath::Abs(phi2[k]-phi1[j]);
      if(dphi>180) dphi = 360-dphi;
      if(dphi>fMaxDeltaPhi) continue; 
      hzv->Fill(z1[j]-(z2[k]-z1[j])/((r2[k]/r1[j])-1));
     }
  }

  // Fitting ...
  Double_t max;
  Int_t binMax;
  Double_t maxcenter;
  
  max = hzv->GetMaximum();
  binMax=hzv->GetMaximumBin();
  maxcenter=hzv->GetBinCenter(binMax);
  Double_t dxy=1.5;
  
  TF1 *fz = new TF1 ("fz","([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])))+[3]",maxcenter-dxy,maxcenter+dxy);  
  fz->SetParameter(0,max);
  fz->SetParameter(1,maxcenter);
  fz->SetParameter(2,0.1);
  fz->SetLineColor(kRed);
  hzv->Fit("fz","RQ0");

  if(TMath::Abs(x0)>0.4 || TMath::Abs(y0)>0.4) {
    Warning("FindVertexForCurrentEvent","AliITSVertexerIonsGeneral finder is not reliable for events with x and y vertex coordinates too far from the origin of the reference system. Their values will be set to 0.\n");
    Double_t position[3]={0,0,fz->GetParameter(1)};
    Double_t resolution[3]={0,0,0};
    Double_t snr[3]={0,0,0};
    Char_t name[30];
    if(fDebug>0)Info("FindVertexForCurrentEvent","Vertex found for event %d",evnumber);
    sprintf(name,"Vertex");
    fCurrentVertex = new AliESDVertex(position,resolution,snr,name);
    return fCurrentVertex;
  }
  
  for(Int_t j=0; j<np1; j++) {
    for(Int_t k=0; k<np2; k++) {
      if(TMath::Abs((z1[j]-(z2[k]-z1[j])/((r2[k]/r1[j])-1))-fz->GetParameter(1))>fMaxDeltaZ) continue;
      if(y2[k]==y1[j]) continue;
      hxv->Fill(x0+(x2[k]-((x2[k]-x1[j])/(y2[k]-y1[j]))*y2[k]));
      if(x2[k]==x1[j]) continue;
      hyv->Fill(y0+(y2[k]-((y2[k]-y1[j])/(x2[k]-x1[j]))*x2[k]));
    }
  }
    
  TF1 *fx = new TF1 ("fx","([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])))+[3]+[4]*x+[5]*x*x",x0-dxy,x0+dxy);  
  fx->SetParameter(0,100);
  Double_t dist=0.3;
  Double_t xapprox=FindMaxAround(x0,hxv,dist);
  if(fDebug>0) cout << "xapprox = " << xapprox << endl;
  fx->SetParameter(1,xapprox);
  Double_t difcentroid=0.07;
  fx->SetParLimits(1,xapprox-difcentroid,xapprox+difcentroid);
  fx->SetParameter(2,0.1);
  fx->SetLineColor(kRed);
  hxv->Fit("fx","RQW0"); 
  
  TF1 *fy = new TF1 ("fy","([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])))+[3]+[4]*x+[5]*x*x",y0-dxy,y0+dxy);  
  fy->SetParameter(0,100);
  Double_t yapprox=FindMaxAround(y0,hyv,dist);
  if(fDebug>0) cout << "yapprox = " << yapprox << endl;
  fy->SetParameter(1,yapprox);
  fy->SetParLimits(1,yapprox-difcentroid,yapprox+difcentroid);
  fy->SetParameter(2,0.1);
  fy->SetLineColor(kRed);
  hyv->Fit("fy","RQW0");
  
  delete [] z1;
  delete [] z2;
  delete [] y1;
  delete [] y2;
  delete [] x1;
  delete [] x2;
  delete [] r1;
  delete [] r2;
  delete [] phi1;
  delete [] phi2;
  delete hxv;
  delete hyv;
  delete hzv;
  
  Double_t position[3]={fx->GetParameter(1),fy->GetParameter(1),fz->GetParameter(1)};
  Double_t resolution[3]={fx->GetParameter(2),fy->GetParameter(2),fz->GetParameter(2)};
  Double_t snr[3]={0,0,0};
  
  Char_t name[30];
  if(fDebug>0)Info("FindVertexForCurrentEvent","Vertex found for event %d",evnumber);
  sprintf(name,"Vertex");
  fCurrentVertex = new AliESDVertex(position,resolution,snr,name);

  return fCurrentVertex;
}

//______________________________________________________________________
void AliITSVertexerIons::PhiFunc(Double_t &x,Double_t &y,Double_t &phi) {
  // Method for the computation of the phi angle given x and y. The result is in degrees. 
  if(y>0 && x>0) phi=(TMath::ATan((Double_t)(y/x))*57.29578);
  if(y>0 && x<0) phi=(TMath::ATan((Double_t)(y/x))*57.29578)+180;
  if(y<0 && x<0) phi=(TMath::ATan((Double_t)(y/x))*57.29578)+180;
  if(y<0 && x>0) phi=(TMath::ATan((Double_t)(y/x))*57.29578)+360;;
}

//______________________________________________________________________
void AliITSVertexerIons::FindVertices(){
  // computes the vertices of the events in the range FirstEvent - LastEvent
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliITSLoader* itsloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  itsloader->LoadRecPoints("read");
  for(Int_t i=fFirstEvent;i<=fLastEvent;i++){
    rl->GetEvent(i);
    FindVertexForCurrentEvent(i);
    if(fCurrentVertex){
      WriteCurrentVertex();
    }
    else {
      if(fDebug>0){
	cout<<"Vertex not found for event "<<i<<endl;
      }
    }
  }
}

//________________________________________________________
void AliITSVertexerIons::PrintStatus() const {
  // Print current status
  cout <<"=======================================================\n";
  cout <<" Debug flag: "<<fDebug<<endl;
  cout<<"First event to be processed "<<fFirstEvent;
  cout<<"\n Last event to be processed "<<fLastEvent<<endl;
  if(fCurrentVertex)fCurrentVertex->PrintStatus();
}

//________________________________________________________
Double_t AliITSVertexerIons::FindMaxAround(Double_t point, TH1F *h, Double_t distance) {
  // It finds the maximum of h within point-distance and distance+point. 
  Int_t maxcontent=0;
  Int_t maxbin=0;
  for(Int_t i=0;i<h->GetNbinsX();i++) {
    Int_t content=(Int_t)h->GetBinContent(i);
    Double_t center=(Double_t)h->GetBinCenter(i);
    if(TMath::Abs(center-point)>distance) continue;
    if(content>maxcontent) {maxcontent=content;maxbin=i;}
  }
  Double_t max=h->GetBinCenter(maxbin);
  return max;
}


