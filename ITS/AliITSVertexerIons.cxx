#include <AliITSVertexerIons.h>
#include "stdlib.h"
#include <TMath.h>
#include <TRandom.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <Riostream.h>

#include "AliRun.h"
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSLoader.h"
#include "AliITSRecPoint.h"
#include "AliGenerator.h"
#include "AliMagF.h"

#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <AliITSVertex.h>
#include <TObjArray.h>
#include <TObject.h>
#include <AliITSVertexerPPZ.h>
//////////////////////////////////////////////////////////////////////
// AliITSVertexerIons is a class for full 3D primary vertex         //
// finding optimized for Ion-Ion interactions                       //
//                                                                  // 
//                                                                  //
//                                                                  //
//                                                                  //
// Written by Giuseppe Lo Re and Francesco Riggi                    //
//
// Giuseppe.Lore@ct.infn.it                                         //
// Franco.Riggi@ct.infn.it                                          //
//                                                                  //
// Release date: May 2001                                           //
//                                                                  //
//                                                                  //       
//////////////////////////////////////////////////////////////////////

ClassImp(AliITSVertexerIons)



  //______________________________________________________________________
  AliITSVertexerIons::AliITSVertexerIons():AliITSVertexer() {
  // Default Constructor

  fITS = 0;
  SetNpThreshold();
}

//______________________________________________________________________
AliITSVertexerIons::AliITSVertexerIons(TString fn):AliITSVertexer(fn) {
  // Standard constructor
  
  fITS = 0;
  SetNpThreshold();
}


//______________________________________________________________________
AliITSVertexerIons::~AliITSVertexerIons() {
  // Default Destructor
  fITS = 0;
}

//______________________________________________________________________
AliITSVertex* AliITSVertexerIons::FindVertexForCurrentEvent(Int_t evnumber){
  // Defines the AliITSVertex for the current event
  fCurrentVertex = 0;
  Double_t position[3];
  Double_t resolution[3];
  Double_t snr[3];
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliITSLoader* itsloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  if(!fITS) { 
    fITS =(AliITS *)gAlice->GetDetector("ITS");
    if(!fITS) {
      Error("FindVertexForCurrentEvent","AliITS object was not found");
      return fCurrentVertex;
    }
  }
  fITS->SetTreeAddress();
  AliITSgeom *g2 = fITS->GetITSgeom(); 
  TClonesArray  *recpoints = fITS->RecPoints();
  AliITSRecPoint *pnt;

  Int_t nopoints1 = 0;                                         
  Int_t nopoints2 = 0;
  Double_t vzero[3];
  Double_t aspar[2];
     
  //------------ Rough Vertex evaluation ---------------------------------   

  Int_t i,npoints,ipoint,j,k,max,binmax;
  Double_t zCentroid;
  Float_t l[3], p[3];

  Double_t mxpiu = 0;
  Double_t mxmeno = 0;
  Double_t mypiu = 0;
  Double_t mymeno = 0;
   
  TH1F *hITSz1         = new TH1F("hITSz1","",100,-14.35,14.35);
  TTree *tr =  itsloader->TreeR(); 
  for(i=g2->GetStartSPD();i<=g2->GetLastSPD();i++) 
    {
      fITS->ResetRecPoints();
      tr->GetEvent(i);    
      npoints = recpoints->GetEntries();
      for (ipoint=0;ipoint<npoints;ipoint++) {
	pnt = (AliITSRecPoint*)recpoints->UncheckedAt(ipoint);
	l[0]=pnt->GetX();
	l[1]=0;
	l[2]=pnt->GetZ();
	g2->LtoG(i, l, p);
	if(i<80 && TMath::Abs(p[2])<14.35) {         
	  if(p[0]>0) mxpiu++; 
	  if(p[0]<0) mxmeno++; 
	  if(p[1]>0) mypiu++; 
	  if(p[1]<0) mymeno++; 
	  hITSz1->Fill(p[2]);       
	}
	if(i>=80 && TMath::Abs(p[2])<14.35) nopoints2++;
      }       
    }  
  
  nopoints1 = (Int_t)(hITSz1->GetEntries()); 
           
  aspar[0] = (mxpiu-mxmeno)/(mxpiu+mxmeno);
  aspar[1] = (mypiu-mymeno)/(mypiu+mymeno); 
   
  vzero[0] = 5.24441*aspar[0]; 
  vzero[1] = 5.24441*aspar[1]; 

  zCentroid=  TMath::Abs(hITSz1->GetMean()); 
  vzero[2] = -5.31040e-02+1.42198*zCentroid+7.44718e-01*TMath::Power(zCentroid,2)
    -5.73426e-01*TMath::Power(zCentroid,3)+2.01500e-01*TMath::Power(zCentroid,4)          
    -3.34118e-02*TMath::Power(zCentroid,5)+2.20816e-03*TMath::Power(zCentroid,6);
   
  if(hITSz1->GetMean()<0) vzero[2] = -vzero[2];
   
  /*cout << "\nXvzero: " << vzero[0] << " cm" << "";
    cout << "\nYvzero: " << vzero[1] << " cm" << "";
    cout << "\nZvzero: " << vzero[2] << " cm" << "\n";*/

  delete hITSz1;

  Double_t dphi,r,deltaZ,a,b;
  Int_t np1=0;
  Int_t np2=0;
  Int_t niter=0;
   
  Double_t deltaPhiZ = 0.08;
   
  Float_t mag[3];
  Float_t origin[3];
  for(Int_t okm=0;okm<3;okm++) origin[okm]=0;
  gAlice->Field()->Field(origin,mag);

  deltaPhiZ = deltaPhiZ*mag[2]/2;
  Double_t deltaPhiXY = 1.0;   

  //   cout << "\ndeltaPhiZ: " << deltaPhiZ << " deg" << "\n";   
  //   cout << "deltaPhiXY: " << deltaPhiXY << " deg" << "\n";   
   
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
        
  deltaZ = 4.91617e-01+2.67567e-02*vzero[2]+1.49626e-02*TMath::Power(vzero[2],2); 
  Float_t multFactorZ = 28000./(Float_t)nopoints1;
  Int_t nbin=(Int_t)((deltaZ/0.005)/multFactorZ);
  Int_t nbinxy=250;
  Int_t *vectorBinZ,*vectorBinXY;
  vectorBinZ=new Int_t[nbin];
  vectorBinXY=new Int_t[nbinxy];
  Float_t f1= 0;
  Float_t f2= 0;
  Double_t sigma,averagebg;
     
  TH1D *hITSZv         = new TH1D("hITSZv","",nbin,vzero[2]-deltaZ,vzero[2]+deltaZ);
  TH1D *hITSXv         = new TH1D("hITSXv","",nbinxy,-3,3);
  TH1D *hITSYv         = new TH1D("hITSYv","",nbinxy,-3,3);

  //   cout << "deltaZeta: " << deltaZ << " cm" << "\n";   
   
   
 start:   

  hITSZv->Add(hITSZv,-1);
  hITSXv->Add(hITSXv,-1);
  hITSYv->Add(hITSYv,-1);

  np1=np2=0;


  for(i=g2->GetStartSPD();i<=g2->GetLastSPD();i++) 
    {
      fITS->ResetRecPoints(); 
      itsloader->TreeR()->GetEvent(i);
      npoints = recpoints->GetEntries();
      for (ipoint=0;ipoint<npoints;ipoint++) {
               
	pnt = (AliITSRecPoint*)recpoints->UncheckedAt(ipoint);
	l[0]=pnt->GetX();
	l[1]=0;
	l[2]=pnt->GetZ();
	g2->LtoG(i, l, p);
              
	if(i<80 && TMath::Abs(p[2])<14.35) {
	  p[0]=p[0]-vzero[0];
	  p[1]=p[1]-vzero[1];
	  r=TMath::Sqrt(TMath::Power(p[0],2)+TMath::Power(p[1],2));
	  y1[np1]=p[1];
	  x1[np1]=p[0];
	  z1[np1]=p[2]; 
	  r1[np1]=r;
	  phi1[np1]=PhiFunc(p);
	  np1++;
	}

	if(i>=80 &&  TMath::Abs(p[2])<14.35) { 
	  p[0]=p[0]-vzero[0];
	  p[1]=p[1]-vzero[1];
	  r=TMath::Sqrt(TMath::Power(p[0],2)+TMath::Power(p[1],2));
	  y2[np2]=p[1];
	  x2[np2]=p[0];
	  z2[np2]=p[2];
	  r2[np2]=r;
	  phi2[np2]=PhiFunc(p);
	  np2++;
	}
                        
      }
    }

  //------------------ Correlation between rec points ----------------------
    
  
  if(np1<fNpThreshold) {
    /*    cout << "Warning: AliITSVertexerIons finder is not reliable for low multiplicity events. May be you have to use AliITSVertexerPPZ.\n";
	  position[0]=vzero[0];
	  position[1]=vzero[1];
	  position[2]=vzero[2];
	  resolution[0]=-123;
	  resolution[1]=-123;
	  resolution[2]=-123;
	  snr[0]=-123;
	  snr[1]=-123;
	  snr[2]=-123;
	  Char_t name[30];
	  sprintf(name,"Vertex");
	  fCurrentVertex = new AliITSVertex(position,resolution,snr,name);
	  return fCurrentVertex;*/

    Warning("FindVertexForCurrentEvent","AliITSVertexerIons finder is not reliable for low multiplicity events. Switching to AliITSVertexerPPZ with default parameters...\n");
    Warning("FindVertexForCurrentEvent","N rec points = %d - Threshold is %d",np1,fNpThreshold);
    AliITSVertexerPPZ *dovert = new AliITSVertexerPPZ("default");
    fCurrentVertex =dovert->FindVertexForCurrentEvent(rl->GetEventNumber());
    delete dovert;
    return fCurrentVertex;
  }


  for(j=0; j<(np2)-1; j++) {
    for(k=0; k<(np1)-1; k++) { 
      dphi=TMath::Abs(phi1[k]-phi2[j]);
      if(dphi>180) dphi = 360-dphi;
      if(dphi<deltaPhiZ && TMath::Abs((z1[k]-(z2[j]-z1[k])/((r2[j]/r1[k])-1))-vzero[2])
	 <deltaZ) hITSZv->Fill(z1[k]-(z2[j]-z1[k])/((r2[j]/r1[k])-1));        
    }
  }
      
  //cout << "\nNumber of used pairs: \n";

  a = vzero[2]-deltaZ; 
  b = vzero[2]+deltaZ; 
  max=(Int_t) hITSZv->GetMaximum();
  binmax=hITSZv->GetMaximumBin();
  sigma=0;
  for(i=0;i<nbin;i++) vectorBinZ[i]=(Int_t)hITSZv->GetBinContent(i);
  for(i=0;i<10;i++) f1=f1+vectorBinZ[i]/10;
  for(i=nbin-10;i<nbin;i++) f2=f2+vectorBinZ[i]/10;
  averagebg=(f1+f2)/2;
  for(i=0;i<nbin;i++) {
    if(vectorBinZ[i]-averagebg>(max-averagebg)*0.4 && 
       vectorBinZ[i]-averagebg<(max-averagebg)*0.7) {
      sigma=hITSZv->GetBinCenter(binmax)-hITSZv->GetBinCenter(i);
      sigma=TMath::Abs(sigma);
      if(sigma==0) sigma=0.05;
    }
  }
  
  
  TF1 *fz = new TF1 ("fz","([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])))+[3]",a,b);  
  
  fz->SetParameter(0,max);
  if(niter==0) {Double_t temp = hITSZv->GetBinCenter(binmax); vzero[2]=temp;}
  fz->SetParameter(1,vzero[2]);
  fz->SetParameter(2,sigma); 
  fz->SetParameter(3,averagebg);  
  fz->SetParLimits(2,0,999);
  fz->SetParLimits(3,0,999);
  
  hITSZv->Fit("fz","RQ0");
  
  snr[2] = fz->GetParameter(0)/fz->GetParameter(3);
  if(snr[2]<0.) { 
    Error("FindVertexForCurrentEvent","\nNegative Signal to noise ratio for z!!!\n");
    Error("FindVertexForCurrentEvent","The algorithm cannot find the z vertex position.\n");
    //    exit(123456789);
    return fCurrentVertex;
  }
  else {
    position[2] = fz->GetParameter(1);
    if(position[2]<0) 
      {
	position[2]=position[2]-TMath::Abs(position[2])*1.11/10000;
      }
    else
      {
	position[2]=position[2]+TMath::Abs(position[2])*1.11/10000;
      }
  }
  resolution[2] = fz->GetParError(1);
  
  for(j=0; j<(np2)-1; j++) {
    for(k=0; k<(np1)-1; k++) { 
      dphi=TMath::Abs(phi1[k]-phi2[j]);
      if(dphi>180) dphi = 360-dphi;
      if(dphi>deltaPhiXY) continue;
      if(TMath::Abs((z1[k]-(z2[j]-z1[k])/((r2[j]/r1[k])-1))-position[2])
	 <4*resolution[2]) {
	hITSXv->Fill(vzero[0]+(x2[j]-((x2[j]-x1[k])/(y2[j]-y1[k]))*y2[j]));
	hITSYv->Fill(vzero[1]+(y2[j]-((y2[j]-y1[k])/(x2[j]-x1[k]))*x2[j]));
      }
    }
  }
  
  TF1 *fx = new TF1("fx","([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])))+[3]",vzero[0]-0.5,vzero[0]+0.5);  
  
  max=(Int_t) hITSXv->GetMaximum();
  binmax=hITSXv->GetMaximumBin();
  sigma=0;
  f1=f2=0;
  for(i=0;i<nbinxy;i++) vectorBinXY[i]=(Int_t)hITSXv->GetBinContent(i);
  for(i=0;i<10;i++) f1=f1+vectorBinXY[i]/10;
  for(i=nbinxy-10;i<nbinxy;i++) f2=f2+vectorBinXY[i]/10;
  averagebg=(f1+f2)/2;
  for(i=0;i<nbinxy;i++) {
    if(vectorBinXY[i]-averagebg>(max-averagebg)*0.4 && 
       vectorBinXY[i]-averagebg<(max-averagebg)*0.7) {
      sigma=hITSXv->GetBinCenter(binmax)-hITSXv->GetBinCenter(i);
      sigma=TMath::Abs(sigma);
      if(sigma==0) sigma=0.05;
    }
  }
  
  
  fx->SetParameter(0,max);
  fx->SetParameter(1,vzero[0]);
  fx->SetParameter(2,sigma); 
  fx->SetParameter(3,averagebg);
  
  hITSXv->Fit("fx","RQ0"); 
  
  snr[0] = fx->GetParameter(0)/fx->GetParameter(3);
  
  if(snr[0]<0.) { 
    Error("FindVertexForCurrentEvent","\nNegative Signal to noise ratio for x!\n");
    Error("FindVertexForCurrentEvent","The algorithm cannot find the x vertex position\n");
    //   exit(123456789);  
    return fCurrentVertex;
  }
  else {
    position[0]=fx->GetParameter(1);
    resolution[0]=fx->GetParError(1);
  }
  
  TF1 *fy = new TF1("fy","([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])))+[3]",vzero[1]-0.5,vzero[1]+0.5);  
  
  max=(Int_t) hITSYv->GetMaximum();
  binmax=hITSYv->GetMaximumBin();
  sigma=0;
  f1=f2=0;
  for(i=0;i<nbinxy;i++) vectorBinXY[i]=(Int_t)hITSYv->GetBinContent(i);
  for(i=0;i<10;i++) f1=f1+vectorBinXY[i]/10;
  for(i=nbinxy-10;i<nbinxy;i++) f2=f2+vectorBinXY[i]/10;
  averagebg=(f1+f2)/2;
  for(i=0;i<nbinxy;i++) {
    if(vectorBinXY[i]-averagebg>(max-averagebg)*0.4 && 
       vectorBinXY[i]-averagebg<(max-averagebg)*0.7) {
      sigma=hITSYv->GetBinCenter(binmax)-hITSYv->GetBinCenter(i);
      sigma=TMath::Abs(sigma);
      if(sigma==0) sigma=0.05;
    }
  }
  
  fy->SetParameter(0,max);
  fy->SetParameter(1,vzero[1]);
  fy->SetParameter(2,sigma); 
  fy->SetParameter(3,averagebg);
  
  hITSYv->Fit("fy","RQ0"); 
  
  snr[1] = fy->GetParameter(0)/fy->GetParameter(3);
  if(snr[1]<0.) {   
    Error("FindVertexForCurrentEvent","\nNegative Signal to noise ratio for y!\n");
    Error("FindVertexForCurrentEvent","The algorithm cannot find the y vertex position.\n");
    //    exit(123456789); 
    return fCurrentVertex;
  }
  else {
    position[1]=fy->GetParameter(1);
    resolution[1]=fy->GetParError(1);
  }
  
    
  /*cout << "iter = " << niter << " -> x = " << position[0] << endl;
    cout << "iter = " << niter << " -> y = " << position[1] << endl;
    cout << "iter = " << niter << " -> z = " << position[2] << endl;
    cout << "---\n";*/

  niter++;
   
  vzero[0] = position[0];
  vzero[1] = position[1];
  vzero[2] = position[2];

  if(niter<3) goto start;  // number of iterations
   

  cout<<"Number of iterations: "<<niter<<endl;
  cout << "\nNo. of Points on the two pixel layers: "<<np1<<" "<<np2<<endl; 
 
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
  delete [] vectorBinZ;
  delete [] vectorBinXY;
   
  delete hITSZv;
  delete hITSXv;
  delete hITSYv;
  Char_t name[30];
  //  sprintf(name,"Vertex_%d",evnumber);
  if(fDebug>0)Info("FindVertexForCurrentEvent","Vertex found for event %d",evnumber);
  sprintf(name,"Vertex");
  fCurrentVertex = new AliITSVertex(position,resolution,snr,name);
  return fCurrentVertex;
}

//______________________________________________________________________
Double_t AliITSVertexerIons::PhiFunc(Float_t p[]) {
  Double_t phi=0;

  if(p[1]>0 && p[0]>0) phi=(TMath::ATan((Double_t)(p[1]/p[0]))*57.29578);
  if(p[1]>0 && p[0]<0) phi=(TMath::ATan((Double_t)(p[1]/p[0]))*57.29578)+180;
  if(p[1]<0 && p[0]<0) phi=(TMath::ATan((Double_t)(p[1]/p[0]))*57.29578)+180;
  if(p[1]<0 && p[0]>0) phi=(TMath::ATan((Double_t)(p[1]/p[0]))*57.29578)+360;
  return phi;
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
