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
//////////////////////////////////////////////////////////////////////
// AliITSVertexerIons is a class for full 3D primary vertex         //
// finding optimized for Ion-Ion interactions                       //
//                                                                  // 
//                                                                  //
//                                                                  //
//                                                                  //
// Written by Giuseppe Lo Re and Francesco Riggi                    //
// Giuseppe.Lore@ct.infn.it                                         //
//                                                                  //
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
}

//______________________________________________________________________
AliITSVertexerIons::AliITSVertexerIons(TString fn):AliITSVertexer(fn) {
  // Standard constructor
  
  fITS = 0;
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
  Double_t Position[3];
  Double_t Resolution[3];
  Double_t SNR[3];
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
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

  Int_t NoPoints1 = 0;                                         
  Int_t NoPoints2 = 0;
  Double_t Vzero[3];
  Double_t AsPar[2];
     
  //------------ Rough Vertex evaluation ---------------------------------   

  Int_t i,npoints,ipoint,j,k,max,BinMax;
  Double_t ZCentroid;
  Float_t l[3], p[3];

  Double_t mxpiu = 0;
  Double_t mxmeno = 0;
  Double_t mypiu = 0;
  Double_t mymeno = 0;
   
  TH1F *hITSz1         = new TH1F("hITSz1","",100,-14.35,14.35);
  TTree *TR =  ITSloader->TreeR(); 
  for(i=g2->GetStartSPD();i<=g2->GetLastSPD();i++) 
    {
      fITS->ResetRecPoints();
      TR->GetEvent(i);    
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
	if(i>=80 && TMath::Abs(p[2])<14.35) NoPoints2++;
      }       
    }  
  
  NoPoints1 = (Int_t)(hITSz1->GetEntries()); 
           
  AsPar[0] = (mxpiu-mxmeno)/(mxpiu+mxmeno);
  AsPar[1] = (mypiu-mymeno)/(mypiu+mymeno); 
   
  Vzero[0] = 5.24441*AsPar[0]; 
  Vzero[1] = 5.24441*AsPar[1]; 

  ZCentroid=  TMath::Abs(hITSz1->GetMean()); 
  Vzero[2] = -5.31040e-02+1.42198*ZCentroid+7.44718e-01*TMath::Power(ZCentroid,2)
    -5.73426e-01*TMath::Power(ZCentroid,3)+2.01500e-01*TMath::Power(ZCentroid,4)          
    -3.34118e-02*TMath::Power(ZCentroid,5)+2.20816e-03*TMath::Power(ZCentroid,6);
   
  if(hITSz1->GetMean()<0) Vzero[2] = -Vzero[2];
   
  /*cout << "\nXvzero: " << Vzero[0] << " cm" << "";
    cout << "\nYvzero: " << Vzero[1] << " cm" << "";
    cout << "\nZvzero: " << Vzero[2] << " cm" << "\n";*/

  delete hITSz1;

  Double_t dphi,r,DeltaZ,a,b;
  Int_t np1=0;
  Int_t np2=0;
  Int_t niter=0;
   
  Double_t DeltaPhiZ = 0.08;
   
  Float_t B[3];
  Float_t origin[3];
  for(Int_t okm=0;okm<3;okm++) origin[okm]=0;
  gAlice->Field()->Field(origin,B);

  DeltaPhiZ = DeltaPhiZ*B[2]/2;
  Double_t DeltaPhiXY = 1.0;   

  //   cout << "\nDeltaPhiZ: " << DeltaPhiZ << " deg" << "\n";   
  //   cout << "DeltaPhiXY: " << DeltaPhiXY << " deg" << "\n";   
   
  Double_t *Z1, *Z2, *Y1, *Y2, *X1, *X2, *phi1, *phi2, *r1, *r2;
  Z1=new Double_t[NoPoints1];
  Z2=new Double_t[NoPoints2];
  Y1=new Double_t[NoPoints1];
  Y2=new Double_t[NoPoints2];
  X1=new Double_t[NoPoints1];
  X2=new Double_t[NoPoints2];
  phi1=new Double_t[NoPoints1];
  phi2=new Double_t[NoPoints2];
  r1=new Double_t[NoPoints1];
  r2=new Double_t[NoPoints2];
        
  DeltaZ = 4.91617e-01+2.67567e-02*Vzero[2]+1.49626e-02*TMath::Power(Vzero[2],2); 
  Float_t MulFactorZ = 28000./(Float_t)NoPoints1;
  Int_t nbin=(Int_t)((DeltaZ/0.005)/MulFactorZ);
  Int_t nbinxy=250;
  Int_t *VectorBinZ,*VectorBinXY;
  VectorBinZ=new Int_t[nbin];
  VectorBinXY=new Int_t[nbinxy];
  Float_t f1= 0;
  Float_t f2= 0;
  Double_t sigma,MediaFondo;
     
  TH1D *hITSZv         = new TH1D("hITSZv","",nbin,Vzero[2]-DeltaZ,Vzero[2]+DeltaZ);
  TH1D *hITSXv         = new TH1D("hITSXv","",nbinxy,-3,3);
  TH1D *hITSYv         = new TH1D("hITSYv","",nbinxy,-3,3);

  //   cout << "DeltaZeta: " << DeltaZ << " cm" << "\n";   
   
   
 start:   

  hITSZv->Add(hITSZv,-1.);
  hITSXv->Add(hITSXv,-1.);
  hITSYv->Add(hITSYv,-1.);

  np1=np2=0;

   

  for(i=g2->GetStartSPD();i<=g2->GetLastSPD();i++) 
    {
      fITS->ResetRecPoints(); 
      ITSloader->TreeR()->GetEvent(i);
      npoints = recpoints->GetEntries();
      for (ipoint=0;ipoint<npoints;ipoint++) {
               
	pnt = (AliITSRecPoint*)recpoints->UncheckedAt(ipoint);
	l[0]=pnt->GetX();
	l[1]=0;
	l[2]=pnt->GetZ();
	g2->LtoG(i, l, p);
              
	if(i<80 && TMath::Abs(p[2])<14.35) {
	  p[0]=p[0]-Vzero[0];
	  p[1]=p[1]-Vzero[1];
	  r=TMath::Sqrt(TMath::Power(p[0],2)+TMath::Power(p[1],2));
	  Y1[np1]=p[1];
	  X1[np1]=p[0];
	  Z1[np1]=p[2]; 
	  r1[np1]=r;
	  phi1[np1]=PhiFunc(p);
	  np1++;
	}

	if(i>=80 &&  TMath::Abs(p[2])<14.35) { 
	  p[0]=p[0]-Vzero[0];
	  p[1]=p[1]-Vzero[1];
	  r=TMath::Sqrt(TMath::Power(p[0],2)+TMath::Power(p[1],2));
	  Y2[np2]=p[1];
	  X2[np2]=p[0];
	  Z2[np2]=p[2];
	  r2[np2]=r;
	  phi2[np2]=PhiFunc(p);
	  np2++;
	}
                        
      }
    }

  //------------------ Correlation between rec points ----------------------
    
  //cout << "\nNo. of Points on the two pixel layers: "<<np1<<" "<<np2<<endl; 

  for(j=0; j<(np2)-1; j++) {
    for(k=0; k<(np1)-1; k++) { 
      dphi=TMath::Abs(phi1[k]-phi2[j]);
      if(dphi>180) dphi = 360-dphi;
      if(dphi<DeltaPhiZ && TMath::Abs((Z1[k]-(Z2[j]-Z1[k])/((r2[j]/r1[k])-1))-Vzero[2])
	 <DeltaZ) hITSZv->Fill(Z1[k]-(Z2[j]-Z1[k])/((r2[j]/r1[k])-1));        
    }
  }
      
  //cout << "\nNumber of used pairs: \n";
  //cout << hITSZv->GetEntries() << '\n' << '\n';      
  a = Vzero[2]-DeltaZ; 
  b = Vzero[2]+DeltaZ; 
  max=(Int_t) hITSZv->GetMaximum();
  BinMax=hITSZv->GetMaximumBin();
  sigma=0;
  for(i=0;i<nbin;i++) VectorBinZ[i]=(Int_t)hITSZv->GetBinContent(i);
  for(i=0;i<10;i++) f1=f1+VectorBinZ[i]/10;
  for(i=nbin-10;i<nbin;i++) f2=f2+VectorBinZ[i]/10;
  MediaFondo=(f1+f2)/2;
  for(i=0;i<nbin;i++) {
    if(VectorBinZ[i]-MediaFondo>(max-MediaFondo)*0.4 && 
       VectorBinZ[i]-MediaFondo<(max-MediaFondo)*0.7) {
      sigma=hITSZv->GetBinCenter(BinMax)-hITSZv->GetBinCenter(i);
      sigma=TMath::Abs(sigma);
      if(sigma==0) sigma=0.05;
    }
  }
   
  /*cout << "f1 " <<f1 <<endl;
    cout << "f2 " <<f2 <<endl;
    cout << "GetMaximumBin " <<hITSZv->GetMaximumBin() <<endl;
    cout << "nbin " <<nbin <<endl;
    cout << "max " << hITSZv->GetBinContent(BinMax)<<endl;
    cout << "sigma " <<sigma <<endl;
    cout << "Fondo " <<   MediaFondo<< endl;*/
   
  TF1 *fz = new TF1 ("fz","([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])))+[3]",a,b);  
   
  fz->SetParameter(0,max);
  if(niter==0) {Double_t temp = hITSZv->GetBinCenter(BinMax); Vzero[2]=temp;}
  fz->SetParameter(1,Vzero[2]);
  fz->SetParameter(2,sigma); 
  fz->SetParameter(3,MediaFondo);  
  fz->SetParLimits(2,0,999);
  fz->SetParLimits(3,0,999);
   
  hITSZv->Fit("fz","RMEQ0");
     
  SNR[2] = fz->GetParameter(0)/fz->GetParameter(3);
  if(SNR[2]<0.) { 
    cout << "\nNegative Signal to noise ratio for z!!!" << endl;
    cout << "The algorithm cannot find the z vertex position." << endl;
    exit(123456789);
  }
  else
    {
      Position[2] = fz->GetParameter(1);
      if(Position[2]<0) 
	{
	  Position[2]=Position[2]-TMath::Abs(Position[2])*1.11/10000;
	}
      else
	{
	  Position[2]=Position[2]+TMath::Abs(Position[2])*1.11/10000;
	}
    }
  Resolution[2] = fz->GetParError(1);
   
   
   
  for(j=0; j<(np2)-1; j++) {
    for(k=0; k<(np1)-1; k++) { 
      dphi=TMath::Abs(phi1[k]-phi2[j]);
      if(dphi>180) dphi = 360-dphi;
      if(dphi>DeltaPhiXY) continue;
      if(TMath::Abs((Z1[k]-(Z2[j]-Z1[k])/((r2[j]/r1[k])-1))-Position[2])
	 <4*Resolution[2]) {
	hITSXv->Fill(Vzero[0]+(X2[j]-((X2[j]-X1[k])/(Y2[j]-Y1[k]))*Y2[j]));
	hITSYv->Fill(Vzero[1]+(Y2[j]-((Y2[j]-Y1[k])/(X2[j]-X1[k]))*X2[j]));
      }
    }
  }
   
  TF1 *fx = new TF1
    ("fx","([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])))+[3]",Vzero[0]-0.5,Vzero[0]+0.5);  

  max=(Int_t) hITSXv->GetMaximum();
  BinMax=hITSXv->GetMaximumBin();
  sigma=0;
  f1=f2=0;
  for(i=0;i<nbinxy;i++) VectorBinXY[i]=(Int_t)hITSXv->GetBinContent(i);
  for(i=0;i<10;i++) f1=f1+VectorBinXY[i]/10;
  for(i=nbinxy-10;i<nbinxy;i++) f2=f2+VectorBinXY[i]/10;
  MediaFondo=(f1+f2)/2;
  for(i=0;i<nbinxy;i++) {
    if(VectorBinXY[i]-MediaFondo>(max-MediaFondo)*0.4 && 
       VectorBinXY[i]-MediaFondo<(max-MediaFondo)*0.7) {
      sigma=hITSXv->GetBinCenter(BinMax)-hITSXv->GetBinCenter(i);
      sigma=TMath::Abs(sigma);
      if(sigma==0) sigma=0.05;
    }
  }
   
  /*cout << "f1 " <<f1 <<endl;
    cout << "f2 " <<f2 <<endl;
    cout << "GetMaximumBin " <<hITSXv->GetMaximumBin() <<endl;
    cout << "max " << hITSXv->GetBinContent(BinMax)<<endl;
    cout << "sigma " <<sigma <<endl;
    cout << "Fondo " <<   MediaFondo<< endl;
    cout << "nbinxy " <<nbinxy <<endl;*/
   
  fx->SetParameter(0,max);
  fx->SetParameter(1,Vzero[0]);
  fx->SetParameter(2,sigma); 
  fx->SetParameter(3,MediaFondo);

  hITSXv->Fit("fx","RMEQ0"); 

  SNR[0] = fx->GetParameter(0)/fx->GetParameter(3);
  if(SNR[0]<0.) { 
    cout << "\nNegative Signal to noise ratio for x!!!" << endl;
    cout << "The algorithm cannot find the x vertex position." << endl;
    exit(123456789);  
  }
  else
    {
      Position[0]=fx->GetParameter(1);
      Resolution[0]=fx->GetParError(1);
    }

  TF1 *fy = new TF1("fy","([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])))+[3]",Vzero[1]-0.5,Vzero[1]+0.5);  

  max=(Int_t) hITSYv->GetMaximum();
  BinMax=hITSYv->GetMaximumBin();
  sigma=0;
  f1=f2=0;
  for(i=0;i<nbinxy;i++) VectorBinXY[i]=(Int_t)hITSYv->GetBinContent(i);
  for(i=0;i<10;i++) f1=f1+VectorBinXY[i]/10;
  for(i=nbinxy-10;i<nbinxy;i++) f2=f2+VectorBinXY[i]/10;
  MediaFondo=(f1+f2)/2;
  for(i=0;i<nbinxy;i++) {
    if(VectorBinXY[i]-MediaFondo>(max-MediaFondo)*0.4 && 
       VectorBinXY[i]-MediaFondo<(max-MediaFondo)*0.7) {
      sigma=hITSYv->GetBinCenter(BinMax)-hITSYv->GetBinCenter(i);
      sigma=TMath::Abs(sigma);
      if(sigma==0) sigma=0.05;
    }
  }
   
  /*cout << "f1 " <<f1 <<endl;
    cout << "f2 " <<f2 <<endl;
    cout << "GetMaximumBin " <<hITSYv->GetMaximumBin() <<endl;
    cout << "max " << hITSYv->GetBinContent(BinMax)<<endl;
    cout << "sigma " <<sigma <<endl;
    cout << "Fondo " <<   MediaFondo<< endl;
    cout << "nbinxy " <<nbinxy <<endl;*/
   
  fy->SetParameter(0,max);
  fy->SetParameter(1,Vzero[1]);
  fy->SetParameter(2,sigma); 
  fy->SetParameter(3,MediaFondo);
   
  hITSYv->Fit("fy","RMEQ0"); 

  SNR[1] = fy->GetParameter(0)/fy->GetParameter(3);
  if(SNR[1]<0.) { 
    cout << "\nNegative Signal to noise ratio for y!!!" << endl;
    cout << "The algorithm cannot find the y vertex position." << endl;
    exit(123456789); 
  }
  else
    {
      Position[1]=fy->GetParameter(1);
      Resolution[1]=fy->GetParError(1);
    }

  niter++;
   
  Vzero[0] = Position[0];
  Vzero[1] = Position[1];
  Vzero[2] = Position[2];

  if(niter<3) goto start;  // number of iterations
   

  cout<<"Number of iterations: "<<niter<<endl;
 
  delete [] Z1;
  delete [] Z2;
  delete [] Y1;
  delete [] Y2;
  delete [] X1;
  delete [] X2;
  delete [] r1;
  delete [] r2;
  delete [] phi1;
  delete [] phi2;
  delete [] VectorBinZ;
  delete [] VectorBinXY;
   
  delete hITSZv;
  delete hITSXv;
  delete hITSYv;
  Char_t name[30];
  //  sprintf(name,"Vertex_%d",evnumber);
  if(fDebug>0)Info("FindVertexForCurrentEvent","Vertex found for event %d",evnumber);
  sprintf(name,"Vertex");
  fCurrentVertex = new AliITSVertex(Position,Resolution,SNR,name);
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
  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  ITSloader->LoadRecPoints("read");
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
