#include <TMath.h>
#include <TRandom.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream.h>

#include "AliRun.h"
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSRecPoint.h"
#include "AliGenerator.h"
#include "AliMagF.h"

#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <AliITSVertex.h>
#include <TObjArray.h>
#include <TObject.h>


ClassImp(AliITSVertex)

////////////////////////////////////////////////////////////////////////
// AliITSVertex is a class for primary vertex finding.
//
// Version: 0
// Written by Giuseppe Lo Re and Francesco Riggi
// Giuseppe.Lore@ct.infn.it
// Franco.Riggi@ct.infn.it
// Marh 2001
//
//
///////////////////////////////////////////////////////////////////////


//______________________________________________________________________________

AliITSVertex::AliITSVertex() {   

   fPosition = new Double_t[3];
   fResolution = new Double_t[3];
   fSNR = new Double_t[3];
  
   AliITS* aliits =(AliITS *)gAlice->GetDetector("ITS");
   AliITSgeom *g2 = ((AliITS*)aliits)->GetITSgeom(); 
   
   TClonesArray  *recpoints = aliits->RecPoints();
   AliITSRecPoint *pnt;
   TRandom rnd;

   Int_t NoPoints1 = 0;						
   Int_t NoPoints2 = 0;
   Double_t Vzero[3];
   Double_t AsPar[2];
     
//------------Rough Vertex------------------------------------------------------   

   Int_t i,npoints,ipoint,j,k;
   Double_t ZCentroid;
   Float_t l[3], p[3];

   TH1F *hITSx1pos      = new TH1F("hITSx1pos","",100,0.,4.2);
   TH1F *hITSx1neg      = new TH1F("hITSx1neg","",100,-4.2,0.);
   TH1F *hITSy1pos      = new TH1F("hITSy1pos","",100,0.,4.2);
   TH1F *hITSy1neg      = new TH1F("hITSy1neg","",100,-4.2,0.);
   TH1F *hITSz1         = new TH1F("hITSz1","",100,-14.35,14.35);
  
   for(i=g2->GetStartSPD();i<=g2->GetLastSPD();i++) 
   {
	aliits->ResetRecPoints(); 
        gAlice->TreeR()->GetEvent(i);	
	
	   npoints = recpoints->GetEntries();
	   for (ipoint=0;ipoint<npoints;ipoint++) {
	       pnt = (AliITSRecPoint*)recpoints->UncheckedAt(ipoint);
	       l[0]=pnt->GetX();
	       l[1]=0;
	       l[2]=pnt->GetZ();
		   g2->LtoG(i, l, p);
		   if(i<80 && TMath::Abs(p[2])<14.35) { 	
		        if(p[0]>0) hITSx1pos->Fill(p[0]); 
		        if(p[0]<0) hITSx1neg->Fill(p[0]); 
		        if(p[1]>0) hITSy1pos->Fill(p[1]); 
		        if(p[1]<0) hITSy1neg->Fill(p[1]); 
		   hITSz1->Fill(p[2]);       
	       }
		   if(i>=80 && TMath::Abs(p[2])<14.35) (NoPoints2)++;
       }       
   }  
  
   NoPoints1 = (Int_t)(hITSz1->GetEntries()); 
   	   
   Double_t mxpiu = hITSx1pos->GetEntries();
   Double_t mxmeno = hITSx1neg->GetEntries();
   Double_t mypiu = hITSy1pos->GetEntries();
   Double_t mymeno = hITSy1neg->GetEntries();

   AsPar[0] = (mxpiu-mxmeno)/(mxpiu+mxmeno);
   AsPar[1] = (mypiu-mymeno)/(mypiu+mymeno); 
   
   Vzero[0] = 5.24434*AsPar[0]; 
   Vzero[1] = 5.24434*AsPar[1]; 

   ZCentroid=  TMath::Abs(hITSz1->GetMean()); 
   Vzero[2] = -5.31040e-02+1.42198*ZCentroid+7.44718e-01*TMath::Power(ZCentroid,2)
   	          -5.73426e-01*TMath::Power(ZCentroid,3)+2.01500e-01*TMath::Power(ZCentroid,4)		
              -3.34118e-02*TMath::Power(ZCentroid,5)+2.20816e-03*TMath::Power(ZCentroid,6);
   
   if(hITSz1->GetMean()<0) Vzero[2] = -Vzero[2];
   
   //cout << "\nCentroid and RMS of hIITSz1: \n";
   //cout << hITSz1->GetMean() << "       "  << hITSz1->GetRMS()  <<endl;
   //cout << "\nAsPar[0]: " << AsPar[0];
   //cout << "\nAsPar[1]: " << AsPar[1];
   //cout << "\nAsPar[2]: " << AsPar[2];
   cout << "\nXvzero: " << Vzero[0] << " cm" << "";
   cout << "\nYvzero: " << Vzero[1] << " cm" << "";
   cout << "\nZvzero: " << Vzero[2] << " cm" << "\n";

   delete hITSx1pos;
   delete hITSx1neg;
   delete hITSy1pos;
   delete hITSy1neg;
   delete hITSz1;

//-----------------------Get points---------------------------------------------

   Double_t dphi,DeltaPhi,DeltaZ,r;
   Int_t np1=0;
   Int_t np2=0;

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
	
   for(i=g2->GetStartSPD();i<=g2->GetLastSPD();i++) 
   {
	aliits->ResetRecPoints(); 
	gAlice->TreeR()->GetEvent(i);
	   npoints = recpoints->GetEntries();
	   for (ipoint=0;ipoint<npoints;ipoint++) {
	       
	        //if(rnd.Integer(10)>2) continue;
		    pnt = (AliITSRecPoint*)recpoints->UncheckedAt(ipoint);
	        l[0]=pnt->GetX();
	        l[1]=0;
	        l[2]=pnt->GetZ();
	        g2->LtoG(i, l, p);
		    p[0]=p[0]-Vzero[0];
		    p[1]=p[1]-Vzero[1];
	        r=TMath::Sqrt(TMath::Power(p[0],2)+TMath::Power(p[1],2));
	      
			if(i<80 && TMath::Abs(p[2])<14.35) { 
			Z1[np1]=p[2]; 
			Y1[np1]=p[1];
			X1[np1]=p[0];
			r1[np1]=r;
			phi1[np1]=PhiFunc(p);
			np1++;
			}

			if(i>=80 &&  TMath::Abs(p[2])<14.35) { 
			Z2[np2]=p[2];
			Y2[np2]=p[1];
			X2[np2]=p[0];
			r2[np2]=r;
			phi2[np2]=PhiFunc(p);
	 		np2++;
	 		}
			
	  }
   }

//-------------------Correlations-----------------------------------------------
    
   //cout << "\nPoints number on the first and second layer: "<<np1<<" "<<np2<<endl;

   DeltaPhi = 0.085-0.1*TMath::Sqrt(Vzero[0]*Vzero[0]+Vzero[1]*Vzero[1]); 
   if(DeltaPhi<0) {
      cout << "\nThe algorith can't find the vertex. " << endl;
	  exit(123456789);
   }
   
   Float_t B[3];
   Float_t origin[3];
   for(Int_t okm=0;okm<3;okm++) origin[okm]=0;
   gAlice->Field()->Field(origin,B);

   DeltaPhi = DeltaPhi*B[2]/2;

   cout << "\nDeltaPhi: " << DeltaPhi << " deg" << "\n";   

   DeltaZ =
   0.2+2.91617e-01+2.67567e-02*Vzero[2]+1.49626e-02*TMath::Power(Vzero[2],2); 
     
   cout << "DeltaZeta: " << DeltaZ << " cm" << "\n";     

   TH1F *hITSZv         = new TH1F("hITSZv","",DeltaZ/0.005,Vzero[2]-DeltaZ,Vzero[2]+DeltaZ);

   for(j=0; j<(np2)-1; j++) {
	 for(k=0; k<(np1)-1; k++) { 
		if(TMath::Abs(Z1[k]-Z2[j])>13.) continue;
		dphi=TMath::Abs(phi1[k]-phi2[j]);
		if(dphi>180) dphi = 360-dphi;
		if(dphi<DeltaPhi) { 
		  if(TMath::Abs((Z1[k]-(Z2[j]-Z1[k])/((r2[j]/r1[k])-1))-Vzero[2])<DeltaZ) 
		      hITSZv->Fill(Z1[k]-(Z2[j]-Z1[k])/((r2[j]/r1[k])-1));
		}
	 }
   }
   
   
   
   cout << "\nNumber of used pairs: \n";
   cout << hITSZv->GetEntries() << '\n' << '\n';	
   cout << "\nCentroid and RMS of hIITSZv: \n";
   cout <<  hITSZv->GetMean() << "       " << hITSZv->GetRMS() << "\n"<< "\n";

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

   Double_t a = Vzero[2]-DeltaZ; 
   Double_t b = Vzero[2]+DeltaZ; 
   
   TF1 *f4 = new TF1 ("f4","([0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2])))+[3]",a,b);  
   f4->SetParameter(0,700.);
   f4->SetParameter(1,Vzero[2]);
   f4->SetParameter(2,0.05); // !!
   f4->SetParameter(3,50.);
   
   hITSZv->Fit("f4","RME0"); 
   
   delete hITSZv;
  
   fSNR[2] = f4->GetParameter(0)/f4->GetParameter(3);
   if(fSNR[2]<1.5) { 
     cout << "\nSignal to noise ratio too small!!!" << endl;
     cout << "The algorithm can't find the z vertex position." << endl;
	 exit(123456789);
   }
   else
   {
     fPosition[2] = (f4->GetParameter(1));
     if(fPosition[2]<0) 
       {
         fPosition[2]=fPosition[2]-TMath::Abs(fPosition[2])*1.11/10000;
       }
     else
       {
         fPosition[2]=fPosition[2]+TMath::Abs(fPosition[2])*1.11/10000;
       }
   }
   fResolution[2] = f4->GetParError(1);
   
   AliGenerator *gener = gAlice->Generator();
   Float_t Vx,Vy,Vz;
   gener->GetOrigin(Vx,Vy,Vz);
 
   fPosition[0]=(Double_t)Vx;
   fPosition[1]=(Double_t)Vy;

   fResolution[0] = 0;
   fResolution[1] = 0;

   fSNR[0] = 0;
   fSNR[1] = 0;
   
}


//______________________________________________________________________________

AliITSVertex::~AliITSVertex() { 
   delete [] fPosition;
   delete [] fResolution;
   delete [] fSNR;
}

//______________________________________________________________________________


Double_t AliITSVertex::PhiFunc(Float_t p[]) {
	 Double_t phi=0;
	 if(p[1]>0 && p[0]>0) phi=(TMath::ATan((Double_t)(p[1]/p[0]))*57.29578);
	 if(p[1]>0 && p[0]<0) phi=(TMath::ATan((Double_t)(p[1]/p[0]))*57.29578)+180;
	 if(p[1]<0 && p[0]<0) phi=(TMath::ATan((Double_t)(p[1]/p[0]))*57.29578)+180;
	 if(p[1]<0 && p[0]>0) phi=(TMath::ATan((Double_t)(p[1]/p[0]))*57.29578)+360;
	 return phi;
}
//
//______________________________________________________________________________

