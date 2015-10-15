/**************************************************************************
 * Copyright(c) 2006-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------
//           AliITSVertexer3DTapan class
//   This is a class for the 3d vertex finding
//    Origin: Tapan Nayak, VECC-CERN, Tapan.Nayak@cern.ch
//-----------------------------------------------------------------

#include <TH1.h>
#include <TTree.h>
#include <TClonesArray.h>

#include "AliITSVertexer3DTapan.h"
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"
#include "AliESDVertex.h"
#include "AliLog.h"

ClassImp(AliITSVertexer3DTapan)

void AliITSVertexer3DTapan::LoadClusters(TTree *cTree) {
  //--------------------------------------------------------------------
  //This function loads the SPD clusters
  //--------------------------------------------------------------------
   TClonesArray dummy("AliITSRecPoint",10000), *clusters=&dummy;
   TBranch *branch=cTree->GetBranch("ITSRecPoints");
   branch->SetAddress(&clusters);

   Int_t nentr=cTree->GetEntries(),nc1=0,nc2=0;
   for (Int_t i=0; i<nentr; i++) {
       if (!cTree->GetEvent(i)) continue;
       //
       //   Below:
       //   "alpha" is the angle from the global X-axis to the
       //    local GEANT X'-axis  ( rot[0]=cos(alpha) and rot[1]=sin(alpha) )
       //    "phi" is the angle from the global X-axis to the
       //         local cluster X"-axis
       //

       Int_t    lay,lad,det; AliITSgeomTGeo::GetModuleId(i,lay,lad,det);

       if (lay>2) break;  //load the SPD clusters only

       Int_t ncl=clusters->GetEntriesFast();
       Float_t hPhi=0.;
       while (ncl--) {
          AliITSRecPoint *c=(AliITSRecPoint*)clusters->UncheckedAt(ncl);
	  Float_t pos[3];
	  c->GetGlobalXYZ(pos);
          if (lay==1) {
	    /*             fX1[nc1]= r*cp - c->GetY()*sp;
             fY1[nc1]= r*sp + c->GetY()*cp;
             fZ1[nc1]= c->GetZ(); */
	    fX1[nc1] = pos[0]; fY1[nc1] = pos[1]; fZ1[nc1] = pos[2];
	     CalculatePhi(fX1[nc1], fY1[nc1], hPhi);
	     fPhi1[nc1]= hPhi;
             nc1++;
          } else {
            /* fX2[nc2]= r*cp - c->GetY()*sp;
             fY2[nc2]= r*sp + c->GetY()*cp;
             fZ2[nc2]= c->GetZ(); */
	    fX2[nc2] = pos[0]; fY2[nc2] = pos[1]; fZ2[nc2] = pos[2];
	     CalculatePhi(fX2[nc2], fY2[nc2], hPhi);
	     fPhi2[nc2]= hPhi;
             nc2++;
          }
       }
   }
   ficlu1 = nc1; ficlu2 = nc2;
   AliInfo(Form("Number of clusters: %d (first layer) and %d (second layer)",ficlu1,ficlu2));
}

AliESDVertex *AliITSVertexer3DTapan::FindVertexForCurrentEvent(TTree *cTree) {
  //
  // This function reconstructs ....
  //
  //
  LoadClusters(cTree);

  Double_t pos[3], postemp[3];
  Double_t sigpos[3]={0.,0.,0.};
  Int_t ncontr, ncontrtemp;
  Float_t cuts[3];
  Int_t vtxstatus=0;

  //....
  pos[0] = 0.;   pos[1] = 0.;   pos[2] = 0.;
  cuts[0]=1.;  cuts[1]=1.;  cuts[2]=20.;
  CalculateVertex3d1(pos, cuts, ncontr);
  if(ncontr==0) {
    pos[0] = 9999.;   pos[1] = 9999.;   pos[2] = 9999.;
    vtxstatus = -1;
  }
  AliInfo(Form("1st step: %d %f %f %f st=%d",ncontr,pos[0],pos[1],pos[2],vtxstatus));

  if(vtxstatus == 0) {
    ncontrtemp = ncontr; postemp[0] = pos[0];   postemp[1] = pos[1];   postemp[2] = pos[2]; 
    cuts[0]=0.3;  cuts[1]=0.3;  cuts[2]=1.;
    CalculateVertex3d1(pos, cuts, ncontr);
    if(ncontr==0) {
      ncontr = ncontrtemp; pos[0] = postemp[0];   pos[1] = postemp[1];   pos[2] = postemp[2];
      vtxstatus = 2;
    }
    AliInfo(Form("2nd step: %d %f %f %f st=%d",ncontr,pos[0],pos[1],pos[2],vtxstatus));
  }

  if(vtxstatus == 0) {
    ncontrtemp = ncontr; postemp[0] = pos[0];   postemp[1] = pos[1];   postemp[2] = pos[2]; 
    cuts[0]=0.25;  cuts[1]=0.25;  cuts[2]=1.0;
    CalculateVertex3d2(pos, cuts, ncontr, sigpos);
    if(ncontr==0) {
      ncontr = ncontrtemp; pos[0] = postemp[0];   pos[1] = postemp[1];   pos[2] = postemp[2];
      vtxstatus = 3;
    }
    AliInfo(Form("3rd step: %d %f %f %f st=%d",ncontr,pos[0],pos[1],pos[2],vtxstatus));
  }

  if(vtxstatus == 0) {
    ncontrtemp = ncontr; postemp[0] = pos[0];   postemp[1] = pos[1];   postemp[2] = pos[2]; 
    cuts[0]=0.1;  cuts[1]=0.1;  cuts[2]=0.2;
    CalculateVertex3d2(pos, cuts, ncontr, sigpos);
    if(ncontr==0) {
      ncontr = ncontrtemp; pos[0] = postemp[0];   pos[1] = postemp[1];   pos[2] = postemp[2];
      vtxstatus = 4;
    }
    AliInfo(Form("4th step: %d %f %f %f st=%d",ncontr,pos[0],pos[1],pos[2],vtxstatus));
  }
  AliInfo(Form("Final step: %d %f %f %f st=%d",ncontr,pos[0],pos[1],pos[2],vtxstatus));

  Double_t covma[6]={0.,0.,0.,0.,0.,0.};
  covma[0]=sigpos[0];
  covma[2]=sigpos[1];
  covma[5]=sigpos[2];
  return new AliESDVertex(pos,covma,(Double_t)vtxstatus,ncontr,"AliITSVertexer3DTapan");

}


void AliITSVertexer3DTapan::CalculateVertex3d1(Double_t pos[3], Float_t cuts[3], Int_t &ncontr) {
  //
  // This function reconstructs first two steps of vertex
  //

  Double_t p1[4], p2[4], p3[4], p4[4];
  Double_t pa[3], pb[3];
  Double_t hphi1, hphi2, hphi3, hphi4;

  ncontr = 0;
  Float_t  phicut = 1.0;
  Double_t distance;  Float_t  distancecut = 1.0;
  Int_t    ibin=20;  Float_t  ilow=-1.; Float_t  ihigh=1.; 
  Int_t    ibinz=400; Float_t  ilowz=-20.; Float_t  ihighz=20.;

  TH1F *hx = new TH1F("hx","",  ibin, ilow, ihigh);
  TH1F *hy = new TH1F("hy","",  ibin, ilow, ihigh);
  TH1F *hz = new TH1F("hz","",  ibinz,ilowz,ihighz);
  
  for (Int_t ip1=0; ip1<ficlu1; ip1++) {
    // Two points on layer1: p1 and p3
    p1[0] = fX1[ip1];   p1[1] = fY1[ip1];    p1[2] = fZ1[ip1];   
    p3[0] = fX1[ip1+1]; p3[1] = fY1[ip1+1];  p3[2] = fZ1[ip1+1]; 
    hphi1 = fPhi1[ip1]; hphi3 = fPhi1[ip1+1];
    
    for (Int_t ip2=0; ip2<ficlu2; ip2++) {
      // Two points on layer 2: p2 and p4
      p2[0] = fX2[ip2];   p2[1] = fY2[ip2];     p2[2] = fZ2[ip2]; 
      p4[0] = fX2[ip2+1]; p4[1] = fY2[ip2+1];   p4[2] = fZ2[ip2+1]; 
      hphi2 = fPhi2[ip2]; hphi4 = fPhi2[ip2+1];

      // First line is formed by p1-p2 and second line by p3-p4
      // We find two points on each line which form the closest distance of the two lines
      // pa[0],pa[1],pa[2]: points on line 1 and pb[0],pb[1],pb[2]: points on line 2
      // Next: Consider x, y and z to be less than cuts[0], cuts[1] and cuts[2], respectively

      if(TMath::Abs(hphi1-hphi2)<phicut && TMath::Abs(hphi3-hphi4)<phicut){
	CalculateLine(p1, p2, p3, p4, pa, pb);

	if (pa[0]>pos[0]-cuts[0] && pa[0]<pos[0]+cuts[0] && pa[1]>pos[1]-cuts[1] && pa[1]<pos[1]+cuts[1] && pa[2]>pos[2]-cuts[2] && pa[2]<pos[2]+cuts[2]){
	  distance = (TMath::Sqrt(pow((pa[0]-pb[0]),2) + pow((pa[1]-pb[1]),2) + pow((pa[2]-pb[2]),2)));
	  if(distance<distancecut){
	    hx->Fill(pa[0]); 		hy->Fill(pa[1]);		hz->Fill(pa[2]);
	    hx->Fill(pb[0]); 		hy->Fill(pb[1]);		hz->Fill(pb[2]);
	    ncontr++;
	  }
	}
      }

      // Third line is formed by p1-p4 and fourth line by p3-p2
      // We find two points on each line which form the closest distance of the two lines
      // pa[0],pa[1],pa[2]: points on line 3 and pb[0],pb[1],pb[2]: points on line 4
      // Next: Consider x, y and z to be less than cuts[0], cuts[1] and cuts[2], respectively      
      if(TMath::Abs(hphi1-hphi4)<phicut && TMath::Abs(hphi3-hphi2)<phicut) {
	CalculateLine(p1, p4, p3, p2, pa, pb);
	if (pa[0]>pos[0]-cuts[0] && pa[0]<pos[0]+cuts[0] && pa[1]>pos[1]-cuts[1] && pa[1]<pos[1]+cuts[1]){
	  distance = (TMath::Sqrt(pow((pa[0]-pb[0]),2) + pow((pa[1]-pb[1]),2) + pow((pa[2]-pb[2]),2)));
	  if(distance<distancecut){
	    hx->Fill(pa[0]); 		hy->Fill(pa[1]);		hz->Fill(pa[2]);
	    hx->Fill(pb[0]); 		hy->Fill(pb[1]);		hz->Fill(pb[2]);
	    ncontr++;
	  }
	}
      }
    }
  }

  Int_t maxbinx = hx->GetMaximumBin(); 
  Int_t maxbiny = hy->GetMaximumBin();
  Int_t maxbinz = hz->GetMaximumBin();
  pos[0] = ilow  + ((ihigh-ilow)/ibin)*maxbinx;
  pos[1] = ilow  + ((ihigh-ilow)/ibin)*maxbiny;
  pos[2] = ilowz + ((ihighz-ilowz)/ibinz)*maxbinz;
  hx->Delete();
  hy->Delete();
  hz->Delete();
}

void AliITSVertexer3DTapan::CalculateVertex3d2(Double_t pos[3], Float_t cuts[3], Int_t &ncontr, Double_t sigpos[3]) {
  //
  // This function reconstructs second two steps of vertex
  //

  Double_t p1[4], p2[4], p3[4], p4[4];
  Double_t pa[3], pb[3];
  Double_t hphi1, hphi2, hphi3, hphi4;

  ncontr = 0;
  Float_t  phicut = 0.3;
  Double_t distance; Float_t  distancecut = 1.0;

  Double_t vertx  =0.;   Double_t verty  =0.;   Double_t vertz  =0.;	 
  Double_t vertx2 =0.;   Double_t verty2 =0.;   Double_t vertz2 =0.;	 

  for (Int_t ip1=0; ip1<ficlu1; ip1++) {
    // Two points on layer1: p1 and p3
    p1[0] = fX1[ip1];   p1[1] = fY1[ip1];    p1[2] = fZ1[ip1];   
    p3[0] = fX1[ip1+1]; p3[1] = fY1[ip1+1];  p3[2] = fZ1[ip1+1]; 
    hphi1 = fPhi1[ip1]; hphi3 = fPhi1[ip1+1];
    
    for (Int_t ip2=0; ip2<ficlu2; ip2++) {
      // Two points on layer 2: p2 and p4
      p2[0] = fX2[ip2];   p2[1] = fY2[ip2];     p2[2] = fZ2[ip2]; 
      p4[0] = fX2[ip2+1]; p4[1] = fY2[ip2+1];   p4[2] = fZ2[ip2+1]; 
      hphi2 = fPhi2[ip2]; hphi4 = fPhi2[ip2+1];

      // First line is formed by p1-p2 and second line by p3-p4
      // We find two points on each line which form the closest distance of the two lines
      // pa[0],pa[1],pa[2] are the points on line 1 and pb[0],pb[1],pb[2] are the points on line 2
      // Next: Consider x, y and z to be less than cuts[0], cuts[1] and cuts[2], respectively

      if(TMath::Abs(hphi1-hphi2)<phicut && TMath::Abs(hphi3-hphi4)<phicut){
	CalculateLine(p1, p2, p3, p4, pa, pb);

	// We consider the points where x, y and z points are less than xcut, ycut and zcut, respectively
	if (pa[0]>pos[0]-cuts[0] && pa[0]<pos[0]+cuts[0] && pa[1]>pos[1]-cuts[1] && pa[1]<pos[1]+cuts[1] && pa[2]>pos[2]-cuts[2] && pa[2]<pos[2]+cuts[2]){
	  distance = (TMath::Sqrt(pow((pa[0]-pb[0]),2) + pow((pa[1]-pb[1]),2) + pow((pa[2]-pb[2]),2)));
	  if(distance<distancecut){
	    ncontr++;
	    vertx  = vertx  + pa[0];       verty  = verty  + pa[1];       vertz  = vertz  + pa[2];
	    vertx2 = vertx2 + pa[0]*pa[0]; verty2 = verty2 + pa[1]*pa[1]; vertz2 = vertz2 + pa[2]*pa[2];
	    ncontr++;
	    vertx  = vertx  + pb[0];       verty  = verty  + pb[1];       vertz  = vertz  + pb[2];
	    vertx2 = vertx2 + pb[0]*pb[0]; verty2 = verty2 + pb[1]*pb[1]; vertz2 = vertz2 + pb[2]*pb[2];
	  }
	}
      }
      
      // Third line is formed by p1-p4 and fourth line by p3-p2
      // We find two points on each line which form the closest distance of the two lines
      // pa[0],pa[1],pa[2] are the points on line 3 and pb[0],pb[1],pb[2] are the points on line 4
      // Next: Consider x, y and z to be less than cuts[0], cuts[1] and cuts[2], respectively      
      if(TMath::Abs(hphi1-hphi4)<phicut && TMath::Abs(hphi3-hphi2)<phicut) {
	
	CalculateLine(p1, p4, p3, p2, pa, pb);
	if (pa[0]>pos[0]-cuts[0] && pa[0]<pos[0]+cuts[0] && pa[1]>pos[1]-cuts[1] && pa[1]<pos[1]+cuts[1] && pa[2]>pos[2]-cuts[2] && pa[2]<pos[2]+cuts[2]){
	  distance = (TMath::Sqrt(pow((pa[0]-pb[0]),2) + pow((pa[1]-pb[1]),2) + pow((pa[2]-pb[2]),2)));
	  if(distance<distancecut){
	    ncontr++;
	    vertx  = vertx  + pa[0];       verty  = verty  + pa[1];       vertz  = vertz  + pa[2];
	    vertx2 = vertx2 + pa[0]*pa[0]; verty2 = verty2 + pa[1]*pa[1]; vertz2 = vertz2 + pa[2]*pa[2];
	    ncontr++;
	    vertx  = vertx  + pb[0];       verty  = verty  + pb[1];       vertz  = vertz  + pb[2];
	    vertx2 = vertx2 + pb[0]*pb[0]; verty2 = verty2 + pb[1]*pb[1]; vertz2 = vertz2 + pb[2]*pb[2];
	  }
	}
      }
    }
  }

  if(ncontr>0){
    pos[0]  = vertx/ncontr;  pos[1] = verty/ncontr;   pos[2] = vertz/ncontr;
    vertx2 = vertx2/ncontr; verty2 = verty2/ncontr; vertz2 = vertz2/ncontr;
    sigpos[0] = TMath::Sqrt(vertx2 - pos[0]*pos[0]); 
    sigpos[1] = TMath::Sqrt(verty2 - pos[1]*pos[1]);
    sigpos[2] = TMath::Sqrt(vertz2 - pos[2]*pos[2]);
  }
  ncontr = ncontr/2;
}

void AliITSVertexer3DTapan::CalculatePhi(Float_t fx, Float_t fy, Float_t & phi)
{
  //calculates phi
  Float_t ybyx, phi1;
  const Float_t kradian = 180./3.141592654;
  
  if(fx==0.)
    {
      if(fy>0.) phi = 90.;
      if(fy<0.) phi = 270.;
    }
  if(fx != 0.)
    {
      ybyx = fy/fx;
      if(ybyx < 0) ybyx = - ybyx;
      phi1 = TMath::ATan(ybyx)*kradian;
      if(fx > 0 && fy > 0) phi = phi1;        // 1st Quadrant
      if(fx < 0 && fy > 0) phi = 180 - phi1;  // 2nd Quadrant
      if(fx < 0 && fy < 0) phi = 180 + phi1;  // 3rd Quadrant
      if(fx > 0 && fy < 0) phi = 360 - phi1;  // 4th Quadrant
      
    }
  phi = phi/kradian;
}

void AliITSVertexer3DTapan::CalculateLine(Double_t p1[4], Double_t p2[4], Double_t p3[4], Double_t p4[4], Double_t pa[3], Double_t pb[3]) const{
  //calculates line
  Double_t p13x, p13y, p13z;
  Double_t p21x, p21y, p21z;
  Double_t p43x, p43y, p43z;
  Double_t d1343, d4321, d1321, d4343, d2121;
  Double_t numer, denom;
  Double_t mua,   mub;
  mua = 0.; mub = 0.;
  
  p13x = p1[0] - p3[0];
  p13y = p1[1] - p3[1];
  p13z = p1[2] - p3[2];

  p21x = p2[0] - p1[0];
  p21y = p2[1] - p1[1];
  p21z = p2[2] - p1[2];

  p43x = p4[0] - p3[0];
  p43y = p4[1] - p3[1];
  p43z = p4[2] - p3[2];

  d1343 = p13x * p43x + p13y * p43y + p13z * p43z;
  d4321 = p43x * p21x + p43y * p21y + p43z * p21z;
  d1321 = p13x * p21x + p13y * p21y + p13z * p21z;
  d4343 = p43x * p43x + p43y * p43y + p43z * p43z;
  d2121 = p21x * p21x + p21y * p21y + p21z * p21z;

  denom = d2121 * d4343 - d4321 * d4321;
  numer = d1343 * d4321 - d1321 * d4343;

  if(denom>0) mua = numer / denom;
  if(d4343>0) mub = (d1343 + d4321 * (mua)) / d4343;

  pa[0] = p1[0] + mua * p21x;
  pa[1] = p1[1] + mua * p21y;
  pa[2] = p1[2] + mua * p21z;

  pb[0] = p3[0] + mub * p43x;
  pb[1] = p3[1] + mub * p43y;
  pb[2] = p3[2] + mub * p43z;
}

