/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/*
$Log$
Revision 1.8  2000/11/12 17:17:03  pcrochet
BuildGeometry of AliMUON for trigger chambers delegated to AliMUONSegmentationTriggerX (same strategy as for tracking chambers)

Revision 1.7  2000/10/03 21:48:07  morsch
Adopt to const declaration of some of the methods in AliSegmentation.

Revision 1.6  2000/10/02 16:58:29  egangler
Cleaning of the code :
-> coding conventions
-> void Streamers
-> some useless includes removed or replaced by "class" statement

Revision 1.5  2000/07/03 11:54:57  morsch
AliMUONSegmentation and AliMUONHitMap have been replaced by AliSegmentation and AliHitMap in STEER
The methods GetPadIxy and GetPadXxy of AliMUONSegmentation have changed name to GetPadI and GetPadC.

Revision 1.4  2000/06/29 12:34:09  morsch
AliMUONSegmentation class has been made independent of AliMUONChamber. This makes
it usable with any other geometry class. The link to the object to which it belongs is
established via an index. This assumes that there exists a global geometry manager
from which the pointer to the parent object can be obtained (in our case gAlice).

Revision 1.3  2000/06/26 10:01:26  pcrochet
global variables removed

Revision 1.2  2000/06/15 07:58:48  morsch
Code from MUON-dev joined

Revision 1.1.2.1  2000/06/09 21:51:03  morsch
Code from AliMUONSegResTrigger.cxx

*/



/*
Old Log:
Revision 1.1.2.4  2000/04/26 12:33:25  morsch
Minor changes in some methods (CP)

Revision 1.1.2.3  2000/03/20 18:14:16  morsch
Missing sector added.

Revision 1.1.2.2  2000/02/20 07:50:49  morsch
Bugs in Dpx, Dpy and ISector methods corrected (P.C.)

Revision 1.1.2.1  2000/02/17 14:33:49  morsch
Draft version from P. Crochet

*/

#include <TMath.h>
#include <TBRIK.h>
#include <TNode.h>
#include <TGeometry.h>

#include "AliMUON.h"
#include "AliMUONSegmentationTriggerX.h"
#include "AliMUONTriggerConstants.h"
#include "TMath.h"
#include "TRandom.h"
#include "TArc.h"
#include "AliMUONChamber.h"
#include "AliRun.h"  // gAlice
#include <iostream.h> 
ClassImp(AliMUONSegmentationTriggerX)

//------------------------------------------------------------------
void AliMUONSegmentationTriggerX::Init(Int_t chamber)
{
// intialize X segmentation 
  cout << "Initialize Trigger Chamber Geometry X " << "\n";    
  AliMUONSegmentationTrigger::Init(chamber);

// calculate x & y position of X strips
  Int_t nModule=AliMUONTriggerConstants::Nmodule();
  for (Int_t imodule=0; imodule<nModule; imodule++) {
    Float_t width=StripSizeX(AliMUONTriggerConstants::ModuleId(imodule));     
    Int_t nStrip=AliMUONTriggerConstants::NstripX(imodule);
    for (Int_t istrip=0; istrip<nStrip; istrip++){    
      fXofxsmin[imodule][istrip] = AliMUONTriggerConstants::XcMin(imodule)*fZscale;
      fXofxsmax[imodule][istrip] = AliMUONTriggerConstants::XcMax(imodule)*fZscale;
      
      fYofxsmin[imodule][istrip] = (fYcmin[imodule]+width*(istrip))*fZscale;
      fYofxsmax[imodule][istrip] = (fYcmin[imodule]+width*(istrip+1))*fZscale;
    }
  }
}

//------------------------------------------------------------------
void AliMUONSegmentationTriggerX::GetPadI(Float_t x,Float_t y,Int_t &ix,Int_t &iy) 
{
//  Returns pad coordinates (ix,iy) for given real coordinates (x,y)
//  x,y = real coordinates; ix = module number , iy = strip number
  ix = 0;    
  iy = 0;
  Int_t nModule=AliMUONTriggerConstants::Nmodule();
  for (Int_t imodule=0; imodule<nModule; imodule++) {
      Int_t nStrip=AliMUONTriggerConstants::NstripX(imodule);
      for (Int_t istrip=0; istrip<nStrip; istrip++){
	  if (x>fXofxsmin[imodule][istrip]&&x<fXofxsmax[imodule][istrip]&&
	      y>fYofxsmin[imodule][istrip]&&y<fYofxsmax[imodule][istrip]){
	      ix = AliMUONTriggerConstants::ModuleId(imodule);
	      iy = istrip;
	  }
      }
  }
}

//------------------------------------------------------------------
void AliMUONSegmentationTriggerX::GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y) 
{
//  Returns real coordinates (x,y) for given pad coordinates (ix,iy)
//  ix = module number , iy = strip number;  x,y = center of strip
  x = 0.;    
  y = 0.; 
  Int_t nModule=AliMUONTriggerConstants::Nmodule();

  for (Int_t imodule=0; imodule<nModule; imodule++) {
    if (AliMUONTriggerConstants::ModuleId(imodule)==ix){
      x=fXofxsmin[imodule][iy]+(fXofxsmax[imodule][iy]-fXofxsmin[imodule][iy])/2.;
      y=fYofxsmin[imodule][iy]+(fYofxsmax[imodule][iy]-fYofxsmin[imodule][iy])/2.;
    }
  }
}

//------------------------------------------------------------------
void AliMUONSegmentationTriggerX::SetPadSize(Float_t p1, Float_t p2)
{
//  Sets the padsize 
//  
  fDpx=p1;
  fDpy=p2;
}

//------------------------------------------------------------------
void AliMUONSegmentationTriggerX::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]){
// Returns list of 10 next neighbours for given X strip (ix, iy)  
// neighbour number 4 in the list -                     
// neighbour number 3 in the list  |                    
// neighbour number 2 in the list  |_ Upper part             
// neighbour number 1 in the list  |            
// neighbour number 0 in the list -           
//      X strip (ix, iy) 
// neighbour number 5 in the list -       
// neighbour number 6 in the list  | _ Lower part
// neighbour number 7 in the list  |
// neighbour number 8 in the list  | 
// neighbour number 9 in the list -
  
  Int_t absiX = TMath::Abs(iX); 
  Int_t numModule = ModuleNumber(absiX);   // module number Id.
  Int_t nStrip = AliMUONTriggerConstants::NstripX(numModule); //numb of strips
  Int_t iCandidateUp, iCandidateDo;
  Int_t j;

  *Nlist = 10;
  for (Int_t i=0; i<10; i++) Xlist[i]=Ylist[i]=0;

  if (iY < nStrip) {

    for (Int_t i=0; i<5; i++) {
      j = i + 5;
      iCandidateUp = iY + (i + 1);
      iCandidateDo = iY - (i + 1);
      if (iCandidateUp < nStrip) { 
	Xlist[i] = iX;
	Ylist[i] = iCandidateUp;  
      }
      if (iCandidateDo >= 0) { 
	Xlist[j] = iX;
	Ylist[j] = iCandidateDo;  
      }
    }
    
  } // iY < nStrip
}

//------------------------------------------------------------------   
void AliMUONSegmentationTriggerX::SetPad(Int_t ix, Int_t iy)
{
  // Sets virtual pad coordinates, needed for evaluating pad response 
  // outside the tracking program
  GetPadC(ix,iy,fX,fY);
  GetPadI(fX,fY,fIx,fIy);
  fSector=Sector(ix,iy);
}

//------------------------------------------------------------------   
Int_t AliMUONSegmentationTriggerX::ISector() 
{ return fSector;}

//------------------------------------------------------------------   
Int_t AliMUONSegmentationTriggerX::Ix()
{ return fIx;}

//------------------------------------------------------------------   

Int_t AliMUONSegmentationTriggerX::Iy()
{ return fIy;}

//------------------------------------------------------------------
Float_t AliMUONSegmentationTriggerX::Dpx(Int_t isec) const
{ 
// returns x size of x strips for sector isec
    
  if (isec==1) {
    return 17.0*fZscale;
  } else if (isec==2) {
    return 34.0*fZscale;
  } else if (isec==3) {
    return 34.0*fZscale;
  } else if (isec==4) {
    return 34.0*fZscale;
  } else if (isec==5) {
    return 34.0*fZscale;
  } else if (isec==6) {
    return 68.0*fZscale;
  } else {
    return 0.;
  }
}

//------------------------------------------------------------------
Float_t AliMUONSegmentationTriggerX::Dpy(Int_t isec) const
{ 
// returns y size of x strips for sector isec

  if (isec==1) {
    return 1.0625*fZscale;
  } else if (isec==2) {
    return 1.0625*fZscale;
  } else if (isec==3) {
    return 1.0625*fZscale;
  } else if (isec==4) {
    return 2.125*fZscale;
  } else if (isec==5) {
    return 4.25*fZscale;
  } else if (isec==6) {
    return 4.25*fZscale;
  } else {
    return 0.;	
  }   
}

//------------------------------------------------------------------   
void AliMUONSegmentationTriggerX::SetHit(Float_t xhit, Float_t yhit)
{ 
// set hit during disIntegration
AliMUONSegmentationTrigger::SetHit(xhit,yhit);
}

//------------------------------------------------------------------   
Int_t AliMUONSegmentationTriggerX::Sector(Int_t ix, Int_t iy) 
{
// Returns sector number for given module
// 
    
  Int_t absix=TMath::Abs(ix);
  Int_t iwidth=Int_t(StripSizeX(absix));

  if (absix==52) {
    return 1;
  } else if (absix==41||absix==61) {
    return 2;
  } else if (iwidth==1) {
    return 3;
  } else if (iwidth==2) {
    return 4;
  } else if ((absix>=11&&absix<17)||(absix>=91&&absix<97)) {
    return 5;
  } else if (iwidth==4) {
    return 6;
  } else {
    return 0;
  }
}

//------------------------------------------------------------------   
void AliMUONSegmentationTriggerX::
IntegrationLimits(Float_t& x1, Float_t& x2, Float_t& x3, Float_t& x4) 
{ 
// returns quantities needed to evaluate neighbour strip response

  Int_t ix,iy;
  Float_t xstrip,ystrip;
  GetPadI(fXhit,fYhit,ix,iy);  
  GetPadC(ix,iy,xstrip,ystrip);  
  x1=fYhit;        // hit y position
  x2=ystrip;       // y coordinate of the main strip
  x3=fY;           // current strip real y coordinate  
  //  width=StripSizeX(ix);   // width of the main strip 

  // find the position of the 2 borders of the current strip
  Float_t ymin = fYofxsmin[ModuleNumber(fIx)][fIy];
  Float_t ymax = fYofxsmax[ModuleNumber(fIx)][fIy];
  
  // dist. between the hit and the closest border of the current strip
  x4 = (TMath::Abs(ymax-x1) > TMath::Abs(ymin-x1)) ? 
    TMath::Abs(ymin-x1):TMath::Abs(ymax-x1);    

}

//------------------------------------------------------------------   
void AliMUONSegmentationTriggerX::Draw(const char* opt) const
{
  
  if (!strcmp(opt,"eventdisplay")) { 
    TNode *node, *nodeS;
    char nameChamber[10], nameNode[10];
    char nameSense1[10], nameSense2[10], nameSense3[10], nameSense4[10];
    
    TNode* top=gAlice->GetGeometry()->GetNode("alice"); 
    sprintf(nameChamber,"C_MUON%d",fId+1);
    new TBRIK(nameChamber,"Mother","void",340.,340.,0.25);
    top->cd();
    sprintf(nameNode,"MUON%d",100+fId+1);
    node = new TNode(nameNode,"Chambernode",nameChamber,0,0,fChamber->Z(),"");
    node->SetLineColor(kBlack);    
    AliMUON *pMUON  = (AliMUON *) gAlice->GetModule("MUON");
    (pMUON->Nodes())->Add(node);
    
    sprintf(nameSense1,"S1_MUON%d",fId+1);
    sprintf(nameSense2,"S2_MUON%d",fId+1);
    sprintf(nameSense3,"S3_MUON%d",fId+1);
    sprintf(nameSense4,"S4_MUON%d",fId+1);
    
    for (Int_t imodule=0; imodule<AliMUONTriggerConstants::Nmodule(); imodule++) {    
      Int_t idModule=AliMUONTriggerConstants::ModuleId(imodule);
      
      if (TMath::Abs(idModule)!=51) {	 
	
	Int_t nStripX=AliMUONTriggerConstants::NstripX(imodule);
	Float_t xmin=fXofxsmin[imodule][0];
	Float_t xmax=fXofxsmax[imodule][nStripX-1];
	Float_t ymin=fYofxsmin[imodule][0];
	Float_t ymax=fYofxsmax[imodule][nStripX-1];
	Float_t xpos=xmin+(xmax-xmin)/2.;
	Float_t ypos=ymin+(ymax-ymin)/2.;
	Float_t halfx=(xmax-xmin)/2.;
	Float_t halfy=(ymax-ymin)/2.;
	
	if (idModule==11) 
	  new TBRIK(nameSense1,"Module","void",halfx,halfy,0.25);   
	if (idModule==17) 
	  new TBRIK(nameSense2,"Module","void",halfx,halfy,0.25);   
	if (idModule==41) 
	  new TBRIK(nameSense3,"Module","void",halfx,halfy,0.25);   
	if (idModule==52) 
	  new TBRIK(nameSense4,"Module","void",halfx,halfy,0.25); 
	node->cd();
	sprintf(nameNode,"S_MUON%d",1000*fId+1+imodule);
	
	if (TMath::Abs(idModule)==41||TMath::Abs(idModule)==61) {
	  nodeS = new TNode(nameNode,"Module",nameSense3,xpos,ypos,0,"");
	} else if (TMath::Abs(idModule)==52) {
	  nodeS = new TNode(nameNode,"Module",nameSense4,xpos,ypos,0,"");
	} else if (TMath::Abs((idModule-Int_t(idModule/10)*10.))!=7) {
	  nodeS = new TNode(nameNode,"Module",nameSense1,xpos,ypos,0,"");
	} else {
	  //	} else if (TMath::Abs((idModule-Int_t(idModule/10)*10.))==7) {
	  nodeS = new TNode(nameNode,"Module",nameSense2,xpos,ypos,0,"");
	}
	nodeS->SetLineColor(kBlue);
	node->cd();
      }
    }
  }
}


