/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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
/////////////////////////////////////////////////////////////////////////
// ALICE DISPLAY CLUSTERS CLASS                                        //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TFile.h>
#include <TPolyMarker3D.h>

#include "AliClusters.h"
#include "AliDisplay2.h"
#include "AliDisplayClusters.h"
#include "AliITS.h"
#include "AliITSLoader.h"
#include "AliITSclusterV2.h"
#include "AliITSgeom.h"
#include "AliModuleInfo.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliTPCLoader.h"
#include "AliTPCParam.h"
#include "AliTPCcluster.h"

ClassImp(AliDisplayClusters)

//_____________________________________________________________
AliDisplayClusters::AliDisplayClusters()
{
  //Default constructor
  fPoints = new TPolyMarker3D[gAliDisplay2->GetNbModules()];
  fName = new char*[gAliDisplay2->GetNbModules()];
  fNb=0;
  for(Int_t i=0;i<gAliDisplay2->GetNbModules();i++){
    fPoints[i].SetMarkerSize(0.2); 
    fPoints[i].SetMarkerColor(2); 
    fPoints[i].SetMarkerStyle(1);
  }
}

//_____________________________________________________________
AliDisplayClusters::~AliDisplayClusters()
{
  // Destructor
  delete [] fPoints;
}

//_____________________________________________________________
Int_t AliDisplayClusters::GetNbClusters()
{
  // Returns the number of clusters
  Int_t r=0;
  for(Int_t i=0;i<fNb;i++){
     if(gAliDisplay2->GetModuleInfo()->IsEnabled(fName[i])) r+=fPoints[i].GetN();
  }
  return r;
}

//_____________________________________________________________
void AliDisplayClusters::LoadClusters(const char *name,Int_t nevent)
{
  // Loads ITS and TPC clusters
  if(strstr(name,"ITS")) LoadITSClusters(nevent);
  if(strstr(name,"TPC")) LoadTPCClusters(nevent);
}

//_____________________________________________________________
void AliDisplayClusters::LoadITSClusters(Int_t nevent)
{
  // Loads ITS clusters
  fName[fNb]=new char[strlen("ITS")];
  strcpy(fName[fNb],"ITS");
  AliRunLoader *rl = AliRunLoader::Open("galice.root");
  if(!rl) {
    cerr<<"Can't open galice.root";
    return;
  }
  AliITSLoader *itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
  AliITS *its  = (AliITS*)gAlice->GetModule("ITS");
  
  rl->GetEvent(nevent);
  itsl->LoadRecPoints();
  TTree *cTree=itsl->TreeR();
  if(!cTree)
  {
    cerr<<"Error occured during ITS clusters load";
    return;
  }

  AliITSgeom *geom=its->GetITSgeom();
  Int_t count = 0;

  TClonesArray *clusters=new TClonesArray("AliITSclusterV2",10000);
  TBranch *branch=cTree->GetBranch("Clusters");
  branch->SetAddress(&clusters);
  Int_t nentr=(Int_t)cTree->GetEntries();
  for (Int_t i=0; i<nentr; i++) {
       if (!cTree->GetEvent(i)) continue;

   	Double_t rot[9];     
	geom->GetRotMatrix(i,rot);
       	Int_t lay,lad,det; 
	geom->GetModuleId(i,lay,lad,det);
  	Float_t tx,ty,tz;  
	geom->GetTrans(lay,lad,det,tx,ty,tz);     

      	Double_t r=-tx*rot[1]+ty*rot[0];          
	if (lay==1) r=-r;
      	Double_t phi=TMath::ATan2(rot[1],rot[0]); 
	if (lay==1) phi-=3.1415927;
      	Double_t cp=TMath::Cos(phi), sp=TMath::Sin(phi);

      	Int_t ncl=clusters->GetEntriesFast();
       	while (ncl--) {
           	AliITSclusterV2 *c=(AliITSclusterV2*)clusters->UncheckedAt(ncl);
           	Double_t g[3];
           	g[0]= r*cp + c->GetY()*sp; //
           	g[1]=-r*sp + c->GetY()*cp; //
	   	g[2]=c->GetZ();
		fPoints[fNb].SetPoint(count,g[0],g[1],g[2]);
		count++;
       	}
   }
  fNb++;
}

//_____________________________________________________________
void AliDisplayClusters::LoadTPCClusters(Int_t nevent)
{ 
  // Loads TPC clusters  
  fName[fNb]=new char[strlen("TPC")];
  strcpy(fName[fNb],"TPC");
  TFile *file = TFile::Open("galice.root");
  AliTPCParam *dig=(AliTPCParam *)file->Get("75x40_100x60_150x60");
  if (!dig) {cerr<<"TPC parameters have not been found !\n";}
  file->Close();

  AliRunLoader *rl = AliRunLoader::Open("galice.root");
  if(!rl) {
    cerr<<"Can't open galice.root";
    return;
  }
  AliTPCLoader *itsl = (AliTPCLoader*)rl->GetLoader("TPCLoader");
  if(!itsl){
    cerr<<"Can't find Loader";
    return;
  }
  
  rl->GetEvent(nevent);
  itsl->LoadRecPoints();
  TTree *cTree=itsl->TreeR();
  if(!cTree)
  {
    cerr<<"Error during TPC clusters load";
    return;
  }

  Int_t count = 0;
  Float_t noiseth = 10;

   AliClusters *clusters=new AliClusters(); 
   clusters->SetClass("AliTPCcluster");

   cTree->SetBranchAddress("Segment",&clusters);

   Int_t nrows=Int_t(cTree->GetEntries());
   for (Int_t n=0; n<nrows; n++) {
       cTree->GetEvent(n);
       Int_t sec,row;
       dig->AdjustSectorRow(clusters->GetID(),sec,row);
       TClonesArray &clrow=*clusters->GetArray();
       Int_t ncl=clrow.GetEntriesFast();
       while (ncl--) {
           AliTPCcluster *cl=(AliTPCcluster*)clrow[ncl];
           Double_t x=dig->GetPadRowRadii(sec,row), y=cl->GetY(), z=cl->GetZ();
	   if (cl->GetQ()<noiseth) continue;
           Float_t cs, sn, tmp;
           dig->AdjustCosSin(sec,cs,sn);
           tmp = x*cs-y*sn; y= x*sn+y*cs; x=tmp;
	   fPoints[fNb].SetPoint(count,x,y,z);
	   count++;
       }
       clrow.Clear();
      
   }
   delete cTree;
   delete dig;
   fNb++;
  
}

//_____________________________________________________________
void AliDisplayClusters::Draw()
{
  // Draws clusters
  for(Int_t i=0;i<fNb;i++){
    if(gAliDisplay2->GetModuleInfo()->IsEnabled(fName[i])) fPoints[i].Draw();
  }
}

