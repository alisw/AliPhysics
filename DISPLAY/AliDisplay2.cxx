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

#define do_mc

//standard modules
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

//ROOT

#include <TGButton.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <TGIcon.h>
#include <TApplication.h>
#include <TGFileDialog.h>
#include <TGMenu.h>
#include <TGTab.h>
#include <TView.h>
#include <TGLayout.h>
#include <TGeometry.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TBranch.h>
#include <TVirtualX.h>
#include <TPolyMarker3D.h>
#include <TGNumberEntry.h>
#include <TSystem.h>
#include <TRootHelpDialog.h>
#include <TParticle.h>
#include <TGShutter.h>
#include <TGWindow.h>
#include <TEnv.h>
#include <TGLabel.h>
#include <TCanvas.h>
#include <TClonesArray.h>

//AliRoot Module
#include "AliModule.h"
#include "AliDetector.h"
#include "AliHeader.h"
#include "AliPoints.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliITSclusterV2.h"
#include "AliITSgeom.h"
#include "AliITSLoader.h"
#include "AliITS.h"
#include "AliTPCcluster.h"
#include "AliTPCLoader.h"
#include "AliTPC.h"
#include "AliClusters.h"
#include "AliTPCParam.h"
#include "AliMC.h"

#ifdef ALI_HLT
//HLT Module 
#include "AliL3Logging.h"
#include "AliL3Display.h"
#include "AliL3Transform.h"
#include "AliL3Track.h"
#include "AliL3TrackArray.h"
#include "AliL3SpacePointData.h"
#include "AliL3MemHandler.h"
#endif

#include "AliDisplay2.h"

 AliDisplay2 *gAliDisplay2;

//extern filetypes;
const char *filetypes[] = {"ROOT files","*.root","All files","*",0,0};
const char *imgtypes[] = {"GIF files","*.gif",0,0};

ClassImp(AliModuleInfo);

//_____________________________________________________________
AliModuleInfo::AliModuleInfo(int n)
{
  fName = new char*[n];
  fId = new Int_t[n];
  fEnabled = new Bool_t[n];
  fNb = 0;
}

//_____________________________________________________________
AliModuleInfo::~AliModuleInfo(){
  delete [] fName;
  delete [] fId;
  delete [] fEnabled;
}

//_____________________________________________________________
void AliModuleInfo::Add(const char *name,Int_t i)
{
  fName[fNb]=new char[strlen(name)];
  strcpy(fName[fNb],name);
  fId[fNb]=i;
  fEnabled[fNb]=kTRUE;
  fNb++;
}

//_____________________________________________________________
void AliModuleInfo::SetId(char *name,Int_t id)
{
  Int_t i=0;
  while((strcmp(name,fName[i])!=0)&&(i!=fNb)) i++;
  if(strcmp(name,fName[i])==0) fId[i]=id;
}

//_____________________________________________________________
char* AliModuleInfo::Name(Int_t id)
{
  Int_t i=0;
  while((fId[i]!=id)&&(i!=fNb)) i++;
  if(fId[i]==id) return fName[i];
  return 0;
}
  
//_____________________________________________________________
Int_t AliModuleInfo::Id(char *name)
{
  Int_t i=0;
  while((strcmp(name,fName[i])!=0)&&(i!=fNb)) i++;
  if(strcmp(name,fName[i])==0) return fId[i];
  return -1;
}

//_____________________________________________________________
Bool_t AliModuleInfo::IsEnabled(Int_t id)
{
  //return the current status of the detector
  Int_t i=0;
  while((fId[i]!=id)&&(i!=fNb)) i++;
  if(fId[i]==id) return fEnabled[i];
  return kFALSE;
}

//_____________________________________________________________
void AliModuleInfo::Disable(Int_t id)
{
  //Disable the detector 
  Int_t i=0;
  while((fId[i]!=id)&&(i!=fNb)) i++;
  if(fId[i]==id) fEnabled[i]=kFALSE;
}

//_____________________________________________________________
void AliModuleInfo::Enable(Int_t id)
{
  //Enable the detector 
  Int_t i=0;
  while((fId[i]!=id)&&(i!=fNb)) i++;
  if(fId[i]==id) fEnabled[i]=kTRUE;
}

//_____________________________________________________________
void AliModuleInfo::Print()
{
  printf("\n***Modules***");
  printf("\nName\tId\tEnabled"); 
  for(Int_t i=0;i<fNb;i++){
    printf("\n%s",fName[i]);
    printf("\t%d",fId[i]);
    printf("\t%d",fEnabled[i]);
  }
}


ClassImp(AliDisplayClusters);

//_____________________________________________________________
AliDisplayClusters::AliDisplayClusters()
{
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
  delete [] fPoints;
}

//_____________________________________________________________
Int_t AliDisplayClusters::GetNbClusters()
{
  Int_t r=0;
  for(Int_t i=0;i<fNb;i++){
     if(gAliDisplay2->GetModuleInfo()->IsEnabled(fName[i])) r+=fPoints[i].GetN();
  }
  return r;
}

//_____________________________________________________________
void AliDisplayClusters::LoadClusters(char *name,Int_t nevent)
{
  if(strstr(name,"ITS")) LoadITSClusters(nevent);
  if(strstr(name,"TPC")) LoadTPCClusters(nevent);
}

//_____________________________________________________________
void AliDisplayClusters::LoadITSClusters(Int_t nevent)
{
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
  for(Int_t i=0;i<fNb;i++){
    if(gAliDisplay2->GetModuleInfo()->IsEnabled(fName[i])) fPoints[i].Draw();
  }
}


ClassImp(AliDisplayHLT);

//_____________________________________________________________
AliDisplayHLT::AliDisplayHLT()
{
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
AliDisplayHLT::~AliDisplayHLT()
{
  delete [] fPoints;
}

//_____________________________________________________________
void AliDisplayHLT::LoadHLT(char *name,Int_t nevent)
{
  if(strstr(name,"TPC")) LoadHLTTPC(nevent);
}

//_____________________________________________________________
void AliDisplayHLT::LoadHLTTPC(Int_t nevent)
{
#ifdef ALI_HLT
  //load TPC Clusters from the raw data
  //raw data must be in the directorie $ALICE_ROOT/raw
  //First we read the data from files
  //do_mc MUST BE DEFINED AND USED FOR RAW DATA GENERATION

  fName[fNb]=new char[strlen("TPC")];
  strcpy(fName[fNb],"TPC");
  Char_t fname[256];
  AliL3MemHandler *clusterfile[36][6];
  AliL3SpacePointData *fClusters[36][6];
  UInt_t fNcl[36][6];
  memset(fClusters,0,36*6*sizeof(AliL3SpacePointData*));
  //  strcpy(path,gSystem->Getenv("ALICE_ROOT"));
  //strcat(path,"/raw");
  //printf("\nRaw data path %s",path);
  char path[128];
  strcpy(path,gAliDisplay2->GetRawDataPath());
  for(Int_t s=0; s<36; s++)
    {
      for(Int_t p=0; p<AliL3Transform::GetNPatches(); p++)
	{
	  Int_t patch;
	  patch=-1;
	  clusterfile[s][p] = new AliL3MemHandler();
	  if(nevent<0)
	    sprintf(fname,"%s/points_%d_%d.raw",path,s,patch);
	  else
	    sprintf(fname,"%s/points_%d_%d_%d.raw",path,nevent,s,patch);
	  if(!clusterfile[s][p]->SetBinaryInput(fname))
	    {
	      LOG(AliL3Log::kError,"AliL3Evaluation::Setup","File Open")
		<<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
	      delete clusterfile[s][p];
              clusterfile[s][p] = 0; 
	      continue;
	    }
	  fClusters[s][p] = (AliL3SpacePointData*)clusterfile[s][p]->Allocate();
	  clusterfile[s][p]->Binary2Memory(fNcl[s][p],fClusters[s][p]);
	  clusterfile[s][p]->CloseBinaryInput();
	  break;
	}
    }

  //Second step: we assign the clusters to the fPoints array
  Int_t nbc=0;
  for(Int_t s=0; s<36; s++){
      for(Int_t p=0;p<6;p++){
	  AliL3SpacePointData *points = fClusters[s][p];
	  if(!points) continue;
	  Float_t xyz[3];
	  for(UInt_t i=0; i<fNcl[s][p]; i++){
	      xyz[0] = points[i].fX;
	      xyz[1] = points[i].fY;
	      xyz[2] = points[i].fZ;
	      fPoints[fNb].SetPoint(i+nbc,xyz[0],xyz[1],xyz[2]);
 	    }
	    nbc += fNcl[s][p];
	}
    }
  fNb++;
#else
  printf("This is event %d\n",nevent);
#endif
}

//_____________________________________________________________
void AliDisplayHLT::Draw()
{
  for(Int_t i=0;i<fNb;i++){
   if(gAliDisplay2->GetModuleInfo()->IsEnabled(fName[i])) fPoints[i].Draw();
  }
}


ClassImp(AliSliderFrame);

//_____________________________________________________________
AliSliderFrame::AliSliderFrame(const TGWindow *p, UInt_t w, UInt_t h)
{
	fMainFrame = new TGCompositeFrame(p, w, h, kVerticalFrame | kRaisedFrame);
	fLayout = new TGLayoutHints( kLHintsBottom | kLHintsRight | kLHintsExpandX,5,5,2,2);

	//Momentum Slider
	fMomentumFrame = new TGCompositeFrame(fMainFrame,0,0, kHorizontalFrame);
	fMomentumLayout = new TGLayoutHints( kLHintsLeft | kLHintsTop | kLHintsExpandX ,5,5,2,2);
	fMomentumSlider = new TGDoubleHSlider(fMomentumFrame, 400, kDoubleScaleBoth, kIdsMOMENTUM);
	fMomentumSlider->Connect("PositionChanged()","AliSliderFrame",this,"DoSlider()");
	fMomentumSlider->Connect("PositionChanged()","AliSliderFrame",this,"DoPositionChanged(Int_t)");
	fMomentumSlider->Connect("Released()","AliSliderFrame",this,"DoReleased()");
	SetMomentumRange(0,2);
	fMomentumSlider->SetPosition(0,2);
	fMomentumLabel = new TGLabel(fMomentumFrame,"Momentum");
	
	fMomentumMinValue = new TGNumberEntry(fMomentumFrame,GetMomentumMin(),7,kIdsMomentumMIN);
	fMomentumMinValue->Connect("ValueChanged(Long_t)","AliSliderFrame",this,"DoField(Long_t)");
	fMomentumMinValue->GetNumberEntry()->Connect("ReturnPressed()","AliSliderFrame",this,"DoField(Long_t)");
	fMomentumMinValue->SetButtonToNum(kFALSE);
	
	fMomentumMaxValue = new TGNumberEntry(fMomentumFrame,GetMomentumMax(),7,kIdsMomentumMAX);
	fMomentumMaxValue->Connect("ValueChanged(Long_t)","AliSliderFrame",this,"DoField(Long_t)");
	fMomentumMaxValue->GetNumberEntry()->Connect("ReturnPressed()","AliSliderFrame",this,"DoField(Long_t)");
	fMomentumMaxValue->SetButtonToNum(kFALSE);
	
	fMomentumFrame->AddFrame(fMomentumLabel,new TGLayoutHints( kLHintsLeft | kLHintsCenterY,5,5,2,2));
	fMomentumFrame->AddFrame(fMomentumMinValue,new TGLayoutHints( kLHintsLeft ,5,5,2,2));
	fMomentumFrame->AddFrame(fMomentumSlider,fMomentumLayout);
	fMomentumFrame->AddFrame(fMomentumMaxValue,new TGLayoutHints( kLHintsRight,5,5,2,2));
	fMomentumFrame->Resize(700,100);


	//Rapidity Slider
	fRapidityFrame = new TGCompositeFrame(fMainFrame,0,0, kHorizontalFrame);
	fRapidityLayout = new TGLayoutHints( kLHintsLeft | kLHintsTop | kLHintsExpandX,5,5,2,2);
	fRapiditySlider = new TGDoubleHSlider(fRapidityFrame, 400, kDoubleScaleBoth, kIdsRAPIDITY);
	fRapiditySlider->Connect("PositionChanged()","AliSliderFrame",this,"DoSlider()");
	fRapiditySlider->Connect("PositionChanged()","AliSliderFrame",this,"DoPositionChanged(Int_t)");
	fRapiditySlider->Connect("Released()","AliSliderFrame",this,"DoReleased()");
	SetRapidityRange(-1.5,1.5);
	fRapiditySlider->SetPosition(-1.5,1.5);
	fRapidityLabel = new TGLabel(fRapidityFrame,"Rapidity    ");
	fRapidityMinValue = new TGNumberEntry(fRapidityFrame,GetRapidityMin(),7,kIdsRapidityMIN);
	fRapidityMinValue->Connect("ValueChanged(Long_t)","AliSliderFrame",this,"DoField(Long_t)");
	fRapidityMinValue->SetButtonToNum(kFALSE);
	fRapidityMinValue->GetNumberEntry()->Connect("ReturnPressed()","AliSliderFrame",this,"DoField(Long_t)");
	
	fRapidityMaxValue = new TGNumberEntry(fRapidityFrame,GetRapidityMax(),7,kIdsRapidityMAX);
	fRapidityMaxValue->Connect("ValueChanged(Long_t)","AliSliderFrame",this,"DoField(Long_t)");
	fRapidityMaxValue->SetButtonToNum(kFALSE);
	fRapidityMaxValue->GetNumberEntry()->Connect("ReturnPressed()","AliSliderFrame",this,"DoField(Long_t)");
	
	fRapidityFrame->AddFrame(fRapidityLabel,new TGLayoutHints( kLHintsLeft | kLHintsCenterY,5,5,2,2));
	fRapidityFrame->AddFrame(fRapidityMinValue,new TGLayoutHints( kLHintsLeft ,5,5,2,2));
	fRapidityFrame->AddFrame(fRapiditySlider,fRapidityLayout);
	fRapidityFrame->AddFrame(fRapidityMaxValue,new TGLayoutHints( kLHintsRight ,5,5,2,2));
	fRapidityFrame->Resize(700,100);
	
	fMainFrame->AddFrame(fMomentumFrame,fLayout);
	fMainFrame->AddFrame(fRapidityFrame,fLayout);
	fMainFrame->MapSubwindows();
	fMainFrame->MapWindow();
	fMainFrame->Resize(700,100);
	LoadFromRC();
}

//_____________________________________________________________
AliSliderFrame::~AliSliderFrame(void)
{
	delete fLayout;
	delete fRapidityLayout;
	delete fMomentumLayout;
	delete fMomentumMinValue;
	delete fMomentumMaxValue;
	delete fMainFrame;
	delete fRapidityFrame;
	delete fMomentumFrame;
	delete fMomentumSlider;
	delete fRapiditySlider;
	delete fMomentumLabel;
	delete fRapidityLabel;
	delete fRapidityMinValue;
	delete fRapidityMaxValue;
}

//_____________________________________________________________
void AliSliderFrame::CloseWindow(void)
{
	delete this;
}

//_____________________________________________________________
void AliSliderFrame::DoSlider(Int_t /*pos*/)
{
	TGFrame *frame = (TGFrame *) gTQSender;
	TGDoubleSlider * ds = (TGDoubleSlider *) frame;
	int id = ds->WidgetId();
	char min[8];
	char max[8];
	switch(id){
	case kIdsMOMENTUM:{
		//sprintf(buf,"momentum min:%f max:%f",GetMomentumMin(),GetMomentumMax());
		//printf("\n%s",buf);
		sprintf(min,"%.4f",GetMomentumMin());
		sprintf(max,"%.4f",GetMomentumMax());
		fMomentumMaxValue->SetText(max);
		fMomentumMinValue->SetText(min);
							}
		break;
	case kIdsRAPIDITY:{
		sprintf(min,"%.4f",GetRapidityMin());
		sprintf(max,"%.4f",GetRapidityMax());
		fRapidityMaxValue->SetText(max);
		fRapidityMinValue->SetText(min);
							}
		break;
	default:break;
	}
}

//_____________________________________________________________
void AliSliderFrame::DoReleased(Int_t /*pos*/)
{
	TGFrame *frame = (TGFrame *) gTQSender;
	TGDoubleSlider * ds = (TGDoubleSlider *) frame;
	int id = ds->WidgetId();
	switch(id){
	case kIdsRAPIDITY:
	case kIdsMOMENTUM:
		gAliDisplay2->Update(kmCUTS);
		break;
	default: break;
	}
}

//_____________________________________________________________
void AliSliderFrame::DoPositionChanged(Int_t /*pos*/)
{
        if(!gAliDisplay2->GetSliderUpdate()) return;
	TGFrame *frame = (TGFrame *) gTQSender;
	TGDoubleSlider * ds = (TGDoubleSlider *) frame;
	int id = ds->WidgetId();
	switch(id){
	case kIdsRapidityMIN:
	case kIdsRapidityMAX:
	case kIdsMomentumMIN:
	case kIdsMomentumMAX:{
		gAliDisplay2->Update(kmCUTS);
	}
		break;
	default: break;
	}
}

//_____________________________________________________________
void AliSliderFrame::DoField(Long_t pos)
{
	TGNumberEntry *ne = (TGNumberEntry *) gTQSender;
	int id = ne->WidgetId();
	char max[8],min[8];
	Int_t sign=0;
	Float_t step=gAliDisplay2->GetSliderStep();
	if((pos/10000)==0){//Up button pressed
		sign = 1;
	}
	else sign = -1;

	switch(id){
	case kIdsMomentumMIN:{
		fMomentumMinValue->SetNumber(fMomentumMinValue->GetNumber()+step*sign);
		sprintf(min,"%.4f",fMomentumMinValue->GetNumber());
		fMomentumMinValue->SetText(min);
								 }
		break;
	case kIdsMomentumMAX:{
		fMomentumMaxValue->SetNumber(fMomentumMaxValue->GetNumber()+step*sign);
		sprintf(max,"%.4f",fMomentumMaxValue->GetNumber());
		fMomentumMaxValue->SetText(max);
								 }
		break;
	case kIdsRapidityMIN:{
		fRapidityMinValue->SetNumber(fRapidityMinValue->GetNumber()+step*sign);
		sprintf(min,"%.4f",fRapidityMinValue->GetNumber());
		fRapidityMinValue->SetText(min);
								 }
		break;
	case kIdsRapidityMAX:{
		fRapidityMaxValue->SetNumber(fRapidityMaxValue->GetNumber()+step*sign);
		sprintf(max,"%.4f",fRapidityMaxValue->GetNumber());
		fRapidityMaxValue->SetText(max);
								 }
		break;
	default:break;
	}

	fMomentumSlider->SetPosition(fMomentumMinValue->GetNumber(),fMomentumMaxValue->GetNumber());
	fRapiditySlider->SetPosition(fRapidityMinValue->GetNumber(),fRapidityMaxValue->GetNumber());
//	gAliDisplay2->Draw();
	gAliDisplay2->Update(kmCUTS);
}

//_____________________________________________________________
void AliSliderFrame::SaveToRC()
{
	TEnv *rc=new TEnv(".alidisplayrc");
	rc->SetValue("AliDisplay.MomentumMin",GetMomentumMin());
	rc->SetValue("AliDisplay.MomentumMax",GetMomentumMax());
	rc->SetValue("AliDisplay.RapidityMin",GetRapidityMin());
	rc->SetValue("AliDisplay.RapidityMax",GetRapidityMax());
	rc->SaveLevel(kEnvLocal);
	rc->Save();
}

//_____________________________________________________________	
void AliSliderFrame::LoadFromRC()
{
  TEnv *rc=new TEnv(".alidisplayrc");
  Float_t a,b;
  a=rc->GetValue("AliDisplay.MomentumMin",0);
  b=rc->GetValue("AliDisplay.MomentumMax",2);
  fMomentumSlider->SetPosition(a,b);
  a=rc->GetValue("AliDisplay.RapidityMin",-1.5);
  b=rc->GetValue("AliDisplay.RapidityMax",1.5);
  fRapiditySlider->SetPosition(a,b);
}


ClassImp(AliDetectorFrame);

int AliDetectorFrame::fgBaseId = 1000;

//_____________________________________________________________
AliDetectorFrame::AliDetectorFrame(const TGWindow *p, Int_t w, Int_t h,UInt_t bgc)
{
	fMainFrame = new TGCompositeFrame(p,w,h,kVerticalFrame,bgc);
	TGLayoutHints *layout = new TGLayoutHints(kLHintsTop | kLHintsLeft,2,2,2,2);
	TGLayoutHints *layout2 = new TGLayoutHints(kLHintsTop | kLHintsRight,2,2,2,2);
	TGLayoutHints *layout3 = new TGLayoutHints(kLHintsTop | kLHintsExpandX,2,2,2,2);
	fCheckButton = new TGCheckButton*[gAliDisplay2->GetNbModules()];
	fCheckedButton = new Bool_t[gAliDisplay2->GetNbModules()];
	fCheckButtonId = new Int_t[gAliDisplay2->GetNbModules()];
	TGCompositeFrame *dframe;
	TGButton 		 *button;
	char			 text[32];
	AliModule 	*mod;
	fCheckedMode = kFALSE;
	for(Int_t i=0;i<gAliDisplay2->GetNbModules();i++){
		mod = dynamic_cast<AliModule*> (gAliDisplay2->GetModules()->At(i));
		if(!mod) continue;
		dframe = new TGCompositeFrame(fMainFrame,150,20,kHorizontalFrame);
		fCheckButton[i] = new TGCheckButton(dframe,mod->GetName(),fgBaseId);
		fCheckButtonId[i]=fgBaseId;
		fCheckButton[i]->Connect("Clicked()","AliDetectorFrame",this,"DoCheckButton(Int_t)");		
		fCheckedButton[i]=kTRUE;
		dframe->AddFrame(fCheckButton[i],layout);
		fCheckButton[i]->SetState(kButtonDown);
		sprintf(text,"Specific %s view",mod->GetName());
		button = new TGTextButton(dframe,"Display",fgBaseId);
		button->SetToolTipText(text);
		button->Connect("Clicked()","AliDetectorFrame",this,"DoSpecific()");
		dframe->AddFrame(button,layout2);
		fMainFrame->AddFrame(dframe,layout3);
		gAliDisplay2->GetModuleInfo()->SetId((char*)mod->GetName(),fgBaseId);
		fgBaseId++;
	}	
	gAliDisplay2->GetModuleInfo()->Print();
	fButtonFrame = new TGCompositeFrame(fMainFrame,w,100,kHorizontalFrame,bgc);
	fButtonAll = new TGTextButton(fButtonFrame,"All",kIdbSelectALL);
	fButtonAll->Connect("Clicked()","AliDetectorFrame",this,"DoButton(Int_t)");
	fButtonFrame->AddFrame(fButtonAll,new TGLayoutHints(kLHintsBottom | kLHintsLeft,2,2,2,2));
	fButtonInvert = new TGTextButton(fButtonFrame,"Invert",kIdbSelectINVERT);
	fButtonInvert->Connect("Clicked()","AliDetectorFrame",this,"DoButton(Int_t)");
	fButtonFrame->AddFrame(fButtonInvert,new TGLayoutHints(kLHintsBottom | kLHintsRight,2,2,2,2));
	fMainFrame->AddFrame(fButtonFrame,new TGLayoutHints(kLHintsBottom | kLHintsLeft|kLHintsExpandX,0,0,2,2));		
	fCheckedMode = kTRUE;
}

//_____________________________________________________________
AliDetectorFrame::~AliDetectorFrame()
{
	for(Int_t i=0;i<gAliDisplay2->GetNbModules();i++){
		delete fCheckButton[i];
	}
	delete fCheckButton;
	delete fCheckedButton;
	delete fCheckButtonId;
	//delete [] fDetectorName;
	delete fButtonFrame;
	delete fMainFrame;
	delete fButtonAll;
	delete fButtonInvert;
}

//_____________________________________________________________
void AliDetectorFrame::DoButton(Int_t /*pos*/)
{
	TGFrame *frame = (TGFrame *) gTQSender;
	TGButton *bu= (TGButton *) frame;
	int id = bu->WidgetId();	
	fCheckedMode = kFALSE;
	AliModule *mo;
	switch(id){
	case kIdbSelectALL:{
			for(Int_t i=0;i<gAliDisplay2->GetNbModules();i++){
				mo = dynamic_cast<AliModule *> (gAliDisplay2->GetModules()->At(i));
				if(!mo) continue;
				fCheckButton[i]->SetState(kButtonDown);
				fCheckedButton[i]=kTRUE;
				gAliDisplay2->EnableDetector(mo->GetName());
			}
							  }
		break;
	case kIdbSelectINVERT:{
			for(Int_t i=0;i<gAliDisplay2->GetNbModules();i++){
				mo = dynamic_cast<AliModule *> (gAliDisplay2->GetModules()->At(i));
				if(!mo) continue;
				if(fCheckedButton[i]==kTRUE) {
					fCheckButton[i]->SetState(kButtonUp);
					fCheckedButton[i]=kFALSE;
					gAliDisplay2->DisableDetector(mo->GetName());
				}
				else if(fCheckedButton[i]==kFALSE)  {
					fCheckButton[i]->SetState(kButtonDown);
					fCheckedButton[i]=kTRUE;
					gAliDisplay2->EnableDetector(mo->GetName());
				}
			}
								  }
		break;
	default:break;
	}
	gAliDisplay2->Update(kmMODULES);
	fCheckedMode = kTRUE;
}

//_____________________________________________________________
void AliDetectorFrame::DoCheckButton(Int_t /*pos*/)
{
	if(fCheckedMode == kFALSE) return;
	TGFrame *frame = (TGFrame *) gTQSender;
	TGCheckButton *bu= (TGCheckButton *) frame;
	Int_t id = bu->WidgetId();
	AliModule *mo;
	for(Int_t i=0;i<gAliDisplay2->GetNbModules();i++){
		mo = dynamic_cast<AliModule *> (gAliDisplay2->GetModules()->At(i));
		if(!mo) continue;
		if(id==fCheckButtonId[i]){
			if(fCheckedButton[i]==kTRUE) {
				fCheckedButton[i]=kFALSE;
				gAliDisplay2->DisableDetector(mo->GetName());
				}
			else {
			fCheckedButton[i]=kTRUE;
			gAliDisplay2->EnableDetector(mo->GetName());
			}
		}
	}
	gAliDisplay2->Update(kmMODULES);
}

//_____________________________________________________________
void AliDetectorFrame::DoSpecific()
{
	TGFrame *frame = (TGFrame *) gTQSender;
	TGCheckButton *bu= (TGCheckButton *) frame;
	Int_t id = bu->WidgetId();
	AliModule *mo;
	for(Int_t i=0;i<gAliDisplay2->GetNbModules();i++){
		mo = dynamic_cast<AliModule *> (gAliDisplay2->GetModules()->At(i));
		if(!mo) continue;
		if(id==fCheckButtonId[i]){
			gAliDisplay2->DrawDetector(mo->GetName());
		}
	}
}


ClassImp(AliShutterItem);

//_____________________________________________________________
AliShutterItem::AliShutterItem(TGShutter *s,char *text, UInt_t id)
{
	fShutterItem = new TGShutterItem(s, new TGHotString(text), id);
	fMainFrame = (TGCompositeFrame *) fShutterItem->GetContainer();
	s->AddItem(fShutterItem);
}

//_____________________________________________________________
AliShutterItem::~AliShutterItem(void)
{
	delete fButton;
	delete fShutterItem;
	delete fMainFrame;
}

//_____________________________________________________________
void AliShutterItem::AddTextButton(char* text,char *tiptext, UInt_t idb)
{
	//Add a TGTextButton in the TGShutterItem. This button will execute the fonction
	fButton = new TGTextButton(fMainFrame,new TGHotString(text),idb);
	fButton->Resize(100,fButton->GetDefaultHeight());
	fButton->Connect("Clicked()","AliShutterItem",this,"DoButton(Int_t)");
	fButton->SetToolTipText(tiptext);
	//fButton->Connect("Clicked()","AliDisplay2",gAliDisplay2,"DoViews(Int_t)");
	fMainFrame->AddFrame(fButton, new TGLayoutHints( kLHintsTop | kLHintsCenterX ,5,5,10,10));
}

//_____________________________________________________________
void AliShutterItem::AddPictureButton(char* file, char *tiptext, UInt_t idb)
{
	//Add a TGPictureButton in the TGShutterItem. The icon file must be in DISPLAY/icons
	TString filename=StrDup(gAliDisplay2->GetIconsPath());
	filename.Append(file);
	TGPicture *picture = (TGPicture *) gClient->GetPicture(filename);
	fButton = new TGPictureButton(fMainFrame,picture,idb);		
	fButton->SetToolTipText(tiptext);
	fButton->Connect("Clicked()","AliShutterItem",this,"DoButton(Int_t)");
	fMainFrame->AddFrame(fButton, new TGLayoutHints( kLHintsTop | kLHintsCenterX ,5,5,10,10));
}

//_____________________________________________________________
void AliShutterItem::AddCheckButton(char *text,Int_t idb)
{
	fButton = new TGCheckButton(fMainFrame,new TGHotString(text),idb);
	fButton->Resize(100,fButton->GetDefaultHeight());
	fButton->Connect("Clicked()","AliShutterItem",this,"DoButton(Int_t)");
	fMainFrame->AddFrame(fButton, new TGLayoutHints( kLHintsTop | kLHintsLeft ,5,5,10,10));
}

//_____________________________________________________________
void AliShutterItem::DoButton(Int_t /*pos*/){
	TGFrame *frame = (TGFrame *) gTQSender;
	TGButton *bu= (TGButton *) frame;
	int id = bu->WidgetId();
	switch(id){
	case kIdbNextEVENT:{
		gAliDisplay2->ShowNextEvent(1);
							  }
		break;
	case kIdbPrevEVENT:{
		gAliDisplay2->ShowNextEvent(-1);
							  }
		break;
	case kIdbCheckHITS:{
	  if(gAliDisplay2->IsEnabled(kHits)) gAliDisplay2->Disable(kHits);
	  else gAliDisplay2->Enable(kHits);	  
		}
		break;
	case kIdbCheckCLUSTERS:{
	  if(gAliDisplay2->IsEnabled(kClusters)) gAliDisplay2->Disable(kClusters);
	  else gAliDisplay2->Enable(kClusters);
	}
		break;
	case kIdbCheckHLT:{
	  if(gAliDisplay2->IsEnabled(kHLT)) gAliDisplay2->Disable(kHLT);
	  else gAliDisplay2->Enable(kHLT);
	}
	  break;
	case kIdbCheckTRACKS:{
	  if(gAliDisplay2->IsEnabled(kTracks)) gAliDisplay2->Disable(kTracks);
	  else gAliDisplay2->Enable(kTracks);
	}
	  break;
	case kIdbSIDEVIEW:{
		gAliDisplay2->DoView(kIdbSIDEVIEW);
							}
		break;
	case kIdbFRONTVIEW:{
		gAliDisplay2->DoView(kIdbFRONTVIEW);
							 }
		break;
	case kIdbTOPVIEW:{
		gAliDisplay2->DoView(kIdbTOPVIEW);
						  }
		break;
	case kIdbALLVIEW:{
		gAliDisplay2->DoView(kIdbALLVIEW);
						  }
		break;
	default:break;
	}
}

ClassImp(AliShutterFrame);

//_____________________________________________________________
AliShutterFrame::AliShutterFrame(TGCompositeFrame *p, UInt_t /*w*/, UInt_t h){
	fShutter = new TGShutter(p,kSunkenFrame);
	fLayout = new TGLayoutHints(kLHintsExpandY | kLHintsTop | kLHintsLeft);
	fMainFrame = (TGCompositeFrame *) fShutter;
	
	//Event Shutter
	AliShutterItem *item = new AliShutterItem(fShutter,"Event",kIdsEVENT);
	
	item->AddPictureButton("next.xpm","Show next event",kIdbNextEVENT);
	item->AddPictureButton("prev.xpm","Show previous event",kIdbPrevEVENT);

	//View Shutter
	item = new AliShutterItem(fShutter,"View",kIdsVIEW);
	item->AddPictureButton("top.xpm","Top view",kIdbTOPVIEW);
	item->AddPictureButton("side.xpm","Side view",kIdbSIDEVIEW);
	item->AddPictureButton("front.xpm","Front view",kIdbFRONTVIEW);
	item->AddPictureButton("four.xpm","Four views",kIdbALLVIEW);

	//Detector Shutter
	item = new AliShutterItem(fShutter,"Detectors",kIdsDETECTORS);
	TGCompositeFrame *frame = item->GetShutterItemFrame();
	fDetectorFrameLayout = new TGLayoutHints( kLHintsTop | kLHintsLeft| kLHintsExpandX | kLHintsCenterX,5,5,5,5);
	fDetectorFrame = new AliDetectorFrame(frame,200,200,item->GetShutterItem()->GetDefaultFrameBackground());
	frame->AddFrame(fDetectorFrame->GetDetectorFrame(),fDetectorFrameLayout);
	
	//Options Shutter
	item = new AliShutterItem(fShutter,"Options",kIdsOPTIONS);
	item->AddCheckButton("Display Hits",kIdbCheckHITS);
	item->AddCheckButton("Display Clusters",kIdbCheckCLUSTERS);
	item->AddCheckButton("Display HLT Clusters",kIdbCheckHLT);
	//	item->AddCheckButton("Display Tracks",kIdbCheckTRACKS);

	fMainFrame->Resize(150,h);
}

//_____________________________________________________________
AliShutterFrame::~AliShutterFrame(void)
{
	delete fLayout;
	delete fShutter;
	delete fMainFrame;
	delete fDetectorFrame;
	delete fDetectorFrameLayout;
}


ClassImp(AliDisplayFrame);

//_____________________________________________________________
AliDisplayFrame::AliDisplayFrame(const TGWindow *p, UInt_t w, UInt_t h)
{
	fClipMin=-20;
	fClipMax=20;
	fPreviousW=0;
	fPreviousH=0;
	fRange = 500;
	fPolyMarkers = new TObjArray(1000);
	
	fMainFrame = new TGCompositeFrame(p,w,h);
	fMainTab = new TGTab(fMainFrame, w, h);
	TGCompositeFrame *fFrame1 = fMainTab->AddTab("Main View");
	fMainEmbeddedCanvas = new TRootEmbeddedCanvas("Main12",fFrame1,w,h,kFixedWidth);
	fFrame1->AddFrame(fMainEmbeddedCanvas,new TGLayoutHints( kLHintsTop | kLHintsLeft|kLHintsExpandX| kLHintsExpandY, 0, 0, 0, 0));
	fMainCanvas = fMainEmbeddedCanvas->GetCanvas();
	fMainCanvas->SetFillColor(1);
	fMainCanvas->SetBorderMode(0);
	fMainCanvas->cd();
	fMainCanvas->SetFixedAspectRatio();
	fMainCanvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","AliDisplayFrame",this,"ExecuteEvent(Int_t,Int_t,Int_t,TObject*)");
	//fView = new TView(1);
	//DoView(kIdbFRONTVIEW);
	
	gAliDisplay2->SetCurrentView(kIdbFRONTVIEW);	
	
	TGCompositeFrame  *fFrame2 = fMainTab->AddTab("No detector");
	fSelectionEmbeddedCanvas = new TRootEmbeddedCanvas("Selection",fFrame2,w,h);
	fSelectionCanvas = fSelectionEmbeddedCanvas->GetCanvas();
	fSelectionCanvas->SetFillColor(1);
	fSelectionCanvas->SetBorderMode(0);
	fSelectionCanvas->cd();
	fFrame2->AddFrame(fSelectionEmbeddedCanvas,new TGLayoutHints( kLHintsTop | kLHintsLeft|kLHintsExpandX| kLHintsExpandY, 0, 0, 0, 0));
	fMainFrame->AddFrame(fMainTab,new TGLayoutHints( kLHintsTop | kLHintsLeft|kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
	fAllViews = kFALSE;
	fMainFrame->MapSubwindows();
	fMainFrame->MapWindow();
}

//_____________________________________________________________
AliDisplayFrame::~AliDisplayFrame(void)
{
	delete fMainTab;
	delete fSelectionEmbeddedCanvas;
	delete fMainEmbeddedCanvas;
	delete fFrame1;
	delete fFrame2;
	delete fMainCanvas;
	delete fSelectionCanvas;
	delete fPoints2;
	delete fPoints;
	delete fModules;
	delete fMainFrame;
	delete [] fActivePoints;
	delete [] fClustersPos;
}

//_____________________________________________________________
void AliDisplayFrame::DoView(Int_t view)
{
	Int_t x,y;
	char vname[16];
	y=fMainFrame->GetDefaultHeight();
	x=fMainFrame->GetDefaultWidth();
	gAliDisplay2->SetCurrentView(view);
	switch(view){
	case kIdbALLVIEW:{
		fAllViews=kTRUE;
		strcpy(vname,"All views");	
		fMainCanvas->cd();
		gPad->Clear();
		fMainCanvas->SetFillColor(15);
		fMainCanvas->Divide(2,2,0.005,0.005,1);
		
		fMainCanvas->cd(1);
		Draw(30,30,0);

		gAliDisplay2->SetCurrentView(kIdbTOPVIEW);
		fMainCanvas->cd(2);
		Draw(90,-90,90);

		gAliDisplay2->SetCurrentView(kIdbSIDEVIEW);
		fMainCanvas->cd(3);		
		Draw(90,0,-90);
		
		gAliDisplay2->SetCurrentView(kIdbFRONTVIEW);
		fMainCanvas->cd(4);
		Draw(0,-90,0);
		
		//fMainCanvas->cd();
		
	}
		break;
	case kIdbTOPVIEW:{
		strcpy(vname,"Top view  ");
		fAllViews=kFALSE;	
		fMainCanvas->cd();											
		gPad->SetFillColor(1);
		gPad->Clear();
		gPad->Draw();
		Draw(90,-90,90);
	}
		break;
	case kIdbSIDEVIEW:{
		strcpy(vname,"Side view");
		fAllViews=kFALSE;	
		fMainCanvas->cd();		
		gPad->SetFillColor(1);
		gPad->Clear();
		gPad->Draw();
		Draw(90,0,-90);
	}
		break;
	case kIdbFRONTVIEW:{
		strcpy(vname,"Front view");
		fAllViews=kFALSE;	
		fMainCanvas->cd();
		gPad->SetFillColor(1);
		gPad->Clear();
		gPad->Draw();
		
		Draw(0,-90,0);
	}
		break;
	default: break;
	}
	(fMainTab->GetTabTab(0))->SetText(new TGString(vname));
}

//_____________________________________________________________
void AliDisplayFrame::DrawDetector(const char *name)
{
	(fMainTab->GetTabTab(1))->SetText(new TGString(name));
}

//_____________________________________________________________
void AliDisplayFrame::EnableDetector(const char *name)
{
	AliModule *module = dynamic_cast<AliModule*>(gAlice->Modules()->FindObject(name));
	if(!module) return;
	gAliDisplay2->GetModuleInfo()->Enable((char*)name);
	module->Enable();
}

//_____________________________________________________________
void AliDisplayFrame::DisableDetector(const char *name)
{
	AliModule *module = dynamic_cast<AliModule*>(gAlice->Modules()->FindObject(name));
	if(!module) return;
	gAliDisplay2->GetModuleInfo()->Disable((char*)name);
	module->Disable();
}

//_____________________________________________________________
void AliDisplayFrame::Draw(Float_t theta, Float_t phi, Float_t psi)
{
  //clock_t t1,t2;
  time_t t1,t2;
  //t1 = clock();
  TStopwatch timer;
  timer.Start();
  time(&t1);
	gPad->SetCursor(kWatch);
	gPad->SetEditable(kTRUE);
	gPad->SetFillColor(1);
	gPad->Clear();
	
	Int_t iret;
	
	TView *view = new TView(1);
	TGDimension dim=((TGCanvas*)fMainEmbeddedCanvas)->GetViewPort()->GetDefaultSize();
	Float_t aspectRatio = dim.fWidth/(Float_t) dim.fHeight;
	//printf("Dimension %d %d",dim.fWidth,dim.fHeight);
	if(gAliDisplay2->GetCurrentView()==kIdbFRONTVIEW){
		view->SetRange(-fRange*aspectRatio,-fRange,-fRange,fRange*aspectRatio,fRange,fRange);
	}
	if(gAliDisplay2->GetCurrentView()==kIdbTOPVIEW){
		view->SetRange(-fRange,-fRange,-fRange*aspectRatio,fRange,fRange,fRange*aspectRatio);
	}
	if(gAliDisplay2->GetCurrentView()==kIdbSIDEVIEW){
		view->SetRange(-fRange,-fRange,-fRange*aspectRatio,fRange,fRange,fRange*aspectRatio);
	}
	
	gAlice->GetGeometry()->Draw("same");
	if(gAliDisplay2->IsEnabled(kHits)) DrawHits();
	if(gAliDisplay2->IsEnabled(kClusters)) fClusters->Draw();
	if(gAliDisplay2->IsEnabled(kHLT)) fHLT->Draw();

	gAliDisplay2->AppendPad();
	view->SetView(phi,theta,psi,iret);
	
	view->ZoomView(gPad,gAliDisplay2->GetZoomFactor());
	//t2 = clock();
	time(&t2);
	//	printf("\nDrawn in....%f sec", ((double)t2-t1)/(10000*CLK_TCK));
	printf("\nDrawn in....%2lf sec", difftime(t2,t1));
	timer.Stop();
	timer.Print("m");
}

//_____________________________________________________________
void AliDisplayFrame::DrawHits()
{
	AliPoints *p;
	if(!fPoints2) return;
	for(Int_t i=0;i<fPoints2->GetEntries();i++){
		if(fActivePoints[i]){
			p=dynamic_cast<AliPoints *>(fPoints2->UncheckedAt(i));	
			if(!p) continue;
			p->Draw();
		}
	}
}

//_____________________________________________________________
void AliDisplayFrame::LoadEnabledModules()
{
	clock_t t1,t2;
	t1=clock(); 
	TIter next(gAlice->Modules());
	AliModule *module;
	fModules = new TObjArray(0,32);
	while((module = dynamic_cast <AliModule*> (next()))){
		if(!module) continue;
		if(!module->IsActive()) continue;
		fModules->AddLast(module);
	}
	t2=clock();
	fNbModules = fModules->GetEntriesFast();
	//	printf("\nModules loaded in.....%f sec", ((double)t2-t1)/(10000*CLK_TCK));
}

//_____________________________________________________________
void AliDisplayFrame::LoadClusters(Int_t nevent)
{
  	fClusters = new AliDisplayClusters();
	fClusters->LoadClusters("ITS TPC",nevent);
}

//_____________________________________________________________
void AliDisplayFrame::LoadHLTClusters(Int_t nevent)
{
  fHLT = new AliDisplayHLT();
  fHLT->LoadHLT("TPC",nevent);
}
	
//_____________________________________________________________
void AliDisplayFrame::LoadHits()
{
	clock_t t1,t2;

	t1=clock(); 
	fPoints2 = new TObjArray(0,1000);
	AliModule *module;
	TObjArray *points;
	for(Int_t i=0;i<fNbModules;i++){
		module = dynamic_cast<AliModule*>(fModules->UncheckedAt(i));
		if(!module) continue;
		points = module->Points();
		if(!points) {
			continue;
		}
		for(Int_t j=0;j<points->GetEntriesFast();j++){
			if(!points->UncheckedAt(j)) continue;
			fPoints2->AddLast((points->UncheckedAt(j)));
		}	
	}
	fActivePoints = new Bool_t[fPoints2->GetEntries()];
	for(Int_t k=0;k<fPoints2->GetEntriesFast();k++){
		fActivePoints[k]=kTRUE;
	}
	printf("\n nb hits %d",fPoints2->GetEntries());
	t2=clock();
	//	printf("\nPoints loaded in....%f sec", ((double)t2-t1)/(10000*CLK_TCK));
}

//_____________________________________________________________
void AliDisplayFrame::ApplyCuts()
{
	clock_t t1,t2;
	t1=clock();
	
	Float_t		*pxyz;
	Float_t		r,theta,eta,cutmin,cutmax,etamin,etamax,pmom,smin,smax;
	Int_t		nbhits=0;
	AliPoints *pm;
	TParticle *particle;
	
	//Get momentum cut
	smin = gAliDisplay2->GetMomentumMin();
	smax = gAliDisplay2->GetMomentumMax();
	cutmin = 2.0*smin;
	if(smax<0.98) 	cutmax = 2.0*smax;
	else 			cutmax = 100000;
	
	//Get rapidity cut
	smax = gAliDisplay2->GetRapidityMax();
	smin = gAliDisplay2->GetRapidityMin();
	//etamin = 1.5*(2*smin-1);
	//etamax = 1.5*(2*smax-1);
	etamin = smin;
	etamax = smax;
	if(smin<-1.46) etamin = -1000;
	if(smax>1.46) etamax = 1000;
	

	if(!fPoints2) return;
	for(Int_t i=0;i<fPoints2->GetEntries();i++){
		pm = dynamic_cast<AliPoints*>(fPoints2->UncheckedAt(i));
		if(!pm) {
			fActivePoints[i]=kFALSE;
			continue;
		}
		particle = pm->GetParticle();
		if(!particle) {
			fActivePoints[i]=kFALSE;
			continue;
		}
		pmom = particle->P();
		if(pmom < cutmin) {
			fActivePoints[i]=kFALSE;
			continue;
		}
		if(pmom > cutmax) {
			fActivePoints[i]=kFALSE;
			continue;
		}
		pxyz = pm->GetP();
		r = TMath::Sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]);
		theta = TMath::ATan2(r,TMath::Abs(pxyz[2]));
		if(theta) eta = -TMath::Log(TMath::Abs(TMath::Tan(0.5*theta)));
		else eta = 1e10;
		if(pxyz[2] < 0) eta = -eta;
		if((eta < etamin) || (eta > etamax)) {
			fActivePoints[i]=kFALSE;
			continue;
		}
		fActivePoints[i]=kTRUE;
		//pm->Draw();
		nbhits += pm->GetN();
	}
	gAliDisplay2->SetNbHits(nbhits);
	t2=clock();
	//	printf("\nCuts applied in....%f sec", ((double)t2-t1)/(10000*CLK_TCK));
	gAliDisplay2->SetNbParticles(GetNbActivePoints()); 
}

//_____________________________________________________________
Int_t AliDisplayFrame::GetNbActivePoints()
{
	Int_t ans=0;
	for(Int_t i=0;i<fPoints2->GetEntries();i++){
		if(fActivePoints[i]) ans++;
	}
	return ans;
}
//_____________________________________________________________
void AliDisplayFrame::DrawX3d()
{
	TPad *pad = dynamic_cast<TPad*>(gPad);
	pad->cd();
	TView *view = pad->GetView();
	if(!view) return;
	pad->x3d();
}

//_____________________________________________________________
void AliDisplayFrame::SavePadGIF(const char *file)
{
  if(!gPad){
    printf("\nThere is no active pad");
    return;
  }
  gPad->SaveAs(file);
}

//_____________________________________________________________
void AliDisplayFrame::DrawGL()
{
	TPad *pad = dynamic_cast<TPad*>(gPad);
	pad->cd();
	TView *view = pad->GetView();
	if(!view) return;
	pad->x3d("OPENGL");
}

//_____________________________________________________________
void AliDisplayFrame::ExecuteEvent(Int_t event, Int_t px,Int_t py,TObject *)
{
	static Float_t x0,y0,x1,y1;
	static Int_t pxold,pyold;
	static Int_t px0,py0;
	static Int_t linedrawn;
	Float_t temp;
	
	
	switch(event){
	case kMouseMotion:{
			
			AliPoints *p=dynamic_cast<AliPoints*> (gPad->GetSelected());
			if(p){
				gAliDisplay2->SetStatusBar(p->GetName(),1);
				gAliDisplay2->SetStatusBar(p->GetDetector()->GetName(),2);
				}
			}
			break;
		default:break;
	}	
	
	if((!gAliDisplay2->GetZoomMode())&&(gPad->GetView())){
		gPad->GetView()->ExecuteRotateView(event,px,py);
		return;
	}
	

	
	if(gAliDisplay2->GetZoomMode()==kTRUE){
		switch(event){
			
		case kButton1Down:{
			gVirtualX->SetLineColor(-1);
			gPad->TAttLine::Modify();
			x0 = gPad->AbsPixeltoX(px);
			y0 = gPad->AbsPixeltoY(py);
			px0 = px;
			py0 = py;
			pxold = px;
			pyold = py;
			linedrawn = 0;
		}
		break;
		case kButton1Motion:{
			if(linedrawn) gVirtualX->DrawBox(px0,py0,pxold,pyold,TVirtualX::kHollow);
			pxold = px;
			pyold = py;
			linedrawn = 1;
			gVirtualX->DrawBox(px0,py0,pxold,pyold,TVirtualX::kHollow);
		}
		break;
		
		case kButton1Up:{
			gPad->GetCanvas()->FeedbackMode(kFALSE);
			if(px == px0) break;
			if(py == py0) break;
			x1 = gPad->AbsPixeltoX(px);
			y1 = gPad->AbsPixeltoY(py);
			if(x1<x0) { 
				temp = x0;
				x0 = x1;
				x1 = temp;
			}
			if(y1<y0) {
				temp = y0;
				y0 = y1;
				y1 = temp;
			}
			printf("\nBox (%f,%f)-(%f,%f)",x0,y0,x1,y1);
			gPad->SetEditable(kTRUE);
			//gPad->Range(x0,y0,x1,y1);
			gPad->SetEditable(kFALSE);
			//gPad->Range(0.5,0.5,1,1);
			//gAliDisplay2->SetZoomFactor(1);
				gPad->Modified(kTRUE);
			gAliDisplay2->Draw();	
			gAliDisplay2->SetZoomMode(kFALSE);
			gPad->SetEditable(kTRUE);
		}
		break; 
		default: break;		
		}		
	}
}

ClassImp(AliInfoFrame);

//_____________________________________________________________
AliInfoFrame::AliInfoFrame(TGCompositeFrame *p, UInt_t w, UInt_t h)
{
	fMainFrame = new TGCompositeFrame(p, w, h, kVerticalFrame);

	fTitleFrame = new TGCompositeFrame(fMainFrame,w,h,kRaisedFrame|kVerticalFrame);
	AddLabel("ALICE",kLHintsTop | kLHintsCenterX);
	AddLabel("Event Display",kLHintsTop | kLHintsCenterX);
	
	TString filename=StrDup(gAliDisplay2->GetIconsPath());
	filename.Append("Alice.xpm");
	TGPicture *alicelogo = (TGPicture *) gClient->GetPicture(filename);
	TGIcon *alice = new TGIcon(fTitleFrame,alicelogo,50,50);
	fTitleFrame->AddFrame(alice,new TGLayoutHints(kLHintsTop | kLHintsCenterX,0,0,0,0));
	
	AddLabel("Powered by",kLHintsTop | kLHintsCenterX);
	AddLabel("AliRoot",kLHintsTop | kLHintsCenterX);
	fMainFrame->AddFrame(fTitleFrame,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,0,0,0,3));

	//Feedback
	fFiguresFrame = new TGCompositeFrame(fMainFrame,w,h,kRaisedFrame|kVerticalFrame);
	TGCompositeFrame *frame = new TGCompositeFrame(fFiguresFrame,w,100, kHorizontalFrame);
	fNbEventLabel = new TGLabel(frame,"");
	TGLabel * label = new TGLabel(frame,"Event number");
	fNbEventLabel->SetText(gAliDisplay2->GetEventNumber());
	frame->AddFrame(label,new TGLayoutHints(kLHintsTop | kLHintsLeft ,10,0,0,0));
	frame->AddFrame(fNbEventLabel,new TGLayoutHints(kLHintsTop | kLHintsRight,5,10,0,0));

	fFiguresFrame->AddFrame(frame,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,0,0,0,0));
	
	frame = new TGCompositeFrame(fFiguresFrame,w,100, kHorizontalFrame);
	label = new TGLabel(frame,"Nb Particles");
	fNbParticuleLabel = new TGLabel(frame,"");
	fNbParticuleLabel->SetText(gAliDisplay2->GetNbParticles());
	frame->AddFrame(label,new TGLayoutHints(kLHintsTop | kLHintsLeft,10,0,0,0));
	frame->AddFrame(fNbParticuleLabel,new TGLayoutHints(kLHintsTop | kLHintsRight,5,10,0,0));

	fFiguresFrame->AddFrame(frame,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,0,0,0,0));

	frame = new TGCompositeFrame(fFiguresFrame,w,100, kHorizontalFrame);
	label = new TGLabel(frame,"Nb Hits");
	fNbHitsLabel = new TGLabel(frame,"");
	fNbHitsLabel->SetText("--");
	frame->AddFrame(label,new TGLayoutHints(kLHintsTop | kLHintsLeft,10,0,0,0));
	frame->AddFrame(fNbHitsLabel,new TGLayoutHints(kLHintsTop | kLHintsRight ,5,10,0,0));
	frame->Layout();
	fFiguresFrame->AddFrame(frame,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,0,0,0,0));

	frame = new TGCompositeFrame(fFiguresFrame,w,100, kHorizontalFrame);
	label = new TGLabel(frame,"Nb Clusters");
	fNbClustersLabel = new TGLabel(frame,"");
	fNbClustersLabel->SetText("--");
	frame->AddFrame(label,new TGLayoutHints(kLHintsTop | kLHintsLeft,10,0,0,0));
	frame->AddFrame(fNbClustersLabel,new TGLayoutHints(kLHintsTop | kLHintsRight ,5,10,0,0));
	frame->Layout();
	fFiguresFrame->AddFrame(frame,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,0,0,0,0));
	fMainFrame->AddFrame(fFiguresFrame,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,0,0,2,0));

	fMainFrame->Layout();
	fMainFrame->MapSubwindows();
	fMainFrame->MapWindow();
};

//_____________________________________________________________
AliInfoFrame::~AliInfoFrame(void){

	delete fMainFrame;
	delete fTitleFrame;
	delete fFiguresFrame;
	delete fNbParticuleLabel;
	delete fNbEventLabel;
	delete fNbHitsLabel;
}

//_____________________________________________________________
void AliInfoFrame::AddLabel(char *text, UInt_t options){
	TGLabel * label = new TGLabel(fTitleFrame,text);
	fTitleFrame->AddFrame(label,new TGLayoutHints(options,0,0,0,0));
}

//_____________________________________________________________
void AliInfoFrame::Update()
{
	fNbParticuleLabel->SetText(gAliDisplay2->GetNbParticles());
	fNbEventLabel->SetText(gAliDisplay2->GetEventNumber());
	if(gAliDisplay2->IsEnabled(kHits))fNbHitsLabel->SetText(gAliDisplay2->GetNbHits());
	else fNbHitsLabel->SetText("--");
	if(gAliDisplay2->IsEnabled(kClusters))fNbClustersLabel->SetText(gAliDisplay2->GetNbClusters());
	else fNbClustersLabel->SetText("--");
	fMainFrame->Layout();
}


ClassImp(AliSettingFrame);

//_____________________________________________________________
AliSettingFrame::AliSettingFrame(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h)
				 :TGTransientFrame(p,main,w,h)
{
	fMainFrame = new TGCompositeFrame((TGWindow *)((TGTransientFrame *)this),w,h,kVerticalFrame);
	
	fZoomStepFrame = new TGCompositeFrame(fMainFrame,w,50,kHorizontalFrame);
	fZoomStepLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft |kLHintsExpandX,5,5,5,5);
	fZoomStepEntry = new TGNumberEntryField(fZoomStepFrame,kIdtZoomSTEP,gAliDisplay2->GetZoomStep());
	fZoomStepEntry->Connect("ReturnPressed()","AliSettingFrame",this,"DoSettings(Int_t)");
	fZoomStepLabel = new TGLabel(fZoomStepFrame,"Zoom step");
	fZoomStepFrame->AddFrame(fZoomStepLabel,new TGLayoutHints(kLHintsTop | kLHintsLeft,0,0,0,0));
	fZoomStepFrame->AddFrame(fZoomStepEntry,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX ,5,5,0,0));
	fMainFrame->AddFrame(fZoomStepFrame,fZoomStepLayout);

	fSliderStepFrame = new TGCompositeFrame(fMainFrame,w,50,kHorizontalFrame);
	fSliderStepLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft |kLHintsExpandX,5,5,5,5);
	fSliderStepEntry = new TGNumberEntryField(fSliderStepFrame,kIdtSliderSTEP,gAliDisplay2->GetSliderStep());
	fSliderStepEntry->Connect("ReturnPressed()","AliSettingFrame",this,"DoSettings(Int_t)");
	fSliderStepLabel = new TGLabel(fSliderStepFrame,"Slider step");
	fSliderStepFrame->AddFrame(fSliderStepLabel,new TGLayoutHints(kLHintsTop | kLHintsLeft,0,0,0,0));
	fSliderStepFrame->AddFrame(fSliderStepEntry,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX ,5,5,0,0));
	fMainFrame->AddFrame(fSliderStepFrame,fSliderStepLayout);

	fSliderUpdateFrame = new TGCompositeFrame(fMainFrame,w,50,kHorizontalFrame);
	fSliderUpdateLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft |kLHintsExpandX,5,5,5,5);
	fSliderUpdateButton = new TGCheckButton(fSliderUpdateFrame,"Update display on slider move",kIdtSliderUPDATE);
	fSliderUpdateButton->Connect("Clicked()","AliSettingFrame",this,"DoSettings(Int_t)");
	fIsLoading = kTRUE;
	if(gAliDisplay2->GetSliderUpdate()) fSliderUpdateButton->SetState(kButtonDown);
	else fSliderUpdateButton->SetState(kButtonUp);
	fIsLoading = kFALSE;
       
	fSliderUpdateFrame->AddFrame(fSliderUpdateButton,new TGLayoutHints(kLHintsTop | kLHintsLeft,0,0,0,0));
	fMainFrame->AddFrame(fSliderUpdateFrame,fSliderUpdateLayout);

	AddFrame(fMainFrame,new TGLayoutHints(kLHintsTop | kLHintsLeft |kLHintsExpandX,0,0,0,0));
	fMainFrame->Layout();
	// position relative to the parent's window
   Window_t wdum;
   int ax, ay;
   gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(),
             (Int_t)(((TGFrame *) main)->GetWidth() - GetWidth()) >> 1,
             (Int_t)(((TGFrame *) main)->GetHeight() - GetHeight()) >> 1,
             ax, ay, wdum);
   Move(ax, ay);

   SetWindowName("Setting frame");
   MapSubwindows();
   MapWindow();
   Layout();
}

//_____________________________________________________________
AliSettingFrame::~AliSettingFrame()
{
	delete fZoomStepLayout;
	delete fZoomStepEntry;
	delete fZoomStepLabel;
	delete fSliderStepLayout;
	delete fSliderStepEntry;
	delete fSliderStepLabel;	

	delete fSliderUpdateLayout;
	delete fSliderUpdateButton;
	
	delete fMainFrame;
	delete fZoomStepFrame;
	delete fSliderUpdateFrame;
	delete fSliderStepLayout;
}

//_____________________________________________________________
void AliSettingFrame::DoSettings(Int_t /*pos*/)
{
	TGNumberEntryField *ne = (TGNumberEntryField *) gTQSender;
	int id = ne->WidgetId();
	switch(id){
	case kIdtZoomSTEP:{
		gAliDisplay2->SetZoomStep(ne->GetNumber());
		}
		break;
	case kIdtSliderSTEP:{
		gAliDisplay2->SetSliderStep(ne->GetNumber());
		}
		break;
	case kIdtSliderUPDATE:{
	  if(fIsLoading) return ;
	  if(gAliDisplay2->GetSliderUpdate()) gAliDisplay2->SetSliderUpdate(kFALSE);
	  else gAliDisplay2->SetSliderUpdate(kTRUE);
	}
	  break;
	default: break;
	}
}


ClassImp(AliMenu);

//_____________________________________________________________
AliMenu::AliMenu(TGCompositeFrame *p, UInt_t w, UInt_t h, UInt_t options)
{
	fMenuBar = new TGMenuBar(p,w,h,options);
	fToolBar = new TGToolBar(p,60,20,options);
	
	fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsExpandX | kLHintsLeft ,0,0,0,0);
	fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft,0,4,0,0);
	
	fMenuFile = new TGPopupMenu(gClient->GetRoot());
	fMenuFile->AddEntry("Open",kIdmOPEN);
	fMenuFile->AddEntry("Save as",kIdmSAVEAS);
	fMenuFile->AddEntry("Close",kIdmCLOSE);
	fMenuFile->AddSeparator();
	fMenuFile->AddEntry("Print",kIdmPRINT);
	fMenuFile->AddEntry("Print setup",kIdmPRINTSETUP);
	fMenuFile->AddSeparator();
	fMenuFile->AddEntry("Exit",kIdmEXIT);
	fMenuFile->DisableEntry(kIdmSAVEAS);
	fMenuFile->Associate(p);
	fMenuBar->AddPopup("File",fMenuFile,fMenuBarItemLayout);
	fMenuFile->Connect("Activated(Int_t)","AliMenu",this,"DoMenu(Int_t)");

	fMenuOptions = new TGPopupMenu(gClient->GetRoot());
	fMenuOptions->AddEntry("Settings",kIdmSETTINGS);
	fMenuOptions->AddEntry("Save settings",kIdmSAVESETTINGS);
	fMenuOptions->Associate(p);
	fMenuBar->AddPopup("Options",fMenuOptions,fMenuBarItemLayout);
	fMenuOptions->Connect("Activated(Int_t)","AliMenu",this,"DoMenu(Int_t)");
	
	fMenuView = new TGPopupMenu(gClient->GetRoot());
	fMenuView->AddEntry("X3d ",kIdmVIEWX3D);
	fMenuView->AddEntry("OpenGL",kIdmVIEWGL);
	fMenuView->Associate(p);
	fMenuBar->AddPopup("View",fMenuView,fMenuBarItemLayout);
	fMenuView->Connect("Activated(Int_t)","AliMenu",this,"DoMenu(Int_t)");
	
	fMenuHelp = new TGPopupMenu(gClient->GetRoot());
	fMenuHelp->AddEntry("Help",kIdmHELP);
	fMenuHelp->AddSeparator();
	fMenuHelp->AddEntry("About",kIdmABOUT);
	fMenuHelp->Associate(p);
	fMenuBar->AddPopup("Help",fMenuHelp,fMenuBarItemLayout);
	fMenuHelp->Connect("Activated(Int_t)","AliMenu",this,"DoMenu(Int_t)");

	p->AddFrame(fMenuBar,fMenuBarLayout);
	fTBD = new ToolBarData_t;
	
	fToolBarLayout = new TGLayoutHints(kLHintsTop | kLHintsExpandX,0,0,0,0);
	AddPictureButton("open.xpm","Open file",kIdmOPEN,5);
	AddPictureButton("save.xpm","Save current pad as gif file",kIdmSAVEAS,0);
	AddPictureButton("settings.xpm","Settings",kIdmSETTINGS,5);
	AddPictureButton("help.xpm","Help",kIdmHELP,5);
	AddPictureButton("quit.xpm","Exit AliDisplay",kIdmEXIT,5);
	AddPictureButton("opengl.xpm","Open GL view",kIdmVIEWGL,5);
	AddPictureButton("x3d.xpm","x3d view",kIdmVIEWX3D,0);
	AddPictureButton("zoomplus16.xpm","Zoom in",kIdbZoomIN,5);
	AddPictureButton("zoommoins16.xpm","Zoom out",kIdbZoomOUT,0);
	AddPictureButton("zoomzone.xpm","Zoom on zone",kIdbZoomZONE,0);
	p->AddFrame(fToolBar,fToolBarLayout);
}

//_____________________________________________________________
AliMenu::~AliMenu()
{
	delete fMenuBarLayout;
	delete fMenuBarItemLayout;
	delete fMenuFile;
	delete fMenuOptions;
	delete fMenuView;
	delete fMenuHelp;
	delete fToolBarLayout;
	delete fToolBar;
	delete fMenuBar;
	delete fTBD;
}

//_____________________________________________________________
void AliMenu::DoMenu(Int_t id)
{
	switch(id){
	case kIdmOPEN:{
		TGFileInfo fi;
		static TString dir(".");
		fi.fFileTypes = filetypes;
		fi.fIniDir = StrDup(dir.Data());
		new TGFileDialog(gClient->GetRoot(),gAliDisplay2->GetMainFrame(),kFDOpen,&fi);
		if(!fi.fFilename) return;
					  }
		break;
	case kIdmEXIT:{
		gApplication->Terminate(0);
					  }
		break;
	case kIdmSAVEAS:{
	  TGFileInfo fi;
	  static TString dir(".");
	  fi.fFileTypes = imgtypes;
	  fi.fIniDir = StrDup(dir.Data());
	  new TGFileDialog(gClient->GetRoot(),gAliDisplay2->GetMainFrame(),kFDSave,&fi);
	  if(!fi.fFilename) return;
	  gAliDisplay2->SavePadGIF(fi.fFilename);
	}
	  break;
	case kIdmSETTINGS:{
		new AliSettingFrame((TGWindow *)gClient->GetRoot(),(TGWindow *)gAliDisplay2->GetMainFrame(),200,150);
							}
		break;
	case kIdmHELP:{
		TRootHelpDialog *hd=new TRootHelpDialog((TGWindow *)gClient->GetRoot(),"Help",300,300);
		hd->SetText(helpTxt);	
        hd->Popup();
		}
		break;
		
	case kIdmSAVESETTINGS:{
		gAliDisplay2->DoSaveSettings();
								 }
		break;
	case kIdmVIEWX3D:{
		gAliDisplay2->DrawX3d();
		}
		break;
	case kIdmVIEWGL:{
		gAliDisplay2->DrawGL();
	}
	  break;
	case kIdbZoomIN:{
		gAliDisplay2->SetZoomFactor(gAliDisplay2->GetZoomFactor()*gAliDisplay2->GetZoomStep());
		gAliDisplay2->Draw();
						  }
		break;
	case kIdbZoomZONE:{
		gAliDisplay2->SetZoomMode(kTRUE);
		gAliDisplay2->SetEditable(kFALSE);
	}
	break;
	case kIdbZoomOUT:{		
		gAliDisplay2->SetZoomFactor(gAliDisplay2->GetZoomFactor()/gAliDisplay2->GetZoomStep());
		gAliDisplay2->Draw();
							}
		break;
	default:break;
	}
}

//_____________________________________________________________
void AliMenu::DoToolBar(Int_t /*id*/)
{
	TGFrame *frame = (TGFrame *) gTQSender;
	TGButton *bu = (TGButton *) frame;
	DoMenu(bu->WidgetId());
}

//_____________________________________________________________
void AliMenu::AddPictureButton(char *fname,char *tiptext,UInt_t id, UInt_t spacing)
{
	TString filename = StrDup(gAliDisplay2->GetIconsPath());
	filename.Append(fname);

	fTBD->fPixmap=filename.Data();
	fTBD->fTipText = tiptext;
	fTBD->fId = id;
	fTBD->fStayDown = kFALSE;

	fToolBar->AddButton(fToolBar,fTBD,spacing);
	if(fTBD->fButton)
		fTBD->fButton->Connect("Clicked()","AliMenu",this,"DoToolBar(Int_t)");
}

ClassImp(AliDisplay2);
//_____________________________________________________________
AliDisplay2::AliDisplay2(const TGWindow *p, UInt_t w, UInt_t h)
			:TObject()
{
	//gAlice->SetDisplay(this);
	gAliDisplay2=this;
	fSliderUpdate = kFALSE;
	fZoomMode = kFALSE;
	fZoomStep = 1.2;
	fZoomFactor = 1.5;
	fNbParticles = 0;
	fEventNumber = 0;
	fNbHits = 0;
	fSliderStep = 0.01;
	fClustersLoaded = kFALSE;
	fHitsLoaded = kFALSE;
	fHLTLoaded = kFALSE;
	fTracksLoaded = kFALSE;
	fMode =0;
	FindModules();

	fIconsPath = new char[32];
	strcpy(fIconsPath,gSystem->Getenv("ALICE_ROOT"));
	strcat(fIconsPath,"/DISPLAY/icons/");
	LoadFromRC();
	fMainFrame = new TGMainFrame(p,w,h,kVerticalFrame);
	fSubFrame = new TGCompositeFrame(fMainFrame,w,h,kHorizontalFrame);
	fLeftFrame = new TGCompositeFrame(fSubFrame,150,h,kVerticalFrame|kFixedWidth);
	fRightFrame = new TGCompositeFrame(fSubFrame,600,h,kVerticalFrame);
	//fMainFrame->Connect("ProcessedEvent(Event_t*)", "AliDisplay2", this,"HandleMouseWheel(Event_t*)");
	fMainFrame->Connect("ProcessedEvent(Event_t*)", "AliDisplay2",this,"HandleResize(Event_t*)");
	//MenuBar
	fMenu = new AliMenu(fMainFrame,1,1,kRaisedFrame|kHorizontalFrame);

	//Slider Frame
	fSliderFrameLayout = new TGLayoutHints( kLHintsBottom| kLHintsRight| kLHintsExpandX | kLHintsCenterX, 2, 2, 2, 2);
	fSliderFrame  = new AliSliderFrame(fRightFrame,600,150);
	fRightFrame->AddFrame(fSliderFrame->GetSliderFrame(),fSliderFrameLayout);

		//Info Frame
	fInfoFrameLayout = new TGLayoutHints( kLHintsTop | kLHintsLeft  | kLHintsExpandX  ,0,0,0,0);
	fInfoFrame = new AliInfoFrame(fLeftFrame,150,200);
	fLeftFrame->AddFrame(fInfoFrame->GetInfoFrame(),fInfoFrameLayout);
	
	
	//Shutter Frame
	fShutterFrameLayout = new TGLayoutHints( kLHintsTop | kLHintsLeft | kLHintsExpandY |kLHintsExpandX, 0, 0, 5, 0);
	fShutterFrame = new AliShutterFrame(fLeftFrame,150,300);
	fLeftFrame->AddFrame(fShutterFrame->GetShutterFrame(),fShutterFrameLayout);

	//Display Frame
	fDisplayFrameLayout = new TGLayoutHints( kLHintsTop | kLHintsRight | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2); 
	fDisplayFrame = new AliDisplayFrame(fRightFrame, w-150,w-110);
	fRightFrame->AddFrame(fDisplayFrame->GetDisplayFrame(),fDisplayFrameLayout);
	fDisplayFrame->GetDisplayFrame()->Connect("ProcessedEvent(Event_t*)", "AliDisplay2", this,"HandleMouseWheel(Event_t*)");
	

	fLeftFrame->Layout();

	fSubFrame->AddFrame(fLeftFrame, new TGLayoutHints( kLHintsBottom | kLHintsLeft | kLHintsExpandY, 5, 5, 2, 2));
	fSubFrame->AddFrame(fRightFrame, new TGLayoutHints( kLHintsBottom | kLHintsRight | kLHintsExpandX | kLHintsExpandY, 5, 5, 2, 2));
	
	Int_t parts[] = {45,45,10};
	fStatusBar = new TGStatusBar(fMainFrame,50,10,kHorizontalFrame);
	fStatusBar->SetParts(parts,3);
	fStatusBar->SetText("AliDisplay v2.0",0);
	fMainFrame->AddFrame(fStatusBar,new TGLayoutHints(kLHintsBottom | kLHintsExpandX,0,0,0,0));
	
	fMainFrame->AddFrame(fSubFrame,new TGLayoutHints( kLHintsBottom | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 2, 2));
	fMainFrame->SetWindowName("Ali Display");
	
	fMainFrame->MapSubwindows();
	fMainFrame->MapWindow();
	LoadSettings();

	
	fMainFrame->Resize(w-10,h);
	fMainFrame->Resize(w,h);
	fMainFrame->SetWMSizeHints(500,500,1280, 1200,1,1);
}

//_____________________________________________________________
AliDisplay2::~AliDisplay2(void)
{
  delete fModules;
  delete [] fEnabledModules;
  delete fModuleInfo;

	delete fSliderFrameLayout;
	delete fSliderFrame;	
	delete fDisplayFrameLayout;
	delete fDisplayFrame;
	//	delete fZoomFrameLayout;
	//	delete fZoomFrame;
	delete fShutterFrameLayout;
	delete fShutterFrame;
	delete fInfoFrameLayout;
	delete fInfoFrame;
	delete fDetectorFrameLayout;
	delete fDetectorFrame;
	
	delete fSubFrame;
	delete fLeftFrame;
	delete fRightFrame;
	delete fMainFrame;
	delete fAliDisplay2rc;

	delete fMenu;
	delete fStatusBar;
}

//_____________________________________________________________
void AliDisplay2::CloseWindow(void)
{
	delete this;
}

//_____________________________________________________________
void AliDisplay2::LoadFromRC()
{
  TEnv *rc=new TEnv(".alidisplayrc");
  SetSliderUpdate(rc->GetValue("AliDisplay.SliderUpdate",kFALSE));
  SetZoomStep(rc->GetValue("AliDisplay.ZoomStep",1.2));
  SetSliderStep(rc->GetValue("AliDisplay.SliderStep",0.01));
  char c[128];
  fRawDataPath = new char[128];
  strcpy(c,gSystem->Getenv("ALICE_ROOT"));
  sprintf(fRawDataPath,"%s%s",c,rc->GetValue("AliDisplay.RawDataPath","/raw"));
  printf("\nRaw data path %s",fRawDataPath);
}

//_____________________________________________________________
void AliDisplay2::SaveToRC()
{
  TEnv *rc=new TEnv(".alidisplayrc");
  rc->SetValue("AliDisplay.SliderUpdate",GetSliderUpdate());
  rc->SetValue("AliDisplay.ZoomStep",GetZoomStep());
  rc->SetValue("AliDisplay.SliderStep",GetSliderStep());
  rc->SetValue("AliDisplay.RawDataPath","/raw");
  rc->SaveLevel(kEnvLocal);
  rc->Save();
}

//_____________________________________________________________
void AliDisplay2::DoSaveSettings(void)
{
  fSliderFrame->SaveToRC();
  SaveToRC();
}

//_____________________________________________________________
void AliDisplay2::LoadSettings()
{
  LoadFromRC();
}

//_____________________________________________________________
void AliDisplay2::Draw(Option_t */*options*/)
{
	fDisplayFrame->DoView(fCurrentView);
}

//_____________________________________________________________
void AliDisplay2::DrawX3d()
{
	fDisplayFrame->DrawX3d();
}

//_____________________________________________________________
void AliDisplay2::DrawGL()
{
	fDisplayFrame->DrawGL();
}

//_____________________________________________________________
void AliDisplay2::ShowNextEvent(Int_t delta)
{
	//Load the next event
	clock_t t1,t2;
	t1=clock();
	Int_t newEvent=0;
	if(delta!=0){
		gAlice->Clear();
		//Int_t currentEvent = gAlice->GetHeader()->GetEvent();
		newEvent = fEventNumber + delta;
		if( newEvent < 0) return;
		gAlice->GetEvent(newEvent);
		fEventNumber += delta;
		//		if(!gAlice->TreeH()) return;
	}
	if(IsEnabled(kHits)) LoadHits();
	if(IsEnabled(kClusters)) LoadClusters(newEvent);
	if(IsEnabled(kHLT)) LoadHLTClusters(newEvent);
	t2=clock();
	//	printf("\nEvent loaded in....%f sec", ((double)t2-t1)/(10000*CLK_TCK));
	Update(kmMODULES);
}

//_____________________________________________________________
void AliDisplay2::FindModules()
{
	//Find the modules used for the simulation and assign these modules to the array fModules
	fModules = new TObjArray;
	TObjArray *modules = gAlice->Modules();
	AliModule *mod;
	Int_t nbm = 0;
	for(Int_t i=0;i<modules->GetEntriesFast();i++){
		mod = (AliModule *) modules->At(i);
		if(!mod) continue;
		const char *avoid = strstr("BODY MAG ABSO DIPO HALL FRAME SHIL PIPE",mod->GetName());
		if(avoid) continue;
		fModules->AddLast(mod);
		nbm++;
	}
	fEnabledModules = new Bool_t[nbm];
	fNbModules = nbm;
	fModuleInfo = new AliModuleInfo(nbm);
	for(Int_t j=0;j<fModules->GetEntriesFast();j++){
	  fModuleInfo->Add(fModules->At(j)->GetName(),j);
	  fEnabledModules[j]=kTRUE;
	}
}

//_____________________________________________________________
void AliDisplay2::LoadHits()
{
	//Load the detected hits from each detector to memory
   gAlice->ResetPoints();
   TIter next(gAlice->Modules());
   AliModule *module;
   Int_t ntracks = gAlice->GetMCApp()->GetNtrack();
   while((module = (AliModule*)next())) 
    {
     AliDetector* detector = dynamic_cast<AliDetector*>(module);
     if(detector) detector->SetTreeAddress();
    }
   next.Reset();
   for (Int_t track=0; track<ntracks;track++) {
      gAlice->ResetHits();
      while((module = (AliModule*)next())) {
         AliDetector* detector = dynamic_cast<AliDetector*>(module);
         if(detector)
           {
             detector->TreeH()->GetEvent(track);
             detector->LoadPoints(track);
           }
      }
      next.Reset();
     }
   fHitsLoaded = kTRUE;
}

//_____________________________________________________________
void AliDisplay2::LoadClusters(Int_t nevent)
{
  //clock_t t1,t2;
	fDisplayFrame->LoadClusters(nevent);
	fClustersLoaded = kTRUE;
	//	printf("\nClusters loaded in....%f sec", ((double)t2-t1)/(10000*CLK_TCK));
}

//_____________________________________________________________
void AliDisplay2::LoadHLTClusters(Int_t nevent)
{
  fDisplayFrame->LoadHLTClusters(nevent);
  fHLTLoaded = kTRUE;
}

//_____________________________________________________________
void AliDisplay2::Enable(Int_t m)
{
  if(m==kHits){
    if((fMode&kHits)==kHits) return;
    fMode = kHits|fMode;
    if(!fHitsLoaded) LoadHits();
    Update(kmPOINTS);
  }
  if(m==kClusters){
    if((fMode&kClusters)==kClusters) return;
    fMode = kClusters|fMode;
    if(!fClustersLoaded) LoadClusters(fEventNumber);
    Update();
  }
  if(m==kHLT){
    if((fMode&kHLT)==kHLT) return;
    fMode = kHLT|fMode;
    if(!fHLTLoaded) {
      LoadHLTClusters(fEventNumber);
    }
    Update();
  }
  if(m==kTracks){
    if((fMode&kTracks)==kTracks) return;
    fMode = kTracks|fMode;
    Update();
  }
};    

//_____________________________________________________________
void AliDisplay2::Disable(Int_t m)
{
  if(m==kHits){
    fMode = fMode|kHits;
    fMode = fMode^kHits;
  }
  if(m==kClusters){
    fMode = fMode|kClusters;
    fMode = fMode^kClusters;
  }
  if(m==kHLT){
    fMode = fMode|kHLT;
    fMode = fMode^kHLT;
  }
  if(m==kTracks){
    fMode = fMode|kTracks;
    fMode = fMode^kTracks;
  }
  Update();
}

//_____________________________________________________________
Bool_t AliDisplay2::IsEnabled(Int_t m)
{
  if(m==kHits){
    if((fMode&kHits)==kHits) return kTRUE;
    return kFALSE;
  }
  if(m==kClusters){
    if((fMode&kClusters)==kClusters) return kTRUE;
    return kFALSE;
  }
  if(m==kHLT){
    if((fMode&kHLT)==kHLT) return kTRUE;
    return kFALSE;
  }
  if(m==kTracks){
    if((fMode&kTracks)==kTracks) return kTRUE;
    return kFALSE;
  }
  return kFALSE;
}

//_____________________________________________________________
void AliDisplay2::HandleMouseWheel(Event_t *event)
{
	if(event->fType != kButtonPress && event->fType != kButtonRelease) return;

	if(event->fCode == kButton4){
		fZoomFactor *= fZoomStep;
		Draw();
	}
	
	if(event->fCode == kButton5){
		fZoomFactor /= fZoomStep;
		Draw();
	}
}

//_____________________________________________________________
void AliDisplay2::HandleResize(Event_t *event)
{
	switch(event->fType){
	case kConfigureNotify:{
		Draw();
	}
	break;
	default:break;
	}
}	

//_____________________________________________________________
void AliDisplay2::Update(Int_t tag)
{
	if(tag==kmMODULES){
		LoadEnabledModules();
		if(((fMode)&kHits)==kHits){
			LoadEnabledHits();
			ApplyCuts();
		}
	}
	if(tag==kmCUTS){
		if(((fMode)&kHits)==kHits)ApplyCuts();
	}
	if(tag==kmPOINTS){
		if(((fMode)&kHits)==kHits){
			LoadEnabledHits();
			ApplyCuts();
		}
	}
	Draw();
	fInfoFrame->Update();
}


