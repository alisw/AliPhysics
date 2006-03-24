#if !defined(__CINT__) || defined(__MAKECINT__)
#include "alles.h"
#include "AliTPCtracker.h"
#include "TView.h"
#include "TPolyMarker3D.h"
#include "AliSimDigits.h"
#include "AliTPCParam.h"
#include "AliRunLoader.h"
#include "AliTPCclusterMI.h"
#endif

/*
  Author:   marian.ivanov@cern.ch

  How to use ?
  .L AliTPCDisplayClustersMI.C+
  AliTPCDisplayClusters disp;
  disp.Init(0,0);   //specify event number, and threshold for the noise
  disp.DisplayClusters();  

*/




class AliTPCDisplayClusters{
public:
  AliTPCDisplayClusters();
  void SetIO(Int_t event);
  void LoadClusters(Int_t noiseth);  
  void DisplayClusters(Int_t first=0, Int_t last=-1);
  void Init(Int_t event, Int_t noiseth){SetIO(event); LoadClusters(noiseth);}
  TObjArray   * fArray;
  AliTPCParam * fParam;
  TTree       * fTree;
  TGeometry   * fGeom;
};


//----------------------------------------------------------------------
AliTPCDisplayClusters::AliTPCDisplayClusters()
{
  fArray = 0;
  fParam = 0;
  fTree  = 0;
  fGeom  = 0;
}
//----------------------------------------------------------------------
void AliTPCDisplayClusters::SetIO(Int_t event)
{
  AliRunLoader* rl = AliRunLoader::Open();
  rl->GetEvent(event);  
  AliLoader* tpcl = (AliLoader*)rl->GetLoader("TPCLoader");
  if (tpcl == 0x0)
    {
      cerr<<"Can not get TPC Loader"<<endl;
      return;
    }  
  rl->CdGAFile();
  fParam=(AliTPCParam *)gDirectory->Get("75x40_100x60_150x60");
  fGeom=(TGeometry*)gDirectory->Get("AliceGeom");
  if(!fParam){
    fParam = new AliTPCParamSR();
    fParam->Update();
  }
  
  if (!fParam) {cerr<<"TPC parameters have not been found !\n"; return ;}  
  tpcl->LoadRecPoints();
  fTree = tpcl->TreeR();
   
}
//----------------------------------------------------------------------
void AliTPCDisplayClusters::LoadClusters(Int_t noiseth)
{
  //
  // load all clusters to memory
  if (fArray) {
    fArray->Delete();
    delete fArray;
  }
  fArray = new TObjArray(fParam->GetNSegmentsTotal());
  //
  AliTPCClustersRow * pclrow = new AliTPCClustersRow;
  pclrow->SetClass("AliTPCclusterMI");
  pclrow->SetArray(0);
  TBranch *br = fTree->GetBranch("Segment");
  br->SetAddress(&pclrow);
  //
  //
  Int_t nrows=Int_t(fTree->GetEntries());
  for (Int_t n=0; n<nrows; n++) {
    //
    pclrow = new AliTPCClustersRow;
    pclrow->SetClass("AliTPCclusterMI");
    pclrow->SetArray(0);
    br->SetAddress(&pclrow);
    //
    br->GetEntry(n);
    AliTPCClustersRow &clrow = *pclrow;
    Int_t ncl=clrow.GetArray()->GetEntriesFast();
    TObjArray * arrrow = new TObjArray(0);
    fArray->AddAt(arrrow,pclrow->GetID());

    // printf("segment\t%d\trow\t%d\tclusters\t%d",n,pclrow->GetID(),ncl);
    while (ncl--) {
      AliTPCclusterMI *cl=(AliTPCclusterMI*)clrow[ncl];
      if (cl->GetQ()>noiseth){
	arrrow->AddLast(cl);    
      }  
    }   
    //    printf("over\t%d\n",arrrow->GetEntries());

  }
}


void AliTPCDisplayClusters::DisplayClusters(Int_t first, Int_t last)
{
  Int_t nrows = fParam->GetNSegmentsTotal();
  // some "constants"
  Int_t markerColorSignal = 5;
  Int_t markerColorBgr = 2;
  Int_t MASK = 10000000;

  TCanvas *c1=new TCanvas("cdisplay", "Cluster display",0,0,700,730);
  TView *v=new TView(1);
  v->SetRange(-330,-360,-330,360,360,1710);
  c1->Clear();
  c1->SetFillColor(1);
  c1->SetTheta(90.);
  c1->SetPhi(0.);
  
  for (Int_t irow=0; irow<nrows; irow++) {

    TObjArray * arr = (TObjArray*)fArray->At(irow);
    if (!arr) continue;
    
    Int_t sec,row;
    fParam->AdjustSectorRow(irow,sec,row);    
    Int_t ncl=arr->GetEntriesFast();

    TPolyMarker3D *pm=new TPolyMarker3D(ncl);
    TPolyMarker3D *pmSignal=new TPolyMarker3D(ncl); // polymarker for signal
    Int_t imarBgr=0;
    Int_t imarSignal=0;
    
    while (ncl--) {
      AliTPCclusterMI *cl=(AliTPCclusterMI*)arr->At(ncl);
      if (cl){
      //
	Double_t x=fParam->GetPadRowRadii(sec,row), y=cl->GetY(), z=cl->GetZ();
	Float_t cs, sn, tmp;
	fParam->AdjustCosSin(sec,cs,sn);
	tmp = x*cs-y*sn; y= x*sn+y*cs; x=tmp; 
	Int_t trackId = cl->GetLabel(0);
	if ( (last>0) &&trackId>last) continue;
	if (trackId<first) continue;
	if (trackId < MASK-1) {
	  pmSignal->SetPoint(imarSignal,x,y,z);
	  imarSignal++;
	} else {
	  pm->SetPoint(imarBgr,x,y,z);
	  imarBgr++;
	}          
      }
    }
        
    // change color for signal
    pm->SetMarkerSize(1); 
    pm->SetMarkerColor(markerColorBgr);
    pm->SetMarkerStyle(1);
    pm->Draw();
    
    pmSignal->SetMarkerSize(1); 
    pmSignal->SetMarkerColor(markerColorSignal);
    pmSignal->SetMarkerStyle(1);
    pmSignal->Draw();      
  }
  
  
  TNode * main = (TNode*)((fGeom->GetListOfNodes())->First());
  TIter next(main->GetListOfNodes());
  TNode  *module=0;
  while((module = (TNode*)next())) {
    char ch[100];
    sprintf(ch,"%s\n",module->GetTitle());
    //printf("%s\n",module->GetTitle());
    if (ch[0]=='T'&&ch[1]=='P' && ch[2]=='C')  //if TPC draw
      module->SetVisibility(3);
    else
      module->SetVisibility(-1);
  }
  
  
  fGeom->Draw("same");
  c1->Modified(); c1->Update(); 
  
  

}
