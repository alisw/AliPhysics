
/*
#include "alles.h"
extern AliTPCParam *gTPCParam;
extern AliTPCClustersArray * gCalcClusters;
extern AliTPCClustersArray * gDifClusters;
AliTPCDigitsArray * GetDigitsArray();
  AliTPCClustersArray * GetExactClustersArray();
AliTPCClustersArray *  GetCalcClustersArray(Bool_t newtree=kFALSE);
AliTPCClustersArray *  GetDifClustersArray(Bool_t newtree=kFALSE);
TClonesArray * CompareClusters(TClonesArray *calcclusters,TClonesArray *exactclusters, 
				TClonesArray *diffclusters, Float_t shx, Float_t shy);
*/

void   TPCDigits2Clusters(AliTPCDigitsArray * digarr, AliTPCClustersArray *calcarr,
			  AliTPCClustersArray *exarr=0, AliTPCClustersArray *difarr=0)
{
  //
  //calculate clusters position from digits information
  //

  
  //initialisation of cluster finder
  AliClusterFinder cf;
  cf.SetThreshold(gTPCParam->GetZeroSup());
  cf.SetDetectorParam(gTPCParam);
  cf.SetBFit(kFALSE);
  cf.GetMinuit()->SetPrintLevel(-1);
  cf.GetMinuit()->SetMaxIterations(20);
  
  if (digarr ==0) digarr = GetDigitsArray();
  if ( (digarr == 0) || (digarr->GetTree()==0)) return;  
  if (exarr==0) exarr =GetExactClustersArray();
  if ((exarr!=0)&&(difarr==0)) difarr= GetDifClustersArray();

  if (calcarr==0) calcarr = GetCalcClustersArray(kTRUE);
  //loop over all writen digits segments
  Int_t nsegment = (Int_t)digarr->GetTree()->GetEntries();
  Int_t segment;
  for (segment = 0; segment<nsegment;segment++){
    //load segment with index (TTree internal)  number segment
    AliSimDigits *digrow= (AliSimDigits*)digarr->LoadEntry(segment);
    Int_t sector,row;
    ((AliTPCParam*)digarr->GetParam())->AdjustSectorRow(digrow->GetID(),sector,row);
    
    AliTPCClustersRow * clrow = calcarr->CreateRow(sector,row);
    AliTPCClustersRow * difrow= difarr->CreateRow(sector,row);    
    AliTPCClustersRow * exrow =  exarr->GetRow(sector,row);
    if (exrow==0) exrow = exarr->LoadRow(sector,row);
    TH2F *his = (TH2F*)digrow->GenerHisto();
    digrow->ExpandTrackBuffer();
    if (his==0) return;

    //set current index for cluster finder
    Int_t  index[3]= {0,sector,row};
    cf.GetHisto(his);
    cf.SetDetectorIndex(index);    
    cf.FindPeaks3(clrow->GetArray());

    //get cluster tracks ID
    Int_t ncl = clrow->GetArray()->GetEntriesFast();
    for (Int_t icl = 0 ;icl<ncl; icl++)
      {
	AliDigitCluster *  cl =(AliDigitCluster*)clrow->GetArray()->UncheckedAt(icl);
	Int_t i = (Int_t)cl->fMaxX;
	Int_t j = (Int_t)cl->fMaxY;		
	cl->fTracks[0]=digrow->GetTrackIDFast(i,j,0);
	cl->fTracks[1]=digrow->GetTrackIDFast(i,j,1);
	cl->fTracks[2]=digrow->GetTrackIDFast(i,j,2);	
      }
    Float_t shx = gTPCParam->GetZOffset()/gTPCParam->GetZWidth()+0.5;
    Float_t shy= 0.5;
    CompareClusters(clrow->GetArray(),exrow->GetArray(),difrow->GetArray(),shx,shy);        
    calcarr->StoreRow(sector,row);
    difarr->StoreRow(sector,row);
    digarr->ClearRow(sector,row);
  }    
  //STORING RESULTS
 char treeName[100];
 sprintf(treeName,"TreeCCalc_%s",gTPCParam->GetTitle());
 calcarr->GetTree()->Write(treeName);
 if (difarr!=0){
   char treeName[100];
   sprintf(treeName,"TreeCDif_%s",gTPCParam->GetTitle());
   difarr->GetTree()->Write(treeName);
 }

}


AliTPCClustersArray *  GetCalcClustersArray(Bool_t newtree=kFALSE, const char* name=0) 
{
  //
  //construct AliTPCClusters Array object
  //
  AliTPCClustersArray * arr=0;
  //  if ( (gAlice!=0) && (gAlice->GetDetector("TPC")!=0) ) {    
  //  arr = ((AliTPC*)gAlice->GetDetector("TPC"))->GetClustersArray();
  //  if (arr!=0) arr->Update();
  //}
  if (arr==0) {    
    arr = new AliTPCClustersArray;
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
    if (file==0) file = new TFile("galice.root","update");   
    arr->SetClusterType("AliDigitCluster");         
    arr->Setup(gTPCParam);
    cout<<"Update status : "<<arr->Update()<<"\n";
    char treeName[100];
    if (name==0)
      sprintf(treeName,"TreeCExact_%s",gTPCParam->GetTitle());
    else
      sprintf(treeName,"TreeCCalc_%s",gTPCParam->GetTitle());    
    if (newtree!=kTRUE){
      cout<<"Connect tree status : "<<arr->ConnectTree(treeName)<<"\n";
    }
    else {
      arr->MakeTree();
    }
  }
  gCalcClusters = arr;
  return arr;
}


AliTPCClustersArray *  GetDifClustersArray(Bool_t newtree=kFALSE, const char* name=0) 
{
  //
  //construct AliTPCClusters Array object
  //
  AliTPCClustersArray * arr=0;
  //  if ( (gAlice!=0) && (gAlice->GetDetector("TPC")!=0) ) {    
  //  arr = ((AliTPC*)gAlice->GetDetector("TPC"))->GetClustersArray();
  //  if (arr!=0) arr->Update();
  //}
  if (arr==0) {    
    arr = new AliTPCClustersArray;
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
    if (file==0) file = new TFile("galice.root","update");   
    arr->SetClusterType("AliDifCluster");         
    arr->Setup(gTPCParam);
    cout<<"Update status : "<<arr->Update()<<"\n";
    char treeName[100];
    if (name==0) 
      sprintf(treeName,"TreeCExact_%s",gTPCParam->GetTitle());
    else
      sprintf(treeName,"TreeCDif_%s",gTPCParam->GetTitle());
    if (newtree!=kTRUE){
      cout<<"Connect tree status : "<<arr->ConnectTree(treeName)<<"\n";
    }
    else {
      arr->MakeTree();
    }
  }
  gDifClusters  = arr; 
  return arr;
}



TClonesArray * CompareClusters(TClonesArray *calcclusters,TClonesArray *exactclusters, 
				TClonesArray *diffclusters, Float_t shx, Float_t shy)
{
  if (calcclusters==0) return 0;
  if (exactclusters==0) return 0;
  if (diffclusters==0) return 0;
  TClonesArray & diff=*diffclusters;
   
   
  Int_t nclusters2 = calcclusters->GetEntriesFast();
  Int_t nclusters1 = exactclusters->GetEntriesFast();
  
  for(Int_t i=0; i<nclusters1; i++){
    AliCluster * cl1 = (AliCluster*)exactclusters->UncheckedAt(i);
    AliDigitCluster *cl2;
    AliDifCluster diffc;
    Int_t index=-1;
    Float_t dx=10000;
    Float_t dy=10000;
    for (Int_t j=0; j<nclusters2;j++){
      cl2 = (AliDigitCluster*)calcclusters->UncheckedAt(j);
      if ( (cl2!=0)&& (cl1!=0)){ 
	Float_t ddx = (cl2->fX-shx)-cl1->fX;
	Float_t ddy = (cl2->fY-shy)-cl1->fY;
	if ((ddx*ddx+ddy*ddy)<(dx*dx+dy*dy)){
	  dx=ddx;
	  dy=ddy;
	  index=j;
	}      
      }
    }
    
    cl2 = (AliDigitCluster*)calcclusters->UncheckedAt(index);
    if (cl2!=0){
      diffc.fDx      =dx;
      diffc.fDy      =dy;
      
      diffc.fAngleX  = cl1->fSigmaX2;
      diffc.fAngleY  = cl1->fSigmaY2; 
      diffc.fTracks[0] = cl2->fTracks[0];
      diffc.fTracks[1] = cl2->fTracks[1];
      diffc.fTracks[2] = cl2->fTracks[2];
      diffc.fGenerTrack= cl1->fTracks[0];

      diffc.fX       = cl2->fX;
      diffc.fY       = cl2->fY;
      diffc.fNx      = cl2->fNx;
      diffc.fNy      = cl2->fNy;
      diffc.fQ       = cl2->fQ;    
      diffc.fOrigQ   = cl1->fQ;
      diffc.fSigmaX2 = cl2->fSigmaX2;
      diffc.fSigmaY2 = cl2->fSigmaY2;
      diffc.fArea    = cl2->fArea;
      diffc.fMax     = cl2->fMax;
    }
    else cout<<"pici[ici/n";
    
    new(diff[i]) AliDifCluster(diffc);
  }
  return diffclusters;
} 
