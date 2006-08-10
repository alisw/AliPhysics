#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>

#include "AliTRDtracker.h"
#include "AliTRDclusterMI.h"
#include "AliTRDhit.h"
#include "AliTRDv1.h"
#include "AliTRDgeometry.h"
#include "AliTRDgeometryDetail.h"
#include "AliTRDparameter.h"
#include "AliTRDclusterCorrection.h"
#include "alles.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "Rtypes.h"
#include "TTree.h"

#include "AliRunLoader.h"
#include "AliStack.h"
#include "TF1.h"
#include "AliTrackReference.h"
#endif    
#include "AliTRDclusterErrors.h"


ClassImp(AliTRDCI)
ClassImp(AliTRDExactPoint)
ClassImp(AliTRDClusterErrAnal)
ClassImp(AliTRDClusterErrDraw)



AliTRDClusterErrAnal::AliTRDClusterErrAnal(Char_t *chloader )
{
  //
  //SET Input loaders
  if (gAlice){
    delete gAlice->GetRunLoader();
    delete gAlice;
    gAlice = 0x0;
  }    
  fRunLoader = AliRunLoader::Open(chloader);
  if (fRunLoader == 0x0){
    cerr<<"Can not open session"<<endl;
    return;
  }
  fTRDLoader = fRunLoader->GetLoader("TRDLoader");  
  if (fTRDLoader == 0x0){
    cerr<<"Can not get TRD Loader"<<endl;
    return ;
  }   
  if (fRunLoader->LoadgAlice()){
    cerr<<"Error occured while l"<<endl;
    return;
  }
  fRunLoader->CdGAFile();
  fTracker = new AliTRDtracker(gFile);
  fParam = (AliTRDparameter*) gFile->Get("TRDparameter");
  fGeometry = new AliTRDgeometryDetail();   
  fTRD      = (AliTRDv1*) gAlice->GetDetector("TRD");
  //
  AliTRDCI * clinfo = new AliTRDCI();
  fFileA  = new TFile("trdclusteranal.root","recreate");
  fFileA->cd();
  fTreeA  = new TTree("trdcl","trdcl");
  fTreeA->Branch("trdcl","AliTRDCI",&clinfo);
  delete clinfo;
}

void AliTRDClusterErrAnal::SetIO(Int_t event)
{
  //
  //set input output for given event
  fRunLoader->SetEventNumber(event);
  fRunLoader->LoadHeader();
  fRunLoader->LoadKinematics();
  fRunLoader->LoadTrackRefs();
  fTRDLoader->LoadHits();
  fTRDLoader->LoadRecPoints("read");
  //
  fStack = fRunLoader->Stack();
  fHitTree = fTRDLoader->TreeH();
  fClusterTree = fTRDLoader->TreeR();
  fReferenceTree = fRunLoader->TreeTR();
  fTracker->LoadClusters(fClusterTree);
  //
}

void AliTRDClusterErrAnal::LoadClusters()
{
  //
  //
  // Load clusters  
  TObjArray *ClusterArray = new TObjArray(400);   
  TObjArray carray(2000);
  TBranch *branch=fClusterTree->GetBranch("TRDcluster");
  if (!branch) {
    Error("ReadClusters","Can't get the branch !");
    return;
  }
  Int_t over5 =0;
  Int_t over10=0;

  branch->SetAddress(&ClusterArray);
  Int_t nentries = (Int_t)fClusterTree->GetEntries();
  for (Int_t i=0;i<nentries;i++){
    fClusterTree->GetEvent(i);
    Int_t nCluster = ClusterArray->GetEntriesFast();
    for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) { 
      AliTRDcluster* c = (AliTRDcluster*)ClusterArray->UncheckedAt(iCluster);
      carray.AddLast(c);
      ClusterArray->RemoveAt(iCluster);
      if (c->GetQ()>5)  over5++;
      if (c->GetQ()>10) over10++;
    }
  }
  Int_t nClusters = carray.GetEntriesFast();
  printf("Total number of clusters %d\t%d\t%d\n", nClusters,over5,over10);
  //
  //
  //SORT clusters
  //
  Int_t all=0;  
  for (Int_t i=0;i<nClusters;i++){
    AliTRDcluster *cl = (AliTRDcluster *) carray.UncheckedAt(i);  
    Int_t plane = fGeometry->GetPlane(cl->GetDetector());
    if (plane>=12) continue;
    Int_t time = cl->GetLocalTimeBin(); 
    if (time>=100) continue;
    Int_t sector = fGeometry->GetSector(cl->GetDetector());
    if (sector>=18){
      printf("problem1\n");
      continue;
    }    
    fClusters[plane][time][sector].AddLast(cl);
    all++;
  }
}

void AliTRDClusterErrAnal::SortReferences()
{
  //
  //
  //
  printf("Sorting references\n");
  fReferences = new TObjArray;
  Int_t ntracks = fStack->GetNtrack();
  fReferences->Expand(ntracks);
  Int_t nentries = (Int_t)fReferenceTree->GetEntries();
  TClonesArray * arr = new TClonesArray("AliTrackReference");
  TBranch * br = fReferenceTree->GetBranch("TRD");
  br->SetAddress(&arr);
  //
  TClonesArray *labarr=0;
  Int_t nreferences=0;
  Int_t nreftracks=0;
  for (Int_t iprim=0;iprim<nentries;iprim++){
    if (br->GetEntry(iprim)){
      for (Int_t iref=0;iref<arr->GetEntriesFast();iref++){
	AliTrackReference *ref =(AliTrackReference*)arr->At(iref);
	if (!ref) continue;
	Int_t lab = ref->GetTrack();
	if ( (lab<0) || (lab>ntracks)) continue;
	//
	if (fReferences->At(lab)==0) {
	  labarr = new TClonesArray("AliTrackReference"); 
	  fReferences->AddAt(labarr,lab);
	  nreftracks++;
	}
	TClonesArray &larr = *labarr;
	new(larr[larr.GetEntriesFast()]) AliTrackReference(*ref);
	nreferences++;
      }
    }
  }
  printf("Total number of references = \t%d\n", nreferences);
  printf("Total number of tracks with references = \t%d\n", nreftracks);
  printf("End - Sorting references\n");
  
}

AliTrackReference * AliTRDClusterErrAnal::FindNearestReference(Int_t lab, Float_t pos[3], Float_t dmax)
{
  //
  //
  //
  if (fReferences->At(lab)==0) return 0;
  AliTrackReference *nearest=0;
  TClonesArray * arr = (TClonesArray *)fReferences->At(lab);
  for (Int_t iref =0;iref<arr->GetEntriesFast();iref++){
    AliTrackReference * ref = ( AliTrackReference *)arr->UncheckedAt(iref);
    if (!ref) continue;
    Float_t delta = (pos[0]-ref->X())*(pos[0]-ref->X());
    delta += (pos[1]-ref->Y())*(pos[1]-ref->Y());
    delta += (pos[2]-ref->Z())*(pos[2]-ref->Z());
    delta = TMath::Sqrt(delta);
    if (delta<dmax){
      dmax=delta;
      nearest = ref;
    }
  }
  return nearest;
}

void AliTRDClusterErrAnal::MakeExactPoints(Int_t trackmax)
{
  //
  //make exact points:-)	
  //
  // 

  fExactPoints.Delete();
  fExactPoints.Expand(fStack->GetNtrack());
  //
  Double_t fSum=0;
  Double_t fSumQ =0;
  Double_t fSumX=0;
  Double_t fSumX2=0;
  Double_t fSumXY=0;
  Double_t fSumXZ=0;
  Double_t fSumY=0;
  Double_t fSumZ=0;
  //
  Int_t entries = Int_t(fHitTree->GetEntries());
  printf("Number of primary  entries\t%d\n",entries);
  entries = TMath::Min(trackmax,entries);
  Int_t nallpoints = 0;

  Int_t nalltracks =0;
  Int_t pointspertrack =0;

  for (Int_t entry=0;entry<entries; entry++){
    gAlice->ResetHits();
    fHitTree->GetEvent(entry);
    Int_t lastlabel    = -1;
    Int_t lastdetector = -1;
    Int_t lasttimebin   = -1;
    Float_t lastpos[3];
    //
    for(AliTRDhit *hit = (AliTRDhit *) fTRD->FirstHit(-1); hit; 
	hit = (AliTRDhit *) fTRD->NextHit()) {
      //
      Int_t label    = hit->Track();
      TParticle * particle = fStack->Particle(label);
      if (!particle) continue;
      if (particle->Pt()<0.05) continue;
      Int_t detector = hit->GetDetector();
      Int_t plane    = fGeometry->GetPlane(detector);
      //
      //
      if (hit->GetCharge()==0) continue;
      Float_t pos[3] = {hit->X(),hit->Y(),hit->Z()};
      Int_t indexes[3];
      fGeometry->Global2Detector(pos,indexes,fParam);
      //
      Float_t rot[3];
      fGeometry->Rotate(detector,pos,rot);
      //rot[0] *=-1;
      //  rot[1] *=-1;
      //
      //
      Float_t  time0    = fParam->GetTime0(plane);
      Int_t    timebin  = Int_t(TMath::Nint(((time0 - rot[0])/fParam->GetTimeBinSize())+ fParam->GetTimeBefore())+0.1); 
      if (timebin<0) continue;
      //
      //
      if (label!=lastlabel || detector != lastdetector || lasttimebin !=timebin){
	//
	if (label!=lastlabel){
	  fExactPoints.AddAt(new TClonesArray("AliTRDExactPoint",0),label);
	  //printf("new particle\t%d\n",label);
	  nalltracks++;
	  //	  printf("particle\t%d- hits\t%d\n",lastlabel, pointspertrack);
	  pointspertrack=0;
	}
	
	if ( (fSum>1) && lasttimebin>=0 && lasttimebin<fParam->GetTimeMax() ){
	  //if we have enough info for given layer time bin - store it
	  AliTrackReference * ref = FindNearestReference(lastlabel,lastpos,4.);
	  Float_t rotmom[3];
	  Float_t rotpos[3];
	  Float_t refangle=0;
	  if (ref){
	    Float_t mom[3] = {ref->Px(),ref->Py(),ref->Pz()};
	    Float_t refpos[3] = {ref->X(),ref->Y(),ref->Z()};
	    fGeometry->Rotate(detector,mom,rotmom);
	    fGeometry->Rotate(detector,refpos,rotpos);
	    refangle = rotmom[1]/rotmom[0];

	  }

	  Double_t ay,by,az,bz;
	  Double_t det = fSum*fSumX2-fSumX*fSumX;
	  if (TMath::Abs(det)> 0.000000000000001) { 
	    by = (fSum*fSumXY-fSumX*fSumY)/det;
	    ay = (fSumX2*fSumY-fSumX*fSumXY)/det;
	    
	  }else{
	    ay =fSumXY/fSumX;
	    by =0;	   
	  }
	  if (TMath::Abs(det)> 0.000000000000001) { 
	    bz = (fSum*fSumXZ-fSumX*fSumZ)/det;
	    az = (fSumX2*fSumZ-fSumX*fSumXZ)/det;	  
	  }else{
	    az =fSumXZ/fSumX;
	    bz =0;	   
	  }
	  //
	  Float_t lastplane = fGeometry->GetPlane(lastdetector);
	  Float_t time0    = fParam->GetTime0(lastplane);
	  Float_t xcenter0 = time0 - (lasttimebin - fParam->GetTimeBefore()+0.5)*fParam->GetTimeBinSize();
	  Float_t xcenter = fTracker->GetX(0,lastplane,lasttimebin);
	  if (TMath::Abs(xcenter-xcenter0)>0.001){
	    printf("problem");
	  }
	  
	  Float_t ty = ay + by * xcenter;
	  Float_t tz = az + bz * xcenter;
	  //

	  TClonesArray * arr = (TClonesArray *) fExactPoints.At(label);
	  TClonesArray & larr= *arr;
	  Int_t arrent = arr->GetEntriesFast();
	  AliTRDExactPoint * point = new (larr[arrent]) AliTRDExactPoint;
	  nallpoints++;

	  if (ref){
	    point->SetReference(ref);
	    point->fTRefAngleY = rotmom[1]/rotmom[0];
	  }
	  point->fTX  = xcenter;
	  point->fTY  = ty;
	  point->fTZ  = tz;
	  point->fTAY = by;
	  point->fTAZ = bz;
	  //
	  point->fGx  = lastpos[0];
	  point->fGy  = lastpos[1];
	  point->fGz  = lastpos[2];

	  //
	  point->fDetector     = lastdetector;
	  point->fLocalTimeBin = lasttimebin;
	  point->fPlane        = fGeometry->GetPlane(lastdetector); 
	  point->fSector      = fGeometry->GetSector(lastdetector); 
	  point->fPlaneMI      = indexes[0]; 
	  //
	  point->fTPrim = fSum;
	  point->fTQ    = fSumQ;	    
	  //
	}
	lastdetector = detector;
	lastlabel    = label;
	lasttimebin  = timebin;
	fSum=fSumQ=fSumX=fSumX2=fSumXY=fSumXZ=fSumY=fSumZ=0.;
      }
      //
      lastpos[0] = hit->X();
      lastpos[1] = hit->Y();
      lastpos[2] = hit->Z();
      fSum++;
      fSumQ  +=hit->GetCharge();
      fSumX  +=rot[0];
      fSumX2 +=rot[0]*rot[0];
      fSumXY +=rot[0]*rot[1];      
      fSumXZ +=rot[0]*rot[2];
      fSumY  +=rot[1];      
      fSumZ  +=rot[2];     
      pointspertrack++;
    }
  }
  //
  printf("Found %d exact points\n",nallpoints); 
}






Int_t AliTRDClusterErrAnal::Analyze(Int_t trackmax) {  

  //
  // comparison works with both cluster types MI and old also
  //dummy cluster to be fill if not cluster info
  AliTRDclusterMI clmi;  
  TClass * classmi =  clmi.IsA();
  //
  //SetOutput
  AliTRDCI * clinfo = new AliTRDCI();
  TBranch * clbr = fTreeA->GetBranch("trdcl");
  clbr->SetAddress(&clinfo);

  SetIO(0);
  SortReferences();
  MakeExactPoints(trackmax);
  LoadClusters();
  //
  trackmax =  fStack->GetNtrack();
  //
  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrack = (Int_t)fExactPoints.GetEntries();
  printf("Found %d charged in TRD in first %d particles", nTrack, trackmax);
  //

  for (Int_t itrack = 0; itrack<trackmax; itrack++){
    TClonesArray *arrpoints = (TClonesArray*)fExactPoints.At(itrack);
    
    if (!arrpoints) continue;
    //printf("new particle\t%d\n",itrack);
    TParticle * particle = fStack->Particle(itrack);
    if (!particle) continue;
    //printf("founded in kine tree \t%d\n",itrack);
    Int_t npoints = arrpoints->GetEntriesFast();
    if (npoints<10) continue;
    //printf("have enough points \t%d\t%d\n",itrack,npoints);

    for (Int_t ipoint=0;ipoint<npoints;ipoint++){
      AliTRDExactPoint * point = (AliTRDExactPoint*)arrpoints->UncheckedAt(ipoint);
      if (!point) continue;
      //
      Int_t sec = fGeometry->GetSector(point->fDetector);
      if (sec>18){
	printf("problem2\n");
      }
      TObjArray & cllocal = fClusters[point->fPlane][point->fLocalTimeBin][sec];       
      Int_t nclusters = cllocal.GetEntriesFast();
      Float_t maxdist = 10;
      AliTRDcluster * nearestcluster =0;
      clinfo->fNClusters=0;
      //find nearest cluster to hit with given label
      for (Int_t icluster =0; icluster<nclusters; icluster++){
	AliTRDcluster * cluster = (AliTRDcluster*)cllocal.UncheckedAt(icluster);
	if (!cluster) continue;
	if ( (cluster->GetLabel(0)!= itrack) &&  (cluster->GetLabel(1)!= itrack)&&(cluster->GetLabel(2)!= itrack))
	  continue;
	Float_t dist = TMath::Abs(cluster->GetY()-point->fTY);
	if (TMath::Abs(cluster->GetZ()-point->fTZ)>5.5 || dist>3.) continue; 
	clinfo->fNClusters++;
	if (dist<maxdist){
	  maxdist = dist;
	  nearestcluster = cluster;
	}
      }      
      //
      clinfo->fEp  = *point;
      clinfo->fP   = *particle;
      if (!nearestcluster) {
	clinfo->fStatus=1;
	clinfo->fCl = clmi;
      }
      else{
	clinfo->fStatus=0;
	if (nearestcluster->IsA()==classmi){
	  clinfo->fCl    =*((AliTRDclusterMI*)nearestcluster);    
	}
	else{
	  clinfo->fCl    = *nearestcluster;
	}
	//     
	Float_t dz      = clinfo->fCl.GetZ()-point->fTZ;
	Double_t h01    = sin(TMath::Pi() / 180.0 * fParam->GetTiltingAngle());
	clinfo->fTDistZ = dz;
	clinfo->fDYtilt = h01*dz*((point->fPlane%2)*2.-1.); 
	//
	clinfo->fNTracks =1;
	if (nearestcluster->GetLabel(1)>=0)  clinfo->fNTracks++;
	if (nearestcluster->GetLabel(2)>=0)  clinfo->fNTracks++;  
	clinfo->Update();
      }
      //
      fTreeA->Fill();
    }
  }
  
  
  fFileA->cd();
  fTreeA->Write();
  fFileA->Close();
  return 0;
}

AliTRDclusterCorrection*   AliTRDClusterErrDraw::MakeCorrection(TTree * tree, Float_t offset)
{
  //
  //
  // make corrections
  AliTRDclusterCorrection * cor = new AliTRDclusterCorrection;
  cor->SetOffsetAngle(offset); 
  for (Int_t iplane=0;iplane<6;iplane++)
    for (Int_t itime=0;itime<15;itime++)
      for (Int_t iangle=0; iangle<20;iangle++){
	Float_t angle = cor->GetAngle(iangle);
	TH1F delta("delta","delta",30,-0.3,0.3);
	char selection[100]="fStatus==0&&fNTracks<2";
	char selectionall[1000];
	sprintf(selectionall,"%s&&abs(fEp.fTAY-%f)<0.2&&fEp.fPlane==%d&&fCl.fTimeBin==%d&&fCl.fQ>20",
		selection,angle,iplane,itime);
	printf("\n%s",selectionall);
	tree->Draw("fEp.fTY-fCl.fY+fDYtilt>>delta",selectionall);
	gPad->Update();
	printf("\nplane\t%d\ttime%d\tangle%f",iplane,itime,angle);
	printf("\tentries%f\tmean\t%f\t%f",delta.GetEntries(),delta.GetMean(),delta.GetRMS());	
	cor->SetCorrection(iplane,itime,angle,delta.GetMean(),delta.GetRMS());
      }
  TFile * f = new TFile("TRDcorrection.root","new");
  if (!f) f = new TFile("TRDcorrection.root","recreate");
  f->cd();
  cor->Write("TRDcorrection");
  f->Close();
  return cor; 
}

TH1F * AliTRDClusterErrDraw::ResDyVsAmp(TTree* tree, const char* selection, Float_t t0, Float_t ampmin, Float_t ampmax)
{
  //
  //
  TH2F hisdy("resy","resy",10,ampmin,ampmax,30,-0.3,0.3);
  char expression[1000];
  sprintf(expression,"fEp.fTY-fCl.fY+fDYtilt+%.4f*fEp.fTAY:fCl.fQ>>resy",t0);
  char selectionall[1000];
  sprintf(selectionall,"fStatus==0&&%s",selection);
  printf("%s\n",expression);
  printf("%s\n",selectionall);
  tree->Draw(expression,selectionall);
  return CreateResHisto(&hisdy);
}


TH1F * AliTRDClusterErrDraw::ResDyVsRelPos(TTree* tree, const char* selection, Float_t t0, Float_t min, Float_t max)
{
  //
  //
  min *=128;
  max *=128;
  TH2F hisdy("resy","resy",10,min,max,30,-0.3,0.3);
  char expression[1000];
  sprintf(expression,"fEp.fTY-fCl.fY+fDYtilt+%.4f*fEp.fTAY:fCl.fRelPos>>resy",t0);
  char selectionall[1000];
  sprintf(selectionall,"fStatus==0&&%s",selection);
  printf("%s\n",expression);
  printf("%s\n",selectionall);
  tree->Draw(expression,selectionall);
  return CreateResHisto(&hisdy);
}


TH1F * AliTRDClusterErrDraw::ResDyVsAngleY(TTree* tree, const char* selection, Float_t t0, Float_t min, Float_t max)
{
  //
  //
  TH2F hisdy("resy","resy",10,min,max,30,-0.3,0.3);

  char expression[1000];
  sprintf(expression,"fEp.fTY-fCl.fY+fDYtilt+%f*fEp.fTAY:fEp.fTAY>>resy",t0);
  char selectionall[1000];
  sprintf(selectionall,"fStatus==0&&%s",selection);

  tree->Draw(expression,selectionall);
  return CreateResHisto(&hisdy);
}

void AliTRDClusterErrDraw::AliLabelAxes(TH1* histo, const char* xAxisTitle, const char* yAxisTitle)
{
  histo->GetXaxis()->SetTitle(xAxisTitle);
  histo->GetYaxis()->SetTitle(yAxisTitle);
}


TH1F* AliTRDClusterErrDraw::CreateEffHisto(TH1F* hGen, TH1F* hRec)
{
  Int_t nBins = hGen->GetNbinsX();
  TH1F* hEff = (TH1F*) hGen->Clone("hEff");
  hEff->SetTitle("");
  hEff->SetStats(kFALSE);
  hEff->SetMinimum(0.);
  hEff->SetMaximum(110.);

  for (Int_t iBin = 0; iBin <= nBins; iBin++) {
    Double_t nGen = hGen->GetBinContent(iBin);
    Double_t nRec = hRec->GetBinContent(iBin);
    if (nGen > 0) {
      Double_t eff = nRec/nGen;
      hEff->SetBinContent(iBin, 100. * eff);
      Double_t error = sqrt((eff*(1.-eff)+0.01) / nGen);      
      //      if (error == 0) error = sqrt(0.1/nGen);
      //
      if (error == 0) error = 0.0001;
      hEff->SetBinError(iBin, 100. * error);
    } else {
      hEff->SetBinContent(iBin, 100. * 0.5);
      hEff->SetBinError(iBin, 100. * 0.5);
    }
  }

  return hEff;
}



TH1F* AliTRDClusterErrDraw::CreateResHisto(TH2F* hRes2, Bool_t draw,  Bool_t drawBinFits, 
		     Bool_t overflowBinFits)
{
  TVirtualPad* currentPad = gPad;
  TAxis* axis = hRes2->GetXaxis();
  Int_t nBins = axis->GetNbins();
  TH1F* hRes, *hMean;
  if (axis->GetXbins()->GetSize()){
    hRes = new TH1F("hRes", "", nBins, axis->GetXbins()->GetArray());
    hMean = new TH1F("hMean", "", nBins, axis->GetXbins()->GetArray());
  }
  else{
    hRes = new TH1F("hRes", "", nBins, axis->GetXmin(), axis->GetXmax());
    hMean = new TH1F("hMean", "", nBins, axis->GetXmin(), axis->GetXmax());

  }
  hRes->SetStats(false);
  hRes->SetOption("E");
  hRes->SetMinimum(0.);
  //
  hMean->SetStats(false);
  hMean->SetOption("E");
 
  // create the fit function
  //TKFitGaus* fitFunc = new TKFitGaus("resFunc");
  TF1 * fitFunc = new TF1("G","[3]+[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",-3,3);
  
  fitFunc->SetLineWidth(2);
  fitFunc->SetFillStyle(0);
  // create canvas for fits
  TCanvas* canBinFits = NULL;
  Int_t nPads = (overflowBinFits) ? nBins+2 : nBins;
  Int_t nx = Int_t(sqrt(nPads-1.));// + 1;
  Int_t ny = (nPads-1) / nx + 1;
  if (drawBinFits) {
    canBinFits = (TCanvas*)gROOT->FindObject("canBinFits");
    if (canBinFits) delete canBinFits;
    canBinFits = new TCanvas("canBinFits", "fits of bins", 200, 100, 500, 700);
    canBinFits->Divide(nx, ny);
  }

  // loop over x bins and fit projection
  Int_t dBin = ((overflowBinFits) ? 1 : 0);
  for (Int_t bin = 1-dBin; bin <= nBins+dBin; bin++) {
    if (drawBinFits) canBinFits->cd(bin + dBin);
    TH1D* hBin = hRes2->ProjectionY("hBin", bin, bin);
    //    
    if (hBin->GetEntries() > 10) {
      fitFunc->SetParameters(hBin->GetMaximum(),hBin->GetMean(),hBin->GetRMS(),0.02*hBin->GetMaximum());
      hBin->Fit(fitFunc,"s");
      Double_t sigma = TMath::Abs(fitFunc->GetParameter(2));

      if (sigma > 0.){
	hRes->SetBinContent(bin, TMath::Abs(fitFunc->GetParameter(2)));
	hMean->SetBinContent(bin, fitFunc->GetParameter(1));	
      }
      else{
	hRes->SetBinContent(bin, 0.);
	hMean->SetBinContent(bin,0);
      }
      hRes->SetBinError(bin, fitFunc->GetParError(2));
      hMean->SetBinError(bin, fitFunc->GetParError(1));
      
      //
      //

    } else {
      hRes->SetBinContent(bin, 0.);
      hRes->SetBinError(bin, 0.);
      hMean->SetBinContent(bin, 0.);
      hMean->SetBinError(bin, 0.);
    }
    

    if (drawBinFits) {
      char name[256];
      if (bin == 0) {
	sprintf(name, "%s < %.4g", axis->GetTitle(), axis->GetBinUpEdge(bin));
      } else if (bin == nBins+1) {
	sprintf(name, "%.4g < %s", axis->GetBinLowEdge(bin), axis->GetTitle());
      } else {
	sprintf(name, "%.4g < %s < %.4g", axis->GetBinLowEdge(bin),
		axis->GetTitle(), axis->GetBinUpEdge(bin));
      }
      canBinFits->cd(bin + dBin);
      hBin->SetTitle(name);
      hBin->SetStats(kTRUE);
      hBin->DrawCopy("E");
      canBinFits->Update();
      canBinFits->Modified();
      canBinFits->Update();
    }
    
    delete hBin;
  }

  delete fitFunc;
  currentPad->cd();
  if (draw){
    currentPad->Divide(1,2);
    currentPad->cd(1);
    hRes->Draw();
    currentPad->cd(2);
    hMean->Draw();
  }
  
  return hRes;
}
