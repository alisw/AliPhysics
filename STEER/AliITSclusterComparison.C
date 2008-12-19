#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>

#include "AliITStrackerMI.h"
#include "AliITSclusterV2.h"

#include "AliITShit.h"
#include "AliITS.h"
#include "AliITSgeom.h"
//#include "AliITSgeometryDetail.h"
//#include "AliITSparameter.h"
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
#include "AliITSclusterComparison.h"


ClassImp(AliITSCI)
ClassImp(AliITSClusterErrAnal)



AliITSClusterErrAnal::AliITSClusterErrAnal(Char_t *chloader )
{
  //
  //SET Input loaders
  if (gAlice){
    delete AliRunLoader::GetRunLoader();
    delete gAlice;
    gAlice = 0x0;
  }    
  fRunLoader = AliRunLoader::Open(chloader);
  if (fRunLoader == 0x0){
    cerr<<"Can not open session"<<endl;
    return;
  }
  fITSLoader = fRunLoader->GetLoader("ITSLoader");  
  if (fITSLoader == 0x0){
    cerr<<"Can not get ITS Loader"<<endl;
    return ;
  }   
  if (fRunLoader->LoadgAlice()){
    cerr<<"Error occured while l"<<endl;
    return;
  }
  fRunLoader->CdGAFile();  
  fITS = (AliITS*)gAlice->GetDetector("ITS");
  if (!fITS) {
    cerr<<"AliITSclusterComparison.C : Can not find the ITS detector !"<<endl;
  }
  AliITSgeom *geom = fITS->GetITSgeom();
  //An instance of the ITS tracker
  fTracker = new AliITStrackerMI(geom);
  //
  //
  AliITSCI * clinfo = new AliITSCI();
  fFileA  = new TFile("itsclusteranal.root","recreate");
  fFileA->cd();
  fTreeA  = new TTree("itscl","itscl");
  fTreeA->Branch("itscl","AliITSCI",&clinfo);

  AliITSclusterV2 * cl = new AliITSclusterV2;
  fTreeB  = new TTree("Clusters","Clusters");
  fTreeB->Branch("cl","AliITSclusterV2",&cl);

  fClusters = 0;
  delete clinfo;
}

void AliITSClusterErrAnal::SetIO(Int_t event)
{
  //
  //set input output for given event
  fRunLoader->SetEventNumber(event);
  fRunLoader->LoadHeader();
  fRunLoader->LoadKinematics();
  fRunLoader->LoadTrackRefs();
  fITSLoader->LoadHits();
  fITSLoader->LoadRecPoints("read");
  //
  fStack = fRunLoader->Stack();
  fHitTree = fITSLoader->TreeH();
  fClusterTree = fITSLoader->TreeR();
  fReferenceTree = fRunLoader->TreeTR();
  //
}

void AliITSClusterErrAnal::LoadClusters()
{
  //
  //
  // Load clusters
  if (fClusters) {
    fClusters->Delete();
    delete fClusters;
  }
  Int_t nparticles = fStack->GetNtrack();
  fClusters = new TObjArray(nparticles);
  //
  TObjArray *ClusterArray = new TClonesArray("AliITSclusterV2");   
  TObjArray carray(2000);
  TBranch *branch=fClusterTree->GetBranch("Clusters");
  if (!branch) {
    Error("ReadClusters","Can't get the branch !");
    return;
  }
  branch->SetAddress(&ClusterArray);
  TBranch *brcl = fTreeB->GetBranch("cl");
  AliITSclusterV2* cluster= new AliITSclusterV2;
  brcl->SetAddress(&cluster);
  Int_t nentries = (Int_t)fClusterTree->GetEntries();
  for (Int_t i=0;i<nentries;i++){
    fClusterTree->GetEvent(i);
    Int_t nCluster = ClusterArray->GetEntriesFast();
    SignDeltas(ClusterArray,-3.8);
    for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) { 
      cluster = (AliITSclusterV2*)ClusterArray->UncheckedAt(iCluster);
      carray.AddLast(new AliITSclusterV2(*cluster));
      fTreeB->Fill();
    }
  }
  //
  Int_t nClusters = carray.GetEntriesFast();
  printf("Total number of clusters %d \n", nClusters);
  fTracker->LoadClusters(fClusterTree);
  //
  //
  //SORT clusters
  //
  Int_t all=0;
  Int_t noisecl =0;
  for (Int_t i=0;i<nClusters;i++){
    AliITSclusterV2 *cl = (AliITSclusterV2 *) carray.UncheckedAt(i);  
    if (cl->GetLabel(0)<0) noisecl++;
    //
    for (Int_t itrack=0; itrack<3;itrack++){
      Int_t lab = cl->GetLabel(itrack);
      if (lab>=0){
	TObjArray * array = (TObjArray*)fClusters->At(lab);
	if (array==0){
	  array = new TObjArray(20);
	  fClusters->AddAt(array,lab);
	}
	array->AddLast(cl);
	all++;
      }     	
    }    
  }
  printf("Total number of track clusters %d \n", all);
  printf("Total number of noise clusters %d \n", noisecl);
}

void AliITSClusterErrAnal::SignDeltas( TObjArray *ClusterArray, Float_t vz)
{
  //
  //  
  Int_t entries = ClusterArray->GetEntriesFast();
  if (entries<4) return;
  AliITSclusterV2* cluster = (AliITSclusterV2*)ClusterArray->At(0);
  Int_t layer = cluster->GetLayer();
  if (layer>1) return;
  Int_t index[10000];
  Int_t ncandidates=0;
  Float_t r = (layer>0)? 7:4;
  // 
  for (Int_t i=0;i<entries;i++){
    AliITSclusterV2* cl0 = (AliITSclusterV2*)ClusterArray->At(i);
    Float_t nz = 1+(cl0->GetZ()-vz)/r;
    if (cl0->GetNy()+cl0->GetNz()<=5+2*layer+nz) continue;
    index[ncandidates] = i;  //candidate to belong to delta electron track
    ncandidates++;
    if (cl0->GetNy()+cl0->GetNz()>9+2*layer+nz) {
      cl0->SetDeltaProbability(1);
    }
  }
  //
  //  
  //
  for (Int_t i=0;i<ncandidates;i++){
    AliITSclusterV2* cl0 = (AliITSclusterV2*)ClusterArray->At(index[i]);
    if (cl0->GetDeltaProbability()>0.8) continue;
    // 
    Int_t ncl = 0;
    Float_t y[100],z[100],sumy,sumz,sumy2, sumyz, sumw;
    sumy=sumz=sumy2=sumyz=sumw=0.0;
    for (Int_t j=0;j<ncandidates;j++){
      if (i==j) continue;
      AliITSclusterV2* cl1 = (AliITSclusterV2*)ClusterArray->At(index[j]);
      //
      Float_t dz = cl0->GetZ()-cl1->GetZ();
      Float_t dy = cl0->GetY()-cl1->GetY();
      if (TMath::Sqrt(dz*dz+dy*dy)<0.2){
	Float_t weight = cl1->GetNy()+cl1->GetNz()-2;
	y[ncl] = cl1->GetY();
	z[ncl] = cl1->GetZ();
	sumy+= y[ncl]*weight;
	sumz+= z[ncl]*weight;
	sumy2+=y[ncl]*y[ncl]*weight;
	sumyz+=y[ncl]*z[ncl]*weight;
	sumw+=weight;
	ncl++;
      }
    }
    if (ncl<4) continue;
    Float_t det = sumw*sumy2  - sumy*sumy;
    Float_t delta=1000;
    if (TMath::Abs(det)>0.01){
      Float_t z0  = (sumy2*sumz - sumy*sumyz)/det;
      Float_t k   = (sumyz*sumw - sumy*sumz)/det;
      delta = TMath::Abs(cl0->GetZ()-(z0+k*cl0->GetY()));
    }
    else{
      Float_t z0  = sumyz/sumy;
      delta = TMath::Abs(cl0->GetZ()-z0);
    }
    if ( delta<0.05) {
      cl0->SetDeltaProbability(1-20.*delta);
    }   
  }
}

void AliITSClusterErrAnal::LoadParticles()
{
  //
  //
  //  fReferences = new TObjArray;
}

void AliITSClusterErrAnal::SortReferences()
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
  TBranch * br = fReferenceTree->GetBranch("ITS");
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
	//if ( (lab<0) || (lab>ntracks)) continue;
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


void AliITSClusterErrAnal::GetNTeor(Int_t layer, Float_t theta, Float_t phi, Float_t &ny, Float_t &nz)
{
  //
  //get "mean shape"
  if (layer==0){
    ny = 1.+TMath::Abs(phi)*3.2;
    nz = 1.+TMath::Abs(theta)*0.34;
    return;
  }
  if (layer==1){
    ny = 1.+TMath::Abs(phi)*3.2;
    nz = 1.+TMath::Abs(theta)*0.28;
    return;
  }
  
  if (layer>3){
    ny = 2.02+TMath::Abs(phi)*1.95;
    nz = 2.02+TMath::Abs(phi)*2.35;
    return;
  }
  ny  = 6.6-2.7*abs(phi);
  nz  = 2.8-3.11*TMath::Abs(phi)+0.45*TMath::Abs(theta);
}


Int_t AliITSClusterErrAnal::GetError(Int_t layer, const AliITSclusterV2*cl, Float_t theta, Float_t phi,Float_t expQ, Float_t &erry, Float_t &errz)
{
  //calculate cluster position error
  //
  Float_t nz,ny;
  GetNTeor(layer, theta,phi,ny,nz);  
  erry   = TMath::Sqrt(cl->GetSigmaY2()); 
  errz   = TMath::Sqrt(cl->GetSigmaZ2()); 
  //
  // PIXELS
  if (layer<2){
    
    if (TMath::Abs(ny-cl->GetNy())>0.6)  {
      if (ny<cl->GetNy()){
	erry*=0.4+TMath::Abs(ny-cl->GetNy());
	errz*=0.4+TMath::Abs(ny-cl->GetNy());
      }else{
	erry*=0.7+0.5*TMath::Abs(ny-cl->GetNy());
	errz*=0.7+0.5*TMath::Abs(ny-cl->GetNy());
      }
    }
    if (TMath::Abs(nz-cl->GetNz())>1.)  {
      erry*=TMath::Abs(nz-cl->GetNz());
      errz*=TMath::Abs(nz-cl->GetNz());	      
    }
    erry*=0.85;
    errz*=0.85;
    erry= TMath::Min(erry,float(0.005));
    errz= TMath::Min(errz,float(0.03));
    return 10;
  }
  //STRIPS
  if (layer>3){ 
    if (cl->GetNy()==100||cl->GetNz()==100){
      erry = 0.004;
      errz = 0.2;
      return 100;
    }
    if (cl->GetNy()+cl->GetNz()>12){
      erry = 0.06;
      errz = 0.57;
      return 100;
    }
    Float_t normq = cl->GetQ()/(TMath::Sqrt(1+theta*theta+phi*phi));
    Float_t chargematch = TMath::Max(double(normq/expQ),2.);
    //
    if (cl->GetType()==1 || cl->GetType()==10 ){     							       
      if (chargematch<1.0 || (cl->GetNy()+cl->GetNz()<nz+ny+0.5)){
	errz = 0.043;
	erry = 0.00094;
	return 101;
      }
      if (cl->GetNy()+cl->GetNz()<nz+ny+1.2){
	errz = 0.06;
	erry =0.0013;
	return 102;
      }
      erry = 0.0027;
      errz = TMath::Min(0.028*(chargematch+cl->GetNy()+cl->GetNz()-nz+ny),0.15);
      return 103;
    }
    if (cl->GetType()==2 || cl->GetType()==11 ){ 
      erry = TMath::Min(0.0010*(1+chargematch+cl->GetNy()+cl->GetNz()-nz+ny),0.05);
      errz = TMath::Min(0.025*(1+chargematch+cl->GetNy()+cl->GetNz()-nz+ny),0.5);
      return 104;
    }
    
    if (cl->GetType()>100 ){     							       
      if ((chargematch+cl->GetNy()+cl->GetNz()-nz-ny<1.5)){
	errz = 0.05;
	erry = 0.00096;
	return 105;
      }
      if (cl->GetNy()+cl->GetNz()-nz-ny<1){
	errz = 0.10;
	erry = 0.0025;
	return 106;
      }

      errz = TMath::Min(0.05*(chargematch+cl->GetNy()+cl->GetNz()-nz-ny),0.4);
      erry = TMath::Min(0.003*(chargematch+cl->GetNy()+cl->GetNz()-nz-ny),0.05);
      return 107;
    }    
    Float_t diff = cl->GetNy()+cl->GetNz()-ny-nz;
    if (diff<1) diff=1;
    if (diff>4) diff=4;
        
    if (cl->GetType()==5||cl->GetType()==6||cl->GetType()==7||cl->GetType()==8){
      errz = 0.14*diff;
      erry = 0.003*diff;
      return 108;
    }  
    erry = 0.04*diff;
    errz = 0.06*diff;
    return 109;
  }

  //DRIFTS
  Float_t normq = cl->GetQ()/(TMath::Sqrt(1+theta*theta+phi*phi));
  Float_t chargematch = normq/expQ;
  Float_t factorz=1;
  Int_t   cnz = cl->GetNz()%10;
  //charge match
  if (cl->GetType()==1){
    if (chargematch<1.25){
      erry =  0.0028*(1.+6./cl->GetQ());  // gold clusters
    }
    else{
      erry = 0.003*chargematch;
      if (cl->GetNz()==3) erry*=1.5;
    }
    if (chargematch<1.0){
      errz =  0.0011*(1.+6./cl->GetQ());
    }
    else{
      errz = 0.002*(1+2*(chargematch-1.));
    }
    if (cnz>nz+0.6) {
      erry*=(cnz-nz+0.5);
      errz*=1.4*(cnz-nz+0.5);
    }
  }
  if (cl->GetType()>1){
    if (chargematch<1){
      erry =  0.00385*(1.+6./cl->GetQ());  // gold clusters
      errz =  0.0016*(1.+6./cl->GetQ());
    }
    else{
      errz = 0.0014*(1+3*(chargematch-1.));
      erry = 0.003*(1+3*(chargematch-1.));
    } 
    if (cnz>nz+0.6) {
      erry*=(cnz-nz+0.5);
      errz*=1.4*(cnz-nz+0.5);
    }
  }

  if (TMath::Abs(cl->GetY())>2.5){
    factorz*=1+2*(TMath::Abs(cl->GetY())-2.5);
  }
  if (TMath::Abs(cl->GetY())<1){
    factorz*=1.+0.5*TMath::Abs(TMath::Abs(cl->GetY())-1.);
  }
  factorz= TMath::Min(factorz,float(4.));  
  errz*=factorz;

  erry= TMath::Min(erry,float(0.05));
  errz= TMath::Min(errz,float(0.05));  
  return 200;
}


AliITSclusterV2 * AliITSClusterErrAnal::FindNearestCluster(AliITSCI * clinfo,  Float_t dmax)
{
  //
  //
  AliTrackReference * ref = &(clinfo->fRef);
  TObjArray * clusters = (TObjArray*)(fClusters->At(ref->GetTrack()));
  if (!clusters) return 0;
  Int_t nclusters = clusters->GetEntriesFast();  
  Float_t maxd = dmax;
  AliITSclusterV2 * cluster=0;
  Float_t radius = ref->R();
  clinfo->fNClusters=0;
  for (Int_t icluster=0;icluster<nclusters;icluster++){
    AliITSclusterV2 * cl = (AliITSclusterV2*)clusters->At(icluster);
    if (!cl) continue;
    Float_t dz = cl->GetZ()-ref->Z();
    //
    if (TMath::Abs(dz)>dmax) continue;        
    Int_t idet=cl->GetDetectorIndex();        
    Double_t r=-1;
    Float_t  dy =0;
    //    for (Int_t i=cl->fLayer;i<=cl->fLayer-1;i++){
    
    AliITStrackerMI::AliITSdetector & det = fTracker->GetDetector(cl->fLayer,idet);
      if (TMath::Abs(radius-det.GetR())<1.5) {
	Bool_t isOK=kFALSE;
	for (Int_t itrack=0;itrack<3;itrack++){
	  if (cl->fTracks[itrack]== clinfo->fRef.GetTrack()) {
	    clinfo->fNClusters++;
	    isOK= kTRUE;
	  }
	}
	if (!isOK) continue;
	r = det.GetR();
	Double_t phicl= det.GetPhi();
	Double_t ca=TMath::Cos(phicl), sa=TMath::Sin(phicl);
	Double_t x =  double(ref->X())*ca +double(ref->Y())*sa;
	Double_t y = -double(ref->X())*sa +double(ref->Y())*ca;
	if (TMath::Abs(x-r-0.012)>0.1) continue;
	//
	Double_t px = ref->Px(), py = ref->Py();
	Double_t lpx =  double(px)*ca +double(py)*sa;
	Double_t lpy =  -double(px)*sa +double(py)*ca;
	Double_t lphi =  TMath::ATan2(lpy,lpx);

	Double_t theta  = ref->Pz()/ref->Pt();
	Double_t dx = x-r;
	dy = y       -cl->GetY()-dx*lphi;
	dz = ref->Z()-cl->GetZ()-dx*theta;
	//dx-=0.012;
	Double_t d = TMath::Sqrt(dy*dy+dz*dz);
	if (d<maxd){
	  cluster=cl;
	  maxd=d;
	  clinfo->Update(fTracker);
	  clinfo->fCl = *cluster;
	  clinfo->fLx    = x;
	  clinfo->fLy    = y;
	  clinfo->fLayer = cl->fLayer;
	  clinfo->fPhiCl = phicl;
	  clinfo->fRCl   = det.GetR();
	  clinfo->fPhiL  = lphi;
	  clinfo->fDz    = dz;
	  clinfo->fDy    = dy;
	  cl->SetSigmaZ2(TMath::Abs(cl->GetSigmaZ2()));
	  cl->SetSigmaY2(TMath::Abs(cl->GetSigmaY2()));
	  clinfo->fPoolZ  = dz/TMath::Sqrt(cl->GetSigmaZ2());
	  clinfo->fPoolY  = dy/TMath::Sqrt(cl->GetSigmaY2());
	  clinfo->fPoolY2  = dy/TMath::Sqrt(cl->GetSigmaY2());	  
	  clinfo->fErrY   = TMath::Sqrt(cl->GetSigmaY2()); 
	  clinfo->fErrZ   = TMath::Sqrt(cl->GetSigmaZ2()); 
	  Float_t erry,errz,ny,nz;
	  clinfo->fErrType = GetError(cl->fLayer,cl,theta,lphi,clinfo->fExpectedQ,erry,errz);
	  GetNTeor(cl->fLayer, theta,lphi,ny,nz);  
	  clinfo->fNormQ = cl->GetQ()/TMath::Sqrt(1+theta*theta+lpy*lpy/(lpx*lpx));
	  clinfo->fTNy  = ny;
	  clinfo->fTNz  = nz;
	  clinfo->fErrY = erry;
	  clinfo->fErrZ = errz;	 
	  clinfo->fPoolY2 = dy/clinfo->fErrY;
	  clinfo->fPoolZ2 = dz/clinfo->fErrZ;
	  clinfo->fLabPos =-1;
	  Double_t xyz[3];
	  det.GetGlobalXYZ(cl,xyz);
	  clinfo->fGx = xyz[0];
	  clinfo->fGy = xyz[1];
	  for (Int_t i=0;i<3;i++) if (cl->GetLabel(i)==clinfo->fRef.GetTrack()) clinfo->fLabPos = i;
	}
      }
    }
  //  }
  return cluster;
}



Int_t AliITSClusterErrAnal::Analyze(Int_t trackmax) {    
  //
  //
  // dummy cluster to be fill if not cluster info
  //
  AliITSCI * clinfo = new AliITSCI;
  AliITSclusterV2 dummy;
  TBranch *  branch = fTreeA->GetBranch("itscl");
  branch->SetAddress(&clinfo);
  Int_t nall =0;
  Int_t withcl =0;
  Int_t ntracks = fReferences->GetEntriesFast();
  ntracks = TMath::Min(trackmax,ntracks);
  for (Int_t itrack=0;itrack<ntracks;itrack++){
    TClonesArray * references = (TClonesArray*) fReferences->At(itrack);
    if (!references) continue;
    Int_t nreferences = references->GetEntriesFast();
    if (nreferences<4) continue;
    TObjArray * clarray = (TObjArray*)fClusters->At(itrack);
    for (Int_t iref=0;iref<nreferences;iref++){
      AliTrackReference* trackref = (AliTrackReference*)references->At(iref);     
      if (trackref&& trackref->P()>0.01){
	TParticle * p = fStack->Particle(trackref->GetTrack());
	if (p) clinfo->fP = *p;
	clinfo->fRef = *trackref;
	clinfo->fStatus = 0;
	clinfo->Update(fTracker);
	nall++;
	AliITSclusterV2 * cluster = FindNearestCluster(clinfo,3.0);
	if (cluster){
	  //clinfo->fCl = *cluster;
	  clinfo->fStatus=1;
	  withcl++;
	}else{
	  clinfo->fCl= dummy;
	}
	clinfo->Update(fTracker);
	fTreeA->Fill();
      }
    }
  }
  
  fFileA->cd();
  fTreeA->Write();
  fTreeB->Write();
  fFileA->Close();
  return 0;
}



void AliITSClusterErrAnal::MakeSeeds(Double_t zv)
{
  //
  //
  // 

  AliITStrackerMI::AliITSlayer & layer3 =   fTracker->GetLayer(3); 
  AliITStrackerMI::AliITSlayer & layer2 =   fTracker->GetLayer(2); 
  AliITStrackerMI::AliITSlayer & layer1 =   fTracker->GetLayer(1); 
  AliITStrackerMI::AliITSlayer & layer0 =   fTracker->GetLayer(0); 
  Int_t ncl3 =  layer3.GetNumberOfClusters();
  Int_t ncl2 =  layer2.GetNumberOfClusters();
  Int_t ncl1 =  layer1.GetNumberOfClusters();
  Int_t ncl0 =  layer0.GetNumberOfClusters();
  TH1F *hc = new TH1F("hc","hc",100,-0.01,0.01);
  TH1F *hrc = new TH1F("hrc","hrc",100,-1,1);
  TH1F *hl = new TH1F("hl","hl",100,-0.05,0.05);
  TH1F *hlr = new TH1F("hlr","hlr",1000,-30.,30.);
  TH2F *hcl = new TH2F("hcl","hcl",1000,-0.01,0.01,1000,-0.03,0.03);
  TH2F *hcc = new TH2F("hcc","hcc",1000,-0.01,0.01,1000,-0.03,0.03);
  TH1F *hr = new TH1F("hr","hr",100,0,2000.);
  //
  TH1F *hpoolc = new TH1F("hpoolc","hpoolc",100,-5,5);
  TH1F *hpooll = new TH1F("hpooll","hpooll",100,-5,5);
  //
  Double_t ratio23 = layer2.GetR()/layer3.GetR();
  Double_t ratio12 = layer1.GetR()/layer2.GetR();
  Double_t referenceR =  (layer0.GetR()+layer1.GetR()+layer2.GetR()+layer3.GetR())*0.25;
  //
  Double_t *clusterxyzr3 = new Double_t[ncl3*6];
  Double_t *clusterxyzr2 = new Double_t[ncl2*6];
  Double_t *clusterxyzr1 = new Double_t[ncl1*6];
  Double_t *clusterxyzr0 = new Double_t[ncl0*6];
  //
  // Get global coordinates 3
  for (Int_t icl3=0;icl3<ncl3;icl3++){
    AliITSclusterV2* cl3 = layer3.GetCluster(icl3);
    Double_t *xyz3 = &clusterxyzr3[6*icl3];
    AliITStrackerMI::AliITSdetector & det = layer3.GetDetector(cl3->GetDetectorIndex());
    det.GetGlobalXYZ(cl3,xyz3);
    xyz3[3] = TMath::Sqrt(xyz3[0]*xyz3[0]+xyz3[1]*xyz3[1]);
    Double_t fi3 = TMath::ATan2(xyz3[1],xyz3[0]);
    xyz3[4] = TMath::Cos(fi3);
    xyz3[5] = TMath::Sin(fi3);
  }  
  //
  // Get global coordinates 2
  for (Int_t icl2=0;icl2<ncl2;icl2++){
    AliITSclusterV2* cl2 = layer2.GetCluster(icl2);
    Double_t *xyz2 = &clusterxyzr2[6*icl2];
    AliITStrackerMI::AliITSdetector & det = layer2.GetDetector(cl2->GetDetectorIndex());
    det.GetGlobalXYZ(cl2,xyz2);
    xyz2[3] = TMath::Sqrt(xyz2[0]*xyz2[0]+xyz2[1]*xyz2[1]);
    Double_t fi2 = TMath::ATan2(xyz2[1],xyz2[0]);
    xyz2[4] = TMath::Cos(fi2);
    xyz2[5] = TMath::Sin(fi2);
  }  
  //
  // Get global coordinates 1
  for (Int_t icl1=0;icl1<ncl1;icl1++){
    AliITSclusterV2* cl1 = layer1.GetCluster(icl1);
    Double_t *xyz1 = &clusterxyzr1[6*icl1];
    AliITStrackerMI::AliITSdetector & det = layer1.GetDetector(cl1->GetDetectorIndex());
    det.GetGlobalXYZ(cl1,xyz1);
    xyz1[3] = TMath::Sqrt(xyz1[0]*xyz1[0]+xyz1[1]*xyz1[1]);
    Double_t fi1 = TMath::ATan2(xyz1[1],xyz1[0]);
    xyz1[4] = TMath::Cos(fi1);
    xyz1[5] = TMath::Sin(fi1);
  }  
  //
  // Get global coordinates 0
  for (Int_t icl0=0;icl0<ncl0;icl0++){
    AliITSclusterV2* cl0 = layer0.GetCluster(icl0);
    Double_t *xyz0 = &clusterxyzr0[6*icl0];
    AliITStrackerMI::AliITSdetector & det = layer0.GetDetector(cl0->GetDetectorIndex());
    det.GetGlobalXYZ(cl0,xyz0);
    xyz0[3] = TMath::Sqrt(xyz0[0]*xyz0[0]+xyz0[1]*xyz0[1]);
    Double_t fi0 = TMath::ATan2(xyz0[1],xyz0[0]);
    xyz0[4] = TMath::Cos(fi0);
    xyz0[5] = TMath::Sin(fi0);

  }  
  //
  //
  Int_t ntracks=0;
  Int_t ntracks0=0;
  Int_t good0=0;
  Int_t accepted0=0;
  //
  Int_t nfoundtracks=0;
  Int_t good1 =0;;
  Int_t accepted1=0;
  Int_t accepted2=0;
  //
  Int_t *trackindex = new Int_t[ncl3*4*300];
  const AliITSclusterV2 **trackclusters = new const AliITSclusterV2*[ncl3*4*300];
  //
  for (Int_t icl3=0;icl3<ncl3;icl3++){
    AliITSclusterV2* cl3 = layer3.GetCluster(icl3);
    Double_t *xyz3 = &clusterxyzr3[icl3*6];
    Double_t &r3   = clusterxyzr3[icl3*6+3];
    ntracks=0;  //number of possible tracks corresponding to given cluster
    //
    //
    // find pairs in layer 2
    Float_t zc   = (xyz3[2]-zv)*ratio23+zv;
    Float_t lphi  = TMath::ATan2(xyz3[1],xyz3[0]);
    layer2.SelectClusters(zc-1,zc+1, layer2.GetR()*lphi-1, layer2.GetR()*lphi+1.);
    const AliITSclusterV2 * cl2=0;
    Int_t npairs2=0;
    Int_t ci2[100];
    while ( (cl2=layer2.GetNextCluster(ci2[npairs2]))!=0){
      Int_t index2 = ci2[npairs2];
      Double_t *xyz2 = &clusterxyzr2[index2*6];
      Double_t &r2   = clusterxyzr2[index2*6+3];
      //      
      if (TMath::Abs((xyz3[2]-zv)*r2/r3-(xyz2[2]-zv))>0.3) continue;
      npairs2++;
    }    
    //
    //
    // find pairs in layer 2
    for (Int_t jcl2 =0;jcl2<npairs2; jcl2++){
      Int_t index2 = ci2[jcl2];
      Double_t *xyz2 = &clusterxyzr2[index2*6];
      Float_t zc   = (xyz2[2]-zv)*ratio12+zv;
      Float_t lphi  = TMath::ATan2(xyz2[1],xyz2[0]);
      cl2 = layer2.GetCluster(index2);
      //
      //      
      layer1.SelectClusters(zc-1,zc+1, layer1.GetR()*lphi-1, layer1.GetR()*lphi+1.);
      const AliITSclusterV2 * cl1=0;
      Int_t ci1;
      //
      Double_t lr2 = (xyz3[0]-xyz2[0])*(xyz3[0]-xyz2[0])+(xyz3[1]-xyz2[1])*(xyz3[1]-xyz2[1]);
      Double_t lr  = TMath::Sqrt(lr2); 
      Double_t sinphi = (xyz3[1]-xyz2[1])/lr;
      Double_t cosphi = (xyz3[0]-xyz2[0])/lr;
      Double_t t1 = (xyz3[2]-xyz2[2])/lr;
      //
      while ( (cl1=layer1.GetNextCluster(ci1))!=0){
	if (cl1->GetQ()==0) continue;
	Int_t index1 = ci1;
	Double_t *xyz1 = &clusterxyzr1[index1*6];
	AliITSclusterV2* cl1 = layer1.GetCluster(index1);
  
	if (cl3->GetLabel(0)==cl2->GetLabel(0) && cl3->GetLabel(0)==cl1->GetLabel(0)) good0++;
	Double_t loc1[2];
	//
	loc1[0] =  (xyz1[0]-xyz2[0])*cosphi+(xyz1[1]-xyz2[1])*sinphi;
	loc1[1] = -(xyz1[0]-xyz2[0])*sinphi+(xyz1[1]-xyz2[1])*cosphi;
	//
	Double_t c1 = 2.*loc1[1]/(loc1[0]*(loc1[0]-lr));   //curvature approximated using 2 derivation
	Double_t c=c1;
	Double_t t2 = (xyz2[2]-xyz1[2])/TMath::Sqrt(loc1[0]*loc1[0]+loc1[1]*loc1[1]); 
	if (TMath::Abs(c)>0.008) continue;
	if (TMath::Abs(t2-t1)>0.02) continue;
	//	
	trackclusters[ntracks*4]    = 0;
	trackclusters[ntracks*4+1]  = cl1;
	trackclusters[ntracks*4+2]  = cl2;
	trackclusters[ntracks*4+3]  = cl3;
	trackindex[ntracks*4+0] = 0;
	trackindex[ntracks*4+1] = index1;
	trackindex[ntracks*4+2] = index2;
	trackindex[ntracks*4+3] = icl3;
	if (ntracks>=300*ncl3){	  
	  break;
	}
	if (cl3->GetLabel(0)==cl2->GetLabel(0) && cl3->GetLabel(0)==cl1->GetLabel(0)) {
	  //printf("%f\n",c);
	  
	  //hc->Fill(c);
	  //hl->Fill(t2-t1);
	  //hr->Fill(1/c);
	  //hcc->Fill(c,TMath::Abs(c1));
	  //hcl->Fill(c1,t2-t1);
	  
	  accepted0++;
	}
	ntracks++;
      }
    }   // end - find pairs in layer 2
    //    
    ntracks0+=ntracks;
    //hcc->Draw();
    Int_t   ntracks2=0;
    Int_t   cindex2[10000*4];
    Float_t pools[10000];
    
    for (Int_t itrack=0;itrack<ntracks;itrack++){
      const AliITSclusterV2 *cl0, *cl1,*cl2,*cl3;
      cl3 = trackclusters[itrack*4+3];
      cl2 = trackclusters[itrack*4+2];
      cl1 = trackclusters[itrack*4+1];
      Int_t index1 = trackindex[itrack*4+1];
      Int_t index2 = trackindex[itrack*4+2];
      Int_t index3 = trackindex[itrack*4+3];
      //
      Double_t *xyz1 = &clusterxyzr1[index1*6],*xyz2 = &clusterxyzr2[index2*6],*xyz3=&clusterxyzr3[index3*6];
      Double_t &r1 = clusterxyzr1[index1*6+3],&r2 = clusterxyzr2[index2*6+3], &r3 = clusterxyzr3[index3*6+3];
      
      Double_t r0 = layer0.GetR();
      Float_t zc   = xyz1[2] + (xyz2[2]-xyz1[2])*(r0-r1)/(r2-r1);
      Float_t lphi  = TMath::ATan2(xyz1[1],xyz1[0]);
      
      layer0.SelectClusters(zc-1,zc+1, layer0.GetR()*lphi-0.3, layer0.GetR()*lphi+0.3);
      Int_t ci0;
      //
      //
      Double_t lr2 = (xyz2[0]-xyz1[0])*(xyz2[0]-xyz1[0])+(xyz2[1]-xyz1[1])*(xyz2[1]-xyz1[1]);
      Double_t lr  = TMath::Sqrt(lr2); 
      Double_t sinphi = (xyz2[1]-xyz1[1])/lr;
      Double_t cosphi = (xyz2[0]-xyz1[0])/lr;
      //
      Double_t loc3[2]; // cluster3 in local coordindate 
      loc3[0] =  (xyz3[0]-xyz1[0])*cosphi+(xyz3[1]-xyz1[1])*sinphi;
      loc3[1] = -(xyz3[0]-xyz1[0])*sinphi+(xyz3[1]-xyz1[1])*cosphi;
      
      Double_t tan1 = (xyz3[2]-xyz2[2])/TMath::Sqrt((xyz3[0]-xyz2[0])*(xyz3[0]-xyz2[0])+(xyz3[1]-xyz2[1])*(xyz3[1]-xyz2[1]));     
      Double_t c3 = 2.*loc3[1]/(loc3[0]*(loc3[0]-lr));   //curvature approximated using 2 derivation
      //
      while ( (cl0=layer0.GetNextCluster(ci0))!=0){ 
	if (cl0->GetQ()==0) continue;
	Double_t *xyz0 = &clusterxyzr0[ci0*6];	
	if (cl3->GetLabel(0)==cl0->GetLabel(0) && cl3->GetLabel(0)==cl2->GetLabel(0) 
	    && cl3->GetLabel(0)==cl1->GetLabel(0)&&cl3->GetLabel(0)==cl0->GetLabel(0)){
	  good1++;
	}      	
	Double_t loc0[2]; // cluster3 in local coordindate 
	loc0[0] =  (xyz0[0]-xyz1[0])*cosphi+(xyz0[1]-xyz1[1])*sinphi;
	loc0[1] = -(xyz0[0]-xyz1[0])*sinphi+(xyz0[1]-xyz1[1])*cosphi;
	//
	Double_t c0 = 2.*loc0[1]/(loc0[0]*(loc0[0]-lr));   //curvature approximated using 2 derivation
	Double_t tan2 = (xyz1[2]-xyz0[2])/TMath::Sqrt((xyz1[0]-xyz0[0])*(xyz1[0]-xyz0[0])+(xyz1[1]-xyz0[1])*(xyz1[1]-xyz0[1]));
	//if (TMath::Abs(tan1-tan2)>0.05) continue;
	//if (TMath::Abs(c3-c0)>0.006) continue;
	Double_t poolphi = (c3-c0)/0.0011; 
	poolphi*=poolphi;
	Double_t poolr = ((c3-c0)/(TMath::Abs(c3+c0)+0.001))/0.14;
	Double_t poolth = (tan1-tan2)/0.0094; 
	poolth*=poolth;
	Double_t poolthr = ((tan1-tan2)/(TMath::Abs(c3+c0)+0.005))/0.83;
	if (poolr*poolr+poolthr*poolthr>40) continue;
	//if ( TMath::Abs((c3-c0)/(c3+c0))>0.6) continue; 
	cindex2[ntracks2*4+0]   = ci0;
	cindex2[ntracks2*4+1]   = index1;
	cindex2[ntracks2*4+2]   = index2;
	cindex2[ntracks2*4+3]   = index3;
	pools[ntracks2] = poolr*poolr+poolthr*poolthr;
	//
	ntracks2++;
	if (cl3->GetLabel(0)==cl0->GetLabel(0) && cl3->GetLabel(0)==cl2->GetLabel(0) 
	    && cl3->GetLabel(0)==cl1->GetLabel(0)){
	  hl->Fill(tan1-tan2);
	  hrc->Fill(poolr);
	  hc->Fill((c3-c0));
	  hlr->Fill(poolthr);
	  hpoolc->Fill(poolr);
	  hpooll->Fill(poolthr);
	  accepted1++;
	}      	
      }
    }
    Int_t sindex[10000];
    TMath::Sort(ntracks2,pools,sindex,kFALSE);
    Float_t minpool = pools[sindex[0]];
    Int_t max = TMath::Min(ntracks2,20);
    for (Int_t itrack=0;itrack<ntracks2;itrack++){
      if (itrack>=max) continue;
      //if (pools[sindex[itrack]]>minpool+5) continue;
      const AliITSclusterV2 *cl0, *cl1,*cl2,*cl3;
      cl3 = layer3.GetCluster(cindex2[sindex[itrack]*4+3]);
      cl2 = layer2.GetCluster(cindex2[sindex[itrack]*4+2]);
      cl1 = layer1.GetCluster(cindex2[sindex[itrack]*4+1]);
      cl0 = layer0.GetCluster(cindex2[sindex[itrack]*4+0]);
      if (cl3->GetLabel(0)==cl0->GetLabel(0) && cl3->GetLabel(0)==cl2->GetLabel(0) 
	  && cl3->GetLabel(0)==cl1->GetLabel(0)){
	accepted2++;
      }
      nfoundtracks++;
    }
  }
  hc->Draw();
  printf("First pass \t%d\t%d\t%d\n",ntracks0,good0,accepted0);  
  printf("Found tracks\t%d\t%d\t%d\t%d\n",nfoundtracks,good1, accepted1,accepted2);
  
}
