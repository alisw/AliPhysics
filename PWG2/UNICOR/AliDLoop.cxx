//Author: Dariusz Miskowiec 2007

//=============================================================================
// simple event loop manager
// First, events are classified according to centrality, reaction plane angle, 
// and z-vertex. Then, single particles and true pairs are processed. Finally, 
// event mixing is done within event classes. For this, the tree and the event 
// passed as arguments are cloned. 
// The class is, at present, deliberately left in a macro-like design, with all 
// the essential things happening inside the function AliDLoop::Run. This is 
// subject to future polishing if the class proves useful. 
//=============================================================================

#include <TROOT.h>
#include <TMath.h>
#include <TTree.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include "AliDLoop.h"
#include "AliDAnalGlobal.h"
#include "AliDAnalSingle.h"
#include "AliDAnalCorrel.h"
#include "AliDAnalPtfluc.h"
#include "AliDEvent.h"

ClassImp(AliDLoop)

//=============================================================================
AliDLoop::AliDLoop(TTree *tr, AliDEvent *ev0, char *outfil): 
  TObject(), 
  fTree0(tr), fTree1(0x0), 
  fEv0(ev0), fEv1(0x0), 
  fOutputFilename(outfil)     
{
  // constructor
  // Clone the tree such that during the event mixing the two events get 
  // populated from two different trees. This costs some memory (37 MB for 
  // ceres) but allows to avoid switching one tree between two different 
  // events, especially difficult for alice because alice trees do not like 
  // to be attached more than once. 

  fTree1 = (TTree*) tr->Clone();
  fEv1 = (AliDEvent*) ev0->Clone();
  fEv0->AttachTree(fTree0);
  fEv1->AttachTree(fTree1);
}
//=============================================================================
AliDLoop::AliDLoop(const AliDLoop &loop) : 
  TObject(loop), fTree0(loop.fTree0), fTree1(loop.fTree1), fEv0(loop.fEv0), 
  fEv1(loop.fEv1), fOutputFilename(loop.fOutputFilename) 
{     
  //copy constructor, shallow copy, just to make the compiler happy
}
//=============================================================================
AliDLoop &AliDLoop::operator=(const AliDLoop &source) 
{
  // substitution operator, shallow copy, just to make the compiler happy

  TObject::operator=(source);
  fTree0 = source.fTree0;
  fTree1 = source.fTree1;
  fEv0 = source.fEv0;
  fEv1 = source.fEv1;
  fOutputFilename = source.fOutputFilename;
  return *this;
}
//=============================================================================
Int_t AliDLoop::Mem() const 
{
  // return memory in MB occupied by this process

  ProcInfo_t info;
  gSystem->GetProcInfo(&info);
  Int_t memmb = info.fMemVirtual/1024;
  return memmb;
}
//=============================================================================
void AliDLoop::Run(int nwish) const 
{
  // process nwish events from fTree0; nwish=-1 (default) means process all

  // control parameters


  //  const int kcen = 10; // number of centrality bins for mixing
  //  const int kphi = 8; // number of phi bins for mixing
  //  const int kzve = 13; // number of z-vertex bins for mixing
  const int kcen = 1; // number of centrality bins for mixing
  const int kphi = 1; // number of phi bins for mixing
  const int kzve = 1; // number of z-vertex bins for mixing
  int mixingFactor = 1; // number of wished event pairs/number of events
  int mixingStep[10] = {2,3,5,7,11,13,17,19,23,29};
  if (mixingFactor>=10) {printf("mixing factor too high\n"); exit(-1);}
  int minn = 20*mixingFactor; // minimum number of events per bin
  int pr = !gROOT->IsBatch(); // more  printing in interactive mode

  // initialize analysis

  Double_t etami = fEv0->Etamin();
  Double_t etama = fEv0->Etamax();
  AliDAnalGlobal *dag = new AliDAnalGlobal("dag");
  AliDAnalSingle *all = new AliDAnalSingle("all",etami,etama,0);
  AliDAnalSingle *pim = new AliDAnalSingle("pim",etami,etama,-211);
  AliDAnalSingle *pip = new AliDAnalSingle("pip",etami,etama, 211);
  AliDAnalCorrel *cnn = new AliDAnalCorrel("cnn",etami,etama,-211,-211);
  AliDAnalCorrel *cpp = new AliDAnalCorrel("cpp",etami,etama, 211, 211);
  AliDAnalPtfluc *ptf = new AliDAnalPtfluc("ptf",-211,-211);

  // initialize number of events etc

  printf("\ndetermining number of events in chain...");
  int nmax = fTree0->GetEntries();
  int nact = nmax;
  if (nwish>-1 && nwish<nmax) nact = nwish;
  printf("analyzing %d events of %d\n\n",nact,nmax);
  TStopwatch sto;

  sto.Reset();
  sto.Start();

  // find groups of of similar events
  // first allocate nact entries to each array; later reduce

  printf("%4d MB; booking   event index arrays\n",Mem());
  TArrayI *ar[kcen][kphi][kzve];
  int nar[kcen][kphi][kzve];
  for (int i=0; i<kcen; i++) 
    for (int j=0; j<kphi; j++)
      for (int k=0; k<kzve; k++) {
	ar[i][j][k] = new TArrayI(nact);
	nar[i][j][k] = 0;
      }
  printf("%4d MB; filling   event index arrays\n",Mem());
  for (int i=0; i<nact; i++) {
    fTree0->GetEvent(i);
    if (!fEv0->Good()) continue;
    int icen = (int) (kcen*fEv0->Centrality());
    int iphi = (int) (kphi*(fEv0->RPphi()/TMath::Pi()+1)/2);
    int izve = (int) (kzve*(fEv0->Zver()+1.0)/2.0);
    if (pr) printf("\revent %5d   %2d%2d%2d",i,icen,iphi,izve);
    if (icen<0 || icen>=kcen) continue;
    if (iphi<0 || iphi>=kphi) continue;
    if (izve<0 || izve>=kzve) continue;
    ar[icen][iphi][izve]->AddAt(i,nar[icen][iphi][izve]);
    nar[icen][iphi][izve]++;
  }
  if (pr) printf("\n");

  printf("%4d MB; shrinking event index arrays\n",Mem());
  int nall = 0;
  int nsup = 0;
  for (int i=0; i<kcen; i++) 
    for (int j=0; j<kphi; j++)
      for (int k=0; k<kzve; k++) {
	nall += nar[i][j][k];
	// suppress scarse bins
	if (nar[i][j][k]<minn) {
	  nsup += nar[i][j][k];
	  nar[i][j][k]=0;
	}
	ar[i][j][k]->Set(nar[i][j][k]);
      }
  printf("%4d MB; %d events sorted into bins\n",Mem(),nall);
  printf("%d events suppressed because bin occupancy was below %d\n",nsup,minn);
  printf("sorting events took %.1f s CPU %.1f s real %.1f ms real per event\n\n",
	 sto.CpuTime(),sto.RealTime(),1000*sto.RealTime()/nall);
  if ((nall-=nsup) < 1) return;

  // actual loop

  printf("%4d MB; starting the main loop\n",Mem());
  sto.Reset();
  sto.Start();
  for (int i=0; i<kcen; i++) 
    for (int j=0; j<kphi; j++)
      for (int k=0; k<kzve; k++) {

	// true pairs

	int nevents = ar[i][j][k]->GetSize(); // number of events in this bin
	for (int l=0; l<nevents; l++) {
	  if (pr) printf("\revent %6d of %6d  %3d%3d%3d",l,nevents,i,j,k);
	  int m = ar[i][j][k]->At(l); 
	  fTree0->GetEvent(m);
	  dag->Process(fEv0);
	  all->Process(fEv0);
	  pim->Process(fEv0);
	  pip->Process(fEv0);
	  cnn->Process(0,fEv0,fEv0,0);
	  cnn->Process(2,fEv0,fEv0,TMath::DegToRad()*180);	  
	  cpp->Process(0,fEv0,fEv0,0);
	  cpp->Process(2,fEv0,fEv0,TMath::DegToRad()*180);	  
	  ptf->Process(0,fEv0,fEv0);	  
	}
	if (pr && nevents) printf("\n");

	// mixed pairs
	// loop somewhat optimized to reduce jumping back and forth

	for (int imix=0; imix<mixingFactor; imix++) {
	  int step = mixingStep[imix];
	  for (int phase=0; phase<step; phase++) for (int l=phase; l<nevents; l+=step) {
	    int m0 = ar[i][j][k]->At(l); 
	    int next = (l+step)%nevents;
	    int m1 = ar[i][j][k]->At(next); 
	    if (pr) printf("\rmixing %5d and %5d",m0,m1);
	    fTree0->GetEvent(m0);
	    fTree1->GetEvent(m1);
	    //printf("   event0 %f   event1 %f\n",fEv0->ParticleP(0),fEv1->ParticleP(0));
	    cnn->Process(1,fEv0,fEv1,0);
	    cpp->Process(1,fEv0,fEv1,0);
	    ptf->Process(1,fEv0,fEv1);	  
	  }
	if (pr && nevents) printf("\r");
	}
      }

  sto.Stop();
  printf("event loop took %.1f s CPU %.1f s real %.1f ms real per event\n\n",
	 sto.CpuTime(),sto.RealTime(),1000*sto.RealTime()/nall);

  // save analysis results

  dag->Save(fOutputFilename.Data(),"recreate");
  all->Save(fOutputFilename.Data());
  pim->Save(fOutputFilename.Data());
  pip->Save(fOutputFilename.Data());
  cnn->Save(fOutputFilename.Data());
  cpp->Save(fOutputFilename.Data());
  ptf->Save(fOutputFilename.Data());
}
//=============================================================================
