#ifndef __CINT__
#include "alles.h"
#include "GetGoodParticles.h"
#include "AliL3Transform.h"
#endif


void GetGoodParticles(Int_t minslice,Int_t maxslice,char *eventfile,char *digitfile)
{

  Int_t good_number = 70;
  
  struct GoodTrack goodtracks[15000];
  Int_t nt=0;
  char *fname = "good_tracks";
  
  ifstream in(fname);
  if(in)
    {
      printf("Delete you old good_tracks file\n");
      return;
    }
  
  TFile *evfile = TFile::Open(eventfile);
  if(!evfile)
    {
      printf("No eventfile\n");
      return;
    }
  gAlice = (AliRun*)evfile->Get("gAlice");
  gAlice->GetEvent(0);

  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");
  AliTPCParam *param = (AliTPCParam*)evfile->Get("75x40_100x60");
  
  TPC->SetParam(param);
  
  Int_t ver = TPC->IsVersion();
  if(ver!=2)
    {
      printf("Wrong TPC version\n");
      return;
    }
  
  Int_t zero=TPC->GetParam()->GetZeroSup();
  
  Int_t np = gAlice->GetNtrack();
  Int_t *good = new Int_t[np];
  for(Int_t ii=0; ii<np; ii++)
    good[ii] = 0;

  AliL3Transform *transform = new AliL3Transform();

  TFile *digfile = TFile::Open(digitfile);
  TTree *tree = (TTree*)digfile->Get("TreeD_75x40_100x60_0");
  AliSimDigits da, *digits=&da;
  tree->GetBranch("Segment")->SetAddress(&digits);
  
  Int_t *count = new Int_t[np]; //np number of particles.
  Int_t i;
  for (i=0; i<np; i++) count[i]=0;
  for(i=0; i<tree->GetEntries(); i++)
    {
      if (!tree->GetEvent(i)) continue;
      Int_t sec,row;
      param->AdjustSectorRow(digits->GetID(),sec,row);
      Int_t sl,padrow;
      transform->Sector2Slice(sl,padrow,sec,row);
      if(sl < minslice) continue;
      if(sl > maxslice) break;
      digits->First();
      do {
	Int_t it=digits->CurrentRow(), ip=digits->CurrentColumn();
	Short_t dig = digits->GetDigit(it,ip);
	Int_t idx0=digits->GetTrackID(it,ip,0); 
	Int_t idx1=digits->GetTrackID(it,ip,1);
	Int_t idx2=digits->GetTrackID(it,ip,2);
	
	if (idx0>=0 && dig>=zero) count[idx0]+=1;
	if (idx1>=0 && dig>=zero) count[idx1]+=1;
	if (idx2>=0 && dig>=zero) count[idx2]+=1;
	  } while (digits->Next());
      
      for (Int_t j=0; j<np; j++) 
	{
	  if (count[j]>1) //at least two digits at this padrow 
	    good[j]++;
	  
	  count[j]=0;
	}
    }
    
  delete[] count;
  
  
  //The following code has been taken from offlinemacro->AliTPCComparison.C
  
  TTree *TH=gAlice->TreeH();
  Int_t npart=(Int_t)TH->GetEntries();
    


  while (npart--) {
    AliTPChit *hit0=0;
    
    TPC->ResetHits();
    TH->GetEvent(npart);
    AliTPChit * hit = (AliTPChit*) TPC->FirstHit(-1);
    while (hit){
      if (hit->fQ==0.) break;
      hit =  (AliTPChit*) TPC->NextHit();
    }
    if (hit) {
      hit0 = new AliTPChit(*hit); //Make copy of hit
      hit = hit0;
    }
    else continue;
    AliTPChit *hit1=(AliTPChit*)TPC->NextHit();       
    if (hit1==0) continue;
    if (hit1->fQ != 0.) continue;
    Int_t i=hit->Track();
    TParticle *p = (TParticle*)gAlice->Particle(i);
    
    if (p->GetFirstMother()>=0) continue;  //secondary particle
    if (good[i] < good_number) continue;
    if (p->Pt() < 0.1) continue;
    if (TMath::Abs(p->Pz()/p->Pt())>0.999) continue;
    printf("checking paricle %d with %d hits \n",i,good[i]);
    
    goodtracks[nt].eta = p->Eta();
    goodtracks[nt].label=i;
    goodtracks[nt].code=p->GetPdgCode();
    goodtracks[nt].px=hit->X(); goodtracks[nt].py=hit->Y(); goodtracks[nt].pz=hit->Z();
    goodtracks[nt].pt = p->Pt();
    goodtracks[nt].nhits = good[i];
    nt++;
    
    if (hit0) delete hit0;
  }
  
  ofstream out(fname);
  if(out) 
    {
      for (Int_t ngd=0; ngd<nt; ngd++)            
	out<<goodtracks[ngd].label<<' '<<goodtracks[ngd].code<<' '<<
	  goodtracks[ngd].px<<' '<<goodtracks[ngd].py<<' '<<goodtracks[ngd].pz<<' '<<
	  goodtracks[ngd].pt<<' '<<goodtracks[ngd].eta<<' '<<goodtracks[ngd].nhits<<endl; 
      
    } 
  else
    printf("cannot open file\n");
  
  out.close();
  
  delete [] good;

  evfile->Close();
  digfile->Close();

  delete transform;

}
