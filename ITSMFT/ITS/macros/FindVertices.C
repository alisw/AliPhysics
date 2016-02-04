#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <Riostream.h>
#include <TSystem.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TBranch.h>
#include <TArrayF.h>
#include <TH1F.h>
#include <TParticle.h>

#include "AliESDVertex.h"
#include "AliRun.h"
#include "AliGenEventHeader.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliLoader.h"
#include "AliStack.h"
#include "AliGeomManager.h"
#include "AliITSUVertexer.h"
#include "AliITSUClusterPix.h"
#include "AliITSURecoDet.h"
#include "AliITSURecoLayer.h"
#include "AliITSUGeomTGeo.h"
#endif

void FindVertices(int nev=-1,int evStart=0){
  gAlice=NULL;
  // phi (rad), z(cm), pair(cm), cluster(cm)
  AliITSUVertexer *vert=new AliITSUVertexer(0.005,0.002,0.04,0.8,10);

  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  runLoader->LoadgAlice();

  gAlice = runLoader->GetAliRun();

  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadRecPoints();

  AliLoader *dl = runLoader->GetDetectorLoader("ITS");

  AliGeomManager::LoadGeometry("geometry.root");
  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE);
  Int_t nLayers = gm->GetNLayers();
  AliITSUClusterPix::SetGeom(gm);
  AliITSURecoDet *its = new AliITSURecoDet(gm, "ITSinterface");
  its->CreateClusterArrays();
 
  TTree * cluTree = 0x0;
  
  int nevTot = (Int_t)runLoader->GetNumberOfEvents();
  printf("N Events : %i \n",nevTot);
  evStart = evStart<nevTot ? evStart : nevTot-1;
  if (evStart<0) evStart = 0;
  //
  int lastEv = nev<0 ? nevTot : evStart+nev;
  if (lastEv > nevTot) lastEv = nevTot;
  //
  TArrayF mcVertex(3); 

  TFile *fint = new TFile("RecVertex.root","recreate");
  TNtuple *nt = new TNtuple("ntvert","vertices","xtrue:ytrue:ztrue:xrec:yrec:zrec:xdiff:ydiff:zdiff:rec:primaryfrac:goodlines:goodphi:lines:linesphi");
  TNtuple *nt1 = new TNtuple("pileup","verticesz","z1:z2:z3:z4:z5:z6:z7:z8:z9:z10:c1:c2:c3:c4:c5:c6:c7:c8:c9:c10");
  TH1F *pthist = new TH1F("pthist","p_{t};p_{t} (GeV/c)",40,0,10);
  TH1F *etahist = new TH1F("etahist","#eta;#eta",40,-10,10);

  UInt_t nophi=0,goodphi=0,nolines1=0,goodlines=0,goodph=0;
  for (Int_t iEvent = evStart; iEvent < lastEv; iEvent++) {
    runLoader->GetEvent(iEvent);
    printf("\n Event %i \n",iEvent);
    
    runLoader->GetHeader()->GenEventHeader()->PrimaryVertex(mcVertex);
    AliStack *stack = runLoader->Stack();

    vert->FindVertexForCurrentEvent(dl->TreeR());
    Int_t noVert=vert->GetNumOfVertices();
    AliESDVertex *recV=vert->GetAllVertices(noVert);
    cout << "Vertexing done, no vert: " << noVert << endl;
    UInt_t nolines;
    #ifdef MC_CHECK
    UInt_t *parr = vert->GetParticleId(nolines);
    #endif
    Int_t cand=0;
    Float_t bestDist=1e27;
    for(UInt_t ii=0;ii<noVert;++ii) {
      Double_t p[3];
      recV[ii].GetXYZ(p);
      Float_t dx=p[0]-mcVertex[0];
      Float_t dy=p[1]-mcVertex[1];
      Float_t dz=p[2]-mcVertex[2];
      Float_t distance=dx*dx+dy*dy+dz*dz;
      if(distance<bestDist) {
        bestDist=distance;
        cand=ii;
      }
    }

    Float_t nt1z[20];for(Int_t iiii=0;iiii<20;++iiii) nt1z[iiii]=-1; 
    Int_t count=0;
    Double_t ref[3];
    
    if(noVert>1) gSystem->Exec(Form("echo \"Event %i\">>pileup",iEvent));

    for(UInt_t ii=0;ii<noVert;++ii) {
      recV[cand].GetXYZ(ref);
      if(ii==cand) continue;
      if(ii>9) break;
      Double_t p[3];
      recV[ii].GetXYZ(p);
      
      Float_t dx=TMath::Abs(p[0]-ref[0]);
      Float_t dy=TMath::Abs(p[1]-ref[1]);
      Float_t dz=TMath::Abs(p[2]-ref[2]);
      nt1z[count+10]=recV[ii].GetNContributors();
      nt1z[count++]=dz;
    }
    nt1->Fill(nt1z);

    #ifdef MC_CHECK
    UInt_t primary=0;
    for(Int_t iPart=0; iPart<nolines; iPart++) {
      TParticle* part = (TParticle*)stack->Particle(parr[iPart]);
      etahist->Fill(part->Eta());
      pthist->Fill(part->Pt());
      if(part->IsPrimary()) primary++;
    }

    nophi=vert->GetLinesPhi();
    goodph=vert->GetGoodLinesPhi();
    nolines1=vert->GetNoLines();
    goodlines=vert->GetGoodLines();
    cout << "Primaries: " << primary <<" Good lines: " << nolines << " Lines: " << vert->GetNoLines() << endl;
    #endif
    Float_t xnt[17];
    xnt[0]=mcVertex[0];
    xnt[1]=mcVertex[1];
    xnt[2]=mcVertex[2];
    if(noVert!=0) {
      xnt[3]=recV[cand].GetX();
      xnt[4]=recV[cand].GetY();
      xnt[5]=recV[cand].GetZ();
      Double_t cov[6];
      recV[cand].GetCovarianceMatrix(cov);
      xnt[11]=cov[0];
      xnt[12]=cov[3];
      xnt[13]=cov[5];
    } else xnt[11]=xnt[12]=xnt[13]=xnt[3]=xnt[4]=xnt[5]=0;
    xnt[6]=mcVertex[0]-xnt[3];
    xnt[7]=mcVertex[1]-xnt[4];
    xnt[8]=mcVertex[2]-xnt[5];
    xnt[9]=noVert;
    #ifdef MC_CHECK
    xnt[10]=nolines-primary;
    #else
    xnt[10]=nolines;
    #endif
    xnt[11]=goodlines;
    xnt[12]=goodph;
    xnt[13]=nolines1;
    xnt[14]=nophi;

    Float_t dx = TMath::Abs(mcVertex[0]-xnt[3])*10000;
    Float_t dy = TMath::Abs(mcVertex[1]-xnt[4])*10000;
    Float_t dz = TMath::Abs(mcVertex[2]-xnt[5])*10000;

    cout << noVert << endl;
    if((dx>=500||dy>=500||dz>=500)&&noVert>0) gSystem->Exec(Form("echo \"Event %i\">>outliers",iEvent));
    
    if(noVert>0) {
      cout << "Contribs: " << recV[cand].GetNContributors() << endl;
      cout << mcVertex[0] << "\t" << xnt[3] << "\t" << dx << endl;
      cout << mcVertex[1] << "\t" << xnt[4] << "\t" << dy << endl;
      cout << mcVertex[2] << "\t" << xnt[5] << "\t" << dz << endl;
    }      
    nt->Fill(xnt);
  }//event loop

  fint->cd();
  pthist->Write();
  etahist->Write();
  nt->Write();
  nt1->Write();
  fint->Close();
}
