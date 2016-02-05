#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TObjArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "AliITSUClusterPix.h"
#include "AliITSURecoLayer.h"
#include "AliITSURecoDet.h"
#include "AliITSUHit.h"
#include "AliITSUGeomTGeo.h"
#include "AliITSMFTSegmentationPix.h"
#include "AliGeomManager.h"
#include "AliStack.h"
#include "AliLoader.h"
#include "AliCDBManager.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGeoMatrix.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TClonesArray.h"

#endif

TObjArray histoArr;
enum {kNPixAll=0,kNPixSPL=1,kDR=0,kDTXodd,kDTXeven,kDTZ, kDTXoddSPL,kDTXevenSPL,kDTZSPL};

TPaveStats* GetStPad(TH1* hst);
TPaveStats* SetStPadPos(TH1* hst,float x1,float x2,float y1,float y2, Int_t stl=-1,Int_t col=-1);
TCanvas* DrawNP(int np, TObjArray* harr=0, TCanvas* cnv=0);
TH1* GetHistoClSize(int npix,int id,TObjArray* harr=0);
void DrawReport(const char* psname, TObjArray* harr=0);

typedef struct {
  Int_t evID;
  Int_t volID;
  Int_t lrID;
  Int_t clID;
  Int_t nPix;
  Int_t nX;
  Int_t nZ;
  Int_t q;
  Float_t pt;
  Float_t eta;
  Float_t phi;
  Float_t xyz[3];
  Float_t dX;
  Float_t dY;
  Float_t dZ;  
  Bool_t split;  
  Bool_t prim;
  Int_t  pdg;
  Int_t  ntr;
  Float_t alpha; // alpha is the angle in y-radius plane in local frame
  Float_t beta;  // beta is the angle in xz plane, taken from z axis, growing counterclockwise
  Int_t nRowPatt;
  Int_t nColPatt;
} clSumm;

void compClusHits(int nev=-1)
{

  const int kSplit=0x1<<22;
  const int kSplCheck=0x1<<23;
  //
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
  gSystem->Load("libITSUpgradeRec");
  gROOT->SetStyle("Plain");

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetSpecificStorage("GRP/GRP/Data",
        Form("local://%s",gSystem->pwd()));
  man->SetSpecificStorage("ITS/Align/Data",
        Form("local://%s",gSystem->pwd()));
  man->SetSpecificStorage("ITS/Calib/RecoParam",
        Form("local://%s",gSystem->pwd()));
  man->SetRun(0);

  gAlice=NULL;
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  runLoader->LoadgAlice();

  gAlice = runLoader->GetAliRun();

  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadRecPoints();
  runLoader->LoadSDigits();
  runLoader->LoadHits();

  AliLoader *dl = runLoader->GetDetectorLoader("ITS");

  AliGeomManager::LoadGeometry("geometry.root");
  TObjArray algITS;
  AliGeomManager::LoadAlignObjsFromCDBSingleDet("ITS",algITS);
  AliGeomManager::ApplyAlignObjsToGeom(algITS);
  //
  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE);
  AliITSUClusterPix::SetGeom(gm);
  //
  AliITSURecoDet *its = new AliITSURecoDet(gm, "ITSinterface");
  its->CreateClusterArrays();
  //
  Double_t xg1,yg1,zg1=0.,xg0,yg0,zg0=0.,tg0;
  Double_t xExit,yExit,zExit,xEnt,yEnt,zEnt,tof1;

  //
  TTree * cluTree = 0x0;
  TTree *hitTree = 0x0;
  TClonesArray *hitList=new TClonesArray("AliITSUHit");
  //
  TObjArray arrMCTracks; // array of hit arrays for each particle
  //
  Float_t xyzClGloF[3];
  Double_t xyzClGlo[3],xyzClTr[3];
  Int_t labels[3];
  int nLab = 0;
  int nlr=its->GetNLayersActive();
  int ntotev = (Int_t)runLoader->GetNumberOfEvents();
  printf("N Events : %i \n",ntotev);
  if (nev>0) ntotev = TMath::Min(nev,ntotev);
  //
  
  // output tree
  TFile* flOut = TFile::Open("clInfo.root","recreate");
  TTree* trOut = new TTree("clitsu","clitsu");
  clSumm cSum;
  trOut->Branch("evID", &cSum.evID ,"evID/I");
  trOut->Branch("volID",&cSum.volID,"volID/I");
  trOut->Branch("lrID", &cSum.lrID ,"lrID/I");  
  trOut->Branch("clID", &cSum.clID ,"clID/I");  
  trOut->Branch("nPix", &cSum.nPix ,"nPix/I");
  trOut->Branch("nX"  , &cSum.nX   ,"nX/I");
  trOut->Branch("nZ"  , &cSum.nZ   ,"nZ/I");
  trOut->Branch("q"   , &cSum.q    ,"q/I");
  trOut->Branch("pt"  , &cSum.pt   ,"pt/F");  
  trOut->Branch("eta"  ,&cSum.eta  ,"eta/F");  
  trOut->Branch("phi"  , &cSum.phi  ,"phi/F");  
  trOut->Branch("xyz",   cSum.xyz,  "xyz[3]/F");  
  trOut->Branch("dX"  , &cSum.dX   ,"dX/F");
  trOut->Branch("dY"  , &cSum.dY   ,"dY/F");
  trOut->Branch("dZ"  , &cSum.dZ   ,"dZ/F");  
  trOut->Branch("split",&cSum.split,"split/O");
  trOut->Branch("prim", &cSum.prim, "prim/O");
  trOut->Branch("pdg",  &cSum.pdg,  "pdg/I");
  trOut->Branch("ntr",  &cSum.ntr,  "ntr/I");
  trOut->Branch("alpha", &cSum.alpha, "alpha/F");
  trOut->Branch("beta", &cSum.beta, "beta/F");
  trOut->Branch("nRowPatt", &cSum.nRowPatt, "nRowPatt/I");
  trOut->Branch("nColPatt", &cSum.nColPatt, "nColPatt/I");
  
  for (Int_t iEvent = 0; iEvent < ntotev; iEvent++) {
    printf("\n Event %i \n",iEvent);
    runLoader->GetEvent(iEvent);
    AliStack *stack = runLoader->Stack();
    cluTree=dl->TreeR();
    hitTree=dl->TreeH();
    hitTree->SetBranchAddress("ITS",&hitList);
    // 
    // read clusters
    for (int ilr=nlr;ilr--;) {
      TBranch* br = cluTree->GetBranch(Form("ITSRecPoints%d",ilr));
      if (!br) {printf("Did not find cluster branch for lr %d\n",ilr); exit(1);}
      br->SetAddress(its->GetLayerActive(ilr)->GetClustersAddress());
    }
    cluTree->GetEntry(0); 
    its->ProcessClusters();
    //
    // read hits
    for(Int_t iEnt=0;iEnt<hitTree->GetEntries();iEnt++){//entries loop degli hits
      hitTree->GetEntry(iEnt);
      int nh = hitList->GetEntries();      
      for(Int_t iHit=0; iHit<nh;iHit++){
        AliITSUHit *pHit = (AliITSUHit*)hitList->At(iHit);
        int mcID = pHit->GetTrack();
        TClonesArray* harr = arrMCTracks.GetEntriesFast()>mcID ? (TClonesArray*)arrMCTracks.At(mcID) : 0;
        if (!harr) {
          harr = new TClonesArray("AliITSUHit"); // 1st encounter of the MC track
          arrMCTracks.AddAtAndExpand(harr,mcID);
        }
        //
        new ( (*harr)[harr->GetEntriesFast()] ) AliITSUHit(*pHit);
      }
    }

    //
    // compare clusters and hits
    //
    printf(" tree entries: %lld\n",cluTree->GetEntries());
    //
    for (int ilr=0;ilr<nlr;ilr++) {
      AliITSURecoLayer* lr = its->GetLayerActive(ilr);
      TClonesArray* clr = lr->GetClusters();
      int nClu = clr->GetEntries();
      printf("Layer %d : %d clusters\n",ilr,nClu);
      //
      for (int icl=0;icl<nClu;icl++) {
        AliITSUClusterPix *cl = (AliITSUClusterPix*)clr->At(icl);
        int modID = cl->GetVolumeId();

        //------------ check if this is a split cluster
        int sInL = modID - gm->GetFirstChipIndex(ilr);
        if (!cl->TestBit(kSplCheck)) {
          cl->SetBit(kSplCheck);
          // check if there is no other cluster with same label on this module
          AliITSURecoSens* sens = lr->GetSensor(sInL);
          int nclSn = sens->GetNClusters();
          int offs = sens->GetFirstClusterId();
          //  printf("To check for %d (mod:%d) N=%d from %d\n",icl,modID,nclSn,offs);
          for (int ics=0;ics<nclSn;ics++) {
            AliITSUClusterPix* clusT = (AliITSUClusterPix*)lr->GetCluster(offs+ics); // access to clusters
            if (clusT==cl) continue;
            for (int ilb0=0;ilb0<3;ilb0++) {
              int lb0 = cl->GetLabel(ilb0); if (lb0<=-1) break;
              for (int ilb1=0;ilb1<3;ilb1++) {
                int lb1 = clusT->GetLabel(ilb1); if (lb1<=-1) break;
                if (lb1==lb0) {
                  cl->SetBit(kSplit);
                  clusT->SetBit(kSplit);
                  /*
                  printf("Discard clusters of module %d:\n",modID);
                  cl->Print();
                  clusT->Print();
                  */
                  break;
                }
              }
            }
          }
        }
        //------------
        const AliITSMFTSegmentationPix* segm = gm->GetSegmentation(ilr);
        //
        cl->GetGlobalXYZ(xyzClGloF);
        int clsize = cl->GetNPix();
        for (int i=3;i--;) xyzClGlo[i] = xyzClGloF[i];
        const TGeoHMatrix* mat = gm->GetMatrixSens(modID);
        if (!mat) {printf("failed to get matrix for module %d\n",cl->GetVolumeId());}
        mat->MasterToLocal(xyzClGlo,xyzClTr);
        //
        int col,row;
        segm->LocalToDet(xyzClTr[0],xyzClTr[2],row,col); // effective col/row
        nLab = 0;
        for (int il=0;il<3;il++) {
          if (cl->GetLabel(il)>=0) labels[nLab++] = cl->GetLabel(il);
          else break;
        }
        // find hit info
        for (int il=0;il<nLab;il++) {
          TClonesArray* htArr = (TClonesArray*)arrMCTracks.At(labels[il]);
          if (!htArr) {printf("did not find MChits for label %d ",labels[il]); cl->Print(); continue;}
          //
          int nh = htArr->GetEntriesFast();
          AliITSUHit *pHit=0;
          for (int ih=nh;ih--;) {
            AliITSUHit* tHit = (AliITSUHit*)htArr->At(ih);
            if (tHit->GetChip()!=modID) continue;
            pHit = tHit;
            break;
          }
          if (!pHit) {
            printf("did not find MChit for label %d on module %d ",il,modID); 
            cl->Print(); 
            htArr->Print();
            continue;
          }
          //
          pHit->GetPositionG(xg1,yg1,zg1);
          pHit->GetPositionG0(xg0,yg0,zg0,tg0);
          //
          double txyzH[3],gxyzH[3] = { (xg1+xg0)/2, (yg1+yg0)/2, (zg1+zg0)/2 };
          mat->MasterToLocal(gxyzH,txyzH);

          double rcl = TMath::Sqrt(xyzClTr[0]*xyzClTr[0]+xyzClTr[1]*xyzClTr[1]);
          double rht = TMath::Sqrt(txyzH[0]*txyzH[0]+txyzH[1]*txyzH[1]);
          //
          //Angles determination

          pHit->GetPositionL(xExit,yExit,zExit,gm);
          pHit->GetPositionL0(xEnt,yEnt,zEnt,tof1,gm);

          Double_t dirHit[3]={(xExit-xEnt),(yExit-yEnt),(zExit-zEnt)};

          /*double PG[3] = {(double)pHit->GetPXG(), (double)pHit->GetPYG(), (double)pHit->GetPZG()}; //Momentum at hit-point in Global Frame
          double PL[3];
          if (TMath::Abs(PG[0])<10e-7 && TMath::Abs(PG[1])<10e-7) {
            pHit->Dump();
            int lb = pHit->GetTrack();
            stack->Particle(lb)->Print();
            continue;
          }
          mat->MasterToLocalVect(PG,PL); //Momentum in local Frame
          //printf(">> %e %e   %e %e   %e %e\n",PG[0],PL[0],PG[1],PL[1],PG[2],PL[2]);*/

          Double_t alpha1 = TMath::ACos(TMath::Abs(dirHit[1])/TMath::Sqrt(dirHit[0]*dirHit[0]+dirHit[1]*dirHit[1]+dirHit[2]*dirHit[2])); //Polar Angle
          Float_t alpha2 = (Float_t) alpha1; //convert to float
          cSum.alpha = alpha2;

          Double_t beta1;
          beta1 = TMath::ATan2(dirHit[0],dirHit[2]); //Azimuthal angle, values from -Pi to Pi
          Float_t beta2 = (Float_t) beta1;
          cSum.beta = beta2;
          
          GetHistoClSize(clsize,kDR,&histoArr)->Fill((rht-rcl)*1e4);
          if (cl->TestBit(kSplit)) {
            if (col%2) GetHistoClSize(clsize,kDTXoddSPL,&histoArr)->Fill((txyzH[0]-xyzClTr[0])*1e4);
            else       GetHistoClSize(clsize,kDTXevenSPL,&histoArr)->Fill((txyzH[0]-xyzClTr[0])*1e4);
            GetHistoClSize(clsize,kDTZSPL,&histoArr)->Fill((txyzH[2]-xyzClTr[2])*1e4);
            GetHistoClSize(0,kNPixSPL,&histoArr)->Fill(clsize);
          }
          if (col%2) GetHistoClSize(clsize,kDTXodd,&histoArr)->Fill((txyzH[0]-xyzClTr[0])*1e4);
          else       GetHistoClSize(clsize,kDTXeven,&histoArr)->Fill((txyzH[0]-xyzClTr[0])*1e4);
          GetHistoClSize(clsize,kDTZ,&histoArr)->Fill((txyzH[2]-xyzClTr[2])*1e4);
          GetHistoClSize(0,kNPixAll,&histoArr)->Fill(clsize);
          //
          cSum.evID = iEvent;
          cSum.volID = cl->GetVolumeId();
          cSum.lrID = ilr;
          cSum.clID = icl;
          cSum.nPix = cl->GetNPix();
          cSum.nX   = cl->GetNx();
          cSum.nZ   = cl->GetNz();
          cSum.q    = cl->GetQ();
          cSum.split = cl->TestBit(kSplit);
          cSum.dX = (txyzH[0]-xyzClTr[0])*1e4;
          cSum.dY = (txyzH[1]-xyzClTr[1])*1e4;
          cSum.dZ = (txyzH[2]-xyzClTr[2])*1e4;
          cSum.nRowPatt = cl-> GetPatternRowSpan();
          cSum.nColPatt = cl-> GetPatternColSpan();
                    
          int label = cl->GetLabel(0);
          TParticle* part = 0;
          if (label>=0 && (part=stack->Particle(label)) ) {
            cSum.pdg = part->GetPdgCode();
            cSum.eta = part->Eta();
            cSum.pt  = part->Pt();
            cSum.phi = part->Phi();
            cSum.prim = stack->IsPhysicalPrimary(label);
          } 
          cSum.ntr = 0;
          for (int ilb=0;ilb<3;ilb++) if (cl->GetLabel(ilb)>=0) cSum.ntr++;
          for (int i=0;i<3;i++) cSum.xyz[i] = xyzClGloF[i];
          //
          trOut->Fill();
          /*
          if (clsize==5) {
            printf("\nL%d(%c) Mod%d, Cl:%d | %+5.1f %+5.1f (%d/%d)|H:%e %e %e | C:%e %e %e\n",ilr,cl->TestBit(kSplit) ? 'S':'N',
             modID,icl,(txyzH[0]-xyzClTr[0])*1e4,(txyzH[2]-xyzClTr[2])*1e4, row,col,
             gxyzH[0],gxyzH[1],gxyzH[2],xyzClGlo[0],xyzClGlo[1],xyzClGlo[2]);
            cl->Print();
            pHit->Print();
            //
            double a0,b0,c0,a1,b1,c1,e0;
            pHit->GetPositionL0(a0,b0,c0,e0,gm);
            pHit->GetPositionL(a1,b1,c1,gm);
            float cloc[3];
            cl->GetLocalXYZ(cloc);
            printf("LocH: %e %e %e | %e %e %e\n",a0,b0,c0,a1,b1,c1);
            printf("LocC: %e %e %e | %e %e %e\n",cloc[0],cloc[1],cloc[2],xyzClTr[0],xyzClTr[1],xyzClTr[2]);
          }
          */
          //
        }
      }
    }
    
    //    layerClus.Clear();
    //
    arrMCTracks.Delete();
  }//event loop
  //
  flOut->cd();
  trOut->Write();
  delete trOut;
  flOut->Close();
  flOut->Delete();
  DrawReport("clinfo.ps",&histoArr);
  //
}

void DrawReport(const char* psname, TObjArray* harr) 
{
  gStyle->SetOptFit(1);
  if (!harr) harr = &histoArr;
  TCanvas* cnv = new TCanvas("cl","cl",900,600);
  //
  TString psnm1 = psname;
  if (psnm1.IsNull()) psnm1 = "clusters.ps";
  TString psnm0 = psnm1.Data(); 
  psnm0 += "[";
  TString psnm2 = psnm1.Data(); 
  psnm2 += "]";
  cnv->Print(psnm0.Data());
  //
  TH1* clall = GetHistoClSize(0,kNPixAll,harr);
  clall->SetLineColor(kRed);
  clall->Draw();
  TH1* clSpl = GetHistoClSize(0,kNPixSPL,harr);
  clSpl->SetLineColor(kBlue);
  clSpl->Draw("sames");
  gPad->Modified();
  gPad->Update();
  SetStPadPos(clall,0.75,0.97,0.8,1.,-1,clall->GetLineColor());
  SetStPadPos(clSpl,0.75,0.97,0.6,0.8,-1,clSpl->GetLineColor());
  gPad->Modified();
  gPad->Update();
  gPad->SetLogy(1);
  //
  cnv->cd();
  cnv->Print(psnm1.Data());
  //
  // plot cluster sized from 1 to 10
  for (int i=1;i<=10;i++) {
    if (clall->GetBinContent(clall->FindBin(i))<100) continue;
    DrawNP(i,harr,cnv);
    cnv->Print(psnm1.Data());
  }
  cnv->Print(psnm2.Data());
}

//------------------------------------------
TH1* GetHistoClSize(int npix,int id,TObjArray* harr)
{
  // book histos
  TH1* h = 0;
  if (!harr) harr = &histoArr;
  //
  if (npix<1) {
    if (harr->GetEntriesFast()>=id && (h=(TH1*)harr->At(id))) return h;
    h = new TH1F("npixAll","npixAll",150,0.5,54.5); 
    h->SetDirectory(0);
    h->SetLineColor(kRed);
    harr->AddAtAndExpand(h, kNPixAll);
    //
    h = new TH1F("npixSpl","npixSpl",150,0.5,54.5);
    h->SetLineColor(kBlue);
    h->SetDirectory(0);
    harr->AddAtAndExpand(h, kNPixSPL);
    //
    h = (TH1*)harr->At(id);
    if (!h) {printf("Unknown histo id=%d\n",id); exit(1);}
    return h;
  }
  //
  int idh = npix*10+id;
  if (harr->GetEntriesFast()>=idh && (h=(TH1*)harr->At(idh))) return h;
  //
  const int nbin=100;
  const double kdiff=80;
  // need to create set of histos
  //
  h = new TH1F(Form("dxy_npix%d",npix),Form("dr_npix%d",npix),nbin,-kdiff,kdiff);
  h->SetDirectory(0);
  harr->AddAtAndExpand(h, npix*10 + kDR);
  //
  h  = new TH1F(Form("dtxODD_npix%d",npix),Form("dtxODD_npix%d",npix),nbin,-kdiff,kdiff);
  h->SetDirectory(0);
  h->SetLineColor(kRed);
  harr->AddAtAndExpand(h, npix*10 + kDTXodd);
  h  = new TH1F(Form("dtxEVN_npix%d",npix),Form("dtxEVN_npix%d",npix),nbin,-kdiff,kdiff);
  h->SetDirectory(0);
  h->SetLineColor(kBlue);
  harr->AddAtAndExpand(h, npix*10 + kDTXeven);
  //
  h  = new TH1F(Form("dtz_npix%d",npix),Form("dtz_npix%d",npix),nbin,-kdiff,kdiff);
  h->SetLineColor(kGreen);
  h->SetDirectory(0);
  harr->AddAtAndExpand(h, npix*10 + kDTZ);
  //
  //
  h  = new TH1F(Form("SPL_dtxODD_npix%d",npix),Form("SPL_dtxODD_npix%d",npix),nbin,-kdiff,kdiff);
  h->SetLineColor(kMagenta);
  h->SetFillColor(kMagenta);
  h->SetFillStyle(3001);  
  h->SetLineStyle(2);
  h->SetDirectory(0);

  harr->AddAtAndExpand(h, npix*10 + kDTXoddSPL);
  h  = new TH1F(Form("SPL_dtxEVN_npix%d",npix),Form("SPL_dtxEVN_npix%d",npix),nbin,-kdiff,kdiff);
  h->SetLineColor(kCyan);
  h->SetFillColor(kCyan);
  h->SetFillStyle(3006);  
  h->SetLineStyle(2);
  h->SetDirectory(0);
  harr->AddAtAndExpand(h, npix*10 + kDTXevenSPL);
  //
  h  = new TH1F(Form("SPL_dtz_npix%d",npix),Form("SPLdtz_npix%d",npix),nbin,-kdiff,kdiff);
  harr->AddAtAndExpand(h, npix*10 + kDTZSPL);
  h->SetDirectory(0);
  //
  h->SetLineColor(kGreen+2);
  h->SetFillColor(kGreen+2);
  h->SetLineStyle(2);
  h->SetFillStyle(3001);
  h = (TH1*)harr->At(idh);
  if (!h) {printf("Unknown histo id=%d\n",idh); exit(1);}
  return h;
}



TCanvas* DrawNP(int np, TObjArray* harr, TCanvas* cnv)
{
  if (!harr) harr = &histoArr;
  if (!cnv) cnv = new TCanvas(Form("cnv%d",np),Form("cnv%d",np),900,700);
  cnv->Clear();
  cnv->Divide(2,1);
  cnv->cd(1);
  //
  TH1* dxodd = (TH1*)harr->At(np*10+kDTXodd);
  TH1* dxevn = (TH1*)harr->At(np*10+kDTXeven);
  TH1* dxoddS =(TH1*)harr->At(np*10+kDTXoddSPL);
  TH1* dxevnS =(TH1*)harr->At(np*10+kDTXevenSPL);
  double max = TMath::Max(dxodd->GetMaximum(),dxevn->GetMaximum());
  dxodd->SetMaximum(1.1*max);
  dxodd->GetXaxis()->SetTitle("#DeltaX, #mum");
  dxodd->SetTitle(Form("#DeltaX for clSize=%d",np));
  dxodd->Fit("gaus","","");
  dxevn->Fit("gaus","","sames");
  //
  dxoddS->Draw("sames");
  dxevnS->Draw("sames");
  //
  gPad->Modified();
  gPad->Update();
  SetStPadPos(dxodd,0.75,0.97,0.8,1., -1,dxodd->GetLineColor());
  SetStPadPos(dxevn,0.75,0.97,0.6,0.8, -1,dxevn->GetLineColor());
  SetStPadPos(dxoddS,0.75,0.97,0.4,0.6, -1,dxoddS->GetLineColor());
  SetStPadPos(dxevnS,0.75,0.97,0.2,0.4, -1,dxevnS->GetLineColor());
  //
  cnv->cd(2);
  TH1* dz  = (TH1*)harr->At(np*10+kDTZ);
  dz->SetTitle(Form("#DeltaZ for clSize=%d",np));
  dz->GetXaxis()->SetTitle("#DeltaZ, #mum");
  dz->Fit("gaus");
  TH1* dzS = (TH1*)harr->At(np*10+kDTZSPL);
  dz->Draw("sames");
  gPad->Modified();
  gPad->Update();
  SetStPadPos(dz,0.75,0.97,0.8,1., -1, dz->GetLineColor());
  SetStPadPos(dzS,0.75,0.97,0.5,0.7, -1, dzS->GetLineColor());
  gPad->Modified();
  gPad->Update();
  //
  cnv->cd();
  return cnv;
}


TPaveStats* GetStPad(TH1* hst)
{
  TList *lst = hst->GetListOfFunctions();
  if (!lst) return 0;
  int nf = lst->GetSize();
  for (int i=0;i<nf;i++) {
    TPaveStats *fnc = (TPaveStats*) lst->At(i);
    if (fnc->InheritsFrom("TPaveStats")) return fnc;
  }
  return 0;
  //
}


TPaveStats* SetStPadPos(TH1* hst,float x1,float x2,float y1,float y2, Int_t stl, Int_t col)
{
  TPaveStats* pad = GetStPad(hst);
  if (!pad) return 0;
  pad->SetX1NDC( x1 );
  pad->SetX2NDC( x2 );
  pad->SetY1NDC( y1 );
  pad->SetY2NDC( y2 );
  if (stl>=0) pad->SetFillStyle(stl);
  if (col>=0) pad->SetTextColor(col);
  pad->SetFillColor(0);
  //
  gPad->Modified();
  return pad;
}
