#if !defined(__CINT__) || defined(__MAKECINT__)
#include<Riostream.h>
#include<TROOT.h>
#include<TArrayI.h>
#include<TBranch.h>
#include<TCanvas.h>
#include<TClassTable.h>
#include<TClonesArray.h>
#include<TFile.h>
#include<TGeoManager.h>
#include<TH1.h>
#include<TH2.h>
#include <TInterpreter.h>
#include<TObject.h>
#include<TObjArray.h>
#include<TTree.h>
#include "AliRun.h"
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSDetTypeRec.h"
#include "AliITSRecPoint.h"
#include "AliITSRecPoint.h"
#include "AliITSdigit.h"
#include "AliITSdigitSSD.h"
#include "AliITShit.h"
#include "AliITSmodule.h" 
#include "AliITSsegmentation.h"
#include "AliITSsegmentationSPD.h" 
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliRunLoader.h"
#include "AliITSLoader.h"
#include "AliHeader.h"
#endif
void GetHitsCoor(TObject *its, Int_t mod, TObjArray & histos, Int_t subd,Bool_t verb);
Int_t GetRecCoor(TObject *ge, TClonesArray *ITSrec, Int_t mod, TH2F *h2, TH1F *h1, Bool_t verb);
Int_t GetClusCoor(TObject *ge, TClonesArray *ITSrec, Int_t mod, TH2F *h2, TH1F *h1, Bool_t verb);
void GetDigits(TObject *tmps,TObject *ge,TClonesArray *ITSdigits, Int_t subd, Int_t mod, Bool_t verbose, TObjArray & histos);

Int_t AliITSGeoPlot (Int_t evesel=0, char *opt="All+ClustersV2", TString filename="galice.root", Int_t isfastpoints = 0) {
  /*******************************************************************
   *  This macro displays geometrical information related to the
   *  hits, digits and rec points (or V2 clusters) in ITS.
   *  There are histograms that are not displayed (i.e. energy
   *  deposition) but that are saved onto a file (see below)
   *
   *  INPUT arguments:
   *
   *  Int_t evesel:  the selected event number. It's =0 by default
   *
   *  Options: Any combination of:
   *    1) subdetector name:  "SPD", "SDD", "SSD", "All" (default)
   *    2) Printouts:        "Verbose" Almost everything is printed out: 
   *                          it is wise to redirect the output onto a file
   *                    e.g.: .x AliITSGeoPlot.C("All+Verbose") > out.log 
   *    3) Rec Points option: "Rec"   ---> Uses Rec. Points (default)
   * 
   *    4) ClustersV2 option: "ClustersV2" ---->Uses ClustersV2
   *                           otherwise ---->uses hits and digits only
   *    Examples:
   *       .x AliITSGeoPlot();  (All subdetectors; no-verbose; no-recpoints)
   *       .x AliITSGeoPlot("SPD+SSD+Verbose+Rec"); 
   *   
   *    filename:   It's "galice.root" by default. 
   *    isfastpoints: integer. It is set to 0 by defaults. This means that
   *                slow reconstruction is assumed. If fast recpoint are in
   *                in use, isfastpoints must be set =1.
   *
   *  OUTPUT: It produces a root file with a list of histograms
   *    
   *  WARNING: spatial information for SSD/DIGITS is obtained by pairing
   *          digits on p and n side originating from the same track, when
   *          possible. This (mis)use of DIGITS is tolerated for debugging 
   *          purposes only !!!!  The pairing in real life should be done
   *          starting from info really available... 
   * 
   *  COMPILATION: this macro can be compiled. 
   *      1)       You need to set your include path with
   * gSystem->SetIncludePath("-I- -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -g");
   *      3)       If you are using root instead of aliroot you need to
   *               execute the macro loadlibs.C in advance
   *      4)       To compile this macro from root (or aliroot):
   *                 ---  .L AliITSGeoPlot.C++
   *                 ---  AliITSGeoPlot();
   *     
   *  M.Masera  14/05/2001 18:30
   *  Last rev. 31/05/2004 14:00 (Clusters V2 added)  E.C.          
   ********************************************************************/

  //Options
  TString choice(opt);
  Bool_t All = choice.Contains("All");
  Bool_t verbose=choice.Contains("Verbose");
  Bool_t userec=choice.Contains("Rec");
  Bool_t useclustersv2=choice.Contains("ClustersV2");
  Int_t retcode=1; //return code
 
  if (gClassTable->GetID("AliRun") < 0) {
    gInterpreter->ExecuteMacro("loadlibs.C");
  }
  else { 
    if(gAlice){
      delete gAlice->GetRunLoader();
      delete gAlice;
      gAlice=0;
    }
  }
  // retrives geometry 
  TString geof(gSystem->DirName(filename));
  geof += "/geometry.root";
  TGeoManager::Import(geof.Data());
  if (!gGeoManager) {
    cout<<"geometry not found\n";
    return -1;
  }
 
  AliRunLoader* rl = AliRunLoader::Open(filename.Data());
  if (rl == 0x0){
    cerr<<"AliITSGeoPlot.C : Can not open session RL=NULL"<< endl;
    return -1;
  }
  Int_t retval = rl->LoadgAlice();
  if (retval){
    cerr<<"AliITSGeoPlot.C : LoadgAlice returned error"<<endl;
    return -1;
  }
  gAlice=rl->GetAliRun();

  retval = rl->LoadHeader();
  if (retval){
    cerr<<"AliITSGeoPlot.C : LoadHeader returned error"<<endl;
    return -1;
  }

  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");

  if(!ITSloader){
    cerr<<"AliITSGeoPlot.C :  ITS loader not found"<<endl;
    return -1;
  }

  ITSloader->LoadHits("read");
  ITSloader->LoadDigits("read");
  if(isfastpoints==1)ITSloader->SetRecPointsFileName("ITS.FastRecPoints.root");
  ITSloader->LoadRecPoints("read");
  rl->GetEvent(evesel);
  Int_t nparticles = rl->GetHeader()->GetNtrack();
  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
  ITS->SetTreeAddress();
  if(verbose) {
    cout<<" "<<endl<<" "<<endl;
    cout<<"******* Event processing started   *******"<<endl;
    cout<<"In the following, hits with 'StatusEntering' flag active"<<endl;
    cout<<"will not be processed"<<endl; 
    cout << "Number of particles=  " << nparticles <<endl;
  }

  // HITS
  TTree *TH = ITSloader->TreeH();
  Stat_t ntracks = TH->GetEntries();
  if(verbose)cout<<"Number of primary tracks= "<<ntracks<<endl;

  // ITS
  Int_t nmodules;
  ITS->InitModules(-1,nmodules);
  cout<<"Number of ITS modules= "<<nmodules<<endl;
  cout<<"Filling modules... It takes a while, now. Please be patient"<<endl;
  ITS->FillModules(0,0,nmodules," "," ");
  cout<<"ITS modules .... DONE!"<<endl;
  
  AliITSDetTypeRec* detTypeRec = new AliITSDetTypeRec(ITSloader);
  detTypeRec->SetDefaults();
  // DIGITS
  TTree *TD = ITSloader->TreeD();

  //RECPOINTS (V2 clusters)
  TTree *TR = ITSloader->TreeR();
  TClonesArray *ITSrec  = detTypeRec->RecPoints();
  TBranch *branch = 0;
  if(userec && TR && ITSrec){
    if(isfastpoints==1){
      branch = ITSloader->TreeR()->GetBranch("ITSRecPointsF");
      cout<<"using fast points\n";
    }
    else {
      branch = ITSloader->TreeR()->GetBranch("ITSRecPoints");
    }
    if(branch)branch->SetAddress(&ITSrec);
  }

  if(userec && (!TR || !ITSrec || !branch)){
    userec = kFALSE;
    cout<<"\n ======================================================= \n";
    cout<<"WARNING: there are no RECPOINTS on this file ! \n";
    cout<<"======================================================= \n \n";
  }
  if(useclustersv2 && TR && ITSrec){
    branch = ITSloader->TreeR()->GetBranch("ITSRecPoints");
    if(branch)branch->SetAddress(&ITSrec);
  }

  if(useclustersv2 && (!TR || !ITSrec || !branch)){
    useclustersv2 = kFALSE;
    cout<<"\n ======================================================= \n";
    cout<<"WARNING: there are no CLUSTERSV2 on this file ! \n";
    cout<<"======================================================= \n \n";
  }


  //local variables
  Int_t mod;   //module number
  Int_t nbytes = 0; 
  Double_t pos[3];  // Global position of the current module
  Float_t ragdet; // Radius of detector (x y plane)
  Int_t first; // first module
  Int_t last; // last module
  Int_t nrecp; //number of RecPoints for a given module

  //List of histograms
  TObjArray histos(26,0);  // contains the pointers to the histograms
  // Book histograms SPD
  TH2F *rspd = new TH2F("rspd","Radii of digits - SPD",50,-10.,10.,50,-10.,10.);
  TH2F *rhspd = new TH2F("rhspd","Radii of hits - SPD",50,-10.,10.,50,-10.,10.);
  TH2F *rmspd = new TH2F("rmspd","Radii of SPD modules",50,-10.,10.,50,-10.,10.);
  TH1F *zspd = new TH1F("zspd","Z of digits - SPD",100,-30.,30.);
  TH1F *zhspd = new TH1F("zhspd","Z of hits - SPD",100,-30.,30.);
  TH1F *zmspd = new TH1F("zmspd","Z of SPD modules",100,-30,30.);

  Char_t title1[50]="";
  Char_t title2[50]="";
  if(userec){ 
    sprintf(title1,"Radii of recpoints - %s","SPD");
    sprintf(title2,"Z of recpoints - %s","SPD");
  }
  if(useclustersv2){
    sprintf(title1,"Radii of clustersV2 - %s","SPD");
    sprintf(title2,"Z of clustersV2 - %s","SPD");
  }
  TH2F *rrspd = new TH2F("rrspd",title1,50,-10.,10.,50,-10.,10.);
  TH1F *zrspd = new TH1F("zrspd",title2,100,-30.,30.);
  TH1F *enespd = new TH1F("enespd","Energy deposition SPD (KeV)",100,0.,1000.);
  histos.AddLast(rspd);  // 0
  histos.AddLast(rhspd); // 1
  histos.AddLast(rmspd); // 2
  histos.AddLast(zspd);  // 3
  histos.AddLast(zhspd); // 4
  histos.AddLast(zmspd); // 5
  histos.AddLast(rrspd); // 6
  histos.AddLast(zrspd); // 7
  histos.AddLast(enespd); // 8
  // Book histograms SDD
  TH2F *rsdd = new TH2F("rsdd","Radii of digits - SDD",50,-40.,40.,50,-40.,40.);
  TH2F *rhsdd = new TH2F("rhsdd","Radii of hits - SDD",50,-40.,40.,50,-40.,40.);
  TH2F *rmsdd = new TH2F("rmsdd","Radii of SDD modules",50,-40.,40.,50,-40.,40.);
  TH1F *zsdd = new TH1F("zsdd","Z of digits - SDD",100,-40.,40.);
  TH1F *zhsdd = new TH1F("zhsdd","Z of hits - SDD",100,-40.,40.);
  TH1F *zmsdd = new TH1F("zmsdd","Z of SDD modules",100,-40,40.);
  Char_t title3[50];
  Char_t title4[50];
  if(userec){ 
    sprintf(title3,"Radii of recpoints - %s","SDD");
    sprintf(title4,"Z of recpoints - %s","SDD");
  }
  if(useclustersv2){
    sprintf(title3,"Radii of clustersV2 - %s","SDD");
    sprintf(title4,"Z of clustersV2 - %s","SDD");
  }
  TH2F *rrsdd = new TH2F("rrsdd",title3,50,-40.,40.,50,-40.,40.);   
  TH1F *zrsdd = new TH1F("zrsdd",title4,100,-40.,40.);
  TH1F *enesdd = new TH1F("enesdd","Energy deposition SDD (KeV)",100,0.,1000.);
  histos.AddLast(rsdd);  // 9
  histos.AddLast(rhsdd); // 10
  histos.AddLast(rmsdd); // 11
  histos.AddLast(zsdd);  // 12
  histos.AddLast(zhsdd); // 13
  histos.AddLast(zmsdd); // 14
  histos.AddLast(rrsdd); // 15
  histos.AddLast(zrsdd); // 16
  histos.AddLast(enesdd); // 17
  // Book histogram SSD
  TH2F *rssd = new TH2F("rssd","Radii of digits - SSD",50,-50.,50.,50,-50.,50.);
  TH2F *rhssd = new TH2F("rhssd","Radii of hits - SSD",50,-50.,50.,50,-50.,50.);
  TH2F *rmssd = new TH2F("rmssd","Radii of SSD modules",50,-50.,50.,50,-50.,50.);
  TH1F *zssd = new TH1F("zssd","Z of digits - SSD",100,-70.,70.);
  TH1F *zhssd = new TH1F("zhssd","Z of hits - SSD",100,-70.,70.);
  TH1F *zmssd = new TH1F("zmssd","Z of SSD modules",100,-70,70.);
  Char_t title5[50];
  Char_t title6[50];
  if(userec){ 
    sprintf(title5,"Radii of recpoints - %s","SSD");
    sprintf(title6,"Z of recpoints - %s","SSD");
  }
  if(useclustersv2){
    sprintf(title5,"Radii of clustersV2 - %s","SSD");
    sprintf(title6,"Z of clustersV2 - %s","SSD");
  }

  TH2F *rrssd = new TH2F("rrssd",title5,50,-50.,50.,50,-50.,50.);
  TH1F *zrssd = new TH1F("zrssd",title6,100,-70.,70.);
  TH1F *enessd = new TH1F("enessd","Energy deposition SSD (KeV)",100,0.,1000.);
  histos.AddLast(rssd);  // 18
  histos.AddLast(rhssd); // 19
  histos.AddLast(rmssd); // 20
  histos.AddLast(zssd);  // 21
  histos.AddLast(zhssd); // 22
  histos.AddLast(zmssd); // 23
  histos.AddLast(rrssd); // 24
  histos.AddLast(zrssd); // 25
  histos.AddLast(enessd); // 26
  //
  // Loop on subdetectors
  // 
  cout<<"CALL GETITSGEOM \n \n \n";
  AliITSgeom *geom = ITS->GetITSgeom();
  TString detna; // subdetector name
  for(Int_t subd=0;subd<3;subd++){
    if(All || (choice.Contains("SPD") && subd==0) || (choice.Contains("SDD") && subd==1) || (choice.Contains("SSD") && subd==2)){
      // Prepare array for the digits
      TClonesArray *ITSdigits  = ITS->DigitsAddress(subd);
      Bool_t usedigits = kTRUE;
      if(!ITSdigits){
        usedigits = kFALSE;
        cout<<"\n ======================================================= \n";
        cout<<"WARNING: there are no DIGITS on this file ! \n";
        cout<<"======================================================= \n \n";
      }
      // Get segmentation model
      if(subd==0)detna="SPD";
      if(subd==1)detna="SDD";
      if(subd==2)detna="SSD";
      AliITSsegmentation *seg=(AliITSsegmentation*)detTypeRec->GetSegmentationModel(subd);
      // Loop on modules
      first = geom->GetStartDet(subd);
      last = geom->GetLastDet(subd);
      if(verbose){
        cout<<"     "<<endl<<"-------------------------------------"<<endl;
        cout<<"Start processing subdetector "<<detna<<endl;
        cout<<detna<<" first module "<<first<<endl;
        cout<<detna<<" last module "<<last<<endl;
        cout<<" "<<endl<<" "<<endl;
      }
      for (mod=first; mod<=last; mod++){
        geom->GetTrans(mod,pos);  // position of the module in the MRS
        ragdet=sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
        // The following 2 histos are a check of the geometry
        TH2F *bidi = (TH2F*)histos.At(2+subd*9);
        TH1F *uni = (TH1F*)histos.At(5+subd*9);
        bidi->Fill(pos[0],pos[1]);
        uni->Fill(pos[2]);
        if(verbose){
          cout<<"=========================================================\n";
          cout<<detna<<" module="<<mod<<endl;
          cout<<"Mod. coordinates: "<<pos[0]<<", "<<pos[1]<<", ";
          cout<<pos[2]<<" Radius "<<ragdet<<endl;
        }

        // Hits
        GetHitsCoor(ITS,mod,histos,subd,verbose);

        //RecPoints     
        if(userec){
          detTypeRec->ResetRecPoints();
          branch->GetEvent(mod);
          TH2F *bidi=(TH2F*)histos.At(6+subd*9);
          TH1F *uni=(TH1F*)histos.At(7+subd*9);
          nrecp=GetRecCoor(geom,ITSrec,mod,bidi,uni,verbose);
        }
        if(useclustersv2){
          detTypeRec->ResetRecPoints();
          branch->GetEvent(mod);
          TH2F *bidi=(TH2F*)histos.At(6+subd*9);
          TH1F *uni=(TH1F*)histos.At(7+subd*9);
          nrecp=GetClusCoor(geom,ITSrec,mod,bidi,uni,verbose);
	  
        }
     
        // Digits
        if(usedigits){
          detTypeRec->ResetDigits();
          nbytes += TD->GetEvent(mod);
          GetDigits(seg,geom,ITSdigits,subd,mod,verbose,histos);
        }

      } // End of loop on the modules
      TH1F *h1tmp;
      TH2F *h2tmp;
      //  Plot the histograms
      TCanvas *current=0; // current Canvas (1--> SPD, 2---> SDD, 3---> SSD)
      if(subd==0){
        // Prepare canvas 1
        TCanvas *c1 = new TCanvas("c1","SPD",10,10,600,900);
        c1->Divide(2,3);
        current=c1;
      }
      if(subd==1){
        // Prepare canvas 2
        TCanvas *c2 = new TCanvas("c2","SDD",40,40,600,900);
        c2->Divide(2,3);
        current=c2;
      }
      if(subd==2){
        // Prepare canvas 3
        TCanvas *c3 = new TCanvas("c3","SSD",70,70,600,900);
        c3->Divide(2,3);
        current=c3;
      }
      current->cd(1);
      h2tmp = (TH2F*)histos.At(9*subd);
      h2tmp->Draw();
      current->cd(2);
      h1tmp=(TH1F*)histos.At(3+subd*9);
      h1tmp->Draw();
      current->cd(3);
      h2tmp=(TH2F*)histos.At(1+9*subd);
      h2tmp->Draw();
      current->cd(4);
      h1tmp=(TH1F*)histos.At(4+subd*9);
      h1tmp->Draw();
   
      if(userec || useclustersv2){
        current->cd(5);
        h2tmp=(TH2F*)histos.At(6+9*subd);
        h2tmp->Draw();
        current->cd(6);
        h1tmp=(TH1F*)histos.At(7+subd*9);
        h1tmp->Draw();
      }
  
      else {
        current->cd(5);
        h2tmp=(TH2F*)histos.At(2+9*subd);
        h2tmp->Draw();
        current->cd(6);
        h2tmp=(TH2F*)histos.At(5+9*subd);
        h2tmp->Draw();
      }
    } // if(All.....
  } // end of loop on subdetectors
  // Save the histograms
  TFile *fh = new TFile("AliITSGeoPlot.root","recreate");
  // The list is written to file as a single entry
  TList *lihi = new TList();
  // copy the pointers to the histograms to a TList object.
  // The histograms concerning recpoints are not copied if
  // 'userec' is false.
  for(Int_t i=0;i<histos.GetEntriesFast();i++){
    if(choice.Contains("All") || (choice.Contains("SPD") && i<8) || (choice.Contains("SDD") && i>7 && i<16) || (choice.Contains("SSD") && i>15)){
      if(!(!userec && ((i+2)%9==0 || (i+1)%9==0)))lihi->Add(histos.At(i));
    }
  }
  lihi->Write("Histograms ITS hits+digits+recpoints",TObject::kSingleKey);
  fh->Close();

  return retcode;
}


void GetHitsCoor(TObject *its, Int_t mod, TObjArray & histos, Int_t subd,Bool_t verb){
  TH2F *h2=(TH2F*)histos.At(1+subd*9);
  TH1F *h1=(TH1F*)histos.At(4+subd*9);
  TH1F *ener=(TH1F*)histos.At(8+subd*9);
  AliITS *ITS= (AliITS*)its;
  AliITSmodule *modu = ITS->GetModule(mod);
  TObjArray *fHits = modu->GetHits();
  Int_t nhits = fHits->GetEntriesFast();
  if(nhits>0){
    if(verb){
      cout<<"-------------------------------------------------------"<<endl;
      cout<<"Number of HITS for module "<<mod<<": "<<nhits<<endl;
    }
    for (Int_t hit=0;hit<nhits;hit++) {
      AliITShit *iHit = (AliITShit*) fHits->At(hit);
      if(!iHit->StatusEntering()){
        Float_t x=iHit->GetXG();
        Float_t y=iHit->GetYG();
        Float_t z=iHit->GetZG();
        Float_t edep=iHit->GetIonization()*1000000;
        h2->Fill(x,y);
        h1->Fill(z);
        ener->Fill(edep);
        if(verb){
          cout<<"hit # "<<hit<<" Global coordinates "<<x<<" "<<y<<" "<<z<<endl;
          cout<<"track # "<<iHit->GetTrack()<<" energy deposition (KeV)= ";
          cout<<edep<<endl;
        }
      }
    }
  }
}


Int_t GetClusCoor(TObject *ge, TClonesArray *ITSrec, Int_t mod, TH2F *h2, TH1F *h1, Bool_t verb){

  AliITSgeom *geom = (AliITSgeom*)ge;
  Int_t nrecp = ITSrec->GetEntries();
  if(nrecp>0){
    Float_t lc[3]; for(Int_t i=0; i<3; i++) lc[i]=0.;
    Float_t gc[3]; for(Int_t i=0; i<3; i++) gc[i]=0.;
    if(verb){
      cout<<"-------------------------------------------------------"<<endl;
      cout<<"Number of CLUSTERS for module "<<mod<<": "<<nrecp<<endl;
    }
    for(Int_t irec=0;irec<nrecp;irec++) {
      AliITSRecPoint *recp = (AliITSRecPoint*)ITSrec->At(irec);
      Double_t rot[9];     
      geom->GetRotMatrix(mod,rot);
      Int_t lay,lad,det;   
      geom->GetModuleId(mod,lay,lad,det);
      Float_t tx,ty,tz;    
      geom->GetTrans(lay,lad,det,tx,ty,tz);     

      Double_t alpha=TMath::ATan2(rot[1],rot[0])+TMath::Pi();
      Double_t phi1=TMath::Pi()/2+alpha;
      if(lay==1) phi1+=TMath::Pi();

      Float_t cp=TMath::Cos(phi1), sp=TMath::Sin(phi1);
      Float_t  r=tx*cp+ty*sp;
      gc[0]= r*cp - recp->GetY()*sp;
      gc[1]= r*sp + recp->GetY()*cp;
      gc[2]=recp->GetZ();
  
      if(verb){
        Float_t r=sqrt(gc[0]*gc[0]+gc[1]*gc[1]);
        cout<<"Global coor. + radius "<<gc[0]<<" "<<gc[1]<<" "<<gc[2]<<" ";
        cout<<r<<endl;
        cout<<"Associated track "<<recp->GetLabel(0)<<endl;
      }
      h2->Fill(gc[0],gc[1]);
      h1->Fill(gc[2]);

    }
  }
  return nrecp;
}
Int_t GetRecCoor(TObject *ge, TClonesArray *ITSrec, Int_t mod, TH2F *h2, TH1F *h1, Bool_t verb){
  AliITSgeom *geom = (AliITSgeom*)ge;
  Int_t nrecp = ITSrec->GetEntries();
  if(nrecp>0){
    Float_t lc[3]; for(Int_t i=0; i<3; i++) lc[i]=0.;
    Float_t gc[3]; for(Int_t i=0; i<3; i++) gc[i]=0.;
    if(verb){
      cout<<"-------------------------------------------------------"<<endl;
      cout<<"Number of REC POINTS for module "<<mod<<": "<<nrecp<<endl;
    }
    for(Int_t irec=0;irec<nrecp;irec++) {
      AliITSRecPoint *recp = (AliITSRecPoint*)ITSrec->At(irec);
      lc[0]=recp->GetDetLocalX();
      lc[2]=recp->GetDetLocalZ();
      geom->LtoG(mod,lc,gc);
      if(verb){
        cout<<"recp # "<<irec<<" local coordinates. lx= "<<lc[0]<<" lz= ";
        cout<<lc[2]<<endl;
        Float_t r=sqrt(gc[0]*gc[0]+gc[1]*gc[1]);
        cout<<"Global coor. + radius "<<gc[0]<<" "<<gc[1]<<" "<<gc[2]<<" ";
        cout<<r<<endl;
        cout<<"Associated track "<<recp->GetLabel(0)<<endl;
      }
      h2->Fill(gc[0],gc[1]);
      h1->Fill(gc[2]);
    }
  }
  return nrecp;
}

void GetDigits(TObject *tmps,TObject *ge,TClonesArray *ITSdigits, Int_t subd, Int_t mod, Bool_t verbose, TObjArray & histos){
  AliITSsegmentation *seg = (AliITSsegmentation*)tmps;
  AliITSgeom *geom = (AliITSgeom*)ge;
  Int_t layer;
  Int_t ladder;
  Int_t detec;
  if(subd==2){
    geom->GetModuleId(mod,layer,ladder,detec);
    seg->SetLayer(layer);
  }
  Float_t lcoor[3]; for(Int_t j=0; j<3; j++) lcoor[j]=0.;  //local coord dig.
  Float_t gcoor[3]; for(Int_t j=0; j<3; j++) gcoor[j]=0.; // global coo. dig.
  Float_t ragdig; // Radius digit
  TArrayI ssdone(5000);  // used to match p and n side digits of strips
  TArrayI pair(5000);    // as above 
  Int_t ndigits = ITSdigits->GetEntries();
  AliITSdigit *digs;
  if(ndigits){ 
    if(verbose){
      cout<<"-------------------------------------------------------"<<endl;
      cout<<"Number of DIGITS for module "<<mod<<": "<<ndigits<<endl;
    }
    // Get the coordinates of the module
    if(subd==2){
      for (Int_t digit=0;digit<ndigits;digit++){
        ssdone[digit]=0;
        pair[digit]=0;
      }
    }
    for (Int_t digit=0;digit<ndigits;digit++) {
      digs = (AliITSdigit*)ITSdigits->UncheckedAt(digit);
      Int_t iz=digs->GetCoord1();  // cell number z
      Int_t ix=digs->GetCoord2();  // cell number x
      // Get local coordinates of the element 
      if(subd<2){
        seg->DetToLocal(ix,iz,lcoor[0],lcoor[2]);
      }
      else{
        // SSD: if iz==0 ---> N side; if iz==1 P side
        if(ssdone[digit]==0){
          ssdone[digit]=1;
          pair[digit]=-1;
          Bool_t pside=(iz==1);
          Bool_t impaired=kTRUE;
          Int_t pstrip=0;
          Int_t nstrip=0;
          if(pside)pstrip=ix;
          if(!pside)nstrip=ix;
          for(Int_t digi2=0;digi2<ndigits;digi2++){
            if(ssdone[digi2]==0 && impaired){
              AliITSdigitSSD *dig2=(AliITSdigitSSD*)ITSdigits->UncheckedAt(digi2);
              if(dig2->GetCoord1() != iz && dig2->GetTrack(0)==digs->GetTrack(0) && dig2->GetTrack(0)>0){
                ssdone[digi2]=2;
                pair[digit]=digi2;
                if(pside)nstrip=dig2->GetCoord2();
                if(!pside)pstrip=dig2->GetCoord2();
                impaired=kFALSE;
              }
            }
          }
          if(!impaired)seg->GetPadCxz(pstrip,nstrip,lcoor[0],lcoor[2]);
        }
      }
      if(subd<2 || (subd==2 && ssdone[digit]==1)){
        Int_t coor1=digs->GetCoord1();
        Int_t coor2=digs->GetCoord2();
        Int_t tra0=digs->GetTrack(0);
        if(verbose){
          cout<<"digit # "<<digit<<" fCoord1= "<<coor1<<" fCoord2= ";
          cout<<coor2<<" track "<<tra0<<" "<<endl;
          if(subd<2)cout<<"local coordinates -- x="<<lcoor[0]<<", z=";
          cout<<lcoor[2]<<endl;
          if(subd==2){
            if(pair[digit]==-1){
              cout<<"(digit  not paired)"<<endl;
            }
            else {
              cout<<"local coordinates -- x="<<lcoor[0]<<", z="<<lcoor[2];
              cout<<endl;
              Int_t dtmp=pair[digit];
              AliITSdigitSSD *dig2=(AliITSdigitSSD*)ITSdigits->UncheckedAt(dtmp);
              Int_t coor1b=dig2->GetCoord1();
              Int_t coor2b=dig2->GetCoord2();
              Int_t tra0b=dig2->GetTrack(0);
              cout<<"(digit paired with digit #"<<dtmp<<endl;
              cout<<"with fCoord1= "<<coor1b<<" fCoord2= "<<coor2b;
              cout<<" track "<<tra0b<<")"<<endl;
            }
          }
        }
        if(subd<2 || (subd==2 && pair[digit]!=-1)){
          // Global coordinates of the element
          //SDD and SPD use cm, SSD microns (GetPadCxz)
          if(subd==2)for(Int_t j=0;j<3;j++)lcoor[j]=lcoor[j]/10000.;
          lcoor[1]=0.;
          geom->LtoG(mod,lcoor,gcoor);  // global coord. in cm
          ragdig=sqrt(gcoor[0]*gcoor[0]+gcoor[1]*gcoor[1]);
          if(verbose){
	    cout<<"global coordinates "<<gcoor[0]<<" "<<gcoor[1];
	    cout<<" "<<gcoor[2]<<" Radius "<<ragdig<<endl;
	  }   
          //Fill histograms
          TH2F *bidi = (TH2F*)histos.At(subd*9);
          TH1F *uni = (TH1F*)histos.At(3+subd*9);
          bidi->Fill(gcoor[0],gcoor[1]);
          uni->Fill(gcoor[2]);
        }
      }
    } // loop on digits for this module
  } // if(ndigits>0....
}

