#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TMath.h>
#include <TStopwatch.h>
#include "AliAlignObjParams.h"
#include "AliTrackPointArray.h"
#include "AliITSAlignMille.h"
#endif
//********************************************************************
//  Example to run ITS alignment using Millepede
//
//  Read tracks from AliTrackPoints.N.root and fill Millepede.
//  Millepede configuration is taken from AliITSAlignMille.conf
//  Results are written as AliAlignObjParams in outfile.
//
//      Origin: M. Lunardon
//********************************************************************

Bool_t SelectLayers(AliTrackPointArray *tpa, int *nreqpts_lay) {
  // selection on layers
  int npts=tpa->GetNPoints();
  int nlay[6]; 
  for (int jj=0;jj<6;jj++) nlay[jj]=0;
  for (int ip=0; ip<npts; ip++) {
    int lay=-1;
    float r=TMath::Sqrt(tpa->GetX()[ip]*tpa->GetX()[ip] + tpa->GetY()[ip]*tpa->GetY()[ip]);
    if (r<5) lay=1;
    else if (r>5 && r<10) lay=2;
    else if (r>10 && r<18) lay=3;
    else if (r>18 && r<30) lay=4;
    else if (r>30 && r<40) lay=5;
    else if (r>40) lay=6;
    if (lay<0) continue;
    nlay[lay-1]++;
  }
  Bool_t isgood=1;
  for (int jj=0; jj<6; jj++) {
    if (nlay[jj]<nreqpts_lay[jj]) isgood=0;
  }
  return isgood;
}
//////////////////////////////////////

int ITSAlignMille(int fromev=1, int toev=200, int maxentries=-1, char *outfile="ITSAlignMille.root") {

  // Read tracks from AliTrackPoints.N.root and fill Millepede.
  // Millepede results are written as AliAlignObjParams in outfile.
  
  int nreqpts=6;
  int nreqpts_lay[6]={0,0,0,0,0,0};

  TFile *fout=new TFile(outfile,"recreate");
  if (!fout->IsOpen()) {
    printf("\nCannot open output file!\n");
    return -1;
  }

  TChain *chain=new TChain("spTree");  
  char dir[100]="AliTrackPoints";
  char st[200];
  
  for (int iev=fromev; iev<=toev; iev++) {
    sprintf(st,"%s/AliTrackPoints.%d.root",dir,iev);
    chain->Add(st);
  }
  
  int nentries=chain->GetEntries();
  printf("There are %d entries in chain\n",nentries);
  
  if (maxentries>0 && maxentries<nentries) nentries=maxentries;
  
  AliTrackPointArray *tpa = 0;
  chain->SetBranchAddress("SP", &tpa);

  ////////////////////////////////////////////
  
  AliITSAlignMille *alig = new AliITSAlignMille("AliITSAlignMille.conf");

  Int_t nmod=alig->GetNModules();

  Double_t *parameters = new Double_t[nmod*6];
  Double_t *errors = new Double_t[nmod*6];
  Double_t *pulls = new Double_t[nmod*6];

  for(Int_t k=0;k<nmod*6;k++) {
    parameters[k]=0.;
    errors[k]=0.;
    pulls[k]=0.;
  }

  Double_t trackParams[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
  alig->InitGlobalParameters(parameters);
  alig->Print();

  ////////////////////////////////////////////////////////////////////

  TStopwatch stimer;
  stimer.Start();

  int iTrackOk=0; // number of good passed track
  // loop on spTree entries
  // one entry = one track
  for (int i=0; i<nentries; i++) {
    chain->GetEvent(i);
    //////////////////////////////////////////////////////
    // track preselection    
    int npts=tpa->GetNPoints();
    if (npts<nreqpts) continue;
    if (!SelectLayers(tpa,nreqpts_lay)) continue;
    //////////////////////////////////////////////////////

    if (!alig->ProcessTrack(tpa)) {
      if (!(iTrackOk%500)) 
	printf("Calling LocalFit n. %d\n",iTrackOk);
      alig->LocalFit(iTrackOk++,trackParams,0);
    }
  }
  printf("\nMillepede fed with %d tracks\n\n",iTrackOk);
  
  alig->GlobalFit(parameters,errors,pulls);
  
  cout << "Done with GlobalFit " << endl;
  
  ////////////////////////////////////////////////////////////
  // output
  

  TClonesArray *array = new TClonesArray("AliAlignObjParams",4000);
  TClonesArray &alobj = *array;

  // first object: ITS
  new(alobj[0]) AliAlignObjParams("ITS", 0, 0., 0., 0., 0., 0., 0., kTRUE);
  
  UShort_t volid;
  const Char_t *symname;
  Double_t dx,dy,dz,dpsi,dtheta,dphi;
  Double_t corr[21];

  // all ITS modules
  for (int idx=0; idx<2198; idx++) {
    volid=alig->GetModuleVolumeID(idx);
    symname = AliGeomManager::SymName(volid);
    for (int jj=0;jj<21;jj++) corr[jj]=0.0;

    if (alig->CheckVolumeID(volid)) { // defined modules
      alig->SetCurrentModule(idx);
      int iidx=alig->GetCurrentModuleInternalIndex();
      dx    = parameters[iidx*6+0];  corr[0] = errors[iidx*6+0]*errors[iidx*6+0];
      dy    = parameters[iidx*6+1];  corr[2] = errors[iidx*6+1]*errors[iidx*6+1];
      dz    = parameters[iidx*6+2];  corr[5] = errors[iidx*6+2]*errors[iidx*6+2];
      dpsi  = parameters[iidx*6+3];  corr[9] = errors[iidx*6+3]*errors[iidx*6+3];
      dtheta= parameters[iidx*6+4];  corr[14]= errors[iidx*6+4]*errors[iidx*6+4];
      dphi  = parameters[iidx*6+5];  corr[20]= errors[iidx*6+5]*errors[iidx*6+5];
    }
    else { // other modules
      dx=0.0; dy=0.0; dz=0.0; dpsi=0.0; dtheta=0.0; dphi=0.0;
    }
    
    new(alobj[idx+1]) AliAlignObjParams(symname, volid, dx, dy, dz, dpsi, dtheta, dphi, kFALSE);
    AliAlignObjParams* alo = (AliAlignObjParams*) array->UncheckedAt(idx+1);
    alo->SetCorrMatrix(corr);
  }
  
  fout->WriteObject(array,"ITSAlignObjs","kSingleKey");
  fout->Close();

  stimer.Stop();
  stimer.Print();
   
  return 0;
}
