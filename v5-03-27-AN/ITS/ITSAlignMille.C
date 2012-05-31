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


int ITSAlignMille(int fromev=1, int toev=1, int fromentry=0, int nentries=-1, char *outfile="ITSAlignMille.root", char *confile="AliITSAlignMille.conf",int nreqpts=3) {

  //AliLog::SetGlobalLogLevel(6);
  
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
  
  int toentry=chain->GetEntries();
  printf("There are %d entries in chain\n",toentry);
  
  if (nentries>0 && nentries<(toentry-fromentry)) toentry=fromentry+nentries;
  
  AliTrackPointArray *tpa = 0;
  chain->SetBranchAddress("SP", &tpa);

  ////////////////////////////////////////////
  
  AliITSAlignMille *alig = new AliITSAlignMille(confile);
  if (!alig->IsConfigured()) return -3;

  Int_t nmod=alig->GetNModules();
  alig->SetMinNPtsPerTrack(nreqpts);

  // require 4 points in SPD (one per layer, up and down)
  if (nreqpts>3) {
    alig->SetRequiredPoint("LAYER",1,1,1);
    alig->SetRequiredPoint("LAYER",1,-1,1);
    alig->SetRequiredPoint("LAYER",2,1,1);
    alig->SetRequiredPoint("LAYER",2,-1,1);
  }

  // correction for SSD bug (1) : inverisione sensor 18 e 19
  // alig->SetBug(1);

  Double_t *parameters = new Double_t[nmod*6];
  Double_t *errors = new Double_t[nmod*6];
  Double_t *pulls = new Double_t[nmod*6];

  for(Int_t k=0;k<nmod*6;k++) {
    parameters[k]=0.;
    errors[k]=0.;
    pulls[k]=0.;
  }

  // need array with size fNLocal*2
  Double_t trackParams[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  alig->InitGlobalParameters(parameters);
  alig->Print();

  Double_t sigmatra=alig->GetParSigTranslations();
  Double_t sigmarot=alig->GetParSigRotations();
  ////////////////////////////////////////////////////////////////////

  TStopwatch stimer;
  stimer.Start();

  /////////////////////

  int iTrackOk=0; // number of good passed track
  // loop on spTree entries
  // one entry = one track
  for (int i=fromentry; i<toentry; i++) {

    chain->GetEvent(i);
    if (!alig->ProcessTrack(tpa)) {
      if (!(iTrackOk%500)) 
	printf("Calling LocalFit n. %d\n",iTrackOk);
      alig->LocalFit(iTrackOk++,trackParams,0);
    }
  }
  printf("\nMillepede fed with %d tracks\n",iTrackOk);
  printf("Total number of rejected points because of bad EqLoc : %d\n\n",alig->GetTotBadLocEqPoints());
  
  stimer.Print();
  stimer.Continue(); 

  // make global fit
  alig->GlobalFit(parameters,errors,pulls);
  
  cout << "Done with GlobalFit " << endl;
  
  ////////////////////////////////////////////////////////////
  // output

  printf("\nProcessed points statistics:\n");
  Int_t maxstat=0;

  printf("\nOutput values and point statistics:\n\n");
  for (int i=0; i<nmod; i++) {
    Int_t statm=alig->GetProcessedPoints()[i];
    if (statm>maxstat) maxstat=statm;
    printf("index: %-4d   stat: %-7d   pars: %f   %f   %f   %f   %f   %f\n",alig->GetModuleIndexArray()[i], statm,
	   parameters[i*6+0],
	   parameters[i*6+1],
	   parameters[i*6+2],
	   parameters[i*6+3],
	   parameters[i*6+4],
	   parameters[i*6+5]);
    
  }
  FILE *pf=fopen("ITSAlignMille.out","w");
  if (pf) {
    fprintf(pf,"# param     dx         dy         dz         dpsi       dtheta     dphi\n");
    for (int i=0; i<nmod; i++) {
      fprintf(pf,"   %-5d %-10f %-10f %-10f %-10f %-10f %-10f\n",i,
	      parameters[i*6+0],
	      parameters[i*6+1],
	      parameters[i*6+2],
	      parameters[i*6+3],
	      parameters[i*6+4],
	      parameters[i*6+5]);
      
    }
    fclose(pf);
  }

  printf("Max statistics = %d\n",maxstat);
  if (maxstat<1) {
    printf("No points for alignment! quitting now...\n");
    return -1;
  }  

  TClonesArray *array = new TClonesArray("AliAlignObjParams",4000);
  TClonesArray &alobj = *array;

  // first object: ITS
  new(alobj[0]) AliAlignObjParams("ITS", 0, 0., 0., 0., 0., 0., 0., kTRUE);
  
  UShort_t volid;
  Char_t *symname;
  Double_t dx,dy,dz,dpsi,dtheta,dphi;
  Double_t corr[21];
  Double_t deltalocal[6];
  Double_t t[3],r[3];

  AliAlignObjParams *tmpa=NULL;

  // quality factor: 0 = NOT ALIGNED or OFF
  //                 N = ALIGNED with statistics = N  
  UInt_t QF; 

  // all ITS sensitive modules
  for (int idx=0; idx<2198; idx++) {
    volid=AliITSAlignMilleModule::GetVolumeIDFromIndex(idx);
    symname = AliGeomManager::SymName(volid);
    // default null misalignment
    for (int jj=0;jj<21;jj++) corr[jj]=0.0;
    for (int jj=0;jj<3;jj++) {t[jj]=0.0;r[jj]=0.0;}
    QF=0;

    if (alig->IsContained(volid)>=0) { // defined modules or inside a supermodule
      alig->SetCurrentModule(idx); // set the supermodule that contains idx
      int iidx=alig->GetCurrentModuleInternalIndex();
      
      // check if all 6 parameters have been evaluated by millepede
      Bool_t isoff=0;
      Double_t parsig=0;
      for (int kk=0; kk<6; kk++) {
	parsig=sigmatra;
	if (kk>2) parsig=sigmarot;
	if (pulls[iidx*6+kk]==0.0 && errors[iidx*6+kk]==parsig) isoff=1;
      }

      // check if module was fixed
      Bool_t isfixed=1;
      for (int kk=0; kk<6; kk++) {
	if (parameters[iidx*6+kk]!=0.0 || pulls[iidx*6+kk]!=0.0 || errors[iidx*6+kk]!=0.0) isfixed=0;
      }

      if (!isoff && !isfixed) { // OK, has been evaluated
	deltalocal[0] = parameters[iidx*6+0];  
	deltalocal[1] = parameters[iidx*6+1]; 
	deltalocal[2] = parameters[iidx*6+2]; 
	deltalocal[3] = parameters[iidx*6+3]; 
	deltalocal[4] = parameters[iidx*6+4];
	deltalocal[5] = parameters[iidx*6+5]; 
	tmpa = alig->GetCurrentModule()->GetSensitiveVolumeMisalignment(volid,deltalocal);
	if (!tmpa) {
	  printf("error transforming millepede parameters for module %d\n",idx);
	  continue;
	}
	tmpa->GetPars(t,r);
	// at the moment sigma of supermodule is given to sensitive modules. to be fixed...
	corr[0] = errors[iidx*6+0]*errors[iidx*6+0];
	corr[2] = errors[iidx*6+1]*errors[iidx*6+1];
	corr[5] = errors[iidx*6+2]*errors[iidx*6+2];
	corr[9] = errors[iidx*6+3]*errors[iidx*6+3];
	corr[14]= errors[iidx*6+4]*errors[iidx*6+4];
	corr[20]= errors[iidx*6+5]*errors[iidx*6+5];
	Int_t statm=alig->GetProcessedPoints()[iidx];
	QF=statm;
      }
      else {
	if (isoff) {
	  printf("Module %d is OFF!\n",idx);
	  QF=0;
	}
	if (isfixed) {
	  printf("Module %d was FIXED!\n",idx);
	  if (alig->GetPreAlignmentQualityFactor(idx)>0) 
	    QF=alig->GetPreAlignmentQualityFactor(idx);
	  else 
	    QF=0;
	}
      }

    }

    new(alobj[idx+1]) AliAlignObjParams(symname,volid,t[0],t[1],t[2],r[0],r[1],r[2],kFALSE);
    AliAlignObjParams* alo = (AliAlignObjParams*) array->UncheckedAt(idx+1);
    alo->SetCorrMatrix(corr);
    alo->SetUniqueID(QF);
  }
  
  fout->WriteObject(array,"ITSAlignObjs","kSingleKey");
  fout->Close();

  stimer.Stop();
  stimer.Print();

  return 0;
}
