#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliAlgMPRecord.h"
#include "Mille.h"
#include <TFile.h>
#include <TChain.h>
#include <TString.h>
#include <TSystem.h>
#include <TClassTable.h>
#include <TMath.h>
//
#endif


// convert MPRecord to Mille format
const char* recBranchName = "mprec";
const char* recTreeName   = "mpTree";
const char* defOutName    = "mpData";
const int defSplit = -200;   // default output chunk size in MB


TChain* LoadMPrecChain(const char* inpData, const char* chName=recTreeName);
int     ConvertAndStore(AliAlgMPRecord* rec, Mille* mille);
Bool_t  ProcessMPRec(AliAlgMPRecord* rec);


void MPRec2Mille
(const char* inpName,                // name of MPRecord file or list of files
 const char* outname = defOutName,   // out file name
 int split           = defSplit      // 0: no split, >0: on N tracks,<0: on size in MB
 )
{
  //
  if (gClassTable->GetID("AliAlgSteer")<0) gSystem->Load("libAlg.so");
  //
  TChain* mprChain = LoadMPrecChain(inpName);
  if (!mprChain) return;
  int nEnt = mprChain->GetEntries();
  //
  TBranch* br = mprChain->GetBranch(recBranchName);
  if (!br) {
    printf("provided tree does not contain branch mprec\n");
    return;
  }
  //
  AliAlgMPRecord* mprec = new AliAlgMPRecord();
  br->SetAddress(&mprec);
  //
  TString mln = outname;
  if (mln.IsNull()) mln = inpName; // use inpname + ".mille"
  if (mln.EndsWith(".mille")>0) mln.Resize(mln.Last('.'));
  printf(">>%s \n<<%s%s%s\n",inpName, mln.Data(),split ? "_XXX":"",".mille");
  if (split) printf("Split on %d %s\n",TMath::Abs(split),split>0 ? "tracks":"MB");
  //
  TString milleName;
  Mille* mille=0;
  int cntTr=0,cntTot=0, cntMille=0;
  double sizeW=0,sizeWTot=0;
  if (split<0) split *= 1000000;
  for (int i=0;i<nEnt;i++) {
    mprChain->GetEntry(i);
    //
    if (!ProcessMPRec(mprec)) continue; // preprocess and skip if needed
    //
    if (!mille || 
	(split>0 && ++cntTr>split) ||
	(split<0 && sizeW>-split) 
	) { // start new mille file
      cntTr = 0;
      sizeW = 0;
      delete mille; 
      if (split) milleName = Form("%s_%03d.%s",mln.Data(),cntMille,"mille");
      else       milleName = Form("%s.%s",mln.Data(),"mille");
      cntMille++;
      printf("Opening output file %s\n",milleName.Data());
      mille = new Mille(milleName.Data());
    }
    cntTot++;
    int nbwr = ConvertAndStore(mprec,mille);
    sizeW += nbwr;
    sizeWTot += nbwr;
  }
  if (mille) delete mille;
  //
  br->SetAddress(0);
  delete mprec;
  delete mprChain;
  //
  printf("converted %d tracks out of %d\n(%ld MB in %d %s%s.mille files written)\n",
	 cntTot,nEnt,long(sizeWTot/1e6),cntMille,mln.Data(),split ? "_XXX":"");
  //
}

//_________________________________________________________
int ConvertAndStore(AliAlgMPRecord* rec, Mille* mille)
{
  // convert and store the record
  static TArrayF buffDLoc;
  int nr = rec->GetNResid(); // number of residual records
  int nloc = rec->GetNVarLoc();
  if (buffDLoc.GetSize()<nloc) buffDLoc.Set(nloc+100);
  float *buffLocV = buffDLoc.GetArray();
  const float *recDGlo = rec->GetArrGlo();
  const float *recDLoc = rec->GetArrLoc();
  const short *recLabLoc = rec->GetArrLabLoc();
  const int   *recLabGlo = rec->GetArrLabGlo();
  //
  for (int ir=0;ir<nr;ir++) {
    memset(buffLocV,0,nloc*sizeof(float));      
    int ndglo = rec->GetNDGlo(ir);
    int ndloc = rec->GetNDLoc(ir);
    // fill 0-suppressed array from MPRecord to non-0-suppressed array of Mille
    for (int l=ndloc;l--;) buffLocV[recLabLoc[l]] = recDLoc[l]; 
    //
    mille->mille(nloc,buffLocV,ndglo,recDGlo,recLabGlo,rec->GetResid(ir),rec->GetResErr(ir));
    //
    recLabGlo += ndglo; // next record
    recDGlo   += ndglo;
    recLabLoc += ndloc;
    recDLoc   += ndloc;
  }
  return mille->end(); // bytes written
  //
}


//____________________________________________________________________
TChain* LoadMPrecChain(const char* inpData, const char* chName)
{
  TChain* chain = new TChain(chName);
  //
  TString inpDtStr = inpData;
  if (inpDtStr.EndsWith(".root")) {
    chain->AddFile(inpData);
  }
  else {
    //
    ifstream inpf(inpData);
    if (!inpf.good()) {
      printf("Failed on input filename %s\n",inpData);
      return 0;
    }
    //
    TString flName;
    flName.ReadLine(inpf);
    while ( !flName.IsNull() ) {
      flName = flName.Strip(TString::kBoth,' ');
      if (flName.BeginsWith("//") || flName.BeginsWith("#")) {flName.ReadLine(inpf); continue;}
      flName = flName.Strip(TString::kBoth,',');
      flName = flName.Strip(TString::kBoth,'"');
      printf("Adding %s\n",flName.Data());
      chain->AddFile(flName.Data());
      flName.ReadLine(inpf);
    }
  }
  //
  int n = chain->GetEntries();
  if (n<1) {
    printf("Obtained chain is empty\n");
    return 0;
  }
  else printf("Opened %s chain with %d entries\n",chName,n);
  return chain;
}

//_________________________________________________________
Bool_t  ProcessMPRec(AliAlgMPRecord* rec)
{
  // put here user code
  //
  return kTRUE;
}

