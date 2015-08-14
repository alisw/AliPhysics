#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TList.h>
#include <TTree.h>
#include <TH1.h>
#include <TAxis.h>
#include <TString.h>
#include "AliAlgMPRecord.h"
#include "AliAlgSteer.h"
#include "Mille.h"
#endif

// Generate fake MP records for DOFs below statistics threshold
// If argument bin==kTRUE, then write directly pede binary file
// otherwise, write MPRecord

TString mpDummy = "mpDummy";
const Float_t kDummyDer = 1e-6;
const Float_t kDummyRes = 0;
const Float_t kDummyErr = 999.;

void MPDummyForLowStat(const char* stfile, int thr=30, int nGen=40,Bool_t bin=kTRUE)
{
  // show degrees of freedom with low stat
  TFile* fl = TFile::Open(stfile);
  if (!fl) {printf("Failed to open %s\n",stfile); return;}
  TList* lst = (TList*)fl->Get("clist");
  if (!lst) {printf("No clist in %s\n",stfile); return;}
  TH1* hstdof = (TH1*)lst->FindObject("DOFstat");
  if (!hstdof) {printf("No DOFstat histo in %s\n",stfile); return;}
  //
  int ndof = hstdof->GetNbinsX();
  TAxis* xax = hstdof->GetXaxis();
  printf("%4s\t%-50s\t%s","cnt"," DOF ID_name","entries");
  Mille* ml = 0;
  AliAlgSteer* algSteer=0;
  AliAlgMPRecord* mpRec=0;
  //
  if (bin) ml = new Mille(Form("%s.%s",mpDummy.Data(),"mille"));
  else {
    algSteer = new AliAlgSteer();
    algSteer->SetMPDatFileName(mpDummy.Data());
    algSteer->SetMPOutType(AliAlgSteer::kMPRec);
    algSteer->InitMPRecOutput();
    mpRec = algSteer->GetMPRecord();
  }
  //
  int   labDum[1] = {0}, cnt=0;
  float locDum[1] = {0}, gloDum[1] = {kDummyDer};
  //
  for (int i=1;i<=ndof;i++) {
    if (hstdof->GetBinContent(i)>thr) continue;
    TString labS = xax->GetBinLabel(i);
    printf("%4d\t%-50s\t%7d\n",cnt++,labS.Data(),(int)hstdof->GetBinContent(i));
    int indL = labS.Index("_");
    if (indL>0) labS.Resize(indL);
    else {
      printf("Failed to extract label from %s\n",labS.Data());
      exit(1);
    }
    labDum[0] = labS.Atoi();
    if (bin) {
      for (int j=nGen;j--;) {
	ml->mille(0, locDum, 1, gloDum, labDum, kDummyRes, kDummyErr);
	ml->end();
      }
    }
    else {
      mpRec->DummyRecord(kDummyRes,kDummyErr,kDummyDer,labDum[0]);
      for (int j=nGen;j--;) algSteer->GetMPRecTree()->Fill();
    }
  }
  //
  if (bin) delete ml;
  else {
    algSteer->CloseMPRecOutput();
    delete algSteer;
  }
  //
  lst->SetOwner();
  delete lst;
  fl->Close();
  delete fl;
}

