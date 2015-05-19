#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TList.h>
#include <TH1.h>
#include <TAxis.h>
#endif

// show DOFs with statistics below given threshold

void showLowStat(const char* stfile, int thr=30)
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
  int cnt=0;
  for (int i=1;i<=ndof;i++) {
    if (hstdof->GetBinContent(i)<thr) {
      printf("%4d\t%-50s\t%7d\n",cnt++,xax->GetBinLabel(i),(int)hstdof->GetBinContent(i));
    }
  }
  //
  lst->SetOwner();
  delete lst;
  fl->Close();
  delete fl;
}
