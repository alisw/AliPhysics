//
// 6-feb-2002 by PAI
//
class TH1F;
class TH2F;

void jetDst(Int_t mode=11, Int_t var=5, Int_t nmax=0);
TString *defineInput(Int_t mode);
void   defineJetFinder(Int_t mode, Int_t var);

void   testDST(char* fname="RES/FILE/jet100_1000ev_4.root"); 
// by PAI - 4-feb-2002
TH2F* bookTrigHist(TH2F *hid);
void  fillTriggerPatch(TH2F *hid, TH2F *htrig, TH1F* hE);

