//
// 8-feb-2002 new version with class AliEMCALJetMicroDst
//
class TH1F;
class TVector3;

void anaDst(Int_t mode=1);
Bool_t defineMicroDst(Int_t mode);

void bookHist1();
void bookHistPartonsAfterHardSc();
Bool_t fillInfoForPartons(TVector3 *vec1, TVector3 *vec2);
void compareJetPartons(TVector3 *vecJet, TVector3 *vec1, TVector3 *vec2);

void drawPartons();
void fitGauss(TH1F* hid, Int_t opt);

void cleanUpEvent();
