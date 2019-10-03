class TPaveStats;

TH1* GetBaseHisto(TPad* pad=0);
TFrame* GetFrame(TPad* pad=0);
TPaveStats* SetStatPad(TH1* hst,float x1,float x2,float y1,float y2,Int_t stl=-1,Int_t col=-1);
TPaveStats* GetStatPad(TH1* hst);
void SetHStyle(TH1* hst,int col=kRed,int mark=20,float mrsize=0.7);
void SetGStyle(TGraph* hst,int col=kRed,int mark=20,float mrsize=0.7);
TH1* Cumulate(TH1* histo, Bool_t doErr=kTRUE, const char* copyName=0);
TLatex* AddLabel(const char*txt,float x=0.1,float y=0.9,int color=kBlack,float size=0.04);
void wAv(double v1,double v2, double err1=0,double err2=0, double* wv=0,double *we=0);
void wSum(double v1,double v2, double err1=0,double err2=0, double* wv=0,double *we=0);

void SaveCanvas(TCanvas* canv,const char* path="canv",const Option_t *option="ecg");
