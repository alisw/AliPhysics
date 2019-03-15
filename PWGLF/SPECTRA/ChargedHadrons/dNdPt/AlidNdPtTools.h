#ifndef AlidNdPtTools_H
#define AlidNdPtTools_H

#include "THnSparse.h"

class AliESDtrackCuts;

class AlidNdPtTools : public TObject
{
    public:
        virtual                 ~AlidNdPtTools() = 0;

    public: 
        static Long64_t        FillHist(THnSparseD* s, Double_t x1, Double_t x2=0, Double_t x3=0, Double_t x4=0, Double_t x5=0, Double_t x6=0, Double_t x7 =0, Double_t x8 =0, Double_t x9 =0, Double_t x10 =0, Double_t x11 =0, Double_t x12 =0);
        static Int_t           AddAxis(const char* label, Int_t nbins, Double_t xmin, Double_t xmax, const char* option = 0);                    // options: <none>
        static Int_t           AddAxis(const char* label, const char* title, Int_t nbins, Double_t xmin, Double_t xmax, const char* option = 0); // options: <none>
        static Int_t           AddAxis(const char* label, Int_t nbins, Double_t* xbins, const char* option = 0);        // options: <none>
        static Int_t           AddAxis(const char* label, const char* title, Int_t nbins, Double_t* xbins, const char* option = 0);        // options: <none>
        static Int_t           AddAxis(const char* label, const char* title, const char* option);                                          // options: pt
        static Int_t           AddAxis(const char* label, const char* option);                                          // options: pt        
        static Int_t           AddAxis(const char* option);                                                             // options: pt
        static THnSparseD*     CreateHist(const char* name);
        static void            ResetHist() { if (fSparseTmp) delete fSparseTmp; }
        static TH1D*           CreateLogHist(const char* name, const char* title);                
        static TH1D*           CreateLogHist(const char* name);         
        static void            Log(TH1D* h, const char* name) { if (h) h->Fill(name,1); }
        
        static AliESDtrackCuts* CreateESDtrackCuts(const char* option); // options
        

    private:
        static THnSparseD*      fSparseTmp;     //! temporary histogram
 
     ClassDef(AlidNdPtTools, 1);
};

#endif
