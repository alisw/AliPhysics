#ifndef ALIHISTOGRAMSUE_H
#define ALIHISTOGRAMSUE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class  to handle histograms for UE analysis
//---------------------------------------------
////////////////////////////////////////////////

class TH1F;
class TH2F;
class TH1I;
class TObjArray;
class TProfile;
class TTree;
class TVector3;

class  AliHistogramsUE : public TObject
  {
  public:
    AliHistogramsUE();
    AliHistogramsUE(TList * list);
    virtual           ~AliHistogramsUE() { }
    AliHistogramsUE(const  AliHistogramsUE &det);
    AliHistogramsUE&   operator=(const  AliHistogramsUE &det);
    
    TObjArray*     CreateCanvas(const Int_t ncanv);
    TObjArray*     GetHistosForPlotting(TString file, TString branches);
    TList*	   CreateHistos(Int_t bins, Double_t min, Double_t max, Double_t etacut);
    void           DrawUE(Int_t debug);  //to draw final plots (normalized)
    void           FillHistogram(const char* name,Double_t fillX); //One dimensional
    void           FillHistogram(const char* name,Int_t fillX); //One dimensional
    void           FillHistogram(const char* name,Double_t fillX, Double_t fillY); //Two dimensional
    void           FillHistogram(const char* name,Double_t fillX, Double_t fillY, Double_t weight); //Two dimensional
    void           FillHistogram(const char* name,Double_t fillX, Int_t fillY, Double_t weight); //Two dimensional
    TList*	   GetListOfHistos();
    TH1F*            GetTrials()       {return fh1Trials;}
    TProfile*        GetXsec()         {return fh1Xsec;}
    void           PlotBranchesUE(TString file, TString branches, Double_t minJetProjection);  //TO BE CALLED BY EXTERNAL MACRO !!!		 
    void           SetStyle(); 
  protected:

  private:
    
    Int_t          fBinsPtInHist;    // Number of pT bins in histograms
    Double_t       fMinJetPtInHist;  // Minimum jet pT in histograms
    Double_t       fMaxJetPtInHist;  // Maximum jet pT in histograms
    Double_t       fTrackEtaCut;     // Track eta cut  
    TList*         fListOfHistos;    //  Output list of histograms
    
   
    // Histograms
    TH1F*  fhNJets;                  //!
    TH1F*  fhEleadingPt;             //!
    
    TH1F*  fhMinRegPtDist;           //!
    TH1F*  fhRegionMultMin;          //!
    TH1F*  fhMinRegAvePt;            //!
    TH1F*  fhMinRegSumPt;            //!
    TH1F*  fhMinRegMaxPtPart;        //!
    TH1F*  fhMinRegSumPtvsMult;      //!
    
    TH2F*  fhdNdEtaPhiDist;          //!
    TH2F*  fhFullRegPartPtDistVsEt;  //!
    TH2F*  fhTransRegPartPtDistVsEt; //!
    
    TH1F*  fhRegionSumPtMaxVsEt;     //!
    TH1I*  fhRegionMultMax;          //!
    TH1F*  fhRegionMultMaxVsEt;      //!
    TH1F*  fhRegionSumPtMinVsEt;     //!
    TH1F*  fhRegionMultMinVsEt;      //!
    TH1F*  fhRegionAveSumPtVsEt;     //!
    TH1F*  fhRegionDiffSumPtVsEt;    //!
    
    TH1F*  fhRegionAvePartPtMaxVsEt; //!
    TH1F*  fhRegionAvePartPtMinVsEt; //!
    TH1F*  fhRegionMaxPartPtMaxVsEt; //!
    
    TH2F*  fhRegForwardMult;         //!
    TH2F*  fhRegForwardSumPtvsMult;  //!
    TH2F*  fhRegBackwardMult;        //!
    TH2F*  fhRegBackwardSumPtvsMult; //!
    TH2F*  fhRegForwardPartPtDistVsEt; //!
    TH2F*  fhRegBackwardPartPtDistVsEt; //!
    TH2F*  fhRegTransMult;         //!
    TH2F*  fhRegTransSumPtVsMult;    //!
    TH2F*  fhMinRegSumPtJetPtBin;    //!
    TH2F*  fhMaxRegSumPtJetPtBin;    //!
    TH1F*  fhVertexMult;             //!
 
    TProfile*  fh1Xsec;		    //! 	
    TH1F*  fh1Trials;               //!

    ClassDef( AliHistogramsUE, 1 ); // Class to manage histograms in UE analysis
  };

#endif

    
