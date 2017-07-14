#ifndef ALIV0READERV1_H
#define ALIV0READERV1_H

#include "AliAnalysisTaskSE.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliConversionPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "AliExternalTrackParam.h"
#include "TObject.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliKFParticle.h"
#include "TParticle.h"
#include <iterator>
#include <vector>
#include "AliESDpid.h"
#include "TF1.h"
#include "TRandom3.h"
#include "AliAnalysisManager.h"

class AliConversionPhotonBase;
class TRandom3;
class TList;
class AliKFConversionPhoton;
class TString;
class TClonesArray;
class TH1F;
class TH2F;
class AliAODConversionPhoton;

using namespace std;

class AliV0ReaderV1 : public AliAnalysisTaskSE {

  public:

    class iterator : public std::iterator<std::bidirectional_iterator_tag, AliConversionPhotonBase> {
    public:
      enum Direction_t {
        kForwardDirection = 0,
        kBackwardDirection = 1
      };
      iterator(const AliV0ReaderV1 *reader, Direction_t dir, int position);
      iterator(const iterator &other);
      iterator &operator=(const iterator &other);
      virtual ~iterator() {}

      bool operator!=(iterator &other) const;
      iterator operator++(int);
      iterator &operator++();
      iterator operator--(int);
      iterator &operator--();
      AliConversionPhotonBase *operator*();

    private:
      const AliV0ReaderV1     *fkData;          ///< V0 reader used to iterate over
      int                      fCurrentIndex;   ///< Index of the current element
      Direction_t              fDirection;      ///< Iterator in forward direction
    };

    AliV0ReaderV1(const char *name="V0ReaderV1");
    virtual                    ~AliV0ReaderV1();                            //virtual destructor

    void                      UserCreateOutputObjects();
    virtual Bool_t            Notify();
    virtual void              UserExec(Option_t *option);
    virtual void              Terminate(Option_t *);
    virtual void              Init();

    Bool_t                    ProcessEvent (AliVEvent *inputEvent, AliMCEvent *mcEvent=NULL);
    Bool_t                    IsEventSelected()                     {return fEventIsSelected;}

    // Return Reconstructed Gammas
    TClonesArray*             GetReconstructedGammas() const        {return fConversionGammas;}
    Int_t                     GetNReconstructedGammas() const       {if(fConversionGammas){return fConversionGammas->GetEntriesFast();} else{ return 0;}}
    AliConversionPhotonBase *operator[](int index) const;

    AliConversionPhotonCuts*  GetConversionCuts()                   {return fConversionCuts;}
    AliConvEventCuts*         GetEventCuts()                        {return fEventCuts;}
    TList*                    GetCutHistograms()                    {if(fConversionCuts) {return fConversionCuts->GetCutHistograms();}
                                                                     return NULL;}
    TList*                    GetEventCutHistograms()               {if(fEventCuts) {return fEventCuts->GetCutHistograms();}
                                                                     return NULL;}
    TString                   GetCurrentFileName()                  {return fCurrentFileName;}
    // Set Options
    void	       SetAddv0sInESDFilter(Bool_t addv0s)	{kAddv0sInESDFilter = addv0s;}
    void               CountTracks();
    void               CountTPCoutTracks();
    void               SetConversionCuts(const TString cut);
    void               SetConversionCuts(AliConversionPhotonCuts *cuts) {fConversionCuts=cuts; return;}
    void               SetEventCuts(const TString cut);
    void               SetEventCuts(AliConvEventCuts *cuts)             {fEventCuts=cuts; return;}

    void               SetUseOwnXYZCalculation(Bool_t flag)             {fUseOwnXYZCalculation=flag; return;}
    void               SetUseConstructGamma(Bool_t flag)                {fUseConstructGamma=flag; return;}
    void               SetUseAODConversionPhoton(Bool_t b)              {if(b){ cout<<"Setting Outputformat to AliAODConversionPhoton "<<endl;}
                                                                         else { cout<<"Setting Outputformat to AliKFConversionPhoton "<<endl;}
                                                                         kUseAODConversionPhoton=b; return;}

    void               SetCreateAODs(Bool_t k=kTRUE)                    {fCreateAOD=k; return;}
    void               SetDeltaAODFilename(TString s)                   {fDeltaAODFilename=s; return;}
    void               SetDeltaAODBranchName(TString string)            {fDeltaAODBranchName = string;
                                                                         fRelabelAODs = kTRUE;
                                                                         AliInfo(Form("Set DeltaAOD BranchName to: %s",fDeltaAODBranchName.Data()));
                                                                         AliInfo(Form("Relabeling of AODs has automatically been switched ON!"));
                                                                         return;}

    void               RelabelAODs(Bool_t relabel=kTRUE)                {fRelabelAODs=relabel; return;}
    Bool_t             AreAODsRelabeled()                               {return fRelabelAODs;}
    Int_t              IsReaderPerformingRelabeling()                   {return fPreviousV0ReaderPerformsAODRelabeling;}
    void               RelabelAODPhotonCandidates(AliAODConversionPhoton *PhotonCandidate);
    void               SetPeriodName(TString name)                      {fPeriodName = name;
                                                                         AliInfo(Form("Set PeriodName to: %s",fPeriodName.Data()));
                                                                         return;}
    TString            GetPeriodName()                                  {return fPeriodName;}
    Int_t              GetPtHardFromFile()                              {return fPtHardBin;}
    Int_t              GetNumberOfPrimaryTracks()                       {return fNumberOfPrimaryTracks;}
    Int_t              GetNumberOfTPCoutTracks()                        {return fNumberOfTPCoutTracks;}
    void               SetUseMassToZero (Bool_t b)                      {if(b){ cout<<"enable set mass to zero for AliAODConversionPhoton"<<endl;}
                                                                         else { cout<<"disable set mass to zero for AliAODConversionPhoton "<<endl;}
                                                                         fUseMassToZero=b; return;}

    void               SetProduceV0FindingEfficiency(Bool_t b)          {fProduceV0findingEffi = b;
                                                                         if(b) AliInfo("Enabled V0finding Efficiency");
                                                                         return;}

    Bool_t             GetProduceV0FindingEfficiency()                  {return fProduceV0findingEffi;}
    TList*             GetV0FindingEfficiencyHistograms()               {return fHistograms;}
    void               SetProduceImpactParamHistograms(Bool_t b)        {fProduceImpactParamHistograms = b;
                                                                         if(b) AliInfo("Producing additional impact parameter histograms");
                                                                         return;}

    Bool_t             GetProduceImpactParamHistograms()                {return fProduceImpactParamHistograms;}
    TList*             GetImpactParamHistograms()                       {return fImpactParamHistograms;}

    Bool_t             ParticleIsConvertedPhoton(AliMCEvent *mcEvent, TParticle *particle, Double_t etaMax, Double_t rMax, Double_t zMax);
    void               CreatePureMCHistosForV0FinderEffiESD();
    void               FillRecMCHistosForV0FinderEffiESD( AliESDv0* currentV0);
    void               FillImpactParamHistograms(AliVTrack *ptrack, AliVTrack* ntrack, AliESDv0 *fCurrentV0, AliKFConversionPhoton *fCurrentMotherKF);
    Bool_t             CheckVectorOnly(vector<Int_t> &vec, Int_t tobechecked);
    Bool_t             CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked);
    void               SetImprovedPsiPair(Int_t p)                      {fImprovedPsiPair=p;return;}
    Int_t              GetImprovedPsiPair()                             {return fImprovedPsiPair;}
  
    iterator           begin() const                                    {return iterator(this, iterator::kForwardDirection, 0);}
    iterator           end() const                                      {return iterator(this, iterator::kForwardDirection, GetNReconstructedGammas());}
    iterator           rbegin() const                                   {return iterator(this, iterator::kBackwardDirection, GetNReconstructedGammas() -1); }
    iterator           rend() const                                     {return iterator(this, iterator::kBackwardDirection, -1);}

  protected:
    // Reconstruct Gammas
    Bool_t                  ProcessESDV0s();
    AliKFConversionPhoton*  ReconstructV0(AliESDv0* fCurrentV0,Int_t currentV0Index);
    void                    FillAODOutput();
    void                    FindDeltaAODBranchName();
    Bool_t                  GetAODConversionGammas();

    // Getter Functions
    const AliExternalTrackParam*   GetExternalTrackParam(AliESDv0 *fCurrentV0, Int_t &tracklabel, Int_t charge);
    const AliExternalTrackParam*   GetExternalTrackParamP(AliESDv0 *fCurrentV0, Int_t &tracklabel) {return GetExternalTrackParam(fCurrentV0,tracklabel,1);}
    const AliExternalTrackParam*   GetExternalTrackParamN(AliESDv0 *fCurrentV0, Int_t &tracklabel) {return GetExternalTrackParam(fCurrentV0,tracklabel,-1);}
    AliKFParticle*                 GetPositiveKFParticle(AliAODv0 *fCurrentV0, Int_t fTrackLabel[2]);
    AliKFParticle*                 GetNegativeKFParticle(AliAODv0 *fCurrentV0, Int_t fTrackLabel[2]);
    AliKFParticle*                 GetPositiveKFParticle(AliESDv0 *fCurrentV0, Int_t fTrackLabel[2]);
    AliKFParticle*                 GetNegativeKFParticle(AliESDv0 *fCurrentV0, Int_t fTrackLabel[2]);

    Bool_t               GetConversionPoint(const AliExternalTrackParam *pparam, const AliExternalTrackParam *nparam, Double_t convpos[3], Double_t dca[2]);
    Bool_t               GetHelixCenter(const AliExternalTrackParam *track, Double_t center[2]);
    Double_t             GetPsiPair(const AliESDv0* v0, const AliExternalTrackParam *positiveparam, const AliExternalTrackParam *negativeparam, const Double_t convpos[3]) const;
    Bool_t 	   kAddv0sInESDFilter; 	          // Add PCM v0s to AOD created in ESD filter
    TBits		     *fPCMv0BitField;	  // Pointer to bitfield of PCM v0s
    AliConversionPhotonCuts  *fConversionCuts;    // Pointer to the ConversionCut Selection
    AliConvEventCuts         *fEventCuts;         // Pointer to the ConversionCut Selection
    TClonesArray             *fConversionGammas;  // TClonesArray holding the reconstructed photons
    Bool_t         fUseImprovedVertex;            // set flag to improve primary vertex estimation by adding photons
    Bool_t         fUseOwnXYZCalculation;         //flag that determines if we use our own calculation of xyz (markus)
    Bool_t         fUseConstructGamma;            //flag that determines if we use ConstructGamma method from AliKF
    Bool_t         kUseAODConversionPhoton;       // set flag to use AOD instead of KF output format for photons
    Bool_t         fCreateAOD;                    // set flag for AOD creation
    TString        fDeltaAODBranchName;           // File where Gamma Conv AOD is located, if not in default AOD
    TString        fDeltaAODFilename;             // set filename for delta/satellite aod
    Bool_t         fRelabelAODs;                  //
    Int_t          fPreviousV0ReaderPerformsAODRelabeling; // 0->not set, meaning V0Reader has not yet determined if it should do AODRelabeling, 1-> V0Reader perfomrs relabeling, 2-> previous V0Reader in list perfomrs relabeling
    Bool_t         fEventIsSelected;
    Int_t          fNumberOfPrimaryTracks;        // Number of Primary Tracks in AOD or ESD
    Int_t          fNumberOfTPCoutTracks;        // Number of TPC Tracks with TPCout flag
    TString        fPeriodName;
    Int_t          fPtHardBin;                    // ptHard bin from file
    Bool_t         fUseMassToZero;                // switch on setting the mass to 0 for AODConversionPhotons
    Bool_t         fProduceV0findingEffi;         // enable histograms for V0finding efficiency
    Bool_t         fProduceImpactParamHistograms; // enable histograms of impact parameters
    Float_t        fCurrentInvMassPair;           // Invariant mass of the pair
    Int_t          fImprovedPsiPair;              // enables the calculation of PsiPair after the precise calculation of R and use of the proper function for propagation
    TList         *fHistograms;                   // list of histograms for V0 finding efficiency
    TList         *fImpactParamHistograms;        // list of histograms of impact parameters
    TH2F          *fHistoMCGammaPtvsR;            // histogram with all converted gammas vs Pt and R (eta < 0.9)
    TH2F          *fHistoMCGammaPtvsPhi;          // histogram with all converted gammas vs Pt and Phi (eta < 0.9)
    TH2F          *fHistoMCGammaPtvsEta;          // histogram with all converted gammas vs Pt and Eta
    TH2F          *fHistoMCGammaRvsPhi;           // histogram with all converted gammas vs R and Phi (eta < 0.9)
    TH2F          *fHistoMCGammaRvsEta;           // histogram with all converted gammas vs R and Eta
    TH2F          *fHistoMCGammaPhivsEta;         // histogram with all converted gammas vs Phi and Eta
    TH2F          *fHistoRecMCGammaPtvsR;         // histogram with all reconstructed converted gammas vs Pt and R (eta < 0.9)
    TH2F          *fHistoRecMCGammaPtvsPhi;       // histogram with all reconstructed converted gammas vs Pt and Phi (eta < 0.9)
    TH2F          *fHistoRecMCGammaPtvsEta;       // histogram with all reconstructed converted gammas vs Pt and Eta
    TH2F          *fHistoRecMCGammaRvsPhi;        // histogram with all reconstructed converted gammas vs R and Phi (eta < 0.9)
    TH2F          *fHistoRecMCGammaRvsEta;        // histogram with all reconstructed converted gammas vs R and Eta
    TH2F          *fHistoRecMCGammaPhivsEta;      // histogram with all reconstructed converted gammas vs Phi and Eta
    TH1F          *fHistoRecMCGammaMultiPt;       // histogram with all at least double counted photons vs Pt (eta < 0.9)
    TH2F          *fHistoRecMCGammaMultiPtvsEta;  // histogram with all at least double counted photons vs Pt vs Eta
    TH1F          *fHistoRecMCGammaMultiR;        // histogram with all at least double counted photons vs R (eta < 0.9)
    TH1F          *fHistoRecMCGammaMultiPhi;      // histogram with all at least double counted photons vs Phi (eta < 0.9)
    TH1F          *fHistoPosTrackImpactParamZ;    //impact parameter z of positive track of V0
    TH1F          *fHistoPosTrackImpactParamY;
    TH1F          *fHistoPosTrackImpactParamX;
    TH2F          *fHistoPosTrackImpactParamZvsPt;
    TH2F          *fHistoPosTrackImpactParamYvsPt;
    TH2F          *fHistoPosTrackImpactParamXvsPt;
    TH1F          *fHistoNegTrackImpactParamZ;
    TH1F          *fHistoNegTrackImpactParamY;
    TH1F          *fHistoNegTrackImpactParamX;
    TH2F          *fHistoNegTrackImpactParamZvsPt;
    TH2F          *fHistoNegTrackImpactParamYvsPt;
    TH2F          *fHistoNegTrackImpactParamXvsPt;
    TH2F          *fHistoImpactParamZvsR;         // conversion point z vs conversion radius
    TH2F          *fHistoImpactParamZvsR2;        // after cuts
    TH1F          *fHistoPt;
    TH1F          *fHistoPt2;                     // Pt after Impact parameter and causality cuts
    TH1F          *fHistoDCAzPhoton;
    TH1F          *fHistoDCAzPhoton2;             // photon dca after impact parameter and causality cuts
    TH1F          *fHistoR;                       // conversion radius
    TH1F          *fHistoRrecalc;                 // recalculated conversion radius
    TH1F          *fHistoRviaAlpha;                       // conversion radius
    TH1F          *fHistoRviaAlphaRecalc;                 // recalculated conversion radius
    TH1F          *fHistoRdiff;                   // difference in R between conflict cluster and conversion radius
    TH1F          *fHistoImpactParameterStudy;    // info about which cut rejected how many V0s
    TTree         *fImpactParamTree;               // tree with y, pt and conversion radius 
   
    vector<Int_t>  fVectorFoundGammas;            // vector with found MC labels of gammas
    TString       fCurrentFileName;               // current file name
    Bool_t        fMCFileChecked;                 // vector with MC file names which are broken
    
  private:
    AliV0ReaderV1(AliV0ReaderV1 &original);
    AliV0ReaderV1 &operator=(const AliV0ReaderV1 &ref);

    ClassDef(AliV0ReaderV1, 16)

};

inline void AliV0ReaderV1::SetConversionCuts(const TString cut){
  if(fConversionCuts != NULL){
  delete fConversionCuts;
    fConversionCuts=NULL;
  }
  if(fConversionCuts == NULL){
    fConversionCuts=new AliConversionPhotonCuts("V0ReaderCuts","V0ReaderCuts");
    fConversionCuts->InitializeCutsFromCutString(cut.Data());
  }
}

inline void AliV0ReaderV1::SetEventCuts(const TString cut){
  if(fEventCuts != NULL){
  delete fEventCuts;
  fEventCuts=NULL;
  }
  if(fEventCuts == NULL){
    fEventCuts=new AliConvEventCuts("V0ReaderEventCuts","V0ReaderEventCuts");
    fEventCuts->InitializeCutsFromCutString(cut.Data());
  }
}

#endif
