/*!
   \file AliMCSpectraWeights.h
   \brief Class for re-weighting particle abundances in MC simulations
   \author Patrick Huhn
   \date 25/10/2019
*/
#ifndef __AliMCSpectraWeights__
#define __AliMCSpectraWeights__

class TParticle;
class AliMCEvent;

#include "TNamed.h"
#include "TObject.h"
#include "THn.h"
#include "TH1D.h"
#include "TString.h"
#include "TF1.h"
#include "TArrayD.h"

class AliMCSpectraWeights : public TNamed {
  public:
    enum ParticleType {kPion=0, kProtons=1, kKaon=2, kSigmaMinus=3, kSigmaPlus=4, kRest=5};/*!< enumerator of different particle types. */
    enum TaskState {kAllEmpty=0, kMCSpectraObtained, kDataFractionLoaded, kMCWeightCalculated};/*!< counter for the status of the task. */
    enum SysFlag {kNominal=0, kPionUp, kPionDown, kProtonUp, kProtonDown, kKaonUp, kKaonDown, kSigmaPlusUp, kSigmaPlusDown, kSigmaMinusUp, kSigmaMinusDown, kBylinkinUpper, kBylinkinLower, kHagedorn, kHagedornUpper, kHagedornLower, kExponential, kExponentialUpper, kExponentialLower, kBlastwave, kBlastwaveUpper, kBlastwaveLower};

    AliMCSpectraWeights();/*!< default root constructor */
    AliMCSpectraWeights(const char *collisionSystem, const char* name, AliMCSpectraWeights::SysFlag flag); /*!< constructor to be used.*/
    ~AliMCSpectraWeights();/*!< default destructor */

    void Init();/* Function to start initalizing after all setters are made. */
    Double_t GetMCSpectraWeight(TParticle* mcGenParticle, Float_t eventMultiplicityOrCentrality);/*main function to use. Will deliver correct weights to re-weight the abundances of different particle species */
    Double_t GetMCSpectraWeight(TParticle* mcGenParticle, AliMCEvent* mcEvent);/* preferable to use this */
    void FillMCSpectra(AliMCEvent* mcEvent);/*function to fill internal mc spectra for calculation of weight factors*/

    //Setter
    void SetBinsPt(TArrayD *bins){if(fBinsPt) delete fBinsPt; fBinsPt = new TArrayD(*bins);}/*!< Set bins in Pt using a TArrayD */
    void SetBinsPt(Int_t nBins, Double_t *binEdges){if(fBinsPt) delete fBinsPt; fBinsPt = new TArrayD(nBins+1,binEdges);}/*!< Set bins in Pt using number of bins and array of bin edges */
    void SetBinsMultCent(TArrayD *bins){if(fBinsMultCent) delete fBinsMultCent; fBinsMultCent = new TArrayD(*bins);}/*!< Set bins in Multiplicity/Centrality using a TArrayD */
    void SetBinsMultCent(Int_t nBins, Double_t *binEdges){if(fBinsMultCent) delete fBinsMultCent; fBinsMultCent = new TArrayD(nBins+1,binEdges);}/*!< Set bins in Multiplicity/Centrality number of bins and array of bin edges */

    void SetMCSpectraFile(const char* file){fstFileMCSpectra=file;}
    void SetDataFractionsFile(const char* file){fstFilePublished=file;}
    void SetCollisionSystem(const char* system) {fstCollisionSystem=system;}
    void SetUseMultiplicity(Bool_t bMult) {fUseMultiplicity = bMult;}
    void SetSavedObjName(const char* name){fstSavedObjName=name;}
    void SetSavedListName(const char* name){fstSavedListName=name;}
    void SetSysFlag(SysFlag flag){fFlag=flag;}

    // Getter
    TArrayD* GetBinsPt() {return fBinsPt;}
    TString* GetParticleTypes() const {return fstPartTypes;}
    Int_t GetNPartTypes() const {return fNPartTypes;}
    Int_t GetTaskStatus() const {return fbTaskStatus;}
    THnF* GetHistMCGenPrimTrackParticles() const {return fHistMCGenPrimTrackParticle;}
    THnF* GetHistDataFraction() const {return fHistDataFractions;}
    THnF* GetHistMCFraction() const {return fHistMCFractions;}
    THnF* GetHistMCWeights() const {return fHistMCWeights;}
    AliMCSpectraWeights::SysFlag GetSysFlag() const {return fFlag;}
    Double_t GetMultOrCent() const {return fMultOrCent;}

    //usefull functions
    Int_t IdentifyMCParticle(TParticle* mcParticle);

    ClassDef(AliMCSpectraWeights,1); // Class for reweighting mc particle abundances

  private:
    void InitHistos();//!
    void LoadMeasuredFractions();//!
    Bool_t LoadFromAliMCSpectraWeight(AliMCSpectraWeights* obj = 0);//!
    Bool_t LoadFromTHnF(const char* histname);//!
    Bool_t CalculateMCWeights();//!
    Bool_t CalcMCFractions();//!
    Bool_t CorrectFractionsforRest();//!
    Double_t GetMultFromCent(int CentBin);//!
    Double_t GetCentFromMult(double dMult);//!
    void CountEventMult();//!
    // histograms
    THnF*       fHistMCGenPrimTrackParticle;// Histogram for MC particle information
    THnF*       fHistDataFractions;// Histogram for particle abundances from published data
    THnF*       fHistMCFractions;// Histogram for particle abundances from MC
    THnF*       fHistMCWeights;// Histogram for weight factors to re-weight MC abundances to data ones.

    TArrayD*    fBinsPt;/* pT binning */
    TArrayD*    fBinsMultCent;/* centrality or multiplicity binning */
    TString     fstCollisionSystem;/* collision system */
    Int_t       fNPartTypes;/* number of particle species */
    TString*    fstPartTypes;/* Array of used particle species */
    int         fNCentralities;
    TString*    fstCentralities;

    // paths to files
    TString     fstFileMCSpectra;/* path to previous train output of MC fractions */
    TString     fstFilePublished;/* path to calculated fractions from published spectra */
    TString     fstSavedObjName;
    TString     fstSavedListName;

    // configs
    Bool_t        fUseMultiplicity; /* switch to use multiplicity instead of centrality */
    Int_t         fbTaskStatus;/* controls internal status of class */
    AliMCSpectraWeights::SysFlag       fFlag;// enum for systematic variation;

    AliMCEvent    *fMCEvent;		//!<! MC event
    Double_t      fMultOrCent;//!
};

#endif /* __AliMCSpectraWeights__ */
