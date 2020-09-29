
/**
 * @file AliMCSpectraWeights.h
 * @brief Class for re-weighting particle abundances in MC simulations
 * @author Patrick Huhn
 * @date 25.10.2019
 */
#ifndef __AliMCSpectraWeights__
#define __AliMCSpectraWeights__

//#define __AliMCSpectraWeights_DEBUG__

class TParticle;
class AliMCEvent;
class TH3F;
class TH1D;
class TF1;
#include <string>
#include <vector>

/**
 * @class AliMCSpectraWeights AliMCSpectraWeights.h "AliMCSpectraWeights.h"
 * @brief TODO
 *
 *  TODO write something up
 *  TODO write full explanation here
 */
class AliMCSpectraWeights {
public:
    enum ParticleType {
        kPion = 0,
        kProtons = 1,
        kKaon = 2,
        kSigmaMinus = 3,
        kSigmaPlus = 4,
        kRest = 5
    }; /*!< enumerator of different particle types. */
    enum TaskState {
        kAllEmpty = 0,
        kMCSpectraObtained,
        kDataFractionLoaded,
        kMCWeightCalculated
    }; /*!< counter for the status of the task. */
    enum SysFlag {
        kNominal = 0,
        kPionUp,
        kPionDown,
        kProtonUp,
        kProtonDown,
        kKaonUp,
        kKaonDown,
        kSigmaPlusUp,
        kSigmaPlusDown,
        kSigmaMinusUp,
        kSigmaMinusDown,
        kBylinkinUpper,
        kBylinkinLower,
        kHagedorn,
        kHagedornUpper,
        kHagedornLower,
        kExponential,
        kExponentialUpper,
        kExponentialLower,
        kBlastwave,
        kBlastwaveUpper,
        kBlastwaveLower
    }; /*!< enumerator for systematic variations */
    
private:
    std::string fstCollisionSystem; /*!< collision system */
    std::string fstFileMCSpectra; /*!< path to previous train output of MC fractions */
    std::string fstFilePublished; /*!< path to calculated fractions from published spectra */
    std::string fstSavedObjName;
    std::string fstSavedListName;
    std::vector<std::string> fstPartTypes; /*!< Array of used particle species */
    std::vector<std::string> fstCentralities;
    std::vector<float>fBinsMultCent;        /*!< centrality or multiplicity binning */
    std::vector<float> fBinsPt; /*!< pT binning */
    TH3F* fHistMCGenPrimTrackParticle;  /*!< Histogram for MC particle information*/
    TH3F*  fHistDataFractions; /*!< Histogram for particle abundances from published data */
    TH3F*
    fHistMCFractions; /*!< Histogram for particle abundances from MC */
    TH3F*     fHistMCWeights;     /*!< Histogram for weight factors to re-weight MCabundances to data ones. */
    AliMCEvent* fMCEvent;   /*!< MC event */
    float fMultOrCent;      /*!< counted multiplicity or centrality class */
    int fNPartTypes;        /*!< number of particle species */
    int fNCentralities;     /*!< number of selected centrality classes */
    TaskState fbTaskStatus; /* controls internal status of class */
    SysFlag fFlag;          /*!< enum flag for systematic variation */
    bool fUseMultiplicity; /*!< switch to use multiplicity instead of centrality*/
    
    // functions
    std::string GetFunctionFromSysFlag(SysFlag flag); //!
    std::string GetSysVarFromSysFlag(SysFlag flag);   //!
    float GetMultFromCent(int CentBin) const;              //!
    std::vector<float> GetMultTupleFromCent(int CentBin) const; //!
    float GetMultFromCent(std::string cent);         //!
    float GetCentFromMult(float dMult);             //!
    void InitHistos();                                //!
    void LoadMeasuredFractions();                     //!
    void CountEventMult();                            //!
    int GetPartTypeNumber(ParticleType type);         //!
    int GetPartTypeNumber(std::string Particle);      //!
    int GetCentFromString(std::string cent);
    bool
    LoadFromAliMCSpectraWeight(AliMCSpectraWeights* obj);                 //!
    bool LoadFromTHnF(const char* histname);                              //!
    bool CalculateMCWeights();                                            //!
    bool CalcMCFractions();                                               //!
    bool CorrectFractionsforRest();                                       //!
     AliMCSpectraWeights(const AliMCSpectraWeights&);
     AliMCSpectraWeights& operator=(const AliMCSpectraWeights&);
public:
    AliMCSpectraWeights(); /*!< default root constructor */
    AliMCSpectraWeights(
                        std::string collisionSystem, std::string stName,
                        AliMCSpectraWeights::SysFlag flag); /*!< constructor to be used.*/

    ~AliMCSpectraWeights();
    
    void
    Init(); /*!< Function to start initalizing after all setters are made. */
    float GetMCSpectraWeight(
                             TParticle* mcGenParticle,
                             float eventMultiplicityOrCentrality); /*!< main function to use. Will deliver correct weights to re-weight the abundances of different particle species */
    float
    GetMCSpectraWeight(TParticle* mcGenParticle,
                       AliMCEvent* mcEvent); /*!< preferable to use this */
    void FillMCSpectra(
                       AliMCEvent* mcEvent); /*!< function to fill internal mc spectra for calculation of weight factors*/
    
    // Setter
    /** @fn void SetBinsPt(std::vector<double> bins)
     *  @brief function to set pt binning of all internal histograms
     *  @param bins a std::vector of floats containing the binning
     */
    void SetBinsPt(std::vector<float> bins);
    
    /** @fn void SetBinsMultCent(std::vector<double> bins)
     *  @brief function to set multiplicity binning of all internal histograms
     *  @param bins a std::vector of floats containing the binning
     */
    void SetBinsMultCent(std::vector<float> bins);
    
    void SetMCSpectraFile(const char* file) { fstFileMCSpectra = file; }
    void SetDataFractionsFile(const char* file) { fstFilePublished = file; }
    void SetCollisionSystem(const char* system) { fstCollisionSystem = system; }
    void SetUseMultiplicity(bool bMult) { fUseMultiplicity = bMult; }
    void SetSavedObjName(const char* name) { fstSavedObjName = name; }
    void SetSavedListName(const char* name) { fstSavedListName = name; }
    void SetSysFlag(SysFlag flag) { fFlag = flag; }
    
    // Getter
    const std::vector<float>& GetBinsPt() const { return fBinsPt; }
    const std::vector<std::string>& GetParticleTypes() const { return fstPartTypes; }
    int GetNPartTypes() const { return fNPartTypes; }
    int GetTaskStatus() const { return fbTaskStatus; }
    TH3F* GetHistMCGenPrimTrackParticles() const {
        return fHistMCGenPrimTrackParticle;
    }
    TH3F* GetHistDataFraction() const { return fHistDataFractions; }
    TH3F* GetHistMCFraction() const { return fHistMCFractions; }
    TH3F* GetHistMCWeights() const { return fHistMCWeights; }
    SysFlag GetSysFlag() const { return fFlag; }
    float GetMultOrCent() const { return fMultOrCent; }
    
    int IdentifyMCParticle(TParticle* mcParticle);
};

#endif /* __AliMCSpectraWeights__ */
