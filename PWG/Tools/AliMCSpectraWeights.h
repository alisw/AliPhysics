
/**
 * @file AliMCSpectraWeights.h
 * @brief Class for re-weighting particle abundances in MC simulations
 * @author Patrick Huhn
 * @date 25.10.2019
 */
#ifndef __AliMCSpectraWeights__
#define __AliMCSpectraWeights__

//#define __AliMCSpectraWeights_DebugPCC__
//#define __AliMCSpectraWeights_DebugTiming__

class TParticle;
class AliMCEvent;
class TH3F;
class TH1D;
class TF1;
#include "TNamed.h"
#include "TRandom3.h"
#include <map>
#include <string>
#include <vector>



/**
 * @class AliMCSpectraWeights AliMCSpectraWeights.h "AliMCSpectraWeights.h"
 * @brief TODO
 *
 *  TODO write something up
 *  TODO write full explanation here
 */
class AliMCSpectraWeights : public TNamed {
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
        kBylinkin,
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
    std::string
        fstFileMCSpectra; /*!< path to previous train output of MC fractions */
    std::string fstFilePublished; /*!< path to calculated fractions from
                                     published spectra */
    std::string fstSavedObjName;
    std::string fstSavedListName;
    std::vector<std::string>
        fstPartTypes; /*!< Array of used particle species */
    std::vector<std::string> fstCentralities;
    std::vector<float> fBinsMultCent; /*!< centrality or multiplicity binning */
    std::vector<float> fBinsPt;       /*!< pT binning */
    std::vector<AliMCSpectraWeights::SysFlag> fAllSystematicFlags;
    std::map<AliMCSpectraWeights::SysFlag, TH3F*>
        fHistMCWeightsSys; /*!< Histograms for systematic variations of weight
                              factors */
    TRandom3 frndGen;
    TH3F* fHistMCGenPrimTrackParticle; /*!< Histogram for MC particle
                                          information*/
    TH3F* fHistDataFractions; /*!< Histogram for particle abundances from
                                 published data */
    TH3F* fHistMCFractions;   /*!< Histogram for particle abundances from MC */
    TH3F* fHistMCWeights;     /*!< Histogram for weight factors to re-weight
                                 MCabundances to data ones. */
    AliMCEvent* fMCEvent;     /*!< MC event */
    int         fNPrimaries;   /*!< MC primaries from last event*/
    float fMultOrCent;        /*!< counted multiplicity or centrality class */
    int fNPartTypes;          /*!< number of particle species */
    int fNCentralities;       /*!< number of selected centrality classes */
    int fNSysFlags;           /*!< number of all systematic flags */
    TaskState fbTaskStatus;   /*!< controls internal status of class */
    SysFlag fFlag;            /*!< enum flag for systematic variation */
    bool fUseMultiplicity;    /*!< switch to use multiplicity instead of
                                 centrality*/
    bool fUseMBFractions;     /*!< switch to use MB fractions instead of mult.
                                 dependent ones*/
    bool fDoSystematics;

    // functions
    // intern getter
    std::string const GetFunctionFromSysFlag(SysFlag flag) const;           //!
    std::string const GetSysVarFromSysFlag(SysFlag flag) const;             //!
    float const GetMultFromCent(int CentBin) const;                   //!
    std::vector<float> GetMultTupleFromCent(int CentBin) const; //!
    float const GetMultFromCent(std::string const& cent) const;                    //!
    float const GetCentFromMult(float const dMult) const;                         //!
    int const GetPartTypeNumber(ParticleType type) const;                   //!
    int const GetPartTypeNumber(std::string const& Particle) const;                //!
    int const GetCentFromString(std::string const& cent) const;
    // main functions
    void InitHistos();                                          //!
    void LoadMeasuredFractions();                               //!
    void CountEventMult();                                      //!
    void SelectRndSysFlagForEvent();                                 //!
    bool CalculateMCWeights();                                 //!
    bool CalcMCFractions();                                    //!
    bool CorrectFractionsforRest();                            //!

    int const CheckAndIdentifyParticle(TParticle* part);
    int const FindBinEntry(float pt, int const part);
    
    // private = to be deleted
    AliMCSpectraWeights(const AliMCSpectraWeights&);//copy
    AliMCSpectraWeights& operator=(const AliMCSpectraWeights&);//copy assign

  public:
    AliMCSpectraWeights(); /*!< default root constructor */
    AliMCSpectraWeights(
        std::string const &collisionSystem, std::string const &stName,
        AliMCSpectraWeights::SysFlag flag); /*!< constructor to be used.*/
//    AliMCSpectraWeights(const AliMCSpectraWeights& org); // copy constructor
//    AliMCSpectraWeights(AliMCSpectraWeights* org); // copy constructor
    ~AliMCSpectraWeights();

    void Init(); /*!< Function to start initalizing after all setters are made. */
    float const
    GetMCSpectraWeight(TParticle* mcGenParticle,
                       AliMCEvent* mcEvent); /*!< old; should not be used */
    float const
    GetMCSpectraWeightNominal(TParticle* mcGenParticle);/*!< main function to use. Will
                                                         deliver correct weights to
                                                         re-weight the abundances of
                                                         different particle species */
    float const
    GetMCSpectraWeightSystematics(TParticle* mcGenParticle);

    void FillMCSpectra(
        AliMCEvent* mcEvent); /*!< function to fill internal mc spectra for
                                 calculation of weight factors*/
    void StartNewEvent();

    // Setter
    void SetMCSpectraFile(std::string const& file) { fstFileMCSpectra = file; }
    void SetDataFractionsFile(std::string const& file) { fstFilePublished = file; }
    void SetCollisionSystem(std::string const& system) { fstCollisionSystem = system; }
    void SetUseMultiplicity(bool bMult) { fUseMultiplicity = bMult; }
    void SetUseMBFractions(bool bMBFractions = true) {fUseMBFractions = bMBFractions;}
    void SetSavedObjName(const char* name) { fstSavedObjName = name; }
    void SetSavedListName(const char* name) { fstSavedListName = name; }
    void SetSysFlag(SysFlag flag) { fFlag = flag; }
    void SetCurrentEvent(AliMCEvent* event) {
        if (event != fMCEvent) {
            fMCEvent = event;
            CountEventMult();
        }
    }
    void SetDoSystematics(bool doSys = true) { fDoSystematics = doSys; }

    // Getter
    std::vector<std::string> const GetParticleTypes() const {return fstPartTypes;}
    int const GetNPartTypes() const { return fNPartTypes; }
    int const GetTaskStatus() const { return fbTaskStatus; }
    TH3F* const GetHistMCGenPrimTrackParticles() const { return fHistMCGenPrimTrackParticle;}
    TH3F* const GetHistDataFraction() const { return fHistDataFractions; }
    TH3F* const GetHistMCFraction() const { return fHistMCFractions; }
    TH3F* const GetHistMCWeights() const { return fHistMCWeights; }
    std::map<SysFlag, TH3F*> const GetHistMCWeightsSys() const {return fHistMCWeightsSys;}
    SysFlag const GetSysFlag() const { return fFlag; }
    float const GetMultOrCent() const { return fMultOrCent; }

    int const IdentifyMCParticle(TParticle* mcParticle);

    ClassDef(AliMCSpectraWeights, 1);
};

struct AliMCSpectraWeightsHandler : public TNamed {
    AliMCSpectraWeightsHandler(); //default for ROOT
    AliMCSpectraWeightsHandler(AliMCSpectraWeights* fMCWeight, const char* name); // to be used
    ~AliMCSpectraWeightsHandler(){}

    AliMCSpectraWeights* fMCSpectraWeight = 0;
private:
    AliMCSpectraWeightsHandler(const AliMCSpectraWeightsHandler&);//copy
    AliMCSpectraWeightsHandler& operator=(const AliMCSpectraWeightsHandler&);//copy assign

    ClassDef(AliMCSpectraWeightsHandler, 1);
};

#endif /* __AliMCSpectraWeights__ */
