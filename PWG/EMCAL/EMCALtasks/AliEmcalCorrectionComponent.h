#ifndef ALIEMCALCORRECTIONCOMPONENT_H
#define ALIEMCALCORRECTIONCOMPONENT_H

#include <map>
#include <string>

class TH1F;
#include <TNamed.h>

class AliMCEvent;
class AliEMCALRecoUtils;
class AliVCaloCells;
class AliVTrack;
class AliVCluster;
class AliVEvent;
#include <AliLog.h>
#include <AliEMCALGeometry.h>
#include "AliYAMLConfiguration.h"
#include "AliEmcalContainerUtils.h"
#include "AliParticleContainer.h"
#include "AliMCParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliClusterContainer.h"
#include "AliEmcalCorrectionEventManager.h"

/**
 * @class AliEmcalCorrectionComponent
 * @ingroup EMCALCORRECTIONFW
 * @brief Base class for correction components in the EMCal correction framework
 *
 * Base class for all correction components in the EMCal Correction Framework. Each correction
 * component corresponds to a correction needed for the EMCal. Creation, configuration, and execution
 * of all of the components is handled by AliEmcalCorrectionTask. Each new component is automatically
 * registered through the AliEmcalCorrectionComponentFactory class, and is thus automatically available
 * to execute. Components are configured through a set of %YAML configuration files stored in
 * AliYAMLConfiguration. For more information about the steering, see AliEmcalCorrectionTask.
 *
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University
 * @date Jul 8, 2016
 */

class AliEmcalCorrectionComponent : public TNamed {
 public:
  AliEmcalCorrectionComponent();
  AliEmcalCorrectionComponent(const char * name);
  virtual ~AliEmcalCorrectionComponent();

  // Virtual functions to be overloaded 
  virtual Bool_t Initialize();
  virtual void UserCreateOutputObjects();
  virtual void ExecOnce();
  virtual Bool_t Run();
  virtual Bool_t UserNotify();
  virtual Bool_t CheckIfRunChanged();
  
  void GetEtaPhiDiff(const AliVTrack *t, const AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
  void UpdateCells();
  void GetPass();
  void FillCellQA(TH1F* h);
  Int_t InitBadChannels();

  // Containers and cells
  AliParticleContainer   *AddParticleContainer(const char *n)                    { return AliEmcalContainerUtils::AddContainer<AliParticleContainer>(n, fParticleCollArray); }
  AliTrackContainer      *AddTrackContainer(const char *n)                       { return AliEmcalContainerUtils::AddContainer<AliTrackContainer>(n, fParticleCollArray); }
  AliMCParticleContainer *AddMCParticleContainer(const char *n)                  { return AliEmcalContainerUtils::AddContainer<AliMCParticleContainer>(n, fParticleCollArray); }
  AliClusterContainer    *AddClusterContainer(const char *n)                     { return AliEmcalContainerUtils::AddContainer<AliClusterContainer>(n, fClusterCollArray); }
  void                    AdoptParticleContainer(AliParticleContainer* cont)     { fParticleCollArray.Add(cont)                        ; }
  void                    AdoptTrackContainer(AliTrackContainer* cont)           { AdoptParticleContainer(cont)                        ; }
  void                    AdoptMCParticleContainer(AliMCParticleContainer* cont) { AdoptParticleContainer(cont)                        ; }
  void                    AdoptClusterContainer(AliClusterContainer* cont)       { fClusterCollArray.Add(cont)                         ; }
  AliParticleContainer   *GetParticleContainer(Int_t i=0)                  const { return AliEmcalContainerUtils::GetContainer<AliParticleContainer>(i, fParticleCollArray); }
  AliParticleContainer   *GetParticleContainer(const char* name)           const { return AliEmcalContainerUtils::GetContainer<AliParticleContainer>(name, fParticleCollArray); }
  AliClusterContainer    *GetClusterContainer(Int_t i=0)                   const { return AliEmcalContainerUtils::GetContainer<AliClusterContainer>(i, fClusterCollArray); }
  AliClusterContainer    *GetClusterContainer(const char* name)            const { return AliEmcalContainerUtils::GetContainer<AliClusterContainer>(name, fClusterCollArray); }
  AliMCParticleContainer *GetMCParticleContainer(Int_t i=0)                const { return dynamic_cast<AliMCParticleContainer*>(GetParticleContainer(i))   ; }
  AliMCParticleContainer *GetMCParticleContainer(const char* name)         const { return dynamic_cast<AliMCParticleContainer*>(GetParticleContainer(name)); }
  AliTrackContainer      *GetTrackContainer(Int_t i=0)                     const { return dynamic_cast<AliTrackContainer*>(GetParticleContainer(i))        ; }
  AliTrackContainer      *GetTrackContainer(const char* name)              const { return dynamic_cast<AliTrackContainer*>(GetParticleContainer(name))     ; }
  void                    RemoveParticleContainer(Int_t i=0)                     { fParticleCollArray.RemoveAt(i)                      ; }
  void                    RemoveClusterContainer(Int_t i=0)                      { fClusterCollArray.RemoveAt(i)                       ; }
  AliEMCALRecoUtils      *GetRecoUtils()  const { return fRecoUtils; }
  AliVCaloCells          *GetCaloCells()  const { return fCaloCells; }
  TList                  *GetOutputList() const { return fOutput; }
  
  void SetCaloCells(AliVCaloCells * cells) { fCaloCells = cells; }
  void SetRecoUtils(AliEMCALRecoUtils *ru) { fRecoUtils = ru; }

  void SetInputEvent(AliVEvent * event) { fEventManager.SetInputEvent(event); }
  void SetMCEvent(AliMCEvent * mcevent) { fMCEvent = mcevent; }
  /**
   * If we are using standard input event then the embedded event should not be used!
   * We store whether the embedding event should be used, so we invert the bool here.
   * Then, if it is only set when we see the standard input event, then any single input
   * with the standard input event will be enough to disable the embedded event.
   */
  void SetUsingInputEvent(bool b = true) { fEventManager.SetUseEmbeddingEvent(!b); }

  void SetEMCALGeometry(AliEMCALGeometry * geometry ) { fGeom = geometry; }
  void SetCentralityBin(Int_t bin) { fCentBin = bin; }
  void SetCentrality(Double_t cent) { fCent = cent; }
  void SetNcentralityBins(Int_t n) { fNcentBins = n; }
  void SetVertex(Double_t * vertex) { fVertex[0] = vertex[0]; fVertex[1] = vertex[1]; fVertex[2] = vertex[2]; }
  void SetIsESD(Bool_t isESD) {fEsdMode = isESD; }
  void SetCustomBadChannels(TString customBC) {fCustomBadChannelFilePath = customBC; }

  /// Set %YAML Configuration
  void SetYAMLConfiguration(PWG::Tools::AliYAMLConfiguration config) { fYAMLConfig = config; }

  /// Retrieve property
  template<typename T> bool GetProperty(std::string propertyName, T & property, bool requiredProperty = true, std::string correctionName = "");
 protected:
  PWG::Tools::AliYAMLConfiguration fYAMLConfig;           ///< Contains the %YAML configuration used to configure the component
  Bool_t                  fCreateHisto;                   ///< Flag to make some basic histograms
  Bool_t                  fLoad1DBadChMap;                ///< Flag to load 1D bad channel map
  Int_t                   fRun;                           //!<! Run number
  TString                 fFilepass;                      ///< Input data pass number
  Bool_t                  fGetPassFromFileName;           ///< Get fFilepass from file name
  AliEmcalCorrectionEventManager fEventManager;           ///< Minimal task which inherits from AliAnalysisTaskSE and manages access to the event
  Bool_t                  fEsdMode;                       ///< flag for ESD
  AliMCEvent             *fMCEvent;                       //!<! MC
  Double_t                fCent;                          //!<! Event centrality
  Int_t                   fNcentBins;                     ///< How many centrality bins (this member copied from AliAnalysisTaskEmcal)
  Int_t                   fCentBin;                       //!<! Event centrality bin
  Int_t                   fNbins;                         ///< No. of pt bins
  Double_t                fMinBinPt;                      ///< Min pt in histograms
  Double_t                fMaxBinPt;                      ///< Max pt in histograms
  Double_t                fVertex[3];                     //!<! Event vertex
  AliEMCALGeometry       *fGeom;                          //!<! Geometry object
  Int_t                   fMinMCLabel;                    ///< Minimum MC label value for the tracks/clusters being considered MC particles
  TObjArray               fClusterCollArray;              ///< Cluster collection array
  TObjArray               fParticleCollArray;             ///< Particle/track collection array
  AliVCaloCells          *fCaloCells;                     //!<! Pointer to CaloCells
  AliEMCALRecoUtils      *fRecoUtils;                     ///<  Pointer to RecoUtils
  TList                  *fOutput;                        //!<! List of output histograms
  
  TString                fBasePath;                       ///< Base folder path to get root files
  TString                fCustomBadChannelFilePath;       ///< Custom path to bad channel map OADB file

 private:
  AliEmcalCorrectionComponent(const AliEmcalCorrectionComponent &);               // Not implemented
  AliEmcalCorrectionComponent &operator=(const AliEmcalCorrectionComponent &);    // Not implemented
  
  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionComponent, 9); // EMCal correction component
  /// \endcond
};

/**
 * Get the requested property from the %YAML configuration. This function is generally used by
 * AliEmcalCorrectionComponent derived tasks as a wrapper around the more complicated functions to
 * retrieve properties.
 *
 * @param[in] propertyName Name of the property to retrieve
 * @param[out] property Containers the retrieved property
 * @param[in] requiredProperty True if the property is required
 * @param[in] correctionName Name of the correction from where the property should be retrieved
 *
 * @return True if the property was set successfully
 */
template<typename T>
bool AliEmcalCorrectionComponent::GetProperty(std::string propertyName, T & property, bool requiredProperty, std::string correctionName)
{
  // Get proper correction name
  if (correctionName == "")
  {
    correctionName = GetName();
  }
  // Add the correction name as a prefix of the property
  return fYAMLConfig.GetProperty(std::vector<std::string>{correctionName, propertyName}, property, requiredProperty);
}

/**
 * @class AliEmcalCorrectionComponentFactory
 * @ingroup EMCALCOREFW
 * @brief Factory for correction components in the EMCal correction framework
 *
 * This class maintains a map between the name of the correction component and a function to create the
 * component. This map can be then be used to create each desired component by passing the name in a string.
 * The benefit of this approach is new components can be automatically registered without changing any of the
 * correction classes. Only the %YAML configuration needs to be changed!
 *
 * The class and associated functions are based on: https://stackoverflow.com/a/582456 , and edited
 * for our purposes.
 *
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University
 * @date Jul 8, 2016
 */

/// Template function for creating a new component. Used to register the component.
template<typename T> AliEmcalCorrectionComponent * createT() { return new T; }

// Factory to create and keep track of new components
class AliEmcalCorrectionComponentFactory
{
 public:
  virtual ~AliEmcalCorrectionComponentFactory() {}

  typedef std::map<std::string, AliEmcalCorrectionComponent*(*)()> map_type;

  /// Creates an instance of an object based on the name if the name is registered in the map.
  static AliEmcalCorrectionComponent * createInstance(std::string const& s)
  {
    map_type::iterator it = getMap()->find(s);
    if(it == getMap()->end())
      return 0;
    // Initializes the function with empty arguments (?)
    return it->second();
  }

 protected:
  /// Creates and access the component map
  static map_type * getMap() {
    // We never delete the map (until program termination) because we cannot guarantee
    // correct destruction order
    if(!componentMap) { componentMap = new map_type;  }
    return componentMap;
  }

 private:
  /// Contains the map to all of the components
  static map_type * componentMap;
};

/**
 * @class RegisterCorrectionComponent
 * @ingroup EMCALCOREFW
 * @brief Registers EMCal correction components in the factory map
 *
 * This class allows EMCal correction components to automatically register in a map, such that new components
 * are automatically available in the correction framework.
 *
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University
 * @date Jul 8, 2016
 */
template<typename T>
class RegisterCorrectionComponent : public AliEmcalCorrectionComponentFactory
{ 
 public:
  /// Registers the name of the component to map to a function that can create the component
  RegisterCorrectionComponent(std::string const& s)
  { 
    getMap()->insert(std::make_pair(s, &createT<T>));
  }
};

#endif /* ALIEMCALCORRECTIONCOMPONENT_H */
