#ifndef ALIEMCALCORRECTIONCOMPONENT_H
#define ALIEMCALCORRECTIONCOMPONENT_H

#include <map>
#include <string>

// CINT can't handle the yaml header!
#if !(defined(__CINT__) || defined(__MAKECINT__))
#include <yaml-cpp/yaml.h>
#endif

class TH1F;
#include <TNamed.h>

class AliMCEvent;
class AliEMCALRecoUtils;
class AliVCaloCells;
class AliVTrack;
class AliVCluster;
class AliVEvent;
#include <AliLog.h>
#include "AliEmcalContainerUtils.h"
#include "AliParticleContainer.h"
#include "AliMCParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliClusterContainer.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalCorrectionEventManager.h"

/**
 * @class AliEmcalCorrectionComponent
 * @ingroup EMCALCOREFW
 * @brief Base class for correction components in the EMCal correction framework
 *
 * Base class for all correction components in the EMCal Correction Framework. Each correction
 * component corresponds to a correction needed for the EMCal. Creation, configuration, and execution
 * of all of the components is handled by AliEmcalCorrectionTask. Each new component is automatically
 * registered through the AliEmcalCorrectionComponentFactory class, and is thus automatically available
 * to execute. Components are configured through a set of YAML configuration files. For more information
 * about the steering, see AliEmcalCorrectionTask.
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

#if !(defined(__CINT__) || defined(__MAKECINT__))
  /// Make copy to ensure that the nodes do not point to each other (?)
  void SetUserConfiguration(YAML::Node & node) { fUserConfiguration = node; }
  void SetDefaultConfiguration(YAML::Node & node) { fDefaultConfiguration = node; }
  
  /// Retrieve property
  template<typename T> bool GetProperty(std::string propertyName, T & property, bool requiredProperty = true, std::string correctionName = "");

  /// Retrieve property driver function. It is static so that it can be used by other classes
  template<typename T> static bool GetProperty(std::string propertyName, T & property, const YAML::Node & userConfiguration, const YAML::Node & defaultConfiguration, bool requiredProperty = true, std::string correctionName = "");
#endif
  static bool IsSharedValue(std::string & value);

 protected:

#if !(defined(__CINT__) || defined(__MAKECINT__))
  template<typename T> static bool GetPropertyFromNodes(const YAML::Node & node, const YAML::Node & sharedParametersNode, std::string propertyName, T & property, const std::string correctionName, const std::string configurationType, int nodesDeep = 0);
  template<typename T> static bool GetPropertyFromNode(const YAML::Node & node, std::string propertyName, T & property);

  template<typename T> static typename std::enable_if<!std::is_arithmetic<T>::value && !std::is_same<T, std::string>::value && !std::is_same<T, bool>::value>::type PrintRetrievedPropertyValue(T & property, std::stringstream & tempMessage);
  template<typename T> static typename std::enable_if<std::is_arithmetic<T>::value || std::is_same<T, std::string>::value || std::is_same<T, bool>::value>::type PrintRetrievedPropertyValue(T & property, std::stringstream & tempMessage);

  YAML::Node              fUserConfiguration;             //!<! User YAML configuration
  YAML::Node              fDefaultConfiguration;          //!<! Default YAML configuration
#endif

  Bool_t                  fCreateHisto;                   ///< Flag to make some basic histograms
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
  Bool_t                  fIsEmbedded;                    ///< Trigger, embedded signal
  Int_t                   fMinMCLabel;                    ///< Minimum MC label value for the tracks/clusters being considered MC particles
  TObjArray               fClusterCollArray;              ///< Cluster collection array
  TObjArray               fParticleCollArray;             ///< Particle/track collection array
  AliVCaloCells          *fCaloCells;                     //!<! Pointer to CaloCells
  AliEMCALRecoUtils      *fRecoUtils;                     ///<  Pointer to RecoUtils
  TList                  *fOutput;                        //!<! List of output histograms
  
  TString                fBasePath;                       ///< Base folder path to get root files

 private:
  AliEmcalCorrectionComponent(const AliEmcalCorrectionComponent &);               // Not implemented
  AliEmcalCorrectionComponent &operator=(const AliEmcalCorrectionComponent &);    // Not implemented
  
  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionComponent, 3); // EMCal correction component
  /// \endcond
};

#if !(defined(__CINT__) || defined(__MAKECINT__))
/**
 * Get the requested property from the YAML configuration. This function is generally used by
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
  bool result = GetProperty(propertyName, property, fUserConfiguration, fDefaultConfiguration, requiredProperty, correctionName);

  return result;
}

/**
 * General function to get a property from a set of YAML configuration files. It first calls checks
 * the user configuration, and then if the property is not found, it then checks the default configuration.
 * If the property is required but it is not found, a fatal error is thrown.
 *
 * @param[in] propertyName Name of the property to retrieve
 * @param[out] property Contains the retrieved property
 * @param[in] userConfiguration YAML user configuration node
 * @param[in] defaultConfiguration YAML default configuration node
 * @param[in] requiredProperty True if the property is required
 * @param[in] correctionName Name of the correction from where the property should be retrieved
 *
 * @return True if the property was set successfully
 */
template<typename T>
bool AliEmcalCorrectionComponent::GetProperty(std::string propertyName, T & property, const YAML::Node & userConfiguration, const YAML::Node & defaultConfiguration, bool requiredProperty, std::string correctionName)
{
  // Remove AliEmcalCorrection if in name
  std::size_t prefixStringLocation = correctionName.find("AliEmcalCorrection");
  if (prefixStringLocation != std::string::npos)
  {
    // AliEmcalCorrection is 18 characters
    correctionName.erase(prefixStringLocation, prefixStringLocation + 18);
  }

  bool setProperty = false;
  // IsNull checks is a node is empty. A node is empty if it is created.
  // IsDefined checks if the node that was requested was not actually created.
  // In this case, the user configuration node is always created. If it IsNull, then we ignore it.
  if (userConfiguration.IsNull() != true)
  {
    AliDebugClass(2, "Looking for parameter in user configuration");
    YAML::Node sharedParameters = userConfiguration["sharedParameters"];
    //std::cout << std::endl << "User Node: " << fUserConfiguration << std::endl;
    setProperty = AliEmcalCorrectionComponent::GetPropertyFromNodes(userConfiguration, sharedParameters, propertyName, property, correctionName, "user");
  }

  if (setProperty != true)
  {
    AliDebugClass(2, "Looking for parameter in default configuration");
    YAML::Node sharedParameters = defaultConfiguration["sharedParameters"];
    //std::cout << std::endl << "Default Node: " << fDefaultConfiguration << std::endl;
    setProperty = AliEmcalCorrectionComponent::GetPropertyFromNodes(defaultConfiguration, sharedParameters, propertyName, property, correctionName, "default");
  }

  if (setProperty != true && requiredProperty == true)
  {
    std::stringstream message;
    message << "Failed to retrieve property \""
        << (correctionName != "" ? correctionName + ":" : "")
        << propertyName << "\" from default config!" << std::endl;
    AliFatalClass(message.str().c_str());
  }

  // Return whether the value was actually set
  return setProperty;
}

/**
 * Actually handles retrieving parameters from YAML nodes.
 *
 * The retrieval procedure is as follows:
 *  - Check if the property is defined in the base of the node
 *  - If the property is found:
 *     - Check whether it requests a shared parameter
 *        - If so, attempt to retrieve the property from the shared parameter.
 *     - If not, attempt to retrieve the parameter from the correction name
 *     - In either case, through a warning (or fatal if requested) if it is not found
 *  - If not found, attempt to find a node named by the correction name:
 *     - If it exists, recurse in this function. It will not go more than two levels deep!
 *     - If it does not exist, attempt to find a substring name by removing all characters after the first "_"
 *         - Recurse if found
 *
 * @param[in] node Main YAML node containing the properties
 * @param[in] sharedParametersNode YAML node containing the shared parameters
 * @param[in] propertyName Name of the property to retrieve
 * @param[out] property Contains the retrieved property
 * @param[in] correctionName Name of the correction from where the property should be retrieved
 * @param[in] configurationType Name of the configuration type. Either "user" or "default"
 * @param[in] nodesDeep Keeps track of how many nodes deep we have recursed. Optional parameter that the user does not need to set.
 *
 * @return True if the property was set successfully
 */
template<typename T>
bool AliEmcalCorrectionComponent::GetPropertyFromNodes(const YAML::Node & node, const YAML::Node & sharedParametersNode, std::string propertyName, T & property, const std::string correctionName, const std::string configurationType, int nodesDeep)
{
  // Used as a buffer for printing complicated messages
  std::stringstream tempMessage;

  bool returnValue = false;
  if (nodesDeep > 2)
  {
    // Ensure that we do not go past two levels
    tempMessage.str("");
    tempMessage << "Went too many levels of recursion. Bailing out on \"" << correctionName << ":" << propertyName << "\" from " << configurationType << " config at level " << nodesDeep;
    AliDebugClass(2, TString::Format("%s", tempMessage.str().c_str()));
    return false;
  }

  // Only want to print once to ensure the user is not overwhelmed!
  if (nodesDeep == 0)
  {
    tempMessage.str("");
    tempMessage << "Retreiving property \""
          << (correctionName != "" ? correctionName + ":" : "")
          << propertyName << "\" from " << configurationType << " config at level " << nodesDeep;
    AliDebugClass(1, TString::Format("%s", tempMessage.str().c_str()));
  }

  if (node.IsDefined() == true)
  {
    //std::cout << "Node: " << node << std::endl;
    if (node[propertyName])
    {
      bool isShared = false;
      // Will contain the name of the property in the sharedParameters section that should be retrieved
      // if it is requested through the user configuration.
      std::string sharedValueName = "";
      // Necessary because it fails on vectors and other complex objects that are YAML sequences.
      if (std::is_arithmetic<T>::value || std::is_same<T, std::string>::value || std::is_same<T, bool>::value)
      {
        // Retrieve value as string to check for shared value
        sharedValueName = node[propertyName].as<std::string>();
        // Check for a shared value
        isShared = AliEmcalCorrectionComponent::IsSharedValue(sharedValueName);
      }

      tempMessage.str("");
      tempMessage << "property \"" 
            << (nodesDeep > 0 ? correctionName + ":" : "")
            << propertyName
            << "\" using " << ( isShared ? "\"sharedParameters:" + sharedValueName + "\" in " : "" )
            << "values from the " << configurationType
            << " configuration at level " << nodesDeep;

      AliDebugClass(2, TString::Format("Retrieveing %s", tempMessage.str().c_str()));

      bool retrievalResult = false;
      if (isShared == true)
      {
        retrievalResult = GetPropertyFromNode(sharedParametersNode, sharedValueName, property);
      }
      else
      {
        retrievalResult = GetPropertyFromNode(node, propertyName, property);
      }

      // Inform about the result
      if (retrievalResult == true)
      {
        // Add the retrieved value to the message (only if trivially printable)
        PrintRetrievedPropertyValue(property, tempMessage);
        AliDebugClassStream(2) << "Succeeded in retrieveing " << tempMessage.str() << std::endl;
        returnValue = true;
      }
      else
      {
        // Only fatal if we have exhausted our last option, the default
        if (configurationType == "default")
        {
          AliFatalClass(TString::Format("Failed to retrieve %s", tempMessage.str().c_str()));
        }
        else
        {
          AliDebugClass(2, TString::Format("Failed to retrieve %s", tempMessage.str().c_str()));
        }
        returnValue = false;
      }
    }
    else
    {
      // Go one node deeper using recursion
      // Must create a new node, since we took the original node by reference
      YAML::Node deeperNode = node[correctionName];
      nodesDeep++;
      AliDebugClass(2, TString::Format("Going a node deeper with \"%s\" to level %i", correctionName.c_str(), nodesDeep));
      returnValue = GetPropertyFromNodes(deeperNode, sharedParametersNode, propertyName, property, correctionName, configurationType, nodesDeep);

      // If we didn't find it, next try a substring
      if (returnValue == false)
      {
        // Don't increment nodesDeep, because we are about to go another node deep anyway
        std::size_t splitLocation = correctionName.find("_");
        // We only want to look at the split if we are only one node deep and
        // if the split is actually meaningful
        if (splitLocation != std::string::npos && nodesDeep == 1)
        {
          // Split from start to underscore
          std::string subCorrectionName = correctionName.substr(0, splitLocation);
          AliDebugClass(2, TString::Format("Attempting to retrieve property \"%s\" from base correction \"%s\" at level %i", propertyName.c_str(), subCorrectionName.c_str(), nodesDeep));
          // Retrieve the base correction node
          // Need to create new node! Otherwise it will assign the correctionName node to subCorrectionName node!
          YAML::Node subCorrectionNode = node[subCorrectionName];
          returnValue = GetPropertyFromNodes(subCorrectionNode, sharedParametersNode, propertyName, property, subCorrectionName, configurationType, nodesDeep);
        }
      }
    }
  }
  else
  {
    tempMessage.str("");
    tempMessage << "Node is undefined for \""
          << (correctionName != "" ? correctionName + ":" : "")
          << propertyName << "\" from " << configurationType << " config at level " << nodesDeep;
    AliDebugClass(2, TString::Format("%s", tempMessage.str().c_str()));

    returnValue = false;
  }

  return returnValue;
}

/**
 * Prints the retrieved property value if the type is easily printable. This function handles arithmetic
 * values, strings, and bools. Further types could be handled by other functions that are selectively enabled
 * by std::enable_if<>.
 *
 * NOTE: This function retunrs void! See: https://stackoverflow.com/a/29044828
 *
 * @param[in] property Property to be printed
 * @param[in, out] tempMessage Stringstream into which the property should be streamed
 */
template<typename T>
typename std::enable_if<std::is_arithmetic<T>::value || std::is_same<T, std::string>::value || std::is_same<T, bool>::value>::type
AliEmcalCorrectionComponent::PrintRetrievedPropertyValue(T & property, std::stringstream & tempMessage)
{
  //AliDebugClassStream(2) << " with value " << property;
  tempMessage << " with value \"" << property << "\"";
}

/**
 * Handles all other cases where the property cannot be printed.
 *
 * NOTE: This function retunrs void! See: https://stackoverflow.com/a/29044828
 *
 * @param[in] property Property to be printed
 * @param[in, out] tempMessage Stringstream into which the property should be streamed
 */
template<typename T>
typename std::enable_if<!std::is_arithmetic<T>::value && !std::is_same<T, std::string>::value && !std::is_same<T, bool>::value>::type
AliEmcalCorrectionComponent::PrintRetrievedPropertyValue(T & property, std::stringstream & tempMessage)
{
  // Cannot easily print these types, so just note that is the case!
  tempMessage << " with a value that cannot be trivially printed";
}

/**
 * Utility function to retrieve the property from a already selected YAML node.
 *
 * @param[in] node YAML node from which the property should be retrieved
 * @param[in] propertyName Name of the property to retrieve
 * @param[out] property Contains the retrieved property
 */
template<typename T>
bool AliEmcalCorrectionComponent::GetPropertyFromNode(const YAML::Node & node, std::string propertyName, T & property)
{
  if (node[propertyName])
  {
    property = node[propertyName].as<T>();
    return true;
  }
  return false;
}

#endif /* Hide yaml from CINT */

/**
 * @class AliEmcalCorrectionComponentFactory
 * @ingroup EMCALCOREFW
 * @brief Factory for correction components in the EMCal correction framework
 *
 * This class maintains a map between the name of the correction component and a function to create the
 * component. This map can be then be used to create each desired component by passing the name in a string.
 * The benefit of this approach is new components can be automatically registered without changing any of the
 * correction classes. Only the YAML configuration needs to be changed!
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
