#ifndef ALIEMCALCORRECTIONCOMPONENT_H
#define ALIEMCALCORRECTIONCOMPONENT_H

// AliEmcalCorrectionComponent
//

#include <map>
// For std::pair
#include <utility>
#include <string>
#include <iostream>

// CINT can't handle the yaml header!
#if !(defined(__CINT__) || defined(__MAKECINT__))
#include <yaml-cpp/yaml.h>
#endif

#include <Rtypes.h>
#include <TNamed.h>
#include <TFile.h>

#include "AliAnalysisTaskEmcal.h"
#include "AliLog.h"
#include "AliEMCALRecoUtils.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODCaloCluster.h"
#include "AliEmcalParticle.h"
#include "AliEMCALGeometry.h"
#include "AliAnalysisManager.h"

class AliEmcalCorrectionComponent : public TNamed {
 public:
  enum containerType {
    kCluster = 0, //!< Cluster container
    kTrack = 1    //!< Track container
  };

  AliEmcalCorrectionComponent();
  AliEmcalCorrectionComponent(const char * name);
  virtual ~AliEmcalCorrectionComponent();

  // Virtual functions to be overloaded 
  virtual Bool_t Initialize();
  virtual void ExecOnce();
  virtual Bool_t Run();
  virtual Bool_t UserNotify();
  
  void GetEtaPhiDiff(const AliVTrack *t, const AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
  void UpdateCells();
  Bool_t RunChanged();
  void GetPass();
  void FillCellQA(TH1F* h);

  AliClusterContainer * GetClusterContainer() { return fClusCont; }
  AliParticleContainer * GetParticleContainer() { return fPartCont; }
  AliVCaloCells * GetCaloCells() { return fCaloCells; }
  TList * GetOutputList() { return fOutput; }
  
  void SetClusterContainer(AliClusterContainer * cont) { fClusCont = cont; }
  void SetParticleContainer(AliParticleContainer * cont) { fPartCont = cont; }
  void SetCaloCells(AliVCaloCells * cells) { fCaloCells = cells; }
  void SetRecoUtils(AliEMCALRecoUtils *ru) { fRecoUtils = ru; }

  void SetEvent(AliVEvent * event) { fEvent = event; }
  void SetMCEvent(AliMCEvent * mcevent) { fMCEvent = mcevent; }

  void SetEMCALGeometry(AliEMCALGeometry * geometry ) { fGeom = geometry; }
  void SetCentralityBin(Int_t bin) { fCentBin = bin; }
  void SetCentrality(Double_t cent) { fCent = cent; }
  void SetNcentralityBins(Int_t n) { fNcentBins = n; }

#if !(defined(__CINT__) || defined(__MAKECINT__))
  /// Make copy to ensure that the nodes do not point to each other (?)
  void SetUserConfiguration(YAML::Node & node) { fUserConfiguration = node; }
  void SetDefaultConfiguration(YAML::Node & node) { fDefaultConfiguration = node; }
  
  /// Retrieve property
  template<typename T> void GetProperty(std::string propertyName, T & property, bool requiredProperty = true, std::string correctionName = "");

  /// Retrieve property driver function. It is static so that it can be used by other classes
  template<typename T> static void GetProperty(std::string propertyName, T & property, const YAML::Node & userConfiguration, const YAML::Node & defaultConfiguration, bool requiredProperty = true, std::string correctionName = "");
#endif
  static bool IsSharedValue(std::string value);

  // Determine branch name using the "usedefault" pattern
  static std::string DetermineUseDefaultName(containerType contType, Bool_t esdMode);

 protected:

  void AddContainer(containerType type);

#if !(defined(__CINT__) || defined(__MAKECINT__))
  template<typename T> static bool GetPropertyFromNodes(const YAML::Node & node, const YAML::Node & sharedParametersNode, std::string propertyName, T & property, const std::string correctionName, const std::string configurationType, int nodesDeep = 0);
  template<typename T> static bool GetPropertyFromNode(const YAML::Node & node, std::string propertyName, T & property);

  YAML::Node              fUserConfiguration;
  YAML::Node              fDefaultConfiguration;
#endif

  Bool_t                  fCreateHisto;                   // flag to make some basic histograms
  Int_t                   fRun;                           //!run number
  TString                 fFilepass;                      // input data pass number
  Bool_t                  fGetPassFromFileName;           // get fFilepass from file name
  AliVEvent              *fEvent;                         // pointer to event
  Bool_t                  fEsdMode;                       // flag for ESD
  AliMCEvent             *fMCEvent;                       // MC
  Double_t                fCent;                          //!<!event centrality
  Int_t                   fNcentBins;                     ///< how many centrality bins (this member copied from AliAnalysisTaskEmcal)
  Int_t                   fCentBin;                       //!<!event centrality bin
  Int_t                   fNbins;                         ///< no. of pt bins
  Double_t                fMinBinPt;                      ///< min pt in histograms
  Double_t                fMaxBinPt;                      ///< max pt in histograms
  Double_t                fVertex[3];                     //!<!event vertex
  AliEMCALGeometry       *fGeom;                          //!<!geometry object
  Bool_t                  fIsEmbedded;                    ///< trigger, embedded signal
  Int_t                   fMinMCLabel;                    ///< minimum MC label value for the tracks/clusters being considered MC particles
  AliClusterContainer    *fClusCont;                      // pointer to the cluster container
  AliParticleContainer   *fPartCont;                      // pointer to the track/particle container
  AliVCaloCells          *fCaloCells;                     // pointer to calo cells
  AliEMCALRecoUtils      *fRecoUtils;                     // pointer to reco utils
  TList                  *fOutput;                        // list of output histograms
  
  TString                fBasePath;                       // base folder path to get root files

 private:
  AliEmcalCorrectionComponent(const AliEmcalCorrectionComponent &);               // Not implemented
  AliEmcalCorrectionComponent &operator=(const AliEmcalCorrectionComponent &);    // Not implemented
  
  ClassDef(AliEmcalCorrectionComponent, 1) // EMCal correction component
};

#if !(defined(__CINT__) || defined(__MAKECINT__))
/**
 *
 *
 */
template<typename T>
void AliEmcalCorrectionComponent::GetProperty(std::string propertyName, T & property, bool requiredProperty, std::string correctionName)
{
  // Get proper correction name
  if (correctionName == "")
  {
    correctionName = GetName();
  }
  GetProperty(propertyName, property, fUserConfiguration, fDefaultConfiguration, requiredProperty, correctionName);
}

/**
 *
 *
 */
template<typename T>
void AliEmcalCorrectionComponent::GetProperty(std::string propertyName, T & property, const YAML::Node & userConfiguration, const YAML::Node & defaultConfiguration, bool requiredProperty, std::string correctionName)
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
    // AliInfo
    std::cout << "INFO: Looking for parameter in user configuration" << std::endl;
    YAML::Node sharedParameters = userConfiguration["sharedParameters"];
    //std::cout << std::endl << "User Node: " << fUserConfiguration << std::endl;
    setProperty = AliEmcalCorrectionComponent::GetPropertyFromNodes(userConfiguration, sharedParameters, propertyName, property, correctionName, "user");
  }

  if (setProperty != true)
  {
    // AliInfo
    std::cout << "INFO: Looking for parameter in default configuration" << std::endl;
    YAML::Node sharedParameters = defaultConfiguration["sharedParameters"];
    //std::cout << std::endl << "Default Node: " << fDefaultConfiguration << std::endl;
    setProperty = AliEmcalCorrectionComponent::GetPropertyFromNodes(defaultConfiguration, sharedParameters, propertyName, property, correctionName, "default");
  }

  if (setProperty != true)
  {
    std::stringstream message;
    message << "Failed to retrieve property \""
        << (correctionName != "" ? correctionName + ":" : "")
        << propertyName << "\" from default config!" << std::endl;
    ::Fatal(correctionName.c_str(), "%s", message.str().c_str());
  }

}

/**
 *
 *
 */
template<typename T>
bool AliEmcalCorrectionComponent::GetPropertyFromNodes(const YAML::Node & node, const YAML::Node & sharedParametersNode, std::string propertyName, T & property, const std::string correctionName, const std::string configurationType, int nodesDeep)
{
  bool returnValue = false;
  if (nodesDeep > 2)
  {
    // Ensure that we do not go past two levels
    // AliFatal
    std::cout << "FATAL: Went too many levels of recursion. Bailing out on \"" << correctionName << ":" << propertyName << "\" from " << configurationType << " config at level " << nodesDeep << std::endl;
    //AliFatal(Form("%s: Went too many levels of recursion. Bailing out on \"%s:%s\" from %s config at lavel %i", GetName(), correctionName, propertyName, configurationType, nodesDeep)); 
    return false;
  }

  // Only want to print once to ensure the user is not overwhelmed!
  if (nodesDeep == 0)
  {
    // AliInfo
    std::cout << "INFO: Retreiving property \""
         << (correctionName != "" ? correctionName + ":" : "")
         << propertyName << "\" from " << configurationType << " config at level " << nodesDeep << std::endl;
  }

  if (node.IsDefined() == true)
  {
    //std::cout << "Node: " << node << std::endl;
    if (node[propertyName])
    {
      bool isShared = false;
      // Necessary because it fails on vectors and other complex objects that are YAML sequences.
      if (std::is_arithmetic<T>::value || std::is_same<T, std::string>::value)
      {
        // Retrieve value as string to check for shared value
        std::string testIfShared = node[propertyName].as<std::string>();
        // Check for a shared value
        //std::cout << "testIfShared: " << testIfShared << std::endl;
        isShared = AliEmcalCorrectionComponent::IsSharedValue(testIfShared);
        //std::cout << "Value shared? " << isShared << std::endl;
      }

      std::stringstream tempMessage;
      tempMessage << "property \"" 
            << (nodesDeep > 0 ? correctionName + ":" : "")
            << propertyName
            << "\" using " << ( isShared ? "shared " : "" )
            << "values from the " << configurationType
            << " configuration at level " << nodesDeep;

      // AliInfo
      std::cout << "INFO: Retrieveing " << tempMessage.str() << std::endl;

      bool retrievalResult = false;
      if (isShared == true)
      {
        retrievalResult = GetPropertyFromNode(sharedParametersNode, propertyName, property);
      }
      else
      {
        retrievalResult = GetPropertyFromNode(node, propertyName, property);
      }

      // Inform about the result
      if (retrievalResult == true)
      {
        // AliInfo
        std::cout << "INFO: Succeeded in retrieving " << tempMessage.str() << std::endl;
        returnValue = true;
      }
      else
      {
        // Only fatal if we have exhausted our last option, the default
        if (configurationType == "default")
        {
          // AliFatal
          std::cout << "FATAL: Failed to retrieve " << tempMessage.str() << std::endl;
        }
        else
        {
          // AliInfo
          std::cout << "INFO: Failed to retrieve " << tempMessage.str() << std::endl;
        }
        returnValue = false;
      }
    }
    else
    {
      // Go one node deeper using recursion
      // Must create a new node, since we took the original node by reference
      // AliInfo
      YAML::Node deeperNode = node[correctionName];
      nodesDeep++;
      std::cout << "INFO: Going a node deeper with " << correctionName << " to level " << nodesDeep << std::endl;
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
          // AliInfo
          std::cout << "INFO: Attempting to retrieve property \"" << propertyName << "\" from base correction \"" << subCorrectionName << "\" at level " << nodesDeep << std::endl;
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
    // AliInfo
    std::cout << "INFO: Node is undefined for \""
         << (correctionName != "" ? correctionName + ":" : "")
         << propertyName << "\" from " << configurationType << " config at level " << nodesDeep << std::endl;

    returnValue = false;
  }

  return returnValue;
}

/**
 *
 *
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

// Below is based on: https://stackoverflow.com/a/582456 , and editted for our purposes.
// 
// Template for creating a new component. Used to register the component.
template<typename T> AliEmcalCorrectionComponent * createT() { return new T; }

// Factory to create and keep track of new components
class AliEmcalCorrectionComponentFactory
{
 public:
  typedef std::map<std::string, AliEmcalCorrectionComponent*(*)()> map_type;

  // Creates an instance of an object based on the name if the name is registered in the map.
  static AliEmcalCorrectionComponent * createInstance(std::string const& s)
  {
    map_type::iterator it = getMap()->find(s);
    if(it == getMap()->end())
      return 0;
    // Initializes the function with empty arguments (?)
    return it->second();
  }

 protected:
  // Creates and access the component map
  static map_type * getMap() {
    // We never delete the map (until program termination) because we cannot guarantee
    // correct destruction order
    if(!componentMap) { componentMap = new map_type;  }
    return componentMap;
  }

 private:
  // Contains the map to all of the components
  static map_type * componentMap;
};

// Allows registration of new components into the factory map.
template<typename T>
class RegisterCorrectionComponent : public AliEmcalCorrectionComponentFactory
{ 
 public:
  RegisterCorrectionComponent(std::string const& s)
  { 
    getMap()->insert(std::make_pair(s, &createT<T>));
  }
};

#endif /* ALIEMCALCORRECTIONCOMPONENT_H */
