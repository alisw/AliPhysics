#ifndef ALIYAMLCONFIGURATION_H
#define ALIYAMLCONFIGURATION_H

// CINT can't handle the yaml header because it has c++11!
#if !(defined(__CINT__) || defined(__MAKECINT__))
#include <yaml-cpp/yaml.h>
#endif

#include <string>
#include <vector>
#include <ostream>

#include <TObject.h>
#include <TString.h>

#include <AliLog.h>

/**
 * @class PWG::Tools::AliYAMLConfiguration
 * @brief %YAML configuration class for AliPhysics.
 *
 * A class to handle generic reading and writing to %YAML files. This can be used to configure tasks, coordinate trains,
 * or many other tasks which require configuration. While yaml-cpp can be used directly, this class handles many details
 * such as accessing files on AliEn, as well as generally simplifying the user experience. The class can also handle
 * multiple configuration files, first looking in the last added file, and then if the value is not found, looking in
 * previous configurations until the value is found. Values are accessed in this order so the values in the initial file
 * can be overridden by values defined in subsequently added files. Values that are used in multiple places can be set
 * together using "Shared Paramaters" (see the section below) or using %YAML anchors.
 *
 * Usage information:
 *
 * Consider the following example %YAML configuration.
 *
 * ~~~
 * hello:
 *   world:
 *     exampleValue: 10
 * example: 12.2
 * importantValues:
 *   - entry1
 *   - entry2
 *   - entry3
 * ~~~
 *
 * To use the class, first add the config class to your task (be certain to stream it to store the
 * configuration with your class!):
 *
 * ~~~{.cxx}
 * PWG::Tools::AliYAMLConfiguration fYAMLConfig; ///< YAML configuration.
 *
 * /// Retrieve the YAML configurations for direct access for printing, etc.
 * PWG::Tools::AliYAMLConfiguration & GetYAMLConfiguration() { return fYAMLConfig; }
 * ~~~
 *
 * Then, to use it in your task:
 *
 * ~~~{.cxx}
 * // Will only be checked if a requested value is not in the second configuration.
 * fYAMLConfig.AddConfiguration(filename, name);
 * // Will be checked first. This is so it can override values in the first added configuration.
 * fYAMLConfig.AddConfiguration(filename2, name2);
 *
 * // Once all configuration is done and the YAML nodes will not be modified any more
 * // (for example, at the end of an AddTask), call Initialize() to lock in the configurations
 * // for streaming to the grid.
 * fYAMLConfig.Initialize();
 * ~~~
 *
 * %YAML objects cannot be streamed due to CINT limitations, so after the configuration class is
 * streamed, it must be re-initialized. This can be done in any function after streaming has been
 * completed, such as `UserCreateOutputObjects()`. It only needs to be performed once.
 *
 * ~~~{.cxx}
 * fYAMLConfig.Reinitialize();
 * ~~~
 *
 * To access or write a value, use the GetProperty(...) or WriteProperty(...) functions. To use
 * them, you must define an object of the desired type that you would like to read or write, and
 * then describe the path to the property. The path consists of the names of YAML nodes, separated
 * by a user specified delimiter (":" is the default).  For the example YAML above, to read
 * "exampleValue", the user would define an int to be set to the read value and the path would
 * be "hello:world:exampleValue". As a convenience, there is a helper function which simplifies
 * specifying the path. It will instead take a std::vector of strings to specify the path,
 * thereby also setting the proper delimiter. So the path would be
 * `{"hello", "world", "exampleValue"}`. Note that you may need to explicitly specify such an
 * initialization as `std::vector<std::string>` if you none of the arguments are already defined
 * as variables.
 *
 * Explicitly, this would look like:
 *
 * ~~~{.cxx}
 * int tempInt = 0;
 * // True specifies that the value must exist in the config
 * fYAMLConfig.GetProperty({"hello", "world", "examplePath"}, tempInt, true);
 * std::cout << tempInt << "\n"; // Returns "10"
 * // Alternatively, you could specify the path explicitly (although not recommended)
 * fYAMLConfig.GetProperty("hello:world:exampleValue", tempInt);
 * // The type specifies what is returned
 * double tempDouble = 0.;
 * fYAMLConfig.GetProperty("example", tempDouble);
 * std::vector <std::string> importantValues;
 * fYAMLConfig.GetProperty("importantValues", importantValues);
 * // The vector "importantValues" now contains three strings
 * ~~~
 *
 * If you are accessing a %YAML file on AliEn (which is supported by default), the file will be
 * copied to the directory where the task is run. Conveniently, this means that the configuration
 * file is also easily accessible in the LEGO test train directory.
 *
 * That's basically all there is to using it. For more information, look at the documentation
 * of the various GetProperty(...) and WriteProperty(...) functions. The interfaces to both
 * functions are fairly similar.
 *
 * Notes on using shared parameters:
 *
 * sharedParameters are inherently limited. It can only retrieve values where the requested
 * type is arithmetic, string, or bool. The retrieved shared parameters value can only be of
 * those same types. Note that the shared parameters correspond to each configuration file.
 * ie. If `sharedParameters:test` is requested in the first configuration file, then it will
 * only look for the sharedParameters value in that configuration. Thus, if a sharedParameter
 * is requested in a later configuration file, the earlier configuration shared parameter
 * values will not be considered.
 *
 * Given the limitations, YAML anchors are recommended for more advanced usage as they can
 * be much more sophisticated.
 *
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University
 * @date Sept 19, 2017
 */

/**
 * Tell yaml-cpp how to encode and decode `TString`. Note that this is basically identical
 * (except for calling `.Data()`) to how std::string is handle in yaml-cpp/node/convert.h.
 * For more information, see: https://github.com/jbeder/yaml-cpp/wiki/Tutorial#converting-tofrom-native-data-types
 */
#if !(defined(__CINT__) || defined(__MAKECINT__))
namespace YAML {
template<>
  struct convert<TString> {
    static Node encode(const TString & str) {
      return Node(str.Data());
    }
    static bool decode(const Node & node, TString & str) {
      if (!node.IsScalar()) {
        return false;
      }
      str = node.Scalar();
      return true;
    }
  };
}
#endif

/**
 * Allow reverse range based iteration. We use this to iterate over the available
 * configurations in reverse to find the requested value. We iterate in reverse
 * because it allows configurations that are added later to override values in
 * configurations that are already added. Otherwise, overriding values could only
 * be performed by rearranging configuration files (which was not possible in the
 * EMCal Correction Framework, for exmaple).
 *
 * From: https://stackoverflow.com/a/21510202
 */
template <typename It>
class Range
{
  It b, e;

 public:
  Range(It b, It e) : b(b), e(e) {}
  It begin() const { return b; }
  It end() const { return e; }
};

template <typename ORange, typename OIt = decltype(std::begin(std::declval<ORange>())),
     typename It = std::reverse_iterator<OIt>>
Range<It> reverse(ORange&& originalRange)
{
  return Range<It>(It(std::end(originalRange)), It(std::begin(originalRange)));
}

// operator<< has to be forward declared carefully to stay in the global namespace so that it works with CINT.
// For generally how to keep the operator in the global namespace, See: https://stackoverflow.com/a/38801633
namespace PWG { namespace Tools { class AliYAMLConfiguration; } }
std::ostream & operator<< (std::ostream &in, const PWG::Tools::AliYAMLConfiguration &myTask);

namespace PWG {
namespace Tools {

class AliYAMLConfiguration : public TObject {
 public:
  AliYAMLConfiguration(const std::string prefixString = "AliEmcalCorrection", const std::string delimiterCharacter = ":");
  virtual ~AliYAMLConfiguration() {}

  /** @{
   * @name Setup the configurations at various points of the analysis.
   */
  bool Initialize();
  bool Reinitialize();
  /** @} */

  /** @{
   * @name Add a new particular configuration under a given name.
   */
  /// Add YAML configuration at configurationFilename to available configurations
  int AddEmptyConfiguration(const std::string & configurationName);
  int AddConfiguration(std::string configurationFilename, std::string configurationName = "");
  #if !(defined(__CINT__) || defined(__MAKECINT__))
  int AddConfiguration(const YAML::Node node, std::string configurationName = "");
  /** @} */

  /** @{
   * @name Get a particular configuration by index or name.
   */
  bool DoesConfigurationExist(const int i)                                              const { return i >= 0 && static_cast<const unsigned int>(i) < fConfigurations.size(); }
  bool DoesConfigurationExist(const std::string & name)                                 const { return DoesConfigurationExist(GetConfigurationIndexFromName(name, fConfigurations)); }
  const std::pair<std::string, YAML::Node> & GetConfiguration(const int i)              const { return fConfigurations.at(i); }
  const std::pair<std::string, YAML::Node> & GetConfiguration(const std::string & name) const { return GetConfiguration(GetConfigurationIndexFromName(name, fConfigurations)); }
  std::pair<std::string, YAML::Node> & GetConfiguration(const int i)                          { return fConfigurations.at(i); }
  std::pair<std::string, YAML::Node> & GetConfiguration(const std::string & name)             { return GetConfiguration(GetConfigurationIndexFromName(name, fConfigurations)); }
  /** @} */

  /** @{
   * @name Translate between configuration name and index.
   */
  std::string GetConfigurationNameFromIndex(const unsigned int i)                       const { return fConfigurations.at(i).first; }
  int GetConfigurationIndexFromName(const std::string & name)                           const { return GetConfigurationIndexFromName(name, fConfigurations); }
  /** @} */

  /** @{
   * @name Remove a particular configuration by index or name.
   */
  bool RemoveConfiguration(const unsigned int i);
  bool RemoveConfiguration(const std::string & name)                                          { return RemoveConfiguration(GetConfigurationIndexFromName(name, fConfigurations)); }
  /** @} */

  /** @{
   * @name Get a property from the available configurations
   */
  // Helper functions to retrieve property values
  template<typename T>
  bool GetProperty(std::vector<std::string> propertyPath, const std::string & propertyName, T & property, const bool requiredProperty) const;
  template<typename T>
  bool GetProperty(const std::vector<std::string> propertyPath, T & property, const bool requiredProperty) const;
  // Retrieve property driver function.
  template<typename T>
  bool GetProperty(std::string propertyName, T & property, const bool requiredProperty = true) const;
  /** @} */

  /** @{
   * @name Write a property to a particular configuration
   */
  // Write property driver function.
  template<typename T>
  bool WriteProperty(std::string propertyName, T & property, std::string configurationName = "");
  /** @} */
  #endif

  /** @{
   * @name Write a particular configuration at a given index or name to a file.
   */
  bool WriteConfiguration(const std::string & filename, const unsigned int i) const;
  bool WriteConfiguration(const std::string & filename, const std::string & configurationName) const;
  /** @} */

  /** @{
   * @name Write a particular configuration at a given index or name to a file.
   */
  bool CompareConfigurations(const int config1, const int config2)                     const;
  bool CompareConfigurations(const int config1, const std::string & config2)           const { return CompareConfigurations(config1, GetConfigurationIndexFromName(config2, fConfigurations)); }
  bool CompareConfigurations(const std::string & config1, const std::string & config2) const { return CompareConfigurations(GetConfigurationIndexFromName(config1, fConfigurations), GetConfigurationIndexFromName(config2, fConfigurations)); }

  /** @} */

  /** @{
   * @name Print a particular configuration at a given index or name to a file.
   */
  std::string toString(const int index = -1) const;
  std::ostream & Print(std::ostream &in, const int index = -1) const;
  std::ostream & Print(std::ostream &in, const std::string & configurationName) const;
  friend std::ostream & ::operator<< (std::ostream &in, const AliYAMLConfiguration &myTask);
  void Print(Option_t* /* opt */ = "") const;
  /** @} */

 protected:

  // Utility functions
  // File utilities
  inline bool DoesFileExist(const std::string & filename) const;
  void SetupReadingConfigurationFilePath(std::string & filename, const std::string & fileIdentifier) const;
  void WriteConfigurationToFilePath(const std::string & localFilename, std::string filename) const;
  #if !(defined(__CINT__) || defined(__MAKECINT__))
  // Printing
  void PrintConfiguration(std::ostream & stream, const std::pair<std::string, YAML::Node> & configPair) const;
  #endif
  // Configuration utilities
  template<typename T>
  unsigned int GetConfigurationIndexFromName(const std::string & name, const std::vector<std::pair<std::string, T>> & configurations) const;

  bool IsSharedValue(std::string & value) const;
  #if !(defined(__CINT__) || defined(__MAKECINT__))
  template<typename T>
  auto PrintRetrievedPropertyValueImpl(std::stringstream & tempMessage, const T & property, int) const -> decltype(tempMessage << property, void());
  template<typename T>
  auto PrintRetrievedPropertyValueImpl(std::stringstream & tempMessage, const std::vector <T> & property, int) const -> decltype(property.begin(), property.end(), tempMessage << std::declval<T>(), void());
  template<typename T>
  void PrintRetrievedPropertyValueImpl(std::stringstream & tempMessage, const T & property, long) const;
  template<typename T>
  void PrintRetrievedPropertyValue(std::stringstream & tempMessage, const T & property) const;

  template<typename T>
  bool GetPropertyFromNode(const YAML::Node & node, std::string propertyName, T & property) const;
  template<typename T>
  bool GetProperty(YAML::Node & node, YAML::Node & sharedParametersNode, const std::string & configurationName, std::string propertyName, T & property) const;

  template<typename T>
  void WriteValue(YAML::Node & node, std::string propertyName, T & proeprty);

  std::vector<std::pair<std::string, YAML::Node> > fConfigurations;         //!<! Contains all YAML configurations. The first element has the highest precedence.
  #endif
  std::vector<std::pair<std::string, std::string> > fConfigurationsStrings; ///<  Contains all YAML configurations as strings so that they can be streamed.

  bool fInitialized;                          ///< True if the configurations have been initialized.
  std::string fPrefixString;                  ///< Contains the prefix of any names base node names which should be removed.
  std::string fDelimiter;                     ///< Delimiter character to separate each level of the request.

  /// \cond CLASSIMP
  ClassDef(AliYAMLConfiguration, 1); // YAML Configuration
  /// \endcond
};

#if !(defined(__CINT__) || defined(__MAKECINT__))
/**
 * Prints the retrieved property value if the type implements operator<<(). For more on how this achieved,
 * see: https://stackoverflow.com/a/22759368 and https://stackoverflow.com/a/11866675.
 * NOTE: This function should not be called directly, but through PrintRetrievedPropertyValue(...)
 *
 * NOTE: This function returns void!
 *
 * @param[in] property Property to be printed
 * @param[in, out] tempMessage Stringstream into which the property should be streamed
 */
template<typename T>
auto AliYAMLConfiguration::PrintRetrievedPropertyValueImpl(std::stringstream & tempMessage, const T & property, int) const -> decltype(tempMessage << property, void())
{
  tempMessage << " with value \"" << property << "\"";
}

/**
 * Prints the retrieved values in a vector if the value is a vector (more specifically, if it implements begin() and end())
 * and if the type of the vector can be streamed in an ostream object.
 * NOTE: This function should not be called directly, but through PrintRetrievedPropertyValue(...)
 *
 * @param[in] property Property to be printed
 * @param[in, out] tempMessage Stringstream into which the property should be streamed
 */
template<typename T>
auto AliYAMLConfiguration::PrintRetrievedPropertyValueImpl(std::stringstream & tempMessage, const std::vector <T> & property, int) const -> decltype(property.begin(), property.end(), tempMessage << std::declval<T>(), void())
{
  tempMessage << " with value(s):";
  for (auto it = property.begin(); it != property.end(); it++) {
    tempMessage << "\n\t- " << *it;
  }
}

/**
 * Handles all other cases where the property is not covered by the other functions.
 *
 * NOTE: This function should not be called directly, but through PrintRetrievedPropertyValue(...)
 *
 * NOTE: Cannot use "..." as a function argument here as a fall back because it causes
 * problems with ROOT dictionary generation...
 *
 * @param[in] property Property to be printed
 * @param[in, out] tempMessage Stringstream into which the property should be streamed
 */
template<typename T>
void AliYAMLConfiguration::PrintRetrievedPropertyValueImpl(std::stringstream & tempMessage, const T & property, long) const
{
  // Cannot easily print these types, so just note that is the case!
  tempMessage << " with a value that cannot be trivially printed";
}

/**
 * Wrapper function to resolve ambiguity in function overloading. By passing 0 as the last argument to the
 * PrintRetrievedPropertyValueImpl(...) functions, it will first attempt to use the matching int function,
 * and then fall back the to the long function. This is effectively a hack to help c++ properly resolve the
 * function to call.
 *
 * Inspired by: https://stackoverflow.com/a/38283990
 *
 * @param[in] property Property to be printed
 * @param[in, out] tempMessage Stringstream into which the property should be streamed
 */
template<typename T>
void AliYAMLConfiguration::PrintRetrievedPropertyValue(std::stringstream & tempMessage, const T & property) const
{
  // By passing zero, it will first lookup the int, and then the long
  AliYAMLConfiguration::PrintRetrievedPropertyValueImpl(tempMessage, property, 0);
}

/**
 * Utility function to retrieve the property from a already selected YAML node.
 *
 * @param[in] node YAML node from which the property should be retrieved
 * @param[in] propertyName Name of the property to retrieve
 * @param[out] property Contains the retrieved property
 *
 * @return True if the property was set successfully
 */
template<typename T>
bool AliYAMLConfiguration::GetPropertyFromNode(const YAML::Node & node, std::string propertyName, T & property) const
{
  if (node[propertyName])
  {
    property = node[propertyName].as<T>();
    return true;
  }
  return false;
}

/**
 * Helper function for main retrieval function. It automatically adds the property name to the end of the property path.
 * By doing so, it allows the user to specify a general property path and then vary the property name between function calls.
 *
 * @param[in] propertyPath Path to the property in the YAML file.
 * @param[in] propertyName Name of the property to be retrieved.
 * @param[out] property Contains the retrieved property.
 * @param[in] requiredProperty True if the property is required
 *
 * @return True if the property was set successfully
 */
template<typename T>
bool AliYAMLConfiguration::GetProperty(std::vector <std::string> propertyPath, const std::string & propertyName, T & property, const bool requiredProperty) const
{
  propertyPath.push_back(propertyName);
  return GetProperty(propertyPath, property, requiredProperty);
}

/**
 * Helper function for main retrieval function. Each value in the property path will be joined together, separated by the delimiter
 * set in the configuration class.
 *
 * @param[in] propertyPath Path to the property in the YAML file.
 * @param[out] property Contains the retrieved property.
 * @param[in] requiredProperty True if the property is required
 *
 * @return True if the property was set successfully
 */
template<typename T>
bool AliYAMLConfiguration::GetProperty(const std::vector <std::string> propertyPath, T & property, const bool requiredProperty) const
{
  // Combine the requested names together
  std::string requestedName = "";
  for (auto & str : propertyPath)
  {
    if (requestedName.length() > 0) {
      requestedName += ":" + str;
    }
    else {
      requestedName = str;
    }
  }

  // Pass on the properly requested call
  return GetProperty(requestedName, property, requiredProperty);
}

/**
 * Main general function to get a property from a set of YAML configuration files. It first calls checks
 * the user configuration, and then if the property is not found, it then checks the default configuration.
 * If the property is required but it is not found, a fatal error is thrown.
 *
 * @param[in] propertyName Name of the property to retrieve
 * @param[out] property Contains the retrieved property
 * @param[in] requiredProperty True if the property is required
 *
 * @return True if the property was set successfully
 */
template<typename T>
bool AliYAMLConfiguration::GetProperty(std::string propertyName, T & property, const bool requiredProperty) const
{
  // Remove the prefix string if in name
  std::size_t prefixStringLocation = propertyName.find(fPrefixString);
  if (prefixStringLocation != std::string::npos)
  {
    // Remove the prefix string
    propertyName.erase(prefixStringLocation, prefixStringLocation + fPrefixString.length());
  }

  bool setProperty = false;
  // Search in reverse so it is possible to override configuration values.
  for (auto configPair : reverse(fConfigurations))
  {
    if (setProperty == true) {
      AliDebugGeneralStream("AliYAMLConfiguration", 2) << "Property \"" << propertyName << "\" found!\n";
      break;
    }

    // IsNull checks is a node is empty. A node is empty if it is created.
    // IsDefined checks if the node that was requested was not actually created.
    if (configPair.second.IsNull() != true)
    {
      AliDebugGeneralStream("AliYAMLConfiguration", 2) << "Looking for parameter \"" << propertyName << "\" in \"" << configPair.first << "\" configuration\n";
      // NOTE: This may not exist, but that is entirely fine.
      YAML::Node sharedParameters = configPair.second["sharedParameters"];
      setProperty = GetProperty(configPair.second, sharedParameters, configPair.first, propertyName, property);
    }
  }

  if (setProperty != true && requiredProperty == true)
  {
    std::stringstream message;
    message << "Failed to retrieve required property \""
        << propertyName << "\" from available configurations!" << std::endl;
    AliFatalGeneral("AliYAMLConfiguration", message.str().c_str());
  }

  // Return whether the value was actually set
  return setProperty;
}

/**
 * Actually handles retrieving parameters from YAML nodes.
 *
 * For propertyName = "hello:example_special:property", the retrieval procedure is as follows:
 *  - Extract node name by searching the propertyName up to the delimiter (in this case, "hello").
 *  - Look up the node name.
 *    - If the node exists, recurse with node[nodeName] and with remaining propertyName = "example_special:property".
 *    - If it doesn't exist, check for the specialization delimiter (default: "_") and take the node name as everything before the delimiter
 *      (for example, "example_special" leads to "example").
 *  - If the delimiter is not found, then atempt to extract the property with the remaining node name.
 *    - If the property is arithmetic, std::string, or bool, it will check if the property is instead a shared parameter. If so, it will retrieve
 *      the shared parameter value.
 *    - The property is stored in the property, and true is returned if the effort was successful.
 *
 * @param[in] node Main YAML node containing the properties
 * @param[in] sharedParametersNode YAML node containing the shared parameters
 * @param[in] configurationName Name of the configuration type.
 * @param[in] propertyName Name of the property to retrieve
 * @param[out] property Contains the retrieved property
 *
 * @return True if the property was set successfully
 */
template<typename T>
bool AliYAMLConfiguration::GetProperty(YAML::Node & node, YAML::Node & sharedParametersNode, const std::string & configurationName, std::string propertyName, T & property) const
{
  // Used as a buffer for printing complicated messages
  std::stringstream tempMessage;

  bool returnValue = false;

  const std::string specializationDelimiter = "_";
  size_t delimiterPosition = 0;
  std::string tempPropertyName = propertyName;

  if ((delimiterPosition = tempPropertyName.find(fDelimiter)) != std::string::npos)
  {
    std::string nodeName = tempPropertyName.substr(0, delimiterPosition);
    tempPropertyName.erase(0, delimiterPosition + fDelimiter.length());

    // Attempt to use the node name
    if (node[nodeName].IsDefined() == true)
    {
      // Retrieve node and then recurse
      YAML::Node tempNode = node[nodeName];
      AliDebugGeneralStream("AliYAMLConfiguration", 2) << "Attempting to retrieving property \"" << tempPropertyName << "\" by going a node deeper with node \"" << nodeName << "\".\n";
      returnValue = GetProperty(tempNode, sharedParametersNode, configurationName, tempPropertyName, property);
    }

    // Check for the specialization if the nodeName is undefined.
    // Alternatively, if the value was not returned successfully, we should also check for the specialization
    //   such as inheritnace for input objects.
    if (node[nodeName].IsDefined() == false || returnValue == false)
    {
      // Check for specialization
      if ((delimiterPosition = nodeName.find(specializationDelimiter)) != std::string::npos)
      {
        std::string specializationNodeName = nodeName.substr(0, delimiterPosition);
        YAML::Node tempNode = node[specializationNodeName];
        AliDebugGeneralStream("AliYAMLConfiguration", 2) << "Attempting to retrieving property \"" << tempPropertyName << "\" by going a node deeper through dropping the specializtion and using node \"" << specializationNodeName << "\".\n";
        returnValue = GetProperty(tempNode, sharedParametersNode, configurationName, tempPropertyName, property);
      }
      else {
        returnValue = false;
      }
    }
  }
  else
  {
    // Check if node exists!
    if (node[propertyName])
    {
      // Handle shared parameters
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
        isShared = IsSharedValue(sharedValueName);
      }

      // Setup printable message
      tempMessage.str("");
      tempMessage << "property \""
            << propertyName
            << "\" using " << ( isShared ? "\"sharedParameters:" + sharedValueName + "\" in " : "" )
            << "values from the " << configurationName << " configuration";

      AliDebugGeneralStream("AliYAMLConfiguration", 2) << "Retrieveing " << tempMessage.str() << "\n";

      // Retrieve result
      bool retrievalResult = false;
      if (isShared == true) {
        // Retrieve from the shared parameter node directly.
        retrievalResult = GetPropertyFromNode(sharedParametersNode, sharedValueName, property);
      }
      else {
        retrievalResult = GetPropertyFromNode(node, propertyName, property);
      }

      // Inform about the result
      if (retrievalResult == true) {
        // Add the retrieved value to the message (only if trivially printable)
        AliYAMLConfiguration::PrintRetrievedPropertyValue(tempMessage, property);
        AliDebugGeneralStream("AliYAMLConfiguration", 2) << "Succeeded in retrieveing " << tempMessage.str() << "\n";
        returnValue = true;
      }
      else {
        returnValue = false;
      }
    }
  }

  return returnValue;
}

/**
 * Write a value to a YAML configuration. Note that the value is written to the YAML configuration, but then
 * the YAML configuration needs to explicitly be written to file if the changes should be persistent.
 *
 * @param[in] propertyName Path to the property in the YAML file.
 * @param[in] property Property to be written.
 * @param[in] configurationName Name of the YAML configuration file to which the property should be written.
 *
 * @return True if the write was successful.
 */
template<typename T>
bool AliYAMLConfiguration::WriteProperty(std::string propertyName, T & property, std::string configurationName)
{
  unsigned int configurationIndex = 0;
  if (configurationName != "")
  {
    configurationIndex = GetConfigurationIndexFromName(configurationName, fConfigurations);
  }

  if (fConfigurations.size() == 0) {
    AliErrorStream() << "No configurations available! Property will not be written!\n";
    return false;
  }

  std::pair<std::string, YAML::Node> & configPair = fConfigurations.at(configurationIndex);

  WriteValue(configPair.second, propertyName, property);
  AliDebugGeneralStream("AliYAMLConfiguration", 2) << "Final Node:\n" << configPair.second << "\n";

  return true;
}

/**
 * Determine where the value should be written and write it to the appropriate YAML node.
 *
 * @param[in, out] node YAML node where the property should (eventually) be written.
 * @param[in] propertyName Name of the property to be written.
 * @param[in] property Property to be written to the YAML node.
 */
template<typename T>
void AliYAMLConfiguration::WriteValue(YAML::Node & node, std::string propertyName, T & property)
{
  // Find the node name we are interested in
  // Derived from: https://stackoverflow.com/a/14266139
  const std::string delimiter = ":";
  size_t delimiterPosition = 0;
  std::string tempPropertyName = propertyName;

  if ((delimiterPosition = tempPropertyName.find(delimiter)) != std::string::npos)
  {
    std::string nodeName = tempPropertyName.substr(0, delimiterPosition);
    tempPropertyName.erase(0, delimiterPosition + delimiter.length());
    AliDebugGeneralStream("AliYAMLConfiguration", 2) << "nodeName: " << nodeName << "\n";
    AliDebugGeneralStream("AliYAMLConfiguration", 2) << "Node Before:\n" << node << "\n";
    if (node[nodeName].IsDefined()) {
      AliDebugGeneralStream("AliYAMLConfiguration", 2) << "Using existing node\n";
      YAML::Node tempNode = node[nodeName];
      WriteValue(tempNode, tempPropertyName, property);
    }
    else {
      AliDebugGeneralStream("AliYAMLConfiguration", 2) << "Creating new node\n";
      YAML::Node tempNode;
      node[nodeName] = tempNode;
      WriteValue(tempNode, tempPropertyName, property);
    }
    AliDebugGeneralStream("AliYAMLConfiguration", 2) << "Node After:\n" << node << "\n";
  }
  else {
    // We have finished the recursion. Now we need to save the actual property.
    node[propertyName] = property;
  }
}

#endif

/**
 * Get the index of the configuration given the configuration name. If the same name is used multiple times, then
 * the first instance will be returned.
 *
 * @param[in] name Name of the configuration.
 * @param[in] configurations Configurations to search through (could contain YAML nodes, or the string copies).
 *
 * @return Index of the configuration, or -1 if not found.
 */
template<typename T>
unsigned int AliYAMLConfiguration::GetConfigurationIndexFromName(const std::string & name, const std::vector<std::pair<std::string, T>> & configurations) const
{
  int index = -1;
  for (const auto & configPair : configurations)
  {
    if (configPair.first == name) {
      // Retrieve index
      // See: https://stackoverflow.com/a/10962459
      index = std::addressof(configPair) - std::addressof(configurations[0]);
      break;
    }
  }

  return index;
}

} // namespace Tools
} // namespace PWG

#endif /* ALIYAMLCONFIGURATION_H */
