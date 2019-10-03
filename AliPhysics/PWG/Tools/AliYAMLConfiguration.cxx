//
//

#include "AliYAMLConfiguration.h"

#include <cstdio>
#include <fstream>

#include <TSystem.h>
#include <TGrid.h>
#include <TFile.h>
#include <TUUID.h>

/// \cond CLASSIMP
ClassImp(PWG::Tools::AliYAMLConfiguration);
/// \endcond

namespace PWG {
namespace Tools {

/**
 * Default constructor.
 *
 * @param[in] prefixString Prefix to remove when gettnig a property.
 * @param[in] delimiterCharacter Character that delimits between each part of the YAML path.
 */
AliYAMLConfiguration::AliYAMLConfiguration(const std::string prefixString, const std::string delimiterCharacter):
  TObject(),
  fConfigurations(),
  fConfigurationsStrings(),
  fInitialized(false),
  fPrefixString(prefixString),
  fDelimiter(delimiterCharacter)
{
}

/**
 * Create an empty YAML node and adds it to the available configurations.
 *
 * @param[in] configurationName Name of the YAML node. The node will be stored under this name.
 *
 * @return Position of configuration if the configuration was added successfully. -1 if unsuccessful.
 */
int AliYAMLConfiguration::AddEmptyConfiguration(const std::string & configurationName)
{
  YAML::Node node;
  AliInfoStream() << "Adding configuration \"" << configurationName << "\" as an empty YAML node.\n";
  return AddConfiguration(node, configurationName);
}

/**
 * Add YAML configuration at a given filename to the available configurations.
 *
 * @param[in] configurationFilename Filename of the YAML configuration file to be added.
 * @param[in] configurationName Name of the YAML node. The node will be stored under this name.
 *
 * @return Position of configuration if the configuration was added successfully. -1 if unsuccessful.
 */
int AliYAMLConfiguration::AddConfiguration(std::string configurationFilename, std::string configurationName)
{
  SetupReadingConfigurationFilePath(configurationFilename, configurationName);

  // Add file
  if (DoesFileExist(configurationFilename) == false) {
    AliErrorStream() << "Configuration filename \"" << configurationFilename << "\" does not exist!\n";
    return -1;
  }

  // Create node from filename
  auto node = YAML::LoadFile(configurationFilename);

  if (node.IsNull() == true) {
    AliErrorStream() << "The node at configuration filename \"" << configurationFilename << "\" is null and will not be added!\n";
    return -1;
  }

  AliInfoStream() << "Adding configuration \"" << configurationName << "\" located at \"" << configurationFilename << "\".\n";
  return AddConfiguration(node, configurationName);
}

/**
 * Add a YAML node to the available configurations.
 *
 * @param[in] node YAML node to be added.
 * @param[in] configurationName Name of the YAML node. The node will be stored under this name.
 *
 * @return Position of configuration if the configuration was added successfully. -1 if unsuccessful.
 */
int AliYAMLConfiguration::AddConfiguration(const YAML::Node node, std::string configurationName)
{
  if (configurationName == "")
  {
    // Attempt to get name from node
    // It is fine if the node is empty
    GetPropertyFromNode(node, "name", configurationName);
  }

  if (configurationName == "") {
    AliErrorStream() << "Could not determine a name for the configuration:\n\"\n" << node << "\n\" Configuration will not be added!\n";
    return -1;
  }

  // Add the configuration
  AliDebugStream(2) << "Adding configuration \"" << configurationName << "\".\n";
  fConfigurations.push_back(std::make_pair(configurationName, node));

  // Return the location of the new configuration
  return fConfigurations.size() - 1;
}

/**
 * Remove a configuration at a given index.
 *
 * @param[in] i Index of the configuration to be rmeoved.
 *
 * @return true if the configuration was removed.
 */
bool AliYAMLConfiguration::RemoveConfiguration(const unsigned int i)
{
  bool returnValue = false;
  if (i < fConfigurations.size())
  {
    fConfigurations.erase(fConfigurations.begin() + i);
    returnValue = true;
  }

  return returnValue;
}

/**
 * Write a YAML configuration node to file.
 *
 * @param[in] index Index of the YAML configuration node.
 * @param[in] filename Filename to write the node to.
 *
 * @return True if write was successful.
 */
bool AliYAMLConfiguration::WriteConfiguration(const std::string & filename, const unsigned int index) const
{
  // Write to a local temp filename
  TUUID tempUUID;
  std::string localFilename = tempUUID.AsString();
  localFilename += ".yaml";

  if (!DoesConfigurationExist(index)) {
    AliWarningStream() << "Requested configuration index " << index << " does not exist - it cannot be written!\n";
    return false;
  }
  auto configPair = fConfigurations.at(index);
  std::ofstream outputFile(localFilename);
  outputFile << configPair.second;
  outputFile.close();

  // Hanldes writing to places like AliEn
  WriteConfigurationToFilePath(localFilename, filename);

  remove(localFilename.c_str());

  return true;
}

/**
 * Write a YAML configuration node to file.
 *
 * @param[in] configurationName Name of the YAML configuration node.
 * @param[in] filename Filename to write the node to.
 *
 * @return True if write was successful.
 */
bool AliYAMLConfiguration::WriteConfiguration(const std::string & filename, const std::string & configurationName) const
{
  return WriteConfiguration(filename, GetConfigurationIndexFromName(configurationName, fConfigurations));
}

/**
 * Compare two configurations to see if they are identical.
 *
 * @param[in] configIndex1 Index of the first configuration.
 * @param[in] configIndex2 Index of the second configuration.
 *
 * @return True if the configurations are the same.
 */
bool AliYAMLConfiguration::CompareConfigurations(const int configIndex1, const int configIndex2) const
{
  bool returnValue = false;

  bool configsExist = true;
  std::vector<int> configs = {configIndex1, configIndex2};
  for (auto config : configs) {
    if (!DoesConfigurationExist(config)) {
      AliErrorStream() << "Configuration at index " << config << " does not exist.\n";
      configsExist = false;
    }
  }

  if (configsExist)
  {
    // Generate YAML nodes for the comparison
    auto configPair1 = GetConfiguration(configIndex1);
    auto configPair2 = GetConfiguration(configIndex2);

    // Need to stream the configuration back to a string to remove the comments
    // since they are not preserved in the YAML node.
    std::stringstream config1SS;
    config1SS << configPair1.second;
    std::stringstream config2SS;
    config2SS << configPair2.second;

    // Compare the nodes. Make the comparison as strings, as the YAML nodes do _not_ match,
    // despite the strings matching. In fact, the YAML nodes will _not_ match even if they
    // are generated from the same string....
    if (config1SS.str() == config2SS.str()) {
      returnValue = true;
    }
    else {
      // Already should be the case, but just to be explicit.
      returnValue = false;

      // Inform the user about the details of the mismatch
      std::stringstream errorMessageSS;
      errorMessageSS << "Configuration mismatch between configuration at index " << configIndex1 << " and at index " << configIndex2 << "\n";
      errorMessageSS << "Config 1:\n";
      Print(errorMessageSS, configIndex1);
      errorMessageSS << "Config 2:\n";
      Print(errorMessageSS, configIndex2);
      AliWarningStream() << errorMessageSS.str();
    }
  }

  return returnValue;
}


/**
 * Checks if a file exists. This is done inline to make it efficient.
 * See: https://stackoverflow.com/a/19841704
 *
 * @param filename String containing the filename of the file to check.
 *
 * @return True if the file exists.
 */
inline bool AliYAMLConfiguration::DoesFileExist(const std::string & filename) const
{
  std::ifstream inFile(filename);
  return inFile.good();
}

/**
 * Handles setting up the configuration file to be opened, including in AliPhysics and on the grid.
 * Cannot just use TFile::Open() because the YAML file is just text as opposed to a root file.
 * In the case of a file on the grid, it is copied locally.
 *
 * @param[in,out] filename Name of the file to be open
 * @param[in] fileIdentifier Additional file identifier to add onto the file name
 */
void AliYAMLConfiguration::SetupReadingConfigurationFilePath(std::string & filename, const std::string & fileIdentifier) const
{
  if (filename != "")
  {
    // Handle if in AliPhysics and includes $ALICE_PHYSICS
    filename = gSystem->ExpandPathName(filename.c_str());

    // Handle grid
    if(filename.find("alien://") != std::string::npos)
    {
      AliDebug(2, TString::Format("Opening file \"%s\" on the grid!", filename.c_str()));
      // Initialize AliEn connection if needed
      if (!gGrid) {
        TGrid::Connect("alien://");
      }

      // Determine the local filename and copy file to local directory
      std::string localFilename = gSystem->BaseName(filename.c_str());
      // Add identifier if it's not an empty string
      if (fileIdentifier != "") {
        localFilename = fileIdentifier + "." + localFilename;
      }
      // Add UUID to ensure there are no conflicts if multiple yaml configs have the same configuration file name
      TUUID tempUUID;
      localFilename = "." + localFilename;
      localFilename = tempUUID.AsString() + localFilename;

      // Copy file
      TFile::Cp(filename.c_str(), localFilename.c_str());

      // yaml-cpp should only open the local file
      filename = localFilename;
    }
  }
}

/**
 * Write a selected YAML configuration to file. Practically, it copies a local file to the desired location]
 * to ensure seamless access to AliEn.
 *
 * @param[in] filename Filename to which the configuration should be written.
 * @param[in] localFilename Filename where the configuration was written locally.
 */
void AliYAMLConfiguration::WriteConfigurationToFilePath(const std::string & localFilename, std::string filename) const
{
  bool cannotWriteFile = false;
  if (localFilename == "") {
    AliErrorStream() << "Local filename is null, so the file cannot be written!\n";
    cannotWriteFile = true;
  }
  if (filename == "") {
    AliErrorStream() << "Filename is null, so the file cannot be written\n";
    cannotWriteFile = true;
  }

  if (cannotWriteFile == false) {
    // Handle if in AliPhysics and includes $ALICE_PHYSICS
    filename = gSystem->ExpandPathName(filename.c_str());

    // Handle grid
    if(filename.find("alien://") != std::string::npos)
    {
      AliDebugStream(2) << "Writing file \"" << filename << "\" on the grid!\n";
      // Initialize AliEn connection if needed
      if (!gGrid) {
        TGrid::Connect("alien://");
      }
    }

    // Copy file
    AliDebugStream(2) << "Copying localFilename \"" << localFilename << "\" to filename \"" << filename << "\".\n";
    TFile::Cp(localFilename.c_str(), filename.c_str());
  }
}

/**
 * Initialize configurations.
 *
 * This includes storing the contents of the YAML configurations
 * into strings so that they can be streamed to the grid.
 * (yaml-cpp objects do not work properly with ROOT streamers).
 *
 * NOTE: fConfigurationInitialized is set to true if the function is successful.
 */
bool AliYAMLConfiguration::Initialize()
{
  // Copy all configurations to their respective strings.
  // That way, the strings will be successfully streamed by ROOT and the YAML nodes can be re-initialized on the grid.
  std::stringstream tempSS;
  for (auto & configPair : fConfigurations)
  {
    tempSS.str("");
    tempSS << configPair.second;
    fConfigurationsStrings.push_back(std::make_pair(configPair.first, tempSS.str()));
  }

  fInitialized = true;

  return fInitialized;
}

/**
 * Reinitialize the configurations from the strings after streaming. This is required because the YAML
 * node objects cannot be streamed.
 *
 * This should be called immediately after the object is streamed on the grid. For example, at
 * UserCreateOutputObjects().
 *
 * @return True if the configurations were re-initialized.
 */
bool AliYAMLConfiguration::Reinitialize()
{
  bool returnValue = false;

  // If they were not streamed then the size of the vectors should not match
  if (fConfigurations.size() != fConfigurationsStrings.size())
  {
    if (fInitialized == false)
    {
      // Not initialized, so the strings do not represent the nodes
      AliFatalGeneral("AliYAMLConfiguration", "Attempted to re-initialize the YAML nodes, but the string based configurations are not available. Did you remember to call initialize?");
    }

    for (const auto & configStrPair : fConfigurationsStrings)
    {
      YAML::Node node = YAML::Load(configStrPair.second);
      fConfigurations.push_back(std::make_pair(configStrPair.first, node));
    }

    returnValue = true;
  }

  return returnValue;
}

/**
 * Check if value is a shared parameter, meaning we should look
 * at another node. Also edits the input string to remove "sharedParameters:"
 * if it exists, making it ready for use.
 *
 * @param[in] value String containing the string value return by the parameter.
 *
 * @return True if the value is shared.
 */
bool AliYAMLConfiguration::IsSharedValue(std::string & value) const
{
  std::string stringToFind = "sharedParameters:";
  std::size_t sharedParameterLocation = value.find(stringToFind);
  if (sharedParameterLocation != std::string::npos)
  {
    value.erase(sharedParameterLocation, sharedParameterLocation + stringToFind.length());
    return true;
  }
  // Return false otherwise
  return false;
}

/**
 * Print a particular configuration.
 *
 * @param[in, out] stream output stream.
 * @param[in] configPair Pair of string and YAML::Node containing the configuration to be printed.
 */
void AliYAMLConfiguration::PrintConfiguration(std::ostream & stream, const std::pair<std::string, YAML::Node> & configPair) const
{
  stream << "\nConfiguration Name: \"" << configPair.first << "\"\n\n";
  stream << configPair.second << "\n";
}

/**
 * Prints information about the YAML configuration(s).
 *
 * If index = -1 (which is the default argument), then all configurations will be printed!
 *
 * @param index Index of the YAML configuration to be printed.
 * @return std::string containing information about the task.
 */
std::string AliYAMLConfiguration::toString(const int index) const
{
  // Reinitialize the YAML configuration nodes for improved visualization. This function
  // only reinitializes from the stream strings if the YAML nodes aren't created.
  // NOTE: We cast away the const here because we want to keep the overload of Print().
  //       This of course usually not a good thing to do, but it's also the only reasonable option
  //       because ROOT 5 cannot handle YAML::Node at all. Also, this basically amounts to
  //       generating the YAML::Node(s) from an existing member, so the modification of the
  //       object isn't modified by much.
  bool reinitialized = false;
  if (fInitialized) {
    // Only attempt to run this if initialized (and therefore would need reinitialization).
    const_cast<AliYAMLConfiguration *>(this)->Reinitialize();
    reinitialized = true;
  }

  std::stringstream tempSS;

  if (index < 0) {
    // Print header information
    tempSS << std::boolalpha;
    tempSS << "AliYAMLConfiguration:\n";
    tempSS << "Initialized: " << fInitialized << "\n";
    // As a note to the user
    if (reinitialized) {
      tempSS << "Reinitialized YAML nodes during print: " << reinitialized << "\n";
    }
    tempSS << "Prefix string: \"" << fPrefixString << "\"\n";
    tempSS << "Delimiter: \"" << fDelimiter << "\"\n";
    tempSS << "Configurations vector length: " << fConfigurations.size() << ", Configurations string vector length: " << fConfigurationsStrings.size() << "\n";
    tempSS << "\nYAML Configurations:\n";

    // Print all configurations
    for (auto configPair : fConfigurations) {
      PrintConfiguration(tempSS, configPair);
    }
  }
  else {
    // Print a particular configuration
    if (static_cast<const unsigned int>(index) < fConfigurations.size()) {
      PrintConfiguration(tempSS, fConfigurations.at(index));
    }
    else {
      tempSS << "Index " << index << " is out of range!\n";
    }
  }

  return tempSS.str();
}

/**
 * Print YAML configuration information on an output stream using the string representation provided by
 * AliYAMLConfiguration::toString(). Used by operator<<.
 *
 * If index = -1 (which is the default argument), then all configurations will be printed!
 *
 * @param in output stream stream
 * @param index Index of the configuration to print.
 * @return reference to the output stream
 */
std::ostream & AliYAMLConfiguration::Print(std::ostream & in, const int index) const {
  in << toString(index);
  return in;
}

/**
 * Print YAML configuration information on an output stream using the string representation provided by
 * AliYAMLConfiguration::toString(). Used by operator<<
 *
 * @param in output stream stream.
 * @param configurationName Name of the configuration to print.
 * @return reference to the output stream
 */
std::ostream & AliYAMLConfiguration::Print(std::ostream & in, const std::string & configurationName) const {
  std::ostream & result = Print(in, GetConfigurationIndexFromName(configurationName, fConfigurations));
  return result;
}

/**
 * Print basic YAML configuration information using the string representation provided by
 * AliYAMLConfiguration::toString()
 */
void AliYAMLConfiguration::Print(Option_t* opt) const
{
  AliInfoStream() << toString();
}

} // namespace Tools
} // namespace PWG

/**
 * Implementation of the output stream operator for AliYAMLConfiguration. Printing
 * basic YAML configuration information provided by the function toString(). Note that
 * since this is friend, it is defined outside the namespace.
 *
 * @param in output stream
 * @param myTask Task which will be printed
 * @return Reference to the output stream
 */
std::ostream & operator<<(std::ostream & in, const PWG::Tools::AliYAMLConfiguration & myTask)
{
  std::ostream & result = myTask.Print(in);
  return result;
}

