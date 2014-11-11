/*
 * AliEMCALJSONReader.h
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */

#ifndef PWG_EMCAL_ALIEMCALJSONREADER_H_
#define PWG_EMCAL_ALIEMCALJSONREADER_H_

#include "AliEMCALConfigurationObject.h"
#include <string>
#include <vector>

class TList;

class AliEMCALJSONSyntaxTreeNode{
public:
  AliEMCALJSONSyntaxTreeNode(const char *name, AliEMCALJSONSyntaxTreeNode *mother);
  ~AliEMCALJSONSyntaxTreeNode();

  void AddEntry(AliEMCALConfigurationObject *entry);
  AliEMCALJSONSyntaxTreeNode *CreateDaughter(const char *name);
  void SetOwner(bool owner = true);
  AliEMCALJSONSyntaxTreeNode *GetMotherNode() const { return fMotherNode; }
  std::vector<AliEMCALJSONSyntaxTreeNode *> &GetDaughters() { return fDaughters; };
  std::vector<AliEMCALConfigurationObject *> &GetEntries() { return fEntries; }
  const char *GetName() const { return fName.c_str(); }

private:
  AliEMCALJSONSyntaxTreeNode(const AliEMCALJSONSyntaxTreeNode &ref);
  AliEMCALJSONSyntaxTreeNode *operator=(const AliEMCALJSONSyntaxTreeNode &ref);

  std::string fName;
  AliEMCALJSONSyntaxTreeNode *fMotherNode;
  std::vector<AliEMCALConfigurationObject *> fEntries;
  std::vector<AliEMCALJSONSyntaxTreeNode *> fDaughters;
  bool fOwner;
};

class AliEMCALJSONReader {
public:
  AliEMCALJSONReader();
  virtual ~AliEMCALJSONReader();

  TList *Decode(const char *jsosnstring) const;

private:
  void AddNodeToList(AliEMCALJSONSyntaxTreeNode *node, TList *consumer) const;
};

#endif /* PWG_EMCAL_ALIEMCALJSONREADER_H_ */
