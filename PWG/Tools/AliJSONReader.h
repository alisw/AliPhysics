/*
 * AliJSONReader.h
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */

#ifndef _ALIJSONREADER_H_
#define _ALIJSONREADER_H_

#include "AliJSONData.h"
#include <string>
#include <vector>

class TList;

class AliJSONSyntaxTreeNode{
public:
  AliJSONSyntaxTreeNode(const char *name, AliJSONSyntaxTreeNode *mother);
  ~AliJSONSyntaxTreeNode();

  void AddEntry(AliJSONData *entry);
  AliJSONSyntaxTreeNode *CreateDaughter(const char *name);
  void SetOwner(bool owner = true);
  AliJSONSyntaxTreeNode *GetMotherNode() const { return fMotherNode; }
  std::vector<AliJSONSyntaxTreeNode *> &GetDaughters() { return fDaughters; };
  std::vector<AliJSONData *> &GetEntries() { return fEntries; }
  const char *GetName() const { return fName.c_str(); }

private:
  AliJSONSyntaxTreeNode(const AliJSONSyntaxTreeNode &ref);
  AliJSONSyntaxTreeNode *operator=(const AliJSONSyntaxTreeNode &ref);

  std::string fName;
  AliJSONSyntaxTreeNode *fMotherNode;
  std::vector<AliJSONData *> fEntries;
  std::vector<AliJSONSyntaxTreeNode *> fDaughters;
  bool fOwner;
};

class AliJSONReader {
public:
  AliJSONReader();
  virtual ~AliJSONReader();

  TList *Decode(const char *jsosnstring) const;

private:
  void AddNodeToList(AliJSONSyntaxTreeNode *node, TList *consumer) const;
};

#endif /* PWG_EMCAL_ALIEMCALJSONREADER_H_ */
