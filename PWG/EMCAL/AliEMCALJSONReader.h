/*
 * AliEMCALJSONReader.h
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */

#ifndef PWG_EMCAL_ALIEMCALJSONREADER_H_
#define PWG_EMCAL_ALIEMCALJSONREADER_H_

#include <string>
#include <vector>

class TList;

class AliEMCALJSONSyntaxTreeNode{

  AliEMCALJSONSyntaxTreeNode(const char *name, AliEMCALJSONSyntaxTreeNode *mother):
    fName(name),
    fMotherNode(mother),
    fEntries(),
    fDaughters(),
    fOwner(false)
  {}

  ~AliEMCALJSONSyntaxTreeNode() {
    if(fOwner){
      for(std::vector<AliEMCALConfigurationObject *>::iterator it = fEntries.begin(); it != fEntries.end(); it++){
        delete *it;
      }
    }
    for(std::vector<AliEMCALJSONSyntaxTreeNode *>::iterator it = fDaughters.begin(); it != fDaughters.end(); it++){
      delete *it;
    }
  }

  void AddEntry(AliEMCALConfigurationObject *entry) { fEntries.push_back(entry); }

  AliEMCALJSONSyntaxTreeNode *CreateDaughter(const char *name){
    AliEMCALJSONSyntaxTreeNode *daughter = new AliEMCALJSONSyntaxTreeNode(name, this);
    fDaughters.push_back(daughter);
    return daughter;
  }

  void SetOwner(bool owner = true) {
    fOwner = owner;
    for(std::vector<AliEMCALJSONSyntaxTreeNode *>::iterator it = fDaughters.begin(); it != fDaughters.end(); it++)
      (*it)->SetOwner(owner);
  }
  AliEMCALJSONSyntaxTreeNode *GetMotherNode() const { return fMotherNode; }
  std::vector<AliEMCALJSONSyntaxTreeNode *> &GetDaughters() const { return fDaughters; };
  std::vector<AliEMCALConfigurationObject *> &GetEntries() const { return fEntries; }
  const char *GetName() const { return fName; }

private:
  AliEMCALJSONSyntaxTreeNode *fMotherNode;
  std::vector<AliEMCALConfigurationObject *> fEntries;
  std::vector<AliEMCALJSONSyntaxTreeNode *> fDaughters;
  bool fOwner;
}

class AliEMCALJSONReader {
public:
  AliEMCALJSONReader() {}
  virtual ~AliEMCALJSONReader() {}

  TList *Decode(const char *jsosnstring) const;

private:
  void AddNodeToList(AliEMCALJSONSyntaxTreeNode *node, TList *consumer) const;
};

#endif /* PWG_EMCAL_ALIEMCALJSONREADER_H_ */
