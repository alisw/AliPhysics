/*
 * AliJSONReader.cxx
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */

#include <TList.h>
#include <TString.h>

#include "AliJSONData.h"
#include "AliJSONReader.h"

AliJSONSyntaxTreeNode::AliJSONSyntaxTreeNode(const char *name, AliJSONSyntaxTreeNode *mother):
    fName(name),
    fMotherNode(mother),
    fEntries(),
    fDaughters(),
    fOwner(false)
{
}

AliJSONSyntaxTreeNode:: ~AliJSONSyntaxTreeNode() {
  if(fOwner){
    for(std::vector<AliJSONData *>::iterator it = fEntries.begin(); it != fEntries.end(); it++){
      delete *it;
    }
  }
  for(std::vector<AliJSONSyntaxTreeNode *>::iterator it = fDaughters.begin(); it != fDaughters.end(); it++){
    delete *it;
  }
}

void AliJSONSyntaxTreeNode::AddEntry(AliJSONData *entry) {
  fEntries.push_back(entry);
}

AliJSONSyntaxTreeNode *AliJSONSyntaxTreeNode::CreateDaughter(const char *name){
  AliJSONSyntaxTreeNode *daughter = new AliJSONSyntaxTreeNode(name, this);
  fDaughters.push_back(daughter);
  return daughter;
}

void AliJSONSyntaxTreeNode::SetOwner(bool owner) {
  fOwner = owner;
  for(std::vector<AliJSONSyntaxTreeNode *>::iterator it = fDaughters.begin(); it != fDaughters.end(); it++)
    (*it)->SetOwner(owner);
}

AliJSONReader::AliJSONReader() {

}

AliJSONReader::~AliJSONReader() {
}

TList* AliJSONReader::Decode(const char* jsonstring) const {
  /*
   * Decode JSON String
   * 1st create the abstract syntax tree
   * 2nd serialise the abstract syntax tree it into a TList
   */
  AliJSONSyntaxTreeNode * ast = new AliJSONSyntaxTreeNode("", NULL),
      *current = ast;

  TString jsontstring(jsonstring);
  jsontstring.ReplaceAll(" ","");

  // First strip away the left and right brackets
  int first(jsontstring.First('{')+1), last(jsontstring.Last('}'));
  jsontstring = jsontstring(first, last-first+1);
  bool hasNext = jsontstring.Length() > 0;
  while(hasNext){ // Create the abstract syntax tree
    if(jsontstring[0] == '}'){
      current = current->GetMotherNode();
      jsontstring = jsontstring(1, jsontstring.Length()-1);
      if(!(jsontstring.Length() || current)) hasNext = false;
      continue;
    }
    // Find the next separator
    int separator = jsontstring.First(':');
    if(separator == kNPOS){
      hasNext = false;
      break;
    } else{
      TString key = jsontstring(0,separator-1);
      key.ReplaceAll("\"", "");
      jsontstring = jsontstring(separator + 1, jsontstring.Length() - (separator + 1));
      if(jsontstring[0] == '{'){
        // Handle complicated object
        current = current->CreateDaughter(key.Data());
        jsontstring = jsontstring(1,jsontstring.Length()-1);
      } else{
        // Handle simple value
        separator = jsontstring.First(',');
    	if(separator == kNPOS){
          separator = jsontstring.First('}');
	}
        TString value = jsontstring(0, separator -1);
        jsontstring = jsontstring(separator+1, jsontstring.Length() - (separator + 1));
        value.ReplaceAll("\"", "");
        current->AddEntry(new AliJSONData(key.Data(), value.Data()));
      }
    }
  }

  ast->SetOwner(false);

  // Serialise it into a TList
  TList *entries = new TList;
  std::vector<AliJSONData *> &rootnodeentries = ast->GetEntries();
  for(std::vector<AliJSONData *>::iterator it = rootnodeentries.begin(); it != rootnodeentries.end(); it++)
    entries->Add(*it);
  std::vector<AliJSONSyntaxTreeNode *> &daughters = ast->GetDaughters();
  for(std::vector<AliJSONSyntaxTreeNode *>::iterator it = daughters.begin(); it != daughters.end(); it++)
    AddNodeToList(*it, entries);
  // Delete the syntax tree after being done
  delete ast;
  return entries;
}

void AliJSONReader::AddNodeToList(AliJSONSyntaxTreeNode* node, TList* consumer) const {
  TList *entries = new TList;
  entries->SetName(node->GetName());
  std::vector<AliJSONData *> &nodeentries = node->GetEntries();
  std::vector<AliJSONSyntaxTreeNode *> &daughters = node->GetDaughters();
  for(std::vector<AliJSONData *>::iterator it = nodeentries.begin(); it != nodeentries.end(); it++)
    entries->Add(*it);
  for(std::vector<AliJSONSyntaxTreeNode *>::iterator it = daughters.begin(); it != daughters.end(); it++)
    AddNodeToList(*it, entries);
  consumer->Add(entries);
}
