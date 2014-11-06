/*
 * AliEMCALJSONReader.cxx
 *
 *  Created on: 06.11.2014
 *      Author: markusfasel
 */
#include <TList.h>
#include <TString.h>

#include <AliEMCALJSONReader.h>

AliEMCALJSONReader::AliEMCALJSONReader() {
  // TODO Auto-generated constructor stub

}

AliEMCALJSONReader::~AliEMCALJSONReader() {
  // TODO Auto-generated destructor stub
}

TList* AliEMCALJSONReader::Decode(const char* jsonstring) const {
  /*
   * Decode JSON String
   * 1st create the abstract syntax tree
   * 2nd serialise the abstract syntax tree it into a TList
   */
  AliEMCALJSONSyntaxTreeNode * ast = new AliEMCALJSONSyntaxTreeNode("", NULL),
      *current = ast;

  TString jsontstring(jsonstring);
  jsontstring.ReplaceAll(" ","");

  // First strip away the left and right brackets
  int first(jsontstring.First('{')+1), last(jsontstring.Last('}')-1);
  jsontstring = jsontstring(first, last-first);
  bool hasNext = jsontstring.Length() > 0;
  while(hasNext){ // Create the abstract syntax tree
    if(jsontstring[0] == '}'){
      current = current->GetMotherNode();
      jsonstring = jsontstring(1, jsontstring.Length()-1);
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
      if(jsontstring[0] == "{"){
        // Handle complicated object
        current = current->CreateDaughter(key.Data());
        jsontstring = jsontstring(1,jsontstring.Length()-1);
      } else{
        // Handle simple value
        separator = jsontstring.First(',');
        TString value = jsontstring(0, separator -1);
        jsontstring = jsontstring(separator+1. jsontstring.Length() - (separator + 1));
        value.ReplaceAll("\"", "");
        current->AddEntry(new AliEMCALConfigurationObject(key.Data(), value.Data()));
      }
    }
  }

  ast->SetOwner(false);

  // Serialise it into a TList
  TList *entries = new TList;
  std::vector<AliEMCALConfigurationObject *> &rootnodeentries = ast->GetEntries();
  for(std::vector<AliEMCALConfigurationObject *>::iterator it = rootnodeentries.begin(); it != rootnodeentries.end(); it++)
    entries->Add(*it);
  std::vector<AliEMCALJSONSyntaxTreeNode *> &daughters = ast->GetDaughters();
  for(std::vector<AliEMCALJSONSyntaxTreeNode *>::iterator it = daughters.begin(); it != daughters.end(); it++)
    AddNodeToList(*it, entries);
  return entries;
}

void AliEMCALJSONReader::AddNodeToList(AliEMCALJSONSyntaxTreeNode* node, TList* consumer) const {
  TList *entries = new TList;
  entries->SetName(node->GetName());
  std::vector<AliEMCALConfigurationObject *> &nodeentries = node->GetEntries();
  std::vector<AliEMCALJSONSyntaxTreeNode *> &daughters = node->GetDaughters();
  for(std::vector<AliEMCALConfigurationObject *>::iterator it = nodeentries.begin(); it != nodeentries.end(); it++)
    entries->Add(*it);
  for(std::vector<AliEMCALJSONSyntaxTreeNode *>::iterator it = daughters.begin(); it != daughters.end(); it++)
    AddNodeToList(*it, entries);
  consumer->Add(entries);
}
