/// \file AliXMLParser.cxx

#include "AliXMLParser.h"

#include <Riostream.h>
#include "TTree.h"
#include "TBranch.h"
#include "TWebFile.h"
#include <iostream>
#include <cstdlib>
#include <TList.h>
#include <TString.h>

using std::cerr;
using std::cout;
using std::endl;

ClassImp(AliXMLParser);

AliXMLParser::AliXMLParser():
   fTreeList(0),
   fTableTag(0),
   fInsideTree(kFALSE),
   fNumTokens(0),
   fEntries(0),
   fVal(0),
   fNumTrees(0),
   fError(kFALSE)
{
   /// Default Contructor

   fEntries = new TList();
   fVal = new TList();
   fTreeList = new TList();
   fInsideTree = kFALSE;
   fNumTrees = 0;
   fNumTokens = -1;
   fError = kFALSE;

   fEntries->SetOwner(kTRUE);
   fVal->SetOwner(kTRUE);
   fTreeList->SetOwner(kTRUE);
}

AliXMLParser::AliXMLParser(const AliXMLParser& obj):
   fTreeList(0),
   fTableTag(0),
   fInsideTree(kFALSE),
   fNumTokens(0),
   fEntries(0),
   fVal(0),
   fNumTrees(0),
   fError(kFALSE)
{
   /// Copy Contructor

   fTreeList = obj.fTreeList;
   fTableTag = obj.fTableTag;
   fInsideTree = obj.fInsideTree;
   fNumTokens = obj.fNumTokens;
   fEntries = obj.fEntries;
   fVal = obj.fVal;
   fNumTrees = obj.fNumTrees;
   fError = obj.fError;
}

AliXMLParser& AliXMLParser::operator=(const AliXMLParser& other)
{
   /// Assignment

   if(this != &other)
   {
      fTreeList = other.fTreeList;
      fTableTag = other.fTableTag;
      fInsideTree = other.fInsideTree;
      fNumTokens = other.fNumTokens;
      fEntries = other.fEntries;
      fVal = other.fVal;
      fNumTrees = other.fNumTrees;
      fError = other.fError;
   }
   return *this;
}

void AliXMLParser::OnStartDocument() {}

int AliXMLParser::GetEntryIndex(TString entry_name) //Finds the attribute index by name
{
   int i;
   for(i=0;i<=fNumTokens;i++)
   {
      if(entry_name.EqualTo(*(TString *)(fEntries->At(i))))
         return i;
   }
   return -1;
}

void AliXMLParser::OnStartElement(const char *name, const TList *attributes) //Stores Parsed XML in a Tree
{
   int n,i;
   TString tree_symbol = "roottreename";
   TString var(name);
   TString temp("");
   TXMLAttr *attr;
   TIter next(attributes);
   TBranch *branch;
   TTree *curr_tree;

   if((attr = (TXMLAttr*) next()))
   {
      if((!fInsideTree) && tree_symbol.EqualTo(attr->GetName())) //If tag found, Enter a new tree
      {
         fTableTag = name;
         fNumTrees++;
         fTreeList->Add((TObject *)(new TTree(attr->GetValue(),attr->GetValue())));
         fInsideTree = kTRUE;
      }
      else if(fInsideTree)
      {
         curr_tree=(TTree *)(fTreeList->At(fNumTrees-1));
         do
         {
            if(!(curr_tree->GetBranch(attr->GetName())))
          
            {
               TString var_type = attr->GetName();
               var_type += "/C";
               fNumTokens++;

               fEntries->Add((TObject *)(new TString(attr->GetName())));

               fVal->Add((TObject *)(new TString(temp)));

               //Adding Attributes as branches of TTree
               branch = curr_tree->Branch(attr->GetName(), (void *)(((TString *)(fVal->At(fNumTokens)))->Data()), var_type);
               n = curr_tree->GetEntries();
               if(n > 0)
               {
                  (*(TString *)(fVal->At(fNumTokens))) = temp;
                  for(i=0;i<n;i++)
                     branch->Fill();
               }

            }

            (*(TString *)(fVal->At(GetEntryIndex(attr->GetName())))) = attr->GetValue();

           
         } while ((attr = (TXMLAttr*) next()));
         curr_tree->Fill(); //Fill the tree with attribute values
         for(i=0;i<=fNumTokens;i++)
         {
            (*(TString *)(fVal->At(i))) = temp;
         }
      }
   }
}

      
   
void AliXMLParser::OnEndElement(const char *name)
{
   /// if tag closes, re-initialize everything for a new tree

   if(fTableTag.EqualTo(name))
   {
      fInsideTree = kFALSE;
      fEntries->Clear("nodelete");
      fVal->Clear("nodelete");
      fNumTokens=-1;
   }
}

void AliXMLParser::OnCharacters(const char *) {} 

void AliXMLParser::OnComment(const char *) {}

void AliXMLParser::OnWarning(const char *warning) 
{
   cout << "Warning in XML: " << warning << endl;
}

void AliXMLParser::OnError(const char *error) 
{ 
   fError = kTRUE;
   cerr << "Error in XML: " << error << endl;
}

void AliXMLParser::OnFatalError(const char *error)
{ 
   fError = kTRUE;
   cerr << "FatalError in XML: " << error << endl;
}

void AliXMLParser::OnCdataBlock(const char *, Int_t ) {}

void AliXMLParser::OnEndDocument() {}

AliXMLParser::~AliXMLParser()
{
   delete fEntries;
   delete fVal;
   delete fTreeList;
}


TList* AliXMLParser::GetTreesFromXML(TString file) //Returns List of Trees by parsing local XML file
{
   TSAXParser *saxParser = new TSAXParser();
   saxParser->ConnectToHandler("AliXMLParser", this);
   if(!TFile::Open(file)) 
   {
      cerr << "File not found" << endl;
      return 0x0;
   }
   saxParser->ParseFile(file);
   delete saxParser;
   if (fError) return 0x0;
   return fTreeList;
}

TList *AliXMLParser::GetTreesFromURL(TString host) //Returns List of Trees by parsing remote XML file
{
   char *buf;
   int size;
   TSAXParser *saxParser = new TSAXParser();
   saxParser->ConnectToHandler("AliXMLParser", this);
   host += "&filetype=raw";
   TFile *file = TFile::Open(host);
   if (!file)
   {
      cerr << "Invalid URL" << endl;
      return 0x0;
   }
   size=file->GetSize();
   buf = (char *)malloc(size+1);
   file->ReadBuffer(buf, size);
   saxParser->ParseBuffer(buf,size);
   delete saxParser;
   free(buf);
   if (fError) return 0x0;
   return fTreeList;
}
