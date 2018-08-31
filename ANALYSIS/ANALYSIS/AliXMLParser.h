/// \class AliXMLParser

#ifndef ALIXMLPARSER_H
#define ALIXMLPARSER_H

#ifndef ROOT_TList
#include <TList.h>
#endif

#ifndef ROOT_TSAXParser
#include <TSAXParser.h>
#endif

#ifndef ROOT_TXMLAttr
#include <TXMLAttr.h>
#endif

#ifndef ROOT_TString
#include <TString.h>
#endif

class AliXMLParser {
public:
   AliXMLParser();
   AliXMLParser(const AliXMLParser& obj);
   virtual ~AliXMLParser();
   TList*   GetTreesFromXML(TString file);
   TList*   GetTreesFromURL(TString host);

//   -------------- Slot Functions -----------
   void     OnStartDocument();
   void     OnEndDocument();
   void     OnStartElement(const char*, const TList*);
   void     OnEndElement(const char*);
   void     OnCharacters(const char*);
   void     OnComment(const char*);
   void     OnWarning(const char*);
   void     OnError(const char*);
   void     OnFatalError(const char*);
   void     OnCdataBlock(const char*, Int_t);
//   -----------------------------------------

private:
 
   int      GetEntryIndex(TString entry_name); //Reverse search for index of entry by name
   TList*   fTreeList; ///< List of Trees made from tables
   TString  fTableTag; ///< Identifier for new table
   Bool_t   fInsideTree; ///< True if table_tag is identified
   Int_t    fNumTokens; ///< Number of Attributes
   TList*   fEntries; ///< List of Entries in a table
   TList*   fVal; ///< Corresponding values of Entries in a Table
   Int_t    fNumTrees; ///< Number of Trees
   Bool_t   fError; ///< True if error in XML is encountered

   ClassDef(AliXMLParser,0);

   AliXMLParser& operator=(const AliXMLParser& other);
};

#endif
