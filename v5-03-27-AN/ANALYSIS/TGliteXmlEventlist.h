#ifndef TGLITEXMLEVENTLIST_H
#define TGLITEXMLEVENTLIST_H

///////////////////////////////////////////////////////////////////////////////////////
// class to read gLite XML collections
// Autor: Dr. A.-J. Peters - CERN 2004/ALICE   Mail-to: Andreas.Joachim.Peters@cern.ch
///////////////////////////////////////////////////////////////////////////////////////

// example: 
//   ---------------------------------------------------------------------------------
//   TGliteXmlEventlist* evlist = new TGliteXmlEventlist("/tmp/fileset.rxml");
//   evlist->Reset()
//   // loop over all events
//   while (evlist->Next()) {
//     printf("URL for file AliESDs.root is: \n", evlist->GetURL("AliESDs.root"));
//   }
//   delete evlist;
//   ---------------------------------------------------------------------------------


#include <TObject.h>
#include <TString.h>
class TList;
class TIter;
class TMap;

class TGliteXmlEventlist : public TObject {
public:
  TGliteXmlEventlist(const char* localfilename);
  virtual ~TGliteXmlEventlist();
  void Reset();
  TMap* Next();
  const char* GetURL(const char* name) const ;
  void        Print(Option_t* opt) const;
private:
  TGliteXmlEventlist(const TGliteXmlEventlist & src);
  TGliteXmlEventlist & operator=(const TGliteXmlEventlist & src);

  TString     fXmlFile;//Andi - please put a comment
  TList*      fEventList;//Andi - please put a comment
  TIter*      fEventListIter;//Andi - please put a comment
  TMap*       fCurrent;//Andi - please put a comment

  virtual void ReadXML();

  ClassDef(TGliteXmlEventlist,0);
};

#endif
