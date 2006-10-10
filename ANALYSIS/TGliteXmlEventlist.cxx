#include "TGliteXmlEventlist.h"

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

#include <TCollection.h>
#include <TMap.h>
#include <TObjString.h>
#include <TXMLEngine.h>
#include <TList.h>


ClassImp(TGliteXmlEventlist)


TGliteXmlEventlist::TGliteXmlEventlist(const char* localfilename):
  TObject(),
  fXmlFile(localfilename),
  fEventList(new TList()),
  fEventListIter(new TIter(fEventList)),
  fCurrent(0)
{
//Andi - please put a comment
  fEventList->SetOwner(kTRUE);
  ReadXML();
}

TGliteXmlEventlist::~TGliteXmlEventlist() {
//Andi - please put a comment
  delete fEventList;
  delete fEventListIter;
}

void
TGliteXmlEventlist::Reset() {
//Andi - please put a comment
  fEventListIter->Reset();
  fCurrent = 0;
}

TMap* TGliteXmlEventlist::Next() {
//Andi - please put a comment
  fCurrent = (TMap*)fEventListIter->Next();
  return fCurrent;
}

void TGliteXmlEventlist::ReadXML() {
//Andi - please put a comment
  TXMLEngine* xml = new TXMLEngine();
#if ROOT_VERSION_CODE < 328704
  xmlDocPointer xdoc = xml->ParseFile(fXmlFile);
  xmlNodePointer xglite = xml->DocGetRootElement(xdoc);
  xmlNodePointer xdtext = xml->GetChild(xglite);
  xmlNodePointer xcollection = xml->GetNext(xdtext);
  
  xmlNodePointer xtext = 0;
  xmlNodePointer xevent = 0;
  xmlNodePointer xtextnext  = 0;
  xmlNodePointer xeventnext = 0;
#else
  XMLDocPointer_t xdoc = xml->ParseFile(fXmlFile);
  XMLNodePointer_t xglite = xml->DocGetRootElement(xdoc);
  XMLNodePointer_t xdtext = xml->GetChild(xglite);
  XMLNodePointer_t xcollection = xml->GetNext(xdtext);

  XMLNodePointer_t xtext = 0;
  XMLNodePointer_t xevent = 0;
  XMLNodePointer_t xtextnext  = 0;
  XMLNodePointer_t xeventnext = 0;
#endif
  Bool_t first_event=kTRUE;
  do {
    if (first_event) {
      xtextnext  = xml->GetChild(xcollection);
      first_event = kFALSE;
    } else {
      xtextnext  = xml->GetNext(xevent);
    }

    if (xtextnext) {
      xeventnext = xml->GetNext(xtextnext);
      xtext = xeventnext;
      xevent = xeventnext;
    } else {
      xtext  = 0;
    }

    if (xevent) {
      if (xml->GetAttr(xevent,"name")) {
	TMap* files = new TMap();
	
	// here is our event
	//	printf("Found xevent: %s\n",xml->GetAttr(xevent,"name"));
	
	Bool_t first_file = kTRUE;
	
#if ROOT_VERSION_CODE < 328704
 	xmlNodePointer xfile = 0;
 	xmlNodePointer xfiletext = 0;

 	xmlNodePointer xfilenext = 0;
 	xmlNodePointer xfiletextnext = 0;
#else
	XMLNodePointer_t xfile = 0;
	XMLNodePointer_t xfiletext = 0;
	
	XMLNodePointer_t xfilenext = 0;
	XMLNodePointer_t xfiletextnext = 0;
#endif
	do {
	  if (first_file) {
	    xfiletextnext = xml->GetChild(xevent);
	    first_file = kFALSE;
	  } else {
	    xfiletextnext = xml->GetNext(xfile);
	  }
	  
	  if (xfiletextnext) {
	    xfilenext = xml->GetNext(xfiletextnext);
	    xfiletext = xfilenext;
	    xfile     = xfilenext;
	  } else {
	    xfile     = 0;
	    xfiletext = 0;
	  }
	  if (xfile) {
	    // here we have an event file  
	    //	    printf("Found file:   %s\n",xml->GetAttr(xfile,"name"));
	    
	    // get the attributes;
	    //	  xml->GetAttr(xfile,"comment");
	    //	  xml->GetAttr(xfile,"date");
	    //	  xml->GetAttr(xfile,"group");
	    //	  xml->GetAttr(xfile,"guid");
	    //	  xml->GetAttr(xfile,"path");
	    //	  xml->GetAttr(xfile,"permissions");
	    //	  xml->GetAttr(xfile,"pfn");
	    //	  xml->GetAttr(xfile,"se");
	    //	  xml->GetAttr(xfile,"size");
	    //	  xml->GetAttr(xfile,"user");
	    
	    Bool_t first_mirror = kTRUE;
	    
#if ROOT_VERSION_CODE < 328704
 	    xmlNodePointer xmirror = 0;
 	    xmlNodePointer xmirrortext = 0;

 	    xmlNodePointer xmirrornext = 0;
 	    xmlNodePointer xmirrortextnext = 0;
#else
	    XMLNodePointer_t xmirror = 0;
	    XMLNodePointer_t xmirrortext = 0;
	    
	    XMLNodePointer_t xmirrornext = 0;
	    XMLNodePointer_t xmirrortextnext = 0;
#endif
	    
	    do {
	      if (first_mirror) {
		xmirrortextnext = xml->GetChild(xfile);
		first_mirror = kFALSE;
	      } else {
		xmirrortextnext = xml->GetNext(xmirror);
	      }
	      
	      if (xmirrortextnext) {
		xmirrornext = xml->GetNext(xmirrortextnext);
		xmirrortext = xmirrornext;
		xmirror     = xmirrornext;
	      } else {
		xmirror     = 0;
		xmirrortext = 0;
	      }
	      if (xmirror) {
		// here we have a file mirror
		xml->GetAttr(xmirror,"name");
		xml->GetAttr(xmirror,"domain");
		xml->GetAttr(xmirror,"latitude");
		xml->GetAttr(xmirror,"longitude");
		xml->GetAttr(xmirror,"location");
		xml->GetAttr(xmirror,"master");
		xml->GetAttr(xmirror,"site");
		xml->GetAttr(xmirror,"rootd");
		const char* master = 0;
		if ( (master = xml->GetAttr(xmirror,"master"))) {
		  if (atoi(master) == 1) {
		    files->Add(new TObjString(xml->GetAttr(xfile,"name")) , new TObjString(xml->GetAttr(xmirror,"rootd")));
		  }
		}
		//		printf("Found File Mirror: %s\n",xml->GetAttr(xmirror,"name"));
	      }
	    } while (xmirror);
	  }
	} while (xfile);
	//	printf("Adding files\n");
	fEventList->Add(files);
      }
    }
  } while ( xtext );
  delete xml;
}    

const char* TGliteXmlEventlist::GetURL(const char* filename) const {
//Andi - please put a comment
  if (fCurrent) {
    TObjString* obj = (TObjString*)fCurrent->GetValue(filename);
    if (obj) {
      if (strlen(obj->GetName())) 
	return (obj->GetName());
      else 
	return 0;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}


void TGliteXmlEventlist::Print(Option_t */*opt*/) const {
//Andi - please put a comment
  printf("Dumping %d elements\n",fEventList->GetSize());
  TIter next(fEventList);
  TMap* filemap;
  Int_t count=0;
  while ( (filemap = (TMap*)next()) ) {
    count++;
    printf("Printing Element %d\n",count);
    filemap->Print();
  }
}



