/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-----------------------------------------------------------------
//           Implementation of the AliXMLCollection class
//   This is the class that creates XML collections after querying the tags
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

//ROOT
#include <Riostream.h>
#include <TString.h>
#include <TEventList.h>

#include "AliXMLCollection.h"

ClassImp(AliXMLCollection)

//___________________________________________________________________________
  AliXMLCollection::AliXMLCollection() :
    TObject(),
    fCollectionName(),
    fout()
{
  //Default constructor
}

//___________________________________________________________________________
AliXMLCollection::~AliXMLCollection() {
  //Destructor
}

//___________________________________________________________________________
Bool_t AliXMLCollection::WriteHeader() {
  //Creates the xml output file

  TString xmlName = fCollectionName;
  xmlName += ".xml";

  TString collectionHeader = "<collection name=";
  collectionHeader += "\"";
  collectionHeader += fCollectionName;
  collectionHeader += "\"";
  collectionHeader += ">";
  
  // Open the output stream
  fout.open(xmlName);
  fout<<"<?xml version=\"1.0\"?>\n";
  fout<<"<tags>\n";
  fout<<"  "<<collectionHeader<<"\n";  

  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliXMLCollection::WriteBody(Int_t counter, const char* guid, const char* turl, TEventList *list) {
  //Writes the body of the xml collection
  TString listline;
  for(Int_t i = 0; i < list->GetN(); i++) {
    listline += list->GetEntry(i);
    listline += ",";
  }  
  listline = listline(0,listline.Length()-1);

  TString line0 = "<event name=\"";
  line0 += counter;
  line0 += "\">";
  
  TString line1 = "<file name=\"AliESDs.root\" ";
  line1 += "guid=\"";
  line1 += guid;
  line1 += "\" ";
  line1 += "turl=\"";
  line1 += turl;
  line1 += "\" ";
  line1 += "evlist=\"";
  line1 += listline;
  line1 += "\"";
  line1 += " />";
  
  fout<<"    "<<line0<<"\n";
  fout<<"      "<<line1<<"\n";
  fout<<"    </event>\n";
  
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliXMLCollection::Export() {
  //Closes the stream
  fout<<"  "<<"</collection>\n";
  fout<<"</tags>\n";

  fout.close();

  return kTRUE;
}
