/***************************************************************************
 *   This macros collects into a common Tree sets of ESD objects 
 *   saved individually in separate files given by "fnames".
 *   Example:
 *        root [0] .L $ALICE_ROOT/STEER/AliESD2Tree.C
 *        root [1] char * names[] = {"rfio:///castor/cern.ch/user/a/aliprod/AliEn-1.0/data/00001/00205.1077678127","rfio:///castor/cern.ch/user/a/aliprod/AliEn-1.0/data/00001/00467.1077678996"}
 *        root [2] AliESD2Tree(names,2)
 *    
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <TError.h>
  #include <Riostream.h>
  #include <TFile.h>
  #include <TTree.h>
  #include <TKey.h>

  #include "AliESD.h"
#endif

static const Char_t *names[]={"AliESDs.root"};
static Int_t n=sizeof(names)/sizeof(Char_t *);

Int_t AliESD2Tree(const Char_t *fnames[]=names, Int_t nf=n) {
  TFile out("AliESDtree.root","recreate");

  TTree *esdTree=new TTree("esdTree","Tree with ESD objects");
  AliESD *event=0;
  esdTree->Branch("ESD","AliESD",&event);  

  for (Int_t i=0; i<nf; i++) {
      TFile * in = TFile::Open(fnames[i]);
      if (!in->IsOpen()) {
         ::Error("AliESD2Tree.C","Can't open file %s",fnames[i]);
         continue;
      }
      TKey *key=0;
      TIter next(in->GetListOfKeys());
      while ((key=(TKey*)next())!=0) {
          event=(AliESD*)key->ReadObj();
          esdTree->Fill();
          delete event;
      }
      in->Close();
  }  
  
  out.cd();
  esdTree->Write();
  delete esdTree;

  return 0;
}
