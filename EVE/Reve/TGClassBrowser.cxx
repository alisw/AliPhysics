
#include "TROOT.h"
#include "TSystem.h"
#include "TGClient.h"
#include "TGListTree.h"
#include "TGLayout.h"
#include "TSystemDirectory.h"
#include "TGMimeTypes.h"
#include "TFile.h"
#include "TKey.h"
#include "TInterpreter.h"
#include "TClass.h"
#include "TDataMember.h"
#include "TMethod.h"
#include "TMethodArg.h"

#include "TGClassBrowser.h"

ClassImp(TGClassBrowser)

//______________________________________________________________________________
TGClassBrowser::TGClassBrowser(const TGWindow *p, UInt_t w, UInt_t h) :
      TGMainFrame(p, w, h)

{
   SetCleanup(kDeepCleanup);
   fCanvas = new TGCanvas(this, 100, 100);
   fListTree = new TGListTree(fCanvas, kHorizontalFrame);
   AddFrame(fCanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop | 
                kLHintsExpandX | kLHintsExpandY));
   
   TGListTreeItem *item = fListTree->AddItem(0,"Classes");
   fListTree->Connect("DoubleClicked(TGListTreeItem *, Int_t)",
      "TGClassBrowser", this, "DoubleClicked(TGListTreeItem *, Int_t)");
   fClassIcon = gClient->GetPicture("class.png");
   fMemberIcon = gClient->GetPicture("member.png");
   fMethodIcon = gClient->GetPicture("method.png");

   MapSubwindows();
   Resize(GetDefaultSize());
   MapWindow();
   TClass *cl;
   TIter nextcl(gROOT->GetListOfClasses());
   while ((cl = (TClass *)nextcl())) {
      TGListTreeItem *itm = fListTree->AddItem(item, cl->GetName(), 
                                               fClassIcon, fClassIcon);
      itm->SetTipText(Form("Decl: %s : %d\nImpl: %s : %d",
                      cl->GetDeclFileName(), cl->GetDeclFileLine(),
                      cl->GetImplFileName(), cl->GetImplFileLine()));
   }
}

//______________________________________________________________________________
TGClassBrowser::~TGClassBrowser()
{
   // Destructor.

   delete fListTree;
   Cleanup();
}

//______________________________________________________________________________
void TGClassBrowser::DisplayClass(TGListTreeItem *, const TString &)
{
   // display details of ROOT class

}

//______________________________________________________________________________
void TGClassBrowser::DoubleClicked(TGListTreeItem *item, Int_t /*btn*/)
{
   // Process double clicks in TGListTree.

   TClass *cl;
   TIter nextcl(gROOT->GetListOfClasses());
   if (item == fListTree->GetFirstItem()) {   
      while ((cl = (TClass *)nextcl())) {
         if (!fListTree->FindChildByName(item, cl->GetName())) {
            TGListTreeItem *itm = fListTree->AddItem(item, cl->GetName(), 
                                                     fClassIcon, fClassIcon);
            itm->SetTipText(Form("Decl: %s : %d\nImpl: %s : %d",
                            cl->GetDeclFileName(), cl->GetDeclFileLine(),
                            cl->GetImplFileName(), cl->GetImplFileLine()));
         }
      }
      fListTree->ClearViewPort();
      return;
   }
   cl = (TClass*)gROOT->GetListOfClasses()->FindObject(item->GetText());
   if (cl) {
      TClass *bc;
      TDataMember *dm;
      TMethod *m;
      
      TIter nextbc(cl->GetListOfBases());
      while ((bc = (TClass *) nextbc())) {
         fListTree->AddItem(item, bc->GetName(), fClassIcon, fClassIcon);
      }
     
      TIter nextdm(cl->GetListOfDataMembers());
      while ((dm = (TDataMember *) nextdm())) {
         if (dm->IsEnum()) continue;
         fListTree->AddItem(item, Form("%s %s",dm->GetFullTypeName(), dm->GetName()), 
                            fMemberIcon, fMemberIcon);
      }

      TIter nextm(cl->GetListOfMethods());
      while ((m = (TMethod *) nextm())) {
         TString method = Form("%s %s(", m->GetReturnTypeName(), m->GetName());
         TIter nextarg(m->GetListOfMethodArgs());
         TMethodArg *a;
         while ((a = (TMethodArg*) nextarg())) {
            method += a->GetFullTypeName();
            method += ", ";
         }
         method.Remove(TString::kTrailing,' ');
         method.Remove(TString::kTrailing,',');
         method += ")";
         fListTree->AddItem(item, method.Data(), fMethodIcon, fMethodIcon);
      }
   }
   fListTree->ClearViewPort();
}

