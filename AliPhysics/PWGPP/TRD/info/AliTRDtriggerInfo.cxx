////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Mergable Trigger List                                                 //
//                                                                        //
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////



#include "TObjArray.h"
#include "TObjString.h"
#include "TH1.h"

#include "AliLog.h"

#include "AliTRDtriggerInfo.h"

ClassImp(AliTRDtriggerInfo)
//____________________________________________
AliTRDtriggerInfo::AliTRDtriggerInfo()
  :TObject()
  ,fTriggerList(NULL)
  ,fHisto(NULL)
{
//  Constructor. Reset all fields.
  memset(fTriggerSel, 0, kTriggerListSize *sizeof(Bool_t));
  memset(fTriggerStat, 0, kTriggerListSize *sizeof(Int_t));
}

//____________________________________________
AliTRDtriggerInfo::~AliTRDtriggerInfo()
{
//  Destructor. 
  if(fHisto) delete fHisto;
  if(fTriggerList){ fTriggerList->Delete(); delete fTriggerList;}
}

//____________________________________________
Int_t AliTRDtriggerInfo::Add(const Char_t *trigger, Int_t nstat, Bool_t select)
{
// Add trigger named "trigger" to the list. If trigger is not on the list it is add.
// The statistics of this trigger is set to "nstat" [default=1]
// On return the position incremented

  Int_t itrig(-1), nt(0);
  if((nt=GetNTriggers()) && (itrig = GetTrigger(trigger))>=0){
    fTriggerSel[itrig]  = select;
    fTriggerStat[itrig]+= nstat;
  } else {
    if(!nt) fTriggerList = new TObjArray;
    fTriggerList->Add(new TObjString(trigger));
    fTriggerStat[nt] = nstat;
    fTriggerSel[nt]  = select;
    nt++;
  }
  return nt;
}

//____________________________________________
void AliTRDtriggerInfo::Draw(Option_t *opt)
{
// Draw trigger statistics. Via parameter "opt" one can select a sublist of trigger names

  Int_t nt(0);
  if(!(nt = GetNTriggers())) return;
  if(fHisto && fHisto->GetNbinsX() != nt){ delete fHisto; fHisto=NULL;}
//
  if(!fHisto) fHisto = new TH1F("hTriggerStat", "Trigger Statistics;TRIGGER;Freq. [%]", nt, 0.5, nt+0.5);
  TAxis *ax(fHisto->GetXaxis());  fHisto->Reset();
  ax->SetTitleOffset(6.5); ax->CenterTitle();
  fHisto->SetFillStyle(3001);fHisto->SetFillColor(kGreen);
  fHisto->SetBarWidth(0.8);fHisto->SetBarOffset(0.1);
//
  Int_t ibin(0);
  if(strcmp(opt, "")!=0){
    TString sopt(opt);
    TObjArray *optTriggers = sopt.Tokenize(" ");
    for(Int_t iopt(0); iopt<optTriggers->GetEntries(); iopt++){
      const Char_t *trigger = ((TObjString*)(*optTriggers)[iopt])->GetName();
      Int_t itrig(-1);
      if((itrig = GetTrigger(trigger))<0) continue;
      ibin++;
      ax->SetBinLabel(ibin, trigger);
      fHisto->SetBinContent(ibin, fTriggerStat[itrig]);
    }
    optTriggers->Delete(); delete optTriggers;
  } else {
    for(Int_t jtrig(nt); jtrig--;){
      ibin++;
      ax->SetBinLabel(ibin, ((TObjString*)(*fTriggerList)[jtrig])->GetName());
      fHisto->SetBinContent(ibin, fTriggerStat[jtrig]);
    }
  }
  ax->SetRange(1, ibin);
  fHisto->Scale(1.e2/fHisto->Integral());
  fHisto->Draw("hbar2");

  return;
}

//____________________________________________
Int_t AliTRDtriggerInfo::GetNTriggers() const
{
// Return no of trigger classes
  return fTriggerList?fTriggerList->GetEntries():0;
}

//____________________________________________
const char* AliTRDtriggerInfo::GetTrigger(Int_t it) const
{
// Return trigger by index. 0 in case it is not found
  if(it<0||it>=GetNTriggers()) return 0;
  return ((TObjString*)(*fTriggerList)[it])->GetName();
}

//____________________________________________
Int_t AliTRDtriggerInfo::GetTrigger(const char *trigger) const
{
// Return trigger index by name. -1 in case trigger not registred

  for(Int_t jtrig(GetNTriggers()); jtrig--;){
    if(((TObjString*)(*fTriggerList)[jtrig])->String().CompareTo(trigger)==0) return jtrig;
  }
  return -1;
}

//____________________________________________
Long64_t AliTRDtriggerInfo::Merge(TCollection* list)
{
// Merge list of trigger statistics. Add new triggers if neccessary

  Int_t imerge(1);
  AliTRDtriggerInfo *o(NULL);
  TIterator *it(list->MakeIterator());
  while((o  = (AliTRDtriggerInfo*)it->Next())){
    for(Int_t jtrig(o->GetNTriggers()); jtrig--;) Add(((TObjString*)(*o->fTriggerList)[jtrig])->GetName(), o->fTriggerStat[jtrig], o->fTriggerSel[jtrig]);
    imerge++;
  }
  return imerge;
}
  
//____________________________________________
void AliTRDtriggerInfo::Print(Option_t *opt) const
{
// Print trigger statistics. If opt="a" detailed trigger/statistics is also provided

  Int_t nt(GetNTriggers());
  Int_t stat(0); for(Int_t jtrig(0); jtrig<nt; jtrig++) stat+=fTriggerStat[jtrig];
  printf("TriggerClasses[%3d] Stat[%5d]\n", nt, stat);
  if(!nt || strcmp(opt,"a")!=0) return;

  for(Int_t jtrig(0); jtrig<nt; jtrig++) printf("  %3d \"%s\" [%5d]\n", jtrig, ((TObjString*)(*fTriggerList)[jtrig])->GetName(), fTriggerStat[jtrig]);
}


