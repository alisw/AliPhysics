////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Trend Value Incapsulation                                             //
//                                                                        //
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliTRDtrendValue.h"

ClassImp(AliTRDtrendValue)


//____________________________________________
AliTRDtrendValue::AliTRDtrendValue() 
  : TNamed("none", "none")
  ,fValue(0.)
  ,fSigma(1.)
  ,fResponsible(NULL)
{
//  Constructor. Reset all fields.
  //memset(fAlarmMessage, 0, kNlevels*sizeof(Char_t*));
  memset(fNotifiable, 0, kNnotifiable*sizeof(TNamed*));
}

//____________________________________________
AliTRDtrendValue::AliTRDtrendValue(const Char_t *n, const Char_t *t) 
  : TNamed("none", t)
  ,fValue(0.)
  ,fSigma(1.)
  ,fResponsible(NULL)
{
//  Constructor. Define name and title for trend variable.
  TString s(n);
  TObjArray *names(s.Tokenize("_"));
  if(names->GetEntriesFast()!=2){
    AliError(Form("Wrong trend value name format. Trend value name should be of the form \"trendClass_trendValue\" with only one \"_\" character."));
  } else SetName(n);

  //memset(fAlarmMessage, 0, kNlevels*sizeof(Char_t*));
  memset(fNotifiable, 0, kNnotifiable*sizeof(TNamed*));
}

//____________________________________________
AliTRDtrendValue::AliTRDtrendValue(const AliTRDtrendValue &ref)
  : TNamed(ref)
  ,fValue(ref.fValue)
  ,fSigma(ref.fSigma)
  ,fResponsible(NULL)
{
  if(ref.fResponsible) fResponsible = new TNamed(*ref.fResponsible);
  //memset(fAlarmMessage, 0, kNlevels*sizeof(Char_t*));
  //for(Int_t ia(0); ia<kNlevels; ia++) SetAlarm(ia, ref.fAlarmMessage[ia]);
  memset(fNotifiable, 0, kNnotifiable*sizeof(TNamed*));
  Int_t in(0);
  while(ref.fNotifiable[in]){
    fNotifiable[in] = new TNamed(*(ref.fNotifiable[in]));
    in++;
  }
}

//____________________________________________
AliTRDtrendValue& AliTRDtrendValue::operator/=(const AliTRDtrendValue &n)
{
  fValue-=n.fValue;
  if(n.fSigma>0.) fValue/=n.fSigma;
  return *this;
}

//____________________________________________
const char* AliTRDtrendValue::GetAlarmMessage(Int_t ns) const
{
// Check if value triggered alarm
  if(ns<0 || ns>kNlevels) return NULL;
  else return "not defined";//fAlarmMessage[ns];
}

//____________________________________________
const char* AliTRDtrendValue::GetClassName() const
{
// Check task to which value belong
  TString s(TNamed::GetName());
  TObjArray *names(s.Tokenize("_"));
  if(names->GetEntriesFast()!=2){
    AliError(Form("Wrong trend value name format."));
    return NULL;
  }

  return ((TObjString*)names->At(0))->String().Data();
}

//____________________________________________
const char* AliTRDtrendValue::GetValueName() const
{
// value name
  TString s(TNamed::GetName());
  TObjArray *names(s.Tokenize("_"));
  if(names->GetEntriesFast()!=2){
    AliError(Form("Wrong trend value name format."));
    return NULL;
  }
  return ((TObjString*)names->At(1))->String().Data();
}

//____________________________________________
const char* AliTRDtrendValue::GetResponsible() const
{
// Get responsible with name and mail
  if(!fResponsible) return NULL;
  return Form("%s <%s>", fResponsible->GetName(), fResponsible->GetTitle());
}

//____________________________________________
const char* AliTRDtrendValue::GetNotifiable(Int_t in) const
{
// Get noticible person "in" with name and mail
  if(in<0||in>kNnotifiable) return NULL;
  if(!fNotifiable[in]) return NULL;
  return Form("%s <%s>", fNotifiable[in]->GetName(), fNotifiable[in]->GetTitle());
}

//____________________________________________
void AliTRDtrendValue::SetAlarm(Int_t level, Char_t */*m*/)
{
// define alarm message for "level"
  if(level<0||level>=kNlevels){
    AliWarning(Form("Alarm level[%d] out of range [0 %d]", level, kNlevels-1));
    return;
  }
  //fAlarmMessage[level] = StrDup(m);
}

//____________________________________________
void AliTRDtrendValue::SetNotifiable(const Char_t *name, const Char_t *mail)
{
// add noticible person to DB
  Int_t n(0); while(GetNotifiable(n)) n++;
  if(n>=kNnotifiable-1){
    AliWarning(Form("Could not add %s for notification. There are already %d persons registered for notification.", name, kNnotifiable-1));
    return;
  }
  fNotifiable[n] = new TNamed(name, mail);
}

//____________________________________________
void AliTRDtrendValue::SetResponsible(const Char_t *name, const Char_t *mail, Option_t *opt) 
{
// set responsible person for trend
  if(fResponsible){
    if(strcmp(opt, "u")!=0){
      AliWarning(Form("Responsible already set %s <%s>", fResponsible->GetName(), fResponsible->GetTitle()));
      return;
    }else{
      AliWarning(Form("Old responsible %s <%s> replaced by %s <%s>", fResponsible->GetName(), fResponsible->GetTitle(), name, mail));
      new(fResponsible) TNamed(name, mail);
    }
  } else fResponsible = new TNamed(name, mail);
}

//____________________________________________
void AliTRDtrendValue::Print(Option_t */*o*/) const
{
//   name - title
//   value - limits
//   alarm level, message
//   responsible

  printf("    %s [%s] - %s\n", GetValueName(), GetClassName(), GetTitle());
  printf("*** %f +- %f\n", fValue, fSigma);
  printf("*** Responsible %s <%s>\n", fResponsible?fResponsible->GetName():"", fResponsible?fResponsible->GetTitle():"");
  printf("*** Notifiable person(s) ***\n");
  Int_t in(0);
  while(fNotifiable[in]){
    printf("        %s <%s>\n", fNotifiable[in]->GetName(), fNotifiable[in]->GetTitle());
    in++;
  }
/*  printf("*** Alarm messages \n");
  for(in=0; in<kNlevels; in++) printf("*** Alarm [%d] : %s\n", in, fAlarmMessage[in]?fAlarmMessage[in]:"not defined");*/
}
