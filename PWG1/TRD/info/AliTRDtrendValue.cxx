#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliTRDtrendValue.h"

ClassImp(AliTRDtrendValue)


//____________________________________________
AliTRDtrendValue::AliTRDtrendValue() 
  : TNamed("none", "none")
  ,fAlarmLevel(0)
  ,fValue(0.)
  ,fResponsible()
  ,fNnotifiable(0)
{
  memset(fLimits, 0, 2*(kNlevels+1)*sizeof(Double_t));
  for(Int_t ilevel(kNlevels); ilevel--; ) sprintf(fAlarmMessage[ilevel], " ");
}

//____________________________________________
AliTRDtrendValue::AliTRDtrendValue(Char_t *n, Char_t *t) 
  : TNamed("none", t)
  ,fAlarmLevel(0)
  ,fValue(0.)
  ,fResponsible()
  ,fNnotifiable(0)
{
  TString s(n);
  TObjArray *names(s.Tokenize("_"));
  if(names->GetEntriesFast()!=2){
    AliError(Form("Wrong trend value name format. Trend value name should be of the form \"trendClass_trendValue\" with only one \"_\" character."));
  } else SetName(n);

  memset(fLimits, 0, 2*(kNlevels+1)*sizeof(Double_t));
  for(Int_t ilevel(kNlevels); ilevel--; ) sprintf(fAlarmMessage[ilevel], " ");
}

//____________________________________________
Int_t AliTRDtrendValue::GetAlarmLevel()
{
  // check value against limits and do some more work
  fAlarmLevel=kNlevels-1;
  for(Int_t ilevel(0); ilevel<kNlevels+1; ilevel++)
    if(fValue<fLimits[2*ilevel+1] &&
       fValue>=fLimits[2*ilevel]){ 
      fAlarmLevel = ilevel;
      break;
    }

  return fAlarmLevel;
}

//____________________________________________
const char* AliTRDtrendValue::GetAlarmMessage() const
{
  if(!fAlarmLevel) return "OK";
  else return fAlarmMessage[fAlarmLevel-1];
}

//____________________________________________
const char* AliTRDtrendValue::GetClassName() const
{
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
  TString s(TNamed::GetName());
  TObjArray *names(s.Tokenize("_"));
  if(names->GetEntriesFast()!=2){
    AliError(Form("Wrong trend value name format."));
    return NULL;
  }
  return ((TObjString*)names->At(1))->String().Data();
}

//____________________________________________
const char* AliTRDtrendValue::GetResponsible(Char_t *n, Char_t *mail) const
{
  if(n) sprintf(n, "%s", fResponsible.name);
  if(mail) sprintf(mail, "%s", fResponsible.mail);
  return Form("%s <%s>", fResponsible.name, fResponsible.mail);
}

//____________________________________________
const char* AliTRDtrendValue::GetNotifiable(Int_t in, Char_t *n, Char_t *mail) const
{
  if(in<0||in>=fNnotifiable) return NULL;
  if(n) sprintf(n, "%s", fNotifiable[in].name);
  if(mail) sprintf(mail, "%s", fNotifiable[in].mail);
  return Form("%s <%s>", fNotifiable[in].name, fNotifiable[in].mail);
}

//____________________________________________
void AliTRDtrendValue::SetNotifiable(const Char_t *name, const Char_t *mail)
{
  if(fNnotifiable==kNnotifiable){
    AliWarning(Form("Could not add %s for notification. Only %d persons can be registered for notification.", name, kNnotifiable));
    return;
  }
  sprintf(fNotifiable[fNnotifiable].name, "%s", name);
  sprintf(fNotifiable[fNnotifiable].mail, "%s", mail);
  fNnotifiable++;
}

//____________________________________________
void AliTRDtrendValue::SetResponsible(const Char_t *name, const Char_t *mail) 
{
  sprintf(fResponsible.name, "%s", name);
  sprintf(fResponsible.mail, "%s", mail);
}

//____________________________________________
void AliTRDtrendValue::Print(Option_t */*o*/) const
{
//   name - title
//   value - limits
//   alarm level, message
//   responsible

  printf("    %s [%s] - %s\n", GetValueName(), GetClassName(), GetTitle());
  printf("*** %f limits[%f %f]\n", fValue, fLimits[0], fLimits[1]);
  if(fAlarmLevel){
    printf("*** Alarm level   : %d limits[%f %f]\n", fAlarmLevel, fLimits[2*fAlarmLevel], fLimits[2*fAlarmLevel+1]);
    printf("*** Alarm message : %s\n", GetAlarmMessage());
  }
  printf("*** Responsible %s <%s>\n", fResponsible.name, fResponsible.mail);
  if(fNnotifiable){
    printf("*** Notifiable person(s) ***\n");
    for(Int_t i(0); i<fNnotifiable; i++)
      printf("        %s <%s>\n", fNotifiable[i].name, fNotifiable[i].mail);
  }
}

//____________________________________________
AliTRDtrendValue::AliTRDtrendValueResponsible::AliTRDtrendValueResponsible(Char_t *n, Char_t *m) 
{
  if(n) sprintf(name, "%s", n); else sprintf(name, " ");
  if(m) sprintf(mail, "%s", m); else sprintf(mail, " ");
}
