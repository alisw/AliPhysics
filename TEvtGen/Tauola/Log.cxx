#include <fstream>
#include "Log.h"
using std::streambuf;
using std::stringstream;
using std::ostream;
using std::cout;
using std::cerr;
using std::endl;

namespace Tauolapp
{

list<Log::Pointer*> *Log::PointerList = NULL;

streambuf   *Log::bCout=cout.rdbuf(),*Log::bCerr=cerr.rdbuf();
ostream     *Log::out=&cout;
stringstream Log::buf;
int  Log::warnLimit=100;
int  Log::decays[4] = {0};
int  Log::dCount =0,Log::dRangeS =65535,Log::dRangeE =65534;
int  Log::faCount=0,Log::faRangeS=65535,Log::faRangeE=65534;
int  Log::iCount =0,Log::wCount =0,Log::eCount =0,Log::asCount=0, Log::asFailedCount=0;
bool Log::iAction=1,Log::wAction=1,Log::eAction=1,Log::asAction=1,Log::rAction=1;

void Log::AddDecay(int type)
{
  decays[type]++;
}

ostream& Log::Debug(unsigned short int code, bool count)
{
  if(count) ++dCount;
  if(code>=dRangeS && code<=dRangeE ) return *out<<"DEBUG("<<code<<") from TAUOLA:"<<endl;
  return buf.seekp(0);
}


ostream& Log::Info(bool count)
{
  if(count) ++iCount;
  if(iAction) return *out<<"INFO from TAUOLA:"<<endl;
  return buf.seekp(0);
}


ostream& Log::Warning(bool count)
{
  if(count) ++wCount;

  if(warnLimit>0 && wCount>=warnLimit)
  {
    if(wAction)
    {
      *out<<"WARNING from TAUOLA:"<<endl<<"Limit reached ("<<warnLimit<<"). Warnings suppressed."<<endl;
      wAction=false;
    }
    return buf.seekp(0);
  }

  if(wAction && count) return *out<<"WARNING from TAUOLA:"<<endl;
  if(wAction)          return *out;
  return buf.seekp(0);
}


ostream& Log::Error(bool count)
{
  if(count) ++eCount;
  if(eAction) return *out<<"ERROR from TAUOLA:"<<endl;
  return buf.seekp(0);
}

void Log::Assert(bool check, char *text)
{
  ++asCount;
  if(check) return;

  ++asFailedCount;
  if(text==NULL)  *out<<"ASSERT from TAUOLA:"<<endl<<"Assertion failed. "<<endl;
  else            *out<<"ASSERT from TAUOLA:"<<endl<<"Assertion failed: "<<text<<endl;

  if(asAction) exit(-1);
}

void Log::Fatal(string text,unsigned short code)
{
  ++faCount;
  if(text.size()==0) *out<<"FATAL ERROR from TAUOLA:"<<endl<<"Terminated by a call to Log::Exit();"<<endl;
  else               *out<<"FATAL ERROR from TAUOLA:"<<endl<<text<<endl;
  if(code<faRangeS || code>faRangeE) exit(-1);
}

void Log::RedirectOutput(void (*func)(), ostream& where)
{

  if(!rAction) { func(); return; }
  cout.rdbuf(where.rdbuf());
  cerr.rdbuf(where.rdbuf());
  where<<endl;
  func();
  cout.rdbuf(bCout);
  cerr.rdbuf(bCerr);
}

void Log::RedirectOutput(ostream& where)
{
  if(!rAction) return;
  cout.rdbuf(where.rdbuf());
  cerr.rdbuf(where.rdbuf());
  where<<endl;
}

void Log::Summary()
{
  *out<<"---------------------------- Tauola Log Summary ------------------------------"<<endl;

  // Debug
  *out<<" Debug:   \t";
  if(dRangeS>dRangeE)  *out<<"(OFF)";
  *out<<"\t\t"<<dCount<<"\t";
  if(dRangeS<=dRangeE) *out<<"Debug range: "<<dRangeS<<" - "<<dRangeE;
  *out<<endl;

  // Info
  *out<<" Info:    \t";
  if(!iAction) *out<<"(OFF)";
  *out<<"\t\t"<<iCount<<"\t"<<endl;

  // Warnings
  *out<<" Warnings:\t";
  if(!wAction)
  {
    if(warnLimit>0 && wCount>warnLimit) *out<<"(SUPP.)";
    else                                *out<<"(OFF)";
  }
  *out<<"\t\t"<<wCount<<"\t"<<endl;

  // Errors
  *out<<" Errors:  \t";
  if(!eAction) *out<<"(OFF)";
  *out<<"\t\t"<<eCount<<"\t"<<endl;

  // Counters
  if(asCount || !asAction || faRangeS<faRangeE) cout<<"-----------------------------------"<<endl;
  if(asCount>0)          *out<<" Asserts:                     "<<asCount<<endl;
  if(!asAction)          *out<<" Failed asserts ignored:      "<<asFailedCount<<endl;
  if(faRangeS<=faRangeE) *out<<" Fatal errors ignored:        "<<faCount<<endl;

  cout<<"-----------------------------------"<<endl;
  if(decays[3]) cout<<" Normal decays:                        "<<decays[3]<<endl;
  if(decays[2]) cout<<" Decays without mother:                "<<decays[2]<<endl;
  if(decays[1]) cout<<" Decays without mother & grandmothers: "<<decays[1]<<endl;
  if(decays[0]) cout<<" Decayed using Tauola gun:             "<<decays[0]<<endl;
  *out<<"------------------------------------------------------------------------------"<<endl;
}

} // namespace Tauolapp
