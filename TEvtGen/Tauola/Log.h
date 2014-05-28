#ifndef __LOG_CLASS_HEADER__
#define __LOG_CLASS_HEADER__

/**
 * This file contains class for logging and filtering output.
 * This header file also includes a debug macro which
 * tracks any possible memory leaks within the program.
 *
 * @author Tomasz Przedzinski
 * @date 14 November 2009
 */

#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <cstring>
#include <list>

using std::stringstream;
using std::string;
using std::streambuf;
using std::ostream;
using std::list;
using std::cout;
using std::endl;

namespace Tauolapp
{

class Log
{
public:
  /** Shows the summary of all messages. */
  static void Summary();

  /** Shows the summary at the end of the program. */
  static void SummaryAtExit()              { atexit(Summary);      }

  /** Adds the decay to the counter. The type is:
      0 - gun, 1 - no mothers & grandmothers, 2 - no mothers, 3 - ok. */
  static void AddDecay(int type);

  /** Four logging entries. Usage:
            Log::Info()<<"Logging some info: "<<8<<" > "<<7.9<<endl;
      Use Log::Info(false) if you don't want the message to be counted.*/
  static ostream& Debug(unsigned short int code=0, bool count=true);
  static ostream& Info(bool count=true);
  static ostream& Warning(bool count=true);
  static ostream& Error(bool count=true);

  /** Turns off or on particular types of messages
      By default, only debugging messages are turned off. */
  static void LogInfo   (bool flag=true)  { iAction=flag;         }
  static void LogWarning(bool flag=true)  { wAction=flag;         }
  static void LogError  (bool flag=true)  { eAction=flag;         }

  static void LogAll    (bool flag=true)  { iAction=wAction=eAction=flag; dRangeS=0; dRangeE=65535; }

        /** Sets the range of debug codes that will be printed.
            By default, the debug messages are turned off. */
  static void LogDebug(unsigned short s=0,unsigned short e=65535)         { dRangeS=s; dRangeE=e;   }

  /** Asserts logical value. If the assertion fails, the default message or 'text'
            will be printed and the program will terminate.
            Program termination can be suppressed by Log::IgnoreFailedAsserts(); */
  static void Assert(bool check, char *text=NULL);

  /** Terminates the program with added default message or 'text'.
            It can be suppressed by Log::IgnoreFatal(); */
  static void Fatal(string text, unsigned short int code=0);
  static void Fatal(unsigned short int code=0)                            { Fatal("",code);       }

  /** Redirects output to log. Redirection can be done for a block of code
      or for one function only. Redirection can be turned off by using
      Log::IgnoreRedirection(); If the target is one of the log streams
      (for example): Log::RedirectOutput( someFunction, Log::Info() );
      You can turn the function's messages off by turning the apropriate
      log entries off. The redirected code will still be executed,
      only messages are redirected. */
  static void RedirectOutput(void (*func)(), ostream& where=*out);
  static void RedirectOutput(ostream& where=*out);
  /** WARNING! If you're redirecting more than one function, do not forget
      to use RevertOutput() afterwards. */
  static void RevertOutput()                      { std::cout.rdbuf(bCout); std::cerr.rdbuf(bCerr); }

  /** Do not exit when Log::Assert() check is false.
      The number of failed asserts will be listed in the summary. */
  static void IgnoreFailedAssert(bool flag=true)                           { asAction=!flag;        }

  /** Ignores redirections of functions' output.
      The function will still be called in a normal way. */
  static void IgnoreRedirection(bool flag=true)                            { rAction=!flag;         }

  /** Do not exit when Log::Fatal() with the code within the provided range is called.
            The number of ignored fatal errors will be listed in the summary. */
  static void IgnoreFatal(unsigned short s=0,unsigned short e=65535) { faRangeS=s; faRangeE=e; }

  /** Change the output of the logged messages.
      Log::SetOutput(cerr);                    //changes the output to cerr
      Log::SetOutput(new ofstream("log.txt")); //changes the output to a file "log.txt" */
  static void SetOutput(ostream *newOut)                                    { out=newOut;           }
  static void SetOutput(ostream &newOut)                                    { out=&newOut;          }

  /** Change the limit of warnings that will be displayed. Set to 0 for no limit. */
  static void SetWarningLimit(int x)                                        { warnLimit=x;          }

protected:
  static streambuf *bCout,*bCerr;
  static ostream *out;
  static stringstream buf;
  static int  warnLimit;
  static int  decays[4];
  static int  dCount,dRangeS,dRangeE,faCount,faRangeS,faRangeE;
  static int  iCount, wCount, eCount, asCount, asFailedCount;
  static bool iAction,wAction,eAction,asAction,rAction;
/**
  Memory leak tracking section. Compile with #define _LOG_DEBUG_MODE_ to turn it on.
  WARNING! Increases execution time significantly. Useful only for debug purposes.
*/
protected:
  typedef struct
  {
    unsigned long address;
    unsigned long size;
    char  file[64];
    unsigned long line;
  } Pointer;
  static list<Pointer*> *PointerList;
public:
#ifdef _LOG_DEBUG_MODE_
  static void NewPointer(unsigned long address,  unsigned long size,  const char *file, unsigned long line)
  {
    if(!PointerList)
    {
      PointerList = new list<Pointer *>();
      atexit(PrintAllocatedPointers);
    }
    Pointer *info = new Pointer();
    info->address = address;
    info->size    = size;
    info->line    = line;
    strncpy(info->file, file, 63);
    PointerList->push_front(info);
  }
  static void DeletePointer(unsigned long address)
  {
    if(!PointerList) return;
    for(list<Pointer*>::iterator i = PointerList->begin(); i!=PointerList->end(); i++)
    {
      if((*i)->address == address)
      {
        PointerList->remove((*i));
        break;
      }
    }
  }
  static bool PointerCompare(Pointer *one, Pointer *two)
  {
    int eq = strcmp(one->file,two->file);
    if(eq<0) return true;
    else if(eq>0) return false;
    return (one->line <= two->line);
  }
  static void PrintAllocatedPointers()
  {
    if(!PointerList) return;
    int pointers=0,buf=0;
    unsigned long total=0;
    char *lastS=" ";
    int lastL=0;
    if(PointerList->size()==0)
    {
      cout<<"----------------------------UNFREED MEMORY POINTERS----------------------------\n";
      cout<<"                                 ... NONE ...\n";
      cout<<"-------------------------------------------------------------------------------\n";
      return;
    }
    PointerList->sort(PointerCompare);
    cout<<"---------------------------UNFREED MEMORY POINTERS---------------------------\n";
    for(list<Pointer*>::iterator i = PointerList->begin(); i!=PointerList->end(); i++)
    {
      total+=(*i)->size;
      ++pointers;
      if(strcmp(lastS,(*i)->file)==0)
      {
        if(lastL==(*i)->line)
        {
          printf("%56s%10lub (%lu)\n"," ",(*i)->size,(*i)->address);
          continue;
        }
      }
      lastS=(*i)->file;
      lastL=(*i)->line;
      printf("%s%n:",(*i)->file,&buf);
      printf("%-*lu%10lub (%lu)\n",55-buf,(*i)->line,(*i)->size,(*i)->address);
    }
    cout<<endl<<total<<"\tbytes"<<endl;
    cout<<pointers<<"\tpointers"<<endl;
    cout<<"-------------------------------------------------------------------------------\n";
  };
#endif //_LOG_DEBUG_MODE_
};

#ifdef _LOG_DEBUG_MODE_

/**
    Redeclare new and delete to use the tracking feature.
    To use __FILE__ and __LINE__ macro efficiently this header file
    should be included in all separately compiled libraries.
*/

inline void* operator new(size_t size, const char *filename, int line)
{
  void *ptr = (void *)malloc(size);
  Log::NewPointer((unsigned long)ptr, size, filename, line);
  return(ptr);
}

inline void  operator delete(void *p)
{
  Log::DeletePointer((unsigned long)p);
  free(p);
}

#define new new(__FILE__, __LINE__)

#endif //_LOG_DEBUG_MODE_

} // namespace Tauolapp
#endif
