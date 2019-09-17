// ProgressLog.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#ifndef PROGRESSLOG_H
#define PROGRESSLOG_H

#include <ctime>
#include <sys/times.h>
#include <unistd.h>

namespace Pythia8 {

using namespace std;

// ProgressLog is a simple helper class to monitor the progress of a run.
// When used in the main event loop, it will with suitably (logarithmically)
// spaced intervals, print out one line with information about the number
// of events generated, two estimates (based on instantaneous and average
// CPU consumption) of when the run will be completed, the host on which
// the program is run and its process number.

class ProgressLog {

public:

  // Create an object for monitoring the progress of a run with NIn iterations.
  ProgressLog( long NIn, int maxsec = 1200) : secstep(maxsec) {
    init(NIn);
  }

  // Intermittently write out a line of progress information, giving
  // the current iteration (in the interval [0:N[ ).
  void operator()(int cnt) {
    tick(cnt + 1, N);
  }

  // Intermittently write out a line of progress information using the
  // internal counter.
  void operator()() {
    tick(++count, N);
  }

  // Intermittently write out a line of progress information giving
  // the current iteration and the total number of iterations.
  void tick(long i, long n) {
    if ( !statusTime(i, n) ) return;

    double fcpui = fclock();
    time_t timei = time(0);
    double ftime0 = time0;
    double ftime1 = time1;
    double ftimei = timei;
    double eff = 1.0;
    if ( ftimei > ftime1 && fcpui > fcpu1 )
      eff = (fcpui-fcpu1)/(ftimei-ftime1);
    if ( eff >= 1.0 ) eff = 0.999999;
    int ieff = 100*eff;
    double eff0 = 1.0;
    if ( ftimei > ftime0 && fcpui > fcpu0 )
      eff0 = (fcpui-fcpu0)/(ftimei-ftime0);
    if ( eff0 >= 1.0 ) eff0 = 0.999999;
    int ieff0 = 100*eff0;
    double fcpun = fcpu0 + (n*(fcpui-fcpu0))/i;
    time_t timen = (time_t)(ftimei + (fcpun-fcpui)/eff + 30.0);
    time_t timen0 = (time_t)(ftimei + (fcpun-fcpui)/eff0 + 30.0);
    char date[1024];
    char daten[1024];
    char daten0[1024];
    strftime(date,1024,"%y.%m.%d %H:%M",localtime(&timei));
    strftime(daten,1024,"%H:%M",localtime(&timen));
    strftime(daten0,1024,"%H:%M",localtime(&timen0));
    long ii = i;
    if ( n - i < n/10 ) ii = i - n;
    time_t dayn = (timen - timei)/86400;
    time_t dayn0 = (timen0 - timei)/86400;

    ostream & os = cout;

    if ( dayn <= 0 && dayn0 <= 0 ) {
      os << date << " " << setw(8) << ii << "/" << setw(9);
      os.setf(ios::left, ios::adjustfield);
      os << n << " etc:   " << daten << "[";
      os.setf(ios::right, ios::adjustfield);
      os << setw(2) << ieff << "%]   " << daten0 << "[" << ieff0 << "%] "
         << host << ":" << pid << endl << flush;
    } else {
      os << date << " " << setw(8) << ii << "/" << setw(9);
      os.setf(ios::left, ios::adjustfield);
      os << n << " etc: " << dayn << "+" << daten << "[";
      os.setf(ios::right, ios::adjustfield);
      os << setw(2) << ieff << "%] "
         << dayn0 << "+" << daten0 << "[" << ieff0 << "%] "
         << host << ":" << pid << endl << flush;
    }

    fcpu1 = fcpui;
    time1 = timei;

  }

  // Interface to the system time information.
  double fclock() {
    struct tms tmsbuf;
    times(&tmsbuf);
    double d =
      tmsbuf.tms_utime+tmsbuf.tms_stime+tmsbuf.tms_cutime+tmsbuf.tms_cstime;
    d /= sysconf(_SC_CLK_TCK);
    return d;
  }

  // Check if this is a good time to print out a status line.
  bool statusTime(long i, long n) const {
    if ( i <= 0 ) return false;
    if ( i == n ) return true;
    if ( i > n/2 ) i = n-i;
    while ( i >= 10 && !(i%10) ) i /= 10;
    if ( i == 1 || i == 2 || i == 5 ) return true;
    if ( secstep > 0 && time(0) > time1 + secstep ) return true;
    return false;
  }

  // Initialise the basic engine.
  void init(long n) {
    N = n;
    count = 0;
    fcpu0 = fcpu1 = fclock();
    time0 = time1 = time(0);
    char name[1024];
    gethostname(name,1024);
    host = name;
    if ( host.find(".") != string::npos )
      host = host.substr(0, host.find("."));
    pid = getpid();
    char date[1024];
    strftime(date,1024,"%y.%m.%d %H:%M",localtime(&time0));
    ostream & os = cout;
    os << date << "        0/" << setw(9);
    os.setf(ios::left, ios::adjustfield);
    os << n;
    os.setf(ios::right, ios::adjustfield);
    os << " Initializing...                "
       << host << ":" << pid << endl << flush;
  }

private:

  // If larger than 0, a status line will be written every secstep
  // second.
  int secstep;

  // The clock when the run was started.
  time_t time0;

  // The cpu clock when the run was started.
  double fcpu0;

  // The clock the last time a status line was written out.
  time_t time1;

  // The cpu clock the last time a status line was written out.
  double fcpu1;

  // The host on which we are running.
  string host;

  // The pid of the current process.
  pid_t pid;

  // The number of iterations
  long N;

  // The number of iterations so far
  long count;

};

}

#endif
