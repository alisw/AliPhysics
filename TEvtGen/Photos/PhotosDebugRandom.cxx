#include <iostream>
#include "PhotosDebugRandom.h"
#include "Log.h"
using std::ios_base;

namespace Photospp
{

int    PhotosDebugRandom::i97_saved      = 0;
int    PhotosDebugRandom::j97_saved      = 0;
double PhotosDebugRandom::uran_saved[97] = { 0 };
double PhotosDebugRandom::cran_saved     = 0;

void PhotosDebugRandom::saveState()
{
  i97_saved=i97;
	j97_saved=j97;
  cran_saved=cran;

	for(int i=0;i<97;i++) uran_saved[i]=uran[i];
}

void PhotosDebugRandom::restoreState()
{
  i97=i97_saved;
	j97=j97_saved;
  cran=cran_saved;

	for(int i=0;i<97;i++) uran[i]=uran_saved[i];
}

void PhotosDebugRandom::setState(int i, int j, double c, double list[97])
{
  i97=i;
  j97=j;
  cran=c;
  for(int i=0;i<97;i++) uran[i]=list[i];
}

void PhotosDebugRandom::setSaveState(int i, int j, double c, double list[97])
{
  i97_saved=i;
  j97_saved=j;
  cran_saved=c;
  for(int i=0;i<97;i++) uran_saved[i]=list[i];
}

void PhotosDebugRandom::print()
{
  int                coutPrec = cout.precision(18);
	ios_base::fmtflags flags    = cout.setf(ios_base::floatfield);

  Log::RedirectOutput(Log::Info());

  cout<<"double uran_state[97] = { ";
  for(int i=0;i<96;i++) cout<<uran[i]<<", ";
  cout<<uran[96]<<" };"<<endl<<endl;
  cout<<"PhotosDebugRandom::setState( "<<i97<<", "<<j97<<", "<<cran<<", uran_state );"<<endl;

  Log::RevertOutput();

	// Revert output stream flags and precision
	cout.precision(coutPrec);
	cout.flags    (flags);
}

} // namespace Photospp
