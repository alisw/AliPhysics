#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include <cctype>
#include <iomanip>

bool GCToBerg()
{
  bool debug = false;

  std::ifstream in("bergToGC.dat");
  if (!in)
    {
      std::cerr << "Could not open bergToGC.dat" << std::endl;
      return false;
    }
  std::map<int,int> gctoberg;

  int berg;
  std::string manu;;
  while ( in >> berg >> manu )
    {
      if (debug) std::cout << berg << " -> " << manu;
      if ( std::isdigit(manu[0]) )
	{
	  int imanu = std::atoi(manu.c_str());
	  if ( debug ) std::cout << " = " << imanu;
	  gctoberg[imanu] = berg;
	}
      if ( debug ) std::cout << std::endl;
    }

  in.close();

  std::map<int,int>::const_iterator it;

  for ( it = gctoberg.begin(); it != gctoberg.end(); ++it )
    {
      std::cout << std::setw(2) << it->first 
		<< std::setw(4) << it->second
		<< std::endl;
    }
  return true;
}
