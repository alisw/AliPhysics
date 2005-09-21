//Given a list of gassiplex (manas) channel ranges, defines
//a motif.dat file.

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <algorithm>

void readGCToBerg(std::vector<int>& gc2berg)
{
  gc2berg.clear();

  std::ifstream in("GCToBerg.dat");
  
  int gc, berg;

  while ( in >> gc >> berg )
    {
      gc2berg.push_back(berg);
    }

  in.close();
}

void defineMotif(const char* motifName)
{
  std::vector<int> gc2berg;

  readGCToBerg(gc2berg);

  char line[80];

  std::ostringstream ofile;

  ofile << "motif" << motifName << ".dat";

  std::ofstream out(ofile.str().c_str());

  out << "# Motif " << motifName << std::endl
      << "#" << std::endl
      << "#connecteur_berg kapton padname not_used" << std::endl
      << "#for slats there's no kapton connector, so it's always 1 (zero make the reader" << std::endl
      << "#abort, so it's not a valid value here)." << std::endl
      << "#" << std::endl;

  int n = 1;
  while ( std::cin.getline(line,80) )
    {
      std::string sline(line);

      std::string::size_type pos = sline.find_first_of('-');
      int i1 = atoi(sline.c_str());
      int i2 = i1;
      if ( pos != std::string::npos )
	{
	  i2 = atoi(sline.substr(pos+1).c_str());
	}
      // assume the first 2 characters only are the gc number,
      // so replace them with berg number.
      if ( i1 < i2 ) 
	{
	  for ( int i = i1; i <= i2; ++i )
	    {
	      int berg = gc2berg[i];
	      out << berg << "\t1\t" << n << "\t-" << std::endl;
	      ++n;
	    }
	}
      else
	{
	  for ( int i = i1; i >= i2; --i )
	    {
	      int berg = gc2berg[i];
	      out << berg << "\t1\t" << n << "\t-" << std::endl;
	      ++n;
	    }
	}
    }
  
  out.close();
}

int main(int argc, char** argv)
{
  defineMotif(argv[1]);
  return 0;
}
