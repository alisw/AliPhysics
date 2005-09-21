// Reads a motif file defined using gassiplex channel and converts
// it (using GCToBerg.dat file) into a motif file using berg pin numbering
// (as used by mapping software).

#include <fstream>
#include <vector>
#include <string>

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

void convertMotif(const char* inputfile, const char* outputfile)
{
  std::vector<int> gc2berg;

  readGCToBerg(gc2berg);

  std::ifstream in(inputfile);
  std::ofstream out(outputfile);

  char line[80];
  int n = 1;
  while ( in.getline(line,80) )
    {
      
      if ( line[0] == '#' )
	{
	  // Comment line. No change.
	  out << line << std::endl;
	}
      else
	{
	  std::string sline(line);
      
	  int gc = atoi(sline.c_str());
	  // assume the first 2 characters only are the gc number,
	  // so replace them with berg number.
	  int berg = gc2berg[gc];
	  out << berg << "\t1\t" << n << "\t-" << std::endl;
	  ++n;
	}
    }
  
  in.close();
  out.close();
}
