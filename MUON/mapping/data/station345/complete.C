#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cctype>

void complete(const char* file)
{
  std::string ofile(file);
  std::string ifile(file);

  ifile += ".old";

  std::ostringstream cmd;

  cmd << "mv " << ofile << " " << ifile;

  system(cmd.str().c_str());

  int n = 1;

  char line[80];

  std::ifstream in(ifile.c_str());
  std::ofstream out(ofile.c_str());

  while ( in.getline(line,80) )
    {
      if ( !std::isdigit(line[0]) )
	{
	  out << line << std::endl;
	}
      else
	{
	  out << line << "\t1\t" << n << "\t-" << std::endl;
	  n++;
	}
    }
  in.close();
  out.close();
}
