#include <vector>
#include <utility>
#include <fstream>
#include <string>

void add(std::vector<std::pair<int,int> >& v, int i1, int i2)
{
  v.push_back(std::make_pair<int,int>(i1,i2));
}

std::vector<std::pair<int,int> > make_pattern(const std::string& what)
{
  std::vector<std::pair<int,int> > cols;

  if ( what == "L5" )
    {
      cols.push_back(std::make_pair<int,int>(0,40)); // starting at zero, length = 40
      cols.push_back(std::make_pair<int,int>(0,24));
    }
  else if ( what == "Z1" || what == "Z5" )
    {
      cols.push_back(std::make_pair<int,int>(24,16));
      cols.push_back(std::make_pair<int,int>(0,40));
      cols.push_back(std::make_pair<int,int>(0,8));
    }
  else if ( what == "O9" || what == "O10" || what == "O11" || what == "O12" 
	    || what == "O17" || what == "O18" || what == "O19" )
    {
      add(cols,0,32);
      add(cols,0,32);
    }
  else if ( what == "Z2" )
    {
      cols.push_back(std::make_pair<int,int>(0,8));
      cols.push_back(std::make_pair<int,int>(0,40));
      cols.push_back(std::make_pair<int,int>(24,16));
    }
  else if ( what == "L6" )
    {
      cols.push_back(std::make_pair<int,int>(0,24));
      cols.push_back(std::make_pair<int,int>(0,40));
    }
  else if ( what == "L7" )
    {
      cols.push_back(std::make_pair<int,int>(0,40));
      cols.push_back(std::make_pair<int,int>(16,24));
    }
  else if ( what == "Z3" )
    {
      cols.push_back(std::make_pair<int,int>(0,16));
      cols.push_back(std::make_pair<int,int>(0,40));
      cols.push_back(std::make_pair<int,int>(32,8));
    }
  else if ( what == "Z4" )
    {
      cols.push_back(std::make_pair<int,int>(32,8));
      cols.push_back(std::make_pair<int,int>(0,40));
      cols.push_back(std::make_pair<int,int>(0,16));
    }
  else if ( what == "L8" )
    {
      cols.push_back(std::make_pair<int,int>(16,24));
      cols.push_back(std::make_pair<int,int>(0,40));
    }
  else if ( what == "L3" )
    {
      for ( int i = 0; i < 4; ++i )  add(cols,0,4);
      for ( int i = 0; i < 16; ++i ) add(cols,1,3);
    }
  else if ( what == "O4" )
    {
      for ( int i = 0; i < 16; ++i ) add(cols,0,4);
    }
  else if ( what == "L4" )
    {
      for ( int i = 0; i < 16; ++i ) add(cols,1,3);
      for ( int i = 0; i <  4; ++i ) add(cols,0,4);
    }
  else if ( what == "P3" )
    {
      for ( int i = 0; i <  4; ++i ) add(cols,0,4);
      for ( int i = 0; i <  8; ++i ) add(cols,0,5);
      for ( int i = 0; i <  2; ++i ) add(cols,0,4);
    }
  else if ( what == "Q3" )
    {
      for ( int i = 0; i <  2; ++i ) add(cols,4,1);
      for ( int i = 0; i <  6; ++i ) add(cols,0,5);
      for ( int i = 0; i <  8; ++i ) add(cols,0,4);
    }
  else if ( what == "Q4" )
    {
      for ( int i = 0; i <  8; ++i ) add(cols,0,4);
      for ( int i = 0; i <  6; ++i ) add(cols,0,5);
      for ( int i = 0; i <  2; ++i ) add(cols,4,1);

    }
  else if ( what == "P4" )
    {
      for ( int i = 0; i <  2; ++i ) add(cols,0,4);
      for ( int i = 0; i <  8; ++i ) add(cols,0,5);
      for ( int i = 0; i <  4; ++i ) add(cols,0,4);
    }
  else if ( what == "O1" || what == "O2" || what == "O13" ||
	    what == "O15" || what == "O16" )
    {
      for ( int i = 0; i < 8; ++i ) add(cols,0,8);
    }
  else if ( what == "O5" || what == "O6" || what == "O7" || what == "O8" )
    {
      for ( int i = 0; i < 28; ++i ) add(cols,0,2);
    }
  else if ( what == "L9" )
    {
      add(cols,0,48);
      add(cols,0,16);
    }
  else if ( what == "L10" )
    {
      add(cols,0,16);
      add(cols,0,48);
    }
  else if ( what == "O14" )
    {
      for ( int i = 0; i < 4; ++i )
	{
	  add(cols,0,16);
	}
    }
  else if ( what == "R1" )
    {
 //      add(cols,0,27);
//       add(cols,3,25);
//       add(cols,6,12);
      add(cols,1,27);
      add(cols,0,25);
      add(cols,10,12);
    }
  else if ( what == "R2" )
    {
//       add(cols,0,11);
//       add(cols,0,14);
//       add(cols,13,4);
//       add(cols,-29,11);
//       add(cols,16,24);
      add(cols,29,11);
      add(cols,26,14);
      add(cols,0,11);
      add(cols,-23,4);
      add(cols,0,24);
    }
  else if ( what == "R3" )
    {
//       add(cols,0,13);
//       add(cols,10,6);
//       add(cols,-40,1);
//       add(cols,13,30);
//       add(cols,16,14);
      add(cols,30,13);
      add(cols,2,1);
      add(cols,-27,6);
      add(cols,0,30);
      add(cols,13,14);

    }
  else if ( what == "R4" )
    {
 //      add(cols,0,10);
//       add(cols,0,13);
//       add(cols,11,5);
//       add(cols,-30,15);
//       add(cols,14,21);
      add(cols,35,10);
      add(cols,32,13);
      add(cols,0,15);
      add(cols,-29,5);
      add(cols,10,21);
    }
  else if ( what == "R5" )
    {
//       add(cols,0,11);
//       add(cols,8,6);
//       add(cols,-35,12);
//       add(cols,12,35);
      add(cols,36,11);
      add(cols,0,12);
      add(cols,-33,6);
      add(cols,0,35);
    }
  else if ( what == "R6" )
    {
 //      add(cols,0,8);
//       add(cols,4,8);
//       add(cols,-47,3);
//       add(cols,10,43);
//       add(cols,50,2);
      add(cols,45,8);
      add(cols,3,3);
      add(cols,-41,8);
      add(cols,0,43);
      add(cols,1,2);
    }
  else if ( what == "R7" )
    {
 //      add(cols,0,4);
//       add(cols,0,10);
//       add(cols,7,43);
//       add(cols,-52,5);
//       add(cols,57,2);
      add(cols,55,4);
      add(cols,49,10);
      add(cols,2,5);
      add(cols,-9,43);
      add(cols,0,2);
    }
  else if ( what == "R8" )
    {
//       add(cols,0,6);
//       add(cols,2,54);
//       add(cols,-58,3);
//       add(cols,59,1);
      add(cols,55,6);
      add(cols,0,3);
      add(cols,-5,54);
      add(cols,1,1);
    }
  else if ( what == "R9" )
    {
 //      add(cols,0,1);
//       add(cols,0,3);
//       add(cols,0,60);
      add(cols,59,1);
      add(cols,57,3);
      add(cols,0,60);
    }
  else if ( what == "R10" )
    {
 //      add(cols,61,3);
//       add(cols,0,61);
      add(cols,0,3);
      add(cols,3,61);
    }
  else if ( what == "R11" )
    {
//       add(cols,61,2);
//       add(cols,0,62);
      add(cols,0,2);
      add(cols,1,63);
    }
  else if ( what == "R12" )
    {
//       add(cols,62,1);
//       add(cols,0,63);
      add(cols,0,1);
      add(cols,0,63);
    }
  else if ( what == "R13" )
    {
//       add(cols,1,3);
//       add(cols,1,10);
//       add(cols,0,17);
//       add(cols,0,17);
//       add(cols,0,17);
      add(cols,13,3);
      add(cols,6,10);
      add(cols,0,17);
      add(cols,0,17);
      add(cols,0,17);
    }
  else if ( what == "R14" )
    {
      for ( int i = 0; i < 5; ++i ) add(cols,1,8);
      for ( int i = 0; i < 3; ++i ) add(cols,0,8);
    }
  else if ( what == "R15" )
    {
      for ( int i = 0; i < 3; ++i ) add(cols,9,1);
      for ( int i = 0; i < 4; ++i ) add(cols,1,9);
      for ( int i = 0; i < 2; ++i ) add(cols,1,8);
      for ( int i = 0; i < 1; ++i ) add(cols,0,9);
    }
  else if ( what == "R16" )
    {
      for ( int i = 0; i < 3; ++i ) add(cols,10,1);
      for ( int i = 0; i < 3; ++i ) add(cols,1,10);
      for ( int i = 0; i < 2; ++i ) add(cols,1,9);
      for ( int i = 0; i < 1; ++i ) add(cols,0,10);
      for ( int i = 0; i < 1; ++i ) add(cols,0,3);
    }
  else if ( what == "R17" )
    {
      for ( int i = 0; i < 3; ++i ) add(cols,11,1);
      add(cols,4,8);
      for ( int i = 0; i < 2; ++i ) add(cols,1,11);
      for ( int i = 0; i < 2; ++i ) add(cols,0,11);
      for ( int i = 0; i < 1; ++i ) add(cols,0,9);        
    }
  else if ( what == "R18" )
    {
      for ( int i = 0; i < 2; ++i ) add(cols,13,1);
      add(cols,11,3);
      add(cols,2,12);
      for ( int i = 0; i < 3; ++i ) add(cols,1,13);
      add(cols,0,8);
    }
  else if ( what == "R19" )
    {
      add(cols,9,6);
      for ( int i = 0; i < 2; ++i ) add(cols,1,14);
      for ( int i = 0; i < 2; ++i ) add(cols,0,15);
    }
  else if ( what == "R20" )
    {
      for ( int i = 0; i < 5; ++i ) add(cols,0,11);
      add(cols,2,9);
    }
  else if ( what == "R21" )
    {
      add(cols,0,2);
      for ( int i = 0; i < 5; ++i ) add(cols,0,11);
      add(cols,4,7);
    }
  else if ( what == "R22" )
    {
      add(cols,1,4);
      for ( int i = 0; i < 2; ++i ) add(cols,1,11);
      for ( int i = 0; i < 3; ++i ) add(cols,0,12);
      add(cols,10,2);
    }
  else if ( what == "R23" )
    {
      add(cols,0,10);
      for ( int i = 0; i < 4; ++i ) add(cols,0,12);
      add(cols,6,6);      
    }
  else if ( what == "R24" )
    {
      add(cols,0,7);
      for ( int i = 0; i < 4; ++i ) add(cols,0,13);
      add(cols,8,5);
    }
  else if ( what == "R25" )
    {
      add(cols,0,9);
      for ( int i = 0; i < 3; ++i ) add(cols,0,14);
      add(cols,1,13);
    }
  else if ( what == "R26" )
    {
      add(cols,1,2);
      for ( int i = 0; i < 2; ++i ) add(cols,1,7);
      for ( int i = 0; i < 6; ++i ) add(cols,0,8);
    }
  else if ( what == "R27" ) 
    {
      add(cols,1,53);
      add(cols,0,11);
    }
  else if ( what == "R28" ) 
    {
      add(cols,12,43);
      add(cols,0,21);
    }
  else if ( what == "R29" ) 
    {
      add(cols,22,34);
      add(cols,0,30);
    }
  else if ( what == "L11" ) 
    {
      add(cols,12,26);
      add(cols,0,38);
    }
  else if ( what == "R31" ) 
    {
      add(cols,0,19);
      add(cols,12,45);
    }
  else if ( what == "R32" ) 
    {
      add(cols,0,15);
      add(cols,11,49);
    }
  else if ( what == "R33" ) 
    {
      add(cols,0,13);
      add(cols,11,51);
    }
  else if ( what == "L12" || what == "L14" || what == "L16" ) 
    {
      add(cols,0,48);
      add(cols,32,16);
    }
  else if ( what == "R34" ) 
    {
      add(cols,0,32);
      add(cols,16,32);
    }
  else if ( what == "L13" || what == "L15" || what == "L17" ) 
    {
      add(cols,0,16);
      add(cols,0,48);
    }
  else if ( what == "R30" ) 
    {
      add(cols,9,15);
      add(cols,5,22);
      add(cols,0,27);
    }
  else if ( what == "R35" )
    {
      for ( int i = 0; i < 3; ++i ) add(cols,9,1);
      for ( int i = 0; i < 2; ++i ) add(cols,8,2);
      add(cols,7,3);
      for ( int i = 0; i < 2; ++i ) add(cols,6,4);
      add(cols,4,6);
      for ( int i = 0; i < 4; ++i ) add(cols,0,10);      
    }
  else if ( what == "R36" )
    {
      for ( int i = 0; i < 6; ++i ) add(cols,0,10);      
      add(cols,6,4);
    }
  else if ( what == "R37" )
    {
      add(cols,0,6);
      for ( int i = 0; i < 4; ++i ) add(cols,0,10);      
      for ( int i = 0; i < 2; ++i ) add(cols,1,9);      
    }
  else if ( what == "R38" )
    {
      for ( int i = 0; i < 3; ++i ) add(cols,0,1);      
      for ( int i = 0; i < 2; ++i ) add(cols,0,2);      
      add(cols,0,3);
      for ( int i = 0; i < 2; ++i ) add(cols,0,4);      
      for ( int i = 0; i < 7; ++i ) add(cols,0,6);      
      add(cols,0,4);
    }
  else if ( what == "R39" )
    {
      add(cols,4,2);
      for ( int i = 0; i < 8; ++i ) add(cols,0,6);            
      for ( int i = 0; i < 2; ++i ) add(cols,0,7);            
    }
  else if ( what == "R40" )
    {     
      add(cols,29,7);
      add(cols,15,21);
      add(cols,0,36);
    }
  else if ( what == "R41" )
    {     
      add(cols,0,7);
      add(cols,0,21);
      add(cols,0,36);      
    }
  else if ( what == "R42" )
    {     
      add(cols,0,32);
      add(cols,16,32);
    }
  else if ( what == "I1" )
    {
      add(cols,0,64);
    }
  else if ( what == "L18" )
    {
      add(cols,0,40);
      add(cols,32,8);
    }
  else if ( what == "L19" )
    {
      add(cols,0,48);
      add(cols,32,16);
    }
  else if ( what == "L20" )
    {
      add(cols,32,16);
      add(cols,0,48);
    }
  else if ( what == "D1" ) // dummy big motif of 16x80 to test
    {
      for ( int i = 0; i < 16; ++i )
	{
	  add(cols,0,80);
	}
    }
  return cols;
}

void makePadPos(const char* padposfile /*="padPosTest.dat"*/, const char* type)
{
  std::ofstream out(padposfile);

  int n = 1;

  std::vector<std::pair<int,int> > cols = make_pattern(type);

  int c = 0;

  for ( size_t i = 0; i < cols.size(); ++i )
    {
      std::pair<int,int> p = cols[i];
      if ( p.first < 0 ) --c;
      for ( int j = abs(p.first); j < abs(p.first)+p.second; ++j )
	{
	  out << n << "\t" << c << "\t" << j << std::endl;
	  ++n;
	}
      ++c;
    }
  out.close();
}

int main(int argc, char** argv)
{
  makePadPos(argv[1],argv[2]);
}
