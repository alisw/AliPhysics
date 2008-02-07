#include <iostream>
#include <iomanip>
#include <fstream>

void
dump10bit()
{
  std::ifstream file("FMD_4096.ddl");
  if (!file) { 
    std::cerr << "No such file" << std::endl;
    return;
  }
  size_t i = 0;
  long long w40;
  std::cout << std::setfill('0');
  while (!file.eof()) { 
    if (i % 5 == 0) 
      std::cout << std::setw(7) << i << " ";

    int w10 =  file.get();
    w40     |= ((0x3ff) & w10) << ((i % 5) * 8);
    std::cout << std::hex << std::setw(2) << w10 << " ";

    if (i % 5 == 4) { 
      std::cout << std::hex << std::setw(5) << w40 << std::endl;
      w40 = 0;
    }
    i++;
  }
}

      
