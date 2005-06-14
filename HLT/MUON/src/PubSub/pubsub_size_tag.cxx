/*
usage: progname -c | -x <filename>

the -c switch inserts the size (in number of 32 bit double words)
into a 32 bit double word at the beginning of the file.
the -x switch removes it.

Author: Gareth de Vaux
Email:  devaux@lhc.phy.uct.ac.za | dhlt@lordcow.org
*/

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <unistd.h>

using namespace std;

int include(char *filename); // writes the size tag
int remove(char *filename);  // removes it


class queue {
  
  class node {
    public:
      char ch;
      node *next;
    };

  node *head, *tail;

  public:
    queue() {
      head = NULL;
      tail = NULL;
      }

    ~queue() {
      node *n = head;
        
      while (n != NULL) {
        head = head->next;
        delete n;
        n = head;
        }
          
      head = NULL;
      tail = NULL;
      }

    void add(char data) {
      node *tmp = new node;
      assert(tmp != NULL);
      tmp->ch = data;
      tmp->next = NULL;
      if (tail != NULL) tail->next = tmp;
      tail = tmp;
      if (head == NULL) head = tail;
      }

    char rem() {
      assert(head != NULL);
      node *tmp = head;
      char tmpdata = head->ch;
      head = head->next;
      if (head == NULL) tail = NULL;
      delete tmp;
      return tmpdata;
      }
  };



int main(int argc, char *argv[])
{
  if (argc!=3) {
    cout << "usage: " << argv[0] << " -c | -x <filename>\n";
    return(1);
    }

  if (!strcmp(argv[1],"-c")) {
    return(include(argv[2]));
    }

  else {
    if (!strcmp(argv[1],"-x")) {
      return(remove(argv[2]));
      }

    else {
      cout << "usage: " << argv[0] << " -c | -x <filename>\n";
      return(1);
      }
    }
}


int include(char *filename)
{
  char ch;
  int i;
  queue q;

  fstream file(filename, ios::in | ios::out | ios::binary | ios::ate);
  
  if (!file) {
    cout << "can't open " << filename << '\n';
    return(1);
    }

  int doubleword = file.tellg()/4;  // find size of file in number of "32 bit double words"
  
  file.seekg(0, ios::beg);  // since file was opened with 'ios::ate' to find the size

  for (i=0; file.get(ch) && i<4; i++) {  // backup the first 4 bytes in the queue
    q.add(ch);
    }

  file.seekp(0, ios::beg);  // point to the beginning of the file again

  for (i=0; i<4; i++) {     // write the "32 bit double word"/size to the beginning of the file
    file.put( doubleword % 0x100 );
    doubleword /= 0x100;
    }

  while(file.get(ch)) {       // now shift all the bytes down the file, using the
    q.add(ch);                // queue as temporary storage
    file.seekp(-1, ios::cur);
    file.put(q.rem());
    }

  for (i=0; i<4; i++) {  // and dump the rest of the queue at the end;
    file.put(q.rem());
    }

  file.close();
  return(0);
}


int remove(char *filename)
{
  char ch;

  fstream file(filename, ios::in | ios::out | ios::binary | ios::ate);
  if (!file) {
    cout << "can't open " << filename << '\n';
    return(1);
    }

  int filesize = file.tellg(); // find size of file in bytes
      
  file.seekg(4, ios::beg); // skip the get pointer past the "32 bit double word"
  
  while(file.get(ch)) {  // shift all the bytes up the file
    file.seekp(-5, ios::cur);
    file.put(ch);
    file.seekg(4, ios::cur);
    }

  file.close();
  truncate(filename,filesize-4);
  return(0);
}
