#include "TMath.h"
#include "TH1.h"
#include "TArrayI.h"
#include "TArrayF.h"

#define MAX_TREE_HT 100

struct MH_Node
{
  int character;
  float frequency;
  struct MH_Node *l, *r;
};
 
 
struct M_Heap
{
  unsigned size;
  unsigned space;
  struct MH_Node **array;
};
 
struct MH_Node* newNode(int character, float frequency)
{
  struct MH_Node* temp = (struct MH_Node*) malloc(sizeof(struct MH_Node));
  temp->l = temp->r = NULL;
  temp->character = character;
  temp->frequency = frequency;
  return temp;
}
 
 
struct M_Heap* createM_Heap(unsigned space)
{
  struct M_Heap* M_Heap = (struct M_Heap*) malloc(sizeof(struct M_Heap));
  M_Heap->size = 0;
  M_Heap->space = space;
  M_Heap->array = (struct MH_Node**)malloc(M_Heap->space * sizeof(struct MH_Node*));
  return M_Heap;
}
 
 
void swapMH_Node(struct MH_Node** a, struct MH_Node** b)
{
  struct MH_Node* t = *a;
  *a = *b;
  *b = t;
}
 
 
void M_Heapify(struct M_Heap* M_Heap, int idx)
{
  // arrange in increasing frequency order
  for (int i=idx;i<M_Heap->size;i++) {
    for (int j=i+1;j<M_Heap->size;j++) {
      if (M_Heap->array[i]->frequency > M_Heap->array[j]->frequency)  
	swapMH_Node(&M_Heap->array[i], &M_Heap->array[j]);
    }
  }
}
 
int isSizeOne(struct M_Heap* M_Heap)
{
  return (M_Heap->size == 1);
}
 
 
struct MH_Node* extractMin(struct M_Heap* M_Heap)
{
  struct MH_Node* temp = M_Heap->array[0];
  M_Heap->array[0] = M_Heap->array[M_Heap->size - 1];
  --M_Heap->size;
  M_Heapify(M_Heap, 0);
  return temp;
}
 
 
void insertM_Heap(struct M_Heap* M_Heap, struct MH_Node* MH_Node)
{
  //  printf("Insert %f\n",MH_Node->frequency);
  int i = M_Heap->size - 1;
  if (i<0) {
    M_Heap->array[M_Heap->size++] = MH_Node;
    return;
  }
  M_Heap->array[M_Heap->size] = MH_Node;
  ++M_Heap->size;
  M_Heapify(M_Heap, 0);
}
 
 
void buildM_Heap(struct M_Heap* M_Heap)
{
  M_Heapify(M_Heap, 0);
}
 
 
void printArr(int arr[], int n)
{
  int i;
  for (i = 0; i < n; ++i) printf("%d", arr[i]);
  printf("\n");
}
 
 
int isLeaf(struct MH_Node* root)
{
  return !(root->l) && !(root->r) ;
}
 
 
struct M_Heap* createAndBuildM_Heap(int character[], float frequency[], int size)
{
  int i;
  struct M_Heap* M_Heap = createM_Heap(size);
  for (i = 0; i < size; ++i) M_Heap->array[i] = newNode(character[i], frequency[i]);
  M_Heap->size = size;
  buildM_Heap(M_Heap);
  return M_Heap;
}
 
struct MH_Node* buildHuffmanTree(int character[], float frequency[], int size)
{
  struct MH_Node *l, *r, *top;
  struct M_Heap* M_Heap = createAndBuildM_Heap(character, frequency, size);
  while (!isSizeOne(M_Heap)) {
    l = extractMin(M_Heap);
    r = extractMin(M_Heap);
    top = newNode('$', l->frequency + r->frequency);
    top->l = l;
    top->r = r;
    insertM_Heap(M_Heap, top);
  }
  return extractMin(M_Heap);
}

float printCodes(struct MH_Node* root, int arr[], int top)
{
  float totsiz = 0;
  if (root->l)
    {
      arr[top] = 0;
      totsiz += printCodes(root->l, arr, top + 1);
    }
  if (root->r)
    {
      arr[top] = 1;
      totsiz += printCodes(root->r, arr, top + 1);
    }
  if (isLeaf(root))
    {
      printf("%d: (%f) ", root->character, root->frequency);
      totsiz += top*root->frequency;
      printArr(arr, top);
    }
  return totsiz;
}

void HuffmanCodes(int character[], float frequency[], int size)
{
  struct MH_Node* root = buildHuffmanTree(character, frequency, size);
  int arr[MAX_TREE_HT], top = 0;
  int defNbits = 1 + TMath::Log2(size);
  double freqNorm = 0;
  for (int i=0;i<size;i++) {
    freqNorm += frequency[i];
  }
  float packSize = printCodes(root, arr, top);
  //
  printf("Brute force: %d vs %f coded\n",defNbits,packSize/freqNorm);

}

void HuffmanCodes(TH1* histo, int maxBin=-1)
{
  int nb = histo->GetNbinsX();
  if (maxBin>1 && nb>maxBin) nb = maxBin;
  TArrayI dat(nb);
  TArrayF frq(nb);  
  for (int i=0;i<nb;i++) {
    dat[i] = i;
    frq[i] = histo->GetBinContent(i+1);
  }
  HuffmanCodes(dat.GetArray(),frq.GetArray(),nb);
  //
}
