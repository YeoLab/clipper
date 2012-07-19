#include <Python.h>
#include <iostream>
#include <strings.h>
#include <string.h>
using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using namespace std;
#include <ctime>
#include <fstream>
using std::ifstream;
#include <cstdlib> // for exit function
#include <vector>



// new in 1.2: now accepts reads from cin (STDIN)
// this program gets the expected frequency of heights of reads
// randomly distributed r times over a gene of length L
// adapted to C++ by Michael Lovci
extern "C" PyMODINIT_FUNC initpeaks(void); /* Forward */

int main(int argc, char **argv)
{
  Py_SetProgramName(argv[0]);
    
  Py_Initialize();
  
  initpeaks();  
}

/*
length: int, effective length of gene
iterations: int, number of iterations to run
reads: list, list of read lengths
cutoff: double, B-H FDR cutoff
stats: boolean, true prints running time stats

Expects length and iterations to be non-zero. If passed [] for reads returns []
*/
extern "C" PyObject *peaks_shuffle(PyObject *self, PyObject *args)
{  
  int L = 2000; //length
  int r = 1000; //number of times to rerun
  int T = 0 ; //0 or 1 show timing info
  float alpha = .05; //FDR
  PyObject *reads; //list of reads
  
  
  //TODO figure out how to set default values and document this
  
  //parse args
  if(!PyArg_ParseTuple(args, "iiifO", &L, &r, &T, &alpha, &reads)) {
      return NULL;
    }
  
  if (L ==0 || r ==0) {
    PyErr_SetString(PyExc_ValueError, "input values should be non-zero");
    return NULL;
  }
  
  //The equation in the initalizer is the estimated depth of the cutoff, initalizes all values to zero
  int num_reads = PyList_Size(reads); 
  std::vector<long> OBS_CUTOFF(((L / num_reads) + 1) * 100, 0L);
  int redone = 0;
  
  
  
  int GENE[L]; //will store the height at every position of the gene

  //This is height distribution, no reason to have length of gene, arbitrary starting depth choosen 
  std::vector<int> HEIGHT(100, 0L); // will have the number of times you see every height across the gene  HEIGHT[3] = 10 means that there are 10 peaks with a height of 3
  
  srand(time(NULL)); // seed random-number generator with the current time
  for (int iteration = 1; iteration <= r; iteration++){

    //re-initalize height and gene on each iteration
    for(int i = 0; i < L; i++) {
      GENE[i] = 0;
    }
    
    for( int i = 0; i < HEIGHT.size(); i++) {
      HEIGHT[i] = 0;
    }

    //for each read assign it randomly
    for (int i=0; i < num_reads; i++){
      long len = PyInt_AsLong(PyList_GetItem(reads, i));

      //error checking for obviously invalid reads
      if (len > L){
	PyErr_SetString(PyExc_ValueError, "there is a read longer than the gene");
	return NULL;
      }
      int ran;
      
      //Pick a random location witin the gene that the read can map to
      ran = rand() % (L - len - 1); //correct possible positions of reads based on size of gene and size of read
      
      if (ran + len >= L) {
	PyErr_SetString(PyExc_ValueError, "read is assigned past end of gene, this is a bug");
	return NULL;
      }



      //Increment the coverage along the length of that read
      for(int position=ran; position < (ran+len); position++)	{
	GENE[position]++;
      }
    }
    
    int total_num_peaks = 0;
    int max_height = 0;
    for (int j = 0; j < L; j++) {
      
      //if not enough read depth, resize
      if (GENE[j] >= HEIGHT.size()) { 
	  HEIGHT.resize(GENE[j] + 100, 0);
	}

      //now add to total height distribution
      HEIGHT[GENE[j]]++; // simulated height distribution
      if (GENE[j] > 0){
	total_num_peaks++;
      }
      if (GENE[j] > max_height){
	max_height = GENE[j];
      }
    }
    
    //here is where you compute what height is significant.
    double PVAL[max_height];
    std::vector<int> sig_heights;
    PVAL[0] = 1; // set pvalue for height 0 = to 1;

    for (int height = 1; height < max_height; height++){
      
      //Initalize p-value to 1 for all heights
      PVAL[height] = 1;

      //counts number of peaks that are heigher than height
      int bigger_peaks = 0; 
      for (int height2=height; height2 < max_height; height2++){
	bigger_peaks += HEIGHT[height2];
      }

      //convert the integers to doubles and divide.
      double b = bigger_peaks;
      double t =  total_num_peaks;
      double a = b/t;
      PVAL[height] = a;

      if (a < 1 && bigger_peaks > 0){
	sig_heights.push_back(height);
      }  
    }

    double N_heights = sig_heights.size();
    double R = 0.0;
    double correctedP[sizeof(PVAL) / sizeof(PVAL[0])];
    std::vector<int> corr_sig_heights;

    //This is janky need help from mike
    for (int ele = sig_heights.size() - 1; ele >= 0; ele--){
      correctedP[sig_heights[ele]] = PVAL[sig_heights[ele]] * N_heights/(N_heights - R);
      R += 1.0;

      if (correctedP[sig_heights[ele]] < alpha){
	corr_sig_heights.push_back(sig_heights[ele]); // add the current height to corrected significant heights list.
      }
    }

    if (corr_sig_heights.size() == 0 ) { // no heights are significant
      iteration--;
      redone++;

      if (redone > (10*r)){
	break;
      }

      continue;
    }

    //Assign height cutoff to output vector
    int height_cutoff = corr_sig_heights[(corr_sig_heights.size() - 1)]; // height cutoff is the smallest height that is still significant.
    if (height_cutoff > OBS_CUTOFF.size()) {
      OBS_CUTOFF.resize(height_cutoff + 100, 0L);
    }
    OBS_CUTOFF[height_cutoff]++;  
  }
  
  PyObject *returnList = PyList_New(0);
  
  //constuct return list
  for (int cut = 0; cut < OBS_CUTOFF.size(); cut++) {
    PyList_Append(returnList, PyInt_FromLong(OBS_CUTOFF[cut]));
  }
  
  return returnList;
}

static PyMethodDef peaks_methods[] = {
    {"shuffle",             peaks_shuffle,      METH_VARARGS,
     "Return the meaning of everything."},
    {NULL, NULL, 0, NULL}           /* sentinel */
};

extern "C" PyMODINIT_FUNC initpeaks(void) 
{
  PyImport_AddModule("peaks");
  Py_InitModule("peaks", peaks_methods);
}

int usage(char *program_name)
{
  fprintf(stderr,"Usage is %s [options] \n", program_name);
  fprintf(stderr,"file_in is a file with a list of the length of aligned reads (only the part that aligns) \n");
  fprintf(stderr,"Options\n");
  fprintf(stderr," -L <int>   Effective Gene Length\n");
  fprintf(stderr," -r <int>   # of iterations \n");
  fprintf(stderr," -f <int>   input file containing read lengths \n");
  fprintf(stderr," -a <float>   B-H FDR cutoff, default(.05) ... gets very slow as alpha gets smaller\n");
  fprintf(stderr," -T 1/0 default(0) print running time statistics \n");
  fprintf(stderr,"Output:\n");
  fprintf(stderr,"Significance Threshold [tab] Iterations with this observed threshold\n");
  exit(8);
  return -1;
}
