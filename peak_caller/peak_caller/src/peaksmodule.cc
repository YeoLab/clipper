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
  if (num_reads == 0) {
        PyErr_SetString(PyExc_ValueError, "List must not be null");
	return NULL;
  }
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
      if (len + 1 >= L){
	continue;
	//PyErr_SetString(PyExc_ValueError, "there is a read longer than the gene");
	//return NULL;
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
  
  PyObject *returnList = PyList_New(OBS_CUTOFF.size());
  
  //constuct return list
  for (int cut = 0; cut < OBS_CUTOFF.size(); cut++) {
    PyList_SetItem(returnList, cut, PyInt_FromLong(OBS_CUTOFF[cut]));
  }

  return returnList;
}

/* Find contigous (within margin) regions that have reads, the area between covered locations is defined as regions without any coverage

Input: data - wiggle track in list form each value is the coverage at that location
margin - distance between section

Output:
A list of strings in the form "start_location|stop_location"

TODO: Modify to allow for thresholded margins"

*/
extern "C" PyObject *peaks_find_sections(PyObject *self, PyObject *args) {
  std::vector<PyObject*> sections(0); //vector of sections because appending to a python list takes a very long time
  PyObject *wiggle; //list of reads
  int margin;
  int start = 0;
  int stop = 0;
  bool in_section = false;
  int gap = 0;
  int loc = 0; //initalize outside because we need to use for the end catch

  //parse args
  if(!PyArg_ParseTuple(args, "Oi", &wiggle, &margin)) {
      return NULL;
    }
  
  int wiggle_size = PyList_Size(wiggle); 
  
  //walk along the wiggle track looking for gaps wider than the margin, when that happens output the region
  for (; loc < wiggle_size; loc++) { 
    double cur_wiggle_value = PyFloat_AsDouble(PyList_GetItem(wiggle, loc));

    //If the current value is non-zero assume you are in a section and if you are not in a section set this location as the start of the section 
    if (cur_wiggle_value > 0) {
      gap = 0;
      
      //if not in section mark this location as the first part of a section
      if(!in_section) {
	start = loc;
      }

      in_section = true;
      //sets the new start after 
      
    } else {
      gap += 1;

      //sets new section if any only if we just left a section 
      if(in_section && gap > margin ) {
	in_section = false;
	stop = loc - gap + 1; //sets the stop to the last location that a real value has been seen
	
	//adds section to list
	PyObject *section = PyTuple_New(2);
	PyTuple_SetItem(section, 0, PyInt_FromLong(start));
	PyTuple_SetItem(section, 1, PyInt_FromLong(stop));
     
	sections.push_back(section);
      }
    }
  }

  //catch last potental section
  if ( in_section ) {
 
    PyObject *section = PyTuple_New(2);
    PyTuple_SetItem(section, 0, PyInt_FromLong(start));
    PyTuple_SetItem(section, 1, PyInt_FromLong(loc - (gap + 1)));
   
    sections.push_back(section);

  }
  
  PyObject *returnList = PyList_New(sections.size());
  
  //transform vector to pylist
  for (int i = 0; i < sections.size(); i++) {
    PyList_SetItem(returnList, i, sections[i]);
  }
  
  return returnList; 
}

/* fast version of reads to wiggle */
extern "C" PyObject *peaks_readsToWiggle_pysam(PyObject *self, PyObject *args) {
  
  
  //Define argunments passed in
  PyObject *reads; //list of reads
  int tx_start;
  int tx_end;
  char *keepstrand;
  char *usePos;
  PyObject *fractional_input;
  bool makeFractional; //flag to return either whole number coverage or coverage normalzied by length of reads
  //parse args
  if(!PyArg_ParseTuple(args, "OiissO", &reads, &tx_start, &tx_end, &keepstrand, &usePos, &fractional_input)) {
      return NULL;
  }
  
  //error checking
  if(!(strcmp(usePos,"center") || strcmp(usePos, "start") || strcmp(usePos, "end"))) {
    	PyErr_SetString(PyExc_NameError, "usePos must be either center, start or end");
	return NULL;
  }
  makeFractional = PyObject_IsTrue(fractional_input);
  //set list sizes to do calculations on (NEVER SET BEFORE parsing the tuple)
  
  std::vector<double> wiggle(tx_end - tx_start + 1, 0);    
  std::vector<int> pos_counts(tx_end - tx_start + 1, 0);
  std::vector<int> lengths;




  PyObject *iterator = PyObject_GetIter(reads);
  PyObject *item;
  PyObject *jxns;
  PyObject *allreads;

  jxns = PyDict_New();
  allreads = PySet_New(NULL);

  while (item = PyIter_Next(iterator)) {

    //skips reads on the wrong strand
    PyObject *is_reverse = PyObject_GetAttrString(item, "is_reverse");
    Py_INCREF(is_reverse); 
    
    if (is_reverse == Py_True && !strcmp(keepstrand, "+")) {
      continue;
    } else if(is_reverse == Py_False && !strcmp(keepstrand, "-")) {
      continue;
    }
    Py_DECREF(is_reverse);
    
    //sets read stard and stop location 
    PyObject *aligned_positions = PyObject_GetAttrString(item, "positions");
    Py_INCREF(aligned_positions);
    int positions_size = PyList_Size(aligned_positions);


      //long read_start = PyLong_AsLong(Pyread_statrtPyList_GetItem(aligned_positions, 0)); //possible bug here
    PyObject *Pyread_start = PyList_GetItem(aligned_positions, 0);
    PyObject *Pyread_stop = PyList_GetItem(aligned_positions, positions_size - 1);

    long read_start = PyLong_AsLong(Pyread_start); //possible bug here
    long read_stop = PyLong_AsLong(Pyread_stop); //get the last item in the list

    //skip if the read aligns past before the tx_start or after the tx_end
    if (read_start < tx_start || read_stop > tx_end) {
      continue;
    }

    //doing the hacky version of getting the read size, just getting all the aligned locations
    lengths.push_back(positions_size);
    
    //choose where to align single reads 
    if (!strcmp(usePos, "center")) {
      pos_counts[(((read_stop + read_start) / 2) - tx_start)] += 1; //assign pos counts 
    } else if (!strcmp(usePos, "start")) {
      if (strcmp(keepstrand, "+")) {
	pos_counts[read_start - tx_start]++;
      } else {
	pos_counts[read_stop - tx_start]++;
      } 
    } else if (!strcmp(usePos, "end")) {
      if (strcmp(keepstrand, "+")) {
	pos_counts[read_stop - tx_start]++;
      } else {
	pos_counts[read_start - tx_start]++;
      } 
    } else {
      return NULL;
    }

    //generate wiggle track from files


    PyObject *read_loc = PyTuple_New(2);
    PyTuple_SetItem(read_loc, 0, Pyread_start);
    PyTuple_SetItem(read_loc, 1, Pyread_stop);
    PySet_Add(allreads, read_loc);

    for(int i = 0; i < positions_size; i++) {
      //cur == pos == genome position (coordinate)
      PyObject *cur  = PyList_GetItem(aligned_positions, i);
      Py_INCREF(cur);
      if (cur == NULL) {
	return NULL;
      }
      long pos = PyLong_AsLong(cur);

      if (i+1 < positions_size){

	PyObject *nextcur  = PyList_GetItem(aligned_positions, (i+1));
	long nextpos = PyLong_AsLong(nextcur);
	if (nextpos > (pos +1)){

	  // the next position is > than this position + 1. this is a junction
	  
	  PyObject *jxn = PyTuple_New(2);
	  PyTuple_SetItem(jxn, 0, cur);
	  PyTuple_SetItem(jxn, 1, nextcur);

	  if (PyDict_Contains(jxns, jxn)){

	    PyObject *jxnCount = PyDict_GetItem(jxns,jxn);
	    long incremented = PyInt_AsLong(jxnCount);
	    incremented++;
	    
	    PyDict_SetItem(jxns, jxn, PyInt_FromLong(incremented));
	  } else{
	    PyDict_SetItem(jxns, jxn, PyInt_FromLong(1));
	  }
	}
      }

      Py_DECREF(cur);
      int wig_index = pos-tx_start;
      if (makeFractional) {
	wiggle[wig_index] += 1.0 / positions_size;
      } else {
	wiggle[wig_index]++;
      }
    }
    //item == read
    

    Py_DECREF(aligned_positions);
    Py_DECREF(item);
    





  }
  Py_DECREF(iterator); // iteration has ended, garbage collet it
  
  //all 3 items have been generated, convert them into PyLists and return them as a tuple 
  
  //convert lengths to pylist
  PyObject *retLengths = PyList_New(lengths.size());
  
  //transform vector to pylist
  for (int i = 0; i < lengths.size(); i++) {
    PyList_SetItem(retLengths, i, PyInt_FromLong(lengths[i]));
  }

  //convert single counts and wiggle as pylist
  PyObject *retWiggle = PyList_New(wiggle.size());
  PyObject *retPos_counts = PyList_New(pos_counts.size());
  
  //transform vector to pylist
  for (int i = 0; i < wiggle.size(); i++) {
    if(makeFractional) {
      PyList_SetItem(retWiggle, i, PyFloat_FromDouble(wiggle[i]));
    } else {
      PyList_SetItem(retWiggle, i, PyInt_FromLong(wiggle[i]));
    }
    
    PyList_SetItem(retPos_counts, i, PyInt_FromLong(pos_counts[i]));
  }
 
  //add 3 lists to tuple and return
  PyObject *return_tuple = PyTuple_New(5);
  PyTuple_SetItem(return_tuple, 0, retWiggle);
  PyTuple_SetItem(return_tuple, 1, jxns);
  PyTuple_SetItem(return_tuple, 2, retPos_counts);
  PyTuple_SetItem(return_tuple, 3, retLengths);
  PyTuple_SetItem(return_tuple, 4, allreads);

  
  return return_tuple;
}
static PyMethodDef peaks_methods[] = {
    {"shuffle",             peaks_shuffle,      METH_VARARGS,
     "Return the meaning of everything."},
    {"find_sections", peaks_find_sections, METH_VARARGS,
    "finds sections given a list and a margin"},
    {"readsToWiggle_pysam", peaks_readsToWiggle_pysam, METH_VARARGS,
    "converts pysam to a wiggle vector and some other stuff"},

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
