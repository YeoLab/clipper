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

//length: int, effective length of gene
//iterations: int, number of iterations to run
//reads: list, list of read lengths
//cutoff: double, B-H FDR cutoff
//stats: boolean, true prints running time stats
extern "C" PyObject shuffle(PyObject *length, PyObject *iterations, PyObject *reads, PyObject *cutoff, PyObject *stats)
{
  long L = PyInt_AsLong(length);
  long r = PyInt_AsLong(iterations);
  bool T = PyInt_AsLong(stats);
  double alpha = PyFloat_AsDouble(cutoff);
  
  //if (L ==0 || r ==0 || readfile ==0) usage(argv[0]);
  //Need to fix this to throw error if L r or the file is null or none
  
  int redone = 0;
  int OBS_CUTOFF[1000]; // more than 1000 observed heights will kill this...
  int size_readslist = 0;
  
  //Initalize all observed hights to zero (maybe need to convert this to pyobject)
  for (int cut = 0; cut < 1000; cut++){
    OBS_CUTOFF[cut] = 0;
  }
  
  int GENE[L]; //will store the height at every position of the gene
  int HEIGHT[L]; // will have the number of times you see every height across the gene  HEIGHT[3] = 10 means that there are 10 peaks with a height of 3
  ifstream indata; // indata is like cin
  ifstream indata2; // indata is like cin

  std::vector<int> readlist;    
  //if (strcasecmp(file, "stdin")==0){
  //  char line[256];

  //  while (!cin.eof()){
  //    cin.getline(line, 256);
  //    readlist.push_back(atoi(line));
  //  }

  //}
  //else{
  //  indata.open(file); // opens the file containing the list of read lengths
  //  if(!indata) { // file couldn't be opened
  //    cout << "Error: file:" << file << "could not be opened" << endl;
  //    exit(1);
  //  }
  //  char line[256];
  //  while ( !indata.eof() ) { // keep reading until end-of-file
  //    indata.getline(line,256);
  //    readlist.push_back(atoi(line));
  //  }
    //}
  //size_readslist = readlist.size();
  
  //Up to here we are just loading the read list into an array for use
  //If I make this an external module I won't have to do that.
  //int reads[size_readslist];
  //for( int i=0; i<(size_readslist);i++){
  //  reads[i]= readlist[i];
  //}
  
  //Reads should take care of everything up to here
  readlist.clear();
  srand(time(NULL)); // seed random-number generator with the current time
  for (int iteration = 1; iteration <= r; iteration++){
    //initialize GENE and HEIGHT
    for (int j = 0; j < L; j++){
      GENE[j] = 0;
    }
    for (int j = 0; j < L; j++){
      HEIGHT[j] = 0;
    }
    //end initialize

    //for each read assign it randomly
    for (int i=0; i<PyList_Size(reads); i++){
      long len = PyInt_AsLong(PyList_GetItem(reads, i));

      //error checking for obviously invalid reads
      if (len > L){
	cout << "read" << i << " is longer than gene (" << len << "), exiting" << endl;
	exit(8);
      }
      int ran;
      
      //Pick a random location witin the gene that the read can map to
      do{           // pick a random number within the length of the gene
	ran = rand()%L;
      }while((ran+len) > L); // pick again if the read will fall outside of the gene
      
      //Increment the coverage along the length of that read
      for(int position=ran; position < (ran+len); position++)	{
	GENE[position]++;
      }
    }
    int total_num_peaks = 0;
    int max_height = 0;
    for (int j = 0; j < L; j++) {
      HEIGHT[GENE[j]]++; // simulated height distribution
      if (GENE[j] > 0){
	total_num_peaks++;
      }
      if (GENE[j] > max_height){
	max_height = GENE[j];
      }
    }
    //cout << " calculating significance " << endl;
    //here is where you compute what height is significant.

    double PVAL[max_height];
    std::vector<int> sig_heights;
    PVAL[0] = 1; // set pvalue for height 0 = to 1;
    for (int height =1; height < max_height; height++){
      PVAL[height] =1;
      int bigger_peaks = 0; //number of peaks that are heigher than height
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
      //      cout << height << "\t" << bigger_peaks << "\t" << total_num_peaks << "\t" << a<< endl;
    }

    //benjamani-hochberg FDR (correction of p-values)
    //    cout << " B-H correction " << endl;

    double N_heights = sig_heights.size();
    double R = 0.0;
    double correctedP[sizeof(PVAL) / sizeof(PVAL[0])];
    std::vector<int> corr_sig_heights;
    //test
    for (int ele = sig_heights.size()-1; ele>=0; ele--){
      correctedP[sig_heights[ele]] = PVAL[sig_heights[ele]] * N_heights/(N_heights - R);
      R += 1.0;
      //      cout << "ele\t" << ele << "\tPVAL:\t" << PVAL[sig_heights[ele]] << "\tN_heights:\t" << N_heights << "\tsig_h:\t" << sig_heights[ele] << "\tcP\t" << correctedP[sig_heights[ele]] << endl; // print out what you wee here.
      if (correctedP[sig_heights[ele]] < alpha){
	corr_sig_heights.push_back(sig_heights[ele]); // add the current height to corrected significant heights list.
      }
    }
    if (corr_sig_heights.size() ==0 ){ // no heights are significant
      //      cout << "redoing iteration # " << iteration << endl;
      iteration--;
      redone++;
      if (redone > (10*r)){
	//	fprintf (stderr, "skipping this gene, there aren't enough reads\n");
	break;
      }
      continue;
    }

    int height_cutoff = corr_sig_heights[(corr_sig_heights.size() - 1)]; // height cutoff is the smallest height that is still significant.
    OBS_CUTOFF[height_cutoff]++;
    //    cout << "height cutoff is " << height_cutoff << endl;
  }

  PyObject* returnList = PyList_New(1000);
  for (int cut =0; cut < 1000; cut++){
    PyList_Append(returnList, PyInt_FromLong(OBS_CUTOFF[cut]));
    //if (OBS_CUTOFF[cut] > 0){
      //cout << cut << "\t" << OBS_CUTOFF[cut] << endl;
      //}
  }

  if (redone > 0){
    //    fprintf(stderr,"Perhaps you should make alpha larger?  I had to re-iterate an extra %i times.\n", redone);
  }
  //calculate runtime
  if (T == 1)    {
    double CPS = CLOCKS_PER_SEC;
    double C = clock();
    double runtime = C/CPS;
    cout << "runtime: " << runtime << " seconds" << endl;
    cout << "number of reads: " << size_readslist << endl;
    cout << "iterations: " << r << endl;
    cout << "gene length: " << L << endl;
  }
  //


  return *returnList;
}

