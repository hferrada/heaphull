#include "HeapCHull.h"
using namespace hch;

bool TRACE = 0;		// true: print all details for console
bool TEST = 0;		// true: apply exhaustive test

// Structure with all globals parameters program
typedef struct {
	float sizeSt;			// size in bytes for the data structure

	ulong n;				// length (numbers of symbols) of sequence
	uint lgn;				// Ceiling of logarithm n

	float *X, *Y;			// array of float points for X and Y coordinates
	HeapCHull *myHull;
} ParProg;

int main(int argc, char *argv[]) {
	ParProg *par = new ParProg();

	if(argc != 2){
		cout << "¡ERROR WITH PARAMETERS! " << endl;
		cout << "./heaphull <n>" << endl;
		exit(1);
	}
	HeapCHull::TRACE = TRACE;
	HeapCHull::TEST = TEST;

	par->n = atol(argv[1]);
	par->sizeSt = 0;

	cout << "heap hull's Parameters..." << endl;
	cout << "[1] n: " << par->n << endl;

	float k = 2.0*par->n*sizeof(float);
	cout << endl <<" Allocating " << k << " Bytes = " << k/(1024.0*1024.0) << " MiB for X[] and Y[] arrays" << endl;

	// random inputs arrays...
	par->X = new float[par->n];
	par->Y = new float[par->n];
	for(ulong i=0; i<par->n; i++){
		par->X[i] = (float) rand()/RAND_MAX;
		par->Y[i] = (float) rand()/RAND_MAX;
	}
	if (TRACE){
		cout << "Input arrays (first 100 points)..." << endl;
		ulong m = par->n;
		if (m > 100) m = 100;

		for(ulong i=0; i<m; i++)
			cout << "(" << par->X[i] << ", " << par->Y[i] << ") " << endl;

		cout << endl;
	}

	cout << " Computing CH... " << endl;
	clock_t t = clock();
	par->myHull = new HeapCHull(par->n, par->X, par->Y);
	t = clock() - t;

	cout << endl << "*** Convex Hull CH[1.." << par->myHull->nCH << "]" << endl;
	if (TRACE){
		for(ulong i=0; i<par->myHull->nCH; i++)
			cout << par->myHull->CH[i] << " " << par->X[par->myHull->CH[i]] << " " << par->Y[par->myHull->CH[i]] << endl;
		cout << endl;
	}

	cout << " CPU construction time....  : " << (t*1000.0)/CLOCKS_PER_SEC << " ms" << endl;;
	cout << " Structure' size            : " << par->myHull->sizeCHull << " Bytes = "<< (float)(par->myHull->sizeCHull)/1024.0 << " KB = " << ((float)(par->myHull->sizeCHull)/1024.0)/1024.0 << " MiB aprox." << endl;
	cout << " nCH (Convex Hull's points) : " << par->myHull->nCH << endl;
	cout << " nIgn (pruned points)       : " << par->myHull->nIgn << " = " << 100.0*float(par->myHull->nIgn)/float(par->n) << " % of n" << endl;

	delete [] par->X;
	delete [] par->Y;

	cout << "#########################################" << endl;
}
