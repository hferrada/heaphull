/*
 * HeapCHull.h
 *
 *  Created on: 5 Apr 2017
 *      Author: ferrada
 */

#ifndef HEAPCHULL_H_
#define HEAPCHULL_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <random>
#include <time.h>
#include <dirent.h>
#include <bits/random.h>

using namespace std;

namespace hch {

#define WITHCH_XY 1				// 1: include arrays CH_X[0..nCH-1] and CH_Y[0..nCH-1] with the floating points of the hull.

#ifndef uchar
#define uchar unsigned char
#endif

//  compute and return the slope between the points (X[i],Y[i]) and (X[j],Y[j])
#define slope(i,j) ((Y[i]-Y[j]) / (X[i]-X[j]))

// return true if the points i,j and k are in counter clockwise order. Condition: X[i] >= X[j] >= X[k]
//#define isCCW_RL(i,j,k) ((Y[i]-Y[j])*(X[j]-X[k]) < (Y[j]-Y[k])*(X[i]-X[j]))

// return true if the points i,j and k are in counter clockwise order. Condition: X[i] <= X[j] <= X[k]
//#define isCCW_LR(i,j,k) ((Y[k]-Y[j])*(X[j]-X[i]) > (Y[j]-Y[i])*(X[k]-X[j]))

const uint W64 = 64;
const uint W64Min1 = 63;
const uint BW64 = 6;		// pow of two for W64
const uint WW64 = 128;
const ulong maskW63 = 0x8000000000000000;
const ulong maskLB = 0x000000000000FFFF;
const uint TWO32m1 = 4294967295;	// it is 2 at 32 minus 1

class HeapCHull {
private:
	void seacrhExtremePoints();

	void searchCHullPointsQ1();
	void searchCHullPointsQ2();
	void searchCHullPointsQ3();
	void searchCHullPointsQ4();

	void searchCHullPointsQ1_40();
	void searchCHullPointsQ2_40();
	void searchCHullPointsQ3_40();
	void searchCHullPointsQ4_40();

	// set the new item at top of the queue Q[1..k] (cells with 32 bits)
	void setTopMinQ(uint *Q, uint k, float *V);

	// set the new item at top of the queue Q[1..k] (cells with 32 bits)
	void setTopMaxQ(uint *Q, uint k, float *V);

	// The same than previous one but now handling queues with cells of 40 bits
	void setTopMinQ40(uint *Q, uchar* _Q, ulong k, float *V);
	void setTopMaxQ40(uint *Q, uchar* _Q, ulong k, float *V);

	// (cells with 32 bits)

	void createMaxHeapQ1();
	void createMaxHeapQ2();
	void createMinHeapQ3();
	void createMinHeapQ4();

	void createMinHeap(uint *Q, uint m, float *V);
	void createMaxHeap(uint *Q, uint m, float *V);
	void createMinHeapQ40(uint *Q, uchar* _Q, ulong m, float *V);
	void createMaxHeapQ40(uint *Q, uchar* _Q, ulong m, float *V);

	void createMinHeapCormen(uint *Q, uint m, float *V);
	void createMaxHeapCormen(uint *Q, uint m, float *V);

	// ####################################################
	// methods not used (queues with 64 bits)
	void setTopMinQ(ulong *Q, ulong k, float *V);
	void setTopMaxQ(ulong *Q, ulong k, float *V);
	void createMinHeap(ulong *Q, ulong m, float *V);
	void createMaxHeap(ulong *Q, ulong m, float *V);
	void seacrhExtremePoints3();
	void seacrhExtremePoints2();
	// ####################################################

	bool isQ40;					// we handle numbers with 40 bits !

	float *X, *Y;			    // array of original float points for X and Y coordinates
	ulong n, n1, n2, n3, n4;		// n1, n2, n3 and n4 are the points in Q1, Q2, Q3 and Q4 respectively
	uint lgn;					// Ceiling of logarithm n

	uint *Q1, *Q2, *Q3, *Q4;	// Queues for each quadrant
	uchar *_Q1, *_Q2, *_Q3, *_Q4;
	ulong ri1, up1; 			// the extreme points for Q1
	ulong up2, le2; 			// the extreme points for Q2
	ulong le3, lo3; 			// the extreme points for Q3
	ulong lo4, ri4; 			// the extreme points for Q4
	ulong c1, c2, c3, c4;		// these are the point more closet to each corner given by the extreme points

	float xri1, xup1, xup2, xle2, xle3, xlo3, xlo4, xri4, xc1, xc2, xc3, xc4;
	float yri1, yup1, yup2, yle2, yle3, ylo3, ylo4, yri4, yc1, yc2, yc3, yc4;

	float m1, m2, m3, m4;		// line-slope for the lines: L1:line(ri1, c1), L2:line(up2, c2), L3:line(le3, c3) and L4:line(lo4, c4) respectively
	float m1b, m2b, m3b, m4b;	// additional line-slope for the lines: L1:line(c1, up1), L2:line(c2, le2), L3:line(c3, lo3) and L4:line(c4, ri4) respectively
	float mh;					// horizontal slope for Lh:line(le2,ri1)

public:
	static bool TRACE;			// true: print all details for console
	static bool TEST;
	static bool LeftRight;		// true: make the final convex hull from the leftmost point to the rightmost one in clockwise order
								// false: from right to left also in clockwise order

	HeapCHull(ulong nPoins, float *xArr, float *yArr);
	virtual ~HeapCHull();

	// create the arrays X[0...nCH] and Y[0..nCH] with the final points of the convex hull form the original arrays  xArr[] and yArr[]
	// It must be called in your code after computing the convex hull. This is equivalent to put the flag 'WITHCH_XY' in 1
	void genXY(float **X, float **Y, float *xArr, float *yArr);

	void inPointsInQ40();
	void inPointsInQ32();

	// set the number x as a bitstring sequence in *A. In the range of bits [ini, .. ini+len-1] of *A. Here x has len bits
	void setNum64(ulong *A, ulong ini, uint len, ulong x);

	// return (in a unsigned long integer) the number in A from bits of position 'ini' to 'ini+len-1'
	ulong getNum64(ulong *A, ulong ini, uint len);

	ulong *CH;					// this is the final Convex Hull in counterclockwise order
	ulong nCH;					// nCH is the number of points in the final Convex Hull CH[1..nCH]
	ulong nIgn;					// Discarded points
	ulong nQs;					// points in queues
	ulong sizeCHull;

	float *CH_X, *CH_Y;			// these will be created only if WITHCH_XY is 1

	// methods and vars for test...
	bool isInCHull(ulong i);
	void testCHull(ulong rep);
	bool isCCW_LR(ulong i, ulong j, ulong k);
	bool isCCW_RL(ulong i, ulong j, ulong k);
	ulong changeP;				// used for test

	// save file chull...
	void createFileChull(char filePoints[300]);
};

} /* namespace hch */

#endif /* HEAPCHULL_H_ */
