/*
 * HeapCHull.cpp
 *
 *  Created on: 5 Apr 2017
 *      Author: ferrada
 */

#include "HeapCHull.h"

namespace hch {
bool HeapCHull::TRACE = false;
bool HeapCHull::TEST = false;
bool HeapCHull::LeftRight = false;	// make the final convex hull from the leftmost point to the rightmost one in clockwise order

HeapCHull::HeapCHull(ulong nPoins, float *xArr, float *yArr) {
	n = nPoins;
	lgn = 1 + log(n)/log(2);
	X = xArr;
	Y = yArr;
	seacrhExtremePoints();

	ulong len = 1+n;
	Q1 = new uint[len];
	Q2 = new uint[len];
	Q3 = new uint[len];
	Q4 = new uint[len];
	sizeCHull = 4*len*sizeof(uint);

	if (len > TWO32m1){
		// we work with cells (of 32 + 8) bits in the queues, that is Qi[k] || _Qi[k]
		_Q1 = new uchar[len];
		_Q2 = new uchar[len];
		_Q3 = new uchar[len];
		_Q4 = new uchar[len];
		sizeCHull += 4*len*sizeof(uchar);
		isQ40 = true;
	}else
		// we work with cells of 32 bits in the queues
		isQ40 = false;

	//if (TRACE) cout << " ** Fixed size for all Queues : " << sizeCHull << " Bytes" << endl;

	// computing main slopes...
	n1=n2=n3=n4=1;
	if (xri1!=xup1){
		if (xc1>xup1 && yc1>yri1){
			m1 = (yri1-yc1)/(xri1-xc1); 		//slope(ri1, c1);
			m1b = (yri1-yup1)/(xri1-xup1);		//slope(ri1, up1);
			if (m1 < m1b){
				m1b = (yc1-yup1)/(xc1-xup1); 	//slope(c1, up1);
				if (isQ40){
					Q1[1]=c1>>8;
					_Q1[1]=c1&maskLB;
				}else Q1[1]=c1;
				n1=2;
			}else{
				m1 = m1b;
				m1b = 0;
				xc1 = yc1 = -1*FLT_MAX;
			}
		}else{
			m1 = (yri1-yup1)/(xri1-xup1);		//slope(ri1, up1);
			m1b = 0;
			xc1 = yc1 = -1*FLT_MAX;
		}
	}else{
		m1 = m1b = 0;
		xc1 = yc1 = -1*FLT_MAX;
	}

	if (xup2!=xle2){
		if (xc2<xup2 && yc2>yle2){
			m2 = (yup2-yc2)/(xup2-xc2);			//slope(up2, c2);
			m2b = (yup2-yle2)/(xup2-xle2);		//slope(up2, le2);
			if (m2 < m2b){
				m2b = (yc2-yle2)/(xc2-xle2);	//slope(c2, le2);
				if (isQ40){
					Q2[1]=c2>>8;
					_Q2[1]=c2&maskLB;
				}else Q2[1]=c2;
				n2=2;
			}else{
				m2 = m2b;
				m2b = 0;
				xc2 = yc2 = -1*FLT_MAX;
			}
		}else{
			m2 = (yup2-yle2)/(xup2-xle2);		//slope(up2, le2);
			m2b = 0;
			xc2 = yc2 = -1*FLT_MAX;
		}
	}else{
		m2 = m2b = 0;
		xc2 = yc2 = -1*FLT_MAX;
	}

	if (xle3!=xlo3){
		if (xc3<xlo3 && yc3<yle3){
			m3 = (yle3-yc3)/(xle3-xc3);			//slope(le3, c3);
			m3b = (yle3-ylo3)/(xle3-xlo3);		//slope(le3, lo3);
			if (m3 < m3b){
				m3b = (ylo3-yc3)/(xlo3-xc3);	//slope(c3, lo3);
				if (isQ40){
					Q3[1]=c3>>8;
					_Q3[1]=c3&maskLB;
				}else Q3[1]=c3;
				n3=2;
			}else{
				m3 = m3b;
				m3b = 0;
				xc3 = yc3 = FLT_MAX;
			}
		}else{
			m3 = (yle3-ylo3)/(xle3-xlo3);		//slope(le3, lo3);
			m3b = 0;
			xc3 = yc3 = FLT_MAX;
		}
	}else{
		m3 = m3b = 0;
		xc3 = yc3 = FLT_MAX;
	}

	if (xlo4!=xri4){
		if (xc4>xlo4 && yc4<yri4){
			m4 = (ylo4-yc4)/(xlo4-xc4);			//slope(lo4, c4);
			m4b = (ylo4-yri4)/(xlo4-xri4);		//slope(lo4, ri4);
			if (m4 < m4b){
				m4b = (yc4-yri4)/(xc4-xri4);	//slope(c4, ri4);
				if (isQ40){
					Q4[1]=c4>>8;
					_Q4[1]=c4&maskLB;
				}else Q4[1]=c4;
				n4=2;
			}else{
				m4 = m4b;
				m4b = 0;
				xc4 = yc4 = FLT_MAX;
			}
		}else{
			m4 = (ylo4-yri4)/(xlo4-xri4);		//slope(lo4, ri4);
			m4b = 0;
			xc4 = yc4 = FLT_MAX;
		}
	}else{
		m4 = m4b = 0;
		xc4 = yc4 = FLT_MAX;
	}

	if (xri1!=xle2)
		mh = (yri1-yle2)/(xri1-xle2);	 		//slope(le2, ri1);
	else
		mh = 0;

	if (isQ40){
		inPointsInQ40();

		nQs = n1+n2+n3+n4;
		nIgn = n-nQs;

		// to make a heap for each Q...
		createMaxHeapQ40(Q1, _Q1, n1, X);
		createMaxHeapQ40(Q2, _Q2, n2, Y);
		createMinHeapQ40(Q3, _Q3, n3, X);
		createMinHeapQ40(Q4, _Q4, n4, Y);
	}else{
		inPointsInQ32();

		// to make a heap for each Q...
		createMaxHeapCormen(Q1, n1, X);
		createMaxHeapCormen(Q2, n2, Y);
		createMinHeapCormen(Q3, n3, X);
		createMinHeapCormen(Q4, n4, Y);

		nQs = n1+n2+n3+n4;
		nIgn = n-nQs;
	}

	if(TRACE){
		ulong i;

		cout << endl;
		cout << "ri1 = " << ri1 << " = (" << xri1 << ", " << yri1 << ")" << endl;
		cout << "c1 = " << c1 << " = (" << xc1 << ", " << yc1 << ")" << endl;
		cout << "up1 = " << up1 << " = (" << xup1 << ", " << yup1 << ")" << endl;
		cout << "up2 = " << up2 << " = (" << xup2 << ", " << yup2 << ")" << endl;
		cout << "c2 = " << c2 << " = (" << xc2 << ", " << yc2 << ")" << endl;
		cout << "le2 = " << le2 << " = (" << xle2 << ", " << yle2 << ")" << endl;
		cout << "le3 = " << le3 << " = (" << xle3 << ", " << yle3 << ")" << endl;
		cout << "c3 = " << c3 << " = (" << xc3 << ", " << yc3 << ")" << endl;
		cout << "lo3 = " << lo3 << " = (" << xlo3 << ", " << ylo3 << ")" << endl;
		cout << "lo4 = " << lo4 << " = (" << xlo4 << ", " << ylo4 << ")" << endl;
		cout << "c4 = " << c4 << " = (" << xc4 << ", " << yc4 << ")" << endl;
		cout << "ri4 = " << ri4 << " = (" << xri4 << ", " << yri4 << ")" << endl << endl;

		if (xc1!=xri1){
			cout << "L1:  (xri1,yri1), (xc1,yc1) = (" << xri1 << ", " << yri1 << "), (" << xc1 << ", " << yc1 << "): m1= " << m1 << endl;
			cout << "L1b: (xc1,yc1), (xup1,yup1) = (" << xc1 << ", " << yc1 << "), (" << xup1 << ", " << yup1 << "): m1= " << m1 << endl;
		}else
			cout << "L1: (xri1,yri1), (xup1,yup1) = (" << xri1 << ", " << yri1 << "), (" << xup1 << ", " << yup1 << "): m1= " << m1 << endl;

		if (xc2!=xup2){
			cout << "L2:  (xup2,yup2), (xc2,yc2) = (" << xup2 << ", " << yup2 << "), (" << xc2 << ", " << yc2 << "): m2= " << m2 << endl;
			cout << "L2b: (xc2,yc2), (xle2,yle2) = (" << xc2 << ", " << yc2 << "), (" << xle2 << ", " << yle2 << "): m2= " << m2 << endl;
		}else
			cout << "L2: (xup2,yup2), (xle2,yle2) = (" << xup2 << ", " << yup2 << "), (" << xle2 << ", " << yle2 << "): m2= " << m2 << endl;

		if (xc3!=xle3){
			cout << "L3:  (xle3,yle3), (xc3,yc3) = (" << xle3 << ", " << yle3 << "), (" << xc3 << ", " << yc3 << "): m3= " << m3 << endl;
			cout << "L3b: (xc3,yc3), (xlo3,ylo3) = (" << xc3 << ", " << yc3 << "), (" << xlo3 << ", " << ylo3 << "): m3= " << m3 << endl;
		}else
			cout << "L3: (xle3,yle3), (xlo3,ylo3) = (" << xle3 << ", " << yle3 << "), (" << xlo3 << ", " << ylo3 << "): m3= " << m3 << endl;

		if (xc4!=xri4){
			cout << "L4:  (xlo4,ylo4), (xc4,yc4) = (" << xlo4 << ", " << ylo4 << "), (" << xc4 << ", " << yc4 << "): m4= " << m4 << endl;
			cout << "L4b: (xc4,yc4), (xri4,yri4) = (" << xc4 << ", " << yc4 << "), (" << xri4 << ", " << yri4 << "): m4= " << m4 << endl;
		}else
			cout << "L4: (xlo4,ylo4), (xri4,yri4) = (" << xlo4 << ", " << ylo4 << "), (" << xri4 << ", " << yri4 << "): m4= " << m4 << endl;

		if (xle2!=xri1)
			cout << "Lh: (xle2,yle2), (xri1,yri1) = (" << xle2 << ", " << yle2 << "), (" << xri1 << ", " << yri1 << "): mh= " << mh << endl;
		else
			cout << "mh = 0" << endl;

		if (!isQ40){
			if (n1>1){
				cout << endl << "Max X[Q1[1..." << n1-1 << "]]:" << endl;
				for(i=1; i<n1; i++)
					//cout << X[Q1[i]] << " ";
					cout << Q1[i] << " ";
				cout << endl;
			}
			if (n2>1){
				cout << "Max Y[Q2[1..." << n2-1 << "]]:" << endl;
				for(i=1; i<n2; i++)
					//cout << Y[Q2[i]] << " ";
					cout << Q2[i] << " ";
				cout << endl;
			}
			if (n3>1){
				cout << "Min X[Q3[1..." << n3-1 << "]]:" << endl;
				for(i=1; i<n3; i++)
					//cout << X[Q3[i]] << " ";
					cout << Q3[i] << " ";
				cout << endl;
			}
			if (n4>1){
				cout << "Min Y[Q4[1..." << n4-1 << "]]:" << endl;
				for(i=1; i<n4; i++)
					//cout << X[Q4[i]] << " ";
					cout << Q4[i] << " ";
				cout << endl;
			}
		}else{
			ulong qi;
			if (n1>1){
				cout << endl << "Max X[Q1[1..." << n1-1 << "]]:" << endl;
				for(i=1; i<n1; i++){
					qi = ((Q1[i]<<8)+_Q1[i]);
					cout << X[qi] << " ";
				}
				cout << endl;
			}
			if (n2>1){
				cout << "Max Y[Q2[1..." << n2-1 << "]]:" << endl;
				for(i=1; i<n2; i++){
					qi = ((Q2[i]<<8)+_Q2[i]);
					cout << Y[qi] << " ";
				}
				cout << endl;
			}
			if (n3>1){
				cout << "Min X[Q3[1..." << n3-1 << "]]:" << endl;
				for(i=1; i<n3; i++){
					qi = ((Q3[i]<<8)+_Q3[i]);
					cout << X[qi] << " ";
				}
				cout << endl;
			}
			if (n4>1){
				cout << "Min Y[Q4[1..." << n4-1 << "]]:" << endl;
				for(i=1; i<n4; i++){
					qi = ((Q4[i]<<8)+_Q4[i]);
					cout << Y[qi] << " ";
				}
				cout << endl;
			}
		}
	}

	if(TEST){
		cout << "Testing queues..." << endl;
		if (!isQ40){
			ulong i,j,lim = n1>>1;
			for(i=1; i<lim; i++){
				j = i<<1;
				if (X[Q1[i]] < X[Q1[j]]){
					cout << "ERROR, for i="<<i<<", X[Q1[i]]=X["<<Q1[i]<<"] = "<<X[Q1[i]]<<" < X[Q1[2i]]=X["<<Q1[j]<<"] = "<<X[Q1[j]]<< endl;
					exit(0);
				}
				if (j+1 < n1 && X[Q1[i]] < X[Q1[j+1]]){
					cout << "ERROR, for i="<<i<<", X[Q1[i]]=X["<<Q1[i]<<"] = "<<X[Q1[i]]<<" < X[Q1[2i+1]]=X["<<Q1[j+1]<<"] = "<<X[Q1[j+1]]<< endl;
					exit(0);
				}
			}

			lim = n2>>1;
			for(i=1; i<lim; i++){
				j = i<<1;
				if (Y[Q2[i]] < Y[Q2[j]]){
					cout << "ERROR, for i="<<i<<", Y[Q2[i]]=Y["<<Q2[i]<<"] = "<<Y[Q2[i]]<<" < Y[Q2[2i]]=Y["<<Q2[j]<<"] = "<<Y[Q2[j]]<< endl;
					exit(0);
				}
				if (j+1 < n2 && Y[Q2[i]] < Y[Q2[j+1]]){
					cout << "ERROR, for i="<<i<<", Y[Q2[i]]=Y["<<Q2[i]<<"] = "<<Y[Q2[i]]<<" < Y[Q2[2i+1]]=Y["<<Q2[j+1]<<"] = "<<Y[Q2[j+1]]<< endl;
					exit(0);
				}
			}

			lim = n3>>1;
			for(i=1; i<lim; i++){
				j = i<<1;
				if (X[Q3[i]] > X[Q3[j]]){
					cout << "ERROR, for i="<<i<<", X[Q3[i]]=X["<<Q3[i]<<"] = "<<X[Q3[i]]<<" > X[Q3[2i]]=X["<<Q3[j]<<"] = "<<X[Q3[j]]<< endl;
					exit(0);
				}
				if (j+1 < n3 && X[Q3[i]] > X[Q3[j+1]]){
					cout << "ERROR, for i="<<i<<", X[Q3[i]]=X["<<Q3[i]<<"] = "<<X[Q3[i]]<<" > X[Q3[2i+1]]=X["<<Q3[j+1]<<"] = "<<X[Q3[j+1]]<< endl;
					exit(0);
				}
			}

			lim = n4>>1;
			for(i=1; i<lim; i++){
				j = i<<1;
				if (Y[Q4[i]] > Y[Q4[j]]){
					cout << "ERROR, for i="<<i<<", Y[Q4[i]]=Y["<<Q4[i]<<"] = "<<Y[Q4[i]]<<" > Y[Q4[2i]]=Y["<<Q4[j]<<"] = "<<Y[Q4[j]]<< endl;
					exit(0);
				}
				if (j+1 < n4 && Y[Q4[i]] > Y[Q4[j+1]]){
					cout << "ERROR, for i="<<i<<", Y[Q4[i]]=Y["<<Q4[i]<<"] = "<<Y[Q4[i]]<<" > Y[Q4[2i+1]]=Y["<<Q4[j+1]<<"] = "<<Y[Q4[j+1]]<< endl;
					exit(0);
				}
			}
		}else{
			ulong i,qi,qj,j,lim = n1>>1;
			for(i=1; i<lim; i++){
				j = i<<1;
				qi = ((Q1[i]<<8)+_Q1[i]);
				qj = ((Q1[j]<<8)+_Q1[j]);
				if (X[qi] < X[qj]){
					cout << "ERROR, for i="<<i<<", X[Q1[i]]=X["<<qi<<"] = "<<X[qi]<<" < X[Q1[2i]]=X["<<qj<<"] = "<<X[qj]<< endl;
					exit(0);
				}
				qj = ((Q1[j+1]<<8)+_Q1[j+1]);
				if (j+1 < n1 && X[qi] < X[qj]){
					cout << "ERROR, for i="<<i<<", X[Q1[i]]=X["<<qi<<"] = "<<X[qi]<<" < X[Q1[2i+1]]=X["<<qj<<"] = "<<X[qj]<< endl;
					exit(0);
				}
			}

			lim = n2>>1;
			for(i=1; i<lim; i++){
				j = i<<1;
				qi = ((Q2[i]<<8)+_Q2[i]);
				qj = ((Q2[j]<<8)+_Q2[j]);
				if (Y[qi] < Y[qj]){
					cout << "ERROR, for i="<<i<<", Y[Q2[i]]=Y["<<qi<<"] = "<<Y[qi]<<" < Y[Q2[2i]]=Y["<<qj<<"] = "<<Y[qj]<< endl;
					exit(0);
				}
				qj = ((Q2[j+1]<<8)+_Q2[j+1]);
				if (j+1 < n2 && Y[qi] < Y[qj]){
					cout << "ERROR, for i="<<i<<", Y[Q2[i]]=Y["<<qi<<"] = "<<Y[qi]<<" < Y[Q2[2i+1]]=Y["<<qj<<"] = "<<Y[qj]<< endl;
					exit(0);
				}
			}

			lim = n3>>1;
			for(i=1; i<lim; i++){
				j = i<<1;
				qi = ((Q3[i]<<8)+_Q3[i]);
				qj = ((Q3[j]<<8)+_Q3[j]);
				if (X[qi] > X[qj]){
					cout << "ERROR, for i="<<i<<", X[Q3[i]]=X["<<qi<<"] = "<<X[qi]<<" > X[Q3[2i]]=X["<<qj<<"] = "<<X[qj]<< endl;
					exit(0);
				}
				qj = ((Q3[j+1]<<8)+_Q3[j+1]);
				if (j+1 < n3 && X[qi] > X[qj]){
					cout << "ERROR, for i="<<i<<", X[Q3[i]]=X["<<qi<<"] = "<<X[qi]<<" > X[Q3[2i+1]]=X["<<qj<<"] = "<<X[qj]<< endl;
					exit(0);
				}
			}
			lim = n4>>1;
			for(i=1; i<lim; i++){
				j = i<<1;
				qi = ((Q4[i]<<8)+_Q4[i]);
				qj = ((Q4[j]<<8)+_Q4[j]);
				if (Y[qi] > Y[qj]){
					cout << "ERROR, for i="<<i<<", Y[Q4[i]]=Y["<<qi<<"] = "<<Y[qi]<<" > Y[Q4[2i]]=Y["<<qj<<"] = "<<Y[qj]<< endl;
					exit(0);
				}
				qj = ((Q4[j+1]<<8)+_Q4[j+1]);
				if (j+1 < n4 && Y[qi] > Y[qj]){
					cout << "ERROR, for i="<<i<<", Y[Q4[i]]=Y["<<qi<<"] = "<<Y[qi]<<" > Y[Q4[2i+1]]=Y["<<qj<<"] = "<<Y[qj]<< endl;
					exit(0);
				}
			}
		}
		//cout << " Test for Queues ok !" << endl;
	}

	// To make the final convex hull from queues (the extreme points are not into the queues)
	nQs+=8;
	len = 1 + nQs;
	CH = new ulong[len];
	sizeCHull += len*sizeof(uint);			// ** SIZE FOR CONVEX HULL WITH MAXIMUM LENGHT

	if (LeftRight==false){
		// make the CH from the rightmost to the leftmost point
		CH[0] = ri1;
		nCH = 1;
		if (ri1 != up1){
			if (n1>1){
				if (isQ40){
					searchCHullPointsQ1_40();
					//delete [] _Q1;
				}
				else searchCHullPointsQ1();
			}
			CH[nCH] = up1;	// the last point in Convex Hull for II-quadrant
			nCH++;
		}
		//delete [] Q1; // be deleted in the destroyer method !!

		if (up1 != up2){
			nIgn--;
			CH[nCH] = up2;
			nCH++;
		}
		changeP = nCH-1; // this is only for test !!
		if (up2 != le2){
			if (n2>1){
				if (isQ40){
					searchCHullPointsQ2_40();
					//delete [] _Q2;
				}
				else searchCHullPointsQ2();
			}
			CH[nCH] = le2;	// the last point in Convex Hull for II-quadrant
			changeP = nCH;
			nCH++;
		}
		//delete [] Q2;

		if (le2 != le3){
			nIgn--;
			CH[nCH] = le3;
			nCH++;
		}
		if (le3 != lo3){
			if (n3>1){
				if (isQ40){
					searchCHullPointsQ3_40();
					//delete [] _Q3;
				}
				else searchCHullPointsQ3();
			}
			if (lo3 != ri1){
				CH[nCH] = lo3;	// the last point in Convex Hull for III-quadrant
				nCH++;
			}
		}
		//delete [] Q3;

		if (lo3 != lo4){
			nIgn--;
			CH[nCH] = lo4;
			nCH++;
		}
		if (lo4 != ri4){
			if (n4>1){
				if (isQ40){
					searchCHullPointsQ4_40();
					//delete [] _Q4;
				}
				else searchCHullPointsQ4();
			}
			if (ri4 != ri1){
				nIgn--;
				CH[nCH] = ri4;	// the last point in Convex Hull for IV-quadrant
				nCH++;
			}
		}
		//delete [] Q4;
	}else{
		// make the CH from the leftmost to the rightmost point
		CH[0] = le3;
		nCH = 1;
		if (le3 != lo3){
			if (n3>1){
				if (isQ40){
					searchCHullPointsQ3_40();
					//delete [] _Q3;
				}
				else searchCHullPointsQ3();
			}
			CH[nCH] = lo3;	// the last point in Convex Hull for III-quadrant
			nCH++;
		}
		//delete [] Q3;

		if (lo3 != lo4){
			CH[nCH] = lo4;
			nCH++;
			nIgn--;
		}
		if (lo4 != ri4){
			if (n4>1){
				if (isQ40){
					searchCHullPointsQ4_40();
					//delete [] _Q4;
				}
				else searchCHullPointsQ4();
			}
			CH[nCH] = ri4;	// the last point in Convex Hull for IV-quadrant
			nCH++;
		}
		//delete [] Q4;
		changeP = nCH-1; // this is only for test !!
		if (ri4 != ri1){
			CH[nCH] = ri1;
			nCH++;
			nIgn--;
			changeP = nCH;
		}
		if (ri1 != up1){
			if (n1>1){
				if (isQ40){
					searchCHullPointsQ1_40();
					//delete [] _Q1;
				}
				else searchCHullPointsQ1();
			}
			CH[nCH] = up1;	// the last point in Convex Hull for II-quadrant
			nCH++;
		}
		//delete [] Q1;

		if (up1 != up2){
			CH[nCH] = up2;
			nCH++;
			nIgn--;
		}
		if (up2 != le2){
			if (n2>1){
				if (isQ40){
					searchCHullPointsQ2_40();
					//delete [] _Q2;
				}
				else searchCHullPointsQ2();
			}
			if (le2 != le3){
				CH[nCH] = le2;	// the last point in Convex Hull for II-quadrant
				nCH++;
				nIgn--;
			}
		}
		//delete [] Q2;
	}
	if (nIgn>n) nIgn=0; // only for statics !!

	if (WITHCH_XY){
		CH_X = new float[nCH];
		CH_Y = new float[nCH];
		for (len=0; len<nCH; len++){
			CH_X[len] = xArr[CH[len]];
			CH_Y[len] = yArr[CH[len]];
		}
	}
}

// create the float arrays X[] and Y[] with the floating convex hull's points
void HeapCHull::genXY(float **X, float **Y, float *xArr, float *yArr){
	ulong i;
	float *XX, *YY;

	*X = XX = new float[nCH];
	*Y = YY = new float[nCH];

	if (WITHCH_XY){
		for (i=0; i<nCH; i++){
			XX[i] = CH_X[i];
			YY[i] = CH_Y[i];
		}
	}else{
		for (i=0; i<nCH; i++){
			XX[i] = xArr[CH[i]];
			YY[i] = yArr[CH[i]];
		}
	}
}


// search the extreme points...
void HeapCHull::seacrhExtremePoints(){
	float xi,yi;
	bool seei=true;

	c1 = c2 = c3 = c4 = n+1;
	xc1 = xc4 = yc1 = yc2 = -1*FLT_MAX;
	xc2 = xc3 = yc3 = yc4 = FLT_MAX;
	xri1 = xup1 = xup2 = xle2 = xle3 = xlo3 = xlo4 = xri4 = X[0];
	yri1 = yup1 = yup2 = yle2 = yle3 = ylo3 = ylo4 = yri4 = Y[0];
	ri1 = up1 = up2 = le2 = le3 = lo3 = lo4 = ri4 = 0;
	for(ulong i=1; i<n; i++, seei=true){
		xi=X[i];yi=Y[i];
		if(xi < xle2){
			if (xc2-yc2 > xle2-yle2){
				c2 = le2;
				xc2=xle2;
				yc2=yle2;
			}
			if (xc3+yc3 > xle3+yle3){
				c3 = le3;
				xc3=xle3;
				yc3=yle3;
			}
			xle2 = xle3 = xi;
			yle2 = yle3 = yi;
			le2 = le3 = i;
			seei=false;
		}else{
			if(xi == xle2){ // no actualizar c2 o c3 pq tiene el mismo x q xi
				if(yi > yle2){
					le2 = i;
					xle2 = xi;
					yle2 = yi;
					seei=false;
				}else{
					if(yi < yle3){
						le3 = i;
						xle3 = xi;
						yle3 = yi;
						seei=false;
					}
				}
			}else{
				if(xi > xri1){
					if (xc1+yc1 < xri1+yri1){
						c1 = ri1;
						xc1=xri1;
						yc1=yri1;
					}
					if (yc4-xc4 > yri4-xri4){
						c4 = ri4;
						xc4=xri4;
						yc4=yri4;
					}
					ri1 = ri4 = i;
					xri1 = xri4 = xi;
					yri1 = yri4 = yi;
					seei=false;
				}else{
					if(xi == xri1){
						if(yi > yri1){
							ri1 = i;
							xri1 = xi;
							yri1 = yi;
							seei=false;
						}else{
							if(yi < yri4){
								ri4 = i;
								xri4 = xi;
								yri4 = yi;
								seei=false;
							}
						}
					}
				}
			}
		}


		if(yi < ylo3){
			if (xc3+yc3 > xlo3+ylo3){
				xc3=xlo3;
				yc3=ylo3;
				c3 = lo3;
			}
			if (yc4-xc4 > ylo4-xlo4){
				xc4=xlo4;
				yc4=ylo4;
				c4 = lo4;
			}

			xlo3 = xlo4 = xi;
			ylo3 = ylo4 = yi;
			lo3 = lo4 = i;
			seei=false;
		}else{
			if(yi == ylo3){
				if(xi < xlo3){
					xlo3 = xi;
					ylo3 = yi;
					lo3 = i;
					seei=false;
				}else{
					if(xi > xlo4){
						xlo4 = xi;
						ylo4 = yi;
						lo4 = i;
						seei=false;
					}
				}
			}else{
				if(yi > yup2){
					if (xc1+yc1 < xup1+yup1){
						xc1=xup1;
						yc1=yup1;
						c1 = up1;
					}
					if (xc2-yc2 > xup2-yup2){
						xc2=xup2;
						yc2=yup2;
						c2 = up2;
					}

					xup2 = xup1 = xi;
					yup2 = yup1 = yi;
					up1 = up2 = i;
					seei=false;
				}else{
					if(yi == yup2){
						if(xi < xup2){
							xup2 = xi;
							yup2 = yi;
							up2 = i;
							seei=false;
						}else{
							if(xi > xup1){
								xup1 = xi;
								yup1 = yi;
								up1 = i;
								seei=false;
							}
						}
					}
				}
			}

		}

		if (seei){
			if (xc1+yc1 < xi+yi){
				c1=i;
				xc1=xi;
				yc1=yi;
			}else{
				if (xc2-yc2 > xi-yi){
					c2=i;
					xc2=xi;
					yc2=yi;
				}else{
					if (xc3+yc3 > xi+yi){
						c3=i;
						xc3=xi;
						yc3=yi;
					}else{
						if (yc4-xc4 > yi-xi){
							c4=i;
							xc4=xi;
							yc4=yi;
						}
					}
				}
			}
		}
	}
}


// looking for convex hull's points in the queue Q4. Quadrant IV: [lo4 --> ri4]
void HeapCHull::searchCHullPointsQ4(){
	ulong i, prv, ini=nCH;
	float mi, mp, mf, xi, yi, auxx, ylast=ylo4, xlast=xlo4;

	mf = (ylo4-yri4)/(xlo4-xri4);							//slope(le2,up2);
	for(n4--; n4; n4--){
		i=Q4[1];
		setTopMinQ(Q4, n4, Y);
		xi = X[i];
		yi = Y[i];

		while (n4>1 && yi==Y[Q4[1]]){
			prv = Q4[1];
			auxx = X[prv];
			if (auxx > xi){
				xi = auxx;
				i = prv;
			}
			n4--;
			setTopMinQ(Q4, n4, Y);
		}

		if(xi>xlast){
			mi = (ylast-yi)/(xlast-xi);						// slope(last,i);
			if (mi<mf){
				// x must be in the current Convex Hull !!
				if (nCH>ini){
					prv = CH[nCH-2];
					mp = (ylast-Y[prv])/(xlast-X[prv]);				//slope(prv,last);
					mi = (yi-Y[prv])/(xi-X[prv]);					//slope(prv,i);
					while(mi<=mp){	// We exclude points in the same line
						nCH--;
						if (nCH>ini){
							xlast=X[prv];
							ylast=Y[prv];
							prv = CH[nCH-2];
							mp = (ylast-Y[prv])/(xlast-X[prv]);		//slope(prv,last);
							mi = (yi-Y[prv])/(xi-X[prv]);			//slope(prv,i);
						}else break;
					}
				}
				xlast=xi;
				ylast=yi;
				mf = (yri4-ylast)/(xri4-xlast);				//slope(ri4,last);
				CH[nCH] = i;
				nCH++;
			}
		}
		//cout << "CH[]: ";for(ulong j=0;j<nCH;cout<<CH[j]<<" ",j++);cout<<endl;
	}
}
void HeapCHull::searchCHullPointsQ4_40(){
	ulong i, prv, ini=nCH;
	float mi, mp, mf, xi, yi, auxx, ylast=ylo4, xlast=xlo4;

	mf = (ylo4-yri4)/(xlo4-xri4);							//slope(le2,up2);
	for(n4--; n4; n4--){
		i = ((Q4[1]<<8)+_Q4[1]);
		setTopMinQ40(Q4, _Q4, n4, Y);
		xi = X[i];
		yi = Y[i];

		if (n4>1){
			prv = ((Q4[1]<<8)+_Q4[1]);
			while (yi==Y[prv]){
				auxx = X[prv];
				if (auxx > xi){
					xi = auxx;
					i = prv;
				}

				n4--;
				setTopMinQ40(Q4, _Q4, n4, Y);
				if (n4==1)
					break;
				prv = ((Q4[1]<<8)+_Q4[1]);
			}
		}

		if(xi>xlast){
			mi = (ylast-yi)/(xlast-xi);						// slope(last,i);
			if (mi<mf){
				// x must be in the current Convex Hull !!
				if (nCH>ini){
					prv = CH[nCH-2];
					mp = (ylast-Y[prv])/(xlast-X[prv]);				//slope(prv,last);
					mi = (yi-Y[prv])/(xi-X[prv]);					//slope(prv,i);
					while(mi<=mp){	// We exclude points in the same line
						nCH--;
						if (nCH>ini){
							xlast=X[prv];
							ylast=Y[prv];
							prv = CH[nCH-2];
							mp = (ylast-Y[prv])/(xlast-X[prv]);		//slope(prv,last);
							mi = (yi-Y[prv])/(xi-X[prv]);			//slope(prv,i);
						}else break;
					}
				}
				xlast=xi;
				ylast=yi;
				mf = (yri4-ylast)/(xri4-xlast);				//slope(ri4,last);
				CH[nCH] = i;
				nCH++;
			}
		}
		//cout << "CH[]: ";for(ulong j=0;j<nCH;cout<<CH[j]<<" ",j++);cout<<endl;
	}
}

// looking for convex hull's points in the queue Q3. Quadrant III: [le3 --> lo3]
void HeapCHull::searchCHullPointsQ3(){
	ulong i, prv, ini=nCH;
	float mi, mp, mf, xi, yi, auxy, ylast=yle3, xlast=xle3;

	mf = (yle3-ylo3)/(xle3-xlo3);							//slope(le2,up2);
	for(n3--; n3; n3--){
		i=Q3[1];
		xi = X[i];
		yi = Y[i];
		setTopMinQ(Q3, n3, X);

		while (n3>1 && xi==X[Q3[1]]){
			prv = Q3[1];
			auxy = Y[prv];
			if (auxy < yi){
				yi = auxy;
				i = prv;
			}
			n3--;
			setTopMinQ(Q3, n3, X);
		}

		if(yi<ylast){
			mi = (ylast-yi)/(xlast-xi);						// slope(last,i);
			if (mi<mf){
				// x must be in the current Convex Hull !!
				if (nCH>ini){
					prv = CH[nCH-2];
					mp = (ylast-Y[prv])/(xlast-X[prv]);				//slope(prv,last);
					mi = (yi-Y[prv])/(xi-X[prv]);					//slope(prv,i);
					while(mi<=mp){// we exclude points in the same line
						nCH--;
						if (nCH>ini){
							xlast=X[prv];
							ylast=Y[prv];
							prv = CH[nCH-2];
							mp = (ylast-Y[prv])/(xlast-X[prv]);		//slope(prv,last);
							mi = (yi-Y[prv])/(xi-X[prv]);			//slope(prv,i);
						}else break;
					}
				}
				xlast=xi;
				ylast=yi;
				mf = (ylo3-ylast)/(xlo3-xlast);				//slope(lo3,last);
				CH[nCH] = i;
				nCH++;
			}
		}
		//cout << "CH[]: ";for(ulong j=0;j<nCH;cout<<CH[j]<<" ",j++);cout<<endl;
	}
}
void HeapCHull::searchCHullPointsQ3_40(){
	ulong i, prv, ini=nCH;
	float mi, mp, mf, xi, yi, auxy, ylast=yle3, xlast=xle3;

	mf = (yle3-ylo3)/(xle3-xlo3);							//slope(le2,up2);
	for(n3--; n3; n3--){
		i = ((Q3[1]<<8)+_Q3[1]);
		setTopMinQ40(Q3, _Q3, n3, X);
		xi = X[i];
		yi = Y[i];

		if (n3>1){
			prv = ((Q3[1]<<8)+_Q3[1]);
			while (xi==X[prv]){
				auxy = Y[prv];
				if (auxy < yi){
					yi = auxy;
					i = prv;
				}

				n3--;
				setTopMinQ40(Q3, _Q3, n3, X);
				if (n3==1)
					break;
				prv = ((Q3[1]<<8)+_Q3[1]);
			}
		}

		if(yi<ylast){
			mi = (ylast-yi)/(xlast-xi);						// slope(last,i);
			if (mi<mf){
				// x must be in the current Convex Hull !!
				if (nCH>ini){
					prv = CH[nCH-2];
					mp = (ylast-Y[prv])/(xlast-X[prv]);				//slope(prv,last);
					mi = (yi-Y[prv])/(xi-X[prv]);					//slope(prv,i);
					while(mi<=mp){// we exclude points in the same line
						nCH--;
						if (nCH>ini){
							xlast=X[prv];
							ylast=Y[prv];
							prv = CH[nCH-2];
							mp = (ylast-Y[prv])/(xlast-X[prv]);		//slope(prv,last);
							mi = (yi-Y[prv])/(xi-X[prv]);			//slope(prv,i);
						}else break;
					}
				}
				xlast=xi;
				ylast=yi;
				mf = (ylo3-ylast)/(xlo3-xlast);				//slope(lo3,last);
				CH[nCH] = i;
				nCH++;
			}
		}
		//cout << "CH[]: ";for(ulong j=0;j<nCH;cout<<CH[j]<<" ",j++);cout<<endl;
	}
}

// looking for convex hull's points in the queue Q2. Quadrant II: [up2 --> le2]
void HeapCHull::searchCHullPointsQ2(){
	ulong i, prv, ini=nCH;
	float mi, mp, mf, xi, yi, auxx, ylast=yup2, xlast=xup2;

	mf = (yle2-yup2)/(xle2-xup2);							//slope(le2,up2);
	for(n2--; n2; n2--){
		i=Q2[1];
		setTopMaxQ(Q2, n2, Y);
		xi = X[i];
		yi = Y[i];

		while (n2>1 && yi==Y[Q2[1]]){
			prv = Q2[1];
			auxx = X[prv];
			if (auxx < xi){
				xi = auxx;
				i = prv;
			}
			n2--;
			setTopMaxQ(Q2, n2, Y);
		}

		if(xi<xlast){
			mi = (ylast-yi)/(xlast-xi);						// slope(last,i);
			if (mi<mf){
				// x must be in the current Convex Hull !!
				if (nCH>ini){
					prv = CH[nCH-2];
					mp = (ylast-Y[prv])/(xlast-X[prv]);				//slope(prv,last);
					mi = (yi-Y[prv])/(xi-X[prv]);					//slope(prv,i);
					while(mi<=mp){// we exclude points in the same line (with '<' these will be include)
						nCH--;
						if (nCH>ini){
							xlast=X[prv];
							ylast=Y[prv];
							prv = CH[nCH-2];
							mp = (ylast-Y[prv])/(xlast-X[prv]);		//slope(prv,last);
							mi = (yi-Y[prv])/(xi-X[prv]);			//slope(prv,i);
						}else break;
					}
				}
				xlast=xi;
				ylast=yi;
				mf = (yle2-ylast)/(xle2-xlast);				//slope(le2,last);
				CH[nCH] = i;
				nCH++;
			}

		}
		//cout << "CH[]: ";for(ulong j=0;j<nCH;cout<<CH[j]<<" ",j++);cout<<endl;
	}
}
void HeapCHull::searchCHullPointsQ2_40(){
	ulong i, prv, ini=nCH;
	float mi, mp, mf, xi, yi, auxx, ylast=yup2, xlast=xup2;

	mf = (yle2-yup2)/(xle2-xup2);							//slope(le2,up2);
	for(n2--; n2; n2--){
		i = ((Q2[1]<<8)+_Q2[1]);
		setTopMaxQ40(Q2, _Q2, n2, Y);
		xi = X[i];
		yi = Y[i];

		if (n2>1){
			prv = ((Q2[1]<<8)+_Q2[1]);
			while (yi==Y[prv]){
				auxx = X[prv];
				if (auxx < xi){
					xi = auxx;
					i = prv;
				}
				n2--;
				setTopMaxQ40(Q2, _Q2, n2, Y);
				if (n2==1)
					break;
				prv = ((Q2[1]<<8)+_Q2[1]);
			}
		}

		if(xi<xlast){
			mi = (ylast-yi)/(xlast-xi);						// slope(last,i);
			if (mi<mf){
				// x must be in the current Convex Hull !!
				if (nCH>ini){
					prv = CH[nCH-2];
					mp = (ylast-Y[prv])/(xlast-X[prv]);				//slope(prv,last);
					mi = (yi-Y[prv])/(xi-X[prv]);					//slope(prv,i);
					while(mi<=mp){// we exclude points in the same line (with '<' these will be include)
						nCH--;
						if (nCH>ini){
							xlast=X[prv];
							ylast=Y[prv];
							prv = CH[nCH-2];
							mp = (ylast-Y[prv])/(xlast-X[prv]);		//slope(prv,last);
							mi = (yi-Y[prv])/(xi-X[prv]);			//slope(prv,i);
						}else break;
					}
				}
				xlast=xi;
				ylast=yi;
				mf = (yle2-ylast)/(xle2-xlast);				//slope(le2,last);
				CH[nCH] = i;
				nCH++;
			}

		}
		//cout << "CH[]: ";for(ulong j=0;j<nCH;cout<<CH[j]<<" ",j++);cout<<endl;
	}
}

// looking for convex hull's points in the queue Q1. Quadrant I: [ri1 --> up1]
void HeapCHull::searchCHullPointsQ1(){
	ulong i, prv, ini=nCH;
	float mi, mp, mf, xi, yi, auxy, ylast=yri1, xlast=xri1;

	mf = (yup1-yri1)/(xup1-xri1);							//slope(up1,last);
	for(n1--; n1; n1--){
		i=Q1[1];
		setTopMaxQ(Q1, n1, X);
		xi = X[i];
		yi = Y[i];

		while (n1>1 && xi==X[Q1[1]]){
			prv = Q1[1];
			auxy = Y[prv];
			if (auxy > yi){
				yi = auxy;
				i = prv;
			}
			n1--;
			setTopMaxQ(Q1, n1, X);
		}

		if(yi>ylast){
			mi = (ylast-yi)/(xlast-xi);						// slope(last,i);
			if (mi<mf){
				// x is above the line (up1)...(last) --> x must be in the current Convex Hull !!
				if (nCH>ini){
					prv = CH[nCH-2];
					mp = (ylast-Y[prv])/(xlast-X[prv]);				//slope(prv,last);
					mi = (yi-Y[prv])/(xi-X[prv]);					//slope(prv,i);
					while(mi<=mp){									// if (mi == mp) --> {i, last, prv} are co-linear points
						nCH--;
						if (nCH>ini){
							xlast=X[prv];
							ylast=Y[prv];
							prv = CH[nCH-2];
							mp = (ylast-Y[prv])/(xlast-X[prv]);		//slope(prv,last);
							mi = (yi-Y[prv])/(xi-X[prv]);			//slope(prv,i);
						}else break;
					}
				}
				xlast=xi;
				ylast=yi;
				mf = (yup1-ylast)/(xup1-xlast);				//slope(up1,last);
				CH[nCH] = i;
				nCH++;
			}
		}
		//cout << "CH[]: ";for(ulong j=0;j<nCH;cout<<CH[j]<<" ",j++);cout<<endl;
	}
}
void HeapCHull::searchCHullPointsQ1_40(){
	ulong i, prv, ini=nCH;
	float mi, mp, mf, xi, yi, auxy, ylast=yri1, xlast=xri1;

	mf = (yup1-yri1)/(xup1-xri1);							//slope(up1,last);
	for(n1--; n1; n1--){
		i = ((Q1[1]<<8)+_Q1[1]);
		setTopMaxQ40(Q1, _Q1, n1, X);
		xi = X[i];
		yi = Y[i];

		if (n1>1){
			prv = ((Q1[1]<<8)+_Q1[1]);
			while (xi==X[prv]){
				auxy = Y[prv];
				if (auxy > yi){
					yi = auxy;
					i = prv;
				}
				n1--;
				setTopMaxQ40(Q1, _Q1, n1, X);
				if (n1==1)
					break;
				prv = ((Q1[1]<<8)+_Q1[1]);
			}
		}

		if(yi>ylast){
			mi = (ylast-yi)/(xlast-xi);						// slope(last,i);
			if (mi<mf){
				// x is above the line (up1)...(last) --> x must be in the current Convex Hull !!
				if (nCH>ini){
					prv = CH[nCH-2];
					mp = (ylast-Y[prv])/(xlast-X[prv]);				//slope(prv,last);
					mi = (yi-Y[prv])/(xi-X[prv]);					//slope(prv,i);
					while(mi<=mp){
						nCH--;
						if (nCH>ini){
							xlast=X[prv];
							ylast=Y[prv];
							prv = CH[nCH-2];
							mp = (ylast-Y[prv])/(xlast-X[prv]);		//slope(prv,last);
							mi = (yi-Y[prv])/(xi-X[prv]);			//slope(prv,i);
						}else break;
					}
				}
				xlast=xi;
				ylast=yi;
				mf = (yup1-ylast)/(xup1-xlast);				//slope(up1,last);
				CH[nCH] = i;
				nCH++;
			}
		}
		//cout << "CH[]: ";for(ulong j=0;j<nCH;cout<<CH[j]<<" ",j++);cout<<endl;
	}
}
// set the new item at top of the queue Q[1..k]
void HeapCHull::setTopMinQ(uint *Q, uint k, float *V){
	uint m=2, i=1, qm, qk=Q[k];
	float num = V[qk];

	while(m<k){
		qm = Q[m];
		if (m+1 < k && V[qm] > V[Q[m+1]]){
			m++;
			qm = Q[m];
		}

		if(num > V[qm]){
			Q[i] = qm;
			i = m;
			m <<= 1;
		}else
			break;
	}
	Q[i] = qk;
}

// set the new item at top of the queue Q[1..k]
void HeapCHull::setTopMaxQ(uint *Q, uint k, float *V){
	uint m=2, i=1, qm, qk=Q[k];
	float num = V[qk];

	while(m<k){
		qm = Q[m];
		if (m+1 < k && V[qm] < V[Q[m+1]]){
			m++;
			qm = Q[m];
		}

		if(num < V[qm]){
			Q[i] = qm;
			i = m;
			m <<= 1;
		}else
			break;
	}
	Q[i] = qk;
}

// set the new item at top of the queue Q[1..k]
void HeapCHull::setTopMinQ40(uint *Q, uchar* _Q, ulong k, float *V){
	ulong m=2, i=1, qm, qmm, qk=((Q[k]<<8)+_Q[k]);
	float num = V[qk];

	while(m<k){
		qm = ((Q[m]<<8)+_Q[m]);
		if (m+1 < k){
			qmm = ((Q[m+1]<<8)+_Q[m+1]);
			if (V[qm] > V[qmm]){
				m++;
				qm=qmm;
			}
		}

		if(num > V[qm]){
			Q[i] = Q[m];
			_Q[i] = _Q[m];
			i = m;
			m <<= 1;
		}else
			break;
	}
	Q[i] = Q[k];
	_Q[i] = _Q[k];
}

// set the new item at top of the queue Q[1..k]
void HeapCHull::setTopMaxQ40(uint *Q, uchar* _Q, ulong k, float *V){
	ulong m=2, i=1, qm, qmm, qk=((Q[k]<<8)+_Q[k]);
	float num = V[qk];

	while(m<k){
		qm = ((Q[m]<<8)+_Q[m]);
		if (m+1 < k){
			qmm = ((Q[m+1]<<8)+_Q[m+1]);
			if (V[qm] < V[qmm]){
				m++;
				qm=qmm;
			}
		}

		if(num < V[qm]){
			Q[i] = Q[m];
			_Q[i] = _Q[m];
			i = m;
			m <<= 1;
		}else
			break;
	}
	Q[i] = Q[k];
	_Q[i] = _Q[k];
}


// to create the minimum priority queue Q[1..m]
void HeapCHull::createMinHeap(uint *Q, uint m, float *V){
	uint i,j,k,t;
	float v;

	for(i=2; i<m; i++){
		j=i;
		k=Q[j];
		v=V[k];
		t=j>>1;
		while(t && V[Q[t]]>v){
			Q[j] = Q[t];
			j=t;
			t>>=1;
		}
		Q[j]=k;
	}
}

// to make maximum priority queue Q1[1..n1]
void HeapCHull::createMaxHeapQ1(){
	uint i,j,k,t,m=0;
	float x;

	for(i=2; i<n1; i++){
		k=Q1[i];
		x=X[k];
		if (x == X[Q1[1]]){		// we select only the highmost y-coordinate when its x-values are equal!
			if (Y[k] > Y[Q1[1]])
				Q1[1]=k;
			m++;
		}else{
			j=i-m;
			t=j>>1;
			while(t && X[Q1[t]]<x){
				Q1[j] = Q1[t];
				j=t;
				t>>=1;
			}
			Q1[j]=k;
		}
	}
	n1-=m;
}

// to make maximum priority queue Q1[1..n1]
void HeapCHull::createMaxHeapQ2(){
	uint i,j,k,t,m=0;
	float y;

	for(i=2; i<n2; i++){
		k=Q2[i];
		y=Y[k];
		if (y == Y[Q2[1]]){		// we select only the highmost y-coordinate when its x-values are equal!
			if (X[k] < X[Q2[1]])
				Q2[1]=k;
			m++;
		}else{
			j=i-m;
			t=j>>1;
			while(t && Y[Q2[t]]<y){
				Q2[j] = Q2[t];
				j=t;
				t>>=1;
			}
			Q2[j]=k;
		}
	}
	n2-=m;
}


// to create the minimum priority queue Q3[1..n3]
void HeapCHull::createMinHeapQ3(){
	uint i,j,k,t,m=0;
	float x;

	for(i=2; i<n3; i++){
		k=Q3[i];
		j=i-m;
		Q3[j]=k;
		x=X[k];
		if (x == X[Q3[1]]){		// we select only the lowest y-coordinate when its x-values are equal!
			if (Y[k] < Y[Q3[1]])
				Q3[1]=k;
			m++;
		}else{
			t=j>>1;
			while(t && X[Q3[t]]>=x){
				if (X[Q3[t]]==x){
					if (Y[k] < Y[Q3[t]]){
						Q3[j] = Q3[t];
						j=t;
						t>>=1;
					}
					m++;
				}else{
					Q3[j] = Q3[t];
					j=t;
					t>>=1;
				}
			}
			Q3[j]=k;
		}
	}
	n3-=m;


	for(i=2; i<n3; i++){
		x=X[Q3[i]];
		if (x == X[Q3[i-1]]){
			cout << "ERRRRRRRRR i = " << i << endl;
			exit(0);
		}

	}
}

// to create the minimum priority queue Q4[1..n4]
void HeapCHull::createMinHeapQ4(){
	uint i,j,k,t,m=0;
	float y;

	for(i=2; i<n4; i++){
		k=Q4[i];
		y=Y[k];
		if (y == Y[Q4[1]]){		// we select only the lowest y-coordinate when its x-values are equal!
			if (X[k] > X[Q4[1]])
				Q4[1]=k;
			m++;
		}else{
			j=i-m;
			t=j>>1;
			while(t && Y[Q4[t]]>y){
				Q4[j] = Q4[t];
				j=t;
				t>>=1;
			}
			Q4[j]=k;
		}
	}
	n4-=m;
}
// to make maximum priority queue Q[1..m]
void HeapCHull::createMaxHeap(uint *Q, uint m, float *V){
	uint i,j,k,t;
	float v;

	for(i=2; i<m; i++){
		j=i;
		k=Q[j];
		v=V[k];
		t=j>>1;
		while(t && V[Q[t]]<v){
			Q[j] = Q[t];
			j=t;
			t>>=1;
		}
		Q[j]=k;
	}
}

// to create the minimum priority queue Q[1..m]
void HeapCHull::createMinHeapQ40(uint *Q, uchar* _Q, ulong m, float *V){
	ulong i,j,k,t,qt;
	float v;

	for(i=2; i<m; i++){
		j=i;
		k=((Q[j]<<8)+_Q[j]);
		v=V[k];
		t=j>>1;
		qt = ((Q[t]<<8)+_Q[t]);
		while(t && V[qt]>v){
			Q[j] = Q[t];
			_Q[j] = _Q[t];
			j=t;
			t>>=1;
			qt = ((Q[t]<<8)+_Q[t]);
		}
		Q[j]=k>>8;
		_Q[j] = k&maskLB;
	}
}

// to make maximum priority queue Q[1..m]
void HeapCHull::createMaxHeapQ40(uint *Q, uchar* _Q, ulong m, float *V){
	uint i,j,k,t,qt;
	float v;

	for(i=2; i<m; i++){
		j=i;
		k=((Q[j]<<8)+_Q[j]);
		v=V[k];
		t=j>>1;
		qt = ((Q[t]<<8)+_Q[t]);
		while(t && V[qt]<v){
			Q[j] = Q[t];
			_Q[j] = _Q[t];
			j=t;
			t>>=1;
			qt = ((Q[t]<<8)+_Q[t]);
		}
		Q[j]=k>>8;
		_Q[j] = k&maskLB;
	}
}

// to create the minimum priority queue Q[1..m]
void HeapCHull::createMinHeapCormen(uint *Q, uint m, float *V){
	uint i,j,k,t,mid=m/2;
	float v;

	for(i=mid; i; i--){
		j=i;
		k=Q[j];
		v=V[k];
		t=j<<1;
		if (t+1<m && V[Q[t+1]] < V[Q[t]])
			t++;

		// bajamos v hasta su posición correcta, intercambiandolo con el menor de sus hijos
		while(t<m && v>V[Q[t]]){
			Q[j] = Q[t];
			j=t;
			t<<=1;
			if (t+1<m && V[Q[t+1]] < V[Q[t]])
				t++;
		}
		Q[j]=k;
	}
}

// to make maximum priority queue Q[1..m]
void HeapCHull::createMaxHeapCormen(uint *Q, uint m, float *V){
	uint i,j,k,t,mid=m/2;
	float v;

	for(i=mid; i; i--){
		j=i;
		k=Q[j];
		v=V[k];
		t=j<<1;
		if (t+1<m && V[Q[t+1]] > V[Q[t]])
			t++;

		// bajamos v hasta su posicion correcta, intercambiandolo con el mayor de sus hijos
		while(t<m && v<V[Q[t]]){
			Q[j] = Q[t];
			j=t;
			t<<=1;
			if (t+1<m && V[Q[t+1]] > V[Q[t]])
				t++;
		}
		Q[j]=k;
	}
}
// Q1 and Q2 are maximum priority queues. Q3 and Q4 are minimum priority queues
void HeapCHull::inPointsInQ40(){
	float m, x, y;

	for(ulong i=0; i<n; i++){
		x=X[i];y=Y[i];
		if(x!=xri1){
			m=(y-yri1)/(x-xri1);									// slope(i, ri1);
			if(m<mh){												// p_i is above horizontal line (le2...ri1)
				if(x>xup1){
					if ((m<m1)||(x<xc1 && (y-yc1)/(x-xc1)<m1b)){
						Q1[n1] = i>>8;
						_Q1[n1] = i&maskLB;
						n1++;
					}
				}else{
					if(x<xup2){
						m=(y-yup2)/(x-xup2);						// slope(i, up2);
						if ((m<m2)||(x<xc2 && (y-yc2)/(x-xc2)<m2b)){
							Q2[n2] = i>>8;
							_Q2[n2] = i&maskLB;
							n2++;
						}
					}
				}
			}else{
				if(x<xlo3){
					m=(y-yle3)/(x-xle3);							//slope(i, le3);
					if ((m<m3)||(x>xc3 && (y-yc3)/(x-xc3)<m3b)){	//slope(i, c3);
						Q3[n3] = i>>8;
						_Q3[n3] = i&maskLB;
						n3++;
					}
				}else{
					if(x>xlo4){
						m=(y-ylo4)/(x-xlo4);						// slope(i, lo4);
						if ((m<m4)||(x>xc4 && (y-yc4)/(x-xc4)<m4b)){// slope(i, c4);
							Q4[n4] = i>>8;
							_Q4[n4] = i&maskLB;
							n4++;
						}
					}
				}
			}
		}
	}

}

// Q1 and Q2 are maximum priority queues. Q3 and Q4 are minimum priority queues
void HeapCHull::inPointsInQ32(){
	float m, x, y;

	for(ulong i=0; i<n; i++){
		x=X[i];y=Y[i];
		if(x!=xri1){
			m=(y-yri1)/(x-xri1);									// slope(i, ri1);
			if(m<mh){												// p_i is above horizontal line (le2...ri1)
				if(x>xup1){
					if ((m<m1)||(x<xc1 && (y-yc1)/(x-xc1)<m1b)){
						Q1[n1] = i;
						n1++;
					}
				}else{
					if(x<xup2){
						m=(y-yup2)/(x-xup2);						// slope(i, up2);
						if ((m<m2)||(x<xc2 && (y-yc2)/(x-xc2)<m2b)){
							Q2[n2] = i;
							n2++;
						}
					}
				}
			}else{
				if(x<xlo3){
					m=(y-yle3)/(x-xle3);							//slope(i, le3);
					if ((m<m3)||(x>xc3 && (y-yc3)/(x-xc3)<m3b)){	//slope(i, c3);
						Q3[n3] = i;
						n3++;
					}
				}else{
					if(x>xlo4){
						m=(y-ylo4)/(x-xlo4);						// slope(i, lo4);
						if ((m<m4)||(x>xc4 && (y-yc4)/(x-xc4)<m4b)){// slope(i, c4);
							Q4[n4] = i;
							n4++;
						}
					}
				}
			}
		}
	}
}


//====================================================================================================================

bool HeapCHull::isInCHull(ulong i){
	if (i == ri1 || i == ri4 || i == le2 || i == le3)
		return true;

	ulong ini, end, m;
	float mi = slope(ri1, i);
	float xi = X[i];

	if (mi < mh){
		// we search in the upper points of the CH
		if (LeftRight){
			ini=changeP;
			if(le2 != le3)
				ini++;
			end=nCH-1;
		}else{
			ini=0;
			end=changeP;
		}

		m = ini + (end-ini)/2;
		while(X[CH[m]] != xi){
			if (X[CH[m]] < xi)
				end = m-1;
			else
				ini = m+1;
			if (ini <= end)
				m = ini + (end-ini)/2;
			else
				break;
		}
	}else{
		if (LeftRight){
			ini=0;
			end=changeP;
		}else{
			ini=changeP;
			if(le2 != le3)
				ini++;
			end=nCH-1;
		}
		m = ini + (end-ini)/2;
		while(X[CH[m]] != xi){
			if (X[CH[m]] < xi)
				ini = m+1;
			else
				end = m-1;
			if (ini <= end)
				m = ini + (end-ini)/2;
			else
				break;
		}
	}

	if(X[CH[m]] == xi && Y[CH[m]] == Y[i])
		return true;

	return false;
}

// return true if the points i,j and k are in counter clockwise order. Condition: X[i] <= X[j] <= X[k]
//#define isCCW_LR(i,j,k) ((Y[k]-Y[j])*(X[j]-X[i]) > (Y[j]-Y[i])*(X[k]-X[j]))
bool HeapCHull::isCCW_LR(ulong i, ulong j, ulong k){
	if ((Y[k]-Y[j])*(X[j]-X[i]) > (Y[j]-Y[i])*(X[k]-X[j]))
		return true;
	return false;
}

// return true if the points i,j and k are in counter clockwise order. Condition: X[i] >= X[j] >= X[k]
//#define isCCW_RL(i,j,k) ((Y[i]-Y[j])*(X[j]-X[k]) < (Y[j]-Y[k])*(X[i]-X[j]))
bool HeapCHull::isCCW_RL(ulong i, ulong j, ulong k){
	if ((Y[i]-Y[j])*(X[j]-X[k]) < (Y[j]-Y[k])*(X[i]-X[j]))
		return true;
	return false;
}

// check the correctness of the convex hull
void HeapCHull::testCHull(ulong rep){
	ulong i, j, k, x, sp, ep;
	float mi, mh2, xi, mij, mik;
	float epsERR = 0.05;

	// I.- test if CH[0..nCh-1] is a convex polygon
	//cout << "  I .- Testing if CH[0..nCh-1] is a convex polygon..." << endl;
	if (LeftRight){
		for(i=0, j=1, k=2; k<=changeP; i++, j++, k++){
			if(!isCCW_LR(CH[i],CH[j],CH[k])){
				mij = slope(CH[i], CH[j]);
				mik = slope(CH[i], CH[k]);
				if (mij == mik) continue;
				if ((mij<mik && (mik-mij)>epsERR) || (mik<mij && (mij-mik)>epsERR)){
					cout << "slope(i, j) = " << mij << endl;
					cout << "slope(i, k) = " << mik << endl;
					cout << "Error in rep(0) = " << rep << ". The convex hull is not a convex polygon !" << endl;
					cout << "Lower lines !isCCW_LR(CH[i],CH[j],CH[k]) to i=" << i << ", j=" << j << ", k=" << k << endl;
					cout << "(CH[i],CH[j],CH[k]) = (" << CH[i] << ", " << CH[j] << "," << CH[k] << ")" << endl;
					cout << "(" << X[CH[i]] << ", " << Y[CH[i]] << ")" << endl;
					cout << "(" << X[CH[j]] << ", " << Y[CH[j]] << ")" << endl;
					cout << "(" << X[CH[k]] << ", " << Y[CH[k]] << ")" << endl;
					exit(1);
				}
			}
		}

		for(i=changeP, j=changeP+1, k=changeP+2; k<nCH; i++, j++, k++){
			if(!isCCW_RL(CH[i],CH[j],CH[k])){
				mij = slope(CH[i], CH[j]);
				mik = slope(CH[i], CH[k]);
				if (mij == mik) continue;
				if ((mij<mik && (mik-mij)>epsERR) || (mik<mij && (mij-mik)>epsERR)){
					cout << "slope(i, j) = " << mij << endl;
					cout << "slope(i, k) = " << mik << endl;
					cout << "Error in rep(0) = " << rep << ". The convex hull is not a convex polygon !" << endl;
					cout << "Upper lines !isCCW_RL(CH[i],CH[j],CH[k]) to i=" << i << ", j=" << j << ", k=" << k << endl;
					cout << "(CH[i],CH[j],CH[k]) = (" << CH[i] << ", " << CH[j] << "," << CH[k] << ")" << endl;
					cout << "(" << X[CH[i]] << ", " << Y[CH[i]] << ")" << endl;
					cout << "(" << X[CH[j]] << ", " << Y[CH[j]] << ")" << endl;
					cout << "(" << X[CH[k]] << ", " << Y[CH[k]] << ")" << endl;
					exit(1);
				}
			}
		}
	}else{
		for(i=0, j=1, k=2; k<=changeP; i++, j++, k++){
			if(!isCCW_RL(CH[i],CH[j],CH[k])){
				mij = slope(CH[i], CH[j]);
				mik = slope(CH[i], CH[k]);
				if (mij == mik) continue;
				if ((mij<mik && (mik-mij)>epsERR) || (mik<mij && (mij-mik)>epsERR)){
					cout << "slope(i, j) = " << mij << endl;
					cout << "slope(i, k) = " << mik << endl;
					cout << "Error in rep(0) = " << rep << ". The convex hull is not a convex polygon !" << endl;
					cout << "Upper lines !isCCW_RL(CH[i],CH[j],CH[k]) to i=" << i << ", j=" << j << ", k=" << k << endl;
					cout << "(CH[i],CH[j],CH[k]) = (" << CH[i] << ", " << CH[j] << "," << CH[k] << ")" << endl;
					cout << "(" << X[CH[i]] << ", " << Y[CH[i]] << ")" << endl;
					cout << "(" << X[CH[j]] << ", " << Y[CH[j]] << ")" << endl;
					cout << "(" << X[CH[k]] << ", " << Y[CH[k]] << ")" << endl;
					exit(1);
				}
			}
		}

		for(i=changeP, j=changeP+1, k=changeP+2; k<nCH; i++, j++, k++){
			if(!isCCW_LR(CH[i],CH[j],CH[k])){
				mij = slope(CH[i], CH[j]);
				mik = slope(CH[i], CH[k]);
				if (mij == mik) continue;
				if ((mij<mik && (mik-mij)>epsERR) || (mik<mij && (mij-mik)>epsERR)){
					cout << "slope(i, j) = " << mij << endl;
					cout << "slope(i, k) = " << mik << endl;
					cout << "Error in rep(0) = " << rep << ". The convex hull is not a convex polygon !" << endl;
					cout << "Lower lines !isCCW_LR(CH[i],CH[j],CH[k]) to i=" << i << ", j=" << j << ", k=" << k << endl;
					cout << "(CH[i],CH[j],CH[k]) = (" << CH[i] << ", " << CH[j] << "," << CH[k] << ")" << endl;
					cout << "(" << X[CH[i]] << ", " << Y[CH[i]] << ")" << endl;
					cout << "(" << X[CH[j]] << ", " << Y[CH[j]] << ")" << endl;
					cout << "(" << X[CH[k]] << ", " << Y[CH[k]] << ")" << endl;
					exit(1);
				}
			}
		}
	}
	//cout << " Test I OK !" << endl;

	mh2 = slope(ri4, le3);
	// II.- test that each point is inside the convex hull...
	//cout << "  II.- Test that each point is inside the convex hull..." << endl;
	for(i=0; i<n; i++){
		xi = X[i];
		if (!isInCHull(i)){
			mi = slope(ri1, i);
			if (mi < mh){
				if (LeftRight){
					// search above: between ri2...le2 --> mh
					sp=changeP; ep=nCH-1;
					x=sp+(ep-sp)/2;
					while(X[CH[x]] < xi || (x+1<nCH && X[CH[x+1]] >= xi)){
						if(X[CH[x]] < xi)
							ep=x-1;
						else
							sp=x+1;
						x=sp+(ep-sp)/2;
					}
				}else{
					// search above: between 0...r12 --> mh
					sp=0; ep=changeP;
					x=ep/2;
					while(X[CH[x]] < xi || (x+1<changeP && X[CH[x+1]] >= xi)){
						if(X[CH[x]] < xi)
							ep=x-1;
						else
							sp=x+1;
						x=sp+(ep-sp)/2;
					}
				}

				if (X[CH[x]] == xi){
					if (Y[CH[x]] < Y[i]){
						cout << "For rep(A) = " << rep << endl;
						cout << "ERROR. Point i=" << i << ": X[i]=" << xi << ", Y[i]=" << Y[i] << " is out of the polygon" << endl;
						cout << "       Between CH[" << x << "]=" << CH[x] << ", and CH[" << x+1 << "]=" << CH[x+1] << endl;
						exit(1);
					}
				}else{
					k=x+1;
					if (k == nCH)
						k=0;

					if(isCCW_RL(CH[x], i, CH[k])){
						mij = slope(CH[x], i);
						mik = slope(CH[x], CH[k]);
						if (mij == mik) continue;
						if ((mij<mik && (mik-mij)>epsERR) || (mik<mij && (mij-mik)>epsERR)){
							cout << "For rep(B) = " << rep << endl;
							cout << "ERROR. Point i=" << i << ": X[i]=" << xi << ", Y[i]=" << Y[i] << " is out of the polygon" << endl;
							cout << "       Between CH[" << x << "]=" << CH[x] << ", and CH[" << k << "]=" << CH[k] << endl;
							exit(1);
						}
					}
				}
			}else{
				mi = slope(i, ri4);
				if (mi > mh2){
					if (LeftRight){
						// search below: between ri2...le2 --> mh
						sp=0; ep=changeP;
						x=sp+(ep-sp)/2;
						while(X[CH[x]] >= xi || (x+1<changeP && X[CH[x+1]] < xi)){
							if(X[CH[x]] >= xi)
								ep=x-1;
							else
								sp=x+1;
							x=sp+(ep-sp)/2;
						}
					}else{
						// search below: between le3...ri4 --> mh2
						sp=changeP; ep=nCH-1;
						x=sp+(ep-sp)/2;
						while(X[CH[x]] >= xi || (x+1<nCH && X[CH[x+1]] < xi)){
							if(X[CH[x]] >= xi)
								ep=x-1;
							else
								sp=x+1;
							x=sp+(ep-sp)/2;
						}
					}
					if (X[CH[x]] == xi){
						if(Y[i] < Y[CH[x]]){
							cout << "For rep(C) = " << rep << endl;
							cout << "ERROR. Point i=" << i << ": X[i]=" << xi << ", Y[i]=" << Y[i] << " is out of the polygon" << endl;
							if (x)
								cout << "       Between CH[" << x-1 << "]=" << CH[x-1] << ", and CH[" << x << "]=" << CH[x] << endl;
							else
								cout << "       Between CH[" << nCH-1 << "]=" << CH[nCH-1] << ", and CH[" << x << "]=" << CH[x] << endl;
							exit(1);
						}
					}else{
						if (x+1 == nCH) x=0;
						if(isCCW_LR(CH[x], i, CH[x+1])){
							mij = slope(CH[x], i);
							mik = slope(CH[x], CH[x+1]);
							if (mij == mik) continue;
							if ((mij<mik && (mik-mij)>epsERR) || (mik<mij && (mij-mik)>epsERR)){
								cout << "For rep(D) = " << rep << endl;
								cout << "ERROR. Point i=" << i << ": X[i]=" << xi << ", Y[i]=" << Y[i] << " is out of the polygon" << endl;
								if (x)
									cout << "       Between CH[" << x-1 << "]=" << CH[x-1] << ", and CH[" << x << "]=" << CH[x] << endl;
								else
									cout << "       Between CH[" << nCH-1 << "]=" << CH[nCH-1] << ", and CH[" << x << "]=" << CH[x] << endl;
								exit(1);
							}
						}
					}
				}
			}
		}
	}
	//cout << " Test II OK !" << endl;
}

void HeapCHull::createFileChull(char filePoints[300]){
	// make extremPointFile...
	char aFile[400];
	strcpy(aFile, "");
	sprintf(aFile, "%s%lu_ext", filePoints, n);

	cout << "extrmePoints File: " << aFile << endl;
	FILE *fExt = fopen(aFile, "w");

	fprintf(fExt, "%lu %f %f\n", ri1, X[ri1], Y[ri1]);
	if (c1 < n){
		fprintf(fExt, "%lu %f %f\n", c1, X[c1], Y[c1]);
		fprintf(fExt, "%lu %f %f\n", up1, X[up1], Y[up1]);
	}else{
		if (ri1 != up1)
			fprintf(fExt, "%lu %f %f\n", up1, X[up1], Y[up1]);
	}

	if (up1 != up2)
		fprintf(fExt, "%lu %f %f\n", up2, X[up2], Y[up2]);
	if (c2 < n){
		fprintf(fExt, "%lu %f %f\n", c2, X[c2], Y[c2]);
		fprintf(fExt, "%lu %f %f\n", le2, X[le2], Y[le2]);
	}else{
		if (up2 != le2)
			fprintf(fExt, "%lu %f %f\n", le2, X[le2], Y[le2]);
	}

	if (le2 != le3)
		fprintf(fExt, "%lu %f %f\n", le3, X[le3], Y[le3]);
	if (c3 < n){
		fprintf(fExt, "%lu %f %f\n", c3, X[c3], Y[c3]);
		fprintf(fExt, "%lu %f %f\n", lo3, X[lo3], Y[lo3]);
	}else{
		if (le3 != lo3)
			fprintf(fExt, "%lu %f %f\n", lo3, X[lo3], Y[lo3]);
	}

	if (lo3 != lo4)
		fprintf(fExt, "%lu %f %f\n", lo4, X[lo4], Y[lo4]);
	if (c4 < n){
		fprintf(fExt, "%lu %f %f\n", c4, X[c4], Y[c4]);
		if (ri1 != ri4)
			fprintf(fExt, "%lu %f %f\n", ri4 ,X[ri4], Y[ri4]);
	}else{
		if (lo4 != ri4 && ri1 != ri4)
			fprintf(fExt, "%lu %f %f\n", ri4 ,X[ri4], Y[ri4]);
	}

	fprintf(fExt, "%lu %f %f\n", ri1, X[ri1], Y[ri1]);
	fclose (fExt);

	// make extremPointFile...
	strcpy(aFile, "");
	sprintf(aFile, "%s%lu_chull", filePoints, n);
	cout << "chull File : " << aFile << endl;
	FILE *fCH = fopen(aFile, "w");
	for(uint k=0; k<nCH; k++){
		fprintf(fCH, "%lu %f %f\n", CH[k] ,X[CH[k]], Y[CH[k]]);
	}
	fprintf(fCH, "%lu %f %f\n", CH[0] ,X[CH[0]], Y[CH[0]]);
	fclose (fCH);
}
//====================================================================================================================

// set the number x as a bitstring sequence in *A. In the range of bits [ini, .. ini+len-1] of *A. Here x has len bits
void HeapCHull::setNum64(ulong *A, ulong ini, uint len, ulong x) {
	ulong i=ini>>BW64, j=ini-(i<<BW64);

	if ((j+len)>W64){
		ulong myMask = ~(~0ul >> j);
		A[i] = (A[i] & myMask) | (x >> (j+len-W64));
		myMask = ~0ul >> (j+len-W64);
		A[i+1] = (A[i+1] & myMask) | (x << (WW64-j-len));
	}else{
		ulong myMask = (~0ul >> j) ^ (~0ul << (W64-j-len)); // XOR: 1^1=0^0=0; 0^1=1^0=1
		A[i] = (A[i] & myMask) | (x << (W64-j-len));
	}
}

// return (in a unsigned long integer) the number in A from bits of position 'ini' to 'ini+len-1'
ulong HeapCHull::getNum64(ulong *A, ulong ini, uint len){
	ulong i=ini>>BW64, j=ini-(i<<BW64);
	ulong result = (A[i] << j) >> (W64-len);

	if (j+len > W64)
		result = result | (A[i+1] >> (WW64-j-len));

	return result;
}

HeapCHull::~HeapCHull() {
	delete [] Q1;
	delete [] Q2;
	delete [] Q3;
	delete [] Q4;

	if (isQ40){
		delete [] _Q1;
		delete [] _Q2;
		delete [] _Q3;
		delete [] _Q4;
	}
	delete [] CH;
	//cout << "Destroyed !" << endl;
}

// ##########################################################################3
// methods not used...
// set the new item at top of the queue Q[1..k]
void HeapCHull::setTopMinQ(ulong *Q, ulong k, float *V){
	ulong pm, qm, qmm, qk, m=2, i=1;
	qk = getNum64(Q,k*lgn,lgn);
	float nqk = V[qk]; 					// nqk = V[Q[k]];

	while(m<k){
		pm = m*lgn;
		qm = getNum64(Q,pm,lgn);
		if (m+1<k){
			qmm = getNum64(Q,pm+lgn,lgn);
			if (V[qm] > V[qmm]){
				m++;
				qm=qmm;
			}
		}

		if(nqk > V[qm]){
			setNum64(Q,i*lgn,lgn,qm);	// Q[i] = qm;
			i = m;
			m <<= 1;
		}else
			break;
	}
	setNum64(Q,i*lgn,lgn,qk);			//Q[i] = Q[k];

}
// set the new item at top of the queue Q[1..k]
void HeapCHull::setTopMaxQ(ulong *Q, ulong k, float *V){
	ulong pm, qm, qmm, qk, m=2, i=1;
	qk = getNum64(Q,k*lgn,lgn);
	float nqk = V[qk]; 					// nqk = V[Q[k]];

	while(m<k){
		pm = m*lgn;
		qm = getNum64(Q,pm,lgn);
		if (m+1<k){
			qmm = getNum64(Q,pm+lgn,lgn);
			if (V[qm] < V[qmm]){
				m++;
				qm=qmm;
			}
		}

		if(nqk < V[qm]){
			setNum64(Q,i*lgn,lgn,qm);	// Q[i] = qm;
			i = m;
			m <<= 1;
		}else
			break;
	}
	setNum64(Q,i*lgn,lgn,qk);			//Q[i] = Q[k];

}
// to create the minimum priority queue Q[1..m]
void HeapCHull::createMinHeap(ulong *Q, ulong m, float *V){
	ulong i,j,k,t,qt;
	float v;

	for(i=2; i<m; i++){
		j=i;
		k=getNum64(Q,j*lgn,lgn);			//k=Q[j];
		v=V[k];
		t=j>>1;
		qt=getNum64(Q,t*lgn,lgn);
		while(t && V[qt]>v){
			setNum64(Q,j*lgn,lgn, qt); 	//Q[j] = Q[t];
			j=t;
			t>>=1;
			qt=getNum64(Q,t*lgn,lgn);
		}
		setNum64(Q,j*lgn,lgn,k); 	//Q[j]=k;
	}
}

// to create the maximum priority queue Q[1..m]
void HeapCHull::createMaxHeap(ulong *Q, ulong m, float *V){
	ulong i,j,k,t,qt;
	float v;

	for(i=2; i<m; i++){
		j=i;
		k=getNum64(Q,j*lgn,lgn);				//k=Q[j];
		v=V[k];
		t=j>>1;
		qt=getNum64(Q,t*lgn,lgn);
		while(t && V[qt]<v){
			setNum64(Q,j*lgn,lgn,qt);	//Q[j] = qt;
			j=t;
			t>>=1;
			qt=getNum64(Q,t*lgn,lgn);
		}
		setNum64(Q,j*lgn,lgn,k); 	//Q[j]=k;
	}
}

// Cambiar y dejar los if con < primeros a los que tienen ==


void HeapCHull::seacrhExtremePoints3(){
	float xi,yi;

	c1 = c2 = c3 = c4 = n+1;
	xc1 = xc4 = -1*FLT_MAX;
	xc2 = xc3 = FLT_MAX;
	yc1 = yc2 = -1*FLT_MAX;
	yc3 = yc4 = FLT_MAX;
	xri1 = xup1 = xup2 = xle2 = xle3 = xlo3 = xlo4 = xri4 = X[0];
	yri1 = yup1 = yup2 = yle2 = yle3 = ylo3 = ylo4 = yri4 = Y[0];
	ri1 = up1 = up2 = le2 = le3 = lo3 = lo4 = ri4 = 0;
	for(ulong i=1; i<n; i++){
		xi=X[i];yi=Y[i];
		if (xi < xle2){
			if (xc2-yc2 > xle2-yle2){
				xc2=xle2;
				yc2=yle2;
				c2 = le2;
			}
			if (xc3+yc3 > xle3+yle3){
				xc3=xle3;
				yc3=yle3;
				c3 = le3;
			}

			xle2 = xle3 = xi;
			yle2 = yle3 = yi;
			le2 = le3 = i;
		}else{

		}
	}

}

void HeapCHull::seacrhExtremePoints2(){
	float xi,yi;

	c1 = c2 = c3 = c4 = n+1;
	xc1 = xc4 = -1*FLT_MAX;
	xc2 = xc3 = FLT_MAX;
	yc1 = yc2 = -1*FLT_MAX;
	yc3 = yc4 = FLT_MAX;
	xri1 = xup1 = xup2 = xle2 = xle3 = xlo3 = xlo4 = xri4 = X[0];
	yri1 = yup1 = yup2 = yle2 = yle3 = ylo3 = ylo4 = yri4 = Y[0];
	ri1 = up1 = up2 = le2 = le3 = lo3 = lo4 = ri4 = 0;
	for(ulong i=1; i<n; i++){
		xi=X[i];yi=Y[i];
		if (xi < xle2){
			if (xc2-yc2 > xle2-yle2){
				xc2=xle2;
				yc2=yle2;
				c2 = le2;
			}
			if (xc3+yc3 > xle3+yle3){
				xc3=xle3;
				yc3=yle3;
				c3 = le3;
			}

			xle2 = xle3 = xi;
			yle2 = yle3 = yi;
			le2 = le3 = i;
		}else{
			if(xi == xle2){ // no actualizar c2 o c3 pq tiene el mismo x q xi
				if(yi > yle2){
					xle2 = xi;
					yle2 = yi;
					le2 = i;
				}else{
					if(yi < yle3){
						xle3 = xi;
						yle3 = yi;
						le3 = i;
					}
				}
			}else{
				if(xi > xri1){
					if (xc1+yc1 < xri1+yri1){
						xc1=xri1;
						yc1=yri1;
						c1 = ri1;
					}
					if (yc4-xc4 > yri4-xri4){
						xc4=xri4;
						yc4=yri4;
						c4 = ri4;
					}
					xri1 = xri4 = xi;
					yri1 = yri4 = yi;
					ri1 = ri4 = i;
				}else{
					if(xi == xri1){
						if(yi > yri1){
							xri1 = xi;
							yri1 = yi;
							ri1 = i;
						}else{
							if(yi < yri4){
								xri4 = xi;
								yri4 = yi;
								ri4 = i;
							}
						}
					}else{
						if(yi < ylo3){
							if (xc3+yc3 > xlo3+ylo3){
								xc3=xlo3;
								yc3=ylo3;
								c3 = lo3;
							}
							if (yc4-xc4 > ylo4-xlo4){
								xc4=xlo4;
								yc4=ylo4;
								c4 = lo4;
							}

							xlo3 = xlo4 = xi;
							ylo3 = ylo4 = yi;
							lo3 = lo4 = i;
						}else{
							if(yi == ylo3){
								if(xi < xlo3){
									xlo3 = xi;
									ylo3 = yi;
									lo3 = i;
								}else{
									if(xi > xlo4){
										xlo4 = xi;
										ylo4 = yi;
										lo4 = i;
									}
								}
							}else{
								if(yi > yup2){
									if (xc1+yc1 < xup1+yup1){
										xc1=xup1;
										yc1=yup1;
										c1 = up1;
									}
									if (xc2-yc2 > xup2-yup2){
										xc2=xup2;
										yc2=yup2;
										c2 = up2;
									}

									xup2 = xup1 = xi;
									yup2 = yup1 = yi;
									up1 = up2 = i;
								}else{
									if(yi == yup2){
										if(xi < xup2){
											xup2 = xi;
											yup2 = yi;
											up2 = i;
										}else{
											if(xi > xup1){
												xup1 = xi;
												yup1 = yi;
												up1 = i;
											}
										}
									}else{
										if (xc1+yc1 < xi+yi){
											xc1=xi;
											yc1=yi;
											c1 = i;
										}else{
											if (xc3+yc3 > xi+yi){
												xc3=xi;
												yc3=yi;
												c3 = i;
											}else{
												if (xi < xc2){
													if (yi>=yc2 || xc2-xi > yc2-yi){
														xc2=xi;
														yc2=yi;
														c2 = i;
													}
												}else{
													if (yi > yc2){
														if (xi<=xc2 || yi-yc2 > xi-xc2){
															xc2=xi;
															yc2=yi;
															c2 = i;
														}
													}
												}
												if (xi > xc4){
													if (yi<=yc4 || xi-xc4 > yi-yc4){
														xc4=xi;
														yc4=yi;
														c4 = i;
													}
												}else{
													if (yi < yc4){
														if (xi>=xc4 || yc4-yi > xc4-xi){
															xc4=xi;
															yc4=yi;
															c4 = i;
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
/*
void HeapCHull::seacrhExtremePoints_old(){
	float xi,yi;
	bool seei=true;

	ri1 = up1 = up2 = le2 = le3 = lo3 = lo4 = ri4 = c1 = c2 = c3 = c4 = 0;
	for(ulong i=1; i<n; i++, seei=true){
		xi=X[i];yi=Y[i];
		if(xi < X[le2]){
			if (X[c2]-Y[c2] > X[le2]-Y[le2])
				c2=le2;
			if (X[c3]+Y[c3] > X[le2]+Y[le2])
				c3=le2;
			le2 = le3 = i;
			seei=false;
		}else{
			if(xi == X[le2]){
				if(yi > Y[le2]){
					le2 = i;
					seei=false;
				}
				if(yi < Y[le3]){
					le3 = i;
					seei=false;
				}
			}
			if(xi > X[ri1]){
				if (X[c1]+Y[c1] < X[ri1]+Y[ri1])
					c1=ri1;
				if (Y[c4]-X[c4] > Y[ri1]-X[ri1])
					c4=ri1;
				ri1 = ri4 = i;
				seei=false;
			}else{
				if(xi == X[ri1]){
					if(yi > Y[ri1]){
						ri1 = i;
						seei=false;
					}
					if(yi < Y[ri4]){
						ri4 = i;
						seei=false;
					}
				}
			}
		}

		if(yi < Y[lo3]){
			if (X[c3]+Y[c3] > X[lo3]+Y[lo3])
				c3=lo3;
			if (Y[c4]-X[c4] > Y[lo3]-X[lo3])
				c4=lo3;

			lo3 = lo4 = i;
			seei=false;
		}else{
			if(yi == Y[lo3]){
				if(xi < X[lo3]){
					lo3 = i;
					seei=false;
				}
				if(xi > X[lo4]){
					lo4 = i;
					seei=false;
				}
			}
			if(yi > Y[up2]){
				if (X[c1]+Y[c1] < X[up2]+Y[up2])
					c1=up2;
				if (X[c2]-Y[c2] > X[up2]-Y[up2])
					c2=up2;

				up2 = up1 = i;
				seei=false;
			}else{
				if(yi == Y[up2]){
					if(xi < X[up2]){
						up2 = i;
						seei=false;
					}
					if(xi > X[up1]){
						up1 = i;
						seei=false;
					}
				}
			}
		}

		if (seei){
			if (X[c1]+Y[c1] < xi+yi)
				c1=i;
			if (X[c2]-Y[c2] > xi-yi)
				c2=i;
			if (X[c3]+Y[c3] > xi+yi)
				c3=i;
			if (Y[c4]-X[c4] > yi-xi)
				c4=i;
		}
	}
}*/

// ##########################################################################3
} /* namespace hch */
