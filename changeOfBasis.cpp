//Copyright 2014 Scott M. Johnson
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include "ChangeOfBasis.h"

namespace cob
{

inline bool matchDirection( int from, int to)
{
	return ((from & 3) == (to & 3));
}

// The caseNumber is a representation of the column matrix MAtoB that transforms column vectors 
// from some frame A to another frame B.
//
// [ Vb ] = [ MAtoB ] . [ Va ]
//
// MAtoB is factored into a permutation Matrix P and a diagonal sign matrix S.
//
// [ MAtoB ] = [ P ] . [ S ]
//                        [ 1  0  0 ] [ -1  0  0 ]
// example : [ MAtoB ] =  [ 0  0  1 ] [  0  1  0 ]
//                        [ 0  1  0 ] [  0  0 -1 ]
//
// Each of the six permutation matrices P is arbitrarily assigned an index from 0 to 5.
// Each of the eight sign matrices is assigned an index from 0 to 7.
//  caseNumber = permutationIndex * 8 + signMatrixIndex
// Once you have a case number you can save it and keep using it later.
// The goal of these functions is to perform the change of basis without adding
// any floating point round off errors.
//
// example:  Prepare to convert rotations or vectors from the BVH Frame to the Unreal 3 Frame
// int caseNumber = getCaseNumber( bvhFrame, Unreal3Frame );

int getCaseNumber( const triple &from, const triple &to)
{
	int p(0), s(0);

	if (matchDirection(from.a, to.a))
	{
		s |= (from.a != to.a) ? 0x04 : 0;
		if (matchDirection(from.b, to.b))
		{
			// a, b, c
			p = 0;
			s |= (from.b != to.b) ? 0x02 : 0;
			s |= (from.c != to.c) ? 0x01 : 0;
		}
		else // matchDirection(from.b, to.c) matches
		{
			// a, c, b
			p = 1;
			s |= (from.b != to.c) ? 0x02 : 0;
			s |= (from.c != to.b) ? 0x01 : 0;
		}
	}
	else if (matchDirection(from.a, to.b))
	{
		s |= (from.a != to.b) ? 0x04 : 0;
		if (matchDirection(from.b, to.a))
		{
			// b, a, c
			p = 2;
			s |= (from.b != to.a) ? 0x02 : 0;
			s |= (from.c != to.c) ? 0x01 : 0;
		}
		else if (matchDirection(from.b, to.c))
		{
			// c, a, b 
			p = 3;
			s |= (from.b != to.c) ? 0x02 : 0;
			s |= (from.c != to.a) ? 0x01 : 0;
		}
	}
	else	// matchDirection(from.a, to.c)
	{
		s |= (from.a != to.c) ? 0x04 : 0;
		if (matchDirection(from.b, to.a))
		{
			// b, c, a
			p = 4;
			s |= (from.b != to.a) ? 0x02 : 0;
			s |= (from.c != to.b) ? 0x01 : 0;
		}
		else 
		{
			// c, b, a 
			p = 5;
			s |= (from.b != to.b) ? 0x02 : 0;
			s |= (from.c != to.a) ? 0x01 : 0;
		}
	}

	return p * 8 + s;
}

void matrixCob3x3(int caseNumber, 
	double &a00, double &a01, double &a02,
	double &a10, double &a11, double &a12,
	double &a20, double &a21, double &a22)
{
	double m00, m01, m02;
	double m10, m11, m12;
	double m20, m21, m22;

	m00 = a00;  m01 = a01; m02 = a02;
	m10 = a10;  m11 = a11; m12 = a12;
	m20 = a20;  m21 = a21; m22 = a22;

	switch(caseNumber)
	{
		case 0 :
			a00 =  m00; a01 =  m01; a02 =  m02;
			a10 =  m10; a11 =  m11; a12 =  m12;
			a20 =  m20; a21 =  m21; a22 =  m22;
			break;
		case 1 :
			a00 =  m00; a01 =  m01; a02 = -m02;
			a10 =  m10; a11 =  m11; a12 = -m12;
			a20 = -m20; a21 = -m21; a22 =  m22;
			break;
		case 2 :
			a00 =  m00; a01 = -m01; a02 =  m02;
			a10 = -m10; a11 =  m11; a12 = -m12;
			a20 =  m20; a21 = -m21; a22 =  m22;
			break;
		case 3 :
			a00 =  m00; a01 = -m01; a02 = -m02;
			a10 = -m10; a11 =  m11; a12 =  m12;
			a20 = -m20; a21 =  m21; a22 =  m22;
			break;
		case 4 :
			a00 =  m00; a01 = -m01; a02 = -m02;
			a10 = -m10; a11 =  m11; a12 =  m12;
			a20 = -m20; a21 =  m21; a22 =  m22;
			break;
		case 5 :
			a00 =  m00; a01 = -m01; a02 =  m02;
			a10 = -m10; a11 =  m11; a12 = -m12;
			a20 =  m20; a21 = -m21; a22 =  m22;
			break;
		case 6 :
			a00 =  m00; a01 =  m01; a02 = -m02;
			a10 =  m10; a11 =  m11; a12 = -m12;
			a20 = -m20; a21 = -m21; a22 =  m22;
			break;
		case 7 :
			a00 =  m00; a01 =  m01; a02 =  m02;
			a10 =  m10; a11 =  m11; a12 =  m12;
			a20 =  m20; a21 =  m21; a22 =  m22;
			break;
		case 8 :
			a00 =  m00; a01 =  m02; a02 =  m01;
			a10 =  m20; a11 =  m22; a12 =  m21;
			a20 =  m10; a21 =  m12; a22 =  m11;
			break;
		case 9 :
			a00 =  m00; a01 = -m02; a02 =  m01;
			a10 = -m20; a11 =  m22; a12 = -m21;
			a20 =  m10; a21 = -m12; a22 =  m11;
			break;
		case 10 :
			a00 =  m00; a01 =  m02; a02 = -m01;
			a10 =  m20; a11 =  m22; a12 = -m21;
			a20 = -m10; a21 = -m12; a22 =  m11;
			break;
		case 11 :
			a00 =  m00; a01 = -m02; a02 = -m01;
			a10 = -m20; a11 =  m22; a12 =  m21;
			a20 = -m10; a21 =  m12; a22 =  m11;
			break;
		case 12 :
			a00 =  m00; a01 = -m02; a02 = -m01;
			a10 = -m20; a11 =  m22; a12 =  m21;
			a20 = -m10; a21 =  m12; a22 =  m11;
			break;
		case 13 :
			a00 =  m00; a01 =  m02; a02 = -m01;
			a10 =  m20; a11 =  m22; a12 = -m21;
			a20 = -m10; a21 = -m12; a22 =  m11;
			break;
		case 14 :
			a00 =  m00; a01 = -m02; a02 =  m01;
			a10 = -m20; a11 =  m22; a12 = -m21;
			a20 =  m10; a21 = -m12; a22 =  m11;
			break;
		case 15 :
			a00 =  m00; a01 =  m02; a02 =  m01;
			a10 =  m20; a11 =  m22; a12 =  m21;
			a20 =  m10; a21 =  m12; a22 =  m11;
			break;
		case 16 :
			a00 =  m11; a01 =  m10; a02 =  m12;
			a10 =  m01; a11 =  m00; a12 =  m02;
			a20 =  m21; a21 =  m20; a22 =  m22;
			break;
		case 17 :
			a00 =  m11; a01 =  m10; a02 = -m12;
			a10 =  m01; a11 =  m00; a12 = -m02;
			a20 = -m21; a21 = -m20; a22 =  m22;
			break;
		case 18 :
			a00 =  m11; a01 = -m10; a02 = -m12;
			a10 = -m01; a11 =  m00; a12 =  m02;
			a20 = -m21; a21 =  m20; a22 =  m22;
			break;
		case 19 :
			a00 =  m11; a01 = -m10; a02 =  m12;
			a10 = -m01; a11 =  m00; a12 = -m02;
			a20 =  m21; a21 = -m20; a22 =  m22;
			break;
		case 20 :
			a00 =  m11; a01 = -m10; a02 =  m12;
			a10 = -m01; a11 =  m00; a12 = -m02;
			a20 =  m21; a21 = -m20; a22 =  m22;
			break;
		case 21 :
			a00 =  m11; a01 = -m10; a02 = -m12;
			a10 = -m01; a11 =  m00; a12 =  m02;
			a20 = -m21; a21 =  m20; a22 =  m22;
			break;
		case 22 :
			a00 =  m11; a01 =  m10; a02 = -m12;
			a10 =  m01; a11 =  m00; a12 = -m02;
			a20 = -m21; a21 = -m20; a22 =  m22;
			break;
		case 23 :
			a00 =  m11; a01 =  m10; a02 =  m12;
			a10 =  m01; a11 =  m00; a12 =  m02;
			a20 =  m21; a21 =  m20; a22 =  m22;
			break;
		case 24 :
			a00 =  m22; a01 =  m20; a02 =  m21;
			a10 =  m02; a11 =  m00; a12 =  m01;
			a20 =  m12; a21 =  m10; a22 =  m11;
			break;
		case 25 :
			a00 =  m22; a01 = -m20; a02 = -m21;
			a10 = -m02; a11 =  m00; a12 =  m01;
			a20 = -m12; a21 =  m10; a22 =  m11;
			break;
		case 26 :
			a00 =  m22; a01 =  m20; a02 = -m21;
			a10 =  m02; a11 =  m00; a12 = -m01;
			a20 = -m12; a21 = -m10; a22 =  m11;
			break;
		case 27 :
			a00 =  m22; a01 = -m20; a02 =  m21;
			a10 = -m02; a11 =  m00; a12 = -m01;
			a20 =  m12; a21 = -m10; a22 =  m11;
			break;
		case 28 :
			a00 =  m22; a01 = -m20; a02 =  m21;
			a10 = -m02; a11 =  m00; a12 = -m01;
			a20 =  m12; a21 = -m10; a22 =  m11;
			break;
		case 29 :
			a00 =  m22; a01 =  m20; a02 = -m21;
			a10 =  m02; a11 =  m00; a12 = -m01;
			a20 = -m12; a21 = -m10; a22 =  m11;
			break;
		case 30 :
			a00 =  m22; a01 = -m20; a02 = -m21;
			a10 = -m02; a11 =  m00; a12 =  m01;
			a20 = -m12; a21 =  m10; a22 =  m11;
			break;
		case 31 :
			a00 =  m22; a01 =  m20; a02 =  m21;
			a10 =  m02; a11 =  m00; a12 =  m01;
			a20 =  m12; a21 =  m10; a22 =  m11;
			break;
		case 32 :
			a00 =  m11; a01 =  m12; a02 =  m10;
			a10 =  m21; a11 =  m22; a12 =  m20;
			a20 =  m01; a21 =  m02; a22 =  m00;
			break;
		case 33 :
			a00 =  m11; a01 = -m12; a02 =  m10;
			a10 = -m21; a11 =  m22; a12 = -m20;
			a20 =  m01; a21 = -m02; a22 =  m00;
			break;
		case 34 :
			a00 =  m11; a01 = -m12; a02 = -m10;
			a10 = -m21; a11 =  m22; a12 =  m20;
			a20 = -m01; a21 =  m02; a22 =  m00;
			break;
		case 35 :
			a00 =  m11; a01 =  m12; a02 = -m10;
			a10 =  m21; a11 =  m22; a12 = -m20;
			a20 = -m01; a21 = -m02; a22 =  m00;
			break;
		case 36 :
			a00 =  m11; a01 =  m12; a02 = -m10;
			a10 =  m21; a11 =  m22; a12 = -m20;
			a20 = -m01; a21 = -m02; a22 =  m00;
			break;
		case 37 :
			a00 =  m11; a01 = -m12; a02 = -m10;
			a10 = -m21; a11 =  m22; a12 =  m20;
			a20 = -m01; a21 =  m02; a22 =  m00;
			break;
		case 38 :
			a00 =  m11; a01 = -m12; a02 =  m10;
			a10 = -m21; a11 =  m22; a12 = -m20;
			a20 =  m01; a21 = -m02; a22 =  m00;
			break;
		case 39 :
			a00 =  m11; a01 =  m12; a02 =  m10;
			a10 =  m21; a11 =  m22; a12 =  m20;
			a20 =  m01; a21 =  m02; a22 =  m00;
			break;
		case 40 :
			a00 =  m22; a01 =  m21; a02 =  m20;
			a10 =  m12; a11 =  m11; a12 =  m10;
			a20 =  m02; a21 =  m01; a22 =  m00;
			break;
		case 41 :
			a00 =  m22; a01 = -m21; a02 = -m20;
			a10 = -m12; a11 =  m11; a12 =  m10;
			a20 = -m02; a21 =  m01; a22 =  m00;
			break;
		case 42 :
			a00 =  m22; a01 = -m21; a02 =  m20;
			a10 = -m12; a11 =  m11; a12 = -m10;
			a20 =  m02; a21 = -m01; a22 =  m00;
			break;
		case 43 :
			a00 =  m22; a01 =  m21; a02 = -m20;
			a10 =  m12; a11 =  m11; a12 = -m10;
			a20 = -m02; a21 = -m01; a22 =  m00;
			break;
		case 44 :
			a00 =  m22; a01 =  m21; a02 = -m20;
			a10 =  m12; a11 =  m11; a12 = -m10;
			a20 = -m02; a21 = -m01; a22 =  m00;
			break;
		case 45 :
			a00 =  m22; a01 = -m21; a02 =  m20;
			a10 = -m12; a11 =  m11; a12 = -m10;
			a20 =  m02; a21 = -m01; a22 =  m00;
			break;
		case 46 :
			a00 =  m22; a01 = -m21; a02 = -m20;
			a10 = -m12; a11 =  m11; a12 =  m10;
			a20 = -m02; a21 =  m01; a22 =  m00;
			break;
		case 47 :
			a00 =  m22; a01 =  m21; a02 =  m20;
			a10 =  m12; a11 =  m11; a12 =  m10;
			a20 =  m02; a21 =  m01; a22 =  m00;
			break;
	}
}

void getAtoBMatrix(int caseNumber, 
	double &m00, double &m01, double &m02,
	double &m10, double &m11, double &m12,
	double &m20, double &m21, double &m22)
{
	switch(caseNumber)
	{
		case 0 :
			m00 = 1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = 1.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = 1.0;
			break;
		case 1 :
			m00 = 1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = 1.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = -1.0;
			break;
		case 2 :
			m00 = 1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = -1.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = 1.0;
			break;
		case 3 :
			m00 = 1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = -1.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = -1.0;
			break;
		case 4 :
			m00 = -1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = 1.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = 1.0;
			break;
		case 5 :
			m00 = -1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = 1.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = -1.0;
			break;
		case 6 :
			m00 = -1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = -1.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = 1.0;
			break;
		case 7 :
			m00 = -1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = -1.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = -1.0;
			break;
		case 8 :
			m00 = 1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = 1.0;
			m20 = 0.0; m21 = 1.0; m22 = 0.0;
			break;
		case 9 :
			m00 = 1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = -1.0;
			m20 = 0.0; m21 = 1.0; m22 = 0.0;
			break;
		case 10 :
			m00 = 1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = 1.0;
			m20 = 0.0; m21 = -1.0; m22 = 0.0;
			break;
		case 11 :
			m00 = 1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = -1.0;
			m20 = 0.0; m21 = -1.0; m22 = 0.0;
			break;
		case 12 :
			m00 = -1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = 1.0;
			m20 = 0.0; m21 = 1.0; m22 = 0.0;
			break;
		case 13 :
			m00 = -1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = -1.0;
			m20 = 0.0; m21 = 1.0; m22 = 0.0;
			break;
		case 14 :
			m00 = -1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = 1.0;
			m20 = 0.0; m21 = -1.0; m22 = 0.0;
			break;
		case 15 :
			m00 = -1.0; m01 = 0.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = -1.0;
			m20 = 0.0; m21 = -1.0; m22 = 0.0;
			break;
		case 16 :
			m00 = 0.0; m01 = 1.0; m02 = 0.0;
			m10 = 1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = 1.0;
			break;
		case 17 :
			m00 = 0.0; m01 = 1.0; m02 = 0.0;
			m10 = 1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = -1.0;
			break;
		case 18 :
			m00 = 0.0; m01 = -1.0; m02 = 0.0;
			m10 = 1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = 1.0;
			break;
		case 19 :
			m00 = 0.0; m01 = -1.0; m02 = 0.0;
			m10 = 1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = -1.0;
			break;
		case 20 :
			m00 = 0.0; m01 = 1.0; m02 = 0.0;
			m10 = -1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = 1.0;
			break;
		case 21 :
			m00 = 0.0; m01 = 1.0; m02 = 0.0;
			m10 = -1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = -1.0;
			break;
		case 22 :
			m00 = 0.0; m01 = -1.0; m02 = 0.0;
			m10 = -1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = 1.0;
			break;
		case 23 :
			m00 = 0.0; m01 = -1.0; m02 = 0.0;
			m10 = -1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = 0.0; m22 = -1.0;
			break;
		case 24 :
			m00 = 0.0; m01 = 0.0; m02 = 1.0;
			m10 = 1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = 1.0; m22 = 0.0;
			break;
		case 25 :
			m00 = 0.0; m01 = 0.0; m02 = -1.0;
			m10 = 1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = 1.0; m22 = 0.0;
			break;
		case 26 :
			m00 = 0.0; m01 = 0.0; m02 = 1.0;
			m10 = 1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = -1.0; m22 = 0.0;
			break;
		case 27 :
			m00 = 0.0; m01 = 0.0; m02 = -1.0;
			m10 = 1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = -1.0; m22 = 0.0;
			break;
		case 28 :
			m00 = 0.0; m01 = 0.0; m02 = 1.0;
			m10 = -1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = 1.0; m22 = 0.0;
			break;
		case 29 :
			m00 = 0.0; m01 = 0.0; m02 = -1.0;
			m10 = -1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = 1.0; m22 = 0.0;
			break;
		case 30 :
			m00 = 0.0; m01 = 0.0; m02 = 1.0;
			m10 = -1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = -1.0; m22 = 0.0;
			break;
		case 31 :
			m00 = 0.0; m01 = 0.0; m02 = -1.0;
			m10 = -1.0; m11 = 0.0; m12 = 0.0;
			m20 = 0.0; m21 = -1.0; m22 = 0.0;
			break;
		case 32 :
			m00 = 0.0; m01 = 1.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = 1.0;
			m20 = 1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 33 :
			m00 = 0.0; m01 = 1.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = -1.0;
			m20 = 1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 34 :
			m00 = 0.0; m01 = -1.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = 1.0;
			m20 = 1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 35 :
			m00 = 0.0; m01 = -1.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = -1.0;
			m20 = 1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 36 :
			m00 = 0.0; m01 = 1.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = 1.0;
			m20 = -1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 37 :
			m00 = 0.0; m01 = 1.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = -1.0;
			m20 = -1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 38 :
			m00 = 0.0; m01 = -1.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = 1.0;
			m20 = -1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 39 :
			m00 = 0.0; m01 = -1.0; m02 = 0.0;
			m10 = 0.0; m11 = 0.0; m12 = -1.0;
			m20 = -1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 40 :
			m00 = 0.0; m01 = 0.0; m02 = 1.0;
			m10 = 0.0; m11 = 1.0; m12 = 0.0;
			m20 = 1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 41 :
			m00 = 0.0; m01 = 0.0; m02 = -1.0;
			m10 = 0.0; m11 = 1.0; m12 = 0.0;
			m20 = 1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 42 :
			m00 = 0.0; m01 = 0.0; m02 = 1.0;
			m10 = 0.0; m11 = -1.0; m12 = 0.0;
			m20 = 1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 43 :
			m00 = 0.0; m01 = 0.0; m02 = -1.0;
			m10 = 0.0; m11 = -1.0; m12 = 0.0;
			m20 = 1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 44 :
			m00 = 0.0; m01 = 0.0; m02 = 1.0;
			m10 = 0.0; m11 = 1.0; m12 = 0.0;
			m20 = -1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 45 :
			m00 = 0.0; m01 = 0.0; m02 = -1.0;
			m10 = 0.0; m11 = 1.0; m12 = 0.0;
			m20 = -1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 46 :
			m00 = 0.0; m01 = 0.0; m02 = 1.0;
			m10 = 0.0; m11 = -1.0; m12 = 0.0;
			m20 = -1.0; m21 = 0.0; m22 = 0.0;
			break;
		case 47 :
			m00 = 0.0; m01 = 0.0; m02 = -1.0;
			m10 = 0.0; m11 = -1.0; m12 = 0.0;
			m20 = -1.0; m21 = 0.0; m22 = 0.0;
			break;
	}
}

void quatCob( int caseNumber, double &qx, double &qy, double &qz, double & /* qw */ )
{
	double tx, ty, tz;
	tx = qx; ty = qy; tz = qz;

	switch(caseNumber)
	{
		case 0:   qx =  tx; qy =  ty; qz =  tz;  break;
		case 1:   qx = -tx; qy = -ty; qz =  tz;  break;
		case 2:   qx = -tx; qy =  ty; qz = -tz;  break;
		case 3:   qx =  tx; qy = -ty; qz = -tz;  break;
		case 4:   qx =  tx; qy = -ty; qz = -tz;  break;
		case 5:   qx = -tx; qy =  ty; qz = -tz;  break;
		case 6:   qx = -tx; qy = -ty; qz =  tz;  break;
		case 7:   qx =  tx; qy =  ty; qz =  tz;  break;
		case 8:   qx = -tx; qy = -tz; qz = -ty;  break;
		case 9:   qx =  tx; qy = -tz; qz =  ty;  break;
		case 10:   qx =  tx; qy =  tz; qz = -ty;  break;
		case 11:   qx = -tx; qy =  tz; qz =  ty;  break;
		case 12:   qx = -tx; qy =  tz; qz =  ty;  break;
		case 13:   qx =  tx; qy =  tz; qz = -ty;  break;
		case 14:   qx =  tx; qy = -tz; qz =  ty;  break;
		case 15:   qx = -tx; qy = -tz; qz = -ty;  break;
		case 16:   qx = -ty; qy = -tx; qz = -tz;  break;
		case 17:   qx =  ty; qy =  tx; qz = -tz;  break;
		case 18:   qx = -ty; qy =  tx; qz =  tz;  break;
		case 19:   qx =  ty; qy = -tx; qz =  tz;  break;
		case 20:   qx =  ty; qy = -tx; qz =  tz;  break;
		case 21:   qx = -ty; qy =  tx; qz =  tz;  break;
		case 22:   qx =  ty; qy =  tx; qz = -tz;  break;
		case 23:   qx = -ty; qy = -tx; qz = -tz;  break;
		case 24:   qx =  tz; qy =  tx; qz =  ty;  break;
		case 25:   qx =  tz; qy = -tx; qz = -ty;  break;
		case 26:   qx = -tz; qy = -tx; qz =  ty;  break;
		case 27:   qx = -tz; qy =  tx; qz = -ty;  break;
		case 28:   qx = -tz; qy =  tx; qz = -ty;  break;
		case 29:   qx = -tz; qy = -tx; qz =  ty;  break;
		case 30:   qx =  tz; qy = -tx; qz = -ty;  break;
		case 31:   qx =  tz; qy =  tx; qz =  ty;  break;
		case 32:   qx =  ty; qy =  tz; qz =  tx;  break;
		case 33:   qx = -ty; qy =  tz; qz = -tx;  break;
		case 34:   qx =  ty; qy = -tz; qz = -tx;  break;
		case 35:   qx = -ty; qy = -tz; qz =  tx;  break;
		case 36:   qx = -ty; qy = -tz; qz =  tx;  break;
		case 37:   qx =  ty; qy = -tz; qz = -tx;  break;
		case 38:   qx = -ty; qy =  tz; qz = -tx;  break;
		case 39:   qx =  ty; qy =  tz; qz =  tx;  break;
		case 40:   qx = -tz; qy = -ty; qz = -tx;  break;
		case 41:   qx = -tz; qy =  ty; qz =  tx;  break;
		case 42:   qx =  tz; qy = -ty; qz =  tx;  break;
		case 43:   qx =  tz; qy =  ty; qz = -tx;  break;
		case 44:   qx =  tz; qy =  ty; qz = -tx;  break;
		case 45:   qx =  tz; qy = -ty; qz =  tx;  break;
		case 46:   qx = -tz; qy =  ty; qz =  tx;  break;
		case 47:   qx = -tz; qy = -ty; qz = -tx;  break;
	}
}

void vectorCob( int caseNumber, double &vx, double &vy, double &vz )
{
	double tx, ty, tz;
	tx = vx; ty = vy; tz = vz;

	switch(caseNumber)
	{
		case 0:   vx =  tx; vy =  ty; vz =  tz;  break;
		case 1:   vx =  tx; vy =  ty; vz = -tz;  break;
		case 2:   vx =  tx; vy = -ty; vz =  tz;  break;
		case 3:   vx =  tx; vy = -ty; vz = -tz;  break;
		case 4:   vx = -tx; vy =  ty; vz =  tz;  break;
		case 5:   vx = -tx; vy =  ty; vz = -tz;  break;
		case 6:   vx = -tx; vy = -ty; vz =  tz;  break;
		case 7:   vx = -tx; vy = -ty; vz = -tz;  break;
		case 8:   vx =  tx; vy =  tz; vz =  ty;  break;
		case 9:   vx =  tx; vy = -tz; vz =  ty;  break;
		case 10:   vx =  tx; vy =  tz; vz = -ty;  break;
		case 11:   vx =  tx; vy = -tz; vz = -ty;  break;
		case 12:   vx = -tx; vy =  tz; vz =  ty;  break;
		case 13:   vx = -tx; vy = -tz; vz =  ty;  break;
		case 14:   vx = -tx; vy =  tz; vz = -ty;  break;
		case 15:   vx = -tx; vy = -tz; vz = -ty;  break;
		case 16:   vx =  ty; vy =  tx; vz =  tz;  break;
		case 17:   vx =  ty; vy =  tx; vz = -tz;  break;
		case 18:   vx = -ty; vy =  tx; vz =  tz;  break;
		case 19:   vx = -ty; vy =  tx; vz = -tz;  break;
		case 20:   vx =  ty; vy = -tx; vz =  tz;  break;
		case 21:   vx =  ty; vy = -tx; vz = -tz;  break;
		case 22:   vx = -ty; vy = -tx; vz =  tz;  break;
		case 23:   vx = -ty; vy = -tx; vz = -tz;  break;
		case 24:   vx =  tz; vy =  tx; vz =  ty;  break;
		case 25:   vx = -tz; vy =  tx; vz =  ty;  break;
		case 26:   vx =  tz; vy =  tx; vz = -ty;  break;
		case 27:   vx = -tz; vy =  tx; vz = -ty;  break;
		case 28:   vx =  tz; vy = -tx; vz =  ty;  break;
		case 29:   vx = -tz; vy = -tx; vz =  ty;  break;
		case 30:   vx =  tz; vy = -tx; vz = -ty;  break;
		case 31:   vx = -tz; vy = -tx; vz = -ty;  break;
		case 32:   vx =  ty; vy =  tz; vz =  tx;  break;
		case 33:   vx =  ty; vy = -tz; vz =  tx;  break;
		case 34:   vx = -ty; vy =  tz; vz =  tx;  break;
		case 35:   vx = -ty; vy = -tz; vz =  tx;  break;
		case 36:   vx =  ty; vy =  tz; vz = -tx;  break;
		case 37:   vx =  ty; vy = -tz; vz = -tx;  break;
		case 38:   vx = -ty; vy =  tz; vz = -tx;  break;
		case 39:   vx = -ty; vy = -tz; vz = -tx;  break;
		case 40:   vx =  tz; vy =  ty; vz =  tx;  break;
		case 41:   vx = -tz; vy =  ty; vz =  tx;  break;
		case 42:   vx =  tz; vy = -ty; vz =  tx;  break;
		case 43:   vx = -tz; vy = -ty; vz =  tx;  break;
		case 44:   vx =  tz; vy =  ty; vz = -tx;  break;
		case 45:   vx = -tz; vy =  ty; vz = -tx;  break;
		case 46:   vx =  tz; vy = -ty; vz = -tx;  break;
		case 47:   vx = -tz; vy = -ty; vz = -tx;  break;
	}
}

void getYawPitchRollRotationOrderFromAxisFrame( const triple &axisFrame, double &yaw, double &pitch, double &roll)
{
	if (matchDirection( UP, axisFrame.a))
	{
		if (UP != axisFrame.a)
		{
			pitch = -pitch;
			roll = -roll;
		}
		if (matchDirection(RIGHT, axisFrame.b))
		{
			if (RIGHT != axisFrame.b)
			{
				yaw = -yaw;
				roll = -roll;
			}
			if (FORWARD != axisFrame.c)
			{
				yaw = -yaw;
				pitch = -pitch;
			}
		}
		else	// matchDirection(RIGHT, axisFrame.c)
		{
			if (RIGHT != axisFrame.c)
			{
				yaw = -yaw;
				roll = -roll;
			}
			if (FORWARD != axisFrame.b)
			{
				yaw = -yaw;
				pitch = -pitch;
			}
		}
	}
	else if (matchDirection( UP, axisFrame.b))
	{
		if (UP != axisFrame.b)
		{
			pitch = -pitch;
			roll = -roll;
		}
		if (matchDirection(RIGHT, axisFrame.a))
		{
			if (RIGHT != axisFrame.a)
			{
				yaw = -yaw;
				roll = -roll;
			}
			if (FORWARD != axisFrame.c)
			{
				yaw = -yaw;
				pitch = -pitch;
			}
		}
		else	// matchDirection(RIGHT, axisFrame.c)
		{
			if (RIGHT != axisFrame.c)
			{
				yaw = -yaw;
				roll = -roll;
			}
			if (FORWARD != axisFrame.a)
			{
				yaw = -yaw;
				pitch = -pitch;
			}
		}
	}
	else // UP matches with axisFrame.c
	{
		if (UP != axisFrame.c)
		{
			pitch = -pitch;
			roll = -roll;
		}
		if (matchDirection(RIGHT, axisFrame.a))
		{
			if (RIGHT != axisFrame.a)
			{
				yaw = -yaw;
				roll = -roll;
			}
			if (FORWARD != axisFrame.b)
			{
				yaw = -yaw;
				pitch = -pitch;
			}
		}
		else	// matchDirection(RIGHT, axisFrame.b)
		{
			if (RIGHT != axisFrame.b)
			{
				yaw = -yaw;
				roll = -roll;
			}
			if (FORWARD != axisFrame.a)
			{
				yaw = -yaw;
				pitch = -pitch;
			}
		}
	}
}

void eulerCob( const triple &fromFrame, const triple &toFrame, double &yaw, double &pitch, double &roll );

int getEulerCaseNumber( const triple &from,	const triple &to)
{
	int eulerCaseNumber = 0;

	double signYaw(1.0), signPitch(1.0), signRoll(1.0);

	eulerCob( from, to, signYaw, signPitch, signRoll);

	if (signYaw < 0.0) eulerCaseNumber |= 0x04;
	if (signPitch < 0.0) eulerCaseNumber |= 0x02;
	if (signRoll < 0.0) eulerCaseNumber |= 0x01;

	return eulerCaseNumber;
}

void eulerCob( int eulerCaseNumber, double &yaw, double &pitch, double &roll)
{
	if (eulerCaseNumber & 0x04) yaw = -yaw;
	if (eulerCaseNumber & 0x02) pitch = -pitch;
	if (eulerCaseNumber & 0x01) roll = -roll;
}

void eulerCob( const triple &fromFrame, const triple &toFrame, double &yaw, double &pitch, double &roll )
{
	// Create normalized versions of the frames that can only contain Forward, Right, and Up
	triple newFromFrame( fromFrame.a & 3, fromFrame.b & 3, fromFrame.c & 3);
	triple newToFrame( toFrame.a & 3, toFrame.b & 3, toFrame.c & 3);

	// Convert to normalized frames by flipping axes if Back, Left, or Down are used.
	// If an axis is flipped, the sign of that rotation is kept but the other two reverse.
	getYawPitchRollRotationOrderFromAxisFrame( fromFrame, yaw, pitch, roll  );
	getYawPitchRollRotationOrderFromAxisFrame( toFrame, yaw, pitch, roll );

	// Now we are only left with six cases between the two normalized frames
	// The conversion between the normalized frames has no effect on yaw,pitch,roll with
	// the exception that if the conversion has a reflection in it then all the angles
	// change sign.
	int p = getCaseNumber( newFromFrame, newToFrame ) >> 3;
	if (p == 1 || p == 2 || p == 5)
	{
		// change between frames has a reflection
		yaw = -yaw;
		pitch = -pitch;
		roll = -roll;
	}
}

} // namespace cob
