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

// Note : The Apache 2.0 license allows free personal and commercial use.

#pragma once

#ifndef CHANGEOFBASIS_H
#define	CHANGEOFBASIS_H

// A Change of Basis is converting a quantity to another reference frame.
// Changing a vector to another reference frame is straightforward but
// changing a rotation to another frame is not as easy and is poorly
// understood by most programmers.  This intends to be a gem that solves
// the change of basis.  It is meant as an add-on to an existing math
// library.  It allows you to change information say, from Maya to Unreal,
// or from a BVH animation file to OpenGL.

namespace cob
{
	const int FORWARD = 0;
	const int RIGHT = 1;
	const int UP = 2;
	const int BACK = 4;		// BACK & 3 == FORWARD
	const int LEFT = 5;		// LEFT & 3 == RIGHT
	const int DOWN = 6;		// DOWN & 3 == UP

	struct triple
	{
		triple(int aa, int bb, int cc)
			: a(aa), b(bb), c(cc)
		{
		}
		int a, b, c;
	};

	// These are sample reference frames.  They list the X, Y, and Z axis from 
	// the point of view of a character.  For instance if a character is looking foward along
	// the positive X axis, Y is to his right and Z is up then the reference frame would
	// be (Forward, Right, Up).
	const triple Unreal3Frame( FORWARD, RIGHT, UP );
	const triple OpenGLFrame( LEFT, UP, FORWARD );
	const triple OculusFrame( RIGHT, UP, BACK );
	const triple BvhFrame( LEFT, UP, FORWARD );
	const triple BvhBlenderFrame( LEFT, FORWARD, UP );
	const triple KinectFrame( RIGHT, UP, BACK );
	const triple PrioVRFrame( RIGHT, UP, FORWARD );

	// To use the Change of Basis on a matrix, vector or a quaternion, call this first
	// with the two frames ("from" and "to").  
	int getCaseNumber( 
		const triple &from,	// make a triple from a permutation of (FORWARD or BACK, LEFT or RIGHT, UP or DOWN )
		const triple &to );	// make a triple from a permutation of (FORWARD or BACK, LEFT or RIGHT, UP or DOWN )

	// Matrix Change of Basis
	// The caseNumber is described above and represents the matrix MAtoB.
	// The matrix MA passed in must be a column matrix of doubles 
	// (ex. The first basis vector is the column 
	// [ MA00 ]
	// [ MA10 ]
	// [ MA20 ]  )
	//
	// Let MA be a rotation in the A basis frame (the "from" frame)
	// Let MB be a rotation in the B basis frame (the "to" frame)
	// To change MA to a rotation MB in the B basis, this function performs this math:
	//
	//  [ MB ] = [ MAtoB ] . [ MA ] . transpose([ MAtoB ])
	//
	void matrixCob3x3(int caseNumber, 
		double &MA00, double &MA01, double &MA02,
		double &MA10, double &MA11, double &MA12,
		double &MA20, double &MA21, double &MA22);

	// Vector Change of Basis
	// Performs [ VB ] = [ MAtoB ] . [ VA ]
	void vectorCob( int caseNumber, double &vAx, double &vAy, double &vAz );

	// Quaternion Change of Basis
	// The caseNumber is described above and represents the matrix MAtoB.
	// This function performs a change of basis on a quaternion by very efficiently
	// changing the quaternion to a matrix, doing a change of basis on the matrix
	// and then converting it back to a quaternion.  That function simplifies to
	// permuting the quaternion components qx, qy, and qz and changing their signs.
	// The calculation is done such that qw does not change.  qw doesn't actually
	// need to be a parameter.
	void quatCob( int caseNumber, double &qx, double &qy, double &qz, double &qw );

	// Euler Angle (yaw, pitch, roll) Change of Basis
	// The Euler Case Number is a piece of data that makes it faster to perform a Change of Basis
	// on multiple sets of euler angles using the same "from" and "to".
	int getEulerCaseNumber( 
		const triple &from,				// make a triple from a permutation of (FORWARD or BACK, LEFT or RIGHT, UP or DOWN )
		const triple &to);				// make a triple from a permutation of (FORWARD or BACK, LEFT or RIGHT, UP or DOWN )

	// Euler Angle (yaw, pitch, roll) Change of Basis
	// The rotation order is Yaw around the Up/Down axis, Pitch around left/right then Roll around Forward/Back
	// The sign convention is based on whether the frames you specify are right or left handed.
	// The eulerCaseNumber is described above.
	// This operation can only change the signs of the yaw, pitch, and/or roll.
	// This function efficiently converts the euler angles to a matrix, then
	// performs a change of basis on the matrix and then converts the matrix back
	// to euler angles.  The net result is just potential change of signs.
	// Since it only changes their signs, this function works on radians and degrees.
	void eulerCob( int eulerCaseNumber, double &yaw, double &pitch, double &roll );

	// Euler Angle (yaw, pitch, roll) Change of Basis
	// Same as above but is slower across multiple calls since it has to the do work
	// to compute the case number over and over again.
	// Since it only changes their signs, this function works on radians and degrees.
	void eulerCob( const triple &from, const triple &to, double &yaw, double &pitch, double &roll );


	// This is for testing or for showing customers what is going on under the hood.
	// You provide the members of a 3x3 column vector and a caseNumber and this sets the matrix
	// elements to mAtoB as mentioned above.
	void getAtoBMatrix(int caseNumber, 
		double &m00, double &m01, double &m02,
		double &m10, double &m11, double &m12,
		double &m20, double &m21, double &m22);

}

#endif // CHANGEOFBASIS_H
