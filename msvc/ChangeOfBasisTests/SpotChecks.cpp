
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
#include "Math.h"
#pragma warning (push, 3)
#include "gtest\gtest.h"
#pragma warning(pop)

#include "CheckAgainstFullMath.h"

using namespace cob;

TEST( QuatToMatrix, Test1)
{
	Quat4d qA = Quat4d::fromAxisAndAngle( 1.0, 2.0, 3.0, RADIANS_PER_DEGREE * 23.0 );
	
	ColumnMatrix3d mA = quat4dToColumnMatrix3d( qA );

	double d = mA.determinant();

	EXPECT_NEAR( d, 1.0, 0.001 );

	Quat4d qAnswer = columnMatrix3dToQuat4d( mA );
	
	bool match = qAnswer.equals( qA );
	if (!match)
	{
		std::cout << "QuatToMatrix Test failed ----" << std::endl;
		std::cout << "Correct ----" << std::endl;
		std::cout << qA << std::endl;

		std::cout << "Wrong " << std::endl;
		std::cout << qAnswer << std::endl;
	}
	EXPECT_TRUE( match );

}

TEST( ChangeOfBasis, SampleUseCaseUsingMatrix )
{
	// The transform from A to B is just a swap of the X and Y axes.
	// But the transform creates a reflection.
	triple InertialLabsFrame( RIGHT, FORWARD, UP);
	int caseNumber = getCaseNumber( InertialLabsFrame, Unreal3Frame );

	// Start with a rotation around X
	ColumnMatrix3d mA = ColumnMatrix3d::fromRotationX( RADIANS_PER_DEGREE * 33.0 );

	// Do the change of basis between the frames
	ColumnMatrix3d mB(mA);
	matrixCob3x3( caseNumber, 
		mB._m[0][0], mB._m[0][1], mB._m[0][2],
		mB._m[1][0], mB._m[1][1], mB._m[1][2],
		mB._m[2][0], mB._m[2][1], mB._m[2][2] );

	// The correct answer is a rotation around Y but the angle is negative because
	// we've switched from a right handed frame to a left handed frame
	ColumnMatrix3d mAnswer = ColumnMatrix3d::fromRotationY( RADIANS_PER_DEGREE * -33.0 );

	bool match = mAnswer.equals( mB );
	if (!match)
	{
		std::cout << "Test failed ----" << std::endl;
		std::cout << "Correct ----" << std::endl;
		std::cout << mAnswer << std::endl;

		std::cout << "Wrong " << std::endl;
		std::cout << mB << std::endl;
	}

	EXPECT_TRUE( match );

}

TEST(ChangeOfBasis, InertialLabsSensorsToUnreal)
{
	triple InertialLabsFrame( RIGHT, FORWARD, UP);

	int caseNumber = getCaseNumber( InertialLabsFrame, Unreal3Frame );

	EXPECT_EQ( caseNumber, 16 );

	// This is computed by hand by looking at a picture
	ColumnMatrix3d answerMAtoB(  0.0,  1.0,  0.0,
								 1.0,  0.0,  0.0,
								 0.0,  0.0,  1.0 );

	ColumnMatrix3d mAtoB;

	getAtoBMatrix( caseNumber, 
		mAtoB._m[0][0], mAtoB._m[0][1], mAtoB._m[0][2],
		mAtoB._m[1][0], mAtoB._m[1][1], mAtoB._m[1][2],
		mAtoB._m[2][0], mAtoB._m[2][1], mAtoB._m[2][2] );

	bool match = answerMAtoB.equals( mAtoB );
	EXPECT_TRUE( match );

	checkAgainstFullMath( InertialLabsFrame , Unreal3Frame );

	double yawIn( 23.0 * RADIANS_PER_DEGREE);
	double pitchIn(33.0 * RADIANS_PER_DEGREE);
	double rollIn( 80.0 * RADIANS_PER_DEGREE);

	double yaw(yawIn), pitch(pitchIn), roll(rollIn);

	int eulerCaseNumber = getEulerCaseNumber( InertialLabsFrame, Unreal3Frame);
	eulerCob( eulerCaseNumber, yaw, pitch, roll);

	EXPECT_DOUBLE_EQ( yaw, -yawIn );
	EXPECT_DOUBLE_EQ( pitch, -pitchIn );
	EXPECT_DOUBLE_EQ( roll, -rollIn );
}

TEST(ChangeOfBasis, KinectToOpenGL)
{
	int caseNumber = getCaseNumber( KinectFrame, OpenGLFrame );

	EXPECT_EQ( caseNumber, 5 );

	ColumnMatrix3d answerMAtoB( -1.0,  0.0,  0.0,
								 0.0,  1.0,  0.0,
								 0.0,  0.0, -1.0 );

	ColumnMatrix3d mAtoB;

	getAtoBMatrix( caseNumber, 
		mAtoB._m[0][0], mAtoB._m[0][1], mAtoB._m[0][2],
		mAtoB._m[1][0], mAtoB._m[1][1], mAtoB._m[1][2],
		mAtoB._m[2][0], mAtoB._m[2][1], mAtoB._m[2][2] );

	bool match = answerMAtoB.equals( mAtoB );
	EXPECT_TRUE( match );

	checkAgainstFullMath( KinectFrame , OpenGLFrame );

	double yawIn( 23.0 * RADIANS_PER_DEGREE);
	double pitchIn(33.0 * RADIANS_PER_DEGREE);
	double rollIn( 80.0 * RADIANS_PER_DEGREE);

	double yaw(yawIn), pitch(pitchIn), roll(rollIn);

	int eulerCaseNumber = getEulerCaseNumber( KinectFrame, OpenGLFrame);
	eulerCob( eulerCaseNumber, yaw, pitch, roll);

	EXPECT_DOUBLE_EQ( yaw, yawIn );
	EXPECT_DOUBLE_EQ( pitch, -pitchIn );
	EXPECT_DOUBLE_EQ( roll, -rollIn );

}