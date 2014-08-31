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

#include "CheckAgainstFullMath.h"

#include "ChangeOfBasis.h"
#include "Math.h"

#pragma warning (push, 3)
#include "gtest\gtest.h"
#pragma warning(pop)

namespace cob
{
	void checkAgainstFullMath( triple const &from, triple const &to )
	{
		int caseNumber = getCaseNumber( from, to );

		ColumnMatrix3d mAtoB;

		getAtoBMatrix( caseNumber, 
			mAtoB._m[0][0], mAtoB._m[0][1], mAtoB._m[0][2],
			mAtoB._m[1][0], mAtoB._m[1][1], mAtoB._m[1][2],
			mAtoB._m[2][0], mAtoB._m[2][1], mAtoB._m[2][2] );

		ColumnMatrix3d mA(  900,  901,  902,
							910,  911,  912,
							920,  921,  922 );

		ColumnMatrix3d mB = mAtoB * mA * mAtoB.transpose();

		ColumnMatrix3d cob(mA);
		matrixCob3x3( caseNumber, 
			cob._m[0][0], cob._m[0][1], cob._m[0][2],
			cob._m[1][0], cob._m[1][1], cob._m[1][2],
			cob._m[2][0], cob._m[2][1], cob._m[2][2] );

		bool match = cob.equals( mB );

		if (!match)
		{
			std::cout << "Test failed ----" << std::endl;
			std::cout << "Correct ----" << std::endl;
			std::cout << mB << std::endl;

			std::cout << "Wrong " << std::endl;
			std::cout << cob << std::endl;
		}

		EXPECT_TRUE( match );

		// Make up some quaternion.  Construction will normalize it.
		Quat4d qA = Quat4d::fromAxisAndAngle( 1.0, 2.0, 3.0, RADIANS_PER_DEGREE * 23.0 );
	
		mA = quat4dToColumnMatrix3d( qA );
		mB = mAtoB * mA * mAtoB.transpose();
		Quat4d qAnswer = columnMatrix3dToQuat4d( mB );

		Quat4d qB(qA);

		quatCob( caseNumber, qB._x, qB._y, qB._z, qB._w );

		match = qAnswer.equals( qB );
		if (!match)
		{
			std::cout << "Test failed ----" << std::endl;
			std::cout << "Correct ----" << std::endl;
			std::cout << qAnswer << std::endl;

			std::cout << "Wrong " << std::endl;
			std::cout << qB << std::endl;
		}
		EXPECT_TRUE( match );
	}
}