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
#pragma warning (push, 3)
#include "gtest\gtest.h"
#pragma warning(pop)

#include "CheckAgainstFullMath.h"

using namespace cob;

// These are not in a meaningful order
static triple allFrames[] =
{
	triple( FORWARD, RIGHT, UP ),
	triple( FORWARD, RIGHT, DOWN ),
	triple( FORWARD, UP, RIGHT ),
	triple( FORWARD, UP, LEFT ),
	triple( FORWARD, LEFT, UP ),
	triple( FORWARD, LEFT, DOWN ),
	triple( FORWARD, DOWN, RIGHT ),
	triple( FORWARD, DOWN, LEFT ),
	triple( RIGHT, FORWARD, UP ),
	triple( RIGHT, FORWARD, DOWN ),
	triple( RIGHT, UP, FORWARD ),
	triple( RIGHT, UP, BACK ),
	triple( RIGHT, BACK, UP ),
	triple( RIGHT, BACK, DOWN ),
	triple( RIGHT, DOWN, FORWARD ),
	triple( RIGHT, DOWN, BACK ),
	triple( UP, FORWARD, RIGHT ),
	triple( UP, FORWARD, LEFT ),
	triple( UP, RIGHT, FORWARD ),
	triple( UP, RIGHT, BACK ),
	triple( UP, BACK, RIGHT ),
	triple( UP, BACK, LEFT ),
	triple( UP, LEFT, FORWARD ),
	triple( UP, LEFT, BACK ),
	triple( BACK, RIGHT, UP ),
	triple( BACK, RIGHT, DOWN ),
	triple( BACK, UP, RIGHT ),
	triple( BACK, UP, LEFT ),
	triple( BACK, LEFT, UP ),
	triple( BACK, LEFT, DOWN ),
	triple( BACK, DOWN, RIGHT ),
	triple( BACK, DOWN, LEFT ),
	triple( LEFT, FORWARD, UP ),
	triple( LEFT, FORWARD, DOWN ),
	triple( LEFT, UP, FORWARD ),
	triple( LEFT, UP, BACK ),
	triple( LEFT, BACK, UP ),
	triple( LEFT, BACK, DOWN ),
	triple( LEFT, DOWN, FORWARD ),
	triple( LEFT, DOWN, BACK ),
	triple( DOWN, FORWARD, RIGHT ),
	triple( DOWN, FORWARD, LEFT ),
	triple( DOWN, RIGHT, FORWARD ),
	triple( DOWN, RIGHT, BACK ),
	triple( DOWN, BACK, RIGHT ),
	triple( DOWN, BACK, LEFT ),
	triple( DOWN, LEFT, FORWARD ),
	triple( DOWN, LEFT, BACK )
};

const triple getFrame( int index )
{
	return allFrames[index];
}

// This checks that the optimized version of the change of basis for a matrix
// or a quaternion matches the full math.
TEST(ChangeOfBasis, CheckMatrixAndQuatCob)
{
	int i, j;

	for (i = 0; i < 48; ++i)
	{
		triple from = getFrame(i);

		for (j = 0; j < 48; ++j)
		{
			triple to = getFrame(j);

			checkAgainstFullMath( from , to );
		}
	}
}


