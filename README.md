# Change of Basis
Supopse you have an orientation, say from 3D Studio Max, and you need to use it in a game engine like Unreal. The orientation could be a matrix, a quaternion or Euler angles. It is not straightforward how you generically can take an orientation in one frame and convert to another frame. The provided code makes this simple without knowing any of the math. The code is designed so that you can add the Change of Basis to an existing math library, as long as you are familiar with the conventions of the math library such as whether it uses column matrices or row matrices.

The code has direct applications when doing work with sensors and Virtual Reality. IMU sensors all report their orientations in space using quaternions. You need to take the quaternion from the sensor and put it into the frame that your application uses. For instance, I frequently use sensors from a vendor that uses X = Right, Y = Forward, Z = Up. I have to use those rotations in Unreal 3 where the global frame is X = Forward, Y = Right, Z = Up. I also have to use the same rotations in an OpenGL application that uses X = Left, Y = Up, Z = Forward. The attached code makes this so simple that you will never think about it again. The code also works without any floating point operations that would change the components of the rotation other than permuting them or changing their sign. This goes with the notion that changing the frame shouldn't change the rotation at all.

The part of the problem that I cannot help solve is that too much code uses reference frames implicitly. In those cases, you will have to reverse engineer what the axes are. When you are lucky, the reference frame you need is explicitly printed in manuals or even shown as an axis frame on the screen with the X axis in red, the Y axis in green and the Z axis in blue. Or you will find that the reference frame convention is not well defined as is the case with the BVH animation format. Some applications use Y = Up and some applications use Z = Up.

## How the Code was Written
The code was mostly generated. No one could write out that many cases without going nuts. The matrix change of basis was derived as in the the attached document. The quaternion version comes from symbolically changing the quaternion to a matrix, performing the matrix change of basis, and then converting back to a matrix. The operation was written for a symbolic math tool called Maxima. Maxima simplified the results into a pattern that I could follow. Then I wrote a code generator to create that pattern. The same method was performed for the Euler change of basis.

## Using the Code
The use case for using the code intends to be very simple. Determine the two reference frames that you need to convert between and get a case number from them. If you see the code, you'll see why the variable is called 'caseNumber'. Then pass in the caseNumber to the function that does the change of basis, either quatCob() or matrixCob3x3().

The code is meant to be added to your exsting code as source. Just put the source files in with your other code and compile.

## Advantage Versus Doing the Math Yourself
This code only permutes and negates your input numbers. The actual math is two matrix multiplies that add round-off error to your numbers. This code keeps the exact precision of your input data, guaranteed.
## Samples
This sample code converts a quaternion from a vendor's frame (RIGHT, FORWARD,UP) to (FORWARD, RIGHT, UP)

```
using namespace cob;

triple sensorVendor( RIGHT, FORWARD, UP );
triple unreal3( FORWARD, RIGHT, UP );

// This sets up the transform from the sensor vendor's frame to the Unreal 3 frame
int caseNumber = getCaseNumber( sensorVendor, unreal3 );

// assume that there is a Quaternion class whose components are qx, qy, qz, and qw;
Quaternion sensor = getSensorRotationFromSensorHardware();

// This takes the sensor quaternion in the sensor vendor's frame and changes it to the Unreal 3 frame
quatCob( caseNumber, sensor.qx, sensor.qy, sensor.qz, sensor.qw );
```
Here is the same example only done with matrices.
```
using namespace cob;

triple sensorVendor( RIGHT, FORWARD, UP );
triple unreal3( FORWARD, RIGHT, UP );

// This sets up the transform from the sensor vendor's frame to the Unreal 3 frame
int caseNumber = getCaseNumber( sensorVendor, unreal3 );

// assume that there is a 3x3 matrix class whose members are
// m00, m01, m02,
// m10, m11, m12,
// m20, m21, m22
Matrix sensor = getSensorRotationFromSensorHardwareAsMatrix3x3();

// This takes the sensor matrix in the sensor vendor's frame and changes it to the Unreal 3 frame
matrixCob3x3( caseNumber,
    sensor.m00, sensor.m01, sensor.m02,
    sensor.m10, sensor.m11, sensor.m12,
    sensor.m20, sensor.m21, sensor.m22);
```
## Using the code
Just cut/paste changeOfBasis.cpp and changeOfBasis.h into your project. The Visual Studio project is just for the test code and you don't need it.
