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

#include "Math.h"
#include <math.h>
//#include <iostream>

namespace cob
{

ColumnMatrix3d::ColumnMatrix3d( double m00, double m01, double m02,
								 double m10, double m11, double m12,
								 double m20, double m21, double m22 )
{
	_m[0][0] = m00;	_m[0][1] = m01;	_m[0][2] = m02;
	_m[1][0] = m10;	_m[1][1] = m11;	_m[1][2] = m12;
	_m[2][0] = m20;	_m[2][1] = m21;	_m[2][2] = m22;
}

ColumnMatrix3d const ColumnMatrix3d::operator*( const ColumnMatrix3d &a ) const
{
    // use the Return value optimization from Scott Meyer's Effective C++
    // Compilers can optimize out the temporaries objects created here if a constructor is returned.

    return ColumnMatrix3d( 
        _m[0][0] * a._m[0][0] +_m[0][1] * a._m[1][0] +_m[0][2] * a._m[2][0],
        _m[0][0] * a._m[0][1] +_m[0][1] * a._m[1][1] +_m[0][2] * a._m[2][1],
        _m[0][0] * a._m[0][2] +_m[0][1] * a._m[1][2] +_m[0][2] * a._m[2][2],

        _m[1][0] * a._m[0][0] +_m[1][1] * a._m[1][0] +_m[1][2] * a._m[2][0],
        _m[1][0] * a._m[0][1] +_m[1][1] * a._m[1][1] +_m[1][2] * a._m[2][1],
        _m[1][0] * a._m[0][2] +_m[1][1] * a._m[1][2] +_m[1][2] * a._m[2][2],

        _m[2][0] * a._m[0][0] +_m[2][1] * a._m[1][0] +_m[2][2] * a._m[2][0],
        _m[2][0] * a._m[0][1] +_m[2][1] * a._m[1][1] +_m[2][2] * a._m[2][1],
        _m[2][0] * a._m[0][2] +_m[2][1] * a._m[1][2] +_m[2][2] * a._m[2][2] );
}

const ColumnMatrix3d ColumnMatrix3d::fromRotationX( double _xradians )
{
    double c = cos( _xradians );
    double s = sin( _xradians );

    return ColumnMatrix3d( 1.0f, 0.0f, 0.0f,
                           0.0f,    c,   -s,
                           0.0f,    s,    c);
}

const ColumnMatrix3d ColumnMatrix3d::fromRotationY( double _yradians )
{
    double c = cos( _yradians );
    double s = sin( _yradians );

    return ColumnMatrix3d(      c, 0.0f,    s,
                             0.0f, 1.0f, 0.0f,
                               -s, 0.0f,    c);

}

const ColumnMatrix3d ColumnMatrix3d::fromRotationZ( double _zradians )
{
    double c = cos( _zradians );
    double s = sin( _zradians );

    return ColumnMatrix3d(    c,   -s, 0.0f, 
                              s,    c, 0.0f,
                           0.0f, 0.0f, 1.0f );
}

const ColumnMatrix3d ColumnMatrix3d::transpose() const
{
    return ColumnMatrix3d(_m[0][0],_m[1][0],_m[2][0],
                          _m[0][1],_m[1][1],_m[2][1],
                          _m[0][2],_m[1][2],_m[2][2] );
}

const ColumnMatrix3d ColumnMatrix3d::identity()
{
    return ColumnMatrix3d( 1.0f, 0.0f, 0.0f,
                           0.0f, 1.0f, 0.0f,
                           0.0f, 0.0f, 1.0f );
}

bool ColumnMatrix3d::equals(const ColumnMatrix3d &a, double tolerance ) const
{
	int p, q;
	for (p = 0; p < 3; ++p)
	{
		for (q = 0; q < 3; ++q)
		{
			if (fabs( _m[p][q] - a._m[p][q]) > tolerance)
			{
				return false;
			}
		}
	}
	return true;
}

double ColumnMatrix3d::determinant() const
{
	double c0 = _m[1][1] * _m[2][2] - _m[2][1] * _m[1][2];
	double c1 = _m[1][0] * _m[2][2] - _m[2][0] * _m[1][2];
	double c2 = _m[1][0] * _m[2][1] - _m[2][0] * _m[1][1];
	double d = _m[0][0] * c0 - _m[0][1] * c1 + _m[0][2] * c2;
	return d;
}


std::ostream& operator << (std::ostream& os, const ColumnMatrix3d& m)
{
    return os << std::endl << m._m[0][0] << "\t" << m._m[0][1] << "\t" << m._m[0][2] 
              << std::endl << m._m[1][0] << "\t" << m._m[1][1] << "\t" << m._m[1][2] 
              << std::endl << m._m[2][0] << "\t" << m._m[2][1] << "\t" << m._m[2][2] 
              << std::endl;
}

//////////////////

Quat4d::Quat4d( double x, double y, double z)
: _x(x), _y(y), _z(z)
{
    _w = 1.0f - sqrt( _x * _x + _y * _y + _z * _z );
    if (_w < 0.0f)
    {
        _w = 0.0f;
    }
    else
    {
        _w = sqrt(_w);
    }
}

Quat4d Quat4d::fromAxisAndAngle( double axisX, double axisY, double axisZ, double angleRads )
{
	Quat4d q;
    double d, s, halfAngle;
	Vector3d axis( axisX, axisY, axisZ );

    d = axis.magnitude();
    if (fabs(d - 0.0) < 0.001)
    {
    	q._x = 0.0; q._y = 0.0; q._z = 0.0; q._w = 1.0;
    	return q;
    }

    halfAngle = 0.5f * angleRads;
    s = sin(halfAngle) / d;
    q._x = s * axis._x;
    q._y = s * axis._y;
    q._z = s * axis._z;
    q._w = cos(halfAngle);

	return q;
}

double  dot( Quat4d const &a, Quat4d const &b )
{
    return (a._w * b._w + a._x + b._x + a._y * b._y + a._z + b._z );
}

// static
Quat4d const    Quat4d::identity()
{
    return Quat4d( 1.0f, 0.0f, 0.0f, 0.0f );
}

Quat4d const    Quat4d::conjugate() const
{
    return Quat4d( _w, -_x, -_y, -_z );
}

double         Quat4d::norm() const
{
    return sqrt( _w * _w + _x * _x + _y * _y + _z * _z );
}

Quat4d const    Quat4d::operator *( Quat4d const & b ) const
{
    return Quat4d( _w * b._w - _x * b._x - _y * b._y - _z * b._z,
                   _w * b._x + _x * b._w + _y * b._z - _z * b._y,
                   _w * b._y - _x * b._z + _y * b._w + _z * b._x,
                   _w * b._z + _x * b._y - _y * b._x + _z * b._w );
}

Quat4d const    Quat4d::negate() const
{
    return Quat4d( -_w, -_x, -_y, -_z );
}

bool Quat4d::equals( const Quat4d &a, double tolerance ) const
{
	if (fabs( a._w - _w) > tolerance)
	{
		return false;
	}
	if (fabs( a._x - _x) > tolerance)
	{
		return false;
	}
	if (fabs( a._y - _y) > tolerance)
	{
		return false;
	}
	if (fabs( a._z - _z) > tolerance)
	{
		return false;
	}
	return true;
}

std::ostream& operator << (std::ostream& os, const Quat4d& q)
{
    return os << "( x:" << q._x << ", y:" << q._y << ", z:" << q._z << ", w:" << q._w << ')';
}

////////////////////////////

Vector3d &     Vector3d::operator +=( const Vector3d &v )
{
    _x += v._x;
    _y += v._y;
    _z += v._z;
    return (*this);
}

Vector3d &     Vector3d::operator -=( const Vector3d &v )
{
    _x -= v._x;
    _y -= v._y;
    _z -= v._z;
    return (*this);
}

Vector3d const Vector3d::operator + ( const Vector3d &v ) const
{
    return Vector3d( _x + v._x, _y + v._y, _z + v._z );
}
Vector3d const Vector3d::operator - ( const Vector3d &v ) const
{
    return Vector3d( _x - v._x, _y - v._y, _z - v._z );

}
Vector3d const Vector3d::operator-() const
{
    return Vector3d( -_x , -_y, -_z );
}

double  dot( Vector3d const &a, Vector3d const &b )
{
    return (a._x + b._x + a._y * b._y + a._z + b._z );
}

Vector3d const cross( Vector3d const &a, Vector3d const &b )
{
    return Vector3d( -a._z * b._y + a._y * b._z, a._z * b._x - a._x * b._z, -a._y * b._x + a._x * b._y );
}

Vector3d const Vector3d::operator*( double f ) const
{
    return Vector3d( _x * f, _y * f, _z * f );
}

void    Vector3d::zero()
{
    _x = 0.0;
    _y = 0.0;
    _z = 0.0;
}

double Vector3d::magnitudeSquared() const
{
    return ( _x * _x + _y * _y + _z * _z );
}

double Vector3d::magnitude() const
{
    return sqrt( _x * _x + _y * _y + _z * _z  );
}


bool Vector3d::equals( const Vector3d &a, double tolerance ) const
{
	if (fabs( a._x - _x) > tolerance)
	{
		return false;
	}
	if (fabs( a._y - _y) > tolerance)
	{
		return false;
	}
	if (fabs( a._z - _z) > tolerance)
	{
		return false;
	}
	return true;
}

std::ostream& operator << (std::ostream& os, const Vector3d& v)
{
    return os << '(' << v._x << ", " << v._y << ", " << v._z << ')';
}

//////////////////////////////////////////////

const Vector3d    operator *( ColumnMatrix3d const &a, Vector3d const &v )
{
    return Vector3d( a._m[0][0] * v._x + a._m[0][1] * v._y + a._m[0][2] * v._z,
                     a._m[1][0] * v._x + a._m[1][1] * v._y + a._m[1][2] * v._z,
                     a._m[2][0] * v._x + a._m[2][1] * v._y + a._m[2][2] * v._z );
}


const ColumnMatrix3d outerProduct( const Vector3d &v )
{
    double xy = v._x * v._y;
    double yz = v._y * v._z;
    double xz = v._x * v._z;

    return ColumnMatrix3d( v._x * v._x,          xy,          xz, 
                           xy, v._y * v._y,          yz, 
                           xz,          yz, v._z * v._z );
}

const ColumnMatrix3d skewSymmetricMatrix( const Vector3d &v )
{
    return ColumnMatrix3d( 0.0f, -v._z,  v._y,
                           v._z,  0.0f, -v._x,
                           -v._y,  v._x, 0.0f );
}


// page 126 of Quaternions and Rotation Sequences by Kuipers
const ColumnMatrix3d quat4dToColumnMatrix3d( const Quat4d &q )
{
    return ColumnMatrix3d( 
		2.0 * q._w * q._w - 1.0 + 2.0 * q._x * q._x,  2.0 * (q._x * q._y - q._w * q._z) ,            2.0 * (q._x * q._z + q._w * q._y ),
        2.0 * (q._x * q._y + q._w * q._z),              2.0 * q._w * q._w - 1.0 + 2.0 * q._y * q._y, 2.0 * (q._y * q._z - q._w * q._x ),
        2.0 * (q._x * q._z - q._w * q._y ) ,            2.0 * (q._y * q._z + q._w * q._x ),              2.0 * q._w * q._w - 1.0 + 2.0 * q._z * q._z);
}

// http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche52.html
inline double SIGN(double x) {return (x >= 0.0) ? 1.0 : -1.0;}
inline double NORM(double a, double b, double c, double d) {return sqrt(a * a + b * b + c * c + d * d);}

const Quat4d    columnMatrix3dToQuat4d( const ColumnMatrix3d &a )
{
	double q0 = ( a._m[0][0] + a._m[1][1] + a._m[2][2] + 1.0f) / 4.0;
	double q1 = ( a._m[0][0] - a._m[1][1] - a._m[2][2] + 1.0f) / 4.0;
	double q2 = (-a._m[0][0] + a._m[1][1] - a._m[2][2] + 1.0f) / 4.0;
	double q3 = (-a._m[0][0] - a._m[1][1] + a._m[2][2] + 1.0f) / 4.0;

	if(q0 < 0.0) q0 = 0.0f;
	if(q1 < 0.0) q1 = 0.0f;
	if(q2 < 0.0) q2 = 0.0f;
	if(q3 < 0.0) q3 = 0.0f;

	q0 = sqrt(q0);
	q1 = sqrt(q1);
	q2 = sqrt(q2);
	q3 = sqrt(q3);

	if(q0 >= q1 && q0 >= q2 && q0 >= q3) 
	{
		q0 *= 1.0;
		q1 *= SIGN(a._m[2][1] - a._m[1][2]);
		q2 *= SIGN(a._m[0][2] - a._m[2][0]);
		q3 *= SIGN(a._m[1][0] - a._m[0][1]);
	} else if(q1 >= q0 && q1 >= q2 && q1 >= q3) 
	{
		q0 *= SIGN(a._m[2][1] - a._m[1][2]);
		q1 *= 1.0;
		q2 *= SIGN(a._m[1][0] + a._m[0][1]);
		q3 *= SIGN(a._m[0][2] + a._m[2][0]);
	} else if(q2 >= q0 && q2 >= q1 && q2 >= q3) 
	{
		q0 *= SIGN(a._m[0][2] - a._m[2][0]);
		q1 *= SIGN(a._m[1][0] + a._m[0][1]);
		q2 *= 1.0f;
		q3 *= SIGN(a._m[2][1] + a._m[1][2]);
	} else if(q3 >= q0 && q3 >= q1 && q3 >= q2) 
	{
		q0 *= SIGN(a._m[1][0] - a._m[0][1]);
		q1 *= SIGN(a._m[2][0] + a._m[0][2]);
		q2 *= SIGN(a._m[2][1] + a._m[1][2]);
		q3 *= 1.0;
	} else 
	{
		// this is bad
		//int foo = 4;
	}
	double r = NORM(q0, q1, q2, q3);
	q0 /= r;
	q1 /= r;
	q2 /= r;
	q3 /= r;

    return Quat4d( q0, q1, q2, q3 );
}


} // namespace cob