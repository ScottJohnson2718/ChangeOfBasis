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

#pragma once

#ifndef COB_MATH
#define COB_MATH

#include <iostream>

// This is a small stub math library that is used for testing purposes.

namespace cob
{
	const double PI = 3.14159265358979643;
	const double PI_OVER_TWO =  0.5 * PI;
	const double TWO_PI = 2.0 * PI;

	const double RADIANS_PER_DEGREE = 0.017456492519946496;	
	const double DEGREES_PER_RADIAN = 57.295779513082641;

	class ColumnMatrix3d
	{
	public:

		explicit ColumnMatrix3d () {}
        
		explicit ColumnMatrix3d( double m00, double m01, double m02,
								 double m10, double m11, double m12,
								 double m20, double m21, double m22 );

		// These are so fundamental that they had to be put in.
		static const ColumnMatrix3d identity();
		static const ColumnMatrix3d fromRotationX( double x_radians );
		static const ColumnMatrix3d fromRotationY( double y_radians );
		static const ColumnMatrix3d fromRotationZ( double z_radians );

		double determinant() const;

		ColumnMatrix3d const  transpose() const;

		ColumnMatrix3d const  operator*( const ColumnMatrix3d &m) const; 
		bool equals (const ColumnMatrix3d &a, double tolerance = 0.001 ) const;
        
	public:
        
		// For rotation matrices...
		// Orthogonal basis vectors are stored in columns {m00, m10, m20}, { m01, m11, m21}, {m02, m12, m22}
		double _m[3][3];
	};

	std::ostream& operator << (std::ostream& os, const ColumnMatrix3d& m);

    class Quat4d
    {
    public:
        explicit Quat4d() {}

        explicit Quat4d( double x, double y, double z );
		static  Quat4d fromAxisAndAngle( double axisX, double axisY, double axisZ, double angleRads );

        explicit Quat4d( double w, double x, double y, double z )
        : _w(w), _x(x), _y(y), _z(z)
        {
        }

        Quat4d const  operator*( const Quat4d &m) const;
        Quat4d const  conjugate() const;
        double       norm() const;
        Quat4d const  negate() const;

        static Quat4d const  identity();
		bool equals( const Quat4d &a, double tolerance = 0.001 ) const;

    public:

        double   _w, _x, _y, _z;
    };

	std::ostream& operator << (std::ostream& os, const Quat4d& q);

    double  dot( Quat4d const &a, Quat4d const &b );

    class Vector3d
    {
    public:
        explicit Vector3d() {}

        explicit Vector3d( double x, double y, double z )
        : _x(x), _y(y), _z(z)
        {
        }

        Vector3d &     operator +=( const Vector3d &v );
        Vector3d &     operator -=( const Vector3d &v );

        Vector3d const operator + ( const Vector3d &v ) const;
        Vector3d const operator - ( const Vector3d &v ) const;
        Vector3d const operator-() const;
        Vector3d const operator * ( double f ) const;

        void        zero();

        double     magnitude() const;
        double     magnitudeSquared() const;

		bool		equals( const Vector3d &a, double tolerance ) const;

    public:

        double   _x, _y, _z;
    };

	std::ostream& operator << (std::ostream& os, const Vector3d& v);
    double  dot( Vector3d const &a, Vector3d const &b );
    Vector3d const cross( Vector3d const &a, Vector3d const &b );

    const ColumnMatrix3d    quat4dToColumnMatrix3d( const Quat4d &q );
    const Quat4d            columnMatrix3dToQuat4d( const ColumnMatrix3d &m );

}

#endif
