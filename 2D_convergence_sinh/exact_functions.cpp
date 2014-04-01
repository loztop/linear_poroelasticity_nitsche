#include "assemble.h"

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "exodusII_io.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "linear_implicit_system.h"
#include "transient_system.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"

// The definition of a geometric element
#include "elem.h"

#include "assemble.h"

#define FAC 1

#define POLY 0
#define SIN 1
#define PIE 3.1415926535897932384626433832

#define TTEST 0
#define FIX_T 1

Number exact_2D_solution_u(const Point& p,
                         const Parameters& parameters,  // parameters, not needed
                         const std::string&, // sys_name, not needed
                         const std::string&) // unk_name, not needed
{ 

  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  
   Real t =  parameters.get<Real>("time");


#if TIME && !THREED
  return cos(x)*sinh(y)*sin(t);
#endif

#if !TIME && !THREED
  return cos(x)*sinh(y);
#endif
  
	#if THREED
  return sin(t)*cos(x)*sinh(y)*sin(z);
#endif
	
}

Number exact_2D_solution_v(const Point& p,
                         const Parameters& parameters,  // parameters, not needed
                         const std::string&, // sys_name, not needed
                         const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");

#if TIME && !THREED
  return sin(x)*cosh(y)*sin(t);
#endif
  
  #if !TIME && !THREED
  return sin(x)*cosh(y);
#endif
	
	
#if THREED
  return sin(x)*cosh(y)*sin(z)*sin(t);
#endif	
	
  
}

Number exact_2D_solution_w(const Point& p,
                         const Parameters& parameters,   // parameters, not needed
                         const std::string&, // sys_name, not needed
                         const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");

#if !THREED
return 0;
#endif	

#if THREED
  return sin(x)*sinh(y)*cos(z)*sin(t);
#endif	
}


// We now define the gradient of the exact solution
Gradient exact_2D_derivative_u(const Point& p,
                             const Parameters& parameters,  // parameters, not needed
                             const std::string&, // sys_name, not needed
                             const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
   Real t =  parameters.get<Real>("time");

#if TTEST
t=FIX_T;
#endif

  // First derivatives to be returned.
  Gradient gradu;

#if TIME && !THREED
  gradu(0) = -sin(x)*sinh(y)*sin(t);
  gradu(1) = cos(x)*cosh(y)*sin(t);
  return gradu;
#endif

  #if !TIME && !THREED
  gradu(0) = -sin(x)*sinh(y);
  gradu(1) = cos(x)*cosh(y);
  return gradu;
#endif
	
	
	#if THREED
  gradu(0) = -sin(x)*sinh(y)*sin(z)*sin(t);
  gradu(1) = cos(x)*cosh(y)*sin(z)*sin(t);
  gradu(2) = cos(x)*sinh(y)*cos(z)*sin(t);
  return gradu;

#endif

}


// We now define the gradient of the exact solution
Gradient exact_2D_derivative_v(const Point& p,
                             const Parameters& parameters,  // parameters, not needed
                             const std::string&, // sys_name, not needed
                             const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
   Real t =  parameters.get<Real>("time");

#if TTEST
t=FIX_T;
#endif

  // First derivatives to be returned.
  Gradient gradv;

#if TIME && !THREED
  gradv(0) = cos(x)*cosh(y)*sin(t);
  gradv(1) = sin(x)*sinh(y)*sin(t);
  return gradv;
#endif
  
#if !TIME && !THREED
  gradv(0) = cos(x)*cosh(y) ;
  gradv(1) = sin(x)*sinh(y) ;
  return gradv;
#endif
  
	
		
#if THREED
  gradv(0) = cos(x)*cosh(y)*sin(z)*sin(t);
  gradv(1) = sin(x)*sinh(y)*sin(z)*sin(t);
  gradv(2) = sin(x)*cosh(y)*cos(z)*sin(t);
  return gradv;
#endif
	
}

// We now define the gradient of the exact solution
Gradient exact_2D_derivative_w(const Point& p,
                             const Parameters& parameters,   // parameters, not needed
                             const std::string&, // sys_name, not needed
                             const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");

  // First derivatives to be returned.
  Gradient gradw;

	#if !THREED 
	  return gradw;
#endif
	
#if THREED 
  gradw(0) = cos(x)*sinh(y)*cos(z)*sin(t);
  gradw(1) = sin(x)*cosh(y)*cos(z)*sin(t);
  gradw(2) = -sin(x)*sinh(y)*sin(z)*sin(t);
  return gradw;
#endif
	
}



Number exact_2D_solution_x(const Point& p,
                         const Parameters& parameters,  // parameters, not needed
                         const std::string&, // sys_name, not needed
                         const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
   Real t =  parameters.get<Real>("time");


#if TIME && !THREED
  return cos(x)*sinh(y)*sin(t);
#endif
  
 #if !TIME && !THREED
  return cos(x)*sinh(y);
#endif 

	#if THREED
  return sin(t)*cos(x)*sinh(y)*sin(z);
#endif
	
}

Number exact_2D_solution_y(const Point& p ,
                         const Parameters& parameters,  // parameters, not needed
                         const std::string&, // sys_name, not needed
                         const std::string&) // unk_name, not needed
{

  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");

#if TIME && !THREED
  return sin(x)*cosh(y)*sin(t);
#endif
  
  #if !TIME && !THREED
  return sin(x)*cosh(y);
#endif
	
	#if THREED
  return sin(x)*cosh(y)*sin(z)*sin(t);
#endif	
	
}

Number exact_2D_solution_z(const Point& p,
                         const Parameters& parameters,   // parameters, not needed
                         const std::string&, // sys_name, not needed
                         const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");

  // analytic solution value

#if !THREED
  return 0;
#endif
	
	#if THREED
  return sin(x)*sinh(y)*cos(z)*sin(t);
#endif	
	
}


// We now define the gradient of the exact solution
Gradient exact_2D_derivative_x(const Point& p,
                             const Parameters& parameters,  // parameters, not needed
                             const std::string&, // sys_name, not needed
                             const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
   Real t =  parameters.get<Real>("time");

#if TTEST
t=FIX_T;
#endif

  // First derivatives to be returned.
  Gradient gradx;

#if TIME && !THREED
  gradx(0) = -sin(x)*sinh(y)*sin(t);
  gradx(1) = cos(x)*cosh(y)*sin(t);
  return gradx;
#endif
  
#if !TIME && !THREED
  gradx(0) = -sin(x)*sinh(y) ;
  gradx(1) = cos(x)*cosh(y) ;
  return gradx;
#endif
	
	
	#if THREED
  gradx(0) = -sin(x)*sinh(y)*sin(z)*sin(t);
  gradx(1) = cos(x)*cosh(y)*sin(z)*sin(t);
  gradx(2) = cos(x)*sinh(y)*cos(z)*sin(t);
	  return gradx;
#endif
	
}


// We now define the gradient of the exact solution
Gradient exact_2D_derivative_y(const Point& p,
                             const Parameters& parameters,  // parameters, not needed
                             const std::string&, // sys_name, not needed
                             const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");

#if TTEST
t=FIX_T;
#endif

  // First derivatives to be returned.
  Gradient grady;

  
#if TIME && !THREED
  grady(0) = cos(x)*cosh(y)*sin(t);
  grady(1) = sin(x)*sinh(y)*sin(t);

  return grady;
#endif
  
  
  #if !TIME && !THREED
  grady(0) = cos(x)*cosh(y);
  grady(1) = sin(x)*sinh(y);

  return grady;
#endif
	
	
		
	#if THREED
  grady(0) = cos(x)*cosh(y)*sin(z)*sin(t);
  grady(1) = sin(x)*sinh(y)*sin(z)*sin(t);
  grady(2) = sin(x)*cosh(y)*cos(z)*sin(t);
  return grady;
#endif
	
  
  
}



// We now define the gradient of the exact solution
Gradient exact_2D_derivative_z(const Point& p,
                             const Parameters& parameters,   // parameters, not needed
                             const std::string&, // sys_name, not needed
                             const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");

  // First derivatives to be returned.
  Gradient gradz;
#if !THREED 
  return gradz;
#endif
#if THREED 
  gradz(0) = cos(x)*sinh(y)*cos(z)*sin(t);
  gradz(1) = sin(x)*cosh(y)*cos(z)*sin(t);
  gradz(2) = -sin(x)*sinh(y)*sin(z)*sin(t);
  return gradz;
#endif
	
}

Number exact_2D_solution_p(const Point& p,
                         const Parameters& parameters,  // parameters, not needed
                         const std::string&, // sys_name, not needed
                         const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");


#if TIME && !THREED
  return (-sin(x)*sinh(y))*sin(t);
#endif
  
 #if !TIME && !THREED
  return (-sin(x)*sinh(y));
#endif
   
	#if THREED
  return (-sin(x)*sinh(y)*sin(z))*sin(t);
#endif
  
}


// We now define the gradient of the exact solution
Gradient exact_2D_derivative_p(const Point& p,
                         const Parameters& parameters,  // parameters, not needed
                             const std::string&, // sys_name, not needed
                             const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
   Real t =  parameters.get<Real>("time");

  // First derivatives to be returned.
  Gradient gradp;


  gradp(0) = -cos(x)*sinh(y)*sin(t);
  gradp(1) = -sin(x)*cosh(y)*sin(t);
  return gradp;

	

}

Number forcing_function_2D(const Point& p, const Parameters& parameters)
{
    const Real x = p(0);
    const Real y = p(1);
    const Real z = p(2);
   Real t =  parameters.get<Real>("time");
  const Real dt =  parameters.get<Real>("dt");

#if !THREED
  return 0   ;
  #endif
 
   
#if THREED
	  return -(cos(t)+sin(t))*(sin(x)*sinh(y)*sin(z)) ;
#endif
   

  }
  
  
  Number forcing_function_2D_u(const Point& p, const Parameters& parameters)
{
    const Real x = p(0);
    const Real y = p(1);
    const Real z = p(2);
   Real t =  parameters.get<Real>("time");
  const Real dt =  parameters.get<Real>("dt");

  
#if TIME && !THREED
    return -cos(x)*sinh(y)*sin(t);
#endif

	
#if !TIME && !THREED
    return -cos(x)*sinh(y);
#endif
		
		#if THREED
		return 0;
#endif
		
	
  }

  
   Number forcing_function_2D_v(const Point& p, const Parameters& parameters)
{
    const Real x = p(0);
    const Real y = p(1);
    const Real z = p(2);
   Real t =  parameters.get<Real>("time");
  const Real dt =  parameters.get<Real>("dt");

#if TIME && !THREED  
    return -sin(x)*cosh(y)*sin(t);
#endif
	

#if !TIME && !THREED  
    return -sin(x)*cosh(y);
#endif
		
		#if THREED
		return 0;
#endif
		
	
  }