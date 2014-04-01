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



  return cos(x)*sinh(y)*sin(2*PIE*t);
  

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

  return sin(x)*cosh(y)*sin(2*PIE*t);

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


return 0;


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


  gradu(0) = -sin(x)*sinh(y)*sin(2*PIE*t);
  gradu(1) = cos(x)*cosh(y)*sin(2*PIE*t);
  return gradu;


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

  gradv(0) = cos(x)*cosh(y)*sin(2*PIE*t);
  gradv(1) = sin(x)*sinh(y)*sin(2*PIE*t);
  return gradv;

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

#if THREED && !TIME
  gradw(0) = ((-1.0/3.0)*cos(2.0*PIE*x)*sin(2.0*PIE*y)*cos(2.0*PIE*z))/FAC;
  gradw(1) = ((-1.0/3.0)*sin(2.0*PIE*x)*cos(2.0*PIE*y)*cos(2.0*PIE*z))/FAC;
  gradw(2) = ((1.0/3.0)*sin(2.0*PIE*x)*sin(2.0*PIE*y)*sin(2.0*PIE*z))/FAC;
  return gradw;
#endif

#if THREED && TIME
  gradw(0) = ((-1.0/3.0)*cos(2.0*PIE*x)*sin(2.0*PIE*y)*cos(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
  gradw(1) = ((-1.0/3.0)*sin(2.0*PIE*x)*cos(2.0*PIE*y)*cos(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
  gradw(2) = ((1.0/3.0)*sin(2.0*PIE*x)*sin(2.0*PIE*y)*sin(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
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



  return cos(x)*sinh(y)*sin(2*PIE*t);

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

  return sin(x)*cosh(y)*sin(2*PIE*t);
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

  return 0;

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


  gradx(0) = -sin(x)*sinh(y)*sin(2*PIE*t);
  gradx(1) = cos(x)*cosh(y)*sin(2*PIE*t);
  return gradx;

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

  grady(0) = cos(x)*cosh(y)*sin(2*PIE*t);
  grady(1) = sin(x)*sinh(y)*sin(2*PIE*t);

  return grady;
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


#if THREED && !TIME
  gradz(0) = (-KPERM*4*PIE*PIE*cos(2*PIE*x)*sin(2*PIE*y)*cos(2*PIE*z))/FAC;
  gradz(1) = (-KPERM*4*PIE*PIE*sin(2*PIE*x)*cos(2*PIE*y)*cos(2*PIE*z))/FAC;
  gradz(2) = (KPERM*4*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z))/FAC;
  return gradz;
#endif


#if THREED && TIME
  gradz(0) = (-KPERM*4*PIE*PIE*cos(2*PIE*x)*sin(2*PIE*y)*cos(2*PIE*z)*sin(2*PIE*t))/FAC;
  gradz(1) = (-KPERM*4*PIE*PIE*sin(2*PIE*x)*cos(2*PIE*y)*cos(2*PIE*z)*sin(2*PIE*t))/FAC;
  gradz(2) = (KPERM*4*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z)*sin(2*PIE*t))/FAC;
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


	
  return (-sin(x)*sinh(y)-(cos(1)-1)*(cosh(1)-1))*sin(2*PIE*t);
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


  gradp(0) = -cos(x)*sinh(y)*sin(2*PIE*t);
  gradp(1) = -sin(x)*cosh(y)*sin(2*PIE*t);
  return gradp;


}

Number forcing_function_2D(const Point& p, const Parameters& parameters)
{
    const Real x = p(0);
    const Real y = p(1);
    const Real z = p(2);
   Real t =  parameters.get<Real>("time");
  const Real dt =  parameters.get<Real>("dt");

#if TTEST
t=FIX_T;
#endif



  return 0   ;
   
   


  }
  
  
  Number forcing_function_2D_u(const Point& p, const Parameters& parameters)
{
    const Real x = p(0);
    const Real y = p(1);
    const Real z = p(2);
   Real t =  parameters.get<Real>("time");
  const Real dt =  parameters.get<Real>("dt");

  
    return -cos(x)*sinh(y)*sin(2*PIE*t);
 

  }

  
   Number forcing_function_2D_v(const Point& p, const Parameters& parameters)
{
    const Real x = p(0);
    const Real y = p(1);
    const Real z = p(2);
   Real t =  parameters.get<Real>("time");
  const Real dt =  parameters.get<Real>("dt");

  
    return -sin(x)*cosh(y)*sin(2*PIE*t);
 

  }