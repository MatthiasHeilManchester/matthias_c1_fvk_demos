//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC//
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC//
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC//
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC//
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC//
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC//
//LIC//====================================================================
#include <fenv.h>

#include <random>

//Generic routines
#include "generic.h"

// The mesh
#include "meshes/triangle_mesh.h"

// The equations
#include "c1_foeppl_von_karman.h"


using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;



// Random number between 0 and 1
namespace Random
{
 double random_between_zero_and_one()
 {
  static std::mt19937 gen(12345); // note this is static so will only be executed once! 
  static std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(gen);
 }
 
}


// hierher move to geom_objects.h

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
// Straight line as geometric object
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


 //=========================================================================
 /// Steady, straight 1D line in 2D space connecting two specified points.
 /// First point reached for zeta = 0; last one for zeta = 1 
 //=========================================================================
class TwoDStraightLineFromTwoPoints : public GeomObject
{
public:
 
 /// Constructor: Pass left and right point.
 TwoDStraightLineFromTwoPoints(const Vector<double>& left,
                               const Vector<double>& right,
                               const double& phi=0.0) 
  : GeomObject(1, 2)
  {
#ifdef PARANOID
   if (left.size() != 2)
    {
     std::ostringstream error_message;
     error_message << "left point should have size 2, not "
                   << left.size() << std::endl;
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   if (right.size() != 2)
    {
     std::ostringstream error_message;
     error_message << "right point should have size 2, not "
                   << right.size() << std::endl;
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

   Left.resize(2);
   Left[0]=left[0]*cos(phi)-left[1]*sin(phi);
   Left[1]=left[1]*cos(phi)+left[0]*sin(phi);
   Right.resize(2);
   Right[0]=right[0]*cos(phi)-right[1]*sin(phi);
   Right[1]=right[1]*cos(phi)+right[0]*sin(phi);
   
  }
 
 /// Broken copy constructor
 TwoDStraightLineFromTwoPoints(const TwoDStraightLineFromTwoPoints& dummy) = delete;
 
 /// Broken assignment operator
 void operator=(const TwoDStraightLineFromTwoPoints&) = delete;
 
 /// Destructor
 ~TwoDStraightLineFromTwoPoints(){}
 
 /// Position Vector at Lagrangian coordinate zeta
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Position Vector
   r[0] = Left[0]+zeta[0]*(Right[0]-Left[0]);
   r[1] = Left[1]+zeta[0]*(Right[1]-Left[1]);
  }
 
 
 /// Parametrised position on object: r(zeta). Evaluated at
 /// previous timestep. t=0: current time; t>0: previous
 /// timestep.
 void position(const unsigned& t,
               const Vector<double>& zeta,
               Vector<double>& r) const
  {
   // Position Vector
   r[0] = Left[0]+zeta[0]*(Right[0]-Left[0]);
   r[1] = Left[1]+zeta[0]*(Right[1]-Left[1]);
  }
 
 
 /// Derivative of position Vector w.r.t. to coordinates:
 /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
 /// Evaluated at current time.
 virtual void dposition(const Vector<double>& zeta,
                        DenseMatrix<double>& drdzeta) const
  {
   // Tangent vector
   drdzeta(0, 0) = Right[0]-Left[0];
   drdzeta(0, 1) = Right[1]-Left[1];
  }
 
 
 /// 2nd derivative of position Vector w.r.t. to coordinates:
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
 /// ddrdzeta(alpha,beta,i). Evaluated at current time.
 virtual void d2position(const Vector<double>& zeta,
                         RankThreeTensor<double>& ddrdzeta) const
  {
   // Derivative of tangent vector
   ddrdzeta(0, 0, 0) = 0.0;
   ddrdzeta(0, 0, 1) = 0.0;
  }
 
 
 /// Posn Vector and its  1st & 2nd derivatives
 /// w.r.t. to coordinates:
 /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
 /// ddrdzeta(alpha,beta,i).
 /// Evaluated at current time.
 virtual void d2position(const Vector<double>& zeta,
                         Vector<double>& r,
                         DenseMatrix<double>& drdzeta,
                         RankThreeTensor<double>& ddrdzeta) const
  {
   // Position Vector
   r[0] = Left[0]+zeta[0]*(Right[0]-Left[0]);
   r[1] = Left[1]+zeta[0]*(Right[1]-Left[1]);
   
   // Tangent vector
   drdzeta(0, 0) = Right[0]-Left[0];
   drdzeta(0, 1) = Right[1]-Left[1];
   
   // Derivative of tangent vector
   ddrdzeta(0, 0, 0) = 0.0;
   ddrdzeta(0, 0, 1) = 0.0;
  }
 
 
 /// How many items of Data does the shape of the object depend on?
 unsigned ngeom_data() const
  {
   return 0;
  }
 
 /// Return pointer to the j-th Data item that the object's
 /// shape depends on
 Data* geom_data_pt(const unsigned& j)
  {
   return 0;
  }

private:
 
 /// Left point 
 Vector<double> Left;

 /// Right point
 Vector<double> Right; 
 
};




/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
// Polynomial approximation to ellipse
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


 //=========================================================================
 /// Polynomial approximation to ellipse
 //=========================================================================
class PolynomialApproxToEllipse : public GeomObject
{
public:
 
 /// Constructor: Pass ellipse and start and end coordinates.
 /// Boundary goes through the start and end points of the ellipse
 /// and uses a polynomial representation of order  m_poly in
 /// between. Used for test of the boundary interpolation.
 /// For m_poly<=3 (or <=5) the boundary must be represented exactly
 /// by curved elements with boundary polynomial order of 3 (or 5). 
 PolynomialApproxToEllipse(Ellipse* ellipse_pt,
                           const double& zeta_start,
                           const double& zeta_end,
                           const unsigned& m_poly)
  : GeomObject(1, 2), Ellipse_pt(ellipse_pt),
    Zeta_start(zeta_start), Zeta_end(zeta_end), M_poly_dev(m_poly-2)
  {

   Vector<double> zeta(1);

   // Get left and right values
   Left.resize(2);
   zeta[0]=zeta_start;
   Ellipse_pt->position(zeta,Left);

   Right.resize(2);
   zeta[0]=zeta_end;
   Ellipse_pt->position(zeta,Right);
   
   // Polynomial coefficients: Fm_minus_2[p][i]
   Fm_minus_2.resize(M_poly_dev);
   for (unsigned p=0;p<M_poly_dev;p++)
    {
     Fm_minus_2[p] = {Random::random_between_zero_and_one(),
                      Random::random_between_zero_and_one()};
    }

   // Choose polynomial so that zeros are located
   // in relevant part of boundary coordinate
   Zeros_in_interval=true;

   // Max of polynomial 
   double product=1.0;
   double eval_point=0.5;
   if ((M_poly_dev+1)%2==0)
    {
     eval_point=0.5-1.0/double(M_poly_dev+2);
    }
   for (unsigned ii=0;ii<M_poly_dev+2;ii++)
    {
     double fract_zero=double(ii)/double(M_poly_dev+1);
     product*=(eval_point-fract_zero);
    }
   Ampl_of_deviation=0.0; // hierher1.0e-2/product; //1.0e-2/product;
  }
 
 /// Broken copy constructor
 PolynomialApproxToEllipse(const PolynomialApproxToEllipse& dummy) = delete;
 
 /// Broken assignment operator
 void operator=(const PolynomialApproxToEllipse&) = delete;
 
 /// Destructor
 ~PolynomialApproxToEllipse(){}
 
 /// Position Vector at Lagrangian coordinate zeta
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   double fract=(zeta[0]-Zeta_start)/(Zeta_end-Zeta_start);
   for (unsigned i=0;i<2;i++)
    {
     r[i]=Left[i]+(Right[i]-Left[i])*fract;
     if (Zeros_in_interval)
      {
       double product=1.0;
       for (unsigned ii=0;ii<M_poly_dev+2;ii++)
        {
         double fract_zero=double(ii)/double(M_poly_dev+1);
         product*=(fract-fract_zero);
        }
       r[i]+=Ampl_of_deviation*product; 
      }
     else
      {
       for (unsigned p=0;p<M_poly_dev;p++)
        {
         r[i]+=fract*(1.0-fract)*Fm_minus_2[p][i]*pow(fract,p);
        }
      }
    }
  }
 
 
 /// Parametrised position on object: r(zeta). Evaluated at
 /// previous timestep. t=0: current time; t>0: previous
 /// timestep.
 void position(const unsigned& t,
               const Vector<double>& zeta,
               Vector<double>& r) const
  {
   position(zeta,r);
  }
 
 
 /// Derivative of position Vector w.r.t. to coordinates:
 /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
 /// Evaluated at current time.
 virtual void dposition(const Vector<double>& zeta,
                        DenseMatrix<double>& drdzeta) const
  {
   double fract=(zeta[0]-Zeta_start)/(Zeta_end-Zeta_start);
   for (unsigned i=0;i<2;i++)
    {
     drdzeta(0,i)=(Right[i]-Left[i])/(Zeta_end-Zeta_start);
     if (Zeros_in_interval)
      {
       double sum=0.0;
       for (unsigned jj=0;jj<M_poly_dev+2;jj++)
        {
         double product=Ampl_of_deviation/(Zeta_end-Zeta_start);
         for (unsigned ii=0;ii<M_poly_dev+2;ii++)
          {
           if (ii!=jj)
            {
             double fract_zero=double(ii)/double(M_poly_dev+1);
             product*=(fract-fract_zero);
            }
          }
         sum+=product;
        }
       drdzeta(0,i)+=sum; 
      }
     else
      {
       for (unsigned p=0;p<M_poly_dev;p++)
        {
         drdzeta(0,i)+=
          1.0/(Zeta_end-Zeta_start)*
          (      (1.0-fract)*Fm_minus_2[p][i]*pow(fract,p)+
                 fract*(   -1.0  )*Fm_minus_2[p][i]*pow(fract,p)
           );
         if (p>0)
          {
           drdzeta(0,i)+=
            1.0/(Zeta_end-Zeta_start)*
            fract*(1.0-fract)*Fm_minus_2[p][i]*p*pow(fract,p-1);
          }
        }
      }
    }
  }
 
 
 /// 2nd derivative of position Vector w.r.t. to coordinates:
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
 /// ddrdzeta(alpha,beta,i). Evaluated at current time.
 virtual void d2position(const Vector<double>& zeta,
                         RankThreeTensor<double>& ddrdzeta) const
  {
   double fract=(zeta[0]-Zeta_start)/(Zeta_end-Zeta_start);
   for (unsigned i=0;i<2;i++)
    {
     ddrdzeta(0,0,i)=0.0;
     if (Zeros_in_interval)
      {
       // Clever trick from recursive definition of polynomial
       // and its derivatives (thanks, ChatGPT)
       double P  = 1.0;
       double P1 = 0.0;
       double P2 = 0.0;
       for (unsigned ii=0;ii<M_poly_dev+2;ii++)
        {
         double fract_zero=double(ii)/double(M_poly_dev+1);
         double d = fract - fract_zero;         
         P2 = d * P2 + 2.0 * P1;
         P1 = d * P1 + P;
         P  = d * P;
        }
       ddrdzeta(0,0,i)+=Ampl_of_deviation/pow((Zeta_end-Zeta_start),2)*P2;
      }
     else
      {
       for (unsigned p=0;p<M_poly_dev;p++)
        {
         double sum=
          (     -1.0)*Fm_minus_2[p][i]*pow(fract,p)+
          (     -1.0)*Fm_minus_2[p][i]*pow(fract,p);
         if (p>0)
          {
           sum+=
            (1.0-fract)*Fm_minus_2[p][i]*p*pow(fract,p-1)+
            fract*(     -1.0)*Fm_minus_2[p][i]*p*pow(fract,p-1)+
            (1.0-fract)*Fm_minus_2[p][i]*p*pow(fract,p-1)+
            fract*(-1.0     )*Fm_minus_2[p][i]*p*pow(fract,p-1);
           if (p>1)
            {
             sum+=
              fract*(1.0-fract)*Fm_minus_2[p][i]*p*(p-1)*pow(fract,p-2);
            }
          }
         ddrdzeta(0,0,i)+=sum/pow((Zeta_end-Zeta_start),2);
        }
      }
    }
  }
 
 
 /// Posn Vector and its  1st & 2nd derivatives
 /// w.r.t. to coordinates:
 /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
 /// ddrdzeta(alpha,beta,i).
 /// Evaluated at current time.
 virtual void d2position(const Vector<double>& zeta,
                         Vector<double>& r,
                         DenseMatrix<double>& drdzeta,
                         RankThreeTensor<double>& ddrdzeta) const
  {
   oomph_info << "hierher broken" << std::endl;
   abort();
  }
 
 
 /// How many items of Data does the shape of the object depend on?
 unsigned ngeom_data() const
  {
   return 0;
  }
 
 /// Return pointer to the j-th Data item that the object's
 /// shape depends on
 Data* geom_data_pt(const unsigned& j)
  {
   return 0;
  }

private:

 /// Left point 
 Vector<double> Left;

 /// Right point
 Vector<double> Right;
 
 /// Pointer to approximated ellipse
 Ellipse* Ellipse_pt;

 /// Start coordinate
 double Zeta_start;
 
 /// End coordinate
 double Zeta_end;

 /// Polynomial order for deviation from straight line
 unsigned M_poly_dev;
 
 /// Polynomial coefficients: Fm_minus_2[p][i]
 Vector<Vector<double>> Fm_minus_2;

 /// Amplitude of deviation
 double Ampl_of_deviation;

 /// Choose polynomial so that zeros are located
 /// in relevant part of boundary coordinate
 bool Zeros_in_interval;
};



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////






//========================================================================
/// Dimensional parameters
//========================================================================
namespace DimensionalParameters
{

 /// Gravity
 double Gravity=9.81; // metres/sec^2

 /// Sheet thickness
 double Thickness=0.8e-3; // metres

 /// Sheet radius
 double Radius=0.2; // this works: 10.0e-2; // metres

 /// Density
 double Density=900.0; // kg/m^3;

 /// Young's modulus
 double Youngs_modulus=1.44e6; // Newton/metre^2

 /// Poisson's ratio
 double Poisson_ratio=0.5;

}





//========================================================================
/// Namespace for problem parameters
//========================================================================
namespace Parameters
{

 /// Number of plot points
 unsigned Nplot=5;

 /// Rotate coordinates on curvilinear boundaries?
 bool Rotate_coordinates_on_all_curvilinear_boundaries=true;

 /// Ellipse half x-axis
 double A = 1.0;
 
 /// Ellipse half y-axis
 double B = 1.0;

 /// Damping constant for damped solves (magnitude sort of irrelevant
 /// since the adaptive timestepping will kick in anyway).
 double Mu =1.0; 


 /// Nondimensional thickness of plate -- dependent parameter compute!
 double Thickness = 0.0;

 
 /// Membrane coupling coefficient (a dependent parameter)
 double Eta = 0.0; // hierher does it have the 1-nu^2 in it?)
                   // 12.0 * (1.0 - Nu * Nu) / (Thickness * Thickness);


 /// Max non-dimensional pressure on bending scale; dependent parameter compute
 double P_max=0.0;

 /// Pressure magnitude
  double P_mag = 0.0;

 /// pressure perturbation
 double P_cos=0.0;

 /// Wavenumber for pressure perturbation
 unsigned N_cos=6;
 
 /// Element area
 double Element_area = 0.5;

 
 /// Pressure depending on the position (x,y)
  void get_pressure(const Vector<double>& x, double& pressure)
  {
   double phi=atan2(x[1],x[0]);
   pressure = P_mag + P_cos*cos(phi*double(N_cos));
  }

  /// In plane forcing (shear stress) depending on the position (x,y)
  void get_in_plane_traction(const Vector<double>& x, Vector<double>& tau)
  {
   // Zero shear stress
   tau[0]=0.0;
   tau[1]=0.0;
  }



 /// Compute/updated dependent non-dimensional parameters
 void update_nondim_parameters()
 {
  // Non-dim thickness
  Thickness=DimensionalParameters::Thickness/DimensionalParameters::Radius;
  
  // FvK parameter
  Eta=12.0*(1.0-DimensionalParameters::Poisson_ratio*
            DimensionalParameters::Poisson_ratio)/(Thickness*Thickness);
  
  // Max. pressure (corresponding to full gravity)
  P_max=DimensionalParameters::Density*
   DimensionalParameters::Gravity*
   DimensionalParameters::Thickness*
   (1.0-DimensionalParameters::Poisson_ratio*
    DimensionalParameters::Poisson_ratio)/
   (sqrt(12.0)*DimensionalParameters::Youngs_modulus)*
   pow(Eta,1.5);
  
  oomph_info
   << "Updated non-dim parameters. \n"
   << "Thickness = " << Thickness << std::endl
   << "Eta       = " << Eta << std::endl
   << "P_max     = " << P_max << std::endl
   << std::endl;
 }


 ///////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////

 
 //===========================================================================
 /// Class to define zero C0 boundary conditions: f=0 for all zeta.
 /// Can be used for in-plane FvK displacements.
 //===========================================================================
 class ZeroC0BoundaryConditions : public virtual BoundaryConditionForC1PlateBending
 {
  
  /// Implement pure virtual function to specify value of the function
  /// (typically a displacement
  /// component) as a function of zeta, the 1D coordinate that parametrises the
  /// boundary
  virtual double f(const double& zeta)
   {return 0.0;}
  
 };
 
 

 ///////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////

 
 //===========================================================================
 /// Class to define zero C1 boundary conditions: f=df/dn=0 for all zeta
 /// Can be used for out-of-plane displacements in FvK and Koiter Steigman.
 //===========================================================================
 class ZeroC1BoundaryConditions : public virtual BoundaryConditionForC1PlateBending
 {
  
  /// Implement pure virtual function to specify value of the function (typically a
  /// displacement component) as a function of zeta, the 1D coordinate that parametrises
  /// the boundary
  virtual double f(const double& zeta)
   {return 0.0;}
  
  
  /// Override broken virtual function to specify the normal derivative of the function
  /// (typically a displacement component) w.r.t zeta, the 1D coordinate that
  /// parametrises the boundary. This is only needed for genuine C1 quantities
  /// (or for large-amplitude problems, e.g. Koiter Steigman, where all three
  /// displacement components need to be C1.
  virtual double dfdn(const double& zeta)
   {return 0.0;}

  
 };

 
} // end parameters namespace




///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////




//==start_of_problem_class============================================
/// Problem definition
//====================================================================
template<class ELEMENT>
class UnstructuredC1PlateProblem : public virtual Problem
{

public:

  /// Constructor
 UnstructuredC1PlateProblem(double const& element_area,
                            const unsigned& m_poly_actual_boundary,
                            const unsigned& boundary_order,
                            bool use_square_domain,
                            const double& phi);

  /// Destructor
  ~UnstructuredC1PlateProblem()
  {}

 /// Overloaded version of the problem's access function to
 /// the mesh. 
 TriangleMesh<ELEMENT>* mesh_pt()
  {
   return Bulk_mesh_pt;
  }
 
  /// Update after solve (empty)
  void actions_after_newton_solve() {}

  /// Update the problem specs before solve: empty
  void actions_before_newton_solve(){}
 
 /// Doc the solution
 void doc_solution(bool steady = true);


 /// Validate all basis functions for curved bell
 void validate_curved_bell_and_bubble_basis_functions(const std::string&
                                                      dir_name_for_output);
 
 /// Plot all basis functions for curved bell
 void plot_curved_bell_and_bubble_basis_functions(const std::string&
                                                  dir_name_for_output);
 

 /// Validate interpolated_x
 void validate_interpolated_x(const std::string&
                              dir_name_for_output,
                              const unsigned& m_poly_actual_boundary,
                              const unsigned& boundary_order);
 
 /// Doc/check boundary coordinates
 void doc_boundary_coords()
  {
   unsigned nb=Bulk_mesh_pt->nboundary();
   for (unsigned b=0;b<nb;b++)
    {
     std::string filename="boundary_coordinate"+to_string(b)+".dat";
     std::ofstream outfile;
     outfile.open(filename.c_str());
     oomph_info << "Checking boundary " << b << std::endl;
     const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);
     oomph_info << "Number of elements on boundary " << b << " : " << nb_element << std::endl;
     for(unsigned e=0;e<nb_element;e++)
      {
       // Get pointer to bulk element adjacent to b
       ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b,e));
       
       unsigned n_node=el_pt->nnode();
       oomph_info << "Element " << e << " has " << n_node << " nodes " << std::endl;
       for (unsigned n = 0; n < n_node; ++n)
        {
         // Get boundary node
         BoundaryNode<Node>* nod_pt =
          dynamic_cast<BoundaryNode<Node>*>(el_pt->node_pt(n));
        if (nod_pt==0)
         {
          oomph_info << "Node n = " << n << " at "
                     << el_pt->node_pt(n)->x(0) << " "
                     << el_pt->node_pt(n)->x(1) << " "
                     << " is not a boundary node" << std::endl;
         }
        else
         {
          std::set<unsigned>* boundaries_pt=0;
          nod_pt->get_boundaries_pt(boundaries_pt);
          if (boundaries_pt==0)
           {
            oomph_info << "Node n = " << n << " at "
                       << nod_pt->x(0) << " "
                       << nod_pt->x(1) << " "
                       << " is not on any boundary" << std::endl;
           }
          else
           {
            oomph_info << "Node n = " << n << " at "
                       << nod_pt->x(0) << " "
                       << nod_pt->x(1) << " "
                       << " is on boundaries:";
            for (auto b : *boundaries_pt)
             {
              oomph_info << b << " ";
             }
            oomph_info << std::endl;
           }
          
          // Check if it is on the boundary
          if (nod_pt->is_on_boundary(b))
           {            
            // We should only have one coordinate on this boundary
            unsigned nzeta=nod_pt->ncoordinates_on_boundary(b);
#ifdef PARANOID
            if (nzeta!=1)
             {
              oomph_info << "Why do we have more than one boundary coordinate?"
                         << std::endl;
              abort();
             }
#endif
            
            Vector<double> zeta(nzeta);
            nod_pt->get_coordinates_on_boundary(b,zeta);            
            outfile << nod_pt->x(0) << " "
                    << nod_pt->x(1) << " "
                    << zeta[0] <<  std::endl;
            
           }
         }
        }
      }
     outfile.close();
    }
  }


 
 
private:


 
 /// Doc boundary elements and faces
 // hierher move this into triangle mesh and describe output and prefix
 void doc_boundary_elements_and_faces(Mesh* mesh_pt, std::string name_prefix="")
  {
   
   // Loop over boundaries
   unsigned nbound=mesh_pt->nboundary();
   std::cout << "Number of boundaries: " << nbound << std::endl;
   for (unsigned b = 0; b < nbound; b++)
    {
     ofstream bulk_file;
     bulk_file.open(name_prefix+"bulk_elements_on_boundary"+to_string(b)+".dat");
     ofstream face_file;
     face_file.open(name_prefix+"face_elements_on_boundary"+to_string(b)+".dat");
     unsigned nel = mesh_pt->nboundary_element(b);
     oomph_info << "Boundary: " << b << " is adjacent to " << nel << " elements"
               << std::endl;
     
     // Loop over elements on given boundary
     for (unsigned e = 0; e < nel; e++)
      {
       FiniteElement* fe_pt = mesh_pt->boundary_element_pt(b,e);
       unsigned nnod_1d=fe_pt->nnode_1d();
       fe_pt->output(bulk_file,nnod_1d);
       unsigned face_index=mesh_pt->face_index_at_boundary(b,e);
       oomph_info << "Boundary element:" << fe_pt
                 << " Face index of boundary is "
                 << face_index
                 << std::endl;
       
       // Build face element (we have a TElement!)
       FiniteElement* face_element_pt;
       switch (nnod_1d)
        {
        case 4:
        {
         face_element_pt = new DummyFaceElement<TElement<2,4>>(fe_pt,face_index);
        }
        break;
        
        case 3:
        {
         face_element_pt = new DummyFaceElement<TElement<2,3>>(fe_pt,face_index);
        }
        break;
        
        case 2:
        {
         face_element_pt = new DummyFaceElement<TElement<2,2>>(fe_pt,face_index);
        }
        break;
        
        default:
         // hierher throw properly
         oomph_info << "Never get here: " << nnod_1d << std::endl;
         abort();
        }
       
       face_file << "ZONE" << std::endl;
       for (unsigned j=0;j<nnod_1d;j++)
        {
         face_file << face_element_pt->node_pt(j)->x(0) << " "
                   << face_element_pt->node_pt(j)->x(1) << " "
                   << std::endl;
        }
       delete face_element_pt;

      }

     bulk_file.close();
     face_file.close();
    }
  }

  /// Pointer to "bulk" mesh
  TriangleMesh<ELEMENT>* Bulk_mesh_pt;

 
  /// Enumeration to keep track of boundary ids
  enum
  {
    Outer_boundary0 = 0,
    Outer_boundary1 = 1,
    Outer_boundary2 = 2,
    Outer_boundary3 = 3
  };

  /// Target element area
  double Element_area;

  /// Pointer to constraint mesh
  Mesh* Constraint_mesh_pt;

  /// Doc info object for labeling output
  DocInfo Doc_info;

}; // end_of_problem_class



//======================================================================
/// Constructor definition
//======================================================================
template<class ELEMENT>
UnstructuredC1PlateProblem<ELEMENT>::UnstructuredC1PlateProblem
(const double& element_area,
 const unsigned& m_poly_actual_boundary,
 const unsigned& boundary_order,
 bool use_square_domain,
 const double& phi)
 : Element_area(element_area)
{
 
 
 // Set directory
 Doc_info.set_directory("RESLT");

 
 // Build the mesh
 //================
 TriangleMeshClosedCurve* outer_boundary_pt=0;
     
 //Outer boundary
 //-------------
 if (use_square_domain)
  {
   // Straight lines
   Vector<double> left(2);
   Vector<double> right(2);
   double zeta_start = 0.0;
   double zeta_end = 1.0;
   unsigned nsegment=4;

    
   // Storage for outer boundaries (for triangle)
   Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(4);

   // Right
   left[0] = 1.0;
   left[1] =-1.0;
   right[0]= 1.0;
   right[1]= 1.0;
   TwoDStraightLineFromTwoPoints* right_line_pt =
    new TwoDStraightLineFromTwoPoints(left,right,phi);
    outer_curvilinear_boundary_pt[0] =
    new TriangleMeshCurviLine(right_line_pt, zeta_start,
                              zeta_end, nsegment, Outer_boundary0);
  
   // Top
   left[0] = 1.0;
   left[1] = 1.0;
   right[0]=-1.0;
   right[1]= 1.0;
   TwoDStraightLineFromTwoPoints* top_line_pt =
    new TwoDStraightLineFromTwoPoints(left,right,phi);
    outer_curvilinear_boundary_pt[1] =
    new TriangleMeshCurviLine(top_line_pt, zeta_start,
                              zeta_end, nsegment, Outer_boundary1);
   // Left
   left[0] =-1.0;
   left[1] = 1.0;
   right[0]=-1.0;
   right[1]=-1.0;
   TwoDStraightLineFromTwoPoints* left_line_pt =
    new TwoDStraightLineFromTwoPoints(left,right,phi);
    outer_curvilinear_boundary_pt[2] =
    new TriangleMeshCurviLine(left_line_pt, zeta_start,
                              zeta_end, nsegment, Outer_boundary2);

   // Bottom
   left[0] =-1.0;
   left[1] =-1.0;
   right[0]= 1.0;
   right[1]=-1.0;
   TwoDStraightLineFromTwoPoints* bottom_line_pt =
    new TwoDStraightLineFromTwoPoints(left,right,phi);
    outer_curvilinear_boundary_pt[3] =
    new TriangleMeshCurviLine(bottom_line_pt, zeta_start,
                              zeta_end, nsegment, Outer_boundary3);
     
   // Combine
   outer_boundary_pt =
    new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);
   
  }
 else
  {
   double A = Parameters::A;
   double B = Parameters::B;
   Ellipse* outer_boundary_ellipse_pt = new Ellipse(A, B);
 
   // Storage for outer boundaries (for triangle)
   Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(4);

   //First bit
   double zeta_start = 0.0;
   double zeta_end = 0.5*MathematicalConstants::Pi;

   // Polynomial approximation for ellipse; matching at the end points
   PolynomialApproxToEllipse* poly_approx_pt=new PolynomialApproxToEllipse
    (outer_boundary_ellipse_pt,zeta_start,zeta_end,m_poly_actual_boundary);
 
   unsigned nsegment = (unsigned)(MathematicalConstants::Pi/sqrt(Element_area));
   outer_curvilinear_boundary_pt[0] = 
    new TriangleMeshCurviLine(poly_approx_pt, zeta_start,
                              zeta_end, nsegment, Outer_boundary0);
 
 
   //Second bit

   // Straight line
   double zeta_start_next_curved=MathematicalConstants::Pi;
   Vector<double> left(2);
   Vector<double> zeta(1);
   zeta[0]=zeta_end;
   outer_boundary_ellipse_pt->position(zeta,left);
   Vector<double> right(2);
   zeta[0]=zeta_start_next_curved;
   outer_boundary_ellipse_pt->position(zeta,right);
   TwoDStraightLineFromTwoPoints* straight_line_pt =
    new TwoDStraightLineFromTwoPoints(left,right);
 
   zeta_start = 0.5*MathematicalConstants::Pi;
   zeta_end = MathematicalConstants::Pi;
   outer_curvilinear_boundary_pt[1] =
    new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
                              zeta_end, nsegment, Outer_boundary1);
  
   
   //Third bit
   zeta_start = zeta_start_next_curved;
   zeta_end = 1.5*MathematicalConstants::Pi;
   outer_curvilinear_boundary_pt[2] = 
    new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
                              zeta_end, nsegment, Outer_boundary2);
 
 
   //Fourth bit
   zeta_start = 1.5*MathematicalConstants::Pi;
   zeta_end = 2.0*MathematicalConstants::Pi;
   outer_curvilinear_boundary_pt[3] =
    new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
                              zeta_end, nsegment, Outer_boundary3);
 
   // Combine
   outer_boundary_pt =
    new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);
   }

 
 
 //Create mesh parameters object
 TriangleMeshParameters mesh_parameters(outer_boundary_pt);

 // Element area
 mesh_parameters.element_area() = Element_area;
 
 // Build an assign bulk mesh
 Bulk_mesh_pt=new TriangleMesh<ELEMENT>(mesh_parameters);
 

  
  // Now upgrade to (potentially) curved C1 boundaries 
  //==================================================
  {

  
   // Let's have a look at the orig mesh
   Bulk_mesh_pt->output(Doc_info.directory()+
                        "/mesh_before_black_box_upgrade.dat");




   // // hierher test what element splitting does ad confirm that it updates
   // the boundary lookup scheme
   // {
   //  // Let's have a look at the new mesh
   //  Bulk_mesh_pt->output("mesh_before_splitting.dat");
   //  doc_boundary_coords();
   //  std::string name_prefix="test_before_splitting_";
   //  doc_boundary_elements_and_faces(Bulk_mesh_pt,name_prefix);

   //  // Split elements that have multiple edges on a boundary
   //  // Note: Sets up the boundary loopup scheme too.
   //  Bulk_mesh_pt->
   //   template split_elements_with_multiple_boundary_edges<ELEMENT>();
   //  // (Split_elements_output_stream);
    
    
    
   //  // Let's have a look at the new mesh
   //  Bulk_mesh_pt->output("mesh_after_splitting.dat");
   //  doc_boundary_coords();
   //  name_prefix="test_after_splitting_";
   //  doc_boundary_elements_and_faces(Bulk_mesh_pt,name_prefix);

   // }





   
   // Create the mesh for the Lagrange multiplier elements that enforce
   // continuity of our smooth solution across different parts of the
   // mesh boundary (only really needed when there are kinks)
   Constraint_mesh_pt = new Mesh();


   // Let's have a look what the black box helper function does:
   C1PlateHelper::Duplicated_node_output_stream.open
    (Doc_info.directory()+"duplicated_nodes.dat");
   C1PlateHelper::Upgraded_to_curved_edge_element_stream.open
    (Doc_info.directory()+"elements_upgraded_to_curved.dat");
   C1PlateHelper::Split_elements_output_stream.open
    (Doc_info.directory()+"split_elements.dat");
   C1PlateHelper::Rotated_node_output_stream.open
    (Doc_info.directory()+"rotated_nodes.dat");
   C1PlateHelper::Rotated_element_output_stream.open
    (Doc_info.directory()+"rotated_elements.dat");
    
   // Rotate coordinates on curvilinear boundaries?
   C1PlateHelper::upgrade_triangle_mesh_for_c1_plate_bending<ELEMENT>(
    Bulk_mesh_pt,
    Constraint_mesh_pt,
    Parameters::Rotate_coordinates_on_all_curvilinear_boundaries,
    boundary_order);

   // Done
   C1PlateHelper::Duplicated_node_output_stream.close();
   C1PlateHelper::Upgraded_to_curved_edge_element_stream.close();
   C1PlateHelper::Split_elements_output_stream.close();
   C1PlateHelper::Rotated_node_output_stream.close();
   C1PlateHelper::Rotated_element_output_stream.close();

   // // Let's have a look at the new mesh
   // Bulk_mesh_pt->output(Doc_info.directory()+"mesh_black_box_upgrade.dat");
   // doc_boundary_coords();
   // std::string name_prefix=Doc_info.directory()+"/test_";
   // doc_boundary_elements_and_faces(Bulk_mesh_pt,name_prefix);

  }

  // Build global mesh
  //==================
  
  //Add submeshes to problem
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Constraint_mesh_pt);

  // Combine submeshes into a single Mesh 
  build_global_mesh();


   
  // Update the corner constraints based on the applied
  // boundary conditions. NOTE: This must be called
  // after all boundary conditions have been applied.
  unsigned n_el = Constraint_mesh_pt->nelement();
  for(unsigned i_el = 0; i_el < n_el; i_el++)
   {
    // hierher rename in src
    dynamic_cast<DuplicateNodeConstraintElement*>
     (Constraint_mesh_pt->element_pt(i_el))
     ->pin_redundant_constraints();
   }
  
  // Assign equation numbers
  oomph_info << "Number of equations: "
             << assign_eqn_numbers()
             << "\n\n\n\n\n";

  oomph_info << "Mesh built with";
  if (!Parameters::Rotate_coordinates_on_all_curvilinear_boundaries)
   {
    oomph_info << "out (!)";
   }
  oomph_info << " rotating coordinates on curvilinear boundaries\n"
             << " and square is rotated by angle " << phi
             << std::endl << std::endl << std::endl;

  
} // end Constructor








//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredC1PlateProblem<ELEMENT>::doc_solution(bool steady)
{
 ofstream some_file,some_file2;
 char filename[100];
 

 oomph_info << "Docing soln" << Doc_info.number()  << ".dat" << std::endl;
 
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file ,Parameters::Nplot);
 some_file.close();

 if (steady)
  {
   sprintf(filename,"%s/steady_soln%i.dat",Doc_info.directory().c_str(),
           Doc_info.number());
   some_file.open(filename);
   Bulk_mesh_pt->output(some_file ,Parameters::Nplot);
   some_file.close();
  }

  // Increment the doc_info number
  Doc_info.number()++;

} // end of doc



//========================================================================
/// Plot all basis functions for curved bell
//========================================================================
template<class ELEMENT>
void UnstructuredC1PlateProblem<ELEMENT>::plot_curved_bell_and_bubble_basis_functions(
 const std::string& dir_name_for_output)
{

 ofstream some_file;
 char filename[100];
 
// Test & plot 'em
 bool plot_em=true;
 if (dir_name_for_output=="") plot_em=false;
 
 // Find a curved element on the outer boundary
 unsigned b=Outer_boundary0;
 {
  unsigned e=0;
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));

   if (plot_em)
    {
     sprintf(filename,"%s/test_curved_element.dat",
             dir_name_for_output.c_str());
     some_file.open(filename);
     
     // Output the lot
     unsigned nplot=30;
     el_pt->full_output(some_file,nplot);
     
     some_file.close();
    }

   // Get map to curvline boundaries of mesh
   std::map<unsigned, C1CurviLine*> c1_curviline_boundary_pt =
    Bulk_mesh_pt->c1_curviline_boundary_pt();
   C1CurviLine* curviline_pt=c1_curviline_boundary_pt[b];
   
    
   // Find the dimension of the element 
   const unsigned dim = el_pt->dim(); //2;
    
   // The number of first derivatives is the dimension of the element
   const unsigned n_deriv = dim;

   // The number of second derivatives is the triangle number of the dimension
   const unsigned n_2deriv = dim * (dim + 1) / 2;
    
   // Find out how many nodes there are for w
   const unsigned n_w_node = el_pt->nw_node(); // 3; 

   // Get the vector of nodes used for each field
   const Vector<unsigned> w_nodes = el_pt->get_w_node_indices(); // {0,1,2}; 

   // Find out how many basis types there are at each node
   const unsigned n_w_nodal_type = el_pt->nw_type_at_each_node(); // 6; 

   // Find out how many basis types there are internally
   unsigned n_w_internal_type =  el_pt->nw_type_internal(); // 3 or 10

   oomph_info << "n_w_internal_type =  " << n_w_internal_type << std::endl;
    
   // Out-of-plane local basis & test functions
   // ------------------------------------------
   // Nodal basis & test functions
   Shape psi_n_w(n_w_node, n_w_nodal_type);
   Shape test_n_w(n_w_node, n_w_nodal_type);
   DShape dpsi_n_wdxi(n_w_node, n_w_nodal_type, n_deriv);
   DShape dtest_n_wdxi(n_w_node, n_w_nodal_type, n_deriv);
   DShape d2psi_n_wdxi2(n_w_node, n_w_nodal_type, n_2deriv);
   DShape d2test_n_wdxi2(n_w_node, n_w_nodal_type, n_2deriv);
   
   // Internal basis & test functions
   Shape psi_i_w(n_w_internal_type);
   Shape test_i_w(n_w_internal_type);
   DShape dpsi_i_wdxi(n_w_internal_type, n_deriv);
   DShape dtest_i_wdxi(n_w_internal_type, n_deriv);
   DShape d2psi_i_wdxi2(n_w_internal_type, n_2deriv);
   DShape d2test_i_wdxi2(n_w_internal_type, n_2deriv);
    
    
   // Plot all basis functions
   if (plot_em)
    {
     // Tecplot header info from some generic triangle element
     TElement<2,2>* aux_el_pt= new TElement<2,2>;
     unsigned nplot=100;

     /// Plot the whole thing or just the edge
     for (unsigned do_curved_edge=0;do_curved_edge<2;do_curved_edge++)
      {
       
       
       // Nodal (curved Bell) basis functions
       Vector<ofstream*> nodal_file_pt;
       unsigned count=0;
       for (unsigned j=0;j<n_w_node;j++)
        { 
         for (unsigned k=0;k<n_w_nodal_type;k++)
          {
           if (do_curved_edge==1)
            {
             sprintf(filename,"%s/test_curved_bell_curved_edge_nodal_basis%i.dat",
                     dir_name_for_output.c_str(),count);
            }
           else
            {
             sprintf(filename,"%s/test_curved_bell_nodal_basis%i.dat",
                     dir_name_for_output.c_str(),count);
            }
           nodal_file_pt.push_back(new ofstream);
           nodal_file_pt[count]->open(filename);
           
           // Tecplot header info
           if (do_curved_edge==0)
            {
             *(nodal_file_pt[count]) << aux_el_pt->tecplot_zone_string(nplot);
            }
           count++;
          }
        }
       
       
       // Internal (bubble) basis
       Vector<ofstream*> internal_file_pt;
       count=0;
       for (unsigned k_type = 0; k_type < n_w_internal_type; k_type++)
        {
         if (do_curved_edge==1)
          {
           sprintf(filename,"%s/test_curved_bell_curved_edge_bubble_basis%i.dat",
                   dir_name_for_output.c_str(),count);
          }
         else
          {
           sprintf(filename,"%s/test_curved_bell_bubble_basis%i.dat",
                   dir_name_for_output.c_str(),count);
          }
         internal_file_pt.push_back(new ofstream);
         internal_file_pt[count]->open(filename);
         
         // Tecplot header info
         if (do_curved_edge==0)
          {
           *(internal_file_pt[count]) << aux_el_pt->tecplot_zone_string(nplot);
          }
         count++;
        }
    
     

       // Only used for curved edge check
       Vector<double> r_from_boundary(2,0.0);
       Vector<double> drdzeta(2,0.0);
       Vector<double> zeta(1);
                             
       // Loop over plot points
       Vector<double> s_plot(2);
       unsigned num_plot_points = aux_el_pt->nplot_points(nplot);
       for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
        {
         // Get local coordinates of plot point
         if (do_curved_edge==0)
          {
           aux_el_pt->get_s_plot(iplot, nplot, s_plot);
          }
         else
          {
           s_plot[0]=double(iplot)/double(num_plot_points-1);
           s_plot[1]=1.0-s_plot[0];
          

           /// Position r as fct of zeta from curvilinear boundary representation
           zeta[0]=el_pt->bernadou_element_basis_pt()->get_s_ubar()+
            s_plot[1]*(el_pt->bernadou_element_basis_pt()->get_s_obar()-
                       el_pt->bernadou_element_basis_pt()->get_s_ubar());
           curviline_pt->position(zeta,r_from_boundary);
   
           /// Derivative of position Vector w.r.t. to zeta:
           curviline_pt->dposition(zeta, drdzeta);
          }

         // Get plot point
         Vector<double> interp_x(dim, 0.0);
         el_pt->interpolated_x(s_plot, interp_x);

       
         // Call the derivatives of the shape and test functions for the out of
         // plane unknown
         //double J =
          el_pt->d2basis_and_d2test_w_eulerian_foeppl_von_karman(s_plot,
                                                                 psi_n_w,
                                                                 psi_i_w,
                                                                 dpsi_n_wdxi,
                                                                 dpsi_i_wdxi,
                                                                 d2psi_n_wdxi2,
                                                                 d2psi_i_wdxi2,
                                                                 test_n_w,
                                                                 test_i_w,
                                                                 dtest_n_wdxi,
                                                                 dtest_i_wdxi,
                                                                 d2test_n_wdxi2,
                                                                 d2test_i_wdxi2);


         count=0;
         for (unsigned j=0;j<n_w_node;j++)
          { 
           for (unsigned k=0;k<n_w_nodal_type;k++)
            {
             if (do_curved_edge==0)
              {
               *(nodal_file_pt[count]) << interp_x[0] << " "
                                       << interp_x[1] << " "
                                       << psi_n_w(j,k) << " "
                                       << dpsi_n_wdxi(j,k,0) << " "
                                       << dpsi_n_wdxi(j,k,1) << " "
                                       << d2psi_n_wdxi2(j,k,0) << " "
                                       << d2psi_n_wdxi2(j,k,1) << " "
                                       << d2psi_n_wdxi2(j,k,2) << " "
                                       << s_plot[0] << " "
                                       << s_plot[1] << " " 
                                       << std::endl;
              }
             else
              {

               // * grad psi \cdot n should be zero for most basis fcts
               //   and vary cubically for the others (as in the Hermite basis, already added); so there
               //   are four nonzero functions.
               //
               //   Normal derivatives should be non-zero only for basis functions 1, 4, 7, 10
               //   ERROR: 0, 2, 5, 6, 8, 11 should be zero but aren't
               //          
               //
               // * psi should be zero for most basis fcts; the others should be a quintic basis
               //   (i.e. interpolants for value, deriv, and second deriv, i.e. six values =
               //    fifth order).
               //
               //   psi should be nonzero only for basis functions 0, 2, 5, 6, 8, 11 **CORRECT!**
               // NOTE WE'RE DOING ROTATED BASES! SO WE'RE INTERPOLATING VALUE v, dv/dn, dv/dt,
               // d^2v/dn^2, d^2v/dndt, d^2/dt^2
               //
               // NOTE: MAY NOT BE THE EXACT HERMITE POLYNOMIALS BECAUSE OF DIFFERENT PARAMETRISATION
               // OF THE DEFORMED ELEMENT (WORK IN PROGRESS, REVISIT)
               //
               double norm=sqrt(drdzeta[0]*drdzeta[0]+
                                drdzeta[1]*drdzeta[1]);
               double s=s_plot[1];
               *(nodal_file_pt[count]) << interp_x[0] << " " // 1
                                       << interp_x[1] << " " // 2
                                       << r_from_boundary[0] << " " // 3 
                                       << r_from_boundary[1] << " " // 4
                                       <<  drdzeta[1]/norm << " " // 5 
                                       << -drdzeta[0]/norm << " " // 6
                                       << s << " "  // 7
                                       << psi_n_w(j,k) << " " // 8
                                       << dpsi_n_wdxi(j,k,0) << " " // 9 (d/dn because we're rotated!)
                                       // << (  dpsi_n_wdxi(j,k,0)*drdzeta[1]
                                       //      -dpsi_n_wdxi(j,k,1)*drdzeta[0])/norm << " " // 9 (dpsi/dn)
                                       << 1.0 - 3.0*s*s + 2.0*s*s*s << " " // 10
                                       << s - 2.0*s*s + s*s*s << " " // 11
                                       << 3.0*s*s - 2.0*s*s*s << " " // 12
                                       << -s*s + s*s*s << " " // 13
                                       << std::endl;
              }
             
             count++;
            }
          }
       
         count=0;
         for (unsigned k_type = 0; k_type < n_w_internal_type; k_type++)
          {
           if (do_curved_edge==0)
            {
             *(internal_file_pt[count]) << interp_x[0] << " "
                                        << interp_x[1] << " "
                                        << psi_i_w(k_type) << " "
                                        << dpsi_i_wdxi(k_type,0) << " "
                                        << dpsi_i_wdxi(k_type,1) << " "
                                        << d2psi_i_wdxi2(k_type,0) << " "
                                        << d2psi_i_wdxi2(k_type,1) << " "
                                        << d2psi_i_wdxi2(k_type,2) << " "
                                        << s_plot[0] << " "
                                        << s_plot[1] << " " 
                                        << std::endl;
            }
           else
            {
             double norm=sqrt(drdzeta[0]*drdzeta[0]+
                              drdzeta[1]*drdzeta[1]);
             double s=s_plot[1];             
              *(internal_file_pt[count])  << interp_x[0] << " " // 1
                                       << interp_x[1] << " " // 2
                                       << r_from_boundary[0] << " " // 3 
                                       << r_from_boundary[1] << " " // 4
                                       <<  drdzeta[1]/norm << " " // 5 
                                       << -drdzeta[0]/norm << " " // 6
                                       << s << " "  // 7
                                       << psi_i_w(k_type) << " " // 8
                                          << dpsi_i_wdxi(k_type,0) << " " // 9 (d/dn because we're rotated!)
                                       // << (  dpsi_i_wdxi(k_type,0)*drdzeta[1]
                                       //      -dpsi_i_wdxi(k_type,1)*drdzeta[0])/norm << " " // 9 (dpsi/dn)
                                       << 1.0 - 3.0*s*s + 2.0*s*s*s << " " // 10
                                       << s - 2.0*s*s + s*s*s << " " // 11
                                       << 3.0*s*s - 2.0*s*s*s << " " // 12
                                       << -s*s + s*s*s << " " // 13
                                       << std::endl;
            }
  
           count++;
          }
        }

     
       // Write tecplot footer (e.g. FE connectivity lists) & close
       count=0;
       for (unsigned j=0;j<n_w_node;j++)
        { 
         for (unsigned k=0;k<n_w_nodal_type;k++)
          {
           if (do_curved_edge==0)
            {
             aux_el_pt->write_tecplot_zone_footer(*(nodal_file_pt[count]), nplot);
            }
           nodal_file_pt[count]->close();
           delete nodal_file_pt[count];
           count++;
          }
        }
       count=0;
       for (unsigned k_type = 0; k_type < n_w_internal_type; k_type++)
        {       
         if (do_curved_edge==0)
          {
           aux_el_pt->write_tecplot_zone_footer(*(internal_file_pt[count]), nplot);
          }
         internal_file_pt[count]->close();
         delete internal_file_pt[count];
         count++;
        }
     
      }
     
     delete aux_el_pt;
     aux_el_pt=0;
     
     oomph_info << "\n\nPlot of curved bell basis functions done! Now do: " << std::endl;
     oomph_info << "cd RESLT" << std::endl;
     oomph_info << "gnuplot -c ../validate_dpsidn_bubble.gp" << std::endl;
     oomph_info << "gnuplot -c ../validate_dpsidn_nodal.gp" << std::endl;
     oomph_info << "gnuplot -c ../validate_psi_nodal.gp" << std::endl;
     oomph_info << "gnuplot -c ../validate_psi_bubble.gp" << std::endl;
     oomph_info << "display validate*png" << std::endl;
     oomph_info << std::endl;
    }
  }
 }
}
   

//========================================================================
/// Validate all basis functions for curved bell
//========================================================================
template<class ELEMENT>
void UnstructuredC1PlateProblem<ELEMENT>::validate_curved_bell_and_bubble_basis_functions(
 const std::string& dir_name_for_output)
{
 
 ofstream some_file;
 char filename[100];
 
// Test & plot 'em
 bool plot_em=true;
 if (dir_name_for_output=="") plot_em=false;
 
 // Find a curved element on the outer boundary
 unsigned b=Outer_boundary0;
 {
  unsigned e=0;
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));
   
   if (plot_em)
    {
     sprintf(filename,"%s/test_curved_element.dat",
             dir_name_for_output.c_str());
     some_file.open(filename);
     
     // Output the lot
     unsigned nplot=30;
     el_pt->full_output(some_file,nplot);
     
     some_file.close();
    }
   
   
   // Find the dimension of the element 
   const unsigned dim = el_pt->dim(); //2;
   
   // The number of first derivatives is the dimension of the element
   const unsigned n_deriv = dim;
   
   // The number of second derivatives is the triangle number of the dimension
   const unsigned n_2deriv = dim * (dim + 1) / 2;
   
   // Find out how many nodes there are for w
   const unsigned n_w_node = el_pt->nw_node(); // 3; 

   // Get the vector of nodes used for each field
   const Vector<unsigned> w_nodes = el_pt->get_w_node_indices(); // {0,1,2}; 

   // Find out how many basis types there are at each node
   const unsigned n_w_nodal_type = el_pt->nw_type_at_each_node(); // 6; 

   // Find out how many basis types there are internally
   unsigned n_w_internal_type =  el_pt->nw_type_internal(); // 3 or 10

   oomph_info << "n_w_internal_type =  " << n_w_internal_type << std::endl;
    
   // Out-of-plane local basis & test functions (bubble and nodal)
   // ------------------------------------------------------------
   // Nodal basis & test functions
   Shape psi_n_w(n_w_node, n_w_nodal_type);
   Shape test_n_w(n_w_node, n_w_nodal_type);
   DShape dpsi_n_wdxi(n_w_node, n_w_nodal_type, n_deriv);
   DShape dtest_n_wdxi(n_w_node, n_w_nodal_type, n_deriv);
   DShape d2psi_n_wdxi2(n_w_node, n_w_nodal_type, n_2deriv);
   DShape d2test_n_wdxi2(n_w_node, n_w_nodal_type, n_2deriv);
   
   // Internal basis & test functions
   Shape psi_i_w(n_w_internal_type);
   Shape test_i_w(n_w_internal_type);
   DShape dpsi_i_wdxi(n_w_internal_type, n_deriv);
   DShape dtest_i_wdxi(n_w_internal_type, n_deriv);
   DShape d2psi_i_wdxi2(n_w_internal_type, n_2deriv);
   DShape d2test_i_wdxi2(n_w_internal_type, n_2deriv);
   
   // Number of interpolation-type tests we have to do:
   unsigned n_interpolation_test=n_w_node*n_w_nodal_type+n_w_internal_type;
   
   oomph_info << "Number of interpolation-type tests to do: "
              << n_interpolation_test << std::endl;


   // Test matrix
   DenseDoubleMatrix test_matrix(n_interpolation_test);
   
   
   // "Points" where interpolation property ought to be satisfied:
   std::map<std::string, // type of interpolation (a,b,...)
            Vector<      // instances of this (e.g. the three vertex points a
             // end up in a Vector of length 3). 
             std::pair<  // for each instance store a pair:
              Vector<double>,     // first  part of pair stores the coordinates;
              Vector<std::string> // second part of pair stores the types of 
              // dofs (value, deriv, ...); in general there
              // are multiple ones (e.g. at the vertices we
              // have 6) so we store them in a vector
              >>> test_point;
   
   // At this stage we have only three vertex nodes ("a") and the internal data ("e")
   test_point["a"].resize(3);
   test_point["a"][0]={{1.0,0.0},{"w","dwdx","dwdy","d2wdx2","d2wdxdy","d2wdy2"}};
   test_point["a"][1]={{0.0,1.0},{"w","dwdx","dwdy","d2wdx2","d2wdxdy","d2wdy2"}};
   test_point["a"][2]={{0.0,0.0},{"w","dwdx","dwdy","d2wdx2","d2wdxdy","d2wdy2"}};
   switch (n_w_internal_type)
    {
    case 3:
     test_point["e"].resize(3);
     test_point["e"][2]={{0.5 ,0.25},{"w"}};
     test_point["e"][0]={{0.25,0.5 },{"w"}};
     test_point["e"][1]={{0.25,0.25},{"w"}};
     
     break;
     
    case 10:
     test_point["e"].resize(10);
     test_point["e"][6]={{1.0/6.0,4.0/6.0},{"w"}};
     test_point["e"][7]={{1.0/6.0,3.0/6.0},{"w"}};
     test_point["e"][8]={{1.0/6.0,2.0/6.0},{"w"}};
     test_point["e"][0]={{1.0/6.0,1.0/6.0},{"w"}};
     
     test_point["e"][1]={{2.0/6.0,1.0/6.0},{"w"}};
     test_point["e"][2]={{3.0/6.0,1.0/6.0},{"w"}};
     test_point["e"][3]={{4.0/6.0,1.0/6.0},{"w"}};
     
     test_point["e"][4]={{3.0/6.0,2.0/6.0},{"w"}};
     test_point["e"][5]={{2.0/6.0,3.0/6.0},{"w"}};
     
     test_point["e"][9]={{2.0/6.0,2.0/6.0},{"w"}};
     
     break;
     
    default:
     throw OomphLibError("Never get here!",
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   
   // Linearly enumerated basis functions themselves
   Shape psi(n_interpolation_test);
   
   // Linearly enumerated first derivs (x,y) or (n,t)
   DShape dpsi(n_interpolation_test,2);
   
   // Linearly enumerated  2nd derivs (xx, xy, yy) or (or nn, nt,tt) 
   DShape d2psi(n_interpolation_test,3);
   
   // Test & plot 'em
   bool plot_em=true;
   if (dir_name_for_output=="") plot_em=false;
   if (plot_em)
    {
     sprintf(filename,"%s/curved_bell_test_points.dat",
             dir_name_for_output.c_str());
     some_file.open(filename);
    }
   
   // output intermediate results to screen
   bool output_to_screen=false;
   
   // Loop over the dofs
   unsigned interpolation_condition_count=0;
   
   // (class of dof: a,[b,d],e)
   for (auto dof_class : test_point)
    {
     
     if (output_to_screen)
      {
       oomph_info << std::fixed << std::setprecision(1);
      }




     // TEST LOGIC:
     //============
     // - Rows of the test matrix loop over all interpolation conditions which specify
     //   the location at which a certain condition (on the value of a basis function or
     //   one of its derivatives being one for exactly one basis function and zero
     //   for all the others) ought to be satisfied.
     // - Columns of the test matrix, then loop over all basis functions and evaluated
     //   the relevant value (value itself or the specified derivative at that point)
     //   for all basis functions in the element.
     // If all goes well, the matrix therefore ought to be a unit matrix.
     
     
     // Loop over location of all dof locations of this class (a1,a2,a3,...)
     // dof_class.second is a vector containing the pairs of location and
     // quantities to be interpolated/checked
     unsigned count=0;
     for (auto dof_location_and_tests : dof_class.second)
      {
       // Tell us what you're doing
       if (output_to_screen)
        {
         oomph_info << dof_class.first << count << " : " << std::endl;
        }
       
       // List coordinates of test point (dof_location_and_tests.first is the vector
       // of coordinates)
       if (plot_em)
        {
         for (unsigned i=0;i<2;i++)
          {
           some_file << (dof_location_and_tests.first)[i] << " ";
          }
        }
       
       
       // Call the derivatives of the shape and test functions 
       //double J =
       el_pt->d2basis_and_d2test_w_eulerian_foeppl_von_karman(dof_location_and_tests.first,
                                                              psi_n_w,
                                                              psi_i_w,
                                                              dpsi_n_wdxi,
                                                              dpsi_i_wdxi,
                                                              d2psi_n_wdxi2,
                                                              d2psi_i_wdxi2,
                                                              test_n_w,
                                                              test_i_w,
                                                              dtest_n_wdxi,
                                                              dtest_i_wdxi,
                                                              d2test_n_wdxi2,
                                                              d2test_i_wdxi2);
       
       
       // Move across into 1D enumeration:
       unsigned counter=0;
       
       // Nodal first...
       for (unsigned j=0;j<n_w_node;j++)
        {
         for (unsigned k=0;k<n_w_nodal_type;k++)
          {
           psi(counter)=psi_n_w(j,k);
           dpsi(counter,0)=dpsi_n_wdxi(j,k,0);
           dpsi(counter,1)=dpsi_n_wdxi(j,k,1);
           d2psi(counter,0)=d2psi_n_wdxi2(j,k,0);
           d2psi(counter,1)=d2psi_n_wdxi2(j,k,1);
           d2psi(counter,2)=d2psi_n_wdxi2(j,k,2);
           counter++;
          }
        }
       // then the internal (bubble) ones:
       for (unsigned j=0;j<n_w_internal_type;j++)
        {
         psi(counter)=psi_i_w(j);
         dpsi(counter,0)=dpsi_i_wdxi(j,0);
         dpsi(counter,1)=dpsi_i_wdxi(j,1);
         d2psi(counter,0)=d2psi_i_wdxi2(j,0);
         d2psi(counter,1)=d2psi_i_wdxi2(j,1);
         d2psi(counter,2)=d2psi_i_wdxi2(j,2);
         counter++;
        }

       
       // Loop over test types: dof_location_and_tests.second
       // is the vector whose strings tell us what quantity
       // we're supposed to interpolate/test
       for (auto test_type : dof_location_and_tests.second)
        {       
         if (output_to_screen)
          {         
           oomph_info << test_type << " " << std::endl;
          }
         
         if (test_type=="w")
          {
           for (unsigned i=0;i<n_interpolation_test;i++)
            {
             test_matrix(interpolation_condition_count,i)=psi[i];
             if (output_to_screen) oomph_info << psi[i] << " ";
            }
          }
         else if (test_type=="dwdx")
          {
           for (unsigned i=0;i<n_interpolation_test;i++)
            {
             test_matrix(interpolation_condition_count,i)=dpsi(i,0);
             if (output_to_screen) oomph_info << dpsi(i,0) << " ";
            }
          }
         else if (test_type=="dwdy")
          {
           for (unsigned i=0;i<n_interpolation_test;i++)
            {
             test_matrix(interpolation_condition_count,i)=dpsi(i,1);
             if (output_to_screen) oomph_info << dpsi(i,1) << " ";
            }
          }
         else if (test_type=="d2wdx2")
          {
           for (unsigned i=0;i<n_interpolation_test;i++)
            {
             test_matrix(interpolation_condition_count,i)=d2psi(i,0);
             if (output_to_screen) oomph_info << d2psi(i,0) << " ";
            }
          }
         else if (test_type=="d2wdxdy")
          {
           for (unsigned i=0;i<n_interpolation_test;i++)
            {
             test_matrix(interpolation_condition_count,i)=d2psi(i,1);
             if (output_to_screen) oomph_info << d2psi(i,1) << " ";
            }
          }
         else if (test_type=="d2wdy2")
          {
           for (unsigned i=0;i<n_interpolation_test;i++)
            {
             test_matrix(interpolation_condition_count,i)=d2psi(i,2);
             if (output_to_screen) oomph_info << d2psi(i,2) << " ";
            }
          }
         else
          {
           throw OomphLibError("Never get here!",
                               OOMPH_CURRENT_FUNCTION,
                               OOMPH_EXCEPTION_LOCATION);
          }
         
         interpolation_condition_count++;
         if (output_to_screen) oomph_info << std::endl;
         
        }
       if (output_to_screen) oomph_info << std::endl;
       if (plot_em)
        {
         some_file << std::endl;
        }
       count++;
       
      }
    }
   if (plot_em)
    {
     some_file.close();
    }
   
   // Test: exactly one unit entry per row
   std::stringstream row_unit_ness_stream;
   std::stringstream col_unit_ness_stream;
   double tol=1.0e-10;
   bool test_passed=true;
   for (unsigned i=0;i<n_interpolation_test;i++)
    {
     unsigned count_one_in_row=0;
     unsigned count_non_one_in_row=0;
     unsigned count_one_in_col=0;
     unsigned count_non_one_in_col=0;
     for (unsigned j=0;j<n_interpolation_test;j++)
      {
       if (std::abs(test_matrix(i,j)    )<tol) count_non_one_in_row++;
       if (std::abs(test_matrix(i,j)-1.0)<tol)
        {
         count_one_in_row++;
         row_unit_ness_stream << "Unit entry in row " << i << " is ";
         if (i==j)
          {
           row_unit_ness_stream << " on diagonal" << std::endl;
          }
         else
          {
           row_unit_ness_stream << " off diagonal, namely in column "
                      << j << std::endl;
          }
          
        }
       if (std::abs(test_matrix(j,i)    )<tol) count_non_one_in_col++;
       if (std::abs(test_matrix(j,i)-1.0)<tol)
        {
         count_one_in_col++;         
         col_unit_ness_stream << "Unit entry in column " << j << " is ";
         if (i==j)
          {
           col_unit_ness_stream << " on diagonal" << std::endl;
          }
         else
          {
           col_unit_ness_stream << " off diagonal, namely in row "
                                << i << std::endl;
          }
        }

      }

     std::stringstream diagnostic;
     diagnostic << "Row "
                 << i << " has "
                 << count_one_in_row << " ones (should be 1) and "
                 << count_non_one_in_row << " non-ones (should be "
                 << n_interpolation_test-1 << ") "
                 << std::endl;
     
     if ((count_one_in_row!=1)||
         (count_one_in_col!=1)||
         (count_non_one_in_row!=(n_interpolation_test-1))||
         (count_non_one_in_col!=(n_interpolation_test-1)))
      {
       test_passed=false;
       oomph_info << "failed: " << diagnostic.str();
      }
     else
      {
       oomph_info << "passed: " << diagnostic.str();
      }
    }


   oomph_info << row_unit_ness_stream.str() << std::endl << std::endl;
   oomph_info << col_unit_ness_stream.str() << std::endl << std::endl;
   
   if (test_passed)
    {
     oomph_info << "Test of curved bell basis functions passed!"
                << std::endl;
    }
   else
    {
     oomph_info << "Test of curved bell basis functions failed!"
                << std::endl;
    }
   
   
// // Plot all basis functions
//    if (plot_em)
//     {
//      // Tecplot header info from some generic triangle element
//      TElement<2,2>* aux_el_pt= new TElement<2,2>;
//      unsigned nplot=100;
//      for (unsigned i=0;i<n_interpolation_test;i++)
//       { 
//        sprintf(filename,"%s/test_basic_basis%i.dat",
//                dir_name_for_output.c_str(),i);
//        some_file.open(filename);
       
//        // Tecplot header info
//        some_file << aux_el_pt->tecplot_zone_string(nplot);
       
//        // Loop over plot points
//        Vector<double> s_plot(2);
//        unsigned num_plot_points = aux_el_pt->nplot_points(nplot);
//        for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
//         {
//          // Get local coordinates of plot point
//          aux_el_pt->get_s_plot(iplot, nplot, s_plot);
         
//          Shape psi(n_interpolation_test);
//          b_pt->full_basic_polynomials(s_plot,psi);
//          DShape dpsi(n_interpolation_test,2); // first derivs
//          b_pt->dfull_basic_polynomials(s_plot,dpsi);
//          DShape d2psi(n_interpolation_test,3); // 2nd derivs xx, xy, yy
//          b_pt->d2full_basic_polynomials(s_plot,d2psi);
         
//          some_file << s_plot[0] << " "
//                    << s_plot[1] << " ";
//          some_file << psi[i] << " ";
//          some_file << dpsi(i,0) << " "
//                    << dpsi(i,1) << " ";
//          some_file << d2psi(i,0) << " "
//                    << d2psi(i,1) << " "
//                    << d2psi(i,2) << " ";
//          some_file << std::endl;
//         }
       
//        // Write tecplot footer (e.g. FE connectivity lists)
//        aux_el_pt->write_tecplot_zone_footer(some_file, nplot);
//        some_file.close();
//       }
     
//      delete aux_el_pt;
//      aux_el_pt=0;
     
//      // oomph_info << "\n\nPlot of basis functions done! Now do: " << std::endl;
//      // oomph_info << "oomph-convert -z test_basic_basis*dat" << std::endl;
//      // oomph_info << "makePvd test_basic_basis test_basic_basis.pvd" << std::endl;
//      // oomph_info << "oomph-convert -p2 test_points.dat " << std::endl;
//      // oomph_info << "paraview --state test_basic_basis.pvsm " << std::endl;
//      // oomph_info << std::endl;
//     }
  }
 }
}



//========================================================================
/// Validate interpolated_x. Should represent third or fifth order polynomials
/// exactly.
//========================================================================
template<class ELEMENT>
void UnstructuredC1PlateProblem<ELEMENT>::validate_interpolated_x(
 const std::string& dir_name_for_output,
 const unsigned& m_poly_actual_boundary,
 const unsigned& boundary_order)
{

 ofstream some_file;
 char filename[100];
 
// Test & plot 'em
 bool plot_em=true;
 if (dir_name_for_output=="") plot_em=false;
 if (plot_em)
  {
   sprintf(filename,"%s/test_poly_%i_bound_%i_interpolated_x.dat",
           dir_name_for_output.c_str(),
           m_poly_actual_boundary,
           boundary_order);
   some_file.open(filename);
  

   // Output the entire mesh
   sprintf(filename,"%s/interpolated_x_mesh_plot_poly_%i_bound_%i_interpolated_x.dat",
           dir_name_for_output.c_str(),
           m_poly_actual_boundary,
           boundary_order);
   std::string fname=filename;
   Bulk_mesh_pt->output(fname);
  }

 
 // Initialise
 double max_error=0.0;
 
 // Loop over outer boundary
 unsigned b=Outer_boundary0;
 {
  // loop over all of them
  unsigned nelem=Bulk_mesh_pt->nboundary_element(b);
  for (unsigned e=0;e<nelem;e++)
   {
    // Get pointer to bulk element adjacent to b
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
     Bulk_mesh_pt->boundary_element_pt(b,e));
    
    // Get map to curvline boundaries of mesh
    std::map<unsigned, C1CurviLine*> c1_curviline_boundary_pt =
     Bulk_mesh_pt->c1_curviline_boundary_pt();
    C1CurviLine* curviline_pt=c1_curviline_boundary_pt[b];
    
    // Only used for curved edge check
    Vector<double> r_from_boundary(2,0.0);
    Vector<double> zeta(1);
    
    // Which edge is the curved one?
    unsigned curved_edge=
     el_pt->bernadou_element_basis_pt()->curved_edge();
          
    // Loop over test points
    Vector<double> s_plot(2);
    unsigned num_plot_points = 100;
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
     {

      // Get local coordinates of plot point
      double s_frac_along_edge=0.0;
      if (curved_edge==C1PlateHelper::CurvedEdgeEnumeration::zero)
       {
        s_plot[0]=0.0;
        s_plot[1]=1.0-double(iplot)/double(num_plot_points-1);
        s_frac_along_edge=1.0-s_plot[1];
       }
      else if (curved_edge==C1PlateHelper::CurvedEdgeEnumeration::one)
       {
        s_plot[0]=double(iplot)/double(num_plot_points-1);
        s_plot[1]=0.0;
        s_frac_along_edge=s_plot[0];
       }
      else if (curved_edge==C1PlateHelper::CurvedEdgeEnumeration::two)
       {
        s_plot[0]=1.0-double(iplot)/double(num_plot_points-1);
        s_plot[1]=1.0-s_plot[0];
        s_frac_along_edge=s_plot[1];
       }
      else
       {
        std::cout << "hierher never get here!" << std::endl;
        abort();
       }
      
      /// Position r as fct of zeta from curvilinear boundary representation
      zeta[0]=el_pt->bernadou_element_basis_pt()->get_s_ubar()+
       s_frac_along_edge*(el_pt->bernadou_element_basis_pt()->get_s_obar()-
                          el_pt->bernadou_element_basis_pt()->get_s_ubar());
      curviline_pt->position(zeta,r_from_boundary);
      
      
      // Get plot point
      Vector<double> interp_x(2, 0.0);
      el_pt->interpolated_x(s_plot, interp_x);
      
      // check
      double error=sqrt(pow(r_from_boundary[0]-interp_x[0],2)+
                        pow(r_from_boundary[1]-interp_x[1],2));
      if (error>max_error) max_error=error;
      
      if (plot_em)
       {
        unsigned prec=17;
        some_file << std::setprecision(prec) << r_from_boundary[0] << " "
                  << std::setprecision(prec) << r_from_boundary[1] << " "
                  << std::setprecision(prec) << interp_x[0] << " "
                  << std::setprecision(prec) << interp_x[1] << " "
                  << std::setprecision(prec) << zeta[0] << " "
                  << std::endl;
       }
     }
   }
 }
 
 oomph_info << "Max. error interpolated_x(): " << max_error << std::endl;
 
 
 if (plot_em)
  {
   some_file.close();
  }
 
 
 // oomph_info << "\n\nPlot of curved bell basis functions done! Now do: " << std::endl;
 // oomph_info << "cd RESLT"
 // oomph_info << "gnuplot -c ../validate_dpsidn_bubble.gp" << std::endl;
 // oomph_info << "gnuplot -c ../validate_dpsidn_nodal.gp" << std::endl;
 // oomph_info << "gnuplot -c ../validate_psi_nodal.gp" << std::endl;
 // oomph_info << "gnuplot -c ../validate_psi_bubble.gp" << std::endl;
 // oomph_info << "display validate*png" << std::endl;
 // oomph_info << std::endl;
 
}

   




/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////



//========================================================================
/// Validate mapping from monomials to 36 [66] basic dofs
//========================================================================
template<unsigned M>
void validate_monomials_to_basic_basis_functions(const std::string&
                                                 dir_name_for_output="")
{


 oomph_info << "\n\n\nTesting Bernadou basic basis functions for M = " << M << "\n"
            << "==================================================" << std::endl;

 // Make Bernadou element
 BernadouElementBasis<M>* b_pt=new BernadouElementBasis<M>;
 unsigned n_basic=b_pt->n_basic_basis_functions();

 // Test matrix
 DenseDoubleMatrix test_matrix(n_basic);
 
 ofstream some_file;
 char filename[100];
 
 // "Points" where interpolation property ought to be satisfied:
 std::map<std::string, // type of interpolation (a,b,...)
          Vector<      // instances of this (e.g. the three vertex points a
                       // end up in a Vector of length 3). 
           std::pair<  // for each instance store a pair:
            Vector<double>,     // first  part of pair stores the coordinates;
            Vector<std::string> // second part of pair stores the types of 
                                // dofs (value, deriv, ...); in general there
                                // are multiple ones (e.g. at the vertices we
                                // have 6) so we store them in a vector
            >>> test_point;
 
 test_point["a"].resize(3);
 test_point["a"][0]={{1.0,0.0},{"w","dwdx","dwdy","d2wdx2","d2wdxdy","d2wdy2"}};
 test_point["a"][1]={{0.0,1.0},{"w","dwdx","dwdy","d2wdx2","d2wdxdy","d2wdy2"}};
 test_point["a"][2]={{0.0,0.0},{"w","dwdx","dwdy","d2wdx2","d2wdxdy","d2wdy2"}};
 test_point["b"].resize(3);
 test_point["b"][0]={{0.0,0.5},{"-dwdx"}};
 test_point["b"][1]={{0.5,0.0},{"-dwdy"}};
 test_point["b"][2]={{0.5,0.5},{"dwdn"}};
 switch (M)
  {
  case 3:
   test_point["d"].resize(6);
   test_point["d"][0]={{0.0,0.75},{"w","-dwdx"}};
   test_point["d"][1]={{0.0,0.25},{"w","-dwdx"}};
   
   test_point["d"][2]={{0.25,0.0},{"w","-dwdy"}};
   test_point["d"][3]={{0.75,0.0},{"w","-dwdy"}};
   
   test_point["d"][4]={{0.75,0.25},{"w","dwdn"}};
   test_point["d"][5]={{0.25,0.75},{"w","dwdn"}};
   
   test_point["e"].resize(3);
   test_point["e"][0]={{0.5 ,0.25},{"w"}};
   test_point["e"][1]={{0.25,0.5 },{"w"}};
   test_point["e"][2]={{0.25,0.25},{"w"}};
   
   break;
   
  case 5:
   test_point["d"].resize(12);
   test_point["d"][0]={{0.0,5.0/6.0},{"w","-dwdx"}};
   test_point["d"][1]={{0.0,4.0/6.0},{"w","-dwdx"}};
   test_point["d"][2]={{0.0,2.0/6.0},{"w","-dwdx"}};
   test_point["d"][3]={{0.0,1.0/6.0},{"w","-dwdx"}};

   test_point["d"][4]={{1.0/6.0,0.0},{"w","-dwdy"}};
   test_point["d"][5]={{2.0/6.0,0.0},{"w","-dwdy"}};
   test_point["d"][6]={{4.0/6.0,0.0},{"w","-dwdy"}};
   test_point["d"][7]={{5.0/6.0,0.0},{"w","-dwdy"}};
  
   test_point["d"][8 ]={{5.0/6.0,1.0/6.0},{"w","dwdn"}};
   test_point["d"][9 ]={{4.0/6.0,2.0/6.0},{"w","dwdn"}};
   test_point["d"][10]={{2.0/6.0,4.0/6.0},{"w","dwdn"}};
   test_point["d"][11]={{1.0/6.0,5.0/6.0},{"w","dwdn"}};
   
   test_point["e"].resize(10);
   test_point["e"][0]={{1.0/6.0,4.0/6.0},{"w"}};
   test_point["e"][1]={{1.0/6.0,3.0/6.0},{"w"}};
   test_point["e"][2]={{1.0/6.0,2.0/6.0},{"w"}};
   test_point["e"][3]={{1.0/6.0,1.0/6.0},{"w"}};

   test_point["e"][4]={{2.0/6.0,1.0/6.0},{"w"}};
   test_point["e"][5]={{3.0/6.0,1.0/6.0},{"w"}};
   test_point["e"][6]={{4.0/6.0,1.0/6.0},{"w"}};
   
   test_point["e"][7]={{3.0/6.0,2.0/6.0},{"w"}};
   test_point["e"][8]={{2.0/6.0,3.0/6.0},{"w"}};
   
   test_point["e"][9]={{2.0/6.0,2.0/6.0},{"w"}};

   break;

  default:
   throw OomphLibError("Never get here!",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 Shape psi(n_basic); // basis functions themselves
 DShape dpsi(n_basic,2); // first derivs
 DShape d2psi(n_basic,3); // 2nd derivs xx, xy, yy 


// Test & plot 'em
bool plot_em=true;
if (dir_name_for_output=="") plot_em=false;
if (plot_em)
 {
  sprintf(filename,"%s/monomial_test_points.dat",
          dir_name_for_output.c_str());
  some_file.open(filename);
 }

// output intermediate results to screen
bool output_to_screen=false;

// Loop over the dofs
unsigned interpolation_condition_count=0;

// (class of dof: a,b,d,e)
for (auto dof_class : test_point)
 {

  if (output_to_screen)
   {
    oomph_info << std::fixed << std::setprecision(1);
   }
  
  // Loop over location of all dof locations of this class (a1,a2,a3,...)
  // dof_class.second is a vector containing the pairs of location and
  // quantities to be interpolated/checked
  unsigned count=0;
  for (auto dof_location_and_tests : dof_class.second)
   {
    // Tell us what you're doing
    if (output_to_screen)
     {
      oomph_info << dof_class.first << count << " : " << std::endl;
     }
    
    // List coordinates of test point (dof_location_and_tests.first is the vector
    // of coordinates)
    if (plot_em)
     {
      for (unsigned i=0;i<2;i++)
       {
        some_file << (dof_location_and_tests.first)[i] << " ";
       }
     }
    
    // Loop over test types: dof_location_and_tests.second
    // is the vector whose strings tell us what quantity
    // we're supposed to interpolate/test
    for (auto test_type : dof_location_and_tests.second)
      {       
       if (output_to_screen)
        {         
         oomph_info << test_type << " " << std::endl;
        }
       
       // Get all the basis functions and derivatives at this point
       b_pt->full_basic_polynomials(dof_location_and_tests.first,psi);
       b_pt->dfull_basic_polynomials(dof_location_and_tests.first,dpsi);
       b_pt->d2full_basic_polynomials(dof_location_and_tests.first,d2psi);

       if (test_type=="w")
        {
         for (unsigned i=0;i<n_basic;i++)
          {
           test_matrix(interpolation_condition_count,i)=psi[i];
           if (output_to_screen) oomph_info << psi[i] << " ";
          }
        }
       else if (test_type=="dwdx")
        {
         for (unsigned i=0;i<n_basic;i++)
          {
           test_matrix(interpolation_condition_count,i)=dpsi(i,0);
           if (output_to_screen) oomph_info << dpsi(i,0) << " ";
          }
        }
       else if (test_type=="dwdy")
        {
         for (unsigned i=0;i<n_basic;i++)
          {
           test_matrix(interpolation_condition_count,i)=dpsi(i,1);
           if (output_to_screen) oomph_info << dpsi(i,1) << " ";
          }
        }
       else if (test_type=="-dwdx")
        {
         for (unsigned i=0;i<n_basic;i++)
          {
           test_matrix(interpolation_condition_count,i)=-dpsi(i,0);
           if (output_to_screen) oomph_info << -dpsi(i,0) << " ";
          }
        }
       else if (test_type=="-dwdy")
        {
         for (unsigned i=0;i<n_basic;i++)
          {
           test_matrix(interpolation_condition_count,i)=-dpsi(i,1);
           if (output_to_screen) oomph_info << -dpsi(i,1) << " ";
          }
        }
       else if (test_type=="dwdn")
        {
         for (unsigned i=0;i<n_basic;i++)
          {
           test_matrix(interpolation_condition_count,i)=
            1.0/sqrt(2.0)*(dpsi(i,0)+dpsi(i,1));
           if (output_to_screen)
            {
             oomph_info << 1.0/sqrt(2.0)*(dpsi(i,0)+dpsi(i,1))
                        << " ";
            }
          }
        }
       else if (test_type=="d2wdx2")
        {
         for (unsigned i=0;i<n_basic;i++)
          {
           test_matrix(interpolation_condition_count,i)=d2psi(i,0);
           if (output_to_screen) oomph_info << d2psi(i,0) << " ";
          }
        }
       else if (test_type=="d2wdxdy")
        {
         for (unsigned i=0;i<n_basic;i++)
          {
           test_matrix(interpolation_condition_count,i)=d2psi(i,1);
           if (output_to_screen) oomph_info << d2psi(i,1) << " ";
          }
        }
       else if (test_type=="d2wdy2")
        {
         for (unsigned i=0;i<n_basic;i++)
          {
           test_matrix(interpolation_condition_count,i)=d2psi(i,2);
           if (output_to_screen) oomph_info << d2psi(i,2) << " ";
          }
        }
       else
        {
         throw OomphLibError("Never get here!",
                             OOMPH_CURRENT_FUNCTION,
                             OOMPH_EXCEPTION_LOCATION);
        }
       
       interpolation_condition_count++;
       if (output_to_screen) oomph_info << std::endl;
       
      }
    if (output_to_screen) oomph_info << std::endl;
    if (plot_em)
     {
      some_file << std::endl;
     }
    count++;
    

   }
 }
if (plot_em)
 {
  some_file.close();
 }

// Test: exactly one unit entry per row
unsigned n_interpolation_test=n_basic;
std::stringstream row_unit_ness_stream;
std::stringstream col_unit_ness_stream;
double tol=1.0e-10;
bool test_passed=true;
for (unsigned i=0;i<n_interpolation_test;i++)
 {
  unsigned count_one_in_row=0;
  unsigned count_non_one_in_row=0;
  unsigned count_one_in_col=0;
  unsigned count_non_one_in_col=0;
  for (unsigned j=0;j<n_interpolation_test;j++)
   {
    if (std::abs(test_matrix(i,j)    )<tol) count_non_one_in_row++;
    if (std::abs(test_matrix(i,j)-1.0)<tol)
     {
      count_one_in_row++;
      row_unit_ness_stream << "Unit entry in row " << i << " is ";
      if (i==j)
       {
        row_unit_ness_stream << " on diagonal" << std::endl;
       }
      else
       {
        row_unit_ness_stream << " off diagonal, namely in column "
                             << j << std::endl;
       }
          
     }
    if (std::abs(test_matrix(j,i)    )<tol) count_non_one_in_col++;
    if (std::abs(test_matrix(j,i)-1.0)<tol)
     {
      count_one_in_col++;         
      col_unit_ness_stream << "Unit entry in column " << j << " is ";
      if (i==j)
       {
        col_unit_ness_stream << " on diagonal" << std::endl;
       }
      else
       {
        col_unit_ness_stream << " off diagonal, namely in row "
                             << i << std::endl;
       }
     }

   }

  std::stringstream diagnostic;
  diagnostic << "Row "
             << i << " has "
             << count_one_in_row << " ones (should be 1) and "
             << count_non_one_in_row << " non-ones (should be "
             << n_interpolation_test-1 << ") "
             << std::endl;
     
  if ((count_one_in_row!=1)||
      (count_one_in_col!=1)||
      (count_non_one_in_row!=(n_interpolation_test-1))||
      (count_non_one_in_col!=(n_interpolation_test-1)))
   {
    test_passed=false;
    oomph_info << "failed: " << diagnostic.str();
   }
  else
   {
    oomph_info << "passed: " << diagnostic.str();
   }
 }


oomph_info << row_unit_ness_stream.str() << std::endl << std::endl;
oomph_info << col_unit_ness_stream.str() << std::endl << std::endl;
   
if (test_passed)
 {
  oomph_info << "Test of basic basis functions passed!"
             << std::endl;
 }
else
 {
  oomph_info << "Test of basic basis functions failed!"
             << std::endl;
 }
 


// Plot all basis functions
if (plot_em)
 {
  // Tecplot header info from some generic triangle element
  TElement<2,2>* aux_el_pt= new TElement<2,2>;
  unsigned nplot=100;
  for (unsigned i=0;i<n_basic;i++)
  { 
   sprintf(filename,"%s/test_basic_basis%i.dat",
           dir_name_for_output.c_str(),i);
   some_file.open(filename);
   
   // Tecplot header info
   some_file << aux_el_pt->tecplot_zone_string(nplot);
   
   // Loop over plot points
   Vector<double> s_plot(2);
   unsigned num_plot_points = aux_el_pt->nplot_points(nplot);
   for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
     // Get local coordinates of plot point
     aux_el_pt->get_s_plot(iplot, nplot, s_plot);
     
     Shape psi(n_basic);
     b_pt->full_basic_polynomials(s_plot,psi);
     DShape dpsi(n_basic,2); // first derivs
     b_pt->dfull_basic_polynomials(s_plot,dpsi);
     DShape d2psi(n_basic,3); // 2nd derivs xx, xy, yy
     b_pt->d2full_basic_polynomials(s_plot,d2psi);
     
     some_file << s_plot[0] << " "
               << s_plot[1] << " ";
     some_file << psi[i] << " ";
     some_file << dpsi(i,0) << " "
               << dpsi(i,1) << " ";
     some_file << d2psi(i,0) << " "
               << d2psi(i,1) << " "
               << d2psi(i,2) << " ";
     some_file << std::endl;
    }
   
   // Write tecplot footer (e.g. FE connectivity lists)
   aux_el_pt->write_tecplot_zone_footer(some_file, nplot);
   some_file.close();
  }
  
  delete aux_el_pt;
  aux_el_pt=0;
  
  oomph_info << "\n\nPlot of basis functions done! Now do: " << std::endl;
  oomph_info << "oomph-convert -z test_basic_basis*dat" << std::endl;
  oomph_info << "makePvd test_basic_basis test_basic_basis.pvd" << std::endl;
  oomph_info << "oomph-convert -p2 test_points.dat " << std::endl;
  oomph_info << "paraview --state test_basic_basis.pvsm " << std::endl;
  oomph_info << std::endl;
 }
  

}



////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


// Run problem with specified boundary order, for square domain (if bool is true)
// and rotation angle (only used for square domain)
void problem_level_test(const unsigned& boundary_order,
                        bool use_square_domain,
                        const double& phi)
{
 
 
 // (only used for non-square domain)
 unsigned m_poly_actual_boundary=5;
 
 // Build problem  
 UnstructuredC1PlateProblem<FoepplVonKarmanC1CurvableBellElement<4>> problem(
  Parameters::Element_area,m_poly_actual_boundary,boundary_order,
  use_square_domain, phi);
  
 
 // Document the initial state
 problem.doc_solution();
 
 // Test curved Bell
 std::string dir_name="RESLT"; // hierher incorporate rotation
 problem.validate_curved_bell_and_bubble_basis_functions(dir_name);
 
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

 

//=======start_of_main========================================
///Driver code 
//============================================================
int main(int argc, char** argv)
{
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);


  // Store command line arguments
  CommandLineArgs::setup(argc, argv);
  
  // Parse command line
  CommandLineArgs::parse_and_assign();

  // Doc what has actually been specified on the command line
  //CommandLineArgs::doc_specified_flags();



  // Test 1: From the very bottom: Monomials are OK
  validate_monomials_to_basic_basis_functions<5>();
  validate_monomials_to_basic_basis_functions<3>();


  // Loop over boundary order
  for (unsigned b=3;b<6;b+=2)
   {
    // Loop over angle
    double phi=0.0;
    for (unsigned i_angle=0;i_angle<2;i_angle++)
     {
      // Loop over rotate dofs
      for (unsigned rotate=0;rotate<2;rotate++)
       {
        if (rotate==1)
         {
          Parameters::Rotate_coordinates_on_all_curvilinear_boundaries=true;
         }
        else
         {
          Parameters::Rotate_coordinates_on_all_curvilinear_boundaries=false;
         }
         
        
        bool use_square_domain=true;      
        problem_level_test(b,
                           use_square_domain,
                           phi);
       }
     }
   }
  
       

  // // Test 3: From the very top: interpolated_x
  // {

  //  // Check accuracy of newton solver when determining local coordiante
  //  // of point on curvilinear boundary
  //  C1CurviLine::Tol_for_get_zeta=1.0e-12;

  //  // Allow massively warped elements
  //  FiniteElement::Accept_negative_jacobian=true;

   
  //  for (unsigned boundary_order=3;boundary_order<=5;boundary_order+=2)
  //   {
  //    for (unsigned m_poly_actual_boundary=2;
  //         m_poly_actual_boundary<=boundary_order+5;
  //         m_poly_actual_boundary++)
  //     {
  //      UnstructuredC1PlateProblem<FoepplVonKarmanC1CurvableBellElement<4>> problem(
  //       Parameters::Element_area,m_poly_actual_boundary,boundary_order);
       
  //      oomph_info << "Testing with m_poly_actual_boundary = " << m_poly_actual_boundary
  //                 << " ; boundary_order = " << boundary_order << " : ";
  //      problem.validate_interpolated_x(".",m_poly_actual_boundary, boundary_order);
  //     }
  //   }
       
       
  //  exit(0);
  // }

} // End of main
