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

//Generic routines
#include "generic.h"

// The mesh
#include "meshes/triangle_mesh.h"

// The equations
#ifdef USE_KS
#include "c1_koiter_steigmann.h"
#else
#include "c1_foeppl_von_karman.h"
#endif

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;



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
                               const Vector<double>& right) 
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
   Left[0]=left[0];
   Left[1]=left[1];
   Right.resize(2);
   Right[0]=right[0];
   Right[1]=right[1];
   
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
 
 /// Enumeration of cases
 enum
 {
  Clamped_validation,
  Pinned_validation,
  Balance_on_edge,
  Free_edges
 };

 /// Rotate coordinates on curvilinear boundaries?
 bool Rotate_coordinates_on_all_curvilinear_boundaries=true;

 /// Which case are we doing
 unsigned Problem_case = Free_edges;
 
 /// Ellipse half x-axis
 double A = 1.0;
 
 /// Ellipse half y-axis
 double B = 1.0;

 /// Damping constant for damped solves (magnitude sort of irrelevant
 /// since the adaptive timestepping will kick in anyway).
 double Mu =1.0; 


 /// Nondimensional thickness of plate -- dependent parameter compute!
 double Thickness = 0.0;

 #ifdef USE_KS

 // hierher update these to make them consistent with fvk
 
 /// Membrane coupling coefficient 
 double Eta_u = 1.0;

  /// hierher what is this?
 double Eta_sigma = 1.0; 

#else
 
 /// Membrane coupling coefficient (a dependent parameter)
 double Eta = 0.0; // hierher does it have the 1-nu^2 in it?)
                   // 12.0 * (1.0 - Nu * Nu) / (Thickness * Thickness);

 #endif

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

 #ifdef USE_KS
 
  /// Traction depending on the position (x,y) and deformation of the sheet
  void get_traction(const Vector<double>& x,
		    const Vector<double>& u,
		    const DenseMatrix<double>& grad_u,
		    const Vector<double>& n,
		    Vector<double>& traction)
  {
    // Metric tensor of deformed surface
    DenseMatrix<double> G(2,2,0.0);
    for(unsigned alpha = 0; alpha < 2; alpha++)
    {
      G(alpha, alpha) += 1.0;
      for (unsigned beta = 0; beta < 2; beta++)
      {
        G(alpha, beta) += grad_u(alpha, beta) + grad_u(beta, alpha);
	for (unsigned i = 0; i < 3; i++)
	{
          G(alpha, beta) += grad_u(i, alpha) * grad_u(i, beta);
	}
      }
    }



    // hierher Aidan: do we really need this conversion? Lagr/Eulerian. why?

    // hierher pressure --> traction in src too

    // hierher: scale KS like FvK otherwise we'll all go insane!
    
    // Find the pressure per undeformed area in terms of the pressure per
    // deformed area
    double p = sqrt(G(0,0)*G(1,1) - G(1,0)*G(0,1)) * P_mag/
     (12.0 * (1.0 - Nu * Nu) / (Thickness * Thickness));
    
    // Assign traction
    traction.resize(3);
    traction[0] = p * n[0];
    traction[1] = p * n[1];
    traction[2] = p * n[2];

    // Dead load
    if (Parameters::Problem_case == Parameters::Balance_on_edge)
     {
      traction[0] = 0.0;
      traction[1] = 0.0;
      traction[2] = p;
     }
  }

#else
 
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

 #endif



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
  UnstructuredC1PlateProblem(double const& element_area = 0.09);

  /// Destructor
  ~UnstructuredC1PlateProblem()
  {
    // Close trace file
    Trace_file.close();
  };

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

 #ifndef USE_KS

 // hierher check dofs for KS; it this type of pinning still correct?
 // and how do we make the equations linear?
 
 /// Make the problem linear (biharmonic) by pinning all in-plane dofs and
 /// setting eta=0; also readjusts the constraints and and reassigns
 /// the equation numbers
 void make_linear()
  {
   // Remove stretching coupling
   Parameters::Eta = 0.0;
   
   // Pin all in-plane displacements
   unsigned n_node = Bulk_mesh_pt->nnode();
   for(unsigned i_node = 0; i_node < n_node; i_node++)
    {
     Bulk_mesh_pt->node_pt(i_node)->pin(0);
     Bulk_mesh_pt->node_pt(i_node)->set_value(0,0.0);
     Bulk_mesh_pt->node_pt(i_node)->pin(1);
     Bulk_mesh_pt->node_pt(i_node)->set_value(1,0.0);
    }
   
   
   // Update the corner constraints based on boundary conditions
   // after changing the boundary conditions
   unsigned n_el = Constraint_mesh_pt->nelement();
   for(unsigned i_el = 0; i_el < n_el; i_el++)
    {
     dynamic_cast<DuplicateNodeConstraintElement*>
      (Constraint_mesh_pt->element_pt(i_el))
      ->pin_redundant_constraints();
    }
   
   // Reassign the equation numbers
   oomph_info << "Reassiging equation numbers after changing BCs. "
              << " ndof = " << assign_eqn_numbers() << std::endl;
   
  } // End make_linear
 
#endif
 
 
 /// Doc the solution
 void doc_solution(bool steady = true);

 

 /// Validate mapping from monomials to 36 [66] basic dofs
 template<unsigned M>
 void validate_monomials_to_basic_basis_functions(const std::string&
                                                  dir_name_for_output="");

 /// Validate all basis functions for curved bell
 void validate_curved_bell_and_bubble_basis_functions(const std::string&
                                                     dir_name_for_output);
 

 // // hierher 
 // /// Validate interpolation of normal derivative along curved edge
 // template<unsigned M>
 // void validate_dpsi_dn_along_edge(const std::string&
 //                                  dir_name_for_output="");
 

 
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



 /// Typedef for the function pointer to the function that allows
 /// documtating the progress of the damped solve
 typedef void (*DampedSolveDocSolutionFctPt)(const unsigned&);

 
// hierher move into base class
/// Use damped solves to get close to a steady solution; when close
/// enough, attempt a steady solve. If that fails, be stricter about the
/// meaning of "close enough" and repeat until a steady solve succeeds.
///
/// Expects:
///   dt_supplied_guess: a guess for a good timestep size
///   epsilon:           an 'error tolerance' for the timestepper to limit
///                        the size of a damped step
///  doc_soln_fct_pt:    function pointer to void function that takes
///                      unsigned (representing the number of the current
///                      damped solve) as arg, This function usually calls
///                      the doc_solution(...) fct of the underlying problem
///                      class. Defaults to null, in which case no doc is
///                      produced.
/// Returns:
///   a double suggesting the timestep dt for the next solve
double damped_solve(
 const double& dt_supplied_guess, 
 const double& epsilon, 
 DampedSolveDocSolutionFctPt doc_soln_fct_pt=0)
  {
   // We are unsteady until a steady solve succeeds
   bool dofs_are_steady = false;

   // This used to be an input parameter but generally it's too wobbly
   // so let's set it to false here; can re-enable if it's ever found to be
   // useful
   bool begin_with_steady_solve=false;
   
   // Max residual of the steady problem before we attempt a steady solve
   double sufficiently_small = 1.0e-2;
   
   // Timestep size
   double dt = dt_supplied_guess;
   
  // Try steady indicates when we should attempt a steady solve
  bool try_steady = begin_with_steady_solve;
  
  // Value to be returned for initial guess for next damped_solve dt
  double suggested_dt_for_next_damped_solve = dt;
  
  // Only set dt_initial_guess once, after the first successful solve
  bool suggested_dt_for_next_damped_solve_is_unset = true;

  // Counter for unsteady solves
  unsigned unsteady_solve_counter=0;
  
  // If we are documenting the damped stage, create an initial state before any
  // unsteady solves have been done.
  if (doc_soln_fct_pt!=0)
   {
    doc_soln_fct_pt(unsteady_solve_counter);
    unsteady_solve_counter++;
   }
  
  
  // Keep looping until all the dofs are steady
  while (!dofs_are_steady)
   {
    //------------------------------------------------------------------------
    // If we are supposed to try a steady solve, do it
    if (try_steady)
     {
      oomph_info << "ATTEMPT A STEADY SOLVE" << std::endl;
      
      // Store the dofs before a steady solve so that they can be put back in
      // case it fails
      store_current_dof_values();

      // Get the max residual in case we need it to adjust sufficiently_small
      DoubleVector res;
      get_residuals(res);
      double max_steady_residual = res.max();
      try
       {
        // <<< Solve >>>
        steady_newton_solve();
        
        // If that worked, we have achieved steady state.
        // Celebrate and take note
        oomph_info << "\nHOORAY\n"
                   << "Steady solve was successful, damped solve complete\n"
                   << std::endl;
        dofs_are_steady = true;
        
        // If we are documenting the unsteady states, add the final solution to
        // the unsteady solution outputs
        if (doc_soln_fct_pt!=0)
         {
          doc_soln_fct_pt(unsteady_solve_counter);
          unsteady_solve_counter++;
         }
       }
      // If the steady solve fails, we need to tidy up before carrying on
      catch (OomphLibError& error)
       {
        // If our tolerance to attempt a steady solve is smaller than the
        // tolerance, then this implies that the initial residual was within
        // tolerance and we still got an error! Odd.
        if (sufficiently_small < newton_solver_tolerance())
         {
          oomph_info << "\nUH OH\n"
                     << "\"sufficiently small\" is now " << sufficiently_small
                     << " which is smaller than the Newton solver tolerance.\n"
                     << "For some reason we still gt an error in the"
                     << " steady Newton solve. \n"
                     << "Giving up on damped solves..." << std::endl;
          throw error;
         } // End of if tolerance is to small
        else
         {
          oomph_info << "\nNOT STEADY ENOUGH.\n"
                     << "\"sufficiently_small\" is insufficiently small \n"
                     << "i.e. we've stopped the unsteady solves too early\n"
                     << "Decreasing it from " << sufficiently_small
                     << " to " << max_steady_residual / 2.0 << ".\n"
                     << "Continuing with damped solves...\n"
                     << std::endl;
          
          // Decrease the threshold for attempting steady solves as this one
          // didn't work
          sufficiently_small = max_steady_residual / 2.0;
          
          // Go back to the state we were in before attempting the steady solve
          restore_dof_values();
          for (unsigned i = 0; i < ntime_stepper(); i++)
           {
            time_stepper_pt(i)->undo_make_steady();
           }
          
          // Stop trying steady solves
          try_steady = false;
          
          // Keep calm and carry on // hierher Aidan: what is this?
          error.disable_error_message();
          
         } // End of else tolerance is not too small
       } // End of catch error
     } // End of if try_steady
    
    //------------------------------------------------------------------------
    // Try get us close to a steady solution by solving the damped version of
    // the equations. When it is time to try a steady solve, break this loop.
    while(!try_steady)
     {
      //----------------------------------------------------------------------
      // Begin by doing a damped solve
      oomph_info << "NEW DAMPED PSEUDO-TIME STEP WITH: dt = "
                 << dt << std::endl;
      double dt_next = adaptive_unsteady_newton_solve(dt, epsilon);
      dt = dt_next;
      
      // If we haven't set the initial guess for the next damped solve dt, then
      // set it now. It should be the recommended timestep after the first
      // successful solve. Assuming the following damped solve will start in a
      // roughly similar state to this one, this is appropriate.
      if (suggested_dt_for_next_damped_solve_is_unset)
       {
        suggested_dt_for_next_damped_solve = dt_next;
        suggested_dt_for_next_damped_solve_is_unset = false;
       }
      
      // If we are documenting the unsteady solutions then do so, else just
      // just increase the unsteady step counter to keep count of damped steps
      if (doc_soln_fct_pt!=0)
       {
        doc_soln_fct_pt(unsteady_solve_counter);
        unsteady_solve_counter++;
       }
      
      //------------------------------------------------------------------------
      // Check how close we are to a steady solution by getting the steady
      // max residual, if it is sufficiently small, try a steady solve.
      // If that doesn't work, restrict what it means to be "sufficiently small"
      // and return to unsteady. We repeat this until the steady solve works,
      // or, we give up.
      
      // First set the timesteppers to steady
      for (unsigned i = 0; i < ntime_stepper(); i++)
       {
        time_stepper_pt(i)->make_steady();
       }
      
      // Then get the residual
      DoubleVector res;
      get_residuals(res);
      double max_steady_residual = res.max();
      oomph_info << std::endl
                 << "The max steady residual is " << max_steady_residual
                 << std::endl;
      
      // If it is "sufficiently small" then try a steady solve
      try_steady = max_steady_residual < sufficiently_small;
      
      // Reset time steppers
      for (unsigned i = 0; i < ntime_stepper(); i++)
       {
        time_stepper_pt(i)->undo_make_steady();
       }
     } // End of while(!try_steady)
   } // End of while(!steady)
  
  
  // Done; return most recent suggestion for timestep
  return suggested_dt_for_next_damped_solve;
  
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

 
 /// Global temporal error norm for pseudo-timestepping
 double global_temporal_error_norm()
  {
#ifdef USE_KS
   oomph_info << "Find w for KS; also fix if statement below" << std::endl;
   abort(); // hierher
#else
   unsigned w_dof_index=2;
#endif
   
   double global_error = 0.0;
   
   //Find out how many nodes there are in the problem
   unsigned n_node = Bulk_mesh_pt->nnode();
   
   //Loop over the nodes and calculate the estimated error in the values
   for(unsigned i=0;i<n_node;i++)
    {
     // Node with only in-plane displacements?
     unsigned nval=Bulk_mesh_pt->node_pt(i)->nvalue();
     if (nval>2)
      {
       // Get error in solution: Difference between predicted and actual
       // value
       double error = Bulk_mesh_pt->node_pt(i)->time_stepper_pt()->
        temporal_error_in_value(Bulk_mesh_pt->node_pt(i),w_dof_index);
       
       //Add the square of the individual error to the global error
       global_error += error*error;
      }
    }
   
   // Divide by the number of nodes
   global_error /= double(n_node);
   
   // Return square root...
   return sqrt(global_error);
   
  } // end of global_temporal_error_norm

 
  /// Pin all displacements and rotation at the centre
  void pin_all_displacements_and_rotation_at_centre_node();

  /// Balance on edge along line
  void pin_for_balance_on_edge();
 
    
  /// Trace file to document norm of solution
  ofstream Trace_file;

  /// Pointer to "bulk" mesh
  TriangleMesh<ELEMENT>* Bulk_mesh_pt;

 
  /// Enumeration to keep track of boundary ids
  enum
  {
    Outer_boundary0 = 0,
    Outer_boundary1 = 1,
    Outer_boundary2 = 2,
    Outer_boundary3 = 3,
    Inner_boundary0 = 4,
    Inner_boundary1 = 5,
    Inner_boundary2 = 6
  };

  /// Target element area
  double Element_area;

  /// Pointer to constraint mesh
  Mesh* Constraint_mesh_pt;

  /// Doc info object for labeling output
  DocInfo Doc_info;

 // The Line Visualiser.
 LineVisualiser* LV_pt;


}; // end_of_problem_class



//======================================================================
/// Constructor definition
//======================================================================
template<class ELEMENT>
UnstructuredC1PlateProblem<ELEMENT>::UnstructuredC1PlateProblem(const double&
                                                                element_area)
 : Element_area(element_area)
{

 // Allocate the timestepper only used in anger for damped solve
 add_time_stepper_pt(new BDF<1>); // (true)); // hierher adaptive


 // Build the mesh
 //================
 
 
 //Outer boundary
 //--------------
 
 double A = Parameters::A;
 double B = Parameters::B;
 Ellipse* outer_boundary_ellipse_pt = new Ellipse(A, B);
 
 // Storage for outer boundaries (for triangle)
 Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(4);

 //First bit
 double zeta_start = 0.0;
 double zeta_end = 0.5*MathematicalConstants::Pi;
 unsigned nsegment = (unsigned)(MathematicalConstants::Pi/sqrt(Element_area));
 outer_curvilinear_boundary_pt[0] = 
  new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
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
 
 bool flatten_one_side=false; 
 if (Parameters::Problem_case==Parameters::Balance_on_edge)
  {
   flatten_one_side=false;
  }
 if (flatten_one_side)
  {
   zeta_start = 0.0;
   zeta_end = 1.0;
   outer_curvilinear_boundary_pt[1] =
    new TriangleMeshCurviLine(straight_line_pt, zeta_start,
                              zeta_end, nsegment, Outer_boundary1);
  }
 else
  {
   zeta_start = 0.5*MathematicalConstants::Pi;
   zeta_end = MathematicalConstants::Pi;
   outer_curvilinear_boundary_pt[1] =
    new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
                              zeta_end, nsegment, Outer_boundary1);
  }

   
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
  TriangleMeshClosedCurve* outer_boundary_pt =
  new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);

  
  // Internal open boundaries
  //-------------------------
  
  // Represent inner boundaries by curvilines?
  bool use_curviline=true;
  if (CommandLineArgs::command_line_flag_has_been_set
      ("--use_polyline_for_internal_boundaries"))
   {
    use_curviline=false;
   }

  
  // We want internal open curves
  Vector<TriangleMeshOpenCurve *> inner_open_boundaries_pt;
  
  // Internal bit 
  
  // Open curve 1
  Vector<Vector<double> > vertices(3,Vector<double>(2,0.0));
  vertices[0][0] =-0.5*Parameters::A;
  vertices[0][1] = 0.0;
  
  vertices[1][0] = 0.0;
  vertices[1][1] = 0.0;
  
  vertices[2][0] = 0.5*Parameters::A;
  vertices[2][1] = 0.0;
  unsigned boundary_id = Inner_boundary0;
  
  TriangleMeshCurveSection* boundary2_pt=0;
  if (use_curviline)
   {
    // Straight curvilinear line
    TwoDStraightLineFromTwoPoints* straight_line_pt =
     new TwoDStraightLineFromTwoPoints(vertices[0],vertices[2]);
    
    double zeta_start=0.0;
    double zeta_end=1.0;
    
    // Two segments to mimick the three-vertex polyline version
    unsigned nsegment=2;
    boundary2_pt =
     new TriangleMeshCurviLine(straight_line_pt, zeta_start,
                               zeta_end, nsegment, boundary_id);
    
   }
  else
   {
    boundary2_pt =
     new TriangleMeshPolyLine(vertices, boundary_id);
   }
  
  // Each internal open curve is defined by a vector of
  // TriangleMeshCurveSections
  Vector<TriangleMeshCurveSection *> internal_curve_section1_pt(1);
  internal_curve_section1_pt[0] = boundary2_pt;
  
  // The open curve that defines this boundary
  inner_open_boundaries_pt.push_back(
   new TriangleMeshOpenCurve(internal_curve_section1_pt));
  
  
  // Make a T-shape in the middle of the domain (to force splitting
  // of elements)
  bool make_t_shape=false;
  if (CommandLineArgs::command_line_flag_has_been_set
      ("--use_t_shape_internal_boundaries"))
   {
    make_t_shape=true;
   }
  if (make_t_shape)
   {
    
    // Open Curve 2
    Vector<Vector<double> > vertices(2,Vector<double>(2,0.0));
    vertices[0][0] = 0.0;
    vertices[0][1] =-0.5;
    
    vertices[1][0] = 0.0;
    vertices[1][1] = 0.0;
    boundary_id = Inner_boundary1;

    TriangleMeshCurveSection* boundary3_pt=0;
    if (use_curviline)
     {
      // Straight curvilinear line
      TwoDStraightLineFromTwoPoints* straight_line_pt =
       new TwoDStraightLineFromTwoPoints(vertices[0],vertices[1]);
       
      double zeta_start=0.0;
      double zeta_end=1.0;

      // Just one segment to mimick the two-vertex polyline version
      unsigned nsegment=1;
      boundary3_pt =
       new TriangleMeshCurviLine(straight_line_pt, zeta_start,
                                 zeta_end, nsegment, boundary_id);
            
      // Connect final vertex on this boundary
      // to middle of the horizontal one:
      double zeta_to_connect_to=0.5;
      boundary3_pt->connect_final_vertex_to_curviline(
       dynamic_cast<TriangleMeshCurviLine*>(boundary2_pt),
       zeta_to_connect_to);
     }
    else
     {
      boundary3_pt =
       new TriangleMeshPolyLine(vertices, boundary_id);
      
      // Connect final vertex on this boundary
      // to middle vertex in the horizontal one:
      unsigned vertex_to_connect_to=1;
      boundary3_pt->connect_final_vertex_to_polyline(
       dynamic_cast<TriangleMeshPolyLine*>(boundary2_pt),
       vertex_to_connect_to);
     }

    // Each internal open curve is defined by a vector of
    // TriangleMeshCurveSections
    Vector<TriangleMeshCurveSection*> internal_curve_section2_pt(1);
    internal_curve_section2_pt[0] = boundary3_pt;
    
    // The open curve that defines this boundary
    inner_open_boundaries_pt.push_back(
     new TriangleMeshOpenCurve(internal_curve_section2_pt));


    // hierher this creates a node that is on three boundaries and (currently overwhelms our lovely
    // little (and limited-scope) black box helper function:
    
    // // Open Curve 3
    // vertices[0][0] = 0.0;
    // vertices[0][1] = 0.5;
    
    // vertices[1][0] = 0.0;
    // vertices[1][1] = 0.0;
    // boundary_id = Inner_boundary2;
    
    // TriangleMeshPolyLine* boundary4_pt =
    //  new TriangleMeshPolyLine(vertices, boundary_id);
     
    // // Connect final vertex on this boundary
    // // to middle vertex in the horizontal one:
    // unsigned vertex_to_connect_to=1;
    // boundary4_pt->connect_final_vertex_to_polyline(
    //  boundary2_pt,
    //  vertex_to_connect_to);
  
    // // Each internal open curve is defined by a vector of
    // // TriangleMeshCurveSections
    // Vector<TriangleMeshCurveSection *> internal_curve_section3_pt(1);
    // internal_curve_section3_pt[0] = boundary4_pt;
    
    // // The open curve that defines this boundary
    // inner_open_boundaries_pt.push_back(
    //  new TriangleMeshOpenCurve(internal_curve_section3_pt));
     
   }


  
  //Create mesh parameters object
  TriangleMeshParameters mesh_parameters(outer_boundary_pt);

  // Element area
  mesh_parameters.element_area() = Element_area;

  // Specify the internal open boundaries
  mesh_parameters.internal_open_curves_pt() = inner_open_boundaries_pt;

  // Build an assign bulk mesh
  Bulk_mesh_pt=new TriangleMesh<ELEMENT>(mesh_parameters,
                                         time_stepper_pt());


  
  // Now upgrade to (potentially) curved C1 boundaries 
  //==================================================
  {

  
   // Let's have a look at the orig mesh
   Bulk_mesh_pt->output("mesh_before_black_box_upgrade.dat");




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
    ("duplicated_nodes.dat");
   C1PlateHelper::Upgraded_to_curved_edge_element_stream.open
    ("elements_upgraded_to_curved.dat");
   C1PlateHelper::Split_elements_output_stream.open
    ("split_elements.dat");
   C1PlateHelper::Rotated_node_output_stream.open
    ("rotated_nodes.dat");
   C1PlateHelper::Rotated_element_output_stream.open
    ("rotated_elements.dat");
 
   
   // Rotate coordinates on curvilinear boundaries?
   C1PlateHelper::upgrade_triangle_mesh_for_c1_plate_bending<ELEMENT>(
    Bulk_mesh_pt,
    Constraint_mesh_pt,
    Parameters::Rotate_coordinates_on_all_curvilinear_boundaries);

   // Done
   C1PlateHelper::Duplicated_node_output_stream.close();
   C1PlateHelper::Upgraded_to_curved_edge_element_stream.close();
   C1PlateHelper::Split_elements_output_stream.close();
   C1PlateHelper::Rotated_node_output_stream.close();
   C1PlateHelper::Rotated_element_output_stream.close();

   // Let's have a look at the new mesh
   Bulk_mesh_pt->output("mesh_black_box_upgrade.dat");
   doc_boundary_coords();
   std::string name_prefix="test_";
   doc_boundary_elements_and_faces(Bulk_mesh_pt,name_prefix);



   
  
  }

  // Build global mesh
  //==================
  
  //Add submeshes to problem
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Constraint_mesh_pt);

  // Combine submeshes into a single Mesh 
  build_global_mesh();


 
  // Complete the build of all elements so they are fully functional
  //================================================================
  unsigned n_element = Bulk_mesh_pt->nelement();
  for(unsigned e=0;e<n_element;e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    //Set the traction and physical constants
#ifdef USE_KS
    
    // hierher: rename pressure --> traction in src
    el_pt->pressure_fct_pt() = &Parameters::get_traction;

    el_pt->mu_pt()=&Parameters::Mu;

    // hierher why do we need thickness and (two!) etas?
    el_pt->thickness_pt() = &Parameters::Thickness;
    el_pt->nu_pt() = &Parameters::Nu;
    el_pt->eta_u_pt() = &Parameters::Eta_u;
    el_pt->eta_sigma_pt() = &Parameters::Eta_sigma;

    // Damping parameter for damped solve
    el_pt->mu_pt() = &Parameters::Mu;

#else
    
    el_pt->pressure_fct_pt() = &Parameters::get_pressure;

    // Damping parameter for damped solve
    el_pt->mu_pt()=&Parameters::Mu;
    
    el_pt->nu_pt() = &DimensionalParameters::Poisson_ratio;
    el_pt->eta_pt() = &Parameters::Eta;

#endif
  }
  

  
  // Set the boundary conditions
  //============================
  
  // Get map of curvline boundaries in the mesh
  std::map<unsigned, TriangleMeshCurviLine*> curviline_boundary_pt =
   Bulk_mesh_pt->curviline_boundary_pt();
  
  
  // Clamp it
  if (Parameters::Problem_case==Parameters::Clamped_validation)
  {
   
   // Create vector of pointers to BoundaryConditionForC1PlateBending objects
   // that specify the imposed displacements along the boundary
   Vector<BoundaryConditionForC1PlateBending*> boundary_values_pt(3);
   boundary_values_pt[0]= new Parameters::ZeroC0BoundaryConditions;
   boundary_values_pt[1]= new Parameters::ZeroC0BoundaryConditions;
   boundary_values_pt[2]= new Parameters::ZeroC1BoundaryConditions;
   
    // Set the boundary conditions on the outer boundaries
   unsigned nbound = 4;
   for(unsigned b = 0; b < nbound; b++)
    {
     const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);
     for(unsigned e=0;e<nb_element;e++)
      {
       // Get pointer to bulk element adjacent to b
       ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b,e));

       // Clamp: i.e. pin the two in-plane displacements, and pin the out-of-plane
       // displacement and its normal derivative. We also apply implied
       // boundary conditions (e.g. specification of dw/dn also implies
       // d^2w/dn/dzeta etc.)
       el_pt->fully_clamp_specified_boundary(b,boundary_values_pt,
                                             curviline_boundary_pt[b]);
      }
    }
  }

  // Pin it 
  else if (Parameters::Problem_case==Parameters::Pinned_validation)
  {
   
   // Create vector of pointers to BoundaryConditionForC1PlateBending objects
   // that specify the imposed displacements along the boundary
   Vector<BoundaryConditionForC1PlateBending*> boundary_values_pt(3);
   boundary_values_pt[0]= new Parameters::ZeroC0BoundaryConditions;
   boundary_values_pt[1]= new Parameters::ZeroC0BoundaryConditions;
   boundary_values_pt[2]= new Parameters::ZeroC0BoundaryConditions;
   
    // Set the boundary conditions on the outer boundaries
   unsigned nbound = 4;
   for(unsigned b = 0; b < nbound; b++)
    {
     const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);
     for(unsigned e=0;e<nb_element;e++)
      {
       // Get pointer to bulk element adjacent to b
       ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b,e));

       // Pin i.e. pin the in-plane and out-of plane displacements
       // We also apply implied
       // boundary conditions (e.g. specification of w also implies
       // dw/dt and d^2w/dt^2 etc.
       el_pt->pin_specified_boundary(b,boundary_values_pt,
                                     curviline_boundary_pt[b]);
      }
    }
  }
  // Balance on edge along line
  else if (Parameters::Problem_case==Parameters::Balance_on_edge)
   {
    pin_for_balance_on_edge();
   }
  // All other cases: simply pin and stop rotation via the centre
  else if (Parameters::Problem_case==Parameters::Free_edges)
   {
    pin_all_displacements_and_rotation_at_centre_node();
   }
  else
   {
    oomph_info << "Never get here" << std::endl;
    abort();
   }
  
   
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
             << assign_eqn_numbers() << '\n';

  
  // Set directory
  Doc_info.set_directory("RESLT");

  // Open trace file
  char filename[100];
  sprintf(filename, "RESLT/trace.dat");
  Trace_file.open(filename);


  
  // Setup sample points for line visualiser
  unsigned  npt=100;
  Vector<Vector<double> > coord_vec(npt);
  coord_vec[0].resize(2);
  coord_vec[0][0]=0.0;
  coord_vec[0][1]=0.0;
  for (unsigned j=1;j<npt;j++)
   {
    coord_vec[j].resize(2);
    coord_vec[j][0]=double(j)/double(npt);
    coord_vec[j][1]=0.0;
   }

  
  // The C1 elements don't provide a standard implementation
  // of dshape(s,...) so locate_zeta has to use finite
  // differencing. Call this before setting up the
  // line visualiser; this is where the locate_zeta happens!
  Locate_zeta_helpers::Evaluate_dzeta_ds_by_fd=true;

  
  // Setup line visualiser
  LV_pt=new LineVisualiser(Bulk_mesh_pt,
                           coord_vec);


  
} // end Constructor






//==start_of_pin_all_displacements_and_rotation_at_centre_node================
/// pin all displacements and rotations in the centre
//============================================================================
template<class ELEMENT>
void UnstructuredC1PlateProblem<ELEMENT>::
pin_all_displacements_and_rotation_at_centre_node()
{
 
 // Choose non-centre node on which we'll supress
 // the rigid body rotation around the z axis.
 double max_x_potentially_pinned_node=-DBL_MAX;
 Node* pinned_rotation_node_pt=0;
 Node* pinned_node_pt=0;
 double min_dist_from_origin=DBL_MAX;
 
 // Pin the node that is at the centre in the domain
 unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(Inner_boundary0);
 for (unsigned inod=0;inod<num_int_nod;inod++)
  {
   // Get node point
   Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(Inner_boundary0,inod);
   
   // Find the node with the largest x coordinate
   if (fabs(nod_pt->x(0))>max_x_potentially_pinned_node)
    {
     max_x_potentially_pinned_node=fabs(nod_pt->x(0));
     pinned_rotation_node_pt=nod_pt;
    }
   
   // If the node is on the other internal boundary too
   double dist=sqrt(pow(nod_pt->x(0),2)+pow(nod_pt->x(1),2));
   if (dist<min_dist_from_origin)
    {
     pinned_node_pt=nod_pt;
     min_dist_from_origin=dist;
    }
  }
 
#ifdef USE_KS
 
 // Get relevant information from first element
 ELEMENT* first_el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0));
 
 // We're setting all dofs to zero
 double value=0.0;
 
 // The three displacement directions
 for (unsigned i_field=0;i_field<3;i_field++)
  {
   const unsigned first_nodal_type_index =
    first_el_pt->first_nodal_type_index_for_field(i_field);
   
   // Types: 0: u; 1: u_n; 2: u_t; 3: u_nn; 4: u_tn; 5: u_tt
   if (i_field<2)
    {
     // In plane: just the value
     unsigned k_type=0;
     pinned_node_pt->pin(first_nodal_type_index + k_type);
     pinned_node_pt->set_value(first_nodal_type_index + k_type, value);

     oomph_info << "Pinning (in-plane) " << first_nodal_type_index + k_type
                << " at "
                << pinned_node_pt->x(0) << " "
                << pinned_node_pt->x(1) << " "
                << std::endl;
    }
   else
    {
     // Out of plane: w, w_n, w_t
     for (unsigned k_type=0;k_type<3;k_type++)
      {
       pinned_node_pt->pin(first_nodal_type_index + k_type);
       pinned_node_pt->set_value(first_nodal_type_index + k_type, value);
       
       oomph_info << "Pinning (oo-plane) " << first_nodal_type_index + k_type
                  << " at "
                  << pinned_node_pt->x(0) << " "
                  << pinned_node_pt->x(1) << " "
                  << std::endl;
      }
    }
    
   // Pin y displacement at node at furthest x distance to suppress rotation about
   // the vertical axis
   if (i_field==1)
    {
     unsigned k_type=0;
     pinned_rotation_node_pt->pin(first_nodal_type_index + k_type);
     pinned_rotation_node_pt->set_value(first_nodal_type_index + k_type, value);
     
     oomph_info << "Pinning (z rot via y displ) " << first_nodal_type_index + k_type
                << " at "
                << pinned_rotation_node_pt->x(0) << " "
                << pinned_rotation_node_pt->x(1) << " "
                << std::endl;
    }
    
  }
 
#else
 
 // Constrain central node which is not rotated (though it doesn't
 // really matter if it was; we can either pin dw/dx and dw/dy or dw/dn
 // and dw/dt (relative to whatever directions the dof has been rotated
 // to)
 // - In-plane dofs are values 0 and 1
 // - Out of plane displacement is value 2;
 // - x and y (or n and t) derivatives of w are values 3 and 4.
 pinned_node_pt->pin(0);
 pinned_node_pt->set_value(0,0.0);
 pinned_node_pt->pin(1);
 pinned_node_pt->set_value(1,0.0);
 pinned_node_pt->pin(2);
 pinned_node_pt->set_value(2,0.0);
 pinned_node_pt->pin(3);
 pinned_node_pt->set_value(3,0.0);
 pinned_node_pt->pin(4);
 pinned_node_pt->set_value(4,0.0);
 
 oomph_info << "Pinning (FvK dofs 0,1,2,3,4) "
            << " at "
            << pinned_node_pt->x(0) << " "
            << pinned_node_pt->x(1) << " "
            << std::endl;
 
 // Pin y displacement at node at furthest x distance to suppress rotation about
 // the vertical axis
 pinned_rotation_node_pt->pin(1);
 
 oomph_info << "Pinning (FvK dofs 1) "
            << " at "
            << pinned_rotation_node_pt->x(0) << " "
            << pinned_rotation_node_pt->x(1) << " "
            << std::endl;
 
#endif
}



//==start_of_pin_for_balance_on_edge=========================================
/// Apply bcs so that sheet is balanced on edge in the middle
//============================================================================
template<class ELEMENT>
void UnstructuredC1PlateProblem<ELEMENT>::pin_for_balance_on_edge()
{
 
 // Pin all nodes along internal boundary "0"
 unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(Inner_boundary0);
 for (unsigned inod=0;inod<num_int_nod;inod++)
  {
   // Get node point
   Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(Inner_boundary0,inod);
   
#ifdef USE_KS
   
   // Get relevant information from first element
   ELEMENT* first_el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0));
   
   // We're setting all dofs to zero
   double value=0.0;
   
   // The three displacement directions
   for (unsigned i_field=0;i_field<3;i_field++)
    {
     const unsigned first_nodal_type_index =
      first_el_pt->first_nodal_type_index_for_field(i_field);
     
     // Types: 0: u; 1: u_n; 2: u_t; 3: u_nn; 4: u_tn; 5: u_tt
     //        0: u; 1: u_x; 2: u_y; 3: u_xx; 4: u_xy; 5: u_yy
     if (i_field<2)
      {
       // In plane: just the value
       unsigned k_type=0;
       nod_pt->pin(first_nodal_type_index + k_type);
       nod_pt->set_value(first_nodal_type_index + k_type, value);
      }
     else
      {
       // Out of plane: w, w_x, w_xx
       for (unsigned k_type=0;k_type<4;k_type++)
        {
         if (k_type!=2)
          {
           nod_pt->pin(first_nodal_type_index + k_type);
           nod_pt->set_value(first_nodal_type_index + k_type, value);
          }
        }
      }
     
     
  }
   
#else
   
   // - In-plane dofs are values 0 and 1
   // - Out of plane displacement is value 2;
   // - x and y (or t and n) derivatives of w are values 3 and 4.


   // oomph_info << "Pinning at: "
   //            << nod_pt->x(0) << " "
   //            << nod_pt->x(1) << " "
   //            << std::endl;

   
   // u
   nod_pt->pin(0);
   nod_pt->set_value(0,0.0);
   // v
   nod_pt->pin(1);
   nod_pt->set_value(1,0.0);

   // Only vertex nodes have w
   if (nod_pt->nvalue()>2)
    {
     // w
     nod_pt->pin(2);
     nod_pt->set_value(2,0.0);
     // w_x
     nod_pt->pin(3);
     nod_pt->set_value(3,0.0);
     
      // w_y (to suppress rotation about edge)
     nod_pt->pin(4);
     nod_pt->set_value(4,0.0);
    }
   
#endif
   
  }
}






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


#ifndef USE_KS
 
 // Full soln (apparently not implemented for KS; hierher add it)
 sprintf(filename,"%s/full_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file2.open(filename);
 unsigned nel=Bulk_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->full_output(some_file2,Parameters::Nplot);
  }
 some_file2.close();

 #endif
 
 // Output line visualiser solution 
 sprintf(filename,"%s/line_soln%i.dat",
         Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 LV_pt->output(some_file);
 some_file.close();



  // Increment the doc_info number
  Doc_info.number()++;

} // end of doc



// // hierher

// //========================================================================
// /// Validate interpolation of normal derivative along curved edge
// //========================================================================
// template<class ELEMENT>
// template<unsigned M>
// void UnstructuredC1PlateProblem<ELEMENT>::validate_dpsi_dn_along_edge(
//  const std::string& dir_name_for_output)
// {

//  ofstream some_file;
//  char filename[100];
 
// // Test & plot 'em
//  bool plot_em=true;
//  if (dir_name_for_output=="") plot_em=false;
//  if (plot_em)
//   {
//    sprintf(filename,"%s/test_curved_element.dat",
//            dir_name_for_output.c_str());
//    some_file.open(filename);
//   }
 
//  // Find a curved element on the outer boundary
//  // hierher (could actually do this for all of them)
//  // unsigned nb=Bulk_mesh_pt->nboundary();
//  // for (unsigned b=0;b<nb;b++)
//  unsigned b=Outer_boundary0;
//  {
//   const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);
//   for(unsigned e=0;e<nb_element;e++)
//    {
//     // Get pointer to bulk element adjacent to b
//     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
//      Bulk_mesh_pt->boundary_element_pt(b,e));
    
//     // Output the lot
//     unsigned nplot=30;
//     el_pt->full_output(some_file,nplot);

//     // Get basis functions
//     basis_w_foeppl_von_karman(const Vector<double>& s,
//                               Shape& psi_n,
//                               Shape& psi_i) const

// //     // All in CurvableBellElement:

 
// //     /// Access function for the Bernadou_element_basis_pt
// //     BernadouElementBasisBase* bernadou_element_basis_pt()
// //     {
// //       // [zdec] Should this throw an error if not upgraded or just return null
// //       // pt?
// //       return Bernadou_element_basis_pt;
// //     }
    
// //  /// Get the physical coordinate
// //     template<unsigned BOUNDARY_ORDER>
// //     void BernadouElementBasis<BOUNDARY_ORDER>::coordinate_x(
// //       const Vector<double>& s, Vector<double>& fk) const
// //     {
// //       Vector<double> s_basic(s);
// //       permute_shape(s_basic);
// //       f_k(s_basic, fk);
// //     }


    
// //     /// Get the Bell/Bernadou basis for the unknowns
// //     virtual void c1_basis(const Vector<double>& s,
// //                           Shape& nodal_basis,
// //                           Shape& bubble_basis) const
// //     {
// //       if (element_is_curved())
// //       {
// //         Bernadou_element_basis_pt->shape(s, nodal_basis, bubble_basis);
// //       }


// //       // hierher should really drop down into the constituent functinos!
      
// //   //======================================================================
// //   /// Out-of-plane basis functions at local coordinate s
// //   //======================================================================
// //   template<unsigned NNODE_1D>
// //   void FoepplVonKarmanC1CurvableBellElement<
// //     NNODE_1D>::basis_w_foeppl_von_karman(const Vector<double>& s,
// //                                          Shape& psi_n,
// //                                          Shape& psi_i) const
// //   {
    
// //    // hierher Aidan: Kill this commented out bit?
// // //     throw OomphLibError("This still needs testing for curved elements.",
// // //                         "void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::shape_and_test_foeppl_von_karman(...)",
// // //                         OOMPH_EXCEPTION_LOCATION);
    
// //     this->c1_basis(s, psi_n, psi_i);
    
// //     // Rotate the degrees of freedom
// //     rotate_shape(psi_n);
// //   }


  
//     break; // hierher
//    }
//   //break; // hierher
//  }

//  if (plot_em)
//   {
//    some_file.close();
//   }
//  exit(0);
// }
 



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

   // Get map to curvline boundaries of mesh
   std::map<unsigned, C1CurviLine*> c1_curviline_boundary_pt =
    Bulk_mesh_pt->c1_curviline_boundary_pt();
   C1CurviLine* curviline_pt=c1_curviline_boundary_pt[b];
   
    
   // Find the dimension of the element [zdec] will this ever not be 2?
   const unsigned dim = el_pt->dim(); //2; // hierher should come from here 
    
   // The number of first derivatives is the dimension of the element
   const unsigned n_deriv = dim;

   // The number of second derivatives is the triangle number of the dimension
   const unsigned n_2deriv = dim * (dim + 1) / 2;
    
   // Find out how many nodes there are for w
   const unsigned n_w_node = el_pt->nw_node(); // 3; // hierher el_pt->nw_node();

   // Get the vector of nodes used for each field
   const Vector<unsigned> w_nodes = el_pt->get_w_node_indices(); // {0,1,2}; // hierher el_pt->get_w_node_indices();

   // Find out how many basis types there are at each node
   const unsigned n_w_nodal_type = el_pt->nw_type_at_each_node(); // 6; // hierher el_pt->nw_type_at_each_node();

   // Find out how many basis types there are internally
   unsigned n_w_internal_type =  el_pt->nw_type_internal(); // 3; // // hierher el_pt->nw_type_internal();
   //if (M==5) n_w_internal_type=10;
    
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
         double J =
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
                                       << dpsi_n_wdxi(j,k,0) << " " // 8
                                       << dpsi_n_wdxi(j,k,1) << " " // 9 
                                       << 1 - 3*s*s + 2*s*s*s << " " // 10
                                       << s - 2*s*s + s*s*s << " " // 11
                                       << 3*s*s - 2*s*s*s << " " // 12
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
              *(internal_file_pt[count]) << interp_x[0] << " " // 1
                                         << interp_x[1] << " " // 2
                                         << r_from_boundary[0] << " " // 3 
                                         << r_from_boundary[1] << " " // 4
                                         <<  drdzeta[1]/norm << " " // 5 
                                         << -drdzeta[0]/norm << " " // 6
                                         << s << " "  // 7
                                         << dpsi_i_wdxi(k_type,0) << " " // 8
                                         << dpsi_i_wdxi(k_type,1) << " " // 9 
                                         << 1 - 3*s*s + 2*s*s*s << " " // 10
                                         << s - 2*s*s + s*s*s << " " // 11
                                         << 3*s*s - 2*s*s*s << " " // 12
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
     exit(0);
     // oomph_info << "oomph-convert -z test_basic_basis*dat" << std::endl;
     // oomph_info << "makePvd test_basic_basis test_basic_basis.pvd" << std::endl;
     // oomph_info << "oomph-convert -p2 test_points.dat " << std::endl;
     // oomph_info << "paraview --state test_basic_basis.pvsm " << std::endl;
     // oomph_info << std::endl;
     
    }
  }
 }
}
   







//========================================================================
/// Validate mapping from monomials to 36 [66] basic dofs
//========================================================================
template<class ELEMENT>
template<unsigned M>
void UnstructuredC1PlateProblem<ELEMENT>::
validate_monomials_to_basic_basis_functions(const std::string&
                                            dir_name_for_output)
{



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
double tol=1.0e-10;
bool test_passed=true;
 for (unsigned i=0;i<n_basic;i++)
  {
   unsigned count_one_in_row=0;
   unsigned count_zero_in_row=0;
   unsigned count_one_in_col=0;
   unsigned count_zero_in_col=0;
   for (unsigned j=0;j<n_basic;j++)
    {
     if (std::abs(test_matrix(i,j)    )<tol) count_zero_in_row++;
     if (std::abs(test_matrix(i,j)-1.0)<tol) count_one_in_row++;
     if (std::abs(test_matrix(j,i)    )<tol) count_zero_in_col++;
     if (std::abs(test_matrix(j,i)-1.0)<tol) count_one_in_col++;
    }
   
   if ((count_one_in_row!=1)||
       (count_one_in_col!=1)||
       (count_zero_in_row!=(n_basic-1))||
       (count_zero_in_col!=(n_basic-1)))
    {
     test_passed=false;
     oomph_info << "Test failed: row/col test: "
                << i << ": "
                <<  count_one_in_row << " "
                <<  count_zero_in_row << " "
                <<  count_one_in_col << " "
                <<  count_zero_in_col << " "
                << std::endl;
     break;
    }
  }

if (test_passed)
 {
  oomph_info << "Test of basic basis functions passed!"
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
  
  exit(0);
 
 
  

}





 

//========================================================================
/// Namespace for function that calls doc_solution() during the damped
/// solves
//========================================================================
namespace DocProgressOfDampedSolutions
{

 /// Pointer to the problem class (to get access the doc solution function
 #ifdef USE_KS
 
  UnstructuredC1PlateProblem<KoiterSteigmannC1CurvableBellElement>*
  Problem_pt=0;

#else

 UnstructuredC1PlateProblem<FoepplVonKarmanC1CurvableBellElement<4>>*
   Problem_pt=0;

#endif

 /// Function to call doc_solution during damped solves
 void doc_solution_during_damped_solve(const unsigned& i_step)
 {
  oomph_info << "Docing solution for damped solve step "
             << i_step << std::endl;

  // bumps up counter by itself.
  bool steady=false;
  Problem_pt->doc_solution(steady);
  
 }

} // end of namespace



//=======start_of_main========================================
///Driver code 
//============================================================
int main(int argc, char** argv)
{
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);


  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  // Define possible command line arguments and parse the ones that
  // were actually specified

  // Clamped boundary conditions?
  CommandLineArgs::specify_command_line_flag("--nplot",&Parameters::Nplot);

  // T-shaped internal boundary
  CommandLineArgs::specify_command_line_flag("--use_t_shape_internal_boundaries");

  // Use polyline for internal boundaries
  CommandLineArgs::specify_command_line_flag("--use_polyline_for_internal_boundaries");

  // Clamped boundary conditions?
  CommandLineArgs::specify_command_line_flag("--use_clamped_bc");

  // Pinned boundary conditions?
  CommandLineArgs::specify_command_line_flag("--use_pinned_bc");

  // Balance on edge boundary conditions?
  CommandLineArgs::specify_command_line_flag("--use_balance_on_edge_bc");
  
  // Rotate coords?
  CommandLineArgs::specify_command_line_flag
   ("--do_not_rotate_coords_on_curved_boundaries");

  // Element area
  CommandLineArgs::specify_command_line_flag("--el_area",
                                             &Parameters::Element_area);
  
  // // Square outer boundary (straight curvilines)
  // CommandLineArgs::specify_command_line_flag
  //  ("--outer_boundary_straight_curved");
  
  // // Square outer boundary (polygonal)
  // CommandLineArgs::specify_command_line_flag
  //  ("--outer_boundary_straight_poly");

  // hierher check that not both are specified

  
  // Test drive damped solve
  CommandLineArgs::specify_command_line_flag
   ("--test_damped_solve");
  
  
  // Parse command line
  CommandLineArgs::parse_and_assign();

  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();

  if (CommandLineArgs::command_line_flag_has_been_set
      ("--do_not_rotate_coords_on_curved_boundaries"))
   {
    Parameters::Rotate_coordinates_on_all_curvilinear_boundaries=false;
   }
  
  
  // Check consistency
  if (CommandLineArgs::command_line_flag_has_been_set("--use_clamped_bc"))
  {
    Parameters::Problem_case = Parameters::Clamped_validation;
    if (!Parameters::Rotate_coordinates_on_all_curvilinear_boundaries)
     {
      oomph_info << "clamped bcs require rotated dofs on boundary" << std::endl;
      abort();
     }
  }
  
  // Check consistency
  if (CommandLineArgs::command_line_flag_has_been_set("--use_pinned_bc"))
  {
    Parameters::Problem_case = Parameters::Pinned_validation;
    if (!Parameters::Rotate_coordinates_on_all_curvilinear_boundaries)
     {
      oomph_info << "pinned bcs require rotated dofs on boundary" << std::endl;
      abort();
     }
  }

  /// Balance on edge
  if (CommandLineArgs::command_line_flag_has_been_set("--use_balance_on_edge_bc"))
   {
    Parameters::Problem_case = Parameters::Balance_on_edge;
  }
  

#ifdef USE_KS
  
  // Create the problem, using FvK elements derived from TElement<2,4>
  // elements (with 4 nodes per element edge and 10 nodes overall).
  UnstructuredC1PlateProblem<KoiterSteigmannC1CurvableBellElement>
    problem(Parameters::Element_area);

#else

  // Build problem 
  UnstructuredC1PlateProblem<FoepplVonKarmanC1CurvableBellElement<4>> problem(
    Parameters::Element_area);

#endif

  // problem.validate_monomials_to_basic_basis_functions<5>();
  // problem.validate_monomials_to_basic_basis_functions<3>();


  std::string dir_name="RESLT";
  problem.validate_curved_bell_and_bubble_basis_functions(dir_name);

  
   //problem.validate_dpsi_dn_along_edge<3>(".");
  exit(0);
  
  // Pass problem pointer to namespace for docing damped solves
  DocProgressOfDampedSolutions::Problem_pt=&problem;


  // Update/set non-dim parameters
  Parameters::update_nondim_parameters();
   
  // Tweak Newton solver parameters
  problem.max_residuals() = 1.0e3;

  // Document the initial state
  problem.doc_solution();

  // Set pressure increment
  unsigned n_step = 100;
  double p_inc = Parameters::P_max/double(n_step); // 1.0e-2;

  // Initialise actual pressure
  Parameters::P_mag = 0.0;

  oomph_info << "Doing nstep = " << n_step << " pressure increments of "
             << p_inc << std::endl;

  // exit(0);

  
  if (CommandLineArgs::command_line_flag_has_been_set("--test_damped_solve"))
   {
    // // 0.1 and 100 steps gives nice animation
    // p_inc=1.0;
    // n_step=10;
    Parameters::P_cos=1.0;
   }

  
  // Overwrite for "Balance on Edge" case
  if (Parameters::Problem_case == Parameters::Balance_on_edge)
   {
    p_inc = 1.0; 
    n_step = 3; 
   }
  


  for( unsigned i = 0; i < n_step; i++ )
  {
   // Bump
   Parameters::P_mag += p_inc;

   if (!CommandLineArgs::command_line_flag_has_been_set("--test_damped_solve"))
    {
     // Solve the system
     problem.newton_solve();
    }
   else
    {
     // initial value for timestep
     double dt=1.0;
     
     // tolerance for adaptive timestepping; somewhat random
     // hierher Aidan: any recommendations?
     double epsilon=0.01; // ten times smaller shows timestepping nicely. 1.0e-3;

     // Damped solve
     double suggested_next_dt=
      problem.damped_solve(dt,epsilon,
                           &DocProgressOfDampedSolutions::doc_solution_during_damped_solve);

     // Can (but don't have to) to use this for next solve
     oomph_info << "Suggested next dt = " << suggested_next_dt << std::endl;
    }

   
   // Document the current solution
   problem.doc_solution();
  }


} // End of main
