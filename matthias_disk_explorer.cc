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



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////



//========================================================================
/// Namespace for problem parameters
//========================================================================
namespace Parameters
{

 /// Enumeration of cases
 enum
 {
  Clamped_validation,
  Pinned_validation,
  Free_edges
 };

 /// Which case are we doing
 unsigned Problem_case = Free_edges;
 
 /// Ellipse half x-axis
 double A = 1.0;
 
 /// Ellipse half y-axis
 double B = 1.0;
 
 /// Poisson ratio
 double Nu = 0.5;
 
 /// Nondimensional thickness of plate
 double Thickness = 0.01;

 #ifdef USE_KS

 // hierher update these to make them consistent with fvk
 
 /// Membrane coupling coefficient 
 double Eta_u = 1.0;

  /// hierher what is this?
 double Eta_sigma = 1.0; 

#else
 
 /// Membrane coupling coefficient (this should really be computed
 /// as a dependent parameter...)
 double Eta = 12.0 * (1.0 - Nu * Nu) / (Thickness * Thickness);

 #endif


 /// Pressure magnitude
  double P_mag = 0.0;

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
  }

#else
 
 /// Pressure depending on the position (x,y)
  void get_pressure(const Vector<double>& x, double& pressure)
  {
   pressure = P_mag;
  }

  /// In plane forcing (shear stress) depending on the position (x,y)
  void get_in_plane_traction(const Vector<double>& x, Vector<double>& tau)
  {
   // Zero shear stress
   tau[0]=0.0;
   tau[1]=0.0;
  }

 #endif




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
      ->validate_and_pin_redundant_constraints();
    }
   
   // Reassign the equation numbers
   oomph_info << "Reassiging equation numbers after changing BCs. "
              << " ndof = " << assign_eqn_numbers() << std::endl;
   
  } // End make_linear
 
#endif
 
 
 /// Doc the solution
 void doc_solution(const std::string& comment="");

 
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
#ifdef PARANOID
            // We should only have one coordinate on this boundary
            unsigned nzeta=nod_pt->ncoordinates_on_boundary(b);
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
 // hierher move this into mesh and describe output and prefix
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
     std::cout << "Boundary: " << b << " is adjacent to " << nel << " elements"
               << std::endl;
     
     // Loop over elements on given boundary
     for (unsigned e = 0; e < nel; e++)
      {
       FiniteElement* fe_pt = mesh_pt->boundary_element_pt(b,e);
       unsigned nnod_1d=fe_pt->nnode_1d();
       fe_pt->output(bulk_file,nnod_1d);
       unsigned face_index=mesh_pt->face_index_at_boundary(b,e);
       std::cout << "Boundary element:" << fe_pt
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
         std::cout << "Never get here: " << nnod_1d << std::endl;
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

 
  /// Pin all displacements and rotation at the centre
  void pin_all_displacements_and_rotation_at_centre_node();

  /// Trace file to document norm of solution
  ofstream Trace_file;

  /// Pointer to "bulk" mesh
  TriangleMesh<ELEMENT>* Bulk_mesh_pt;

 
  /// Enumeration to keep track of boundary ids
  enum
  {
    Outer_boundary0 = 0,
    Outer_boundary1 = 1,
    Inner_boundary0 = 2
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
UnstructuredC1PlateProblem<ELEMENT>::UnstructuredC1PlateProblem(const double&
                                                                element_area)
 : Element_area(element_area)
{

 // Build the mesh
 //================
 
 Vector<double> zeta(1);
 Vector<double> posn(2);
 
 //Outer boundary
 //--------------
 
 double A = Parameters::A;
 double B = Parameters::B;
 Ellipse* outer_boundary_ellipse_pt = new Ellipse(A, B);


 // hierher MH break up into N randomly enumerated bits; smooth meets
 // smooth; straight lines, etc.
 
 // Storage for outer boundaries (for triangle)
 Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(2);

 //First bit
 double zeta_start = 0.0;
 double zeta_end = MathematicalConstants::Pi;
 unsigned nsegment = (unsigned)(MathematicalConstants::Pi/sqrt(Element_area));
 outer_curvilinear_boundary_pt[0] = 
  new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
                            zeta_end, nsegment, Outer_boundary0);
 
 
 //Second bit
 zeta_start = MathematicalConstants::Pi;
 zeta_end = 2.0*MathematicalConstants::Pi;
 nsegment = (int)(MathematicalConstants::Pi/sqrt(Element_area));
 outer_curvilinear_boundary_pt[1] =
  new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
                            zeta_end, nsegment, Outer_boundary1);
 
 // Combine
  TriangleMeshClosedCurve* outer_boundary_pt =
  new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);
 
 // Internal open boundaries
 //-------------------------
 // Total number of open curves in the domain
 unsigned n_open_curves = 1;
 
 // We want internal open curves
 Vector<TriangleMeshOpenCurve *> inner_open_boundaries_pt(n_open_curves);


 // hierher MH: optinoally re-introduce cross in the middle to check
 // element splitting
 
 
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
 
 TriangleMeshPolyLine* boundary2_pt =
  new TriangleMeshPolyLine(vertices, boundary_id);
  
  // Each internal open curve is defined by a vector of
  // TriangleMeshCurveSections
  Vector<TriangleMeshCurveSection *> internal_curve_section1_pt(1);
  internal_curve_section1_pt[0] = boundary2_pt;
    
  // The open curve that defines this boundary
  inner_open_boundaries_pt[0] =
   new TriangleMeshOpenCurve(internal_curve_section1_pt);

  //Create mesh parameters object
  TriangleMeshParameters mesh_parameters(outer_boundary_pt);

  // Element area
  mesh_parameters.element_area() = Element_area;

  // Specify the internal open boundaries
  mesh_parameters.internal_open_curves_pt() = inner_open_boundaries_pt;

  // Build an assign bulk mesh
  Bulk_mesh_pt=new TriangleMesh<ELEMENT>(mesh_parameters);


  
  // Now upgrade to (potentially) curved C1 boundaries 
  //==================================================
  {

  
   // Create the mesh for the Lagrange multiplier elements that enforce
   // continuity of our smooth solution across different parts of the
   // mesh boundary (only really needed when there are kinks)
   Constraint_mesh_pt = new Mesh();

   // hierher inside this helper function issue warning if any of
   // the boundaries are not curvilines
   
   // hierher explain
   C1Helper::upgrade_triangle_mesh_for_c1_plate_bending<ELEMENT>(
    Bulk_mesh_pt,
    Constraint_mesh_pt);

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
    
    // hierher: pressure --> traction in src
    el_pt->pressure_fct_pt() = &Parameters::get_traction;

    // hierher why do we need thickness and (two!) etas?
    el_pt->thickness_pt() = &Parameters::Thickness;
    el_pt->nu_pt() = &Parameters::Nu;
    el_pt->eta_u_pt() = &Parameters::Eta_u;
    el_pt->eta_sigma_pt() = &Parameters::Eta_sigma;

    
    // hierher need an example that uses this!
    // el_pt->mu_pt() = &Parameters::Mu;

#else
    
    el_pt->pressure_fct_pt() = &Parameters::get_pressure;
    el_pt->nu_pt() = &Parameters::Nu;
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
   
    // Set the boundary conditions on the two outer boundaries
   unsigned nbound = 2;
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
   
    // Set the boundary conditions on the two outer boundaries
   unsigned nbound = 2;
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
  // boundary conditions
  unsigned n_el = Constraint_mesh_pt->nelement();
  for(unsigned i_el = 0; i_el < n_el; i_el++)
   {
    // hierher rename in src
    dynamic_cast<DuplicateNodeConstraintElement*>
     (Constraint_mesh_pt->element_pt(i_el))
     ->validate_and_pin_redundant_constraints();
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
 unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(2);
 for (unsigned inod=0;inod<num_int_nod;inod++)
  {
   // Get node point
   Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(2,inod);
   
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
    }
   else
    {
     // Out of plane: w, w_n, w_t
     for (unsigned k_type=0;k_type<3;k_type++)
      {
       pinned_node_pt->pin(first_nodal_type_index + k_type);
       pinned_node_pt->set_value(first_nodal_type_index + k_type, value);
      }
    }
    
   // Pin y displacement at node at furthest x distance to suppress rotation about
   // the vertical axis
   if (i_field==1)
    {
     unsigned k_type=0;
     pinned_rotation_node_pt->pin(first_nodal_type_index + k_type);
     pinned_rotation_node_pt->set_value(first_nodal_type_index + k_type, value);
    }
    
  }
 
#else
 
 // Constrain central node which is not rotated (though it doesn't
 // really matter if it was; we can either pin dw/dx and dw/dy or dw/dn
 // and dw/dt (relative to whatever directions the dof has been rotated
 // to)
 // - In-plane dofs are values 0 and 1
 // - Out of plane displacement is value 2;
 // - x and y (or t and n) derivatives of w are values 3 and 4.
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
 
 // Pin y displacement at node at furthest x distance to suppress rotation about
 // the vertical axis
 pinned_rotation_node_pt->pin(1);
 
#endif
}





//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredC1PlateProblem<ELEMENT>::doc_solution(
 const std::string& comment)
{
  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts = 5; // 50;

  sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file,npts);
  some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
  << comment << "\"\n";
  some_file.close();

  // Increment the doc_info number
  Doc_info.number()++;

} // end of doc



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
  CommandLineArgs::specify_command_line_flag("--use_clamped_bc");

  // Clamped boundary conditions?
  CommandLineArgs::specify_command_line_flag("--use_pinned_bc");

  
  // Rotate dofs?
  CommandLineArgs::specify_command_line_flag("--rotate_dofs_on_boundary");

  // Parse command line
  CommandLineArgs::parse_and_assign();

  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();

  // Check consistency
  if (CommandLineArgs::command_line_flag_has_been_set("--use_clamped_bc"))
  {
    Parameters::Problem_case = Parameters::Clamped_validation;
    if (!CommandLineArgs::command_line_flag_has_been_set("--rotate_dofs_on_boundary"))
     {
      oomph_info << "clamped bcs require rotated dofs on boundary" << std::endl;
      abort();
     }
  }
  
  // Check consistency
  if (CommandLineArgs::command_line_flag_has_been_set("--use_pinned_bc"))
  {
    Parameters::Problem_case = Parameters::Pinned_validation;
    if (!CommandLineArgs::command_line_flag_has_been_set("--rotate_dofs_on_boundary"))
     {
      oomph_info << "pinned bcs require rotated dofs on boundary" << std::endl;
      abort();
     }
  }
 

#ifdef USE_KS
  
  // Create the problem, using FvK elements derived from TElement<2,4>
  // elements (with 4 nodes per element edge and 10 nodes overall).
  UnstructuredC1PlateProblem<KoiterSteigmannC1CurvableBellElement>
    problem(Parameters::Element_area);

#else

  // Build problem // hierher what's the 4 for? What else can I do
  UnstructuredC1PlateProblem<FoepplVonKarmanC1CurvableBellElement<4>> problem(
    Parameters::Element_area);

#endif

  
  // Tweak Newton solver parameters
  problem.max_residuals() = 1.0e3;
  problem.max_newton_iterations() = 100;
  problem.newton_solver_tolerance() = 1.0e-11;

  // Document the initial state
  problem.doc_solution();

  // Set the Poisson ratio
  Parameters::Nu = 0.5;
  
  // Do we want to solve the linear problem?
  // problem.make_linear();

  // Set pressure and incrementation
  Parameters::P_mag = 0.0;
  double p_inc = 1.0e-3;
  unsigned n_step = 10;
  for( unsigned i = 0; i < n_step; i++ )
  {
   // Bump
   Parameters::P_mag += p_inc;
   
   // Solve the system
   problem.newton_solve();
   
   // Document the current solution
   problem.doc_solution();
  }


} // End of main
