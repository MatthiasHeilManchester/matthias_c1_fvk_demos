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

// The equations
#include "c1_foeppl_von_karman.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;

//==============================================================================
/// Namespace to deal update triangle meshes to deal with C1 elements
// hierher this will move into C1_helper.h in src/generic
//==============================================================================
namespace C1Helper
{
 
 
//==============================================================================
// hierher update
/// Duplicate nodes at corners in order to properly apply boundary
/// conditions from each edge. Also adds (8) Lagrange multiplier dofs to the
/// problem in order to constrain continuous interpolation here across its (8)
/// vertex dofs. (Note "corner" here refers to the meeting point of any two
/// sub-boundaries in the closed external boundary)
//==============================================================================
 void duplicate_corner_nodes(Mesh* bulk_mesh_pt, 
                             std::map<unsigned,C1CurviLine*> c1_curviline_pt,
                             Mesh* constraint_mesh_pt)
 {

  // hierher check if mesh is distributed!
  
  // Collection of nodes that occupy two boundaries together with the boundary IDs
  // (ordered: first < second)
  std::map<Node*,std::pair<unsigned,unsigned>> boundaries_of_boundary_node_pt;

  // Loop over the curvilinear parts of the outer boundary
  for (const auto& [i_bound, para] : c1_curviline_pt)
  {
   unsigned n_b_node = bulk_mesh_pt->nboundary_node(i_bound);
   for(unsigned i_b_node = 0; i_b_node < n_b_node; i_b_node++)
    {
     // Store the node we are checking
     Node* node_pt = bulk_mesh_pt->boundary_node_pt(i_bound,i_b_node);

     // Pointer to set that contains the boundaries we're on
     std::set<unsigned>* boundaries_pt=0;
     node_pt->get_boundaries_pt(boundaries_pt);
     if (boundaries_pt!=0)
      {
       if (boundaries_pt->size()==2)
        {
         unsigned b_min=UINT_MAX;
         unsigned b_max=0;
         for (unsigned b : (*boundaries_pt))
          {
           oomph_info << "Node " << node_pt << " is on boundary " << b << std::endl;
           if (b<b_min) b_min=b;
           if (b>b_max) b_max=b;           
          }
         // Ordered!
         boundaries_of_boundary_node_pt[node_pt].first=b_min;
         boundaries_of_boundary_node_pt[node_pt].second=b_max;
        }
      }
    }
  }
  
  // Here are the nodes that need to be duplicated. We duplicate them
  // on the lower of its two boundaries (this is stored first)
  for (auto a : boundaries_of_boundary_node_pt)
   {
    Node* node_to_be_duplicated_pt=a.first;
    unsigned boundary_on_which_node_is_duplicated=a.second.first;
    unsigned boundary_on_which_node_is_left=a.second.second;
    
    oomph_info << "Node " <<  node_to_be_duplicated_pt
               << " is duplicated on boundary "
               << boundary_on_which_node_is_duplicated << " and kept on boundary"
               << boundary_on_which_node_is_left
               << std::endl;
   
    // Find the boundary element that contains the node to be duplicated on
    // the boundary where the node is to be duplicated
    FiniteElement* el_where_node_is_to_be_duplicated_pt=0;
    unsigned n_b_el = bulk_mesh_pt->nboundary_element(boundary_on_which_node_is_duplicated);
    for (unsigned i_b_el = 0; i_b_el < n_b_el; i_b_el++)
    {
      // Get the element pointer
      FiniteElement* el_pt = bulk_mesh_pt->boundary_element_pt
       (boundary_on_which_node_is_duplicated, i_b_el);
      // If the corner node pt is in the element we have found the right
      // element
      if (el_pt->get_node_number(node_to_be_duplicated_pt) != -1)
       {
        el_where_node_is_to_be_duplicated_pt = el_pt;
        break;
      }
    }

    oomph_info << "Boundary element that contains that node: "
               << el_where_node_is_to_be_duplicated_pt << std::endl;
    
    // Now we need to create a new node and substitute the element's
    // old corner node for this new one
    Node* new_node_pt = el_where_node_is_to_be_duplicated_pt->construct_boundary_node(
     el_where_node_is_to_be_duplicated_pt->get_node_number(node_to_be_duplicated_pt));
    
    // Copy the position and other info from the old node into the new node
    new_node_pt->x(0)=node_to_be_duplicated_pt->x(0);
    new_node_pt->x(1)=node_to_be_duplicated_pt->x(1);

    // Then we add this node to the mesh
    bulk_mesh_pt->add_node_pt(new_node_pt);

    // Then replace the old node for the new one on the boundary
    bulk_mesh_pt->remove_boundary_node(boundary_on_which_node_is_duplicated,node_to_be_duplicated_pt);
    bulk_mesh_pt->   add_boundary_node(boundary_on_which_node_is_duplicated,new_node_pt);

    // The final job is to constrain this duplication using the specialised
    // Lagrange multiplier elements which enforce equality of displacement and
    // its derivatives either side of this corner.
    C1CurviLine* left_parametrisation_pt  = c1_curviline_pt[boundary_on_which_node_is_left];
    C1CurviLine* right_parametrisation_pt = c1_curviline_pt[boundary_on_which_node_is_duplicated];

    // Get the coordinates on each node on their respective boundaries
    Vector<double> left_boundary_coordinate =
     {left_parametrisation_pt->get_zeta(node_to_be_duplicated_pt->position())};
    Vector<double> right_boundary_coordinate =
     {right_parametrisation_pt->get_zeta(new_node_pt->position())};

    // Create the constraining element
    DuplicateNodeConstraintElement* constraint_element_pt =
     new DuplicateNodeConstraintElement(node_to_be_duplicated_pt,
                                        new_node_pt,
                                        left_parametrisation_pt,
                                        right_parametrisation_pt,
                                        left_boundary_coordinate,
                                        right_boundary_coordinate);

    // Add the constraining element to the mesh
    constraint_mesh_pt->add_element_pt(constraint_element_pt);
   }   
 }



//==============================================================================
/// A function that upgrades straight sided elements to be curved. This involves
/// Setting up the parametric boundary, F(s) and the first derivative F'(s)
/// We also need to set the edge number of the upgraded element and the positions
/// of the nodes j and k (defined below) and set which edge (k) is to be exterior
///            @ k               
///           /(                 
///          /. \                
///         /._._)               
///      i @     @ j             
/// For RESTING or FREE boundaries we need to have a C2 CONTINUOUS boundary
/// representation. That is we need to have a continuous 2nd derivative defined
/// too. This is well discussed in by [Zenisek 1981] (Aplikace matematiky ,
/// Vol. 26 (1981), No. 2, 121--141). This results in the necessity for F''(s)
/// as well.
//=============================================================================
 void upgrade_edge_elements_to_curved_boundaries(
  Mesh* bulk_mesh_pt, 
  std::map<unsigned,C1CurviLine*> c1_curviline_pt) 
 {
  
  // Loop over the curvilinear parts of the outer boundary
  for (const auto& [ibound, c1_curve_pt] : c1_curviline_pt)
   {
    
    // Loop over the bulk elements adjacent to boundary ibound
    const unsigned n_els=bulk_mesh_pt->nboundary_element(ibound);
    for(unsigned e=0; e<n_els; e++)
     {
      // Get pointer to bulk element adjacent to b
      FiniteElement* bulk_el_pt =
       bulk_mesh_pt->boundary_element_pt(ibound,e);
      
      // hierher Aidan what is that? why "My"?
      // Initialise enum for the curved edge
      MyC1CurvedElements::Edge edge(MyC1CurvedElements::none);
      
      // Loop over all (three) vertex nodes of the element and
      // identify single node that is interior (i.e. not on any
      // of the outer boundaries
      unsigned index_of_interior_node = 3;
      unsigned nnode_not_on_any_outer_boundary = 0;
      const unsigned nnode = 3;
      Vector<Vector<double> > xn(nnode,Vector<double>(2,0.0));
      for(unsigned n=0;n<nnode;++n)
       {
        Node* nod_pt = bulk_el_pt->node_pt(n);
        xn[n][0]=nod_pt->x(0);
        xn[n][1]=nod_pt->x(1);
        
        // Check if it is on any of the outer boundaries
        bool node_is_on_some_outer_boundary=false;
        for (const auto& [b, dummy_c1_curve_pt] : c1_curviline_pt)
         {
          if (nod_pt->is_on_boundary(b))
           {
            node_is_on_some_outer_boundary=true;
            break;
           }
         }
        if (!node_is_on_some_outer_boundary)
         {
          index_of_interior_node = n;
          nnode_not_on_any_outer_boundary++;
         }
       }// end record boundary nodes
      
      // hierher shouldn't these be called zeta (everywhere; sigh)
      // boundary coordinate at the next (cyclic) node after interior
      const double s_ubar =
       c1_curve_pt->get_zeta(xn[(index_of_interior_node+1) % 3]);
      
      // boundary coordinate at the previous (cyclic) node before interior
      const double s_obar =
       c1_curve_pt->get_zeta(xn[(index_of_interior_node+2) % 3]);
      
      // Assign edge case
      edge = static_cast<MyC1CurvedElements::Edge>(index_of_interior_node);
      
#ifdef PARANOID
      // Check nnode_on_neither_boundary
      if (nnode_not_on_any_outer_boundary == 0)
       {
        throw OomphLibError(
         "No interior nodes. One node per CurvedElement must be interior.",
         OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
       }
      else if (nnode_not_on_any_outer_boundary> 1)
       {
        throw OomphLibError(
         "Multiple interior nodes. Only one node per CurvedElement can be interior.",
         OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
       }
      
      // Check for inverted elements
      if (s_ubar>s_obar)
       {
        throw OomphLibError(
         "Decreasing parametric coordinate. Parametric coordinate must increase as the edge is traversed anti-clockwise.",
         OOMPH_CURRENT_FUNCTION,
         OOMPH_EXCEPTION_LOCATION);
       } // end checks
#endif
      
      // Upgrade it
      TemplateFreeCurvableBellElement* curv_el_pt=
       dynamic_cast<TemplateFreeCurvableBellElement*>(bulk_el_pt);
#ifdef PARANOID
      if (curv_el_pt==0)
       {
        throw OomphLibError(
         "Cast to TemplateFreeCurvableBellElement failed",
         OOMPH_CURRENT_FUNCTION,
         OOMPH_EXCEPTION_LOCATION);
       }
#endif

      // hierher shouldn't we hard code this to only allow specific options
      // this can't be any number, right?
      unsigned boundary_order=5;
      curv_el_pt->upgrade_element_to_curved(edge, s_ubar, s_obar,
                                            c1_curve_pt,
                                            boundary_order);
     }

   } // end of loop over outer boundaries
  
 } // end_upgrade_elements


 
//======================================================================
/// Function to set up rotated nodes on the boundary: necessary if we want to set
/// up physical boundary conditions on a curved boundary with Hermite type dofs.
/// For example if we know w(n,t) = f(t) (where n and t are the
/// normal and tangent to a boundary) we ALSO know dw/dt and d2w/dt2.
/// NB no rotation is needed if the edges are completely free!
//======================================================================
 void rotate_edge_degrees_of_freedom(
  Mesh* bulk_mesh_pt, 
  std::map<unsigned,C1CurviLine*> c1_curviline_pt) 
{
 
 // Loop over the bulk elements: Yes, really because we also need to deal with those that only
 // have a single node on the boundary!
 unsigned n_element = bulk_mesh_pt-> nelement();
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element 
   FiniteElement* el_pt = bulk_mesh_pt->finite_element_pt(e);
   
   // Loop over the curvilinear parts of the outer boundary
   for (const auto& [b, c1_curve_pt] : c1_curviline_pt)
    {
     // local node numbers of nodes on external boundaries
     Vector<unsigned> boundary_node;
     
     // Boundary coordinates of nodes on the external boundaries
     Vector<double> boundary_coordinate_of_node;
     
     // Loop over vertex nodes (they come first)
     const unsigned nnode=3;
     for (unsigned n=0; n<nnode;++n)
      {
       // If on external boundary b
       if (el_pt->node_pt(n)->is_on_boundary(b))
        {
         boundary_node.push_back(n);
         double coord = c1_curve_pt->get_zeta(el_pt->node_pt(n)->position());
         boundary_coordinate_of_node.push_back(coord);
        }
      }
     
     // If the element has nodes on the boundary, rotate the Hermite dofs
     if(!boundary_node.empty())
      {
       // Rotate the nodes by passing the index of the nodes and the
       // normal / tangent vectors to the element
       
       // Upgrade it
       TemplateFreeCurvableBellElement* curv_el_pt=
        dynamic_cast<TemplateFreeCurvableBellElement*>(el_pt);
#ifdef PARANOID
       if (curv_el_pt==0)
        {
         throw OomphLibError(
          "Cast to TemplateFreeCurvableBellElement failed",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
        }
#endif
       
       curv_el_pt->
        rotated_boundary_helper_pt()->
        set_nodal_boundary_parametrisation(boundary_node,
                                           boundary_coordinate_of_node,
                                           c1_curve_pt);
      }
    }
  }
 
} // end rotate_edge_degrees_of_freedom



} // end namespace


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
 
 /// Membrane coupling coefficient (this should really be computed
 /// as a dependent parameter...)
 double Eta = 12.0 * (1.0 - Nu * Nu) / (Thickness * Thickness);
 
 /// Pressure magnitude
  double P_mag = 0.0;

 /// Element area
 double Element_area = 0.5;
 
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

  /// Get the null function for applying homogenous BCs
  void null_fct(const Vector<double>& x, double& exact_w)
  {
    exact_w = 0.0;
  }


}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Problem definition
//====================================================================
template<class ELEMENT>
class UnstructuredFvKProblem : public virtual Problem
{

public:

  /// Constructor
  UnstructuredFvKProblem(double const& element_area = 0.09);

  /// Destructor
  ~UnstructuredFvKProblem()
  {
    // Close trace file
    Trace_file.close();
  };

  /// Update after solve (empty)
  void actions_after_newton_solve() {}

  /// Update the problem specs before solve: empty
  void actions_before_newton_solve(){}
 
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

  } // End make_linear()

 
  /// Doc the solution
  void doc_solution(const std::string& comment="");

  /// Overloaded version of the problem's access function to
  /// the mesh. 
  TriangleMesh<ELEMENT>* mesh_pt()
  {
    return Bulk_mesh_pt;
  }


private:

  /// Pin all displacements and rotation (dofs 0-4) at the centre
  void pin_all_displacements_and_rotation_at_centre_node();

  /// Pin all in-plane displacements in the domain
  /// (for solving the linear problem)
  void pin_all_in_plane_displacements();

  /// Trace file to document norm of solution
  ofstream Trace_file;

  /// Pointer to "bulk" mesh
  TriangleMesh<ELEMENT>* Bulk_mesh_pt;

  /// Enumeration to keep track of boundary ids
  enum
  {
    Outer_boundary0 = 0,
    Outer_boundary1 = 1,
    Inner_boundary0 = 2,
    Inner_boundary1 = 3
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
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(const double&
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

   // Map as "sparse vector" for curvilines associated with outer boundaries
   std::map<unsigned,C1CurviLine*> c1_curviline_pt;

   // hierher there should be an interface in the triangle mesh class that
   // returns al of these automatically
   c1_curviline_pt[Outer_boundary0]
    = new C1CurviLine(dynamic_cast<TriangleMeshCurviLine*>(
                       outer_curvilinear_boundary_pt[0]));
   c1_curviline_pt[Outer_boundary1]
    = new C1CurviLine(dynamic_cast<TriangleMeshCurviLine*>(
                       outer_curvilinear_boundary_pt[1])); 
   
   // Split elements that have two boundary edges
   TimeStepper* time_stepper_pt = Bulk_mesh_pt->Time_stepper_pt;
   Bulk_mesh_pt->
    template split_elements_with_multiple_boundary_edges<ELEMENT>(time_stepper_pt);
   
   // Create the mesh for the Lagrange multiplier elements that enforce
   // continuity of our smooth solution across different parts of the
   // mesh boundary (only really needed when there are kinks)
   Constraint_mesh_pt = new Mesh();


   // New general helper function
   C1Helper::duplicate_corner_nodes(Bulk_mesh_pt,
                                    c1_curviline_pt,
                                    Constraint_mesh_pt);

   // New general helper function
   C1Helper::upgrade_edge_elements_to_curved_boundaries(
    Bulk_mesh_pt,
    c1_curviline_pt);

   // Rotate degrees of freedom (only needed for clamped bcs
   if (CommandLineArgs::command_line_flag_has_been_set("--rotate_dofs_on_boundary"))
    {
     C1Helper::rotate_edge_degrees_of_freedom(
      Bulk_mesh_pt,
      c1_curviline_pt);
    }
   
  }
  // End upgrade C1 boundaries

  
  //Add submeshes to problem
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Constraint_mesh_pt);

  // Combine submeshes into a single Mesh 
  build_global_mesh();

  // Complete the build of all elements so they are fully functional
  unsigned n_element = Bulk_mesh_pt->nelement();
  for(unsigned e=0;e<n_element;e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    //Set the pressure function pointers and the physical constants
    el_pt->pressure_fct_pt() = &Parameters::get_pressure;
    el_pt->in_plane_forcing_fct_pt() = &Parameters::get_in_plane_traction;

    // Assign the parameter pointers for the element
    el_pt->nu_pt() = &Parameters::Nu;
    el_pt->eta_pt() = &Parameters::Eta;
  }

  // Set the boundary conditions
  //============================

  // Clamp it
  if (Parameters::Problem_case==Parameters::Clamped_validation)
  {
    // Set the boundary conditions
    unsigned nbound = 2;
    for(unsigned b = 0; b < nbound; b++)
    {
      const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);
      for(unsigned e=0;e<nb_element;e++)
      {
        // Get pointer to bulk element adjacent to b
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b,e));

        // A true clamp, so we set everything except the second normal to zero
        for(unsigned idof=0; idof<6; ++idof)
        {
          // Pin both in plane displacements
          // Cannot set second normal derivative
          if(idof<2)
          {
            el_pt->fix_in_plane_displacement_dof(idof,b,Parameters::null_fct);
          }

          // Cannot set second normal derivative
          if(idof!=3)
          {
            el_pt->fix_out_of_plane_displacement_dof(idof,b,Parameters::null_fct);
          }
        }
      }
    }
  }
  // All other cases: simply pin and stop rotation via the centre
  else
  {
    pin_all_displacements_and_rotation_at_centre_node();
  }
  
  // Update the corner constraints based on the applied
  // boundary conditions
  unsigned n_el = Constraint_mesh_pt->nelement();
  for(unsigned i_el = 0; i_el < n_el; i_el++)
   {
    // hierher rename and unify FvK/KS
    dynamic_cast<DuplicateNodeConstraintElement*>
     (Constraint_mesh_pt->element_pt(i_el))
     ->validate_and_pin_redundant_constraints();
   }
  
  // Assign equation numbers
  oomph_info << "Number of equations: " << assign_eqn_numbers() << '\n';

  
  // Set directory
  Doc_info.set_directory("RESLT");

  // Open trace file
  char filename[100];
  sprintf(filename, "RESLT/trace.dat");
  Trace_file.open(filename);

} // end Constructor






//==start_of_pin_all_displacements_and_rotation_at_centre_node======================
/// pin all displacements and rotations in the centre
//==============================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
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

  // Pin central node:
  // - In-plane dofs are always 0 and 1
  // - Out of plane displacement is 2, x and y derivatives are 3 and 4.
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
  
  // Pin y displacement at max distance from origin
  pinned_rotation_node_pt->pin(1);
}





//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(
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


  // Find the solution at r=0
  // ------------------------

  // should really pre-compute this since it doesn't change...
  MeshAsGeomObject Mesh_as_geom_obj(Bulk_mesh_pt);
  Vector<double> s(2);
  GeomObject* geom_obj_pt=0;
  Vector<double> r(2,0.0);
  Mesh_as_geom_obj.locate_zeta(r,geom_obj_pt,s);

  // Compute the interpolated displacement vector
  Vector<double> u_0(3,0.0);
  u_0=dynamic_cast<ELEMENT*>(geom_obj_pt)->interpolated_fvk_disp(s);
  oomph_info << "w in the middle: " << std::setprecision(15)
             << u_0[2] << std::endl;
  Trace_file << Parameters::P_mag << " " << u_0[2] << '\n';

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

  // Rotate dofs?
  CommandLineArgs::specify_command_line_flag("--rotate_dofs_on_boundary");

  // Parse command line
  CommandLineArgs::parse_and_assign();

  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();

  // Constant pressure for validation case
  if (CommandLineArgs::command_line_flag_has_been_set("--use_clamped_bc"))
  {
    Parameters::Problem_case = Parameters::Clamped_validation;
    if (!CommandLineArgs::command_line_flag_has_been_set("--rotate_dofs_on_boundary"))
     {
      oomph_info << "clamped bcs require rotated dofs on boundary" << std::endl;
      abort();
     }
  }

  // Build problem
  UnstructuredFvKProblem<FoepplVonKarmanC1CurvableBellElement<4>> problem(
    Parameters::Element_area);

  // Tweak Newton solver parameters
  problem.max_residuals() = 1.0e3;
  problem.max_newton_iterations() = 30;
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
