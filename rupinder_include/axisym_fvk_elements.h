#ifndef OOMPH_AXISYM_FVK_ELEMENTS_HEADER
#define OOMPH_AXISYM_FVK_ELEMENTS_HEADER

// Generic oomph-lib routines
#include "generic.h"

// Include the singular integration scheme
#include "singular_integration.h"

namespace oomph
{
  //=============================================================
  /// A class for the element that solves the axisymmetric
  /// Foeppl-von Karman equations and the perturbation
  /// eigenvalue problem.
  //=============================================================
  template<unsigned NNODE_1D>
  class AxisymFvkElement : public QElement<1, NNODE_1D>,
                           public QHermiteElement<1>
  {
  public:
    /// Function pointer to the pressure function
    typedef double (*AxisymFvKPressureFctPt)(const double& r);

    /// Function pointer to the singular pressure function
    typedef double (*AxisymFvKSingularPressureFctPt)(const double& r);

    /// Function pointer to the forcing term in the Cauchy equation
    typedef double (*AxisymFvKForcingFctPt)(const double& r);

    /// Constructor (Initialise all function pointers to 0)
    AxisymFvkElement()
      : Pressure_fct_pt(0), Singular_pressure_fct_pt(0), Forcing_fct_pt(0)
    {
      // Set the number of nodes
      this->set_n_node(NNODE_1D);

      // Set the elemental and nodal dimensions
      this->set_dimension(1);

      // Set the number of interpolated position types (1 for the Lagrangian
      // basis and shape functions)
      this->set_nnodal_position_type(1);

      // Assign pointer to the default integration scheme
      this->set_integration_scheme(&Default_integration_scheme);

      // Assign pointer to the singular integration scheme
      this->set_singular_integration_scheme(&Singular_integration_scheme);
    }

    /// Return the number of Hermite nodes
    unsigned nhermite_node() const
    {
      return 2;
    }

    /// Return a pointer to the nth Hermite node
    /// There are only 2 Hermite nodes in this element so n = 0 calls the first
    /// node and n = 1 calls the last node
    Node*& hermite_node_pt(const unsigned& n)
    {
      // Call the first node if n = 0
      if (n == 0)
      {
        return node_pt(0);
      }
      else if (n == 1)
      {
        // Call the first node if n = 1
        return node_pt(nnode() - 1);
      }

      // We should only get here if n != 0, 1
      std::string error_message =
        "There are only two nodes in this element that contain "
        "Hermite information, you have called node " +
        std::to_string(n) + ".";

      throw OomphLibError(
       error_message,
       OOMPH_EXCEPTION_LOCATION,
       OOMPH_CURRENT_FUNCTION);
    }

    /// Return a pointer to the nth Hermite node (const version)
    Node* const& hermite_node_pt(const unsigned& n) const
    {
      // Call the first node if n = 0
      if (n == 0)
      {
        return node_pt(0);
      }
      else if (n == 1)
      {
        // Call the first node if n = 1
        return node_pt(nnode() - 1);
      }

      // We should only get here if n != 0, 1
      std::string error_message =
        "There are only two nodes in this element that contain "
        "Hermite information, you have called node " +
        std::to_string(n) + ".";

      throw OomphLibError(
       error_message,
       OOMPH_EXCEPTION_LOCATION,
       OOMPH_CURRENT_FUNCTION);
      
    }

    /// Returns true if the nth node is Hermite, false otherwise
    bool is_hermite_node(const unsigned& n) const
    {
      // Only the first and last nodes in the element are Hermite
      if (n == 0 || n == nnode() - 1)
      {
        return true;
      }

      return false;
    }

    /// Required # of `values' (pinned or dofs) at node n
    /// Store values at nodes
    unsigned required_nvalue(const unsigned& n) const override
    {
      // Store u, w, w', \tilde{u}_{r}, \tilde{u}_{\theta}, \tilde{w},
      // \tilde{w}' at the Hermite nodes
      if (n == 0 || n == nnode() - 1)
      {
        return 7;
      }

      // Store u, \tilde{u}_{r}, \tilde{u}_{\theta} at all other nodes
      return 3;
    }

    /// Return the index at which the radial displacement is stored
    unsigned u_index_fvk() const
    {
      return 0;
    }

    /// Return the index at which the radial displacement perturbation is stored
    unsigned u_pert_index_fvk() const
    {
      return 1;
    }

    /// Return the index at which the azimuthal displacement perturbation is
    /// stored
    unsigned u_theta_pert_index_fvk() const
    {
      return 2;
    }

    /// Return the index at which the out of plane displacement is stored
    unsigned w_index_fvk() const
    {
      return 3;
    }

    /// Return the index at which the local derivative of the out of plane
    /// displacement is stored
    unsigned dwds_index_fvk() const
    {
      return 4;
    }

    /// Return the index at which the out of plane displacement perturbation is
    /// stored
    unsigned w_pert_index_fvk() const
    {
      return 5;
    }

    /// Return the index at which the local derivative of the out of plane
    /// displacement perturbation is stored
    unsigned dwds_pert_index_fvk() const
    {
      return 6;
    }

    /// Return the number of types of basis function required to interpolate u.
    unsigned nlagrange_type() const
    {
      // There is only 1 type for Lagrange elements
      return 1;
    }

    /// Return the number of types of basis function required to interpolate w
    unsigned nhermite_type() const
    {
      // There are two types of basis functions for Hermite elements
      // the type with zero slope and unit value,
      // the type with unit slope and zero value.
      return 2;
    }

    /// Return the interpolated value of u at a given local coordinate s
    double interpolated_u(const Vector<double>& s) const
    {
      // Store the result
      double interpolated_u = 0.0;

      // Find out how many nodes there are
      unsigned n_node = nnode();

      // Set up storage for shape functions
      Shape psi(n_node);

      // Get the values of the shape functions at s
      shape(s, psi);

      // Sum up the contributions from the shape function at each node
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += raw_nodal_value(l, u_index_fvk()) * psi(l);
      }

      return interpolated_u;
    }

    /// Return the interpolated value of the perturbed radial displacement at a
    /// given local coordinate s
    double interpolated_u_pert(const Vector<double>& s) const
    {
      // Store the result
      double interpolated_u_pert = 0.0;

      // Find out how many nodes there are
      unsigned n_node = nnode();

      // Set up storage for shape functions
      Shape psi(n_node);

      // Get the values of the shape functions at s
      shape(s, psi);

      // Sum up the contributions from the shape function at each node
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u_pert += raw_nodal_value(l, u_pert_index_fvk()) * psi(l);
      }

      return interpolated_u_pert;
    }

    /// Return the interpolated value of the perturbed azimuthal displacement at
    /// a given local coordinate s
    double interpolated_u_theta_pert(const Vector<double>& s) const
    {
      // Store the result
      double interpolated_u_theta_pert = 0.0;

      // Find out how many nodes there are
      unsigned n_node = nnode();

      // Set up storage for shape functions
      Shape psi(n_node);

      // Get the values of the shape functions at s
      shape(s, psi);

      // Sum up the contributions from the shape function at each node
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u_theta_pert +=
          raw_nodal_value(l, u_theta_pert_index_fvk()) * psi(l);
      }

      return interpolated_u_theta_pert;
    }

    /// Return the interpolated value of w at a given local coordinate s
    double interpolated_w(const Vector<double>& s) const
    {
      // Store the result
      double interpolated_w = 0.0;

      // There are two basis functions for interpolating w
      unsigned n_nodal_position_type = nhermite_type();

      // Get the number of Hermite nodes
      unsigned n_hermite_node = nhermite_node();

      // Set up storage for shape functions
      Shape psi(n_hermite_node, n_nodal_position_type);

      // Get the value of the shape function at s
      QHermiteElement<1>::shape(s, psi);

      // Sum up the contributions from the two shape functions at the first and
      // last node
      for (unsigned l = 0; l < n_hermite_node; l++)
      {
        for (unsigned k = 0; k < n_nodal_position_type; k++)
        {
          // Get the appropriate raw nodal value from the node for each type of
          // function
          // When k = 0, raw_value returns the value of w
          // When k = 1, raw_value returns the value of dwds
          interpolated_w +=
            hermite_node_pt(l)->raw_value(w_index_fvk() + k) * psi(l, k);
        }
      }

      return interpolated_w;
    }

    /// Return the interpolated value of w_pert at a given local coordinate s
    double interpolated_w_pert(const Vector<double>& s) const
    {
      // Store the result
      double interpolated_w_pert = 0.0;

      // There are two basis functions for interpolating w
      unsigned n_nodal_position_type = nhermite_type();

      // Get the number of Hermite nodes
      unsigned n_hermite_node = nhermite_node();

      // Set up storage for shape functions
      Shape psi(n_hermite_node, n_nodal_position_type);

      // Get the value of the shape function at s
      QHermiteElement<1>::shape(s, psi);

      // Sum up the contributions from the two shape functions at the first and
      // last node
      for (unsigned l = 0; l < n_hermite_node; l++)
      {
        for (unsigned k = 0; k < n_nodal_position_type; k++)
        {
          // Get the appropriate raw nodal value from the node for each type of
          // function
          // When k = 0, raw_value returns the value of w_pert
          // When k = 1, raw_value returns the value of dw_pertds
          interpolated_w_pert +=
            hermite_node_pt(l)->raw_value(w_pert_index_fvk() + k) * psi(l, k);
        }
      }

      return interpolated_w_pert;
    }

    /// Return the interpolated value of du/dr at a given local coordinate s
    double interpolated_dudr(const Vector<double>& s) const
    {
      // Store the result
      double interpolated_dudr = 0.0;

      // Find out how many nodes there are
      unsigned n_node = nnode();

      // Set up storage for shape functions
      Shape psi(n_node);
      DShape dpsidr(n_node, 1);

      // Get the values of the shape functions at s
      dshape_eulerian(s, psi, dpsidr);

      // Sum up from the contributions from the shape function at each node
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_dudr += raw_nodal_value(l, u_index_fvk()) * dpsidr(l, 0);
      }

      return interpolated_dudr;
    }

    /// Return the interpolated radial derivative of the perturbation in the
    /// radial displacement
    double interpolated_dudr_pert(const Vector<double>& s) const
    {
      // Store the result
      double interpolated_dudr_pert = 0.0;

      // Find out how many nodes there are
      unsigned n_node = nnode();

      // Set up storage for shape functions
      Shape psi(n_node);
      DShape dpsidr(n_node, 1);

      // Get the values of the shape functions at s
      dshape_eulerian(s, psi, dpsidr);

      // Sum up from the contributions from the shape function at each node
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_dudr_pert +=
          raw_nodal_value(l, u_pert_index_fvk()) * dpsidr(l, 0);
      }

      return interpolated_dudr_pert;
    }

    /// Return the interpolated radial derivative of the perturbed azimuthal
    /// displacement at a given local coordinate s
    double interpolated_dudr_theta_pert(const Vector<double>& s) const
    {
      // Store the result
      double interpolated_dudr_theta_pert = 0.0;

      // Find out how many nodes there are
      unsigned n_node = nnode();

      // Set up storage for shape functions
      Shape psi(n_node);
      DShape dpsidr(n_node, 1);

      // Get the values of the shape functions at s
      dshape_eulerian(s, psi, dpsidr);

      // Sum up the contributions from the shape function at each node
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_dudr_theta_pert +=
          raw_nodal_value(l, u_theta_pert_index_fvk()) * dpsidr(l, 0);
      }

      return interpolated_dudr_theta_pert;
    }


    /// Return the interpolated value of dw/dr at a given local coordinate s
    double interpolated_dwdr(const Vector<double>& s) const
    {
      // Store the result
      double interpolated_dwdr = 0.0;

      // Get the number of Hermite nodes
      unsigned n_hermite_node = nhermite_node();

      // There are two shape functions at each node
      unsigned n_nodal_position_type = nhermite_type();

      // Get the index of the coordinate we are taking the derivative with
      // respect to
      unsigned derivative_coordinate_index = 0;

      // Get the number of dimensions we are in (1D)
      unsigned n_dim = 1;

      // Set up storage for shape functions
      Shape psi(n_hermite_node, n_nodal_position_type);
      DShape dpsidr(n_hermite_node, n_nodal_position_type, n_dim);

      // Get the values of the shape functions at s
      dshape_eulerian_hermite(s, psi, dpsidr);

      // Sum up the contributions from the two shape functions at each node
      for (unsigned l = 0; l < n_hermite_node; l++)
      {
        for (unsigned k = 0; k < n_nodal_position_type; k++)
        {
          // Get the appropriate raw nodal value from the node for each type of
          // function
          // When k = 0, raw_value returns the value of w
          // When k = 1, raw_value returns the value of dwds
          interpolated_dwdr +=
            hermite_node_pt(l)->raw_value(w_index_fvk() + k) *
            dpsidr(l, k, derivative_coordinate_index);
        }
      }

      return interpolated_dwdr;
    }

    /// Return the interpolated value of dw_pert/dr at a given local coordinate
    /// s
    double interpolated_dwdr_pert(const Vector<double>& s) const
    {
      // Store the result
      double interpolated_dwdr_pert = 0.0;

      // Get the number of Hermite nodes
      unsigned n_hermite_node = nhermite_node();

      // There are two shape functions at each node
      unsigned n_nodal_position_type = nhermite_type();

      // Get the index of the coordinate we are taking the derivative with
      // respect to
      unsigned derivative_coordinate_index = 0;

      // Get the number of dimensions we are in (1D)
      unsigned n_dim = 1;

      // Set up storage for shape functions
      Shape psi(n_hermite_node, n_nodal_position_type);
      DShape dpsidr(n_hermite_node, n_nodal_position_type, n_dim);

      // Get the values of the shape functions at s
      dshape_eulerian_hermite(s, psi, dpsidr);

      // Sum up the contributions from the two shape functions at each node
      for (unsigned l = 0; l < n_hermite_node; l++)
      {
        for (unsigned k = 0; k < n_nodal_position_type; k++)
        {
          // Get the appropriate raw nodal value from the node for each type of
          // function
          // When k = 0, raw_value returns the value of w_pert
          // When k = 1, raw_value returns the value of dw_pertds
          interpolated_dwdr_pert +=
            hermite_node_pt(l)->raw_value(w_pert_index_fvk() + k) *
            dpsidr(l, k, derivative_coordinate_index);
        }
      }

      return interpolated_dwdr_pert;
    }

    /// Return the interpolated value of d^2w/dr^2 at a given local coordinate s
    double interpolated_d2wdr(const Vector<double>& s) const
    {
      // Store the result
      double interpolated_d2wdr = 0.0;

      // Get the number of Hermite nodes
      unsigned n_hermite_node = nhermite_node();

      // There are two shape functions at each node
      unsigned n_nodal_position_type = nhermite_type();

      // Get the index of the coordinate we are taking the derivative with
      // respect to
      unsigned derivative_coordinate_index = 0;

      // Get the number of dimensions we are in (1)
      unsigned n_dim = 1;

      // Set up storage for shape functions
      Shape psi(n_hermite_node, n_nodal_position_type);
      DShape dpsidr(n_hermite_node, n_nodal_position_type, n_dim);
      DShape d2psidr(n_hermite_node, n_nodal_position_type, n_dim);

      // Get the values of the shape functions at s
      d2shape_eulerian_hermite(s, psi, dpsidr, d2psidr);

      // Sum up the contributions from the two shape functions at each node
      for (unsigned l = 0; l < n_hermite_node; l++)
      {
        for (unsigned k = 0; k < n_nodal_position_type; k++)
        {
          // Get the appropriate raw nodal value from the node for each type of
          // function
          // When k = 0, raw_value returns the value of w
          // When k = 1, raw_value returns the value of dwds
          interpolated_d2wdr +=
            hermite_node_pt(l)->raw_value(w_index_fvk() + k) *
            d2psidr(l, k, derivative_coordinate_index);
        }
      }

      return interpolated_d2wdr;
    }

    /// Return the interpolated value of d^2w_pert/dr^2 at a given local
    /// coordinate s
    double interpolated_d2wdr_pert(const Vector<double>& s) const
    {
      // Store the result
      double interpolated_d2wdr_pert = 0.0;

      // Get the number of Hermite nodes
      unsigned n_hermite_node = nhermite_node();

      // There are two shape functions at each node
      unsigned n_nodal_position_type = nhermite_type();

      // Get the index of the coordinate we are taking the derivative with
      // respect to
      unsigned derivative_coordinate_index = 0;

      // Get the number of dimensions we are in (1)
      unsigned n_dim = 1;

      // Set up storage for shape functions
      Shape psi(n_hermite_node, n_nodal_position_type);
      DShape dpsidr(n_hermite_node, n_nodal_position_type, n_dim);
      DShape d2psidr(n_hermite_node, n_nodal_position_type, n_dim);

      // Get the values of the shape functions at s
      d2shape_eulerian_hermite(s, psi, dpsidr, d2psidr);

      // Sum up the contributions from the two shape functions at each node
      for (unsigned l = 0; l < n_hermite_node; l++)
      {
        for (unsigned k = 0; k < n_nodal_position_type; k++)
        {
          // Get the appropriate raw nodal value from the node for each type of
          // function
          // When k = 0, raw_value returns the value of w_pert
          // When k = 1, raw_value returns the value of dw_pertds
          interpolated_d2wdr_pert +=
            hermite_node_pt(l)->raw_value(w_pert_index_fvk() + k) *
            d2psidr(l, k, derivative_coordinate_index);
        }
      }

      return interpolated_d2wdr_pert;
    }

    /// Return the von Karman strain in the radial direction
    double radial_strain(const Vector<double>& s) const
    {
      return (interpolated_dudr(s) + 0.5 * std::pow(interpolated_dwdr(s), 2));
    }

    /// Return the azimuthal von Karman strain
    double hoop_strain(const Vector<double>& s) const
    {
      return (interpolated_u(s) / interpolated_x(s, 0));
    }

    /// Return the linearised radial von Karman strain of the perturbation
    double radial_strain_pert(const Vector<double>& s) const
    {
      return interpolated_dudr_pert(s) +
             interpolated_dwdr_pert(s) * interpolated_dwdr(s);
    }

    /// Return the linearised hoop von Karman strain of the perturbation
    double hoop_strain_pert(const Vector<double>& s) const
    {
      return (
        (interpolated_u_pert(s) + wavemode() * interpolated_u_theta_pert(s)) /
        interpolated_x(s, 0));
    }

    /// Return the off diagonal term of the in-plane strain tensor of the
    /// perturbation
    double inplane_shear_strain_pert(const Vector<double>& s) const
    {
      return 0.5 *
             ((-wavemode() * interpolated_u_pert(s) / interpolated_x(s, 0)) +
              interpolated_dudr_theta_pert(s) -
              (wavemode() * interpolated_w_pert(s) * interpolated_dwdr(s) /
               interpolated_x(s, 0)));
    }

    // Return the radial stress
    double radial_stress(const Vector<double>& s) const
    {
      return (1.0 / (1.0 - nu() * nu())) *
             (radial_strain(s) + nu() * hoop_strain(s));
    }

    // Return the radial stress of the perturbation
    double radial_stress_pert(const Vector<double>& s) const
    {
      return (1.0 / (1.0 - nu() * nu())) *
             (radial_strain_pert(s) + nu() * hoop_strain_pert(s));
    }

    // Return the azimuthal stress
    double hoop_stress(const Vector<double>& s) const
    {
      return (1.0 / (1.0 - nu() * nu())) *
             (hoop_strain(s) + nu() * radial_strain(s));
    }

    // Return the azimuthal stress of the perturbation
    double hoop_stress_pert(const Vector<double>& s) const
    {
      return (1.0 / (1.0 - nu() * nu())) *
             (hoop_strain_pert(s) + nu() * radial_strain_pert(s));
    }

    // Return the inplane shear stress of the perturbation
    double inplane_shear_stress_pert(const Vector<double>& s) const
    {
      return (1.0 / (1.0 + nu())) * (inplane_shear_strain_pert(s));
    }

    /// Integrate the pressure over the element
    double integrate_pressure()
    {
      // We are in one dimension
      unsigned n_dim = 1;

      // Storage for variables
      double result = 0.0;
      double pressure = 0.0;
      double w = 0.0;
      double J = 0.0;
      double interpolated_r = 0.0;
      Vector<double> s(1, 0.0);

      // Get the shape functions
      const unsigned n_node = this->nnode();
      Shape psi(n_node);
      DShape dpsidr(n_node, n_dim);

      // Get the number of integration points
      unsigned n_intpt = this->integral_pt()->nweight();

      // Integrate the regular part of the pressure
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the weight of the ipt-th integration point
        w = this->integral_pt()->weight(ipt);

        // Get the Jacobian and shape function at the ipt-th knot
        J = this->dshape_eulerian(s, psi, dpsidr);

        // Get the coordinates of the ipt-th knot
        s[0] = this->integral_pt()->knot(ipt, 0);
        interpolated_r = this->interpolated_x(s, 0);

        // Get the pressure at the ipt-th knot
        pressure = this->get_pressure_fvk(interpolated_r);

        // The pressure is singular on the edge of the domain, if this element
        // is not on the edge then use the regular integration scheme to
        // integrate the pressure.
        if (!element_is_on_edge())
        {
          pressure += this->get_singular_pressure_fvk(interpolated_r);
        }

        // Apply the quadrature rule (Pressure is axisymmetric so multiply by 2
        // pi)
        result +=
          2.0 * MathematicalConstants::Pi * interpolated_r * pressure * J * w;
      }

      // If this element is on the edge of the domain, integrate the singular
      // pressure with the singular integration scheme
      if (element_is_on_edge())
      {
        // Get the number of integration points
        n_intpt = this->singular_integral_pt()->nweight();

        // Loop over the integration points
        for (unsigned ipt = 0; ipt < n_intpt; ipt++)
        {
          // Get the weight of the ipt-th integration point
          w = this->singular_integral_pt()->weight(ipt);

          // Get the local coordinates of the ipt-th knot
          s[0] = this->singular_integral_pt()->knot(ipt, 0);
          interpolated_r = this->interpolated_x(s, 0);

          // Get the Jacobian and shape function at the ipt-th knot
          J = this->dshape_eulerian(s, psi, dpsidr);

          // Get the singular pressure at the ipt-th knot
          pressure = this->get_singular_pressure_fvk(interpolated_r);

          // Apply the quadrature rule (Pressure is axisymmetric so multiply by
          // 2 pi)
          result +=
            2.0 * MathematicalConstants::Pi * interpolated_r * pressure * J * w;
        }
      }
      return result;
    }

    /// Add the element's contribution to its residual vector
    void fill_in_contribution_to_residuals(Vector<double>& residuals) override
    {
      // Get the number of nodes
      unsigned n_node = nnode();

      // Get the number of Hermite nodes
      unsigned n_hermite_node = nhermite_node();

      // Store the number of basis functions at each node
      unsigned n_hermite_type = nhermite_type();

      // Store the number of dimensions of the problem
      unsigned n_dim = 1;

      // Store the index of the Hermite node disregarding non-Hermite nodes
      // i.e. if the n-th node is the second Hermite node, the
      // hermite_node_index value is 1.
      unsigned hermite_node_index = 0;

      // Set up memory for the shape functions and their derivatives
      Shape hermite_psi(n_hermite_node, n_hermite_type);
      DShape hermite_dpsidr(n_hermite_node, n_hermite_type, n_dim);
      DShape hermite_d2psidr(n_hermite_node, n_hermite_type, n_dim);
      Shape lagrange_psi(n_node);
      DShape lagrange_dpsidr(n_node, n_dim);

      // Set the number of integration points
      unsigned n_intpt = integral_pt()->nweight();

      // Store the local equation number
      int local_eqn = 0;

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the integral weight
        double w = integral_pt()->weight(ipt);

        // Get the integral knot position
        Vector<double> s(1, integral_pt()->knot(ipt, 0));

        // Call the derivatives of the Hermite basis functions
        double J = d2shape_eulerian_hermite(
          s, hermite_psi, hermite_dpsidr, hermite_d2psidr);

        // Call the derivatives of the Lagrange basis functions
        dshape_eulerian_at_knot(ipt, lagrange_psi, lagrange_dpsidr);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Store values interpolated at the current integration point
        double interpolated_r = this->interpolated_x(s, 0);
        double interpolated_dwdr = this->interpolated_dwdr(s);
        double interpolated_d2wdr = this->interpolated_d2wdr(s);

        // Assemble the contributions to the residual
        // Loop over the Lagrange test functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Get the local equation number
          local_eqn = nodal_local_eqn(l, u_index_fvk());

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            // The steady Cauchy equation with a forcing term
            // (Strong form: \nabla \cdot \sigma = F)
            // (Weak form: \int (\sigma : \nabla \psi + \psi F) r dr)
            residuals[local_eqn] +=
              (lagrange_dpsidr(l, 0) * radial_stress(s) * interpolated_r +
               lagrange_psi(l) * hoop_stress(s) +
               lagrange_psi(l) * get_forcing_fvk(interpolated_r) *
                 interpolated_r) *
              W;
          }

          // If looping over Hermite test functions
          if (is_hermite_node(l))
          {
            // The first and last nodes correspond to the zeroth and first shape
            // functions respectively
            if (l == 0)
            {
              hermite_node_index = 0;
            }
            else
            {
              hermite_node_index = 1;
            }

            // Loop over the types of Hermite test function
            for (unsigned k = 0; k < n_hermite_type; k++)
            {
              // Get the local equation number
              local_eqn = nodal_local_eqn(l, w_index_fvk() + k);

              // If it's not a boundary condition
              if (local_eqn >= 0)
              {
                // The axisymmetric Laplacian operator
                residuals[local_eqn] +=
                  ((interpolated_dwdr / interpolated_r) + interpolated_d2wdr) *
                  ((hermite_dpsidr(hermite_node_index, k, 0) / interpolated_r) +
                   hermite_d2psidr(hermite_node_index, k, 0)) *
                  interpolated_r * W;

                // Inplane stress terms
                residuals[local_eqn] +=
                  eta() * radial_stress(s) *
                  hermite_dpsidr(hermite_node_index, k, 0) * interpolated_dwdr *
                  interpolated_r * W;

                // Regular pressure residual
                residuals[local_eqn] -= get_pressure_fvk(interpolated_r) *
                                        hermite_psi(hermite_node_index, k) *
                                        interpolated_r * W;

                // If not on the edge, use the same integral scheme to integrate
                // the singular pressure field
                if (!element_is_on_edge())
                {
                  residuals[local_eqn] -=
                    get_singular_pressure_fvk(interpolated_r) *
                    hermite_psi(hermite_node_index, k) * interpolated_r * W;
                }
              }
            }
          }
        }
      }

      // Use the singular integration scheme to add the residuals from the
      // singular pressure if on the edge
      //------------------------------------------------------------------

      if (element_is_on_edge())
      {
        // Set the number of integration points
        n_intpt = singular_integral_pt()->nweight();

        // Loop over the integration points
        for (unsigned ipt = 0; ipt < n_intpt; ipt++)
        {
          // Get the integral weight
          double w = singular_integral_pt()->weight(ipt);

          // Get the integral knot position
          Vector<double> s(1, singular_integral_pt()->knot(ipt, 0));

          // Call the derivatives of the Hermite basis functions
          double J = dshape_eulerian_hermite(s, hermite_psi, hermite_dpsidr);

          // Premultiply the weights and the Jacobian
          double W = w * J;

          // Store the interpolated position at the integration point
          double interpolated_r = this->interpolated_x(s, 0);

          // Assemble the contributions to the residual
          // Loop over the Hermite test functions
          for (unsigned l = 0; l < n_node; l++)
          {
            // Do nothing if looping over a Lagrange test function
            if (is_hermite_node(l))
            {
              // The first and last nodes correspond to the zeroth and first
              // shape functions respectively
              if (l == 0)
              {
                hermite_node_index = 0;
              }
              else
              {
                hermite_node_index = 1;
              }

              // Loop over the types of Hermite test function
              for (unsigned k = 0; k < n_hermite_type; k++)
              {
                // Get the local equation number
                local_eqn = nodal_local_eqn(l, w_index_fvk() + k);

                // If it's not a boundary condition
                if (local_eqn >= 0)
                {
                  // Singular pressure residual
                  residuals[local_eqn] -=
                    get_singular_pressure_fvk(interpolated_r) *
                    hermite_psi(hermite_node_index, k) * interpolated_r * W;
                }
              }
            }
          }
        }
      }
    }

    // Fill in the contributions to the Jacobian for use in the eigenvalue
    // problem, the mass matrix is the identity matrix
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix) override
    {
      // Find out how many nodes there are
      unsigned n_node = nnode();

      // Get the number of Hermite nodes
      unsigned n_hermite_node = nhermite_node();

      // Store the number of types of shape function for the Hermite
      // interpolated values
      unsigned n_hermite_type = nhermite_type();

      // Store the number of dimensions of the problem
      unsigned n_dim = 1;

      // Store the index of the Hermite node disregarding non-Hermite nodes
      // i.e. if the n-th node is the second Hermite node, the
      // hermite_node_index value is 1.
      unsigned hermite_node_index = 0;
      unsigned hermite_node_index2 = 0;

      // Set up memory for the shape functions and their derivatives
      Shape lagrange_psi(n_node);
      DShape lagrange_dpsidr(n_node, n_dim);
      Shape hermite_psi(n_hermite_node, n_hermite_type);
      DShape hermite_dpsidr(n_hermite_node, n_hermite_type, n_dim);
      DShape hermite_d2psidr(n_hermite_node, n_hermite_type, n_dim);

      // Set the number of integration points
      unsigned n_intpt = integral_pt()->nweight();

      // Integers to store the local equation and unknown numbers
      int local_eqn = 0;
      int local_unknown = 0;

      // Set the mass matrix as the identity matrix
      // Loop through every entry of the mass matrix and set to zero, unless on
      // the diagonal when we set the entry to 1

      // Loop through every degree of freedom in every node
      for (unsigned l = 0; l < n_node; l++)
      {
        unsigned ndof = required_nvalue(l);
        for (unsigned dof = 0; dof < ndof; dof++)
        {
          // Get the local equation number associated with the current degree of
          // freedom
          local_eqn = nodal_local_eqn(l, dof);

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            // Loop through every degree of freedom in every node again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              unsigned ndof2 = required_nvalue(l2);
              for (unsigned dof2 = 0; dof2 < ndof2; dof2++)
              {
                // Get the local unknown number associated with the current
                // degree of freedom
                local_unknown = nodal_local_eqn(l2, dof2);

                // If at a non-zero degree of freedom
                if (local_unknown >= 0)
                {
                  // If the local equation number and local unknown number
                  // match, we are on the diagonal of the mass matrix so set the
                  // entry to 1
                  if (local_eqn == local_unknown)
                  {
                    // This is the local mass matrix, when the full mass matrix
                    // is assembled, the nodes shared between two elements have
                    // the entries from the local mass matrices corresponding to
                    // both elements are summed together. Therefore we set the
                    // local mass matrix terms to 0.5 if the current node
                    // belongs to two elements.

                    // Check if the current node could be shared between two
                    // elements
                    if (l == 0 || l == n_node - 1)
                    {
                      // Check if the node is on the edge of the domain and
                      // therefore not shared between two elements. (This is
                      // dependent on the mesh being between 0 and 1)
                      if (abs(node_pt(l)->x(0)) < 1.0e-8 ||
                          abs(node_pt(l)->x(0) - 1.0) < 1.0e-8)
                      {
                        mass_matrix(local_eqn, local_unknown) = 1.0;
                      }
                      else
                      {
                        mass_matrix(local_eqn, local_unknown) = 0.5;
                      }
                    }
                    // The current node is not shared between elements so set
                    // the entry to 1
                    else
                    {
                      mass_matrix(local_eqn, local_unknown) = 1.0;
                    }
                  }
                  else
                  {
                    mass_matrix(local_eqn, local_unknown) = 0.0;
                  }
                }
              }
            }
          }
        }
      }


      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the integral weight
        double w = integral_pt()->weight(ipt);

        // Get the integral knot position
        Vector<double> s(1, integral_pt()->knot(ipt, 0));

        // Get the radial coordinate
        double interpolated_r = interpolated_x(s, 0);

        // Call the derivatives of the Hermite basis functions
        double J = d2shape_eulerian_hermite(
          s, hermite_psi, hermite_dpsidr, hermite_d2psidr);

        // Call the derivatives of the Lagrange basis functions
        dshape_eulerian_at_knot(ipt, lagrange_psi, lagrange_dpsidr);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Assemble the contributions to the Jacobian and mass matrices
        // Loop over the test functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Assemble Cauchy equation section of the Jacobian
          //-------------------------------------------------

          // Get the local equation number
          local_eqn = this->nodal_local_eqn(l, u_pert_index_fvk());

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Assemble the part of the Jacobian corresponding to the first
              // equation of the Cauchy equations
              //-------------------------------------------------------------

              // Get the local unknown number
              local_unknown = nodal_local_eqn(l2, u_pert_index_fvk());

              // If at a non-zero degree of freedom, add in the entry
              if (local_unknown >= 0)
              {
                // The radial displacement part of the radial stress term
                jacobian(local_eqn, local_unknown) +=
                  (1.0 / (1.0 - nu() * nu())) *
                  (lagrange_dpsidr(l, 0) *
                   (interpolated_r * lagrange_dpsidr(l2, 0) +
                    nu() * lagrange_psi(l2))) *
                  W;

                // The radial displacement part of the inplane shear stress term
                jacobian(local_eqn, local_unknown) +=
                  ((wavemode() * wavemode()) / (2.0 * (1.0 + nu()))) *
                  ((lagrange_psi(l) * lagrange_psi(l2)) / interpolated_r) * W;

                // The radial displacement part of the hoop stress term
                jacobian(local_eqn, local_unknown) +=
                  (1.0 / (1.0 - nu() * nu())) * lagrange_psi(l) *
                  (nu() * lagrange_dpsidr(l2, 0) +
                   lagrange_psi(l2) / interpolated_r) *
                  W;
              }

              // Get the local unknown number
              local_unknown = nodal_local_eqn(l2, u_theta_pert_index_fvk());

              // If at a non-zero degree of freedom, add in the entry
              if (local_unknown >= 0)
              {
                // The azimuthal displacement part of the radial stress term
                jacobian(local_eqn, local_unknown) +=
                  ((wavemode() * nu()) / (1.0 - nu() * nu())) *
                  lagrange_psi(l2) * lagrange_dpsidr(l, 0) * W;

                // The azimuthal displacement part of the inplane shear stress
                // term
                jacobian(local_eqn, local_unknown) -=
                  (wavemode() / (2.0 * (1.0 + nu()))) * lagrange_psi(l) *
                  lagrange_dpsidr(l2, 0) * W;

                // The azimuthal displacement part of the hoop stress term
                jacobian(local_eqn, local_unknown) +=
                  (wavemode() / (1.0 - nu() * nu())) *
                  (lagrange_psi(l) * lagrange_psi(l2) / interpolated_r) * W;
              }

              if (is_hermite_node(l2))
              {
                // The first and last nodes correspond to the zeroth and first
                // shape functions respectively
                if (l2 == 0)
                {
                  hermite_node_index2 = 0;
                }
                else
                {
                  hermite_node_index2 = 1;
                }

                // Loop over the Hermite terms
                for (unsigned k2 = 0; k2 < n_hermite_type; k2++)
                {
                  // Get the local unknown number
                  local_unknown = nodal_local_eqn(l2, w_pert_index_fvk() + k2);

                  // If at a non-zero degree of freedom, add in the entry
                  if (local_unknown >= 0)
                  {
                    // The transverse displacement part of the radial stress
                    // term
                    jacobian(local_eqn, local_unknown) +=
                      (1.0 / (1.0 - nu() * nu())) *
                      (interpolated_dwdr(s) * lagrange_dpsidr(l, 0) *
                       hermite_dpsidr(hermite_node_index2, k2, 0)) *
                      interpolated_r * W;

                    // The transverse displacement part of the inplane shear
                    // stress term
                    jacobian(local_eqn, local_unknown) +=
                      ((wavemode() * wavemode()) / (2.0 * (1.0 + nu()))) *
                      (lagrange_psi(l) * interpolated_dwdr(s) *
                       hermite_psi(hermite_node_index2, k2) / interpolated_r) *
                      W;

                    // The transverse displacement part of the hoop stress term
                    jacobian(local_eqn, local_unknown) +=
                      (nu() / (1.0 - nu() * nu())) *
                      (interpolated_dwdr(s) * lagrange_psi(l) *
                       hermite_dpsidr(hermite_node_index2, k2, 0)) *
                      W;
                  }
                }
              }
            }
          }

          // Assemble the part of the Jacobian corresponding to the second
          // equation of the Cauchy equations
          //--------------------------------------------------------------

          // Get the local equation number
          local_eqn = this->nodal_local_eqn(l, u_theta_pert_index_fvk());

          // If it's not a boundary condition
          if (local_eqn >= 0)
          {
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              // Get the local unknown number
              local_unknown = nodal_local_eqn(l2, u_pert_index_fvk());

              // If at a non-zero degree of freedom, add in the entry
              if (local_unknown >= 0)
              {
                // The two radial displacement parts of the inplane shear stress
                // term
                jacobian(local_eqn, local_unknown) +=
                  (wavemode() / (2.0 * (1.0 + nu()))) *
                  (lagrange_psi(l) * lagrange_psi(l2) / interpolated_r) * W;

                jacobian(local_eqn, local_unknown) -=
                  (wavemode() / (2.0 * (1.0 + nu()))) *
                  (lagrange_dpsidr(l, 0) * lagrange_psi(l2)) * W;

                // The radial displacement part of the hoop stress term
                jacobian(local_eqn, local_unknown) +=
                  (wavemode() / (1.0 - nu() * nu())) * lagrange_psi(l) *
                  (nu() * lagrange_dpsidr(l2, 0) +
                   lagrange_psi(l2) / interpolated_r) *
                  W;
              }

              // Get the local unknown number
              local_unknown = nodal_local_eqn(l2, u_theta_pert_index_fvk());

              // If at a non-zero degree of freedom, add in the entry
              if (local_unknown >= 0)
              {
                // The two azimuthal displacement parts of the inplane shear
                // stress term
                jacobian(local_eqn, local_unknown) -=
                  (1.0 / (2.0 * (1.0 + nu()))) * lagrange_psi(l) *
                  lagrange_dpsidr(l2, 0) * W;

                jacobian(local_eqn, local_unknown) +=
                  (1.0 / (2.0 * (1.0 + nu()))) *
                  (lagrange_dpsidr(l, 0) * lagrange_dpsidr(l2, 0)) *
                  interpolated_r * W;

                // The azimuthal displacement part of the hoop stress term
                jacobian(local_eqn, local_unknown) +=
                  (wavemode() * wavemode() / (1.0 - nu() * nu())) *
                  (lagrange_psi(l) * lagrange_psi(l2) / interpolated_r) * W;
              }

              if (is_hermite_node(l2))
              {
                // The first and last nodes correspond to the zeroth and first
                // shape
                // functions respectively
                if (l2 == 0)
                {
                  hermite_node_index2 = 0;
                }
                else
                {
                  hermite_node_index2 = 1;
                }
                // Loop over the Hermite terms
                for (unsigned k2 = 0; k2 < n_hermite_type; k2++)
                {
                  // Get the local unknown number
                  local_unknown = nodal_local_eqn(l2, w_pert_index_fvk() + k2);

                  // If at a non-zero degree of freedom, add in the entry
                  if (local_unknown >= 0)
                  {
                    // The two out of plane displacement parts of the inplane
                    // shear stress term
                    jacobian(local_eqn, local_unknown) +=
                      (wavemode() / (2.0 * (1.0 + nu()))) *
                      (lagrange_psi(l) * interpolated_dwdr(s) *
                       hermite_psi(hermite_node_index2, k2) / interpolated_r) *
                      W;

                    jacobian(local_eqn, local_unknown) -=
                      (wavemode() / (2.0 * (1.0 + nu()))) *
                      (lagrange_dpsidr(l, 0) * interpolated_dwdr(s) *
                       hermite_psi(hermite_node_index2, k2)) *
                      W;

                    // The out of plane displacement part of the hoop stress
                    // term
                    jacobian(local_eqn, local_unknown) +=
                      (nu() * wavemode() / (1.0 - nu() * nu())) *
                      lagrange_psi(l) * interpolated_dwdr(s) *
                      hermite_dpsidr(hermite_node_index2, k2, 0) * W;
                  }
                }
              }
            }
          }

          // Assemble bending equation section of the Jacobian
          //--------------------------------------------------

          if (is_hermite_node(l))
          {
            // The first and last nodes correspond to the zeroth and first
            // shape
            // functions respectively
            if (l == 0)
            {
              hermite_node_index = 0;
            }
            else
            {
              hermite_node_index = 1;
            }

            for (unsigned k = 0; k < n_hermite_type; k++)
            {
              // Get the local equation number
              local_eqn = this->nodal_local_eqn(l, w_pert_index_fvk() + k);

              // If it's not a boundary condition
              if (local_eqn >= 0)
              {
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  // Get the local unknown number
                  local_unknown = nodal_local_eqn(l2, u_pert_index_fvk());

                  // If at a non-zero degree of freedom, add in the entry
                  if (local_unknown >= 0)
                  {
                    // The radial displacement part of the radial stress term
                    jacobian(local_eqn, local_unknown) +=
                      (eta() / (1.0 - nu() * nu())) * interpolated_dwdr(s) *
                      hermite_dpsidr(hermite_node_index, k, 0) *
                      (lagrange_dpsidr(l2, 0) * interpolated_r +
                       nu() * lagrange_psi(l2)) *
                      W;

                    // The radial displacement part of the inplane shear stress
                    // term
                    jacobian(local_eqn, local_unknown) +=
                      ((eta() * wavemode() * wavemode()) /
                       (2.0 * (1.0 + nu()))) *
                      hermite_psi(hermite_node_index, k) *
                      interpolated_dwdr(s) *
                      (lagrange_psi(l2) / interpolated_r) * W;
                  }

                  // Get the local unknown number
                  local_unknown = nodal_local_eqn(l2, u_theta_pert_index_fvk());

                  // If at a non-zero degree of freedom, add in the entry
                  if (local_unknown >= 0)
                  {
                    // The azimuthal displacement part of the radial stress term
                    jacobian(local_eqn, local_unknown) +=
                      (eta() * wavemode() * nu() / (1.0 - nu() * nu())) *
                      interpolated_dwdr(s) *
                      hermite_dpsidr(hermite_node_index, k, 0) *
                      lagrange_psi(l2) * W;

                    // The azimuthal displacement part of the inplane shear
                    // stress term
                    jacobian(local_eqn, local_unknown) -=
                      (eta() * wavemode() / (2.0 * (1.0 + nu()))) *
                      interpolated_dwdr(s) *
                      hermite_psi(hermite_node_index, k) *
                      lagrange_dpsidr(l2, 0) * W;
                  }

                  if (is_hermite_node(l2))
                  {
                    // The first and last nodes correspond to the zeroth and
                    // first shape functions respectively
                    if (l2 == 0)
                    {
                      hermite_node_index2 = 0;
                    }
                    else
                    {
                      hermite_node_index2 = 1;
                    }

                    for (unsigned k2 = 0; k2 < n_hermite_type; k2++)
                    {
                      // Get the local unknown number
                      local_unknown =
                        nodal_local_eqn(l2, w_pert_index_fvk() + k2);

                      // If at a non-zero degree of freedom, add in the entry
                      if (local_unknown >= 0)
                      {
                        // The biharmonic term
                        jacobian(local_eqn, local_unknown) +=
                          ((hermite_dpsidr(hermite_node_index, k, 0) /
                            interpolated_r) +
                           hermite_d2psidr(hermite_node_index, k, 0) -
                           (wavemode() * wavemode() /
                            (interpolated_r * interpolated_r)) *
                             hermite_psi(hermite_node_index, k)) *
                          ((hermite_dpsidr(hermite_node_index2, k2, 0) /
                            interpolated_r) +
                           hermite_d2psidr(hermite_node_index2, k2, 0) -
                           (wavemode() * wavemode() /
                            (interpolated_r * interpolated_r)) *
                             hermite_psi(hermite_node_index2, k2)) *
                          interpolated_r * W;

                        // The following terms come from the nonlinear term in
                        // the FvK equation
                        jacobian(local_eqn, local_unknown) +=
                          (eta() / (1.0 - nu() * nu())) * interpolated_dwdr(s) *
                          interpolated_dwdr(s) *
                          hermite_dpsidr(hermite_node_index, k, 0) *
                          hermite_dpsidr(hermite_node_index2, k2, 0) *
                          interpolated_r * W;

                        jacobian(local_eqn, local_unknown) +=
                          ((eta() * wavemode() * wavemode()) /
                           (2.0 * (1.0 + nu()))) *
                          interpolated_dwdr(s) * interpolated_dwdr(s) *
                          (hermite_psi(hermite_node_index, k) *
                           hermite_psi(hermite_node_index2, k2) /
                           interpolated_r) *
                          W;

                        jacobian(local_eqn, local_unknown) +=
                          eta() * radial_stress(s) *
                          hermite_dpsidr(hermite_node_index, k, 0) *
                          hermite_dpsidr(hermite_node_index2, k2, 0) *
                          interpolated_r * W;

                        jacobian(local_eqn, local_unknown) +=
                          (eta() * wavemode() * wavemode()) * hoop_stress(s) *
                          (hermite_psi(hermite_node_index, k) *
                           hermite_psi(hermite_node_index2, k2) /
                           interpolated_r) *
                          W;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    /// Compute the Hermite basis functions and also first derivatives w.r.t
    /// global coordinates at local coordinate s;
    /// Returns Jacobian of mapping from global to local coordinates.
    double dshape_eulerian_hermite(const Vector<double>& s,
                                   Shape& psi,
                                   DShape& dpsids) const
    {
      // Find the element dimension
      const unsigned el_dim = dim();

      // Get the number of nodes
      const unsigned n_node = nnode();

      // Get the values of the basis functions and their local derivatives
      dshape_local_hermite(s, psi, dpsids);

      // Allocate memory for the Jacobian and inverse Jacobian
      DenseMatrix<double> jacobian(el_dim);
      DenseMatrix<double> inverse_jacobian(el_dim);

      // Get the shape functions
      Shape shape_psi(n_node);
      DShape shape_dpsids(n_node, el_dim);

      // Get the values of the shape functions and their local derivatives
      dshape_local(s, shape_psi, shape_dpsids);

      // Now calculate the inverse jacobian
      const double det =
        local_to_eulerian_mapping(shape_dpsids, jacobian, inverse_jacobian);

      // Now set the values of the derivatives to be dpsidx
      transform_derivatives(inverse_jacobian, dpsids);

      // Return the determinant of the jacobian
      return det;
    }

    /// Compute the Hermite basis functions and also the first and second
    /// derivatives w.r.t global coordinates at local coordinate s; Returns
    /// Jacobian of mapping from global to local coordinates.
    double d2shape_eulerian_hermite(const Vector<double>& s,
                                    Shape& psi,
                                    DShape& dpsids,
                                    DShape& d2psids) const
    {
      // Find the element dimension
      const unsigned el_dim = dim();

      // Get the number of nodes
      const unsigned n_node = nnode();

      // Get the number of second derivatives required (Only 1 in 1D)
      const unsigned n_deriv = 1;

      // Get the values of the basis functions and their local derivatives
      d2shape_local_hermite(s, psi, dpsids, d2psids);

      // Allocate memory for the Jacobian and inverse Jacobian
      DenseMatrix<double> jacobian(el_dim);
      DenseMatrix<double> inverse_jacobian(el_dim);

      // Get the shape functions
      Shape shape_psi(n_node);
      DShape shape_dpsids(n_node, el_dim);
      DShape shape_d2psids(n_node, el_dim);

      // Get the values of the shape functions and their local derivatives
      d2shape_local(s, shape_psi, shape_dpsids, shape_d2psids);

      // Now calculate the inverse jacobian
      const double det =
        local_to_eulerian_mapping(shape_dpsids, jacobian, inverse_jacobian);

      // Allocate memory for the Jacobian of second derivatives
      DenseMatrix<double> jacobian2(n_deriv, el_dim);

      // Assemble the Jacobian of second derivatives
      assemble_local_to_eulerian_jacobian2(shape_d2psids, jacobian2);

      // Now set the values of the derivatives to be dpsidx
      transform_second_derivatives(
        jacobian, inverse_jacobian, jacobian2, dpsids, d2psids);

      // Return the determinant of the jacobian
      return det;
    }

    /// Local derivatives of the 1D Hermite basis functions
    void dshape_local_hermite(const Vector<double>& s,
                              Shape& psi,
                              DShape& dpsids) const
    {
      // Local storage
      double Psi[2][2], DPsi[2][2];

      // Call the OneDimensional Shape functions
      OneDimHermite::shape(s[0], Psi);
      OneDimHermite::dshape(s[0], DPsi);

      // Loop over number of nodes
      for (unsigned l = 0; l < 2; l++)
      {
        // Loop over number of dofs
        for (unsigned k = 0; k < 2; k++)
        {
          psi(l, k) = Psi[l][k];
          dpsids(l, k, 0) = DPsi[l][k];
        }
      }
    }

    /// Local first and second local derivatives of the 1D Hermite basis
    /// functions
    void d2shape_local_hermite(const Vector<double>& s,
                               Shape& psi,
                               DShape& dpsids,
                               DShape& d2psids) const
    {
      // Local storage
      double Psi[2][2], DPsi[2][2], D2Psi[2][2];

      // Call the OneDimensional Shape functions
      OneDimHermite::shape(s[0], Psi);
      OneDimHermite::dshape(s[0], DPsi);
      OneDimHermite::d2shape(s[0], D2Psi);

      // Loop over number of nodes
      for (unsigned l = 0; l < 2; l++)
      {
        // Loop over number of dofs
        for (unsigned k = 0; k < 2; k++)
        {
          psi(l, k) = Psi[l][k];
          dpsids(l, k, 0) = DPsi[l][k];
          d2psids(l, k, 0) = D2Psi[l][k];
        }
      }
    }

    //======================================================================
    /// Output exact solution
    ///
    /// Solution is provided via function pointer.
    /// Plot at a given number of plot points.
    ///
    /// r, u_exact, w_exact
    //======================================================================
    void output_fct(
      std::ostream& outfile,
      const unsigned& nplot,
      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt) override
    {
      // Vector of local coordinates
      Vector<double> s(1);

      // Vector for coordintes
      Vector<double> r(1);

      // Exact solution Vector
      Vector<double> exact_soln(2);

      // Loop over plot points
      unsigned num_plot_points = nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        get_s_plot(iplot, nplot, s);

        // Get x position as Vector
        interpolated_x(s, r);

        // Get exact solution at this point
        (*exact_soln_pt)(r, exact_soln);

        // Output r, u_exact, w_exact
        outfile << r[0] << " " << exact_soln[0] << " " << exact_soln[1]
                << std::endl;
      }

      // Write tecplot footer (e.g. FE connectivity lists)
      write_tecplot_zone_footer(outfile, nplot);
    }

    /// Output the global coordinate then value of w for nplot plot points in
    /// the element
    void output(std::ostream& outfile, const unsigned& nplot) override
    {
      // Vector of local coordinates (1 dimensional)
      Vector<double> s(1);

      // Vector of global coordinates (1 dimensional)
      Vector<double> x(1);

      // Vector of solution (s, x, u, du/dr, w, dw/dr, d^2w/dr^2, u_pert,
      // dudr_pert, u_theta_pert, dudr_theta_pert, w_pert, dw_pert/dr)
      Vector<double> soln(14);

      // Loop over plot points
      unsigned num_plot_points = nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        get_s_plot(iplot, nplot, s);

        // Get the coordinates of x
        interpolated_x(s, x);

        // Get the local coordinate
        soln[0] = s[0];

        // Get the radial coordinate
        soln[1] = x[0];

        // Get the value of u
        soln[2] = interpolated_u(s);

        // Get the value of du/dr
        soln[3] = interpolated_dudr(s);

        // Get the value of u_pert
        soln[4] = interpolated_u_pert(s);

        // Get the radial derivative of u_pert
        soln[5] = interpolated_dudr_pert(s);

        // Get the value of u_theta_pert
        soln[6] = interpolated_u_theta_pert(s);

        // Get the radial derivative of u_theta_pert
        soln[7] = interpolated_dudr_theta_pert(s);

        // Get the value of w
        soln[8] = interpolated_w(s);

        // Get the value of dw/dr
        soln[9] = interpolated_dwdr(s);

        // Get the value of d^2w/dr^2
        soln[10] = interpolated_d2wdr(s);

        // Get the value of w_pert
        soln[11] = interpolated_w_pert(s);

        // Get the value of dw_pert/dr
        soln[12] = interpolated_dwdr_pert(s);

        // Get the value of d^2_pert/dr^2
        soln[13] = interpolated_d2wdr_pert(s);

        // Output s, r, u, du/dr, u_pert, du_pert/dr, u_theta_pert,
        // du_theta_pertdr w, dw/dr, d^2w/dr^2, w_pert, dw_pertdr, d^2w_pertdr^2
        outfile << soln[0] << " " << soln[1] << " " << soln[2] << " " << soln[3]
                << " " << soln[4] << " " << soln[5] << " " << soln[6] << " "
                << soln[7] << " " << soln[8] << " " << soln[9] << " "
                << soln[10] << " " << soln[11] << " " << soln[12] << " "
                << soln[13] << std::endl;
      }
    }

    /// Get the pressure at radial position r
    double get_pressure_fvk(const double& r) const
    {
      // If no pressure function has been set, return zero
      if (Pressure_fct_pt == 0)
      {
        return 0.0;
      }
      else
      {
        // Get the pressure
        return (*Pressure_fct_pt)(r);
      }
    }

    /// Get the singular pressure at radial position r
    double get_singular_pressure_fvk(const double& r) const
    {
      // If no singular pressure function has been set, return zero
      if (Singular_pressure_fct_pt == 0)
      {
        return 0.0;
      }
      else
      {
        // Get the singular pressure
        return (*Singular_pressure_fct_pt)(r);
      }
    }

    /// Get the forcing term in the Cauchy equation (Used simply for validation)
    double get_forcing_fvk(const double& r) const
    {
      // If no forcing term has been set, return zero
      if (Forcing_fct_pt == 0)
      {
        return 0.0;
      }
      else
      {
        // Get the forcing term
        return (*Forcing_fct_pt)(r);
      }
    }

    // Pointer to Eta
    double*& eta_pt()
    {
      return Eta_pt;
    }

    // Access function to Eta
    const double& eta() const
    {
      return *Eta_pt;
    }

    // Pointer to Nu
    double*& nu_pt()
    {
      return Nu_pt;
    }

    // Access function to Nu
    const double& nu() const
    {
      return *Nu_pt;
    }

    // Pointer to Wavemode
    unsigned*& wavemode_pt()
    {
      return Wavemode_pt;
    }

    // Access function to wavemode
    const unsigned& wavemode() const
    {
      return *Wavemode_pt;
    }

    /// Access function: Pointer to the pressure function
    AxisymFvKPressureFctPt& pressure_fct_pt()
    {
      return Pressure_fct_pt;
    }

    /// Access function: Pointer to the pressure function (const version)
    AxisymFvKPressureFctPt pressure_fct_pt() const
    {
      return Pressure_fct_pt;
    }

    /// Access function: Pointer to the singular pressure function
    AxisymFvKSingularPressureFctPt& singular_pressure_fct_pt()
    {
      return Singular_pressure_fct_pt;
    }

    /// Access function: Pointer to the singular pressure function (const
    /// version)
    AxisymFvKSingularPressureFctPt singular_pressure_fct_pt() const
    {
      return Singular_pressure_fct_pt;
    }

    /// Access function: Pointer to the forcing term
    AxisymFvKForcingFctPt& forcing_fct_pt()
    {
      return Forcing_fct_pt;
    }

    /// Access function: Pointer to the forcing term (const version)
    AxisymFvKForcingFctPt forcing_fct_pt() const
    {
      return Forcing_fct_pt;
    }

    /// Set the singular integration scheme
    void set_singular_integration_scheme(Integral* const& singular_integral_pt)
    {
      Singular_integral_pt = singular_integral_pt;
    }

    /// Return the pointer to the singular integration scheme (const)
    Integral* const& singular_integral_pt() const
    {
      return Singular_integral_pt;
    }

    /// Access function: Pointer to boolean (true if the element is on the edge)
    bool& element_is_on_edge()
    {
      return Element_is_on_edge;
    }

    /// Access function: Pointer to boolean (true if the element is on the edge)
    /// (const)
    bool element_is_on_edge() const
    {
      return Element_is_on_edge;
    }

    ////////////////////////////////////////////////////////////////////////////
    ///////////// Start of doubly inherited member functions ///////////////////
    ////////////////////////////////////////////////////////////////////////////

    // The current element inherits from QElement and QHermiteElement which both
    // contain member functions of the same name. We are therefore required to
    // define these member functions instead of inheriting them. The following
    // member functions are all chosen to be inherited from QElement. The doubly
    // inherited member functions which don't copy the member functions from
    // parent classes are left outside this 'section' marked by 'Start of doubly
    // inherited member functions'

    /// Overload the template-free interface for the calculation of
    /// inverse jacobian matrix. This is a one-dimensional element, so
    /// use the 1D version.
    double invert_jacobian_mapping(
      const DenseMatrix<double>& jacobian,
      DenseMatrix<double>& inverse_jacobian) const override
    {
      return QElement<1, NNODE_1D>::invert_jacobian_mapping(jacobian,
                                                            inverse_jacobian);
    }

    /// Check whether the local coordinate are valid or not
    bool local_coord_is_valid(const Vector<double>& s) override
    {
      return QElement<1, NNODE_1D>::local_coord_is_valid(s);
    }

    /// Adjust local coordinates so that they're located inside
    /// the element
    void move_local_coord_back_into_element(Vector<double>& s) const override
    {
      QElement<1, NNODE_1D>::move_local_coord_back_into_element(s);
    }

    /// Get local coordinates of node j in the element; vector sets its own size
    void local_coordinate_of_node(const unsigned& j,
                                  Vector<double>& s) const override
    {
      QElement<1, NNODE_1D>::local_coordinate_of_node(j, s);
    }

    /// Get the local fraction of node j in the element
    void local_fraction_of_node(const unsigned& j,
                                Vector<double>& s_fraction) override
    {
      QElement<1, NNODE_1D>::local_fraction_of_node(j, s_fraction);
    }

    /// This function returns the local fraction of all nodes at the n-th
    /// position in a one dimensional expansion along the i-th local coordinate
    double local_one_d_fraction_of_node(const unsigned& n1d,
                                        const unsigned& i) override
    {
      return QElement<1, NNODE_1D>::local_one_d_fraction_of_node(n1d, i);
    }

    /// Calculate the geometric shape functions at local coordinate s
    void shape(const Vector<double>& s, Shape& psi) const override
    {
      QElement<1, NNODE_1D>::shape(s, psi);
    }

    /// Compute the geometric shape functions and derivatives w.r.t. local
    /// coordinates at local coordinate s
    void dshape_local(const Vector<double>& s,
                      Shape& psi,
                      DShape& dpsids) const override
    {
      QElement<1, NNODE_1D>::dshape_local(s, psi, dpsids);
    }

    /// Compute the geometric shape functions, derivatives and
    /// second derivatives w.r.t. local coordinates at local coordinate s
    /// d2psids(i,0) = \f$ d^2 \psi_j / d s^2 \f$
    void d2shape_local(const Vector<double>& s,
                       Shape& psi,
                       DShape& dpsids,
                       DShape& d2psids) const override
    {
      QElement<1, NNODE_1D>::d2shape_local(s, psi, dpsids, d2psids);
    }

    /// Number of nodes along each element edge
    unsigned nnode_1d() const override
    {
      return QElement<1, NNODE_1D>::nnode_1d();
    }

    /// Min. value of local coordinate
    double s_min() const override
    {
      return QElement<1, NNODE_1D>::s_min();
    }

    /// Max. value of local coordinate
    double s_max() const override
    {
      return QElement<1, NNODE_1D>::s_max();
    }

    /// Output
    void output(std::ostream& outfile) override
    {
      QElement<1, NNODE_1D>::output(outfile);
    }

    /// C-style output
    void output(FILE* file_pt) override
    {
      QElement<1, NNODE_1D>::output(file_pt);
    }

    /// C-style output at n_plot points
    void output(FILE* file_pt, const unsigned& n_plot) override
    {
      QElement<1, NNODE_1D>::output(file_pt, n_plot);
    }

    /// Get vector of local coordinates of plot point i (when plotting
    /// nplot points in each "coordinate direction).
    void get_s_plot(const unsigned& i,
                    const unsigned& nplot,
                    Vector<double>& s,
                    const bool& shifted_to_interior = false) const override
    {
      QElement<1, NNODE_1D>::get_s_plot(i, nplot, s, shifted_to_interior);
    }

    /// Return string for tecplot zone header (when plotting
    /// nplot points in each "coordinate direction)
    std::string tecplot_zone_string(const unsigned& nplot) const override
    {
      return QElement<1, NNODE_1D>::tecplot_zone_string(nplot);
    }

    /// Return total number of plot points (when plotting
    /// nplot points in each "coordinate direction)
    unsigned nplot_points(const unsigned& nplot) const override
    {
      return QElement<1, NNODE_1D>::nplot_points(nplot);
    }

    ////////////////////////////////////////////////////////////////////////////
    /////////////// End of doubly inherited member functions ///////////////////
    ////////////////////////////////////////////////////////////////////////////

  private:
    /// Store pointer to Eta
    double* Eta_pt;

    /// Store Poisson's ratio
    double* Nu_pt;

    /// Store the azimuthal wavemode
    unsigned* Wavemode_pt;

    /// Pointer to the pressure function
    AxisymFvKPressureFctPt Pressure_fct_pt;

    /// Pointer to the singular pressure function
    AxisymFvKSingularPressureFctPt Singular_pressure_fct_pt;

    /// Pointer to the forcing term
    AxisymFvKForcingFctPt Forcing_fct_pt;

    /// Default integration rule: Gaussian integration of the same order as the
    /// Lagrangian shape+basis functions
    /// Integration is exact for the shape function, u-basis and w-basis
    /// functions
    // hierher: Determine correct integration scheme order
    static Gauss<1, 4> Default_integration_scheme;

    /// The integration scheme dealing with the square root singularity
    static SingularIntegration<1, 3> Singular_integration_scheme;

    /// Pointer to the integration scheme that integrates square root
    /// singularities
    Integral* Singular_integral_pt;

    /// Boolean that returns true if the element is on the edge
    bool Element_is_on_edge;
  };

  //=======================================================================
  /// Assign the static Default_integration_scheme
  //=======================================================================
  // hierher: Determine the correct Gauss scheme order
  template<unsigned NNODE_1D>
  Gauss<1, 4> AxisymFvkElement<NNODE_1D>::Default_integration_scheme;

  //=======================================================================
  /// Assign the static Singular_integration_scheme
  //=======================================================================
  // The term to be integrated with the singular scheme will contain a fifth
  // order polynomial multiplied by the weight function (This is assuming the
  // pressure function varies quadratically (ignoring the singular weight
  // function)).
  // The singular Gauss integration scheme with 3 knots is sufficient to
  // integrate a polynomial of order 5.
  // hierher: Double check the order of the singular integration scheme is
  // sufficient
  template<unsigned NNODE_1D>
  SingularIntegration<1, 3>
    AxisymFvkElement<NNODE_1D>::Singular_integration_scheme;
} // namespace oomph

#endif
