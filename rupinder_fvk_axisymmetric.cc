// Generic oomph-lib routines
#include "generic.h"

// Include AxisymFvkEleme
#include "rupinder_include/axisym_fvk_elements.h"
#include "rupinder_include/singular_integration.cc"


using namespace oomph;

// Namespace for storing problem parameters
namespace GlobalParameters
{
  // Specify Eta
  double Eta = 1.0e5;

  // Specify Poisson's ratio
  double Nu = 0.5;

  // Specify the azimuthal wavemode
  unsigned Wavemode = 16;

  // Specify the number of elements being used
  unsigned NElements = 10;

  // Not the pressure function but simply sets the magnitude of the pressure
  // function
  double PressureMagnitude = 200.0;

  // This is the solution that have been chosen that obey the validation
  // boundary conditions
  void validation_solution(const Vector<double>& x, Vector<double>& u)
  {
    u.resize(2);
    u[0] = x[0] * exp(-(Nu + 1.0) * x[0]);
    u[1] = pow(x[0], 4) * pow(x[0] - 1.0, 3);
  }

  // We obtain the validation pressure function and validation forcing function
  // by plugging the validation solution into the FvK equations and seeing what
  // the pressure and forcing terms must be in order to satisfy the equations.

  // The pressure function in the validation case
  double validation_pressure_function(const double& r)
  {
    // Obtained from Mathematica
    return (2.0 * pow(-1.0 + r, 2) * pow(r, 2) * Eta * (1.0 + Nu) *
              (-16.0 + r * (67.0 + 4.0 * Nu +
                            r * (-67.0 - 11.0 * Nu + 7.0 * r * (1.0 + Nu)))) +
            exp(r * (1.0 + Nu)) *
              (-128 * (-1.0 + pow(Nu, 2)) +
               r * (1350 * (-1.0 + pow(Nu, 2)) +
                    r * (pow(4.0 - 7.0 * r, 2) * pow(-1.0 + r, 5) * pow(r, 6) *
                           (40.0 + r * (-155.0 + 133.0 * r)) * Eta -
                         3456.0 * (-1.0 + pow(Nu, 2)) +
                         2450.0 * r * (-1.0 + pow(Nu, 2)))))) /
           (2.0 * exp(r * (1.0 + Nu)) * (-1.0 + pow(Nu, 2)));
  }


  double validation_forcing_function(const double& r)
  {
    // Obtained from Mathematica
    return (pow(-1.0 + r, 3) * pow(r, 5) * (-4.0 + 7.0 * r) *
              (r * (107.0 + 7.0 * r * (-13.0 + Nu) - 11.0 * Nu) +
               4.0 * (-7.0 + Nu)) -
            (2.0 * (1.0 + Nu) * (-3.0 + r + r * Nu)) / exp(r * (1.0 + Nu))) /
           (2.0 * (-1.0 + pow(Nu, 2)));
  }

  // Specify the transverse pressure
  double pressure_function(const double& r)
  {
    return PressureMagnitude;
  }

  // Specify the singular transverse pressure
  double singular_pressure_function(const double& r)
  {
    return 0.0 / sqrt(1.0 - r * r);
  }
} // namespace GlobalParameters

// The mesh used in the AxisymFvKProblem, almost the same as the OneDMesh
// templated with AxisymFvKElements. Except we also define functions to count,
// refer to, or check Hermite nodes and our mesh is only between 0 and 1.
template<unsigned NNODE_1D>
class AxisymFvkMesh : public OneDMesh<AxisymFvkElement<NNODE_1D>>
{
public:
  /// The constructor for the AxisymFvkMesh calls the OneDMesh constructor, we
  /// also check that the nodes are equally spaced apart.
  // hierher: We do this because the local derivative of w instead of the global
  // derivative is stored at the Hermite nodes, meaning that if adjacent
  // elements are different sizes, the local derivative is not continuous. To
  // remove this problem, we store the global derivative at nodes instead.
  AxisymFvkMesh(const unsigned& n_element,
                TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
    : OneDMesh<AxisymFvkElement<NNODE_1D>>(n_element, 0.0, 1.0, time_stepper_pt)
  {
    // Get the number of nodes in the mesh
    unsigned n_node = this->nnode();

    // Check that the spacing between nodes is equal up to a tolerance
    double tolerance = 1.0e-8;

    // Storage for the distance between successive nodes
    double spacing_1 = 0.0;
    double spacing_2 = 0.0;

    // Loop over the nodes and check that the distance between nodes l and l+1
    // is the same as between nodes l+1 and l+2.
    for (unsigned l = 0; l < n_node - 2; l++)
    {
      spacing_1 = abs(this->node_pt(l)->x(0) - this->node_pt(l + 1)->x(0));
      spacing_2 = abs(this->node_pt(l + 1)->x(0) - this->node_pt(l + 2)->x(0));

      // Throw an error if the spacings are not equal within the tolerance
      if (abs(spacing_2 - spacing_1) > tolerance)
       {
        std::string error_message=
         "The nodes are not equally spaced in the mesh.";
        throw OomphLibError(
         error_message,
         OOMPH_EXCEPTION_LOCATION,
         OOMPH_CURRENT_FUNCTION);
      }
    }
  }

  /// Each element has two Hermite nodes, these nodes are shared between
  /// neighbouring elements, so there are N+1 Hermite nodes for N elements.
  unsigned nhermite_node() const
  {
    return Mesh::nelement() + 1;
  }

  /// Return pointer to the nth global node with Hermite data
  Node*& hermite_node_pt(const unsigned& n)
  {
    unsigned n_node = NNODE_1D;
    return Mesh::node_pt(n * (n_node - 1));
  }

  /// Return pointer to the nth global node with Hermite data (const version)
  Node* hermite_node_pt(const unsigned& n) const
  {
    unsigned n_node = NNODE_1D;
    return Mesh::node_pt(n * (n_node - 1));
  }

  /// Returns true if the chosen node is Hermite
  bool is_hermite_node(const unsigned& n) const
  {
    if (n % (NNODE_1D - 1) == 0)
    {
      return true;
    }

    return false;
  }
};


/// The problem class for the axisymmetric Foeppl-von Karman equations and the
/// eigenvalue problem created by adding a perturbation.
template<unsigned NNODE_1D, class EIGEN_SOLVER>
class AxisymFvkProblem : public Problem
{
public:
  /// Return the index at which the radial displacement is stored
  unsigned u_index_fvk() const
  {
    // Get the index from the element
    return dynamic_cast<AxisymFvkElement<NNODE_1D>*>(
             this->mesh_pt()->element_pt(0))
      ->u_index_fvk();
  }

  /// Return the index at which the radial displacement perturbation is stored
  unsigned u_pert_index_fvk() const
  {
    // Get the index from the element
    return dynamic_cast<AxisymFvkElement<NNODE_1D>*>(
             this->mesh_pt()->element_pt(0))
      ->u_pert_index_fvk();
  }

  /// Return the index at which the azimuthal displacement perturbation is
  /// stored
  unsigned u_theta_pert_index_fvk() const
  {
    // Get the index from the element
    return dynamic_cast<AxisymFvkElement<NNODE_1D>*>(
             this->mesh_pt()->element_pt(0))
      ->u_theta_pert_index_fvk();
  }

  /// Return the index at which the out of plane displacement is stored
  unsigned w_index_fvk() const
  {
    // Get the index from the element
    return dynamic_cast<AxisymFvkElement<NNODE_1D>*>(
             this->mesh_pt()->element_pt(0))
      ->w_index_fvk();
  }

  /// Return the index at which the local derivative of the out of plane
  /// displacement is stored
  unsigned dwds_index_fvk() const
  {
    // Get the index from the element
    return dynamic_cast<AxisymFvkElement<NNODE_1D>*>(
             this->mesh_pt()->element_pt(0))
      ->dwds_index_fvk();
  }

  /// Return the index at which the out of plane displacement perturbation is
  /// stored
  unsigned w_pert_index_fvk() const
  {
    // Get the index from the element
    return dynamic_cast<AxisymFvkElement<NNODE_1D>*>(
             this->mesh_pt()->element_pt(0))
      ->w_pert_index_fvk();
  }

  /// Return the index at which the local derivative of the out of plane
  /// displacement perturbation is stored
  unsigned dwds_pert_index_fvk() const
  {
    // Get the index from the element
    return dynamic_cast<AxisymFvkElement<NNODE_1D>*>(
             this->mesh_pt()->element_pt(0))
      ->dwds_pert_index_fvk();
  }

  // This deconstructor deletes the eigen solver and mesh that are stored on the
  // heap.
  ~AxisymFvkProblem()
  {
    delete mesh_pt();
    delete this->eigen_solver_pt();
  }

  /// The problem constructor
  AxisymFvkProblem(const unsigned& n_element)
  {
    // Create the eigen solver
    this->eigen_solver_pt() = new EIGEN_SOLVER;

    // Build the mesh and store pointer in Problem
    mesh_pt() = new AxisymFvkMesh<NNODE_1D>(n_element);

    // Set a flag signifying that this element is on the edge
    auto el_pt = dynamic_cast<AxisymFvkElement<NNODE_1D>*>(
      mesh_pt()->element_pt(n_element - 1));
    el_pt->element_is_on_edge() = true;

    // Set external data for all the elements
    for (unsigned i = 0; i < n_element; i++)
    {
      // Get a pointer to the element
      el_pt =
        dynamic_cast<AxisymFvkElement<NNODE_1D>*>(mesh_pt()->element_pt(i));

      // Set Eta in the element
      el_pt->eta_pt() = &GlobalParameters::Eta;

      // Set Poisson's ratio in the element
      el_pt->nu_pt() = &GlobalParameters::Nu;

      // Set the pressure in the element
      el_pt->pressure_fct_pt() = &GlobalParameters::pressure_function;

      // Set the wavemode in the element
      el_pt->wavemode_pt() = &GlobalParameters::Wavemode;

      // Set the singular pressure in the element
      el_pt->singular_pressure_fct_pt() =
        &GlobalParameters::singular_pressure_function;

      if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
      {
        // Set the pressure in the element
        el_pt->pressure_fct_pt() =
          &GlobalParameters::validation_pressure_function;

        // Set the forcing in the element
        el_pt->forcing_fct_pt() =
          &GlobalParameters::validation_forcing_function;

        // Set the singular pressure in the element
        el_pt->singular_pressure_fct_pt() = 0;
      }
    }

    // The problem initially requires solving the base axisymmetric equations so
    // we pin the dofs in the perturbation problem and set the boundary
    // conditions
    pin_pert_solution();
    set_boundary_conditions();

    // Setup equation numbering scheme
    assign_eqn_numbers();
  }

  // Set the boundary conditions
  void set_boundary_conditions()
  {
    // Get the boundary nodes
    auto centre_boundary_node = mesh_pt()->boundary_node_pt(0, 0);
    auto outer_boundary_node = mesh_pt()->boundary_node_pt(1, 0);

    // Base solution boundary conditions
    //----------------------------------

    // Centre boundary conditions
    centre_boundary_node->pin(u_index_fvk());
    centre_boundary_node->pin(dwds_index_fvk());

    // Outer boundary conditions
    outer_boundary_node->pin(w_index_fvk());
    outer_boundary_node->pin(dwds_index_fvk());

    // Perturbation boundary conditions
    //---------------------------------

    // Centre boundary conditions
    centre_boundary_node->pin(u_pert_index_fvk());
    centre_boundary_node->pin(dwds_pert_index_fvk());

    // Outer boundary conditions
    // outer_boundary_node->pin(u_theta_pert_index_fvk());
    outer_boundary_node->pin(w_pert_index_fvk());
    outer_boundary_node->pin(dwds_pert_index_fvk());
  }

  // Pin the base solution dofs at every node
  void pin_base_solution()
  {
    // Loop over each node and pin the base solution data
    unsigned n_node = this->mesh_pt()->nnode();
    auto mesh_pt = dynamic_cast<AxisymFvkMesh<NNODE_1D>*>(this->mesh_pt());
    for (unsigned l = 0; l < n_node; l++)
    {
      // Pin the radial displacement
      mesh_pt->node_pt(l)->pin(u_index_fvk());

      // Check if the node has Hermite data to pin
      if (mesh_pt->is_hermite_node(l))
      {
        // Pin the transverse displacement
        mesh_pt->node_pt(l)->pin(w_index_fvk());
        // Pin the transverse displacement local derivative
        mesh_pt->node_pt(l)->pin(dwds_index_fvk());
      }
    }
  }

  // Unpin the base solution dofs at every node
  void unpin_base_solution()
  {
    // Loop over each node and unpin the base solution data
    unsigned n_node = this->mesh_pt()->nnode();
    auto mesh_pt = dynamic_cast<AxisymFvkMesh<NNODE_1D>*>(this->mesh_pt());
    for (unsigned l = 0; l < n_node; l++)
    {
      // Unpin the radial displacement
      mesh_pt->node_pt(l)->unpin(u_index_fvk());

      // Check if the node has Hermite data to unpin
      if (mesh_pt->is_hermite_node(l))
      {
        // Unpin transverse displacement
        mesh_pt->node_pt(l)->unpin(w_index_fvk());
        // Unpin transverse displacement local derivative
        mesh_pt->node_pt(l)->unpin(dwds_index_fvk());
      }
    }
  }

  // Pin the perturbation dofs at every node
  void pin_pert_solution()
  {
    // Loop over each node and pin the perturbation solution data
    unsigned n_node = this->mesh_pt()->nnode();
    auto mesh_pt = dynamic_cast<AxisymFvkMesh<NNODE_1D>*>(this->mesh_pt());
    for (unsigned l = 0; l < n_node; l++)
    {
      // Pin the radial displacement
      mesh_pt->node_pt(l)->pin(u_pert_index_fvk());
      // Pin the azimuthal displacement
      mesh_pt->node_pt(l)->pin(u_theta_pert_index_fvk());

      // Check if the node has Hermite data to pin
      if (mesh_pt->is_hermite_node(l))
      {
        // Pin the transverse displacement
        mesh_pt->node_pt(l)->pin(w_pert_index_fvk());
        // Pin the transverse displacement derivative
        mesh_pt->node_pt(l)->pin(dwds_pert_index_fvk());
      }
    }
  }

  // Unpin the perturbation dofs at every node
  void unpin_pert_solution()
  {
    // Loop over each node and unpin the perturbation solution data
    unsigned n_node = this->mesh_pt()->nnode();
    auto mesh_pt = dynamic_cast<AxisymFvkMesh<NNODE_1D>*>(this->mesh_pt());
    for (unsigned l = 0; l < n_node; l++)
    {
      // Unpin the radial displacement
      mesh_pt->node_pt(l)->unpin(u_pert_index_fvk());
      // Unpin the azimuthal displacement
      mesh_pt->node_pt(l)->unpin(u_theta_pert_index_fvk());

      // Check if the node has Hermite data to unpin
      if (mesh_pt->is_hermite_node(l))
      {
        // Unpin the transverse displacement
        mesh_pt->node_pt(l)->unpin(w_pert_index_fvk());
        // Unpin the transverse displacement derivative
        mesh_pt->node_pt(l)->unpin(dwds_pert_index_fvk());
      }
    }
  }

  /// Document the solution
  void doc_solution(std::string filename, const unsigned& npts) const
  {
    // Output the computed solution
    std::ofstream outfile(filename);
    mesh_pt()->output(outfile, npts);
    outfile.close();

    // If we run the validation case, output the exact solution
    if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
    {
      // Output exact solution
      outfile.open("RESLT/exact_soln.dat");
      mesh_pt()->output_fct(
        outfile, npts, GlobalParameters::validation_solution);
      outfile.close();
    }
  }

  /// Integrate the pressure over the domain
  double integrate_pressure()
  {
    // Storage for the answer
    double result = 0.0;

    // Loop over the elements in the mesh
    unsigned n_element = this->mesh_pt()->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      // Get the pressure integrals over each element
      result += dynamic_cast<AxisymFvkElement<NNODE_1D>*>(
                  this->mesh_pt()->element_pt(e))
                  ->integrate_pressure();
    }

    return result;
  }

  // This struct stores associated eigenvalues and vectors together, useful so
  // that we can sort a vector of EigenData by an element of the struct.
  struct EigenData
  {
    std::complex<double> eigenvalue;
    DoubleVector eigenvector_real;
    DoubleVector eigenvector_imag;
  };

  /// Solve the eigenproblem for at least n_eval eigenvalues/vectors, sort the
  /// eigenvalues/vectors by the real part of the eigenvalues, keep the first
  /// n_eval eigenvectors/values and assign the first eigenvector to the dofs.
  Vector<std::complex<double>> eigensolve(const unsigned& n_eval)
  {
    // Storage for the eigenvalues
    Vector<std::complex<double>> eigenvalues;

    // Storage for the eigenvectors
    Vector<DoubleVector> eigenvectors_real;
    Vector<DoubleVector> eigenvectors_imag;

    // Solve the eigenvalue problem
    this->solve_eigenproblem(
      n_eval, eigenvalues, eigenvectors_real, eigenvectors_imag);

    // Store the eigenvalues and eigenvectors in a vector of EigenData structs
    // for easy sorting
    Vector<EigenData> eigen_data(eigenvalues.size());
    for (unsigned i = 0; i < eigenvalues.size(); i++)
    {
      eigen_data[i].eigenvalue = eigenvalues[i];
      eigen_data[i].eigenvector_real = eigenvectors_real[i];
      eigen_data[i].eigenvector_imag = eigenvectors_imag[i];
    }

    // Sort the eigen data by the real component of the eigenvalues
    // The third argument is a predicate that defines how an element of
    // eigen_data is 'bigger' than another
    sort(eigen_data.begin(),
         eigen_data.end(),
         [](const EigenData& a, const EigenData& b) -> bool
         { return abs(a.eigenvalue.real()) < abs(b.eigenvalue.real()); });

    // Unpack the sorted vector of eigendata back into the eigenvalue and
    // eigenvector of eigenvalues and eigenvectors
    for (unsigned i = 0; i < n_eval; i++)
    {
      eigenvalues[i] = eigen_data[i].eigenvalue;
      eigenvectors_real[i] = eigen_data[i].eigenvector_real;
      eigenvectors_imag[i] = eigen_data[i].eigenvector_imag;
    }

    // Truncate the additional eigenvalues and vectors
    eigenvalues.resize(n_eval);
    eigenvectors_real.resize(n_eval);
    eigenvectors_imag.resize(n_eval);

    // Assign the first eigenvector to dofs so that the displacement can be seen
    assign_eigenvector_to_dofs(eigenvectors_real[0]);

    return eigenvalues;
  }

  // Solve the base axisymmetric problem then return the smallest eigenvalue of
  // the Jacobian in the eigenproblem
  double solve_base_problem_and_find_eigenvalue()
  {
    // Pin variables in the perturbation problem, unpin variables in the base
    // problem and apply boundary conditions
    unpin_base_solution();
    pin_pert_solution();
    set_boundary_conditions();
    assign_eqn_numbers();

    // Solve the axisymmetric base problem
    newton_solve();

    // Pin variables in the base problem, unpin variables in the
    // perturbation problem and apply boundary conditions
    pin_base_solution();
    unpin_pert_solution();
    set_boundary_conditions();
    assign_eqn_numbers();

    // Find the smallest eigenvalue of the Jacobian in the eigenproblem
    return eigensolve(1)[0].real();
  }
};


int main(int argc, char** argv)
{
  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  // Choose between pure bending or fvk model
  CommandLineArgs::specify_command_line_flag("--validation");

  // Set the wavemode from the command line
  CommandLineArgs::specify_command_line_flag("--wavemode",
                                             &GlobalParameters::Wavemode);

  // Set the number of elements from the command line
  CommandLineArgs::specify_command_line_flag("--n_elements",
                                             &GlobalParameters::NElements);

  // Parse command line
  CommandLineArgs::parse_and_assign();

  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();

  // If running the validation case, run this section of code otherwise, run
  // everything else.
  if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
    // Assign storage for the number of elements
    unsigned n_element = 0;

    // Set the number of output points per element
    unsigned npts = 5;

    // Run the validation case for the following varying number of elements
    std::vector<unsigned> n_element_list = {2, 5, 10, 100};

    for (unsigned i = 0; i < n_element_list.size(); i++)
    {
      n_element = n_element_list[i];

      // Output the results to the file of this name
      std::string filename = "RESLT/validation_solution_" +
                             std::to_string(n_element_list[i]) +
                             "_elements.dat";

      // Create the problem, solve the axisymmetric problem and document it.
      AxisymFvkProblem<3, LAPACK_QZ> problem(n_element);

      problem.max_residuals() = 1.0e100;
      problem.max_newton_iterations() = 100;

      problem.newton_solve();

      problem.doc_solution(filename, npts);
    }
  }
  else
  {
    ////////////////////////////////////////////////////////////////////////////
    // Start of setting parameters
    ////////////////////////////////////////////////////////////////////////////
    // Set the number of output points per element
    unsigned npts = 5;

    // Solution doc file
    std::string doc_solution_filename = "RESLT/soln.dat";
    std::string doc_critical_pressure_filename = "RESLT/critical_pressure.dat";

    // How small the eigenvalue must be to accept the critical pressure value
    double eigenvalue_tolerance = 1.0e-6;

    // The nondimensional measure aspect ratio
    GlobalParameters::Eta = 1.0e5;

    // Magnitude of the transverse pressure
    GlobalParameters::PressureMagnitude = 200.0;

    // Set the range of pressures that we guess within
    double minimum_pressure_guess = 0.0;
    double maximum_pressure_guess = 1000.0;

    ////////////////////////////////////////////////////////////////////////////
    // End of setting parameters
    ////////////////////////////////////////////////////////////////////////////

    // Create the problem
    AxisymFvkProblem<3, LAPACK_QZ> problem(GlobalParameters::NElements);

    // Allow the Newton solver to have a bad guess, use more steps and stop with
    // a larger residual
    problem.max_residuals() = 1.0e100;
    problem.max_newton_iterations() = 1000;
    problem.newton_solver_tolerance() = 1.0e-7;

    // Temp values
    double min_guess_eigenvalue = 0.0;
    double max_guess_eigenvalue = 0.0;
    double current_guess_eigenvalue = 0.0;

    // Flag to say if the eigenvalue increases when the pressure increases
    bool eigenvalue_is_increasing = false;

    // First, check that the smallest eigenvalue when using the minimum and
    // maximum pressures have opposite signs and therefore there may be a root
    // to find within the range provided
    GlobalParameters::PressureMagnitude = minimum_pressure_guess;
    min_guess_eigenvalue = problem.solve_base_problem_and_find_eigenvalue();

    GlobalParameters::PressureMagnitude = maximum_pressure_guess;
    max_guess_eigenvalue = problem.solve_base_problem_and_find_eigenvalue();

    if (max_guess_eigenvalue * min_guess_eigenvalue > 0.0)
    {
      std::cout << "Invalid range to search" << std::endl;
      return 1;
    }

    // If the eigenvalue from the maximum pressure is positive then the
    // eigenvalues must be increasing
    if (max_guess_eigenvalue > 0)
    {
      eigenvalue_is_increasing = true;
    }
    else
    {
      eigenvalue_is_increasing = false;
    }

    // Do a binary search to find the pressure corresponding to a zero
    // eigenvalue
    current_guess_eigenvalue = 1.0;
    while (abs(current_guess_eigenvalue) > eigenvalue_tolerance)
    {
      // Choose the average pressure as a guess
      GlobalParameters::PressureMagnitude =
        (minimum_pressure_guess + maximum_pressure_guess) / 2.0;

      // Compute the smallest eigenvalue
      current_guess_eigenvalue =
        problem.solve_base_problem_and_find_eigenvalue();

      // Set the average pressure as the new min or max pressure guess such that
      // the range of pressures contains the root
      if ((current_guess_eigenvalue > 0.0) == (eigenvalue_is_increasing))
      {
        maximum_pressure_guess = GlobalParameters::PressureMagnitude;
      }
      else
      {
        minimum_pressure_guess = GlobalParameters::PressureMagnitude;
      }
    }

    // Document the solution
    problem.doc_solution(doc_solution_filename, npts);

    // Document the critical pressure
    std::ofstream output_file(doc_critical_pressure_filename);

    output_file.precision(16);
    output_file << GlobalParameters::NElements << " "
                << GlobalParameters::Wavemode << " "
                << GlobalParameters::PressureMagnitude << " "
                << current_guess_eigenvalue << std::endl;

    output_file.close();
  }

  return 0;
}
