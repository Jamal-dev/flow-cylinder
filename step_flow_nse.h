// step_flow_nse.h

#ifndef STEP_FLOW_NSE_H
#define STEP_FLOW_NSE_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
//#include <deal.II/grid/tria_boundary_lib.h> // deprecated
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
//#include <deal.II/lac/constraint_matrix.h>  // deprecated
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>

#include <fstream>
#include <sstream>

using namespace dealii;

#include "tensors.h"
#include "initial_values.h"
#include "boundary_values_step.h"


template <int dim>
class StepFlowNSE
{
public:
  StepFlowNSE(const unsigned int degree);
  ~StepFlowNSE();
  void run();

private:
  void set_runtime_parameters();
  void setup_system();
  void assemble_system_matrix();
  void assemble_system_rhs();
  void set_initial_bc();
  void set_newton_bc();
  void solve();
  void newton_iteration();
  void output_results(const unsigned int refinement_cycle,
                      const BlockVector<double> solution_1) const;
  double compute_point_value(Point<dim> p,
                             const unsigned int component) const;
  void compute_drag_lift_tensor();
  void compute_functional_values();

  const unsigned int degree;

  Triangulation<dim> triangulation;
  FESystem<dim> fe;
  DoFHandler<dim> dof_handler;

  AffineConstraints<double> constraints;

  BlockSparsityPattern sparsity_pattern;
  BlockSparseMatrix<double> system_matrix;

  BlockVector<double> solution, newton_update, old_timestep_solution, old_old_timestep_solution;
  BlockVector<double> system_rhs;

  TimerOutput timer;

  // Global variables for timestepping scheme
  unsigned int timestep_number;
  unsigned int max_no_timesteps;
  double timestep, theta, time;
  std::string time_stepping_scheme;

  double force_fluid_x, force_fluid_y;

  double volume_source, traction_x, traction_y, traction_x_biot, traction_y_biot;

  // Fluid parameters
  double density_fluid, kinematic_viscosity;

  SparseDirectUMFPACK A_direct;

  std::string test_case;

  std::string filename_basis;
};

template <int dim>
StepFlowNSE<dim>::StepFlowNSE(const unsigned int degree)
    : degree(degree),
      triangulation(Triangulation<dim>::maximum_smoothing),
      fe(FE_Q<dim>(degree + 1), dim, // velocity (vector-valued)
         FE_Q<dim>(degree + 1), dim, // displacements (vector-valued)
         FE_Q<dim>(degree), 1),      // pressure (scalar-valued)
      dof_handler(triangulation),
      timer(std::cout, TimerOutput::summary, TimerOutput::cpu_times)
{}

// Destructor
template <int dim>
StepFlowNSE<dim>::~StepFlowNSE()
{}

// Include implementation files
#include "step_flow_nse_impl/set_runtime_parameters.cc"
#include "step_flow_nse_impl/setup_system.cc"
#include "step_flow_nse_impl/assemble_system_matrix.cc"
#include "step_flow_nse_impl/assemble_system_rhs.cc"
#include "step_flow_nse_impl/set_initial_bc.cc"
#include "step_flow_nse_impl/set_newton_bc.cc"
#include "step_flow_nse_impl/solve.cc"
#include "step_flow_nse_impl/newton_iteration.cc"
#include "step_flow_nse_impl/output_results.cc"
#include "step_flow_nse_impl/compute_point_value.cc"
#include "step_flow_nse_impl/compute_drag_lift_tensor.cc"
#include "step_flow_nse_impl/compute_functional_values.cc"
#include "step_flow_nse_impl/run.cc"

#endif // STEP_FLOW_NSE_H
