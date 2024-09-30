// set_runtime_parameters.cc

#ifndef SET_RUNTIME_PARAMETERS_CC
#define SET_RUNTIME_PARAMETERS_CC

#include "../step_flow_nse.h"

template <int dim>
void StepFlowNSE<dim>::set_runtime_parameters()
{
    std::string grid_name;

    // Fluid parameters for all test cases
    if (test_case == "Channel_NSE")
    {
        density_fluid = 1.0;
        kinematic_viscosity = 1.0e-2;

        // Timestepping schemes
        // BE, CN, CN_shifted
        time_stepping_scheme = "CN";

        // Timestep size:
        timestep = 1.0;

        // Max number of time steps
        max_no_timesteps = 1000;

        grid_name = "geometry_files/rectangle_mandel.inp";
        // grid_name  = "unit_square_1.inp";

        GridIn<dim> grid_in;
        grid_in.attach_triangulation(triangulation);
        std::ifstream input_file(grid_name.c_str());
        Assert(dim == 2, ExcInternalError());
        grid_in.read_ucd(input_file);

        triangulation.refine_global(3);

        filename_basis = "output_files/solution_channel_NSE_";
    }
    else if (test_case == "2D-1_NSE")
    {
        // Fluid parameters for all test cases
        density_fluid = 1.0;
        kinematic_viscosity = 1.0e-3;

        // Timestepping schemes
        // BE, CN, CN_shifted
        time_stepping_scheme = "BE";

        // Timestep size:
        timestep = 1.0;

        // Max number of time steps
        max_no_timesteps = 25;

        grid_name = "geometry_files/nsbench4_original.inp";

        GridIn<dim> grid_in;
        grid_in.attach_triangulation(triangulation);
        std::ifstream input_file(grid_name.c_str());
        Assert(dim == 2, ExcInternalError());
        grid_in.read_ucd(input_file);

        // TODO: only when ns bench
        Point<dim> p(0.2, 0.2);
        // double radius = 0.05;
        const SphericalManifold<dim> boundary(p);
        triangulation.set_all_manifold_ids_on_boundary(80, 8);
        // triangulation.set_all_manifold_ids_on_boundary(81,9);
        triangulation.set_manifold(8, boundary);
        // triangulation.set_manifold (9, boundary);

        triangulation.refine_global(2);

        filename_basis = "output_files/solution_2D_1_NSE_";
    }
    else if (test_case == "Cavity_IV")
    {
        // Fluid parameters for all test cases
        density_fluid = 1.0;
        kinematic_viscosity = 1.0e-3;

        // Timestepping schemes
        // BE, CN, CN_shifted
        time_stepping_scheme = "BE";

        // Timestep size:
        timestep = 1.0;

        // Max number of time steps
        max_no_timesteps = 10000; // 500;

        grid_name = "geometry_files/cavity.inp";
        // grid_name  = "unit_square_1.inp";

        GridIn<dim> grid_in;
        grid_in.attach_triangulation(triangulation);
        std::ifstream input_file(grid_name.c_str());
        Assert(dim == 2, ExcInternalError());
        grid_in.read_ucd(input_file);

        triangulation.refine_global(3);

        filename_basis = "output_files/solution_cavity_";
    }
    else if (test_case == "L_shaped_IV")
    {
        // Fluid parameters for all test cases
        density_fluid = 1.0;
        kinematic_viscosity = 1.0e-3;

        // Timestepping schemes
        // BE, CN, CN_shifted
        time_stepping_scheme = "BE";

        // Timestep size:
        timestep = 0.5;

        // Max number of time steps
        max_no_timesteps = 500;

        grid_name = "geometry_files/l_shaped.inp";
        // grid_name  = "unit_square_1.inp";

        GridIn<dim> grid_in;
        grid_in.attach_triangulation(triangulation);
        std::ifstream input_file(grid_name.c_str());
        Assert(dim == 2, ExcInternalError());
        grid_in.read_ucd(input_file);

        triangulation.refine_global(3);

        filename_basis = "output_files/solution_l_shaped_IV_";
    }
    else if (test_case == "Venturi_IV")
    {
        // Fluid parameters for all test cases
        density_fluid = 1.0;
        kinematic_viscosity = 1.0e-3;

        // Timestepping schemes
        // BE, CN, CN_shifted
        time_stepping_scheme = "CN";

        // Timestep size:
        timestep = 0.5;

        // Max number of time steps
        max_no_timesteps = 100;

        grid_name = "geometry_files/venturi.inp";
        // grid_name  = "unit_square_1.inp";

        GridIn<dim> grid_in;
        grid_in.attach_triangulation(triangulation);
        std::ifstream input_file(grid_name.c_str());
        Assert(dim == 2, ExcInternalError());
        grid_in.read_ucd(input_file);

        triangulation.refine_global(4);

        filename_basis = "output_files/solution_venturi_IV_";
    }
    else
    {
        std::cout << "Aborting in set runtime parameters." << std::endl;
        abort();
    }

    // Right hand side flow equation, e.g., gravity
    force_fluid_x = 0.0;
    force_fluid_y = 0.0;

    // Traction - not in use here
    traction_x = 0.0;
    traction_y = 0.0;

    // A variable to count the number of time steps
    timestep_number = 0;

    // Counts total time
    time = 0;

    // Here, we choose a time-stepping scheme that
    // is based on finite differences:
    // BE         = backward Euler scheme
    // CN         = Crank-Nicolson scheme
    // CN_shifted = time-shifted Crank-Nicolson scheme
    // For further properties of these schemes,
    // we refer to standard literature.
    if (time_stepping_scheme == "BE")
        theta = 1.0;
    else if (time_stepping_scheme == "CN")
        theta = 0.5;
    else if (time_stepping_scheme == "CN_shifted")
        theta = 0.5 + 0.1; // timestep;
    else
        std::cout << "No such timestepping scheme" << std::endl;
}

#endif // SET_RUNTIME_PARAMETERS_CC
