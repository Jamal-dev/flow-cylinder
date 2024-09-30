// output_results.cc

#ifndef OUTPUT_RESULTS_CC
#define OUTPUT_RESULTS_CC

#include "../step_flow_nse.h"

// This function is known from almost all other
// tutorial steps in deal.II.
template <int dim>
void StepFlowNSE<dim>::output_results(const unsigned int refinement_cycle,
                                      const BlockVector<double> output_vector_1) const
{

    std::vector<std::string> solution_names_1;
    solution_names_1.push_back("x_velo");
    solution_names_1.push_back("y_velo");
    solution_names_1.push_back("x_dis");
    solution_names_1.push_back("y_dis");
    solution_names_1.push_back("press");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation_1(dim + dim + 1, DataComponentInterpretation::component_is_scalar);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

    data_out.add_data_vector(output_vector_1, solution_names_1,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation_1);

    // INFO: Build higher order patches
    data_out.build_patches(2);

    std::ostringstream filename_vtk;
    std::ostringstream filename_ucd;

    std::cout << "------------------" << std::endl;
    std::cout << "Write solution" << std::endl;
    std::cout << "------------------" << std::endl;
    std::cout << std::endl;

    filename_vtk << filename_basis
                 << Utilities::int_to_string(refinement_cycle, 5)
                 << ".vtk";

    filename_ucd << filename_basis
                 << Utilities::int_to_string(refinement_cycle, 5)
                 << ".ucd";

    std::ofstream output_vtk(filename_vtk.str().c_str());
    data_out.write_vtk(output_vtk);
}

#endif // OUTPUT_RESULTS_CC
