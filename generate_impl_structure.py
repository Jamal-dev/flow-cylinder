import os

# List of member functions to create implementation files for
functions = [
    'set_runtime_parameters',
    'setup_system',
    'assemble_system_matrix',
    'assemble_system_rhs',
    'set_initial_bc',
    'set_newton_bc',
    'solve',
    'newton_iteration',
    'output_results'
]

# Create the 'step_flow_nse_impl' directory if it doesn't exist
os.makedirs('step_flow_nse_impl', exist_ok=True)

# Generate implementation files for each function
for func in functions:
    # Construct the filename
    filename = f'step_flow_nse_impl/{func}.cc'
    
    # Define the include guard macro
    guard_macro = func.upper() + '_CC'
    
    # Open the file for writing
    with open(filename, 'w') as f:
        # Write the file contents
        f.write(f'// {func}.cc\n\n')
        f.write(f'#ifndef {guard_macro}\n')
        f.write(f'#define {guard_macro}\n\n')
        f.write('#include "../step_flow_nse.h"\n\n')
        
        # Function signature
        f.write('template <int dim>\n')
        if func == 'output_results':
            # output_results has parameters
            f.write(f'void StepFlowNSE<dim>::{func}(const unsigned int refinement_cycle,\n')
            f.write(f'                              const BlockVector<double> solution_1) const\n')
        else:
            f.write(f'void StepFlowNSE<dim>::{func}()\n')
        f.write('{\n')
        f.write('    // Your implementation here\n')
        f.write('}\n\n')
        f.write(f'#endif // {guard_macro}\n')
