// compute_point_value.cc

#ifndef COMPOUTE_POINT_VALUE_CC
#define COMPOUTE_POINT_VALUE_CC

#include "../step_flow_nse.h"

template <int dim>
double StepFlowNSE<dim>::compute_point_value(Point<dim> p,
                                             const unsigned int component) const
{

    Vector<double> tmp_vector(dim + dim + 1);
    VectorTools::point_value(dof_handler,
                             solution,
                             p,
                             tmp_vector);

    return tmp_vector(component);
}

#endif // COMPOUTE_POINT_VALUE_CC
