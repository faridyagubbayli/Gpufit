#ifndef GPUFIT_SIMPLIFIED_GAMMA_VARIATE_CUH_INCLUDED
#define GPUFIT_SIMPLIFIED_GAMMA_VARIATE_CUH_INCLUDED

#include <stdio.h>

__device__ void calculate_simplified_gamma_variate(
    REAL const * parameters,
    int const n_fits,
    int const n_points,
    REAL * value,
    REAL * derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char * user_info,
    std::size_t const user_info_size)
{
    // indices

    REAL * user_info_float = (REAL*)user_info;
    REAL x = 0;
    if (!user_info_float)
    {
        x = point_index;
    }
    else if (user_info_size / sizeof(REAL) == n_points)
    {
        x = user_info_float[point_index];
    }
    else if (user_info_size / sizeof(REAL) > n_points)
    {
        int const chunk_begin = chunk_index * n_fits * n_points;
        int const fit_begin = fit_index * n_points;
        x = user_info_float[chunk_begin + fit_begin + point_index];
    }

    // parameters

    REAL const * p = parameters;
	// p[0] => t0
	// p[1] => t_max
	// p[2] => y_max
	// p[3] => alpha

	// value

	if (x - p[0] >= 0)
	{
		REAL t_prime = (x - p[0]) / (p[1] - p[0]);
		value[point_index] = exp(log(p[2]) + p[3] * (1 + log(t_prime + 0.001) - t_prime));
	}
	else
	{
		value[point_index] = 0;
	}

	// derivative

	REAL* current_derivative = derivative + point_index;

	if (x - p[0] >= 0)
	{
		REAL t_prime = (x - p[0]) / (p[1] - p[0]);

		current_derivative[0 * n_points] = 0;
		current_derivative[1 * n_points] = 0;
		current_derivative[2 * n_points] = 0;
		current_derivative[3 * n_points] = value[point_index] * (1 + log(t_prime + 0.001) - t_prime);
	}
	else
	{
		current_derivative[0 * n_points] = 0.1;
		current_derivative[1 * n_points] = 0;
		current_derivative[2 * n_points] = 0;
		current_derivative[3 * n_points] = 0;
	}




    
    /*// value

	if (x - p[3] >= 0)
	{
		value[point_index] = p[0] * powf(x - p[3], p[1]) * exp(-(x - p[3]) / p[2]);
	}
	else
	{
		value[point_index] = 0;
	}

    // derivative

    REAL * current_derivative = derivative + point_index;

	if (x - p[3] >= 0)
	{
		current_derivative[0 * n_points] = powf(x - p[3], p[1]) * exp(-(x - p[3]) / p[2]);
		current_derivative[1 * n_points] = p[0] * powf(x - p[3], p[1]) * exp(-(x - p[3]) / p[2]) * log(x - p[3] + 0.001);
		current_derivative[2 * n_points] = (p[0] * exp(-(x - p[3]) / p[2]) * powf(x - p[3], p[1]) * (x - p[3])) / (p[2] * p[2]);
		current_derivative[3 * n_points] = ((p[0] * exp(-(x - p[3]) / p[2]) * powf(x - p[3], p[1])) / p[2]) - (p[0] * p[1] * exp(-(x - p[3]) / p[2]) * powf(x - p[3], p[1] - 1));
	}
	else
	{
		current_derivative[0 * n_points] = 0.1;
		current_derivative[1 * n_points] = 0;
		current_derivative[2 * n_points] = 0;
		current_derivative[3 * n_points] = 0.5;
	}*/

	
	/*printf("X: %f\n", x);
	printf("Out val: %f\n", value[point_index]);
	printf("Grad A: %f\n", current_derivative[0 * n_points]);
	printf("Grad alpha: %f\n", current_derivative[1 * n_points]);
	printf("Grad beta: %f\n", current_derivative[2 * n_points]);
	printf("Grad t0: %f\n", current_derivative[3 * n_points]);*/

}

#endif
