#include "../gpufit.h"
#include <iostream>
#include <vector>

#include <array>
#include <cmath>
#include <time.h>

template<std::size_t n_points, std::size_t n_parameters>
void generate_gauss_1d(
    std::array< REAL, n_points >& values,
    std::array< REAL, n_parameters > const & parameters,
	std::size_t const n_fits)
{
    REAL const A = parameters[ 0 ];
    REAL const alpha = parameters[ 1 ];
    REAL const beta = parameters[ 2 ];
    REAL const t0 = parameters[ 3 ];

	std::size_t points_per_fit = n_points / n_fits;

	for (int f = 0; f < n_fits; f++)
	{
		for (int point_index = 0; point_index < points_per_fit; point_index++)
		{
			REAL const t = point_index;
			//printf("%d, %d, %d\n", (f * points_per_fit) + point_index, f, point_index);
			if (t - t0 >= 0)
			{
				values[(f * points_per_fit) + point_index] = A * powf(t - t0, alpha) * exp(-(t - t0) / beta);
			}
			else
			{
				values[(f * points_per_fit) + point_index] = 0;
			}
		}
	}
    
}

std::array< REAL, 520 * 520 * 36 > data{};
std::array< REAL, 520 * 520 * 4 > initial_parameters{};
std::array< REAL, 520 * 520 * 4 > output_parameters;

std::array< int, 520 * 520 > output_states;
std::array< REAL, 520 * 520 > output_chi_square;
std::array< int, 520 * 520 > output_n_iterations;

void simplified_gamma_variate_fit()
{
    /*
    Performs a single fit using the SIMPLIFIED_GAMMA_VARIATE model.
    - Doesn't use user_info or weights.
    - No noise is added.
    - Checks fitted parameters equalling the true parameters.
    */

	std::size_t const n_fits{ 520 * 520 };
    std::size_t const n_points{ 36 };
    std::size_t const n_parameters{ 4 };

    std::array< REAL, n_parameters > const true_parameters{ { 2.0f, 1.8f, 1.2f, 3.0f } };

    
    generate_gauss_1d(data, true_parameters, n_fits);

	/*printf("Data points:\n");
	for (int i = 0; i < n_points; i++)
	{
		printf("%d >>\t\t %.3f\n", i, data[i]);
	}*/

    
	for (int i = 0; i < n_fits; i++)
	{
		initial_parameters[i * n_parameters + 0] = 3.0f;
		initial_parameters[i * n_parameters + 1] = 5.16f;
		initial_parameters[i * n_parameters + 2] = 1.3154f;
		initial_parameters[i * n_parameters + 3] = 0.1f;
	}

    REAL tolerance{ 0.0001f };

    int max_n_iterations{ 2000 };

    std::array< int, n_parameters > parameters_to_fit{ { 0, 0, 0, 1} };

	clock_t tStart = clock();

    int const status
        = gpufit
        (
            n_fits,
            n_points,
            data.data(),
            0,
			SIMPLIFIED_GAMMA_VARIATE,
            initial_parameters.data(),
            tolerance,
            max_n_iterations,
            parameters_to_fit.data(),
            LSE,
            0,
            0,
			output_parameters.data(),
			output_states.data(),
			output_chi_square.data(),
			output_n_iterations.data()
        );
	clock_t tEnd = clock();

	printf("Params:\n");
	for (int i = 0; i < 40; i++)
	{
		printf("%.3f\n", output_parameters[i * n_parameters + 3]);
	}
	//printf("Output Chi Square %.4f:\n", output_chi_square);

	printf("Time taken: %.2fs\n", (double)(tEnd - tStart) / CLOCKS_PER_SEC);
}

int main(int argc, char* argv[])
{
	simplified_gamma_variate_fit();

	std::cout << std::endl << "Simplified Gamma-Variate Example completed!" << std::endl;
	std::cout << "Press ENTER to exit" << std::endl;
	std::getchar();

	return 0;
}
