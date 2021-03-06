
package algorithms;

public class LevenbergMarquardt
{
	/**
	 *	Constants
	 */
	public static final int GAUSSIAN_1D_PARAMETERS = 4;
	public static final int PARAM_BACKGROUND = 0;
	public static final int PARAM_AMPLITUDE = 1;
	public static final int PARAM_MU = 2;
	public static final int PARAM_SIGMA = 3;
	
	public static final double LAMBDA_FACTOR = 10.0f;
	public static final double LAMBDA_FACTOR_INV = 0.1f;
	
	public static final double MIN_EPSILON = 1e-30f;
	
	/**
	 *	Members
	 */
	//public double[][] 
	
	/**
	 *	Constructor
	 */
	public LevenbergMarquardt()
	{
		/* empty */
	}
	
	/**
	 *
	 */
	public double gaussian1D(double x, double[] parameters)
	{
		return parameters[PARAM_BACKGROUND] + parameters[PARAM_AMPLITUDE] * Math.exp(-0.5f * ((x-parameters[PARAM_MU])*(x-parameters[PARAM_MU])/(parameters[PARAM_SIGMA]*parameters[PARAM_SIGMA])));
	}
	
	/*public double gaussian2Dn(double x, double y, double[] parameters)
	{
		return (1/(parameters[PARAM_SIGMA]*Math.sqrt(2*Math.PI))) * Math.exp(-0.5f * ((x-parameters[PARAM_MU])*(x-parameters[PARAM_MU])/(parameters[PARAM_SIGMA]*parameters[PARAM_SIGMA])));
	}*/

	/**
	 *
	 */
	public double calculateChi2(double[] x_data, double[] y_data, int data_points, double[] parameters)
	{
		double chi2 = 0.0;
		
		// loop over all pixels in image to sum the product of difference between image and model
		for(int p = 0; p < data_points; ++p)
		{
			// calculate sum of squared difference
			double dy = y_data[p] - gaussian1D(x_data[p], parameters); // RSLV: inline Gaussian1D function?
			// NOTE: no normalisation factor necessary because it is constant between iterations and would only add to the computational complexity of the algorithm
			chi2 += dy * dy;
		}
		
		// store result
		return chi2;
	}

	/**
	 *
	 */
	public double[][] calculateAlphaMatrix(double[] x_data, double[] y_data, int data_points, double[] parameters, double lambda)

	{
		double[][] alpha_matrix = new double[GAUSSIAN_1D_PARAMETERS][GAUSSIAN_1D_PARAMETERS];
		
		// diagonal
		//double alpha_res_00 = 0.0f; // sum of ones
		double alpha_res_00 = data_points; // OPTIMISATION: sum of ones
		double alpha_res_11 = 0.0f;
		double alpha_res_22 = 0.0f;
		double alpha_res_33 = 0.0f;
		
		// triangular matrix (triangular, since matrix is symmetrical)
		double alpha_res_01 = 0.0f;
		double alpha_res_02 = 0.0f;
		double alpha_res_03 = 0.0f;
		
		double alpha_res_12 = 0.0f;
		double alpha_res_13 = 0.0f;
		
		double alpha_res_23 = 0.0f;
		
		// loop over all data points, computating all intermediate values (only triangular matrix, since matrix is symmetrical)
		for(int p = 0; p < data_points; ++p)
		{
			
			// calculate derivates for position
			double xmxpos = x_data[p]-parameters[PARAM_MU];
			double xmxpos2 = xmxpos * xmxpos;
			double xsig2 = parameters[PARAM_SIGMA] * parameters[PARAM_SIGMA];
			double xsig3 = xsig2 * parameters[PARAM_SIGMA];
			double xmxpos2dxsig2 = xmxpos2 / xsig2;
			
			double common_exp = Math.exp(-0.5f * xmxpos2dxsig2);
			double amp_common_exp = parameters[PARAM_AMPLITUDE] * common_exp;
			
			//double dy_p0 = 1.0f; // optimisation
			double dy_p1 = common_exp;
			double dy_p2 = amp_common_exp * (xmxpos / xsig2);
			double dy_p3 = amp_common_exp * (xmxpos2 / xsig3);
			
			// sum product
			//alpha_res_00 += dy_p0 * dy_p0; // OPTIMISATION
			alpha_res_11 += dy_p1 * dy_p1;
			alpha_res_22 += dy_p2 * dy_p2;
			alpha_res_33 += dy_p3 * dy_p3;
			
			alpha_res_01 += dy_p1; // optimisation // dy_p0 (=1) * dy_p1;
			alpha_res_02 += dy_p2; // optimisation // dy_p0 (=1) * dy_p2;
			alpha_res_03 += dy_p3; // optimisation // dy_p0 (=1) * dy_p3;
			
			alpha_res_12 += dy_p1 * dy_p2;
			alpha_res_13 += dy_p1 * dy_p3;
			
			alpha_res_23 += dy_p2 * dy_p3;
		}
		
		// marquardt's lambda addition
		double lambda1 = lambda + 1.0f;
		alpha_res_00 *= lambda1;
		alpha_res_11 *= lambda1;
		alpha_res_22 *= lambda1;
		alpha_res_33 *= lambda1;
		
		// store results (expand computation from triangular matrix to full matrix)
		alpha_matrix[0][0] = alpha_res_00;
		alpha_matrix[0][1] = alpha_res_01;
		alpha_matrix[0][2] = alpha_res_02;
		alpha_matrix[0][3] = alpha_res_03;
		
		alpha_matrix[1][0] = alpha_res_01;
		alpha_matrix[1][1] = alpha_res_11;
		alpha_matrix[1][2] = alpha_res_12;
		alpha_matrix[1][3] = alpha_res_13;
		
		alpha_matrix[2][0] = alpha_res_02;
		alpha_matrix[2][1] = alpha_res_12;
		alpha_matrix[2][2] = alpha_res_22;
		alpha_matrix[2][3] = alpha_res_23;
		
		alpha_matrix[3][0] = alpha_res_03;
		alpha_matrix[3][1] = alpha_res_13;
		alpha_matrix[3][2] = alpha_res_23;
		alpha_matrix[3][3] = alpha_res_33;
		
		return alpha_matrix;
	}

	/**
	 *
	 */
	public double[] calculateBetaVector(double[] x_data, double[] y_data, int data_points, double[] parameters)
	{
		double[] beta_vector = new double[GAUSSIAN_1D_PARAMETERS];
		
		// loop over all pixels in image
		for(int p = 0; p < data_points; ++p)
		{
			double xmxpos = x_data[p]-parameters[PARAM_MU];
			double xmxpos2 = xmxpos * xmxpos;
			double xsig2 = parameters[PARAM_SIGMA] * parameters[PARAM_SIGMA];
			double xsig3 = xsig2 * parameters[PARAM_SIGMA];
			double xmxpos2dxsig2 = xmxpos2 / xsig2;
			
			double common_exp = Math.exp(-0.5f * xmxpos2dxsig2);
			double amp_common_exp = parameters[PARAM_AMPLITUDE] * common_exp;
			
			// compute intermediate shared result (difference)
			double px_dg = y_data[p] - gaussian1D(x_data[p], parameters);
			beta_vector[0] +=  px_dg;
			beta_vector[1] +=  px_dg * common_exp;
			beta_vector[2] +=  px_dg * amp_common_exp * (xmxpos / xsig2);
			beta_vector[3] +=  px_dg * amp_common_exp * (xmxpos2 / xsig3);
		}
		
		return beta_vector;
	}

	/**
	 *	Choleksy decomposition of matrix[R][C] done in place
	 */
	public double[][] choleskyDecomposition(double[][] matrix)
	{
		// loop over diagonals
		for(int j = 0; j < GAUSSIAN_1D_PARAMETERS; ++j)
		{
			double l_jj = matrix[j][j];
			for(int k = 0; k < j; ++k)
			{
				l_jj -= matrix[k][j] * matrix[k][j];
			}
			
			// RSLV: check if l_jj > 0 for positive-definite
			
			// take square root of l_jj
			l_jj = Math.sqrt(l_jj);
			
			// set new diagonal values (in place!)
			matrix[j][j] = l_jj;
			
			// calculate inverse of l_jj;
			double l_jj_inv = 1.0 / l_jj;
			
			// STEP: calculate values of submatrix
			for(int i = j+1; i < GAUSSIAN_1D_PARAMETERS; ++i)
			{
				double l_ji = matrix[j][i];
				for(int k = 0; k < j; ++k)
				{
					l_ji -= matrix[k][i] * matrix[k][j];
				}
				l_ji *= l_jj_inv;
				matrix[j][i] = l_ji;
			}
		}
		
		// reset upper triangle of matrix to zero
		for(int j = 1; j < GAUSSIAN_1D_PARAMETERS; ++j)
		{
			for(int i = 0; i < j; ++i)
			{
				matrix[j][i] = 0.0;
			}
		}
		
		// return matrix
		return matrix;
	}
	
	/**
	 *	Forward substitution of matrix
	 */
	public double[][] forwardSubstitute(double[][] l_matrix)
	{
		double[][] temp_matrix = new double[GAUSSIAN_1D_PARAMETERS][GAUSSIAN_1D_PARAMETERS];
		
		// initialize temp matrix with identity matrix
		/*for(int i = 0; i < GAUSSIAN_1D_PARAMETERS; ++i)
		{
			for(int j = i+1; j < GAUSSIAN_1D_PARAMETERS; ++j)
			{
				temp_matrix[i][j] = 0.0;
				temp_matrix[j][i] = 0.0;
			}
		}*/
		for(int d = 0; d < GAUSSIAN_1D_PARAMETERS; ++d)
		{
			temp_matrix[d][d] = 1.0f;
		}
		
		// keep track if inverse exists (i.e. matrix is non-singular)
		// boolean inverse_exists = true;
		
		// perform forward substitution
		for(int i = 0; i < GAUSSIAN_1D_PARAMETERS; ++i)
		{
			double pivot_value = l_matrix[i][i];
			//if(pivot_value == 0.0)
			//{
			//	inverse_exists = false;
			//	return null;
			//}
			
			// echelon elimination on succeeding (+1) rows
			for(int j = i+1; j < GAUSSIAN_1D_PARAMETERS; ++j)
			{
				double pivot_factor = l_matrix[i][j] / pivot_value;
				
				// subtract pivot factor multiple of pivot row from current row
				for(int k = 0; k < j; ++k)
				{
					l_matrix[k][j] -= pivot_factor * l_matrix[k][i];
					temp_matrix[k][j] -= pivot_factor * temp_matrix[k][i];
				}
			}
			
			// divide row by pivot value
			for(int k = 0; k <= i; ++k)
			{
				l_matrix[k][i] /= pivot_value;
				temp_matrix[k][i] /= pivot_value;
			}
		}
		
		return temp_matrix;
	}
	
	/**
	 *
	 */
	public double[][] multiplyMatrixSelf(double[][] matrix)
	{
		double[][] result_matrix = new double[GAUSSIAN_1D_PARAMETERS][GAUSSIAN_1D_PARAMETERS];
		
		// multiply the two inverse triangular matrices (in this case a special function that transposes the alpha matrix inline for the second matrix) to get the inverse alpha matrix
		for(int i = 0; i < GAUSSIAN_1D_PARAMETERS; ++i)
		{
			for(int j = 0; j < GAUSSIAN_1D_PARAMETERS; ++j)
			//for(__public int j = 0; j <= i; ++j) // RSLV: optimisation?
			{
				double mat_mult = 0.0f;
				for(int k = 0; k < GAUSSIAN_1D_PARAMETERS; ++k)
				{
					mat_mult += matrix[j][k] * matrix[i][k]; // NOTE in last term k,j is transposed to j,k!
				}
				result_matrix[i][j] = mat_mult;
			}
		}
		
		return result_matrix;
	}

	/**
	 *
	 */
	public double[] calculateDaVector(double[][] inverse_alpha_matrix, double[] beta_vector)
	{
		double[] da_vector = new double[GAUSSIAN_1D_PARAMETERS];
		
		// calculate da = inverse_alpha matrix * beta vector
		for(int i = 0; i < GAUSSIAN_1D_PARAMETERS; ++i)
		{
			double da_res = 0.0f;
			for(int j = 0; j < GAUSSIAN_1D_PARAMETERS; ++j)
			{
				da_res += beta_vector[j] * inverse_alpha_matrix[i][j];
			}
			
			// store result
			da_vector[i] = da_res;
		}
		
		return da_vector;
	}
	
	/**
	 *
	 */
	public double[] calculateUpdatedParameters(double[] fitted_parameters, double[] da_vector)
	{
		double[] updated_parameters = new double[GAUSSIAN_1D_PARAMETERS];
		
		for(int i = 0; i < GAUSSIAN_1D_PARAMETERS; ++i)
		{
			updated_parameters[i] = fitted_parameters[i] + da_vector[i];
		}
		
		return updated_parameters;
	}
	
	/**
	 *
	 */
	public double[] calculateStandardErrors(double[] x_data, double[] y_data, int data_points, double[] parameters)
	{
		double[] standard_errors = new double[GAUSSIAN_1D_PARAMETERS];
		
		// calculate alpha matrix
		double lambda = 1.0;
		double[][] alpha_matrix = calculateAlphaMatrix(x_data, y_data, data_points, parameters, lambda);
		
		// invert alpha matrix
		double[][] inverse_alpha_matrix = choleskyDecomposition(alpha_matrix);
		inverse_alpha_matrix = forwardSubstitute(inverse_alpha_matrix);
		inverse_alpha_matrix = multiplyMatrixSelf(inverse_alpha_matrix);
		
		for(int d = 0; d < GAUSSIAN_1D_PARAMETERS; ++d)
		{
			standard_errors[d] = Math.sqrt(inverse_alpha_matrix[d][d]);
		}
		
		return standard_errors;
	}
	
	/**
	 *
	 */
	public double[] run(double[] x_data, double[] y_data, int data_points, double[] initial_parameters, int num_iterations, double lambda)
	{
		double[] fitted_parameters = initial_parameters; //new double[GAUSSIAN_1D_PARAMETERS];
		
		// calculate chi2 value of current model
		double current_chi2 = calculateChi2(x_data, y_data, data_points, initial_parameters);
		
		// loop for fixed number of iterations
		int iteration_count = 0;
		while(iteration_count < num_iterations)
		{
			// calculate alpha matrix
			double[][] alpha_matrix = calculateAlphaMatrix(x_data, y_data, data_points, fitted_parameters, lambda);
			
			// calculate beta vector
			double[] beta_vector = calculateBetaVector(x_data, y_data, data_points, fitted_parameters);
			
			// invert alpha matrix
			double[][] inverse_alpha_matrix = choleskyDecomposition(alpha_matrix);
			inverse_alpha_matrix = forwardSubstitute(inverse_alpha_matrix);
			inverse_alpha_matrix = multiplyMatrixSelf(inverse_alpha_matrix);
			
			// calculate da vector
			double[] da_vector = calculateDaVector(inverse_alpha_matrix, beta_vector);
			
			// calculate updated parameters
			double[] updated_parameters = calculateUpdatedParameters(fitted_parameters, da_vector);
			
			// calculate updated chi2
			double updated_chi2 = calculateChi2(x_data, y_data, data_points, updated_parameters);
			
			// update parameters
			if(updated_chi2 <= current_chi2)
			{
				fitted_parameters = updated_parameters;
				current_chi2 = updated_chi2;
				lambda *= LAMBDA_FACTOR_INV; // decrease step size
			}
			else
			{
				lambda *= LAMBDA_FACTOR; // increase step size
			}
			
			// increase number of iterations
			++iteration_count;
		}
		
		return fitted_parameters;
	}
}
