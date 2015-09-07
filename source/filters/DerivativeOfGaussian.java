
package filters;

// import ImageJ classes
import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import ij.IJ.*;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;
//import ij.ImagePlus;


/**
 *	
 */
public class DerivativeOfGaussian
{
	public static final double DEFAULT_SIGMA = 1.0;
	public static final float DEFAULT_THETA = 0.0f;
	//public static final double DEFAULT_ANISOTROPY_RATIO = 1.0;
	
	/**
	 *	Constructor
	 */
	public DerivativeOfGaussian()
	{
		/* do nothing */
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static ImageProcessor derivativeX(ImageProcessor ip)
	{
		return derivativeX(ip, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor derivativeX(ImageProcessor ip, double sigma)
	{
		// convolve image with DoG_dx kernel
		//double[][] kernel = computeKernelGaussian2D_du(sigma);
		//return convolve(ip, kernel);
		return derivativeX(ip, sigma, DEFAULT_THETA);
	}
	
	public static ImageProcessor derivativeX(ImageProcessor ip, double sigma, float theta)
	{
		// convolve image with DoG_dx kernel
		//double[][] kernel = computeKernelGaussian2D_du(sigma);
		//return convolve(ip, kernel);
		return derivativeX(ip, sigma, sigma, theta);
	}
	
	public static ImageProcessor derivativeX(ImageProcessor ip, double sigma_x, double sigma_y)
	{
		return derivativeX(ip, sigma_x, sigma_y, DEFAULT_THETA);
	}
	
	public static ImageProcessor derivativeX(ImageProcessor ip, double sigma_x, double sigma_y, float theta)
	{
		// convolve image with DoG_dx kernel
		double[][] kernel = computeKernelGaussian2D_du(sigma_x, sigma_y, theta);
		return convolve(ip, kernel);
	}
	
	public static ImageProcessor derivativeY(ImageProcessor ip)
	{
		return derivativeY(ip, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor derivativeY(ImageProcessor ip, double sigma)
	{
		// convolve image with DoG_dy kernel
		//double[][] kernel = computeKernelGaussian2D_dv(sigma);
		//return convolve(ip, kernel);
		return derivativeY(ip, sigma, DEFAULT_THETA);
	}
	
	public static ImageProcessor derivativeY(ImageProcessor ip, double sigma, float theta)
	{
		// convolve image with DoG_dy kernel
		//double[][] kernel = computeKernelGaussian2D_dv(sigma);
		//return convolve(ip, kernel);
		return derivativeY(ip, sigma, sigma, theta);
	}
	
	public static ImageProcessor derivativeY(ImageProcessor ip, double sigma_x, double sigma_y)
	{
		return derivativeY(ip, sigma_x, sigma_y, DEFAULT_THETA);
	}
	
	public static ImageProcessor derivativeY(ImageProcessor ip, double sigma_x, double sigma_y, float theta)
	{
		// convolve image with DoG_dy kernel
		double[][] kernel = computeKernelGaussian2D_dv(sigma_x, sigma_y, theta);
		return convolve(ip, kernel);
	}
	
	public static ImageProcessor derivativeXX(ImageProcessor ip)
	{
		return derivativeXX(ip, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor derivativeXX(ImageProcessor ip, double sigma)
	{
		// convolve image with DoG_dx_dx kernel
		//double[][] kernel = computeKernelGaussian2D_du_du(sigma);
		//return convolve(ip, kernel);
		return derivativeXX(ip, sigma, DEFAULT_THETA);
	}
	
	public static ImageProcessor derivativeXX(ImageProcessor ip, double sigma, float theta)
	{
		// convolve image with DoG_dx_dx kernel
		//double[][] kernel = computeKernelGaussian2D_du_du(sigma);
		//return convolve(ip, kernel);
		return derivativeXX(ip, sigma, sigma, theta);
	}
	
	public static ImageProcessor derivativeXX(ImageProcessor ip, double sigma_x, double sigma_y)
	{
		return derivativeXX(ip, sigma_x, sigma_y, DEFAULT_THETA);
	}
	
	public static ImageProcessor derivativeXX(ImageProcessor ip, double sigma_x, double sigma_y, float theta)
	{
		// convolve image with DoG_dx_dx kernel
		double[][] kernel = computeKernelGaussian2D_du_du(sigma_x, sigma_y, theta);
		return convolve(ip, kernel);
	}
	
	public static ImageProcessor derivativeYY(ImageProcessor ip)
	{
		return derivativeYY(ip, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor derivativeYY(ImageProcessor ip, double sigma)
	{
		// convolve image with DoG_dy_dy kernel
		//double[][] kernel = computeKernelGaussian2D_dv_dv(sigma);
		//return convolve(ip, kernel);
		return derivativeYY(ip, sigma, DEFAULT_THETA);
	}
	
	public static ImageProcessor derivativeYY(ImageProcessor ip, double sigma, float theta)
	{
		// convolve image with DoG_dy_dy kernel
		//double[][] kernel = computeKernelGaussian2D_dv_dv(sigma);
		//return convolve(ip, kernel);
		return derivativeYY(ip, sigma, sigma, theta);
	}
	
	public static ImageProcessor derivativeYY(ImageProcessor ip, double sigma_x, double sigma_y)
	{
		return derivativeYY(ip, sigma_x, sigma_y, DEFAULT_THETA);
	}
	
	public static ImageProcessor derivativeYY(ImageProcessor ip, double sigma_x, double sigma_y, float theta)
	{
		// convolve image with DoG_dy_dy kernel
		double[][] kernel = computeKernelGaussian2D_dv_dv(sigma_x, sigma_y, theta);
		return convolve(ip, kernel);
	}
	
	public static ImageProcessor derivativeXY(ImageProcessor ip)
	{
		return derivativeXY(ip, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor derivativeXY(ImageProcessor ip, double sigma)
	{
		// convolve image with DoG_dx_dy kernel
		//double[][] kernel = computeKernelGaussian2D_du_dv(sigma);
		//return convolve(ip, kernel);
		return derivativeXY(ip, sigma, DEFAULT_THETA);
	}
	
	public static ImageProcessor derivativeXY(ImageProcessor ip, double sigma, float theta)
	{
		// convolve image with DoG_dx_dy kernel
		//double[][] kernel = computeKernelGaussian2D_du_dv(sigma);
		//return convolve(ip, kernel);
		return derivativeXY(ip, sigma, sigma, theta);
	}
	
	public static ImageProcessor derivativeXY(ImageProcessor ip, double sigma_x, double sigma_y)
	{
		return derivativeXY(ip, sigma_x, sigma_y, DEFAULT_THETA);
	}
	
	public static ImageProcessor derivativeXY(ImageProcessor ip, double sigma_x, double sigma_y, float theta)
	{
		// convolve image with DoG_dx_dy kernel
		double[][] kernel = computeKernelGaussian2D_du_dv(sigma_x, sigma_y, theta);
		return convolve(ip, kernel);
	}
	
	public static ImageProcessor derivativeYX(ImageProcessor ip)
	{
		return derivativeYX(ip, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor derivativeYX(ImageProcessor ip, double sigma)
	{
		// convolve image with DoG_dy_dx kernel
		//double[][] kernel = computeKernelGaussian2D_dv_du(sigma);
		//return convolve(ip, kernel);
		return derivativeYX(ip, sigma, DEFAULT_THETA);
	}
	
	public static ImageProcessor derivativeYX(ImageProcessor ip, double sigma, float theta)
	{
		// convolve image with DoG_dy_dx kernel
		//double[][] kernel = computeKernelGaussian2D_dv_du(sigma);
		//return convolve(ip, kernel);
		return derivativeYX(ip, sigma, sigma, theta);
	}
	
	public static ImageProcessor derivativeYX(ImageProcessor ip, double sigma_x, double sigma_y)
	{
		return derivativeYX(ip, sigma_x, sigma_y, DEFAULT_THETA);
	}
	
	public static ImageProcessor derivativeYX(ImageProcessor ip, double sigma_x, double sigma_y, float theta)
	{
		// convolve image with DoG_dy_dx kernel
		double[][] kernel = computeKernelGaussian2D_dv_du(sigma_x, sigma_y, theta);
		return convolve(ip, kernel);
	}
	
	public static double[][][] get_line_points(ImageProcessor ip, double sigma)
	{
		int image_width = ip.getWidth();
		int image_height = ip.getHeight();
		
		//final line points
		double[][][] line_points = new double[image_width][image_height][8];
		// store results of eigendecomposition
		//	[px][py][0] = lambda1_magnitude		n(t)
		//	[px][py][1] = lambda1_direction_x	n_x(t)
		//	[px][py][2] = lambda1_direction_y	n_y(t)
		//	[px][py][3] = lambda2_magnitude		s(t)
		//	[px][py][4] = lambda2_direction_x	s_x(t)
		//	[px][py][5] = lambda2_direction_y	s_y(t)
		//	[px][py][6] = super-resolved_x		t_x, or dlpx
		//	[px][py][7] = super-resolved_y		t_y, or dlpy		

		// calculate derivatives of gaussian from image
		ImageProcessor dx = derivativeX(ip, sigma);
		ImageProcessor dy = derivativeY(ip, sigma);
		ImageProcessor dxdx = derivativeXX(ip, sigma);
		ImageProcessor dxdy = derivativeXY(ip, sigma);
		ImageProcessor dydx = dxdy;//DerivativeOfGaussian.derivativeYX(ip, SIGMA);
		ImageProcessor dydy = derivativeYY(ip, sigma);

		//calculate Eigen decomposition of Hessian matrix
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				Matrix m = new Matrix(2, 2, 0); // 2x2 RC matrix with zeros
				m.set(0, 0, dxdx.getf(px, py));
				m.set(0, 1, dxdy.getf(px, py));
				m.set(1, 0, dydx.getf(px, py));
				m.set(1, 1, dydy.getf(px, py));
				
				// compute eigenvalues and eigenvectors
				EigenvalueDecomposition evd = m.eig();
				Matrix d = evd.getD();
				Matrix v = evd.getV();
				
				// determine first and second eigenvalue and eigenvector
				double first_eigenvalue = 0.0; // |n(t)|
				double first_eigenvector_x = 0.0; // n(t) -> perpendicular to s(t)
				double first_eigenvector_y = 0.0; // n(t) -> perpendicular to s(t)
				double second_eigenvalue = 0.0;
				double second_eigenvector_x = 0.0;
				double second_eigenvector_y = 0.0;
				
				if(d.get(0,0) <= d.get(1,1))
				{
					// d(0,0) is most negative minimum eigenvalue
					first_eigenvalue = d.get(0,0); // L1
					first_eigenvector_x = v.get(0,0); // V1x
					first_eigenvector_y = v.get(1,0); // V1y
					second_eigenvalue = d.get(1,1); // L2
					second_eigenvector_x = v.get(0,1); // V2x
					second_eigenvector_y = v.get(1,1); // V2y
				}
				else
				{
					// d(1,1) is most negative minimum eigenvalue
					first_eigenvalue = d.get(1,1); // L1
					first_eigenvector_x = v.get(0,1); // V1x
					first_eigenvector_y = v.get(1,1); // V1y
					second_eigenvalue = d.get(0,0); // L2
					second_eigenvector_x = v.get(0,0); // V2x
					second_eigenvector_y = v.get(1,0); // V2y
				}
				
				// reorient all eigenvecors in same direction
				if(first_eigenvector_y < 0)
				{
					first_eigenvector_x = -first_eigenvector_x;
					first_eigenvector_y = -first_eigenvector_y;
				}
				
				if(second_eigenvector_y < 0)
				{
					second_eigenvector_x = -second_eigenvector_x;
					second_eigenvector_y = -second_eigenvector_y;
				}
				
				// store eigenvalues and eigenvector for new optimum
				line_points[px][py][0] = first_eigenvalue;
				line_points[px][py][1] = first_eigenvector_x;
				line_points[px][py][2] = first_eigenvector_y;
				line_points[px][py][3] = second_eigenvalue;
				line_points[px][py][4] = second_eigenvector_x;
				line_points[px][py][5] = second_eigenvector_y;
				
				// calculate position of peak in second order Taylor polynomial from Steger's algorithm
				double t = (dx.getf(px,py)*first_eigenvector_x + dy.getf(px,py)*first_eigenvector_y)/(dxdx.getf(px,py)*first_eigenvector_x*first_eigenvector_x + 2.0*dxdy.getf(px,py)*first_eigenvector_x*first_eigenvector_y + dydy.getf(px,py)*first_eigenvector_y*first_eigenvector_y);
				double dlpx = t*first_eigenvector_x;
				double dlpy = t*first_eigenvector_y;
				
				// store line point
				line_points[px][py][6] = dlpx;
				line_points[px][py][7] = dlpy;
			}
		}
		
		
		return line_points;
	}
	
	/** the same as line point above, only in "Steger's" format 
	    used for verification and debug */
	public static double[][] get_line_points_stegers(ImageProcessor ip, double sigma, boolean mode)
	{
		int image_width = ip.getWidth();
		int image_height = ip.getHeight();
		int l;
		double a;
		double b;
	
		
		//final line points
		double[][] line_points = new double[7][image_width*image_height];
		// store results of eigendecomposition
		//	[0] = lambda1_magnitude		n(t)
		//	[1] = lambda1_direction_x	n_x(t)
		//	[2] = lambda1_direction_y	n_y(t)
		//	[3] = derivative of image in y
		//	[4] = derivative of image in x
		//	[5] = super-resolved_x		t_y, or dlpy
		//	[6] = super-resolved_y		t_x, or dlpx		


		// calculate derivatives of gaussian from image
		ImageProcessor dx = derivativeX(ip, sigma);
		ImageProcessor dy = derivativeY(ip, sigma);
		ImageProcessor dxdx = derivativeXX(ip, sigma);
		ImageProcessor dxdy = derivativeXY(ip, sigma);
		ImageProcessor dydy = derivativeYY(ip, sigma);

		//calculate Eigen decomposition of Hessian matrix
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{

				l = LINCOOR(py,px,image_width);
				
				// determine first and second eigenvalue and eigenvector
				double first_eigenvalue = 0.0; // |n(t)|
				double first_eigenvector_x = 0.0; // n(t) -> perpendicular to s(t)
				double first_eigenvector_y = 0.0; // n(t) -> perpendicular to s(t)							
				double t;
				
				
				//*
				double dfdrc = dxdy.getf(px, py);
				double dfdrr = dxdx.getf(px, py);
				double dfdcc = dydy.getf(px, py);

				double theta, e1,e2,n1,n2,c,s;
				
				 if (dfdrc != 0.0) 
				 {
					    theta = 0.5*(dfdcc-dfdrr)/dfdrc;
					    t = 1.0/(Math.abs(theta)+Math.sqrt(theta*theta+1.0));
					    if (theta < 0.0) t = -t;
					    c = 1.0/Math.sqrt(t*t+1.0);
					    s = t*c;
					    e1 = dfdrr-t*dfdrc;
					    e2 = dfdcc+t*dfdrc;
					  
				 } else 
				 {
					    c = 1.0;
					    s = 0.0;
					    e1 = dfdrr;
					    e2 = dfdcc;
					  
				 }
				 n1 = c;
				 n2 = -s;
				 //actually since hessian is a symmetric matrix
				 //eigen vectors are perpendicular, so no need to store another vector
				
				 
				  if (Math.abs(e1) > Math.abs(e2)) 
				  {
					  first_eigenvalue = e1;
					    //eigval[1] = e2;
					    first_eigenvector_x = n1;
					    first_eigenvector_y = n2;
					    //eigvec[1][0] = -n2;
					    //eigvec[1][1] = n1;
				  } 
				  else if (Math.abs(e1) < Math.abs(e2))
				  {
						  first_eigenvalue = e2;
					    //eigval[1] = e1;
						  first_eigenvector_x = -n2;
						  first_eigenvector_y = n1;
					    //eigvec[1][0] = n1;
					    //eigvec[1][1] = n2;
				   } 
				   else 
				   {
					    if (e1 < e2) 
					    {
					    	first_eigenvalue = e1;
					      //eigval[1] = e2;
					    	first_eigenvector_x = n1;
					    	first_eigenvector_y = n2;
					      //eigvec[1][0] = -n2;
					      //eigvec[1][1] = n1;
					    } else 
					    {
					    	first_eigenvalue = e2;
					      //eigval[1] = e1;
					    	first_eigenvector_x = -n2;
					    	first_eigenvector_y = n1;
					      //eigvec[1][0] = n1;
					      //eigvec[1][1] = n2;
					    }
					}
					//*/ 
				
				
				
				// store eigenvalues and eigenvector for new optimum
				  if (mode)
					  line_points[0][l] = first_eigenvalue*(-1);
				  else
					  line_points[0][l] = first_eigenvalue;
				line_points[1][l] = first_eigenvector_y;
				line_points[2][l] = first_eigenvector_x;				
				line_points[3][l] = dy.getf(px,py);
				line_points[4][l] = dx.getf(px,py);
				
				// calculate position of peak in second order Taylor polynomial from Steger's algorithm
				b = (dx.getf(px,py)*first_eigenvector_x + dy.getf(px,py)*first_eigenvector_y);
				a = (dxdx.getf(px,py)*first_eigenvector_x*first_eigenvector_x + 2.0*dxdy.getf(px,py)*first_eigenvector_x*first_eigenvector_y + dydy.getf(px,py)*first_eigenvector_y*first_eigenvector_y);
				//a = (dxdx.getf(px,py)*first_eigenvector_x*first_eigenvector_x + dxdy.getf(px,py)*dxdy.getf(px,py)*first_eigenvector_x*first_eigenvector_y + dydy.getf(px,py)*first_eigenvector_y*first_eigenvector_y);
				t = (+1)*b/a;
				double dlpx = t*first_eigenvector_x;
				double dlpy = t*first_eigenvector_y;
				
				// store line point
				line_points[5][l] = dlpy+py;
				line_points[6][l] = dlpx+px;								
								
			}
		}
		
		
		return line_points;
	}
	
	
	
	
	public static double[][][] frangi_measures(double[][][] line_points, int image_height, int image_width, double frangi_beta, double frangi_c)
	{
		// store Frangi measures on eigenvalues; NOTE: |L1| <= |L2|
		// beta is control parameter, set at 0.5
		// c dependens on image bit depth, about half maximum Hessian matrix norm
		//	[0] = frangi L1
		//	[1] = frangi L2
		//	[2] = blobness (eccentricity), L1 / L2; note: keep sign!
		//	[3] = second order structureness, RSS of eigenvalues, or Frobius norm
		//	[4] = vesselness = exp(-[2]^2/2*FRANGI_BETA^2)(1-exp(-[3]^2/2*FRANGI_C^2))

		double[][][] frangi_measures = new double[image_width][image_height][5];
		double frangi_l1;
		double frangi_l2;
		double first_eigenvalue, second_eigenvalue; 
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				// calculate Frangi measures
				first_eigenvalue = line_points[px][py][0];
				second_eigenvalue = line_points[px][py][3];

				if(Math.abs(first_eigenvalue) > Math.abs(second_eigenvalue))
				{
					frangi_l1 = second_eigenvalue;
					frangi_l2 = first_eigenvalue;
				}
				else
				{
					frangi_l1 = first_eigenvalue;
					frangi_l2 = second_eigenvalue;
				}
				
				double eccentricity = frangi_l1 / (frangi_l2 + 1e20); // NOTE: beware of division by zero!
				double structureness = Math.sqrt(frangi_l1*frangi_l1+frangi_l2*frangi_l2); // ro Frobius norm
				double fm = Math.exp(-((eccentricity*eccentricity)/(2*frangi_beta*frangi_beta))) * (1-Math.exp(-((structureness*structureness))/(2*frangi_c*frangi_c)));
				if(frangi_l2 > 0)
				{
					fm = 0.0;
				}
				// store Frangi measures
				frangi_measures[px][py][0] = frangi_l1;
				frangi_measures[px][py][1] = frangi_l2;
				frangi_measures[px][py][2] = eccentricity;
				frangi_measures[px][py][3] = structureness;
				frangi_measures[px][py][4] = fm;
				
			}			
		}
		return frangi_measures;
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static double[][] computeKernelGaussian2D()
	{
		return computeKernelGaussian2D(DEFAULT_SIGMA);
	}
	
	public static double[][] computeKernelGaussian2D(double sigma)
	{
		return computeKernelGaussian2D(sigma, DEFAULT_THETA);
	}
	
	public static double[][] computeKernelGaussian2D(double sigma, float theta)
	{
		return computeKernelGaussian2D(sigma, sigma, theta);
	}
	
	public static double[][] computeKernelGaussian2D(double sigma_x, double sigma_y)
	{
		return computeKernelGaussian2D(sigma_x, sigma_y, DEFAULT_THETA);
	}
	
	public static double[][] computeKernelGaussian2D(double sigma_x, double sigma_y, float theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*Math.max(sigma_x, sigma_y)); // RSLV: use floor instead of round?
		int kernel_size = 1+2*kernel_radius;
		
		// compute kernel
		double[][] kernel = new double[kernel_size][kernel_size];
		for(int ky = 0; ky < kernel_size; ++ky)
		{
			int y = ky - kernel_radius;
			for(int kx = 0; kx < kernel_size; ++kx)
			{
				int x = kx - kernel_radius;
				double u = x * Math.cos(theta) - y * Math.sin(theta);
				double v = x * Math.sin(theta) + y * Math.cos(theta);
				kernel[kx][ky] = gaussian2D(u, v, sigma_x, sigma_y);
			}
		}
		
		// normalize kernel
		kernel = normalize_kernel(kernel);
		
		// return kernel
		return kernel;
	}
	
	public static double[][] computeKernelGaussian2D_du()
	{
		return computeKernelGaussian2D_du(DEFAULT_SIGMA);
	}
	
	public static double[][] computeKernelGaussian2D_du(double sigma)
	{
		return computeKernelGaussian2D_du(sigma, DEFAULT_THETA);
	}
	
	public static double[][] computeKernelGaussian2D_du(double sigma, float theta)
	{
		return computeKernelGaussian2D_du(sigma, sigma, theta);
	}
	
	public static double[][] computeKernelGaussian2D_du(double sigma_x, double sigma_y)
	{
		return computeKernelGaussian2D_du(sigma_x, sigma_y, DEFAULT_THETA);
	}
	
	public static double[][] computeKernelGaussian2D_du(double sigma_x, double sigma_y, float theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*Math.max(sigma_x, sigma_y)); // RSLV: use floor instead of round?
		int kernel_size = 1+2*kernel_radius;
		
		// compute kernel
		double[][] kernel = new double[kernel_size][kernel_size];
		for(int ky = 0; ky < kernel_size; ++ky)
		{
			int y = ky - kernel_radius;
			for(int kx = 0; kx < kernel_size; ++kx)
			{
				int x = kx - kernel_radius;
				double u = x * Math.cos(theta) - y * Math.sin(theta);
				double v = x * Math.sin(theta) + y * Math.cos(theta);
				kernel[kx][ky] = gaussian2D_dx(u, v, sigma_x, sigma_y);
			}
		}
		
		// normalize kernel
		kernel = normalize_kernel(kernel);
		
		// return kernel
		return kernel;
	}
	
	public static double[][] computeKernelGaussian2D_dv()
	{
		return computeKernelGaussian2D_dv(DEFAULT_SIGMA);
	}
	
	public static double[][] computeKernelGaussian2D_dv(double sigma)
	{
		return computeKernelGaussian2D_dv(sigma, DEFAULT_THETA);
	}
	
	public static double[][] computeKernelGaussian2D_dv(double sigma, float theta)
	{
		return computeKernelGaussian2D_dv(sigma, sigma, theta);
	}
	
	public static double[][] computeKernelGaussian2D_dv(double sigma_x, double sigma_y)
	{
		return computeKernelGaussian2D_dv(sigma_x, sigma_y, DEFAULT_THETA);
	}
	
	public static double[][] computeKernelGaussian2D_dv(double sigma_x, double sigma_y, float theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*Math.max(sigma_x, sigma_y)); // RSLV: use floor instead of round?
		int kernel_size = 1+2*kernel_radius;
		
		// compute kernel
		double[][] kernel = new double[kernel_size][kernel_size];
		for(int ky = 0; ky < kernel_size; ++ky)
		{
			int y = ky - kernel_radius;
			for(int kx = 0; kx < kernel_size; ++kx)
			{
				int x = kx - kernel_radius;
				double u = x * Math.cos(theta) - y * Math.sin(theta);
				double v = x * Math.sin(theta) + y * Math.cos(theta);
				kernel[kx][ky] = gaussian2D_dy(u, v, sigma_x, sigma_y);
			}
		}
		
		// normalize kernel
		kernel = normalize_kernel(kernel);
		
		// return kernel
		return kernel;
	}
	
	public static double[][] computeKernelGaussian2D_du_du()
	{
		return computeKernelGaussian2D_du_du(DEFAULT_SIGMA);
	}
	
	public static double[][] computeKernelGaussian2D_du_du(double sigma)
	{
		return computeKernelGaussian2D_du_du(sigma, DEFAULT_THETA);
	}
	
	public static double[][] computeKernelGaussian2D_du_du(double sigma, float theta)
	{
		return computeKernelGaussian2D_du_du(sigma, sigma, theta);
	}
	
	public static double[][] computeKernelGaussian2D_du_du(double sigma_x, double sigma_y)
	{
		return computeKernelGaussian2D_du_du(sigma_x, sigma_y, DEFAULT_THETA);
	}
	
	public static double[][] computeKernelGaussian2D_du_du(double sigma_x, double sigma_y, float theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*Math.max(sigma_x, sigma_y)); // RSLV: use floor instead of round?
		int kernel_size = 1+2*kernel_radius;
		
		// compute kernel
		double[][] kernel = new double[kernel_size][kernel_size];
		for(int ky = 0; ky < kernel_size; ++ky)
		{
			int y = ky - kernel_radius;
			for(int kx = 0; kx < kernel_size; ++kx)
			{
				int x = kx - kernel_radius;
				double u = x * Math.cos(theta) - y * Math.sin(theta);
				double v = x * Math.sin(theta) + y * Math.cos(theta);
				kernel[kx][ky] = gaussian2D_dx_dx(u, v, sigma_x, sigma_y);
			}
		}
		
		// normalize kernel
		kernel = normalize_kernel(kernel);
		
		// return kernel
		return kernel;
	}
	
	public static double[][] computeKernelGaussian2D_dv_dv()
	{
		return computeKernelGaussian2D_dv_dv(DEFAULT_SIGMA);
	}
	
	public static double[][] computeKernelGaussian2D_dv_dv(double sigma)
	{
		return computeKernelGaussian2D_dv_dv(sigma, DEFAULT_THETA);
	}
	
	public static double[][] computeKernelGaussian2D_dv_dv(double sigma, float theta)
	{
		return computeKernelGaussian2D_dv_dv(sigma, sigma, theta);
	}
	
	public static double[][] computeKernelGaussian2D_dv_dv(double sigma_x, double sigma_y)
	{
		return computeKernelGaussian2D_dv_dv(sigma_x, sigma_y, DEFAULT_THETA);
	}
	
	public static double[][] computeKernelGaussian2D_dv_dv(double sigma_x, double sigma_y, float theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*Math.max(sigma_x, sigma_y)); // RSLV: use floor instead of round?
		int kernel_size = 1+2*kernel_radius;
		
		// compute kernel
		double[][] kernel = new double[kernel_size][kernel_size];
		for(int ky = 0; ky < kernel_size; ++ky)
		{
			int y = ky - kernel_radius;
			for(int kx = 0; kx < kernel_size; ++kx)
			{
				int x = kx - kernel_radius;
				double u = x * Math.cos(theta) - y * Math.sin(theta);
				double v = x * Math.sin(theta) + y * Math.cos(theta);
				kernel[kx][ky] = gaussian2D_dy_dy(u, v, sigma_x, sigma_y);
			}
		}
		
		// normalize kernel
		kernel = normalize_kernel(kernel);
		
		// return kernel
		return kernel;
	}
	
	public static double[][] computeKernelGaussian2D_du_dv()
	{
		return computeKernelGaussian2D_du_dv(DEFAULT_SIGMA);
	}
	
	public static double[][] computeKernelGaussian2D_du_dv(double sigma)
	{
		return computeKernelGaussian2D_du_dv(sigma, DEFAULT_THETA);
	}
	
	public static double[][] computeKernelGaussian2D_du_dv(double sigma, float theta)
	{
		return computeKernelGaussian2D_du_dv(sigma, sigma, theta);
	}
	
	public static double[][] computeKernelGaussian2D_du_dv(double sigma_x, double sigma_y)
	{
		return computeKernelGaussian2D_du_dv(sigma_x, sigma_y, DEFAULT_THETA);
	}
	
	public static double[][] computeKernelGaussian2D_du_dv(double sigma_x, double sigma_y, float theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*Math.max(sigma_x, sigma_y)); // RSLV: use floor instead of round?
		int kernel_size = 1+2*kernel_radius;
		
		// compute kernel
		double[][] kernel = new double[kernel_size][kernel_size];
		for(int ky = 0; ky < kernel_size; ++ky)
		{
			int y = ky - kernel_radius;
			for(int kx = 0; kx < kernel_size; ++kx)
			{
				int x = kx - kernel_radius;
				double u = x * Math.cos(theta) - y * Math.sin(theta);
				double v = x * Math.sin(theta) + y * Math.cos(theta);
				kernel[kx][ky] = gaussian2D_dx_dy(u, v, sigma_x, sigma_y);
			}
		}
		
		// normalize kernel
		kernel = normalize_kernel(kernel);
		
		// return kernel
		return kernel;
	}
	
	public static double[][] computeKernelGaussian2D_dv_du()
	{
		return computeKernelGaussian2D_dv_du(DEFAULT_SIGMA);
	}
	
	public static double[][] computeKernelGaussian2D_dv_du(double sigma)
	{
		return computeKernelGaussian2D_dv_du(sigma, DEFAULT_THETA);
	}
	
	public static double[][] computeKernelGaussian2D_dv_du(double sigma, float theta)
	{
		return computeKernelGaussian2D_dv_du(sigma, sigma, theta);
	}
	
	public static double[][] computeKernelGaussian2D_dv_du(double sigma_x, double sigma_y)
	{
		return computeKernelGaussian2D_dv_du(sigma_x, sigma_y, DEFAULT_THETA);
	}
	
	public static double[][] computeKernelGaussian2D_dv_du(double sigma_x, double sigma_y, float theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*Math.max(sigma_x, sigma_y)); // RSLV: use floor instead of round?
		int kernel_size = 1+2*kernel_radius;
		
		// compute kernel
		double[][] kernel = new double[kernel_size][kernel_size];
		for(int ky = 0; ky < kernel_size; ++ky)
		{
			int y = ky - kernel_radius;
			for(int kx = 0; kx < kernel_size; ++kx)
			{
				int x = kx - kernel_radius;
				double u = x * Math.cos(theta) - y * Math.sin(theta);
				double v = x * Math.sin(theta) + y * Math.cos(theta);
				kernel[kx][ky] = gaussian2D_dy_dx(u, v, sigma_x, sigma_y);
			}
		}
		
		// normalize kernel
		kernel = normalize_kernel(kernel);
		
		// return kernel
		return kernel;
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static double gaussian2D(double x, double y)
	{
		return gaussian2D(x, y, DEFAULT_SIGMA);
	}
	
	public static double gaussian2D(double x, double y, double sigma)
	{
		//return (1/(2*Math.PI*Math.pow(sigma, 2)))*Math.exp(-(x*x+y*y)/(2*sigma*sigma));
		return gaussian2D(x, y, sigma, sigma);
	}
	
	public static double gaussian2D(double x, double y, double sigma_x, double sigma_y)
	{
		return (1/(2*Math.PI*sigma_x*sigma_y))*Math.exp(-0.5*((x*x)/(sigma_x*sigma_x)+(y*y)/(sigma_y*sigma_y)));
	}
	
	public static double gaussian2D_dx(double x, double y)
	{
		return gaussian2D_dx(x, y, DEFAULT_SIGMA);
	}
	
	public static double gaussian2D_dx(double x, double y, double sigma)
	{
		//return ((-x)/(2*Math.PI*Math.pow(sigma, 4)))*Math.exp(-(x*x+y*y)/(2*sigma*sigma));
		return gaussian2D_dx(x, y, sigma, sigma);
	}
	
	public static double gaussian2D_dx(double x, double y, double sigma_x, double sigma_y)
	{
		return ((-x)/(2*Math.PI*Math.pow(sigma_x, 3)*sigma_y))*Math.exp(-0.5*((x*x)/(sigma_x*sigma_x)+(y*y)/(sigma_y*sigma_y)));
	}
	
	public static double gaussian2D_dy(double x, double y)
	{
		return gaussian2D_dy(x, y, DEFAULT_SIGMA);
	}
	
	public static double gaussian2D_dy(double x, double y, double sigma)
	{
		//return gaussian2D_dx(y, x, sigma); // NOTE x and y are swapped
		//return ((-y)/(2*Math.PI*Math.pow(sigma, 4)))*Math.exp(-(x*x+y*y)/(2*sigma*sigma));
		return gaussian2D_dy(y, x, sigma, sigma); // NOTE x and y are swapped
	}
	
	public static double gaussian2D_dy(double x, double y, double sigma_x, double sigma_y)
	{
		//return gaussian2D_dx(y, x, sigma); // NOTE x and y are swapped
		return ((-y)/(2*Math.PI*sigma_x*Math.pow(sigma_y, 3)))*Math.exp(-0.5*((x*x)/(sigma_x*sigma_x)+(y*y)/(sigma_y*sigma_y)));
	}
	
	public static double gaussian2D_dx_dx(double x, double y)
	{
		return gaussian2D_dx_dx(x, y, DEFAULT_SIGMA);
	}
	
	public static double gaussian2D_dx_dx(double x, double y, double sigma)
	{
		//return ((x*x-sigma*sigma)/(2*Math.PI*Math.pow(sigma, 6)))*Math.exp(-(x*x+y*y)/(2*sigma*sigma));
		return gaussian2D_dx_dx(x, y, sigma, sigma);
	}
	
	public static double gaussian2D_dx_dx(double x, double y, double sigma_x, double sigma_y)
	{
		return ((x*x-sigma_x*sigma_x)/(2*Math.PI*Math.pow(sigma_x, 5)*sigma_y))*Math.exp(-0.5*((x*x)/(sigma_x*sigma_x)+(y*y)/(sigma_y*sigma_y)));
	}
	
	public static double gaussian2D_dy_dy(double x, double y)
	{
		return gaussian2D_dy_dy(x, y, DEFAULT_SIGMA);
	}
	
	public static double gaussian2D_dy_dy(double x, double y, double sigma)
	{
		//return gaussian2D_dx_dx(y, x, sigma); // NOTE x and y are swapped
		return gaussian2D_dy_dy(x, y, sigma, sigma);
	}
	
	public static double gaussian2D_dy_dy(double x, double y, double sigma_x, double sigma_y)
	{
		//return gaussian2D_dx_dx(y, x, sigma_x, sigma_y); // NOTE x and y are swapped
		return ((y*y-sigma_y*sigma_y)/(2*Math.PI*Math.pow(sigma_y, 5)*sigma_x))*Math.exp(-0.5*((x*x)/(sigma_x*sigma_x)+(y*y)/(sigma_y*sigma_y)));
	}
	
	public static double gaussian2D_dx_dy(double x, double y)
	{
		return gaussian2D_dx_dy(x, y, DEFAULT_SIGMA);
	}
	
	public static double gaussian2D_dx_dy(double x, double y, double sigma)
	{
		// return ((x*y)/(2*Math.PI*Math.pow(sigma, 6)))*Math.exp(-(x*x+y*y)/(2*sigma*sigma));
		return gaussian2D_dx_dy(x, y, sigma, sigma);
	}
	
	public static double gaussian2D_dx_dy(double x, double y, double sigma_x, double sigma_y)
	{
		return ((x*y)/(2*Math.PI*Math.pow(sigma_x, 3)*Math.pow(sigma_y, 3)))*Math.exp(-0.5*((x*x)/(sigma_x*sigma_x)+(y*y)/(sigma_y*sigma_y)));
	}
	
	public static double gaussian2D_dy_dx(double x, double y)
	{
		return gaussian2D_dy_dx(x, y, DEFAULT_SIGMA);
	}
	
	public static double gaussian2D_dy_dx(double x, double y, double sigma)
	{
		//return gaussian2D_dx_dy(y, x, sigma); // NOTE x and y are swapped (although not necessary since dx_dy == dy_dx)
		return gaussian2D_dy_dx(x, y, sigma, sigma);
	}
	
	public static double gaussian2D_dy_dx(double x, double y, double sigma_x, double sigma_y)
	{
		//return gaussian2D_dx_dy(y, x, sigma_x, sigma_y); // NOTE x and y are swapped (although not necessary since dx_dy == dy_dx)
		return ((x*y)/(2*Math.PI*Math.pow(sigma_x, 3)*Math.pow(sigma_y, 3)))*Math.exp(-0.5*((x*x)/(sigma_x*sigma_x)+(y*y)/(sigma_y*sigma_y)));
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static double[][] normalize_kernel(double[][] kernel)
	{
		// calculate sum of components
		double sum = 0.0;
		for(int kx = 0; kx < kernel.length; ++kx)
		{
			for(int ky = 0; ky < kernel[kx].length; ++ky)
			{
				sum += Math.abs(kernel[kx][ky]); // NOTE: use abs to normalize symmetrical kernel with a positive and negative lobe
			}
		}
		
		// avoid division by zero
		if(sum == 0.0) { return kernel; }
		
		// calculate scale factor
		double scale_factor = 1 / sum;
		
		// scale components
		for(int kx = 0; kx < kernel.length; ++kx)
		{
			for(int ky = 0; ky < kernel[kx].length; ++ky)
			{
				kernel[kx][ky] *= scale_factor;
			}
		}
		
		// return normalized kernel
//		ImageProcessor kernel_ip = kernelToImage(kernel);
//		ImagePlus kernel_imp = new ImagePlus("normalized kernel", kernel_ip);
//		kernel_imp.resetDisplayRange();
//		kernel_imp.show();
		return kernel;
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	/*
	public static ImageProcessor convolve(ImageProcessor ip, float[] kernel, int kernel_width, int kernel_height)
	{
		ij.plugin.filter.Convolver conv = new ij.plugin.filter.Convolver();
		conv.setNormalize(false);
		
		ImageProcessor ip_dup = ip.duplicate();
		
		conv.convolve(ip_dup, kernel, kernel_width, kernel_height);
		
		return ip_dup;
		
	}
	*/
	
	public static ImageProcessor convolve(ImageProcessor ip, double[][] kernel)
	{
		// get image and kernel sizes
		int image_width = ip.getWidth();
		int image_height = ip.getHeight();
		
		int kernel_width = kernel.length;
		int kernel_height = kernel_width; // NOTE: assume square kernel
		int kernel_half_width = (int)Math.floor(0.5 * kernel_width);
		int kernel_half_height = kernel_half_width;
		
		// convert input image processor to float
		ImageProcessor ip_inp = ip.convertToFloat(); // RSLV: duplicate first?
		
		// create new empty output float processor
		ImageProcessor ip_res = new FloatProcessor(image_width, image_height);
		
		// convolve input image with kernel
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				double kernel_product = 0.0;
				for(int ky = 0; ky < kernel_height; ++ky)
				{
					int ppy = py + ky - kernel_half_height;
					if(ppy < 0) ppy = 0; // clamp at border
					if(ppy >= image_height) ppy = image_height - 1; // clamp at border
					for(int kx = 0; kx < kernel_width; ++kx)
					{
						int ppx = px + kx - kernel_half_width;
						if(ppx < 0) ppx = 0; // clamp at border
						if(ppx >= image_width) ppx = image_width - 1; // clamp at border
						kernel_product += ip_inp.getf(ppx, ppy) * kernel[kx][ky];
					}
				}
				ip_res.setf(px, py, (float)kernel_product);
			}
		}
		
		// return 
		ip_res.resetMinAndMax();
		return ip_res;
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static ImageProcessor directionalDerivative(ImageProcessor ip, float theta)
	{
		return directionalDerivative(ip, theta, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor directionalDerivative(ImageProcessor ip, float theta, double sigma)
	{
		return directionalDerivative(ip, theta, sigma, sigma);
	}
	
	// NOTE: not sure if this trick also works for anistropic kernels
	public static ImageProcessor directionalDerivative(ImageProcessor ip, float theta, double sigma_x, double sigma_y)
	{
		ImageProcessor ip_dx = derivativeX(ip, sigma_x, sigma_y);
		ImageProcessor ip_dy = derivativeY(ip, sigma_x, sigma_y);
		return directionalDerivative(ip_dx, ip_dy, theta);
	}
	
	public static ImageProcessor directionalDerivative(ImageProcessor ip_dx, ImageProcessor ip_dy, float theta)
	{
		// assert: images ip_dx and ip_dy are of same size
		ImageProcessor ip_output = new FloatProcessor(ip_dx.getWidth(), ip_dy.getHeight());
		for(int y = 0; y < ip_dx.getHeight(); ++y)
		{
			for(int x = 0; x < ip_dx.getWidth(); ++x)
			{
				ip_output.setf(x, y, (float)(Math.cos(theta) * ip_dx.getf(x, y) + Math.sin(theta) * ip_dy.getf(x, y)));
			}
		}
		ip_output.resetMinAndMax();
		return ip_output;
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static ImageProcessor kernelToImage(double[][] kernel)
	{
		FloatProcessor ip = new FloatProcessor(kernel.length, kernel.length); // NOTE: assume sqaure kernel
		for(int i = 0; i < kernel.length; ++i)
		{
			for(int j = 0; j < kernel.length; ++j)
			{
				ip.setf(i, j, (float)kernel[i][j]);
			}
		}
		return ip;
	}
	public static int LINCOOR(int row, int col,int width)
	{ 
		return (row)*(width)+(col);				
	}
}
