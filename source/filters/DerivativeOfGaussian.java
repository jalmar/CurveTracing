
package filters;

// import ImageJ classes
import ij.IJ.*;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;


/**
 *	
 */
public class DerivativeOfGaussian
{
	public static final double DEFAULT_SIGMA = 1.0;
	public static final double DEFAULT_THETA = 0.0;
	
	/**
	 *	Constructor
	 */
	public DerivativeOfGaussian()
	{
		/* do nothing */
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	/*public static ImageProcessor firstOrder(ImageProcessor ip)
	{
		return firstOrder(ip, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor firstOrder(ImageProcessor ip, double sigma)
	{
		return firstOrder(ip, sigma, DEFAULT_THETA);
	}
	
	public static ImageProcessor firstOrder(ImageProcessor ip, double sigma, double theta)
	{
		// convolve image with DoG_dx kernel
		double[][] kernel = computeKernelGaussian2D_du(sigma);
		return convolve(ip, kernel);
	}*/
	
	public static ImageProcessor derivativeX(ImageProcessor ip)
	{
		return derivativeX(ip, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor derivativeX(ImageProcessor ip, double sigma)
	{
		// convolve image with DoG_dx kernel
		double[][] kernel = computeKernelGaussian2D_du(sigma);
		return convolve(ip, kernel);
	}
	
	public static ImageProcessor derivativeY(ImageProcessor ip)
	{
		return derivativeY(ip, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor derivativeY(ImageProcessor ip, double sigma)
	{
		// convolve image with DoG_dy kernel
		double[][] kernel = computeKernelGaussian2D_dv(sigma);
		return convolve(ip, kernel);
	}
	
	public static ImageProcessor derivativeXX(ImageProcessor ip)
	{
		return derivativeXX(ip, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor derivativeXX(ImageProcessor ip, double sigma)
	{
		// convolve image with DoG_dx_dx kernel
		double[][] kernel = computeKernelGaussian2D_du_du(sigma);
		return convolve(ip, kernel);
	}
	
	public static ImageProcessor derivativeYY(ImageProcessor ip)
	{
		return derivativeYY(ip, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor derivativeYY(ImageProcessor ip, double sigma)
	{
		// convolve image with DoG_dy_dy kernel
		double[][] kernel = computeKernelGaussian2D_dv_dv(sigma);
		return convolve(ip, kernel);
	}
	
	public static ImageProcessor derivativeXY(ImageProcessor ip)
	{
		return derivativeXY(ip, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor derivativeXY(ImageProcessor ip, double sigma)
	{
		// convolve image with DoG_dx_dy kernel
		double[][] kernel = computeKernelGaussian2D_du_dv(sigma);
		return convolve(ip, kernel);
	}
	
	public static ImageProcessor derivativeYX(ImageProcessor ip)
	{
		return derivativeYX(ip, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor derivativeYX(ImageProcessor ip, double sigma)
	{
		// convolve image with DoG_dy_dx kernel
		double[][] kernel = computeKernelGaussian2D_dv_du(sigma);
		return convolve(ip, kernel);
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
	
	public static double[][] computeKernelGaussian2D(double sigma, double theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*sigma); // RSLV: use floor instead of round?
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
				kernel[kx][ky] = gaussian2D(u, v, sigma);
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
	
	public static double[][] computeKernelGaussian2D_du(double sigma, double theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*sigma); // RSLV: use floor instead of round?
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
				kernel[kx][ky] = gaussian2D_dx(u, v, sigma);
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
	
	public static double[][] computeKernelGaussian2D_dv(double sigma, double theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*sigma); // RSLV: use floor instead of round?
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
				kernel[kx][ky] = gaussian2D_dy(u, v, sigma);
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
	
	public static double[][] computeKernelGaussian2D_du_du(double sigma, double theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*sigma); // RSLV: use floor instead of round?
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
				kernel[kx][ky] = gaussian2D_dx_dx(u, v, sigma);
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
	
	public static double[][] computeKernelGaussian2D_dv_dv(double sigma, double theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*sigma); // RSLV: use floor instead of round?
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
				kernel[kx][ky] = gaussian2D_dy_dy(x, y, sigma);
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
	
	public static double[][] computeKernelGaussian2D_du_dv(double sigma, double theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*sigma); // RSLV: use floor instead of round?
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
				kernel[kx][ky] = gaussian2D_dx_dy(x, y, sigma);
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
	
	public static double[][] computeKernelGaussian2D_dv_du(double sigma, double theta)
	{
		// calculate required kernel size (2*3*sigma~99%)
		int kernel_radius = (int)Math.round(3*sigma); // RSLV: use floor instead of round?
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
				kernel[kx][ky] = gaussian2D_dy_dx(x, y, sigma);
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
		return (1/(2*Math.PI*Math.pow(sigma, 2)))*Math.exp(-(x*x+y*y)/(2*sigma*sigma));
	}
	
	public static double gaussian2D_dx(double x, double y)
	{
		return gaussian2D_dx(x, y, DEFAULT_SIGMA);
	}
	
	public static double gaussian2D_dx(double x, double y, double sigma)
	{
		return ((-x)/(2*Math.PI*Math.pow(sigma, 4)))*Math.exp(-(x*x+y*y)/(2*sigma*sigma));
	}
	
	public static double gaussian2D_dy(double x, double y)
	{
		return gaussian2D_dy(x, y, DEFAULT_SIGMA);
	}
	
	public static double gaussian2D_dy(double x, double y, double sigma)
	{
		return gaussian2D_dx(y, x, sigma); // NOTE x and y are swapped
	}
	
	public static double gaussian2D_dx_dx(double x, double y)
	{
		return gaussian2D_dx_dx(x, y, DEFAULT_SIGMA);
	}
	
	public static double gaussian2D_dx_dx(double x, double y, double sigma)
	{
		return ((x*x-sigma*sigma)/(2*Math.PI*Math.pow(sigma, 6)))*Math.exp(-(x*x+y*y)/(2*sigma*sigma));
	}
	
	public static double gaussian2D_dy_dy(double x, double y)
	{
		return gaussian2D_dy_dy(x, y, DEFAULT_SIGMA);
	}
	
	public static double gaussian2D_dy_dy(double x, double y, double sigma)
	{
		return gaussian2D_dx_dx(y, x, sigma); // NOTE x and y are swapped
	}
	
	public static double gaussian2D_dx_dy(double x, double y)
	{
		return gaussian2D_dx_dy(x, y, DEFAULT_SIGMA);
	}
	
	public static double gaussian2D_dx_dy(double x, double y, double sigma)
	{
		return ((x*y)/(2*Math.PI*Math.pow(sigma, 6)))*Math.exp(-(x*x+y*y)/(2*sigma*sigma));
	}
	
	public static double gaussian2D_dy_dx(double x, double y)
	{
		return gaussian2D_dy_dx(x, y, DEFAULT_SIGMA);
	}
	
	public static double gaussian2D_dy_dx(double x, double y, double sigma)
	{
		return gaussian2D_dx_dy(y, x, sigma); // NOTE x and y are swapped (although not necessary since dx_dy == dy_dx)
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
	
	public static ImageProcessor directionalDerivative(ImageProcessor ip, double direction)
	{
		return directionalDerivative(ip, direction, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor directionalDerivative(ImageProcessor ip, double theta, double sigma)
	{
		ImageProcessor ip_dx = derivativeX(ip, sigma);
		ImageProcessor ip_dy = derivativeY(ip, sigma);
		return directionalDerivative(ip_dx, ip_dy, theta);
	}
	
	public static ImageProcessor directionalDerivative(ImageProcessor ip_dx, ImageProcessor ip_dy, double theta)
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
}
