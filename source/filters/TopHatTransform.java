
package filters;

// import ImageJ classes
import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

/**
 *	
 */
public class TopHatTransform
{
	/**
	 *
	 */
	public static final int DEFAULT_KERNEL_RADIUS = 3;
	public static final int KERNEL_ON = 255;
	public static final int KERNEL_OFF = 0;
	
	/**
	 *	Constructor
	 */
	public TopHatTransform()
	{
		/* do nothing */
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static ImageProcessor run(ImageProcessor ip)
	{
		return run(ip, DEFAULT_KERNEL_RADIUS);
	}
	
	public static ImageProcessor run(ImageProcessor ip, boolean whiteTopHatTransform)
	{
		return run(ip, DEFAULT_KERNEL_RADIUS, whiteTopHatTransform);
	}
	
	public static ImageProcessor run(ImageProcessor ip, double sigma)
	{
		return run(ip, sigma, true);
	}
	
	public static ImageProcessor run(ImageProcessor ip, double sigma, boolean whiteTopHatTransform)
	{
		// construct structure element (disc kernel with radius 3*sigma)
		int kernel_radius = (int)Math.round(3*sigma);
		int kernel_size = 1+2*kernel_radius;
		int[][] se_kernel = createDiscKernel(kernel_radius);
		
		// white top hat := f - opening(f,se); opening -> min-max
		// black top hat := f - closing(f,se); closing -> max-min
		ImageProcessor result1 = new FloatProcessor(ip.getWidth(), ip.getHeight());
		for(int py = 0; py < ip.getHeight(); ++py)
		{
			for(int px = 0; px < ip.getWidth(); ++px)
			{
				double new_value = (whiteTopHatTransform ? Double.MAX_VALUE : Double.MIN_VALUE);
				for(int kx = 0; kx < se_kernel.length; ++kx)
				{
					int x = px - kernel_radius + kx;
					for(int ky = 0; ky < se_kernel.length; ++ky)
					{
						int y = py - kernel_radius + ky;
						if(se_kernel[kx][ky] == KERNEL_ON && y >= 0 && y < ip.getHeight() && x >= 0 && x < ip.getWidth())
						{
							new_value = (whiteTopHatTransform ? Math.min(new_value, ip.getf(x, y)) : Math.max(new_value, ip.getf(x, y)));
						}
					}
				}
				result1.setf(px, py, (float)new_value);
			}
		}
		
		ImageProcessor result2 = new FloatProcessor(ip.getWidth(), ip.getHeight());
		for(int py = 0; py < result1.getHeight(); ++py)
		{
			for(int px = 0; px < result1.getWidth(); ++px)
			{
				double new_value = (whiteTopHatTransform ? Double.MIN_VALUE : Double.MAX_VALUE);
				for(int kx = 0; kx < se_kernel.length; ++kx)
				{
					int x = px - kernel_radius + kx;
					for(int ky = 0; ky < se_kernel.length; ++ky)
					{
						int y = py - kernel_radius + ky;
						if(se_kernel[kx][ky] == KERNEL_ON && y >= 0 && y < result1.getHeight() && x >= 0 && x < result1.getWidth())
						{
							new_value = (whiteTopHatTransform ? Math.max(new_value, result1.getf(x, y)) : Math.min(new_value, result1.getf(x, y)));
						}
					}
				}
				result2.setf(px, py, (float)new_value);
			}
		}
		
		// return result
		return result2;
	}
	
	public static int[][] createDiscKernel()
	{
		return createDiscKernel(1);
	}
	
	public static int[][] createDiscKernel(int kernel_radius)
	{
		int kernel_size = 1+2*kernel_radius;
		int kernel[][] = new int[kernel_size][kernel_size];
		
		for(int ky = 0; ky < kernel_size; ++ky)
		{
			int y = ky - kernel_radius;
			for(int kx = 0; kx < kernel_size; ++kx)
			{
				int x = kx - kernel_radius;
				if(Math.sqrt(x*x + y*y) <= kernel_radius)
				{
					kernel[kx][ky] = KERNEL_ON;
				}
				//else
				//{
				//	kernel[kx][ky] = KERNEL_OFF;
				//}
			}
		}
		return kernel;
	}
}
