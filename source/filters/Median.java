
package filters;

// import ImageJ classes
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

/**
 *	
 */
public class Median
{
	/**
	 *	Constants
	 */
	public static final int DEFAULT_KERNEL_SIZE = 3;
	
	/**
	 *	Constructor
	 */
	public Median()
	{
		/* do nothing */
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static ImagePlus run(ImagePlus imp)
	{
		return run(imp, DEFAULT_KERNEL_SIZE);
	}
	
	public static ImagePlus run(ImagePlus imp, int size)
	{
		// duplicate image
		ImagePlus imp_dup = imp;//.duplicate();
		
		// process stack
		ImageStack stack = imp_dup.getStack();
		for(int slice = 1; slice <= stack.getSize(); ++slice)
		{
			ImageProcessor slice_ip = stack.getProcessor(slice);
			ImageProcessor slice_result = run(slice_ip, size);
			stack.setProcessor(slice_result, slice);
		}
		imp_dup.setStack(stack);
		
		// return image
		return imp_dup;
	}
	
	public static ImageProcessor run(ImageProcessor ip)
	{
		return run(ip, DEFAULT_KERNEL_SIZE);
	}
	
	public static ImageProcessor run(ImageProcessor ip, int size)
	{
		// assert size is odd, size > 1
		int half_size = (int)Math.round(0.5*size);
		//int square_size = size*size;
		ImageProcessor ip_out = new FloatProcessor(ip.getWidth(), ip.getHeight());
		for(int py = 0; py < ip_out.getHeight(); ++py)
		{
			for(int px = 0; px < ip_out.getWidth(); ++px)
			{
				//double[] values = new double[square_size];
				java.util.Vector<Double> values = new java.util.Vector<Double>();
				for(int ky = -half_size; ky <= half_size; ++ky)
				{
					for(int kx = -half_size; kx <= half_size; ++kx)
					{
						int y = py + ky;
						int x = px + kx;
						if(y >= 0 && y < ip.getHeight() && x >= 0 && x < ip.getWidth())
						{
							values.add(new Double(ip.getf(x, y)));
						}
					}
				}
				Double[] values_arr = values.toArray(new Double[]{});
				java.util.Arrays.sort(values_arr);
				//System.err.println(java.util.Arrays.toString(values_arr));
				double median_value = values_arr[(int)Math.round(0.5*values_arr.length)];
				ip_out.setf(px, py, (float)median_value);
			}
		}
		return ip_out;
	}
}
