
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
public class MaskFilter
{
	/**
	 *	Constants
	 */
	//public static final int DEFAULT_KERNEL_SIZE = 3;
	public static final boolean DEFAULT_FILTER_MODE_PASSTHROUGH = true; //
	
	/**
	 *	Constructor
	 */
	public MaskFilter()
	{
		/* do nothing */
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static ImagePlus run(ImagePlus imp, ImagePlus mask)
	{
		return run(imp, mask, DEFAULT_FILTER_MODE_PASSTHROUGH);
	}
	
	public static ImagePlus run(ImagePlus imp, ImagePlus mask, boolean filter_mode_passthrough)
	{
		// duplicate image
		ImagePlus imp_dup = imp;//.duplicate();
		
		// process stack
		ImageStack stack = imp_dup.getStack();
		for(int slice = 1; slice <= stack.getSize(); ++slice)
		{
			ImageProcessor slice_ip = stack.getProcessor(slice);
			ImageProcessor slice_result = run(slice_ip, mask.getProcessor(), filter_mode_passthrough);
			stack.setProcessor(slice_result, slice);
		}
		imp_dup.setStack(stack);
		
		// return image
		return imp_dup;
	}
	
	public static ImageProcessor run(ImageProcessor ip, ImageProcessor mask)
	{
		return run(ip, mask, DEFAULT_FILTER_MODE_PASSTHROUGH);
	}
	
	public static ImageProcessor run(ImageProcessor ip, ImageProcessor mask, boolean filter_mode_passthrough)
	{
		ImageProcessor filter_passthrough = ip.duplicate();
		ImageProcessor filter_block = ip.duplicate();
		for(int py = 0; py < ip.getHeight(); ++py)
		{
			for(int px = 0; px < ip.getWidth(); ++px)
			{
				if(mask.get(px, py) == 0)
				{
					filter_passthrough.set(px, py, 0);
				}
				else
				{
					filter_block.set(px, py, 0);
				}
			}
		}
		
		return (filter_mode_passthrough ? filter_passthrough : filter_block);
	}
}
