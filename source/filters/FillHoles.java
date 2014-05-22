
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
public class FillHoles
{
	/**
	 *	Constants
	 */
	public static final int BACKGROUND = 0;
	public static final int FOREGROUND = 255;
	
	public static final int DEFAULT_VOTING_THRESHOLD = 5; // x out of 8
	
	/**
	 *	Constructor
	 */
	public FillHoles()
	{
		/* do nothing */
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static ImagePlus run(ImagePlus imp)
	{
		return run(imp, DEFAULT_VOTING_THRESHOLD);
	}
	
	public static ImagePlus run(ImagePlus imp, int voting_threshold)
	{
		// duplicate image
		ImagePlus imp_dup = imp;//.duplicate();
		
		// process stack
		ImageStack stack = imp_dup.getStack();
		for(int slice = 1; slice <= stack.getSize(); ++slice)
		{
			ImageProcessor slice_ip = stack.getProcessor(slice);
			ImageProcessor slice_result = run(slice_ip, voting_threshold);
			stack.setProcessor(slice_result, slice);
		}
		imp_dup.setStack(stack);
		
		// return image
		return imp_dup;
	}
	
	public static ImageProcessor run(ImageProcessor ip)
	{
		return run(ip, DEFAULT_VOTING_THRESHOLD);
	}
	
	public static ImageProcessor run(ImageProcessor ip, int voting_threshold)
	{
		int image_width = ip.getWidth();
		int image_height = ip.getHeight();
		boolean image_changed = true; // convergence test
		while(image_changed)
		{
			image_changed = false;
			for(int py = 0; py < image_height; ++py)
			{
				for(int px = 0; px < image_width; ++px)
				{
					if(ip.get(px, py) == BACKGROUND)
					{
						// check neighbourhood
						int voting_count = 0;
						
						// count neighbourhood votes
						for(int ky = -1; ky <= 1; ++ky)
						{
							for(int kx = -1; kx <= 1; ++kx)
							{
								if(ip.getPixel(px+kx, py+ky) == FOREGROUND)
								{
									++voting_count;
								}
							}
						}
						
						// update image
						if(voting_count >= voting_threshold)
						{
							ip.set(px, py, FOREGROUND);
							image_changed = true;
						}
					}
				}
			}
		}
		
		return ip;
	}
}
