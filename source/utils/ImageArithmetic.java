
package utils;

// import ImageJ classes
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

/**
 *	
 */
public class ImageArithmetic
{
	/**
	 *	Constructor
	 */
	public ImageArithmetic()
	{
		/* do nothing */
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static ImageProcessor add(ImageProcessor lh, ImageProcessor rh)
	{
		// assert lh.getWidth() == rh.getWidth()
		// assert lh.getHeight() == rh.getHeight()
		ImageProcessor res = new FloatProcessor(lh.getWidth(), lh.getHeight());
		for(int py = 0; py < res.getHeight(); ++py)
		{
			for(int px = 0; px < res.getWidth(); ++px)
			{
				res.setf(px, py, (float)(lh.getf(px, py) + rh.getf(px, py)));
			}
		}
		return res;
	}
	
	public static ImageProcessor subtract(ImageProcessor lh, ImageProcessor rh)
	{
		// assert lh.getWidth() == rh.getWidth()
		// assert lh.getHeight() == rh.getHeight()
		ImageProcessor res = new FloatProcessor(lh.getWidth(), lh.getHeight());
		for(int py = 0; py < res.getHeight(); ++py)
		{
			for(int px = 0; px < res.getWidth(); ++px)
			{
				res.setf(px, py, (float)(lh.getf(px, py) - rh.getf(px, py)));
			}
		}
		return res;
	}
}
