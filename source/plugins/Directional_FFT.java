
package plugins;

// import ImageJ classes
import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.ByteProcessor;
import ij.process.FHT;

import ij.plugin.PlugIn;
import ij.gui.GenericDialog;

/**
 *	
 */
public class Directional_FFT implements PlugIn
{
	/**
	 *	Constants
	 */
	public static final int MASK_ON = 255;
	public static final int MASK_OFF = 0;
	
	/**
	 *	Constructor
	 */
	//public Directional_FFT()
	//{
	//	/* do nothing */
	//}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public void run(String arg)
	{
		// get active image
		ImagePlus imp = IJ.getImage();
		if(null == imp) return;
		
		// ask for parameters
		GenericDialog gd = new GenericDialog("Directional FFT filter");
		gd.addNumericField("min theta:", 0, 0);
		gd.addNumericField("max theta:", 0, 0);
		gd.showDialog();
		if(gd.wasCanceled()) return;
		
		// retrieve parameters
		double min_theta = gd.getNextNumber();
		double max_theta = gd.getNextNumber();
		
		// convert degrees to radians
		min_theta = Math.toRadians(min_theta);
		max_theta = Math.toRadians(max_theta);
		
		// swap min/max if necessary
		/*if(min_theta > max_theta)
		{
			double tmp = min_theta;
			min_theta = max_theta;
			max_theta = tmp;
		}*/
		
		System.err.println(min_theta);
		System.err.println(max_theta);
		
		// execute filter
		ImagePlus result = exec(imp, min_theta, max_theta);
		
		// show image
		if(null != result)
		{
			result.show();
		}
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public ImagePlus exec(ImagePlus imp, double min_theta, double max_theta)
	{
		// TODO: check arguments
		return directionalFFT(imp, min_theta, max_theta);
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static ImagePlus directionalFFT(ImagePlus imp, double min_theta, double max_theta)
	{
		// forward FFT
		ImageProcessor ip = imp.getProcessor();
		FHT fourier = new FHT(ip);
		fourier.transform();
		fourier = maskFourierImage(fourier, min_theta, max_theta);
		fourier.inverseTransform();
		// TODO: convert inverse transformed to original type (8, 16, or 32 bit)
		fourier.resetMinAndMax();
		return new ImagePlus(imp.getTitle(), fourier);
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	private static FHT maskFourierImage(FHT ip, double min_theta, double max_theta)
	{
		// RSLV: reverse angular direction from clockwise to counter-clockwise?
		// min_theta = -min_theta;
		// max_theta = -max_theta;
		
		// helper variables
		int width = ip.getWidth();
		int height = ip.getHeight();
		double center_x = 0.5 * width;
		double center_y = 0.5 * height;
		
		double min_theta_unit_vector_x = Math.cos(min_theta);
		double min_theta_unit_vector_y = Math.sin(min_theta);
		double max_theta_unit_vector_x = Math.cos(max_theta);
		double max_theta_unit_vector_y = Math.sin(max_theta);
		
		// mask fourier image
		for(int py = 0; py < height; ++py)
		{
			for(int px = 0; px < width; ++px)
			{
				// calculate pixel vector from center
				double pixel_vector_x = px - center_x;
				double pixel_vector_y = py - center_y;
				
				// check if pixel is within band using cross product
				double cross_pixel_min_theta = pixel_vector_x*min_theta_unit_vector_y - min_theta_unit_vector_x*pixel_vector_y;
				double cross_pixel_max_theta = pixel_vector_x*max_theta_unit_vector_y - max_theta_unit_vector_x*pixel_vector_y;
				if(Math.signum(cross_pixel_min_theta) == Math.signum(cross_pixel_max_theta) && cross_pixel_min_theta != 0) // outside mask (exclude zero)
				{
					ip.set(px, py, MASK_OFF);
				}
			}
		}
		// return the mask
		return ip;
	}
}
