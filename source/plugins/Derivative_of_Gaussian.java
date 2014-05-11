
package plugins;

// import ImageJ classes
import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;

import ij.plugin.PlugIn;
import ij.gui.GenericDialog;

import ij.Prefs;

// import own classes
import filters.DerivativeOfGaussian;


/**
 *	
 */
public class Derivative_of_Gaussian implements PlugIn
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
		GenericDialog gd = new GenericDialog("Derivative of Gaussian filter");
		gd.addChoice("Derivative", new String[]{"d/dx", "d/dy", "d^2/dxdx", "d^2/dxdy", "d^2/dydx", "d^2/dydy"}, Prefs.get("DoG_filter.derivative_s", "d/dx"));
		gd.addNumericField("sigma(X)", Prefs.get("DoG_filter.sigma_x", 1.0), 1);
		gd.addNumericField("sigma(Y)", Prefs.get("DoG_filter.sigma_y", 1.0), 1);
		gd.addNumericField("orientation", Prefs.get("DoG_filter.orientation", 0.00), 2);
		
		gd.showDialog();
		if(gd.wasCanceled()) return;
		
		// retrieve parameters
		String derivative_s = gd.getNextChoice();
		double sigma_x = gd.getNextNumber();
		double sigma_y = gd.getNextNumber();
		double orientation = gd.getNextNumber();
		
		// save preferences
		Prefs.set("DoG_filter.derivative_s", derivative_s);
		Prefs.set("DoG_filter.sigma_x", sigma_x);
		Prefs.set("DoG_filter.sigma_y", sigma_y);
		Prefs.set("DoG_filter.orientation", orientation);
		
		// convert degrees to radians
		orientation = Math.toRadians(orientation);
		
		// execute filter
		ImagePlus result = exec(imp, derivative_s, sigma_x, sigma_y, orientation);
		
		// show image
		if(null != result)
		{
			result.show();
		}
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public ImagePlus exec(ImagePlus imp, String derivative, double sigma_x, double sigma_y, double orientation)
	{
		// TODO: check arguments
		return derivativeOfGaussian(imp, derivative, sigma_x, sigma_y, orientation);
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static ImagePlus derivativeOfGaussian(ImagePlus imp, String derivative, double sigma_x, double sigma_y, double orientation)
	{
		// forward FFT
		ImageProcessor ip = imp.getProcessor();
		if(derivative.equals("d/dx"))
		{
			ip = DerivativeOfGaussian.derivativeX(ip, sigma_x, sigma_y, (float)orientation);
		}
		else if(derivative.equals("d/dy"))
		{
			ip = DerivativeOfGaussian.derivativeY(ip, sigma_x, sigma_y, (float)orientation);
		}
		else if(derivative.equals("d^2/dxdx"))
		{
			ip = DerivativeOfGaussian.derivativeXX(ip, sigma_x, sigma_y, (float)orientation);
		}
		else if(derivative.equals("d^2/dxdy"))
		{
			ip = DerivativeOfGaussian.derivativeXY(ip, sigma_x, sigma_y, (float)orientation);
		}
		else if(derivative.equals("d^2/dydx"))
		{
			ip = DerivativeOfGaussian.derivativeYX(ip, sigma_x, sigma_y, (float)orientation);
		}
		else if(derivative.equals("d^2/dydy"))
		{
			ip = DerivativeOfGaussian.derivativeYY(ip, sigma_x, sigma_y, (float)orientation);
		}
		else 
		{
			// NOT supported
			IJ.error("Wrong derivative specified: " + derivative);
		}
		return new ImagePlus("DoG_" + imp.getTitle(), ip);
	}
}
