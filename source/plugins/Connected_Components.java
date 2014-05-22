
package plugins;

// import ImageJ classes
import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;

import ij.plugin.PlugIn;
import ij.gui.GenericDialog;

import ij.Prefs;

// import own classes
import filters.ConnectedComponents;


/**
 *	
 */
public class Connected_Components implements PlugIn
{
	/**
	 *	Constants
	 */
	public static final int MASK_ON = 255;
	public static final int MASK_OFF = 0;
	
	public static final int DEFAULT_VOTING_THRESHOLD = 5;
	
	/**
	 *	Constructor
	 */
	//public Fill_Holes()
	//{
	//	/* do nothing */
	//}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public void run(String arg)
	{
		// get active image
		ImagePlus imp = IJ.getImage();
		if(null == imp) return;
		
		// ask for parameters; TODO: add 4 or 8 connectivity choice
		/*GenericDialog gd = new GenericDialog("Fill Holes filter");
		gd.addNumericField("voting threshold", Prefs.get("Fill_Holes_filter.voting_threshold", DEFAULT_VOTING_THRESHOLD), 0);
		
		gd.showDialog();
		if(gd.wasCanceled()) return;
		
		// retrieve parameters
		int voting_threshold = (int)gd.getNextNumber();
		
		// save preferences
		Prefs.set("Fill_Holes_filter.voting_threshold", voting_threshold);
		*/
		
		// execute filter
		ImagePlus result = exec(imp);
		
		// show image
		if(null != result)
		{
			result.show();
		}
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public ImagePlus exec(ImagePlus imp)
	{
		// TODO: check arguments
		return connectedComponents(imp);
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static ImagePlus connectedComponents(ImagePlus imp)
	{
		ImageProcessor ip = imp.getProcessor();
		ip = ConnectedComponents.run(ip);
		return new ImagePlus("Connected_Components_" + imp.getTitle(), ip);
	}
}
