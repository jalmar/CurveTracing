
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
	//public static final int BACKGROUND = 0; // NOTE: best if kept 0
	//public static final int FOREGROUND = 255;
	
	public static final String DEFAULT_CONNECTIVITY_S = "8-connectivity";
	
	public static final ConnectedComponents.Connectivity DEFAULT_CONNECTIVITY = ConnectedComponents.Connectivity.EIGHT_CONNECTIVITY;
	
	public static final int DEFAULT_COMPONENT_MIN_AREA_SIZE_THRESHOLD = 1;
	public static final int DEFAULT_COMPONENT_MAX_AREA_SIZE_THRESHOLD = Integer.MAX_VALUE;
	
	/**
	 *	Constructor
	 */
	//public Connected_Components()
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
		GenericDialog gd = new GenericDialog("Connected Components filter");
		gd.addChoice("Connectivity", new String[]{"4-connectivity", "8-connectivity"}, Prefs.get("Connected_Components_filter.connectivity_s", DEFAULT_CONNECTIVITY_S));
		
		gd.addNumericField("Minimum_component_area_size_threshold", Prefs.get("Connected_Components_filter.min_area_threshold", DEFAULT_COMPONENT_MIN_AREA_SIZE_THRESHOLD), 0);
		gd.addNumericField("Maximum_component_area_size_threshold", Prefs.get("Connected_Components_filter.max_area_threshold", DEFAULT_COMPONENT_MAX_AREA_SIZE_THRESHOLD), 0);
		
		gd.showDialog();
		if(gd.wasCanceled()) return;
		
		// retrieve parameters
		String connectivity_s = gd.getNextChoice();
		
		int min_area = (int)gd.getNextNumber();
		int max_area = (int)gd.getNextNumber();
		
		ConnectedComponents.Connectivity connectivity = DEFAULT_CONNECTIVITY;
		if(connectivity_s.equals("4-connectivity"))
		{
			connectivity = ConnectedComponents.Connectivity.FOUR_CONNECTIVITY;
		}
		else if(connectivity_s.equals("8-connectivity"))
		{
			connectivity = ConnectedComponents.Connectivity.EIGHT_CONNECTIVITY;
		}
		
		// save preferences
		Prefs.set("Connected_Components_filter.connectivity_s", connectivity_s);
		Prefs.set("Connected_Components_filter.min_area_threshold", min_area);
		Prefs.set("Connected_Components_filter.max_area_threshold", max_area);
		
		// execute filter
		ImagePlus result = exec(imp, connectivity, min_area, max_area);
		
		// show image
		if(null != result)
		{
			result.show();
		}
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public ImagePlus exec(ImagePlus imp, ConnectedComponents.Connectivity connectivity, int min_area, int max_area)
	{
		// TODO: check arguments
		return connectedComponents(imp, connectivity, min_area, max_area);
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static ImagePlus connectedComponents(ImagePlus imp, ConnectedComponents.Connectivity connectivity, int min_area, int max_area)
	{
		ImageProcessor ip = imp.getProcessor();
		ip = ConnectedComponents.run(ip, connectivity, min_area, max_area);
		return new ImagePlus("Connected_Components_" + imp.getTitle(), ip);
	}
}
