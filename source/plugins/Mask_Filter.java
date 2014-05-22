
package plugins;

// import ImageJ classes
import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;

import ij.WindowManager;

import ij.plugin.PlugIn;
import ij.gui.GenericDialog;

import ij.Prefs;

// import own classes
import filters.MaskFilter;


/**
 *	
 */
public class Mask_Filter implements PlugIn
{
	/**
	 *	Constants
	 */
	public static final int MASK_ON = 255;
	public static final int MASK_OFF = 0;
	
	public static final boolean DEFAULT_FILTER_MODE_PASSTHROUGH = true;
	
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
		// get list of open images
		int[] image_ids = WindowManager.getIDList();
		int image_count = WindowManager.getImageCount();
		String[] image_names = new String[image_count];
		for(int i = 0; i < image_count; ++i)
		{
			int image_id = image_ids[i];
			ImagePlus img = WindowManager.getImage(image_id);
			image_names[i] = img.getTitle();
		}
		
		if(image_count < 2) return; // requires at least two open images
		
		// ask for parameters
		GenericDialog gd = new GenericDialog("Mask filter");
		gd.addChoice("Source image", image_names, image_names[0]);
		gd.addChoice("Mask image", image_names, image_names[0]);
		gd.addCheckbox("Passthrough", Prefs.get("Mask_filter.passthrough", DEFAULT_FILTER_MODE_PASSTHROUGH));
		
		gd.showDialog();
		if(gd.wasCanceled()) return;
		
		// retrieve parameters
		//String source_img_title = gd.getNextChoice();
		//String mask_img_title = gd.getNextChoice();
		int source_img_index = gd.getNextChoiceIndex();
		int mask_img_index = gd.getNextChoiceIndex();
		boolean filter_mode_passthrough = gd.getNextBoolean();
		
		// save preferences
		Prefs.set("Mask_filter.filter_mode_passthrough", filter_mode_passthrough);
		
		// execute filter
		ImagePlus source_img = WindowManager.getImage(image_ids[source_img_index]);
		ImagePlus mask_img = WindowManager.getImage(image_ids[mask_img_index]);
		ImagePlus result = exec(source_img, mask_img, filter_mode_passthrough);
		
		// show image
		if(null != result)
		{
			result.show();
		}
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public ImagePlus exec(ImagePlus imp, ImagePlus mask, boolean filter_mode_passthrough)
	{
		// TODO: check arguments
		return maskFilter(imp, mask, filter_mode_passthrough);
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static ImagePlus maskFilter(ImagePlus imp, ImagePlus mask, boolean filter_mode_passthrough)
	{
		ImageProcessor source_ip = imp.getProcessor();
		ImageProcessor mask_ip = mask.getProcessor();
		ImageProcessor result_ip = MaskFilter.run(source_ip, mask_ip, filter_mode_passthrough);
		return new ImagePlus("Mask_Filter_" + imp.getTitle(), result_ip);
	}
}
