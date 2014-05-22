
package filters;

// import ImageJ classes
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

/**
 *
 */
public class ConnectedComponents
{
	/**
	 *	Constants
	 */
	public static final int BACKGROUND = 0; // NOTE: best if kept 0
	public static final int FOREGROUND = 255;
	
	/**
	 *	Members
	 */
	//public double[][] 
	
	/**
	 *	Constructor
	 */
	public ConnectedComponents()
	{
		/* empty */
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static ImagePlus run(ImagePlus imp, int connectivity)
	{
		// duplicate image
		ImagePlus imp_dup = imp;//.duplicate();
		
		// process stack
		ImageStack stack = imp_dup.getStack();
		for(int slice = 1; slice <= stack.getSize(); ++slice)
		{
			ImageProcessor slice_ip = stack.getProcessor(slice);
			ImageProcessor slice_result = run(slice_ip);
			stack.setProcessor(slice_result, slice);
		}
		imp_dup.setStack(stack);
		
		// return image
		return imp_dup;
	}
	
	public static ImageProcessor run(ImageProcessor ip)
	{
		int img_width = ip.getWidth();
		int img_height = ip.getHeight();
		
		// first pass: classify each pixel
		int label = 0;
		//int pixel_count_label = 0;
		ImageProcessor labels_ip = new ShortProcessor(img_width, img_height);
		for(int py = 0; py < img_height; ++py)
		{
			for(int px = 0; px < img_width; ++px)
			{
				// get pixel value
				int pv = ip.get(px, py);
				int lx = labels_ip.get(px, py);
				
				// check pixels
				if(pv == FOREGROUND && lx == BACKGROUND)
				{
					// start new label
					label += 1;
					
					int area = labelRegion(ip, labels_ip, px, py, label);
					System.err.println("label=" + label + ", area=" + area);
					
					/*
					// get neighbouring pixel values
					int lp = class_labels_ip.getPixel(px - 1, py - 1);
					int lq = class_labels_ip.getPixel(px    , py - 1);
					int lr = class_labels_ip.getPixel(px + 1, py - 1);
					int ls = class_labels_ip.getPixel(px - 1, py    );
				
					// check if all neighbours are background => new label
					if(lp == BACKGROUND && lq == BACKGROUND && lr == BACKGROUND && ls == BACKGROUND)
					{
						// assign new label
						class_labels_ip.set(px, py, ++last_label);
						//++pixel_count_label;
					}
					else
					{
						// TODO: check label of neighbours
						if(ls != BACKGROUND)
						{
							class_labels_ip.set(px, py, ls);
						}
						else if(lr != BACKGROUND)
						{
							class_labels_ip.set(px, py, lr);
						}
						else if(lq != BACKGROUND)
						{
							class_labels_ip.set(px, py, lq);
						}
						else if(lp != BACKGROUND)
						{
							class_labels_ip.set(px, py, lp);
						}
						else
						{
							// RSLV: panic!! this should never happen
						}
						
						// TODO: resolve conflicts; i.e. register class equivalence
						
					}
					*/
				}
			}
		}
		
		return labels_ip;
	}
	
	public static int labelRegion(ImageProcessor ip, ImageProcessor labels_ip, int x, int y, int label)
	{
		int area = 0;
		int pv = ip.getPixel(x, y);
		int pl = labels_ip.getPixel(x, y);
		
		if(pv == FOREGROUND && pl == BACKGROUND) // NOTE: latter condition prevent the algorithm from retracking on itself
		{
			// add pixel to region
			labels_ip.set(x, y, label);
			area += 1;
			
			// continue search to neighbours (recursive function call)
			area += labelRegion(ip, labels_ip, x - 1, y - 1, label); // 8
			area += labelRegion(ip, labels_ip, x    , y - 1, label); // 4
			area += labelRegion(ip, labels_ip, x + 1, y - 1, label); // 8
			area += labelRegion(ip, labels_ip, x - 1, y    , label); // 4
			area += labelRegion(ip, labels_ip, x + 1, y    , label); // 4
			area += labelRegion(ip, labels_ip, x - 1, y + 1, label); // 8
			area += labelRegion(ip, labels_ip, x    , y + 1, label); // 4
			area += labelRegion(ip, labels_ip, x + 1, y + 1, label); // 8
		}
		return area;
	}
	
	/*
	public static int relabelRegion(ImageProcessor ip, int x, int y, int find, int replace)
	{
		int area = 0;
		int pv = ip.get(x, y);
		if(pv == find)
		{
			// replace and continue search in neighbourhood
			ip.set(x, y, replace);
			area += 1;
			area += relabelRegion(ip, x - 1, y - 1, find, replace); // 8
			area += relabelRegion(ip, x    , y - 1, find, replace); // 4
			area += relabelRegion(ip, x + 1, y - 1, find, replace); // 8
			area += relabelRegion(ip, x - 1, y    , find, replace); // 4
			area += relabelRegion(ip, x + 1, y    , find, replace); // 4
			area += relabelRegion(ip, x - 1, y + 1, find, replace); // 8
			area += relabelRegion(ip, x    , y + 1, find, replace); // 4
			area += relabelRegion(ip, x + 1, y + 1, find, replace); // 8
		}
		return area;
	}*/
}
