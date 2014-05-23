
package filters;

// import ImageJ classes
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import ij.process.LUT;

import ij.measure.ResultsTable;

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
	
	public static final int DEFAULT_COMPONENT_MIN_AREA_SIZE_THRESHOLD = 1;
	public static final int DEFAULT_COMPONENT_MAX_AREA_SIZE_THRESHOLD = Integer.MAX_VALUE;
	
	public static enum Connectivity { FOUR_CONNECTIVITY, EIGHT_CONNECTIVITY }; // RSLV: add constructor to allow FOUR_CONNECTIVITY[4], EIGHT_CONNECTIVITY[8]; OR drop enum and use integer {4, 8}
	public static final Connectivity DEFAULT_CONNECTIVITY = Connectivity.EIGHT_CONNECTIVITY;
	
	public static final LUT cc_lut = getConnectedComponentLUT();
	
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
	
	public static ImagePlus run(ImagePlus imp, Connectivity connectivity, int min_area, int max_area)
	{
		// duplicate image
		ImagePlus imp_dup = imp;//.duplicate();
		
		// process stack
		ImageStack stack = imp_dup.getStack();
		for(int slice = 1; slice <= stack.getSize(); ++slice)
		{
			ImageProcessor slice_ip = stack.getProcessor(slice);
			ImageProcessor slice_result = run(slice_ip, connectivity, min_area, max_area);
			stack.setProcessor(slice_result, slice);
		}
		imp_dup.setStack(stack);
		
		// return image
		return imp_dup;
	}
	
	public static ImageProcessor run(ImageProcessor ip, Connectivity connectivity, int min_area, int max_area)
	{
		// get image dimensions
		int img_width = ip.getWidth();
		int img_height = ip.getHeight();
		
		// duplicate image
		ImageProcessor ip_dup = ip.duplicate();
		
		// first pass: classify each pixel
		int label = 0;
		ImageProcessor labels_ip = new ShortProcessor(img_width, img_height);
//		ResultsTable rt = ResultsTable.getResultsTable();
//		rt.reset(); // clear table
		for(int py = 0; py < img_height; ++py)
		{
			for(int px = 0; px < img_width; ++px)
			{
				// get pixel value
				int pv = ip_dup.get(px, py);
				int lx = labels_ip.get(px, py);
				
				// check pixels
				if(pv == FOREGROUND && lx == BACKGROUND)
				{
					// label connected component with a new label value
					int area = labelRegion(ip_dup, labels_ip, px, py, ++label, connectivity);
					
					// filter too small and too large components (and indirectly relabels the components afterwards)
					if(area < min_area || area > max_area)
					{
						relabelRegion(ip_dup, px, py, FOREGROUND, BACKGROUND, connectivity); // remove component from duplicate input
						relabelRegion(labels_ip, px, py, label--, BACKGROUND, connectivity);
					}
//					else
//					{
//						// add area to ResultsTable
//						rt.incrementCounter();
//						rt.addValue("Label", label);
//						rt.addValue("Area", area);
//					}
				}
			}
		}
		
//		// show results table
//		rt.show("Results");
		
		// return labeled image
		labels_ip.setLut(cc_lut);
		return labels_ip;
	}
	
	/**
	 *	Recursively label connected pixels with label value
	 */
	public static int labelRegion(ImageProcessor ip, ImageProcessor labels_ip, int x, int y, int label, Connectivity connectivity)
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
			area += labelRegion(ip, labels_ip, x    , y - 1, label, connectivity); // 4
			area += labelRegion(ip, labels_ip, x - 1, y    , label, connectivity); // 4
			area += labelRegion(ip, labels_ip, x + 1, y    , label, connectivity); // 4
			area += labelRegion(ip, labels_ip, x    , y + 1, label, connectivity); // 4
			if(connectivity == Connectivity.EIGHT_CONNECTIVITY)
			{
				area += labelRegion(ip, labels_ip, x - 1, y - 1, label, connectivity); // 8
				area += labelRegion(ip, labels_ip, x + 1, y - 1, label, connectivity); // 8
				area += labelRegion(ip, labels_ip, x - 1, y + 1, label, connectivity); // 8
				area += labelRegion(ip, labels_ip, x + 1, y + 1, label, connectivity); // 8
			}
		}
		return area;
	}
	
	public static int relabelRegion(ImageProcessor ip, int x, int y, int find, int replace, Connectivity connectivity)
	{
		int area = 0;
		int pv = ip.getPixel(x, y); // NOTE: use getPixel to avoid out of bounds
		if(pv == find)
		{
			// replace and continue search in neighbourhood
			ip.set(x, y, replace);
			area += 1;
			
			// continue search to neighbours (recursive function call)
			area += relabelRegion(ip, x    , y - 1, find, replace, connectivity); // 4
			area += relabelRegion(ip, x - 1, y    , find, replace, connectivity); // 4
			area += relabelRegion(ip, x + 1, y    , find, replace, connectivity); // 4
			area += relabelRegion(ip, x    , y + 1, find, replace, connectivity); // 4
			if(connectivity == Connectivity.EIGHT_CONNECTIVITY)
			{
				area += relabelRegion(ip, x - 1, y - 1, find, replace, connectivity); // 8
				area += relabelRegion(ip, x + 1, y - 1, find, replace, connectivity); // 8
				area += relabelRegion(ip, x - 1, y + 1, find, replace, connectivity); // 8
				area += relabelRegion(ip, x + 1, y + 1, find, replace, connectivity); // 8
			}
		}
		return area;
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static LUT getConnectedComponentLUT()
	{
		byte[] cc_lut_r = new byte[256];
		byte[] cc_lut_g = new byte[256];
		byte[] cc_lut_b = new byte[256];
		
		// 254 / 6 = 42.333
		// R: 255  -   0   0   +  255
		// G:  +  255 255  -   0   0
		// B:  0   0   +  255 255  -
		byte r = (byte)255;
		byte g = (byte)0;
		byte b = (byte)0;
		for(int i = 1; i < 256; ++i)
		{
			int section = (int)(i / 42.333);
			//int index = (int)(section % 42.333);
			switch(section)
			{
				case 0:
					r = (byte)255;
					g += 6;
					b = (byte)0;
					break;
				case 1:
					r -= 6;
					g = (byte)255;
					b = (byte)0;
					break;
				case 2:
					r = (byte)0;
					g = (byte)255;
					b += 6;
					break;
				case 3:
					r = (byte)0;
					g -= 6;
					b = (byte)255;
					break;
				case 4:
					r += 6;
					g = (byte)0;
					b = (byte)255;
					break;
				case 5:
					r = (byte)255;
					g = (byte)0;
					b -= 6;
					break;
			}
			cc_lut_r[i] = r;
			cc_lut_g[i] = g;
			cc_lut_b[i] = b;
		}
		
		return new LUT(cc_lut_r, cc_lut_g, cc_lut_b);
	}
}
