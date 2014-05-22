
package plugins;

// import ImageJ classes
import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;

import ij.measure.ResultsTable;

import ij.plugin.PlugIn;
import ij.gui.GenericDialog;

import ij.Prefs;

// import own classes
import filters.MaskFilter;


/**
 *	
 */
public class NonZero_Statistics implements PlugIn
{
	/**
	 *	Constants
	 */
	public static enum Measures {AREA, MIN, MAX, MEAN, STDEV};
	public static final int MEASURE_COUNT = Measures.values().length;
	public static final ResultsTable results_table = ResultsTable.getResultsTable(); //new ResultsTable();
	
	/**
	 *	Constructor
	 */
	//public NonZero_Statistics()
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
		//GenericDialog gd = new GenericDialog("Mask filter");
		//gd.addChoice("Source image", image_names, image_names[0]);
		
		//gd.showDialog();
		//if(gd.wasCanceled()) return;
		
		// retrieve parameters
		//int source_img_index = gd.getNextChoiceIndex();
		
		// save preferences
		//Prefs.set("Mask_filter.filter_mode_passthrough", filter_mode_passthrough);
		
		// execute filter
		//ImagePlus result = exec(source_img, mask_img, filter_mode_passthrough);
		double[] results = nonZeroStatistics(imp);
		
		// show image
		//if(null != result)
		//{
		//	result.show();
		//}
		ResultsTable rt = results_table;
		rt.reset();
		rt.incrementCounter();
		rt.addValue("Image", imp.getTitle());
		rt.addValue("Area", results[Measures.AREA.ordinal()]);
		rt.addValue("Min", results[Measures.MIN.ordinal()]);
		rt.addValue("Max", results[Measures.MAX.ordinal()]);
		rt.addValue("Mean", results[Measures.MEAN.ordinal()]);
		rt.addValue("Stdev", results[Measures.STDEV.ordinal()]);
		//rt.show("Non-zero Statistics Results");
		rt.show("Results");
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public double[] exec(ImagePlus imp)
	{
		// TODO: check arguments
		return nonZeroStatistics(imp);
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static double[] nonZeroStatistics(ImagePlus imp)
	{
		ImageProcessor ip = imp.getProcessor();
		double[] results = new double[MEASURE_COUNT]; // AREA, MIN, MAX, MEAN, STDEV
		results[Measures.MIN.ordinal()] = 65535;
		
		for(int py = 0; py < ip.getHeight(); ++py)
		{
			for(int px = 0; px < ip.getWidth(); ++px)
			{
				int pv = ip.get(px, py);
				if(pv != 0)
				{
					results[Measures.AREA.ordinal()] += 1;
					results[Measures.MIN.ordinal()] = Math.min(results[Measures.MIN.ordinal()], pv);
					results[Measures.MAX.ordinal()] = Math.max(results[Measures.MAX.ordinal()], pv);
					results[Measures.MEAN.ordinal()] += pv;
					results[Measures.STDEV.ordinal()] += pv*pv;
				}
			}
		}
		
		// calculate actual mean and stdev measures
		if(results[Measures.AREA.ordinal()] > 0)
		{
			results[Measures.MEAN.ordinal()] = results[Measures.MEAN.ordinal()] / results[Measures.AREA.ordinal()];
			results[Measures.STDEV.ordinal()] = Math.sqrt((results[Measures.STDEV.ordinal()] / results[Measures.AREA.ordinal()]) - (results[Measures.MEAN.ordinal()] * results[Measures.MEAN.ordinal()]));
		}
		
		// return results
		return results;
	}
}
