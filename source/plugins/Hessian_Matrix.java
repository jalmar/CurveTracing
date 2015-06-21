
package plugins;

// import Java classes
import java.util.Collections;
import java.util.Comparator;
import java.util.Vector;
import java.util.Map;
import java.util.HashMap;
import java.util.HashSet;

import java.awt.Color;

// import ImageJ classes
import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.LUT;

import ij.plugin.PlugIn;
import ij.gui.GenericDialog;
import ij.gui.DialogListener;
import ij.Prefs;

import ij.plugin.frame.RoiManager;

import ij.gui.Overlay;
import ij.gui.Roi;
//import ij.gui.Line; // LineRoi // NOTE: will conflict with core.Line
import ij.gui.OvalRoi;
import ij.gui.PolygonRoi; // for POLYLINE

import ij.measure.ResultsTable;

// import Jama classes
import Jama.Matrix;
import Jama.EigenvalueDecomposition;

// import own classes
//import algorithms.Stegers;
import algorithms.LevenbergMarquardt;
import filters.DerivativeOfGaussian;
import filters.ConnectedComponents; // for rainbow LUT
import core.Point;
import core.Line;
import utils.Profiling;
import utils.Tuple;


/**
 *	
 */
public class Hessian_Matrix implements PlugIn
{
	/**
	 *	Parameters and constants
	 */
	public static double SIGMA = 1.8;
	public static double LINEPOINT_THRESHOLD = 0.5;
	
	public static boolean USE_GAUSSIAN_FIT_PEAKS_BEFORE = false;
	public static boolean USE_GAUSSIAN_FIT_PEAKS_AFTER = false;
	public static final int LMA_NUM_ITERATIONS = 5;
	public static final double LMA_DEFAULT_LAMBDA = 0.001;
	
	public static boolean USE_PRESET_USER_THRESHOLDS = false;
	public static double UPPER_THRESHOLD = 5000.0;
	public static double LOWER_THRESHOLD = 1000.0;
	
	public static boolean FILTER_MULTIRESPONSE = true;
	public static double MULTIRESPONSE_THRESHOLD = 360.0 / 18;
	
	public static boolean PREVENT_SPLINTER_FRAGMENTS = true;
	public static int SPLINTER_FRAGMENT_THRESHOLD = 5;
	
	public static boolean INCLUDE_JUNCTIONS_IN_LINE = true;
	public static boolean INCLUDE_FIRST_PROCESSED_IN_LINE = true;
	
	public static double COST_FUNCTION_WEIGHT = 1.0;
	
	public static boolean ENABLE_LINKING = true;
	public static double MAX_BENDING_ANGLE = 30.0;
	
	public static boolean FILTER_SHORT_LINE_SEGMENTS = true;
	public static int LINE_LENGTH_THRESHOLD = 10;
	
	public static boolean DEBUG_MODE_ENABLED = false;
	
	// ////////////////////////////////////////////////////////////////////////
	
	/**
	 *	Constructor
	 */
	public Hessian_Matrix()
	{
		/* do nothing */
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public void run(String arg)
	{
		// get active image
		ImagePlus imp = IJ.getImage();
		if(null == imp) return;
		
		// ask for parameters
		GenericDialog gd = new GenericDialog("Steger's algorithm");
		
		gd.addNumericField("Sigma", Prefs.get("stegers.sigma", SIGMA), 2);
		gd.addNumericField("Linepoint_threshold", Prefs.get("stegers.linepoint_threshold", LINEPOINT_THRESHOLD), 3);
		gd.addCheckbox("Use_gaussian_fit_peaks_before", Prefs.get("stegers.use_gaussian_fit_peaks_before", USE_GAUSSIAN_FIT_PEAKS_BEFORE));
		gd.addCheckbox("Use_gaussian_fit_peaks_after", Prefs.get("stegers.use_gaussian_fit_peaks_after", USE_GAUSSIAN_FIT_PEAKS_AFTER));
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Use_preset_user_thresholds", Prefs.get("stegers.preset_user_thresholds", USE_PRESET_USER_THRESHOLDS));
		gd.addNumericField("Upper_threshold", Prefs.get("stegers.upper_threshold", UPPER_THRESHOLD), 0);
		gd.addNumericField("Lower_threshold", Prefs.get("stegers.lower_threshold", LOWER_THRESHOLD), 0);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Filter_multiresponse", Prefs.get("stegers.filter_multiresponse", FILTER_MULTIRESPONSE));
		gd.addNumericField("Multireponse_angle_threshold", Prefs.get("stegers.multiresponse_threshold", MULTIRESPONSE_THRESHOLD), 2);
		
		gd.addCheckbox("Prevent_splinter_fragments", Prefs.get("stegers.prevent_splinter_fragments", PREVENT_SPLINTER_FRAGMENTS));
		gd.addNumericField("Splinter_fragment_threshold", Prefs.get("stegers.splinter_fragment_threshold", SPLINTER_FRAGMENT_THRESHOLD), 0);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Include_junctions", Prefs.get("stegers.include_junctions", INCLUDE_JUNCTIONS_IN_LINE));
		gd.addCheckbox("Include_first_processed", Prefs.get("stegers.include_first_processed", INCLUDE_FIRST_PROCESSED_IN_LINE));
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addNumericField("Cost_function_weight", Prefs.get("stegers.cost_function_weight", COST_FUNCTION_WEIGHT), 2);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Enable_linking", Prefs.get("stegers.enable_linking", ENABLE_LINKING));
		gd.addNumericField("Max_bending_angle", Prefs.get("stegers.max_bending_angle", MAX_BENDING_ANGLE), 0);
		
		gd.addCheckbox("Filter_short_line_segments", Prefs.get("stegers.filter_short_line_segments", FILTER_SHORT_LINE_SEGMENTS));
		gd.addNumericField("Line_length_threshold", Prefs.get("stegers.line_length_threshold", LINE_LENGTH_THRESHOLD), 0);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Enable_debug_mode", Prefs.get("stegers.debug_mode", DEBUG_MODE_ENABLED));
		
		gd.setOKLabel("Trace");
		
		gd.showDialog();
		if(gd.wasCanceled()) return;
		
		// retrieve parameters
		SIGMA = gd.getNextNumber();
		LINEPOINT_THRESHOLD = gd.getNextNumber();
		
		USE_GAUSSIAN_FIT_PEAKS_BEFORE = gd.getNextBoolean();
		USE_GAUSSIAN_FIT_PEAKS_AFTER = gd.getNextBoolean();
		
		USE_PRESET_USER_THRESHOLDS = gd.getNextBoolean();
		UPPER_THRESHOLD = gd.getNextNumber();
		LOWER_THRESHOLD = gd.getNextNumber();
		
		FILTER_MULTIRESPONSE = gd.getNextBoolean();
		MULTIRESPONSE_THRESHOLD = gd.getNextNumber();
		
		PREVENT_SPLINTER_FRAGMENTS = gd.getNextBoolean();
		SPLINTER_FRAGMENT_THRESHOLD = (int)gd.getNextNumber();
		
		INCLUDE_JUNCTIONS_IN_LINE = gd.getNextBoolean();
		INCLUDE_FIRST_PROCESSED_IN_LINE = gd.getNextBoolean();
		
		COST_FUNCTION_WEIGHT = gd.getNextNumber();
		
		ENABLE_LINKING = gd.getNextBoolean();
		MAX_BENDING_ANGLE = gd.getNextNumber();
		
		FILTER_SHORT_LINE_SEGMENTS = gd.getNextBoolean();
		LINE_LENGTH_THRESHOLD = (int)gd.getNextNumber();
		
		DEBUG_MODE_ENABLED = gd.getNextBoolean();
		
		// store parameters in preferences
		Prefs.set("stegers.sigma", SIGMA);
		Prefs.set("stegers.linepoint_threshold", LINEPOINT_THRESHOLD);
		Prefs.set("stegers.use_gaussian_fit_peaks_before", USE_GAUSSIAN_FIT_PEAKS_BEFORE);
		Prefs.set("stegers.use_gaussian_fit_peaks_after", USE_GAUSSIAN_FIT_PEAKS_AFTER);
		Prefs.set("stegers.preset_user_thresholds", USE_PRESET_USER_THRESHOLDS);
		Prefs.set("stegers.upper_threshold", UPPER_THRESHOLD);
		Prefs.set("stegers.lower_threshold", LOWER_THRESHOLD);
		Prefs.set("stegers.filter_multiresponse", FILTER_MULTIRESPONSE);
		Prefs.set("stegers.multiresponse_threshold", MULTIRESPONSE_THRESHOLD);
		Prefs.set("stegers.prevent_splinter_fragments", PREVENT_SPLINTER_FRAGMENTS);
		Prefs.set("stegers.splinter_fragment_threshold", SPLINTER_FRAGMENT_THRESHOLD);
		Prefs.set("stegers.include_junctions", INCLUDE_JUNCTIONS_IN_LINE);
		Prefs.set("stegers.include_first_processed", INCLUDE_FIRST_PROCESSED_IN_LINE);
		Prefs.set("stegers.cost_function_weight", COST_FUNCTION_WEIGHT);
		Prefs.set("stegers.enable_linking", ENABLE_LINKING);
		Prefs.set("stegers.max_bending_angle", MAX_BENDING_ANGLE);
		Prefs.set("stegers.filter_short_line_segments", FILTER_SHORT_LINE_SEGMENTS);
		Prefs.set("stegers.line_length_threshold", LINE_LENGTH_THRESHOLD);
		Prefs.set("stegers.debug_mode", DEBUG_MODE_ENABLED);
		
		// execute filter
		exec(imp, SIGMA, LINEPOINT_THRESHOLD);
		/*ImagePlus result = exec(imp, sigma, threshold);
		
		// show image
		if(null != result)
		{
			result.show();
		}*/
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public void exec(ImagePlus imp, double sigma, double linepoint_threshold)
	{
		// TODO: check arguments
		stegersAlgorithm(imp, sigma, linepoint_threshold);
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public static void stegersAlgorithm(ImagePlus imp)
	{
		stegersAlgorithm(imp, SIGMA);
	}
	
	public static void stegersAlgorithm(ImagePlus imp, double sigma)
	{
		stegersAlgorithm(imp, sigma, LINEPOINT_THRESHOLD);
	}
	
	public static void stegersAlgorithm(ImagePlus imp, double sigma, double linepoint_threshold)
	{
		run(imp.getProcessor(), sigma, linepoint_threshold); // Stegers.run()
	}
	
	public static ImageProcessor run(ImageProcessor ip, double sigma, double linepoint_threshold)
	{
		// image properties
		int image_width = ip.getWidth();
		int image_height = ip.getHeight();
		
		// *********************************************************************
		
		// STEP: calculate Gaussian derivatives, Hessian matrix, eigenvalues and eigenvector, determine Taylor polynomical response function peak position
		// 
		Profiling.tic();
		IJ.showStatus("Calculating Eigen decomposition of Hessian matrix");				
		double[][][] line_points;
		line_points = DerivativeOfGaussian.get_line_points(ip, SIGMA);
		Profiling.toc("Calculating Eigen decomposition of Hessian matrix");
		
		// =====================================================================
		
		Profiling.tic();
		ImagePlus normals_imp = IJ.createImage("normals", image_width, image_height, 4, 32);
		int[] indices = new int[]{1, 2, 4, 5};
		for(int i = 0; i < indices.length; ++i)
		{
			ImageProcessor slice_ip = normals_imp.getImageStack().getProcessor(i+1);
			for(int px = 0; px < image_width; ++px)
			{
				for(int py = 0; py < image_height; ++py)
				{
					slice_ip.setf(px, py, (float)line_points[px][py][indices[i]]);
				}
			}
		}
		normals_imp.show();
		Profiling.toc("Generating results table");
		
		return null;
	}
}
