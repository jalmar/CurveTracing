
package plugins;

// import Java classes
import java.util.Collections;
import java.util.Comparator;
import java.util.Vector;
import java.util.Map;
import java.util.HashMap;
import java.util.HashSet;
import java.awt.Color;

//import java.awt.event.MouseListener; 
import java.awt.event.MouseAdapter; // more convenient class than MouseListener interface
import java.awt.event.MouseEvent;

// import ImageJ classes
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.LUT;

import ij.plugin.PlugIn;
import ij.gui.GenericDialog;
import ij.gui.DialogListener;
import ij.Prefs;

import ij.gui.Overlay;
import ij.gui.Roi;
//import ij.gui.Line; // LineRoi // NOTE: will conflict with core.Line
import ij.gui.OvalRoi;
import ij.gui.PolygonRoi; // for POLYLINE

import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import ij.gui.Plot;
import ij.gui.PlotWindow;

import ij.measure.ResultsTable;

import ij.plugin.frame.RoiManager;

// import Jama classes
import Jama.Matrix;
import Jama.EigenvalueDecomposition;

// import own classes
import algorithms.LevenbergMarquardt;
//import algorithms.Stegers;
import filters.DerivativeOfGaussian;
//import filters.Laplacian;
import filters.Median;
import filters.TopHatTransform;
import filters.FillHoles;
import filters.ConnectedComponents;
import utils.ImageArithmetic;
import utils.Profiling;
import utils.Tuple;
import core.Point;
import core.Line;

/**
 *	Microtubule tracing algorithm
 */
public class Trace_Microtubules_Modified_Stegers implements PlugIn
{
	/**
	 *	Parameters
	 */
	
	// sigma
	public static double SIGMA = 2.0;
	
	// noise reduction and background suppression
	public static boolean REDUCE_NOISE = true;
	public static int MEDIAN_FILTER_SIZE = 3;
	public static boolean SUPPRESS_BACKGROUND = true;
	
	// Frangi vesselness measure
	//public static double FRANGI_ALPHA = 0.5;
	public static double FRANGI_BETA = 0.5;
	public static double FRANGI_C = 0.5*65535; // NOTE: assumes 16-bit images
	
	// profile fitting
	public static int SAMPLE_RATE = 1;
	public static double SAMPLE_RATE_INV = 1.0 / (double)SAMPLE_RATE;
	public static int INTERPOLATION_METHOD = ImageProcessor.BILINEAR;
	public static final int LMA_NUM_ITERATIONS = 5;
	public static final double LMA_DEFAULT_LAMBDA = 0.001;
	
	// filter criteria
	public static int AMPLITUDE_LOWER_THRESHOLD = 0;
	public static int AMPLITUDE_UPPER_THRESHOLD = 2*65535; // NOTE: assumes 16-bit images
	public static double MU_THRESHOLD = 0.707;
	public static double SIGMA_RATIO_LOWER_THRESHOLD = 0.67;
	public static double SIGMA_RATIO_UPPER_THRESHOLD = 1.50;
	public static double R_SQUARED_THRESHOLD = 0.50;
	
	public static boolean FILTER_ON_AMPLITUDE = true;
	public static boolean FILTER_ON_MU = true;
	public static boolean FILTER_ON_SIGMA = true;
	public static boolean FILTER_ON_R_SQUARED = true;
	
	// additional filtering
	public static boolean FILTER_HIT_MAP = true;
	public static int VALID_LINE_POINT_THRESHOLD = 2;
//	public static boolean SKELETONIZE_VALID_LINE_POINTS_MAP = true;
//	public static boolean FILTER_JUNCTION_POINTS = true;
//	public static boolean FILTER_SINGLE_PIXEL_ISLETS = true;
	
	// tracing
	public static boolean USE_GAUSSIAN_FIT_PEAKS = true;
	public static double WEIGHT_HIT_COUNT = 2.0;
	public static double WEIGHT_R_SQUARED = 3.0;
	public static double WEIGHT_DISTANCE = 1.0;
	public static double WEIGHT_ORIENTATION = 2.0;
	
	// linking
	public static boolean ENABLE_LINKING = true;
	public static double SEARCH_DISTANCE = (int)(2*Math.ceil(6*SIGMA));
	public static int BACKTRACK_DISTANCE = 3;
	public static int AVERAGING_DISTANCE = 5;
	
//	// ODF map
//	public static double ANGULAR_RESOLUTION = 6.0; // degrees per orientation step; NOTE: should be a divisor of 180.0
//	public static int ODF_STEPS = (int)(180.0 / ANGULAR_RESOLUTION);
//	public static double ANISOTROPY_FACTOR = 5;
//	public static int ODF_PEAK_WIDTH = 5;
	
	// additional filtering
	public static int LINE_LENGTH_THRESHOLD = 5;
	
	// debug mode
	public static boolean DEBUG_MODE_ENABLED = true;
	public static boolean USE_WEIGHTED_AVERAGES = true;
	public static boolean SHOW_RESULTS_TABLE = false;
	public static boolean SHOW_VECTOR_OVERLAY = true;
//	public static boolean SHOW_ODF_OVERLAY = true;
	public static boolean USE_VECTOR_NORMALS = true;
	public static boolean USE_TRANSPARENCY = true;
	
	public static int SLICE_INDEX = -1;
	
	/**
	 *
	 */
	public Trace_Microtubules_Modified_Stegers()
	{
		/* nothing */
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public void run(String arg)
	{
		// get the current image
		ImagePlus imp = IJ.getImage();
		if(imp == null) return;
		
		// show dialog with options
		GenericDialog gd = new GenericDialog("Trace microtubules modified Stegers");
		gd.addNumericField("PSF_sigma", Prefs.get("mt_trace.psf_sigma", SIGMA), 2);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Reduce_noise", Prefs.get("mt_trace.reduce_noise", REDUCE_NOISE));
		gd.addCheckbox("Suppress_background", Prefs.get("mt_trace.suppress_background", SUPPRESS_BACKGROUND));
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Filter_on_amplitude", Prefs.get("mt_trace.filter_amplitude", true));
		gd.addNumericField("Amplitude_lower_threshold", Prefs.get("mt_trace.filter_amplitude_lower_threshold", AMPLITUDE_LOWER_THRESHOLD), 0);
		gd.addNumericField("Amplitude_upper_threshold", Prefs.get("mt_trace.filter_amplitude_upper_threshold", AMPLITUDE_UPPER_THRESHOLD), 0);
		
		gd.addCheckbox("Filter_on_mu", Prefs.get("mt_trace.filter_mu", true));
		gd.addNumericField("Mu_absolute_threshold", Prefs.get("mt_trace.filter_mu_threshold", MU_THRESHOLD), 3);
		
		gd.addCheckbox("Filter_on_sigma", Prefs.get("mt_trace.filter_sigma", true));
		gd.addNumericField("Sigma_ratio_lower_threshold", Prefs.get("mt_trace.filter_sigma_lower_threshold", SIGMA_RATIO_LOWER_THRESHOLD), 2);
		gd.addNumericField("Sigma_ratio_upper_threshold", Prefs.get("mt_trace.filter_sigma_upper_threshold", SIGMA_RATIO_UPPER_THRESHOLD), 2);
		
		gd.addCheckbox("Filter_on_R-squared", Prefs.get("mt_trace.filter_r_squared", true));
		gd.addNumericField("R-squared_threshold", Prefs.get("mt_trace.filter_r_squared_threshold", R_SQUARED_THRESHOLD), 2);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Filter_hit_map", Prefs.get("mt_trace.filter_hit_map", FILTER_HIT_MAP));
		gd.addNumericField("Valid_line_point_threshold", Prefs.get("mt_trace.valid_line_point_threshold", VALID_LINE_POINT_THRESHOLD), 0);
//		gd.addCheckbox("Skeletonize_valid_line_points_map", Prefs.get("mt_trace.skeletonize_valid_line_points_map", SKELETONIZE_VALID_LINE_POINTS_MAP));
//		gd.addCheckbox("Filter_junction_points", Prefs.get("mt_trace.filter_junction_points", FILTER_JUNCTION_POINTS));
//		gd.addCheckbox("Filter_single_pixel_islets", Prefs.get("mt_trace.filter_single_pixel_islets", FILTER_SINGLE_PIXEL_ISLETS));
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Use_Gaussian_fit_peaks", Prefs.get("mt_trace.use_gaussian_fit_peaks", USE_GAUSSIAN_FIT_PEAKS));
		gd.addNumericField("Hit_count_weight", Prefs.get("mt_trace.weight_hit_count", WEIGHT_HIT_COUNT), 3);
		gd.addNumericField("R_squared_weight", Prefs.get("mt_trace.weight_r_squared", WEIGHT_R_SQUARED), 3);
		gd.addNumericField("Distance_weight", Prefs.get("mt_trace.weight_distance", WEIGHT_DISTANCE), 3);
		gd.addNumericField("Orientation_weight", Prefs.get("mt_trace.weight_orientation", WEIGHT_ORIENTATION), 3);
		
//		gd.setInsets(10, 20, 0); // seperate parameter groups
//		
//		gd.addCheckbox("Enable_linking", Prefs.get("mt_trace.enable_linking", //ENABLE_LINKING));
//		gd.addNumericField("Search_distance", Prefs.get("mt_trace.search_distance", SEARCH_DISTANCE), 1);
//		gd.addNumericField("Backtrack_distance", Prefs.get("mt_trace.backtrack_distance", BACKTRACK_DISTANCE), 0);
//		gd.addNumericField("Averaging_distance", Prefs.get("mt_trace.averaging_distance", AVERAGING_DISTANCE), 0);
		
//		gd.setInsets(10, 20, 0); // seperate parameter groups
//		
//		gd.addNumericField("Angular_resolution", Prefs.get("mt_trace.angular_resolution", ANGULAR_RESOLUTION), 1);
//		gd.addNumericField("Anisotropy_factor", Prefs.get("mt_trace.anisotropy_factor", ANISOTROPY_FACTOR), 1);
//		gd.addNumericField("ODF_peak_width", Prefs.get("mt_trace.odf_peak_width", ODF_PEAK_WIDTH), 0);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addNumericField("Line_length_lower_threshold", Prefs.get("mt_trace.line_length_lower_threshold", LINE_LENGTH_THRESHOLD), 0);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Enable_debug_mode", Prefs.get("mt_trace.debug_mode", DEBUG_MODE_ENABLED));
		gd.addCheckbox("Use_weighted_averages", Prefs.get("mt_trace.use_weighted_averages", USE_WEIGHTED_AVERAGES));
		gd.addCheckbox("Show_results_table", Prefs.get("mt_trace.show_results_table", SHOW_RESULTS_TABLE));
		gd.addCheckbox("Show_vector_overlay", Prefs.get("mt_trace.show_vector_overlay", SHOW_VECTOR_OVERLAY));
//		gd.addCheckbox("Show_ODF_overlay", Prefs.get("mt_trace.show_odf_overlay", SHOW_ODF_OVERLAY));
		gd.addCheckbox("Use_vector_normals", Prefs.get("mt_trace.use_vector_normals", USE_VECTOR_NORMALS));
		gd.addCheckbox("Use_transparency", Prefs.get("mt_trace.use_transparency", USE_TRANSPARENCY));
		
		gd.setOKLabel("Trace");
		
		gd.showDialog();
		if(gd.wasCanceled()) return;
		
		// retrieve parameters
		SIGMA = gd.getNextNumber();
		REDUCE_NOISE = gd.getNextBoolean();
		SUPPRESS_BACKGROUND = gd.getNextBoolean();
		
		FILTER_ON_AMPLITUDE = gd.getNextBoolean();
		AMPLITUDE_LOWER_THRESHOLD = (int)gd.getNextNumber();
		AMPLITUDE_UPPER_THRESHOLD = (int)gd.getNextNumber();
		
		FILTER_ON_MU = gd.getNextBoolean();
		MU_THRESHOLD = gd.getNextNumber();
		
		FILTER_ON_SIGMA = gd.getNextBoolean();
		SIGMA_RATIO_LOWER_THRESHOLD = gd.getNextNumber();
		SIGMA_RATIO_UPPER_THRESHOLD = gd.getNextNumber();
		
		FILTER_ON_R_SQUARED = gd.getNextBoolean();
		R_SQUARED_THRESHOLD = gd.getNextNumber();
		
		FILTER_HIT_MAP = gd.getNextBoolean();
		VALID_LINE_POINT_THRESHOLD = (int)gd.getNextNumber();
//		SKELETONIZE_VALID_LINE_POINTS_MAP = gd.getNextBoolean();
//		FILTER_JUNCTION_POINTS = gd.getNextBoolean();
//		FILTER_SINGLE_PIXEL_ISLETS = gd.getNextBoolean();
		
		USE_GAUSSIAN_FIT_PEAKS = gd.getNextBoolean();
		WEIGHT_HIT_COUNT = gd.getNextNumber();
		WEIGHT_R_SQUARED = gd.getNextNumber();
		WEIGHT_DISTANCE = gd.getNextNumber();
		WEIGHT_ORIENTATION = gd.getNextNumber();
		
//		ENABLE_LINKING = gd.getNextBoolean();
//		SEARCH_DISTANCE = gd.getNextNumber();
//		BACKTRACK_DISTANCE = (int)gd.getNextNumber();
//		AVERAGING_DISTANCE = (int)gd.getNextNumber();
		
//		ANGULAR_RESOLUTION = gd.getNextNumber();
//		ODF_STEPS = (int)(180.0 / ANGULAR_RESOLUTION); // NOTE: also update ODF_STEPS
//		ANISOTROPY_FACTOR = gd.getNextNumber();
//		ODF_PEAK_WIDTH = (int)gd.getNextNumber();
		
		LINE_LENGTH_THRESHOLD = (int)gd.getNextNumber();
		
		DEBUG_MODE_ENABLED = gd.getNextBoolean();
		USE_WEIGHTED_AVERAGES = gd.getNextBoolean();
		SHOW_RESULTS_TABLE = gd.getNextBoolean();
		SHOW_VECTOR_OVERLAY = gd.getNextBoolean();
//		SHOW_ODF_OVERLAY = gd.getNextBoolean();
		USE_VECTOR_NORMALS = gd.getNextBoolean();
		USE_TRANSPARENCY = gd.getNextBoolean();
		
		// store parameters in preferences
		Prefs.set("mt_trace.psf_sigma", SIGMA);
		Prefs.set("mt_trace.reduce_noise", REDUCE_NOISE);
		Prefs.set("mt_trace.suppress_background", SUPPRESS_BACKGROUND);
		
		Prefs.set("mt_trace.filter_amplitude", FILTER_ON_AMPLITUDE);
		Prefs.set("mt_trace.filter_amplitude_lower_threshold", AMPLITUDE_LOWER_THRESHOLD);
		Prefs.set("mt_trace.filter_amplitude_upper_threshold", AMPLITUDE_UPPER_THRESHOLD);
		
		Prefs.set("mt_trace.filter_mu", FILTER_ON_MU);
		Prefs.set("mt_trace.filter_mu_threshold", MU_THRESHOLD);
		
		Prefs.set("mt_trace.filter_sigma", FILTER_ON_SIGMA);
		Prefs.set("mt_trace.filter_sigma_lower_threshold", SIGMA_RATIO_LOWER_THRESHOLD);
		Prefs.set("mt_trace.filter_sigma_upper_threshold", SIGMA_RATIO_UPPER_THRESHOLD);
		
		Prefs.set("mt_trace.filter_r_squared", FILTER_ON_R_SQUARED);
		Prefs.set("mt_trace.filter_r_squared_threshold", R_SQUARED_THRESHOLD);
		
		Prefs.set("mt_trace.filter_hit_map", FILTER_HIT_MAP);
		Prefs.set("mt_trace.valid_line_point_threshold", VALID_LINE_POINT_THRESHOLD);
//		Prefs.set("mt_trace.skeletonize_valid_line_points_map", SKELETONIZE_VALID_LINE_POINTS_MAP);
//		Prefs.set("mt_trace.filter_junction_points", FILTER_JUNCTION_POINTS);
//		Prefs.set("mt_trace.filter_single_pixel_islets", FILTER_SINGLE_PIXEL_ISLETS);
		
		Prefs.set("mt_trace.use_gaussian_fit_peaks", USE_GAUSSIAN_FIT_PEAKS);
		Prefs.set("mt_trace.weight_hit_count", WEIGHT_HIT_COUNT);
		Prefs.set("mt_trace.weight_r_squared", WEIGHT_R_SQUARED);
		Prefs.set("mt_trace.weight_distance", WEIGHT_DISTANCE);
		Prefs.set("mt_trace.weight_orientation", WEIGHT_ORIENTATION);
		
//		Prefs.set("mt_trace.enable_linking", ENABLE_LINKING);
//		Prefs.set("mt_trace.search_distance", SEARCH_DISTANCE);
//		Prefs.get("mt_trace.backtrack_distance", BACKTRACK_DISTANCE);
//		Prefs.get("mt_trace.averaging_distance", AVERAGING_DISTANCE);
		
//		Prefs.set("mt_trace.angular_resolution", ANGULAR_RESOLUTION);
//		Prefs.set("mt_trace.anisotropy_factor", ANISOTROPY_FACTOR);
//		Prefs.set("mt_trace.odf_peak_width", ODF_PEAK_WIDTH);
		
		Prefs.set("mt_trace.line_length_lower_threshold", LINE_LENGTH_THRESHOLD);
		
		Prefs.set("mt_trace.debug_mode", DEBUG_MODE_ENABLED);
		Prefs.set("mt_trace.use_weighted_averages", USE_WEIGHTED_AVERAGES);
		Prefs.set("mt_trace.show_results_table", SHOW_RESULTS_TABLE);
		Prefs.set("mt_trace.show_vector_overlay", SHOW_VECTOR_OVERLAY);
//		Prefs.set("mt_trace.show_odf_overlay", SHOW_ODF_OVERLAY);
		Prefs.set("mt_trace.use_vector_normals", USE_VECTOR_NORMALS);
		Prefs.set("mt_trace.use_transparency", USE_TRANSPARENCY);
		
		// trace microtubules
		run(imp, SIGMA, REDUCE_NOISE, SUPPRESS_BACKGROUND);
		/*ImagePlus result = run(imp, sigma, reduce_noise, suppress_background);
		
		// show result
		if (result != null)
		{
			result.show();
		}*/
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	/**
	 *
	 */
	public static ImagePlus run(ImagePlus imp, double sigma, boolean reduce_noise, boolean suppress_background)
	{
		// create new image for output
		ImagePlus imp_out = IJ.createImage("Microtubule traces of " + imp.getTitle(), imp.getWidth(), imp.getHeight(), imp.getNSlices(), 8);
		
		// process slices in image stack
		ImageStack stack_in = imp.getStack();
		ImageStack stack_out = imp_out.getStack();
		for(int slice = 1; slice <= stack_in.getSize(); ++slice)
		{
			SLICE_INDEX = slice;
			ImageProcessor slice_ip = stack_in.getProcessor(slice);
			ImageProcessor slice_result = run(slice_ip, sigma, reduce_noise, suppress_background);
			if(slice_result != null)
			{
				stack_out.setProcessor(slice_result, slice);
			}
			else
			{
				// RSLV: add empty slice?
			}
		}
		imp_out.setStack(stack_out);
		imp_out.resetDisplayRange();
		imp_out.show();
		
		// return image
		return imp_out;
	}
	
	public static ImageProcessor run(ImageProcessor ip, double psf_sigma, boolean reduce_noise, boolean suppress_background)
	{
		// duplicate image processor
		ImageProcessor ip_original = ip;
		ImageProcessor ip_dup = ip_original.duplicate(); // NOTE: duplicate!
		int image_width = ip.getWidth();
		int image_height = ip.getHeight();
		
		// get image statistics
		int original_image_max_intensity = 0;
		int original_image_min_intensity = 65535; // NOTE: assumes 16-bit images!
		double original_image_avg_intensity = 0.0;
		for(int px = 0; px < image_width; ++px)
		{
			for(int py = 0; py < image_height; ++py)
			{
				int pv = ip_original.get(px, py);
				if(pv > original_image_max_intensity) original_image_max_intensity = pv;
				if(pv < original_image_min_intensity) original_image_min_intensity = pv;
				original_image_avg_intensity += pv;
			}
		}
		original_image_avg_intensity /= image_width * image_height;
		
		// ---------------------------------------------------------------------
		
		// Step 1: median filtering to reduce camera shot noise
		ImageProcessor ip_step_1 = ip_dup;
		if(reduce_noise)
		{
			Profiling.tic();
			
			IJ.showStatus("Removing shot noise");
			//ip_dup.medianFilter(); // NOTE: only works on 8-bit or RGB images
			ip_step_1 = Median.run(ip_dup, MEDIAN_FILTER_SIZE);
			
			if(DEBUG_MODE_ENABLED)
			{
				//ip_step_1.resetMinAndMax();
				ImagePlus debug_imp = new ImagePlus("DEBUG: median image", ip_step_1);
				debug_imp.resetDisplayRange();
				debug_imp.show();
			}
			
			Profiling.toc("Step 1: Removing shot noise");
		}
		else
		{
			System.err.println("Skipping step 1: removing shot noise");
		}
		
		// ---------------------------------------------------------------------
		
		// Step 2: white top hat transform to suppress background noise
		ImageProcessor ip_step_2 = ip_step_1.duplicate(); // NOTE: duplicate!
		if(suppress_background)
		{
			Profiling.tic();
			IJ.showStatus("Removing background noise");
			ImageProcessor ip_wth = TopHatTransform.run(ip_step_1, SIGMA); // RSLV: slightly larger than expected sigma?
			ip_step_2 = ImageArithmetic.subtract(ip_step_1, ip_wth);
			
			if(DEBUG_MODE_ENABLED)
			{
				//ip_wth.resetMinAndMax();
				ImagePlus debug_imp = new ImagePlus("DEBUG: top hat transform image", ip_wth);
				debug_imp.resetDisplayRange();
				debug_imp.show();
				
				//ip_step_2.resetMinAndMax();
				debug_imp = new ImagePlus("DEBUG: background removed image", ip_step_2);
				debug_imp.resetDisplayRange();
				debug_imp.show();
			}
			
			Profiling.toc("Step 2: Removing background");
		}
		else
		{
			System.err.println("Skipping step 2: removing background");
		}
		
		// get image statistics
		// RLSV: use ImageStatistics.getStatistics(ImageProcessor ip, int mOptions, Calibration cal)
		int preprocessed_image_max_intensity = 0;
		int preprocessed_image_min_intensity = 65535; // NOTE: assumes 16-bit images!
		double preprocessed_image_avg_intensity = 0.0;
		for(int px = 0; px < image_width; ++px)
		{
			for(int py = 0; py < image_height; ++py)
			{
				int pv = ip_step_2.get(px, py);
				if(pv > preprocessed_image_max_intensity) preprocessed_image_max_intensity = pv;
				if(pv < preprocessed_image_min_intensity) preprocessed_image_min_intensity = pv;
				preprocessed_image_avg_intensity += pv;
			}
		}
		preprocessed_image_avg_intensity /= image_width * image_height;
		
		// RSLV: updating Frangi C?
		FRANGI_C = preprocessed_image_avg_intensity;
		System.err.println("DEBUG: automatically updating FRANGI_C to average intensity of preprocessed image, FRANGI_C=" + FRANGI_C);
		
		// ---------------------------------------------------------------------
		
		// Step 3a: eigenvalue/vector decomposition of Hessian matrix
		
		double[][][] results_step_3;
		
		// store Frangi measures on eigenvalues; NOTE: |L1| <= |L2|
		double[][][] frangi_measures;
				
		Profiling.tic();
		IJ.showStatus("Step 3b: Calculating Eigen decomposition of Hessian matrix");

		//get line points
		results_step_3 = DerivativeOfGaussian.get_line_points(ip_step_2, SIGMA);
		
		frangi_measures = DerivativeOfGaussian.frangi_measures( results_step_3, image_height, image_width, FRANGI_BETA, FRANGI_C); 
		
		Profiling.toc("Step 3b: Calculating Eigen decomposition of Hessian matrix");

		
		// show intermediate images
		if(DEBUG_MODE_ENABLED)
		{
			// eigenvalues and eigenvectors for debug images
			ImageProcessor first_eigenvalues_ip = new FloatProcessor(image_width, image_height);
			ImageProcessor second_eigenvalues_ip = new FloatProcessor(image_width, image_height);
			
			ImageProcessor first_eigenvectors_ip = new ByteProcessor(image_width*5, image_height*5);
			ImageProcessor second_eigenvectors_ip = new ByteProcessor(image_width*5, image_height*5);
		
			ImageProcessor first_eigenvectors_theta_ip = new FloatProcessor(image_width, image_height);
			ImageProcessor second_eigenvectors_theta_ip = new FloatProcessor(image_width, image_height);
			
			ImageProcessor first_theta_direction_ip = new FloatProcessor(image_width, image_height);
			ImageProcessor second_theta_direction_ip = new FloatProcessor(image_width, image_height);
			
			ImageProcessor frangi_blobness_ip = new FloatProcessor(image_width, image_height);
			ImageProcessor frangi_structureness_ip = new FloatProcessor(image_width, image_height);
			ImageProcessor frangi_vesselness_ip = new FloatProcessor(image_width, image_height);
			
			// fill images with data
			for(int py = 0; py < image_height; ++py)
			{
				for( int px = 0; px < image_width; ++px)
				{
					// DEBUG: store eigenvalues in image processor
					first_eigenvalues_ip.setf(px, py, (float)results_step_3[px][py][0]);
					second_eigenvalues_ip.setf(px, py, (float)results_step_3[px][py][3]);
					
					// DEBUG: store eigenvectors in image processor
					int cx = px*5+2;
					int cy = py*5+2;
					first_eigenvectors_ip.set((int)cx, (int)cy, 255);
					first_eigenvectors_ip.set((int)Math.floor(cx-results_step_3[px][py][1]), (int)Math.floor(cy-results_step_3[px][py][2]), 255);
					first_eigenvectors_ip.set((int)Math.ceil(cx+results_step_3[px][py][1]), (int)Math.ceil(cy+results_step_3[px][py][2]), 255);
					
					second_eigenvectors_ip.set((int)cx, (int)cy, 255);
					second_eigenvectors_ip.set((int)Math.floor(cx-results_step_3[px][py][4]), (int)Math.floor(cy-results_step_3[px][py][5]), 255);
					second_eigenvectors_ip.set((int)Math.ceil(cx+results_step_3[px][py][4]), (int)Math.ceil(cy+results_step_3[px][py][5]), 255);
					
					// store orientation of eigenvectors
					double first_theta = Math.atan2(results_step_3[px][py][2], results_step_3[px][py][1]);
					double second_theta = Math.atan2(results_step_3[px][py][5], results_step_3[px][py][4]);
					
					first_eigenvectors_theta_ip.setf(px, py, (float)first_theta);
					second_eigenvectors_theta_ip.setf(px, py, (float)second_theta);
					
					// store direction (bi-directional range)?
					// NOTE: map to CCW coordinate system with 0 degree pointing right
					if(first_theta < 0) first_theta = Math.PI + first_theta; // map negative to positive
					if(second_theta < 0) second_theta = Math.PI + second_theta; // map negative to positive
					
					first_theta = first_theta % Math.PI; // just to be sure
					second_theta = second_theta % Math.PI; // just to be sure
					
					first_theta = first_theta * (180.0 / Math.PI); // radians to degrees
					second_theta = second_theta * (180.0 / Math.PI); // radians to degrees
					
					first_theta_direction_ip.setf(px, py, (float)first_theta);
					second_theta_direction_ip.setf(px, py, (float)second_theta);
					
					// store Frangi measures
					frangi_blobness_ip.setf(px, py, (float)frangi_measures[px][py][2]);
					frangi_structureness_ip.setf(px, py, (float)frangi_measures[px][py][3]);
					frangi_vesselness_ip.setf(px, py, (float)frangi_measures[px][py][4]);
				}
			}
			
			// first and second eigenvalues
			ImagePlus first_eigenvalues_imp = new ImagePlus("DEBUG: first eigenvalues", first_eigenvalues_ip);
			first_eigenvalues_imp.resetDisplayRange();
			first_eigenvalues_imp.show();
			
			ImagePlus second_eigenvalues_imp = new ImagePlus("DEBUG: second eigenvalues", second_eigenvalues_ip);
			second_eigenvalues_imp.resetDisplayRange();
			second_eigenvalues_imp.show();
			
			// first and second eigenvectors
			ImagePlus first_eigenvectors_imp = new ImagePlus("DEBUG: first eigenvectors", first_eigenvectors_ip);
			first_eigenvectors_imp.resetDisplayRange();
//			first_eigenvectors_imp.show();
			
			ImagePlus second_eigenvectors_imp = new ImagePlus("DEBUG: second eigenvectors", second_eigenvectors_ip);
			second_eigenvectors_imp.resetDisplayRange();
//			second_eigenvectors_imp.show();
			
			// theta orientation
			ImagePlus first_eigenvectors_theta_imp = new ImagePlus("DEBUG: theta of first eigenvectors", first_eigenvectors_theta_ip);
			first_eigenvectors_theta_imp.resetDisplayRange();
//			first_eigenvectors_theta_imp.show();
			
			ImagePlus second_eigenvectors_theta_imp = new ImagePlus("DEBUG: theta of second eigenvectors", second_eigenvectors_theta_ip);
			second_eigenvectors_theta_imp.resetDisplayRange();
//			second_eigenvectors_theta_imp.show();
			
			// direction
			ImagePlus first_theta_direction_imp = new ImagePlus("DEBUG: direction of theta of first eigenvectors", first_theta_direction_ip);
			first_theta_direction_imp.resetDisplayRange();
//			first_theta_direction_imp.show();
			
			ImagePlus second_theta_direction_imp = new ImagePlus("DEBUG: direction of theta of second eigenvectors", second_theta_direction_ip);
			second_theta_direction_imp.resetDisplayRange();
//			second_theta_direction_imp.show();
			
			// Frangi measures
			ImagePlus frangi_blobness_imp = new ImagePlus("DEBUG: Frangi blobness measure", frangi_blobness_ip);
			frangi_blobness_imp.resetDisplayRange();
//			frangi_blobness_imp.show();
			
			ImagePlus frangi_structureness_imp = new ImagePlus("DEBUG: Frangi structureness measure", frangi_structureness_ip);
			frangi_structureness_imp.resetDisplayRange();
//			frangi_structureness_imp.show();
			
			ImagePlus frangi_vesselness_imp = new ImagePlus("DEBUG: Frangi vesselness measure", frangi_vesselness_ip);
			frangi_vesselness_imp.resetDisplayRange();
			frangi_vesselness_imp.show();
		}
		
		// ---------------------------------------------------------------------
		
		// Step 4: perform line profile fitting to data
		Profiling.tic();
		IJ.showStatus("Fitting Gaussian line profiles to pixels");
		
		//ImageProcessor fitting_ip = ip_original;
		double[][][] fitting_results = new double[image_width][image_height][4]; // [x][y][bg=0|amp=1|mu=2|sigma=3]
		double[][][] standard_error_fit_results = new double[image_width][image_height][4]; // [x][y][bg=0|amp=1|mu=2|sigma=3]
		double[][] chi_squared_fit_results = new double[image_width][image_height];
		double[][] r_squared_fit_results = new double[image_width][image_height];
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				// get center pixel and vector orientation from (previous) step 3
				double cx = px;
				double cy = py;
				double nx = results_step_3[px][py][1];
				double ny = results_step_3[px][py][2];
				
				// extract line profile data from *original* image
				ip.setInterpolationMethod(INTERPOLATION_METHOD);
				int line_profile_width = (int)Math.ceil(3*SIGMA); // RSLV: 3*sigma minimum? Currently using ceil!
				int data_points = 2*SAMPLE_RATE*line_profile_width+1;
				double[] x_data = new double[data_points];
				double[] y_data = new double[data_points];
				double min_value = Double.MAX_VALUE; // keep track of minimum value
				double max_value = Double.MIN_VALUE; // keep track of maximum value
				double mean_value = 0.0;
				double ii = -line_profile_width;
				for(int i = 0; i < data_points; ++i)
				{
					// interpolated x,y
					double ix = cx + ii * nx;
					double iy = cy + ii * ny;
					double pv = ip.getPixelInterpolated(ix, iy); // NOTE: *original* image!
					x_data[i] = ii; // NOTE use relative x-coordinate!!
					y_data[i] = pv;
					
					// update min/max value
					min_value = Math.min(min_value, pv);
					max_value = Math.max(max_value, pv);
					mean_value += pv;
					
					// increment ii
					ii += SAMPLE_RATE_INV;
				}
				mean_value /= data_points;
				
				// initial parameter estimation
				double background = min_value;
				double amplitude = max_value - min_value;
				double mu = 0; // 0 = center in relative coordinate system
				double sig = SIGMA;
				double[] initial_parameters = new double[]{background, amplitude, mu, sig};
				
				// set up new LMA instance
				LevenbergMarquardt lma = new LevenbergMarquardt();
				
				// run LMA fitting procedure
				double[] fitted_parameters = lma.run(x_data, y_data, data_points, initial_parameters, LMA_NUM_ITERATIONS, LMA_DEFAULT_LAMBDA);
				double[] standard_errors = lma.calculateStandardErrors(x_data, y_data, data_points, fitted_parameters);
				double chi_squared = lma.calculateChi2(x_data, y_data, data_points, fitted_parameters);
				
				// calculate R^2 measure
				double sstot = 0.0;
				for(int i = 0; i < data_points; ++i)
				{
					double dd = y_data[i] - mean_value; // difference with mean
					sstot += dd*dd; // sum of sqaure
				}
				double r_squared = 1.0 - (chi_squared / sstot);
				
				// store result of fitting
				fitting_results[px][py] = fitted_parameters;
				standard_error_fit_results[px][py] = standard_errors;
				chi_squared_fit_results[px][py] = chi_squared;
				r_squared_fit_results[px][py] = r_squared;
				
				// override Taylor peak estimation
				if(USE_GAUSSIAN_FIT_PEAKS)
				{
					results_step_3[px][py][6] = fitted_parameters[2]*nx;
					results_step_3[px][py][7] = fitted_parameters[2]*ny;
				}
			}
		}
		Profiling.toc("Step 4: Fitting line profiles");
		
		// DEBUG: show intermediate images
		if(DEBUG_MODE_ENABLED)
		{
			// generate images
			ImageProcessor background_fit_ip = new FloatProcessor(image_width, image_height);
			ImageProcessor amplitude_fit_ip = new FloatProcessor(image_width, image_height);
			ImageProcessor mu_fit_ip = new FloatProcessor(image_width, image_height);
			ImageProcessor mu_squared_fit_ip = new FloatProcessor(image_width, image_height);
			ImageProcessor sigma_fit_ip = new FloatProcessor(image_width, image_height);
			ImageProcessor chi_squared_fit_ip = new FloatProcessor(image_width, image_height);
			ImageProcessor log_chi_squared_fit_ip = new FloatProcessor(image_width, image_height);
			ImageProcessor r_squared_fit_ip = new FloatProcessor(image_width, image_height);
			for(int py = 0; py < image_height; ++py)
			{
				for( int px = 0; px < image_width; ++px)
				{
					background_fit_ip.setf(px, py, (float)fitting_results[px][py][0]);
					amplitude_fit_ip.setf(px, py, (float)fitting_results[px][py][1]);
					mu_fit_ip.setf(px, py, (float)fitting_results[px][py][2]);
					mu_squared_fit_ip.setf(px, py, (float)(fitting_results[px][py][2]*fitting_results[px][py][2]));
					sigma_fit_ip.setf(px, py, (float)Math.abs(fitting_results[px][py][3]));
					chi_squared_fit_ip.setf(px, py, (float)chi_squared_fit_results[px][py]);
					log_chi_squared_fit_ip.setf(px, py, (float)Math.log(chi_squared_fit_results[px][py]));
					r_squared_fit_ip.setf(px, py, (float)r_squared_fit_results[px][py]);
				}
			}
			
			// smallest/largest [absolute] eigenvalues
//			ImagePlus background_fit_imp = new ImagePlus("DEBUG: background fit", background_fit_ip);
//			background_fit_imp.resetDisplayRange();
//			background_fit_imp.show();
			
//			ImagePlus amplitude_fit_imp = new ImagePlus("DEBUG: amplitude fit", amplitude_fit_ip);
//			amplitude_fit_imp.resetDisplayRange();
//			amplitude_fit_imp.show();
			
//			ImagePlus mu_fit_imp = new ImagePlus("DEBUG: mu fit", mu_fit_ip);
//			mu_fit_imp.resetDisplayRange();
//			mu_fit_imp.setDisplayRange(-SIGMA, SIGMA);
//			mu_fit_imp.show();
			
			ImagePlus mu_squared_fit_imp = new ImagePlus("DEBUG: mu squared fit", mu_squared_fit_ip);
			mu_squared_fit_imp.resetDisplayRange();
			mu_squared_fit_imp.setDisplayRange(0, SIGMA*SIGMA);
			mu_squared_fit_imp.show();
			
			ImagePlus sigma_fit_imp = new ImagePlus("DEBUG: sigma fit", sigma_fit_ip);
			//sigma_fit_imp.resetDisplayRange();
			sigma_fit_imp.setDisplayRange(0, 5*SIGMA); // RSLV: 5 * max scale space sigma?
			sigma_fit_imp.show();
			
//			ImagePlus chi_fit_imp = new ImagePlus("DEBUG: chi sqaured fit", chi_fit_ip);
//			chi_fit_imp.resetDisplayRange();
//			chi_fit_imp.show();
			
//			ImagePlus log_chi_squared_fit_imp = new ImagePlus("DEBUG: log chi squared fit", log_chi_squared_fit_ip);
//			log_chi_squared_fit_imp.resetDisplayRange();
//			log_chi_squared_fit_imp.show();
			
			ImagePlus r_squared_fit_imp = new ImagePlus("DEBUG: r squared fit", r_squared_fit_ip);
			r_squared_fit_imp.setDisplayRange(0, 1);
			r_squared_fit_imp.show();
		}
		
		// *********************************************************************
		
		// Step 5a: additional filtering of point (determine which pixels are microtubules and which belong to the background) based on selection criteria
		Profiling.tic();
		IJ.showStatus("Filter line points and creating projection maps");
		
		// create hit map
		ImageProcessor hit_map_ip = new ByteProcessor(image_width, image_height);
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				// apply filter criteria
				boolean filtered = false;
				
				// filter on amplitude
				if(FILTER_ON_AMPLITUDE && (fitting_results[px][py][1] < AMPLITUDE_LOWER_THRESHOLD || fitting_results[px][py][1] > AMPLITUDE_UPPER_THRESHOLD))
				{
					// probably not a valid center line pixel; skip
					filtered = true; //continue;
				}
				
				// filter on mu
				if(FILTER_ON_MU && Math.abs(fitting_results[px][py][2]) > MU_THRESHOLD)
				{
					// probably not a valid center line pixel; skip
					filtered = true; //continue;
				}
				
				// filter on sigma
				if(FILTER_ON_SIGMA && (fitting_results[px][py][3] < SIGMA_RATIO_LOWER_THRESHOLD * SIGMA || fitting_results[px][py][3] > SIGMA_RATIO_UPPER_THRESHOLD * SIGMA))
				{
					// probably not a valid center line pixel; skip
					filtered = true; //continue;
				}
				
				// filter on r-squared
				if(FILTER_ON_R_SQUARED && (Double.isNaN(r_squared_fit_results[px][py]) || r_squared_fit_results[px][py] < R_SQUARED_THRESHOLD))
				{
					// probably not a valid center line pixel; skip
					filtered = true; //continue;
				}
				
				// set pixel in hit images
				if(!filtered)
				{
					// set hit image pixel
					hit_map_ip.set(px, py, 255);
				}
			}
		}
		Profiling.toc("Step 5a: Filter line points and creating projection maps");
		
		// ---------------------------------------------------------------------
		
		// Step 5b: [optional] filter out smaller particles using opening and closing
		// RSLV: use connected components with minimum size filter instead? Benefits of not changing morphology of remaining components!
		
		ImageProcessor hit_map_filtered_ip = hit_map_ip.duplicate();
		if(FILTER_HIT_MAP)
		{
			Profiling.tic();
			((ByteProcessor)hit_map_filtered_ip).erode(3, 0);
			((ByteProcessor)hit_map_filtered_ip).dilate(3, 0);
			((ByteProcessor)hit_map_filtered_ip).dilate(3, 0);
			((ByteProcessor)hit_map_filtered_ip).erode(3, 0);
			Profiling.toc("Step 5b: Filter hit map");
		}
		
		// RSLV: filter using connected components instead
		// NOTE: because of recursion it can run out of stack memory with large segments!
//		ImageProcessor hit_map_filtered_ip = hit_map_ip.duplicate();
//		ConnectedComponents.run(hit_map_filtered_ip, ConnectedComponents.Connectivity.EIGHT_CONNECTIVITY, 10, 1000000);
		
		// ---------------------------------------------------------------------
		
		// Step 5c: create projection maps from hit map
		
		// TODO: average measures!!!!
		
		// store averaged results of eigendecomposition
		//	[px][py][0] = lambda1_magnitude		n(t)
		//	[px][py][1] = lambda1_direction_x	n_x(t)
		//	[px][py][2] = lambda1_direction_y	n_y(t)
		//	[px][py][3] = lambda2_magnitude		s(t)
		//	[px][py][4] = lambda2_direction_x	s_x(t)
		//	[px][py][5] = lambda2_direction_y	s_y(t)
		//	[px][py][6] = super-resolved_x		t_x, or dlpx
		//	[px][py][7] = super-resolved_y		t_y, or dlpy
		double[][][] avg_results_step_3 = new double[image_width][image_height][8];
		// RSLV: how to average n_x, n_y and s_x, s_y? Using average theta
		
		// store Frangi measures on eigenvalues; NOTE: |L1| <= |L2|
		// beta is control parameter, set at 0.5
		// c dependens on image bit depth, about half maximum Hessian matrix norm
		//	[0] = frangi L1
		//	[1] = frangi L2
		//	[2] = blobness (eccentricity), L1 / L2; note: keep sign!
		//	[3] = second order structureness, RSS of all elements, or Frobius norm
		//	[4] = vesselness = exp(-[2]^2/2*FRANGI_BETA^2)(1-exp(-[3]^2/2*FRANGI_C^2))
//		double[][][] avg_frangi_measures = new double[image_width][image_height][5];
		// RSLV: condition |L1| <= |L2| may change by averging?
		
		double[][][] avg_fitting_results = new double[image_width][image_height][4]; // [x][y][bg=0|amp=1|mu=2|sigma=3]
//		double[][][] avg_standard_error_fit_results = new double[image_width][image_height][4]; // [x][y][bg=0|amp=1|mu=2|sigma=3]
//		double[][] avg_chi_squared_fit_results = new double[image_width][image_height];
		double[][] avg_r_squared_fit_results = new double[image_width][image_height];
		
		
		// project all valid line points in hit map
		Profiling.tic();
		ImageProcessor hit_count_map_ip = new ByteProcessor(image_width, image_height); // count of hits
		ImageProcessor weight_map_ip = new FloatProcessor(image_width, image_height); // sum of weights
		
		int hit_count_max = 0; // keep track of maximum hit count
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				if(hit_map_filtered_ip.get(px, py) != 0)
				{
					// correct for peak position in coordinate
					int cpx = (int)(px+0.5+fitting_results[px][py][2]*results_step_3[px][py][1]);
					int cpy = (int)(py+0.5+fitting_results[px][py][2]*results_step_3[px][py][2]);
					
					// bounds checking
					if(cpx >= 0 && cpx < image_width && cpy >= 0 && cpy < image_height)
					{
						// update hit count
						int cc = 1 + hit_count_map_ip.get(cpx, cpy);
						hit_count_map_ip.set(cpx, cpy, cc);
						
						// retain maximum hit count value
						if(cc > hit_count_max)
						{
							hit_count_max = cc;
						}
						
						// calculate weight
						double projection_distance = fitting_results[px][py][2];
						double weight = 1.0;
						if(USE_WEIGHTED_AVERAGES)
						{
							weight = Math.abs(gaussianPDF(SIGMA, projection_distance));
						}
						weight_map_ip.setf(cpx, cpy, (float)(weight_map_ip.getf(cpx, cpy) + weight));
						
						// update sums for averages
						avg_results_step_3[cpx][cpy][0] += weight * results_step_3[px][py][0]; // L1
						avg_results_step_3[cpx][cpy][1] += weight * avg_results_step_3[px][py][1]; // V1_x; RSLV: average n_x;
						avg_results_step_3[cpx][cpy][2] += weight * avg_results_step_3[px][py][2]; // V1_y; RSLV: average n_y
						avg_results_step_3[cpx][cpy][3] += weight * results_step_3[px][py][3]; // L2
						avg_results_step_3[cpx][cpy][4] += weight * avg_results_step_3[px][py][4]; // V2_x; RSLV: average s_x
						avg_results_step_3[cpx][cpy][5] += weight * avg_results_step_3[px][py][5]; // V2_y; RSLV: average s_y
						avg_results_step_3[cpx][cpy][6] += weight * results_step_3[px][py][6]; // dlpx; RSLV: no weight?
						avg_results_step_3[cpx][cpy][7] += weight * results_step_3[px][py][7]; // dlpy; RSLV: no weight?
/*						
						avg_frangi_measures[cpx][cpy][0] = weight * frangi_measures[px][py][0];
						avg_frangi_measures[cpx][cpy][1] = weight * frangi_measures[px][py][1];
						avg_frangi_measures[cpx][cpy][2] = weight * frangi_measures[px][py][2];
						avg_frangi_measures[cpx][cpy][3] = weight * frangi_measures[px][py][3];
						avg_frangi_measures[px][py][4] = weight * frangi_measures[px][py][4];
*/						
						avg_fitting_results[cpx][cpy][0] += weight * fitting_results[px][py][0]; // BG
						avg_fitting_results[cpx][cpy][1] += weight * fitting_results[px][py][1]; // AMP
						avg_fitting_results[cpx][cpy][2] += weight * fitting_results[px][py][2]; // MU; RSLV: no weight?
						avg_fitting_results[cpx][cpy][3] += weight * fitting_results[px][py][3]*fitting_results[px][py][3]; // SIGMA; NOTE: average of sigma is sqrt sum of variance (sigma squared)!
/*						
						avg_standard_error_fit_results[cpx][cpy][0] += weight * standard_error_fit_results[px][py][0];
						avg_standard_error_fit_results[cpx][cpy][1] += weight * standard_error_fit_results[px][py][1];
						avg_standard_error_fit_results[cpx][cpy][2] += weight * standard_error_fit_results[px][py][2];
						avg_standard_error_fit_results[cpx][cpy][3] += weight * standard_error_fit_results[px][py][3];
						
						avg_chi_squared_fit_results[cpx][cpy] += weight * chi_squared_fit_results[px][py];
*/						
						avg_r_squared_fit_results[cpx][cpy] += weight * r_squared_fit_results[px][py]; // RSLV: average R sqaured is just linear sum?
					}
				}
			}
		}
		
		// create average projection maps
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				double weight = weight_map_ip.getf(px, py);
				if(weight > 0) // NOTE: avoid division by zero
				{
					double inv_weight = 1.0 / weight;
					
					avg_results_step_3[px][py][0] *= inv_weight; // L1
					avg_results_step_3[px][py][1] *= inv_weight; // V1_x; RSLV: average n_x
					avg_results_step_3[px][py][2] *= inv_weight; // V1_y; RSLV: average n_y
					avg_results_step_3[px][py][3] *= inv_weight; // L2
					avg_results_step_3[px][py][4] *= inv_weight; // V2_x; RSLV: average s_x
					avg_results_step_3[px][py][5] *= inv_weight; // V2_y; RSLV: average s_y
					avg_results_step_3[px][py][6] *= inv_weight; // dlpx
					avg_results_step_3[px][py][7] *= inv_weight; // dlpy
/*					
					avg_frangi_measures[px][py][0] /= hc;
					avg_frangi_measures[px][py][1] /= hc;
					avg_frangi_measures[px][py][2] /= hc;
					avg_frangi_measures[px][py][3] /= hc;
					avg_frangi_measures[px][py][4] /= hc;
*/					
					avg_fitting_results[px][py][0] *= inv_weight; // BG
					avg_fitting_results[px][py][1] *= inv_weight; // AMP
					avg_fitting_results[px][py][2] *= inv_weight; // MU
					avg_fitting_results[px][py][3] *= inv_weight; // SIGMA; NOTE: average of sigma is sqrt sum of variance (sigma squared)!
					//avg_fitting_results[px][py][3] = Math.sqrt(avg_fitting_results[px][py][3]);
/*					
					avg_standard_error_fit_results[px][py][0] /= hc;
					avg_standard_error_fit_results[px][py][1] /= hc;
					avg_standard_error_fit_results[px][py][2] /= hc;
					avg_standard_error_fit_results[px][py][3] /= hc;
					
					avg_chi_squared_fit_results[px][py] /= hc;
*/					
					avg_r_squared_fit_results[px][py] *= inv_weight; // RSLV: average R sqaured is just linear sum?
				}
			}
		}
		Profiling.toc("Step 5c: Create (averaged) projection maps");
		
		// DEBUG: show intermediate images
		if(DEBUG_MODE_ENABLED)
		{
			// generate images
			// TODO: create images for other projected/averaged parameters
			ImageProcessor r_squared_avg_projection_map_ip = new FloatProcessor(image_width, image_height);
			for(int py = 0; py < image_height; ++py)
			{
				for(int px = 0; px < image_width; ++px)
				{
					r_squared_avg_projection_map_ip.setf(px, py, (float)avg_r_squared_fit_results[px][py]);
				}
			}
			
			// -----------------------------------------------------------------
			
			// show hit images
			ImagePlus hit_map_imp = new ImagePlus("DEBUG: hit map", hit_map_ip);
			hit_map_imp.resetDisplayRange();
			hit_map_imp.show();
			
			ImagePlus hit_map_filtered_imp = new ImagePlus("DEBUG: hit map filtered", hit_map_filtered_ip);
			hit_map_filtered_imp.resetDisplayRange();
			hit_map_filtered_imp.show();
			
			ImagePlus hit_count_map_imp = new ImagePlus("DEBUG: hit count map", hit_count_map_ip);
			hit_count_map_imp.setDisplayRange(0, hit_count_max); //resetDisplayRange();
			hit_count_map_imp.show();
			
			// show R^2 projection map
			ImagePlus r_squared_avg_projection_map_imp = new ImagePlus("DEBUG: r-squared avg projection map", r_squared_avg_projection_map_ip);
			r_squared_avg_projection_map_imp.setDisplayRange(0, 1); //resetDisplayRange();
			r_squared_avg_projection_map_imp.show();
		}
		
		// ---------------------------------------------------------------------
		
		// Step 5d: filter hit count map
		
		ImageProcessor valid_line_points_map_ip = new ByteProcessor(image_width, image_height);
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				// filter pixels with less than VALID_LINE_POINT_THRESHOLD hits
				if(hit_count_map_ip.get(px, py) >= VALID_LINE_POINT_THRESHOLD)
				{
					valid_line_points_map_ip.set(px, py, 255);
				}
				//else
				//{
				//	valid_line_points_mask_ip.set(px, py, 0);
				//}
			}
		}
		
/*
		// skeletonize valid line point map
		ImageProcessor valid_line_points_skeleton_ip = valid_line_points_map_ip.duplicate();
		if(SKELETONIZE_VALID_LINE_POINTS_MAP)
		{
			valid_line_points_skeleton_ip.invert();
			((ByteProcessor)valid_line_points_skeleton_ip).skeletonize();
			valid_line_points_skeleton_ip.invert();
		}
		
		// RSLV: after skeletonization, should be possible to pick up on junctions; they are points with more than two neighours?
		ImageProcessor filter_junctions_ip = valid_line_points_skeleton_ip.duplicate();
		if(FILTER_JUNCTION_POINTS)
		{
			for(int py = 0; py < image_height; ++py)
			{
				for(int px = 0; px < image_width; ++px)
				{
					int sum = 0;
					for(int ky = -1; ky <= 1; ++ky)
					{
						for(int kx = -1; kx <= 1; ++kx)
						{
							sum += valid_line_points_skeleton_ip.getPixel(px+kx, py+ky); // NOTE: getPixel for bounds checking!
							// RSLV: could use filter_jucntions_ip to avoid filtering too many junction points; i.e. after removing one junction, second point may no longer be regarded as a junction?
						}
					}
					
					if(sum >= 4*255) // NOTE: assumes binary mask with [BG=0|FG=255]
					{
						filter_junctions_ip.set(px, py, 0);
					}
				}
			}
		}
		
		// remove single pixel islets
		// RSLV: could attempt removing everything with less than three neighbours? Would filter out two-pixel islets AND line end points. However, it filters out quite a few true points in SNR <= 3
		ImageProcessor valid_line_points_filtered_ip = filter_junctions_ip.duplicate();
		if(FILTER_SINGLE_PIXEL_ISLETS)
		{
			for(int py = 0; py < image_height; ++py)
			{
				for(int px = 0; px < image_width; ++px)
				{
					if(valid_line_points_filtered_ip.get(px, py) == 255)
					{
						int sum = 0;
						for(int ky = -1; ky <= 1; ++ky)
						{
							for(int kx = -1; kx <= 1; ++kx)
							{
								sum += filter_junctions_ip.getPixel(px+kx, py+ky); // NOTE: getPixel for bounds checking!
							}
						}
						
						// filter single islets
						if(sum <= 255) // NOTE: assumes binary mask with [BG=0|FG=255]
						{
							valid_line_points_filtered_ip.set(px, py, 0);
						}
					}
				}
			}
		}
*/
		
		// DEBUG: show intermediate images
		if(DEBUG_MODE_ENABLED)
		{
			// show valid line points map image
			ImagePlus valid_line_points_map_imp = new ImagePlus("DEBUG: valid line points map", valid_line_points_map_ip);
			valid_line_points_map_imp.resetDisplayRange();
			valid_line_points_map_imp.show();
			
/*
			if(SKELETONIZE_VALID_LINE_POINTS_MAP)
			{
				ImagePlus valid_line_points_skeleton_imp = new ImagePlus("DEBUG: valid line points skeleton", valid_line_points_skeleton_ip);
				valid_line_points_skeleton_imp.resetDisplayRange();
				valid_line_points_skeleton_imp.show();
			}
			
			if(FILTER_JUNCTION_POINTS)
			{
				ImagePlus filter_junctions_imp = new ImagePlus("DEBUG: junction points filtered", filter_junctions_ip);
				filter_junctions_imp.resetDisplayRange();
				filter_junctions_imp.show();
			}
			
			if(FILTER_SINGLE_PIXEL_ISLETS)
			{
				ImagePlus valid_line_points_filtered_imp = new ImagePlus("DEBUG: valid line points filtered", valid_line_points_filtered_ip);
				valid_line_points_filtered_imp.resetDisplayRange();
				valid_line_points_filtered_imp.show();
			}
*/
		}
		
		// *********************************************************************
		
		// STEP 3a: first sort all seeding points in decending order
		Vector<Vector<Double> > seeding_points = new Vector<Vector<Double> >();
		for(int px = 0; px < image_width; ++px)
		{
			for(int py = 0; py < image_height; ++py)
			{
				// use only valid line points
				if(valid_line_points_map_ip.get(px, py) == 255)
				{
					// create point vector
					Vector<Double> p = new Vector<Double>(3);
					p.add(new Double(px));
					p.add(new Double(py));
					p.add(new Double((double)(hit_count_map_ip.get(px, py)) / (double)(hit_count_max))); // hit count
					//p.add(new Double(0.0)); // RSLV: avg mu?
					p.add(new Double(avg_r_squared_fit_results[px][py])); // avg R^2
					seeding_points.add(p);
				}
			}
		}
		System.err.println("Number of seeding points = " + seeding_points.size());
		Collections.sort(seeding_points, new Comparator<Vector<Double> >() {
			public int compare(Vector<Double> o1, Vector<Double> o2)
			{
				Double cost_o1 = new Double(WEIGHT_HIT_COUNT * o1.get(2) + WEIGHT_R_SQUARED * o1.get(3));
				Double cost_o2 = new Double(WEIGHT_HIT_COUNT * o2.get(2) + WEIGHT_R_SQUARED * o2.get(3));
				return cost_o2.compareTo(cost_o1); // descending order
			}
		});
		
		// trace lines, starting from highest seeding point
		Line[][] point_to_line_map = new Line[image_width][image_height];
		ImageProcessor processing_map_ip = new ByteProcessor(image_width, image_height);
		int AVAILABLE = 0x00;
		int LIMITED = 0x04;
		int LINEPOINT = 0x02;
		int JUNCTION = 0x08;
		Vector<Line> lines = new Vector<Line>();
		for(Vector<Double> seeding_point : seeding_points)
		{
			// get current pixel information
			int cp_px = seeding_point.get(0).intValue();
			int cp_py = seeding_point.get(1).intValue();
			double cp_dlpx = results_step_3[cp_px][cp_py][6];
			double cp_dlpy = results_step_3[cp_px][cp_py][7];
			
			// get orientation of point
			double cp_nx = results_step_3[cp_px][cp_py][1];
			double cp_ny = results_step_3[cp_px][cp_py][2];
			double cp_na = getAngle(cp_nx, cp_ny);
			int cp_ni = getOrientationIndex(cp_na);
			double cp_sx = results_step_3[cp_px][cp_py][4];
			double cp_sy = results_step_3[cp_px][cp_py][5];
			double cp_sa = getAngle(cp_sx, cp_sy);
			int cp_si = getOrientationIndex(cp_sa);
			
			// check if not already processed
			if(processing_map_ip.get(cp_px, cp_py) == AVAILABLE)
			{
				// start local processing map for line
				ImageProcessor local_processing_map_ip = new ByteProcessor(image_width, image_height);
				
				// add point to line
				Line line = new Line();
				line.add(new Point(cp_px, cp_py, cp_dlpx, cp_dlpy));
				point_to_line_map[cp_px][cp_py] = line;
				processing_map_ip.set(cp_px, cp_py, LINEPOINT);
				local_processing_map_ip.set(cp_px, cp_py, 255);
				lines.add(line); // @#$%@#%@#%
				
				// trace in both directions
				double[] trace_directions = new double[]{0, 180};
				for(int tdi = 0; tdi < trace_directions.length; ++tdi)
				{
					double trace_direction = trace_directions[tdi];
					System.err.println("Trace in direction " + trace_direction);
					
					// trace line in first direction
					int tp_px = cp_px;
					int tp_py = cp_py;
					
					double tp_sao = (cp_sa+trace_direction)%360.0;
					
					boolean tracing = true;
					while(tracing)
					{
						// get information on current trace point
						double tp_dlpx = results_step_3[tp_px][tp_py][6];
						double tp_dlpy = results_step_3[tp_px][tp_py][7];
						
						// get orientation of point
						double tp_nx = results_step_3[tp_px][tp_py][1];
						double tp_ny = results_step_3[tp_px][tp_py][2];
						double tp_na = getAngle(tp_nx, tp_ny);
						int tp_ni = getOrientationIndex(tp_na);
						double tp_sx = results_step_3[tp_px][tp_py][4];
						double tp_sy = results_step_3[tp_px][tp_py][5];
						double tp_sa = getAngle(tp_sx, tp_sy);
						int tp_si = getOrientationIndex(tp_sa);
						
						// reorient vector information in global tracing direction
						System.err.println("angle="+tp_sa);
						System.err.println("index="+tp_si);
						if(Math.abs(tp_sa - tp_sao) > 90.0)
						{
							System.err.println("Reorienting s(t) vector");
							tp_sx = -tp_sx;
							tp_sy = -tp_sy;
							tp_sa = getAngle(tp_sx, tp_sy);
							tp_si = getOrientationIndex(tp_sa);
							System.err.println("corrected angle="+tp_sa);
							System.err.println("corrected index="+tp_si);
						}
						
						// determine next best neighbour
						int np_px = -1;
						int np_py = -1;
						double last_cost = Double.MAX_VALUE;
						
						System.err.println("tp_px="+tp_px + "    tp_py="+tp_py);
						
						// check all three neighbours
						for(int i = -1; i <= 1; ++i)
						{
							// calculate new orientation index within range [0..7]
							int tp_sic = tp_si + i;
							if(tp_sic < 0) tp_sic = 8 + tp_sic;
							if(tp_sic > 7) tp_sic = tp_sic % 8;
							
							// determine next neighbour position
							int[] dp = getNextPixel(tp_sic);
							//int[] dp = getNextPixel((tp_ni + i) % 8);
							System.err.println("  dp[0]=" + dp[0] + "    dp[1]=" + dp[1]);
							int dp_px = tp_px + dp[0];
							int dp_py = tp_py + dp[1];
							System.err.println("  dp_px=" + dp_px + "    dp_py=" + dp_py);
							
							// skip if outside of image window
							if(dp_px < 0 || dp_px >= image_width || dp_py < 0 || dp_py >= image_height)
							{
								System.err.println("  Next point lies outside image window");
								continue;
							}
							
							// skip if not in hit count map
							if(hit_count_map_ip.get(dp_px, dp_py) == 0)
							{
								System.err.println("  Next point does not have a registered hit");
								continue;
							}
							
							// get next neighbour information
							double dp_dlpx = results_step_3[dp_px][dp_py][6];
							double dp_dlpy = results_step_3[dp_px][dp_py][7];
							
							// get orientation of point
							double dp_nx = results_step_3[dp_px][dp_py][1];
							double dp_ny = results_step_3[dp_px][dp_py][2];
							double dp_na = getAngle(dp_nx, dp_ny); //Math.atan2(dp_ny, dp_nx);
							int dp_ni = getOrientationIndex(dp_na);
							double dp_sx = results_step_3[dp_px][dp_py][4];
							double dp_sy = results_step_3[dp_px][dp_py][5];
							double dp_sa = getAngle(dp_sx, dp_sy); //Math.atan2(dp_sy, dp_sx);
							int dp_si = getOrientationIndex(dp_sa);
							
							// calclate cost function
							double ddist_x = (dp_px + 0.5 + dp_dlpx) - (tp_px + 0.5 + tp_dlpx);
							double ddist_y = (dp_py + 0.5 + dp_dlpy) - (tp_py + 0.5 + tp_dlpy);
							double ddist = Math.sqrt(ddist_x * ddist_x + ddist_y * ddist_y);
							double dtheta = Math.abs(dp_sa - tp_sa); // NOTE: use s(t), not n(t), although should not matter
							System.err.println("dtheta="+dtheta);
							if(dtheta > 180.0) // NOTE: bigger than, not equal!
							{
								dtheta = 360.0 - dtheta;
								System.err.println("corrected dtheta="+dtheta);
							}
							if(dtheta > 90.0) // NOTE: bigger than, not equal!
							{
								dtheta = 180.0 - dtheta;
								System.err.println("corrected dtheta="+dtheta);
							}
							dtheta = Math.toRadians(dtheta);
							
							double hit_count = (double)(hit_count_map_ip.get(dp_px, dp_py)) / (double)(hit_count_max);
							//double peak_offset = 0.0; // RLSV: use avg mu?
							double r_squared = avg_r_squared_fit_results[dp_px][dp_py];
							double cost = WEIGHT_DISTANCE * ddist + WEIGHT_ORIENTATION * dtheta + WEIGHT_R_SQUARED * (1.0 - r_squared) + WEIGHT_HIT_COUNT * hit_count;
							
							System.err.println("  cost=w1*distance+w2*theta+w3*r_squared+w4*hit_count" + WEIGHT_DISTANCE + "*" + ddist + "+" + WEIGHT_ORIENTATION + "*" + dtheta + "+" + WEIGHT_R_SQUARED + "*" + r_squared + "+" + WEIGHT_HIT_COUNT + "*" + hit_count + "=" + cost);
							System.err.println("  last_cost="+last_cost);
							
							// check if current neighour is better than previous
							if(cost < last_cost)
							{
								System.err.println("    Point is better than previous point");
								// update best neighbour
								last_cost = cost;
								np_px = dp_px;
								np_py = dp_py;
							}
						}
						
						// *****************************************************
						
						// skip if no best match found (outside of image window or no valid line point)
						if(np_px == -1 || np_py == -1)
						{
							System.err.println("Did not find any suitable neighbour");
							tracing = false;
							continue; //break;
						}
						
						System.err.println("Best neighbour at ("+np_px+","+np_py+")");
						
						// get next neighbour information
						double np_dlpx = results_step_3[np_px][np_py][6];
						double np_dlpy = results_step_3[np_px][np_py][7];
						
						// get orientation of point
						double np_nx = results_step_3[np_px][np_py][1];
						double np_ny = results_step_3[np_px][np_py][2];
						double np_na = getAngle(np_nx, np_ny); //Math.atan2(cp_ny, cp_nx);
						int np_ni = getOrientationIndex(np_na);
						double np_sx = results_step_3[np_px][np_py][4];
						double np_sy = results_step_3[np_px][np_py][5];
						double np_sa = getAngle(np_sx, np_sy); //Math.atan2(cp_sy, cp_sx);
						int np_si = getOrientationIndex(np_sa);
						
						// reorient vector information in global tracing direction
						if(Math.abs(np_sa - tp_sao) > 90.0 && Math.abs(np_sa - tp_sao) < 270.0)
						{
							System.err.println("Reorienting s(t) vector");
							np_sx = -np_sx;
							np_sy = -np_sy;
							np_sa = getAngle(np_sx, np_sy);
							np_si = getOrientationIndex(np_sa);
						}
						
						// determine action
						System.err.println("Determine best action");
						if(processing_map_ip.get(np_px, np_py) == AVAILABLE || processing_map_ip.get(np_px, np_py) == LIMITED) // i.e. not a linepoint nor a junction
						{
							System.err.println("  Point is available");
							System.err.println("    status = " + processing_map_ip.get(np_px, np_py));
							
							// add point to line
							if(trace_direction == 0)
							{
								line.addFirst(new Point(np_px, np_py, np_dlpx, np_dlpy));
							}
							else
							{
								line.addLast(new Point(np_px, np_py, np_dlpx, np_dlpy));
							}
							point_to_line_map[np_px][np_py] = line;
							processing_map_ip.set(np_px, np_py, LINEPOINT);
							local_processing_map_ip.set(np_px, np_py, 255);
						}
						else if(processing_map_ip.get(np_px, np_py) == LINEPOINT)
						{
							System.err.println("  Point is part of line, mark as junction, and stop tracing");
							
							// find current point on line
							Point jp = new Point(np_px, np_py, np_dlpx, np_dlpy);
							System.err.println("jp = " + jp);
							Line line_a = point_to_line_map[np_px][np_py];
							System.err.println("line_a = " + line_a);
							System.err.println("line_a.size() = " + line_a.size());
							int jp_pos = line_a.indexOf(jp);
							System.err.println("jp_pos = " + jp_pos);
							
							// check if not the same line
							if(line_a != line)
							{
								// split the line
								System.err.println("splicing line a");
								System.err.println("line_a = " + line_a);
								System.err.println("line_a.size() = " + line_a.size());
								System.err.println("point jp = " + jp);
								System.err.println("position = " + jp_pos);
								
								Line line_b = line_a.splice(jp_pos);
								System.err.println("line_a = " + line_a);
								System.err.println("line_a.size() = " + line_a.size());
								
								System.err.println("line_b = " + line_b);
								System.err.println("line_b.size() = " + line_b.size());
								
								// RSLV: add junction point to line?
								line_a.addLast(new Point(np_px, np_py, np_dlpx, np_dlpy));
								
								lines.add(line_b);
								
								// link points to second line
								for(Point pb : line_b)
								{
									point_to_line_map[pb.px][pb.py] = line_b;
								}
								
								// mark point as junction
								processing_map_ip.set(np_px, np_py, JUNCTION);
								local_processing_map_ip.set(np_px, np_py, 255);
								point_to_line_map[np_px][np_py] = null;
							}
							
							// stop tracing
							tracing = false;
						}
						else if(processing_map_ip.get(np_px, np_py) == JUNCTION)
						{
							System.err.println("  Point is junction, stop tracing");
							
							// RLSV: include junction point in new line?
							if(trace_direction == 0)
							{
								line.addFirst(new Point(np_px, np_py, np_dlpx, np_dlpy));
							}
							else
							{
								line.addLast(new Point(np_px, np_py, np_dlpx, np_dlpy));
							}
							//point_to_line_map[np_px][np_py] = line;
							//processing_map_ip.set(np_px, np_py, LINEPOINT);
							//local_processing_map_ip.set(np_px, np_py, 255);
							
							// stop tracing
							tracing = false;
						}
						
						// continue tracing from neighbour
						tp_px = np_px;
						tp_py = np_py;
						tp_sao = np_sa;
					}
					
					System.err.println("Done tracing in direction " + trace_direction);
				}
				
				// mark surrounding of line as processed
				local_processing_map_ip.erode();
				local_processing_map_ip.erode();
				
				for(int lpy = 0; lpy < image_height; ++lpy)
				{
					for(int lpx = 0; lpx < image_width; ++lpx)
					{
						if(local_processing_map_ip.get(lpx, lpy) == 255 && processing_map_ip.get(lpx, lpy) == 0)
						{
							processing_map_ip.set(lpx, lpy, LIMITED);
						}
					}
				}
			}
		}
		
		// ---------------------------------------------------------------------
		
		if(DEBUG_MODE_ENABLED)
		{
			ImagePlus processing_map_imp = new ImagePlus("Processing map", processing_map_ip);
			processing_map_imp.setDisplayRange(0,8); // NOTE: keep up to date with number of processing types
			processing_map_imp.show();
		}
		
		// *********************************************************************
		
		// connect lines segments into larger lines
//		if(ENABLE_LINKING)
//		{
//			if(DEBUG_MODE_ENABLED)
//			{
				// find all initial line end points
				double MAX_BENDING_ANGLE = 2*22.5; // bending angle is actually FOV
				double SEARCH_DISTANCE = 2*Math.ceil(6*SIGMA); // pixels
				int BACKTRACK_DISTANCE = 3; // number of point to backtrack
				int AVERAGING_DISTANCE = 10;
				ImageProcessor search_window_overlay_ip = ip.duplicate();
				Overlay search_window_overlay = new Overlay(); // TMP: DEBUG
				
				// for each line
				for(int i = 0; i < lines.size(); ++i)
				{
					// get line
					Line l = lines.get(i);
					
					// skip if less than some size
					if(l.size() < 2+AVERAGING_DISTANCE+2*BACKTRACK_DISTANCE) continue;
					
					// for both end points
					int wx1 = l.get(0+BACKTRACK_DISTANCE).px;
					int wy1 = l.get(0+BACKTRACK_DISTANCE).py;
					double wsx1 = l.get(0+BACKTRACK_DISTANCE).sx;
					double wsy1 = l.get(0+BACKTRACK_DISTANCE).sy;
					int wdx1 = l.get(0+BACKTRACK_DISTANCE).px - l.get(1+BACKTRACK_DISTANCE).px;
					int wdy1 = l.get(0+BACKTRACK_DISTANCE).py - l.get(1+BACKTRACK_DISTANCE).py;
					double wdsx1 = (l.get(0+BACKTRACK_DISTANCE).px + l.get(0+BACKTRACK_DISTANCE).sx) - (l.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).px + l.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).sx);
					double wdsy1 = (l.get(0+BACKTRACK_DISTANCE).py + l.get(0+BACKTRACK_DISTANCE).sy) - (l.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).py + l.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).sy);
					double wsa1 = Math.atan2(wdsy1, wdsx1);
					double lp1 = 0.0; // TODO: calculate persistence length
					
					int wx2 = l.get(l.size()-BACKTRACK_DISTANCE-1).px;
					int wy2 = l.get(l.size()-BACKTRACK_DISTANCE-1).py;
					double wsx2 = l.get(l.size()-BACKTRACK_DISTANCE-1).sx;
					double wsy2 = l.get(l.size()-BACKTRACK_DISTANCE-1).sy;
					int wdx2 = l.get(l.size()-BACKTRACK_DISTANCE-1).px - l.get(l.size()-BACKTRACK_DISTANCE-2).px;
					int wdy2 = l.get(l.size()-BACKTRACK_DISTANCE-1).py - l.get(l.size()-BACKTRACK_DISTANCE-2).py;
					double wdsx2 = (l.get(l.size()-BACKTRACK_DISTANCE-1).px + l.get(l.size()-BACKTRACK_DISTANCE-1).sx) - (l.get(l.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).px + l.get(l.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).sx);
					double wdsy2 = (l.get(l.size()-BACKTRACK_DISTANCE-1).py + l.get(l.size()-BACKTRACK_DISTANCE-1).sy) - (l.get(l.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).py + l.get(l.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).sy);
					double wsa2 = Math.atan2(wdsy2, wdsx2);
					double lp2 = 0.0; // TODO: calculate persistence length
					
					// define search window
					double wtlx1 = wx1 + Math.cos(wsa1-Math.PI/2)*0.5*1.0 + 0.5 + wsx1; // NOTE: half pex offset
					double wtly1 = wy1 + Math.sin(wsa1-Math.PI/2)*0.5*1.0 + 0.5 + wsy1; // NOTE: half pex offset
					double wtrx1 = wx1 + Math.cos(wsa1+Math.PI/2)*0.5*1.0 + 0.5 + wsx1; // NOTE: half pex offset
					double wtry1 = wy1 + Math.sin(wsa1+Math.PI/2)*0.5*1.0 + 0.5 + wsy1; // NOTE: half pex offset
					
					double wtlx2 = wx2 + Math.cos(wsa2-Math.PI/2)*0.5*1.0 + 0.5 + wsx2; // NOTE: half pex offset
					double wtly2 = wy2 + Math.sin(wsa2-Math.PI/2)*0.5*1.0 + 0.5 + wsy2; // NOTE: half pex offset
					double wtrx2 = wx2 + Math.cos(wsa2+Math.PI/2)*0.5*1.0 + 0.5 + wsx2; // NOTE: half pex offset
					double wtry2 = wy2 + Math.sin(wsa2+Math.PI/2)*0.5*1.0 + 0.5 + wsy2; // NOTE: half pex offset
					
					// DEBUG: draw search windows in overlay
					double persistence_length = Math.toRadians(MAX_BENDING_ANGLE);
					
					Roi w1 = new PolygonRoi(new float[]{(float)wtlx1, (float)(wtlx1+Math.cos(wsa1-(0.5*persistence_length))*SEARCH_DISTANCE), (float)(wtrx1+Math.cos(wsa1+(0.5*persistence_length))*SEARCH_DISTANCE), (float)wtrx1}, new float[]{(float)wtly1, (float)(wtly1+Math.sin(wsa1-(0.5*persistence_length))*SEARCH_DISTANCE), (float)(wtry1+Math.sin(wsa1+(0.5*persistence_length))*SEARCH_DISTANCE), (float)wtry1}, PolygonRoi.POLYGON);
					w1.setStrokeWidth(0.0);
					w1.setStrokeColor(Color.YELLOW);
					w1.setPosition(1);
					search_window_overlay.add(w1);
					Roi w1c = new OvalRoi(wx1+0.5+wsx1-0.125, wy1+0.5+wsy1-0.125, 0.250, 0.250);
					w1c.setStrokeWidth(0.0);
					w1c.setStrokeColor(Color.YELLOW);
					w1c.setPosition(1);
					search_window_overlay.add(w1c);
					float[] w1_line_xs = new float[BACKTRACK_DISTANCE+AVERAGING_DISTANCE];
					float[] w1_line_ys = new float[BACKTRACK_DISTANCE+AVERAGING_DISTANCE];
					for(int j = 0; j < BACKTRACK_DISTANCE+AVERAGING_DISTANCE; ++j)
					{
						w1_line_xs[j] = (float)(l.get(j).px + 0.5 + l.get(j).sx);
						w1_line_ys[j] = (float)(l.get(j).py + 0.5 + l.get(j).sy);
					}
					Roi w1l = new PolygonRoi(w1_line_xs, w1_line_ys, PolygonRoi.POLYLINE);
					w1l.setStrokeWidth(0.0);
					w1l.setStrokeColor(Color.YELLOW);
					w1l.setPosition(1);
					search_window_overlay.add(w1l);
					
					Roi w2 = new PolygonRoi(new float[]{(float)wtlx2, (float)(wtlx2+Math.cos(wsa2-(0.5*persistence_length))*SEARCH_DISTANCE), (float)(wtrx2+Math.cos(wsa2+(0.5*persistence_length))*SEARCH_DISTANCE), (float)wtrx2}, new float[]{(float)wtly2, (float)(wtly2+Math.sin(wsa2-(0.5*persistence_length))*SEARCH_DISTANCE), (float)(wtry2+Math.sin(wsa2+(0.5*persistence_length))*SEARCH_DISTANCE), (float)wtry2}, PolygonRoi.POLYGON);
					w2.setStrokeWidth(0.0);
					w2.setStrokeColor(Color.ORANGE);
					w2.setPosition(1);
					search_window_overlay.add(w2);
					Roi w2c = new OvalRoi(wx2+0.5+wsx2-0.125, wy2+0.5+wsy2-0.125, 0.250, 0.250);    
					w2c.setStrokeWidth(0.0);
					w2c.setStrokeColor(Color.ORANGE);
					w2c.setPosition(1);
					search_window_overlay.add(w2c);
					float[] w2_line_xs = new float[BACKTRACK_DISTANCE+AVERAGING_DISTANCE];
					float[] w2_line_ys = new float[BACKTRACK_DISTANCE+AVERAGING_DISTANCE];
					for(int j = 0; j < BACKTRACK_DISTANCE+AVERAGING_DISTANCE; ++j)
					{
						w2_line_xs[j] = (float)(l.get(l.size()-j-1).px + 0.5 + l.get(l.size()-j-1).sx);
						w2_line_ys[j] = (float)(l.get(l.size()-j-1).py + 0.5 +  l.get(l.size()-j-1).sy);
					}
					Roi w2l = new PolygonRoi(w2_line_xs, w2_line_ys, PolygonRoi.POLYLINE);
					w2l.setStrokeWidth(0.0);
					w2l.setStrokeColor(Color.ORANGE);
					w2l.setPosition(1);
					search_window_overlay.add(w2l);
				}
				
				ImagePlus search_window_overlay_imp = new ImagePlus("Search windows", search_window_overlay_ip);
				search_window_overlay_imp.setOverlay(search_window_overlay); // TMP: DEBUG
				//search_window_overlay_imp.updateAndRepaintWindow();
////			search_window_overlay_imp.show();
//			}
			
			// -----------------------------------------------------------------
			
			// start linking
			boolean changed = true;
			while(changed)
			{
				changed = false;
				
//				double SEARCH_DISTANCE = 2*Math.ceil(6*SIGMA); // pixels
//				int BACKTRACK_DISTANCE = 3; // number of point to backtrack
//				int AVERAGING_DISTANCE = 10;
				HashMap<Tuple<Line, Integer>, HashSet<Tuple<Line, Integer> > > candidate_matches = new HashMap<Tuple<Line, Integer>, HashSet<Tuple<Line, Integer> > >();
				
				// for each line
				for(int i = 0; i < lines.size(); ++i)
				{
					// get line
					Line l = lines.get(i);
					
					// skip if less than some size
					if(l.size() < 2+AVERAGING_DISTANCE+2*BACKTRACK_DISTANCE) continue;
					
					// for both end points
					int wx1 = l.get(0+BACKTRACK_DISTANCE).px;
					int wy1 = l.get(0+BACKTRACK_DISTANCE).py;
					double wsx1 = l.get(0+BACKTRACK_DISTANCE).sx;
					double wsy1 = l.get(0+BACKTRACK_DISTANCE).sy;
					int wdx1 = l.get(0+BACKTRACK_DISTANCE).px - l.get(1+BACKTRACK_DISTANCE).px;
					int wdy1 = l.get(0+BACKTRACK_DISTANCE).py - l.get(1+BACKTRACK_DISTANCE).py;
					double wdsx1 = (l.get(0+BACKTRACK_DISTANCE).px + l.get(0+BACKTRACK_DISTANCE).sx) - (l.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).px + l.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).sx);
					double wdsy1 = (l.get(0+BACKTRACK_DISTANCE).py + l.get(0+BACKTRACK_DISTANCE).sy) - (l.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).py + l.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).sy);
					double wsa1 = Math.atan2(wdsy1, wdsx1);
					double lp1 = 0.0; // TODO: calculate persistence length
					
					int wx2 = l.get(l.size()-BACKTRACK_DISTANCE-1).px;
					int wy2 = l.get(l.size()-BACKTRACK_DISTANCE-1).py;
					double wsx2 = l.get(l.size()-BACKTRACK_DISTANCE-1).sx;
					double wsy2 = l.get(l.size()-BACKTRACK_DISTANCE-1).sy;
					int wdx2 = l.get(l.size()-BACKTRACK_DISTANCE-1).px - l.get(l.size()-BACKTRACK_DISTANCE-2).px;
					int wdy2 = l.get(l.size()-BACKTRACK_DISTANCE-1).py - l.get(l.size()-BACKTRACK_DISTANCE-2).py;
					double wdsx2 = (l.get(l.size()-BACKTRACK_DISTANCE-1).px + l.get(l.size()-BACKTRACK_DISTANCE-1).sx) - (l.get(l.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).px + l.get(l.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).sx);
					double wdsy2 = (l.get(l.size()-BACKTRACK_DISTANCE-1).py + l.get(l.size()-BACKTRACK_DISTANCE-1).sy) - (l.get(l.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).py + l.get(l.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).sy);
					double wsa2 = Math.atan2(wdsy2, wdsx2);
					double lp2 = 0.0; // TODO: calculate persistence length
					
					// define search window
					double wtlx1 = wx1 + Math.cos(wsa1-Math.PI/2)*0.5*1.0 + 0.5 + wsx1; // NOTE: half pex offset
					double wtly1 = wy1 + Math.sin(wsa1-Math.PI/2)*0.5*1.0 + 0.5 + wsy1; // NOTE: half pex offset
					double wtrx1 = wx1 + Math.cos(wsa1+Math.PI/2)*0.5*1.0 + 0.5 + wsx1; // NOTE: half pex offset
					double wtry1 = wy1 + Math.sin(wsa1+Math.PI/2)*0.5*1.0 + 0.5 + wsy1; // NOTE: half pex offset
					
					double wtlx2 = wx2 + Math.cos(wsa2-Math.PI/2)*0.5*1.0 + 0.5 + wsx2; // NOTE: half pex offset
					double wtly2 = wy2 + Math.sin(wsa2-Math.PI/2)*0.5*1.0 + 0.5 + wsy2; // NOTE: half pex offset
					double wtrx2 = wx2 + Math.cos(wsa2+Math.PI/2)*0.5*1.0 + 0.5 + wsx2; // NOTE: half pex offset
					double wtry2 = wy2 + Math.sin(wsa2+Math.PI/2)*0.5*1.0 + 0.5 + wsy2; // NOTE: half pex offset
					
					// DEBUG: draw search windows in overlay
					double persistence_length = Math.toRadians(MAX_BENDING_ANGLE);
					
					Roi w1 = new PolygonRoi(new float[]{(float)wtlx1, (float)(wtlx1+Math.cos(wsa1-(0.5*persistence_length))*SEARCH_DISTANCE), (float)(wtrx1+Math.cos(wsa1+(0.5*persistence_length))*SEARCH_DISTANCE), (float)wtrx1}, new float[]{(float)wtly1, (float)(wtly1+Math.sin(wsa1-(0.5*persistence_length))*SEARCH_DISTANCE), (float)(wtry1+Math.sin(wsa1+(0.5*persistence_length))*SEARCH_DISTANCE), (float)wtry1}, PolygonRoi.POLYGON);
					
					Roi w2 = new PolygonRoi(new float[]{(float)wtlx2, (float)(wtlx2+Math.cos(wsa2-(0.5*persistence_length))*SEARCH_DISTANCE), (float)(wtrx2+Math.cos(wsa2+(0.5*persistence_length))*SEARCH_DISTANCE), (float)wtrx2}, new float[]{(float)wtly2, (float)(wtly2+Math.sin(wsa2-(0.5*persistence_length))*SEARCH_DISTANCE), (float)(wtry2+Math.sin(wsa2+(0.5*persistence_length))*SEARCH_DISTANCE), (float)wtry2}, PolygonRoi.POLYGON);
					
					// find all possible matches
					for(int j = 0; j < lines.size(); ++j)
					{
						// skip self
						if(j == i) continue;
						
						// get other line
						Line ol = lines.get(j);
						
						// again, skip if less than some size
						if(ol.size() < 2+AVERAGING_DISTANCE+2*BACKTRACK_DISTANCE) continue;
						
						// get coordinates of points
						int olx1 = ol.get(0+BACKTRACK_DISTANCE).px;
						int oly1 = ol.get(0+BACKTRACK_DISTANCE).py;
						int olx2 = ol.get(ol.size()-BACKTRACK_DISTANCE-1).px;
						int oly2 = ol.get(ol.size()-BACKTRACK_DISTANCE-1).py;
						
						// find a match
						int source_window = -1;
						int destination_window = -1;
						if(w1.contains(olx1, oly1))
						{
							System.err.println("found possible match between W1 (" + wx1 + "," + wy1 + ") and OL1 (" + olx1 + "," + oly1 + ")");
							source_window = 1;
							destination_window = 1;
						}
						if(w1.contains(olx2, oly2))
						{
							System.err.println("found possible match between W1 (" + wx1 + "," + wy1 + ") and OL2 (" + olx2 + "," + oly2 + ")");
							source_window = 1;
							destination_window = 2;
						}
						if(w2.contains(olx1, oly1))
						{
							System.err.println("found possible match between W2 (" + wx2 + "," + wy2 + ") and OL1 (" + olx1 + "," + oly1 + ")");
							source_window = 2;
							destination_window = 1;
						}
						if(w2.contains(olx2, oly2))
						{
							System.err.println("found possible match between W2 (" + wx2 + "," + wy2 + ") and OL2 (" + olx2 + "," + oly2 + ")");
							source_window = 2;
							destination_window = 2;
						}
						
						// add to candidate matches if match was found
						if(source_window != -1 && destination_window != -1)
						{
							Tuple<Line, Integer> source = new Tuple<Line, Integer>(l, source_window);
							Tuple<Line, Integer> destination = new Tuple<Line, Integer>(ol, destination_window);
							
							// first direction
							if(candidate_matches.containsKey(source))
							{
								// add to existing set
								HashSet<Tuple<Line, Integer> > source_set = candidate_matches.get(source);
								source_set.add(destination);
								candidate_matches.put(source, source_set);
							}
							else
							{
								// create new entry
								HashSet<Tuple<Line, Integer> > source_set = new HashSet<Tuple<Line, Integer> >();
								source_set.add(destination);
								candidate_matches.put(source, source_set);
							}
							
							// bi-directional
							if(candidate_matches.containsKey(destination))
							{
								// add to existing set
								HashSet<Tuple<Line, Integer> > destination_set = candidate_matches.get(destination);
								destination_set.add(source);
								candidate_matches.put(destination, destination_set);
							}
							else
							{
								// create new entry
								HashSet<Tuple<Line, Integer> > destination_set = new HashSet<Tuple<Line, Integer> >();
								destination_set.add(source);
								candidate_matches.put(destination, destination_set);
							}
						}
					}
				}
				
				// -------------------------------------------------------------
				
				// straight forward linking of matches
				for(Map.Entry<Tuple<Line, Integer>, HashSet<Tuple<Line, Integer> > > entry : candidate_matches.entrySet())
				{
					// get source
					Tuple<Line, Integer> source = entry.getKey();
					
//					// skip if already matched
//					if(matched.contains(source))
//					{
//						continue;
//					}
					
					// get details on source
					Line source_line = source.first;
					Integer source_window = source.second;
					int source_wx = -1, source_wy = -1;
					double source_wsx = -1, source_wsy = -1;
					int source_wdx = -1, source_wdy = -1;
					double source_wdsx = -1, source_wdsy = -1;
					double source_wsa = -1;
					if(source_window % 2 == 1) // front
					{
						source_wx = source_line.get(0+BACKTRACK_DISTANCE).px;
						source_wy = source_line.get(0+BACKTRACK_DISTANCE).py;
						source_wsx = source_line.get(0+BACKTRACK_DISTANCE).sx;
						source_wsy = source_line.get(0+BACKTRACK_DISTANCE).sy;
						source_wdx = source_line.get(0+BACKTRACK_DISTANCE).px - source_line.get(1+BACKTRACK_DISTANCE).px;
						source_wdy = source_line.get(0+BACKTRACK_DISTANCE).py - source_line.get(1+BACKTRACK_DISTANCE).py;
						source_wdsx = (source_line.get(0+BACKTRACK_DISTANCE).px + source_line.get(0+BACKTRACK_DISTANCE).sx) - (source_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).px + source_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).sx);
						source_wdsy = (source_line.get(0+BACKTRACK_DISTANCE).py + source_line.get(0+BACKTRACK_DISTANCE).sy) - (source_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).py + source_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).sy);
						source_wsa = Math.atan2(source_wdsy, source_wdsx);
					}
					else if(source_window % 2 == 0) // back
					{
						source_wx = source_line.get(source_line.size()-BACKTRACK_DISTANCE-1).px;
						source_wy = source_line.get(source_line.size()-BACKTRACK_DISTANCE-1).py;
						source_wsx = source_line.get(source_line.size()-BACKTRACK_DISTANCE-1).sx;
						source_wsy = source_line.get(source_line.size()-BACKTRACK_DISTANCE-1).sy;
						source_wdx = source_line.get(source_line.size()-BACKTRACK_DISTANCE-1).px - source_line.get(source_line.size()-BACKTRACK_DISTANCE-2).px;
						source_wdy = source_line.get(source_line.size()-BACKTRACK_DISTANCE-1).py - source_line.get(source_line.size()-BACKTRACK_DISTANCE-2).py;
						source_wdsx = (source_line.get(source_line.size()-BACKTRACK_DISTANCE-1).px + source_line.get(source_line.size()-BACKTRACK_DISTANCE-1).sx) - (source_line.get(source_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).px + source_line.get(source_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).sx);
						source_wdsy = (source_line.get(source_line.size()-BACKTRACK_DISTANCE-1).py + source_line.get(source_line.size()-BACKTRACK_DISTANCE-1).sy) - (source_line.get(source_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).py + source_line.get(source_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).sy);
						source_wsa = Math.atan2(source_wdsy, source_wdsx);
					}
					
					// get set of destinations
					HashSet<Tuple<Line, Integer> > destinations = entry.getValue();
					
					// skip if more than one match
	//				if(destinations.size() > 1)
	//				{
	//					continue;
	//				}
					
					// find best match
					Tuple<Line, Integer> best_matching_destination = null;
					Line best_matching_line = null;
					int best_matching_window = -1;
					double best_matching_cost = Double.MAX_VALUE;
					for(Tuple<Line, Integer> destination : destinations)
					{
//						// skip if already matched
//						if(matched.contains(destination))
//						{
//							continue;
//						}
						
						// get details on source
						Line destination_line = destination.first;
						Integer destination_window = destination.second;
						int destination_wx = -1, destination_wy = -1;
						double destination_wsx = -1, destination_wsy = -1;
						int destination_wdx = -1, destination_wdy = -1;
						double destination_wdsx = -1, destination_wdsy = -1;
						double destination_wsa = -1;
						if(destination_window % 2 == 1) // front
						{
							destination_wx = destination_line.get(0+BACKTRACK_DISTANCE).px;
							destination_wy = destination_line.get(0+BACKTRACK_DISTANCE).py;
							destination_wsx = destination_line.get(0+BACKTRACK_DISTANCE).sx;
							destination_wsy = destination_line.get(0+BACKTRACK_DISTANCE).sx;
							destination_wdx = destination_line.get(0+BACKTRACK_DISTANCE).px - destination_line.get(1+BACKTRACK_DISTANCE).px;
							destination_wdy = destination_line.get(0+BACKTRACK_DISTANCE).py - destination_line.get(1+BACKTRACK_DISTANCE).py;
							destination_wdsx = (destination_line.get(0+BACKTRACK_DISTANCE).px + destination_line.get(0+BACKTRACK_DISTANCE).sx) - (destination_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).px + destination_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).sx);
							destination_wdsy = (destination_line.get(0+BACKTRACK_DISTANCE).py + destination_line.get(0+BACKTRACK_DISTANCE).sy) - (destination_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).py + destination_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).sy);
							destination_wsa = Math.atan2(destination_wdsy, destination_wdsx);
						}
						else if(destination_window % 2 == 0) // back
						{
							destination_wx = destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-1).px;
							destination_wy = destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-1).py;
							destination_wsx = destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-1).sx;
							destination_wsy = destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-1).sy;
							destination_wdx = destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-1).px - destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-2).px;
							destination_wdy = destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-1).py - destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-2).py;
							destination_wdsx = (destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-1).px + destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-1).sx) - (destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).px + destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).sx);
							destination_wdsy = (destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-1).py + destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-1).sy) - (destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).py + destination_line.get(destination_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).sy);
							destination_wsa = Math.atan2(destination_wdsy, destination_wdsx);
						}
						
						// calculate distance between points
						double dspx = (source_wx + source_wsx) - (destination_wx + destination_wsx);
						double dspy = (source_wy + source_wsy) - (destination_wy + destination_wsy);
						double diff_distance = Math.sqrt(dspx*dspx+dspy*dspy);
							
						// calculate difference in angle of traces
						System.err.println("source_wsa = " + source_wsa*180/Math.PI);
						System.err.println("destination_wsa = " + destination_wsa*180/Math.PI);
						double diff_angle = Math.abs(source_wsa - destination_wsa);
						System.err.println("diff_angle = " + diff_angle*180/Math.PI);
						if(diff_angle > Math.PI)
						{
							diff_angle = (2*Math.PI) - diff_angle;
							System.err.println("new diff_angle = " + diff_angle*180/Math.PI);
						}
						diff_angle = Math.PI - diff_angle; // NOTE: angle difference of 180 is optimal, 0 is -> flip angle so can be used as low cost
						
						// RSLV: move up to candidate match?
						if(diff_angle > Math.toRadians(MAX_BENDING_ANGLE))
						{
							continue; // skip match
						}
						
						double matching_cost = diff_angle;//+0.1*diff_distance; //diff_distance+5*3*SIGMA*diff_angle;//diff_distance+COST_FUNCTION_WEIGHT*diff_angle; // TODO: define cost function based on distance and orientation
						
						// check for better match
						if(matching_cost < best_matching_cost)
						{
							best_matching_destination = destination;
							best_matching_cost = matching_cost;
							best_matching_line = destination_line;
							best_matching_window = destination.second;
						}
					}
					if(best_matching_destination != null && best_matching_line != null && best_matching_window != -1)
					{
						// skip if already matched
	//					if(matched.contains(destination))
	//					{
	//						continue;
	//					}
						
						// update destination:windows2 to source:window2
						for(Map.Entry<Tuple<Line, Integer>, HashSet<Tuple<Line, Integer> > > e : candidate_matches.entrySet())
						{
							Tuple<Line, Integer> k = e.getKey();
							HashSet<Tuple<Line, Integer> > v = e.getValue();
							
							for(Tuple<Line, Integer> d : v)
							{
								if(d.first == best_matching_line && d.second == (best_matching_window == 1 ? 2 : 1))
								{
									d.first = source_line;
									d.second = source_window + 2;
								}
							}
						}
						//source_window += 2;
						
						// connect the two lines together
						for(int i = 0; i < best_matching_line.size(); ++i)
						{
							if(source_window % 2 == 1 && best_matching_window % 2 == 1)
							{
								source_line.addFirst(best_matching_line.get(i));
							}
							else if(source_window % 2 == 1 && best_matching_window % 2 == 0)
							{
								source_line.addFirst(best_matching_line.get(best_matching_line.size()-i-1));
							}
							else if(source_window % 2 == 0 && best_matching_window % 2 == 1)
							{
								source_line.addLast(best_matching_line.get(i));
							}
							else if(source_window % 2 == 0 && best_matching_window % 2 == 0)
							{
								source_line.addLast(best_matching_line.get(best_matching_line.size()-i-1));
							}
						}
						
						// remove line
						lines.remove(best_matching_line);
						
//						// add to list of matched end points
//						matched.add(source);
//						matched.add(best_matching_destination);
						
						changed = true;
						break;
					}
					
					if(changed)
					{
						break;
					}
				}
			}
//		}
		
		// *********************************************************************
		
		// Present interactive vector overlay tool
		if(SHOW_VECTOR_OVERLAY)
		{
			// duplicate image for vector overlay
			Profiling.tic();
			ip_original.setInterpolationMethod(ImageProcessor.NONE); // NONE, NEAREST_NEIGHBOR, BILINEAR, BICUBIC
			ImageProcessor vector_overlay_ip = ip_original.duplicate();
			Overlay eigenvectors_overlay = new Overlay();
			
			// create vector overlay of primary [and secondary] eigenvector on top of *original* image
			for(int py = 0; py < image_height; ++py)
			{
				for(int px = 0; px < image_width; ++px)
				{
					// check filter status of pixel
					boolean filtered = (hit_map_filtered_ip.get(px, py) == 0); // RSLV: hit_count_map_ip.get(px, py) >= LINE_THREHSOLD?
					
					// DEBUG: overlay vector on scaled image
					double cx = px+0.5;
					double cy = py+0.5;
					
					Roi second_eigenvector_roi = new ij.gui.Line(cx-0.4*results_step_3[px][py][4], cy-0.4*results_step_3[px][py][5], cx+0.4*results_step_3[px][py][4], cy+0.4*results_step_3[px][py][5]);
					
					second_eigenvector_roi.setStrokeColor(filtered ? Color.DARK_GRAY : new Color(255, 255, 0, (USE_TRANSPARENCY ? (int)(255*Math.max(0, r_squared_fit_results[px][py])) : 255)));
					second_eigenvector_roi.setStrokeWidth(0.0);
					second_eigenvector_roi.setPosition(0); // NOTE: only applicable to single frame image!
					eigenvectors_overlay.add(second_eigenvector_roi);
					
					if(!filtered)
					{
						Roi second_eigenvector_corrected_roi = new ij.gui.Line(cx-0.4*results_step_3[px][py][4]+fitting_results[px][py][2]*results_step_3[px][py][1], cy-0.4*results_step_3[px][py][5]+fitting_results[px][py][2]*results_step_3[px][py][2], cx+0.4*results_step_3[px][py][4]+fitting_results[px][py][2]*results_step_3[px][py][1], cy+0.4*results_step_3[px][py][5]+fitting_results[px][py][2]*results_step_3[px][py][2]);
				
						second_eigenvector_corrected_roi.setStrokeColor(new Color(0, 255, 0, (USE_TRANSPARENCY ? (int)(255*Math.max(0, r_squared_fit_results[px][py])) : 255)));
						second_eigenvector_corrected_roi.setStrokeWidth(0.0);
						second_eigenvector_corrected_roi.setPosition(0); // RSLV: only applicable to single frame image
						eigenvectors_overlay.add(second_eigenvector_corrected_roi);
					}
					
					/*Roi first_eigenvector_roi = new ij.gui.Line(cx-0.4*results_step_3[px][py][1], cy-.4*results_step_3[px][py][2], cx+0.4*results_step_3[px][py][1], cy+0.4*results_step_3[px][py][2]);
					
					first_eigenvector_roi.setStrokeWidth(0.0);
					first_eigenvector_roi.setStrokeColor(Color.ORANGE);
					first_eigenvector_roi.setPosition(0); // NOTE: redundant in single frame instance
					eigenvectors_overlay.add(first_eigenvector_roi);*/
				}
			}
			
			// show original image with eigenvectors overlay
			final ImagePlus vector_overlay_imp = new ImagePlus("DEBUG: eigenvectors overlay", vector_overlay_ip); // TMP: pass as model to custom MouseAdapter
			//vector_overlay_imp.resetDisplayRange();
			vector_overlay_imp.setOverlay(eigenvectors_overlay);
			//vector_overlay_imp.updateAndRepaintWindow();
			vector_overlay_imp.show();
			
			// RSLV: if this interaction tool will be used in the release; create subclass of MouseAdapter to accept the model as parameters, rather than using final class from parent class!
			final ImageWindow img_wnd = vector_overlay_imp.getWindow(); // final for inner class access
			final ImageCanvas img_cnv = img_wnd.getCanvas(); // final for inner class access
//			final double[][][] odf_map_tmp = odf_map; // TMP: pass as model to custom MouseAdapter
//			final double[][][] odf_filtered_map_tmp = odf_filtered_map; // TMP: pass as model to custom MouseAdapter
//			final double[][][] odf_min_max_map_tmp = odf_min_max_map; // TMP: pass as model to custom MouseAdapter
//			final int[][][] odf_peaks_tmp = odf_peaks; // TMP: pass as model to custom MouseAdapter
//			final double[][][] odf_map2_tmp = odf_map2; // TMP: pass as model to custom MouseAdapter
//			final double[][][] odf_filtered_map2_tmp = odf_filtered_map2; // TMP: pass as model to custom MouseAdapter
//			final double[][][] odf_min_max_map2_tmp = odf_min_max_map2; // TMP: pass as model to custom MouseAdapter
//			final int[][][] odf_peaks2_tmp = odf_peaks2; // TMP: pass as model to custom MouseAdapter
			final double[][][] hessian_results_tmp = results_step_3; // TMP: pass as model to custom MouseAdapter
			final double[][][] fitting_results_tmp = fitting_results; // TMP: pass as model to custom MouseAdapter
			final double[][][] standard_error_fit_results_tmp = standard_error_fit_results; // TMP: pass as model to custom MouseAdapter
			final double[][] chi_squared_fit_results_tmp = chi_squared_fit_results; // TMP: pass as model to custom MouseAdapter
			final double[][] r_squared_fit_results_tmp = r_squared_fit_results; // TMP: pass as model to custom MouseAdapter
			final ImageProcessor ip_tmp = ip; // TMP: pass as model to custom MouseAdapter
			
			//img_cnv.setMagnification(24.0); // 2400% magnification; RSLV: cannot move view after setMagnification?!
			
			img_cnv.addMouseListener(new MouseAdapter(){
			
				/**
				 *	Private members
				 */
				private int previous_x;
				private int previous_y;
				private int current_x;
				private int current_y;
				
				private Roi clicked_pixel = null;
				private ij.gui.Line perpendicular_line = null;
//				private OvalRoi odf_range_one = null;
//				private OvalRoi odf_range_three = null;
				
				private PlotWindow profile_plot_wnd = null;
				private Plot profile_plot = null;
				
				private PolygonRoi trace = null;
				
				/**
				 *	User clicked a pixel in the image
				 */
				@Override
				public void mouseClicked(MouseEvent e)
				{
					// backup previous coordinate
					previous_x = current_x;
					previous_y = current_y;
					
					// update current coordinate
					current_x = img_cnv.offScreenX(e.getX());
					current_y = img_cnv.offScreenY(e.getY());
					
					// draw border around selected pixel
					if(clicked_pixel != null)
					{
						// remove previous ROI first
						vector_overlay_imp.getOverlay().remove(clicked_pixel);
					}
					clicked_pixel = new Roi(current_x, current_y, 1, 1);
					clicked_pixel.setStrokeColor(Color.CYAN);
					clicked_pixel.setStrokeWidth(0.0);
					vector_overlay_imp.getOverlay().add(clicked_pixel);
					
//					// draw ODF ranges
//					if(odf_range_one != null)
//					{
//						// remove previous ROI first
//						vector_overlay_imp.getOverlay().remove(odf_range_one);
//					}
//					odf_range_one = new OvalRoi(current_x - 0.5*ANISOTROPY_FACTOR*SIGMA + 0.5, current_y - 0.5*ANISOTROPY_FACTOR*SIGMA + 0.5, ANISOTROPY_FACTOR*SIGMA, ANISOTROPY_FACTOR*SIGMA);
//					odf_range_one.setStrokeColor(Color.BLUE);
//					odf_range_one.setStrokeWidth(0.0);
//					vector_overlay_imp.getOverlay().add(odf_range_one);
//					
//					if(odf_range_three != null)
//					{
//						// remove previous ROI first
//						vector_overlay_imp.getOverlay().remove(odf_range_three);
//					}
//					odf_range_three = new OvalRoi(current_x - 3*0.5*ANISOTROPY_FACTOR*SIGMA + 0.5, current_y - 3*0.5*ANISOTROPY_FACTOR*SIGMA + 0.5, 3*ANISOTROPY_FACTOR*SIGMA, 3*ANISOTROPY_FACTOR*SIGMA);
//					odf_range_three.setStrokeColor(Color.BLUE);
//					odf_range_three.setStrokeWidth(0.0);
//					vector_overlay_imp.getOverlay().add(odf_range_three);
					
					//vector_overlay_imp.updateAndRepaintWindow(); // force update of window
					
					// determine action:
					// [1] new pixel selected -> show information
					// [2] same pixel selected -> show trace menu
					if(previous_x == current_x && previous_y == current_y)
					{
						// show trace line options dialog
						GenericDialog gd = new GenericDialog("Trace line from pixel x=" + current_x + ", y=" + current_y);
						gd.addNumericField("Max trace length", 100, 0);
						//gd.addNumericField("Step size", 1, 2);
						//gd.hideCancelButton();
						gd.showDialog();
						if(gd.wasOKed()) // !gd.wasCanceled()
						{
							// trace line from selected pixel
							final int MAX_TRACE_LENGTH = (int)gd.getNextNumber();
							//final double STEP_SIZE = gd.getNextNumber();
							
							double trace_length = 0.0;
							Vector<Double> trace_xs_vec = new Vector<Double>();
							Vector<Double> trace_ys_vec = new Vector<Double>();
							
							// first position and direction
							double tx = current_x + 0.5;
							double ty = current_y + 0.5;
							
							// mu correction (first eigenvector)
							// RSLV: repeat m correction until converged to same pixel?
							tx += hessian_results_tmp[(int)tx][(int)ty][1] * fitting_results_tmp[(int)tx][(int)ty][2];
							ty += hessian_results_tmp[(int)tx][(int)ty][2] * fitting_results_tmp[(int)tx][(int)ty][2];
							
							// add current point to trace
							trace_xs_vec.add(tx);
							trace_ys_vec.add(ty);
							
							// trace path along centerline pixels
							boolean tracing = true;
							double STEP_SIZE = 0.1;
							
							while(tracing && trace_length < MAX_TRACE_LENGTH)
							{
								// get second eigenvector of pixel
								int cx = (int)tx;
								int cy = (int)ty;
								double sx = hessian_results_tmp[cx][cy][4]; // TODO: solve ArrayIndexOutOfBounds when running off the image!
								double sy = hessian_results_tmp[cx][cy][5]; // TODO: solve ArrayIndexOutOfBounds when running off the image!
								double sa = Math.atan2(sy, sx);
								double asa = Math.abs(sa);
								
								// get linear interpolated vector field
								double isa = 0.0;
								double wsa = 0;
								for(int kx = -1; kx <= 1; ++kx)
								{
									for(int ky = -1; ky <= 1; ++ky)
									{
										// get neighbour pixel
										int nx = cx + kx;
										int ny = cy + ky;
										if(nx < 0 || nx >= ip_tmp.getWidth() || ny < 0 || ny >= ip_tmp.getHeight())
										{
											continue; // skip out of image
										}
										
										// get neighbour distance
										double dx = tx - (nx+0.5);
										double dy = ty - (ny+0.5);
										if(Math.abs(dx) >= 1.0 || Math.abs(dy) >= 1.0)
										{
											continue; // skip out of range
										}
										
										// get neighbour information
										double nsx = hessian_results_tmp[nx][ny][4];
										double nsy = hessian_results_tmp[nx][ny][5];
										double nsa = Math.atan2(nsy, nsx);
										double ansa = Math.abs(nsa);
										
										// add sum to interpolated value
										double idx = 1.0 - dx;
										double idy = 1.0 - dy;
										double w = Math.sqrt(idx*idx+idy*idy);
										
										// add to sum
										isa += w * ansa;
										wsa += w;
									}
								}
								isa /= wsa;
								
								// take step in interpolated direction
								tx += STEP_SIZE * Math.cos(isa);
								ty += STEP_SIZE * Math.sin(isa);
								
								if(tx >= 0 && tx < ip_tmp.getWidth() && ty >= 0 && ty < ip_tmp.getHeight())
								{
									// add point to line
									trace_xs_vec.add(tx);
									trace_ys_vec.add(ty);
									
									// increase trace length
									trace_length += STEP_SIZE;
								}
								else
								{
									// stop tracing (outside image window)
									System.err.println("Stop tracing line, outside image window");
									tracing = false;
								}
							}
							
							/*
							boolean tracing = true;
							while(tracing && trace_length < MAX_TRACE_LENGTH)
							{
								// get second eigenvector of pixel
								double sx = hessian_results_tmp[(int)tx][(int)ty][4]; // TODO: solve ArrayIndexOutOfBounds when running off the image!
								double sy = hessian_results_tmp[(int)tx][(int)ty][5]; // TODO: solve ArrayIndexOutOfBounds when running off the image!
								double sa = Math.atan2(sy, sx);
								
								System.err.println("sa = " + sa + ", " + (180.0/Math.PI*sa));
									
								// TODO: reorient s(t) vector in direction of trace
								
								// determine distance to next pixels
								double dx = Math.floor(tx)-tx;
								double dy = Math.floor(ty)-ty;
								if(Math.abs(sa) < Math.PI/2)
								{
									System.err.println("true");
									dx = dx + 1.0;
								}
								if(Math.abs(sa) > 0)
								{
									System.err.println("false");
									dy = dy + 1.0;
								}
								
								System.err.println("tx = " + tx + ", ty = " + ty);
								System.err.println("dx = " + dx + ", dy = " + dy);
								
								// calculate remainging sides of triangles
								double dxp = dy*(Math.PI/2-sa);
								double dyp = dx*sa;
								
								System.err.println("dxp = " + dxp + ", dyp = " + dyp);
								
								// take triangle with smallest opposing side
								if(dxp < dyp)
								{
									System.err.println("true");
									dyp = dy;
								}
								else
								{
									System.err.println("false");
									dxp = dx;
								}
								
								System.err.println("dxp = " + dxp + ", dyp = " + dyp);
								
								// backup previous position
								//double prev_tx = tx;
								//double prev_ty = ty;
								
								// perform step
								tx += dxp;
								ty += dyp;
								
								System.err.println("ntx = " + tx + ", nty = " + ty);
								
								// check if still within image window
								if(tx >= 0 && tx < ip_tmp.getWidth() && ty >= 0 && ty < ip_tmp.getWidth())
								{
									System.err.println("true");
									// add point to line
									trace_xs_vec.add(tx);
									trace_ys_vec.add(ty);
									//trace_length += Math.sqrt(dxp*dxp+dyp*dyp);
									
								}
								else
								{
									System.err.println("false");
									// stop tracing (outside image window)
									System.err.println("Stop tracing line, outside image window");
									tracing = false;
								}
								
								trace_length += 1; // DEBUG; TODO: remove
							}
							*/
							
							// TODO: opposite direction
							
							// -------------------------------------------------
							
							// manually convert Float collection to array of floats
							float[] trace_xs = new float[trace_xs_vec.size()];
							float[] trace_ys = new float[trace_ys_vec.size()];
							for(int i = 0; i < trace_xs_vec.size(); ++i)
							{
								if(trace_xs_vec.get(i) != null && trace_ys_vec.get(i) != null)
								{
									trace_xs[i] = trace_xs_vec.get(i).floatValue();
									trace_ys[i] = trace_ys_vec.get(i).floatValue();
								}
								else
								{
									// RSLV: how to deal with null pointers if they exist
								}
							}
							
							// display trace line in image overlay
							if(trace != null)
							{
								// remove previous trace first
								vector_overlay_imp.getOverlay().remove(trace);
							}
							trace = new PolygonRoi(trace_xs, trace_ys, Roi.POLYLINE);
							trace.setStrokeColor(Color.RED);
							trace.setStrokeWidth(0.0);
							vector_overlay_imp.getOverlay().add(trace);
							//vector_overlay_imp.updateAndRepaintWindow(); // force update of window
						}
					}
					else
					{
						// display line profile perpendicular to curve
						if(perpendicular_line != null)
						{
							// remove previous trace first
							vector_overlay_imp.getOverlay().remove(perpendicular_line);
						}
						double nx = hessian_results_tmp[current_x][current_y][1];
						double ny = hessian_results_tmp[current_x][current_y][2];
						perpendicular_line = new ij.gui.Line(current_x+0.5-nx*3*SIGMA, current_y+0.5-ny*3*SIGMA, current_x+0.5+nx*3*SIGMA, current_y+0.5+ny*3*SIGMA);
						perpendicular_line.setStrokeColor(Color.RED);
						perpendicular_line.setStrokeWidth(0.0);
						vector_overlay_imp.getOverlay().add(perpendicular_line);
						//vector_overlay_imp.updateAndRepaintWindow(); // force update of window
						
						// TODO: display information about pixel
						// raw data [intensity, Hessian, fitting]
						// image region and parameters
						
						// populate line profile and Gaussian curve with data (NOTE: adapted from code of fitting procedure)
						ip_tmp.setInterpolationMethod(INTERPOLATION_METHOD); // NONE, NEAREST_NEIGHBOR, BILINEAR, BICUBIC, 
						int line_profile_width = (int)Math.ceil(3*SIGMA);
						int data_points = 2*SAMPLE_RATE*line_profile_width+1;
						double[] plot_xs = new double[data_points];
						double[] profile_ys = new double[data_points];
						double[] gaussian_ys = new double[data_points];
						double min_value = Double.MAX_VALUE; // keep track of minimum value
						double max_value = Double.MIN_VALUE; // keep track of maximum value
						double first_order_moment = 0.0;
						double second_order_moment = 0.0;
						double[] fit_param = fitting_results_tmp[current_x][current_y];
						double ii = -line_profile_width;
						for(int i = 0; i < data_points; ++i)
						{
							// interpolated x,y
							double ix = current_x + ii * nx;
							double iy = current_y + ii * ny;
							double pv = ip_tmp.getPixelInterpolated(ix, iy);
							plot_xs[i] = ii; // NOTE use relative x-coordinate!!
							profile_ys[i] = pv;
							
							double gv = fit_param[0] + fit_param[1] * Math.exp(-0.5 * ((fit_param[2] - ii) * (fit_param[2] - ii))/(fit_param[3] * fit_param[3])); // RSLV: normalised Gaussian? (1/(Math.sqrt(2*Math.PI)*fit_param[3])) *
							gaussian_ys[i] = gv;
							
							// update min/max value
							min_value = Math.min(min_value, pv);
							max_value = Math.max(max_value, pv);
							min_value = Math.min(min_value, gv);
							max_value = Math.max(max_value, gv);
							
							first_order_moment += pv;
							second_order_moment += pv*pv;
							
							// increment ii
							ii += SAMPLE_RATE_INV;
						}
						
						double mean_value = first_order_moment / data_points;
						double variance_value = (second_order_moment / data_points) - (mean_value * mean_value);
						
						double chi_square = 0.0;
						double chi_square_variance = 0.0;
						double chi_square_sum_square = 0.0;
						double chi_square_expected = 0.0;
						for(int i = 0; i < data_points; ++i)
						{
							double dy = profile_ys[i] - gaussian_ys[i];
							dy = dy * dy;
							chi_square += dy;
							chi_square_variance += dy / variance_value;
							chi_square_sum_square += dy / second_order_moment;
							chi_square_expected += dy / gaussian_ys[i];
						}
						
						// add fit and profile data to plot
						profile_plot = new Plot("Profile plot @ (x=" + current_x + ", y=" + current_y + ")", "Relative X-position", "Image intensity", plot_xs, gaussian_ys);
						profile_plot.setLimits(-line_profile_width-1, line_profile_width+1, Math.min(0.0, min_value*1.05), max_value*1.05); // RSLV: bottom-y to zero?
						profile_plot.setSize(600, 600);
						profile_plot.setFrameSize(600, 600);
						
						profile_plot.setColor(Color.RED);
						profile_plot.addPoints(plot_xs, profile_ys, Plot.CIRCLE); // LINE, DOT, CROSS, CIRCLE, BOX, TRIANGLE
						
						// add labels (reversed order)
						profile_plot.setColor(Color.RED);
						profile_plot.addLabel(0.02, 0.05, " o   Line profile");
						profile_plot.setColor(Color.BLUE);
						profile_plot.addLabel(0.02, 0.08, "---  Gaussian fit\n       bg  = " + String.format("%.2f", fit_param[0]) + "\n       amp = " + String.format("%.2f", fit_param[1]) + "\n       mu  = " + String.format("%.4f", fit_param[2]) + "\n       sig = " + String.format("%.4f", fit_param[3]) + "\n       chi = " + String.format("%.1f", chi_squared_fit_results_tmp[current_x][current_y]) + "\n       log = " + String.format("%.5f", Math.log(chi_squared_fit_results_tmp[current_x][current_y])) + "\n       R^2 = " + String.format("%.4f", r_squared_fit_results_tmp[current_x][current_y]) + "\n\n       ecc = " + String.format("%.2f", Math.sqrt((hessian_results_tmp[current_x][current_y][3]*hessian_results_tmp[current_x][current_y][3]) - (hessian_results_tmp[current_x][current_y][0]*hessian_results_tmp[current_x][current_y][0]))) + "\n       mean = " + String.format("%.4f", mean_value) + "\n       var = " + String.format("%.4f", variance_value)+ "\n       lvar = " + String.format("%.4f", Math.log(variance_value)) + "\n       chi = " + String.format("%.1f", chi_square)+ "\n       chi_log = " + String.format("%.4f", Math.log(chi_square)) + "\n       chi_var = " + String.format("%.4f", chi_square_variance)+ "\n       chi_sum = " + String.format("%.4f", chi_square_sum_square) + "\n       chi_exp = " + String.format("%.4f", chi_square_expected));
						
						// create new plot window or draw in existing plot window
						if(profile_plot_wnd != null && !profile_plot_wnd.isClosed())
						{
							profile_plot_wnd.drawPlot(profile_plot);
							profile_plot_wnd.setTitle("Profile plot @ (x=" + current_x + ", y=" + current_y + ")");
							// RSLV: set focus on plot window?
						}
						else
						{
							profile_plot_wnd = profile_plot.show();
							profile_plot_wnd.setLocation(700, 80);
						}
					}
					
					// force update on window
					vector_overlay_imp.updateAndRepaintWindow();
				}
			});
			
			Profiling.toc("DEBUG: drawing eigenvector overlay on duplicate of original image");
		}
		
		// *********************************************************************
		
		// Display table with results
		if(SHOW_RESULTS_TABLE)
		{
			Profiling.tic();
			ResultsTable raw_data_table = new ResultsTable();
			raw_data_table.setPrecision(5);
			raw_data_table.showRowNumbers(false);
			//raw_data_table.reset(); // to clear the table
			
			// TODO: fill table
			for(int py = 0; py < image_height; ++py)
			{
				for(int px = 0; px < image_width; ++px)
				{
					// new row of results
					raw_data_table.incrementCounter();
					
					// pixel coordinate and intensity value
					raw_data_table.addValue("px", px);
					raw_data_table.addValue("py", py);
					raw_data_table.addValue("pv", ip_original.get(px, py));
					
					// gradients
					/*
					raw_data_table.addValue("dx", dx.getf(px, py));
					raw_data_table.addValue("dy", dy.getf(px, py));
					raw_data_table.addValue("dxdx", dxdx.getf(px, py));
					raw_data_table.addValue("dxdy", dxdy.getf(px, py));
					raw_data_table.addValue("dydx", dydx.getf(px, py));
					raw_data_table.addValue("dydy", dydy.getf(px, py));*/
					
					// Hessian matrix
					raw_data_table.addValue("L1", results_step_3[px][py][0]); // L1
					raw_data_table.addValue("V1x", results_step_3[px][py][1]); // V1x
					raw_data_table.addValue("V1y", results_step_3[px][py][2]); // V1y
					raw_data_table.addValue("L2", results_step_3[px][py][3]); // L2
					raw_data_table.addValue("V2x", results_step_3[px][py][4]); // V2x
					raw_data_table.addValue("V2y", results_step_3[px][py][5]); // V2y
					
					// Frangi measures
					raw_data_table.addValue("blobness", frangi_measures[px][py][2]);
					raw_data_table.addValue("structureness", frangi_measures[px][py][3]);
					raw_data_table.addValue("vesselness", frangi_measures[px][py][4]);
					
					// fitted parameters of line profile fit
					raw_data_table.addValue("bg", fitting_results[px][py][0]); // bg
					raw_data_table.addValue("amp", fitting_results[px][py][1]); // amp
					raw_data_table.addValue("mu", fitting_results[px][py][2]); // mu
					raw_data_table.addValue("sig", fitting_results[px][py][3]); // sigma
					
					// standard error of line profile fit
					raw_data_table.addValue("SE_bg", standard_error_fit_results[px][py][0]); // bg
					raw_data_table.addValue("SE_amp", standard_error_fit_results[px][py][1]); // amp
					raw_data_table.addValue("SE_mu", standard_error_fit_results[px][py][2]); // mu
					raw_data_table.addValue("SE_sig", standard_error_fit_results[px][py][3]); // sigma
					
					// goodness of line profile fit
					raw_data_table.addValue("chi2", chi_squared_fit_results[px][py]); // chi^2
					raw_data_table.addValue("R2", r_squared_fit_results[px][py]); // R^2
					
					// filter status
					raw_data_table.addValue("hit", hit_map_filtered_ip.get(px, py)); // hit map
					raw_data_table.addValue("hit_count", hit_count_map_ip.get(px, py)); // hit count map
				}
			}
			
			raw_data_table.show("Result of segmentation");
			Profiling.toc("Generating results table");
		}
		
		// *********************************************************************
		
		// filter short line segments
//		if(FILTER_SHORT_LINE_SEGMENTS)
//		{
			Vector<Line> new_lines = new Vector<Line>();
			for(Line l : lines)
			{
				if(l != null && l.contourLength() >= LINE_LENGTH_THRESHOLD)
				{
					new_lines.add(l);
				}
			}
			lines = new_lines;
//		}
		
		// *********************************************************************
		
		// draw overlay on duplicate of original image
		RoiManager roi_manager = RoiManager.getInstance();
		if(roi_manager == null) roi_manager = new RoiManager();
		//roi_manager.show(); // .setVisible(true)?
		ImageProcessor overlay_ip = ip.duplicate();
//		ImageProcessor overlay2_ip = valid_line_points_magnitude_ip.duplicate();
//		overlay2_ip.abs(); // RSLV: to abs() or not?
		Overlay lines_overlay = new Overlay();
		
		ImageProcessor pixelated_traces_ip = new ByteProcessor(image_width, image_height);
		
		for(int px = 0; px < image_width; ++px)
		{
			for(int py = 0; py < image_height; ++py)
			{
				if(valid_line_points_map_ip.get(px, py) == 255)
				{
					ij.gui.Line s_vector = new ij.gui.Line(px + 0.5 - 0.3 * results_step_3[px][py][4], py + 0.5 - 0.3 * results_step_3[px][py][5], px + 0.5 + 0.3 * results_step_3[px][py][4], py + 0.5 + 0.3 * results_step_3[px][py][5]);
					s_vector.setStrokeColor(Color.GRAY);
					s_vector.setStrokeWidth(0.0);
					lines_overlay.add(s_vector);
				}
			}
		}
		int line_color_index = 1;
		for(Line l : lines)
		{
			// skip empty lines (RSLV: or less than 3 line points!)
			if(l == null || l.size() < 3) continue;
			
			// add origin marker to overlay
			OvalRoi origin_p = new OvalRoi(l.getFirst().px + 0.375, l.getFirst().py + 0.375, 0.25, 0.25);
			origin_p.setStrokeColor(Color.RED);
			origin_p.setStrokeWidth(0.0);
			lines_overlay.add(origin_p);
			
			OvalRoi origin_s = new OvalRoi(l.getFirst().px + 0.5 + l.getFirst().sx - 0.125, l.getFirst().py + 0.5 + l.getFirst().sy - 0.125, 0.25, 0.25);
			origin_s.setStrokeColor(Color.BLUE);
			origin_s.setStrokeWidth(0.0);
			lines_overlay.add(origin_s);
			
			// add trace path to overlay
			float[] pxs = new float[l.size()];
			float[] pys = new float[l.size()];
			float[] sxs = new float[l.size()];
			float[] sys = new float[l.size()];
			float[] dxs = new float[l.size()];
			float[] dys = new float[l.size()];
			
			int index = 0;
			for(Point p : l)
			{
				pixelated_traces_ip.set(p.px, p.py, line_color_index);
				pxs[index] = (float)(p.px + 0.5);
				pys[index] = (float)(p.py + 0.5);
				dxs[index] = (float)(p.px + 0.5 + results_step_3[p.px][p.py][6]); // NOTE: use Steger' super-resolved coordinate
				dys[index] = (float)(p.py + 0.5 + results_step_3[p.px][p.py][7]); // NOTE: use Steger' super-resolved coordinate
				sxs[index] = (float)(p.px + 0.5 + p.sx); // NOTE: use super-resolved coordinate
				sys[index] = (float)(p.py + 0.5 + p.sy); // NOTE: use super-resolved coordinate
				++index;
			}
			
			PolygonRoi polyline_p = new PolygonRoi(pxs, pys, Roi.POLYLINE);
			polyline_p.setStrokeColor(Color.YELLOW); // TODO: rainbow colors :D
			polyline_p.setStrokeWidth(0.0);
			lines_overlay.add(polyline_p);
			
			PolygonRoi polyline_d = new PolygonRoi(dxs, dys, Roi.POLYLINE);
			polyline_d.setStrokeColor(Color.ORANGE); // TODO: rainbow colors :D
			polyline_d.setStrokeWidth(0.0);
			lines_overlay.add(polyline_d);
			
			PolygonRoi polyline_s = new PolygonRoi(sxs, sys, Roi.POLYLINE);
			polyline_s.setStrokeColor(Color.RED); // TODO: rainbow colors :D
			polyline_s.setStrokeWidth(0.0);
			lines_overlay.add(polyline_s);
			
////		roi_manager.addRoi(polyline_s); // RSLV: use add(imp, roi, index)?
			
			++line_color_index;
		}
		
//		if(DEBUG_MODE_ENABLED)
//		{
			ImagePlus overlay_imp = new ImagePlus("Detected line traces", overlay_ip);
			overlay_imp.setOverlay(lines_overlay);
			//overlay.updateAndRepaintWindow();
//			overlay_imp.show();
			
//			ImagePlus overlay2_imp = new ImagePlus("Steger's algorithm line traces", overlay2_ip);
//			overlay2_imp.setOverlay(lines_overlay);
//			//overlay.updateAndRepaintWindow();
//			overlay2_imp.show();
			
			LUT rainbow_lut = ConnectedComponents.getConnectedComponentLUT();
			pixelated_traces_ip.setLut(rainbow_lut);
			ImagePlus pixelated_traces_imp = new ImagePlus("Detected line traces", pixelated_traces_ip);
		pixelated_traces_imp.show();
//		}
		
		// show ROI manager
		roi_manager.setVisible(true);
		roi_manager.toFront();
		
		// calculate measures
		ResultsTable raw_data_table = ResultsTable.getResultsTable();
		raw_data_table.setPrecision(5);
		raw_data_table.showRowNumbers(false);
		System.err.println("===== LINE MEASURES =====");
		for(int i = 0; i < lines.size(); ++i)
		{
			Line l = lines.get(i);
			System.err.println("Line " + i);
			System.err.println("L_c = " + l.contourLength());
			System.err.println("R-R = " + l.endToEndDistance());
			System.err.println("L_p = " + l.persistenceLength());
			
			raw_data_table.incrementCounter();
			raw_data_table.addValue("S#", SLICE_INDEX);
			raw_data_table.addValue("L#", i);
			raw_data_table.addValue("Lc", l.contourLength());
			raw_data_table.addValue("Lr", l.endToEndDistance());
			raw_data_table.addValue("Lp", l.persistenceLength());
		}
		raw_data_table.show("Results");
		
		return pixelated_traces_ip;
	}
	
	// *************************************************************************
	
	public static double getAngle(double dx, double dy)
	{
		//return Math.atan2(dy, dx) + Math.PI; // RADIAN MODE
		double angle = (180.0/Math.PI)*Math.atan2(dy, dx);
		if(angle < 0)
		{
			angle = 360.0 + angle;
		}
		return angle; // DEGREE MODE
	}
	
	public static int getOrientationIndex(double theta)
	{
		/* RADIAN MODE
		double theta_c = (theta + (Math.PI / 8));
		if(theta_c >= 2 * Math.PI)
		{
			theta_c -= 2 * Math.PI;
		}
		int orientation_index = (int)Math.floor(theta_c / (Math.PI / 4));
		return orientation_index;
		*/
		
		/* DEGREE MODE */
		double theta_c = (theta + (180.0 / 8));
		if(theta_c >= 2 * 180.0)
		{
			theta_c -= 2 * 180.0;
		}
		int orientation_index = (int)Math.floor(theta_c / (180.0 / 4));
		return orientation_index;
	}
	
	public static int[] getNextPixel(int orientation_index)
	{
		// RSLV: find a mathematical expression for the switch statement?
		int[] dp = new int[2]; // [0]=dx, [1]=dy
		switch(orientation_index)
		{
			case 0:
				dp[0] = 1;
				dp[1] = 0;
				break;
			case 1:
				dp[0] = 1;
				dp[1] = 1;
				break;
			case 2:
				dp[0] = 0;
				dp[1] = 1;
				break;
			case 3:
				dp[0] = -1;
				dp[1] = 1;
				break;
			case 4:
				dp[0] = -1;
				dp[1] = 0;
				break;
			case 5:
				dp[0] = -1;
				dp[1] = -1;
				break;
			case 6:
				dp[0] = 0;
				dp[1] = -1;
				break;
			case 7:
				dp[0] = 1;
				dp[1] = -1;
				break;
		}
		return dp;
	}
	
	public static double gaussianPDF(double mu, double sigma, double x)
	{
		return (1/(sigma*Math.sqrt(2*Math.PI))) * Math.exp(-0.5f * ((x-mu)*(x-mu)/(sigma*sigma)));
	}
	
	public static double gaussianPDF(double sigma, double x)
	{
		return (1/(sigma*Math.sqrt(2*Math.PI))) * Math.exp(-0.5f * ((x*x)/(sigma*sigma)));
	}
}
