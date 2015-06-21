
package plugins;

// import Java classes
import java.util.Collections;
import java.util.Comparator;
import java.util.Vector;
import java.util.HashMap;

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
//import filters.FillHoles;
import filters.ConnectedComponents;
import utils.ImageArithmetic;
import utils.Profiling;

import core.Point;
import core.Line;

//import utils.Tuple;
//import utils.Triple;
import utils.Quadruple;


/**
 *	Microtubule tracing algorithm
 */
public class Trace_Microtubules_ODFs implements PlugIn
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
	public static int VALID_LINE_POINT_THRESHOLD = 2;
	public static boolean SKELETONIZE_POINTS_MASK = true;
	public static boolean FILTER_POINTS_MASK = true;
	
	// ODF map
	public static double ANGULAR_RESOLUTION = 6.0; // degrees per orientation step; NOTE: should be a divisor of 180.0
	public static int ODF_STEPS = (int)(180.0 / ANGULAR_RESOLUTION);
	public static double ANISOTROPY_FACTOR = 5;
	public static int ODF_PEAK_WIDTH = 5;
	
	// debug mode
	public static boolean DEBUG_MODE_ENABLED = true;
	public static boolean SHOW_RESULTS_TABLE = false;
	public static boolean SHOW_VECTOR_OVERLAY = true;
	public static boolean SHOW_ODF_OVERLAY = true;
	public static boolean USE_VECTOR_NORMALS = true;
	public static boolean USE_TRANSPARENCY = true;
	
	/**
	 *
	 */
	public Trace_Microtubules_ODFs()
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
		GenericDialog gd = new GenericDialog("Trace microtubules ODFs");
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
		
		gd.addNumericField("Valid_line_point_threshold", Prefs.get("mt_trace.valid_line_point_threshold", VALID_LINE_POINT_THRESHOLD), 0);
		gd.addCheckbox("Skeletonize_point_map", Prefs.get("mt_trace.skeletonize_point_mask", SKELETONIZE_POINTS_MASK));
		gd.addCheckbox("Filter_point_map", Prefs.get("mt_trace.filter_point_mask", FILTER_POINTS_MASK));
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addNumericField("Angular_resolution", Prefs.get("mt_trace.angular_resolution", ANGULAR_RESOLUTION), 1);
		gd.addNumericField("Anisotropy_factor", Prefs.get("mt_trace.anisotropy_factor", ANISOTROPY_FACTOR), 1);
		gd.addNumericField("ODF_peak_width", Prefs.get("mt_trace.odf_peak_width", ODF_PEAK_WIDTH), 0);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Enable_debug_mode", Prefs.get("mt_trace.debug_mode", DEBUG_MODE_ENABLED));
		gd.addCheckbox("Show_results_table", Prefs.get("mt_trace.show_results_table", SHOW_RESULTS_TABLE));
		gd.addCheckbox("Show_vector_overlay", Prefs.get("mt_trace.show_vector_overlay", SHOW_VECTOR_OVERLAY));
		gd.addCheckbox("Show_ODF_overlay", Prefs.get("mt_trace.show_odf_overlay", SHOW_ODF_OVERLAY));
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
		
		VALID_LINE_POINT_THRESHOLD = (int)gd.getNextNumber();
		SKELETONIZE_POINTS_MASK = gd.getNextBoolean();
		FILTER_POINTS_MASK = gd.getNextBoolean();
		
		ANGULAR_RESOLUTION = gd.getNextNumber();
		ODF_STEPS = (int)(180.0 / ANGULAR_RESOLUTION); // NOTE: also update ODF_STEPS
		ANISOTROPY_FACTOR = gd.getNextNumber();
		ODF_PEAK_WIDTH = (int)gd.getNextNumber();
		
		DEBUG_MODE_ENABLED = gd.getNextBoolean();
		SHOW_RESULTS_TABLE = gd.getNextBoolean();
		SHOW_VECTOR_OVERLAY = gd.getNextBoolean();
		SHOW_ODF_OVERLAY = gd.getNextBoolean();
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
		
		Prefs.set("mt_trace.valid_line_point_threshold", VALID_LINE_POINT_THRESHOLD);
		Prefs.set("mt_trace.skeletonize_point_map", SKELETONIZE_POINTS_MASK);
		Prefs.set("mt_trace.filter_point_map", FILTER_POINTS_MASK);
		
		Prefs.set("mt_trace.angular_resolution", ANGULAR_RESOLUTION);
		Prefs.set("mt_trace.anisotropy_factor", ANISOTROPY_FACTOR);
		Prefs.set("mt_trace.odf_peak_width", ODF_PEAK_WIDTH);
		
		Prefs.set("mt_trace.debug_mode", DEBUG_MODE_ENABLED);
		Prefs.set("mt_trace.show_results_table", SHOW_RESULTS_TABLE);
		Prefs.set("mt_trace.show_vector_overlay", SHOW_VECTOR_OVERLAY);
		Prefs.set("mt_trace.show_odf_overlay", SHOW_ODF_OVERLAY);
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
		ImagePlus imp_out = IJ.createImage("Microtubule traces of " + imp.getTitle(), imp.getWidth(), imp.getHeight(), imp.getNSlices(), 32);
		
		// process slices in image stack
		ImageStack stack_in = imp.getStack();
		ImageStack stack_out = imp_out.getStack();
		for(int slice = 1; slice <= stack_in.getSize(); ++slice)
		{
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
		// RLSV: use ImageStatistics.getStatistics(ImageProcessor ip, int mOptions, Calibration cal)
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
		//if(FIT_ON_RAW_IMAGE) // TODO: implement choice
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
				//double theta = Math.atan2(ny, nx); // RSLV: add rotation as parameter for more optimal fit than Hessian can produce?
				
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
		
		// Step 5a: additional filtering of point (determine which pixels are microtubules and which belong to the background
		Profiling.tic();
		IJ.showStatus("Filter line points and creating projection maps");
		
		// projection maps
		//	[0] = count
		//	[1] = r-squared sum
		//	[2] = r-squared avg
		//	[3] = 
		//double[][][] projection_maps = new double[image_width][image_height][6];
		
		// projection maps; TODO: average other measures, e.g. eigenvectors?
		ImageProcessor hit_map_ip = new ByteProcessor(image_width, image_height);
		ImageProcessor hit_count_map_ip = new ByteProcessor(image_width, image_height); // count of hits
		ImageProcessor r_squared_sum_projection_map_ip = new FloatProcessor(image_width, image_height); // NOTE: sum of R^2 measures
		
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
		double[][][] avg_frangi_measures = new double[image_width][image_height][5];
		// RSLV: condition |L1| <= |L2| may change by averging?
		
		
		double[][][] avg_fitting_results = new double[image_width][image_height][4]; // [x][y][bg=0|amp=1|mu=2|sigma=3]
		double[][][] avg_standard_error_fit_results = new double[image_width][image_height][4]; // [x][y][bg=0|amp=1|mu=2|sigma=3]
		double[][] avg_chi_squared_fit_results = new double[image_width][image_height];
		double[][] avg_r_squared_fit_results = new double[image_width][image_height];
		
		
		int hit_count_max = 0; // keep track of maximum hit count
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
				if(FILTER_ON_SIGMA && (fitting_results[px][py][3] < SIGMA_RATIO_LOWER_THRESHOLD * SIGMA || fitting_results[px][py][3] > SIGMA_RATIO_UPPER_THRESHOLD * SIGMA)) // RSLV: use scale space sigma? Or fixed SIGMA?
				{
					// probably not a valid center line pixel; skip
					filtered = true; //continue;
				}
				
				// filter on r-squared
				if(FILTER_ON_R_SQUARED && (Double.isNaN(r_squared_fit_results[px][py]) || Math.abs(r_squared_fit_results[px][py]) < R_SQUARED_THRESHOLD))
				{
					// probably not a valid center line pixel; skip
					filtered = true; //continue;
				}
				
				// set pixel in hit images
				if(!filtered)
				{
					// set hit image pixel
					hit_map_ip.set(px, py, 255);
					
					// correct for peak position in coordinate
					int cpx = (int)(px+0.5+fitting_results[px][py][2]*results_step_3[px][py][1]);
					int cpy = (int)(py+0.5+fitting_results[px][py][2]*results_step_3[px][py][2]);
					
					// bounds checking
					if(cpx >= 0 && cpx < image_width && cpy >= 0 && cpy < image_height)
					{
						// update hit count
						int cc = 1 + hit_count_map_ip.get(cpx, cpy); // NOTE: getPixel for bounds checking!
						hit_count_map_ip.set(cpx, cpy, cc);
						
						// retain maximum hit count value
						if(cc > hit_count_max)
						{
							hit_count_max = cc;
						}
						
						// calculate tensor orientation, range [0..180] degree
						double theta_n = getAngle(results_step_3[px][py][1], results_step_3[px][py][2]);
						double theta_s = getAngle(results_step_3[px][py][4], results_step_3[px][py][5]);
						theta_n = theta_n % 180.0;
						theta_s = theta_s % 180.0;
						theta_n = theta_n * (Math.PI / 180.0);
						theta_s = theta_s * (Math.PI / 180.0);
						
						// update sums for averages
						avg_results_step_3[cpx][cpy][0] += results_step_3[px][py][0];
						avg_results_step_3[cpx][cpy][1] += Math.cos(theta_n); // RSLV: average n_x
						avg_results_step_3[cpx][cpy][2] += Math.sin(theta_n); // RSLV: average n_y
						avg_results_step_3[cpx][cpy][3] += results_step_3[px][py][3];
						avg_results_step_3[cpx][cpy][4] += Math.cos(theta_s); // RSLV: average s_x
						avg_results_step_3[cpx][cpy][5] += Math.sin(theta_s); // RSLV: average s_y
						avg_results_step_3[cpx][cpy][6] += results_step_3[px][py][6];
						avg_results_step_3[cpx][cpy][7] += results_step_3[px][py][7];
						
						avg_frangi_measures[cpx][cpy][0] = frangi_measures[px][py][0];
						avg_frangi_measures[cpx][cpy][1] = frangi_measures[px][py][1];
						avg_frangi_measures[cpx][cpy][2] = frangi_measures[px][py][2];
						avg_frangi_measures[cpx][cpy][3] = frangi_measures[px][py][3];
						avg_frangi_measures[px][py][4] = frangi_measures[px][py][4];
						
						avg_fitting_results[cpx][cpy][0] += fitting_results[px][py][0]; // BG
						avg_fitting_results[cpx][cpy][1] += fitting_results[px][py][1]; // AMP
						avg_fitting_results[cpx][cpy][2] += fitting_results[px][py][2]; // MU
						avg_fitting_results[cpx][cpy][3] += fitting_results[px][py][3]*fitting_results[px][py][3]; // SIGMA; NOTE: average of sigma is sqrt sum of variance (sigma squared)!
						
						avg_standard_error_fit_results[cpx][cpy][0] += standard_error_fit_results[px][py][0];
						avg_standard_error_fit_results[cpx][cpy][1] += standard_error_fit_results[px][py][1];
						avg_standard_error_fit_results[cpx][cpy][2] += standard_error_fit_results[px][py][2];
						avg_standard_error_fit_results[cpx][cpy][3] += standard_error_fit_results[px][py][3];
						
						avg_chi_squared_fit_results[cpx][cpy] += chi_squared_fit_results[px][py];
						
						avg_r_squared_fit_results[cpx][cpy] += r_squared_fit_results[px][py]; // RSLV: average R sqaured is just linear sum?
					}
				}
			}
		}
		
		// create average projection maps
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				int hc = hit_count_map_ip.get(px, py);
				if(hc > 0) // just to make sure
				{
					avg_results_step_3[px][py][0] /= hc;
					avg_results_step_3[px][py][1] /= hc; // RSLV: average n_x
					avg_results_step_3[px][py][2] /= hc; // RSLV: average n_y
					avg_results_step_3[px][py][3] /= hc;
					avg_results_step_3[px][py][4] /= hc; // RSLV: average s_x
					avg_results_step_3[px][py][5] /= hc; // RSLV: average s_y
					avg_results_step_3[px][py][6] /= hc;
					avg_results_step_3[px][py][7] /= hc;
					
					avg_frangi_measures[px][py][0] /= hc;
					avg_frangi_measures[px][py][1] /= hc;
					avg_frangi_measures[px][py][2] /= hc;
					avg_frangi_measures[px][py][3] /= hc;
					avg_frangi_measures[px][py][4] /= hc;
					
					avg_fitting_results[px][py][0] /= hc; // BG
					avg_fitting_results[px][py][1] /= hc; // AMP
					avg_fitting_results[px][py][2] /= hc; // MU
					avg_fitting_results[px][py][3] /= hc; // SIGMA; NOTE: average of sigma is sqrt sum of variance (sigma squared)!
					avg_fitting_results[px][py][3] = Math.sqrt(avg_fitting_results[px][py][3]);
					
					avg_standard_error_fit_results[px][py][0] /= hc;
					avg_standard_error_fit_results[px][py][1] /= hc;
					avg_standard_error_fit_results[px][py][2] /= hc;
					avg_standard_error_fit_results[px][py][3] /= hc;
					
					avg_chi_squared_fit_results[px][py] /= hc;
					
					avg_r_squared_fit_results[px][py] /= hc; // RSLV: average R sqaured is just linear sum?
				}
			}
		}
		Profiling.toc("Step 5a: Filter line points and creating projection maps");
		
		// ---------------------------------------------------------------------
		
		// Step 5b: additional filtering of the valid line points
		
		// RSLV: filter small particles using opening (erode+dilate) and closing (dilate+erode)
		ImageProcessor hit_map_filtered_ip = hit_map_ip.duplicate();
		((ByteProcessor)hit_map_filtered_ip).erode(3, 0);
		((ByteProcessor)hit_map_filtered_ip).dilate(3, 0);
		((ByteProcessor)hit_map_filtered_ip).dilate(3, 0);
		((ByteProcessor)hit_map_filtered_ip).erode(3, 0);
		
		// RSLV: filter using minimum connected components
		// NOTE: because of recursion it can run out of stack memory with large segments!
//		ImageProcessor hit_map_filtered_ip = hit_map_ip.duplicate();
//		ConnectedComponents.run(hit_map_filtered_ip, ConnectedComponents.Connectivity.EIGHT_CONNECTIVITY, 10, 1000000);
		
		// TODO: seperate creation of projection maps?
		
		// TODO: DEBUG: show intermediate images
		ImageProcessor r_squared_avg_projection_map_ip = new FloatProcessor(image_width, image_height);
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				r_squared_avg_projection_map_ip.setf(px, py, (float)avg_r_squared_fit_results[px][py]);
			}
		}
		
		// DEBUG: show intermediate images
		if(DEBUG_MODE_ENABLED)
		{
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
		
		//if(true)
		//	return null;
		
		// *********************************************************************
		
		// construct ODF map (through scale-space search)
		Profiling.tic();
		double[][][] odf_map = new double[image_width][image_height][ODF_STEPS];
		double[][][] odf_min_max_map = new double[image_width][image_height][2]; // [0] = min, [1] = max
		double[][][] odf_map2 = new double[image_width][image_height][ODF_STEPS];
		double[][][] odf_min_max_map2 = new double[image_width][image_height][2]; // [0] = min, [1] = max
		
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				odf_min_max_map[px][py][0] = Double.MAX_VALUE;
				odf_min_max_map[px][py][1] = Double.MIN_VALUE;
				odf_min_max_map2[px][py][0] = Double.MAX_VALUE;
				odf_min_max_map2[px][py][1] = Double.MIN_VALUE;
			}
		}
		
		// find response for all orientation
		ImagePlus odf_maps_imp = IJ.createImage("ODF maps n(t)", image_width, image_height, ODF_STEPS, 32);
		ImagePlus odf_maps2_imp = IJ.createImage("ODF maps s(t)", image_width, image_height, ODF_STEPS, 32);
		
		for(int orientation_step = 0; orientation_step < ODF_STEPS; ++orientation_step)
		{
			// calculate derivatives of Gaussian from image in direction
			float current_orientation = (float)(orientation_step * ANGULAR_RESOLUTION);
			IJ.showStatus("Calculating derivative of Gaussian at angle=" + current_orientation + " anisotropy=" + ANISOTROPY_FACTOR);
			current_orientation = (float)Math.toRadians(current_orientation);
			
			// NOTE: rotate in reverse direction
			ImageProcessor odf_dx = DerivativeOfGaussian.derivativeX(ip_step_2, 0.5*SIGMA, ANISOTROPY_FACTOR*SIGMA, -current_orientation);
			ImageProcessor odf_dy = DerivativeOfGaussian.derivativeY(ip_step_2, 0.5*SIGMA, ANISOTROPY_FACTOR*SIGMA, -current_orientation);
			ImageProcessor odf_dxdx = DerivativeOfGaussian.derivativeXX(ip_step_2, 0.5*SIGMA, ANISOTROPY_FACTOR*SIGMA, -current_orientation);
			ImageProcessor odf_dxdy = DerivativeOfGaussian.derivativeXY(ip_step_2, 0.5*SIGMA, ANISOTROPY_FACTOR*SIGMA, -current_orientation);
			ImageProcessor odf_dydx = DerivativeOfGaussian.derivativeYX(ip_step_2, 0.5*SIGMA, ANISOTROPY_FACTOR*SIGMA, -current_orientation);
			ImageProcessor odf_dydy = DerivativeOfGaussian.derivativeYY(ip_step_2, 0.5*SIGMA, ANISOTROPY_FACTOR*SIGMA, -current_orientation);
			
			ImageProcessor odf_map_ip = new FloatProcessor(image_width, image_height);
			odf_maps_imp.getImageStack().setProcessor(odf_map_ip, orientation_step+1); // NOTE: 1-offset index
			ImageProcessor odf_map2_ip = new FloatProcessor(image_width, image_height);
			odf_maps2_imp.getImageStack().setProcessor(odf_map2_ip, orientation_step+1); // NOTE: 1-offset index
			
			for(int py = 0; py < image_height; ++py)
			{
				for(int px = 0; px < image_width; ++px)
				{
					// skip if point has already been filtered out by the previous step (to save computation time); NOTE: does NOT save computation time
//					if(hit_map_filtered_ip.get(px, py) == 0)
//					{
//						odf_map[px][py][orientation_step] = 10.0;
//						continue;
//					}
					
					// calculate Hessian matrix
					Matrix m = new Matrix(2, 2, 0); // 2x2 RC matrix with zeros
					m.set(0, 0, odf_dxdx.getf(px, py));
					m.set(0, 1, odf_dxdy.getf(px, py));
					m.set(1, 0, odf_dydx.getf(px, py));
					m.set(1, 1, odf_dydy.getf(px, py));
					
					// compute eigenvalues and eigenvectors
					EigenvalueDecomposition evd = m.eig();
					Matrix d = evd.getD();
					Matrix v = evd.getV();
					
					// determine first and second eigenvalue and eigenvector
					double first_eigenvalue = 0.0; // |n(t)|
					double first_eigenvector_x = 0.0; // n(t) -> perpendicular to s(t)
					double first_eigenvector_y = 0.0; // n(t) -> perpendicular to s(t)
					double second_eigenvalue = 0.0;
					double second_eigenvector_x = 0.0;
					double second_eigenvector_y = 0.0;
					
					if(d.get(0,0) <= d.get(1,1))
					{
						// d(0,0) is most negative minimum eigenvalue
						first_eigenvalue = d.get(0,0); // L1
						first_eigenvector_x = v.get(0,0); // V1x
						first_eigenvector_y = v.get(1,0); // V1y
						second_eigenvalue = d.get(1,1); // L2
						second_eigenvector_x = v.get(0,1); // V2x
						second_eigenvector_y = v.get(1,1); // V2y
					}
					else
					{
						// d(1,1) is most negative minimum eigenvalue
						first_eigenvalue = d.get(1,1); // L1
						first_eigenvector_x = v.get(0,1); // V1x
						first_eigenvector_y = v.get(1,1); // V1y
						second_eigenvalue = d.get(0,0); // L2
						second_eigenvector_x = v.get(0,0); // V2x
						second_eigenvector_y = v.get(1,0); // V2y
					}
					
					// RSLV: store eigenvalues and eigenvector for all orientations?
					//odf_map[px][py][orientation_step][0] = first_eigenvalue;
					//odf_map[px][py][orientation_step][1] = first_eigenvector_x;
					//odf_map[px][py][orientation_step][2] = first_eigenvector_y;
					//odf_map[px][py][orientation_step][3] = second_eigenvalue;
					//odf_map[px][py][orientation_step][4] = second_eigenvector_x;
					//odf_map[px][py][orientation_step][5] = second_eigenvector_y;
					
					// RSLV: or store only response (looking for normal -> min value)
					odf_map[px][py][orientation_step] = first_eigenvalue;
					odf_map_ip.setf(px, py, (float)first_eigenvalue);
					if(first_eigenvalue < odf_min_max_map[px][py][0])
					{
						odf_min_max_map[px][py][0] = first_eigenvalue;
					}
					if(first_eigenvalue > odf_min_max_map[px][py][1])
					{
						odf_min_max_map[px][py][1] = first_eigenvalue;
					}
					
					// RSLV: or store only response (looking for s(t) -> max value?)
					odf_map2[px][py][orientation_step] = second_eigenvalue;
					odf_map2_ip.setf(px, py, (float)second_eigenvalue);
					if(second_eigenvalue < odf_min_max_map2[px][py][0])
					{
						odf_min_max_map2[px][py][0] = second_eigenvalue;
					}
					if(second_eigenvalue > odf_min_max_map2[px][py][1])
					{
						odf_min_max_map2[px][py][1] = second_eigenvalue;
					}
				}
			}
		}
		Profiling.toc("Step 3a: calculating ODF map");
		
		// show intermediate results
		if(DEBUG_MODE_ENABLED)
		{
			odf_maps_imp.resetDisplayRange();
			odf_maps_imp.show();
			odf_maps2_imp.resetDisplayRange();
			odf_maps2_imp.show();
		}
		
		// STEP 3b: finding the two peak orientations
		Profiling.tic();
		
		// estimate background of side lobes for peak extraction
		double[][][] odf_background = new double[image_width][image_height][ODF_STEPS];
		double[][][] odf_background2 = new double[image_width][image_height][ODF_STEPS];
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				for(int i = 0; i < ODF_STEPS; ++i)
				{
					// find minimum for opening, maximum for closing
					double min_val = Double.MAX_VALUE;
					double max_val2 = Double.MIN_VALUE;
					for(int k = (int)(-(0.5*ODF_PEAK_WIDTH)); k <= (int)(0.5*ODF_PEAK_WIDTH); ++k)
					{
						// RSLV: use wrap around?
						int ik = (i + k + ODF_STEPS) % ODF_STEPS;
						//if(i + k >= 0 && i + k < ODF_STEPS)
						//{
							if(odf_map[px][py][ik] < min_val)
							{
								min_val = odf_map[px][py][ik];
							}
							if(odf_map2[px][py][ik] > max_val2)
							{
								max_val2 = odf_map2[px][py][ik];
							}
						//}
					}
					odf_background[px][py][i] = min_val;
					odf_background2[px][py][i] = max_val2;
				}
			}
		}
		
		// remove side lobe background and TODO: reset min/max values
		double[][][] odf_filtered_map = new double[image_width][image_height][ODF_STEPS];
		double[][][] odf_filtered_map2 = new double[image_width][image_height][ODF_STEPS];
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				// find first peak
				for(int i = 0; i < ODF_STEPS; ++i)
				{
					odf_filtered_map[px][py][i] = odf_map[px][py][i] - odf_background[px][py][i];
					odf_filtered_map2[px][py][i] = odf_map2[px][py][i] - odf_background2[px][py][i];
				}
			}
		}
		
		// find peaks
		int[][][] odf_peaks = new int[image_width][image_height][2];
		int[][][] odf_peaks2 = new int[image_width][image_height][2];
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				// find first peak
				double first_peak_value = Double.MIN_VALUE;
				double first_peak_value2 = Double.MAX_VALUE;
				int first_peak_position = -1;
				int first_peak_position2 = -1;
				for(int i = 0; i < ODF_STEPS; ++i)
				{
					if(odf_filtered_map[px][py][i] > first_peak_value)
					{
						first_peak_value = odf_filtered_map[px][py][i];
						first_peak_position = i;
					}
					
					if(odf_filtered_map2[px][py][i] < first_peak_value2)
					{
						first_peak_value2 = odf_filtered_map2[px][py][i];
						first_peak_position2 = i;
					}
				}
				odf_peaks[px][py][0] = first_peak_position;
				odf_peaks2[px][py][0] = first_peak_position2;
				
				// calculate range that should be disabled for second peak
				// NOTE: upper and lower bound may wrap around; i.e. upper > lower!
				int lower_bound = (odf_peaks[px][py][0]+ODF_STEPS-(int)(0.5*ODF_PEAK_WIDTH))%ODF_STEPS;
				int upper_bound = (odf_peaks[px][py][0]+ODF_STEPS+(int)(0.5*ODF_PEAK_WIDTH))%ODF_STEPS;
				
				int lower_bound2 = (odf_peaks2[px][py][0]+ODF_STEPS-(int)(0.5*ODF_PEAK_WIDTH))%ODF_STEPS;
				int upper_bound2 = (odf_peaks2[px][py][0]+ODF_STEPS+(int)(0.5*ODF_PEAK_WIDTH))%ODF_STEPS;
				
				// find second peak
				double second_peak_value = Double.MIN_VALUE;
				double second_peak_value2 = Double.MAX_VALUE;
				int second_peak_position = -1;
				int second_peak_position2 = -1;
				for(int i = 0; i < ODF_STEPS; ++i)
				{
					// select only second peak away from first peak
					if((lower_bound < upper_bound && (i <= lower_bound || i >= upper_bound)) || (upper_bound < lower_bound && (i <= lower_bound && i >= upper_bound)))
					{
						if(odf_filtered_map[px][py][i] > second_peak_value)
						{
							second_peak_value = odf_filtered_map[px][py][i];
							second_peak_position = i;
						}
					}
					
					// select second peak away from first peak
					if((lower_bound2 < upper_bound2 && (i <= lower_bound2 || i >= upper_bound2)) || (upper_bound2 < lower_bound2 && (i <= lower_bound2 && i >= upper_bound2)))
					{
						if(odf_filtered_map2[px][py][i] < second_peak_value2)
						{
							second_peak_value2 = odf_filtered_map2[px][py][i];
							second_peak_position2 = i;
						}
					}
				}
				odf_peaks[px][py][1] = second_peak_position;
				odf_peaks2[px][py][1] = second_peak_position2;
			}
		}
		Profiling.toc("Step 3a: calculating ODF map");
		
		// *********************************************************************
		
		if(SHOW_ODF_OVERLAY)
		{
			// duplicate image for vector overlay
			Profiling.tic();
			ip_original.setInterpolationMethod(ImageProcessor.NONE); // NONE, NEAREST_NEIGHBOR, BILINEAR, BICUBIC
			ImageProcessor odf_overlay_ip = ip_original.duplicate();
			Overlay odf_overlay = new Overlay();
			
			// create vector overlay of primary [and secondary] eigenvector on top of *original* image
			for(int py = 0; py < image_height; ++py)
			{
				for(int px = 0; px < image_width; ++px)
				{
					// check filter status of pixel
					boolean filtered = (hit_map_ip.get(px, py) == 0); // RSLV: hit_count_map_ip.get(px, py) >= LINE_THREHSOLD?
					
					// DEBUG: overlay vector on duplicate of image
					double cx = px+0.5;
					double cy = py+0.5;
					
					// draw first peak
					double fpx = 0.4*Math.cos(Math.toRadians(odf_peaks[px][py][0]*ANGULAR_RESOLUTION));
					double fpy = 0.4*Math.sin(Math.toRadians(odf_peaks[px][py][0]*ANGULAR_RESOLUTION));
					
					Roi first_odf_peak_roi = null;
					if(USE_VECTOR_NORMALS)
					{
						// keep using normals
						first_odf_peak_roi = new ij.gui.Line(cx-fpx, cy-fpy, cx+fpx, cy+fpy);
					}
					else
					{
						// flip normals to vectors
						first_odf_peak_roi = new ij.gui.Line(cx-fpy, cy+fpx, cx+fpy, cy-fpx);
					}
					
					first_odf_peak_roi.setStrokeWidth(0.0);
					first_odf_peak_roi.setStrokeColor(filtered ? Color.DARK_GRAY : new Color((int)(255*Math.abs(Math.cos(Math.toRadians(odf_peaks[px][py][0]*ANGULAR_RESOLUTION+(USE_VECTOR_NORMALS ? 0 : 90))))), (int)(255*Math.abs(Math.sin(Math.toRadians(odf_peaks[px][py][0]*ANGULAR_RESOLUTION+(USE_VECTOR_NORMALS ? 0 : 90))))), 0, (USE_TRANSPARENCY ? (int)(255*Math.max(0, r_squared_fit_results[px][py])) : 255)));
					first_odf_peak_roi.setPosition(0); // NOTE: redundant in single frame instance
					odf_overlay.add(first_odf_peak_roi);
					
					// draw second peak
					double spx = 0.4*Math.cos(Math.toRadians(odf_peaks[px][py][1]*ANGULAR_RESOLUTION));
					double spy = 0.4*Math.sin(Math.toRadians(odf_peaks[px][py][1]*ANGULAR_RESOLUTION));
					
					Roi second_odf_peak_roi = null;
					if(USE_VECTOR_NORMALS)
					{
						// keep using normals
						second_odf_peak_roi = new ij.gui.Line(cx-spx, cy-spy, cx+spx, cy+spy);
					}
					else
					{
						// flip normals to vectors
						second_odf_peak_roi = new ij.gui.Line(cx-spy, cy+spx, cx+spy, cy-spx);
					}
					
					second_odf_peak_roi.setStrokeWidth(0.0);
					second_odf_peak_roi.setStrokeColor(filtered ? Color.DARK_GRAY : new Color((int)(255*Math.abs(Math.cos(Math.toRadians(odf_peaks[px][py][1]*ANGULAR_RESOLUTION+(USE_VECTOR_NORMALS ? 0 : 90))))), (int)(255*Math.abs(Math.sin(Math.toRadians(odf_peaks[px][py][1]*ANGULAR_RESOLUTION+(USE_VECTOR_NORMALS ? 0 : 90))))), 0, (USE_TRANSPARENCY ? (int)(255*Math.max(0, r_squared_fit_results[px][py])) : 255)));
					second_odf_peak_roi.setPosition(0); // NOTE: redundant in single frame instance
					odf_overlay.add(second_odf_peak_roi);
					
/*					// also add vector normal to overlay
					double npx = 0.4*results_step_3[px][py][1];
					double npy = 0.4*results_step_3[px][py][2];
					
					Roi vector_normal_roi = null;
					if(USE_VECTOR_NORMALS)
					{
						// keep using normals
						vector_normal_roi = new ij.gui.Line(cx-npx, cy-npy, cx+npx, cy+npy);
					}
					else
					{
						// flip normals to vectors
						vector_normal_roi = new ij.gui.Line(cx-npy, cy+npx, cx+npy, cy-npx);
					}
					
					vector_normal_roi.setStrokeWidth(0.0);
					vector_normal_roi.setStrokeColor(Color.GRAY);
					vector_normal_roi.setPosition(0); // NOTE: redundant in single frame instance
					odf_overlay.add(vector_normal_roi);
*/				}
			}
			
			// show scaled image with ODF overlay
			final ImagePlus odf_overlay_imp = new ImagePlus("DEBUG: ODF overlay", odf_overlay_ip); // TMP: pass as model to custom MouseAdapter
			//odf_overlay_imp.resetDisplayRange();
			odf_overlay_imp.setOverlay(odf_overlay);
			//odf_overlay_imp.updateAndRepaintWindow();
			odf_overlay_imp.show();
			
			// -----------------------------------------------------------------
			
			// RSLV: if this interaction tool will be used in the release; create subclass of MouseAdapter to accept the model as parameters, rather than using final class from parent class!
			final ImageWindow img_wnd = odf_overlay_imp.getWindow(); // final for inner class access
			final ImageCanvas img_cnv = img_wnd.getCanvas(); // final for inner class access
			final double[][][] odf_map_tmp = odf_map; // TMP: pass as model to custom MouseAdapter
			final double[][][] odf_filtered_map_tmp = odf_filtered_map; // TMP: pass as model to custom MouseAdapter
			final double[][][] odf_min_max_map_tmp = odf_min_max_map; // TMP: pass as model to custom MouseAdapter
			final int[][][] odf_peaks_tmp = odf_peaks; // TMP: pass as model to custom MouseAdapter
			final double[][][] odf_map2_tmp = odf_map2; // TMP: pass as model to custom MouseAdapter
			final double[][][] odf_filtered_map2_tmp = odf_filtered_map2; // TMP: pass as model to custom MouseAdapter
			final double[][][] odf_min_max_map2_tmp = odf_min_max_map2; // TMP: pass as model to custom MouseAdapter
			final int[][][] odf_peaks2_tmp = odf_peaks2; // TMP: pass as model to custom MouseAdapter
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
				private OvalRoi odf_range_one = null;
				private OvalRoi odf_range_three = null;
				
				private ij.gui.Line odf_first_peak = null;
				private ij.gui.Line odf_second_peak = null;
				
				private PolygonRoi trace = null;
				
				private PlotWindow classical_plot_wnd = null;
				private PlotWindow odf_map_plot_wnd = null;
				private Plot odf_map_plot = null;
				
				private PlotWindow classical_plot2_wnd = null;
				private PlotWindow odf_map_plot2_wnd = null;
				private Plot odf_map_plot2 = null;
				
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
						odf_overlay_imp.getOverlay().remove(clicked_pixel);
					}
					clicked_pixel = new Roi(current_x, current_y, 1, 1);
					clicked_pixel.setStrokeColor(Color.CYAN);
					clicked_pixel.setStrokeWidth(0.0);
					odf_overlay_imp.getOverlay().add(clicked_pixel);
					
					// draw ODF ranges
					if(odf_range_one != null)
					{
						// remove previous ROI first
						odf_overlay_imp.getOverlay().remove(odf_range_one);
					}
					odf_range_one = new OvalRoi(current_x - 0.5*ANISOTROPY_FACTOR*SIGMA + 0.5, current_y - 0.5*ANISOTROPY_FACTOR*SIGMA + 0.5, ANISOTROPY_FACTOR*SIGMA, ANISOTROPY_FACTOR*SIGMA);
					odf_range_one.setStrokeColor(Color.BLUE);
					odf_range_one.setStrokeWidth(0.0);
					odf_overlay_imp.getOverlay().add(odf_range_one);
					
					if(odf_range_three != null)
					{
						// remove previous ROI first
						odf_overlay_imp.getOverlay().remove(odf_range_three);
					}
					odf_range_three = new OvalRoi(current_x - 3*0.5*ANISOTROPY_FACTOR*SIGMA + 0.5, current_y - 3*0.5*ANISOTROPY_FACTOR*SIGMA + 0.5, 3*ANISOTROPY_FACTOR*SIGMA, 3*ANISOTROPY_FACTOR*SIGMA);
					odf_range_three.setStrokeColor(Color.BLUE);
					odf_range_three.setStrokeWidth(0.0);
					odf_overlay_imp.getOverlay().add(odf_range_three);
					
					// TODO: draw enlarged vectors / normals of selected pixel
					// DEBUG: overlay vector on duplicate of image
					double cx = current_x+0.5;
					double cy = current_y+0.5;
					
					// draw first peak
					double fpx = 0.4*ANISOTROPY_FACTOR*Math.cos(Math.toRadians(odf_peaks_tmp[current_x][current_y][0]*ANGULAR_RESOLUTION));
					double fpy = 0.4*ANISOTROPY_FACTOR*Math.sin(Math.toRadians(odf_peaks_tmp[current_x][current_y][0]*ANGULAR_RESOLUTION));
					
					// first remove existing roi
					if(odf_first_peak != null)
					{
						// remove previous ROI first
						odf_overlay_imp.getOverlay().remove(odf_first_peak);
					}
					
					if(USE_VECTOR_NORMALS)
					{
						// keep using normals
						odf_first_peak = new ij.gui.Line(cx-fpx, cy-fpy, cx+fpx, cy+fpy);
					}
					else
					{
						// flip normals to vectors
						odf_first_peak = new ij.gui.Line(cx-fpy, cy+fpx, cx+fpy, cy-fpx);
					}
					
					odf_first_peak.setStrokeWidth(0.5);
					odf_first_peak.setStrokeColor(new Color((int)(255*Math.abs(Math.cos(Math.toRadians(odf_peaks_tmp[current_x][current_y][0]*ANGULAR_RESOLUTION+(USE_VECTOR_NORMALS ? 0 : 90))))), (int)(255*Math.abs(Math.sin(Math.toRadians(odf_peaks_tmp[current_x][current_y][0]*ANGULAR_RESOLUTION+(USE_VECTOR_NORMALS ? 0 : 90))))), 0, 255));
					odf_first_peak.setPosition(0); // NOTE: redundant in single frame instance
					odf_overlay_imp.getOverlay().add(odf_first_peak);
					
					// draw second peak
					double spx = 0.4*ANISOTROPY_FACTOR*Math.cos(Math.toRadians(odf_peaks_tmp[current_x][current_y][1]*ANGULAR_RESOLUTION));
					double spy = 0.4*ANISOTROPY_FACTOR*Math.sin(Math.toRadians(odf_peaks_tmp[current_x][current_y][1]*ANGULAR_RESOLUTION));
					
					// first remove existing roi
					if(odf_second_peak != null)
					{
						// remove previous ROI first
						odf_overlay_imp.getOverlay().remove(odf_second_peak);
					}
					
					if(USE_VECTOR_NORMALS)
					{
						// keep using normals
						odf_second_peak = new ij.gui.Line(cx-spx, cy-spy, cx+spx, cy+spy);
					}
					else
					{
						// flip normals to vectors
						odf_second_peak = new ij.gui.Line(cx-spy, cy+spx, cx+spy, cy-spx);
					}
					
					odf_second_peak.setStrokeWidth(0.5);
					odf_second_peak.setStrokeColor(new Color((int)(255*Math.abs(Math.cos(Math.toRadians(odf_peaks_tmp[current_x][current_y][1]*ANGULAR_RESOLUTION+(USE_VECTOR_NORMALS ? 0 : 90))))), (int)(255*Math.abs(Math.sin(Math.toRadians(odf_peaks_tmp[current_x][current_y][1]*ANGULAR_RESOLUTION+(USE_VECTOR_NORMALS ? 0 : 90))))), 0, 255));
					odf_second_peak.setPosition(0); // NOTE: redundant in single frame instance
					odf_overlay_imp.getOverlay().add(odf_second_peak);
					
					//odf_overlay_imp.updateAndRepaintWindow(); // force update of window
					
					// determine action:
					// [1] new pixel selected -> show information
					// [2] same pixel selected -> show trace menu
					if(false) // TODO: implement tracing for ODF map!
					//if(previous_x == current_x && previous_y == current_y)
					{
						// show trace line options dialog
						GenericDialog gd = new GenericDialog("Trace line from pixel x=" + current_x + ", y=" + current_y);
						gd.addNumericField("Max steps", 100, 0);
						gd.addNumericField("Step size", 1, 2);
						//gd.hideCancelButton();
						gd.showDialog();
						if(gd.wasOKed()) // !gd.wasCanceled()
						{
							// trace line from selected pixel
							final int MAX_STEPS = (int)gd.getNextNumber();
							final double STEP_SIZE = gd.getNextNumber();
							
							// TODO: implement tracing
							System.err.println("TODO: implement tracing for ODF map!");
						}
					}
					else
					{
						// create classical plot
						double[] classical_filtered_xs = new double[ODF_STEPS];
						for(int i = 0; i < ODF_STEPS; ++i)
						{
							classical_filtered_xs[i] = i * ANGULAR_RESOLUTION;
						}
						
						Plot classical_plot = new Plot("Classical plot @ (x=" + current_x + ", y=" + current_y + ")", "Orientation", "Magnitude", classical_filtered_xs, odf_filtered_map_tmp[current_x][current_y]);
						classical_plot.setSize(300, 300);
						classical_plot.setFrameSize(300, 300);
						
						// marker for first peak
						if(odf_peaks_tmp[current_x][current_y][0] != -1)
						{
							classical_plot.setColor(Color.RED);
							classical_plot.setLineWidth(2);
							classical_plot.addPoints(new double[]{odf_peaks_tmp[current_x][current_y][0]*ANGULAR_RESOLUTION}, new double[]{odf_filtered_map_tmp[current_x][current_y][odf_peaks_tmp[current_x][current_y][0]]}, Plot.CIRCLE);
						}
						
						// marker for second peak
						if(odf_peaks_tmp[current_x][current_y][1] != -1)
						{
							classical_plot.setColor(Color.GREEN);
							classical_plot.setLineWidth(2);
							classical_plot.addPoints(new double[]{odf_peaks_tmp[current_x][current_y][1]*ANGULAR_RESOLUTION}, new double[]{odf_filtered_map_tmp[current_x][current_y][odf_peaks_tmp[current_x][current_y][1]]}, Plot.CIRCLE);
						}
						
						// set line plot style
						classical_plot.setColor(Color.DARK_GRAY);
						classical_plot.setLineWidth(1);
						
						if(classical_plot_wnd != null && !classical_plot_wnd.isClosed())
						{
							classical_plot_wnd.drawPlot(classical_plot);
							classical_plot_wnd.setTitle("Classical plot @ (x=" + current_x + ", y=" + current_y + ")");
							// RSLV: set focus on plot window?
						}
						else
						{
							classical_plot_wnd = classical_plot.show();
							classical_plot_wnd.setLocation(660, 30);
						}
						
						// plot ODF as polar plot
						double[] odf_filtered_xs = new double[2*ODF_STEPS+1];
						double[] odf_filtered_ys = new double[2*ODF_STEPS+1];
						double odf_filtered_x_range = 0;
						double odf_filtered_y_range = 0;
						for(int i = 0; i < ODF_STEPS; ++i)
						{
							double vx = Math.cos(Math.toRadians(i*ANGULAR_RESOLUTION));
							double vy = Math.sin(Math.toRadians(i*ANGULAR_RESOLUTION));
							double px = vx * (200 + odf_filtered_map_tmp[current_x][current_y][i]);
							double py = vy * (200 + odf_filtered_map_tmp[current_x][current_y][i]);
							odf_filtered_xs[i] = px;
							odf_filtered_ys[i] = -py; // NOTE: flip y-axis
							odf_filtered_xs[ODF_STEPS+i] = -odf_filtered_xs[i];
							odf_filtered_ys[ODF_STEPS+i] = -odf_filtered_ys[i];
							if(px > odf_filtered_x_range) odf_filtered_x_range = px;
							if(py > odf_filtered_y_range) odf_filtered_y_range = py;
						}
						odf_filtered_xs[2*ODF_STEPS] = odf_filtered_xs[0]; // closed polyline
						odf_filtered_ys[2*ODF_STEPS] = odf_filtered_ys[0]; // closed polyline
						double odf_filtered_max_range = Math.max(odf_filtered_x_range, odf_filtered_y_range);
						
						odf_map_plot = new Plot("ODF plot @ (x=" + current_x + ", y=" + current_y + ")", "Polar plot", "Polar plot", odf_filtered_xs, odf_filtered_ys);
						odf_map_plot.setLimits(-1.05*odf_filtered_max_range, 1.05*odf_filtered_max_range, -1.05*odf_filtered_max_range, 1.05*odf_filtered_max_range);
						odf_map_plot.setSize(300, 300);
						odf_map_plot.setFrameSize(300, 300);
						
						// add current vector normal direction
						odf_map_plot.setColor(Color.YELLOW);
						odf_map_plot.setLineWidth(2);
						odf_map_plot.drawLine(hessian_results_tmp[current_x][current_y][1]*0.5*odf_filtered_max_range, -hessian_results_tmp[current_x][current_y][2]*0.5*odf_filtered_max_range, -hessian_results_tmp[current_x][current_y][1]*0.5*odf_filtered_max_range, hessian_results_tmp[current_x][current_y][2]*0.5*odf_filtered_max_range); // NOTE: flip y-axis
						
						// draw first peak
						odf_map_plot.setColor(Color.RED);
						odf_map_plot.setLineWidth(2);
						if(odf_peaks_tmp[current_x][current_y][0] != -1)
						{
							odf_map_plot.drawLine(Math.cos(Math.toRadians(odf_peaks_tmp[current_x][current_y][0]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map_tmp[current_x][current_y][odf_peaks_tmp[current_x][current_y][0]], -Math.sin(Math.toRadians(odf_peaks_tmp[current_x][current_y][0]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map_tmp[current_x][current_y][odf_peaks_tmp[current_x][current_y][0]], -Math.cos(Math.toRadians(odf_peaks_tmp[current_x][current_y][0]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map_tmp[current_x][current_y][odf_peaks_tmp[current_x][current_y][0]], Math.sin(Math.toRadians(odf_peaks_tmp[current_x][current_y][0]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map_tmp[current_x][current_y][odf_peaks_tmp[current_x][current_y][0]]); // NOTE: flip y-axis
						}
						
						// draw second peak
						odf_map_plot.setColor(Color.GREEN);
						odf_map_plot.setLineWidth(2);
						if(odf_peaks_tmp[current_x][current_y][1] != -1)
						{
							odf_map_plot.drawLine(Math.cos(Math.toRadians(odf_peaks_tmp[current_x][current_y][1]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map_tmp[current_x][current_y][odf_peaks_tmp[current_x][current_y][1]], -Math.sin(Math.toRadians(odf_peaks_tmp[current_x][current_y][1]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map_tmp[current_x][current_y][odf_peaks_tmp[current_x][current_y][1]], -Math.cos(Math.toRadians(odf_peaks_tmp[current_x][current_y][1]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map_tmp[current_x][current_y][odf_peaks_tmp[current_x][current_y][1]], Math.sin(Math.toRadians(odf_peaks_tmp[current_x][current_y][1]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map_tmp[current_x][current_y][odf_peaks_tmp[current_x][current_y][1]]); // NOTE: flip y-axis
						}
						
						// set color and line width for ODF contour
						odf_map_plot.setColor(Color.DARK_GRAY);
						odf_map_plot.setLineWidth(1);
						
						// create new plot window or draw in existing plot window
						if(odf_map_plot_wnd != null && !odf_map_plot_wnd.isClosed())
						{
							odf_map_plot_wnd.drawPlot(odf_map_plot);
							odf_map_plot_wnd.setTitle("ODF @ (x=" + current_x + ", y=" + current_y + ")");
							// RSLV: set focus on plot window?
						}
						else
						{
							odf_map_plot_wnd = odf_map_plot.show();
							odf_map_plot_wnd.setLocation(660, 450);
						}
						
						// -----------------------------------------------------
						
						// create classical plot
						double[] classical_filtered_xs2 = new double[ODF_STEPS];
						for(int i = 0; i < ODF_STEPS; ++i)
						{
							classical_filtered_xs2[i] = i * ANGULAR_RESOLUTION;
						}
						
						Plot classical_plot2 = new Plot("Classical plot2 @ (x=" + current_x + ", y=" + current_y + ")", "Orientation", "Magnitude", classical_filtered_xs2, odf_filtered_map2_tmp[current_x][current_y]);
						classical_plot2.setSize(300, 300);
						classical_plot2.setFrameSize(300, 300);
						
						// marker for first peak
						if(odf_peaks2_tmp[current_x][current_y][0] != -1)
						{
							classical_plot2.setColor(Color.RED);
							classical_plot2.setLineWidth(2);
							classical_plot2.addPoints(new double[]{odf_peaks2_tmp[current_x][current_y][0]*ANGULAR_RESOLUTION}, new double[]{odf_filtered_map2_tmp[current_x][current_y][odf_peaks2_tmp[current_x][current_y][0]]}, Plot.CIRCLE);
						}
						
						// marker for second peak
						if(odf_peaks_tmp[current_x][current_y][1] != -1)
						{
							classical_plot2.setColor(Color.GREEN);
							classical_plot2.setLineWidth(2);
							classical_plot2.addPoints(new double[]{odf_peaks2_tmp[current_x][current_y][1]*ANGULAR_RESOLUTION}, new double[]{odf_filtered_map2_tmp[current_x][current_y][odf_peaks2_tmp[current_x][current_y][1]]}, Plot.CIRCLE);
						}
						
						// set line plot style
						classical_plot2.setColor(Color.DARK_GRAY);
						classical_plot2.setLineWidth(1);
						
						if(classical_plot2_wnd != null && !classical_plot2_wnd.isClosed())
						{
							classical_plot2_wnd.drawPlot(classical_plot2);
							classical_plot2_wnd.setTitle("Classical plot2 @ (x=" + current_x + ", y=" + current_y + ")");
							// RSLV: set focus on plot window?
						}
						else
						{
							classical_plot2_wnd = classical_plot2.show();
							classical_plot2_wnd.setLocation(1050, 30);
						}
						
						// plot ODF as polar plot
						double[] odf_filtered_xs2 = new double[2*ODF_STEPS+1];
						double[] odf_filtered_ys2 = new double[2*ODF_STEPS+1];
						double odf_filtered_x_range2 = 0;
						double odf_filtered_y_range2 = 0;
						for(int i = 0; i < ODF_STEPS; ++i)
						{
							double vx2 = Math.cos(Math.toRadians(i*ANGULAR_RESOLUTION));
							double vy2 = Math.sin(Math.toRadians(i*ANGULAR_RESOLUTION));
							double px2 = vx2 * (200 + Math.abs(odf_filtered_map2_tmp[current_x][current_y][i]));
							double py2 = vy2 * (200 + Math.abs(odf_filtered_map2_tmp[current_x][current_y][i]));
							odf_filtered_xs2[i] = px2;
							odf_filtered_ys2[i] = -py2; // NOTE: flip y-axis
							odf_filtered_xs2[ODF_STEPS+i] = -odf_filtered_xs2[i];
							odf_filtered_ys2[ODF_STEPS+i] = -odf_filtered_ys2[i];
							if(px2 > odf_filtered_x_range2) odf_filtered_x_range2 = px2;
							if(py2 > odf_filtered_y_range2) odf_filtered_y_range2 = py2;
						}
						odf_filtered_xs2[2*ODF_STEPS] = odf_filtered_xs2[0]; // closed polyline
						odf_filtered_ys2[2*ODF_STEPS] = odf_filtered_ys2[0]; // closed polyline
						double odf_filtered_max_range2 = Math.max(odf_filtered_x_range2, odf_filtered_y_range2);
						
						odf_map_plot2 = new Plot("ODF2 @ (x=" + current_x + ", y=" + current_y + ")", "Polar plot", "Polar plot", odf_filtered_xs2, odf_filtered_ys2);
						odf_map_plot2.setLimits(-1.05*odf_filtered_max_range2, 1.05*odf_filtered_max_range2, -1.05*odf_filtered_max_range2, 1.05*odf_filtered_max_range2);
						odf_map_plot2.setSize(300, 300);
						odf_map_plot2.setFrameSize(300, 300);
						
						// add current vector direction
						odf_map_plot2.setColor(Color.YELLOW);
						odf_map_plot2.setLineWidth(2);
						odf_map_plot2.drawLine(hessian_results_tmp[current_x][current_y][1]*0.5*odf_filtered_max_range2, -hessian_results_tmp[current_x][current_y][2]*0.5*odf_filtered_max_range2, -hessian_results_tmp[current_x][current_y][1]*0.5*odf_filtered_max_range2, hessian_results_tmp[current_x][current_y][2]*0.5*odf_filtered_max_range2); // NOTE: flip y-axis
						
						// draw first peak
						odf_map_plot2.setColor(Color.RED);
						odf_map_plot2.setLineWidth(2);
						if(odf_peaks2_tmp[current_x][current_y][0] != -1)
						{
							odf_map_plot2.drawLine(Math.cos(Math.toRadians(odf_peaks2_tmp[current_x][current_y][0]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map2_tmp[current_x][current_y][odf_peaks2_tmp[current_x][current_y][0]], -Math.sin(Math.toRadians(odf_peaks2_tmp[current_x][current_y][0]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map2_tmp[current_x][current_y][odf_peaks2_tmp[current_x][current_y][0]], -Math.cos(Math.toRadians(odf_peaks2_tmp[current_x][current_y][0]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map2_tmp[current_x][current_y][odf_peaks2_tmp[current_x][current_y][0]], Math.sin(Math.toRadians(odf_peaks2_tmp[current_x][current_y][0]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map2_tmp[current_x][current_y][odf_peaks2_tmp[current_x][current_y][0]]); // NOTE: flip y-axis
						}
						
						// draw second peak
						odf_map_plot2.setColor(Color.GREEN);
						odf_map_plot2.setLineWidth(2);
						if(odf_peaks2_tmp[current_x][current_y][1] != -1)
						{
							odf_map_plot2.drawLine(Math.cos(Math.toRadians(odf_peaks2_tmp[current_x][current_y][1]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map2_tmp[current_x][current_y][odf_peaks2_tmp[current_x][current_y][1]], -Math.sin(Math.toRadians(odf_peaks2_tmp[current_x][current_y][1]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map2_tmp[current_x][current_y][odf_peaks2_tmp[current_x][current_y][1]], -Math.cos(Math.toRadians(odf_peaks2_tmp[current_x][current_y][1]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map2_tmp[current_x][current_y][odf_peaks2_tmp[current_x][current_y][1]], Math.sin(Math.toRadians(odf_peaks2_tmp[current_x][current_y][1]*ANGULAR_RESOLUTION))*0.5*odf_filtered_map2_tmp[current_x][current_y][odf_peaks2_tmp[current_x][current_y][1]]); // NOTE: flip y-axis
						}
						
						// set color and line width for ODF contour
						odf_map_plot2.setColor(Color.DARK_GRAY);
						odf_map_plot2.setLineWidth(1);
						
						// create new plot window or draw in existing plot window
						if(odf_map_plot2_wnd != null && !odf_map_plot2_wnd.isClosed())
						{
							odf_map_plot2_wnd.drawPlot(odf_map_plot2);
							odf_map_plot2_wnd.setTitle("ODF2 @ (x=" + current_x + ", y=" + current_y + ")");
							// RSLV: set focus on plot window?
						}
						else
						{
							odf_map_plot2_wnd = odf_map_plot2.show();
							odf_map_plot2_wnd.setLocation(1050, 450);
						}
					}
					
					// force update on window
					odf_overlay_imp.updateAndRepaintWindow();
				}
			});
			
			Profiling.toc("DEBUG: drawing ODF overlay on duplicate of original image");
		}
		
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
					boolean filtered = (hit_map_ip.get(px, py) == 0); // RSLV: hit_count_map_ip.get(px, py) >= LINE_THREHSOLD?
					
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
			final double[][][] odf_map_tmp = odf_map; // TMP: pass as model to custom MouseAdapter
			final double[][][] odf_filtered_map_tmp = odf_filtered_map; // TMP: pass as model to custom MouseAdapter
			final double[][][] odf_min_max_map_tmp = odf_min_max_map; // TMP: pass as model to custom MouseAdapter
			final int[][][] odf_peaks_tmp = odf_peaks; // TMP: pass as model to custom MouseAdapter
			final double[][][] odf_map2_tmp = odf_map2; // TMP: pass as model to custom MouseAdapter
			final double[][][] odf_filtered_map2_tmp = odf_filtered_map2; // TMP: pass as model to custom MouseAdapter
			final double[][][] odf_min_max_map2_tmp = odf_min_max_map2; // TMP: pass as model to custom MouseAdapter
			final int[][][] odf_peaks2_tmp = odf_peaks2; // TMP: pass as model to custom MouseAdapter
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
				private OvalRoi odf_range_one = null;
				private OvalRoi odf_range_three = null;
				
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
					
					// draw ODF ranges
					if(odf_range_one != null)
					{
						// remove previous ROI first
						vector_overlay_imp.getOverlay().remove(odf_range_one);
					}
					odf_range_one = new OvalRoi(current_x - 0.5*ANISOTROPY_FACTOR*SIGMA + 0.5, current_y - 0.5*ANISOTROPY_FACTOR*SIGMA + 0.5, ANISOTROPY_FACTOR*SIGMA, ANISOTROPY_FACTOR*SIGMA);
					odf_range_one.setStrokeColor(Color.BLUE);
					odf_range_one.setStrokeWidth(0.0);
					vector_overlay_imp.getOverlay().add(odf_range_one);
					
					if(odf_range_three != null)
					{
						// remove previous ROI first
						vector_overlay_imp.getOverlay().remove(odf_range_three);
					}
					odf_range_three = new OvalRoi(current_x - 3*0.5*ANISOTROPY_FACTOR*SIGMA + 0.5, current_y - 3*0.5*ANISOTROPY_FACTOR*SIGMA + 0.5, 3*ANISOTROPY_FACTOR*SIGMA, 3*ANISOTROPY_FACTOR*SIGMA);
					odf_range_three.setStrokeColor(Color.BLUE);
					odf_range_three.setStrokeWidth(0.0);
					vector_overlay_imp.getOverlay().add(odf_range_three);
					
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
					raw_data_table.addValue("hit", hit_map_ip.get(px, py)); // hit map
					raw_data_table.addValue("hit_count", hit_count_map_ip.get(px, py)); // hit count map
				}
			}
			
			raw_data_table.show("Result of segmentation");
			Profiling.toc("Generating results table");
		}
		
		// *********************************************************************
		
		// Step Nx: prepare maps for tracing of individual lines
		
		// TEMP: linking datastructures between algorithms
		ImageProcessor valid_line_points_magnitude_ip = r_squared_avg_projection_map_ip;
		double[][][] line_points = results_step_3;
		
		//boolean[][] valid_line_points_mask = new boolean[image_width][image_height];
		ImageProcessor valid_line_points_mask_ip = new ByteProcessor(image_width, image_height);
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				// filter pixels with less than VALID_LINE_POINT_THRESHOLD hits
				if(hit_count_map_ip.get(px, py) >= VALID_LINE_POINT_THRESHOLD)
				{
					//valid_line_points_mask[px][py] = true;
					valid_line_points_mask_ip.set(px, py, 255);
				}
				else
				{
					//valid_line_points_mask[px][py] = false;
					valid_line_points_mask_ip.set(px, py, 0);
				}
			}
		}
		
		// skeletonize valid line point map
		ImageProcessor valid_line_points_skeleton_ip = valid_line_points_mask_ip.duplicate();
		
		if(SKELETONIZE_POINTS_MASK)
		{
			valid_line_points_skeleton_ip.invert();
			((ByteProcessor)valid_line_points_skeleton_ip).skeletonize();
			valid_line_points_skeleton_ip.invert();
		}
		
		// remove single islets
		ImageProcessor valid_line_points_filtered_ip = valid_line_points_skeleton_ip.duplicate();
		if(FILTER_POINTS_MASK)
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
								sum += valid_line_points_filtered_ip.getPixel(px+kx, py+ky); // NOTE: getPixel for bounds checking!
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
		
		/*
		ImageProcessor junctions_ip = new ByteProcessor(image_width, image_height);
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				int sum = 0;
				for(int ky = -1; ky <= 1; ++ky)
				{
					for(int kx = -1; kx <= 1; ++kx)
					{
						sum += junctions_ip.getPixel(px+kx, py+ky); // NOTE: getPixel for bounds checking!
					}
				}
				
				if(sum >= 4*255) // NOTE: assumes binary mask with [BG=0|FG=255]
				{
					junctions_ip.set(px, py, 255);
				}
			}
		}
		*/
		
		// create final valid line point mask for line tracing
		boolean[][] valid_line_points_mask = new boolean[image_width][image_height];
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
			{
				// filter pixels with less than VALID_LINE_POINT_THRESHOLD hits
				valid_line_points_mask[px][py] = valid_line_points_filtered_ip.get(px, py) > 0;
			}
		}
		
		// show intermediate images
		if(DEBUG_MODE_ENABLED)
		{
			ImagePlus valid_line_points_mask_imp = new ImagePlus("valid line points mask", valid_line_points_mask_ip);
			valid_line_points_mask_imp.resetDisplayRange();
			valid_line_points_mask_imp.show();
			
			if(SKELETONIZE_POINTS_MASK)
			{
				ImagePlus valid_line_points_skeleton_imp = new ImagePlus("valid line points skeleton", valid_line_points_skeleton_ip);
				valid_line_points_skeleton_imp.resetDisplayRange();
				valid_line_points_skeleton_imp.show();
			}
			
			if(FILTER_POINTS_MASK)
			{
				ImagePlus valid_line_points_filtered_imp = new ImagePlus("valid line points filtered", valid_line_points_filtered_ip);
				valid_line_points_filtered_imp.resetDisplayRange();
				valid_line_points_filtered_imp.show();
			}
			
//			ImagePlus junctions_imp = new ImagePlus("junctions", junctions_ip);
//			junctions_imp.resetDisplayRange();
//			junctions_imp.show();
			
		}
		
		// *********************************************************************
		
		// create ordered list of seeding points, selecting high R^2 and high count first: {px, py, R^2, hit count}
		Vector<Quadruple<Integer, Integer, Double, Integer> > seeding_points = new Vector<Quadruple<Integer, Integer, Double, Integer> >(); // RSLV: I know the number of valid line points in advance, set initial capacity?
		for(int px = 0; px < image_width; ++px)
		{
			for(int py = 0; py < image_height; ++py)
			{
				// use only valid line points
				if(valid_line_points_filtered_ip.get(px, py) != 0)
				{
					// create point vector
					Quadruple<Integer, Integer, Double, Integer> q = new Quadruple<Integer, Integer, Double, Integer>(px, py, avg_r_squared_fit_results[px][py], hit_count_map_ip.get(px, py));
					seeding_points.add(q);
				}
			}
		}
		Collections.sort(seeding_points, new Comparator<Quadruple<Integer, Integer, Double, Integer> >() {
			public int compare(Quadruple<Integer, Integer, Double, Integer> o1, Quadruple<Integer, Integer, Double, Integer> o2)
			{
				int comp = o2.getThird().compareTo(o1.getThird()); // sort on R^2
				//if(comp == 0)
				//{
				//	comp = o2.getFourth().compareTo(o1.getFourth()); // sort on hit count
				//}
				return comp;
			}
		});
		
		for(Quadruple<Integer, Integer, Double, Integer> q : seeding_points)
		{
			System.err.println(q.toString());
		}
		System.err.println("Number of seeding points = " + seeding_points.size());
		
		// start tracing lines
		int[][] processing_map = new int[image_width][image_height];
		ImageProcessor processing_map_ip = new ByteProcessor(image_width, image_height);
		
		Vector<Line> lines = new Vector<Line>();
		for(int pi = 0; pi < seeding_points.size(); ++pi)
		{
			// get current pixel information
			Quadruple<Integer, Integer, Double, Integer> q = seeding_points.get(pi);
			int cp_px = q.getFirst();
			int cp_py = q.getSecond();
			double cp_dlpx = line_points[cp_px][cp_py][6];
			double cp_dlpy = line_points[cp_px][cp_py][7];
			
			// get orientation of point
			double cp_nx = line_points[cp_px][cp_py][1];
			double cp_ny = line_points[cp_px][cp_py][2];
			double cp_na = getAngle(cp_nx, cp_ny);
			int cp_ni = getOrientationIndex(cp_na);
			double cp_sx = line_points[cp_px][cp_py][4];
			double cp_sy = line_points[cp_px][cp_py][5];
			double cp_sa = getAngle(cp_sx, cp_sy);
			int cp_si = getOrientationIndex(cp_sa);
			
			// RSLV: additional selection criteria, skip if normal of single does not come close to any of the two peak; i.e., if normal is between two peak, it probably is a point on an intersection
			
			// check if point has not alread been used
			if(processing_map[cp_px][cp_py] != 0)
			{
				continue; // already used, skip to next
			}
			
			// start new line from this point
			Line line = new Line();
			line.add(new Point(cp_px, cp_py, cp_dlpx, cp_dlpy));
			//point_to_line_map[cp_px][cp_py] = line;
			processing_map[cp_px][cp_py] = 255;
			processing_map_ip.set(cp_px, cp_py, 255);
			lines.add(line); // @#$%@#%@#%
			
			// *****************************************************************
		}
		
		
		// DEBUG: early exit; TODO: remove
		if(true)
		{
			return null;
		}
		
		// *********************************************************************
		
		// TODO: connect lines at/near junctions [5x5 search window]
		
		// *********************************************************************
		
		// TODO: connect lines at a slightly more global scale
		// for each line
		//   check both end point
		//      extend window into direction of line (say 15 pixels long, 10 pixels wide)
		//      if another lines endpoint is within reach
		//         connect these two lines together
		//         option: smooth out connection (i.e. drop first couple of points from endpoint)
		
		// *********************************************************************
		
		/*
		// filter short line segments
		if(FILTER_SHORT_LINE_SEGMENTS)
		{
			Vector<Line> new_lines = new Vector<Line>();
			for(Line l : lines)
			{
				if(l != null && l.size() >= LINE_LENGTH_THRESHOLD)
				{
					new_lines.add(l);
				}
			}
			lines = new_lines;
		}
		
		// *********************************************************************
		
		// draw overlay on duplicate of original image
		RoiManager roi_manager = RoiManager.getInstance();
		if(roi_manager == null) roi_manager = new RoiManager();
		//roi_manager.show(); // .setVisible(true)?
		ImageProcessor overlay_ip = ip.duplicate();
		ImageProcessor overlay2_ip = valid_line_points_magnitude_ip.duplicate();
		overlay2_ip.abs(); // RSLV: to abs() or not?
		Overlay lines_overlay = new Overlay();
		
		ImageProcessor pixelated_traces_ip = new ByteProcessor(image_width, image_height);
		
		for(int px = 0; px < image_width; ++px)
		{
			for(int py = 0; py < image_height; ++py)
			{
				if(valid_line_points_mask[px][py])
				{
					ij.gui.Line s_vector = new ij.gui.Line(px + 0.5 - 0.3 * line_points[px][py][4], py + 0.5 - 0.3 * line_points[px][py][5], px + 0.5 + 0.3 * line_points[px][py][4], py + 0.5 + 0.3 * line_points[px][py][5]);
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
			int index = 0;
			for(Point p : l)
			{
				pixelated_traces_ip.set(p.px, p.py, line_color_index);
				pxs[index] = (float)(p.px + 0.5);
				pys[index] = (float)(p.py + 0.5);
				sxs[index] = (float)(p.px + 0.5 + p.sx); // NOTE: use super-resolved coordinate
				sys[index] = (float)(p.py + 0.5 + p.sy); // NOTE: use super-resolved coordinate
				++index;
			}
			
			PolygonRoi polyline_p = new PolygonRoi(pxs, pys, Roi.POLYLINE);
			polyline_p.setStrokeColor(Color.YELLOW); // TODO: rainbow colors :D
			polyline_p.setStrokeWidth(0.0);
			lines_overlay.add(polyline_p);
			
			PolygonRoi polyline_s = new PolygonRoi(sxs, sys, Roi.POLYLINE);
			polyline_s.setStrokeColor(Color.GREEN); // TODO: rainbow colors :D
			polyline_s.setStrokeWidth(0.0);
			lines_overlay.add(polyline_s);
			
			//roi_manager.addRoi(polyline_s); // RSLV: use add(imp, roi, index)
			
			++line_color_index;
		}
		
		if(DEBUG_MODE_ENABLED)
		{
			ImagePlus overlay_imp = new ImagePlus("Steger's algorithm line traces", overlay_ip);
			overlay_imp.setOverlay(lines_overlay);
			//overlay.updateAndRepaintWindow();
			overlay_imp.show();
			
			ImagePlus overlay2_imp = new ImagePlus("Steger's algorithm line traces", overlay2_ip);
			overlay2_imp.setOverlay(lines_overlay);
			//overlay.updateAndRepaintWindow();
			overlay2_imp.show();
			
			LUT rainbow_lut = ConnectedComponents.getConnectedComponentLUT();
			pixelated_traces_ip.setLut(rainbow_lut);
			ImagePlus pixelated_traces_imp = new ImagePlus("Steger's algorithm line traces", pixelated_traces_ip);
			pixelated_traces_imp.show();
		}*/
		
		// RSLV: return image with overlay?
		//roi_manager.setVisible(true);
		//roi_manager.toFront();
		return null;
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
}
