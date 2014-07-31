
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

// import Jama classes
import Jama.Matrix;
import Jama.EigenvalueDecomposition;

// import own classes
//import algorithms.Stegers;
import filters.DerivativeOfGaussian;
import filters.ConnectedComponents; // for rainbow LUT
import core.Point;
import core.Line;
import utils.Tuple;


/**
 *	
 */
public class Stegers_Algorithm implements PlugIn
{
	/**
	 *	Constants
	 */
	public enum Mode { ABS_MAX, MAX, ABS_MIN, MIN };
	public static final String[] MODES_S = new String[]{"Abs. max", "Max", "Abs. min", "Min"};
	
	public static final Mode DEFAULT_MODE = Mode.ABS_MAX;
	public static final int DEFAULT_MODE_I = 0;
	public static final String DEFAULT_MODE_S = MODES_S[DEFAULT_MODE_I];
	
	public static Mode MODE = DEFAULT_MODE;
	public static int MODE_I = DEFAULT_MODE_I;
	public static String MODE_S = DEFAULT_MODE_S;
	
	public static final double DEFAULT_SIGMA = 1.8;
	public static final double DEFAULT_SIGMA_RANGE = 0.0;
	public static final int DEFAULT_SIGMA_STEPS = 1;
	public static final double DEFAULT_LINEPOINT_THRESHOLD = 0.5;
	
	public static double SIGMA = DEFAULT_SIGMA;
	public static double SIGMA_RANGE = DEFAULT_SIGMA_RANGE;
	public static int SIGMA_STEPS = DEFAULT_SIGMA_STEPS;
	public static double LINEPOINT_THRESHOLD = DEFAULT_LINEPOINT_THRESHOLD;
	
	public static final boolean DEFAULT_USE_PRESET_USER_THRESHOLDS = false;
	public static final double DEFAULT_UPPER_THRESHOLD = 5000.0; // RSLV: convert to percentage of histogram?
	public static final double DEFAULT_LOWER_THRESHOLD = 1000.0; // RSLV: convert to percentage of histogram?
	
	public static boolean USE_PRESET_USER_THRESHOLDS = DEFAULT_USE_PRESET_USER_THRESHOLDS;
	public static double UPPER_THRESHOLD = DEFAULT_UPPER_THRESHOLD;
	public static double LOWER_THRESHOLD = DEFAULT_LOWER_THRESHOLD;
	
	public static final boolean DEFAULT_FILTER_MULTIRESPONSE = true;
	public static final double DEFAULT_MULTIRESPONSE_THRESHOLD = 360.0 / 18.0; // about 20 degrees parallelism
	
	public static boolean FILTER_MULTIRESPONSE = DEFAULT_FILTER_MULTIRESPONSE;
	public static double MULTIRESPONSE_THRESHOLD = DEFAULT_MULTIRESPONSE_THRESHOLD;
	
	public static final boolean DEFAULT_PREVENT_SPLINTER_FRAGMENTS = true;
	public static final int DEFAULT_SPLINTER_FRAGMENT_THRESHOLD = 3;
	
	public static boolean PREVENT_SPLINTER_FRAGMENTS = DEFAULT_PREVENT_SPLINTER_FRAGMENTS;
	public static int SPLINTER_FRAGMENT_THRESHOLD = DEFAULT_SPLINTER_FRAGMENT_THRESHOLD;
	
	public static final boolean DEFAULT_FILTER_SHORT_LINE_SEGMENTS = true;
	public static final int DEFAULT_LINE_LENGTH_THRESHOLD = 10;
	
	public static boolean FILTER_SHORT_LINE_SEGMENTS = DEFAULT_FILTER_SHORT_LINE_SEGMENTS;
	public static int LINE_LENGTH_THRESHOLD = DEFAULT_LINE_LENGTH_THRESHOLD;
	
	public static final boolean DEFAULT_INCLUDE_JUNCTIONS_IN_LINE = true;
	
	public static boolean INCLUDE_JUNCTIONS_IN_LINE = DEFAULT_INCLUDE_JUNCTIONS_IN_LINE;
	
	public static final boolean DEFAULT_INCLUDE_FIRST_PROCESSED_IN_LINE = true;
	
	public static boolean INCLUDE_FIRST_PROCESSED_IN_LINE = DEFAULT_INCLUDE_FIRST_PROCESSED_IN_LINE;
	
	public static final double DEFAULT_COST_FUNCTION_WEIGHT = 1.0;
	
	public static double COST_FUNCTION_WEIGHT = DEFAULT_COST_FUNCTION_WEIGHT;
	
	public static final boolean DEFAULT_USE_WEINGARTEN_MATRIX = false;
	
	public static boolean USE_WEINGARTEN_MATRIX = DEFAULT_USE_WEINGARTEN_MATRIX;
	
	public static final boolean DEFAULT_DEBUG_MODE_ENABLED = true;
	
	public static boolean DEBUG_MODE_ENABLED = DEFAULT_DEBUG_MODE_ENABLED;
	
	public static final boolean DEFAULT_ENABLE_LINKING = true;
	
	public static boolean ENABLE_LINKING = DEFAULT_ENABLE_LINKING;
	
	public static final double DEFAULT_MAX_BENDING_ANGLE = 30.0;
	
	public static double MAX_BENDING_ANGLE = DEFAULT_MAX_BENDING_ANGLE;
	
	// ////////////////////////////////////////////////////////////////////////
	
	/**
	 *	Constructor
	 */
	public Stegers_Algorithm()
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
		gd.addChoice("Mode", MODES_S, Prefs.get("stegers.mode_s", DEFAULT_MODE_S));
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addNumericField("Sigma", Prefs.get("stegers.sigma", DEFAULT_SIGMA), 2);
//		gd.addNumericField("Sigma_range", Prefs.get("stegers.sigma_range", DEFAULT_SIGMA_RANGE), 2);
//		gd.addNumericField("Sigma_steps", Prefs.get("stegers.sigma_steps", DEFAULT_SIGMA_STEPS), 0);
		gd.addNumericField("Linepoint_threshold", Prefs.get("stegers.linepoint_threshold", DEFAULT_LINEPOINT_THRESHOLD), 2);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Use_preset_user_thresholds", Prefs.get("stegers.preset_user_thresholds", DEFAULT_USE_PRESET_USER_THRESHOLDS));
		gd.addNumericField("Upper_threshold", Prefs.get("stegers.upper_threshold", DEFAULT_UPPER_THRESHOLD), 0);
		gd.addNumericField("Lower_threshold", Prefs.get("stegers.lower_threshold", DEFAULT_LOWER_THRESHOLD), 0);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Filter_multiresponse", Prefs.get("stegers.filter_multiresponse", DEFAULT_FILTER_MULTIRESPONSE));
		gd.addNumericField("Multireponse_angle_threshold", Prefs.get("stegers.multiresponse_threshold", DEFAULT_MULTIRESPONSE_THRESHOLD), 2);
		
		gd.addCheckbox("Prevent_splinter_fragments", Prefs.get("stegers.prevent_splinter_fragments", DEFAULT_PREVENT_SPLINTER_FRAGMENTS));
		gd.addNumericField("Splinter_fragment_threshold", Prefs.get("stegers.splinter_fragment_threshold", DEFAULT_SPLINTER_FRAGMENT_THRESHOLD), 0);
		
		gd.addCheckbox("Filter_short_line_segments", Prefs.get("stegers.filter_short_line_segments", DEFAULT_FILTER_SHORT_LINE_SEGMENTS));
		gd.addNumericField("Line_length_threshold", Prefs.get("stegers.line_length_threshold", DEFAULT_LINE_LENGTH_THRESHOLD), 0);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Include_junctions", Prefs.get("stegers.include_junctions", DEFAULT_INCLUDE_JUNCTIONS_IN_LINE));
		gd.addCheckbox("Include_first_processed", Prefs.get("stegers.include_first_processed", DEFAULT_INCLUDE_FIRST_PROCESSED_IN_LINE));
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addNumericField("Cost_funtion_weight", Prefs.get("stegers.cost_function_weight", DEFAULT_COST_FUNCTION_WEIGHT), 2);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
//		gd.addCheckbox("Use_Weingarten_matrix", Prefs.get("stegers.weingarten_matrix", DEFAULT_USE_WEINGARTEN_MATRIX));
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Enable_debug_mode", Prefs.get("stegers.debug_mode", DEFAULT_DEBUG_MODE_ENABLED));
		
		gd.addCheckbox("Enable_linking", Prefs.get("stegers.enable_linking", DEFAULT_ENABLE_LINKING));
		gd.addNumericField("Max_bending_angle", Prefs.get("stegers.max_bending_angle", DEFAULT_MAX_BENDING_ANGLE), 0);
		
		gd.setOKLabel("Trace");
		
		gd.showDialog();
		if(gd.wasCanceled()) return;
		
		// retrieve parameters
		MODE_I = gd.getNextChoiceIndex();
		MODE_S = MODES_S[MODE_I];
		MODE = Mode.values()[MODE_I];
		
		SIGMA = gd.getNextNumber();
//		SIGMA_RANGE = gd.getNextNumber();
//		SIGMA_STEPS = (int)gd.getNextNumber();
		LINEPOINT_THRESHOLD = gd.getNextNumber();
		
		USE_PRESET_USER_THRESHOLDS = gd.getNextBoolean();
		UPPER_THRESHOLD = gd.getNextNumber();
		LOWER_THRESHOLD = gd.getNextNumber();
		
		FILTER_MULTIRESPONSE = gd.getNextBoolean();
		MULTIRESPONSE_THRESHOLD = gd.getNextNumber();
		
		PREVENT_SPLINTER_FRAGMENTS = gd.getNextBoolean();
		SPLINTER_FRAGMENT_THRESHOLD = (int)gd.getNextNumber();
		
		FILTER_SHORT_LINE_SEGMENTS = gd.getNextBoolean();
		LINE_LENGTH_THRESHOLD = (int)gd.getNextNumber();
		
		INCLUDE_JUNCTIONS_IN_LINE = gd.getNextBoolean();
		INCLUDE_FIRST_PROCESSED_IN_LINE = gd.getNextBoolean();
		
		COST_FUNCTION_WEIGHT = gd.getNextNumber();
		
//		USE_WEINGARTEN_MATRIX = gd.getNextBoolean();
		
		DEBUG_MODE_ENABLED = gd.getNextBoolean();
		
		ENABLE_LINKING = gd.getNextBoolean();
		
		MAX_BENDING_ANGLE = gd.getNextNumber();
		
		// store parameters in preferences
		Prefs.set("stegers.mode_s", MODE_S);
		Prefs.set("stegers.sigma", SIGMA);
//		Prefs.set("stegers.sigma_range", SIGMA_RANGE);
//		Prefs.set("stegers.sigma_steps", SIGMA_STEPS);
		Prefs.set("stegers.linepoint_threshold", LINEPOINT_THRESHOLD);
		Prefs.set("stegers.preset_user_thresholds", USE_PRESET_USER_THRESHOLDS);
		Prefs.set("stegers.upper_threshold", UPPER_THRESHOLD);
		Prefs.set("stegers.lower_threshold", LOWER_THRESHOLD);
		Prefs.set("stegers.filter_multiresponse", FILTER_MULTIRESPONSE);
		Prefs.set("stegers.multiresponse_threshold", MULTIRESPONSE_THRESHOLD);
		Prefs.set("stegers.prevent_splinter_fragments", PREVENT_SPLINTER_FRAGMENTS);
		Prefs.set("stegers.splinter_fragment_threshold", SPLINTER_FRAGMENT_THRESHOLD);
		Prefs.set("stegers.filter_short_line_segments", FILTER_SHORT_LINE_SEGMENTS);
		Prefs.set("stegers.line_length_threshold", LINE_LENGTH_THRESHOLD);
		Prefs.set("stegers.include_junctions", INCLUDE_JUNCTIONS_IN_LINE);
		Prefs.set("stegers.include_first_processed", INCLUDE_FIRST_PROCESSED_IN_LINE);
		Prefs.set("stegers.cost_function_weight", COST_FUNCTION_WEIGHT);
//		Prefs.set("stegers.weingarten_matrix", USE_WEINGARTEN_MATRIX);
		Prefs.set("stegers.debug_mode", DEBUG_MODE_ENABLED);
		Prefs.set("stegers.enable_linking", ENABLE_LINKING);
		Prefs.set("stegers.max_bending_angle", MAX_BENDING_ANGLE);
		
		MAX_BENDING_ANGLE = Math.toRadians(MAX_BENDING_ANGLE);
		
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
		stegersAlgorithm(imp, DEFAULT_SIGMA);
	}
	
	public static void stegersAlgorithm(ImagePlus imp, double sigma)
	{
		stegersAlgorithm(imp, sigma, DEFAULT_LINEPOINT_THRESHOLD);
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
		
		// STEP 1: calculate Gaussian derivatives, Hessian matrix, eigenvalues and eigenvector, determine Taylor polynomical response function peak position
		
		// store intermediate results in table
		//	[px][py][0] = lambda1_magnitude		n(t)
		//	[px][py][1] = lambda1_direction_x	n_x(t)
		//	[px][py][2] = lambda1_direction_y	n_y(t)
		//	[px][py][3] = lambda2_magnitude		s(t)
		//	[px][py][4] = lambda2_direction_x	s_x(t)
		//	[px][py][5] = lambda2_direction_y	s_y(t)
		//	[px][py][6] = super-resolved_x		t_x, or dlpx
		//	[px][py][7] = super-resolved_y		t_y, or dlpy
		double[][][] line_points = new double[image_width][image_height][8];
		
		ImageProcessor dx = new FloatProcessor(image_width, image_height);
		ImageProcessor dy = new FloatProcessor(image_width, image_height);
		ImageProcessor dxdx = new FloatProcessor(image_width, image_height);
		ImageProcessor dxdy = new FloatProcessor(image_width, image_height);
		ImageProcessor dydx = new FloatProcessor(image_width, image_height);
		ImageProcessor dydy = new FloatProcessor(image_width, image_height);
		
		double[][] sigma_map = new double[image_width][image_height];
		ImageProcessor sigma_map_ip = new FloatProcessor(image_width, image_height);
		
		// search scale space: maximize selected eigenvalue's response
		double lower_sigma = sigma - SIGMA_RANGE;
		double upper_sigma = sigma + SIGMA_RANGE;
		double step_sigma = (upper_sigma - lower_sigma) / (SIGMA_STEPS > 1 ? SIGMA_STEPS - 1 : SIGMA_STEPS); // 2 * SIGMA_RANGE / SIGMA_STEPS
		for(int k = 0; k < SIGMA_STEPS; ++k) //(double current_sigma = lower_sigma; current_sigma <= upper_sigma; current_sigma += step_sigma)
		{
			// calculate current sigma
			double current_sigma = lower_sigma + k * step_sigma;
			System.err.println("current_sigma = " + current_sigma);
			
			// calculate derivatives of gaussian from image
			ImageProcessor dx_t = DerivativeOfGaussian.derivativeX(ip, current_sigma);
			ImageProcessor dy_t = DerivativeOfGaussian.derivativeY(ip, current_sigma);
			ImageProcessor dxdx_t = DerivativeOfGaussian.derivativeXX(ip, current_sigma);
			ImageProcessor dxdy_t = DerivativeOfGaussian.derivativeXY(ip, current_sigma);
			ImageProcessor dydx_t = dxdy_t;//DerivativeOfGaussian.derivativeYX(ip, current_sigma);
			ImageProcessor dydy_t = DerivativeOfGaussian.derivativeYY(ip, current_sigma);
			
			// calculate line points from eigenvalues and eigenvectors based on Hessian matrix
			for(int py = 0; py < image_height; ++py)
			{
				for(int px = 0; px < image_width; ++px)
				{
					// construct Hessian matrix in Jama
					Matrix m = null; // NOTE: beware of null-pointer exceptions!
					if(USE_WEINGARTEN_MATRIX)
					{
						double dx_squared = dx_t.getf(px, py) * dx_t.getf(px, py);
						double dy_squared = dy_t.getf(px, py) * dy_t.getf(px, py);
						double dx_times_dy = dx_t.getf(px, py) * dy_t.getf(px, py);
						double dy_times_dx = dx_times_dy; // dy.getf(px, py) * dx.getf(px, py);
						
						Matrix f1 = new Matrix(2, 2, 0); // 2x2 RC matrix with zeros
						f1.set(0, 0, 1 + dx_squared);
						f1.set(0, 1, dx_times_dy);
						f1.set(1, 0, dy_times_dx);
						f1.set(1, 1, 1 + dy_squared);
						
						Matrix f2 = new Matrix(2, 2, 0); // 2x2 RC matrix with zeros
						f2.set(0, 0, dxdx_t.getf(px, py));
						f2.set(0, 1, dxdy_t.getf(px, py));
						f2.set(1, 0, dydx_t.getf(px, py));
						f2.set(1, 1, dydy_t.getf(px, py));
						
						m = f2.times(f1.inverse());
						//m.timesEquals(-1 / Math.sqrt(1 + dx_squared + dy_squared)); // inplace scalar multiplication: same as m = m.times(-1 / Math.sqrt(1 + dx_squared + dy_squared));
					}
					else
					{
						// use Hessian matrix
						m = new Matrix(2, 2, 0); // 2x2 RC matrix with zeros
						m.set(0, 0, dxdx_t.getf(px, py));
						m.set(0, 1, dxdy_t.getf(px, py));
						m.set(1, 0, dydx_t.getf(px, py));
						m.set(1, 1, dydy_t.getf(px, py));
						//System.err.println("Cond="+h.cond());
					}
					
					// compute eigenvalues and eigenvectors
					EigenvalueDecomposition evd = m.eig();
					Matrix d = evd.getD();
					Matrix v = evd.getV();
					
					// determine largest absolute (perpendicular -> n(t)) and smallest absolute (parallel -> s(t)) eigenvalue and corresponding eigenvector
					double selected_eigenvalue = 0.0; // |n(t)|
					double nx = 0.0; // n(t) -> perpendicular to s(t)
					double ny = 0.0; // n(t) -> perpendicular to s(t)
					double second_eigenvalue = 0.0;
					double sx = 0.0;
					double sy = 0.0;
					if((MODE == Mode.ABS_MAX && Math.abs(d.get(0,0)) >= Math.abs(d.get(1,1))) // Stegers*: absolute maximum
					|| (MODE == Mode.MAX && d.get(0,0) >= d.get(1,1)) // real maximum
					|| (MODE == Mode.ABS_MIN && Math.abs(d.get(0,0)) <= Math.abs(d.get(1,1))) // absolute minimum
					|| (MODE == Mode.MIN && d.get(0,0) <= d.get(1,1))) // real minimum
					{
						selected_eigenvalue = d.get(0,0);
						nx = v.get(0,0);
						ny = v.get(1,0);
						second_eigenvalue = d.get(1,1);
						sx = v.get(0,1);
						sy = v.get(1,1);
					}
					else
					{
						selected_eigenvalue = d.get(1,1);
						nx = v.get(0,1);
						ny = v.get(1,1);
						second_eigenvalue = d.get(0,0);
						sx = v.get(0,0);
						sy = v.get(1,0);
					}
					
					// reorient all eigenvecors in same direction
					if(ny < 0)
					{
						nx = -nx;
						ny = -ny;
					}
					
					if(sy < 0)
					{
						sx = -sx;
						sy = -sy;
					}
					
					// calculate position of line point and filter line points
					double t = ((dx_t.getf(px,py)*nx + dy_t.getf(px,py)*ny) / (dxdx_t.getf(px,py)*nx*nx + dxdy_t.getf(px,py)*nx*ny + dydx_t.getf(px,py)*ny*nx + dydy_t.getf(px,py)*ny*ny)); // NOTE: removed '-' from start of equation!!
					double dlpx = t*nx;
					double dlpy = t*ny;
					
					// store line point information (only if larger than previous scale space search)
					// RSLV: always absolute value?
					if(Math.abs(selected_eigenvalue) > Math.abs(line_points[px][py][0]))
					{
						line_points[px][py][0] = selected_eigenvalue;
						line_points[px][py][1] = nx;
						line_points[px][py][2] = ny;
						line_points[px][py][3] = second_eigenvalue;
						line_points[px][py][4] = sx;
						line_points[px][py][5] = sy;
						line_points[px][py][6] = dlpx;
						line_points[px][py][7] = dlpy;
						
						dx.setf(px, py, dx_t.getf(px, py));
						dy.setf(px, py, dy_t.getf(px, py));
						dxdx.setf(px, py, dxdx_t.getf(px, py));
						dxdy.setf(px, py, dxdy_t.getf(px, py));
						dydx.setf(px, py, dydx_t.getf(px, py));
						dydy.setf(px, py, dydy_t.getf(px, py));
						
						sigma_map[px][py] = current_sigma;
						sigma_map_ip.setf(px, py, (float)current_sigma);
					}
				}
			}
		} // end scale space search
		
		if(DEBUG_MODE_ENABLED)
		{
			ImagePlus sigma_map_imp = new ImagePlus("Sigma map", sigma_map_ip);
			sigma_map_imp.resetDisplayRange();
//			sigma_map_imp.show();
		}
		
		// *********************************************************************
		
		// STEP 2: filter line points using Taylor polynomial response function fit and distance criterion of peak from pixel center
		boolean[][] valid_line_points_mask = new boolean[image_width][image_height];
		double[][] valid_line_points_magnitude = new double[image_width][image_height];
		
		// TEMP: debug images
		ImageProcessor valid_line_points_mask_ip = new ByteProcessor(line_points.length, line_points[0].length); // use width, height from double array
		ImageProcessor valid_line_points_magnitude_ip = new FloatProcessor(line_points.length, line_points[0].length); // use width, height from double array
		
		for(int px = 0; px < image_width; ++px)
		{
			for(int py = 0; py < image_height; ++py)
			{
				if(Math.abs(line_points[px][py][6]) <= linepoint_threshold && Math.abs(line_points[px][py][7]) <= linepoint_threshold)
				{
					valid_line_points_mask[px][py] = true;
					valid_line_points_magnitude[px][py] = line_points[px][py][0];
					
					// TEMP: debug images
					valid_line_points_mask_ip.set(px, py, 255);
					valid_line_points_magnitude_ip.setf(px, py, (float)line_points[px][py][0]);
				}
				else
				{
					valid_line_points_mask[px][py] = false;
					valid_line_points_magnitude[px][py] = 0.0;
					
					// TEMP: debug images
					valid_line_points_mask_ip.setf(px, py, 0);
					valid_line_points_magnitude_ip.setf(px, py, 0.0f);
				}
			}
		}
		
		// show debug images
		if(DEBUG_MODE_ENABLED)
		{
			ImagePlus valid_line_points_mask_image = new ImagePlus("Valid line points mask", valid_line_points_mask_ip);
			valid_line_points_mask_image.resetDisplayRange();
			valid_line_points_mask_image.show();
			
			// TEMP: debug images
			ImagePlus valid_line_points_magnitude_image = new ImagePlus("Valid line points magnitude", valid_line_points_magnitude_ip);
			valid_line_points_magnitude_image.resetDisplayRange();
			valid_line_points_magnitude_image.show();
		}
		
		// *********************************************************************
		
		// STEP 3: trace lines; requires two user specified thresholds
		
		// get user threshold levels (interactive)
		if(!USE_PRESET_USER_THRESHOLDS)
		{
			ImageProcessor threshold_tmp_ip = valid_line_points_magnitude_ip.duplicate();
			threshold_tmp_ip.abs();
			//threshold_tmp_ip.setMinAndMax(0, 65535);
			threshold_tmp_ip.resetMinAndMax();
			double threshold_tmp_ip_min = threshold_tmp_ip.getMin();
			double threshold_tmp_ip_max = threshold_tmp_ip.getMax();
			System.err.println("ip.getMin()="+threshold_tmp_ip_min);
			System.err.println("ip.getMax()="+threshold_tmp_ip_max);
			final double threshold_map_scale_factor = threshold_tmp_ip.getMax() / 256;
			final ImageProcessor threshold_ip = threshold_tmp_ip.convertToByteProcessor(); // final to make accessible in anonymous inner class
			final ImagePlus threshold_imp = new ImagePlus("Threshold map", threshold_ip); // final to make accessible in anonymous inner class
			threshold_imp.resetDisplayRange();
			threshold_imp.show();
			
			GenericDialog thresholds_gd = new GenericDialog("User threshold levels"); // RSLV: make NonBlockingGenericDialog();
			double min_intensity_value = threshold_tmp_ip_min; // 0
			double max_intensity_value = threshold_tmp_ip_max; //(int)Math.pow(2, ip.getBitDepth())-1;
			thresholds_gd.addSlider("Upper threshold", min_intensity_value, max_intensity_value, UPPER_THRESHOLD); // RSLV: limit?
			thresholds_gd.addSlider("Lower threshold", min_intensity_value, max_intensity_value, LOWER_THRESHOLD); // RSLV: limit?
			
			DialogListener dlg_listener = new DialogListener(){
				@Override
				public boolean dialogItemChanged(GenericDialog gd, java.awt.AWTEvent e)
				{
					// get threshold parameters
					double ut = gd.getNextNumber();
					double lt = gd.getNextNumber();
					
					// make sure upper threshold is always greater than or equal to lower threshold
					if(lt > ut)
					{
						// TODO: update sliders in dialog; how?
						return false;
					}
					
					// use luts to highlight regions
					// NOTE: LUTs are 8-bit, so map is only an approximation
					int but = (int)(ut / threshold_map_scale_factor);
					int blt = (int)(lt / threshold_map_scale_factor);
					
					byte[] threshold_lut_r = new byte[256];
					byte[] threshold_lut_g = new byte[256];
					byte[] threshold_lut_b = new byte[256];
					
					// create lut
					for(int i = 0; i < 256; ++i)
					{
						if(i < blt)
						{
							// retain gray scale
							threshold_lut_r[i] = (byte)i;
							threshold_lut_g[i] = (byte)i;
							threshold_lut_b[i] = (byte)i;
						}
						else if(i < but)
						{
							// set to lower threshold colour
							threshold_lut_r[i] = (byte)255;
							threshold_lut_g[i] = (byte)255;
							threshold_lut_b[i] = (byte)0;
						}
						else
						{
							// set to upper threshold colour
							threshold_lut_r[i] = (byte)255;
							threshold_lut_g[i] = (byte)0;
							threshold_lut_b[i] = (byte)0;
						}
					}
					LUT threshold_lut = new LUT(threshold_lut_r, threshold_lut_g, threshold_lut_b);
					
					// set LUT and update image
					threshold_ip.setLut(threshold_lut);
					threshold_imp.updateAndDraw();
					
					return true; // true == accepted values, false == incorrect values
				}
			};
			thresholds_gd.addDialogListener(dlg_listener);
			dlg_listener.dialogItemChanged(thresholds_gd, null); // force update of lut
			
			thresholds_gd.setOKLabel("Continue tracing");
			thresholds_gd.setCancelLabel("Cancel tracing");
			
			// focus on threshold imp window, (sometimes get lost behind other debug image windows)
			threshold_imp.getWindow().setVisible(true);
			threshold_imp.getWindow().toFront();
			
			thresholds_gd.showDialog();
			
			if(thresholds_gd.wasCanceled())
			{
				return null;
			}
			
			// get user specified threshold
			UPPER_THRESHOLD = thresholds_gd.getNextNumber();
			LOWER_THRESHOLD = thresholds_gd.getNextNumber();
			System.err.println("UPPER_THRESHOLD="+UPPER_THRESHOLD);
			System.err.println("LOWER_THRESHOLD="+LOWER_THRESHOLD);
			// store user thresholds
			Prefs.set("stegers.upper_threshold", UPPER_THRESHOLD);
			Prefs.set("stegers.lower_threshold", LOWER_THRESHOLD);
			
			// close intermediate image
			if(!DEBUG_MODE_ENABLED)
			{
				threshold_imp.close();
			}
		}
		
		// STEP 3a: first find all valid lines above threshold in descending order
		Vector<Vector<Double> > high_threshold_points = new Vector<Vector<Double> >();
		for(int px = 0; px < image_width; ++px)
		{
			for(int py = 0; py < image_height; ++py)
			{
				// use only valid line points
				if(Math.abs(valid_line_points_magnitude[px][py]) >= UPPER_THRESHOLD)
				{
					// create point vector
					Vector<Double> p = new Vector<Double>(3);
					p.add(new Double(px));
					p.add(new Double(py));
					p.add(new Double(Math.abs(valid_line_points_magnitude[px][py])));
					high_threshold_points.add(p);
				}
			}
		}
		Collections.sort(high_threshold_points, new Comparator<Vector<Double> >() {
			public int compare(Vector<Double> o1, Vector<Double> o2)
			{
				return o2.get(2).compareTo(o1.get(2)); // descending order
			}
		});
		
		// intermediate containers and flags
		HashMap<Integer, Vector<Line> > junctions_map = new HashMap<Integer, Vector<Line> >(); // RSLV: create tuple class for key; Tuple<Integer, Integer> for x,y-coordinate?
		Line[][] point_to_line_map = new Line[image_width][image_height]; // NOTE: for computational efficiency of finding which line belongs to a point
		int[][] processing_map = new int[image_width][image_height];
		int DEFAULT = 0x00;
		int AVAILABLE = 0x00; // NOTE: same as DEFAULT
		int LIMITED = 0x01;
		int PROCESSED = 0x02;
		int LINEPOINT = 0x04;
		int JUNCTION = 0x08;
		
		ImageProcessor processing_map_ip = new ByteProcessor(image_width, image_height);
		
		// STEP 3c: start tracing from high threshold line points
		Vector<Line> lines = new Vector<Line>();
		for(int pi = 0; pi < high_threshold_points.size(); ++pi)
		{
			// get current pixel information
			int cp_px = high_threshold_points.get(pi).get(0).intValue();
			int cp_py = high_threshold_points.get(pi).get(1).intValue();
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
			
			// check if point has not alread been used
			if(processing_map[cp_px][cp_py] != AVAILABLE)
			{
				continue; // already used, skip to next
			}
			
			// start new line from this point
			Line line = new Line();
			line.add(new Point(cp_px, cp_py, cp_dlpx, cp_dlpy));
			point_to_line_map[cp_px][cp_py] = line;
			processing_map[cp_px][cp_py] = LINEPOINT;
			processing_map_ip.set(cp_px, cp_py, LINEPOINT);
			lines.add(line); // @#$%@#%@#%
			
			// *****************************************************************
			
			//for(int trace_direction = 0; trace_direction != 180; trace_direction = 180)
			double[] trace_directions = new double[]{0, 180};
			for(int tdi = 0; tdi < trace_directions.length; ++tdi)
			{
				double trace_direction = trace_directions[tdi];
				System.err.println("Trace in direction " + trace_direction);
				
				// trace line in first direction
				int tp_px = cp_px;
				int tp_py = cp_py;
				
				//double tp_sxo = cp_sx;
				//double tp_syo = cp_sy;
				double tp_sao = (cp_sa+trace_direction)%360.0;
				//int tp_sio = cp_si;
				
				boolean tracing = true;
				while(tracing)
				{
					// get information on current trace point
					double tp_dlpx = line_points[tp_px][tp_py][6];
					double tp_dlpy = line_points[tp_px][tp_py][7];
					
					// get orientation of point
					double tp_nx = line_points[tp_px][tp_py][1];
					double tp_ny = line_points[tp_px][tp_py][2];
					double tp_na = getAngle(tp_nx, tp_ny);
					int tp_ni = getOrientationIndex(tp_na);
					double tp_sx = line_points[tp_px][tp_py][4];
					double tp_sy = line_points[tp_px][tp_py][5];
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
						
						// NOTE: added check for lower threshold here, rather then a couple of lines below
						// only process valid line point
						// and above (or equal to) lower user-specified threshold
						if(valid_line_points_mask[dp_px][dp_py] && Math.abs(line_points[dp_px][dp_py][0]) >= LOWER_THRESHOLD)
						{
							// get next neighbour information
							double dp_dlpx = line_points[dp_px][dp_py][6];
							double dp_dlpy = line_points[dp_px][dp_py][7];
							
							// get orientation of point
							double dp_nx = line_points[dp_px][dp_py][1];
							double dp_ny = line_points[dp_px][dp_py][2];
							double dp_na = getAngle(dp_nx, dp_ny); //Math.atan2(dp_ny, dp_nx);
							int dp_ni = getOrientationIndex(dp_na);
							double dp_sx = line_points[dp_px][dp_py][4];
							double dp_sy = line_points[dp_px][dp_py][5];
							double dp_sa = getAngle(dp_sx, dp_sy); //Math.atan2(dp_sy, dp_sx);
							int dp_si = getOrientationIndex(dp_sa);
							
							// calclate cost function
							double ddist_x = (dp_px + dp_dlpx) - (tp_px + tp_dlpx);
							double ddist_y = (dp_py + dp_dlpy) - (tp_py + tp_dlpy);
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
							dtheta *= Math.PI / 180.0; // deg to rad
							
							if(dtheta > MAX_BENDING_ANGLE)
							{
								System.err.println("  Exceeding maximum bending angle, skipping point");
								continue;
							}
							
							double cost = ddist + COST_FUNCTION_WEIGHT * dtheta;
							
							System.err.println("  cost=" + ddist + "+" +  COST_FUNCTION_WEIGHT + "*" + dtheta + "=" + cost);
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
						else
						{
							System.err.println("  Next point is not a valid line point");
							System.err.println("  Or below lower user-specified threshold");
						}
					}
					
					// *********************************************************
					
					// skip if no best match found (outside of image window or no valid line point)
					if(np_px == -1 || np_py == -1)
					{
						System.err.println("Did not find any suitable neighbour");
						tracing = false;
						continue; //break;
					}
					
					System.err.println("Best neighbour at ("+np_px+","+np_py+")");
					
					// get next neighbour information
					double np_dlpx = line_points[np_px][np_py][6];
					double np_dlpy = line_points[np_px][np_py][7];
					
					// get orientation of point
					double np_nx = line_points[np_px][np_py][1];
					double np_ny = line_points[np_px][np_py][2];
					double np_na = getAngle(np_nx, np_ny); //Math.atan2(cp_ny, cp_nx);
					int np_ni = getOrientationIndex(np_na);
					double np_sx = line_points[np_px][np_py][4];
					double np_sy = line_points[np_px][np_py][5];
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
					if(processing_map[np_px][np_py] == AVAILABLE || processing_map[np_px][np_py] == LIMITED)
					{
						System.err.println("  Point is available");
						
						// NOTE: check for below lower threshold is moved up where selecting the best neighbour; this avoids having two lines that meet, but not connected
//						// skip if next neighbour below lowest threshold level
//						double np_pv = line_points[np_px][np_py][0];
//						if(Math.abs(np_pv) < LOWER_THRESHOLD)
//						{
//							System.err.println("  Next point is below lower threshold");
//							// stop tracing
//							tracing = false;
//						}
//						else
//						{
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
							processing_map[np_px][np_py] = LINEPOINT;
							processing_map_ip.set(np_px, np_py, LINEPOINT);
							
							if(FILTER_MULTIRESPONSE)
							{
								// mark perpendicular neighbours with similar orientation as processed
								int[] pp = getNextPixel(np_ni);
								int pp1_dx = pp[0];
								int pp1_dy = pp[1];
								int pp2_dx = -pp1_dx; // NOTE: pp2_dx = -pp1_dx
								int pp2_dy = -pp1_dy; // NOTE: pp2_dy = -pp1_dy
								
								// parallel point 1
								int pp1_px = np_px + pp1_dx;
								int pp1_py = np_py + pp1_dy;
								if(pp1_px >= 0 && pp1_px < image_width && pp1_py >= 0 && pp1_py < image_height)
								{
									double pp1_nx = line_points[pp1_px][pp1_py][1];
									double pp1_ny = line_points[pp1_px][pp1_py][2];
									double pp1_na = getAngle(pp1_nx, pp1_ny); //Math.atan2(pp1_ny, pp1_nx);
									int pp1_ni = getOrientationIndex(pp1_na);
									double pp1_sx = line_points[pp1_px][pp1_py][4];
									double pp1_sy = line_points[pp1_px][pp1_py][5];
									double pp1_sa = getAngle(pp1_sx, pp1_sy); //Math.atan2(pp1_sy, pp1_sx);
									int pp1_si = getOrientationIndex(pp1_sa);
									
									// mark as processed if roughly the same orientation as current point
									if((processing_map[pp1_px][pp1_py] == AVAILABLE || processing_map[pp1_px][pp1_py] == LIMITED) && Math.abs((pp1_na % 180.0) - (np_na % 180.0)) <= MULTIRESPONSE_THRESHOLD)
									{
										System.err.println("Checkpoint D");
										processing_map[pp1_px][pp1_py] = PROCESSED;
										processing_map_ip.set(pp1_px, pp1_py, PROCESSED);
									}
								}
								
								// parallel point 2
								int pp2_px = np_px + pp2_dx;
								int pp2_py = np_py + pp2_dy;
								if(pp2_px >= 0 && pp2_px < image_width && pp2_py >= 0 && pp2_py < image_height)
								{
									double pp2_nx = line_points[pp2_px][pp2_py][1];
									double pp2_ny = line_points[pp2_px][pp2_py][2];
									double pp2_na = getAngle(pp2_nx, pp2_ny); //Math.atan2(pp2_ny, pp2_nx);
									int pp2_ni = getOrientationIndex(pp2_na);
									double pp2_sx = line_points[pp2_px][pp2_py][4];
									double pp2_sy = line_points[pp2_px][pp2_py][5];
									double pp2_sa = getAngle(pp2_sx, pp2_sy); //Math.atan2(pp2_sy, pp2_sx);
									int pp2_si = getOrientationIndex(pp2_sa);
									
									// mark as processed if roughly the same orientation as current point
									if((processing_map[pp2_px][pp2_py] == AVAILABLE || processing_map[pp2_px][pp2_py] == LIMITED) && Math.abs((pp2_na % 180.0) - (np_na % 180.0)) <= MULTIRESPONSE_THRESHOLD)
									{
										System.err.println("Checkpoint E");
										processing_map[pp2_px][pp2_py] = PROCESSED;
										processing_map_ip.set(pp2_px, pp2_py, PROCESSED);
									}
								}
							}
//						}
					}
					else if(processing_map[np_px][np_py] == PROCESSED)
					{
						System.err.println("  Point is already processed, stop tracing");
						
						// add point
						if(INCLUDE_FIRST_PROCESSED_IN_LINE)
						{
							if(trace_direction == 0)
							{
								line.addFirst(new Point(np_px, np_py, np_dlpx, np_dlpy));
							}
							else
							{
								line.addLast(new Point(np_px, np_py, np_dlpx, np_dlpy));
							}
						}
						
						// RSLV: mark as line point?
						processing_map[np_px][np_py] = LINEPOINT;
						
						// add to point to line map
						point_to_line_map[np_px][np_py] = line;
						
						// stop tracing
						tracing = false;
					}
					else if(processing_map[np_px][np_py] == LINEPOINT)
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
						
						if(line_a == line)
						{
							System.err.println("line_a == line!");
							// TODO: remove junction point from line?
						}
						else
						{
							if(PREVENT_SPLINTER_FRAGMENTS && (jp_pos <= SPLINTER_FRAGMENT_THRESHOLD || (line_a.size() - jp_pos) <= SPLINTER_FRAGMENT_THRESHOLD))
							{
								System.err.println("Not splitting line to avoid fragmentation");
							}
							else
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
								
								// RSLV: add junction point to line (see also next else if statement)
								line_a.addLast(new Point(np_px, np_py, np_dlpx, np_dlpy));
								
								lines.add(line_b);
								
								// link points to second line
								for(Point pb : line_b)
								{
									point_to_line_map[pb.px][pb.py] = line_b;
								}
							}
						}
						
						// mark point as junction
						processing_map[np_px][np_py] = JUNCTION;
						processing_map_ip.set(np_px, np_py, JUNCTION);
						point_to_line_map[np_px][np_py] = null;
						
						// add junction point to current line trace
						if(INCLUDE_JUNCTIONS_IN_LINE)
						{
							if(trace_direction == 0)
							{
								line.addFirst(new Point(np_px, np_py, np_dlpx, np_dlpy));
							}
							else
							{
								line.addLast(new Point(np_px, np_py, np_dlpx, np_dlpy));
							}
						}
						
						// stop tracing
						tracing = false;
					}
					else if(processing_map[np_px][np_py] == JUNCTION)
					{
						System.err.println("  Point is junction, stop tracing");
						
						// add line point to current trace
						if(INCLUDE_JUNCTIONS_IN_LINE)
						{
							if(trace_direction == 0)
							{
								line.addFirst(new Point(np_px, np_py, np_dlpx, np_dlpy));
							}
							else
							{
								line.addLast(new Point(np_px, np_py, np_dlpx, np_dlpy));
							}
						}
						
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
			
			// prevent small fragments
			if(PREVENT_SPLINTER_FRAGMENTS && line.size() <= SPLINTER_FRAGMENT_THRESHOLD)
			{
				System.err.println("Rejecting small line fragment");
				lines.remove(line); // TODO: what if line was split from previous line?
				
				for(Point p : line)
				{
					point_to_line_map[p.px][p.py] = null;
					processing_map[p.px][p.py] = LIMITED;
				}
			}
		}
		
		// ------
		
		if(DEBUG_MODE_ENABLED)
		{
			ImagePlus processing_map_imp = new ImagePlus("Processing map", processing_map_ip);
			processing_map_imp.setDisplayRange(0,4); // NOTE: keep up to date with number of processing types
			processing_map_imp.show();
		}
		
		// *********************************************************************
		
		// RSLV: connect lines at/near junctions [5x5 search window]
		
		// *********************************************************************
		
		// TODO: connect lines at a slightly more global scale
		// for each line
		//   check both end point
		//      extend window into direction of line (say 15 pixels long, 10 pixels wide)
		//      if another lines endpoint is within reach
		//         connect these two lines together
		//         option: smooth out connection (i.e. drop first couple of points from endpoint)
		if(ENABLE_LINKING)
		{
			double SEARCH_DISTANCE = 2*Math.ceil(6*SIGMA); // pixels
			int BACKTRACK_DISTANCE = 3; // number of point to backtrack
			int AVERAGING_DISTANCE = 10;
			ImageProcessor search_window_overlay_ip = ip.duplicate();
			Overlay search_window_overlay = new Overlay(); // TMP: DEBUG
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
				
				int wx2 = l.get(l.size()-BACKTRACK_DISTANCE-1).px;
				int wy2 = l.get(l.size()-BACKTRACK_DISTANCE-1).py;
				double wsx2 = l.get(l.size()-BACKTRACK_DISTANCE-1).sx;
				double wsy2 = l.get(l.size()-BACKTRACK_DISTANCE-1).sy;
				int wdx2 = l.get(l.size()-BACKTRACK_DISTANCE-1).px - l.get(l.size()-BACKTRACK_DISTANCE-2).px;
				int wdy2 = l.get(l.size()-BACKTRACK_DISTANCE-1).py - l.get(l.size()-BACKTRACK_DISTANCE-2).py;
				double wdsx2 = (l.get(l.size()-BACKTRACK_DISTANCE-1).px + l.get(l.size()-BACKTRACK_DISTANCE-1).sx) - (l.get(l.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).px + l.get(l.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).sx);
				double wdsy2 = (l.get(l.size()-BACKTRACK_DISTANCE-1).py + l.get(l.size()-BACKTRACK_DISTANCE-1).sy) - (l.get(l.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).py + l.get(l.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).sy);
				double wsa2 = Math.atan2(wdsy2, wdsx2);
				
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
				double MAX_BENDING_ANGLE = 30.0;
				double persistence_length = MAX_BENDING_ANGLE * (Math.PI / 180.0);
				
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
		
//			if(DEBUG_MODE_ENABLED)
//			{
				ImagePlus search_window_overlay_imp = new ImagePlus("Search windows", search_window_overlay_ip);
				search_window_overlay_imp.setOverlay(search_window_overlay); // TMP: DEBUG
				//search_window_overlay_imp.updateAndRepaintWindow();
				search_window_overlay_imp.show();
//			}
			
/*			//HashMap<Tuple<Line, Integer>, HashSet<Tuple<Line, Integer> > > candidate_matches
			for(Map.Entry<Tuple<Line, Integer>, HashSet<Tuple<Line, Integer> > > entry : candidate_matches.entrySet())
			{
				// get source
				Tuple<Line, Integer> source = entry.getKey();
				
				// get details on source
				Line source_line = source.first;
				int source_wx, source_wy;
				double source_wsx, source_wsy;
				int source_wdx, source_wdy;
				double source_wdsx, source_wdsy;
				double source_wsa;
				if(source.second == 1) // front
				{
					source_wx = source_line.get(0+BACKTRACK_DISTANCE).px;
					source_wy = source_line.get(0+BACKTRACK_DISTANCE).py;
					source_wsx = source_line.get(0+BACKTRACK_DISTANCE).sx;
					source_wsy = source_line.get(0+BACKTRACK_DISTANCE).sx;
					source_wdx = source_line.get(0+BACKTRACK_DISTANCE).px - source_line.get(1+BACKTRACK_DISTANCE).px;
					source_wdy = source_line.get(0+BACKTRACK_DISTANCE).py - source_line.get(1+BACKTRACK_DISTANCE).py;
					source_wdsx = (source_line.get(0+BACKTRACK_DISTANCE).px + source_line.get(0+BACKTRACK_DISTANCE).sx) - (source_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).px + source_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).sx);
					source_wdsy = (source_line.get(0+BACKTRACK_DISTANCE).py + source_line.get(0+BACKTRACK_DISTANCE).sy) - (source_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).py + source_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).sy);
					source_wsa = Math.atan2(source_wdsy, source_wdsx);
				}
				else if(source.second == 2) // back
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
				
				// find best match
				Line best_matching_line = null;
				int best_matching_window = -1;
				double best_matching_cost = Double.MAX_VALUE;
				for(Tuple<Line, Integer> destination : destinations)
				{
					// get details on source
					Line destination_line = destination.first;
					int destination_wx, destination_wy;
					double destination_wsx, destination_wsy;
					int destination_wdx, destination_wdy;
					double destination_wdsx, destination_wdsy;
					double destination_wsa;
					if(destination.second == 1) // front
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
					else if(destination.second == 2) // back
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
					// RSLV: redundant with triangular field of views?
					//if(diff_angle > MAX_BENDING_ANGLE)
					//{
					//	continue; // skip match
					//}
					
					double matching_cost = diff_angle; //diff_distance+5*3*SIGMA*diff_angle;//diff_distance+COST_FUNCTION_WEIGHT*diff_angle; // TODO: define cost function based on distance and orientation
					
					// check for better match
					if(matching_cost < best_matching_cost)
					{
						best_matching_cost = matching_cost;
						best_matching_line = destination_line;
						best_matching_window = destination.second;
					}
				}
				
				// if match found
				if(best_matching_line != null && best_matching_window != -1)
				{
					// TODO: connect two lines together
					
					// TODO: remove key destination:window
					// NOTE: iterator invalidated?!
					
					// TODO: replace all instances of destination:window with source:window
					
					// remove from candidate list (both source and destination:source?)
				}
			}
*/
			
/*  LINKING OF LINES 
				// find matching lines
				Vector<Line> lmatches = new Vector<Line>();
				Vector<Integer> wmatches = new Vector<Integer>();
				Vector<Integer> pmatches = new Vector<Integer>();
				
				// try to match to all lines
				for(int j = 0; j < lines.size(); ++j) // RSLV: i+1 may be safer option?
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
					if(olx1 >= wtlx1 && olx1 <= wtlx1+SEARCH_DISTANCE && oly1 >= wtly1 && oly1 <= wtly1+SEARCH_DISTANCE)
					{
						System.err.println("found match between W1 and OL1");
						lmatches.add(ol);
						wmatches.add(new Integer(1));
						pmatches.add(new Integer(1));
					}
					if(olx1 >= wtlx2 && olx1 <= wtlx2+SEARCH_DISTANCE && oly1 >= wtly2 && oly1 <= wtly2+SEARCH_DISTANCE)
					{
						System.err.println("found match between W2 and OL1");
						lmatches.add(ol);
						wmatches.add(new Integer(2));
						pmatches.add(new Integer(1));
					}
					if(olx2 >= wtlx1 && olx2 <= wtlx1+SEARCH_DISTANCE && oly2 >= wtly1 && oly2 <= wtly1+SEARCH_DISTANCE)
					{
						System.err.println("found match between W1 and OL2");
						lmatches.add(ol);
						wmatches.add(new Integer(1));
						pmatches.add(new Integer(2));
					}
					if(olx2 >= wtlx2 && olx2 <= wtlx2+SEARCH_DISTANCE && oly2 >= wtly2 && oly2 <= wtly2+SEARCH_DISTANCE)
					{
						System.err.println("found match between W2 and OL2");
						lmatches.add(ol);
						wmatches.add(new Integer(2));
						pmatches.add(new Integer(2));
					}
				}
				
				// find best matching line
				double best_matching_cost = Double.MAX_VALUE;
				int best_matching_index = -1;
				for(int k = 0; k < lmatches.size(); ++k)
				{
					// get details on match
					Line best_line = lmatches.get(k);
					int best_window = wmatches.get(k);
					int best_point = pmatches.get(k);
					
					// calculate matching cost
					double diff_distance = Double.MAX_VALUE;
					double diff_angle = Math.PI;
					
//				// TMP: copy
//				int wx1 = l.getFirst().px;
//				int wy1 = l.getFirst().py;
//				double wsx1 = l.getFirst().sx;
//				double wsy1 = l.getFirst().sx;
//				int wdx1 = l.get(1).px - l.get(2).px;
//				int wdy1 = l.get(1).py - l.get(2).py;
//				double wdsx1 = (l.get(1).px + l.get(1).sx) - (l.get(2).px + l.get(2).sx);
//				double wdsy1 = (l.get(1).py + l.get(1).sy) - (l.get(2).py + l.get(2).sy);
//				double wsa1 = Math.atan2(wdsy1, wdsx1);
					
					if(best_window == 1 && best_point == 1)
					{
						// calculate distance difference
						double dspx = (wx1 + wsx1) - (best_line.get(0+BACKTRACK_DISTANCE).px + best_line.get(0+BACKTRACK_DISTANCE).sx);
						double dspy = (wy1 + wsy1) - (best_line.get(0+BACKTRACK_DISTANCE).py + best_line.get(0+BACKTRACK_DISTANCE).sy);
						
						diff_distance = Math.sqrt(dspx*dspx+dspy*dspy);
						
						// calculate angle difference
						double ddsx = (best_line.get(0+BACKTRACK_DISTANCE).px + best_line.get(0+BACKTRACK_DISTANCE).sx) - (best_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).px + best_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).sx);
						double ddsy = (best_line.get(0+BACKTRACK_DISTANCE).py + best_line.get(0+BACKTRACK_DISTANCE).sy) - (best_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).py + best_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).sy);
						double dwsa = Math.atan2(ddsy, ddsx);
						
						System.err.println("wsa1 = " + wsa1*180/Math.PI);
						System.err.println("dwsa = " + dwsa*180/Math.PI);
						diff_angle = Math.abs(wsa1 - dwsa);
						System.err.println("diff_angle = " + diff_angle*180/Math.PI);
						if(diff_angle > Math.PI)
						{
							diff_angle = (2*Math.PI) - diff_angle;
							System.err.println("new diff_angle = " + diff_angle*180/Math.PI);
						}
						diff_angle = Math.PI - diff_angle; // flip
					}
					else if(best_window == 1 && best_point == 2)
					{
						// calculate distance difference
						double dspx = (wx1 + wsx1) - (best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).px + best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).sx);
						double dspy = (wy1 + wsy1) - (best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).py + best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).sy);
						
						diff_distance = Math.sqrt(dspx*dspx+dspy*dspy);
						
						// calculate angle difference
						double ddsx = (best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).px + best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).sx) - (best_line.get(best_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).px + best_line.get(best_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).sx);
						double ddsy = (best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).py + best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).sy) - (best_line.get(best_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).py + best_line.get(best_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).sy);
						double dwsa = Math.atan2(ddsy, ddsx);
						
						System.err.println("wsa1 = " + wsa1*180/Math.PI);
						System.err.println("dwsa = " + dwsa*180/Math.PI);
						diff_angle = Math.abs(wsa1 - dwsa);
						System.err.println("diff_angle = " + diff_angle*180/Math.PI);
						if(diff_angle > Math.PI)
						{
							diff_angle = (2*Math.PI) - diff_angle;
							System.err.println("new diff_angle = " + diff_angle*180/Math.PI);
						}
						diff_angle = Math.PI - diff_angle; // flip
					}
					else if(best_window == 2 && best_point == 1)
					{
						// calculate distance difference
						double dspx = (wx2 + wsx2) - (best_line.get(0+BACKTRACK_DISTANCE).px + best_line.get(0+BACKTRACK_DISTANCE).sx);
						double dspy = (wy2 + wsy2) - (best_line.get(0+BACKTRACK_DISTANCE).py + best_line.get(0+BACKTRACK_DISTANCE).sy);
						
						diff_distance = Math.sqrt(dspx*dspx+dspy*dspy);
						
						// calculate angle difference
						double ddsx = (best_line.get(0+BACKTRACK_DISTANCE).px + best_line.get(0+BACKTRACK_DISTANCE).sx) - (best_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).px + best_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).sx);
						double ddsy = (best_line.get(0+BACKTRACK_DISTANCE).py + best_line.get(0+BACKTRACK_DISTANCE).sy) - (best_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).py + best_line.get(1+BACKTRACK_DISTANCE+AVERAGING_DISTANCE).sy);
						double dwsa = Math.atan2(ddsy, ddsx);
						
						System.err.println("wsa2 = " + wsa2*180/Math.PI);
						System.err.println("dwsa = " + dwsa*180/Math.PI);
						diff_angle = Math.abs(wsa2 - dwsa);
						System.err.println("diff_angle = " + diff_angle*180/Math.PI);
						if(diff_angle > Math.PI)
						{
							diff_angle = (2*Math.PI) - diff_angle;
							System.err.println("new diff_angle = " + diff_angle*180/Math.PI);
						}
						diff_angle = Math.PI - diff_angle; // flip
					}
					else if(best_window == 2 && best_point == 2)
					{
						// calculate distance difference
						double dspx = (wx2 + wsx2) - (best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).px + best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).sx);
						double dspy = (wy2 + wsy2) - (best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).py + best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).sy);
						
						diff_distance = Math.sqrt(dspx*dspx+dspy*dspy);
						
						// calculate angle difference
						double ddsx = (best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).px + best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).sx) - (best_line.get(best_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).px + best_line.get(best_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).sx);
						double ddsy = (best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).py + best_line.get(best_line.size()-BACKTRACK_DISTANCE-1).sy) - (best_line.get(best_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).py + best_line.get(best_line.size()-BACKTRACK_DISTANCE-AVERAGING_DISTANCE-2).sy);
						double dwsa = Math.atan2(ddsy, ddsx);
						
						System.err.println("wsa2 = " + wsa2*180/Math.PI);
						System.err.println("dwsa = " + dwsa*180/Math.PI);
						diff_angle = Math.abs(wsa2 - dwsa);
						System.err.println("diff_angle = " + diff_angle*180/Math.PI);
						if(diff_angle > Math.PI)
						{
							diff_angle = (2*Math.PI) - diff_angle;
							System.err.println("new diff_angle = " + diff_angle*180/Math.PI);
						}
						diff_angle = Math.PI - diff_angle; // flip
					}
					
					//if(diff_angle > MAX_BENDING_ANGLE)
					//{
					//	continue; // skip match
					//}
					
					double matching_cost = diff_distance+5*3*SIGMA*diff_angle;//diff_distance+COST_FUNCTION_WEIGHT*diff_angle; // TODO: define cost function based on distance and orientation
					
					// check for better match
					if(matching_cost < best_matching_cost)
					{
						best_matching_cost = matching_cost;
						best_matching_index = k;
					}
				}
				
				// check if match found
				if(best_matching_index != -1)
				{
					// link lines together
					System.err.println("best match found between W" + wmatches.get(best_matching_index) + " and OL" + pmatches.get(best_matching_index));
					
					// get details on match
					Line best_line = lmatches.get(best_matching_index);
					int best_window = wmatches.get(best_matching_index);
					int best_point = pmatches.get(best_matching_index);
					
					// get both lines
					Line l1 = l; // lines.get(i);
					Line l2 = best_line;
					
					// check how to merge
					if(best_window == 1 && best_point == 1)
					{
						for(int k = 0; k < l2.size(); ++k)
						{
							l1.addFirst(l2.get(k)); // RSLV: backtrack?
						}
					}
					else if(best_window == 1 && best_point == 2)
					{
						for(int k = l2.size()-1; k >= 0; --k)
						{
							l1.addFirst(l2.get(k)); // RSLV: backtrack?
						}
					}
					else if(best_window == 2 && best_point == 1)
					{
						for(int k = 0; k < l2.size(); ++k)
						{
							l1.addLast(l2.get(k)); // RSLV: backtrack?
						}
					}
					else if(best_window == 2 && best_point == 2)
					{
						for(int k = l2.size()-1; k >= 0; --k)
						{
							l1.addLast(l2.get(k)); // RSLV: backtrack?
						}
					}
					
					// remove second line from list
					lines.remove(l2);
					
					// set status to updated
					updated = true;
					
					// check again for self (continue merging consequtive line pieces) NOTE: not very neat programming, I know :|
					--i;
				}
			}
		}
		END OF LINKING */
		
		}
		
		// *********************************************************************
		
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
			// skip empty lines
			if(l == null || l.size() == 0) continue;
			
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
			
			roi_manager.addRoi(polyline_s); // RSLV: use add(imp, roi, index)
			
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
		}
		
		// RSLV: return image with overlay?
		roi_manager.setVisible(true);
		roi_manager.toFront();
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
