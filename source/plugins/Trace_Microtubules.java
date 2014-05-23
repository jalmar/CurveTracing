
package plugins;

// import Java classes
//import java.io.File;
import java.awt.Color;

//import java.awt.event.MouseListener; 
import java.awt.event.MouseAdapter; // more convenient class than MouseListener interface
import java.awt.event.MouseEvent;

import java.util.Vector;

// import ImageJ classes
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import ij.plugin.PlugIn;
import ij.gui.GenericDialog;
import ij.Prefs;

import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.Line; // LineRoi
import ij.gui.PolygonRoi; // for POLYLINE

import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import ij.gui.Plot;
import ij.gui.PlotWindow;

import ij.measure.ResultsTable;

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


/**
 *	Microtubule tracing algorithm
 */
public class Trace_Microtubules implements PlugIn
{
	/**
	 *	Constants
	 */
	public static final double DEFAULT_SIGMA = 2.0;
	public static final int MEDIAN_FILTER_SIZE = 3; // direct neighbourhood
	
	public static final String DEFAULT_INTERPOLATION_METHOD_S = "Bilinear"; // None, Nearest Neighbor, Bilinear, Bicubic
	public static String INTERPOLATION_METHOD_S = "Bilinear";
	public static int INTERPOLATION_METHOD = ImageProcessor.BILINEAR; // NONE, NEAREST_NEIGHBOR, BILINEAR, BICUBIC
	public static final int DEFAULT_SAMPLE_RATE = 3;
	public static int SAMPLE_RATE = 3; // NOTE: int for now
	public static double SAMPLE_RATE_INV = 1.0 / (double)SAMPLE_RATE;
	public static double SAMPLE_OFFSET = 0.0;
	
	public static final int DEFAULT_BACKGROUND_LOWER_THRESHOLD = 0;
	public static final int DEFAULT_BACKGROUND_UPPER_THRESHOLD = 65535;
	public static int BACKGROUND_LOWER_THRESHOLD = 0;
	public static int BACKGROUND_UPPER_THRESHOLD = 65535;
	public static final int DEFAULT_AMPLITUDE_LOWER_THRESHOLD = 21845; // 1/3
	public static final int DEFAULT_AMPLITUDE_UPPER_THRESHOLD = 2*65535;
	public static int AMPLITUDE_LOWER_THRESHOLD = 21845; // 1/3
	public static int AMPLITUDE_UPPER_THRESHOLD = 2*65535;
	public static final double DEFAULT_MU_THRESHOLD = 0.707;
	public static double MU_THRESHOLD = 0.707; // sqrt(2)
	public static final double DEFAULT_SIGMA_RATIO_LOWER_THRESHOLD = 0.67;
	public static final double DEFAULT_SIGMA_RATIO_UPPER_THRESHOLD = 1.50;
	public static double SIGMA_RATIO_LOWER_THRESHOLD = 0.67;
	public static double SIGMA_RATIO_UPPER_THRESHOLD = 1.50;
	public static final double DEFAULT_R_SQUARED_THRESHOLD = 0.80;
	public static double R_SQUARED_THRESHOLD = 0.80;
	
	public static boolean FILTER_ON_BACKGROUND = true;
	public static boolean FILTER_ON_AMPLITUDE = true;
	public static boolean FILTER_ON_MU = true;
	public static boolean FILTER_ON_SIGMA = true;
	public static boolean FILTER_ON_R_SQUARED = true;
	
	public static boolean APPLY_FILL_HOLES_ALGORITHM = true;
	public static final int DEFAULT_VOTING_THRESHOLD = 5; // out of eight; RSLV: convert to unity?
	public static int VOTING_THRESHOLD = DEFAULT_VOTING_THRESHOLD;
	
	public static boolean FILTER_ON_COMPONENT_MIN_AREA_SIZE = true;
	public static final int DEFAULT_COMPONENT_MIN_AREA_SIZE_THRESHOLD = 4;
	public static int COMPONENT_MIN_AREA_SIZE_THRESHOLD = DEFAULT_COMPONENT_MIN_AREA_SIZE_THRESHOLD;
	
	public static final boolean DEFAULT_USE_ABSOLUTE_EIGENVALUES = false;
	public static boolean USE_ABSOLUTE_EIGENVALUES = DEFAULT_USE_ABSOLUTE_EIGENVALUES;
	public static final boolean DEFAULT_DEBUG_MODE_ENABLED = false;
	public static boolean DEBUG_MODE_ENABLED = DEFAULT_DEBUG_MODE_ENABLED;
	
	public static final boolean DEFAULT_SHOW_VECTOR_OVERLAY = true;
	public static boolean SHOW_VECTOR_OVERLAY = DEFAULT_SHOW_VECTOR_OVERLAY;
	public static final boolean DEFAULT_SHOW_RESULTS_TABLE = false;
	public static boolean SHOW_RESULTS_TABLE = DEFAULT_SHOW_RESULTS_TABLE;
	
	/**
	 *
	 */
	public Trace_Microtubules()
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
		GenericDialog gd = new GenericDialog("Trace microtubules");
		gd.addNumericField("PSF_sigma", Prefs.get("mt_trace.psf_sigma", DEFAULT_SIGMA), 1);
		gd.addCheckbox("Reduce_noise", Prefs.get("mt_trace.reduce_noise", true));
		gd.addCheckbox("Suppress_background", Prefs.get("mt_trace.suppress_background", true));
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addNumericField("Sample_rate", Prefs.get("mt_trace.sample_rate", DEFAULT_SAMPLE_RATE), 0);
		gd.addChoice("Interpolation", new String[]{"None", "Nearest Neighbor", "Bilinear", "Bicubic"}, Prefs.get("mt_trace.interpolation_method_s", DEFAULT_INTERPOLATION_METHOD_S));
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Filter_on_background", Prefs.get("mt_trace.filter_background", true));
		gd.addNumericField("Background_lower_threshold", Prefs.get("mt_trace.filter_background_lower_threshold", DEFAULT_BACKGROUND_LOWER_THRESHOLD), 0);
		gd.addNumericField("Background_upper_threshold", Prefs.get("mt_trace.filter_background_threshold", DEFAULT_BACKGROUND_UPPER_THRESHOLD), 0);
		
		gd.addCheckbox("Filter_on_amplitude", Prefs.get("mt_trace.filter_amplitude", true));
		gd.addNumericField("Amplitude_lower_threshold", Prefs.get("mt_trace.filter_amplitude_lower_threshold", DEFAULT_AMPLITUDE_LOWER_THRESHOLD), 0);
		gd.addNumericField("Amplitude_upper_threshold", Prefs.get("mt_trace.filter_amplitude_upper_threshold", DEFAULT_AMPLITUDE_UPPER_THRESHOLD), 0);
		
		gd.addCheckbox("Filter_on_mu", Prefs.get("mt_trace.filter_mu", true));
		gd.addNumericField("Mu_absolute_threshold", Prefs.get("mt_trace.filter_mu_threshold", DEFAULT_MU_THRESHOLD), 3);
		
		gd.addCheckbox("Filter_on_sigma", Prefs.get("mt_trace.filter_sigma", true));
		gd.addNumericField("Sigma_ratio_lower_threshold", Prefs.get("mt_trace.filter_sigma_lower_threshold", DEFAULT_SIGMA_RATIO_LOWER_THRESHOLD), 2);
		gd.addNumericField("Sigma_ratio_upper_threshold", Prefs.get("mt_trace.filter_sigma_upper_threshold", DEFAULT_SIGMA_RATIO_UPPER_THRESHOLD), 2);
		
		gd.addCheckbox("Filter_on_R-squared", Prefs.get("mt_trace.filter_r_squared", true));
		gd.addNumericField("R-squared_threshold", Prefs.get("mt_trace.filter_r_squared_threshold", DEFAULT_R_SQUARED_THRESHOLD), 2);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Apply_hole_filling", Prefs.get("mt_trace.fill_holes", true));
		gd.addNumericField("Voting_threshold", Prefs.get("mt_trace.fill_holes_voting_threshold", DEFAULT_VOTING_THRESHOLD), 0);
		
		gd.addCheckbox("Filter_on_component_min_area_size", Prefs.get("mt_trace.filter_component_min_area_size", true));
		gd.addNumericField("Component_min_area_size_threshold", Prefs.get("mt_trace.filter_component_min_area_size_threshold", DEFAULT_COMPONENT_MIN_AREA_SIZE_THRESHOLD), 0);
		
		gd.setInsets(10, 20, 0); // seperate parameter groups
		
		gd.addCheckbox("Use_absolute_eigenvalues", Prefs.get("mt_trace.absolute_eigenvalues", DEFAULT_USE_ABSOLUTE_EIGENVALUES));
		gd.addCheckbox("Enable_debug_mode", Prefs.get("mt_trace.debug_mode", DEFAULT_DEBUG_MODE_ENABLED));
		gd.addCheckbox("Show_vector_overlay", Prefs.get("mt_trace.vector_overlay", DEFAULT_SHOW_VECTOR_OVERLAY));
		gd.addCheckbox("Show_results_table", Prefs.get("mt_trace.results_table", DEFAULT_SHOW_RESULTS_TABLE));
		
		gd.setOKLabel("Trace");
		
		gd.showDialog();
		if(gd.wasCanceled()) return;
		
		// retrieve parameters
		double sigma = gd.getNextNumber();
		boolean reduce_noise = gd.getNextBoolean();
		boolean suppress_background = gd.getNextBoolean();
		
		SAMPLE_RATE = (int)gd.getNextNumber();
		SAMPLE_RATE_INV = 1.0 / SAMPLE_RATE;
		String INTERPOLATION_METHOD_S = gd.getNextChoice();
		if(INTERPOLATION_METHOD_S.equals("None"))
		{
			INTERPOLATION_METHOD = ImageProcessor.NONE;
		}
		else if(INTERPOLATION_METHOD_S.equals("Nearest Neighbor"))
		{
			INTERPOLATION_METHOD = ImageProcessor.NEAREST_NEIGHBOR;
		}
		else if(INTERPOLATION_METHOD_S.equals("Bilinear"))
		{
			INTERPOLATION_METHOD = ImageProcessor.BILINEAR;
		}
		else if(INTERPOLATION_METHOD_S.equals("Bicubic"))
		{
			INTERPOLATION_METHOD = ImageProcessor.BICUBIC;
		}
		
		FILTER_ON_BACKGROUND = gd.getNextBoolean();
		BACKGROUND_LOWER_THRESHOLD = (int)gd.getNextNumber();
		BACKGROUND_UPPER_THRESHOLD = (int)gd.getNextNumber();
		
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
		
		APPLY_FILL_HOLES_ALGORITHM = gd.getNextBoolean();
		VOTING_THRESHOLD = (int)gd.getNextNumber();
		
		FILTER_ON_COMPONENT_MIN_AREA_SIZE = gd.getNextBoolean();
		COMPONENT_MIN_AREA_SIZE_THRESHOLD = (int)gd.getNextNumber();
		
		USE_ABSOLUTE_EIGENVALUES = gd.getNextBoolean();
		DEBUG_MODE_ENABLED = gd.getNextBoolean();
		SHOW_VECTOR_OVERLAY = gd.getNextBoolean();
		SHOW_RESULTS_TABLE = gd.getNextBoolean();
		
		// store parameters in preferences
		Prefs.set("mt_trace.psf_sigma", sigma);
		Prefs.set("mt_trace.reduce_noise", reduce_noise);
		Prefs.set("mt_trace.suppress_background", suppress_background);
		Prefs.set("mt_trace.sample_rate", SAMPLE_RATE);
		Prefs.set("mt_trace.interpolation_method_s", INTERPOLATION_METHOD_S);
		Prefs.set("mt_trace.filter_background", FILTER_ON_BACKGROUND);
		Prefs.set("mt_trace.filter_background_lower_threshold", BACKGROUND_LOWER_THRESHOLD);
		Prefs.set("mt_trace.filter_background_upper_threshold", BACKGROUND_UPPER_THRESHOLD);
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
		Prefs.set("mt_trace.fill_holes", APPLY_FILL_HOLES_ALGORITHM);
		Prefs.set("mt_trace.fill_holes_voting_threshold", VOTING_THRESHOLD);
		Prefs.set("mt_trace.filter_component_min_area_size", FILTER_ON_COMPONENT_MIN_AREA_SIZE);
		Prefs.set("mt_trace.filter_component_min_area_size_threshold", COMPONENT_MIN_AREA_SIZE_THRESHOLD);
		Prefs.set("mt_trace.absolute_eigenvalues", USE_ABSOLUTE_EIGENVALUES);
		Prefs.set("mt_trace.debug_mode", DEBUG_MODE_ENABLED);
		Prefs.set("mt_trace.vector_overlay", SHOW_VECTOR_OVERLAY);
		Prefs.set("mt_trace.results_table", SHOW_RESULTS_TABLE);
		
		// trace microtubules
		run(imp, sigma, reduce_noise, suppress_background);
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
	public static ImagePlus run(ImagePlus imp)
	{
		return run(imp, DEFAULT_SIGMA, true, true);
	}
	
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
	
	public static ImageProcessor run(ImageProcessor ip)
	{
		return run(ip, DEFAULT_SIGMA, true, true);
	}
	
	public static ImageProcessor run(ImageProcessor ip, double sigma, boolean reduce_noise, boolean suppress_background)
	{
		// duplicate image processor
		ImageProcessor ip_original = ip;
		ImageProcessor ip_dup = ip_original.duplicate();
		int original_image_width = ip.getWidth();
		int original_image_height = ip.getHeight();
		
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
		
		// ---------------------------------------------------------------------
		
		// Step 2: white top hat transform to suppress background noise
		ImageProcessor ip_step_2 = ip_step_1.duplicate(); // RSLV: to duplicate or not?
		if(suppress_background)
		{
			Profiling.tic();
			
			IJ.showStatus("Removing background noise");
			ImageProcessor ip_wth = TopHatTransform.run(ip_step_1, (sigma+1)*DEFAULT_SIGMA);
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
		
		// ---------------------------------------------------------------------
		
		// Step 3: eigenvalue/vector decomposition of Hessian matrix
		
		// calculate derivatives of gaussian from image
		Profiling.tic();
		ImageProcessor dx = DerivativeOfGaussian.derivativeX(ip_step_2, sigma);
		ImageProcessor dy = DerivativeOfGaussian.derivativeY(ip_step_2, sigma);
		ImageProcessor dxdx = DerivativeOfGaussian.derivativeXX(ip_step_2, sigma);
		ImageProcessor dxdy = DerivativeOfGaussian.derivativeXY(ip_step_2, sigma);
		ImageProcessor dydx = dxdy;//DerivativeOfGaussian.derivativeYX(ip_step_2, sigma);
		ImageProcessor dydy = DerivativeOfGaussian.derivativeYY(ip_step_2, sigma);
		Profiling.toc("Step 3a: Calculating image derivatives");
		
		// store results [x][y][lambda1=0|lambda2=1|V1x=2|V1y=3|V2x=4|V2y=5]
		double[][][] results_step_3 = new double[original_image_width][original_image_height][6];
		
		// TMP: eigenvalues and eigenvectors in seperate image
		ImageProcessor smallest_eigenvalues_ip = new FloatProcessor(original_image_width, original_image_height);
		ImageProcessor largest_eigenvalues_ip = new FloatProcessor(original_image_width, original_image_height);
		
		ImageProcessor smallest_eigenvectors_ip = new ByteProcessor(original_image_width*5, original_image_height*5);
		ImageProcessor largest_eigenvectors_ip = new ByteProcessor(original_image_width*5, original_image_height*5);
		
		ImageProcessor smallest_theta_ip = new FloatProcessor(original_image_width, original_image_height);
		ImageProcessor largest_theta_ip = new FloatProcessor(original_image_width, original_image_height);
		
		// calculate line points from eigenvalues and eigenvectors based on Hessian matrix
		Profiling.tic();
		for(int py = 0; py < original_image_height; ++py)
		{
			for(int px = 0; px < original_image_width; ++px)
			{
				// construct Hessian matrix in Jama
				Matrix h = new Matrix(2, 2, 0); // 2x2 matrix with zeros
				h.set(0, 0, dxdx.getf(px, py));
				h.set(0, 1, dxdy.getf(px, py));
				h.set(1, 0, dydx.getf(px, py));
				h.set(1, 1, dydy.getf(px, py));
				//System.err.println("Cond="+h.cond());
				
				// compute eigenvalues and eigenvectors
				EigenvalueDecomposition evd = h.eig();
				Matrix d = evd.getD();
				Matrix v = evd.getV();
				
				// determine largest [absolute] (perpendicular -> n(t)) and smallest [absolute] (parallel -> s(t)) eigenvalue and corresponding eigenvector
				double d_00 = d.get(0,0);
				double d_11 = d.get(1,1);
				if(USE_ABSOLUTE_EIGENVALUES)
				{
					d_00 = Math.abs(d_00);
					d_11 = Math.abs(d_11);
				}
				
				// store smallest/largest eigenvalues and corresponding eigenvectors
				if(d_00 >= d_11)
				{
					// d(0,0) is largest (absolute) eigenvalue
					results_step_3[px][py][0] = d_00; //d.get(0,0); // lambda1
					results_step_3[px][py][1] = d_11; //d.get(1,1); // lambda2
					results_step_3[px][py][2] = v.get(0,0); // V1x
					results_step_3[px][py][3] = v.get(1,0); // V1y
					results_step_3[px][py][4] = v.get(0,1); // V2x
					results_step_3[px][py][5] = v.get(1,1); // V2y
				}
				else
				{
					// d(1,1) is largest (absolute) eigenvalue
					results_step_3[px][py][0] = d_11; //d.get(1,1); // lambda1
					results_step_3[px][py][1] = d_00; //d.get(0,0); // lambda2
					results_step_3[px][py][2] = v.get(0,1); // V1x
					results_step_3[px][py][3] = v.get(1,1); // V1y
					results_step_3[px][py][4] = v.get(0,0); // V2x
					results_step_3[px][py][5] = v.get(1,0); // V2y
				}
				
				// DEBUG: store eigenvalues in image processor
				smallest_eigenvalues_ip.setf(px, py, (float)results_step_3[px][py][1]);
				largest_eigenvalues_ip.setf(px, py, (float)results_step_3[px][py][0]);
				
				smallest_theta_ip.setf(px, py, (float)Math.atan2(results_step_3[px][py][4], results_step_3[px][py][5]));
				largest_theta_ip.setf(px, py, (float)Math.atan2(results_step_3[px][py][2], results_step_3[px][py][3]));
				
				// DEBUG: store eigenvectors in image processor
				int cx = px*5+2;
				int cy = py*5+2;
				smallest_eigenvectors_ip.set((int)cx, (int)cy, 255);
				smallest_eigenvectors_ip.set((int)Math.floor(cx-results_step_3[px][py][4]), (int)Math.floor(cy-results_step_3[px][py][5]), 255);
				smallest_eigenvectors_ip.set((int)Math.ceil(cx+results_step_3[px][py][4]), (int)Math.ceil(cy+results_step_3[px][py][5]), 255);
				
				largest_eigenvectors_ip.set((int)cx, (int)cy, 255);
				largest_eigenvectors_ip.set((int)Math.floor(cx-results_step_3[px][py][2]), (int)Math.floor(cy-results_step_3[px][py][3]), 255);
				largest_eigenvectors_ip.set((int)Math.ceil(cx+results_step_3[px][py][2]), (int)Math.ceil(cy+results_step_3[px][py][3]), 255);
				
				/* old code
				double selected_eigenvalue = 0.0; // |n(t)|
				double nx = 0.0; // n(t) -> perpendicular to s(t)
				double ny = 0.0; // n(t) -> perpendicular to s(t)
				
				//if(Math.abs(d.get(0,0)) <= Math.abs(d.get(1,1)))
				if(d.get(0,0) >= d.get(1,1)) // Stegers*
				{
					//selected_eigenvalue = Math.abs(d.get(1,1));
					selected_eigenvalue = d.get(1,1); // Stegers*
					nx = v.get(0,0);
					ny = v.get(1,0);
				}
				else
				{
					//selected_eigenvalue = Math.abs(d.get(0,0));
					selected_eigenvalue = d.get(0,0); // Stegers*
					nx = v.get(0,1);
					ny = v.get(1,1);
				}
				
				// calculate position of line point
				double t = -(dx.getf(px,py)*nx + dy.getf(px,py)*ny)/(dxdx.getf(px,py)*nx*nx + dxdy.getf(px,py)*dydx.getf(px,py)*nx*ny + dydy.getf(px,py)*ny*ny);
				double dlpx = t*nx;
				double dlpy = t*ny;
				
				// store line point
				results_step_3[px][py][0] = dlpx;
				results_step_3[px][py][1] = dlpy;
				results_step_3[px][py][2] = selected_eigenvalue;
				*/
				
				// TEMP: directional image
				//direction_ip.setf(px, py, (float)Math.atan2(ny, nx));
			}
		}
		Profiling.toc("Step 3b: Calculating Hessian eigenvalues and eigenvectors");
		
		// DEBUG: show intermediate images
		if(DEBUG_MODE_ENABLED)
		{
			// smallest/largest [absolute] eigenvalues
			ImagePlus smallest_eigenvalues_imp = new ImagePlus("DEBUG: smallest" + (USE_ABSOLUTE_EIGENVALUES ? " absolute " : " ") + "eigenvalues", smallest_eigenvalues_ip);
			smallest_eigenvalues_imp.resetDisplayRange();
			smallest_eigenvalues_imp.show();
			
			ImagePlus largest_eigenvalues_imp = new ImagePlus("DEBUG: largest" + (USE_ABSOLUTE_EIGENVALUES ? " absolute " : " ") + "eigenvalues", largest_eigenvalues_ip);
			largest_eigenvalues_imp.resetDisplayRange();
			largest_eigenvalues_imp.show();
			
			// smallest/largest [absolute] eigenvectors
			ImagePlus smallest_eigenvectors_imp = new ImagePlus("DEBUG: smallest eigenvectors", smallest_eigenvectors_ip);
			smallest_eigenvectors_imp.resetDisplayRange();
			smallest_eigenvectors_imp.show();
			
			ImagePlus largest_eigenvectors_imp = new ImagePlus("DEBUG: largest eigenvectors", largest_eigenvectors_ip);
			largest_eigenvectors_imp.resetDisplayRange();
			largest_eigenvectors_imp.show();
			
			// theta direction
			ImagePlus smallest_theta_imp = new ImagePlus("DEBUG: theta of smallest eigenvector", smallest_theta_ip);
			smallest_theta_imp.resetDisplayRange();
			smallest_theta_imp.show();
			
			ImagePlus largest_theta_imp = new ImagePlus("DEBUG: theta of largest eigenvector", largest_theta_ip);
			largest_theta_imp.resetDisplayRange();
			largest_theta_imp.show();
		}
		
		// TEMP: directional image
		//ImagePlus direction_img = new ImagePlus("Directional image", direction_ip);
		//direction_img.resetDisplayRange();
		//direction_img.show();
		
		/* using seperate algorithm */
//		ImageProcessor ip_step_3_magnitude = Stegers.run(ip_step_2, sigma);
//		ip_step_3_magnitude.invert();
//		
//		if(DEBUG_MODE_ENABLED)
//		{
//			//ip_step_3.resetMinAndMax();
//			ImagePlus debug_imp = new ImagePlus("DEBUG: Steger's algorithm magnitude image", ip_step_3_magnitude);
//			debug_imp.resetDisplayRange();
//			debug_imp.show();
//		}
//		// TODO: get more than just the magnitude image from Steger's algorithm
		
		// ---------------------------------------------------------------------
		
		// Step 4: perform line profile fitting to data
		Profiling.tic();
		//ImageProcessor fitting_ip = ip_original;
		//if(FIT_ON_RAW_IMAGE) // TODO: implement choice
		double[][][] fitting_results = new double[original_image_width][original_image_height][4]; // [x][y][bg=0|amp=1|mu=2|sigma=3]
		double[][][] standard_error_fit_results = new double[original_image_width][original_image_height][4]; // [x][y][bg=0|amp=1|mu=2|sigma=3]
		double[][] chi_squared_fit_results = new double[original_image_width][original_image_height];
		double[][] r_squared_fit_results = new double[original_image_width][original_image_height];
		for(int py = 0; py < original_image_height; ++py)
		{
			for(int px = 0; px < original_image_width; ++px)
			{
				// get center pixel and vector orientation from (previous) step 3
				double cx = px; // RSLV: use Steger's estimation?
				double cy = py; // RSLV: use Steger's estimation?
				double nx = results_step_3[px][py][4];
				double ny = results_step_3[px][py][5];
				//double theta = atan(ny, nx); // RSLV: add rotation as parameter for more optimal fit than Hessian can produce?
				
				// extract line profile data from *original* image
				ip.setInterpolationMethod(INTERPOLATION_METHOD);
				int line_profile_width = (int)Math.ceil(3*sigma); // RSLV: 3*sigma minimum?
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
					double ix = cx + ii * nx + SAMPLE_OFFSET; // RSLV: no half pex offset?
					double iy = cy + ii * ny + SAMPLE_OFFSET; // RSLV: no half pex offset?
					double pv = ip.getPixelInterpolated(ix, iy); //or ip.getInterpolatedPixel(ix, iy), or ip.getInterpolatedValue(ix, iy);
					//double pv = ip.getPixelValue((int)ix, (int)iy); // no interpolation
					//double pv = ip.getPixelValue((int)Math.round(ix), (int)Math.round(iy)); // nearest neighbor interpolation
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
				//double sigma = sigma;
				double[] initial_parameters = new double[]{background, amplitude, mu, sigma};
				
				// set up new LMA instance
				LevenbergMarquardt lma = new LevenbergMarquardt();
				
				// run LMA fitting procedure
				double[] fitted_parameters = lma.run(x_data, y_data, data_points, initial_parameters, 10, 0.001);
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
		Profiling.toc("Step 3c: Fitting line profiles");
		
		// DEBUG: show intermediate images
		if(DEBUG_MODE_ENABLED)
		{
			// generate images
			ImageProcessor background_fit_ip = new FloatProcessor(original_image_width, original_image_height);
			ImageProcessor amplitude_fit_ip = new FloatProcessor(original_image_width, original_image_height);
			ImageProcessor mu_fit_ip = new FloatProcessor(original_image_width, original_image_height);
			ImageProcessor mu_squared_fit_ip = new FloatProcessor(original_image_width, original_image_height);
			ImageProcessor sigma_fit_ip = new FloatProcessor(original_image_width, original_image_height);
			ImageProcessor chi_squared_fit_ip = new FloatProcessor(original_image_width, original_image_height);
			ImageProcessor log_chi_squared_fit_ip = new FloatProcessor(original_image_width, original_image_height);
			ImageProcessor r_squared_fit_ip = new FloatProcessor(original_image_width, original_image_height);
			for(int py = 0; py < original_image_height; ++py)
			{
				for( int px = 0; px < original_image_width; ++px)
				{
					background_fit_ip.setf(px, py, (float)fitting_results[px][py][0]);
					amplitude_fit_ip.setf(px, py, (float)fitting_results[px][py][1]);
					mu_fit_ip.setf(px, py, (float)fitting_results[px][py][2]);
					mu_squared_fit_ip.setf(px, py, (float)(fitting_results[px][py][2]*fitting_results[px][py][2]));
					sigma_fit_ip.setf(px, py, (float)Math.abs(fitting_results[px][py][3]));
					chi_squared_fit_ip.setf(px, py, (float)chi_squared_fit_results[px][py]);
					log_chi_squared_fit_ip.setf(px, py, (float)Math.log(chi_squared_fit_results[px][py]));
					r_squared_fit_ip.setf(px, py, (float)Math.log(r_squared_fit_results[px][py]));
				}
			}
			
			// smallest/largest [absolute] eigenvalues
			ImagePlus background_fit_imp = new ImagePlus("DEBUG: background fit", background_fit_ip);
			background_fit_imp.resetDisplayRange();
			background_fit_imp.show();
			
			ImagePlus amplitude_fit_imp = new ImagePlus("DEBUG: amplitude fit", amplitude_fit_ip);
			amplitude_fit_imp.resetDisplayRange();
			amplitude_fit_imp.show();
			
//			ImagePlus mu_fit_imp = new ImagePlus("DEBUG: mu fit", mu_fit_ip);
//			mu_fit_imp.resetDisplayRange();
//			mu_fit_imp.show();
			
			ImagePlus mu_squared_fit_imp = new ImagePlus("DEBUG: mu squared fit", mu_squared_fit_ip);
			mu_squared_fit_imp.resetDisplayRange();
			mu_squared_fit_imp.setDisplayRange(0, sigma);
			mu_squared_fit_imp.show();
			
			ImagePlus sigma_fit_imp = new ImagePlus("DEBUG: sigma fit", sigma_fit_ip);
			//sigma_fit_imp.resetDisplayRange();
			sigma_fit_imp.setDisplayRange(0, 5*sigma);
			sigma_fit_imp.show();
			
//			ImagePlus chi_fit_imp = new ImagePlus("DEBUG: chi sqaured fit", chi_fit_ip);
//			chi_fit_imp.resetDisplayRange();
//			chi_fit_imp.show();
			
			ImagePlus log_chi_squared_fit_imp = new ImagePlus("DEBUG: log chi squared fit", log_chi_squared_fit_ip);
			log_chi_squared_fit_imp.resetDisplayRange();
			log_chi_squared_fit_imp.show();
			
			ImagePlus r_squared_fit_imp = new ImagePlus("DEBUG: r squared fit", r_squared_fit_ip);
			r_squared_fit_imp.resetDisplayRange();
			r_squared_fit_imp.show();
		}
		
		// *********************************************************************
		
		// Step 5: additional filtering of point (determine which pixels are microtubules and which belong to the background
		double[][] filtered_pixels = new double[original_image_width][original_image_height];
		ImageProcessor hit_ip = new ByteProcessor(original_image_width, original_image_height);
		ImageProcessor hit_count_ip = new ByteProcessor(original_image_width, original_image_height);
		for(int py = 0; py < original_image_height; ++py)
		{
			for(int px = 0; px < original_image_width; ++px)
			{
				// filter on background
				boolean filtered = false;
				if(FILTER_ON_BACKGROUND && (fitting_results[px][py][0] < BACKGROUND_LOWER_THRESHOLD || fitting_results[px][py][0] > BACKGROUND_UPPER_THRESHOLD))
				{
					// probably not a valid center line pixel; skip
					filtered = true; //continue;
				}
				
				
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
				if(FILTER_ON_SIGMA && (fitting_results[px][py][3] < SIGMA_RATIO_LOWER_THRESHOLD * sigma || fitting_results[px][py][3] > SIGMA_RATIO_UPPER_THRESHOLD * sigma))
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
					hit_ip.set(px, py, 255);
					hit_count_ip.putPixel((int)(px+0.5+fitting_results[px][py][2]*results_step_3[px][py][4]), (int)(py+0.5+fitting_results[px][py][2]*results_step_3[px][py][5]), hit_count_ip.getPixel((int)(px+0.5+fitting_results[px][py][2]*results_step_3[px][py][4]), (int)(py+0.5+fitting_results[px][py][2]*results_step_3[px][py][5])) + 1); // NOTE: include mu correction
				}
			}
		}
		
		// Step 5b: relax filtering criteria based on neighbourhood
		// using an iterative binary voting algorithm to test for membership
		ImageProcessor hit_filled_ip = hit_ip.duplicate();
		ImageProcessor hit_filled_count_ip = hit_count_ip.duplicate();
		if(APPLY_FILL_HOLES_ALGORITHM)
		{
			// run fill holes algorithm
			hit_filled_ip = FillHoles.run(hit_filled_ip, VOTING_THRESHOLD);
			
			// create relaxed hit count image; RSLV: find a way to avoid second pass?
			hit_filled_count_ip = new ByteProcessor(original_image_width, original_image_height); // clear image (fill with zeros)
			for(int py = 0; py < original_image_height; ++py)
			{
				for(int px = 0; px < original_image_width; ++px)
				{
					if(hit_filled_ip.get(px, py) > 0)
					{
						hit_filled_count_ip.putPixel((int)(px+0.5+fitting_results[px][py][2]*results_step_3[px][py][4]), (int)(py+0.5+fitting_results[px][py][2]*results_step_3[px][py][5]), hit_filled_count_ip.getPixel((int)(px+0.5+fitting_results[px][py][2]*results_step_3[px][py][4]), (int)(py+0.5+fitting_results[px][py][2]*results_step_3[px][py][5])) + 1); // NOTE: include mu correction
					}
				}
			}
		}
		
		// Step 5c: filter components on minimum area size
		ImageProcessor hit_filled_filtered_ip = hit_filled_ip.duplicate();
		ImageProcessor hit_filled_filtered_count_ip = hit_filled_count_ip.duplicate();
		if(FILTER_ON_COMPONENT_MIN_AREA_SIZE)
		{
			hit_filled_filtered_ip = ConnectedComponents.run(hit_filled_filtered_ip, ConnectedComponents.Connectivity.EIGHT_CONNECTIVITY, COMPONENT_MIN_AREA_SIZE_THRESHOLD, Integer.MAX_VALUE);
			
			// create relaxed hit count image; RSLV: find a way to avoid second pass!!
			hit_filled_filtered_count_ip = new ByteProcessor(original_image_width, original_image_height); // clear image (fill with zeros)
			for(int py = 0; py < original_image_height; ++py)
			{
				for(int px = 0; px < original_image_width; ++px)
				{
					if(hit_filled_filtered_ip.get(px, py) > 0)
					{
						hit_filled_filtered_count_ip.putPixel((int)(px+0.5+fitting_results[px][py][2]*results_step_3[px][py][4]), (int)(py+0.5+fitting_results[px][py][2]*results_step_3[px][py][5]), hit_filled_filtered_count_ip.getPixel((int)(px+0.5+fitting_results[px][py][2]*results_step_3[px][py][4]), (int)(py+0.5+fitting_results[px][py][2]*results_step_3[px][py][5])) + 1); // NOTE: include mu correction
					}
				}
			}
		}
		
		// DEBUG: show intermediate images
//		if(DEBUG_MODE_ENABLED)
//		{
			// show hit images
			ImagePlus hit_imp = new ImagePlus("DEBUG: hit image", hit_ip);
			ImagePlus hit_count_imp = new ImagePlus("DEBUG: hit count image", hit_count_ip);
			hit_imp.resetDisplayRange();
			hit_count_imp.resetDisplayRange();
			hit_imp.show();
			hit_count_imp.show();
			
			// show relaxed hit images
			if(APPLY_FILL_HOLES_ALGORITHM)
			{
				ImagePlus hit_filled_imp = new ImagePlus("DEBUG: hit filled image", hit_filled_ip);
				ImagePlus hit_filled_count_imp = new ImagePlus("DEBUG: hit filled count image", hit_filled_count_ip);
				hit_filled_imp.resetDisplayRange();
				hit_filled_count_imp.resetDisplayRange();
				hit_filled_imp.show();
				hit_filled_count_imp.show();
			}
			
			// show hit images filtered on minimum component size
			if(FILTER_ON_COMPONENT_MIN_AREA_SIZE)
			{
				ImagePlus hit_filled_filtered_imp = new ImagePlus("DEBUG: hit filled filtered image", hit_filled_filtered_ip);
				//ImagePlus hit_filled_filtered_count_imp = new ImagePlus("DEBUG: hit filled filtered count image", hit_filled_filtered_count_ip);
				hit_filled_filtered_imp.resetDisplayRange();
				//hit_filled_filtered_count_imp.resetDisplayRange();
				hit_filled_filtered_imp.show();
				//hit_filled_filtered_count_imp.show();
			}
//		}
		
		// *********************************************************************
		
		// DEBUG: show intermediate images
		if(SHOW_VECTOR_OVERLAY)
		{
			// TMP: create scaled image for vector overlay
			Profiling.tic();
			final int SCALE_FACTOR = 1;
			ip_original.setInterpolationMethod(ImageProcessor.NONE); // NONE, NEAREST_NEIGHBOR, BILINEAR, BICUBIC
			ImageProcessor scaled_ip = ip_original.resize(original_image_width*SCALE_FACTOR, original_image_height*SCALE_FACTOR);//duplicate();
			Overlay eigenvectors_overlay = new Overlay();
			
			// create vector overlay of primary [and secondary] eigenvector on top of scaled *original* image
			for(int py = 0; py < original_image_height; ++py)
			{
				for(int px = 0; px < original_image_width; ++px)
				{
					// check filter status of pixel
					boolean filtered = (hit_filled_filtered_ip.get(px, py) == 0); // RSLV: probability; e.g 0.5 for pixels picked up by the hole filling algorithm?
					
					// DEBUG: overlay vector on scaled image
					double cx = px*SCALE_FACTOR+0.5*SCALE_FACTOR;
					double cy = py*SCALE_FACTOR+0.5*SCALE_FACTOR;
					
					Roi largest_eigenvector_roi = new Line(cx-0.4*SCALE_FACTOR*results_step_3[px][py][2], cy-0.4*SCALE_FACTOR*results_step_3[px][py][3], cx+0.4*SCALE_FACTOR*results_step_3[px][py][2], cy+0.4*SCALE_FACTOR*results_step_3[px][py][3]);
					
					largest_eigenvector_roi.setStrokeColor(filtered ? Color.DARK_GRAY : Color.YELLOW);
					largest_eigenvector_roi.setStrokeWidth(0.0);
					largest_eigenvector_roi.setPosition(0); // RSLV: only applicable to single frame image
					eigenvectors_overlay.add(largest_eigenvector_roi);
					
					if(!filtered)
					{
						Roi largest_eigenvector_corrected_roi = new Line(cx-0.4*SCALE_FACTOR*results_step_3[px][py][2]+fitting_results[px][py][2]*results_step_3[px][py][4], cy-0.4*SCALE_FACTOR*results_step_3[px][py][3]+fitting_results[px][py][2]*results_step_3[px][py][5], cx+0.4*SCALE_FACTOR*results_step_3[px][py][2]+fitting_results[px][py][2]*results_step_3[px][py][4], cy+0.4*SCALE_FACTOR*results_step_3[px][py][3]+fitting_results[px][py][2]*results_step_3[px][py][5]);
						
						largest_eigenvector_corrected_roi.setStrokeColor(Color.GREEN); // keep yellow for now
						largest_eigenvector_corrected_roi.setStrokeWidth(0.0);
						largest_eigenvector_corrected_roi.setPosition(0); // RSLV: only applicable to single frame image
						eigenvectors_overlay.add(largest_eigenvector_corrected_roi);
					}
					
					/*Roi smallest_eigenvector_roi = new Line(cx-1.0*results_step_3[px][py][4], cy-1.0*results_step_3[px][py][5], cx+1.0*results_step_3[px][py][4], cy+1.0*results_step_3[px][py][5]);
					
					smallest_eigenvector_roi.setStrokeWidth(0.0);
					smallest_eigenvector_roi.setStrokeColor(Color.orange);
					smallest_eigenvector_roi.setPosition(0); // NOTE: redundant in single frame instance
					eigenvectors_overlay.add(smallest_eigenvector_roi);*/
				}
			}
			
			
			// TEMP: outside debug
			// show scaled image with eigenvectors overlay
			final ImagePlus scaled_imp = new ImagePlus("DEBUG: eigenvectors overlay", scaled_ip); // TMP: pass as model to custom MouseAdapter
			//scaled_imp.resetDisplayRange();
			scaled_imp.setOverlay(eigenvectors_overlay);
			//scaled_imp.updateAndRepaintWindow();
			scaled_imp.show();
			
			// RSLV: if this interaction tool will be used in the release; create subclass of MouseAdapter to accept the model as parameters, rather than using final class from parent class!
			final ImageWindow img_wnd = scaled_imp.getWindow(); // final for inner class access
			final ImageCanvas img_cnv = img_wnd.getCanvas(); // final for inner class access
			final double[][][] hessian_results_tmp = results_step_3; // TMP: pass as model to custom MouseAdapter
			final double[][][] fitting_results_tmp = fitting_results; // TMP: pass as model to custom MouseAdapter
			final double[][][] standard_error_fit_results_tmp = standard_error_fit_results; // TMP: pass as model to custom MouseAdapter
			final double[][] chi_squared_fit_results_tmp = chi_squared_fit_results; // TMP: pass as model to custom MouseAdapter
			final double[][] r_squared_fit_results_tmp = r_squared_fit_results; // TMP: pass as model to custom MouseAdapter
			final double sigma_tmp = sigma; // TMP: pass as model to custom MouseAdapter
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
				private Line perpendicular_line = null;
				
				private PlotWindow plot_wnd = null;
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
						scaled_imp.getOverlay().remove(clicked_pixel);
					}
					clicked_pixel = new Roi(current_x, current_y, 1, 1);
					clicked_pixel.setStrokeColor(Color.BLUE);
					clicked_pixel.setStrokeWidth(0.0);
					scaled_imp.getOverlay().add(clicked_pixel);
					scaled_imp.updateAndRepaintWindow(); // force update of window
					
					// determine action:
					// [1] new pixel selected -> show information
					// [2] same pixel selected -> show trace menu
					if(previous_x == current_x && previous_y == current_y)
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
							
							int segment_count = 0;
							Vector<Double> trace_xs_vec = new Vector<Double>();
							Vector<Double> trace_ys_vec = new Vector<Double>();
							
							// first position and direction
							double tx = current_x + 0.5;
							double ty = current_y + 0.5;
							
							// mu correction
							tx += hessian_results_tmp[(int)tx][(int)ty][4] * fitting_results_tmp[(int)tx][(int)ty][2];
							ty += hessian_results_tmp[(int)tx][(int)ty][5] * fitting_results_tmp[(int)tx][(int)ty][2];
							
							if((int)tx != current_x || (int)ty != current_y)
							{
								tx = (int)tx + 0.5;
								ty = (int)ty + 0.5;
								tx += hessian_results_tmp[(int)tx][(int)ty][4] * fitting_results_tmp[(int)tx][(int)ty][2];
								ty += hessian_results_tmp[(int)tx][(int)ty][5] * fitting_results_tmp[(int)tx][(int)ty][2];
							}
							
							/*
							double prev_tx = tx;
							double prev_ty = ty;
							// search for centerline pixel
							do
							{
								// backup coordinates
								prev_tx = tx;
								prev_ty = ty;
								
								// pixel center
								tx = (int)tx + 0.5;
								ty = (int)tx + 0.5;
								
								
								
							}
							while((int)tx != (int)prev_tx || (int)ty != (int)prev_ty); // WARN: could generate a deadlock!
							*/
							
							trace_xs_vec.add(tx);
							trace_ys_vec.add(ty);
							
							// trace path along centerline pixels
							while(segment_count < MAX_STEPS)
							{
								// backup previous position
								double prev_tx = tx;
								double prev_ty = ty;
								
								// get largest eigenvector of pixel
								double sx = hessian_results_tmp[(int)tx][(int)ty][2]; // TODO: solve ArrayIndexOutOfBounds when running off the image!
								double sy = hessian_results_tmp[(int)tx][(int)ty][3]; // TODO: solve ArrayIndexOutOfBounds when running off the image!
								tx += STEP_SIZE * sx;
								ty += STEP_SIZE * sy;
								
								// mu correction
								tx += hessian_results_tmp[(int)tx][(int)ty][4] * fitting_results_tmp[(int)tx][(int)ty][2];
								ty += hessian_results_tmp[(int)tx][(int)ty][5] * fitting_results_tmp[(int)tx][(int)ty][2];
								
								// repeat if in next pixel
								if((int)tx != (int)prev_tx || (int)ty != (int)prev_ty)
								{
									sx = hessian_results_tmp[(int)tx][(int)ty][2]; // TODO: solve ArrayIndexOutOfBounds when running off the image!
									sy = hessian_results_tmp[(int)tx][(int)ty][3]; // TODO: solve ArrayIndexOutOfBounds when running off the image!
									tx = (int)tx + 0.5; //+= STEP_SIZE * sx;
									ty = (int)ty + 0.5; //+= STEP_SIZE * sy;
									
									// mu correction
									tx += hessian_results_tmp[(int)tx][(int)ty][4] * fitting_results_tmp[(int)tx][(int)ty][2];
									ty += hessian_results_tmp[(int)tx][(int)ty][5] * fitting_results_tmp[(int)tx][(int)ty][2];
								}
								
								trace_xs_vec.add(tx);
								trace_ys_vec.add(ty);
								++segment_count;
							}
							
							// TODO: opposite direction
							//tx = img_cnv.offScreenX(e.getX()) + 0.5;
							//ty = img_cnv.offScreenY(e.getY()) + 0.5;
							//segment = 0;
							
							// manually convert Float collection to array of floats
							float[] trace_xs = new float[segment_count];
							float[] trace_ys = new float[segment_count];
							for(int i = 0; i < segment_count; ++i)
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
								scaled_imp.getOverlay().remove(trace);
							}
							trace = new PolygonRoi(trace_xs, trace_ys, Roi.POLYLINE);
							trace.setStrokeColor(Color.RED);
							trace.setStrokeWidth(0.0);
							scaled_imp.getOverlay().add(trace);
							scaled_imp.updateAndRepaintWindow(); // force update of window
						}
					}
					else
					{
						// display line profile perpendicular to curve
						if(perpendicular_line != null)
						{
							// remove previous trace first
							scaled_imp.getOverlay().remove(perpendicular_line);
						}
						double nx = hessian_results_tmp[current_x][current_y][4];
						double ny = hessian_results_tmp[current_x][current_y][5];
						perpendicular_line = new Line(current_x+0.5-nx*3*sigma_tmp, current_y+0.5-ny*3*sigma_tmp, current_x+0.5+nx*3*sigma_tmp, current_y+0.5+ny*3*sigma_tmp);
						perpendicular_line.setStrokeColor(Color.RED);
						perpendicular_line.setStrokeWidth(0.0);
						scaled_imp.getOverlay().add(perpendicular_line);
						scaled_imp.updateAndRepaintWindow(); // force update of window
						
						// TODO: display information about pixel
						// raw data [intensity, Hessian, fitting]
						// image region and parameters
						
						// populate line profile and Gaussian curve with data (NOTE: adapted from code of fitting procedure)
						ip_tmp.setInterpolationMethod(INTERPOLATION_METHOD); // NONE, NEAREST_NEIGHBOR, BILINEAR, BICUBIC, 
						int line_profile_width = (int)Math.ceil(3*sigma_tmp);
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
							double ix = current_x + ii * nx + SAMPLE_OFFSET; // RSLV: no half pex offset?
							double iy = current_y + ii * ny + SAMPLE_OFFSET; // RSLV: no half pex offset?
							double pv = ip_tmp.getPixelInterpolated(ix, iy); //or ip_tmp.getInterpolatedPixel(ix, iy), or ip_tmp.getInterpolatedValue(ix, iy);
							//double pv = ip_tmp.getPixelValue((int)ix, (int)iy); // no interpolation
							//double pv = ip_tmp.getPixelValue((int)Math.round(ix), (int)Math.round(iy)); // nearest neighbor interpolation
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
						
						profile_plot.setColor(Color.RED);
						profile_plot.addPoints(plot_xs, profile_ys, Plot.CIRCLE); // LINE, DOT, CROSS, CIRCLE, BOX, TRIANGLE
						
						// add labels (reversed order)
						profile_plot.setColor(Color.RED);
						profile_plot.addLabel(0.02, 0.05, " o   Line profile");
						profile_plot.setColor(Color.BLUE);
						profile_plot.addLabel(0.02, 0.08, "---  Gaussian fit\n       bg  = " + String.format("%.2f", fit_param[0]) + "\n       amp = " + String.format("%.2f", fit_param[1]) + "\n       mu  = " + String.format("%.4f", fit_param[2]) + "\n       sig = " + String.format("%.4f", fit_param[3]) + "\n       chi = " + String.format("%.1f", chi_squared_fit_results_tmp[current_x][current_y]) + "\n       log = " + String.format("%.5f", Math.log(chi_squared_fit_results_tmp[current_x][current_y])) + "\n       R^2 = " + String.format("%.4f", r_squared_fit_results_tmp[current_x][current_y]) + "\n\n       ecc = " + String.format("%.2f", Math.sqrt((hessian_results_tmp[current_x][current_y][1]*hessian_results_tmp[current_x][current_y][1]) - (hessian_results_tmp[current_x][current_y][0]*hessian_results_tmp[current_x][current_y][0]))) + "\n       mean = " + String.format("%.4f", mean_value) + "\n       var = " + String.format("%.4f", variance_value)+ "\n       lvar = " + String.format("%.4f", Math.log(variance_value)) + "\n       chi = " + String.format("%.1f", chi_square)+ "\n       chi_log = " + String.format("%.4f", Math.log(chi_square)) + "\n       chi_var = " + String.format("%.4f", chi_square_variance)+ "\n       chi_sum = " + String.format("%.4f", chi_square_sum_square) + "\n       chi_exp = " + String.format("%.4f", chi_square_expected));
						
						// create new plot window or draw in existing plot window
						if(plot_wnd != null && !plot_wnd.isClosed())
						{
							plot_wnd.drawPlot(profile_plot);
							plot_wnd.setTitle("Profile plot @ (x=" + current_x + ", y=" + current_y + ")");
							// RSLV: set focus on plot window?
						}
						else
						{
							plot_wnd = profile_plot.show();
						}
					}
				}
			});
			
			Profiling.toc("DEBUG: drawing eigenvector overlay on scaled image");
		}
		
		// *********************************************************************
		
		// ---------------------------------------------------------------------
		
		
		
		// Future steps: trace individual lines
		
		// Display table with results
		if(SHOW_RESULTS_TABLE)
		{
			Profiling.tic();
			ResultsTable raw_data_table = new ResultsTable();
			raw_data_table.setPrecision(5);
			raw_data_table.showRowNumbers(false);
			//raw_data_table.reset(); // to clear the table
			
			// TODO: fill table
			for(int py = 0; py < original_image_height; ++py)
			{
				for(int px = 0; px < original_image_width; ++px)
				{
					// new row of results
					raw_data_table.incrementCounter();
					
					// pixel coordinate
					raw_data_table.addValue("px", px);
					raw_data_table.addValue("py", py);
					
					// gradients
					raw_data_table.addValue("dx", dx.getf(px, py));
					raw_data_table.addValue("dy", dy.getf(px, py));
					raw_data_table.addValue("dxdx", dxdx.getf(px, py));
					raw_data_table.addValue("dxdy", dxdy.getf(px, py));
					raw_data_table.addValue("dydx", dydx.getf(px, py));
					raw_data_table.addValue("dydy", dydy.getf(px, py));
					
					// Hessian matrix
					raw_data_table.addValue("L1", results_step_3[px][py][0]); // L1
					raw_data_table.addValue("L2", results_step_3[px][py][1]); // L2
					raw_data_table.addValue("V1x", results_step_3[px][py][2]); // V1x
					raw_data_table.addValue("V1y", results_step_3[px][py][3]); // V1y
					raw_data_table.addValue("V2x", results_step_3[px][py][4]); // V2x
					raw_data_table.addValue("V2y", results_step_3[px][py][5]); // V2y
					
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
					raw_data_table.addValue("chi^2", chi_squared_fit_results[px][py]); // chi^2
					raw_data_table.addValue("R^2", r_squared_fit_results[px][py]); // R^2
					
					// filter status
					raw_data_table.addValue("hit", hit_ip.get(px, py)); // hit
					raw_data_table.addValue("hit_count", hit_count_ip.get(px, py)); // hit count
					raw_data_table.addValue("hit_filled", hit_filled_ip.get(px, py)); // hit filled
					raw_data_table.addValue("hit_filled_count", hit_filled_count_ip.get(px, py)); // hit filled count
				}
			}
			
			raw_data_table.show("Result of tracing");
			Profiling.toc("Generating results table");
		}
		
		// Finally, return the result (should eventually be a list of coordinates for each line trace in the image)
		return null;//ip_step_2; // TODO: keep up to date with final processing step
	}
}
