
package algorithms;

// import Java classes
import java.util.Collections;
import java.util.Comparator;
import java.util.Vector;

import java.awt.Color;

// import ImageJ classes
import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.LUT;

import ij.gui.Overlay;
import ij.gui.Roi;
//import ij.gui.Line; // LineRoi // NOTE: will conflict with core.Line
import ij.gui.OvalRoi;
import ij.gui.PolygonRoi; // for POLYLINE

// import Jama classes
import Jama.Matrix;
import Jama.EigenvalueDecomposition;

// import own classes
import filters.DerivativeOfGaussian;
import filters.ConnectedComponents; // for rainbow LUT
import core.Point;
import core.Line;


/**
 *	Steger's curve tracing algorithm
 */
public class Stegers
{
	/**
	 *
	 */
	public static final double DEFAULT_SIGMA = 1.8;
	public static final double DEFAULT_THRESHOLD = 0.5;
	public enum Mode { ABS_MAX, MAX, ABS_MIN, MIN };
	
	public static Mode current_mode = Mode.ABS_MAX;
	
	public static double COST_FUNCTION_WEIGHT = 1.0;
	public static double ROUGHLY_PARALLEL_THRESHOLD = 360.0 / 32; // about 11.25 degree
	
	/**
	 *
	 */
	//public Stegers()
	//{
	//	/* nothing */
	//}
	
	// ////////////////////////////////////////////////////////////////////////
	
	/**
	 *
	 */
	public static void setMode(Mode m)
	{
		current_mode = m;
	}
	
	public static ImageProcessor run(ImageProcessor ip)
	{
		return run(ip, DEFAULT_SIGMA);
	}
	
	public static ImageProcessor run(ImageProcessor ip, double sigma)
	{
		return run(ip, sigma, DEFAULT_THRESHOLD);
	}
	
	public static ImageProcessor run(ImageProcessor ip, double sigma, double threshold)
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
		
		// calculate derivatives of gaussian from image
		ImageProcessor dx = DerivativeOfGaussian.derivativeX(ip, sigma);
		ImageProcessor dy = DerivativeOfGaussian.derivativeY(ip, sigma);
		ImageProcessor dxdx = DerivativeOfGaussian.derivativeXX(ip, sigma);
		ImageProcessor dxdy = DerivativeOfGaussian.derivativeXY(ip, sigma);
		ImageProcessor dydx = dxdy;//DerivativeOfGaussian.derivativeYX(ip, sigma);
		ImageProcessor dydy = DerivativeOfGaussian.derivativeYY(ip, sigma);
		
		// calculate line points from eigenvalues and eigenvectors based on Hessian matrix
		for(int py = 0; py < image_height; ++py)
		{
			for(int px = 0; px < image_width; ++px)
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
				
				// determine largest absolute (perpendicular -> n(t)) and smallest absolute (parallel -> s(t)) eigenvalue and corresponding eigenvector
				double selected_eigenvalue = 0.0; // |n(t)|
				double nx = 0.0; // n(t) -> perpendicular to s(t)
				double ny = 0.0; // n(t) -> perpendicular to s(t)
				double second_eigenvalue = 0.0;
				double sx = 0.0;
				double sy = 0.0;
				if((current_mode == Mode.ABS_MAX && Math.abs(d.get(0,0)) >= Math.abs(d.get(1,1))) // Stegers*: absolute maximum
				|| (current_mode == Mode.MAX && d.get(0,0) >= d.get(1,1)) // real maximum
				|| (current_mode == Mode.ABS_MIN && Math.abs(d.get(0,0)) <= Math.abs(d.get(1,1))) // absolute minimum
				|| (current_mode == Mode.MIN && d.get(0,0) <= d.get(1,1))) // real minimum
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
				
				// calculate position of line point and filter line points
				double t = ((dx.getf(px,py)*nx + dy.getf(px,py)*ny) / (dxdx.getf(px,py)*nx*nx + dxdy.getf(px,py)*nx*ny + dydx.getf(px,py)*ny*nx + dydy.getf(px,py)*ny*ny)); // NOTE: removed '-' from start of equation!!
				double dlpx = t*nx;
				double dlpy = t*ny;
				
				// store line point information
				line_points[px][py][0] = selected_eigenvalue;
				line_points[px][py][1] = nx;
				line_points[px][py][2] = ny;
				line_points[px][py][3] = second_eigenvalue;
				line_points[px][py][4] = sx;
				line_points[px][py][5] = sy;
				line_points[px][py][6] = dlpx;
				line_points[px][py][7] = dlpy;
			}
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
				if(Math.abs(line_points[px][py][6]) <= threshold && Math.abs(line_points[px][py][7]) <= threshold)
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
		
		// TEMP: debug images
		ImagePlus valid_line_points_mask_image = new ImagePlus("Valid line points mask", valid_line_points_mask_ip);
		valid_line_points_mask_image.resetDisplayRange();
		valid_line_points_mask_image.show();
		
		// TEMP: debug images
		ImagePlus valid_line_points_magnitude_image = new ImagePlus("Valid line points magnitude", valid_line_points_magnitude_ip);
		valid_line_points_magnitude_image.resetDisplayRange();
		valid_line_points_magnitude_image.show();
		
		// *********************************************************************
		
		// STEP 3: trace lines; requires two thresholds
		
		// get tracing parameters
		double upper_threshold = 5000.0; // RSLV: how to get user threshold?
		double lower_threshold = 1000.0; // RSLV: how to get user threshold?
		
		//upper_threshold = 1000.0; // RSLV: how to get user threshold?
		//lower_threshold = 650.0; // RSLV: how to get user threshold?
		
		// STEP 3a: first find all valid lines above threshold in descending order
		Vector<Vector<Double> > high_threshold_points = new Vector<Vector<Double> >();
		for(int px = 0; px < image_width; ++px)
		{
			for(int py = 0; py < image_height; ++py)
			{
				// RSLV: do not use abs in non-abs modes?
				// RSLV: use only valid line points?
				if(Math.abs(valid_line_points_magnitude[px][py]) >= upper_threshold)
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
		Line[][] point_to_line_map = new Line[image_width][image_height]; // NOTE: for computational efficiency of finding which line belongs to a point
		int[][] processing_map = new int[image_width][image_height];
		int DEFAULT = 0x00;
		int AVAILABLE = 0x00;
		int PROCESSED = 0x01;
		int LINEPOINT = 0x02;
		int JUNCTION = 0x04;
		
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
			lines.add(line);
			
			// *****************************************************************
			
			System.err.println("Trace in first direction");
			
			// trace line in first direction
			int tp_px = cp_px;
			int tp_py = cp_py;
			
			//double tp_sxo = cp_sx;
			//double tp_syo = cp_sy;
			double tp_sao = cp_sa;
			//int tp_sio = cp_si;
			
			boolean tracing = true;
			//int debug_index = 0;
			//while(tracing && debug_index++ < 1000)
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
				// RSLV: redudant because of check in neighbour orientation?
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
					
					// RSLV: only process valid line point
					if(valid_line_points_mask[dp_px][dp_py])
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
						if(dtheta > 90.0) // NOTE: bigger than, not equal!
						{
							dtheta = 180.0 - dtheta; // RSLV: check if correct?
							System.err.println("corrected dtheta="+dtheta);
						}
						dtheta *= Math.PI / 180.0; // deg to rad
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
					}
				}
				
				// *************************************************************
				
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
				if(Math.abs(np_sa - tp_sao) > 90.0)
				{
					System.err.println("Reorienting s(t) vector");
					np_sx = -np_sx;
					np_sy = -np_sy;
					np_sa = getAngle(np_sx, np_sy);
					np_si = getOrientationIndex(np_sa);
				}
				
				// determine action
				System.err.println("Determine best action");
				if(processing_map[np_px][np_py] == AVAILABLE)
				{
					System.err.println("  Point is available, add point to line");
					
					// skip if next neighbour below lowest threshold level
					double np_pv = line_points[np_px][np_py][0];
					if(Math.abs(np_pv) >= lower_threshold) // RSLV: do not use absolute value in non-absolute modes
					{
						// add point to line
						line.addFirst(new Point(np_px, np_py, np_dlpx, np_dlpy));
						point_to_line_map[np_px][np_py] = line;
						processing_map[np_px][np_py] = LINEPOINT;
						processing_map_ip.set(np_px, np_py, LINEPOINT);
						
						// mark perpendicular neighbours with similar orientation as processed
						int[] pp = getNextPixel(np_ni);
						int pp1_dx = pp[0];
						int pp1_dy = pp[1];
						int pp2_dx = -pp1_dx; // NOTE: pp2_dx = -pp1_dx
						int pp2_dy = -pp1_dy; // NOTE: pp2_dy = -pp1_dy
						/*switch(np_ni) // NOTE: np_ni -> normals in perpendicular direction
						{
							case 0:
								pp1_dx = 1;
								pp1_dy = 0;
								pp2_dx = -1;
								pp2_dy = 0;
								break;
							case 1:
								pp1_dx = 1;
								pp1_dy = 1;
								pp2_dx = -1;
								pp2_dy = -1;
								break;
							case 2:
								pp1_dx = 1;
								pp1_dy = 0;
								pp2_dx = -1;
								pp2_dy = 0;
								break;
							case 3:
								pp1_dx = -1;
								pp1_dy = 1;
								pp2_dx = 1;
								pp2_dy = -1;
								break;
							case 4:
								pp1_dx = -1;
								pp1_dy = 0;
								pp2_dx = 1;
								pp2_dy = 0;
								break;
							case 5:
								pp1_dx = -1;
								pp1_dy = -1;
								pp2_dx = 1;
								pp2_dy = 1;
								break;
							case 6:
								pp1_dx = 0;
								pp1_dy = -1;
								pp2_dx = 0;
								pp2_dy = 1;
								break;
							case 7:
								pp1_dx = 1;
								pp1_dy = -1;
								pp2_dx = -1;
								pp2_dy = 1;
								break;
						}*/
						
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
							if(processing_map[pp1_px][pp1_py] == AVAILABLE && Math.abs((pp1_na % 180.0) - (np_na % 180.0)) <= ROUGHLY_PARALLEL_THRESHOLD)
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
							if(processing_map[pp2_px][pp2_py] == AVAILABLE && Math.abs((pp2_na % 180.0) - (np_na % 180.0)) <= ROUGHLY_PARALLEL_THRESHOLD)
							{
								System.err.println("Checkpoint E");
								processing_map[pp2_px][pp2_py] = PROCESSED;
								processing_map_ip.set(pp2_px, pp2_py, PROCESSED);
							}
						}
					}
					else
					{
						// stop tracing
						tracing = false;
					}
				}
				else if(processing_map[np_px][np_py] == PROCESSED)
				{
					System.err.println("  Point is already processed, stop tracing");
					// stop tracing
					tracing = false;
				}
				else if(processing_map[np_px][np_py] == LINEPOINT)
				{
					System.err.println("  Point is part of line, (TODO) mark as junction, and stop tracing");
					// TODO: split the other line
					// TODO: mark point as junction
					
					// stop tracing
					tracing = false;
				}
				else if(processing_map[np_px][np_py] == JUNCTION)
				{
					System.err.println("  Point is junction, stop tracing");
					// RSLV: new case, how to continue?
					tracing = false;
				}
				
				// continue tracing from neighbour
				tp_px = np_px;
				tp_py = np_py;
				tp_sao = np_sa;
			}
			
			System.err.println("Done tracing in first direction");
			
			// *****************************************************************
			
			System.err.println("Trace in second direction");
			
			// trace line in second direction
			tp_px = cp_px; //int tp_px = cp_px;
			tp_py = cp_py; //int tp_py = cp_py;
			
			tp_sao = (cp_sa + 180.0) % 360.0; // NOTE: reverse direction of search!
			
			//boolean tracing = true;
			tracing = true;
			//debug_index = 0;
			//while(tracing && debug_index++ < 1000)
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
				// RSLV: redudant because of check in neighbour orientation?
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
					
					// RSLV: only process valid line point
					if(valid_line_points_mask[dp_px][dp_py])
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
						if(dtheta > 90.0) // NOTE: bigger than, not equal!
						{
							dtheta = 180.0 - dtheta; // RSLV: check if correct?
							System.err.println("corrected dtheta="+dtheta);
						}
						dtheta *= Math.PI / 180.0; // deg to rad
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
					}
				}
				
				// *************************************************************
				
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
				if(Math.abs(np_sa - tp_sao) > 90.0)
				{
					System.err.println("Reorienting s(t) vector");
					np_sx = -np_sx;
					np_sy = -np_sy;
					np_sa = getAngle(np_sx, np_sy);
					np_si = getOrientationIndex(np_sa);
				}
				
				// determine action
				System.err.println("Determine best action");
				if(processing_map[np_px][np_py] == AVAILABLE)
				{
					System.err.println("  Point is available, add point to line");
					
					// skip if next neighbour below lowest threshold level
					double np_pv = line_points[np_px][np_py][0];
					if(Math.abs(np_pv) >= lower_threshold) // RSLV: do not use absolute value in non-absolute modes
					{
						// add point to line
						line.addLast(new Point(np_px, np_py, np_dlpx, np_dlpy));
						point_to_line_map[np_px][np_py] = line;
						processing_map[np_px][np_py] = LINEPOINT;
						processing_map_ip.set(np_px, np_py, LINEPOINT);
						
						// mark perpendicular neighbours with similar orientation as processed
						int[] pp = getNextPixel(np_ni);
						int pp1_dx = pp[0];
						int pp1_dy = pp[1];
						int pp2_dx = -pp1_dx; // NOTE: pp2_dx = -pp1_dx
						int pp2_dy = -pp1_dy; // NOTE: pp2_dy = -pp1_dy
						/*switch(np_ni) // NOTE: np_ni -> normals in perpendicular direction
						{
							case 0:
								pp1_dx = 1;
								pp1_dy = 0;
								pp2_dx = -1;
								pp2_dy = 0;
								break;
							case 1:
								pp1_dx = 1;
								pp1_dy = 1;
								pp2_dx = -1;
								pp2_dy = -1;
								break;
							case 2:
								pp1_dx = 1;
								pp1_dy = 0;
								pp2_dx = -1;
								pp2_dy = 0;
								break;
							case 3:
								pp1_dx = -1;
								pp1_dy = 1;
								pp2_dx = 1;
								pp2_dy = -1;
								break;
							case 4:
								pp1_dx = -1;
								pp1_dy = 0;
								pp2_dx = 1;
								pp2_dy = 0;
								break;
							case 5:
								pp1_dx = -1;
								pp1_dy = -1;
								pp2_dx = 1;
								pp2_dy = 1;
								break;
							case 6:
								pp1_dx = 0;
								pp1_dy = -1;
								pp2_dx = 0;
								pp2_dy = 1;
								break;
							case 7:
								pp1_dx = 1;
								pp1_dy = -1;
								pp2_dx = -1;
								pp2_dy = 1;
								break;
						}*/
						
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
							System.err.println("pp1_na="+pp1_na+"    np_na="+np_na);
							if(processing_map[pp1_px][pp1_py] == AVAILABLE && Math.abs((pp1_na % 180.0) - (np_na % 180.0)) <= ROUGHLY_PARALLEL_THRESHOLD)
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
							System.err.println("pp2_na="+pp2_na+"    np_na="+np_na);
							if(processing_map[pp2_px][pp2_py] == AVAILABLE && Math.abs((pp2_na % 180.0) - (np_na % 180.0)) <= ROUGHLY_PARALLEL_THRESHOLD)
							{
								System.err.println("Checkpoint E");
								processing_map[pp2_px][pp2_py] = PROCESSED;
								processing_map_ip.set(pp2_px, pp2_py, PROCESSED);
							}
						}
					}
					else
					{
						// stop tracing
						tracing = false;
					}
				}
				else if(processing_map[np_px][np_py] == PROCESSED)
				{
					System.err.println("  Point is already processed, stop tracing");
					// stop tracing
					tracing = false;
				}
				else if(processing_map[np_px][np_py] == LINEPOINT)
				{
					System.err.println("  Point is part of line, (TODO) mark as junction, and stop tracing");
					
					// TODO: split the other line
					
					// TODO: mark point as junction
					
					// stop tracing
					tracing = false;
				}
				else if(processing_map[np_px][np_py] == JUNCTION)
				{
					System.err.println("  Point is junction, stop tracing");
					// RSLV: new case, how to continue?
					tracing = false;
				}
				
				// continue tracing from neighbour
				tp_px = np_px;
				tp_py = np_py;
				tp_sao = np_sa;
			}
			
			System.err.println("Done tracing second direction");
		}
		
		// ------
		
		ImagePlus processing_map_imp = new ImagePlus("Processing map", processing_map_ip);
		processing_map_imp.show();
		
		// *********************************************************************
		
		// draw overlay on duplicate of original image
		ImageProcessor overlay_ip = ip.duplicate();
		ImageProcessor overlay2_ip = valid_line_points_magnitude_ip.duplicate();
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
			
			++line_color_index;
		}
		
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
		
		// RSLV: return image with overlay?
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
			/*case 0:
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
				break;*/
			
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
			
			/*case 0:
				dp[0] = 1;
				dp[1] = 0;
				break;
			case 1:
				dp[0] = 1;
				dp[1] = 1;
				break;
			case 2:
				dp[0] = 1;
				dp[1] = 0;
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
				break;*/
		}
		return dp;
	}
}
