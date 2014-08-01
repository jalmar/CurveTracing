
package core;

// import Java classes
import java.util.List;
import java.util.LinkedList;
import java.util.ListIterator;


/**
 *	Line class, contains an ordered list of connected points
 */
public class Line extends LinkedList<Point>
{
	/**
	 *
	 */
	
	// ////////////////////////////////////////////////////////////////////////
	
	/**
	 *
	 */
	public Line()
	{
		super();
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public Line splice(int position)
	{
		if(position < 0 || position >= this.size()) return null;
		
		Line second_half = new Line();
		while(this.size() > position)
		{
			Point p = this.pollLast();
			second_half.addFirst(p); // Note: retain order
		}
		return second_half;
	}
	
	// reverse
	public void reverse()
	{
		for(int i = 0; i < this.size(); ++i)
		{
			Point p = this.pollLast();
			this.addFirst(p);
		}
	}
	
	// contour length
	public double contourLength()
	{
		return contourLength(0, this.size()-1);
	}
	
	// NOTE: start and end are inclusive
	public double contourLength(int start, int end)
	{
		double result = 0.0;
		for(int i = start; i <= end - 1; ++i)
		{
			result += euclideanDistance(this.get(i), this.get(i+1));
		}
		return result;
	}
	
	// end-to-end distance
	public double endToEndDistance()
	{
		return euclideanDistance(this.get(0), this.get(this.size()-1));
	}
	
	// euclidean distance
	public double euclideanDistance(Point p1, Point p2)
	{
		return euclideanDistance(p1.px+p1.sx, p1.py+p1.sy, p1.pz+p1.sz, p2.px+p2.sx, p2.py+p2.sy, p2.pz+p2.sz);
	}
	
	public double euclideanDistance(double x1, double y1, double x2, double y2)
	{
		return euclideanDistance(x1, y1, x2, y2, 0.0, 0.0);
	}
	
	public double euclideanDistance(double x1, double y1, double z1, double x2, double y2, double z2)
	{
		return Math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
	}
	
	// tangent
	// NOTE: only in 2D!
/*	public double tangent(Point p1, Point p2)
	{
		double dx = (p1.px+p1.sx) - (p2.px+p2.sx);
		double dy = (p1.py+p1.sy) - (p2.py+p2.sy);
		//double dz = (p1.pz+p1.sz) - (p2.pz+p2.sz);
		
		return Math.atan2(dy, dx);
	}
*/	
	// persistence length
	public double persistenceLength()
	{
		return persistenceLength(0, this.size()-1);
	}
	
	// NOTE: start and end are inclusive
	public double persistenceLength(int start, int end)
	{
		// < cos ( theta ) > = exp ( - L_c / L_p )
		// L_p = -L_c / ln( < cos ( theta ) > ); RSLV: cos ( theta ) can be negative?
		// RSLv: prefer L_c / L_p
		double theta_sum = 0.0;
		double theta_squared_sum = 0.0;
		int count_n = 0;
		for(int i = start; i <= end - 2; ++i)
		{
			// dot product
			double[] v1 = segmentVector(i, true); // NOTE: normalized
			double[] v2 = segmentVector(i+1, true); // NOTE: normalized
			double theta_cos = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
			// RLSV: theta_cos = Math.abs(theta_cos)?
			theta_sum += theta_cos;
			theta_squared_sum += theta_cos * theta_cos;
			++count_n
		}
		double theta_mean = theta_sum / count_n; // RSLV: count_n = start - end - (2|1)
		double theta_var = (theta_squared_sum / count_n) - theta_mean*theta_mean;
		//double theta_stdev = Math.sqrt(theta_var);
		double lc = contourLength(start, end);
		double ls = lc / ls; // approximately / average segment length
		//return (-lc / Math.log(theta_mean));
		return (2*ls) / theta_var; // stdev = Math.sqrt(2 * ls / lp) => lp = (2 * ls) / stdev^2 = 2*ls / var
	}
	
	// NOTE: get vector of segment (index,index+1)
	// NOTE: range of index is 0 to size()-2
	public double[] segmentVector(int index)
	{
		return segmentVector(index, false);
	}
	
	public double[] segmentVector(int index, boolean normalize)
	{
		Point p1 = this.get(index);
		Point p2 = this.get(index+1);
		double dx = (p2.px+p2.sx) - (p1.px+p1.sx);
		double dy = (p2.py+p2.sy) - (p1.py+p1.sy);
		double dz = (p2.pz+p2.sz) - (p1.pz+p1.sz);
		if(normalize)
		{
			double dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
			if(dist != 0)
			{
				dx /= dist;
				dy /= dist;
				dz /= dist;
			}
			else
			{
				System.err.println("WARNING: dividion by zero in Line::segmentVector()");
			}
		}
		return new double[]{dx, dy, dz};
	}
	
}
