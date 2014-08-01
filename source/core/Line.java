
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
	
	public double contourLength(int start, int end)
	{
		double result = 0.0;
		for(int i = start; i <= end; ++i)
		{
			result += euclideanDistance(this.get(i), this.get(i+1));
		}
		return result;
	}
	
	// end-to-end distance
	public endToEndDistance()
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
	public double tangent(Point p1, Point p2)
	{
		double dx = (p1.px+p1.sx) - (p2.px+p2.sx);
		double dy = (p1.py+p1.sy) - (p2.py+p2.sy);
		//double dz = (p1.pz+p1.sz) - (p2.pz+p2.sz);
		
		return Math.atan2(dy, dx);
	}
	
	// persistence length
	public double persistenceLength()
	{
		return persistenceLength(0, this.size()-1);
	}
	
	public double persistenceLength(int start, int end)
	{
		// < cos ( theta ) > = exp ( - L_c / L_p )
		// L_p = -L_c / ln( < cos ( theta ) > ); RSLV: cos ( theta ) can be negative?
		double t1 = tangent(this.get(start), this.get(start+1));
		double t2 = tangent(this.get(end-1), this.get(end));
		double theta = 0.0;
		double lc = contourLength(start, end);
		return -lc / Math.ln(Math.cos(theta)); // RSLV: factor 2 because we are working in 2D?
	}
}
