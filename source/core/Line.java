
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
}
