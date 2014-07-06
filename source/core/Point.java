
package core;

// import Java classes


/**
 *	Point class, contains 3D pixel coordinate data, 3D super-resolved coordinate data, and flag attributes
 */
public class Point
{
	/**
	 *
	 */
	public int px, py, pz;
	public double sx, sy, sz; // (0,0,0) is center of pixel
	private int flags;
	
	// ////////////////////////////////////////////////////////////////////////
	
	/**
	 *
	 */
	public Point()
	{
		this(0, 0, 0);
	}
	
	/**
	 *
	 */
	public Point(int npx, int npy)
	{
		this(npx, npy, 0);
	}
	
	/**
	 *
	 */
	public Point(int npx, int npy, int npz)
	{
		this(npx, npy, npz, 0.0, 0.0, 0.0); // RSLV: default to nsx=npx, nsy=npy, nsz=npz?
	}
	
	/**
	 *
	 */
	public Point(int npx, int npy, double nsx, double nsy)
	{
		this(npx, npy, 0, nsx, nsy, 0.0);
	}
	
	/**
	 *
	 */
	public Point(int npx, int npy, int npz, double nsx, double nsy, double nsz)
	{
		this(npx, npy, npz, nsx, nsy, nsz, 0x00);
	}
	
	/**
	 *
	 */
	public Point(int npx, int npy, int npz, double nsx, double nsy, double nsz, int nflags)
	{
		px = npx; py = npy; pz = npz;
		sx = nsx; sy = nsy; sz = nsz;
		flags = nflags;
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	public boolean checkFlag(int flag)
	{
		return (flags & flag) != 0;
	}
	
	public void setFlag(int flag)
	{
		flags |= flag;
	}
	
	public void clearFlag(int flag)
	{
		flags &= ~flag;
	}
	
	public void clearFlags()
	{
		flags = 0x00;
	}
	
	// ////////////////////////////////////////////////////////////////////////
	
	@Override
	public boolean equals(Object o)
	{
		Point p = (Point)o;
		return p.px == px && p.py == py && p.pz == pz;// && p.sx == sx && p.sy == sy && p.sz == sz; // NOTE: no flags are checked
	}
	
	public String toString()
	{
		return "{" + px + "," + py + "," + pz + ";" + sx + "," + sy + "," + sz + "}";
	}
}
