
package utils;

/**
 *	
 */
public class Profiling
{
	/**
	 *
	 */
	public static long last_tic = 0l;
	public static long last_toc = 0l;
	/**
	 *	Constructor
	 */
	public Profiling()
	{
		/* do nothing */
	}
	
	// ////////////////////////////////////////////////////////////////////////

	public static long tic()
	{
		last_tic = System.nanoTime();
		return last_tic;
	}
	
	public static long toc()
	{
		last_toc = System.nanoTime();
		return last_toc;
	}
	
	public static long toc(String message)
	{
		toc();
		print(message);
		return last_toc;
	}
	
	public static void print()
	{
		print("Process");
	}
	
	public static void print(String message)
	{
		System.err.println(message + " took " + format_time((last_toc - last_tic)));
	}
	
	public static String format_time(double time)
	{
		final String[] units = new String[]{"ns", "us", "ms", "sec", "min", "hours", "too long"};
		int unit = 0;
		while(time >= 1000 && unit < 3)
		{
			time /= 1000;
			++unit;
		}
		while(time >= 60 && unit >= 3 && unit < 5)
		{
			time /= 60;
			++unit;
		}
		return String.format("%.2f %s", time, units[unit]);
	}
}
