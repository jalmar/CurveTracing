
package utils;

/**
 *	Note: class implementation adapted from http://stackoverflow.com/questions/156275/what-is-the-equivalent-of-the-c-pairl-r-in-java and http://stackoverflow.com/questions/2670982/using-pairs-or-2-tuples-in-java
 */
public class Triple<A, B, C>
{
	// members
	public A first;
	public B second;
	public C third;
	
	// *************************************************************************
	
	// constructor
	public Triple()
	{
		super();
		this.first = null;
		this.second = null;
		this.third = null;
	}
	
	// NOTE: cannot do single argument constructor in case type A == type B
	
	public Triple(A first, B second, C third)
	{
		super();
		this.first = first;
		this.second = second;
		this.third = third;
	}
	
	// *************************************************************************
	
	public A getFirst()
	{
		return this.first;
	}
	
	public void setFirst(A first)
	{
		this.first = first;
	}
	
	public B getSecond()
	{
		return this.second;
	}
	
	public void setSecond(B second)
	{
		this.second = second;
	}
	
	public C getThird()
	{
		return this.third;
	}
	
	public void setThird(C third)
	{
		this.third = third;
	}
	
	// *************************************************************************
	
	@Override
	public int hashCode() {
		int hashFirst = this.first != null ? this.first.hashCode() : 0;
		int hashSecond = this.second != null ? this.second.hashCode() : 0;
		int hashThird = this.third != null ? this.third.hashCode() : 0;
		
		return (hashFirst + hashSecond + hashThird) * hashThird + (hashFirst + hashSecond) * hashSecond + hashFirst;
		
		//final int prime = 31;
		//int result = 1;
		//result = prime * result + ((x == null) ? 0 : x.hashCode());
		//result = prime * result + ((y == null) ? 0 : y.hashCode());
		//return result;
	}
	
	@Override
	public boolean equals(Object other)
	{
		if (other instanceof Triple) {
			Triple otherTriple = (Triple) other;
			return 
			((  this.first == otherTriple.first ||
				( this.first != null && otherTriple.first != null &&
				  this.first.equals(otherTriple.first))) &&
			 (	this.second == otherTriple.second ||
				( this.second != null && otherTriple.second != null &&
				  this.second.equals(otherTriple.second))) &&
			 (	this.third == otherTriple.third ||
				( this.third != null && otherTriple.third != null &&
				  this.third.equals(otherTriple.third))) );
		}
		
		return false;
	}
	
	public String toString()
	{ 
		return "(" + this.first + ", " + this.second + ", " + this.third + ")"; 
	}
}
